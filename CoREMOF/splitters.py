"""Deterministic, parent-aware train/validation/test splitting.

This module intentionally has no NumPy, pandas, scikit-learn, or chemistry
imports.  It can be used on a login node or in a lightweight downstream ML
environment as long as release metadata have already been loaded.

The project-defined default ``priority_main`` is the conflict-aware hierarchy
over exact RAC5, then MOFid-v2, then MOFid-v1 groups; conflicting lower groups
never merge stronger components, and missing evidence becomes a singleton or
is explicitly excluded.  It does not use optional reference methods.
``main_union`` is the broader split-leakage guard, not a parent or explanatory
method and not proof of identity or parentage: before filtering it takes
transitive connected components over the complete release from exact full CIF
SHA-256, database-namespaced source, RAC5, MOFid-v2, and MOFid-v1 edges.
``parent_only`` is the narrower guard that uses only groups from the selected
explanatory parent relation as split blocks.  ``auto`` selects ``main_union``
for ``priority_main`` and ``parent_only`` for every explicitly selected direct
or reference method.
"""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass, field, fields
from decimal import Decimal
import hashlib
import json
import math
from numbers import Real
import os
from pathlib import Path
import re
import shutil
import tempfile
from types import MappingProxyType
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

from CoREMOF import __version__
from CoREMOF._authority import (
    AuthorityStateError,
    IdentitySealRegistry,
    reject_sealed_copy,
    state_fingerprint,
)
from CoREMOF.dataset import (
    _is_authenticated_official_checker_view,
    _release_receipt_state,
    _validate_classified_generation_if_present,
    _validate_dataset_generation_if_present,
)
from CoREMOF.labels import CHECKER_PRESETS, LABELS, resolve_checker_view
from CoREMOF.parents import (
    LEAKAGE_GUARD_CHOICES,
    SELECTABLE_PARENT_METHODS,
    ParentResolver,
    generated_output_terminology_definitions,
    leakage_guard_definition,
    parent_method_definition,
    resolve_leakage_guard,
)


_SPLIT_NAMES = ("train", "validation", "test")
SPLIT_API_VERSION = "0.2.0"
_SPLIT_RESULT_FACTORY_TOKEN = object()
_SPLIT_RESULTS = IdentitySealRegistry("split result generation")


class OfficialSplitUnavailableError(RuntimeError):
    """Raised when an official split is requested without an audited manifest."""


def _deep_freeze(value):
    if isinstance(value, Mapping):
        result = {}
        for key, item in value.items():
            if not isinstance(key, str) or not key or key != key.strip():
                raise TypeError("mapping keys must be exact nonblank strings")
            result[key] = _deep_freeze(item)
        return MappingProxyType(result)
    if isinstance(value, (list, tuple, set, frozenset)):
        return tuple(_deep_freeze(item) for item in value)
    return value


def _jsonable(value):
    if isinstance(value, Mapping):
        result = {}
        for key, item in value.items():
            if not isinstance(key, str) or not key or key != key.strip():
                raise TypeError("mapping keys must be exact nonblank strings")
            result[key] = _jsonable(item)
        return result
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_jsonable(item) for item in value]
    return value


_BASE_IMPLEMENTATION_FILES = (
    "_authority.py",
    "dataset.py",
    "labels.py",
    "parents.py",
    "splitters.py",
)
_TARGET_IMPLEMENTATION_FILES = (
    "_authority.py",
    "dataset.py",
    "labels.py",
    "targets.py",
)


def _current_base_implementation_hashes() -> Dict[str, str]:
    package_root = Path(__file__).resolve().parent
    return {
        filename: hashlib.sha256(
            (package_root / filename).read_bytes()
        ).hexdigest()
        for filename in _BASE_IMPLEMENTATION_FILES
    }


_IMPORTED_BASE_IMPLEMENTATION_HASHES = MappingProxyType(
    _current_base_implementation_hashes()
)


def _target_receipt_implementation_hashes(
    target_receipt: Mapping[str, object],
) -> Dict[str, str]:
    implementation = target_receipt.get("implementation")
    if not isinstance(implementation, Mapping):
        raise ValueError("target merge receipt has no implementation mapping")
    source_hashes = implementation.get("source_sha256")
    if not isinstance(source_hashes, Mapping):
        raise ValueError("target merge receipt has no source_sha256 mapping")
    observed = {str(filename) for filename in source_hashes}
    if observed != set(_TARGET_IMPLEMENTATION_FILES):
        raise ValueError(
            "target merge implementation closure must contain exactly {}".format(
                ", ".join(_TARGET_IMPLEMENTATION_FILES)
            )
        )
    result = {}
    for filename in _TARGET_IMPLEMENTATION_FILES:
        digest = source_hashes.get(filename)
        if (
            not isinstance(digest, str)
            or len(digest) != 64
            or any(character not in "0123456789abcdef" for character in digest)
        ):
            raise ValueError(
                "target merge receipt has an invalid SHA-256 for {}".format(
                    filename
                )
            )
        result[filename] = digest
    return result


def _implementation_hashes(
    include_targets: bool = False,
    target_receipt: Optional[Mapping[str, object]] = None,
) -> Mapping[str, str]:
    """Return import-bound split sources and the frozen target-merge source."""

    current = _current_base_implementation_hashes()
    imported = dict(_IMPORTED_BASE_IMPLEMENTATION_HASHES)
    if current != imported:
        changed = sorted(
            filename
            for filename in set(current).union(imported)
            if current.get(filename) != imported.get(filename)
        )
        raise ValueError(
            "CoREMOF split implementation source changed after module import: {}".format(
                ", ".join(changed)
            )
        )
    result = dict(imported)
    if not include_targets:
        return result
    if target_receipt is None:
        raise ValueError(
            "target-aware implementation hashing requires the target merge receipt"
        )
    target_hashes = _target_receipt_implementation_hashes(target_receipt)
    for filename in ("_authority.py", "dataset.py", "labels.py"):
        if target_hashes[filename] != imported[filename]:
            raise ValueError(
                "{} changed between target merge and split import".format(filename)
            )
    result["targets.py"] = target_hashes["targets.py"]
    return result


def _as_optional_tuple(values, name: str) -> Optional[Tuple[str, ...]]:
    if values is None:
        return None
    if isinstance(values, str):
        result = (values,)
    else:
        if not isinstance(values, (list, tuple)):
            raise TypeError(
                "%s must be an exact string or ordered list/tuple of strings"
                % name
            )
        result = tuple(values)
    if not result:
        raise ValueError("%s must not be empty; use None for no filter" % name)
    if any(not isinstance(value, str) for value in result):
        raise TypeError("%s values must be strings" % name)
    if any(not value or value != value.strip() for value in result):
        raise ValueError("%s values must be exact nonblank strings" % name)
    return result


def _column_name_tuple(value: object, name: str) -> Tuple[str, ...]:
    """Validate one optional externally supplied column-name sequence."""

    if value is None:
        return ()
    if isinstance(value, str) or not isinstance(value, (list, tuple)):
        raise TypeError("%s must be a list or tuple of strings when present" % name)
    result = tuple(value)
    if any(
        not isinstance(item, str) or not item or item != item.strip()
        for item in result
    ):
        raise ValueError("%s values must be exact nonblank strings" % name)
    if len(set(result)) != len(result):
        raise ValueError("%s must not contain duplicate names" % name)
    return result


def _stable_digest(*values: object) -> str:
    text = "\0".join(str(value) for value in values)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _split_source_payload(classified: object) -> Tuple[object, ...]:
    """Return exact generic-protocol state on which a split can depend."""

    dataset = getattr(classified, "dataset", None)
    if dataset is None:
        raise AuthorityStateError("split source has no dataset")
    metadata_rows = getattr(dataset, "metadata_rows", None)
    if not isinstance(metadata_rows, (list, tuple)):
        raise AuthorityStateError("split source metadata_rows must be a sequence")
    if any(not isinstance(row, Mapping) for row in metadata_rows):
        raise AuthorityStateError("split source metadata rows must be mappings")
    labels = getattr(classified, "label_by_id", None)
    if labels is None:
        labels = getattr(classified, "labels", None)
    if not isinstance(labels, (Mapping, list, tuple)):
        raise AuthorityStateError("split source labels must be a mapping or sequence")
    official = getattr(classified, "checker_view_official", False)
    if type(official) is not bool:
        raise AuthorityStateError("checker_view_official must be a boolean")
    return (
        "coremof-generic-split-source/1",
        getattr(classified, "checker_view", None),
        official,
        tuple(getattr(classified, "structure_ids", ())),
        labels,
        getattr(classified, "selection_filters", None),
        tuple(metadata_rows),
        getattr(dataset, "parent_by_id", None),
        getattr(dataset, "cif_hashes", None),
        getattr(dataset, "input_hashes", None),
        getattr(dataset, "dataset_info", None),
        getattr(dataset, "parent_group_methods", None),
        getattr(dataset, "dataset_version", None),
        getattr(dataset, "cif_files_verified", None),
        getattr(dataset, "target_columns", None),
        getattr(dataset, "feature_columns", None),
        getattr(dataset, "target_input_receipt", None),
    )


def _split_result_identity_snapshot(result: object) -> Tuple[object, ...]:
    snapshot = []
    for definition in fields(result):
        name = definition.name
        if name == "_authority_factory_token":
            continue
        value = getattr(result, name)
        if name == "_authority_source":
            snapshot.append((name, "identity", id(value)))
        elif isinstance(value, Mapping):
            snapshot.append((name, "mapping", id(value)))
        elif type(value) is tuple:
            snapshot.append((name, "tuple", id(value)))
        else:
            snapshot.append((name, type(value).__module__, type(value).__qualname__, value))
    return tuple(snapshot)


def _split_result_seal_fingerprint(result: object) -> str:
    """Bind every public dataclass field, including nested receipt evidence."""

    payload = tuple(
        (definition.name, getattr(result, definition.name))
        for definition in fields(result)
        if definition.name not in {"_authority_factory_token", "_authority_source"}
    )
    return state_fingerprint(("coremof-split-result-generation/2", payload))


@dataclass(frozen=True)
class SplitResult:
    """A complete split assignment and its reproducibility receipt.

    The object can be unpacked like MOFDescribe's splitters; unpacking yields
    ``train_indices, validation_indices, test_indices``.  Structure IDs remain
    available explicitly and should be preferred when persisting a split.
    """

    train_ids: Tuple[str, ...]
    validation_ids: Tuple[str, ...]
    test_ids: Tuple[str, ...]
    train_indices: Tuple[int, ...]
    validation_indices: Tuple[int, ...]
    test_indices: Tuple[int, ...]
    assignments: Mapping[str, str]
    exclusions: Mapping[str, str]
    labels: Mapping[str, str]
    parent_groups: Mapping[str, str]
    parent_diagnostics: Mapping[str, str]
    parent_conflicts: Tuple[Mapping[str, object], ...]
    leakage_groups: Mapping[str, str]
    index_by_id: Mapping[str, int]
    view_index_by_id: Mapping[str, int]
    dataset_version: Optional[str]
    checker_view: Optional[str]
    checker_view_official: bool
    parent_method: str
    requested_leakage_guard: str
    leakage_guard: str
    missing_parent: str
    fractions: Tuple[float, float, float]
    random_state: str
    stratify_by: Tuple[str, ...]
    filters: Mapping[str, object]
    input_hashes: Mapping[str, str]
    implementation_package_version: str
    implementation_split_api_version: str
    implementation_hashes: Mapping[str, str]
    cif_files_verified: bool
    parent_input_status: Optional[str]
    dataset_input_status: Optional[str]
    provisional_input: bool
    target_data: Optional[Mapping[str, object]] = None
    official_split: bool = False
    _authority_factory_token: object = field(
        default=None, repr=False, compare=False
    )
    _authority_source: object = field(default=None, repr=False, compare=False)

    def __post_init__(self) -> None:
        for field_name in (
            "train_ids",
            "validation_ids",
            "test_ids",
            "train_indices",
            "validation_indices",
            "test_indices",
            "fractions",
            "stratify_by",
        ):
            object.__setattr__(self, field_name, tuple(getattr(self, field_name)))
        for field_name in (
            "assignments",
            "exclusions",
            "labels",
            "parent_groups",
            "parent_diagnostics",
            "leakage_groups",
            "index_by_id",
            "view_index_by_id",
            "input_hashes",
            "implementation_hashes",
        ):
            object.__setattr__(
                self,
                field_name,
                MappingProxyType(dict(getattr(self, field_name))),
            )
        object.__setattr__(self, "filters", _deep_freeze(self.filters))
        if self.target_data is not None:
            object.__setattr__(self, "target_data", _deep_freeze(self.target_data))
        object.__setattr__(
            self,
            "parent_conflicts",
            tuple(_deep_freeze(conflict) for conflict in self.parent_conflicts),
        )

    def __copy__(self):
        reject_sealed_copy(self, "copied")
        raise TypeError("unsealed split results cannot be copied")

    def __deepcopy__(self, memo):
        reject_sealed_copy(self, "deep-copied")
        raise TypeError("unsealed split results cannot be deep-copied")

    def __reduce_ex__(self, protocol):
        reject_sealed_copy(self, "pickled")
        raise TypeError("split results are point-in-time non-serializable objects")

    def __iter__(self):
        return iter((self.train_indices, self.validation_indices, self.test_indices))

    @property
    def valid_ids(self) -> Tuple[str, ...]:
        """MOFDescribe-style alias for ``validation_ids``."""

        return self.validation_ids

    @property
    def valid_indices(self) -> Tuple[int, ...]:
        """MOFDescribe-style alias for ``validation_indices``."""

        return self.validation_indices

    @property
    def train_release_indices(self) -> Tuple[int, ...]:
        return tuple(self.index_by_id[value] for value in self.train_ids)

    @property
    def validation_release_indices(self) -> Tuple[int, ...]:
        return tuple(self.index_by_id[value] for value in self.validation_ids)

    @property
    def test_release_indices(self) -> Tuple[int, ...]:
        return tuple(self.index_by_id[value] for value in self.test_ids)

    @property
    def counts(self) -> Mapping[str, int]:
        return {
            "train": len(self.train_ids),
            "validation": len(self.validation_ids),
            "test": len(self.test_ids),
            "excluded": len(self.exclusions),
        }

    @property
    def achieved_fractions(self) -> Mapping[str, float]:
        """Fractions among assigned (non-excluded) structures."""

        assigned = len(self.train_ids) + len(self.validation_ids) + len(self.test_ids)
        if assigned == 0:
            return {split: 0.0 for split in _SPLIT_NAMES}
        return {
            "train": len(self.train_ids) / assigned,
            "validation": len(self.validation_ids) / assigned,
            "test": len(self.test_ids) / assigned,
        }

    @property
    def label_counts_by_split(self) -> Mapping[str, Mapping[str, int]]:
        """Return the achieved checker-label distribution in each partition."""

        ids_by_split = {
            "train": self.train_ids,
            "validation": self.validation_ids,
            "test": self.test_ids,
        }
        return {
            split: {
                label: sum(self.labels.get(structure_id) == label for structure_id in ids)
                for label in sorted(
                    {
                        self.labels[structure_id]
                        for structure_id in ids
                        if structure_id in self.labels
                    }
                )
            }
            for split, ids in ids_by_split.items()
        }

    @property
    def warnings(self) -> Tuple[str, ...]:
        """Stable machine-readable caveats for this result."""

        warnings = []
        if self.provisional_input:
            warnings.append("PROVISIONAL_PARENT_INPUT")
        assigned = len(self.train_ids) + len(self.validation_ids) + len(self.test_ids)
        if assigned == 0:
            warnings.append("NO_ELIGIBLE_STRUCTURES")
        else:
            achieved = self.achieved_fractions
            integer_tolerance = 1.0 / assigned
            if any(
                abs(achieved[split] - requested) > integer_tolerance
                for split, requested in zip(_SPLIT_NAMES, self.fractions)
            ):
                warnings.append("GROUP_CONSTRAINED_FRACTIONS")
            for split, requested in zip(_SPLIT_NAMES, self.fractions):
                if requested > 0.0 and self.counts[split] == 0:
                    warnings.append("EMPTY_{}_PARTITION".format(split.upper()))
        if any(reason == "PARENT_METHOD_CONFLICT" for reason in self.exclusions.values()):
            warnings.append("PARENT_METHOD_CONFLICTS_EXCLUDED")
        if any(
            reason == "PARENT_METHOD_CONFLICT"
            for reason in self.parent_diagnostics.values()
        ):
            warnings.append("PARENT_METHOD_CONFLICTS_GUARDED")
        if self.parent_conflicts:
            warnings.append("PARENT_METHOD_CONFLICT_LEDGER_PRESENT")
        return tuple(warnings)

    @property
    def leakage_audit(self) -> Mapping[str, object]:
        splits_by_block: Dict[str, Set[str]] = {}
        sizes: Dict[str, int] = {}
        for structure_id, split in self.assignments.items():
            block = self.leakage_groups[structure_id]
            splits_by_block.setdefault(block, set()).add(split)
            sizes[block] = sizes.get(block, 0) + 1
        cross_split = sorted(
            block for block, splits in splits_by_block.items() if len(splits) > 1
        )
        return {
            "block_count": len(splits_by_block),
            "cross_split_block_count": len(cross_split),
            "cross_split_blocks": cross_split,
            "max_block_size": max(sizes.values()) if sizes else 0,
            "passed": not cross_split,
        }

    def assignment_rows(self) -> List[Mapping[str, object]]:
        """Return deterministic rows used by :meth:`to_csv`."""

        _require_split_result(self)
        rows: List[Mapping[str, object]] = []
        for structure_id, index in sorted(self.index_by_id.items(), key=lambda item: item[1]):
            rows.append(
                {
                    "structure_id": structure_id,
                    "release_index": index,
                    "view_index": self.view_index_by_id.get(structure_id, ""),
                    "split": self.assignments.get(structure_id, "EXCLUDED"),
                    "label": self.labels.get(structure_id, ""),
                    "parent_group": self.parent_groups.get(structure_id, ""),
                    "parent_diagnostic": self.parent_diagnostics.get(structure_id, ""),
                    "leakage_group": self.leakage_groups.get(structure_id, ""),
                    "exclusion_reason": self.exclusions.get(structure_id, ""),
                }
            )
        return rows

    def _assignment_sha256(self) -> str:
        lines = []
        for structure_id in sorted(self.index_by_id):
            lines.append(
                "\t".join(
                    (
                        structure_id,
                        self.assignments.get(structure_id, "EXCLUDED"),
                        self.labels.get(structure_id, ""),
                        self.parent_groups.get(structure_id, ""),
                        self.parent_diagnostics.get(structure_id, ""),
                        self.leakage_groups.get(structure_id, ""),
                        self.exclusions.get(structure_id, ""),
                    )
                )
            )
        return hashlib.sha256("\n".join(lines).encode("utf-8")).hexdigest()

    def receipt(self) -> Mapping[str, object]:
        """Return a JSON-serializable reproducibility receipt."""

        _require_split_result(self)
        partitions = {
            "train": {
                "ids": list(self.train_ids),
                "indices": list(self.train_indices),
                "release_indices": list(self.train_release_indices),
            },
            "validation": {
                "ids": list(self.validation_ids),
                "indices": list(self.validation_indices),
                "release_indices": list(self.validation_release_indices),
            },
            "test": {
                "ids": list(self.test_ids),
                "indices": list(self.test_indices),
                "release_indices": list(self.test_release_indices),
            },
        }
        receipt = {
            "schema_version": "coremof-split-receipt/1.1",
            "implementation": {
                "package": "CoREMOF-tools",
                "package_version": self.implementation_package_version,
                "split_api_version": self.implementation_split_api_version,
                "source_sha256": dict(self.implementation_hashes),
            },
            "algorithm": {
                "name": "deterministic_greedy_group_stratification",
                "version": "1.0",
                "group_constraint_precedes_ratio_balance": True,
            },
            "dataset_version": self.dataset_version,
            "checker_view": self.checker_view,
            "checker_view_kind": (
                "OFFICIAL_RELEASE_VIEW"
                if self.checker_view_official
                else "USER_DEFINED"
            ),
            "contract_definitions": _jsonable(
                generated_output_terminology_definitions()
            ),
            "parent_input_status": self.parent_input_status,
            "dataset_input_status": self.dataset_input_status,
            "provisional_input": self.provisional_input,
            "official_split": self.official_split,
            "parent_method": self.parent_method,
            "parent_method_definition": _jsonable(
                parent_method_definition(self.parent_method)
            ),
            "requested_leakage_guard": self.requested_leakage_guard,
            "requested_leakage_guard_definition": _jsonable(
                leakage_guard_definition(self.requested_leakage_guard)
            ),
            "leakage_guard": self.leakage_guard,
            "leakage_guard_definition": _jsonable(
                leakage_guard_definition(self.leakage_guard)
            ),
            "missing_parent": self.missing_parent,
            "fractions": dict(zip(_SPLIT_NAMES, self.fractions)),
            "achieved_fractions": dict(self.achieved_fractions),
            "label_counts_by_split": {
                split: dict(counts)
                for split, counts in self.label_counts_by_split.items()
            },
            "random_state": self.random_state,
            "stratify_by": list(self.stratify_by),
            "filters": _jsonable(self.filters),
            "counts": dict(self.counts),
            "input_hashes": dict(sorted(self.input_hashes.items())),
            "integrity": {
                "cif_files_verified": self.cif_files_verified,
                "manifest_full_sha256_required_by_main_union": (
                    self.leakage_guard == "main_union"
                ),
            },
            "assignment_sha256": self._assignment_sha256(),
            "leakage_audit": dict(self.leakage_audit),
            "warnings": list(self.warnings),
            "parent_diagnostics": [
                {"structure_id": structure_id, "diagnostic": diagnostic}
                for structure_id, diagnostic in sorted(self.parent_diagnostics.items())
            ],
            "parent_conflict_count": len(self.parent_conflicts),
            "parent_conflicts": _jsonable(self.parent_conflicts),
            "partitions": partitions,
            "exclusions": [
                {
                    "structure_id": structure_id,
                    "release_index": self.index_by_id[structure_id],
                    "view_index": self.view_index_by_id.get(structure_id),
                    "reason": reason,
                }
                for structure_id, reason in sorted(
                    self.exclusions.items(), key=lambda item: self.index_by_id[item[0]]
                )
            ],
        }
        if self.target_data is not None:
            receipt["target_data"] = _jsonable(self.target_data)
        return receipt

    @staticmethod
    def _atomic_target(path, overwrite: bool) -> Tuple[Path, Path]:
        if type(overwrite) is not bool:
            raise TypeError("overwrite must be a boolean")
        output = Path(path)
        output.parent.mkdir(parents=True, exist_ok=True)
        if output.exists() and not overwrite:
            raise FileExistsError("Refusing to overwrite existing split output: %s" % output)
        descriptor, temporary_name = tempfile.mkstemp(
            prefix=".%s." % output.name, suffix=".tmp", dir=str(output.parent)
        )
        os.close(descriptor)
        return output, Path(temporary_name)

    @staticmethod
    def _publish_temporary(temporary: Path, output: Path, overwrite: bool) -> None:
        if overwrite:
            os.replace(str(temporary), str(output))
        else:
            # Linking within the destination directory is an atomic
            # create-if-absent operation.  Unlike os.replace(), it cannot
            # overwrite a file created by a concurrent writer after preflight.
            os.link(str(temporary), str(output))
            temporary.unlink()

    def to_csv(self, path, overwrite: bool = False) -> Path:
        """Write assignment rows, including explicit exclusions."""

        _require_split_result(self)
        output, temporary = self._atomic_target(path, overwrite)
        fieldnames = (
            "structure_id",
            "release_index",
            "view_index",
            "split",
            "label",
            "parent_group",
            "parent_diagnostic",
            "leakage_group",
            "exclusion_reason",
        )
        try:
            with temporary.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(self.assignment_rows())
                handle.flush()
                os.fsync(handle.fileno())
            self._publish_temporary(temporary, output, overwrite)
        finally:
            if temporary.exists():
                temporary.unlink()
        return output

    def to_json(self, path, overwrite: bool = False) -> Path:
        """Write the point-in-time split receipt."""

        _require_split_result(self)
        output, temporary = self._atomic_target(path, overwrite)
        try:
            with temporary.open("w", encoding="utf-8") as handle:
                json.dump(self.receipt(), handle, indent=2, sort_keys=True)
                handle.write("\n")
                handle.flush()
                os.fsync(handle.fileno())
            self._publish_temporary(temporary, output, overwrite)
        finally:
            if temporary.exists():
                temporary.unlink()
        return output

    @staticmethod
    def _quarantine_published_generation(
        target: Path, identity: Path, quarantine: Path
    ) -> bool:
        """Remove only the published generation identified by ``identity``.

        Moving the public pathname into the private staging directory before
        comparing inodes closes the race between ``samefile`` and ``unlink``.
        A foreign generation is restored create-if-absent; if a newer writer
        already owns the public pathname, the detached generation remains in
        staging and the caller preserves that directory for recovery.
        """

        try:
            os.replace(str(target), str(quarantine))
        except FileNotFoundError:
            return True
        try:
            published_by_us = identity.exists() and os.path.samefile(
                str(identity), str(quarantine)
            )
        except FileNotFoundError:
            published_by_us = False
        if published_by_us:
            quarantine.unlink()
            return True
        try:
            os.link(str(quarantine), str(target))
        except OSError:
            return False
        quarantine.unlink()
        return True

    def write(
        self,
        output_directory,
        stem: str = "coremof_split",
        overwrite: bool = False,
    ) -> Tuple[Path, Path]:
        """Write ``<stem>.csv`` and ``<stem>.json`` and return both paths.

        Default create-if-absent rollback never deletes a concurrent writer's
        replacement.  Since two ordinary files cannot be replaced as one
        filesystem transaction, ``overwrite=True`` has explicit single-writer
        semantics; callers must serialize writers for the same stem.
        """

        _require_split_result(self)
        if type(overwrite) is not bool:
            raise TypeError("overwrite must be a boolean")
        if (
            not isinstance(stem, str)
            or not stem
            or stem in {".", ".."}
            or Path(stem).name != stem
        ):
            raise ValueError("stem must be a non-empty filename stem without directories")
        directory = Path(output_directory)
        directory.mkdir(parents=True, exist_ok=True)
        csv_path = directory / (stem + ".csv")
        json_path = directory / (stem + ".json")
        if not overwrite:
            existing = [path for path in (csv_path, json_path) if path.exists()]
            if existing:
                raise FileExistsError(
                    "Refusing to overwrite existing split output: %s"
                    % ", ".join(str(path) for path in existing)
                )

        # Render the complete pair before touching either public destination.
        # Publishing two ordinary files cannot be one filesystem transaction,
        # so preserve/restore prior files if the second rename fails.
        staging_directory = Path(
            tempfile.mkdtemp(prefix=".{}.".format(stem), dir=str(directory))
        )
        staged_csv = staging_directory / csv_path.name
        staged_json = staging_directory / json_path.name
        targets = (csv_path, json_path)
        staged = (staged_csv, staged_json)
        backups = {}
        published_identities = {}
        published = []
        remove_staging = True
        try:
            self.to_csv(staged_csv)
            self.to_json(staged_json)
            if overwrite:
                for target in targets:
                    if target.exists():
                        backup = staging_directory / (target.name + ".previous")
                        try:
                            os.link(str(target), str(backup))
                        except OSError:
                            shutil.copy2(str(target), str(backup))
                        backups[target] = backup
            for source, target in zip(staged, targets):
                if overwrite:
                    if target not in backups:
                        identity = staging_directory / (target.name + ".published")
                        os.link(str(source), str(identity))
                        published_identities[target] = identity
                    os.replace(str(source), str(target))
                else:
                    os.link(str(source), str(target))
                published.append(target)
        except BaseException:
            for target in reversed(targets):
                backup = backups.get(target)
                if backup is not None and backup.exists():
                    os.replace(str(backup), str(target))
                elif target in published:
                    identity = published_identities.get(
                        target, staged[targets.index(target)]
                    )
                    quarantine = staging_directory / (target.name + ".rollback")
                    if not self._quarantine_published_generation(
                        target, identity, quarantine
                    ):
                        remove_staging = False
            raise
        finally:
            if remove_staging:
                shutil.rmtree(staging_directory, ignore_errors=True)
        return csv_path, json_path


def _validate_split_result_contract(result: SplitResult, source: object) -> None:
    """Reverse-check every public result/receipt claim against its source."""

    if result._authority_factory_token is not _SPLIT_RESULT_FACTORY_TOKEN:
        raise AuthorityStateError(
            "split result was not produced by ParentGroupSplitter"
        )
    for name in (
        "checker_view_official",
        "cif_files_verified",
        "provisional_input",
        "official_split",
    ):
        if type(getattr(result, name)) is not bool:
            raise AuthorityStateError("{} must be an exact boolean".format(name))
    if result.official_split:
        raise AuthorityStateError(
            "official_split requires a separately authenticated assignment manifest; "
            "none is available"
        )
    for name in (
        "parent_method",
        "requested_leakage_guard",
        "leakage_guard",
        "missing_parent",
        "random_state",
        "implementation_package_version",
        "implementation_split_api_version",
    ):
        value = getattr(result, name)
        if not isinstance(value, str) or not value or value != value.strip():
            raise AuthorityStateError("{} must be an exact nonblank string".format(name))
    for name in ("dataset_version", "checker_view", "parent_input_status", "dataset_input_status"):
        value = getattr(result, name)
        if value is not None and (
            not isinstance(value, str)
            or not value
            or value != value.strip()
        ):
            raise AuthorityStateError(
                "{} must be an exact nonblank string or None".format(name)
            )
    if (
        result.parent_method not in SELECTABLE_PARENT_METHODS
        or result.requested_leakage_guard not in LEAKAGE_GUARD_CHOICES
        or result.leakage_guard
        != resolve_leakage_guard(result.requested_leakage_guard, result.parent_method)
        or result.missing_parent not in {"singleton", "exclude"}
    ):
        raise AuthorityStateError("split parent/leakage policy is invalid")
    if result.implementation_package_version != __version__ or (
        result.implementation_split_api_version != SPLIT_API_VERSION
    ):
        raise AuthorityStateError("split implementation version claim is invalid")

    partition_ids = (
        tuple(result.train_ids),
        tuple(result.validation_ids),
        tuple(result.test_ids),
    )
    flat_ids = tuple(value for partition in partition_ids for value in partition)
    if any(
        not isinstance(value, str) or not value or value != value.strip()
        for value in flat_ids
    ) or len(set(flat_ids)) != len(flat_ids):
        raise AuthorityStateError("split partitions contain invalid or duplicate IDs")
    assignment_ids = set(result.assignments)
    if assignment_ids != set(flat_ids):
        raise AuthorityStateError("split assignments differ from partition IDs")
    for split, ids in zip(_SPLIT_NAMES, partition_ids):
        if any(result.assignments.get(structure_id) != split for structure_id in ids):
            raise AuthorityStateError("split assignment label differs from its partition")
    if set(result.exclusions).intersection(assignment_ids):
        raise AuthorityStateError("assigned and excluded structure IDs overlap")
    if set(result.index_by_id) != assignment_ids.union(result.exclusions):
        raise AuthorityStateError("split result does not cover its release universe")
    if any(
        not isinstance(key, str)
        or not key
        or key != key.strip()
        or not isinstance(value, int)
        or isinstance(value, bool)
        or value < 0
        for key, value in result.index_by_id.items()
    ):
        raise AuthorityStateError("split release index mapping is invalid")
    if sorted(result.index_by_id.values()) != list(range(len(result.index_by_id))):
        raise AuthorityStateError("split release indices must be contiguous and unique")
    for ids, indices in zip(
        partition_ids,
        (result.train_indices, result.validation_indices, result.test_indices),
    ):
        if tuple(result.view_index_by_id[value] for value in ids) != tuple(indices):
            raise AuthorityStateError("split view indices differ from partition IDs")
    if result.leakage_audit.get("passed") is not True:
        raise AuthorityStateError("split leakage audit failed")
    if len(result.fractions) != 3 or any(
        type(value) is not float or not math.isfinite(value) or value < 0.0
        for value in result.fractions
    ) or not math.isclose(sum(result.fractions), 1.0, rel_tol=0.0, abs_tol=1e-9):
        raise AuthorityStateError("split fractions are invalid")
    if any(
        not isinstance(value, str) or not value or value != value.strip()
        for value in result.stratify_by
    ):
        raise AuthorityStateError("split strata names are invalid")

    classified = source
    dataset = getattr(classified, "dataset", None)
    if dataset is None:
        raise AuthorityStateError("split result source dataset is unavailable")
    expected_official = _is_authenticated_official_checker_view(classified)
    if result.checker_view_official is not expected_official:
        raise AuthorityStateError(
            "split official checker claim differs from its authenticated source"
        )
    expected_checker_view = getattr(classified, "checker_view", None)
    expected_dataset_version = getattr(dataset, "dataset_version", None)
    if result.checker_view != expected_checker_view or (
        result.dataset_version != expected_dataset_version
    ):
        raise AuthorityStateError("split source identity claims changed")
    expected_hashes = getattr(dataset, "input_hashes", None)
    if expected_hashes is None:
        expected_hashes = {}
    if not isinstance(expected_hashes, Mapping) or dict(result.input_hashes) != dict(
        expected_hashes
    ):
        raise AuthorityStateError("split input hashes differ from its source")
    expected_cif = getattr(dataset, "cif_files_verified", False)
    if type(expected_cif) is not bool or result.cif_files_verified is not expected_cif:
        raise AuthorityStateError("split CIF verification claim differs from its source")
    try:
        (
            expected_parent_status,
            expected_dataset_status,
            expected_provisional,
        ) = _release_receipt_state(dataset)
    except (TypeError, ValueError) as exc:
        raise AuthorityStateError(
            "split source release-state contract is invalid"
        ) from exc
    if result.parent_input_status != expected_parent_status or (
        result.dataset_input_status != expected_dataset_status
    ):
        raise AuthorityStateError("split release-status claims differ from source")
    if result.provisional_input is not expected_provisional:
        raise AuthorityStateError("split provisional-input claim is invalid")
    target_receipt = getattr(dataset, "target_input_receipt", None)
    if _jsonable(result.target_data) != _jsonable(target_receipt):
        raise AuthorityStateError("split target receipt differs from source")
    if target_receipt is not None:
        # TargetMergedDataset receipt() revalidates its own factory generation.
        receipt_method = getattr(dataset, "receipt", None)
        if not callable(receipt_method) or _jsonable(receipt_method()) != _jsonable(
            target_receipt
        ):
            raise AuthorityStateError("split target generation is unauthenticated")
    expected_hash_closure = _implementation_hashes(
        include_targets=target_receipt is not None,
        target_receipt=target_receipt,
    )
    if dict(result.implementation_hashes) != dict(expected_hash_closure):
        raise AuthorityStateError("split implementation hash closure is invalid")


def _register_split_result(result: SplitResult, source: object) -> SplitResult:
    if result._authority_factory_token is not _SPLIT_RESULT_FACTORY_TOKEN:
        raise AuthorityStateError("cannot register a non-factory split result")
    object.__setattr__(result, "_authority_source", source)
    _validate_split_result_contract(result, source)
    official = _validate_classified_generation_if_present(source)
    source_fingerprint = None
    if not official:
        source_fingerprint = state_fingerprint(_split_source_payload(source))
    fingerprint = _split_result_seal_fingerprint(result)
    _SPLIT_RESULTS.register(
        result,
        fingerprint,
        {
            "source": source,
            "source_official": official,
            "source_fingerprint": source_fingerprint,
            "identity_snapshot": _split_result_identity_snapshot(result),
        },
    )
    return result


def _require_split_result(result: SplitResult) -> None:
    current = _SPLIT_RESULTS.entry(result)
    if current is None:
        raise AuthorityStateError(
            "split result was not produced by the current ParentGroupSplitter factory"
        )
    expected_fingerprint, context = current
    if not isinstance(context, Mapping):  # pragma: no cover - private invariant
        raise AuthorityStateError("split result authority context is malformed")
    if _split_result_identity_snapshot(result) != context.get("identity_snapshot"):
        raise AuthorityStateError(
            "split result changed after construction; recompute the split"
        )
    if _split_result_seal_fingerprint(result) != expected_fingerprint:
        raise AuthorityStateError(
            "split result fingerprint changed after construction; recompute the split"
        )
    source = context.get("source")
    if context.get("source_official") is True:
        if not _validate_classified_generation_if_present(source):
            raise AuthorityStateError("official split source lost authentication")
    else:
        try:
            source_fingerprint = state_fingerprint(_split_source_payload(source))
        except (AttributeError, KeyError, TypeError, ValueError) as exc:
            raise AuthorityStateError(
                "generic split source changed after result construction"
            ) from exc
        if source_fingerprint != context.get("source_fingerprint"):
            raise AuthorityStateError(
                "generic split source changed after result construction"
            )
    _validate_split_result_contract(result, source)


class ParentGroupSplitter:
    """Split a classified CoRE-MOF dataset while keeping parents together.

    ``priority_main`` is the project-defined, conflict-aware explanatory
    hierarchy RAC5, then MOFid v2, then MOFid v1; it is not a row-wise
    first-nonmissing fallback.  ``leakage_guard="auto"`` resolves to the
    broader full-release ``main_union`` graph for ``priority_main`` and to
    ``parent_only`` for every other explanatory method.  See
    :func:`CoREMOF.parents.parent_method_definition` and
    :func:`CoREMOF.parents.leakage_guard_definition` for exact contracts.
    """

    def __init__(
        self,
        classified_dataset,
        parent_method: str = "priority_main",
        leakage_guard: str = "auto",
        missing_parent: str = "singleton",
        random_state=42,
        stratify_by: Sequence[str] = ("label",),
        labels: Optional[Sequence[str]] = ("CR", "NCR"),
        sources: Optional[Sequence[str]] = None,
        variants: Optional[Sequence[str]] = None,
        metals: Optional[Sequence[str]] = None,
        structure_ids: Optional[Sequence[str]] = None,
        required_targets: Optional[Sequence[str]] = None,
        required_target_mode: str = "all",
        official: bool = False,
    ):
        self.classified_dataset = classified_dataset
        _validate_classified_generation_if_present(classified_dataset)
        self.dataset = getattr(classified_dataset, "dataset", None)
        if self.dataset is None and hasattr(classified_dataset, "metadata_rows"):
            self.dataset = classified_dataset
        if self.dataset is None:
            raise TypeError("classified_dataset must expose its source as .dataset")
        _validate_dataset_generation_if_present(self.dataset)
        if not isinstance(parent_method, str):
            raise TypeError("parent method must be a string")
        if parent_method not in SELECTABLE_PARENT_METHODS:
            raise ValueError(
                "Unknown parent method %r; choose one of %s"
                % (parent_method, ", ".join(SELECTABLE_PARENT_METHODS))
            )
        if not isinstance(leakage_guard, str):
            raise TypeError("leakage guard must be a string")
        if leakage_guard not in LEAKAGE_GUARD_CHOICES:
            raise ValueError(
                "leakage_guard must be 'auto', 'parent_only', or 'main_union'"
            )
        requested_leakage_guard = leakage_guard
        leakage_guard = resolve_leakage_guard(leakage_guard, parent_method)
        if missing_parent not in ("exclude", "singleton"):
            raise ValueError("missing_parent must be 'exclude' or 'singleton'")
        if type(official) is not bool:
            raise TypeError("official must be a boolean")
        if official:
            raise OfficialSplitUnavailableError(
                "No audited official split manifest is available for this release. "
                "Generate an exploratory split with official=False; do not label it "
                "as an official CoRE-MOF benchmark."
            )

        self.parent_method = parent_method
        self.requested_leakage_guard = requested_leakage_guard
        self.leakage_guard = leakage_guard
        self.missing_parent = missing_parent
        if isinstance(random_state, bool) or not isinstance(random_state, (int, str)):
            raise TypeError("random_state must be a non-boolean integer or string")
        if isinstance(random_state, str) and (
            not random_state or random_state != random_state.strip()
        ):
            raise ValueError("random_state string must be exact and non-empty")
        self.random_state = str(random_state)
        if isinstance(stratify_by, str):
            self.stratify_by = (stratify_by,)
        else:
            if not isinstance(stratify_by, (list, tuple)):
                raise TypeError(
                    "stratify_by must be a string or iterable of strings; only "
                    "an exact string or ordered list/tuple is accepted"
                )
            self.stratify_by = tuple(stratify_by)
        if any(
            not isinstance(value, str)
            or not value
            or value != value.strip()
            for value in self.stratify_by
        ):
            raise ValueError("stratify_by values must be exact nonblank strings")
        self.label_filter = _as_optional_tuple(labels, "labels")
        if self.label_filter is not None:
            invalid_labels = {
                value.upper() for value in self.label_filter
            }.difference(LABELS)
            if invalid_labels:
                raise ValueError(
                    "unknown split label(s): %s" % ", ".join(sorted(invalid_labels))
                )
        self.source_filter = _as_optional_tuple(sources, "sources")
        self.variant_filter = _as_optional_tuple(variants, "variants")
        self.metal_filter = _as_optional_tuple(metals, "metals")
        self.structure_id_filter = _as_optional_tuple(structure_ids, "structure_ids")
        self._structure_id_filter_set = (
            set(self.structure_id_filter)
            if self.structure_id_filter is not None
            else None
        )
        self.required_targets = _as_optional_tuple(required_targets, "required_targets")
        if self.required_targets is not None:
            if not self.required_targets or any(not value for value in self.required_targets):
                raise ValueError("required_targets must contain non-empty target names")
            if len(set(self.required_targets)) != len(self.required_targets):
                raise ValueError("required_targets contains duplicates")
        if required_target_mode not in ("all", "any"):
            raise ValueError("required_target_mode must be 'all' or 'any'")
        self.required_target_mode = required_target_mode

        self._rows = tuple(getattr(self.dataset, "metadata_rows"))
        self._row_by_id: Dict[str, Mapping[str, object]] = {}
        self._index_by_id: Dict[str, int] = {}
        for index, row in enumerate(self._rows):
            if not isinstance(row, Mapping):
                raise TypeError("Every metadata row must be a mapping")
            raw_structure_id = row.get("structure_id")
            if (
                not isinstance(raw_structure_id, str)
                or not raw_structure_id
                or raw_structure_id != raw_structure_id.strip()
            ):
                raise ValueError(
                    "Every metadata row must have an exact nonblank string structure_id"
                )
            structure_id = raw_structure_id
            if structure_id in self._row_by_id:
                raise ValueError("Duplicate structure_id in metadata_rows: %s" % structure_id)
            self._row_by_id[structure_id] = row
            self._index_by_id[structure_id] = index
        attached_targets = _column_name_tuple(
            getattr(self.dataset, "target_columns", None), "target_columns"
        )
        self._attached_targets = attached_targets
        if self.required_targets is not None:
            if not attached_targets:
                raise ValueError(
                    "required_targets needs a dataset created by merge_targets()"
                )
            unknown_targets = set(self.required_targets).difference(attached_targets)
            if unknown_targets:
                raise ValueError(
                    "unknown required target(s): %s"
                    % ", ".join(sorted(unknown_targets))
                )
        label_values = getattr(classified_dataset, "label_by_id", None)
        if label_values is None:
            label_values = getattr(classified_dataset, "labels", {})
        self._labels = self._coerce_labels(label_values)
        selected_values = getattr(classified_dataset, "structure_ids", None)
        if selected_values is None:
            selected_values = tuple(self._labels)
        self._view_ids = tuple(selected_values)
        if any(
            not isinstance(value, str)
            or not value
            or value != value.strip()
            for value in self._view_ids
        ):
            raise ValueError(
                "classified_dataset structure IDs must be exact nonblank strings"
            )
        if len(set(self._view_ids)) != len(self._view_ids):
            raise ValueError("classified_dataset contains duplicate structure IDs")
        self._view_index_by_id = {
            structure_id: index
            for index, structure_id in enumerate(self._view_ids)
        }
        self._preselected_ids = set(self._view_ids)
        unknown_selected = self._preselected_ids.difference(self._row_by_id)
        if unknown_selected:
            raise ValueError(
                "classified_dataset contains IDs absent from its release: %s"
                % ", ".join(sorted(unknown_selected))
            )
        if self._structure_id_filter_set is not None:
            unknown_filter = self._structure_id_filter_set.difference(self._row_by_id)
            if unknown_filter:
                raise KeyError(
                    "Unknown structure_id filter value(s): %s"
                    % ", ".join(sorted(unknown_filter)[:10])
                )
        available_strata = set().union(
            *(set(row) for row in self._rows)
        ) if self._rows else set()
        unknown_strata = [
            key
            for key in self.stratify_by
            if key != "label" and key not in available_strata
        ]
        if unknown_strata:
            raise ValueError(
                "Unknown stratify_by field(s): %s"
                % ", ".join(sorted(set(unknown_strata)))
            )
        selection_filters = getattr(classified_dataset, "selection_filters", None)
        if selection_filters is None:
            selection_filters = {}
        elif not isinstance(selection_filters, Mapping):
            raise TypeError("selection_filters must be a mapping when present")
        self._selection_filters = self._selection_json_safe(selection_filters)
        target_receipt = getattr(self.dataset, "target_input_receipt", None)
        self._target_input_receipt = (
            self._json_safe(target_receipt) if target_receipt is not None else None
        )

    @classmethod
    def _json_safe(cls, value):
        if isinstance(value, Mapping):
            result = {}
            for key, item in value.items():
                if not isinstance(key, str) or not key or key != key.strip():
                    raise TypeError("receipt mapping keys must be exact nonblank strings")
                result[key] = cls._json_safe(item)
            return result
        if isinstance(value, (list, tuple)):
            return [cls._json_safe(item) for item in value]
        if value is None or type(value) in (str, int, bool):
            return value
        if type(value) is float and math.isfinite(value):
            return value
        raise TypeError("receipt values must be finite JSON data without coercion")

    @classmethod
    def _selection_json_safe(cls, value):
        """Copy selection provenance without coercing keys or unknown values."""

        if isinstance(value, Mapping):
            result = {}
            for key, item in value.items():
                if not isinstance(key, str) or not key or key != key.strip():
                    raise TypeError(
                        "selection_filters keys must be exact nonblank strings"
                    )
                result[key] = cls._selection_json_safe(item)
            return result
        if isinstance(value, (list, tuple)):
            return [cls._selection_json_safe(item) for item in value]
        if value is None or isinstance(value, (str, int, bool)):
            return value
        if isinstance(value, float) and math.isfinite(value):
            return value
        raise TypeError("selection_filters must contain only finite JSON values")

    def _coerce_labels(self, labels) -> Dict[str, str]:
        if isinstance(labels, Mapping):
            result = {}
            for structure_id, value in labels.items():
                if (
                    not isinstance(structure_id, str)
                    or not structure_id
                    or structure_id != structure_id.strip()
                ):
                    raise ValueError(
                        "classification label keys must be exact nonblank structure IDs"
                    )
                if isinstance(value, Mapping) and "label" in value:
                    value = value["label"]
                elif hasattr(value, "label"):
                    value = value.label
                if value is not None:
                    if not isinstance(value, str):
                        raise ValueError(
                            "classification labels must be strings or null"
                        )
                    canonical = value.upper()
                    if canonical not in LABELS:
                        raise ValueError(
                            "unknown classification label for {}: {!r}".format(
                                structure_id, value
                            )
                        )
                    result[structure_id] = canonical
            return result
        values = tuple(labels)
        if len(values) != len(self._rows):
            raise ValueError("A label sequence must align one-to-one with metadata_rows")
        result = {}
        for row, label in zip(self._rows, values):
            if label is None:
                continue
            if not isinstance(label, str):
                raise ValueError("classification labels must be strings or null")
            canonical = label.upper()
            if canonical not in LABELS:
                raise ValueError(
                    "unknown classification label for {}: {!r}".format(
                        row["structure_id"], label
                    )
                )
            result[row["structure_id"]] = canonical
        return result

    @staticmethod
    def _validate_fractions(fractions: Sequence[float]) -> Tuple[float, float, float]:
        if not isinstance(fractions, (list, tuple)):
            raise TypeError("fractions must be an ordered list/tuple of numeric values")
        if any(isinstance(value, bool) for value in fractions):
            raise ValueError("fractions must be numeric values, not booleans")
        if any(not isinstance(value, (Real, Decimal)) for value in fractions):
            raise TypeError(
                "fractions must contain finite non-boolean numeric values"
            )
        values = tuple(float(value) for value in fractions)
        if len(values) != 3:
            raise ValueError("fractions must contain train, validation, and test values")
        if any(not math.isfinite(value) or value < 0.0 for value in values):
            raise ValueError("fractions must be finite and non-negative")
        if not math.isclose(sum(values), 1.0, rel_tol=0.0, abs_tol=1e-9):
            raise ValueError("fractions must sum to 1.0")
        if not any(value > 0.0 for value in values):
            raise ValueError("at least one fraction must be positive")
        return values  # type: ignore[return-value]

    def train_valid_test_split(
        self, fractions: Sequence[float] = (0.8, 0.1, 0.1)
    ) -> SplitResult:
        """Return a deterministic parent-safe three-way split."""

        fractions_tuple = self._validate_fractions(fractions)
        resolver = ParentResolver(self.dataset, missing_parent=self.missing_parent)
        primary = resolver.resolve(self.parent_method)
        if self.leakage_guard == "main_union":
            broad = resolver.resolve("main_union", missing_parent="singleton")
            # ``main_union`` was already computed before any experiment
            # filtering and includes bridges through rows that may later be
            # excluded.  Rebuilding it on the eligible subset would leak.
            block_by_id = dict(broad.group_by_id)
        else:
            block_by_id = dict(primary.group_by_id)

        exclusions: Dict[str, str] = {}
        parent_diagnostics: Dict[str, str] = {}
        candidates: List[str] = []
        for structure_id in self._row_by_id:
            filter_reason = self._filter_reason(structure_id)
            if filter_reason is not None:
                exclusions[structure_id] = filter_reason
            elif structure_id in primary.exclusions:
                parent_reason = primary.exclusions[structure_id]
                if (
                    parent_reason == "PARENT_METHOD_CONFLICT"
                    and self.leakage_guard == "main_union"
                    and structure_id in block_by_id
                ):
                    # The explanatory hierarchy deliberately leaves this row
                    # unresolved, but the broader graph still gives it a safe,
                    # indivisible split block.  Keep the data and retain the
                    # scientific conflict as a diagnostic.
                    parent_diagnostics[structure_id] = parent_reason
                    candidates.append(structure_id)
                else:
                    exclusions[structure_id] = parent_reason
            elif structure_id not in block_by_id:
                exclusions[structure_id] = "MISSING_LEAKAGE_GROUP"
            else:
                candidates.append(structure_id)

        stratum_by_id = {
            structure_id: self._stratum(structure_id) for structure_id in candidates
        }
        assignment_by_block = self._assign_blocks(
            candidates, block_by_id, stratum_by_id, fractions_tuple
        )
        assignments = {
            structure_id: assignment_by_block[block_by_id[structure_id]]
            for structure_id in candidates
        }
        splits_by_block: Dict[str, Set[str]] = {}
        for structure_id, split in assignments.items():
            splits_by_block.setdefault(block_by_id[structure_id], set()).add(split)
        if any(len(splits) > 1 for splits in splits_by_block.values()):
            raise RuntimeError("Internal error: a leakage block crossed split partitions")
        ids_by_split = {
            split: tuple(
                structure_id
                for structure_id in self._row_by_id
                if assignments.get(structure_id) == split
            )
            for split in _SPLIT_NAMES
        }

        dataset_version = getattr(self.dataset, "dataset_version", None)
        if dataset_version is not None and (
            not isinstance(dataset_version, str)
            or not dataset_version
            or dataset_version != dataset_version.strip()
        ):
            raise TypeError(
                "dataset_version must be an exact nonblank string or None"
            )
        checker_view = getattr(self.classified_dataset, "checker_view", None)
        if checker_view is not None and (
            not isinstance(checker_view, str)
            or not checker_view
            or checker_view != checker_view.strip()
        ):
            raise TypeError("checker_view must be an exact nonblank string or None")
        marker = object()
        checker_view_official_value = getattr(
            self.classified_dataset, "checker_view_official", marker
        )
        if checker_view_official_value is marker:
            checker_view_official = False
        elif type(checker_view_official_value) is not bool:
            raise TypeError("checker_view_official must be a boolean when present")
        else:
            checker_view_official = checker_view_official_value
        if checker_view_official and checker_view not in CHECKER_PRESETS:
            raise ValueError(
                "checker_view_official=True requires a canonical official "
                "checker preset"
            )
        if checker_view_official and not _is_authenticated_official_checker_view(
            self.classified_dataset
        ):
            raise ValueError(
                "checker_view_official=True requires an internally authenticated "
                "recomputed release view"
            )

        input_hashes_value = getattr(self.dataset, "input_hashes", None)
        if input_hashes_value is None:
            input_hashes = {}
        elif not isinstance(input_hashes_value, Mapping):
            raise TypeError("input_hashes must be a mapping when present")
        else:
            input_hashes = input_hashes_value
        if any(
            not isinstance(key, str)
            or not key
            or key != key.strip()
            or not isinstance(value, str)
            or not value
            or value != value.strip()
            for key, value in input_hashes.items()
        ):
            raise TypeError(
                "input_hashes keys and values must be exact nonblank strings"
            )

        cif_verified_value = getattr(self.dataset, "cif_files_verified", marker)
        if cif_verified_value is marker:
            cif_files_verified = False
        elif type(cif_verified_value) is not bool:
            raise TypeError("cif_files_verified must be a boolean when present")
        else:
            cif_files_verified = cif_verified_value

        (
            parent_input_status,
            dataset_input_status,
            provisional_input,
        ) = _release_receipt_state(self.dataset)
        filters = {
            "preselection": {
                "active": self._preselected_ids != set(self._row_by_id),
                "selected_count": len(self._preselected_ids),
                "release_count": len(self._row_by_id),
                "structure_ids_sha256": hashlib.sha256(
                    "\n".join(sorted(self._preselected_ids)).encode("utf-8")
                ).hexdigest(),
                "selection_filters": self._selection_filters,
            },
            "labels": list(self.label_filter) if self.label_filter is not None else None,
            "sources": list(self.source_filter) if self.source_filter is not None else None,
            "variants": list(self.variant_filter) if self.variant_filter is not None else None,
            "metals": list(self.metal_filter) if self.metal_filter is not None else None,
            "structure_ids": (
                list(self.structure_id_filter)
                if self.structure_id_filter is not None
                else None
            ),
        }
        if self._target_input_receipt is not None or self.required_targets is not None:
            required = tuple(self.required_targets or ())
            available_counts = {
                target: sum(
                    self._target_is_available(row.get(target)) for row in self._rows
                )
                for target in required
            }
            eligible_count = sum(
                self._passes_target_requirement(row) for row in self._rows
            ) if required else len(self._rows)
            filters["targets"] = {
                "attached_columns": list(self._attached_targets),
                "required": list(required),
                "mode": self.required_target_mode,
                "available_counts": available_counts,
                "eligible_release_count": eligible_count,
                "excluded_release_count": len(self._rows) - eligible_count,
                "filter_precedes_assignment": True,
                "leakage_blocks_use_full_release_universe": True,
            }
        result = SplitResult(
            train_ids=ids_by_split["train"],
            validation_ids=ids_by_split["validation"],
            test_ids=ids_by_split["test"],
            train_indices=tuple(
                self._view_index_by_id[value] for value in ids_by_split["train"]
            ),
            validation_indices=tuple(
                self._view_index_by_id[value] for value in ids_by_split["validation"]
            ),
            test_indices=tuple(
                self._view_index_by_id[value] for value in ids_by_split["test"]
            ),
            assignments=assignments,
            exclusions={
                structure_id: exclusions[structure_id]
                for structure_id in self._row_by_id
                if structure_id in exclusions
            },
            labels={
                structure_id: self._labels[structure_id]
                for structure_id in self._row_by_id
                if structure_id in self._labels
            },
            parent_groups=dict(primary.group_by_id),
            parent_diagnostics=parent_diagnostics,
            parent_conflicts=tuple(
                {
                    "lower_method": conflict.lower_method,
                    "lower_group": conflict.lower_group,
                    "stronger_components": conflict.stronger_components,
                    "member_ids": conflict.member_ids,
                    "unresolved_ids": conflict.unresolved_ids,
                    "component_members": {
                        component: members
                        for component, members in conflict.component_members
                    },
                }
                for conflict in primary.conflicts
            ),
            leakage_groups=block_by_id,
            index_by_id=dict(self._index_by_id),
            view_index_by_id=dict(self._view_index_by_id),
            dataset_version=dataset_version,
            checker_view=checker_view,
            checker_view_official=checker_view_official,
            parent_method=self.parent_method,
            requested_leakage_guard=self.requested_leakage_guard,
            leakage_guard=self.leakage_guard,
            missing_parent=self.missing_parent,
            fractions=fractions_tuple,
            random_state=self.random_state,
            stratify_by=self.stratify_by,
            filters=filters,
            input_hashes=dict(input_hashes),
            implementation_package_version=__version__,
            implementation_split_api_version=SPLIT_API_VERSION,
            implementation_hashes=dict(
                _implementation_hashes(
                    include_targets=self._target_input_receipt is not None,
                    target_receipt=self._target_input_receipt,
                )
            ),
            cif_files_verified=cif_files_verified,
            parent_input_status=parent_input_status,
            dataset_input_status=dataset_input_status,
            provisional_input=provisional_input,
            target_data=self._target_input_receipt,
            official_split=False,
            _authority_factory_token=_SPLIT_RESULT_FACTORY_TOKEN,
        )
        return _register_split_result(result, self.classified_dataset)

    def _filter_reason(self, structure_id: str) -> Optional[str]:
        row = self._row_by_id[structure_id]
        if structure_id not in self._preselected_ids:
            return "PRESELECTION_FILTER"
        if (
            self._structure_id_filter_set is not None
            and structure_id not in self._structure_id_filter_set
        ):
            return "STRUCTURE_ID_FILTER"
        if self.required_targets is not None and not self._passes_target_requirement(row):
            return "MISSING_REQUIRED_TARGET"
        label = self._labels.get(structure_id)
        if label is None:
            return "LABEL_NOT_AVAILABLE"
        if self.label_filter is not None:
            allowed = {value.upper() for value in self.label_filter}
            if label.upper() not in allowed:
                return "LABEL_FILTER"
        if self.source_filter is not None:
            allowed = {value.upper() for value in self.source_filter}
            source = row.get("source_database", "")
            if not isinstance(source, str) or source != source.strip():
                raise TypeError("source_database must be an exact string")
            if source.upper() not in allowed:
                return "SOURCE_FILTER"
        if self.variant_filter is not None:
            allowed = {value.upper() for value in self.variant_filter}
            variant = row.get("structure_variant", "")
            if not isinstance(variant, str) or variant != variant.strip():
                raise TypeError("structure_variant must be an exact string")
            if variant.upper() not in allowed:
                return "VARIANT_FILTER"
        if self.metal_filter is not None:
            requested = {value.casefold() for value in self.metal_filter}
            if not requested.intersection(self._metal_elements(row)):
                return "METAL_FILTER"
        return None

    @staticmethod
    def _target_is_available(value: object) -> bool:
        return value is not None

    def _passes_target_requirement(self, row: Mapping[str, object]) -> bool:
        required = self.required_targets
        if required is None:
            return True
        available = [self._target_is_available(row.get(target)) for target in required]
        if self.required_target_mode == "all":
            return all(available)
        return any(available)

    @staticmethod
    def _metal_elements(row: Mapping[str, object]) -> Set[str]:
        value = row.get("metal_elements")
        if value is None and isinstance(row.get("metals"), Mapping):
            value = row["metals"].get("elements")  # type: ignore[index]
        if value is None:
            return set()
        if isinstance(value, str):
            if value != value.strip():
                raise TypeError("metal_elements must be an exact string")
            values = tuple(item for item in re.split(r"[,;|\s]+", value) if item)
        elif isinstance(value, (list, tuple)):
            values = tuple(value)
            if any(
                not isinstance(item, str)
                or not item
                or item != item.strip()
                for item in values
            ):
                raise TypeError(
                    "metal_elements sequences must contain exact nonblank strings"
                )
        else:
            raise TypeError("metal_elements must be a string or string sequence")
        return {item.casefold() for item in values}

    def _stratum(self, structure_id: str) -> Tuple[str, ...]:
        row = self._row_by_id[structure_id]
        values: List[str] = []
        for key in self.stratify_by:
            if key == "label":
                value = self._labels.get(structure_id)
            else:
                value = row.get(key)
            if value is None or value == "":
                values.append("<MISSING>")
            elif not isinstance(value, str) or value != value.strip():
                raise TypeError(
                    "stratification values must be exact strings or missing"
                )
            else:
                values.append(value)
        return tuple(values)

    def _assign_blocks(
        self,
        structure_ids: Sequence[str],
        block_by_id: Mapping[str, str],
        stratum_by_id: Mapping[str, Tuple[str, ...]],
        fractions: Tuple[float, float, float],
    ) -> Dict[str, str]:
        members_by_block: Dict[str, List[str]] = {}
        for structure_id in structure_ids:
            members_by_block.setdefault(block_by_id[structure_id], []).append(structure_id)
        if not members_by_block:
            return {}

        total = len(structure_ids)
        active_splits = [
            split for split, fraction in zip(_SPLIT_NAMES, fractions) if fraction > 0.0
        ]
        fraction_by_split = dict(zip(_SPLIT_NAMES, fractions))
        target_total = {
            split: fraction_by_split[split] * total for split in active_splits
        }
        stratum_total = Counter(stratum_by_id[value] for value in structure_ids)
        strata = tuple(sorted(stratum_total))
        target_stratum = {
            split: {
                stratum: fraction_by_split[split] * count
                for stratum, count in stratum_total.items()
            }
            for split in active_splits
        }
        current_total = {split: 0 for split in active_splits}
        current_stratum = {split: Counter() for split in active_splits}

        block_strata: Dict[str, Counter] = {}
        for block, members in members_by_block.items():
            block_strata[block] = Counter(
                stratum_by_id[structure_id] for structure_id in members
            )

        blocks = sorted(
            members_by_block,
            key=lambda block: (
                -len(members_by_block[block]),
                -max(block_strata[block].values()),
                _stable_digest(self.random_state, "block-order", block),
            ),
        )
        assignment: Dict[str, str] = {}
        for block_index, block in enumerate(blocks):
            unfilled = [split for split in active_splits if current_total[split] == 0]
            remaining = len(blocks) - block_index
            candidates = active_splits
            if unfilled and remaining <= len(unfilled):
                candidates = unfilled
            scored = []
            for split in candidates:
                score = self._projected_loss(
                    split,
                    len(members_by_block[block]),
                    block_strata[block],
                    current_total,
                    current_stratum,
                    target_total,
                    target_stratum,
                    len(strata),
                )
                scored.append(
                    (
                        score,
                        _stable_digest(self.random_state, "split-tie", block, split),
                        split,
                    )
                )
            chosen = min(scored)[2]
            assignment[block] = chosen
            current_total[chosen] += len(members_by_block[block])
            for stratum, count in block_strata[block].items():
                current_stratum[chosen][stratum] += count
        return assignment

    @staticmethod
    def _projected_loss(
        candidate: str,
        block_size: int,
        block_strata: Mapping[Tuple[str, ...], int],
        current_total: Mapping[str, int],
        current_stratum: Mapping[str, Mapping[Tuple[str, ...], int]],
        target_total: Mapping[str, float],
        target_stratum: Mapping[str, Mapping[Tuple[str, ...], float]],
        stratum_count: int,
    ) -> float:
        # The loss of every untouched split/stratum is identical for all
        # candidate placements.  Compare only the sparse delta introduced by
        # this block, keeping high-cardinality metadata strata linear in N.
        old_total = current_total[candidate]
        new_total = old_total + block_size
        total_target = target_total[candidate]
        total_scale = max(total_target, 1.0)
        total_delta = (
            ((new_total - total_target) / total_scale) ** 2
            - ((old_total - total_target) / total_scale) ** 2
        )
        overflow_delta = (
            (max(0.0, new_total - total_target) / total_scale) ** 2
            - (max(0.0, old_total - total_target) / total_scale) ** 2
        )
        stratum_delta = 0.0
        for stratum, addition in sorted(block_strata.items()):
            old_value = current_stratum[candidate].get(stratum, 0)
            new_value = old_value + addition
            target = target_stratum[candidate][stratum]
            scale = max(target, 1.0)
            stratum_delta += (
                ((new_value - target) / scale) ** 2
                - ((old_value - target) / scale) ** 2
            )
        if stratum_count:
            stratum_delta /= stratum_count
        return total_delta + stratum_delta + (2.0 * overflow_delta)


# Concise alias used in examples and compatible with the MOFDescribe naming style.
ParentSplitter = ParentGroupSplitter


def split_release(
    release_root,
    checkers=None,
    fractions: Sequence[float] = (0.8, 0.1, 0.1),
    verify_cif_files: bool = False,
    **splitter_options
) -> SplitResult:
    """Load, classify, and split a release in one call.

    ``release_root`` is normally the extracted release directory.  An already
    loaded :class:`~CoREMOF.dataset.CoREMOFDataset` or an already classified
    dataset is also accepted, which avoids redundant I/O in notebooks.
    Use :func:`CoREMOF.parents.parent_method_definition` and
    :func:`CoREMOF.parents.leakage_guard_definition` to inspect the exact
    project-defined parent and leakage semantics recorded in the receipt.
    """

    if type(verify_cif_files) is not bool:
        raise TypeError("verify_cif_files must be a boolean")
    requested_official = splitter_options.get("official", False)
    if type(requested_official) is not bool:
        raise TypeError("official must be a boolean")
    if requested_official:
        raise OfficialSplitUnavailableError(
            "No audited official split manifest is available for this release. "
            "Generate an exploratory split with official=False; do not label it "
            "as an official CoRE-MOF benchmark."
        )

    if hasattr(release_root, "label_by_id"):
        if verify_cif_files:
            verified = getattr(
                getattr(release_root, "dataset", None),
                "cif_files_verified",
                False,
            )
            if type(verified) is not bool:
                raise TypeError("cif_files_verified must be a boolean when present")
            if not verified:
                raise ValueError(
                    "verify_cif_files=True requires a release path or a dataset "
                    "previously loaded with verify_cif_files=True"
                )
        if checkers is not None:
            requested_view, _, _ = resolve_checker_view(checkers)
            existing_view = getattr(release_root, "checker_view", None)
            if existing_view is None:
                raise ValueError(
                    "an explicitly requested checker view cannot be verified on "
                    "this preclassified dataset"
                )
            if (
                not isinstance(existing_view, str)
                or not existing_view
                or existing_view != existing_view.strip()
            ):
                raise TypeError(
                    "preclassified checker_view must be an exact nonblank string"
                )
            if existing_view != requested_view:
                raise ValueError(
                    "preclassified dataset uses {!r}, but checkers={!r} requests {!r}".format(
                        existing_view, checkers, requested_view
                    )
                )
        classified_dataset = release_root
    elif hasattr(release_root, "classify") and hasattr(release_root, "metadata_rows"):
        if verify_cif_files:
            verified = getattr(release_root, "cif_files_verified", False)
            if type(verified) is not bool:
                raise TypeError("cif_files_verified must be a boolean when present")
            if not verified:
                raise ValueError(
                    "the supplied dataset was not loaded with verify_cif_files=True"
                )
        classified_dataset = release_root.classify(
            "5checker" if checkers is None else checkers
        )
    else:
        # Lazy import preserves the lightweight import contract for users who
        # construct ParentGroupSplitter from another dataset implementation.
        from CoREMOF.dataset import CoREMOFDataset

        classified_dataset = CoREMOFDataset.from_release(
            release_root, verify_cif_files=verify_cif_files
        ).classify("5checker" if checkers is None else checkers)
    return ParentGroupSplitter(
        classified_dataset, **splitter_options
    ).train_valid_test_split(fractions=fractions)


__all__ = [
    "ParentGroupSplitter",
    "ParentSplitter",
    "OfficialSplitUnavailableError",
    "SPLIT_API_VERSION",
    "SplitResult",
    "split_release",
]
