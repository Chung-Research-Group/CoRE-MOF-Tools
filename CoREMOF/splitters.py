"""Deterministic, parent-aware train/validation/test splitting.

This module intentionally has no NumPy, pandas, scikit-learn, or chemistry
imports.  It can be used on a login node or in a lightweight downstream ML
environment as long as release metadata have already been loaded.
"""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import tempfile
from types import MappingProxyType
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

from CoREMOF import __version__
from CoREMOF.labels import CHECKER_PRESETS, LABELS, resolve_checker_view
from CoREMOF.parents import PARENT_METHODS, ParentResolver


_SPLIT_NAMES = ("train", "validation", "test")
SPLIT_API_VERSION = "0.1.0"


class OfficialSplitUnavailableError(RuntimeError):
    """Raised when an official split is requested without an audited manifest."""


def _deep_freeze(value):
    if isinstance(value, Mapping):
        return MappingProxyType(
            {str(key): _deep_freeze(item) for key, item in value.items()}
        )
    if isinstance(value, (list, tuple, set, frozenset)):
        return tuple(_deep_freeze(item) for item in value)
    return value


def _jsonable(value):
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_jsonable(item) for item in value]
    return value


def _implementation_hashes() -> Mapping[str, str]:
    package_root = Path(__file__).resolve().parent
    result = {}
    for filename in ("dataset.py", "labels.py", "parents.py", "splitters.py"):
        result[filename] = hashlib.sha256(
            (package_root / filename).read_bytes()
        ).hexdigest()
    return result


def _as_optional_tuple(values) -> Optional[Tuple[str, ...]]:
    if values is None:
        return None
    if isinstance(values, str):
        return (values,)
    return tuple(str(value) for value in values)


def _stable_digest(*values: object) -> str:
    text = "\0".join(str(value) for value in values)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


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
    official_split: bool = False

    def __post_init__(self) -> None:
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
        object.__setattr__(
            self,
            "parent_conflicts",
            tuple(_deep_freeze(conflict) for conflict in self.parent_conflicts),
        )

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
        return {
            "schema_version": "coremof-split-receipt/1.0",
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
            "parent_input_status": self.parent_input_status,
            "dataset_input_status": self.dataset_input_status,
            "provisional_input": self.provisional_input,
            "official_split": self.official_split,
            "parent_method": self.parent_method,
            "leakage_guard": self.leakage_guard,
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

    @staticmethod
    def _atomic_target(path, overwrite: bool) -> Tuple[Path, Path]:
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

    def write(
        self,
        output_directory,
        stem: str = "coremof_split",
        overwrite: bool = False,
    ) -> Tuple[Path, Path]:
        """Write ``<stem>.csv`` and ``<stem>.json`` and return both paths."""

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
        existed = {target: target.exists() for target in targets}
        backups = {}
        published = []
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
                    os.replace(str(source), str(target))
                else:
                    os.link(str(source), str(target))
                published.append(target)
        except BaseException:
            for target in reversed(targets):
                backup = backups.get(target)
                if backup is not None and backup.exists():
                    os.replace(str(backup), str(target))
                elif target in published and not existed[target] and target.exists():
                    source = staged[targets.index(target)]
                    # A concurrent writer must never be removed during
                    # rollback.  Delete only the hard link we published.
                    if source.exists() and os.path.samefile(str(source), str(target)):
                        target.unlink()
            raise
        finally:
            shutil.rmtree(staging_directory, ignore_errors=True)
        return csv_path, json_path


class ParentGroupSplitter:
    """Split a classified CoRE-MOF dataset while keeping parents together.

    ``priority_main`` is the recommended parent method.  It uses RAC5 first,
    then MOFid v2, then MOFid v1.  Its automatic leakage guard is the broader
    ``main_union`` graph.  Explicit direct/reference methods use only their
    own parent groups unless ``leakage_guard="main_union"`` is requested.
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
        official: bool = False,
    ):
        self.classified_dataset = classified_dataset
        self.dataset = getattr(classified_dataset, "dataset", None)
        if self.dataset is None and hasattr(classified_dataset, "metadata_rows"):
            self.dataset = classified_dataset
        if self.dataset is None:
            raise TypeError("classified_dataset must expose its source as .dataset")
        if not isinstance(parent_method, str):
            raise TypeError("parent method must be a string")
        parent_method = parent_method.strip().lower()
        splitter_parent_methods = tuple(
            method for method in PARENT_METHODS if method != "main_union"
        )
        if parent_method not in splitter_parent_methods:
            raise ValueError(
                "Unknown parent method %r; choose one of %s"
                % (parent_method, ", ".join(splitter_parent_methods))
            )
        if leakage_guard not in ("auto", "parent_only", "main_union"):
            raise ValueError("leakage_guard must be 'auto', 'parent_only', or 'main_union'")
        if leakage_guard == "auto":
            leakage_guard = "main_union" if parent_method == "priority_main" else "parent_only"
        if missing_parent not in ("exclude", "singleton"):
            raise ValueError("missing_parent must be 'exclude' or 'singleton'")
        if official:
            raise OfficialSplitUnavailableError(
                "No audited official split manifest is available for this release. "
                "Generate an exploratory split with official=False; do not label it "
                "as an official CoRE-MOF benchmark."
            )

        self.parent_method = parent_method
        self.leakage_guard = leakage_guard
        self.missing_parent = missing_parent
        self.random_state = str(random_state)
        if isinstance(stratify_by, str):
            self.stratify_by = (stratify_by,)
        else:
            self.stratify_by = tuple(str(value) for value in (stratify_by or ()))
        self.label_filter = _as_optional_tuple(labels)
        if self.label_filter is not None:
            invalid_labels = {
                value.upper() for value in self.label_filter
            }.difference(LABELS)
            if invalid_labels:
                raise ValueError(
                    "unknown split label(s): %s" % ", ".join(sorted(invalid_labels))
                )
        self.source_filter = _as_optional_tuple(sources)
        self.variant_filter = _as_optional_tuple(variants)
        self.metal_filter = _as_optional_tuple(metals)
        self.structure_id_filter = _as_optional_tuple(structure_ids)
        self._structure_id_filter_set = (
            set(self.structure_id_filter)
            if self.structure_id_filter is not None
            else None
        )

        self._rows = tuple(getattr(self.dataset, "metadata_rows"))
        self._row_by_id: Dict[str, Mapping[str, object]] = {}
        self._index_by_id: Dict[str, int] = {}
        for index, row in enumerate(self._rows):
            structure_id = str(row.get("structure_id", "")).strip()
            if not structure_id:
                raise ValueError("Every metadata row must have a non-empty structure_id")
            if structure_id in self._row_by_id:
                raise ValueError("Duplicate structure_id in metadata_rows: %s" % structure_id)
            self._row_by_id[structure_id] = row
            self._index_by_id[structure_id] = index
        label_values = getattr(classified_dataset, "label_by_id", None)
        if label_values is None:
            label_values = getattr(classified_dataset, "labels", {})
        self._labels = self._coerce_labels(label_values)
        selected_values = getattr(classified_dataset, "structure_ids", None)
        if selected_values is None:
            selected_values = tuple(self._labels)
        self._view_ids = tuple(str(value) for value in selected_values)
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
        selection_filters = getattr(classified_dataset, "selection_filters", {}) or {}
        self._selection_filters = self._json_safe(selection_filters)

    @classmethod
    def _json_safe(cls, value):
        if isinstance(value, Mapping):
            return {str(key): cls._json_safe(item) for key, item in value.items()}
        if isinstance(value, (set, frozenset)):
            return [cls._json_safe(item) for item in sorted(value, key=lambda item: str(item))]
        if isinstance(value, (list, tuple)):
            return [cls._json_safe(item) for item in value]
        if value is None or isinstance(value, (str, int, float, bool)):
            return value
        return str(value)

    def _coerce_labels(self, labels) -> Dict[str, str]:
        if isinstance(labels, Mapping):
            result = {}
            for structure_id, value in labels.items():
                if isinstance(value, Mapping) and "label" in value:
                    value = value["label"]
                elif hasattr(value, "label"):
                    value = value.label
                if value is not None:
                    canonical = str(value).upper()
                    if canonical not in LABELS:
                        raise ValueError(
                            "unknown classification label for {}: {!r}".format(
                                structure_id, value
                            )
                        )
                    result[str(structure_id)] = canonical
            return result
        values = tuple(labels)
        if len(values) != len(self._rows):
            raise ValueError("A label sequence must align one-to-one with metadata_rows")
        result = {}
        for row, label in zip(self._rows, values):
            if label is None:
                continue
            canonical = str(label).upper()
            if canonical not in LABELS:
                raise ValueError(
                    "unknown classification label for {}: {!r}".format(
                        row["structure_id"], label
                    )
                )
            result[str(row["structure_id"])] = canonical
        return result

    @staticmethod
    def _validate_fractions(fractions: Sequence[float]) -> Tuple[float, float, float]:
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
        checker_view = getattr(self.classified_dataset, "checker_view", None)
        checker_view_official = bool(
            getattr(
                self.classified_dataset,
                "checker_view_official",
                checker_view in CHECKER_PRESETS,
            )
        )
        input_hashes = getattr(self.dataset, "input_hashes", {}) or {}
        cif_files_verified = bool(
            getattr(self.dataset, "cif_files_verified", False)
        )
        parent_methods = getattr(self.dataset, "parent_group_methods", {}) or {}
        if isinstance(parent_methods, Mapping):
            parent_input_status_value = parent_methods.get("release_status")
        else:
            parent_input_status_value = None
        parent_input_status = (
            str(parent_input_status_value)
            if parent_input_status_value is not None
            else None
        )
        dataset_info = getattr(self.dataset, "dataset_info", {}) or {}
        if isinstance(dataset_info, Mapping):
            dataset_input_status_value = dataset_info.get("release_status")
        else:
            dataset_input_status_value = None
        dataset_input_status = (
            str(dataset_input_status_value)
            if dataset_input_status_value is not None
            else None
        )
        # Publication state is deliberately closed.  A label such as
        # FINAL_CANDIDATE or NOT_FINAL must never be guessed to mean final.
        provisional_input = not (
            parent_input_status == "FINAL" and dataset_input_status == "FINAL"
        )
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
        return SplitResult(
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
            dataset_version=str(dataset_version) if dataset_version is not None else None,
            checker_view=str(checker_view) if checker_view is not None else None,
            checker_view_official=checker_view_official,
            parent_method=self.parent_method,
            leakage_guard=self.leakage_guard,
            missing_parent=self.missing_parent,
            fractions=fractions_tuple,
            random_state=self.random_state,
            stratify_by=self.stratify_by,
            filters=filters,
            input_hashes={str(key): str(value) for key, value in input_hashes.items()},
            implementation_package_version=__version__,
            implementation_split_api_version=SPLIT_API_VERSION,
            implementation_hashes=dict(_implementation_hashes()),
            cif_files_verified=cif_files_verified,
            parent_input_status=parent_input_status,
            dataset_input_status=dataset_input_status,
            provisional_input=provisional_input,
            official_split=False,
        )

    def _filter_reason(self, structure_id: str) -> Optional[str]:
        row = self._row_by_id[structure_id]
        if structure_id not in self._preselected_ids:
            return "PRESELECTION_FILTER"
        if (
            self._structure_id_filter_set is not None
            and structure_id not in self._structure_id_filter_set
        ):
            return "STRUCTURE_ID_FILTER"
        label = self._labels.get(structure_id)
        if label is None:
            return "LABEL_NOT_AVAILABLE"
        if self.label_filter is not None:
            allowed = {value.upper() for value in self.label_filter}
            if label.upper() not in allowed:
                return "LABEL_FILTER"
        if self.source_filter is not None:
            allowed = {value.upper() for value in self.source_filter}
            if str(row.get("source_database", "")).upper() not in allowed:
                return "SOURCE_FILTER"
        if self.variant_filter is not None:
            allowed = {value.upper() for value in self.variant_filter}
            if str(row.get("structure_variant", "")).upper() not in allowed:
                return "VARIANT_FILTER"
        if self.metal_filter is not None:
            requested = {value.casefold() for value in self.metal_filter}
            if not requested.intersection(self._metal_elements(row)):
                return "METAL_FILTER"
        return None

    @staticmethod
    def _metal_elements(row: Mapping[str, object]) -> Set[str]:
        value = row.get("metal_elements")
        if value is None and isinstance(row.get("metals"), Mapping):
            value = row["metals"].get("elements")  # type: ignore[index]
        if value is None:
            return set()
        if isinstance(value, str):
            values = re.split(r"[,;|\s]+", value.strip("[](){} "))
        else:
            try:
                values = list(value)  # type: ignore[arg-type]
            except TypeError:
                values = [value]
        return {str(item).strip("'\" ").casefold() for item in values if str(item).strip()}

    def _stratum(self, structure_id: str) -> Tuple[str, ...]:
        row = self._row_by_id[structure_id]
        values: List[str] = []
        for key in self.stratify_by:
            if key == "label":
                value = self._labels.get(structure_id)
            else:
                value = row.get(key)
            if isinstance(value, (dict, list, tuple, set)):
                value = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)
            values.append("<MISSING>" if value is None or value == "" else str(value))
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
    """

    if splitter_options.get("official"):
        raise OfficialSplitUnavailableError(
            "No audited official split manifest is available for this release. "
            "Generate an exploratory split with official=False; do not label it "
            "as an official CoRE-MOF benchmark."
        )

    if hasattr(release_root, "label_by_id"):
        if verify_cif_files and not getattr(
            getattr(release_root, "dataset", None), "cif_files_verified", False
        ):
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
            if str(existing_view) != requested_view:
                raise ValueError(
                    "preclassified dataset uses {!r}, but checkers={!r} requests {!r}".format(
                        existing_view, checkers, requested_view
                    )
                )
        classified_dataset = release_root
    elif hasattr(release_root, "classify") and hasattr(release_root, "metadata_rows"):
        if verify_cif_files and not getattr(release_root, "cif_files_verified", False):
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
