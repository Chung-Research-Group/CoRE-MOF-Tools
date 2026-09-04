"""Deferred target attachment for frozen split and benchmark assignments.

Attachment is a left join by exact current public ``structure_id`` or an
explicit hash-bound alias accepted by :class:`CoREMOF.targets.AliasRegistry`.
``keep`` preserves every frozen selected ID and leaves unavailable targets as
null. ``error`` requires every attached target for every selected ID. ``drop``
creates only a derived filtered view: it never refills a cohort, changes an
existing partition, rebalances, or resplits.  The attachment receipt references
the original assignment digest.  Target input/value hashes remain confined to
the attachment receipt and never modify the frozen split receipt.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field, fields, is_dataclass
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
from types import MappingProxyType
from typing import Dict, Mapping, Optional, Sequence, Tuple

from CoREMOF import __version__
from CoREMOF._authority import (
    AuthorityStateError,
    IdentitySealRegistry,
    state_fingerprint,
)


TARGET_ATTACHMENT_API_VERSION = "0.1.0"


class TargetAttachmentError(ValueError):
    """Raised when targets cannot be attached without changing the assignment."""


_MISSING_POLICY_DEFINITIONS = MappingProxyType(
    {
        "keep": "preserve every frozen selected ID and retain null targets",
        "error": "require every target for every frozen selected ID",
        "drop": "derive a filtered view only; never refill, rebalance, or resplit",
    }
)
_RESERVED_EXPORT_COLUMNS = frozenset({"structure_id", "partition", "run_key"})
_FROZEN_ASSIGNMENT_FACTORY_TOKEN = object()
_ATTACHED_RESULT_FACTORY_TOKEN = object()
_FROZEN_ASSIGNMENTS = IdentitySealRegistry("frozen assignment verification")
_ATTACHED_RESULTS = IdentitySealRegistry("target attachment generation")
_ACCEPTED_ASSIGNMENT_RECEIPT_SCHEMAS = frozenset(
    {
        "coremof-data-split-receipt/1.0",
        "coremof-cr-ncr-benchmark-suite-receipt/1.0",
        "coremof-split-receipt/1.1",
    }
)


def _freeze(value: object) -> object:
    if isinstance(value, Mapping):
        result = {}
        for key, item in value.items():
            if not isinstance(key, str) or not key or key != key.strip():
                raise TypeError("mapping keys must be exact nonblank strings")
            result[key] = _freeze(item)
        return MappingProxyType(result)
    if isinstance(value, (list, tuple, set, frozenset)):
        return tuple(_freeze(item) for item in value)
    return value


def _jsonable(value: object) -> object:
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_jsonable(item) for item in value]
    return value


def _canonical_json(value: object) -> str:
    return json.dumps(
        _jsonable(value),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )


def _digest(value: object) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _seal_projection(value: object) -> object:
    """Convert nested frozen dataclass state to authority-fingerprint values."""

    if is_dataclass(value) and not isinstance(value, type):
        return (
            "dataclass",
            "{}.{}".format(type(value).__module__, type(value).__qualname__),
            tuple(
                (definition.name, _seal_projection(getattr(value, definition.name)))
                for definition in fields(value)
                if definition.name
                not in {
                    "_factory_token",
                    "_authority_factory_token",
                    "_authority_generation_marker",
                }
            ),
        )
    if isinstance(value, Mapping):
        return (
            "mapping",
            tuple(
                (_seal_projection(key), _seal_projection(item))
                for key, item in value.items()
            ),
        )
    if type(value) in (tuple, list):
        return tuple(_seal_projection(item) for item in value)
    return value


def _frozen_assignment_fingerprint(value: object) -> str:
    return state_fingerprint(
        ("coremof-frozen-assignment-verification/1", _seal_projection(value))
    )


def _register_frozen_assignment(value: object):
    if getattr(value, "_factory_token", None) is not _FROZEN_ASSIGNMENT_FACTORY_TOKEN:
        raise TargetAttachmentError("cannot register a non-factory frozen assignment")
    _FROZEN_ASSIGNMENTS.register(
        value,
        _frozen_assignment_fingerprint(value),
        {"source_schema": getattr(value, "_source_schema", None)},
    )
    return value


def _require_frozen_assignment(value: object) -> None:
    try:
        context = _FROZEN_ASSIGNMENTS.require(
            value, _frozen_assignment_fingerprint(value)
        )
    except AuthorityStateError as error:
        raise TargetAttachmentError(str(error)) from error
    if not isinstance(context, Mapping) or context.get("source_schema") != getattr(
        value, "_source_schema", None
    ):
        raise TargetAttachmentError("frozen assignment authority context changed")
    _source_hashes()


def _attached_result_fingerprint(value: object) -> str:
    return state_fingerprint(
        ("coremof-target-attachment-generation/1", _seal_projection(value))
    )


def _register_attached_result(value: object):
    if (
        getattr(value, "_authority_factory_token", None)
        is not _ATTACHED_RESULT_FACTORY_TOKEN
    ):
        raise TargetAttachmentError("cannot register a non-factory target attachment")
    run_views = getattr(value, "run_views", None)
    if isinstance(run_views, Mapping):
        for view in run_views.values():
            _require_attached_result(view)
    _ATTACHED_RESULTS.register(value, _attached_result_fingerprint(value))
    return value


def _require_attached_result(value: object) -> None:
    try:
        _ATTACHED_RESULTS.require(value, _attached_result_fingerprint(value))
    except AuthorityStateError as error:
        raise TargetAttachmentError(str(error)) from error
    run_views = getattr(value, "run_views", None)
    if isinstance(run_views, Mapping):
        for view in run_views.values():
            _require_attached_result(view)
    _source_hashes()


def _sha256(value: object, name: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise TargetAttachmentError(
            "{} must be a lowercase SHA-256 digest".format(name)
        )
    return value


def _assignment_digest(
    assignments: Mapping[str, str], exclusions: Mapping[str, str]
) -> str:
    lines = []
    for structure_id in sorted(set(assignments).union(exclusions)):
        lines.append(
            "\t".join(
                (
                    structure_id,
                    assignments.get(structure_id, "EXCLUDED"),
                    exclusions.get(structure_id, ""),
                )
            )
        )
    return hashlib.sha256("\n".join(lines).encode("utf-8")).hexdigest()


def _release_binding(dataset: object) -> Mapping[str, object]:
    base = getattr(dataset, "base_dataset", None)
    if base is not None and getattr(dataset, "target_columns", ()):
        dataset = base
    structure_ids = tuple(getattr(dataset, "structure_ids", ()))
    input_hashes = getattr(dataset, "input_hashes", None)
    if not structure_ids or not isinstance(input_hashes, Mapping) or not input_hashes:
        raise TargetAttachmentError(
            "target attachment requires a release with structure IDs and input hashes"
        )
    normalized_hashes = {}
    for name, digest in input_hashes.items():
        if not isinstance(name, str) or not name or name != name.strip():
            raise TargetAttachmentError(
                "release input-hash names must be exact strings"
            )
        normalized_hashes[name] = _sha256(digest, "release input hash")
    return MappingProxyType(
        {
            "dataset_version": getattr(dataset, "dataset_version", None),
            "structure_count": len(structure_ids),
            "ordered_structure_ids_sha256": hashlib.sha256(
                "\n".join(structure_ids).encode("utf-8")
            ).hexdigest(),
            "structure_ids_sha256": hashlib.sha256(
                "\n".join(sorted(structure_ids)).encode("utf-8")
            ).hexdigest(),
            "input_hashes": MappingProxyType(dict(sorted(normalized_hashes.items()))),
        }
    )


def _release_binding_matches(
    expected: Mapping[str, object], observed: Mapping[str, object]
) -> bool:
    if expected.get("dataset_version") != observed.get("dataset_version"):
        return False
    if "structure_count" in expected and expected.get(
        "structure_count"
    ) != observed.get("structure_count"):
        return False
    for key in ("ordered_structure_ids_sha256", "structure_ids_sha256"):
        if key in expected and expected.get(key) != observed.get(key):
            return False
    expected_hashes = expected.get("input_hashes")
    observed_hashes = observed.get("input_hashes")
    if not isinstance(expected_hashes, Mapping) or not isinstance(
        observed_hashes, Mapping
    ):
        return False
    # ``release_binding`` is the base-release identity, not the broader split
    # input closure.  Accepting a subset would let a caller omit one frozen
    # release input while retaining the same version and ID universe.
    return dict(expected_hashes) == dict(observed_hashes)


def _receipt_release_binding(receipt: Mapping[str, object]) -> object:
    binding = receipt.get("release_binding")
    if isinstance(binding, Mapping):
        return binding
    cohort_receipt = receipt.get("cohort_receipt")
    if isinstance(cohort_receipt, Mapping):
        binding = cohort_receipt.get("release_binding")
        if isinstance(binding, Mapping):
            return binding
    return None


def _receipt_official_split(receipt: Mapping[str, object]) -> object:
    value = receipt.get("official_split")
    if value is None:
        cohort_receipt = receipt.get("cohort_receipt")
        if isinstance(cohort_receipt, Mapping):
            value = cohort_receipt.get("official_split")
    return value


def _validate_live_assignment_authority(
    assignment: object, observed_binding: Mapping[str, object]
) -> Optional[bool]:
    """Revalidate mutable live objects against their frozen source receipt."""

    if isinstance(assignment, FrozenAssignmentManifest):
        assignment._validate_integrity()
        if not _release_binding_matches(assignment.release_binding, observed_binding):
            raise TargetAttachmentError(
                "loaded release does not match the frozen assignment receipt"
            )
        return assignment.official_split

    receipt_method = getattr(assignment, "receipt", None)
    receipt = receipt_method() if callable(receipt_method) else None
    if not isinstance(receipt, Mapping):
        if isinstance(assignment, (Mapping, list, tuple)):
            return None
        raise TargetAttachmentError(
            "live assignment objects require a known sealed receipt; use a raw "
            "ID list, partition mapping, or verified FrozenAssignmentManifest"
        )

    schema = receipt.get("schema_version")
    # A receipt's text is not authority by itself.  Only the exact package
    # result type paired with its one accepted schema can delegate a live
    # assignment.  The result's own receipt() method performs its identity-
    # seal check before this dispatch is reached.
    from CoREMOF.benchmarks import CRNCRBenchmarkSuite, DataSplitResult
    from CoREMOF.splitters import SplitResult

    schema_by_type = {
        DataSplitResult: "coremof-data-split-receipt/1.0",
        CRNCRBenchmarkSuite: "coremof-cr-ncr-benchmark-suite-receipt/1.0",
        SplitResult: "coremof-split-receipt/1.1",
    }
    expected_schema = schema_by_type.get(type(assignment))
    if expected_schema is None:
        raise TargetAttachmentError(
            "live assignment type is not an accepted sealed package result"
        )
    if schema != expected_schema:
        raise TargetAttachmentError(
            "live assignment receipt schema does not match its sealed result type"
        )

    expected_binding = _receipt_release_binding(receipt)
    if schema in {
        "coremof-data-split-receipt/1.0",
        "coremof-cr-ncr-benchmark-suite-receipt/1.0",
    } and not isinstance(expected_binding, Mapping):
        raise TargetAttachmentError(
            "source assignment receipt is missing its release binding"
        )
    if isinstance(expected_binding, Mapping) and not _release_binding_matches(
        expected_binding, observed_binding
    ):
        raise TargetAttachmentError(
            "live assignment release differs from its frozen receipt binding"
        )

    expected_official = _receipt_official_split(receipt)
    if expected_official is not None and type(expected_official) is not bool:
        raise TargetAttachmentError(
            "source receipt official_split must be boolean or null"
        )
    if expected_official is True:
        raise TargetAttachmentError(
            "official assignment attachment is unavailable without a separately "
            "audited official manifest loader"
        )
    observed_official = getattr(assignment, "official_split", None)
    if observed_official != expected_official:
        raise TargetAttachmentError(
            "live assignment official_split differs from its frozen receipt"
        )

    if schema == "coremof-data-split-receipt/1.0":
        assignments = getattr(assignment, "assignments", None)
        exclusions = getattr(assignment, "exclusions", None)
        if not isinstance(assignments, Mapping) or not isinstance(exclusions, Mapping):
            raise TargetAttachmentError("live split lacks its frozen assignment maps")
        computed = _assignment_digest(
            {str(key): str(value) for key, value in assignments.items()},
            {str(key): str(value) for key, value in exclusions.items()},
        )
        declared = getattr(assignment, "assignment_digest", None)
        expected = receipt.get("assignment_sha256")
        if computed != expected or declared != expected:
            raise TargetAttachmentError(
                "live assignment digest differs from its frozen receipt"
            )
    elif schema == "coremof-cr-ncr-benchmark-suite-receipt/1.0":
        runs = getattr(assignment, "runs", None)
        declared_suite = getattr(assignment, "assignment_digest", None)
        receipt_runs = receipt.get("runs")
        if not isinstance(runs, (list, tuple)) or not isinstance(receipt_runs, list):
            raise TargetAttachmentError(
                "live benchmark suite lacks frozen run evidence"
            )
        expected_by_key = {
            item.get("run_key"): item.get("assignment_sha256")
            for item in receipt_runs
            if isinstance(item, Mapping)
        }
        recomputed = []
        for run in runs:
            run_key = getattr(run, "run_key", None)
            assignments = getattr(run, "assignments", None)
            if not isinstance(run_key, str) or not isinstance(assignments, Mapping):
                raise TargetAttachmentError("live benchmark run is malformed")
            digest = _assignment_digest(
                {str(key): str(value) for key, value in assignments.items()}, {}
            )
            if (
                expected_by_key.get(run_key) != digest
                or getattr(run, "assignment_digest", None) != digest
            ):
                raise TargetAttachmentError(
                    "live benchmark run differs from its frozen receipt"
                )
            recomputed.append((run_key, digest))
        suite_digest = _digest(recomputed)
        if (
            suite_digest != receipt.get("suite_assignment_sha256")
            or declared_suite != suite_digest
            or len(recomputed) != len(expected_by_key)
        ):
            raise TargetAttachmentError(
                "live benchmark suite differs from its frozen receipt"
            )
    return expected_official


@dataclass(frozen=True)
class FrozenAssignmentManifest:
    """Receipt-verified assignment used by the persistence CLI."""

    assignments: Mapping[str, str]
    assignment_digest: str
    release_binding: Mapping[str, object]
    official_split: Optional[bool] = None
    exclusions: Mapping[str, str] = field(default_factory=dict)
    _source_schema: str = field(default="", repr=False, compare=False)
    _factory_token: object = field(default=None, repr=False, compare=False)
    _integrity_digest: str = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if self._factory_token is not _FROZEN_ASSIGNMENT_FACTORY_TOKEN:
            raise TargetAttachmentError(
                "FrozenAssignmentManifest must be created by "
                "frozen_assignment_manifest() from a verified manifest and receipt"
            )
        assignments, _ = _partition_mapping(dict(self.assignments))
        object.__setattr__(self, "assignments", assignments)
        object.__setattr__(
            self,
            "assignment_digest",
            _sha256(self.assignment_digest, "assignment_digest"),
        )
        if not isinstance(self.release_binding, Mapping):
            raise TargetAttachmentError("frozen assignment requires a release binding")
        object.__setattr__(self, "release_binding", _freeze(self.release_binding))
        if not isinstance(self.exclusions, Mapping):
            raise TargetAttachmentError(
                "frozen assignment exclusions must be a mapping"
            )
        exclusions = {}
        for structure_id, reason in self.exclusions.items():
            if (
                not isinstance(structure_id, str)
                or not structure_id
                or structure_id != structure_id.strip()
                or not isinstance(reason, str)
            ):
                raise TargetAttachmentError(
                    "frozen assignment exclusions require exact IDs and string reasons"
                )
            if structure_id in assignments:
                raise TargetAttachmentError(
                    "frozen assignment ID appears as both selected and excluded"
                )
            exclusions[structure_id] = reason
        object.__setattr__(self, "exclusions", _freeze(exclusions))
        if self._source_schema not in _ACCEPTED_ASSIGNMENT_RECEIPT_SCHEMAS:
            raise TargetAttachmentError(
                "frozen assignment receipt schema is not accepted"
            )
        if self.official_split is not None and type(self.official_split) is not bool:
            raise TargetAttachmentError(
                "frozen assignment official_split must be boolean or null"
            )
        if self.official_split is True:
            raise TargetAttachmentError(
                "official assignment attachment is unavailable without a separately "
                "audited official manifest loader"
            )
        if self._source_schema != "coremof-split-receipt/1.1":
            observed = _assignment_digest(assignments, exclusions)
            if observed != self.assignment_digest:
                raise TargetAttachmentError(
                    "frozen assignment digest differs from its selected/excluded rows"
                )
        object.__setattr__(self, "_integrity_digest", self._current_integrity_digest())

    def _current_integrity_digest(self) -> str:
        return _digest(
            {
                "assignments": self.assignments,
                "exclusions": self.exclusions,
                "assignment_digest": self.assignment_digest,
                "release_binding": self.release_binding,
                "official_split": self.official_split,
                "source_schema": self._source_schema,
            }
        )

    def _validate_integrity(self) -> None:
        _require_frozen_assignment(self)
        if self._factory_token is not _FROZEN_ASSIGNMENT_FACTORY_TOKEN:
            raise TargetAttachmentError("frozen assignment factory binding changed")
        if self.official_split is True:
            raise TargetAttachmentError(
                "official assignment attachment is unavailable without a separately "
                "audited official manifest loader"
            )
        if self._current_integrity_digest() != self._integrity_digest:
            raise TargetAttachmentError(
                "frozen assignment changed after receipt verification"
            )


def _current_source_hashes() -> Dict[str, str]:
    root = Path(__file__).resolve().parent
    names = (
        "_authority.py",
        "attachments.py",
        "dataset.py",
        "targets.py",
        "_transactions.py",
    )
    return {
        name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in names
    }


_IMPORTED_SOURCE_HASHES = MappingProxyType(_current_source_hashes())


def _source_hashes() -> Mapping[str, str]:
    """Return import-bound source hashes and reject an on-disk code drift."""

    current = _current_source_hashes()
    if current != dict(_IMPORTED_SOURCE_HASHES):
        changed = sorted(
            name
            for name in set(current).union(_IMPORTED_SOURCE_HASHES)
            if current.get(name) != _IMPORTED_SOURCE_HASHES.get(name)
        )
        raise TargetAttachmentError(
            "CoREMOF implementation source changed after module import: {}".format(
                ", ".join(changed)
            )
        )
    return _IMPORTED_SOURCE_HASHES


def _exact_ids(values: object, name: str = "structure IDs") -> Tuple[str, ...]:
    if isinstance(values, str) or not isinstance(values, (list, tuple)):
        raise TypeError("{} must be an ordered list/tuple of strings".format(name))
    result = tuple(values)
    if any(
        not isinstance(value, str) or not value or value != value.strip()
        for value in result
    ):
        raise TargetAttachmentError(
            "{} must contain exact nonblank strings".format(name)
        )
    if len(set(result)) != len(result):
        raise TargetAttachmentError("{} contains duplicates".format(name))
    return result


def _dataset_from_assignment(value: object, explicit: object = None) -> object:
    embedded = None
    dataset = getattr(value, "_dataset", None)
    if dataset is not None:
        embedded = dataset
    if embedded is None:
        source = getattr(value, "_authority_source", None)
        if source is not None:
            dataset = getattr(source, "dataset", None)
            if dataset is not None:
                embedded = dataset
    if embedded is None:
        source = getattr(value, "dataset", None)
        if source is not None and hasattr(source, "metadata_rows"):
            embedded = source
    if explicit is not None:
        explicit_dataset = getattr(explicit, "dataset", explicit)
        frozen_binding = getattr(value, "release_binding", None)
        if isinstance(frozen_binding, Mapping) and not _release_binding_matches(
            frozen_binding, _release_binding(explicit_dataset)
        ):
            raise TargetAttachmentError(
                "loaded release does not match the frozen assignment receipt"
            )
        if embedded is not None and _canonical_json(
            _release_binding(embedded)
        ) != _canonical_json(_release_binding(explicit_dataset)):
            raise TargetAttachmentError(
                "dataset= does not match the frozen assignment release binding"
            )
        return explicit_dataset
    if embedded is not None:
        return embedded
    raise TargetAttachmentError(
        "dataset= is required for raw ID lists and partition mappings"
    )


def _config_like(path: Path) -> bool:
    if path.suffix.casefold() != ".json":
        return False
    try:
        with path.open("r", encoding="utf-8") as handle:
            value = json.load(handle)
    except (OSError, UnicodeError, json.JSONDecodeError):
        return False
    return isinstance(value, Mapping) and "sources" in value


def _merged_targets(
    dataset: object,
    target_data: object,
    alias_registry: object = None,
    feature_tables: Sequence[str] = (),
    verify_cif_files: bool = False,
):
    from CoREMOF.targets import (
        TargetMergedDataset,
        merge_targets,
        merge_targets_from_config,
    )

    if isinstance(target_data, TargetMergedDataset):
        merged = target_data
        try:
            merged.receipt()
        except Exception as error:
            raise TargetAttachmentError(
                "premerged target generation failed its sealed integrity check"
            ) from error
        target_base = getattr(merged, "base_dataset", merged)
        if _canonical_json(_release_binding(target_base)) != _canonical_json(
            _release_binding(dataset)
        ):
            raise TargetAttachmentError(
                "premerged target data does not match the assignment release binding"
            )
        return merged
    if isinstance(target_data, (str, Path)) and _config_like(Path(target_data)):
        if alias_registry is not None or feature_tables:
            raise TargetAttachmentError(
                "a target config already declares aliases/features; do not also pass overrides"
            )
        return merge_targets_from_config(
            dataset,
            target_data,
            verify_cif_files=verify_cif_files,
        )
    return merge_targets(
        dataset,
        target_data,
        alias_registry=alias_registry,
        feature_tables=feature_tables,
        verify_cif_files=verify_cif_files,
    )


def _partition_mapping(value: object) -> Tuple[Mapping[str, str], str]:
    assignments = getattr(value, "assignments", None)
    if isinstance(assignments, Mapping):
        result = {str(key): str(item) for key, item in assignments.items()}
        digest = getattr(value, "assignment_digest", None)
        if not isinstance(digest, str):
            method = getattr(value, "_assignment_sha256", None)
            digest = method() if callable(method) else None
        if not isinstance(digest, str):
            receipt = getattr(value, "receipt", None)
            receipt_value = receipt() if callable(receipt) else {}
            digest = (
                receipt_value.get("assignment_sha256")
                if isinstance(receipt_value, Mapping)
                else None
            )
        if not isinstance(digest, str):
            digest = _digest(sorted(result.items()))
        return MappingProxyType(result), _sha256(digest, "assignment digest")
    if isinstance(value, Mapping):
        if all(isinstance(item, (list, tuple)) for item in value.values()):
            result = {}
            for partition, ids in value.items():
                if (
                    not isinstance(partition, str)
                    or not partition
                    or partition != partition.strip()
                ):
                    raise TargetAttachmentError(
                        "partition names must be exact nonblank strings"
                    )
                for structure_id in _exact_ids(ids, "partition IDs"):
                    if structure_id in result:
                        raise TargetAttachmentError(
                            "structure ID appears in multiple partitions: {}".format(
                                structure_id
                            )
                        )
                    result[structure_id] = partition
        elif all(isinstance(item, str) for item in value.values()):
            result = {}
            for structure_id, partition in value.items():
                if (
                    not isinstance(structure_id, str)
                    or not structure_id
                    or structure_id != structure_id.strip()
                    or not isinstance(partition, str)
                    or not partition
                    or partition != partition.strip()
                ):
                    raise TargetAttachmentError(
                        "ID-to-partition mappings require exact nonblank strings"
                    )
                result[structure_id] = partition
        else:
            raise TargetAttachmentError(
                "partition mapping must be partition->ID-list or ID->partition"
            )
        return MappingProxyType(result), _digest(sorted(result.items()))
    ids = _exact_ids(value)
    result = {structure_id: "selected" for structure_id in ids}
    return MappingProxyType(result), _digest(sorted(result.items()))


def frozen_assignment_manifest(
    rows: Sequence[Mapping[str, object]], receipt: Mapping[str, object]
) -> FrozenAssignmentManifest:
    """Verify one persisted assignment CSV against its paired receipt.

    New exploratory manifests, historical splitter manifests, and a single
    per-run benchmark manifest are accepted.  The returned wrapper retains the
    receipt's frozen digest even though excluded rows are not target-joined.
    """

    if not isinstance(rows, (list, tuple)) or not rows:
        raise TargetAttachmentError("assignment manifest must contain at least one row")
    if not isinstance(receipt, Mapping):
        raise TargetAttachmentError("assignment receipt must be a JSON object")
    partitions_seen = set()
    run_keys = set()
    assignments: Dict[str, str] = {}
    exclusions: Dict[str, str] = {}
    normalized_rows = []
    for row_number, row in enumerate(rows, start=2):
        if not isinstance(row, Mapping):
            raise TargetAttachmentError(
                "assignment manifest row {} is not a mapping".format(row_number)
            )
        structure_id = row.get("structure_id")
        partition_column = (
            "partition" if "partition" in row else "split" if "split" in row else None
        )
        if partition_column is None:
            raise TargetAttachmentError(
                "assignment manifest requires a partition or split column"
            )
        partition = row.get(partition_column)
        if (
            not isinstance(structure_id, str)
            or not structure_id
            or structure_id != structure_id.strip()
            or not isinstance(partition, str)
            or not partition
            or partition != partition.strip()
        ):
            raise TargetAttachmentError(
                "assignment manifest row {} has an invalid ID or partition".format(
                    row_number
                )
            )
        if structure_id in assignments or structure_id in exclusions:
            raise TargetAttachmentError(
                "assignment manifest repeats structure_id {}".format(structure_id)
            )
        reason = row.get("exclusion_reason", "")
        if reason is None:
            reason = ""
        if not isinstance(reason, str):
            raise TargetAttachmentError(
                "assignment manifest exclusion reasons must be strings"
            )
        if partition == "EXCLUDED":
            exclusions[structure_id] = reason
        else:
            assignments[structure_id] = partition
            partitions_seen.add(partition)
        run_key = row.get("run_key")
        if run_key not in (None, ""):
            if not isinstance(run_key, str) or run_key != run_key.strip():
                raise TargetAttachmentError(
                    "assignment manifest has an invalid run_key"
                )
            run_keys.add(run_key)
        normalized_rows.append((structure_id, partition, reason, row))

    schema = receipt.get("schema_version")
    if schema not in _ACCEPTED_ASSIGNMENT_RECEIPT_SCHEMAS:
        raise TargetAttachmentError(
            "assignment receipt schema is not accepted: {!r}".format(schema)
        )
    if schema == "coremof-data-split-receipt/1.0":
        if run_keys:
            raise TargetAttachmentError(
                "a data-split manifest must not declare benchmark run_key values"
            )
        expected = receipt.get("assignment_sha256")
        digest_kind = "new"
    elif schema == "coremof-cr-ncr-benchmark-suite-receipt/1.0":
        if "assignment_sha256" in receipt:
            raise TargetAttachmentError(
                "a benchmark-suite receipt must not declare a top-level "
                "assignment_sha256"
            )
        receipt_runs = receipt.get("runs")
        if not isinstance(receipt_runs, list):
            raise TargetAttachmentError(
                "a benchmark-suite receipt requires its per-run digest list"
            )
        if len(run_keys) != 1:
            raise TargetAttachmentError(
                "a benchmark suite manifest must contain exactly one run_key"
            )
        run_key = next(iter(run_keys))
        matches = [
            item
            for item in receipt_runs
            if isinstance(item, Mapping) and item.get("run_key") == run_key
        ]
        if len(matches) != 1:
            raise TargetAttachmentError(
                "assignment receipt has no unique entry for run_key {}".format(run_key)
            )
        expected = matches[0].get("assignment_sha256")
        digest_kind = "new"
    else:
        # Historical SplitResult receipts bind all diagnostic CSV fields.
        if run_keys:
            raise TargetAttachmentError(
                "a historical split manifest must not declare benchmark run_key values"
            )
        expected = receipt.get("assignment_sha256")
        digest_kind = "legacy"
    expected = _sha256(expected, "receipt assignment_sha256")

    if digest_kind == "legacy":
        lines = []
        for structure_id, partition, reason, row in sorted(normalized_rows):
            lines.append(
                "\t".join(
                    (
                        structure_id,
                        partition,
                        str(row.get("label", "") or ""),
                        str(row.get("parent_group", "") or ""),
                        str(row.get("parent_diagnostic", "") or ""),
                        str(row.get("leakage_group", "") or ""),
                        reason,
                    )
                )
            )
        observed = hashlib.sha256("\n".join(lines).encode("utf-8")).hexdigest()
    else:
        observed = _assignment_digest(assignments, exclusions)
    if observed != expected:
        raise TargetAttachmentError(
            "assignment manifest does not match receipt assignment_sha256"
        )

    receipt_partitions = receipt.get("partitions")
    if isinstance(receipt_partitions, Mapping):
        declared = {}
        for partition, value in receipt_partitions.items():
            ids = value.get("ids") if isinstance(value, Mapping) else value
            if isinstance(ids, (list, tuple)):
                for structure_id in ids:
                    declared[str(structure_id)] = str(partition)
        if declared and declared != assignments:
            raise TargetAttachmentError(
                "assignment manifest partitions differ from the paired receipt"
            )
    receipt_exclusions = receipt.get("exclusions")
    if isinstance(receipt_exclusions, Mapping):
        declared_exclusions = {
            str(key): str(value) for key, value in receipt_exclusions.items()
        }
    elif isinstance(receipt_exclusions, list):
        declared_exclusions = {
            str(item.get("structure_id")): str(item.get("reason", ""))
            for item in receipt_exclusions
            if isinstance(item, Mapping)
        }
    else:
        declared_exclusions = {}
    if declared_exclusions and declared_exclusions != exclusions:
        raise TargetAttachmentError(
            "assignment manifest exclusions differ from the paired receipt"
        )
    release_binding = receipt.get("release_binding")
    if not isinstance(release_binding, Mapping):
        cohort_receipt = receipt.get("cohort_receipt")
        if isinstance(cohort_receipt, Mapping):
            release_binding = cohort_receipt.get("release_binding")
    if (
        not isinstance(release_binding, Mapping)
        and schema != "coremof-split-receipt/1.1"
    ):
        raise TargetAttachmentError(
            "assignment receipt is missing its exact release binding"
        )
    if not isinstance(release_binding, Mapping):
        receipt_hashes = receipt.get("input_hashes")
        if not isinstance(receipt_hashes, Mapping) or not receipt_hashes:
            raise TargetAttachmentError(
                "assignment receipt has no release input-hash binding"
            )
        all_ids = tuple(
            sorted(structure_id for structure_id, _, _, _ in normalized_rows)
        )
        release_binding = {
            "dataset_version": receipt.get("dataset_version"),
            "structure_count": len(all_ids),
            "structure_ids_sha256": hashlib.sha256(
                "\n".join(all_ids).encode("utf-8")
            ).hexdigest(),
            "input_hashes": dict(receipt_hashes),
        }
    official_split = receipt.get("official_split")
    if official_split is None:
        cohort_receipt = receipt.get("cohort_receipt")
        if isinstance(cohort_receipt, Mapping):
            official_split = cohort_receipt.get("official_split")
    if official_split is not None and type(official_split) is not bool:
        raise TargetAttachmentError(
            "assignment receipt official_split must be boolean when present"
        )
    if official_split is True:
        raise TargetAttachmentError(
            "official assignment attachment is unavailable without a separately "
            "audited official manifest loader"
        )
    return _register_frozen_assignment(
        FrozenAssignmentManifest(
            assignments,
            expected,
            release_binding,
            official_split=official_split,
            exclusions=exclusions,
            _source_schema=str(schema),
            _factory_token=_FROZEN_ASSIGNMENT_FACTORY_TOKEN,
        )
    )


def _coverage(
    ids: Sequence[str],
    assignments: Mapping[str, str],
    target_columns: Sequence[str],
    values_by_id: Mapping[str, Mapping[str, object]],
) -> Tuple[Mapping[str, object], Mapping[str, object]]:
    partitions = tuple(sorted(set(assignments[value] for value in ids)))
    coverage = {}
    missing = {}
    for partition in partitions:
        partition_ids = tuple(value for value in ids if assignments[value] == partition)
        by_target = {
            target: sum(
                values_by_id[value].get(target) is not None for value in partition_ids
            )
            for target in target_columns
        }
        complete = sum(
            all(
                values_by_id[value].get(target) is not None for target in target_columns
            )
            for value in partition_ids
        )
        coverage[partition] = {
            "selected_structure_count": len(partition_ids),
            "complete_target_structure_count": complete,
            "available_by_target": by_target,
        }
        missing[partition] = {
            target: len(partition_ids) - by_target[target] for target in target_columns
        }
    return _freeze(coverage), _freeze(missing)


@dataclass(frozen=True)
class TargetAttachedView:
    """Immutable target view over one unchanged frozen assignment."""

    structure_ids: Tuple[str, ...]
    original_structure_ids: Tuple[str, ...]
    assignments: Mapping[str, str]
    target_columns: Tuple[str, ...]
    values_by_id: Mapping[str, Mapping[str, object]]
    provenance_by_id: Mapping[str, Mapping[str, object]]
    target_definitions: Mapping[str, Mapping[str, object]]
    coverage_by_partition: Mapping[str, object]
    missing_counts_by_partition: Mapping[str, object]
    derived_coverage_by_partition: Mapping[str, object]
    derived_missing_counts_by_partition: Mapping[str, object]
    missing_policy: str
    dropped_ids: Tuple[str, ...]
    original_assignment_digest: str
    derived_view_digest: str
    target_input_receipt: Mapping[str, object]
    official_split: Optional[bool]
    _receipt: Mapping[str, object] = field(repr=False, compare=False)
    _authority_factory_token: object = field(default=None, repr=False, compare=False)

    def __post_init__(self) -> None:
        for name in (
            "structure_ids",
            "original_structure_ids",
            "target_columns",
            "dropped_ids",
        ):
            object.__setattr__(self, name, tuple(getattr(self, name)))
        for name in (
            "assignments",
            "values_by_id",
            "provenance_by_id",
            "target_definitions",
            "coverage_by_partition",
            "missing_counts_by_partition",
            "derived_coverage_by_partition",
            "derived_missing_counts_by_partition",
            "target_input_receipt",
            "_receipt",
        ):
            object.__setattr__(self, name, _freeze(getattr(self, name)))

    def receipt(self) -> Mapping[str, object]:
        _require_attached_result(self)
        return _jsonable(self._receipt)  # type: ignore[return-value]

    def rows(self) -> Tuple[Mapping[str, object], ...]:
        _require_attached_result(self)
        return tuple(
            {
                "structure_id": structure_id,
                "partition": self.assignments[structure_id],
                **{
                    target: self.values_by_id[structure_id].get(target)
                    for target in self.target_columns
                },
            }
            for structure_id in self.structure_ids
        )

    @staticmethod
    def _csv_cell(value: object) -> object:
        if value is None:
            return ""
        if isinstance(value, bool):
            return "true" if value else "false"
        if isinstance(value, (Mapping, list, tuple)):
            return _canonical_json(value)
        return value

    def write(
        self,
        output_directory: object,
        stem: str = "coremof_attached_targets",
        overwrite: bool = False,
    ) -> Tuple[Path, Path, Path]:
        """Transactionally write target rows, provenance, and receipt."""

        _require_attached_result(self)
        if (
            not isinstance(stem, str)
            or not stem
            or stem in {".", ".."}
            or Path(stem).name != stem
        ):
            raise TargetAttachmentError(
                "stem must be a filename stem without directories"
            )
        if type(overwrite) is not bool:
            raise TypeError("overwrite must be a boolean")
        directory = Path(output_directory)
        directory.mkdir(parents=True, exist_ok=True)
        targets = (
            directory / (stem + ".csv"),
            directory / (stem + ".provenance.jsonl"),
            directory / (stem + ".json"),
        )
        if not overwrite and any(path.exists() for path in targets):
            raise FileExistsError("Refusing to overwrite target-attached output")
        staging = Path(tempfile.mkdtemp(prefix=".{}-".format(stem), dir=str(directory)))
        staged = tuple(staging / path.name for path in targets)
        preserve_staging = False
        try:
            fields = ("structure_id", "partition") + self.target_columns
            with staged[0].open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                for row in self.rows():
                    writer.writerow(
                        {key: self._csv_cell(row.get(key)) for key in fields}
                    )
                handle.flush()
                os.fsync(handle.fileno())
            with staged[1].open("w", encoding="utf-8") as handle:
                for structure_id in self.structure_ids:
                    by_target = self.provenance_by_id.get(structure_id, {})
                    for target in sorted(by_target):
                        observations = by_target[target]
                        if not isinstance(observations, (list, tuple)):
                            observations = (observations,)
                        for observation in observations:
                            handle.write(
                                _canonical_json(
                                    {
                                        "structure_id": structure_id,
                                        "partition": self.assignments[structure_id],
                                        "target": target,
                                        "observation": observation,
                                    }
                                )
                                + "\n"
                            )
                handle.flush()
                os.fsync(handle.fileno())
            with staged[2].open("w", encoding="utf-8") as handle:
                json.dump(
                    self.receipt(), handle, indent=2, sort_keys=True, allow_nan=False
                )
                handle.write("\n")
                handle.flush()
                os.fsync(handle.fileno())
            from CoREMOF._transactions import publish_file_bundle

            publish_file_bundle(staged, targets, overwrite=overwrite)
        except BaseException as error:
            preserve_staging = hasattr(error, "coremof_preserved_staging_directory")
            raise
        finally:
            if not preserve_staging:
                shutil.rmtree(staging, ignore_errors=True)
        return targets


def _attach_one(
    assignment: object,
    merged: object,
    missing: str,
    release_binding: Mapping[str, object],
    official_split: Optional[bool],
) -> TargetAttachedView:
    if missing not in {"keep", "error", "drop"}:
        raise TargetAttachmentError("missing must be 'keep', 'error', or 'drop'")
    assignments, assignment_digest = _partition_mapping(assignment)
    original_ids = tuple(sorted(assignments))
    release_ids = set(getattr(merged, "structure_ids", ()))
    unknown = set(original_ids).difference(release_ids)
    if unknown:
        raise TargetAttachmentError(
            "assignment IDs are absent from the target release: {}".format(
                ", ".join(sorted(unknown)[:5])
            )
        )
    target_columns = tuple(getattr(merged, "target_columns"))
    collisions = sorted(_RESERVED_EXPORT_COLUMNS.intersection(target_columns))
    if collisions:
        raise TargetAttachmentError(
            "target columns collide with reserved attachment export columns: {}".format(
                ", ".join(collisions)
            )
        )
    values_source = getattr(merged, "target_values_by_id")
    values = {
        structure_id: {
            target: values_source[structure_id].get(target) for target in target_columns
        }
        for structure_id in original_ids
    }
    incomplete = tuple(
        structure_id
        for structure_id in original_ids
        if not all(
            values[structure_id].get(target) is not None for target in target_columns
        )
    )
    if missing == "error" and incomplete:
        by_target = {
            target: sum(values[value].get(target) is None for value in original_ids)
            for target in target_columns
        }
        raise TargetAttachmentError(
            "missing='error' requires a completed target set: {} of {} selected "
            "structures are incomplete; missing counts by target: {}".format(
                len(incomplete),
                len(original_ids),
                ", ".join(
                    "{}={}".format(target, count) for target, count in by_target.items()
                ),
            )
        )
    coverage, missing_counts = _coverage(
        original_ids,
        assignments,
        target_columns,
        values,
    )
    dropped = incomplete if missing == "drop" else ()
    dropped_set = set(dropped)
    selected = tuple(value for value in original_ids if value not in dropped_set)
    selected_assignments = {
        structure_id: assignments[structure_id] for structure_id in selected
    }
    selected_values = {structure_id: values[structure_id] for structure_id in selected}
    derived_coverage, derived_missing_counts = _coverage(
        selected,
        selected_assignments,
        target_columns,
        selected_values,
    )
    provenance_source = getattr(merged, "target_provenance_by_id", {})
    provenance = {
        structure_id: provenance_source.get(structure_id, {})
        for structure_id in selected
    }
    derived_digest = _digest(
        {
            "original_assignment_sha256": assignment_digest,
            "missing_policy": missing,
            "selected_assignments": sorted(selected_assignments.items()),
            "target_values": selected_values,
        }
    )
    target_receipt = getattr(merged, "target_input_receipt")
    receipt = {
        "schema_version": "coremof-target-attachment-receipt/1.0",
        "target_attachment_api_version": TARGET_ATTACHMENT_API_VERSION,
        "implementation": {
            "package": "CoREMOF-tools",
            "package_version": __version__,
            "source_sha256": dict(_source_hashes()),
        },
        "join_policy_definition": (
            "left join by exact current public structure_id or an explicitly "
            "hash-bound alias; no fuzzy matching, unit/condition inference, "
            "imputation, refill, rebalance, or resplit"
        ),
        "missing_policy": missing,
        "missing_policy_definition": _MISSING_POLICY_DEFINITIONS[missing],
        "missing_policy_definitions": dict(_MISSING_POLICY_DEFINITIONS),
        "original_assignment_sha256": assignment_digest,
        "original_selected_structure_count": len(original_ids),
        "derived_selected_structure_count": len(selected),
        "dropped_structure_count": len(dropped),
        "dropped_ids": list(dropped),
        "target_columns": list(target_columns),
        "target_definitions": _jsonable(getattr(merged, "target_definitions")),
        "coverage_by_partition": _jsonable(coverage),
        "missing_counts_by_partition": _jsonable(missing_counts),
        "derived_coverage_by_partition": _jsonable(derived_coverage),
        "derived_missing_counts_by_partition": _jsonable(derived_missing_counts),
        "release_binding": _jsonable(release_binding),
        "release_binding_sha256": _digest(release_binding),
        "official_split": official_split,
        "official_split_definition": (
            "source frozen assignment status; null only for a raw ID list or "
            "historical receipt that did not record this field"
        ),
        "target_input_receipt": _jsonable(target_receipt),
        "derived_view_sha256": derived_digest,
        "original_split_receipt_changed": False,
    }
    return _register_attached_result(
        TargetAttachedView(
            structure_ids=selected,
            original_structure_ids=original_ids,
            assignments=selected_assignments,
            target_columns=target_columns,
            values_by_id=selected_values,
            provenance_by_id=provenance,
            target_definitions=getattr(merged, "target_definitions"),
            coverage_by_partition=coverage,
            missing_counts_by_partition=missing_counts,
            derived_coverage_by_partition=derived_coverage,
            derived_missing_counts_by_partition=derived_missing_counts,
            missing_policy=missing,
            dropped_ids=dropped,
            original_assignment_digest=assignment_digest,
            derived_view_digest=derived_digest,
            target_input_receipt=target_receipt,
            official_split=official_split,
            _receipt=receipt,
            _authority_factory_token=_ATTACHED_RESULT_FACTORY_TOKEN,
        )
    )


@dataclass(frozen=True)
class TargetAttachedBenchmarkSuite:
    """Targets attached independently to every frozen benchmark run."""

    run_views: Mapping[str, TargetAttachedView]
    target_columns: Tuple[str, ...]
    missing_policy: str
    original_assignment_digest: str
    target_input_receipt: Mapping[str, object]
    official_split: Optional[bool]
    _receipt: Mapping[str, object] = field(repr=False, compare=False)
    _authority_factory_token: object = field(default=None, repr=False, compare=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "run_views", _freeze(self.run_views))
        object.__setattr__(self, "target_columns", tuple(self.target_columns))
        object.__setattr__(
            self, "target_input_receipt", _freeze(self.target_input_receipt)
        )
        object.__setattr__(self, "_receipt", _freeze(self._receipt))

    def receipt(self) -> Mapping[str, object]:
        _require_attached_result(self)
        return _jsonable(self._receipt)  # type: ignore[return-value]

    def write(
        self,
        output_directory: object,
        stem: str = "coremof_cr_ncr_targets",
        overwrite: bool = False,
    ) -> Path:
        _require_attached_result(self)
        if (
            not isinstance(stem, str)
            or not stem
            or stem in {".", ".."}
            or Path(stem).name != stem
        ):
            raise TargetAttachmentError(
                "stem must be a filename stem without directories"
            )
        if type(overwrite) is not bool:
            raise TypeError("overwrite must be a boolean")
        parent = Path(output_directory)
        parent.mkdir(parents=True, exist_ok=True)
        target = parent / stem
        if target.exists() and not overwrite:
            raise FileExistsError("Refusing to overwrite attached benchmark output")
        staging = Path(tempfile.mkdtemp(prefix=".{}-".format(stem), dir=str(parent)))
        try:
            fields = ("run_key", "structure_id", "partition") + self.target_columns
            with (staging / "attached_membership.csv").open(
                "w", encoding="utf-8", newline=""
            ) as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                for run_key in sorted(self.run_views):
                    view = self.run_views[run_key]
                    for row in view.rows():
                        writer.writerow(
                            {
                                "run_key": run_key,
                                **{
                                    key: TargetAttachedView._csv_cell(row.get(key))
                                    for key in fields
                                    if key != "run_key"
                                },
                            }
                        )
            with (staging / "provenance.jsonl").open("w", encoding="utf-8") as handle:
                for run_key in sorted(self.run_views):
                    view = self.run_views[run_key]
                    for structure_id in view.structure_ids:
                        by_target = view.provenance_by_id.get(structure_id, {})
                        for target_name in sorted(by_target):
                            observations = by_target[target_name]
                            if not isinstance(observations, (list, tuple)):
                                observations = (observations,)
                            for observation in observations:
                                handle.write(
                                    _canonical_json(
                                        {
                                            "run_key": run_key,
                                            "structure_id": structure_id,
                                            "partition": view.assignments[structure_id],
                                            "target": target_name,
                                            "observation": observation,
                                        }
                                    )
                                    + "\n"
                                )
            with (staging / "coverage.json").open("w", encoding="utf-8") as handle:
                json.dump(
                    {
                        run_key: {
                            "coverage_by_partition": _jsonable(
                                view.coverage_by_partition
                            ),
                            "missing_counts_by_partition": _jsonable(
                                view.missing_counts_by_partition
                            ),
                            "derived_coverage_by_partition": _jsonable(
                                view.derived_coverage_by_partition
                            ),
                            "derived_missing_counts_by_partition": _jsonable(
                                view.derived_missing_counts_by_partition
                            ),
                            "dropped_structure_count": len(view.dropped_ids),
                            "dropped_ids": list(view.dropped_ids),
                            "original_assignment_sha256": view.original_assignment_digest,
                            "derived_view_sha256": view.derived_view_digest,
                        }
                        for run_key, view in sorted(self.run_views.items())
                    },
                    handle,
                    indent=2,
                    sort_keys=True,
                    allow_nan=False,
                )
                handle.write("\n")
            with (staging / "receipt.json").open("w", encoding="utf-8") as handle:
                json.dump(
                    self.receipt(), handle, indent=2, sort_keys=True, allow_nan=False
                )
                handle.write("\n")
            files = sorted(path for path in staging.iterdir() if path.is_file())
            with (staging / "SHA256SUMS").open("w", encoding="utf-8") as handle:
                for path in files:
                    handle.write(
                        "{}  {}\n".format(
                            hashlib.sha256(path.read_bytes()).hexdigest(), path.name
                        )
                    )
            from CoREMOF._transactions import publish_directory

            return publish_directory(staging, target, overwrite=overwrite)
        except BaseException:
            shutil.rmtree(staging, ignore_errors=True)
            raise


def attach_targets(
    split_or_id_list: object,
    target_data: object,
    dataset: object = None,
    missing: str = "keep",
    alias_registry: object = None,
    feature_tables: Sequence[str] = (),
    verify_cif_files: bool = False,
):
    """Attach targets to a split, suite, partition mapping, or raw ID list."""

    _source_hashes()
    if type(verify_cif_files) is not bool:
        raise TypeError("verify_cif_files must be a boolean")
    source = _dataset_from_assignment(split_or_id_list, explicit=dataset)
    release_binding = _release_binding(source)
    official_split = _validate_live_assignment_authority(
        split_or_id_list, release_binding
    )
    merged = _merged_targets(
        source,
        target_data,
        alias_registry=alias_registry,
        feature_tables=feature_tables,
        verify_cif_files=verify_cif_files,
    )
    runs = getattr(split_or_id_list, "runs", None)
    suite_digest = getattr(split_or_id_list, "assignment_digest", None)
    if isinstance(runs, (list, tuple)) and isinstance(suite_digest, str):
        views = {
            run.run_key: _attach_one(
                run, merged, missing, release_binding, official_split
            )
            for run in runs
        }
        receipt = {
            "schema_version": "coremof-target-attached-benchmark-suite-receipt/1.0",
            "target_attachment_api_version": TARGET_ATTACHMENT_API_VERSION,
            "implementation": {
                "package": "CoREMOF-tools",
                "package_version": __version__,
                "source_sha256": dict(_source_hashes()),
            },
            "missing_policy": missing,
            "missing_policy_definition": _MISSING_POLICY_DEFINITIONS[missing],
            "missing_policy_definitions": dict(_MISSING_POLICY_DEFINITIONS),
            "original_suite_assignment_sha256": suite_digest,
            "run_count": len(views),
            "target_columns": list(getattr(merged, "target_columns")),
            "target_input_receipt": _jsonable(getattr(merged, "target_input_receipt")),
            "coverage_by_run_and_partition": {
                run_key: _jsonable(view.coverage_by_partition)
                for run_key, view in sorted(views.items())
            },
            "missing_counts_by_run_and_partition": {
                run_key: _jsonable(view.missing_counts_by_partition)
                for run_key, view in sorted(views.items())
            },
            "derived_coverage_by_run_and_partition": {
                run_key: _jsonable(view.derived_coverage_by_partition)
                for run_key, view in sorted(views.items())
            },
            "derived_missing_counts_by_run_and_partition": {
                run_key: _jsonable(view.derived_missing_counts_by_partition)
                for run_key, view in sorted(views.items())
            },
            "dropped_by_run": {
                run_key: {
                    "count": len(view.dropped_ids),
                    "ids": list(view.dropped_ids),
                }
                for run_key, view in sorted(views.items())
            },
            "release_binding": _jsonable(release_binding),
            "release_binding_sha256": _digest(release_binding),
            "official_split": official_split,
            "official_split_definition": (
                "source frozen benchmark-suite status; null only when the "
                "source did not record this field"
            ),
            "original_split_receipts_changed": False,
        }
        return _register_attached_result(
            TargetAttachedBenchmarkSuite(
                run_views=views,
                target_columns=tuple(getattr(merged, "target_columns")),
                missing_policy=missing,
                original_assignment_digest=suite_digest,
                target_input_receipt=getattr(merged, "target_input_receipt"),
                official_split=official_split,
                _receipt=receipt,
                _authority_factory_token=_ATTACHED_RESULT_FACTORY_TOKEN,
            )
        )
    return _attach_one(
        split_or_id_list,
        merged,
        missing,
        release_binding,
        official_split,
    )


__all__ = [
    "TARGET_ATTACHMENT_API_VERSION",
    "FrozenAssignmentManifest",
    "TargetAttachedBenchmarkSuite",
    "TargetAttachedView",
    "TargetAttachmentError",
    "attach_targets",
    "frozen_assignment_manifest",
]
