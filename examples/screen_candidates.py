#!/usr/bin/env python3
"""Rank finite CoRE-MOF screening candidates and optionally split them safely.

This example intentionally uses only the Python standard library plus the
public :mod:`CoREMOF` APIs.  Checker consensus, release validation, target
joining, parent resolution, and leakage-safe splitting remain package-owned.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from dataclasses import dataclass, field
from decimal import Decimal, InvalidOperation
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile
from types import MappingProxyType
from typing import Dict, Iterable, Mapping, Optional, Sequence, Tuple, Union

# Make the checked-out example directly runnable before an editable install.
# An installed copy is otherwise imported normally.
_REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if (_REPOSITORY_ROOT / "CoREMOF").is_dir():
    repository_text = str(_REPOSITORY_ROOT)
    if repository_text not in sys.path:
        sys.path.insert(0, repository_text)

from CoREMOF import __version__
from CoREMOF._authority import (
    AuthorityStateError,
    IdentitySealRegistry,
    reject_sealed_copy,
    state_fingerprint,
)
import CoREMOF.dataset as _coremof_dataset_module
from CoREMOF.dataset import (
    CoREMOFDataset,
    ReleaseValidationError,
    _is_authenticated_official_checker_view,
    _release_receipt_state,
    _validate_classified_generation_if_present,
    _validate_dataset_generation_if_present,
)
from CoREMOF.labels import CHECKER_PRESETS
from CoREMOF.parents import LEAKAGE_GUARD_CHOICES, SELECTABLE_PARENT_METHODS
from CoREMOF.targets import TargetDataError, merge_targets_from_config


SCHEMA_VERSION = "coremof-screening-receipt/1.0"
_COREMOF_PACKAGE_ROOT = Path(_coremof_dataset_module.__file__).resolve().parent
_CLASSIFICATION_SOURCE_FILES = ("_authority.py", "dataset.py", "labels.py")
_TARGET_SOURCE_FILES = _CLASSIFICATION_SOURCE_FILES + ("targets.py",)
_SPLIT_SOURCE_FILES = _CLASSIFICATION_SOURCE_FILES + (
    "parents.py",
    "splitters.py",
)
_REQUIRED_TARGET_CSV_PREFIX = "required_target:"
_MISSING_TEXT = frozenset({""})
_SCREENING_RESULT_FACTORY_TOKEN = object()
_SCREENING_RESULTS = IdentitySealRegistry("screening result generation")
_NONFINITE_TEXT = frozenset(
    {
        "nan",
        "inf",
        "+inf",
        "-inf",
        "infinity",
        "+infinity",
        "-infinity",
    }
)


class ScreeningError(ValueError):
    """Raised when a requested screening operation cannot be audited safely."""


def _exact_nonblank_string(value: object, name: str) -> str:
    """Return one already-canonical string without coercion or repair."""

    if (
        not isinstance(value, str)
        or not value
        or value != value.strip()
    ):
        raise ScreeningError("{} must be an exact nonblank string".format(name))
    return value


def _optional_exact_string(value: object, name: str) -> Optional[str]:
    """Validate an optional exact nonblank string receipt claim."""

    if value is None:
        return None
    return _exact_nonblank_string(value, name)


def _screen_release_receipt_state(
    dataset: object,
) -> Tuple[Optional[str], Optional[str], bool]:
    """Translate the package's closed release-state contract to screen errors."""

    try:
        return _release_receipt_state(dataset)
    except (TypeError, ValueError) as exc:
        raise ScreeningError(str(exc)) from exc


def _deep_freeze(value: object) -> object:
    """Copy one result value into a recursively immutable representation."""

    if isinstance(value, Mapping):
        result = {}
        for key, item in value.items():
            if not isinstance(key, str) or not key or key != key.strip():
                raise ScreeningError("mapping keys must be exact nonblank strings")
            result[key] = _deep_freeze(item)
        return MappingProxyType(result)
    if isinstance(value, (list, tuple)):
        return tuple(_deep_freeze(item) for item in value)
    if isinstance(value, (set, frozenset)):
        raise ScreeningError("unordered sets are not valid screening result data")
    return value


@dataclass(frozen=True)
class ScreeningResult:
    """In-memory ranked result before files are published."""

    rows: Tuple[Mapping[str, object], ...]
    selected_ids: Tuple[str, ...]
    receipt: Mapping[str, object]
    classified_dataset: object
    filters: Mapping[str, object]
    _authority_factory_token: object = field(
        default=None, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        """Detach result evidence from mutable caller-owned containers."""

        selected_ids = tuple(self.selected_ids)
        for structure_id in selected_ids:
            _exact_nonblank_string(structure_id, "selected structure_id")
        object.__setattr__(
            self,
            "rows",
            tuple(_deep_freeze(row) for row in self.rows),
        )
        object.__setattr__(
            self,
            "selected_ids",
            selected_ids,
        )
        object.__setattr__(self, "receipt", _deep_freeze(self.receipt))
        object.__setattr__(self, "filters", _deep_freeze(self.filters))

    def __copy__(self):
        reject_sealed_copy(self, "copied")
        raise TypeError("unsealed screening results cannot be copied")

    def __deepcopy__(self, memo):
        reject_sealed_copy(self, "deep-copied")
        raise TypeError("unsealed screening results cannot be deep-copied")

    def __reduce_ex__(self, protocol):
        reject_sealed_copy(self, "pickled")
        raise TypeError("screening results are point-in-time non-serializable objects")


RankingNumber = Union[int, float, Decimal]


def _as_tuple(values: Optional[Iterable[str]]) -> Optional[Tuple[str, ...]]:
    if values is None:
        return None
    if isinstance(values, str):
        result = (values,)
    else:
        if not isinstance(values, (list, tuple)):
            raise ScreeningError(
                "filter values must be an exact string or ordered list/tuple "
                "of strings"
            )
        result = tuple(values)
    if not result or any(not isinstance(value, str) for value in result):
        raise ScreeningError("filter values must be non-empty strings")
    if any(not value or value != value.strip() for value in result):
        raise ScreeningError("filter values must be exact non-empty strings")
    if len(set(result)) != len(result):
        raise ScreeningError("filter values must not contain duplicates")
    return result


def _optional_mapping_attribute(value: object, name: str) -> Mapping[str, object]:
    """Return an optional mapping attribute without truthiness coercion."""

    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise ScreeningError("{} must be a mapping when present".format(name))
    return value


def _input_hash_mapping(dataset) -> Mapping[str, str]:
    values = _optional_mapping_attribute(
        getattr(dataset, "input_hashes", None), "input_hashes"
    )
    if any(
        not isinstance(key, str)
        or not key
        or key != key.strip()
        or not isinstance(value, str)
        or not value
        or value != value.strip()
        for key, value in values.items()
    ):
        raise ScreeningError(
            "input_hashes keys and values must be exact non-empty strings"
        )
    return values  # type: ignore[return-value]


def _column_name_tuple(value: object, name: str) -> Tuple[str, ...]:
    """Validate an optional dataset-like column-name sequence."""

    if value is None:
        return ()
    if isinstance(value, str) or not isinstance(value, (list, tuple)):
        raise ScreeningError(
            "{} must be a list or tuple of strings when present".format(name)
        )
    result = tuple(value)
    if any(
        not isinstance(item, str) or not item or item != item.strip()
        for item in result
    ):
        raise ScreeningError(
            "{} values must be exact nonblank strings".format(name)
        )
    if len(set(result)) != len(result):
        raise ScreeningError("{} must not contain duplicate names".format(name))
    return result


def _boolean_attribute(obj: object, name: str, default: bool) -> bool:
    marker = object()
    value = getattr(obj, name, marker)
    if value is marker:
        return default
    if type(value) is not bool:
        raise ScreeningError("{} must be a boolean when present".format(name))
    return value


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def _current_coremof_source_hashes(filenames: Sequence[str]) -> Dict[str, str]:
    """Hash a closed, named CoREMOF source-file set."""

    return {
        filename: _sha256_file(_COREMOF_PACKAGE_ROOT / filename)
        for filename in filenames
    }


_EXAMPLE_SCRIPT_PATH = Path(__file__).resolve()
_IMPORTED_EXAMPLE_SCRIPT_SHA256 = _sha256_file(_EXAMPLE_SCRIPT_PATH)
_IMPORTED_CLASSIFICATION_SOURCE_SHA256 = _current_coremof_source_hashes(
    _CLASSIFICATION_SOURCE_FILES
)


def _verify_imported_screening_sources() -> Dict[str, str]:
    """Fail if executing screening sources drifted after module import."""

    current_script = _sha256_file(_EXAMPLE_SCRIPT_PATH)
    if current_script != _IMPORTED_EXAMPLE_SCRIPT_SHA256:
        raise ScreeningError(
            "screening example script changed after module import; restart the run"
        )
    current_classification = _current_coremof_source_hashes(
        _CLASSIFICATION_SOURCE_FILES
    )
    for filename in _CLASSIFICATION_SOURCE_FILES:
        if (
            current_classification[filename]
            != _IMPORTED_CLASSIFICATION_SOURCE_SHA256[filename]
        ):
            raise ScreeningError(
                "CoREMOF source {} changed after module import; restart the run".format(
                    filename
                )
            )
    return dict(_IMPORTED_CLASSIFICATION_SOURCE_SHA256)


def _receipt_source_hashes(
    receipt: Mapping[str, object],
    expected_filenames: Sequence[str],
    description: str,
) -> Dict[str, str]:
    """Read and validate one package-owned implementation hash contract."""

    implementation = receipt.get("implementation")
    if not isinstance(implementation, Mapping):
        raise ScreeningError(
            "{} receipt has no implementation mapping".format(description)
        )
    source_hashes = implementation.get("source_sha256")
    if not isinstance(source_hashes, Mapping):
        raise ScreeningError(
            "{} receipt has no source_sha256 mapping".format(description)
        )
    expected = set(expected_filenames)
    observed = {str(filename) for filename in source_hashes}
    if observed != expected:
        raise ScreeningError(
            "{} source closure differs from the expected files: {}".format(
                description, ", ".join(sorted(observed)) or "<empty>"
            )
        )
    result = {}
    for filename in expected_filenames:
        digest = source_hashes.get(filename)
        if (
            not isinstance(digest, str)
            or len(digest) != 64
            or any(character not in "0123456789abcdef" for character in digest)
        ):
            raise ScreeningError(
                "{} has an invalid SHA-256 for {}".format(description, filename)
            )
        result[filename] = digest
    return result


def _screening_source_hashes(
    target_receipt: Optional[Mapping[str, object]],
) -> Dict[str, str]:
    """Freeze the exact classification and optional target-merge sources."""

    classification = _verify_imported_screening_sources()
    if target_receipt is None:
        return classification
    target_hashes = _receipt_source_hashes(
        target_receipt, _TARGET_SOURCE_FILES, "target-merge"
    )
    for filename, digest in classification.items():
        if target_hashes[filename] != digest:
            raise ScreeningError(
                "CoREMOF source {} changed between target merge and screening".format(
                    filename
                )
            )
    return target_hashes


def _verify_result_implementation(
    receipt: Mapping[str, object],
) -> None:
    """Require one result receipt to match the still-loaded implementation."""

    classification = _verify_imported_screening_sources()
    implementation = receipt.get("implementation")
    if not isinstance(implementation, Mapping):
        raise ScreeningError("screening receipt has no implementation mapping")
    if implementation.get("example_script_sha256") != _IMPORTED_EXAMPLE_SCRIPT_SHA256:
        raise ScreeningError(
            "screening result was created by a different example script"
        )
    source_hashes = implementation.get("coremof_source_sha256")
    if not isinstance(source_hashes, Mapping):
        raise ScreeningError("screening receipt has no CoREMOF source hashes")
    with_targets = receipt.get("target_merge_receipt") is not None
    expected_files = _TARGET_SOURCE_FILES if with_targets else _CLASSIFICATION_SOURCE_FILES
    if {str(filename) for filename in source_hashes} != set(expected_files):
        raise ScreeningError("screening receipt has an invalid CoREMOF source closure")
    current_hashes = _current_coremof_source_hashes(expected_files)
    for filename in expected_files:
        if source_hashes.get(filename) != current_hashes[filename]:
            raise ScreeningError(
                "screening result was created by a different {}".format(filename)
            )
    for filename, digest in classification.items():
        if current_hashes[filename] != digest:
            raise ScreeningError(
                "screening source {} changed after module import".format(filename)
            )
    if with_targets:
        target_receipt = receipt.get("target_merge_receipt")
        if not isinstance(target_receipt, Mapping):
            raise ScreeningError("screening target receipt must be a mapping")
        target_hashes = _receipt_source_hashes(
            target_receipt, _TARGET_SOURCE_FILES, "target-merge"
        )
        if target_hashes != current_hashes:
            raise ScreeningError(
                "target-merge implementation differs from screening sources"
            )


def _split_source_hashes(
    screening_receipt: Mapping[str, object],
    split_receipt: Mapping[str, object],
) -> Dict[str, str]:
    """Validate and return the implementation closure used by one split."""

    with_targets = screening_receipt.get("target_merge_receipt") is not None
    expected = _SPLIT_SOURCE_FILES + (("targets.py",) if with_targets else ())
    split_hashes = _receipt_source_hashes(split_receipt, expected, "split")
    implementation = screening_receipt.get("implementation")
    if not isinstance(implementation, Mapping):
        raise ScreeningError("screening receipt has no implementation mapping")
    screening_hashes = implementation.get("coremof_source_sha256")
    if not isinstance(screening_hashes, Mapping):
        raise ScreeningError("screening receipt has no CoREMOF source hashes")
    for filename, digest in screening_hashes.items():
        if split_hashes.get(str(filename)) != digest:
            raise ScreeningError(
                "CoREMOF source {} changed between screening and splitting".format(
                    filename
                )
            )
    return split_hashes


def _required_target_csv_column(target: str) -> str:
    return _REQUIRED_TARGET_CSV_PREFIX + target


def _ids_sha256(values: Iterable[str]) -> str:
    return hashlib.sha256("\n".join(values).encode("utf-8")).hexdigest()


def _json_safe(value: object) -> object:
    if isinstance(value, Mapping):
        result = {}
        for key, item in value.items():
            if not isinstance(key, str) or not key or key != key.strip():
                raise ScreeningError("receipt mapping keys must be exact nonblank strings")
            result[key] = _json_safe(item)
        return result
    if isinstance(value, (set, frozenset)):
        raise ScreeningError("unordered sets are not valid receipt values")
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, Decimal):
        if not value.is_finite():
            raise ScreeningError("non-finite decimal cannot be written to a receipt")
        return str(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ScreeningError("non-finite value cannot be written to a receipt")
        return value
    raise ScreeningError(
        "receipt value type {} is unsupported without coercion".format(
            type(value).__name__
        )
    )


def _target_available(value: object) -> bool:
    """Return true only for an explicitly present, finite target value."""

    if value is None:
        return False
    if isinstance(value, Decimal):
        return value.is_finite()
    if isinstance(value, float):
        return math.isfinite(value)
    return True


def _finite_number(
    value: object,
) -> Tuple[Optional[RankingNumber], Optional[str]]:
    """Normalize one ranking value without losing integer/text precision.

    Exact integers remain integers, numeric text becomes :class:`Decimal`, and
    native floats remain floats.  Sorting converts those accepted forms to
    exact decimal comparison values separately.
    """

    if value is None:
        return None, "MISSING_RANK_VALUE"
    if isinstance(value, bool):
        return None, "NON_NUMERIC_RANK_VALUE"
    if isinstance(value, int):
        return value, None
    if isinstance(value, float):
        if not math.isfinite(value):
            return None, "NONFINITE_RANK_VALUE"
        if value == 0.0:
            value = 0.0
        return value, None
    if isinstance(value, Decimal):
        if not value.is_finite():
            return None, "NONFINITE_RANK_VALUE"
        if value.is_zero():
            value = Decimal(0)
        return value, None
    if isinstance(value, str):
        text = value.strip()
        normalized = text.casefold()
        if normalized in _MISSING_TEXT:
            return None, "MISSING_RANK_VALUE"
        if normalized in _NONFINITE_TEXT:
            return None, "NONFINITE_RANK_VALUE"
        try:
            number = Decimal(text)
        except (InvalidOperation, ValueError):
            return None, "NON_NUMERIC_RANK_VALUE"
        if not number.is_finite():
            return None, "NONFINITE_RANK_VALUE"
        if number.is_zero():
            number = Decimal(0)
        return number, None
    return None, "NON_NUMERIC_RANK_VALUE"


def _comparison_decimal(value: RankingNumber) -> Decimal:
    """Return the exact mathematical value used by the ranking comparator."""

    if isinstance(value, Decimal):
        return value
    if isinstance(value, int):
        return Decimal(value)
    return Decimal.from_float(value)


def _passes_required_targets(
    metadata: Mapping[str, object], required: Sequence[str], mode: str
) -> bool:
    available = [_target_available(metadata.get(target)) for target in required]
    return all(available) if mode == "all" else any(available)


def load_dataset(
    release: Path,
    target_config: Optional[Path] = None,
    verify_cif_files: bool = False,
):
    """Load a validated release, optionally attaching configured targets/features."""

    if target_config is None:
        return CoREMOFDataset.from_release(
            release, verify_cif_files=verify_cif_files
        )
    return merge_targets_from_config(
        release,
        target_config,
        verify_cif_files=verify_cif_files,
    )


def screen_dataset(
    dataset,
    *,
    rank_by: str,
    order: str = "descending",
    checkers: str = "5checker",
    labels: Optional[Iterable[str]] = ("CR",),
    sources: Optional[Iterable[str]] = None,
    variants: Optional[Iterable[str]] = None,
    metals: Optional[Iterable[str]] = None,
    required_targets: Optional[Iterable[str]] = None,
    required_target_mode: str = "all",
    limit: Optional[int] = None,
) -> ScreeningResult:
    """Filter a dataset first, then rank only finite numeric values.

    The returned object is write-free.  Scientific classification and the
    source/variant/metal filters are delegated to ``CoREMOFDataset``.
    """

    _validate_dataset_generation_if_present(dataset)
    if (
        not isinstance(rank_by, str)
        or not rank_by
        or rank_by != rank_by.strip()
    ):
        raise ScreeningError("rank_by must be a non-empty column name")
    if order not in ("ascending", "descending"):
        raise ScreeningError("order must be 'ascending' or 'descending'")
    if required_target_mode not in ("all", "any"):
        raise ScreeningError("required_target_mode must be 'all' or 'any'")
    if limit is not None and (
        isinstance(limit, bool) or not isinstance(limit, int) or limit <= 0
    ):
        raise ScreeningError("limit must be a positive integer")

    labels_tuple = _as_tuple(labels)
    sources_tuple = _as_tuple(sources)
    variants_tuple = _as_tuple(variants)
    metals_tuple = _as_tuple(metals)
    required_tuple = _as_tuple(required_targets) or ()

    (
        parent_input_status,
        dataset_input_status,
        _provisional_input,
    ) = _screen_release_receipt_state(dataset)
    input_hashes = _input_hash_mapping(dataset)
    cif_files_verified = _boolean_attribute(
        dataset, "cif_files_verified", False
    )
    dataset_version = _exact_nonblank_string(
        getattr(dataset, "dataset_version", None), "dataset_version"
    )

    target_columns = _column_name_tuple(
        getattr(dataset, "target_columns", None), "target_columns"
    )
    feature_columns = _column_name_tuple(
        getattr(dataset, "feature_columns", None), "feature_columns"
    )
    unknown_required = set(required_tuple).difference(target_columns)
    if unknown_required:
        if not target_columns:
            raise ScreeningError(
                "required targets need --target-config; no targets are attached"
            )
        raise ScreeningError(
            "unknown required target(s): {}".format(
                ", ".join(sorted(unknown_required))
            )
        )

    available_columns = (
        set(dataset.metadata_rows[0]) if dataset.metadata_rows else set()
    )
    if rank_by not in available_columns:
        raise ScreeningError("unknown ranking column {!r}".format(rank_by))

    if rank_by in target_columns:
        ranking_kind = "TARGET"
    elif rank_by in feature_columns:
        ranking_kind = "JOINED_FEATURE"
    else:
        ranking_kind = "RELEASE_METADATA"

    classified = dataset.classify(checkers)
    checker_view = _exact_nonblank_string(
        getattr(classified, "checker_view", None), "checker_view"
    )
    classified_dataset_version = _exact_nonblank_string(
        getattr(classified, "dataset_version", None),
        "classified dataset_version",
    )
    if classified_dataset_version != dataset_version:
        raise ScreeningError(
            "classified dataset_version differs from its source dataset"
        )
    checker_view_official_marker = object()
    checker_view_official_value = getattr(
        classified, "checker_view_official", checker_view_official_marker
    )
    if (
        checker_view_official_value is not checker_view_official_marker
        and type(checker_view_official_value) is not bool
    ):
        raise ScreeningError(
            "checker_view_official must be a boolean when present"
        )
    if checker_view_official_value is True and checker_view not in CHECKER_PRESETS:
        raise ScreeningError(
            "checker_view_official=True requires a canonical official checker preset"
        )
    if (
        checker_view_official_value is True
        and not _is_authenticated_official_checker_view(classified)
    ):
        raise ScreeningError(
            "checker_view_official=True requires an internally authenticated "
            "recomputed release view"
        )
    filtered = classified.filter(
        labels=labels_tuple,
        sources=sources_tuple,
        variants=variants_tuple,
        metals=metals_tuple,
    )

    required_ids = []
    required_excluded = []
    for record in filtered:
        structure_id = _exact_nonblank_string(
            getattr(record, "structure_id", None), "screening structure_id"
        )
        if not required_tuple or _passes_required_targets(
            record.metadata, required_tuple, required_target_mode
        ):
            required_ids.append(structure_id)
        else:
            required_excluded.append(structure_id)

    required_view = filtered.filter(structure_ids=required_ids)
    numeric = []
    ranking_exclusions: Dict[str, str] = {}
    for record in required_view:
        value, reason = _finite_number(record.metadata.get(rank_by))
        if reason is not None:
            ranking_exclusions[record.structure_id] = reason
            continue
        if value is None:  # pragma: no cover - guarded by the result contract
            raise ScreeningError("finite-number parser returned no value")
        numeric.append((record, value))

    if not numeric:
        raise ScreeningError(
            "no eligible structure has a finite numeric value for {!r}".format(
                rank_by
            )
        )
    # A stable two-pass sort keeps the structure-ID tie breaker ascending even
    # for descending values, without negating/rounding a Decimal sort key.
    numeric.sort(key=lambda item: item[0].structure_id)
    numeric.sort(
        key=lambda item: _comparison_decimal(item[1]),
        reverse=(order == "descending"),
    )
    finite_count = len(numeric)
    if limit is not None:
        numeric = numeric[:limit]

    rows = []
    for rank, (record, value) in enumerate(numeric, start=1):
        row = {
            "rank": rank,
            "structure_id": record.structure_id,
            "ranking_column": rank_by,
            "ranking_value": value,
            "checker_view": checker_view,
            "classification_label": record.label,
            "source_database": record.source_database,
            "source_id": record.metadata.get("source_id", ""),
            "structure_variant": record.structure_variant,
            "metal_elements": ";".join(record.metal_elements),
            "cif_file": record.metadata.get("cif_file", ""),
        }
        for target in sorted(required_tuple):
            row[_required_target_csv_column(target)] = record.metadata.get(target)
        rows.append(row)

    selected_ids = tuple(
        _exact_nonblank_string(row["structure_id"], "screening structure_id")
        for row in rows
    )
    exclusion_counts = Counter(ranking_exclusions.values())
    checker_view_official = _boolean_attribute(
        classified,
        "checker_view_official",
        False,
    )
    if checker_view_official and checker_view not in CHECKER_PRESETS:
        raise ScreeningError(
            "checker_view_official=True requires a canonical official checker preset"
        )
    if checker_view_official and not _is_authenticated_official_checker_view(
        classified
    ):
        raise ScreeningError(
            "checker_view_official=True requires an internally authenticated "
            "recomputed release view"
        )
    target_receipt = getattr(dataset, "target_input_receipt", None)
    coremof_source_hashes = _screening_source_hashes(target_receipt)
    filters = {
        "labels": list(labels_tuple) if labels_tuple is not None else None,
        "sources": list(sources_tuple) if sources_tuple is not None else None,
        "variants": list(variants_tuple) if variants_tuple is not None else None,
        "metals": list(metals_tuple) if metals_tuple is not None else None,
        "required_targets": list(required_tuple),
        "required_target_mode": required_target_mode,
    }
    receipt = {
        "schema_version": SCHEMA_VERSION,
        "implementation": {
            "package": "CoREMOF-tools",
            "package_version": __version__,
            "source_capture": "MODULE_IMPORT_BOUND",
            "example_script_sha256": _IMPORTED_EXAMPLE_SCRIPT_SHA256,
            "coremof_source_sha256": dict(
                sorted(coremof_source_hashes.items())
            ),
        },
        "dataset_version": dataset_version,
        "dataset_input_status": dataset_input_status,
        "parent_input_status": parent_input_status,
        "cif_files_verified": cif_files_verified,
        "checker_view": checker_view,
        "checker_view_kind": (
            "OFFICIAL_RELEASE_VIEW"
            if checker_view_official
            else "USER_DEFINED"
        ),
        "filters": filters,
        "ranking": {
            "column": rank_by,
            "kind": ranking_kind,
            "order": order,
            "tie_breaker": "structure_id_ascending",
            "numeric_policy": {
                "integer": "EXACT",
                "numeric_text": "EXACT_DECIMAL",
                "native_float": "PRESERVE_IEEE_754_VALUE",
                "mixed_comparison": "EXACT_DECIMAL_VALUE",
            },
            "null_policy": "EXCLUDE",
            "non_numeric_policy": "EXCLUDE",
            "nonfinite_policy": "EXCLUDE",
            "limit": limit,
        },
        "counts": {
            "release": len(dataset),
            "after_checker_and_identity_filters": len(filtered),
            "excluded_missing_required_target": len(required_excluded),
            "after_required_target_filter": len(required_view),
            "excluded_invalid_ranking_value": len(ranking_exclusions),
            "finite_rankable": finite_count,
            "emitted": len(rows),
        },
        "ranking_exclusion_counts": dict(sorted(exclusion_counts.items())),
        "eligible_before_ranking_ids_sha256": _ids_sha256(
            sorted(required_view.structure_ids)
        ),
        "emitted_ids_sha256": _ids_sha256(selected_ids),
        "input_hashes": dict(sorted(input_hashes.items())),
        "target_merge_receipt": (
            _json_safe(target_receipt) if target_receipt is not None else None
        ),
        "csv_output": {
            "required_target_columns": {
                target: _required_target_csv_column(target)
                for target in sorted(required_tuple)
            },
            "required_targets_are_flat_columns": True,
        },
        "policies": {
            "filter_precedes_ranking": True,
            "missing_values_imputed": False,
            "nonfinite_values_ranked": False,
            "ranking_changes_scientific_labels": False,
        },
    }
    result = ScreeningResult(
        rows=tuple(rows),
        selected_ids=selected_ids,
        receipt=receipt,
        classified_dataset=classified,
        filters=filters,
        _authority_factory_token=_SCREENING_RESULT_FACTORY_TOKEN,
    )
    return _register_screening_result(result)


def _screen_source_payload(classified: object) -> Tuple[object, ...]:
    dataset = getattr(classified, "dataset", None)
    if dataset is None:
        raise AuthorityStateError("screening source has no dataset")
    metadata_rows = getattr(dataset, "metadata_rows", None)
    records = getattr(classified, "records", None)
    if not isinstance(metadata_rows, (list, tuple)) or not isinstance(
        records, (list, tuple)
    ):
        raise AuthorityStateError("screening source rows must be sequences")
    if any(not isinstance(row, Mapping) for row in metadata_rows):
        raise AuthorityStateError("screening metadata rows must be mappings")
    record_payload = []
    for record in records:
        structure_id = _exact_nonblank_string(
            getattr(record, "structure_id", None), "classified structure_id"
        )
        metadata = getattr(record, "metadata", None)
        if not isinstance(metadata, Mapping):
            raise AuthorityStateError("classified record metadata must be a mapping")
        record_payload.append(
            (
                structure_id,
                getattr(record, "checker_view", None),
                getattr(record, "label", None),
                getattr(record, "checker_statuses", None),
                metadata,
            )
        )
    return (
        "coremof-generic-screen-source/1",
        getattr(classified, "checker_view", None),
        getattr(classified, "checker_view_official", False),
        tuple(record_payload),
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


def _screening_identity_snapshot(result: ScreeningResult) -> Tuple[object, ...]:
    return (
        id(result.rows),
        id(result.selected_ids),
        id(result.receipt),
        id(result.classified_dataset),
        id(result.filters),
        result._authority_factory_token is _SCREENING_RESULT_FACTORY_TOKEN,
        tuple(
            (
                id(row),
                row.get("structure_id"),
                row.get("rank"),
                type(row.get("ranking_value")).__name__,
                row.get("ranking_value"),
            )
            for row in result.rows
        ),
    )


def _screening_seal_fingerprint(result: ScreeningResult) -> str:
    """Bind all ranked rows, filters, and receipt claims without coercion."""

    return state_fingerprint(
        (
            "coremof-screening-result-generation/2",
            result.rows,
            result.selected_ids,
            result.receipt,
            result.filters,
        )
    )


def _register_screening_result(result: ScreeningResult) -> ScreeningResult:
    if result._authority_factory_token is not _SCREENING_RESULT_FACTORY_TOKEN:
        raise AuthorityStateError("cannot register a non-factory screening result")
    _verify_result_implementation(result.receipt)
    _validate_screening_result(result)
    source = result.classified_dataset
    official = _validate_classified_generation_if_present(source)
    source_fingerprint = None
    if not official:
        source_fingerprint = state_fingerprint(_screen_source_payload(source))
    fingerprint = _screening_seal_fingerprint(result)
    _SCREENING_RESULTS.register(
        result,
        fingerprint,
        {
            "source": source,
            "source_official": official,
            "source_fingerprint": source_fingerprint,
            "identity_snapshot": _screening_identity_snapshot(result),
        },
    )
    return result


def _require_screening_result(result: ScreeningResult) -> None:
    current = _SCREENING_RESULTS.entry(result)
    if current is None:
        raise AuthorityStateError(
            "screening result was not produced by screen_dataset()"
        )
    expected_fingerprint, context = current
    if not isinstance(context, Mapping):  # pragma: no cover
        raise AuthorityStateError("screening authority context is malformed")
    if _screening_identity_snapshot(result) != context.get("identity_snapshot"):
        raise AuthorityStateError(
            "screening result changed after construction; recompute it"
        )
    if _screening_seal_fingerprint(result) != expected_fingerprint:
        raise AuthorityStateError(
            "screening result fingerprint changed after construction; recompute it"
        )
    source = context.get("source")
    if context.get("source_official") is True:
        if not _validate_classified_generation_if_present(source):
            raise AuthorityStateError("official screening source lost authentication")
    else:
        try:
            source_fingerprint = state_fingerprint(_screen_source_payload(source))
        except (AttributeError, KeyError, TypeError, ValueError) as exc:
            raise AuthorityStateError(
                "generic screening source changed after result construction"
            ) from exc
        if source_fingerprint != context.get("source_fingerprint"):
            raise AuthorityStateError(
                "generic screening source changed after result construction"
            )
    _verify_result_implementation(result.receipt)
    _validate_screening_result(result)


def make_parent_aware_split(
    result: ScreeningResult,
    *,
    parent_method: str = "priority_main",
    leakage_guard: str = "auto",
    missing_parent: str = "singleton",
    fractions: Sequence[float] = (0.8, 0.1, 0.1),
    random_state: int = 42,
    stratify_by: Sequence[str] = ("label",),
):
    """Split ranked candidates with explicitly defined parent and leakage policies.

    The default ``priority_main`` is a conflict-aware hierarchy over exact RAC5,
    then complete MOFid-v2, then complete MOFid-v1 release triads.  A lower group
    never merges two stronger components: it records ``PARENT_METHOD_CONFLICT``.
    Missing evidence becomes a singleton unless exclusion is requested.
    ``priority_main`` does not use and explicitly excludes Zeo++, topology,
    source IDs, common names, CIF hashes, provisional source-ID/MOFid transitive
    groups, and StructureMatcher evidence from that explanatory hierarchy.  Here
    priority means parent-evidence precedence; it does not rank, schedule, or
    recalculate failed scientific features.

    ``leakage_guard="auto"`` selects ``main_union`` for ``priority_main`` and
    ``parent_only`` otherwise.  ``main_union`` is a full-release transitive
    leakage block over exact CIF SHA-256, source sibling, RAC5, MOFid-v2, and
    MOFid-v1 relations before candidate filtering; it is not an explanatory
    parent claim or proof of structural identity. ``parent_only`` uses only the
    selected parent method's groups as split blocks.
    """

    _require_screening_result(result)
    filters = result.filters
    return result.classified_dataset.train_valid_test_split(
        parent_method=parent_method,
        leakage_guard=leakage_guard,
        missing_parent=missing_parent,
        fractions=fractions,
        random_state=random_state,
        stratify_by=stratify_by,
        labels=filters["labels"],
        sources=filters["sources"],
        variants=filters["variants"],
        metals=filters["metals"],
        structure_ids=result.selected_ids,
        required_targets=(
            filters["required_targets"] or None
        ),
        required_target_mode=filters["required_target_mode"],
    )


def _csv_value(value: object) -> object:
    if isinstance(value, Decimal):
        if not value.is_finite():
            raise ScreeningError("refusing to write a non-finite decimal value")
        return str(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ScreeningError("refusing to write a non-finite CSV value")
        return format(value, ".17g")
    if isinstance(value, (Mapping, list, tuple)):
        return json.dumps(
            _json_safe(value), sort_keys=True, separators=(",", ":"), allow_nan=False
        )
    if value is None:
        return ""
    return value


def _validate_stem(stem: str) -> str:
    if (
        not isinstance(stem, str)
        or not stem
        or stem in (".", "..")
        or Path(stem).name != stem
    ):
        raise ScreeningError("stem must be a filename stem without directories")
    return stem


def _publish_bundle(
    staged_to_final: Sequence[Tuple[Path, Path]], overwrite: bool
) -> None:
    """Publish a staged bundle with race-safe create-if-absent rollback.

    ``overwrite=False`` is safe against another writer replacing a pathname
    during publication: rollback removes only a hard link still pointing to
    this call's staged inode.  ``overwrite=True`` is intentionally a
    single-writer replacement operation; callers must serialize concurrent
    overwrite publishers for the same output stem.
    """

    final_paths = tuple(final for _, final in staged_to_final)
    staged_by_final = {final: staged for staged, final in staged_to_final}
    if not overwrite:
        existing = [path for path in final_paths if path.exists()]
        if existing:
            raise FileExistsError(
                "refusing to overwrite existing output: {}".format(
                    ", ".join(str(path) for path in existing)
                )
            )
    staging = staged_to_final[0][0].parent
    backups: Dict[Path, Path] = {}
    published_identities: Dict[Path, Path] = {}
    published = []

    def quarantine_published_generation(final: Path, identity: Path) -> None:
        """Detach and remove only the generation identified by ``identity``."""

        quarantine = staging / (final.name + ".rollback")
        try:
            os.replace(str(final), str(quarantine))
        except FileNotFoundError:
            return
        try:
            published_by_this_call = identity.exists() and os.path.samefile(
                str(identity), str(quarantine)
            )
        except FileNotFoundError:
            published_by_this_call = False
        if published_by_this_call:
            quarantine.unlink()
            return

        try:
            os.link(str(quarantine), str(final))
        except OSError:
            # A newer writer may already own ``final``.  Preserve the detached
            # foreign generation outside the soon-to-be-removed staging tree;
            # losing it would be worse than leaving a clearly named recovery
            # artifact for manual reconciliation.
            descriptor, preserved_name = tempfile.mkstemp(
                prefix=".{}-rollback-preserved-".format(final.name),
                dir=str(final.parent),
            )
            os.close(descriptor)
            os.replace(str(quarantine), preserved_name)
            return
        quarantine.unlink()

    try:
        if overwrite:
            for final in final_paths:
                if final.exists():
                    backup = staging / (final.name + ".previous")
                    try:
                        os.link(str(final), str(backup))
                    except OSError:
                        shutil.copy2(str(final), str(backup))
                    backups[final] = backup
        for staged, final in staged_to_final:
            if overwrite:
                if final not in backups:
                    identity = staging / (final.name + ".published")
                    os.link(str(staged), str(identity))
                    published_identities[final] = identity
                os.replace(str(staged), str(final))
            else:
                # A hard link is an atomic create-if-absent operation.  It
                # closes the race between the preflight existence check and
                # publication without needing a platform-specific lock.
                os.link(str(staged), str(final))
            published.append(final)
        directory_fd = os.open(str(final_paths[0].parent), os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    except BaseException:
        for final in reversed(final_paths):
            backup = backups.get(final)
            if backup is not None and backup.exists():
                os.replace(str(backup), str(final))
            elif final in published:
                quarantine_published_generation(
                    final,
                    published_identities.get(final, staged_by_final[final]),
                )
        raise


def _validate_split_correspondence(
    result: ScreeningResult,
    split_result: object,
    split_receipt: Mapping[str, object],
) -> None:
    """Require an optional split to describe exactly this screening universe."""

    receipt = result.receipt
    if split_receipt.get("dataset_version") != receipt.get("dataset_version"):
        raise ScreeningError("split dataset version differs from screening result")
    if split_receipt.get("checker_view") != receipt.get("checker_view"):
        raise ScreeningError("split checker view differs from screening result")
    if _json_safe(split_receipt.get("input_hashes")) != _json_safe(
        receipt.get("input_hashes")
    ):
        raise ScreeningError("split inputs differ from screening result inputs")

    split_filters = split_receipt.get("filters")
    if not isinstance(split_filters, Mapping):
        raise ScreeningError("split receipt has no filter mapping")
    expected_filters = result.filters
    for key in ("labels", "sources", "variants", "metals"):
        if _json_safe(split_filters.get(key)) != _json_safe(expected_filters.get(key)):
            raise ScreeningError(
                "split {} filter differs from screening result".format(key)
            )
    if _json_safe(split_filters.get("structure_ids")) != list(result.selected_ids):
        raise ScreeningError(
            "split structure-ID selection differs from ranked screening rows"
        )
    target_filters = split_filters.get("targets")
    required = list(expected_filters.get("required_targets") or ())
    if target_filters is None:
        if required:
            raise ScreeningError("split omits required-target filtering")
    elif not isinstance(target_filters, Mapping):
        raise ScreeningError("split target filters must be a mapping")
    elif (
        _json_safe(target_filters.get("required")) != required
        or target_filters.get("mode") != expected_filters.get("required_target_mode")
        or target_filters.get("filter_precedes_assignment") is not True
        or target_filters.get("leakage_blocks_use_full_release_universe") is not True
    ):
        raise ScreeningError("split target filters differ from screening result")

    classified_ids = tuple(
        _exact_nonblank_string(
            getattr(record, "structure_id", None), "classified structure_id"
        )
        for record in result.classified_dataset
    )
    preselection = split_filters.get("preselection")
    if not isinstance(preselection, Mapping):
        raise ScreeningError("split receipt has no full-universe preselection")
    if (
        preselection.get("release_count") != len(classified_ids)
        or preselection.get("selected_count") != len(classified_ids)
        or preselection.get("structure_ids_sha256")
        != _ids_sha256(sorted(classified_ids))
    ):
        raise ScreeningError("split full-release universe differs from screening result")

    try:
        assignment_ids = {
            _exact_nonblank_string(value, "split assignment structure_id")
            for value in split_result.assignments
        }
        exclusion_ids = {
            _exact_nonblank_string(value, "split exclusion structure_id")
            for value in split_result.exclusions
        }
        index_ids = {
            _exact_nonblank_string(value, "split index structure_id")
            for value in split_result.index_by_id
        }
    except (AttributeError, TypeError) as exc:
        raise ScreeningError("split result lacks auditable assignment mappings") from exc
    release_ids = set(classified_ids)
    selected_ids = set(result.selected_ids)
    if index_ids != release_ids or assignment_ids.union(exclusion_ids) != release_ids:
        raise ScreeningError("split assignment universe differs from screening result")
    if assignment_ids.difference(selected_ids):
        raise ScreeningError("split assigns structures absent from ranked screening rows")
    if not selected_ids.issubset(assignment_ids.union(exclusion_ids)):
        raise ScreeningError("split omits ranked screening structures")


def _validate_screening_result(result: ScreeningResult) -> None:
    """Reverse-check ranked rows against their immutable screening receipt."""

    rows = tuple(result.rows)
    selected_ids = tuple(result.selected_ids)
    for structure_id in selected_ids:
        _exact_nonblank_string(structure_id, "selected structure_id")
    row_ids = tuple(
        _exact_nonblank_string(row.get("structure_id"), "screening row structure_id")
        for row in rows
    )
    if row_ids != selected_ids:
        raise ScreeningError(
            "screening rows and selected_ids differ; rebuild the screening result"
        )
    if len(set(selected_ids)) != len(selected_ids) or any(not value for value in selected_ids):
        raise ScreeningError(
            "screening result contains duplicate or empty selected structure IDs"
        )

    receipt = result.receipt
    expected_receipt_keys = {
        "schema_version",
        "implementation",
        "dataset_version",
        "dataset_input_status",
        "parent_input_status",
        "cif_files_verified",
        "checker_view",
        "checker_view_kind",
        "filters",
        "ranking",
        "counts",
        "ranking_exclusion_counts",
        "eligible_before_ranking_ids_sha256",
        "emitted_ids_sha256",
        "input_hashes",
        "target_merge_receipt",
        "csv_output",
        "policies",
    }
    if set(receipt) != expected_receipt_keys:
        raise ScreeningError("screening receipt has an invalid top-level contract")
    if receipt.get("schema_version") != SCHEMA_VERSION:
        raise ScreeningError("screening receipt has an invalid schema version")
    receipt_dataset_version = _exact_nonblank_string(
        receipt.get("dataset_version"), "receipt dataset_version"
    )
    receipt_checker_view = _exact_nonblank_string(
        receipt.get("checker_view"), "receipt checker_view"
    )
    _optional_exact_string(
        receipt.get("dataset_input_status"), "receipt dataset_input_status"
    )
    _optional_exact_string(
        receipt.get("parent_input_status"), "receipt parent_input_status"
    )
    implementation = receipt.get("implementation")
    if not isinstance(implementation, Mapping) or set(implementation) != {
        "package",
        "package_version",
        "source_capture",
        "example_script_sha256",
        "coremof_source_sha256",
    }:
        raise ScreeningError("screening receipt has an invalid implementation contract")
    if (
        implementation.get("package") != "CoREMOF-tools"
        or implementation.get("package_version") != __version__
        or implementation.get("source_capture") != "MODULE_IMPORT_BOUND"
    ):
        raise ScreeningError("screening receipt implementation claims are invalid")
    counts = receipt.get("counts")
    ranking = receipt.get("ranking")
    receipt_filters = receipt.get("filters")
    if not isinstance(counts, Mapping) or not isinstance(ranking, Mapping):
        raise ScreeningError("screening receipt is missing counts or ranking metadata")
    if _json_safe(receipt_filters) != _json_safe(result.filters):
        raise ScreeningError(
            "screening filters differ from their receipt; rebuild the screening result"
        )
    if counts.get("emitted") != len(rows):
        raise ScreeningError("screening receipt emitted count differs from ranked rows")
    if receipt.get("emitted_ids_sha256") != _ids_sha256(selected_ids):
        raise ScreeningError("screening receipt emitted-ID digest differs from ranked rows")

    ranking_column = ranking.get("column")
    order = ranking.get("order")
    if not isinstance(ranking_column, str) or order not in ("ascending", "descending"):
        raise ScreeningError("screening receipt has an invalid ranking contract")
    parsed = []
    for expected_rank, row in enumerate(rows, start=1):
        if row.get("rank") != expected_rank:
            raise ScreeningError("screening rows have non-contiguous ranks")
        if row.get("ranking_column") != ranking_column:
            raise ScreeningError("screening row ranking column differs from its receipt")
        if row.get("checker_view") != receipt.get("checker_view"):
            raise ScreeningError("screening row checker view differs from its receipt")
        value, reason = _finite_number(row.get("ranking_value"))
        if reason is not None or value is None:
            raise ScreeningError("screening result contains an invalid ranking value")
        parsed.append(
            (
                _exact_nonblank_string(
                    row["structure_id"], "screening row structure_id"
                ),
                value,
            )
        )

    expected_order = sorted(parsed, key=lambda item: item[0])
    expected_order.sort(
        key=lambda item: _comparison_decimal(item[1]),
        reverse=(order == "descending"),
    )
    if tuple(structure_id for structure_id, _ in expected_order) != selected_ids:
        raise ScreeningError("screening rows do not follow their declared ranking order")

    integer_count_keys = (
        "release",
        "after_checker_and_identity_filters",
        "excluded_missing_required_target",
        "after_required_target_filter",
        "excluded_invalid_ranking_value",
        "finite_rankable",
        "emitted",
    )
    if any(
        not isinstance(counts.get(key), int)
        or isinstance(counts.get(key), bool)
        or counts.get(key) < 0
        for key in integer_count_keys
    ):
        raise ScreeningError("screening receipt contains invalid count values")
    if (
        counts["after_checker_and_identity_filters"]
        - counts["excluded_missing_required_target"]
        != counts["after_required_target_filter"]
        or counts["after_required_target_filter"]
        - counts["excluded_invalid_ranking_value"]
        != counts["finite_rankable"]
        or counts["emitted"] > counts["finite_rankable"]
        or counts["after_checker_and_identity_filters"] > counts["release"]
    ):
        raise ScreeningError("screening receipt counts are internally inconsistent")
    limit = ranking.get("limit")
    expected_emitted = (
        min(counts["finite_rankable"], limit)
        if isinstance(limit, int) and not isinstance(limit, bool)
        else counts["finite_rankable"]
    )
    if counts["emitted"] != expected_emitted:
        raise ScreeningError("screening receipt limit differs from emitted row count")

    expected_filter_keys = {
        "labels",
        "sources",
        "variants",
        "metals",
        "required_targets",
        "required_target_mode",
    }
    if set(result.filters) != expected_filter_keys:
        raise ScreeningError("screening result contains an invalid filter contract")
    if result.filters["required_target_mode"] not in ("all", "any"):
        raise ScreeningError("screening result has an invalid required-target mode")
    if limit is not None and (
        isinstance(limit, bool) or not isinstance(limit, int) or limit <= 0
    ):
        raise ScreeningError("screening receipt has an invalid ranking limit")

    classified = result.classified_dataset
    try:
        classified_records = tuple(classified)
        classified_by_id = {
            _exact_nonblank_string(
                getattr(record, "structure_id", None),
                "classified structure_id",
            ): record
            for record in classified_records
        }
        dataset = classified.dataset
    except (AttributeError, TypeError) as exc:
        raise ScreeningError("screening result lacks an auditable classified dataset") from exc
    # Validate every external receipt-claim source before comparing derived
    # fields, so a malformed container or boolean cannot be masked by an
    # earlier secondary mismatch.
    input_hashes = _input_hash_mapping(dataset)
    (
        expected_parent_status,
        expected_dataset_status,
        _expected_provisional,
    ) = _screen_release_receipt_state(dataset)
    expected_cif_files_verified = _boolean_attribute(
        dataset, "cif_files_verified", False
    )
    checker_view_official = _boolean_attribute(
        classified,
        "checker_view_official",
        False,
    )
    if len(classified_by_id) != len(classified_records):
        raise ScreeningError("classified dataset contains duplicate structure IDs")
    if len(classified_records) != counts["release"]:
        raise ScreeningError("screening release count differs from classified dataset")
    classified_checker_view = _exact_nonblank_string(
        getattr(classified, "checker_view", None), "classified checker_view"
    )
    if checker_view_official and classified_checker_view not in CHECKER_PRESETS:
        raise ScreeningError(
            "checker_view_official=True requires a canonical official checker preset"
        )
    if checker_view_official and not _is_authenticated_official_checker_view(
        classified
    ):
        raise ScreeningError(
            "checker_view_official=True requires an internally authenticated "
            "recomputed release view"
        )
    if classified_checker_view != receipt_checker_view:
        raise ScreeningError("screening checker view differs from classified dataset")
    classified_dataset_version = _exact_nonblank_string(
        getattr(classified, "dataset_version", None),
        "classified dataset_version",
    )
    if classified_dataset_version != receipt_dataset_version:
        raise ScreeningError("screening dataset version differs from classified dataset")
    if _json_safe(input_hashes) != _json_safe(
        receipt.get("input_hashes")
    ):
        raise ScreeningError("screening input hashes differ from classified dataset")
    if _json_safe(getattr(dataset, "target_input_receipt", None)) != _json_safe(
        receipt.get("target_merge_receipt")
    ):
        raise ScreeningError("screening target receipt differs from classified dataset")
    if receipt.get("dataset_input_status") != expected_dataset_status:
        raise ScreeningError("screening dataset status differs from classified dataset")
    if receipt.get("parent_input_status") != expected_parent_status:
        raise ScreeningError("screening parent status differs from classified dataset")
    if receipt.get("cif_files_verified") is not expected_cif_files_verified:
        raise ScreeningError("screening CIF-verification claim differs from dataset")
    expected_view_kind = (
        "OFFICIAL_RELEASE_VIEW"
        if checker_view_official
        else "USER_DEFINED"
    )
    if receipt.get("checker_view_kind") != expected_view_kind:
        raise ScreeningError("screening checker-view kind differs from classified dataset")

    target_columns = _column_name_tuple(
        getattr(dataset, "target_columns", None), "target_columns"
    )
    feature_columns = _column_name_tuple(
        getattr(dataset, "feature_columns", None), "feature_columns"
    )
    expected_ranking_kind = (
        "TARGET"
        if ranking_column in target_columns
        else "JOINED_FEATURE"
        if ranking_column in feature_columns
        else "RELEASE_METADATA"
    )
    expected_ranking_contract = {
        "column": ranking_column,
        "kind": expected_ranking_kind,
        "order": order,
        "tie_breaker": "structure_id_ascending",
        "numeric_policy": {
            "integer": "EXACT",
            "numeric_text": "EXACT_DECIMAL",
            "native_float": "PRESERVE_IEEE_754_VALUE",
            "mixed_comparison": "EXACT_DECIMAL_VALUE",
        },
        "null_policy": "EXCLUDE",
        "non_numeric_policy": "EXCLUDE",
        "nonfinite_policy": "EXCLUDE",
        "limit": limit,
    }
    if _json_safe(ranking) != expected_ranking_contract:
        raise ScreeningError("screening ranking contract differs from dataset replay")
    expected_csv_contract = {
        "required_target_columns": {
            target: _required_target_csv_column(target)
            for target in sorted(result.filters["required_targets"] or ())
        },
        "required_targets_are_flat_columns": True,
    }
    if _json_safe(receipt.get("csv_output")) != expected_csv_contract:
        raise ScreeningError("screening CSV contract differs from requested targets")
    if _json_safe(receipt.get("policies")) != {
        "filter_precedes_ranking": True,
        "missing_values_imputed": False,
        "nonfinite_values_ranked": False,
        "ranking_changes_scientific_labels": False,
    }:
        raise ScreeningError("screening scientific policies are invalid")

    try:
        filtered = classified.filter(
            labels=result.filters["labels"],
            sources=result.filters["sources"],
            variants=result.filters["variants"],
            metals=result.filters["metals"],
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise ScreeningError("screening filters cannot be replayed") from exc
    required_targets = tuple(result.filters["required_targets"] or ())
    required_records = tuple(
        record
        for record in filtered
        if not required_targets
        or _passes_required_targets(
            record.metadata,
            required_targets,
            str(result.filters["required_target_mode"]),
        )
    )
    numeric_records = []
    replay_exclusions: Dict[str, str] = {}
    for record in required_records:
        replay_value, replay_reason = _finite_number(
            record.metadata.get(ranking_column)
        )
        if replay_reason is not None:
            replay_exclusions[record.structure_id] = replay_reason
            continue
        if replay_value is None:  # pragma: no cover - finite parser invariant
            raise ScreeningError("ranking replay produced no finite value")
        numeric_records.append((record, replay_value))
    numeric_records.sort(key=lambda item: item[0].structure_id)
    numeric_records.sort(
        key=lambda item: _comparison_decimal(item[1]),
        reverse=(order == "descending"),
    )
    finite_count = len(numeric_records)
    emitted_records = numeric_records[:limit] if limit is not None else numeric_records
    replay_ids = tuple(record.structure_id for record, _ in emitted_records)
    if replay_ids != selected_ids:
        raise ScreeningError(
            "ranked screening rows differ from a replay of the classified dataset"
        )

    replay_counts = {
        "release": len(classified_records),
        "after_checker_and_identity_filters": len(filtered),
        "excluded_missing_required_target": len(filtered) - len(required_records),
        "after_required_target_filter": len(required_records),
        "excluded_invalid_ranking_value": len(replay_exclusions),
        "finite_rankable": finite_count,
        "emitted": len(emitted_records),
    }
    if dict(counts) != replay_counts:
        raise ScreeningError("screening receipt counts differ from dataset replay")
    if _json_safe(receipt.get("ranking_exclusion_counts")) != dict(
        sorted(Counter(replay_exclusions.values()).items())
    ):
        raise ScreeningError("ranking exclusion counts differ from dataset replay")
    if receipt.get("eligible_before_ranking_ids_sha256") != _ids_sha256(
        sorted(record.structure_id for record in required_records)
    ):
        raise ScreeningError("eligible-ID digest differs from dataset replay")

    expected_row_keys = {
        "rank",
        "structure_id",
        "ranking_column",
        "ranking_value",
        "checker_view",
        "classification_label",
        "source_database",
        "source_id",
        "structure_variant",
        "metal_elements",
        "cif_file",
    }.union(_required_target_csv_column(target) for target in required_targets)
    for expected_rank, ((record, expected_value), row) in enumerate(
        zip(emitted_records, rows), start=1
    ):
        expected_row = {
            "rank": expected_rank,
            "structure_id": record.structure_id,
            "ranking_column": ranking_column,
            "ranking_value": expected_value,
            "checker_view": classified.checker_view,
            "classification_label": record.label,
            "source_database": record.source_database,
            "source_id": record.metadata.get("source_id", ""),
            "structure_variant": record.structure_variant,
            "metal_elements": ";".join(record.metal_elements),
            "cif_file": record.metadata.get("cif_file", ""),
        }
        for target in required_targets:
            expected_row[_required_target_csv_column(target)] = record.metadata.get(
                target
            )
        if set(row) != expected_row_keys:
            raise ScreeningError("screening row has an invalid public column contract")
        for key, expected_value in expected_row.items():
            actual_value = row.get(key)
            if type(actual_value) is not type(expected_value) or actual_value != expected_value:
                raise ScreeningError(
                    "screening row {!r} differs from classified dataset field {!r}".format(
                        record.structure_id, key
                    )
                )


def write_screening(
    result: ScreeningResult,
    output_directory: Path,
    *,
    stem: str = "coremof_screening",
    split_result=None,
    overwrite: bool = False,
) -> Mapping[str, Path]:
    """Write the ranked table, receipt, and optional split as one bundle.

    ``overwrite=True`` has explicit single-writer semantics for a given output
    stem.  The default create-if-absent mode protects another writer's
    replacement inode during rollback.
    """

    if type(overwrite) is not bool:
        raise TypeError("overwrite must be a boolean")
    _require_screening_result(result)
    split_receipt = None
    split_hashes = None
    if split_result is not None:
        try:
            split_receipt = split_result.receipt()
        except (AttributeError, TypeError, ValueError) as exc:
            raise ScreeningError("split result has no valid receipt") from exc
        if not isinstance(split_receipt, Mapping):
            raise ScreeningError("split receipt must be a mapping")
        _validate_split_correspondence(result, split_result, split_receipt)
        split_hashes = _split_source_hashes(result.receipt, split_receipt)
    stem = _validate_stem(stem)
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    final = {
        "ranked_csv": output_directory / (stem + ".csv"),
        "receipt_json": output_directory / (stem + ".receipt.json"),
    }
    if split_result is not None:
        final.update(
            {
                "split_csv": output_directory / (stem + ".split.csv"),
                "split_json": output_directory / (stem + ".split.json"),
            }
        )
    if not overwrite:
        existing = [path for path in final.values() if path.exists()]
        if existing:
            raise FileExistsError(
                "refusing to overwrite existing output: {}".format(
                    ", ".join(str(path) for path in existing)
                )
            )

    staging = Path(
        tempfile.mkdtemp(prefix=".{}-".format(stem), dir=str(output_directory))
    )
    staged = {name: staging / path.name for name, path in final.items()}
    try:
        output_receipt = dict(_json_safe(result.receipt))
        if split_result is not None:
            if split_receipt is None or split_hashes is None:  # pragma: no cover
                raise ScreeningError("split validation did not produce a receipt")
            output_receipt["implementation"]["coremof_source_sha256"] = dict(
                sorted(split_hashes.items())
            )

        required_target_fields = tuple(
            _required_target_csv_column(target)
            for target in sorted(result.filters["required_targets"])
        )
        fields = (
            "rank",
            "structure_id",
            "ranking_column",
            "ranking_value",
            "checker_view",
            "classification_label",
            "source_database",
            "source_id",
            "structure_variant",
            "metal_elements",
            "cif_file",
        ) + required_target_fields
        with staged["ranked_csv"].open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for row in result.rows:
                writer.writerow({key: _csv_value(row.get(key)) for key in fields})
            handle.flush()
            os.fsync(handle.fileno())

        if split_result is not None:
            split_result.to_csv(staged["split_csv"])
            split_result.to_json(staged["split_json"])

        output_receipt["outputs"] = {
            "ranked_csv": {
                "file": final["ranked_csv"].name,
                "sha256": _sha256_file(staged["ranked_csv"]),
                "row_count": len(result.rows),
            }
        }
        if split_receipt is not None:
            output_receipt["split"] = {
                "enabled": True,
                "contract_definitions": split_receipt["contract_definitions"],
                "assignment_csv": final["split_csv"].name,
                "receipt_json": final["split_json"].name,
                "assignment_csv_sha256": _sha256_file(staged["split_csv"]),
                "receipt_json_sha256": _sha256_file(staged["split_json"]),
                "assignment_sha256": split_receipt["assignment_sha256"],
                "parent_method": split_receipt["parent_method"],
                "parent_method_definition": split_receipt[
                    "parent_method_definition"
                ],
                "requested_leakage_guard": split_receipt[
                    "requested_leakage_guard"
                ],
                "requested_leakage_guard_definition": split_receipt[
                    "requested_leakage_guard_definition"
                ],
                "leakage_guard": split_receipt["leakage_guard"],
                "leakage_guard_definition": split_receipt[
                    "leakage_guard_definition"
                ],
                "counts": split_receipt["counts"],
                "leakage_audit": split_receipt["leakage_audit"],
                "official_split": split_receipt["official_split"],
                "provisional_input": split_receipt["provisional_input"],
                "implementation_source_sha256": dict(
                    sorted(split_hashes.items())
                ),
            }
        else:
            output_receipt["split"] = {"enabled": False}
        _verify_result_implementation(result.receipt)
        with staged["receipt_json"].open("w", encoding="utf-8") as handle:
            json.dump(
                output_receipt, handle, indent=2, sort_keys=True, allow_nan=False
            )
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())

        ordered_names = ["ranked_csv"]
        if split_result is not None:
            ordered_names.extend(("split_csv", "split_json"))
        ordered_names.append("receipt_json")
        _verify_result_implementation(result.receipt)
        _publish_bundle(
            [(staged[name], final[name]) for name in ordered_names], overwrite
        )
    finally:
        shutil.rmtree(staging, ignore_errors=True)
    return final


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Filter CoRE-MOF rows before ranking a finite numeric target or "
            "feature, with an optional parent-aware split whose project-defined "
            "terms are expanded in its receipts"
        )
    )
    parser.add_argument("release", type=Path, help="extracted CoRE-MOF release root")
    parser.add_argument("--target-config", type=Path)
    parser.add_argument("--rank-by", required=True)
    parser.add_argument(
        "--order", choices=("ascending", "descending"), default="descending"
    )
    parser.add_argument(
        "--checkers",
        choices=("3checker", "4checker", "5checker"),
        default="5checker",
    )
    parser.add_argument(
        "--label",
        action="append",
        dest="labels",
        help="classification label to retain; repeat as needed (default: CR)",
    )
    parser.add_argument("--source", action="append", dest="sources")
    parser.add_argument("--variant", action="append", dest="variants")
    parser.add_argument("--metal", action="append", dest="metals")
    parser.add_argument(
        "--require-target", action="append", dest="required_targets"
    )
    parser.add_argument(
        "--required-target-mode", choices=("all", "any"), default="all"
    )
    parser.add_argument("--limit", type=int)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--stem", default="coremof_screening")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--verify-cifs", action="store_true")

    parser.add_argument(
        "--split",
        action="store_true",
        help="also create parent-aware train/validation/test assignments",
    )
    parser.add_argument(
        "--parent-method",
        choices=SELECTABLE_PARENT_METHODS,
        default="priority_main",
        help=(
            "explanatory relation for --split. Each non-control choice reads a "
            "validated status/group/size triad: MATCHED means available size>=2, "
            "UNMATCHED available size=1, and NOT_AVAILABLE no edge; nulls never match. "
            "priority_main is the project-defined full-release hierarchy over rac, "
            "mofid2, and mofid1 triads: RAC5 anchors, then MOFid v2, then MOFid v1; a lower "
            "group attaches unresolved rows to at most one stronger component, never "
            "merges multiple stronger components, and records PARENT_METHOD_CONFLICT. "
            "Missing evidence is a unique singleton or explicit exclusion; Zeo++, "
            "topology, source ID, common name, CIF hash, provisional "
            "source-ID/MOFid transitive groups, and "
            "StructureMatcher are excluded. Here priority means parent-evidence "
            "precedence; it does not rank, schedule, or recalculate failed scientific "
            "features. rac5 uses all 264 finite binary64 values "
            "after -0.0 to +0.0 and float.hex with rtol=atol=0; mofid_v2/mofid_v1 "
            "compare complete canonicalized strings. rac5_zeo combines RAC5 and zeo; "
            "zeo uses exact float.hex equality for 13 intensive N2/He fields (radii "
            "1.655/1.32 A), equal N2 channel dimension, and available equal framework "
            "periodic dimension. source_id uses the namespaced database/ID pair, "
            "common_name exact canonicalized text, and none one singleton per row. "
            "identity_union selects the separate project-defined provisional "
            "source-ID/MOFid transitive groups read from identity "
            "status/group/size: "
            "each named release is freshly recomputed over all its current rows by "
            "transitively unioning exact database-namespaced source ID plus eligible "
            "complete MOFid-v2 or MOFid-v1 text edges with no precedence and without "
            "importing an earlier component. Each group and identity_size count "
            "one transitive connected component of structures, not edges or identifiers. "
            "Missing identifiers and unsuccessful MOFid statuses add no edge. "
            "It is not proof of structural identity and does not enter main_union. "
            "Canonicalized text means convert to text, collapse Unicode whitespace, "
            "trim, reject case-insensitive empty, -, nan, none, null, n/a, na, "
            "unknown, missing, timeout, timed out, error, failed, fail, fail process, "
            "failed process, or process failed whole fields, then apply Unicode NFKC and casefold, "
            "so those declared whole-field placeholders never compare; "
            "with no fuzzy match. This text processing only compares identifiers and "
            "does not modify a CIF, atom, bond, occupancy, coordinate, chemistry, or "
            "unit cell. Eligible MOFid statuses are SUCCESS, "
            "SUCCESS_TOPOLOGY_UNKNOWN, SUCCESS_TOPOLOGY_ERROR, and "
            "SUCCESS_TOPOLOGY_TIMEOUT; unresolved reconciliation, ambiguous node, "
            "timeout, error, no-MOF, unmatched-node, and decomposition-error values "
            "never match through a shared null. "
            "rac5_crystalnets (RT-) reads rac_crystalnets_status/group/size and "
            "requires exact equality of all 264 finite RAC5 values plus the complete "
            "successful current CrystalNets fingerprint; "
            "mofid_v2_crystalnets (M2T-) reads mofid2_crystalnets_status/group/size. "
            "They require exact complete RAC5/MOFid plus successful "
            "release-authorized current CrystalNets network/count/net/agreement and "
            "every subnet node status/dimension/key/name/genome field; counts equal "
            "subnet count and valid heterogeneous summary nulls are retained. Missing, "
            "nonfinite, partial, timed-out, failed, or otherwise incomplete input adds "
            "no edge. M2T is provisional whenever MOFid-v2 is and must be rebuilt "
            "before use if release-authorized MOFid-v2 values change. Its eligible statuses "
            "are exactly SUCCESS, SUCCESS_TOPOLOGY_UNKNOWN, SUCCESS_TOPOLOGY_ERROR, "
            "and SUCCESS_TOPOLOGY_TIMEOUT; the latter two are successful calculated "
            "identifiers with embedded topology qualifiers, while every other status "
            "adds no edge. "
            "structure_matcher_strict (SM-) is a "
            "Python 3.9/pymatgen 2024.2.8/NumPy 1.26.4 component view of pairs "
            "exhaustively blocked by parsed ElementComparator "
            "fractional composition and tested in both directions with "
            "pymatgen ElementComparator fits with ltol=stol=0.001, angle_tol=0.01, "
            "primitive_cell/attempt_supercell true, scale/allow_subset false, "
            "supercell_size=num_sites, no ignored species, and symmetric=true. Its "
            "parser expands symmetry, uses site/frac tolerances 1e-4, checks occupancy, "
            "sorts, and preserves disorder without manual repair, occupancy selection, "
            "atom deletion, or chemistry editing; parser, timeout, OOM, matcher, and "
            "asymmetric cases are NOT_AVAILABLE rather than unmatched. Directional "
            "displacement/(V/Nsites)^(1/3) is dimensionless, not angstrom RMSD; SM is "
            "a convenience component view rather than duplicate proof or an all-pairs "
            "claim, and direct edges are authoritative. RT/M2T/SM "
            "prefix digests are criterion-bound length-delimited UTF-8 SHA-256 with "
            "at least eight uppercase hex characters, extended only on collision. Receipts "
            "store exact per-method semantics (default: priority_main)"
        ),
    )
    parser.add_argument(
        "--leakage-guard",
        choices=LEAKAGE_GUARD_CHOICES,
        default="auto",
        help=(
            "split-block policy. main_union is the project-defined leakage guard, not "
            "an explanatory parent relation or proof of structural identity. It is a "
            "transitive union "
            "over the complete release of full manifest CIF SHA-256 and parent-table "
            "source-sibling/RAC5/MOFid-v2/MOFid-v1 relations before filters; missing "
            "CIF hashes fail, missing optional evidence adds no edge, and all edges "
            "have equal union status. Loaded releases require the source triad; only a "
            "manual object lacking it uses strip/upper database plus strip/casefold/"
            "whitespace-collapse source-ID fallback without NFKC. parent_only is the project-defined narrow guard "
            "that uses only selected explanatory groups and adds no cross-method edge. "
            "auto is only a selector: main_union for priority_main, parent_only "
            "otherwise. A priority conflict is excluded under parent_only but may remain "
            "assigned and diagnosed under main_union (default: auto)"
        ),
    )
    parser.add_argument(
        "--missing-parent",
        choices=("singleton", "exclude"),
        default="singleton",
        help=(
            "unavailable parent evidence: singleton creates one unique "
            "SINGLETON:<structure_id>; exclude records MISSING_PARENT_EVIDENCE and "
            "assigns no partition (default: singleton)"
        ),
    )
    parser.add_argument(
        "--fractions", type=float, nargs=3, default=(0.8, 0.1, 0.1)
    )
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--stratify-by", action="append", dest="stratify_by")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    labels = tuple(args.labels) if args.labels else ("CR",)
    try:
        dataset = load_dataset(
            args.release,
            target_config=args.target_config,
            verify_cif_files=args.verify_cifs,
        )
        result = screen_dataset(
            dataset,
            rank_by=args.rank_by,
            order=args.order,
            checkers=args.checkers,
            labels=labels,
            sources=args.sources,
            variants=args.variants,
            metals=args.metals,
            required_targets=args.required_targets,
            required_target_mode=args.required_target_mode,
            limit=args.limit,
        )
        split_result = None
        if args.split:
            split_result = make_parent_aware_split(
                result,
                parent_method=args.parent_method,
                leakage_guard=args.leakage_guard,
                missing_parent=args.missing_parent,
                fractions=tuple(args.fractions),
                random_state=args.random_state,
                stratify_by=(
                    tuple(args.stratify_by) if args.stratify_by else ("label",)
                ),
            )
        outputs = write_screening(
            result,
            args.output_directory,
            stem=args.stem,
            split_result=split_result,
            overwrite=args.overwrite,
        )
    except (
        FileExistsError,
        FileNotFoundError,
        KeyError,
        OSError,
        ReleaseValidationError,
        ScreeningError,
        TargetDataError,
        ValueError,
    ) as error:
        print("screen_candidates: error: {}".format(error), file=sys.stderr)
        return 2

    summary = {
        "dataset_version": result.receipt["dataset_version"],
        "checker_view": result.receipt["checker_view"],
        "ranking": _json_safe(result.receipt["ranking"]),
        "counts": _json_safe(result.receipt["counts"]),
        "outputs": {name: str(path) for name, path in outputs.items()},
        "split": (
            {
                "counts": dict(split_result.counts),
                "parent_method": split_result.parent_method,
                "requested_leakage_guard": getattr(
                    split_result,
                    "requested_leakage_guard",
                    split_result.leakage_guard,
                ),
                "leakage_guard": split_result.leakage_guard,
                "leakage_audit": dict(split_result.leakage_audit),
            }
            if split_result is not None
            else None
        ),
    }
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
