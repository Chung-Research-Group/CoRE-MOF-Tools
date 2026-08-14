"""Audited target-data joins for CoRE-MOF modelling workflows.

This module joins user-supplied response variables to validated release
metadata and optional public feature tables.  It intentionally performs no
fuzzy name matching, unit inference, condition inference, or value
imputation.  Earlier database identifiers are accepted only through an
explicit, hash-bound alias registry supplied by the user.

The implementation uses only the Python standard library.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
import hashlib
import io
import json
import math
import os
from pathlib import Path
import shutil
import stat
import tempfile
from types import MappingProxyType
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

from CoREMOF import __version__

from ._authority import state_fingerprint
from .dataset import (
    CoREMOFDataset,
    StructureRecord,
    _register_dataset_generation,
    _validate_dataset_generation_if_present,
)


TARGET_API_VERSION = "0.1.0"
_TARGET_MERGE_FACTORY_TOKEN = object()

CURRENT_FEATURE_TABLES = MappingProxyType(
    {
        "rac5": "features/rac5_features.csv",
        "zeo": "features/zeo_features.csv",
        "zeo_zero_probe": "features/zeo_zero_probe_features.csv",
        "topology": "features/topology_features.csv",
    }
)

_VALUE_TYPES = frozenset({"string", "float", "int", "bool", "json"})
_FORMATS = frozenset({"csv", "json", "jsonl"})


class TargetDataError(ValueError):
    """Raised when target inputs cannot be joined without ambiguity."""


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@dataclass(frozen=True)
class _FileSnapshot:
    """One stable regular-file byte generation and its receipt values."""

    path: Path
    data: bytes
    sha256: str
    size_bytes: int


def _stat_signature(value: os.stat_result) -> Tuple[int, int, int, int, int]:
    """Return identity and mutation-sensitive fields for one file stat."""

    return (
        int(value.st_dev),
        int(value.st_ino),
        int(value.st_size),
        int(getattr(value, "st_mtime_ns", int(value.st_mtime * 1_000_000_000))),
        int(getattr(value, "st_ctime_ns", int(value.st_ctime * 1_000_000_000))),
    )


def _read_snapshot_bytes(handle) -> bytes:
    """Read one input while keeping the open descriptor available for audit."""

    return handle.read()


def _capture_file(path: Union[str, Path]) -> _FileSnapshot:
    """Capture and hash one stable path generation exactly once.

    Parsing always consumes the returned bytes.  Path and descriptor stat
    signatures before and after the read reject concurrent replacement or
    in-place mutation instead of binding a receipt hash to different bytes.
    """

    path = Path(path).expanduser()
    path_before = path.stat()
    if not stat.S_ISREG(path_before.st_mode):
        raise TargetDataError("input is not a regular file: {}".format(path))
    with path.open("rb") as handle:
        descriptor_before = os.fstat(handle.fileno())
        if _stat_signature(path_before) != _stat_signature(descriptor_before):
            raise TargetDataError(
                "input changed or was replaced while being opened: {}".format(path)
            )
        data = _read_snapshot_bytes(handle)
        snapshot_sha256 = hashlib.sha256(data).hexdigest()
        handle.seek(0)
        verification_digest = hashlib.sha256()
        verification_size = 0
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            verification_digest.update(block)
            verification_size += len(block)
        descriptor_after = os.fstat(handle.fileno())
    path_after = path.stat()
    signatures = {
        _stat_signature(path_before),
        _stat_signature(descriptor_before),
        _stat_signature(descriptor_after),
        _stat_signature(path_after),
    }
    if (
        len(signatures) != 1
        or len(data) != descriptor_after.st_size
        or verification_size != len(data)
        or verification_digest.hexdigest() != snapshot_sha256
    ):
        raise TargetDataError(
            "input changed or was replaced during byte capture: {}".format(path)
        )
    return _FileSnapshot(
        path=path,
        data=data,
        sha256=snapshot_sha256,
        size_bytes=len(data),
    )


def _current_implementation_hashes() -> Dict[str, str]:
    package_root = Path(__file__).resolve().parent
    return {
        filename: _sha256_file(package_root / filename)
        for filename in ("_authority.py", "dataset.py", "labels.py", "targets.py")
    }


_IMPORTED_IMPLEMENTATION_HASHES = MappingProxyType(
    _current_implementation_hashes()
)


def _implementation_hashes() -> Mapping[str, str]:
    """Return import-bound sources, failing if their paths have since drifted."""

    current = _current_implementation_hashes()
    if current != dict(_IMPORTED_IMPLEMENTATION_HASHES):
        changed = sorted(
            filename
            for filename in set(current).union(_IMPORTED_IMPLEMENTATION_HASHES)
            if current.get(filename) != _IMPORTED_IMPLEMENTATION_HASHES.get(filename)
        )
        raise TargetDataError(
            "CoREMOF implementation source changed after module import: {}".format(
                ", ".join(changed)
            )
        )
    return _IMPORTED_IMPLEMENTATION_HASHES


def _canonical_data(value: object) -> object:
    """Thaw immutable mappings/tuples without accepting non-JSON containers."""

    if isinstance(value, Mapping):
        result = {}
        for key, item in value.items():
            if not isinstance(key, str) or not key or key != key.strip():
                raise TargetDataError(
                    "JSON contract mapping keys must be exact nonblank strings"
                )
            result[key] = _canonical_data(item)
        return result
    if isinstance(value, (list, tuple)):
        return [_canonical_data(item) for item in value]
    if isinstance(value, (set, frozenset)):
        raise TargetDataError("unordered sets are not valid target contract data")
    return value


def _canonical_json(value: object) -> str:
    try:
        return json.dumps(
            _canonical_data(value),
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
    except (TypeError, ValueError) as error:
        raise TargetDataError("value is not finite JSON data: {}".format(error))


def _deep_freeze(value: object) -> object:
    if isinstance(value, Mapping):
        result = {}
        for key, item in value.items():
            if not isinstance(key, str) or not key or key != key.strip():
                raise TargetDataError("mapping keys must be exact nonblank strings")
            result[key] = _deep_freeze(item)
        return MappingProxyType(result)
    if isinstance(value, (list, tuple)):
        return tuple(_deep_freeze(item) for item in value)
    if isinstance(value, (set, frozenset)):
        raise TargetDataError("unordered sets are not valid target contract data")
    return value


def _jsonable(value: object) -> object:
    if isinstance(value, Mapping):
        result = {}
        for key, item in value.items():
            if not isinstance(key, str) or not key or key != key.strip():
                raise TargetDataError("mapping keys must be exact nonblank strings")
            result[key] = _jsonable(item)
        return result
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, (set, frozenset)):
        raise TargetDataError("unordered sets are not valid target contract data")
    return value


def _clean_mapping(value: Optional[Mapping[str, object]], name: str) -> Dict[str, object]:
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise TypeError("{} must be a mapping".format(name))
    result = {}
    for key, item in value.items():
        if not isinstance(key, str) or not key or key != key.strip():
            raise TargetDataError("{} contains an empty target name".format(name))
        _canonical_json(item)
        result[key] = item
    return result


def _clean_columns(values: Optional[Iterable[str]], name: str) -> Optional[Tuple[str, ...]]:
    if values is None:
        return None
    if isinstance(values, str):
        values = (values,)
    elif not isinstance(values, (list, tuple)):
        raise TypeError("{} must be a string, list, or tuple".format(name))
    result = tuple(values)
    if not result or any(
        not isinstance(value, str)
        or not value
        or value != value.strip()
        for value in result
    ):
        raise TargetDataError("{} must contain non-empty column names".format(name))
    if len(set(result)) != len(result):
        raise TargetDataError("{} contains duplicate column names".format(name))
    return result


def _exact_name_tuple(values: object, name: str) -> Tuple[str, ...]:
    if isinstance(values, str) or not isinstance(values, (list, tuple)):
        raise TypeError("{} must be a list or tuple of strings".format(name))
    result = tuple(values)
    if any(
        not isinstance(value, str)
        or not value
        or value != value.strip()
        for value in result
    ):
        raise TargetDataError(
            "{} values must be exact nonblank strings".format(name)
        )
    if len(set(result)) != len(result):
        raise TargetDataError("{} must not contain duplicate names".format(name))
    return result


def _infer_format(path: Path, declared: Optional[str]) -> str:
    if declared is not None:
        if (
            not isinstance(declared, str)
            or not declared
            or declared != declared.strip()
        ):
            raise TargetDataError("target format must be an exact nonblank string")
        value = declared
    else:
        value = path.suffix.lower().lstrip(".")
        if value == "ndjson":
            value = "jsonl"
    if value not in _FORMATS:
        raise TargetDataError(
            "cannot determine target format for {}; use csv, json, or jsonl".format(
                path.name
            )
        )
    return value


@dataclass(frozen=True)
class TargetSource:
    """One target result file and its explicit scientific declarations.

    CSV values remain strings unless ``value_types`` explicitly requests a
    conversion.  Empty CSV cells are null by default.  JSON/JSONL native
    scalar types and JSON null are preserved.
    """

    path: Union[str, Path]
    name: Optional[str] = None
    id_column: str = "structure_id"
    target_columns: Optional[Sequence[str]] = None
    target_names: Mapping[str, str] = field(default_factory=dict)
    units: Mapping[str, object] = field(default_factory=dict)
    conditions: Mapping[str, object] = field(default_factory=dict)
    value_types: Mapping[str, str] = field(default_factory=dict)
    null_values: Sequence[str] = ("",)
    format: Optional[str] = None

    def __post_init__(self) -> None:
        path = Path(self.path).expanduser()
        name = self.name if self.name is not None else path.stem
        if not isinstance(name, str) or not name or name != name.strip():
            raise TargetDataError("target source name may not be empty")
        if (
            not isinstance(self.id_column, str)
            or not self.id_column
            or self.id_column != self.id_column.strip()
        ):
            raise TargetDataError("target id_column may not be empty")
        id_column = self.id_column
        columns = _clean_columns(self.target_columns, "target_columns")
        raw_target_names = _clean_mapping(self.target_names, "target_names")
        target_names = {}
        for input_column, raw_name in raw_target_names.items():
            if (
                not isinstance(raw_name, str)
                or not raw_name
                or raw_name != raw_name.strip()
            ):
                raise TargetDataError(
                    "target_names contains an empty canonical target name"
                )
            target_names[input_column] = raw_name
        if columns is not None:
            unknown_names = set(target_names).difference(columns)
            if unknown_names:
                raise TargetDataError(
                    "target_names refers to unselected input column(s): {}".format(
                        ", ".join(sorted(unknown_names))
                    )
                )
        units = _clean_mapping(self.units, "units")
        conditions = _clean_mapping(self.conditions, "conditions")
        raw_types = _clean_mapping(self.value_types, "value_types")
        value_types = {}
        for target, raw_type in raw_types.items():
            if not isinstance(raw_type, str) or raw_type != raw_type.strip():
                raise TargetDataError("value_types values must be exact strings")
            value_type = raw_type
            if value_type not in _VALUE_TYPES:
                raise TargetDataError(
                    "unknown value type {!r} for {}; choose one of {}".format(
                        raw_type, target, ", ".join(sorted(_VALUE_TYPES))
                    )
                )
            value_types[target] = value_type
        if isinstance(self.null_values, str) or not isinstance(
            self.null_values, (list, tuple)
        ):
            raise TypeError("null_values must be a list or tuple of strings")
        null_values = tuple(self.null_values)
        if any(not isinstance(value, str) for value in null_values):
            raise TargetDataError("null_values must contain only strings")
        object.__setattr__(self, "path", path)
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "id_column", id_column)
        object.__setattr__(self, "target_columns", columns)
        object.__setattr__(self, "target_names", MappingProxyType(target_names))
        object.__setattr__(self, "units", _deep_freeze(units))
        object.__setattr__(self, "conditions", _deep_freeze(conditions))
        object.__setattr__(self, "value_types", MappingProxyType(value_types))
        object.__setattr__(self, "null_values", null_values)
        object.__setattr__(self, "format", _infer_format(path, self.format))


@dataclass(frozen=True)
class AliasRegistry:
    """Explicit mapping from earlier identifiers to current public IDs."""

    path: Union[str, Path]
    current_id_column: str = "structure_id"
    alias_columns: Sequence[str] = ()
    format: Optional[str] = None
    null_values: Sequence[str] = ("",)

    def __post_init__(self) -> None:
        path = Path(self.path).expanduser()
        current = self.current_id_column
        aliases = _clean_columns(self.alias_columns, "alias_columns")
        if (
            not isinstance(current, str)
            or not current
            or current != current.strip()
        ):
            raise TargetDataError("alias current_id_column may not be empty")
        if not aliases:
            raise TargetDataError(
                "alias_columns must explicitly name every accepted earlier-ID column"
            )
        if current in aliases:
            raise TargetDataError("current_id_column cannot also be an alias column")
        object.__setattr__(self, "path", path)
        object.__setattr__(self, "current_id_column", current)
        object.__setattr__(self, "alias_columns", aliases)
        object.__setattr__(self, "format", _infer_format(path, self.format))
        if isinstance(self.null_values, str) or not isinstance(
            self.null_values, (list, tuple)
        ):
            raise TypeError("null_values must be a list or tuple of strings")
        object.__setattr__(self, "null_values", tuple(self.null_values))
        if any(not isinstance(value, str) for value in self.null_values):
            raise TargetDataError("null_values must contain only strings")


def _read_csv_records(
    snapshot: _FileSnapshot,
) -> Tuple[Tuple[str, ...], List[Tuple[int, Dict[str, object]]]]:
    path = snapshot.path
    with io.StringIO(snapshot.data.decode("utf-8-sig"), newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise TargetDataError("CSV has no header: {}".format(path))
        fields = tuple(reader.fieldnames)
        if len(set(fields)) != len(fields):
            raise TargetDataError("CSV has duplicate header names: {}".format(path))
        records = []
        for row_number, row in enumerate(reader, start=2):
            if None in row:
                raise TargetDataError(
                    "{}:{} contains more values than its header".format(path, row_number)
                )
            if any(value is None for value in row.values()):
                raise TargetDataError(
                    "{}:{} contains fewer values than its header".format(path, row_number)
                )
            records.append((row_number, dict(row)))
    return fields, records


def _read_json_records(
    snapshot: _FileSnapshot,
    id_column: str,
    target_columns: Optional[Sequence[str]],
) -> Tuple[Tuple[str, ...], List[Tuple[int, Dict[str, object]]]]:
    path = snapshot.path
    try:
        value = json.loads(snapshot.data.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise TargetDataError("cannot read JSON {}: {}".format(path, error))
    rows: List[Dict[str, object]] = []
    if isinstance(value, list):
        for index, item in enumerate(value, start=1):
            if not isinstance(item, Mapping):
                raise TargetDataError("{} JSON item {} is not an object".format(path, index))
            rows.append(dict(item))
    elif isinstance(value, Mapping):
        for raw_id, item in value.items():
            if isinstance(item, Mapping):
                row = dict(item)
            else:
                if target_columns is None or len(target_columns) != 1:
                    raise TargetDataError(
                        "JSON scalar maps require exactly one declared target column"
                    )
                row = {target_columns[0]: item}
            if not isinstance(raw_id, str) or not raw_id or raw_id != raw_id.strip():
                raise TargetDataError("JSON keyed identifiers must be exact strings")
            if id_column in row and (
                not isinstance(row[id_column], str)
                or row[id_column] != raw_id
            ):
                raise TargetDataError(
                    "{} keyed ID {!r} conflicts with embedded {!r}".format(
                        path, raw_id, row[id_column]
                    )
                )
            row[id_column] = raw_id
            rows.append(row)
    else:
        raise TargetDataError("{} must contain a JSON array or object".format(path))
    fields: List[str] = []
    for row in rows:
        for key in row:
            if not isinstance(key, str) or not key or key != key.strip():
                raise TargetDataError("JSON field names must be exact strings")
            if key not in fields:
                fields.append(key)
    return tuple(fields), [(index, row) for index, row in enumerate(rows, start=1)]


def _read_jsonl_records(
    snapshot: _FileSnapshot,
) -> Tuple[Tuple[str, ...], List[Tuple[int, Dict[str, object]]]]:
    path = snapshot.path
    fields: List[str] = []
    records = []
    try:
        text = snapshot.data.decode("utf-8")
    except UnicodeDecodeError as error:
        raise TargetDataError("cannot read JSONL {}: {}".format(path, error))
    with io.StringIO(text) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                item = json.loads(line)
            except json.JSONDecodeError as error:
                raise TargetDataError(
                    "cannot read JSONL {}:{}: {}".format(path, line_number, error)
                )
            if not isinstance(item, Mapping):
                raise TargetDataError(
                    "{}:{} JSONL record is not an object".format(path, line_number)
                )
            row = {str(key): value for key, value in item.items()}
            for key in row:
                if key not in fields:
                    fields.append(key)
            records.append((line_number, row))
    return tuple(fields), records


def _read_records(
    snapshot: _FileSnapshot,
    file_format: str,
    id_column: str,
    target_columns: Optional[Sequence[str]],
) -> Tuple[Tuple[str, ...], List[Tuple[int, Dict[str, object]]]]:
    if file_format == "csv":
        return _read_csv_records(snapshot)
    if file_format == "json":
        return _read_json_records(snapshot, id_column, target_columns)
    return _read_jsonl_records(snapshot)


def _convert_value(raw: object, source: TargetSource, target: str) -> object:
    if raw is None:
        return None
    if source.format == "csv" and isinstance(raw, str) and raw in source.null_values:
        return None
    value_type = source.value_types.get(target)
    if value_type is None:
        value = raw
    elif value_type == "string":
        value = str(raw)
    elif value_type == "float":
        if isinstance(raw, bool):
            raise TargetDataError("boolean cannot be converted to float")
        try:
            value = float(raw)
        except (TypeError, ValueError):
            raise TargetDataError("cannot convert {!r} to float".format(raw))
    elif value_type == "int":
        if isinstance(raw, bool):
            raise TargetDataError("boolean cannot be converted to int")
        if isinstance(raw, int):
            value = raw
        elif isinstance(raw, str) and raw.strip() and raw.strip().lstrip("+-").isdigit():
            value = int(raw.strip())
        else:
            raise TargetDataError("cannot convert {!r} to int exactly".format(raw))
    elif value_type == "bool":
        if isinstance(raw, bool):
            value = raw
        elif isinstance(raw, str) and raw.strip().casefold() in ("true", "false"):
            value = raw.strip().casefold() == "true"
        else:
            raise TargetDataError("cannot convert {!r} to bool".format(raw))
    else:  # explicit JSON conversion
        if isinstance(raw, str):
            try:
                value = json.loads(raw)
            except json.JSONDecodeError as error:
                raise TargetDataError("cannot parse JSON value {!r}: {}".format(raw, error))
        else:
            value = raw
    if isinstance(value, float) and not math.isfinite(value):
        raise TargetDataError("target values must be finite; got {!r}".format(value))
    if (
        value_type != "json"
        and value is not None
        and not isinstance(value, (str, int, float, bool))
    ):
        raise TargetDataError(
            "target {!r} is not a scalar JSON value".format(target)
        )
    _canonical_json(value)
    return value


def _read_alias_registry(
    registry: AliasRegistry, current_ids: Sequence[str]
) -> Tuple[Dict[str, str], Mapping[str, object]]:
    snapshot = _capture_file(registry.path)
    fields, rows = _read_records(
        snapshot,
        str(registry.format),
        registry.current_id_column,
        registry.alias_columns,
    )
    required = (registry.current_id_column,) + tuple(registry.alias_columns)
    missing = [column for column in required if column not in fields]
    if missing:
        raise TargetDataError(
            "alias registry is missing column(s): {}".format(", ".join(missing))
        )
    current_set = set(current_ids)
    lookup = {structure_id: structure_id for structure_id in current_ids}
    alias_count = 0
    for row_number, row in rows:
        current = row[registry.current_id_column]
        if not isinstance(current, str) or not current or current != current.strip():
            raise TargetDataError(
                "alias registry identifiers must be exact nonblank strings"
            )
        if current not in current_set:
            raise TargetDataError(
                "alias registry row {} refers to unknown current ID {!r}".format(
                    row_number, current
                )
            )
        for column in registry.alias_columns:
            raw = row[column]
            if raw is None or (isinstance(raw, str) and raw in registry.null_values):
                continue
            if not isinstance(raw, str) or not raw or raw != raw.strip():
                raise TargetDataError(
                    "alias registry aliases must be exact nonblank strings"
                )
            alias = raw
            existing = lookup.get(alias)
            if existing is not None and existing != current:
                raise TargetDataError(
                    "ambiguous alias {!r} maps to both {} and {}".format(
                        alias, existing, current
                    )
                )
            if existing is None:
                alias_count += 1
            lookup[alias] = current
    receipt = {
        "file_name": snapshot.path.name,
        "format": registry.format,
        "sha256": snapshot.sha256,
        "size_bytes": snapshot.size_bytes,
        "current_id_column": registry.current_id_column,
        "alias_columns": list(registry.alias_columns),
        "row_count": len(rows),
        "unique_alias_count": alias_count,
    }
    return lookup, receipt


def _read_feature_table(
    dataset: CoREMOFDataset, name: str
) -> Tuple[Tuple[str, ...], Dict[str, Dict[str, str]], Mapping[str, object]]:
    if name not in CURRENT_FEATURE_TABLES:
        raise TargetDataError(
            "unknown feature table {!r}; choose one of {}".format(
                name, ", ".join(sorted(CURRENT_FEATURE_TABLES))
            )
        )
    relative = CURRENT_FEATURE_TABLES[name]
    path = dataset.release_root / relative
    snapshot = _capture_file(path)
    fields, rows = _read_csv_records(snapshot)
    if "structure_id" not in fields:
        raise TargetDataError("{} has no structure_id column".format(relative))
    feature_columns = tuple(field for field in fields if field != "structure_id")
    indexed: Dict[str, Dict[str, str]] = {}
    for row_number, row in rows:
        structure_id = row["structure_id"]
        if (
            not isinstance(structure_id, str)
            or not structure_id
            or structure_id != structure_id.strip()
        ):
            raise TargetDataError("{}:{} has an empty structure_id".format(path, row_number))
        if structure_id in indexed:
            raise TargetDataError("{} has duplicate structure_id {}".format(path, structure_id))
        indexed[structure_id] = {column: row[column] for column in feature_columns}
    expected = set(dataset.structure_ids)
    observed = set(indexed)
    if expected != observed:
        raise TargetDataError(
            "{} ID set differs from release metadata: missing={} extra={}".format(
                relative,
                sorted(expected.difference(observed))[:5],
                sorted(observed.difference(expected))[:5],
            )
        )
    receipt = {
        "name": name,
        "release_path": relative,
        "sha256": snapshot.sha256,
        "size_bytes": snapshot.size_bytes,
        "column_count": len(feature_columns),
        "row_count": len(rows),
    }
    return feature_columns, indexed, receipt


def _stable_target_digest(
    structure_ids: Sequence[str],
    target_columns: Sequence[str],
    values_by_id: Mapping[str, Mapping[str, object]],
) -> str:
    digest = hashlib.sha256()
    for structure_id in sorted(structure_ids):
        payload = {
            target: values_by_id.get(structure_id, {}).get(target)
            for target in sorted(target_columns)
        }
        digest.update(structure_id.encode("utf-8"))
        digest.update(b"\0")
        digest.update(_canonical_json(payload).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


class TargetMergedDataset(CoREMOFDataset):
    """Validated release rows joined to feature columns and user targets."""

    def __init__(
        self,
        base_dataset: CoREMOFDataset,
        records: Sequence[StructureRecord],
        target_columns: Sequence[str],
        feature_columns: Sequence[str],
        target_values_by_id: Mapping[str, Mapping[str, object]],
        target_provenance_by_id: Mapping[str, Mapping[str, Sequence[Mapping[str, object]]]],
        target_definitions: Mapping[str, Mapping[str, object]],
        target_receipt: Mapping[str, object],
        input_hashes: Mapping[str, str],
        *,
        _factory_token: object = None,
    ) -> None:
        if _factory_token is not _TARGET_MERGE_FACTORY_TOKEN:
            raise TargetDataError(
                "TargetMergedDataset must be created by merge_targets()"
            )
        base_is_official = _validate_dataset_generation_if_present(base_dataset)
        if not isinstance(base_dataset, CoREMOFDataset):
            raise TypeError("base_dataset must be a CoREMOFDataset")
        target_columns_tuple = _exact_name_tuple(target_columns, "target_columns")
        feature_columns_tuple = _exact_name_tuple(feature_columns, "feature_columns")
        if set(target_columns_tuple).intersection(feature_columns_tuple):
            raise TargetDataError("target and feature columns must be disjoint")
        current_ids = tuple(base_dataset.structure_ids)
        if len(records) != len(current_ids) or tuple(
            record.structure_id for record in records
        ) != current_ids:
            raise TargetDataError(
                "target-merged records must preserve the base release order and universe"
            )
        if set(target_values_by_id) != set(current_ids):
            raise TargetDataError(
                "target values must contain exactly the base release structure IDs"
            )

        frozen_records = []
        for base_record, record in zip(base_dataset.records, records):
            if not isinstance(record, StructureRecord) or not isinstance(
                record.metadata, Mapping
            ):
                raise TargetDataError(
                    "target-merged rows must be StructureRecord values with mappings"
                )
            if any(
                key not in record.metadata
                or type(record.metadata[key]) is not type(value)
                or record.metadata[key] != value
                for key, value in base_record.metadata.items()
            ):
                raise TargetDataError(
                    "target merge changed release metadata for {}".format(
                        base_record.structure_id
                    )
                )
            expected_columns = set(base_record.metadata).union(
                feature_columns_tuple, target_columns_tuple
            )
            if set(record.metadata) != expected_columns:
                raise TargetDataError(
                    "target-merged metadata has missing or unexpected columns"
                )
            values = target_values_by_id[record.structure_id]
            if not isinstance(values, Mapping) or set(values) != set(target_columns_tuple):
                raise TargetDataError(
                    "each target-value row must contain exactly the target columns"
                )
            for target in target_columns_tuple:
                if _canonical_json(record.metadata[target]) != _canonical_json(
                    values[target]
                ):
                    raise TargetDataError(
                        "joined target value differs from its target table"
                    )
            frozen_records.append(
                StructureRecord(
                    structure_id=record.structure_id,
                    metadata=_deep_freeze(record.metadata),
                    parent_groups=_deep_freeze(base_record.parent_groups),
                    cif_manifest=(
                        _deep_freeze(base_record.cif_manifest)
                        if base_record.cif_manifest is not None
                        else None
                    ),
                )
            )
        raw_parent_by_id = base_dataset.parent_by_id
        if not isinstance(raw_parent_by_id, Mapping):
            raise TypeError("base parent_by_id must be a mapping")
        normalized_parent_by_id = {}
        for structure_id in current_ids:
            row = raw_parent_by_id.get(structure_id, {})
            if row is None:
                row = {}
            if not isinstance(row, Mapping):
                raise TypeError("base parent_by_id rows must be mappings")
            normalized_parent_by_id[structure_id] = _deep_freeze(row)
        super().__init__(
            release_root=base_dataset.release_root,
            records=frozen_records,
            dataset_info=_deep_freeze(base_dataset.dataset_info),
            parent_group_methods=_deep_freeze(base_dataset.parent_group_methods),
            parent_by_id=MappingProxyType(normalized_parent_by_id),
            input_hashes=_deep_freeze(input_hashes),
            cif_files_verified=base_dataset.cif_files_verified,
        )
        self.base_dataset = base_dataset
        self.target_columns = target_columns_tuple
        self.feature_columns = feature_columns_tuple
        self.target_values_by_id = _deep_freeze(
            {
                structure_id: dict(values)
                for structure_id, values in target_values_by_id.items()
            }
        )
        self.target_provenance_by_id = _deep_freeze(target_provenance_by_id)
        self.target_definitions = _deep_freeze(target_definitions)
        self.target_input_receipt = _deep_freeze(target_receipt)
        self._validate_target_generation_contract()
        self._authority_extra_state = (
            "coremof-target-transformation/1",
            state_fingerprint(self.target_input_receipt),
            state_fingerprint(self.target_values_by_id),
            state_fingerprint(self.target_provenance_by_id),
            state_fingerprint(self.target_definitions),
            self.target_columns,
            self.feature_columns,
        )
        self._authority_extra_bindings = MappingProxyType(
            {
                "base_dataset": self.base_dataset,
                "target_columns": self.target_columns,
                "feature_columns": self.feature_columns,
                "target_values_by_id": self.target_values_by_id,
                "target_provenance_by_id": self.target_provenance_by_id,
                "target_definitions": self.target_definitions,
                "target_input_receipt": self.target_input_receipt,
            }
        )
        _register_dataset_generation(
            self,
            kind="target_merged",
            official_release_source=base_is_official,
            base_dataset=base_dataset,
        )

    def _validate_target_generation_contract(self) -> None:
        receipt = self.target_input_receipt
        expected_keys = {
            "schema_version",
            "target_api_version",
            "implementation",
            "dataset_version",
            "release_structure_count",
            "target_columns",
            "target_definitions",
            "sources",
            "alias_registry",
            "feature_tables",
            "counts",
            "target_values_sha256",
            "policies",
        }
        if "config" in receipt:
            expected_keys.add("config")
        if not isinstance(receipt, Mapping) or set(receipt) != expected_keys:
            raise TargetDataError("target merge receipt has an invalid closed contract")
        if (
            receipt.get("schema_version") != "coremof-target-merge-receipt/1.0"
            or receipt.get("target_api_version") != TARGET_API_VERSION
            or receipt.get("dataset_version") != self.dataset_version
            or receipt.get("release_structure_count") != len(self)
            or tuple(receipt.get("target_columns", ())) != self.target_columns
            or _canonical_json(receipt.get("target_definitions"))
            != _canonical_json(self.target_definitions)
            or receipt.get("target_values_sha256")
            != _stable_target_digest(
                self.structure_ids,
                self.target_columns,
                self.target_values_by_id,
            )
        ):
            raise TargetDataError(
                "target merge receipt differs from the exact joined generation"
            )
        implementation = receipt.get("implementation")
        if not isinstance(implementation, Mapping) or set(implementation) != {
            "package",
            "package_version",
            "target_api_version",
            "source_sha256",
        }:
            raise TargetDataError("target merge receipt implementation is invalid")
        if (
            implementation.get("package") != "CoREMOF-tools"
            or implementation.get("package_version") != __version__
            or implementation.get("target_api_version") != TARGET_API_VERSION
            or dict(implementation.get("source_sha256", {}))
            != dict(_implementation_hashes())
        ):
            raise TargetDataError(
                "target merge receipt implementation does not match executing sources"
            )
        if not isinstance(self.target_provenance_by_id, Mapping):
            raise TargetDataError("target provenance must be a mapping")
        unknown_ids = set(self.target_provenance_by_id).difference(self.structure_ids)
        if unknown_ids:
            raise TargetDataError("target provenance contains unknown structure IDs")
        for structure_id, by_target in self.target_provenance_by_id.items():
            if not isinstance(by_target, Mapping) or set(by_target).difference(
                self.target_columns
            ):
                raise TargetDataError("target provenance contains unknown target columns")
            for target, observations in by_target.items():
                if not isinstance(observations, (list, tuple)):
                    raise TargetDataError("target provenance observations must be a sequence")
                if any(not isinstance(item, Mapping) for item in observations):
                    raise TargetDataError("target provenance observations must be mappings")

    def target_values(self, structure_id: str) -> Mapping[str, object]:
        """Return every target column for one current public structure ID."""

        if structure_id not in self._by_id:
            raise KeyError(structure_id)
        values = self.target_values_by_id.get(structure_id, {})
        return MappingProxyType(
            {target: values.get(target) for target in self.target_columns}
        )

    def receipt(self) -> Mapping[str, object]:
        """Return the hash-bound target/feature merge receipt."""

        _validate_dataset_generation_if_present(self)
        return _jsonable(self.target_input_receipt)  # type: ignore[return-value]

    @staticmethod
    def _csv_cell(value: object) -> object:
        if value is None:
            return ""
        if isinstance(value, bool):
            return "true" if value else "false"
        if isinstance(value, (Mapping, list, tuple)):
            return _canonical_json(_jsonable(value))
        return value

    @staticmethod
    def _quarantine_published_generation(
        target: Path, identity: Path, quarantine: Path
    ) -> bool:
        """Remove only the published generation identified by ``identity``.

        The final pathname is first renamed into the private staging directory.
        This makes the subsequent identity check and deletion operate on a name
        that a same-stem writer cannot replace between ``samefile`` and
        ``unlink``.  If the detached generation is foreign, restore it with an
        atomic create-if-absent link rather than overwriting a newer final.

        Return ``False`` only when a foreign generation must remain quarantined
        for safety; the caller then preserves the staging directory.
        """

        try:
            os.replace(str(target), str(quarantine))
        except FileNotFoundError:
            return True

        try:
            published_by_us = identity.exists() and os.path.samefile(
                str(identity), str(quarantine)
            )
        except OSError:
            published_by_us = False
        if published_by_us:
            quarantine.unlink()
            return True

        try:
            os.link(str(quarantine), str(target))
        except OSError:
            # Do not destroy a foreign inode merely because it could not be
            # restored without overwriting a newer final.  Leaving the private
            # staging directory is preferable to losing either generation.
            return False
        quarantine.unlink()
        return True

    def write(
        self,
        output_directory: Union[str, Path],
        stem: str = "coremof_targets",
        overwrite: bool = False,
    ) -> Tuple[Path, Path, Path]:
        """Write merged CSV, provenance JSONL, and receipt JSON as one bundle.

        With ``overwrite=False``, create-if-absent publication and rollback are
        safe against another writer replacing a file during this call.  Since
        three ordinary files cannot be replaced as one filesystem transaction,
        ``overwrite=True`` has explicit single-writer semantics: callers must
        serialize writers for the same output directory and stem (for example
        with an external lock).
        """

        if type(overwrite) is not bool:
            raise TypeError("overwrite must be a boolean")
        _validate_dataset_generation_if_present(self)
        if (
            not isinstance(stem, str)
            or not stem
            or stem in {".", ".."}
            or Path(stem).name != stem
        ):
            raise ValueError("stem must be a non-empty filename stem without directories")
        directory = Path(output_directory)
        directory.mkdir(parents=True, exist_ok=True)
        targets = (
            directory / (stem + ".csv"),
            directory / (stem + ".provenance.jsonl"),
            directory / (stem + ".json"),
        )
        if not overwrite and any(path.exists() for path in targets):
            raise FileExistsError(
                "Refusing to overwrite existing target output: {}".format(
                    ", ".join(str(path) for path in targets if path.exists())
                )
        )
        staging = Path(tempfile.mkdtemp(prefix=".{}-".format(stem), dir=str(directory)))
        staged = tuple(staging / path.name for path in targets)
        backups: Dict[Path, Path] = {}
        published_identities: Dict[Path, Path] = {}
        published: List[Path] = []
        remove_staging = True
        try:
            fieldnames = tuple(self.records[0].metadata) if self.records else ("structure_id",)
            with staged[0].open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames)
                writer.writeheader()
                for record in self.records:
                    writer.writerow(
                        {
                            field: self._csv_cell(record.metadata.get(field))
                            for field in fieldnames
                        }
                    )
                handle.flush()
                os.fsync(handle.fileno())
            with staged[1].open("w", encoding="utf-8") as handle:
                for structure_id in sorted(self.target_provenance_by_id):
                    by_target = self.target_provenance_by_id[structure_id]
                    for target in sorted(by_target):
                        for observation in by_target[target]:
                            payload = {
                                "structure_id": structure_id,
                                "target": target,
                                **_jsonable(observation),
                            }
                            handle.write(_canonical_json(payload) + "\n")
                handle.flush()
                os.fsync(handle.fileno())
            with staged[2].open("w", encoding="utf-8") as handle:
                json.dump(self.receipt(), handle, indent=2, sort_keys=True, allow_nan=False)
                handle.write("\n")
                handle.flush()
                os.fsync(handle.fileno())
            if overwrite:
                for target in targets:
                    if target.exists():
                        backup = staging / (target.name + ".previous")
                        try:
                            os.link(str(target), str(backup))
                        except OSError:
                            shutil.copy2(str(target), str(backup))
                        backups[target] = backup
            for source, target in zip(staged, targets):
                if overwrite:
                    if target not in backups:
                        identity = staging / (target.name + ".published")
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
                    quarantine = staging / (target.name + ".rollback")
                    if not self._quarantine_published_generation(
                        target, identity, quarantine
                    ):
                        remove_staging = False
            raise
        finally:
            if remove_staging:
                shutil.rmtree(staging, ignore_errors=True)
        return targets


def merge_targets(
    dataset: Union[str, Path, CoREMOFDataset],
    sources: Sequence[Union[str, Path, TargetSource]],
    alias_registry: Optional[AliasRegistry] = None,
    feature_tables: Sequence[str] = (),
    verify_cif_files: bool = False,
    *,
    _config_receipt: Optional[Mapping[str, object]] = None,
    _config_factory_token: object = None,
) -> TargetMergedDataset:
    """Join one or more target files to a complete release universe.

    Exact current public IDs are always accepted.  Any other identifier is an
    error unless it appears unambiguously in ``alias_registry``.  Target names
    may not collide with release metadata or selected feature columns.
    """

    if type(verify_cif_files) is not bool:
        raise TypeError("verify_cif_files must be a boolean")
    if isinstance(dataset, (str, Path)):
        base = CoREMOFDataset.from_release(dataset, verify_cif_files=verify_cif_files)
    elif isinstance(dataset, CoREMOFDataset):
        base = dataset
        if verify_cif_files and not base.cif_files_verified:
            raise ValueError(
                "the supplied dataset was not loaded with verify_cif_files=True"
            )
    else:
        raise TypeError("dataset must be a release path or CoREMOFDataset")
    if isinstance(base, TargetMergedDataset):
        raise TargetDataError("target data are already attached to this dataset")
    if isinstance(sources, (str, Path, TargetSource)):
        sources = (sources,)
    elif not isinstance(sources, (list, tuple)):
        raise TypeError("sources must be a target source or a list/tuple of sources")
    normalized_sources = tuple(
        source if isinstance(source, TargetSource) else TargetSource(source)
        for source in sources
    )
    if not normalized_sources:
        raise TargetDataError("at least one target source is required")
    source_names = [source.name for source in normalized_sources]
    if len(set(source_names)) != len(source_names):
        raise TargetDataError("target source names must be unique")

    current_ids = base.structure_ids
    lookup = {structure_id: structure_id for structure_id in current_ids}
    alias_receipt = None
    if alias_registry is not None:
        lookup, alias_receipt = _read_alias_registry(alias_registry, current_ids)

    if isinstance(feature_tables, str):
        feature_tables = (feature_tables,)
    elif not isinstance(feature_tables, (list, tuple)):
        raise TypeError("feature_tables must be a string, list, or tuple")
    feature_names = tuple(feature_tables)
    if any(
        not isinstance(name, str)
        or not name
        or name != name.strip()
        for name in feature_names
    ):
        raise TargetDataError("feature_tables must contain exact nonblank strings")
    if len(set(feature_names)) != len(feature_names):
        raise TargetDataError("feature_tables contains duplicates")
    feature_columns: List[str] = []
    feature_rows: Dict[str, Dict[str, str]] = {
        structure_id: {} for structure_id in current_ids
    }
    feature_receipts = []
    base_columns = set(base.metadata_rows[0]) if base.metadata_rows else {"structure_id"}
    for name in feature_names:
        columns, rows, receipt = _read_feature_table(base, name)
        collisions = base_columns.intersection(columns).union(
            set(feature_columns).intersection(columns)
        )
        if collisions:
            raise TargetDataError(
                "feature table {} collides with existing column(s): {}".format(
                    name, ", ".join(sorted(collisions))
                )
            )
        feature_columns.extend(columns)
        feature_receipts.append(receipt)
        for structure_id in current_ids:
            feature_rows[structure_id].update(rows[structure_id])

    definitions: Dict[str, Mapping[str, object]] = {}
    target_order: List[str] = []
    values_by_id: Dict[str, Dict[str, object]] = {}
    provenance_by_id: Dict[str, Dict[str, List[Mapping[str, object]]]] = {}
    source_receipts = []
    exact_count = 0
    alias_count = 0
    observation_count = 0
    null_observation_count = 0
    for source in normalized_sources:
        path = Path(source.path).resolve()
        snapshot = _capture_file(path)
        source_sha = snapshot.sha256
        fields, rows = _read_records(
            snapshot, str(source.format), source.id_column, source.target_columns
        )
        if source.id_column not in fields:
            raise TargetDataError(
                "{} is missing ID column {!r}".format(path.name, source.id_column)
            )
        input_targets = source.target_columns
        if input_targets is None:
            input_targets = tuple(field for field in fields if field != source.id_column)
        if not input_targets:
            raise TargetDataError("{} contains no target columns".format(path.name))
        if source.id_column in input_targets:
            raise TargetDataError(
                "{} ID column cannot also be a target column".format(path.name)
            )
        missing = [target for target in input_targets if target not in fields]
        if missing:
            raise TargetDataError(
                "{} is missing target column(s): {}".format(path.name, ", ".join(missing))
            )
        unknown_names = set(source.target_names).difference(input_targets)
        if unknown_names:
            raise TargetDataError(
                "{} target_names refers to unselected input column(s): {}".format(
                    path.name, ", ".join(sorted(unknown_names))
                )
            )
        canonical_targets = tuple(
            source.target_names.get(input_target, input_target)
            for input_target in input_targets
        )
        if len(set(canonical_targets)) != len(canonical_targets):
            raise TargetDataError(
                "{} maps multiple input columns to the same canonical target".format(
                    path.name
                )
            )
        undeclared = set(source.units).union(source.conditions).union(source.value_types).difference(canonical_targets)
        if undeclared:
            raise TargetDataError(
                "{} declarations refer to unselected target(s): {}".format(
                    path.name, ", ".join(sorted(undeclared))
                )
            )
        for target in canonical_targets:
            if target == "structure_id" or target in base_columns or target in feature_columns:
                raise TargetDataError(
                    "target column {!r} collides with release/feature metadata".format(target)
                )
            definition = {
                "unit": source.units.get(target),
                "conditions": source.conditions.get(target),
                "value_type": source.value_types.get(target),
            }
            if target in definitions and _canonical_json(definitions[target]) != _canonical_json(definition):
                raise TargetDataError(
                    "target {!r} has conflicting unit/condition/type declarations".format(target)
                )
            if target not in definitions:
                definitions[target] = definition
                target_order.append(target)
        source_exact = 0
        source_alias = 0
        source_null = 0
        for row_number, row in rows:
            raw_id = row.get(source.id_column)
            if raw_id is None:
                raise TargetDataError(
                    "{}:{} has a null structure identifier".format(path.name, row_number)
                )
            if (
                not isinstance(raw_id, str)
                or not raw_id
                or raw_id != raw_id.strip()
            ):
                raise TargetDataError(
                    "{}:{} has an empty structure identifier".format(path.name, row_number)
                )
            input_id = raw_id
            structure_id = lookup.get(input_id)
            if structure_id is None:
                raise TargetDataError(
                    "{}:{} identifier {!r} is neither a current public ID nor an audited alias".format(
                        path.name, row_number, input_id
                    )
                )
            resolution = "CURRENT" if input_id == structure_id else "ALIAS"
            if resolution == "CURRENT":
                exact_count += 1
                source_exact += 1
            else:
                alias_count += 1
                source_alias += 1
            for input_target, target in zip(input_targets, canonical_targets):
                try:
                    value = _convert_value(row.get(input_target), source, target)
                except TargetDataError as error:
                    raise TargetDataError(
                        "{}:{} target {}: {}".format(path.name, row_number, target, error)
                    )
                observation_count += 1
                if value is None:
                    null_observation_count += 1
                    source_null += 1
                current_values = values_by_id.setdefault(structure_id, {})
                if target in current_values:
                    old = current_values[target]
                    if old is not None and value is not None and _canonical_json(old) != _canonical_json(value):
                        raise TargetDataError(
                            "conflicting duplicate target {} for {}: {} versus {}".format(
                                target,
                                structure_id,
                                _canonical_json(old),
                                _canonical_json(value),
                            )
                        )
                    if old is None and value is not None:
                        current_values[target] = value
                else:
                    current_values[target] = value
                provenance_by_id.setdefault(structure_id, {}).setdefault(target, []).append(
                    {
                        "source": source.name,
                        "source_file": path.name,
                        "source_sha256": source_sha,
                        "row_number": row_number,
                        "input_column": input_target,
                        "input_structure_id": input_id,
                        "id_resolution": resolution,
                        "value": value,
                    }
                )
        source_receipts.append(
            {
                "name": source.name,
                "file_name": path.name,
                "format": source.format,
                "sha256": source_sha,
                "size_bytes": snapshot.size_bytes,
                "id_column": source.id_column,
                "input_target_columns": list(input_targets),
                "target_names": {
                    input_target: target
                    for input_target, target in zip(input_targets, canonical_targets)
                },
                "target_columns": list(canonical_targets),
                "units": _jsonable(source.units),
                "conditions": _jsonable(source.conditions),
                "value_types": dict(source.value_types),
                "row_count": len(rows),
                "current_id_row_count": source_exact,
                "alias_id_row_count": source_alias,
                "null_observation_count": source_null,
            }
        )

    target_columns = tuple(target_order)
    complete_values = {
        structure_id: {
            target: values_by_id.get(structure_id, {}).get(target)
            for target in target_columns
        }
        for structure_id in current_ids
    }
    records = []
    for record in base.records:
        merged = dict(record.metadata)
        merged.update(feature_rows[record.structure_id])
        merged.update(
            {
                target: _deep_freeze(value)
                for target, value in complete_values[record.structure_id].items()
            }
        )
        records.append(
            StructureRecord(
                structure_id=record.structure_id,
                metadata=MappingProxyType(merged),
                parent_groups=record.parent_groups,
                cif_manifest=record.cif_manifest,
            )
        )
    target_digest = _stable_target_digest(current_ids, target_columns, complete_values)
    with_any = sum(
        any(value is not None for value in complete_values[structure_id].values())
        for structure_id in current_ids
    )
    with_all = sum(
        all(complete_values[structure_id][target] is not None for target in target_columns)
        for structure_id in current_ids
    )
    receipt = {
        "schema_version": "coremof-target-merge-receipt/1.0",
        "target_api_version": TARGET_API_VERSION,
        "implementation": {
            "package": "CoREMOF-tools",
            "package_version": __version__,
            "target_api_version": TARGET_API_VERSION,
            "source_sha256": dict(_implementation_hashes()),
        },
        "dataset_version": base.dataset_version,
        "release_structure_count": len(base),
        "target_columns": list(target_columns),
        "target_definitions": _jsonable(definitions),
        "sources": source_receipts,
        "alias_registry": alias_receipt,
        "feature_tables": feature_receipts,
        "counts": {
            "source_files": len(normalized_sources),
            "observations": observation_count,
            "null_observations": null_observation_count,
            "current_id_rows": exact_count,
            "alias_id_rows": alias_count,
            "structures_with_any_target": with_any,
            "structures_with_all_targets": with_all,
        },
        "target_values_sha256": target_digest,
        "policies": {
            "identifier_matching": "EXACT_CURRENT_OR_EXPLICIT_ALIAS",
            "fuzzy_matching": False,
            "units_inferred": False,
            "conditions_inferred": False,
            "missing_targets_imputed": False,
            "conflicting_duplicates": "ERROR",
        },
    }
    if _config_receipt is not None:
        if _config_factory_token is not _TARGET_MERGE_FACTORY_TOKEN:
            raise TargetDataError(
                "target config receipt can be attached only by merge_targets_from_config()"
            )
        receipt["config"] = dict(_config_receipt)
    input_hashes = dict(base.input_hashes)
    for item in feature_receipts:
        input_hashes[str(item["release_path"])] = str(item["sha256"])
    for item in source_receipts:
        input_hashes["targets/{}".format(item["name"])] = str(item["sha256"])
    if alias_receipt is not None:
        input_hashes["targets/alias_registry"] = str(alias_receipt["sha256"])
    if _config_receipt is not None:
        input_hashes["targets/config"] = _config_receipt["sha256"]
    return TargetMergedDataset(
        base_dataset=base,
        records=records,
        target_columns=target_columns,
        feature_columns=feature_columns,
        target_values_by_id=complete_values,
        target_provenance_by_id=provenance_by_id,
        target_definitions=definitions,
        target_receipt=receipt,
        input_hashes=input_hashes,
        _factory_token=_TARGET_MERGE_FACTORY_TOKEN,
    )


def load_target_config(path: Union[str, Path]) -> Mapping[str, object]:
    """Read a JSON target configuration and resolve paths beside that file."""

    config_path = Path(path).expanduser().resolve()
    snapshot = _capture_file(config_path)
    try:
        config = json.loads(snapshot.data.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise TargetDataError("cannot read target config {}: {}".format(config_path, error))
    if not isinstance(config, Mapping):
        raise TargetDataError("target config must contain a JSON object")
    allowed = {"sources", "alias_registry", "feature_tables"}
    unknown = set(config).difference(allowed)
    if unknown:
        raise TargetDataError(
            "unknown target config key(s): {}".format(", ".join(sorted(unknown)))
        )
    raw_sources = config.get("sources")
    if not isinstance(raw_sources, list) or not raw_sources:
        raise TargetDataError("target config sources must be a non-empty array")
    sources = []
    source_keys = {
        "path", "name", "id_column", "target_columns", "target_names", "units", "conditions",
        "value_types", "null_values", "format",
    }
    for index, raw in enumerate(raw_sources, start=1):
        if not isinstance(raw, Mapping):
            raise TargetDataError("target config source {} must be an object".format(index))
        unknown_source = set(raw).difference(source_keys)
        if unknown_source:
            raise TargetDataError(
                "target config source {} has unknown key(s): {}".format(
                    index, ", ".join(sorted(unknown_source))
                )
            )
        if "path" not in raw:
            raise TargetDataError("target config source {} has no path".format(index))
        values = dict(raw)
        raw_source_path = values["path"]
        if (
            not isinstance(raw_source_path, str)
            or not raw_source_path
            or raw_source_path != raw_source_path.strip()
        ):
            raise TargetDataError(
                "target config source {} path must be an exact nonblank string".format(
                    index
                )
            )
        source_path = Path(raw_source_path).expanduser()
        if not source_path.is_absolute():
            source_path = config_path.parent / source_path
        values["path"] = source_path
        sources.append(TargetSource(**values))
    raw_alias = config.get("alias_registry")
    alias = None
    if raw_alias is not None:
        if not isinstance(raw_alias, Mapping):
            raise TargetDataError("alias_registry must be an object")
        alias_keys = {"path", "current_id_column", "alias_columns", "format", "null_values"}
        unknown_alias = set(raw_alias).difference(alias_keys)
        if unknown_alias:
            raise TargetDataError(
                "alias_registry has unknown key(s): {}".format(", ".join(sorted(unknown_alias)))
            )
        if "path" not in raw_alias:
            raise TargetDataError("alias_registry has no path")
        values = dict(raw_alias)
        raw_alias_path = values["path"]
        if (
            not isinstance(raw_alias_path, str)
            or not raw_alias_path
            or raw_alias_path != raw_alias_path.strip()
        ):
            raise TargetDataError(
                "alias_registry path must be an exact nonblank string"
            )
        alias_path = Path(raw_alias_path).expanduser()
        if not alias_path.is_absolute():
            alias_path = config_path.parent / alias_path
        values["path"] = alias_path
        alias = AliasRegistry(**values)
    raw_features = config.get("feature_tables", ())
    if isinstance(raw_features, str) or not isinstance(raw_features, (list, tuple)):
        raise TargetDataError("feature_tables must be an array")
    feature_tables = tuple(raw_features)
    if any(
        not isinstance(value, str)
        or not value
        or value != value.strip()
        for value in feature_tables
    ):
        raise TargetDataError(
            "feature_tables values must be exact nonblank strings"
        )
    if len(set(feature_tables)) != len(feature_tables):
        raise TargetDataError("feature_tables must not contain duplicates")
    return MappingProxyType(
        {
            "sources": tuple(sources),
            "alias_registry": alias,
            "feature_tables": feature_tables,
            "config_file_name": config_path.name,
            "config_sha256": snapshot.sha256,
            "config_size_bytes": snapshot.size_bytes,
        }
    )


def merge_targets_from_config(
    dataset: Union[str, Path, CoREMOFDataset],
    config_path: Union[str, Path],
    verify_cif_files: bool = False,
) -> TargetMergedDataset:
    """Load a JSON configuration and return its merged modelling dataset."""

    config = load_target_config(config_path)
    config_receipt = {
        "file_name": config["config_file_name"],
        "sha256": config["config_sha256"],
        "size_bytes": config["config_size_bytes"],
    }
    return merge_targets(
        dataset,
        config["sources"],  # type: ignore[arg-type]
        alias_registry=config["alias_registry"],  # type: ignore[arg-type]
        feature_tables=config["feature_tables"],  # type: ignore[arg-type]
        verify_cif_files=verify_cif_files,
        _config_receipt=config_receipt,
        _config_factory_token=_TARGET_MERGE_FACTORY_TOKEN,
    )


__all__ = [
    "AliasRegistry",
    "CURRENT_FEATURE_TABLES",
    "TARGET_API_VERSION",
    "TargetDataError",
    "TargetMergedDataset",
    "TargetSource",
    "load_target_config",
    "merge_targets",
    "merge_targets_from_config",
]
