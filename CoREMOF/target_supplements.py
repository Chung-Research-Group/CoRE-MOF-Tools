"""Explicit attachment of separately packaged, hash-bound CoRE-MOF targets.

Only the three response values are attached. Full calculation statuses,
diagnostics and provenance remain in the verified optional supplement. This
module performs no imputation, unit conversion, filtering or resampling.
"""

import csv
import io
import json
import math
from pathlib import Path, PurePosixPath
import re
from types import MappingProxyType

from .targets import TargetDataError, TargetSource, _capture_file, _sha256_file, merge_targets


SUPPLEMENT_KIND = "coremof-target-supplement/1.0"
_SUPPLEMENT_FACTORY_TOKEN = object()


def _current_implementation_hashes():
    return {"target_supplements.py": _sha256_file(Path(__file__).resolve())}


_IMPORTED_IMPLEMENTATION_HASHES = MappingProxyType(_current_implementation_hashes())


def _implementation_hashes():
    """Bind the validator that was imported, rejecting later source drift."""

    if _current_implementation_hashes() != dict(_IMPORTED_IMPLEMENTATION_HASHES):
        raise TargetDataError(
            "CoREMOF implementation source changed after module import: "
            "target_supplements.py"
        )
    return _IMPORTED_IMPLEMENTATION_HASHES


_CORE_FILES = (
    "dataset_info.json",
    "metadata/metadata.csv",
    "manifests/cif_manifest.csv",
)
_REQUIRED_FILES = (
    "targets.csv", "targets.jsonl", "targets_schema.json", "target_details.json"
)
_ENDPOINT_DEFINITIONS = {
    "ch4_loading_298K_65bar_mol_kg_framework": {
        "unit": "mol/kg-framework",
        "conditions": {"temperature_K": 298, "pressure_Pa": 6500000},
    },
    "co2_n2_henry_selectivity_298K": {
        "unit": "dimensionless",
        "conditions": {"temperature_K": 298, "pressure_Pa": None},
    },
    "h2_loading_77K_100bar_mol_kg_framework": {
        "unit": "mol/kg-framework",
        "conditions": {"temperature_K": 77, "pressure_Pa": 10000000},
    },
}


def _json(data, name):
    from .targets import TargetDataError

    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise TargetDataError("{} contains duplicate JSON key {!r}".format(name, key))
            result[key] = value
        return result

    def reject_constant(value):
        raise TargetDataError("{} contains non-finite JSON value {}".format(name, value))

    def finite_float(text):
        value = float(text)
        if not math.isfinite(value):
            raise TargetDataError("{} contains an overflowed non-finite JSON number".format(name))
        return value

    try:
        return json.loads(
            data, object_pairs_hook=pairs, parse_constant=reject_constant,
            parse_float=finite_float,
        )
    except (ValueError, UnicodeError) as exc:
        raise TargetDataError("invalid JSON in {}: {}".format(name, exc)) from exc


def _digest(value, name):
    from .targets import TargetDataError

    if not isinstance(value, str) or re.fullmatch(r"[0-9a-f]{64}", value) is None:
        raise TargetDataError("{} must be a lowercase SHA-256 digest".format(name))
    return value


def _bound_path(root, name):
    from .targets import TargetDataError

    if not isinstance(name, str) or not name or "\\" in name:
        raise TargetDataError("supplement file names must be relative POSIX paths")
    relative = PurePosixPath(name)
    if relative.as_posix() != name or relative.is_absolute() or any(part in {".", ".."} for part in relative.parts):
        raise TargetDataError("supplement path escapes its root: {}".format(name))
    candidate = (root / name).resolve()
    if root not in candidate.parents or candidate == root:
        raise TargetDataError("supplement path escapes its root: {}".format(name))
    return candidate


def _parse_csv(data, targets):
    from .targets import TargetDataError

    expected_fields = {"structure_id"}
    expected_fields.update(targets)
    expected_fields.update(name + "__status" for name in targets)
    try:
        reader = csv.DictReader(io.StringIO(data.decode("utf-8-sig"), newline=""))
        fields = reader.fieldnames
        if fields is None or len(fields) != len(set(fields)) or set(fields) != expected_fields:
            raise TargetDataError("targets.csv must contain exactly the ID, value and status columns")
        rows = {}
        for number, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise TargetDataError("targets.csv:{} has an inconsistent row width".format(number))
            identity = row["structure_id"]
            if not identity or identity != identity.strip() or identity in rows:
                raise TargetDataError("targets.csv contains empty, padded or duplicate structure IDs")
            values = {}
            for name in targets:
                text = row[name]
                if text == "":
                    value = None
                else:
                    try:
                        value = float(text)
                    except (ValueError, OverflowError) as exc:
                        raise TargetDataError("targets.csv has a nonnumeric {} value".format(name)) from exc
                    if not math.isfinite(value):
                        raise TargetDataError("targets.csv has a non-finite {} value".format(name))
                values[name] = (value, row[name + "__status"])
            rows[identity] = values
        return rows
    except UnicodeError as exc:
        raise TargetDataError("targets.csv is not UTF-8") from exc


def _parse_jsonl(data, targets, definitions):
    from .targets import TargetDataError

    rows = {}
    for number, line in enumerate(data.splitlines(), start=1):
        if not line.strip():
            continue
        row = _json(line, "targets.jsonl:{}".format(number))
        if not isinstance(row, dict) or set(row) != {"structure_id", "targets"}:
            raise TargetDataError("targets.jsonl rows must contain structure_id and targets")
        identity = row["structure_id"]
        if not isinstance(identity, str) or not identity or identity != identity.strip() or identity in rows:
            raise TargetDataError("targets.jsonl contains empty, padded or duplicate structure IDs")
        payloads = row["targets"]
        if not isinstance(payloads, dict) or set(payloads) != set(targets):
            raise TargetDataError("targets.jsonl has an inconsistent target universe")
        values = {}
        for name in targets:
            payload = payloads[name]
            if not isinstance(payload, dict) or "value" not in payload or "execution_status" not in payload:
                raise TargetDataError("targets.jsonl target payload lacks value or execution_status")
            value = payload["value"]
            if value is not None:
                try:
                    valid_value = type(value) in (int, float) and math.isfinite(value)
                except OverflowError:
                    valid_value = False
                if not valid_value:
                    raise TargetDataError("targets.jsonl has a nonnumeric or non-finite response value")
            status = payload["execution_status"]
            if not isinstance(status, str) or not status or status != status.strip():
                raise TargetDataError("targets.jsonl has an invalid execution_status")
            if status not in {"SUCCESS", "ERROR", "NOT_AVAILABLE"}:
                raise TargetDataError("targets.jsonl has an unsupported canonical target status")
            if type(payload.get("available")) is not bool or payload["available"] != (value is not None):
                raise TargetDataError("targets.jsonl availability disagrees with its response value")
            if (status == "SUCCESS") != (value is not None):
                raise TargetDataError("targets.jsonl execution_status disagrees with its response value")
            if payload.get("unit") != definitions[name]["unit"] or payload.get("conditions") != definitions[name]["conditions"]:
                raise TargetDataError("targets.jsonl unit/conditions disagree with target definitions")
            values[name] = (value, status)
        rows[identity] = values
    return rows


def attach_target_supplement(dataset, supplement_root, *, expected_sha256):
    """Validate a supplement completely, then reuse the audited left join.

    ``expected_sha256`` binds the bytes of the supplied manifest. The merge
    receipt names that manifest digest and records the exact CSV and imported
    supplement-validator source digests.
    Compatibility is defined by the loaded core's original input hashes,
    not by labels, filenames, a shared database version alone, or old aliases.
    """
    from .dataset import CoREMOFDataset, _LEGACY_TARGET_VALUE_COLUMNS

    if not isinstance(dataset, CoREMOFDataset):
        raise TypeError("dataset must be a CoREMOFDataset")
    _implementation_hashes()
    expected_sha256 = _digest(expected_sha256, "expected_sha256")
    root = Path(supplement_root).expanduser().resolve()
    if not root.is_dir():
        raise TargetDataError("target supplement directory does not exist")
    manifest_snapshot = _capture_file(_bound_path(root, "manifest.json"))
    if manifest_snapshot.sha256 != expected_sha256:
        raise TargetDataError("target supplement manifest checksum does not match expected_sha256")
    manifest = _json(manifest_snapshot.data, "manifest.json")
    if not isinstance(manifest, dict) or manifest.get("kind") != SUPPLEMENT_KIND or manifest.get("schema_version") != SUPPLEMENT_KIND:
        raise TargetDataError("unsupported target supplement kind")
    if manifest.get("dataset_version") != dataset.dataset_version:
        raise TargetDataError("target supplement dataset_version does not match the core")
    if type(manifest.get("structure_count")) is not int or manifest["structure_count"] != len(dataset):
        raise TargetDataError("target supplement structure_count does not match the core")
    compatibility = manifest.get("compatible_core_sha256")
    if not isinstance(compatibility, dict) or set(compatibility) != set(_CORE_FILES):
        raise TargetDataError("target supplement must bind all three compatible core input files")
    for name, digest in compatibility.items():
        if _digest(digest, name) != dataset.input_hashes.get(name):
            raise TargetDataError("target supplement is incompatible with core input {}".format(name))
    targets = _LEGACY_TARGET_VALUE_COLUMNS
    definitions = manifest.get("target_definitions")
    if not isinstance(definitions, dict) or set(definitions) != set(targets):
        raise TargetDataError("target supplement must declare exactly the three supported responses")
    for name, definition in definitions.items():
        if not isinstance(definition, dict) or set(definition) != {"unit", "conditions"}:
            raise TargetDataError("target definitions must contain explicit unit and conditions")
        if not isinstance(definition["unit"], str) or not definition["unit"] or not isinstance(definition["conditions"], dict):
            raise TargetDataError("target unit/conditions declarations are invalid")
        expected = _ENDPOINT_DEFINITIONS[name]
        conditions = definition["conditions"]
        if definition["unit"] != expected["unit"] or set(conditions) != set(expected["conditions"]):
            raise TargetDataError("target unit/conditions do not match the fixed response endpoint")
        for condition, expected_value in expected["conditions"].items():
            value = conditions[condition]
            if expected_value is None:
                valid = value is None
            else:
                valid = type(value) in (int, float) and value == expected_value
            if not valid:
                raise TargetDataError("target condition {} does not match the fixed response endpoint".format(condition))
    files = manifest.get("files")
    if not isinstance(files, dict) or not set(_REQUIRED_FILES).issubset(files):
        raise TargetDataError("target supplement manifest lacks required bound files")
    snapshots = {}
    for name, binding in files.items():
        if not isinstance(binding, dict) or set(binding) != {"sha256", "bytes"}:
            raise TargetDataError("supplement file bindings require sha256 and bytes")
        digest = _digest(binding["sha256"], name)
        size = binding["bytes"]
        if type(size) is not int or size < 0:
            raise TargetDataError("supplement file byte counts must be nonnegative integers")
        snapshot = _capture_file(_bound_path(root, name))
        if snapshot.sha256 != digest or snapshot.size_bytes != size:
            raise TargetDataError("target supplement file checksum/size mismatch: {}".format(name))
        snapshots[name] = snapshot
    # Required JSON companions must be valid JSON as well as hash-bound. The
    # full schema/protocol payload is retained, not reinterpreted as model data.
    for name in ("targets_schema.json", "target_details.json"):
        if not isinstance(_json(snapshots[name].data, name), dict):
            raise TargetDataError("{} must contain a JSON object".format(name))
    csv_rows = _parse_csv(snapshots["targets.csv"].data, targets)
    jsonl_rows = _parse_jsonl(snapshots["targets.jsonl"].data, targets, definitions)
    expected_ids = set(dataset.structure_ids)
    if set(csv_rows) != expected_ids or set(jsonl_rows) != expected_ids:
        raise TargetDataError("target supplement must contain the exact core structure-ID universe")
    if csv_rows != jsonl_rows:
        raise TargetDataError("target supplement CSV/JSONL response values or statuses disagree")
    source = TargetSource(
        path=snapshots["targets.csv"].path,
        name="target-supplement:" + manifest_snapshot.sha256,
        target_columns=targets,
        value_types={name: "float" for name in targets},
        units={name: definitions[name]["unit"] for name in targets},
        conditions={name: definitions[name]["conditions"] for name in targets},
    )
    attached = merge_targets(
        dataset, (source,), _supplement_factory_token=_SUPPLEMENT_FACTORY_TOKEN
    )
    source_receipts = attached.receipt().get("sources", ())
    if len(source_receipts) != 1 or source_receipts[0].get("sha256") != snapshots["targets.csv"].sha256:
        raise TargetDataError("target supplement CSV changed between validation and attachment")
    return attached
