#!/usr/bin/env python3
"""Build a deterministic, fill-only CoRE-MOF adsorption-target snapshot.

The builder combines three evidence classes without changing cohort or split
membership:

* reusable historical values classified ``EXISTING`` by a frozen completion
  manifest;
* collector-validated results for keys classified ``MISSING`` by that same
  manifest; and
* collector-validated results for later release additions.

Rows classified ``EXCLUDED`` by the completion manifest remain explicit nulls.
No value is imputed, averaged, clipped, or allowed to overwrite a historical
assignment.  The output is one immutable as-of-cutoff attachment dataset, not
a claim that the calculation campaign or dataset publication is complete.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import shutil
import sys
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


SCHEMA_VERSION = "coremof-combined-adsorption-targets/1.0"
RECEIPT_SCHEMA_VERSION = "coremof-combined-adsorption-targets-receipt/1.0"
HISTORICAL_SCIENTIFIC_NULL_DEFINITION = (
    "The frozen source status is EXISTING with explicit_null=true and a nonempty "
    "diagnostic; the structure remains eligible, but this endpoint is unavailable "
    "and native-null and cannot be filled because current results may fill only "
    "source keys marked MISSING."
)
SOURCE_DIAGNOSTIC_DEFINITION = (
    "The frozen source diagnostic retained verbatim for audit; for the current "
    "historical scientific-null record it is ZERO_DENOMINATOR."
)

ENDPOINTS: Mapping[str, Mapping[str, Any]] = {
    "CH4": {
        "short": "ch4",
        "target": "ch4_loading_298K_65bar_mol_kg_framework",
        "semantics": "absolute_loading_mol_kg_framework",
        "unit": "mol/kg-framework",
        "temperature_K": 298.0,
        "pressure_Pa": 6_500_000.0,
        "method": "grand-canonical Monte Carlo",
    },
    "H2": {
        "short": "h2",
        "target": "h2_loading_77K_100bar_mol_kg_framework",
        "semantics": "absolute_loading_mol_kg_framework",
        "unit": "mol/kg-framework",
        "temperature_K": 77.0,
        "pressure_Pa": 10_000_000.0,
        "method": "grand-canonical Monte Carlo",
    },
    "WIDOM_CO2_N2": {
        "short": "widom",
        "target": "co2_n2_widom_rosenbluth_weight_ratio_298K_1bar",
        "semantics": "widom_rosenbluth_weight_ratio_co2_over_n2",
        "unit": "dimensionless",
        "temperature_K": 298.0,
        "pressure_Pa": 100_000.0,
        "method": "four-component Widom insertion",
    },
}
SHORT_TO_ENDPOINT = {str(value["short"]): key for key, value in ENDPOINTS.items()}

LONG_FIELDS = (
    "schema_version",
    "snapshot_id",
    "release_version",
    "release_phase",
    "structure_id",
    "cif_sha256",
    "endpoint",
    "target_name",
    "target_semantics",
    "eligibility_status",
    "eligibility_reasons",
    "assignment_status",
    "available",
    "value",
    "reported_error",
    "unit",
    "temperature_K",
    "pressure_Pa",
    "source_class",
    "source_alias_conflict",
    "source_diagnostic",
    "source_dataset_version",
    "source_logical_name",
    "source_sha256",
    "source_record_id",
    "historical_candidate_present",
    "protocol_charge_method",
    "cif_preprocessing_action",
    "raspa_version",
    "task_id",
    "task_payload_sha256",
    "result_sha256",
    "receipt_file_sha256",
    "receipt_payload_sha256",
    "source_batch",
    "job_id",
    "validation_scope",
)


class BuildError(RuntimeError):
    """Raised when an input or merge invariant fails closed."""


def _reject_constant(value: str) -> None:
    raise BuildError("non-finite JSON constant is forbidden: {}".format(value))


def _unique_object(pairs: Sequence[Tuple[str, Any]]) -> Dict[str, Any]:
    output: Dict[str, Any] = {}
    for key, value in pairs:
        if key in output:
            raise BuildError("duplicate JSON key: {}".format(key))
        output[key] = value
    return output


def _loads_json(text: str, label: str) -> Any:
    try:
        return json.loads(
            text,
            object_pairs_hook=_unique_object,
            parse_constant=_reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError, TypeError, ValueError) as exc:
        raise BuildError("invalid JSON in {}: {}".format(label, exc)) from exc


def _read_json(path: Path) -> Tuple[Any, bytes]:
    data = path.read_bytes()
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise BuildError("{} is not strict UTF-8".format(path.name)) from exc
    return _loads_json(text, path.name), data


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _verify_payload_hash(
    value: Mapping[str, Any],
    field: str,
    label: str,
) -> None:
    expected = value.get(field)
    payload = dict(value)
    payload.pop(field, None)
    actual = _sha256(_canonical_json_bytes(payload))
    if expected != actual:
        raise BuildError("{} payload hash mismatch".format(label))


def _strict_bool(value: str, label: str) -> bool:
    if value == "true":
        return True
    if value == "false":
        return False
    raise BuildError("{} must be literal true or false".format(label))


def _finite_float(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise BuildError("{} is not numeric".format(label))
    try:
        converted = float(value)
    except (TypeError, ValueError) as exc:
        raise BuildError("{} is not numeric".format(label)) from exc
    if not math.isfinite(converted):
        raise BuildError("{} is not finite".format(label))
    return converted


def _float_text(value: Optional[float]) -> str:
    return "" if value is None else repr(float(value))


def _json_text(value: Any) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def _canonical_json_bytes(value: Any) -> bytes:
    return (_json_text(value) + "\n").encode("utf-8")


def _pretty_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
            indent=2,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def _write_bytes(path: Path, data: bytes) -> None:
    with path.open("xb") as handle:
        handle.write(data)
        handle.flush()
        os.fsync(handle.fileno())


def _write_json(path: Path, value: Any) -> None:
    _write_bytes(path, _pretty_json_bytes(value))


def _write_csv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
        handle.flush()
        os.fsync(handle.fileno())


def _write_gzip_csv(
    path: Path,
    fields: Sequence[str],
    rows: Iterable[Mapping[str, Any]],
) -> None:
    with path.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text_handle:
                writer = csv.DictWriter(
                    text_handle,
                    fieldnames=list(fields),
                    lineterminator="\n",
                )
                writer.writeheader()
                for row in rows:
                    writer.writerow({field: row.get(field, "") for field in fields})
                text_handle.flush()
        raw.flush()
        os.fsync(raw.fileno())


def _load_release_manifest(path: Path) -> Tuple[Dict[str, Dict[str, str]], str, int]:
    data = path.read_bytes()
    try:
        text = data.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise BuildError("release manifest is not UTF-8") from exc
    reader = csv.DictReader(io.StringIO(text, newline=""))
    required = {"structure_id", "sha256"}
    if reader.fieldnames is None or not required.issubset(reader.fieldnames):
        raise BuildError("release manifest lacks structure_id or sha256")
    rows: Dict[str, Dict[str, str]] = {}
    for index, row in enumerate(reader, start=2):
        structure_id = str(row["structure_id"]).strip()
        sha = str(row["sha256"]).strip().lower()
        if not structure_id or len(sha) != 64 or any(c not in "0123456789abcdef" for c in sha):
            raise BuildError("invalid release manifest row {}".format(index))
        if structure_id in rows:
            raise BuildError("duplicate release structure_id {}".format(structure_id))
        rows[structure_id] = dict(row)
    return rows, _sha256(data), len(data)


def _load_source_manifest(
    path: Path,
    release_ids: Mapping[str, Mapping[str, str]],
    expected_phase: str,
) -> Tuple[Dict[Tuple[str, str], Dict[str, Any]], str, int, Counter]:
    digest = hashlib.sha256()
    output: Dict[Tuple[str, str], Dict[str, Any]] = {}
    status_counts: Counter = Counter()
    size_bytes = 0
    with path.open("rb") as handle:
        for line_number, raw in enumerate(handle, start=1):
            digest.update(raw)
            size_bytes += len(raw)
            if not raw.endswith(b"\n"):
                raise BuildError("{} line {} lacks newline".format(path.name, line_number))
            try:
                text = raw.decode("utf-8")
            except UnicodeDecodeError as exc:
                raise BuildError("{} line {} is not UTF-8".format(path.name, line_number)) from exc
            row = _loads_json(text, "{}:{}".format(path.name, line_number))
            if not isinstance(row, dict):
                raise BuildError("{} line {} is not an object".format(path.name, line_number))
            structure_id = row.get("structure_id")
            endpoint = row.get("endpoint")
            if not isinstance(structure_id, str) or structure_id not in release_ids:
                raise BuildError("unknown structure at {}:{}".format(path.name, line_number))
            if endpoint not in ENDPOINTS:
                raise BuildError("unknown endpoint at {}:{}".format(path.name, line_number))
            if row.get("manifest_schema_version") != "coremof-gcmc-manifest/1.0":
                raise BuildError("wrong source schema at {}:{}".format(path.name, line_number))
            if row.get("release_phase") != expected_phase:
                raise BuildError("wrong release phase at {}:{}".format(path.name, line_number))
            expected_release_version = (
                "v26.0.1" if expected_phase == "v26.0.1_base" else "v26.0.2"
            )
            if row.get("release_version") != expected_release_version:
                raise BuildError("wrong release version at {}:{}".format(path.name, line_number))
            if not isinstance(row.get("alias_conflict"), bool):
                raise BuildError(
                    "source alias-conflict flag is not Boolean at {}:{}".format(
                        path.name, line_number
                    )
                )
            expected_cif_sha = str(release_ids[structure_id]["sha256"]).lower()
            if str(row.get("cif_sha256", "")).lower() != expected_cif_sha:
                raise BuildError("CIF hash mismatch at {}:{}".format(path.name, line_number))
            key = (structure_id, str(endpoint))
            if key in output:
                raise BuildError("duplicate source assignment {} {}".format(*key))
            status = row.get("status")
            if status not in {"EXISTING", "MISSING", "EXCLUDED"}:
                raise BuildError("invalid source status at {}:{}".format(path.name, line_number))
            if status == "EXISTING":
                if row.get("result_semantics") != ENDPOINTS[str(endpoint)]["semantics"]:
                    raise BuildError(
                        "historical result semantics mismatch at {}:{}".format(
                            path.name, line_number
                        )
                    )
                if row.get("has_existing_result") is not True or not isinstance(
                    row.get("explicit_null"), bool
                ):
                    raise BuildError(
                        "inconsistent EXISTING row at {}:{}".format(path.name, line_number)
                    )
                if row["explicit_null"] is True:
                    if row.get("result_value") is not None or not row.get(
                        "result_diagnostic"
                    ):
                        raise BuildError(
                            "invalid historical scientific null at {}:{}".format(
                                path.name, line_number
                            )
                        )
                else:
                    try:
                        _finite_float(row.get("result_value"), "historical result_value")
                    except BuildError as exc:
                        raise BuildError(
                            "{} at {}:{} for {} {}".format(
                                exc,
                                path.name,
                                line_number,
                                structure_id,
                                endpoint,
                            )
                        ) from exc
            elif status == "MISSING":
                if (
                    not isinstance(row.get("explicit_null"), bool)
                    or not isinstance(row.get("has_existing_result"), bool)
                    or row.get("result_value") is not None
                ):
                    raise BuildError(
                        "inconsistent MISSING row at {}:{}".format(path.name, line_number)
                    )
            output[key] = row
            status_counts[(str(endpoint), str(status))] += 1
    structure_ids = {key[0] for key in output}
    expected_keys = {
        (structure_id, endpoint)
        for structure_id in structure_ids
        for endpoint in ENDPOINTS
    }
    if set(output) != expected_keys:
        raise BuildError("{} does not contain all three endpoints per structure".format(path.name))
    for structure_id in structure_ids:
        rows = [output[(structure_id, endpoint)] for endpoint in ENDPOINTS]
        excluded = {row["status"] == "EXCLUDED" for row in rows}
        reasons = {
            _json_text(row.get("exclusion_reasons", []))
            for row in rows
        }
        if len(excluded) != 1 or len(reasons) != 1:
            raise BuildError(
                "{} has endpoint-inconsistent scientific eligibility for {}".format(
                    path.name, structure_id
                )
            )
    return output, digest.hexdigest(), size_bytes, status_counts


def _verify_source_summary(
    summary_path: Path,
    manifest_path: Path,
    manifest_sha: str,
    rows: Mapping[Tuple[str, str], Mapping[str, Any]],
    status_counts: Counter,
) -> Tuple[Mapping[str, Any], str, int]:
    summary, data = _read_json(summary_path)
    if not isinstance(summary, dict):
        raise BuildError("{} is not an object".format(summary_path.name))
    artifact_sha = (
        summary.get("artifacts", {}).get("jsonl", {}).get("sha256")
        if isinstance(summary.get("artifacts"), dict)
        else None
    )
    if artifact_sha != manifest_sha:
        raise BuildError("{} does not bind {}".format(summary_path.name, manifest_path.name))
    if summary.get("record_count") != len(rows):
        raise BuildError("{} record_count mismatch".format(summary_path.name))
    structure_ids = {key[0] for key in rows}
    excluded_ids = {
        structure_id
        for structure_id in structure_ids
        if rows[(structure_id, next(iter(ENDPOINTS)))]["status"] == "EXCLUDED"
    }
    alias_conflict_rows = sum(row.get("alias_conflict") is True for row in rows.values())
    alias_conflict_ids = {
        structure_id
        for (structure_id, _), row in rows.items()
        if row.get("alias_conflict") is True
    }
    explicit_null_rows = sum(row.get("explicit_null") is True for row in rows.values())
    expected_summary_scalars = {
        "structure_count": len(structure_ids),
        "eligible_structure_count": len(structure_ids) - len(excluded_ids),
        "excluded_structure_count": len(excluded_ids),
        "alias_conflict_record_count": alias_conflict_rows,
        "alias_conflict_structure_count": len(alias_conflict_ids),
        "explicit_null_record_count": explicit_null_rows,
    }
    for key, expected_value in expected_summary_scalars.items():
        if summary.get(key) != expected_value:
            raise BuildError("{} {} mismatch".format(summary_path.name, key))
    expected_overall_status = summary.get("status_counts")
    if not isinstance(expected_overall_status, dict):
        raise BuildError("{} lacks overall status counts".format(summary_path.name))
    overall_status_counts = {}
    for status in ("EXCLUDED", "EXISTING", "MISSING"):
        count = sum(
            value
            for (_, observed), value in status_counts.items()
            if observed == status
        )
        if count or status in expected_overall_status:
            overall_status_counts[status] = count
    if expected_overall_status != overall_status_counts:
        raise BuildError("{} overall status counts mismatch".format(summary_path.name))
    expected_endpoint_counts = summary.get("endpoint_status_counts")
    if not isinstance(expected_endpoint_counts, dict):
        raise BuildError("{} lacks endpoint_status_counts".format(summary_path.name))
    for endpoint in ENDPOINTS:
        expected = expected_endpoint_counts.get(endpoint)
        if not isinstance(expected, dict):
            raise BuildError("{} lacks counts for {}".format(summary_path.name, endpoint))
        actual = {
            status: status_counts[(endpoint, status)]
            for status in ("EXCLUDED", "EXISTING", "MISSING")
            if status_counts[(endpoint, status)] or status in expected
        }
        if actual != expected:
            raise BuildError(
                "{} endpoint counts mismatch for {}".format(summary_path.name, endpoint)
            )
    return summary, _sha256(data), len(data)


def _load_current_evidence(
    path: Path,
    release_ids: Mapping[str, Mapping[str, str]],
    source_rows: Mapping[Tuple[str, str], Mapping[str, Any]],
) -> Tuple[Dict[Tuple[str, str], Dict[str, str]], str, int, Counter]:
    compressed = path.read_bytes()
    try:
        decompressed = gzip.decompress(compressed).decode("utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise BuildError("current evidence is not deterministic UTF-8 gzip CSV") from exc
    reader = csv.DictReader(io.StringIO(decompressed, newline=""))
    required = {
        "dataset_version",
        "release_phase",
        "structure_id",
        "cif_sha256",
        "endpoint",
        "target_name",
        "target_semantics",
        "status",
        "available",
        "value",
        "error",
        "unit",
        "temperature_kelvin",
        "pressure_pascal",
    }
    if reader.fieldnames is None or not required.issubset(reader.fieldnames):
        raise BuildError("current evidence lacks required columns")
    output: Dict[Tuple[str, str], Dict[str, str]] = {}
    counts: Counter = Counter()
    for line_number, row in enumerate(reader, start=2):
        structure_id = str(row["structure_id"]).strip()
        short_endpoint = str(row["endpoint"]).strip()
        endpoint = SHORT_TO_ENDPOINT.get(short_endpoint)
        if structure_id not in release_ids or endpoint is None:
            raise BuildError("unknown current-evidence key at line {}".format(line_number))
        key = (structure_id, endpoint)
        if key in output:
            raise BuildError("duplicate current-evidence key {} {}".format(*key))
        source = source_rows.get(key)
        if source is None or source.get("status") != "MISSING":
            raise BuildError("current result is not fill-only for {} {}".format(*key))
        expected_phase = str(source.get("release_phase", ""))
        expected_version = "v26.0.1" if expected_phase == "v26.0.1_base" else "v26.0.2"
        if row["release_phase"] != expected_phase or row["dataset_version"] != expected_version:
            raise BuildError(
                "current-evidence release identity mismatch at line {}".format(line_number)
            )
        if str(row["cif_sha256"]).lower() != str(release_ids[structure_id]["sha256"]).lower():
            raise BuildError("current-evidence CIF hash mismatch at line {}".format(line_number))
        contract = ENDPOINTS[endpoint]
        status = str(row["status"])
        if status not in {"SUCCESS", "SUCCESS_NULL", "ERROR"}:
            raise BuildError("unsupported current status {}".format(status))
        if row.get("target_name") != contract["target"]:
            raise BuildError(
                "target-name mismatch at current-evidence line {}".format(line_number)
            )
        if row.get("target_semantics") != contract["semantics"]:
            raise BuildError(
                "target-semantics mismatch at current-evidence line {}".format(line_number)
            )
        observed_unit = str(row["unit"])
        if observed_unit != str(contract["unit"]) and not (
            status == "ERROR" and observed_unit == ""
        ):
            raise BuildError(
                "unit mismatch at current-evidence line {} for {} {}: {!r} != {!r}".format(
                    line_number,
                    structure_id,
                    endpoint,
                    observed_unit,
                    contract["unit"],
                )
            )
        if _finite_float(row["temperature_kelvin"], "temperature") != contract["temperature_K"]:
            raise BuildError(
                "temperature mismatch at current-evidence line {}".format(line_number)
            )
        if _finite_float(row["pressure_pascal"], "pressure") != contract["pressure_Pa"]:
            raise BuildError("pressure mismatch at current-evidence line {}".format(line_number))
        available = _strict_bool(str(row["available"]), "available")
        if status == "SUCCESS":
            if not available:
                raise BuildError("SUCCESS is unavailable at line {}".format(line_number))
            _finite_float(row["value"], "current result value")
            if row.get("error", "") != "":
                _finite_float(row["error"], "current result error")
        elif status in {"SUCCESS_NULL", "ERROR"}:
            if available or row["value"] != "":
                raise BuildError("non-success exposes value at line {}".format(line_number))
        output[key] = dict(row)
        counts[(short_endpoint, status)] += 1
        counts[(str(row["dataset_version"]), status)] += 1
    return output, _sha256(compressed), len(compressed), counts


def _logical_source_name(path_text: Any) -> str:
    if not isinstance(path_text, str) or not path_text:
        return ""
    return Path(path_text).name


def _build_assignment(
    snapshot_id: str,
    structure_id: str,
    endpoint: str,
    release_row: Mapping[str, str],
    source: Mapping[str, Any],
    current: Optional[Mapping[str, str]],
) -> Dict[str, Any]:
    contract = ENDPOINTS[endpoint]
    source_status = str(source["status"])
    exclusions = source.get("exclusion_reasons", [])
    if not isinstance(exclusions, list) or any(not isinstance(item, str) for item in exclusions):
        raise BuildError("invalid exclusion reasons for {} {}".format(structure_id, endpoint))
    eligible = source_status != "EXCLUDED"
    row: Dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "snapshot_id": snapshot_id,
        "release_version": "v26.0.2",
        "release_phase": str(source.get("release_phase", "")),
        "structure_id": structure_id,
        "cif_sha256": str(release_row["sha256"]).lower(),
        "endpoint": str(contract["short"]),
        "target_name": str(contract["target"]),
        "target_semantics": str(contract["semantics"]),
        "eligibility_status": "ELIGIBLE" if eligible else "EXCLUDED",
        "eligibility_reasons": _json_text(exclusions),
        "assignment_status": "",
        "available": "false",
        "value": "",
        "reported_error": "",
        "unit": str(contract["unit"]),
        "temperature_K": _float_text(float(contract["temperature_K"])),
        "pressure_Pa": _float_text(float(contract["pressure_Pa"])),
        "source_class": "NONE",
        "source_alias_conflict": "true" if source.get("alias_conflict") else "false",
        "source_diagnostic": str(source.get("result_diagnostic", "")),
        "source_dataset_version": "",
        "source_logical_name": "",
        "source_sha256": "",
        "source_record_id": "",
        "historical_candidate_present": (
            "true" if source.get("has_existing_result") is True else "false"
        ),
        "protocol_charge_method": "",
        "cif_preprocessing_action": "",
        "raspa_version": "",
        "task_id": "",
        "task_payload_sha256": "",
        "result_sha256": "",
        "receipt_file_sha256": "",
        "receipt_payload_sha256": "",
        "source_batch": "",
        "job_id": "",
        "validation_scope": "",
    }
    if not eligible:
        row["assignment_status"] = "EXCLUDED_SCIENTIFIC_ELIGIBILITY"
        row["source_class"] = "EXCLUDED"
        if source.get("has_existing_result") is True:
            row["source_logical_name"] = _logical_source_name(source.get("selected_source_file"))
            row["source_sha256"] = str(source.get("selected_source_sha256", ""))
            row["source_record_id"] = str(source.get("selected_legacy_id", ""))
        return row

    if source_status == "EXISTING":
        row.update(
            {
                "source_class": "HISTORICAL_EXISTING",
                "source_dataset_version": str(source.get("release_version", "v26.0.1")),
                "source_logical_name": _logical_source_name(source.get("selected_source_file")),
                "source_sha256": str(source.get("selected_source_sha256", "")),
                "source_record_id": str(source.get("selected_legacy_id", "")),
                "validation_scope": "FROZEN_COMPLETION_MANIFEST_ACCEPTED_HISTORICAL",
            }
        )
        if source.get("explicit_null") is True:
            row["assignment_status"] = "HISTORICAL_SCIENTIFIC_NULL"
            return row
        value = _finite_float(source.get("result_value"), "historical result")
        error_value = source.get("result_error")
        error = None if error_value is None else _finite_float(error_value, "historical error")
        row.update(
            {
                "assignment_status": "SUCCESS",
                "available": "true",
                "value": _float_text(value),
                "reported_error": _float_text(error),
            }
        )
        return row

    if source_status != "MISSING":
        raise BuildError("unexpected eligible source status")
    if current is None:
        row["assignment_status"] = "MISSING_NO_VALIDATED_RESULT_AS_OF_CUTOFF"
        return row

    status = str(current["status"])
    source_class = (
        "CURRENT_VALIDATED_V2602"
        if current.get("dataset_version") == "v26.0.2"
        else "CURRENT_VALIDATED_V2601"
    )
    row.update(
        {
            "assignment_status": status,
            "source_class": source_class,
            "source_dataset_version": str(current.get("dataset_version", "")),
            "source_logical_name": str(current.get("source_batch", "")),
            "protocol_charge_method": str(current.get("protocol_charge_method", "")),
            "cif_preprocessing_action": str(current.get("cif_preprocessing_action", "")),
            "raspa_version": str(current.get("raspa_version", "")),
            "task_id": str(current.get("task_id", "")),
            "task_payload_sha256": str(current.get("task_payload_sha256", "")),
            "result_sha256": str(current.get("result_sha256", "")),
            "receipt_file_sha256": str(current.get("receipt_file_sha256", "")),
            "receipt_payload_sha256": str(current.get("receipt_payload_sha256", "")),
            "source_batch": str(current.get("source_batch", "")),
            "job_id": str(current.get("job_id", "")),
            "validation_scope": str(current.get("validation_scope", "")),
        }
    )
    if status == "SUCCESS":
        row["available"] = "true"
        row["value"] = _float_text(_finite_float(current["value"], "current value"))
        if current.get("error", "") != "":
            row["reported_error"] = _float_text(
                _finite_float(current["error"], "current error")
            )
    return row


def _coverage(assignments: Sequence[Mapping[str, Any]], release_count: int) -> Dict[str, Any]:
    by_endpoint: Dict[str, Any] = {}
    available_sets: Dict[str, set] = {}
    total_finite = 0
    for endpoint in ENDPOINTS:
        short = str(ENDPOINTS[endpoint]["short"])
        rows = [row for row in assignments if row["endpoint"] == short]
        if len(rows) != release_count:
            raise BuildError("endpoint {} does not span release".format(short))
        available_ids = {str(row["structure_id"]) for row in rows if row["available"] == "true"}
        eligible = sum(row["eligibility_status"] == "ELIGIBLE" for row in rows)
        available_sets[short] = available_ids
        total_finite += len(available_ids)
        phase_counts: Dict[str, int] = {}
        for phase in sorted({str(row["release_phase"]) for row in rows}):
            phase_counts[phase] = sum(
                row["available"] == "true" and row["release_phase"] == phase for row in rows
            )
        source_counts = Counter(
            str(row["source_class"]) for row in rows if row["available"] == "true"
        )
        status_counts = Counter(str(row["assignment_status"]) for row in rows)
        excluded_historical = sum(
            row["eligibility_status"] == "EXCLUDED"
            and row["historical_candidate_present"] == "true"
            for row in rows
        )
        alias_conflict_counts = Counter(
            str(row["assignment_status"])
            for row in rows
            if row["source_alias_conflict"] == "true"
        )
        by_endpoint[short] = {
            "release_structures": release_count,
            "eligible_structures": eligible,
            "excluded_structures": release_count - eligible,
            "finite_unique_structures": len(available_ids),
            "eligible_without_finite_target": eligible - len(available_ids),
            "full_release_coverage_fraction": len(available_ids) / release_count,
            "eligible_coverage_fraction": len(available_ids) / eligible,
            "finite_by_release_phase": phase_counts,
            "finite_by_source_class": dict(sorted(source_counts.items())),
            "assignment_status_counts": dict(sorted(status_counts.items())),
            "excluded_rows_with_historical_candidate": excluded_historical,
            "source_alias_conflict_rows_by_assignment_status": dict(
                sorted(alias_conflict_counts.items())
            ),
        }
    any_ids = set().union(*available_sets.values())
    all_ids = set.intersection(*available_sets.values())
    pairwise = {}
    shorts = [str(ENDPOINTS[key]["short"]) for key in ENDPOINTS]
    for left_index, left in enumerate(shorts):
        for right in shorts[left_index + 1 :]:
            pairwise["{}_and_{}".format(left, right)] = len(
                available_sets[left].intersection(available_sets[right])
            )
    return {
        "schema_version": SCHEMA_VERSION,
        "release_structure_count": release_count,
        "finite_endpoint_assignments": total_finite,
        "structures_with_any_finite_target": len(any_ids),
        "structures_with_all_three_finite_targets": len(all_ids),
        "structures_with_no_finite_target": release_count - len(any_ids),
        "pairwise_finite_intersections": pairwise,
        "endpoint": by_endpoint,
    }


def _wide_rows(
    release_ids: Sequence[str],
    assignments_by_key: Mapping[Tuple[str, str], Mapping[str, Any]],
) -> Tuple[List[str], List[Dict[str, Any]]]:
    fields = ["structure_id", "cif_sha256"]
    for endpoint in ENDPOINTS:
        target = str(ENDPOINTS[endpoint]["target"])
        fields.extend(
            [
                target,
                target + "__reported_error",
                target + "__status",
                target + "__source_class",
                target + "__eligible",
            ]
        )
    rows: List[Dict[str, Any]] = []
    for structure_id in release_ids:
        first = assignments_by_key[(structure_id, next(iter(ENDPOINTS)))]
        row: Dict[str, Any] = {
            "structure_id": structure_id,
            "cif_sha256": first["cif_sha256"],
        }
        for endpoint in ENDPOINTS:
            assignment = assignments_by_key[(structure_id, endpoint)]
            target = str(ENDPOINTS[endpoint]["target"])
            row[target] = assignment["value"]
            row[target + "__reported_error"] = assignment["reported_error"]
            row[target + "__status"] = assignment["assignment_status"]
            row[target + "__source_class"] = assignment["source_class"]
            row[target + "__eligible"] = (
                "true" if assignment["eligibility_status"] == "ELIGIBLE" else "false"
            )
        rows.append(row)
    return fields, rows


def _target_config() -> Dict[str, Any]:
    target_columns = [str(ENDPOINTS[key]["target"]) for key in ENDPOINTS]
    return {
        "sources": [
            {
                "path": "targets_for_attachment.csv",
                "name": "combined_current_available_adsorption_targets",
                "format": "csv",
                "id_column": "structure_id",
                "target_columns": target_columns,
                "value_types": {target: "float" for target in target_columns},
                "units": {
                    str(ENDPOINTS[key]["target"]): str(ENDPOINTS[key]["unit"])
                    for key in ENDPOINTS
                },
                "conditions": {
                    str(ENDPOINTS[key]["target"]): {
                        "method": ENDPOINTS[key]["method"],
                        "quantity": ENDPOINTS[key]["semantics"],
                        "temperature_K": ENDPOINTS[key]["temperature_K"],
                        "pressure_Pa": ENDPOINTS[key]["pressure_Pa"],
                    }
                    for key in ENDPOINTS
                },
                "null_values": [""],
            }
        ],
    }


def _parse_expected(values: Sequence[str]) -> Dict[str, int]:
    output: Dict[str, int] = {}
    for value in values:
        if "=" not in value:
            raise BuildError("--expected-finite requires endpoint=count")
        endpoint, raw_count = value.split("=", 1)
        if endpoint not in SHORT_TO_ENDPOINT or endpoint in output:
            raise BuildError("unknown or duplicate expected endpoint {}".format(endpoint))
        try:
            count = int(raw_count)
        except ValueError as exc:
            raise BuildError("invalid expected count {}".format(raw_count)) from exc
        if count < 0:
            raise BuildError("expected count cannot be negative")
        output[endpoint] = count
    return output


def _parse_expected_sha(values: Sequence[str]) -> Dict[str, str]:
    output: Dict[str, str] = {}
    for value in values:
        if "=" not in value:
            raise BuildError("--expected-input-sha requires logical_name=sha256")
        logical_name, digest = value.split("=", 1)
        digest = digest.lower()
        if logical_name in output:
            raise BuildError("duplicate expected input hash {}".format(logical_name))
        if len(digest) != 64 or any(character not in "0123456789abcdef" for character in digest):
            raise BuildError("invalid expected SHA-256 for {}".format(logical_name))
        output[logical_name] = digest
    return output


def build(args: argparse.Namespace) -> Path:
    output = args.output_directory.resolve()
    if output.exists():
        raise BuildError("refusing to overwrite existing output {}".format(output))
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = output.with_name(".{}.staging.{}".format(output.name, os.getpid()))
    if staging.exists():
        raise BuildError("staging path already exists {}".format(staging))

    release, release_sha, release_size = _load_release_manifest(args.release_manifest)
    if len(release) != args.expected_release_count:
        raise BuildError(
            "release count {} != {}".format(len(release), args.expected_release_count)
        )

    base_rows, base_sha, base_size, base_status = _load_source_manifest(
        args.base_manifest,
        release,
        "v26.0.1_base",
    )
    additions_rows, additions_sha, additions_size, additions_status = _load_source_manifest(
        args.additions_manifest,
        release,
        "v26.0.2_additions",
    )
    base_summary, base_summary_sha, base_summary_size = _verify_source_summary(
        args.base_summary,
        args.base_manifest,
        base_sha,
        base_rows,
        base_status,
    )
    additions_summary, additions_summary_sha, additions_summary_size = _verify_source_summary(
        args.additions_summary,
        args.additions_manifest,
        additions_sha,
        additions_rows,
        additions_status,
    )
    if base_summary.get("structure_count") != args.expected_base_count:
        raise BuildError("base summary structure count mismatch")
    if additions_summary.get("structure_count") != args.expected_additions_count:
        raise BuildError("additions summary structure count mismatch")
    source_rows = dict(base_rows)
    overlap = set(source_rows).intersection(additions_rows)
    if overlap:
        raise BuildError("base/additions assignment overlap")
    source_rows.update(additions_rows)
    expected_assignment_count = len(release) * len(ENDPOINTS)
    if len(source_rows) != expected_assignment_count:
        raise BuildError(
            "source manifests span {} assignments; expected {}".format(
                len(source_rows), expected_assignment_count
            )
        )

    current_receipt, current_receipt_bytes = _read_json(args.current_build_receipt)
    if not isinstance(current_receipt, dict):
        raise BuildError("current build receipt is not an object")
    if current_receipt.get("status") != "PASS_CURRENT_FINISHED_EXPLORATORY_TARGET_SNAPSHOT":
        raise BuildError("current build receipt is not passing")
    _verify_payload_hash(
        current_receipt,
        "receipt_payload_sha256",
        "current build receipt",
    )
    current_nonpublication_flags = (
        "campaign_complete",
        "official_split",
        "promotion_performed",
        "publication_authorized",
    )
    if any(current_receipt.get(flag) is not False for flag in current_nonpublication_flags):
        raise BuildError("current build receipt overclaims completion or publication")
    current_audit, current_audit_bytes = _read_json(args.current_independent_audit)
    if not isinstance(current_audit, dict):
        raise BuildError("current independent audit is not an object")
    if current_audit.get("status") != "PASS_INDEPENDENT_CURRENT_FINISHED_TARGET_SNAPSHOT":
        raise BuildError("current independent audit is not passing")
    _verify_payload_hash(
        current_audit,
        "audit_payload_sha256",
        "current independent audit",
    )
    current, current_sha, current_size, current_counts = _load_current_evidence(
        args.current_evidence,
        release,
        source_rows,
    )
    audit_hashes = current_audit.get("input_hashes")
    if not isinstance(audit_hashes, dict):
        audit_hashes = {}
    expected_current_sha = audit_hashes.get("target_evidence.csv.gz")
    if expected_current_sha != current_sha:
        raise BuildError("current independent audit does not bind target evidence")
    if audit_hashes.get("BUILD_RECEIPT.json") != _sha256(current_receipt_bytes):
        raise BuildError("current independent audit does not bind build receipt")
    if audit_hashes.get("receipt_payload_sha256") != current_receipt.get(
        "receipt_payload_sha256"
    ):
        raise BuildError("current independent audit does not bind receipt payload")

    current_output_counts = current_receipt.get("output_counts")
    if not isinstance(current_output_counts, dict):
        raise BuildError("current build receipt lacks output counts")
    current_finite_by_endpoint = {
        short: current_counts[(short, "SUCCESS")] for short in SHORT_TO_ENDPOINT
    }
    if current_output_counts.get("worker_evidence_rows") != len(current):
        raise BuildError("current receipt evidence-row count mismatch")
    if current_output_counts.get("finite_target_values") != sum(
        current_finite_by_endpoint.values()
    ):
        raise BuildError("current receipt finite-value count mismatch")
    if current_output_counts.get("finite_target_values_by_endpoint") != (
        current_finite_by_endpoint
    ):
        raise BuildError("current receipt endpoint count mismatch")

    actual_input_hashes = {
        "release_manifest": release_sha,
        "base_manifest": base_sha,
        "base_summary": base_summary_sha,
        "additions_manifest": additions_sha,
        "additions_summary": additions_summary_sha,
        "current_evidence": current_sha,
        "current_build_receipt": _sha256(current_receipt_bytes),
        "current_independent_audit": _sha256(current_audit_bytes),
    }
    for logical_name, expected_sha in _parse_expected_sha(args.expected_input_sha).items():
        if logical_name not in actual_input_hashes:
            raise BuildError("unknown expected input hash name {}".format(logical_name))
        if actual_input_hashes[logical_name] != expected_sha:
            raise BuildError("input SHA-256 mismatch for {}".format(logical_name))

    release_ids = sorted(release)
    assignments: List[Dict[str, Any]] = []
    assignments_by_key: Dict[Tuple[str, str], Dict[str, Any]] = {}
    for structure_id in release_ids:
        for endpoint in ENDPOINTS:
            key = (structure_id, endpoint)
            assignment = _build_assignment(
                args.snapshot_id,
                structure_id,
                endpoint,
                release[structure_id],
                source_rows[key],
                current.get(key),
            )
            assignments.append(assignment)
            assignments_by_key[key] = assignment

    coverage = _coverage(assignments, len(release))
    expected_finite = _parse_expected(args.expected_finite)
    for endpoint, expected in expected_finite.items():
        actual = coverage["endpoint"][endpoint]["finite_unique_structures"]
        if actual != expected:
            raise BuildError(
                "finite {} count {} != expected {}".format(endpoint, actual, expected)
            )
    if args.expected_finite_total is not None:
        if coverage["finite_endpoint_assignments"] != args.expected_finite_total:
            raise BuildError("total finite assignment count mismatch")

    wide_fields, wide = _wide_rows(release_ids, assignments_by_key)
    try:
        staging.mkdir(mode=0o750)
        _write_csv(staging / "targets_for_attachment.csv", wide_fields, wide)
        _write_gzip_csv(staging / "target_assignments.csv.gz", LONG_FIELDS, assignments)
        _write_json(staging / "targets.json", _target_config())
        coverage.update(
            {
                "snapshot_id": args.snapshot_id,
                "data_cutoff_utc": args.data_cutoff_utc,
                "status": "PASS_COMBINED_CURRENT_AVAILABLE_TARGET_DATASET",
                "campaign_complete": False,
                "official_split": False,
                "promotion_performed": False,
                "publication_authorized": False,
                "redistribution_scope": "PRIVATE_INTERNAL_LOCAL_BUILD",
            }
        )
        _write_json(staging / "coverage_summary.json", coverage)
        readme = """# Combined current-available CoRE-MOF v26.0.2 adsorption targets

In this dataset, **final as of the recorded cutoff** means an immutable,
fill-only union of every accepted input bound by this receipt at that cutoff.
It does not mean that the calculation campaign is complete, that every
published structure has a target, or that publication/redistribution is
authorized.

The canonical attachment table has one row for each of the 42,574 published
v26.0.2 structure identifiers.  It combines reusable historical values marked
`EXISTING` by the frozen completion manifests with collector-validated results
for keys those same manifests marked `MISSING`.  A current result may fill only
its corresponding missing key; it cannot overwrite an existing value.

Rows excluded because current MOFChecker reports a lone molecule or current
Chen--Manz reports an isolated structure retain explicit exclusion status and
native null targets.  Missing, error, physical-null, held, failed-before-science,
and not-yet-calculated cases also remain null.  No value is imputed, averaged,
clipped, or inferred.

`HISTORICAL_SCIENTIFIC_NULL` means the frozen source status is `EXISTING` with
`explicit_null=true` and a nonempty diagnostic.  The structure remains eligible,
but that endpoint is unavailable and native-null; it cannot be filled because
current results may fill only source keys marked `MISSING`.  The long-form
`source_diagnostic` field retains that upstream explanation verbatim (currently
`ZERO_DENOMINATOR` for the one historical Widom record in this snapshot).

The three endpoint contracts are methane absolute loading at 298 K and 65 bar,
hydrogen absolute loading at 77 K and 100 bar, and the dimensionless raw CO2/N2
Rosenbluth-weight ratio from four-component Widom insertion at 298 K and a
1-bar simulation state.  The ratio is not a Henry coefficient or adsorption
selectivity.

Files:

- `targets_for_attachment.csv`: 42,574-row exact-ID left-join input.
- `target_assignments.csv.gz`: complete 42,574 x 3 assignment/status/provenance
  accounting.
- `targets.json`: typed CoREMOF-tools deferred-attachment configuration.
- `coverage_summary.json`: exact finite, null, exclusion, source, and endpoint
  intersection counts.
- `BUILD_RECEIPT.json`: immutable input, method, policy, and output bindings.
- `SHA256SUMS`: checksum ledger for every file in this directory except the
  ledger itself.

Target values and their availability never participate in checker consensus,
parent grouping, leakage blocks, diversity balancing, cohort construction, or
train/validation/test assignment.  Attach this table only after those choices
are frozen.  This structure-resolved internal dataset can contain CSD-derived
values and is not cleared for a public Git repository.
"""
        _write_bytes(staging / "README.md", readme.encode("utf-8"))

        output_hashes = {
            name: {
                "sha256": _sha256_file(staging / name),
                "size_bytes": (staging / name).stat().st_size,
            }
            for name in (
                "README.md",
                "coverage_summary.json",
                "target_assignments.csv.gz",
                "targets.json",
                "targets_for_attachment.csv",
            )
        }
        receipt: Dict[str, Any] = {
            "schema_version": RECEIPT_SCHEMA_VERSION,
            "snapshot_id": args.snapshot_id,
            "status": "PASS_COMBINED_CURRENT_AVAILABLE_TARGET_DATASET",
            "created_at_utc": args.data_cutoff_utc,
            "data_cutoff_utc": args.data_cutoff_utc,
            "final_as_of_cutoff_definition": (
                "Immutable fill-only union of all accepted inputs bound by this receipt "
                "at the cutoff; not campaign completion or publication authorization."
            ),
            "release_version": "v26.0.2",
            "release_structure_count": len(release),
            "merge_policy": {
                "historical_existing_reused": True,
                "current_results_fill_only_manifest_missing_keys": True,
                "overwrite_existing": False,
                "conflicts": "FAIL_CLOSED",
                "imputation": False,
                "excluded_values_exposed_as_targets": False,
                "target_independent_assignments_unchanged": True,
            },
            "terminology": {
                "HISTORICAL_SCIENTIFIC_NULL": (
                    HISTORICAL_SCIENTIFIC_NULL_DEFINITION
                ),
                "source_diagnostic": SOURCE_DIAGNOSTIC_DEFINITION,
            },
            "endpoint_contracts": {key: dict(value) for key, value in ENDPOINTS.items()},
            "input_bindings": {
                "release_manifest": {
                    "logical_name": args.release_manifest.name,
                    "sha256": release_sha,
                    "size_bytes": release_size,
                    "row_count": len(release),
                },
                "base_manifest": {
                    "logical_name": args.base_manifest.name,
                    "sha256": base_sha,
                    "size_bytes": base_size,
                    "row_count": len(base_rows),
                },
                "base_summary": {
                    "logical_name": args.base_summary.name,
                    "sha256": base_summary_sha,
                    "size_bytes": base_summary_size,
                },
                "additions_manifest": {
                    "logical_name": args.additions_manifest.name,
                    "sha256": additions_sha,
                    "size_bytes": additions_size,
                    "row_count": len(additions_rows),
                },
                "additions_summary": {
                    "logical_name": args.additions_summary.name,
                    "sha256": additions_summary_sha,
                    "size_bytes": additions_summary_size,
                },
                "current_evidence": {
                    "logical_name": args.current_evidence.name,
                    "sha256": current_sha,
                    "size_bytes": current_size,
                    "row_count": len(current),
                },
                "current_build_receipt": {
                    "logical_name": args.current_build_receipt.name,
                    "sha256": _sha256(current_receipt_bytes),
                    "size_bytes": len(current_receipt_bytes),
                },
                "current_independent_audit": {
                    "logical_name": args.current_independent_audit.name,
                    "sha256": _sha256(current_audit_bytes),
                    "size_bytes": len(current_audit_bytes),
                },
            },
            "source_summary_bindings": {
                "base": {
                    "eligible_structure_count": base_summary["eligible_structure_count"],
                    "excluded_structure_count": base_summary["excluded_structure_count"],
                    "alias_conflict_record_count": base_summary[
                        "alias_conflict_record_count"
                    ],
                    "alias_conflict_structure_count": base_summary[
                        "alias_conflict_structure_count"
                    ],
                    "endpoint_status_counts": base_summary["endpoint_status_counts"],
                },
                "additions": {
                    "eligible_structure_count": additions_summary["eligible_structure_count"],
                    "excluded_structure_count": additions_summary["excluded_structure_count"],
                    "alias_conflict_record_count": additions_summary[
                        "alias_conflict_record_count"
                    ],
                    "alias_conflict_structure_count": additions_summary[
                        "alias_conflict_structure_count"
                    ],
                    "endpoint_status_counts": additions_summary["endpoint_status_counts"],
                },
                "current_status_counts": {
                    "{}|{}".format(key[0], key[1]): value
                    for key, value in sorted(current_counts.items())
                },
            },
            "output_counts": coverage,
            "output_bindings": output_hashes,
            "campaign_complete": False,
            "official_split": False,
            "promotion_performed": False,
            "publication_authorized": False,
            "redistribution_scope": "PRIVATE_INTERNAL_LOCAL_BUILD",
            "builder": {
                "logical_name": Path(__file__).name,
                "sha256": _sha256_file(Path(__file__).resolve()),
            },
        }
        payload = dict(receipt)
        receipt["receipt_payload_sha256"] = _sha256(_canonical_json_bytes(payload))
        _write_json(staging / "BUILD_RECEIPT.json", receipt)

        ledger_names = sorted(path.name for path in staging.iterdir() if path.name != "SHA256SUMS")
        ledger = "".join(
            "{}  {}\n".format(_sha256_file(staging / name), name) for name in ledger_names
        )
        _write_bytes(staging / "SHA256SUMS", ledger.encode("ascii"))
        os.replace(str(staging), str(output))
    except BaseException:
        if staging.exists():
            shutil.rmtree(staging)
        raise
    return output


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-manifest", type=Path, required=True)
    parser.add_argument("--base-manifest", type=Path, required=True)
    parser.add_argument("--base-summary", type=Path, required=True)
    parser.add_argument("--additions-manifest", type=Path, required=True)
    parser.add_argument("--additions-summary", type=Path, required=True)
    parser.add_argument("--current-evidence", type=Path, required=True)
    parser.add_argument("--current-build-receipt", type=Path, required=True)
    parser.add_argument("--current-independent-audit", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--snapshot-id", required=True)
    parser.add_argument("--data-cutoff-utc", required=True)
    parser.add_argument("--expected-release-count", type=int, default=42_574)
    parser.add_argument("--expected-base-count", type=int, default=36_628)
    parser.add_argument("--expected-additions-count", type=int, default=5_946)
    parser.add_argument(
        "--expected-finite",
        action="append",
        default=[],
        metavar="ENDPOINT=COUNT",
        help="Fail unless a short endpoint (ch4, h2, widom) has this finite-ID count.",
    )
    parser.add_argument("--expected-finite-total", type=int)
    parser.add_argument(
        "--expected-input-sha",
        action="append",
        default=[],
        metavar="LOGICAL_NAME=SHA256",
        help="Fail unless a named input has this exact SHA-256.",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    try:
        args = _parser().parse_args(argv)
        output = build(args)
    except (BuildError, OSError) as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 2
    print(str(output))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
