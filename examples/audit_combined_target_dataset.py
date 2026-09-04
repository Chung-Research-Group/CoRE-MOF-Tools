#!/usr/bin/env python3
"""Independently audit a combined CoRE-MOF adsorption-target snapshot."""

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
import re
import sys
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple


ENDPOINTS = {
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
SHORT_TO_ENDPOINT = {value["short"]: key for key, value in ENDPOINTS.items()}
EXPECTED_MEMBERS = {
    "BUILD_RECEIPT.json",
    "README.md",
    "coverage_summary.json",
    "target_assignments.csv.gz",
    "targets.json",
    "targets_for_attachment.csv",
}
EXPECTED_TERMINOLOGY = {
    "HISTORICAL_SCIENTIFIC_NULL": (
        "The frozen source status is EXISTING with explicit_null=true and a nonempty "
        "diagnostic; the structure remains eligible, but this endpoint is unavailable "
        "and native-null and cannot be filled because current results may fill only "
        "source keys marked MISSING."
    ),
    "source_diagnostic": (
        "The frozen source diagnostic retained verbatim for audit; for the current "
        "historical scientific-null record it is ZERO_DENOMINATOR."
    ),
}


class AuditError(RuntimeError):
    """Raised when an independent validation fails."""


def _reject_constant(value: str) -> None:
    raise AuditError("non-finite JSON constant {}".format(value))


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise AuditError("duplicate JSON key {}".format(key))
        result[key] = value
    return result


def _loads(text: str, label: str) -> Any:
    try:
        return json.loads(
            text,
            object_pairs_hook=_unique_object,
            parse_constant=_reject_constant,
        )
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        raise AuditError("invalid JSON in {}: {}".format(label, exc)) from exc


def _read_json(path: Path) -> Tuple[Any, bytes]:
    data = path.read_bytes()
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise AuditError("{} is not UTF-8".format(path)) from exc
    return _loads(text, path.name), data


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def _canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def _verify_payload_hash(
    value: Mapping[str, Any],
    field: str,
    label: str,
) -> None:
    expected = value.get(field)
    payload = dict(value)
    payload.pop(field, None)
    if expected != _sha256_bytes(_canonical_json_bytes(payload)):
        raise AuditError("{} payload hash mismatch".format(label))


def _finite(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise AuditError("{} is Boolean".format(label))
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise AuditError("{} is not numeric".format(label)) from exc
    if not math.isfinite(result):
        raise AuditError("{} is not finite".format(label))
    return result


def _parse_expected(values: Sequence[str]) -> Dict[str, int]:
    result: Dict[str, int] = {}
    for value in values:
        if "=" not in value:
            raise AuditError("--expected-finite requires endpoint=count")
        endpoint, raw_count = value.split("=", 1)
        if endpoint not in SHORT_TO_ENDPOINT or endpoint in result:
            raise AuditError("unknown or duplicate expected endpoint {}".format(endpoint))
        try:
            count = int(raw_count)
        except ValueError as exc:
            raise AuditError("invalid expected count {}".format(raw_count)) from exc
        if count < 0:
            raise AuditError("expected count cannot be negative")
        result[endpoint] = count
    if set(result) != set(SHORT_TO_ENDPOINT):
        raise AuditError("expected finite counts must cover ch4, h2, and widom")
    return result


def _parse_expected_sha(values: Sequence[str]) -> Dict[str, str]:
    result: Dict[str, str] = {}
    for value in values:
        if "=" not in value:
            raise AuditError("--expected-input-sha requires logical_name=sha256")
        logical_name, digest = value.split("=", 1)
        digest = digest.lower()
        if logical_name in result:
            raise AuditError("duplicate expected input hash {}".format(logical_name))
        if len(digest) != 64 or any(character not in "0123456789abcdef" for character in digest):
            raise AuditError("invalid expected SHA-256 for {}".format(logical_name))
        result[logical_name] = digest
    return result


def _load_release(path: Path) -> Dict[str, str]:
    data = path.read_bytes()
    try:
        text = data.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise AuditError("release manifest is not UTF-8") from exc
    reader = csv.DictReader(io.StringIO(text, newline=""))
    if reader.fieldnames is None or not {"structure_id", "sha256"}.issubset(reader.fieldnames):
        raise AuditError("invalid release manifest header")
    result: Dict[str, str] = {}
    for row in reader:
        structure_id = str(row["structure_id"]).strip()
        if structure_id in result:
            raise AuditError("duplicate release ID {}".format(structure_id))
        result[structure_id] = str(row["sha256"]).strip().lower()
    return result


def _load_sources(
    sources: Sequence[Tuple[Path, str]],
    release: Mapping[str, str],
) -> Dict[Tuple[str, str], Dict[str, Any]]:
    result: Dict[Tuple[str, str], Dict[str, Any]] = {}
    for path, expected_phase in sources:
        phase_rows: Dict[Tuple[str, str], Dict[str, Any]] = {}
        with path.open("rb") as handle:
            for line_number, raw in enumerate(handle, start=1):
                if not raw.endswith(b"\n"):
                    raise AuditError("source JSONL line lacks newline")
                try:
                    row = _loads(raw.decode("utf-8"), "{}:{}".format(path.name, line_number))
                except UnicodeDecodeError as exc:
                    raise AuditError("source JSONL is not UTF-8") from exc
                if not isinstance(row, dict):
                    raise AuditError("source JSONL row is not an object")
                if row.get("manifest_schema_version") != "coremof-gcmc-manifest/1.0":
                    raise AuditError("source schema mismatch")
                expected_version = (
                    "v26.0.1" if expected_phase == "v26.0.1_base" else "v26.0.2"
                )
                if (
                    row.get("release_phase") != expected_phase
                    or row.get("release_version") != expected_version
                ):
                    raise AuditError("source release identity mismatch")
                if not isinstance(row.get("alias_conflict"), bool):
                    raise AuditError("source alias-conflict flag is not Boolean")
                structure_id = row.get("structure_id")
                endpoint = row.get("endpoint")
                if structure_id not in release or endpoint not in ENDPOINTS:
                    raise AuditError("unknown source assignment")
                if str(row.get("cif_sha256", "")).lower() != release[structure_id]:
                    raise AuditError("source CIF hash mismatch")
                key = (str(structure_id), str(endpoint))
                if key in result:
                    raise AuditError("duplicate source assignment")
                if row.get("status") not in {"EXISTING", "MISSING", "EXCLUDED"}:
                    raise AuditError("invalid source status")
                if row.get("status") == "EXISTING":
                    if row.get("result_semantics") != ENDPOINTS[str(endpoint)]["semantics"]:
                        raise AuditError("historical source semantics mismatch")
                    if row.get("has_existing_result") is not True or not isinstance(
                        row.get("explicit_null"), bool
                    ):
                        raise AuditError("historical source state mismatch")
                    if row["explicit_null"] is True:
                        if row.get("result_value") is not None or not row.get(
                            "result_diagnostic"
                        ):
                            raise AuditError("historical scientific null is invalid")
                    else:
                        _finite(row.get("result_value"), "historical source value")
                elif row.get("status") == "MISSING" and row.get("result_value") is not None:
                    raise AuditError("missing source row exposes a value")
                result[key] = row
                phase_rows[key] = row
        phase_ids = {key[0] for key in phase_rows}
        expected_phase_keys = {
            (structure_id, endpoint)
            for structure_id in phase_ids
            for endpoint in ENDPOINTS
        }
        if set(phase_rows) != expected_phase_keys:
            raise AuditError("source phase lacks complete endpoint triplets")
        for structure_id in phase_ids:
            triplet = [phase_rows[(structure_id, endpoint)] for endpoint in ENDPOINTS]
            excluded = {row["status"] == "EXCLUDED" for row in triplet}
            reasons = {
                json.dumps(
                    row.get("exclusion_reasons", []),
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                    allow_nan=False,
                )
                for row in triplet
            }
            if len(excluded) != 1 or len(reasons) != 1:
                raise AuditError("source phase has inconsistent scientific eligibility")
    if len(result) != len(release) * 3:
        raise AuditError("source assignments do not span the release")
    return result


def _load_current(
    path: Path,
    source: Mapping[Tuple[str, str], Mapping[str, Any]],
    release: Mapping[str, str],
) -> Dict[Tuple[str, str], Dict[str, str]]:
    compressed = path.read_bytes()
    try:
        text = gzip.decompress(compressed).decode("utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise AuditError("invalid current evidence gzip") from exc
    reader = csv.DictReader(io.StringIO(text, newline=""))
    result: Dict[Tuple[str, str], Dict[str, str]] = {}
    for row in reader:
        structure_id = str(row.get("structure_id", ""))
        endpoint = SHORT_TO_ENDPOINT.get(str(row.get("endpoint", "")))
        if structure_id not in release or endpoint is None:
            raise AuditError("unknown current assignment")
        key = (structure_id, endpoint)
        if key in result or source[key].get("status") != "MISSING":
            raise AuditError("current evidence is duplicate or not fill-only")
        expected_phase = str(source[key].get("release_phase", ""))
        expected_version = "v26.0.1" if expected_phase == "v26.0.1_base" else "v26.0.2"
        if (
            row.get("release_phase") != expected_phase
            or row.get("dataset_version") != expected_version
        ):
            raise AuditError("current evidence release identity mismatch")
        if str(row.get("cif_sha256", "")).lower() != release[structure_id]:
            raise AuditError("current evidence CIF hash mismatch")
        contract = ENDPOINTS[endpoint]
        status = row.get("status")
        if status not in {"SUCCESS", "SUCCESS_NULL", "ERROR"}:
            raise AuditError("invalid current status")
        if (
            row.get("target_name") != contract["target"]
            or row.get("target_semantics") != contract["semantics"]
        ):
            raise AuditError("current evidence target identity mismatch")
        if row.get("unit") != contract["unit"] and not (
            status == "ERROR" and row.get("unit") == ""
        ):
            raise AuditError("current evidence unit mismatch")
        if _finite(row.get("temperature_kelvin"), "temperature") != contract["temperature_K"]:
            raise AuditError("current evidence temperature mismatch")
        if _finite(row.get("pressure_pascal"), "pressure") != contract["pressure_Pa"]:
            raise AuditError("current evidence pressure mismatch")
        if status == "SUCCESS":
            if row.get("available") != "true":
                raise AuditError("current SUCCESS is unavailable")
            _finite(row.get("value"), "current value")
        elif status in {"SUCCESS_NULL", "ERROR"}:
            if row.get("available") != "false" or row.get("value") != "":
                raise AuditError("current non-success exposes a value")
        result[key] = row
    return result


def _verify_ledger(root: Path) -> Dict[str, str]:
    ledger_path = root / "SHA256SUMS"
    lines = ledger_path.read_text(encoding="ascii").splitlines()
    entries: Dict[str, str] = {}
    pattern = re.compile(r"^([0-9a-f]{64})  ([A-Za-z0-9_.-]+)$")
    for line in lines:
        match = pattern.fullmatch(line)
        if match is None:
            raise AuditError("malformed checksum line")
        digest, name = match.groups()
        if name in entries:
            raise AuditError("duplicate checksum member")
        entries[name] = digest
    if set(entries) != EXPECTED_MEMBERS:
        raise AuditError("checksum membership mismatch")
    actual_members = {path.name for path in root.iterdir() if path.name != "SHA256SUMS"}
    if actual_members != EXPECTED_MEMBERS:
        raise AuditError("output directory membership mismatch")
    for name, digest in entries.items():
        if _sha256_file(root / name) != digest:
            raise AuditError("checksum mismatch for {}".format(name))
    return entries


def _verify_deterministic_rebuild(
    root: Path,
    comparison_root: Path,
    primary_ledger: Mapping[str, str],
) -> str:
    if not comparison_root.is_dir() or comparison_root.is_symlink():
        raise AuditError("comparison dataset is absent, not a directory, or a symlink")
    comparison_ledger = _verify_ledger(comparison_root)
    if dict(primary_ledger) != comparison_ledger:
        raise AuditError("deterministic rebuild checksum ledgers differ")
    if (root / "SHA256SUMS").read_bytes() != (comparison_root / "SHA256SUMS").read_bytes():
        raise AuditError("deterministic rebuild ledger bytes differ")
    for name in sorted(EXPECTED_MEMBERS):
        if (root / name).read_bytes() != (comparison_root / name).read_bytes():
            raise AuditError("deterministic rebuild bytes differ for {}".format(name))
    return _sha256_file(comparison_root / "SHA256SUMS")


def _load_long(root: Path) -> Tuple[Dict[Tuple[str, str], Dict[str, str]], int]:
    data = (root / "target_assignments.csv.gz").read_bytes()
    try:
        text = gzip.decompress(data).decode("utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise AuditError("invalid target assignments gzip") from exc
    reader = csv.DictReader(io.StringIO(text, newline=""))
    result: Dict[Tuple[str, str], Dict[str, str]] = {}
    prior = None
    for row in reader:
        short = row.get("endpoint")
        endpoint = SHORT_TO_ENDPOINT.get(str(short))
        structure_id = str(row.get("structure_id", ""))
        if endpoint is None:
            raise AuditError("unknown output endpoint")
        key = (structure_id, endpoint)
        sort_key = (structure_id, list(ENDPOINTS).index(endpoint))
        if prior is not None and sort_key <= prior:
            raise AuditError("target assignments are not strictly sorted")
        prior = sort_key
        if key in result:
            raise AuditError("duplicate output assignment")
        result[key] = row
    return result, len(result)


def _expected_status(
    source: Mapping[str, Any],
    current: Optional[Mapping[str, str]],
) -> Tuple[str, str, Optional[float], Optional[float]]:
    status = source.get("status")
    if status == "EXCLUDED":
        return "EXCLUDED_SCIENTIFIC_ELIGIBILITY", "false", None, None
    if status == "EXISTING":
        if source.get("explicit_null") is True:
            return "HISTORICAL_SCIENTIFIC_NULL", "false", None, None
        value = _finite(source.get("result_value"), "historical value")
        raw_error = source.get("result_error")
        error = None if raw_error is None else _finite(raw_error, "historical error")
        return "SUCCESS", "true", value, error
    if current is None:
        return "MISSING_NO_VALIDATED_RESULT_AS_OF_CUTOFF", "false", None, None
    if current["status"] == "SUCCESS":
        error = (
            None
            if current.get("error", "") == ""
            else _finite(current["error"], "current error")
        )
        return "SUCCESS", "true", _finite(current["value"], "current value"), error
    return current["status"], "false", None, None


def _assert_optional_float(actual: str, expected: Optional[float], label: str) -> None:
    if expected is None:
        if actual != "":
            raise AuditError("{} should be null".format(label))
    elif _finite(actual, label) != expected:
        raise AuditError("{} value mismatch".format(label))


def _verify_assignments(
    rows: Mapping[Tuple[str, str], Mapping[str, str]],
    source: Mapping[Tuple[str, str], Mapping[str, Any]],
    current: Mapping[Tuple[str, str], Mapping[str, str]],
    release: Mapping[str, str],
    snapshot_id: str,
) -> Dict[str, Any]:
    if set(rows) != set(source):
        raise AuditError("output assignment key set mismatch")
    available: Dict[str, set] = {value["short"]: set() for value in ENDPOINTS.values()}
    source_counts: Counter = Counter()
    status_counts: Counter = Counter()
    eligible_counts: Counter = Counter()
    excluded_historical: Counter = Counter()
    alias_conflicts: Counter = Counter()
    phase_available: Counter = Counter()
    for key in sorted(source, key=lambda item: (item[0], list(ENDPOINTS).index(item[1]))):
        structure_id, endpoint = key
        output = rows[key]
        source_row = source[key]
        current_row = current.get(key)
        contract = ENDPOINTS[endpoint]
        status, is_available, value, error = _expected_status(source_row, current_row)
        if (
            output.get("schema_version") != "coremof-combined-adsorption-targets/1.0"
            or output.get("snapshot_id") != snapshot_id
            or output.get("release_version") != "v26.0.2"
            or output.get("release_phase") != source_row.get("release_phase")
        ):
            raise AuditError("output release or snapshot identity mismatch")
        if output.get("cif_sha256") != release[structure_id]:
            raise AuditError("output CIF hash mismatch")
        if (
            output.get("target_name") != contract["target"]
            or output.get("target_semantics") != contract["semantics"]
            or output.get("unit") != contract["unit"]
            or _finite(output.get("temperature_K"), "output temperature")
            != contract["temperature_K"]
            or _finite(output.get("pressure_Pa"), "output pressure")
            != contract["pressure_Pa"]
        ):
            raise AuditError("output target contract mismatch")
        if output.get("assignment_status") != status or output.get("available") != is_available:
            raise AuditError("output assignment status mismatch")
        expected_eligibility = "EXCLUDED" if source_row.get("status") == "EXCLUDED" else "ELIGIBLE"
        if output.get("eligibility_status") != expected_eligibility:
            raise AuditError("output eligibility mismatch")
        expected_reasons = json.dumps(
            source_row.get("exclusion_reasons", []),
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        if output.get("eligibility_reasons") != expected_reasons:
            raise AuditError("output eligibility reasons mismatch")
        expected_historical = "true" if source_row.get("has_existing_result") is True else "false"
        if output.get("historical_candidate_present") != expected_historical:
            raise AuditError("output historical-candidate flag mismatch")
        expected_alias_conflict = "true" if source_row.get("alias_conflict") else "false"
        if output.get("source_alias_conflict") != expected_alias_conflict:
            raise AuditError("output source alias-conflict flag mismatch")
        if output.get("source_diagnostic") != str(
            source_row.get("result_diagnostic", "")
        ):
            raise AuditError("output source diagnostic mismatch")
        if source_row.get("status") == "EXCLUDED":
            expected_source_class = "EXCLUDED"
        elif source_row.get("status") == "EXISTING":
            expected_source_class = "HISTORICAL_EXISTING"
        elif current_row is None:
            expected_source_class = "NONE"
        elif current_row.get("dataset_version") == "v26.0.2":
            expected_source_class = "CURRENT_VALIDATED_V2602"
        else:
            expected_source_class = "CURRENT_VALIDATED_V2601"
        if output.get("source_class") != expected_source_class:
            raise AuditError("output source class mismatch")
        if current_row is not None:
            current_pairs = {
                "source_dataset_version": "dataset_version",
                "protocol_charge_method": "protocol_charge_method",
                "cif_preprocessing_action": "cif_preprocessing_action",
                "raspa_version": "raspa_version",
                "task_id": "task_id",
                "task_payload_sha256": "task_payload_sha256",
                "result_sha256": "result_sha256",
                "receipt_file_sha256": "receipt_file_sha256",
                "receipt_payload_sha256": "receipt_payload_sha256",
                "source_batch": "source_batch",
                "job_id": "job_id",
                "validation_scope": "validation_scope",
            }
            for output_field, current_field in current_pairs.items():
                if output.get(output_field, "") != current_row.get(current_field, ""):
                    raise AuditError("output current provenance mismatch")
            if output.get("source_logical_name", "") != current_row.get("source_batch", ""):
                raise AuditError("output current logical-source mismatch")
        elif source_row.get("status") == "EXISTING":
            historical_expected = {
                "source_dataset_version": str(source_row.get("release_version", "")),
                "source_logical_name": Path(
                    str(source_row.get("selected_source_file", ""))
                ).name,
                "source_sha256": str(source_row.get("selected_source_sha256", "")),
                "source_record_id": str(source_row.get("selected_legacy_id", "")),
                "validation_scope": "FROZEN_COMPLETION_MANIFEST_ACCEPTED_HISTORICAL",
            }
            for field, expected_value in historical_expected.items():
                if output.get(field, "") != expected_value:
                    raise AuditError("output historical provenance mismatch")
        elif source_row.get("status") == "EXCLUDED":
            candidate_expected = {
                "source_logical_name": (
                    Path(str(source_row.get("selected_source_file", ""))).name
                    if source_row.get("has_existing_result") is True
                    else ""
                ),
                "source_sha256": (
                    str(source_row.get("selected_source_sha256", ""))
                    if source_row.get("has_existing_result") is True
                    else ""
                ),
                "source_record_id": (
                    str(source_row.get("selected_legacy_id", ""))
                    if source_row.get("has_existing_result") is True
                    else ""
                ),
            }
            for field, expected_value in candidate_expected.items():
                if output.get(field, "") != expected_value:
                    raise AuditError("output excluded provenance mismatch")
        else:
            blank_fields = (
                "source_dataset_version",
                "source_logical_name",
                "source_sha256",
                "source_record_id",
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
            if any(output.get(field, "") != "" for field in blank_fields):
                raise AuditError("missing output exposes unsupported provenance")
        _assert_optional_float(str(output.get("value", "")), value, "target")
        _assert_optional_float(str(output.get("reported_error", "")), error, "reported error")
        short = str(contract["short"])
        status_counts[(short, status)] += 1
        if source_row.get("alias_conflict") is True:
            alias_conflicts[(short, status)] += 1
        if expected_eligibility == "ELIGIBLE":
            eligible_counts[short] += 1
        if is_available == "true":
            available[short].add(structure_id)
            source_counts[(short, str(output.get("source_class")))] += 1
            phase_available[(short, str(output.get("release_phase")))] += 1
        if (
            source_row.get("status") == "EXCLUDED"
            and source_row.get("has_existing_result") is True
        ):
            excluded_historical[short] += 1
    all_sets = list(available.values())
    any_ids = set().union(*all_sets)
    all_ids = set.intersection(*all_sets)
    release_phases = sorted(
        {str(row.get("release_phase", "")) for row in source.values()}
    )
    endpoint_summary = {}
    for short in sorted(available):
        endpoint_summary[short] = {
            "finite_unique_structures": len(available[short]),
            "eligible_structures": eligible_counts[short],
            "source_counts": {
                source_class: count
                for (endpoint, source_class), count in sorted(source_counts.items())
                if endpoint == short
            },
            "status_counts": {
                status: count
                for (endpoint, status), count in sorted(status_counts.items())
                if endpoint == short
            },
            "phase_counts": {
                phase: phase_available[(short, phase)] for phase in release_phases
            },
            "excluded_historical_candidate_count": excluded_historical[short],
            "alias_conflict_status_counts": {
                status: count
                for (endpoint, status), count in sorted(alias_conflicts.items())
                if endpoint == short
            },
        }
    pairwise = {}
    shorts = [str(ENDPOINTS[key]["short"]) for key in ENDPOINTS]
    for left_index, left in enumerate(shorts):
        for right in shorts[left_index + 1 :]:
            pairwise["{}_and_{}".format(left, right)] = len(
                available[left].intersection(available[right])
            )
    return {
        "assignment_count": len(rows),
        "finite_endpoint_assignments": sum(len(value) for value in available.values()),
        "structures_with_any_finite_target": len(any_ids),
        "structures_with_all_three_finite_targets": len(all_ids),
        "structures_with_no_finite_target": len(release) - len(any_ids),
        "pairwise_finite_intersections": pairwise,
        "endpoint": endpoint_summary,
    }


def _verify_wide(
    root: Path,
    release: Mapping[str, str],
    long_rows: Mapping[Tuple[str, str], Mapping[str, str]],
) -> None:
    with (root / "targets_for_attachment.csv").open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        observed_ids = []
        for row in reader:
            structure_id = str(row.get("structure_id", ""))
            observed_ids.append(structure_id)
            if structure_id not in release or row.get("cif_sha256") != release[structure_id]:
                raise AuditError("wide release identity mismatch")
            for endpoint, contract in ENDPOINTS.items():
                assignment = long_rows[(structure_id, endpoint)]
                target = str(contract["target"])
                pairs = {
                    target: "value",
                    target + "__reported_error": "reported_error",
                    target + "__status": "assignment_status",
                    target + "__source_class": "source_class",
                }
                for wide_field, long_field in pairs.items():
                    if row.get(wide_field, "") != assignment.get(long_field, ""):
                        raise AuditError("wide/long mismatch for {}".format(wide_field))
                expected_eligible = (
                    "true"
                    if assignment["eligibility_status"] == "ELIGIBLE"
                    else "false"
                )
                if row.get(target + "__eligible") != expected_eligible:
                    raise AuditError("wide eligibility mismatch")
    if observed_ids != sorted(release):
        raise AuditError("wide IDs are not the complete sorted release")


def _verify_coverage(
    root: Path,
    summary: Mapping[str, Any],
    release_count: int,
    receipt: Mapping[str, Any],
) -> None:
    coverage, _ = _read_json(root / "coverage_summary.json")
    if not isinstance(coverage, dict):
        raise AuditError("coverage summary is not an object")
    expected_endpoint = {}
    for short, endpoint in summary["endpoint"].items():
        finite = endpoint["finite_unique_structures"]
        eligible = endpoint["eligible_structures"]
        expected_endpoint[short] = {
            "release_structures": release_count,
            "eligible_structures": eligible,
            "excluded_structures": release_count - eligible,
            "finite_unique_structures": finite,
            "eligible_without_finite_target": eligible - finite,
            "full_release_coverage_fraction": finite / release_count,
            "eligible_coverage_fraction": finite / eligible,
            "finite_by_release_phase": endpoint["phase_counts"],
            "finite_by_source_class": endpoint["source_counts"],
            "assignment_status_counts": endpoint["status_counts"],
            "excluded_rows_with_historical_candidate": endpoint[
                "excluded_historical_candidate_count"
            ],
            "source_alias_conflict_rows_by_assignment_status": endpoint[
                "alias_conflict_status_counts"
            ],
        }
    expected_coverage = {
        "schema_version": "coremof-combined-adsorption-targets/1.0",
        "release_structure_count": release_count,
        "finite_endpoint_assignments": summary["finite_endpoint_assignments"],
        "structures_with_any_finite_target": summary[
            "structures_with_any_finite_target"
        ],
        "structures_with_all_three_finite_targets": summary[
            "structures_with_all_three_finite_targets"
        ],
        "structures_with_no_finite_target": summary["structures_with_no_finite_target"],
        "pairwise_finite_intersections": summary["pairwise_finite_intersections"],
        "endpoint": expected_endpoint,
        "snapshot_id": receipt["snapshot_id"],
        "data_cutoff_utc": receipt["data_cutoff_utc"],
        "status": "PASS_COMBINED_CURRENT_AVAILABLE_TARGET_DATASET",
        "campaign_complete": False,
        "official_split": False,
        "promotion_performed": False,
        "publication_authorized": False,
        "redistribution_scope": "PRIVATE_INTERNAL_LOCAL_BUILD",
    }
    if coverage != expected_coverage:
        for key in sorted(set(coverage).union(expected_coverage)):
            if coverage.get(key) != expected_coverage.get(key):
                raise AuditError(
                    "coverage {} differs from independent recomputation: {} != {}".format(
                        key,
                        coverage.get(key),
                        expected_coverage.get(key),
                    )
                )
        raise AuditError("coverage summary differs from independent recomputation")


def audit(args: argparse.Namespace) -> Dict[str, Any]:
    root = args.dataset.resolve()
    if not root.is_dir() or root.is_symlink():
        raise AuditError("dataset root is absent, not a directory, or a symlink")
    ledger = _verify_ledger(root)
    receipt, receipt_bytes = _read_json(root / "BUILD_RECEIPT.json")
    if not isinstance(receipt, dict):
        raise AuditError("build receipt is not an object")
    payload_hash = receipt.get("receipt_payload_sha256")
    payload = dict(receipt)
    payload.pop("receipt_payload_sha256", None)
    if payload_hash != _sha256_bytes(_canonical_json_bytes(payload)):
        raise AuditError("build receipt payload hash mismatch")
    if receipt.get("status") != "PASS_COMBINED_CURRENT_AVAILABLE_TARGET_DATASET":
        raise AuditError("build receipt status is not passing")
    if (
        receipt.get("release_version") != "v26.0.2"
        or receipt.get("release_structure_count") != args.expected_release_count
        or not isinstance(receipt.get("snapshot_id"), str)
        or not receipt.get("snapshot_id")
    ):
        raise AuditError("build receipt release or snapshot identity mismatch")
    nonpublication_flags = (
        "campaign_complete",
        "official_split",
        "promotion_performed",
        "publication_authorized",
    )
    if any(receipt.get(flag) is not False for flag in nonpublication_flags):
        raise AuditError("build receipt overclaims completion or publication")
    if receipt.get("terminology") != EXPECTED_TERMINOLOGY:
        raise AuditError("build receipt terminology definitions mismatch")
    receipt_output_bindings = receipt.get("output_bindings")
    expected_output_binding_names = EXPECTED_MEMBERS.difference(
        {"BUILD_RECEIPT.json"}
    )
    if (
        not isinstance(receipt_output_bindings, dict)
        or set(receipt_output_bindings) != expected_output_binding_names
    ):
        raise AuditError("build receipt output-binding membership mismatch")
    for name in sorted(expected_output_binding_names):
        binding = receipt_output_bindings[name]
        if (
            not isinstance(binding, dict)
            or binding.get("sha256") != ledger[name]
            or binding.get("size_bytes") != (root / name).stat().st_size
        ):
            raise AuditError("build receipt output binding mismatch for {}".format(name))

    expected_inputs = {
        "release_manifest": args.release_manifest,
        "base_manifest": args.base_manifest,
        "base_summary": args.base_summary,
        "additions_manifest": args.additions_manifest,
        "additions_summary": args.additions_summary,
        "current_evidence": args.current_evidence,
        "current_build_receipt": args.current_build_receipt,
        "current_independent_audit": args.current_independent_audit,
    }
    input_bindings = receipt.get("input_bindings")
    if not isinstance(input_bindings, dict):
        raise AuditError("build receipt lacks input bindings")
    input_hashes = {}
    for name, path in expected_inputs.items():
        digest = _sha256_file(path)
        input_hashes[name] = digest
        binding = input_bindings.get(name)
        if not isinstance(binding, dict) or binding.get("sha256") != digest:
            raise AuditError("input binding mismatch for {}".format(name))
    for logical_name, expected_sha in _parse_expected_sha(args.expected_input_sha).items():
        if logical_name not in input_hashes:
            raise AuditError("unknown expected input hash name {}".format(logical_name))
        if input_hashes[logical_name] != expected_sha:
            raise AuditError("expected input SHA-256 mismatch for {}".format(logical_name))

    current_receipt, current_receipt_bytes = _read_json(args.current_build_receipt)
    current_audit, _ = _read_json(args.current_independent_audit)
    if not isinstance(current_receipt, dict) or not isinstance(current_audit, dict):
        raise AuditError("current receipt or audit is not an object")
    if current_receipt.get("status") != "PASS_CURRENT_FINISHED_EXPLORATORY_TARGET_SNAPSHOT":
        raise AuditError("current build receipt is not passing")
    if current_audit.get("status") != "PASS_INDEPENDENT_CURRENT_FINISHED_TARGET_SNAPSHOT":
        raise AuditError("current independent audit is not passing")
    _verify_payload_hash(current_receipt, "receipt_payload_sha256", "current receipt")
    _verify_payload_hash(current_audit, "audit_payload_sha256", "current audit")
    current_audit_inputs = current_audit.get("input_hashes")
    if not isinstance(current_audit_inputs, dict):
        raise AuditError("current audit lacks input hashes")
    if (
        current_audit_inputs.get("BUILD_RECEIPT.json")
        != _sha256_bytes(current_receipt_bytes)
        or current_audit_inputs.get("target_evidence.csv.gz")
        != input_hashes["current_evidence"]
        or current_audit_inputs.get("receipt_payload_sha256")
        != current_receipt.get("receipt_payload_sha256")
    ):
        raise AuditError("current audit input binding mismatch")

    release = _load_release(args.release_manifest)
    if len(release) != args.expected_release_count:
        raise AuditError(
            "release contains {} structures, expected {}".format(
                len(release), args.expected_release_count
            )
        )
    source = _load_sources(
        (
            (args.base_manifest, "v26.0.1_base"),
            (args.additions_manifest, "v26.0.2_additions"),
        ),
        release,
    )
    current = _load_current(args.current_evidence, source, release)
    current_finite_by_endpoint = {
        short: sum(
            row["endpoint"] == short and row["status"] == "SUCCESS"
            for row in current.values()
        )
        for short in SHORT_TO_ENDPOINT
    }
    current_output_counts = current_receipt.get("output_counts")
    if (
        not isinstance(current_output_counts, dict)
        or current_output_counts.get("worker_evidence_rows") != len(current)
        or current_output_counts.get("finite_target_values")
        != sum(current_finite_by_endpoint.values())
        or current_output_counts.get("finite_target_values_by_endpoint")
        != current_finite_by_endpoint
    ):
        raise AuditError("current receipt counts do not match evidence")
    long_rows, long_count = _load_long(root)
    if long_count != args.expected_release_count * len(ENDPOINTS):
        raise AuditError("long assignment count mismatch")
    summary = _verify_assignments(
        long_rows,
        source,
        current,
        release,
        str(receipt["snapshot_id"]),
    )
    _verify_wide(root, release, long_rows)
    _verify_coverage(root, summary, len(release), receipt)
    coverage_document, _ = _read_json(root / "coverage_summary.json")
    if receipt.get("output_counts") != coverage_document:
        raise AuditError("build receipt output counts do not match coverage summary")

    targets, _ = _read_json(root / "targets.json")
    expected_targets = [ENDPOINTS[key]["target"] for key in ENDPOINTS]
    expected_types = {target: "float" for target in expected_targets}
    expected_units = {
        ENDPOINTS[key]["target"]: ENDPOINTS[key]["unit"] for key in ENDPOINTS
    }
    expected_conditions = {
        ENDPOINTS[key]["target"]: {
            "method": ENDPOINTS[key]["method"],
            "quantity": ENDPOINTS[key]["semantics"],
            "temperature_K": ENDPOINTS[key]["temperature_K"],
            "pressure_Pa": ENDPOINTS[key]["pressure_Pa"],
        }
        for key in ENDPOINTS
    }
    expected_target_configuration = {
        "sources": [
            {
                "path": "targets_for_attachment.csv",
                "name": "combined_current_available_adsorption_targets",
                "format": "csv",
                "id_column": "structure_id",
                "target_columns": expected_targets,
                "value_types": expected_types,
                "units": expected_units,
                "conditions": expected_conditions,
                "null_values": [""],
            }
        ],
    }
    if targets != expected_target_configuration:
        raise AuditError("target configuration differs from the exact contract")

    expected_finite = _parse_expected(
        args.expected_finite
        or ("ch4=28979", "h2=28974", "widom=28944")
    )
    actual_finite = {
        short: value["finite_unique_structures"]
        for short, value in summary["endpoint"].items()
    }
    if actual_finite != expected_finite:
        raise AuditError("published-cutoff finite counts mismatch: {}".format(actual_finite))
    if summary["finite_endpoint_assignments"] != args.expected_finite_total:
        raise AuditError("published-cutoff finite assignment total mismatch")

    comparison_ledger_sha256 = None
    if args.comparison_dataset is not None:
        comparison_ledger_sha256 = _verify_deterministic_rebuild(
            root,
            args.comparison_dataset.resolve(),
            ledger,
        )

    checks = {
        "checksum_membership_and_hashes": True,
        "receipt_payload_hash": True,
        "receipt_input_bindings": True,
        "release_exact_expected_ids_and_hashes": True,
        "complete_release_times_three_assignment_key_set": True,
        "historical_existing_values_exact": True,
        "current_results_fill_only_missing_keys": True,
        "excluded_targets_remain_null": True,
        "nonfinite_values_rejected": True,
        "wide_long_exact_parity": True,
        "coverage_recomputed": True,
        "target_configuration_exact": True,
        "completion_and_publication_flags_false": True,
        "byte_identical_deterministic_rebuild": args.comparison_dataset is not None,
    }
    result: Dict[str, Any] = {
        "schema_version": "coremof-combined-adsorption-targets-independent-audit/1.0",
        "audit_id": args.audit_id,
        "audited_dataset": root.name,
        "status": "PASS_INDEPENDENT_COMBINED_CURRENT_AVAILABLE_TARGET_DATASET",
        "checks": checks,
        "input_hashes": input_hashes,
        "dataset_ledger_sha256": _sha256_file(root / "SHA256SUMS"),
        "dataset_receipt_sha256": _sha256_bytes(receipt_bytes),
        "dataset_files": ledger,
        "comparison_dataset_ledger_sha256": comparison_ledger_sha256,
        "counts": summary,
        "redistribution_scope": "PRIVATE_INTERNAL_LOCAL_BUILD",
        "publication_authorized": False,
        "auditor": {
            "logical_name": Path(__file__).name,
            "sha256": _sha256_file(Path(__file__).resolve()),
        },
    }
    result["audit_payload_sha256"] = _sha256_bytes(_canonical_json_bytes(result))
    return result


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--release-manifest", type=Path, required=True)
    parser.add_argument("--base-manifest", type=Path, required=True)
    parser.add_argument("--base-summary", type=Path, required=True)
    parser.add_argument("--additions-manifest", type=Path, required=True)
    parser.add_argument("--additions-summary", type=Path, required=True)
    parser.add_argument("--current-evidence", type=Path, required=True)
    parser.add_argument("--current-build-receipt", type=Path, required=True)
    parser.add_argument("--current-independent-audit", type=Path, required=True)
    parser.add_argument("--audit-output", type=Path, required=True)
    parser.add_argument("--audit-id", required=True)
    parser.add_argument("--comparison-dataset", type=Path)
    parser.add_argument("--expected-release-count", type=int, default=42_574)
    parser.add_argument(
        "--expected-finite",
        action="append",
        default=None,
        metavar="ENDPOINT=COUNT",
    )
    parser.add_argument("--expected-finite-total", type=int, default=86_897)
    parser.add_argument(
        "--expected-input-sha",
        action="append",
        default=[],
        metavar="LOGICAL_NAME=SHA256",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    try:
        args = _parser().parse_args(argv)
        if args.audit_output.exists():
            raise AuditError("refusing to overwrite {}".format(args.audit_output))
        result = audit(args)
        data = (
            json.dumps(result, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False)
            + "\n"
        ).encode("utf-8")
        args.audit_output.parent.mkdir(parents=True, exist_ok=True)
        with args.audit_output.open("xb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
    except (AuditError, OSError) as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 2
    print(str(args.audit_output.resolve()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
