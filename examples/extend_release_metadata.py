"""Add accepted targets to a new private release copy, without changing structures.

Requires a complete, checksum-bound source release, the separately verified
combined-target snapshot, and jsonschema for record validation. The output is
always a private candidate. This does not authorize publication, alter labels
or groups, regenerate features, or change any benchmark assignments.
"""
from __future__ import annotations

import argparse
from collections import Counter
import copy
import csv
import hashlib
import io
import json
import math
from pathlib import Path, PurePosixPath
import re
import tempfile

from CoREMOF._transactions import publish_directory
from CoREMOF.dataset import CoREMOFDataset


BUILDER_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


ENDPOINTS = {
    "ch4": ("ch4_loading_298K_65bar_mol_kg_framework", "mol/kg-framework", 298.0, 6500000.0),
    "h2": ("h2_loading_77K_100bar_mol_kg_framework", "mol/kg-framework", 77.0, 10000000.0),
    "widom": ("co2_n2_henry_selectivity_298K", "dimensionless", 298.0, None),
}
STATUS = {
    "SUCCESS": ("SUCCESS", None),
    "EXCLUDED_SCIENTIFIC_ELIGIBILITY": ("NOT_AVAILABLE", "EXCLUDED_BY_TARGET_ELIGIBILITY_POLICY"),
    "MISSING_NO_VALIDATED_RESULT_AS_OF_CUTOFF": ("NOT_AVAILABLE", "NO_ACCEPTED_TARGET"),
    "SUCCESS_NULL": ("NOT_AVAILABLE", "UNDEFINED_SCIENTIFIC_RESULT"),
    "HISTORICAL_SCIENTIFIC_NULL": ("NOT_AVAILABLE", "UNDEFINED_SCIENTIFIC_RESULT"),
    "ERROR": ("ERROR", "CALCULATION_ERROR"),
}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def canonical(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def strict_json(data):
    def reject(value):
        raise ValueError(f"Invalid JSON constant: {value}")
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"Duplicate JSON key: {key}")
            result[key] = value
        return result
    return json.loads(data, object_pairs_hook=unique, parse_constant=reject)


def relative_path(name):
    path = PurePosixPath(name)
    if not name or path.is_absolute() or ".." in path.parts or str(path) != name or "\\" in name:
        raise ValueError(f"Unsafe release-relative path: {name!r}")
    return path


def read_bound(root, name, expected):
    path = root / relative_path(name)
    if path.is_symlink() or not path.is_file() or any(p.is_symlink() for p in path.parents):
        raise ValueError(f"Not a regular bound input: {path}")
    data = path.read_bytes()
    if digest(data) != expected:
        raise ValueError(f"Input checksum mismatch: {path}")
    return data


def ledger(root, expected):
    data = read_bound(root, "SHA256SUMS", expected)
    result = {}
    for line in data.decode().splitlines():
        hash_value, name = line.split("  ", 1)
        relative_path(name)
        if len(hash_value) != 64 or any(c not in "0123456789abcdef" for c in hash_value):
            raise ValueError("Malformed checksum")
        if name in result or name == "SHA256SUMS":
            raise ValueError("Duplicate or self-referential checksum entry")
        result[name] = hash_value
    return result


def csv_rows(data):
    reader = csv.DictReader(io.StringIO(data.decode(), newline=""))
    if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
        raise ValueError("Missing or duplicate CSV header")
    rows = list(reader)
    if any(None in row or None in row.values() for row in rows):
        raise ValueError("Malformed CSV row")
    return reader.fieldnames, rows


def csv_bytes(fields, rows):
    out = io.StringIO(newline="")
    writer = csv.DictWriter(out, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return out.getvalue().encode()


def number(text):
    if not isinstance(text, str):
        raise ValueError("Expected an explicitly declared CSV scalar")
    if text == "":
        return None
    result = float(text)
    if not math.isfinite(result):
        raise ValueError("Non-finite target or condition")
    return result


def public_target(row):
    """Project scientific fields only, never private paths, IDs or worker details."""
    endpoint = row["endpoint"]
    name, unit, temperature, pressure = ENDPOINTS[endpoint]
    if (row["target_name"] != name or row["unit"] != unit
            or number(row["temperature_K"]) != temperature
            or number(row["pressure_Pa"]) != pressure):
        raise ValueError("Target endpoint or conditions do not match the declared contract")
    status, reason = STATUS[row["assignment_status"]]
    value, error = number(row["value"]), number(row["reported_error"])
    available = status == "SUCCESS"
    if row["available"] not in ("true", "false") or (row["available"] == "true") != available:
        raise ValueError("Target availability differs from status")
    if available != (value is not None) or (error is not None and (not available or error < 0)):
        raise ValueError("Target status, value and reported uncertainty are inconsistent")
    method = {}
    for key in ("protocol_charge_method", "cif_preprocessing_action", "raspa_version"):
        method_value = row.get(key, "")
        if method_value and (len(method_value) > 120 or re.fullmatch(r"[A-Za-z0-9 ._()+-]+", method_value) is None
                      or re.search(r"[0-9a-fA-F]{64}", method_value)):
            raise ValueError("Method metadata contains an unsupported or private value")
        method[key] = method_value or None
    return {
        "available": available, "execution_status": status,
        "value": value, "reported_error": error, "unit": unit,
        "conditions": {"temperature_K": temperature, "pressure_Pa": pressure},
        "calculation_method": method,
        "diagnostic": None if available else {
            "code": reason, "category": "target_availability",
            "message": {
                "EXCLUDED_BY_TARGET_ELIGIBILITY_POLICY": "No accepted target under the recorded target-eligibility policy.",
                "NO_ACCEPTED_TARGET": "No validated target in this metadata revision.",
                "UNDEFINED_SCIENTIFIC_RESULT": "The completed calculation did not define a finite target.",
                "CALCULATION_ERROR": "The recorded target calculation was unsuccessful.",
            }[reason],
            "retry_action": "NONE_AUTOMATIC",
        },
    }


def target_schema():
    diagnostic = {
        "type": ["object", "null"], "additionalProperties": False,
        "required": ["code", "category", "message", "retry_action"],
        "properties": {key: {"type": "string"} for key in ("code", "category", "message", "retry_action")},
    }
    properties = {}
    for name, unit, temperature, pressure in ENDPOINTS.values():
        properties[name] = {
            "type": "object", "additionalProperties": False,
            "required": ["available", "execution_status", "value", "reported_error", "unit", "conditions", "diagnostic", "calculation_method"],
            "properties": {
                "available": {"type": "boolean"},
                "execution_status": {"enum": ["SUCCESS", "ERROR", "NOT_AVAILABLE"]},
                "value": {"type": ["number", "null"]},
                "reported_error": {"type": ["number", "null"], "minimum": 0},
                "unit": {"const": unit}, "diagnostic": diagnostic,
                "calculation_method": {"type": "object", "additionalProperties": False,
                    "required": ["protocol_charge_method", "cif_preprocessing_action", "raspa_version"],
                    "properties": {key: {"type": ["string", "null"], "maxLength": 120}
                        for key in ("protocol_charge_method", "cif_preprocessing_action", "raspa_version")}},
                "conditions": {"type": "object", "additionalProperties": False,
                    "required": ["temperature_K", "pressure_Pa"],
                    "properties": {"temperature_K": {"const": temperature}, "pressure_Pa": {"const": pressure}}},
            },
            "allOf": [{"if": {"properties": {"available": {"const": True}}},
                "then": {"properties": {"value": {"type": "number"}, "execution_status": {"const": "SUCCESS"}, "diagnostic": {"type": "null"}}},
                "else": {"properties": {"value": {"type": "null"}, "reported_error": {"type": "null"}, "execution_status": {"enum": ["ERROR", "NOT_AVAILABLE"]}, "diagnostic": {"type": "object"}}}}],
        }
    return {"type": "object", "additionalProperties": False, "required": list(properties), "properties": properties}


def build(source, target_root, output, source_ledger_hash, target_ledger_hash, *, revision_id=None):
    """Return an independent private copy with all original scientific fields retained."""
    from jsonschema import Draft202012Validator

    source, target_root, output = Path(source).resolve(), Path(target_root).resolve(), Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if any(parent.is_symlink() for parent in output.parents):
        raise ValueError("Output parents must not be symlinks")
    output = output.resolve()
    revision_id = output.name if revision_id is None else revision_id
    if not isinstance(revision_id, str) or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,79}", revision_id) is None:
        raise ValueError("Use a short explicit metadata revision identifier")
    if any(root == output or root in output.parents or output in root.parents for root in (source, target_root)):
        raise ValueError("Output and inputs must be separate directory trees")
    source_files, target_files = ledger(source, source_ledger_hash), ledger(target_root, target_ledger_hash)
    required_files = {"dataset_info.json", "README.md", "manifests/cif_manifest.csv",
        "metadata/metadata.csv", "metadata/metadata.jsonl", "metadata/structure_record_schema.json",
        "manifests/structure_json_manifest.csv", "manifests/validation.json", "manifests/metadata_revision_validation.json"}
    if not required_files.issubset(source_files):
        raise ValueError("Source release omits required metadata or validation files")
    info = strict_json(read_bound(source, "dataset_info.json", source_files["dataset_info.json"]))
    # This exporter never grants a missing evidence or redistribution permission.
    if "metadata_revision" not in info:
        raise ValueError("A validated metadata revision is required")
    _, cif_rows = csv_rows(read_bound(source, "manifests/cif_manifest.csv", source_files["manifests/cif_manifest.csv"]))
    cif_hashes = {row["structure_id"]: row["sha256"] for row in cif_rows}
    if len(cif_hashes) != len(cif_rows) or len(cif_rows) != info["structure_count"]:
        raise ValueError("CIF membership is not unique and complete")
    for sid, expected in cif_hashes.items():
        if source_files.get(f"cifs/{sid}.cif") != expected:
            raise ValueError("CIF manifest and file ledger disagree")
    _, rows = csv_rows(read_bound(target_root, "target_assignments.csv", target_files["target_assignments.csv"]))
    targets, seen = {sid: {} for sid in cif_hashes}, set()
    for row in rows:
        key = row["structure_id"], row["endpoint"]
        if key in seen:
            raise ValueError("Duplicate target endpoint")
        seen.add(key)
        if row["structure_id"] not in cif_hashes:
            continue  # The superset's complete target table may cover both releases.
        if row["cif_sha256"] != cif_hashes[row["structure_id"]]:
            raise ValueError("Target refers to different structure bytes")
        name = ENDPOINTS[row["endpoint"]][0]
        targets[row["structure_id"]][name] = public_target(row)
    expected_names = {item[0] for item in ENDPOINTS.values()}
    if any(set(value) != expected_names for value in targets.values()):
        raise ValueError("Every structure must have every target endpoint, including unavailable targets")
    original_schema = strict_json(read_bound(source, "metadata/structure_record_schema.json", source_files["metadata/structure_record_schema.json"]))
    schema = copy.deepcopy(original_schema)
    schema["$id"] = "urn:coremof:schema:structure-record:1.1"
    schema["properties"]["schema_version"] = {"const": "coremof-structure-record/1.1"}
    schema["properties"]["targets"] = target_schema()
    schema["required"].append("targets")
    Draft202012Validator.check_schema(schema)
    validator = Draft202012Validator(schema)
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".target-metadata-", dir=output.parent))
    produced, changed, counts = {}, [], Counter()

    def emit(name, data):
        path = staging / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        if path.read_bytes() != data:
            raise ValueError(f"Output readback differs: {name}")
        produced[name] = digest(data)

    try:
        json_manifest = []
        for name, expected in sorted(source_files.items()):
            data = read_bound(source, name, expected)
            if name.startswith("metadata/structures/") and name.endswith(".json"):
                record = strict_json(data)
                sid = record["structure_id"]
                if name != f"metadata/structures/{sid}.json" or sid not in targets or "targets" in record:
                    raise ValueError("Invalid structure record membership or existing targets")
                original = copy.deepcopy(record)
                record["schema_version"] = "coremof-structure-record/1.1"
                record["targets"] = targets[sid]
                validator.validate(record)
                restored = {k: v for k, v in record.items() if k != "targets"}
                restored["schema_version"] = original["schema_version"]
                if restored != original:
                    raise ValueError("Original scientific record changed")
                data = canonical(record)
                json_manifest.append({"structure_id": sid, "json_file": name, "size_bytes": len(data), "sha256": digest(data)})
            elif name == "metadata/metadata.csv":
                fields, metadata = csv_rows(data)
                if len(metadata) != len(cif_hashes) or {row["structure_id"] for row in metadata} != set(cif_hashes):
                    raise ValueError("Metadata CSV membership is duplicated or incomplete")
                fields += [field for key in sorted(expected_names) for field in (key, key + "__status")]
                if len(fields) != len(set(fields)):
                    raise ValueError("Target column collides with existing metadata")
                for row in metadata:
                    for key, target in targets[row["structure_id"]].items():
                        row[key] = "" if target["value"] is None else repr(target["value"])
                        row[key + "__status"] = target["execution_status"]
                data = csv_bytes(fields, metadata)
            elif name == "metadata/metadata.jsonl":
                enriched, metadata_ids = [], set()
                for line in data.splitlines():
                    row = strict_json(line)
                    if row["structure_id"] in metadata_ids or row["structure_id"] not in targets:
                        raise ValueError("Metadata JSONL membership is duplicated or unknown")
                    metadata_ids.add(row["structure_id"])
                    if "targets" in row:
                        raise ValueError("Target object already exists")
                    row["targets"] = targets[row["structure_id"]]
                    enriched.append(json.dumps(row, sort_keys=True, allow_nan=False, separators=(",", ":")))
                if metadata_ids != set(cif_hashes):
                    raise ValueError("Metadata JSONL membership is incomplete")
                data = ("\n".join(enriched) + "\n").encode()
            elif name == "metadata/structure_record_schema.json":
                data = canonical(schema)
            elif name == "metadata/analysis_metadata_schema.json":
                analysis = strict_json(data)
                key = "metadata/structure_annotations.jsonl"
                if key in analysis.get("formats", {}):
                    analysis["formats"][key] = analysis["formats"][key].replace(
                        "joins structure-record/1.0 without changing it",
                        "joins structure-record/1.1, with unchanged checker/group annotations")
                data = canonical(analysis)
            elif name in ("dataset_info.json", "README.md", "manifests/structure_json_manifest.csv", "manifests/validation.json", "manifests/metadata_revision_validation.json"):
                continue  # Regenerated below, never presented as a current old audit.
            if digest(data) != expected:
                changed.append(name)
            emit(name, data)
        if len(json_manifest) != len(cif_hashes) or {row["structure_id"] for row in json_manifest} != set(cif_hashes):
            raise ValueError("Per-structure JSON membership differs from CIF membership")
        emit("manifests/structure_json_manifest.csv", csv_bytes(["structure_id", "json_file", "size_bytes", "sha256"], json_manifest))
        for records in targets.values():
            for name, target in records.items():
                counts[(name, target["execution_status"])] += 1
        summary = {name: {status: counts[(name, status)] for status in ("SUCCESS", "NOT_AVAILABLE", "ERROR")} for name in sorted(expected_names)}
        emit("metadata/targets_schema.json", canonical(target_schema()))
        emit("metadata/targets.jsonl", ("\n".join(json.dumps({"structure_id": sid, "targets": value}, sort_keys=True, allow_nan=False) for sid, value in sorted(targets.items())) + "\n").encode())
        binding = {"schema_version": "coremof-target-metadata-inputs/1.0",
            "source_release_ledger_sha256": source_ledger_hash,
            "accepted_targets_ledger_sha256": target_ledger_hash,
            "accepted_assignments_sha256": target_files["target_assignments.csv"],
            "builder_sha256": BUILDER_SHA256,
            "source_record_schema_sha256": source_files["metadata/structure_record_schema.json"],
            "inherited_scientific_evidence_changed": False}
        emit("manifests/target_metadata_inputs.json", canonical(binding))
        validation = {
            "schema_version": "coremof-target-metadata-validation/1.0",
            "validation_status": "PASS_STAGED_NOT_PUBLISHED", "publication_authorized": False,
            "official_split": False, "dataset_version": info["dataset_version"],
            "structure_count": len(cif_hashes), "original_scientific_fields_unchanged": True,
            "source_files_read_and_verified": len(source_files), "per_structure_schema_validations": len(json_manifest),
            "targets": summary, "source_cifs_changed": False, "checker_or_grouping_changes": False,
            "frozen_benchmarks_changed": False, "publication_gate": "Inherited MOFid and redistribution gates remain unresolved.",
            "input_bindings": "manifests/target_metadata_inputs.json",
        }
        emit("manifests/validation.json", canonical(validation))
        emit("manifests/metadata_revision_validation.json", canonical(validation))
        info["metadata_revision"] = {**info["metadata_revision"],
            "revision_id": revision_id, "core_structure_record_schema_unchanged": False,
            "targets_modified": True, "frozen_splits_modified": False,
            "validation_manifest": "manifests/metadata_revision_validation.json"}
        info["metadata_revision"].pop("targets_or_frozen_splits_modified", None)
        info["target_metadata"] = {"schema": "metadata/targets_schema.json", "jsonl": "metadata/targets.jsonl", "counts": summary,
            "target_magnitudes_used_for_grouping": False, "targets_are_simulation_results": True,
            "henry_selectivity": "CO2/N2 infinite-dilution Henry selectivity at 298 K; dimensionless. No finite-pressure uptake or individual gas Henry constants are inferred."}
        info["per_structure_json"].update(schema_version="coremof-structure-record/1.1", record_count=len(json_manifest))
        # Preserve the source's exact, loader-validated release status. This
        # target-only export grants no new authority. Candidate/publication and
        # split declarations belong in the validation receipt, not in new
        # dataset_info authority fields that the loader deliberately rejects.
        # Refresh tabular descriptors only for changed existing tables.
        for name in ("metadata/metadata.csv", "metadata/metadata.jsonl"):
            if name in info.get("tabular_files", {}):
                entry = info["tabular_files"][name]
                entry["size_bytes"] = (staging / name).stat().st_size
                if name.endswith(".csv"):
                    with (staging / name).open() as stream:
                        entry["columns"] = next(csv.reader(stream))
        info.setdefault("tabular_files", {})["metadata/targets.jsonl"] = {
            "row_count": len(cif_hashes), "size_bytes": (staging / "metadata/targets.jsonl").stat().st_size}
        emit("dataset_info.json", canonical(info))
        target_readme = ("# CoRE-MOF-COD target-enriched metadata candidate\n\n"
            "This is a private, staged copy, not an authorized public release.\n"
            "CIFs, checker evidence, feature values and related-structure groups are unchanged.\n"
            "Accepted targets are included in metadata.csv, metadata.jsonl and every structure JSON.\n"
            "Join all tables by structure_id. Missing targets remain null with an explicit reason.\n"
            "CH4 and H2 uptake use mol/kg-framework. CO2/N2 Henry selectivity is dimensionless.\n"
            "Each target records its temperature and, for uptake, pressure.\n"
            "The per-structure record schema is now 1.1, adding targets to the unchanged 1.0 scientific fields.\n"
            "Do not substitute this metadata into a frozen benchmark.\n"
            "Final MOFid admissibility and asset-level sharing permissions are still required before public release.\n"
            "In particular, CSD-derived material is licence-gated and SI material requires rights review.\n"
            "Input checksums and the builder identity are in manifests/target_metadata_inputs.json.\n"
            "calculation_method retains recorded protocol fields; null means unrecorded, not a guessed method.\n"
            "The protocol charge method is the electrostatics setting, not evidence of the charge-prediction method.\n")
        target_readme += "Checker/group annotation revision labels retain their recorded value because only targets were added.\n"
        source_readme = read_bound(source, "README.md", source_files["README.md"]).decode()
        marker = "## Checker information"
        if marker in source_readme:
            retained = source_readme[source_readme.index(marker):]
            retained = retained.replace("The established `coremof-structure-record/1.0` schema is unchanged.",
                "The `coremof-structure-record/1.1` schema adds targets while retaining the existing scientific fields.")
            retained = re.sub(r"Revision identifier: `[^`]+`\.", f"Revision identifier: `{revision_id}`.", retained)
            target_readme += "\n" + retained
        emit("README.md", target_readme.encode())
        emit("SHA256SUMS", "".join(f"{hash_value}  {name}\n" for name, hash_value in sorted(produced.items())).encode())
        # Prove all emitted bytes once more independently of their original write.
        for name, expected in produced.items():
            read_bound(staging, name, expected)
        loaded = CoREMOFDataset.from_release(staging, verify_cif_files=True)
        if len(loaded) != len(cif_hashes) or loaded.dataset_version != info["dataset_version"]:
            raise ValueError("The public loader did not preserve release identity and membership")
        publish_directory(staging, output, overwrite=False)
        return {**validation, "source_release_ledger_sha256": source_ledger_hash,
            "accepted_targets_ledger_sha256": target_ledger_hash, "output_ledger_sha256": produced["SHA256SUMS"],
            "changed_structure_records": sum(name.startswith("metadata/structures/") for name in changed),
            "changed_other_payloads": [name for name in changed if not name.startswith("metadata/structures/")],
            "builder_sha256": BUILDER_SHA256}
    except BaseException as exc:
        # Keep failed staging for diagnosis; never expose it as a valid release.
        setattr(exc, "coremof_preserved_staging_directory", str(staging))
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--targets", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--source-ledger-sha256", required=True)
    parser.add_argument("--target-ledger-sha256", required=True)
    parser.add_argument("--revision-id", required=True)
    args = parser.parse_args()
    receipt = build(args.source, args.targets, args.output, args.source_ledger_sha256, args.target_ledger_sha256,
        revision_id=args.revision_id)
    print(json.dumps(receipt, sort_keys=True, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
