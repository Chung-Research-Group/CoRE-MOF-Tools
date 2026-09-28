"""Extend an accepted target snapshot using saved, validated collector results.

No scheduler command or scientific program is called. Existing values are
immutable. The output is a private candidate, not a redistribution decision or
a new benchmark assignment. Full artifact validation is reused from the frozen
collector reports; the exact result/receipt bytes are verified again here.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import shutil
import tempfile

from CoREMOF._transactions import publish_directory


TARGETS = {
    "ch4": "ch4_loading_298K_65bar_mol_kg_framework",
    "h2": "h2_loading_77K_100bar_mol_kg_framework",
    "widom": "co2_n2_henry_selectivity_298K",
}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def canonical(data):
    return (json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()


def strict_json(data):
    def reject(value):
        raise ValueError(f"Non-finite JSON constant: {value}")
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"Duplicate JSON key: {key}")
            result[key] = value
        return result
    return json.loads(data, parse_constant=reject, object_pairs_hook=unique)


def capture(path, expected=None):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"Expected a regular input file: {path}")
    data = path.read_bytes()
    digest = sha(data)
    if expected is not None and digest != expected:
        raise ValueError(f"Input checksum mismatch: {path}")
    return data, {"path": str(path.resolve()), "sha256": digest, "size_bytes": len(data)}


def csv_rows(data):
    reader = csv.DictReader(io.StringIO(data.decode("utf-8"), newline=""))
    if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
        raise ValueError("Missing or duplicate CSV header")
    rows = list(reader)
    if any(None in row or None in row.values() for row in rows):
        raise ValueError("Malformed CSV row")
    return rows


def finite(value):
    return type(value) in (int, float) and math.isfinite(value)


def verify_saved_observation(task, record):
    """Verify exact task/result/receipt binding without repeating science."""
    for key in ("task_id", "structure_id", "endpoint"):
        if task[key] != record[key]:
            raise ValueError(f"Collector {key} differs from planned task")
    task_digest = sha(canonical(task))
    if record["task_payload_sha256"] != task_digest:
        raise ValueError("Task payload checksum differs from collector")
    if record["classification"] not in ("TERMINAL_SUCCESS", "TERMINAL_NULL"):
        return None
    if (record["errors"] or record["validated_result"] is not True
            or record["validated_receipt"] is not True):
        raise ValueError("Collector did not validate this observation")
    bundle = Path(record["bundle_path"])
    if bundle.is_symlink():
        raise ValueError("Worker bundle must not be a symlink")
    result_bytes, _ = capture(bundle / "result.json", record["result_sha256"])
    receipt_bytes, _ = capture(bundle / "receipt.json", record["receipt_file_sha256"])
    result, receipt = strict_json(result_bytes), strict_json(receipt_bytes)
    receipt_payload = dict(receipt)
    stored = receipt_payload.pop("receipt_payload_sha256")
    if stored != record["receipt_payload_sha256"] or sha(canonical(receipt_payload)) != stored:
        raise ValueError("Worker receipt payload checksum mismatch")
    if receipt["result_sha256"] != record["result_sha256"]:
        raise ValueError("Worker receipt refers to a different result")
    if receipt["task_payload_sha256"] != task_digest or receipt["cif_sha256"] != task["cif_sha256"]:
        raise ValueError("Worker receipt task/CIF mismatch")
    for payload in (result, receipt):
        for key in ("task_id", "structure_id", "endpoint"):
            if payload[key] != task[key]:
                raise ValueError(f"Worker {key} differs from planned task")
    expected_status = "SUCCESS" if record["classification"] == "TERMINAL_SUCCESS" else "SUCCESS_NULL"
    if result["status"] != expected_status or receipt["status"] != expected_status:
        raise ValueError("Worker status differs from collector")
    if result.get("failure") is not None or receipt.get("failure") is not None:
        raise ValueError("Successful worker observation contains a failure")
    if result["task_provenance"] != {k: v for k, v in task.items() if k != "cif_path"}:
        raise ValueError("Worker result does not match the exact task")
    if result["target"] != record["target"]:
        raise ValueError("Collector target differs from canonical worker target")
    process = receipt["process"]
    if process["returncode"] != 0 or process["timed_out"] is not False:
        raise ValueError("Successful observation has unsuccessful process evidence")
    if not isinstance(receipt.get("protocol"), dict):
        raise ValueError("Missing calculation protocol")
    return result


def apply_observation(wide, long, task, record, batch):
    """Fill a missing eligible cell, never replace an accepted value or null."""
    sid, endpoint = task["structure_id"], task["endpoint"]
    column = TARGETS[endpoint]
    if sid not in wide or wide[sid]["cif_sha256"] != task["cif_sha256"]:
        raise ValueError("Unknown structure or changed CIF")
    if record["classification"] not in ("TERMINAL_SUCCESS", "TERMINAL_NULL"):
        return record["classification"]
    if task["status"] != "MISSING":
        raise ValueError("New observation is not a source-manifest MISSING key")
    target = record["target"]
    value = target["value"]
    if record["classification"] == "TERMINAL_SUCCESS" and not finite(value):
        raise ValueError("Successful target must be finite and numeric")
    if record["classification"] == "TERMINAL_NULL" and (endpoint != "widom" or value is not None):
        raise ValueError("Invalid scientific null")
    error = target.get("error")
    if error is not None and (not finite(error) or error < 0):
        raise ValueError("Target uncertainty must be finite and nonnegative")
    expected_unit = "dimensionless" if endpoint == "widom" else "mol/kg-framework"
    if target["unit"] != expected_unit:
        raise ValueError("Incompatible target units")
    row = wide[sid]
    if row[column] != "":
        if value is None or float(row[column]) != value:
            raise ValueError(f"Conflicting accepted target: {sid} {endpoint}")
        return "ALREADY_ACCEPTED_IDENTICAL"
    if row[column + "__eligible"] != "true":
        raise ValueError("Cannot fill an excluded structure")
    if row[column + "__status"] in ("HISTORICAL_SCIENTIFIC_NULL", "SUCCESS_NULL"):
        if value is None:
            return "ALREADY_ACCEPTED_NULL"
        raise ValueError("Cannot replace an accepted scientific null")
    row[column] = "" if value is None else repr(value)
    row[column + "__reported_error"] = "" if error is None else repr(error)
    row[column + "__status"] = "SUCCESS" if value is not None else "SUCCESS_NULL"
    row[column + "__source_class"] = "COLLECTOR_VALIDATED_NUMBERED_BATCH"
    long[(sid, endpoint)].update(
        assignment_status=row[column + "__status"], available=str(value is not None).lower(),
        value=row[column], reported_error=row[column + "__reported_error"],
        source_class="COLLECTOR_VALIDATED_NUMBERED_BATCH", source_batch=str(batch),
        task_id=task["task_id"], task_payload_sha256=record["task_payload_sha256"],
        result_sha256=record["result_sha256"], receipt_file_sha256=record["receipt_file_sha256"],
        receipt_payload_sha256=record["receipt_payload_sha256"],
        source_diagnostic=target.get("diagnostic") or target.get("null_reason") or "",
        validation_scope="frozen collector artifact validation and current exact result/receipt hash verification",
    )
    return record["classification"]


def recorded_method(record, result):
    """Project only explicitly recorded, hash-verified calculation settings."""
    data, _ = capture(Path(record["bundle_path"]) / "receipt.json", record["receipt_file_sha256"])
    protocol = strict_json(data)["protocol"]
    systems = protocol.get("simulation", {}).get("Systems", [])
    if not isinstance(systems, list) or len(systems) > 1:
        raise ValueError("Cannot select a calculation method from ambiguous systems")
    method = {
        "protocol_charge_method": systems[0].get("ChargeMethod") if systems else None,
        "cif_preprocessing_action": protocol.get("cif_preprocessing", {}).get("action"),
        "raspa_version": result.get("output_observations", {}).get("raspa_version"),
    }
    for value in method.values():
        if value is not None and (not isinstance(value, str) or not value):
            raise ValueError("Recorded calculation settings must be explicit nonempty strings")
    return {key: value if value is not None else "" for key, value in method.items()}


def extend(baseline, input_audit, output):
    baseline, input_audit, output = Path(baseline), Path(input_audit), Path(output)
    if output.exists():
        raise FileExistsError(f"Output already exists: {output}")
    ledger_bytes, ledger_binding = capture(baseline / "SHA256SUMS")
    ledger = {}
    for line in ledger_bytes.decode().splitlines():
        digest, name = line.split(maxsplit=1)
        name = name.lstrip("*")
        if name in ledger or Path(name).is_absolute() or ".." in Path(name).parts:
            raise ValueError("Invalid or duplicate checksum-ledger member")
        ledger[name] = digest
    inputs = {"baseline_ledger": ledger_binding, "baseline": {}}
    captured = {}
    required = ("targets.csv", "target_assignments.csv", "targets.json", "merge_receipt.json", "coverage.json")
    prefixes = [prefix for prefix in ("targets/", "") if all(prefix + name in ledger for name in required)]
    if len(prefixes) != 1:
        raise ValueError("Baseline must contain exactly one complete accepted-target layout")
    for filename in required:
        name = prefixes[0] + filename
        data, record = capture(baseline / name, ledger[name])
        inputs["baseline"][name] = record
        captured["targets/" + filename] = data
    previous_receipt = strict_json(captured["targets/merge_receipt.json"])
    previous_inputs = previous_receipt["inputs"]
    cif_binding = previous_inputs.get("release_manifest", previous_inputs.get("cif_manifest"))
    if not isinstance(cif_binding, dict):
        raise ValueError("Accepted baseline does not bind its original release CIF manifest")
    cif_bytes, inputs["cif_manifest"] = capture(cif_binding["path"], cif_binding["sha256"])
    cif_rows = csv_rows(cif_bytes)
    cif_ids = {r["structure_id"] for r in cif_rows}
    wide_rows = csv_rows(captured["targets/targets.csv"])
    wide = {r["structure_id"]: r for r in wide_rows}
    if len(wide) != len(wide_rows) or len(cif_ids) != len(cif_rows) or set(wide) != cif_ids:
        raise ValueError("Master target membership differs from the release")
    for row in cif_rows:
        if wide[row["structure_id"]]["cif_sha256"] != row["sha256"]:
            raise ValueError("Baseline targets refer to different CIF bytes")
    long_rows = csv_rows(captured["targets/target_assignments.csv"])
    for row in long_rows:
        for key in ("protocol_charge_method", "cif_preprocessing_action", "raspa_version"):
            row.setdefault(key, "")
    long = {(r["structure_id"], r["endpoint"]): r for r in long_rows}
    if len(long) != len(long_rows) or set(long) != {(sid, ep) for sid in wide for ep in TARGETS}:
        raise ValueError("Invalid endpoint accounting")
    previous = {(sid, col): row[col] for sid, row in wide.items() for col in TARGETS.values() if row[col] != ""}
    if any(not math.isfinite(float(value)) for value in previous.values()):
        raise ValueError("Non-finite value in accepted baseline")
    audit_bytes, inputs["input_audit"] = capture(input_audit)
    audit = strict_json(audit_bytes)
    seen = set()
    counts, additions, unavailable = {}, Counter(), []
    observations, input_batches = [], []
    for batch in audit["batches"]:
        if "collection" not in batch:
            unavailable.append({"batch": batch["batch"], "reason": "NO_COLLECTOR_REPORT"})
            continue
        collection_bytes, collection_binding = capture(batch["collection"]["path"], batch["collection"]["sha256"])
        plan_bytes, plan_binding = capture(batch["plan"]["path"], batch["plan"]["sha256"])
        collection, plan = strict_json(collection_bytes), strict_json(plan_bytes)
        tasks = {r["task_id"]: r for r in plan["tasks"]}
        records = {r["task_id"]: r for r in collection["records"]}
        if (len(tasks) != len(plan["tasks"]) or len(records) != len(collection["records"])
                or set(tasks) != set(records)):
            raise ValueError("Collection/plan membership mismatch")
        if Counter(r["classification"] for r in records.values()) != Counter(collection["counts"]):
            raise ValueError("Collection count mismatch")
        outcomes = Counter()
        for task_id in sorted(tasks):
            task, record = tasks[task_id], records[task_id]
            key = task["structure_id"], task["endpoint"]
            if key in seen:
                raise ValueError("Repeated endpoint task in the update")
            seen.add(key)
            result = verify_saved_observation(task, record)
            outcome = apply_observation(wide, long, task, record, batch["batch"])
            outcomes[outcome] += 1
            if outcome == "TERMINAL_SUCCESS":
                additions[task["endpoint"]] += 1
            if outcome in ("TERMINAL_SUCCESS", "TERMINAL_NULL"):
                long[key].update(recorded_method(record, result))
            if outcome not in ("TERMINAL_SUCCESS", "ALREADY_ACCEPTED_IDENTICAL"):
                unavailable.append({"batch": batch["batch"], "structure_id": key[0], "endpoint": key[1], "reason": outcome})
            observations.append({"batch": batch["batch"], "task": task, "audit": record, "merge_outcome": outcome})
        counts[str(batch["batch"])] = dict(outcomes)
        input_batches.append({"batch": batch["batch"], "plan": plan_binding, "collection": collection_binding,
                              "time": collection["time"]})
    if any(wide[sid][col] != value for (sid, col), value in previous.items()):
        raise ValueError("Accepted target value changed")
    for row in long.values():
        row["snapshot_id"] = output.name
    coverage = {ep: sum(row[col] != "" for row in wide.values()) for ep, col in TARGETS.items()}
    summary = {"release_structures": len(wide), "finite_targets": coverage,
               "all_three_targets": sum(all(row[col] != "" for col in TARGETS.values()) for row in wide.values()),
               "new_finite_targets": dict(additions), "previous_finite_values_unchanged": len(previous),
               "batch_merge_outcomes": counts, "unavailable": unavailable,
               "data_cutoff_utc": max((x["time"] for x in input_batches), default=None),
               "publication_authorized": False, "benchmark_assignments_changed": False,
               "science_launched": False}
    inputs["batches"] = input_batches
    inputs["implementation"] = capture(Path(__file__))[1]
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".target-update-", dir=output.parent) as temporary:
        stage = Path(temporary) / output.name
        stage.mkdir()
        for name, rows in (("targets.csv", [wide[k] for k in sorted(wide)]),
                           ("target_assignments.csv", [long[k] for k in sorted(long)])):
            with (stage / name).open("x", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
        config = strict_json(captured["targets/targets.json"])
        config["sources"][0]["path"] = "targets.csv"
        config["sources"][0]["name"] = output.name
        for name, value in (("targets.json", config), ("coverage.json", summary),
                            ("merge_receipt.json", {"inputs": inputs, "summary": summary, "policy": "immutable fill-only accepted union"})):
            (stage / name).write_bytes(canonical(value))
        with (stage / "new_observations.jsonl").open("x", encoding="utf-8") as stream:
            for item in observations:
                stream.write(json.dumps(item, sort_keys=True, allow_nan=False) + "\n")
        shutil.copy2(__file__, stage / "extend_collected_targets.py")
        (stage / "SHA256SUMS").write_text("".join(f"{sha(path.read_bytes())}  {path.name}\n" for path in sorted(stage.iterdir())))
        if output.exists():
            raise FileExistsError(output)
        publish_directory(stage, output, overwrite=False)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True, help="Frozen snapshot with SHA256SUMS and targets/, or a previous flat output of this example")
    parser.add_argument("--input-audit", type=Path, required=True, help="Hash-bound batch collection/plan inventory")
    parser.add_argument("--output", type=Path, required=True, help="New private target metadata directory")
    options = parser.parse_args()
    print(json.dumps(extend(options.baseline, options.input_audit, options.output), indent=2, allow_nan=False))
