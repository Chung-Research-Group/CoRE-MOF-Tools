import copy
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


EXAMPLE = Path(__file__).resolve().parents[1] / "examples/extend_collected_targets.py"
SPEC = importlib.util.spec_from_file_location("extend_collected_targets_example", EXAMPLE)
workflow = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(workflow)


class CollectedTargetUpdateTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.sid = "2020[Cu][nan]3[FSR]1"
        self.cif_hash = hashlib.sha256(b"synthetic CIF").hexdigest()
        self.task = {"structure_id": self.sid, "endpoint": "ch4", "task_id": "test-ch4",
                     "status": "MISSING", "cif_sha256": self.cif_hash}
        bundle = self.root / "result bundle"
        bundle.mkdir()
        target = {"value": 0.0, "error": 0.0, "unit": "mol/kg-framework"}
        self.result = {"status": "SUCCESS", "target": target, "failure": None,
                       "task_provenance": dict(self.task),
                       **{k: self.task[k] for k in ("structure_id", "endpoint", "task_id")}}
        result_bytes = workflow.canonical(self.result)
        receipt = {"status": "SUCCESS", "failure": None, "cif_sha256": self.cif_hash,
                   "task_payload_sha256": workflow.sha(workflow.canonical(self.task)),
                   "result_sha256": workflow.sha(result_bytes), "protocol": {"test_fixture": True},
                   "process": {"returncode": 0, "timed_out": False},
                   **{k: self.task[k] for k in ("structure_id", "endpoint", "task_id")}}
        receipt["receipt_payload_sha256"] = workflow.sha(workflow.canonical(receipt))
        receipt_bytes = workflow.canonical(receipt)
        (bundle / "result.json").write_bytes(result_bytes)
        (bundle / "receipt.json").write_bytes(receipt_bytes)
        self.record = {"classification": "TERMINAL_SUCCESS", "target": target, "errors": [],
                       "validated_result": True, "validated_receipt": True, "bundle_path": str(bundle),
                       "result_sha256": workflow.sha(result_bytes), "receipt_file_sha256": workflow.sha(receipt_bytes),
                       "receipt_payload_sha256": receipt["receipt_payload_sha256"],
                       "task_payload_sha256": receipt["task_payload_sha256"],
                       **{k: self.task[k] for k in ("structure_id", "endpoint", "task_id")}}
        self.wide = {self.sid: {"structure_id": self.sid, "cif_sha256": self.cif_hash}}
        self.long = {}
        fields = ("assignment_status", "available", "value", "reported_error", "source_class",
                  "source_batch", "task_id", "task_payload_sha256", "result_sha256",
                  "receipt_file_sha256", "receipt_payload_sha256", "source_diagnostic",
                  "validation_scope", "snapshot_id", "raspa_version")
        for endpoint, column in workflow.TARGETS.items():
            self.wide[self.sid].update({column: "", column + "__eligible": "true",
                                       column + "__status": "MISSING", column + "__reported_error": "",
                                       column + "__source_class": "MISSING"})
            self.long[(self.sid, endpoint)] = {k: "" for k in fields}
            self.long[(self.sid, endpoint)].update(structure_id=self.sid, endpoint=endpoint)

    def test_complete_saved_binding_and_physical_zero(self):
        self.assertEqual(workflow.verify_saved_observation(self.task, self.record), self.result)
        outcome = workflow.apply_observation(self.wide, self.long, self.task, self.record, 1)
        self.assertEqual(outcome, "TERMINAL_SUCCESS")
        self.assertEqual(self.wide[self.sid][workflow.TARGETS["ch4"]], "0.0")

    def test_recorded_method_is_not_inferred_and_receipt_hash_is_checked(self):
        self.assertEqual(workflow.recorded_method(self.record, self.result), {
            "protocol_charge_method": "", "cif_preprocessing_action": "", "raspa_version": ""})
        path = Path(self.record["bundle_path"]) / "receipt.json"
        receipt = json.loads(path.read_text())
        receipt["protocol"] = {
            "simulation": {"Systems": [{"ChargeMethod": "Ewald"}]},
            "cif_preprocessing": {"action": "UNCHANGED_COMPLETE_CHARGES"},
        }
        payload = workflow.canonical(receipt)
        path.write_bytes(payload)
        with self.assertRaises(ValueError):
            workflow.recorded_method(self.record, self.result)
        record = {**self.record, "receipt_file_sha256": workflow.sha(payload)}
        result = {**self.result, "output_observations": {"raspa_version": "3.1.0"}}
        self.assertEqual(workflow.recorded_method(record, result), {
            "protocol_charge_method": "Ewald", "cif_preprocessing_action": "UNCHANGED_COMPLETE_CHARGES",
            "raspa_version": "3.1.0"})
        receipt["protocol"]["simulation"]["Systems"].append({"ChargeMethod": "None"})
        payload = workflow.canonical(receipt)
        path.write_bytes(payload)
        record["receipt_file_sha256"] = workflow.sha(payload)
        with self.assertRaisesRegex(ValueError, "ambiguous systems"):
            workflow.recorded_method(record, result)

    def test_incomplete_or_altered_receipts_are_rejected(self):
        for field, value in (("validated_result", False), ("validated_receipt", False),
                             ("task_payload_sha256", "0" * 64), ("result_sha256", "0" * 64),
                             ("receipt_file_sha256", "0" * 64), ("target", {"value": 99})):
            record = copy.deepcopy(self.record)
            record[field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                workflow.verify_saved_observation(self.task, record)

    def test_nonfinite_boolean_units_and_uncertainty_fail(self):
        for field, value in (("value", True), ("value", float("nan")), ("value", float("inf")),
                             ("unit", "cm3/g"), ("error", -1), ("error", float("nan"))):
            record = copy.deepcopy(self.record)
            record["target"][field] = value
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                workflow.apply_observation(self.wide, self.long, self.task, record, 1)

    def test_accepted_value_and_accepted_null_cannot_be_replaced(self):
        column = workflow.TARGETS["ch4"]
        self.wide[self.sid][column] = "5.1250"
        with self.assertRaisesRegex(ValueError, "Conflicting accepted target"):
            workflow.apply_observation(self.wide, self.long, self.task, self.record, 1)
        self.assertEqual(self.wide[self.sid][column], "5.1250")
        record = copy.deepcopy(self.record)
        record["target"]["value"] = 5.125
        self.assertEqual(workflow.apply_observation(self.wide, self.long, self.task, record, 1),
                         "ALREADY_ACCEPTED_IDENTICAL")
        self.assertEqual(self.wide[self.sid][column], "5.1250")
        self.wide[self.sid][column] = ""
        self.wide[self.sid][column + "__status"] = "HISTORICAL_SCIENTIFIC_NULL"
        with self.assertRaisesRegex(ValueError, "accepted scientific null"):
            workflow.apply_observation(self.wide, self.long, self.task, self.record, 1)

    def test_excluded_changed_cif_and_nonmissing_task_fail(self):
        self.wide[self.sid][workflow.TARGETS["ch4"] + "__eligible"] = "false"
        with self.assertRaisesRegex(ValueError, "excluded"):
            workflow.apply_observation(self.wide, self.long, self.task, self.record, 1)
        for field, value in (("cif_sha256", "0" * 64), ("status", "EXISTING")):
            task = {**self.task, field: value}
            with self.assertRaises(ValueError):
                workflow.apply_observation(self.wide, self.long, task, self.record, 1)

    def test_end_to_end_output_and_hashes_without_changing_baseline(self):
        baseline = self.root / "baseline"
        (baseline / "targets").mkdir(parents=True)

        def write_csv(path, rows):
            with path.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)

        write_csv(baseline / "targets/targets.csv", list(self.wide.values()))
        write_csv(baseline / "targets/target_assignments.csv", list(self.long.values()))
        cif_manifest = self.root / "cif_manifest.csv"
        write_csv(cif_manifest, [{"structure_id": self.sid, "sha256": self.cif_hash}])
        (baseline / "targets/merge_receipt.json").write_bytes(workflow.canonical({"inputs": {
            "release_manifest": {"path": str(cif_manifest), "sha256": workflow.sha(cif_manifest.read_bytes())}}}))
        (baseline / "targets/coverage.json").write_text("{}\n")
        (baseline / "targets/targets.json").write_bytes(workflow.canonical({"sources": [{"path": "targets.csv"}]}))
        files = {p: p.read_bytes() for p in (baseline / "targets").iterdir()}
        (baseline / "SHA256SUMS").write_text("".join(
            f"{workflow.sha(data)}  {path.relative_to(baseline)}\n" for path, data in files.items()))
        collection = self.root / "collection.json"
        collection.write_bytes(workflow.canonical({"records": [self.record], "counts": {"TERMINAL_SUCCESS": 1},
                                                    "time": "2026-09-22T00:00:00+00:00"}))
        plan = self.root / "plan.json"
        plan.write_bytes(workflow.canonical({"tasks": [self.task]}))
        audit = self.root / "audit.json"
        audit.write_bytes(workflow.canonical({"batches": [
            {"batch": 1, "collection": {"path": str(collection), "sha256": workflow.sha(collection.read_bytes())},
             "plan": {"path": str(plan), "sha256": workflow.sha(plan.read_bytes())}},
            {"batch": 2}]}))
        output = self.root / "new targets"
        summary = workflow.extend(baseline, audit, output)
        self.assertEqual(summary["finite_targets"]["ch4"], 1)
        self.assertFalse(summary["benchmark_assignments_changed"])
        self.assertIn({"batch": 2, "reason": "NO_COLLECTOR_REPORT"}, summary["unavailable"])
        for path, data in files.items():
            self.assertEqual(path.read_bytes(), data)
        for line in (output / "SHA256SUMS").read_text().splitlines():
            expected, filename = line.split(maxsplit=1)
            self.assertEqual(workflow.sha((output / filename).read_bytes()), expected)
        # A saved extension is itself an accepted baseline for the next batch.
        # Replaying the same collector entry cannot count or alter it twice.
        next_output = self.root / "next targets"
        next_summary = workflow.extend(output, audit, next_output)
        self.assertEqual(next_summary["new_finite_targets"], {})
        self.assertEqual(next_summary["previous_finite_values_unchanged"], 1)
        self.assertEqual(next_summary["batch_merge_outcomes"]["1"], {"ALREADY_ACCEPTED_IDENTICAL": 1})
        self.assertEqual((next_output / "targets.csv").read_bytes(), (output / "targets.csv").read_bytes())
        with self.assertRaises(FileExistsError):
            workflow.extend(baseline, audit, output)

    def test_strict_json_rejects_ambiguous_values(self):
        for payload in ('{"a":1,"a":2}', '{"a":NaN}', '{"a":Infinity}'):
            with self.assertRaises(ValueError):
                workflow.strict_json(payload)


if __name__ == "__main__":
    unittest.main()
