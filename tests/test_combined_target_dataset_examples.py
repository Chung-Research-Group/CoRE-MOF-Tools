import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
BUILD_SCRIPT = REPOSITORY_ROOT / "examples" / "build_combined_target_dataset.py"
AUDIT_SCRIPT = REPOSITORY_ROOT / "examples" / "audit_combined_target_dataset.py"

ENDPOINTS = {
    "CH4": {
        "short": "ch4",
        "target": "ch4_loading_298K_65bar_mol_kg_framework",
        "unit": "mol/kg-framework",
        "temperature": "298.0",
        "pressure": "6500000.0",
    },
    "H2": {
        "short": "h2",
        "target": "h2_loading_77K_100bar_mol_kg_framework",
        "unit": "mol/kg-framework",
        "temperature": "77.0",
        "pressure": "10000000.0",
    },
    "WIDOM_CO2_N2": {
        "short": "widom",
        "target": "co2_n2_widom_rosenbluth_weight_ratio_298K_1bar",
        "unit": "dimensionless",
        "temperature": "298.0",
        "pressure": "100000.0",
    },
}

CURRENT_FIELDS = (
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


def _json_bytes(value):
    return (
        json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False)
        + "\n"
    ).encode("utf-8")


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _csv_bytes(fieldnames, rows):
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue().encode("utf-8")


def _source_row(structure_id, cif_sha256, endpoint, phase, status):
    row = {
        "manifest_schema_version": "coremof-gcmc-manifest/1.0",
        "structure_id": structure_id,
        "cif_sha256": cif_sha256,
        "endpoint": endpoint,
        "release_phase": phase,
        "release_version": "v26.0.1" if phase == "v26.0.1_base" else "v26.0.2",
        "status": status,
        "alias_conflict": False,
        "exclusion_reasons": [],
        "explicit_null": status != "EXISTING",
        "has_existing_result": status == "EXISTING",
    }
    if status == "EXISTING":
        is_scientific_null = (
            structure_id == "BASE_EXISTING" and endpoint == "WIDOM_CO2_N2"
        )
        row.update(
            {
                "explicit_null": is_scientific_null,
                "result_value": (
                    None
                    if is_scientific_null
                    else 10.0 + list(ENDPOINTS).index(endpoint)
                ),
                "result_error": None if is_scientific_null else 0.1,
                "result_diagnostic": "ZERO_DENOMINATOR" if is_scientific_null else "",
                "result_semantics": (
                    "absolute_loading_mol_kg_framework"
                    if endpoint in {"CH4", "H2"}
                    else "widom_rosenbluth_weight_ratio_co2_over_n2"
                ),
                "selected_source_file": "synthetic/historical.csv",
                "selected_source_sha256": "e" * 64,
                "selected_legacy_id": structure_id + "-legacy",
            }
        )
    elif status == "EXCLUDED":
        # An excluded source can have a historical candidate, but it must still
        # produce a native null in the attachment dataset.
        row.update(
            {
                "exclusion_reasons": ["MOFCHECKER_LONE_MOLECULE"],
                "has_existing_result": True,
                "result_value": 99.0,
                "selected_source_file": "synthetic/excluded-candidate.csv",
                "selected_source_sha256": "d" * 64,
                "selected_legacy_id": structure_id + "-excluded-candidate",
            }
        )
    return row


def _current_row(structure_id, cif_sha256, endpoint, phase, status, value=""):
    contract = ENDPOINTS[endpoint]
    return {
        "dataset_version": "v26.0.1" if phase == "v26.0.1_base" else "v26.0.2",
        "release_phase": phase,
        "structure_id": structure_id,
        "cif_sha256": cif_sha256,
        "endpoint": contract["short"],
        "target_name": contract["target"],
        "target_semantics": (
            "absolute_loading_mol_kg_framework"
            if endpoint in {"CH4", "H2"}
            else "widom_rosenbluth_weight_ratio_co2_over_n2"
        ),
        "status": status,
        "available": "true" if status == "SUCCESS" else "false",
        "value": str(value) if status == "SUCCESS" else "",
        "error": "0.02" if status == "SUCCESS" else "",
        "unit": contract["unit"],
        "temperature_kelvin": contract["temperature"],
        "pressure_pascal": contract["pressure"],
        "protocol_charge_method": "PACMAN",
        "cif_preprocessing_action": "NONE",
        "raspa_version": "3.0-test",
        "task_id": structure_id + "-" + contract["short"],
        "task_payload_sha256": "a" * 64,
        "result_sha256": "b" * 64,
        "receipt_file_sha256": "c" * 64,
        "receipt_payload_sha256": "f" * 64,
        "source_batch": "synthetic-test",
        "job_id": "test-job",
        "validation_scope": "COLLECTOR_VALIDATED_SYNTHETIC",
    }


def _write_current_bundle(paths, rows):
    current_bytes = gzip.compress(_csv_bytes(CURRENT_FIELDS, rows), mtime=0)
    paths["current_evidence"].write_bytes(current_bytes)
    finite_by_endpoint = {
        contract["short"]: sum(
            row["endpoint"] == contract["short"] and row["status"] == "SUCCESS"
            for row in rows
        )
        for contract in ENDPOINTS.values()
    }
    receipt = {
        "status": "PASS_CURRENT_FINISHED_EXPLORATORY_TARGET_SNAPSHOT",
        "campaign_complete": False,
        "official_split": False,
        "promotion_performed": False,
        "publication_authorized": False,
        "output_counts": {
            "worker_evidence_rows": len(rows),
            "finite_target_values": sum(finite_by_endpoint.values()),
            "finite_target_values_by_endpoint": finite_by_endpoint,
        },
    }
    compact_receipt = (
        json.dumps(
            receipt,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")
    receipt["receipt_payload_sha256"] = hashlib.sha256(compact_receipt).hexdigest()
    paths["current_build_receipt"].write_bytes(_json_bytes(receipt))
    audit = {
        "status": "PASS_INDEPENDENT_CURRENT_FINISHED_TARGET_SNAPSHOT",
        "input_hashes": {
            "target_evidence.csv.gz": hashlib.sha256(current_bytes).hexdigest(),
            "BUILD_RECEIPT.json": _sha256(paths["current_build_receipt"]),
            "receipt_payload_sha256": receipt["receipt_payload_sha256"],
        },
    }
    compact_audit = (
        json.dumps(
            audit,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")
    audit["audit_payload_sha256"] = hashlib.sha256(compact_audit).hexdigest()
    paths["current_independent_audit"].write_bytes(_json_bytes(audit))


def _write_fixture(root):
    inputs = root / "inputs"
    inputs.mkdir()
    paths = {
        "release_manifest": inputs / "cif_manifest.csv",
        "base_manifest": inputs / "v26.0.1_base.manifest.jsonl",
        "base_summary": inputs / "v26.0.1_base.summary.json",
        "additions_manifest": inputs / "v26.0.2_additions.manifest.jsonl",
        "additions_summary": inputs / "v26.0.2_additions.summary.json",
        "current_evidence": inputs / "target_evidence.csv.gz",
        "current_build_receipt": inputs / "BUILD_RECEIPT.json",
        "current_independent_audit": inputs / "current.independent_audit.json",
    }
    structure_ids = (
        "BASE_EXISTING",
        "BASE_MISSING",
        "ADD_CURRENT",
        "ADD_EXCLUDED",
    )
    cif_hashes = {
        structure_id: hashlib.sha256(structure_id.encode("ascii")).hexdigest()
        for structure_id in structure_ids
    }
    paths["release_manifest"].write_bytes(
        _csv_bytes(
            ("structure_id", "sha256"),
            [
                {"structure_id": structure_id, "sha256": cif_hashes[structure_id]}
                for structure_id in structure_ids
            ],
        )
    )

    phases = (
        (
            "base",
            "v26.0.1_base",
            (("BASE_EXISTING", "EXISTING"), ("BASE_MISSING", "MISSING")),
        ),
        (
            "additions",
            "v26.0.2_additions",
            (("ADD_CURRENT", "MISSING"), ("ADD_EXCLUDED", "EXCLUDED")),
        ),
    )
    for label, phase, structures in phases:
        manifest_path = paths[label + "_manifest"]
        rows = [
            _source_row(structure_id, cif_hashes[structure_id], endpoint, phase, status)
            for structure_id, status in structures
            for endpoint in ENDPOINTS
        ]
        manifest_path.write_bytes(
            b"".join(
                json.dumps(
                    row,
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                    allow_nan=False,
                ).encode("utf-8")
                + b"\n"
                for row in rows
            )
        )
        status_counts = {
            endpoint: {
                status: sum(
                    row["endpoint"] == endpoint and row["status"] == status
                    for row in rows
                )
                for status in ("EXCLUDED", "EXISTING", "MISSING")
            }
            for endpoint in ENDPOINTS
        }
        excluded_count = sum(status == "EXCLUDED" for _, status in structures)
        summary = {
            "structure_count": len(structures),
            "eligible_structure_count": len(structures) - excluded_count,
            "excluded_structure_count": excluded_count,
            "record_count": len(rows),
            "alias_conflict_record_count": 0,
            "alias_conflict_structure_count": 0,
            "explicit_null_record_count": sum(row["explicit_null"] for row in rows),
            "status_counts": {
                status: sum(row["status"] == status for row in rows)
                for status in ("EXCLUDED", "EXISTING", "MISSING")
            },
            "endpoint_status_counts": status_counts,
            "artifacts": {"jsonl": {"sha256": _sha256(manifest_path)}},
        }
        paths[label + "_summary"].write_bytes(_json_bytes(summary))

    current_rows = [
        _current_row(
            "BASE_MISSING",
            cif_hashes["BASE_MISSING"],
            "CH4",
            "v26.0.1_base",
            "SUCCESS",
            4.2,
        ),
        _current_row(
            "BASE_MISSING",
            cif_hashes["BASE_MISSING"],
            "H2",
            "v26.0.1_base",
            "SUCCESS_NULL",
        ),
        _current_row(
            "BASE_MISSING",
            cif_hashes["BASE_MISSING"],
            "WIDOM_CO2_N2",
            "v26.0.1_base",
            "ERROR",
        ),
    ]
    for index, endpoint in enumerate(ENDPOINTS, start=1):
        current_rows.append(
            _current_row(
                "ADD_CURRENT",
                cif_hashes["ADD_CURRENT"],
                endpoint,
                "v26.0.2_additions",
                "SUCCESS",
                index + 0.25,
            )
        )
    _write_current_bundle(paths, current_rows)
    return paths, cif_hashes, current_rows


class CombinedTargetDatasetExampleTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.paths, self.cif_hashes, self.current_rows = _write_fixture(self.root)

    def tearDown(self):
        self.temporary_directory.cleanup()

    def _build(self, output, expected_success=True):
        command = [
            sys.executable,
            str(BUILD_SCRIPT),
            "--release-manifest",
            str(self.paths["release_manifest"]),
            "--base-manifest",
            str(self.paths["base_manifest"]),
            "--base-summary",
            str(self.paths["base_summary"]),
            "--additions-manifest",
            str(self.paths["additions_manifest"]),
            "--additions-summary",
            str(self.paths["additions_summary"]),
            "--current-evidence",
            str(self.paths["current_evidence"]),
            "--current-build-receipt",
            str(self.paths["current_build_receipt"]),
            "--current-independent-audit",
            str(self.paths["current_independent_audit"]),
            "--output-directory",
            str(output),
            "--snapshot-id",
            "synthetic_combined_targets_v1",
            "--data-cutoff-utc",
            "2026-09-04T00:00:00Z",
            "--expected-release-count",
            "4",
            "--expected-base-count",
            "2",
            "--expected-additions-count",
            "2",
            "--expected-finite",
            "ch4=3",
            "--expected-finite",
            "h2=2",
            "--expected-finite",
            "widom=1",
            "--expected-finite-total",
            "6",
        ]
        result = subprocess.run(command, check=False, capture_output=True, text=True)
        if expected_success and result.returncode != 0:
            self.fail("builder failed:\nstdout={}\nstderr={}".format(result.stdout, result.stderr))
        return result

    def _audit(self, dataset, audit_output, comparison=None):
        command = [
            sys.executable,
            str(AUDIT_SCRIPT),
            "--dataset",
            str(dataset),
            "--release-manifest",
            str(self.paths["release_manifest"]),
            "--base-manifest",
            str(self.paths["base_manifest"]),
            "--base-summary",
            str(self.paths["base_summary"]),
            "--additions-manifest",
            str(self.paths["additions_manifest"]),
            "--additions-summary",
            str(self.paths["additions_summary"]),
            "--current-evidence",
            str(self.paths["current_evidence"]),
            "--current-build-receipt",
            str(self.paths["current_build_receipt"]),
            "--current-independent-audit",
            str(self.paths["current_independent_audit"]),
            "--audit-output",
            str(audit_output),
            "--audit-id",
            "synthetic-independent-audit-v1",
            "--expected-release-count",
            "4",
            "--expected-finite",
            "ch4=3",
            "--expected-finite",
            "h2=2",
            "--expected-finite",
            "widom=1",
            "--expected-finite-total",
            "6",
        ]
        if comparison is not None:
            command.extend(("--comparison-dataset", str(comparison)))
        return subprocess.run(command, check=False, capture_output=True, text=True)

    def test_fill_only_merge_is_complete_deterministic_and_auditable(self):
        first = self.root / "combined-first"
        second = self.root / "combined-second"
        self._build(first)
        self._build(second)

        self.assertEqual(
            sorted(path.name for path in first.iterdir()),
            sorted(path.name for path in second.iterdir()),
        )
        for first_path in first.iterdir():
            self.assertEqual(first_path.read_bytes(), (second / first_path.name).read_bytes())

        with (first / "targets_for_attachment.csv").open(
            "r", encoding="utf-8", newline=""
        ) as handle:
            wide = list(csv.DictReader(handle))
        self.assertEqual(len(wide), 4)
        self.assertEqual(
            [row["structure_id"] for row in wide],
            sorted(self.cif_hashes),
        )
        by_id = {row["structure_id"]: row for row in wide}

        excluded = by_id["ADD_EXCLUDED"]
        for contract in ENDPOINTS.values():
            target = contract["target"]
            self.assertEqual(excluded[target], "")
            self.assertEqual(
                excluded[target + "__status"],
                "EXCLUDED_SCIENTIFIC_ELIGIBILITY",
            )
            self.assertEqual(excluded[target + "__eligible"], "false")

        base_missing = by_id["BASE_MISSING"]
        self.assertEqual(
            base_missing[ENDPOINTS["CH4"]["target"]],
            "4.2",
        )
        self.assertEqual(
            base_missing[ENDPOINTS["H2"]["target"]],
            "",
        )
        self.assertEqual(
            base_missing[ENDPOINTS["H2"]["target"] + "__status"],
            "SUCCESS_NULL",
        )
        self.assertEqual(
            base_missing[ENDPOINTS["WIDOM_CO2_N2"]["target"] + "__status"],
            "ERROR",
        )
        self.assertEqual(
            by_id["BASE_EXISTING"][ENDPOINTS["WIDOM_CO2_N2"]["target"]],
            "",
        )
        self.assertEqual(
            by_id["BASE_EXISTING"][
                ENDPOINTS["WIDOM_CO2_N2"]["target"] + "__status"
            ],
            "HISTORICAL_SCIENTIFIC_NULL",
        )

        with gzip.open(first / "target_assignments.csv.gz", "rt", encoding="utf-8") as handle:
            long_rows = list(csv.DictReader(handle))
        self.assertEqual(len(long_rows), 12)
        self.assertEqual(
            sum(row["source_class"] == "HISTORICAL_EXISTING" for row in long_rows),
            3,
        )
        self.assertEqual(
            sum(row["available"] == "true" for row in long_rows),
            6,
        )

        coverage = json.loads((first / "coverage_summary.json").read_text(encoding="utf-8"))
        self.assertEqual(coverage["release_structure_count"], 4)
        self.assertEqual(coverage["finite_endpoint_assignments"], 6)
        self.assertEqual(
            {
                endpoint: summary["finite_unique_structures"]
                for endpoint, summary in coverage["endpoint"].items()
            },
            {"ch4": 3, "h2": 2, "widom": 1},
        )

        from CoREMOF.targets import load_target_config

        loaded_config = load_target_config(first / "targets.json")
        self.assertEqual(len(loaded_config["sources"]), 1)
        self.assertEqual(
            tuple(loaded_config["sources"][0].target_columns),
            tuple(contract["target"] for contract in ENDPOINTS.values()),
        )

        audit_path = self.root / "combined-first.independent-audit.json"
        audit_result = self._audit(first, audit_path, comparison=second)
        if audit_result.returncode != 0:
            self.fail(
                "auditor failed:\nstdout={}\nstderr={}".format(
                    audit_result.stdout, audit_result.stderr
                )
            )
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        self.assertEqual(
            audit["status"],
            "PASS_INDEPENDENT_COMBINED_CURRENT_AVAILABLE_TARGET_DATASET",
        )
        self.assertEqual(audit["counts"]["assignment_count"], 12)
        self.assertEqual(audit["counts"]["finite_endpoint_assignments"], 6)

    def test_builder_refuses_to_fill_an_existing_assignment(self):
        rows = list(self.current_rows)
        rows.append(
            _current_row(
                "BASE_EXISTING",
                self.cif_hashes["BASE_EXISTING"],
                "CH4",
                "v26.0.1_base",
                "SUCCESS",
                123.0,
            )
        )
        _write_current_bundle(self.paths, rows)
        output = self.root / "must-not-exist"
        result = self._build(output, expected_success=False)
        self.assertEqual(result.returncode, 2)
        self.assertIn("not fill-only", result.stderr)
        self.assertFalse(output.exists())

    def test_builder_rejects_a_nonfinite_current_value(self):
        rows = [dict(row) for row in self.current_rows]
        rows[0]["value"] = "nan"
        _write_current_bundle(self.paths, rows)
        output = self.root / "must-not-exist"
        result = self._build(output, expected_success=False)
        self.assertEqual(result.returncode, 2)
        self.assertIn("not finite", result.stderr)
        self.assertFalse(output.exists())

    def test_independent_auditor_rejects_a_tampered_dataset(self):
        output = self.root / "combined"
        self._build(output)
        with (output / "targets_for_attachment.csv").open("a", encoding="utf-8") as handle:
            handle.write("\n")
        audit_path = self.root / "tampered.independent-audit.json"
        result = self._audit(output, audit_path)
        self.assertEqual(result.returncode, 2)
        self.assertIn("checksum mismatch", result.stderr)
        self.assertFalse(audit_path.exists())


if __name__ == "__main__":
    unittest.main()
