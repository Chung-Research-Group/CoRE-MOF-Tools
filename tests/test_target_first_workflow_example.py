"""Execute the current target-first example against a small validated release."""

import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from CoREMOF.parents import ZEO_NUMERIC_FINGERPRINT_DEFINITION
from test_benchmarks import EXACT_BENCHMARK_BACKEND
from test_dataset_labels import (
    CHECKER_COLUMNS, _declare_crystalnets_reference_contracts, _make_release,
    _metadata_rows, _parent_rows, _write_csv,
)

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "build_target_first_benchmark.py"
SPEC = importlib.util.spec_from_file_location("target_first_workflow_example", EXAMPLE)
workflow = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(workflow)
ENDPOINTS = ("ch4_loading", "h2_loading", "widom_ratio")


def _release_and_targets(root):
    release = root / "release"
    release.mkdir()
    _make_release(release)
    metadata, parents, manifest = [], [], []
    for index in range(60):
        sid = "ASR-COD-2026-{:04d}".format(index + 1)
        label = "CR" if index < 40 else "NCR" if index < 56 else "AMBIGUOUS" if index < 58 else "UNCHECKED"
        votes = (("PASS",) * 5 if label == "CR" else ("FAIL",) * 5 if label == "NCR"
                 else ("PASS", "FAIL", "PASS", "FAIL", "PASS") if label == "AMBIGUOUS"
                 else ("NOT_AVAILABLE",) * 5)
        row = dict(_metadata_rows()[0], structure_id=sid, cif_file="cifs/{}.cif".format(sid),
                   source_id="SOURCE-{}".format(index), mofid_v2="MOF-{}".format(index))
        row.update(dict(zip(CHECKER_COLUMNS, votes)))
        row.update({"label_{}checker".format(n): label for n in (3, 4, 5)})
        metadata.append(row)
        parent = {"structure_id": sid}
        for key in _parent_rows()[0]:
            if not key.endswith("_group"):
                continue
            prefix = key[:-6]
            group_prefix = _parent_rows()[0][key].split("-")[0]
            group_index = 1 if prefix == "rac" and index == 2 else index
            matched = prefix == "rac" and index in (1, 2)
            parent.update({prefix + "_group": "{}-{:08X}".format(group_prefix, group_index),
                           prefix + "_size": "2" if matched else "1",
                           prefix + "_status": "MATCHED" if matched else "UNMATCHED"})
        for prefix, group_prefix in (("rac_crystalnets", "RT"), ("mofid2_crystalnets", "M2T")):
            matched = prefix == "rac_crystalnets" and index in (1, 2)
            group_index = 1 if matched else index
            parent.update({prefix + "_group": "{}-{:08X}".format(group_prefix, group_index),
                           prefix + "_size": "2" if matched else "1",
                           prefix + "_status": "MATCHED" if matched else "UNMATCHED"})
        parents.append(parent)
        # 0--1 by CIF, 1--2 by RAC/RT: row 1 has a missing target.
        # 3--40 links opposite labels; row 40 also has a missing target.
        hash_key = 0 if index in (0, 1) else 3 if index in (3, 40) else index
        manifest.append({"structure_id": sid, "cif_file": row["cif_file"],
                         "size_bytes": "1", "sha256": hashlib.sha256(str(hash_key).encode()).hexdigest()})
    _write_csv(release / "metadata/metadata.csv", tuple(metadata[0]), metadata)
    _write_csv(release / "parent_groups/parent_groups.csv", tuple(parents[0]), parents)
    _write_csv(release / "manifests/cif_manifest.csv", tuple(manifest[0]), manifest)
    info_path = release / "dataset_info.json"
    info = json.loads(info_path.read_text())
    info["structure_count"] = len(metadata)
    info_path.write_text(json.dumps(info), encoding="utf-8")
    _declare_crystalnets_reference_contracts(release)
    ids = [row["structure_id"] for row in metadata]
    observations = [{"structure_id": sid, **{key: float(index) for key in ENDPOINTS}}
                    for index, sid in enumerate(ids)]
    for index in (1, 40):
        observations[index][ENDPOINTS[0]] = None
    observations[4][ENDPOINTS[0]] = True
    observations[5][ENDPOINTS[0]] = "1.25"
    (root / "targets.json").write_text(json.dumps(observations), encoding="utf-8")
    source = {"path": "targets.json", "id_column": "structure_id", "target_columns": list(ENDPOINTS),
              "units": {key: "dimensionless" if key == ENDPOINTS[2] else "mol/kg" for key in ENDPOINTS},
              "conditions": {key: {"temperature_K": 298} for key in ENDPOINTS}}
    config = root / "target_config.json"
    config.write_text(json.dumps({"sources": [source]}), encoding="utf-8")
    # Native JSON values preserve bool/text for eligibility to reject. CSV users
    # should declare numeric types; the merge rejects invalid conversions itself.
    return release, config, ids


def _feature_tables(release, ids):
    rac_fields = tuple("rac_{:03d}".format(i) for i in range(264))
    rac, zeo, topology = [], [], []
    zeo_fields = tuple(ZEO_NUMERIC_FINGERPRINT_DEFINITION["numeric_fields"]) + (
        "n2_channel_dimension", "structure_periodic_dimension")
    for index, sid in enumerate(ids):
        value_index = 1 if index == 2 else index
        rac.append({"structure_id": sid, "rac5_available": "true" if index < 30 else "false",
                    **{name: str((value_index + 1) * (j % 7 + 1)) if index < 30 else ""
                       for j, name in enumerate(rac_fields)}})
        zeo.append({"structure_id": sid, "n2_he_available": "true" if 30 <= index < 50 else "false",
                    "periodicity_available": "true" if 30 <= index < 50 else "false",
                    **{name: str((index + 1) / (j + 2)) if 30 <= index < 50 else ""
                       for j, name in enumerate(zeo_fields)}})
        topology.append({"structure_id": sid, "topology_available": "true", "network_dimension": "3",
                         "single_node_net": "pcu", "all_node_net": "pcu", "single_all_agree": "true"})
    for name, rows in (("rac5", rac), ("zeo", zeo), ("topology", topology)):
        _write_csv(release / "features" / (name + "_features.csv"), tuple(rows[0]), rows)


class TargetFirstWorkflowExampleTests(unittest.TestCase):
    def test_finite_rule_preserves_zero_and_rejects_invalid_values(self):
        for value in (0, -0.0, 1.2, 10 ** 400):
            self.assertIsNone(workflow.finite_target_reason(value))
        for value in (None, True, False, "1.0", float("nan"), float("inf"), -float("inf")):
            self.assertIsNotNone(workflow.finite_target_reason(value))

    def check_workflow(self, root, diversity):
        release, config, ids = _release_and_targets(root)
        if diversity == "representative":
            _feature_tables(release, ids)
        output = root / "output"
        args = [str(release), "--target-config", str(config), "--output-directory", str(output),
                "--seeds", "912", "913", "--diversity", diversity]
        for endpoint in ENDPOINTS:
            args.extend(("--require-target", endpoint))
        interpreter_flags = ["-S"] if sys.flags.no_site else []
        completed = subprocess.run([sys.executable, *interpreter_flags, str(EXAMPLE), *args], text=True,
                                   capture_output=True, timeout=120)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        receipt = json.loads((output / "workflow_receipt.json").read_text())
        eligible = receipt["eligibility"]["eligible_ids"]
        self.assertIn(ids[0], eligible)  # zero is available
        self.assertEqual(set(ids).difference(eligible), {ids[1], ids[4], ids[5], ids[40]})
        self.assertEqual(receipt["eligibility"]["eligible_ids_sha256"], workflow.canonical_digest(eligible))
        suite_path = output / receipt["suite_receipt"]["path"]
        suite = json.loads(suite_path.read_text())
        self.assertEqual(suite["cohort_receipt"]["eligibility_filter"]["eligible_ids_sha256"],
                         receipt["eligibility"]["eligible_ids_sha256"])
        self.assertEqual(suite["cohort_receipt"]["eligible_pool_counts"], {"C_CR": 36, "M_NCR": 15})
        self.assertEqual(receipt["suite_assignment_sha256"], receipt["target_attachment_assignment_sha256"])
        attached_path = output / receipt["target_attachment_receipt"]["path"]
        attached = json.loads(attached_path.read_text())
        self.assertEqual(attached["original_suite_assignment_sha256"], receipt["suite_assignment_sha256"])
        self.assertEqual(workflow.canonical_digest(receipt["required_endpoint_definitions"]),
                         receipt["required_endpoint_definitions_sha256"])
        for relative in ("target_eligibility.csv",):
            self.assertTrue((output / relative).is_file())
        tests = set()
        by_seed = {}
        for path in sorted((suite_path.parent / "runs").glob("*.csv")):
            with path.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            assignment = {row["structure_id"]: row["partition"] for row in rows}
            self.assertNotIn(ids[3], assignment)  # missing opposite label cannot hide impurity
            self.assertNotIn(ids[40], assignment)
            self.assertEqual(assignment.get(ids[0]), assignment.get(ids[2]))
            tests.add(tuple(sorted(sid for sid, part in assignment.items() if part == "test")))
            by_seed.setdefault(rows[0]["seed"], []).append(assignment)
        self.assertEqual(len(tests), 1)
        for assignments in by_seed.values():
            for first, second in zip(assignments, assignments[1:]):
                self.assertTrue(all(first[sid] == second[sid] for sid in set(first).intersection(second)))
        for line in (output / "SHA256SUMS").read_text().splitlines():
            digest, path = line.split("  ", 1)
            self.assertEqual(hashlib.sha256((output / path).read_bytes()).hexdigest(), digest)
        before = (output / "workflow_receipt.json").read_bytes()
        with self.assertRaises(FileExistsError):
            workflow.build_target_first_benchmark(release, config, output, required_targets=ENDPOINTS)
        self.assertEqual((output / "workflow_receipt.json").read_bytes(), before)
        if diversity == "none":
            # Changing required endpoints creates a separately receipted dataset;
            # it cannot silently reinterpret the earlier eligibility snapshot.
            second_path = workflow.build_target_first_benchmark(
                release, config, root / "second", required_targets=ENDPOINTS[1:],
                diversity="none", seeds=(912,), group_criteria=("RT", "M2T"))
            second = json.loads(second_path.read_text())
            self.assertEqual(second["eligibility"]["eligible_count"], 60)
            self.assertNotEqual(second["eligibility"]["eligible_ids_sha256"],
                                receipt["eligibility"]["eligible_ids_sha256"])
            self.assertNotEqual(second["required_endpoint_definitions_sha256"],
                                receipt["required_endpoint_definitions_sha256"])
            self.assertEqual(second["target_merge_receipt_sha256"], receipt["target_merge_receipt_sha256"])
            self.assertEqual((output / "workflow_receipt.json").read_bytes(), before)
        return receipt

    def test_loader_merge_target_first_export_and_receipts(self):
        with tempfile.TemporaryDirectory() as temporary:
            self.check_workflow(Path(temporary), "none")

    @unittest.skipUnless(EXACT_BENCHMARK_BACKEND, "requires the pinned benchmark extra")
    def test_representative_rt_m2t_workflow_executes(self):
        with tempfile.TemporaryDirectory() as temporary:
            receipt = self.check_workflow(Path(temporary), "representative")
            suite = json.loads((Path(temporary) / "output" / receipt["suite_receipt"]["path"]).read_text())
            self.assertEqual(suite["cohort_receipt"]["diversity_profile"]["backend"]["numpy"], "1.26.4")

    def test_changed_required_endpoint_is_not_silently_accepted(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            release, config, _ = _release_and_targets(root)
            with self.assertRaisesRegex(ValueError, "unknown required"):
                workflow.build_target_first_benchmark(release, config, root / "bad",
                    required_targets=("wrong_endpoint",), diversity="none")
            self.assertFalse((root / "bad").exists())


if __name__ == "__main__":
    unittest.main()
