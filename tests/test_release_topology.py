"""Release topology tests that need neither Julia nor scientific dependencies."""
import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

from CoREMOF import release_topology as topology


SID = "FSR-COD-2016-0106"
SHA = "a" * 64


def record():
    raw = topology._unavailable(SID, SHA, "fixture", "fixture", 0.0)
    part = {"status": "SUCCESS", "dimension": 3, "topology": "pcu", "is_named_net": True}
    raw.update(execution_status="SUCCESS", error_type=None, error_message=None,
               interpenetrated_subnet_count=1, catenation_degree=1,
               network_dimension=3, single_node_net="pcu", all_node_net="pcu",
               all_subnets_single_all_agree=True,
               subnets=[{"subnet_index": 1, "single_node": dict(part), "all_node": dict(part),
                         "single_all_agree": True}])
    return copy.deepcopy(raw)


class ReleaseTopologyRecordTests(unittest.TestCase):
    def test_standard_library_import(self):
        process = subprocess.run([sys.executable, "-B", "-S", "-c",
            "from CoREMOF import release_topology; import sys; "
            "assert not ({'juliacall','numpy','pymatgen'} & set(sys.modules))"],
            cwd=Path(__file__).resolve().parents[1], capture_output=True, text=True)
        self.assertEqual(process.returncode, 0, process.stderr)

    def test_vendored_scientific_protocol_bytes(self):
        topology._check(topology.projection.__file__, topology.PROJECTION_SHA256)
        root = Path(topology.__file__).with_name("_topology_profile")
        for name, digest in topology.PROFILE_FILES.items():
            topology._check(root / name, digest)
        for name, digest in topology.PROJECT_FILES.items():
            topology._check(root / name, digest)

    def test_complete_both_modes(self):
        result = topology._validate(record(), SID, SHA)
        self.assertTrue(result["topology_available"])
        self.assertEqual(result["execution_status"], "SUCCESS")
        self.assertEqual(result["single_node_net"], "pcu")
        self.assertTrue(result["single_all_agree"])

    def test_mode_disagreement_keeps_complete_custom_genome(self):
        raw = record()
        raw["subnets"][0]["all_node"]["topology"] = "3-abcdefghijkl (3 1 1 0 0 1)"
        result = topology._validate(raw, SID, SHA)
        self.assertTrue(result["topology_available"])
        self.assertFalse(result["single_all_agree"])
        self.assertEqual(result["all_node_net"], "3-abcdefghijkl")
        self.assertEqual(result["subnets"][0]["all_node"]["topological_genome"], "3 1 1 0 0 1")

    def test_failed_component_keeps_other_mode_and_is_unavailable_for_matching(self):
        raw = record()
        raw["subnets"][0]["single_node"].update(dimension=0, topology="FAILED with: AssertionError(test)")
        result = topology._validate(raw, SID, SHA)
        self.assertEqual(result["execution_status"], "PARTIAL")
        self.assertFalse(result["topology_available"])
        self.assertIsNone(result["network_dimension"])
        self.assertIsNone(result["single_node_net"])
        self.assertIsNone(result["all_node_net"])
        self.assertEqual(result["subnets"][0]["all_node"]["topology_key"], "pcu")
        self.assertEqual(result["subnets"][0]["single_node"]["status"], "ERROR")

    def test_zero_dimensional_and_multiple_subnets_retained(self):
        raw = record()
        raw["subnets"].append(copy.deepcopy(raw["subnets"][0]))
        raw["subnets"][1]["subnet_index"] = 2
        raw.update(interpenetrated_subnet_count=2, catenation_degree=2)
        for subnet in raw["subnets"]:
            for mode in ("single_node", "all_node"):
                subnet[mode].update(dimension=0, topology="0-dimensional")
        result = topology._validate(raw, SID, SHA)
        self.assertEqual(result["network_dimension"], 0)
        self.assertEqual(result["catenation_degree"], 2)
        self.assertEqual(len(result["subnets"]), 2)

    def test_empty_success_remains_partial(self):
        raw = record()
        raw.update(subnets=[], interpenetrated_subnet_count=0, catenation_degree=0)
        result = topology._validate(raw, SID, SHA)
        self.assertEqual(result["execution_status"], "PARTIAL")
        self.assertFalse(result["topology_available"])

    def test_different_identity_version_method_and_bad_numbers_rejected(self):
        for key, value in (("structure_id", "wrong"), ("cif_sha256", "b" * 64),
                           ("software", {}), ("method", {}), ("runtime_seconds", float("nan")),
                           ("catenation_degree", True), ("execution_status", "UNKNOWN")):
            raw = record()
            raw[key] = value
            with self.subTest(key=key), self.assertRaises(topology.ReleaseTopologyError):
                topology._validate(raw, SID, SHA)

    def test_inconsistent_order_or_dimensions_rejected(self):
        raw = record()
        raw["subnets"][0]["subnet_index"] = 2
        with self.assertRaises(topology.ReleaseTopologyError):
            topology._validate(raw, SID, SHA)
        raw = record()
        raw["subnets"][0]["all_node"]["dimension"] = float("nan")
        with self.assertRaises(topology.ReleaseTopologyError):
            topology._validate(raw, SID, SHA)

    def test_error_cannot_keep_scientific_values(self):
        raw = topology._unavailable(SID, SHA, "TIMEOUT", "timeout", 2.0)
        result = topology._validate(raw, SID, SHA)
        self.assertFalse(result["topology_available"])
        raw["network_dimension"] = 3
        with self.assertRaises(topology.ReleaseTopologyError):
            topology._validate(raw, SID, SHA)


class ReleaseTopologyRuntimeTests(unittest.TestCase):
    def test_unreadable_tree_cannot_produce_an_incomplete_manifest(self):
        with tempfile.TemporaryDirectory() as directory:
            def broken_walk(root, *, followlinks, onerror):
                onerror(PermissionError("fixture unreadable directory"))
                return iter(())
            with patch.object(topology.os, "walk", side_effect=broken_walk), \
                 self.assertRaisesRegex(topology.ReleaseTopologyError, "Unreadable runtime tree"):
                topology._tree(Path(directory))

    def test_runtime_bytes_and_membership_are_bound(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            julia = root / "julia/bin/julia"
            julia.parent.mkdir(parents=True)
            julia.write_text("fixture executable")
            depot = root / "depot"
            for name in ("packages", "compiled", "artifacts"):
                (depot / name).mkdir(parents=True)
            manifest = root / "runtime.json"
            digest = topology.capture_runtime_manifest(julia=julia, depot=depot, output_path=manifest)
            topology._verify_runtime(julia, depot, manifest, digest)
            with self.assertRaises(FileExistsError):
                topology.capture_runtime_manifest(julia=julia, depot=depot, output_path=manifest)
            julia.write_text("changed executable")
            with self.assertRaisesRegex(topology.ReleaseTopologyError, "Runtime tree differs"):
                topology._verify_runtime(julia, depot, manifest, digest)
            julia.write_text("fixture executable")
            extra = depot / "compiled/unrecorded"
            extra.write_text("unrecorded code")
            with self.assertRaisesRegex(topology.ReleaseTopologyError, "Runtime tree differs"):
                topology._verify_runtime(julia, depot, manifest, digest)

    def test_runtime_symlink_escape_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inside = root / "runtime"
            inside.mkdir()
            (root / "outside").write_text("keep")
            (inside / "link").symlink_to(root / "outside")
            with self.assertRaisesRegex(topology.ReleaseTopologyError, "Escaping/broken"):
                topology._tree(inside)


class ReleaseTopologyIsolationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.cif = self.root / "original.cif"
        self.cif.write_text("data_example\n")
        self.project = self.root / "project"
        self.project.mkdir()
        for name in topology.PROJECT_FILES:
            (self.project / name).write_text("fixture")
        self.options = dict(cif_path=self.cif, structure_id=SID,
            output_dir=self.root / "output", julia=Path(sys.executable), project=self.project,
            depot=self.root, runtime_manifest=self.cif, runtime_manifest_sha256=SHA, timeout_seconds=0.1)

    def test_existing_output_preserved_before_runtime_or_science(self):
        destination = self.options["output_dir"]
        destination.mkdir()
        with patch.object(topology, "_verify_runtime") as verify, \
             patch.object(topology.subprocess, "Popen") as launch, \
             self.assertRaises(FileExistsError):
            topology.calculate_release_topology(**self.options)
        verify.assert_not_called()
        launch.assert_not_called()

    def patches(self):
        assets = Path(topology.__file__).with_name("_topology_profile")
        return (patch.object(topology, "_profile", return_value=assets),
                patch.object(topology, "_verify_runtime"),
                patch.object(topology, "_check", side_effect=lambda path, digest: Path(path)))

    def test_child_isolation_and_failure_receipt(self):
        fake = types.SimpleNamespace(wait=lambda timeout: None, returncode=1)
        def launch(command, **kwargs):
            isolated_cif = Path(command[-3])
            self.assertNotEqual(isolated_cif, self.cif)
            self.assertEqual(isolated_cif.read_bytes(), self.cif.read_bytes())
            self.assertIn("--startup-file=no", command)
            self.assertIn("--history-file=no", command)
            self.assertEqual(kwargs["env"]["JULIA_PKG_OFFLINE"], "true")
            self.assertNotIn("LD_PRELOAD", kwargs["env"])
            self.assertNotIn("PYTHONPATH", kwargs["env"])
            self.assertTrue(kwargs["start_new_session"])
            return fake
        profile, verify, check = self.patches()
        with profile, verify, check, patch.object(topology.subprocess, "Popen", side_effect=launch):
            result = topology.calculate_release_topology(**self.options)
        self.assertFalse(result["topology_available"])
        self.assertEqual(result["error"]["type"], "PROCESS_ERROR")
        self.assertEqual(self.cif.read_text(), "data_example\n")
        receipt = json.loads((self.options["output_dir"] / "receipt.json").read_text())
        self.assertFalse(receipt["release_metadata_promoted"])

    def test_timeout_stops_only_owned_process_and_keeps_nulls(self):
        fake = types.SimpleNamespace(returncode=-15,
            wait=lambda timeout: (_ for _ in ()).throw(subprocess.TimeoutExpired("julia", 1)))
        profile, verify, check = self.patches()
        with profile, verify, check, patch.object(topology.subprocess, "Popen", return_value=fake), \
             patch.object(topology, "_stop") as stop:
            result = topology.calculate_release_topology(**self.options)
        stop.assert_called_once_with(fake)
        self.assertEqual(result["error"]["type"], "TIMEOUT")
        self.assertIsNone(result["single_node_net"])
        self.assertEqual(result["subnets"], [])

    def test_bad_id_and_timeout_rejected(self):
        for change in ({"structure_id": "../bad"}, {"timeout_seconds": 0},
                       {"timeout_seconds": float("inf")}, {"timeout_seconds": True}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                topology.calculate_release_topology(**{**self.options, **change})


if __name__ == "__main__":
    unittest.main()
