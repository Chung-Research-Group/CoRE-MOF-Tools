"""Release-method tests, independent of optional scientific installations."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

from CoREMOF import _release_mofid_protocol as protocol
from CoREMOF import release_mofid as replay


class ReleaseMOFidContractTests(unittest.TestCase):
    def test_vendored_protocol_is_byte_identical(self):
        self.assertEqual(hashlib.sha256(Path(protocol.__file__).read_bytes()).hexdigest(),
                         replay.PROTOCOL_SHA256)

    def test_standard_library_import(self):
        done = subprocess.run([sys.executable, "-S", "-B", "-c",
            "from CoREMOF import release_mofid; import sys; "
            "assert not ({'ase','numpy','pymatgen','mofid'} & set(sys.modules))"],
            cwd=Path(__file__).resolve().parents[1], capture_output=True, text=True)
        self.assertEqual(done.returncode, 0, done.stderr)

    def test_path_like_and_mismatched_ids_rejected(self):
        for name in ("../ASR-COD-2000-0001", "/tmp/example", "ASR-COD-2000-0001/..",
                     "FSR-COD-2000-0001", ""):
            with self.subTest(name=name), self.assertRaises(ValueError):
                replay._safe_id(name, "ASR")
        replay._safe_id("ASR-COD-UNKN-0001", "ASR")

    def test_profile_modified_manifest_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "method.json"
            path.write_text("{}")
            with self.assertRaisesRegex(replay.ReleaseMOFidError, "SHA-256"):
                replay._read_profile(path, path)

    def test_topology_failure_is_not_execution_failure(self):
        for token in ("UNKNOWN", "ERROR", "TIMEOUT"):
            value, status = protocol.normalize_mofid_v1(
                f"[Zn].[O] MOFid-v1.{token}.cat0.NO_REF;record")
            self.assertEqual(value, f"[Zn].[O] MOFid-v1.{token}.cat0")
            self.assertEqual(status, f"SUCCESS_TOPOLOGY_{token}")
        self.assertEqual(protocol.normalize_mofid_v1("* MOFid-v1.NA.NO_REF"),
                         (None, "NOT_AVAILABLE_NO_MOF"))

    def test_pinned_coordinate_convention_is_preserved(self):
        builder = object.__new__(protocol.PinnedV2Builder)
        builder.Lattice = lambda cell: cell
        captured = {}
        def structure(*args, **kwargs):
            captured.update(args=args, kwargs=kwargs)
            return "structure"
        builder.Structure = structure
        atoms = types.SimpleNamespace(cell="cell", get_chemical_symbols=lambda: ["Zn"],
                                      get_positions=lambda: [[1, 2, 3]])
        self.assertEqual(builder._to_pymatgen(atoms), "structure")
        self.assertEqual(captured["kwargs"], {})
        self.assertEqual(captured["args"], ("cell", ["Zn"], [[1, 2, 3]]))

    def test_frozen_order_first_match_and_all_matches_recorded(self):
        builder = object.__new__(protocol.PinnedV2Builder)
        builder.library = {"Zn1": [
            {"node_label": "Zn1_Type-9", "archive_index": 2, "path": "nine"},
            {"node_label": "Zn1_Type-1", "archive_index": 8, "path": "one"},
        ]}
        options = {}
        def matcher(**kwargs):
            options.update(kwargs)
            return types.SimpleNamespace(fit=lambda left, right: True)
        builder.StructureMatcher = matcher
        builder.ElementComparator = lambda: "elements"
        builder._to_pymatgen = lambda atoms: atoms
        builder._remove_pbc_cuts = lambda atoms: atoms
        builder.ase_read = lambda path: path
        result = builder.match_node("Zn1", "query")
        self.assertEqual(result["selected_node_label"], "Zn1_Type-9")
        self.assertEqual(result["matching_node_labels"], ["Zn1_Type-9", "Zn1_Type-1"])
        self.assertTrue(result["multiple_structural_matches"])
        self.assertEqual((options["ltol"], options["stol"]), (0.25, 1.5))
        self.assertFalse(options["scale"])

    def test_formula_retains_explicit_ones(self):
        self.assertEqual(protocol.formula_from_symbols(["Zn", "O", "O", "C"]), "C1O2Zn1")

    def test_timeout_does_not_revalidate_existing_v1(self):
        row = {"structure_id": "FSR-COD-2000-0001", "structure_variant": "FSR",
               "cif_file": "input.cif", "existing_mofid_v1": "[Zn] MOFid-v1.pcu.cat0"}
        result = protocol.timeout_record(row, "a"*64, "b"*64, "c"*64, "d"*64, 30)
        protocol.validate_result_record(result, row["structure_id"])
        self.assertEqual(result["v1_comparison"], "NOT_REVALIDATED_TIMEOUT")
        self.assertEqual(result["mofid_v2_status"], "TIMEOUT")
        self.assertEqual(result["mofid_v2_scope"], "COREMOF_FSR_EXTENSION")


class ReleaseMOFidIsolationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.cif = self.root / "source.cif"
        self.cif.write_bytes(b"data_test\n")
        self.options = dict(cif_path=self.cif, structure_id="ASR-COD-2000-0001",
            structure_variant="ASR", output_dir=self.root / "result", python=sys.executable,
            method_manifest=self.cif, node_manifest=self.cif, node_root=self.root,
            source_root=self.root, pinned_site=self.root, mofid_site=self.root)

    def test_existing_output_is_not_touched(self):
        output = self.options["output_dir"]
        output.mkdir()
        (output / "sentinel").write_text("keep")
        with patch.object(replay, "_read_profile", return_value={}), \
             patch.object(replay.subprocess, "Popen") as process, \
             self.assertRaises(FileExistsError):
            replay.calculate_release_mofid(**self.options)
        process.assert_not_called()
        self.assertEqual((output / "sentinel").read_text(), "keep")

    def test_worker_failure_preserves_source_and_publishes_nothing(self):
        fake = types.SimpleNamespace(wait=lambda timeout: None, returncode=1)
        def launch(command, **kwargs):
            request = json.loads(Path(command[-1]).read_text())
            copied = Path(request['private_root']) / 'input' / 'ASR-COD-2000-0001.cif'
            self.assertEqual(copied.read_bytes(), self.cif.read_bytes())
            self.assertNotEqual(copied.resolve(), self.cif.resolve())
            return fake
        with patch.object(replay, "_read_profile", return_value={}), \
             patch.object(replay.subprocess, "Popen", side_effect=launch), \
             self.assertRaises(replay.ReleaseMOFidError):
            replay.calculate_release_mofid(**self.options)
        self.assertFalse(self.options["output_dir"].exists())
        self.assertEqual(self.cif.read_bytes(), b"data_test\n")
        self.assertEqual(sorted(path.name for path in self.root.iterdir()), ["source.cif"])

    def test_timeout_publishes_explicit_unverified_timeout(self):
        fake = types.SimpleNamespace(
            wait=lambda timeout: (_ for _ in ()).throw(subprocess.TimeoutExpired("worker", 1)),
            returncode=-15)
        with patch.object(replay, "_read_profile", return_value={}), \
             patch.object(replay.subprocess, "Popen", return_value=fake) as process, \
             patch.object(replay, "_stop_own_group") as stop:
            result = replay.calculate_release_mofid(**self.options)
        stop.assert_called_once_with(fake)
        self.assertTrue(process.call_args.kwargs["start_new_session"])
        self.assertEqual(result["mofid_v2_status"], "TIMEOUT")
        receipt = json.loads((self.options["output_dir"] / "receipt.json").read_text())
        self.assertFalse(receipt["runtime_verified"])
        self.assertFalse(receipt["release_promotion"])

    def test_invalid_timeout_rejected_before_io(self):
        for value in (0, -1, float("nan"), float("inf"), True):
            with self.subTest(value=value), self.assertRaises(ValueError):
                replay.calculate_release_mofid(**self.options, timeout_seconds=value)


if __name__ == "__main__":
    unittest.main()
