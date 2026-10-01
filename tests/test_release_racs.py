"""Exact-schema and isolation checks without installing the external runtime."""
import hashlib
import json
import re
from pathlib import Path
import subprocess
import stat
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

from CoREMOF import release_racs as racs


class ReleaseRACSchemaTests(unittest.TestCase):
    def test_release_guide_heading_underlines(self):
        root = Path(__file__).resolve().parents[1] / "docs/source"
        for name in ("release_mofid_replay.rst", "release_racs_replay.rst", "release_topology_replay.rst", "release_zeopp_replay.rst"):
            lines = (root / name).read_text().splitlines()
            for index, line in enumerate(lines[1:], 1):
                if re.fullmatch(r"[=-]{3,}", line):
                    with self.subTest(file=name, line=index + 1):
                        self.assertGreaterEqual(len(line), len(lines[index - 1].strip()))

    def test_standard_library_import(self):
        result = subprocess.run([sys.executable, "-S", "-B", "-c",
            "from CoREMOF import release_racs; import sys; "
            "assert not ({'numpy','molSimplify','pymatgen'} & set(sys.modules))"],
            cwd=Path(__file__).resolve().parents[1], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_frozen_176_and_264_column_schemas(self):
        for depth, count in ((3, 176), (5, 264)):
            with self.subTest(depth=depth):
                schema = racs._schema(depth)
                self.assertEqual(len(schema), count)
                self.assertEqual(list(schema), sorted(schema))
                digest = hashlib.sha256(json.dumps(sorted(schema.values()), separators=(",", ":")).encode()).hexdigest()
                self.assertEqual(digest, racs.SCHEMA_SHA256[depth])
        self.assertTrue(set(racs._schema(3)) < set(racs._schema(5)))

    def test_unauthorized_depth_rejected(self):
        for depth in (0, 4, 6, True, 5.0, "5"):
            with self.subTest(depth=depth), self.assertRaises(ValueError):
                racs._schema(depth)

    def test_input_order_normalized_without_rounding(self):
        schema = racs._schema(5)
        names = list(schema)[::-1]
        precise = float.fromhex("0x1.23456789abcdep-7")
        values = [precise] * len(names)
        values[0] = -0.0
        values[1] = 0.0
        result = racs._validated_values(names, values, 5)
        self.assertEqual(list(result), list(schema.values()))
        self.assertEqual(result[schema[names[2]]].hex(), precise.hex())
        self.assertEqual(result[schema[names[0]]].hex(), "-0x0.0p+0")
        self.assertEqual(result[schema[names[1]]].hex(), "0x0.0p+0")

    def test_partial_duplicate_extra_and_nonfinite_vectors_rejected(self):
        names = list(racs._schema(5))
        for bad_names, values in (
            (names[:-1], [0.0] * 263),
            (names[:-1] + names[:1], [0.0] * 264),
            (names + ["extra"], [0.0] * 265),
            (names, [0.0] * 263),
            (names, [float("nan")] + [0.0] * 263),
            (names, [float("inf")] + [0.0] * 263),
        ):
            with self.subTest(length=len(bad_names)), self.assertRaises(racs.ReleaseRACError):
                racs._validated_values(bad_names, values, 5)

    def test_unsafe_archive_paths_rejected(self):
        for name in ("../outside", "/etc/config", "a/../../b", "a//b"):
            with self.subTest(name=name), self.assertRaises(racs.ReleaseRACError):
                racs._relative_path(name)

    def test_failed_result_cannot_contain_a_partial_vector(self):
        result = self.result()
        result.update(execution_status="ERROR", available=False, descriptors={"partial": 0},
                      values_float_hex=None, error={"type": "TestError"})
        with self.assertRaisesRegex(racs.ReleaseRACError, "null entire vector"):
            racs._validate_result(result, "2000[Cu][nan]3[ASR]1", "a" * 64, 5)

    def test_success_float_hex_and_identity_checked(self):
        result = self.result()
        racs._validate_result(result, "2000[Cu][nan]3[ASR]1", "a" * 64, 5)
        result["values_float_hex"][next(iter(result["descriptors"]))] = "0x1.0p+0"
        with self.assertRaisesRegex(racs.ReleaseRACError, "representation differ"):
            racs._validate_result(result, "2000[Cu][nan]3[ASR]1", "a" * 64, 5)
        with self.assertRaisesRegex(racs.ReleaseRACError, "identity/schema"):
            racs._validate_result(result, "2000[Cu][nan]3[ASR]2", "a" * 64, 5)

    @staticmethod
    def result():
        schema = racs._schema(5)
        return {"schema": "coremof-release-racs/1.0", "profile": racs.PROFILE,
            "structure_id": "2000[Cu][nan]3[ASR]1", "cif_sha256": "a" * 64,
            "depth": 5, "schema_sha256": racs.SCHEMA_SHA256[5],
            "metadata_names": list(schema.values()), "call_arguments": racs.CALL_ARGUMENTS,
            "execution_status": "SUCCESS", "available": True,
            "descriptors": {name: 0.0 for name in schema.values()},
            "values_float_hex": {name: "0x0.0p+0" for name in schema.values()}}


class ReleaseRACIsolationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.cif = self.root / "source.cif"
        self.cif.write_bytes(b"data_test\n")
        self.options = dict(cif_path=self.cif, structure_id="2000[Cu][nan]3[ASR]1",
            output_dir=self.root / "result", python=sys.executable, source_archive=self.cif,
            environment_manifest=self.cif, timeout_seconds=1)

    @staticmethod
    def checked(path, *args, **kwargs):
        return Path(path)

    def test_existing_output_preserved(self):
        output = self.options["output_dir"]
        output.mkdir()
        (output / "sentinel").write_text("keep")
        with patch.object(racs, "_check", side_effect=self.checked), \
             patch.object(racs.subprocess, "Popen") as process, \
             self.assertRaises(FileExistsError):
            racs.calculate_release_racs(**self.options)
        process.assert_not_called()
        self.assertEqual((output / "sentinel").read_text(), "keep")

    def test_child_has_clean_environment_and_isolated_input(self):
        fake = types.SimpleNamespace(wait=lambda timeout: None, returncode=1)
        def launch(command, **kwargs):
            request = json.loads(Path(command[-1]).read_text())
            copied = Path(request['private_root']) / 'input' / '2000[Cu][nan]3[ASR]1.cif'
            self.assertEqual(copied.read_bytes(), self.cif.read_bytes())
            self.assertNotEqual(copied.resolve(), self.cif.resolve())
            self.assertNotIn("LD_PRELOAD", kwargs["env"])
            self.assertNotIn("PYTHONPATH", kwargs["env"])
            self.assertEqual(kwargs["env"]["PYTHONHASHSEED"], "0")
            self.assertTrue(kwargs["start_new_session"])
            self.assertIn("-S", command)
            return fake
        with patch.object(racs, "_check", side_effect=self.checked), \
             patch.object(racs.subprocess, "Popen", side_effect=launch), \
             self.assertRaises(racs.ReleaseRACError):
            racs.calculate_release_racs(**self.options)
        self.assertEqual(self.cif.read_bytes(), b"data_test\n")
        self.assertFalse(self.options["output_dir"].exists())
        self.assertEqual([p.name for p in self.root.iterdir()], ["source.cif"])

    def test_timeout_is_explicit_null_vector_not_imputation(self):
        fake = types.SimpleNamespace(returncode=-15,
            wait=lambda timeout: (_ for _ in ()).throw(subprocess.TimeoutExpired("worker", 1)))
        with patch.object(racs, "_check", side_effect=self.checked), \
             patch.object(racs.subprocess, "Popen", return_value=fake), \
             patch.object(racs, "_stop") as stop:
            result = racs.calculate_release_racs(**self.options)
        stop.assert_called_once_with(fake)
        self.assertEqual(result["execution_status"], "TIMEOUT")
        self.assertIsNone(result["descriptors"])
        self.assertIsNone(result["values_float_hex"])
        receipt = json.loads((self.options["output_dir"] / "receipt.json").read_text())
        self.assertFalse(receipt["runtime_verified"])
        self.assertFalse(receipt["release_promotion"])


class ReleaseRACEnvironmentTests(unittest.TestCase):
    def test_changed_bytes_and_unrecorded_files_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prefix = root / "environment"
            prefix.mkdir()
            source = prefix / "module.py"
            source.write_bytes(b"value = 1\n")
            manifest = root / "manifest.json"
            entries = [
                {"path": ".", "type": "directory",
                 "mode": format(stat.S_IMODE(prefix.stat().st_mode), "04o")},
                {"path": "module.py", "type": "file", "link_count": 1,
                 "mode": format(stat.S_IMODE(source.stat().st_mode), "04o"),
                 "size_bytes": source.stat().st_size, "sha256": racs._sha(source)},
            ]
            manifest.write_text(json.dumps({"entries": entries, "entry_count": 2,
                "regular_file_bytes": source.stat().st_size, "tree_merkle_sha256": "fixture"}))
            with patch.object(racs, "TREE_SHA256", racs._sha(manifest)):
                result = racs._verify_environment_tree(prefix, manifest)
                self.assertEqual(result["entry_count"], 2)
                source.write_bytes(b"value = 2\n")
                with self.assertRaisesRegex(racs.ReleaseRACError, "SHA-256 mismatch"):
                    racs._verify_environment_tree(prefix, manifest)
                source.write_bytes(b"value = 1\n")
                (prefix / "unexpected.py").write_text("extra")
                with self.assertRaisesRegex(racs.ReleaseRACError, "added or missing"):
                    racs._verify_environment_tree(prefix, manifest)


if __name__ == "__main__":
    unittest.main()
