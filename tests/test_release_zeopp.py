"""Frozen Zeo++ profile and isolated worker tests without the external binary."""
import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

from CoREMOF import release_zeopp as zeo

SID = "FSR-COD-2016-0106"
SHA = "a" * 64


def probe_record():
    return {"structure_id": SID, "input": {"cif_sha256": SHA},
            "protocol_id": zeo.probe.PROTOCOL_ID, "execution_status": "SUCCESS",
            "features": {"intrinsic_props": dict.fromkeys(zeo.probe.INTRINSIC_FEATURES, 0),
                         "N2_probe_props": dict.fromkeys(zeo.probe.N2_FEATURES, 0),
                         "He_probe_props": {"AV_VF": 0}},
            "N2_channel_topology": {"channel_count": 0, "channel_dimensions": [], "maximum_channel_dimension": 0}}


class ReleaseZeoppProfileTests(unittest.TestCase):
    def test_standard_library_import(self):
        result = subprocess.run([sys.executable, "-B", "-S", "-c",
            "from CoREMOF import release_zeopp; import sys; "
            "assert not ({'numpy','ase','pymatgen'} & set(sys.modules))"],
            cwd=Path(__file__).resolve().parents[1], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_recorded_science_sources_and_policy_are_byte_identical(self):
        zeo._check(zeo.probe.__file__, zeo.PROBE_SHA256)
        zeo._check(zeo.framework.__file__, zeo.FRAMEWORK_SHA256)
        zeo._check(Path(zeo.__file__).with_name("data") / "zeopp_probe_namespace_v2.json", zeo.POLICY_SHA256)
        self.assertEqual(zeo.probe.N2_RADIUS_A, 1.655)
        self.assertEqual(zeo.probe.HE_RADIUS_A, 1.32)
        self.assertEqual(zeo.probe.SURFACE_SAMPLES_PER_ATOM, 5000)
        self.assertEqual(zeo.probe.VOLUME_SAMPLES_TOTAL, 5000)

    def test_wrong_binary_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            network = Path(directory) / "network"
            network.write_text("not the recorded binary")
            with self.assertRaisesRegex(zeo.ReleaseZeoppError, "SHA-256"):
                zeo._profile(network)

    def test_complete_zero_features_are_retained(self):
        raw = probe_record()
        zeo._validate(raw, SID, SHA, "n2_he")
        self.assertEqual(raw["features"]["He_probe_props"]["AV_VF"], 0)

    def test_partial_nonfinite_and_negative_values_rejected(self):
        for value in (float("nan"), float("inf"), -1, True):
            raw = probe_record()
            raw["features"]["intrinsic_props"]["LCD_A"] = value
            with self.subTest(value=value), self.assertRaises(zeo.ReleaseZeoppError):
                zeo._validate(raw, SID, SHA, "n2_he")
        raw = probe_record()
        raw["features"]["intrinsic_props"].pop("LCD_A")
        with self.assertRaises(zeo.ReleaseZeoppError):
            zeo._validate(raw, SID, SHA, "n2_he")

    def test_fraction_channel_and_identity_mismatch_rejected(self):
        raw = probe_record()
        raw["features"]["He_probe_props"]["AV_VF"] = 1.1
        with self.assertRaises(zeo.ReleaseZeoppError):
            zeo._validate(raw, SID, SHA, "n2_he")
        raw = probe_record()
        raw["N2_channel_topology"]["channel_count"] = 2
        with self.assertRaises(zeo.ReleaseZeoppError):
            zeo._validate(raw, SID, SHA, "n2_he")
        with self.assertRaises(zeo.ReleaseZeoppError):
            zeo._validate(probe_record(), "ASR-COD-2000-0001", SHA, "n2_he")

    def test_failed_record_cannot_keep_features(self):
        raw = probe_record()
        raw["execution_status"] = "ERROR"
        with self.assertRaises(zeo.ReleaseZeoppError):
            zeo._validate(raw, SID, SHA, "n2_he")
        raw["features"] = None
        raw["N2_channel_topology"] = None
        zeo._validate(raw, SID, SHA, "n2_he")

    def test_framework_matches_full_raw_component_counts(self):
        line = ("MOF C2 3 segments: 2 framework(s) (1D/2D/3D 1 0 1 ) and 1 molecule(s). "
                "Identified dimensionality of framework(s): 1 3")
        raw = {"structure_id": SID, "input": {"cif_sha256": SHA},
               "protocol_id": zeo.framework.PROTOCOL_ID, "execution_status": "SUCCESS",
               "framework_dimension": zeo.framework.parse_strinfo(line)}
        zeo._validate(raw, SID, SHA, "framework")
        raw["framework_dimension"]["molecule_count"] = 0
        with self.assertRaises(zeo.ReleaseZeoppError):
            zeo._validate(raw, SID, SHA, "framework")


class ReleaseZeoppIsolationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.cif = self.root / "input.cif"
        self.cif.write_bytes(b"data_fixture\n")
        self.options = dict(cif_path=self.cif, structure_id=SID, output_dir=self.root / "out",
                            network=sys.executable, timeout_seconds=1)

    def test_existing_output_is_untouched(self):
        self.options["output_dir"].mkdir()
        with patch.object(zeo, "_profile") as profile, self.assertRaises(FileExistsError):
            zeo.calculate_release_zeopp(**self.options)
        profile.assert_not_called()

    def test_invalid_options_never_execute(self):
        for change in ({"structure_id": "../../x"}, {"timeout_seconds": 0},
                       {"timeout_seconds": True}, {"timeout_seconds": 2.5}):
            with self.subTest(change=change), patch.object(zeo.subprocess, "Popen") as launch, \
                 self.assertRaises(ValueError):
                zeo.calculate_release_zeopp(**{**self.options, **change})
            launch.assert_not_called()

    def _patched_call(self, launch):
        with patch.object(zeo, "_profile", return_value=self.cif), \
             patch.object(zeo, "_check", side_effect=lambda path, digest: Path(path)), \
             patch.object(zeo.subprocess, "Popen", side_effect=launch):
            return zeo.calculate_release_zeopp(**self.options)

    def test_isolated_manifest_clean_environment_and_error_outputs(self):
        def launch(command, **kwargs):
            import csv
            path = Path(command[command.index("--manifest") + 1])
            with path.open() as stream:
                row = next(csv.DictReader(stream))
            copied = Path(row["canonical_cif_path"])
            self.assertNotEqual(copied, self.cif)
            self.assertEqual(copied.read_bytes(), self.cif.read_bytes())
            self.assertIn("-S", command)
            self.assertTrue(kwargs["start_new_session"])
            self.assertNotIn("LD_PRELOAD", kwargs["env"])
            return types.SimpleNamespace(returncode=2, wait=lambda timeout: None)
        result = self._patched_call(launch)
        self.assertEqual(result["execution_status"], "ERROR")
        self.assertIsNone(result["features"])
        self.assertIsNone(result["framework_dimension"])
        self.assertEqual(self.cif.read_bytes(), b"data_fixture\n")
        self.assertFalse(json.loads((self.options["output_dir"] / "receipt.json").read_text())["release_metadata_promoted"])

    def test_outer_timeout_terminates_only_owned_groups(self):
        fake = types.SimpleNamespace(returncode=-15,
            wait=lambda timeout: (_ for _ in ()).throw(subprocess.TimeoutExpired("worker", 1)))
        with patch.object(zeo, "_stop") as stop:
            result = self._patched_call(lambda *args, **kwargs: fake)
        self.assertEqual(stop.call_count, 2)
        self.assertTrue(all(x["error_type"] == "TIMEOUT" for x in result["components"].values()))

    def test_component_failure_does_not_discard_independent_success(self):
        def launch(command, **kwargs):
            if "--namespace-policy" not in command:
                return types.SimpleNamespace(returncode=2, wait=lambda timeout: None)
            raw = copy.deepcopy(probe_record())
            raw["input"]["cif_sha256"] = zeo._sha(self.cif)
            output = Path(command[command.index("--output-root") + 1]) / "records"
            output.mkdir()
            (output / (SID + ".json")).write_text(json.dumps(raw))
            return types.SimpleNamespace(returncode=0, wait=lambda timeout: None)
        result = self._patched_call(launch)
        self.assertEqual(result["execution_status"], "PARTIAL")
        self.assertIsNotNone(result["features"])
        self.assertIsNone(result["framework_dimension"])


if __name__ == "__main__":
    unittest.main()
