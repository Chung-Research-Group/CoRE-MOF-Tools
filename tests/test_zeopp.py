import os
from pathlib import Path
import stat
import tempfile
import unittest
from unittest.mock import patch

from CoREMOF.calculation import Zeopp


FAKE_NETWORK = """#!/usr/bin/env python3
import pathlib
import sys

args = sys.argv[1:]
output = pathlib.Path(args[-2])
if '-fail' in args:
    print('intentional failure', file=sys.stderr)
    raise SystemExit(7)
if '-chan' in args:
    text = 'Channel dimensionality 3\\n'
elif '-strinfo' in args:
    if args[-2] != '-strinfo':
        raise SystemExit('strinfo accepts only the CIF, not an output argument')
    output = pathlib.Path(args[-1]).with_suffix('.strinfo')
    text = ('MOF formula 6 segments: 6 framework(s) (1D/2D/3D 1 2 3 ) '
            'and 0 molecule(s). Identified dimensionality of framework(s): 1 2 2 3 3 3\\n')
elif '-res' in args:
    text = 'MOF 12.5 4.25 8.0\\n'
elif '-sa' in args:
    text = ('ASA_A^2: 1 ASA_m^2/cm^3: 2 ASA_m^2/g: 3 '
            'NASA_A^2: 4 NASA_m^2/cm^3: 5 NASA_m^2/g: 6\\n')
elif '-volpo' in args:
    text = ('POAV_A^3: 1 PONAV_A^3: 2 POAV_cm^3/g: 3 PONAV_cm^3/g: 4 '
            'POAV_Volume_fraction: 0.25 PONAV_Volume_fraction: 0.10\\n')
else:
    raise SystemExit(2)
output.write_text(text)
"""


class ZeoppTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.executable = self.root / "fake network"
        self.executable.write_text(FAKE_NETWORK, encoding="utf-8")
        self.executable.chmod(self.executable.stat().st_mode | stat.S_IXUSR)
        self.structure = self.root / "structure with spaces.cif"
        self.structure.write_text("data_test\n", encoding="utf-8")
        self.environment = patch.dict(
            os.environ, {"COREMOF_NETWORK_EXECUTABLE": str(self.executable)}
        )
        self.environment.start()

    def tearDown(self):
        self.environment.stop()
        self.temporary_directory.cleanup()

    def test_all_parsers_and_paths_with_spaces(self):
        prefix = str(self.root / "parallel run")
        self.assertEqual(Zeopp.ChanDim(self.structure, prefix=prefix)["Dimension"], 3)
        self.assertEqual(Zeopp.FrameworkDim(self.structure, prefix=prefix)["N_2D"], 2)
        self.assertEqual(Zeopp.PoreDiameter(self.structure, prefix=prefix)["PLD"], 4.25)
        self.assertEqual(Zeopp.SurfaceArea(self.structure, prefix=prefix)["ASA"], [1, 2, 3])
        self.assertEqual(Zeopp.PoreVolume(self.structure, prefix=prefix)["NVF"], 0.10)
        self.assertEqual(list(self.root.glob("parallel run_*.txt")), [])

    def test_missing_input_is_reported_before_execution(self):
        with self.assertRaisesRegex(FileNotFoundError, "CIF file does not exist"):
            Zeopp.PoreDiameter(self.root / "missing.cif")

    def test_nonzero_exit_includes_zeopp_diagnostic(self):
        with self.assertRaisesRegex(RuntimeError, "intentional failure"):
            Zeopp._run_network(self.structure, ["-fail"], str(self.root / "failure"))

    def test_nonfinite_and_negative_output_is_unavailable_not_a_value(self):
        surface_labels = (
            "ASA_A^2:", "ASA_m^2/cm^3:", "ASA_m^2/g:",
            "NASA_A^2:", "NASA_m^2/cm^3:", "NASA_m^2/g:",
        )
        volume_labels = (
            "POAV_A^3:", "PONAV_A^3:", "POAV_cm^3/g:", "PONAV_cm^3/g:",
            "POAV_Volume_fraction:", "PONAV_Volume_fraction:",
        )
        for parser, labels in ((Zeopp.SurfaceArea, surface_labels), (Zeopp.PoreVolume, volume_labels)):
            for label in labels:
                for bad in ("nan", "inf", "-Infinity", "-1"):
                    with self.subTest(parser=parser.__name__, label=label, bad=bad):
                        output = " ".join(f"{key} {bad if key == label else '0'}" for key in labels)
                        with patch.object(Zeopp, "_run_network", return_value=output):
                            with self.assertRaises(ValueError):
                                parser(self.structure)
        for index in range(3):
            for bad in ("nan", "inf", "-Infinity", "-1"):
                fields = ["0", "0", "0"]
                fields[index] = bad
                with self.subTest(diameter=index, bad=bad):
                    with patch.object(Zeopp, "_run_network", return_value="MOF " + " ".join(fields)):
                        with self.assertRaises(ValueError):
                            Zeopp.PoreDiameter(self.structure)

    def test_zero_is_preserved(self):
        with patch.object(Zeopp, "_run_network", return_value="MOF 0 0 0"):
            self.assertEqual(Zeopp.PoreDiameter(self.structure)["PLD"], 0.0)
        with patch.object(Zeopp, "_run_network", return_value="Channel dimensionality 0"):
            self.assertEqual(Zeopp.ChanDim(self.structure, probe_radius=0)["Dimension"], 0)
        self.assertEqual(Zeopp._labelled_float("ASA_A^2: 0", "ASA_A^2:"), 0.0)

    def test_invalid_dimensions_and_counts_are_rejected(self):
        for value in (-1, 4):
            with patch.object(Zeopp, "_run_network", return_value=f"Channel dimensionality {value}"):
                with self.assertRaises(ValueError):
                    Zeopp.ChanDim(self.structure)
        for counts in ("1 2 3 4", "1 2 3 -1", "-1 2 3 3", "1 -2 3 3", "1 2 -3 3"):
            with patch.object(Zeopp, "_run_network", return_value="a b c d e f g " + counts):
                with self.assertRaises(ValueError):
                    Zeopp.FrameworkDim(self.structure)

    def test_invalid_options_never_launch_network(self):
        with patch.object(Zeopp, "_run_network") as run:
            for radius in (-1, float("nan"), float("inf"), True, "1.655"):
                for parser, kwargs in (
                    (Zeopp.ChanDim, {"probe_radius": radius}),
                    (Zeopp.SurfaceArea, {"chan_radius": radius}),
                    (Zeopp.SurfaceArea, {"probe_radius": radius}),
                    (Zeopp.PoreVolume, {"chan_radius": radius}),
                    (Zeopp.PoreVolume, {"probe_radius": radius}),
                ):
                    with self.subTest(parser=parser.__name__, kwargs=kwargs):
                        with self.assertRaises(ValueError):
                            parser(self.structure, **kwargs)
            for count in (None, 0, -1, True, 2.5, "5000"):
                for parser in (Zeopp.SurfaceArea, Zeopp.PoreVolume):
                    with self.assertRaises(ValueError):
                        parser(self.structure, num_samples=count)
            with self.assertRaises(ValueError):
                Zeopp.PoreDiameter(self.structure, high_accuracy="False")
            run.assert_not_called()

    def test_all_channel_dimensions_and_zero_channels(self):
        for text, expected in (("0 channels identified of dimensionality", 0),
                               ("3 channels identified of dimensionality 1 3 2", 3)):
            with patch.object(Zeopp, "_run_network", return_value=text):
                self.assertEqual(Zeopp.ChanDim(self.structure)["Dimension"], expected)
        for text in ("0 channels identified of dimensionality 1", "2 channels identified of dimensionality 3",
                     "1 channels identified of dimensionality 0", "1 channels identified of dimensionality 4"):
            with patch.object(Zeopp, "_run_network", return_value=text), self.assertRaises(ValueError):
                Zeopp.ChanDim(self.structure)

    def test_discrete_molecules_are_valid_zero_framework_dimension(self):
        text = ("MOF C2 2 segments: 0 framework(s) (1D/2D/3D 0 0 0 ) and 2 molecule(s). "
                "Identified dimensionality of framework(s):")
        with patch.object(Zeopp, "_run_network", return_value=text):
            self.assertEqual(Zeopp.FrameworkDim(self.structure),
                             {"unit": "nan", "Dimension": 0, "N_1D": 0, "N_2D": 0, "N_3D": 0})

    def test_surface_labels_cannot_match_inside_nasa(self):
        self.assertEqual(Zeopp._labelled_float("NASA_A^2: 99 ASA_A^2: 2", "ASA_A^2:"), 2)
        for text in ("NASA_A^2: 99", "ASA_A^2: 1 ASA_A^2: 2"):
            with self.assertRaises(ValueError):
                Zeopp._labelled_float(text, "ASA_A^2:")

    def test_source_cif_and_siblings_are_not_modified(self):
        before = self.structure.read_bytes()
        sibling = self.structure.with_suffix(".strinfo")
        sibling.write_text("keep original")
        Zeopp.FrameworkDim(self.structure, prefix=str(self.root / "isolated"))
        self.assertEqual(self.structure.read_bytes(), before)
        self.assertEqual(sibling.read_text(), "keep original")
        self.assertEqual(list(self.root.glob("isolated_*")), [])

    def test_void_fraction_larger_than_one_is_rejected(self):
        text = ("POAV_A^3: 1 PONAV_A^3: 2 POAV_cm^3/g: 3 PONAV_cm^3/g: 4 "
                "POAV_Volume_fraction: 1.01 PONAV_Volume_fraction: 0.1")
        with patch.object(Zeopp, "_run_network", return_value=text), self.assertRaises(ValueError):
            Zeopp.PoreVolume(self.structure)


if __name__ == "__main__":
    unittest.main()
