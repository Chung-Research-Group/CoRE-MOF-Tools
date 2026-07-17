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
    text = 'a b c d e f g 1 2 3 3\\n'
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


if __name__ == "__main__":
    unittest.main()
