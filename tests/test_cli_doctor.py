import contextlib
import io
import unittest
from unittest.mock import patch

from CoREMOF import cli


class DoctorDependencyTests(unittest.TestCase):
    def test_mofid_requires_its_actual_python_module(self):
        self.assertIn("mofid", cli.FEATURES["MOFid"]["modules"])
        output = io.StringIO()
        with patch.object(cli, "_module_available", side_effect=lambda name: name != "mofid"):
            with contextlib.redirect_stdout(output):
                result = cli.doctor()
        self.assertEqual(result, 0)
        self.assertIn("[MISSING] MOFid: mofid", output.getvalue())
        self.assertNotIn("[OK]      MOFid", output.getvalue())


if __name__ == "__main__":
    unittest.main()
