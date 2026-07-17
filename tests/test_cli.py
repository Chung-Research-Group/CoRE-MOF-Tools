from contextlib import redirect_stdout
from io import StringIO
import unittest

from CoREMOF.cli import main


class CliTests(unittest.TestCase):
    def test_doctor_runs_without_optional_dependencies(self):
        output = StringIO()
        with redirect_stdout(output):
            return_code = main(["doctor"])
        self.assertEqual(return_code, 0)
        self.assertIn("Database lookup", output.getvalue())
        self.assertIn("Missing optional components", output.getvalue())


if __name__ == "__main__":
    unittest.main()
