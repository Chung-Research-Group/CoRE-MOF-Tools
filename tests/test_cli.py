from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from CoREMOF.cli import main


class CliTests(unittest.TestCase):
    def test_doctor_runs_without_optional_dependencies(self):
        output = StringIO()
        with redirect_stdout(output):
            return_code = main(["doctor"])
        self.assertEqual(return_code, 0)
        self.assertIn("Database lookup", output.getvalue())
        self.assertIn("Missing optional components", output.getvalue())

    def test_split_subcommand_uses_safe_scientific_defaults(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_directory = Path(temporary_directory) / "split"
            result = SimpleNamespace(
                dataset_version="v26.test",
                checker_view="5checker",
                parent_method="priority_main",
                leakage_guard="main_union",
                provisional_input=True,
                cif_files_verified=False,
                counts={"train": 8, "validation": 1, "test": 1, "excluded": 2},
                leakage_audit={"passed": True, "cross_split_block_count": 0},
                write=lambda directory, stem, overwrite: (
                    Path(directory) / (stem + ".csv"),
                    Path(directory) / (stem + ".json"),
                ),
            )
            with patch("CoREMOF.splitters.split_release", return_value=result) as mocked:
                output = StringIO()
                with redirect_stdout(output):
                    return_code = main(
                        [
                            "split",
                            "/release",
                            "--output-directory",
                            str(output_directory),
                        ]
                    )
            self.assertEqual(return_code, 0)
            self.assertTrue(mocked.called)
            call = mocked.call_args
            self.assertEqual(call.kwargs["checkers"], "5checker")
            self.assertEqual(call.kwargs["parent_method"], "priority_main")
            self.assertEqual(call.kwargs["leakage_guard"], "auto")
            self.assertEqual(call.kwargs["missing_parent"], "singleton")
            self.assertEqual(call.kwargs["labels"], ("CR", "NCR"))
            self.assertFalse(call.kwargs["verify_cif_files"])
            self.assertFalse(call.kwargs["official"])
            self.assertIn('"provisional_input": true', output.getvalue())
            self.assertIn('"cif_files_verified": false', output.getvalue())

    def test_official_cli_request_fails_cleanly(self):
        from CoREMOF.splitters import OfficialSplitUnavailableError

        error = OfficialSplitUnavailableError("no audited manifest")
        with patch("CoREMOF.splitters.split_release", side_effect=error):
            stderr = StringIO()
            with redirect_stderr(stderr):
                return_code = main(
                    [
                        "split",
                        "/release",
                        "--official",
                        "--output-directory",
                        "/output",
                    ]
                )
        self.assertEqual(return_code, 2)
        self.assertIn("no audited manifest", stderr.getvalue())

    def test_expected_split_validation_error_fails_without_traceback(self):
        with patch(
            "CoREMOF.splitters.split_release",
            side_effect=ValueError("bad release evidence"),
        ):
            stderr = StringIO()
            with redirect_stderr(stderr):
                return_code = main(
                    ["split", "/release", "--output-directory", "/output"]
                )
        self.assertEqual(return_code, 2)
        self.assertIn("bad release evidence", stderr.getvalue())
        self.assertNotIn("Traceback", stderr.getvalue())

    def test_unexpected_programming_error_is_not_hidden(self):
        with patch(
            "CoREMOF.splitters.split_release",
            side_effect=RuntimeError("unexpected internal failure"),
        ):
            with self.assertRaisesRegex(RuntimeError, "unexpected internal failure"):
                main(["split", "/release", "--output-directory", "/output"])


if __name__ == "__main__":
    unittest.main()
