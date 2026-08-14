from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from CoREMOF.cli import build_parser, main
from CoREMOF.parents import LEAKAGE_GUARD_CHOICES, SELECTABLE_PARENT_METHODS


class CliTests(unittest.TestCase):
    def test_split_help_defines_project_parent_and_leakage_terms(self):
        output = StringIO()
        with redirect_stdout(output), self.assertRaises(SystemExit) as raised:
            main(["split", "--help"])
        self.assertEqual(raised.exception.code, 0)
        help_text = " ".join(output.getvalue().split()).replace(
            "project- defined", "project-defined"
        )
        self.assertIn("priority_main is the project-defined", help_text)
        self.assertIn("RAC5 then MOFid v2 then MOFid v1", help_text)
        self.assertIn("not a row-wise first-nonmissing fallback", help_text)
        self.assertIn("priority means parent-evidence precedence", help_text)
        self.assertIn("does not rank, schedule, or recalculate", help_text)
        self.assertIn("identity_union selects the separate project-defined", help_text)
        self.assertIn("source-ID/MOFid transitive groups", help_text)
        self.assertIn("identity_size count one transitive connected component", help_text)
        self.assertIn("Unicode NFKC and casefold", help_text)
        self.assertIn("each named release is freshly recomputed", help_text)
        self.assertIn("without importing an earlier component", help_text)
        self.assertIn("unsuccessful MOFid statuses add no edge", help_text)
        self.assertIn("SUCCESS_TOPOLOGY_UNKNOWN", help_text)
        self.assertIn("rac5_crystalnets (RT-)", help_text)
        self.assertIn("mofid_v2_crystalnets (M2T-)", help_text)
        self.assertIn("release-authorized MOFid-v2 values change", help_text)
        self.assertIn("rebuild M2T before use", help_text)
        self.assertNotIn("rac5_topology", help_text)
        self.assertNotIn("mofid_v2_topology", help_text)
        self.assertIn("structure_matcher_strict (SM-)", help_text)
        self.assertIn("main_union is the project-defined leakage guard", help_text)
        self.assertIn("parent_only is the project-defined guard", help_text)
        self.assertIn("auto is only a project-defined selector", help_text)
        self.assertIn("main_union for priority_main", help_text)
        self.assertIn("transitive union over the complete release", help_text)
        self.assertIn("before filters", help_text)
        self.assertIn("missing CIF hashes fail", help_text)
        self.assertIn("adds no cross-method edge", help_text)
        self.assertIn("MATCHED means available size>=2", help_text)
        self.assertIn("NOT_AVAILABLE supplies no edge", help_text)
        self.assertIn("all 264 finite binary64 values", help_text)
        self.assertIn("float.hex with rtol=atol=0", help_text)
        self.assertIn("not proof of structural identity", help_text)
        self.assertIn("(V/Nsites)^(1/3) is dimensionless", help_text)
        self.assertIn("manual compatibility object", help_text)
        for choice in SELECTABLE_PARENT_METHODS:
            self.assertIn(choice, help_text)
        for guard in LEAKAGE_GUARD_CHOICES:
            self.assertIn(guard, help_text)

    def test_split_parser_exposes_every_defined_parent_and_guard_choice(self):
        parser = build_parser()
        for method in SELECTABLE_PARENT_METHODS:
            with self.subTest(parent_method=method):
                args = parser.parse_args(
                    [
                        "split",
                        "/release",
                        "--parent-method",
                        method,
                        "--output-directory",
                        "/output",
                    ]
                )
                self.assertEqual(args.parent_method, method)
        for guard in LEAKAGE_GUARD_CHOICES:
            with self.subTest(leakage_guard=guard):
                args = parser.parse_args(
                    [
                        "split",
                        "/release",
                        "--leakage-guard",
                        guard,
                        "--output-directory",
                        "/output",
                    ]
                )
                self.assertEqual(args.leakage_guard, guard)

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
                requested_leakage_guard="auto",
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
            self.assertIsNone(call.kwargs["required_targets"])
            self.assertEqual(call.kwargs["required_target_mode"], "all")
            self.assertIn('"provisional_input": true', output.getvalue())
            self.assertIn('"cif_files_verified": false', output.getvalue())
            self.assertIn('"requested_leakage_guard": "auto"', output.getvalue())

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

    def test_split_cli_accepts_optional_reference_parent_methods(self):
        result = SimpleNamespace(
            dataset_version="v26.test",
            checker_view="5checker",
            parent_method="structure_matcher_strict",
            leakage_guard="parent_only",
            provisional_input=True,
            cif_files_verified=False,
            counts={"train": 2, "validation": 1, "test": 1, "excluded": 0},
            leakage_audit={"passed": True, "cross_split_block_count": 0},
            write=lambda directory, stem, overwrite: (
                Path(directory) / (stem + ".csv"),
                Path(directory) / (stem + ".json"),
            ),
        )
        with patch("CoREMOF.splitters.split_release", return_value=result) as mocked:
            with redirect_stdout(StringIO()):
                return_code = main(
                    [
                        "split",
                        "/release",
                        "--parent-method",
                        "structure_matcher_strict",
                        "--leakage-guard",
                        "parent_only",
                        "--output-directory",
                        "/output",
                    ]
                )
        self.assertEqual(return_code, 0)
        self.assertEqual(
            mocked.call_args.kwargs["parent_method"], "structure_matcher_strict"
        )
        self.assertEqual(mocked.call_args.kwargs["leakage_guard"], "parent_only")

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

    def test_target_config_is_merged_before_required_target_split(self):
        target_dataset = object()
        result = SimpleNamespace(
            dataset_version="v26.test",
            checker_view="5checker",
            parent_method="priority_main",
            leakage_guard="main_union",
            provisional_input=True,
            cif_files_verified=False,
            counts={"train": 2, "validation": 0, "test": 0, "excluded": 2},
            leakage_audit={"passed": True, "cross_split_block_count": 0},
            write=lambda directory, stem, overwrite: (
                Path(directory) / (stem + ".csv"),
                Path(directory) / (stem + ".json"),
            ),
        )
        with patch(
            "CoREMOF.targets.merge_targets_from_config",
            return_value=target_dataset,
        ) as merge_mock, patch(
            "CoREMOF.splitters.split_release", return_value=result
        ) as split_mock:
            with redirect_stdout(StringIO()):
                return_code = main(
                    [
                        "split",
                        "/release",
                        "--target-config",
                        "/inputs/targets.json",
                        "--require-target",
                        "uptake",
                        "--require-target",
                        "selectivity",
                        "--required-target-mode",
                        "all",
                        "--output-directory",
                        "/output",
                    ]
                )
        self.assertEqual(return_code, 0)
        self.assertIs(split_mock.call_args.args[0], target_dataset)
        self.assertEqual(
            split_mock.call_args.kwargs["required_targets"],
            ("uptake", "selectivity"),
        )
        self.assertEqual(split_mock.call_args.kwargs["required_target_mode"], "all")
        merge_mock.assert_called_once()

    def test_required_target_without_config_fails_cleanly(self):
        stderr = StringIO()
        with redirect_stderr(stderr):
            return_code = main(
                [
                    "split",
                    "/release",
                    "--require-target",
                    "uptake",
                    "--output-directory",
                    "/output",
                ]
            )
        self.assertEqual(return_code, 2)
        self.assertIn("requires --target-config", stderr.getvalue())

    def test_merge_targets_subcommand_writes_three_outputs(self):
        class Merged:
            dataset_version = "v26.test"
            target_columns = ("uptake", "selectivity")
            feature_columns = ("rac_a",)

            def __len__(self):
                return 4

            def receipt(self):
                return {"target_values_sha256": "a" * 64}

            def write(self, directory, stem, overwrite):
                return (
                    Path(directory) / (stem + ".csv"),
                    Path(directory) / (stem + ".provenance.jsonl"),
                    Path(directory) / (stem + ".json"),
                )

        merged = Merged()
        with patch(
            "CoREMOF.targets.merge_targets_from_config", return_value=merged
        ):
            stdout = StringIO()
            with redirect_stdout(stdout):
                return_code = main(
                    [
                        "merge-targets",
                        "/release",
                        "--config",
                        "/inputs/targets.json",
                        "--output-directory",
                        "/output",
                    ]
                )
        self.assertEqual(return_code, 0)
        summary = stdout.getvalue()
        self.assertIn('"target_columns"', summary)
        self.assertIn('"provenance_jsonl"', summary)


if __name__ == "__main__":
    unittest.main()
