from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
import hashlib
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from CoREMOF.cli import build_parser, main
from CoREMOF.parents import LEAKAGE_GUARD_CHOICES, SELECTABLE_PARENT_METHODS


class CliTests(unittest.TestCase):
    def test_benchmark_subcommand_calls_package_api_with_frozen_defaults(self):
        classified = object()
        dataset = SimpleNamespace(classify=lambda checker: classified)

        class Suite:
            dataset_version = "v26.test"
            checker_view = "5checker"
            runs = tuple(range(30))
            fixed_test_ids = ("A", "B")
            group_criteria = ("priority_main",)
            diversity_index_hash = "d" * 64
            assignment_digest = "a" * 64
            official_split = False

            def write(self, directory, stem, overwrite):
                return Path(directory) / stem

        suite = Suite()
        with patch(
            "CoREMOF.dataset.CoREMOFDataset.from_release", return_value=dataset
        ), patch(
            "CoREMOF.benchmarks.build_cr_ncr_benchmark", return_value=suite
        ) as build_mock:
            stdout = StringIO()
            with redirect_stdout(stdout):
                return_code = main(
                    [
                        "benchmark-cr-ncr",
                        "/release",
                        "--output-directory",
                        "/output",
                    ]
                )
        self.assertEqual(return_code, 0)
        self.assertIs(build_mock.call_args.args[0], classified)
        self.assertEqual(
            build_mock.call_args.kwargs["ncr_pool_fractions"],
            (0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
        )
        self.assertEqual(build_mock.call_args.kwargs["seeds"], (42, 43, 44, 45, 46))
        self.assertEqual(build_mock.call_args.kwargs["group_criteria"], ("priority_main",))
        self.assertIsNone(build_mock.call_args.kwargs["cohort_eligibility"])
        self.assertIn('"run_count": 30', stdout.getvalue())
        self.assertIn('"official_split": false', stdout.getvalue())

    def test_benchmark_feasibility_error_is_concise(self):
        from CoREMOF.benchmarks import BenchmarkFeasibilityError

        dataset = SimpleNamespace(classify=lambda checker: object())
        with patch(
            "CoREMOF.dataset.CoREMOFDataset.from_release", return_value=dataset
        ), patch(
            "CoREMOF.benchmarks.build_cr_ncr_benchmark",
            side_effect=BenchmarkFeasibilityError("C=9143 CR and M=9769 NCR, so M>C"),
        ):
            stderr = StringIO()
            with redirect_stderr(stderr):
                return_code = main(
                    [
                        "benchmark-cr-ncr",
                        "/release",
                        "--output-directory",
                        "/output",
                    ]
                )
        self.assertEqual(return_code, 2)
        self.assertIn("C=9143 CR and M=9769 NCR", stderr.getvalue())
        self.assertNotIn("Traceback", stderr.getvalue())

    def test_attach_targets_subcommand_reads_assignment_manifest_and_calls_api(self):
        class Attached:
            missing_policy = "keep"
            structure_ids = ("A", "B")
            dropped_ids = ()
            target_columns = ("uptake",)
            original_assignment_digest = "a" * 64

            def write(self, directory, stem, overwrite):
                root = Path(directory)
                return (
                    root / (stem + ".csv"),
                    root / (stem + ".provenance.jsonl"),
                    root / (stem + ".json"),
                )

        with tempfile.TemporaryDirectory() as temporary_directory:
            manifest = Path(temporary_directory) / "split.csv"
            manifest.write_text(
                "structure_id,partition\nA,train\nB,test\nC,EXCLUDED\n",
                encoding="utf-8",
            )
            digest = hashlib.sha256(
                "A\ttrain\t\nB\ttest\t\nC\tEXCLUDED\t".encode("utf-8")
            ).hexdigest()
            manifest.with_suffix(".json").write_text(
                json.dumps(
                    {
                        "schema_version": "coremof-data-split-receipt/1.0",
                        "assignment_sha256": digest,
                        "partitions": {"train": ["A"], "validation": [], "test": ["B"]},
                        "exclusions": {"C": ""},
                        "release_binding": {
                            "dataset_version": "vtest",
                            "structure_count": 3,
                            "structure_ids_sha256": "f" * 64,
                            "input_hashes": {"metadata/metadata.csv": "b" * 64},
                        },
                    }
                ),
                encoding="utf-8",
            )
            dataset = object()
            with patch(
                "CoREMOF.dataset.CoREMOFDataset.from_release", return_value=dataset
            ), patch(
                "CoREMOF.attachments.attach_targets", return_value=Attached()
            ) as attach_mock:
                with redirect_stdout(StringIO()):
                    return_code = main(
                        [
                            "attach-targets",
                            "/release",
                            "--manifest",
                            str(manifest),
                            "--config",
                            "/targets.json",
                            "--output-directory",
                            "/output",
                        ]
                    )
        self.assertEqual(return_code, 0)
        frozen = attach_mock.call_args.args[0]
        self.assertEqual(dict(frozen.assignments), {"A": "train", "B": "test"})
        self.assertEqual(frozen.assignment_digest, digest)
        self.assertIs(attach_mock.call_args.kwargs["dataset"], dataset)
        self.assertEqual(attach_mock.call_args.kwargs["missing"], "keep")

    def test_attach_targets_rejects_a_manifest_not_bound_to_its_receipt(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            manifest = Path(temporary_directory) / "split.csv"
            manifest.write_text("structure_id,partition\nA,train\n", encoding="utf-8")
            manifest.with_suffix(".json").write_text(
                json.dumps(
                    {
                        "schema_version": "coremof-data-split-receipt/1.0",
                        "assignment_sha256": "0" * 64,
                        "release_binding": {
                            "dataset_version": "vtest",
                            "structure_count": 1,
                            "structure_ids_sha256": "f" * 64,
                            "input_hashes": {"metadata/metadata.csv": "b" * 64},
                        },
                    }
                ),
                encoding="utf-8",
            )
            stderr = StringIO()
            with patch("CoREMOF.attachments.attach_targets") as attach_mock, redirect_stderr(
                stderr
            ):
                return_code = main(
                    [
                        "attach-targets",
                        "/release",
                        "--manifest",
                        str(manifest),
                        "--config",
                        "/targets.json",
                        "--output-directory",
                        "/output",
                    ]
                )
        self.assertEqual(return_code, 2)
        self.assertIn("does not match receipt assignment_sha256", stderr.getvalue())
        attach_mock.assert_not_called()

    def test_attach_targets_help_defines_receipt_and_release_verification(self):
        output = StringIO()
        with redirect_stdout(output), self.assertRaises(SystemExit) as raised:
            main(["attach-targets", "--help"])
        self.assertEqual(raised.exception.code, 0)
        text = " ".join(output.getvalue().split())
        self.assertIn("paired receipt digest", text)
        self.assertIn("release version/universe/input-hash binding", text)
        self.assertIn("default: MANIFEST with a .json suffix", text)

    def test_benchmark_help_defines_fixed_size_and_leakage_terms(self):
        output = StringIO()
        with redirect_stdout(output), self.assertRaises(SystemExit) as raised:
            main(["benchmark-cr-ncr", "--help"])
        self.assertEqual(raised.exception.code, 0)
        text = " ".join(output.getvalue().split()).replace("per- structure", "per-structure")
        self.assertIn("Strict CR has five available PASS votes", text)
        self.assertIn("strict NCR has five available FAIL votes", text)
        self.assertIn("round-half-up(q*M)", text)
        self.assertIn("q=1 uses every eligible NCR row", text)
        self.assertIn("Raw strict-pool and policy-excluded counts", text)
        self.assertIn("fixed_pure_cr", text)
        self.assertIn("main_union plus every ordered selected criterion", text)
        self.assertIn("not identity or parent proof", text)
        self.assertIn("per-structure singleton", text)
        self.assertIn("median/IQR scaling without imputation", text)
        self.assertIn("RT is exact complete 264-value finite RAC5", text)
        self.assertIn("M2T", text)
        self.assertIn("SM", text)
        self.assertIn("complete_release_label_pure_effective_blocks", text)

    def test_attach_targets_help_defines_all_missing_policies(self):
        output = StringIO()
        with redirect_stdout(output), self.assertRaises(SystemExit) as raised:
            main(["attach-targets", "--help"])
        self.assertEqual(raised.exception.code, 0)
        text = " ".join(output.getvalue().split())
        self.assertIn("keep performs a left join and retains nulls", text)
        self.assertIn("error requires every target", text)
        self.assertIn("drop creates only a filtered derived view", text)
        self.assertIn("never refills, rebalances, or resplits", text)

    def test_new_command_defaults_match_the_benchmark_contract(self):
        parser = build_parser()
        args = parser.parse_args(
            [
                "benchmark-cr-ncr",
                "/release",
                "--output-directory",
                "/output",
            ]
        )
        self.assertEqual(args.ncr_pool_fractions, (0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
        self.assertEqual(args.seeds, (42, 43, 44, 45, 46))
        self.assertEqual(args.group_criteria, ("priority_main",))
        self.assertIsNone(args.cohort_eligibility)
        self.assertEqual(args.diversity, "representative")
        self.assertEqual(args.fractions, (0.8, 0.1, 0.1))
        self.assertFalse(args.no_full_cr_diagnostic)

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

    def test_doctor_reports_exact_benchmark_version_drift(self):
        expected = {
            "numpy": "1.26.4",
            "scikit-learn": "9.9.9",
            "scipy": "1.13.1",
            "joblib": "1.5.3",
            "threadpoolctl": "3.6.0",
        }
        output = StringIO()
        with patch(
            "CoREMOF.cli.importlib_metadata.version",
            side_effect=lambda name: expected[name],
        ), redirect_stdout(output):
            return_code = main(["doctor"])
        self.assertEqual(return_code, 0)
        self.assertIn(
            "[MISSING] Representative benchmark: "
            "scikit-learn==1.5.0 (found 9.9.9)",
            output.getvalue(),
        )

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
