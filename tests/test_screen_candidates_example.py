import contextlib
import csv
from decimal import Decimal
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
from types import MappingProxyType
import unittest
from unittest.mock import patch

from CoREMOF.dataset import CoREMOFDataset, StructureRecord
from CoREMOF.parents import LEAKAGE_GUARD_CHOICES, SELECTABLE_PARENT_METHODS
from CoREMOF.targets import TargetSource, merge_targets


SCRIPT_PATH = (
    Path(__file__).resolve().parents[1] / "examples" / "screen_candidates.py"
)
SPEC = importlib.util.spec_from_file_location(
    "coremof_screen_candidates_example", SCRIPT_PATH
)
if SPEC is None or SPEC.loader is None:  # pragma: no cover - import machinery guard
    raise RuntimeError("could not load examples/screen_candidates.py")
SCREEN = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = SCREEN
SPEC.loader.exec_module(SCREEN)


IDS = (
    "ASR-COD-2026-0001",  # eligible, lower score
    "ASR-COD-2026-0002",  # eligible, higher score
    "ASR-COD-2026-0003",  # missing required target
    "ASR-COD-2026-0004",  # non-finite ranking value
    "ASR-CSD-2026-0001",  # wrong source
    "FSR-COD-2026-0001",  # wrong variant
    "ASR-COD-2026-0005",  # NCR
    "ASR-COD-2026-0006",  # wrong metal
)


def _write_csv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _base_dataset(root, screening_scores=None):
    records = []
    parents = {}
    for index, structure_id in enumerate(IDS):
        variant, source, _, _ = structure_id.split("-")
        is_ncr = structure_id == "ASR-COD-2026-0005"
        metal = "Zn" if structure_id == "ASR-COD-2026-0006" else "Cu"
        status = "FAIL" if is_ncr else "PASS"
        metadata = {
            "structure_id": structure_id,
            "cif_file": "cifs/{}.cif".format(structure_id),
            "source_database": source,
            "source_id": "SRC-{:02d}".format(index),
            "structure_variant": variant,
            "metal_elements": metal,
            "cell_volume_A3": str(index + 1),
            "mofclassifier_status": status,
            "mofchecker_status": status,
            "chen_manz_status": status,
            "mosaec_status": status,
            "setc_gat_status": status,
            "label_3checker": "NCR" if is_ncr else "CR",
            "label_4checker": "NCR" if is_ncr else "CR",
            "label_5checker": "NCR" if is_ncr else "CR",
        }
        records.append(
            StructureRecord(
                structure_id=structure_id,
                metadata=MappingProxyType(metadata),
                parent_groups=MappingProxyType({}),
                cif_manifest=MappingProxyType(
                    {
                        "structure_id": structure_id,
                        "cif_file": metadata["cif_file"],
                        "size_bytes": "1",
                        "sha256": hashlib.sha256(
                            structure_id.encode("utf-8")
                        ).hexdigest(),
                    }
                ),
            )
        )
        if structure_id in IDS[:2]:
            parents[structure_id] = {
                "rac_status": "MATCHED",
                "rac_group": "R-ABCD0001",
                "rac_size": "2",
            }
        else:
            parents[structure_id] = {
                "rac_status": "UNMATCHED",
                "rac_group": "R-{:08x}".format(index),
                "rac_size": "1",
            }

    dataset = CoREMOFDataset(
        release_root=root,
        records=records,
        dataset_info={"dataset_version": "v-test", "release_status": "FINAL"},
        parent_group_methods={"release_status": "FINAL"},
        parent_by_id=parents,
        input_hashes={"metadata/metadata.csv": "a" * 64},
        cif_files_verified=False,
    )

    scores = {
        IDS[0]: "2.0",
        IDS[1]: "9.0",
        IDS[2]: "100.0",
        IDS[3]: "Infinity",
        IDS[4]: "50.0",
        IDS[5]: "40.0",
        IDS[6]: "30.0",
        IDS[7]: "20.0",
    }
    if screening_scores is not None:
        scores.update(screening_scores)
    _write_csv(
        root / "features" / "rac5_features.csv",
        ("structure_id", "rac5_available", "screening_score"),
        (
            {
                "structure_id": structure_id,
                "rac5_available": "true",
                "screening_score": scores[structure_id],
            }
            for structure_id in IDS
        ),
    )
    _write_csv(
        root / "uptake.csv",
        ("structure_id", "uptake"),
        (
            {
                "structure_id": structure_id,
                "uptake": "" if structure_id == IDS[2] else str(index + 0.5),
            }
            for index, structure_id in enumerate(IDS)
        ),
    )
    return merge_targets(
        dataset,
        (
            TargetSource(
                root / "uptake.csv",
                target_columns=("uptake",),
                value_types={"uptake": "float"},
                units={"uptake": "mol/kg"},
            ),
        ),
        feature_tables=("rac5",),
    )


class ScreeningExampleTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.dataset = _base_dataset(self.root)

    def tearDown(self):
        self.temporary_directory.cleanup()

    def _screen(self, **overrides):
        options = {
            "rank_by": "screening_score",
            "order": "descending",
            "checkers": "5checker",
            "labels": ("CR",),
            "sources": ("COD",),
            "variants": ("ASR",),
            "metals": ("Cu",),
            "required_targets": ("uptake",),
        }
        options.update(overrides)
        return SCREEN.screen_dataset(self.dataset, **options)

    def test_parser_help_defines_parent_leakage_and_auto_terms(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output), self.assertRaises(SystemExit) as raised:
            SCREEN.build_parser().parse_args(("--help",))
        self.assertEqual(raised.exception.code, 0)
        help_text = " ".join(output.getvalue().split())
        self.assertIn("priority_main is the project-defined", help_text)
        self.assertIn("RAC5 then MOFid v2 then MOFid v1", help_text)
        self.assertIn("does not simply choose the first available value on each row", help_text)
        self.assertIn("auto becomes main_union for priority_main", help_text)
        self.assertIn("full-release transitive union", help_text)
        self.assertIn("constructed before filters", help_text)

    def test_parser_exposes_every_defined_parent_and_guard_choice(self):
        parser = SCREEN.build_parser()
        common = (
            "/release",
            "--rank-by",
            "score",
            "--output-directory",
            "/output",
        )
        for method in SELECTABLE_PARENT_METHODS:
            with self.subTest(parent_method=method):
                args = parser.parse_args(common + ("--parent-method", method))
                self.assertEqual(args.parent_method, method)
        for guard in LEAKAGE_GUARD_CHOICES:
            with self.subTest(leakage_guard=guard):
                args = parser.parse_args(common + ("--leakage-guard", guard))
                self.assertEqual(args.leakage_guard, guard)

    def test_filters_before_ranking_and_excludes_nonfinite_fail_closed(self):
        result = self._screen()
        self.assertEqual(result.selected_ids, (IDS[1], IDS[0]))
        self.assertEqual(
            [row["ranking_value"] for row in result.rows],
            [Decimal("9.0"), Decimal("2.0")],
        )
        self.assertEqual(
            result.receipt["counts"],
            {
                "release": 8,
                "after_checker_and_identity_filters": 4,
                "excluded_missing_required_target": 1,
                "after_required_target_filter": 3,
                "excluded_invalid_ranking_value": 1,
                "finite_rankable": 2,
                "emitted": 2,
            },
        )
        self.assertEqual(
            result.receipt["ranking_exclusion_counts"],
            {"NONFINITE_RANK_VALUE": 1},
        )
        self.assertTrue(result.receipt["policies"]["filter_precedes_ranking"])
        self.assertFalse(result.receipt["policies"]["missing_values_imputed"])
        self.assertEqual(
            result.receipt["ranking"]["numeric_policy"],
            {
                "integer": "EXACT",
                "numeric_text": "EXACT_DECIMAL",
                "native_float": "PRESERVE_IEEE_754_VALUE",
                "mixed_comparison": "EXACT_DECIMAL_VALUE",
            },
        )
        source_hashes = result.receipt["implementation"][
            "coremof_source_sha256"
        ]
        self.assertEqual(
            set(source_hashes), {"dataset.py", "labels.py", "targets.py"}
        )
        package_root = Path(__file__).resolve().parents[1] / "CoREMOF"
        for filename, digest in source_hashes.items():
            self.assertEqual(
                digest,
                hashlib.sha256((package_root / filename).read_bytes()).hexdigest(),
            )

    def test_screening_result_is_recursively_immutable_and_revalidated(self):
        result = self._screen()
        with self.assertRaises(TypeError):
            result.rows[0]["structure_id"] = IDS[0]
        with self.assertRaises(TypeError):
            result.receipt["counts"]["emitted"] = 999
        with self.assertRaises(TypeError):
            result.filters["sources"] = ("CSD",)
        with self.assertRaises(TypeError):
            result.receipt["target_merge_receipt"]["sources"][0]["name"] = "changed"

        forged_rows = [dict(row) for row in result.rows]
        forged_rows[0]["structure_id"] = IDS[0]
        forged = SCREEN.ScreeningResult(
            rows=tuple(forged_rows),
            selected_ids=result.selected_ids,
            receipt=SCREEN._json_safe(result.receipt),
            classified_dataset=result.classified_dataset,
            filters=SCREEN._json_safe(result.filters),
        )
        output = self.root / "forged-screening"
        with self.assertRaisesRegex(
            SCREEN.ScreeningError, "rows and selected_ids differ"
        ):
            SCREEN.write_screening(forged, output)
        self.assertFalse(output.exists())

        forged_rows = [dict(row) for row in result.rows]
        forged_rows[0]["classification_label"] = "NCR"
        forged = SCREEN.ScreeningResult(
            rows=tuple(forged_rows),
            selected_ids=result.selected_ids,
            receipt=SCREEN._json_safe(result.receipt),
            classified_dataset=result.classified_dataset,
            filters=SCREEN._json_safe(result.filters),
        )
        with self.assertRaisesRegex(
            SCREEN.ScreeningError, "differs from classified dataset field"
        ):
            SCREEN.write_screening(forged, output)
        self.assertFalse(output.exists())

        receipt_forgeries = (
            (
                "dataset status",
                lambda receipt: receipt.__setitem__("dataset_input_status", "FORGED"),
            ),
            (
                "parent status",
                lambda receipt: receipt.__setitem__("parent_input_status", "FORGED"),
            ),
            (
                "CIF-verification claim",
                lambda receipt: receipt.__setitem__("cif_files_verified", True),
            ),
            (
                "checker-view kind",
                lambda receipt: receipt.__setitem__(
                    "checker_view_kind", "USER_DEFINED"
                ),
            ),
            (
                "ranking contract",
                lambda receipt: receipt["ranking"].__setitem__("kind", "TARGET"),
            ),
            (
                "scientific policies",
                lambda receipt: receipt["policies"].__setitem__(
                    "missing_values_imputed", True
                ),
            ),
        )
        for expected_error, mutate in receipt_forgeries:
            with self.subTest(receipt_claim=expected_error):
                forged_receipt = SCREEN._json_safe(result.receipt)
                mutate(forged_receipt)
                forged = SCREEN.ScreeningResult(
                    rows=tuple(dict(row) for row in result.rows),
                    selected_ids=result.selected_ids,
                    receipt=forged_receipt,
                    classified_dataset=result.classified_dataset,
                    filters=SCREEN._json_safe(result.filters),
                )
                with self.assertRaisesRegex(SCREEN.ScreeningError, expected_error):
                    SCREEN.write_screening(forged, output)
                self.assertFalse(output.exists())

    def test_python_api_rejects_non_integer_and_non_positive_limits(self):
        for invalid in (True, False, 1.2, "2", 0, -1):
            with self.subTest(limit=invalid):
                with self.assertRaisesRegex(
                    SCREEN.ScreeningError, "positive integer"
                ):
                    self._screen(limit=invalid)

    def test_plain_release_receipt_hashes_only_classification_sources(self):
        result = SCREEN.screen_dataset(
            self.dataset.base_dataset,
            rank_by="cell_volume_A3",
            labels=("CR",),
            sources=("COD",),
            variants=("ASR",),
            metals=("Cu",),
        )
        self.assertEqual(
            set(result.receipt["implementation"]["coremof_source_sha256"]),
            {"dataset.py", "labels.py"},
        )
        self.assertIsNone(result.receipt["target_merge_receipt"])
        self.assertEqual(
            result.receipt["implementation"]["source_capture"],
            "MODULE_IMPORT_BOUND",
        )

    def test_import_bound_source_drift_refuses_screening_and_publication(self):
        original_hash = SCREEN._sha256_file

        def drift_script(path):
            if Path(path).resolve() == SCREEN._EXAMPLE_SCRIPT_PATH:
                return "0" * 64
            return original_hash(path)

        with patch.object(SCREEN, "_sha256_file", side_effect=drift_script):
            with self.assertRaisesRegex(
                SCREEN.ScreeningError, "script changed after module import"
            ):
                self._screen()

        result = self._screen()

        forged_receipt = SCREEN._json_safe(result.receipt)
        forged_receipt["implementation"]["coremof_source_sha256"][
            "unexpected.py"
        ] = "0" * 64
        forged = SCREEN.ScreeningResult(
            rows=tuple(dict(row) for row in result.rows),
            selected_ids=result.selected_ids,
            receipt=forged_receipt,
            classified_dataset=result.classified_dataset,
            filters=SCREEN._json_safe(result.filters),
        )
        with self.assertRaisesRegex(SCREEN.ScreeningError, "source closure"):
            SCREEN.write_screening(
                forged,
                self.root / "forged-source-closure",
                stem="must-not-publish",
            )
        self.assertFalse((self.root / "forged-source-closure").exists())

        def drift_dataset(path):
            if Path(path).name == "dataset.py":
                return "f" * 64
            return original_hash(path)

        with patch.object(SCREEN, "_sha256_file", side_effect=drift_dataset):
            with self.assertRaisesRegex(
                SCREEN.ScreeningError,
                "dataset.py changed after module import",
            ):
                SCREEN.write_screening(
                    result,
                    self.root / "drift-output",
                    stem="must-not-publish",
                )
        self.assertFalse((self.root / "drift-output").exists())

    def test_numeric_target_can_rank_ascending_before_deterministic_limit(self):
        result = self._screen(rank_by="uptake", order="ascending", limit=2)
        self.assertEqual(result.selected_ids, (IDS[0], IDS[1]))
        self.assertEqual(
            [row["ranking_value"] for row in result.rows], [0.5, 1.5]
        )
        self.assertEqual(result.receipt["ranking"]["kind"], "TARGET")
        self.assertEqual(result.receipt["counts"]["finite_rankable"], 3)
        self.assertEqual(result.receipt["counts"]["emitted"], 2)

    def test_exact_integer_and_decimal_text_ordering_never_rounds_through_float(self):
        precision_root = self.root / "precision"
        dataset = _base_dataset(
            precision_root,
            {
                IDS[0]: "9007199254740992",
                IDS[1]: "9007199254740993",
            },
        )
        result = SCREEN.screen_dataset(
            dataset,
            rank_by="screening_score",
            order="descending",
            checkers="5checker",
            labels=("CR",),
            sources=("COD",),
            variants=("ASR",),
            metals=("Cu",),
            required_targets=("uptake",),
        )
        self.assertEqual(result.selected_ids, (IDS[1], IDS[0]))
        self.assertEqual(
            [row["ranking_value"] for row in result.rows],
            [Decimal("9007199254740993"), Decimal("9007199254740992")],
        )
        first_paths = SCREEN.write_screening(
            result, precision_root / "first-output", stem="exact-ranking"
        )
        second_paths = SCREEN.write_screening(
            result, precision_root / "second-output", stem="exact-ranking"
        )
        with first_paths["ranked_csv"].open(
            "r", encoding="utf-8", newline=""
        ) as handle:
            written = list(csv.DictReader(handle))
        self.assertEqual(
            [row["ranking_value"] for row in written],
            ["9007199254740993", "9007199254740992"],
        )
        self.assertEqual(
            first_paths["ranked_csv"].read_bytes(),
            second_paths["ranked_csv"].read_bytes(),
        )
        self.assertEqual(
            first_paths["receipt_json"].read_bytes(),
            second_paths["receipt_json"].read_bytes(),
        )

        number, reason = SCREEN._finite_number(9007199254740993)
        self.assertIsNone(reason)
        self.assertEqual(number, 9007199254740993)
        self.assertIsInstance(number, int)

        huge, reason = SCREEN._finite_number("1e999")
        self.assertIsNone(reason)
        self.assertEqual(huge, Decimal("1e999"))

    def test_exact_numeric_ties_use_structure_id_and_invalid_values_fail_closed(self):
        tie_root = self.root / "ties"
        dataset = _base_dataset(
            tie_root,
            {IDS[0]: "1.000000000000000000", IDS[1]: "1.0"},
        )
        result = SCREEN.screen_dataset(
            dataset,
            rank_by="screening_score",
            order="descending",
            labels=("CR",),
            sources=("COD",),
            variants=("ASR",),
            metals=("Cu",),
            required_targets=("uptake",),
        )
        self.assertEqual(result.selected_ids, (IDS[0], IDS[1]))

        for value in (
            True,
            False,
            None,
            "NaN",
            "Infinity",
            Decimal("NaN"),
            Decimal("Infinity"),
            float("nan"),
            float("inf"),
        ):
            with self.subTest(value=value):
                number, reason = SCREEN._finite_number(value)
                self.assertIsNone(number)
                self.assertIn(
                    reason,
                    {
                        "MISSING_RANK_VALUE",
                        "NON_NUMERIC_RANK_VALUE",
                        "NONFINITE_RANK_VALUE",
                    },
                )

    def test_declared_string_na_is_available_but_not_numeric(self):
        token_path = self.root / "tokens.csv"
        _write_csv(
            token_path,
            ("structure_id", "screening_token"),
            (
                {
                    "structure_id": structure_id,
                    "screening_token": (
                        "NA"
                        if structure_id == IDS[0]
                        else "10"
                        if structure_id == IDS[1]
                        else ""
                    ),
                }
                for structure_id in IDS
            ),
        )
        dataset = merge_targets(
            self.dataset.base_dataset,
            (
                TargetSource(
                    token_path,
                    target_columns=("screening_token",),
                    value_types={"screening_token": "string"},
                    units={"screening_token": "user-defined"},
                ),
            ),
            feature_tables=("rac5",),
        )
        common = {
            "checkers": "5checker",
            "labels": ("CR",),
            "sources": ("COD",),
            "variants": ("ASR",),
            "metals": ("Cu",),
            "required_targets": ("screening_token",),
        }
        eligibility = SCREEN.screen_dataset(
            dataset, rank_by="screening_score", **common
        )
        self.assertEqual(eligibility.selected_ids, (IDS[1], IDS[0]))
        self.assertEqual(
            eligibility.receipt["counts"]["excluded_missing_required_target"],
            2,
        )

        numeric = SCREEN.screen_dataset(
            dataset, rank_by="screening_token", **common
        )
        self.assertEqual(numeric.selected_ids, (IDS[1],))
        self.assertEqual(
            numeric.receipt["ranking_exclusion_counts"],
            {"NON_NUMERIC_RANK_VALUE": 1},
        )

    def test_json_required_targets_write_canonical_objects_and_arrays(self):
        json_targets = self.root / "json-targets.csv"
        _write_csv(
            json_targets,
            ("structure_id", "profile", "steps"),
            (
                {
                    "structure_id": structure_id,
                    "profile": '{"b":[2,1],"a":%d}' % (index + 1),
                    "steps": "[1,2]",
                }
                for index, structure_id in enumerate(IDS[:2])
            ),
        )
        dataset = merge_targets(
            self.dataset.base_dataset,
            (
                TargetSource(
                    json_targets,
                    target_columns=("profile", "steps"),
                    value_types={"profile": "json", "steps": "json"},
                ),
            ),
            feature_tables=("rac5",),
        )
        result = SCREEN.screen_dataset(
            dataset,
            rank_by="screening_score",
            labels=("CR",),
            sources=("COD",),
            variants=("ASR",),
            metals=("Cu",),
            required_targets=("profile", "steps"),
        )
        paths = SCREEN.write_screening(
            result, self.root / "json-screening", stem="json_targets"
        )
        with paths["ranked_csv"].open(
            "r", encoding="utf-8", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(rows[0]["structure_id"], IDS[1])
        self.assertEqual(
            rows[0]["required_target:profile"], '{"a":2,"b":[2,1]}'
        )
        self.assertEqual(rows[0]["required_target:steps"], "[1,2]")

    def test_parent_aware_split_and_atomic_outputs_are_auditable(self):
        result = self._screen()
        split = SCREEN.make_parent_aware_split(
            result,
            parent_method="rac5",
            leakage_guard="parent_only",
            fractions=(0.5, 0.0, 0.5),
            random_state=17,
        )
        self.assertEqual(split.assignments[IDS[0]], split.assignments[IDS[1]])
        self.assertTrue(split.leakage_audit["passed"])

        output_directory = self.root / "screening"
        paths = SCREEN.write_screening(
            result,
            output_directory,
            stem="cu_candidates",
            split_result=split,
        )
        self.assertEqual(
            set(paths), {"ranked_csv", "receipt_json", "split_csv", "split_json"}
        )
        with paths["ranked_csv"].open(
            "r", encoding="utf-8", newline=""
        ) as handle:
            reader = csv.DictReader(handle)
            fields = tuple(reader.fieldnames or ())
            rows = list(reader)
        self.assertEqual([row["structure_id"] for row in rows], [IDS[1], IDS[0]])
        self.assertIn("required_target:uptake", fields)
        self.assertNotIn("required_targets", fields)
        self.assertEqual(
            [row["required_target:uptake"] for row in rows], ["1.5", "0.5"]
        )
        receipt_text = paths["receipt_json"].read_text(encoding="utf-8")
        self.assertNotIn("NaN", receipt_text)
        receipt = json.loads(receipt_text)
        self.assertTrue(receipt["split"]["enabled"])
        self.assertTrue(receipt["split"]["leakage_audit"]["passed"])
        self.assertEqual(receipt["split"]["parent_method"], "rac5")
        self.assertEqual(
            receipt["split"]["parent_method_definition"]["criterion"], "rac5"
        )
        self.assertEqual(
            receipt["split"]["requested_leakage_guard"], "parent_only"
        )
        self.assertEqual(receipt["split"]["leakage_guard"], "parent_only")
        self.assertEqual(
            receipt["split"]["requested_leakage_guard_definition"]["identifier"],
            "parent_only",
        )
        self.assertEqual(
            receipt["split"]["leakage_guard_definition"]["identifier"],
            "parent_only",
        )
        self.assertEqual(receipt["outputs"]["ranked_csv"]["row_count"], 2)
        self.assertEqual(
            receipt["csv_output"]["required_target_columns"],
            {"uptake": "required_target:uptake"},
        )
        split_source_hashes = receipt["implementation"][
            "coremof_source_sha256"
        ]
        self.assertEqual(
            set(split_source_hashes),
            {
                "dataset.py",
                "labels.py",
                "parents.py",
                "splitters.py",
                "targets.py",
            },
        )
        self.assertEqual(
            receipt["split"]["implementation_source_sha256"],
            split_source_hashes,
        )
        self.assertEqual(
            receipt["outputs"]["ranked_csv"]["sha256"],
            hashlib.sha256(paths["ranked_csv"].read_bytes()).hexdigest(),
        )
        with self.assertRaises(FileExistsError):
            SCREEN.write_screening(
                result,
                output_directory,
                stem="cu_candidates",
                split_result=split,
            )

    def test_unrelated_split_is_rejected_before_publication(self):
        full_result = self._screen()
        one_row_result = self._screen(limit=1)
        unrelated_split = SCREEN.make_parent_aware_split(
            full_result,
            parent_method="rac5",
            leakage_guard="parent_only",
        )
        output = self.root / "unrelated-split"
        with self.assertRaisesRegex(
            SCREEN.ScreeningError, "structure-ID selection differs"
        ):
            SCREEN.write_screening(
                one_row_result,
                output,
                split_result=unrelated_split,
            )
        self.assertFalse(output.exists())

    def test_nonoverwrite_rollback_preserves_concurrent_replacement_inodes(self):
        staging = self.root / "race-staging"
        output = self.root / "race-output"
        staging.mkdir()
        output.mkdir()
        staged_first = staging / "first.csv"
        staged_second = staging / "second.json"
        final_first = output / staged_first.name
        final_second = output / staged_second.name
        staged_first.write_text("our first\n", encoding="utf-8")
        staged_second.write_text("our second\n", encoding="utf-8")

        real_link = SCREEN.os.link
        call_count = 0

        def replace_between_links(source, destination):
            nonlocal call_count
            call_count += 1
            if call_count == 2:
                final_first.unlink()
                final_first.write_text("other writer first\n", encoding="utf-8")
                final_second.write_text("other writer second\n", encoding="utf-8")
            return real_link(source, destination)

        with patch.object(SCREEN.os, "link", side_effect=replace_between_links):
            with self.assertRaises(FileExistsError):
                SCREEN._publish_bundle(
                    (
                        (staged_first, final_first),
                        (staged_second, final_second),
                    ),
                    overwrite=False,
                )

        self.assertEqual(
            final_first.read_text(encoding="utf-8"), "other writer first\n"
        )
        self.assertEqual(
            final_second.read_text(encoding="utf-8"), "other writer second\n"
        )
        self.assertTrue(staged_first.exists())
        self.assertTrue(staged_second.exists())

    def test_nonoverwrite_rollback_cannot_unlink_post_identity_replacement(self):
        staging = self.root / "post-check-staging"
        output = self.root / "post-check-output"
        staging.mkdir()
        output.mkdir()
        staged_first = staging / "first.csv"
        staged_second = staging / "second.json"
        final_first = output / staged_first.name
        final_second = output / staged_second.name
        staged_first.write_text("our first\n", encoding="utf-8")
        staged_second.write_text("our second\n", encoding="utf-8")

        real_link = SCREEN.os.link
        real_replace = SCREEN.os.replace
        real_samefile = SCREEN.os.path.samefile
        publication_count = 0
        identity_checked = False

        def fail_second_publication(source, destination):
            nonlocal publication_count
            publication_count += 1
            if publication_count == 2:
                raise OSError("deterministic second publication failure")
            return real_link(source, destination)

        def replace_final_after_identity_check(left, right):
            nonlocal identity_checked
            answer = real_samefile(left, right)
            if answer and not identity_checked:
                identity_checked = True
                replacement = output / ".post-check-replacement"
                replacement.write_text("other writer\n", encoding="utf-8")
                real_replace(str(replacement), str(final_first))
            return answer

        with patch.object(
            SCREEN.os, "link", side_effect=fail_second_publication
        ), patch.object(
            SCREEN.os.path,
            "samefile",
            side_effect=replace_final_after_identity_check,
        ):
            with self.assertRaisesRegex(
                OSError, "deterministic second publication failure"
            ):
                SCREEN._publish_bundle(
                    (
                        (staged_first, final_first),
                        (staged_second, final_second),
                    ),
                    overwrite=False,
                )

        self.assertTrue(identity_checked)
        self.assertEqual(final_first.read_text(encoding="utf-8"), "other writer\n")
        self.assertFalse(final_second.exists())

    def test_overwrite_rollback_removes_new_and_restores_existing_output(self):
        staging = self.root / "overwrite-failure-staging"
        output = self.root / "overwrite-failure-output"
        staging.mkdir()
        output.mkdir()
        staged_first = staging / "first.csv"
        staged_second = staging / "second.json"
        final_first = output / staged_first.name
        final_second = output / staged_second.name
        staged_first.write_text("new first\n", encoding="utf-8")
        staged_second.write_text("new second\n", encoding="utf-8")
        final_second.write_text("old second\n", encoding="utf-8")

        real_replace = SCREEN.os.replace
        replacement_count = 0

        def fail_second_replacement(source, destination):
            nonlocal replacement_count
            replacement_count += 1
            if replacement_count == 2:
                raise OSError("deterministic second replacement failure")
            return real_replace(source, destination)

        with patch.object(
            SCREEN.os, "replace", side_effect=fail_second_replacement
        ):
            with self.assertRaisesRegex(
                OSError, "deterministic second replacement failure"
            ):
                SCREEN._publish_bundle(
                    (
                        (staged_first, final_first),
                        (staged_second, final_second),
                    ),
                    overwrite=True,
                )

        self.assertFalse(final_first.exists())
        self.assertEqual(final_second.read_text(encoding="utf-8"), "old second\n")

    def test_overwrite_bundle_is_explicitly_single_writer_and_replaces(self):
        self.assertIn("single-writer", SCREEN._publish_bundle.__doc__)
        self.assertIn("single-writer", SCREEN.write_screening.__doc__)

        result = self._screen()
        output = self.root / "single-writer-overwrite"
        paths = SCREEN.write_screening(result, output, stem="candidate")
        first_bytes = {name: path.read_bytes() for name, path in paths.items()}
        replaced = SCREEN.write_screening(
            result,
            output,
            stem="candidate",
            overwrite=True,
        )
        self.assertEqual(
            first_bytes,
            {name: path.read_bytes() for name, path in replaced.items()},
        )

    def test_screen_receipt_preserves_auto_request_and_main_union_resolution(self):
        result = self._screen()
        split = SCREEN.make_parent_aware_split(
            result,
            parent_method="priority_main",
            leakage_guard="auto",
            fractions=(0.5, 0.0, 0.5),
            random_state=23,
        )
        paths = SCREEN.write_screening(
            result,
            self.root / "auto-screening",
            stem="auto_contract",
            split_result=split,
        )
        receipt = json.loads(paths["receipt_json"].read_text(encoding="utf-8"))
        split_summary = receipt["split"]
        self.assertEqual(split_summary["parent_method"], "priority_main")
        self.assertEqual(
            split_summary["parent_method_definition"]["priority_order"],
            ["rac5", "mofid_v2", "mofid_v1"],
        )
        self.assertEqual(split_summary["requested_leakage_guard"], "auto")
        self.assertEqual(
            split_summary["requested_leakage_guard_definition"]["identifier"],
            "auto",
        )
        self.assertEqual(split_summary["leakage_guard"], "main_union")
        self.assertEqual(
            split_summary["leakage_guard_definition"]["identifier"], "main_union"
        )

    def test_cli_selects_target_config_and_returns_error_when_nothing_rankable(self):
        output_directory = self.root / "cli-output"
        release = self.root / "release"
        config = self.root / "targets.json"
        with patch.object(SCREEN, "load_dataset", return_value=self.dataset) as load:
            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                return_code = SCREEN.main(
                    (
                        str(release),
                        "--target-config",
                        str(config),
                        "--rank-by",
                        "screening_score",
                        "--source",
                        "COD",
                        "--variant",
                        "ASR",
                        "--metal",
                        "Cu",
                        "--require-target",
                        "uptake",
                        "--output-directory",
                        str(output_directory),
                    )
                )
        self.assertEqual(return_code, 0)
        load.assert_called_once_with(
            release, target_config=config, verify_cif_files=False
        )
        summary = json.loads(stdout.getvalue())
        self.assertEqual(summary["checker_view"], "5checker")
        self.assertEqual(summary["counts"]["emitted"], 2)

        with self.assertRaisesRegex(SCREEN.ScreeningError, "no eligible structure"):
            self._screen(rank_by="screening_score", metals=("Ag",))


if __name__ == "__main__":
    unittest.main()
