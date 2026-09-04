import csv
from collections import Counter
import ctypes
from dataclasses import replace
import errno
import hashlib
import inspect
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import MappingProxyType
import unittest
from unittest.mock import patch

import CoREMOF.dataset as dataset_module
import CoREMOF.attachments as attachments_module
import CoREMOF.benchmarks as benchmarks_module

from CoREMOF import (
    CRNCRBenchmarkSuite as PublicCRNCRBenchmarkSuite,
    DataSplitResult as PublicDataSplitResult,
    TargetAttachedView as PublicTargetAttachedView,
    attach_targets,
)
from CoREMOF.attachments import (
    FrozenAssignmentManifest,
    TargetAttachedBenchmarkSuite,
    TargetAttachmentError,
    frozen_assignment_manifest,
)
from CoREMOF.benchmarks import (
    BenchmarkDependencyError,
    BenchmarkError,
    BenchmarkFeasibilityError,
    BENCHMARK_JOBLIB_VERSION,
    BENCHMARK_NUMPY_VERSION,
    BENCHMARK_SCIPY_VERSION,
    BENCHMARK_SKLEARN_VERSION,
    BENCHMARK_THREADPOOLCTL_VERSION,
    LABEL_PURE_BLOCK_ELIGIBILITY,
    available_group_criteria,
    build_cr_ncr_benchmark,
    build_cr_ncr_cohorts,
    build_diversity_index,
    data_split,
    _fixed_test_blocks,
    _cluster_tier,
    _nested_exact_block_memberships,
    normalize_group_criteria,
)
from CoREMOF.dataset import (
    CoREMOFDataset,
    StructureRecord,
)
from CoREMOF.parents import ZEO_NUMERIC_FINGERPRINT_DEFINITION
from CoREMOF.targets import AliasRegistry, TargetSource, merge_targets


try:
    import joblib as _benchmark_joblib
    import numpy as _benchmark_numpy
    import scipy as _benchmark_scipy
    import sklearn as _benchmark_sklearn
    import threadpoolctl as _benchmark_threadpoolctl
except (ImportError, ModuleNotFoundError):
    EXACT_BENCHMARK_BACKEND = False
else:
    EXACT_BENCHMARK_BACKEND = (
        _benchmark_numpy.__version__ == BENCHMARK_NUMPY_VERSION
        and _benchmark_sklearn.__version__ == BENCHMARK_SKLEARN_VERSION
        and _benchmark_scipy.__version__ == BENCHMARK_SCIPY_VERSION
        and _benchmark_joblib.__version__ == BENCHMARK_JOBLIB_VERSION
        and _benchmark_threadpoolctl.__version__ == BENCHMARK_THREADPOOLCTL_VERSION
    )


class _Dataset:
    def __init__(self, ids, parents=None, hashes=None, rows=None):
        if rows is None:
            rows = [
                {
                    "structure_id": structure_id,
                    "source_database": "COD",
                    "source_id": "SOURCE-{}".format(structure_id),
                    "structure_variant": "ASR" if index % 2 == 0 else "FSR",
                    "metal_elements": "Cu" if index % 3 else "Zn",
                }
                for index, structure_id in enumerate(ids)
            ]
        self.metadata_rows = tuple(rows)
        self.parent_by_id = parents or {}
        self.cif_hashes = hashes or {
            structure_id: hashlib.sha256(structure_id.encode("utf-8")).hexdigest()
            for structure_id in ids
        }
        self.input_hashes = {"metadata/metadata.csv": "a" * 64}
        self.dataset_version = "vtest"
        self.parent_group_methods = {}
        self.structure_ids = tuple(row["structure_id"] for row in self.metadata_rows)


class _Classified:
    def __init__(self, dataset, labels, checker_view="5checker", ids=None):
        self.dataset = dataset
        self.structure_ids = tuple(dataset.structure_ids if ids is None else ids)
        self.label_by_id = {
            structure_id: labels[structure_id] for structure_id in self.structure_ids
        }
        self.checker_view = checker_view


def _authenticated_classified(
    root,
    cr=20,
    ncr=8,
    ambiguous=0,
    unchecked=0,
    metadata_hash="b" * 64,
    input_hashes=None,
    cif_hashes=None,
):
    records = []
    for index in range(cr + ncr + ambiguous + unchecked):
        structure_id = "ID-{:03d}".format(index)
        if index < cr:
            statuses = ("PASS",) * 5
        elif index < cr + ncr:
            statuses = ("FAIL",) * 5
        elif index < cr + ncr + ambiguous:
            statuses = ("PASS", "FAIL", "PASS", "FAIL", "PASS")
        else:
            statuses = ("NOT_AVAILABLE", "FAIL", "FAIL", "FAIL", "FAIL")
        metadata = MappingProxyType(
            {
                "structure_id": structure_id,
                "source_database": "COD",
                "source_id": "SOURCE-{}".format(structure_id),
                "structure_variant": "ASR" if index % 2 == 0 else "FSR",
                "metal_elements": "Cu" if index % 3 else "Zn",
                "mofclassifier_status": statuses[0],
                "mofchecker_status": statuses[1],
                "chen_manz_status": statuses[2],
                "mosaec_status": statuses[3],
                "setc_gat_status": statuses[4],
            }
        )
        sha = (cif_hashes or {}).get(
            structure_id,
            hashlib.sha256(structure_id.encode("utf-8")).hexdigest(),
        )
        records.append(
            StructureRecord(
                structure_id=structure_id,
                metadata=metadata,
                parent_groups=MappingProxyType({}),
                cif_manifest=MappingProxyType(
                    {
                        "structure_id": structure_id,
                        "cif_file": "cifs/{}.cif".format(structure_id),
                        "size_bytes": "1",
                        "sha256": sha,
                    }
                ),
            )
        )
    dataset = CoREMOFDataset(
        release_root=root,
        records=records,
        dataset_info={"dataset_version": "vtest"},
        parent_group_methods={},
        parent_by_id={},
        input_hashes=(
            {"metadata/metadata.csv": metadata_hash}
            if input_hashes is None
            else input_hashes
        ),
        cif_files_verified=False,
    )
    dataset_module._register_dataset_generation(
        dataset,
        kind="validated_release",
        official_release_source=True,
    )
    return dataset, dataset.classify("5checker")


def _classified(cr=20, ncr=8, ambiguous=0):
    return _authenticated_classified(
        Path("/synthetic-coremof-release"), cr=cr, ncr=ncr, ambiguous=ambiguous
    )[1]


def _write_csv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _real_dataset(root, cr=10, ncr=4):
    return _authenticated_classified(root, cr=cr, ncr=ncr)


class CriterionAndGeneralSplitTests(unittest.TestCase):
    def test_published_parent_registry_accepts_only_the_exact_release_status(self):
        dataset_module._reject_retired_or_reserved_keys(
            {"release_status": "PROVISIONAL_LATEST_AUDITED_SNAPSHOT"},
            "parent_group_methods",
        )
        with self.assertRaisesRegex(Exception, "retired or reserved"):
            dataset_module._reject_retired_or_reserved_keys(
                {"release_status": "PUBLICATION_AUTHORIZED"},
                "parent_group_methods",
            )

    def test_result_types_are_exported_from_the_package_root(self):
        result = data_split(_classified(), group_criteria="none", diversity="none")
        self.assertIsInstance(result, PublicDataSplitResult)
        suite = build_cr_ncr_benchmark(
            _classified(),
            ncr_pool_fractions=(0, 1),
            seeds=(42,),
            group_criteria="none",
            diversity="none",
        )
        self.assertIsInstance(suite, PublicCRNCRBenchmarkSuite)

    def test_classified_dataset_exposes_additive_convenience_methods(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            _, classified = _real_dataset(Path(temporary_directory))
            split = classified.data_split(group_criteria="none", diversity="none")
            self.assertTrue(split.leakage_audit["passed"])
            cohorts = classified.build_cr_ncr_cohorts(
                ncr_pool_fractions=(0, 1),
                seeds=(42,),
                group_criteria="none",
                diversity="none",
            )
            self.assertEqual(len(cohorts.cohorts), 2)

    def test_classified_convenience_method_signatures_expose_public_options(self):
        classified = _classified()
        for bound, standalone in (
            (classified.data_split, data_split),
            (classified.build_cr_ncr_cohorts, build_cr_ncr_cohorts),
            (classified.build_cr_ncr_benchmark, build_cr_ncr_benchmark),
        ):
            self.assertEqual(
                tuple(inspect.signature(bound).parameters.values()),
                tuple(inspect.signature(standalone).parameters.values())[1:],
            )
            self.assertTrue(
                all(
                    parameter.kind is not inspect.Parameter.VAR_KEYWORD
                    for parameter in inspect.signature(bound).parameters.values()
                )
            )

    def test_aliases_normalize_to_canonical_receipt_names(self):
        self.assertEqual(
            normalize_group_criteria(("RT", "M2T", "SM")),
            (
                "rac5_crystalnets",
                "mofid_v2_crystalnets",
                "structure_matcher_strict",
            ),
        )
        self.assertEqual(
            normalize_group_criteria("rac5_topology"),
            ("rac5_crystalnets",),
        )
        with self.assertRaisesRegex(BenchmarkError, "duplicate"):
            normalize_group_criteria(("RT", "rac5_crystalnets"))

    def test_available_criteria_explains_missing_optional_release_evidence(self):
        classified = _classified()
        report = available_group_criteria(classified.dataset)
        self.assertTrue(report["priority_main"]["available"])
        self.assertTrue(report["rac5"]["available"])
        self.assertFalse(report["rac5_crystalnets"]["available"])
        self.assertIn("RT", report["rac5_crystalnets"]["aliases"])
        self.assertIn("required_evidence", report["rac5_crystalnets"])
        with self.assertRaises(TypeError):
            report["priority_main"] = {}

    def test_hidden_bridge_is_retained_before_label_filtering(self):
        ids = ("A", "HIDDEN", "C", "D")
        rows = [
            {
                "structure_id": "A",
                "source_database": "COD",
                "source_id": "A",
                "structure_variant": "ASR",
                "metal_elements": "Cu",
            },
            {
                "structure_id": "HIDDEN",
                "source_database": "COD",
                "source_id": "BRIDGE",
                "structure_variant": "ASR",
                "metal_elements": "Cu",
            },
            {
                "structure_id": "C",
                "source_database": "COD",
                "source_id": "BRIDGE",
                "structure_variant": "ASR",
                "metal_elements": "Cu",
            },
            {
                "structure_id": "D",
                "source_database": "COD",
                "source_id": "D",
                "structure_variant": "ASR",
                "metal_elements": "Cu",
            },
        ]
        hashes = {
            "A": "1" * 64,
            "HIDDEN": "1" * 64,
            "C": "2" * 64,
            "D": "3" * 64,
        }
        dataset = _Dataset(ids, hashes=hashes, rows=rows)
        labels = {"A": "CR", "HIDDEN": "AMBIGUOUS", "C": "NCR", "D": "CR"}
        result = data_split(
            _Classified(dataset, labels),
            group_criteria="none",
            diversity="none",
        )
        self.assertEqual(
            result.effective_leakage_blocks["A"],
            result.effective_leakage_blocks["C"],
        )
        self.assertEqual(result.assignments["A"], result.assignments["C"])
        self.assertEqual(result.exclusions["HIDDEN"], "LABEL_FILTER")
        self.assertTrue(result.leakage_audit["passed"])
        self.assertEqual(
            result.receipt()["effective_leakage_blocks"][
                "selected_mixed_label_block_count"
            ],
            1,
        )

    def test_general_split_rejects_unknown_label_filter_values(self):
        classified = _classified()
        with self.assertRaisesRegex(BenchmarkError, "unknown labels: CCR"):
            data_split(
                classified,
                labels=("CR", "CCR"),
                group_criteria="none",
                diversity="none",
            )

    def test_selected_criteria_are_unioned_and_missing_values_are_singletons(self):
        ids = ("A", "B", "C", "D")
        parents = {
            "A": {"rac_status": "MATCHED", "rac_group": "R-X", "rac_size": "2"},
            "C": {"rac_status": "MATCHED", "rac_group": "R-X", "rac_size": "2"},
        }
        dataset = _Dataset(ids, parents=parents)
        labels = {value: "CR" for value in ids}
        result = data_split(
            _Classified(dataset, labels),
            group_criteria=("rac5", "none"),
            diversity="none",
        )
        self.assertEqual(result.group_criteria, ("rac5", "none"))
        self.assertEqual(
            result.effective_leakage_blocks["A"],
            result.effective_leakage_blocks["C"],
        )
        self.assertNotEqual(
            result.criterion_groups["rac5"]["B"],
            result.criterion_groups["rac5"]["D"],
        )

    def test_row_order_does_not_change_assignments_or_digest(self):
        original = _classified(cr=18, ncr=6)
        reversed_rows = tuple(reversed(original.dataset.metadata_rows))
        reversed_dataset = _Dataset(
            tuple(row["structure_id"] for row in reversed_rows),
            hashes=original.dataset.cif_hashes,
            rows=reversed_rows,
        )
        reversed_view = _Classified(reversed_dataset, original.label_by_id)
        left = data_split(original, group_criteria="none", diversity="none")
        right = data_split(reversed_view, group_criteria="none", diversity="none")
        self.assertEqual(dict(left.assignments), dict(right.assignments))
        self.assertEqual(left.assignment_digest, right.assignment_digest)

    def test_large_indivisible_block_reports_partition_deviation(self):
        ids = tuple("X{}".format(index) for index in range(10))
        hashes = {
            value: (
                "f" * 64 if index < 9 else hashlib.sha256(value.encode()).hexdigest()
            )
            for index, value in enumerate(ids)
        }
        dataset = _Dataset(ids, hashes=hashes)
        result = data_split(
            _Classified(dataset, {value: "CR" for value in ids}),
            group_criteria="none",
            diversity="none",
        )
        deviations = result.diagnostics["partition_count_deviations"]
        self.assertTrue(any(abs(value) > 0 for value in deviations.values()))
        self.assertTrue(result.leakage_audit["passed"])

    def test_representative_backend_fails_closed_without_pinned_environment(self):
        if EXACT_BENCHMARK_BACKEND:
            self.skipTest("the exact benchmark backend is installed")
        classified = _classified()
        with self.assertRaisesRegex(BenchmarkDependencyError, "pinned benchmark extra"):
            data_split(classified, group_criteria="none", diversity="representative")

    def test_data_split_writer_restores_existing_pair_if_first_publish_fails(self):
        result = data_split(_classified(), group_criteria="none", diversity="none")
        with tempfile.TemporaryDirectory() as temporary_directory:
            csv_path, json_path = result.write(temporary_directory)
            original_csv = csv_path.read_bytes()
            original_json = json_path.read_bytes()
            real_replace = os.replace
            failed = []

            def fail_first_publish(source, target):
                if target == str(csv_path) and not failed:
                    failed.append((source, target))
                    raise OSError("injected first-publication failure")
                return real_replace(source, target)

            with patch(
                "CoREMOF._transactions.os.replace", side_effect=fail_first_publish
            ):
                with self.assertRaisesRegex(OSError, "injected"):
                    result.write(temporary_directory, overwrite=True)
            self.assertEqual(csv_path.read_bytes(), original_csv)
            self.assertEqual(json_path.read_bytes(), original_json)

    def test_data_split_writer_requires_an_exact_boolean_overwrite_flag(self):
        result = data_split(_classified(), group_criteria="none", diversity="none")
        with tempfile.TemporaryDirectory() as temporary_directory:
            result.write(temporary_directory)
            for value in (1, "yes", [], None):
                with (
                    self.subTest(value=value),
                    self.assertRaisesRegex(TypeError, "overwrite must be a boolean"),
                ):
                    result.write(temporary_directory, overwrite=value)

    def test_data_split_writer_never_removes_a_concurrent_foreign_generation(self):
        result = data_split(_classified(), group_criteria="none", diversity="none")
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            csv_path = directory / "race.csv"
            json_path = directory / "race.json"
            real_link = os.link
            real_replace = os.replace
            injected = []

            def inject_before_second_publication(source, target):
                if target == str(json_path) and not injected:
                    injected.append(target)
                    replacement = directory / "foreign.csv"
                    replacement.write_text("foreign-csv\n", encoding="utf-8")
                    real_replace(str(replacement), str(csv_path))
                    json_path.write_text("foreign-json\n", encoding="utf-8")
                return real_link(source, target)

            with (
                patch(
                    "CoREMOF._transactions.os.link",
                    side_effect=inject_before_second_publication,
                ),
                self.assertRaises(FileExistsError),
            ):
                result.write(directory, stem="race")
            self.assertEqual(csv_path.read_text(encoding="utf-8"), "foreign-csv\n")
            self.assertEqual(json_path.read_text(encoding="utf-8"), "foreign-json\n")

    def test_claimed_official_generation_is_revalidated(self):
        classified = _classified()
        object.__setattr__(classified, "checker_view", "4checker")
        with self.assertRaisesRegex(ValueError, "changed|differs"):
            data_split(classified, group_criteria="none", diversity="none")


class CRNCRBenchmarkTests(unittest.TestCase):
    def test_nested_solver_jointly_avoids_a_greedy_dead_end(self):
        ids = ("A", "B1", "B2", "C1", "C2", "C3")
        blocks = {
            "A": "BLOCK-1",
            "B1": "BLOCK-2",
            "B2": "BLOCK-2",
            "C1": "BLOCK-3",
            "C2": "BLOCK-3",
            "C3": "BLOCK-3",
        }
        strata = {structure_id: ("all",) for structure_id in ids}
        memberships = _nested_exact_block_memberships(
            ids,
            blocks,
            strata,
            target_counts=(3, 5),
            seed=42,
            namespace="joint-regression",
            label="synthetic",
        )
        self.assertEqual(len(memberships[3]), 3)
        self.assertEqual(len(memberships[5]), 5)
        self.assertTrue(set(memberships[3]).issubset(memberships[5]))
        for selected in memberships.values():
            selected = set(selected)
            for block in set(blocks.values()):
                members = {value for value in ids if blocks[value] == block}
                self.assertTrue(
                    not selected.intersection(members) or members <= selected
                )

    def test_nested_solver_exact_fallback_recovers_a_multibin_solution(self):
        sizes = (3, 3, 2, 2, 2)
        ids = tuple(
            "B{}-{}".format(block_index, member_index)
            for block_index, size in enumerate(sizes)
            for member_index in range(size)
        )
        blocks = {
            "B{}-{}".format(block_index, member_index): "BLOCK-{}".format(block_index)
            for block_index, size in enumerate(sizes)
            for member_index in range(size)
        }
        strata = {structure_id: ("all",) for structure_id in ids}
        memberships = _nested_exact_block_memberships(
            ids,
            blocks,
            strata,
            target_counts=(3, 7, 12),
            seed=42,
            namespace="joint-fallback-regression",
            label="synthetic",
        )
        self.assertEqual([len(memberships[value]) for value in (3, 7, 12)], [3, 7, 12])
        self.assertTrue(
            set(memberships[3]) <= set(memberships[7]) <= set(memberships[12])
        )

    def test_label_impure_blocks_require_explicit_sensitivity_policy(self):
        shared_cr_ncr = "1" * 64
        shared_cr_ambiguous = "2" * 64
        shared_ncr_unchecked = "3" * 64
        hashes = {
            "ID-000": shared_cr_ncr,
            "ID-006": shared_cr_ncr,
            "ID-001": shared_cr_ambiguous,
            "ID-008": shared_cr_ambiguous,
            "ID-007": shared_ncr_unchecked,
            "ID-009": shared_ncr_unchecked,
        }
        _, classified = _authenticated_classified(
            Path("/synthetic-impure-release"),
            cr=6,
            ncr=2,
            ambiguous=1,
            unchecked=1,
            cif_hashes=hashes,
        )
        with self.assertRaisesRegex(
            BenchmarkFeasibilityError,
            "cannot be selected as whole effective leakage blocks.*2 CR and 2 NCR",
        ):
            build_cr_ncr_cohorts(
                classified,
                ncr_pool_fractions=(0, 1),
                seeds=(42,),
                train=0.8,
                val=0.2,
                test=0.0,
                group_criteria="none",
                diversity="none",
            )
        cohorts = build_cr_ncr_cohorts(
            classified,
            ncr_pool_fractions=(0, 1),
            seeds=(42,),
            train=0.8,
            val=0.2,
            test=0.0,
            group_criteria="none",
            cohort_eligibility=LABEL_PURE_BLOCK_ELIGIBILITY,
            diversity="none",
        )
        self.assertEqual(cohorts.raw_pool_counts, {"CR": 6, "NCR": 2})
        self.assertEqual(cohorts.pool_counts, {"CR": 4, "NCR": 0})
        self.assertEqual(
            cohorts.receipt()["excluded_strict_membership"],
            {
                "CR": 2,
                "NCR": 2,
                "reason": "BLOCK_NOT_SINGLE_STRICT_LABEL",
                "membership_sha256": cohorts.receipt()["excluded_strict_membership"][
                    "membership_sha256"
                ],
                "companion_manifest": "label_accounting_manifest.csv",
            },
        )
        self.assertTrue(
            all(
                cohort.whole_block_membership_audit["passed"]
                for cohort in cohorts.cohorts
            )
        )

    def test_nested_cohorts_add_and_remove_only_complete_blocks(self):
        hashes = {}
        for first, second, digest in (
            ("ID-000", "ID-001", "1" * 64),
            ("ID-002", "ID-003", "2" * 64),
            ("ID-008", "ID-009", "3" * 64),
            ("ID-010", "ID-011", "4" * 64),
        ):
            hashes[first] = digest
            hashes[second] = digest
        _, classified = _authenticated_classified(
            Path("/synthetic-whole-block-release"),
            cr=8,
            ncr=4,
            cif_hashes=hashes,
        )
        cohorts = build_cr_ncr_cohorts(
            classified,
            ncr_pool_fractions=(0, 0.5, 1),
            seeds=(42,),
            train=0.8,
            val=0.2,
            test=0.0,
            group_criteria="none",
            diversity="none",
        )
        self.assertEqual(
            [cohort.selected_ncr_count for cohort in cohorts.cohorts],
            [0, 2, 4],
        )
        for cohort in cohorts.cohorts:
            self.assertEqual(
                cohort.whole_block_membership_audit["partially_selected_block_count"],
                0,
            )

    def test_unreachable_formula_count_fails_without_splitting_a_block(self):
        hashes = {}
        for offset, digest in ((0, "1" * 64), (2, "2" * 64), (4, "3" * 64)):
            hashes["ID-{:03d}".format(8 + offset)] = digest
            hashes["ID-{:03d}".format(9 + offset)] = digest
        _, classified = _authenticated_classified(
            Path("/synthetic-unreachable-block-release"),
            cr=8,
            ncr=6,
            cif_hashes=hashes,
        )
        with self.assertRaisesRegex(
            BenchmarkFeasibilityError,
            "exact whole-block NCR addition.*requested 3.*nearest reachable",
        ):
            build_cr_ncr_cohorts(
                classified,
                ncr_pool_fractions=(0, 0.5, 1),
                seeds=(42,),
                train=0.8,
                val=0.2,
                test=0.0,
                group_criteria="none",
                diversity="none",
            )

    def test_fixed_size_nested_paired_suite_and_common_clean_test(self):
        classified = _classified(cr=20, ncr=8, ambiguous=2)
        suite = build_cr_ncr_benchmark(
            classified,
            ncr_pool_fractions=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
            seeds=(42, 43, 44, 45, 46),
            group_criteria="none",
            diversity="none",
        )
        self.assertEqual(len(suite.runs), 30)
        self.assertFalse(suite.official_split)
        cohort_receipt = suite.receipt()["cohort_receipt"]
        self.assertEqual(
            cohort_receipt["input_hashes"],
            dict(sorted(classified.dataset.input_hashes.items())),
        )
        self.assertEqual(
            cohort_receipt["checker_names"],
            ["MOFClassifier", "MOFChecker", "Chen-Manz", "MOSAEC", "SETC-GAT"],
        )
        self.assertIn("effective_leakage_policy_definition", cohort_receipt)
        self.assertIn("test_policy_definition", cohort_receipt)
        self.assertIn("test_policy_definition", suite.receipt())
        self.assertEqual(
            cohort_receipt["complete_release_label_accounting"]["counts"],
            {"CR": 20, "NCR": 8, "AMBIGUOUS": 2, "UNCHECKED": 0},
        )
        self.assertEqual(len(cohort_receipt["effective_leakage_blocks"]["sha256"]), 64)
        self.assertTrue(suite.fixed_test_ids)
        self.assertTrue(all(run.test_ids == suite.fixed_test_ids for run in suite.runs))
        self.assertTrue(all(len(run.assignments) == 20 for run in suite.runs))
        self.assertTrue(all(run.leakage_audit["passed"] for run in suite.runs))
        for seed in suite.seeds:
            runs = [run for run in suite.runs if run.seed == seed]
            runs.sort(key=lambda run: float(run.requested_ncr_pool_fraction))
            self.assertEqual(
                set(runs[-1].ncr_ids), set(classified.dataset.structure_ids[20:28])
            )
            for previous, current in zip(runs, runs[1:]):
                self.assertTrue(set(previous.ncr_ids).issubset(current.ncr_ids))
                self.assertTrue(set(current.cr_ids).issubset(previous.cr_ids))
                common = set(previous.assignments).intersection(current.assignments)
                self.assertTrue(
                    all(
                        previous.assignments[value] == current.assignments[value]
                        for value in common
                    )
                )
            self.assertTrue(all(run.achieved_counts["validation"] > 0 for run in runs))
            endpoint = runs[-1]
            self.assertEqual(endpoint.label_counts_by_partition["validation"]["NCR"], 1)
            for mixed in runs[3:]:
                self.assertEqual(
                    mixed.label_counts_by_partition["validation"],
                    {"CR": 1, "NCR": 1},
                )
            self.assertIn(
                "label_count_deviations_by_partition",
                suite.receipt()["runs"][suite.runs.index(endpoint)],
            )
            for run in runs:
                for partition in ("train", "validation", "test"):
                    self.assertEqual(
                        sum(run.label_counts_by_partition[partition].values()),
                        run.achieved_counts[partition],
                    )
        diagnostic = suite.runs[0].full_cr_diagnostic
        self.assertEqual(diagnostic["structure_count"], 20)
        self.assertGreater(diagnostic["exact_training_overlap_count"], 0)

    def test_strict_benchmark_rejects_an_unauthenticated_label_container(self):
        ids = ("A", "B")
        classified = _Classified(_Dataset(ids), {"A": "CR", "B": "NCR"})
        with self.assertRaisesRegex(BenchmarkError, "authenticated ClassifiedDataset"):
            build_cr_ncr_cohorts(
                classified,
                ncr_pool_fractions=(0,),
                seeds=(42,),
                train=0.5,
                val=0.0,
                test=0.5,
                group_criteria="none",
                diversity="none",
            )

    def test_fixed_test_never_adds_a_block_that_worsens_count_distance(self):
        ids = tuple("A{:02d}".format(index) for index in range(21))
        blocks = {
            structure_id: "small" if index < 6 else "large"
            for index, structure_id in enumerate(ids)
        }
        selected = _fixed_test_blocks(
            ids,
            {structure_id: "CR" for structure_id in ids},
            blocks,
            {structure_id: ("CR", "none") for structure_id in ids},
            desired_count=10,
        )
        self.assertEqual(len(selected), 6)
        tied_ids = tuple("T{:02d}".format(index) for index in range(20))
        tied = _fixed_test_blocks(
            tied_ids,
            {structure_id: "CR" for structure_id in tied_ids},
            {
                structure_id: "below" if index < 8 else "above"
                for index, structure_id in enumerate(tied_ids)
            },
            {structure_id: ("CR", "none") for structure_id in tied_ids},
            desired_count=10,
        )
        self.assertEqual(len(tied), 8)

    def test_fixed_test_finds_the_globally_closest_whole_block_count(self):
        ids = tuple("G{:03d}".format(index) for index in range(100))
        boundaries = (
            (0, 8, "eight"),
            (8, 15, "seven"),
            (15, 18, "three"),
            (18, 100, "rest"),
        )
        blocks = {
            structure_id: next(
                block for start, stop, block in boundaries if start <= index < stop
            )
            for index, structure_id in enumerate(ids)
        }
        selected = _fixed_test_blocks(
            ids,
            {structure_id: "CR" for structure_id in ids},
            blocks,
            {structure_id: ("CR", "none") for structure_id in ids},
            desired_count=10,
        )
        self.assertEqual(len(selected), 10)
        self.assertEqual({blocks[value] for value in selected}, {"seven", "three"})

    def test_generated_benchmark_receipt_passes_first_use_terminology_audit(self):
        auditor = (
            Path(__file__).resolve().parents[2]
            / "CoREMOF-COD"
            / "dataset_split"
            / "scripts"
            / "audit_first_use_terminology.py"
        )
        if not auditor.is_file():
            self.skipTest("CoREMOF-COD terminology auditor is not available")
        suite = build_cr_ncr_benchmark(
            _classified(),
            ncr_pool_fractions=(0, 1),
            seeds=(42,),
            group_criteria="none",
            diversity="none",
        )
        with tempfile.TemporaryDirectory() as temporary_directory:
            receipt_path = Path(temporary_directory) / "receipt.json"
            receipt_path.write_text(
                json.dumps(suite.receipt(), indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            for no_site in (False, True):
                command = [sys.executable]
                if no_site:
                    command.append("-S")
                command.extend(["-B", str(auditor), str(receipt_path)])
                completed = subprocess.run(
                    command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=False,
                )
                self.assertEqual(completed.returncode, 0, completed.stdout)

    def test_cohorts_are_inspectable_before_partitioning(self):
        classified = _classified(cr=12, ncr=4)
        cohorts = build_cr_ncr_cohorts(
            classified,
            ncr_pool_fractions=(0, 0.5, 1),
            seeds=(42,),
            group_criteria="none",
            diversity="none",
        )
        self.assertEqual(len(cohorts.cohorts), 3)
        self.assertEqual(cohorts.pool_counts, {"CR": 12, "NCR": 4})
        suite = cohorts.data_split()
        self.assertEqual(len(suite.runs), 3)
        with self.assertRaisesRegex(BenchmarkError, "rebuild the cohorts"):
            cohorts.data_split(train=0.7, val=0.2, test=0.1)
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = cohorts.write(temporary_directory)
            checksum_path = root / "SHA256SUMS"
            self.assertTrue(checksum_path.is_file())
            for line in checksum_path.read_text(encoding="utf-8").splitlines():
                digest, relative = line.split("  ", 1)
                self.assertEqual(
                    hashlib.sha256((root / relative).read_bytes()).hexdigest(),
                    digest,
                )

    def test_round_half_up_formula(self):
        classified = _classified(cr=12, ncr=5)
        cohorts = build_cr_ncr_cohorts(
            classified,
            ncr_pool_fractions=(0.1, 0.5, 1.0),
            seeds=(42,),
            group_criteria="none",
            diversity="none",
        )
        self.assertEqual(
            [cohort.selected_ncr_count for cohort in cohorts.cohorts],
            [1, 3, 5],
        )

    def test_m_greater_than_c_fails_before_diversity_backend(self):
        classified = _classified(cr=4, ncr=5)
        with self.assertRaisesRegex(
            BenchmarkFeasibilityError, r"C=4 CR and M=5 NCR.*M>C"
        ):
            build_cr_ncr_benchmark(classified, diversity="representative")

    def test_clean_test_capacity_failure_reports_maximum_fraction(self):
        classified = _classified(cr=10, ncr=10)
        with self.assertRaisesRegex(
            BenchmarkFeasibilityError,
            r"after reserving.*C=10 CR, M=10 NCR.*maximum feasible",
        ):
            build_cr_ncr_cohorts(
                classified,
                group_criteria="none",
                diversity="none",
            )

    def test_zero_test_fraction_reserves_no_test_block(self):
        suite = build_cr_ncr_benchmark(
            _classified(cr=10, ncr=10),
            ncr_pool_fractions=(0, 1),
            seeds=(42,),
            train=0.9,
            val=0.1,
            test=0.0,
            group_criteria="none",
            diversity="none",
        )
        self.assertEqual(suite.fixed_test_ids, ())
        self.assertTrue(all(run.test_ids == () for run in suite.runs))
        self.assertTrue(all(len(run.assignments) == 10 for run in suite.runs))

    def test_zero_validation_fraction_never_receives_a_block(self):
        suite = build_cr_ncr_benchmark(
            _classified(cr=12, ncr=4),
            ncr_pool_fractions=(0, 1),
            seeds=(42,),
            train=0.9,
            val=0.0,
            test=0.1,
            group_criteria="none",
            diversity="none",
        )
        self.assertTrue(all(run.validation_ids == () for run in suite.runs))
        self.assertTrue(
            all("validation" not in set(run.assignments.values()) for run in suite.runs)
        )

    def test_benchmark_writer_emits_all_manifests_and_valid_checksums(self):
        suite = build_cr_ncr_benchmark(
            _classified(cr=12, ncr=4),
            ncr_pool_fractions=(0, 0.5, 1),
            seeds=(42, 43),
            group_criteria="none",
            diversity="none",
        )
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = suite.write(temporary_directory)
            for name in (
                "suite_index.json",
                "membership_manifest.csv",
                "fixed_test_manifest.csv",
                "full_cr_diagnostic_manifest.csv",
                "label_accounting_manifest.csv",
                "receipt.json",
                "SHA256SUMS",
            ):
                self.assertTrue((root / name).is_file())
            lines = (root / "SHA256SUMS").read_text(encoding="utf-8").splitlines()
            self.assertTrue(lines)
            for line in lines:
                digest, relative = line.split("  ", 1)
                self.assertEqual(
                    hashlib.sha256((root / relative).read_bytes()).hexdigest(),
                    digest,
                )
            first_run = suite.runs[0]
            with (root / "runs" / (first_run.run_key + ".csv")).open(
                "r", encoding="utf-8", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle))
            verified = frozen_assignment_manifest(rows, suite.receipt())
            self.assertEqual(verified.assignment_digest, first_run.assignment_digest)
            self.assertEqual(dict(verified.assignments), dict(first_run.assignments))
            with self.assertRaises(FileExistsError):
                suite.write(temporary_directory)

    def test_label_accounting_eligibility_does_not_depend_on_requested_q(self):
        suite = build_cr_ncr_benchmark(
            _classified(cr=10, ncr=4),
            ncr_pool_fractions=(0.5,),
            seeds=(42,),
            train=0.8,
            val=0.2,
            test=0.0,
            group_criteria="none",
            diversity="none",
        )
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = suite.write(temporary_directory)
            with (root / "label_accounting_manifest.csv").open(
                "r", encoding="utf-8", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle))
        strict_rows = [row for row in rows if row["label"] in {"CR", "NCR"}]
        self.assertEqual(len(strict_rows), 14)
        self.assertTrue(all(row["cohort_eligible"] == "true" for row in strict_rows))
        self.assertTrue(all(row["exclusion_reason"] == "" for row in strict_rows))

    def test_disabled_full_cr_diagnostic_omits_the_manifest(self):
        suite = build_cr_ncr_benchmark(
            _classified(cr=12, ncr=4),
            ncr_pool_fractions=(0, 1),
            seeds=(42,),
            group_criteria="none",
            diversity="none",
            include_full_cr_diagnostic=False,
        )
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = suite.write(temporary_directory)
            self.assertFalse((root / "full_cr_diagnostic_manifest.csv").exists())

    def test_directory_publication_is_atomic_create_if_absent(self):
        suite = build_cr_ncr_benchmark(
            _classified(cr=12, ncr=4),
            ncr_pool_fractions=(0, 1),
            seeds=(42,),
            group_criteria="none",
            diversity="none",
        )
        from CoREMOF import _transactions

        real_publish = _transactions._rename_directory_noreplace
        injected = []

        def concurrent_directory(source, target):
            if not injected:
                injected.append(target)
                target.mkdir()
                (target / "foreign.txt").write_text("foreign\n", encoding="utf-8")
            return real_publish(source, target)

        with tempfile.TemporaryDirectory() as temporary_directory:
            with (
                patch(
                    "CoREMOF._transactions._rename_directory_noreplace",
                    side_effect=concurrent_directory,
                ),
                self.assertRaises(FileExistsError),
            ):
                suite.write(temporary_directory, stem="race")
            target = Path(temporary_directory) / "race"
            self.assertEqual(
                (target / "foreign.txt").read_text(encoding="utf-8"), "foreign\n"
            )

    def test_directory_publication_uses_locked_fallback_when_noreplace_is_unsupported(
        self,
    ):
        from CoREMOF import _transactions

        class UnsupportedRename:
            argtypes = None
            restype = None

            def __call__(self, *args):
                ctypes.set_errno(errno.EINVAL)
                return -1

        class FakeLibc:
            renameat2 = UnsupportedRename()

        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            staging = directory / ".staging"
            staging.mkdir()
            (staging / "payload.txt").write_text("complete\n", encoding="utf-8")
            target = directory / "published"
            with patch("CoREMOF._transactions.ctypes.CDLL", return_value=FakeLibc()):
                _transactions._rename_directory_noreplace(staging, target)
            self.assertEqual(
                (target / "payload.txt").read_text(encoding="utf-8"), "complete\n"
            )

            second_staging = directory / ".second-staging"
            second_staging.mkdir()
            with (
                patch("CoREMOF._transactions.ctypes.CDLL", return_value=FakeLibc()),
                self.assertRaises(FileExistsError),
            ):
                _transactions._rename_directory_noreplace(second_staging, target)
            self.assertTrue(second_staging.is_dir())

    def test_feasible_fraction_receipt_is_limited_to_the_q_domain(self):
        cohorts = build_cr_ncr_cohorts(
            _classified(cr=20, ncr=8),
            ncr_pool_fractions=(0, 1),
            seeds=(42,),
            group_criteria="none",
            diversity="none",
        )
        self.assertEqual(
            cohorts.receipt()["feasibility_audit"][
                "maximum_feasible_NCR_pool_fraction"
            ],
            "1",
        )
        self.assertIn(
            "targets.py", cohorts.receipt()["implementation"]["source_sha256"]
        )


@unittest.skipUnless(
    EXACT_BENCHMARK_BACKEND,
    "requires exact NumPy/scikit-learn benchmark extra",
)
class RepresentativeDiversityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary_directory = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temporary_directory.name)
        cls.dataset, cls.classified = _real_dataset(cls.root, cr=50, ncr=10)
        ids = cls.dataset.structure_ids
        rac_fields = tuple("rac_{:03d}".format(index) for index in range(264))
        rac_rows = []
        for row_index, structure_id in enumerate(ids):
            row = {
                "structure_id": structure_id,
                "rac5_available": "true" if row_index <= 30 else "false",
            }
            for field_index, field in enumerate(rac_fields):
                if row_index <= 30:
                    row[field] = (
                        "7.0"
                        if field_index == 263
                        else str((row_index + 1) * ((field_index % 7) + 1))
                    )
                else:
                    row[field] = ""
            if row_index == 30:
                row[rac_fields[0]] = ""
            rac_rows.append(row)
        _write_csv(
            cls.root / "features" / "rac5_features.csv",
            ("structure_id", "rac5_available") + rac_fields,
            rac_rows,
        )

        zeo_fields = tuple(ZEO_NUMERIC_FINGERPRINT_DEFINITION["numeric_fields"]) + (
            "n2_channel_dimension",
            "structure_periodic_dimension",
        )
        zeo_rows = []
        for row_index, structure_id in enumerate(ids):
            available = 30 <= row_index <= 50
            row = {
                "structure_id": structure_id,
                "n2_he_available": "true" if available else "false",
                "periodicity_available": "true" if available else "false",
            }
            for field_index, field in enumerate(zeo_fields):
                row[field] = (
                    str((row_index + 1) / (field_index + 2)) if available else ""
                )
            if row_index == 50:
                row[zeo_fields[0]] = ""
            zeo_rows.append(row)
        _write_csv(
            cls.root / "features" / "zeo_features.csv",
            (
                "structure_id",
                "n2_he_available",
                "periodicity_available",
            )
            + zeo_fields,
            zeo_rows,
        )

        topology_rows = [
            {
                "structure_id": structure_id,
                "topology_available": "true",
                "network_dimension": "3",
                "single_node_net": "pcu",
                "all_node_net": "pcu",
                "single_all_agree": "true",
            }
            for structure_id in ids
        ]
        _write_csv(
            cls.root / "features" / "topology_features.csv",
            (
                "structure_id",
                "topology_available",
                "network_dimension",
                "single_node_net",
                "all_node_net",
                "single_all_agree",
            ),
            topology_rows,
        )

    @classmethod
    def tearDownClass(cls):
        cls.temporary_directory.cleanup()

    def test_tiers_retain_missing_rows_without_imputation_and_are_deterministic(self):
        first = build_diversity_index(self.dataset)
        reversed_dataset = CoREMOFDataset(
            release_root=self.dataset.release_root,
            records=tuple(reversed(self.dataset.records)),
            dataset_info=self.dataset.dataset_info,
            parent_group_methods=self.dataset.parent_group_methods,
            parent_by_id=self.dataset.parent_by_id,
            input_hashes=self.dataset.input_hashes,
            cif_files_verified=self.dataset.cif_files_verified,
        )
        second = build_diversity_index(reversed_dataset)
        self.assertEqual(dict(first.strata_by_id), dict(second.strata_by_id))
        self.assertEqual(first.digest, second.digest)
        counts = Counter(first.tier_by_id.values())
        self.assertEqual(counts, {"rac5": 30, "zeo": 20, "no_numeric": 10})
        self.assertEqual(first.tier_by_id[self.dataset.structure_ids[30]], "zeo")
        self.assertEqual(first.tier_by_id[self.dataset.structure_ids[50]], "no_numeric")
        self.assertFalse(first.profile["scientific_feature_imputation"])
        backend = first.profile["backend"]
        self.assertEqual(backend["scipy"], BENCHMARK_SCIPY_VERSION)
        self.assertEqual(backend["joblib"], BENCHMARK_JOBLIB_VERSION)
        self.assertEqual(backend["threadpoolctl"], BENCHMARK_THREADPOOLCTL_VERSION)
        self.assertEqual(backend["thread_limit"], 1)
        self.assertFalse(backend["cross_architecture_bit_identity_guaranteed"])
        self.assertTrue(
            all("filepath" not in record for record in backend["threadpool_libraries"])
        )
        self.assertGreaterEqual(
            first.profile["tier_diagnostics"]["rac5"]["zero_iqr_field_count"],
            1,
        )

    def test_zero_iqr_unit_divisor_retains_outlier_deviation(self):
        captured = {}

        class CapturingKMeans:
            def __init__(self, **kwargs):
                pass

            def fit_predict(self, matrix):
                captured["matrix"] = matrix.copy()
                return _benchmark_numpy.arange(matrix.shape[0])

        ids = tuple("Z{}".format(index) for index in range(5))
        vectors = {
            structure_id: ((10.0,) if index == 4 else (0.0,))
            for index, structure_id in enumerate(ids)
        }
        _, details = _cluster_tier(
            ids,
            vectors,
            "zeo",
            _benchmark_numpy,
            object,
            CapturingKMeans,
        )
        self.assertEqual(details["zero_iqr_field_count"], 1)
        self.assertEqual(captured["matrix"][:, 0].tolist(), [0.0, 0.0, 0.0, 0.0, 10.0])

    def test_split_retains_every_tier_and_balances_proportional_coverage(self):
        result = data_split(
            self.classified,
            group_criteria="none",
            diversity="representative",
            train=0.6,
            val=0.2,
            test=0.2,
            random_state=42,
        )
        self.assertEqual(len(result.assignments), len(self.dataset.structure_ids))
        self.assertEqual(result.receipt()["target_columns_consumed"], [])
        tier_totals = Counter(result.diversity_tiers.values())
        fractions = {"train": 0.6, "validation": 0.2, "test": 0.2}
        for tier, total in tier_totals.items():
            observed = Counter(
                result.assignments[structure_id]
                for structure_id, value in result.diversity_tiers.items()
                if value == tier
            )
            for partition, fraction in fractions.items():
                self.assertLessEqual(
                    abs(observed[partition] - fraction * total),
                    3.0,
                    msg="{} {} coverage".format(tier, partition),
                )
        self.assertTrue(result.leakage_audit["passed"])

    def test_paired_benchmark_uses_the_same_reusable_diversity_index(self):
        index = build_diversity_index(self.dataset)
        suite = build_cr_ncr_benchmark(
            self.classified,
            ncr_pool_fractions=(0.0, 0.5, 1.0),
            seeds=(42, 43),
            group_criteria="none",
            diversity="representative",
        )
        self.assertEqual(suite.diversity_index_hash, index.digest)
        self.assertEqual(len(suite.runs), 6)
        self.assertTrue(all(run.leakage_audit["passed"] for run in suite.runs))
        self.assertEqual(
            set(suite.input_hashes).intersection(
                {
                    "features/rac5_features.csv",
                    "features/zeo_features.csv",
                    "features/topology_features.csv",
                }
            ),
            {
                "features/rac5_features.csv",
                "features/zeo_features.csv",
                "features/topology_features.csv",
            },
        )


class DeferredTargetAttachmentTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.dataset, self.classified = _real_dataset(self.root)
        self.target_path = self.root / "targets.csv"
        _write_csv(
            self.target_path,
            ("structure_id", "uptake"),
            [
                {
                    "structure_id": structure_id,
                    "uptake": "" if index in (1, 12) else str(index + 0.5),
                }
                for index, structure_id in enumerate(self.dataset.structure_ids)
            ],
        )
        self.source = TargetSource(
            self.target_path,
            target_columns=("uptake",),
            value_types={"uptake": "float"},
            units={"uptake": "mol/kg"},
            conditions={"uptake": {"temperature_K": 298}},
        )

    def tearDown(self):
        self.temporary_directory.cleanup()

    def test_keep_error_and_drop_preserve_original_assignment_digest(self):
        split = data_split(self.classified, group_criteria="none", diversity="none")
        original_receipt = json.dumps(split.receipt(), sort_keys=True)
        kept = split.attach_targets(self.source, missing="keep")
        self.assertIsInstance(kept, PublicTargetAttachedView)
        self.assertEqual(kept.original_structure_ids, tuple(sorted(split.assignments)))
        self.assertIsNone(kept.values_by_id[self.dataset.structure_ids[1]]["uptake"])
        self.assertEqual(kept.original_assignment_digest, split.assignment_digest)
        self.assertIs(kept.official_split, False)
        self.assertIs(kept.receipt()["official_split"], False)
        with self.assertRaisesRegex(TargetAttachmentError, "completed target set"):
            split.attach_targets(self.source, missing="error")
        dropped = split.attach_targets(self.source, missing="drop")
        self.assertIn(self.dataset.structure_ids[1], dropped.dropped_ids)
        self.assertNotIn(self.dataset.structure_ids[1], dropped.assignments)
        missing_partition = split.assignments[self.dataset.structure_ids[1]]
        self.assertGreaterEqual(
            dropped.missing_counts_by_partition[missing_partition]["uptake"], 1
        )
        self.assertEqual(
            dropped.derived_missing_counts_by_partition.get(missing_partition, {}).get(
                "uptake", 0
            ),
            0,
        )
        self.assertEqual(
            dropped.receipt()["release_binding"]["input_hashes"],
            dict(self.dataset.input_hashes),
        )
        for structure_id, partition in dropped.assignments.items():
            self.assertEqual(partition, split.assignments[structure_id])
        self.assertEqual(json.dumps(split.receipt(), sort_keys=True), original_receipt)

    def test_live_assignment_is_revalidated_against_its_frozen_receipt(self):
        other_dataset, _ = _authenticated_classified(
            self.root / "mutated-release",
            cr=10,
            ncr=4,
            metadata_hash="c" * 64,
        )
        split = data_split(self.classified, group_criteria="none", diversity="none")
        object.__setattr__(split, "_dataset", other_dataset)
        with self.assertRaisesRegex(BenchmarkError, "dataset binding changed"):
            split.attach_targets(self.source)

        split = data_split(self.classified, group_criteria="none", diversity="none")
        object.__setattr__(split, "assignment_digest", "f" * 64)
        with self.assertRaisesRegex(BenchmarkError, "changed after authentication"):
            split.attach_targets(self.source)

        split = data_split(self.classified, group_criteria="none", diversity="none")
        object.__setattr__(split, "official_split", True)
        with self.assertRaisesRegex(BenchmarkError, "changed after authentication"):
            split.attach_targets(self.source)

        split = data_split(self.classified, group_criteria="none", diversity="none")
        stripped_receipt = split.receipt()
        stripped_receipt.pop("release_binding")
        object.__setattr__(split, "_receipt", stripped_receipt)
        with self.assertRaisesRegex(BenchmarkError, "changed after authentication"):
            split.attach_targets(self.source)

        suite = build_cr_ncr_benchmark(
            self.classified,
            ncr_pool_fractions=(0, 1),
            seeds=(42,),
            train=0.8,
            val=0.2,
            test=0.0,
            group_criteria="none",
            diversity="none",
        )
        object.__setattr__(suite.runs[0], "assignment_digest", "e" * 64)
        with self.assertRaisesRegex(BenchmarkError, "changed after authentication"):
            suite.attach_targets(self.source)

    def test_premerged_target_generation_is_revalidated_before_attachment(self):
        split = data_split(self.classified, group_criteria="none", diversity="none")
        for field_name in (
            "target_values_by_id",
            "target_definitions",
            "target_input_receipt",
        ):
            with self.subTest(field_name=field_name):
                merged = merge_targets(self.dataset, self.source)
                object.__setattr__(merged, field_name, {})
                with self.assertRaisesRegex(
                    TargetAttachmentError, "sealed integrity check"
                ):
                    split.attach_targets(merged)

    def test_receipt_source_hashes_are_bound_at_module_import(self):
        with (
            patch.object(benchmarks_module, "_current_source_hashes", return_value={}),
            self.assertRaisesRegex(BenchmarkError, "changed after module import"),
        ):
            benchmarks_module._source_hashes()
        with (
            patch.object(attachments_module, "_current_source_hashes", return_value={}),
            self.assertRaisesRegex(
                TargetAttachmentError, "changed after module import"
            ),
        ):
            attachments_module._source_hashes()

    def test_frozen_assignment_requires_verified_factory_and_known_schema(self):
        binding = {
            "dataset_version": "vtest",
            "structure_count": 1,
            "structure_ids_sha256": "f" * 64,
            "input_hashes": {"metadata/metadata.csv": "b" * 64},
        }
        with self.assertRaisesRegex(TargetAttachmentError, "must be created"):
            FrozenAssignmentManifest(
                {"A": "train"},
                "0" * 64,
                binding,
                official_split=False,
            )

        split = data_split(self.classified, group_criteria="none", diversity="none")
        unknown = split.receipt()
        unknown["schema_version"] = "unrecognized-assignment-receipt/1.0"
        with self.assertRaisesRegex(TargetAttachmentError, "schema is not accepted"):
            frozen_assignment_manifest(split.assignment_rows(), unknown)

        claimed_official = split.receipt()
        claimed_official["official_split"] = True
        with self.assertRaisesRegex(TargetAttachmentError, "official assignment"):
            frozen_assignment_manifest(split.assignment_rows(), claimed_official)

    def test_verified_frozen_assignment_detects_post_factory_mutation(self):
        split = data_split(self.classified, group_criteria="none", diversity="none")
        frozen = frozen_assignment_manifest(split.assignment_rows(), split.receipt())
        object.__setattr__(frozen, "assignments", {"FORGED": "test"})
        with self.assertRaisesRegex(TargetAttachmentError, "changed after"):
            attach_targets(frozen, self.source, dataset=self.dataset)

        frozen = frozen_assignment_manifest(split.assignment_rows(), split.receipt())
        object.__setattr__(frozen, "official_split", True)
        with self.assertRaisesRegex(
            TargetAttachmentError, "changed after authentication"
        ):
            attach_targets(frozen, self.source, dataset=self.dataset)

    def test_authority_does_not_transfer_through_dataclass_replace(self):
        split = data_split(self.classified, group_criteria="none", diversity="none")
        with self.assertRaisesRegex(BenchmarkError, "not produced"):
            replace(split).receipt()

        cohorts = build_cr_ncr_cohorts(
            self.classified,
            ncr_pool_fractions=(0, 1),
            seeds=(42,),
            train=0.8,
            val=0.2,
            test=0.0,
            group_criteria="none",
            diversity="none",
        )
        with self.assertRaisesRegex(BenchmarkError, "not produced"):
            replace(cohorts).data_split()

        suite = cohorts.data_split()
        with self.assertRaisesRegex(BenchmarkError, "not produced"):
            replace(suite).receipt()

        frozen = frozen_assignment_manifest(split.assignment_rows(), split.receipt())
        with self.assertRaisesRegex(TargetAttachmentError, "not produced"):
            attach_targets(replace(frozen), self.source, dataset=self.dataset)

        view = split.attach_targets(self.source)
        with self.assertRaisesRegex(TargetAttachmentError, "not produced"):
            replace(view).rows()

        attached_suite = suite.attach_targets(self.source)
        with self.assertRaisesRegex(TargetAttachmentError, "not produced"):
            replace(attached_suite).receipt()

    def test_live_receipt_requires_an_exact_known_result_type(self):
        split = data_split(self.classified, group_criteria="none", diversity="none")

        class PretendSplit:
            assignments = split.assignments
            official_split = False

            @staticmethod
            def receipt():
                return split.receipt()

        with self.assertRaisesRegex(TargetAttachmentError, "type is not an accepted"):
            attach_targets(PretendSplit(), self.source, dataset=self.dataset)

    def test_frozen_benchmark_receipt_dispatches_strictly_by_schema(self):
        suite = build_cr_ncr_benchmark(
            self.classified,
            ncr_pool_fractions=(0,),
            seeds=(42,),
            train=0.8,
            val=0.2,
            test=0.0,
            group_criteria="none",
            diversity="none",
        )
        run = suite.runs[0]
        rows = [
            {
                "run_key": run.run_key,
                "structure_id": structure_id,
                "partition": partition,
            }
            for structure_id, partition in sorted(run.assignments.items())
        ]
        malformed = suite.receipt()
        malformed["assignment_sha256"] = run.assignment_digest
        with self.assertRaisesRegex(TargetAttachmentError, "top-level"):
            frozen_assignment_manifest(rows, malformed)

    def test_receipt_boundaries_recheck_imported_source_hashes(self):
        split = data_split(self.classified, group_criteria="none", diversity="none")
        view = split.attach_targets(self.source)
        with (
            patch.object(benchmarks_module, "_current_source_hashes", return_value={}),
            self.assertRaisesRegex(BenchmarkError, "changed after module import"),
        ):
            split.receipt()
        with (
            patch.object(attachments_module, "_current_source_hashes", return_value={}),
            self.assertRaisesRegex(
                TargetAttachmentError, "changed after module import"
            ),
        ):
            view.receipt()

    def test_reserved_export_target_columns_are_rejected(self):
        split = data_split(self.classified, group_criteria="none", diversity="none")
        for target_name in ("partition", "run_key"):
            with self.subTest(target_name=target_name):
                path = self.root / (target_name + ".csv")
                _write_csv(
                    path,
                    ("structure_id", target_name),
                    [
                        {"structure_id": structure_id, target_name: "value"}
                        for structure_id in self.dataset.structure_ids
                    ],
                )
                source = TargetSource(
                    path,
                    target_columns=(target_name,),
                    value_types={target_name: "string"},
                )
                with self.assertRaisesRegex(TargetAttachmentError, "reserved"):
                    split.attach_targets(source)

    def test_premerged_targets_and_frozen_receipts_cannot_cross_releases(self):
        split = data_split(self.classified, group_criteria="none", diversity="none")
        other_dataset, _ = _authenticated_classified(
            self.root / "other",
            cr=10,
            ncr=4,
            metadata_hash="c" * 64,
        )
        premerged = merge_targets(other_dataset, self.source)
        with self.assertRaisesRegex(TargetAttachmentError, "release binding"):
            split.attach_targets(premerged)

        frozen = frozen_assignment_manifest(split.assignment_rows(), split.receipt())
        self.assertIsInstance(frozen, FrozenAssignmentManifest)
        with self.assertRaisesRegex(TargetAttachmentError, "frozen assignment receipt"):
            attach_targets(
                frozen,
                self.source,
                dataset=other_dataset,
            )

    def test_frozen_release_binding_rejects_an_observed_hash_subset(self):
        full_hashes = {
            "metadata/metadata.csv": "b" * 64,
            "parent_groups/parent_groups.csv": "d" * 64,
        }
        full_dataset, full_classified = _authenticated_classified(
            self.root / "full-binding",
            cr=10,
            ncr=4,
            input_hashes=full_hashes,
        )
        split = data_split(full_classified, group_criteria="none", diversity="none")
        frozen = frozen_assignment_manifest(split.assignment_rows(), split.receipt())
        subset_dataset, _ = _authenticated_classified(
            self.root / "subset-binding",
            cr=10,
            ncr=4,
            input_hashes={"metadata/metadata.csv": "b" * 64},
        )
        subset_target = self.root / "subset_targets.csv"
        _write_csv(
            subset_target,
            ("structure_id", "uptake"),
            [
                {"structure_id": value, "uptake": "1.0"}
                for value in subset_dataset.structure_ids
            ],
        )
        with self.assertRaisesRegex(TargetAttachmentError, "frozen assignment receipt"):
            attach_targets(
                frozen,
                TargetSource(
                    subset_target,
                    target_columns=("uptake",),
                    value_types={"uptake": "float"},
                ),
                dataset=subset_dataset,
            )

    def test_historical_split_manifest_receipt_is_verified(self):
        legacy = self.classified.train_valid_test_split(
            parent_method="none",
            leakage_guard="parent_only",
        )
        frozen = frozen_assignment_manifest(legacy.assignment_rows(), legacy.receipt())
        self.assertEqual(
            frozen.assignment_digest, legacy.receipt()["assignment_sha256"]
        )
        self.assertEqual(dict(frozen.assignments), dict(legacy.assignments))

    def test_target_merged_input_is_unwrapped_before_selection_and_receipting(self):
        base_split = data_split(
            self.classified,
            group_criteria="none",
            diversity="none",
            random_state="target-independence",
        )
        merged = merge_targets(self.dataset, self.source)
        labels = dict(self.classified.label_by_id)
        merged_view = _Classified(merged, labels)
        merged_split = data_split(
            merged_view,
            group_criteria="none",
            diversity="none",
            random_state="target-independence",
        )
        self.assertEqual(dict(base_split.assignments), dict(merged_split.assignments))
        self.assertEqual(base_split.assignment_digest, merged_split.assignment_digest)
        self.assertEqual(
            base_split.receipt()["input_hashes"],
            merged_split.receipt()["input_hashes"],
        )
        receipt_text = json.dumps(merged_split.receipt(), sort_keys=True)
        self.assertNotIn(str(self.target_path), receipt_text)
        self.assertEqual(merged_split.receipt()["target_columns_consumed"], [])

    def test_partition_mapping_raw_ids_and_aliases(self):
        aliases_path = self.root / "aliases.csv"
        _write_csv(
            aliases_path,
            ("structure_id", "old_id"),
            [
                {"structure_id": value, "old_id": "OLD-{}".format(index)}
                for index, value in enumerate(self.dataset.structure_ids)
            ],
        )
        alias_target = self.root / "alias_targets.jsonl"
        alias_target.write_text(
            "\n".join(
                json.dumps({"name": "OLD-{}".format(index), "score": index})
                for index in range(len(self.dataset.structure_ids))
            )
            + "\n",
            encoding="utf-8",
        )
        source = TargetSource(
            alias_target,
            id_column="name",
            target_columns=("score",),
            value_types={"score": "int"},
        )
        aliases = AliasRegistry(aliases_path, alias_columns=("old_id",))
        partitioned = attach_targets(
            {
                "train": self.dataset.structure_ids[:2],
                "test": self.dataset.structure_ids[2:4],
            },
            source,
            dataset=self.dataset,
            alias_registry=aliases,
        )
        self.assertEqual(
            partitioned.assignments[self.dataset.structure_ids[0]], "train"
        )
        raw = attach_targets(
            list(self.dataset.structure_ids[:3]),
            self.source,
            dataset=self.dataset,
        )
        self.assertEqual(set(raw.assignments.values()), {"selected"})
        self.assertIsNone(raw.official_split)
        self.assertIsNone(raw.receipt()["official_split"])

    def test_suite_attachment_reports_coverage_per_run_and_partition(self):
        suite = build_cr_ncr_benchmark(
            self.classified,
            ncr_pool_fractions=(0, 0.5, 1),
            seeds=(42, 43),
            group_criteria="none",
            diversity="none",
        )
        original = json.dumps(suite.receipt(), sort_keys=True)
        attached = suite.attach_targets(self.source, missing="keep")
        self.assertIsInstance(attached, TargetAttachedBenchmarkSuite)
        self.assertEqual(len(attached.run_views), 6)
        self.assertEqual(attached.original_assignment_digest, suite.assignment_digest)
        self.assertIs(attached.official_split, False)
        self.assertIs(attached.receipt()["official_split"], False)
        self.assertEqual(json.dumps(suite.receipt(), sort_keys=True), original)
        with tempfile.TemporaryDirectory() as output:
            root = attached.write(output)
            self.assertTrue((root / "attached_membership.csv").is_file())
            self.assertTrue((root / "provenance.jsonl").is_file())
            self.assertTrue((root / "coverage.json").is_file())
            self.assertTrue((root / "receipt.json").is_file())
            for line in (root / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
                digest, name = line.split("  ", 1)
                self.assertEqual(
                    hashlib.sha256((root / name).read_bytes()).hexdigest(),
                    digest,
                )

        dropped = suite.attach_targets(self.source, missing="drop")
        dropped_receipt = dropped.receipt()
        self.assertIn("missing_policy_definition", dropped_receipt)
        self.assertEqual(
            set(dropped_receipt["missing_policy_definitions"]),
            {"keep", "error", "drop"},
        )
        self.assertEqual(set(dropped_receipt["dropped_by_run"]), set(dropped.run_views))
        self.assertTrue(
            any(item["count"] for item in dropped_receipt["dropped_by_run"].values())
        )


if __name__ == "__main__":
    unittest.main()
