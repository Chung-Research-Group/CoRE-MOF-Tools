"""Regression tests for eligibility-before-cohort selection."""
import tempfile
import unittest
from pathlib import Path
from test_benchmarks import _authenticated_classified
from CoREMOF.benchmarks import build_cr_ncr_benchmark, BenchmarkError


class TargetFirstTests(unittest.TestCase):
    def build(self, classified, **kwargs):
        return build_cr_ncr_benchmark(
            classified, ncr_pool_fractions=(0, .5, 1), seeds=(912, 913),
            diversity="none", include_full_cr_diagnostic=False, **kwargs)

    def test_missing_target_member_does_not_break_full_release_group(self):
        with tempfile.TemporaryDirectory() as directory:
            same = "a" * 64
            _, classified = _authenticated_classified(
                Path(directory), cif_hashes={"ID-000": same, "ID-001": same,
                                             "ID-002": same})
            ids = tuple(v for v in classified.structure_ids if v != "ID-001")
            result = self.build(classified, eligible_structure_ids=ids)
            self.assertEqual(len(result.eligible_cr_ids), 19)
            for run in result.runs:
                self.assertNotIn("ID-001", run.assignments)
                self.assertEqual(run.assignments.get("ID-000"),
                                 run.assignments.get("ID-002"))
                self.assertTrue(run.leakage_audit["passed"])
            self.assertIn("ID-001", result.effective_leakage_blocks)
            self.assertIn("NOT_IN_ELIGIBLE_STRUCTURE_IDS",
                          result.cohort_eligibility_exclusions["ID-001"])

    def test_order_independence_and_legacy_defaults(self):
        with tempfile.TemporaryDirectory() as directory:
            _, classified = _authenticated_classified(Path(directory))
            baseline = self.build(classified)
            explicit = self.build(classified, eligible_structure_ids=classified.structure_ids)
            self.assertEqual([dict(r.assignments) for r in baseline.runs],
                             [dict(r.assignments) for r in explicit.runs])
            ids = classified.structure_ids[1:]
            first = self.build(classified, eligible_structure_ids=ids)
            second = self.build(classified, eligible_structure_ids=tuple(reversed(ids)))
            self.assertEqual([r.assignment_digest for r in first.runs],
                             [r.assignment_digest for r in second.runs])
            self.assertNotIn("eligibility_filter", baseline.receipt()["cohort_receipt"])

    def test_invalid_ids_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            _, classified = _authenticated_classified(Path(directory))
            for ids in [("not-a-release-id",), ("ID-000", "ID-000"), ()]:
                with self.assertRaises(BenchmarkError):
                    self.build(classified, eligible_structure_ids=ids)
            with self.assertRaises(TypeError):
                self.build(classified, eligible_structure_ids="ID-000")

    def test_filter_cannot_hide_opposite_label(self):
        with tempfile.TemporaryDirectory() as directory:
            same = "a" * 64
            _, classified = _authenticated_classified(Path(directory),
                cif_hashes={"ID-000": same, "ID-020": same})
            ids = tuple(v for v in classified.structure_ids if v != "ID-020")
            result = self.build(classified, eligible_structure_ids=ids,
                cohort_eligibility="complete_release_label_pure_effective_blocks")
            self.assertNotIn("ID-000", result.eligible_cr_ids)
            self.assertIn("BLOCK_NOT_SINGLE_STRICT_LABEL",
                          result.cohort_eligibility_exclusions["ID-000"])

    def test_transition_balancing_keeps_partition_sizes_and_membership(self):
        with tempfile.TemporaryDirectory() as directory:
            _, classified = _authenticated_classified(Path(directory), cr=200, ncr=75)
            legacy = self.build(classified)
            result = self.build(classified, partition_strategy="transition_balanced")
            for old, run in zip(legacy.runs, result.runs):
                self.assertEqual(set(old.assignments), set(run.assignments))
                self.assertEqual(dict(run.achieved_counts),
                                 {"train":160,"validation":20,"test":20,"total":200})
                self.assertTrue(run.leakage_audit["passed"])
                if run.ncr_ids:
                    self.assertGreater(run.label_counts_by_partition["validation"]["NCR"], 0)
            self.assertTrue(result.receipt()["paired_partition_assignment_profile"]
                            ["constant_partition_counts_across_q"])

    def test_transition_balancing_retains_groups_with_missing_members(self):
        with tempfile.TemporaryDirectory() as directory:
            _, classified = _authenticated_classified(Path(directory), cr=80, ncr=24,
                cif_hashes={"ID-000":"a"*64,"ID-001":"a"*64,"ID-002":"a"*64,
                            "ID-081":"b"*64,"ID-082":"b"*64})
            ids=tuple(v for v in classified.structure_ids if v != "ID-001")
            result=self.build(classified, eligible_structure_ids=ids,
                              partition_strategy="transition_balanced")
            for run in result.runs:
                self.assertEqual(run.assignments.get("ID-000"), run.assignments.get("ID-002"))
                self.assertEqual(run.assignments.get("ID-081"), run.assignments.get("ID-082"))
                self.assertNotIn("ID-001", run.assignments)
            for seed in result.seeds:
                self.assertEqual(len({tuple(sorted(r.achieved_counts.items()))
                    for r in result.runs if r.seed==seed}),1)


if __name__ == "__main__":
    unittest.main()
