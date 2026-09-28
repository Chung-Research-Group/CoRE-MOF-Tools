"""Optional integration checks against extracted CoRE-MOF release trees."""

import csv
import gc
import hashlib
import json
import os
from pathlib import Path
import unittest

from CoREMOF.dataset import CoREMOFDataset
from CoREMOF.benchmarks import (
    LABEL_PURE_BLOCK_ELIGIBILITY,
    build_cr_ncr_cohorts,
)


V2601 = os.environ.get("COREMOF_V2601_RELEASE")
V2602 = os.environ.get("COREMOF_V2602_RELEASE")


@unittest.skipUnless(V2601 and V2602, "set both COREMOF_V2601_RELEASE and COREMOF_V2602_RELEASE")
class ReleaseIntegrationTests(unittest.TestCase):
    @staticmethod
    def _csv_digests(path, relative_path):
        result = {}
        with Path(path, relative_path).open(
            "r", encoding="utf-8-sig", newline=""
        ) as handle:
            for row in csv.DictReader(handle):
                structure_id = row["structure_id"]
                payload = json.dumps(
                    row, sort_keys=True, separators=(",", ":")
                ).encode("utf-8")
                result[structure_id] = hashlib.sha256(payload).hexdigest()
        return result

    def test_release_membership_and_shared_metadata(self):
        base = self._csv_digests(Path(V2601), "metadata/metadata.csv")
        self.assertEqual(len(base), 36628)
        superset = self._csv_digests(Path(V2602), "metadata/metadata.csv")
        self.assertEqual(len(superset), 42574)
        self.assertTrue(set(base).issubset(superset))
        self.assertEqual(
            base,
            {structure_id: superset[structure_id] for structure_id in base},
        )
        base_manifest = self._csv_digests(
            Path(V2601), "manifests/cif_manifest.csv"
        )
        superset_manifest = self._csv_digests(
            Path(V2602), "manifests/cif_manifest.csv"
        )
        self.assertEqual(set(base_manifest), set(base))
        self.assertEqual(
            base_manifest,
            {
                structure_id: superset_manifest[structure_id]
                for structure_id in base_manifest
            },
        )

    def test_every_official_checker_view_recomputes_declared_counts(self):
        superset = CoREMOFDataset.from_release(Path(V2602))
        declared = superset.dataset_info["label_counts"]
        for view in ("3checker", "4checker", "5checker"):
            classified = superset.classify(checkers=view)
            observed = dict(classified.label_counts())
            self.assertEqual(observed, declared["label_{}".format(view)])
            del classified
            gc.collect()

    def test_v2602_published_strict_and_label_pure_block_pools(self):
        superset = CoREMOFDataset.from_release(Path(V2602))
        classified = superset.classify(checkers="5checker")
        counts = dict(classified.label_counts())
        self.assertEqual(counts["CR"], 6294)
        self.assertEqual(counts["NCR"], 2299)
        self.assertEqual(counts["AMBIGUOUS"], 7367)
        self.assertEqual(counts["UNCHECKED"], 26614)
        cohorts = build_cr_ncr_cohorts(
            classified,
            cohort_eligibility=LABEL_PURE_BLOCK_ELIGIBILITY,
            diversity="none",
        )
        self.assertEqual(cohorts.raw_pool_counts, {"CR": 6294, "NCR": 2299})
        self.assertEqual(cohorts.pool_counts, {"CR": 4693, "NCR": 1727})
        self.assertTrue(
            cohorts.receipt()["all_zero_partially_selected_effective_blocks"]
        )
        full = [
            cohort
            for cohort in cohorts.cohorts
            if cohort.requested_ncr_pool_fraction == "1.0"
        ]
        self.assertTrue(full)
        self.assertTrue(
            all(set(cohort.ncr_ids) == set(cohorts.full_ncr_ids) for cohort in full)
        )

    def test_open_scope_priority_split_is_complete_deterministic_and_leakage_free(self):
        superset = CoREMOFDataset.from_release(Path(V2602))
        classified = superset.classify(checkers="5checker")
        parameters = {
            "parent_method": "priority_main",
            "leakage_guard": "auto",
            "labels": ("CR", "NCR"),
            "sources": ("COD", "SI"),
            "random_state": "integration-seed",
        }
        first = classified.train_valid_test_split(**parameters)
        first_assignments = dict(first.assignments)
        first_digest = first.receipt()["assignment_sha256"]
        first_conflict_count = len(first.parent_conflicts)
        first_counts = dict(first.counts)
        first_diagnostics = dict(first.parent_diagnostics)
        first_audit = dict(first.leakage_audit)
        first_provisional = first.provisional_input
        first_official = first.official_split
        expected_eligible = {
            record.structure_id
            for record in classified
            if record.label in ("CR", "NCR")
            and record.source_database in ("COD", "SI")
        }
        conflict_member_ids = {
            structure_id
            for conflict in first.parent_conflicts
            for structure_id in conflict["member_ids"]
        }
        del first
        gc.collect()
        second = classified.train_valid_test_split(**parameters)
        self.assertEqual(first_assignments, second.assignments)
        self.assertEqual(first_digest, second.receipt()["assignment_sha256"])
        self.assertEqual(set(first_assignments), expected_eligible)
        self.assertEqual(
            sum(first_counts[split] for split in ("train", "validation", "test")),
            len(expected_eligible),
        )
        self.assertEqual(first_audit["cross_split_block_count"], 0)
        self.assertTrue(first_provisional)
        self.assertFalse(first_official)
        self.assertTrue(
            all(
                value == "PARENT_METHOD_CONFLICT"
                for value in first_diagnostics.values()
            )
        )
        self.assertTrue(set(first_diagnostics).issubset(first_assignments))
        self.assertTrue(set(first_diagnostics).issubset(conflict_member_ids))
        self.assertEqual(
            first_conflict_count, second.receipt()["parent_conflict_count"]
        )


if __name__ == "__main__":
    unittest.main()
