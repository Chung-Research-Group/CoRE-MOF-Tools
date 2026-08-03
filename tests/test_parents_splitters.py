import json
import hashlib
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from CoREMOF.parents import ParentResolver
from CoREMOF.splitters import (
    OfficialSplitUnavailableError,
    ParentGroupSplitter,
    split_release,
)


class _Dataset:
    def __init__(self, rows, parents=None, hashes=None):
        self.metadata_rows = tuple(rows)
        self.parent_by_id = parents or {}
        if hashes is None:
            self.cif_hashes = {
                row["structure_id"]: hashlib.sha256(
                    row["structure_id"].encode("utf-8")
                ).hexdigest()
                for row in self.metadata_rows
            }
        else:
            self.cif_hashes = hashes
        self.dataset_version = "v-test"
        self.input_hashes = {"metadata.csv": "abc123"}
        self.parent_group_methods = {
            "release_status": "PROVISIONAL_LATEST_AUDITED_SNAPSHOT"
        }


class _Classified:
    def __init__(self, dataset, labels):
        self.dataset = dataset
        self.label_by_id = labels
        self.checker_view = "3checker"


def _rows(ids):
    return [
        {
            "structure_id": structure_id,
            "source_database": "COD",
            "source_id": "SRC-%s" % structure_id,
            "structure_variant": "ASR",
            "metal_elements": "Cu;Zn" if index % 2 else "Cu",
        }
        for index, structure_id in enumerate(ids)
    ]


def _unique_hashes(ids):
    return {
        structure_id: hashlib.sha256(structure_id.encode("utf-8")).hexdigest()
        for structure_id in ids
    }


class ParentResolverTests(unittest.TestCase):
    def test_missing_direct_value_is_unique_singleton_by_default(self):
        dataset = _Dataset(
            _rows(("a", "b", "c")),
            {"a": {"rac5": "R1"}, "b": {"rac5": "R1"}, "c": {}},
        )
        resolution = ParentResolver(dataset).resolve("rac5")
        self.assertEqual(resolution.groups["a"], resolution.groups["b"])
        self.assertEqual(resolution.groups["c"], "SINGLETON:c")
        self.assertFalse(resolution.exclusions)
        with self.assertRaises(TypeError):
            resolution.groups["a"] = "changed"

        strict = ParentResolver(dataset, missing_parent="exclude").resolve("rac5")
        self.assertEqual(strict.exclusions["c"], "MISSING_PARENT_EVIDENCE")

    def test_release_not_available_group_is_not_treated_as_evidence(self):
        dataset = _Dataset(
            _rows(("missing", "x", "y")),
            {
                "missing": {
                    "rac_group": "R-SYNTHETIC",
                    "rac_status": "NOT_AVAILABLE",
                    "rac_size": "1",
                },
                "x": {"rac_group": "R1", "rac_status": "MATCHED", "rac_size": "2"},
                "y": {"rac_group": "R1", "rac_status": "MATCHED", "rac_size": "2"},
            },
        )
        resolution = ParentResolver(dataset).resolve("rac5")
        self.assertEqual(resolution.groups["missing"], "SINGLETON:missing")
        self.assertEqual(resolution.evidence_by_id["missing"], ())
        self.assertEqual(resolution.groups["x"], "R1")

    def test_release_status_and_size_are_closed_and_validated(self):
        bad_status = _Dataset(
            _rows(("a",)),
            {"a": {"rac_group": "R1", "rac_status": "AVAILABLE", "rac_size": "1"}},
        )
        with self.assertRaisesRegex(ValueError, "MATCHED, UNMATCHED, or NOT_AVAILABLE"):
            ParentResolver(bad_status).resolve("rac5")
        bad_size = _Dataset(
            _rows(("a",)),
            {"a": {"rac_group": "R1", "rac_status": "MATCHED", "rac_size": "1"}},
        )
        with self.assertRaisesRegex(ValueError, "size >= 2"):
            ParentResolver(bad_size).resolve("rac5")

    def test_priority_conflict_preserves_anchors_and_excludes_lower_only_row(self):
        ids = ("a", "b", "c", "d", "e", "f", "g")
        parents = {
            "a": {"rac5": "R1", "mofid_v2": "X"},
            "b": {"rac5": "R2", "mofid_v2": "X"},
            "c": {"mofid_v2": "X", "mofid_v1": "Z"},
            "d": {"mofid_v2": "Y", "mofid_v1": "Q"},
            "e": {"mofid_v2": "Y"},
            "f": {"mofid_v1": "Q"},
            "g": {},
        }
        resolution = ParentResolver(_Dataset(_rows(ids), parents)).resolve("priority_main")
        self.assertNotEqual(resolution.groups["a"], resolution.groups["b"])
        self.assertEqual(resolution.exclusions["c"], "PARENT_METHOD_CONFLICT")
        self.assertNotIn("c", resolution.groups)
        self.assertEqual(resolution.groups["d"], resolution.groups["e"])
        self.assertEqual(resolution.groups["d"], resolution.groups["f"])
        self.assertEqual(resolution.groups["g"], "SINGLETON:g")
        self.assertEqual(len(resolution.conflicts), 1)
        conflict = resolution.conflicts[0]
        self.assertEqual(conflict.lower_method, "mofid_v2")
        self.assertEqual(conflict.lower_group, "X")
        self.assertEqual(conflict.stronger_components, ("RAC5:R1", "RAC5:R2"))
        self.assertEqual(conflict.member_ids, ("a", "b", "c"))
        self.assertEqual(conflict.unresolved_ids, ("c",))

    def test_priority_records_anchor_only_conflicts_without_excluding_rows(self):
        parents = {
            "a": {"rac5": "R1", "mofid_v2": "M-X"},
            "b": {"rac5": "R2", "mofid_v2": "M-X"},
        }
        resolution = ParentResolver(_Dataset(_rows(("a", "b")), parents)).resolve(
            "priority_main"
        )
        self.assertFalse(resolution.exclusions)
        self.assertEqual(len(resolution.conflicts), 1)
        conflict = resolution.conflicts[0]
        self.assertEqual(conflict.unresolved_ids, ())
        self.assertEqual(
            conflict.component_members,
            (("RAC5:R1", ("a",)), ("RAC5:R2", ("b",))),
        )

    def test_main_union_keeps_transitive_full_universe_bridges(self):
        ids = ("a", "b", "c", "d", "e", "f")
        rows = _rows(ids)
        # b-c is a source bridge; c-d is a RAC bridge; d-e is MOFid v2.
        for row in rows:
            if row["structure_id"] in ("b", "c"):
                row["source_id"] = "SHARED-SOURCE"
        parents = {
            "c": {"rac5": "R-CD"},
            "d": {"rac5": "R-CD", "mofid_v2": "M-DE"},
            "e": {"mofid_v2": "M-DE"},
        }
        hashes = _unique_hashes(ids)
        hashes["a"] = hashes["b"] = "a" * 64
        resolution = ParentResolver(_Dataset(rows, parents, hashes)).resolve("main_union")
        self.assertEqual(len({resolution.groups[value] for value in ("a", "b", "c", "d", "e")}), 1)
        self.assertNotEqual(resolution.groups["f"], resolution.groups["a"])

    def test_main_union_rejects_truncated_or_placeholder_cif_hashes(self):
        dataset = _Dataset(_rows(("a",)), hashes={"a": "abc123"})
        with self.assertRaisesRegex(ValueError, "full lowercase SHA-256"):
            ParentResolver(dataset).resolve("main_union")

    def test_main_union_requires_a_hash_for_every_structure(self):
        dataset = _Dataset(_rows(("a", "b")), hashes={})
        with self.assertRaisesRegex(ValueError, "requires a full CIF SHA-256"):
            ParentResolver(dataset).resolve("main_union")

    def test_not_available_source_group_is_not_resurrected_from_metadata(self):
        rows = _rows(("a", "b"))
        for row in rows:
            row["source_id"] = "SAME-DISPLAY-ID"
        parents = {
            value: {
                "source_group": "S-SYNTHETIC-%s" % value,
                "source_status": "NOT_AVAILABLE",
                "source_size": "1",
            }
            for value in ("a", "b")
        }
        resolution = ParentResolver(_Dataset(rows, parents)).resolve("main_union")
        self.assertNotEqual(resolution.groups["a"], resolution.groups["b"])

    def test_mofid_v1_conflict_cannot_merge_two_mofid_v2_components(self):
        parents = {
            "a": {"mofid_v2": "M2-A", "mofid_v1": "M1-X"},
            "b": {"mofid_v2": "M2-B", "mofid_v1": "M1-X"},
            "c": {"mofid_v1": "M1-X"},
        }
        resolution = ParentResolver(_Dataset(_rows(("a", "b", "c")), parents)).resolve(
            "priority_main"
        )
        self.assertNotEqual(resolution.groups["a"], resolution.groups["b"])
        self.assertEqual(resolution.exclusions["c"], "PARENT_METHOD_CONFLICT")


class ParentGroupSplitterTests(unittest.TestCase):
    def _singleton_classified(self, count=100, order=None):
        ids = ["s%03d" % index for index in range(count)]
        if order is not None:
            ids = list(order)
        dataset = _Dataset(_rows(ids))
        labels = {
            structure_id: ("CR" if int(structure_id[1:]) % 2 else "NCR")
            for structure_id in ids
        }
        return _Classified(dataset, labels)

    def test_singleton_split_hits_requested_fractions_and_is_unpackable(self):
        result = ParentGroupSplitter(
            self._singleton_classified(), parent_method="none", random_state="stable-seed"
        ).train_valid_test_split((0.6, 0.2, 0.2))
        self.assertEqual(result.counts["train"], 60)
        self.assertEqual(result.counts["validation"], 20)
        self.assertEqual(result.counts["test"], 20)
        train, valid, test = result
        self.assertEqual(train, result.train_indices)
        self.assertEqual(valid, result.valid_indices)
        self.assertEqual(test, result.test_indices)
        self.assertEqual(result.leakage_audit["cross_split_block_count"], 0)
        self.assertEqual(
            result.achieved_fractions,
            {"train": 0.6, "validation": 0.2, "test": 0.2},
        )
        self.assertEqual(result.label_counts_by_split["train"], {"CR": 30, "NCR": 30})

    def test_fraction_contract_accepts_zero_and_rejects_invalid_values(self):
        splitter = ParentGroupSplitter(
            self._singleton_classified(count=10), parent_method="none"
        )
        all_train = splitter.train_valid_test_split((1.0, 0.0, 0.0))
        self.assertEqual(all_train.counts["train"], 10)
        self.assertEqual(all_train.counts["validation"], 0)
        for fractions in (
            (0.8, 0.2),
            (0.8, 0.3, -0.1),
            (0.8, 0.3, 0.0),
            (float("nan"), 0.5, 0.5),
            (0.0, 0.0, 0.0),
        ):
            with self.subTest(fractions=fractions):
                with self.assertRaises(ValueError):
                    splitter.train_valid_test_split(fractions)

    def test_assignment_is_independent_of_metadata_row_order(self):
        ids = ["s%03d" % index for index in range(80)]
        shuffled = list(ids)
        random.Random(991).shuffle(shuffled)
        first = ParentGroupSplitter(
            self._singleton_classified(order=ids),
            parent_method="none",
            random_state="order-test",
        ).train_valid_test_split()
        second = ParentGroupSplitter(
            self._singleton_classified(order=shuffled),
            parent_method="none",
            random_state="order-test",
        ).train_valid_test_split()
        self.assertEqual(first.assignments, second.assignments)
        self.assertEqual(
            first.receipt()["assignment_sha256"],
            second.receipt()["assignment_sha256"],
        )

    def test_unknown_stratification_field_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "Unknown stratify_by"):
            ParentGroupSplitter(
                self._singleton_classified(),
                parent_method="none",
                stratify_by=("label", "sorce_database"),
            )

    def test_external_label_maps_use_the_closed_canonical_vocabulary(self):
        classified = self._singleton_classified(count=4)
        classified.label_by_id["s000"] = "MAYBE"
        with self.assertRaisesRegex(ValueError, "unknown classification label"):
            ParentGroupSplitter(classified, parent_method="none")

    def test_main_union_is_a_guard_not_an_explanatory_parent_method(self):
        classified = self._singleton_classified(count=4)
        with self.assertRaisesRegex(TypeError, "parent method must be a string"):
            ParentResolver(classified.dataset).resolve(None)
        with self.assertRaisesRegex(TypeError, "parent method must be a string"):
            ParentGroupSplitter(classified, parent_method=None)
        with self.assertRaisesRegex(ValueError, "Unknown parent method"):
            ParentGroupSplitter(
                classified, parent_method="main_union"
            )

    def test_explicit_structure_id_filter_is_reported(self):
        classified = self._singleton_classified(count=10)
        result = ParentGroupSplitter(
            classified,
            parent_method="none",
            structure_ids=("s001", "s003", "s005"),
        ).train_valid_test_split((1.0, 0.0, 0.0))
        self.assertEqual(set(result.train_ids), {"s001", "s003", "s005"})
        self.assertEqual(result.exclusions["s000"], "STRUCTURE_ID_FILTER")
        with self.assertRaisesRegex(KeyError, "Unknown structure_id"):
            ParentGroupSplitter(
                classified,
                parent_method="none",
                structure_ids=("absent",),
            )

    def test_assignment_is_independent_of_python_hash_seed(self):
        script = r'''
import json
from types import SimpleNamespace
from CoREMOF.splitters import ParentGroupSplitter
ids = ["s%03d" % i for i in range(40)]
rows = [{"structure_id": x, "source_database": "COD", "source_id": x,
         "structure_variant": "ASR", "metal_elements": "Cu"} for x in ids]
dataset = SimpleNamespace(metadata_rows=rows, parent_by_id={}, cif_hashes={},
                          dataset_version="v", input_hashes={}, parent_group_methods={})
classified = SimpleNamespace(dataset=dataset,
    label_by_id={x: ("CR" if i % 2 else "NCR") for i, x in enumerate(ids)},
    checker_view="3checker")
result = ParentGroupSplitter(classified, parent_method="none",
                             random_state="hash-seed").train_valid_test_split()
print(json.dumps(dict(result.assignments), sort_keys=True))
'''
        outputs = []
        for seed in ("1", "918273"):
            environment = dict(os.environ)
            environment["PYTHONHASHSEED"] = seed
            outputs.append(
                subprocess.check_output(
                    [sys.executable, "-c", script],
                    cwd=str(Path(__file__).resolve().parents[1]),
                    env=environment,
                    text=True,
                )
            )
        self.assertEqual(outputs[0], outputs[1])

    def test_main_union_guard_blocks_mofid_bridge_across_rac_anchors(self):
        ids = ("a", "b", "c", "d", "e", "f")
        parents = {
            "a": {"rac5": "R1", "mofid_v2": "M-SHARED"},
            "b": {"rac5": "R2", "mofid_v2": "M-SHARED"},
            "c": {"rac5": "R3"},
            "d": {"rac5": "R4"},
            "e": {"rac5": "R5"},
            "f": {"rac5": "R6"},
        }
        classified = _Classified(_Dataset(_rows(ids), parents), {value: "CR" for value in ids})
        result = ParentGroupSplitter(
            classified,
            parent_method="rac5",
            leakage_guard="main_union",
            random_state=7,
        ).train_valid_test_split((0.5, 0.25, 0.25))
        self.assertEqual(result.assignments["a"], result.assignments["b"])
        self.assertEqual(result.leakage_groups["a"], result.leakage_groups["b"])
        self.assertTrue(result.leakage_audit["passed"])

    def test_priority_conflict_row_is_retained_when_main_union_guards_it(self):
        ids = ("a", "b", "bridge", "d", "e")
        parents = {
            "a": {"rac5": "R1", "mofid_v2": "M-SHARED"},
            "b": {"rac5": "R2", "mofid_v2": "M-SHARED"},
            "bridge": {"mofid_v2": "M-SHARED"},
            "d": {"rac5": "R3"},
            "e": {"rac5": "R4"},
        }
        classified = _Classified(
            _Dataset(_rows(ids), parents), {value: "CR" for value in ids}
        )
        result = ParentGroupSplitter(
            classified,
            parent_method="priority_main",
            leakage_guard="main_union",
        ).train_valid_test_split((0.6, 0.2, 0.2))
        self.assertIn("bridge", result.assignments)
        self.assertEqual(
            result.parent_diagnostics["bridge"], "PARENT_METHOD_CONFLICT"
        )
        self.assertNotIn("bridge", result.exclusions)
        self.assertEqual(result.assignments["a"], result.assignments["bridge"])
        self.assertIn("PARENT_METHOD_CONFLICTS_GUARDED", result.warnings)
        self.assertIn("PARENT_METHOD_CONFLICT_LEDGER_PRESENT", result.warnings)
        self.assertEqual(result.receipt()["parent_conflict_count"], 1)
        ledger = result.receipt()["parent_conflicts"][0]
        self.assertEqual(ledger["lower_method"], "mofid_v2")
        self.assertEqual(ledger["unresolved_ids"], ["bridge"])

    def test_filtered_middle_row_still_bridges_selected_rows(self):
        ids = ("a", "bridge", "c", "d", "e", "f")
        parents = {
            "bridge": {"mofid_v2": "M-BC"},
            "c": {"mofid_v2": "M-BC"},
        }
        hashes = _unique_hashes(ids)
        hashes["a"] = hashes["bridge"] = "b" * 64
        labels = {value: "CR" for value in ids}
        labels["bridge"] = "AMBIGUOUS"
        classified = _Classified(_Dataset(_rows(ids), parents, hashes), labels)
        result = ParentGroupSplitter(
            classified,
            parent_method="priority_main",
            leakage_guard="main_union",
            labels=("CR",),
        ).train_valid_test_split((0.5, 0.25, 0.25))
        self.assertEqual(result.assignments["a"], result.assignments["c"])
        self.assertEqual(result.exclusions["bridge"], "LABEL_FILTER")

    def test_prefiltered_classification_absence_is_an_explicit_exclusion(self):
        dataset = _Dataset(_rows(("a", "b", "c")))
        classified = _Classified(dataset, {"a": "CR", "b": "NCR"})
        result = split_release(classified, parent_method="none")
        self.assertEqual(result.exclusions["c"], "PRESELECTION_FILTER")
        self.assertTrue(result.filters["preselection"]["active"])
        self.assertTrue(
            all(
                index < 2
                for index in (
                    *result.train_indices,
                    *result.valid_indices,
                    *result.test_indices,
                )
            )
        )
        self.assertEqual(result.view_index_by_id, {"a": 0, "b": 1})
        self.assertEqual(result.index_by_id["c"], 2)
        self.assertEqual(result.filters["preselection"]["selected_count"], 2)

    def test_split_release_honours_or_rejects_explicit_preclassified_view(self):
        classified = self._singleton_classified(count=10)
        same = split_release(classified, checkers="3checker", parent_method="none")
        self.assertEqual(same.checker_view, "3checker")
        with self.assertRaisesRegex(ValueError, "preclassified dataset uses"):
            split_release(classified, checkers="5checker", parent_method="none")

    def test_split_release_loads_classifies_and_splits_a_path(self):
        from tests.test_dataset_labels import _make_release

        with tempfile.TemporaryDirectory() as directory:
            release_root = Path(directory, "release")
            release_root.mkdir()
            _make_release(release_root)
            result = split_release(
                release_root,
                checkers=3,
                parent_method="none",
                fractions=(0.5, 0.25, 0.25),
            )
        self.assertEqual(result.dataset_version, "vtest")
        self.assertEqual(result.checker_view, "3checker")
        self.assertEqual(sum(result.counts[name] for name in ("train", "validation", "test")), 2)

    def test_receipt_is_atomic_provisional_and_refuses_overwrite(self):
        result = ParentGroupSplitter(
            self._singleton_classified(count=20), parent_method="none"
        ).train_valid_test_split()
        with tempfile.TemporaryDirectory() as directory:
            csv_path, json_path = result.write(directory, stem="example")
            self.assertTrue(csv_path.exists())
            receipt = json.loads(Path(json_path).read_text(encoding="utf-8"))
            self.assertTrue(receipt["provisional_input"])
            self.assertFalse(receipt["official_split"])
            self.assertTrue(receipt["leakage_audit"]["passed"])
            self.assertEqual(receipt["algorithm"]["version"], "1.0")
            self.assertEqual(receipt["implementation"]["split_api_version"], "0.1.0")
            self.assertFalse(receipt["integrity"]["cif_files_verified"])
            self.assertEqual(set(receipt["implementation"]["source_sha256"]), {
                "dataset.py", "labels.py", "parents.py", "splitters.py"
            })
            self.assertIn("PROVISIONAL_PARENT_INPUT", receipt["warnings"])
            self.assertIn("achieved_fractions", receipt)
            self.assertIn("label_counts_by_split", receipt)
            self.assertEqual(receipt["parent_conflict_count"], 0)
            with self.assertRaises(FileExistsError):
                result.write(directory, stem="example")
            result.write(directory, stem="example", overwrite=True)
            with self.assertRaisesRegex(ValueError, "filename stem"):
                result.write(directory, stem="../escape")
        with self.assertRaises(TypeError):
            result.assignments["s000"] = "test"
        with self.assertRaises(TypeError):
            result.filters["labels"] = ["AMBIGUOUS"]

    def test_receipt_freezes_implementation_hashes_at_split_construction(self):
        first_hashes = {name: "a" * 64 for name in (
            "dataset.py", "labels.py", "parents.py", "splitters.py"
        )}
        second_hashes = {name: "b" * 64 for name in first_hashes}
        with patch("CoREMOF.splitters._implementation_hashes", return_value=first_hashes):
            result = ParentGroupSplitter(
                self._singleton_classified(count=10), parent_method="none"
            ).train_valid_test_split()
        with patch("CoREMOF.splitters._implementation_hashes", return_value=second_hashes):
            receipt = result.receipt()
        self.assertEqual(receipt["implementation"]["source_sha256"], first_hashes)

    def test_pair_write_rolls_back_if_second_render_fails(self):
        result = ParentGroupSplitter(
            self._singleton_classified(count=20), parent_method="none"
        ).train_valid_test_split()
        with tempfile.TemporaryDirectory() as directory:
            original_to_json = result.to_json
            try:
                object.__setattr__(
                    result,
                    "to_json",
                    lambda *args, **kwargs: (_ for _ in ()).throw(
                        RuntimeError("synthetic render failure")
                    ),
                )
                with self.assertRaisesRegex(RuntimeError, "synthetic render failure"):
                    result.write(directory, stem="pair")
            finally:
                object.__setattr__(result, "to_json", original_to_json)
            self.assertFalse(Path(directory, "pair.csv").exists())
            self.assertFalse(Path(directory, "pair.json").exists())

    def test_official_split_request_fails_closed(self):
        with self.assertRaisesRegex(
            OfficialSplitUnavailableError, "No audited official split manifest"
        ):
            ParentGroupSplitter(
                self._singleton_classified(count=10),
                parent_method="none",
                official=True,
            )

    def test_custom_checker_view_is_marked_user_defined_in_receipt(self):
        classified = self._singleton_classified(count=10)
        classified.checker_view = "custom:MOFClassifier+MOSAEC"
        classified.checker_view_official = False
        result = ParentGroupSplitter(
            classified, parent_method="none"
        ).train_valid_test_split()
        self.assertEqual(
            result.receipt()["checker_view_kind"], "USER_DEFINED"
        )

    def test_final_input_state_uses_a_closed_two_artifact_contract(self):
        classified = self._singleton_classified(count=6)
        classified.dataset.parent_group_methods = {"release_status": "FINAL_CANDIDATE"}
        classified.dataset.dataset_info = {"release_status": "FINAL"}
        candidate = ParentGroupSplitter(
            classified, parent_method="none"
        ).train_valid_test_split()
        self.assertTrue(candidate.provisional_input)

        classified.dataset.parent_group_methods = {"release_status": "FINAL"}
        final = ParentGroupSplitter(
            classified, parent_method="none"
        ).train_valid_test_split()
        self.assertFalse(final.provisional_input)


if __name__ == "__main__":
    unittest.main()
