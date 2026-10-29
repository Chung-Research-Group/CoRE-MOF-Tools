"""Synthetic consumer-helper checks, without scientific inputs or calculations."""

import importlib.util
from pathlib import Path
from types import SimpleNamespace
import unittest


_PATH = Path(__file__).resolve().parents[1] / "examples" / "benchmark_inputs.py"
_SPEC = importlib.util.spec_from_file_location("benchmark_inputs_example", _PATH)
example = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(example)


class ThreeEndpointHelperTests(unittest.TestCase):
    def target_table(self, rows, columns=("a", "b", "c"), ids=None):
        return SimpleNamespace(
            structure_ids=tuple(rows) if ids is None else ids,
            target_columns=columns,
            target_values=rows.__getitem__,
        )

    def test_finite_values_preserve_zero_and_release_order(self):
        rows = {
            "synthetic-a": {"a": 0, "b": 1.5, "c": -2},
            "synthetic-b": {"a": 1, "b": 2, "c": 3},
        }
        release = SimpleNamespace(structure_ids=("synthetic-b", "synthetic-a"))
        self.assertEqual(
            example.complete_target_ids(release, self.target_table(rows)),
            release.structure_ids,
        )
        self.assertEqual(rows["synthetic-a"]["a"], 0)

    def test_unavailable_and_nonnumeric_values_are_excluded_without_changes(self):
        values = (None, True, "1.5", float("inf"), float("-inf"), float("nan"))
        for value in values:
            with self.subTest(value=value):
                rows = {"synthetic-a": {"a": value, "b": 2, "c": 3}}
                release = SimpleNamespace(structure_ids=("synthetic-a",))
                self.assertEqual(
                    example.complete_target_ids(release, self.target_table(rows)), (),
                )
                self.assertIs(rows["synthetic-a"]["a"], value)

    def test_exactly_three_target_columns_are_required(self):
        release = SimpleNamespace(structure_ids=("synthetic-a",))
        rows = {"synthetic-a": {"a": 1, "b": 2, "c": 3}}
        for columns in ((), ("a",), ("a", "b"), ("a", "b", "c", "d")):
            with self.subTest(columns=columns), self.assertRaisesRegex(
                ValueError, "exactly three"
            ):
                example.complete_target_ids(release, self.target_table(rows, columns))

    def test_mismatched_id_universe_is_rejected(self):
        release = SimpleNamespace(structure_ids=("synthetic-a",))
        rows = {"synthetic-b": {"a": 1, "b": 2, "c": 3}}
        with self.assertRaisesRegex(ValueError, "same ID universe"):
            example.complete_target_ids(release, self.target_table(rows))

    def test_defaults_are_illustrative_options_not_assignments(self):
        self.assertEqual(example.options, {
            "total_size": "full_cr", "train": 0.8, "val": 0.1, "test": 0.1,
            "cohort_eligibility": "complete_release_label_pure_effective_blocks",
            "diversity": "representative", "test_policy": "fixed_pure_cr",
            "include_full_cr_diagnostic": False,
        })


if __name__ == "__main__":
    unittest.main()
