"""Illustrative target-first API helpers, not a frozen benchmark rebuilder.

Declare exactly three numeric endpoint columns in targets.json. Full-release
grouping and label-purity checks are performed inside the benchmark API before
eligible_structure_ids is applied. Intersect these IDs with certified model-input
availability before reproducing a common-input benchmark. Never refill a saved
assignment using this example.
"""
import math


def complete_target_ids(db, targets):
    columns = tuple(targets.target_columns)
    if len(columns) != 3:
        raise ValueError("Declare exactly three required numeric target columns")
    if set(db.structure_ids) != set(targets.structure_ids):
        raise ValueError("Target table and release must have the same ID universe")
    def finite(value):
        return (not isinstance(value, bool) and isinstance(value, (int, float))
                and math.isfinite(value))
    return tuple(sid for sid in db.structure_ids
                 if all(finite(targets.target_values(sid)[key]) for key in columns))


options = dict(
    total_size="full_cr",
    train=0.8, val=0.1, test=0.1,
    cohort_eligibility="complete_release_label_pure_effective_blocks",
    diversity="representative",
    test_policy="fixed_pure_cr",
    include_full_cr_diagnostic=False,
)
