# Grouped dataset recipes

These copyable recipes use saved release evidence and do not run calculations
or training. Use a new output directory for each analysis. For a new benchmark
that requires complete targets, use the [target-first builder](build_target_first_benchmark.py)
instead. To reproduce the paper dataset, use the [frozen replay](../docs/source/frozen_assignment_replay.rst).

The general recipes below do not select on target availability and do not
regenerate the paper's frozen 3,737/466/468 assignments. Install the pinned
`benchmark` extra for representative diversity. `diversity="none"` is only an
explicit non-representative test or sensitivity profile.

## Checker views

```python
# workflow: checker_views
from itertools import combinations
from CoREMOF.labels import CHECKER_COLUMNS

def checker_views(dataset):
    checkers = tuple(CHECKER_COLUMNS)
    for size in (3, 4, 5):
        for selected in combinations(checkers, size):
            yield selected, dataset.classify(checkers=selected)
```

This yields ten 3-of-5, five 4-of-5 and one 5-of-5 combinations. CR means all
selected votes PASS, NCR all FAIL, a complete mixture AMBIGUOUS, and unavailable
evidence UNCHECKED. The paired constructor uses the named `"5checker"` context.

## General grouped split

```python
# workflow: general_split
from CoREMOF import available_group_criteria

def write_general_split(dataset, output_directory, *, checkers="5checker",
                        group_criteria=("RT", "M2T"),
                        diversity="representative", seed=42):
    availability = available_group_criteria(dataset)
    classified = dataset.classify(checkers=checkers)
    split = classified.data_split(
        group_criteria=group_criteria,
        train=0.8, val=0.1, test=0.1,
        leakage_guard="main_union_plus_criteria",
        diversity=diversity, random_state=seed,
    )
    paths = split.write(output_directory)
    return availability, split, paths
```

`RT` uses exact complete RAC5 and CrystalNets evidence. `M2T` uses the complete
canonical MOFid-v2 identifier and independent CrystalNets evidence. The union
with the full-release `main_union` guard forms indivisible related-structure
groups before filtering. Missing evidence adds no matching edge. See the
[criterion definitions](../README_DATASET_SPLITTING.md).

## Paired target-independent cohorts

```python
# workflow: frozen_benchmark
def freeze_benchmark(dataset, output_directory, *,
                     group_criteria=("RT", "M2T"),
                     diversity="representative"):
    classified = dataset.classify(checkers="5checker")
    cohorts = classified.build_cr_ncr_cohorts(
        ncr_pool_fractions=(0.0, 0.5, 1.0),
        seeds=(912, 913, 914, 915), total_size="full_cr",
        train=0.8, val=0.1, test=0.1,
        group_criteria=group_criteria,
        cohort_eligibility="complete_release_label_pure_effective_blocks",
        diversity=diversity, test_policy="fixed_pure_cr",
    )
    suite = cohorts.data_split(include_full_cr_diagnostic=True)
    suite.write(output_directory)
    return cohorts, suite
```

The explicit label-purity policy excludes an effective group with any other
checker label anywhere in the full release. Inspect raw, excluded and eligible
counts. Increasing `q` adds eligible NCR and removes the same number of CR
structures without moving persistent structures between partitions. Whole-group
constraints may make a requested design infeasible, which is reported rather
than silently repaired. `q=1` includes all eligible NCR, not a pure-NCR cohort.
The common test is pure CR. `full_cr_diagnostic` includes seen structures and
is supplementary, not an independent test. Outputs use `official_split=false`.

The equivalent one-stage CLI is:

```bash
coremof benchmark-cr-ncr /secure/release/CoRE-MOF-COD \
  --group-criteria RT M2T \
  --cohort-eligibility complete_release_label_pure_effective_blocks \
  --ncr-pool-fractions 0.0 0.5 1.0 \
  --seeds 912 913 914 915 --fractions 0.8 0.1 0.1 \
  --output-directory /secure/work/benchmark
```

## Explicit attachment to frozen assignments

```python
# workflow: attach_frozen
def attach_frozen(suite, target_data, output_directory):
    before_digest = suite.assignment_digest
    before_receipt = suite.receipt()
    attached = suite.attach_targets(target_data, missing="keep")
    assert suite.assignment_digest == before_digest
    assert suite.receipt() == before_receipt
    assert attached.original_assignment_digest == before_digest
    attached.write(output_directory)
    return attached
```

Use declared `TargetSource` objects, a verified `targets.json` configuration,
or a premerged table bound to this exact release. `keep` preserves nulls,
`error` requires complete declared targets, and `drop` creates only a filtered
derived view without refill, rebalancing or resplitting.

For an already saved suite, attach one run at a time:

```bash
coremof attach-targets /secure/release/CoRE-MOF-COD \
  --manifest /secure/work/benchmark/coremof_cr_ncr_benchmark/runs/seed912_q0p0.csv \
  --receipt /secure/work/benchmark/coremof_cr_ncr_benchmark/receipt.json \
  --config /secure/targets/targets.json --missing keep \
  --output-directory /secure/work/attached_seed912_q0p0
```

The repeated-ID membership accounting table is not an attachment manifest.
Preserve the original suite and keep attached outputs separately. For the
explicit optional-supplement API, see [target supplements](../docs/source/target_supplements.rst).
