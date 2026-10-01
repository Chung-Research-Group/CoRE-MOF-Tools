# Reproducible CoRE-MOF-Tools workflows

Scope: the `0.4.0.dev0` source baseline recorded in [README.md](README.md).
These recipes use the package APIs; they do not submit calculations or rebuild
release evidence. Output directories must be new and appropriately restricted.

For the **current target-first workflow**, use the runnable
[target-first example](../examples/build_target_first_benchmark.py) and
[guide](../docs/source/target_first_benchmark.rst). It merges required targets,
selects finite eligible IDs, constructs full-release related-structure groups
before filtering, and saves a receipt connecting eligibility and assignments.
The 30-run paired construction in sections 5–6 below is retained as the
historical target-independent workflow. It does not reproduce the current
paper's frozen, model-input-certified dataset by itself.

## 1. Choose a workflow and freeze the inputs

For historical structure lookup, curation, descriptors, and pretrained models,
use the [feature guide](../docs/source/features.rst) and
[examples](../examples/README.md). These functions can invoke external tools,
download resources, or write derived CIFs. They are separate from the frozen
release recipes below. Do not run a charge predictor or curation function on
the immutable release CIF path; use an isolated derived-structure workflow.

For dataset work, obtain the exact extracted release selected by the project
handoff, verify its supplied integrity ledger, and record its version. Obtain
target files separately. In the current target-first workflow, required target
availability determines eligibility before cohort selection. In the historical
target-independent workflow below, availability is not a selection criterion.
Target magnitudes never select cohorts or partitions in either workflow.
A compact loader-complete release can omit CIF bytes;
it cannot support a coordinate-based model without the separately authorized
structures. `verify_cif_files=True` requires actual CIF files.

Install the selected source in an isolated Python 3.9–3.11 environment:

```bash
python -m pip install ".[benchmark]"
python -m CoREMOF doctor
```

Run this inside the checkout of the commit named by the handoff. The exact
benchmark stack is NumPy 1.26.4, scikit-learn 1.5.0, SciPy 1.13.1, joblib
1.5.3, and threadpoolctl 3.6.0. Require the representative-backend diagnostic
to pass; no silent numerical fallback is supported. The examples below accept
`diversity="none"` only as an explicit non-representative test/sensitivity
option. Model training should use a separate environment if its requirements
conflict with these pins.

## 2. Define the experiment before using short criterion names

The five checker names are MOFClassifier, original MOFChecker, Chen–Manz,
MOSAEC, and SETC-GAT. Computation-ready (CR) means all selected votes PASS;
non-computation-ready (NCR) means all selected votes FAIL. A complete mixture
is AMBIGUOUS, and any unavailable vote is UNCHECKED. A timeout is not FAIL.

`priority_main` is the complete-release conflict-aware explanatory hierarchy:
exact available depth-five revised autocorrelation (RAC5) groups anchor
first, then MOFid-v2 and MOFid-v1 groups attach unresolved rows. A lower group
touching zero stronger components creates one; touching one attaches unresolved
members; touching several records `PARENT_METHOD_CONFLICT` without merging
them. Missing evidence leaves unique singletons. It excludes Zeo++, topology,
source IDs, CIF hashes, common names, and StructureMatcher.

`main_union` is the separate conservative leakage guard, not a parent claim:
before any filter, it forms transitive connected components over the complete
release from exact full CIF SHA-256, database-namespaced source siblings, and
release-authorized RAC5, MOFid-v2, and MOFid-v1 groups. A missing optional key
adds no edge; missing required CIF hashes fail closed. The additive
`main_union_plus_criteria` policy adds co-membership edges from every chosen
criterion, takes connected-component closure, and keeps the resulting blocks
indivisible. Two missing values never match.

The latest combined reference workflow selects `RT` and `M2T` together.
`RT` means exact equality of all 264 finite RAC5 values plus a complete
successful current CrystalNets fingerprint. `M2T` means exact equality of
complete canonical MOFid-v2 plus that same fingerprint. Canonical text
conversion collapses Unicode whitespace, trims, rejects empty/whole-field
placeholders, applies Unicode NFKC, then case-folds; it does not change a CIF
or its chemistry. MOFid-v2 eligibility is exactly `SUCCESS`,
`SUCCESS_TOPOLOGY_UNKNOWN`, `SUCCESS_TOPOLOGY_ERROR`, or
`SUCCESS_TOPOLOGY_TIMEOUT`; every other status or incomplete input adds no
edge. This reference criterion is provisional if its loaded MOFid evidence is
provisional. The fingerprint includes complete SingleNodes/AllNodes subnet
status, dimension, topology key/name/genome and agreement, with network,
subnet/catenation-count and net/agreement summaries. Missing/error/partial
CrystalNets evidence adds no match. RAC5 equality uses binary64 values with
negative zero mapped to positive zero and zero tolerance, without scaling.

Together these criteria mean **either available relation can add an edge**,
followed by transitive closure with the base guard. They do not require every
structure to have both features, and do not change the definitions of
`priority_main` or `main_union`. A row missing both remains in the general API
by default and may still be connected by the base guard. To reproduce a prior
explicit evidence-complete subset, carry that subset's exclusion manifest;
do not silently adopt it as a general default.

Representative strata balance distributions, not structural identity: complete
264-value RAC5 is used first, otherwise 13 selected intensive Zeo++ fields
plus channel/framework dimensions, otherwise an explicit no-numeric tier.
There is no imputation. The target-free median/interquartile scaling, RAC5
principal-component reduction, and deterministic clustering are defined in
the [manuscript](manuscript.md). Neither these strata nor a plotting embedding
may divide leakage blocks.

## 3. Inspect all 16 checker combinations

Copy the following function into a Python session. It returns each ordered
checker tuple and its classified full-release view, without target filtering.

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

Use `dict(view.label_counts())` for counts and `view.checker_view` for the
exact custom identifier. There are ten 3-of-5, five 4-of-5, and one 5-of-5
combinations. A custom five-checker list is not the named `"5checker"`
authority context required by the paired constructor. For that constructor,
always call `dataset.classify(checkers="5checker")` separately.

## 4. General combined-criterion splits

This function supports any preset or explicit ordered checker combination.
It constructs the full-release guard before selecting CR/NCR rows. The
80/10/10 fractions are requested sizes, not permission to divide blocks.

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

The availability report describes supported criteria and required evidence;
it is not an instruction to silently substitute another criterion. Inspect
the receipt, exclusions, effective blocks, achieved partition counts, and
zero-crossing audit. For several checker views, give each call a distinct
output directory. Independent calls do not promise a common test across views.

## 5. Historical target-independent paired five-checker benchmark

The explicit label-pure sensitivity policy below excludes a whole effective
block if any of its complete-release members has another label. It retains
raw and excluded counts. This is the documented feasible sensitivity design,
not the default complete-strict-pool definition; remove the option only when
requesting a fail-closed audit of that original definition.

```python
# workflow: frozen_benchmark
def freeze_benchmark(dataset, output_directory, *,
                     group_criteria=("RT", "M2T"),
                     diversity="representative"):
    classified = dataset.classify(checkers="5checker")
    cohorts = classified.build_cr_ncr_cohorts(
        ncr_pool_fractions=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
        seeds=(42, 43, 44, 45, 46), total_size="full_cr",
        train=0.8, val=0.1, test=0.1,
        group_criteria=group_criteria,
        cohort_eligibility="complete_release_label_pure_effective_blocks",
        diversity=diversity, test_policy="fixed_pure_cr",
    )
    # Inspect raw_pool_counts, pool_counts, and cohort exclusions here.
    suite = cohorts.data_split(include_full_cr_diagnostic=True)
    suite.write(output_directory)
    return cohorts, suite
```

`C` is the eligible CR count and `M` the eligible NCR count. For NCR-pool
fraction `q`, selected NCR = `round_half_up(q*M)` and selected CR = `C` minus
that number, giving a constant total `C`. Failures such as `M > C`,
`M > C - test_count`, or an infeasible nested whole-block ladder must be
resolved as a new declared experiment, not by capping or duplication. At
`q=1`, all **eligible** NCR are used; the cohort is not necessarily 100% NCR.

The common `fixed_pure_cr` test is pure CR and shared across all 30 runs. The
separate `full_cr_diagnostic` covers the entire raw CR pool, including seen
structures, and is supplementary. Do not reuse the recorded 4,693/1,727
eligible counts from the earlier `priority_main` configuration as the counts
for a newly selected criterion tuple.

Equivalent one-stage CLI, after inspecting the definitions above:

```bash
coremof benchmark-cr-ncr /secure/release/CoREMOF-COD \
  --group-criteria RT M2T \
  --cohort-eligibility complete_release_label_pure_effective_blocks \
  --ncr-pool-fractions 0.0 0.2 0.4 0.6 0.8 1.0 \
  --seeds 42 43 44 45 46 --fractions 0.8 0.1 0.1 \
  --output-directory /secure/work/benchmark
```

The output subdirectory `coremof_cr_ncr_benchmark` contains per-run CSVs,
membership and fixed-test manifests, the full-CR diagnostic, eligibility
accounting, receipt, suite index, and checksum ledger. Use the returned writer
path as authoritative. All outputs have `official_split=false`.

## 6. Historical workflow: attach the combined targets after freezing

The [combined-target guide](../COMBINED_TARGET_DATASET.md) documents deterministic
construction and independent auditing of the historical-plus-current target
snapshot. Its `targets.json` contains explicit types, units, conditions, and
input paths. A newer snapshot is a separate input; it must not overwrite the
old one or change frozen partitions.

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

For this historical workflow, now pass `target_data` as the path to the validated `targets.json`, or as
declared `TargetSource` objects. `keep` preserves all IDs and nulls; `error`
requires completeness; `drop` is a filtered derived view with no refill,
rebalance, or resplit. Store original and attached outputs separately.

For an already saved suite, this CLI attaches **one run** at a time:

```bash
coremof attach-targets /secure/release/CoREMOF-COD \
  --manifest /secure/work/benchmark/coremof_cr_ncr_benchmark/runs/seed42_q0p0.csv \
  --receipt /secure/work/benchmark/coremof_cr_ncr_benchmark/receipt.json \
  --config /secure/targets/targets.json --missing keep \
  --output-directory /secure/work/attached_seed42_q0p0
```

Use a separate output directory for each run. The repeated-ID membership
accounting table is not an attachment manifest. Use `suite.attach_targets`
before serialization to attach all runs in one operation.

## 7. Model benchmarking: planned downstream work

Freeze each model's preprocessing configuration; record success/failure by
structure ID and representation. Report finite targets separately from
MatDeepLearn graph and PMTransformer graph/grid success. Fit model scalers,
feature selection, and learned transformations on training data only.
Keep the same clean-test assignments across models; report actual evaluable
IDs and any common intersection rather than refilling the test.

Report raw target units, endpoint-specific nulls, partition counts, actual NCR
ratios after availability filtering, paired seed results, and exact-ID and
same-block overlap for the supplementary full-CR diagnostic. No training,
preprocessing certification, or prediction result is produced by this recipe.

## 8. Reproduction gates

The runnable Python blocks above are compiled and smoke-tested with synthetic
data by `tests/test_manuscript_workflows.py`. Those small tests explicitly use
non-representative diversity and do not replace the numerical or full-release
integration tests. Consult [verification.json](verification.json) for what was
actually rerun. Before reporting a new production suite, additionally verify
input hashes, numerical pins, whole-block/nesting audits, a common pure-CR
test, exact formula counts, constant assignments within each seed, and target
attachment's unchanged original digest. Transfer restricted files only through
an approved institutional channel; Git contains documentation and aggregate
counts, not the new structure-resolved datasets.
