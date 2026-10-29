# CoRE-MOF-COD benchmark handoff

Choose between reproducing frozen assignments and constructing a new benchmark.
The two workflows are not interchangeable. Use the dataset, checker-evidence
view, software identity and input hashes named by the supplied handoff.

## Reproduce the paper's frozen adsorption benchmark

The `coremof_cod_common_input_equal_size_20260912_v1` experiment has twelve
shared assignments: NCR-pool fractions `q=0, 0.5, 1` and seeds `912–915`.
Each case contains **train/val/test: 3,737/466/468** structures. The same
pure-CR test is shared across all cases. `q` is the fraction of the eligible
NCR pool, not the NCR percentage in the training set. Validation includes
CR and NCR when `q>0`.

The frozen design combines `RT` and `M2T` with the full-release leakage guard.
Related-structure groups cannot cross partitions, but every structure retains
its own targets and statistical weight. `RT` requires exact equality of all
264 finite RAC5 descriptors and the complete successful CrystalNets
fingerprint. `M2T` requires the complete canonical MOFid-v2 identifier and that
independent CrystalNets fingerprint. Missing evidence never matches another
missing value. Full criterion definitions are in the
[splitting handbook](README_DATASET_SPLITTING.md).

Use the [frozen-assignment replay](docs/source/frozen_assignment_replay.rst):

```bash
python examples/replay_common_input_benchmark.py \
  --handoff-archive /private/approved_workflow_handoff.tar.gz \
  --archive-sha256 VERIFIED_ARCHIVE_SHA256 \
  --output /private/new-assignment-replay
```

This validates and copies the recorded assignments. It does not resample,
recalculate features or train models. A translated CoRE-ID export has its own
file and assignment hashes. Preserve its verified relationship to the frozen
membership instead of expecting the original bytes to retain their hashes.
New release metadata and newly available targets do not replace the experiment's
bound checker evidence, input certification or assignments.

## Build a new target-complete benchmark

Load target-free core metadata, explicitly merge the requested targets and
select finite-target eligibility before cohort construction. The grouping graph
and checker-label purity still use the complete release, including rows lacking
targets. Target values do not influence grouping, diversity or partitioning.
Model graph, grid and descriptor eligibility requires separate certification.

Use [the target-first guide](docs/source/target_first_benchmark.rst) and its
executable example:

```bash
python examples/build_target_first_benchmark.py /private/validated_release \
  --target-config /private/targets.json \
  --require-target ch4_loading \
  --require-target h2_loading \
  --require-target henry_selectivity \
  --output-directory /private/new-benchmark
```

Replace endpoint names with those in the verified target configuration. The
example selects `RT` and `M2T`, the explicit
`complete_release_label_pure_effective_blocks` policy, representative diversity,
`q=0, 0.5, 1`, seeds `912–915`, and `transition_balanced` partitioning. It writes
eligibility, grouping, target-source and assignment bindings in
`workflow_receipt.json`. These settings do not guarantee the paper's population
sizes on another release, and they do not change legacy package defaults.

Strict CR means all five selected checkers are available and PASS. Strict NCR
means all five are available and FAIL. An unavailable result is a non-vote, not
FAIL, and produces UNCHECKED. A complete mixed set of PASS/FAIL votes is
AMBIGUOUS. The full-release label-purity policy excludes mixed-label groups
before sampling and reports raw, excluded and eligible counts separately.

If the eligible CR pool has size `C` and NCR pool size `M`, the constant-size
cohort uses `round_half_up(q*M)` NCR structures and the remaining `C` rows as
CR. Thus `q=1` means all eligible NCR, not a pure-NCR cohort. Infeasible pool
sizes or whole-group constraints fail explicitly rather than silently capping,
duplicating or splitting a related group. All generated assignments remain
exploratory, with `official_split=false`.

The target-independent builder remains available for applications that do not
select on target availability. Its API and deferred attachment examples are in
the splitting handbook and [grouped recipes](examples/grouped_workflow_recipes.md).
Neither route overrides a supplied frozen experiment.

## Install and verify the selected software

Use the exact source revision or wheel hash recorded with the handoff or release
catalogue. A shared version string alone does not identify an installed build.
Use Python 3.9–3.11. Representative diversity requires the pinned benchmark
extra:

```bash
python -m pip install "/path/to/verified/CoRE-MOF-Tools[benchmark]"
coremof doctor
```

The required numerical stack is NumPy 1.26.4, scikit-learn 1.5.0, SciPy 1.13.1,
joblib 1.5.3 and threadpoolctl 3.6.0. Missing packages and version drift fail
explicitly. Numerical libraries use a one-thread limit, with runtime identity
recorded. Cross-architecture bit identity is not assumed. Compare receipts.

Verify the archive checksum before use, then the ledgers named by its receiver
instructions. A compact metadata-only release can load with
`verify_cif_files=False`, but cannot support graph models without separately
authorized CIFs. Never combine similarly named releases by timestamp.

## Targets, evaluation and transfer

For explicit optional-supplement attachment, use
`dataset.attach_target_supplement(root, expected_sha256=trusted_manifest_hash)`.
See [target supplements](docs/source/target_supplements.rst) and
[combined target construction](COMBINED_TARGET_DATASET.md). Default loading
does not discover or attach targets. Frozen attachment preserves assignments.
`missing="drop"` is only a derived view and never refills or resplits.

Report CH4/H2 uptake in mol/kg-framework and CO2/N2 Henry selectivity as
dimensionless `S`, with the frozen experiment's `log10(1+S)` primary evaluation
and raw scale secondary. Preserve zero and null values. Average seed-level
metrics and report sample standard deviation and actual seed counts. Do not
pool repeated test predictions, substitute screening scores, fill incomplete
cases with another experiment or confuse standardized outputs with evaluated
predictions. A full-CR prediction view containing training-related structures
is diagnostic, not an independent test.

This repository supplies code and documentation, not third-party checker code
or structure-resolved benchmark data. Obtain inputs through the approved
distribution route, applying source-specific permissions and recipient access
requirements. CSD CIFs and other licence-gated payloads do not belong in a
public code repository. A valid checksum does not grant redistribution rights.
See [database access](README_DATABASE_ACCESS.md) and [upload guide](UPLOAD_GUIDE.md).
