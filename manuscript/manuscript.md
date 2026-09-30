# CoRE-MOF-Tools: traceable workflows from crystal curation to leakage-aware machine-learning datasets

Author-editable draft · 2026-09-07 · software `0.4.0.dev0`

Authors, affiliations, corresponding author, funding, and journal: to be
confirmed. Source baseline and verification scope are recorded in the
[workspace index](README.md) and [evidence map](evidence.md). This draft does
not announce a release, an official benchmark split, or completed predictive
benchmark results.

## Abstract

Computational studies of metal–organic frameworks require more than a table
of structures and properties: curation choices, unavailable calculations,
related structures, and changing target coverage must remain traceable.
CoRE-MOF-Tools provides a Python interface connecting database access, crystal
curation, structural checks, descriptor calculation, pretrained property
models, and reproducible dataset preparation. Its development version adds a
lightweight release loader, strict checker-consensus classification,
criterion-dependent structural grouping, target-independent partitioning,
and deferred attachment of typed target data. Leakage protection is evaluated
over the complete release before experiment-specific filtering. A paired
benchmark constructor varies the fraction of the eligible non-computation-ready
pool while preserving cohort size and a common computation-ready test set,
subject to indivisible-block feasibility. An application snapshot over 42,574
CoRE-MOF v26.0.2 structures contains finite methane, hydrogen, and raw Widom
ratio targets for 28,979, 28,974, and 28,944 structures, respectively, as of
4 September 2026. These counts describe available targets, not completed
model benchmarks. Machine-readable receipts separate frozen assignments from
later target updates, enabling auditable evaluation as calculation coverage
increases.

Keywords: metal–organic frameworks; scientific software; data provenance;
structural curation; dataset splitting; machine learning.

## 1. Introduction

The CoRE MOF database supports computational screening of experimentally
reported metal–organic frameworks. The 2025 database study combined curation,
calculated and machine-learned properties, and material–process screening.
That study is the scientific context for the software described here; it is
not a publication of the present v26.0.2 snapshot.
[CoRE MOF DB, Matter (2025)](https://doi.org/10.1016/j.matt.2025.102140).

A reusable software workflow must make several distinctions explicit. A
failed calculation is not a negative scientific prediction. Equality under
one structural fingerprint is not proof of chemical identity. A dataset with
an available target is not necessarily a valid model input. Finally,
structures connected through records excluded from a particular experiment
can still create train–test overlap. These distinctions motivate separate
layers for evidence ingestion, classification, structural relationships,
partitioning, target attachment, and model-specific preprocessing.

This work documents the current CoRE-MOF-Tools implementation and its
reproducibility boundaries. The contribution is the integration and audit
contract, not a claim to have invented the external chemistry algorithms or
to have demonstrated superior model performance. The package retains its
historical scientific interfaces while adding dataset APIs that can operate
without importing the scientific calculation stack.

## 2. Software architecture

The distribution name is `CoREMOF_tools`; Python imports use `CoREMOF` and
the command-line entry point is `coremof`. The inspected source version is
`0.4.0.dev0` with Python 3.9–3.11 declared in packaging metadata. The base
installation has no mandatory third-party dependencies. The `full` extra
contains the historical scientific dependency set, whereas `benchmark`
pins the numerical environment used for representative sampling. Installing
an extra does not supply a CSD licence or certify external executables and
model weights.

| Layer | Implemented entry points | Main boundary |
|---|---|---|
| Database and structure access | `structure.information`, SI/CSD download helpers | Licensed and asset-specific access requirements remain explicit |
| CIF curation | `curate.preprocess`, `clean` | May write derived structures; not invoked by the release splitter |
| Saved checker results | `CoREMOFDataset.classify`, `examples/read_checker_results.py` | Reads recorded findings and combines selected votes; external checker execution is not distributed |
| Descriptors and identifiers | Zeo++ wrappers, `mof_features.RACs`, topology/OMS helpers, MOFid wrappers | Method-specific environments and evidence; no implied equivalence between historical and release-pinned runs |
| Pretrained predictions | `prediction.pacman`, `stability`, `cp` | Model/backend availability and domain of applicability require separate verification |
| Release datasets | `CoREMOFDataset.from_release`, `classify`, `filter` | Exact release membership, validated evidence contracts, explicit exclusions |
| Partitions and cohorts | `train_valid_test_split`, `data_split`, `build_cr_ncr_cohorts`, `build_cr_ncr_benchmark` | New target-first cohorts screen target availability before selection; historical deferred-target cohorts retain their frozen assignments |
| Targets | `merge_targets`, `attach_targets`, combined-target example builder/auditor | Exact-ID or declared-alias matching, explicit types/units/conditions |

Zeo++ supplies pore-geometry analysis; revised autocorrelation descriptors
encode graph-based chemical environments; CrystalNets identifies network
topology; and MOFid supplies framework identifiers. The exact depth-five
264-value schema and current-release contracts used here are project-specific
choices, not universal defaults of those tools.
[Zeo++ methods](https://www.maciejharanczyk.info/Zeopp/docs.html),
[RACs](https://doi.org/10.1021/acs.jpca.7b08750),
[CrystalNets](https://doi.org/10.21468/SciPostChem.1.2.005),
[MOFid](https://doi.org/10.1021/acs.cgd.9b01050).

The private release-curation and adsorption controllers are companion project
workflows, not packaged `coremof` job-submission commands. Their independently
validated outputs can become dataset inputs. This manuscript does not treat a
Slurm exit status as proof of scientific target availability.

## 3. Dataset methods

### 3.1 Release universe and strict checker views

The published v26.0.2 membership used by the application contains 42,574
unique structure IDs: the 36,628-member v26.0.1 base and 5,946 additions.
These are released structures, not unique hypothetical chemical parents;
different released solvent-removal variants retain their own IDs. Joins use
`structure_id`, never row position or a filename guessed from an alias.

Computation-ready (CR) and non-computation-ready (NCR) are strict consensus
labels under an explicitly selected checker set. The five canonical checker
names, in order, are MOFClassifier, original MOFChecker, Chen–Manz, MOSAEC,
and SETC-GAT. All selected results must be available and PASS for CR, or
available and FAIL for NCR. A complete mixture of PASS and FAIL is AMBIGUOUS;
any unavailable result makes the view UNCHECKED. Execution failures are
non-votes, not scientific FAIL votes, and majority voting is not used.

Named presets use the first three, first four, or all five checkers.
Explicit ordered combinations also permit all ten three-of-five and five
four-of-five comparisons. The complete combination series therefore has
16 views when the five-of-five view is included. Custom views have distinct
identifiers even when their checker list coincides with a named preset.
The paired contamination constructor currently requires the named, recomputed
five-checker view; support for other views in general splitting does not imply
support in that constructor.

### 3.2 Explanatory relationships and conservative leakage protection

An explanatory relationship and an indivisible partition block answer
different questions. The project-defined `priority_main` is a conflict-aware
explanatory hierarchy constructed across the complete release. Exact
release-authorized depth-five revised autocorrelation (RAC5) groups anchor
first, followed by MOFid-v2 and MOFid-v1 groups. A lower-priority group touching
no stronger component creates a component; one touching exactly one attaches
only unresolved members. A group touching multiple stronger components
records `PARENT_METHOD_CONFLICT` and does not merge them. Missing evidence
leaves structure-specific singletons. This hierarchy excludes Zeo++,
CrystalNets, source IDs, CIF hashes, common names, and StructureMatcher.

Separately, `main_union` is a leakage guard, not a parent or identity claim.
Before filtering the complete release, it takes transitive connected
components over exact full CIF SHA-256, database-namespaced source siblings,
and release-authorized RAC5, MOFid-v2, and MOFid-v1 relations. A missing
optional relation adds no edge; two nulls never match. A missing required CIF
hash fails validation. Building this graph before selecting source, checker
label, or experiment membership preserves connections through hidden bridge
rows.

The additive policy `main_union_plus_criteria` unions that graph with every
ordered user-selected criterion and takes connected-component closure. The
resulting effective leakage blocks cannot cross partitions. This guarantees
zero crossings relative to the registered evidence, not the absence of every
possible unmeasured structural relationship.

Two complementary optional reference criteria support the combined grouping
workflow. `RT` (`rac5_crystalnets`) requires exact equality of all 264 finite
RAC5 values and a complete successful current CrystalNets fingerprint.
`M2T` (`mofid_v2_crystalnets`) requires exact complete canonical MOFid-v2 text
and the same fingerprint. Canonical text conversion collapses Unicode
whitespace, trims, rejects whole-field placeholders, applies Unicode NFKC,
then case-folds; it does not change a CIF or its chemistry. Eligible MOFid-v2
statuses are `SUCCESS`, `SUCCESS_TOPOLOGY_UNKNOWN`,
`SUCCESS_TOPOLOGY_ERROR`, and `SUCCESS_TOPOLOGY_TIMEOUT`; every other status
or incomplete input adds no edge. The latter two statuses describe a
successful identifier with an embedded topology qualifier, not a successful
MOFid execution inferred from a timeout. The criterion remains provisional
when its release-authorized MOFid input is provisional.

The CrystalNets fingerprint includes network dimension, subnet/catenation
counts, SingleNodes and AllNodes net summaries and agreement, and every
complete subnet's status, dimension, topology key/name/genome, and agreement.
It retains duplicate subnets while removing their arbitrary ordering. Exact
RAC5 equality uses finite binary64 values with negative zero mapped to positive
zero and zero numerical tolerance; it uses no rounding or scaling.

Selecting `group_criteria=("RT", "M2T")` adds edges from either available
criterion, not only pairs satisfying both. Neither criterion changes
`priority_main` or the pre-existing `main_union` by itself. A structure
missing both stays eligible by default, adds no optional edge, and may still
be connected by the base leakage guard. An evidence-complete sensitivity
dataset can explicitly exclude it while retaining an exclusion ledger;
unavailability must not be described as proof of uniqueness or as NCR.

### 3.3 Representative target-free diversity

The representative profile first uses a complete finite 264-value RAC5
vector, otherwise a complete vector of 13 contract-selected intensive Zeo++
fields plus channel and bonded-framework dimensions. Structures lacking both
remain in an explicit no-numeric tier. No scientific feature is imputed.
Within each numerical tier, features are centered by the median and divided
by the interquartile range; a zero interquartile range uses a unit divisor,
retaining centered deviations. RAC5 is reduced to at most 32 principal
components with a full singular-value decomposition. A completely zero
scaled matrix has a defined constant-coordinate path.

MiniBatchKMeans strata are generated from sorted IDs with profile seed 2602,
ten initializations, batch size `max(1024, 3*k)`, no cluster reassignment,
and `k = min(n, 256, max(16, ceil(sqrt(n))))`. The requested and effective
number of strata can differ for constant or repeated vectors. Source,
structure variant, current topology category, and feature-availability tier
are balanced alongside these strata. Numerical routines use scikit-learn;
the exact five-package version contract is given in the workflow guide.
[scikit-learn](https://www.jmlr.org/papers/v12/pedregosa11a.html).

These full-universe, target-free strata are part of experiment design, not
model preprocessing or structural identity. They never divide an effective
leakage block. The design uses knowledge of the frozen release's unlabeled
feature distribution and should not be described as a wholly inductive
train-only representation-learning pipeline. Model scalers, feature
selection, and learned encoders must be fit on training data only.

### 3.4 Paired fixed-size CR/NCR cohorts

The paired constructor starts from the unfiltered named five-checker view.
If a strict-label structure shares a complete-release effective block with
another label, the default complete-pool request fails closed. The explicit
sensitivity option `complete_release_label_pure_effective_blocks` retains
only blocks whose complete-release members all have the same strict label.
It changes eligibility, not checker labels. Mixed-label and unavailable-label
blocks are accounted for separately.

Let `C` and `M` be the eligible CR and NCR structure counts under the declared
policy, and `q` the fraction of that eligible NCR pool to include. Each cohort
contains exactly `C` structures, with

```text
NCR_count(q) = round_half_up(q * M)
CR_count(q)  = C - NCR_count(q)
actual_NCR_ratio(q) = NCR_count(q) / C
```

Thus `q=1` uses all eligible NCR structures; it does not mean a 100%-NCR
composition. `M > C` means that the NCR pool is larger than the entire allowed
cohort. Reserving a common CR test further requires `M <= C - test_count`.
The implementation reports actionable counts when these constraints fail.
Whole-block subset feasibility is an additional requirement: satisfying the
count inequalities alone does not guarantee a feasible ladder. Counts are
not silently capped, duplicated, or resized.

Six pool fractions from 0 to 1 in increments of 0.2 and seeds 42–46 give
30 paired runs. A common approximately 10% pure-CR test is selected jointly
with a feasible CR-removal ladder from whole CR-only effective blocks. Its
IDs are identical across every ratio and seed and have zero block overlap
with training or validation. Within each seed, increasing `q` adds nested NCR
memberships and removes reverse-nested CR memberships; persistent structures
retain their partition. Requested train/validation sizes are approximate
because blocks are indivisible, and each deviation is reported. Validation
is composition-balanced with the mixed development set; it is not forced to
be pure CR.

The separate `full_cr_diagnostic` predicts over the entire raw CR pool and
reports exact-ID and same-block overlap with training. It includes seen
structures and is supplementary, not a substitute for the independent test.
All assignments remain exploratory with `official_split=false`.

### 3.5 Frozen assignments and deferred targets

The benchmark workflow freezes structural evidence, diversity, memberships,
test selection, and assignments before opening target values or availability.
`attach_targets` then reuses the existing typed CSV/JSON/JSONL target parser,
exact current IDs, explicitly declared aliases, units, conditions, and
conflict checks. Fuzzy matching and implicit unit conversion are not used.

The default `missing="keep"` is a left join retaining every selected ID and
native null. `missing="error"` validates completeness. `missing="drop"`
creates a derived filtered view without refilling, rebalancing, or resplitting.
Target receipts reference the original assignment digest; target hashes never
alter the frozen split receipt. Persisted manifests are checked against their
paired receipts and the release binding before attachment. Immutable objects
and transactional writers reduce accidental partial or inconsistent exports;
checksums establish integrity, not a substitute for scientific validation or
an adversarial security boundary.

## 4. Application snapshots and verification

### 4.1 Denominator-controlled target coverage

The combined snapshot `v2602_combined_available_targets_20260904_v3` merges
accepted historical results with validated new results at cutoff
`2026-09-04T05:43:23Z`. The denominator below is the complete 42,574-structure
published release, not only eligible simulations or strict checker labels.
Counts are reproduced from the repository's
[public-safe aggregate record](../V2602_COMBINED_TARGET_COVERAGE_20260904.json).

| Endpoint | Finite unique structure IDs | Full-release coverage |
|---|---:|---:|
| CH4 absolute loading, 298 K / 65 bar, mol kg−1 framework | 28,979 | 68.0674% |
| H2 absolute loading, 77 K / 100 bar, mol kg−1 framework | 28,974 | 68.0556% |
| Raw CO2/N2 Widom Rosenbluth-weight ratio, 298 K / 1 bar, dimensionless | 28,944 | 67.9852% |

The snapshot contains 86,897 finite structure–endpoint assignments, 28,983
structures with at least one finite target, and 28,937 with all three. A
historical Widom record with a zero denominator remains explicitly null; the
fill-only merge does not overwrite accepted existing records. Earlier
completion-only counts omitted accepted historical targets and therefore
must not be presented as total available-data coverage. This fixed snapshot
does not include later calculation batches and does not imply campaign
completion or redistribution approval.

### 4.2 Recorded eligibility example, not a new grouping result

The maintained real-release integration fixture records 6,294 raw strict CR,
2,299 raw strict NCR, 7,367 AMBIGUOUS, and 26,614 UNCHECKED structures. Under
the earlier explanatory-hierarchy benchmark configuration, the explicit
label-pure-block policy excludes 1,601 CR and 572 NCR structures, leaving
`C=4,693` and `M=1,727`. Its `q=1` cohort therefore contains 2,966 CR and
1,727 NCR structures, an actual NCR composition of approximately 36.80%.
These policy-specific counts are not asserted for the combined optional
criteria workflow above; that workflow must recompute and receipt its own
eligibility counts. See the [evidence map](evidence.md) for configuration and
test provenance.

A fresh read-only count on 7 September 2026 independently recomputed these
four raw labels from the published metadata and confirmed 42,574 unique IDs
with no duplicates. This count-only check did not rebuild the structural
relations or rerun the full-release cohort integration; the eligible counts
above remain the explicitly identified recorded integration example.

### 4.3 Software verification scope

Maintained tests exercise strict labels, release validation, criterion
aliases, hidden bridge rows, missing evidence, zero block crossings,
whole-block feasibility, nesting, stable paired assignments, target-parser
parity, conflict rejection, immutable receipts, and transactional exports.
The historical splitter has separate compatibility regressions. The
[verification record](verification.json) distinguishes checks rerun for this
documentation update from existing integration assertions and optional checks
not rerun. Passing API tests does not validate every optional chemistry
backend, certify every model input, or establish a predictive advantage.

## 5. Planned evaluation and limitations

The next machine-learning evaluation should consume frozen assignments and
the selected target snapshot, then produce a per-model input-availability
ledger. Finite labels, successful MatDeepLearn graphs, and successful
PMTransformer graphs/grids are different denominators. Their certified input
counts and predictive metrics are pending, not inferred from target coverage.
For each endpoint, report test size after target and model-input availability,
error metrics with units, paired seed differences, and the exact clean-test
intersection used across models. Preserve missing cases; do not fill a gap by
moving another structure into the test set.

Coverage visualization should report exact-ID coverage, effective blocks
touched, and the number of release structures reachable through those blocks
as separate quantities. Companion analyses can overlay checker categories and
source categories on a frozen target-free embedding, but a two-dimensional
map does not establish duplicate identity or complete chemical-space
coverage. Uniform Manifold Approximation and Projection (UMAP) is an optional
visualization in the companion analysis workflow, not the splitter's
clustering backend. Probe-accessible-only and full-textural Zeo++ spaces must
be fitted and reported separately.
[UMAP method](https://arxiv.org/abs/1802.03426).

Strict consensus can exclude useful structures and can disagree across
checker combinations; it is an operational definition rather than ground
truth for all simulations. Missing features do not establish uniqueness.
Conservative transitive leakage blocks can overconnect structures and reduce
sample sizes. The label-pure sensitivity cohort changes the population being
evaluated. Model-specific preprocessing failure and missing targets can
further bias the evaluable subset. These losses require accounting, not
silent data repair. Numerical package pins support reproducibility but do
not promise bit-identical clustering across architectures and numerical
libraries.

## 6. Conclusion

CoRE-MOF-Tools connects existing MOF analysis capabilities with explicit
release validation, checker semantics, structural relationships, and
target-independent dataset design. Its frozen-assignment/late-target
separation supports incremental growth in property coverage without changing
the original experiment. The current application demonstrates traceable data
preparation and available-target coverage; completed predictive comparisons
remain a separate, explicitly pending step.

## Code and data availability draft

Source is maintained at the
[CoRE-MOF-Tools repository](https://github.com/Chung-Research-Group/CoRE-MOF-Tools).
The inspected development commit is recorded in the workspace index; final
submission requires a selected release commit and archival identifier.
Packaging declares CC-BY-4.0 for the project; dependency, model, and data-asset
licences require their own review. Package installation does not confer rights
to CSD data or licensed checker outputs. Structure-resolved target files and
restricted checksum manifests are not included in this documentation update.
COD, SI, and CSD redistribution must follow their separately reviewed asset
permissions. All current benchmark assignments are exploratory, and data
publication approval remains distinct from software availability.

## Author completion items

Confirm authorship and scope; complete model-input and prediction results;
select receipt-matched figures; finish method citations; review availability
and rights statements; freeze the submission source/data versions. Supporting
workflow instructions, captions, evidence mapping, and primary-source links
are supplied as separate files in this workspace.
