# Evidence map and scope of claims

This map binds the draft to source commit
`fad92f7cc8991325aff780da4dbeb6d3b0d537f6` (`0.4.0.dev0`), inspected
2026-09-07. “Implemented” means code exists, not that every optional backend
was executed during this documentation update. Fresh test results and skipped
checks are recorded separately in [verification.json](verification.json).

## Implementation map

| Capability | Source | Maintained verification | Status and caveat |
|---|---|---|---|
| Database lookup and download | [structure.py](../CoREMOF/structure.py) | [test_structure.py](../tests/test_structure.py) | Existing; licensed/asset access and actual remote downloads not revalidated here |
| CIF input resolution and curation | [inputs.py](../CoREMOF/inputs.py), [curate.py](../CoREMOF/curate.py) | [test_inputs.py](../tests/test_inputs.py) | Existing; legacy scientific execution separate from release dataset loading |
| Pore geometry and RAC descriptors | [Zeopp.py](../CoREMOF/calculation/Zeopp.py), [mof_features.py](../CoREMOF/calculation/mof_features.py) | [test_zeopp.py](../tests/test_zeopp.py), [test_racs.py](../tests/test_racs.py) | Existing wrappers; using them does not by itself reproduce the release-pinned methods |
| Property-model wrappers | [prediction.py](../CoREMOF/prediction.py) | [test_predictions.py](../tests/test_predictions.py), [test_cp_utils.py](../tests/test_cp_utils.py) | Existing; no new accuracy or applicability-domain result |
| Release loading and checker consensus | [dataset.py](../CoREMOF/dataset.py), [labels.py](../CoREMOF/labels.py) | [test_dataset_labels.py](../tests/test_dataset_labels.py), [test_split_release_integration.py](../tests/test_split_release_integration.py) | Implemented; unknown statuses fail instead of being guessed |
| Structural relations and historical partitions | [parents.py](../CoREMOF/parents.py), [splitters.py](../CoREMOF/splitters.py) | [test_parents_splitters.py](../tests/test_parents_splitters.py) | Implemented; keep legacy defaults/receipts unchanged |
| Additive splits, numerical strata, paired cohorts | [benchmarks.py](../CoREMOF/benchmarks.py) | [test_benchmarks.py](../tests/test_benchmarks.py) | Implemented; strict five-checker paired constructor, explicit feasibility errors |
| Typed target merging and late attachment | [targets.py](../CoREMOF/targets.py), [attachments.py](../CoREMOF/attachments.py) | [test_targets.py](../tests/test_targets.py), [test_benchmarks.py](../tests/test_benchmarks.py) | Implemented; derived target view retains original assignment digest |
| Combined available-target snapshots | [builder](../examples/build_combined_target_dataset.py), [independent auditor](../examples/audit_combined_target_dataset.py) | [test_combined_target_dataset_examples.py](../tests/test_combined_target_dataset_examples.py) | Implemented; accepted historical values plus fill-only validated current evidence |
| CLI and executable documentation | [cli.py](../CoREMOF/cli.py), [workflow guide](workflows.md) | [test_cli.py](../tests/test_cli.py), [test_notebook.py](../tests/test_notebook.py), [test_manuscript_workflows.py](../tests/test_manuscript_workflows.py) | Implemented; this documentation adds no scientific execution controller |
| Distribution and documentation | [setup.py](../setup.py), [CI](../.github/workflows/tests.yml) | [test_handbook.py](../tests/test_handbook.py) | Development source; CI configuration is not evidence of a fresh release/build pass |

## Count provenance: never merge these denominators

The published-release integration fixture and target aggregate are separate
evidence sources. A fresh read-only count of the published CoREMOF-COD metadata
on 2026-09-07 reconfirmed 42,574 unique IDs, no duplicates, and the four raw
five-checker counts below. It was not a CIF/adapter or whole-cohort audit.
Computation-ready (CR) means all selected checker votes
PASS; non-computation-ready (NCR) means all selected votes FAIL. All votes must
be available. A complete mixed vote is AMBIGUOUS; an unavailable vote gives
UNCHECKED. The named five-checker view uses MOFClassifier, original MOFChecker,
Chen–Manz, MOSAEC, and SETC-GAT.

The recorded benchmark configuration is `priority_main`: the complete-release
conflict-aware explanatory hierarchy of exact RAC5, then MOFid-v2, then
MOFid-v1 groups. Lower groups do not merge multiple stronger components;
conflicts are recorded, missing evidence leaves singletons, and it excludes
Zeo++, topology, source IDs, CIF hashes, and StructureMatcher. Its separate
`main_union` leakage guard is the transitive connected-component closure of
full CIF SHA-256, database-namespaced source siblings, and release-authorized
RAC5/MOFid-v2/MOFid-v1 relations over the complete release before filtering;
it is not a parent claim. Effective blocks add the selected explanatory edges
to that guard. The label-pure policy retains only whole blocks with one strict
label across their complete-release membership.

| Quantity | Value | Denominator / provenance |
|---|---:|---|
| CoREMOF-COD structures | 42,574 | Unique published structure IDs; release integration fixture |
| inherited base cohort / additions | 36,628 / 5,946 | Exact membership-preservation integration assertions |
| Raw strict five-checker CR / NCR | 6,294 / 2,299 | Full-release integration fixture, before optional grouping eligibility |
| Raw AMBIGUOUS / UNCHECKED | 7,367 / 26,614 | Same named five-checker full-release fixture |
| Label-pure eligible CR / NCR | 4,693 / 1,727 | Recorded explanatory-hierarchy configuration above, not arbitrary selected criteria |
| Policy-excluded CR / NCR | 1,601 / 572 | Mixed-label complete-release blocks; original labels unchanged |
| Finite CH4 / H2 / Widom targets | 28,979 / 28,974 / 28,944 | All 42,574 IDs at the September 4 target cutoff |
| Any / all-three finite targets | 28,983 / 28,937 | Union/intersection of endpoint-specific finite-ID sets |
| Finite structure–endpoint assignments | 86,897 | Sum of endpoint counts; not number of structures |
| Simulation-eligible / excluded | 32,988 / 9,586 | Snapshot-specific simulation screening, not strict CR/NCR consensus |

The target source is
[CoREMOF-COD_COMBINED_TARGET_COVERAGE_20260904.json](../CoREMOF-COD_COMBINED_TARGET_COVERAGE_20260904.json),
snapshot `coremof_cod_combined_available_targets_20260904_v3`, cutoff
`2026-09-04T05:43:23Z`. It explicitly records `campaign_complete=false`,
`publication_authorized=false`, and `official_split=false`. The new documentation
does not re-audit private structure-resolved target rows. Source hashes and
restricted transfer ledgers must remain in the private evidence bundle.

## Companion artifacts, not new package outputs from this task

The project workflow maintains separate checker-combination, source-coverage,
and target-ready visualization bundles. Their scope differs from a packaged
API and from a new benchmark generation. The figure plan names the bundles
that an authorized local analyst should inspect; selecting plots does not
authorize copying their structure-resolved data into Git. No existing dataset
or figure is moved, overwritten, or recomputed by this documentation change.

## Claims not supported yet

- A stable 0.4 release, tagged archive, or audited official assignment manifest.
- Completed MatDeepLearn or PMTransformer preprocessing coverage and performance.
- A model accuracy gain, speedup, calibrated interval, or causal effect of NCR
  contamination; each requires a completed, controlled model comparison.
- The same eligible pool counts after changing the grouping criteria or
  source/feature-availability filters.
- Universal duplicate detection, proof of chemical parentage, or absence of
  leakage through relationships missing from the registered evidence.
- Automatic redistribution rights for CSD/SI or licensed checker-derived rows.
- A fully completed adsorption campaign or inclusion of post-cutoff batches
  in the September 4 snapshot.

## Updating this workspace

For a new software baseline, update the source identifier, inspect changed APIs,
run the workflow regression and relevant maintained tests, and create a new
verification record. For new data, freeze a separate snapshot and update all
dependent table captions/denominators together. Do not edit the old target
snapshot or frozen assignments just to make a coverage number larger.
