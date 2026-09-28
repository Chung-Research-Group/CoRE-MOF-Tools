# Figure, table, and analysis plan

This file organizes manuscript assets; it does not claim that model results
have been generated. The one new figure is an editable workflow schematic.
Existing scientific plots must be selected from their audited snapshot, not
redrawn from a newer live table under an old caption.

## Plotting contract

- Save each scientific panel in its own file. Do not export a multi-panel
  montage as the only PNG. Keep PNG (at least 300 dpi), PDF, and editable-text
  SVG versions when generating numerical figures.
- Prefer about 3.4–4.5 inches wide for a single panel; use up to 6.9 × 4.25
  inches when the relationships need more width. Text must remain at least
  9 pt at the intended final placement size, including ticks and legends.
- Use fixed label colors across the series, clear axis units, concise titles,
  and locally defined categories. Do not use color as the only distinction.
  Keep legends outside dense point clouds when possible.
- Include a data table and caption per figure. State snapshot ID/cutoff,
  checker combination, inclusion/exclusion rules, finite and missing counts,
  counting unit, denominator, feature fields/transforms, and statistical
  settings. Do not place all this information in a tiny title.
- For overlapping subsets, use a count-preserving flow or explicit overlap
  table. A structure is not counted twice within a mutually exclusive flow
  stage. Missing, mixed, and unavailable states must remain visible.

## New workflow figure

![Target-independent assignment and deferred targets](figures/workflow.svg)

**Figure 1. Separation of experiment design from target availability.**
Frozen release IDs, checker evidence, and structural features feed a
complete-release relation graph before filtering. Its connected components
are indivisible partition blocks, not chemical-identity claims. Target-free
diversity can balance distributions but cannot split a block. Only after
cohorts and assignments are frozen are accepted historical and validated
current targets joined by exact ID. Unavailable targets remain null. Model
preprocessing and evaluation are downstream work, not results of this
schematic. The two lanes explicitly show that target availability does not
enter cohort selection. The figure is separately saved as
[editable SVG](figures/workflow.svg), [PDF](figures/workflow.pdf), and
[320-dpi PNG](figures/workflow.png) at 6.9 × 4.25 inches; all visible text is
at least 10.08 pt at that size. Use this overview at its stated full-width
placement; shrinking it requires a new typography/layout check.

## Scientific figures to select or produce

Computation-ready (CR) means all selected checker votes are available and
PASS; non-computation-ready (NCR) means all are available and FAIL. A complete
mixture is AMBIGUOUS; any unavailable vote gives UNCHECKED. Every plot must
name the selected checker combination, not merely say “three checkers.”

For the latest combined reference analysis, `RT` means exact equality of all
264 finite depth-five revised autocorrelation (RAC5) values plus a complete
successful current CrystalNets fingerprint. `M2T` means exact complete
canonical MOFid-v2 text plus that fingerprint. Canonical text conversion
collapses Unicode whitespace, trims, rejects placeholders, applies Unicode
NFKC, then case-folds; it does not change a CIF or its chemistry. MOFid-v2
eligibility is exactly `SUCCESS`, `SUCCESS_TOPOLOGY_UNKNOWN`,
`SUCCESS_TOPOLOGY_ERROR`, or `SUCCESS_TOPOLOGY_TIMEOUT`; every other status
or incomplete input adds no edge. The reference remains provisional with
provisional MOFid evidence. The fingerprint includes complete
SingleNodes/AllNodes subnet status, dimension, key/name/genome and agreement,
with network, subnet/catenation-count and net/agreement summaries. Missing
inputs add no match; RAC5 equality is binary64 exact with signed-zero
canonicalization and no tolerance.

These optional criteria enter neither `priority_main` nor `main_union` by
themselves. `priority_main` is the complete-release conflict-aware explanatory
hierarchy of exact RAC5, then MOFid-v2, then MOFid-v1: lower groups cannot
merge multiple stronger components, conflicts are recorded, missing rows
remain singletons, and it excludes Zeo++, topology, source IDs, CIF hashes,
and StructureMatcher. `main_union` is a separate leakage guard, not a parent
claim; before filtering it forms transitive connected components across the
complete release from exact full CIF SHA-256, database-namespaced source
siblings, and release-authorized RAC5/MOFid-v2/MOFid-v1 groups. The latest
combined analysis adds both optional edge sets to this guard and takes
connected-component closure. State whether an analysis additionally requires
either optional criterion to be available: that explicit filter is not the
general API's default.

| Proposed individual stem | Question and caption requirement | Status |
|---|---|---|
| `target_coverage_by_endpoint` | Finite unique IDs / 42,574, with the September 4 cutoff and endpoint units | Aggregate data available; select snapshot-matched existing plot |
| `checker_composition__<view>__<cohort>` | CR/NCR/AMBIGUOUS/UNCHECKED structure counts; distinguish any-target and all-three-target cohorts | Existing companion analysis; select one file per view/cohort |
| `source_to_checker_to_target__<view>` | Sankey: mutually exclusive source → checker label → target-ready state; all missing states included | Existing companion analysis; verify each stage sums to its declared universe |
| `umap_probe_accessible__<view>__<cohort>` | Probe-accessible-only space; show zero-accessibility and unavailable populations explicitly | Existing frozen embedding; select matching overlay |
| `umap_full_textural__<view>__<cohort>` | Full Zeo++ textural space including density; same cohort colors as accessible-space view | Separate frozen embedding; not interchangeable with previous row |
| `distribution__<feature>__<view>__<cohort>` | One feature per file, common bins/axes, finite n and missing n; distinguish density from counts | Existing companion analysis; choose a compact representative set |
| `source_coverage__<metric>` | OA, SI, COD, modified-CSD and unmodified-CSD subset coverage under the declared relation | Keep exact-ID, blocks-touched, and expanded-ID metrics in separate files |
| `block_label_purity__<view>` | How many structures/blocks are excluded by strict label-pure eligibility? | Existing companion analysis; match grouping policy to caption |
| `paired_error__<endpoint>__<model>` | Error versus actual NCR composition; paired seeds and clean-test n | Pending model training; do not fabricate curves |
| `model_input_coverage__<model>` | Finite label → certified graph/grid → evaluable split, with failure categories | Pending certified preprocessing counts |

An example UMAP caption should spell out Uniform Manifold Approximation and
Projection (UMAP), identify all numerical inputs and transforms, state that
checker labels/targets are overlay annotations only, and say that visual
proximity is not structural identity or quantitative coverage proof.
Do not refit an embedding on each target subset when comparing overlays on
one frozen reference space. If a new fit is necessary, identify it as a new
space and record package versions, seed, neighbors, minimum distance, metric,
thread count, and complete-case exclusions.

The maintained September 3 reference analysis used two different numerical
spaces: eight probe-accessible fields and eighteen full-textural fields.
Those are plotting profiles, not the splitter's 13-plus-two-field fallback.
Retrieve exact field lists, transforms, and excluded-ID ledgers from the
selected figure receipt. Retain zero accessible volume as a scientifically
meaningful state; do not turn it into a missing observation or a positive
logarithmic value by imputation. Do not infer that “OA” is a redistribution
licence, or infer modified-CSD status from a filename.

## Existing companion bundle catalogue

These are recorded project artifact identifiers, not paths in this package
and not a claim of fresh checksum verification:

| Recorded artifact | Intended use | Scope warning |
|---|---|---|
| `published_v2602_checker_combinations_rt_m2t_parent_v3` (2026-09-02) | All 16 checker views and combined-reference grouping dataset | Explicit optional-evidence eligibility; counts differ from raw pools |
| `v2602_latest_group_coverage_zeopp_umaps_v8` (2026-09-03) | Source coverage and both frozen Zeo++ spaces | Full-release analysis, not necessarily a target-ready cohort |
| `v2602_combined_available_targets_20260904_v3` | Canonical historical-plus-current targets at the declared cutoff | Later calculations are absent; private row-level files |
| `target_ready_cr_ncr_analysis_v1` (2026-09-06) | Checker composition, distributions, overlays, Sankey and block-purity views of that target snapshot | Plot creation date is not a new target-data cutoff |

The coordinator's local handoff provides exact paths and private integrity
receipts. Verify those before selecting files. Keep the large plot archive
outside Git and select individual files through the authorized data channel.

For a selected structure set S in release U, with effective blocks B, report
these separately: exact-ID coverage is `|S| / |U|`; block coverage is the
fraction of blocks touched by S; expanded structure coverage is the number of
members in those touched blocks divided by `|U|`. Expansion does not create a
new target value or transfer a source's provenance or redistribution rights.

## Tables

1. Software capability and dependency matrix: manuscript Section 2.
2. Target coverage with units and common denominator: Section 4.1.
3. Raw, excluded, eligible, target-ready, and model-input-ready counts: the
   [evidence table](evidence.md) supplies only the verified stages; model-input
   rows remain pending.
4. Per-run audit table: cohort size, requested pool fraction, actual NCR ratio,
   partition deviations, fixed-test digest, zero crossed blocks, nesting, and
   target coverage. Populate from the selected suite, never from hand-rounded
   intended counts.

Further useful analyses are target-missingness by source/checker/feature tier,
endpoint intersections, within-block label disagreement, and comparisons of
raw versus label-pure populations. These should be descriptive sensitivity
analyses; none alone establishes a causal effect of data curation.
