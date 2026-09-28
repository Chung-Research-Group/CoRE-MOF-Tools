# CoRE-MOF-Tools manuscript and workflow workspace

This is an editable, journal-neutral software-paper draft, with methods that
can also be reused in the main CoRE-MOF database paper. It documents the
implemented checkout, not an announced software release or a completed model
comparison.

Source baseline: `0.4.0.dev0`, commit
`fad92f7cc8991325aff780da4dbeb6d3b0d537f6`, branch
`agent/dataset-splitting-api`, inspected on 2026-09-07. The source worktree was
clean before this documentation work. Database version `v26.0.2` is a separate
identifier from the software version. These documentation additions follow the
baseline commit; they are not represented as already committed or published.

## Reading order

| File | Purpose |
|---|---|
| [manuscript.md](manuscript.md) | Title, abstract, introduction, implementation, methods, example results, limitations, and availability statements |
| [workflows.md](workflows.md) | Practical Python/CLI workflow, checker combinations, assignment freeze, and later target attachment |
| [evidence.md](evidence.md) | Claim-to-code/test map, count denominators, and implemented versus pending work |
| [figures_and_tables.md](figures_and_tables.md) | Individual figure plan, captions, plotting requirements, and existing-artifact selection rules |
| [references.md](references.md) | Primary-source references and citation tasks still requiring author review |
| [figures/workflow.svg](figures/workflow.svg), [PDF](figures/workflow.pdf), [PNG](figures/workflow.png) | Editable/vector/320-dpi workflow schematic; no scientific measurements are encoded |
| [verification.json](verification.json) | Scope and results of this documentation verification, including skipped checks |

Existing detailed guides remain in place:
[dataset-splitting handbook](../README_DATASET_SPLITTING.md),
[benchmark handoff](../ML_BENCHMARK_HANDOFF.md),
[combined target dataset](../COMBINED_TARGET_DATASET.md), and
[example notebook](../examples/CoREMOF_dataset_splitting_quickstart.ipynb).
This workspace organizes and connects them; it does not replace their detailed
contracts or move existing datasets.

## Before manuscript submission

1. Confirm paper scope, author list/order, affiliations, funding, and journal.
2. Select and freeze the exact software commit and data snapshot to report.
   A development checkout is not evidence that these APIs are on stable PyPI.
3. Complete the model-training comparison, including model-input success
   counts and effective evaluation denominators. No predictive performance,
   speedup, or calibrated uncertainty result is claimed in this draft.
4. Select individual figures from the audited analysis bundles; match each
   caption to its own checker view, eligibility policy, feature space, and
   snapshot cutoff. Do not silently replace an old snapshot with live results.
5. Complete reference, software/dependency licence, and asset-level data-rights
   review. Existing tracked legacy assets mean the whole repository must not
   be described as an already sanitized data-free distribution.
6. Perform the final release/build/reproduction audit and obtain publication
   approval. All current benchmark assignments have `official_split=false`.

No scheduler action, target recalculation, release promotion, Git push, or
data transfer is part of preparing this documentation workspace.
