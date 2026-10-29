---
name: coremof-dataset-use
description: Use CoRE-MOF-Tools to read released metadata, select structures from saved checker results, attach optional targets, or run the documented dataset examples.
---

# Use CoRE-MOF-COD with CoRE-MOF-Tools

This skill helps users operate the existing package with released or otherwise
authorized inputs. It does not describe how the database was constructed or
how to operate its internal curation and calculation campaigns.

## Choose the relevant guide

Read only the guide for the requested task. Paths below are relative to this
skill in the source checkout or extracted source distribution.

| Task | Guide and runnable example |
| --- | --- |
| Install the package | [Installation](../../../docs/source/installation.rst) |
| Load metadata and locate permitted CIF downloads | [Database access](../../../README_DATABASE_ACCESS.md), [read metadata](../../../examples/read_release_metadata.py) |
| Inspect saved checker results and select CR/NCR criteria | [Dataset handbook](../../../README_DATASET_SPLITTING.md), [read checker results](../../../examples/read_checker_results.py) |
| Attach explicitly requested target data | [Target supplements](../../../docs/source/target_supplements.rst), [attachment example](../../../examples/attach_target_supplement.py) |
| Create a new grouped dataset with the existing API | [Target-first example](../../../docs/source/target_first_benchmark.rst), [grouped recipes](../../../examples/grouped_workflow_recipes.md) |
| Read an existing benchmark's assignments without resampling | [Assignment replay](../../../docs/source/frozen_assignment_replay.rst) |

Run example scripts from the repository root and inspect their `--help` for
required inputs. Use documentation matching the installed source revision.

## Interpret inputs and outputs correctly

- Core metadata is target-free. Attach the optional target supplement only
  when requested, using its verified manifest hash and exact CoRE IDs.
- CR means computation-ready and NCR means non-computation-ready under the
  selected checker policy. Unavailable checker results are not FAIL votes.
- Read source, variant and access conditions from metadata, not guesses from
  filenames. Use the package's CoRE-ID parser for identifiers.
- Preserve zero and missing values. Do not fill missing results or infer new
  scientific values from names or neighbouring records.
- Use the documented grouping options. Replaying a benchmark preserves its
  recorded groups, assignments and targets. A new split is a separate dataset.
- Reading metadata and checker results does not run checker software. No
  third-party checker implementation or licensed CIF payload is supplied here.

Summarize the requested outputs and relevant missing inputs. Do not initiate
feature calculations, curation, training or publication as a side effect of
reading a dataset or running a read-only example.
