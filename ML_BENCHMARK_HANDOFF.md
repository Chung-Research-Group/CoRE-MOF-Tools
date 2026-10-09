# Using the CoRE-MOF-Tools benchmark APIs

Use the documentation shipped with your installed version of CoRE-MOF-Tools.
This page links to supported user workflows. It does not distribute the
project's internal curation, calculation, or workstation instructions.

## Choose the workflow for your analysis

- [Dataset-splitting handbook](README_DATASET_SPLITTING.md): checker-based
  selection, related-structure grouping, partitioning, benchmark construction,
  and target attachment.
- [Target-first benchmarking](docs/source/target_first_benchmark.rst): construct
  eligible datasets using the documented target-first API.
- [Frozen assignment replay](docs/source/frozen_assignment_replay.rst): reuse
  existing assignments rather than creating a different experiment.
- [Installation](docs/source/installation.rst): installation options and
  supported environments for this version.

Choose one documented workflow explicitly. Target attachment to an existing
split and target-first cohort construction are different operations. Do not
substitute one for the other when reproducing an experiment.

## Use your authorized inputs

Supply the release, checker criteria, target table, grouping criteria, and
assignment configuration required by your analysis. Example settings are not
an official database split. Record the settings and receipts produced by the
API, and check the reported exclusions and achieved partition sizes.

Availability of source code does not grant redistribution rights to source
CIFs or other third-party data. Obtain required inputs through their documented
access routes and respect the terms that apply to each input.

For agents helping users run these APIs, the repository includes only the
[dataset-use guide](.agents/skills/coremof-dataset-use/SKILL.md). Internal
development and curation skills are not part of the public distribution.
