# Portable ML benchmark workflow selection

Read the parent `SKILL.md` for the shared scientific and transfer boundaries.
Choose the guide matching the requested experiment, not the newest directory
name or a previous chat's counts.

- For a new target-first workflow, read
  [target_first_benchmark.rst](../../../../docs/source/target_first_benchmark.rst)
  and use `examples/build_target_first_benchmark.py`. Required finite targets
  determine eligibility before cohort selection. Grouping and checker-label
  purity still use the complete release, including target-missing structures.
  The example writes the source/endpoint/eligibility/assignment bindings.
- For the historical target-independent experiment, read repository-root
  `ML_BENCHMARK_HANDOFF.md`. This is the dated command, data-layout,
  deferred-target and restricted-transfer guide, not the current target-first
  protocol. Keep its frozen assignments and recorded settings unchanged.
- For an already transferred experiment, its exact manifest and receiver
  prompt define which workflow and metadata view apply. Neither guide grants
  permission to overwrite earlier results or start training on a new dataset.

For combined target construction and the audited as-of-cutoff coverage, also
read repository-root `COMBINED_TARGET_DATASET.md` as the historical baseline.
`examples/extend_collected_targets.py` can extend a compatible accepted target
snapshot with saved collector evidence, preserving values and scientific nulls
and writing a new private candidate. This is not a new split or public release.
Treat completion-only counts as source contributions, never total availability.
Do not substitute an original-host curation runbook or infer production state
on the receiving machine.
