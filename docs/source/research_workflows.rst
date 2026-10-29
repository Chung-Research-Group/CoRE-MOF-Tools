Research workflows
==================

The examples connect target-free structure metadata, saved checker results,
related-structure grouping and explicit target attachment. They do not submit
calculations or train models.

Choose the analysis
-------------------

* Read metadata and recorded checker outcomes with :doc:`database_access`.
* Attach a requested optional target package with :doc:`target_supplements`.
* Construct a new target-complete benchmark with :doc:`target_first_benchmark`.
  Required finite targets determine eligibility before sampling, while grouping
  and checker-label purity use the complete release.
* Reproduce the paper's fixed assignments with :doc:`frozen_assignment_replay`.
  Do not rebuild its groups or replace the bound metadata with a newer release.
* Use :doc:`splitting` for general target-independent splits and later target
  attachment. These supported API recipes do not regenerate the paper dataset.

Executable examples
--------------------

* :download:`Grouped workflow recipes <../../examples/grouped_workflow_recipes.md>`
* :download:`Target-first benchmark builder <../../examples/build_target_first_benchmark.py>`
* :download:`Frozen-assignment replay <../../examples/replay_common_input_benchmark.py>`
* :download:`Optional target attachment <../../examples/attach_target_supplement.py>`
* :download:`Benchmark handoff guide <../../ML_BENCHMARK_HANDOFF.md>`

The grouped recipes are compiled and exercised on synthetic data by
``tests/test_manuscript_workflows.py``. The target-first and replay examples
have separate regression tests. These checks verify software behavior, not
scientific feature quality, model performance or public data permission.

Record source hashes, endpoint definitions, eligibility decisions and assignment
digests with each new analysis. Preserve native nulls and zeros. All generated
assignments are exploratory unless an independently audited official assignment
manifest establishes otherwise. Structure-resolved inputs and results follow
the distribution conditions in :doc:`database_access`.
