Replaying the manuscript's frozen assignments
=============================================

The common-input experiment ``coremof_cod_common_input_equal_size_20260912_v1``
has twelve shared assignments: seeds 912–915 and NCR-pool fractions
``q=0, 0.5, 1``. Every assignment has 3,737 training, 466 validation and
468 test structures. The same pure-CR test is used throughout. Here ``q``
is the fraction of the eligible NCR pool, not the training NCR percentage.

The executable example ``examples/replay_common_input_benchmark.py``
reconstructs these assignments with the currently installed package's exact
group sampler. It requires the private September 13 workflow handoff and its
independently verified archive checksum::

    python examples/replay_common_input_benchmark.py \
      --handoff-archive /private/coremof_cod_workflow_handoff_20260913_v1.tar.gz \
      --archive-sha256 VERIFIED_ARCHIVE_SHA256 \
      --output /private/new-assignment-replay

Run from an installed environment, or set ``PYTHONPATH`` to the checkout when
using ``python -B -S``. The output directory must not exist. Structure-resolved
inputs and results remain private and are not included in the code repository.

What is reproduced
-------------------

The example reuses the archived full-release related-structure groups,
representative-diversity strata, original test IDs and recorded common-input
availability decisions. It excludes complete affected eligible groups,
verifies full-release group label purity, then reruns nested NCR selection,
equal-sized CR removal and group-preserving train/validation assignment.
Each structure remains an individual row. No target or prediction table is
read, no archived Python code is executed, and no files are extracted.

Both the canonical assignment digest and every exported membership row must
match the frozen experiment. The required assignment SHA-256 is
``9e72992970518d039f9631b1f45b516f4ff3603f7945bdcb82dbad283529dcbd``.
Output includes the reproduced CSV, a receipt binding the inputs and current
code, partition counts, a group-crossing audit and a checksum ledger.
Any mismatch fails without publishing an output directory.

What is not reproduced
-----------------------

This replay does not recalculate descriptors, rebuild the original groups or
diversity index, certify graph/grid preprocessing, regenerate targets or
train models. Those stages have separate scientific/runtime contracts.
The frozen checker view and model-input eligibility are not replaced with a
later release's metadata or a newly completed target table. A passing replay
is not authorization to publish licence-gated data or call the assignment an
official database split.

For a genuinely new target-first experiment, use the target-first benchmark
example instead. It is a separate workflow, not a replacement for this replay.
