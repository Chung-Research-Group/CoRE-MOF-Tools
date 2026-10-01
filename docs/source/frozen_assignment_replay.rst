Replaying the manuscript's frozen assignments
=============================================

The common-input experiment ``coremof_cod_common_input_equal_size_20260912_v1``
has twelve shared assignments: seeds 912–915 and NCR-pool fractions
``q=0, 0.5, 1``. Every assignment has 3,737 training, 466 validation and
468 test structures. The same pure-CR test is used throughout. Here ``q``
is the fraction of the eligible NCR pool, not the training NCR percentage.

The executable example ``examples/replay_common_input_benchmark.py``
validates and copies these frozen assignments without rerunning a group
sampler. It requires the approved workflow handoff and its
independently verified archive checksum::

    python examples/replay_common_input_benchmark.py \
      --handoff-archive /private/approved_workflow_handoff.tar.gz \
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
verifies full-release group label purity, nested NCR additions, equal-sized
CR removals, constant partition sizes, and stable within-seed partitions.
It preserves archived row order and exact assignments. This is essential when
identifiers change, because sorting new names and resampling would create a
different experiment.
Each structure remains an individual row. No target or prediction table is
read, no archived Python code is executed, and no files are extracted.

The canonical typed assignment digest must match
``suite_assignment_sha256`` in the receipt bound by the verified archive.
An identifier-translated export has its own digest, not the original file's
digest. Translation does not authorize changing membership or scientific data.
Output includes the reproduced CSV, a receipt binding the inputs and current
code, partition counts, a group-crossing audit and a checksum ledger.
Any mismatch fails without publishing an output directory.

What is not reproduced
-----------------------

This replay does not resample, recalculate descriptors, rebuild the original groups or
diversity index, certify graph/grid preprocessing, regenerate targets or
train models. Those stages have separate scientific/runtime contracts.
The frozen checker view and model-input eligibility are not replaced with a
later release's metadata or a newly completed target table. A passing replay
is not authorization to publish licence-gated data or call the assignment an
official database split.

For a genuinely new target-first experiment, use the target-first benchmark
example instead. It is a separate workflow, not a replacement for this replay.
