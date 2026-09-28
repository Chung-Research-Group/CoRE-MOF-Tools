Target-complete CR/NCR benchmark datasets
================================================

When a benchmark requires complete labels, merge target data and select eligible
IDs before building cohorts. The grouping graph must still use the complete
release, including missing-target rows that connect retained structures.

.. code-block:: python

   import math
   from CoREMOF.dataset import CoREMOFDataset
   from CoREMOF.targets import merge_targets_from_config

   release = CoREMOFDataset.from_release("/path/to/validated_release")
   merged = merge_targets_from_config(release, "targets.json")
   required = tuple(merged.target_columns)  # explicitly choose the endpoints
   eligible = tuple(
       sid for sid in release.structure_ids
       if all(
           isinstance(merged.target_values(sid)[name], (float, int))
           and not isinstance(merged.target_values(sid)[name], bool)
           and math.isfinite(merged.target_values(sid)[name])
           for name in required
       )
   )
   suite = release.classify("5checker").build_cr_ncr_benchmark(
       eligible_structure_ids=eligible,
       group_criteria=("RT", "M2T"),
       cohort_eligibility="complete_release_label_pure_effective_blocks",
       ncr_pool_fractions=(0, 0.5, 1),
       seeds=(912, 913, 914, 915),
       partition_strategy="transition_balanced",
       include_full_cr_diagnostic=False,
   )
   suite.write("splits")
   suite.attach_targets(merged, missing="error").write("model_inputs")

``RT`` requires exact equality of all 264 finite depth-5 RAC descriptors and
the complete successful CrystalNets fingerprint. ``M2T`` requires the complete
canonical MOFid-v2 string, including chemistry, topology and catenation, and
that independent CrystalNets fingerprint. These criteria must be declared in
the input release. Their union is combined with the full-release ``main_union``
leakage guard. Missing evidence never matches another missing value.

The optional ID filter leaves the original grouping IDs intact. It does not
refit grouping on a smaller release. Under the explicit label-purity policy,
an eligible structure in a group with opposite, ambiguous or unchecked labels
is excluded even when those other members have no target. All eligible members
of each selected group are indivisible during sampling and splitting.

The NCR-pool fraction refers to the eligible NCR pool after filtering, not the
fraction of all database NCR structures or the final NCR composition. Save the
merge receipt, required-target list and eligible-ID digest beside the suite
receipt. A changed target-availability set requires a new version of the
benchmark dataset. Never overwrite a prior assignment or reuse its predictions
as results for the rebuilt dataset.

Omitting ``eligible_structure_ids`` preserves the original target-independent
cohort behavior. Target values are never used for grouping or balancing.
Missing numerical descriptors are not imputed. Model graph/grid/feature
eligibility must be checked separately before training.

The opt-in ``transition_balanced`` partition strategy pairs equal achievable
validation counts for each CR-removal/NCR-addition transition. This keeps
partition sizes constant across the NCR ladder without moving persistent
structures. It requires q=0 as the first level. Diversity balance operates
within transitions. Whole-group rounding deviations are reported. Omitting
this option retains the original bounded partition optimizer.

Executable workflow and provenance
----------------------------------

The checkout supplies :download:`build_target_first_benchmark.py
<../../examples/build_target_first_benchmark.py>`, an executable version of
the current target-first workflow. It loads an extracted local release,
merges the declared targets, checks finite eligibility, builds the paired
cohorts, attaches their individual targets, and writes a workflow receipt:

.. code-block:: bash

   python examples/build_target_first_benchmark.py /path/to/validated_release \
     --target-config targets.json \
     --require-target ch4_loading \
     --require-target h2_loading \
     --require-target widom_ratio \
     --output-directory target_first_outputs

Replace the three example names with the exact endpoint names in your target
configuration. Declare each required endpoint's units and conditions. For CSV,
declare ``float`` or ``int`` value types; native JSON numbers are also accepted.
Finite zero remains eligible. Null, boolean and nonnumeric values are excluded
from eligibility, and the target parser rejects nonfinite numeric observations
without rewriting them. Additional configured endpoints may remain optional.

The default example uses RT and M2T, the full-release leakage guard and label
purity, representative diversity, NCR-pool fractions 0/0.5/1, seeds 912--915,
requested 80/10/10 partitions, a common pure-CR test and the explicit
``transition_balanced`` strategy. Install the pinned benchmark extra first.
``--diversity none`` is available only as an explicit non-representative test
or sensitivity choice. These defaults belong to this example and do not change
the package API's existing defaults.

The new output directory contains the materialized target merge, a full-release
target-eligibility table, split and target-attached suites, ``SHA256SUMS`` and
``workflow_receipt.json``. The workflow receipt binds required endpoint
definitions, the merge receipt, eligible IDs and their digest, full-release
grouping and label purity, requested settings, assignment digests and every
written artifact. It verifies that attachment preserves assignments and that
all required values remain finite. It is written last; a directory without it
is incomplete. Existing output directories are rejected.

This example creates new exploratory assignments. It does not reproduce an
existing paper dataset or certify model graph, grid or descriptor inputs.
