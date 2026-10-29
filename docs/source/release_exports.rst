Core metadata and explicit target exports
============================================================

Current core releases are target-free. For the separately versioned optional
target supplement and explicit attachment, see :doc:`target_supplements`.
The enrichment examples below remain available for deliberate private analysis
and historical experiment reproduction. Their combined output is not the
current core metadata layout and must not replace a target-free release.

The source distribution includes two executable examples for collecting
existing accepted targets and adding them to a new private metadata copy.
Neither example runs a simulation, trains a model, modifies a source CIF or
changes a frozen train/validation/test assignment.

These examples do not grant publication or redistribution permission. The
final MOFid evidence gate and any asset-specific restrictions remain in force.
A successful metadata export is not an approved Zenodo release.

1. Extend the accepted target table
-----------------------------------

``examples/extend_collected_targets.py`` reads an accepted target snapshot and
an explicit inventory of completed collector reports. It verifies the saved
task, result and receipt hashes, then fills eligible missing cells. Existing
accepted values and accepted scientific nulls are not replaced. Running jobs,
unvalidated results and conflicting observations cannot supply new values.
Use ``--help`` for the required baseline and inventory arguments.

The baseline may be the original snapshot with a ``targets/`` subdirectory
or a previous flat output of this example. Select one unambiguous layout.
Repeated collector observations are checked against the accepted record and
are not counted again. Where a verified worker receipt records the RASPA
version, electrostatics setting or CIF-preprocessing action, these are carried
forward. Missing method information is left unavailable.

The output retains all structure IDs. A missing value is not zero. Failed,
excluded and mathematically undefined targets remain distinguishable, and
no cohort is refilled or rebalanced.

2. Add accepted targets to release metadata
--------------------------------------------------

Install the release-validation dependency from a source checkout::

    python -m pip install ".[release]"

Then run the exporter using hashes from the independently verified input
receipts, not hashes copied from an untrusted download::

    python examples/extend_release_metadata.py \
      --source /path/to/validated/source-release \
      --targets /path/to/accepted-target-snapshot \
      --output /path/to/new/private-release \
      --source-ledger-sha256 SOURCE_LEDGER_SHA256 \
      --target-ledger-sha256 TARGET_LEDGER_SHA256 \
      --revision-id YYYYMMDD-target-metadata-v1

The destination must not exist. All source payloads are checked against the
source checksum ledger. Each target is matched by exact structure ID and
full CIF hash. The source release may be a base release and the accepted
target table its complete superset. Every retained structure must still have
all three declared endpoint records, including explicit unavailable records.

The new directory contains:

* unchanged CIFs, checker findings, descriptor values and grouping tables;
* target values and statuses in ``metadata/metadata.csv``;
* typed targets in ``metadata/metadata.jsonl`` and
  ``metadata/targets.jsonl``;
* one schema-validated ``metadata/structures/<structure_id>.json`` per CIF;
* the extended ``coremof-structure-record/1.1`` schema;
* updated manifests, an input-checksum record, and ``SHA256SUMS``.

Historical schema 1.1 adds a ``targets`` object to the unchanged scientific fields of
schema 1.0. This is an explicit schema revision, not an in-place modification
of an older release. Existing files and frozen experiments remain intact.
On failure, the destination is not published. A private temporary staging
directory is retained for diagnosis.

Source-separated review bundles
--------------------------------

An archive containing ``source_bundle.json`` instead of ``dataset_info.json``
is a **review slice**, not a standalone release-loader input. Its parent group
sizes still describe the complete release. The loader rejects this layout
without changing the files. Do not manufacture a release contract or reduce
group sizes to make validation pass.

For source-specific analysis with the existing API, load a complete release
that you are authorized to access. Filter through the public API so the
underlying full-release relationships remain available:

.. code-block:: python

   from CoREMOF.dataset import CoREMOFDataset

   dataset = CoREMOFDataset.from_release("/path/to/complete-local-release")
   classified = dataset.classify("5checker")
   cod = classified.filter(sources=("COD",))
   print(cod.label_counts())

   split = classified.data_split(
       sources=("COD",),
       labels=None,  # retain every checker label for this example
       group_criteria="priority_main",
       diversity="none",  # grouping-only example, not the paper benchmark
       random_state=42,
   )

``diversity="none"`` is explicit here to demonstrate grouping without the
optional numerical backend. For representative feature balancing, use the
pinned benchmark dependencies and ``diversity="representative"``. This
example does not recreate or replace a frozen benchmark.

The full-release leakage groups are constructed before source filtering.
Two COD records may be linked through an omitted SI or CSD record, so
rebuilding groups from COD rows alone is not equivalent. The opaque
``full_release_split_guards.csv`` in a review slice preserves the two declared
grouping profiles, but is not currently a package loader contract for arbitrary
criteria or benchmark construction. Combining separately projected criteria
can also miss paths through omitted records. In a CR/NCR benchmark, purity
must be assessed using the labels of **all** full-release group members, not
only the selected source.

Standalone source-projection contracts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``CoREMOFDataset.from_projection`` accepts a review slice together with a
separate, hash-verified contract generated from the complete release. The
slice alone remains insufficient. This is an additive loader, not a relaxed
``from_release`` validator. It does not rewrite any source files or reduce
their full-release parent group sizes.

The producer, who must have authorized access to the complete release, runs:

.. code-block:: python

   from CoREMOF import export_source_projection
   from CoREMOF.dataset import CoREMOFDataset

   complete = CoREMOFDataset.from_release("/path/to/complete-local-release")
   receipt = export_source_projection(
       complete,
       "/path/to/unchanged-COD-review-slice",
       "/path/to/new/COD_projection.json",
       sources=("COD",),
       group_profiles=(("priority_main",), ("RT", "M2T")),
       diversity="representative",
   )
   print(receipt["sha256"])

The destination and its checksum companion must not already exist. The
producer checks exact selected metadata, grouping, feature and CIF-manifest
rows against the complete release. Numeric strata require the pinned
``benchmark`` environment at export time. Use ``diversity="none"`` explicitly
if only grouping is required. An executable wrapper is supplied as
``examples/export_source_projection.py``.

The receiver needs only the selected-source files and the new contract:

.. code-block:: python

   selected = CoREMOFDataset.from_projection(
       "/path/to/COD-review-slice",
       "/path/to/COD_projection.json",
       expected_sha256="EXPECTED_SHA256_FROM_TRUSTED_HANDOFF",
       verify_cif_files=True,
   )
   split = selected.classify("5checker").data_split(
       labels=None,
       group_criteria=("RT", "M2T"),
       diversity="representative",
       random_state=42,
   )

The contract retains the exact **ordered, combined** grouping profiles,
including paths through omitted records, and all-member five-checker counts
for each represented group. It contains no omitted structure IDs or omitted
CIFs. Individual selected records retain their original checker labels.
Strict CR/NCR cohort purity is evaluated from the saved complete-group
counts. Pool sizes and prediction diagnostics describe the selected sources,
not the complete release. Default mixed-group and infeasibility errors still
apply. Any explicit label-pure cohort selection remains recorded in receipts.

Saved representative strata are reused, never fitted again on the source
subset. Reuse needs no numerical backend on the receiver and is explicitly
identified as a cached complete-release index in the receipt. This is not a
fallback: requesting representative diversity when it was not exported fails.
Missing scientific features remain in their original availability tier.

Use ``available_group_criteria(selected)`` to inspect the available ordered
profiles. A criterion appearing within a combined profile is not separately
usable unless that separate profile was exported. Unsupported combinations
fail rather than rebuilding incomplete relationships. Legacy parent resolver
and legacy split calls cannot reconstruct omitted-source links and therefore
reject projection datasets explicitly. Existing complete-release calls are
unchanged. Target merging and frozen target attachment use the selected IDs
and preserve their assignments.

A checksum verifies integrity, not a digital signature. Obtain the expected
value from a trusted handoff rather than calculating it from an untrusted
download. All projection outputs remain exploratory. Neither the contract
nor its successful validation grants redistribution permission or resolves
the release's MOFid gate. Frozen earlier benchmarks are not replaced.

Target meanings and missing data
----------------------------------------

The accepted endpoint contract is:

.. list-table::
   :header-rows: 1
   :widths: 50 25 25

   * - Target
     - Conditions
     - Unit
   * - CH4 absolute uptake
     - 298 K, 65 bar
     - mol/kg-framework
   * - H2 absolute uptake
     - 77 K, 100 bar
     - mol/kg-framework
   * - CO2/N2 Henry selectivity
     - 298 K, infinite dilution
     - dimensionless

Henry selectivity has no finite-pressure condition in the target record.
The exporter does not infer volumetric uptake or individual gas Henry
constants. The target table supplies accepted simulation results, not
experimental measurements or standardized model predictions.

Each endpoint retains availability, execution status, its value, a reported
uncertainty when recorded, conditions and a compact missing-result reason.
The ``calculation_method`` object preserves the recorded RASPA version,
electrostatics setting and CIF-preprocessing action. A null method field means
that this input table did not record the setting. In particular, ``Ewald`` is
an electrostatics method, not the method used to generate atomic charges.
No missing field is guessed, imputed or filled from a different structure.

Validation scope
----------------

The exporter verifies source hashes, exact target/CIF correspondence,
CSV/JSONL membership, schema validity and preservation of the original
scientific record. Its integrity manifest records the source-release ledger,
accepted-target ledger and assignment-table hashes, plus the builder hash.
It does not independently rerun the source release's scientific calculations
or adjudicate unresolved MOFid cases. Those are separate release requirements.

Before public distribution, run the full release and cross-release audits,
verify extracted archives, and resolve the release's scientific and sharing
gates. Do not publish licence-gated structure-resolved material merely because
it is present in a validated local candidate.
