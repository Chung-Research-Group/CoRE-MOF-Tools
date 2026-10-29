Target-free metadata and optional targets
=========================================

Core metadata contains structure identity, structural features, checker results
and related-structure information. Simulation targets are a separate, optional
versioned supplement. Loading a core release never discovers, downloads or
attaches that supplement automatically.

Explicit attachment
-------------------

Obtain the checksum of the supplement's ``manifest.json`` from the trusted
delivery record, then attach it explicitly:

.. code-block:: python

   from CoREMOF.dataset import CoREMOFDataset

   core = CoREMOFDataset.from_release("/path/to/CoRE-MOF-COD")
   with_targets = core.attach_target_supplement(
       "/path/to/CoRE-MOF-COD_targets",
       expected_sha256="SUPPLEMENT_MANIFEST_SHA256",
   )

The loader checks the package kind/version, exact core input hashes, complete
structure-ID universe, file hashes, CSV/JSONL agreement, units and conditions.
It joins only the three response values. Full execution statuses, diagnostics,
reported uncertainties and calculation methods remain in ``targets.jsonl``.
They are not additional prediction targets. Zeros remain zero and missing
values remain null. Attachment does not filter, rebalance, group or split data.

The responses are CH4 uptake at 298 K/65 bar and H2 uptake at 77 K/100 bar,
both in mol/kg-framework, and dimensionless CO2/N2 Henry selectivity at 298 K.
No volumetric uptake or individual gas Henry constants are inferred.

User-supplied targets continue to use the existing ``merge_targets`` and
``TargetSource`` interfaces. A supplement is neither a CIF download nor a
redistribution permission. Source-specific access and release-admission
conditions apply separately.

Historical experiment reproduction
----------------------------------

Earlier combined metadata files remain unchanged. Ordinary loading hides
their target values, statuses, summaries and target-file declarations. A
workflow explicitly reproducing that historical combined layout can opt in:

.. code-block:: python

   historical = CoREMOFDataset.from_release(
       "/path/to/pinned/combined-release", include_legacy_targets=True
   )

This compatibility option does not make an old archive the current release,
and does not allow embedded targets in a new ``target_free_core/1.0`` release.
Keep frozen benchmark inputs, target sources and assignments pinned. New
target-first eligibility workflows attach their requested targets explicitly
before cohort selection. Target magnitudes never define structural grouping.
