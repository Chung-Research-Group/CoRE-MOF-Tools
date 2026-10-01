Separate database releases
============================================================

CoRE-MOF-Tools provides APIs; the separate CoREMOF-COD data repository provides
release catalogues and schemas. Bulk metadata and permitted CIF archives belong
in versioned data deposits, not the Python package. Confirmed URLs and data DOIs
will be linked at publication, rather than guessed here.

Metadata-only use
-----------------

Verify an authorized archive against its independently obtained checksum and
extract to a new directory. Complete metadata can then be loaded without CIFs:

.. code-block:: python

   from CoREMOF.dataset import CoREMOFDataset

   dataset = CoREMOFDataset.from_release(
       "/authorized/path/CoREMOF-COD", verify_cif_files=False
   )
   view = dataset.classify("5checker")
   print(dict(view.label_counts()))

This retains grouping and classification contracts but does not verify absent
CIF bytes or authorize publication. The example
``examples/read_release_metadata.py`` checks an explicit version and optionally
a trusted metadata checksum ledger. Full-release verification should run on a
compute node when inspecting tens of thousands of records.

Source-specific inputs
----------------------

Use an authenticated ``CoREMOFDataset.from_projection`` contract, exported from
the authorized complete release, for source-only input. Neither a JSONL slice
nor a CIF-only overlay reconstructs links through omitted structures. See
:doc:`release_exports` for the existing export API.

Archive and file-level catalogues
------------------------------------------------------------

The database publication catalogue describes archive assets, licences, approval
status, DOIs and download mirrors. Its archive downloader is separate from the
existing package's ``coremof-release-catalog/1.0`` file-level retrieval contract.
Use the appropriate format and exact version; do not pass one to the other.
The file-level API remains unchanged. See :doc:`retrieval`.

Reproducibility and permissions
------------------------------------------------------------

New target-complete benchmarks attach targets by exact ID before missing-target
eligibility filtering and cohort splitting. They do not regenerate old frozen
experiments. See :doc:`target_first_benchmark` and
:doc:`frozen_assignment_replay`. Current candidate admission/publication gates
remain unchanged, and new splits are exploratory.

Unmodified CSD CIFs remain licence-gated. Source names, ASR/FSR/ION variants and
the scientific OA category do not themselves prove redistribution rights.
Precomputed checker results can be read without redistributing checker engines.
The repository's ``README_DATABASE_ACCESS.md`` explains the rights separation
and the future publication sequence in more detail.
