Separate database releases
============================================================

CoRE-MOF-Tools provides APIs; the separate CoRE-MOF-COD data repository provides
release catalogues and schemas. Full metadata and COD/SI CIF archives are
prepared for Zenodo, not bundled in the Python package. Modified CSD and
unmodified CSD CIFs have been built as separate local review packages for subsequent CCDC
handoff. Both CSD categories stay out of GitHub, GitHub Release assets and
Zenodo.

Publication destinations
------------------------

The public project branch is
https://github.com/Chung-Research-Group/CoRE-MOF-Tools/tree/2026-core-mof-cod.
After advisor approval, the planned workflow folder is ``2026-CoRE-MOF-COD``
in https://github.com/Chung-Research-Group/reproducible-workflows. A separate
public database repository is planned under Chung-Research-Group; its exact
URL is not yet assigned. The advisor-review repositories remain private.
The new database repository, workflow contribution and Zenodo deposit are
not claimed completed.

The current branch does not replace historically recorded fit implementations,
environments or checkpoints. Authorized reviewer access, exact computational
assets, recipient replay, source permissions and production MOFid admission
remain separate requirements. SI-CIF deposition requires asset-level rights
clearance. CSD CIF payloads and private mappings stay out of public deposits.

CoRE IDs
--------

Release structure names follow
``YYYY[elements][topology]dimension[variant]serial``. For example,
``2013[Cu][nan]3[ASR]5`` has framework dimensionality 3 and serial number 5.
The element segment includes metals or metalloids. Unknown publication years
use ``0000`` and no unambiguous named topology is represented by ``nan``.
Established IDs remain persistent by default; newer metadata alone does not
trigger renaming. A previously unavailable ``[nan]`` topology may be updated
only with explicit user approval and a validated, unambiguous named
CrystalNets SingleNodes result. Allocate the destination serial from the
maintained registry for its publication year, elements, topology, dimension
and variant, checking that the complete ID is collision-free. Update current
identifier-bearing release artifacts consistently. Retired aliases and
old-to-current crosswalks remain private. Frozen benchmark translations must
preserve structure membership and assignments; updated names alone do not
authorize a new split or scientific recalculation.

Use ``CoREMOF.parse_core_id`` and ``CoREMOF.format_core_id`` for exact parsing
and formatting. They do not allocate serials, infer source/access rights,
or construct related-structure groups. Read ``source_database`` and access
categories from release metadata. CoRE IDs are not chemical MOFid strings.
The compact topology token is not a substitute for the complete CrystalNets
evidence used by grouping. Quote bracket-containing shell arguments.

Release files and links
-----------------------

.. list-table:: Database distribution
   :header-rows: 1
   :widths: 50 20 30

   * - Resource
     - Structure records
     - Distribution route
   * - Full metadata/JSON, without CIF bytes
     - 42,574
     - Zenodo
   * - COD CIF archive
     - 19,598
     - Zenodo
   * - SI CIF archive
     - 5,727
     - Zenodo
   * - Modified CSD CIF package
     - 13,001
     - Local review, then CCDC
   * - Unmodified CSD CIF package
     - 4,248
     - Local review, then licensed CCDC access

The two CSD categories together cover 17,249 structures and retain the
release-recorded classification, not a new source-CIF comparison or permission
grant. Each local package has a separate membership manifest, checksum ledger
and access/licensing document. Category counts and hashes come from those manifests, not from
ASR/FSR/ION or CR/NCR labels. A CIF-only package does not reconstruct the
complete release-loader contract.

The built local archives are
``CoRE-MOF-COD_CSD_modified_cifs_20261001.zip`` and
``CoRE-MOF-COD_CSD_unmodified_cifs_20261001.zip``. Their membership records are
``manifests/ccdc_package_manifest.json`` and
``manifests/classification_manifest.csv`` inside each archive. The public
catalogue records aggregate counts and hashes, not structure-resolved CSD
membership.

The following links and DOIs are deliberately blank for the authors to fill
when released. Packaging does not claim an upload or waive the scientific
admission and permissions checks. Historical software/data DOI citations still
identify their own earlier releases.

.. list-table:: Current release links
   :header-rows: 1
   :widths: 70 15 15

   * - Resource
     - Link
     - DOI
   * - CoRE-MOF-COD metadata and COD/SI CIF deposit on Zenodo
     -
     -
   * - CoRE-MOF-Tools software archive on Zenodo
     -
     -
   * - Modified CSD CIF collection at CCDC
     -
     -
   * - Unmodified CSD CIF collection at CCDC
     -
     -

Metadata-only use
-----------------

Verify an authorized archive against its independently obtained checksum and
extract to a new directory. Complete metadata can then be loaded without CIFs:

.. code-block:: python

   from CoREMOF.dataset import CoREMOFDataset

   dataset = CoREMOFDataset.from_release(
       "/authorized/path/CoRE-MOF-COD", verify_cif_files=False
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

Unmodified CSD CIFs require the applicable valid CSD licence. A modified CSD
classification records an evidenced curation change, such as solvent removal
under the agreed classification, but does not itself grant redistribution
rights. Confirm the modified collection's own access terms rather than
inheriting a grant from an earlier collection. Both categories remain excluded
from GitHub/Zenodo CIF archives. Source names, ASR/FSR/ION variants, CR/NCR
labels and the scientific OA category do not prove redistribution rights.
Precomputed checker results can be read without redistributing checker engines.
The repository's ``README_DATABASE_ACCESS.md`` explains the rights separation
and the future publication sequence in more detail.
