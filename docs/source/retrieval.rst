Version-selected release retrieval
==================================

``CoREMOF.retrieval.fetch_release`` retrieves one exact release version from a
caller-supplied, checksum-pinned JSON catalog. It downloads each declared file,
checks its byte count and SHA-256, validates the staged release with
``CoREMOFDataset.from_release``, and publishes a new local directory only after
all checks pass. This API does not bundle a hosted catalog, assert that a
particular version is remotely available, or grant access to licensed data.
``CoREMOFDataset.from_release`` remains the local release loader. Historical
structure downloads are separate and have the output protections described below.

Retrieve and open a release
--------------------------------

Obtain a catalog path or HTTP(S) URL and its full lowercase SHA-256 from a
trusted provider. A hash downloaded alongside an untrusted catalog is not an
independent authenticity guarantee. Keep the pinned catalog for provenance.
Use the exact version key present in that catalog; there is no implicit
``latest`` resolution or automatic fallback to a different version.

.. code-block:: python

   from CoREMOF.dataset import CoREMOFDataset
   from CoREMOF.retrieval import fetch_release

   release_root = fetch_release(
       catalog_path_or_url,
       requested_version,
       "retrieved-release",
       catalog_sha256=trusted_catalog_sha256,
       verify_cif_files=True,
       max_total_bytes=20_000_000_000,
   )
   dataset = CoREMOFDataset.from_release(release_root, verify_cif_files=True)

The checkout includes the executable equivalent:

.. code-block:: console

   python examples/fetch_release.py CATALOG VERSION NEW_DIRECTORY \
       --catalog-sha256 FULL_SHA256 --verify-cifs --max-total-bytes 20000000000

The capitalized arguments above are caller-provided values, not published
release locations. HTTP(S), local JSON paths, and local ``file:`` URLs are
supported. Resource URLs may be relative to the catalog's final URL after any
redirect. URLs containing username/password credentials are rejected; signed
query URLs are accepted but never copied to the retrieval receipt. The API
does not manage provider login, retry failed requests, or extract archives.

Every file listed in the selected catalog entry is always size- and
checksum-verified. ``verify_cif_files=True`` additionally requires every CIF
named by the release's own CIF manifest and checks its release-declared hash.
With the default ``False``, a valid metadata-only release is permitted; this
does **not** assert the presence or verification of every CIF. The flag and
actual verification status are recorded in the receipt.

Catalog contract
----------------

A catalog is UTF-8 JSON with exactly ``schema_version`` and ``releases`` keys.
Its schema version is ``coremof-release-catalog/1.0``. ``releases`` maps exact
dataset-version strings to objects containing exactly a nonempty ``files``
list. Each file contains exactly:

* ``path``: canonical relative POSIX destination path within the release;
* ``url``: absolute HTTP(S)/local-file URL or URL relative to the catalog;
* ``size_bytes``: nonnegative integer, not a Boolean;
* ``sha256``: the complete lowercase 64-character SHA-256 of the file bytes.

The selected entry must include a complete loadable release layout, including
``dataset_info.json``, ``metadata/metadata.csv``,
``parent_groups/parent_groups.csv``, and
``parent_groups/parent_group_methods.json``. Include the release's manifests,
method receipts, and any other files its declared contracts require. Include
the CIF files when providing a complete structure release. The downloaded
``dataset_info.json`` version must equal the selected catalog version exactly.

For a provider preparing a catalog from an existing, authorized local release,
the following creates a working local-file catalog without assuming a public
host. ``release_root`` must name the intended release directory, not a parent
directory containing unrelated files. Do not redistribute licensed files
without permission.

.. code-block:: python

   import hashlib
   import json
   from pathlib import Path
   from CoREMOF.dataset import CoREMOFDataset
   from CoREMOF.retrieval import CATALOG_SCHEMA, RECEIPT_PATH

   release_root = Path("authorized-local-release").resolve()
   dataset = CoREMOFDataset.from_release(release_root, verify_cif_files=True)
   files = []
   for path in sorted(release_root.rglob("*")):
       if path.is_symlink():
           raise ValueError("Catalog provider must resolve symlinks explicitly")
       if not path.is_file():
           continue
       logical = path.relative_to(release_root).as_posix()
       if logical == RECEIPT_PATH:
           continue  # Retrieval writes its own reserved receipt.
       digest = hashlib.sha256()
       size = 0
       with path.open("rb") as handle:
           for block in iter(lambda: handle.read(1024 * 1024), b""):
               size += len(block)
               digest.update(block)
       files.append({"path": logical, "url": path.as_uri(),
                     "size_bytes": size, "sha256": digest.hexdigest()})
   catalog = {"schema_version": CATALOG_SCHEMA,
              "releases": {dataset.dataset_version: {"files": files}}}
   catalog_bytes = (json.dumps(catalog, indent=2, sort_keys=True) + "\n").encode("utf-8")
   Path("local-release-catalog.json").write_bytes(catalog_bytes)
   print(hashlib.sha256(catalog_bytes).hexdigest())

That local catalog contains local source paths; keep it private when those
paths are sensitive. A remote provider can instead supply relative URLs to
the exact same immutable file bytes. Trust and retain the complete catalog
hash, not a mutable catalog URL alone. JSON duplicate keys, nonfinite values,
unsafe paths, case/Unicode-normalized path collisions, file/directory conflicts,
and the reserved ``manifests/retrieval_receipt.json`` path are rejected. Catalogs
are bounded to 64 MiB, and ``max_total_bytes`` can bound the selected payload.

Failure handling and provenance
-------------------------------

An existing destination is never intentionally replaced. Staging happens in a
temporary sibling directory and is removed on failure. Publication uses the
package's atomic directory writer; on filesystems without a no-replace rename,
the fallback serializes cooperating package writers. Avoid a destination
concurrently modified by unrelated processes.

``manifests/retrieval_receipt.json`` records the selected version, pinned
catalog hash and byte count, every downloaded logical path/hash/byte count,
package and retrieval-implementation versions or hashes, release-loader input
hashes, structure count, and actual CIF-verification status. It does not embed
the source catalog, source URLs, signed credentials, or private local source
paths. The receipt binds what was retrieved; it does not independently certify
the provider's scientific claims or data-access rights.

Historical SI and CSD downloads
-------------------------------

``CoREMOF.structure.download_from_SI`` and ``download_from_CSD`` retain their
historical destinations and successful output content. They do not select a
CoRE-MOF-COD release version. SI/metadata cache URLs refer to the repository's
mutable branch, while CSD retrieval uses the user's installed, licensed CSD.
Use the checksum-pinned catalog workflow above for version-selected releases.
Access to an asset does not itself grant permission to redistribute it.

Existing output files now raise ``FileExistsError`` by default. To intentionally
replace only the selected files, pass ``overwrite=True``. Unrelated files in
the destination are retained. This is an explicit safety change from older
versions that silently overwrote files. Positional arguments, default output
folder and successful return values remain compatible. ``pathlib.Path`` output
folders are supported.

.. code-block:: python

   from CoREMOF.structure import download_from_SI, download_from_CSD

   # Use separate, new destinations. CSD requires the user's licensed API.
   download_from_SI("historical-si")
   download_from_CSD("ABCDEF", "licensed-csd")  # Replace with an actual refcode.

The command-line example requires an explicit source and destination:

.. code-block:: console

   python examples/retrieve_legacy_structures.py si historical-si
   python examples/retrieve_legacy_structures.py csd licensed-csd --refcode ABCDEF

Both SI archives are validated and completely staged before output publication.
Absolute/traversing paths, duplicate names, symlinks, special files, CRC errors
and file/directory conflicts are rejected. Existing output symlinks are rejected.
Downloaded cache files are validated before publication and a concurrently
created cache is not replaced. Existing cache files are not silently repaired.
CSD refcodes must be safe filename tokens. The reader is closed on success and
failure, and the exported CIF text is preserved without chemical changes.

Publication is per-file, with rollback on errors, not a simultaneous multi-file
transaction visible to all readers. Do not concurrently modify output directories
with unrelated processes. If rollback cannot safely restore an existing file,
its private staging directory is retained and its path is attached to the
exception as ``coremof_preserved_staging_directory``. Inspect that directory
before cleanup. Offline regression tests use synthetic ZIP/CIF fixtures and a
fake CSD reader, not unlicensed access to CSD or a scientific retrieval audit.

.. automodule:: CoREMOF.retrieval
   :members: fetch_release, RetrievalError
