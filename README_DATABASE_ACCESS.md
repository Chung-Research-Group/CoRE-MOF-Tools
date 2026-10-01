# Separate database releases and Python tools

CoRE-MOF-Tools contains the API and examples. The separate **CoRE-MOF-COD**
database repository contains release documentation, schemas and versioned data
catalogues. Full metadata and the COD/SI CIF archives are prepared for a
Zenodo deposit, not this Python source tree. All CSD CIFs stay out
of GitHub, GitHub Release assets and Zenodo. Modified and unmodified CSD CIFs
have been built as separate local packages for review and subsequent CCDC handoff.

## CoRE IDs

Current records use `YYYY[elements][topology]dimension[variant]serial`, for
example `2013[Cu][nan]3[ASR]5`. The fields are publication year, metals or
metalloids, topology, bonded-framework dimensionality, curation variant and
distinguishing serial number. `0000` denotes an unknown publication year and
`nan` means that no unambiguous named topology was assigned. The dimension is
not pore-channel dimensionality. Established published IDs are persistent,
so their descriptive tokens need not change when newer metadata is available.

```python
from CoREMOF import parse_core_id

identity = parse_core_id("2013[Cu][nan]3[ASR]5")
print(identity.dimension, identity.variant)  # 3 ASR
```

The source database (`COD`, `CSD`, or `SI`) and access category are explicit
metadata fields, not encoded in a CoRE ID. Join records by the exact
`structure_id`. CoRE IDs are not the chemical identifiers MOFid-v1/v2 and
their short topology field is not the full CrystalNets evidence used for
grouping. Quote bracket-containing filenames in shell commands and percent-
encode them in URLs. Historical-to-current crosswalks are private and are
not part of the public release.

## Release files and links

| Resource | Structure records | Distribution route |
| --- | ---: | --- |
| Full metadata/JSON archive, without CIF bytes | 42,574 | Zenodo |
| COD CIF archive | 19,598 | Zenodo |
| SI CIF archive | 5,727 | Zenodo |
| Modified CSD CIF package | 13,001 | Local review, then CCDC |
| Unmodified CSD CIF package | 4,248 | Local review, then licensed CCDC access |

The two CSD categories together account for the 17,249 CSD-derived structures.
They preserve the release-recorded modified/unmodified classification, not a
new source-CIF comparison or permission grant. Category counts and package
hashes come from their verified manifests, not ASR/FSR/ION labels. Each local CSD
package has its own membership manifest, checksum ledger and access/licensing
document. It must not be presented as a
COD/SI overlay or a complete loader-ready release.

The local archives are `CoRE-MOF-COD_CSD_modified_cifs_20261001.zip` and
`CoRE-MOF-COD_CSD_unmodified_cifs_20261001.zip`. Their private membership records
are `manifests/ccdc_package_manifest.json` and
`manifests/classification_manifest.csv` inside each archive. The public catalogue
records aggregate package counts and hashes, not structure-resolved CSD membership.

Current release links and DOIs below are intentionally blank for the authors to
fill after release. A local package is not an uploaded or published collection.

| Current release resource | Link | DOI |
| --- | --- | --- |
| CoRE-MOF-COD Zenodo deposit | | |
| CoRE-MOF-Tools Zenodo software archive | | |
| Modified CSD collection at CCDC | | |
| Unmodified CSD collection at CCDC | | |

Publication still requires the existing scientific admission and
asset-specific permissions checks. A blank URL is not permission to guess a
provider endpoint or promote a prepared catalogue to a published release.

## Metadata-only use

After verifying an authorized archive against its independently supplied hash,
extract it into a new private directory. The complete metadata-only layout can
be loaded without CIF bytes:

```python
from CoREMOF.dataset import CoREMOFDataset

dataset = CoREMOFDataset.from_release(
    "/authorized/path/CoRE-MOF-COD", verify_cif_files=False
)
view = dataset.classify("5checker")
print(dict(view.label_counts()))
```

Use `examples/read_release_metadata.py` for explicit version checking, a
read-only summary and optional verification against a trusted metadata-ledger
hash. Loading metadata does not verify absent CIF bytes, execute a checker,
change targets or promote a staged release. `verify_cif_files=True` requires
all CIFs in the loaded manifest.

## Source-only data

Do not treat a source JSONL slice or source-only CIF ZIP as a full release.
Use the existing `examples/export_source_projection.py` to create a separately
authenticated contract from the authorized complete release. Then use:

```python
dataset = CoREMOFDataset.from_projection(
    "/authorized/path/source-root", "/authorized/path/source-projection.json",
    expected_sha256="<independently-received-contract-SHA256>",
    verify_cif_files=False,
)
```

The contract retains full-release grouping relationships through omitted
sources. It does not grant source-data redistribution permission.

## Download catalogue formats

The database repository's **archive publication catalogue** describes ZIP
assets, licences, release/permission status, DOIs and mirror URLs. Its own
standard-library downloader verifies and safely extracts an approved asset.
Prepared/unpublished assets are not network-downloadable through that route.

The existing `CoREMOF.retrieval.fetch_release` API and
`examples/fetch_release.py` use a **different, file-level catalogue** with schema
`coremof-release-catalog/1.0`. It declares each file's path, URL, size and hash
and requires an independently supplied catalogue SHA-256. It remains available
for authorized file-level providers. Do not pass the archive publication
catalogue into this API or invent a hosted catalogue. No existing retrieval
API, input contract or defaults are changed by these access examples.

## New analyses versus frozen experiments

For a new target-complete benchmark, use `examples/build_target_first_benchmark.py`:
join by exact structure ID, preserve zeros and remove missing required targets
before cohort construction and grouped splitting. Magnitudes never guide
grouping, diversity or assignments. Build full-release relationships before
filtering, even when only one source will be used.

The paper's frozen adsorption experiment has assignment SHA-256
`9e72992970518d039f9631b1f45b516f4ff3603f7945bdcb82dbad283529dcbd`
and train/val/test = 3,737/466/468. Use its original handoff/checker/feature
revision for `examples/replay_common_input_benchmark.py`; a newer descriptive
metadata catalogue must not replace frozen evidence or saved predictions.

Current data candidates still record provisional membership, `STAGE_ONLY`
MOFid and a blocked publication gate. They can be inspected in authorized
private analyses, but documentation edits do not make them final releases.
New grouping/cohort outputs remain exploratory (`official_split=false`).

## Structure permissions

COD, SI, modified CSD and unmodified CSD structures need separate permission
records. A modified CSD classification records a supported curation change,
such as solvent removal under the agreed classification, but is not itself a
redistribution licence. Unmodified CSD access requires the applicable valid CSD
licence. The access terms for the modified collection must be confirmed for
that collection, rather than inherited from an earlier collection. Neither
ASR/FSR/ION nor CR/NCR nor the scientific OA category proves redistribution
permission. Both CSD categories are excluded from the GitHub/Zenodo CIF
archives. Curation/descriptor calculations need the actual CIFs and the user's
applicable software/source permissions. Reading precomputed checker results
needs neither checker engines nor a CCDC installation.

This explanation belongs in access/repository documentation, not as an
implementation-status discussion in the manuscript.
