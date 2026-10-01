# Separate database releases and Python tools

CoRE-MOF-Tools contains the API and examples. The separate **CoREMOF-COD**
database repository contains release documentation, schemas and versioned data
catalogues. Large metadata/JSON/CIF bundles belong in approved Zenodo deposits
and, optionally, matching GitHub Release assets, not in this Python source tree.
The repository URL and data DOI will be linked when confirmed.

## Metadata-only use

After verifying an authorized archive against its independently supplied hash,
extract it into a new private directory. The complete metadata-only layout can
be loaded without CIF bytes:

```python
from CoREMOF.dataset import CoREMOFDataset

dataset = CoREMOFDataset.from_release(
    "/authorized/path/CoREMOF-COD", verify_cif_files=False
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

COD, SI, previously licensed modified-CSD packages and newly processed CSD
structures need separate permission records. Neither ASR/FSR/ION nor the OA
scientific category proves redistribution permission. Unmodified CSD CIFs are
not publicly bundled. Curation/descriptor calculations need the actual CIFs
and the user's applicable software/source permissions. Reading precomputed
checker results needs neither checker engines nor a CCDC installation.

This explanation belongs in access/repository documentation, not as an
implementation-status discussion in the manuscript.
