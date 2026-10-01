# Local code-review and upload guide

This source snapshot contains the updated CoRE-MOF-Tools `0.4.0.dev0` code,
examples, tests and documentation from the verified September 26 results-only
checkpoint. Its scientific source files were compared with the maintained
checkout. Original edits, Git history and prior artifacts remain untouched.

The package reads the five checker results and combines user-selected votes.
It does not distribute the third-party checker algorithms, reference tables or
execution workers. Compatibility entry points report the removal explicitly.
The optional own-group MOFClassifier integration remains separate.

## Separate software and database releases

This repository is the **code-only Python tools** contribution. Its package
version is `0.4.0.dev0`; a local source/wheel build is not a published PyPI or
Zenodo release. Database archives and the manuscript/plot-input transfer do not
belong inside the Python package.

The database publication uses full metadata for 42,574 structures, without CIF
bytes, and separate COD/SI CIF archives for 19,598/5,727 structures on Zenodo.
The separate modified CSD and unmodified CSD local packages have been built
with 13,001 and 4,248 structures, respectively, for later CCDC handoff. Their
archive names are `CoREMOF-COD_CSD_modified_cifs_20261001.zip` and
`CoREMOF-COD_CSD_unmodified_cifs_20261001.zip`. Neither CSD package may be added to this repository,
GitHub Release assets or Zenodo. Use each package's manifest and confirmed
access terms, not ASR/FSR/ION or CR/NCR labels, to identify its category.

The authors will fill these deliberately blank release links and DOIs later.
Existing historical DOI citations are not replacements for the current links.

| Current release resource | Link | DOI |
| --- | --- | --- |
| CoRE-MOF-Tools software archive on Zenodo | | |
| CoREMOF-COD metadata and COD/SI CIF deposit on Zenodo | | |
| Modified CSD CIF collection at CCDC | | |
| Unmodified CSD CIF collection at CCDC | | |

Keep all scientific admission, licence and upload checks in effect. Preparing
the local packages does not claim that any collection is deposited, published
or authorized for wider distribution. See `README_DATABASE_ACCESS.md` for the
metadata-only loader, source-projection and separate catalogue contracts.

## Installation and examples

```bash
python -m pip install .
python -m CoREMOF doctor
python examples/read_checker_results.py /authorized/path/CoREMOF-COD
python examples/build_target_first_benchmark.py --help
python examples/replay_common_input_benchmark.py --help
```

The latter two workflows are different: one creates a new target-complete
exploratory experiment; the other reproduces frozen paper assignments without
recalculating scientific features. The code supports both without changing
legacy API defaults. See the matching documentation under `docs/source/`.

## Local assets versus Git

The original local preparation retains the legacy lookup tables, node archive
and historical predictor assets separately for authorized local use. They are
not included in this code-only fork. `LOCAL_ASSETS.json` records their identities,
not a promise that the files exist in a clone. Previously built asset-bearing
wheels and source archives also remain outside this contribution.

A code-only clone supports the lightweight release-loading/classification/split
API. Legacy table lookup and optional predictors need their separately obtained
authorized assets. Their absence must not be replaced by fabricated tables or
different model weights. Building a source archive or wheel from the local
asset-bearing tree can include those assets, so the resulting distributions are
**not** cleared for upload merely because the code passed `verify_upload.py`.

The existing licence and third-party notices are retained. Read
`THIRD_PARTY_NOTICES.md` before distributing optional data/model assets. The
approved missing-MOFid policy is documented in `README.md`: unresolved values
remain unavailable, and never create a grouping match. The policy is not proof
that a candidate release has been promoted or that all data may be redistributed.

No cluster-specific agent skill, Git history, cache, credential or old checker
recovery directory is copied into this Git surface. The reviewed export contains
no Git metadata and does not configure a remote. In a fork checkout, retain the
existing fork remote. Run `python3 verify_upload.py` before adding files and
inspect the staged diff. The check includes both new Git-visible files and
already tracked files, even if `.gitignore` now excludes them. It also compares
the current code surface with `GIT_UPLOAD_MANIFEST.json`.

When updating a fork, copy the contents of this folder into the fork checkout,
not a nested `02_CoRE-MOF-Tools/` directory. Preserve unrelated work and inspect
removed files as well as additions. Copying files on top of an old checkout
does not remove retired checker workers or untrack private assets. Do not push
the original development checkout wholesale: its existing Git index includes
legacy database tables, structure archives and predictor weights. This check
inspects current files and the Git index, not earlier commits. Creating a new
branch or deleting a file does not remove it from existing history. Review
the fork's existing history before publishing an asset-bearing repository.
Do not rewrite that history or delete unrelated work as part of copying this
code-only export.

For intended edits after this review, inspect them first, then refresh the
manifest with `python3 verify_upload.py --write-manifest` and rerun the default
check. Do not refresh it to bypass a reported private-data or checker-code
violation. Temporary environments, build outputs and local assets stay out of
the contribution.
The existing portable benchmark skill is retained under
`.agents/skills/coremof-release-curation/`, with its local-asset note aligned
to this layout. It does not carry scheduler or licensed-node instructions.
