# Local code-review and upload guide

This source snapshot contains the updated CoRE-MOF-Tools `0.4.0.dev0` code,
examples, tests and documentation from the verified September 26 results-only
checkpoint. Its scientific source files were compared with the maintained
checkout. Original edits, Git history and prior artifacts remain untouched.

The package reads the five checker results and combines user-selected votes.
It does not distribute the third-party checker algorithms, reference tables or
execution workers. Compatibility entry points report the removal explicitly.
The optional own-group MOFClassifier integration remains separate.

## Installation and examples

```bash
python -m pip install .
python -m CoREMOF doctor
python examples/read_checker_results.py /authorized/path/coremof_v26.0.2
python examples/build_target_first_benchmark.py --help
python examples/replay_common_input_benchmark.py --help
```

The latter two workflows are different: one creates a new target-complete
exploratory experiment; the other reproduces frozen paper assignments without
recalculating scientific features. The code supports both without changing
legacy API defaults. See the matching documentation under `docs/source/`.

## Local assets versus Git

The local copy includes the recorded legacy lookup tables, node archive and
historical predictor assets needed to preserve the source snapshot. These are
explicitly ignored by Git pending asset-level permissions. `LOCAL_ASSETS.json`
lists them. `local/artifacts/` retains the previously tested wheel and source
archive unchanged, for authorized local installation only.

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
recovery directory is copied into this Git surface. No remote is configured.
Run `python3 verify_upload.py` before adding files and inspect the staged diff.
The existing portable benchmark skill is retained under
`.agents/skills/coremof-release-curation/`, with its local-asset note aligned
to this layout. It does not carry scheduler or licensed-node instructions.
