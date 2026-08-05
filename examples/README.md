# Examples

Run `coremof doctor` before using an example. The notebooks are grouped by the external requirements they need.

The dataset tools have two standard-library-only entry points:

- [`CoREMOF_dataset_splitting_quickstart.ipynb`](CoREMOF_dataset_splitting_quickstart.ipynb)
  is the interactive target/feature merge, CR/NCR classification, and
  parent-aware splitting tutorial.
- [`screen_candidates.py`](screen_candidates.py) is an automatic,
  command-line high-throughput screen. It applies eligibility filters before
  ranking, excludes null, non-numeric, and non-finite ranking values, and
  writes a ranked CSV plus a hash-bound JSON receipt.

Ranking is precision-safe: integer inputs remain exact integers, and numeric
text is compared as an exact decimal rather than being converted through a
binary64 float. Native float inputs retain their original IEEE-754 value.
This means, for example, `9007199254740993` correctly ranks above
`9007199254740992`; exact numeric ties are resolved by ascending
`structure_id`. For exact integer target columns, declare `value_types` as
`"int"` in the target configuration. For exact decimal text, leave the target
untyped or declare it as `"string"`; declaring it as `"float"` intentionally
uses normal binary64 parsing. Booleans, nulls, NaN, and infinities are never
ranked.

When `--split` is used below, `priority_main`, `main_union`, `auto`, and
`parent_only` are project-defined CoREMOF-tools identifiers rather than
standard scientific terms. A **release-authorized parent triad** is a
validated `status/group/size` triple declared by the release. `MATCHED` means
available with at least two observed members, `UNMATCHED` means available with
one member, and `NOT_AVAILABLE` supplies no edge; unavailable rows never match.
Here **exact RAC5** means all 264 ordered finite binary64 values match after
`-0.0` is mapped to `+0.0` and `float.hex()` is applied, with
`rtol=atol=0` and no scaling/imputation. **Exact MOFid** means the complete
string matches after Unicode-whitespace collapse, trim, whole-field
missing/execution-placeholder rejection, Unicode NFKC, and case-folding, in
that order; the case-insensitive rejected tokens are empty, `-`, `nan`,
`none`, `null`, `n/a`, `na`, `unknown`, `missing`, `timeout`, `timed out`,
`error`, `failed`, `fail`, `fail process`, `failed process`, and
`process failed`. It is not partial/fuzzy matching and does not alter a CIF.

`priority_main` builds an explanatory hierarchy over the full release: it
reads `rac_status/rac_group/rac_size`,
`mofid2_status/mofid2_group/mofid2_size`, and
`mofid1_status/mofid1_group/mofid1_size` from
`parent_groups/parent_groups.csv`; exact release-authorized RAC5 groups anchor
components, then exact MOFid-v2 and MOFid-v1 groups are processed in that
order. A lower group touching no stronger component creates a component; one
touching exactly one attaches only its unresolved rows; one touching two or
more never merges them, records `PARENT_METHOD_CONFLICT`, and leaves
lower-only rows unresolved. Missing evidence becomes a unique singleton by
default; this is not row-wise first-nonmissing selection.

`main_union` is the separate leakage guard, not a parent claim. It reads full
CIF SHA-256 values from `manifests/cif_manifest.csv` and database-namespaced
source-sibling/RAC5/MOFid
group/status/size columns from the parent table; it forms transitive components
from those five exact relations before filtering. Missing CIF hashes fail
closed, missing optional evidence adds no edge, and all available edges are
unioned without precedence. `--leakage-guard auto` is only a selector: it
chooses `main_union` for `priority_main` and `parent_only` for an explicit
direct/reference method. `parent_only` uses only that selected explanatory
grouping as blocks and adds no cross-method edges. Neither guard adds Zeo++,
topology, common-name, provisional source-ID/MOFid transitive groups, or
StructureMatcher relations.

Rank a numeric release-metadata field directly:

```bash
python examples/screen_candidates.py /path/to/coremof_v26.0.2 \
  --rank-by cell_volume_A3 \
  --source COD \
  --metal Cu \
  --output-directory screening/copper
```

To rank a target or joined feature, pass the same target-merge configuration
used by the notebook. Required-target filtering happens before ranking and,
when requested, before split assignment:

```bash
python examples/screen_candidates.py /path/to/coremof_v26.0.2 \
  --target-config targets.json \
  --rank-by xe_uptake \
  --require-target xe_uptake \
  --checkers 5checker \
  --label CR \
  --source COD \
  --variant ASR \
  --metal Cu \
  --order descending \
  --limit 1000 \
  --split \
  --parent-method priority_main \
  --leakage-guard auto \
  --output-directory screening/xe_uptake
```

Repeat `--source`, `--variant`, `--metal`, `--label`, or `--require-target` to
select multiple values. `--required-target-mode any` accepts rows with any one
of several required targets; the default is `all`. With `--split`, the script
uses the package's full-release parent graph and emits assignment CSV/JSON
alongside the ranking outputs. Such user-generated splits remain exploratory,
not official release splits. Each required target is a separate deterministic
CSV column named `required_target:<target-name>`; it is not hidden inside a
JSON cell. The receipt maps target names to those columns and hashes the exact
example script, CoREMOF classification, target-merge, and optional
parent-splitting source files used by the run. It also binds the release
inputs, target/config/alias/feature inputs through the embedded merge receipt,
ranked CSV, and optional split outputs. When splitting is enabled, the
top-level receipt also includes the explanatory parent definition plus both
the requested and resolved leakage guard and the resolved guard definition.
The script and loaded classification sources are frozen at module import and
rechecked before output publication; if any changes during the run, the script
fails and must be restarted instead of writing a misleading receipt.
Create-if-absent publication is race-safe. `--overwrite` is deliberately a
single-writer operation; serialize concurrent writers for the same directory
and stem with an external lock, or use a new immutable stem for each run.

The historical local files named `henry.txt` contain dimensionless average
Rosenbluth weights. They are not Henry coefficients, uptake, or selectivity;
ranking them in descending order is a workflow demonstration only.

| Directory | Purpose | Main requirements |
|---|---|---|
| `checker/` | Structural validation examples | MOFChecker, MOFClassifier; CSD API for MOSAEC |
| `curation/` | Database download, lookup, and CIF curation | ASE, pymatgen, gemmi; CSD API for CSD downloads |
| `features/` | Geometric, topology, OMS, RAC, and heat-capacity features | Zeo++, juliacall/CrystalNets, molSimplify as applicable |
| `ion_pacman/` | PACMAN charge-based ion-preserving curation | PACMAN-charge |

Examples write outputs relative to the current working directory. Copy an example to a separate working directory if you want to preserve the checked-in reference outputs.
