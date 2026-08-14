# CoREMOF-tools dataset classification and splitting handbook

This handbook explains how to use the lightweight CoREMOF-tools API to:

- load and validate a CoRE-MOF release;
- obtain CR, NCR, AMBIGUOUS, and UNCHECKED labels from a selected checker set;
- filter structures by source, variant, metal, label, or public ID;
- merge one or more target files with current metadata and feature tables;
- inspect criterion-specific parent groups; and
- create deterministic, parent-aware train/validation/test splits.

The examples target CoREMOF-tools `0.4.0.dev0` and the v26 release layout.
The API uses only the Python standard library. It reads completed release
metadata; it does not run crystallographic checkers or calculate descriptors.
Supported Python versions are 3.9, 3.10, and 3.11 (`>=3.9,<3.12`).

> **Development-build notice:** version `0.4.0.dev0` is currently provided by
> the validated local wheel/source checkout described below. The public PyPI
> package may not yet contain this API. Verify `coremof --version` before using
> the handbook.

The companion notebook
[examples/CoREMOF_dataset_splitting_quickstart.ipynb](examples/CoREMOF_dataset_splitting_quickstart.ipynb)
runs the main workflow end to end and writes a verified assignment CSV and
reproducibility receipt.

## Contents

1. [What the package does](#1-what-the-package-does)
2. [Installation](#2-installation)
3. [Required release files](#3-required-release-files)
4. [Five-minute quick start](#4-five-minute-quick-start)
5. [Loading and validating a release](#5-loading-and-validating-a-release)
6. [Reading structure records](#6-reading-structure-records)
7. [CR/NCR classification](#7-crncr-classification)
8. [Filtering structures](#8-filtering-structures)
9. [Parent-group criteria](#9-parent-group-criteria)
10. [Leakage guards](#10-leakage-guards)
11. [Merging target results with release features](#11-merging-target-results-with-release-features)
12. [Creating train/validation/test splits](#12-creating-trainvalidationtest-splits)
13. [Understanding `SplitResult`](#13-understanding-splitresult)
14. [Writing CSV and JSON outputs](#14-writing-csv-and-json-outputs)
15. [Command-line usage](#15-command-line-usage)
16. [Common workflows](#16-common-workflows)
17. [Reproducibility and scientific interpretation](#17-reproducibility-and-scientific-interpretation)
18. [Licensing and redistribution](#18-licensing-and-redistribution)
19. [Troubleshooting](#19-troubleshooting)
20. [Current validation status and limitations](#20-current-validation-status-and-limitations)
21. [Developer checks](#21-developer-checks)

## 1. What the package does

The splitting package consumes an extracted CoRE-MOF release. It independently
recomputes checker-consensus labels from the public checker status columns and
uses the release's audited parent-group tables to prevent related structures
from being separated unintentionally.

It provides three related operations:

1. **Classification:** decide whether a structure is CR, NCR, AMBIGUOUS, or
   UNCHECKED under a selected checker view.
2. **Dataset splitting:** place eligible structures in train, validation, and
   test partitions while keeping a selected structural relation—or a broader
   leakage component—inside one partition.
3. **Target preparation:** join user endpoints and current feature tables,
   preserve declarations/provenance, and explicitly exclude missing targets
   before split assignment.

The package does **not**:

- run MOFClassifier, MOFChecker, Chen-Manz, MOSAEC, or SETC-GAT;
- calculate RACs, MOFid, Zeo++, topology, or open-metal-site properties;
- repair or rewrite CIF files;
- impute missing checker or parent information;
- convert an execution error into a checker FAIL; or
- create an official CoRE-MOF benchmark assignment.

These boundaries are intentional. The package consumes validated scientific
evidence without silently changing the method that produced it.

Keep five choices conceptually separate:

| Choice | Question it answers |
|---|---|
| Checker view | How is CR/NCR/AMBIGUOUS/UNCHECKED recomputed? |
| Parent method | Which relation explains structural grouping? |
| Leakage guard | Which known relations must not cross partitions? |
| Required targets | Which rows contain the endpoints needed by this model? |
| Filters | Which labelled rows are eligible for this experiment? |

### Project-defined identifiers used in this handbook

`priority_main`, `main_union`, `auto`, and `parent_only` are CoREMOF-tools API
identifiers created for this workflow; they are not community-standard
crystallographic terms.

At their first use below, **release-authorized and available** mean that a
criterion's `status/group/size` triad is declared by
`parent_group_methods.json`, appears in `parent_groups.csv`, and passes loader
validation. `MATCHED` means available with an observed group size of at least
two; `UNMATCHED` means available with size one; `NOT_AVAILABLE` means the
criterion supplies no edge and its size-one group exists only to preserve the
table shape. The recorded size must equal the number of rows carrying the
group. Missing/unavailable rows never match and become unique singletons or
explicit exclusions.

At their first use below, **exact RAC5** means equality of all 264 ordered
finite values after IEEE-754 binary64 parsing, mapping `-0.0` to `+0.0`, and
`float.hex()` serialization, with `rtol=atol=0` and no scaling, feature
deletion, imputation, or rounding. **Exact MOFid** means equality of the
complete release-authorized current string after converting to text;
collapsing each
Unicode-whitespace run to one ASCII space; trimming; case-insensitively
rejecting an empty field or `-`, `nan`, `none`, `null`, `n/a`, `na`,
`unknown`, `missing`, `timeout`, `timed out`, `error`, `failed`, `fail`,
`fail process`, `failed process`, or `process failed`; applying Unicode NFKC;
and case-folding, in that order. It is never partial/fuzzy matching and never
changes CIF bytes, atoms, bonds, occupancies, coordinates, chemistry, or unit
cells.

`priority_main` is the conflict-aware **explanatory parent hierarchy**. It is
computed on the complete release from the parent table's
`rac_status/rac_group/rac_size`,
`mofid2_status/mofid2_group/mofid2_size`, and
`mofid1_status/mofid1_group/mofid1_size` fields in this exact order:

1. every available release-authorized exact RAC5 group seeds a component;
2. each exact MOFid-v2 group attaches its unresolved rows to zero or one
   stronger component, or forms a new MOFid-v2 component if it touches none;
3. exact MOFid-v1 groups apply the same rule last; and
4. a structure still lacking all three inputs becomes its own unique singleton
   unless `missing_parent="exclude"` was explicitly selected.

The word **priority** in `priority_main` means the RAC5 → MOFid-v2 → MOFid-v1
precedence used to explain parent groups. It does not rank, schedule, or
recalculate failed scientific features; calculation-retry scheduling is a
separate curation-operation queue outside this package.

If a MOFid group touches two or more stronger components, those components are
never merged. The package records `PARENT_METHOD_CONFLICT`; a lower-only row is
unresolved under the explanatory hierarchy. `priority_main` does **not** use
Zeo++, CrystalNets topology, source ID, common name, CIF hash, provisional
source-ID/MOFid transitive groups, or StructureMatcher evidence.

`main_union` is a distinct conservative **leakage guard**, not an explanatory
parent claim. It forms transitive connected components over the complete
unfiltered release from the full `sha256` values in
`manifests/cif_manifest.csv` and source/RAC5/MOFid group/status/size columns in
`parent_groups.csv`.
Its five relations are exact full CIF SHA-256 equality, database-namespaced
source siblings, and available release-authorized RAC5, MOFid-v2, and MOFid-v1
group edges. The source group compares the ordered `(source_database,
source_id)` pair after applying the text procedure above to each field
separately, so equal IDs in different databases do not match. A missing full CIF hash fails closed. Missing optional evidence
adds no edge, common nulls never match, and there is no priority or conflict
resolution: every available listed edge is unioned. It can place two different `priority_main` groups in the same
indivisible split block without claiming they have the same parent.
Normal release loading requires the source triad, so it never needs metadata
fallback. A manually constructed compatibility object with no source-criterion
fields for a row instead strips and uppercases `source_database` (using
`UNKNOWN` only if blank), and strips, case-folds, Unicode-whitespace-splits,
and ASCII-space-joins `source_id`; it does not apply NFKC or reject a nonblank
placeholder. An explicit source `NOT_AVAILABLE` never falls back.
`leakage_guard="auto"` is only a selector: it chooses `main_union` for
`priority_main` and `parent_only` for an explicitly selected direct/reference
parent method. `parent_only` uses only that selected explanatory grouping as
split blocks and adds no cross-method edges. An unresolved `priority_main`
conflict is excluded with `parent_only`; with `main_union`, it can remain safely
assigned inside the broader block while retaining `PARENT_METHOD_CONFLICT`.
All label, source, variant, metal, ID, and target filters are applied only
after the full-release grouping step.

## 2. Installation

### 2.1 Install a development build

Until version `0.4` is published, clone this repository and install it in a
clean environment:

```bash
python3 -m venv coremof-split-env
source coremof-split-env/bin/activate
git clone https://github.com/Chung-Research-Group/CoRE-MOF-Tools.git
cd CoRE-MOF-Tools
python -m pip install .
coremof --version
```

If an audited wheel has been supplied separately, install it with
`python -m pip install /path/to/coremof_tools-0.4.0.dev0-py3-none-any.whl`.

The distribution name is `CoREMOF-tools`, the Python import namespace is
`CoREMOF`, and the console command is `coremof`. The base installation has no
third-party runtime dependencies.

### 2.2 Install from the source checkout

For development:

```bash
cd /path/to/CoRE-MOF-Tools
python -m pip install -e .
```

### 2.3 Optional historical scientific features

The release loader and splitter do not need NumPy, pandas, scikit-learn,
CCDC, Zeo++, MOFid, or CrystalNets. Other historical CoREMOF-tools features do
need scientific dependencies. Until `0.4` is published, install them from the
source checkout so an older PyPI release is not substituted:

```bash
cd /path/to/CoRE-MOF-Tools
python -m pip install -e ".[full]"
```

Run the feature-level installation check with:

```bash
coremof doctor
```

Missing optional components do not prevent release loading, classification,
or splitting.

## 3. Required release files

Pass the **extracted release root**—the directory containing
`dataset_info.json`—to the loader. Do not pass a `.tar`, `.zip`, or the outer
directory that merely contains the extracted release. A normal release tree
looks like:

```text
coremof_v26.0.2/
├── dataset_info.json
├── cifs/
│   ├── ASR-COD-2000-0001.cif
│   └── ...
├── metadata/
│   └── metadata.csv
├── parent_groups/
│   ├── parent_groups.csv
│   └── parent_group_methods.json
└── manifests/
    └── cif_manifest.csv
```

The loader requires:

- `dataset_info.json`;
- `metadata/metadata.csv`;
- `parent_groups/parent_groups.csv`; and
- `parent_groups/parent_group_methods.json`.

`manifests/cif_manifest.csv` is optional for basic metadata inspection, but it
is effectively required for the recommended project-defined `priority_main`
split because its automatic project-defined `main_union` leakage guard uses
full CIF SHA-256 identity as one of the five edge types defined above.

The loader validates the exact structure-ID sets across files, naming/source
consistency, checker label recomputation, parent group status/size consistency,
CIF paths, manifest sizes, and full SHA-256 syntax. It fails closed when the
release contract is inconsistent.

## 4. Five-minute quick start

```python
from CoREMOF.dataset import CoREMOFDataset

release = "/path/to/coremof_v26.0.2"

# Load and validate release tables.
dataset = CoREMOFDataset.from_release(release)
print(dataset.dataset_version, len(dataset))

# Recompute the strict five-checker consensus.
classified = dataset.classify(checkers="5checker")
print(dict(classified.label_counts()))
print("CR:", len(classified.cr_ids))
print("NCR:", len(classified.ncr_ids))

# Build a COD+SI CR/NCR modelling split.
split = classified.train_valid_test_split(
    parent_method="priority_main",
    leakage_guard="auto",
    labels=("CR", "NCR"),
    sources=("COD", "SI"),
    fractions=(0.8, 0.1, 0.1),
    random_state=42,
)

print(dict(split.counts))
print(dict(split.leakage_audit))
split.write("model_splits", stem="cod_si_5checker_seed42")
```

The output directory will contain:

```text
model_splits/
├── cod_si_5checker_seed42.csv
└── cod_si_5checker_seed42.json
```

The CSV is the assignment table. The JSON is the reproducibility receipt.

### Configure extracted release roots

Set the variables to the extracted release directories on your machine:

```bash
export COREMOF_V2601_RELEASE=/path/to/coremof_v26.0.1
export COREMOF_V2602_RELEASE=/path/to/coremof_v26.0.2
```

Then run, for example:

```bash
coremof split "$COREMOF_V2602_RELEASE" \
  --checkers 5checker \
  --sources COD SI \
  --output-directory "$HOME/coremof_splits" \
  --stem cod_si_5checker_seed42
```

## 5. Loading and validating a release

### 5.1 Normal validation

```python
from CoREMOF.dataset import CoREMOFDataset

dataset = CoREMOFDataset.from_release("/path/to/coremof_v26.0.2")
```

Normal loading verifies table contents and manifest declarations, but it does
not prove that every CIF currently exists or matches the declared bytes. This
is the appropriate default for repeated interactive work with a previously
audited local release. Use strict verification for a newly obtained, copied,
or extracted tree.

### 5.2 Strict CIF-byte verification

Use strict verification after downloading, copying, extracting, or archiving a
release:

```python
dataset = CoREMOFDataset.from_release(
    "/path/to/coremof_v26.0.2",
    verify_cif_files=True,
)
assert dataset.cif_files_verified
```

This reads every CIF, checks its byte size, calculates its full SHA-256, and
compares both with the manifest. It is more I/O-intensive but catches damaged,
missing, or substituted files.

`split_release(..., verify_cif_files=True)` can perform the same check when
given a release path. It cannot retroactively verify a dataset/classified
object loaded without strict verification; load that object with
`CoREMOFDataset.from_release(..., verify_cif_files=True)` first.

### 5.3 Basic dataset properties

```python
print(dataset.dataset_version)
print(len(dataset))
print(dataset.release_root)
print(dataset.structure_ids[:5])
print(dict(dataset.input_hashes))
```

Important attributes include:

| Attribute | Meaning |
|---|---|
| `dataset_version` | Release version declared by `dataset_info.json` |
| `release_root` | Resolved local release path |
| `structure_ids` | Structure IDs in release metadata order |
| `records` | Immutable sequence of joined structure records |
| `metadata_rows` | Immutable metadata row views |
| `dataset_info` | Parsed release-level metadata |
| `parent_group_methods` | Parsed parent-method contract |
| `parent_by_id` | Raw parent CSV rows indexed by structure ID |
| `cif_hashes` | Full manifest SHA-256 values indexed by structure ID |
| `input_hashes` | Hashes of release tables used by the loader |
| `cif_files_verified` | Whether CIF bytes were rehashed during loading |

## 6. Reading structure records

Access a structure by public ID or by its release position:

```python
structure_id = dataset.structure_ids[0]
record = dataset[structure_id]
first_record = dataset[0]
first_ten = dataset[:10]
```

Read identity and metadata fields:

```python
print(record.structure_id)
print(record.source_database)
print(record.structure_variant)
print(record.metal_elements)
print(record["source_id"])
print(record.get("doi"))
print(dict(record.metadata))
```

Metadata CSV values are exposed as strings. Convert numeric fields explicitly
when needed:

```python
year_text = record.get("publication_year")
publication_year = int(year_text) if year_text and year_text != "UNKN" else None
```

Inspect one criterion-specific parent group:

```python
rac_parent = record.parent_group("rac5")
print(rac_parent.status)    # MATCHED, UNMATCHED, or NOT_AVAILABLE
print(rac_parent.group_id)  # for example R-1A2B3C4D
print(rac_parent.size)
print(rac_parent.available)
print(rac_parent.matched)
```

`StructureRecord.parent_group()` reads current criteria stored in the loaded release table, such as
`rac5`, `mofid_v2`, or `source_id`. Computed `priority_main` and `main_union`
relations are not per-record table fields; inspect them with `ParentResolver`
or `SplitResult`.

Records and release mappings are read-only. Create your own dictionary if you
need a mutable working copy:

```python
working_row = dict(record.metadata)
```

## 7. CR/NCR classification

### 7.1 Official checker views

The package provides three official views:

| View | Included checkers |
|---|---|
| `3checker` | MOFClassifier, original MOFChecker, Chen-Manz |
| `4checker` | 3-checker set plus MOSAEC |
| `5checker` | 4-checker set plus SETC-GAT |

Use either the canonical name or the corresponding integer:

```python
three = dataset.classify(checkers="3checker")
four = dataset.classify(checkers=4)
five = dataset.classify(checkers="5checker")
```

### 7.2 Strict consensus rule

Only `PASS` and `FAIL` are scientific votes.

| Included checker results | Consensus label |
|---|---|
| Every checker is `PASS` | `CR` |
| Every checker is `FAIL` | `NCR` |
| Complete mixture of `PASS` and `FAIL` | `AMBIGUOUS` |
| Any included result is unavailable/non-voting | `UNCHECKED` |

There is no majority vote. An error, timeout, missing input, unsupported
representation, or missing upstream MOSAEC output cannot become NCR.

### 7.3 Reading labels

```python
print(dict(five.label_counts()))
print(five.cr_ids[:5])
print(five.ncr_ids[:5])
print(five.ambiguous_ids[:5])
print(five.unchecked_ids[:5])

one = five[five.structure_ids[0]]
print(one.label)
print(dict(one.checker_statuses))
```

Convenience attributes are:

- `cr_ids`;
- `ncr_ids`;
- `ambiguous_ids`;
- `unchecked_ids`; and
- `ids_for_label("CR")` for an explicit label.

### 7.4 Custom checker views

For a user-defined sensitivity analysis, pass an ordered sequence of canonical
checker names:

```python
custom = dataset.classify(
    checkers=("MOFClassifier", "MOFChecker", "MOSAEC")
)
print(custom.checker_view)
# custom:MOFClassifier+MOFChecker+MOSAEC
```

Accepted checker names are exactly:

- `MOFClassifier`;
- `MOFChecker`;
- `Chen-Manz`;
- `MOSAEC`; and
- `SETC-GAT`.

Custom order is retained in the view identifier. Duplicates, sets, unknown
names, and empty sequences are rejected. A custom result uses the same strict
consensus rule but is marked `USER_DEFINED`, not an official release view.
The CLI accepts only the three official presets; custom checker sequences are
available through Python.

## 8. Filtering structures

Filters can be applied during classification:

```python
cu_cod = dataset.classify(
    checkers="5checker",
    labels=("CR", "NCR"),
    sources=("COD",),
    variants=("ASR", "FSR"),
    metals=("Cu",),
)
```

Or applied to an existing classified view:

```python
classified = dataset.classify("5checker")
cu_cod = classified.filter(
    labels=("CR", "NCR"),
    sources=("COD",),
    variants=("ASR", "FSR"),
    metals=("Cu",),
)
```

`classified.select(...)` is an alias for `classified.filter(...)`.

### 8.1 Filter semantics

| Filter | Behavior |
|---|---|
| `labels` | Exact canonical labels after case normalization |
| `sources` | Case-insensitive match; normal release values are COD/CSD/SI |
| `variants` | Case-insensitive match; normal values are ASR/FSR/ION |
| `metals` | Keep a structure containing at least one requested element |
| `structure_ids` | Exact public-ID membership |

Different filter categories are combined with AND. Values inside the metal
filter are combined with OR:

```python
# CR or NCR, from COD, containing Cu or Zn.
subset = classified.filter(
    labels=("CR", "NCR"),
    sources=("COD",),
    metals=("Cu", "Zn"),
)
```

Unknown structure IDs raise rather than disappearing silently. A filtered view
remembers its selection steps so a later split can report preselection in its
receipt.

At the Python API boundary, each filter accepts either one exact string or an
ordered `list`/`tuple` of exact strings. Sets, frozensets, mappings, byte
strings, generators, and other one-shot or unordered containers are rejected;
the package does not repair their order before hashing a receipt.

Filtering does not recompute or modify labels. It also does not erase a parent
bridge from the complete release when the recommended leakage guard is used.

## 9. Parent-group criteria

Parent criteria describe structural relatedness. They do not change checker
labels.

### 9.1 Recommended hierarchy

Use `parent_method="priority_main"` as the normal explanatory policy:

1. RAC5 is the strongest available anchor.
2. MOFid v2 attaches otherwise unresolved structures when it does not merge
   different stronger components.
3. MOFid v1 is the final fallback under the same conflict rule.

This is a conflict-aware hierarchy, not a naive first-nonmissing value applied
independently to every row.

### 9.2 Available parent methods

In this section, **canonicalized identifier text** has one narrow, exact
meaning. For release-authorized current fields, the release builder converts
the value to text, collapses each Unicode-whitespace run to one ASCII space,
trims it, rejects a whole-field missing/execution placeholder, applies Unicode NFKC,
and then applies Unicode case-folding. The rejected case-insensitive whole
fields are empty text, `-`, `nan`, `none`, `null`, `n/a`, `na`,
`unknown`, `missing`, `timeout`, `timed out`, `error`, `failed`, `fail`,
`fail process`, `failed process`, and `process failed`. A source key is the
exact ordered pair of separately canonicalized `source_database` and
`source_id`, so equal text from different databases cannot match. A MOFid key
is the complete canonicalized release-authorized current MOFid string; there
is no prefix, substring, or fuzzy match. MOFid equality is eligible only when
its release status is `SUCCESS`, `SUCCESS_TOPOLOGY_UNKNOWN`,
`SUCCESS_TOPOLOGY_ERROR`, or `SUCCESS_TOPOLOGY_TIMEOUT`.

The release builder recomputes these keys and their connected components
freshly over every current row of each named release. v26.0.1 is recomputed
over its 36,628 rows; v26.0.2 is recomputed independently over all 42,574
superset rows. It does not seed v26.0.2 with a prior v26.0.1 component or
import an earlier source-ID/MOFid edge. Missing, null, placeholder,
unresolved-reconciliation, ambiguous-node, timeout, error, no-MOF,
unmatched-node, and decomposition-error MOFid values add no edge and never
match through a shared null. This text processing does not modify atoms,
bonds, occupancies, coordinates, unit cells, general chemical punctuation, or
CIF bytes. The splitter consumes the resulting release group/status/size
columns; it does not perform this canonicalization or recompute the relation.

The **RAC5 fingerprint** is exact equality of all 264 finite descriptor columns
listed, in order, at `criteria.rac5.ordered_descriptors` in
`parent_group_methods.json`. Each value is parsed as binary64, must be finite,
has `-0.0` mapped to `+0.0`, and is serialized with `float.hex()` before exact
ordered equality (`rtol=atol=0`); there is no scaling, feature deletion,
imputation, or rounding. The **selected Zeo++ fingerprint** is exact equality of
`density_g_cm3`, `largest_cavity_diameter_A`, `pore_limiting_diameter_A`,
`largest_free_path_diameter_A`, the accessible and nonaccessible N₂ surface
areas in m²/cm³ and m²/g, the accessible and nonaccessible N₂ volumes in
cm³/g and as fractions, and `he_void_fraction`, plus exact N₂ channel dimension
and available bonded-framework periodic dimension. Its numeric tolerance is
zero; its 13 numeric values use the same finite-binary64, `-0.0` to `+0.0`,
`float.hex()` representation; N₂ and He probe radii are 1.655 Å and 1.32 Å.
One invalid numeric value or hard gate makes the fingerprint unavailable. Unit-cell-extensive areas
and volumes, OMS, topology labels, and framework-component counts are excluded.

The **current CrystalNets scientific fingerprint**, where current means the
topology evidence authorized by the loaded release rather than a search for a
newer runtime result, requires `SUCCESS`,
`topology_available=true`, `error=null`, and nonempty complete SingleNodes and
AllNodes subnets. It contains network dimension, interpenetrated-subnet count,
catenation degree, subnet count, top-level single-node/all-node nets and their
agreement, and—for every canonically sorted subnet—single/all agreement plus
each node view's status, dimension, topology key, topology name, and
topological genome. It is canonical sorted-key JSON with complete subnet
projections sorted while duplicate subnets are retained, then SHA-256 hashed.
Interpenetrated-subnet count and catenation degree must each be integers at
least one and equal the observed subnet count; subnet indices must be unique
and contiguous before sorting. Top-level dimension, nets, or agreement may be
null for heterogeneous subnets and that null is retained. Every node view must
have `SUCCESS`, dimension 0--3, and a nonblank topology key; topology name and
genome may be null.
Runtime, CIF paths/hashes, diagnostics, software boilerplate, and original
subnet order are excluded. An incomplete or unsuccessful result supplies no
edge; this is exact fingerprint equality, not a topology-similarity tolerance.

| Method | Role | Interpretation |
|---|---|---|
| `priority_main` | Recommended | Conflict-aware RAC5 → MOFid v2 → MOFid v1 hierarchy |
| `rac5` | Main/direct | Exact equality of the release-authorized depth-5 RAC fingerprint |
| `mofid_v2` | Main/direct | Exact equality of the complete canonicalized release-authorized MOFid v2 value; missing values never match |
| `mofid_v1` | Main/direct | Exact equality of the complete canonicalized release-authorized MOFid v1 value; missing values never match |
| `rac5_zeo` | Reference | Exact combined RAC5 and selected Zeo++ fingerprint |
| `rac5_crystalnets` | Optional non-decisive reference | Reads `rac_crystalnets_status/group/size`; an `RT-` group requires exact equality of the full 264-value finite RAC5 fingerprint plus a complete successful current CrystalNets scientific fingerprint; missing, nonfinite, partial, timed-out, or failed input creates no evidence |
| `mofid_v2_crystalnets` | Optional non-decisive reference | Reads `mofid2_crystalnets_status/group/size`; an `M2T-` group requires exact equality of complete canonicalized MOFid v2 with status exactly `SUCCESS`, `SUCCESS_TOPOLOGY_UNKNOWN`, `SUCCESS_TOPOLOGY_ERROR`, or `SUCCESS_TOPOLOGY_TIMEOUT`, plus that complete CrystalNets fingerprint. The latter two are successful calculated identifiers whose embedded topology qualifier is ERROR or TIMEOUT, not MOFid execution failures. Every other MOFid-v2 status and every incomplete CrystalNets input adds no edge; this reference is provisional whenever the release-authorized MOFid-v2 input is provisional. If the release-authorized MOFid-v2 values change, rebuild the M2T groups before use |
| `structure_matcher_strict` | Optional reference | `SM-` is a convenience connected component of exhaustive composition-compatible pairs whose forward and reverse pinned strict pymatgen 2024.2.8 `ElementComparator` fits both pass with `fit(..., symmetric=True)`, `ltol=stol=0.001`, `angle_tol=0.01`, `primitive_cell=true`, `scale=false`, `attempt_supercell=true`, `allow_subset=false`, `supercell_size=num_sites`, and no ignored species. Direct symmetric edges are authoritative; a component is not duplicate proof. Parser failures, timeouts, and execution errors are `NOT_AVAILABLE` rather than unmatched; this method enters neither `priority_main` nor `main_union` |
| `zeo` | Reference | Exact selected Zeo++ fingerprint |
| `source_id` | Reference | Exact canonicalized `(source_database, source_id)` sibling key; no cross-database match |
| `common_name` | Reference | Exact NFKC/whitespace/case-folded common-name match; sparse and non-unique |
| `identity_union` | Reference | Provisional source-ID/MOFid transitive groups: freshly recompute v26.0.1 over its 36,628 current rows and v26.0.2 independently over its 42,574 current rows from exact canonicalized database-namespaced source-ID, eligible complete MOFid-v2, or eligible complete MOFid-v1 edges, then take connected-component closure; no prior component import and no RAC5, Zeo++, topology, CIF hash, common name, or StructureMatcher input |
| `none` | Control | Treat every structure as an independent singleton |

RAC/Zeo equality is fingerprint equivalence, not proof of a common synthetic
parent. The project-defined API key `identity_union` means exactly the
provisional source-ID/MOFid transitive groups stated in the table. Each group,
and therefore each value counted by `identity_size`, is one transitive
connected component of structures joined by those identifier-equality edges;
it is not a count of edges or identifiers. All three edge types have equal
union status—there is no priority or conflict rule—and one type can bridge
components formed by another. Each named release is rebuilt independently
from its current rows; no earlier component or MOFid edge is imported. A null,
rejected placeholder, or non-success MOFid status adds no edge, so missing or
unresolved identifiers never match. The criterion is not proof that members
are chemically or geometrically identical, and it is excluded from both
`priority_main` and `main_union`. A relation built from stage-only MOFid
evidence remains a staged candidate and must be rebuilt if the authorized
MOFid evidence changes.

The prefixes are criterion labels, not scientific results: `RT-` means an
exact RAC5-plus-current-CrystalNets group, `M2T-` means an exact
MOFid-v2-plus-current-CrystalNets group, and `SM-` means a connected component
of strict direct StructureMatcher edges. Text after each prefix is only a
compact deterministic group digest; it is not a topology name, MOFid, RMSD,
or similarity score. The release builder hashes the criterion name and its
complete comparison key with length-delimited UTF-8 SHA-256, starts with eight
uppercase hexadecimal characters after the prefix, and extends every actually
colliding prefix one character at a time until unique in the v26.0.2 superset.

The optional criteria are accepted only when the release declares and
validates their columns. They are sensitivity relations and do not alter
`priority_main` or `main_union`. A missing, partial, timed-out, or failed
CrystalNets result is not group evidence. `mofid_v2_crystalnets` is also
provisional whenever its MOFid v2 input bundle is provisional; the release
method file is authoritative. If release-authorized MOFid-v2 values change,
rebuild the M2T groups before use. These optional non-decisive reference alternatives do
not enter `priority_main`, whose explanatory precedence remains exact RAC5,
then complete MOFid v2, then complete MOFid v1, and they add no `main_union`
edge.

`structure_matcher_strict` is derived from an audited direct-pair edge ledger
using Python 3.9, pymatgen 2024.2.8, and NumPy 1.26.4. Candidate blocking uses
the parsed structure's `ElementComparator` fractional-composition hash, and
every pair within an equal block is attempted. Routine formula, DOI, source,
RAC5, topology, site count, and cell volume are not exclusionary blocks.
Both directional calls use `fit(..., symmetric=True)`; an undirected edge is
added only when both calls succeed under the settings in the table.
Because a nonzero-tolerance match relation need not be transitive, its group is
a connected component for inspection or an explicitly requested sensitivity
split. It is not proof that every member directly matches every other member.
Consult the release duplicate-evidence ledger for directional normalized RMS
and maximum displacement values, clique status, parser state, and exact frozen
matcher settings. The historical relaxed matcher is intentionally not exposed
as a parent method.

The strict v2 ledger records forward and reverse `normalized_rms` and
`normalized_max_displacement`; each direction divides periodic site
displacement by that direction's `(V/Nsites)^(1/3)`. They are dimensionless and
not distances in angstroms. This is pymatgen's periodic lattice/site matcher;
it does not call the separate `charnley/rmsd` molecular Kabsch program. Its
pinned direct `CifParser` expands declared symmetry, merges generated sites at
`site_tolerance=1e-4`, rounds coordinates near 1/3 or 2/3 at
`frac_tolerance=1e-4`, checks occupancy, sorts the periodic structure, and
preserves disorder. It performs no manual repair, occupancy choice, atom
deletion, or chemistry edit. Parser, timeout, OOM, or matcher errors are
`NOT_AVAILABLE`, not nonmatches. A directional fit disagreement is
`ASYMMETRIC_RESULT`, remains unavailable, and cannot create an edge.

The loader fails closed on `sm_*` columns unless
`parent_group_methods.json` contains the exact
`criteria.structure_matcher_strict` contract. That declaration fixes method
ID `pymatgen_structure_matcher_strict_v2`, method schema
`coremof-structure-matcher-method/2.0`, direct symmetric strict-match edges as
the authoritative evidence, connected components as a non-transitive
convenience view, and exclusion from both `priority_main` and `main_union`.
Its public status policy is exactly `MATCHED`, `UNMATCHED`, or
`NOT_AVAILABLE`; any component touched by an unresolved pair must be projected
as structure-specific `NOT_AVAILABLE` singletons.

The declaration must bind
`parent_groups/structure_matcher_strict_evidence_receipt.json` by full
lowercase SHA-256. The receipt must be a passing
`coremof-structure-matcher-release-adapter-receipt/1.0` record bound to the
exact `parent_groups.csv` bytes and a full strict-pair-ledger SHA-256. It also
records candidate, successful, unresolved, direct-edge, structure, and
`NOT_AVAILABLE` counts. The receipt and method declaration must both state
that historical relaxed output was neither executed nor exposed. The verified receipt becomes a
split input hash, so a split receipt is bound to this optional evidence.

The required criterion object is:

```json
{
  "role": "OPTIONAL_REFERENCE",
  "method_id": "pymatgen_structure_matcher_strict_v2",
  "method_schema_version": "coremof-structure-matcher-method/2.0",
  "authoritative_evidence": "DIRECT_SYMMETRIC_STRICT_MATCH_EDGES",
  "fit_policy": "FIT_SYMMETRIC_TRUE_REQUIRED",
  "component_semantics": "CONNECTED_COMPONENT_CONVENIENCE_DIRECT_EDGES_AUTHORITATIVE",
  "component_completeness_policy": "INCOMPLETE_COMPONENTS_NOT_AVAILABLE_UNIQUE_SINGLETON",
  "public_status_policy": ["MATCHED", "UNMATCHED", "NOT_AVAILABLE"],
  "included_in_priority_main": false,
  "included_in_main_union": false,
  "historical_relaxed_executed": false,
  "historical_relaxed_exposed": false,
  "evidence_receipt": {
    "file": "parent_groups/structure_matcher_strict_evidence_receipt.json",
    "sha256": "<64 lowercase hexadecimal characters>"
  }
}
```

The bound receipt has schema
`coremof-structure-matcher-release-adapter-receipt/1.0` and exactly these
additional fields: `status` (`PASS`), `dataset_version`, `method_id`,
`method_schema_version`, `parent_groups_sha256`,
`strict_pair_ledger_sha256`, `structure_count`, `candidate_pair_count`,
`successful_pair_count`, `unresolved_pair_count`,
`strict_direct_match_edge_count`, `not_available_structure_count`, and
`historical_relaxed_executed` (`false`) and
`historical_relaxed_exposed` (`false`). Counts must be non-negative;
successful plus unresolved pairs must equal candidate pairs, and the receipt's
structure and unavailable counts must agree with the release table.

### 9.3 Directly inspect a parent resolution

```python
from CoREMOF.parents import ParentResolver

resolver = ParentResolver(dataset)
resolution = resolver.resolve("priority_main")

structure_id = dataset.structure_ids[0]
print(resolution.groups[structure_id])
print(resolution.exclusions.get(structure_id))
print(len(resolution.conflicts))
```

To inspect one method without splitting:

```python
rac = resolver.resolve("rac5")
members = rac.members_by_group()
```

`ParentResolver.resolve()` computes the relation on the full release before an
optional `structure_ids=` subset is returned. This prevents a requested subset
from changing the underlying relation.

### 9.4 Missing parent evidence

The default is:

```python
missing_parent="singleton"
```

Each missing value becomes `SINGLETON:<structure_id>`. Missing structures never
match each other merely because they share a null value.

For a complete-case sensitivity analysis:

```python
missing_parent="exclude"
```

The structure is then retained in the output assignment table with the
explicit exclusion reason `MISSING_PARENT_EVIDENCE`.

### 9.5 Parent conflicts

A conflict occurs when a lower-priority MOFid group spans more than one
stronger RAC/MOFid component. `priority_main` does not merge those stronger
anchors.

The result provides two views of the issue:

- `parent_diagnostics`: row-level unresolved conflicts relevant to assignment;
- `parent_conflicts`: structured group-level ledger, including anchor-only
  conflicts that require no row exclusion.

Each ledger entry records:

- the lower method and group;
- stronger components;
- members already assigned to each stronger component;
- all affected members; and
- unresolved members.

## 10. Leakage guards

`parent_method` explains the selected parent relation. `leakage_guard` decides
which known relations must stay together during splitting.

| Guard | Behavior |
|---|---|
| `auto` | Use `main_union` for `priority_main`; otherwise use `parent_only` |
| `main_union` | Keep each full transitive CIF/source/RAC5/MOFid2/MOFid1 component together |
| `parent_only` | Keep only the selected explanatory method's groups together |

The recommended combination is:

```python
parent_method="priority_main"
leakage_guard="auto"
```

`auto` resolves to `main_union` for this case. The union contains exact edges
from:

1. full CIF SHA-256 identity;
2. database-namespaced source siblings;
3. available release-authorized RAC5 groups;
4. available release-authorized MOFid v2 groups; and
5. available release-authorized MOFid v1 groups.

The union is constructed on the complete COD+CSD+SI release **before** label,
source, variant, metal, or structure-ID filters are applied. A filtered-out CSD
or UNCHECKED row can therefore remain an invisible bridge between two selected
COD CR/NCR rows. This is intentional leakage protection.

Use `parent_only` only when the scientific purpose is explicitly to test the
selected relation without the broader conservative guard:

```python
rac_only = classified.train_valid_test_split(
    parent_method="rac5",
    leakage_guard="parent_only",
)
```

## 11. Merging target results with release features

Join uptake, selectivity, stability, or another user-supplied endpoint **before**
creating a split. This makes target availability an explicit eligibility rule
and avoids assigning rows that cannot be used for the intended model. The
target layer still builds leakage components from the complete release, so an
untargeted or filtered structure can remain an invisible bridge between two
eligible structures.

The target API uses only the standard library and accepts CSV, JSON, and JSONL.
It can combine multiple files with different target columns:

```python
from CoREMOF.dataset import CoREMOFDataset
from CoREMOF.targets import AliasRegistry, TargetSource

dataset = CoREMOFDataset.from_release("/path/to/coremof_v26.0.2")

uptake = TargetSource(
    "/path/to/uptake.csv",
    name="xe_uptake_298K_1bar",
    id_column="structure_name",
    target_columns=("uptake",),
    target_names={"uptake": "xe_uptake"},
    value_types={"xe_uptake": "float"},
    units={"xe_uptake": "mol/kg"},
    conditions={
        "xe_uptake": {"temperature_K": 298, "pressure_bar": 1.0}
    },
)
selectivity = TargetSource(
    "/path/to/selectivity.jsonl",
    name="xe_kr_selectivity_298K_1bar",
    id_column="structure_id",
    target_columns=("xe_kr_selectivity",),
    units={"xe_kr_selectivity": "dimensionless"},
    conditions={
        "xe_kr_selectivity": {
            "temperature_K": 298,
            "pressure_bar": 1.0,
            "feed": {"Xe": 0.2, "Kr": 0.8},
        }
    },
)

modelling_data = dataset.merge_targets(
    (uptake, selectivity),
    feature_tables=("rac5", "zeo", "topology"),
)
print(modelling_data.target_columns)
print(modelling_data.target_values(modelling_data.structure_ids[0]))
```

`feature_tables` selects current release tables by name. Supported names are
`rac5`, `zeo`, `zeo_zero_probe`, and `topology`. Historical feature artifacts
are not selected implicitly. Omitting `feature_tables` joins targets only to
the main public metadata and avoids loading large descriptor tables.

### 11.1 Identifier matching and earlier database names

Current public IDs are matched exactly after trimming surrounding whitespace.
There is no fuzzy, case-insensitive, basename, chemical-name, or substring
matching. Unknown IDs fail the complete merge.

Earlier public IDs, earlier CIF names, or v11 identifiers are accepted only
through an explicit local alias registry. The user must name the authorized
alias columns; the package never searches migration fields automatically:

```python
aliases = AliasRegistry(
    "/path/to/audited_structure_name_registry.csv",
    current_id_column="structure_id",
    alias_columns=("previous_public_id", "v11_structure_id", "v11_cif_file"),
)

modelling_data = dataset.merge_targets(
    (uptake, selectivity),
    alias_registry=aliases,
    feature_tables=("rac5", "zeo"),
)
```

The registry file, selected columns, size, and SHA-256 are bound into the merge
receipt. If one alias maps to two current IDs—or collides with a different
current public ID—the merge fails. Output rows use only current public IDs;
the earlier-ID columns are not copied into the modelling table.

### 11.2 Nulls, types, units, conditions, and duplicate values

- JSON `null` remains Python `None`; an empty CSV cell is null by default.
- CSV non-null values remain strings unless `value_types` explicitly requests
  `float`, `int`, `bool`, `string`, or `json` conversion.
- Declared `json` targets may be JSON scalars, objects, or arrays. Objects and
  arrays are retained as recursively read-only values in memory, serialized
  with sorted object keys in CSV/provenance output, and compared for duplicate
  equality by that canonical JSON representation; object key order therefore
  does not create a false conflict.
- `target_names` explicitly maps a generic input column such as `uptake` to a
  canonical endpoint such as `xe_uptake`; declarations use the canonical name
  and provenance retains the input column.
- Units and experimental/simulation conditions are never inferred. Declare
  them per target, or leave them explicitly unspecified.
- A target name denotes one endpoint contract. Reusing it with different
  units, conditions, or declared types fails; use distinct target names for
  scientifically different endpoints.
- Equal duplicate non-null observations are accepted and retain both
  provenance entries. Conflicting non-null values fail the entire merge.
- A null duplicate does not overwrite a non-null observation. No missing
  target is imputed or converted to zero.

Every merged dataset retains per-observation provenance: source label, input
filename and SHA-256, row number, original input identifier, current/alias
resolution, and supplied value.

Preserve the scientific meaning of each source column. For example, a
Rosenbluth weight is a dimensionless observation, not a Henry coefficient or
an uptake. Declare and name it accordingly; the package does not reinterpret
or convert it:

```python
rosenbluth = TargetSource(
    "/path/to/rosenbluth.csv",
    id_column="structure_name",
    target_columns=("rosenbluth",),
    target_names={"rosenbluth": "rosenbluth_weight"},
    value_types={"rosenbluth_weight": "float"},
    units={"rosenbluth_weight": "dimensionless"},
    conditions={"rosenbluth_weight": {"temperature_K": 298}},
)
```

If a later workflow derives Henry or uptake values, publish those as separate
canonical targets with their own method, unit, condition, and provenance.

### 11.3 Require targets before assignment

Classify the merged dataset, then require all selected endpoints:

```python
classified = modelling_data.classify("5checker")
split = classified.train_valid_test_split(
    parent_method="priority_main",
    leakage_guard="auto",
    labels=("CR", "NCR"),
    required_targets=("xe_uptake", "xe_kr_selectivity"),
    required_target_mode="all",
    fractions=(0.8, 0.1, 0.1),
    random_state=42,
)
```

Rows lacking at least one required endpoint are written as `EXCLUDED` with
`MISSING_REQUIRED_TARGET`. Set `required_target_mode="any"` only when a row
with any one of the listed endpoints is useful. The receipt records each
target's full-release availability count, the combined eligible count, target
definitions and source hashes, and confirms that the target filter preceded
assignment while leakage blocks used the full release universe.

### 11.4 Deterministic merged outputs

Write the modelling table before or alongside the split:

```python
csv_path, provenance_path, receipt_path = modelling_data.write(
    "model_inputs",
    stem="xe_kr_targets",
)
```

This creates:

- `xe_kr_targets.csv`: release metadata, requested current feature columns,
  and every target column for every release structure;
- `xe_kr_targets.provenance.jsonl`: one record per supplied observation; and
- `xe_kr_targets.json`: target definitions, source/registry/feature hashes,
  counts, policies, and the deterministic target-value digest.

### 11.5 Automatic high-throughput candidate screening

The standard-library example
[`examples/screen_candidates.py`](examples/screen_candidates.py) automates a
validated screen without changing checker labels or filling missing values:

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

The operation order is fixed: validate the release, recompute the selected
checker view, apply identity filters, require declared targets, exclude null,
non-numeric, and non-finite ranking values, then rank. Equal numeric values
are ordered by ascending `structure_id`. Integer inputs and decimal text are
ranked without converting through binary64, so values above `2**53` do not
collapse into false ties; native Python floats retain their IEEE-754 values.
Booleans, nulls, NaN, and infinities are never rankable. A non-null string such
as `NA` remains available according to the declared target contract but is
excluded if used as a numeric ranking field. Each required target is a flat
CSV column named `required_target:<target-name>`.

With `--split`, only the emitted candidates receive train/validation/test
assignments, but leakage components are still constructed from the complete
release. The output remains exploratory and `official_split=false`. The
screening receipt hashes the script and exact CoREMOF source closure, embeds
the target/config/alias/feature receipt, binds release inputs, hashes the
ranked CSV, and binds both optional split files.

Do not reinterpret a ranking column. In particular, the historical local
`henry.txt` files contain dimensionless average Rosenbluth weights—not Henry
coefficients, uptake, or selectivity. Descending Rosenbluth-weight ranking is
only a workflow demonstration.

## 12. Creating train/validation/test splits

### 12.1 Method on a classified dataset

```python
split = classified.train_valid_test_split(
    parent_method="priority_main",
    fractions=(0.8, 0.1, 0.1),
    random_state=42,
    labels=("CR", "NCR"),
    sources=("COD", "SI"),
    variants=None,
    metals=None,
    structure_ids=None,
    missing_parent="singleton",
    leakage_guard="auto",
    stratify_by=("label",),
)
```

### 12.2 One-call convenience function

```python
from CoREMOF.splitters import split_release

split = split_release(
    "/path/to/coremof_v26.0.2",
    checkers="5checker",
    parent_method="priority_main",
    leakage_guard="auto",
    fractions=(0.8, 0.1, 0.1),
    random_state=42,
    labels=("CR", "NCR"),
    sources=("COD", "SI"),
)
```

For a release path or unclassified `CoREMOFDataset`, omitted `checkers`
defaults to `5checker`. For an already classified object, omitted `checkers`
preserves its existing checker view. An explicitly conflicting checker request
raises rather than being silently ignored.

### 12.3 Split parameters

| Parameter | Default | Meaning |
|---|---|---|
| `parent_method` | `priority_main` | Explanatory parent relation |
| `fractions` | `(0.8, 0.1, 0.1)` | Train, validation, and test fractions |
| `random_state` | `42` | Deterministic integer or text tie-break seed |
| `labels` | `("CR", "NCR")` | Eligible checker labels; `None` means no label filter |
| `sources` | `None` | Eligible source databases |
| `variants` | `None` | Eligible ASR/FSR/ION variants |
| `metals` | `None` | At-least-one metal selection |
| `structure_ids` | `None` | Optional exact public-ID selection |
| `missing_parent` | `singleton` | `singleton` or `exclude` |
| `leakage_guard` | `auto` | `auto`, `main_union`, or `parent_only` |
| `stratify_by` | `("label",)` | Label and/or public metadata fields to balance approximately |
| `official` | `False` | Request an audited official manifest; currently fails closed |

For an explicit MOFDescribe-style splitter object:

```python
from CoREMOF.splitters import ParentGroupSplitter

splitter = ParentGroupSplitter(
    classified,
    parent_method="priority_main",
    leakage_guard="auto",
    labels=("CR", "NCR"),
    sources=("COD", "SI"),
    random_state=42,
)
split = splitter.train_valid_test_split((0.8, 0.1, 0.1))
```

`ParentSplitter` is an alias for `ParentGroupSplitter`, and
`classified.split(...)` is an alias for
`classified.train_valid_test_split(...)`.

Fractions must be an ordered `list`/`tuple` of exactly three finite,
non-boolean numeric values, be non-negative, sum to one within numerical
tolerance, and contain at least one positive value. Numeric text such as
`"0.8"` is rejected instead of being coerced. Zero validation or test
fractions are allowed deliberately:

```python
train_only = classified.train_valid_test_split(
    parent_method="priority_main",
    fractions=(1.0, 0.0, 0.0),
)
```

### 12.4 Stratification

The default splitter approximately balances checker labels while respecting
indivisible leakage blocks:

```python
stratify_by=("label",)
```

Additional public metadata fields can be included:

```python
split = classified.train_valid_test_split(
    stratify_by=("label", "source_database", "structure_variant"),
)
```

Every non-`label` key must exist in the metadata table. More or smaller strata
can make exact proportions harder to achieve. Leakage protection always takes
priority over ratio or stratum balance. Like the filters, `stratify_by` accepts
one exact string or an ordered `list`/`tuple` of exact strings; unordered and
one-shot containers are rejected so receipt hashes cannot depend on iteration
order.

## 13. Understanding `SplitResult`

### 13.1 IDs and indices

```python
print(split.train_ids[:5])
print(split.validation_ids[:5])
print(split.test_ids[:5])

train_idx, validation_idx, test_idx = split
```

The unpacked values are positions in the classified view supplied to the
splitter. The result exposes both index systems:

| Field | Index space |
|---|---|
| `train_indices`, `validation_indices`, `test_indices` | Classified view order |
| `train_release_indices`, `validation_release_indices`, `test_release_indices` | Original release metadata order |

`valid_ids` and `valid_indices` are concise aliases for the validation fields.

Persist structure IDs, not integer positions. Positions can change when the
view or metadata ordering changes.

### 13.2 Counts and balance

```python
print(dict(split.counts))
print(dict(split.achieved_fractions))
print({k: dict(v) for k, v in split.label_counts_by_split.items()})
```

`counts` contains `train`, `validation`, `test`, and `excluded`. Every release
row is either assigned or explicitly excluded.

Exact requested fractions may be impossible when a leakage block is large.
Inspect the achieved fractions rather than assuming exact 80/10/10 counts.

### 13.3 Leakage audit

```python
audit = split.leakage_audit
assert audit["passed"]
assert audit["cross_split_block_count"] == 0
print(audit["block_count"], audit["max_block_size"])
```

The audit reports:

- number of assigned leakage blocks;
- any block found in more than one partition;
- maximum assigned block size; and
- a final pass/fail flag.

### 13.4 Assignments, exclusions, and diagnostics

```python
structure_id = split.train_ids[0]
print(split.assignments[structure_id])
print(split.labels[structure_id])
print(split.parent_groups.get(structure_id))
print(split.leakage_groups[structure_id])

for structure_id, reason in list(split.exclusions.items())[:5]:
    print(structure_id, reason)
```

Common exclusion reasons include:

- `PRESELECTION_FILTER`;
- `STRUCTURE_ID_FILTER`;
- `LABEL_NOT_AVAILABLE`;
- `LABEL_FILTER`;
- `SOURCE_FILTER`;
- `VARIANT_FILTER`;
- `METAL_FILTER`;
- `MISSING_PARENT_EVIDENCE`;
- `PARENT_METHOD_CONFLICT`; and
- `MISSING_LEAKAGE_GROUP`.

With the default `priority_main`/`main_union` combination, an unresolved
hierarchy conflict can remain assigned safely inside the broader leakage block.
It then appears in `parent_diagnostics` instead of `exclusions`.

### 13.5 Warnings

```python
print(split.warnings)
```

Possible machine-readable warnings include:

- `PROVISIONAL_PARENT_INPUT`;
- `NO_ELIGIBLE_STRUCTURES`;
- `GROUP_CONSTRAINED_FRACTIONS`;
- `EMPTY_TRAIN_PARTITION`, `EMPTY_VALIDATION_PARTITION`, or
  `EMPTY_TEST_PARTITION`;
- `PARENT_METHOD_CONFLICTS_EXCLUDED`;
- `PARENT_METHOD_CONFLICTS_GUARDED`; and
- `PARENT_METHOD_CONFLICT_LEDGER_PRESENT`.

## 14. Writing CSV and JSON outputs

### 14.1 Write both files transactionally

```python
csv_path, receipt_path = split.write(
    "model_splits",
    stem="cod_si_5checker_seed42",
)
```

Existing files are not overwritten unless requested:

```python
split.write(
    "model_splits",
    stem="cod_si_5checker_seed42",
    overwrite=True,
)
```

The pair is staged before publication. If the second render fails, the method
restores or removes the first output instead of leaving a misleading half-pair.
The stem must be a simple filename stem and cannot contain directory traversal.

### 14.2 Write files separately

```python
split.to_csv("assignments.csv")
split.to_json("split_receipt.json")
```

### 14.3 Assignment CSV columns

The CSV contains one row for every structure in the release:

| Column | Meaning |
|---|---|
| `structure_id` | Stable public structure ID |
| `release_index` | Position in release metadata |
| `view_index` | Position in the classified view; blank outside that view |
| `split` | `train`, `validation`, `test`, or `EXCLUDED` |
| `label` | Recomputed checker-consensus label |
| `parent_group` | Selected explanatory parent group |
| `parent_diagnostic` | Row-level parent conflict diagnostic when present |
| `leakage_group` | Actual indivisible split block |
| `exclusion_reason` | Explicit reason for non-assignment |

A row removed before splitting through `classified.filter(...)` remains in the
CSV as `PRESELECTION_FILTER`; its `view_index` and label can be blank because
it was absent from the supplied classified view.

Example reading with the standard library:

```python
import csv

with open("model_splits/cod_si_5checker_seed42.csv", newline="") as handle:
    rows = list(csv.DictReader(handle))

train_ids = [row["structure_id"] for row in rows if row["split"] == "train"]
```

### 14.4 JSON receipt

The JSON receipt records:

- dataset and checker view;
- package/API/algorithm versions;
- implementation source hashes;
- release input hashes;
- requested filters, fractions, seed, and stratification;
- achieved counts and fractions;
- labels and partition IDs;
- parent method plus its machine-readable definition;
- requested and resolved leakage guard plus its machine-readable definition;
- parent conflicts and diagnostics;
- exclusions;
- assignment digest;
- CIF verification state;
- provisional/official state; and
- zero-cross-partition leakage audit.

Access the same object before writing:

```python
receipt = split.receipt()
print(receipt["assignment_sha256"])
print(receipt["checker_view_kind"])
print(receipt["parent_conflict_count"])
print(receipt["parent_method_definition"]["summary"])
print(receipt["requested_leakage_guard"], receipt["leakage_guard"])
```

Implementation hashes are frozen when `SplitResult` is constructed. Editing
source code later cannot silently change the receipt of an already-created
in-memory result.

## 15. Command-line usage

### 15.1 Help and installation check

```bash
coremof --version
coremof doctor
coremof split --help
```

### 15.2 Recommended COD+SI split

```bash
coremof split /path/to/coremof_v26.0.2 \
  --checkers 5checker \
  --parent-method priority_main \
  --leakage-guard auto \
  --labels CR NCR \
  --sources COD SI \
  --fractions 0.8 0.1 0.1 \
  --random-state 42 \
  --output-directory model_splits \
  --stem cod_si_5checker_seed42
```

The command prints a concise JSON summary to standard output and writes the
assignment CSV and receipt JSON.

### 15.3 Strict CIF verification from the CLI

```bash
coremof split /path/to/coremof_v26.0.2 \
  --verify-cifs \
  --output-directory verified_split
```

### 15.4 Other useful CLI filters

```bash
coremof split /path/to/coremof_v26.0.2 \
  --checkers 3checker \
  --parent-method rac5 \
  --leakage-guard parent_only \
  --labels CR NCR AMBIGUOUS \
  --sources COD \
  --variants ASR FSR \
  --metals Cu Zn \
  --stratify-by label structure_variant \
  --output-directory sensitivity_split
```

### 15.5 Merge and require targets from a configuration file

Keep multi-file schemas and scientific declarations in JSON. Relative paths
are resolved beside the configuration file:

```json
{
  "sources": [
    {
      "path": "results/xe_uptake.csv",
      "name": "xe_uptake_298K_1bar",
      "id_column": "structure_name",
      "target_columns": ["uptake"],
      "target_names": {"uptake": "xe_uptake"},
      "value_types": {"xe_uptake": "float"},
      "units": {"xe_uptake": "mol/kg"},
      "conditions": {
        "xe_uptake": {"temperature_K": 298, "pressure_bar": 1.0}
      }
    },
    {
      "path": "results/xe_kr_selectivity.jsonl",
      "name": "xe_kr_selectivity_298K_1bar",
      "id_column": "structure_id",
      "target_columns": ["selectivity"],
      "target_names": {"selectivity": "xe_kr_selectivity"},
      "units": {"xe_kr_selectivity": "dimensionless"},
      "conditions": {
        "xe_kr_selectivity": {
          "temperature_K": 298,
          "pressure_bar": 1.0,
          "feed": {"Xe": 0.2, "Kr": 0.8}
        }
      }
    }
  ],
  "alias_registry": {
    "path": "private/audited_structure_name_registry.csv",
    "current_id_column": "structure_id",
    "alias_columns": ["previous_public_id", "v11_structure_id"]
  },
  "feature_tables": ["rac5", "zeo", "topology"]
}
```

Create the joined modelling table:

```bash
coremof merge-targets /path/to/coremof_v26.0.2 \
  --config targets.json \
  --output-directory model_inputs \
  --stem xe_kr_targets
```

Or merge and split in one command:

```bash
coremof split /path/to/coremof_v26.0.2 \
  --target-config targets.json \
  --require-target xe_uptake \
  --require-target xe_kr_selectivity \
  --required-target-mode all \
  --parent-method priority_main \
  --output-directory model_splits \
  --stem xe_kr_split
```

CLI options are:

| Option | Default / purpose |
|---|---|
| positional `RELEASE` | Extracted release root |
| `--output-directory DIR` | Required output directory |
| `--checkers` | `5checker`; official presets only |
| `--parent-method` | `priority_main` |
| `--leakage-guard` | `auto` |
| `--missing-parent` | `singleton` |
| `--fractions TRAIN VALIDATION TEST` | `0.8 0.1 0.1` |
| `--random-state` | `42` |
| `--labels ...` | `CR NCR` |
| `--sources ...` | No source filter |
| `--variants ...` | No variant filter |
| `--metals ...` | No metal filter |
| `--structure-ids ...` | No exact-ID filter |
| `--target-config FILE` | No target join; JSON config when supplied |
| repeatable `--require-target NAME` | No target-availability filter |
| `--required-target-mode` | `all`; `any` is an explicit alternative |
| `--stratify-by ...` | `label` |
| `--stem` | `coremof_split` |
| `--overwrite` | Permit replacement of the output pair |
| `--verify-cifs` | Rehash all CIF bytes before splitting |
| `--official` | Require an official manifest; currently fails closed |

CLI validation, missing-path, integrity, and output-collision errors return
exit code 2 with a concise error message. Existing outputs require either a new
stem/directory or `--overwrite`.

Target-merge and screening bundles use race-safe create-if-absent publication
by default. Their explicit overwrite mode is sequential **single-writer**
replacement: concurrent processes writing the same directory and stem must be
serialized with an external lock or, preferably, use distinct immutable stems.

The Python API uses `labels=None` to include all labels. The CLI has no null
value; list all four explicitly:

```bash
--labels CR NCR AMBIGUOUS UNCHECKED
```

## 16. Common workflows

### 16.1 Export CR and NCR IDs without splitting

```python
from pathlib import Path
from CoREMOF.dataset import CoREMOFDataset

dataset = CoREMOFDataset.from_release("/path/to/coremof_v26.0.2")
classified = dataset.classify("5checker", sources=("COD",))

Path("cod_cr_ids.txt").write_text("\n".join(classified.cr_ids) + "\n")
Path("cod_ncr_ids.txt").write_text("\n".join(classified.ncr_ids) + "\n")
```

### 16.2 Compare the three official checker views

```python
for view in ("3checker", "4checker", "5checker"):
    classified = dataset.classify(view)
    print(view, dict(classified.label_counts()))
```

This compares definitions; it does not overwrite one view with another.

### 16.3 Compare main parent assumptions

```python
classified = dataset.classify("5checker")
results = {}

for method in ("rac5", "mofid_v2", "mofid_v1"):
    results[method] = classified.train_valid_test_split(
        parent_method=method,
        leakage_guard="parent_only",
        labels=("CR", "NCR"),
        sources=("COD", "SI"),
        random_state=42,
    )
    print(method, dict(results[method].counts))
```

These are sensitivity analyses under different scientific relations. They are
not interchangeable replicas.

### 16.4 Select a metal-containing subset

```python
cu_zn = dataset.classify(
    "5checker",
    labels=("CR", "NCR"),
    metals=("Cu", "Zn"),
)

split = cu_zn.train_valid_test_split(
    parent_method="priority_main",
    leakage_guard="auto",
)
```

The broader leakage graph is still built from the full release; filtering first
does not remove hidden bridges.

### 16.5 Split an exact list of structures

```python
requested = tuple(dataset.structure_ids[:3])

split = dataset.classify("5checker").train_valid_test_split(
    structure_ids=requested,
    labels=None,
    parent_method="priority_main",
)
```

Unrequested release rows remain in the assignment CSV as explicit exclusions.

### 16.6 Use pandas only in downstream analysis

Pandas is not required by the package, but it can read the output normally:

```python
import pandas as pd

assignments = pd.read_csv("model_splits/cod_si_5checker_seed42.csv")
train = assignments.loc[assignments["split"] == "train"]
```

## 17. Reproducibility and scientific interpretation

### 17.1 What is deterministic

For the same validated release, package source, parameters, and
`random_state`, the `structure_id -> partition` assignment and assignment
digest are deterministic. Python hash randomization and metadata row order do
not change that mapping.

Changing `random_state` can change which complete blocks enter each partition.
It cannot change labels, parent groups, or leakage-block membership.

### 17.2 Why requested fractions may not be exact

The splitter assigns whole leakage blocks. It never breaks a block merely to
hit an exact percentage. A large block can therefore make an exact 80/10/10
split impossible. This is reported through achieved fractions, maximum block
size, and `GROUP_CONSTRAINED_FRACTIONS`.

### 17.3 Use IDs rather than row numbers

Integer indices are convenience values tied to a particular release or view.
Persist:

- `structure_id`;
- the assignment CSV;
- the JSON receipt; and
- the exact release files used.

### 17.4 Current versus official splits

Current package outputs are reproducible exploratory splits. They contain:

```text
official_split = false
```

Passing `official=True` raises `OfficialSplitUnavailableError` because the v26
release currently contains no audited assignment manifest.

An official future v26.0.2 split must freeze v26.0.1 assignments. If a new
structure connects base groups already frozen in different partitions, the new
row must be reported as a bridge conflict rather than moving old rows silently.
That frozen-extension feature is not yet implemented.

### 17.5 Provisional parent inputs

`provisional_input` is derived only from a package-authenticated closed release
contract. Public or caller-supplied `release_status` text is descriptive and
cannot authenticate publication or clear the provisional flag. The current
closed RT/M2T contract is staged, non-decisive, and explicitly
`publication_authorized=false`; an arbitrary `FINAL` or `FINAL_CANDIDATE`
token is rejected rather than treated as authority.

Current v26 parent tables remain provisional pending the final approved MOFid
evidence and parent rebuild. The package is usable now for exploratory work,
but every MOFid-dependent or `priority_main` result must retain that caveat.

### 17.6 Recording and citing a split

For a reproducible publication or model release, report at least:

- CoRE-MOF dataset version and its official citation/DOI when available;
- CoREMOF-tools package version;
- checker view;
- parent method and resolved leakage guard;
- source/label/variant/metal filters;
- train/validation/test fractions and random seed; and
- the receipt's `assignment_sha256` and release input hashes.

Publish the assignment CSV and JSON receipt with the model when licensing
allows. Do not invent a citation for an unpublished local snapshot; cite the
official release record supplied by the data distributor.

## 18. Licensing and redistribution

Loading and splitting existing release tables do not require a CCDC licence.
The package does not invoke CSD software, MOSAEC, SETC-GAT, Zeo++, molSimplify,
MOFid, or CrystalNets during these operations.

However, selecting a structure does not grant redistribution rights:

- **COD:** default open-data publication base.
- **SI:** redistribution remains subject to an asset-level rights review.
- **CSD:** CIFs and structure-resolved CSD-derived outputs remain
  licence-gated/local-build unless CCDC grants written permission for the
  exact material.

The recommended COD or COD+SI split uses leakage blocks constructed from the
complete COD+CSD+SI release. A standalone open bundle cannot reconstruct a
CSD-mediated bridge unless it includes an audited opaque projection of the
private full-universe leakage groups. Do not rebuild the guard on COD-only rows
and call it equivalent.

Assignment CSVs produced from a full licensed release contain one row per
release structure, including excluded CSD rows. Review the output before
redistribution; source filtering is a modelling selection, not a licence
sanitizer.

## 19. Troubleshooting

| Problem | Likely cause | Action |
|---|---|---|
| `release directory does not exist` | Wrong path or archive not extracted | Pass the directory containing `dataset_info.json` |
| Missing metadata/parent file | Incomplete release tree | Re-extract or obtain the complete release |
| `ReleaseValidationError` | ID, label, group size/status, path, count, or manifest inconsistency | Treat the release as damaged/stale; do not bypass validation |
| Full SHA-256 required by `main_union` | Missing or truncated CIF manifest hash | Use the complete audited manifest or choose an explicit narrower experiment |
| CIF size/hash mismatch | Corrupted or changed CIF | Restore the exact release CIF; do not silently edit it |
| Unknown checker status | Misspelled or unsupported token | Use canonical public PASS/FAIL/NOT_AVAILABLE values |
| Unknown checker preset | Invalid official view | Use `3checker`, `4checker`, `5checker`, or an ordered canonical sequence |
| Unknown structure ID | Typo or ID from another release | Check `dataset.structure_ids` and release version |
| Unknown `stratify_by` field | Column absent from metadata | Use `label` or an existing metadata column |
| Fractions rejected | Wrong length, negative/NaN, or sum not equal to one | Supply train/validation/test fractions such as `(0.8, 0.1, 0.1)` |
| Output already exists | Safe overwrite protection | Prefer a new immutable stem/directory; overwrite is single-writer and same-stem concurrent writers require an external lock |
| Empty partition warning | Too few indivisible blocks for the requested nonzero partition | Use more data or set the fraction to zero deliberately |
| Requested ratios not reached | Large parent/leakage groups | Inspect `achieved_fractions` and `max_block_size`; do not split the group manually |
| `official=True` fails | No audited official assignment manifest exists | Generate an exploratory split and retain `official_split=false` |
| Preclassified checker conflict | Explicit `checkers` differs from the existing view | Omit `checkers` or pass the matching view |
| High memory use on full release | Complete metadata, conflict, and leakage ledgers are retained | Run one full split at a time and release unused objects between sensitivity runs |
| Unknown target structure ID | Input name is not a current ID and no alias maps it | Correct the name or supply an explicit audited alias registry |
| Ambiguous alias | One alias maps to multiple current IDs | Repair/audit the registry; the package will not choose one |
| Conflicting duplicate target | Two non-null observations disagree | Reconcile the source data or use scientifically distinct canonical target names |
| Conflicting target definition | One canonical name has different units/conditions/types | Rename the endpoints with `target_names` or correct the declaration |
| Missing required target | Target is null/absent for that row | Keep the explicit exclusion, choose `any`, or change the required endpoint list |

Expected validation errors are concise in the CLI. In Python, catch specific
exceptions only when you can handle them meaningfully:

```python
from CoREMOF.dataset import ReleaseValidationError
from CoREMOF.splitters import OfficialSplitUnavailableError

try:
    dataset = CoREMOFDataset.from_release("/path/to/release")
except (FileNotFoundError, ReleaseValidationError) as error:
    print("Cannot use this release:", error)
```

## 20. Current validation status and limitations

The package and this handbook are validated at release-candidate freeze with
the full focused package/notebook/target/screen/parent suite under both normal
Python and `python -S` without site-packages, dependency-complete discovery,
and opt-in real-release integration tests. Exact interpreter versions, test
counts, hashes, timings, and generated-output comparisons belong in the
immutable validation receipt produced for that freeze; this handbook avoids
embedding counts that become stale whenever a fail-closed regression is added.
The validation also includes:

- exact loading of 36,628 v26.0.1 and 42,574 v26.0.2 structures from the
  terminal staged pair;
- strict byte-level CIF-manifest verification;
- the full target/feature join, all notebook code cells, and deterministic
  target-aware split and screening examples against that same staged pair;
- byte-for-byte comparisons of scientific CSV outputs plus semantic
  comparisons of regenerated JSON receipts, with zero crossed leakage blocks;
- empty-environment wheel installation and `pip check`;
- fail-closed synthetic tests for ambiguous aliases, duplicate conflicts,
  canonical target renaming, target nulls, hidden leakage bridges, unrelated
  split bundles, forged screening rows, concurrent-publication rollback, and
  JSON object/array target canonicalization; and
- fail-closed StructureMatcher contract tests for method/version/semantics,
  optional-only hierarchy membership, parent-table and ledger hashes,
  receipt counts, public completeness policy, and relaxed-output exclusion.

For the local v26.0.2 provisional audited snapshot captured on 2026-08-03, a
five-checker COD+SI CR/NCR run has 5,902 eligible structures and produced
4,722/590/590 train/validation/test rows with zero crossed leakage blocks.
These counts are dated validation evidence, not permanent API constants.

As a full-release loading/classification sanity check for that same snapshot,
the five-checker counts are:

| Label | Count |
|---|---:|
| CR | 6,294 |
| NCR | 2,299 |
| AMBIGUOUS | 7,367 |
| UNCHECKED | 26,614 |

That same dated 2026-08-03 parent audit recorded 4,128
lower-versus-stronger group conflicts:

- 1,991 MOFid-v2 conflicts;
- 2,137 MOFid-v1 conflicts;
- 3,871 anchor-only conflicts; and
- 257 conflicts containing 654 unresolved rows.

They remain leakage-safe under the default full union and are now explicit in
the receipt.

At the current release-candidate freeze, strict documentation rendering and
the Python 3.9/3.10/3.11 package matrices are completed checks rather than
open gates. Remaining authorization or publication gates are:

- independent red-team verification of the frozen source, artifact, and
  generated-output hashes before any downstream terminal replay;
- a fresh terminal target/split/screen run after that independent GO;
- an explicitly audited assignment manifest before any split can be called
  official (the package intentionally emits exploratory splits otherwise);
- separate release-level publication authorization and publication-eligible
  MOFid evidence: the currently validated MOFid projection remains
  `STAGE_ONLY`, so it is a non-published candidate. The RT/M2T CrystalNets
  relations remain staged optional non-decisive references and are not
  consumed by `priority_main` or `main_union`; and
- asset-level redistribution clearance for any SI archives. The distributable
  wheel and source archive exclude those archives while that clearance is not
  established; the internal checkout remains unchanged.

## 21. Developer checks

From the source checkout:

```bash
cd /path/to/CoRE-MOF-Tools

# Lightweight focused suite.
python -m unittest \
  tests.test_handbook \
  tests.test_notebook \
  tests.test_dataset_labels \
  tests.test_parents_splitters \
  tests.test_targets \
  tests.test_cli \
  tests.test_screen_candidates_example

# Prove the new API does not depend on site-packages.
python -S -m unittest \
  tests.test_handbook \
  tests.test_notebook \
  tests.test_dataset_labels \
  tests.test_parents_splitters \
  tests.test_targets \
  tests.test_cli \
  tests.test_screen_candidates_example
```

Run opt-in real-release tests by supplying both roots:

```bash
COREMOF_V2601_RELEASE=/path/to/coremof_v26.0.1 \
COREMOF_V2602_RELEASE=/path/to/coremof_v26.0.2 \
python -m unittest tests.test_split_release_integration
```

Build a dependency-free wheel:

```bash
python -m pip wheel . --no-deps --no-build-isolation --wheel-dir dist
```

Before calling a package or split official, rerun all checks against the final
release inputs and verify the resulting assignment receipt independently.

## Further project documentation

- `docs/source/splitting.rst` — concise Sphinx user guide.
- `examples/CoREMOF_dataset_splitting_quickstart.ipynb` — executable end-to-end
  notebook.
- `README.md` — package overview and installation summary.

The loader does not download releases. External users should obtain an
authorized CoRE-MOF archive from the project/data distributor, preserve its
published checksum, extract it, and pass the extracted directory containing
`dataset_info.json`. CSD-containing material remains subject to the licence
rules in [Licensing and redistribution](#18-licensing-and-redistribution).
