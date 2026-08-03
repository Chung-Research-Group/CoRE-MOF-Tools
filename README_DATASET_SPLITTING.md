# CoREMOF-tools dataset classification and splitting handbook

This handbook explains how to use the lightweight CoREMOF-tools API to:

- load and validate a CoRE-MOF release;
- obtain CR, NCR, AMBIGUOUS, and UNCHECKED labels from a selected checker set;
- filter structures by source, variant, metal, label, or public ID;
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
11. [Creating train/validation/test splits](#11-creating-trainvalidationtest-splits)
12. [Understanding `SplitResult`](#12-understanding-splitresult)
13. [Writing CSV and JSON outputs](#13-writing-csv-and-json-outputs)
14. [Command-line usage](#14-command-line-usage)
15. [Common workflows](#15-common-workflows)
16. [Reproducibility and scientific interpretation](#16-reproducibility-and-scientific-interpretation)
17. [Licensing and redistribution](#17-licensing-and-redistribution)
18. [Troubleshooting](#18-troubleshooting)
19. [Current validation status and limitations](#19-current-validation-status-and-limitations)
20. [Developer checks](#20-developer-checks)

## 1. What the package does

The splitting package consumes an extracted CoRE-MOF release. It independently
recomputes checker-consensus labels from the public checker status columns and
uses the release's audited parent-group tables to prevent related structures
from being separated unintentionally.

It provides two related operations:

1. **Classification:** decide whether a structure is CR, NCR, AMBIGUOUS, or
   UNCHECKED under a selected checker view.
2. **Dataset splitting:** place eligible structures in train, validation, and
   test partitions while keeping a selected structural relation—or a broader
   leakage component—inside one partition.

The package does **not**:

- run MOFClassifier, MOFChecker, Chen-Manz, MOSAEC, or SETC-GAT;
- calculate RACs, MOFid, Zeo++, topology, or open-metal-site properties;
- repair or rewrite CIF files;
- impute missing checker or parent information;
- convert an execution error into a checker FAIL; or
- create an official CoRE-MOF benchmark assignment.

These boundaries are intentional. The package consumes validated scientific
evidence without silently changing the method that produced it.

Keep four choices conceptually separate:

| Choice | Question it answers |
|---|---|
| Checker view | How is CR/NCR/AMBIGUOUS/UNCHECKED recomputed? |
| Parent method | Which relation explains structural grouping? |
| Leakage guard | Which known relations must not cross partitions? |
| Filters | Which labelled rows are eligible for this experiment? |

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
is effectively required for the recommended `priority_main` split because the
default `main_union` leakage guard uses full CIF SHA-256 identity.

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

`StructureRecord.parent_group()` reads published table criteria such as
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

| Method | Role | Interpretation |
|---|---|---|
| `priority_main` | Recommended | Conflict-aware RAC5 → MOFid v2 → MOFid v1 hierarchy |
| `rac5` | Main/direct | Exact equality of the release-authorized depth-5 RAC fingerprint |
| `mofid_v2` | Main/direct | Exact normalized release-authorized MOFid v2 group when available |
| `mofid_v1` | Main/direct | Exact normalized MOFid v1 group when available |
| `rac5_zeo` | Reference | Exact combined RAC5 and selected Zeo++ fingerprint |
| `zeo` | Reference | Exact selected Zeo++ fingerprint |
| `source_id` | Reference | Normalized, database-namespaced source sibling |
| `common_name` | Reference | Normalized common-name match |
| `identity_union` | Reference | Audited identity-union screening relation |
| `none` | Control | Treat every structure as an independent singleton |

RAC/Zeo equality is fingerprint equivalence, not proof of a common synthetic
parent. `identity_union` is a screening relation, not proof that every member
is chemically identical.

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

## 11. Creating train/validation/test splits

### 11.1 Method on a classified dataset

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

### 11.2 One-call convenience function

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

### 11.3 Split parameters

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

Fractions must be finite, non-negative, contain exactly three values, sum to
one within numerical tolerance, and contain at least one positive value. Zero
validation or test fractions are allowed deliberately:

```python
train_only = classified.train_valid_test_split(
    parent_method="priority_main",
    fractions=(1.0, 0.0, 0.0),
)
```

### 11.4 Stratification

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
priority over ratio or stratum balance.

## 12. Understanding `SplitResult`

### 12.1 IDs and indices

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

### 12.2 Counts and balance

```python
print(dict(split.counts))
print(dict(split.achieved_fractions))
print({k: dict(v) for k, v in split.label_counts_by_split.items()})
```

`counts` contains `train`, `validation`, `test`, and `excluded`. Every release
row is either assigned or explicitly excluded.

Exact requested fractions may be impossible when a leakage block is large.
Inspect the achieved fractions rather than assuming exact 80/10/10 counts.

### 12.3 Leakage audit

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

### 12.4 Assignments, exclusions, and diagnostics

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

### 12.5 Warnings

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

## 13. Writing CSV and JSON outputs

### 13.1 Write both files transactionally

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

### 13.2 Write files separately

```python
split.to_csv("assignments.csv")
split.to_json("split_receipt.json")
```

### 13.3 Assignment CSV columns

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

### 13.4 JSON receipt

The JSON receipt records:

- dataset and checker view;
- package/API/algorithm versions;
- implementation source hashes;
- release input hashes;
- requested filters, fractions, seed, and stratification;
- achieved counts and fractions;
- labels and partition IDs;
- parent method, leakage guard, conflicts, and diagnostics;
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
```

Implementation hashes are frozen when `SplitResult` is constructed. Editing
source code later cannot silently change the receipt of an already-created
in-memory result.

## 14. Command-line usage

### 14.1 Help and installation check

```bash
coremof --version
coremof doctor
coremof split --help
```

### 14.2 Recommended COD+SI split

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

### 14.3 Strict CIF verification from the CLI

```bash
coremof split /path/to/coremof_v26.0.2 \
  --verify-cifs \
  --output-directory verified_split
```

### 14.4 Other useful CLI filters

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
| `--stratify-by ...` | `label` |
| `--stem` | `coremof_split` |
| `--overwrite` | Permit replacement of the output pair |
| `--verify-cifs` | Rehash all CIF bytes before splitting |
| `--official` | Require an official manifest; currently fails closed |

CLI validation, missing-path, integrity, and output-collision errors return
exit code 2 with a concise error message. Existing outputs require either a new
stem/directory or `--overwrite`.

The Python API uses `labels=None` to include all labels. The CLI has no null
value; list all four explicitly:

```bash
--labels CR NCR AMBIGUOUS UNCHECKED
```

## 15. Common workflows

### 15.1 Export CR and NCR IDs without splitting

```python
from pathlib import Path
from CoREMOF.dataset import CoREMOFDataset

dataset = CoREMOFDataset.from_release("/path/to/coremof_v26.0.2")
classified = dataset.classify("5checker", sources=("COD",))

Path("cod_cr_ids.txt").write_text("\n".join(classified.cr_ids) + "\n")
Path("cod_ncr_ids.txt").write_text("\n".join(classified.ncr_ids) + "\n")
```

### 15.2 Compare the three official checker views

```python
for view in ("3checker", "4checker", "5checker"):
    classified = dataset.classify(view)
    print(view, dict(classified.label_counts()))
```

This compares definitions; it does not overwrite one view with another.

### 15.3 Compare main parent assumptions

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

### 15.4 Select a metal-containing subset

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

### 15.5 Split an exact list of structures

```python
requested = tuple(dataset.structure_ids[:3])

split = dataset.classify("5checker").train_valid_test_split(
    structure_ids=requested,
    labels=None,
    parent_method="priority_main",
)
```

Unrequested release rows remain in the assignment CSV as explicit exclusions.

### 15.6 Use pandas only in downstream analysis

Pandas is not required by the package, but it can read the output normally:

```python
import pandas as pd

assignments = pd.read_csv("model_splits/cod_si_5checker_seed42.csv")
train = assignments.loc[assignments["split"] == "train"]
```

## 16. Reproducibility and scientific interpretation

### 16.1 What is deterministic

For the same validated release, package source, parameters, and
`random_state`, the `structure_id -> partition` assignment and assignment
digest are deterministic. Python hash randomization and metadata row order do
not change that mapping.

Changing `random_state` can change which complete blocks enter each partition.
It cannot change labels, parent groups, or leakage-block membership.

### 16.2 Why requested fractions may not be exact

The splitter assigns whole leakage blocks. It never breaks a block merely to
hit an exact percentage. A large block can therefore make an exact 80/10/10
split impossible. This is reported through achieved fractions, maximum block
size, and `GROUP_CONSTRAINED_FRACTIONS`.

### 16.3 Use IDs rather than row numbers

Integer indices are convenience values tied to a particular release or view.
Persist:

- `structure_id`;
- the assignment CSV;
- the JSON receipt; and
- the exact release files used.

### 16.4 Current versus official splits

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

### 16.5 Provisional parent inputs

`provisional_input` becomes false only when both the dataset and parent-method
artifacts declare the exact release status `FINAL`. Values such as
`FINAL_CANDIDATE`, missing status, or any other token remain provisional.

Current v26 parent tables remain provisional pending the final approved MOFid
evidence and parent rebuild. The package is usable now for exploratory work,
but every MOFid-dependent or `priority_main` result must retain that caveat.

### 16.6 Recording and citing a split

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

## 17. Licensing and redistribution

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

## 18. Troubleshooting

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
| Output already exists | Safe overwrite protection | Choose a new stem/directory or set `overwrite=True` / `--overwrite` |
| Empty partition warning | Too few indivisible blocks for the requested nonzero partition | Use more data or set the fraction to zero deliberately |
| Requested ratios not reached | Large parent/leakage groups | Inspect `achieved_fractions` and `max_block_size`; do not split the group manually |
| `official=True` fails | No audited official assignment manifest exists | Generate an exploratory split and retain `official_split=false` |
| Preclassified checker conflict | Explicit `checkers` differs from the existing view | Omit `checkers` or pass the matching view |
| High memory use on full release | Complete metadata, conflict, and leakage ledgers are retained | Run one full split at a time and release unused objects between sensitivity runs |

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

## 19. Current validation status and limitations

The package and this handbook have been validated with:

- 61 focused tests under normal Python;
- the same 61 focused tests under `python -S` without site-packages;
- 75 passing dependency-complete tests, with the three opt-in real-release
  integrations also executed separately and passing;
- exact loading of 36,628 v26.0.1 and 42,574 v26.0.2 structures;
- strict byte-level rehashing of all 42,574 v26.0.2 CIFs;
- empty-environment wheel installation and `pip check`; and
- an independent implementation review with no release-blocking package defect.

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

The current parent audit records 4,128 lower-versus-stronger group conflicts:

- 1,991 MOFid-v2 conflicts;
- 2,137 MOFid-v1 conflicts;
- 3,871 anchor-only conflicts; and
- 257 conflicts containing 654 unresolved rows.

They remain leakage-safe under the default full union and are now explicit in
the receipt.

Remaining publication gates are:

- final approved MOFid evidence and rebuilt parent tables;
- rerunning the release integrations after that rebuild;
- audited official/frozen-base assignment-manifest support;
- Sphinx documentation rendering; and
- Python 3.10 and 3.11 package test matrices.

## 20. Developer checks

From the source checkout:

```bash
cd /path/to/CoRE-MOF-Tools

# Lightweight focused suite.
python -m unittest \
  tests.test_handbook \
  tests.test_notebook \
  tests.test_dataset_labels \
  tests.test_parents_splitters \
  tests.test_cli

# Prove the new API does not depend on site-packages.
python -S -m unittest \
  tests.test_handbook \
  tests.test_notebook \
  tests.test_dataset_labels \
  tests.test_parents_splitters \
  tests.test_cli
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
rules in [Licensing and redistribution](#17-licensing-and-redistribution).
