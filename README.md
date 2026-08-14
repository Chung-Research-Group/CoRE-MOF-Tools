<img src="https://raw.githubusercontent.com/Chung-Research-Group/CoRE-MOF-Tools/main/logo.png" alt="CoRE MOF Tools logo" width="500"/>

# CoRE MOF Tools

[![Documentation](https://img.shields.io/badge/docs-Read_the_Docs-blue?logo=readthedocs)](https://core-mof-tools.readthedocs.io/)
[![PyPI](https://img.shields.io/pypi/v/CoREMOF-tools?logo=pypi)](https://pypi.org/project/CoREMOF-tools/)
[![Python](https://img.shields.io/badge/Python-3.9--3.11-blue.svg?logo=python)](https://python.org/downloads/)
[![License](https://img.shields.io/github/license/Chung-Research-Group/CoRE-MOF-Tools)](https://github.com/Chung-Research-Group/CoRE-MOF-Tools/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15055758.svg)](https://doi.org/10.5281/zenodo.15055758)

CoRE MOF Tools is the Python interface accompanying the CoRE MOF database. It provides database lookup and structure download, CIF curation and validation, geometric and chemical descriptors, and pretrained property models.

The package combines pure-Python functions with optional external software. Run the installation check before starting:

```bash
coremof doctor
# or: python -m CoREMOF doctor
```

## Installation

Python 3.9–3.11 is supported. The release-loading, classification, and
dataset-splitting API uses only the Python standard library. Version
`0.4.0.dev0` is not yet the stable PyPI release; install this checkout to use
the new API:

```bash
conda create -n coremof python=3.11
conda activate coremof
git clone https://github.com/Chung-Research-Group/CoRE-MOF-Tools.git
cd CoRE-MOF-Tools
python -m pip install .
coremof doctor
```

Install the historical scientific feature set with the `full` extra. A clean
conda environment is recommended because several of these dependencies contain
compiled extensions:

```bash
python -m pip install ".[full]"
coremof doctor
```

Version 0.4 changes the installation contract: bare installation is the
lightweight dataset/splitting API, while `[full]` preserves the dependency set
previously installed by default.

To reproduce the repository environment, use:

```bash
conda env create -f env.yaml
conda activate coremof_tools
```

Some features require additional software:

| Feature | Additional requirement |
|---|---|
| Zeo++ pore geometry | `conda install -c conda-forge zeopp-lsmo` |
| CSD structure download and MOSAEC | Licensed CSD software and CSD Python API |
| MOFChecker | `pip install git+https://github.com/Au-4/mofchecker_2.0.git@main` |
| MOFid v1/v2 | [MOFid installation](https://snurr-group.github.io/mofid/compiling/#installation) and Open Babel |
| Crystal topology | Julia/CrystalNets through `juliacall`; the first call may install Julia packages |
| Heat-capacity prediction | Full repository checkout containing the ~1.3 GB ensemble model directory |

## Quick start

For a complete practical guide to release loading, CR/NCR definitions,
parent-group choices, leakage guards, Python/CLI examples, output schemas,
licensing, and troubleshooting, see
[the dataset-splitting handbook](https://github.com/Chung-Research-Group/CoRE-MOF-Tools/blob/main/README_DATASET_SPLITTING.md).
For a complete executable walkthrough, open
[the companion notebook](https://github.com/Chung-Research-Group/CoRE-MOF-Tools/blob/main/examples/CoREMOF_dataset_splitting_quickstart.ipynb).

Current v26 status: the null-unresolved MOFid projection is explicitly
`STAGE_ONLY`. Parent relations built from it are non-published candidates and
cannot be promoted by the publication command; splits from current live or
staged inputs remain exploratory rather than official benchmark assignments.

### Build CR/NCR train/validation/test splits

Four project-defined API identifiers are used below; none is a standard
crystallographic term. A **release-authorized parent triad** means one
``status/group/size`` triple declared by
`parent_groups/parent_group_methods.json`, stored in
`parent_groups/parent_groups.csv`, and validated by the loader. `MATCHED`
means the criterion is available and the observed group has at least two
members; `UNMATCHED` means it is available and has one member;
`NOT_AVAILABLE` means it supplies no scientific edge and its size-one group is
only a table-shape placeholder. Two unavailable rows never match.

Here **exact RAC5** means equality of all 264 ordered finite depth-5 values
after binary64 parsing, conversion of `-0.0` to `+0.0`, and `float.hex()`
serialization (`rtol=atol=0`), with no scaling, deletion, imputation, or
rounding. **Exact MOFid** means equality of the complete release-authorized
current string after converting it to text, collapsing each Unicode-whitespace
run to one ASCII space, trimming, rejecting an empty or declared whole-field
missing/execution placeholder, applying Unicode NFKC, and then case-folding;
it is never prefix, substring, or fuzzy matching and never edits a CIF or its
chemistry. Placeholder comparison is case-insensitive and rejects `-`, `nan`,
`none`, `null`, `n/a`, `na`, `unknown`, `missing`, `timeout`, `timed out`,
`error`, `failed`, `fail`, `fail process`, `failed process`, and
`process failed`.

- `priority_main` is the explanatory parent hierarchy. Across the complete
  release it reads `rac_status/rac_group/rac_size`,
  `mofid2_status/mofid2_group/mofid2_size`, and
  `mofid1_status/mofid1_group/mofid1_size` from
  `parent_groups/parent_groups.csv`. It
  creates exact RAC5 anchor components, then uses exact MOFid-v2 groups and
  finally exact MOFid-v1 groups. At each lower step, a group touching
  no stronger component creates a component; one touching exactly one stronger
  component attaches only its unresolved rows; one touching two or more never
  merges those stronger components, records `PARENT_METHOD_CONFLICT`, and
  leaves lower-only rows unresolved. Missing evidence becomes one unique
  singleton per structure unless exclusion is explicitly requested. This is
  not a row-wise first-nonmissing fallback. Here priority means parent-evidence
  precedence, not a queue: it does not rank, schedule, or recalculate failed
  scientific features. Zeo++, topology, source IDs, CIF
  hashes, common names, provisional source-ID/MOFid transitive groups, and
  StructureMatcher do not enter this hierarchy.
- `main_union` is the separate conservative leakage guard, not a parent claim.
  It reads the full CIF SHA-256 values from the `sha256` column in
  `manifests/cif_manifest.csv` and the database-namespaced
  source-sibling/RAC5/MOFid group/status/size columns in `parent_groups.csv`.
  It forms transitive
  connected components over the full unfiltered release from those five exact
  relations. The release source group is exact equality of the ordered
  `(source_database, source_id)` pair after applying the text procedure above
  to each field separately, so IDs never match across databases. A missing CIF SHA-256 fails closed; a missing optional relation
  adds no edge, and nulls never match. There is no evidence precedence or
  conflict resolution: every available listed edge is unioned.
- `leakage_guard="auto"` is only a selector: it chooses `main_union` for
  `priority_main` and `parent_only` for an explicitly chosen direct/reference
  method. `parent_only` uses only that selected explanatory grouping as split
  blocks; it adds none of `main_union`'s cross-method edges. With
  `priority_main`, an unresolved lower-method conflict is excluded under
  `parent_only` but can remain safely assigned and diagnosed under
  `main_union`.

Load an extracted release, recompute a checker view, and keep related
structures in the same partition:

```python
from CoREMOF.dataset import CoREMOFDataset

dataset = CoREMOFDataset.from_release("/path/to/coremof_v26.0.2")
classified = dataset.classify(checkers="5checker")
print(len(classified.cr_ids), len(classified.ncr_ids))
split = classified.train_valid_test_split(
    parent_method="priority_main",  # RAC5, then MOFid v2, then MOFid v1
    leakage_guard="auto",
    labels=("CR", "NCR"),
    sources=("COD", "SI"),
    fractions=(0.8, 0.1, 0.1),
    random_state=42,
)

print(split.train_ids[:3])
split.write("my_split")
```

The default `priority_main` policy is not a row-by-row first-nonmissing
fallback, and its explanatory groups are distinct from the `main_union`
partition blocks. Both are computed before applying experiment filters.
Direct methods such as `rac5`, `mofid_v2`, and
`mofid_v1` are separately selectable for sensitivity studies. Missing parent
evidence becomes a unique singleton and never causes missing rows to match.
Releases may also expose three optional non-decisive reference criteria.
`rac5_crystalnets` (`RT-`) reads the public
`rac_crystalnets_status/group/size` triad and means exact equality of all 264
finite RAC5 values plus a complete
successful **current CrystalNets fingerprint**, where current means the
topology evidence authorized by the loaded release rather than a runtime
search for newer output. The fingerprint requires `SUCCESS`,
`topology_available=true`, `error=null`, nonempty complete SingleNodes and
AllNodes subnets, and count/catenation values equal to the subnet count. It
contains network/count/net/agreement fields and every subnet node's status,
dimension, key, name, and genome in sorted-key JSON; subnet projections are
sorted, duplicate subnets are retained, and SHA-256 is applied. Null top-level
dimension/net/agreement summaries are retained for heterogeneous subnets, and
node topology name/genome may be null; runtime, paths, hashes, diagnostics,
software text, and original subnet order are excluded.
`mofid_v2_crystalnets` (`M2T-`) reads the public
`mofid2_crystalnets_status/group/size` triad and replaces RAC5 with exact
complete canonicalized MOFid-v2 text: convert to text, collapse each
Unicode-whitespace run to one
ASCII space, trim, reject empty/whole-field missing or execution placeholders,
apply Unicode NFKC, and then case-fold. This text processing does not modify a
CIF, atoms, bonds, occupancies, coordinates, chemistry, or unit cell and is
not fuzzy matching. The eligible MOFid-v2 statuses are exactly `SUCCESS`,
`SUCCESS_TOPOLOGY_UNKNOWN`, `SUCCESS_TOPOLOGY_ERROR`, and
`SUCCESS_TOPOLOGY_TIMEOUT`. The latter two are successful calculated
identifiers whose embedded topology qualifier is ERROR or TIMEOUT, not MOFid
execution failures; every other MOFid-v2 status and every incomplete
CrystalNets input adds no edge. This reference method remains provisional
whenever the loaded release's MOFid-v2 input is provisional. If the
release-authorized MOFid-v2 values change, rebuild the M2T groups before use.
`structure_matcher_strict` (`SM-`) is a convenience connected component of
exhaustive composition-compatible pairs whose symmetric forward and reverse
pymatgen 2024.2.8
`ElementComparator` fits both pass with `ltol=stol=0.001`,
`angle_tol=0.01`, `primitive_cell/attempt_supercell=true`,
`scale/allow_subset=false`, `supercell_size=num_sites`, and no ignored species;
direct edges are authoritative; the component is not proof that every pair in
it directly matches or that its members are duplicates. Parser failures,
timeouts, and execution errors are `NOT_AVAILABLE` rather than unmatched and
add no optional edge. Missing or failed input adds no optional edge. The digest after each prefix is only a
group label, not a topology, MOFid, RMSD, or score. These methods do not change
`priority_main` or `main_union`: the recommended explanatory hierarchy remains
exact RAC5, then complete MOFid v2, then complete MOFid v1. The loader accepts `sm_*`
columns only with the exact optional-reference method declaration and its
hash-verified release-adapter receipt; incomplete components must project to
unique `NOT_AVAILABLE` singletons, and relaxed evidence cannot be executed or exposed.
The strict parser expands declared symmetry, merges generated sites at
`site_tolerance=1e-4`, rounds fractions near 1/3 or 2/3 at
`frac_tolerance=1e-4`, checks occupancy, sorts the periodic structure, and
preserves disorder; parser, timeout, OOM, matcher, or directional-disagreement
cases are unavailable rather than nonmatches. Its directional normalized RMS
and maximum displacements divide periodic site displacement by
`(V/Nsites)^(1/3)` for that direction, so they are dimensionless and are not
angstrom RMSD or output from the separate `charnley/rmsd` package.
The COD/SI filter still constructs its leakage blocks over the complete
COD+CSD+SI release. A future standalone open-data bundle must therefore ship
an audited full-universe block projection before it can reproduce
CSD-mediated bridges without exposing CSD rows.

The same workflow is available from the command line:

```bash
coremof split /path/to/coremof_v26.0.2 \
  --checkers 5checker \
  --sources COD SI \
  --output-directory my_split
```

Every output includes explicit exclusions, release-input hashes, the selected
parent method, the requested and resolved leakage guard, machine-readable
definitions of both project-defined policies, a structured
lower-versus-stronger parent conflict ledger, and a zero-cross-partition
leakage audit.
`provisional_input` is derived from a package-authenticated closed release
contract, never from caller-supplied status strings. The current authenticated
RT/M2T reference contract is staged, non-decisive, and explicitly not
publication-authorized, so current outputs remain provisional; an arbitrary
`FINAL` token is rejected and cannot clear that flag. A user-generated split
is never labelled as an official CoRE-MOF split, and `official=True` currently
fails closed because no audited official assignment manifest exists.

### Join target results before splitting

Combine one or more uptake, selectivity, or other endpoint files with current
release metadata and selected feature tables before assigning partitions:

```python
from CoREMOF.targets import TargetSource

targets = dataset.merge_targets(
    (
        TargetSource(
            "uptake.csv",
            id_column="structure_name",
            target_columns=("uptake",),
            target_names={"uptake": "xe_uptake"},
            value_types={"xe_uptake": "float"},
            units={"xe_uptake": "mol/kg"},
            conditions={"xe_uptake": {"temperature_K": 298, "pressure_bar": 1}},
        ),
    ),
    feature_tables=("rac5", "zeo"),
)
split = targets.classify("5checker").train_valid_test_split(
    required_targets=("xe_uptake",),
    parent_method="priority_main",
)
```

Matching uses exact current IDs. Earlier IDs require an explicit audited alias
registry; fuzzy matching, unit inference, target imputation, and silent
conflict resolution are never performed. See the handbook for multi-file,
alias-registry, provenance, configuration-file, and CLI examples.

### Screen candidates automatically

[`examples/screen_candidates.py`](examples/screen_candidates.py) turns the
same validated release or target configuration into a deterministic ranked
CSV and hash-bound receipt:

```bash
python examples/screen_candidates.py /path/to/coremof_v26.0.2 \
  --target-config targets.json \
  --rank-by xe_uptake \
  --require-target xe_uptake \
  --checkers 5checker \
  --label CR \
  --source COD \
  --order descending \
  --limit 1000 \
  --split \
  --parent-method priority_main \
  --leakage-guard auto \
  --output-directory screening/xe_uptake
```

Eligibility and required-target filters run before ranking. Null,
non-numeric, and non-finite values are excluded without imputation; ties use
ascending `structure_id`. Integers and decimal text are ranked without
binary64 rounding, while native floats retain their IEEE-754 values. The
optional split still builds leakage components from the full release and is
exploratory, not an official CoRE-MOF split.
See the [executable notebook](examples/CoREMOF_dataset_splitting_quickstart.ipynb)
for portable and guarded real-target examples.

### Query a CoRE MOF record

```python
from CoREMOF.structure import information

record = information("CR-ASR", "2020[Cu][sql]2[ASR]1")
print(record)
```

Valid datasets are `CR-ASR`, `CR-FSR`, `CR-Ion`, and `NCR`. Invalid dataset or entry names now raise an error with valid choices or close matches.

### Calculate pore geometry with Zeo++

```python
from CoREMOF.calculation import Zeopp

diameters = Zeopp.PoreDiameter("my_mof.cif")
print(diameters["LCD"], diameters["PLD"], diameters["unit"])

surface_area = Zeopp.SurfaceArea(
    "my_mof.cif",
    chan_radius=1.86,
    probe_radius=1.86,
    num_samples=10_000,
)
print(surface_area["ASA"])
```

Zeo++ calls use unique temporary files and are safe for paths containing spaces and for concurrent workflows. Set `COREMOF_NETWORK_EXECUTABLE` if the executable is not named `network`.

### Precheck and standardize a CIF

```python
from CoREMOF.curate import preprocess

job = preprocess("my_mof.cif", output_folder="result_curation")
print(job.result_check)
```

The curation classes write CIF/JSON/CSV outputs to the requested output directory. They do not return the processed CIF in memory.

### Predict properties

```python
from CoREMOF.prediction import pacman, stability

charge_result = pacman("my_mof.cif", output_folder="result_pacman")
stability_result = stability("my_mof.cif")
```

`stability()` additionally calls Zeo++ and RAC descriptors. Check all requirements with `coremof doctor` before running it.

## Feature map

| Task | Module | Main entry point |
|---|---|---|
| Database lookup | `CoREMOF.structure` | `information()` |
| SI/CSD structure download | `CoREMOF.structure` | `download_from_SI`, `download_from_CSD()` |
| Zeo++ geometry | `CoREMOF.calculation.Zeopp` | `PoreDiameter()`, `SurfaceArea()`, `PoreVolume()` |
| Basic descriptors | `CoREMOF.calculation.mof_features` | `SpaceGroup()`, `Mass()`, `Volume()`, `n_atom()` |
| Topology, OMS, RACs | `CoREMOF.calculation.mof_features` | `topology()`, `get_oms_file()`, `RACs()` |
| CIF curation | `CoREMOF.curate` | `preprocess`, `clean`, `clean_pacman` |
| Structural validation | `CoREMOF.curate` | `mof_check`, `run_MOSAEC()`, `run_mofclassifier()` |
| Property models | `CoREMOF.prediction` | `pacman()`, `cp()`, `stability()` |
| MOF identifiers | `CoREMOF.get_mofid` | `run_v1()`, `run_v2()` |
| Dataset classification and splitting | `CoREMOF.dataset`, `CoREMOF.splitters` | `CoREMOFDataset`, `split_release()` |
| Feature/target joining | `CoREMOF.targets` | `TargetSource`, `AliasRegistry`, `merge_targets()` |

Full guides and API documentation are available on [Read the Docs](https://core-mof-tools.readthedocs.io/). Executable notebooks and CIF examples are in [`examples/`](https://github.com/Chung-Research-Group/CoRE-MOF-Tools/tree/main/examples).

## Citation

If you use the database or toolkit, cite:

> Zhao G., Brabson L., Chheda S., Huang J., Kim H., Liu K., et al. “CoRE MOF DB: a curated experimental metal–organic framework database with machine-learned properties for integrated material-process screening.” *Matter* 8 (2025), 102140. [https://doi.org/10.1016/j.matt.2025.102140](https://doi.org/10.1016/j.matt.2025.102140)

Also cite the underlying method used in your workflow (for example Zeo++, PACMAN-charge, MOFClassifier, MOFChecker, MOSAEC, CrystalNets, or the relevant stability model). A method-by-method list is provided in the [documentation](https://core-mof-tools.readthedocs.io/).

## Support and development

- Report reproducible bugs through [GitHub Issues](https://github.com/Chung-Research-Group/CoRE-MOF-Tools/issues).
- Include `coremof doctor` output, Python version, operating system, a minimal code example, and—when shareable—the failing CIF.
- Developed by [Guobin Zhao](https://github.com/sxm13) at MTAP, Pusan National University.
