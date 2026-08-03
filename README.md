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

### Build CR/NCR train/validation/test splits

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

The default `priority_main` policy is conflict-aware; it is not a row-by-row
first-nonmissing fallback. Its automatic leakage guard builds full-release
components from exact CIF, source-sibling, RAC5, MOFid-v2, and MOFid-v1 edges
before applying filters. Direct methods such as `rac5`, `mofid_v2`, and
`mofid_v1` are separately selectable for sensitivity studies. Missing parent
evidence becomes a unique singleton and never causes missing rows to match.
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
parent method and leakage guard, a structured lower-versus-stronger parent
conflict ledger, and a zero-cross-partition leakage audit.
`provisional_input` becomes false only when both the dataset and parent-method
release statuses are exactly `FINAL`; missing values, `FINAL_CANDIDATE`, and
all other tokens remain provisional. A user-generated split is never labelled
as an official CoRE-MOF split, and `official=True` currently fails closed
because no audited official assignment manifest exists.

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

Full guides and API documentation are available on [Read the Docs](https://core-mof-tools.readthedocs.io/). Executable notebooks and CIF examples are in [`examples/`](https://github.com/Chung-Research-Group/CoRE-MOF-Tools/tree/main/examples).

## Citation

If you use the database or toolkit, cite:

> Zhao G., Brabson L., Chheda S., Huang J., Kim H., Liu K., et al. “CoRE MOF DB: a curated experimental metal–organic framework database with machine-learned properties for integrated material-process screening.” *Matter* 8 (2025), 102140. [https://doi.org/10.1016/j.matt.2025.102140](https://doi.org/10.1016/j.matt.2025.102140)

Also cite the underlying method used in your workflow (for example Zeo++, PACMAN-charge, MOFClassifier, MOFChecker, MOSAEC, CrystalNets, or the relevant stability model). A method-by-method list is provided in the [documentation](https://core-mof-tools.readthedocs.io/).

## Support and development

- Report reproducible bugs through [GitHub Issues](https://github.com/Chung-Research-Group/CoRE-MOF-Tools/issues).
- Include `coremof doctor` output, Python version, operating system, a minimal code example, and—when shareable—the failing CIF.
- Developed by [Guobin Zhao](https://github.com/sxm13) at MTAP, Pusan National University.
