<img src="https://raw.githubusercontent.com/Chung-Research-Group/CoRE-MOF-Tools/main/logo.png" alt="CoRE MOF Tools logo" width="500"/>

# CoRE MOF Tools

[![Documentation](https://img.shields.io/badge/docs-Read_the_Docs-blue?logo=readthedocs)](https://coremof-tools.readthedocs.io/)
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

Python 3.9–3.11 is supported. A clean conda environment is recommended because several scientific dependencies contain compiled extensions.

```bash
conda create -n coremof python=3.11
conda activate coremof
pip install CoREMOF-tools
coremof doctor
```

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

Full guides and API documentation are available on [Read the Docs](https://coremof-tools.readthedocs.io/). Executable notebooks and CIF examples are in [`examples/`](https://github.com/Chung-Research-Group/CoRE-MOF-Tools/tree/main/examples).

## Citation

If you use the database or toolkit, cite:

> Zhao G., Brabson L., Chheda S., Huang J., Kim H., Liu K., et al. “CoRE MOF DB: a curated experimental metal–organic framework database with machine-learned properties for integrated material-process screening.” *Matter* 8 (2025), 102140. [https://doi.org/10.1016/j.matt.2025.102140](https://doi.org/10.1016/j.matt.2025.102140)

Also cite the underlying method used in your workflow (for example Zeo++, PACMAN-charge, MOFClassifier, MOFChecker, MOSAEC, CrystalNets, or the relevant stability model). A method-by-method list is provided in the [documentation](https://coremof-tools.readthedocs.io/).

## Support and development

- Report reproducible bugs through [GitHub Issues](https://github.com/Chung-Research-Group/CoRE-MOF-Tools/issues).
- Include `coremof doctor` output, Python version, operating system, a minimal code example, and—when shareable—the failing CIF.
- Developed by [Guobin Zhao](https://github.com/sxm13) at MTAP, Pusan National University.
