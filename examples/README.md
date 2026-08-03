# Examples

Run `coremof doctor` before using an example. The notebooks are grouped by the external requirements they need.

For the standard-library-only CR/NCR classification and parent-aware dataset
split workflow, start with
[`CoREMOF_dataset_splitting_quickstart.ipynb`](CoREMOF_dataset_splitting_quickstart.ipynb).

| Directory | Purpose | Main requirements |
|---|---|---|
| `checker/` | Structural validation examples | MOFChecker, MOFClassifier; CSD API for MOSAEC |
| `curation/` | Database download, lookup, and CIF curation | ASE, pymatgen, gemmi; CSD API for CSD downloads |
| `features/` | Geometric, topology, OMS, RAC, and heat-capacity features | Zeo++, juliacall/CrystalNets, molSimplify as applicable |
| `ion_pacman/` | PACMAN charge-based ion-preserving curation | PACMAN-charge |

Examples write outputs relative to the current working directory. Copy an example to a separate working directory if you want to preserve the checked-in reference outputs.
