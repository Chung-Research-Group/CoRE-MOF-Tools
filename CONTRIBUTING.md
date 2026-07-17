# Contributing

Bug fixes, documentation corrections, tests, and focused feature additions are welcome.

## Development setup

```bash
git clone https://github.com/Chung-Research-Group/CoRE-MOF-Tools.git
cd CoRE-MOF-Tools
conda env create -f env.yaml
conda activate coremof_tools
pip install -e '.[test,docs]'
```

Run the checks used by continuous integration:

```bash
ruff check CoREMOF tests setup.py --select F
python -m compileall -q CoREMOF tests
python -m unittest discover -s tests -v
python -m sphinx -W --keep-going -b html docs/source docs/build/html
```

## Pull requests

- Keep changes focused and explain the user-visible effect.
- Add a regression test for a bug fix when the affected component can run without licensed software.
- Preserve existing public function names and output keys unless the change is explicitly documented as breaking.
- Use temporary directories for intermediate files; do not introduce fixed filenames in the working directory.
- Do not commit `__pycache__`, notebook checkpoints, generated documentation, local environments, or proprietary CSD data.
- State which external programs and model files were used for scientific validation.

## Reporting scientific discrepancies

Include the exact function call, package version, dependency versions, complete traceback, and the smallest shareable CIF. For numerical disagreements, also include probe radii, sampling count, external-program version, and expected units. “Runs without error” is not sufficient validation of a scientific result.
