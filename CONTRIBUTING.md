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

## Project-defined terminology

Do not expose a coined identifier, policy name, abbreviation, or convenience
method without defining it. Before a project-defined term appears in a public
API, command-line option, notebook, metadata field, or receipt:

- explain it in plain language at first use and mark it as project-defined;
- provide one centralized machine-readable contract containing its purpose,
  exact ordered inputs and algorithm, conflict and missing-data rules, excluded
  inputs, output-ID grammar, and how it differs from adjacent concepts;
- make command-line help and user documentation point to the same semantics;
- record both the identifier and its expanded definition in generated receipts;
- for selectors such as `auto`, record the requested selector and the concrete
  resolved policy separately; and
- test every exposed choice, its definition lookup, and receipt serialization.

For example, `priority_main` is an explanatory parent hierarchy while
`main_union` is a broader leakage-block graph. A contribution must not present
either token without preserving that distinction.

## Reporting scientific discrepancies

Include the exact function call, package version, dependency versions, complete traceback, and the smallest shareable CIF. For numerical disagreements, also include probe radii, sampling count, external-program version, and expected units. “Runs without error” is not sufficient validation of a scientific result.
