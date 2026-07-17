"""Sphinx configuration for CoRE MOF Tools."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from CoREMOF import __version__

project = "CoRE MOF Tools"
copyright = "2025–2026, MTAP at Pusan National University"
author = "Guobin Zhao and contributors"
version = __version__
release = __version__

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx_autodoc_typehints",
]

autodoc_mock_imports = [
    "PACMANCharge",
    "MOFClassifier",
    "ase",
    "ccdc",
    "cloudpickle",
    "gemmi",
    "juliacall",
    "keras",
    "matminer",
    "matplotlib",
    "mendeleev",
    "mofchecker",
    "mofid",
    "molSimplify",
    "networkx",
    "numpy",
    "openbabel",
    "pandas",
    "phonopy",
    "pymatgen",
    "scipy",
    "selfies",
    "sklearn",
    "tensorflow",
    "torch",
    "xgboost",
    "yaml",
]
autodoc_member_order = "bysource"
autoclass_content = "both"
napoleon_google_docstring = True
napoleon_numpy_docstring = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}

templates_path = ["_templates"]
exclude_patterns = []
html_theme = "sphinx_rtd_theme"
html_title = f"CoRE MOF Tools {release}"
