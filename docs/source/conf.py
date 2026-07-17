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
    "sphinx_design",
]

autodoc_mock_imports = [
    "PACMANCharge",
    "MOFClassifier",
    "ase",
    "ccdc",
    "cloudpickle",
    "gemmi",
    "joblib",
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
html_theme = "pydata_sphinx_theme"
html_title = f"CoRE MOF Tools {release}"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_meta = {
    "description": (
        "Documentation for CoRE MOF Tools: database access, structure curation, "
        "validation, descriptors, and property prediction for metal-organic frameworks."
    ),
}
html_context = {
    "github_user": "Chung-Research-Group",
    "github_repo": "CoRE-MOF-Tools",
    "github_version": "main",
    "doc_path": "docs/source",
}
html_theme_options = {
    "navbar_align": "left",
    "navigation_depth": 3,
    "show_toc_level": 2,
    "header_links_before_dropdown": 5,
    "show_nav_level": 2,
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/Chung-Research-Group/CoRE-MOF-Tools",
            "icon": "fa-brands fa-github",
        },
    ],
    "secondary_sidebar_items": ["page-toc", "sourcelink"],
    "footer_start": ["copyright"],
    "footer_center": ["sphinx-version"],
    "footer_end": ["theme-version"],
}
