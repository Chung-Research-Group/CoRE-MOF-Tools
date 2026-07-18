"""Sphinx configuration for CoRE MOF Tools."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from CoREMOF import __version__

project = "CoRE MOF Tools"
copyright = "2026, MTAP at Pusan National University"
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
html_theme = "furo"
html_title = f"CoRE MOF Tools {release}"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_meta = {
    "description": (
        "Documentation for CoRE MOF Tools: database access, structure curation, "
        "validation, descriptors, and property prediction for metal-organic frameworks."
    ),
}
html_theme_options = {
    "sidebar_hide_name": False,
    "navigation_with_keys": True,
    "top_of_page_buttons": ["view", "edit"],
    "source_repository": "https://github.com/Chung-Research-Group/CoRE-MOF-Tools/",
    "source_branch": "main",
    "source_directory": "docs/source/",
    "light_css_variables": {
        "color-brand-primary": "#2177a9",
        "color-brand-content": "#176f9f",
        "color-link": "#176f9f",
        "color-link--hover": "#c77400",
        "color-sidebar-background": "#f5f8fa",
        "color-sidebar-background-border": "#dbe5ec",
    },
    "dark_css_variables": {
        "color-brand-primary": "#69bee9",
        "color-brand-content": "#69bee9",
        "color-link": "#69bee9",
        "color-link--hover": "#ffbd4a",
        "color-sidebar-background": "#132636",
        "color-sidebar-background-border": "#32495c",
    },
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/Chung-Research-Group/CoRE-MOF-Tools",
            "html": """
                <svg stroke="currentColor" fill="currentColor" stroke-width="0"
                     viewBox="0 0 16 16" aria-hidden="true"
                     height="1em" width="1em"
                     xmlns="http://www.w3.org/2000/svg">
                  <path d="M8 0C3.58 0 0 3.64 0 8.13c0 3.59 2.29 6.64 5.47 7.71.4.08.55-.18.55-.39 0-.19-.01-.83-.01-1.5-2.01.38-2.53-.5-2.69-.96-.09-.23-.48-.96-.82-1.15-.28-.15-.68-.53-.01-.54.63-.01 1.08.59 1.23.83.72 1.23 1.87.88 2.33.67.07-.53.28-.88.51-1.08-1.78-.21-3.64-.9-3.64-4.02 0-.89.31-1.62.82-2.19-.08-.2-.36-1.04.08-2.16 0 0 .67-.22 2.2.84A7.48 7.48 0 0 1 8 3.9c.68 0 1.36.09 2 .27 1.53-1.06 2.2-.84 2.2-.84.44 1.12.16 1.96.08 2.16.51.57.82 1.3.82 2.19 0 3.13-1.87 3.81-3.65 4.02.29.25.54.73.54 1.49 0 1.08-.01 1.95-.01 2.22 0 .21.15.47.55.39A8.15 8.15 0 0 0 16 8.13C16 3.64 12.42 0 8 0Z"></path>
                </svg>
            """,
            "class": "",
        },
    ],
}
pygments_style = "sphinx"
pygments_dark_style = "monokai"
