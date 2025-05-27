# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os, sys
sys.path.insert(0, os.path.abspath('../..')) 

project = 'CoREMOF'
copyright = '2025, MTAP @ Pusan National University'
author = 'Guobin Zhao'
release = '0.1.8'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']


extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx_autodoc_typehints'
]
autodoc_mock_imports = [
    "tensorflow", "torch", "zeopp", "pymatgen", "ase", "molSimplify","optree",
    "PACMAN_charge", "mofchecker", "gemmi", "phonopy", "xgboost", 'juliacall', 'mofid','keras',
    'cloudpickle','scikit-learn==1.3.2','networkx', 'selfies', 'mendeleev', 'requests'
]

