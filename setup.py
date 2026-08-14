from setuptools import setup, find_packages
from pathlib import Path
import re

here = Path(__file__).parent
readme_path = here / "README.md"
long_description = readme_path.read_text(encoding="utf-8") if readme_path.exists() else ""
init_text = (here / "CoREMOF" / "__init__.py").read_text(encoding="utf-8")
version = re.search(r'^__version__ = "([^\"]+)"', init_text, re.MULTILINE).group(1)

LEGACY_SCIENTIFIC_REQUIREMENTS = [
    'pymatgen',
    'ase',
    'juliacall',
    'molSimplify',
    'PACMAN-charge',
    'cloudpickle',
    'joblib',
    'matminer',
    'numpy',
    'pandas',
    'PyYAML',
    'scipy',
    'xgboost',
    'scikit-learn==1.5.0',
    'gemmi==0.7.0',
    'phonopy',
    'networkx',
    'selfies',
    'mendeleev',
    'requests',
    'MOFClassifier==0.1.1',
]

setup(
    name='CoREMOF_tools',
    version=version,
    author='Guobin Zhao',
    author_email='sxmzhaogb@gmail.com',
    description='Python API for CoRE MOF DB',
    license="CC-BY-4.0",
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'CoREMOF': [
            'calculation/juliapkg.json',
            'data/*.json',
            'data/mofid/*.zip',
            'data/mosaec/*.csv',
            'models/stability/*',
        ],
    },
    # The supporting-information structure archives remain usable from a
    # source checkout and through the existing on-demand downloader, but are
    # not redistributed in wheel or source-distribution artifacts without an
    # independently documented asset-level grant.
    exclude_package_data={'CoREMOF': ['data/SI/*.zip']},
    # The release loader, checker labels, parent resolver, and splitter are
    # standard-library-only. Historical scientific features remain available
    # through the compatibility-preserving ``full`` extra.
    install_requires=[],
    extras_require={
        'full': LEGACY_SCIENTIFIC_REQUIREMENTS,
        'openbabel': LEGACY_SCIENTIFIC_REQUIREMENTS + ['openbabel-wheel'],
        'docs': LEGACY_SCIENTIFIC_REQUIREMENTS + [
            'sphinx',
            'furo',
            'sphinx-design',
            'sphinx-autodoc-typehints',
        ],
        'test': LEGACY_SCIENTIFIC_REQUIREMENTS + ['ruff'],
    },
    entry_points={
        'console_scripts': [
            'coremof=CoREMOF.cli:main',
        ],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    license_files = ("LICENSE",),
    python_requires='>=3.9, <3.12',
    project_urls={
        "Homepage": "https://core-mof-tools.readthedocs.io/",
        "Repository": "https://github.com/Chung-Research-Group/CoRE-MOF-Tools",
        "Dataset splitting handbook": (
            "https://github.com/Chung-Research-Group/CoRE-MOF-Tools/"
            "blob/main/README_DATASET_SPLITTING.md"
        ),
        "Issues": "https://github.com/Chung-Research-Group/CoRE-MOF-Tools/issues",
        "PyPI": "https://pypi.org/project/CoREMOF-tools/",
    },
)
