from setuptools import setup, find_packages
from pathlib import Path
import re

here = Path(__file__).parent
readme_path = here / "README.md"
long_description = readme_path.read_text(encoding="utf-8") if readme_path.exists() else ""
init_text = (here / "CoREMOF" / "__init__.py").read_text(encoding="utf-8")
version = re.search(r'^__version__ = "([^\"]+)"', init_text, re.MULTILINE).group(1)

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
            'data/SI/*.zip',
            'data/mofid/*.zip',
            'data/mosaec/*.csv',
            'models/stability/*',
        ],
    },
    install_requires=[
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
        'scikit-learn==1.3.2',
        'gemmi==0.7.0',
        'phonopy',
        'networkx',
        'selfies',
        'mendeleev',
        'requests',
        "MOFClassifier==0.1.1"
    ],
    extras_require={
        'openbabel': ['openbabel-wheel'],
        'docs': ['sphinx', 'sphinx-rtd-theme', 'sphinx-autodoc-typehints'],
        'test': ['ruff'],
    },
    entry_points={
        'console_scripts': [
            'coremof=CoREMOF.cli:main',
        ],
    },
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    license_files = ("LICENSE",),
    python_requires='>=3.9, <3.12',
    project_urls={
        "Homepage": "https://coremof-tools.readthedocs.io/",
        "Repository": "https://github.com/Chung-Research-Group/CoRE-MOF-Tools",
        "Issues": "https://github.com/Chung-Research-Group/CoRE-MOF-Tools/issues",
        "PyPI": "https://pypi.org/project/CoREMOF-tools/",
    },
)
