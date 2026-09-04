from setuptools import setup, find_packages
from setuptools.command.build_py import build_py as _build_py
from pathlib import Path
import re
import shutil

here = Path(__file__).parent
readme_path = here / "README.md"
long_description = (
    readme_path.read_text(encoding="utf-8") if readme_path.exists() else ""
)
init_text = (here / "CoREMOF" / "__init__.py").read_text(encoding="utf-8")
version = re.search(r'^__version__ = "([^\"]+)"', init_text, re.MULTILINE).group(1)

LEGACY_SCIENTIFIC_REQUIREMENTS = [
    "pymatgen",
    "ase",
    "juliacall",
    "molSimplify",
    "PACMAN-charge",
    "cloudpickle",
    "joblib",
    "matminer",
    "numpy",
    "pandas",
    "PyYAML",
    "scipy",
    "xgboost",
    "scikit-learn==1.5.0",
    "gemmi==0.7.0",
    "phonopy",
    "networkx",
    "selfies",
    "mendeleev",
    "requests",
    "MOFClassifier==0.1.1",
]


class _RestrictedAssetSafeBuildPy(_build_py):
    """Rebuild package staging and exclude non-redistributable SI ZIPs."""

    def run(self):
        # ``build_py`` is incremental and can otherwise carry deleted source,
        # bytecode, or restricted assets from an older build/lib generation
        # into a new wheel. Reconstruct only this package's staging subtree.
        source_package = (here / "CoREMOF").resolve()
        staged_package = (Path(self.build_lib) / "CoREMOF").resolve()
        if (
            staged_package == source_package
            or staged_package in source_package.parents
            or source_package in staged_package.parents
        ):
            raise RuntimeError(
                "refusing to clean a build path that contains the source package: {}".format(
                    staged_package
                )
            )
        if staged_package.exists():
            shutil.rmtree(staged_package)
        super().run()
        staged_si = staged_package / "data" / "SI"
        if staged_si.is_dir():
            for archive in staged_si.glob("*.zip"):
                archive.unlink()


setup(
    name="CoREMOF_tools",
    version=version,
    author="Guobin Zhao",
    author_email="sxmzhaogb@gmail.com",
    description="Python API for CoRE MOF DB",
    license="CC-BY-4.0",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "CoREMOF": [
            "calculation/juliapkg.json",
            "data/*.json",
            "data/mofid/*.zip",
            "data/mosaec/*.csv",
            "models/stability/*",
        ],
    },
    # The supporting-information structure archives remain usable from a
    # source checkout and through the existing on-demand downloader, but are
    # not redistributed in wheel or source-distribution artifacts without an
    # independently documented asset-level grant.
    exclude_package_data={"CoREMOF": ["data/SI/*.zip"]},
    # The release loader, checker labels, parent resolver, and splitter are
    # standard-library-only. Historical scientific features remain available
    # through the compatibility-preserving ``full`` extra.
    install_requires=[],
    extras_require={
        "benchmark": [
            "numpy==1.26.4",
            "scikit-learn==1.5.0",
            "scipy==1.13.1",
            "joblib==1.5.3",
            "threadpoolctl==3.6.0",
        ],
        "full": LEGACY_SCIENTIFIC_REQUIREMENTS,
        "openbabel": LEGACY_SCIENTIFIC_REQUIREMENTS + ["openbabel-wheel"],
        "docs": LEGACY_SCIENTIFIC_REQUIREMENTS
        + [
            "sphinx",
            "furo",
            "sphinx-design",
            "sphinx-autodoc-typehints",
        ],
        "test": LEGACY_SCIENTIFIC_REQUIREMENTS + ["ruff"],
    },
    entry_points={
        "console_scripts": [
            "coremof=CoREMOF.cli:main",
        ],
    },
    cmdclass={"build_py": _RestrictedAssetSafeBuildPy},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    license_files=("LICENSE",),
    python_requires=">=3.9, <3.12",
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
