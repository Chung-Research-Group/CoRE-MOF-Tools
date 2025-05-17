from setuptools import setup, find_packages
from pathlib import Path

here = Path(__file__).parent
readme_path = here / "README.md"
long_description = readme_path.read_text(encoding="utf-8") if readme_path.exists() else ""

setup(
    name='CoREMOF_tools',
    version='0.1.6',
    author='Guobin Zhao',
    author_email='sxmzhaogb@gmail.com',
    description='Python API for CoRE MOF DB',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=find_packages(),
    install_requires=[
        'pymatgen',
        'ase',
        'juliacall',
        'molSimplify',
        'PACMAN-charge',
        'cloudpickle',
        'matminer',
        'xgboost',
        'scikit-learn==1.3.2',
        'mofchecker',
        'gemmi==0.7.0',
        'phonopy',
        'networkx',
        'selfies',
        'mendeleev',
        'requests'
    ],
    extras_require={
        'zeopp': ['zeopp-lsmo'],
        'openbabel': ['openbabel-wheel']
    },
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.9',
    ],
    python_requires='>=3.9, <4',
    project_urls={
        "Homepage": "https://coremof-tools.readthedocs.io/",
        "Repository": "https://github.com/sxm13/CoREMOF_tools",
        "Issues": "https://github.com/mtap-research/CoRE-MOF-Tools/issues",
        "PyPI": "https://pypi.org/project/CoREMOF-tools/",
    },
)
