from setuptools import setup, find_packages
from pathlib import Path

setup(
        name='CoREMOF_tools',
        version='0.1.3',
        author='Guobin Zhao',
        author_email='sxmzhaogb@gmail.com',
        description='Python API for CoRE MOF 2024 DB',
        readme_file = Path(__file__).parent / "README.md"
        long_description = readme_file.read_text() if readme_file.exists() else ""
        long_description_content_type='text/markdown',
        # url='https://github.com/sxm13/CoREMOF_tools',
        packages=find_packages(),
        # package_data={
        # 'CoREMOF': [
        #     'data/*.json',
        #     'data/CSD/*.json',
        #     'data/SI/*.zip',
        #     'models/cp_app/ensemble_models_smallML_120_100/300/*',
        #     'models/cp_app/ensemble_models_smallML_120_100/350/*',
        #     'models/cp_app/ensemble_models_smallML_120_100/400/*',
        #     'models/stability/*'
        # ],
        # },
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
            'mendeleev'
        ],
        extras_require={
            'zeopp': ['zeopp-lsmo'],
            "openbabel": ["openbabel-wheel"]
        },
        classifiers=[
            'Development Status :: 6 - Mature',
            'Intended Audience :: Developers',
            'Topic :: Software Development :: Build Tools',
            'License :: OSI Approved :: MIT License',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.9',
        ],
        python_requires='>=3.9, <4',
        # entry_points={
        #     'console_scripts': [
        #         'coremof=CoREMOF:curate',
            # ],
        # },
        project_urls={
                        "Homepage": "https://coremof-tools.readthedocs.io/en/latest/index.html#",
                        "Documentation": "https://coremof-tools.readthedocs.io/en/latest/index.html#",
                        "Repository": "https://github.com/sxm13/CoREMOF_tools",
                        "Issues": "https://github.com/mtap-research/CoRE-MOF-Tools/issues",
                        "PyPI": "https://pypi.org/project/CoREMOF-tools/",
                        },
    )
