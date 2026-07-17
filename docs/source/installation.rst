Installation
============

Supported Python versions
-------------------------

CoRE MOF Tools 0.3.1 supports Python 3.9–3.11. Python 3.11 is recommended for a
new environment. The upper bound is intentional because the pinned
scikit-learn model dependency is not compatible with every newer Python release.

Recommended installation
------------------------

Create an isolated conda environment, install the package, and run the diagnostic
command:

.. code-block:: bash

   conda create -n coremof python=3.11
   conda activate coremof
   pip install CoREMOF-tools
   coremof doctor

For repository development or exact environment reproduction:

.. code-block:: bash

   git clone https://github.com/Chung-Research-Group/CoRE-MOF-Tools.git
   cd CoRE-MOF-Tools
   conda env create -f env.yaml
   conda activate coremof_tools
   coremof doctor

Optional software by feature
----------------------------

Zeo++ geometry
~~~~~~~~~~~~~~

Install Zeo++ from conda-forge and confirm that ``network`` is on ``PATH``:

.. code-block:: bash

   conda install -c conda-forge zeopp-lsmo
   network

If your executable has another name or location, set
``COREMOF_NETWORK_EXECUTABLE`` to its path.

CSD and MOSAEC
~~~~~~~~~~~~~~

Install the licensed CSD software and its Python API using the CCDC instructions.
The package cannot supply or activate a CSD licence.

MOFChecker
~~~~~~~~~~

.. code-block:: bash

   pip install git+https://github.com/Au-4/mofchecker_2.0.git@main

MOFid
~~~~~

Follow the external MOFid compilation guide and install Open Babel. Confirm that
MOFid works independently before calling :mod:`CoREMOF.get_mofid`.

Heat-capacity ensemble
~~~~~~~~~~~~~~~~~~~~~~

The PyPI wheel does not include the approximately 1.3 GB heat-capacity ensemble.
Use a full GitHub checkout when calling :func:`CoREMOF.prediction.cp`. The function
checks the requested temperature directories and reports which models are absent.
