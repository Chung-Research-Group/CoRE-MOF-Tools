Installation
============

Supported Python versions
-------------------------

CoRE MOF Tools 0.4 supports Python 3.9–3.11. Python 3.11 is recommended for a
new environment. The upper bound is retained while the complete scientific
feature set is validated against newer Python releases.

Recommended installation
------------------------

Version ``0.4.0.dev0`` is not yet the stable PyPI release. Create an isolated
conda environment and install this checkout to use the lightweight release
loader, checker classification, and dataset splitter:

.. code-block:: bash

   conda create -n coremof python=3.11
   conda activate coremof
   git clone https://github.com/Chung-Research-Group/CoRE-MOF-Tools.git
   cd CoRE-MOF-Tools
   python -m pip install .
   coremof doctor

Install the historical scientific feature set with the ``full`` extra:

.. code-block:: bash

   python -m pip install ".[full]"
   coremof doctor

Version 0.4 changes the clean-install dependency contract. The base is
standard-library-only; ``[full]`` preserves the dependencies installed by
default in 0.3. Existing environments normally retain already installed
packages, but new scientific-workflow environments should request ``[full]``.

The target-independent ``representative`` diversity profile used by
``data_split()`` and ``benchmark-cr-ncr`` has a narrower reproducibility extra:

.. code-block:: bash

   python -m pip install ".[benchmark]"

This installs exactly NumPy 1.26.4, scikit-learn 1.5.0, SciPy 1.13.1,
joblib 1.5.3, and threadpoolctl 3.6.0. The profile uses
complete scientific vectors without imputation, median/interquartile-range
scaling, at most 32 RAC5 principal components, and deterministic
MiniBatchKMeans strata. Missing dependencies or version drift raise an error;
the package never silently switches to a different numerical backend.

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

CSD retrieval
~~~~~~~~~~~~~

Install the licensed CSD software and its Python API using the CCDC instructions.
The package cannot supply or activate a CSD licence.

Checker results
~~~~~~~~~~~~~~~

Reading precomputed checker results and combining votes require no checker
installation or CCDC licence. The external-checker implementations and replay
workers are not distributed. Use the original software separately if new
calculations are needed. See :doc:`release_checkers_replay`.

MOFid
~~~~~

Follow the external MOFid compilation guide and install Open Babel. Confirm that
MOFid works independently before calling :mod:`CoREMOF.get_mofid`.

Heat-capacity ensemble
~~~~~~~~~~~~~~~~~~~~~~

The PyPI wheel does not include the approximately 1.3 GB heat-capacity ensemble.
Use a trusted full checkout or supply its ensemble directory to
:func:`CoREMOF.prediction.cp`. Install ``.[heat-capacity]`` in a separate
environment, not alongside ``.[full]`` or ``.[benchmark]``. The supplied
serialized models use scikit-learn 1.4.2 and XGBoost 2.0.3. The function checks
the selected complete 100-model ensembles, their hashes and model-library
versions before inference. See :doc:`features` for units and an example.

Historical stability models
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use ``.[historical-stability]`` in a separate environment for the original
bundled predictors. It deliberately differs from the ``full``, ``benchmark``
and ``heat-capacity`` dependencies. Zeo++, the original seven model/scaler
assets and the pinned external molSimplify source are also required. See
:doc:`historical_stability` for feature settings, usage and the limits of
cross-version reproduction. These models are not the newer MIT benchmark.
