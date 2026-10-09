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
loader, checker classification, and dataset splitter. Use the audited checkout
supplied with your handoff; a fresh clone must use the exact fork and checkpoint
named by that handoff. The upstream default branch and an unversioned PyPI
installation are not a substitute for that checkpoint:

.. code-block:: bash

   conda create -n coremof python=3.11
   conda activate coremof
   cd /path/to/audited/CoRE-MOF-Tools
   python -m pip install .
   coremof doctor

Alternatively, install an independently verified wheel supplied by the project:

.. code-block:: bash

   python -m pip install /path/to/coremof_tools-0.4.0.dev0-py3-none-any.whl
   coremof doctor

The wheel installs the Python API and package resources. The human guides,
example scripts, companion notebook, and portable agent skill are available in
the source checkout or an extracted source distribution, not in the wheel's
installed files. Open ``examples/CoREMOF_dataset_splitting_quickstart.ipynb``
from that source tree. Agents can use
``.agents/skills/coremof-dataset-use/SKILL.md`` there. This consumer guide links
to the existing metadata, checker-result and dataset examples. Internal
development and curation skills are not distributed. Installing
a source distribution also does not copy those guides into ``site-packages``.

For a fresh clone, the verified public fork API checkpoint is:

.. code-block:: bash

   git clone https://github.com/DrakeChan/CoRE-MOF-Tools.git
   cd CoRE-MOF-Tools
   git checkout --detach c66796b77b0e86775f43b9021c6bcb3ebd93abfa
   python -m pip install .
   coremof doctor

That commit provides the ``0.4.0.dev0`` API. Later documentation follow-ups in
a supplied working tree are not part of this commit until separately published;
use the exact source or wheel hashes recorded by their handoff. The hosted
Read the Docs site and stable PyPI package describe earlier published versions.
Read the guides in the source tree matching your installed checkpoint.

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

For the historical scientific environment recorded in ``env.yaml``:

.. code-block:: bash

   cd /path/to/audited/CoRE-MOF-Tools
   conda env create -f env.yaml
   conda activate coremof_tools
   python -m pip install --no-deps -e .
   coremof doctor

This environment file does not define the pinned ``benchmark`` or predictor
environments above. Use the extra matching the requested workflow.

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
historical predictor interfaces. It deliberately differs from the ``full``,
``benchmark`` and ``heat-capacity`` dependencies. Obtain the original seven
model/scaler assets separately; the API verifies their hashes before loading.
Zeo++ and the pinned external molSimplify source are also required. See
:doc:`historical_stability` for feature settings, usage and the limits of
cross-version reproduction. These models are not the newer MIT benchmark.
