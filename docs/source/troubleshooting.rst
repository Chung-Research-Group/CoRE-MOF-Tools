Troubleshooting
===============

``network`` was not found
-------------------------

Install ``zeopp-lsmo`` from conda-forge. If Zeo++ is installed under a custom
name, set ``COREMOF_NETWORK_EXECUTABLE`` to the executable path.

Heat-capacity models are missing
--------------------------------

The PyPI wheel intentionally excludes the approximately 1.3 GB ensemble. Clone
the full repository and run from that checkout. Confirm that directories named
``300``, ``350``, and/or ``400`` exist below
``CoREMOF/models/cp_app/ensemble_models_smallML_120_100``.

CSD or MOSAEC import failure
----------------------------

The CSD Python API must be installed in the same Python environment and a valid
licence must be available. A normal PyPI installation cannot satisfy this
requirement.

The structure name is absent from a feature CSV
------------------------------------------------

Use the exact CIF filename, including ``.cif``, as ``structure_name``. Current
feature files contain a dedicated ``structure_name`` column; legacy files with
an ``Unnamed: 0`` site-index column are also accepted.

A CIF path contains spaces
--------------------------

Current Zeo++ wrappers pass arguments without a shell and support spaces. If a
failure persists, run ``coremof doctor`` and include the full Zeo++ diagnostic
reported by the exception.

Read the Docs import warnings
-----------------------------

The documentation mocks optional scientific packages and does not execute
scientific calculations. A successful documentation build validates API
structure and cross-references, not numerical functionality.
