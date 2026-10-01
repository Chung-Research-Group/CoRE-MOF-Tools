Feature guide
=============

Database and structures
-----------------------

Use :func:`CoREMOF.structure.information` for metadata lookup. SI archives are
extracted by instantiating :class:`CoREMOF.structure.download_from_SI`. CSD
downloads require the licensed CSD Python API and use
:func:`CoREMOF.structure.download_from_CSD`.

Geometry and descriptors
------------------------

The :mod:`CoREMOF.calculation.Zeopp` functions return dictionaries with an
explicit ``unit`` field. ``SurfaceArea`` returns three values for each area in
the order Å² per unit cell, m²/cm³, and m²/g. ``PoreVolume`` returns Å³ per unit
cell and cm³/g. Probe and channel radii are in ångström.

Basic crystallographic descriptors are in
:mod:`CoREMOF.calculation.mof_features`. Topology, open-metal-site, and RAC
calculations have additional Julia, OMS, or molSimplify requirements.

Curation and validation
-----------------------

The classes in :mod:`CoREMOF.curate` execute their workflow during
initialization and write results to ``output_folder``. Keep the returned object
when you need in-memory status such as ``preprocess.result_check``.

For checker-based dataset selection, load precomputed release results with
``CoREMOFDataset.from_release`` and call ``dataset.classify``. This reads
recorded votes without running external checkers. Detailed findings and raw
scores remain in the release's ``metadata/checker_findings.csv`` or ``.jsonl``.
Missing results remain unavailable, not FAIL.

``mof_check``, ``run_MOSAEC`` and the external-checker replay interfaces are
retired, results-only migration notices. Their former implementations are not
distributed. See :doc:`release_checkers_replay` and the repository README.
The optional ``run_mofclassifier`` integration uses our separately installed
MOFClassifier and is independent of reading existing checker results.

Predictions
-----------

:func:`CoREMOF.prediction.pacman` writes a charge-annotated CIF.
:func:`CoREMOF.prediction.stability` combines pretrained models with Zeo++ and
RAC descriptors. :func:`CoREMOF.prediction.cp` uses temperature-specific model
ensembles from the full repository.

Heat capacity
~~~~~~~~~~~~~

Install ``CoREMOF-tools[heat-capacity]`` in a separate Python 3.9–3.11
environment. Its supplied models record scikit-learn 1.4.2 and XGBoost 2.0.3, whereas the
``benchmark`` and legacy ``full`` extras select 1.5.0. Do not combine these
extras in one environment. The model files are not included in the wheel.
Use the complete trusted local repository ensembles, or provide their root
with ``cp(cif, T=[300], model_directory=...)``. The root contains ``300``,
``350`` and ``400`` directories, each with ``model_0`` through ``model_99``.
The API verifies the selected complete ensembles against shipped size/SHA-256
records before model loading. It does not fetch models or accept a partial
ensemble. These records identify the supplied repository assets, not a proven
byte-level match to the original publication's model files.

The method predicts harmonic constant-volume heat capacity, as described by
`Moosavi et al. (2022) <https://doi.org/10.1038/s41563-022-01374-3>`_.
The legacy function name ``cp`` and return-unit strings are unchanged.
The first returned value is in J/g/K. The second, labelled J/mol/K, is
normalized per mole of atoms, not per mole of framework formula units.
For each model, atomic contributions are summed and divided by the sum of
atomic weights or the atom count, respectively. The returned mean and
population standard deviation use all 100 models at that temperature.
No constant-pressure correction is calculated.

Only whole-kelvin temperatures with supplied models are accepted. An empty,
duplicate, fractional or non-finite temperature selection is rejected.
Missing/non-finite atomic features, invalid atomic weights and incomplete
model predictions raise errors rather than being filled or silently omitted.
The source CIF is preserved and featurization uses an isolated copy.
See ``examples/predict_heat_capacity.py`` for a local-model invocation.

Other predictors
~~~~~~~~~~~~~~~~

Each predictor loads only its own optional dependencies. Importing
``CoREMOF.prediction`` does not initialize Keras, PACMAN or the heat-capacity
featurizer. Missing packages still prevent the corresponding calculation,
and no alternative model is substituted.

``stability()`` is the compatibility interface to the historical
models and their original scalers, obtained separately and hash-verified before
use. It does not load the later CoREMOF-COD
multi-seed benchmarking models. Reproducing those experiments requires their
separately recorded inputs, assignments, fitted models and environments.
Use the separate ``historical-stability`` extra and the original hash-verified
assets. Probe radii, descriptor order and cross-version limits are documented
in :doc:`historical_stability`.
``pacman()`` isolates the input CIF and refuses an existing output, but its
generic result is not a certification for a release or simulation protocol.
Validate derived structures and recorded charges before such use.

The reported ensemble ``std`` is model-to-model dispersion. It is not, by
itself, a calibrated prediction interval or a guarantee of applicability to an
out-of-distribution structure.

The low-level ``prediction.download_file(url, path, expected_sha256=...)``
helper streams to a temporary file and publishes only after download and
checksum checks succeed. It never replaces an existing file and validates an
existing cache when a hash is supplied. Obtain that hash from a separately
trusted record. Omitting it preserves legacy usage but does not verify model
identity. Pickle/joblib model files can execute code when loaded, so use only
trusted model sources. The helper itself never loads or executes a model.

Concurrency and temporary files
-------------------------------

Zeo++, OMS, RAC, and heat-capacity workflows use unique temporary paths. This
prevents one process from overwriting another in a shared working directory.
Output directories supplied by the user remain the user's responsibility.

Frozen manuscript workflows
----------------------------

For the exact common-input benchmark assignment replay, see
:doc:`frozen_assignment_replay`. This is separate from recomputing descriptors
or running the historical property predictors with separately obtained,
hash-verified model assets.
