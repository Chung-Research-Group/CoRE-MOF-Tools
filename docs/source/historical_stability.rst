Historical stability predictors
===============================

``CoREMOF.prediction.stability`` applies the original bundled solvent-removal
and thermal neural networks and water-stability random forest. It is **not**
the later CoREMOF-COD multi-seed MIT benchmark, and does not reproduce those
benchmark predictions. It does not train a model or replace release metadata.
The original studies are `Nandy et al. (2021)
<https://doi.org/10.1021/jacs.1c07217>`_ and the `water-stability study (2024)
<https://doi.org/10.1021/jacs.4c05879>`_.

Installation and use
--------------------

Use a separate environment with ``pip install ".[historical-stability]"``.
Do not combine this extra with ``full``, ``benchmark`` or ``heat-capacity``:
their model-library and pandas versions differ. The tested inference stack
uses TensorFlow CPU 2.18.0, Keras 3.8.0, scikit-learn 1.3.0, NumPy 1.26.4,
pandas 1.5.3 and pymatgen 2024.2.8. The extra does not supply Zeo++ or the
external molSimplify source.

The compatibility audit used unmodified molSimplify source commit
``60211676a71039f37adce57a60505f2bdaf7f184`` (version 1.7.3). Install that
trusted source without changing its code or dependency versions. This version
is not available as a PyPI 1.7.3 distribution. Its dataframe append operation
requires pandas 1.x. Set ``COREMOF_NETWORK_EXECUTABLE`` to a tested Zeo++
``network`` binary. Changing the descriptor implementation or probe settings
changes model inputs, even when the output remains a plausible number.

.. code-block:: python

   from CoREMOF.prediction import stability

   result = stability("example.cif", model_directory="original_stability_models")
   print(result["thermal stability"])  # degrees Celsius
   print(result["solvent removal probability"])
   print(result["water probability"])

The optional directory contains the seven original model/scaler files.
Omit it when the trusted checkout already contains them. Missing or changed
assets raise an error before deserialization. The API verifies their SHA-256
hashes, copies the verified bytes and the input CIF privately, and never
substitutes later models, edits the source CIF or imputes missing features.
The legacy result keys and ``"nan, °C, nan"`` unit string remain unchanged.
The two probabilities are dimensionless, not checker PASS/FAIL classifications.

``examples/predict_historical_stability.py`` provides a CPU-only command-line
example. On an HPC system, run it inside an appropriately bounded allocation.
If a global ``LD_LIBRARY_PATH`` points to another TensorFlow installation,
start a fresh process with that conflicting entry removed. Do not change
system libraries or remove paths required by unrelated running work.

Exact feature contract
----------------------

Both neural models use 134 selected revised autocorrelation (RAC) descriptors
through depth three and 14 geometric values. ``RACs(cif, depth=3)`` retains its
historical four-decimal rounding. The solvent-removal model places geometry
first, whereas the thermal model places RAC descriptors first. Geometry is
LCD, PLD, LFPD, accessible/inaccessible/total gravimetric pore volumes,
gravimetric accessible surface area, accessible pore volume and fraction,
inaccessible pore volume and fraction, their summed fraction, volumetric
accessible surface area and unit-cell volume, in that order.

Surface/volume sampling uses a 1.86 Å channel/probe radius, 10,000 samples
and high accuracy. The water model uses 11 selected RAC descriptors and
accessible gravimetric surface area at 1.4 Å, also with 10,000 samples and
high accuracy. These are **not** the release N2-probe settings at 1.655 Å.
The neural input widths are 148 and the water-model input width is 12.
All inputs, scaled inputs and predictions must be finite. Probabilities must
lie in [0, 1]. Zeros are retained. Thermal predictions retain the original
one-decimal NumPy rounding without changing the model output dtype.

Reproducibility limits
----------------------

The neural HDF5 files record Keras 2.3.0. Their three scalers record
scikit-learn 0.22.1, whereas the water model/scaler record 1.3.0. Thus a single
modern environment is a compatibility environment, not proof of the original
training runtime. The `scikit-learn persistence guidance
<https://scikit-learn.org/stable/model_persistence.html>`_ does not guarantee
cross-version loading. Its version warnings remain visible. Keras models are
loaded with ``compile=False`` for inference, preserving their weights without
reconstructing unused training metrics or optimizer state, as supported by
the `Keras loading API
<https://keras.io/api/models/model_saving_apis/model_saving_and_loading/>`_.

The archived example notebook and historical metadata already differ for
the same old SI identifier. Preserve those distinct references. A successful
compatibility replay must not be described as exact reproduction of every
historical release value, nor may it overwrite accepted metadata or frozen
benchmark results.
