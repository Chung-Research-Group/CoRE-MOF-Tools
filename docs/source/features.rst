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

``mof_check`` combines the Chen–Manz checks with MOFChecker. ``run_MOSAEC`` and
``run_mofclassifier`` are batch functions and validate that their external
dependencies and input CIFs are present before starting.

Predictions
-----------

:func:`CoREMOF.prediction.pacman` writes a charge-annotated CIF.
:func:`CoREMOF.prediction.stability` combines pretrained models with Zeo++ and
RAC descriptors. :func:`CoREMOF.prediction.cp` uses temperature-specific model
ensembles from the full repository.

The reported ensemble ``std`` is model-to-model dispersion. It is not, by
itself, a calibrated prediction interval or a guarantee of applicability to an
out-of-distribution structure.

Concurrency and temporary files
-------------------------------

Zeo++, OMS, RAC, and heat-capacity workflows use unique temporary paths. This
prevents one process from overwriting another in a shared working directory.
Output directories supplied by the user remain the user's responsibility.
