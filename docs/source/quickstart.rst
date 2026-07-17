Quick start
===========

Check the environment
---------------------

.. code-block:: bash

   coremof doctor

Missing components are reported by feature. They do not prevent unrelated
features from working.

Query database metadata
-----------------------

.. code-block:: python

   from CoREMOF.structure import information

   record = information("CR-ASR", "2020[Cu][sql]2[ASR]1")
   print(record)

The allowed dataset names are ``CR-ASR``, ``CR-FSR``, ``CR-Ion``, and ``NCR``.
Use ``show_units=True`` to print the CR dataset unit metadata.

Calculate pore diameters
------------------------

.. code-block:: python

   from CoREMOF.calculation import Zeopp

   result = Zeopp.PoreDiameter("my_mof.cif")
   print(f"LCD = {result['LCD']} Å")
   print(f"PLD = {result['PLD']} Å")

Calculate surface area and pore volume
--------------------------------------

.. code-block:: python

   surface = Zeopp.SurfaceArea(
       "my_mof.cif",
       chan_radius=1.86,
       probe_radius=1.86,
       num_samples=10_000,
   )
   volume = Zeopp.PoreVolume(
       "my_mof.cif",
       chan_radius=1.86,
       probe_radius=1.86,
       num_samples=10_000,
   )

   # ASA order: Å² per unit cell, m²/cm³, m²/g
   print(surface["ASA"])
   # PV order: Å³ per unit cell, cm³/g
   print(volume["PV"])

Precheck a CIF
--------------

.. code-block:: python

   from CoREMOF.curate import preprocess

   job = preprocess("my_mof.cif", output_folder="result_curation")
   print(job.result_check)

The precheck writes standardized CIF files and a ``*_precheck.json`` report.

Handle a failed run
-------------------

When reporting a problem, include:

.. code-block:: bash

   coremof doctor
   python --version

Also include the minimal Python command, complete traceback, operating system,
and failing CIF when its licence and publication status allow sharing.
