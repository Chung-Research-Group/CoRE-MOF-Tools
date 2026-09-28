Recorded open-metal-site calculation
====================================

``CoREMOF.release_oms.calculate_release_oms`` reproduces the release's Zeo++
open-metal-site (OMS) diagnostic. It is **not** the same method as the legacy
``get_oms_file`` / ``get_oms_folder`` functions, which use the separate
``open_metal_detector`` implementation. Those APIs remain available with their
existing behaviour. Do not substitute their results for the release OMS field.

.. code-block:: python

   from CoREMOF.release_oms import calculate_release_oms

   result = calculate_release_oms(
       "structure.cif",
       structure_id="FSR-COD-2016-0106",
       output_dir="new_oms_result",
       network="/path/to/recorded/zeopp/bin/network",
   )
   print(result["execution_status"], result["open_metal_site_props"])

The exact invocation is ``network -oms INPUT.cif``. It uses neither a probe
radius nor an added ``-ha`` flag. The pinned binary is
``zeopp-lsmo=0.4.7=h27087fc_0``, using its built-in CCDC radius table, without
requiring the licensed CCDC software. The routine verifies binary, worker,
contract and CIF hashes and operates on a private copy in a bounded process.
The historical identities of all linked system libraries are not established.

Outputs are ``open_metal_site_count``, ``has_open_metal_sites`` and the
nullable ``surface_definition_A`` distance reported by this binary. Zero is
a successful scientific result. Missing, malformed or failed output is
unavailable, not zero. Counts are structure-preparation-dependent geometry
diagnostics, not proof of experimental accessibility or catalytic activity.
They do not independently assign CR/NCR or authorize release promotion.

The result directory must not exist. It contains ``record.json``, a receipt,
the full raw record, and execution logs. An equivalent command is
``python -m CoREMOF.release_oms`` or ``examples/replay_release_oms.py`` with
``--structure-id``, ``--output-dir`` and ``--network``. No backend installation,
CIF modification, production-data replacement or network download occurs.
