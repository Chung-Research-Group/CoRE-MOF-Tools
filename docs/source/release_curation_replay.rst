Recorded solvent-removal workflow
=================================

``CoREMOF.release_curation.curate_cif_v41`` reproduces the recorded July 2026
v4.1 CSD-addition curation method on one source CIF. It is an explicit route,
separate from the historical ``CoREMOF.curate.clean`` API. The latter keeps
its original defaults and independent cleaning branches. This method is not
a promise to recreate the older database's curation. The separate COD
workstation route is described in :doc:`cod_curation_replay`.

Method and scope
----------------

Free-solvent removal (FSR) removes selected disconnected neutral components.
All-solvent removal (ASR) additionally applies the recorded coordinated-solvent
selection, starting from the **final FSR**, not independently from the source.
Confirmed counterions are retained. An ionic candidate keeps ``ION_FSR`` and
skips ASR. Identical FSR/ASR atom sets share one physical output with both roles.
Public naming, including mapping ``ION_FSR`` to ``ION``, is a later release step.

The initial **total pair margin** is 0.25 angstrom. ASE adds its ``skin`` to
each atom's radius, so this method passes ``skin=0.125`` angstrom. If a
metal-containing component would be removed, the total margin increases by
0.05 angstrom per iteration, up to 40 iterations. Each attempt has a separate
directory. This does not change the legacy API's per-atom skin default.

The byte-identical recorded protocol includes preprocessing, atom-subset and
metal-count checks, the interim neutral-solvent formula list, and optional
neutralized PACMAN DDEC6 charging. The formula list is not a general chemical
solvent recognizer. Unknown removals, ambiguous atom mapping and other failed
invariants retain explicit ``REVIEW`` or ``ERROR`` outcomes. Neither means NCR.
Preprocessing can change the derived representation, including symmetry and
site metadata. The original source is never overwritten or presented as repaired.

Runtime and usage
-----------------

Supply the recorded external Python 3.9.23 environment and the frozen
CoREMOF source tree used by that environment. The adapter checks the Python
binary, cleaner, ion and radius tables, ASE neighbour implementation, PACMAN
code and recorded package versions before calculation. It records current
PACMAN source/model hashes. A full historical byte identity for every dependency
and model is **not** established by those checks. Scientific comparisons must
therefore be tied to the relevant saved reference cases.

The required recorded versions include ASE 3.23.0, pymatgen 2024.8.9,
gemmi 0.7.0, PACMAN-charge 1.4.2 and PyTorch 2.7.0+cu118. The frozen cleaner
also imports its historical optional packages. This route is not installed
by the dependency-minimal package. Models and dependencies must already be
present. Network downloads, installation and fallback runtimes are disabled.

.. code-block:: python

   from CoREMOF.release_curation import curate_cif_v41

   result = curate_cif_v41(
       "source.cif",
       source_id="SOURCE1",  # source record, before release IDs are assigned
       output_dir="new_curation_result",
       python="/path/to/recorded/environment/bin/python",
       runtime_coremof_root="/path/to/frozen/CoRE-MOF-Tools",
       timeout_seconds=600,
   )
   print(result["execution_status"], result["review_reasons"])

The equivalent command is ``python -m CoREMOF.release_curation`` or
``examples/replay_release_curation.py`` with ``--source-id``, ``--output-dir``,
``--python`` and ``--runtime-coremof-root``. The destination must not exist.
Use ``--skip-charges`` for a solvent-removal inspection only. Such candidates
are explicitly **not curation-stage eligible**. Even a charged, complete
curation result is not release eligible without the independent release gates.

Outputs and verification
------------------------

``record.json`` gives status, candidate roles and artifact paths relative to
the result directory. ``protocol_record.json`` retains the complete method
record, mapping, iteration history and review reasons. ``runtime.json`` and
``receipt.json`` record the method/runtime/input identity. ``artifacts/`` retains
derived CIFs and the original protocol's audit files, while ``stdout.txt`` and
``stderr.txt`` retain execution diagnostics. Original nested audit files retain
their temporary execution paths. Use the normalized top-level protocol record
to locate exported candidate files. Keep licensed raw CIFs and structure-resolved
results private unless redistribution is explicitly permitted.

The process runs on one CPU with hidden GPUs, fixed PyTorch seed 0,
deterministic algorithms and a bounded timeout. It does not assign CR/NCR,
modify a production registry or overwrite a release CIF. Four private saved
controls are used to check distinct/identical outputs, counterion retention,
adaptive margins and unknown-removal review. These controls are not a
release-wide chemical validation.
