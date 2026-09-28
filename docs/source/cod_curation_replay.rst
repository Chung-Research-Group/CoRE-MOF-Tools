Recorded COD curation workflow
==============================

``CoREMOF.cod_curation.curate_cod_cif_v41`` runs the recorded COD-workstation
v4.1 method on one privately copied raw source CIF. This route is separate
from the CSD-addition ``curate_cif_v41`` method and the legacy ``curate.clean``
API. None of these methods is a promise to reconstruct the older base release.
Retrieve archived CIFs when exact released bytes are required.

Method
------

The method parses all structural data blocks and records occupancy and atom-site
evidence. It constructs a primitive-cell candidate in P1, identifies periodic
framework components from lattice translations and removes only recognized
neutral free solvent. Its total pair margin is 0.25 angstrom, corresponding
to ASE's per-atom ``skin=0.125`` angstrom.

All-solvent removal (ASR) starts from the final free-solvent-removal (FSR)
candidate. Recorded O-, N- and S-donor rules protect ambiguous, bridging and
multidentate branches. Confirmed counterions remain and prevent ASR. Identical
FSR/ASR atom sets have one CIF with both roles. Unknown components, disorder
and failed invariants remain explicit review outcomes. Review, parser failure
and exclusion are not NCR classifications. Independent checkers are not run.

Use
---

.. code-block:: python

   from CoREMOF.cod_curation import curate_cod_cif_v41

   result = curate_cod_cif_v41(
       "source.cif",
       cod_id="7135365",
       python="/path/to/recorded/python",
       workflow_root="/path/to/COD_new",
       output_dir="new_cod_result",
       timeout_seconds=600,
   )
   print(result["execution_status"])

The equivalent CLI is ``python -m CoREMOF.cod_curation``. The runnable example
is ``examples/replay_cod_curation.py``. Supply the transferred COD workflow's
``scripts/``, policy/configuration files, vendored radius/ion tables and
``v4_1/runtime/`` model assets. The package profile binds 26 file hashes.
Missing or changed files fail before execution. All PACMAN model files are
required for charging because the historical module checks for them at import,
including models not selected by DDEC6. They are never downloaded by this API.

Use a separate Python 3.9 environment with ASE 3.25.0, NumPy 1.26.4,
pymatgen 2024.8.9, spglib 2.6.0, gemmi 0.7.0, SciPy 1.13.1,
NetworkX 3.2.1, PyCifRW 4.4.6 and PyTorch 2.7.0+cu118. The recorded
workstation used Python 3.9.21. The local audit used 3.9.23 and records that
difference. Equal package versions alone do not prove equal numerical-library
bytes or results. No dependency installation or silent fallback occurs.

The API uses one CPU thread, hides GPUs, disables external connections and
limits worker address space to 6 GiB. Its timeout covers the full worker process
group. PACMAN uses the recorded DDEC6 model, seed 0, neutralization, atom-type
averaging and ten decimal places, only after structural selection. Source CIFs
and external model files are never passed as writable work products.
``skip_charges=True`` produces inspection-only proposals, never accepted
charged candidates. The destination must not already exist.

Outputs and scientific limits
-----------------------------

``record.json`` and ``protocol_record.json`` retain statuses, child structures,
solvent/ion decisions, atom membership, review reasons and relative output paths.
``artifacts/`` retains derived CIFs and the recorded worker's supporting files.
``runtime.json``, ``execution_fingerprint.json``, ``input_manifest.jsonl`` and
``receipt.json`` bind the current run. Raw nested records retain temporary
execution paths, so use the normalized top-level protocol record to locate
exported files. No output is automatically release eligible.

Five raw-CIF controls retain all compared curation decisions and atom selections,
including counterion retention and disorder review. A returned review CIF is
byte-identical to its archive. Accepted CIFs are not all byte-identical:
coordinate roundoff and PACMAN charge differences remain recorded. In the
initial raw-CIF replay the largest charge difference is 0.0012394934 electron.

The frozen PACMAN neighbour routine uses zero numerical tolerance when excluding
self-pairs. Tests with saved and replayed geometries show near-zero self-pairs
and changes in neighbour counts under coordinate roundoff. This is a numerical
sensitivity of the recorded method, not proof that every historical target is
wrong or that CPU/GPU differences alone explain the discrepancy. A changed
neighbour rule would require a separately named method and scientific validation.
This API retains the recorded rule and never replaces historical charged CIFs,
accepted targets or frozen experiments.
