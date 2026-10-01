Reproduce the release Zeo++ calculations
=========================================

``CoREMOF.release_zeopp.calculate_release_zeopp`` runs the recorded Zeo++
0.4.7 N2/He and bonded-framework calculations through the byte-identical
scientific workers. The external ``network`` binary is supplied explicitly
and checked against its recorded SHA-256. No executable is downloaded and
no arbitrary executable on ``PATH`` is silently accepted.

Recorded scientific settings
-------------------------------

The calculations use high accuracy (``-ha``) and the binary's built-in atomic
radii. Pore diameters are probe-independent. N2-accessible surface area,
probe-occupiable volume and channel dimension use a radius of 1.655 angstrom.
Helium uses 1.32 angstrom and contributes its accessible void fraction only.
Surface-area sampling uses 5,000 samples per atom, while volume sampling uses
5,000 samples in total. The recorded CLI does not expose Monte Carlo seed
control. The implementation does not invent a seed or promise mathematical
determinism across arbitrary builds or systems.

The probe output has eighteen scalar values: four intrinsic fields, thirteen
N2 fields and one He field. Bonded-framework dimensionality is calculated
separately with ``-ha -strinfo``. It is not the dimension of probe-accessible
void channels. Retain all framework counts and channel dimensions, using
their maximum as the respective summary. A valid zero means no accessible
channel or no periodic framework, not a missing calculation.

Usage and output
-------------------

Run on a compute node with the recorded binary installed::

    from CoREMOF.release_zeopp import calculate_release_zeopp

    result = calculate_release_zeopp(
        "structure.cif",
        structure_id="2016[Co][sqc27]3[FSR]3",
        source_database="COD",  # from release metadata, not the CoRE ID
        output_dir="new_result",
        network="/runtime/zeopp/bin/network",
        timeout_seconds=300,
    )

The output parent must exist. Existing output directories are never
overwritten. Each call copies the CIF into an isolated work directory,
clears inherited Python/library overrides, bounds every external command
and its worker process group, and verifies that the source CIF is unchanged.
No atoms, charges, occupancies or coordinates are modified.

The transactional output includes ``record.json``, ``receipt.json``, separate
raw probe/framework records and execution logs. A component that fails has
null values and an explicit diagnostic. Independently successful components
are retained, with an overall ``PARTIAL`` status. No missing feature is
imputed, substituted from another probe or used to change a CR/NCR label.
Open-metal-site calculations are a separate workflow, not a Zeo++ output.

The binary and scientific workers are pinned. The historical identity of
every linked system library is not established by these records. Two real
COD fixtures reproduce all saved probe values, channel lists and bonded
framework counts. This is fixture-level evidence, not full-release numerical
certification. No public evidence bundle or benchmark assignment is changed.

The unchanged public function signatures in ``CoREMOF.calculation.Zeopp``
also remain available. Their general-purpose defaults are not the release
profile: explicitly pass the required N2/He radii. The channel parser now
handles every reported channel, and the framework command uses Zeo++'s
single-CIF output convention. Auxiliary files stay in a disposable directory
instead of beside the source CIF. Return keys and units are unchanged.

See ``examples/replay_release_zeopp.py`` or
``python -m CoREMOF.release_zeopp --help`` for command-line use.
