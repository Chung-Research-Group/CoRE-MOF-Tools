Replaying full-precision release RACs
========================================

Use ``CoREMOF.release_racs.calculate_release_racs`` to reproduce the frozen
molSimplify 1.7.3 calculation. This explicit route is separate from the
general-purpose ``CoREMOF.calculation.RACs`` function, whose depth-3 default
and four-decimal output are retained for compatibility. Selecting depth 5
in the generic function alone does not reproduce the release method.

RACs are revised autocorrelation descriptors computed from the molecular
graph. Here, depth is the maximum number of bonds separating the paired
atoms, not a distance in angstroms. The pinned schema contains 44 descriptor
families at depths 0 through 5, giving 264 RAC5 values. The same pinned
implementation can calculate the depth-3, 176-value schema.

Frozen inputs and settings
----------------------------------------

The numerical calculation is imported from the exact molSimplify source
archive at commit ``60211676a71039f37adce57a60505f2bdaf7f184``. This includes
its atomic property tables. Source-archive SHA-256 is
``e8fe1bda764ae107f976fc3b989f4aab29949343dbf89fe87f0cca2ecdc0f883``.

The adapter checks the recorded Python 3.9.23 executable and the complete
sealed numerical environment against its supplied, hash-verified manifest.
The environment includes NumPy 1.26.4, SciPy 1.13.1, pandas 1.5.3,
NetworkX 2.8.8 and pymatgen 2023.9.10. Matching version labels without matching
runtime bytes is insufficient for this reproduction route.

The descriptor call retains these settings::

    depth=5
    graph_provided=False
    wiggle_room=1
    max_num_atoms=6000
    get_sbu_linker_bond_info=False
    surrounded_sbu_file_generation=False
    detect_1D_rod_sbu=False

No graph is supplied or repaired. Each process has isolated configuration
and temporary directories, a fixed Python hash seed and one numerical-library
thread. The original CIF is copied without changing its bytes and is checked
again afterwards. Ambient Python paths, dynamic-library overrides and user
molSimplify configuration are not inherited by the scientific process.

Usage
-----

The external source archive, sealed environment and environment manifest must
already be available. This function never downloads or installs them::

    from CoREMOF.release_racs import calculate_release_racs

    result = calculate_release_racs(
        "FSR-COD-2016-0106.cif",
        structure_id="FSR-COD-2016-0106",
        output_dir="new_rac5_replay",  # must not already exist
        python=pinned_python,
        source_archive=molsimplify_source_archive,
        environment_manifest=sealed_environment_manifest,
        depth=5,
        timeout_seconds=1200,
    )
    if result["available"]:
        vector = [result["descriptors"][name] for name in result["metadata_names"]]

``examples/replay_release_racs.py`` provides the equivalent command-line
example. Run this strict reproduction workflow on an allocated compute node.
The complete environment check reads roughly 1 GB across 44,204 recorded
filesystem entries and is repeated for each invocation. It is not suitable
for repeatedly launching from a shared login host. Runtime verification and
calculation together are covered by the wall-time limit.

Outputs and interpretation
----------------------------------------

The result records the descriptor names in canonical order, full-precision
values, and their exact hexadecimal floating-point representation. The latter
allows bit-level comparison without decimal-formatting ambiguity. Values are
never rounded, clipped or imputed, and both signs of zero are preserved.

The entire vector must validate. A scientific error or timeout leaves both
``descriptors`` and ``values_float_hex`` null, with a diagnostic. Changed or
missing runtime files raise ``ReleaseRACError`` before calculation. A timeout
during verification does not mean the descriptor calculation started.

The private output directory contains the result, calculation/worker logs,
runtime verification and a receipt binding input, source, environment, schema
and adapter hashes. Existing directories are not overwritten. This route
does not promote results, alter registered metadata, or modify benchmark
features and assignments. Agreement on selected examples does not establish
reproduction for every release structure or for separately trained models.
