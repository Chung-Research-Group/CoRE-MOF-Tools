Reproduce the release topology calculation
==========================================

The explicit ``CoREMOF.release_topology.calculate_release_topology`` API
replays the CrystalNets 1.2.0 calculation with Julia 1.12.6. It is separate
from the unchanged legacy ``topology(..., node_type="single")`` function,
whose Julia configuration requests CrystalNets 1.0 and returns only one mode.

The release uses ``StructureType.MOF`` with **both SingleNodes and AllNodes**.
These clustering modes reduce atomic structures to nets differently. Retain
both results, every subnet and any mode disagreement. The recorded
``catenation_degree`` is the number of subnets in the CrystalNets result,
not a new calculation of physical entanglement. Unknown nets retain their
complete topological genome, not just a shortened display label.

Required external inputs
-------------------------

Supply the existing Julia executable, the recorded ``Project.toml`` and
``Manifest.toml``, and the CrystalNets depot. The package verifies the exact
worker, projection and project-file hashes. It does not install or download
Julia, packages, archives or atomic structures.

The exact project files are also included in the installed package at
``CoREMOF/_topology_profile/``. This directory may serve as ``project`` once
the recorded dependencies are available in the selected depot. For a new
installation, provision Julia and those locked dependencies separately,
then record and review its runtime snapshot before calculation.

Record a runtime snapshot once, on the compute node where it will be used::

    from CoREMOF.release_topology import capture_runtime_manifest

    digest = capture_runtime_manifest(
        julia="/runtime/julia-1.12.6/bin/julia",
        depot="/runtime/crystalnets-depot",
        output_path="runtime.json",
    )

Keep the returned SHA-256 separately with the run inputs. Each calculation
checks full membership and bytes of the Julia installation and depot package,
artifact and compiled-cache trees against this snapshot. This reads about
1.4 GB for the tested installation, so use a bounded compute allocation.
Changes require review and a distinct snapshot, never silent manifest repair.
The snapshot binds the **current** runtime. It is not retrospective proof
that every binary was identical during the original release calculation.
Exact comparison with saved scientific results supplies a separate replay
check. The current audit matches two COD successes and one partial SI result.
It is not a full-database rerun or release-promotion decision.

Calculate one structure
-------------------------

Create the output parent first and choose a new, nonexistent output directory::

    from CoREMOF.release_topology import calculate_release_topology

    result = calculate_release_topology(
        "structure.cif",
        structure_id="FSR-COD-2016-0106",
        output_dir="new_result",
        julia="/runtime/julia-1.12.6/bin/julia",
        project="/runtime/crystalnets",
        depot="/runtime/crystalnets-depot",
        runtime_manifest="runtime.json",
        runtime_manifest_sha256=digest,
        timeout_seconds=300,
    )

The CIF, project and workers run from private copies. Julia startup scripts,
history, automatic precompilation and online package access are disabled.
The timeout includes Julia startup and terminates only the new worker's
process group. No CIF repair or scientific fallback is performed.

Outputs and unavailable results
--------------------------------

The output directory is published transactionally without overwriting existing
results. It contains ``raw.json``, ``record.json``, ``receipt.json`` and the
worker's standard-output/error logs. Source CIF bytes are rechecked afterward.

``record.json`` uses the byte-identical release projection. A failed component
leaves successful subnet information intact but marks the structure
``PARTIAL`` with ``topology_available=false``. It must not supply a complete
topology fingerprint for related-structure matching. Process errors and
timeouts remain unavailable with explicit diagnostics, not valid unmatched
structures. Zero-dimensional results are retained as zero-dimensional.

No metadata, parent grouping, frozen assignment or trained model is updated
by this function. Scientific output can be compared with existing records
before a separate release-evidence decision.

``examples/replay_release_topology.py`` and
``python -m CoREMOF.release_topology --help`` provide the same command-line
interface. Julia is needed only when calculating, not when importing this
standard-library-only module or inspecting a saved record.
