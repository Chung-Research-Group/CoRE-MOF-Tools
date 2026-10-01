Replaying the release MOFid method
========================================

``CoREMOF.release_mofid.calculate_release_mofid`` runs one CIF through the
frozen release protocol in a separate, time-limited process. It is an
explicit reproduction route, not a change to the general-purpose
``CoREMOF.get_mofid`` API. The latter must not be described as release-equivalent.
Neither route repairs a structure or promotes a result into release metadata.

What is frozen
--------------

The packaged ``_release_mofid_protocol.py`` is byte-identical to the project's
recorded scientific protocol (SHA-256
``54481294fbb78e94aae9f0db4baaaa879ad321309b18aecf9fe51d3c00061bd9``).
The adapter requires the original method manifest and ordered 1,182-node
library manifest, verifies library/runtime bytes in the calculating process,
and requires Python 3.9, ASE 3.22.1, pymatgen 2024.2.8, NetworkX 3.2.1,
SELFIES 2.1.1 and NumPy 1.26.4. The pinned MOFid commit is
``5c1b7d3345fc7f3aca1bb346244962ed70634fc7``. Modified Open Babel and Java
are checked against the recorded runtime, not a substitute executable.

The node matcher uses ``ltol=0.25``, ``stol=1.5``, an angle tolerance of
5 degrees, no primitive-cell reduction and no scaling. It selects the first
matching node in the frozen official archive order and records all matches.
It preserves the reference coordinate-conversion convention exactly, including
passing the ASE positions without a Cartesian flag. This is historical method
reproduction, not an endorsement of that convention for new algorithms.
The generic wrapper uses a different coordinate conversion, tolerances and
multiple-match policy. These paths must remain explicitly distinguished.

Usage
-----

Supply paths from the preserved method/runtime bundle to
``examples/replay_release_mofid.py``. Its ``--help`` lists all required inputs.
No path is taken from a particular workstation by default. For example::

    from CoREMOF.release_mofid import calculate_release_mofid

    record = calculate_release_mofid(
        "2016[Co][sqc27]3[FSR]3.cif",
        structure_id="2016[Co][sqc27]3[FSR]3",
        structure_variant="FSR",
        source_database="COD",  # from release metadata, not the CoRE ID
        output_dir="new_replay",  # must not already exist
        python=pinned_python,
        method_manifest=method_manifest,
        node_manifest=node_manifest,
        node_root=node_root,
        source_root=mofid_source,
        pinned_site=pinned_python_packages,
        mofid_site=mofid_python_packages,
        library_paths=runtime_library_paths,
        timeout_seconds=300,
    )

Run inside an appropriately limited compute allocation. The adapter limits
wall time and numerical-library threads, but it does not allocate cluster
resources. The original CIF is copied to an isolated temporary directory and
is rechecked afterwards. Existing output directories are never overwritten.
The result folder contains the result, input identity, worker log and receipt.
The receipt binds the protocol, adapter, method, node library and input hashes.

Interpretation and limits
-------------------------

* A topology ``UNKNOWN``, ``ERROR`` or ``TIMEOUT`` embedded in an otherwise
  calculated identifier is distinct from a failed calculation.
* Missing/unmatched nodes, execution errors and timeouts remain explicit.
  They are not repaired, imputed, or converted into checker FAIL labels.
* A supplied existing MOFid-v1 is retained under the frozen policy. Inspect
  ``v1_comparison``: ``SUBSTANTIVE_MISMATCH`` is not a successful revalidation.
* FSR records retain the ``COREMOF_FSR_EXTENSION`` scope.
* Runtime verification may itself time out. Such a timeout does not prove
  that the scientific calculation began.
* A successful replay of selected cases is not proof that all release records
  reproduce, and it does not clear the separate MOFid release-promotion gate.

The external binaries and node library are not bundled or downloaded. If the
recorded runtime is absent or differs, this route fails with an explicit
diagnostic instead of falling back to the generic wrapper.
