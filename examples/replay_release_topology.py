"""Replay the recorded two-mode topology method with explicit external paths.

Run on a compute node, not a shared login host. The external Julia installation
and its recorded project/depot must already exist. No dependencies download.

Example::

    python examples/replay_release_topology.py structure.cif \
      --structure-id "2016[Co][sqc27]3[FSR]3" --output-dir new_result \
      --julia /runtime/julia-1.12.6/bin/julia --project /runtime/crystalnets \
      --depot /runtime/crystalnets-depot --runtime-manifest runtime.json \
      --runtime-manifest-sha256 RECORDED_SHA256
"""
from CoREMOF.release_topology import main


if __name__ == "__main__":
    raise SystemExit(main())
