"""Replay the recorded Zeo++ N2/He and framework-dimension calculations.

Example::

    python examples/replay_release_zeopp.py structure.cif \
      --structure-id FSR-COD-2016-0106 --output-dir new_zeopp_result \
      --network /runtime/zeopp/bin/network

The output parent must exist, but the destination must not. The recorded
binary is required. No installation, CIF repair or release promotion occurs.
"""
from CoREMOF.release_zeopp import main


if __name__ == "__main__":
    raise SystemExit(main())
