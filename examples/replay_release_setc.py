#!/usr/bin/env python3
"""Retired SETC calculation entry point; see read_checker_results.py.

Run ``python replay_release_setc.py --help`` for the migration notice. This
example performs no training, download, CCDC activation or feature generation.
"""
from CoREMOF.release_setc import main


if __name__ == "__main__":
    raise SystemExit(main())
