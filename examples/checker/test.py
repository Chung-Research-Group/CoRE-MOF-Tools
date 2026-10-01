"""Read saved checker results through the current results-only example.

Checker engines are not distributed or launched here. Supply an explicit
authorized release directory; no old example CIF or checker output is assumed.
New calculations require the original checker software obtained separately.
"""
import argparse
from pathlib import Path
import runpy


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("release_root", help="Extracted authorized release directory")
    parser.add_argument("--structure-id")
    parser.add_argument("--checkers", nargs="+")
    args = parser.parse_args(argv)
    forwarded = [args.release_root]
    if args.structure_id:
        forwarded.extend(("--structure-id", args.structure_id))
    if args.checkers:
        forwarded.extend(("--checkers", *args.checkers))
    # Resolve package imports only after parsing, so --help needs no installation.
    reader = runpy.run_path(str(Path(__file__).resolve().parents[1] / "read_checker_results.py"))
    return reader["main"](forwarded)


if __name__ == "__main__":
    raise SystemExit(main())
