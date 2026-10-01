#!/usr/bin/env python3
"""Fetch a version from an explicitly supplied, checksum-pinned release catalog."""

import argparse
from pathlib import Path
import sys

_ROOT = Path(__file__).resolve().parents[1]
if (_ROOT / "CoREMOF").is_dir() and str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from CoREMOF.retrieval import fetch_release


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, epilog=(
        "No hosted catalog or release URL is bundled. Obtain the manifest and its "
        "trusted SHA-256 from the provider entitled to supply your selected data."))
    parser.add_argument("catalog", help="local catalog JSON path or HTTP(S)/file URL")
    parser.add_argument("version", help="exact catalog version; no implicit latest")
    parser.add_argument("destination", type=Path, help="new local release directory")
    parser.add_argument("--catalog-sha256", required=True)
    parser.add_argument("--verify-cifs", action="store_true")
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--max-total-bytes", type=int)
    args = parser.parse_args(argv)
    try:
        path = fetch_release(args.catalog, args.version, args.destination,
                             catalog_sha256=args.catalog_sha256, verify_cif_files=args.verify_cifs,
                             timeout=args.timeout, max_total_bytes=args.max_total_bytes)
    except (OSError, TypeError, ValueError) as error:
        parser.exit(2, "error: {}\n".format(error))
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
