"""Run the frozen release method once, without downloading or publishing data."""
import argparse
import json

from CoREMOF.release_mofid import calculate_release_mofid


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ("cif-path", "structure-id", "structure-variant", "output-dir",
                   "python", "method-manifest", "node-manifest", "node-root",
                   "source-root", "pinned-site", "mofid-site"):
        parser.add_argument("--" + option, required=True)
    parser.add_argument("--library-path", action="append", default=[])
    parser.add_argument("--existing-mofid-v1")
    parser.add_argument("--timeout-seconds", type=float, default=300)
    options = vars(parser.parse_args())
    options["library_paths"] = options.pop("library_path")
    record = calculate_release_mofid(**options)
    print(json.dumps({key: record[key] for key in (
        "structure_id", "mofid_v1_status", "mofid_v2_status", "v1_comparison")}, indent=2))


if __name__ == "__main__":
    main()
