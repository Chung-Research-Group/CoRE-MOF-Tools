"""Run one full-precision RAC3/RAC5 calculation with the frozen runtime."""
import argparse
import json

from CoREMOF.release_racs import calculate_release_racs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("cif-path", "structure-id", "output-dir", "python",
                 "source-archive", "environment-manifest"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--depth", type=int, choices=(3, 5), default=5)
    parser.add_argument("--timeout-seconds", type=int, default=1200)
    record = calculate_release_racs(**vars(parser.parse_args()))
    print(json.dumps({name: record[name] for name in (
        "structure_id", "depth", "execution_status", "available")}, indent=2))


if __name__ == "__main__":
    main()
