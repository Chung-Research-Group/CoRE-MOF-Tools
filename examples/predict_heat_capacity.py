"""Predict heat capacity with complete, locally supplied repository ensembles."""

import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cif", type=Path)
    parser.add_argument("--models", type=Path, required=True,
                        help="Trusted ensemble root containing 300/350/400 directories")
    parser.add_argument("--temperatures", nargs="+", type=int, default=[300, 350, 400],
                        choices=[300, 350, 400])
    args = parser.parse_args()
    from CoREMOF.prediction import cp
    result = cp(args.cif, T=args.temperatures, model_directory=args.models)
    print(json.dumps({"normalization": ["per gram", "per mole of atoms"],
                      "heat_capacity": result}, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
