"""Use the historical stability interface with separately supplied model assets.

Their original hashes are verified before loading. These models are separate
from the later MIT benchmark.
"""
import argparse
import json
import os
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cif", type=Path)
    parser.add_argument("--models", type=Path, help="Original seven-asset model directory")
    args = parser.parse_args()
    # Local to this example process. No server environment or GPU job changes.
    os.environ["CUDA_VISIBLE_DEVICES"] = ""
    for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                 "TF_NUM_INTRAOP_THREADS", "TF_NUM_INTEROP_THREADS"):
        os.environ[name] = "1"
    from CoREMOF.prediction import stability
    result = stability(args.cif, model_directory=args.models)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
