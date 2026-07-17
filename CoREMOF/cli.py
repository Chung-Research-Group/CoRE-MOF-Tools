"""Command-line utilities for checking a CoREMOF installation."""

from __future__ import annotations

import argparse
from importlib.util import find_spec
from pathlib import Path
import shutil

from CoREMOF import __version__


FEATURES = {
    "Database lookup": {
        "modules": ["requests", "gemmi"],
    },
    "Zeo++ geometry": {
        "executables": ["network"],
    },
    "CIF curation": {
        "modules": ["ase", "pymatgen", "scipy", "gemmi"],
    },
    "CSD / MOSAEC": {
        "modules": ["ccdc"],
    },
    "PACMAN charges": {
        "modules": ["PACMANCharge"],
    },
    "MOFClassifier": {
        "modules": ["MOFClassifier"],
    },
    "Stability prediction": {
        "modules": ["cloudpickle", "keras", "molSimplify"],
        "executables": ["network"],
        "paths": ["models/stability/water_model.pkl"],
    },
    "Heat capacity": {
        "modules": ["pandas", "matminer", "joblib"],
        "paths": [
            "models/cp_app/ensemble_models_smallML_120_100/300",
            "models/cp_app/ensemble_models_smallML_120_100/350",
            "models/cp_app/ensemble_models_smallML_120_100/400",
        ],
    },
    "MOFid": {
        "modules": ["ase", "pymatgen", "networkx", "selfies", "openbabel"],
    },
}


def _module_available(name: str) -> bool:
    try:
        return find_spec(name) is not None
    except (ImportError, ModuleNotFoundError, ValueError):
        return False


def doctor() -> int:
    """Print feature-level dependency status and return zero."""

    package_root = Path(__file__).resolve().parent
    print(f"CoREMOF {__version__}")
    print(f"Package: {package_root}")
    print()
    for feature, requirements in FEATURES.items():
        missing = [
            name
            for name in requirements.get("modules", [])
            if not _module_available(name)
        ]
        missing.extend(
            name
            for name in requirements.get("executables", [])
            if shutil.which(name) is None
        )
        missing.extend(
            f"{relative_path} (data)"
            for relative_path in requirements.get("paths", [])
            if not (package_root / relative_path).exists()
        )
        if missing:
            print(f"[MISSING] {feature}: {', '.join(missing)}")
        else:
            print(f"[OK]      {feature}")
    print()
    print("Missing optional components affect only the listed feature.")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="coremof",
        description="CoRE MOF database and analysis utilities",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("doctor", help="check optional dependencies by feature")
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "doctor":
        return doctor()
    return 2
