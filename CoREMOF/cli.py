"""Command-line utilities for checking a CoREMOF installation."""

from __future__ import annotations

import argparse
from importlib.util import find_spec
import json
from pathlib import Path
import shutil
import sys

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


def split_release(args) -> int:
    """Build and write one audited train/validation/test split."""

    # Keep the normal installation check lightweight.  Dataset/split modules
    # are imported only when this subcommand is used.
    from CoREMOF.splitters import (
        OfficialSplitUnavailableError,
        split_release as build_split,
    )
    from CoREMOF.dataset import ReleaseValidationError

    try:
        result = build_split(
            args.release,
            checkers=args.checkers,
            fractions=tuple(args.fractions),
            parent_method=args.parent_method,
            leakage_guard=args.leakage_guard,
            missing_parent=args.missing_parent,
            random_state=args.random_state,
            stratify_by=tuple(args.stratify_by),
            labels=tuple(args.labels) if args.labels else None,
            sources=tuple(args.sources) if args.sources else None,
            variants=tuple(args.variants) if args.variants else None,
            metals=tuple(args.metals) if args.metals else None,
            structure_ids=tuple(args.structure_ids) if args.structure_ids else None,
            verify_cif_files=args.verify_cifs,
            official=args.official,
        )
        csv_path, receipt_path = result.write(
            args.output_directory,
            stem=args.stem,
            overwrite=args.overwrite,
        )
    except (
        OfficialSplitUnavailableError,
        ReleaseValidationError,
        FileNotFoundError,
        FileExistsError,
        KeyError,
        OSError,
        ValueError,
    ) as error:
        print("coremof split: error: {}".format(error), file=sys.stderr)
        return 2
    summary = {
        "assignment_csv": str(csv_path),
        "receipt_json": str(receipt_path),
        "dataset_version": result.dataset_version,
        "checker_view": result.checker_view,
        "parent_method": result.parent_method,
        "leakage_guard": result.leakage_guard,
        "provisional_input": result.provisional_input,
        "cif_files_verified": result.cif_files_verified,
        "counts": dict(result.counts),
        "leakage_audit": dict(result.leakage_audit),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="coremof",
        description="CoRE MOF database and analysis utilities",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("doctor", help="check optional dependencies by feature")

    split_parser = subparsers.add_parser(
        "split",
        help="make a deterministic parent-aware dataset split",
    )
    split_parser.add_argument(
        "release",
        type=Path,
        help="extracted CoRE-MOF release root (the directory containing dataset_info.json)",
    )
    split_parser.add_argument(
        "--checkers",
        choices=("3checker", "4checker", "5checker"),
        default="5checker",
        help="official checker-consensus view (default: 5checker)",
    )
    split_parser.add_argument(
        "--parent-method",
        choices=(
            "priority_main",
            "rac5",
            "mofid_v2",
            "mofid_v1",
            "rac5_zeo",
            "zeo",
            "source_id",
            "common_name",
            "identity_union",
            "none",
        ),
        default="priority_main",
        help="parent relation used to explain the split (default: priority_main)",
    )
    split_parser.add_argument(
        "--leakage-guard",
        choices=("auto", "main_union", "parent_only"),
        default="auto",
        help="split-block policy (default: auto)",
    )
    split_parser.add_argument(
        "--missing-parent",
        choices=("singleton", "exclude"),
        default="singleton",
        help="handling of unavailable selected parent evidence (default: singleton)",
    )
    split_parser.add_argument(
        "--fractions",
        nargs=3,
        type=float,
        metavar=("TRAIN", "VALIDATION", "TEST"),
        default=(0.8, 0.1, 0.1),
    )
    split_parser.add_argument(
        "--random-state",
        default="42",
        help="deterministic integer or text seed (default: 42)",
    )
    split_parser.add_argument(
        "--labels",
        nargs="+",
        choices=("CR", "NCR", "AMBIGUOUS", "UNCHECKED"),
        default=("CR", "NCR"),
    )
    split_parser.add_argument("--sources", nargs="+", choices=("COD", "CSD", "SI"))
    split_parser.add_argument("--variants", nargs="+", choices=("ASR", "FSR", "ION"))
    split_parser.add_argument("--metals", nargs="+")
    split_parser.add_argument(
        "--structure-ids",
        nargs="+",
        help="optional exact public structure IDs to include",
    )
    split_parser.add_argument(
        "--stratify-by",
        nargs="+",
        default=("label",),
        help="label and/or metadata columns used for approximate balance",
    )
    split_parser.add_argument(
        "--output-directory",
        type=Path,
        required=True,
        help="directory for the assignment CSV and JSON receipt",
    )
    split_parser.add_argument("--stem", default="coremof_split")
    split_parser.add_argument("--overwrite", action="store_true")
    split_parser.add_argument(
        "--verify-cifs",
        action="store_true",
        help="hash every CIF and compare it with the release manifest before splitting",
    )
    split_parser.add_argument(
        "--official",
        action="store_true",
        help="require an audited official manifest (fails closed when unavailable)",
    )
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "doctor":
        return doctor()
    if args.command == "split":
        return split_release(args)
    return 2
