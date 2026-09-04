"""Command-line utilities for checking a CoREMOF installation."""

from __future__ import annotations

import argparse
import csv
from importlib import metadata as importlib_metadata
from importlib.util import find_spec
import json
from pathlib import Path
import shutil
import sys

from CoREMOF import __version__
from CoREMOF.parents import LEAKAGE_GUARD_CHOICES, SELECTABLE_PARENT_METHODS


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
    "Representative benchmark": {
        "versions": {
            "numpy": "1.26.4",
            "scikit-learn": "1.5.0",
            "scipy": "1.13.1",
            "joblib": "1.5.3",
            "threadpoolctl": "3.6.0",
        },
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
        for distribution, expected in requirements.get("versions", {}).items():
            try:
                found = importlib_metadata.version(distribution)
            except importlib_metadata.PackageNotFoundError:
                missing.append("{}=={} (not installed)".format(distribution, expected))
            else:
                if found != expected:
                    missing.append(
                        "{}=={} (found {})".format(distribution, expected, found)
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
    from CoREMOF.targets import TargetDataError, merge_targets_from_config

    try:
        release_input = args.release
        if args.target_config is not None:
            release_input = merge_targets_from_config(
                args.release,
                args.target_config,
                verify_cif_files=args.verify_cifs,
            )
        elif args.require_target:
            raise TargetDataError(
                "--require-target requires --target-config so targets are joined before splitting"
            )
        result = build_split(
            release_input,
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
            required_targets=(
                tuple(args.require_target) if args.require_target else None
            ),
            required_target_mode=args.required_target_mode,
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
        TargetDataError,
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
        "requested_leakage_guard": getattr(
            result, "requested_leakage_guard", result.leakage_guard
        ),
        "leakage_guard": result.leakage_guard,
        "provisional_input": result.provisional_input,
        "cif_files_verified": result.cif_files_verified,
        "counts": dict(result.counts),
        "leakage_audit": dict(result.leakage_audit),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


def merge_target_data(args) -> int:
    """Join configured target files to release metadata and feature tables."""

    from CoREMOF.dataset import ReleaseValidationError
    from CoREMOF.targets import TargetDataError, merge_targets_from_config

    try:
        merged = merge_targets_from_config(
            args.release,
            args.config,
            verify_cif_files=args.verify_cifs,
        )
        csv_path, provenance_path, receipt_path = merged.write(
            args.output_directory,
            stem=args.stem,
            overwrite=args.overwrite,
        )
    except (
        FileExistsError,
        FileNotFoundError,
        OSError,
        ReleaseValidationError,
        TargetDataError,
        ValueError,
    ) as error:
        print("coremof merge-targets: error: {}".format(error), file=sys.stderr)
        return 2
    summary = {
        "merged_csv": str(csv_path),
        "provenance_jsonl": str(provenance_path),
        "receipt_json": str(receipt_path),
        "dataset_version": merged.dataset_version,
        "release_structure_count": len(merged),
        "target_columns": list(merged.target_columns),
        "feature_column_count": len(merged.feature_columns),
        "target_values_sha256": merged.receipt()["target_values_sha256"],
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


def benchmark_cr_ncr(args) -> int:
    """Build and transactionally write the paired exploratory benchmark."""

    from CoREMOF.benchmarks import (
        BenchmarkError,
        build_cr_ncr_benchmark,
    )
    from CoREMOF.dataset import CoREMOFDataset, ReleaseValidationError

    try:
        dataset = CoREMOFDataset.from_release(
            args.release, verify_cif_files=args.verify_cifs
        )
        classified = dataset.classify("5checker")
        suite = build_cr_ncr_benchmark(
            classified,
            ncr_pool_fractions=tuple(args.ncr_pool_fractions),
            seeds=tuple(args.seeds),
            total_size="full_cr",
            train=args.fractions[0],
            val=args.fractions[1],
            test=args.fractions[2],
            group_criteria=tuple(args.group_criteria),
            cohort_eligibility=args.cohort_eligibility,
            diversity=args.diversity,
            test_policy="fixed_pure_cr",
            include_full_cr_diagnostic=not args.no_full_cr_diagnostic,
        )
        output_root = suite.write(
            args.output_directory,
            stem=args.stem,
            overwrite=args.overwrite,
        )
    except (
        BenchmarkError,
        FileExistsError,
        FileNotFoundError,
        OSError,
        ReleaseValidationError,
        TypeError,
        ValueError,
    ) as error:
        print("coremof benchmark-cr-ncr: error: {}".format(error), file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                "output_root": str(output_root),
                "dataset_version": suite.dataset_version,
                "checker_view": suite.checker_view,
                "run_count": len(suite.runs),
                "fixed_test_count": len(suite.fixed_test_ids),
                "group_criteria": list(suite.group_criteria),
                "diversity_index_sha256": suite.diversity_index_hash,
                "suite_assignment_sha256": suite.assignment_digest,
                "official_split": suite.official_split,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


def attach_target_data(args) -> int:
    """Attach configured targets to one persisted assignment manifest."""

    from CoREMOF.attachments import (
        TargetAttachmentError,
        attach_targets,
        frozen_assignment_manifest,
    )
    from CoREMOF.dataset import CoREMOFDataset, ReleaseValidationError
    from CoREMOF.targets import TargetDataError

    try:
        with args.manifest.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            fields = tuple(reader.fieldnames or ())
            if "structure_id" not in fields:
                raise TargetAttachmentError(
                    "assignment manifest has no structure_id column"
                )
            rows = list(reader)
        receipt_path = args.receipt or args.manifest.with_suffix(".json")
        with receipt_path.open("r", encoding="utf-8") as handle:
            assignment_receipt = json.load(handle)
        frozen_assignment = frozen_assignment_manifest(rows, assignment_receipt)
        dataset = CoREMOFDataset.from_release(
            args.release, verify_cif_files=args.verify_cifs
        )
        attached = attach_targets(
            frozen_assignment,
            args.config,
            dataset=dataset,
            missing=args.missing,
            verify_cif_files=args.verify_cifs,
        )
        csv_path, provenance_path, receipt_path = attached.write(
            args.output_directory,
            stem=args.stem,
            overwrite=args.overwrite,
        )
    except (
        FileExistsError,
        FileNotFoundError,
        OSError,
        ReleaseValidationError,
        TargetAttachmentError,
        TargetDataError,
        TypeError,
        ValueError,
    ) as error:
        print("coremof attach-targets: error: {}".format(error), file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                "attached_csv": str(csv_path),
                "provenance_jsonl": str(provenance_path),
                "receipt_json": str(receipt_path),
                "missing_policy": attached.missing_policy,
                "selected_structure_count": len(attached.structure_ids),
                "dropped_structure_count": len(attached.dropped_ids),
                "target_columns": list(attached.target_columns),
                "original_assignment_sha256": attached.original_assignment_digest,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="coremof",
        description="CoRE MOF database and analysis utilities",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("doctor", help="check optional dependencies by feature")

    target_parser = subparsers.add_parser(
        "merge-targets",
        help="join user target files to release metadata and selected features",
    )
    target_parser.add_argument(
        "release",
        type=Path,
        help="extracted CoRE-MOF release root",
    )
    target_parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="JSON file declaring target sources, units, conditions, and aliases",
    )
    target_parser.add_argument("--output-directory", type=Path, required=True)
    target_parser.add_argument("--stem", default="coremof_targets")
    target_parser.add_argument("--overwrite", action="store_true")
    target_parser.add_argument(
        "--verify-cifs",
        action="store_true",
        help="hash every CIF before joining targets",
    )

    attach_parser = subparsers.add_parser(
        "attach-targets",
        help="attach typed targets after a frozen split without resplitting",
        description=(
            "Verify the assignment CSV against its paired receipt digest and exact "
            "release version/universe/input-hash binding, then attach targets after "
            "selection by an exact current public structure ID "
            "or a hash-bound alias declared in the target config. keep performs a "
            "left join and retains nulls; error requires every target; drop creates "
            "only a filtered derived view and never refills, rebalances, or resplits."
        ),
    )
    attach_parser.add_argument(
        "release", type=Path, help="extracted CoRE-MOF release root"
    )
    attach_parser.add_argument(
        "--manifest",
        type=Path,
        required=True,
        help="assignment CSV with structure_id and partition (or split) columns",
    )
    attach_parser.add_argument(
        "--receipt",
        type=Path,
        help=(
            "paired assignment receipt JSON used to verify the frozen digest "
            "(default: MANIFEST with a .json suffix)"
        ),
    )
    attach_parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="JSON target config using the same typed units/conditions/alias contract as merge-targets",
    )
    attach_parser.add_argument(
        "--missing",
        choices=("keep", "error", "drop"),
        default="keep",
        help="null-target policy defined above (default: keep)",
    )
    attach_parser.add_argument("--output-directory", type=Path, required=True)
    attach_parser.add_argument("--stem", default="coremof_attached_targets")
    attach_parser.add_argument("--overwrite", action="store_true")
    attach_parser.add_argument("--verify-cifs", action="store_true")

    benchmark_parser = subparsers.add_parser(
        "benchmark-cr-ncr",
        help="build paired fixed-size strict five-checker CR/NCR cohorts",
        description=(
            "Build exploratory target-independent strict five-checker cohorts. Strict "
            "CR has five available PASS votes; strict NCR has five available FAIL "
            "votes, while NOT_AVAILABLE is a non-vote. After the declared cohort-"
            "eligibility policy, full_cr fixes total size to the eligible CR count C: "
            "at NCR-pool fraction q, round-half-up(q*M) of the eligible NCR pool M "
            "replaces the same number of eligible CR rows, so q=1 uses every eligible "
            "NCR row rather than making a 100%-NCR composition. Raw strict-pool and "
            "policy-excluded counts are reported separately. "
            "fixed_pure_cr reserves one common approximately requested-size test from "
            "complete-release effective leakage blocks containing only CR. An effective "
            "block is the connected-component closure of main_union plus every ordered "
            "selected criterion before filtering. main_union is a leakage guard, not a "
            "parent, not an explanatory relation, and not identity proof; it transitively combines "
            "full CIF SHA-256, database-namespaced source siblings, release-authorized "
            "RAC5 groups, MOFid-v2 groups, and MOFid-v1 groups over the complete release; "
            "it is a partition guard, not identity or parent proof, and missing criterion "
            "evidence is a per-structure singleton. "
            "representative diversity uses no target: complete 264-value finite RAC5, "
            "otherwise 13 intensive N2/He Zeo++ fields plus channel/framework dimensions, "
            "otherwise an explicit no-numeric tier; it uses median/IQR scaling without "
            "imputation, at most 32 RAC5 principal components, and deterministic "
            "MiniBatchKMeans from the pinned benchmark extra."
        ),
    )
    benchmark_parser.add_argument(
        "release", type=Path, help="finalized extracted release root"
    )
    benchmark_parser.add_argument(
        "--ncr-pool-fractions",
        nargs="+",
        type=float,
        default=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    )
    benchmark_parser.add_argument(
        "--seeds", nargs="+", type=int, default=(42, 43, 44, 45, 46)
    )
    benchmark_parser.add_argument(
        "--group-criteria",
        nargs="+",
        default=("priority_main",),
        help=(
            "ordered canonical criteria or their documented short aliases. priority_main is the "
            "full-release conflict-aware RAC5 then MOFid-v2 then MOFid-v1 explanatory "
            "hierarchy: a lower group touching zero stronger components creates one, "
            "one attaches unresolved members, and two or more records "
            "PARENT_METHOD_CONFLICT without merging them; missing rows are singletons. "
            "It excludes Zeo++, CrystalNets, source ID/name, CIF hash, provisional "
            "identity unions, and StructureMatcher. RT is exact complete 264-value "
            "finite RAC5 plus a complete current-success CrystalNets fingerprint: "
            "network dimension; interpenetrated-subnet, catenation, and subnet counts; "
            "top-level single/all-node nets and agreement; and every SingleNodes and "
            "AllNodes subnet status, dimension, topology key/name, topological genome, "
            "and agreement. M2T is exact complete release-authorized MOFid-v2 plus that "
            "fingerprint. In order, its text conversion collapses Unicode whitespace to "
            "ASCII spaces, trims, rejects empty/whole-field placeholders, applies NFKC, "
            "then case-folds; it changes no structure. Eligible statuses are SUCCESS, "
            "SUCCESS_TOPOLOGY_UNKNOWN, SUCCESS_TOPOLOGY_ERROR, and "
            "SUCCESS_TOPOLOGY_TIMEOUT; every other status, including other execution "
            "or NOT_AVAILABLE states, adds no edge. M2T remains provisional with "
            "provisional MOFid. SM is the strict StructureMatcher convenience component "
            "of authoritative direct pairs whose "
            "forward/reverse fit(..., symmetric=True) pass under Python 3.9, pymatgen "
            "2024.2.8, NumPy 1.26.4, ElementComparator, ltol=stol=0.001, "
            "angle_tol=0.01, primitive_cell=true, scale=false, attempt_supercell=true, "
            "allow_subset=false, supercell_size=num_sites, and no ignored species. Its "
            "direct CifParser uses site/coordinate tolerances 0.0001, declared-symmetry "
            "expansion, within-tolerance merging, one-third/two-thirds rounding, "
            "occupancy checks, sorted periodic Structures, native disorder, and no "
            "manual repair/deletion/chemistry edit. Direct edges are authoritative; "
            "parser, timeout, OOM, matcher, asymmetric, or execution failures are "
            "NOT_AVAILABLE, not unmatched. Incomplete evidence adds no match; RT, M2T, "
            "and SM enter neither priority_main nor main_union."
        ),
    )
    benchmark_parser.add_argument(
        "--diversity", choices=("representative", "none"), default="representative"
    )
    benchmark_parser.add_argument(
        "--cohort-eligibility",
        choices=("complete_release_label_pure_effective_blocks",),
        default=None,
        help=(
            "explicit sensitivity policy required when any strict CR or NCR row "
            "shares a complete-release effective leakage block with another label; "
            "it excludes that entire block before applying the exact ladder formula"
        ),
    )
    benchmark_parser.add_argument(
        "--fractions",
        nargs=3,
        type=float,
        metavar=("TRAIN", "VALIDATION", "TEST"),
        default=(0.8, 0.1, 0.1),
    )
    benchmark_parser.add_argument(
        "--no-full-cr-diagnostic",
        action="store_true",
        help=(
            "omit the supplementary full-CR prediction view, which otherwise reports "
            "exact-ID and same-effective-block training overlap and is not the independent test"
        ),
    )
    benchmark_parser.add_argument("--output-directory", type=Path, required=True)
    benchmark_parser.add_argument("--stem", default="coremof_cr_ncr_benchmark")
    benchmark_parser.add_argument("--overwrite", action="store_true")
    benchmark_parser.add_argument("--verify-cifs", action="store_true")

    split_parser = subparsers.add_parser(
        "split",
        help="make a deterministic parent-aware dataset split",
        description=(
            "Make a deterministic split using separate explanatory-parent and "
            "leakage-block policies. Project-defined terms are expanded in the "
            "JSON receipt."
        ),
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
        choices=SELECTABLE_PARENT_METHODS,
        default="priority_main",
        help=(
            "explanatory parent relation. Every non-control key reads its validated "
            "status/group/size triad: MATCHED means available size>=2, UNMATCHED "
            "means available size=1, and NOT_AVAILABLE supplies no edge; nulls never "
            "match. priority_main is the project-defined, conflict-aware full-release "
            "hierarchy over rac, mofid2, and mofid1 triads: RAC5 then MOFid v2 then "
            "MOFid v1. A lower "
            "group may attach unresolved rows to at most one stronger component "
            "but never merges multiple stronger components. It is not a row-wise "
            "first-nonmissing fallback; missing evidence becomes a per-structure "
            "singleton or is explicitly excluded, and Zeo++, topology, source ID, "
            "common name, CIF hash, provisional source-ID/MOFid transitive groups, "
            "and StructureMatcher "
            "are excluded. Here priority means parent-evidence precedence; it does "
            "not rank, schedule, or recalculate failed scientific features. rac5 "
            "uses all 264 finite binary64 values after -0.0 to "
            "+0.0 and float.hex with rtol=atol=0. mofid_v2 and mofid_v1 compare the "
            "complete canonicalized strings. rac5_zeo combines RAC5 with zeo; zeo "
            "uses exact binary64 float.hex equality for 13 intensive fields, N2/He "
            "radii 1.655/1.32 A, equal N2 channel dimension, and available equal "
            "framework periodic dimension. source_id uses the exact namespaced "
            "(source_database, source_id) pair; common_name uses exact canonicalized "
            "text; none makes one independent singleton per structure. "
            "identity_union selects the separate project-defined provisional "
            "source-ID/MOFid transitive groups read from "
            "identity_status/group/size: each named release is freshly recomputed over "
            "all its current rows by transitively unioning exact database-namespaced "
            "source ID plus eligible complete MOFid-v2 or MOFid-v1 text-equality edges "
            "with no precedence and without importing an earlier component. Each group "
            "and identity_size count one transitive "
            "connected component of structures, not edges or identifiers. Missing "
            "identifiers and unsuccessful MOFid statuses add no edge. It is not proof of "
            "structural identity and does not enter main_union. Here canonicalized "
            "text means convert to text, collapse Unicode whitespace, trim, reject "
            "case-insensitive empty, -, nan, none, null, n/a, na, unknown, missing, "
            "timeout, timed out, error, failed, fail, fail process, failed process, "
            "or process failed whole fields, then apply Unicode NFKC and casefold; it does not "
            "alter structures or use fuzzy matching. Eligible MOFid statuses are "
            "SUCCESS, SUCCESS_TOPOLOGY_UNKNOWN, SUCCESS_TOPOLOGY_ERROR, and "
            "SUCCESS_TOPOLOGY_TIMEOUT; unresolved reconciliation, ambiguous node, "
            "timeout, error, no-MOF, unmatched-node, and decomposition-error values "
            "never match through a shared null. "
            "rac5_crystalnets (RT-) reads rac_crystalnets_status/group/size and "
            "requires exact equality of all 264 finite RAC5 values plus the complete "
            "successful current CrystalNets fingerprint; "
            "mofid_v2_crystalnets (M2T-) reads mofid2_crystalnets_status/group/size. "
            "They require complete exact RAC5/MOFid plus successful "
            "release-authorized current CrystalNets fingerprints over "
            "network/count/net/agreement and every subnet node "
            "status/dimension/key/name/genome field; counts equal subnet count, "
            "heterogeneous top summaries and node name/genome may be null, and "
            "runtime/path/diagnostics are excluded. Missing, nonfinite, partial, "
            "timed-out, failed, or otherwise incomplete input adds no edge. M2T "
            "remains provisional whenever its MOFid-v2 input is provisional; if "
            "release-authorized MOFid-v2 values change, rebuild M2T before use. Its "
            "eligible MOFid-v2 statuses are exactly SUCCESS, "
            "SUCCESS_TOPOLOGY_UNKNOWN, SUCCESS_TOPOLOGY_ERROR, and "
            "SUCCESS_TOPOLOGY_TIMEOUT; the latter two are successful calculated "
            "identifiers with embedded topology qualifiers, and every other status "
            "adds no edge. "
            "structure_matcher_strict (SM-) is a "
            "Python 3.9/pymatgen 2024.2.8/NumPy 1.26.4 connected-component view of "
            "pairs exhaustively blocked only by parsed "
            "ElementComparator fractional-composition hash, with "
            "bidirectional pymatgen ElementComparator fits with ltol=stol=0.001, "
            "angle_tol=0.01, primitive_cell/attempt_supercell true, scale/allow_subset "
            "false, supercell_size=num_sites, no ignored species, and symmetric=true. "
            "Its pinned parser expands symmetry, uses site/frac tolerances 1e-4, "
            "checks occupancy, sorts, and preserves disorder without manual repair, "
            "occupancy selection, atom deletion, or chemistry editing; parser, timeout, OOM, "
            "matcher, or asymmetric cases are NOT_AVAILABLE rather than unmatched. "
            "Reported directional displacement divided by (V/Nsites)^(1/3) is "
            "dimensionless, not angstrom RMSD. SM components are convenience views, "
            "not duplicate proof or all-pairs claims; direct edges are authoritative. "
            "RT/M2T/SM prefixes use a criterion-bound "
            "length-delimited UTF-8 SHA-256 key, at least eight uppercase hex "
            "characters, extended only on collision. RT, M2T, and SM remain outside "
            "priority_main and main_union. Every receipt stores the selected "
            "method's exact fields, "
            "algorithm, missing/conflict behavior, exclusions, and policy relation "
            "(default: priority_main)"
        ),
    )
    split_parser.add_argument(
        "--leakage-guard",
        choices=LEAKAGE_GUARD_CHOICES,
        default="auto",
        help=(
            "split-block policy. main_union is the project-defined leakage guard, not "
            "an explanatory parent relation or proof of structural identity. It reads "
            "full manifest SHA-256 plus parent-table source/RAC5/MOFid groups and forms the "
            "transitive union over the complete release from those five exact "
            "relations before "
            "filters; missing CIF hashes fail, missing optional evidence adds no edge, "
            "and all edges have equal union status. Loaded releases require the source "
            "triad; only a manual compatibility object with no source fields falls back "
            "to database strip/upper plus source-ID strip/casefold/whitespace collapse, "
            "without NFKC or nonblank-placeholder rejection. parent_only is the project-defined "
            "guard that uses only the selected "
            "explanatory groups and adds no cross-method edge; an unresolved "
            "priority_main conflict is excluded with parent_only but may remain "
            "assigned and diagnosed with main_union. auto is only a project-defined "
            "selector: it chooses main_union for priority_main and parent_only for "
            "every other parent method and reads no scientific evidence itself "
            "(default: auto)"
        ),
    )
    split_parser.add_argument(
        "--missing-parent",
        choices=("singleton", "exclude"),
        default="singleton",
        help=(
            "handling of unavailable selected parent evidence: singleton creates "
            "SINGLETON:<structure_id>, so nulls never match; exclude records "
            "MISSING_PARENT_EVIDENCE and assigns no partition (default: singleton)"
        ),
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
        "--target-config",
        type=Path,
        help="join target files described by this JSON config before splitting",
    )
    split_parser.add_argument(
        "--require-target",
        action="append",
        help="require a non-null target before assignment; repeat for multiple targets",
    )
    split_parser.add_argument(
        "--required-target-mode",
        choices=("all", "any"),
        default="all",
        help="require all or any named targets (default: all)",
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
    if args.command == "merge-targets":
        return merge_target_data(args)
    if args.command == "benchmark-cr-ncr":
        return benchmark_cr_ncr(args)
    if args.command == "attach-targets":
        return attach_target_data(args)
    return 2
