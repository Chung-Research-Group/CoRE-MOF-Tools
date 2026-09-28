#!/usr/bin/env python3
"""Build and export a target-complete CR/NCR benchmark using public APIs.

This executable example does not certify model inputs or reproduce a frozen
paper dataset. It builds new exploratory assignments in a new output directory.
The workflow receipt is written last; an interrupted directory without it is
incomplete and must not be treated as a successful dataset.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
from collections.abc import Mapping

_ROOT = Path(__file__).resolve().parents[1]
if (_ROOT / "CoREMOF").is_dir() and str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from CoREMOF import __version__
from CoREMOF.dataset import CoREMOFDataset
from CoREMOF.targets import merge_targets_from_config


def _plain(value):
    if isinstance(value, Mapping):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(item) for item in value]
    return value


def canonical_digest(value):
    """SHA-256 of UTF-8 sorted, compact JSON, without ASCII or NaN coercion."""
    data = json.dumps(_plain(value), sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False, allow_nan=False).encode("utf-8")
    return hashlib.sha256(data).hexdigest()


def finite_target_reason(value):
    """Return None for a finite native number, otherwise an exclusion reason."""
    if value is None:
        return "MISSING"
    if isinstance(value, bool):
        return "BOOLEAN"
    # Python integers are finite, including integers too large for binary64.
    if isinstance(value, int):
        return None
    if isinstance(value, float):
        return None if math.isfinite(value) else "NONFINITE"
    return "NON_NUMERIC"


def target_eligibility(dataset, targets, required_targets):
    """Return sorted finite-target IDs and full-release exclusion accounting."""
    if not isinstance(required_targets, (list, tuple)):
        raise TypeError("required_targets must be an ordered sequence of exact names")
    required = tuple(required_targets)
    if (not required or any(not isinstance(name, str) or not name
                            or name != name.strip() for name in required)
            or len(set(required)) != len(required)):
        raise ValueError("required_targets must contain distinct nonblank names")
    unknown = set(required).difference(targets.target_columns)
    if unknown:
        raise ValueError("unknown required targets: " + ", ".join(sorted(unknown)))
    if set(dataset.structure_ids) != set(targets.structure_ids):
        raise ValueError("target table and release must have the same ID universe")
    eligible, rows = [], []
    for sid in sorted(dataset.structure_ids):
        values = targets.target_values(sid)
        reasons = {name: reason for name in required
                   if (reason := finite_target_reason(values[name])) is not None}
        if not reasons:
            eligible.append(sid)
        rows.append({"structure_id": sid, "target_complete": not reasons,
                     "excluded_targets": reasons})
    return tuple(eligible), tuple(rows)


def _file_binding(path, root):
    data = path.read_bytes()
    return {"path": path.relative_to(root).as_posix(), "size_bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest()}


def build_target_first_benchmark(
    release_path, target_config, output_directory, *, required_targets,
    group_criteria=("RT", "M2T"), ncr_pool_fractions=(0, 0.5, 1),
    seeds=(912, 913, 914, 915), fractions=(0.8, 0.1, 0.1),
    diversity="representative", verify_cif_files=False,
):
    """Load, merge, select finite eligibility, build, attach, and receipt outputs.

    Grouping and label purity use the entire validated release inside the
    benchmark API. Only target availability determines eligibility here; model
    graph/grid/descriptor certification remains a separate workflow.
    """
    output = Path(output_directory)
    if output.exists() or output.is_symlink():
        raise FileExistsError("output directory must be new: {}".format(output))
    dataset = CoREMOFDataset.from_release(release_path, verify_cif_files=verify_cif_files)
    targets = merge_targets_from_config(dataset, target_config)
    eligible, eligibility_rows = target_eligibility(dataset, targets, required_targets)
    required = tuple(required_targets)
    definitions = {name: _plain(targets.target_definitions[name]) for name in required}
    for name, definition in definitions.items():
        if definition["unit"] is None or definition["conditions"] is None:
            raise ValueError("required target {} must declare units and conditions".format(name))
        if definition["value_type"] not in (None, "float", "int"):
            raise ValueError("required target {} must use numeric or native JSON values".format(name))
    if len(fractions) != 3:
        raise ValueError("fractions must contain train, validation, and test")
    requested = {
        "checkers": "5checker", "required_targets": list(required),
        "group_criteria": [group_criteria] if isinstance(group_criteria, str) else list(group_criteria),
        "ncr_pool_fractions": list(ncr_pool_fractions), "seeds": list(seeds),
        "total_size": "full_cr", "train": fractions[0], "val": fractions[1], "test": fractions[2],
        "cohort_eligibility": "complete_release_label_pure_effective_blocks",
        "diversity": diversity, "test_policy": "fixed_pure_cr",
        "include_full_cr_diagnostic": False, "partition_strategy": "transition_balanced",
        "verify_cif_files": verify_cif_files,
    }
    suite = dataset.classify("5checker").build_cr_ncr_benchmark(
        eligible_structure_ids=eligible, group_criteria=group_criteria,
        ncr_pool_fractions=ncr_pool_fractions, seeds=seeds,
        train=fractions[0], val=fractions[1], test=fractions[2],
        cohort_eligibility=requested["cohort_eligibility"], diversity=diversity,
        partition_strategy="transition_balanced", include_full_cr_diagnostic=False,
    )
    suite_receipt = suite.receipt()
    cohort_receipt = suite_receipt["cohort_receipt"]
    eligibility_digest = canonical_digest(eligible)
    if cohort_receipt["eligibility_filter"]["eligible_ids_sha256"] != eligibility_digest:
        raise RuntimeError("benchmark eligibility digest differs from the workflow")
    before = canonical_digest(suite_receipt)
    # Extra configured endpoints remain optional; every requested endpoint was
    # checked above and is rechecked in the attached views below.
    missing_policy = "error" if set(required) == set(targets.target_columns) else "keep"
    attached = suite.attach_targets(targets, missing=missing_policy)
    if (attached.original_assignment_digest != suite.assignment_digest
            or canonical_digest(suite.receipt()) != before):
        raise RuntimeError("target attachment changed the frozen assignments")
    for run in suite.runs:
        view = attached.run_views[run.run_key]
        if dict(view.assignments) != dict(run.assignments):
            raise RuntimeError("target attachment changed a run partition")
        if any(finite_target_reason(view.values_by_id[sid][name]) is not None
               for sid in view.structure_ids for name in required):
            raise RuntimeError("an assigned structure lacks a required finite target")

    # Exclusive creation prevents reruns from overwriting a frozen comparison.
    # Public writers provide their own transactional output generations.
    output.mkdir(parents=True, exist_ok=False)
    targets.write(output / "targets", stem="merged_targets")
    with (output / "target_eligibility.csv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("structure_id", "target_complete", "excluded_targets"))
        writer.writeheader()
        for row in eligibility_rows:
            writer.writerow({**row, "target_complete": str(row["target_complete"]).lower(),
                             "excluded_targets": json.dumps(row["excluded_targets"], sort_keys=True)})
    suite_path = suite.write(output / "splits")
    attached_path = attached.write(output / "model_inputs")
    artifacts = [_file_binding(path, output) for path in sorted(output.rglob("*")) if path.is_file()]
    receipt = {
        "schema_version": "coremof-target-first-workflow-receipt/1.0",
        "status": "PASS", "package_version": __version__,
        "example_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "model_input_certification_performed": False,
        "official_split": False, "publication_authorized": False,
        "requested_configuration": requested,
        "required_endpoint_definitions": definitions,
        "required_endpoint_definitions_sha256": canonical_digest(definitions),
        "target_merge_receipt": targets.receipt(),
        "target_merge_receipt_sha256": canonical_digest(targets.receipt()),
        "release_binding": cohort_receipt["release_binding"],
        "eligibility": {
            "rule": "all required native numeric values finite; null, bool and nonnumeric excluded; zero retained",
            "required_targets": list(required), "eligible_ids": list(eligible),
            "eligible_ids_sha256": eligibility_digest, "eligible_count": len(eligible),
            "release_structure_count": len(dataset.structure_ids),
            "manifest": "target_eligibility.csv", "target_magnitudes_used_for_assignment": False,
        },
        "grouping": {
            "group_criteria": cohort_receipt["group_criteria"],
            "group_criterion_definitions": cohort_receipt["group_criterion_definitions"],
            "criterion_groups_sha256": cohort_receipt["criterion_groups_sha256"],
            "effective_leakage_blocks": cohort_receipt["effective_leakage_blocks"],
            "leakage_guard": cohort_receipt["leakage_guard"],
            "leakage_guard_definition": cohort_receipt["effective_leakage_policy_definition"],
            "label_purity": cohort_receipt["cohort_eligibility_policy"],
            "grouping_and_label_purity_universe": "complete release before target eligibility",
        },
        "suite_assignment_sha256": suite.assignment_digest,
        "suite_receipt": _file_binding(suite_path / "receipt.json", output),
        "target_attachment_assignment_sha256": attached.original_assignment_digest,
        "target_attachment_receipt": _file_binding(attached_path / "receipt.json", output),
        "target_attachment_missing_policy": missing_policy,
        "all_required_targets_finite_after_attachment": True,
        "original_assignments_unchanged": True,
        "artifacts": artifacts,
    }
    receipt_bytes = (json.dumps(_plain(receipt), sort_keys=True, indent=2,
                               ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")
    receipt_binding = {"path": "workflow_receipt.json",
                       "sha256": hashlib.sha256(receipt_bytes).hexdigest()}
    with (output / "SHA256SUMS").open("x", encoding="utf-8") as handle:
        for item in artifacts + [receipt_binding]:
            handle.write("{}  {}\n".format(item["sha256"], item["path"]))
    with (output / "workflow_receipt.json").open("xb") as handle:
        handle.write(receipt_bytes)
    return output / "workflow_receipt.json"


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog=("RT matches all 264 finite RAC5 descriptors plus the complete successful "
                "CrystalNets fingerprint; M2T matches the complete canonical MOFid-v2 "
                "string plus that fingerprint. Both are release-authorized reference "
                "criteria, not proof of chemical identity. Their groups are combined "
                "with the full-release CIF-hash/source/RAC5/MOFid leakage guard before "
                "eligibility. The common test contains only strictly five-checker CR "
                "structures, and assignments remain exploratory."),
    )
    parser.add_argument("release", type=Path, help="extracted validated local release")
    parser.add_argument("--target-config", type=Path, required=True)
    parser.add_argument("--require-target", action="append", required=True,
                        help="exact required endpoint name; repeat for every endpoint")
    parser.add_argument("--output-directory", type=Path, required=True, help="new directory only")
    parser.add_argument("--group-criteria", nargs="+", default=["RT", "M2T"])
    parser.add_argument("--ncr-pool-fractions", nargs="+", type=float, default=[0, 0.5, 1])
    parser.add_argument("--seeds", nargs="+", type=int, default=[912, 913, 914, 915])
    parser.add_argument("--fractions", nargs=3, type=float, default=[0.8, 0.1, 0.1])
    parser.add_argument("--diversity", choices=("representative", "none"), default="representative",
                        help="none is an explicit non-representative sensitivity/test profile")
    parser.add_argument("--verify-cifs", action="store_true")
    args = parser.parse_args(argv)
    try:
        receipt = build_target_first_benchmark(
            args.release, args.target_config, args.output_directory,
            required_targets=args.require_target, group_criteria=args.group_criteria,
            ncr_pool_fractions=args.ncr_pool_fractions, seeds=args.seeds,
            fractions=args.fractions, diversity=args.diversity, verify_cif_files=args.verify_cifs,
        )
    except (OSError, TypeError, ValueError) as error:
        parser.exit(2, "error: {}\n".format(error))
    print(receipt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
