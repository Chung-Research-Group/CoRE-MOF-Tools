#!/usr/bin/env python3
"""Run one hash-bound Zeo++ bonded-framework dimensionality calculation.

This is deliberately separate from the probe-dependent void-channel
dimensionality emitted by ``run_zeopp17_candidate_one.py``.  It runs Zeo++
``-ha -strinfo`` and records the number of 1D, 2D, and 3D bonded frameworks
and discrete molecular components.  Results are diagnostic metadata; this
runner does not convert any result directly into a CR/NCR label.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import stat
import subprocess
import tempfile
import time
import traceback
import uuid
from collections import Counter
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence


PROTOCOL_ID = "coremof-zeopp-framework-dimension-0.4.7-v1"
RECORD_SCHEMA_VERSION = "1.0"
MANIFEST_SCHEMA_VERSION = "1.0"
MANIFEST_FIELDS = (
    "manifest_schema_version",
    "row_index",
    "structure_id",
    "original_source_label",
    "canonical_cif_version",
    "canonical_cif_path",
    "cif_basename",
    "cif_size_bytes",
    "cif_sha256",
    "v10_fallback_path",
    "source_refcode_raw",
    "normalized_refcode",
    "source_family",
    "parent_group",
    "parent_group_full_size",
    "target_manifest_path",
    "target_manifest_sha256",
)


class ZeoppFrameworkDimensionError(RuntimeError):
    """A framework-dimensionality input or output contract was violated."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(".{}.{}.tmp".format(path.name, uuid.uuid4().hex))
    try:
        with temporary.open("x", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(str(temporary), str(path))
    finally:
        if temporary.exists():
            temporary.unlink()


def load_row(manifest: Path, expected_hash: str, row_index: int) -> Mapping[str, str]:
    if sha256_file(manifest) != expected_hash:
        raise ZeoppFrameworkDimensionError("manifest SHA-256 mismatch")
    selected: Optional[Mapping[str, str]] = None
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != MANIFEST_FIELDS:
            raise ZeoppFrameworkDimensionError("manifest header is not exact")
        for index, row in enumerate(reader):
            if row.get("manifest_schema_version") != MANIFEST_SCHEMA_VERSION:
                raise ZeoppFrameworkDimensionError("manifest schema version is not exact")
            if row.get("row_index") != str(index):
                raise ZeoppFrameworkDimensionError("manifest row indices are not contiguous")
            if index == row_index:
                selected = dict(row)
    if selected is None:
        raise ZeoppFrameworkDimensionError("row index is outside manifest")
    return selected


def read_bound(path: Path, expected_size: int, expected_hash: str) -> bytes:
    before = path.lstat()
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode):
        raise ZeoppFrameworkDimensionError("CIF is not a regular non-symlink file")
    payload = path.read_bytes()
    after = path.lstat()
    identity_before = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
    identity_after = (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
    if identity_before != identity_after:
        raise ZeoppFrameworkDimensionError("CIF changed while being read")
    if len(payload) != expected_size or hashlib.sha256(payload).hexdigest() != expected_hash:
        raise ZeoppFrameworkDimensionError("CIF size/SHA-256 binding failed")
    return payload


def parse_strinfo(text: str) -> Mapping[str, Any]:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if len(lines) != 1:
        raise ZeoppFrameworkDimensionError("Zeo++ strinfo output is not exactly one nonempty line")
    line = lines[0]
    pattern = re.compile(
        r"^(?P<prefix>.*?)\s+"
        r"(?P<segments>\d+)\s+segments:\s+"
        r"(?P<frameworks>\d+)\s+framework\(s\)\s+"
        r"\(1D/2D/3D\s+(?P<d1>\d+)\s+(?P<d2>\d+)\s+(?P<d3>\d+)\s*\)\s+"
        r"and\s+(?P<molecules>\d+)\s+molecule\(s\)\.\s+"
        r"Identified dimensionality of framework\(s\):\s*(?P<dimensions>[123\s]*)$"
    )
    match = pattern.match(line)
    if not match:
        raise ZeoppFrameworkDimensionError("malformed Zeo++ strinfo output")
    segments = int(match.group("segments"))
    frameworks = int(match.group("frameworks"))
    molecules = int(match.group("molecules"))
    counts = {"1": int(match.group("d1")), "2": int(match.group("d2")), "3": int(match.group("d3"))}
    dimensions = [int(item) for item in match.group("dimensions").split()]
    if frameworks != sum(counts.values()):
        raise ZeoppFrameworkDimensionError("framework total differs from 1D/2D/3D counts")
    if segments != frameworks + molecules:
        raise ZeoppFrameworkDimensionError("segment total differs from framework plus molecule counts")
    if Counter(dimensions) != Counter({int(key): value for key, value in counts.items() if value}):
        raise ZeoppFrameworkDimensionError("listed framework dimensions differ from dimension counts")
    return {
        "segment_count": segments,
        "framework_count": frameworks,
        "framework_1d_count": counts["1"],
        "framework_2d_count": counts["2"],
        "framework_3d_count": counts["3"],
        "molecule_count": molecules,
        "framework_dimensions": dimensions,
        "maximum_framework_dimension": max(dimensions) if dimensions else 0,
        "has_periodic_framework": frameworks > 0,
        "all_components_are_discrete_molecules": frameworks == 0,
        "raw_line": line,
    }


def execute(args: argparse.Namespace) -> Mapping[str, Any]:
    started = time.time()
    row = load_row(args.manifest, args.manifest_sha256, args.row_index)
    structure_id = row["structure_id"]
    record_path = args.output_root / "records" / (structure_id + ".json")
    if record_path.exists():
        existing = json.loads(record_path.read_text(encoding="utf-8"))
        if (
            existing.get("protocol_id") == PROTOCOL_ID
            and existing.get("manifest", {}).get("sha256") == args.manifest_sha256
            and existing.get("manifest_row_index") == args.row_index
            and existing.get("input", {}).get("cif_sha256") == row["cif_sha256"]
        ):
            return {
                "structure_id": structure_id,
                "execution_status": existing.get("execution_status"),
                "resumed": True,
            }
        raise ZeoppFrameworkDimensionError("existing record does not match this run")

    runner_path = Path(__file__).resolve()
    if sha256_file(runner_path) != args.runner_sha256:
        raise ZeoppFrameworkDimensionError("runner SHA-256 mismatch")
    if sha256_file(args.network) != args.network_sha256:
        raise ZeoppFrameworkDimensionError("Zeo++ binary SHA-256 mismatch")

    cif_path = Path(row["canonical_cif_path"])
    cif_size = int(row["cif_size_bytes"])
    cif_hash = row["cif_sha256"]
    record: dict[str, Any] = {
        "record_schema_version": RECORD_SCHEMA_VERSION,
        "protocol_id": PROTOCOL_ID,
        "candidate_metadata": True,
        "curation_role": "DIAGNOSTIC_ONLY_NO_AUTOMATIC_CR_NCR_EXCLUSION",
        "manifest": {"path": str(args.manifest.resolve()), "sha256": args.manifest_sha256},
        "manifest_row_index": args.row_index,
        "structure_id": structure_id,
        "source_family": row["source_family"],
        "canonical_cif_version": row["canonical_cif_version"],
        "input": {"path": str(cif_path), "cif_size_bytes": cif_size, "cif_sha256": cif_hash},
        "configuration": {
            "zeopp_operation": "strinfo",
            "high_accuracy": True,
            "atomic_radii": "Zeo++ built-in default CCDC radii",
            "distinguish_from_channel_dimension": True,
            "timeout_seconds": args.timeout_seconds,
        },
        "provenance": {
            "runner_path": str(runner_path),
            "runner_sha256": args.runner_sha256,
            "network_path": str(args.network.resolve()),
            "network_sha256": args.network_sha256,
            "conda_package": "zeopp-lsmo=0.4.7=h27087fc_0",
            "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
            "slurm_array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID", ""),
            "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID", ""),
        },
        "command": None,
        "framework_dimension": None,
        "execution_status": "RUNNING",
        "error_type": None,
        "error_message": None,
        "traceback": None,
    }
    try:
        payload = read_bound(cif_path, cif_size, cif_hash)
        args.private_root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="zeopp-strdim-{}-".format(structure_id), dir=str(args.private_root)) as raw_tmp:
            tmp = Path(raw_tmp)
            private_cif = tmp / (structure_id + ".cif")
            private_cif.write_bytes(payload)
            output = private_cif.with_suffix(".strinfo")
            command = [str(args.network), "-ha", "-strinfo", str(private_cif)]
            command_started = time.time()
            completed = subprocess.run(
                command,
                cwd=str(tmp),
                capture_output=True,
                text=True,
                timeout=args.timeout_seconds,
                check=True,
            )
            if not output.is_file():
                raise ZeoppFrameworkDimensionError("Zeo++ did not create strinfo output")
            record["command"] = {
                "argv": command,
                "elapsed_seconds": time.time() - command_started,
                "returncode": completed.returncode,
                "stdout_size_bytes": len(completed.stdout.encode("utf-8")),
                "stdout_tail": completed.stdout[-4000:],
                "stderr_size_bytes": len(completed.stderr.encode("utf-8")),
                "stderr_tail": completed.stderr[-4000:],
                "output_sha256": sha256_file(output),
            }
            record["framework_dimension"] = parse_strinfo(output.read_text(encoding="utf-8", errors="replace"))
        record["execution_status"] = "SUCCESS"
    except BaseException as exc:
        record["execution_status"] = "ERROR"
        record["framework_dimension"] = None
        record["error_type"] = type(exc).__name__
        record["error_message"] = str(exc)[:2000]
        record["traceback"] = traceback.format_exc()[-8000:]
    record["elapsed_seconds"] = time.time() - started
    atomic_json(record_path, record)
    return {"structure_id": structure_id, "execution_status": record["execution_status"], "resumed": False}


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--row-index", required=True, type=int)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--private-root", required=True, type=Path)
    parser.add_argument("--network", required=True, type=Path)
    parser.add_argument("--network-sha256", required=True)
    parser.add_argument("--runner-sha256", required=True)
    parser.add_argument("--timeout-seconds", type=int, default=600)
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    try:
        result = execute(parse_args(argv))
        print(json.dumps(result, sort_keys=True))
        return 0
    except BaseException as exc:
        print(
            json.dumps(
                {"execution_status": "FATAL", "error_type": type(exc).__name__, "error_message": str(exc)[:2000]},
                sort_keys=True,
            ),
            file=os.sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
