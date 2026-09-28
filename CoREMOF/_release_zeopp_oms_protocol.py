#!/usr/bin/env python3
"""Run one provenance-bound Zeo++ open-metal-site calculation.

The result is diagnostic metadata.  A zero or nonzero OMS count does not by
itself assign a CR/NCR label, and the detector does not establish experimental
accessibility or catalytic activity.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import stat
import subprocess
import tempfile
import time
import traceback
import uuid
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence


PROTOCOL_ID = "coremof-zeopp-open-metal-sites-0.4.7-v1"
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


class ZeoppOMSError(RuntimeError):
    """An OMS input, output, or provenance contract was violated."""


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
        raise ZeoppOMSError("manifest SHA-256 mismatch")
    selected: Optional[Mapping[str, str]] = None
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != MANIFEST_FIELDS:
            raise ZeoppOMSError("manifest header is not exact")
        for index, row in enumerate(reader):
            if row.get("manifest_schema_version") != MANIFEST_SCHEMA_VERSION:
                raise ZeoppOMSError("manifest schema version is not exact")
            if row.get("row_index") != str(index):
                raise ZeoppOMSError("manifest row indices are not contiguous")
            if index == row_index:
                selected = dict(row)
    if selected is None:
        raise ZeoppOMSError("row index is outside manifest")
    return selected


def read_bound(path: Path, expected_size: int, expected_hash: str) -> bytes:
    before = path.lstat()
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode):
        raise ZeoppOMSError("CIF is not a regular non-symlink file")
    payload = path.read_bytes()
    after = path.lstat()
    identity_before = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
    identity_after = (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
    if identity_before != identity_after:
        raise ZeoppOMSError("CIF changed while being read")
    if len(payload) != expected_size or hashlib.sha256(payload).hexdigest() != expected_hash:
        raise ZeoppOMSError("CIF size/SHA-256 binding failed")
    return payload


def parse_oms(text: str) -> Mapping[str, Any]:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if len(lines) != 1:
        raise ZeoppOMSError("Zeo++ OMS output is not exactly one nonempty line")
    match = re.fullmatch(r"(?P<prefix>.*?)\s+#OMS=\s*(?P<count>\d+)\s*", lines[0])
    if not match:
        raise ZeoppOMSError("malformed Zeo++ OMS output")
    count = int(match.group("count"))
    return {
        "open_metal_site_count": count,
        "has_open_metal_sites": count > 0,
        "raw_line": lines[0],
    }


def parse_surface_definition(stdout: str) -> Optional[float]:
    values = [
        float(value)
        for value in re.findall(r"\bSurface definition\s*=\s*([-+0-9.eE]+)", stdout)
    ]
    if not values:
        return None
    if any(not math.isfinite(value) or value < 0 for value in values):
        raise ZeoppOMSError("invalid OMS surface definition in Zeo++ stdout")
    if any(value != values[0] for value in values[1:]):
        raise ZeoppOMSError("inconsistent OMS surface definitions in Zeo++ stdout")
    return values[0]


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
        raise ZeoppOMSError("existing record does not match this run")

    runner_path = Path(__file__).resolve()
    if sha256_file(runner_path) != args.runner_sha256:
        raise ZeoppOMSError("runner SHA-256 mismatch")
    if sha256_file(args.network) != args.network_sha256:
        raise ZeoppOMSError("Zeo++ binary SHA-256 mismatch")
    if sha256_file(args.contract) != args.contract_sha256:
        raise ZeoppOMSError("OMS metadata contract SHA-256 mismatch")
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    if contract.get("contract_id") != PROTOCOL_ID:
        raise ZeoppOMSError("OMS metadata contract ID mismatch")

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
        "contract": {"path": str(args.contract.resolve()), "sha256": args.contract_sha256},
        "configuration": {
            "zeopp_operation": "oms",
            "command_flags": ["-oms"],
            "probe_dependent": False,
            "atomic_radii": "Zeo++ built-in default CCDC radii",
            "requires_ccdc_licence": False,
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
            "hostname": os.uname().nodename,
        },
        "command": None,
        "open_metal_site_props": None,
        "execution_status": "RUNNING",
        "error_type": None,
        "error_message": None,
        "traceback": None,
    }
    try:
        payload = read_bound(cif_path, cif_size, cif_hash)
        args.private_root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix="zeopp-oms-{}-".format(structure_id), dir=str(args.private_root)
        ) as raw_tmp:
            tmp = Path(raw_tmp)
            private_cif = tmp / (structure_id + ".cif")
            private_cif.write_bytes(payload)
            output = private_cif.with_suffix(".oms")
            command = [str(args.network), "-oms", str(private_cif)]
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
                raise ZeoppOMSError("Zeo++ did not create OMS output")
            props = dict(parse_oms(output.read_text(encoding="utf-8", errors="replace")))
            props["surface_definition_A"] = parse_surface_definition(completed.stdout)
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
            record["open_metal_site_props"] = props
        record["execution_status"] = "SUCCESS"
    except BaseException as exc:
        record["execution_status"] = "ERROR"
        record["open_metal_site_props"] = None
        record["error_type"] = type(exc).__name__
        record["error_message"] = str(exc)[:2000]
        record["traceback"] = traceback.format_exc()[-8000:]
    record["elapsed_seconds"] = time.time() - started
    atomic_json(record_path, record)
    return {
        "structure_id": structure_id,
        "execution_status": record["execution_status"],
        "resumed": False,
    }


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
    parser.add_argument("--contract", required=True, type=Path)
    parser.add_argument("--contract-sha256", required=True)
    parser.add_argument("--timeout-seconds", type=int, default=300)
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    try:
        result = execute(parse_args(argv))
        print(json.dumps(result, sort_keys=True))
        return 0
    except BaseException as exc:
        print(
            json.dumps(
                {
                    "execution_status": "FATAL",
                    "error_type": type(exc).__name__,
                    "error_message": str(exc)[:2000],
                },
                sort_keys=True,
            ),
            file=os.sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
