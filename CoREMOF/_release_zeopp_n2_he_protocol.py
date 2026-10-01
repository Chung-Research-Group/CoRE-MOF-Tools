#!/usr/bin/env python3
"""Calculate one hash-bound Zeo++ N2/He metadata record.

The protocol keeps probe-independent pore diameters separate from probe-
dependent properties.  Nitrogen uses 1.655 A as both channel and probe radius
for surface area and probe-occupiable volume, and as the probe radius for
channel dimensionality.  Helium uses 1.32 A and contributes only the
probe-occupiable accessible void fraction requested for metadata.

Legacy 0 A results are deliberately not overwritten or copied into this
protocol.  They remain a separate ``zero_probe_props`` evidence layer.
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


PROTOCOL_ID = "coremof-zeopp-n2-1p655-he-1p32-0.4.7-v2.1"
RECORD_SCHEMA_VERSION = "2.0"
N2_RADIUS_A = 1.655
HE_RADIUS_A = 1.32
SURFACE_SAMPLES_PER_ATOM = 5000
VOLUME_SAMPLES_TOTAL = 5000

INTRINSIC_FEATURES = (
    "LCD_A",
    "PLD_A",
    "LFPD_A",
    "density_g_cm3",
)
N2_FEATURES = (
    "ASA_A2",
    "ASA_m2_cm3",
    "ASA_m2_g",
    "NASA_A2",
    "NASA_m2_cm3",
    "NASA_m2_g",
    "AV_A3",
    "AV_cm3_g",
    "AV_VF",
    "NAV_A3",
    "NAV_cm3_g",
    "NAV_VF",
    "channel_dimension",
)
HE_FEATURES = ("AV_VF",)

COMMANDS = {
    "intrinsic_pore_diameter": ["-ha", "-res"],
    "n2_surface_area": [
        "-ha", "-sa", str(N2_RADIUS_A), str(N2_RADIUS_A),
        str(SURFACE_SAMPLES_PER_ATOM),
    ],
    "n2_pore_volume": [
        "-ha", "-volpo", str(N2_RADIUS_A), str(N2_RADIUS_A),
        str(VOLUME_SAMPLES_TOTAL),
    ],
    "n2_channel_dimension": ["-ha", "-chan", str(N2_RADIUS_A)],
    "he_pore_volume": [
        "-ha", "-volpo", str(HE_RADIUS_A), str(HE_RADIUS_A),
        str(VOLUME_SAMPLES_TOTAL),
    ],
}


class ZeoppN2HeError(RuntimeError):
    """The N2/He Zeo++ input, output, or provenance contract failed."""


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
        raise ZeoppN2HeError("manifest SHA-256 mismatch")
    selected: Optional[Mapping[str, str]] = None
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        for index, row in enumerate(csv.DictReader(handle)):
            if row.get("row_index") != str(index):
                raise ZeoppN2HeError("manifest row indices are not contiguous")
            if index == row_index:
                selected = dict(row)
    if selected is None:
        raise ZeoppN2HeError("row index is outside manifest")
    return selected


def read_bound(path: Path, expected_size: int, expected_hash: str) -> bytes:
    before = path.lstat()
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode):
        raise ZeoppN2HeError("CIF is not a regular non-symlink file")
    payload = path.read_bytes()
    after = path.lstat()
    identity_before = (
        before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns,
    )
    identity_after = (
        after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns,
    )
    if identity_before != identity_after:
        raise ZeoppN2HeError("CIF changed while being read")
    if len(payload) != expected_size or hashlib.sha256(payload).hexdigest() != expected_hash:
        raise ZeoppN2HeError("CIF size/SHA-256 binding failed")
    return payload


def key_float(line: str, key: str) -> float:
    match = re.search(re.escape(key) + r":\s*([-+0-9.eE]+)", line)
    if not match:
        raise ZeoppN2HeError("missing Zeo++ output key {}".format(key))
    value = float(match.group(1))
    if not math.isfinite(value):
        raise ZeoppN2HeError("non-finite Zeo++ output key {}".format(key))
    return value


def run_command(
    binary: Path,
    arguments: list[str],
    output: Path,
    cif: Path,
    timeout: int,
) -> Mapping[str, Any]:
    command = [str(binary)] + arguments + [str(output), str(cif)]
    started = time.time()
    completed = subprocess.run(
        command,
        cwd=str(output.parent),
        capture_output=True,
        text=True,
        timeout=timeout,
        check=True,
    )
    if not output.is_file():
        raise ZeoppN2HeError("Zeo++ did not create {}".format(output.name))
    return {
        "argv": command,
        "elapsed_seconds": time.time() - started,
        "returncode": completed.returncode,
        "stdout_size_bytes": len(completed.stdout.encode("utf-8")),
        "stdout_tail": completed.stdout[-2000:],
        "stderr_size_bytes": len(completed.stderr.encode("utf-8")),
        "stderr_tail": completed.stderr[-2000:],
        "output_sha256": sha256_file(output),
    }


def validate_fraction(name: str, value: float) -> None:
    if value < 0.0 or value > 1.0:
        raise ZeoppN2HeError("{} is outside [0, 1]".format(name))


def parse_channel_topology(text: str) -> Mapping[str, Any]:
    """Parse every accessible channel dimension and return their maximum."""
    matches = []
    for line in text.splitlines():
        match = re.search(
            r"\b(\d+)\s+channels?\s+identified\s+of\s+dimensionality(?P<dimensions>(?:\s+[0-3])*)\s*$",
            line.strip(),
        )
        if match:
            dimensions = [int(value) for value in match.group("dimensions").split()]
            matches.append((int(match.group(1)), dimensions))
    if len(matches) != 1:
        raise ZeoppN2HeError("malformed N2 channel-dimensionality output")
    channel_count, dimensions = matches[0]
    if channel_count == 0:
        if dimensions:
            raise ZeoppN2HeError("zero N2 channels unexpectedly report dimensions")
        maximum = 0
    else:
        if len(dimensions) != channel_count or any(value not in (1, 2, 3) for value in dimensions):
            raise ZeoppN2HeError("N2 channel count/dimension list is inconsistent")
        maximum = max(dimensions)
    return {
        "channel_count": channel_count,
        "channel_dimensions": dimensions,
        "maximum_channel_dimension": maximum,
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
            and existing.get("manifest_row_index") == args.row_index
            and existing.get("manifest", {}).get("sha256") == args.manifest_sha256
            and existing.get("input", {}).get("cif_sha256") == row["cif_sha256"]
        ):
            return {
                "structure_id": structure_id,
                "execution_status": existing.get("execution_status"),
                "resumed": True,
            }
        raise ZeoppN2HeError("existing record does not match this run")

    runner_path = Path(__file__).resolve()
    runner_hash = sha256_file(runner_path)
    if runner_hash != args.runner_sha256:
        raise ZeoppN2HeError("runner SHA-256 mismatch")
    if sha256_file(args.network) != args.network_sha256:
        raise ZeoppN2HeError("Zeo++ binary SHA-256 mismatch")
    if sha256_file(args.namespace_policy) != args.namespace_policy_sha256:
        raise ZeoppN2HeError("probe namespace policy SHA-256 mismatch")

    cif_path = Path(row["canonical_cif_path"])
    cif_size = int(row["cif_size_bytes"])
    cif_hash = row["cif_sha256"]
    record: dict[str, Any] = {
        "record_schema_version": RECORD_SCHEMA_VERSION,
        "protocol_id": PROTOCOL_ID,
        "candidate_metadata": True,
        "historical_zeopp_equivalence_claimed": False,
        "parent_edge_authorized": False,
        "legacy_zero_probe_policy": {
            "namespace": "zero_probe_props",
            "status": "PRESERVED_AS_SEPARATE_LEGACY_EVIDENCE_NOT_RECOMPUTED_HERE",
            "legacy_protocol_id": "coremof-zeopp17-0.4.7-candidate-v1",
        },
        "namespace_policy": {
            "path": str(args.namespace_policy.resolve()),
            "sha256": args.namespace_policy_sha256,
        },
        "manifest": {
            "path": str(args.manifest.resolve()),
            "sha256": args.manifest_sha256,
        },
        "manifest_row_index": args.row_index,
        "structure_id": structure_id,
        "source_family": row["source_family"],
        "input": {
            "path": str(cif_path),
            "cif_size_bytes": cif_size,
            "cif_sha256": cif_hash,
        },
        "configuration": {
            "high_accuracy": True,
            "atomic_radii": "Zeo++ built-in default CCDC radii",
            "intrinsic_pore_diameter_probe_radius_A": None,
            "n2_channel_radius_A": N2_RADIUS_A,
            "n2_probe_radius_A": N2_RADIUS_A,
            "n2_surface_samples_per_atom": SURFACE_SAMPLES_PER_ATOM,
            "n2_volume_samples_total": VOLUME_SAMPLES_TOTAL,
            "n2_channel_dimension_probe_radius_A": N2_RADIUS_A,
            "he_channel_radius_A": HE_RADIUS_A,
            "he_probe_radius_A": HE_RADIUS_A,
            "he_volume_samples_total": VOLUME_SAMPLES_TOTAL,
            "he_output_scope": "POAV_Volume_fraction_ONLY",
            "monte_carlo_seed_control": "NOT_EXPOSED_BY_THIS_ZEOPP_CLI_PROTOCOL",
        },
        "provenance": {
            "runner_path": str(runner_path),
            "runner_sha256": runner_hash,
            "network_path": str(args.network.resolve()),
            "network_sha256": args.network_sha256,
            "conda_package": "zeopp-lsmo=0.4.7=h27087fc_0",
            "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
            "slurm_array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID", ""),
            "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID", ""),
            "runtime_node": os.environ.get("SLURMD_NODENAME", ""),
        },
        "feature_schema": {
            "intrinsic_props": list(INTRINSIC_FEATURES),
            "N2_probe_props": list(N2_FEATURES),
            "He_probe_props": list(HE_FEATURES),
        },
        "features": None,
        "N2_channel_topology": None,
        "commands": None,
        "execution_status": "RUNNING",
        "error_type": None,
        "error_message": None,
        "traceback": None,
    }

    try:
        payload = read_bound(cif_path, cif_size, cif_hash)
        args.private_root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix="zeopp-n2-he-{}-".format(structure_id),
            dir=str(args.private_root),
        ) as raw_tmp:
            tmp = Path(raw_tmp)
            private_cif = tmp / (structure_id + ".cif")
            private_cif.write_bytes(payload)
            paths = {
                "intrinsic_pore_diameter": tmp / "intrinsic.res",
                "n2_surface_area": tmp / "n2.sa",
                "n2_pore_volume": tmp / "n2.volpo",
                "n2_channel_dimension": tmp / "n2.chan",
                "he_pore_volume": tmp / "he.volpo",
            }
            command_records: dict[str, Any] = {}
            for name, arguments in COMMANDS.items():
                command_records[name] = run_command(
                    args.network,
                    arguments,
                    paths[name],
                    private_cif,
                    args.timeout_seconds,
                )
            record["commands"] = command_records

            res = paths["intrinsic_pore_diameter"].read_text(
                encoding="utf-8", errors="replace",
            ).split()
            if len(res) < 4:
                raise ZeoppN2HeError("malformed pore-diameter output")
            lcd, pld, lfpd = map(float, res[-3:])
            n2_sa_line = paths["n2_surface_area"].read_text(
                encoding="utf-8", errors="replace",
            ).splitlines()[0]
            n2_pv_line = paths["n2_pore_volume"].read_text(
                encoding="utf-8", errors="replace",
            ).splitlines()[0]
            he_pv_line = paths["he_pore_volume"].read_text(
                encoding="utf-8", errors="replace",
            ).splitlines()[0]
            n2_chan_line = paths["n2_channel_dimension"].read_text(
                encoding="utf-8", errors="replace",
            ).splitlines()[0]
            n2_channel_topology = parse_channel_topology(n2_chan_line)

            densities = {
                key_float(n2_sa_line, "Density"),
                key_float(n2_pv_line, "Density"),
                key_float(he_pv_line, "Density"),
            }
            if len(densities) != 1:
                raise ZeoppN2HeError("Zeo++ density differs between commands")
            density = next(iter(densities))
            n2_av_vf = key_float(n2_pv_line, "POAV_Volume_fraction")
            n2_nav_vf = key_float(n2_pv_line, "PONAV_Volume_fraction")
            he_av_vf = key_float(he_pv_line, "POAV_Volume_fraction")
            validate_fraction("N2 POAV volume fraction", n2_av_vf)
            validate_fraction("N2 PONAV volume fraction", n2_nav_vf)
            validate_fraction("He POAV volume fraction", he_av_vf)

            features = {
                "intrinsic_props": {
                    "LCD_A": lcd,
                    "PLD_A": pld,
                    "LFPD_A": lfpd,
                    "density_g_cm3": density,
                },
                "N2_probe_props": {
                    "ASA_A2": key_float(n2_sa_line, "ASA_A^2"),
                    "ASA_m2_cm3": key_float(n2_sa_line, "ASA_m^2/cm^3"),
                    "ASA_m2_g": key_float(n2_sa_line, "ASA_m^2/g"),
                    "NASA_A2": key_float(n2_sa_line, "NASA_A^2"),
                    "NASA_m2_cm3": key_float(n2_sa_line, "NASA_m^2/cm^3"),
                    "NASA_m2_g": key_float(n2_sa_line, "NASA_m^2/g"),
                    "AV_A3": key_float(n2_pv_line, "POAV_A^3"),
                    "AV_cm3_g": key_float(n2_pv_line, "POAV_cm^3/g"),
                    "AV_VF": n2_av_vf,
                    "NAV_A3": key_float(n2_pv_line, "PONAV_A^3"),
                    "NAV_cm3_g": key_float(n2_pv_line, "PONAV_cm^3/g"),
                    "NAV_VF": n2_nav_vf,
                    "channel_dimension": n2_channel_topology["maximum_channel_dimension"],
                },
                "He_probe_props": {"AV_VF": he_av_vf},
            }
            if tuple(features["intrinsic_props"]) != INTRINSIC_FEATURES:
                raise ZeoppN2HeError("intrinsic feature schema mismatch")
            if tuple(features["N2_probe_props"]) != N2_FEATURES:
                raise ZeoppN2HeError("N2 feature schema mismatch")
            if tuple(features["He_probe_props"]) != HE_FEATURES:
                raise ZeoppN2HeError("He feature schema mismatch")
            for namespace in features.values():
                if not all(math.isfinite(float(value)) for value in namespace.values()):
                    raise ZeoppN2HeError("non-finite Zeo++ feature")
            record["features"] = features
            record["N2_channel_topology"] = n2_channel_topology
        record["execution_status"] = "SUCCESS"
    except BaseException as exc:
        record["execution_status"] = "ERROR"
        record["features"] = None
        record["N2_channel_topology"] = None
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
    parser.add_argument("--namespace-policy", required=True, type=Path)
    parser.add_argument("--namespace-policy-sha256", required=True)
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
