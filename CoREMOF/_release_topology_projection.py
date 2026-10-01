#!/usr/bin/env python3
"""Audit CrystalNets shards and materialize compact/rich topology features."""

from __future__ import annotations

import argparse
import collections
import csv
import datetime as dt
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple


PUBLIC_COLUMNS = (
    "structure_id",
    "topology_available",
    "execution_status",
    "network_dimension",
    "interpenetrated_subnet_count",
    "catenation_degree",
    "single_node_net",
    "all_node_net",
    "single_all_agree",
)
CUSTOM_TOPOLOGY = re.compile(
    r"^(?P<key>[0-3]-[a-z]{12}) \((?P<genome>.*)\)$"
)


class AuditError(RuntimeError):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_topology(value: Mapping[str, Any]) -> Dict[str, Any]:
    raw = str(value.get("topology") or "")
    dimension = value.get("dimension")
    if raw.startswith("FAILED with:"):
        return {
            "status": "ERROR",
            "dimension": None,
            "topology_name": None,
            "topology_key": None,
            "topological_genome": None,
            "error_message": raw[len("FAILED with:") :].strip(),
        }
    if raw == "0-dimensional":
        return {
            "status": "SUCCESS",
            "dimension": 0,
            "topology_name": None,
            "topology_key": "0-dimensional",
            "topological_genome": None,
        }
    match = CUSTOM_TOPOLOGY.fullmatch(raw)
    if match:
        return {
            "status": "SUCCESS",
            "dimension": dimension,
            "topology_name": None,
            "topology_key": match.group("key"),
            "topological_genome": match.group("genome"),
        }
    return {
        "status": str(value.get("status") or "NOT_AVAILABLE"),
        "dimension": dimension,
        "topology_name": raw or None,
        "topology_key": raw or None,
        "topological_genome": None,
    }


def consensus(values: Sequence[Any]) -> Any:
    return values[0] if values and all(value == values[0] for value in values) else None


def normalized_record(raw: Mapping[str, Any]) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    status = str(raw.get("execution_status") or "")
    subnets = []
    if status == "SUCCESS":
        for expected_index, source in enumerate(raw.get("subnets") or [], start=1):
            if source.get("subnet_index") != expected_index:
                raise AuditError(
                    f"{raw.get('structure_id')}: nonsequential subnet indices"
                )
            single = parse_topology(source.get("single_node") or {})
            all_node = parse_topology(source.get("all_node") or {})
            subnets.append(
                {
                    "subnet_index": expected_index,
                    "single_node": single,
                    "all_node": all_node,
                    "single_all_agree": (
                        single["status"] == "SUCCESS"
                        and all_node["status"] == "SUCCESS"
                        and single["dimension"] == all_node["dimension"]
                        and single["topology_key"] == all_node["topology_key"]
                    ),
                }
            )
        count = len(subnets)
        if (
            raw.get("interpenetrated_subnet_count") != count
            or raw.get("catenation_degree") != count
        ):
            raise AuditError(
                f"{raw.get('structure_id')}: inconsistent subnet count"
            )
        successful = bool(subnets) and all(
            subnet["single_node"]["status"] == "SUCCESS"
            and subnet["all_node"]["status"] == "SUCCESS"
            for subnet in subnets
        )
        public_status = "SUCCESS" if successful else "PARTIAL"
        dimensions = [
            subnet["single_node"]["dimension"] for subnet in subnets
        ] + [subnet["all_node"]["dimension"] for subnet in subnets]
        keys_single = [subnet["single_node"]["topology_key"] for subnet in subnets]
        keys_all = [subnet["all_node"]["topology_key"] for subnet in subnets]
        dimension = consensus(dimensions) if successful else None
        single_key = consensus(keys_single) if successful else None
        all_key = consensus(keys_all) if successful else None
        agree = (
            all(subnet["single_all_agree"] for subnet in subnets)
            if successful
            else None
        )
    else:
        public_status = status or "ERROR"
        successful = False
        count = None
        dimension = None
        single_key = None
        all_key = None
        agree = None
        if raw.get("subnets"):
            raise AuditError(
                f"{raw.get('structure_id')}: error record has subnets"
            )
    public = {
        "structure_id": raw.get("structure_id"),
        "topology_available": successful,
        "execution_status": public_status,
        "network_dimension": dimension,
        "interpenetrated_subnet_count": count,
        "catenation_degree": count,
        "single_node_net": single_key,
        "all_node_net": all_key,
        "single_all_agree": agree,
    }
    rich = {
        **public,
        "cif_sha256": raw.get("cif_sha256"),
        "runtime_seconds": raw.get("runtime_seconds"),
        "software": raw.get("software"),
        "method": raw.get("method"),
        "subnets": subnets,
        "error": (
            {
                "type": raw.get("error_type"),
                "message": raw.get("error_message"),
            }
            if status != "SUCCESS"
            else None
        ),
    }
    return public, rich


def load_manifest(path: Path) -> Tuple[list[str], Dict[str, str]]:
    order = []
    hashes = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            structure_id = row["structure_id"]
            if structure_id in hashes:
                raise AuditError(f"duplicate manifest ID: {structure_id}")
            order.append(structure_id)
            hashes[structure_id] = row["cif_sha256"]
    return order, hashes


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--manifest", type=Path, required=True)
    result.add_argument("--results-root", type=Path, required=True)
    result.add_argument("--legacy-csv", type=Path)
    result.add_argument("--expected-shards", type=int, default=96)
    result.add_argument("--output-root", type=Path, required=True)
    return result


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parser().parse_args(argv)
    order, expected_hashes = load_manifest(args.manifest)
    temporary = sorted(args.results_root.glob("*.tmp.*"))
    if temporary:
        raise AuditError(f"{len(temporary)} incomplete temporary shard files")
    shards = sorted(args.results_root.glob("shard_*.jsonl"))
    if len(shards) != args.expected_shards:
        raise AuditError(
            f"completed shards {len(shards)} != {args.expected_shards}"
        )
    raw_by_id = {}
    for shard in shards:
        with shard.open("r", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                raw = json.loads(line)
                structure_id = str(raw.get("structure_id") or "")
                if structure_id in raw_by_id:
                    raise AuditError(f"duplicate result: {structure_id}")
                if expected_hashes.get(structure_id) != raw.get("cif_sha256"):
                    raise AuditError(
                        f"{shard}:{line_number}: CIF binding mismatch"
                    )
                serialized = json.dumps(raw, ensure_ascii=False)
                if any(token in serialized for token in ("/home/", "/scratch/", "/old-home/")):
                    raise AuditError(
                        f"{structure_id}: private path leaked into result"
                    )
                raw_by_id[structure_id] = raw
    missing = sorted(set(order).difference(raw_by_id))
    unexpected = sorted(set(raw_by_id).difference(order))
    if missing or unexpected:
        raise AuditError(
            f"membership mismatch: missing={missing[:10]}, unexpected={unexpected[:10]}"
        )

    legacy = None
    if args.legacy_csv is not None:
        legacy = {
            row["structure_id"]: row
            for row in csv.DictReader(
                args.legacy_csv.open("r", encoding="utf-8", newline="")
            )
        }
        if set(legacy) != set(order):
            raise AuditError(
                "legacy comparison membership differs from manifest"
            )

    args.output_root.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_root / "topology_features.csv"
    json_path = args.output_root / "topology_features.jsonl"
    summary_path = args.output_root / "summary.json"
    status_counts = collections.Counter()
    dimension_counts = collections.Counter()
    error_types = collections.Counter()
    legacy_differences = collections.Counter()
    runtimes = []
    with csv_path.open("w", encoding="utf-8", newline="") as csv_handle, (
        json_path.open("w", encoding="utf-8")
    ) as json_handle:
        writer = csv.DictWriter(csv_handle, fieldnames=PUBLIC_COLUMNS)
        writer.writeheader()
        for structure_id in order:
            raw = raw_by_id[structure_id]
            public, rich = normalized_record(raw)
            writer.writerow(public)
            json_handle.write(
                json.dumps(
                    rich,
                    sort_keys=True,
                    separators=(",", ":"),
                    ensure_ascii=False,
                    allow_nan=False,
                )
                + "\n"
            )
            status_counts[public["execution_status"]] += 1
            dimension_counts[str(public["network_dimension"])] += 1
            if raw.get("error_type"):
                error_types[str(raw["error_type"])] += 1
            if isinstance(raw.get("runtime_seconds"), (int, float)):
                runtimes.append(float(raw["runtime_seconds"]))
            if legacy is not None:
                old = legacy[structure_id]
                for key, old_key in (
                    ("network_dimension", "network_dimension"),
                    ("catenation_degree", "catenation_degree"),
                    ("single_node_net", "single_node_net"),
                    ("all_node_net", "all_node_net"),
                ):
                    current = (
                        "" if public[key] is None else str(public[key])
                    )
                    if current != str(old.get(old_key) or ""):
                        legacy_differences[key] += 1

    summary = {
        "schema_version": "crystalnets-topology-audit/1.0",
        "audit_status": "PASS",
        "audited_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "manifest": {
            "file_name": args.manifest.name,
            "row_count": len(order),
            "sha256": sha256_file(args.manifest),
        },
        "shard_count": len(shards),
        "result_count": len(raw_by_id),
        "status_counts": dict(sorted(status_counts.items())),
        "network_dimension_counts": dict(sorted(dimension_counts.items())),
        "error_type_counts": dict(sorted(error_types.items())),
        "runtime_seconds": {
            "count": len(runtimes),
            "sum": sum(runtimes),
            "maximum": max(runtimes) if runtimes else None,
        },
        "software": {
            "julia_version": "1.12.6",
            "crystalnets_version": "1.2.0",
        },
        "method": {
            "structure_type": "MOF",
            "clusterings": ["SingleNodes", "AllNodes"],
            "catenation_degree_definition": (
                "number of subnets in InterpenetratedTopologyResult"
            ),
            "single_all_disagreement_is_preserved": True,
            "unknown_nets_keep_compact_key_in_CSV_and_genome_in_JSONL": True,
        },
        "legacy_comparison": {
            "performed": legacy is not None,
            "interpretation": (
                "differences are reported, not treated as audit failures, "
                "because legacy software version/options are unavailable"
            ),
            "field_difference_counts": dict(
                sorted(legacy_differences.items())
            ),
        },
        "outputs": {
            "topology_features.csv": {
                "row_count": len(order),
                "sha256": sha256_file(csv_path),
            },
            "topology_features.jsonl": {
                "row_count": len(order),
                "sha256": sha256_file(json_path),
            },
        },
    }
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
