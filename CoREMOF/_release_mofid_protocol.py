#!/usr/bin/env python3
"""Deterministic, non-mutating MOFid-v1/v2 calculation protocol.

The scientific implementation is intentionally separated from scheduler and
release code.  Heavy scientific imports are delayed until a record is
calculated so that manifests and audits can be inspected without activating
the MOFid runtime.
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import json
import os
import re
import shutil
import time
import traceback
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


RESULT_SCHEMA = "coremof-pinned-mofid-result/1.0"
METHOD_SCHEMA = "coremof-pinned-mofid-method/1.0"
MOFID_SOURCE_COMMIT = "5c1b7d3345fc7f3aca1bb346244962ed70634fc7"
FINAL_DATA_ARCHIVE_SHA256 = (
    "c3d74379dbf26026badf5be1f1452eea3719c0d2c5354842fc29678f968679d1"
)
NODE_ARCHIVE_SHA256 = (
    "eb3ddf256ab3fef25fa4b538453632d56f2135ae4443ecd86072fbc2cc6ce049"
)
EXPECTED_NODE_COUNT = 1182
NEIGHBOR_SKIN_ANGSTROM = 0.3
JAVA_TOOL_OPTIONS = (
    "-Xms32m -Xmx768m -XX:MaxMetaspaceSize=256m "
    "-Xss512k -XX:+UseSerialGC"
)
STRUCTURE_MATCHER = {
    "ltol": 0.25,
    "stol": 1.5,
    "angle_tol": 5,
    "primitive_cell": False,
    "scale": False,
    "comparator": "ElementComparator",
}
PINNED_PYTHON_PACKAGES = {
    "python": "3.9",
    "ase": "3.22.1",
    "pymatgen": "2024.2.8",
    "networkx": "3.2.1",
    "selfies": "2.1.1",
    "numpy": "1.26.4",
}

V1_STATUSES = {
    "SUCCESS",
    "SUCCESS_TOPOLOGY_UNKNOWN",
    "SUCCESS_TOPOLOGY_ERROR",
    "SUCCESS_TOPOLOGY_TIMEOUT",
    "NOT_AVAILABLE_NO_MOF",
    "ERROR",
    "TIMEOUT",
}
V2_STATUSES = V1_STATUSES | {
    "NOT_AVAILABLE_UNMATCHED_NODE",
    "NOT_AVAILABLE_AMBIGUOUS_NODE",
    "ERROR_DECOMPOSITION",
}
V2_SCOPES = {"PUBLISHED_METHOD", "COREMOF_FSR_EXTENSION"}
NODE_FILE_RE = re.compile(r"^(?P<formula>.+)_Type-(?P<type_number>[0-9]+)\.xyz$")
COMMIT_SUFFIX_RE = re.compile(r"\.(?:NO_REF|[0-9a-fA-F]{8})$")


class MofidProtocolError(RuntimeError):
    """The pinned scientific contract could not be satisfied."""


class NodeMatchUnavailable(MofidProtocolError):
    """A v2 node cannot be assigned from the frozen official library."""

    def __init__(
        self,
        status: str,
        formula: str,
        candidates: Sequence[str],
        matches: Sequence[str],
        message: str,
    ) -> None:
        super().__init__(message)
        self.status = status
        self.formula = formula
        self.candidates = list(candidates)
        self.matches = list(matches)


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def tree_receipt(root: Path, files: Iterable[Path]) -> dict[str, Any]:
    root = root.resolve()
    rows: list[dict[str, Any]] = []
    tree_digest = hashlib.sha256()
    for path in sorted({item.resolve() for item in files}, key=str):
        try:
            relative = path.relative_to(root)
        except ValueError as error:
            raise MofidProtocolError(
                f"tree member escapes {root}: {path}"
            ) from error
        if not path.is_file():
            raise MofidProtocolError(f"tree member is not a file: {path}")
        size = path.stat().st_size
        digest = sha256_file(path)
        relative_text = relative.as_posix()
        tree_digest.update(
            f"{relative_text}\0{size}\0{digest}\n".encode("utf-8")
        )
        rows.append(
            {
                "path": relative_text,
                "size_bytes": size,
                "sha256": digest,
            }
        )
    if not rows:
        raise MofidProtocolError(f"tree receipt would be empty: {root}")
    return {
        "root": str(root),
        "file_count": len(rows),
        "size_bytes": sum(row["size_bytes"] for row in rows),
        "tree_sha256": tree_digest.hexdigest(),
        "files": rows,
    }


def validate_tree_receipt(receipt: Mapping[str, Any]) -> None:
    root = Path(str(receipt.get("root") or "")).resolve()
    rows = receipt.get("files")
    if not isinstance(rows, list) or not rows:
        raise MofidProtocolError(f"tree receipt has no files: {root}")
    rebuilt = tree_receipt(
        root,
        [root / str(row.get("path") or "") for row in rows],
    )
    for key in ("file_count", "size_bytes", "tree_sha256", "files"):
        if rebuilt[key] != receipt.get(key):
            raise MofidProtocolError(
                f"tree receipt changed for {root}: field={key}"
            )


def atomic_write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, path)


def clean_text(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    if text.casefold() in {"", "none", "null", "nan", "na", "n/a"}:
        return None
    return text


def strip_record_suffix(value: Any) -> str | None:
    text = clean_text(value)
    if text is None:
        return None
    return text.split(";", 1)[0].strip()


def normalize_mofid_v1(value: Any) -> tuple[str | None, str]:
    """Return public MOFid-v1 plus its terminal scientific status."""
    primary = strip_record_suffix(value)
    if primary is None:
        return None, "ERROR"
    primary = COMMIT_SUFFIX_RE.sub("", primary).strip()
    if re.match(r"^\*?\s*MOFid-v1\.NA(?:\.|$)", primary, re.IGNORECASE):
        return None, "NOT_AVAILABLE_NO_MOF"
    if "MOFid-v1." not in primary:
        return None, "ERROR"
    topology_match = re.search(r"\bMOFid-v1\.([^. ;]+)", primary)
    topology = topology_match.group(1).upper() if topology_match else ""
    topology_tokens = {token.strip() for token in topology.split(",")}
    if "ERROR" in topology_tokens:
        return primary, "SUCCESS_TOPOLOGY_ERROR"
    if "TIMEOUT" in topology_tokens:
        return primary, "SUCCESS_TOPOLOGY_TIMEOUT"
    if "UNKNOWN" in topology_tokens:
        return primary, "SUCCESS_TOPOLOGY_UNKNOWN"
    if topology == "NA":
        return None, "NOT_AVAILABLE_NO_MOF"
    return primary, "SUCCESS"


def canonical_v1_comparison(value: Any) -> str | None:
    primary, _status = normalize_mofid_v1(value)
    if primary is None:
        return None
    # The archived public data contain one lowercase topology spelling.
    primary = re.sub(
        r"(\bMOFid-v1\.)(unknown|error|timeout)(?=\.)",
        lambda match: match.group(1) + match.group(2).upper(),
        primary,
        flags=re.IGNORECASE,
    )
    return re.sub(r"\s+", " ", primary).strip()


def status_for_v2_topology(topology: Any) -> str:
    tokens = {
        token.strip()
        for token in str(topology or "").upper().split(",")
    }
    if "ERROR" in tokens:
        return "SUCCESS_TOPOLOGY_ERROR"
    if "TIMEOUT" in tokens:
        return "SUCCESS_TOPOLOGY_TIMEOUT"
    if "UNKNOWN" in tokens:
        return "SUCCESS_TOPOLOGY_UNKNOWN"
    if tokens.issubset({"", "NA", "NONE"}):
        return "NOT_AVAILABLE_NO_MOF"
    return "SUCCESS"


def formula_from_symbols(symbols: Iterable[str]) -> str:
    counts = Counter(str(symbol) for symbol in symbols)
    return "".join(f"{element}{counts[element]}" for element in sorted(counts))


def parse_node_file_name(file_name: str) -> tuple[str, int, str]:
    match = NODE_FILE_RE.fullmatch(file_name)
    if match is None:
        raise MofidProtocolError(f"invalid official node filename: {file_name}")
    formula = match.group("formula")
    type_number = int(match.group("type_number"))
    return formula, type_number, file_name[:-4]


def load_node_manifest(
    path: Path,
    node_root: Path,
    verify_files: bool = True,
) -> dict[str, list[dict[str, Any]]]:
    import csv

    required = {
        "archive_index",
        "node_label",
        "formula",
        "type_number",
        "file",
        "size_bytes",
        "sha256",
    }
    grouped: dict[str, list[dict[str, Any]]] = {}
    seen_labels: set[str] = set()
    seen_archive_indices: set[int] = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or set(reader.fieldnames) != required:
            raise MofidProtocolError(
                f"node manifest columns differ from contract: {reader.fieldnames}"
            )
        for row in reader:
            label = row["node_label"]
            if label in seen_labels:
                raise MofidProtocolError(f"duplicate node label: {label}")
            seen_labels.add(label)
            archive_index = int(row["archive_index"])
            if archive_index in seen_archive_indices:
                raise MofidProtocolError(
                    f"duplicate official archive index: {archive_index}"
                )
            seen_archive_indices.add(archive_index)
            file_path = node_root / row["file"]
            if verify_files:
                if not file_path.is_file():
                    raise MofidProtocolError(f"node file is missing: {file_path}")
                if file_path.stat().st_size != int(row["size_bytes"]):
                    raise MofidProtocolError(f"node size mismatch: {file_path}")
                if sha256_file(file_path) != row["sha256"]:
                    raise MofidProtocolError(f"node hash mismatch: {file_path}")
            normalized = {
                **row,
                "archive_index": archive_index,
                "type_number": int(row["type_number"]),
                "path": file_path,
            }
            grouped.setdefault(row["formula"], []).append(normalized)
    if len(seen_labels) != EXPECTED_NODE_COUNT:
        raise MofidProtocolError(
            f"official node count is {len(seen_labels)}, expected {EXPECTED_NODE_COUNT}"
        )
    for candidates in grouped.values():
        candidates.sort(key=lambda row: row["archive_index"])
    return grouped


class PinnedV2Builder:
    """Convert a successful v1 decomposition into a frozen-library v2 ID."""

    def __init__(
        self,
        node_root: Path,
        node_manifest: Path,
        verify_library_files: bool = False,
    ) -> None:
        from ase import neighborlist
        from ase.build import sort as ase_sort
        from ase.io import read as ase_read
        import ase
        import networkx
        import selfies
        from pymatgen.analysis.structure_matcher import (
            ElementComparator,
            StructureMatcher,
        )
        from pymatgen.core.lattice import Lattice
        from pymatgen.core.structure import Structure

        self.ase = ase
        self.ase_read = ase_read
        self.ase_sort = ase_sort
        self.neighborlist = neighborlist
        self.networkx = networkx
        self.selfies = selfies
        self.ElementComparator = ElementComparator
        self.StructureMatcher = StructureMatcher
        self.Lattice = Lattice
        self.Structure = Structure
        self.node_root = node_root
        self.node_manifest = node_manifest
        self.library = load_node_manifest(
            node_manifest,
            node_root,
            verify_files=verify_library_files,
        )

    def _remove_pbc_cuts(self, atoms: Any) -> Any:
        import collections
        import numpy as np

        try:
            cutoffs = [1.4 * value for value in self.neighborlist.natural_cutoffs(atoms)]
            left, right, displacement = self.neighborlist.neighbor_list(
                "ijD", atoms, cutoff=cutoffs
            )
            neighbors: list[list[tuple[int, Any]]] = [[] for _ in atoms]
            for i, j, delta in zip(left, right, displacement):
                neighbors[int(i)].append((int(j), delta))
            visited = [False for _ in atoms]
            queue: collections.deque[tuple[int, Any]] = collections.deque()
            queue.append((0, np.array([0.0, 0.0, 0.0])))
            positions: dict[int, Any] = {}
            while queue:
                index, position = queue.pop()
                if visited[index]:
                    continue
                visited[index] = True
                positions[index] = position
                for neighbor, delta in neighbors[index]:
                    if not visited[neighbor]:
                        queue.append((neighbor, position + delta))
            if len(positions) != len(atoms):
                return atoms
            center = np.sum(atoms.get_cell(), axis=0) * 0.5
            centroid = sum(positions.values()) / len(positions)
            shifted = [
                positions[index] - centroid + center for index in range(len(atoms))
            ]
            rebuilt = self.ase.Atoms(
                symbols=atoms.get_chemical_symbols(),
                positions=shifted,
                pbc=True,
                cell=atoms.get_cell(),
            )
            spans = [
                float(rebuilt.positions[:, axis].max() - rebuilt.positions[:, axis].min())
                for axis in range(3)
            ]
            side = max(spans) + 2.0
            rebuilt.set_cell([side, side, side, 90, 90, 90])
            rebuilt.positions += (
                rebuilt.cell.cellpar()[0:3] / 2 - rebuilt.get_center_of_mass()
            )
            return rebuilt
        except Exception:
            return atoms

    def _to_pymatgen(self, atoms: Any) -> Any:
        # Preserve the published/reference implementation exactly: the ASE
        # position array is passed to Structure without the cartesian flag.
        return self.Structure(
            self.Lattice(atoms.cell),
            atoms.get_chemical_symbols(),
            atoms.get_positions(),
        )

    def split_nodes(self, nodes_cif: Path) -> list[tuple[str, Any]]:
        try:
            atoms = self.ase_read(nodes_cif)
        except Exception as error:
            raise MofidProtocolError(
                f"cannot read AllNode nodes CIF: {nodes_cif}: {error}"
            ) from error
        cutoffs = self.neighborlist.natural_cutoffs(atoms)
        neighbors = self.neighborlist.NeighborList(
            cutoffs,
            self_interaction=False,
            bothways=True,
            skin=NEIGHBOR_SKIN_ANGSTROM,
        )
        neighbors.update(atoms)
        graph = self.networkx.Graph()
        graph.add_nodes_from(range(len(atoms)))
        for index in range(len(atoms)):
            for neighbor in neighbors.get_neighbors(index)[0]:
                graph.add_edge(index, int(neighbor))
        components = sorted(
            self.networkx.connected_components(graph),
            key=lambda component: (-len(component), min(component)),
        )
        unique: list[tuple[str, Any]] = []
        seen_formulas: set[str] = set()
        for component in components:
            indices = sorted(int(index) for index in component)
            fragment = self.ase_sort(atoms[indices])
            formula = formula_from_symbols(fragment.get_chemical_symbols())
            if formula not in seen_formulas:
                seen_formulas.add(formula)
                unique.append((formula, fragment))
        if not unique:
            raise MofidProtocolError("AllNode decomposition produced no node components")
        return unique

    def match_node(self, formula: str, atoms: Any) -> dict[str, Any]:
        candidates = self.library.get(formula, [])
        labels = [str(candidate["node_label"]) for candidate in candidates]
        if not candidates:
            raise NodeMatchUnavailable(
                "NOT_AVAILABLE_UNMATCHED_NODE",
                formula,
                [],
                [],
                f"formula {formula} is absent from the frozen official library",
            )
        matcher = self.StructureMatcher(
            ltol=STRUCTURE_MATCHER["ltol"],
            stol=STRUCTURE_MATCHER["stol"],
            angle_tol=STRUCTURE_MATCHER["angle_tol"],
            primitive_cell=STRUCTURE_MATCHER["primitive_cell"],
            scale=STRUCTURE_MATCHER["scale"],
            comparator=self.ElementComparator(),
        )
        query = self._to_pymatgen(self._remove_pbc_cuts(atoms))
        matches: list[str] = []
        selected: dict[str, Any] | None = None
        for candidate in candidates:
            reference_atoms = self.ase_read(candidate["path"])
            reference = self._to_pymatgen(self._remove_pbc_cuts(reference_atoms))
            if matcher.fit(query, reference):
                matches.append(str(candidate["node_label"]))
                if selected is None:
                    selected = candidate
        if selected is None:
            raise NodeMatchUnavailable(
                "NOT_AVAILABLE_UNMATCHED_NODE",
                formula,
                labels,
                [],
                f"formula {formula} has no structural match in the official library",
            )
        # The final archive has a scientifically meaningful, frozen member
        # order.  Its reference implementation selects the first match.  We
        # bind and expose that order explicitly rather than relying on readdir
        # or glob order.  A non-unique archive rank would be genuinely
        # ambiguous and is rejected by the manifest loader.
        return {
            "formula": formula,
            "selected_node_label": selected["node_label"],
            "selected_archive_index": selected["archive_index"],
            "candidate_node_labels": labels,
            "matching_node_labels": matches,
            "multiple_structural_matches": len(matches) > 1,
            "selection_rule": "FIRST_MATCH_IN_FROZEN_OFFICIAL_ARCHIVE_ORDER",
        }

    def build(
        self,
        identifiers: Mapping[str, Any],
        output_path: Path,
        refname: str,
    ) -> tuple[str, str, list[dict[str, Any]], list[str]]:
        topology = str(identifiers.get("topology") or "")
        cat = identifiers.get("cat")
        if not topology or topology.upper() == "NA" or cat is None:
            raise NodeMatchUnavailable(
                "NOT_AVAILABLE_NO_MOF",
                "",
                [],
                [],
                "v1 decomposition did not produce topology and catenation",
            )
        node_components = self.split_nodes(output_path / "AllNode/nodes.cif")
        node_matches = [
            self.match_node(formula, atoms) for formula, atoms in node_components
        ]
        linker_selfies = [
            self.selfies.encoder(str(linker))
            for linker in list(identifiers.get("smiles_linkers") or [])
        ]
        nodes_part = ".".join(
            f"[{match['selected_node_label']}]" for match in node_matches
        )
        linkers_part = ".".join(linker_selfies)
        raw = (
            f"{nodes_part}.{linkers_part} "
            f"MOFid-v2.{topology}.cat{cat};{refname}"
        )
        public = raw.rsplit(";", 1)[0]
        return raw, public, node_matches, linker_selfies


def scope_for_variant(variant: Any) -> str:
    return (
        "COREMOF_FSR_EXTENSION"
        if str(variant or "").upper() == "FSR"
        else "PUBLISHED_METHOD"
    )


def _base_record(
    row: Mapping[str, Any],
    cif_sha256: str,
    input_manifest_sha256: str,
    node_library_manifest_sha256: str,
    method_manifest_sha256: str,
) -> dict[str, Any]:
    scope = scope_for_variant(row.get("structure_variant"))
    return {
        "schema_version": RESULT_SCHEMA,
        "structure_id": str(row["structure_id"]),
        "source_database": str(row.get("source_database") or ""),
        "source_id": str(row.get("source_id") or ""),
        "structure_variant": str(row.get("structure_variant") or ""),
        "cif_file": str(row["cif_file"]),
        "cif_sha256": cif_sha256,
        "input_manifest_sha256": input_manifest_sha256,
        "node_library_manifest_sha256": node_library_manifest_sha256,
        "method_manifest_sha256": method_manifest_sha256,
        "existing_mofid_v1": clean_text(row.get("existing_mofid_v1")),
        "existing_mofid_v2": clean_text(row.get("existing_mofid_v2")),
        "final_mofid_v1": None,
        "mofid_v1_status": "ERROR",
        "final_mofid_v2": None,
        "mofid_v2_status": "ERROR",
        "mofid_v2_scope": scope,
        "v1_comparison": "NOT_RUN",
        "calculation": {},
    }


def timeout_record(
    row: Mapping[str, Any],
    cif_sha256: str,
    input_manifest_sha256: str,
    node_library_manifest_sha256: str,
    method_manifest_sha256: str,
    timeout_seconds: int,
) -> dict[str, Any]:
    record = _base_record(
        row,
        cif_sha256,
        input_manifest_sha256,
        node_library_manifest_sha256,
        method_manifest_sha256,
    )
    existing = clean_text(row.get("existing_mofid_v1"))
    existing_public, existing_status = normalize_mofid_v1(existing)
    record["final_mofid_v1"] = existing_public
    record["mofid_v1_status"] = existing_status if existing_public else "TIMEOUT"
    record["mofid_v2_status"] = "TIMEOUT"
    record["v1_comparison"] = (
        "NOT_REVALIDATED_TIMEOUT" if existing_public else "NOT_AVAILABLE_TIMEOUT"
    )
    record["calculation"] = {
        "execution_status": "TIMEOUT",
        "timeout_seconds": timeout_seconds,
    }
    return record


def calculate_record(
    row: Mapping[str, Any],
    cif_path: Path,
    work_root: Path,
    node_root: Path,
    node_manifest: Path,
    input_manifest_sha256: str,
    node_library_manifest_sha256: str,
    method_manifest_sha256: str,
) -> dict[str, Any]:
    """Calculate one record.  Call this inside a killable worker process."""
    started = time.monotonic()
    expected_cif_sha256 = str(row["cif_sha256"])
    actual_cif_sha256 = sha256_file(cif_path)
    if actual_cif_sha256 != expected_cif_sha256:
        raise MofidProtocolError(
            f"CIF hash mismatch for {row['structure_id']}: "
            f"{actual_cif_sha256} != {expected_cif_sha256}"
        )
    record = _base_record(
        row,
        actual_cif_sha256,
        input_manifest_sha256,
        node_library_manifest_sha256,
        method_manifest_sha256,
    )
    structure_id = str(row["structure_id"])
    output_path = work_root / structure_id / "Output"
    shutil.rmtree(output_path.parent, ignore_errors=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    stdout = io.StringIO()
    stderr = io.StringIO()
    try:
        from mofid.run_mofid import cif2mofid

        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            identifiers = cif2mofid(str(cif_path), str(output_path))
        raw_v1 = str(identifiers.get("mofid") or "")
        calculated_v1, calculated_v1_status = normalize_mofid_v1(raw_v1)
        existing_v1 = clean_text(row.get("existing_mofid_v1"))
        if existing_v1 is not None:
            existing_public, existing_status = normalize_mofid_v1(existing_v1)
            record["final_mofid_v1"] = existing_public
            record["mofid_v1_status"] = existing_status
            if canonical_v1_comparison(existing_v1) == canonical_v1_comparison(raw_v1):
                record["v1_comparison"] = "MATCH"
            else:
                record["v1_comparison"] = "SUBSTANTIVE_MISMATCH"
        else:
            record["final_mofid_v1"] = calculated_v1
            record["mofid_v1_status"] = calculated_v1_status
            record["v1_comparison"] = (
                "FILLED_MISSING"
                if calculated_v1 is not None
                else "REMAINS_NOT_AVAILABLE"
            )
        calculation: dict[str, Any] = {
            "execution_status": "SUCCESS",
            "raw_mofid_v1": raw_v1,
            "calculated_mofid_v1": calculated_v1,
            "calculated_mofid_v1_status": calculated_v1_status,
            "mofkey": clean_text(identifiers.get("mofkey")),
            "smiles_nodes": list(identifiers.get("smiles_nodes") or []),
            "smiles_linkers": list(identifiers.get("smiles_linkers") or []),
            "topology": identifiers.get("topology"),
            "catenation_degree": identifiers.get("cat"),
        }
        if calculated_v1 is None:
            if calculated_v1_status == "NOT_AVAILABLE_NO_MOF":
                record["mofid_v2_status"] = "NOT_AVAILABLE_NO_MOF"
                calculation["v2_execution_status"] = "NOT_AVAILABLE_NO_MOF"
            else:
                record["mofid_v2_status"] = "ERROR_DECOMPOSITION"
                calculation["v2_execution_status"] = "ERROR_DECOMPOSITION"
                calculation["v2_error"] = {
                    "type": "InvalidV1Decomposition",
                    "message": (
                        "v1 execution returned without a valid MOFid or a "
                        "genuine no-MOF marker"
                    ),
                }
        else:
            builder = PinnedV2Builder(node_root, node_manifest)
            try:
                raw_v2, public_v2, node_matches, linker_selfies = builder.build(
                    identifiers, output_path, structure_id
                )
                record["final_mofid_v2"] = public_v2
                record["mofid_v2_status"] = status_for_v2_topology(
                    identifiers.get("topology")
                )
                calculation.update(
                    {
                        "v2_execution_status": "SUCCESS",
                        "raw_mofid_v2": raw_v2,
                        "calculated_mofid_v2": public_v2,
                        "node_matches": node_matches,
                        "linker_selfies": linker_selfies,
                    }
                )
            except NodeMatchUnavailable as error:
                record["mofid_v2_status"] = error.status
                calculation.update(
                    {
                        "v2_execution_status": error.status,
                        "v2_error": {
                            "type": type(error).__name__,
                            "message": str(error),
                            "formula": error.formula,
                            "candidate_node_labels": error.candidates,
                            "matching_node_labels": error.matches,
                        },
                    }
                )
            except Exception as error:
                record["mofid_v2_status"] = "ERROR_DECOMPOSITION"
                calculation.update(
                    {
                        "v2_execution_status": "ERROR_DECOMPOSITION",
                        "v2_error": {
                            "type": type(error).__name__,
                            "message": str(error),
                            "traceback": traceback.format_exc(),
                        },
                    }
                )
        calculation["captured_stdout"] = stdout.getvalue()
        calculation["captured_stderr"] = stderr.getvalue()
        calculation["runtime_seconds"] = round(time.monotonic() - started, 6)
        record["calculation"] = calculation
        return record
    except Exception as error:
        existing = clean_text(row.get("existing_mofid_v1"))
        existing_public, existing_status = normalize_mofid_v1(existing)
        record["final_mofid_v1"] = existing_public
        record["mofid_v1_status"] = existing_status if existing_public else "ERROR"
        record["v1_comparison"] = (
            "NOT_REVALIDATED_ERROR" if existing_public else "NOT_AVAILABLE_ERROR"
        )
        record["mofid_v2_status"] = "ERROR"
        record["calculation"] = {
            "execution_status": "ERROR",
            "error": {
                "type": type(error).__name__,
                "message": str(error),
                "traceback": traceback.format_exc(),
            },
            "captured_stdout": stdout.getvalue(),
            "captured_stderr": stderr.getvalue(),
            "runtime_seconds": round(time.monotonic() - started, 6),
        }
        return record
    finally:
        shutil.rmtree(output_path.parent, ignore_errors=True)


def validate_result_record(
    record: Mapping[str, Any],
    expected_structure_id: str | None = None,
) -> None:
    if record.get("schema_version") != RESULT_SCHEMA:
        raise MofidProtocolError("record schema version is invalid")
    structure_id = str(record.get("structure_id") or "")
    if not structure_id:
        raise MofidProtocolError("record structure_id is empty")
    if expected_structure_id is not None and structure_id != expected_structure_id:
        raise MofidProtocolError(
            f"record ID mismatch: {structure_id} != {expected_structure_id}"
        )
    for key in (
        "cif_sha256",
        "input_manifest_sha256",
        "node_library_manifest_sha256",
        "method_manifest_sha256",
    ):
        if re.fullmatch(r"[0-9a-f]{64}", str(record.get(key) or "")) is None:
            raise MofidProtocolError(f"invalid {key} in {structure_id}")
    if record.get("mofid_v1_status") not in V1_STATUSES:
        raise MofidProtocolError(
            f"invalid v1 status for {structure_id}: {record.get('mofid_v1_status')}"
        )
    v1_available = clean_text(record.get("final_mofid_v1")) is not None
    if str(record.get("mofid_v1_status")).startswith("SUCCESS"):
        if not v1_available:
            raise MofidProtocolError(f"successful v1 is null for {structure_id}")
    elif v1_available:
        raise MofidProtocolError(
            f"unavailable/error v1 is non-null for {structure_id}"
        )
    if record.get("mofid_v2_status") not in V2_STATUSES:
        raise MofidProtocolError(
            f"invalid v2 status for {structure_id}: {record.get('mofid_v2_status')}"
        )
    if record.get("mofid_v2_scope") not in V2_SCOPES:
        raise MofidProtocolError(
            f"invalid v2 scope for {structure_id}: {record.get('mofid_v2_scope')}"
        )
    if record.get("mofid_v2_status", "").startswith("SUCCESS"):
        if clean_text(record.get("final_mofid_v2")) is None:
            raise MofidProtocolError(f"successful v2 is null for {structure_id}")
    elif record.get("final_mofid_v2") is not None:
        raise MofidProtocolError(f"unavailable v2 is non-null for {structure_id}")
