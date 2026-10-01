#!/usr/bin/env python3
"""Run pinned, fresh MOFClassifier core inference without touching source CIFs.

The runner accepts a small, hash-bound manifest and writes one atomic JSON
record per selected row.  CIF parsing is performed only on a private copy.
Existing records are resumable only after their complete input, method, asset,
runtime, and result contracts have been revalidated.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import io
import json
import math
import os
import platform
import random
import re
import shutil
import stat
import tempfile
import types
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple


TOOL_VERSION = "1.0.0"
CONFIG_SCHEMA_VERSION = "1.0"
MANIFEST_SCHEMA_VERSION = "1.0"
RECORD_SCHEMA_VERSION = "1.0"
EXPECTED_BAG_COUNT = 100
EXPECTED_THRESHOLD = 0.6
EXPECTED_MEAN_AGGREGATION = "float64_math_fsum_divide_by_100"
EXPECTED_MAX_ROWS_PER_INVOCATION = 64
EXPECTED_CPU_INTRAOP_THREADS = 1
EXPECTED_CPU_INTEROP_THREADS = 1
EXPECTED_RANDOM_SEED = 20260719
EXPECTED_DETERMINISTIC_ALGORITHMS = True
EXPECTED_RUNTIME_APPROVAL_STATE = "UNAPPROVED_REQUIRES_RELEASE_REVIEW"
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
UID_RE = re.compile(r"^coremof:[a-z0-9][a-z0-9._-]*:[^:]+$")
try:
    from .identifiers import parse_core_id
except ImportError:  # Isolated worker loads this module directly by file path.
    from identifiers import parse_core_id

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = REPO_ROOT / "dataset_split/config/mofclassifier_fresh_core_v1.json"

MANIFEST_FIELDS = (
    "manifest_schema_version",
    "row_index",
    "persistent_uid",
    "structure_id",
    "canonical_cif_path",
    "canonical_cif_size_bytes",
    "canonical_cif_sha256",
)

RECORD_FIELDS = (
    "record_schema_version",
    "runner_tool_version",
    "checker_name",
    "method_id",
    "fresh_model_inference",
    "mechanistic_hard_fail_eligible",
    "method_config_path",
    "method_config_sha256",
    "asset_bundle_sha256",
    "package_distribution",
    "package_version",
    "package_root",
    "source_module_path",
    "source_module_sha256",
    "atom_initializer_path",
    "atom_initializer_sha256",
    "model_family",
    "checkpoint_bag_count",
    "execution_device",
    "data_loader_num_workers",
    "data_loader_shuffle",
    "data_loader_pin_memory",
    "positive_class_index",
    "pass_threshold",
    "mean_aggregation",
    "maximum_rows_per_invocation",
    "cpu_intraop_threads",
    "cpu_interop_threads",
    "random_seed",
    "deterministic_algorithms",
    "runtime_approval_state",
    "runtime",
    "runtime_fingerprint_sha256",
    "source_manifest_path",
    "source_manifest_sha256",
    "manifest_row_index",
    "persistent_uid",
    "structure_id",
    "canonical_cif_path",
    "canonical_cif_size_bytes",
    "canonical_cif_sha256",
    "canonical_pre_run_sha256",
    "canonical_post_run_sha256",
    "canonical_unchanged",
    "private_copy_used",
    "private_copy_initial_sha256",
    "private_copy_model_input_sha256",
    "private_copy_rewritten",
    "execution_status",
    "completed_bag_count",
    "bag_scores",
    "mean_score",
    "operational_vote",
    "error_stage",
    "error_type",
    "error_message",
    "failed_bag_index",
)


class MofclassifierRunError(ValueError):
    """An input or existing output violated the fail-closed run contract."""


def _reject_duplicate_keys(pairs: Sequence[Tuple[str, Any]]) -> Dict[str, Any]:
    result: Dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise MofclassifierRunError(f"duplicate JSON key: {key!r}")
        result[key] = value
    return result


def load_json(path: Path) -> Any:
    path = path if path.is_absolute() else (Path.cwd() / path).absolute()
    payload, _size, _digest = read_regular_file_bytes(path)
    try:
        return json.loads(
            payload.decode("utf-8"), object_pairs_hook=_reject_duplicate_keys
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise MofclassifierRunError(f"invalid JSON {path}: {exc}") from exc


def read_regular_file_bytes(path: Path) -> Tuple[bytes, int, str]:
    """Read exactly the bytes that are hashed, rejecting an in-read mutation."""
    if not path.is_absolute():
        raise MofclassifierRunError(f"path must be absolute: {path}")
    try:
        before = path.lstat()
        if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode):
            raise MofclassifierRunError(f"file must be a regular non-symlink: {path}")
        payload = path.read_bytes()
        after = path.lstat()
    except OSError as exc:
        raise MofclassifierRunError(f"file is unavailable: {path}: {exc}") from exc
    identity_before = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
    identity_after = (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
    if identity_before != identity_after or len(payload) != after.st_size:
        raise MofclassifierRunError(f"file changed while it was being read: {path}")
    digest = hashlib.sha256(payload).hexdigest()
    return payload, len(payload), digest


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_json_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        allow_nan=False,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def atomic_write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}.{uuid.uuid4().hex}")
    try:
        with temporary.open("x", encoding="utf-8") as handle:
            json.dump(value, handle, allow_nan=False, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _strict_hash(value: Any, context: str) -> str:
    token = str(value)
    if not SHA256_RE.fullmatch(token):
        raise MofclassifierRunError(f"{context} must be a lowercase SHA-256 digest")
    return token


def _canonical_positive_int(value: Any, context: str) -> int:
    token = str(value)
    try:
        number = int(token)
    except ValueError as exc:
        raise MofclassifierRunError(f"{context} must be a positive integer") from exc
    if number <= 0 or str(number) != token:
        raise MofclassifierRunError(f"{context} must be a canonical positive integer")
    return number


def observe_regular_file(path: Path) -> Tuple[int, str]:
    """Hash a non-symlink regular file and reject mutation during hashing."""
    if not path.is_absolute():
        raise MofclassifierRunError(f"path must be absolute: {path}")
    try:
        before = path.lstat()
    except OSError as exc:
        raise MofclassifierRunError(f"file is unavailable: {path}: {exc}") from exc
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode):
        raise MofclassifierRunError(f"file must be a regular non-symlink: {path}")
    digest = sha256_file(path)
    after = path.lstat()
    identity_before = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
    identity_after = (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
    if identity_before != identity_after:
        raise MofclassifierRunError(f"file changed while it was being hashed: {path}")
    return int(after.st_size), digest


@dataclass(frozen=True)
class Asset:
    role: str
    relative_path: str
    path: Path
    size_bytes: int
    sha256: str
    bag_index: Optional[int] = None

    def fingerprint_value(self) -> Dict[str, Any]:
        return {
            "role": self.role,
            "bag_index": self.bag_index,
            "relative_path": self.relative_path,
            "size_bytes": self.size_bytes,
            "sha256": self.sha256,
        }


@dataclass(frozen=True)
class MethodContract:
    config_path: Path
    config_sha256: str
    method_id: str
    package_distribution: str
    package_version: str
    package_root: Path
    model_family: str
    threshold: float
    positive_class_index: int
    mean_aggregation: str
    maximum_rows_per_invocation: int
    cpu_intraop_threads: int
    cpu_interop_threads: int
    random_seed: int
    deterministic_algorithms: bool
    runtime_approval_state: str
    source_module: Asset
    atom_initializer: Asset
    checkpoints: Tuple[Asset, ...]
    asset_bundle_sha256: str


def _asset_from_config(
    package_root: Path,
    raw: Mapping[str, Any],
    role: str,
    bag_index: Optional[int] = None,
) -> Asset:
    expected = {"relative_path", "size_bytes", "sha256"}
    if set(raw) != expected:
        raise MofclassifierRunError(
            f"{role} asset fields differ from {sorted(expected)}: {sorted(raw)}"
        )
    relative = str(raw["relative_path"])
    relative_path = Path(relative)
    if relative_path.is_absolute() or ".." in relative_path.parts or relative_path.as_posix() != relative:
        raise MofclassifierRunError(f"unsafe/non-canonical {role} relative path: {relative!r}")
    path = (package_root / relative_path).resolve()
    try:
        path.relative_to(package_root)
    except ValueError as exc:
        raise MofclassifierRunError(f"{role} asset escapes package root: {path}") from exc
    size = _canonical_positive_int(raw["size_bytes"], f"{role} size_bytes")
    digest = _strict_hash(raw["sha256"], f"{role} sha256")
    observed_size, observed_digest = observe_regular_file(path)
    if observed_size != size or observed_digest != digest:
        raise MofclassifierRunError(
            f"{role} asset mismatch for {path}: expected size/hash {size}/{digest}, "
            f"observed {observed_size}/{observed_digest}"
        )
    return Asset(role, relative, path, size, digest, bag_index)


def load_method_contract(config_path: Path, package_root_override: Optional[Path] = None) -> MethodContract:
    config_path = config_path.resolve()
    config_payload, _config_size, config_sha256 = read_regular_file_bytes(config_path)
    try:
        raw = json.loads(
            config_payload.decode("utf-8"), object_pairs_hook=_reject_duplicate_keys
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise MofclassifierRunError(f"invalid method config JSON: {config_path}: {exc}") from exc
    if not isinstance(raw, Mapping):
        raise MofclassifierRunError("method config must be a JSON object")
    expected_fields = {
        "config_schema_version",
        "method_id",
        "package_distribution",
        "package_version",
        "package_root_hint",
        "model_family",
        "bag_count",
        "pass_threshold",
        "positive_class_index",
        "mean_aggregation",
        "maximum_rows_per_invocation",
        "cpu_intraop_threads",
        "cpu_interop_threads",
        "random_seed",
        "deterministic_algorithms",
        "runtime_approval_state",
        "execution_device",
        "data_loader_num_workers",
        "data_loader_shuffle",
        "data_loader_pin_memory",
        "source_module",
        "atom_initializer",
        "checkpoints",
    }
    if set(raw) != expected_fields:
        raise MofclassifierRunError(
            f"method config fields differ from contract: missing={sorted(expected_fields-set(raw))}, "
            f"extra={sorted(set(raw)-expected_fields)}"
        )
    if raw["config_schema_version"] != CONFIG_SCHEMA_VERSION:
        raise MofclassifierRunError("unsupported method config schema")
    if raw["model_family"] != "core" or raw["execution_device"] != "cpu":
        raise MofclassifierRunError("only the pinned core CPU method is supported")
    if raw["bag_count"] != EXPECTED_BAG_COUNT:
        raise MofclassifierRunError("method must contain exactly 100 bags")
    if raw["pass_threshold"] != EXPECTED_THRESHOLD:
        raise MofclassifierRunError("method threshold must be exactly 0.6")
    if raw["positive_class_index"] != 1:
        raise MofclassifierRunError("positive class index must be 1")
    if raw["mean_aggregation"] != EXPECTED_MEAN_AGGREGATION:
        raise MofclassifierRunError(
            f"mean aggregation must be {EXPECTED_MEAN_AGGREGATION!r}"
        )
    deterministic_settings = {
        "maximum_rows_per_invocation": EXPECTED_MAX_ROWS_PER_INVOCATION,
        "cpu_intraop_threads": EXPECTED_CPU_INTRAOP_THREADS,
        "cpu_interop_threads": EXPECTED_CPU_INTEROP_THREADS,
        "random_seed": EXPECTED_RANDOM_SEED,
        "deterministic_algorithms": EXPECTED_DETERMINISTIC_ALGORITHMS,
        "runtime_approval_state": EXPECTED_RUNTIME_APPROVAL_STATE,
    }
    mismatched_settings = {
        key: (raw[key], expected)
        for key, expected in deterministic_settings.items()
        if raw[key] != expected
    }
    if mismatched_settings:
        raise MofclassifierRunError(
            f"method deterministic/runtime settings mismatch: {mismatched_settings}"
        )
    if raw["data_loader_num_workers"] != 0:
        raise MofclassifierRunError("data loader num_workers must be 0")
    if raw["data_loader_shuffle"] is not False or raw["data_loader_pin_memory"] is not False:
        raise MofclassifierRunError("data loader shuffle and pin_memory must both be false")
    package_root = (
        package_root_override.resolve()
        if package_root_override is not None
        else Path(str(raw["package_root_hint"])).resolve()
    )
    if not package_root.is_dir():
        raise MofclassifierRunError(f"MOFClassifier package root is missing: {package_root}")

    source = _asset_from_config(package_root, raw["source_module"], "source_module")
    atom = _asset_from_config(package_root, raw["atom_initializer"], "atom_initializer")
    checkpoint_raw = raw["checkpoints"]
    if not isinstance(checkpoint_raw, list) or len(checkpoint_raw) != EXPECTED_BAG_COUNT:
        raise MofclassifierRunError("checkpoints must be a 100-entry array")
    checkpoints: List[Asset] = []
    for expected_bag, entry in enumerate(checkpoint_raw, start=1):
        if not isinstance(entry, Mapping) or set(entry) != {"bag_index", "relative_path", "size_bytes", "sha256"}:
            raise MofclassifierRunError(f"invalid checkpoint entry for bag {expected_bag}")
        if entry["bag_index"] != expected_bag:
            raise MofclassifierRunError(
                f"checkpoint bag order mismatch: expected {expected_bag}, found {entry['bag_index']!r}"
            )
        asset_entry = {key: entry[key] for key in ("relative_path", "size_bytes", "sha256")}
        asset = _asset_from_config(package_root, asset_entry, "checkpoint", expected_bag)
        expected_relative = f"models/checkpoint_bag_{expected_bag}.pth.tar"
        if asset.relative_path != expected_relative:
            raise MofclassifierRunError(
                f"bag {expected_bag} checkpoint path must be {expected_relative!r}"
            )
        checkpoints.append(asset)
    all_assets = [source, atom] + checkpoints
    bundle = canonical_json_sha256([asset.fingerprint_value() for asset in all_assets])
    return MethodContract(
        config_path=config_path,
        config_sha256=config_sha256,
        method_id=str(raw["method_id"]),
        package_distribution=str(raw["package_distribution"]),
        package_version=str(raw["package_version"]),
        package_root=package_root,
        model_family="core",
        threshold=EXPECTED_THRESHOLD,
        positive_class_index=1,
        mean_aggregation=EXPECTED_MEAN_AGGREGATION,
        maximum_rows_per_invocation=EXPECTED_MAX_ROWS_PER_INVOCATION,
        cpu_intraop_threads=EXPECTED_CPU_INTRAOP_THREADS,
        cpu_interop_threads=EXPECTED_CPU_INTEROP_THREADS,
        random_seed=EXPECTED_RANDOM_SEED,
        deterministic_algorithms=EXPECTED_DETERMINISTIC_ALGORITHMS,
        runtime_approval_state=EXPECTED_RUNTIME_APPROVAL_STATE,
        source_module=source,
        atom_initializer=atom,
        checkpoints=tuple(checkpoints),
        asset_bundle_sha256=bundle,
    )


def load_manifest(path: Path, expected_sha256: str) -> Tuple[List[Dict[str, str]], str]:
    path = path if path.is_absolute() else (Path.cwd() / path).absolute()
    expected_sha256 = _strict_hash(expected_sha256, "manifest sha256")
    payload, observed_size, observed_sha256 = read_regular_file_bytes(path)
    if observed_size <= 0 or observed_sha256 != expected_sha256:
        raise MofclassifierRunError(
            f"manifest hash mismatch: expected {expected_sha256}, observed {observed_sha256}"
        )
    try:
        text = payload.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise MofclassifierRunError(f"manifest is not valid UTF-8: {path}: {exc}") from exc
    try:
        reader = csv.reader(io.StringIO(text, newline=""))
        try:
            header = next(reader)
        except StopIteration as exc:
            raise MofclassifierRunError("manifest is empty") from exc
        if tuple(header) != MANIFEST_FIELDS:
            raise MofclassifierRunError(f"manifest header mismatch: {header!r}")
        rows: List[Dict[str, str]] = []
        for line_number, values in enumerate(reader, start=2):
            if len(values) != len(header):
                raise MofclassifierRunError(f"manifest row width mismatch at line {line_number}")
            rows.append(dict(zip(header, values)))
    except csv.Error as exc:
        raise MofclassifierRunError(f"invalid manifest CSV {path}: {exc}") from exc
    seen_uid = set()
    seen_id = set()
    seen_path = set()
    for expected_index, row in enumerate(rows):
        if row["manifest_schema_version"] != MANIFEST_SCHEMA_VERSION:
            raise MofclassifierRunError(f"manifest schema mismatch at row {expected_index}")
        if row["row_index"] != str(expected_index):
            raise MofclassifierRunError(f"non-contiguous row_index at position {expected_index}")
        if not UID_RE.fullmatch(row["persistent_uid"]):
            raise MofclassifierRunError(f"invalid persistent_uid at row {expected_index}")
        try:
            parse_core_id(row["structure_id"])
        except ValueError as exc:
            raise MofclassifierRunError(f"invalid structure_id at row {expected_index}") from exc
        cif_path = Path(row["canonical_cif_path"])
        if not cif_path.is_absolute() or str(cif_path) != row["canonical_cif_path"]:
            raise MofclassifierRunError(f"canonical CIF path is not exact absolute syntax at row {expected_index}")
        _canonical_positive_int(row["canonical_cif_size_bytes"], "canonical CIF size")
        _strict_hash(row["canonical_cif_sha256"], "canonical CIF sha256")
        folded = row["structure_id"].casefold()
        if row["persistent_uid"] in seen_uid or folded in seen_id or row["canonical_cif_path"] in seen_path:
            raise MofclassifierRunError(f"manifest identity/path collision at row {expected_index}")
        seen_uid.add(row["persistent_uid"])
        seen_id.add(folded)
        seen_path.add(row["canonical_cif_path"])
    return rows, observed_sha256


def verify_canonical_input(row: Mapping[str, str]) -> Tuple[int, str]:
    path = Path(row["canonical_cif_path"])
    size, digest = observe_regular_file(path)
    expected_size = int(row["canonical_cif_size_bytes"])
    expected_digest = row["canonical_cif_sha256"]
    if size != expected_size or digest != expected_digest:
        raise MofclassifierRunError(
            f"canonical CIF binding mismatch for {row['structure_id']}: expected "
            f"{expected_size}/{expected_digest}, observed {size}/{digest}"
        )
    return size, digest


class PinnedMofclassifierBackend:
    """Thin CPU-only adapter around the hash-verified upstream CLscore source."""

    def __init__(self, contract: MethodContract):
        try:
            distribution = importlib.metadata.distribution(contract.package_distribution)
        except importlib.metadata.PackageNotFoundError as exc:
            raise MofclassifierRunError(
                f"distribution {contract.package_distribution!r} is not installed in this Python"
            ) from exc
        if distribution.version != contract.package_version:
            raise MofclassifierRunError(
                f"MOFClassifier version mismatch: expected {contract.package_version}, "
                f"observed {distribution.version}"
            )
        located_root = Path(distribution.locate_file("MOFClassifier")).resolve()
        if located_root != contract.package_root:
            raise MofclassifierRunError(
                f"installed distribution root {located_root} differs from pinned root {contract.package_root}"
            )
        # CLscore's import-time fallback downloads missing qsp/h models.  Refuse
        # import unless those directories already exist and are non-empty.
        for name in ("models", "models_qsp", "models_h"):
            directory = contract.package_root / name
            if not directory.is_dir() or not any(directory.iterdir()):
                raise MofclassifierRunError(
                    f"refusing CLscore import because it could mutate/download into {directory}"
                )
        source_payload, source_size, source_sha256 = read_regular_file_bytes(
            contract.source_module.path
        )
        if (
            source_size != contract.source_module.size_bytes
            or source_sha256 != contract.source_module.sha256
        ):
            raise MofclassifierRunError("CLscore source changed between contract validation and import")
        thread_environment = {
            "OMP_NUM_THREADS": str(contract.cpu_intraop_threads),
            "MKL_NUM_THREADS": str(contract.cpu_intraop_threads),
            "OPENBLAS_NUM_THREADS": str(contract.cpu_intraop_threads),
            "NUMEXPR_NUM_THREADS": str(contract.cpu_intraop_threads),
        }
        os.environ.update(thread_environment)
        module = types.ModuleType("_coremof_pinned_mofclassifier_clscore_v1")
        module.__file__ = str(contract.source_module.path)
        module.__package__ = ""
        exec(
            compile(source_payload, str(contract.source_module.path), "exec"),
            module.__dict__,
        )
        self.module = module
        self.torch = module.torch
        self.contract = contract
        self.torch.set_num_threads(contract.cpu_intraop_threads)
        try:
            self.torch.set_num_interop_threads(contract.cpu_interop_threads)
        except RuntimeError as exc:
            if self.torch.get_num_interop_threads() != contract.cpu_interop_threads:
                raise MofclassifierRunError(
                    "PyTorch inter-op thread count was already fixed to a conflicting value"
                ) from exc
        random.seed(contract.random_seed)
        try:
            module.np.random.seed(contract.random_seed)
        except AttributeError as exc:
            raise MofclassifierRunError("pinned CLscore source does not expose NumPy as np") from exc
        self.torch.manual_seed(contract.random_seed)
        self.torch.use_deterministic_algorithms(contract.deterministic_algorithms)
        observed_intraop = self.torch.get_num_threads()
        observed_interop = self.torch.get_num_interop_threads()
        observed_deterministic = self.torch.are_deterministic_algorithms_enabled()
        if (
            observed_intraop != contract.cpu_intraop_threads
            or observed_interop != contract.cpu_interop_threads
            or observed_deterministic is not contract.deterministic_algorithms
        ):
            raise MofclassifierRunError(
                "PyTorch deterministic CPU controls differ from the pinned method contract"
            )
        self.runtime = {
            "python": platform.python_version(),
            "torch": str(self.torch.__version__),
            "numpy": importlib.metadata.version("numpy"),
            "ase": importlib.metadata.version("ase"),
            "pymatgen": importlib.metadata.version("pymatgen"),
            "mofclassifier": distribution.version,
            "execution_device": "cpu",
            "cpu_intraop_threads": str(observed_intraop),
            "cpu_interop_threads": str(observed_interop),
            "random_seed": str(contract.random_seed),
            "deterministic_algorithms": str(observed_deterministic).lower(),
            "omp_num_threads": thread_environment["OMP_NUM_THREADS"],
            "mkl_num_threads": thread_environment["MKL_NUM_THREADS"],
            "openblas_num_threads": thread_environment["OPENBLAS_NUM_THREADS"],
            "numexpr_num_threads": thread_environment["NUMEXPR_NUM_THREADS"],
            "runtime_approval_state": contract.runtime_approval_state,
        }

    def preprocess(self, private_cif: Path, atom_initializer: Path, num_workers: int) -> Any:
        if num_workers != 0:
            raise MofclassifierRunError("internal error: num_workers is not zero")
        preload = self.module.preprocess(
            root_cif=str(private_cif), atom_init_file=str(atom_initializer)
        )
        loader = self.torch.utils.data.DataLoader(
            [preload],
            batch_size=1,
            shuffle=False,
            num_workers=0,
            collate_fn=self.module.collate_pool,
            pin_memory=False,
        )
        try:
            graph, identifiers = next(iter(loader))
        except StopIteration as exc:
            raise MofclassifierRunError("MOFClassifier produced no graph batch") from exc
        if len(identifiers) != 1:
            raise MofclassifierRunError("MOFClassifier graph batch does not contain exactly one CIF")
        return graph

    def load_model(self, checkpoint: Asset, reference_graph: Any) -> Any:
        payload, size, digest = read_regular_file_bytes(checkpoint.path)
        if size != checkpoint.size_bytes or digest != checkpoint.sha256:
            raise MofclassifierRunError(
                f"bag {checkpoint.bag_index} changed between contract validation and model load"
            )
        checkpoint_value = self.torch.load(
            io.BytesIO(payload), weights_only=False, map_location=self.torch.device("cpu")
        )
        if not isinstance(checkpoint_value, Mapping):
            raise MofclassifierRunError(f"bag {checkpoint.bag_index} checkpoint is not a mapping")
        required = {"args", "state_dict", "normalizer"}
        if not required.issubset(checkpoint_value):
            raise MofclassifierRunError(
                f"bag {checkpoint.bag_index} checkpoint lacks {sorted(required-set(checkpoint_value))}"
            )
        args = checkpoint_value["args"]
        if not isinstance(args, Mapping):
            raise MofclassifierRunError(f"bag {checkpoint.bag_index} model args are not a mapping")
        for key in ("atom_fea_len", "n_conv", "h_fea_len", "n_h"):
            if key not in args:
                raise MofclassifierRunError(f"bag {checkpoint.bag_index} lacks model arg {key}")
        model = self.module.CrystalGraphConvNet(
            reference_graph[0].shape[-1],
            reference_graph[1].shape[-1],
            atom_fea_len=args["atom_fea_len"],
            n_conv=args["n_conv"],
            h_fea_len=args["h_fea_len"],
            n_h=args["n_h"],
            classification=True,
        )
        model.load_state_dict(checkpoint_value["state_dict"])
        model.eval()
        return model

    def predict(self, model: Any, graph: Any) -> float:
        with self.torch.no_grad():
            output = model(graph[0], graph[1], graph[2], graph[3])
            probability = self.torch.exp(output.detach().cpu())
        if tuple(probability.shape) != (1, 2):
            raise MofclassifierRunError(f"unexpected probability tensor shape: {tuple(probability.shape)}")
        score = float(probability[0, self.contract.positive_class_index].item())
        if not math.isfinite(score) or score < 0.0 or score > 1.0:
            raise MofclassifierRunError(f"invalid MOFClassifier probability: {score!r}")
        return score


@dataclass
class PendingRecord:
    row: Dict[str, str]
    output_path: Path
    base: Dict[str, Any]
    graph: Any
    scores: List[float]
    error_stage: str = ""
    error_type: str = ""
    error_message: str = ""
    failed_bag_index: Optional[int] = None


def _record_output_path(output_root: Path, row: Mapping[str, str]) -> Path:
    return output_root / "records" / f"{int(row['row_index']):08d}_{row['structure_id']}.json"


def _base_record(
    row: Mapping[str, str],
    manifest_path: Path,
    manifest_sha256: str,
    contract: MethodContract,
    runtime: Mapping[str, str],
) -> Dict[str, Any]:
    return {
        "record_schema_version": RECORD_SCHEMA_VERSION,
        "runner_tool_version": TOOL_VERSION,
        "checker_name": "MOFClassifier",
        "method_id": contract.method_id,
        "fresh_model_inference": True,
        "mechanistic_hard_fail_eligible": False,
        "method_config_path": str(contract.config_path),
        "method_config_sha256": contract.config_sha256,
        "asset_bundle_sha256": contract.asset_bundle_sha256,
        "package_distribution": contract.package_distribution,
        "package_version": contract.package_version,
        "package_root": str(contract.package_root),
        "source_module_path": str(contract.source_module.path),
        "source_module_sha256": contract.source_module.sha256,
        "atom_initializer_path": str(contract.atom_initializer.path),
        "atom_initializer_sha256": contract.atom_initializer.sha256,
        "model_family": contract.model_family,
        "checkpoint_bag_count": len(contract.checkpoints),
        "execution_device": "cpu",
        "data_loader_num_workers": 0,
        "data_loader_shuffle": False,
        "data_loader_pin_memory": False,
        "positive_class_index": contract.positive_class_index,
        "pass_threshold": contract.threshold,
        "mean_aggregation": contract.mean_aggregation,
        "maximum_rows_per_invocation": contract.maximum_rows_per_invocation,
        "cpu_intraop_threads": contract.cpu_intraop_threads,
        "cpu_interop_threads": contract.cpu_interop_threads,
        "random_seed": contract.random_seed,
        "deterministic_algorithms": contract.deterministic_algorithms,
        "runtime_approval_state": contract.runtime_approval_state,
        "runtime": dict(runtime),
        "runtime_fingerprint_sha256": canonical_json_sha256(dict(runtime)),
        "source_manifest_path": str(manifest_path),
        "source_manifest_sha256": manifest_sha256,
        "manifest_row_index": int(row["row_index"]),
        "persistent_uid": row["persistent_uid"],
        "structure_id": row["structure_id"],
        "canonical_cif_path": row["canonical_cif_path"],
        "canonical_cif_size_bytes": int(row["canonical_cif_size_bytes"]),
        "canonical_cif_sha256": row["canonical_cif_sha256"],
        "canonical_pre_run_sha256": row["canonical_cif_sha256"],
        "canonical_post_run_sha256": "",
        "canonical_unchanged": False,
        "private_copy_used": False,
        "private_copy_initial_sha256": "",
        "private_copy_model_input_sha256": "",
        "private_copy_rewritten": False,
        "execution_status": "",
        "completed_bag_count": 0,
        "bag_scores": [],
        "mean_score": None,
        "operational_vote": None,
        "error_stage": "",
        "error_type": "",
        "error_message": "",
        "failed_bag_index": None,
    }


def _exact_static_record_fields(expected: Mapping[str, Any]) -> Tuple[str, ...]:
    dynamic = {
        "canonical_post_run_sha256",
        "canonical_unchanged",
        "private_copy_used",
        "private_copy_initial_sha256",
        "private_copy_model_input_sha256",
        "private_copy_rewritten",
        "execution_status",
        "completed_bag_count",
        "bag_scores",
        "mean_score",
        "operational_vote",
        "error_stage",
        "error_type",
        "error_message",
        "failed_bag_index",
    }
    return tuple(field for field in RECORD_FIELDS if field not in dynamic)


def validate_existing_record(path: Path, expected: Mapping[str, Any], retry_errors: bool) -> str:
    value = load_json(path)
    if not isinstance(value, Mapping) or set(value) != set(RECORD_FIELDS):
        raise MofclassifierRunError(f"existing record schema mismatch: {path}")
    for field in _exact_static_record_fields(expected):
        if value[field] != expected[field]:
            raise MofclassifierRunError(f"existing record {field} mismatch: {path}")
    if value["canonical_post_run_sha256"] != expected["canonical_cif_sha256"]:
        raise MofclassifierRunError(f"existing record lacks an exact post-run CIF hash: {path}")
    if value["canonical_unchanged"] is not True:
        raise MofclassifierRunError(f"existing record did not preserve its canonical CIF: {path}")
    if value["private_copy_used"] is True:
        if value["private_copy_initial_sha256"] != expected["canonical_cif_sha256"]:
            raise MofclassifierRunError(f"existing record private-copy hash mismatch: {path}")
        if value["private_copy_model_input_sha256"]:
            _strict_hash(value["private_copy_model_input_sha256"], "private model-input hash")
    elif value["private_copy_used"] is False:
        if (
            value["private_copy_initial_sha256"]
            or value["private_copy_model_input_sha256"]
            or value["private_copy_rewritten"] is not False
        ):
            raise MofclassifierRunError(f"existing record has inconsistent private-copy fields: {path}")
    else:
        raise MofclassifierRunError(f"existing record private_copy_used is not boolean: {path}")
    status = value["execution_status"]
    if status == "SUCCESS":
        if value["private_copy_used"] is not True or not value["private_copy_model_input_sha256"]:
            raise MofclassifierRunError(f"existing success record lacks a complete private-copy proof: {path}")
        scores = value["bag_scores"]
        if not isinstance(scores, list) or len(scores) != EXPECTED_BAG_COUNT:
            raise MofclassifierRunError(f"existing success record lacks 100 bag scores: {path}")
        for score in scores:
            if isinstance(score, bool) or not isinstance(score, (int, float)):
                raise MofclassifierRunError(f"existing record has a non-numeric bag score: {path}")
            if not math.isfinite(float(score)) or not 0.0 <= float(score) <= 1.0:
                raise MofclassifierRunError(f"existing record has an invalid bag score: {path}")
        mean = math.fsum(float(score) for score in scores) / EXPECTED_BAG_COUNT
        if value["completed_bag_count"] != EXPECTED_BAG_COUNT or value["mean_score"] != mean:
            raise MofclassifierRunError(f"existing record mean/completed count mismatch: {path}")
        vote = "PASS" if mean >= EXPECTED_THRESHOLD else "FAIL"
        if value["operational_vote"] != vote:
            raise MofclassifierRunError(f"existing record vote mismatch: {path}")
        if any(value[field] not in ("", None) for field in ("error_stage", "error_type", "error_message", "failed_bag_index")):
            raise MofclassifierRunError(f"existing success record carries an error: {path}")
        return "RESUME_SUCCESS"
    if status == "ERROR":
        # Records produced before the 2026-07-24 contract alignment used an
        # empty string for an explicit non-vote.  Accept those records for
        # resume compatibility, while all newly written error records use the
        # unambiguous JSON null representation.
        if value["operational_vote"] not in ("", None) or value["mean_score"] is not None:
            raise MofclassifierRunError(f"existing error record has a nonblank vote/mean: {path}")
        scores = value["bag_scores"]
        if not isinstance(scores, list):
            raise MofclassifierRunError(f"existing error record bag_scores is not an array: {path}")
        for score in scores:
            if isinstance(score, bool) or not isinstance(score, (int, float)):
                raise MofclassifierRunError(
                    f"existing error record has a non-numeric partial bag score: {path}"
                )
            if not math.isfinite(float(score)) or not 0.0 <= float(score) <= 1.0:
                raise MofclassifierRunError(
                    f"existing error record has an invalid partial bag score: {path}"
                )
        count = value["completed_bag_count"]
        if (
            not isinstance(count, int)
            or isinstance(count, bool)
            or not 0 <= count <= EXPECTED_BAG_COUNT
            or count != len(scores)
        ):
            raise MofclassifierRunError(f"existing error record completed count is invalid: {path}")
        stage = value["error_stage"]
        if stage not in {
            "PRIVATE_COPY_OR_PREPROCESS",
            "MODEL_LOAD",
            "INFERENCE",
            "CANONICAL_POSTCHECK",
        } or not isinstance(value["error_type"], str) or not value["error_type"]:
            raise MofclassifierRunError(f"existing error record lacks error provenance: {path}")
        if not isinstance(value["error_message"], str):
            raise MofclassifierRunError(f"existing error record error_message is not text: {path}")
        failed_bag = value["failed_bag_index"]
        if stage in {"MODEL_LOAD", "INFERENCE"}:
            if (
                not isinstance(failed_bag, int)
                or isinstance(failed_bag, bool)
                or not 1 <= failed_bag <= EXPECTED_BAG_COUNT
                or count != failed_bag - 1
                or value["private_copy_used"] is not True
                or not value["private_copy_model_input_sha256"]
            ):
                raise MofclassifierRunError(
                    f"existing error record has inconsistent failed-bag provenance: {path}"
                )
        elif failed_bag is not None:
            raise MofclassifierRunError(
                f"existing non-bag error record carries a failed bag index: {path}"
            )
        if stage == "PRIVATE_COPY_OR_PREPROCESS" and count != 0:
            raise MofclassifierRunError(
                f"existing preprocess error record carries completed bag scores: {path}"
            )
        return "RETRY_ERROR" if retry_errors else "RESUME_ERROR"
    raise MofclassifierRunError(f"existing record has unsupported status {status!r}: {path}")


def _prepare_one(
    row: Dict[str, str],
    output_path: Path,
    base: Dict[str, Any],
    backend: Any,
    contract: MethodContract,
    private_root: Path,
) -> PendingRecord:
    private_root.mkdir(parents=True, exist_ok=True)
    pending = PendingRecord(row, output_path, base, graph=None, scores=[])
    try:
        with tempfile.TemporaryDirectory(
            prefix=f"mofclassifier_{int(row['row_index']):08d}_", dir=str(private_root)
        ) as directory:
            private_cif = Path(directory) / f"{row['structure_id']}.cif"
            private_atom = Path(directory) / "atom_init.json"
            shutil.copyfile(row["canonical_cif_path"], private_cif)
            shutil.copyfile(contract.atom_initializer.path, private_atom)
            private_size, private_initial = observe_regular_file(private_cif.resolve())
            if private_size != int(row["canonical_cif_size_bytes"]) or private_initial != row["canonical_cif_sha256"]:
                raise MofclassifierRunError("private copy does not match the canonical CIF binding")
            base["private_copy_used"] = True
            base["private_copy_initial_sha256"] = private_initial
            atom_size, atom_hash = observe_regular_file(private_atom.resolve())
            if (
                atom_size != contract.atom_initializer.size_bytes
                or atom_hash != contract.atom_initializer.sha256
            ):
                raise MofclassifierRunError("private atom initializer does not match its pinned asset")
            # Detect a source-file race before giving any path to upstream code.
            verify_canonical_input(row)
            try:
                pending.graph = backend.preprocess(
                    private_cif.resolve(), private_atom.resolve(), num_workers=0
                )
            finally:
                if private_cif.is_file():
                    _private_size, model_input_hash = observe_regular_file(private_cif.resolve())
                    base["private_copy_model_input_sha256"] = model_input_hash
                    base["private_copy_rewritten"] = model_input_hash != private_initial
        if pending.graph is None:
            raise MofclassifierRunError("preprocessing returned no graph")
    except Exception as exc:  # per-record execution defect; never a chemical FAIL
        pending.error_stage = "PRIVATE_COPY_OR_PREPROCESS"
        pending.error_type = type(exc).__name__
        pending.error_message = str(exc)
    return pending


def _finalize_record(pending: PendingRecord, contract: MethodContract) -> Dict[str, Any]:
    record = pending.base
    try:
        _size, post_hash = verify_canonical_input(pending.row)
        record["canonical_post_run_sha256"] = post_hash
        record["canonical_unchanged"] = post_hash == pending.row["canonical_cif_sha256"]
    except Exception as exc:
        record["canonical_post_run_sha256"] = ""
        record["canonical_unchanged"] = False
        pending.error_stage = "CANONICAL_POSTCHECK"
        pending.error_type = type(exc).__name__
        pending.error_message = str(exc)
        pending.failed_bag_index = None
    if pending.error_type:
        if len(pending.scores) > EXPECTED_BAG_COUNT:
            raise MofclassifierRunError("internal error: error record has too many bag scores")
        for score in pending.scores:
            if (
                isinstance(score, bool)
                or not isinstance(score, (int, float))
                or not math.isfinite(float(score))
                or not 0.0 <= float(score) <= 1.0
            ):
                raise MofclassifierRunError(
                    "internal error: error record has an invalid partial bag score"
                )
        record["execution_status"] = "ERROR"
        record["completed_bag_count"] = len(pending.scores)
        record["bag_scores"] = list(pending.scores)
        record["mean_score"] = None
        record["operational_vote"] = None
        record["error_stage"] = pending.error_stage
        record["error_type"] = pending.error_type
        record["error_message"] = pending.error_message[:4000]
        record["failed_bag_index"] = pending.failed_bag_index
    else:
        if len(pending.scores) != len(contract.checkpoints):
            raise MofclassifierRunError("internal error: successful record lacks all bag scores")
        mean = math.fsum(pending.scores) / len(contract.checkpoints)
        record["execution_status"] = "SUCCESS"
        record["completed_bag_count"] = len(pending.scores)
        record["bag_scores"] = pending.scores
        record["mean_score"] = mean
        record["operational_vote"] = "PASS" if mean >= contract.threshold else "FAIL"
        record["error_stage"] = ""
        record["error_type"] = ""
        record["error_message"] = ""
        record["failed_bag_index"] = None
    if tuple(record) != RECORD_FIELDS:
        raise MofclassifierRunError("internal record schema/order drift")
    atomic_write_json(pending.output_path, record)
    return record


def normalize_selection(
    selected_indices: Optional[Sequence[int]],
    row_count: int,
    maximum_rows_per_invocation: int,
) -> List[int]:
    if selected_indices is None:
        raise MofclassifierRunError("explicit row selection is required")
    if any(not isinstance(value, int) or isinstance(value, bool) for value in selected_indices):
        raise MofclassifierRunError("row selection values must be integers")
    selected = sorted(set(selected_indices))
    if not selected:
        raise MofclassifierRunError("row selection is empty")
    if len(selected) > maximum_rows_per_invocation:
        raise MofclassifierRunError(
            f"row selection exceeds the pinned {maximum_rows_per_invocation}-row invocation limit"
        )
    if selected[0] < 0 or selected[-1] >= row_count:
        raise MofclassifierRunError("row selection is outside the manifest")
    return selected


def run(
    manifest_path: Path,
    manifest_sha256: str,
    config_path: Path,
    output_root: Path,
    selected_indices: Optional[Sequence[int]],
    *,
    retry_errors: bool = False,
    package_root_override: Optional[Path] = None,
    private_root: Optional[Path] = None,
    backend_factory: Callable[[MethodContract], Any] = PinnedMofclassifierBackend,
) -> Dict[str, Any]:
    contract = load_method_contract(config_path, package_root_override)
    rows, observed_manifest_sha256 = load_manifest(manifest_path, manifest_sha256)
    manifest_path = manifest_path.resolve()
    output_root = output_root.resolve()
    selected = normalize_selection(
        selected_indices, len(rows), contract.maximum_rows_per_invocation
    )
    backend = backend_factory(contract)
    runtime = dict(backend.runtime)
    if not runtime:
        raise MofclassifierRunError("backend runtime provenance is empty")
    expected_runtime_controls = {
        "execution_device": "cpu",
        "cpu_intraop_threads": str(contract.cpu_intraop_threads),
        "cpu_interop_threads": str(contract.cpu_interop_threads),
        "random_seed": str(contract.random_seed),
        "deterministic_algorithms": str(contract.deterministic_algorithms).lower(),
        "omp_num_threads": str(contract.cpu_intraop_threads),
        "mkl_num_threads": str(contract.cpu_intraop_threads),
        "openblas_num_threads": str(contract.cpu_intraop_threads),
        "numexpr_num_threads": str(contract.cpu_intraop_threads),
        "runtime_approval_state": contract.runtime_approval_state,
    }
    runtime_mismatches = {
        key: (runtime.get(key), expected)
        for key, expected in expected_runtime_controls.items()
        if runtime.get(key) != expected
    }
    if runtime_mismatches:
        raise MofclassifierRunError(
            f"backend runtime controls differ from the pinned method contract: {runtime_mismatches}"
        )
    private_root = (private_root or (output_root / ".private_work")).resolve()
    pending: List[PendingRecord] = []
    resumed_success = 0
    resumed_error = 0
    for index in selected:
        row = rows[index]
        verify_canonical_input(row)
        base = _base_record(row, manifest_path, observed_manifest_sha256, contract, runtime)
        output_path = _record_output_path(output_root, row)
        if output_path.exists() or output_path.is_symlink():
            disposition = validate_existing_record(output_path, base, retry_errors)
            if disposition == "RESUME_SUCCESS":
                resumed_success += 1
                continue
            if disposition == "RESUME_ERROR":
                resumed_error += 1
                continue
        pending.append(_prepare_one(row, output_path, base, backend, contract, private_root))

    active = [item for item in pending if not item.error_type]
    for checkpoint in contract.checkpoints:
        if not active:
            break
        try:
            model = backend.load_model(checkpoint, active[0].graph)
        except Exception as exc:
            for item in active:
                item.error_stage = "MODEL_LOAD"
                item.error_type = type(exc).__name__
                item.error_message = str(exc)
                item.failed_bag_index = checkpoint.bag_index
            active = []
            break
        still_active: List[PendingRecord] = []
        for item in active:
            try:
                score = backend.predict(model, item.graph)
                if isinstance(score, bool) or not isinstance(score, (int, float)):
                    raise MofclassifierRunError(f"backend returned non-numeric score {score!r}")
                score = float(score)
                if not math.isfinite(score) or not 0.0 <= score <= 1.0:
                    raise MofclassifierRunError(f"backend returned invalid score {score!r}")
                item.scores.append(score)
                still_active.append(item)
            except Exception as exc:
                item.error_stage = "INFERENCE"
                item.error_type = type(exc).__name__
                item.error_message = str(exc)
                item.failed_bag_index = checkpoint.bag_index
        active = still_active
        del model

    statuses: Dict[str, int] = {"SUCCESS": 0, "ERROR": 0}
    for item in pending:
        record = _finalize_record(item, contract)
        statuses[record["execution_status"]] += 1
    return {
        "run_summary_schema_version": "1.0",
        "runner_tool_version": TOOL_VERSION,
        "method_id": contract.method_id,
        "method_config_path": str(contract.config_path),
        "method_config_sha256": contract.config_sha256,
        "asset_bundle_sha256": contract.asset_bundle_sha256,
        "source_manifest_path": str(manifest_path),
        "source_manifest_sha256": observed_manifest_sha256,
        "manifest_row_count": len(rows),
        "selected_row_indices": selected,
        "selected_row_count": len(selected),
        "selection_covers_full_manifest": selected == list(range(len(rows))),
        "maximum_rows_per_invocation": contract.maximum_rows_per_invocation,
        "runtime_approval_state": contract.runtime_approval_state,
        "new_record_counts": statuses,
        "resumed_success_count": resumed_success,
        "resumed_error_count": resumed_error,
        "output_root": str(output_root),
    }


def _parse_row_range(value: str) -> List[int]:
    match = re.fullmatch(r"(0|[1-9][0-9]*):(0|[1-9][0-9]*)", value)
    if not match:
        raise argparse.ArgumentTypeError("row range must be START:STOP with non-negative integers")
    start, stop = int(match.group(1)), int(match.group(2))
    if stop <= start:
        raise argparse.ArgumentTypeError("row range STOP must exceed START")
    return list(range(start, stop))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--package-root", type=Path)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--private-root", type=Path)
    parser.add_argument("--retry-errors", action="store_true")
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--row-index", type=int, action="append")
    selection.add_argument("--row-range", type=_parse_row_range)
    parser.add_argument("--summary-output", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    selected = args.row_index if args.row_index is not None else args.row_range
    summary = run(
        args.manifest,
        args.manifest_sha256,
        args.config,
        args.output_root,
        selected,
        retry_errors=args.retry_errors,
        package_root_override=args.package_root,
        private_root=args.private_root,
    )
    if args.summary_output:
        atomic_write_json(args.summary_output.resolve(), summary)
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
