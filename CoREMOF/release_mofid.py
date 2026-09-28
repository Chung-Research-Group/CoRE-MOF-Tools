"""Explicit replay of the frozen CoRE-MOF release MOFid method.

This is separate from the general-purpose :mod:`CoREMOF.get_mofid` wrapper.
It needs the recorded external runtime and official node library, never
downloads dependencies, and never promotes results into a release.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import tempfile

from . import _release_mofid_protocol as protocol


PROTOCOL_SHA256 = "54481294fbb78e94aae9f0db4baaaa879ad321309b18aecf9fe51d3c00061bd9"
METHOD_SHA256 = "5f1d01bd45d3bde59c575672b1be82f8aaac29b586bf4c5128c5dba168356d32"
NODE_MANIFEST_SHA256 = "9f230f8f028907829c255c47f38326ba222bed4502154b737d7a408fb64f664f"
PROFILE = "coremof-release-mofid-20260730"


class ReleaseMOFidError(RuntimeError):
    """The frozen method, inputs or runtime could not be verified."""


def _checked_file(path, digest, size=None):
    path = Path(path)
    if not path.is_file() or (size is not None and path.stat().st_size != size):
        raise ReleaseMOFidError(f"Missing file or size mismatch: {path}")
    if protocol.sha256_file(path) != digest:
        raise ReleaseMOFidError(f"SHA-256 mismatch: {path}")
    return path


def _read_profile(method_manifest, node_manifest):
    _checked_file(protocol.__file__, PROTOCOL_SHA256)
    method = json.loads(_checked_file(method_manifest, METHOD_SHA256).read_text())
    _checked_file(node_manifest, NODE_MANIFEST_SHA256)
    # The fixed hashes bind every setting, including official archive order.
    # Do not accept an edited manifest by silently trusting its new hash.
    return method


def _safe_id(structure_id, variant):
    if variant not in {"ASR", "FSR", "ION"} or re.fullmatch(
        r"(?:ASR|FSR|ION)-(?:COD|CSD|SI)-(?:[0-9]{4}|UNKN)-[0-9]{4,}",
        structure_id,
    ) is None or not structure_id.startswith(variant + "-"):
        raise ValueError("Use a public CoRE-MOF structure ID and its matching variant")


def _runtime_check(request, method):
    """Check bytes and imports in the process that will calculate the record."""
    versions = {"python": ".".join(map(str, sys.version_info[:2]))}
    for name in protocol.PINNED_PYTHON_PACKAGES:
        if name != "python":
            versions[name] = importlib.metadata.version(name)
    if versions != protocol.PINNED_PYTHON_PACKAGES:
        raise ReleaseMOFidError(f"Pinned Python/package versions differ: {versions}")
    from mofid import paths
    from mofid import run_mofid

    source = Path(request["source_root"]).resolve()
    package = Path(request["mofid_site"]).resolve() / "mofid"
    if Path(paths.mofid_path).resolve() != source:
        raise ReleaseMOFidError("Imported MOFid points at another source tree")
    if Path(run_mofid.__file__).resolve().parent != package:
        raise ReleaseMOFidError("Imported MOFid is outside the requested installation")
    runtime = method["runtime_validation"]
    for name, receipt in runtime["runtime_files"].items():
        _checked_file(source / name, receipt["sha256"], receipt["size_bytes"])
    roots = {
        "imported_mofid_package": package,
        "openbabel_shared_libraries_and_plugins": source / "openbabel/build/lib",
        "openbabel_data": source / "openbabel/data",
    }
    for name, receipt in runtime["runtime_trees"].items():
        protocol.validate_tree_receipt({**receipt, "root": str(roots[name])})
    for command, expected in (
        (["java", "-version"], runtime["java"]),
        ([str(source / "openbabel/build/bin/obabel"), "-V"], runtime["openbabel"]),
    ):
        completed = subprocess.run(command, capture_output=True, text=True, timeout=20)
        if completed.returncode or expected not in completed.stdout + completed.stderr:
            raise ReleaseMOFidError(f"Runtime version differs: {command[0]}")
    # Check all 1,182 library files before the first scientific operation.
    protocol.load_node_manifest(
        Path(request["node_manifest"]), Path(request["node_root"]), verify_files=True
    )
    return versions


def _worker(request_path):
    request = json.loads(Path(request_path).read_text())
    row = request["row"]
    _safe_id(row["structure_id"], row["structure_variant"])
    method = _read_profile(request["method_manifest"], request["node_manifest"])
    versions = _runtime_check(request, method)
    # The parent creates this directory exclusively. The protocol's historical
    # cleanup can therefore affect only this attempt's disposable descendants.
    private = Path(request["private_root"])
    cif_path = private / "input" / (row["structure_id"] + ".cif")
    _checked_file(cif_path, row["cif_sha256"])
    record = protocol.calculate_record(
        row, cif_path, private / "work", Path(request["node_root"]),
        Path(request["node_manifest"]), request["input_sha256"],
        NODE_MANIFEST_SHA256, METHOD_SHA256,
    )
    _checked_file(cif_path, row["cif_sha256"])
    protocol.validate_result_record(record, row["structure_id"])
    protocol.atomic_write_json(private / "record.json", record)
    protocol.atomic_write_json(private / "runtime.json", {"versions": versions})


def _stop_own_group(process):
    """Terminate only the new worker session and any of its descendants."""
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        process.wait(timeout=3)
    except subprocess.TimeoutExpired:
        pass
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait()


def calculate_release_mofid(
    cif_path, *, structure_id, structure_variant, output_dir, python,
    method_manifest, node_manifest, node_root, source_root, pinned_site,
    mofid_site, library_paths=(), existing_mofid_v1=None, timeout_seconds=300,
):
    """Calculate one release-method record in an isolated, bounded process.

    ``output_dir`` must not exist. ``python`` is the pinned Python 3.9 runtime,
    not necessarily the caller's interpreter. External paths are explicit.
    Timeouts, decomposition errors and substantive v1 mismatches remain
    visible in the returned record. A completed record is not a promotion
    decision. The original CIF and existing output directories are untouched.
    """
    if os.name != "posix":
        raise ReleaseMOFidError("The frozen external runtime requires POSIX")
    _safe_id(structure_id, structure_variant)
    if isinstance(timeout_seconds, bool) or not math.isfinite(timeout_seconds) or timeout_seconds <= 0:
        raise ValueError("timeout_seconds must be finite and positive")
    _read_profile(method_manifest, node_manifest)
    source_cif = Path(cif_path).resolve(strict=True)
    cif_bytes = source_cif.read_bytes()
    cif_sha = hashlib.sha256(cif_bytes).hexdigest()
    destination = Path(output_dir).absolute()
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(f"Output already exists: {destination}")
    if not destination.parent.is_dir():
        raise FileNotFoundError(f"Create the output parent first: {destination.parent}")
    row = {
        "structure_id": structure_id, "structure_variant": structure_variant,
        "source_database": structure_id.split("-")[1], "source_id": "",
        "cif_file": source_cif.name, "cif_sha256": cif_sha,
        "existing_mofid_v1": existing_mofid_v1, "existing_mofid_v2": None,
    }
    input_bytes = (json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n").encode()
    input_sha = hashlib.sha256(input_bytes).hexdigest()
    with tempfile.TemporaryDirectory(prefix=".coremof-mofid-", dir=destination.parent) as temporary:
        private = Path(temporary)
        (private / "input").mkdir()
        # MOFid includes the CIF stem in its raw record suffix. Preserve the
        # public structure name even though the bytes live in an isolated copy.
        (private / "input" / (structure_id + ".cif")).write_bytes(cif_bytes)
        request = {
            "row": row, "private_root": str(private), "input_sha256": input_sha,
            **{name: str(Path(value).resolve(strict=True)) for name, value in (
                ("method_manifest", method_manifest), ("node_manifest", node_manifest),
                ("node_root", node_root), ("source_root", source_root),
                ("pinned_site", pinned_site), ("mofid_site", mofid_site),
            )},
        }
        protocol.atomic_write_json(private / "request.json", request)
        environment = os.environ.copy()
        environment.update({
            "PYTHONPATH": os.pathsep.join((str(Path(__file__).resolve().parent.parent),
                request["pinned_site"], request["mofid_site"])),
            "PYTHONNOUSERSITE": "1", "PYTHONDONTWRITEBYTECODE": "1",
            "JAVA_TOOL_OPTIONS": protocol.JAVA_TOOL_OPTIONS,
            "LD_LIBRARY_PATH": os.pathsep.join(
                [str(Path(path).resolve(strict=True)) for path in library_paths]
                + [str(Path(request["source_root"]) / "openbabel/build/lib")]
            ),
            "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
        })
        timed_out = False
        command = [str(Path(python).resolve(strict=True)), "-B", "-m",
                   "CoREMOF.release_mofid", "--worker", str(private / "request.json")]
        with (private / "worker.log").open("wb") as log:
            process = subprocess.Popen(command, cwd=private, env=environment,
                stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            try:
                process.wait(timeout=timeout_seconds)
            except subprocess.TimeoutExpired:
                timed_out = True
                _stop_own_group(process)
            except BaseException:
                _stop_own_group(process)
                raise
        _checked_file(source_cif, cif_sha)
        if timed_out:
            # Verification may itself have timed out. Do not assert that the
            # scientific calculation started or that runtime validation passed.
            record = protocol.timeout_record(row, cif_sha, input_sha,
                NODE_MANIFEST_SHA256, METHOD_SHA256, timeout_seconds)
            record["calculation"]["scope"] = "runtime_verification_or_calculation"
        elif process.returncode:
            message = (private / "worker.log").read_text(errors="replace")[-6000:]
            raise ReleaseMOFidError(f"Frozen runtime worker failed ({process.returncode}): {message}")
        else:
            record = json.loads((private / "record.json").read_text())
        protocol.validate_result_record(record, structure_id)
        payload = private / "publish"
        payload.mkdir(mode=0o700)
        protocol.atomic_write_json(payload / "record.json", record)
        (payload / "input.json").write_bytes(input_bytes)
        shutil.copyfile(private / "worker.log", payload / "worker.log")
        receipt = {
            "schema": "coremof-release-mofid-replay/1.0", "profile": PROFILE,
            "cif_sha256": cif_sha, "input_sha256": input_sha,
            "method_manifest_sha256": METHOD_SHA256,
            "node_library_manifest_sha256": NODE_MANIFEST_SHA256,
            "scientific_protocol_sha256": PROTOCOL_SHA256,
            "adapter_sha256": protocol.sha256_file(Path(__file__)),
            "runtime_verified": (private / "runtime.json").exists(),
            "timeout_seconds": timeout_seconds, "process_returncode": process.returncode,
            "record_sha256": protocol.sha256_file(payload / "record.json"),
            "release_promotion": False,
        }
        protocol.atomic_write_json(payload / "receipt.json", receipt)
        # Reserve the destination exclusively. Another process cannot be
        # overwritten even if it created this directory after the initial check.
        destination.mkdir(mode=0o700)
        try:
            for child in payload.iterdir():
                shutil.move(str(child), str(destination / child.name))
        except BaseException:
            shutil.rmtree(destination)
            raise
        return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", required=True, help=argparse.SUPPRESS)
    _worker(parser.parse_args().worker)
