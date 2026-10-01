"""Reproduce the recorded v4.1 source-CIF curation without altering inputs.

This explicit historical-method route is separate from ``curate.clean``. It
does not promise to recreate the historical base release or the separate COD
workstation implementation. No checker labels or release promotion are made.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import tempfile

from . import _release_curation_v41_protocol as protocol
from . import _release_curation_worker as worker
from ._transactions import publish_directory


PROFILE = "coremof-curation-v4.1-20260716"


class ReleaseCurationError(RuntimeError):
    """The explicit historical curation input/runtime/output contract failed."""


def _stop(process):
    for sig in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.killpg(process.pid, sig)
        except ProcessLookupError:
            pass
        try:
            process.wait(timeout=3)
        except subprocess.TimeoutExpired:
            pass
    process.wait()


def _json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def _relative_records(value, root):
    if isinstance(value, dict):
        return {key: _relative_records(item, root) for key, item in value.items()}
    if isinstance(value, list):
        return [_relative_records(item, root) for item in value]
    if isinstance(value, str) and value.startswith(str(root) + os.sep):
        return str(Path(value).relative_to(root))
    return value


def _validate(record, source_id, source_sha, skip_charges, root):
    if (record.get("refcode") != source_id or record.get("source_sha256") != source_sha
            or record.get("pipeline_version") != protocol.PIPELINE_VERSION
            or record.get("status") not in {"complete", "review", "failed"}):
        raise ReleaseCurationError("Curation record has an unexpected identity, method or status")
    if record.get("release_eligible") not in (False, None):
        raise ReleaseCurationError("Curation alone cannot certify release eligibility")
    if record["status"] == "failed":
        if not record.get("errors"):
            raise ReleaseCurationError("Failed curation has no diagnostic")
        return
    if record["status"] == "complete":
        if (record.get("review_reasons") or not record.get("invariants", {}).get("all_required_invariants_passed")
                or not record.get("deliverables")):
            raise ReleaseCurationError("Complete curation lacks passing chemistry/mapping checks")
    elif not record.get("review_reasons"):
        raise ReleaseCurationError("Review curation has no review reason")
    if record.get("category") == "ION":
        if (record.get("curation", {}).get("asr", {}).get("status") != "skipped"
                or [item.get("roles") for item in record["deliverables"]] != [["ION_FSR"]]):
            raise ReleaseCurationError("Ionic curation must retain ION_FSR and skip ASR")
    for item in record.get("deliverables", []):
        if item.get("release_eligible") is not False:
            raise ReleaseCurationError("A curation candidate was marked release eligible")
        pairs = [(item["uncharged_cif"], item["uncharged_sha256"])]
        charge = item.get("pacman", {})
        if skip_charges or record["status"] == "review":
            if item.get("curation_stage_eligible") is not False or charge.get("status") == "complete":
                raise ReleaseCurationError("Uncharged/review candidate was marked ready")
        elif record["status"] == "complete":
            if (charge.get("status") != "complete" or charge.get("net_charge_check") != "passed"
                    or item.get("curation_stage_eligible") is not True):
                raise ReleaseCurationError("Complete candidate lacks validated charges")
            pairs.append((charge["charged_cif"], charge["charged_sha256"]))
        for name, digest in pairs:
            path = Path(name).resolve()
            if root.resolve() not in path.parents:
                raise ReleaseCurationError("Curation output escaped its private directory")
            worker.check(path, digest)


def curate_cif_v41(cif_path, *, source_id, output_dir, python, runtime_coremof_root,
                   skip_charges=False, timeout_seconds=600):
    """Run the recorded v4.1 method on one private copy of a raw source CIF.

    ``source_id`` is the source record name, before public ASR/FSR/ION IDs are
    assigned. The external runtime must have the exact recorded versions and
    frozen CoREMOF source. Missing models are errors, never downloaded. The
    destination must not exist. REVIEW/ERROR/TIMEOUT never mean NCR. Setting
    ``skip_charges=True`` yields inspection-only, non-eligible candidates.
    """
    if os.name != "posix":
        raise ReleaseCurationError("The historical curation runtime requires POSIX")
    if not isinstance(source_id, str) or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,99}", source_id) is None:
        raise ValueError("source_id must be a safe source record name")
    if type(timeout_seconds) is not int or timeout_seconds < 1 or type(skip_charges) is not bool:
        raise ValueError("Use a positive integer timeout and a boolean skip_charges")
    destination = Path(output_dir).absolute()
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    if not destination.parent.is_dir():
        raise FileNotFoundError(destination.parent)
    python = worker.check(Path(python).absolute(), worker.PYTHON_SHA256)
    runtime_root = Path(runtime_coremof_root).resolve(strict=True)
    for name, digest in worker.SOURCE_HASHES.items():
        if name.startswith("CoREMOF."):
            worker.check(runtime_root.joinpath(*name.split(".")).with_suffix(".py"), digest)
    profile = worker.check(protocol.__file__, worker.PROTOCOL_SHA256)
    source = Path(cif_path).resolve(strict=True)
    data = source.read_bytes()
    source_sha = hashlib.sha256(data).hexdigest()
    with tempfile.TemporaryDirectory(prefix=".curation-v41-", dir=destination.parent) as temporary:
        root = Path(temporary)
        for name in ("input", "home", "tmp", "work", "output"):
            (root / name).mkdir(mode=0o700)
        output = root / "output"
        copied = root / "input" / (source_id + ".cif")
        copied.write_bytes(data)
        private_protocol = root / "protocol.py"
        private_protocol.write_bytes(profile.read_bytes())
        private_worker = root / "worker.py"
        private_worker.write_bytes(Path(worker.__file__).read_bytes())
        record_path = output / "protocol_record.json"
        request = {"input": str(copied), "input_sha256": source_sha,
                   "protocol": str(private_protocol), "output": str(output / "artifacts"),
                   "record": str(record_path), "runtime_receipt": str(output / "runtime.json"),
                   "skip_charges": skip_charges}
        request_path = root / "request.json"
        _json(request_path, request)
        env = {"PATH": str(python.parent) + ":/usr/bin:/bin", "HOME": str(root / "home"),
               "TMPDIR": str(root / "tmp"), "MPLCONFIGDIR": str(root / "tmp" / "matplotlib"),
               "PYTHONPATH": str(runtime_root), "PYTHONNOUSERSITE": "1", "PYTHONDONTWRITEBYTECODE": "1",
               "PYTHONHASHSEED": "0", "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1",
               "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
               "LANG": "C.UTF-8", "LC_ALL": "C.UTF-8", "TZ": "UTC",
               "LD_LIBRARY_PATH": str(python.parent.parent / "lib")}
        with (output / "stdout.txt").open("w") as stdout, (output / "stderr.txt").open("w") as stderr:
            process = subprocess.Popen([str(python), "-B", str(private_worker), str(request_path)],
                                       cwd=root / "work", env=env, stdout=stdout, stderr=stderr,
                                       start_new_session=True)
            timed_out = False
            try:
                process.wait(timeout=timeout_seconds)
            except subprocess.TimeoutExpired:
                timed_out = True
                _stop(process)
            except BaseException:
                _stop(process)
                raise
        worker.check(source, source_sha)
        worker.check(private_protocol, worker.PROTOCOL_SHA256)
        if process.returncode == 0 and record_path.is_file() and not timed_out:
            record = json.loads(record_path.read_text())
            _validate(record, source_id, source_sha, skip_charges, output)
            status = {"complete": "COMPLETE", "review": "REVIEW", "failed": "ERROR"}[record["status"]]
            record = _relative_records(record, output)
            # The source copy is not an output artifact. Point to the original
            # immutable input in the private receipt, not a vanished tmp path.
            record["source_cif"] = str(source)
            _json(record_path, record)
        else:
            status = "TIMEOUT" if timed_out else "ERROR"
            record = {"refcode": source_id, "source_sha256": source_sha, "status": "failed",
                      "errors": [{"type": status, "message": "See stderr.txt for the worker diagnostic"}],
                      "deliverables": [], "curation_stage_eligible": False, "release_eligible": False}
            _json(record_path, record)
        result = {"schema_version": "coremof-curation-v41-replay/1.0", "profile": PROFILE,
                  "source_id": source_id, "input_sha256": source_sha, "execution_status": status,
                  "category": record.get("category"), "deliverables": record.get("deliverables", []),
                  "curation_stage_eligible": record.get("curation_stage_eligible", False),
                  "release_eligible": False, "review_reasons": record.get("review_reasons", []),
                  "errors": record.get("errors", [])}
        receipt = {"profile": PROFILE, "protocol_sha256": worker.PROTOCOL_SHA256,
                   "worker_sha256": worker.sha(private_worker), "input_sha256": source_sha,
                   "original_cif_modified": False, "release_metadata_promoted": False,
                   "skip_charges": skip_charges, "timeout_seconds": timeout_seconds,
                   "returncode": process.returncode, "runtime_coremof_root": str(runtime_root),
                   "historical_full_runtime_byte_identity_proven": False,
                   "scope": "Recorded CSD v4.1 method, not base-release or COD-workstation equivalence"}
        _json(output / "record.json", result)
        _json(output / "receipt.json", receipt)
        publish_directory(output, destination, overwrite=False)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cif_path")
    for name in ("source-id", "output-dir", "python", "runtime-coremof-root"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--skip-charges", action="store_true")
    parser.add_argument("--timeout-seconds", type=int, default=600)
    result = curate_cif_v41(**vars(parser.parse_args(argv)))
    print(json.dumps({key: result[key] for key in ("source_id", "execution_status")}))
    return 0 if result["execution_status"] == "COMPLETE" else 1


if __name__ == "__main__":
    raise SystemExit(main())
