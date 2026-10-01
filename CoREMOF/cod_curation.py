"""Run the recorded COD v4.1 curation on private copies, without downloads.

This is distinct from both legacy ``curate.clean`` and CSD ``curate_cif_v41``.
It preserves the recorded scientific method, not a promise of bitwise equality
with workstation coordinates or charges. No release CIF or metadata is replaced.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

from . import _cod_curation_runtime as runtime_helper
from ._transactions import publish_directory
from .release_curation import _relative_records, _stop

PROFILE_PATH = Path(__file__).with_name("data") / "cod_curation_v41_profile.json"
METHOD = "cod_curation_v4_1_20260716"


class CODCurationError(RuntimeError):
    """The recorded input, dependency, asset or result contract did not pass."""


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def _canonical_sha(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def _limit_memory():
    import resource
    resource.setrlimit(resource.RLIMIT_AS, (6 * 1024**3, 6 * 1024**3))


def _check_asset(path, spec):
    if not path.is_file() or path.stat().st_size != spec["bytes"] or _sha(path) != spec["sha256"]:
        raise CODCurationError("Recorded COD asset missing or changed: " + str(path))


def _stage_assets(workflow_root, destination, profile, skip_charges):
    recorded = {}
    for name, spec in profile["assets"].items():
        relative = Path(name)
        if relative.is_absolute() or ".." in relative.parts:
            raise CODCurationError("Unsafe profile asset path")
        if skip_charges and spec["required_for"] == "charging":
            continue
        source = workflow_root / relative
        _check_asset(source, spec)
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        with source.open("rb") as reader, target.open("xb") as writer:
            shutil.copyfileobj(reader, writer, length=1024 * 1024)
        _check_asset(target, spec)
        recorded[name] = dict(spec)
    return recorded


def _validate_record(record, source_id, source_sha, artifact_root, skip_charges):
    if (record.get("pipeline_id") != METHOD or record.get("source_id") != source_id
            or record.get("source", {}).get("sha256") != source_sha
            or record.get("status") not in {"complete", "review", "excluded", "failed"}):
        raise CODCurationError("Unexpected COD result identity, method or status")
    children = record.get("children", [])
    if record["status"] == "complete" and (skip_charges or not children or record.get("review_reasons")):
        raise CODCurationError("Complete COD result lacks accepted charged children")
    for child in children:
        state = child.get("status")
        if state not in {"complete", "review", "excluded", "failed"}:
            raise CODCurationError("Invalid COD child status")
        if state == "complete":
            if (skip_charges or child.get("review_reasons") or not child.get("deliverables")
                    or child.get("invariants", {}).get("all_passed") is not True):
                raise CODCurationError("Accepted COD child lacks complete invariants")
            if child.get("category") == "ION_FSR" and child.get("asr", {}).get("executed") is not False:
                raise CODCurationError("Counterion-containing COD child executed ASR")
        elif state == "review" and not child.get("review_reasons"):
            raise CODCurationError("COD review result lacks a reason")
        if record["status"] == "complete" and state not in {"complete", "excluded"}:
            raise CODCurationError("Complete parent contains an unaccepted child")
        for item in child.get("deliverables", []):
            path = Path(item["path"]).resolve(strict=True)
            if artifact_root.resolve() not in path.parents:
                raise CODCurationError("COD artifact escaped its private output")
            _check_asset(path, item["identity"])
            if state == "complete":
                charge = item.get("pacman", {})
                value = charge.get("charge_sum")
                if (item.get("charged") is not True or charge.get("charge_count") != item.get("atom_count")
                        or type(value) not in (int, float) or not math.isfinite(value)
                        or abs(value) > 1e-5 or charge.get("used_for_structural_selection") is not False):
                    raise CODCurationError("Accepted COD artifact lacks valid post-curation charges")


def curate_cod_cif_v41(cif_path, *, cod_id, output_dir, python, workflow_root,
                      skip_charges=False, timeout_seconds=600):
    """Apply the recorded COD policy to a raw source CIF, without overwriting it.

    Supply the transferred COD workflow root with its exact scripts, policy
    tables and models, and a separate Python 3.9 runtime with the pinned package
    versions. All required assets are verified and copied privately. No model
    download, fallback dependency, CR/NCR classification or release promotion
    occurs. Uncharged proposals remain REVIEW. Historical charge equality is
    not guaranteed across numerical environments, even at equal package versions.
    """
    if os.name != "posix":
        raise CODCurationError("The recorded COD runtime requires POSIX")
    if not isinstance(cod_id, str) or re.fullmatch(r"[0-9]{1,12}", cod_id) is None:
        raise ValueError("cod_id must be a numeric COD identifier string")
    if type(skip_charges) is not bool or type(timeout_seconds) is not int or timeout_seconds < 1:
        raise ValueError("Use a boolean skip_charges and positive integer timeout_seconds")
    destination = Path(output_dir).absolute()
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    if not destination.parent.is_dir():
        raise FileNotFoundError(destination.parent)
    source = Path(cif_path).resolve(strict=True)
    if not source.is_file() or source.stat().st_size > 64 * 1024**2:
        raise ValueError("Expected one regular source CIF no larger than 64 MiB")
    executable = Path(python).absolute()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise CODCurationError("External Python is missing or not executable")
    workflow_root = Path(workflow_root).resolve(strict=True)
    profile = json.loads(PROFILE_PATH.read_text())
    if profile.get("profile") != METHOD:
        raise CODCurationError("Unexpected packaged COD profile")
    source_sha = _sha(source)
    source_id = "cod_" + cod_id
    with tempfile.TemporaryDirectory(prefix=".cod-curation-v41-", dir=destination.parent) as temporary:
        root = Path(temporary)
        for name in ("workflow", "input", "output", "home", "tmp", "guard"):
            (root / name).mkdir(mode=0o700)
        output = root / "output"
        local_workflow = root / "workflow"
        assets = _stage_assets(workflow_root, local_workflow, profile, skip_charges)
        copied = root / "input" / (source_id + ".cif")
        shutil.copyfile(source, copied)
        if _sha(copied) != source_sha:
            raise CODCurationError("Source changed while copying")
        guard = root / "guard" / "sitecustomize.py"
        shutil.copyfile(runtime_helper.__file__, guard)
        local_profile = root / "profile.json"
        _json(local_profile, profile)
        environment = {"PATH": str(executable.parent) + ":/usr/bin:/bin",
                       "HOME": str(root / "home"), "TMPDIR": str(root / "tmp"),
                       "PYTHONPATH": str(root / "guard"), "PYTHONNOUSERSITE": "1",
                       "PYTHONDONTWRITEBYTECODE": "1", "PYTHONHASHSEED": "0",
                       "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                       "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1", "LC_ALL": "C.UTF-8",
                       "LANG": "C.UTF-8", "TZ": "UTC", "MPLCONFIGDIR": str(root / "tmp/matplotlib")}
        probe = subprocess.run([str(executable), "-B", str(guard), str(local_profile), "1" if skip_charges else "0"],
                               env=environment, capture_output=True, text=True, timeout=min(timeout_seconds, 45))
        if probe.returncode:
            raise CODCurationError("COD runtime preflight failed: " + probe.stderr[-4000:])
        runtime = json.loads(probe.stdout)
        _json(output / "runtime.json", runtime)
        entry = {"source_id": source_id, "cod_id": cod_id, "source_sha256": source_sha,
                 "source_cif_relative_path": copied.name, "metadata_status": "not_supplied"}
        manifest = output / "input_manifest.jsonl"
        manifest.write_text(json.dumps(entry, sort_keys=True, separators=(",", ":")) + "\n")
        fingerprint_payload = {"pipeline_id": METHOD, "record_schema_version": 1,
                               "scientific_policy": profile["scientific_policy"], "assets": assets,
                               "input_manifest_sha256": _sha(manifest), "runtime": runtime,
                               "reference_cpu_fingerprint_sha256": profile["original_cpu_fingerprint_sha256"]}
        fingerprint = output / "execution_fingerprint.json"
        fingerprint_sha = _canonical_sha(fingerprint_payload)
        _json(fingerprint, {"payload": fingerprint_payload, "sha256": fingerprint_sha})
        artifact_root = output / "artifacts"
        command = [str(executable), "-B", str(local_workflow / "scripts/curate_cod_one_v4_1.py"),
                   "--source", str(copied), "--source-id", source_id, "--cod-id", cod_id,
                   "--expected-source-sha256", source_sha, "--manifest", str(manifest),
                   "--manifest-sha256", _sha(manifest), "--pipeline-fingerprint", str(fingerprint),
                   "--pipeline-fingerprint-sha256", fingerprint_sha, "--output-root", str(artifact_root),
                   "--attempt-id", "recorded_method", "--charge-mode", "cpu-reference",
                   "--skip-validators", "--keep-work"]
        if skip_charges:
            command.append("--skip-charge")
        timed_out = False
        with (output / "stdout.txt").open("w") as stdout, (output / "stderr.txt").open("w") as stderr:
            process = subprocess.Popen(command, cwd=root, env=environment, stdout=stdout, stderr=stderr,
                                       start_new_session=True, preexec_fn=_limit_memory)
            try:
                process.wait(timeout=timeout_seconds)
            except subprocess.TimeoutExpired:
                timed_out = True
                _stop(process)
            except BaseException:
                _stop(process)
                raise
        if _sha(source) != source_sha or _sha(copied) != source_sha:
            raise CODCurationError("Original or private source was modified")
        for name, spec in assets.items():
            _check_asset(local_workflow / name, spec)
        record_path = artifact_root / f"state/attempts/{cod_id[:3]}/{source_id}/recorded_method.json"
        protocol = json.loads(record_path.read_text()) if record_path.is_file() else {}
        if not timed_out and process.returncode == 0:
            _validate_record(protocol, source_id, source_sha, artifact_root, skip_charges)
            status = {"complete": "COMPLETE", "review": "REVIEW", "excluded": "EXCLUDED", "failed": "ERROR"}[protocol["status"]]
        else:
            status = "TIMEOUT" if timed_out else "ERROR"
        normalized = _relative_records(protocol, output)
        if normalized.get("source"):
            normalized["source"]["path"] = str(source)
        _json(output / "protocol_record.json", normalized)
        result = {"schema_version": "coremof-cod-curation-v41-replay/1.0", "profile": METHOD,
                  "source_id": source_id, "cod_id": cod_id, "input_sha256": source_sha,
                  "execution_status": status, "children": normalized.get("children", []),
                  "review_reasons": normalized.get("review_reasons", []),
                  "errors": normalized.get("errors", []) or ([{"type": status, "message": "See stderr.txt and protocol_record.json"}]
                                                            if status in {"TIMEOUT", "ERROR"} else []),
                  "release_eligible": False, "exact_historical_equality_proven": False}
        receipt = {"profile": METHOD, "profile_sha256": _sha(PROFILE_PATH), "assets": assets,
                   "runtime_helper_sha256": _sha(guard), "input_sha256": source_sha,
                   "execution_fingerprint_sha256": fingerprint_sha, "returncode": process.returncode,
                   "timeout_seconds": timeout_seconds, "memory_limit_gib": 6, "skip_charges": skip_charges,
                   "original_cif_modified": False, "release_metadata_promoted": False,
                   "external_connections": "disabled", "models_downloaded": False,
                   "scope": "Recorded COD method, separate from CSD additions and historical base curation",
                   "numerical_limit": "Historical coordinates and PACMAN charges may differ across numerical runtimes; do not replace frozen CIFs."}
        _json(output / "record.json", result)
        _json(output / "receipt.json", receipt)
        publish_directory(output, destination, overwrite=False)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cif_path")
    for name in ("cod-id", "output-dir", "python", "workflow-root"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--skip-charges", action="store_true")
    parser.add_argument("--timeout-seconds", type=int, default=600)
    result = curate_cod_cif_v41(**vars(parser.parse_args(argv)))
    print(json.dumps({key: result[key] for key in ("source_id", "execution_status")}))
    return 0 if result["execution_status"] == "COMPLETE" else 1


if __name__ == "__main__":
    raise SystemExit(main())
