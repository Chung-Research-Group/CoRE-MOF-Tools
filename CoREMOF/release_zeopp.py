"""Isolated replay of the recorded N2/He and framework-dimension Zeo++ methods."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import tempfile

from . import _release_zeopp_n2_he_protocol as probe
from . import _release_zeopp_framework_protocol as framework
from ._transactions import publish_directory


PROFILE = "coremof-release-zeopp-n2-he-framework-0.4.7"
NETWORK_SHA256 = "6f55b8c36e2e03a027b30c2350f190c614d1adc83211e429266478f99128bc06"
PROBE_SHA256 = "396e5a756815271d92ae40e9bf2320d723aa715741724a167ff7096abdfd8630"
FRAMEWORK_SHA256 = "9683b7202c030994b8c9e9e607e6c0c155e4dea0aca4b65486edeafb08b50668"
POLICY_SHA256 = "65583539c80f11ac19d32ed601151bf3a01fb5a9fd85f8347fde08b437296bbc"


class ReleaseZeoppError(RuntimeError):
    """The input, external binary, method or output contract was not met."""


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _check(path, digest):
    path = Path(path)
    if not path.is_file() or _sha(path) != digest:
        raise ReleaseZeoppError(f"Missing or SHA-256-mismatched input: {path}")
    return path


def _profile(network):
    _check(network, NETWORK_SHA256)
    _check(probe.__file__, PROBE_SHA256)
    _check(framework.__file__, FRAMEWORK_SHA256)
    return _check(Path(__file__).with_name("data") / "zeopp_probe_namespace_v2.json", POLICY_SHA256)


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


def _validate(record, structure_id, cif_sha, kind):
    expected_protocol = probe.PROTOCOL_ID if kind == "n2_he" else framework.PROTOCOL_ID
    if (record.get("structure_id") != structure_id or record.get("protocol_id") != expected_protocol
            or record.get("input", {}).get("cif_sha256") != cif_sha
            or record.get("execution_status") not in {"SUCCESS", "ERROR"}):
        raise ReleaseZeoppError("Unexpected result identity, protocol or status")
    values = record.get("features" if kind == "n2_he" else "framework_dimension")
    if record["execution_status"] == "ERROR":
        if values is not None:
            raise ReleaseZeoppError("Failed calculation contains a scientific value")
        return
    if kind == "n2_he":
        schemas = {"intrinsic_props": probe.INTRINSIC_FEATURES,
                   "N2_probe_props": probe.N2_FEATURES, "He_probe_props": probe.HE_FEATURES}
        if not isinstance(values, dict) or set(values) != set(schemas):
            raise ReleaseZeoppError("Incomplete Zeo++ feature namespaces")
        for name, keys in schemas.items():
            if set(values[name]) != set(keys):
                raise ReleaseZeoppError("Incomplete Zeo++ feature schema")
            for key, value in values[name].items():
                if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                    raise ReleaseZeoppError("Invalid Zeo++ feature")
                if key.endswith("VF") and value > 1:
                    raise ReleaseZeoppError("Invalid void fraction")
        channel = record.get("N2_channel_topology") or {}
        dims = channel.get("channel_dimensions")
        if (not isinstance(dims, list) or any(type(x) is not int or x not in (1, 2, 3) for x in dims)
                or channel.get("channel_count") != len(dims)
                or channel.get("maximum_channel_dimension") != max(dims, default=0)
                or values["N2_probe_props"]["channel_dimension"] != max(dims, default=0)):
            raise ReleaseZeoppError("Inconsistent channel topology")
    else:
        if not isinstance(values, dict) or "raw_line" not in values:
            raise ReleaseZeoppError("Missing framework topology")
        if framework.parse_strinfo(values["raw_line"]) != values:
            raise ReleaseZeoppError("Framework values differ from raw output")


def calculate_release_zeopp(cif_path, *, structure_id, output_dir, network, timeout_seconds=300):
    """Calculate N2/He features and bonded-framework dimensionality separately.

    Use the recorded binary, not an arbitrary executable on PATH. The output
    directory must not exist. Input CIFs are copied, not repaired. The timeout
    bounds each external command, with an additional outer process-group bound.
    Explicit component errors remain unavailable and never trigger a fallback.
    """
    if os.name != "posix":
        raise ReleaseZeoppError("The recorded external binary requires POSIX")
    if not isinstance(structure_id, str) or re.fullmatch(
        r"(?:ASR|FSR|ION)-(?:COD|CSD|SI)-(?:[0-9]{4}|UNKN)-[0-9]{4,}", structure_id
    ) is None:
        raise ValueError("Use a public CoRE-MOF structure ID")
    if type(timeout_seconds) is not int or timeout_seconds <= 0:
        raise ValueError("timeout_seconds must be a positive integer")
    destination = Path(output_dir).absolute()
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    if not destination.parent.is_dir():
        raise FileNotFoundError(destination.parent)
    network = Path(network).resolve(strict=True)
    policy = _profile(network)
    source = Path(cif_path).resolve(strict=True)
    payload = source.read_bytes()
    cif_sha = hashlib.sha256(payload).hexdigest()
    with tempfile.TemporaryDirectory(prefix=".zeopp-replay-", dir=destination.parent) as temporary:
        private = Path(temporary)
        for name in ("input", "home", "tmp", "work", "output"):
            (private / name).mkdir(mode=0o700)
        output = private / "output"
        cif = private / "input" / (structure_id + ".cif")
        cif.write_bytes(payload)
        row = dict.fromkeys(framework.MANIFEST_FIELDS, "")
        row.update(manifest_schema_version="1.0", row_index="0", structure_id=structure_id,
                   canonical_cif_version="user-supplied-exact-bytes", canonical_cif_path=str(cif),
                   cif_basename=cif.name, cif_size_bytes=str(len(payload)), cif_sha256=cif_sha,
                   source_family=structure_id.split("-")[1])
        manifest = private / "manifest.csv"
        with manifest.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=framework.MANIFEST_FIELDS)
            writer.writeheader()
            writer.writerow(row)
        env = {"PATH": "/usr/bin:/bin", "HOME": str(private / "home"),
               "TMPDIR": str(private / "tmp"), "LANG": "C.UTF-8", "LC_ALL": "C.UTF-8",
               "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "PYTHONDONTWRITEBYTECODE": "1"}
        records, attempts = {}, []
        for kind, module, digest, calls in (("n2_he", probe, PROBE_SHA256, 5),
                                          ("framework", framework, FRAMEWORK_SHA256, 1)):
            work = private / "work" / kind
            work.mkdir()
            run_output = output / kind
            run_output.mkdir()
            command = [sys.executable, "-B", "-S", str(_check(module.__file__, digest)),
                       "--manifest", str(manifest), "--manifest-sha256", _sha(manifest),
                       "--row-index", "0", "--output-root", str(run_output),
                       "--private-root", str(work), "--network", str(network),
                       "--network-sha256", NETWORK_SHA256, "--runner-sha256", digest,
                       "--timeout-seconds", str(timeout_seconds)]
            if kind == "n2_he":
                command.extend(["--namespace-policy", str(policy), "--namespace-policy-sha256", POLICY_SHA256])
            with (run_output / "stdout.txt").open("w") as stdout, (run_output / "stderr.txt").open("w") as stderr:
                process = subprocess.Popen(command, cwd=work, env=env, stdout=stdout,
                                           stderr=stderr, start_new_session=True)
                timed_out = False
                try:
                    process.wait(timeout=calls * timeout_seconds + 30)
                except subprocess.TimeoutExpired:
                    timed_out = True
                    _stop(process)
                except BaseException:
                    _stop(process)
                    raise
            path = run_output / "records" / (structure_id + ".json")
            if process.returncode == 0 and path.is_file() and not timed_out:
                record = json.loads(path.read_text())
            else:
                record = {"structure_id": structure_id, "input": {"cif_sha256": cif_sha},
                          "protocol_id": module.PROTOCOL_ID, "execution_status": "ERROR",
                          "error_type": "TIMEOUT" if timed_out else "PROCESS_ERROR",
                          "error_message": f"Worker return code {process.returncode}; see stderr.txt",
                          "features": None, "N2_channel_topology": None, "framework_dimension": None}
                path.parent.mkdir(exist_ok=True)
                path.write_text(json.dumps(record, indent=2) + "\n")
            _validate(record, structure_id, cif_sha, kind)
            records[kind] = record
            attempts.append({"kind": kind, "returncode": process.returncode,
                             "execution_status": record["execution_status"], "record_sha256": _sha(path)})
        _check(source, cif_sha)
        _check(cif, cif_sha)
        _check(network, NETWORK_SHA256)
        successful = sum(x["execution_status"] == "SUCCESS" for x in records.values())
        features = records["n2_he"].get("features")
        periodicity = records["framework"].get("framework_dimension")
        result = {"schema_version": "coremof-release-zeopp-replay/1.0", "profile": PROFILE,
                  "structure_id": structure_id, "cif_sha256": cif_sha,
                  "execution_status": "SUCCESS" if successful == 2 else "PARTIAL" if successful else "ERROR",
                  "features": features, "N2_channel_topology": records["n2_he"].get("N2_channel_topology"),
                  "framework_dimension": ({k: v for k, v in periodicity.items() if k != "raw_line"}
                                          if periodicity is not None else None),
                  "components": {kind: {"execution_status": rec["execution_status"],
                                        "error_type": rec.get("error_type"), "error_message": rec.get("error_message")}
                                 for kind, rec in records.items()}}
        receipt = {"profile": PROFILE, "network_sha256": NETWORK_SHA256,
                   "probe_protocol_sha256": PROBE_SHA256, "framework_protocol_sha256": FRAMEWORK_SHA256,
                   "namespace_policy_sha256": POLICY_SHA256, "cif_sha256": cif_sha,
                   "timeout_seconds_per_command": timeout_seconds, "attempts": attempts,
                   "python_version": sys.version, "source_cif_modified": False, "release_metadata_promoted": False,
                   "historical_full_runtime_byte_identity_proven": False,
                   "monte_carlo_seed_control": "not exposed by the recorded CLI protocol"}
        for name, value in (("record.json", result), ("receipt.json", receipt)):
            (output / name).write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n")
        publish_directory(output, destination, overwrite=False)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cif_path", type=Path)
    for name in ("structure-id", "output-dir", "network"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--timeout-seconds", type=int, default=300)
    result = calculate_release_zeopp(**vars(parser.parse_args(argv)))
    print(json.dumps({key: result[key] for key in ("structure_id", "execution_status")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
