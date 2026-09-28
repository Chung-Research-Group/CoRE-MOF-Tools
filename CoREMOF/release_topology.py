"""Explicit, isolated replay of the release CrystalNets 1.2.0 calculation.

The legacy topology API is unchanged. This module uses the byte-identical
release worker and result projection, both topology modes, and an explicitly
hash-bound external Julia installation. It never installs software or promotes
its output into a database release.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import stat
import subprocess
import tempfile
import time

from . import _release_topology_projection as projection
from ._transactions import publish_directory


PROFILE = "coremof-release-crystalnets-1.2.0"
VERSIONS = {"julia_version": "1.12.6", "crystalnets_version": "1.2.0"}
METHOD = {
    "structure_type": "MOF", "clusterings": ["SingleNodes", "AllNodes"],
    "exports_enabled": False, "warnings_enabled": False,
    "interpenetration_definition": "number of subnets in InterpenetratedTopologyResult",
}
PROFILE_FILES = {
    "run_crystalnets_one.jl": "182a213606ef2c997856efcadba596ebd6a2660e793eedb828800f67befc8d13",
    "run_crystalnets_shard.jl": "a994be05514a9b2db4ce1954a9b43f26e27c62bea2c673c526ad6d3b5ed09d92",
}
PROJECT_FILES = {
    "Project.toml": "bc413dc8352fd1abe87d3ca85db281b86eb3d1b5ca8b2f87ce0199dd1b25eecf",
    "Manifest.toml": "ee0ced50c1d7d62ecdc0785afbb1160265f95f583864e3d22e43e5dd055809a7",
}
PROJECTION_SHA256 = "3854774aec5cef76d322f4dd9270edd6d9a372f61d0ffdd957262574f6d66527"
RUNTIME_SCHEMA = "coremof-crystalnets-runtime/1.0"


class ReleaseTopologyError(RuntimeError):
    """The pinned method, runtime, inputs or result could not be verified."""


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _check(path, digest):
    path = Path(path)
    if not path.is_file() or _sha(path) != digest:
        raise ReleaseTopologyError(f"Missing or SHA-256-mismatched file: {path}")
    return path


def _profile(project):
    assets = Path(__file__).with_name("_topology_profile")
    _check(projection.__file__, PROJECTION_SHA256)
    for name, digest in PROFILE_FILES.items():
        _check(assets / name, digest)
    for name, digest in PROJECT_FILES.items():
        _check(Path(project) / name, digest)
    return assets


def _roots(julia, depot):
    julia = Path(julia).absolute()
    if julia.name != "julia" or julia.parent.name != "bin" or julia.is_symlink():
        raise ValueError("julia must name the installation's regular bin/julia executable")
    depot = Path(depot).resolve(strict=True)
    return {"julia": julia.parent.parent.resolve(strict=True),
            **{f"depot_{name}": depot / name for name in ("packages", "artifacts", "compiled")}}


def _tree(root):
    """Inventory loaded runtime trees without following symlink directories."""
    root = Path(root)
    if not root.is_dir() or root.is_symlink():
        raise ReleaseTopologyError(f"Missing or symlinked runtime root: {root}")
    entries = []
    def unreadable(error):
        raise ReleaseTopologyError(f"Unreadable runtime tree: {error}") from error

    for folder, dirs, files in os.walk(root, followlinks=False, onerror=unreadable):
        for name in sorted(dirs + files):
            path = Path(folder) / name
            info = path.lstat()
            entry = {"path": path.relative_to(root).as_posix()}
            if stat.S_ISLNK(info.st_mode):
                target = os.readlink(path)
                # Runtime symlinks must stay within the selected tree.
                try:
                    path.resolve(strict=True).relative_to(root.resolve())
                except (ValueError, OSError, RuntimeError) as exc:
                    raise ReleaseTopologyError(f"Escaping/broken runtime link: {path}") from exc
                entry.update(kind="symlink", target=target)
            elif stat.S_ISDIR(info.st_mode):
                entry.update(kind="directory")
            elif stat.S_ISREG(info.st_mode):
                entry.update(kind="file", size=info.st_size, sha256=_sha(path))
                after = path.stat()
                if (info.st_size, info.st_mtime_ns, info.st_ino) != (after.st_size, after.st_mtime_ns, after.st_ino):
                    raise ReleaseTopologyError(f"Runtime changed while hashing: {path}")
            else:
                raise ReleaseTopologyError(f"Non-regular runtime entry: {path}")
            entries.append(entry)
    return sorted(entries, key=lambda item: item["path"])


def capture_runtime_manifest(*, julia, depot, output_path):
    """Record the current external runtime, without asserting historical identity.

    Run on a compute node. This reads the full Julia installation plus depot
    packages, artifacts and compiled caches, about 1.4 GB for the tested profile.
    It does not execute Julia, install dependencies or authorize a release.
    """
    roots = _roots(julia, depot)
    result = {"schema_version": RUNTIME_SCHEMA, "profile": PROFILE,
              "scope": "current runtime snapshot, not retrospective proof of historical bytes",
              "trees": {name: _tree(root) for name, root in roots.items()}}
    destination = Path(output_path)
    with destination.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, sort_keys=True, indent=2)
        stream.write("\n")
    return _sha(destination)


def _verify_runtime(julia, depot, manifest, digest):
    if not isinstance(digest, str) or re.fullmatch(r"[0-9a-f]{64}", digest) is None:
        raise ValueError("Supply the recorded runtime_manifest_sha256")
    value = json.loads(_check(manifest, digest).read_text())
    roots = _roots(julia, depot)
    if (value.get("schema_version") != RUNTIME_SCHEMA or value.get("profile") != PROFILE
            or set(value.get("trees", {})) != set(roots)):
        raise ReleaseTopologyError("Unexpected runtime manifest schema/profile/trees")
    for name, root in roots.items():
        # Compare full membership, not just listed files. Additional source or
        # compiled-cache files must not become invisible runtime replacements.
        if value["trees"][name] != _tree(root):
            raise ReleaseTopologyError(f"Runtime tree differs from snapshot: {name}")
    return value


def _validate(raw, structure_id, cif_sha):
    if (raw.get("structure_id") != structure_id or raw.get("cif_sha256") != cif_sha
            or raw.get("schema_version") != "crystalnets-topology-result/1.0"
            or raw.get("software") != VERSIONS or raw.get("method") != METHOD):
        raise ReleaseTopologyError("Worker identity, software or method differs")
    if raw.get("execution_status") not in {"SUCCESS", "ERROR"}:
        raise ReleaseTopologyError("Unknown raw execution status")
    runtime = raw.get("runtime_seconds")
    if isinstance(runtime, bool) or not isinstance(runtime, (float, int)) or not math.isfinite(runtime) or runtime < 0:
        raise ReleaseTopologyError("Invalid worker duration")
    if not isinstance(raw.get("subnets"), list):
        raise ReleaseTopologyError("Missing subnet list")
    for index, subnet in enumerate(raw["subnets"], 1):
        if type(subnet.get("subnet_index")) is not int or subnet["subnet_index"] != index:
            raise ReleaseTopologyError("Invalid subnet order")
        for mode in ("single_node", "all_node"):
            part = subnet.get(mode)
            if not isinstance(part, dict) or part.get("status") not in {"SUCCESS", "NOT_AVAILABLE"}:
                raise ReleaseTopologyError("Invalid subnet mode/status")
            if part["status"] == "SUCCESS":
                if (not isinstance(part.get("topology"), str) or not part["topology"]
                        or type(part.get("dimension")) is not int or part["dimension"] not in range(4)):
                    raise ReleaseTopologyError("Invalid topology label or dimension")
    if raw["execution_status"] == "SUCCESS":
        for key in ("interpenetrated_subnet_count", "catenation_degree"):
            if type(raw.get(key)) is not int or raw[key] != len(raw["subnets"]):
                raise ReleaseTopologyError("Inconsistent subnet count")
    elif raw["subnets"] or any(raw.get(key) is not None for key in (
        "network_dimension", "interpenetrated_subnet_count", "catenation_degree",
        "single_node_net", "all_node_net", "all_subnets_single_all_agree",
    )):
        raise ReleaseTopologyError("An unavailable result contains scientific values")
    try:
        return projection.normalized_record(raw)[1]
    except (projection.AuditError, KeyError, TypeError) as exc:
        raise ReleaseTopologyError(str(exc)) from exc


def _unavailable(structure_id, cif_sha, kind, message, runtime):
    return {"schema_version": "crystalnets-topology-result/1.0",
            "structure_id": structure_id, "cif_sha256": cif_sha,
            "execution_status": "ERROR", "software": VERSIONS, "method": METHOD,
            "runtime_seconds": runtime, "error_type": kind, "error_message": message,
            "network_dimension": None, "interpenetrated_subnet_count": None,
            "catenation_degree": None, "single_node_net": None, "all_node_net": None,
            "all_subnets_single_all_agree": None, "subnets": []}


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


def calculate_release_topology(
    cif_path, *, structure_id, output_dir, julia, project, depot,
    runtime_manifest, runtime_manifest_sha256, timeout_seconds=300,
):
    """Replay both release topology modes on an isolated copy of one CIF.

    ``output_dir`` must not exist. The timeout bounds Julia, including startup.
    The full external runtime is verified before launch. Preserve scientific
    component errors as PARTIAL results and process failures as unavailable.
    No source CIF, installed dependency, existing result or release is edited.
    """
    if os.name != "posix":
        raise ReleaseTopologyError("The recorded Julia runtime requires POSIX")
    if not isinstance(structure_id, str) or re.fullmatch(
        r"(?:ASR|FSR|ION)-(?:COD|CSD|SI)-(?:[0-9]{4}|UNKN)-[0-9]{4,}", structure_id
    ) is None:
        raise ValueError("Use a public CoRE-MOF structure ID")
    if isinstance(timeout_seconds, bool) or not math.isfinite(timeout_seconds) or timeout_seconds <= 0:
        raise ValueError("timeout_seconds must be finite and positive")
    destination = Path(output_dir).absolute()
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    if not destination.parent.is_dir():
        raise FileNotFoundError(destination.parent)
    assets = _profile(project)
    julia, project, depot = (Path(path).resolve(strict=True) for path in (julia, project, depot))
    _verify_runtime(julia, depot, runtime_manifest, runtime_manifest_sha256)
    source = Path(cif_path).resolve(strict=True)
    payload = source.read_bytes()
    cif_sha = hashlib.sha256(payload).hexdigest()
    with tempfile.TemporaryDirectory(prefix=".topology-replay-", dir=destination.parent) as temporary:
        private = Path(temporary)
        for name in ("input", "work", "home", "tmp", "profile", "project", "output"):
            (private / name).mkdir(mode=0o700)
        isolated_cif = private / "input" / (structure_id + ".cif")
        isolated_cif.write_bytes(payload)
        for name, digest in PROFILE_FILES.items():
            (private / "profile" / name).write_bytes(_check(assets / name, digest).read_bytes())
        for name, digest in PROJECT_FILES.items():
            (private / "project" / name).write_bytes(_check(project / name, digest).read_bytes())
        output = private / "output"
        raw_path = output / "raw.json"
        env = {"PATH": "/usr/bin:/bin", "HOME": str(private / "home"),
               "TMPDIR": str(private / "tmp"), "LANG": "C.UTF-8", "LC_ALL": "C.UTF-8",
               "JULIA_DEPOT_PATH": str(depot), "JULIA_LOAD_PATH": "@:@stdlib",
               "JULIA_PKG_OFFLINE": "true", "JULIA_PKG_PRECOMPILE_AUTO": "0",
               "JULIA_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"}
        command = [str(julia), "--startup-file=no", "--history-file=no",
                   "--compiled-modules=existing", "--pkgimages=existing",
                   f"--project={private / 'project'}", "-e",
                   'using CrystalNets; @assert VERSION == v"1.12.6"; '
                   '@assert Base.pkgversion(CrystalNets) == v"1.2.0"; include(popfirst!(ARGS))',
                   str(private / "profile/run_crystalnets_one.jl"),
                   structure_id, str(isolated_cif), cif_sha, str(raw_path)]
        started = time.monotonic()
        timed_out = False
        with (output / "stdout.txt").open("w") as stdout, (output / "stderr.txt").open("w") as stderr:
            process = subprocess.Popen(command, cwd=private / "work", env=env,
                                       stdout=stdout, stderr=stderr, start_new_session=True)
            try:
                process.wait(timeout=timeout_seconds)
            except subprocess.TimeoutExpired:
                timed_out = True
                _stop(process)
            except BaseException:
                _stop(process)
                raise
        elapsed = time.monotonic() - started
        _check(source, cif_sha)
        _check(isolated_cif, cif_sha)
        if timed_out or process.returncode != 0:
            raw = _unavailable(structure_id, cif_sha, "TIMEOUT" if timed_out else "PROCESS_ERROR",
                               f"Julia return code {process.returncode}; see stderr.txt", elapsed)
        else:
            if not raw_path.is_file():
                raise ReleaseTopologyError("Julia returned no result")
            raw = json.loads(raw_path.read_text())
        rich = _validate(raw, structure_id, cif_sha)
        receipt = {"schema_version": "coremof-release-topology-replay/1.0", "profile": PROFILE,
                   "structure_id": structure_id, "cif_sha256": cif_sha,
                   "runtime_manifest_sha256": runtime_manifest_sha256,
                   "runtime_scope": "hash-bound current snapshot; historical bytes not independently attested",
                   "worker_sha256": PROFILE_FILES, "project_sha256": PROJECT_FILES,
                   "projection_sha256": PROJECTION_SHA256, "software": VERSIONS,
                   "execution_status": rich["execution_status"], "returncode": process.returncode,
                   "timeout_seconds": timeout_seconds, "elapsed_seconds": elapsed,
                   "source_cif_modified": False, "release_metadata_promoted": False}
        for name, value in (("raw.json", raw), ("record.json", rich), ("receipt.json", receipt)):
            (output / name).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
        publish_directory(output, destination, overwrite=False)
    return rich


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cif", type=Path)
    for name in ("structure-id", "output-dir", "julia", "project", "depot", "runtime-manifest", "runtime-manifest-sha256"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--timeout-seconds", type=float, default=300)
    args = vars(parser.parse_args(argv))
    args["cif_path"] = args.pop("cif")
    result = calculate_release_topology(**args)
    print(json.dumps({key: result[key] for key in ("structure_id", "execution_status", "topology_available")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
