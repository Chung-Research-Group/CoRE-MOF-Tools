"""Full-precision RAC descriptors using the frozen release calculation.

The generic ``RACs`` API keeps its legacy rounding and defaults. This module
requires the sealed external runtime and molSimplify source archive, validates
their bytes, and performs no installation, CIF repair or release promotion.
"""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path, PurePosixPath
import signal
import stat
import subprocess
import sys
import tarfile
import tempfile
import time
import traceback

from ._transactions import publish_directory


PROFILE = "coremof-racs-ms173-pmg2023"
SOURCE_COMMIT = "60211676a71039f37adce57a60505f2bdaf7f184"
SOURCE_SHA256 = "e8fe1bda764ae107f976fc3b989f4aab29949343dbf89fe87f0cca2ecdc0f883"
SOURCE_SIZE = 1751040
PYTHON_SHA256 = "34f89896a80b3edb926a2402fa8c883d3dc938a7d790ca01005b751929cacaf9"
TREE_SHA256 = "3cff3cafeb7fa7865b0e9e421dfc1e343feb0aa614295ec340db86706831a6f2"
SCHEMA_SHA256 = {
    3: "1faf4ba9986b7c522929d21c166a04706c4a19d9e42340cc5b460e2c9d0e251a",
    5: "5e5576c2495a101dbc1eeda65c2989855f6099be52e4f3563ecfdd19df5a162e",
}
VERSIONS = {
    "monty": "2024.12.10", "networkx": "2.8.8", "numpy": "1.26.4",
    "pandas": "1.5.3", "pymatgen": "2023.9.10", "PyYAML": "6.0.2",
    "scikit-learn": "1.3.2", "scipy": "1.13.1", "spglib": "2.6.0",
}
CALL_ARGUMENTS = {
    "graph_provided": False, "wiggle_room": 1, "max_num_atoms": 6000,
    "get_sbu_linker_bond_info": False, "surrounded_sbu_file_generation": False,
    "detect_1D_rod_sbu": False,
}
RUNTIME_ENV = {
    "PYTHONNOUSERSITE": "1", "PYTHONDONTWRITEBYTECODE": "1",
    "PYTHONHASHSEED": "0", "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
    "LC_ALL": "C.UTF-8", "LANG": "C.UTF-8", "TZ": "UTC",
}


class ReleaseRACError(RuntimeError):
    """An input, runtime or output violates the frozen RAC contract."""


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for data in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(data)
    return digest.hexdigest()


def _check(path, digest, size=None):
    path = Path(path)
    if not path.is_file() or (size is not None and path.stat().st_size != size):
        raise ReleaseRACError(f"Missing file or size mismatch: {path}")
    if _sha(path) != digest:
        raise ReleaseRACError(f"SHA-256 mismatch: {path}")
    return path


def _schema(depth):
    if type(depth) is not int or depth not in SCHEMA_SHA256:
        raise ValueError("The frozen RAC profile supports only integer depth 3 or 5")
    templates = []
    for prefix in ("D_func", "func"):
        for prop in ("I", "S", "T", "Z", "alpha", "chi"):
            templates.append(("Function-group", prefix + "-" + prop + "-{depth}-all"))
    for prefix, props, suffix in (
        ("D_lc", ("I", "S", "T", "Z", "alpha", "chi"), "-all"),
        ("f-lig", ("I", "S", "T", "Z", "chi"), ""),
        ("lc", ("I", "S", "T", "Z", "alpha", "chi"), "-all"),
    ):
        for prop in props:
            templates.append(("Linker", prefix + "-" + prop + "-{depth}" + suffix))
    for prefix in ("D_mc", "f", "mc"):
        for prop in ("I", "S", "T", "Z", "chi"):
            templates.append(("Metal", prefix + "-" + prop + "-{depth}-all"))
    mapping = {template.format(depth=d): "RACs_" + group + "_" + template.format(depth=d)
               for group, template in templates for d in range(depth + 1)}
    digest = hashlib.sha256(json.dumps(sorted(mapping.values()), separators=(",", ":")).encode()).hexdigest()
    if len(templates) != 44 or digest != SCHEMA_SHA256[depth]:
        raise ReleaseRACError("Internal RAC schema differs from the frozen schema")
    return {key: mapping[key] for key in sorted(mapping)}


def _validated_values(names, values, depth):
    if hasattr(names, "tolist"):
        names = names.tolist()
    if hasattr(values, "tolist"):
        values = values.tolist()
    expected = _schema(depth)
    if not isinstance(names, (list, tuple)) or not isinstance(values, (list, tuple)):
        raise ReleaseRACError("molSimplify returned non-sequence descriptors")
    names = [str(name) for name in names]
    if len(names) != len(expected) or len(set(names)) != len(names) or set(names) != set(expected):
        raise ReleaseRACError(f"RAC schema is not the exact {len(expected)}-name set")
    if len(values) != len(names):
        raise ReleaseRACError("RAC name/value lengths differ")
    by_name = {}
    for name, value in zip(names, values):
        number = float(value)
        if not math.isfinite(number):
            raise ReleaseRACError(f"Non-finite RAC descriptor: {name}")
        by_name[name] = number
    # Preserve full float precision and both signs of zero. No rounding or
    # filling of missing descriptors is permitted on this route.
    return {column: by_name[raw] for raw, column in expected.items()}


def _relative_path(value):
    path = PurePosixPath(value)
    if path.is_absolute() or ".." in path.parts or str(path) != value:
        raise ReleaseRACError(f"Unsafe archive/manifest path: {value}")
    return path


def _verify_environment_tree(prefix, manifest):
    """Verify the complete sealed numerical environment, not just versions."""
    tree = json.loads(_check(manifest, TREE_SHA256).read_text())
    entries = tree["entries"]
    expected_paths = {item["path"] for item in entries}
    if len(expected_paths) != len(entries) or len(entries) != tree["entry_count"]:
        raise ReleaseRACError("Environment manifest has duplicate/missing paths")
    prefix = Path(prefix).resolve()
    for item in entries:
        path = prefix / _relative_path(item["path"])
        info = path.lstat()
        kind = item["type"]
        if format(stat.S_IMODE(info.st_mode), "04o") != item["mode"]:
            raise ReleaseRACError(f"Sealed environment mode changed: {path}")
        if kind == "symlink":
            if not stat.S_ISLNK(info.st_mode) or os.readlink(path) != item["target"]:
                raise ReleaseRACError(f"Environment symlink changed: {path}")
        elif kind == "directory":
            if not stat.S_ISDIR(info.st_mode):
                raise ReleaseRACError(f"Environment directory changed: {path}")
        elif kind == "file":
            if not stat.S_ISREG(info.st_mode) or info.st_nlink != item["link_count"]:
                raise ReleaseRACError(f"Environment file identity changed: {path}")
            _check(path, item["sha256"], item["size_bytes"])
        else:
            raise ReleaseRACError(f"Unknown environment entry type: {kind}")
    observed = {"."}
    for parent, directories, files in os.walk(prefix, followlinks=False):
        for name in directories + files:
            observed.add((Path(parent) / name).relative_to(prefix).as_posix())
    if observed != expected_paths:
        raise ReleaseRACError("Sealed environment contains added or missing entries")
    return {"manifest_sha256": TREE_SHA256, "entry_count": len(entries),
            "regular_file_bytes": tree["regular_file_bytes"],
            "tree_merkle_sha256": tree["tree_merkle_sha256"]}


def _extract_source(archive, destination):
    _check(archive, SOURCE_SHA256, SOURCE_SIZE)
    seen = set()
    total = 0
    with tarfile.open(archive, "r:") as bundle:
        for member in bundle.getmembers():
            # Tar directory names conventionally end with a slash.
            relative = _relative_path(member.name.rstrip("/"))
            if relative.as_posix() in seen or len(seen) >= 1000:
                raise ReleaseRACError("Duplicate or excessive source archive members")
            seen.add(relative.as_posix())
            target = Path(destination) / relative
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
            elif member.isfile():
                total += member.size
                if total > 16 * 1024 * 1024:
                    raise ReleaseRACError("Source archive exceeds its size bound")
                target.parent.mkdir(parents=True, exist_ok=True)
                with bundle.extractfile(member) as stream:
                    payload = stream.read()
                if len(payload) != member.size:
                    raise ReleaseRACError("Truncated source archive member")
                target.write_bytes(payload)
            else:
                raise ReleaseRACError("Non-file source archive member")
    return Path(destination)


def _json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def _validate_result(result, structure_id, cif_sha, depth):
    schema = _schema(depth)
    for name, expected in (("schema", "coremof-release-racs/1.0"), ("profile", PROFILE),
        ("structure_id", structure_id), ("cif_sha256", cif_sha), ("depth", depth),
        ("schema_sha256", SCHEMA_SHA256[depth]), ("metadata_names", list(schema.values())),
        ("call_arguments", CALL_ARGUMENTS)):
        if result.get(name) != expected:
            raise ReleaseRACError(f"Worker record identity/schema mismatch: {name}")
    if result.get("execution_status") == "SUCCESS":
        values, hex_values = result.get("descriptors"), result.get("values_float_hex")
        if result.get("available") is not True or not isinstance(values, dict) or set(values) != set(schema.values()):
            raise ReleaseRACError("Successful RAC output is incomplete")
        if not isinstance(hex_values, dict) or set(hex_values) != set(values):
            raise ReleaseRACError("RAC float.hex vector is incomplete")
        for name, value in values.items():
            if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value):
                raise ReleaseRACError("RAC vector contains a non-finite/nonnumeric value")
            if float(value).hex() != hex_values[name]:
                raise ReleaseRACError("RAC values and exact float.hex representation differ")
    elif result.get("execution_status") in {"ERROR", "TIMEOUT"}:
        if result.get("available") is not False or result.get("descriptors") is not None or result.get("values_float_hex") is not None:
            raise ReleaseRACError("Unavailable RAC result must have a null entire vector")
        if not isinstance(result.get("error"), dict) or not result["error"].get("type"):
            raise ReleaseRACError("Unavailable RAC result lacks a diagnostic")
    else:
        raise ReleaseRACError("Unknown RAC execution status")


def _worker(request_path):
    request = json.loads(Path(request_path).read_text())
    private = Path(request["private_root"])
    prefix = Path(sys.executable).resolve().parent.parent
    _check(sys.executable, PYTHON_SHA256)
    if tuple(sys.version_info[:3]) != (3, 9, 23) or hash("coremof-rac5-hash-seed-v1") != -4169351043222113396:
        raise ReleaseRACError("Python version or hash seed differs from the release runtime")
    tree = _verify_environment_tree(prefix, request["environment_manifest"])
    versions = {name: importlib.metadata.version(name) for name in VERSIONS}
    if versions != VERSIONS:
        raise ReleaseRACError(f"Pinned dependency versions differ: {versions}")
    source = _extract_source(request["source_archive"], private / "source")
    if "molSimplify" in sys.modules:
        raise ReleaseRACError("molSimplify was imported before source verification")
    sys.path.insert(0, str(source))
    from molSimplify.Informatics.MOF.MOF_descriptors import get_MOF_descriptors
    module = Path(sys.modules[get_MOF_descriptors.__module__].__file__).resolve()
    if source.resolve() not in module.parents:
        raise ReleaseRACError("molSimplify imported outside the frozen source archive")
    _json(private / "runtime.json", {"tree": tree, "versions": versions,
                                    "python_sha256": PYTHON_SHA256})
    cif = private / "input" / (request["structure_id"] + ".cif")
    _check(cif, request["cif_sha256"])
    result = {"schema": "coremof-release-racs/1.0", "profile": PROFILE,
        "structure_id": request["structure_id"], "cif_sha256": request["cif_sha256"],
        "depth": request["depth"], "schema_sha256": SCHEMA_SHA256[request["depth"]],
        "metadata_names": list(_schema(request["depth"]).values()),
        "call_arguments": dict(CALL_ARGUMENTS), "execution_status": "ERROR",
        "available": False, "descriptors": None, "values_float_hex": None}
    work = private / "calculation"
    work.mkdir()
    started = time.monotonic()
    with (private / "calculation.log").open("w") as log:
        try:
            with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                names, values = get_MOF_descriptors(data=str(cif), depth=request["depth"],
                    path=str(work), xyzpath=str(work / (request["structure_id"] + ".xyz")),
                    **CALL_ARGUMENTS)
            descriptors = _validated_values(names, values, request["depth"])
            result.update(execution_status="SUCCESS", available=True, descriptors=descriptors,
                values_float_hex={name: value.hex() for name, value in descriptors.items()})
        except (Exception, SystemExit) as error:
            result["error"] = {"type": type(error).__name__, "message": str(error),
                               "traceback": traceback.format_exc()}
    _check(cif, request["cif_sha256"])
    result["elapsed_calculation_seconds"] = time.monotonic() - started
    _json(private / "record.json", result)


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


def calculate_release_racs(cif_path, *, structure_id, output_dir, python,
                           source_archive, environment_manifest, depth=5,
                           timeout_seconds=1200):
    """Calculate exact, unrounded release RAC3/RAC5 with the sealed runtime.

    The timeout covers runtime verification and calculation. The full
    environment check is intentionally explicit and can be costly on shared
    storage. Run this reproduction route on a compute node, not a login host.
    Scientific errors have a null *entire* vector and a diagnostic; missing or
    changed runtime files raise ``ReleaseRACError``. Outputs are private and
    never replace existing directories or registered release evidence.
    """
    if os.name != "posix":
        raise ReleaseRACError("The sealed release runtime requires POSIX")
    schema = _schema(depth)
    from .identifiers import parse_core_id
    parse_core_id(structure_id)
    if type(timeout_seconds) is not int or not 1 <= timeout_seconds <= 3600:
        raise ValueError("timeout_seconds must be an integer in 1..3600")
    interpreter = _check(Path(python).resolve(strict=True), PYTHON_SHA256)
    _check(source_archive, SOURCE_SHA256, SOURCE_SIZE)
    _check(environment_manifest, TREE_SHA256)
    original = Path(cif_path).resolve(strict=True)
    cif_bytes = original.read_bytes()
    cif_sha = hashlib.sha256(cif_bytes).hexdigest()
    destination = Path(output_dir).absolute()
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    if not destination.parent.is_dir():
        raise FileNotFoundError(destination.parent)
    with tempfile.TemporaryDirectory(prefix=".coremof-racs-", dir=destination.parent) as temporary:
        private = Path(temporary)
        (private / "input").mkdir()
        (private / "input" / (structure_id + ".cif")).write_bytes(cif_bytes)
        environment = {**RUNTIME_ENV, "PATH": str(interpreter.parent) + ":/usr/bin:/bin"}
        for variable, child in (
            ("HOME", "home"), ("XDG_CONFIG_HOME", "xdg_config"),
            ("XDG_CACHE_HOME", "xdg_cache"), ("XDG_DATA_HOME", "xdg_data"),
            ("XDG_STATE_HOME", "xdg_state"), ("TMPDIR", "tmp"),
        ):
            path = private / child
            path.mkdir(mode=0o700)
            environment[variable] = str(path)
        request = {"structure_id": structure_id, "private_root": str(private),
            "source_archive": str(Path(source_archive).resolve()),
            "environment_manifest": str(Path(environment_manifest).resolve()),
            "cif_sha256": cif_sha, "depth": depth}
        _json(private / "request.json", request)
        prefix = interpreter.parent.parent
        paths = [str(Path(__file__).resolve().parent.parent),
                 str(prefix / "lib/python3.9"), str(prefix / "lib/python3.9/lib-dynload"),
                 str(prefix / "lib/python3.9/site-packages")]
        launcher = ("import sys; sys.path[:] = " + repr(paths) + "; "
                    "from CoREMOF.release_racs import _worker; _worker(sys.argv[1])")
        command = [str(interpreter), "-S", "-s", "-B", "-u", "-c", launcher,
                   str(private / "request.json")]
        timed_out = False
        with (private / "worker.log").open("wb") as log:
            process = subprocess.Popen(command, cwd=private, env=environment,
                stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            try:
                process.wait(timeout=timeout_seconds)
            except subprocess.TimeoutExpired:
                timed_out = True
                _stop(process)
            except BaseException:
                _stop(process)
                raise
        _check(original, cif_sha)
        if timed_out:
            result = {"schema": "coremof-release-racs/1.0", "profile": PROFILE,
                "structure_id": structure_id, "cif_sha256": cif_sha, "depth": depth,
                "schema_sha256": SCHEMA_SHA256[depth], "metadata_names": list(schema.values()),
                "call_arguments": dict(CALL_ARGUMENTS), "execution_status": "TIMEOUT",
                "available": False, "descriptors": None, "values_float_hex": None,
                "error": {"type": "TimeoutError", "message": "Runtime verification or calculation timed out"}}
        elif process.returncode:
            tail = (private / "worker.log").read_text(errors="replace")[-6000:]
            raise ReleaseRACError(f"Sealed RAC worker failed ({process.returncode}): {tail}")
        else:
            result = json.loads((private / "record.json").read_text())
        _validate_result(result, structure_id, cif_sha, depth)
        payload = private / "publish"
        payload.mkdir(mode=0o700)
        _json(payload / "record.json", result)
        for name in ("worker.log", "calculation.log", "runtime.json"):
            if (private / name).is_file():
                (payload / name).write_bytes((private / name).read_bytes())
        _json(payload / "receipt.json", {
            "profile": PROFILE, "cif_sha256": cif_sha, "source_archive_sha256": SOURCE_SHA256,
            "molsimplify_commit": SOURCE_COMMIT, "environment_tree_sha256": TREE_SHA256,
            "python_sha256": PYTHON_SHA256, "adapter_sha256": _sha(__file__),
            "schema_sha256": SCHEMA_SHA256[depth], "depth": depth,
            "record_sha256": _sha(payload / "record.json"),
            "runtime_verified": (private / "runtime.json").is_file(),
            "timeout_seconds": timeout_seconds, "release_promotion": False,
        })
        publish_directory(payload, destination, overwrite=False)
        return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("cif-path", "structure-id", "output-dir", "python", "source-archive", "environment-manifest"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--depth", type=int, choices=(3, 5), default=5)
    parser.add_argument("--timeout-seconds", type=int, default=1200)
    record = calculate_release_racs(**vars(parser.parse_args()))
    print(json.dumps({key: record[key] for key in ("structure_id", "depth", "execution_status", "available")}))
