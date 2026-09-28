"""Fresh-process entry point for the recorded v4.1 curation protocol."""
from __future__ import annotations

import hashlib
import importlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import socket
import sys


PROTOCOL_SHA256 = "7e9b35f74530b75ba5e2d81313d82ba52fc6498733cdb17600ed715968ef797a"
PYTHON_SHA256 = "9024c0445314fe47972ac371dd460c4949aa9bb211ba8934f682c608bd209176"
SOURCE_HASHES = {
    "CoREMOF.curate": "72baca42a5238dc64268928a0b77f9d9b250ad0b9cd53a8cb4cb9c9517cd7628",
    "CoREMOF.utils.atoms_definitions": "3f987856e8a2204a49f41018aaef4676e3dd6015cb9923c7f7772b4a1dec7ea6",
    "CoREMOF.utils.ions_list": "11cd071d713bb2bc8a4fa4b20c972fb66b569b6b4ad2bde7ca161933ad1d6653",
    "ase.neighborlist": "9517884e9d488b820af452988e2a7c93d3c3ec31feb71d275c3a8dfbc35d1f8f",
    "PACMANCharge.pmcharge": "403c0dc5bec584bfb745f3aeae1bc944c29147fe052b79ece9888d287c02c485",
}
VERSIONS = {"ase": "3.23.0", "gemmi": "0.7.0", "pymatgen": "2024.8.9",
            "PACMAN-charge": "1.4.2", "torch": "2.7.0+cu118"}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def check(path, expected):
    path = Path(path)
    if not path.is_file() or sha(path) != expected:
        raise RuntimeError(f"Recorded curation source/runtime differs: {path}")
    return path


def deny_network(*args, **kwargs):
    raise RuntimeError("Curation replay does not download dependencies or models")


def runtime_evidence():
    check(sys.executable, PYTHON_SHA256)
    versions = {name: importlib.metadata.version(name) for name in VERSIONS}
    if versions != VERSIONS:
        raise RuntimeError(f"Curation runtime versions differ: {versions}")
    files = {}
    # Resolve and check code before importing the scientific modules. Their
    # package initializers remain covered by the no-network execution guard.
    for name, expected in SOURCE_HASHES.items():
        spec = importlib.util.find_spec(name)
        if spec is None or spec.origin is None:
            raise RuntimeError(f"Missing curation dependency: {name}")
        path = check(spec.origin, expected)
        files[name] = {"path": str(path), "sha256": expected}
    pacman = Path(files["PACMANCharge.pmcharge"]["path"]).parent
    # PACMAN 1.4.2 attempts downloads for any missing entry at import time.
    # Require the installed files and forbid network access, including when
    # charging is skipped. Only DDEC6 weights are used by the v4.1 protocol.
    for stem in ("cm5", "bader", "ddec", "repeat", "pbe", "bandgap"):
        for suffix in (".pth", ".pkl"):
            if not (pacman / (stem + suffix)).is_file():
                raise RuntimeError(f"Install PACMAN assets explicitly first: {stem}{suffix}")
    assets = {}
    for path in sorted(pacman.rglob("*")):
        if path.is_file() and (path.suffix in (".py", ".json") or path.name in ("ddec.pth", "ddec.pkl")):
            assets[path.relative_to(pacman).as_posix()] = {"bytes": path.stat().st_size, "sha256": sha(path)}
    return {"python_sha256": PYTHON_SHA256, "python_version": sys.version,
            "recorded_package_versions": versions, "source_files": files,
            "pacman_current_assets": assets,
            "additional_current_versions": {name: importlib.metadata.version(name)
                                            for name in ("numpy", "scipy", "pandas", "spglib")},
            "historical_full_runtime_byte_identity_proven": False}


def main():
    request = json.loads(Path(sys.argv[1]).read_text())
    socket.socket.connect = deny_network
    socket.create_connection = deny_network
    evidence = runtime_evidence()
    protocol_path = check(request["protocol"], PROTOCOL_SHA256)
    spec = importlib.util.spec_from_file_location("frozen_curation_v41", protocol_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    source = check(request["input"], request["input_sha256"])
    record = module.run_one(source, Path(request["output"]), {},
                            skip_pacman=request["skip_charges"])
    # A parser/standardizer works on the private input only. The caller checks
    # the original source separately and records any copied-input change.
    evidence["private_input_sha256_after"] = sha(source)
    for item in evidence["source_files"].values():
        check(item["path"], item["sha256"])
    pacman = Path(evidence["source_files"]["PACMANCharge.pmcharge"]["path"]).parent
    for name, item in evidence["pacman_current_assets"].items():
        check(pacman / name, item["sha256"])
    Path(request["runtime_receipt"]).write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
    Path(request["record"]).write_text(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
