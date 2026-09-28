"""Private runtime probe and offline guard for the recorded COD method."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
from pathlib import Path
import socket
import sys


def deny_network(*args, **kwargs):
    raise RuntimeError("External connections are disabled for COD curation replay")


def disable_network():
    socket.socket.connect = deny_network
    socket.socket.connect_ex = deny_network
    socket.create_connection = deny_network
    socket.getaddrinfo = deny_network


def inspect_runtime(profile, skip_charges=False):
    if sys.version_info[:2] != (3, 9):
        raise RuntimeError("The recorded COD method requires Python 3.9")
    required = dict(profile["packages"])
    if skip_charges:
        for name in ("torch", "PyCifRW"):
            required.pop(name)
    versions = {name: importlib.metadata.version(name) for name in required}
    differences = {name: {"expected": expected, "observed": versions[name]}
                   for name, expected in required.items() if versions[name] != expected}
    if differences:
        raise RuntimeError("Recorded COD dependency mismatch: " + json.dumps(differences, sort_keys=True))
    executable = Path(sys.executable).resolve(strict=True)
    return {"python_version": sys.version, "python_executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
            "packages": versions, "original_python_version": profile["original_python_version"],
            "original_python_patch_matches": sys.version.split()[0] == profile["original_python_version"],
            "full_historical_runtime_byte_identity_proven": False,
            "external_connections": "disabled", "charge_mode": "cpu-reference",
            "threads": 1, "gpus_visible": False}


if __name__ == "sitecustomize":
    disable_network()
elif __name__ == "__main__":
    disable_network()
    profile = json.loads(Path(sys.argv[1]).read_text())
    print(json.dumps(inspect_runtime(profile, skip_charges=sys.argv[2] == "1"), sort_keys=True))
