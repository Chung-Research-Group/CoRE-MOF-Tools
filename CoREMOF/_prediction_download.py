"""Staged, no-overwrite download of an explicitly requested model asset."""

import hashlib
import os
from pathlib import Path
import re
import tempfile


def _check_file(path, expected_sha256):
    if path.is_symlink() or any(parent.is_symlink() for parent in path.parents):
        raise ValueError("Model cache paths must not contain symlinks")
    if os.path.lexists(path):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Model cache is not a nonempty regular file: {path}")
        if expected_sha256 is not None:
            digest = hashlib.sha256()
            with path.open("rb") as stream:
                for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                    digest.update(chunk)
            if digest.hexdigest() != expected_sha256:
                raise ValueError(f"Model checksum mismatch: {path}")
        return True
    return False


def download_model_file(url, save_path, *, expected_sha256=None):
    """Return whether a new asset was published; never replace an existing one.

    A checksum must come from a separately trusted record. Without one, this
    function ensures complete staged publication, not scientific model identity.
    It never deserializes the downloaded file.
    """
    if expected_sha256 is not None:
        if not isinstance(expected_sha256, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", expected_sha256):
            raise ValueError("expected_sha256 must contain 64 hexadecimal characters")
        expected_sha256 = expected_sha256.lower()
    destination = Path(save_path).absolute()
    if _check_file(destination, expected_sha256):
        return False
    import requests

    destination.parent.mkdir(parents=True, exist_ok=True)
    if _check_file(destination, expected_sha256):
        return False
    descriptor, name = tempfile.mkstemp(prefix=".coremof-model-", dir=destination.parent)
    temporary = Path(name)
    try:
        digest = hashlib.sha256()
        size = 0
        with os.fdopen(descriptor, "wb") as stream:
            response = requests.get(url, stream=True, timeout=60)
            try:
                response.raise_for_status()
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        stream.write(chunk)
                        digest.update(chunk)
                        size += len(chunk)
                if size == 0:
                    raise ValueError("Downloaded model is empty")
                if expected_sha256 is not None and digest.hexdigest() != expected_sha256:
                    raise ValueError("Downloaded model checksum mismatch")
                stream.flush()
                os.fsync(stream.fileno())
            finally:
                response.close()
        # Hard-link publication cannot replace a concurrent writer's file.
        try:
            os.link(temporary, destination)
        except FileExistsError:
            if not _check_file(destination, expected_sha256):
                raise FileNotFoundError("Concurrent model cache disappeared before verification")
            return False
        return True
    finally:
        temporary.unlink()
