"""Retrieve an exact release from a caller-supplied, checksum-pinned catalog.

No hosted catalog or dataset access entitlement is bundled with this API.
HTTP(S) and local file resources use the same size and SHA-256 checks. The
existing release loader validates the staged data before it is published.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import tempfile
import unicodedata
from urllib.error import HTTPError, URLError
from urllib.parse import urljoin, urlsplit
from urllib.request import Request, urlopen

from . import __version__
from ._transactions import publish_directory
from .dataset import CoREMOFDataset


CATALOG_SCHEMA = "coremof-release-catalog/1.0"
RECEIPT_PATH = "manifests/retrieval_receipt.json"
MAX_CATALOG_BYTES = 64 * 1024 * 1024


class RetrievalError(ValueError):
    """The requested release could not be retrieved and verified."""


def _digest(value, label):
    if not isinstance(value, str) or re.fullmatch(r"[0-9a-f]{64}", value) is None:
        raise RetrievalError("{} must be a lowercase full SHA-256".format(label))
    return value


def _name(value, label):
    if not isinstance(value, str) or not value or value != value.strip():
        raise RetrievalError("{} must be an exact nonblank string".format(label))
    return value


def _url(value):
    _name(value, "resource URL")
    try:
        parsed = urlsplit(value)
        if parsed.scheme not in ("http", "https", "file"):
            raise RetrievalError("resource URLs must use HTTP, HTTPS or file")
        if parsed.username is not None or parsed.password is not None:
            raise RetrievalError("credentials embedded in resource URLs are unsupported")
        if parsed.fragment:
            raise RetrievalError("resource URL fragments are unsupported")
        if parsed.scheme == "file":
            if parsed.netloc not in ("", "localhost") or parsed.query:
                raise RetrievalError("file resources must be local and have no query")
        elif not parsed.hostname:
            raise RetrievalError("HTTP(S) resource URL requires a hostname")
    except ValueError as error:
        if isinstance(error, RetrievalError):
            raise
        raise RetrievalError("invalid resource URL") from None
    return value


def _catalog_url(source):
    if isinstance(source, os.PathLike):
        return Path(source).expanduser().resolve().as_uri()
    _name(source, "catalog path or URL")
    if urlsplit(source).scheme:
        return _url(source)
    return Path(source).expanduser().resolve().as_uri()


def _open_resource(url, timeout):
    # Exceptions deliberately omit URLs: signed query strings can be secrets.
    try:
        response = urlopen(Request(_url(url), headers={"Accept-Encoding": "identity"}), timeout=timeout)
        try:
            _url(response.geturl())
        except ValueError:
            response.close()
            raise
        return response
    except HTTPError as error:
        code = error.code
        error.close()
        raise RetrievalError("resource request failed with HTTP {}".format(code)) from None
    except (URLError, OSError):
        raise RetrievalError("resource could not be read") from None


def _parse_catalog(data):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise RetrievalError("duplicate catalog JSON key")
            result[key] = value
        return result

    def nonfinite(value):
        raise RetrievalError("nonfinite catalog JSON value")

    try:
        catalog = json.loads(data.decode("utf-8"), object_pairs_hook=unique, parse_constant=nonfinite)
    except (UnicodeError, json.JSONDecodeError):
        raise RetrievalError("catalog must be UTF-8 JSON") from None
    if (not isinstance(catalog, dict) or set(catalog) != {"schema_version", "releases"}
            or catalog["schema_version"] != CATALOG_SCHEMA
            or not isinstance(catalog["releases"], dict) or not catalog["releases"]):
        raise RetrievalError("unsupported release catalog schema")
    for version in catalog["releases"]:
        _name(version, "catalog release version")
    return catalog


def _selected_files(catalog, version, base_url):
    if version not in catalog["releases"]:
        raise RetrievalError("requested version is absent from the pinned catalog")
    release = catalog["releases"][version]
    if (not isinstance(release, dict) or set(release) != {"files"}
            or not isinstance(release["files"], list) or not release["files"]):
        raise RetrievalError("selected release must declare a nonempty files list")
    files, paths = [], {}
    for item in release["files"]:
        if not isinstance(item, dict) or set(item) != {"path", "url", "sha256", "size_bytes"}:
            raise RetrievalError("each file requires path, url, sha256 and size_bytes")
        name = _name(item["path"], "release file path")
        path = PurePosixPath(name)
        if (path.is_absolute() or path.as_posix() != name or ".." in path.parts
                or name == "." or "\\" in name or ":" in name or "\x00" in name):
            raise RetrievalError("release file paths must be safe canonical relative paths")
        key = unicodedata.normalize("NFKC", name).casefold()
        if key == RECEIPT_PATH.casefold() or key in paths:
            raise RetrievalError("duplicate, ambiguous or reserved release file path")
        paths[key] = name
        size = item["size_bytes"]
        if type(size) is not int or size < 0:
            raise RetrievalError("file size_bytes must be a nonnegative integer")
        source = _name(item["url"], "file URL")
        files.append({"path": name, "url": _url(urljoin(base_url, source)),
                      "sha256": _digest(item["sha256"], "file sha256"), "size_bytes": size})
    # Reject a file whose name would also have to be a directory, including the
    # reserved receipt's parent. No symlinks or archive members are extracted.
    all_paths = set(paths).union({RECEIPT_PATH.casefold()})
    for name in all_paths:
        if any(parent.as_posix() in all_paths for parent in PurePosixPath(name).parents
               if parent.as_posix() != "."):
            raise RetrievalError("release file paths conflict with required directories")
    return tuple(sorted(files, key=lambda item: item["path"]))


def _download_file(item, destination, timeout):
    destination.parent.mkdir(parents=True, exist_ok=True)
    digest, length = hashlib.sha256(), 0
    with _open_resource(item["url"], timeout) as response, destination.open("xb") as stream:
        while True:
            chunk = response.read(min(1024 * 1024, item["size_bytes"] - length + 1))
            if not chunk:
                break
            length += len(chunk)
            if length > item["size_bytes"]:
                raise RetrievalError("download exceeds declared file size: {}".format(item["path"]))
            digest.update(chunk)
            stream.write(chunk)
    if length != item["size_bytes"] or digest.hexdigest() != item["sha256"]:
        raise RetrievalError("download size or checksum mismatch: {}".format(item["path"]))


def fetch_release(catalog, version, destination, *, catalog_sha256,
                  verify_cif_files=False, timeout=60.0, max_total_bytes=None):
    """Retrieve an exact catalog version into a new validated local directory.

    ``catalog`` is a local JSON path or HTTP(S)/file URL. Its expected SHA-256
    must come from a trusted channel; there is no implicit latest version or
    built-in hosted catalog. File URLs can be relative to the catalog URL.
    Every listed file is checked against its declared byte length and full
    SHA-256. The package release loader checks schemas and the exact requested
    dataset version before publication. ``verify_cif_files=True`` additionally
    requires and hashes every CIF named by the release's CIF manifest.

    Return the published Path. Existing destinations are never intentionally
    replaced. Publication uses the package's shared atomic directory writer;
    its NFS fallback serializes cooperating writers. Failed staging is removed.
    The receipt contains logical file paths and hashes, but no source URLs,
    signed query strings or private local catalog paths.
    """
    expected = _digest(catalog_sha256, "catalog_sha256")
    _name(version, "version")
    if type(verify_cif_files) is not bool:
        raise TypeError("verify_cif_files must be a boolean")
    if (type(timeout) not in (int, float) or not math.isfinite(timeout) or timeout <= 0):
        raise ValueError("timeout must be a finite positive number")
    if max_total_bytes is not None and (type(max_total_bytes) is not int or max_total_bytes < 0):
        raise ValueError("max_total_bytes must be a nonnegative integer or None")
    target = Path(destination).expanduser()
    if os.path.lexists(target):
        raise FileExistsError("release destination already exists")
    source_url = _catalog_url(catalog)
    with _open_resource(source_url, timeout) as response:
        resolved_url = response.geturl()
        data = response.read(MAX_CATALOG_BYTES + 1)
    if len(data) > MAX_CATALOG_BYTES:
        raise RetrievalError("catalog exceeds the 64 MiB size limit")
    if hashlib.sha256(data).hexdigest() != expected:
        raise RetrievalError("catalog checksum mismatch")
    files = _selected_files(_parse_catalog(data), version, resolved_url)
    total = sum(item["size_bytes"] for item in files)
    if max_total_bytes is not None and total > max_total_bytes:
        raise RetrievalError("selected release exceeds max_total_bytes")
    implementation_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    target.parent.mkdir(parents=True, exist_ok=True)
    target = target.parent.resolve() / target.name
    with tempfile.TemporaryDirectory(prefix=".coremof-retrieval-", dir=target.parent) as temporary:
        staged = Path(temporary) / "release"
        staged.mkdir()
        for item in files:
            _download_file(item, staged / item["path"], timeout)
        dataset = CoREMOFDataset.from_release(staged, verify_cif_files=verify_cif_files)
        if dataset.dataset_version != version:
            raise RetrievalError("downloaded dataset version differs from the selected catalog version")
        receipt = {
            "schema_version": "coremof-release-retrieval-receipt/1.0",
            "package_version": __version__, "retrieval_source_sha256": implementation_hash,
            "catalog_sha256": expected, "catalog_size_bytes": len(data), "selected_version": version,
            "file_count": len(files), "downloaded_bytes": total,
            "structure_count": len(dataset), "cif_files_verified": dataset.cif_files_verified,
            "release_input_hashes": dict(dataset.input_hashes),
            "files": [{key: item[key] for key in ("path", "size_bytes", "sha256")} for item in files],
        }
        receipt_path = staged / RECEIPT_PATH
        receipt_path.parent.mkdir(parents=True, exist_ok=True)
        receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n",
                                encoding="utf-8")
        publish_directory(staged, target, overwrite=False)
    return target


__all__ = ["CATALOG_SCHEMA", "RECEIPT_PATH", "RetrievalError", "fetch_release"]
