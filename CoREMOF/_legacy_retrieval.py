"""Filesystem guards for the historical SI/CSD retrieval interfaces.

No optional chemistry libraries, network requests or licensed imports belong
here. This module does not establish redistribution rights or release identity.
"""

from __future__ import annotations

from contextlib import ExitStack
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import stat
import tempfile
import zipfile

from ._transactions import publish_file_bundle


def validate_refcode(refcode):
    """Validate a safe filename token, not the complete CSD identifier grammar."""
    if not isinstance(refcode, str) or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", refcode) is None:
        raise ValueError("refcode must be a nonempty filename token containing letters, digits, '-' or '_'")
    return refcode


def _member_path(info):
    name = info.filename
    trimmed = name[:-1] if name.endswith("/") else name
    parts = trimmed.split("/")
    kind = stat.S_IFMT(info.external_attr >> 16)
    expected_kind = stat.S_IFDIR if info.is_dir() else stat.S_IFREG
    if (
        info.orig_filename != name
        or not trimmed
        or any(part in ("", ".", "..") for part in parts)
        or "\\" in name
        or ":" in name
        or any(ord(char) < 32 for char in name)
        or kind not in (0, expected_kind)
        or (info.is_dir() and info.file_size != 0)
    ):
        raise ValueError("Unsafe or unsupported ZIP member: {!r}".format(name))
    return Path(*PurePosixPath(trimmed).parts)


def _validated_members(archive):
    members = []
    seen = set()
    files = set()
    directories = set()
    for info in archive.infolist():
        relative = _member_path(info)
        if relative in seen:
            raise ValueError("Duplicate ZIP member: {}".format(relative))
        seen.add(relative)
        (directories if info.is_dir() else files).add(relative)
        members.append((info, relative))
        directories.update(tuple(relative.parents)[:-1])
    if files & directories:
        raise ValueError("ZIP member is both a file and a directory")
    return members


def validate_cache(path, file_name):
    """Reject malformed downloaded JSON/ZIP content before caching it."""
    if file_name.endswith(".json"):
        with Path(path).open(encoding="utf-8") as stream:
            value = json.load(stream)
        if not isinstance(value, dict):
            raise ValueError("Legacy metadata must contain a JSON object")
    elif file_name.endswith(".zip"):
        with zipfile.ZipFile(path) as archive:
            _validated_members(archive)
            failed = archive.testzip()
            if failed is not None:
                raise zipfile.BadZipFile("CRC failure in {}".format(failed))
    else:
        raise ValueError("Unsupported legacy cache format: {}".format(file_name))


def _directory_path(path):
    # Do not resolve symlinks before checking them, including a dangling link.
    path = Path(os.path.abspath(os.fspath(path)))
    for candidate in (path, *path.parents):
        if candidate.is_symlink():
            raise ValueError("Refusing a symlink output directory: {}".format(candidate))
        if candidate.exists() and not candidate.is_dir():
            raise NotADirectoryError(str(candidate))
    return path


def _check_target(path, overwrite):
    _directory_path(path.parent)
    if path.is_symlink():
        raise ValueError("Refusing a symlink output file: {}".format(path))
    if path.exists():
        if not path.is_file():
            raise ValueError("Output is not an ordinary file: {}".format(path))
        if not overwrite:
            raise FileExistsError("Refusing to overwrite {}; use overwrite=True explicitly".format(path))


def check_output_file(path, *, overwrite):
    if type(overwrite) is not bool:
        raise TypeError("overwrite must be a boolean")
    path = Path(os.path.abspath(os.fspath(path)))
    _check_target(path, overwrite)
    return path


def _publish_writers(writers, directories, *, overwrite):
    """Stage every file, then publish only the selected files with rollback."""
    if type(overwrite) is not bool:
        raise TypeError("overwrite must be a boolean")
    directories = set(directories)
    for target, _ in writers:
        _check_target(target, overwrite)
        directories.add(target.parent)
    for directory in directories:
        _directory_path(directory)
    if not writers and not directories:
        return
    staging_parent = min(directories, key=lambda path: len(path.parts))
    while not staging_parent.exists():
        staging_parent = staging_parent.parent
    staging = Path(tempfile.mkdtemp(prefix=".coremof-download-", dir=staging_parent))
    created = []
    preserve = False
    try:
        staged = []
        targets = []
        for index, (target, writer) in enumerate(writers):
            source = staging / (str(index) + ".data")
            with source.open("xb") as stream:
                writer(stream)
            staged.append(source)
            targets.append(target)
        # No destination files/directories are created until every source has
        # been read successfully (including the ZIP CRC check on stream close).
        expanded = set(directories)
        for directory in directories:
            expanded.update(directory.parents)
        for directory in sorted(expanded, key=lambda path: (len(path.parts), str(path))):
            _directory_path(directory)
            if not directory.exists():
                directory.mkdir()
                created.append(directory)
        for target in targets:
            _check_target(target, overwrite)
        if staged:
            publish_file_bundle(staged, targets, overwrite=overwrite)
    except BaseException as exc:
        preserve = bool(getattr(exc, "coremof_preserved_staging_directory", None))
        for directory in reversed(created):
            try:
                directory.rmdir()  # Only our newly created, still-empty directories.
            except OSError:
                pass
        raise
    finally:
        if not preserve:
            shutil.rmtree(staging)


def extract_archives(sources, output_folder, *, overwrite=False, ensure_directories=()):
    """Extract whole archives or named entries as one protected file bundle.

    ``sources`` contains (ZIP path, entry or None) pairs. A missing selected entry
    retains the legacy no-op behavior. Unrelated destination files are untouched.
    """
    if type(overwrite) is not bool:
        raise TypeError("overwrite must be a boolean")
    destination = _directory_path(output_folder)
    directories = {destination / item for item in ensure_directories}
    writers = []
    file_paths = set()
    with ExitStack() as stack:
        for path, entry in sources:
            archive = stack.enter_context(zipfile.ZipFile(path))
            for info, relative in _validated_members(archive):
                if entry is not None and entry != info.filename:
                    continue
                target = destination / relative
                if info.is_dir():
                    directories.add(target)
                    continue
                if target in file_paths:
                    raise ValueError("Duplicate output across ZIP archives: {}".format(relative))
                file_paths.add(target)
                directories.update(
                    parent for parent in target.parents
                    if parent == destination or destination in parent.parents
                )

                def write_member(stream, archive=archive, info=info):
                    with archive.open(info) as source:
                        shutil.copyfileobj(source, stream)

                writers.append((target, write_member))
        if file_paths & directories:
            raise ValueError("ZIP outputs have a file/directory conflict")
        _publish_writers(writers, directories, overwrite=overwrite)


def publish_text(path, text, *, overwrite=False):
    """Publish the exact exported text without altering scientific content."""
    target = check_output_file(path, overwrite=overwrite)
    content = text.encode("utf-8")
    _publish_writers([(target, lambda stream: stream.write(content))], {target.parent}, overwrite=overwrite)
