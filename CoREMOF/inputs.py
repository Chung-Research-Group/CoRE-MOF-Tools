"""Shared input-path handling for public CoREMOF workflows."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Tuple, Union
import unicodedata


PathInput = Union[str, os.PathLike]


def _is_cif(path: Path) -> bool:
    return path.suffix.casefold() == ".cif"


def _basename_key(path: Path) -> str:
    """Return a portable comparison key without changing the returned path."""

    return unicodedata.normalize("NFKC", path.stem).casefold()


def resolve_cif_inputs(input_path: PathInput) -> Tuple[Path, ...]:
    """Resolve one CIF file or a directory of direct CIF children.

    A file must be a regular, non-symlink path whose suffix is ``.cif``
    (case-insensitive). A directory is scanned non-recursively; its matching
    children are returned in lexical filename order. Matching symlinks and
    special files are rejected rather than followed or silently skipped.

    The return type is always a tuple, including for a single file. Basenames
    are compared after Unicode NFKC normalization and case folding solely to
    prevent ambiguous output identifiers; the input paths and CIF bytes are
    not modified.
    """

    path = Path(input_path).expanduser()
    if path.is_symlink():
        raise ValueError(f"CIF input must not be a symlink: {path}")
    if not path.exists():
        raise FileNotFoundError(f"CIF input does not exist: {path}")

    if path.is_file():
        if not _is_cif(path):
            raise ValueError(f"CIF input file must have a .cif suffix: {path}")
        return (path,)

    if not path.is_dir():
        raise ValueError(f"CIF input must be a regular file or directory: {path}")

    matches = []
    for child in path.iterdir():
        if not _is_cif(child):
            continue
        if child.is_symlink():
            raise ValueError(f"CIF directory contains a symlink: {child}")
        if not child.is_file():
            raise ValueError(f"CIF directory contains a special CIF path: {child}")
        matches.append(child)

    matches.sort(key=lambda item: item.name)
    if not matches:
        raise ValueError(f"CIF input directory contains no direct CIF files: {path}")

    seen = {}
    for child in matches:
        key = _basename_key(child)
        previous = seen.get(key)
        if previous is not None:
            raise ValueError(
                "CIF input basenames collide after Unicode NFKC normalization "
                f"and case folding: {previous.name!r}, {child.name!r}"
            )
        seen[key] = child

    return tuple(matches)


def collect_cifs(input_path: PathInput) -> Tuple[Path, ...]:
    """Alias for :func:`resolve_cif_inputs`."""

    return resolve_cif_inputs(input_path)


__all__ = ["PathInput", "collect_cifs", "resolve_cif_inputs"]
