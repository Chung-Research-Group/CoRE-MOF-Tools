"""Small filesystem publication primitives shared by additive workflows."""

from __future__ import annotations

import ctypes
import errno
import fcntl
import os
from pathlib import Path
import shutil
import sys
from typing import Sequence


def _same_file(left: Path, right: Path) -> bool:
    try:
        return (
            left.exists() and right.exists() and os.path.samefile(str(left), str(right))
        )
    except OSError:
        return False


def _restore_detached(detached: Path, target: Path) -> bool:
    """Restore without overwriting a newer generation; retain on failure."""

    try:
        os.link(str(detached), str(target))
    except OSError:
        return False
    detached.unlink()
    return True


def publish_file_bundle(
    staged: Sequence[Path],
    targets: Sequence[Path],
    *,
    overwrite: bool,
) -> None:
    """Publish ordinary files with create-if-absent and identity-safe rollback."""

    if type(overwrite) is not bool:
        raise TypeError("overwrite must be a boolean")
    staged = tuple(staged)
    targets = tuple(targets)
    if len(staged) != len(targets) or not staged:
        raise ValueError(
            "staged and target file sequences must have equal nonzero length"
        )
    staging = staged[0].parent
    if any(path.parent != staging for path in staged):
        raise ValueError("all staged files must share one private directory")
    if not overwrite and any(path.exists() for path in targets):
        raise FileExistsError("refusing to overwrite an existing output bundle")

    backups = {}
    identities = {}
    published = []
    preserve_staging = False
    try:
        if overwrite:
            for target in targets:
                if target.exists():
                    backup = staging / (target.name + ".previous")
                    try:
                        os.link(str(target), str(backup))
                    except OSError:
                        shutil.copy2(str(target), str(backup))
                    backups[target] = backup
        for source, target in zip(staged, targets):
            identity = staging / (target.name + ".published")
            os.link(str(source), str(identity))
            identities[target] = identity
            if overwrite:
                os.replace(str(source), str(target))
            else:
                # Hard-link publication is an atomic create-if-absent operation.
                os.link(str(source), str(target))
            published.append(target)
    except BaseException:
        for target in reversed(published):
            detached = staging / (target.name + ".rollback")
            try:
                os.replace(str(target), str(detached))
            except FileNotFoundError:
                detached = None
            if detached is not None:
                if _same_file(detached, identities[target]):
                    detached.unlink()
                    backup = backups.get(target)
                    if backup is not None:
                        try:
                            os.link(str(backup), str(target))
                        except OSError:
                            preserve_staging = True
                elif not _restore_detached(detached, target):
                    preserve_staging = True
        if preserve_staging:
            # Retaining a private backup/quarantine is preferable to deleting a
            # generation that could not be restored without overwriting another.
            setattr(
                sys.exc_info()[1],
                "coremof_preserved_staging_directory",
                str(staging),
            )
        raise


def _rename_directory_noreplace(source: Path, target: Path) -> None:
    """Atomically rename a directory only when the target does not exist."""

    if sys.platform.startswith("linux"):
        libc = ctypes.CDLL(None, use_errno=True)
        renameat2 = getattr(libc, "renameat2", None)
        if renameat2 is None:
            raise OSError(errno.ENOSYS, "renameat2 is unavailable")
        renameat2.argtypes = (
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_uint,
        )
        renameat2.restype = ctypes.c_int
        result = renameat2(
            -100,
            os.fsencode(str(source)),
            -100,
            os.fsencode(str(target)),
            1,  # RENAME_NOREPLACE
        )
        if result == 0:
            return
        error_number = ctypes.get_errno()
        if error_number in (errno.EEXIST, errno.ENOTEMPTY):
            raise FileExistsError(error_number, os.strerror(error_number), str(target))
        if error_number in (
            errno.EINVAL,
            errno.ENOSYS,
            errno.ENOTSUP,
            errno.EOPNOTSUPP,
        ):
            # Some shared filesystems (including deployed NFS variants) reject
            # RENAME_NOREPLACE even though ordinary same-directory rename is
            # atomic.  Serialize cooperating writers with a crash-releasing
            # advisory lock, recheck the destination while holding it, and
            # then use that atomic rename.  A non-cooperating process can still
            # race on such a filesystem; a pre-existing destination is never
            # deliberately replaced.
            _locked_directory_rename_noreplace(source, target)
            return
        raise OSError(error_number, os.strerror(error_number), str(target))
    if sys.platform == "darwin":  # pragma: no cover - Linux CI/HPC path
        libc = ctypes.CDLL(None, use_errno=True)
        renamex_np = getattr(libc, "renamex_np", None)
        if renamex_np is None:
            raise OSError(errno.ENOSYS, "renamex_np is unavailable")
        renamex_np.argtypes = (ctypes.c_char_p, ctypes.c_char_p, ctypes.c_uint)
        renamex_np.restype = ctypes.c_int
        if (
            renamex_np(os.fsencode(str(source)), os.fsencode(str(target)), 0x00000004)
            == 0
        ):
            return
        error_number = ctypes.get_errno()
        if error_number in (errno.EEXIST, errno.ENOTEMPTY):
            raise FileExistsError(error_number, os.strerror(error_number), str(target))
        raise OSError(error_number, os.strerror(error_number), str(target))
    # On Windows os.rename fails when the destination exists.
    if os.name == "nt":  # pragma: no cover - Linux CI/HPC path
        os.rename(str(source), str(target))
        return
    raise OSError(
        errno.ENOTSUP,
        "atomic no-replace directory publication is unsupported on this platform",
        str(target),
    )


def _locked_directory_rename_noreplace(source: Path, target: Path) -> None:
    """NFS-compatible no-replace fallback serialized across our writers."""

    lock_path = target.parent / (".{}.publish.lock".format(target.name))
    descriptor = os.open(str(lock_path), os.O_CREAT | os.O_RDWR, 0o600)
    try:
        fcntl.flock(descriptor, fcntl.LOCK_EX)
        if os.path.lexists(str(target)):
            raise FileExistsError(errno.EEXIST, os.strerror(errno.EEXIST), str(target))
        os.rename(str(source), str(target))
    finally:
        os.close(descriptor)


def publish_directory(staging: Path, target: Path, *, overwrite: bool) -> Path:
    """Publish one fully staged directory without a no-overwrite check/rename race."""

    if type(overwrite) is not bool:
        raise TypeError("overwrite must be a boolean")
    if not overwrite:
        _rename_directory_noreplace(staging, target)
        return target
    backup = None
    if target.exists():
        backup = target.parent / (".{}.previous.{}".format(target.name, os.getpid()))
        if backup.exists():
            raise FileExistsError("stale output backup exists: {}".format(backup))
        os.replace(str(target), str(backup))
    try:
        os.replace(str(staging), str(target))
    except BaseException:
        if backup is not None and backup.exists() and not target.exists():
            os.replace(str(backup), str(target))
        raise
    if backup is not None:
        shutil.rmtree(backup)
    return target


__all__ = ["publish_directory", "publish_file_bundle"]
