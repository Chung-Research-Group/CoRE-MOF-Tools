"""Private, identity-bound state seals for authoritative package objects.

Python objects are mutable even when their public containers are tuples or
mapping proxies: ``object.__setattr__`` and ``dataclasses.replace`` can still
manufacture a new object carrying copied private attributes.  Authority
therefore never follows an attribute.  It is held in a module-private
identity registry and is accepted only after a type-sensitive fingerprint of
the complete authority-relevant state is recomputed at the consuming
boundary.

This module is deliberately internal.  The seals authenticate objects made
by validated package factories; they are not a signature format and are not
serializable delegation credentials.
"""

from __future__ import annotations

import hashlib
import math
from decimal import Decimal
from pathlib import PurePath
import struct
from typing import Mapping
import weakref


class AuthorityStateError(ValueError):
    """Raised when a factory-authenticated object no longer matches its seal."""


_SEALED_OBJECT_MARKER = object()


def _type_name(value: object) -> str:
    kind = type(value)
    return "{}.{}".format(kind.__module__, kind.__qualname__)


def _feed_bytes(digest, value: bytes) -> None:
    digest.update(struct.pack(">Q", len(value)))
    digest.update(value)


def _feed_value(digest, value: object) -> None:
    """Feed one value without coercing its type, order, or representation."""

    if value is None:
        digest.update(b"N;")
    elif type(value) is bool:
        digest.update(b"B1;" if value else b"B0;")
    elif type(value) is int:
        digest.update(b"I")
        digest.update(str(value).encode("ascii"))
        digest.update(b";")
    elif type(value) is float:
        if not math.isfinite(value):
            raise AuthorityStateError(
                "authority-relevant state contains a non-finite float"
            )
        digest.update(b"F")
        digest.update(value.hex().encode("ascii"))
        digest.update(b";")
    elif type(value) is Decimal:
        if not value.is_finite():
            raise AuthorityStateError(
                "authority-relevant state contains a non-finite decimal"
            )
        sign, digits, exponent = value.as_tuple()
        digest.update(b"D")
        digest.update(str(sign).encode("ascii"))
        digest.update(b":")
        digest.update("".join(str(digit) for digit in digits).encode("ascii"))
        digest.update(b":")
        digest.update(str(exponent).encode("ascii"))
        digest.update(b";")
    elif type(value) is str:
        digest.update(b"S")
        _feed_bytes(digest, value.encode("utf-8"))
    elif type(value) is bytes:
        digest.update(b"Y")
        _feed_bytes(digest, value)
    elif isinstance(value, PurePath):
        # Require the exact path class as recorded above and preserve the path
        # spelling; callers separately enforce resolved release roots.
        digest.update(b"P")
        _feed_bytes(digest, _type_name(value).encode("utf-8"))
        _feed_bytes(digest, str(value).encode("utf-8"))
    elif isinstance(value, Mapping):
        digest.update(b"M")
        _feed_bytes(digest, _type_name(value).encode("utf-8"))
        digest.update(struct.pack(">Q", len(value)))
        for key, item in value.items():
            _feed_value(digest, key)
            _feed_value(digest, item)
    elif type(value) in (tuple, list):
        digest.update(b"T" if type(value) is tuple else b"L")
        digest.update(struct.pack(">Q", len(value)))
        for item in value:
            _feed_value(digest, item)
    else:
        raise AuthorityStateError(
            "authority-relevant state contains unsupported type {}".format(
                _type_name(value)
            )
        )


def state_fingerprint(value: object) -> str:
    """Return a type-, order-, and content-sensitive SHA-256 fingerprint."""

    digest = hashlib.sha256()
    digest.update(b"coremof-authority-state/1\0")
    _feed_value(digest, value)
    return digest.hexdigest()


class IdentitySealRegistry:
    """Keep non-transferable seals keyed by exact live object identity."""

    def __init__(self, name: str) -> None:
        self._name = name
        self._entries = {}

    def register(self, value: object, fingerprint: str, context=None) -> None:
        key = id(value)

        def discard(reference, *, registry=self, object_id=key) -> None:
            current = registry._entries.get(object_id)
            if current is not None and current[0] is reference:
                registry._entries.pop(object_id, None)

        reference = weakref.ref(value, discard)
        self._entries[key] = (reference, fingerprint, context)
        object.__setattr__(value, "_authority_generation_marker", _SEALED_OBJECT_MARKER)

    def entry(self, value: object):
        current = self._entries.get(id(value))
        if current is not None and current[0]() is value:
            return current[1], current[2]
        try:
            own_state = object.__getattribute__(value, "__dict__")
        except (AttributeError, TypeError):
            own_state = {}
        if own_state.get("_authority_generation_marker") is _SEALED_OBJECT_MARKER:
            raise AuthorityStateError(
                "{} authority does not transfer through copying, pickling, "
                "replacement, or manual construction".format(self._name)
            )
        return None

    def require(self, value: object, fingerprint: str):
        current = self.entry(value)
        if current is None:
            raise AuthorityStateError(
                "{} was not produced by its authenticated package factory".format(
                    self._name
                )
            )
        expected, context = current
        if type(fingerprint) is not str or fingerprint != expected:
            raise AuthorityStateError(
                "{} changed after authentication; rebuild it from validated "
                "inputs".format(self._name)
            )
        return context

    def is_registered(self, value: object) -> bool:
        return self.entry(value) is not None


def reject_sealed_copy(value: object, operation: str) -> None:
    """Reject copy/pickle operations for authority-bearing objects."""

    try:
        own_state = object.__getattribute__(value, "__dict__")
    except (AttributeError, TypeError):
        own_state = {}
    if own_state.get("_authority_generation_marker") is _SEALED_OBJECT_MARKER:
        raise TypeError(
            "authenticated CoREMOF objects cannot be {}; rebuild from the "
            "validated factory instead".format(operation)
        )
