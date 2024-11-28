# -*- coding: utf-8 -*-
"""Compute hashes for symmetrized structures based on the Wyckoff letters."""
import base64
import hashlib
from collections import Counter

from pymatgen.symmetry.structure import SymmetrizedStructure


def make_sha256_hash(tupl: tuple) -> str:
    """Based on https://stackoverflow.com/a/42151923."""
    hasher = hashlib.sha256()
    hasher.update(repr(tupl).encode())
    return base64.b64encode(hasher.digest()).decode()


def make_hashable(counter: Counter) -> tuple:
    """Based on https://stackoverflow.com/a/42151923."""
    return tuple(sorted((key, value) for key, value in counter.items()))


def hash_symmetrized_structure(
    symmetrized_structure: SymmetrizedStructure, tight: bool = False
) -> str:
    """Run the hashing.

    Args:
        symmetrized_structure (SymmetrizedStructure): A structure object
            that has the Wyckoff letters as property/attribute
        tight (bool): If True, also consider the ordering
            of the Wyckoff letters. Defaults to False.

    Returns:
        str: hash
    """
    if tight:
        return "".join(symmetrized_structure.wyckoff_letters) + str(
            symmetrized_structure.spacegroup.int_number
        )
    wyckoff_letter_counts = tuple(set(symmetrized_structure.wyckoff_letters))
    return make_sha256_hash(wyckoff_letter_counts) + str(
        symmetrized_structure.spacegroup.int_number
    )
