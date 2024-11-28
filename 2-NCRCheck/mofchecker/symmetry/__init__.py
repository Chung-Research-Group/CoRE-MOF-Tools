# -*- coding: utf-8 -*-
"""Analyze symmetrized structures and hash them."""
import functools
from typing import Union

from pymatgen.core import IStructure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.structure import SymmetrizedStructure

from .hash import hash_symmetrized_structure


@functools.lru_cache(maxsize=2, typed=False)
def get_symmetrized_structure(structure: IStructure) -> SymmetrizedStructure:
    """Construct a SymmetrizedStructure.

    That is. a structure where the spacegroup and symmetry operations are defined.

    Args:
        structure (IStructure): structure to symmetrize

    Returns:
        SymmetrizedStructure: symmetrized structure
    """
    return SpacegroupAnalyzer(structure).get_symmetrized_structure()


def symmetrize_if_not_symmetrized(
    structure: Union[IStructure, SymmetrizedStructure]
) -> SymmetrizedStructure:
    """Return a symmetrized structure."""
    if not isinstance(structure, SymmetrizedStructure):
        structure = get_symmetrized_structure(structure)
    return structure


def get_spacegroup_symbol_and_number(structure: Union[IStructure, SymmetrizedStructure]) -> dict:
    """Return a dict with spacegroup symbol and number."""
    structure = symmetrize_if_not_symmetrized(structure)
    return {
        "symbol": structure.spacegroup.int_symbol,
        "number": structure.spacegroup.int_number,
    }


def get_symmetry_hash(
    structure: Union[IStructure, SymmetrizedStructure], tight: bool = False
) -> str:
    """Hashes the symmetrical positions of the SymmetrizedStructure.

    The tight setting also considers the ordering,
    otherwise only the number and identity
    of the elements is considered.

    Args:
        structure (Union[Structure, IStructure, SymmetrizedStructure]):
            A structure for which the symmetry hash is calculated
        tight (bool): If True, also consider the ordering. Defaults to False.

    Returns:
        str: hash
    """
    structure = symmetrize_if_not_symmetrized(structure)
    return hash_symmetrized_structure(structure, tight)
