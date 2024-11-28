# -*- coding: utf-8 -*-
"""Utility function for getting the indices for certain atoms in the structure."""
import functools
from typing import Union

import pymatgen
from pymatgen.core import IStructure, Structure

from ..data import _get_vdw_radius
from ...definitions import METALS


def _vdw_radius_neighbors(structure, site_index, tolerance: float = 1.5):
    elem = str(structure[site_index].specie)
    radius = _get_vdw_radius(elem)
    return structure.get_neighbors(structure[site_index], tolerance * radius)


def is_metal(site: pymatgen.core.Site) -> bool:
    """Return True if the site is a metal.

    Considers transition metal, lanthanide, actinide,
    or Al, Ga, In, Tl, Ge, Sn, Pb, Sb, Bi, Po

    Args:
        site: pymatgen.core.Site

    Returns:
        bool: True if the site is a metal
    """
    if str(site.specie) in METALS:
        return True
    return False


@functools.lru_cache(maxsize=2, typed=False)
def _get_indices(immutable_structure: IStructure) -> dict:
    return {
        "c": _get_c_indices(immutable_structure),
        "h": _get_h_indices(immutable_structure),
        "n": _get_n_indices(immutable_structure),
        "o": _get_o_indices(immutable_structure),
        "metal": _get_metal_indices(immutable_structure),
        "rare_earth": _get_rare_earth_indices(immutable_structure),
        "alkali_alkaline": _get_alkali_alkaline_indices(immutable_structure),
        "halogen": _get_halogen_indices(immutable_structure),
    }


def _get_c_indices(structure):
    return [i for i, species in enumerate(structure.species) if str(species) == "C"]


def _get_h_indices(structure):
    return [i for i, species in enumerate(structure.species) if str(species) == "H"]


def _get_n_indices(structure):
    return [i for i, species in enumerate(structure.species) if str(species) == "N"]


def _get_o_indices(structure):
    return [i for i, site in enumerate(structure) if str(site.specie) == "O"]


def _get_metal_indices(structure):
    return [i for i, site in enumerate(structure) if is_metal(site)]


def _get_rare_earth_indices(structure):
    return [i for i, site in enumerate(structure) if site.specie.is_rare_earth_metal]


def _get_alkali_alkaline_indices(structure):
    return [
        i for i, site in enumerate(structure) if site.specie.is_alkali or site.specie.is_alkaline
    ]


def _get_halogen_indices(structure):
    return [i for i, site in enumerate(structure) if site.specie.is_halogen]


def get_h_indices(structure):
    """Get the indices of all H."""
    return get_indices(structure)["h"]


def get_c_indices(structure):
    """Get the indices of all C."""
    return get_indices(structure)["c"]


def get_n_indices(structure):
    """Get the indices of all N."""
    return get_indices(structure)["n"]


def get_o_indices(structure):
    """Get the indices of all O."""
    return get_indices(structure)["o"]


def get_metal_indices(structure):
    """Get the indices of all metals."""
    return get_indices(structure)["metal"]


def get_rare_earth_indices(structure):
    """Get the indices of all rare-earth metals."""
    return get_indices(structure)["rare_earth"]


def get_halogen_indices(structure):
    """Get the indices of all halogens."""
    return get_indices(structure)["halogen"]


def get_alkali_alkaline_indices(structure):
    """Get the indices of all alkali and alkaline earth metals."""
    return get_indices(structure)["alkali_alkaline"]


def get_indices(structure: Union[Structure, IStructure]) -> dict:
    """Get all the relevant indices."""
    if isinstance(structure, Structure):
        structure = IStructure.from_sites(structure)
    return _get_indices(structure)


def _is_any_neighbor_metal(neighbors):
    return any(is_metal(neighbor.site) for neighbor in neighbors)
