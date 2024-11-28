# -*- coding: utf-8 -*-
"""Checks if there are atomic overlaps, based on dist < min(covr 1, covr 2)."""
import warnings

import numpy as np
from pymatgen.core import Structure
from scipy import sparse

from ...types import StructureIStructureType

from ..check_base import AbstractIndexCheck
from ..data import _get_covalent_radius


class AtomicOverlapCheck(AbstractIndexCheck):
    """Checks if there are atomic overlaps, based on dist < min(covr 1, covr 2)."""

    def __init__(self, structure: StructureIStructureType):
        """Initialize the check.

        Args:
            structure (StructureIStructureType): The structure to check.
        """
        self.structure = structure
        self.indices = None

    @property
    def name(self):
        """Return the name of the check."""
        return "Atomic overlaps"

    def _run_check(self):
        overlaps = _get_overlaps(self.structure)
        return len(overlaps) == 0, overlaps

    @property
    def description(self):
        """Return a description of the check."""
        return "True, if there are no atomic overlaps, based on dist < min(covr 1, covr 2)"


def _compute_overlap_matrix(distance_matrix: np.array, allatomtypes: list, tolerance: float = 1.0):
    """Find atomic overlap based on pairwise distance and Ccvalent radii.

    Criterion: if dist < min (covr 1, covr 2) -> overlap
        (this function is used in molsimplify)

    Args:
        distance_matrix (np.array): pairwise distance matrix
        allatomtypes (list): list of atom types
        tolerance (float): tolerance for overlap

    Returns:
        overlap_matrix (sparse matrix): overlap matrix
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("once")  # only warn once for missing radius data

        overlap_matrix = np.zeros(distance_matrix.shape)
        for i, elem_1 in enumerate(allatomtypes[:-1]):
            for j, elem_2 in enumerate(allatomtypes[i + 1 :]):
                dist = distance_matrix[i, i + j + 1]
                # check for atomic overlap:
                if dist < tolerance * min(
                    _get_covalent_radius(elem_1), _get_covalent_radius(elem_2)
                ):
                    overlap_matrix[i, i + j + 1] = 1
                    overlap_matrix[i + j + 1, i] = 1
    return sparse.csr_matrix(overlap_matrix)


def _get_overlaps(s: Structure) -> list:
    """Find overlapping atoms in a structure."""
    distance_matrix = s.distance_matrix
    atomtypes = [str(species) for species in s.species]
    overlap_matrix = _compute_overlap_matrix(distance_matrix, atomtypes)
    overlap_atoms = []
    for atom in set(sparse.find(overlap_matrix)[0]):
        overlap_atoms.append(atom.item())
    return overlap_atoms
