# -*- coding: utf-8 -*-
"""Flagging overcoordinated hydrogens."""
from pymatgen.analysis.graphs import StructureGraph

from ...types import StructureIStructureType

from .base_coordination_check import BaseCoordinationCheck
from ..data import _get_vdw_radius
from ..utils.get_indices import get_h_indices


class OverCoordinatedHydrogenCheck(BaseCoordinationCheck):
    """Flagging overcoordinated hydrogens."""

    def __init__(self, structure: StructureIStructureType, structure_graph: StructureGraph):
        """Initialize the OverCoordinatedHydrogenCheck check.

        Args:
            structure (StructureIStructureType): The structure to check.
            structure_graph (StructureGraph): The structure graph to use for the check.
        """
        self.structure = structure
        self.h_indices = get_h_indices(self.structure)
        self.structure_graph = structure_graph

    @property
    def name(self):
        """Return the name of the check."""
        return "Overcoordinated hydrogen"

    @property
    def description(self):
        """Return a description of the check."""
        return "Checks, using geometric heuristics,\
             if there are any hydrogen that are likely overcoordinated (i.e., CN>1)."

    def _run_check(self):
        overcoordinated_hydrogens = self._get_overcoordinated_hydrogens()
        return len(overcoordinated_hydrogens) == 0, overcoordinated_hydrogens

    def _get_overcoordinated_hydrogens(self):
        """Check for all H if CN>1, ignore metal bonds."""
        overcoordinated_hydrogens = []

        for site_index in self.h_indices:
            neighbors = self.structure.get_neighbors(
                self.structure[site_index], _get_vdw_radius("H")
            )  # pylint:disable=invalid-name
            if len(neighbors) > 1:
                overcoordinated_hydrogens.append(site_index)
        return overcoordinated_hydrogens
