# -*- coding: utf-8 -*-
"""Check if there are carbons with more neighbors than expected."""
from pymatgen.analysis.graphs import StructureGraph

from ...types import StructureIStructureType

from .base_coordination_check import BaseCoordinationCheck
from ..utils.get_indices import _is_any_neighbor_metal, get_c_indices


class OverCoordinatedCarbonCheck(BaseCoordinationCheck):
    """Check if there are carbons with more neighbors than expected."""

    def __init__(self, structure: StructureIStructureType, structure_graph: StructureGraph):
        """Initialize the OverCoordinatedCarbonCheck check.

        Args:
            structure (StructureIStructureType): The structure to check.
            structure_graph (StructureGraph): The structure graph to use for the check.
        """
        self.structure = structure
        self.c_indices = get_c_indices(self.structure)
        self.structure_graph = structure_graph

    @property
    def name(self):
        """Return the name of the check."""
        return "Overcoordinated carbons"

    @property
    def description(self):
        """Return a description of the check."""
        return "Checks, using geometric heuristics,\
             if there are any carbons that are likely overcoordinated (i.e., CN>4)."

    def _run_check(self):
        overcoordinated_carbons = self._get_overcoordinated_carbons()
        return len(overcoordinated_carbons) == 0, overcoordinated_carbons

    def _get_overcoordinated_carbons(self):
        """Check for all C if CN>4, ignore metal bonds."""
        overcoordinated_carbons = []

        for site_index in self.c_indices:
            cn = self.get_cn(site_index)  # pylint:disable=invalid-name
            if cn > 4 and not _is_any_neighbor_metal(self.get_connected_sites(site_index)):
                overcoordinated_carbons.append(site_index)

        return overcoordinated_carbons
