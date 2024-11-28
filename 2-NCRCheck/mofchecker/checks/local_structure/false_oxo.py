# -*- coding: utf-8 -*-
"""Check if there is a site with an unexpected Oxo group.

Original idea and implementation idea contributed by Andrew Rosen.
https://github.com/kjappelbaum/mofchecker/issues/122
"""

from pymatgen.analysis.graphs import StructureGraph

from ...types import StructureIStructureType

from .base_coordination_check import BaseCoordinationCheck
from ..utils.get_indices import get_metal_indices

NO_TERMINAL_OXO = [
    "Li",
    "Na",
    "K",
    "Rb",
    "Cs",
    "Fr",
    "Be",
    "Mg",
    "Ca",
    "Sr",
    "Ba",
    "Ra",
    "Sc",
    "Y",
    "La",
    "Ac",
    "Ti",
    "Zr",
    "Hf",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Ag",
    "Zn",
    "Cd",
    "Al",
    "Ga",
    "In",
    "Tl",
]


class FalseOxoCheck(BaseCoordinationCheck):
    """Check if there is a metal with oxo group for which such a group is unexpected."""

    def __init__(self, structure: StructureIStructureType, structure_graph: StructureGraph):
        """Initialize the FalseOxoCheck check.

        Args:
            structure (StructureIStructureType): The structure to check.
            structure_graph (StructureGraph): The structure graph to use for the check.
        """
        self.structure = structure
        self.metal_indices = get_metal_indices(self.structure)
        self.structure_graph = structure_graph

    @property
    def name(self):
        """Return the name of the check."""
        return "Unexpected oxo groups"

    @property
    def description(self):
        """Return a description of the check."""
        return "Checks if there is a metal with oxo group,\
             for which such a group is unexpected."

    def _run_check(self):
        wrong_oxo = self._get_wrong_oxo()
        return len(wrong_oxo) == 0, wrong_oxo

    def _get_wrong_oxo(self):
        """Check for all metals if there are unexpected oxo group."""
        wrong_oxo = []

        for site_index in self.metal_indices:
            if str(self.structure[site_index].specie) in NO_TERMINAL_OXO:

                neighbors = self.get_connected_sites(site_index)
                for neighbor in neighbors:
                    neighbor_neighbors = self.get_connected_sites(neighbor.index)
                    if len(neighbor_neighbors) == 1 and str(neighbor.site.specie) == "O":
                        wrong_oxo.append(neighbor_neighbors[0].index)

        return wrong_oxo
