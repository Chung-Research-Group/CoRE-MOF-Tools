# -*- coding: utf-8 -*-
"""Check if there are any lanthanides/actinides that are likely undercoordinated (i.e., CN<4)."""
from pymatgen.analysis.graphs import StructureGraph

from .base_coordination_check import BaseCoordinationCheck
from ..utils.get_indices import get_rare_earth_indices
from ...types import StructureIStructureType


class UnderCoordinatedRareEarthCheck(BaseCoordinationCheck):
    """Check if there are any lanthanides/actinides that are likely undercoordinated (i.e., CN<4)."""

    def __init__(self, structure: StructureIStructureType, structure_graph: StructureGraph):
        """Initialize the UnderCoordinatedRareEarthCheck check.

        Args:
            structure (StructureIStructureType): The structure to check.
            structure_graph (StructureGraph): The structure graph to use for the check.
        """
        self.structure = structure
        self.rare_earth_indices = get_rare_earth_indices(structure)
        self.structure_graph = structure_graph

    @property
    def name(self):
        """Return the name of the check."""
        return "Undercoordinated rare earth metal"

    @property
    def description(self):
        """Return a description of the check."""
        return "Check if there are any lanthanides/actinides\
            that are likely undercoordinated (i.e., CN<4)."

    def _run_check(self):
        undercoordinated_rare_earth_metals = self._get_undercoordinated_rare_earth_metals()
        return (
            len(undercoordinated_rare_earth_metals) == 0,
            undercoordinated_rare_earth_metals,
        )

    def _get_undercoordinated_rare_earth_metals(self):
        """Check for all rare earth metals if CN < 4."""
        undercoordinated_rare_earth_metals = []

        for site_index in self.rare_earth_indices:
            cn = self.get_cn(site_index)  # pylint:disable=invalid-name
            if cn < 4:
                undercoordinated_rare_earth_metals.append(site_index)

        return undercoordinated_rare_earth_metals
