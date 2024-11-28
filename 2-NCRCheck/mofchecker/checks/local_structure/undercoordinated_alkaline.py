# -*- coding: utf-8 -*-
"""Check if there are any alkali/alkaline earth metals that are likely undercoordinated (i.e., CN<4)."""
from pymatgen.analysis.graphs import StructureGraph

from .base_coordination_check import BaseCoordinationCheck
from ..utils.get_indices import get_alkali_alkaline_indices
from ...types import StructureIStructureType


class UnderCoordinatedAlkaliAlkaline(BaseCoordinationCheck):
    """Check if there are any alkali/alkaline earth metals that are likely undercoordinated (i.e., CN<4)."""

    def __init__(self, structure: StructureIStructureType, structure_graph: StructureGraph):
        """Initialize the UnderCoordinatedAlkaliAlkaline check.

        Args:
            structure (StructureIStructureType): The structure to check.
            structure_graph (StructureGraph): The structure graph to use for the check.
        """
        self.structure = structure
        self.alkali_alkaline_indices = get_alkali_alkaline_indices(structure)
        self.structure_graph = structure_graph

    @property
    def name(self):
        """Return the name of the check."""
        return "Undercoordinated alkali/alkaline earth metal"

    @property
    def description(self):
        """Return a description of the check."""
        return "Check if there are any alkali/alkaline earth metals\
            that are likely undercoordinated (i.e., CN<4)."

    def _run_check(self):
        undercoordinated_alkali_alkaline_earth_metals = self._get_undercoordinated_alkali_alkaline()
        return (
            len(undercoordinated_alkali_alkaline_earth_metals) == 0,
            undercoordinated_alkali_alkaline_earth_metals,
        )

    def _get_undercoordinated_alkali_alkaline(self):
        """Check for all alkali/alkaline earth metals of CN < 4."""
        undercoordinated_alkali_alkaline_earth_metals = []

        for site_index in self.alkali_alkaline_indices:
            cn = self.get_cn(site_index)  # pylint:disable=invalid-name
            if cn < 4:
                undercoordinated_alkali_alkaline_earth_metals.append(site_index)

        return undercoordinated_alkali_alkaline_earth_metals
