# -*- coding: utf-8 -*-
"""Check if there are any metals that are sterically exposed."""
from pymatgen.analysis.graphs import StructureGraph

from .base_coordination_check import BaseCoordinationCheck
from ..utils.geometry import get_open_angle
from ..utils.get_indices import (
    get_alkali_alkaline_indices,
    get_metal_indices,
    get_rare_earth_indices,
)
from ...types import StructureIStructureType


class GeometricallyExposedMetal(BaseCoordinationCheck):
    """Check if there are any metals that are likely geometrically exposed.

    We consider alkali/alkaline earth or rare earth metals.

    That is, which form a small cone angle with their binding partners.
    """

    def __init__(
        self,
        structure: StructureIStructureType,
        structure_graph: StructureGraph,
        tight: bool = True,
    ):
        """Construct a GeometricallyExposedMetal check.

        Args:
            structure (StructureIStructureType): structure to check
            structure_graph (StructureGraph): structure graph of the structure
            tight (bool): whether to use a tight metal set of test all metals
        """
        self.structure = structure
        if not tight:
            self.relevant_metals = get_alkali_alkaline_indices(structure) + get_rare_earth_indices(
                structure
            )
        else:
            self.relevant_metals = get_metal_indices(structure)
        self.structure_graph = structure_graph
        self.threshold = 150

    @property
    def name(self):
        """Return the name of the check."""
        return "Geometrically exposed metal."

    @property
    def description(self):
        """Return a description of the check."""
        return "Check if there are any metals (alkali/alkaline earth or rare earth) that are likely \
            geometrically exposed, i.e. which form a small cone angle \
                with their binding partners"

    def _run_check(self):
        exposed_metals = self._get_exposed_metals()
        return (
            len(exposed_metals) == 0,
            exposed_metals,
        )

    def _get_exposed_metals(self):
        """Check for all geometrically exposed metals."""
        geometrically_exposed_metals = []

        for site_index in self.relevant_metals:
            angle = get_open_angle(self.structure_graph, site_index)
            # print(angle, self.get_cn(site_index), site_index)
            if angle > self.threshold:
                if self.get_cn(site_index) < 6:
                    geometrically_exposed_metals.append(site_index)

        return geometrically_exposed_metals
