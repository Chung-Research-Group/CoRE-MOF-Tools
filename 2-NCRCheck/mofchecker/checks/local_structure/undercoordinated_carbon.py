# -*- coding: utf-8 -*-
"""Check for undercoordinated carbons."""
import numpy as np
from pymatgen.analysis.graphs import StructureGraph

from ...types import StructureIStructureType

from .base_missing_check import BaseMissingCheck
from .geometry import _maximum_angle, add_sp2_hydrogen, add_sp3_hydrogens_on_cn1
from ..utils.get_indices import get_c_indices


class UnderCoordinatedCarbonCheck(BaseMissingCheck):
    """Check for undercoordinated carbons."""

    def __init__(self, structure: StructureIStructureType, structure_graph: StructureGraph):
        """Initialize the check.

        Args:
            structure (StructureIStructureType): The structure to check]
            structure_graph (StructureGraph): The structure graph of the structure
        """
        self.structure = structure
        self.c_indices = get_c_indices(self.structure)
        self.structure_graph = structure_graph
        self._position_candidates = None

    @property
    def name(self):
        """Return the name of the check."""
        return "Undercoordinated carbon"

    @property
    def description(self):
        """Return a description of the check."""
        return "Checks, using geometric heuristics,\
             if there are any carbons that are likely undercoordinated."

    def _run_check(self):
        (
            undercoordinated_carbons,
            candidate_positions,
        ) = self._get_undercoordinated_carbons()
        assert len(undercoordinated_carbons) == len(candidate_positions), "Unexpected check error"
        return (
            len(undercoordinated_carbons) == 0,
            undercoordinated_carbons,
            candidate_positions,
        )

    def _get_undercoordinated_carbons(self, tolerance: float = 10):
        """Return a list of undercoordinated carbons and a list of candidate positions.

        Idea is that carbon should at least have three neighbors if it is not sp1.
        In sp1 case it is linear. So we can just check if there are carbons with
        non-linear coordination with less than three neighbors. An example in CoRE
        MOF would be AHOKIR. In principle this should also flag the quite common
        case of benzene rings with missing hydrogens.

        Args:
            tolerance (float): The tolerance for the angle between the neighbors of the carbon.

        Returns:
            List[int], np.typing.ArrayLike: The list of undercoordinated carbons and a list of candidate positions.
        """
        undercoordinated_carbons = []
        h_positions = []  # output must be list of lists to allow for filtering

        for site_index in self.c_indices:
            cn = self.get_cn(site_index)  # pylint:disable=invalid-name
            neighbors = self.get_connected_sites(site_index)
            if cn == 1:
                # this will fail for alkine
                undercoordinated_carbons.append(site_index)
                # make it sp3
                h_positions.append(add_sp3_hydrogens_on_cn1(self.structure[site_index], neighbors))
            if cn == 2:
                angle = _maximum_angle(
                    self.structure.get_angle(site_index, neighbors[0].index, neighbors[1].index)
                )
                if np.abs(180 - angle) > tolerance:
                    undercoordinated_carbons.append(site_index)

                    h_positions.append(add_sp2_hydrogen(self.structure[site_index], neighbors))

            # i wond't catch CN3 as this would need careful evaluation of the bond order

        return undercoordinated_carbons, h_positions
