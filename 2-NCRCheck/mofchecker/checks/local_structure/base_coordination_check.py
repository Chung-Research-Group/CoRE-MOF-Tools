# -*- coding: utf-8 -*-
"""Base class for checks on coordination numbers/environments."""
import abc

from pymatgen.analysis.graphs import StructureGraph
from structuregraph_helpers.analysis import get_cn

from ...types import StructureIStructureType

from ..check_base import AbstractIndexCheck


class BaseCoordinationCheck(AbstractIndexCheck):
    """Base class for checks on coordination numbers/environments."""

    @abc.abstractmethod
    def __init__(self, structure: StructureIStructureType, structure_graph: StructureGraph):
        """Intialize the check.

        Args:
            structure (StructureIStructureType): The structure to check
            structure_graph (StructureGraph): The structure graph of the structure
        """
        self.structure = structure
        self.structure_graph = structure_graph

    def get_cn(self, index):
        """Get coordination number of index."""
        return get_cn(self.structure_graph, index)

    def get_connected_sites(self, index):
        """Get sites connected to index."""
        return self.structure_graph.get_connected_sites(index)

    @classmethod
    def from_mofchecker(cls, mofchecker):
        """Create a checker instance from a mofchecker instance."""
        checker = cls(mofchecker.structure, mofchecker.graph)
        checker.get_cn = mofchecker.get_cn
        checker.get_connected_sites = mofchecker.get_connected_sites
        return checker
