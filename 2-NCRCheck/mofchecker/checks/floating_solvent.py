# -*- coding: utf-8 -*-
"""Find connected components in the cell that do not cross PBC."""
import networkx as nx
from pymatgen.analysis.graphs import StructureGraph
from structuregraph_helpers.subgraph import get_subgraphs_as_molecules

from .check_base import AbstractIndexCheck


class FloatingSolventCheck(AbstractIndexCheck):
    """Check if there is any non-periodic connected component in the cell.

    This might be a solvent, charge compensating counter ion, etc.
    """

    def __init__(self, structure_graph: StructureGraph):
        """Create a floating solvent check instance.

        Args:
            structure_graph (StructureGraph): The structure graph to check.
        """
        self.structure_graph = structure_graph
        nx.set_node_attributes(
            self.structure_graph.graph,
            name="idx",
            values=dict(zip(range(len(structure_graph)), range(len(structure_graph)))),
        )

    @property
    def name(self):
        """Return the name of the check."""
        return "Floating atom or molecule"

    @classmethod
    def from_mofchecker(cls, mofchecker):
        """Initialize a checker from a mofchecker instance."""
        checker = cls(mofchecker.graph)
        return checker

    def _run_check(self):
        _, _, idx, _, _ = get_subgraphs_as_molecules(self.structure_graph, return_unique=False)
        return len(idx) == 0, idx

    @property
    def description(self):
        """Return a description of the check."""
        return "Checks if there is any non-periodic connected component in the cell,\
             which can be a solvent, charge compensating counter ion, etc."
