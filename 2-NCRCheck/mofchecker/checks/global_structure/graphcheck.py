"""Checks operating on the structure graph."""
from pymatgen.analysis.dimensionality import get_dimensionality_larsen

from ..check_base import AbstractCheck


class BaseStructureGraphCheck(AbstractCheck):
    """Base class for checks operating on the structure graph."""

    def __init__(self, structure_graph):
        """Initialize the check.

        Args:
            structure_graph (StructureGraph): The structure graph of the structure
        """
        self.structure_graph = structure_graph

    @classmethod
    def from_mofchecker(cls, mofchecker):
        """Initialize a checker from a mofchecker instance."""
        checker = cls(mofchecker.graph)
        return checker


class IsThreeDimensional(BaseStructureGraphCheck):
    """Check if the structure is 3D."""

    def _run_check(self):
        return bool(get_dimensionality_larsen(self.structure_graph) == 3)

    @property
    def name(self):
        """Return the name of the check."""
        return "3D connected structure graph."

    @property
    def description(self):
        """Return a description of the check."""
        return "Check if the structure graph is 3D connected."
