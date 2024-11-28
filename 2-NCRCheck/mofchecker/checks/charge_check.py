# -*- coding: utf-8 -*-
"""Check that the charges of the structure are reasonable."""
import warnings
from tempfile import NamedTemporaryFile

import numpy as np

from .check_base import AbstractCheck
from ..types import StructureIStructureType


class ChargeCheck(AbstractCheck):
    """Check that the charges of the structure are reasonable."""

    def __init__(self, structure: StructureIStructureType):
        """Create a charge check instance.

        Args:
            structure (StructureIStructureType): The structure to check.
        """
        self.structure = structure
        self.threshold = 3

    @property
    def name(self) -> str:
        """Return the name of the check."""
        return "High charges"

    @property
    def description(self) -> str:
        """Return a description of the check."""
        return f"Check that the charges of the structure are reasonable\
             (abs(charge) not higher than {self.threshold})."

    def _run_check(self):
        try:
            from pyeqeq.main import run_on_cif

            with NamedTemporaryFile("w", suffix=".cif") as file:
                self.structure.to(fmt="cif", filename=file.name)
                charges = run_on_cif(file.name, verbose=False)
                has_high_charges = np.sum(np.abs(charges) > self.threshold)

            return not has_high_charges
        except ImportError:
            warnings.warn("Install the eqeq extra to run the charge check")
            return None
