# -*- coding: utf-8 -*-
"""Running basic pore analysis with zeo++."""
import os
import subprocess
import warnings
from tempfile import TemporaryDirectory
from typing import Union

import numpy as np
from pymatgen.core import Structure

from .check_base import AbstractCheck
from .utils import is_tool
from ..types import StructureIStructureType

ZEOPP_BASE_COMMAND = ["network", "-ha", "-res"]
NO_ZEOPP_WARNING = "Did not find the zeo++ network binary in the path. \
            Can not run pore analysis."

__all__ = ["check_if_porous"]


def run_zeopp(structure: Structure) -> dict:
    """Run zeopp with network -ha -res to find the pore diameters.

    Args:
        structure (Structure): pymatgen Structure object

    Returns:
        dict: pore analysis results
    """
    if is_tool("network"):
        with TemporaryDirectory() as tempdir:
            structure_path = os.path.join(tempdir, "structure.cif")
            result_path = os.path.join(tempdir, "result.res")
            structure.to(fmt="cif", filename=structure_path)
            cmd = ZEOPP_BASE_COMMAND + [str(result_path), str(structure_path)]
            _ = subprocess.run(
                cmd,
                universal_newlines=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
            )

            with open(result_path, "r") as handle:
                results = handle.read()

            zeopp_results = _parse_zeopp(results)

            return zeopp_results
    else:
        warnings.warn(NO_ZEOPP_WARNING)
        return {
            "lis": np.nan,  # largest included sphere
            "lifs": np.nan,  # largest free sphere
            "lifsp": np.nan,  # largest included sphere along free sphere path
        }


def _parse_zeopp(filecontent: str) -> dict:
    """Parse the results line of a network call to zeopp.

    Args:
        filecontent (str): results file

    Returns:
        dict: largest included sphere, largest free sphere,
            largest included sphera along free sphere path
    """
    first_line = filecontent.split("\n")[0]
    parts = first_line.split()

    results = {
        "lis": float(parts[1]),  # largest included sphere
        "lifs": float(parts[2]),  # largest free sphere
        "lifsp": float(parts[3]),  # largest included sphere along free sphere path
    }

    return results


def check_if_porous(structure: Structure, threshold: float = 2.4) -> Union[bool, None]:
    """Run zeo++ to check if structure is porous.

    We consider the CoRE-MOF definition (PLD > 2.4, https://pubs.acs.org/doi/10.1021/acs.jced.9b00835)

    Args:
        structure (Structure): MOF structure to check
        threshold (float): Threshold on the sphere diameter in Angstrom.
            Defaults to 2.4.

    Returns:
        bool: True if porous.
    """
    if is_tool("network"):
        zeopp_results = run_zeopp(structure)
        if zeopp_results["lifs"] >= threshold:
            return True
        return False

    warnings.warn(NO_ZEOPP_WARNING)
    return None


class PorosityCheck(AbstractCheck):
    """Use zeo++ to check if the structure is porous."""

    def __init__(self, structure: StructureIStructureType):
        """Create a porosity check instance."""
        self.structure = structure
        self.threshold = 2.4

    @property
    def name(self):
        """Return the name of the check."""
        return "Porosity"

    @property
    def description(self):
        """Return a description of the check."""
        return f"Check if the pore limiting diameter is greater than {self.threshold}."

    def _run_check(self):
        return check_if_porous(self.structure)
