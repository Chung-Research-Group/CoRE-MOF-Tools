"""Helpers for analysing the geometry of a structure."""
import numpy as np
from element_coder.encode import encode_many
from libconeangle import cone_angle
from numpy.linalg import matrix_rank
from pymatgen.analysis.graphs import StructureGraph


def are_coplanar(coords: np.typing.ArrayLike) -> bool:
    """Check if the given coordinates are coplanar.

    Args:
        coords (np.typing.ArrayLike): The coordinates to check.

    Returns:
        bool: True if the coordinates are coplanar, False otherwise.
    """
    coords = np.unique(np.array(coords), axis=0)
    coords -= coords.mean(axis=0)
    if matrix_rank(coords, tol=0.1) <= 2:
        return True
    return False


def _get_coords_and_elements_of_neighbors(graph, index):
    neighbors = graph.get_connected_sites(index)
    coords = []
    species = []
    coords.append(graph.structure.cart_coords[index])
    species.append(str(graph.structure[index].specie))
    for neighbor in neighbors:
        coords.append(neighbor.site.coords)
        species.append(str(neighbor.site.specie))
    return coords, species


def get_open_angle(graph: StructureGraph, index: int) -> float:
    """Get 360 - cone angle of the site with the given index.

    Args:
        graph (StructureGraph): The StructureGraph to analyse.
        index (int): The index of the site to analyse.

    Returns:
        float: The open angle of the site.
    """
    coords, species = _get_coords_and_elements_of_neighbors(graph, index)
    encodings = encode_many(species, "van_der_waals_radius")
    # encodings = [0.5] * len(encodings)
    try:
        angle, _, _ = cone_angle(coords, encodings, 0)
        return 360 - angle
    except ValueError:
        coords = np.unique(np.array(coords), axis=0)
        coords -= coords.mean(axis=0)
        if matrix_rank(coords) <= 2:
            return 180
        return np.nan


def has_open_angle(graph: StructureGraph, index: int, threshold: float = 80) -> bool:
    """Check if the site with the given index has an open angle.

    Nans are treated as False.

    Args:
        graph (StructureGraph): The StructureGraph to analyse.
        index (int): The index of the site to analyse.
        threshold (float): The threshold for the open angle. Defaults to 80.

    Returns:
        bool: True if the site has an open angle, False otherwise.
    """
    angle = get_open_angle(graph, index)

    if angle > threshold:
        return True
    return False
