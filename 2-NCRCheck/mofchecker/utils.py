# -*- coding: utf-8 -*-
"""Helper functions for the MOFChecker."""
import functools
import json
import pickle
import warnings
from types import FunctionType

import pymatgen
from backports.cached_property import cached_property

from .types import PathType


def deprecated(func: FunctionType) -> FunctionType:
    """Mark function as deprecated using a decorator.

    It will result in a warning being emitted
    when the function is used.

    Args:
        func (FunctionType): function to be decorated

    Returns:
        FunctionType: decorated function
    """

    @functools.wraps(func)
    def new_func(*args, **kwargs):
        warnings.simplefilter("always", DeprecationWarning)  # turn off filter
        warnings.warn(
            "Call to deprecated function {}.".format(func.__name__),
            category=DeprecationWarning,
            stacklevel=2,
        )
        warnings.simplefilter("default", DeprecationWarning)  # reset filter
        return func(*args, **kwargs)

    return new_func


def read_pickle(file: PathType):
    """Read a pickle file."""
    with open(file, "rb") as handle:
        return pickle.load(handle)


def _check_metal_coordination(site, coordination_number: int) -> bool:
    # Lanthanides like to have many neighbors
    # Low coordinatio number is usually only
    # possible with really bulky ligands
    # for this reason, a Lanthanide with low
    # coordination number, e.g., <= 4 can be considered "interesting"
    if (
        (site.specie.is_lanthanoid)
        or (site.specie.is_actinoid)
        or (site.specie.symbol in ("Mo", "Cr", "Hf", "Mb"))
    ) and coordination_number <= 4:
        return True
    return False


def print_dict(dictionary: dict) -> None:
    """Print a dictionary to stdout line by line."""
    for k, v in dictionary.items():
        print(k, v)  # noqa:T201


def _check_if_ordered(structure):
    if not structure.is_ordered:
        raise NotImplementedError(
            "Support of unordered structures with partial occupancies \
                is not implemented (yet)."
        )
    for site in structure:
        if site.specie.atomic_radius is None:
            raise NotImplementedError(
                f"Pymatgen currently does not support this {str(site.specie)} element"
            )


class IStructure(pymatgen.core.structure.IStructure):
    """pymatgen IStructure with faster equality comparison.

    This dramatically speeds up lookups in the LRU cache when an object
    with the same __hash__ is already in the cache.
    """

    __hash__ = pymatgen.core.structure.IStructure.__hash__

    def __eq__(self, other):
        """Use specific, yet performant hash for equality comparison."""
        return self._dict_hash == other._dict_hash

    @cached_property
    def _dict_hash(self):
        """Specific, yet performant hash."""
        return hash(json.dumps(self.as_dict(), sort_keys=True))
