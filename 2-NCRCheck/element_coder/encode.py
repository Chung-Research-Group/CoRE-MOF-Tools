# -*- coding: utf-8 -*-
"""Numerically encode an Element."""
from typing import Iterable, List, Union

import numpy as np
from pymatgen.core import Element

from .data.coding_data import get_coding_dict

__all__ = ("encode",)


def _get_element_symbol(element):
    if isinstance(element, int):
        element = Element.from_Z(element).symbol
    if isinstance(element, Element):
        element = element.symbol
    if not isinstance(element, str):
        raise ValueError(
            f"`element` must be a string or a `pymatgen` `Element`, not {type(element)}"
        )

    return element


def encode(element: Union[Element, str, int], property: str) -> Union[float, int, np.ndarray]:
    """Numerically encode an element.

    Args:
        element (Union[Element, str, int]): Input element, e.g. "Fe",  or Element("Fe"),
            or atomic number (Z) as int.
        property (str): Property that was used for encoding, e.g. "mod_pettifor"

    Raises:
        ValueError: If element is not of type `str` or `pymatgen.core.Element`

    Returns:
        Union[float, int, np.ndarray]: Numerical encoding of element.

    Examples:
        >>> encode("Fe", "mod_pettifor")
        70
    """
    element = _get_element_symbol(element)
    if not isinstance(element, str):
        raise ValueError(
            f"`element` must be a string or a `pymatgen` `Element`, not {type(element)}"
        )

    return get_coding_dict(property)[element]


def encode_many(elements: Iterable[Union[Element, str, int]], property: str) -> List:
    """Numerically encode a collection of elements.

    Note, however, that this method brings the largest computational benefits if
    element symbols are provided as strings. If you provide, atomic numbers or `Element`
    instances, there is still an iterative conversion into strings.

    Args:
        elements (Iterable[Union[Element, str, int]]): Input element, e.g. "Fe",  or Element("Fe"),
            or atomic number (Z) as int.
        property (str): Property that was used for encoding, e.g. "mod_pettifor"

    Raises:
        ValueError: If element is not of type `str` or `pymatgen.core.Element` # noqa: DAR402

    Returns:
        List: Numerical encoding of elements.

    Examples:
        >>> encode_many(["Fe", 'H'], "mod_pettifor")
        [70, 102]
    """
    if not all(isinstance(item, str) for item in elements):
        elements = [_get_element_symbol(item) for item in elements]

    return [get_coding_dict(property)[element] for element in elements]
