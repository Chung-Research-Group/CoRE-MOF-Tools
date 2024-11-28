# -*- coding: utf-8 -*-
"""Utilities for `element_coder`."""
import sys
from typing import Iterable, List, Union

from loguru import logger
from pymatgen.core import Element

from element_coder.encode import encode

__all__ = ["enable_logging"]


def get_range(elements: Iterable[Union[str, Element, int]], property: str = "Z"):
    """Get the range of a property for a list of elements.

    Args:
        elements (Iterable[Union[str, Element, int]]): List of elements.
        property (str): Property that was used for encoding, e.g. "mod_pettifor"

    Returns:
        Tuple[float, int, np.ndarray]: Numerical encoding of element.

    Examples:
        >>> get_range(["Fe", "Tl", "Pb", "U", "Zn"], "mod_pettifor")
        (70, 74)
    """
    encodings = [encode(element, property) for element in elements]
    return min(encodings), max(encodings)


def enable_logging() -> List[int]:
    """Set up the element_coder logging with sane defaults."""
    logger.enable("element_coder")

    config = dict(
        handlers=[
            dict(
                sink=sys.stderr,
                format="<green>{time:YYYY-MM-DD HH:mm:ss.SSS Z UTC}</>"
                " <red>|</> <lvl>{level}</> <red>|</> <cyan>{name}:{function}:{line}</>"
                " <red>|</> <lvl>{message}</>",
                level="INFO",
            ),
            dict(
                sink=sys.stderr,
                format="<red>{time:YYYY-MM-DD HH:mm:ss.SSS Z UTC} | {level} | {name}:{function}:{line} | {message}</>",
                level="WARNING",
            ),
        ]
    )
    return logger.configure(**config)
