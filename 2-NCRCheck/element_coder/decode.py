# -*- coding: utf-8 -*-
"""Decode an elemental encoding."""
from typing import Callable, Union

import numpy as np
from loguru import logger
from numpy.typing import ArrayLike
from scipy.spatial.distance import cdist

from .data.coding_data import get_coding_dict


def _get_coding_data(coding_dict):
    try:
        coding_data = np.fromiter(coding_dict.values(), dtype=float)
    except ValueError:
        coding_data = np.vstack(coding_dict.values())
    coding_shape = coding_data.shape
    if len(coding_shape) == 1:
        coding_data = coding_data.reshape(-1, 1)
    return coding_data


def decode(
    encoding: Union[int, float, np.ndarray, list, tuple],
    property: str,
    metric: Union[str, Callable] = "euclidean",
) -> str:
    """Decode an elemental encoding.

    Args:
        encoding (Union[int, float, np.ndarray]): Numerical encoding of an element
        property (str): Property that was used for encoding, e.g. "mod_pettifor"
        metric (Union[str, callable] optional): Metric to measure distance between (noisy) input encoding
            and tabulated encodings.
            If a string, the distance function can be 'braycurtis', 'canberra', 'chebyshev',
            'cityblock', 'correlation', 'cosine', 'dice', 'euclidean', 'hamming', 'jaccard',
            'jensenshannon', 'kulsinski', 'kulczynski1', 'mahalanobis', 'matching', 'minkowski',
            'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath',
            'sqeuclidean', 'yule'. Defaults to "euclidean".

    Returns:
        str: Element symbol

    Examples:
        >>> decode(80, "mod_pettifor")
        'Tl'
    """
    if isinstance(encoding, (int, float)):
        encoding = np.array([encoding])
    elif isinstance(encoding, (list, tuple)):
        encoding = np.array(encoding)

    coding_dict = get_coding_dict(property)
    coding_data = _get_coding_data(coding_dict)
    distance = cdist(coding_data, encoding.reshape(1, -1), metric=metric)
    matching_index = distance.argmin()
    logger.debug(
        f"Matching distance: {distance[matching_index]}, mean distance: {distance.mean()}, std: {distance.std()}"
    )

    return list(coding_dict.keys())[matching_index]


def decode_many(
    encoding: ArrayLike,
    property: str,
    metric: Union[str, Callable] = "euclidean",
) -> np.array:
    """Decode a collecgtio of elemental encodings.

    Args:
        encoding (ArrayLike): Numerical encodings of elements
        property (str): Property that was used for encoding, e.g. "mod_pettifor"
        metric (Union[str, callable] optional): Metric to measure distance between (noisy) input encoding
            and tabulated encodings.
            If a string, the distance function can be 'braycurtis', 'canberra', 'chebyshev',
            'cityblock', 'correlation', 'cosine', 'dice', 'euclidean', 'hamming', 'jaccard',
            'jensenshannon', 'kulsinski', 'kulczynski1', 'mahalanobis', 'matching', 'minkowski',
            'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath',
            'sqeuclidean', 'yule'. Defaults to "euclidean".

    Returns:
        np.array: Numpy array of element symbols

    Examples:
        >>> decode([80, 90], "mod_pettifor")
        ['Tl', 'Sb']
    """
    encoding = np.array(encoding)

    coding_dict = get_coding_dict(property)
    coding_data = _get_coding_data(coding_dict)
    distance = cdist(coding_data, encoding.reshape(-1, coding_data.shape[1]), metric=metric)
    matching_index = distance.argmin(axis=0)
    logger.debug(
        f"Matching distance: {distance[matching_index]}, mean distance: {distance.mean()}, std: {distance.std()}"
    )

    keys = list(coding_dict.keys())

    return [keys[i] for i in matching_index]
