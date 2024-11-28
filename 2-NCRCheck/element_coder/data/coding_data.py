# -*- coding: utf-8 -*-
"""We represent each possible elemental encoding as a `sciris` dict. This allows for more convenient bi-directional lookup.
"""

import json
import os
from collections import OrderedDict
from difflib import get_close_matches
from functools import lru_cache

import numpy as np
from loguru import logger

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_CODING_DATA_DIR = os.path.join(_THIS_DIR, "raw")

CODING_FILES = {
    "x": "x.json",
    "electronegativity": "x.json",
    "atomic": "atomic.json",
    "atomic_radius_calculated": "atomic_radius_calculated.json",
    "boiling_point": "boiling_point.json",
    "brinell_hardness": "brinell_hardness.json",
    "bulk_modulus": "bulk_modulus.json",
    "cgcnn": "cgcnn.json",
    "coefficient_of_linear_thermal_expansion": "coefficient_of_linear_thermal_expansion.json",
    "critical_temperature": "critical_temperature.json",
    "density_of_solid": "density_of_solid.json",
    "electrical_resistivity": "electrical_resistivity.json",
    "electron_affinity": "electron_affinity.json",
    "elemnet": "elemnet.json",
    "jarvis": "jarvis.json",
    "jarvis_sc": "jarvis_sc.json",
    "magpie": "magpie.json",
    "magpie_sc": "magpie_sc.json",
    "mat2vec": "mat2vec.json",
    "matscholar": "matscholar.json",
    "megnet16": "megnet16.json",
    "melting_point": "melting_point.json",
    "mendeleev": "mendeleev.json",
    "mineral_hardness": "mineral_hardness.json",
    "mod_petti": "mod_petti.json",
    "molar_volume": "molar_volume.json",
    "oliynyk": "oliynyk.json",
    "oliynyk_sc": "oliynyk_sc.json",
    "petti": "petti.json",
    "poissons_ratio": "poissons_ratio.json",
    "reflectivity": "reflectivity.json",
    "refractive_index": "refractive_index.json",
    "superconduction_temperature": "superconduction_temperature.json",
    "thermal_conductivity": "thermal_conductivity.json",
    "van_der_waals_radius": "van_der_waals_radius.json",
    "velocity_of_sound": "velocity_of_sound.json",
    "vickers_hardness": "vickers_hardness.json",
    "youngs_modulus": "youngs_modulus.json",
}

_PROPERTY_KEYS = set(list(CODING_FILES.keys()))

__all__ = ("get_coding_dict",)


def _load_json(file):
    with open(file, "r") as f:
        return json.load(f)


def _load_coding_data(property_key: str) -> OrderedDict:
    """Load the coding data for a given property key."""
    file = os.path.join(_CODING_DATA_DIR, CODING_FILES[property_key])
    return OrderedDict(_load_json(file))


@lru_cache()
def get_coding_dict(key: str) -> dict:
    """Get the coding dict for the given key.
    If no exact match is found, it performs a fuzzy search to find a close match.

    Args:
        key (str): property name

    Returns:
        dict: property dictionary
    """
    input_key = key
    key = key.lower()
    if not key in _PROPERTY_KEYS:
        key = get_close_matches(key, _PROPERTY_KEYS, n=1)[0]
        logger.warning(
            f"No matching key found for {input_key}. Falling back to closest match {key}."
        )
        cd = _load_coding_data(key)
    cd = _load_coding_data(key)

    if len(cd) != len(np.unique(cd.values())):
        logger.warning(
            f"This coding is not unique for certain elements. This will cause problems when decoding."
        )

    return cd
