# -*- coding: utf-8 -*-
"""Thresholds on the order parameters for OMS detection.

Manually tuned and hence probably not super general.

original reference:
Zimmermann, N. E. R.; Jain, A.
Local Structure Order Parameters and Site Fingerprints
for Quantification of Coordination Environment and Crystal Structure Similarity.
RSC Adv. 2020, 10 (10), 6063–6081. https://doi.org/10.1039/C9RA07755C.
"""


OP_DEF = {
    4: {
        "names": ["sq_plan", "sq", "see_saw_rect", "tet", "tri_pyr"],
        "weights": [0.2, 0.1, 0.1, 0.5, 0.5],
        "open": [0, 1, 2, 4],
    },
    5: {
        "names": ["pent_plan", "sq_pyr", "tri_bipyr"],
        "weights": [1, 0.5, 0.5],
        "open": [0, 1],
    },
    6: {"names": ["pent_pyr", "oct"], "weights": [0.3, 0.7], "open": [0]},
    7: {
        "names": ["hex_pyr", "pent_bipyr"],
        "weights": [0.7, 0.3],
        "open": [0],
    },
    8: {"names": ["hex_bipyr"], "weights": [1], "open": None},
}
