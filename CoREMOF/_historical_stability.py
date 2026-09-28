"""Integrity and finite-value checks for the original bundled predictors.

This contract does not describe the later CoREMOF-COD benchmark models.
No optional scientific dependency is imported here.
"""
from __future__ import annotations

import hashlib
import math
from numbers import Real
from pathlib import Path


ASSET_SHA256 = {
    "final_model_T_few_epochs.h5": "f3e7fa57290ed5ac48d3c0af2d454c4b0853b5bf6b46a8b286505049a4897c71",
    "final_model_flag_few_epochs.h5": "0a51d3a7cd44dc598567ef9dbf1ef75ee8ef95b156acc5f5ad82c1037a38a409",
    "solvent_scaler.pkl": "021e0058335d8fd2f310ffcfcda3a7db5e9c67912887c0c5a150327b9e216a1c",
    "thermal_x_scaler.pkl": "b4c73a29307da8e3e9a066e136ae9dde9024f15c6eb4566f524bdb62fdf4810b",
    "thermal_y_scaler.pkl": "e84b7df385e4f07a23b46a1679ea70bb8aea91719692d9d0ff650cc20e85c2db",
    "water_model.pkl": "4733aa5342362f528656e325846f2af285e283733c4c7a3596f132a5be1ff3a4",
    "water_scaler.pkl": "c8062d44e46576dab6ca44fdaf822eb1afe957f64539545f23e3a7ec0e8f033f",
}
MAX_ASSET_BYTES = 4 * 1024 * 1024  # every original asset is below 2.5 MB


def copy_verified_models(source, destination):
    """Copy only trusted bytes before any pickle/Keras deserialization."""
    source, destination = Path(source), Path(destination)
    missing = [name for name in ASSET_SHA256 if not (source / name).is_file()]
    if missing:
        raise FileNotFoundError(
            "Historical stability models are missing: " + ", ".join(missing)
            + ". Supply the original asset directory with model_directory=. "
            "The later MIT benchmark weights are not replacements."
        )
    destination.mkdir(exist_ok=False)
    for name, expected in ASSET_SHA256.items():
        with (source / name).open("rb") as stream:
            payload = stream.read(MAX_ASSET_BYTES + 1)
        if len(payload) > MAX_ASSET_BYTES:
            raise ValueError(f"Historical stability asset exceeds size bound: {name}")
        if hashlib.sha256(payload).hexdigest() != expected:
            raise ValueError(f"Historical stability asset SHA-256 mismatch: {name}")
        (destination / name).write_bytes(payload)
    return destination


def validate_values(values, count, label, *, probability=False):
    """Validate without rounding, clipping, imputation, or dtype promotion."""
    if len(values) != count:
        raise ValueError(f"{label} requires exactly {count} values")
    for value in values:
        if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
            raise ValueError(f"{label} requires finite numeric values")
        if probability and not 0 <= value <= 1:
            raise ValueError(f"{label} probability is outside [0, 1]")
