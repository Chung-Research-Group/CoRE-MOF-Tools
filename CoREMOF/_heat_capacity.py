"""Standard-library preflight for the bundled heat-capacity predictor."""

import hashlib
import importlib.metadata
import json
import math
from numbers import Real
from pathlib import Path


def normalize_temperatures(temperatures):
    """Use exact integer kelvin directories, never rounded temperatures."""
    result = []
    for temperature in temperatures:
        if isinstance(temperature, bool) or not isinstance(temperature, Real):
            raise ValueError("Temperatures must be positive whole numbers in kelvin")
        value = float(temperature)
        if not math.isfinite(value) or value <= 0 or not value.is_integer():
            raise ValueError("Temperatures must be positive whole numbers in kelvin")
        result.append(int(value))
    if not result:
        raise ValueError("At least one prediction temperature is required")
    if len(set(result)) != len(result):
        raise ValueError("Prediction temperatures must be distinct")
    return result


def validate_ensemble(directory, temperatures):
    """Verify the complete known ensemble before unpickling any model.

    These hashes identify the supplied repository assets. They do not establish
    that those assets are byte-identical to the original paper's model files.
    """
    root = Path(directory)
    manifest = json.loads(
        (Path(__file__).parent / "models" / "heat_capacity_assets.json").read_text()
    )
    for temperature in temperatures:
        expected = manifest["temperatures"].get(str(temperature))
        if expected is None:
            raise ValueError(
                f"No bundled heat-capacity model at {temperature} K; use 300, 350 or 400 K"
            )
        folder = root / str(temperature)
        if not folder.is_dir():
            raise FileNotFoundError(
                f"Heat-capacity ensemble models are missing for {temperature} K under {root}"
            )
        if folder.is_symlink() or any(p.is_symlink() for p in folder.parents):
            raise ValueError("Heat-capacity model directories must not be symlinks")
        names = {path.name for path in folder.iterdir()}
        if names != set(expected):
            raise ValueError(
                f"Expected the complete {len(expected)}-model ensemble at {temperature} K; "
                f"missing={sorted(set(expected) - names)}, unexpected={sorted(names - set(expected))}"
            )
        for name, record in expected.items():
            path = folder / name
            if path.is_symlink() or not path.is_file():
                raise ValueError(f"Heat-capacity model is not a regular file: {path}")
            digest = hashlib.sha256()
            size = 0
            with path.open("rb") as stream:
                for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                    size += len(chunk)
                    digest.update(chunk)
            if size != record["bytes"] or digest.hexdigest() != record["sha256"]:
                raise ValueError(f"Heat-capacity model checksum mismatch: {path}")
    for package, key in (("scikit-learn", "scikit_learn_version"), ("xgboost", "xgboost_version")):
        try:
            version = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError as exc:
            raise ImportError("Install CoREMOF-tools[heat-capacity] in a separate environment") from exc
        if version != manifest[key]:
            raise RuntimeError(
                f"The supplied heat-capacity models require {package} "
                f"{manifest[key]}, found {version}. "
                "Install CoREMOF-tools[heat-capacity] in a separate environment."
            )
