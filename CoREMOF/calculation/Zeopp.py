"""Geometric-property calculations backed by the Zeo++ ``network`` binary.

Install Zeo++ independently, for example with
``conda install -c conda-forge zeopp-lsmo``.  The executable can be overridden
with the ``COREMOF_NETWORK_EXECUTABLE`` environment variable.
"""

from __future__ import annotations

import math
from numbers import Real
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from typing import Iterable

from .. import _release_zeopp_n2_he_protocol as _probe_protocol
from .. import _release_zeopp_framework_protocol as _framework_protocol


def _run_network(structure: str | os.PathLike, arguments: Iterable[object], prefix: str) -> str:
    """Run Zeo++ safely and return the generated output text.

    A unique output file is used for every invocation.  This is important for
    high-throughput workflows where several structures may be analysed in the
    same working directory at once.
    """

    structure_path = Path(structure)
    if not structure_path.is_file():
        raise FileNotFoundError(f"CIF file does not exist: {structure_path}")

    executable = os.environ.get("COREMOF_NETWORK_EXECUTABLE", "network")
    if not shutil.which(executable):
        raise FileNotFoundError(
            f"Zeo++ executable '{executable}' was not found. Install zeopp-lsmo "
            "or set COREMOF_NETWORK_EXECUTABLE."
        )

    prefix_path = Path(prefix)
    output_dir = prefix_path.parent if str(prefix_path.parent) != "." else Path.cwd()
    if not output_dir.is_dir():
        raise FileNotFoundError(f"Temporary-output directory does not exist: {output_dir}")

    arguments = [str(value) for value in arguments]
    # -strinfo accepts only the CIF and writes CIF_BASENAME.strinfo. Other
    # operations accept an explicit output path. Keep both the input and all
    # auxiliary files in this call's disposable directory.
    with tempfile.TemporaryDirectory(prefix=f"{prefix_path.name}_", dir=output_dir) as temporary:
        private = Path(temporary).resolve()
        isolated_cif = private / "input.cif"
        isolated_cif.write_bytes(structure_path.read_bytes())
        framework = "-strinfo" in arguments
        output_path = isolated_cif.with_suffix(".strinfo") if framework else private / "output.txt"
        command = [str(Path(shutil.which(executable)).resolve()), *arguments]
        if not framework:
            command.append(str(output_path))
        command.append(str(isolated_cif))
        completed = subprocess.run(command, cwd=private, capture_output=True, text=True, check=False)
        if completed.returncode != 0:
            detail = completed.stderr.strip() or completed.stdout.strip() or "no diagnostic output"
            raise RuntimeError(
                f"Zeo++ failed with exit code {completed.returncode}: {detail}"
            )
        if not output_path.is_file():
            raise RuntimeError("Zeo++ completed without creating its output file")
        return output_path.read_text(encoding="utf-8")


def _arguments(high_accuracy: bool, *values: object) -> list[object]:
    if type(high_accuracy) is not bool:
        raise ValueError("high_accuracy must be a Boolean")
    return (["-ha"] if high_accuracy else []) + list(values)


def _validate_probe_options(*radii: Real, num_samples=5000) -> None:
    for radius in radii:
        if isinstance(radius, bool) or not isinstance(radius, Real):
            raise ValueError("Probe and channel radii must be finite nonnegative numbers")
        if not math.isfinite(radius) or radius < 0:
            raise ValueError("Probe and channel radii must be finite nonnegative numbers")
    if type(num_samples) is not int or num_samples <= 0:
        raise ValueError("num_samples must be a positive integer, not a Boolean")


def _finite_nonnegative(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number < 0:
        raise ValueError("Zeo++ returned a non-finite or negative geometric property")
    return number


def ChanDim(structure, probe_radius=0, high_accuracy=True, prefix="tmp_chan"):
    """Return the maximum accessible-channel dimension, or zero for no channel."""

    _validate_probe_options(probe_radius)
    text = _run_network(
        structure, _arguments(high_accuracy, "-chan", probe_radius), prefix
    )
    first_line = text.splitlines()[0] if text.splitlines() else ""
    try:
        # Retain support for the historical compact representation, while
        # validating every count/dimension in the actual Zeo++ representation.
        compact = re.fullmatch(r"Channel dimensionality ([0-3])", first_line.strip())
        dimension = int(compact.group(1)) if compact else _probe_protocol.parse_channel_topology(text)["maximum_channel_dimension"]
    except (IndexError, ValueError, _probe_protocol.ZeoppN2HeError) as exc:
        raise ValueError(f"Could not parse Zeo++ channel output: {first_line!r}") from exc
    return {"unit": "nan", "Dimension": dimension}


def FrameworkDim(structure, high_accuracy=True, prefix="tmp_strinfo"):
    """Return framework dimensionality and counts of 1D, 2D, and 3D parts."""

    text = _run_network(structure, _arguments(high_accuracy, "-strinfo"), prefix)
    try:
        parsed = _framework_protocol.parse_strinfo(text)
    except _framework_protocol.ZeoppFrameworkDimensionError as exc:
        raise ValueError(f"Could not parse Zeo++ framework output: {text!r}") from exc
    return {
        "unit": "nan",
        "Dimension": parsed["maximum_framework_dimension"],
        "N_1D": parsed["framework_1d_count"],
        "N_2D": parsed["framework_2d_count"],
        "N_3D": parsed["framework_3d_count"],
    }


def PoreDiameter(structure, high_accuracy=True, prefix="tmp_pd"):
    """Return largest-cavity, pore-limiting, and largest-free-pore diameters."""

    text = _run_network(structure, _arguments(high_accuracy, "-res"), prefix)
    fields = text.splitlines()[0].split() if text.splitlines() else []
    try:
        if len(fields) < 4:
            raise ValueError("Missing pore-diameter fields")
        lcd, pld, lfpd = map(_finite_nonnegative, fields[-3:])
    except (IndexError, ValueError) as exc:
        raise ValueError(f"Could not parse Zeo++ pore-diameter output: {' '.join(fields)!r}") from exc
    return {"unit": "angstrom, Å", "LCD": lcd, "PLD": pld, "LFPD": lfpd}


def _labelled_float(line: str, label: str) -> float:
    try:
        matches = re.findall(r"(?:^|\s)" + re.escape(label) + r"\s*(\S+)", line)
        if len(matches) != 1:
            raise ValueError("Missing or duplicate output label")
        return _finite_nonnegative(matches[0])
    except (IndexError, ValueError) as exc:
        raise ValueError(f"Could not parse Zeo++ field {label!r} from: {line!r}") from exc


def SurfaceArea(
    structure,
    chan_radius=1.655,
    probe_radius=1.655,
    num_samples=5000,
    high_accuracy=True,
    prefix="tmp_sa",
):
    """Return accessible and non-accessible surface areas."""

    _validate_probe_options(chan_radius, probe_radius, num_samples=num_samples)
    text = _run_network(
        structure,
        _arguments(high_accuracy, "-sa", chan_radius, probe_radius, num_samples),
        prefix,
    )
    line = text.splitlines()[0] if text.splitlines() else ""
    asa = _labelled_float(line, "ASA_A^2:")
    vsa = _labelled_float(line, "ASA_m^2/cm^3:")
    gsa = _labelled_float(line, "ASA_m^2/g:")
    nasa = _labelled_float(line, "NASA_A^2:")
    nvsa = _labelled_float(line, "NASA_m^2/cm^3:")
    ngsa = _labelled_float(line, "NASA_m^2/g:")
    return {
        "unit": "Å^2, m^2/cm^3, m^2/g",
        "ASA": [asa, vsa, gsa],
        "NASA": [nasa, nvsa, ngsa],
    }


def PoreVolume(
    structure,
    chan_radius=0,
    probe_radius=0,
    num_samples=5000,
    high_accuracy=True,
    prefix="tmp_pv",
):
    """Return accessible/non-accessible pore volumes and void fractions."""

    _validate_probe_options(chan_radius, probe_radius, num_samples=num_samples)
    text = _run_network(
        structure,
        _arguments(high_accuracy, "-volpo", chan_radius, probe_radius, num_samples),
        prefix,
    )
    line = text.splitlines()[0] if text.splitlines() else ""
    poav = _labelled_float(line, "POAV_A^3:")
    ponav = _labelled_float(line, "PONAV_A^3:")
    gpoav = _labelled_float(line, "POAV_cm^3/g:")
    gponav = _labelled_float(line, "PONAV_cm^3/g:")
    poav_fraction = _labelled_float(line, "POAV_Volume_fraction:")
    ponav_fraction = _labelled_float(line, "PONAV_Volume_fraction:")
    if poav_fraction > 1 or ponav_fraction > 1:
        raise ValueError("Zeo++ void fractions must be in [0, 1]")
    return {
        "unit": "PV: Å^3, cm^3/g; VF: nan",
        "PV": [poav, gpoav],
        "NPV": [ponav, gponav],
        "VF": poav_fraction,
        "NVF": ponav_fraction,
    }
