"""Geometric-property calculations backed by the Zeo++ ``network`` binary.

Install Zeo++ independently, for example with
``conda install -c conda-forge zeopp-lsmo``.  The executable can be overridden
with the ``COREMOF_NETWORK_EXECUTABLE`` environment variable.
"""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Iterable


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

    handle = tempfile.NamedTemporaryFile(
        prefix=f"{prefix_path.name}_", suffix=".txt", dir=output_dir, delete=False
    )
    output_path = Path(handle.name)
    handle.close()
    output_path.unlink()  # Zeo++ creates the output itself.

    command = [executable, *(str(value) for value in arguments), str(output_path), str(structure_path)]
    try:
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
        if completed.returncode != 0:
            detail = completed.stderr.strip() or completed.stdout.strip() or "no diagnostic output"
            raise RuntimeError(
                f"Zeo++ failed with exit code {completed.returncode}: {detail}"
            )
        if not output_path.is_file():
            raise RuntimeError("Zeo++ completed without creating its output file")
        return output_path.read_text(encoding="utf-8")
    finally:
        output_path.unlink(missing_ok=True)


def _arguments(high_accuracy: bool, *values: object) -> list[object]:
    return (["-ha"] if high_accuracy else []) + list(values)


def ChanDim(structure, probe_radius=0, high_accuracy=True, prefix="tmp_chan"):
    """Return the dimensionality of channels accessible to a probe."""

    text = _run_network(
        structure, _arguments(high_accuracy, "-chan", probe_radius), prefix
    )
    first_line = text.splitlines()[0] if text.splitlines() else ""
    try:
        dimension = int(first_line.split("dimensionality", 1)[1].split()[0])
    except (IndexError, ValueError) as exc:
        raise ValueError(f"Could not parse Zeo++ channel output: {first_line!r}") from exc
    return {"unit": "nan", "Dimension": dimension}


def FrameworkDim(structure, high_accuracy=True, prefix="tmp_strinfo"):
    """Return framework dimensionality and counts of 1D, 2D, and 3D parts."""

    text = _run_network(structure, _arguments(high_accuracy, "-strinfo"), prefix)
    fields = text.splitlines()[0].split() if text.splitlines() else []
    try:
        dimension = int(fields[-1])
        one_dim, two_dim, three_dim = map(int, fields[7:10])
    except (IndexError, ValueError) as exc:
        raise ValueError(f"Could not parse Zeo++ framework output: {' '.join(fields)!r}") from exc
    return {
        "unit": "nan",
        "Dimension": dimension,
        "N_1D": one_dim,
        "N_2D": two_dim,
        "N_3D": three_dim,
    }


def PoreDiameter(structure, high_accuracy=True, prefix="tmp_pd"):
    """Return largest-cavity, pore-limiting, and largest-free-pore diameters."""

    text = _run_network(structure, _arguments(high_accuracy, "-res"), prefix)
    fields = text.splitlines()[0].split() if text.splitlines() else []
    try:
        lcd, pld, lfpd = map(float, fields[1:4])
    except (IndexError, ValueError) as exc:
        raise ValueError(f"Could not parse Zeo++ pore-diameter output: {' '.join(fields)!r}") from exc
    return {"unit": "angstrom, Å", "LCD": lcd, "PLD": pld, "LFPD": lfpd}


def _labelled_float(line: str, label: str) -> float:
    try:
        return float(line.split(label, 1)[1].split()[0])
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
    return {
        "unit": "PV: Å^3, cm^3/g; VF: nan",
        "PV": [poav, gpoav],
        "NPV": [ponav, gponav],
        "VF": poav_fraction,
        "NVF": ponav_fraction,
    }
