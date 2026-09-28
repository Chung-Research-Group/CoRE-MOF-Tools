#!/usr/bin/env python3
"""Audit-first v4 candidate curation for the 2026 CSD MOF update.

This driver intentionally does not replace ``curate_csd_mof_update.py`` or any
v2/v3 result.  It corrects three defects discovered in the documented-clean
baseline:

* ASE adds ``skin`` to every atomic cutoff, so a requested total pair margin
  must be passed to ASE as ``skin = total_pair_margin / 2``.
* Every adaptive iteration is written to its own deterministic directory and
  the path returned by that iteration is used directly; file mtimes never
  select scientific state.
* FSR is generated first.  Non-ionic ASR is then generated from the final FSR,
  making it possible to enforce ``ASR subset FSR subset preprocessed input``.

Confirmed ion structures retain only the ion-preserving FSR candidate.  If an
ion is detected only after the ASR bond-cutting pass, the structure is routed
to REVIEW rather than silently changing category.  Distinct non-ionic FSR and
ASR structures are both retained, as described by the CoRE MOF manuscript;
identical roles share one physical deliverable.

PACMAN is run only on final candidate deliverables, with ``neutral=True``, a
fixed PyTorch seed, deterministic algorithms, one CPU thread, and CUDA hidden.
The explicit seed is required because PACMAN-charge 1.4.2 constructs a random
auxiliary encoder during prediction.  Every charged CIF must contain one
charge per atom, have a net charge within the configured tolerance, and
preserve the uncharged atom set.  ``--skip-pacman`` is provided for
solvent-removal pilots; such records are never marked as curation-stage
eligible.

Dimensionality, porosity, Chen--Manz, MOFChecker, MOSAEC, and MOFClassifier are
recorded as explicit deferred gates.  This script does not assign final CR/NCR
labels and must write to a new isolated output directory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import re
import shutil
import sys
import tempfile
import traceback
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable


PIPELINE_VERSION = (
    "2026.07.16-4.1-corrected-margin-fsr-derived-asr-seeded-pacman-candidate"
)
SOLVENT_CLEANER = "CoREMOF.curate.clean (v4 controlled adapter)"
INITIAL_TOTAL_PAIR_MARGIN_ANGSTROM = 0.25
ADAPTIVE_TOTAL_MARGIN_INCREMENT_ANGSTROM = 0.05
PACMAN_DIGITS = 10
PACMAN_TORCH_SEED = 0
PACMAN_DEVICE_POLICY = "cpu"
REFERENCE_NUM_THREADS = 1

# Interim, deliberately small registry for automatic removal by the legacy
# component selector.  Formula-only recognition is not sufficient for release;
# unknown removals are routed to review, and the immutable/preprocessed inputs
# remain available for component-level re-curation.
NEUTRAL_SOLVENT_REGISTRY_REVISION = "2026-07-16.1-interim"
KNOWN_NEUTRAL_SOLVENT_FORMULAS = frozenset({
    "C2H3N",     # acetonitrile
    "C2H6O",     # ethanol
    "C2H6OS",    # dimethyl sulfoxide
    "C3H6O",     # acetone
    "C3H7NO",    # dimethylformamide
    "C4H10O",    # diethyl ether
    "C4H8O",     # tetrahydrofuran
    "C5H11NO",   # diethylformamide
    "C5H5N",     # pyridine
    "C6H6",      # benzene
    "C7H8",      # toluene
    "CH2Cl2",    # dichloromethane
    "CH4O",      # methanol
    "CHCl3",     # chloroform
    "H2O",       # water
})


class AdaptiveMarginExhausted(RuntimeError):
    """Raised when removed metal-containing components never disappear."""


class AtomMappingError(RuntimeError):
    """Raised when a cleaner output cannot be mapped to its parent atoms."""


def configure_reference_runtime() -> dict[str, Any]:
    """Configure the deterministic CPU runtime before scientific imports.

    PACMAN-charge chooses its device when ``pmcharge`` is imported and creates
    a randomly initialized auxiliary GCN inside every prediction.  Therefore
    device visibility and the PyTorch seed are part of the scientific method,
    not merely performance settings.
    """
    for variable in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[variable] = str(REFERENCE_NUM_THREADS)
    if PACMAN_DEVICE_POLICY == "cpu":
        os.environ["CUDA_VISIBLE_DEVICES"] = ""

    import torch

    if PACMAN_DEVICE_POLICY == "cpu" and torch.cuda.is_available():
        raise RuntimeError(
            "reference CPU policy failed: torch still reports CUDA available; "
            "launch a fresh process with CUDA_VISIBLE_DEVICES=''"
        )
    torch.set_num_threads(REFERENCE_NUM_THREADS)
    try:
        torch.set_num_interop_threads(REFERENCE_NUM_THREADS)
    except RuntimeError:
        if torch.get_num_interop_threads() != REFERENCE_NUM_THREADS:
            raise
    torch.use_deterministic_algorithms(True)
    if hasattr(torch.backends, "cudnn"):
        torch.backends.cudnn.benchmark = False
        torch.backends.cudnn.deterministic = True
    torch.manual_seed(PACMAN_TORCH_SEED)

    return {
        "pacman_device_policy": PACMAN_DEVICE_POLICY,
        "torch_device": "cpu",
        "torch_version": torch.__version__,
        "torch_cuda_build": torch.version.cuda,
        "torch_cuda_available": torch.cuda.is_available(),
        "torch_seed": PACMAN_TORCH_SEED,
        "torch_deterministic_algorithms": (
            torch.are_deterministic_algorithms_enabled()
        ),
        "torch_num_threads": torch.get_num_threads(),
        "torch_num_interop_threads": torch.get_num_interop_threads(),
        "environment": {
            key: os.environ.get(key)
            for key in (
                "CUDA_VISIBLE_DEVICES",
                "OMP_NUM_THREADS",
                "MKL_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "NUMEXPR_NUM_THREADS",
                "PYTHONHASHSEED",
                "LC_ALL",
                "LANG",
            )
        },
    }


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json_atomic(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent,
        prefix=f".{path.name}.", delete=False,
    ) as out:
        json.dump(data, out, indent=2, sort_keys=True, default=str)
        out.write("\n")
        temporary = Path(out.name)
    os.replace(temporary, path)


def copy_atomic(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=destination.parent, prefix=f".{destination.name}."
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        shutil.copy2(source, temporary)
        os.replace(temporary, destination)
    finally:
        if temporary.exists():
            temporary.unlink()


def reset_directory(path: Path) -> None:
    """Reset only a refcode-local v4 work directory."""
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def load_csd_metadata(path: Path | None) -> dict[str, dict[str, Any]]:
    if path is None:
        return {}
    result: dict[str, dict[str, Any]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                item = json.loads(line)
                refcode = item.get("refcode")
                if refcode:
                    result[refcode] = item
    return result


def cif_inputs(root: Path) -> list[Path]:
    return sorted(path for path in root.rglob("*.cif") if path.is_file())


def record_path(state_dir: Path, refcode: str) -> Path:
    return state_dir / "records" / f"{refcode}.json"


def checker_result(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"status": "missing"}
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def package_version(*names: str) -> str | None:
    for name in names:
        try:
            return importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            continue
    return None


def git_commit_from_path(source: Path) -> str | None:
    """Read a nearby Git HEAD without invoking Git or changing repository state."""
    for directory in (source.resolve(), *source.resolve().parents):
        git_dir = directory / ".git"
        if not git_dir.is_dir():
            continue
        head = git_dir / "HEAD"
        if not head.exists():
            return None
        value = head.read_text(encoding="utf-8").strip()
        if value.startswith("ref: "):
            ref_path = git_dir / value[5:]
            if ref_path.exists():
                return ref_path.read_text(encoding="utf-8").strip()
            packed = git_dir / "packed-refs"
            if packed.exists():
                target = value[5:]
                for line in packed.read_text(encoding="utf-8").splitlines():
                    if line and not line.startswith(("#", "^")):
                        commit, ref = line.split(" ", 1)
                        if ref == target:
                            return commit
            return None
        return value
    return None


def scientific_provenance(coremof_curate_path: Path, ase_neighborlist_path: Path,
                          ion_list_path: Path, radii_path: Path) -> dict[str, Any]:
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "packages": {
            "ase": package_version("ase"),
            "coremof_tools": package_version("coremof-tools", "CoREMOF"),
            "gemmi": package_version("gemmi"),
            "pacman_charge": package_version("PACMAN-charge", "PACMANCharge"),
            "pymatgen": package_version("pymatgen"),
        },
        "coremof_git_commit": git_commit_from_path(coremof_curate_path.parent),
        "source_files": {
            "coremof_curate": {
                "path": str(coremof_curate_path),
                "sha256": sha256(coremof_curate_path),
            },
            "coremof_ions_list": {
                "path": str(ion_list_path),
                "sha256": sha256(ion_list_path),
            },
            "coremof_atomic_definitions": {
                "path": str(radii_path),
                "sha256": sha256(radii_path),
            },
            "ase_neighborlist": {
                "path": str(ase_neighborlist_path),
                "sha256": sha256(ase_neighborlist_path),
            },
            "v4_driver": {
                "path": str(Path(__file__).resolve()),
                "sha256": sha256(Path(__file__).resolve()),
            },
        },
    }


def corrected_cleaner() -> tuple[Any, dict[str, Path]]:
    """Return an upstream cleaner adapter with corrected pair-margin semantics."""
    import ase.neighborlist as ase_neighborlist_module
    import CoREMOF.curate as coremof_curate_module
    import CoREMOF.utils.atoms_definitions as atoms_definitions_module
    import CoREMOF.utils.ions_list as ions_list_module
    from ase.neighborlist import NeighborList
    from CoREMOF.curate import clean as UpstreamClean
    from CoREMOF.utils.atoms_definitions import COVALENTRADII, METAL
    from CoREMOF.utils.ions_list import ALLIONS

    class CorrectedMarginClean(UpstreamClean):
        """Use total pair margin/2 as ASE's per-atom skin."""

        def __init__(self) -> None:
            # Deliberately do not call UpstreamClean.__init__, which would run
            # both independent branches and overwrite iteration outputs.
            self.cambridge_radii = COVALENTRADII
            self.metal_list = [element for element, is_metal in METAL.items() if is_metal]
            self.ions_list = set(ALLIONS)

        def build_ASE_neighborlist(self, cif: Any, skin: float) -> Any:
            radii = [self.cambridge_radii[symbol]
                     for symbol in cif.get_chemical_symbols()]
            neighborlist = NeighborList(
                radii,
                self_interaction=False,
                bothways=True,
                # ASE adds skin to EACH atomic cutoff.  Halving here yields
                # r_i + r_j + requested_total_pair_margin.
                skin=skin / 2.0,
            )
            neighborlist.update(cif)
            return neighborlist

    paths = {
        "coremof_curate": Path(coremof_curate_module.__file__).resolve(),
        "ase_neighborlist": Path(ase_neighborlist_module.__file__).resolve(),
        "ions_list": Path(ions_list_module.__file__).resolve(),
        "atomic_definitions": Path(atoms_definitions_module.__file__).resolve(),
    }
    return CorrectedMarginClean(), paths


def formula_elements(formula: str) -> list[str]:
    return [match[0] for match in re.findall(r"([A-Z][a-z]?)(\d*)", formula)]


def removed_formula_has_metal(formulas: list[str], metals: set[str]) -> bool:
    return any(any(element in metals for element in formula_elements(formula))
               for formula in formulas)


def metal_counts(atoms: Any, metals: set[str]) -> dict[str, int]:
    """Return deterministic element counts for atoms considered metals."""
    counts = Counter(
        symbol for symbol in atoms.get_chemical_symbols() if symbol in metals
    )
    return dict(sorted(counts.items()))


def margin_token(value: float) -> str:
    return f"{value:.3f}".replace(".", "p")


def unrecognized_removed_formulas(stage_result: dict[str, Any]) -> list[str]:
    """Return final-iteration removals outside the interim solvent registry."""
    final_index = int(stage_result["final_iteration"])
    formulas = stage_result["iterations"][final_index][
        "removed_component_formulas"
    ]
    return sorted(
        formula for formula in formulas
        if formula not in KNOWN_NEUTRAL_SOLVENT_FORMULAS
    )


def run_adaptive_stage(
    cleaner: Any,
    input_cif: Path,
    stage_root: Path,
    stage: str,
    initial_total_margin: float,
    max_iterations: int,
    ase_read: Callable[..., Any],
) -> dict[str, Any]:
    """Run one branch with explicit deterministic iteration directories."""
    if stage not in {"FSR", "ASR"}:
        raise ValueError(f"unknown cleaning stage {stage!r}")
    method = cleaner.free_clean if stage == "FSR" else cleaner.all_clean
    input_stem = str(input_cif.with_suffix(""))
    iterations: list[dict[str, Any]] = []
    metals = set(cleaner.metal_list)

    for index in range(max_iterations):
        total_margin = round(
            initial_total_margin
            + index * ADAPTIVE_TOTAL_MARGIN_INCREMENT_ANGSTROM,
            10,
        )
        iteration_dir = (
            stage_root
            / f"iter_{index:03d}_total_margin_{margin_token(total_margin)}"
        )
        iteration_dir.mkdir(parents=True, exist_ok=False)
        result = method(
            input_stem,
            str(iteration_dir),
            cleaner.ions_list,
            total_margin,
        )
        if not isinstance(result, tuple) or len(result) != 2:
            raise RuntimeError(
                f"upstream {stage} cleaner returned no structured result at "
                f"total margin {total_margin:.3f} A"
            )
        _, raw_formulas = result
        removed_formulas = sorted(str(value) for value in (raw_formulas or []))
        candidates = sorted(iteration_dir.glob("*.cif"))
        if len(candidates) != 1:
            raise RuntimeError(
                f"{stage} iteration {index} produced {len(candidates)} CIFs; expected one"
            )
        candidate = candidates[0]
        if stage == "FSR" and not (
            candidate.stem.endswith("_FSR")
            or candidate.stem.endswith("_ION_FSR")
        ):
            raise RuntimeError(f"unexpected FSR filename: {candidate.name}")
        if stage == "ASR" and not (
            candidate.stem.endswith("_ASR")
            or candidate.stem.endswith("_ION_ASR")
        ):
            raise RuntimeError(f"unexpected ASR filename: {candidate.name}")

        atoms = ase_read(str(candidate))
        has_removed_metal = removed_formula_has_metal(removed_formulas, metals)
        iteration_record = {
            "iteration": index,
            "requested_total_pair_margin_angstrom": total_margin,
            "ase_per_atom_skin_angstrom": total_margin / 2.0,
            "effective_pair_threshold": "r_i + r_j + requested_total_pair_margin",
            "candidate_cif": str(candidate),
            "candidate_sha256": sha256(candidate),
            "candidate_atom_count": len(atoms),
            "candidate_formula": atoms.get_chemical_formula(),
            "removed_component_formulas": removed_formulas,
            "removed_formula_contains_metal": has_removed_metal,
        }
        iterations.append(iteration_record)
        if not has_removed_metal:
            return {
                "stage": stage,
                "input_cif": str(input_cif),
                "initial_total_pair_margin_angstrom": initial_total_margin,
                "adaptive_total_margin_increment_angstrom": (
                    ADAPTIVE_TOTAL_MARGIN_INCREMENT_ANGSTROM
                ),
                "iterations": iterations,
                "final_iteration": index,
                "final_total_pair_margin_angstrom": total_margin,
                "final_cif": str(candidate),
                "final_has_ion_suffix": f"_ION_{stage}" in candidate.stem,
            }

    raise AdaptiveMarginExhausted(
        f"{stage} still removed a metal-containing component after "
        f"{max_iterations} iterations"
    )


def map_child_atoms(parent_path: Path, child_path: Path, tolerance_angstrom: float,
                    ase_read: Callable[..., Any]) -> dict[str, Any]:
    """Map every child atom to one parent atom using element and periodic position."""
    import numpy as np
    from ase.geometry import find_mic

    parent = ase_read(str(parent_path))
    child = ase_read(str(child_path))
    parent_cell = np.asarray(parent.cell.array, dtype=float)
    child_cell = np.asarray(child.cell.array, dtype=float)
    if not np.allclose(parent_cell, child_cell, rtol=1e-8, atol=1e-5):
        raise AtomMappingError(
            f"cell changed between {parent_path.name} and {child_path.name}"
        )
    if list(parent.get_pbc()) != list(child.get_pbc()):
        raise AtomMappingError(
            f"PBC changed between {parent_path.name} and {child_path.name}"
        )
    if len(child) > len(parent):
        raise AtomMappingError(
            f"child {child_path.name} has {len(child)} atoms; parent has {len(parent)}"
        )

    parent_symbols = parent.get_chemical_symbols()
    child_symbols = child.get_chemical_symbols()
    unused = set(range(len(parent)))
    mapping: list[int] = []
    ambiguous: list[dict[str, Any]] = []
    distances: list[float] = []
    for child_index, (symbol, position) in enumerate(
        zip(child_symbols, child.positions)
    ):
        matches: list[tuple[float, int]] = []
        for parent_index in sorted(unused):
            if parent_symbols[parent_index] != symbol:
                continue
            delta, distance = find_mic(
                position - parent.positions[parent_index],
                parent.cell,
                pbc=parent.get_pbc(),
            )
            del delta
            numeric_distance = float(distance)
            if numeric_distance <= tolerance_angstrom:
                matches.append((numeric_distance, parent_index))
        if not matches:
            raise AtomMappingError(
                f"cannot map child atom {child_index} ({symbol}) from "
                f"{child_path.name} into {parent_path.name} within "
                f"{tolerance_angstrom:g} A"
            )
        matches.sort()
        distance, selected = matches[0]
        if len(matches) > 1:
            ambiguous.append({
                "child_index": child_index,
                "symbol": symbol,
                "candidate_parent_indices": [item[1] for item in matches],
                "candidate_distances_angstrom": [item[0] for item in matches],
                "selected_parent_index": selected,
            })
        mapping.append(selected)
        distances.append(distance)
        unused.remove(selected)

    return {
        "parent_cif": str(parent_path),
        "child_cif": str(child_path),
        "parent_atom_count": len(parent),
        "child_atom_count": len(child),
        "child_to_parent_indices": mapping,
        "unique_parent_indices": len(set(mapping)) == len(mapping),
        "ambiguous_matches": ambiguous,
        "mapping_unambiguous": not ambiguous,
        "maximum_mapping_distance_angstrom": max(distances, default=0.0),
        "tolerance_angstrom": tolerance_angstrom,
    }


def numeric_cif_value(value: str) -> float:
    match = re.match(
        r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?",
        value.strip(),
    )
    if not match:
        raise ValueError(f"not a numeric CIF value: {value!r}")
    return float(match.group(0))


def charge_statistics(cif_path: Path, expected_atoms: int) -> dict[str, Any]:
    import gemmi

    document = gemmi.cif.read_file(str(cif_path))
    block = document.sole_block()
    values = list(block.find_values("_atom_site_charge"))
    if len(values) != expected_atoms:
        raise RuntimeError(
            f"charged CIF has {len(values)} _atom_site_charge values; "
            f"expected {expected_atoms}"
        )
    charges = [numeric_cif_value(value) for value in values]
    return {
        "charge_tag": "_atom_site_charge",
        "charge_count": len(charges),
        "net_charge_e": math.fsum(charges),
        "minimum_charge_e": min(charges),
        "maximum_charge_e": max(charges),
    }


def run_neutral_pacman(
    uncharged: Path,
    destination: Path,
    charge_work_dir: Path,
    mapping_tolerance: float,
    net_charge_tolerance: float,
    ase_read: Callable[..., Any],
) -> dict[str, Any]:
    import torch
    from PACMANCharge import pmcharge

    if str(pmcharge.device) != PACMAN_DEVICE_POLICY:
        raise RuntimeError(
            f"PACMAN selected device {pmcharge.device!s}; expected "
            f"{PACMAN_DEVICE_POLICY!r}"
        )
    # PACMAN-charge 1.4.2 constructs an auxiliary GCN without loading its
    # weights from the checkpoint.  Resetting the seed immediately before each
    # prediction makes each deliverable reproducible and independent of the
    # number/order of other deliverables processed by the same interpreter.
    torch.manual_seed(PACMAN_TORCH_SEED)
    charge_work_dir.mkdir(parents=True, exist_ok=False)
    charge_input = charge_work_dir / uncharged.name
    copy_atomic(uncharged, charge_input)
    pmcharge.predict(
        cif_file=str(charge_input),
        charge_type="DDEC6",
        digits=PACMAN_DIGITS,
        atom_type=True,
        neutral=True,
        # Preserve the candidate input; PyCifRW adds the charge column to the
        # generated output without pymatgen rewriting the source structure.
        keep_connect=True,
    )
    generated = charge_input.with_name(f"{charge_input.stem}_pacman.cif")
    if not generated.exists():
        raise RuntimeError("PACMAN did not create the charged candidate CIF")

    uncharged_atoms = ase_read(str(uncharged))
    charged_atoms = ase_read(str(generated))
    if uncharged_atoms.get_chemical_symbols() != charged_atoms.get_chemical_symbols():
        raise RuntimeError("PACMAN output does not preserve the atom sequence")
    charged_mapping = map_child_atoms(
        uncharged, generated, mapping_tolerance, ase_read
    )
    if charged_mapping["child_atom_count"] != charged_mapping["parent_atom_count"]:
        raise RuntimeError("PACMAN output changed the atom count")
    statistics = charge_statistics(generated, len(uncharged_atoms))
    net_charge_passed = abs(statistics["net_charge_e"]) <= net_charge_tolerance
    if not net_charge_passed:
        raise RuntimeError(
            f"neutralized PACMAN net charge {statistics['net_charge_e']:.12g} e "
            f"exceeds tolerance {net_charge_tolerance:g} e"
        )
    copy_atomic(generated, destination)
    return {
        "status": "complete",
        "charge_type": "DDEC6",
        "neutral": True,
        "digits": PACMAN_DIGITS,
        "atom_type": True,
        "keep_connect": True,
        "torch_seed": PACMAN_TORCH_SEED,
        "torch_device": str(pmcharge.device),
        "torch_deterministic_algorithms": (
            torch.are_deterministic_algorithms_enabled()
        ),
        "torch_num_threads": torch.get_num_threads(),
        "charged_cif": str(destination),
        "charged_sha256": sha256(destination),
        "atom_mapping": charged_mapping,
        "net_charge_tolerance_e": net_charge_tolerance,
        "net_charge_check": "passed",
        **statistics,
    }


def build_deliverable_specs(
    category: str,
    fsr_cif: Path,
    asr_cif: Path | None,
    fsr_input_indices: list[int],
    asr_input_indices: list[int] | None,
    roles_identical: bool,
) -> list[dict[str, Any]]:
    if category == "ION":
        return [{
            "token": "ION_FSR",
            "roles": ["ION_FSR"],
            "source_candidate": fsr_cif,
            "input_atom_indices": fsr_input_indices,
            "deduplicated_identical_roles": False,
        }]
    if asr_cif is None or asr_input_indices is None:
        return [{
            "token": "FSR_REVIEW",
            "roles": ["FSR"],
            "source_candidate": fsr_cif,
            "input_atom_indices": fsr_input_indices,
            "deduplicated_identical_roles": False,
        }]
    if roles_identical:
        return [{
            "token": "FSR_ASR_IDENTICAL",
            "roles": ["FSR", "ASR"],
            "source_candidate": fsr_cif,
            "input_atom_indices": fsr_input_indices,
            "deduplicated_identical_roles": True,
        }]
    return [
        {
            "token": "FSR",
            "roles": ["FSR"],
            "source_candidate": fsr_cif,
            "input_atom_indices": fsr_input_indices,
            "deduplicated_identical_roles": False,
        },
        {
            "token": "ASR",
            "roles": ["ASR"],
            "source_candidate": asr_cif,
            "input_atom_indices": asr_input_indices,
            "deduplicated_identical_roles": False,
        },
    ]


def run_one(
    source: Path,
    output: Path,
    metadata: dict[str, dict[str, Any]],
    overwrite: bool = False,
    skip_pacman: bool = False,
    initial_total_margin: float = INITIAL_TOTAL_PAIR_MARGIN_ANGSTROM,
    max_adaptive_iterations: int = 40,
    mapping_tolerance: float = 1.0e-3,
    net_charge_tolerance: float = 1.0e-5,
) -> dict[str, Any]:
    """Process one immutable source CIF and return a complete audit record."""
    refcode = source.stem
    state_dir = output / "state"
    done_path = record_path(state_dir, refcode)
    if done_path.exists() and not overwrite:
        with done_path.open(encoding="utf-8") as handle:
            prior = json.load(handle)
        prior_version = prior.get("pipeline_version")
        if prior_version != PIPELINE_VERSION:
            raise RuntimeError(
                f"output contains {prior_version!r}; use a new output directory "
                f"for {PIPELINE_VERSION!r}"
            )
        if prior.get("status") in {"complete", "review"}:
            return prior

    # Configure device visibility, threading, and RNG before any package can
    # import PACMAN/torch.  Heavy imports remain local so --help and manifest
    # inspection work in a lightweight shell environment.
    runtime_provenance = configure_reference_runtime()

    # Heavy imports remain local so --help and manifest inspection work in a
    # lightweight shell environment.
    from ase.io import read as ase_read
    from CoREMOF.curate import preprocess

    cleaner, source_paths = corrected_cleaner()

    work = output / "work" / refcode
    candidate_ref_dir = output / "candidate_cifs" / refcode
    primary_ref_dir = output / "primary_cifs" / refcode
    reset_directory(work)
    reset_directory(candidate_ref_dir)
    reset_directory(primary_ref_dir)
    preprocessed_dir = work / "preprocessed"
    fsr_iterations_dir = work / "iterations" / "fsr"
    asr_iterations_dir = work / "iterations" / "asr_from_final_fsr"
    charge_work_root = work / "pacman"
    for directory in (
        preprocessed_dir,
        fsr_iterations_dir,
        asr_iterations_dir,
        charge_work_root,
    ):
        directory.mkdir(parents=True, exist_ok=True)

    source_meta = metadata.get(refcode)
    record: dict[str, Any] = {
        "pipeline_version": PIPELINE_VERSION,
        "refcode": refcode,
        "source_cif": str(source),
        "source_sha256": sha256(source),
        "ccdc_metadata": source_meta,
        "started_at": utcnow(),
        "status": "running",
        "category": None,
        "curation": {
            "solvent_cleaner": SOLVENT_CLEANER,
            "neighbor_policy": {
                "radii": "CoREMOF COVALENTRADII",
                "initial_requested_total_pair_margin_angstrom": initial_total_margin,
                "ase_per_atom_skin": "requested_total_pair_margin / 2",
                "effective_pair_threshold": (
                    "r_i + r_j + requested_total_pair_margin"
                ),
                "adaptive_total_margin_increment_angstrom": (
                    ADAPTIVE_TOTAL_MARGIN_INCREMENT_ANGSTROM
                ),
                "max_adaptive_iterations": max_adaptive_iterations,
            },
            "branch_policy": (
                "build FSR from preprocessed input; detect ion from final FSR; "
                "skip ASR for confirmed ion; otherwise build ASR from final FSR"
            ),
            "ion_policy": (
                "final ION_FSR suffix confirms ion; ASR-only ion evidence routes "
                "to REVIEW"
            ),
            "deliverable_policy": (
                "retain distinct non-ion FSR and ASR; one physical deliverable "
                "represents both roles when atom sets are identical"
            ),
            "component_removal_policy": {
                "neutral_solvent_registry_revision": (
                    NEUTRAL_SOLVENT_REGISTRY_REVISION
                ),
                "known_neutral_solvent_formulas": sorted(
                    KNOWN_NEUTRAL_SOLVENT_FORMULAS
                ),
                "unknown_removal_action": (
                    "route to REVIEW; do not charge or release the legacy-cleaner "
                    "candidate; immutable and preprocessed inputs remain authoritative"
                ),
                "limitation": (
                    "formula-only registry is an interim guard, not final chemical "
                    "component perception"
                ),
            },
            "pacman_policy": {
                "charge_type": "DDEC6",
                "neutral": True,
                "digits": PACMAN_DIGITS,
                "device_policy": PACMAN_DEVICE_POLICY,
                "torch_seed": PACMAN_TORCH_SEED,
                "reference_num_threads": REFERENCE_NUM_THREADS,
                "deterministic_algorithms": True,
                "net_charge_tolerance_e": net_charge_tolerance,
                "skipped_by_request": skip_pacman,
            },
        },
        "runtime_provenance": runtime_provenance,
        "software_provenance": scientific_provenance(
            source_paths["coremof_curate"],
            source_paths["ase_neighborlist"],
            source_paths["ions_list"],
            source_paths["atomic_definitions"],
        ),
        "eligibility": {
            "framework_dimensionality": {
                "status": "deferred",
                "policy_placeholder": "manuscript 2D/3D framework gate",
            },
            "porosity": {
                "status": "deferred",
                "policy_placeholder": "manuscript PLD > 2.4 A gate",
                "source_csd_pore_analyser_num_percolated_dimensions": (
                    (source_meta or {}).get(
                        "pore_analyser_num_percolated_dimensions"
                    )
                ),
                "source_csd_pore_limiting_diameter_angstrom": (
                    (source_meta or {}).get("pore_analyser_pore_limiting_diameter")
                ),
            },
        },
        "checkers": {
            "Chen_Manz": {"status": "deferred"},
            "MOFChecker": {"status": "deferred"},
            "MOSAEC": {"status": "deferred", "requires_csd_license": True},
            "MOFClassifier": {"status": "deferred"},
            "final_cr_ncr_policy": {"status": "deferred_to_v12_merge_policy"},
        },
        "outputs": {},
        "atom_mapping": {},
        "invariants": {},
        "review_reasons": [],
        "deliverables": [],
        "errors": [],
    }

    try:
        try:
            preprocess(str(source), output_folder=str(preprocessed_dir))
        except Exception as exc:
            if type(exc).__name__ != "SymmetryUndetermined":
                raise
            fallback = preprocessed_dir / source.name
            shutil.copy2(source, fallback)
            record["precheck_fallback"] = {
                "reason": "SymmetryUndetermined",
                "preprocessed_cif": str(fallback),
            }
        precheck_path = preprocessed_dir / f"{refcode}_precheck.json"
        record["precheck"] = checker_result(precheck_path)
        prepared = sorted(preprocessed_dir.glob(f"{refcode}*.cif"))
        if not prepared:
            raise RuntimeError("preprocess produced no CIF")
        if len(prepared) != 1:
            raise RuntimeError(
                f"preprocess produced {len(prepared)} structures; split review required"
            )
        prepared_cif = prepared[0]
        prepared_atoms = ase_read(str(prepared_cif))
        record["outputs"]["preprocessed_cif"] = str(prepared_cif)
        record["outputs"]["preprocessed_sha256"] = sha256(prepared_cif)
        record["validation"] = {
            "preprocessed_ase_parse": "passed",
            "preprocessed_atom_count": len(prepared_atoms),
            "preprocessed_formula": prepared_atoms.get_chemical_formula(),
        }

        fsr_run = run_adaptive_stage(
            cleaner,
            prepared_cif,
            fsr_iterations_dir,
            "FSR",
            initial_total_margin,
            max_adaptive_iterations,
            ase_read,
        )
        record["curation"]["fsr"] = fsr_run
        fsr_cif = Path(fsr_run["final_cif"])
        fsr_atoms = ase_read(str(fsr_cif))
        fsr_mapping = map_child_atoms(
            prepared_cif, fsr_cif, mapping_tolerance, ase_read
        )
        record["atom_mapping"]["fsr_to_preprocessed"] = fsr_mapping
        fsr_input_indices = fsr_mapping["child_to_parent_indices"]
        fsr_is_ion = bool(fsr_run["final_has_ion_suffix"])
        fsr_unrecognized_removals = unrecognized_removed_formulas(fsr_run)
        if fsr_unrecognized_removals:
            record["review_reasons"].append(
                "FSR legacy cleaner removed unrecognized component formula(s): "
                + ", ".join(fsr_unrecognized_removals)
            )
        record["outputs"]["final_fsr_cif"] = str(fsr_cif)

        asr_cif: Path | None = None
        asr_input_indices: list[int] | None = None
        asr_is_ion = False
        asr_mapping: dict[str, Any] | None = None
        asr_atoms: Any | None = None
        if fsr_is_ion:
            record["curation"]["asr"] = {
                "status": "skipped",
                "reason": "confirmed ion in final FSR",
            }
            category = "ION"
        else:
            asr_run = run_adaptive_stage(
                cleaner,
                fsr_cif,
                asr_iterations_dir,
                "ASR",
                initial_total_margin,
                max_adaptive_iterations,
                ase_read,
            )
            record["curation"]["asr"] = asr_run
            asr_cif = Path(asr_run["final_cif"])
            asr_atoms = ase_read(str(asr_cif))
            record["outputs"]["final_asr_cif"] = str(asr_cif)
            asr_is_ion = bool(asr_run["final_has_ion_suffix"])
            asr_mapping = map_child_atoms(
                fsr_cif, asr_cif, mapping_tolerance, ase_read
            )
            record["atom_mapping"]["asr_to_fsr"] = asr_mapping
            asr_to_fsr = asr_mapping["child_to_parent_indices"]
            asr_input_indices = [fsr_input_indices[index] for index in asr_to_fsr]
            record["atom_mapping"]["asr_to_preprocessed_indices"] = (
                asr_input_indices
            )
            if asr_is_ion:
                category = "AMBIGUOUS"
                record["review_reasons"].append(
                    "ion evidence appears only after ASR bond cutting"
                )
            else:
                category = "FSR_ASR"
            asr_unrecognized_removals = unrecognized_removed_formulas(asr_run)
            if asr_unrecognized_removals:
                record["review_reasons"].append(
                    "ASR legacy cleaner removed unrecognized component formula(s): "
                    + ", ".join(asr_unrecognized_removals)
                )

        fsr_unique = (
            fsr_mapping["unique_parent_indices"]
            and len(fsr_input_indices) == len(set(fsr_input_indices))
        )
        fsr_subset = all(
            0 <= index < fsr_mapping["parent_atom_count"]
            for index in fsr_input_indices
        )
        metals = set(cleaner.metal_list)
        preprocessed_metal_counts = metal_counts(prepared_atoms, metals)
        fsr_metal_counts = metal_counts(fsr_atoms, metals)
        fsr_preserves_metal_counts = (
            fsr_metal_counts == preprocessed_metal_counts
        )
        if asr_mapping is None:
            asr_unique = None
            asr_subset = None
            asr_count_order = None
            roles_identical = False
            asr_preserves_metal_counts = None
            asr_metal_counts = None
        else:
            asr_to_fsr = asr_mapping["child_to_parent_indices"]
            asr_unique = (
                asr_mapping["unique_parent_indices"]
                and len(asr_to_fsr) == len(set(asr_to_fsr))
            )
            asr_subset = all(
                0 <= index < asr_mapping["parent_atom_count"]
                for index in asr_to_fsr
            )
            asr_count_order = (
                asr_mapping["child_atom_count"]
                <= asr_mapping["parent_atom_count"]
                <= fsr_mapping["parent_atom_count"]
            )
            roles_identical = (
                asr_mapping["child_atom_count"]
                == asr_mapping["parent_atom_count"]
                and set(asr_to_fsr)
                == set(range(asr_mapping["parent_atom_count"]))
            )
            asr_metal_counts = metal_counts(asr_atoms, metals)
            asr_preserves_metal_counts = asr_metal_counts == fsr_metal_counts

        invariant_checks = {
            "fsr_unique_atom_mapping": fsr_unique,
            "fsr_subset_of_preprocessed": fsr_subset,
            "fsr_atom_count_not_greater_than_preprocessed": (
                fsr_mapping["child_atom_count"]
                <= fsr_mapping["parent_atom_count"]
            ),
            "fsr_mapping_unambiguous": fsr_mapping["mapping_unambiguous"],
            "fsr_preserves_metal_counts": fsr_preserves_metal_counts,
            "asr_unique_atom_mapping": asr_unique,
            "asr_subset_of_fsr": asr_subset,
            "asr_fsr_input_atom_count_order": asr_count_order,
            "asr_mapping_unambiguous": (
                asr_mapping["mapping_unambiguous"]
                if asr_mapping is not None else None
            ),
            "asr_preserves_metal_counts": asr_preserves_metal_counts,
            "fsr_asr_roles_identical": roles_identical,
            "ion_skips_asr": (not fsr_is_ion) or asr_cif is None,
        }
        required_invariants = [
            value for key, value in invariant_checks.items()
            if key != "fsr_asr_roles_identical" and value is not None
        ]
        invariants_passed = all(required_invariants)
        invariant_checks["all_required_invariants_passed"] = invariants_passed
        record["invariants"] = invariant_checks
        record["validation"]["metal_counts"] = {
            "preprocessed": preprocessed_metal_counts,
            "fsr": fsr_metal_counts,
            "asr": asr_metal_counts,
            "policy": (
                "any metal-count change is conservatively routed to review; "
                "framework-vs-guest metal assignment requires later component evidence"
            ),
        }
        if not invariants_passed:
            record["review_reasons"].append(
                "one or more atom-mapping/subset invariants did not pass"
            )

        record["category"] = category
        specs = build_deliverable_specs(
            category,
            fsr_cif,
            asr_cif,
            fsr_input_indices,
            asr_input_indices,
            roles_identical,
        )
        for position, spec in enumerate(specs):
            token = spec["token"]
            uncharged = candidate_ref_dir / f"{refcode}_{token}.cif"
            copy_atomic(spec["source_candidate"], uncharged)
            atoms = ase_read(str(uncharged))
            deliverable: dict[str, Any] = {
                "id": token,
                "roles": spec["roles"],
                "deduplicated_identical_roles": (
                    spec["deduplicated_identical_roles"]
                ),
                "uncharged_cif": str(uncharged),
                "uncharged_sha256": sha256(uncharged),
                "atom_count": len(atoms),
                "formula": atoms.get_chemical_formula(),
                "input_atom_indices": spec["input_atom_indices"],
                "release_eligible": False,
            }
            if record["review_reasons"]:
                deliverable["pacman"] = {
                    "status": "not_run_review",
                    "reason": "; ".join(record["review_reasons"]),
                }
                deliverable["curation_stage_eligible"] = False
            elif skip_pacman:
                deliverable["pacman"] = {
                    "status": "skipped_by_request",
                    "neutral_policy_required_for_production": True,
                }
                deliverable["curation_stage_eligible"] = False
            else:
                charged = primary_ref_dir / f"{refcode}_{token}_pacman.cif"
                pacman_result = run_neutral_pacman(
                    uncharged,
                    charged,
                    charge_work_root / f"deliverable_{position:02d}_{token}",
                    mapping_tolerance,
                    net_charge_tolerance,
                    ase_read,
                )
                deliverable["pacman"] = pacman_result
                deliverable["curation_stage_eligible"] = True
            record["deliverables"].append(deliverable)

        if record["review_reasons"]:
            record["status"] = "review"
        else:
            record["status"] = "complete"
        record["curation_stage_eligible"] = (
            record["status"] == "complete"
            and not skip_pacman
            and all(
                item.get("curation_stage_eligible", False)
                for item in record["deliverables"]
            )
        )
        record["release_eligible"] = False
        record["release_gate_reason"] = (
            "dimensionality, porosity, chemistry checkers, manual regression, "
            "and v12 merge policy remain deferred"
        )
    except Exception as exc:
        record["status"] = "failed"
        record["errors"].append({
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        })
    record["finished_at"] = utcnow()
    write_json_atomic(done_path, record)
    return record


def merge_records(output: Path, report_dir: Path, label: str) -> None:
    records_dir = output / "state" / "records"
    records: list[dict[str, Any]] = []
    for path in sorted(records_dir.glob("*.json")):
        with path.open(encoding="utf-8") as handle:
            records.append(json.load(handle))
    report_dir.mkdir(parents=True, exist_ok=True)
    jsonl_path = report_dir / f"{label}_curation_records.jsonl"
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=report_dir,
        prefix=f".{jsonl_path.name}.", delete=False,
    ) as out:
        for record in records:
            out.write(json.dumps(record, sort_keys=True, default=str) + "\n")
        temporary = Path(out.name)
    os.replace(temporary, jsonl_path)

    versions = sorted({record.get("pipeline_version", "unknown") for record in records})
    categories: dict[str, int] = {}
    for record in records:
        category = str(record.get("category") or "unknown")
        categories[category] = categories.get(category, 0) + 1
    summary = {
        "pipeline_version": versions[0] if len(versions) == 1 else "mixed",
        "record_pipeline_versions": versions,
        "created_at": utcnow(),
        "records": len(records),
        "complete": sum(record.get("status") == "complete" for record in records),
        "review": sum(record.get("status") == "review" for record in records),
        "failed": sum(record.get("status") == "failed" for record in records),
        "curation_stage_eligible": sum(
            bool(record.get("curation_stage_eligible")) for record in records
        ),
        "deliverables": sum(len(record.get("deliverables") or []) for record in records),
        "categories": categories,
        "records_jsonl": str(jsonl_path),
    }
    write_json_atomic(report_dir / f"{label}_curation_summary.json", summary)

    csv_path = report_dir / f"{label}_curation_summary.csv"
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=report_dir,
        prefix=f".{csv_path.name}.", delete=False,
    ) as out:
        writer = csv.DictWriter(out, fieldnames=[
            "refcode", "status", "category", "deliverable_roles",
            "uncharged_cifs", "charged_cifs", "curation_stage_eligible",
            "deposition_date", "doi", "review_reason", "error",
        ])
        writer.writeheader()
        for record in records:
            source_meta = record.get("ccdc_metadata") or {}
            bibliography = source_meta.get("bibliography") or {}
            deliverables = record.get("deliverables") or []
            errors = record.get("errors") or []
            writer.writerow({
                "refcode": record.get("refcode"),
                "status": record.get("status"),
                "category": record.get("category"),
                "deliverable_roles": ";".join(
                    "+".join(item.get("roles") or []) for item in deliverables
                ),
                "uncharged_cifs": ";".join(
                    item.get("uncharged_cif", "") for item in deliverables
                ),
                "charged_cifs": ";".join(
                    (item.get("pacman") or {}).get("charged_cif", "")
                    for item in deliverables
                ),
                "curation_stage_eligible": record.get("curation_stage_eligible"),
                "deposition_date": source_meta.get("deposition_date"),
                "doi": bibliography.get("doi") or source_meta.get("doi"),
                "review_reason": "; ".join(record.get("review_reasons") or []),
                "error": errors[0].get("message") if errors else "",
            })
        temporary = Path(out.name)
    os.replace(temporary, csv_path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, help="recursive source-CIF directory")
    parser.add_argument(
        "--source-cif", type=Path,
        help="process exactly this CIF (recommended for fresh-process execution)",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report-dir", type=Path, required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--metadata", type=Path, help="CSD-export metadata JSONL")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--refcode", action="append")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--skip-pacman", action="store_true",
        help="generate and validate solvent-removal candidates without charges",
    )
    parser.add_argument(
        "--initial-total-margin", type=float,
        default=INITIAL_TOTAL_PAIR_MARGIN_ANGSTROM,
        help="requested total pair margin in A (ASE receives half per atom)",
    )
    parser.add_argument("--max-adaptive-iterations", type=int, default=40)
    parser.add_argument("--mapping-tolerance", type=float, default=1.0e-3)
    parser.add_argument("--net-charge-tolerance", type=float, default=1.0e-5)
    parser.add_argument("--merge-only", action="store_true")
    parser.add_argument("--no-merge", action="store_true")
    args = parser.parse_args()

    if args.merge_only:
        merge_records(args.output, args.report_dir, args.label)
        return 0
    if args.input is None and args.source_cif is None:
        raise SystemExit("--input or --source-cif is required unless --merge-only is used")
    if args.initial_total_margin < 0:
        raise SystemExit("--initial-total-margin must be non-negative")
    if args.max_adaptive_iterations < 1:
        raise SystemExit("--max-adaptive-iterations must be positive")
    if args.mapping_tolerance <= 0 or args.net_charge_tolerance < 0:
        raise SystemExit("mapping tolerance must be positive and charge tolerance non-negative")

    source_metadata = load_csd_metadata(args.metadata)
    inputs = [args.source_cif] if args.source_cif is not None else cif_inputs(args.input)
    if args.shard_count < 1 or not 0 <= args.shard_index < args.shard_count:
        raise SystemExit("shard-index must be in [0, shard-count)")
    if args.refcode:
        wanted = set(args.refcode)
        inputs = [path for path in inputs if path.stem in wanted]
    inputs = inputs[args.shard_index::args.shard_count]
    if args.limit is not None:
        inputs = inputs[:args.limit]
    if not inputs:
        raise SystemExit("No input CIFs selected")
    if len(inputs) != 1:
        raise SystemExit(
            f"Selected {len(inputs)} CIFs, but v4 requires one fresh Python "
            "process per CIF. Use --source-cif from a launcher or Slurm array."
        )

    args.output.mkdir(parents=True, exist_ok=True)
    source = inputs[0]
    result = run_one(
        source,
        args.output,
        source_metadata,
        overwrite=args.overwrite,
        skip_pacman=args.skip_pacman,
        initial_total_margin=args.initial_total_margin,
        max_adaptive_iterations=args.max_adaptive_iterations,
        mapping_tolerance=args.mapping_tolerance,
        net_charge_tolerance=args.net_charge_tolerance,
    )
    print(f"{result['refcode']}\t{result['status']}\t{result.get('category')}")
    if not args.no_merge:
        merge_records(args.output, args.report_dir, args.label)
    return 0 if result["status"] in {"complete", "review"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
