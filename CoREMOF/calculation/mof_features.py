"""Analysis topology, open metal sites, revised autocorrelation and so on.
"""

import math
import numbers
import os
import tempfile

from ase.io import read
from pymatgen.core.structure import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def SpaceGroup(structure):

    """Analysis space group of structure.

    Args:
        structure (str): path to your CIF.
       
    Returns:
        Dictionary:
            -   unit by ["unit"], always nan
            -   hall symbol by ["hall_symbol"]
            -   space group number by ["space_group_number"]
            -   crystal system by ["crystal_system"]
    """
           
    result_sg = {}
    result_sg["unit"]="nan"
    atoms = read(structure)
    structure_ = AseAtomsAdaptor.get_structure(atoms)
    result_ = SpacegroupAnalyzer(structure_, symprec=0.01, angle_tolerance=5)
    hall_symbol = result_.get_hall()
    space_group_number = result_.get_space_group_number()
    crystal_system = result_.get_crystal_system()
    result_sg["hall_symbol"]=hall_symbol
    result_sg["space_group_number"]=space_group_number
    result_sg["crystal_system"]=crystal_system

    return result_sg

def Mass(structure):

    """Analysis total mass of structure.

    Args:
        structure (str): path to your CIF.
       
    Returns:
        Dictionary:
            -   unit by ["unit"], always amu
            -   total mass by ["total_mass"]
    """

    result_m = {}
    result_m["unit"]="amu"
    atoms = read(structure)
    total_mass = atoms.get_masses().sum()
    result_m["total_mass"]=float(total_mass)

    return result_m

def Volume(structure):

    """Analysis total volume of structure.

    Args:
        structure (str): path to your CIF.
       
    Returns:
        Dictionary:
            -   unit by ["unit"], always Å^3
            -   total volume by ["total_volume"]
    """

    result_v = {}
    result_v["unit"]="Å^3"
    with open(structure, "r", encoding="utf-8") as f:
        cif_file = f.read()
    atoms = Structure.from_str(cif_file, fmt="cif")
    total_volume = atoms.volume
    result_v["total_volume"]=total_volume

    return result_v

def n_atom(structure):

    """Analysis number of atoms of structure.

    Args:
        structure (str): path to your CIF.
       
    Returns:
        Dictionary:
            -   unit by ["unit"], always nan
            -   number of atoms by ["number_atoms"]
    """

    result_na = {}
    result_na["unit"]="nan"
    atoms = read(structure)
    number_atoms = len(atoms)
    result_na["number_atoms"]=number_atoms

    return result_na

def topology(structure, node_type="single"):

    """Analysis topology of structure by CrystalNets.jl (https://github.com/coudertlab/CrystalNets.jl?tab=readme-ov-file).

    Args:
        structure (str): path to your CIF.
        node_type (str): the clustering algorithm used to group atoms into vertices. single: each already-defined cluster is mapped to a vertex; all: keep points of extension for organic clusters.
       
    Returns:
        Dictionary:
            -   dimension by ["dimension"]
            -   topology by ["topology"]
            -   catenation by ["catenation"]
    """

    try:
        import juliacall
    except ImportError as exc:
        raise ImportError(
            "Topology analysis requires juliacall and CrystalNets.jl. "
            "Install juliacall, then retry while connected to the internet."
        ) from exc

    jl = juliacall.newmodule("topo")
    jl.seval("using CrystalNets")

    if node_type == "single":
        clustering = jl.Clustering.SingleNodes
    elif node_type == "all":
        clustering = jl.Clustering.AllNodes
    else:
        raise ValueError("node_type should be single or all")
    options = jl.CrystalNets.Options(structure=jl.StructureType.MOF, clusterings=[clustering])
    result = jl.determine_topology(structure, options)
    
    result_tp = {}
    result_tp["dimension"] = []
    result_tp["topology"] = []
    result_tp["catenation"] = []
    interpenetration = jl.CrystalNets.total_interpenetration(result, clustering)
    for x in result:
        info = x[0][clustering]
        result_tp["dimension"].append(jl.ndims(info.genome))
        result_tp["topology"].append(str(info))
        result_tp["catenation"].append(interpenetration[info])

    return result_tp


def get_oms_file(structure):

    """Analysis open metal site of structure from CoRE MOF 2019 (https://github.com/emmhald/open_metal_detector).

    Args:
        structure (str): path to your CIF.
       
    Returns:
        Dictionary:
            -   all types of metal by ["Metal Types"]
            -    has OMS or not by ["Has OMS"], -> True or False
            -    of type of OMS if has by ["OMS Types"]
    """
        
    from CoREMOF.calculation.mof_collection import MofCollection

    with tempfile.TemporaryDirectory(prefix="coremof_oms_") as analysis_folder:
        a_mof_collection = MofCollection(
            path_list=[structure], analysis_folder=analysis_folder
        )
        a_mof_collection.analyse_mofs(num_batches=1, overwrite=False)
        name = os.path.splitext(os.path.basename(structure))[0]
        oms_result = {
            "Metal Types": a_mof_collection.mof_oms_df["Metal Types"][name],
            "Has OMS": a_mof_collection.mof_oms_df["Has OMS"][name],
            "OMS Types": a_mof_collection.mof_oms_df["OMS Types"][name],
        }

    return oms_result


def get_oms_folder(input_folder, n_batch = 1):

    """Analysis open metal site of folder with structures from CoRE MOF 2019 (https://github.com/emmhald/open_metal_detector).

    Args:
        input_folder (str): path to your folder.
        n_batch (int): number of batches.
       
    Returns:
        Dictionary:
            -   all types of metal of each structure by [structure]["Metal Types"]
            -   has OMS or not of each structure by [structure]["Has OMS"], -> True or False
            -   type of OMS if has of each structure by [structure]["OMS Types"]
    """

    from CoREMOF.calculation.mof_collection import MofCollection

    with tempfile.TemporaryDirectory(prefix="coremof_oms_") as analysis_folder:
        mof_collection = MofCollection.from_folder(
            collection_folder=input_folder, analysis_folder=analysis_folder
        )
        mof_collection.analyse_mofs(num_batches=n_batch, overwrite=False)
        oms_result = {}
        for name, row in mof_collection.mof_oms_df.iterrows():
            oms_result[name] = {
                "Metal Types": row.iloc[0],
                "Has OMS": row.iloc[1],
                "OMS Types": row.iloc[2],
            }

    return oms_result


def _rac_family_names(prefix, properties, maximum_depth, suffix="-all"):
    return [
        f"{prefix}-{property_name}-{current_depth}{suffix}"
        for property_name in properties
        for current_depth in range(maximum_depth + 1)
    ]


def _rac_names_by_group(maximum_depth):
    common_properties = ("I", "S", "T", "Z", "chi")
    linker_properties = ("I", "S", "T", "Z", "alpha", "chi")
    return {
        "Metal": (
            _rac_family_names("D_mc", common_properties, maximum_depth)
            + _rac_family_names("f", common_properties, maximum_depth)
            + _rac_family_names("mc", common_properties, maximum_depth)
        ),
        "Linker": (
            _rac_family_names("D_lc", linker_properties, maximum_depth)
            + _rac_family_names("lc", linker_properties, maximum_depth)
            + _rac_family_names(
                "f-lig", common_properties, maximum_depth, suffix=""
            )
        ),
        "Function-group": (
            _rac_family_names("D_func", linker_properties, maximum_depth)
            + _rac_family_names("func", linker_properties, maximum_depth)
        ),
    }


def RACs(structure, depth=3):

    """Revised Autocorrelation features (https://github.com/hjkgrp/molSimplify).

    Args:
        structure (str or path-like): path to one CIF.
        depth (int): maximum autocorrelation depth. The historical public
            default is 3. Setting this to 5 changes the descriptor depth but
            does not by itself reproduce the sealed CoRE-MOF v26 RAC5 method.
       
    Returns:
        Dictionary:
            -   metal by ["Metal"]
            -   linker by ["Linker"]
            -   function group by ["Function-group"]
    """

    if isinstance(depth, bool) or not isinstance(depth, numbers.Integral):
        raise TypeError("depth must be a non-bool integer")
    if depth < 0:
        raise ValueError("depth must be greater than or equal to zero")
    depth = int(depth)
    expected_names = _rac_names_by_group(depth)

    result_rac = {}
    result_rac["Metal"] = {}
    result_rac["Linker"] = {}
    result_rac["Function-group"] = {}

    try:
        from molSimplify.Informatics.MOF.MOF_descriptors import get_MOF_descriptors
    except ImportError as exc:
        raise ImportError(
            "RAC descriptors require molSimplify. Install it with "
            "'pip install molSimplify'."
        ) from exc

    name = os.path.splitext(os.path.basename(structure))[0]
    with tempfile.TemporaryDirectory(prefix="coremof_rac_") as workdir:
        full_names, full_descriptors = get_MOF_descriptors(
            structure,
            depth,
            path=workdir,
            xyzpath=os.path.join(workdir, f"{name}.xyz"),
            max_num_atoms=6000,
        )
                                                        
    names = list(full_names)
    descriptors = list(full_descriptors)
    if len(names) != len(descriptors):
        raise ValueError("molSimplify returned different RAC name/value counts")
    if len(names) != len(set(names)):
        raise ValueError("molSimplify returned duplicate RAC names")
    descriptor_data = dict(zip(names, descriptors))

    for group, group_names in expected_names.items():
        for descriptor_name in group_names:
            if descriptor_name not in descriptor_data:
                raise ValueError(
                    "molSimplify RAC output is missing expected descriptor "
                    f"{descriptor_name!r} at depth {depth}"
                )
            try:
                value = float(descriptor_data[descriptor_name])
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"molSimplify RAC descriptor {descriptor_name!r} is not numeric"
                ) from exc
            if not math.isfinite(value):
                raise ValueError(
                    f"molSimplify RAC descriptor {descriptor_name!r} is not finite"
                )
            result_rac[group][descriptor_name] = round(value, 4)

    return result_rac
