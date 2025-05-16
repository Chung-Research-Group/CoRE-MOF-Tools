"""Analysis topology, open metal sites, revised autocorrelation and so on.
"""

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import read
from pymatgen.core.structure import Structure

import os,juliacall
from CoREMOF.calculation.mof_collection import MofCollection
from CoREMOF.utils import remove

from molSimplify.Informatics.MOF.MOF_descriptors import get_MOF_descriptors


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
    """

    package_directory = os.path.abspath(__file__).replace("mof_features.py","")
    os.environ["JULIA_DEPOT_PATH"] = package_directory

    juliacall.Main.seval('import Pkg; Pkg.add("CrystalNets")')

    jl = juliacall.newmodule("topo")
    jl.seval("using CrystalNets")

    options = jl.CrystalNets.Options(structure=jl.StructureType.MOF)
    result = jl.determine_topology(structure, options)
    result_tp = {}
    result_save = []
    for x in result:
        if node_type=="single":
            result_ = x[0][jl.Clustering.SingleNodes]
            result_save.append(result_)
        elif node_type=="all":
            result_ = x[0][jl.Clustering.AllNodes]
            result_save.append(result_)

    result_tp["dimension"]=[]
    result_tp["topology"] = []
    for info in result_save:
        
        result_tp["dimension"].append(jl.ndims(info.genome))
        result_tp["topology"].append(str(info))

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
        
    a_mof_collection = MofCollection(path_list = [structure], 
                                 analysis_folder="tmp_oms")
    a_mof_collection.analyse_mofs(num_batches=1,overwrite=False)
    oms_result = {
                    "Metal Types": a_mof_collection.mof_oms_df["Metal Types"][structure.replace(".cif","").split("/")[-1]],
                    "Has OMS": a_mof_collection.mof_oms_df["Has OMS"][structure.replace(".cif","").split("/")[-1]],
                    "OMS Types": a_mof_collection.mof_oms_df["OMS Types"][structure.replace(".cif","").split("/")[-1]],
                }

    remove.remove_dir_with_permissions("tmp_oms")

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

    mof_collection = MofCollection.from_folder(collection_folder = input_folder, 
                                                analysis_folder="tmp_oms")
    mof_collection.analyse_mofs(num_batches=n_batch,overwrite=False)
    oms_result = {}
    for name, row in mof_collection.mof_oms_df.iterrows():
        
        oms_result[name] = {
                            "Metal Types": row[0],
                            "Has OMS": row[1],
                            "OMS Types": row[2]
                        }

    remove.remove_dir_with_permissions("tmp_oms")

    return oms_result


def RACs(structure):

    """Revised Autocorrelation features (https://github.com/hjkgrp/molSimplify).

    Args:
        input_folder (str): path to your folder.
        n_batch (int): number of batches.
       
    Returns:
        Dictionary:
            -   metal by ["Metal"]
            -   linker by ["Linker"]
            -   function group by ["Function-group"]
    """

    metal_fnames = ['D_mc-I-0-all', 'D_mc-I-1-all', 'D_mc-I-2-all', 'D_mc-I-3-all', 'D_mc-S-0-all', 'D_mc-S-1-all', 'D_mc-S-2-all',
    'D_mc-S-3-all', 'D_mc-T-0-all', 'D_mc-T-1-all', 'D_mc-T-2-all', 'D_mc-T-3-all', 'D_mc-Z-0-all', 'D_mc-Z-1-all', 'D_mc-Z-2-all',
    'D_mc-Z-3-all', 'D_mc-chi-0-all', 'D_mc-chi-1-all', 'D_mc-chi-2-all', 'D_mc-chi-3-all', 'f-I-0-all', 'f-I-1-all', 'f-I-2-all',
    'f-I-3-all', 'f-S-0-all', 'f-S-1-all', 'f-S-2-all', 'f-S-3-all', 'f-T-0-all', 'f-T-1-all', 'f-T-2-all', 'f-T-3-all', 'f-Z-0-all',
    'f-Z-1-all', 'f-Z-2-all', 'f-Z-3-all', 'f-chi-0-all', 'f-chi-1-all', 'f-chi-2-all', 'f-chi-3-all', 'mc-I-0-all', 'mc-I-1-all',
    'mc-I-2-all', 'mc-I-3-all', 'mc-S-0-all', 'mc-S-1-all', 'mc-S-2-all', 'mc-S-3-all', 'mc-T-0-all', 'mc-T-1-all', 'mc-T-2-all',
    'mc-T-3-all', 'mc-Z-0-all', 'mc-Z-1-all', 'mc-Z-2-all', 'mc-Z-3-all', 'mc-chi-0-all', 'mc-chi-1-all', 'mc-chi-2-all', 'mc-chi-3-all']

    linker_fnames = ['D_lc-I-0-all', 'D_lc-I-1-all', 'D_lc-I-2-all', 'D_lc-I-3-all', 'D_lc-S-0-all', 'D_lc-S-1-all', 'D_lc-S-2-all',
    'D_lc-S-3-all', 'D_lc-T-0-all', 'D_lc-T-1-all', 'D_lc-T-2-all', 'D_lc-T-3-all', 'D_lc-Z-0-all', 'D_lc-Z-1-all', 'D_lc-Z-2-all',
    'D_lc-Z-3-all', 'D_lc-alpha-0-all', 'D_lc-alpha-1-all', 'D_lc-alpha-2-all', 'D_lc-alpha-3-all', 'D_lc-chi-0-all', 'D_lc-chi-1-all',
    'D_lc-chi-2-all', 'D_lc-chi-3-all', 'lc-I-0-all', 'lc-I-1-all', 'lc-I-2-all', 'lc-I-3-all', 'lc-S-0-all', 'lc-S-1-all', 'lc-S-2-all',
    'lc-S-3-all', 'lc-T-0-all', 'lc-T-1-all', 'lc-T-2-all', 'lc-T-3-all', 'lc-Z-0-all', 'lc-Z-1-all', 'lc-Z-2-all', 'lc-Z-3-all',
    'lc-alpha-0-all', 'lc-alpha-1-all', 'lc-alpha-2-all', 'lc-alpha-3-all', 'lc-chi-0-all', 'lc-chi-1-all', 'lc-chi-2-all', 'lc-chi-3-all',
    'f-lig-I-0', 'f-lig-I-1', 'f-lig-I-2', 'f-lig-I-3', 'f-lig-S-0', 'f-lig-S-1', 'f-lig-S-2', 'f-lig-S-3', 'f-lig-T-0', 'f-lig-T-1',
    'f-lig-T-2', 'f-lig-T-3', 'f-lig-Z-0', 'f-lig-Z-1', 'f-lig-Z-2', 'f-lig-Z-3', 'f-lig-chi-0', 'f-lig-chi-1', 'f-lig-chi-2', 'f-lig-chi-3']

    fg_fnames = ['D_func-I-0-all', 'D_func-I-1-all', 'D_func-I-2-all', 'D_func-I-3-all', 'D_func-S-0-all', 'D_func-S-1-all',
    'D_func-S-2-all', 'D_func-S-3-all', 'D_func-T-0-all', 'D_func-T-1-all', 'D_func-T-2-all', 'D_func-T-3-all', 'D_func-Z-0-all',
    'D_func-Z-1-all', 'D_func-Z-2-all', 'D_func-Z-3-all', 'D_func-alpha-0-all', 'D_func-alpha-1-all', 'D_func-alpha-2-all',
    'D_func-alpha-3-all', 'D_func-chi-0-all', 'D_func-chi-1-all', 'D_func-chi-2-all', 'D_func-chi-3-all', 'func-I-0-all',
    'func-I-1-all', 'func-I-2-all', 'func-I-3-all', 'func-S-0-all', 'func-S-1-all', 'func-S-2-all', 'func-S-3-all', 'func-T-0-all',
    'func-T-1-all', 'func-T-2-all', 'func-T-3-all', 'func-Z-0-all', 'func-Z-1-all', 'func-Z-2-all', 'func-Z-3-all', 'func-alpha-0-all',
    'func-alpha-1-all', 'func-alpha-2-all', 'func-alpha-3-all', 'func-chi-0-all', 'func-chi-1-all', 'func-chi-2-all', 'func-chi-3-all']

    result_rac = {}
    result_rac["Metal"] = {}
    result_rac["Linker"] = {}
    result_rac["Function-group"] = {}

    os.makedirs("tmp_rac", exist_ok=True)
    
    name = os.path.basename(structure).replace(".cif", "")
    
    full_names, full_descriptors = get_MOF_descriptors(
                                                        structure,
                                                        3,
                                                        path='tmp_rac',
                                                        xyzpath=f'tmp_rac/{name}.xyz',
                                                        max_num_atoms=6000
                                                        )
                                                        
    descriptor_data = dict(zip(full_names, full_descriptors))

    for metal in metal_fnames:
        result_rac["Metal"][metal] =  round(float(descriptor_data[metal]), 4)
    for linker in linker_fnames:
        result_rac["Linker"][linker] =  round(float(descriptor_data[linker]), 4)
    for fg in fg_fnames:
        result_rac["Function-group"][fg] =  round(float(descriptor_data[fg]), 4)

    remove.remove_dir_with_permissions("tmp_rac")

    return result_rac
