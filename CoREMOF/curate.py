"""Process your CIF to "CoRE MOF" CIF.
"""

import collections
import csv
import functools
import itertools
import json
import math
import os
import re
import shutil
import tempfile
import warnings
from ase.io import read, write

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from CoREMOF.utils.atoms_definitions import METAL, COVALENTRADII
from CoREMOF.utils.ions_list import ALLIONS
import numpy as np
import pandas as pd

from ase.neighborlist import NeighborList
from scipy.sparse.csgraph import connected_components
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

from gemmi import cif as CIF
try:
    from PACMANCharge import pmcharge
except ImportError:
    pmcharge = None

def ensure_data(structure):
    """Precheck your CIF.

    Args:
        structure (str): path to your CIF.

    Returns:
        cif:
            -   added "data_struc" CIF
    """
        
    if not os.path.isfile(structure):
        raise FileNotFoundError(f"CIF file does not exist: {structure}")
    with open(structure, 'r', encoding='utf-8') as file:
        lines = file.readlines()
    first_content = next((i for i, line in enumerate(lines) if line.strip()), 0)
    if not lines or not lines[first_content].strip().lower().startswith('data_'):
        lines.insert(first_content, 'data_struc\n')
        with open(structure, 'w', encoding='utf-8') as file:
            file.writelines(lines)
        return True
    return False

def ase_format(structure):
    """try to read CIF and convert to ASE format.

    Args:
        structure (str): path to your CIF.

    Returns:
        cif:
            -   ASE format CIF.
    """
        
    if not os.path.isfile(structure):
        raise FileNotFoundError(f"CIF file does not exist: {structure}")
    errors = []
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            mof_temp = Structure.from_file(structure)
            mof_temp.to(filename=structure, fmt="cif")
            struc = read(structure)
            write(structure, struc)
            return structure
    except Exception as exc:
        errors.append(exc)
    try:
        struc = read(structure)
        write(structure, struc)
        return structure
    except Exception as exc:
        errors.append(exc)

    ensure_data(structure)
    try:
        struc = read(structure)
        write(structure, struc)
        return structure
    except Exception as exc:
        errors.append(exc)
        raise ValueError(
            f"Could not parse CIF {structure!r} with pymatgen or ASE: {errors[-1]}"
        ) from exc

class preprocess():

    """Precheck your CIF.

    Args:
        structure (str): path to your CIF.
        output_folder (str): the path to save processed CIF.

    Returns:
        Dictionary & cif:
            -   result of pre-check, has metal and carbon, has multi-structures.
            -   CIF by spliting, making primitive and making P1.
    """

    def __init__(self, structure, output_folder="result_curation"):
        self.structure = structure
        self.output = output_folder + os.sep
        os.makedirs(self.output, exist_ok=True)
        self.result_check=self.process()
                
    def process(self):
        result_check = self.split_pri_p1(self.structure, self.output)
        with open(self.output + os.path.basename(self.structure).replace(".cif","") + "_precheck.json", "w") as f:
            json.dump(result_check,f,indent=2)
        return result_check


    def split_pri_p1(self, structure, output_folder):
        result_check = {}
        
        structures = read(structure, index=':')
        n_struc = len(structures)
        result_check["N_structures"] = n_struc

        if n_struc > 1:
            if not os.path.exists(output_folder):
                os.makedirs(output_folder, exist_ok=True)
            print(structure, "with more than one crystal structures")
            for i, atoms in enumerate(structures):
                
                structure_= AseAtomsAdaptor.get_structure(atoms)
                sga = SpacegroupAnalyzer(structure_)
                structure_prep = sga.get_primitive_standard_structure(international_monoclinic=True, keep_site_properties=False)
                struct_name = os.path.basename(structure).replace(".cif","")+"_"+str(i+1)
                structure_prep.to(filename=os.path.join(output_folder, f"{struct_name}.cif"))

                has_metal = any(METAL.get(atom.symbol) for atom in atoms)
                has_carbon = any(atom.symbol == 'C' for atom in atoms)
                
                if has_metal:
                    if has_carbon:
                        result_check[struct_name] = "has metal and carbon"
                    else:
                        result_check[struct_name] = "missing carbon"
                else:
                    if has_carbon:
                        result_check[struct_name] = "missing metal"
                    else:
                        result_check[struct_name] = "missing metal and carbon"
        else:
            atoms = structures[0]
            struct_name = os.path.basename(structure).replace(".cif","")

            has_metal = any(METAL.get(atom.symbol) for atom in atoms)
            has_carbon = any(atom.symbol == 'C' for atom in atoms)
            
            if has_metal:
                if has_carbon:
                    result_check[struct_name] = "has metal and carbon"
                else:
                    result_check[struct_name] = "missing carbon"
            else:
                if has_carbon:
                    result_check[struct_name] = "missing metal"
                else:
                    result_check[struct_name] = "missing metal and carbon"

            structure_= AseAtomsAdaptor.get_structure(atoms)
            sga = SpacegroupAnalyzer(structure_)
            structure_prep = sga.get_primitive_standard_structure(international_monoclinic=True, keep_site_properties=False)
            temp_path = os.path.join(output_folder, f"{struct_name}.cif")
            structure_prep.to(filename=temp_path)
            
        return result_check
    

class clean():

    """Removing free solvent and coordinated solvent but keep ions based on a list.         
    
    Args:
        structure (str): path to your CIF.
        initial_skin (float): ASE skin in angstrom, added to each atom's
            covalent radius. The total pair margin is twice this value.
            The legacy default of 0.25 therefore gives a 0.50 angstrom
            pair margin. Recorded release-curation methods have separate
            settings and do not change this default.
        output_folder (str): the path to save processed CIF.
        saveto (bool or str): the name of csv file with clean result.

    Returns:
        CSV or cif:
            -   result of curating (name, skin, removed solvent).
            -   CIF by curating.
    """

    def __init__(self, structure, initial_skin = 0.25, output_folder="result_curation", saveto: str="clean_result.csv")-> pd.DataFrame:
        
        self.cambridge_radii = COVALENTRADII
        self.metal_list = [element for element, is_metal in METAL.items() if is_metal]
        self.ions_list = set(ALLIONS)
        
        self.structure = structure
        self.initial_skin = initial_skin
        self.output = output_folder
        os.makedirs(self.output, exist_ok=True)

        if saveto:
            self.csv_path = os.path.join(self.output, saveto)
        self.process()

    def process(self):

        """start to run curation.
        """
            
        try:
            with open(self.csv_path, mode='w', newline='') as csv_file:
                csv_writer = csv.writer(csv_file)
                csv_writer.writerow(["Name", "Skin_FSR", "Skin_ASR", "Removed_FSR", "Removed_ASR"])

                fsr_skin, fsr_results = self.run_fsr(self.structure, self.output, self.initial_skin, self.metal_list, self.ions_list)
                print(f"FSR results for {self.structure}: {fsr_results}")
                
                asr_skin, asr_results = self.run_asr(self.structure, self.output, self.initial_skin, self.metal_list, self.ions_list)
                print(f"ASR results for {self.structure}: {asr_results}")

                csv_writer.writerow([
                                        os.path.basename(self.structure),
                                        str(fsr_skin) if fsr_skin else "none",
                                        str(asr_skin) if asr_skin else "none",
                                        str(fsr_results) if fsr_results else "none",
                                        str(asr_results) if asr_results else "none"
                                    ])
        except:
            fsr_skin, fsr_results = self.run_fsr(self.structure, self.output, self.initial_skin, self.metal_list, self.ions_list)
            print(f"FSR results for {self.structure}: {fsr_results}")
            
            asr_skin, asr_results = self.run_asr(self.structure, self.output, self.initial_skin, self.metal_list, self.ions_list)
            print(f"ASR results for {self.structure}: {asr_results}")

    def run_fsr(self, mof, save_folder, initial_skin, metal_list, ions_list):

        """free solvent function.
        """

        m = mof.replace(".cif", "")
        skin = initial_skin
        try:
            while True:
                skin, printed_formulas = self.free_clean(m, save_folder, ions_list, skin)
                has_metals = False
                for e_s in printed_formulas:
                    split_formula = re.findall(r'([A-Z][a-z]?)(\d*)', e_s)
                    elements = [match[0] for match in split_formula]
                    if any(e in metal_list for e in elements):
                        has_metals = True
                        skin += 0.05
                if not has_metals:
                    break
            return skin, printed_formulas
        except Exception as e:
            print(m, str(e))


    def run_asr(self, mof, save_folder, initial_skin, metal_list, ions_list):

        """all solvent function.
        """

        m = mof.replace(".cif", "")
        skin = initial_skin
        try:
            while True:
                skin, printed_formulas = self.all_clean(m, save_folder, ions_list, skin)
                has_metals = False
                for e_s in printed_formulas:
                    split_formula = re.findall(r'([A-Z][a-z]?)(\d*)', e_s)
                    elements = [match[0] for match in split_formula]
                    if any(e in metal_list for e in elements):
                        has_metals = True
                        skin += 0.05
                if not has_metals:
                    break
            return skin, printed_formulas
        except Exception as e:
            print(f"{m} Fail: {e}")

    def build_ASE_neighborlist(self, cif, skin):

        """get list of neighbor.
        """

        radii = [self.cambridge_radii[i] for i in cif.get_chemical_symbols()]
        ASE_neighborlist = NeighborList(radii, self_interaction=False, bothways=True, skin=skin)
        ASE_neighborlist.update(cif)
        return ASE_neighborlist

    def find_clusters(self, adjacency_matrix, atom_count):

        """define cluster by the connected components of a sparse graph.
        """

        clusters = []
        cluster_count, clusterIDs = connected_components(adjacency_matrix, directed=True)
        for n in range(cluster_count):
            clusters.append([i for i in range(atom_count) if clusterIDs[i] == n])
        return clusters

    def find_metal_connected_atoms(self, structure, neighborlist):

        """get the atom connected with metal atom.
        """

        metal_connected_atoms = []
        metal_atoms = []
        for i, elem in enumerate(structure.get_chemical_symbols()):
            if elem in self.metal_list:
                neighbors, _ = neighborlist.get_neighbors(i)
                metal_connected_atoms.append(neighbors)
                metal_atoms.append(i)
        return metal_connected_atoms, metal_atoms, structure

    def CustomMatrix(self, neighborlist, atom_count):

        """convert to matrix.
        """

        matrix = np.zeros((atom_count, atom_count), dtype=int)
        for i in range(atom_count):
            neighbors, _ = neighborlist.get_neighbors(i)
            for j in neighbors:
                matrix[i][j] = 1
        return matrix

    def mod_adjacency_matrix(self, adj_matrix, MetalConAtoms, MetalAtoms, atom_count, struct):

        """modify matrix by breakdown the bond that atom connect with metal atom.
        """

        clusters = self.find_clusters(adj_matrix, atom_count)
        for i, element_1 in enumerate(MetalAtoms):
            for j, element_2 in enumerate(MetalConAtoms[i]):
                if struct[element_2].symbol == "O":
                    tmp = len(self.find_clusters(adj_matrix, atom_count))
                    adj_matrix[element_2][element_1] = 0
                    adj_matrix[element_1][element_2] = 0
                    new_clusters = self.find_clusters(adj_matrix, atom_count)
                    if tmp == len(new_clusters):
                        adj_matrix[element_2][element_1] = 1
                        adj_matrix[element_1][element_2] = 1
                    for ligand in new_clusters:
                        if ligand not in clusters:
                            tmp3 = struct[ligand].get_chemical_symbols()
                            if "O" and "H" in tmp3 and len(tmp3) == 2:
                                adj_matrix[element_2][element_1] = 1
                                adj_matrix[element_1][element_2] = 1
        return adj_matrix

    def cmp(self, x, y):
        return (x > y) - (x < y)

    def cluster_to_formula(self, cluster, cif):

        """convert to chemical formula.
        """

        symbols = [cif[i].symbol for i in cluster]
        count = collections.Counter(symbols)
        formula = ''.join([atom + (str(count[atom]) if count[atom] > 1 else '') for atom in sorted(count)])
        return formula

    def free_clean(self, input_file, save_folder, ions_list, skin):

        """workflow of removing free solvent.
        """

        try:
            print(input_file+".cif")
            cif = read(input_file+".cif")
            refcode = input_file.split("/")[-1]
            atom_count = len(cif.get_chemical_symbols())
            ASE_neighborlist = self.build_ASE_neighborlist(cif,skin)
            a = self.CustomMatrix(ASE_neighborlist,atom_count)
            b = self.find_clusters(a,atom_count)
            b.sort(key=functools.cmp_to_key(lambda x,y: self.cmp(len(x), len(y))))
            b.reverse()
            cluster_length=[]
            solvated_cluster = []
            ions_cluster = []
            printed_formulas = []
            iii=False
            print("using",skin,"as skin")
            for index, _ in enumerate(b):
                cluster_formula = self.cluster_to_formula(b[index], cif) 
                if cluster_formula in ions_list:
                    print(cluster_formula, "is ion")
                    ions_cluster.append(b[index])
                    iii=True
                else:
                    tmp = len(b[index])
                    if len(cluster_length) > 0:
                        if tmp > max(cluster_length):
                            cluster_length = []
                            solvated_cluster = []
                            solvated_cluster.append(b[index])
                            cluster_length.append(tmp)
                        elif tmp > 0.5 * max(cluster_length):
                            solvated_cluster.append(b[index])
                            cluster_length.append(tmp)
                        else:
                            formula = self.cluster_to_formula(b[index], cif)
                            if formula not in printed_formulas:
                                printed_formulas.append(formula)
                    else:
                        solvated_cluster.append(b[index])
                        cluster_length.append(tmp)
            solvated_cluster = solvated_cluster + ions_cluster
            solvated_merged = list(itertools.chain.from_iterable(solvated_cluster))
            atom_count = len(cif[solvated_merged].get_chemical_symbols())

            if iii:
                refcode=refcode.replace("_FSR","")
                new_fn = refcode + '_ION_FSR.cif'
            else:
                refcode=refcode.replace("_FSR","")
                new_fn = refcode + '_FSR.cif'
            write(os.path.join(save_folder, new_fn), cif[solvated_merged])
        
            return skin, printed_formulas
        except Exception as e:
            print(f"{input_file} Fail: {e}")

    def all_clean(self, input_file, save_folder, ions_list, skin):

        """workflow of removing all solvent.
        """

        try:
            fn = input_file
            cif = read(input_file+".cif")
            refcode = fn.split("/")[-1]
            atom_count = len(cif.get_chemical_symbols())
            ASE_neighborlist = self.build_ASE_neighborlist(cif,skin)
            a = self.CustomMatrix(ASE_neighborlist,atom_count)
            b = self.find_clusters(a,atom_count)
            b.sort(key=functools.cmp_to_key(lambda x,y: self.cmp(len(x), len(y))))
            b.reverse()
            cluster_length=[]
            solvated_cluster = []

            printed_formulas = []

            ions_cluster = []
            iii=False
            print("using",skin,"as skin")

            for index, _ in enumerate(b):
                
                cluster_formula = self.cluster_to_formula(b[index], cif) 
                if cluster_formula in ions_list:
                    print(cluster_formula, "is ion")
                    ions_cluster.append(b[index])
                    solvated_cluster.append(b[index])
                    iii=True

                else:
                    tmp = len(b[index])
                    if len(cluster_length) > 0:
                        if tmp > max(cluster_length):
                            cluster_length = []
                            solvated_cluster = []
                            solvated_cluster.append(b[index])
                            cluster_length.append(tmp)
                        if tmp > 0.5 * max(cluster_length):
                            solvated_cluster.append(b[index])
                            cluster_length.append(tmp)
                        else:
                            formula = self.cluster_to_formula(b[index], cif)
                            if formula not in printed_formulas:
                                printed_formulas.append(formula)
                    else:
                        solvated_cluster.append(b[index])
                        cluster_length.append(tmp)
                    
            solvated_merged = list(itertools.chain.from_iterable(solvated_cluster))
            
            atom_count = len(cif[solvated_merged].get_chemical_symbols())
            
            newASE_neighborlist = self.build_ASE_neighborlist(cif[solvated_merged],skin)
            MetalCon, MetalAtoms, struct = self.find_metal_connected_atoms(cif[solvated_merged], newASE_neighborlist)
            c = self.CustomMatrix(newASE_neighborlist,atom_count)
            d = self.mod_adjacency_matrix(c, MetalCon, MetalAtoms,atom_count,struct)
            solvated_clusters2 = self.find_clusters(d,atom_count)
            solvated_clusters2.sort(key=functools.cmp_to_key(lambda x,y: self.cmp(len(x), len(y))))
            solvated_clusters2.reverse()
            cluster_length=[]
            final_clusters = []

            ions_cluster2 = []
            for index, _ in enumerate(solvated_clusters2):

                cluster_formula2 = self.cluster_to_formula(solvated_clusters2[index], struct) 
                if cluster_formula2 in ions_list:
                    final_clusters.append(solvated_clusters2[index])
                    iii=True
                else:
                    tmp = len(solvated_clusters2[index])
                    if len(cluster_length) > 0:
                        if tmp > max(cluster_length):
                            cluster_length = []
                            final_clusters = []
                            final_clusters.append(solvated_clusters2[index])
                            cluster_length.append(tmp)
                        if tmp > 0.5 * max(cluster_length):
                            final_clusters.append(solvated_clusters2[index])
                            cluster_length.append(tmp)
                        else:
                            formula = self.cluster_to_formula(solvated_clusters2[index], struct)
                            if formula not in printed_formulas:
                                printed_formulas.append(formula)
                    else:
                        final_clusters.append(solvated_clusters2[index])
                        cluster_length.append(tmp)
            if iii:
                final_clusters = final_clusters+ions_cluster2
            else:
                final_clusters = final_clusters
            final_merged = list(itertools.chain.from_iterable(final_clusters))
            tmp = struct[final_merged].get_chemical_symbols()
            tmp.sort()
            if iii:
                new_fn = refcode + '_ION_ASR.cif'
            else:
                new_fn = refcode + "_ASR.cif"
            write(os.path.join(save_folder, new_fn), struct[final_merged])
            return skin, printed_formulas
    
        except Exception as e:
            print(f"{input_file} Fail: {e}")

class mof_check:
    """Retired checker execution. Use the release's precomputed results."""

    def __init__(self, *args, **kwargs):
        from ._checker_execution import results_only
        results_only("Chen-Manz and MOFChecker")

    def check(self, *args, **kwargs):
        from ._checker_execution import results_only
        return results_only("Chen-Manz and MOFChecker")

    def Chen_Manz(self, *args, **kwargs):
        from ._checker_execution import results_only
        return results_only("Chen-Manz")

    def mof_checker(self, *args, **kwargs):
        from ._checker_execution import results_only
        return results_only("MOFChecker")




def _validated_pacman_charges(source_path, charged_path, atoms):
    """Require complete charges on unchanged, ordered, fully occupied CIF sites.

    This deliberately rejects symmetry expansion, disorder and reordered rows;
    it does not guess a bijection between a CIF atom loop and ASE atoms.
    """
    source = CIF.read_file(os.fspath(source_path)).sole_block()
    charged = CIF.read_file(os.fspath(charged_path)).sole_block()
    atom_count = len(atoms)
    atom_tags = {
        "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x",
        "_atom_site_fract_y", "_atom_site_fract_z",
    }
    for block, anchor, required in (
        (source, "_atom_site_label", atom_tags),
        (charged, "_atom_site_charge", atom_tags | {"_atom_site_charge"}),
    ):
        loop = block.find_loop(anchor).get_loop()
        if loop is None or not required.issubset({tag.lower() for tag in loop.tags}):
            raise ValueError("PACMAN charges and atom identities must share one atom-site loop")

    def column(block, tag):
        values = list(block.find_loop(tag))
        if len(values) != atom_count or not values:
            raise ValueError(
                f"PACMAN atom mapping requires {atom_count} values for {tag}"
            )
        return values

    def numbers(values, label):
        parsed = [float(CIF.as_number(value)) for value in values]
        if not all(math.isfinite(value) for value in parsed):
            raise ValueError(f"PACMAN {label} contains missing or non-finite values")
        return parsed

    charges = numbers(column(charged, "_atom_site_charge"), "charge vector")
    labels = column(source, "_atom_site_label")
    if len(set(labels)) != atom_count or labels != column(charged, "_atom_site_label"):
        raise ValueError("PACMAN atom labels are duplicated or reordered")
    symbols = column(source, "_atom_site_type_symbol")
    if symbols != column(charged, "_atom_site_type_symbol"):
        raise ValueError("PACMAN changed atom types or their order")
    if symbols != atoms.get_chemical_symbols():
        raise ValueError("CIF atom rows do not match the ordered ASE atoms")

    coordinates = []
    for axis in "xyz":
        tag = "_atom_site_fract_" + axis
        original = numbers(column(source, tag), tag)
        if original != numbers(column(charged, tag), tag):
            raise ValueError("PACMAN changed fractional coordinates or atom order")
        coordinates.append(original)
    # Only numerical roundoff in ASE's Cartesian/fractional conversion is
    # tolerated here. Source/derived CIF coordinates above must be equal.
    parsed_coordinates = list(atoms.get_scaled_positions(wrap=False))
    if len(parsed_coordinates) != atom_count or any(len(row) != 3 for row in parsed_coordinates):
        raise ValueError("ASE did not return one coordinate triplet per atom")
    if any(
        not math.isclose(original, float(parsed), rel_tol=0, abs_tol=1e-12)
        for original_row, parsed_row in zip(zip(*coordinates), parsed_coordinates)
        for original, parsed in zip(original_row, parsed_row)
    ):
        raise ValueError("CIF atom rows require an unsupported ASE atom mapping")

    for block in (source, charged):
        occupancy = list(block.find_loop("_atom_site_occupancy"))
        if occupancy and (
            len(occupancy) != atom_count
            or any(value != 1.0 for value in numbers(occupancy, "occupancies"))
        ):
            raise ValueError("PACMAN curation requires complete unit occupancies")
    for tag in (
        "_cell_length_a", "_cell_length_b", "_cell_length_c",
        "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma",
    ):
        if numbers([source.find_value(tag)], tag) != numbers([charged.find_value(tag)], tag):
            raise ValueError("PACMAN changed the unit cell")
    for tags in (
        ("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz"),
        ("_space_group_name_H-M_alt", "_symmetry_space_group_name_H-M"),
        ("_space_group_IT_number", "_symmetry_Int_Tables_number"),
    ):
        def symmetry_values(block):
            for tag in tags:
                values = tuple(block.find_values(tag))
                if values:
                    return values
            return ()
        if symmetry_values(source) != symmetry_values(charged):
            raise ValueError("PACMAN changed the declared symmetry")
    return charges


class clean_pacman():

    """Remove free solvent while retaining PACMAN-predicted charged components.

    Only FSR is implemented. Coordinated-solvent removal (ASR) is unavailable;
    no ASR-labelled CIF is produced. This generic heuristic is not the sealed
    CoRE-MOF release curation method.
    
    Args:
        structure (str): path to your CIF.
        initial_skin (float): ASE skin in angstrom, added to each atom's
            covalent radius. The total pair margin is twice this value.
            The legacy default of 0.25 therefore gives a 0.50 angstrom
            pair margin. Recorded release-curation methods have separate
            settings and do not change this default.
        output_folder (str): the path to save processed CIF.
        saveto (bool or str): the name of csv file with clean result.

    Returns:
        CSV or cif:
            -   result of curating (name, skin, removed solvent, charge of solvent).
            -   CIF by curating.
    """
    def __init__(self, structure, initial_skin=0.25, output_folder="result_curation", saveto: str="clean_result.csv") -> pd.DataFrame:
        self.structure = structure
        self.initial_skin = initial_skin
        self.output = output_folder
        self.saveto = saveto
        self.cambridge_radii = COVALENTRADII
        self.metal_list = [element for element, is_metal in METAL.items() if is_metal]

        os.makedirs(self.output, exist_ok=True)
        if self.saveto:
            self.csv_path = os.path.join(self.output, self.saveto)

        self.run_pacman()
        self.process()
        
    def run_pacman(self):
        if pmcharge is None:
            raise ImportError(
                "PACMAN-charge is required for clean_pacman. "
                "Install it with 'pip install PACMAN-charge'."
            )
        source = os.path.abspath(os.fspath(self.structure))
        with tempfile.TemporaryDirectory(prefix="coremof_pacman_curation_") as workdir:
            isolated = os.path.join(workdir, "input.cif")
            shutil.copy2(source, isolated)
            pmcharge.predict(
                cif_file=isolated,
                charge_type="DDEC6",
                digits=10,
                atom_type=True,
                neutral=False,
                keep_connect=False,
            )
            predicted = os.path.join(workdir, "input_pacman.cif")
            if not os.path.isfile(predicted):
                raise RuntimeError("PACMAN did not create the expected charged CIF")
            _validated_pacman_charges(source, predicted, read(source))
            destination = os.path.join(
                self.output, os.path.splitext(os.path.basename(source))[0] + "_pacman.cif"
            )
            shutil.copy2(predicted, destination)

    def process(self):
        fsr_skin, fsr_solvent, fsr_ion, fsr_ion_charge = self.run_fsr()
        self.asr_status = "NOT_AVAILABLE"
        self.asr_diagnostic = "ASR_UNSUPPORTED: coordinated-solvent removal is not implemented"
        warnings.warn(self.asr_diagnostic, RuntimeWarning, stacklevel=2)
        if self.saveto:
            mode = 'a' if os.path.exists(self.csv_path) else 'w'
            with open(self.csv_path, mode=mode, newline='') as f:
                writer = csv.writer(f)
                if mode == 'w':
                    writer.writerow(["Name", "Skin_FSR", "Skin_ASR", "FSR_Solvent", "ASR_Solvent", "FSR_Ion", "ASR_Ion", "FSR_Ion_Charge", "ASR_Ion_Charge"])
                writer.writerow([
                    os.path.basename(self.structure), str(fsr_skin), "ASR_UNSUPPORTED",
                    fsr_solvent, "ASR_UNSUPPORTED", fsr_ion, "ASR_UNSUPPORTED",
                    fsr_ion_charge, "ASR_UNSUPPORTED",
                ])

    def run_fsr(self):
        return self.run_clean(mode="FSR")

    def run_asr(self):
        return self.run_clean(mode="ASR")

    def run_clean(self, mode="FSR"):
        if mode == "ASR":
            raise NotImplementedError("PACMAN coordinated-solvent removal (ASR) is not implemented")
        if mode != "FSR":
            raise ValueError("mode must be FSR or ASR")
        file_prefix = os.path.splitext(os.fspath(self.structure))[0]
        skin = self.initial_skin
        clean_func = self.free_clean if mode == "FSR" else self.all_clean

        while True:
            result = clean_func(file_prefix, self.output, skin)
            if result is None:
                raise RuntimeError("PACMAN FSR curation did not produce a result")
            cleaned_skin, solvents, ions, ion_charges = result
            has_metals = any(
                any(e in self.metal_list for e in re.findall(r'([A-Z][a-z]?)\d*', formula))
                for formula in solvents
            )
            if has_metals:
                skin += 0.05
            else:
                break
        return skin, solvents, ions, ion_charges

    def build_ASE_neighborlist(self, cif, skin):
        radii = [self.cambridge_radii[i] for i in cif.get_chemical_symbols()]
        neighborlist = NeighborList(radii, self_interaction=False, bothways=True, skin=skin)
        neighborlist.update(cif)
        return neighborlist

    def find_clusters(self, adjacency_matrix, atom_count):
        _, labels = connected_components(adjacency_matrix, directed=True)
        return [[i for i in range(atom_count) if labels[i] == n] for n in set(labels)]

    def CustomMatrix(self, neighborlist, atom_count):
        mat = np.zeros((atom_count, atom_count), dtype=int)
        for i in range(atom_count):
            neighbors, _ = neighborlist.get_neighbors(i)
            for j in neighbors:
                mat[i][j] = 1
        return mat

    def cluster_to_formula(self, cluster, atoms):
        symbols = [atoms[i].symbol for i in cluster]
        count = collections.Counter(symbols)
        return ''.join([el + (str(count[el]) if count[el] > 1 else '') for el in sorted(count)])

    def free_clean(self, input_file, save_folder, skin):
        cif = read(input_file + ".cif")
        charges = _validated_pacman_charges(
            input_file + ".cif",
            os.path.join(save_folder, os.path.basename(input_file) + "_pacman.cif"),
            cif,
        )

        try:
            neighborlist = self.build_ASE_neighborlist(cif, skin)
            matrix = self.CustomMatrix(neighborlist, len(cif))
            clusters = sorted(self.find_clusters(matrix, len(cif)), key=lambda x: len(x), reverse=True)

            main_clusters, ions, solvents = [], [], []
            ion_formulas, ion_charges = [], []

            for cl in clusters:
                formula = self.cluster_to_formula(cl, cif)
                cluster_charge = math.fsum(charges[i] for i in cl)
                if not math.isfinite(cluster_charge):
                    raise ValueError("PACMAN component charge is not finite")

                if not main_clusters:
                    main_clusters.append(cl)
                elif len(cl) > 0.5 * len(main_clusters[0]):
                    main_clusters.append(cl)
                elif abs(cluster_charge) > 0.1:
                    ions.append(cl)
                    ion_formulas.append(formula)
                    ion_charges.append(cluster_charge)
                else:
                    solvents.append(formula)

            final_atoms = list(itertools.chain.from_iterable(main_clusters + ions))
            suffix = "_FSR_ION.cif" if ions else "_FSR.cif"
            write(os.path.join(save_folder, os.path.basename(input_file) + suffix), cif[final_atoms])
            return skin, solvents, ion_formulas, ion_charges
        except Exception as e:
            raise RuntimeError(f"PACMAN FSR curation failed for {input_file}: {e}") from e

    def all_clean(self, input_file, save_folder, skin):
        raise NotImplementedError("PACMAN coordinated-solvent removal (ASR) is not implemented")

def run_MOSAEC(*args, **kwargs):
    """Retired execution interface, retained only for a clear migration error."""
    from ._checker_execution import results_only
    return results_only("MOSAEC")



def run_mofclassifier(cif_folder, save_path="./mofclassifier_results.json", model="core", batch_size=64, *, overwrite=False):
    """Check MOF by MOFClassifier: https://github.com/Chung-Research-Group/MOFClassifier. Ref: https://doi.org/10.1021/jacs.5c10126        
    
    Args:
        cif_folder (str): path to the folder including all CIFs.
        save_path (str): path to save the predictions.
        model (str): the name of model used for predictions.
        batch_size (int): batch size for predicting.
        overwrite (bool): explicitly replace an existing result file.

    Source CIFs are never given to the upstream mutable parser. Models must
    already be installed; importing this module does not download them.
    The upstream batch mean is retained. To reproduce the recorded release's
    exact 100-bag CPU mean, use ``release_mofclassifier`` instead.

    Returns:
        dict:
            -   results of MOFClassifier.
    """
    from ._mofclassifier import predict_directory
    return predict_directory(cif_folder, save_path, model, batch_size, overwrite=overwrite)
