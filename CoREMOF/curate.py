"""Process your CIF to "CoRE MOF" CIF.
"""

import os, re, csv, json, requests, functools, warnings, itertools, collections, subprocess
from ase.io import read, write

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from CoREMOF.utils.atoms_definitions import METAL, COVALENTRADII
from CoREMOF.utils.ions_list import ALLIONS
import numpy as np
import pandas as pd

from ase.neighborlist import NeighborList
from scipy.sparse.csgraph import connected_components
warnings.filterwarnings('ignore')


from pymatgen.io.ase import AseAtomsAdaptor
from mofchecker import MOFChecker
from CoREMOF.utils.atoms_definitions import ATR, Coef_A, Coef_C #, BO_list, metals4check

from gemmi import cif as CIF
from PACMANCharge import pmcharge


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
        initial_skin (float): skin distance is added to the sum of vdW radii of two atoms.
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

class mof_check():

    """Not Computation-Ready (NCR) MOFs classification.

    Args:
        structure (str): path to your CIF.
        output_folder (str): the path to save checking result.

    Returns:
        Dictionary:
            -   result of NCR classification.
    """
    
    def __init__(self, structure, output_folder="result_curation"):
        self.output = output_folder + os.sep
        self.structure = structure
        os.makedirs(self.output, exist_ok=True)
        self.check()

    def check(self):

        """run checking.
        """
        
        result_check = {}

        chen_manz_result = self.Chen_Manz(self.structure)
        mof_checker_result = self.mof_checker(self.structure)

        result_check["Chen_Manz"] = chen_manz_result
        result_check["mofchecker"] = mof_checker_result

        with open(self.output + os.path.basename(self.structure).replace(".cif","") + "_Chen_Manz_mofchecker.json", "w") as f:
            json.dump(result_check,f,indent=2)

    def Chen_Manz(self, structure):

        """checking MOF by Chen and Manz method: RSC Adv., 2019,9, 36492-36507. https://doi.org/10.1039/C9RA07327B.
        """

        try:
            
            has_problem =[]

            atoms = read(structure)
            sym = atoms.get_chemical_symbols()
            
            for a in range(len(atoms)):
                H_connected = []
                nl = []
                for b in range(len(atoms)):
                    if a == b:
                        continue
                    d = atoms.get_distance(a, b, mic = True)
                    if sym[a] == 'H':
                        if d <= (0.3 + ATR[sym[a]] + ATR[sym[b]]):
                            H_connected.append(b)
                    if d < 0.5*(ATR[sym[a]] + ATR[sym[b]]):
                        has_problem.append("overlapping")
                    if d <= (ATR[sym[a]] + ATR[sym[b]]):
                        nl.append(b)
                if sym[a] == 'C':
                    bonded_ele = [sym[e] for e in nl if sym[e] not in Coef_A]
                    if len(bonded_ele) != 0:
                        pass
                    else:
                        BO = []
                        for bidx in range(len(nl)):
                            b = nl[bidx]
                            d = atoms.get_distance(a, b, mic = True)
                            BO_ab = 10**(Coef_A[sym[b]]*d + Coef_C[sym[b]])
                            if sym[b] == 'H':
                                if BO_ab > 1.25:
                                    BO_ab = 1.25
                            BO.append(BO_ab)
                        sum_BO = sum(BO)
                        if sum_BO < 3.3:
                            has_problem.append("under_carbon")
                        elif sum_BO >= 5.5:
                            has_problem.append("over_carbon")
                if len(nl) == 0:
                    has_problem.append("isolated")

            if len(has_problem) > 0:
                return list(set(has_problem))
            else:
                return ["no issues"]
        except Exception as e:
            return ["unknown"]
        

    def mof_checker(self, structure):
        
        """checking MOF by mofchecker 1.0: https://github.com/lamalab-org/mofchecker.
        """

        try:
            atoms = read(structure)
            structure_ = AseAtomsAdaptor.get_structure(atoms)
            checker = MOFChecker(structure_)
            check_result = checker.get_mof_descriptors()

            has_problem = []
            problem_keys_true = [
                "has_atomic_overlaps", "has_overcoordinated_c", "has_overcoordinated_n",
                "has_overcoordinated_h", "has_suspicious_terminal_oxo",
                "has_undercoordinated_c", "has_undercoordinated_n",
                "has_undercoordinated_rare_earth", "has_undercoordinated_alkali_alkaline",
                "has_geometrically_exposed_metal", 
                "has_lone_molecule"
            ]
            problem_keys_false = ["has_metal", "has_carbon"]

            for key in problem_keys_true:
                if check_result.get(key, False):
                    has_problem.append(key)

            for key in problem_keys_false:
                if not check_result.get(key, True):
                    has_problem.append(key)

            if len(has_problem) > 0:
                return list(set(has_problem))
            else:
                return ["no issues"]
        except Exception as e:
            return ["unknown"]



class clean_pacman():

    """Removing free solvent and coordinated solvent but keep ions based on PACMAN-charge.         
    
    Args:
        structure (str): path to your CIF.
        initial_skin (float): skin distance is added to the sum of vdW radii of two atoms.
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
        pmcharge.predict(
            cif_file=self.structure,
            charge_type="DDEC6",
            digits=10,
            atom_type=True,
            neutral=False,
            keep_connect=False
        )
        src = self.structure.replace(".cif", "_pacman.cif")
        dst = os.path.join(self.output, os.path.basename(src))
        if os.path.exists(src):
            os.rename(src, dst)

    def process(self):
        try:
            fsr_skin, fsr_solvent, fsr_ion, fsr_ion_charge = self.run_fsr()
            asr_skin, asr_solvent, asr_ion, asr_ion_charge = self.run_asr()

            if self.saveto:
                mode = 'a' if os.path.exists(self.csv_path) else 'w'
                with open(self.csv_path, mode=mode, newline='') as f:
                    writer = csv.writer(f)
                    if mode == 'w':
                        writer.writerow(["Name", "Skin_FSR", "Skin_ASR", "FSR_Solvent", "ASR_Solvent", "FSR_Ion", "ASR_Ion", "FSR_Ion_Charge", "ASR_Ion_Charge"])
                    writer.writerow([
                        os.path.basename(self.structure), str(fsr_skin), str(asr_skin),
                        fsr_solvent, asr_solvent, fsr_ion, asr_ion, fsr_ion_charge, asr_ion_charge
                    ])
            os.remove(os.path.join(self.output, os.path.basename(self.structure.replace(".cif","")) + "_pacman.cif"))
        except Exception as e:
            print("[clean_pacman.process] Error:", e)

    def run_fsr(self):
        return self.run_clean(mode="FSR")

    def run_asr(self):
        return self.run_clean(mode="ASR")

    def run_clean(self, mode="FSR"):
        file_prefix = self.structure.replace(".cif", "")
        skin = self.initial_skin
        clean_func = self.free_clean if mode == "FSR" else self.all_clean

        while True:
            result = clean_func(file_prefix, self.output, skin)
            if result is None:
                break
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
        try:
            cif = read(input_file + ".cif")
            charges = list(CIF.read_file(os.path.join(save_folder, os.path.basename(input_file) + "_pacman.cif")).sole_block().find_loop('_atom_site_charge'))

            neighborlist = self.build_ASE_neighborlist(cif, skin)
            matrix = self.CustomMatrix(neighborlist, len(cif))
            clusters = sorted(self.find_clusters(matrix, len(cif)), key=lambda x: len(x), reverse=True)

            main_clusters, ions, solvents = [], [], []
            ion_formulas, ion_charges = [], []

            for cl in clusters:
                formula = self.cluster_to_formula(cl, cif)
                cluster_charge = sum([float(charges[i]) for i in cl if i < len(charges)])

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
            print("[free_clean]", input_file, "failed:", e)

    def all_clean(self, input_file, save_folder, skin):
        try:
            cif = read(input_file + ".cif")
            charges = list(CIF.read_file(os.path.join(save_folder, os.path.basename(input_file) + "_pacman.cif")).sole_block().find_loop('_atom_site_charge'))

            neighborlist = self.build_ASE_neighborlist(cif, skin)
            matrix = self.CustomMatrix(neighborlist, len(cif))
            clusters = sorted(self.find_clusters(matrix, len(cif)), key=lambda x: len(x), reverse=True)

            main_clusters, ions, solvents = [], [], []
            ion_formulas, ion_charges = [], []

            for cl in clusters:
                formula = self.cluster_to_formula(cl, cif)
                cluster_charge = sum([float(charges[i]) for i in cl if i < len(charges)])

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
            suffix = "_ASR_ION.cif" if ions else "_ASR.cif"
            write(os.path.join(save_folder, os.path.basename(input_file) + suffix), cif[final_atoms])
            return skin, solvents, ion_formulas, ion_charges
        except Exception as e:
            print("[all_clean]", input_file, "failed:", e)

try:
    from ccdc import io
except:
    print("Before using MOSAEC to check your structure, please install CSD Python API with license")


class MOSAEC():
    """Check MOF by Metal Oxidation State Automated Error Checker [MOSAEC](https://github.com/uowoolab/MOSAEC).         
    
    Args:
        folder (str): path to your folder contain MOFs.
        saveto (bool or str): the name of csv file with clean result.

    Returns:
        CSV:
            -   result of MOSAEC.
    """

    FILES_TO_DOWNLOAD = {
        "Ionization_Energies.csv": "https://raw.githubusercontent.com/uowoolab/MOSAEC/main/Ionization_Energies.csv",
        "KnownON.csv": "https://raw.githubusercontent.com/uowoolab/MOSAEC/main/KnownON.csv",
        "Oxidation_Probabilities.csv": "https://raw.githubusercontent.com/uowoolab/MOSAEC/main/Oxidation_Probabilities.csv",
        "mosaec.py": "https://raw.githubusercontent.com/uowoolab/MOSAEC/main/mosaec.py"
    }

    def __init__(self, folder: str, saveto: str = "mosaec.csv"):
        self.folder = folder
        self.saveto = saveto
        os.makedirs(self.folder, exist_ok=True)

        for filename, url in self.FILES_TO_DOWNLOAD.items():
            path = os.path.join(self.folder, filename)
            if not os.path.exists(path):
                response = requests.get(url)
                if response.status_code == 200:
                    with open(path, "wb") as f:
                        f.write(response.content)
                    print(f"Downloaded: {filename}")
                else:
                    print(f"Failed to download {filename}: {response.status_code}")
            else:
                print(f"Skipped (already exists): {filename}")
        self.run_mosaec()
        if saveto:
            os.rename(self.folder + "/OxStatesOutput.csv", self.folder + "/" + saveto)
    def run_mosaec(self):
        script_path = "mosaec.py"
        print(f"Running {script_path}...")

        result = subprocess.run(
            ["python", script_path],
            cwd=self.folder,
            capture_output=True,
            text=True
        )

        if result.returncode == 0:
            print(result.stdout)
        else:
            print(result.stderr)

    def check_result(self, saveto):
        df = pd.read_csv(saveto)
        status_cols = df.columns[-10:]

        def merge_group(group):
            result = {}
            for col in status_cols:
                if all(group[col] == 'GOOD'):
                    result[col] = 'GOOD'
                else:
                    result[col] = 'NOT_GOOD'
            return pd.Series(result)


        merged_df = df.groupby('CIF').apply(merge_group).reset_index()

        binary_df = merged_df.copy()
        for col in status_cols:
            binary_df[col] = binary_df[col].apply(lambda x: True if x == 'GOOD' else False)

        binary_df.to_csv(saveto, index=False)
