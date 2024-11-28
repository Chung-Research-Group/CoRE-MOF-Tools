import os, re, csv, glob, functools, warnings, itertools, collections
import numpy as np
from .utils.atoms_definitions import METAL, COVALENTRADII
from .utils.ions_list import ALLIONS
from ase.io import read, write
from ase.neighborlist import NeighborList
from scipy.sparse.csgraph import connected_components
warnings.filterwarnings('ignore')

cambridge_radii = COVALENTRADII
metal_list = [element for element, is_metal in METAL.items() if is_metal]
ions_list = set(ALLIONS)

class run():
    def __init__(self, input_folder, output_folder):
        self.input_folder = input_folder
        self.input_files = glob.glob(os.path.join(input_folder, "*.cif"))
        self.output = output_folder
        self.csv_path = os.path.join(output_folder, "clean_log.csv")
        self.process()

    def process(self):
        with open(self.csv_path, mode='w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(["Structure", "FSR Printed Formulas", "ASR Printed Formulas"])

            for mof in self.input_files:
                fsr_results = self.run_fsr(mof, self.output)
                print(f"FSR results for {mof}: {fsr_results}")
                
                asr_results = self.run_asr(mof, self.output)
                print(f"ASR results for {mof}: {asr_results}")

                csv_writer.writerow([os.path.basename(mof), "; ".join(fsr_results), "; ".join(asr_results)])

    def run_fsr(self, mof, save_folder, initial_skin=0.1, metal_list=metal_list):
        m = mof.replace(".cif", "")
        skin = initial_skin
        try:
            while True:
                printed_formulas = self.free_clean(m, save_folder, ions_list, skin)
                has_metals = False
                for e_s in printed_formulas:
                    split_formula = re.findall(r'([A-Z][a-z]?)(\d*)', e_s)
                    elements = [match[0] for match in split_formula]
                    if any(e in metal_list for e in elements):
                        has_metals = True
                        skin += 0.05
                if not has_metals:
                    break
        except Exception as e:
            print(m, str(e))

        return printed_formulas

    def run_asr(self, mof, save_folder, initial_skin=0.1, metal_list=metal_list):
        m = mof.replace(".cif", "")
        skin = initial_skin
        try:
            while True:
                printed_formulas = self.all_clean(m, save_folder, ions_list, skin)
                has_metals = False
                for e_s in printed_formulas:
                    split_formula = re.findall(r'([A-Z][a-z]?)(\d*)', e_s)
                    elements = [match[0] for match in split_formula]
                    if any(e in metal_list for e in elements):
                        has_metals = True
                        skin += 0.05
                if not has_metals:
                    break
            return printed_formulas
        except Exception as e:
            print(f"{m} Fail: {e}")

    def build_ASE_neighborlist(self, cif, skin):
        radii = [cambridge_radii[i] for i in cif.get_chemical_symbols()]
        ASE_neighborlist = NeighborList(radii, self_interaction=False, bothways=True, skin=skin)
        ASE_neighborlist.update(cif)
        return ASE_neighborlist

    def find_clusters(self, adjacency_matrix, atom_count):
        clusters = []
        cluster_count, clusterIDs = connected_components(adjacency_matrix, directed=True)
        for n in range(cluster_count):
            clusters.append([i for i in range(atom_count) if clusterIDs[i] == n])
        return clusters

    def find_metal_connected_atoms(self, struct, neighborlist):
        metal_connected_atoms = []
        metal_atoms = []
        for i, elem in enumerate(struct.get_chemical_symbols()):
            if elem in metal_list:
                neighbors, _ = neighborlist.get_neighbors(i)
                metal_connected_atoms.append(neighbors)
                metal_atoms.append(i)
        return metal_connected_atoms, metal_atoms, struct

    def CustomMatrix(self, neighborlist, atom_count):
        matrix = np.zeros((atom_count, atom_count), dtype=int)
        for i in range(atom_count):
            neighbors, _ = neighborlist.get_neighbors(i)
            for j in neighbors:
                matrix[i][j] = 1
        return matrix

    def mod_adjacency_matrix(self, adj_matrix, MetalConAtoms, MetalAtoms, atom_count, struct):
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
        symbols = [cif[i].symbol for i in cluster]
        count = collections.Counter(symbols)
        formula = ''.join([atom + (str(count[atom]) if count[atom] > 1 else '') for atom in sorted(count)])
        return formula

    def free_clean(self, input_file, save_folder, ions_list, skin):
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
                new_fn = refcode + '_ION_F.cif'
            else:
                refcode=refcode.replace("_FSR","")
                new_fn = refcode + '_FSR.cif'
            write(save_folder+new_fn, cif[solvated_merged])
        
            return printed_formulas
        except Exception as e:
            print(f"{input_file} Fail: {e}")

    def all_clean(self, input_file, save_folder, ions_list, skin):
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
                new_fn = refcode + '_ION_A.cif'
            else:
                new_fn = refcode + "_ASR.cif"
            write(save_folder+new_fn, struct[final_merged])
            return printed_formulas
    
        except Exception as e:
            print(f"{input_file} Fail: {e}")
