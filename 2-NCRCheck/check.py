import os, csv, glob, warnings
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from .mofchecker import MOFChecker
import warnings
warnings.filterwarnings('ignore')

ATR = {'H':0.38, 'Li':0.86, 'Be':0.53, 'B':1.01, 'C':0.88, 'N':0.86, 'O':0.89, 'F':0.82, 'Na':1.15, 'Mg':1.28, 'Al':1.53, 'Si':1.38, 'P':1.28,
      'S':1.20, 'Cl':1.17, 'K':1.44, 'Ca':1.17, 'Sc':1.62, 'Ti':1.65, 'V':1.51, 'Cr':1.53, 'Mn':1.53, 'Fe':1.43, 'Co':1.31, 'Ni':1.33, 'Cu':1.31,
      'Zn':1.41, 'Ga':1.40, 'Ge':1.35, 'As':1.39, 'Se':1.40, 'Br':1.39, 'Rb':1.65, 'Sr':1.30, 'Y':1.84, 'Zr':1.73, 'Nb':1.66, 'Mo':1.57,
      'Ru':1.58, 'Rh':1.63, 'Pd':1.68, 'Ag':1.56, 'Cd':1.56, 'In':1.53, 'Sn':1.64, 'Sb':1.64, 'Te':1.65, 'I':1.58, 'Cs':1.85, 'Ba':1.52,
      'La':1.91, 'Ce':1.98, 'Pr':1.75, 'Nd':1.92, 'Sm':1.89, 'Eu':1.83, 'Gd':1.79, 'Tb':1.82, 'Dy':1.79, 'Ho':1.63, 'Er':1.80, 'Tm':1.84,
      'Yb':1.80, 'Lu':1.86, 'Hf':1.73, 'W':1.33, 'Re':1.29, 'Ir':1.50, 'Pt':1.66, 'Au':1.68, 'Hg':1.88, 'Pb':1.72, 'Bi':1.72, 'Th':1.97,
      'U':1.76, 'Np':1.73, 'Pu':1.71}
Coef_A = {'H':-0.6093, 'B':-2.2011, 'C':-1.2685, 'N':-1.2680, 'O':-1.0525, 'Cl':-0.7621, 'Br':-0.8003}
Coef_C = {'H':0.5927, 'B':3.4380, 'C':1.8855, 'N':1.8401, 'O':1.5189, 'Cl':1.3723, 'Br':1.5272}
BO_list = ['H','B','C','N','O','Cl','Br']
metals = ['Li','Be','Na','Mg','Al','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Rb','Sr','Y','Zr','Nb','Mo','Ru','Rh',
         'Pd','Ag','Cd','In','Sn','Cs','Ba','La','Ce','Pr','Nd','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','W','Re','Re','Ir',
         'Pt','Au','Hg','Pb','Bi','Th','U','Np','Pu'] # do not include metalloids

class run():
    def __init__(self, input_folder, output_folder):
        self.input_folder = input_folder
        self.input_files = glob.glob(os.path.join(input_folder, "*.cif"))
        self.output = output_folder
        self.process()

    def process(self):
        chen_manz_csv_path = os.path.join(self.output, "chen_manz_checker_log.csv")
        mof_checker_csv_path = os.path.join(self.output, "mof_checker_log.csv")

        with open(chen_manz_csv_path, mode="w", newline="") as chen_manz_file:
            chen_manz_writer = csv.writer(chen_manz_file)
            chen_manz_writer.writerow(["Name", "Isolated Atom", "Overlapping Atoms", "Undercoordinated Carbon", "Overcoordinated Carbon"])

            with open(mof_checker_csv_path, mode="w", newline="") as mof_checker_file:
                mof_checker_writer = csv.writer(mof_checker_file)
                mof_checker_writer.writerow(["Name", "MOF Checker Issues"])

                for cif_file in self.input_files:
                    chen_manz_result = self.Chen_Manz(cif_file)
                    chen_manz_writer.writerow(chen_manz_result)
                    mof_checker_result = self.mof_checker(cif_file)
                    mof_name = chen_manz_result[0]
                    mof_issues = mof_checker_result.get(mof_name, "no issues")
                    mof_checker_writer.writerow([mof_name, mof_issues])

    def Chen_Manz(self,cif_file):
        try:
            name = cif_file.split(".cif")[0].split("/")[-1]
            atoms = read(cif_file)
            sym = atoms.get_chemical_symbols()
            isolated = False
            overlapping = False
            # misplaced_hydro = False
            under_carbon = False
            over_carbon = False
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
                        overlapping = True
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
                            under_carbon = True
                        elif sum_BO >= 5.5:
                            over_carbon = True
                if len(nl) == 0:
                    isolated = True
                # if sym[a] == 'H':
                    # if ('N' in H_connected) or ('O' in H_connected):
                        # common_metals = [e for e in H_connected if e in metals]
                        # if len(common_metals) != 0:
                        #     misplaced_hydro = True
            return [name,isolated, overlapping,  under_carbon, over_carbon] # misplaced_hydro
        except:
            return [name,"unknown", "unknown",  "unknown", "unknown"] # "unknown", 
        

    def mof_checker(self,mof):
        name = mof.split("/")[-1].replace(".cif", "")
        results = {}
        try:
            atoms = read(mof)
            structure = AseAtomsAdaptor.get_structure(atoms)
            checker = MOFChecker(structure)
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
                results[name] = has_problem
            else:
                results[name] = "no issues"
        except Exception as e:
            results[name] = "fail"
            print(f"{name} failed due to {str(e)}")

        return results
