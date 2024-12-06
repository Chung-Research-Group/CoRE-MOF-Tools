from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.structure_matcher import StructureMatcher
import pandas as pd

def check_same(cif_1, cif_2):
    atoms1 = read(cif_1)
    stru_1 = AseAtomsAdaptor.get_structure(atoms1)
    atoms2 = read(cif_2)
    stru_2 = AseAtomsAdaptor.get_structure(atoms2)
    checker = StructureMatcher(
        ltol=0.2, stol=0.3, angle_tol=5, primitive_cell=True, scale=True, 
        attempt_supercell=True, allow_subset=False, comparator=None, 
        supercell_size='volume', ignored_species=()
    )
    result = checker.fit(stru_1, stru_2)
    return result

data_name = pd.read_csv("ASR_data.csv")["refcode"].iloc[:]

for refcode in data_name:
    name = refcode.replace("_ASR_pacman", "")
    cif_1 = f"dataset-ASR/{refcode}.cif"
    cif_2 = f"ori/{name}.cif"
    result = check_same(cif_1, cif_2)
    print(name,result)
