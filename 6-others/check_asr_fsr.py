import fireducks.pandas as pd

from tqdm import tqdm

import shutil

from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.structure_matcher import StructureMatcher

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


ASR_data = pd.read_csv("../all/ASR_data_20250124_internal.csv")
FSR_data = pd.read_csv("../all/FSR_data_20250124_internal.csv")

same_list = []

for i, id in tqdm(enumerate(ASR_data["coreid"][:])):
    refcode_asr = ASR_data["refcode"][i]
    refcode_asr2fsr = refcode_asr.replace("_ASR","_FSR")
    indx_FSR = FSR_data.index[FSR_data["refcode"] == refcode_asr2fsr].tolist()
    if len(FSR_data["coreid"][indx_FSR])>0:
        FSR_id = FSR_data["coreid"][indx_FSR].iloc[0]
        result = check_same("CoREMOF2024DB/ASR/"+id+".cif","CoREMOF2024DB/FSR/"+FSR_id+".cif")
        if result:
            print(id, result)
            shutil.move("CoREMOF2024DB/FSR/"+FSR_id+".cif", "./FSR_remove/")
            same_list.append([id, FSR_id])
        else:
            print(id, result)
    else:
        print(id, "unique")
df_same_list = pd.DataFrame(same_list, columns = ["ASR_id", "FSR_id"])
df_same_list.to_csv("check_results.csv",index=False)