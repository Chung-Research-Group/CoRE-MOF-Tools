import os
from fragment_MOFs_for_pormake import make_MOF_fragments
import pandas as pd
import warnings
from ase.io import read, write
import shutil
import numpy as np
import glob
from pymatgen.core.structure import Structure
from pymatgen.analysis import structure_matcher
import pymatgen.core as mg
from tqdm import tqdm

names = pd.read_excel("./mofs.xlsx")

def ase_format(root_dir):
    cifs = os.listdir(root_dir)
    cifs.sort()
    for cif in cifs:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                mof = read(os.path.join(root_dir, cif))
                write(os.path.join(root_dir, cif), mof)
                print('ASE reading: ' + cif)
        except:
            print('ASE reading fail: ' + cif)

def conv2qmof(root_dir):
    cifs = os.listdir(root_dir)
    cifs.sort()
    mofs = []
    for cif in cifs:
        mof_temp = Structure.from_file(os.path.join(root_dir, cif),primitive=True)
        mofs.append(mof_temp)
        print('trans to primituve cell: ' + cif)
    sm = structure_matcher.StructureMatcher(primitive_cell=True)
    groups = sm.group_structures(mofs)
    for group in groups:
        mof_temp = group[0]
        mof_temp.to(filename=os.path.join(root_dir,cif))



fragmentation_directory = './coremof/FSR_clean/'

ase_format(fragmentation_directory+'cif/')
conv2qmof(fragmentation_directory+'cif/')



log_list = []
for cif_file in tqdm(os.listdir(fragmentation_directory)): # assumes that all of the cifs are in a directory called cif inside of the parent directory.
    # print('now on', cif_file)  
    try:
        return_code = make_MOF_fragments(fragmentation_directory+cif_file,
                                         path=fragmentation_directory+'/',
                                         xyzpath=fragmentation_directory+'/xyz/'+cif_file.replace('cif','xyz'))
        if return_code == None:
            return_code = 'good'
        log_list.append({'cif':cif_file,'return_code':return_code})
    except:
        print(cif_file,"fail")

df = pd.DataFrame(log_list)
df.to_csv(fragmentation_directory+'/fragment_status.csv')



def should_delete_line(previous_line, current_line, line_index):
    if line_index <= 2:
        return False
    if len(previous_line.split()) == 4 and len(current_line.split()) <= 3:
        return True
    return False


periodic_table_symbols = [
    'X','H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg',
    'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr',
    'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
    'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po',
    'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm',
    'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs',
    'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'
    ]


names = ["ABIWAL_FSR"]


for name in qmofid["name"]:
# for name in names:
    lines_to_keep = []
    previous_line = None

    try:
        # linker data
        lines_to_keep = []
        line_index = 1
        with open("./qmofid/linkers/" + name + '_linker_0.xyz', 'r') as file:
            for line in file:
                if previous_line and should_delete_line(previous_line, line, line_index):
                    continue 
                else:
                    lines_to_keep.append(line) 
                line_index += 1
                previous_line = line 
        with open("./qmofid/linkers/" + name + '_linker_0_4ase.xyz', 'w') as file:
            file.writelines(lines_to_keep)
        file.close()
        atoms1 = read("./qmofid/linkers/" + name + '_linker_0_4ase.xyz')
        pos1 = atoms1.positions
        atom1 = atoms1.get_chemical_symbols()
        linkerinfo = []
        for (i,ele1) in enumerate(atom1):
            atom_index1 = periodic_table_symbols.index(ele1)
            x1 = pos1[i][0]
            y1 = pos1[i][1]
            z1 = pos1[i][2]
            atom_per1 = [atom_index1,x1,y1,z1]
            linkerinfo.append(atom_per1)
        np.save("./qmofid/input/" + name + "_linker.npy", linkerinfo)
    
        # sbu data
        lines_to_keep = []
        line_index = 1
        with open("./qmofid/sbus/" + name + '_sbu_0.xyz', 'r') as file:
            for line in file:
                if previous_line and should_delete_line(previous_line, line, line_index):
                    continue 
                else:
                    lines_to_keep.append(line) 
                line_index += 1
                previous_line = line 
        with open("./qmofid/sbus/" + name + '_sbu_0_4ase.xyz', 'w') as file:
            file.writelines(lines_to_keep)
        file.close()
        atoms2 = read("./qmofid/sbus/" + name + '_sbu_0_4ase.xyz')
        pos2 = atoms2.positions
        atom2 = atoms2.get_chemical_symbols()
        sbuinfo = []
        for (i,ele2) in enumerate(atom2):
            atom_index2 = periodic_table_symbols.index(ele2)
            x2 = pos2[i][0]
            y2 = pos2[i][1]
            z2 = pos2[i][2]
            atom_per2 = [atom_index2,x2,y2,z2]
            sbuinfo.append(atom_per2)
        np.save("./qmofid/input/" + name + "_sbu.npy", sbuinfo)
    
        #topology data
        topo_name = qmofid[qmofid["name"]==name]["info.mofid.topology"]
        topo_struc = "./qmofid/topo/" + topo_name.values[0] + ".cif"
        atoms3 = read(topo_struc)
        pos3 = atoms3.positions
        atom3 = atoms3.get_chemical_symbols()
        topoinfo = []
        for (i,ele3) in enumerate(atom3):
            atom_index3 = periodic_table_symbols.index(ele3)
            x3 = pos3[i][0]
            y3 = pos3[i][1]
            z3 = pos3[i][2]
            atom_per3 = [atom_index3,x3,y3,z3]
            topoinfo.append(atom_per3)
        np.save("./qmofid/input/" + name + "_topo.npy", topoinfo)
    
        # atom data
        atoms4 = read("./qmofid/cif/" + name +".cif")
        pos4 = atoms4.positions
        atom4 = atoms4.get_chemical_symbols()
        atominfo = []
        for (i,ele4) in enumerate(atom4):
            atom_index4 = periodic_table_symbols.index(ele4)
            x4 = pos4[i][0]
            y4 = pos4[i][1]
            z4 = pos4[i][2]
            atom_per4 = [atom_index4,x4,y4,z4]
            atominfo.append(atom_per4)
        np.save("./qmofid/input/" + name + "_atom.npy", atominfo)
    
        # tar: pos & cell
        atoms5 = read("../qmof/relaxed/" + name +".cif")
        pos5 = atoms5.positions
        atom5 = atoms5.get_chemical_symbols()
        postarinfo = []
        for (i,ele5) in enumerate(atom5):
            x5 = pos5[i][0]
            y5 = pos5[i][1]
            z5 = pos5[i][2]
            atom_per5 = [x5,y5,z5]
            postarinfo.append(atom_per5)
        np.save("./qmofid/output/" + name + "_pos.npy", postarinfo)
        
        structure = mg.Structure.from_file("../qmof/relaxed/" + name +".cif")
        lattice = structure.lattice.matrix
        np.save("./qmofid/output/" + name + "_cell.npy", lattice)
    except:
        print(name)


np.load("./qmofid/input/" + name + "_linker.npy")