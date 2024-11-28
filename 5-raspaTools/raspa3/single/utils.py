import numpy as np
from PACMANCharge import pmcharge
import re

def pre_charge(cif_file,charge_type="DDEC6",digits=10,atom_type=True,neutral=True):
    try:
        pmcharge.predict(cif_file,charge_type=charge_type,digits=digits,atom_type=atom_type,neutral=neutral)
    except:
        print(cif_file,"fail")

def unit_cell(cif_file,cutoff =14):
    with open(cif_file, 'r') as cif:
        mof_cif= cif.read()
    for line in mof_cif.split("\n"):
        if "_cell_length_a" in line:
            length_a = line.split()[1]
            length_a =float(length_a)
        if "_cell_length_b" in line:
            length_b = line.split()[1]
            length_b = float(length_b)
        if "_cell_length_c" in line:
            length_c= line.split()[1]
            length_c= float(length_c)
        if "_cell_angle_alpha" in line:
            alpha = line.split()[1]
            alpha = float(alpha)
        if "_cell_angle_beta" in line:
            beta= line.split()[1]
            beta= float(beta)
        if "_cell_angle_gamma" in line:
            gamma = line.split()[1]
            gamma = float(gamma)
    
    ax = length_a
    ay = 0.0
    az = 0.0
    bx = length_b * np.cos(gamma * np.pi / 180.0)
    by = length_b * np.sin(gamma * np.pi / 180.0)
    bz = 0.0
    cx = length_c * np.cos(beta * np.pi / 180.0)
    cy = (length_c * length_b * np.cos(alpha * np.pi /180.0) - bx * cx) / by
    cz = (length_c ** 2 - cx ** 2 - cy ** 2) ** 0.5
    
    unit_cell =  np.asarray([[ax, ay, az],[bx, by, bz], [cx, cy, cz]])
    A = unit_cell[0]
    B = unit_cell[1]
    C = unit_cell[2]

    Wa = np.divide(np.linalg.norm(np.dot(np.cross(B,C),A)), np.linalg.norm(np.cross(B,C)))
    Wb = np.divide(np.linalg.norm(np.dot(np.cross(C,A),B)), np.linalg.norm(np.cross(C,A)))
    Wc = np.divide(np.linalg.norm(np.dot(np.cross(A,B),C)), np.linalg.norm(np.cross(A,B)))
    
    uc_x = int(np.ceil(cutoff/(0.5*Wa)))
    uc_y = int(np.ceil(cutoff/(0.5*Wb)))
    uc_z = int(np.ceil(cutoff/(0.5*Wc)))
    
    return [uc_x,uc_y,uc_z]

def fix_atom_site_labels(struc_name):
    second_loop_found = False
    start_processing = False
    new_lines = []

    with open(struc_name+".cif", 'r') as file:
        for line in file:
            if 'loop_' in line:
                if second_loop_found:
                    start_processing = True
                else:
                    second_loop_found = True
            elif start_processing:
                if 'loop_' in line:
                    start_processing = False
                else:
                    parts = line.split()
                    if len(parts) > 1:
                        parts[1] = re.sub(r'\d+', '', parts[1]) 
                    new_lines.append("  ".join(parts) + "\n")
            else:
                new_lines.append(line)
  
    return new_lines
