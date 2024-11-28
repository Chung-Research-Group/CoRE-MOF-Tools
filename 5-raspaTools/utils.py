import numpy as np
from PACMANCharge import pmcharge


def unit_cell(cif_file,cutoff):
    with open(cif_file+".cif", 'r') as cif:
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
    
    return str(uc_x) + "\t" + str(uc_y) + "\t" + str(uc_z)


def predict_charge(cif_file):
    pmcharge.predict(cif_file+".cif",charge_type="DDEC6",digits=10,atom_type=True,neutral=True,keep_connect=False)


def load_lists(*keys, kwargs):
    lists = []
    for key in keys:
        if key in kwargs:
            with open(kwargs[key], "r") as f:
                lists.append([line.strip() for line in f])
        else:
            lists.append(None)
    return tuple(lists)