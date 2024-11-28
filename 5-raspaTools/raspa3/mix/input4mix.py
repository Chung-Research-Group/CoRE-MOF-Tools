import os
import json
import shutil
import argparse
from utils import unit_cell,fix_atom_site_labels,get_mol_frac


def main():

    parser = argparse.ArgumentParser(description="Generate input files for RASPA3 software")

    # Add common arguments
    parser.add_argument('--cifs_folder', default="./", type=str, help='Folder of cifs')
    parser.add_argument('--structure_list', default="list", type=str, help='Structure list used for your simulation')
    parser.add_argument('--temperature_list', default="t_list", type=str, help='Temperature list used for your simulation')
    parser.add_argument('--pressure_list', default="p_list", type=str, help='Pressure list used for your simulation')
    parser.add_argument('--gas_list', default="g_list", type=str, help='Gas list used for your simulation')

    parser.add_argument('--have_water', default=False, type=bool, help='Mix has water')
    
    
    parser.add_argument('--fix_atomnumber', default=False, type=bool, help='Remove the number from atom label')
    parser.add_argument('--pacman', default=False, type=bool, help='Get atomic charges for your framework')

    args, _ = parser.parse_known_args()
    cifs_floder = args.cifs_folder
    structure_list = args.structure_list
    temperature_list = args.temperature_list
    pressure_list = args.pressure_list
    gas_list = args.gas_list
    pacman = args.pacman
    fix_atomnumber = args.fix_atomnumber
    have_water = args.have_water
    with open(structure_list) as structures:
        S =  structures.readlines()

    with open(gas_list) as gases:
        G =  gases.readlines()

    with open(pressure_list) as pressures:
        P =  pressures.readlines()

    with open(temperature_list) as temperatures:
        T =  temperatures.readlines()

    if not have_water:
        parser.add_argument('--mol_frac_list', default="mol_frac_list", type=str, help='Mol fraction used for your simulation')
        args = parser.parse_args()
        mol_frac_list = args.mol_frac_list
    else:
        parser.add_argument('--RH', default=0.2, type=float, help='Relative Humidity')
        parser.add_argument('--others', nargs='+', type=float, default=[0.1, 0.9], help='Mol frac for components except water')
        args = parser.parse_args()
        others = args.others
        RH = args.RH
    

    print(args)

    for g in G:
        g = g.replace("\n","")

    for s in S:
        s = s.replace("\n","")
        if not os.path.exists(s):
            os.makedirs(s)
            
        if pacman:
            from utils import pre_charge
            pre_charge(cif_file=cifs_floder+s+".cif")
            s=s+"_pacman"

        super_cell=unit_cell(cif_file=cifs_floder+s+".cif",cutoff =14)

        for t in T:
            t = t.replace("\n","")
            if not os.path.exists(s+"/"+t):
                os.makedirs(s+"/"+t)

            for p in P:
                p = p.replace("\n","")

                if have_water:
                    
                    path_folder = s+"/"+t+"/"+p+"/"

                    if not os.path.exists(path_folder):
                        os.makedirs(path_folder, exist_ok=True)

                    shutil.copy("force_field.json",path_folder)

                    if fix_atomnumber:
                        lines = fix_atom_site_labels(cifs_floder+s)
                        with open(path_folder+s+".cif",'w') as f:
                            f.write(''.join(lines))
                    else:
                        shutil.copy(cifs_floder+s+".cif",path_folder)
                    
                    M =  get_mol_frac(tot_p=p, RH=float(RH), others=others)

                    input ={}
                    input["SimulationType"]="MonteCarlo"
                    input["NumberOfCycles"]=1000000
                    input["NumberOfInitializationCycles"]=1000000
                    input["PrintEvery"]=10000
                    
                    input["ForceField"]="."
                    
                    input["Systems"]=[]
                    input["Components"]=[]
                    
                    systems_paras={}
                    systems_paras["Type"]="Framework"                     
                    systems_paras["Name"]=s                          # structure
                    systems_paras["CutOff"]=14
                    systems_paras["NumberOfUnitCells"]=super_cell    # unitcell
                    systems_paras["ChargeMethod"]="Ewald"
                    systems_paras["EwaldPrecision"]=1e-6
                    systems_paras["UseChargesFromCIFFile"]=True
                    systems_paras["ExternalTemperature"]=float(t)           # temperature
                    systems_paras["ExternalPressure"]=float(p)              # pressure
                    systems_paras["OutputPDBMovie"]=False              

                    input["Systems"].append(systems_paras)

                    for i in range(len(G)):
                        
                        shutil.copy(G[i].replace("\n","")+".json",path_folder)

                        component_i={}
                        component_i["Name"]=G[i].replace("\n","")          # gas
                        component_i["MoleculeDefinition"]="."
                        component_i["MolFraction"] = round(float(M[i]), 8)         # mol fraction

                        # component_i["FugacityCoefficient"]=1.0
                        # component_i["TranslationProbability"]=0.5
                        # component_i["RotationProbability"]=0.5
                        # component_i["ReinsertionProbability"]=0.5
                        # component_i["SwapProbability"]=1.0
                        
                        component_i["TranslationProbability"]=1.0
                        component_i["RotationProbability"]=1.0
                        component_i["ReinsertionProbability"]=1.0
                        component_i["SwapProbability"]=1.0
                        component_i["IdentityChangeProbability"]=1.0
                        component_i["NumberOfIdentityChanges"]=int(len(G))
                        component_i["IdentityChangesList"]= [n for n in range(len(G))]  # index of gas list
                        component_i["CreateNumberofMolecules"]=0

                        input["Components"].append(component_i)

                    with open(os.path.join(path_folder, "simulation.json"),"w") as sl_file:
                        json.dump(input,sl_file,indent=2)
                    sl_file.close()
                else:

                    path_folder = s+"/"+t+"/"+p+"/"

                    if not os.path.exists(path_folder):
                        os.makedirs(path_folder, exist_ok=True)
                        
                    shutil.copy("force_field.json",path_folder)

                    if fix_atomnumber:
                        lines = fix_atom_site_labels(cifs_floder+s)
                        with open(path_folder+s+".cif",'w') as f:
                            f.write(''.join(lines))
                    else:
                        shutil.copy(cifs_floder+s+".cif",path_folder)
                        
                    input ={}
                    input["SimulationType"]="MonteCarlo"
                    input["NumberOfCycles"]=1000000
                    input["NumberOfInitializationCycles"]=1000000

                    input["ForceField"]="."
                    input["RemoveAtomNumberCodeFromLabel"]=True

                    input["Systems"]=[]
                    input["Components"]=[]
                    
                    systems_paras={}
                    systems_paras["Type"]="Framework"                     
                    systems_paras["Name"]=s                          # structure
                    systems_paras["CutOff"]=14
                    systems_paras["NumberOfUnitCells"]=super_cell    # unitcell
                    systems_paras["ChargeMethod"]="Ewald"
                    systems_paras["EwaldPrecision"]=1e-6
                    systems_paras["UseChargesFromCIFFile"]=True
                    systems_paras["ExternalTemperature"]=float(t)           # temperature
                    systems_paras["ExternalPressure"]=float(p)              # pressure
                    systems_paras["OutputPDBMovie"]=False              

                    input["Systems"].append(systems_paras)

                    with open(mol_frac_list) as mol:
                        M =  mol.readlines()

                    for i in range(len(G)):
                        
                        shutil.copy(G[i].replace("\n","")+".json",path_folder)

                        component_i={}
                        component_i["Name"]=G[i].replace("\n","")                   # gas
                        component_i["MoleculeDefinition"]="."
                        component_i["MolFraction"] = round(float(M[i]), 8)         # mol fraction

                        # component_i["FugacityCoefficient"]=1.0
                        # component_i["TranslationProbability"]=0.5
                        # component_i["RotationProbability"]=0.5
                        # component_i["ReinsertionProbability"]=0.5
                        # component_i["SwapProbability"]=1.0

                        component_i["TranslationProbability"]=1.0
                        component_i["RotationProbability"]=1.0
                        component_i["ReinsertionProbability"]=1.0
                        component_i["SwapProbability"]=1.0
                        component_i["IdentityChangeProbability"]=1.0
                        component_i["IdentityChangeProbability"]=1.0
                        component_i["NumberOfIdentityChanges"]=int(len(G))
                        component_i["IdentityChangesList"]= [n for n in range(len(G))]  # index of gas list
                        component_i["CreateNumberofMolecules"]=0

                        input["Components"].append(component_i)

                    with open(os.path.join(path_folder, "simulation.json"),"w") as sl_file:
                        json.dump(input,sl_file,indent=2) # 
                    sl_file.close()

if __name__ == "__main__":
    main()
