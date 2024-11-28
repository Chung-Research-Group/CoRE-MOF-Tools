import os
import json
import shutil
import argparse
from utils import unit_cell,fix_atom_site_labels


def main():

    parser = argparse.ArgumentParser(description="generate input files for RASPA3 software")
    parser.add_argument('--cifs_folder', default="./", type=str, help='Folder of cifs')
    parser.add_argument('--structure_list', default="list", type=str, help='structure list used for you simulation')
    parser.add_argument('--temperature_list', default="t_list", type=str, help='temperature list used for you simulation')
    parser.add_argument('--pressure_list', default="p_list", type=str, help='pressure list used for you simulation')
    parser.add_argument('--gas_type', default="H2O", type=str, help='gas used for you simulation')
    parser.add_argument('--fix_atomnumber', default=False, type=bool, help='remove the number from atom lable')
    parser.add_argument('--pacman', default=False, type=bool, help='get atomic charges for you framework')

    args = parser.parse_args()

    cifs_floder = args.cifs_folder
    structure_list = args.structure_list
    temperature_list = args.temperature_list
    pressure_list = args.pressure_list
    gas_type = args.gas_type
    pacman = args.pacman
    fix_atomnumber = args.fix_atomnumber

    with open(structure_list) as structures:
        S =  structures.readlines()

    with open(pressure_list) as pressures:
        P =  pressures.readlines()

    with open(temperature_list) as temperatures:
        T =  temperatures.readlines()
    

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
                systems_paras["ExternalTemperature"]=float(t)    # temperature
                systems_paras["ExternalPressure"]=float(p)       # pressure
                systems_paras["OutputPDBMovie"]=False              

                input["Systems"].append(systems_paras)

                shutil.copy(gas_type+".json",path_folder)

                component={}
                component["Name"]=gas_type                       # gas
                component["MoleculeDefinition"]="."
                component["TranslationProbability"]=1.0
                component["RotationProbability"]=1.0
                component["ReinsertionProbability"]=1.0
                component["SwapProbability"]=1.0
                component["CreateNuMberofmolecules"]=0

                input["Components"].append(component)

                with open(os.path.join(path_folder, "simulation.json"),"w") as sl_file:
                    json.dump(input,sl_file,indent=2) # 
                sl_file.close()

if __name__ == "__main__":
    main()
