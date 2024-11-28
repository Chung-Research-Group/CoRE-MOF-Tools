import json, yaml
import pandas as pd
import argparse

def main():

    parser = argparse.ArgumentParser(description="generate force field files for RASPA3 software")
    parser.add_argument('--mol_csv', nargs='+', default="./CH4.csv", type=str, help='define molecule CSV files')
    parser.add_argument('--force_field_csv', default="./FF_para.csv", type=str, help='parameter such LJ')
    args = parser.parse_args()

    for i in range(len(args.mol_csv)):
        gas = str(args.mol_csv[i].split("/")[-1].replace(".csv",""))
        info=pd.read_csv(args.mol_csv[i])
        gas_def = {}
        
        gas_def["CriticalTemperature"]=float(info["CriticalTemperature"].dropna().values[0])
        gas_def["CriticalPressure"]=float(info["CriticalPressure"].dropna().values[0])
        gas_def["AcentricFactor"]=float(info["AcentricFactor"].dropna().values[0])
        gas_def["Type"]=str(info["Type"][0])
        
        gas_def["pseudoAtoms"]=[]
        for i in range(len(info["atom"].dropna().values)):
            gas_def["pseudoAtoms"].append([info["atom"][i],
                                        [float(info["posx"][i]),
                                            float(info["posy"][i]),
                                            float(info["posz"][i])]])

        if len(info["Bonds"].dropna().values)>0:
            gas_def["Bonds"]=[]
            for j in range(len(info["Bonds"].dropna().values)):
                gas_def["Bonds"].append([int(info["Bonds"][j]),int(info["Bonds0"][j])])
        else:
            pass

        with open(gas+".json","w") as mol_file:
            json.dump(gas_def,mol_file,indent=2)
        mol_file.close()



    ff = {}
    ff["PseudoAtoms"]=[]
    ff["SelfInteractions"]=[]

    para=pd.read_csv(args.force_field_csv)

    for i,name in enumerate(para["name"]):
        paras={}
        paras["name"]=name
        paras["print_to_output"]=True if para["print_to_output"][i]=="yes" else False
        paras["element"]=para["element"][i]
        paras["print_as"]=para["print_as"][i]
        paras["mass"]=para["mass"][i]
        paras["charge"]=para["charge"][i]
        ff["PseudoAtoms"].append(paras)

    for i,name in enumerate(para["name"]):
        paras={}
        paras["name"]=name
        paras["type"]=para["type"][i]
        paras["parameters"]=[para["parameter1"][i],para["parameter2"][i]]
        paras["source"]=para["source"][i]
        ff["SelfInteractions"].append(paras)

    ff["MixingRule"]="Lorentz-Berthelot"
    ff["TruncationMethod"]="shifted"
    ff["TailCorrections"]=False

    with open("force_field.json","w") as FF_file:
        json.dump(ff,FF_file,indent=2)
    FF_file.close()


if __name__ == "__main__":
    main()