import juliacall
import glob,os
import pandas as pd
from tqdm import tqdm

jl = juliacall.newmodule("topo")
jl.seval("using CrystalNets")


def check_unique_topology(result):
    singlenodes = result[jl.Clustering.SingleNodes]
    allnodes = result[jl.Clustering.AllNodes]
    if singlenodes != allnodes:
        raise Exception("SingleNodes result "+str(singlenodes)+" != AllNodes result "+str(allnodes))
    return (jl.ndims(singlenodes.genome), str(singlenodes))
def identify_topology(cif):
    options = jl.CrystalNets.Options(structure=jl.StructureType.MOF)
    result = jl.determine_topology(cif, options)
    return [check_unique_topology(x[0]) for x in result]

def run(input_folder, saveto: str="topo.csv")-> pd.DataFrame:
    dim_topo = []
    input_files = glob.glob(os.path.join(input_folder, "*.cif"))
    for name in tqdm(input_files):
        try:
            print(identify_topology(name))
            dim_topo.append([name,identify_topology(name)])
        except:
            dim_topo.append([name,"unknown"])


    treat_topo = []
    for i in dim_topo:
        if i[1] == "unknown":
            print(i[0], "unknown")
            treat_topo.append((i[0], "unknown"))
        elif i[1][0][1].split(" ")[0] == 'UNKNOWN':
            print(i[0], "unknown")
            treat_topo.append((i[0], "unknown"))
        elif i[1][0][1].split(" ")[0] == 'unstable':
            print(i[0], "unstable")
            treat_topo.append((i[0], "unstable"))
        elif i[1][0][1].split(" ")[0] == 'FAILED':
            print(i[0], "FAILED")
            treat_topo.append((i[0], "FAILED"))
        else:
            print(i[0], i[1][0][1])
            treat_topo.append((i[0], i[1][0][1]))

    df_treat_topo = pd.DataFrame(treat_topo, columns=["cif", "topologies"])
    df_treat_topo.to_csv(saveto, index=False)