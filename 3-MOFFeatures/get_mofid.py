import os
from mofid.run_mofid import cif2mofid
import re
import pandas as pd
import numpy as np

def get_mofid(root_dir, save_csv):
    
    files_name = os.listdir(root_dir)
    all_names = []
    for cif_name in files_name:
        if os.path.splitext(cif_name)[1] == '.cif':
            all_names.append(os.path.splitext(cif_name)[0])
        else:
            pass
    
    all_id = []
    for name in all_names:
        cif_path = root_dir + name + '.cif'
        mofid = cif2mofid(cif_path)['mofid']
        
        
        id_mof = re.split(r"[;]+", mofid)[0]
        all_id.append([name,id_mof])
        
    df_all_id = pd.DataFrame(all_id, columns = ["name","mofid"])
    
    df_all_id.to_csv(save_csv,index=False)
    
    return all_id

def get_mofkey(root_dir, save_csv):
    
    files_name = os.listdir(root_dir)
    all_names = []
    for cif_name in files_name:
        if os.path.splitext(cif_name)[1] == '.cif':
            all_names.append(os.path.splitext(cif_name)[0])
        else:
            pass
    
    all_key = []
    for name in all_names:
        cif_path = root_dir + name + '.cif'
        mofkey = cif2mofid(cif_path)['mofkey']
        
        key_mof = re.split(r"[;]+", mofkey)[0]
        all_key.append([name,key_mof])
        np.savetxt("mofkey.txt",all_key, fmt='%s')
    df_all_id = pd.DataFrame(all_key, columns = ["name","mofkey"])
    
    df_all_id.to_csv(save_csv,index=False)
    
    return all_key