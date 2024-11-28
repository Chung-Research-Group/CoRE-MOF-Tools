# HTS for RASPA

## RASPA2

* Install by Conda
```sh
conda install conda-forge::raspa2                
```
*  Usage
  
```sh
from raspa2 import Generate, Check, Result

Generate(task="Widom",use_pacman=True,
         s_list="s_list",t_list="t_list", #,p_list="p_list",g_list="g_list"
         template_file="./data/template_RASPA2/2_Widom.input",
         def_folder="./data/template_RASPA2/",
         cif_folder="./data/structures/",
         job_script="./data/template_RASPA2/raspa.slurm")

Check(task="Widom",use_pacman=True,s_list="s_list",t_list="t_list")
Result(task="Widom",use_pacman=True, s_list="s_list",t_list="t_list",out_folder="results") #unit="molecules/unit"
```
  
## gRASPA                      
            
It will be done soon.

## RASPA3

It will be done soon.                            

