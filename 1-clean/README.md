### preprocess.py                                      
                              
-  split multi structures from a single CIF file to different CIF files <img src="/figs/split.png" alt="Schematic for multi split" width="400" height="100">
-  check that the structure contains carbon and metal atoms to make sure it is MOF
-  make primitive cell for MOF
-  make P1, remove symmetry

### clean.py                             
                                
-  remove free solvent
-  remove Coordinating solvent
-  keep ion for charge balance
  <img src="/figs/clean.png" alt="Schematic for clean" width="600" height="300">

## Usage                       
                                                                      
```
from preprocess import run
run(input_folder="./example/origin/", output_folder="./example/out_preprocess/")

from clean import run
run(input_folder="./example/out_preprocess/", output_folder="./example/out_clean/")
```
