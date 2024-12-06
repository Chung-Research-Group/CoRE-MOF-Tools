Note that the probes used for the Zeo++ computed features we provide in the CoRE MOF 2024 database tables are 1.655 angstroms (for ASA computations, 0 for the others).The probes for zeopp.py are the ones that are used in the training of the stability model in 4-MLFeatures.

## Usage                                
                                   
```
import oms, zeopp,rac, topology,others

topology.run(input_folder="./example/",output_folder="./example/topology/")
oms.run(input_folder="./example/")
others.run(input_folder="./example/",output_folder="./example/others/")
zeopp.run(input_folder="./example/",output_folder="./example/zeo/",zeo_path="/mnt/d/zeo++-0.3/")
rac.run(input_folder="./example/",output_folder="./example/RACs/")
```

