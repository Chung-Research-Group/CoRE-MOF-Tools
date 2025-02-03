## Usage                               
                                    
```
import stability, pacman, cp

pacman.run(input_folder="./example/",output_folder="./example/pacman/")
cp.run(input_folder="./example/",output_folder="./example/heat_capacity/")
stability.run(zeo_input="./example/zeopp_log.csv", rac_input="./example/rac_log.csv", output_folder="./example/stability/")
```

final_model_flag_few_epochs.h5 is the solvent removal stability model.

final_model_T_few_epochs.h5 is the thermal stability model, pulled from MOFSimplify.

Models are from /Users/gianmarcoterrones/Research/MOFs/water_stability_work/current_work/CoRE_data/model/burtch-labels_v2/current_model/combined_data_sets/combined_train/with_hyperopt/random_forest_RACs_and_Zeo++/models
	seed 9 used for train/test splits

