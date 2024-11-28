import os, glob
import pandas as pd
from .cp_app.descriptors import cv_features
from .cp_app.featurizer import featurize_structure
from .cp_app.predictions import predict_Cv_ensemble_structure_multitemperatures


class run:
    def __init__(self, input_folder, output_folder):
        self.input_folder = input_folder
        self.input_files = glob.glob(os.path.join(input_folder, "*.cif"))
        self.output = output_folder
        self.csv_path = os.path.join(output_folder, "cp_log.csv")
        self.process()

    def process(self):
        os.makedirs(self.output, exist_ok=True)
        for mof in self.input_files:
            name = os.path.basename(mof).replace(".cif", "")
            try:
                self.pre_cp(mof, self.output)
            except Exception as e:
                print(f"Failed to process {name}: {e}")

        all_cp_files = glob.glob(os.path.join(self.output, "*_cp.csv"))
        if all_cp_files:
            combined_df = pd.concat([pd.read_csv(f) for f in all_cp_files], ignore_index=True)
            combined_df.to_csv(self.csv_path, index=False)
        combined_df = pd.read_csv(self.csv_path)
        mass_columns = [col for col in combined_df.columns if "Cv_gravimetric" in col]
        combined_df[mass_columns] = combined_df[mass_columns] * 1000
        new_column_names = {col: col.replace("Cv_gravimetric", "Cv_gravimetric_kg") for col in mass_columns}
        combined_df.rename(columns=new_column_names, inplace=True)
        combined_df.to_csv(self.csv_path, index=False)


    def pre_cp(self, cif_file, save_folder):
        name = cif_file.split("/")[-1].split(".")[0]

        featurize_structure(cif_file, verbos=False, saveto=os.path.join(save_folder, f"{name}_features.csv"))
        current_path = os.getcwd()
        predict_Cv_ensemble_structure_multitemperatures(
            path_to_models=current_path+"/features/ML/cp_app/ensemble_models_smallML_120_100",
            structure_name=name + ".cif",
            features_file=os.path.join(save_folder, f"{name}_features.csv"), 
            FEATURES=cv_features,
            temperatures=[300, 350, 400],
            save_to=os.path.join(save_folder, f"{name}_cp.csv")
        )