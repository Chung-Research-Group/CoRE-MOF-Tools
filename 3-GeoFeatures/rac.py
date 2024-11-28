import os, glob
import pandas as pd
from molSimplify.Informatics.MOF.MOF_descriptors import get_MOF_descriptors

class run:
    def __init__(self, input_folder, output_folder):
        self.input_folder = input_folder
        self.input_files = glob.glob(os.path.join(input_folder, "*.cif"))
        self.output = output_folder
        self.csv_path = os.path.join(output_folder, "rac_log.csv")
        self.results = []
        self.process()

    def process(self):
        os.makedirs(self.output, exist_ok=True)
        for mof in self.input_files:
            name = os.path.basename(mof).replace(".cif", "")
            try:
                full_names, full_descriptors = get_MOF_descriptors(
                    mof,
                    3,
                    path=self.output,
                    xyzpath=f'{self.output}/{name}.xyz',
                    max_num_atoms=6000
                )
                
                descriptor_data = dict(zip(full_names, full_descriptors))
                descriptor_data['name'] = name
                self.results.append(descriptor_data)

            except Exception as e:
                print(f"Failed to process {name}: {e}")

        if self.results:
            final_df = pd.DataFrame(self.results)
            final_df.to_csv(self.csv_path, index=False)
