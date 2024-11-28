from PACMANCharge import pmcharge
import os, glob
import pandas as pd


class run:
    def __init__(self, input_folder, output_folder):
        self.input_folder = input_folder
        self.input_files = glob.glob(os.path.join(input_folder, "*.cif"))
        self.output = output_folder
        self.csv_path = os.path.join(output_folder, "pacman_log.csv")
        self.results = []
        self.process()

    def process(self):
        os.makedirs(self.output, exist_ok=True)
        for mof in self.input_files:
            name = os.path.basename(mof).replace(".cif", "")
            try:
                pmcharge.predict(
                    cif_file=mof,
                    charge_type="DDEC6",
                    digits=10,
                    atom_type=True,
                    neutral=True,
                    keep_connect=True
                )
                
                pbe, bandgap = pmcharge.Energy(cif_file=mof)
                
                self.results.append({
                    'name': name,
                    'pbe': pbe,
                    'bandgap': bandgap
                })

            except Exception as e:
                print(f"Failed to process {name}: {e}")

        if self.results:
            final_df = pd.DataFrame(self.results)
            final_df.to_csv(self.csv_path, index=False)
