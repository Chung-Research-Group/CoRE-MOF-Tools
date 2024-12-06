import os, glob, subprocess
import numpy as np
import pandas as pd


class run:
    def __init__(self, input_folder, output_folder, zeo_path):
        self.input_folder = input_folder
        self.input_files = glob.glob(os.path.join(input_folder, "*.cif"))
        self.output = output_folder
        self.csv_path = os.path.join(output_folder, "zeopp_log.csv")
        self.zeo_path = zeo_path
        self.results = []
        self.process()

    def process(self):
        os.makedirs(self.output, exist_ok=True)
        for mof in self.input_files:
            name = os.path.basename(mof).replace(".cif", "")
            try:
                self.run_zeopp_commands(mof, name)

                geo_dict = self.parse_geometric_descriptors(mof, name)
                if geo_dict:
                    self.results.append(geo_dict)

            except Exception as e:
                print(f"Failed to process {name}: {e}")

        if self.results:
            final_df = pd.DataFrame(self.results)
            final_df.to_csv(self.csv_path, index=False)

    def run_zeopp_commands(self, mof, name):
        commands = [
            f'{self.zeo_path}/network -ha -res {self.output}/{name}_pd.txt {mof}',
            f'{self.zeo_path}/network -sa 1.4 1.4 10000 {self.output}/{name}_sa_1_4.txt {mof}',
            f'{self.zeo_path}/network -volpo 1.4 1.4 10000 {self.output}/{name}_pov_1_4.txt {mof}',
            f'{self.zeo_path}/network -sa 1.86 1.86 10000 {self.output}/{name}_sa_1_86.txt {mof}',
            f'{self.zeo_path}/network -volpo 1.86 1.86 10000 {self.output}/{name}_pov_1_86.txt {mof}',
        ]
        for cmd in commands:
            subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

    def parse_geometric_descriptors(self, mof, name):
        try:
            geo_dict = {
                'name': name,
                'cif_file': mof,
                'Di': np.nan,
                'Df': np.nan,
                'Dif': np.nan,
                'unit_cell_volume_1_4': np.nan,
                'VSA_1_4': np.nan,
                'GSA_1_4': np.nan,
                'VPOV_1_4': np.nan,
                'GPOV_1_4': np.nan,
                'density_1_4': np.nan,
                'unit_cell_volume_1_86': np.nan,
                'VSA_1_86': np.nan,
                'GSA_1_86': np.nan,
                'VPOV_1_86': np.nan,
                'GPOV_1_86': np.nan,
                'density_1_86': np.nan,
            }

            pd_file = os.path.join(self.output, f"{name}_pd.txt")
            sa_file_1_4 = os.path.join(self.output, f"{name}_sa_1_4.txt")
            pov_file_1_4 = os.path.join(self.output, f"{name}_pov_1_4.txt")
            sa_file_1_86 = os.path.join(self.output, f"{name}_sa_1_86.txt")
            pov_file_1_86 = os.path.join(self.output, f"{name}_pov_1_86.txt")

            if not all(map(os.path.exists, [pd_file, sa_file_1_4, pov_file_1_4, sa_file_1_86, pov_file_1_86])):
                print(f"One or more required files are missing for {name}")
                return None
            
            with open(pd_file) as f:
                line = f.readline().split()
                geo_dict['Di'], geo_dict['Df'], geo_dict['Dif'] = map(float, line[1:4])

            geo_dict['unit_cell_volume_1_4'], geo_dict['density_1_4'], geo_dict['VSA_1_4'], geo_dict['GSA_1_4'] = self.sa_process(sa_file_1_4)
            geo_dict['density_1_4'], geo_dict['POAV_1_4'], geo_dict['PONAV_1_4'], geo_dict['GPOAV_1_4'], geo_dict['GPONAV_1_4'], \
            geo_dict['POAV_vol_frac_1_4'], geo_dict['PONAV_vol_frac_1_4'], geo_dict['VPOV_1_4'], geo_dict['GPOV_1_4'] = self.pov_process(pov_file_1_4)

            geo_dict['unit_cell_volume_1_86'], geo_dict['density_1_86'], geo_dict['VSA_1_86'], geo_dict['GSA_1_86'] = self.sa_process(sa_file_1_86)
            geo_dict['density_1_86'], geo_dict['POAV_1_86'], geo_dict['PONAV_1_86'], geo_dict['GPOAV_1_86'], geo_dict['GPONAV_1_86'], \
            geo_dict['POAV_vol_frac_1_86'], geo_dict['PONAV_vol_frac_1_86'], geo_dict['VPOV_1_86'], geo_dict['GPOV_1_86'] = self.pov_process(pov_file_1_86)

            return geo_dict

        except Exception as e:
            print(f"Error parsing geometric descriptors for {name}: {e}")
            return None

    def sa_process(self, file_path):
        with open(file_path) as f:
            for i, row in enumerate(f):
                if i == 0:
                    unit_cell_volume = float(row.split('Unitcell_volume:')[1].split()[0])
                    crystal_density = float(row.split('Density:')[1].split()[0])
                    VSA = float(row.split('ASA_m^2/cm^3:')[1].split()[0])
                    GSA = float(row.split('ASA_m^2/g:')[1].split()[0])
        return unit_cell_volume, crystal_density, VSA, GSA

    def pov_process(self, file_path):
        with open(file_path) as f:
            for i, row in enumerate(f):
                if i == 0:
                    density = float(row.split('Density:')[1].split()[0])
                    POAV = float(row.split('POAV_A^3:')[1].split()[0])
                    PONAV = float(row.split('PONAV_A^3:')[1].split()[0])
                    GPOAV = float(row.split('POAV_cm^3/g:')[1].split()[0])
                    GPONAV = float(row.split('PONAV_cm^3/g:')[1].split()[0])
                    POAV_volume_fraction = float(row.split('POAV_Volume_fraction:')[1].split()[0])
                    PONAV_volume_fraction = float(row.split('PONAV_Volume_fraction:')[1].split()[0])
                    VPOV = POAV_volume_fraction + PONAV_volume_fraction
                    GPOV = VPOV / density
        return density, POAV, PONAV, GPOAV, GPONAV, POAV_volume_fraction, PONAV_volume_fraction, VPOV, GPOV