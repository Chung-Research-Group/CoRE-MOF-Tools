import os, sys
import pandas as pd
from tqdm import tqdm

data_file = sys.argv[1]
datasets = ["ASR","FSR","Ion"]

avogadro_number = 6.022e23

for dataset in datasets[:]:
    data = pd.read_excel(data_file, sheet_name=dataset)
    # data_info = pd.read_csv("../../CoREMOF2024/"+dataset + "_data_20241010.csv")
    for i, name in tqdm(enumerate(data["coreid"][:])):
        data_sample = data[data["coreid"] == name]
        try:
            n_unitcells = data_sample["NUC_simulation"].values[0]
            V = "{:.2e}".format(data_sample["Volume (cm3)"].values[0])
            status = data_sample["Sim_Status"].values[0]
            id_code = f"{dataset}_{i}_GEMC_{status}"
            
            pressure_raw = data_sample["Pressure (Pa)"].values[0]
            loading_raw = data_sample["Loading (Molecules/Supercell)"].values[0]

            # super_cells = data_sample["NUC_simulation"].values[0]
            # amu_unitcell = data_info[data_info["coreid"]==name]["average_atomic_mass"].values[0] * data_info[data_info["coreid"]==name]["natoms"].values[0]
            # unit_cell_mass_g = amu_unitcell * 1.66054e-2

            if not pressure_raw:
                print(f"Skipping {name} in {dataset} due to empty pressure data.")
                continue

            if isinstance(pressure_raw, str):
                pressure = pressure_raw.split(',')
            else:
                pressure = [pressure_raw]

            if isinstance(loading_raw, str):
                loading = loading_raw.split(',')
            else:
                loading = [loading_raw]
            p0 = "4540"
            aif_info = (
                f"# Gibbs Ensemble Monte Carlo Simulation Software: MCCCS‒MN (Monte Carlo for Complex Chemical Systems‒Minnesota)\n"
                f"# Setting: cutoff=9 Å, simulation cycles: 5E3 + 1.5E4 + 2.5E4\n"
                f"data_aif\n"
                f"_adsnt_material_id '{name}'\n"
                f"_exptl_temperature 298\n"
                f"_simltn_code {id_code}\n"
                f"_simltn_date 2024-11-19\n"
                f"_simltn_forcefield_adsorbent UFF4MOF\n"
                f"_simltn_forcefield_adsorptive TIP4P\n"
                f"_simltn_size 'N of unitcells: {n_unitcells}, Volume (cm3): {V}'\n"
                f"_units_loading Molecules/Supercell\n"
                f"_units_pressure Pa\n"
                f"_units_temperature K\n"
                f"loop_\n"
                f"_adsorp_pressure\n"
                f"_adsorp_p0\n"
                f"_adsorp_amount\n"
            )

            # for p, l in zip(pressure, loading):
            #     p = p.strip() if isinstance(p, str) else p
            #     l = float(l.strip()) if isinstance(l, str) else l

            #     mol_per_supercell = l / avogadro_number
            #     loading_mmol_per_g = (mol_per_supercell * 1000) / (super_cells * unit_cell_mass_g)

            #     aif_info += f"{p} {p0} {loading_mmol_per_g:.4e}\n"

            # file_path = f"water_data/{dataset}"
            # os.makedirs(file_path, exist_ok=True)
            # file_name = f"{file_path}/{name}.aif"
            # with open(file_name, "w") as file:
            #     file.write(aif_info.strip())


            for p, l in zip(pressure, loading):
                p = p.strip() if isinstance(p, str) else p
                l = l.strip() if isinstance(l, str) else l
                aif_info += f"{p} {p0} {l}\n"
            file_path = f"GEMC/{dataset}"
            os.makedirs(file_path, exist_ok=True)
            file_name = f"{file_path}/{name}.aif"
            with open(file_name, "w") as file:
                file.write(aif_info.strip())


        except KeyError as e:
            print(f"KeyError encountered in {dataset} for {name}: {e}")
        except IndexError as e:
            print(f"IndexError encountered in {dataset} for {name}: {e}")
        except Exception as e:
            print(f"An error occurred while processing {name} in {dataset}: {e}")
