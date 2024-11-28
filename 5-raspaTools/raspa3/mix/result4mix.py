import os
import re
import pandas as pd
import argparse
import numpy as np


def main():
    parser = argparse.ArgumentParser(description="get results for RASPA3 software")
    parser.add_argument('--structure_list', type=str, default="list", help='structure list used for your simulation')
    parser.add_argument('--temperature_list', type=str, default="t_list", help='temperature list used for your simulation')
    parser.add_argument('--pressure_list', type=str, default="p_list", help='pressure list used for your simulation')
    parser.add_argument('--gas_list', type=str, default="g_list", help='gas list used for your simulation')
    parser.add_argument('--unit', type=str, default="mol/kg", choices=['mol/kg', 'mg/g', 'molecules'], help='the unit of uptake')
    parser.add_argument('--adsorption_type', type=str, default="absolute", choices=['absolute', 'excess'], help='the type of uptake')
    parser.add_argument('--output_csv', type=str, default="CO2_N2", help='the output CSV file name')
    args = parser.parse_args()

    with open(args.structure_list) as file:
        structures = [line.strip() for line in file]
    with open(args.gas_list) as file:
        gases = [line.strip() for line in file]
    with open(args.pressure_list) as file:
        pressures = [line.strip() for line in file]
    with open(args.temperature_list) as file:
        temperatures = [line.strip() for line in file]

    keyword = "Abs. loading average" if args.adsorption_type == "absolute" else "Excess loading average"
    unit_key = args.unit + "-framework" if args.unit != "molecules" else "molecules/cell"

    for structure in structures:
        s = structure.replace("\n","")
        for t in temperatures:
            all_data = []
            for p in pressures:
                single_data = [p]
                path_folder = os.path.join(s, t, p, "output")
                if not os.path.exists(path_folder):
                    continue
                data_found = False
                for f in os.listdir(path_folder):
                    if f.endswith(".data"):
                        with open(os.path.join(path_folder, f), "r") as file:
                            content = file.read()
                        if keyword in content:
                            with open(os.path.join(path_folder, f), "r") as file:
                                for line in file:
                                    if keyword in line and unit_key in line:
                                        uptake_values = re.split(r'\s+', line)
                                        single_data.extend([float(uptake_values[-5]), float(uptake_values[-3])])
                                        data_found = True
                        if data_found:
                            break
                if not data_found:
                    single_data.extend(np.zeros(len(gases)*2))

                all_data.append(single_data)
                # print(all_data)
            columns = ["Pressure"] + [g + suffix for g in gases for suffix in ["", "_error"]]
            df_all_data = pd.DataFrame(all_data, columns=columns)
            output_path = s+ "_" + args.output_csv + "_" + t + ".csv"
            df_all_data.to_csv(output_path, index=False)

if __name__ == "__main__":
    main()
