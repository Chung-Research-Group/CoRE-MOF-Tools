import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description="combine single result for RASPA3 software")
    parser.add_argument('--structure_list', type=str, default="list", help='structure list used for your simulation')
    parser.add_argument('--input_csv', type=str, default="_CO2_N2_298", help='the output CSV file name')
    parser.add_argument('--output_csv', type=str, default="result_400ppm_CO2_N2", help='the output CSV file name')
    args = parser.parse_args()


    with open(args.structure_list) as file:
        structures = [line.strip() for line in file]

    all_data = pd.DataFrame()

    for structure in structures:
        s = structure.replace("\n","")
        df = pd.read_csv(s+args.input_csv+".csv")
        df['name'] = s

        all_data = pd.concat([all_data, df], ignore_index=True)

    all_data.to_csv(args.output_csv+".csv",index=False)

if __name__ == "__main__":
    main()