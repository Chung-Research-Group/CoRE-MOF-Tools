import os, argparse
import pandas as pd

def find_data(filename, unit_marker, gas_marker, uptake_marker):

    start_found = False
    second_found = False

    with open(filename, 'r') as file:
        lines = file.readlines()

    for line in lines:
        if not start_found and unit_marker in line:
            start_found = True
        elif start_found and not second_found and gas_marker in line:
            second_found = True
        elif second_found and uptake_marker in line:
            return line.strip()

    return None

def main():
    parser = argparse.ArgumentParser(description="get results for RASPA3 software")
    parser.add_argument('--RH', type=str, default="rh_list", help='pressure list used for your simulation')
    parser.add_argument('--unit', type=str, default="mol/kg", help='unit for block averages')
    parser.add_argument('--output_csv', type=str, default="H2O_CO2_N2", help='the output CSV file name')
    args = parser.parse_args()

    with open(args.RH) as file:
        rh_list = [line.strip() for line in file]

    block_averages_keyword = "BLOCK AVERAGES (LOADING: " + args.unit + ")"
    overall_keyword = "Overall: Average:"
    columns = ["RH","H2O","H2O_error","CO2","CO2_error","N2","N2_error"]

    results = []
    for rh in rh_list:
        single_data = [rh]
        path_file = os.path.join(rh, "result.txt")
        if not os.path.isfile(path_file):
            continue
        
        for gas in ["H2O","CO2","N2"]:
            result = find_data(path_file, block_averages_keyword, gas, overall_keyword)
            try:
                single_data.append(float(result.split(" ")[2].replace(",","")))
                single_data.append(float(result.split(" ")[4]))
            except:
                single_data.append(0)
                single_data.append(0)
            
        results.append(single_data)
    print(results)
    df = pd.DataFrame(results,columns=columns)
    df.to_csv("result/" + args.output_csv + ".csv", index=False)

if __name__ == "__main__":
    main()