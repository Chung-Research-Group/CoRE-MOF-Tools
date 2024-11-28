import os
import re
import argparse

def main():
    parser = argparse.ArgumentParser(description="Generate input files for RASPA3 software")
    parser.add_argument('--structure_list', default="list", type=str, help='Structure list used for your simulation')
    parser.add_argument('--temperature_list', default="t_list", type=str, help='Temperature list used for your simulation')
    parser.add_argument('--pressure_list', default="p_list", type=str, help='Pressure list used for your simulation')
    args = parser.parse_args()

    structure_list = args.structure_list
    temperature_list = args.temperature_list
    pressure_list = args.pressure_list

    with open(structure_list) as structures:
        S = [s.strip() for s in structures.readlines()]
    with open(pressure_list) as pressures:
        P = [p.strip() for p in pressures.readlines()]
    with open(temperature_list) as temperatures:
        T = [t.strip() for t in temperatures.readlines()]

    total_jobs = len(S) * len(T) * len(P)
    i = 0
    j = 0
    k = 0
    l = 0

    header_format = "{:<15} {:<10} {:<10} {:<10} {:<15} {:<15}"
    row_format = "{:<15} {:<10} {:<10} {:<10} {:<15} {:<15}"

    print(header_format.format("Status", "CompletedCycles", "TotalCycles", "Structure", "Temperature", "Pressure"))
    print("-" * 75)

    for s in S:
        for t in T:
            for p in P:
                path = os.path.join(s, t, p, "output")
                key_word = "Current cycle:"
                key_check_comp = 'Loadings'
                finish = False
                last_line = None

                if os.path.isdir(path):
                    for f in os.listdir(path):
                        if f.endswith(".data"):
                            with open(os.path.join(path, f), "r") as file:
                                for line in file:
                                    if key_check_comp in line:
                                        finish = True
                                        i += 1
                                        print(row_format.format("Finish", "-", "-", s, f"{t} K", f"{p} Pa"))
                                        break

                                if not finish:
                                    file.seek(0)
                                    for line in file:
                                        if key_word in line:
                                            last_line = line

                                    if last_line:
                                        info = re.split(r'\s+', last_line)
                                        status = "Initialization" if info[0] == "Initialization:" else "Production"
                                        finished = info[3] if status == "Initialization" else info[2]
                                        total = info[6] if status == "Initialization" else info[5]
                                        if status == "Initialization":
                                            k += 1
                                        else:
                                            j += 1
                                        print(row_format.format(status, finished, total, s, f"{t} K", f"{p} Pa"))
                                    else:
                                        print(row_format.format("No Data", '-', '-', s, f"{t} K", f"{p} Pa"))
                else:
                    l += 1
                    print(row_format.format("Not started", '-', '-', s, f"{t} K", f"{p} Pa"))

    print("-" * 75)
    print("Finish / Total:",i,"/",total_jobs)
    print("On Initialization:",k)
    print("On Production:",j)
    print("Not Start:",l)
if __name__ == "__main__":
    main()
