import os
import re
import glob
import shutil
import subprocess
import logging
import pandas as pd
from GZTOOLs.utils import unit_cell, predict_charge, load_lists

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def submit_job(job_script, folder):
    try:
        cwd_dir = os.getcwd()
        os.chdir(folder)
        result = subprocess.run(["sbatch", job_script], capture_output=True, text=True, check=True)
        os.chdir(cwd_dir)
        output = result.stdout.strip()
        logging.info(f"Job submission successful: {output}")
        job_id = output.split()[-1]
        return job_id
    except subprocess.CalledProcessError as e:
        logging.error(f"Error submitting job in folder {folder}: {e.stderr}")
        return None
    finally:
        os.chdir(cwd_dir)


class Generate:
    def __init__(self, task, use_pacman, **kwargs):
        self.task = task
        self.use_pacman = use_pacman
        self.kwargs = kwargs
        self.validate_inputs()
        self.run_task()

    def validate_inputs(self):
        required_keys = ["def_folder", "cif_folder", "template_file", "job_script"]
        for key in required_keys:
            if key not in self.kwargs:
                raise ValueError(f"Missing required parameter: {key}")
        logging.info("All required inputs are validated.")

    def run_task(self):
        if hasattr(self, self.task):
            logging.info(f"Starting task: {self.task}")
            getattr(self, self.task)()
        else:
            raise ValueError(f"Unsupported task: {self.task}")

    def prepare_files(self, s, t=None, p=None, g=None):
        def_folder = self.kwargs["def_folder"]
        cif_folder = self.kwargs["cif_folder"]
        template_file = self.kwargs["template_file"]
        job_script = self.kwargs["job_script"]
        folder_parts = [s, g, t, p]
        folder_parts = [str(part) for part in folder_parts if part is not None]
        folder = os.path.join(*folder_parts)
        os.makedirs(folder, exist_ok=True)
        for file in glob.glob(os.path.join(def_folder, "*def")):
            shutil.copy(file, folder)
        shutil.copy(os.path.join(cif_folder, f"{s}.cif"), folder)
        with open(template_file, "r") as f:
            input_template = f.read()
        cutoff = self.extract_cutoff(input_template)
        UC = unit_cell(os.path.join(cif_folder, s), cutoff)
        input_file = [
            f"FrameworkName\t{s}" if "FrameworkName" in line else
            f"ExternalTemperature\t{t}" if "ExternalTemperature" in line else
            f"ExternalPressure\t{p}" if "ExternalPressure" in line else
            f"UnitCells\t{UC}" if "UnitCells" in line else
            f"Component 0 MoleculeName\t{g}" if g and "Component" in line else
            line
            for line in input_template.splitlines()
        ]
        with open(os.path.join(folder, "simulation.input"), "w") as f:
            f.write("\n".join(input_file))

        with open(job_script, "r") as f:
            job_template = f.read()
        job_file = job_template.replace("my_job", f"{s}_{g}_{t}_{p}")
        with open(os.path.join(folder, "submit.sh"), "w") as f:
            f.write(job_file)
        return folder

    @staticmethod
    def extract_cutoff(template):
        try:
            return float(re.search(r"CutOff\s+(\d+)", template).group(1))
        except AttributeError:
            raise ValueError("CutOff value not found in template file.")

    def SingleAdsorption(self):
        S, G, T, P = load_lists("s_list", "g_list", "t_list", "p_list", kwargs=self.kwargs)
        for s in S:
            if self.use_pacman:
                predict_charge(os.path.join(self.kwargs["cif_folder"], s))
                s += "_pacman"
            for g in G:
                for t in T:
                    for p in P:
                        folder = self.prepare_files(s, t=t, p=p, g=g)
                        submit_job("submit.sh", folder)

    def MixtureAdsorption(self):
        S, T, P = load_lists("s_list", "t_list", "p_list", kwargs=self.kwargs)
        for s in S:
            if self.use_pacman:
                predict_charge(os.path.join(self.kwargs["cif_folder"], s))
                s += "_pacman"
            for t in T:
                for p in P:
                    folder = self.prepare_files(s, t=t, p=p)
                    submit_job("submit.sh", folder)


class Result:
    def __init__(self, task, use_pacman, out_folder, **kwargs):
        self.task = task
        self.use_pacman = use_pacman
        self.out_folder = out_folder
        self.kwargs = kwargs
        self.s_list = kwargs.get("s_list")
        self.p_list = kwargs.get("p_list")
        self.t_list = kwargs.get("t_list")
        self.g_list = kwargs.get("g_list") if task in ["SingleAdsorption"] else None
        if self.task == "SingleAdsorption":
            self._uptake_single()
        elif self.task == "MixtureAdsorption":
            self._uptake_mixture()
        if self.task == "Widom":
            self._henry_heat()

    def load_data(self):
        with open(self.s_list) as structures:
            S = [s.strip() for s in structures.readlines()]
        if hasattr(self, 'p_list') and self.p_list:
            with open(self.p_list) as pressures:
                P = [p.strip() for p in pressures.readlines()]
        else:
            P = []
        if hasattr(self, 't_list') and self.t_list:
            with open(self.t_list) as temperatures:
                T = [t.strip() for t in temperatures.readlines()]
        else:
            T = [] 
        G = []
        if hasattr(self, 'g_list') and self.g_list:
            with open(self.g_list) as gases:
                G = [g.strip() for g in gases.readlines()]
        if hasattr(self, 'use_pacman') and self.use_pacman:
            S = [f"{s}_pacman" for s in S]
        return S, P, T, G
    
    def _uptake_single(self):
        S, P, T, G = self.load_data()
        for s in S:
            for g in G:
                for t in T:
                    T_data = []
                    for p in P:
                        P_data = [p]
                        output = os.path.join(s, g, str(t), str(p), "Output", "System_0")
                        if not os.path.isdir(output):
                            continue
                        for file in os.listdir(output):
                            if file.endswith(".data"):
                                data_path = os.path.join(output, file)
                            with open(data_path, "r") as f:
                                data_lines = f.readlines()
                            key_word = "Average loading excess [" + self.kwargs["unit"] # molecules/unit, mol/kg, milligram/gram, cm^3 (STP)/gr, cm^3 (STP)/cm^3
                            data_line = [line for line in data_lines if key_word in line]
                            parts = re.split(r'\s+', data_line[0])
                            uptake = float(parts[6])
                            error = float(parts[8])
                            P_data.extend([uptake, error])
                        T_data.append(P_data)
                    pd.DataFrame(T_data, columns=["P","N","E"]).to_csv(f"{self.out_folder}/{s}_{g}_{t}.csv",index=False)
                    
    def _uptake_mixture(self):
        S, P, T, _, = self.load_data()
        for s in S:
            for t in T:
                T_data = []
                for p in P:
                    P_data = [p]
                    output = os.path.join(s, str(t), str(p), "Output", "System_0")
                    if not os.path.isdir(output):
                        continue
                    uptake_error_pairs = []
                    for file in os.listdir(output):
                        if file.endswith(".data"):
                            data_path = os.path.join(output, file)
                            with open(data_path, "r") as f:
                                data_lines = f.readlines()
                            key_word = "Average loading excess [" + self.kwargs["unit"]  # molecules/unit, mol/kg, etc.
                            data_line = [line for line in data_lines if key_word in line]
                            for i in range(len(data_line)):
                                parts = re.split(r'\s+', data_line[i])
                                uptake = float(parts[6])
                                error = float(parts[8])
                                uptake_error_pairs.append((uptake, error))
                    for pair in uptake_error_pairs:
                        P_data.extend(pair)
                    T_data.append(P_data)

                max_n = max(len(row) for row in T_data) - 1
                column_names = ["P"] + [f"{label}{i}" for i in range(max_n // 2) for label in ("N", "E")]
                df = pd.DataFrame(T_data, columns=column_names)
                df.to_csv(f"{self.out_folder}/{s}_{t}.csv", index=False)


class Check:
    def __init__(self, task, use_pacman, **kwargs):
        self.task = task
        self.use_pacman = use_pacman
        self.kwargs = kwargs
        self.s_list = kwargs.get("s_list")
        self.p_list = kwargs.get("p_list")
        self.t_list = kwargs.get("t_list")
        self.g_list = kwargs.get("g_list") if task in ["SingleAdsorption"] else None
        self.check_status()

    def load_data(self):
        with open(self.s_list) as structures:
            S = [s.strip() for s in structures.readlines()]
        if hasattr(self, 'p_list') and self.p_list:
            with open(self.p_list) as pressures:
                P = [p.strip() for p in pressures.readlines()]
        else:
            P = []
        if hasattr(self, 't_list') and self.t_list:
            with open(self.t_list) as temperatures:
                T = [t.strip() for t in temperatures.readlines()]
        else:
            T = [] 
        G = []
        if hasattr(self, 'g_list') and self.g_list:
            with open(self.g_list) as gases:
                G = [g.strip() for g in gases.readlines()]
        if hasattr(self, 'use_pacman') and self.use_pacman:
            S = [f"{s}_pacman" for s in S]
        return S, P, T, G

    def format_output(self, status, finished, total, structure, temperature, pressure, gas=None):
        row_format = "{:<15} {:<10} {:<10} {:<20} {:<15} {:<15} {:<15}"
        gas_str = gas if gas else "-"
        print(row_format.format(status, finished, total, structure, f"{temperature}", f"{pressure}", gas_str))

    def check_directory(self, path, key_check_comp, key_word):
        finish = False
        last_line = None

        if not os.path.isdir(path):
            return "Not started", "-", "-", None

        for f in os.listdir(path):
            if f.endswith(".data"):
                with open(os.path.join(path, f), "r") as file:
                    lines = file.readlines()
                    if any(key_check_comp in line for line in lines):
                        finish = True
                        return "Finished", "-", "-", None
                    for line in lines:
                        if key_word in line:
                            last_line = line
                            break

        if last_line:
            info = re.split(r'\s+', last_line)
            status = "Initialization" if info[0] == "[Init]" else "Production"
            finished = info[3] if status == "Initialization" else info[5]
            total = info[6] if status == "Initialization" else info[8]
            return status, finished, total, last_line

        return "No Data", "-", "-", None

    def check_status(self):
        S, P, T, G = self.load_data()
        total_jobs = len(S) * (len(P) if P else 1) * len(T) * (len(G) if G else 1)
        uncom_names = []

        print("{:<15} {:<10} {:<10} {:<20} {:<15} {:<15} {:<15}".format(
            "Status", "CompletedCycles", "TotalCycles", "Structure", "Temperature", "Pressure", "Gas"
        ))
        print("-" * 95)

        for s in S:
            gases = G if G else [None]
            for g in gases:
                for t in T:
                    Ps = P if P else [None]
                    for p in Ps:
                        path = os.path.join(s, g if g else "", t, p if p else "", "Output/System_0/")
                        status, finished, total, last_line = self.check_directory(
                            path, "Simulation finished", "Current cycle:"
                        )

                        if status == "Not started":
                            self.format_output(status, "-", "-", s, t, p, g)
                        elif status == "Finished":
                            self.format_output(status, "-", "-", s, t, p, g)
                        else:
                            uncom_names.append(s)
                            self.format_output(status, finished, total, s, t, p, g)
        print("-" * 95)
        print(f"Running / Total: {len([s for s in uncom_names if s])} / {total_jobs}")
        with open("running_list", "w") as unfinish:
            unfinish.write("\n".join(uncom_names))