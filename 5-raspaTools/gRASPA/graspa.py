import os
import re
import glob
import shutil
import subprocess
import logging
import pandas as pd
from utils import unit_cell, predict_charge, load_lists

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