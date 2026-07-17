"""Download structures and query information of CoRE MOF Database.
"""

import json
import os
from difflib import get_close_matches
from functools import lru_cache
from pathlib import Path
import tempfile
import zipfile

import requests

from gemmi import cif

package_directory = os.path.abspath(__file__).replace("structure.py","")

files_to_download = {
                    'data/CR.json': 'https://raw.githubusercontent.com/Chung-Research-Group/CoRE-MOF-Tools/main/CoREMOF/data/info/CR.json',
                    'data/NCR.json': 'https://raw.githubusercontent.com/Chung-Research-Group/CoRE-MOF-Tools/main/CoREMOF/data/info/NCR.json',
                    'data/SI/CR.zip': 'https://raw.githubusercontent.com/Chung-Research-Group/CoRE-MOF-Tools/main/CoREMOF/data/SI/CR.zip',
                    'data/SI/NCR.zip': 'https://raw.githubusercontent.com/Chung-Research-Group/CoRE-MOF-Tools/main/CoREMOF/data/SI/NCR.zip'
                    }


def _ensure_data_file(file_name):
    """Return a local data file, downloading it atomically only when required."""

    if file_name not in files_to_download:
        raise ValueError(f"Unknown CoRE MOF data file: {file_name}")
    file_path = Path(package_directory, file_name)
    if file_path.is_file():
        return file_path

    file_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with requests.get(files_to_download[file_name], timeout=60, stream=True) as response:
            response.raise_for_status()
            with tempfile.NamedTemporaryFile(dir=file_path.parent, delete=False) as tmp:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        tmp.write(chunk)
                temporary_path = Path(tmp.name)
        temporary_path.replace(file_path)
    except Exception:
        if "temporary_path" in locals():
            temporary_path.unlink(missing_ok=True)
        raise
    return file_path


@lru_cache(maxsize=4)
def _load_json_data(file_name):
    path = _ensure_data_file(file_name)
    with open(path, "r", encoding="utf-8") as stream:
        return json.load(stream)

class download_from_SI():

    """download structures that we got from supporting information.

    Args:
        output_folder (str): path to save structures.

    Returns:
        cif:
            CoRE MOF SI dataset.   
    """

    def __init__(self, output_folder="./CoREMOF2024DB"):
        
        self.SI_path = Path(package_directory, 'data', 'SI')
        self.output = output_folder
        self.run()

    def run(self):

        """start to run. 
        """
            
        cr_zip = _ensure_data_file('data/SI/CR.zip')
        ncr_zip = _ensure_data_file('data/SI/NCR.zip')
        CR_files = self.list_zip(cr_zip)
        NCR_files = self.list_zip(ncr_zip)
     
        os.makedirs(self.output+"/CR/", exist_ok=True)
        os.makedirs(self.output+"/NCR/", exist_ok=True)

        for file in CR_files[:]:
            self.get_from_SI(cr_zip, file, self.output)
        for file in NCR_files[:]:
            self.get_from_SI(ncr_zip, file, self.output)

    def list_zip(self, zip_path):

        """list of files from a ZIP.

        Args:
            zip_path (str): path to ZIP.

        Returns:
            List:
                name list from a ZIP.  
        """
            
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            file_list = zip_ref.namelist()
            return file_list
    
    def get_from_SI(self, zip_path, entry, output_folder):

        """unzip files from a ZIP.

        Args:
            zip_path (str): path to ZIP.
            entry (str): name of structure.
            output_folder (str): path to save structures. 
        """
                
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            file_list = zip_ref.namelist()
            if entry in file_list:
                destination = Path(output_folder).resolve()
                member_path = (destination / entry).resolve()
                if destination not in member_path.parents and member_path != destination:
                    raise ValueError(f"Unsafe path in ZIP archive: {entry}")
                zip_ref.extract(entry, destination)
            
def download_from_CSD(refcode, output_folder="./CoREMOF2024DB"):

    """download structures from CSD, you need to install [CSD python API](https://downloads.ccdc.cam.ac.uk/documentation/API/installation_notes.html) with licence.

    Args:
        refcode (str): CSD refcode.
        output_folder (str): path to save structures.

    Returns:
        cif:
            downloading CIF.  
    """

    try:
        from ccdc import io
    except ImportError as exc:
        raise ImportError(
            "The licensed CSD Python API is required to download CSD structures."
        ) from exc

    csd_reader = io.EntryReader('CSD')
    cryst = csd_reader.crystal(refcode)
    data = cryst.to_string('cif')
    os.makedirs(output_folder, exist_ok=True)
    with open(os.path.join(output_folder, refcode+'.cif'), 'w', encoding='utf-8') as f:
        f.write(data)


def information(dataset, entry, show_units=False):

    """get information of CoRE MOF database.

    Args:
        dataset (str): name of subset.
        entry (str): name of structure.
        show_units (bool): print the dataset unit metadata when ``True``.

    Returns:
        Dictionary:
            properties, DOI, issues and so on. 
    """     

    valid_datasets = {'CR-ASR', 'CR-FSR', 'CR-Ion', 'NCR'}
    if dataset not in valid_datasets:
        raise ValueError(
            f"Unknown dataset {dataset!r}; expected one of {sorted(valid_datasets)}"
        )

    if dataset == "CR-ASR":
        CR_data = _load_json_data('data/CR.json')
        query_data = CR_data["ASR"]
    elif dataset == "CR-FSR":
        CR_data = _load_json_data('data/CR.json')
        query_data = CR_data["FSR"]
    elif dataset == "CR-Ion":
        CR_data = _load_json_data('data/CR.json')
        query_data = CR_data["Ion"]
    else:
        query_data = _load_json_data('data/NCR.json')

    if entry not in query_data:
        suggestions = get_close_matches(entry, query_data.keys(), n=3)
        hint = f" Close matches: {', '.join(suggestions)}" if suggestions else ""
        raise KeyError(f"Entry {entry!r} was not found in {dataset}.{hint}")
    if show_units and dataset.startswith("CR-"):
        print("unit:\n", CR_data["unit"])
    return query_data[entry]

def read_aif(GEMC_data):

    """get adsorption amount of water from GEMC.

    Args:
        GEMC_data (list): from detail_of_CR.json, for example, information("CR-ASR", "2020[Cu][sql]2[ASR]1")["GEMC"].

    Returns:
        Dictionary:
            -   information,by ["info"] always "('_units_loading', 'Molecules/Supercell')".
            -   pressure by ["pressure"].
            -   uptake by ["uptake"].
    """    
        
    if isinstance(GEMC_data, str):
        aif_text = GEMC_data
    else:
        aif_text = "".join(GEMC_data)
    data = cif.read_string(aif_text)
    block = data.sole_block()

    item = block.find_pair_item('_units_loading')

    adsorption_data = {}
    adsorption_data["info"]=item.pair
    adsorption_data["pressure"]=list(block.find_loop('_adsorp_pressure'))
    adsorption_data["uptake"]=list(block.find_loop('_adsorp_amount'))

    return adsorption_data
