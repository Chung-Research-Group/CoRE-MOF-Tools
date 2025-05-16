"""Download structures and query information of CoRE MOF Database.
"""

try:
    from ccdc import io
    csd_reader = io.EntryReader('CSD')
except:
    print("You need to install CSD software with the license if you want to download all structures")

import os, json, requests, zipfile

from gemmi import cif

package_directory = os.path.abspath(__file__).replace("structure.py","")

files_to_download = {
                    'data/CR.json': 'https://raw.githubusercontent.com/sxm13/CoREMOF_tools/main/CoREMOF/data/info/CR.json',
                    'data/NCR.json': 'https://raw.githubusercontent.com/sxm13/CoREMOF_tools/main/CoREMOF/data/info/NCR.json',
                    'data/SI/CR.zip': 'https://raw.githubusercontent.com/sxm13/CoREMOF_tools/main/data/SI/CR.zip',
                    'data/SI/NCR.zip': 'https://raw.githubusercontent.com/sxm13/CoREMOF_tools/main/data/SI/NCR.zip'
                    }

for file_name, url in files_to_download.items():
    
    file_path = os.path.join(package_directory, file_name)
    directory = os.path.dirname(file_path) 

    if not os.path.exists(directory):
        os.makedirs(directory)
        
    if not os.path.exists(file_path):
        response = requests.get(url)
        if response.status_code == 200:
            with open(file_path, 'wb') as f:
                f.write(response.content)
            print(f"Downloaded {file_name} to {file_path}")
        else:
            print(f"Failed to download {file_name} from {url}")
    else:
        pass

class download_from_SI():

    """download structures that we got from supporting information.

    Args:
        output_folder (str): path to save structures.

    Returns:
        cif:
            CoRE MOF SI dataset.   
    """

    def __init__(self, output_folder="./CoREMOF2024DB"):
        
        self.SI_path = package_directory+'/data/SI/'
        self.output = output_folder
        self.run()

    def run(self):

        """start to run. 
        """
            
        CR_files = self.list_zip(self.SI_path+"CR.zip")
        NCR_files = self.list_zip(self.SI_path+"NCR.zip")
     
        os.makedirs(self.output+"/CR/", exist_ok=True)
        os.makedirs(self.output+"/NCR/", exist_ok=True)

        for file in CR_files[:]:
            self.get_from_SI(self.SI_path+"CR.zip", file, self.output)
        for file in NCR_files[:]:
            self.get_from_SI(self.SI_path+"NCR.zip", file, self.output)

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
                zip_ref.extract(entry, output_folder)
            
def download_from_CSD(refcode, output_folder="./CoREMOF2024DB"):

    """download structures from CSD, you need to install [CSD python API](https://downloads.ccdc.cam.ac.uk/documentation/API/installation_notes.html) with licence.

    Args:
        refcode (str): CSD refcode.
        output_folder (str): path to save structures.

    Returns:
        cif:
            downloading CIF.  
    """

    cryst = csd_reader.crystal(refcode)
    data = cryst.to_string('cif')
    f = open(os.path.join(output_folder, refcode+'.cif'),'w')
    f.write(data)
    f.close()


def information(dataset, entry):

    """get information of CoRE MOF database.

    Args:
        dataset (str): name of subset.
        entry (str): name of structure

    Returns:
        Dictionary:
            properties, DOI, issues and so on. 
    """     

    CR_data_path = package_directory+'/data/CR.json'
    NCR_data_path = package_directory+'/data/NCR.json'

    with open (CR_data_path, "r") as CR_f:
        CR_data = json.load(CR_f)
        
    with open (NCR_data_path, "r") as NCR_f:
        NR_data = json.load(NCR_f)
    
    if dataset == "CR-ASR":
        print("unit:\n", CR_data["unit"])
        query_data = CR_data["ASR"]
    elif dataset == "CR-FSR":
        print("unit:\n", CR_data["unit"])
        query_data = CR_data["FSR"]
    elif dataset == "CR-Ion":
        print("unit:\n", CR_data["unit"])
        query_data = CR_data["Ion"]
    elif dataset == "NCR":
        query_data = NR_data

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
        
    with open("temp_gemc.aif", "w") as f:
        f.write("".join(GEMC_data) )
    data = cif.read_file("temp_gemc.aif")
    os.remove("temp_gemc.aif")
    block = data.sole_block()

    item = block.find_pair_item('_units_loading')

    adsorption_data = {}
    adsorption_data["info"]=item.pair
    adsorption_data["pressure"]=list(block.find_loop('_adsorp_pressure'))
    adsorption_data["uptake"]=list(block.find_loop('_adsorp_amount'))

    return adsorption_data