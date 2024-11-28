from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor
import os, csv, glob
from ase.io import read


class run:
    def __init__(self, input_folder, output_folder):
        self.input_folder = input_folder
        self.input_files = glob.glob(os.path.join(input_folder, "*.cif"))
        self.output = output_folder
        self.csv_path = os.path.join(output_folder, "others_log.csv")
        self.process()

    def process(self):
        with open(self.csv_path, mode='w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(["File", "hall_symbol", "space_group_number", "crystal_system", "total_mass", "number_atoms"])

            for mof in self.input_files:
                name = os.path.basename(mof).replace(".cif", "")
                try:
                    atoms = read(mof)
                    structure = AseAtomsAdaptor.get_structure(atoms)
                    result = SpacegroupAnalyzer(structure, symprec=0.01, angle_tolerance=5)
                    total_mass = atoms.get_masses().sum()
                    number_atoms = len(atoms)
                    
                    hall_symbol = result.get_hall()
                    space_group_number = result.get_space_group_number()
                    crystal_system = result.get_crystal_system()

                    csv_writer.writerow([name, hall_symbol, space_group_number, crystal_system, total_mass, number_atoms])

                except Exception as e:
                    csv_writer.writerow([name, "Error", "Error", "Error", "Error", "Error"])
                    print(f"Error processing {name}: {e}")

