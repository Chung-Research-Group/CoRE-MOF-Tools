import os, csv, glob, warnings
from ase.io import read, write
from pymatgen.core import Structure
from pymatgen.io.cif import CifWriter
from .utils.atoms_definitions import METAL
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
warnings.filterwarnings("ignore")

class run:
    def __init__(self, input_folder, output_folder):
        self.input_folder = input_folder
        self.input_files = glob.glob(os.path.join(input_folder, "*.cif"))
        self.output = output_folder
        self.csv_path = os.path.join(output_folder, "preprocess_log.csv")
        self.process()

    def process(self):
        with open(self.csv_path, mode='w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(["File", "Multi Structure Count", "Has Metal", "Has Carbon", "P1 Status"])
            
            for mof in self.input_files:
                name = os.path.basename(mof).replace(".cif", "")
                try:
                    structures = read(mof, index=':')
                    multi_count = len(structures)
                    print(name, "has", multi_count, "multi structures")

                    for i, structure in enumerate(structures):
                        struct_name = f"{name}_{i+1}" if multi_count > 1 else name
                        temp_path = os.path.join(self.output, f"{struct_name}.cif")
                        write(temp_path, structure, format='cif')

                        has_metal, has_carbon = self._check_structure(structure, struct_name)
                        
                        primitive_path = self.make_primitive(temp_path)
                        if primitive_path:
                            p1_status = self.make_p1(primitive_path)
                        else:
                            p1_status = self.make_p1(temp_path)
                            p1_status = "Primitive fail"

                        csv_writer.writerow([struct_name, multi_count, has_metal, has_carbon, p1_status])

                except Exception as e:
                    print(name, "fail read:", e)
                    csv_writer.writerow([name, 0, "N/A", "N/A", "fail read"])

    def _check_structure(self, structure, name):
        has_metal = any(METAL.get(atom.symbol) for atom in structure)
        has_carbon = any(atom.symbol == 'C' for atom in structure)
        if not has_metal:
            print(name, "missing Metal")
        if not has_carbon:
            print(name, "missing Carbon")
        return has_metal, has_carbon

    def make_primitive(self, mof_path):
        try:
            structure = Structure.from_file(mof_path, primitive=True)
            primitive_path = os.path.join(self.output, f"{os.path.basename(mof_path)}")
            structure.to(filename=primitive_path)
            return primitive_path
        except Exception as e:
            print(f"Primitive conversion failed for {mof_path}: {e}")
            return None

    def make_p1(self, primitive_path):
        try:
            structure = Structure.from_file(primitive_path)
            sga = SpacegroupAnalyzer(structure)
            structure_p1 = sga.get_primitive_standard_structure(international_monoclinic=True, keep_site_properties=False)
            p1_path = os.path.join(self.output, f"{os.path.basename(primitive_path)}")
            writer = CifWriter(structure_p1)
            writer.write_file(p1_path)
            return "P1 complete"
        except Exception as e:
            print(f"P1 conversion failed for {primitive_path}: {e}")
            return "P1 fail"
