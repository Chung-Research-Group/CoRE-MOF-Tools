import os,re,glob
from ase.io import read,write
from ase import Atoms


def ddec2cif(input_folder, output_folder):
    structures = glob.glob(os.path.join(input_folder))
    t=len(structures)
    for name in structures:
        try:
            name = name.strip()
            atoms = read(input_folder+name+'/DDEC6_even_tempered_net_atomic_charges.xyz')
            with open(input_folder+name+"/DDEC6_even_tempered_net_atomic_charges.xyz",'r') as chg:
                lines=chg.readlines()
            unitcell = lines[1]
            match = re.search(r"unitcell \[\{(.+?)\}\]", unitcell)
            if match:
                cell_params_str = match.group(1).split("}, {")
                cell_params = [list(map(float, re.findall(r"-?\d+\.\d+", vector))) for vector in cell_params_str]
            atoms = Atoms(atoms,cell=cell_params)
            write(output_folder+name.split("/")[-1]+'_ddec.cif', atoms)
            with open(input_folder+name+"/DDEC6_even_tempered_net_atomic_charges.xyz",'r') as chg:
                lines=chg.readlines()
            total_atoms = int(lines[0])
            total_lines = total_atoms+2
            charges = []
            for line in lines[2:total_lines]: 
                parts = line.split()
                try:
                    charge = parts[4:5][0]
                    charges.append(charge)
                except:
                    pass
            with open(output_folder+name.split("/")[-1]+'_ddec.cif', 'r') as f:
                lines = f.readlines()
            lines[0] = "# ddec2cif by CoRE MOF 2024 Tools (https://github.com/sxm13/CoRE-MOF-2024-Tools)\n"
            lines[1] = "data_"+name+"\n"
            for i, line in enumerate(lines):
                if '_atom_site_occupancy' in line:
                    lines.insert(i + 1, "  _atom_site_charge\n")
                    break
            charge_index = 0
            for j in range(i + 2, len(lines)):
                if charge_index < len(charges):
                    lines[j] = lines[j].strip() + " " + str(charges[charge_index]) + "\n"
                    charge_index += 1
                else:
                    break
            with open(output_folder+"/"+name.split("/")[-1]+'_ddec.cif', 'w') as file:
                file.writelines(lines)
            file.close()
            with open(output_folder+"/"+name.split("/")[-1]+'_ddec.cif', 'r') as file:
                content = file.read()
            file.close()
            new_content = content.replace('_space_group_name_H-M_alt', '_symmetry_space_group_name_H-M')
            new_content = new_content.replace('_space_group_IT_number', '_symmetry_Int_Tables_number')
            new_content = new_content.replace('_space_group_symop_operation_xyz', '_symmetry_equiv_pos_as_xyz')
            with open(output_folder+"/"+name.split("/")[-1]+'_ddec.cif', 'w') as file:
                file.write(new_content)
            file.close
        except:
            print(name)
    return t