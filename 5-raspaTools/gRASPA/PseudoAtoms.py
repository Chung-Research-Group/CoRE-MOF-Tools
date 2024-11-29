import pandas as pd


DEFAULT_SETTINGS = {
                    "FF_csv": "",
                    "out_folder": ""
                    }


class Generate:
    def __init__(self, **kwargs):
        config = DEFAULT_SETTINGS.copy()
        config.update(kwargs)
        for key, value in config.items():
            setattr(self, key, value)
        self.PseudoAtoms()

    def PseudoAtoms(self):
        paras = pd.read_csv(self.FF_csv)
        padef = ""
        padef += f"#number of pseudo atoms\n"
        padef += str(len(paras))
        padef += "\n#type\tprint\tas\tchem\toxidation\tmass\tcharge\tpolarization\tB-factor\tradii\tconnectivity\tanisotropic\tanisotropic-type\ttinker-type\n"

        for _, row in paras.iterrows():
            padef += (
                f"{row['name']}\t" +
                f"{row['print_to_output']}\t" +
                f"{row['element']}\t" +
                f"{row['print_as']}\t" +
                f"{row['oxidation']}\t" +
                f"{float(row['mass']):.4f}\t" +
                f"{float(row['charge']):.4f}\t" +
                f"{row['polarization']}\t1\t" +
                f"{float(row['parameter2']):.4f}\t0\t0\trelative\t0\n"
            )

        with open(f"{self.out_folder}/pseudo_atoms.def", "w") as file:
            file.write(padef)
