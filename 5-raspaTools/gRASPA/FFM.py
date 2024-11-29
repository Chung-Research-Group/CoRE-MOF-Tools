import pandas as pd


DEFAULT_SETTINGS = {
                    "FF_csv": "",
                    "tail_correction": True,
                    "shifted": False,
                    "mixing_rule": "Lorentz-Berthelot",
                    "out_folder": ""
                    }


class Generate:
    def __init__(self, **kwargs):
        config = DEFAULT_SETTINGS.copy()
        config.update(kwargs)
        for key, value in config.items():
            setattr(self, key, value)
        self.SelfInteractions()

    def SelfInteractions(self):
        paras = pd.read_csv(self.FF_csv)
        ffmdef = ""
        ffmdef += f"# general rule for shifted vs truncated\n{'shifted' if self.shifted else 'truncated'}\n"
        ffmdef += f"# general rule tailcorrections\n{'yes' if self.tail_correction else 'no'}\n# number of defined interactions\n"
        ffmdef += str(len(paras))
        ffmdef += "\n# type interaction\n"
        for _, row in paras.iterrows():
            ffmdef += (
                        row["name"] + "\t" +
                        row["type"] + "\t" +
                        f"{row['parameter1']:.4f}" + "\t" +
                        f"{row['parameter2']:.4f}" + "\t" +
                        "#" +
                        row["source"] +
                        "\n"
                    )
        ffmdef += "# general mixing rule\n" + self.mixing_rule + "\n"
        with open(f"{self.out_folder}/force_field_mixing_rules.def", "w") as file:
            file.write(ffmdef)