import os
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
        os.makedirs(self.out_folder, exist_ok=True)
        self.BinaryInteractions()

    def BinaryInteractions(self):
        paras = pd.read_csv(self.FF_csv)
        ffdef = ""
        ffdef += f"# rules to overwrite\n"
        ffdef += f"{len(paras)}\n"
        ffdef += "# pair truncated/shifted tailcorrections\n"
        for _, row in paras.iterrows():
            ffdef += (
                f"{row['name1']}\t{row['name2']}\t"
                + ("shifted\t" if row['shifted'] else "truncated\t")
                + ("yes\n" if row['tailcorrections'] else "no\n")
            )
        valid_rows = paras[~pd.isna(paras["type"]) & ~pd.isna(paras["paras"]) & (paras["type"] != "") & (paras["paras"] != "")]
        ffdef += f"# number of defined interactions\n{len(valid_rows)}\n"
        for _, row in paras.iterrows():
            if pd.isna(row["type"]) or pd.isna(row["paras"]) or row["type"] == "" or row["paras"] == "":
                continue
            ffdef += f"{row['name1']}\t{row['name2']}\t{row['type']}\t{row['paras']}\n"
        ffdef += "# mixing rules to overwrite\n0\n"
        with open(f"{self.out_folder}/force_field.def", "w") as file:
            file.write(ffdef)
