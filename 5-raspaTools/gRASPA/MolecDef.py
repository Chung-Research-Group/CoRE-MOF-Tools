import pandas as pd
import os

DEFAULT_SETTINGS = {
    "Molec_csv": [],
    "out_folder": "",
}


class Generate:
    def __init__(self, **kwargs):
        config = DEFAULT_SETTINGS.copy()
        config.update(kwargs)
        for key, value in config.items():
            setattr(self, key, value)
        self.MoleculeDefinitions()

    def MoleculeDefinitions(self):
        for csv_file in self.Molec_csv:
            try:
                MoleData = pd.read_csv(csv_file)
            except FileNotFoundError:
                print(f"File not found: {csv_file}")
                continue
            FileName = os.path.basename(csv_file).replace(".csv", "")
            moldef = ""
            moldef += "# critical constants: Temperature [T], Pressure [Pa], and Acentric factor [-]\n"
            moldef += f"{MoleData['CriticalTemperature'].iloc[0]}\n"
            moldef += f"{MoleData['CriticalPressure'].iloc[0]}\n"
            moldef += f"{MoleData['AcentricFactor'].iloc[0]}\n"
            moldef += f"# Number Of Atoms\n{len(MoleData)}\n# Number of groups\n1\n# group_type\nrigid\n# atomic positions\n"
            for i, row in MoleData.iterrows():
                if pd.isna(row["atom"]) or pd.isna(row["posx"]) or pd.isna(row["posy"]) or pd.isna(row["posz"]):
                    continue
                moldef += (
                    f"{i}\t{row['atom']}\t{int(row['posx'])}\t{int(row['posy'])}\t{float(row['posz']):.2f}\n"
                )
            moldef += "# Chiral\tcenters\tBond\tBondDipoles\tBend\tUrayBradley\tInvBend\tTorsion\tImp.\tTorsion\tBond/Bond\tStretch/Bend\tBend/Bend\tStretch/Torsion\tBend/Torsion\tIntraVDW\tIntraCoulomb\n"
            moldef += f"0\t{len(MoleData['Bonds0'].dropna().values)}\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"
            moldef += "#\tBond\tstretch:\tatom\tn1-n2,\ttype,\tparameters\n"
            for i, row in MoleData.iterrows():
                if pd.isna(row["Bonds0"]) or pd.isna(row["Bonds1"]):
                    continue
                moldef += (
                    f"{i}\t{int(row['Bonds0'])}\t{int(row['Bonds1'])}\tRIGID_BOND\n"
                )
            moldef += "# Number of config moves\n0\n"
            with open(f"{self.out_folder}/{FileName}.def", "w") as file:
                file.write(moldef)

