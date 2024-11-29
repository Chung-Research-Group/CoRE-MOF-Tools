DEFAULT_SETTINGS = {
                    "task": "SingleAdsorption",
                    "adsorbate": [],
                    "MolecFile": "Local",
                    "Fraction": [],
                    "WIG": 1,
                    "Fugacity": 1,
                    "TranslationProbability": 1.0,
                    "RotationProbability": 1.0,
                    "ReinsertionProbability": 1.0,
                    "SwapProbability": 1.0,
                    "IdentityChangeProbability": 1.0,
                    "WidomProbability": 1.0,
                    "CreateNumberOfMolecules": 0,
                }


class Generate:
    def __init__(self, **kwargs):
        config = DEFAULT_SETTINGS.copy()
        config.update(kwargs)
        for key, value in config.items():
            setattr(self, key, value)

    def generate_config(self):
        if self.task == "SingleAdsorption":
            return self.SingleComponent()
        elif self.task == "MixtureAdsorption":
            return self.MixtureComponent()
        else:
            raise ValueError(f"Unsupported task type: {self.task}")

    def SingleComponent(self):
        content = (
            f"Component\t0\n"
            f"MoleculeName\t{self.adsorbate[0]}\n"
            f"MoleculeDefinition\t{self.MolecFile}\n"
            f"TranslationProbability\t{self.TranslationProbability}\n"
            f"RotationProbability\t{self.RotationProbability}\n"
            f"ReinsertionProbability\t{self.ReinsertionProbability}\n"
            f"SwapProbability\t{self.SwapProbability}\n"
            f"CreateNumberOfMolecules\t{self.CreateNumberOfMolecules}\n"
        )
        if self.WIG != 1:
            content += f"IdealRosenbluthValue\t{self.WIG}\n"
        if self.Fugacity != 1:
            content += f"FugacityCoefficient\t{self.Fugacity}\n"
        return content

    def MixtureComponent(self):
        content = ""
        for i in range(len(self.adsorbate)):
            adsorbate_list = " ".join(map(str, range(len(self.adsorbate))))
            content += (
                f"Component\t{i}\n"
                f"MoleculeName\t{self.adsorbate[i]}\n"
                f"MoleculeDefinition\t{self.MolecFile}\n"
                f"MolFraction\t{self.Fraction[i]}\n"
                f"TranslationProbability\t{self.TranslationProbability}\n"
                f"RotationProbability\t{self.RotationProbability}\n"
                f"ReinsertionProbability\t{self.ReinsertionProbability}\n"
                f"SwapProbability\t{self.SwapProbability}\n"
                f"IdentityChangeProbability\t{self.IdentityChangeProbability}\n"
                f"NumberOfIdentityChanges\t{len(self.adsorbate)}\n"
                f"IdentityChangesList\t{adsorbate_list}\n"
                f"CreateNumberOfMolecules\t{self.CreateNumberOfMolecules}\n\n"
            )
        if self.WIG != 1:
            content += f"IdealRosenbluthValue\t{self.WIG}\n"
        if self.Fugacity != 1:
            content += f"FugacityCoefficient\t{self.Fugacity}\n"
        return content