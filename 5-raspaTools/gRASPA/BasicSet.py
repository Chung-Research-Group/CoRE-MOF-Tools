DEFAULT_SETTINGS = {
                    "cycles": 20000,
                    "intiacycles": 10000
                    }


class Generate:
    def __init__(self, **kwargs):
        config = DEFAULT_SETTINGS.copy()
        config.update(kwargs)
        for key, value in config.items():
            setattr(self, key, value)

    def BasicSetting(self):

        return (
            f"SimulationType\tMonteCarlo\n"
            f"NumberOfProductionCycles\t{self.cycles}\n"
            f"NumberOfEquilibrationCycles\t0\n"
            f"NumberOfInitializationCycles\t{self.intiacycles}\n"
            f"RestartFile\tno\n"
            f"RandomSeed\t0\n"
            f"NumberOfTrialPositions\t10\n"
            f"NumberOfTrialOrientations\t10\n"
            f"NumberOfBlocks\t1\n"
            f"AdsorbateAllocateSpace\t10240\n"
            f"NumberOfSimulations\t1\n"
            f"SingleSimulation\tyes\n"
            f"InputFileType\tcif\n"
            f"UseDNNforHostGuest\tno\n"
            f"UseMaxStep\t0\n\n"
        )

        
    def generate_config(self):
        return self.BasicSetting()