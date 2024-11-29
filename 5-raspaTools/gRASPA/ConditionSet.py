DEFAULT_SETTINGS = {
                    "T": 298,
                    "P": 100000,
                    }


class Generate:
    def __init__(self, **kwargs):
        config = DEFAULT_SETTINGS.copy()
        config.update(kwargs)
        for key, value in config.items():
            setattr(self, key, value)

    def ConditionSetting(self):
        return (
            f"Temperature\t{self.T}\n"
            f"Pressure\t{self.P}\n\n"
        )
    
    def generate_config(self):
        return self.ConditionSetting()