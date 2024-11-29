from utils import unit_cell, predict_charge


DEFAULT_SETTINGS = {
                    "FF": "Local",
                    "cutoff": 14,
                    "structure": "",
                    "has_charge": True,
                    "charge_from_cif": True,
                    "EwaldPrecision": 1e-6,
                }


class Generate:
    def __init__(self, **kwargs):
        config = DEFAULT_SETTINGS.copy()
        config.update(kwargs)
        for key, value in config.items():
            setattr(self, key, value)
        if self.has_charge:
            predict_charge(self.structure)
            self.part_charge = self.Charge()
            self.structure += "_pacman"
            self.UC = unit_cell(self.structure, self.cutoff)
            self.UC = f"0\t{self.UC}"
            self.structure = self.structure.split("/")[-1]
        else:
            self.part_charge = "ChargeMethod\nNone\n\n"
            self.UC = unit_cell(self.structure, self.cutoff)
            self.structure = self.structure.split("/")[-1]
        
    def ForceField(self):
        return f"Forcefield\t{self.FF}\nCutOffVDW\t{str(self.cutoff)}\n\n"
    
    def Structure(self):
        return f"FrameworkName\t{self.structure}\nUnitCells\t{self.UC}\n\n"
    
    def Charge(self):
        if self.charge_from_cif:
            return (
                f"ChargeMethod\tEwald\n"
                f"UseChargesFromCIFFile\tyes\n"
                f"EwaldPrecision\t{self.EwaldPrecision}\n\n"
            )
        else:
            return (
                f"ChargeMethod\tEwald\n"
                f"UseChargesFromCIFFile\tno\n"
                f"EwaldPrecision\t{self.EwaldPrecision}\n\n"
            )
    
    def generate_config(self):
        return self.ForceField() + self.Structure() + self.part_charge
