"""Parameters of atom used for open metal site analysis."""

class Atom:

    """A class to hold atomic information, and check bonds. Covalent radii taken from DOI: Covalent radii revisited Beatriz Cordero, Verónica Gómez, Ana E. Platero-Prats, Marc Revés, Jorge Echeverría,  Eduard Cremades, Flavia Barragána and Santiago Alvarez Dalton Trans., 2008, 2832-2838. DOI: 10.1039/B801115J.
    
    Args:
        element: Element of the Atom object. This determines it's properties.
    """

    def __init__(self, element):
  
        self._co_all = {'H': 0.31,
                        'D': 0.31,
                        'He': 0.28,
                        'Li': 1.28,
                        'Be': 0.96,
                        'B': 0.84,
                        'C': 0.73,
                        'N': 0.71,
                        'O': 0.66,  # 0.8,
                        'F': 0.57,
                        'Ne': 0.58,
                        'Na': 1.66,
                        'Mg': 1.41,
                        'Al': 1.21,
                        'Si': 1.11,
                        'P': 1.07,
                        'S': 1.05,
                        'Cl': 1.02,
                        'Ar': 1.06,
                        'K': 2.03,
                        'Ca': 1.76,
                        'Sc': 1.7,
                        'Ti': 1.6,
                        'V': 1.53,
                        'Cr': 1.39,
                        'Mn': 1.5,
                        'Fe': 1.42,
                        'Co': 1.38,
                        'Ni': 1.24,
                        'Cu': 1.32,
                        'Zn': 1.22,
                        'Ga': 1.22,
                        'Ge': 1.2,
                        'As': 1.19,
                        'Se': 1.2,
                        'Br': 1.2,
                        'Kr': 1.16,
                        'Rb': 2.2,
                        'Sr': 1.95,
                        'Y': 1.9,
                        'Zr': 1.75,
                        'Nb': 1.64,
                        'Mo': 1.54,
                        'Tc': 1.47,
                        'Ru': 1.46,
                        'Rh': 1.42,
                        'Pd': 1.39,
                        'Ag': 1.45,
                        'Cd': 1.44,
                        'In': 1.42,
                        'Sn': 1.39,
                        'Sb': 1.39,
                        'Te': 1.38,
                        'I': 1.39,
                        'Xe': 1.4,
                        'Cs': 2.44,
                        'Ba': 2.15, #1.80
                        'La': 2.07,
                        'Ce': 2.04,
                        'Pr': 2.03,
                        'Nd': 2.01,
                        'Pm': 1.99,
                        'Sm': 1.98,
                        'Eu': 1.98,
                        'Gd': 1.96,
                        'Tb': 1.94,
                        'Dy': 1.92,
                        'Ho': 1.92,
                        'Er': 1.89,
                        'Tm': 1.9,
                        'Yb': 1.87,
                        'Lu': 1.87,
                        'Hf': 1.75,
                        'Ta': 1.7,
                        'W': 1.62,
                        'Re': 1.51,
                        'Os': 1.44,
                        'Ir': 1.41,
                        'Pt': 1.36,
                        'Au': 1.36,
                        'Hg': 1.32,
                        'Tl': 1.45,
                        'Pb': 1.46,
                        'Bi': 1.48,
                        'Po': 1.4,
                        'At': 1.5,
                        'Rn': 1.5,
                        'Fr': 2.6,
                        'Ra': 2.21,
                        'Ac': 2.15,
                        'Th': 2.06,
                        'Pa': 2,
                        'U': 1.96,
                        'Np': 1.9,
                        'Pu': 1.87,
                        'Am': 1.8,
                        'Cm': 1.69}

        elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
                    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
                    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
                    'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
                    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
                    'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
                    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
                    'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm']

        self._list_of_non_metals = ['H', 'D', 'B', 'C', 'N', 'O', 'F',
                                    'P', 'S', 'Cl', 'Se', 'Br', 'I']

        self.element = element
        if self.element == 'D':
            self.atomic_number = 1
        else:
            self.atomic_number = elements.index(self.element) + 1
        self._co = self._co_all[self.element]

    @property
    def co(self):
        """Coordination radius."""
        return self._co

    @property
    def is_metal(self):
        """Check if atom is metal or not."""
        return self.element not in self._list_of_non_metals

    @property
    def is_lanthanide_or_actinide(self):
        """Check if atom is a lanthanide or actinide."""
        return self.is_lanthanide or self.is_actinide

    @property
    def is_lanthanide(self):
        """Check if atom is a lanthanide."""
        return 72 > self.atomic_number > 56

    @property
    def is_actinide(self):
        """Check if atom is a actinide."""
        return 97 > self.atomic_number > 88

    def bond_tolerance(self, ele2):
        """Determine if atom is a actinide."""
        if self._check_if_heavy_metal_bond(ele2):
            return 0.2
        else:
            return 0.5  # 0.4

    def _check_if_heavy_metal_bond(self, ele2):
        """Determine if atom is a actinide."""
        return self.is_heavy_metal or Atom(ele2).is_heavy_metal

    @property
    def is_heavy_metal(self):  # \m/
        """Determine if atom has a covelant radii larger than 1.95."""
        return self.co > 1.95

    def check_bond(self, ele2, dist, bond_tol=None):
        """Check if the atom is bonded with a given atom"""
        max_bond = self.max_bond(ele2, bond_tol)
        return dist < max_bond

    def max_bond(self, ele2, bond_tol=None):

        """Get the maximum possible distance between the Atom object and another atom of type ele2.

        Args:
            ele2: Element of the atom that the max_bond corresponds to.
            bond_tol: Bold tolerance to use, if not set then the default will be used.

        Returns:    
            The some of covelant radii for Atom and Atom(ele2) plus their covalant radii.
        """

        if bond_tol is None:
            bond_tol = self.bond_tolerance(ele2)
        return Atom(ele2).co + self.co + bond_tol
