# -*- coding: utf-8 -*-
"""Compare hashes with the ones computed for databased of nanoporous materials."""
# ToDo: potentially replace pickle with something else
# that depends less on the Python version
import os
from typing import List

from backports.cached_property import cached_property

from mofchecker.utils import read_pickle

THIS_DIR = os.path.dirname(os.path.realpath(__file__))

__all__ = ["MOFCheckerDB"]


class MOFCheckerDB:
    """Wrapper class for lookup in all databases.

    Databases will be loaded only on first call.
    The output of all lookups is a list of matching entries, and also contains
    information about the database and the mofchecker version that was used
    to compute the hash
    """

    @cached_property
    def graph_hash_dict(self) -> dict:
        """Load the dict of graph hashes."""
        return read_pickle(os.path.join(THIS_DIR, "graph_hash_dict.pkl"))

    @cached_property
    def scaffold_hash_dict(self) -> dict:
        """Load the dict of scaffold hashes."""
        return read_pickle(os.path.join(THIS_DIR, "scaffold_hash_dict.pkl"))

    @cached_property
    def symmetry_hash_dict(self) -> dict:
        """Load the dict of symmetry hashes."""
        return read_pickle(os.path.join(THIS_DIR, "scaffold_hash_dict.pkl"))

    @cached_property
    def composition_dict(self) -> dict:
        """Load the dict of compositions."""
        return read_pickle(os.path.join(THIS_DIR, "composition_dict.pkl"))

    def lookup_graph_hash(self, hash_string: str) -> List[str]:
        """Look up a graph hash_string in the database."""
        return self.graph_hash_dict.get(hash_string, [])

    def lookup_scaffold_hash(self, hash_string: str) -> List[str]:
        """Look up a scaffold hash_string in the database."""
        return self.scaffold_hash_dict.get(hash_string, [])

    def lookup_symmetry_hash(self, hash_string: str) -> List[str]:
        """Loop up symmetry hash_string in the database."""
        return self.symmetry_hash_dict.get(hash_string, [])

    def lookup_composition(self, composition: str) -> List[str]:
        """Look up composition in the database."""
        return self.composition_dict.get(composition, [])
