"""Types reused across the package."""
from pathlib import Path
from typing import Union

from pymatgen.core import IStructure, Structure
from typing_extensions import TypeAlias

PathType: TypeAlias = Union[str, Path]
StructureIStructureType: TypeAlias = Union[Structure, IStructure]
