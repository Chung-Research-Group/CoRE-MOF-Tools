"""Process your CIF to get mofid v1 and v2.
"""

from ase.io import read as ase_read
from ase.io import write as ase_write
from ase import neighborlist
import networkx as nx
from ase.build import sort
import ase
import os
from pathlib import Path
import tempfile
import numpy as np
import collections

from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice

from mofid.id_constructor import extract_fragments
from collections import Counter
from mofid.run_mofid import cif2mofid

from pymatgen.analysis.structure_matcher import ElementComparator, StructureMatcher
import selfies as sf

from openbabel import openbabel as ob


metals  = ['Li', 'Na', 'K', 'Rb', 'Cs', 'Fr', 'Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra',
          'Al', 'Ga', 'Ge', 'In', 'Sn', 'Sb', 'Tl', 'Pb', 'Bi', 'Po',
          'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
          'Y' , 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
          'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
          'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'U', 'Tm', 'Yb', 'Lu',
          'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr'
         ]

def run_v1(structure, *, output_path=None):

    """Converting CIF to mofid-v1, see https://snurr-group.github.io/mofid/ for additional installation. Tip: Please check CMAKE, JAVA, etc. before installing.

    Args:
        structure (str): path to your CIF.
        output_path (str, optional): new directory for supporting outputs.
            If omitted, temporary outputs are isolated and removed after use.

    Returns:
        dict: MOFid-v1 identifier and its fragment/topology information.
    """

    structure_path = Path(structure).resolve(strict=True)
    if not structure_path.is_file():
        raise ValueError("MOFid input must be a CIF file")
    if output_path is not None:
        destination = Path(output_path).absolute()
        destination.mkdir(parents=True, exist_ok=False)
        return cif2mofid(str(structure_path), output_path=str(destination))
    with tempfile.TemporaryDirectory(prefix="coremof-mofid-v1-") as directory:
        return cif2mofid(str(structure_path), output_path=directory)

def dict2str(dct):
    """Convert symbol-to-number dict to str.
    """
    return ''.join(symb + (str(n)) for symb, n in dct.items())

def split_nodes_from_cif(structure, prefix='Output'):
    """Split nodes (CIF) to single XYZ.

    Args:
        structure (str): path to your CIF.
        prefix (str): the path to save processed XYZ.

    Returns:
        int:
            -   0: sucess; 1: fail.
    """
    try:
        atoms = ase_read(structure)
    except:
        print('Error with reading CIF in {}'.format(structure))
        return 1
    cutOff = neighborlist.natural_cutoffs(atoms)
    neighborList = neighborlist.NeighborList(cutOff, self_interaction=False, bothways=True, skin=0.3)
    neighborList.update(atoms)
    G = nx.Graph()
    for k in range(len(atoms)):
        tup = (k, {"element":"{}".format(atoms.get_chemical_symbols()[k]), "pos": atoms.get_positions()[k]})
        G.add_nodes_from([tup])
    for k in range(len(atoms)):
        for i in neighborList.get_neighbors(k)[0]:
            G.add_edge(k, i)
    Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
    form_dicts = []
    for index, g in enumerate(Gcc):
        g = list(g)
        fragment = atoms[g]
        fragment = sort(fragment)
        form_dict = fragment.symbols.formula.count()
        form_dicts.append(dict2str(form_dict))
    #nodes = [atoms[list(Gcc[0])]]
    #unique_formdicts = [form_dicts[0]]
    nodes = []
    unique_formdicts = []
    if len(form_dicts) > 1:
        for index, form_dict in enumerate(form_dicts):
            if form_dict not in unique_formdicts:
                nodes.append(atoms[list(Gcc[index])])
                unique_formdicts.append(form_dict)
    elif len(form_dicts) == 1:
        nodes.append(atoms[list(Gcc[0])])
        unique_formdicts.append(form_dicts[0])
    for index, _ in enumerate(nodes):
        ase_write('{}/node{}.xyz'.format(prefix, index), nodes[index]) # /AllNode/nodes.cif
    return 0

def split_linkers_from_cif(structure, prefix='Output'):
    """Split linkers (CIF) to single XYZ.

    Args:
        structure (str): path to your CIF.
        prefix (str): the path to save processed XYZ.

    Returns:
        int:
            -   0: sucess; 1: fail.
    """
    try:
        atoms = ase_read(structure)
    except:
        print('Error with reading CIF in {}'.format(structure))
        return 1
    cutOff = neighborlist.natural_cutoffs(atoms)
    neighborList = neighborlist.NeighborList(cutOff, self_interaction=False, bothways=True, skin=0.3)
    neighborList.update(atoms)
    G = nx.Graph()
    for k in range(len(atoms)):
        tup = (k, {"element":"{}".format(atoms.get_chemical_symbols()[k]), "pos": atoms.get_positions()[k]})
        G.add_nodes_from([tup])
    for k in range(len(atoms)):
        for i in neighborList.get_neighbors(k)[0]:
            G.add_edge(k, i)
    Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
    form_dicts = []
    for index, g in enumerate(Gcc):
        g = list(g)
        fragment = atoms[g]
        fragment = sort(fragment)
        form_dict = fragment.symbols.formula.count()
        form_dicts.append(dict2str(form_dict))
    nodes = []
    unique_formdicts = []
    if len(form_dicts) > 1:
        for index, form_dict in enumerate(form_dicts):
            if form_dict not in unique_formdicts:
                nodes.append(atoms[list(Gcc[index])])
                unique_formdicts.append(form_dict)
    elif len(form_dicts) == 1:
        nodes.append(atoms[list(Gcc[0])])
        unique_formdicts.append(form_dicts[0])
    for index, _ in enumerate(nodes):
        ase_write('{}_MetalOxolinker{}.xyz'.format(prefix, index), nodes[index]) # /MetalOxo/linkers.cif
    return 0

def get_node_linker_files(structure, prefix='Output'):
    """Split MOF to sbu + linker.

    Args:
        structure (str): path to your CIF.
        prefix (str): the path to save processed XYZ.

    Returns:
        files:
            -   structure, node, linker, ...
    """
    os.makedirs(prefix, exist_ok=True)
    extract_fragments(structure, prefix)

def xyz2fomula(xyzpath):
    """from XYZ to chemical fomula based on A-Z.

    Args:
        xyzpath (str): path to your XYZ.

    Returns:
        str:
            -   fomula.
    """
    atoms = ase_read(xyzpath)
    symbols = atoms.get_chemical_symbols()
    counts = Counter(symbols)
    sorted_items = sorted(counts.items(), key=lambda x: x[0])
    # formula = ''.join(f"{el}{n if n > 1 else ''}" for el, n in sorted_items)
    formula = ''.join(f"{el}{n}" for el, n in sorted_items)
    return formula

def convert_ase_pymat(ase_objects):
    """convert to ase to pymatgen atoms.

    Args:
        ase_objects (ase.Atoms): ase-type atoms.

    Returns:
        pymatgen.core.structure.Structure:
            -   pymatgen-type structure.
    """
    structure_lattice = Lattice(ase_objects.cell)
    structure_species = ase_objects.get_chemical_symbols()
    structure_positions = ase_objects.get_positions()
    # ASE positions are Cartesian coordinates in angstrom.  Pymatgen assumes
    # fractional coordinates unless ``coords_are_cartesian`` is explicitly set.
    return Structure(
        structure_lattice,
        structure_species,
        structure_positions,
        coords_are_cartesian=True,
    )

def remove_pbc_cuts(atoms):
    """Remove building block cuts due to periodic boundary conditions. After the
    removal, the atoms object is centered at the center of the unit cell.

    Args:
        atoms (ase.Atoms): aase.Atoms.

    Returns:
        ase.Atoms:
            -   The processed atoms object.
    """

    try:
        # setting cuttoff parameter
        scale = 1.4
        cutoffs = ase.neighborlist.natural_cutoffs(atoms)
        cutoffs = [scale * c for c in cutoffs]
        #making neighbor_list
        #graphs of single metal is not constructed because there is no neighbor.
        I, J, D = ase.neighborlist.neighbor_list("ijD",atoms,cutoff=cutoffs)
        nl = [[] for _ in atoms]
        for i, j, d in zip(I, J, D):
            nl[i].append((j, d))
        visited = [False for _ in atoms]
        q = collections.deque()
        # Center of the unit cell.
        abc_half = np.sum(atoms.get_cell(), axis=0) * 0.5
        positions = {}
        q.append((0, np.array([0.0, 0.0, 0.0])))
        while q:
            i, pos = q.pop()
            visited[i] = True
            positions[i] = pos
            for j, d in nl[i]:
                if not visited[j]:
                    q.append((j, pos + d))
                    visited[j] = True
        centroid = np.array([0.0, 0.0, 0.0])
        for v in positions.values():
            centroid += v
        centroid /= len(positions)
        syms = [None for _ in atoms]
        poss = [None for _ in atoms]
        for i in range(len(atoms)):
            syms[i] = atoms.symbols[i]
            poss[i] = positions[i] - centroid + abc_half
        atoms = ase.Atoms(
            symbols=syms, positions=poss, pbc=True, cell=atoms.get_cell()
        )
        # resize of cell
        cell_x = np.max(atoms.positions[:,0]) - np.min(atoms.positions[:,0])
        cell_y = np.max(atoms.positions[:,1]) - np.min(atoms.positions[:,1])
        cell_z = np.max(atoms.positions[:,2]) - np.min(atoms.positions[:,2])
        cell = max([cell_x,cell_y,cell_z])
        atoms.set_cell([cell+2,cell+2,cell+2, 90,90,90])
        center_mass = atoms.get_center_of_mass()
        cell_half  = atoms.cell.cellpar()[0:3]/2
        atoms.positions = atoms.positions - center_mass + cell_half    
        return atoms
    except:
        return atoms

def run_v2(structure, nodes_dataset, refname):
    """Run the generic MOFid-v2 wrapper with a read-only node library.

    Args:
        structure (str): path to your MOF.
        nodes_dataset (str): path to your node dataset (download from https://github.com/sxm13/CoREMOF_tools/tree/main/CoREMOF/data/mofid/nodes.zip).
        refname (str): file name of your structure or define by yourself.

    Returns:
        str:
            -   mofid-v2.

    Raises:
        ValueError: a node is unmatched or has multiple library matches.
        RuntimeError: fragment extraction or MOFid-v1 evidence is unavailable.

    This generic compatibility method retains ``ltol=0.3`` and ``stol=2``.
    It is not the pinned CoRE-MOF-COD release calculation (``ltol=0.25``, ``stol=1.5``).
    Unavailable nodes never produce invented Type identifiers. Intermediate
    files use one private temporary directory; the caller's ``Output`` and
    the reference node library are never deleted or modified.
    """
    structure_path = Path(structure).resolve(strict=True)
    library = Path(nodes_dataset).resolve(strict=True)
    if not structure_path.is_file() or not library.is_dir():
        raise ValueError("MOFid requires a CIF file and a node-library directory")
    nodes_type = sorted(library.glob("*.xyz"))
    if not nodes_type:
        raise ValueError("The node library contains no XYZ files")
    if not isinstance(refname, str) or not refname or "\n" in refname or "\r" in refname:
        raise ValueError("refname must be a nonempty single-line string")
    with tempfile.TemporaryDirectory(prefix="coremof-mofid-v2-") as directory:
        output = Path(directory) / "fragments"
        mofidv1 = run_v1(structure_path, output_path=output)
        if not isinstance(mofidv1, dict):
            raise RuntimeError("MOFid-v1 returned no fragment evidence")
        linkers = mofidv1.get("smiles_linkers")
        topology, cat = mofidv1.get("topology"), mofidv1.get("cat")
        if (not isinstance(linkers, list) or any(not isinstance(item, str) for item in linkers)
                or not isinstance(topology, str) or not topology or cat is None):
            raise RuntimeError("MOFid-v1 returned incomplete fragment/topology evidence")
        all_linkers = [sf.encoder(linker) for linker in linkers]
        if split_nodes_from_cif(str(output / "AllNode/nodes.cif"), str(output)) != 0:
            raise RuntimeError("MOFid-v2 node extraction is unavailable")
        all_nodes_xyz = sorted(output.glob("node*.xyz"))
        if not all_nodes_xyz:
            raise RuntimeError("MOFid-v2 node extraction produced no nodes")
        all_nodes_part = []
        matcher = StructureMatcher(ltol=0.3, stol=2, angle_tol=5,
                                   primitive_cell=False, scale=False,
                                   comparator=ElementComparator())
        for node_xyz in all_nodes_xyz:
            node_formula = xyz2fomula(node_xyz)
            known_nodes = [path for path in nodes_type if path.stem.split("_", 1)[0] == node_formula]
            a = convert_ase_pymat(remove_pbc_cuts(ase_read(node_xyz)))
            matches = []
            for known_node in known_nodes:
                b = convert_ase_pymat(remove_pbc_cuts(ase_read(known_node)))
                if matcher.fit(a, b):
                    matches.append(known_node.stem)
            if not matches:
                raise ValueError(f"NOT_AVAILABLE_UNMATCHED_NODE: {node_formula}")
            if len(matches) != 1:
                raise ValueError(f"NOT_AVAILABLE_AMBIGUOUS_NODE: {node_formula}: {matches}")
            all_nodes_part.append(matches[0])

        linkers_part = ".".join(all_linkers)
        nodes_part = ".".join(f"[{node}]" for node in all_nodes_part)
        return f"{nodes_part}.{linkers_part} MOFid-v2.{topology}.cat{cat};{refname}"

def are_identical_smiles(smiles1, smiles2):
    obConversion = ob.OBConversion()
    obConversion.SetInFormat("smi")
    obConversion.SetOutFormat("can")  # Set output to canonical SMILES
    mol1 = ob.OBMol()
    mol2 = ob.OBMol()
    # Read the SMILES strings into molecules
    obConversion.ReadString(mol1, smiles1)
    obConversion.ReadString(mol2, smiles2)
    # Convert molecules to canonical SMILES
    can_smiles1 = obConversion.WriteString(mol1).strip()
    can_smiles2 = obConversion.WriteString(mol2).strip()
    # Compare canonical SMILES
    return can_smiles1 == can_smiles2
