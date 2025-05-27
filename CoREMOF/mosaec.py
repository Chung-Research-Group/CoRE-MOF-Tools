#!/usr/bin/env python3
from __future__ import annotations

# handle differing python versions (< 3.8, & > 3.9) of CSD python API
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

import os
import math

import pandas as pd

import mendeleev

from ccdc import io
from ccdc import crystal
from ccdc import molecule

from ccdc.crystal import Crystal
from ccdc.molecule import Atom, Bond, Molecule

from ccdc import descriptors


def readentry(input_cif: str) -> Crystal:
    """
    Reads a CIF file containing structure data and converts it to a
    standard atom labeling convention using the ccdc.crystal module.

    Parameters:
        input_cif (str): filename (.CIF) containing crystal structure data.

    Returns:
        newcif (ccdc.crystal.Crystal): Crystal object containing structural data
                                       in the standard atom labeling convention.
    """
    # read in the cif to a crystal object
    with io.CrystalReader(input_cif, format="cif") as readcif:
        cif = readcif[0]
    readcif.close()

    # to remove duplicate atoms, need the empirical formula
    formula = cif.formula
    elamnt = formula.split(" ")

    # now convert to standard labelling convention and identify
    # duplicate atoms to be removed
    with open(input_cif, "r") as file:
        file.seek(0)
        newstring = str()
        lines = file.readlines()
        loop_pos = 0
        start = 0
        end = 0
        columncount = 0
        type_pos = 0
        label_pos = 0
        x_pos = 0
        y_pos = 0
        z_pos = 0
        for i, line in enumerate(lines):
            lines[i] = lines[i].lstrip()
        for i, line in enumerate(lines):
            # locate atom type and site label columns
            if "loop_" in line:
                loop_pos = i
            if ("_atom" in line) and (not "_geom" in line) and (not "_aniso" in line):
                start = loop_pos + 1
                end = i + 1
        for i in range(start, end):
            if "atom_site_type_symbol" in lines[i]:
                type_pos = columncount
            if "atom_site_label" in lines[i]:
                label_pos = columncount
            columncount += 1
        counting = {}
        cutoff = {}
        to_remove = []
        for i in range(end, len(lines)):
            if "loop_" in lines[i]:
                break
            # lines with atom information will contain a ., so only look at these
            if "." in lines[i]:
                # split lines by whitespace
                col = lines[i].split()
                # keep count of how many of each element type
                if not col[type_pos] in counting:
                    counting[col[type_pos]] = 1
                elif col[type_pos] in counting:
                    counting[col[type_pos]] += 1
                # new atom labels
                newlabel = f"{col[type_pos]}{counting[col[type_pos]]}"
                lines[i] = lines[i].replace(col[label_pos], newlabel)
                # cutoff repeated atoms
                if newlabel in elamnt:
                    cutoff[col[type_pos]] = counting[col[type_pos]]
                if col[type_pos] in cutoff:
                    if counting[col[type_pos]] > cutoff[col[type_pos]]:
                        to_remove.append(lines[i])
        # remove unnecessary atoms
        for i in to_remove:
            lines.remove(i)
        # combine to new string
        for i in lines:
            newstring += i
        # read into new crystal object and assign bonds
        newcif = crystal.Crystal.from_string(newstring, format="cif")
        newcif.assign_bonds()
        file.close()
    return newcif


def readSBU(input_mol2: str) -> Crystal:
    """
    Reads a MOL2 file containing SBU/metal complex structural data
    and converts it to a standard atom labeling convention using
    the ccdc.crystal module.

    Parameters:
        input_mol2 (str): filename (.mol2) containing SBU/metal complex
                          structural data.

    Returns:
        mol (ccdc.crystal.Crystal): Crystal object containing structural data
                                    in the standard atom labeling convention.
    """
    # First have to convert the connection points (unknown atom types)
    # to H. Has to be done as a string so CSD reader can
    # interpret bonding
    with open(input_mol2, "r") as file:
        file.seek(0)
        newstring = str()
        lines = file.readlines()
        for i, line in enumerate(lines):
            if "*" in line:
                lines[i] = lines[i].replace("Du", "H")
            newstring += lines[i]
        cif = crystal.Crystal.from_string(newstring, format="mol2")
        cif.assign_bonds()
        file.close()

    mol = cif.molecule
    # MOSAEC needs atoms to have unique labels, so make them unique
    count = 1
    for atom in mol.atoms:
        atom.label = f"{atom.label}{count}"
        count += 1

    return mol


def read_CSD_entry(input_refcode: str) -> Crystal:
    """
    Read entries directly from the CSD CrystalReader according to CSD refcode.

    Parameters:
        input_refcode (str): string used to identify materials in the CSD.

    Returns:
        cif (ccdc.crystal.Crystal): Crystal object containing structural data
                                    in the standard atom labeling convention.
    """
    # read in the cif to a crystal object
    csd_crystal_reader = io.CrystalReader("CSD")
    cif = csd_crystal_reader.crystal(input_refcode)
    cif.assign_bonds()
    csd_crystal_reader.close()
    return cif


def get_no_metal_molecule(inputmolecule: Molecule) -> Molecule:
    """
    Remove metal atoms from the input Molecule object.

    Parameters:
        inputmolecule (ccdc.molecule.Molecule): original Molecule object.

    Returns:
        workingmol (ccdc.molecule.Molecule): Molecule object with all metal
                                             atoms removed.
    """
    workingmol = inputmolecule.copy()
    for atom in workingmol.atoms:
        if atom.is_metal:
            workingmol.remove_atom(atom)
    workingmol.assign_bond_types(which="All")
    return workingmol


def get_unique_sites(mole: Molecule, asymmole: Molecule) -> list[Atom]:
    """
    Get the unique atoms in a structure belonging to the asymmetric unit.

    Parameters:
        mole (ccdc.molecule.Molecule): original structure Molecule object.
        asymmole (ccdc.molecule.Molecule): asymmetric unit of the structure.

    Returns:
        uniquesites (list[ccdc.molecule.Atom]): list of unique atoms in the structure
                                                that belong to the asymmetric unit.
    """
    # blank list for unique sites
    uniquesites = []
    labels = []
    asymmcoords = []
    molecoords = []
    duplicates = []
    for atom in asymmole.atoms:
        asymmcoords.append(atom.coordinates)
    for atom in mole.atoms:
        if atom.coordinates in asymmcoords:
            if not atom.coordinates in molecoords:
                if not atom.label in labels:
                    uniquesites.append(atom)
                    molecoords.append(atom.coordinates)
                    labels.append(atom.label)
                else:
                    duplicates.append(atom)
            else:
                duplicates.append(atom)
    if len(duplicates) >= 1:
        for datom in duplicates:
            for atom in uniquesites:
                if any(
                    [
                        (datom.coordinates == atom.coordinates),
                        (datom.label == atom.label),
                    ]
                ):
                    if datom.atomic_symbol == atom.atomic_symbol:
                        if len(datom.neighbours) > len(atom.neighbours):
                            uniquesites.remove(atom)
                            uniquesites.append(datom)
                    if not datom.label in labels:
                        uniquesites.append(datom)
                        labels.append(datom.label)
    return uniquesites


def get_metal_sites(sites: list[Atom]) -> list[Atom]:
    """
    Get the metal sites in a structure belonging to the asymmetric unit.

    Parameters:
        sites (list[ccdc.molecule.Atom]): list of unique atoms in the structure
                                          that belong to the asymmetric unit.

    Returns:
        metalsites (list[ccdc.molecule.Atom]): list of metal sites in the structure
                                               that belong to the asymmetric unit.
    """
    metalsites = []
    for site in sites:
        if site.is_metal == True:
            metalsites.append(site)
    return metalsites


def get_ligand_sites(
    metalsites: list[Atom], sites: list[Atom]
) -> dict[Atom, list[Atom]]:
    """
    Get the ligand sites binding each metal atom in a structure.

    Parameters:
        metalsites (list[ccdc.molecule.Atom]): list of metal sites in the structure
                                               that belong to the asymmetric unit.
        sites (list[ccdc.molecule.Atom]):  list of unique atoms in the structure
                                           that belong to the asymmetric unit.

    Returns:
        metal_sphere (dict[ccdc.molecule.Atom, list[ccdc.molecule.Atom]]):
                        dictionary with metal Atom object as keys and the the list
                        of ligand atoms which bind them as values.
    """
    metal_sphere = {}
    for metal in metalsites:
        sphere1 = []
        for ligand in metal.neighbours:
            if not ligand.is_metal == True:
                for site in sites:
                    if ligand.label == site.label:
                        sphere1.append(site)
        metal_sphere[metal] = sphere1
    return metal_sphere


def get_binding_sites(metalsites: list[Atom], uniquesites: list[Atom]) -> list[Atom]:
    """
    Get the binding sites in a structure, given the list of unique metal atoms
    and all unique atoms.

    Parameters:
        metalsites (list[ccdc.molecule.Atom]): list of unique metal atoms.
        uniquesites (list[ccdc.molecule.Atom]): list of unique atoms.

    Returns:
        binding_sites (list[ccdc.molecule.Atom]): list of binding sites connecting
                                                  metal atoms and ligands.
    """
    binding_sites = set()
    for metal in metalsites:
        for ligand in metal.neighbours:
            for site in uniquesites:
                if ligand.label == site.label:
                    binding_sites.add(site)
    return binding_sites


def ringVBOs(mole: Molecule) -> dict[int, int]:
    """
    Calculates the VBO (valence bond order) for each atom in the structure.

    Parameters:
        mole (ccdc.molecule.Molecule): Molecule object representing the structure.

    Returns:
        ringVBO (dict[int, int]): dictionary with each atom's index in mole.atoms
                                  as keys and VBO (valence bond order) as values.
    """
    ringVBO = {}
    unassigned = mole.atoms
    ringcopy = mole.copy()
    oncycle_atoms = []
    offcycle_atoms = []
    oncycle_labels = []
    offcycle_labels = []
    cyclic_periodic = []
    cyclic_periodic_labels = []
    offcycle_periodic = []

    # remove all the metals, this
    # prevents metal-containing rings (i.e. pores)
    # from interfering
    for atom in ringcopy.atoms:
        if atom.is_metal:
            ringcopy.remove_atom(atom)

    # collect all the cyclic atoms
    for atom in ringcopy.atoms:
        if atom.is_cyclic:
            if not atom in oncycle_atoms:
                oncycle_atoms.append(atom)
                oncycle_labels.append(atom.label)

    # we also need everything that the cyclic atoms are bound to
    for atom in oncycle_atoms:
        for neighbour in atom.neighbours:
            if not neighbour in oncycle_atoms:
                if not neighbour in offcycle_atoms:
                    offcycle_atoms.append(neighbour)
                    offcycle_labels.append(neighbour.label)

    # combine cyclic atoms and 1st coordination sphere
    cyclicsystem = oncycle_atoms + offcycle_atoms

    # initialize ringVBO dictionary
    for atom in unassigned:
        if atom.label in oncycle_labels:
            ringVBO[atom] = 0

    # CSD doesn't do periodic boundary conditions, need a workaround
    # check for any periodic copies of cyclic atoms
    for atom in ringcopy.atoms:
        if all([(not atom in oncycle_atoms), (atom.label in oncycle_labels)]):
            if not atom in cyclic_periodic:
                cyclic_periodic.append(atom)
                cyclic_periodic_labels.append(atom.label)
    for atom in cyclic_periodic:
        # print (atom.neighbours)
        for neighbour in atom.neighbours:
            if not neighbour in (offcycle_periodic + cyclic_periodic):
                if not neighbour.label in (oncycle_labels):
                    offcycle_periodic.append(neighbour)

    # remove every atom that isn't part of or directly bound to a cycle
    for atom in ringcopy.atoms:
        if not atom in (cyclicsystem + cyclic_periodic):
            ringcopy.remove_atom(atom)

    # find all non-cyclic bonds
    # single bonds between cycles, break and cap with H
    for bond in ringcopy.bonds:
        if not bond.is_cyclic:
            # bonds between cycles
            if all(
                [
                    (all((member.label in oncycle_labels for member in bond.atoms))),
                    (
                        all(
                            (
                                not member.label in cyclic_periodic_labels
                                for member in bond.atoms
                            )
                        )
                    ),
                ]
            ):
                member1 = bond.atoms[0]
                member2 = bond.atoms[1]
                Hcap1 = molecule.Atom("H", coordinates=member1.coordinates)
                Hcap2 = molecule.Atom("H", coordinates=member2.coordinates)
                Hcap1_id = ringcopy.add_atom(Hcap1)
                Hcap2_id = ringcopy.add_atom(Hcap2)
                ringcopy.add_bond(bond.bond_type, Hcap1_id, member2)
                ringcopy.add_bond(bond.bond_type, Hcap2_id, member1)
                ringcopy.remove_bond(bond)

    # cap off-cycle atoms
    for offatom in offcycle_atoms + offcycle_periodic:
        # get the VBO for each off-cycle atom
        # (VBO with respect to cyclic atoms)
        offVBO = 0

        # quick check for delocalized systems in the ring
        # if there are any, get the delocalised bond orders
        if any(bond.bond_type == "Delocalised" for bond in offatom.bonds):
            offdVBO = delocalisedLBO(ringcopy)
        # add the non-delocalized bond orders
        for bond in offatom.bonds:
            # only interested in bonds to cyclic atoms
            if any(batom.label in oncycle_labels for batom in bond.atoms):
                # Each bond contributes to Ligand Bond Order according to its type
                if bond.bond_type == "Single":
                    offVBO += 1
                elif bond.bond_type == "Double":
                    offVBO += 2
                elif bond.bond_type == "Triple":
                    offVBO += 3
                elif bond.bond_type == "Quadruple":
                    offVBO += 4
                elif bond.bond_type == "Delocalised":
                    offVBO += offdVBO[offatom]
                elif bond.bond_type == "Aromatic":
                    offVBO += 0
                    print("impossible Aromatic bond")
        # cap with appropriate element for VBO
        if offVBO == 1:
            offatom.atomic_symbol = "H"
        elif offVBO == 2:
            offatom.atomic_symbol = "O"
        elif offVBO == 3:
            offatom.atomic_symbol = "N"
        elif offVBO == 4:
            offatom.atomic_symbol = "C"
        elif offVBO == 5:
            offatom.atomic_symbol = "P"
        elif offVBO == 6:
            offatom.atomic_symbol = "S"
        elif offVBO > 6:
            print("no, that's too many")

    # for each cyclic system, reassign bonds, kekulize, and get VBO
    # the bond and atom pruning we did above ensures that fused cycles
    # will be treated as a single system
    # while non-fused cycles that are connected via bonding are treated
    # as seperate systems
    for cyclesys in ringcopy.components:
        # reassign bonds and kekulize
        cyclesys.assign_bond_types()
        cyclesys.kekulize()

        # porhpyrins and similar molecules are misassigned, we will code a hard fix
        # first identify and isolate the inner porphyrin(/like) atoms
        # (these atoms determine overall charge)
        # store these in a dictionary for later
        joining_atoms = dict()
        joining_rings = dict()
        subring_labels = dict()
        porphyrinatoms = dict()
        porphyrin_to_correct = set()
        # ring by ring
        for subring in cyclesys.rings:
            subring_labels[subring] = []
            # get a list of all atom labels in each subring
            for sratom in subring.atoms:
                subring_labels[subring].append(sratom.label)
            # check atom by atom
            for sratom in subring.atoms:
                # check each atom neighbour
                for srneighbour in sratom.neighbours:
                    srn_label = srneighbour.label
                    # if the neighbour is not part of the current ring
                    if not srn_label in subring_labels[subring]:
                        # if the nighbour IS a cyclic atom
                        # consider this a "joining atom"
                        if srn_label in oncycle_labels:
                            try:
                                joining_atoms[srn_label].append(subring)
                            except KeyError:
                                joining_atoms[srn_label] = [subring]
                            try:
                                joining_rings[subring].append(srn_label)
                            except KeyError:
                                joining_rings[subring] = [srn_label]
        for jring in joining_rings:
            if all([(len(jring) == 16), (jring.is_fully_conjugated)]):
                for patom in jring.atoms:
                    plabel = patom.label
                    if not plabel in joining_atoms:
                        ncyclicbonds = 0
                        for pbond in patom.bonds:
                            if pbond.is_cyclic:
                                ncyclicbonds += 1
                        if ncyclicbonds == 2:
                            try:
                                porphyrinatoms[jring].append(patom)
                            except KeyError:
                                porphyrinatoms[jring] = [patom]
        for porph in porphyrinatoms:
            if all(i.atomic_symbol == "N" for i in porphyrinatoms[porph]):
                protonated = 0
                for patom in porphyrinatoms[porph]:
                    if any(i.atomic_symbol == "H" for i in patom.neighbours):
                        protonated += 1
                if protonated == 0:
                    for patom in porphyrinatoms[porph]:
                        porphyrin_to_correct.add(patom.label)

        # quick check for delocalized systems in the ring
        # if there are any, get the delocalised bond orders
        if any(bond.bond_type == "Delocalised" for bond in cyclesys.bonds):
            rdVBO = delocalisedLBO(cyclesys)

        # assign VBO for each on-cycle atom
        for ratom in cyclesys.atoms:
            rVBO = 0
            if ratom.label in oncycle_labels:
                if ratom.label in porphyrin_to_correct:
                    rVBO -= 0.5
                for rbond in ratom.bonds:
                    # Each bond contributes to Ligand Bond Order
                    # according to its type except periodic copies
                    if any(
                        [
                            (rbond.is_cyclic),
                            (
                                not all(
                                    (
                                        mem.label in cyclic_periodic_labels
                                        for mem in rbond.atoms
                                    )
                                )
                            ),
                        ]
                    ):
                        if rbond.bond_type == "Single":
                            rVBO += 1
                        elif rbond.bond_type == "Double":
                            rVBO += 2
                        elif rbond.bond_type == "Triple":
                            rVBO += 3
                        elif rbond.bond_type == "Quadruple":
                            rVBO += 4
                        elif rbond.bond_type == "Delocalised":
                            rVBO += rdVBO[ratom]
                        elif rbond.bond_type == "Aromatic":
                            rVBO += 0
                            print("impossible Aromatic bond")

                # the VBOs are currently associated to atom objects
                # in molecule objects that we have modified
                # we need these to be associated to atom objects in
                # the parent (unmodified) molecule object
                for matom in unassigned:
                    if matom.label == ratom.label:
                        ringVBO[matom] += rVBO
                        # unassigned.remove(matom)
    return ringVBO


def assign_VBS(atom: Atom, rVBO: dict[int, int], dVBO: dict[int, float]) -> int:
    """
    Assigns a Valence-Bond-Sum (VBS) to an atom.

    Parameters:
        atom (ccdc.molecule.Atom): Atom object.
        rVBO (dict[int, int]): dictionary with each atom's index in mole.atoms
                               as keys and VBO (valence bond order) as values.
        dVBO (dict[int, float]): dictionary with delocalized bond-possessing
                                 atom's index in mole.atoms as keys and their
                                 corresponding (delocalized-only) VBS.

    Returns:
        VBO (int): valence bond sum value.
    """
    VBO = 0
    if atom.is_metal:
        return 0
    if atom in rVBO:
        VBO = rVBO[atom]
    else:
        for bond in atom.bonds:
            if any(batom.is_metal for batom in bond.atoms):
                VBO += 0
            # Each bond contributes to Ligand Bond Order according to its type
            elif bond.bond_type == "Single":
                VBO += 1
            elif bond.bond_type == "Double":
                VBO += 2
            elif bond.bond_type == "Triple":
                VBO += 3
            elif bond.bond_type == "Quadruple":
                VBO += 4
            elif bond.bond_type == "Delocalised":
                VBO += dVBO[atom]
            elif bond.bond_type == "Aromatic":
                # necessary? aVBO not defined
                # VBO += aVBO[atom]
                VBO += dVBO[atom]
    return VBO


def delocalisedLBO(molecule: Molecule) -> dict[int, float]:
    """
    Writes a dictionary of all atoms in the molecule with delocalized bonds
    and their (delocalized-only) valence bond sum (VBS).

    Parameters:
        molecule (ccdc.molecule.Molecule): Molecule object.

    Returns:
        delocal_dict (dict[int, float]): dictionary with delocalized bond-possessing
                                        atom's index in mole.atoms as keys and their
                                        corresponding (delocalized-only) VBS.
    """

    def TerminusCounter(atomlist: list[Atom]) -> int:
        """
        Counts the number of termini in the input delocalized bond system.

        Parameters:
            atomlist (list[ccdc.molecule.Atom]): list of atoms in delocalised system.

        Returns:
            NTerminus (int): number of termini in delocalized bond system.
        """
        NTerminus = 0
        for member in atomlist:
            connectivity = 0
            for bond in member.bonds:
                if bond.bond_type == "Delocalised":
                    connectivity += 1
            if connectivity == 1:
                NTerminus += 1
        return NTerminus

    def delocal_crawl(atomlist: list[Atom]) -> list[Atom]:
        """
        Recursively searches for atoms in delocalised bond systems starting from
        an input list containing at least one delocalised bonding atom.

        Parameters:
            atomlist (list[ccdc.molecule.Atom)]: list of atoms in delocalised system.

        Returns:
            atomlist (list[ccdc.molecule.Atom]): modified list of atoms in
                                                 delocalised system.
        """
        for delocatom in atomlist:
            for bond in delocatom.bonds:
                if bond.bond_type == "Delocalised":
                    for member in bond.atoms:
                        if not member in atomlist:
                            atomlist.append(member)
                            return delocal_crawl(atomlist)
        return atomlist

    delocal_dict = {}
    for atom in molecule.atoms:
        if all(
            [
                (any(bond.bond_type == "Delocalised" for bond in atom.bonds)),
                (not atom in delocal_dict),
            ]
        ):
            delocal_dict[atom] = []
            delocal_system = delocal_crawl([atom])
            NTerminus = TerminusCounter(delocal_system)
            for datom in delocal_system:
                connectivity = 0
                delocLBO = 0
                for neighbour in datom.neighbours:
                    if neighbour in delocal_system:
                        connectivity += 1
                if connectivity == 1:
                    # terminus
                    delocLBO = (NTerminus + 1) / NTerminus
                if connectivity > 1:
                    # node
                    delocLBO = (connectivity + 1) / connectivity
                delocal_dict[datom] = delocLBO
    return delocal_dict


def iVBS_FormalCharge(atom: Atom) -> int:
    """
    Determines the formal charge of an atom NOT involved in any aromatic or
    delocalized bonding system.

    Parameters:
        atom (ccdc.molecule.Atom): Atom object

    Returns:
        charge (int): formal charge of the input atom.
    """
    VBO = 0
    if atom.is_metal:
        return VBO
    CN = 0
    for neighbour in atom.neighbours:
        if not neighbour.is_metal:
            CN += 1
    valence = valence_e(atom)
    charge = 0
    for bond in atom.bonds:
        if any(batom.is_metal for batom in bond.atoms):
            VBO += 0
        # Each bond contributes to Ligand Bond Order according to its type
        elif bond.bond_type == "Single":
            VBO += 1
        elif bond.bond_type == "Double":
            VBO += 2
        elif bond.bond_type == "Triple":
            VBO += 3
        elif bond.bond_type == "Quadruple":
            VBO += 4
    # need the unpaired electrons
    unpaired_e = 4 - abs(4 - valence)
    # expanded valences require special handling
    if VBO <= (unpaired_e):
        charge = VBO - unpaired_e
    # Expanded (2e) valences:
    elif (VBO > unpaired_e) and (VBO < valence):
        diff = VBO - unpaired_e
        if diff <= 2:
            UPE = valence - unpaired_e - 2
        elif diff <= 4:
            UPE = valence - unpaired_e - 4
        elif diff <= 6:
            UPE = valence - unpaired_e - 6
        elif diff <= 8:
            UPE = valence - unpaired_e - 8
        charge = valence - (VBO + UPE)
    elif VBO >= (valence):
        charge = valence - VBO
    return charge


def get_CN(atom: Atom) -> int:
    """
    Determines the coordination number of the input atom.

    Parameters:
        atom (ccdc.molecule.Atom): Atom object.

    Returns:
        coord_number (int): Atom's coordination number.
    """
    CN = 0
    for neighbour in atom.neighbours:
        if not neighbour.is_metal:
            CN += 1
    return CN


def valence_e(elmnt: Atom) -> int:
    """
    Determines the number of valence electrons of an atom/element.

    Parameters:
        elmnt (ccdc.molecule.Atom): Atom object.

    Returns:
        valence (int): Atom's valence electron count.
    """
    atom = mendeleev.element(elmnt.atomic_symbol)
    if atom.block == "s":
        valence = atom.group_id
    if atom.block == "p":
        valence = atom.group_id - 10
    if atom.block == "d":
        valence = atom.group_id
    if atom.block == "f":
        if atom.atomic_number in range(56, 72):
            valence = atom.atomic_number - 57 + 3
        elif atom.atomic_number in range(88, 104):
            valence = atom.atomic_number - 89 + 3
    if atom.group_id == 18:
        valence = 8
    if atom.symbol == "He":
        valence = 2
    return valence


def carbocation_check(atom: Atom) -> Literal["tetrahedral", "trigonal"]:
    """
    Check carbocation/carbanion geometry according to bond angles.

    Parameters:
        atom (ccdc.molecule.Atom): Atom object.

    Returns:
        Literal["tetrahedral", "trigonal"]: geometry at input atom.
    """
    abc = []
    # get atom neighbours
    for neighbours in atom.neighbours:
        if not neighbours.is_metal:
            abc.append(neighbours)
    # get all three relevant bond angles
    angle1 = descriptors.MolecularDescriptors.atom_angle(abc[0], atom, abc[1])
    angle2 = descriptors.MolecularDescriptors.atom_angle(abc[0], atom, abc[2])
    angle3 = descriptors.MolecularDescriptors.atom_angle(abc[1], atom, abc[2])
    # average the angels
    AVGangle = abs(angle1 + angle2 + angle3) / 3
    # take the difference between the averaged bond angles and
    # ideal trigonal planar/tetrahedral bond angles
    tet = abs(AVGangle - 109.5)
    trig = abs(AVGangle - 120)
    if tet < trig:
        return "tetrahedral"
    if trig < tet:
        return "trigonal"


def carbene_type(atom: Atom) -> Literal["singlet", "triplet"]:
    """
    Distinguishes between singlet and triplet carbenes.

    Parameters:
        atom (ccdc.molecule.Atom): Atom object(s) suspected of belonging to a
                                   carbene (2-coordinate carbon II).

    Returns:
        Literal["singlet", "triplet"]: carbene type at input atom.
    """
    # get alpha-atoms
    alpha = atom.neighbours
    alpha_type = []
    # get element symbols for alpha atoms
    for a in alpha:
        if not a.is_metal:
            alpha_type.append(a.atomic_symbol)
    # if any alpha atom is a heteroatom, return "singlet"
    # these are Fischer carbenes
    for a in alpha_type:
        if not any([(a == "C"), (a == "H")]):
            return "singlet"
    # if the carbene C is in a heterocycle,
    # return "singlet"
    # there are Arduengo carbenes (NHCs, CAACs)
    if atom.is_cyclic == True:
        for ring in atom.rings:
            for species in ring.atoms:
                if not species.atomic_symbol == "C":
                    return "singlet"
    # for all other carbenes, return "triplet"
    # these are Schrock carbenes
    return "triplet"


def hapticity(atom: Atom, metalsite: list[Atom]) -> bool:
    """
    Determines if a ligand binding site possesses hapticity (any n-hapto).

    Parameters:
        atom (ccdc.molecule.Atom): Atom object.
        metalsites (list[ccdc.molecule.Atom]): list of metal sites in the structure
                                               that belong to the asymmetric unit.

    Returns:
        bool: whether the the input ligand is hapto-.
    """
    for atom2 in atom.neighbours:
        if not atom2.is_metal:
            if any(n2.label == metalsite.label for n2 in atom2.neighbours):
                return True
    return False


def bridging(atom: Atom) -> int:
    """
    Determines how many metal atoms the input atom binds to search for
    bridging sites.

    Parameters:
        atom (ccdc.molecule.Atom): binding site Atom object.

    Returns:
       bridge (int): number of metal atoms bound to the atom.
    """
    bridge = 0
    for n in atom.neighbours:
        if n.is_metal:
            bridge += 1
    return bridge


def iVBS_Oxidation_Contrib(
    unique_atoms: list[Atom], rVBO: dict[int, int], dVBO: dict[int, float]
) -> dict[Atom, float]:
    """
    Determines the oxidation state contribution of all unique atoms.

    Parameters:
        unique_atoms (list[ccdc.molecule.Atom]): unique atoms belonging to the
                                                 asymmetric unit.
        rVBO (dict[int, int]): dictionary with each atom's index in mole.atoms
                               as keys and VBO (valence bond order) as values.
        dVBO (dict[int, float]): dictionary with delocalized bond-possessing atom's
                                index in mole.atoms as keys and their corresponding
                                (delocalized-only) VBS.

    Returns:
        oxi_contrib (dict[ccdc.molecule.Atom, float)]: dictionary with Atom object
                         as keys and their oxidation state contribution as values.
    """
    VBS = 0
    CN = 0
    valence = 0
    oxi_contrib = {}
    # for each unique atom
    for atom in unique_atoms:
        # assign valence-bond-sum
        VBS = assign_VBS(atom, rVBO, dVBO)
        # determine coordination number
        CN = get_CN(atom)
        #  determine number of valence electrons
        valence = valence_e(atom)
        # get number of unpaired electrons in the free element
        unpaired_e = 4 - abs(4 - valence)

        #  metals do not contribute:
        if atom.is_metal:
            oxi_contrib[atom] = 0
        # Normal valences:
        elif VBS <= (unpaired_e):
            oxi_contrib[atom] = unpaired_e - VBS
        # Expanded (2e) valences:
        elif (VBS > unpaired_e) and (VBS < valence):
            diff = VBS - unpaired_e
            if diff <= 2:
                UPE = valence - unpaired_e - 2
            elif diff <= 4:
                UPE = valence - unpaired_e - 4
            elif diff <= 6:
                UPE = valence - unpaired_e - 6
            elif diff <= 8:
                UPE = valence - unpaired_e - 8
            oxi_contrib[atom] = VBS + UPE - valence
        elif VBS >= (valence):
            oxi_contrib[atom] = VBS - valence

        # need to check for 3-coordinate carbocations,
        # 3-coordinate carbanions, carbenes, and heavier
        # homologues (these are not immediately detectable)
        if any(
            [
                (atom.atomic_symbol == "C"),
                (atom.atomic_symbol == "Si"),
                (atom.atomic_symbol == "Ge"),
                (atom.atomic_symbol == "Pb"),
            ]
        ):
            if not atom in rVBO:
                # 3 coordinate and VBS 3 could be
                # carbanion or carbocation
                if VBS == 3 and CN == 3:
                    geom = carbocation_check(atom)
                    if geom == "trigonal":
                        oxi_contrib[atom] = -1
                    if geom == "tetrahedral":
                        oxi_contrib[atom] = 1
            # VBS 2 and 2 coordinate is carbene,
            # but singlet or triplet?
            if VBS == 2 and CN == 2:
                carbene = carbene_type(atom)
                if carbene == "singlet":
                    oxi_contrib[atom] = 0
                if carbene == "triplet":
                    oxi_contrib[atom] = 2

        # Nitro groups frequently have both N-O bonds assigned
        # as double bonds, giving incorrect VBS of 5
        # and oxidation contribution of -2
        # this block catches this and applies a fix
        if all(
            [
                (atom.atomic_symbol == "N"),
                (VBS == 5 and CN == 3),
            ]
        ):
            N_sphere1 = atom.neighbours
            O_count = 0
            for neighbour in N_sphere1:
                if neighbour.atomic_symbol == "O":
                    O_count += 1
            geom = carbocation_check(atom)
            if O_count == 2 and geom == "trigonal":
                oxi_contrib[atom] = 0

    return oxi_contrib


def redundantAON(AON: dict[Atom, float], molecule: Molecule) -> dict[Atom, float]:
    """
    Maps the oxidation contributions of unique atom sites to the redundant atom
    sites according to their shared atom labels.

    Parameters:
        AON (dict[ccdc.molecule.Atom, float]): dictionary with Atom object as keys
                        and their oxidation state contribution as values for unique
                        Atom objects.
        molecule (ccdc.molecule.Molecule): Molecule object.

    Returns:
        redAON (dict[ccdc.molecule.Atom, float]): dictionary with Atom object as keys
                        and their oxidation state contribution as values for all
                        (including redundant) Atom objects.
    """
    redAON = {}
    for rsite1 in molecule.atoms:
        for usite1 in AON:
            redAON[usite1] = AON[usite1]
            if rsite1.label == usite1.label:
                redAON[rsite1] = AON[usite1]
    return redAON


def binding_domain(
    binding_sites: list[Atom],
    AON: dict[Atom, float],
    molecule: Molecule,
    usites: list[Atom],
) -> dict[Atom, list[Atom]]:
    """
    Builds bonding domains within the crystal structure to determine which
    metal binding sites (Atom objects directly bonded to a metal) are connected
    via conjugation. Function accounts for the inconsistent assignment of
    delocalized bonds, by using the bonding domains (see methodology section
    for details on the implementation and validation).

    Parameters:
        binding_sites (list[ccdc.molecule.Atom]): list of binding sites connecting
                                                  metal atoms and ligands.
        AON (dict[ccdc.molecule.Atom, float]): dictionary with Atom object as keys
                        and their oxidation state contribution as values for unique
                        Atom objects.
        molecule (ccdc.molecule.Molecule): Molecule object.
        uniquesites (list[ccdc.molecule.Atom]): list of unique atoms.

    Returns:
        sitedomain (dict[ccdc.molecule.Atom, list[ccdc.molecule.Atom]]):
                        dictionary with Atom object as keys and a list of Atoms
                        connected through bonding that form a binding domain
                        as values.
    """

    def arom_domains(
        site: Atom, usites: list[Atom], aromlist: list[Atom], bondset: list[Bond]
    ) -> list[Atom]:
        """
        Recursively generate aromatic binding domains.
        """
        for bond in site.bonds:
            bondset.add(bond)
        for bond in bondset:
            for member in bond.atoms:
                if all(
                    [
                        (not member in aromlist),
                        (not member.is_metal),
                        (any(mbond.bond_type == "Aromatic" for mbond in member.bonds)),
                    ]
                ):
                    aromlist.append(member)
                    for mbond in member.bonds:
                        bondset.add(mbond)
                    return arom_domains(site, usites, aromlist, bondset)
        # aromlist currently contains non-unique instances of atoms
        # this will cause problems further down the line, so correct
        for index, member in enumerate(aromlist):
            aromlist[index] = usites[member.label]
        return aromlist

    def deloc_domains(
        site: Atom,
        usites: list[Atom],
        AON: dict[Atom, float],
        molecule: Molecule,
        deloclist: list[Atom],
        bondset: list[Bond],
        checked_bonds: list[Bond],
    ) -> list[Atom]:
        """
        Recursively generate delocalised binding domains.
        """
        for bond in site.bonds:
            if not bond in bondset:
                bondset.add(bond)
        for bond in bondset:
            if not bond in checked_bonds:
                for member in bond.atoms:
                    if all(
                        [
                            (not member in deloclist),
                            (not member.is_metal),
                            (
                                not any(
                                    mbond.bond_type == "Aromatic"
                                    for mbond in member.bonds
                                )
                            ),
                            (
                                any(
                                    [
                                        (
                                            len(
                                                molecule.shortest_path_bonds(
                                                    site, member
                                                )
                                            )
                                            <= 2
                                        ),
                                        (bond.bond_type == "Delocalised"),
                                        (bond.is_conjugated),
                                        (
                                            all(
                                                [
                                                    (bond.bond_type == "Single"),
                                                    (not AON[member] == 0),
                                                ]
                                            )
                                        ),
                                        (
                                            all(
                                                [
                                                    (
                                                        not any(
                                                            mbond.bond_type == "Single"
                                                            for mbond in member.bonds
                                                        )
                                                    ),
                                                    (
                                                        not any(
                                                            mbond.bond_type
                                                            == "Aromatic"
                                                            for mbond in member.bonds
                                                        )
                                                    ),
                                                    (
                                                        not any(
                                                            mbond.bond_type
                                                            == "Delocalised"
                                                            for mbond in member.bonds
                                                        )
                                                    ),
                                                ]
                                            )
                                        ),
                                    ]
                                )
                            ),
                        ]
                    ):
                        deloclist.append(member)
                        for mbond in member.bonds:
                            bondset.add(mbond)
                checked_bonds.add(bond)
                return deloc_domains(
                    site, usites, AON, molecule, deloclist, bondset, checked_bonds
                )
        # deloclist currently contains non-unique instances of atoms
        # this will cause problems further down the line, so correct
        for index, member in enumerate(deloclist):
            deloclist[index] = usites[member.label]
        return deloclist

    sitedomain = {}
    for site in binding_sites:
        if not site.is_metal == True:
            if any(sbond.bond_type == "Aromatic" for sbond in site.bonds):
                sitedomain[site] = arom_domains(
                    site, usites, aromlist=[site], bondset=set()
                )
            if not any(sbond.bond_type == "Aromatic" for sbond in site.bonds):
                sitedomain[site] = deloc_domains(
                    site,
                    usites,
                    AON,
                    molecule,
                    deloclist=[site],
                    bondset=set(),
                    checked_bonds=set(),
                )

    for site in sitedomain:
        olapset = set()
        for site2 in sitedomain:
            for member in sitedomain[site]:
                if member in sitedomain[site2]:
                    olapset.add(site2)
        for olap in olapset:
            sitedomain[site] = list(set(sitedomain[site]) | set(sitedomain[olap]))
            sitedomain[olap] = sitedomain[site]
    return sitedomain


def binding_contrib(
    binding_sphere: dict[Atom, list[Atom]],
    binding_sites: list[Atom],
    AON: dict[Atom, float],
) -> dict[Atom, float]:
    """
    Redistributes oxidation state contributions within a binding domain.
    Equal distribution is assumed across connected binding sites in each domain.

    Parameters:
        binding_sphere (dict[ccdc.molecule.Atom, list[ccdc.molecule.Atom]]):
                        dictionary with Atom object as keys and a list of Atoms
                        connected through bonding that form a binding domain as values.
        binding_sites (list[ccdc.molecule.Atom]): list of binding sites connecting
                                                  metal atoms and ligands.
        AON (dict[ccdc.molecule.Atom, float]): dictionary with Atom object as keys
                        and their oxidation state contribution as values for unique Atoms.

    Returns:
        site_contrib (dict[ccdc.molecule.Atom, float]): dictionary with Atom object as keys
                        and their updated oxidation state contribution as values accounting
                        for distribution within the binding domain.
    """
    site_contrib = {}
    for site in binding_sphere:
        site_contrib[site] = 0
        nbinding = 0
        for member in binding_sphere[site]:
            if member in binding_sites:
                nbinding += 1
            site_contrib[site] += AON[member]
        site_contrib[site] /= nbinding
    return site_contrib


def outer_sphere_domain(
    uniquesites: list[Atom], binding_domains: dict[Atom, list[Atom]]
) -> list[Atom]:
    """
    Identifies sites outside of the binding domains which must be checked for
    outer sphere charge contributions.

    Parameters:
        uniquesites (list[ccdc.molecule.Atom]): list of unique atoms in the structure
                                                belonging to the asymmetric unit.
        binding_domains (dict[ccdc.molecule.Atom, list[ccdc.molecule.Atom]]):
                        dictionary with Atom object as keys and a list of Atoms
                        connected through bonding that form a binding domain as values.

    Returns:
        outer_sphere (list[ccdc.molecule.Atom]): list of unique, non-metal atoms
                                                 outside of binding domains.
    """
    outer_sphere = []
    for site in uniquesites:
        if all(
            [
                (
                    not any(
                        site in binding_domains[domain] for domain in binding_domains
                    )
                ),
                (not site.is_metal),
            ]
        ):
            outer_sphere.append(site)
    return outer_sphere


def outer_sphere_contrib(outer_sphere: list[Atom], AON: dict[Atom, float]) -> int:
    """
    Calculates the total oxidation state contribution of the outer sphere atoms as
    the sum of their formal charge/contributions.

    Parameters:
        outer_sphere (list[ccdc.molecule.Atom]): list of unique, non-metal atoms
                                                 outside of binding domains.
        AON (dict[ccdc.molecule.Atom, float]): dictionary with Atom object as keys
                        and their oxidation state contribution as values for unique Atoms.

    Returns:
        contrib (int): sum of outer sphere charge contributions.
    """
    contrib = 0
    for site in outer_sphere:
        contrib += AON[site]
    return contrib


def get_metal_networks(
    ligand_sites: dict[Atom, list[Atom]],
    binding_sphere: dict[Atom, list[Atom]],
    bindingAON: dict[Atom, float],
) -> dict[Atom, list[Atom]]:
    """
    Determines the metal atoms that are connected through binding domains and
    charged ligands. Any connections through neutral ligands are ignored as they
    do not contribute to the charge accounting.

    Parameters:
        ligand_sites (dict[ccdc.molecule.Atom, list[ccdc.molecule.Atom]]):
                            dictionary with metal Atom object as key and the
                            list of ligand atoms which bind them as values.
        binding_sphere (dict[ccdc.molecule.Atom, list[ccdc.molecule.Atom]]):
                            dictionary with Atom object as keys and a list of Atoms
                            connected through bonding that form a binding domain
                            as values.
        bindingAON (dict[ccdc.molecule.Atom, float]): dictionary with Atom object as
                            keys and their updated oxidation state contribution as
                            values accounting for distribution within the binding domain.

    Returns:
        network_dict (dict[ccdc.molecule.Atom, list[ccdc.molecule.Atom]]):
                            dictionary with as metal Atom objects as keys and a list
                            of other metal Atom objects connected through binding
                            domains/charged ligands as values. Ignores neutral ligand
                            connections.
    """

    def network_crawl(
        ligand_sites: dict[Atom, list[Atom]],
        binding_sphere: dict[Atom, list[Atom]],
        bindingAON: dict[Atom, float],
        metal_networks: list[Atom],
        checked_sites: list[Atom],
        group: list[Atom],
    ) -> list[Atom]:
        """
        Recursively crawl through bonds to identify metals connected through direct bondings
        or delocalised/conjugated systems.
        """
        for metal in group:
            # This block will find all metals connected to an input metal by metal-metal bonds
            checked_sites.append(metal)
            for neighbour in metal.neighbours:
                if neighbour.is_metal:
                    if not neighbour in checked_sites:
                        checked_sites.append(neighbour)
                        for site in ligand_sites:
                            if neighbour.label == site.label:
                                if not site in group:
                                    group.append(site)
                        return network_crawl(
                            ligand_sites,
                            binding_sphere,
                            bindingAON,
                            metal_networks,
                            checked_sites,
                            group,
                        )
            # this block will find all metals connected to an input metal by
            # conjugation and delocalized charge ligands
            # metals connected through NEUTRAL ligands will be ignored
            for site in ligand_sites[metal]:
                if all([(not bindingAON[site] == 0), (not site in checked_sites)]):
                    for dsite in binding_sphere[site]:
                        if all(
                            [(not dsite in checked_sites), (dsite in binding_sphere)]
                        ):
                            checked_sites.append(dsite)
                            for environ in dsite.neighbours:
                                if environ.is_metal:
                                    if environ in ligand_sites:
                                        if all(
                                            [
                                                (
                                                    all(
                                                        not environ in network
                                                        for network in metal_networks
                                                    )
                                                ),
                                                (not environ in group),
                                            ]
                                        ):
                                            group.append(environ)
                                    else:
                                        for umetal in ligand_sites:
                                            if all(
                                                [
                                                    (umetal.label == environ.label),
                                                    (
                                                        all(
                                                            not umetal in network
                                                            for network in metal_networks
                                                        )
                                                    ),
                                                    (not umetal in group),
                                                ]
                                            ):
                                                group.append(umetal)
                    return network_crawl(
                        ligand_sites,
                        binding_sphere,
                        bindingAON,
                        metal_networks,
                        checked_sites,
                        group,
                    )
        return group

    metal_networks = []
    for metal in ligand_sites:
        if all(not metal in network for network in metal_networks):
            metal_networks.append(
                network_crawl(
                    ligand_sites,
                    binding_sphere,
                    bindingAON,
                    metal_networks,
                    checked_sites=[],
                    group=[metal],
                )
            )

    network_dict = {}
    for network in metal_networks:
        for metal in network:
            network_dict[metal] = network
    return network_dict


def distribute_ONEC(
    sONEC: dict[Atom, list[float, float]],
    metal_networks: dict[Atom, list[Atom]],
    IEs: dict[str, list[float]],
    ONP: dict[str, list[float]],
    highest_known_ON: dict[str, int],
    metal_CN: dict[Molecule, int],
    most_probable_ON: dict[str, int],
) -> dict[Atom, list[float, float]]:
    """
    Redistributes the oxidation state contributions across all metal atoms in
    the structure according to their metal networks (fully local distribution)
    & calculates their associated electron counts. Features utilizing electron
    counts are minimally implemented at this time.

    Parameters:
        sONEC (dict[ccdc.molecule.Atom, list[float, float]]): dictionary with
                        metal Atom object as keys and lists containing the initial
                        oxidation state and electron count implied by only the
                        equal splitting of binding domain charges as values.
        metal_networks (dict[ccdc.molecule.Atom, list[ccdc.molecule.Atom]]):
                        dictionary with as metal Atom objects as keys and a list of
                        other metal Atom objects connected through binding domains/
                        charged ligands as values. Ignores neutral ligand connections.
        IEs (dict[str, list[float]]): dictionary with metal element symbols as keys
                        and a list of their reported ionization energies as values.
        ONP (dict[str, list[float]]): dictionary with metal element symbols as keys
                        and a list of the probability at the relevant oxidation states
                        as values.
        highest_known_ON (dict[str, int]) : dictionary with metal element symbols as
                        keys and a their highest known oxidation state as values.
        metal_CN (dict[ccdc.molecule.Molecule, int]): dictionary with as metal Atom
                        objects as keys and their effective coordination number as values.
        most_probable_ON (dict[str, int]) : dictionary with metal element symbols as
                        keys and a their oxidation state with the highest probability
                        as values.

    Returns:
        distributed_ONEC (dict[ccdc.molecule.Atom, list[float, float]]): dictionary
                        with metal Atom object as keys and lists containing their
                        redistributed oxidation state and electron count as values.
    """

    def recursive_distributor_single_network(
        iONEC: dict[Atom, list[float, float]],
        available_charge: int,
        sorted_metals: dict[str, Atom],
        IEs: dict[str, list[float]],
        ONP: dict[str, list[float]],
        highest_known_ON: dict[str, int],
    ) -> dict[Atom, list[float, float]]:
        """
        Distribute network charge according to ionization energy and probability
        until all charge is distributed. Performed after tallying available
        network charge and sorting network metals by element type.
        """
        # initialize working dictionary
        dONEC = {}
        dONEC = dict(iONEC)

        # positive contribution?
        if available_charge > 0:

            # get list of improbable and improbable next oxidations
            prob_metal_type = []
            improb_metal_type = []
            for metal_type in sorted_metals:
                try:
                    prob = float(
                        100
                        * ONP[metal_type][
                            math.floor(dONEC[sorted_metals[metal_type][0]][0]) + 1
                        ]
                    )
                except IndexError:
                    prob = 0
                if prob >= 1:
                    prob_metal_type.append(metal_type)
                else:
                    improb_metal_type.append(metal_type)

            # if only one metal type has a probable next oxidation state, do that
            if len(prob_metal_type) == 1:
                lowestMetal = prob_metal_type[0]

            # if more than one metal type has a probable next oxidation state,
            # determine next lowest ionization energy among probable next
            # oxidation states
            elif len(prob_metal_type) > 1:
                # find lowest next ionization energy
                for metal_type in prob_metal_type:
                    # metal in a negative oxidation state? Use IE = 0.
                    if dONEC[sorted_metals[metal_type][0]][0] < 0:
                        currentIE = 0
                    # metal oxidation state at or higher than highest known? Set IE arbitrarily high.
                    elif (
                        dONEC[sorted_metals[metal_type][0]][0]
                        >= highest_known_ON[metal_type]
                    ):
                        currentIE = 9999
                    # otherwise, use the appropriate IE.
                    else:
                        currentIE = float(
                            IEs[metal_type][
                                math.floor(dONEC[sorted_metals[metal_type][0]][0])
                            ]
                        )
                    if not "lowestIE" in locals():
                        lowestIE = currentIE
                        lowestMetal = metal_type
                    else:
                        if currentIE < lowestIE:
                            lowestIE = currentIE
                            lowestMetal = metal_type

            # if there is no probable next oxidation state available,
            # determine lowest ionization energy among improbable next oxidation states
            elif len(prob_metal_type) == 0:
                # find lowest next ionization energy
                for metal_type in improb_metal_type:
                    # metal in a negative oxidation state? Use IE = 0.
                    if dONEC[sorted_metals[metal_type][0]][0] < 0:
                        currentIE = 0
                    # metal oxidation state at or higher than highest known? Set IE arbitrarily high.
                    elif (
                        dONEC[sorted_metals[metal_type][0]][0]
                        >= highest_known_ON[metal_type]
                    ):
                        currentIE = 9999
                    # otherwise, use the appropriate IE.
                    else:
                        currentIE = float(
                            IEs[metal_type][
                                math.floor(dONEC[sorted_metals[metal_type][0]][0])
                            ]
                        )
                    if not "lowestIE" in locals():
                        lowestIE = currentIE
                        lowestMetal = metal_type
                    else:
                        if currentIE < lowestIE:
                            lowestIE = currentIE
                            lowestMetal = metal_type

            # distribute one ionization energy level worth of charge
            if available_charge >= len(sorted_metals[lowestMetal]):
                for metal in sorted_metals[lowestMetal]:
                    dONEC[metal][0] += 1
                    available_charge -= 1
            elif available_charge < len(sorted_metals[lowestMetal]):
                for metal in sorted_metals[lowestMetal]:
                    dONEC[metal][0] += available_charge / (
                        len(sorted_metals[lowestMetal])
                    )
                available_charge = 0

        # negative contribution?
        if available_charge < 0:
            # get list of improbable and improbable next oxidations
            prob_metal_type = []
            improb_metal_type = []
            for metal_type in sorted_metals:
                try:
                    prob = float(
                        100
                        * ONP[metal_type][
                            math.floor(dONEC[sorted_metals[metal_type][0]][0]) - 1
                        ]
                    )
                except IndexError:
                    prob = 0
                if prob >= 1:
                    prob_metal_type.append(metal_type)
                else:
                    improb_metal_type.append(metal_type)

            # if only one metal type has a probable next oxidation state, do that
            if len(prob_metal_type) == 1:
                highestMetal = prob_metal_type[0]

            # if more than one metal type has a probable next oxidation state,
            # determine next highest ionization energy among probable next
            # oxidation states
            elif len(prob_metal_type) > 1:
                for metal_type in prob_metal_type:
                    # metal in a negative oxidation state? Use IE = 0.
                    if dONEC[sorted_metals[metal_type][0]][0] <= 0:
                        currentIE = 0
                    # metal oxidation state higher than highest known? Set IE arbitrarily high.
                    elif (
                        dONEC[sorted_metals[metal_type][0]][0]
                        > highest_known_ON[metal_type]
                    ):
                        currentIE = 9999
                    # otherwise, use the appropriate IE.
                    else:
                        if not any(
                            [
                                (
                                    math.isclose(
                                        (dONEC[sorted_metals[metal_type][0]][0] % 1),
                                        0,
                                        abs_tol=0.0001,
                                    )
                                ),
                                (
                                    math.isclose(
                                        (dONEC[sorted_metals[metal_type][0]][0] % 1),
                                        1,
                                        abs_tol=0.0001,
                                    )
                                ),
                            ]
                        ):
                            currentIE = float(
                                IEs[metal_type][
                                    math.floor(dONEC[sorted_metals[metal_type][0]][0])
                                ]
                            )
                        else:
                            currentIE = float(
                                IEs[metal_type][
                                    math.floor(dONEC[sorted_metals[metal_type][0]][0])
                                    - 1
                                ]
                            )
                    if not "highestIE" in locals():
                        highestIE = currentIE
                        highestMetal = metal_type
                    else:
                        if currentIE > highestIE:
                            highestIE = currentIE
                            highestMetal = metal_type

            # if no probable next oxidation states are available,
            # determine next highest ionization energy among probable next
            # oxidation states
            elif len(improb_metal_type) > 0:
                for metal_type in improb_metal_type:
                    # metal in a negative oxidation state? Use IE = 0.
                    if dONEC[sorted_metals[metal_type][0]][0] <= 0:
                        currentIE = 0
                    # metal oxidation state higher than highest known? Set IE arbitrarily high.
                    elif (
                        dONEC[sorted_metals[metal_type][0]][0]
                        > highest_known_ON[metal_type]
                    ):
                        currentIE = 9999
                    # otherwise, use the appropriate IE.
                    else:
                        if not any(
                            [
                                (
                                    math.isclose(
                                        (dONEC[sorted_metals[metal_type][0]][0] % 1),
                                        0,
                                        abs_tol=0.0001,
                                    )
                                ),
                                (
                                    math.isclose(
                                        (dONEC[sorted_metals[metal_type][0]][0] % 1),
                                        1,
                                        abs_tol=0.0001,
                                    )
                                ),
                            ]
                        ):
                            currentIE = float(
                                IEs[metal_type][
                                    math.floor(dONEC[sorted_metals[metal_type][0]][0])
                                ]
                            )
                        else:
                            currentIE = float(
                                IEs[metal_type][
                                    math.floor(dONEC[sorted_metals[metal_type][0]][0])
                                    - 1
                                ]
                            )
                    if not "highestIE" in locals():
                        highestIE = currentIE
                        highestMetal = metal_type
                    else:
                        if currentIE > highestIE:
                            highestIE = currentIE
                            highestMetal = metal_type

            # distribute one ionization energy level worth of charge
            if (-1 * available_charge) >= len(sorted_metals[highestMetal]):
                for metal in sorted_metals[highestMetal]:
                    dONEC[metal][0] -= 1
                    available_charge += 1
            elif (-1 * available_charge) < len(sorted_metals[highestMetal]):
                for metal in sorted_metals[highestMetal]:
                    dONEC[metal][0] += available_charge / (
                        len(sorted_metals[highestMetal])
                    )
                available_charge = 0

        # if all charge has been distributed, we're done, otherwise, roll again
        if available_charge == 0:
            return dONEC
        else:
            return recursive_distributor_single_network(
                dONEC, available_charge, sorted_metals, IEs, ONP, highest_known_ON
            )

    # operate on each network individually
    distributed_ONEC = {}
    for network in metal_networks:

        # sort metals by element type
        sorted_metals = {}
        for metal in metal_networks[network]:
            sorted_metals[metal.atomic_symbol] = []
        for metal in metal_networks[network]:
            sorted_metals[metal.atomic_symbol].append(metal)

        # tally up network charge to be distributed
        # and initialize metals to most probable ON
        # (adjust network charge accordingly)
        network_charge = 0
        for metal in metal_networks[network]:
            network_charge += sONEC[metal][0]
            distributed_ONEC[metal] = [most_probable_ON[metal.atomic_symbol]]
            network_charge -= most_probable_ON[metal.atomic_symbol]
            distributed_ONEC[metal].append(int(sONEC[metal][1]))

        # if the most probable oxidation distribution has already balanced the charge, we're done
        # if not, recursively distribute network charge according to ionization energy
        if not (math.isclose(network_charge, 0, abs_tol=0.0001)):
            distributed_ONEC = recursive_distributor_single_network(
                distributed_ONEC,
                network_charge,
                sorted_metals,
                IEs,
                ONP,
                highest_known_ON,
            )

        # finally, adjust electron count to new oxidation state (OiL RiG)
        for metal in metal_networks[network]:
            distributed_ONEC[metal][1] = (
                valence_e(metal) + (2 * metal_CN[metal]) - distributed_ONEC[metal][0]
            )

    return distributed_ONEC


def distribute_OuterSphere(
    sONEC: dict[Atom, list[float, float]],
    outer_sphere_charge: int,
    IEs: dict[Atom, list[Atom]],
    ONP: dict[str, list[float]],
    highest_known_ON: dict[str, int],
    metal_CN: dict[Molecule, int],
) -> dict[Atom, list[float, float]]:
    """
    Redistributes the oxidation state contributions across all metal atoms in
    the structure according to the outer sphere charge contribution (partially
    local distribution) & calculates their associated electron counts.
    Features utilizing electron counts are minimally implemented at this time.

    Parameters:
        sONEC (dict[ccdc.molecule.Atom, list[float, float]]): dictionary with
                        metal Atom object as keys and lists containing the initial
                        oxidation state and electron count implied by only the equal
                        splitting of binding domain charges as values.
        outer_sphere_charge (int): sum of outer sphere charge contributions.
        IEs (dict[str, list(float)]): dictionary with metal element symbols as keys
                        and a list of  their reported ionization energies as values.
        ONP (dict[str, list(float)]): dictionary with metal element symbols as keys
                        and a list of the probability at the relevant oxidation
                        states as values.
        highest_known_ON (dict[str, int]) : dictionary with metal element symbols
                        as keys and a their highest known oxidation state as values.
        metal_CN (dict[ccdc.molecule.Molecule, int]): dictionary with as metal Atom
                        objects as keys and their  effective coordination number
                        as values.

    Returns:
        distributed_ONEC (dict[ccdc.molecule.Atom, list[float, float]]): dictionary
                        with metal Atom object as keys and lists containing their
                        redistributed oxidation state and electron count as values.
    """

    def recursive_distributor(
        iONEC: dict[Atom, list[float, float]],
        available_charge: int,
        IEs: dict[str, list[float]],
        ONP: dict[str, list[float]],
        highest_known_ON: dict[str, int],
    ) -> dict[Atom, list[float, float]]:
        """
        Distribute network charge according to ionization energy and highest
        allowable oxidation state until all charge is distributed. Performed after
        tallying available network charge and sorting network metals by element type.
        """
        # initialize working dictionary
        dONEC = {}
        dONEC = dict(iONEC)

        # positive contribution?
        if available_charge > 0:
            # get list of probable and improbable next oxidations
            prob_metals = []
            improb_metals = []
            for metal in dONEC:
                try:
                    prob = float(
                        100 * ONP[metal.atomic_symbol][math.floor(dONEC[metal][0]) + 1]
                    )
                except IndexError:
                    prob = 0
                if prob >= 1:
                    prob_metals.append(metal)
                else:
                    improb_metals.append(metal)

            if len(prob_metals) == 1:
                lowestMetals = prob_metals
            elif len(prob_metals) > 1:
                for metal in prob_metals:
                    # metal in a negative oxidation state? Use IE = 0.
                    if dONEC[metal][0] < 0:
                        currentIE = 0
                    # metal oxidation state at or higher than highest known? Set IE arbitrarily high.
                    elif dONEC[metal][0] >= highest_known_ON[metal.atomic_symbol]:
                        currentIE = 9999
                    # otherwise, use the appropriate IE.
                    else:
                        currentIE = float(
                            IEs[metal.atomic_symbol][math.floor(dONEC[metal][0])]
                        )
                    if not "lowestIE" in locals():
                        lowestIE = currentIE
                        lowestMetals = [metal]
                    else:
                        if currentIE == lowestIE:
                            lowestMetals.append(metal)
                        if currentIE < lowestIE:
                            lowestIE = currentIE
                            lowestMetals = [metal]
            elif len(prob_metals) == 0:
                for metal in improb_metals:
                    # metal in a negative oxidation state? Use IE = 0.
                    if dONEC[metal][0] < 0:
                        currentIE = 0
                    # metal oxidation state at or higher than highest known? Set IE arbitrarily high.
                    elif dONEC[metal][0] >= highest_known_ON[metal.atomic_symbol]:
                        currentIE = 9999
                    # otherwise, use the appropriate IE.
                    else:
                        currentIE = float(
                            IEs[metal.atomic_symbol][math.floor(dONEC[metal][0])]
                        )
                    if not "lowestIE" in locals():
                        lowestIE = currentIE
                        lowestMetals = [metal]
                    else:
                        if currentIE == lowestIE:
                            lowestMetals.append(metal)
                        if currentIE < lowestIE:
                            lowestIE = currentIE
                            lowestMetals = [metal]
            # distribute one ionization energy level worth of charge
            if available_charge >= len(lowestMetals):
                for metal in lowestMetals:
                    dONEC[metal][0] += 1
                    available_charge -= 1
            elif available_charge < len(lowestMetals):
                for metal in lowestMetals:
                    dONEC[metal][0] += available_charge / (len(lowestMetals))
                available_charge = 0

        # negative contribution?
        if available_charge < 0:
            # get list of improbable and improbable next oxidations
            prob_metals = []
            improb_metals = []
            for metal in dONEC:
                try:
                    prob = float(
                        100 * ONP[metal.atomic_symbol][math.floor(dONEC[metal][0]) - 1]
                    )
                except IndexError:
                    prob = 0
                if prob >= 1:
                    prob_metals.append(metal)
                else:
                    improb_metals.append(metal)

            if len(prob_metals) == 1:
                highestMetals = prob_metals
            elif len(prob_metals) > 1:
                for metal in prob_metals:
                    # metal in a negative oxidation state? Use IE = 0.
                    if dONEC[metal][0] <= 0:
                        currentIE = 0
                    # metal oxidation state higher than highest known? Set IE arbitrarily high.
                    elif dONEC[metal][0] > highest_known_ON[metal.atomic_symbol]:
                        currentIE = 9999
                    # otherwise, use the appropriate IE.
                    else:
                        if not any(
                            [
                                (
                                    math.isclose(
                                        (dONEC[metal][0] % 1), 0, abs_tol=0.0001
                                    )
                                ),
                                (
                                    math.isclose(
                                        (dONEC[metal][0] % 1), 1, abs_tol=0.0001
                                    )
                                ),
                            ]
                        ):
                            currentIE = float(
                                IEs[metal.atomic_symbol][math.floor(dONEC[metal][0])]
                            )
                        else:
                            currentIE = float(
                                IEs[metal.atomic_symbol][
                                    math.floor(dONEC[metal][0]) - 1
                                ]
                            )
                    if not "highestIE" in locals():
                        highestIE = currentIE
                        highestMetals = [metal]
                    else:
                        if currentIE == highestIE:
                            highestMetals.append(metal)
                        if currentIE > highestIE:
                            highestIE = currentIE
                            highestMetals = [metal]

            else:
                for metal in improb_metals:
                    # metal in a negative oxidation state? Use IE = 0.
                    if dONEC[metal][0] <= 0:
                        currentIE = 0
                    # metal oxidation state higher than highest known? Set IE arbitrarily high.
                    elif dONEC[metal][0] > highest_known_ON[metal.atomic_symbol]:
                        currentIE = 9999
                    # otherwise, use the appropriate IE.
                    else:
                        if not any(
                            [
                                (
                                    math.isclose(
                                        (dONEC[metal][0] % 1), 0, abs_tol=0.0001
                                    )
                                ),
                                (
                                    math.isclose(
                                        (dONEC[metal][0] % 1), 1, abs_tol=0.0001
                                    )
                                ),
                            ]
                        ):
                            currentIE = float(
                                IEs[metal.atomic_symbol][math.floor(dONEC[metal][0])]
                            )
                        else:
                            currentIE = float(
                                IEs[metal.atomic_symbol][
                                    math.floor(dONEC[metal][0]) - 1
                                ]
                            )
                    if not "highestIE" in locals():
                        highestIE = currentIE
                        highestMetals = [metal]
                    else:
                        if currentIE == highestIE:
                            highestMetals.append(metal)
                        if currentIE > highestIE:
                            highestIE = currentIE
                            highestMetals = [metal]

            # distribute one ionization energy level worth of charge
            if (-1 * available_charge) >= len(highestMetals):
                for metal in highestMetals:
                    dONEC[metal][0] -= 1
                    available_charge += 1
            elif (-1 * available_charge) < len(highestMetals):
                for metal in highestMetals:
                    dONEC[metal][0] += available_charge / (len(highestMetals))
                available_charge = 0

        # if all charge has been distributed, we're done, otherwise, roll again
        if available_charge == 0:
            return dONEC
        else:
            return recursive_distributor(
                dONEC, available_charge, IEs, ONP, highest_known_ON
            )

    if outer_sphere_charge == 0:
        return sONEC

    distributed_ONEC = {}

    # initialize dictionary for charge distribution
    for metal in sONEC:
        distributed_ONEC[metal] = [sONEC[metal][0]]
        distributed_ONEC[metal].append(int(sONEC[metal][1]))

    # recursively distribute network charge according to ionization energy
    distributed_ONEC = recursive_distributor(
        distributed_ONEC, outer_sphere_charge, IEs, ONP, highest_known_ON
    )

    # finally, adjust electron count to new oxidation state (OiL RiG)
    for metal in sONEC:
        distributed_ONEC[metal][1] = (
            valence_e(metal) + (2 * metal_CN[metal]) - distributed_ONEC[metal][0]
        )

    return distributed_ONEC


def global_charge_distribution(
    metalONdict: dict[Atom, list[float, float]],
    IEs: dict[str, list[float]],
    ONP: dict[str, list[float]],
    highest_known_ON: dict[str, int],
    metal_CN: dict[Molecule, int],
    most_probable_ON: dict[str, int],
) -> dict[Atom, list[float, float]]:
    """
    Redistributes the oxidation state contributions across all metal atoms in
    the structure according to full/global shating (fully delocalized distribution)
    & calculates their associated electron counts. Features utilizing electron
    counts are minimally implemented at this time.

    Parameters:
        metalONdict (dict[ccdc.molecule.Atom, list[float, float]]): dictionary with
                        metal Atom object as keys and lists containing the initial
                        oxidation state and electron count implied by only the
                        equal splitting of binding domain charges as values.
        IEs (dict[str, list[float]]): dictionary with metal element symbols as keys
                        and a list of their reported ionization energies as values.
        ONP (dict[str, list[float]]): dictionary with metal element symbols as keys
                        and a list of the probability at the relevant oxidation
                        states as values.
        highest_known_ON (dict[str, int]) : dictionary with metal element symbols
                        as keys and a their highest known oxidation state as values.
        metal_CN (dict[ccdc.molecule.Molecule, int]): dictionary with as metal Atom
                        objects as keys and their effective coordination number
                        as values.
        most_probable_ON (dict[str, int]) : dictionary with metal element symbols
                        as keys and a their oxidation state with the highest
                        probability as values.

    Returns:
        global_ONEC (dict[ccdc.molecule.Atom, list[float, float]]): dictionary with
                        metal Atom object  as keys and lists containing their
                        redistributed oxidation state and electron count as values.
    """
    global_ONEC = {}

    def recursive_distributor_global(
        iONEC: dict[Atom, list[float, float]],
        available_charge: int,
        sorted_metals: dict[str, Atom],
        IEs: dict[str, list[float]],
        ONP: dict[str, list[float]],
        highest_known_ON: dict[str, int],
    ) -> dict[Atom, list[float, float]]:
        """
        Distribute network charge according to ionization energy and probability
        until all charge is distributed. Performed after tallying available network
        charge and sorting network metals by element type.
        """
        # initialize working dictionary
        dONEC = {}
        dONEC = dict(iONEC)

        # positive contribution?
        if available_charge > 0:

            # get list of improbable and improbable next oxidations
            prob_metal_type = []
            improb_metal_type = []
            for metal_type in sorted_metals:
                try:
                    prob = float(
                        100
                        * ONP[metal_type][
                            math.floor(dONEC[sorted_metals[metal_type][0]][0]) + 1
                        ]
                    )
                except IndexError:
                    prob = 0
                if prob >= 1:
                    prob_metal_type.append(metal_type)
                else:
                    improb_metal_type.append(metal_type)

            # if only one metal type has a probable next oxidation state, do that
            if len(prob_metal_type) == 1:
                lowestMetal = prob_metal_type[0]

            # if more than one metal type has a probable next oxidation state,
            # determine next lowest ionization energy among probable next
            # oxidation states
            elif len(prob_metal_type) > 1:
                # find lowest next ionization energy
                for metal_type in prob_metal_type:
                    # metal in a negative oxidation state? Use IE = 0.
                    if dONEC[sorted_metals[metal_type][0]][0] < 0:
                        currentIE = 0
                    # metal oxidation state at or higher than highest known? Set IE arbitrarily high.
                    elif (
                        dONEC[sorted_metals[metal_type][0]][0]
                        >= highest_known_ON[metal_type]
                    ):
                        currentIE = 9999
                    # otherwise, use the appropriate IE.
                    else:
                        currentIE = float(
                            IEs[metal_type][
                                math.floor(dONEC[sorted_metals[metal_type][0]][0])
                            ]
                        )
                    if not "lowestIE" in locals():
                        lowestIE = currentIE
                        lowestMetal = metal_type
                    else:
                        if currentIE < lowestIE:
                            lowestIE = currentIE
                            lowestMetal = metal_type

            # if there is no probable next oxidation state available,
            # determine lowest ionization energy among improbable next oxidation states
            elif len(prob_metal_type) == 0:
                # find lowest next ionization energy
                for metal_type in improb_metal_type:
                    # metal in a negative oxidation state? Use IE = 0.
                    if dONEC[sorted_metals[metal_type][0]][0] < 0:
                        currentIE = 0
                    # metal oxidation state at or higher than highest known? Set IE arbitrarily high.
                    elif (
                        dONEC[sorted_metals[metal_type][0]][0]
                        >= highest_known_ON[metal_type]
                    ):
                        currentIE = 9999
                    # otherwise, use the appropriate IE.
                    else:
                        currentIE = float(
                            IEs[metal_type][
                                math.floor(dONEC[sorted_metals[metal_type][0]][0])
                            ]
                        )
                    if not "lowestIE" in locals():
                        lowestIE = currentIE
                        lowestMetal = metal_type
                    else:
                        if currentIE < lowestIE:
                            lowestIE = currentIE
                            lowestMetal = metal_type

            # distribute one ionization energy level worth of charge
            if available_charge >= len(sorted_metals[lowestMetal]):
                for metal in sorted_metals[lowestMetal]:
                    dONEC[metal][0] += 1
                    available_charge -= 1
            elif available_charge < len(sorted_metals[lowestMetal]):
                for metal in sorted_metals[lowestMetal]:
                    dONEC[metal][0] += available_charge / (
                        len(sorted_metals[lowestMetal])
                    )
                available_charge = 0

        # negative contribution?
        if available_charge < 0:
            # get list of improbable and improbable next oxidations
            prob_metal_type = []
            improb_metal_type = []
            for metal_type in sorted_metals:
                try:
                    prob = float(
                        100
                        * ONP[metal_type][
                            math.floor(dONEC[sorted_metals[metal_type][0]][0]) - 1
                        ]
                    )
                except IndexError:
                    prob = 0
                if prob >= 1:
                    prob_metal_type.append(metal_type)
                else:
                    improb_metal_type.append(metal_type)

            # if only one metal type has a probable next oxidation state, do that
            if len(prob_metal_type) == 1:
                highestMetal = prob_metal_type[0]

            # if more than one metal type has a probable next oxidation state,
            # determine next highest ionization energy among probable next
            # oxidation states
            elif len(prob_metal_type) > 1:
                for metal_type in prob_metal_type:
                    # metal in a negative oxidation state? Use IE = 0.
                    if dONEC[sorted_metals[metal_type][0]][0] <= 0:
                        currentIE = 0
                    # metal oxidation state higher than highest known? Set IE arbitrarily high.
                    elif (
                        dONEC[sorted_metals[metal_type][0]][0]
                        > highest_known_ON[metal_type]
                    ):
                        currentIE = 9999
                    # otherwise, use the appropriate IE.
                    else:
                        if not any(
                            [
                                (
                                    math.isclose(
                                        (dONEC[sorted_metals[metal_type][0]][0] % 1),
                                        0,
                                        abs_tol=0.0001,
                                    )
                                ),
                                (
                                    math.isclose(
                                        (dONEC[sorted_metals[metal_type][0]][0] % 1),
                                        1,
                                        abs_tol=0.0001,
                                    )
                                ),
                            ]
                        ):
                            currentIE = float(
                                IEs[metal_type][
                                    math.floor(dONEC[sorted_metals[metal_type][0]][0])
                                ]
                            )
                        else:
                            currentIE = float(
                                IEs[metal_type][
                                    math.floor(dONEC[sorted_metals[metal_type][0]][0])
                                    - 1
                                ]
                            )
                    if not "highestIE" in locals():
                        highestIE = currentIE
                        highestMetal = metal_type
                    else:
                        if currentIE > highestIE:
                            highestIE = currentIE
                            highestMetal = metal_type

            # if no probable next oxidation states are available,
            # determine next highest ionization energy among probable next
            # oxidation states
            elif len(improb_metal_type) > 0:
                for metal_type in improb_metal_type:
                    # metal in a negative oxidation state? Use IE = 0.
                    if dONEC[sorted_metals[metal_type][0]][0] <= 0:
                        currentIE = 0
                    # metal oxidation state higher than highest known? Set IE arbitrarily high.
                    elif (
                        dONEC[sorted_metals[metal_type][0]][0]
                        > highest_known_ON[metal_type]
                    ):
                        currentIE = 9999
                    # otherwise, use the appropriate IE.
                    else:
                        if not any(
                            [
                                (
                                    math.isclose(
                                        (dONEC[sorted_metals[metal_type][0]][0] % 1),
                                        0,
                                        abs_tol=0.0001,
                                    )
                                ),
                                (
                                    math.isclose(
                                        (dONEC[sorted_metals[metal_type][0]][0] % 1),
                                        1,
                                        abs_tol=0.0001,
                                    )
                                ),
                            ]
                        ):
                            currentIE = float(
                                IEs[metal_type][
                                    math.floor(dONEC[sorted_metals[metal_type][0]][0])
                                ]
                            )
                        else:
                            currentIE = float(
                                IEs[metal_type][
                                    math.floor(dONEC[sorted_metals[metal_type][0]][0])
                                    - 1
                                ]
                            )
                    if not "highestIE" in locals():
                        highestIE = currentIE
                        highestMetal = metal_type
                    else:
                        if currentIE > highestIE:
                            highestIE = currentIE
                            highestMetal = metal_type

            # distribute one ionization energy level worth of charge
            if (-1 * available_charge) >= len(sorted_metals[highestMetal]):
                for metal in sorted_metals[highestMetal]:
                    dONEC[metal][0] -= 1
                    available_charge += 1
            elif (-1 * available_charge) < len(sorted_metals[highestMetal]):
                for metal in sorted_metals[highestMetal]:
                    dONEC[metal][0] += available_charge / (
                        len(sorted_metals[highestMetal])
                    )
                available_charge = 0

        # if all charge has been distributed, we're done, otherwise, roll again
        if available_charge == 0:
            return dONEC
        else:
            return recursive_distributor_global(
                dONEC, available_charge, sorted_metals, IEs, ONP, highest_known_ON
            )

    # sort metals by element type
    sorted_metals = {}
    for metal in metalONdict:
        sorted_metals[metal.atomic_symbol] = []
    for metal in metalONdict:
        sorted_metals[metal.atomic_symbol].append(metal)

    # tally up the global charge to be distributed
    # and initialize global ON to most probable for all metals
    global_charge = 0
    for metal in metalONdict:
        global_charge += metalONdict[metal][0]
        global_ONEC[metal] = [most_probable_ON[metal.atomic_symbol]]
        global_charge -= most_probable_ON[metal.atomic_symbol]
        global_ONEC[metal].append(int(metalONdict[metal][1]))

    # recursively distribute network charge according to ionization energy
    if math.isclose(global_charge, 0, abs_tol=0.0001):
        distributed_ONEC = global_ONEC
    else:
        distributed_ONEC = recursive_distributor_global(
            global_ONEC, global_charge, sorted_metals, IEs, ONP, highest_known_ON
        )

    # finally, adjust electron count to new oxidation state (OiL RiG)
    for metal in metalONdict:
        global_ONEC[metal][1] = (
            valence_e(metal) + (2 * metal_CN[metal]) - distributed_ONEC[metal][0]
        )
    return global_ONEC


def KnownONs() -> dict[str, list[int]]:
    """
    Reads in the known oxidation states for each metal element.


    Returns:
        KONs (dict[str, list(int)]) : dictionary with metal element symbols as keys
                              and a list of their known oxidation states as values.
    """
    
    data_str = """
            H	-1	1										
            He												
            Li	1											
            Be	0	1	2									
            B	-5	-1	0	1	2	3						
            C	-4	-3	-2	-1	0	1	2	3	4			
            N	-3	-2	-1	1	2	3	4	5				
            O	-2	-1	0	1	2							
            F	-1	0										
            Ne												
            Na	-1	1										
            Mg	0	1	2									
            Al	-2	-1	1	2	3							
            Si	-4	-3	-2	-1	0	1	2	3	4			
            P	-3	-2	-1	0	1	2	3	4	5			
            S	-2	-1	0	1	2	3	4	5	6			
            Cl	-1	1	2	3	4	5	6	7				
            Ar	0											
            K	-1	1										
            Ca	0	1	2									
            Sc	0	1	2	3								
            Ti	-2	-1	0	1	2	3	4					
            V	-3	1	0	1	2	3	4	5				
            Cr	-4	-2	-1	0	1	2	3	4	5	6		
            Mn	-3	-2	-1	0	1	2	3	4	5	6	7	
            Fe	-4	-2	-1	0	1	2	3	4	5	6	7	
            Co	-3	-1	0	1	2	3	4	5				
            Ni	-2	-1	0	1	2	3	4					
            Cu	-2	0	1	2	3	4						
            Zn	-2	0	1	2								
            Ga	-5	-4	-3	-2	-1	1	2	3				
            Ge	-3	-2	-1	0	1	2	3	4				
            As	-3	-2	-1	0	1	2	3	4	5			
            Se	-2	-1	1	2	3	4	5	6				
            Br	-1	1	3	4	5	7						
            Kr	0	1	2									
            Rb	-1	1										
            Sr	0	1	2									
            Y	0	1	2	3								
            Zr	-2	0	1	2	3	4						
            Nb	-3	-1	0	1	2	3	4	5				
            Mo	-4	-2	-1	0	1	2	3	4	5	6		
            Tc	-3	-1	0	1	2	3	4	5	6	7		
            Ru	-4	-2	0	1	2	3	4	5	6	7	8	
            Rh	-3	-1	0	1	2	3	4	5	6			
            Pd	0	1	2	3	4							
            Ag	-2	-1	1	2	3							
            Cd	-2	1	2									
            In	-5	-2	-2	1	2	3						
            Sn	-4	-3	-2	-1	0	1	2	3	4			
            Sb	-3	-2	-1	0	1	2	3	4	5			
            Te	-2	-1	1	3	4	5	6	7				
            I	-1	1	3	4	5	6	7					
            Xe	0	2	4	6	8							
            Cs	-1	1										
            Ba	0	1	2									
            La	0	1	2	3								
            Ce	0	2	3	4								
            Pr	0	1	2	3	4	5						
            Nd	0	2	3	4								
            Pm	2	3										
            Sm	0	2	3									
            Eu	0	2	3									
            Gd	0	1	2	3								
            Tb	0	1	2	3	4							
            Dy	0	2	3	4								
            Ho	0	2	3									
            Er	0	2	3									
            Tm	0	2	3									
            Yb	0	2	3									
            Lu	0	2	3									
            Hf	-2	0	1	2	3	4						
            Ta	-3	-1	0	1	2	3	4	5				
            W	-4	-2	-1	0	1	2	3	4	5	6		
            Re	-3	-1	0	1	2	3	4	5	6	7		
            Os	-4	-2	-1	0	1	2	3	4	5	6	7	8
            Ir	-3	-1	0	1	2	3	4	5	6	7	8	9
            Pt	-3	-2	-1	0	1	2	3	4	5	6		
            Au	-3	-2	-1	0	1	2	3	5				
            Hg	-2	1	2									
            Tl	-5	-2	-1	1	2	3						
            Pb	-4	-2	-1	1	2	3	4					
            Bi	-3	-2	-1	1	2	3	4	5				
            Po	-2	2	4	5	6							
            At	-1	1	3	5	7							
            Rn	2	6										
            Fr	1											
            Ra	2											
            Ac	2	3										
            Th	1	2	3	4								
            Pa	3	4	5									
            U	1	2	3	4	5	6						
            Np	2	3	4	5	6	7						
            Pu	2	3	4	5	6	7	8					
            Am	2	3	4	5	6	7						
            Cm	3	4	5	6								
            Bk	2	3	4	5								
            Cf	2	3	4	5								
            Es	2	3	4									
            Fm	2	3										
            Md	2	3										
            No	2	3										
            Lr	3											
            Rf	4											
            Db	5											
            Sg	0	6										
            Bh	7											
            Hs	8											
            Mt												
            Ds												
            Rg												
            Cn	2											
            Nh												
            Fl												
            Mc												
            Lv												
            Ts												
            Og												
            """.strip()

    KONs: dict[str, list[int]] = {}
    for line in data_str.splitlines():
        parts = line.split()
        elem = parts[0]
        # only parse lines with at least one oxidation state
        if len(parts) > 1:
            # map all remaining tokens to integers
            states = [int(x) for x in parts[1:]]
            KONs[elem] = states
        else:
            # no data: use empty list (or choose a default if you prefer)
            KONs[elem] = []

    return KONs

def IonizationEnergies() -> dict[str, list[float]]:
    """
    Reads in the reported ionization energies for each metal element.

    Returns:
        KIEs (dict[str, list[float]]): dictionary with metal element symbols as keys
                          and a list of their reported ionization energies as values.
    """
    data_str = """
            0	H	13.5984346
            0	He	24.58738901
            1	He	54.41776549
            0	Li	5.391714996
            1	Li	75.640097
            2	Li	122.4543591
            0	Be	9.322699
            1	Be	18.21115
            2	Be	153.896205
            3	Be	217.7185861
            0	B	8.298019
            1	B	25.15483
            2	B	37.93059
            3	B	259.3715
            4	B	340.2260229
            0	C	11.260288
            1	C	24.383154
            2	C	47.88778
            3	C	64.49352
            4	C	392.090518
            5	C	489.993198
            0	N	14.53413
            1	N	29.60125
            2	N	47.4453
            3	N	77.4735
            4	N	97.8901
            5	N	552.06733
            6	N	667.046121
            0	O	13.618055
            1	O	35.12112
            2	O	54.93554
            3	O	77.4135
            4	O	113.899
            5	O	138.1189
            6	O	739.32683
            7	O	871.409883
            0	F	17.42282
            1	F	34.97081
            2	F	62.70798
            3	F	87.175
            4	F	114.249
            5	F	157.16311
            6	F	185.1868
            7	F	953.89805
            8	F	1103.11748
            0	Ne	21.564541
            1	Ne	40.96297
            2	Ne	63.4233
            3	Ne	97.19
            4	Ne	126.247
            5	Ne	157.934
            6	Ne	207.271
            7	Ne	239.097
            8	Ne	1195.80784
            9	Ne	1362.19916
            0	Na	5.13907696
            1	Na	47.28636
            2	Na	71.62
            3	Na	98.936
            4	Na	138.404
            5	Na	172.23
            6	Na	208.504
            7	Na	264.192
            8	Na	299.856
            9	Na	1465.134502
            0	Mg	7.646236
            1	Mg	15.035271
            2	Mg	80.1436
            3	Mg	109.2654
            4	Mg	141.33
            5	Mg	186.76
            6	Mg	225.02
            7	Mg	265.924
            8	Mg	327.99
            9	Mg	367.489
            0	Al	5.985769
            1	Al	18.82855
            2	Al	28.447642
            3	Al	119.9924
            4	Al	153.8252
            5	Al	190.49
            6	Al	241.76
            7	Al	284.64
            8	Al	330.21
            9	Al	398.65
            0	Si	8.15168
            1	Si	16.34585
            2	Si	33.493
            3	Si	45.14179
            4	Si	166.767
            5	Si	205.279
            6	Si	246.57
            7	Si	303.59
            8	Si	351.28
            9	Si	401.38
            0	P	10.486686
            1	P	19.76949
            2	P	30.20264
            3	P	51.44387
            4	P	65.02511
            5	P	220.43
            6	P	263.57
            7	P	309.6
            8	P	372.31
            9	P	424.4
            0	S	10.36001
            1	S	23.33788
            2	S	34.86
            3	S	47.222
            4	S	72.5945
            5	S	88.0529
            6	S	280.954
            7	S	328.794
            8	S	379.84
            9	S	447.7
            0	Cl	12.967633
            1	Cl	23.81364
            2	Cl	39.8
            3	Cl	53.24
            4	Cl	67.68
            5	Cl	96.94
            6	Cl	114.2013
            7	Cl	348.306
            8	Cl	400.851
            9	Cl	456.7
            0	Ar	15.7596119
            1	Ar	27.62967
            2	Ar	40.735
            3	Ar	59.58
            4	Ar	74.84
            5	Ar	91.29
            6	Ar	124.41
            7	Ar	143.4567
            8	Ar	422.6
            9	Ar	479.76
            0	K	4.34066373
            1	K	31.625
            2	K	45.8031
            3	K	60.917
            4	K	82.66
            5	K	99.44
            6	K	117.56
            7	K	154.87
            8	K	175.8174
            9	K	503.67
            0	Ca	6.11315547
            1	Ca	11.871719
            2	Ca	50.91316
            3	Ca	67.2732
            4	Ca	84.34
            5	Ca	108.78
            6	Ca	127.21
            7	Ca	147.24
            8	Ca	188.54
            9	Ca	211.275
            0	Sc	6.56149
            1	Sc	12.79977
            2	Sc	24.756839
            3	Sc	73.4894
            4	Sc	91.95
            5	Sc	110.68
            6	Sc	137.99
            7	Sc	158.08
            8	Sc	180.03
            9	Sc	225.18
            0	Ti	6.82812
            1	Ti	13.5755
            2	Ti	27.49171
            3	Ti	43.26717
            4	Ti	99.299
            5	Ti	119.533
            6	Ti	140.68
            7	Ti	170.5
            8	Ti	192.1
            9	Ti	215.92
            0	V	6.746187
            1	V	14.634
            2	V	29.3111
            3	V	46.709
            4	V	65.28165
            5	V	128.125
            6	V	150.72
            7	V	173.55
            8	V	206
            9	V	230.5
            0	Cr	6.76651
            1	Cr	16.486305
            2	Cr	30.959
            3	Cr	49.16
            4	Cr	69.46
            5	Cr	90.6349
            6	Cr	160.29
            7	Cr	184.76
            8	Cr	209.5
            9	Cr	244.5
            0	Mn	7.434038
            1	Mn	15.63999
            2	Mn	33.668
            3	Mn	51.21
            4	Mn	72.41
            5	Mn	95.604
            6	Mn	119.203
            7	Mn	195.5
            8	Mn	221.89
            9	Mn	248.6
            0	Fe	7.9024681
            1	Fe	16.19921
            2	Fe	30.651
            3	Fe	54.91
            4	Fe	75
            5	Fe	98.985
            6	Fe	124.976
            7	Fe	151.06
            8	Fe	233.6
            9	Fe	262.1
            0	Co	7.88101
            1	Co	17.0844
            2	Co	33.5
            3	Co	51.27
            4	Co	79.5
            5	Co	102
            6	Co	128.9
            7	Co	157.8
            8	Co	186.14
            9	Co	275.4
            0	Ni	7.639878
            1	Ni	18.168838
            2	Ni	35.187
            3	Ni	54.92
            4	Ni	76.06
            5	Ni	108
            6	Ni	132
            7	Ni	162
            8	Ni	193.2
            9	Ni	224.7
            0	Cu	7.72638
            1	Cu	20.29239
            2	Cu	36.841
            3	Cu	57.38
            4	Cu	79.8
            5	Cu	103
            6	Cu	139
            7	Cu	166
            8	Cu	198
            9	Cu	232.2
            0	Zn	9.394197
            1	Zn	17.96439
            2	Zn	39.7233
            3	Zn	59.573
            4	Zn	82.6
            5	Zn	108
            6	Zn	133.9
            7	Zn	173.9
            8	Zn	203
            9	Zn	238
            0	Ga	5.999302
            1	Ga	20.51514
            2	Ga	30.72576
            3	Ga	63.241
            4	Ga	86.01
            5	Ga	112.7
            6	Ga	140.8
            7	Ga	169.9
            8	Ga	211
            9	Ga	244
            0	Ge	7.899435
            1	Ge	15.93461
            2	Ge	34.0576
            3	Ge	45.7155
            4	Ge	90.5
            5	Ge	115.9
            6	Ge	144.9
            7	Ge	176.4
            8	Ge	212.5
            9	Ge	252.1
            0	As	9.78855
            1	As	18.5892
            2	As	28.349
            3	As	50.15
            4	As	62.77
            5	As	121.19
            6	As	147
            7	As	180
            8	As	213
            9	As	247
            0	Se	9.752392
            1	Se	21.196
            2	Se	31.697
            3	Se	42.947
            4	Se	68.3
            5	Se	81.83
            6	Se	155.327
            7	Se	184
            8	Se	219
            9	Se	255
            0	Br	11.81381
            1	Br	21.591
            2	Br	34.871
            3	Br	47.782
            4	Br	59.595
            5	Br	87.39
            6	Br	103.03
            7	Br	192.61
            8	Br	224
            9	Br	261
            0	Kr	13.9996055
            1	Kr	24.35984
            2	Kr	35.838
            3	Kr	50.85
            4	Kr	64.69
            5	Kr	78.49
            6	Kr	109.13
            7	Kr	125.802
            8	Kr	233
            9	Kr	268
            0	Rb	4.1771281
            1	Rb	27.28954
            2	Rb	39.247
            3	Rb	52.2
            4	Rb	68.44
            5	Rb	82.9
            6	Rb	98.67
            7	Rb	132.79
            8	Rb	150.628
            9	Rb	277.12
            0	Sr	5.69486745
            1	Sr	11.0302765
            2	Sr	42.88353
            3	Sr	56.28
            4	Sr	70.7
            5	Sr	88
            6	Sr	104
            7	Sr	121.21
            8	Sr	158.33
            9	Sr	177.3
            0	Y	6.21726
            1	Y	12.2236
            2	Y	20.52441
            3	Y	60.6072
            4	Y	75.35
            5	Y	91.39
            6	Y	110.02
            7	Y	128.12
            8	Y	145.64
            9	Y	185.7
            0	Zr	6.634126
            1	Zr	13.13
            2	Zr	23.17
            3	Zr	34.41836
            4	Zr	80.348
            5	Zr	96.38
            6	Zr	112
            7	Zr	133.7
            8	Zr	153
            9	Zr	172.02
            0	Nb	6.75885
            1	Nb	14.32
            2	Nb	25.04
            3	Nb	37.611
            4	Nb	50.5728
            5	Nb	102.069
            6	Nb	119.1
            7	Nb	136
            8	Nb	159.2
            9	Nb	180
            0	Mo	7.09243
            1	Mo	16.16
            2	Mo	27.13
            3	Mo	40.33
            4	Mo	54.417
            5	Mo	68.82704
            6	Mo	125.638
            7	Mo	143.6
            8	Mo	164.12
            9	Mo	186.3
            0	Tc	7.11938
            1	Tc	15.26
            2	Tc	29.55
            3	Tc	41
            4	Tc	57
            5	Tc	72
            6	Tc	88
            7	Tc	150
            8	Tc	169
            9	Tc	189.9
            0	Ru	7.3605
            1	Ru	16.76
            2	Ru	28.47
            3	Ru	45
            4	Ru	59
            5	Ru	76
            6	Ru	93
            7	Ru	110
            8	Ru	178.41
            9	Ru	198
            0	Rh	7.4589
            1	Rh	18.08
            2	Rh	31.06
            3	Rh	42
            4	Rh	63
            5	Rh	80
            6	Rh	97
            7	Rh	115.1
            8	Rh	135
            9	Rh	207.51
            0	Pd	8.336839
            1	Pd	19.43
            2	Pd	32.93
            3	Pd	46
            4	Pd	61
            5	Pd	84.1
            6	Pd	101
            7	Pd	120
            8	Pd	141
            9	Pd	159.9
            0	Ag	7.576234
            1	Ag	21.4844
            2	Ag	34.8
            3	Ag	49
            4	Ag	65
            5	Ag	82
            6	Ag	106
            7	Ag	125
            8	Ag	145.1
            9	Ag	167
            0	Cd	8.99382
            1	Cd	16.908313
            2	Cd	37.468
            3	Cd	51
            4	Cd	67.9
            5	Cd	87
            6	Cd	105
            7	Cd	130.1
            8	Cd	150
            9	Cd	173
            0	In	5.7863558
            1	In	18.87041
            2	In	28.04415
            3	In	55.45
            4	In	69.3
            5	In	90
            6	In	109
            7	In	130.1
            8	In	156
            9	In	178
            0	Sn	7.343918
            1	Sn	14.63307
            2	Sn	30.506
            3	Sn	40.74
            4	Sn	77.03
            5	Sn	94
            6	Sn	112.9
            7	Sn	135
            8	Sn	156
            9	Sn	184
            0	Sb	8.608389
            1	Sb	16.626
            2	Sb	25.3235
            3	Sb	43.804
            4	Sb	55
            5	Sb	99.51
            6	Sb	117
            7	Sb	139
            8	Sb	162
            9	Sb	185
            0	Tb	9.009808
            1	Tb	18.6
            2	Tb	27.84
            3	Tb	37.4155
            4	Tb	59.3
            5	Tb	69.1
            6	Tb	124.2
            7	Tb	143
            8	Tb	167
            9	Tb	191.1
            0	I	10.45126
            1	I	19.13126
            2	I	29.57
            3	I	40.357
            4	I	51.52
            5	I	74.4
            6	I	87.61
            7	I	150.81
            8	I	171
            9	I	197
            0	Xe	12.1298437
            1	Xe	20.975
            2	Xe	31.05
            3	Xe	42.2
            4	Xe	54.1
            5	Xe	66.703
            6	Xe	91.6
            7	Xe	105.9778
            8	Xe	179.84
            9	Xe	202
            0	Cs	3.893905727
            1	Cs	23.15745
            2	Cs	33.195
            3	Cs	43
            4	Cs	56
            5	Cs	69.1
            6	Cs	82.9
            7	Cs	110.1
            8	Cs	125.61
            9	Cs	213.3
            0	Ba	5.2116646
            1	Ba	10.003826
            2	Ba	35.8438
            3	Ba	47
            4	Ba	58
            5	Ba	71
            6	Ba	86
            7	Ba	101
            8	Ba	130.5
            9	Ba	146.52
            0	La	5.5769
            1	La	11.18496
            2	La	19.1773
            3	La	49.95
            4	La	61.6
            5	La	74
            6	La	88
            7	La	105
            8	La	119
            9	La	151.4
            0	Ce	5.5386
            1	Ce	10.956
            2	Ce	20.1974
            3	Ce	36.906
            4	Ce	65.55
            5	Ce	77.6
            6	Ce	91
            7	Ce	106
            8	Ce	125
            9	Ce	140
            0	Pr	5.4702
            1	Pr	10.631
            2	Pr	21.6237
            3	Pr	38.981
            4	Pr	57.53
            5	Pr	82
            6	Pr	97
            7	Pr	112
            8	Pr	131
            9	Pr	148
            0	Nd	5.525
            1	Nd	10.783
            2	Nd	22.09
            3	Nd	40.6
            4	Nd	60
            5	Nd	84
            6	Nd	99
            7	Nd	114
            8	Nd	136
            9	Nd	152
            0	Pm	5.58187
            1	Pm	10.938
            2	Pm	22.44
            3	Pm	41.17
            4	Pm	61.7
            5	Pm	85
            6	Pm	101
            7	Pm	116
            8	Pm	138
            9	Pm	155
            0	Sm	5.643722
            1	Sm	11.078
            2	Sm	23.55
            3	Sm	41.64
            4	Sm	62.7
            5	Sm	87
            6	Sm	103
            7	Sm	118
            8	Sm	141
            9	Sm	158
            0	Eu	5.670385
            1	Eu	11.24
            2	Eu	24.84
            3	Eu	42.94
            4	Eu	63.2
            5	Eu	89
            6	Eu	105
            7	Eu	120
            8	Eu	144
            9	Eu	161
            0	Gd	6.1498
            1	Gd	12.076
            2	Gd	20.54
            3	Gd	44.44
            4	Gd	64.8
            5	Gd	89
            6	Gd	106
            7	Gd	123
            8	Gd	144
            9	Gd	165
            0	Tb	5.8638
            1	Tb	11.513
            2	Tb	21.82
            3	Tb	39.33
            4	Tb	66.5
            5	Tb	90
            6	Tb	108
            7	Tb	125
            8	Tb	143
            9	Tb	168
            0	Dy	5.93905
            1	Dy	11.647
            2	Dy	22.89
            3	Dy	41.23
            4	Dy	62.1
            5	Dy	93
            6	Dy	110
            7	Dy	127
            8	Dy	152
            9	Dy	170
            0	Ho	6.0215
            1	Ho	11.781
            2	Ho	22.79
            3	Ho	42.52
            4	Ho	63.9
            5	Ho	95
            6	Ho	112
            7	Ho	129
            8	Ho	155
            9	Ho	173
            0	Er	6.1077
            1	Er	11.916
            2	Er	22.7
            3	Er	42.42
            4	Er	65.1
            5	Er	96
            6	Er	114
            7	Er	131
            8	Er	158
            9	Er	177
            0	Tm	6.18431
            1	Tm	12.065
            2	Tm	23.66
            3	Tm	42.41
            4	Tm	65.4
            5	Tm	98
            6	Tm	116
            7	Tm	133
            8	Tm	160
            9	Tm	180
            0	Yb	6.25416
            1	Yb	12.179185
            2	Yb	25.053
            3	Yb	43.61
            4	Yb	65.6
            5	Yb	99
            6	Yb	117
            7	Yb	135
            8	Yb	163
            9	Yb	182
            0	Lu	5.425871
            1	Lu	14.13
            2	Lu	20.9594
            3	Lu	45.249
            4	Lu	66.8
            5	Lu	98
            6	Lu	117
            7	Lu	136
            8	Lu	159
            9	Lu	185
            0	Hf	6.82507
            1	Hf	14.61
            2	Hf	22.55
            3	Hf	33.37
            4	Hf	68.37
            5	Hf	98
            6	Hf	118
            7	Hf	137
            8	Hf	157
            9	Hf	187
            0	Ta	7.549571
            1	Ta	16.2
            2	Ta	23.1
            3	Ta	35
            4	Ta	48.272
            5	Ta	94.01
            6	Ta	119
            7	Ta	139
            8	Ta	159
            9	Ta	180
            0	W	7.86403
            1	W	16.37
            2	W	26
            3	W	38.2
            4	W	51.6
            5	W	64.77
            6	W	122.01
            7	W	141.2
            8	W	160.2
            9	W	179
            0	Re	7.83352
            1	Re	16.6
            2	Re	27
            3	Re	39.1
            4	Re	51.9
            5	Re	67
            6	Re	82.71
            7	Re	144.4
            8	Re	165
            9	Re	187
            0	Os	8.43823
            1	Os	17
            2	Os	25
            3	Os	41
            4	Os	55
            5	Os	70.1
            6	Os	85.1
            7	Os	102.02
            8	Os	168.7
            9	Os	190
            0	Ir	8.96702
            1	Ir	17
            2	Ir	28
            3	Ir	40
            4	Ir	57
            5	Ir	72
            6	Ir	89
            7	Ir	105
            8	Ir	122.7
            9	Ir	194.8
            0	Pt	8.95883
            1	Pt	18.56
            2	Pt	29
            3	Pt	43
            4	Pt	56
            5	Pt	75
            6	Pt	91
            7	Pt	109
            8	Pt	126
            9	Pt	144.9
            0	Au	9.225554
            1	Au	20.203
            2	Au	30
            3	Au	45
            4	Au	60
            5	Au	74
            6	Au	94
            7	Au	112
            8	Au	130.1
            9	Au	149
            0	Hg	10.437504
            1	Hg	18.75688
            2	Hg	34.49
            3	Hg	48.55
            4	Hg	61.2
            5	Hg	76.6
            6	Hg	93
            7	Hg	113.9
            8	Hg	134
            9	Hg	153
            0	Tl	6.1082873
            1	Tl	20.4283
            2	Tl	29.852
            3	Tl	51.14
            4	Tl	62.6
            5	Tl	80
            6	Tl	97.9
            7	Tl	116
            8	Tl	135
            9	Tl	158
            0	Pb	7.4166799
            1	Pb	15.032499
            2	Pb	31.9373
            3	Pb	42.33256
            4	Pb	68.8
            5	Pb	82.9
            6	Pb	100.1
            7	Pb	120
            8	Pb	138
            9	Pb	158
            0	Bi	7.285516
            1	Bi	16.703
            2	Bi	25.57075
            3	Bi	45.37
            4	Bi	54.856
            5	Bi	88.4
            6	Bi	103
            7	Bi	122
            8	Bi	143
            9	Bi	161.1
            0	Po	8.41807
            1	Po	19.3
            2	Po	27.3
            3	Po	36
            4	Po	57
            5	Po	69.1
            6	Po	108
            7	Po	125
            8	Po	146.1
            9	Po	166
            0	At	9.31751
            1	At	17.88
            2	At	26.58
            3	At	39.65
            4	At	50.39
            5	At	72
            6	At	85.1
            7	At	130.1
            8	At	149
            9	At	169
            0	Rn	10.7485
            1	Rn	21.4
            2	Rn	29.4
            3	Rn	36.9
            4	Rn	52.9
            5	Rn	64
            6	Rn	88
            7	Rn	102
            8	Rn	154
            9	Rn	173.9
            0	Fr	4.0727411
            1	Fr	22.4
            2	Fr	33.5
            3	Fr	39.1
            4	Fr	50
            5	Fr	67
            6	Fr	80
            7	Fr	106
            8	Fr	120
            9	Fr	179
            0	Ra	5.2784239
            1	Ra	10.14718
            2	Ra	31
            3	Ra	41
            4	Ra	52.9
            5	Ra	64
            6	Ra	82
            7	Ra	97
            8	Ra	124
            9	Ra	140
            0	Ac	5.380226
            1	Ac	11.75
            2	Ac	17.431
            3	Ac	44.8
            4	Ac	55
            5	Ac	67
            6	Ac	79
            7	Ac	98.9
            8	Ac	113.9
            9	Ac	143.9
            0	Th	6.3067
            1	Th	12.1
            2	Th	18.32
            3	Th	28.648
            4	Th	58
            5	Th	69.1
            6	Th	82
            7	Th	95
            8	Th	118
            9	Th	133
            0	Pa	5.89
            1	Pa	11.9
            2	Pa	18.6
            3	Pa	30.9
            4	Pa	44.3
            5	Pa	72
            6	Pa	85.1
            7	Pa	98.9
            8	Pa	111
            9	Pa	137
            0	U	6.19405
            1	U	11.6
            2	U	19.8
            3	U	36.7
            4	U	46
            5	U	62
            6	U	89
            7	U	101
            8	U	116
            9	U	128.9
            0	Np	6.26554
            1	Np	11.5
            2	Np	19.7
            3	Np	33.8
            4	Np	48
            5	Np	65
            6	Np	92
            7	Np	107
            8	Np	121
            9	Np	136
            0	Pu	6.02576
            1	Pu	11.5
            2	Pu	21.1
            3	Pu	35
            4	Pu	49
            5	Pu	80
            6	Pu	95
            7	Pu	109
            8	Pu	124
            9	Pu	139
            0	Am	5.97381
            1	Am	11.7
            2	Am	21.7
            3	Am	36.8
            4	Am	50
            5	Am	67.9
            6	Am	95
            7	Am	110
            8	Am	125
            9	Am	141
            0	Cm	5.99141
            1	Cm	12.4
            2	Cm	20.1
            3	Cm	37.7
            4	Cm	51
            5	Cm	69.1
            6	Cm	97
            7	Cm	112
            8	Cm	128
            9	Cm	144
            0	Bk	6.19785
            1	Bk	11.9
            2	Bk	21.6
            3	Bk	36
            4	Bk	56
            5	Bk	70.1
            6	Bk	90
            7	Bk	114
            8	Bk	130
            9	Bk	147
            0	Cf	6.28166
            1	Cf	12
            2	Cf	22.4
            3	Cf	37.7
            4	Cf	51.9
            5	Cf	75
            6	Cf	91
            7	Cf	112.9
            8	Cf	133
            9	Cf	152
            0	Es	6.36758
            1	Es	12.2
            2	Es	22.7
            3	Es	38.8
            4	Es	54.1
            5	Es	71
            6	Es	97
            7	Es	112.9
            8	Es	137
            9	Es	157
            0	Fm	6.5
            1	Fm	12.4
            2	Fm	23.2
            3	Fm	39.3
            4	Fm	55
            5	Fm	74
            6	Fm	93
            7	Fm	120
            8	Fm	136
            9	Fm	162
            0	Md	6.58
            1	Md	12.4
            2	Md	24.3
            3	Md	40
            4	Md	54.1
            5	Md	76
            6	Md	96
            7	Md	115.1
            8	Md	143.9
            9	Md	162
            0	No	6.62621
            1	No	12.93
            2	No	25.8
            3	No	41.5
            4	No	60
            5	No	74
            6	No	97
            7	No	119
            8	No	140
            9	No	170
            0	Lr	4.96
            1	Lr	14.54
            2	Lr	21.8
            3	Lr	43.6
            4	Lr	56
            5	Lr	80
            6	Lr	96
            7	Lr	121
            8	Lr	143
            9	Lr	165
            0	Rf	6.02
            1	Rf	14.35
            2	Rf	23.84
            3	Rf	31.87
            4	Rf	64
            5	Rf	77
            6	Rf	102
            7	Rf	119
            8	Rf	146.1
            9	Rf	169
            0	Db	6.8
            1	Db	14
            2	Db	23.1
            3	Db	33
            4	Db	43
            5	Db	86
            6	Db	98.9
            7	Db	126
            8	Db	145.1
            9	Db	172
            0	Sg	7.8
            1	Sg	17.1
            2	Sg	25.8
            3	Sg	35.5
            4	Sg	47.2
            5	Sg	59.3
            6	Sg	109
            7	Sg	122
            8	Sg	152
            9	Sg	170
            0	Bh	7.7
            1	Bh	17.5
            2	Bh	26.7
            3	Bh	37.3
            4	Bh	49
            5	Bh	62.1
            6	Bh	74.9
            7	Bh	134
            8	Bh	148
            9	Bh	178
            0	Hs	7.6
            1	Hs	18.2
            2	Hs	29.3
            3	Hs	37.7
            4	Hs	51.2
            5	Hs	64
            6	Hs	78.1
            7	Hs	91.7
            8	Hs	159.9
            9	Hs	173.9
            4	Mt	50
            5	Mt	
            6	Mt	
            7	Mt	94
            8	Mt	109
            9	Mt	187
            5	Ds	65
            6	Ds	
            7	Ds	
            8	Ds	112.9
            9	Ds	128
            """.strip()

    KIEs: dict[str, list[float]] = {}
    for line in data_str.splitlines():
        parts = line.split()
        if len(parts) < 3:
            continue     # skip blank or malformed lines
        _, elem, energy = parts
        KIEs.setdefault(elem, []).append(float(energy))

    return KIEs


def HighestKnownONs() -> dict[str, int]:
    """
    Determines the highest known oxidation states for each metal element.

    Returns:
        HKONs (dict[str, int]) : dictionary with metal element symbols as keys
                              and their highest known oxidation state as values.
    """
    
    data_str = """
        H	-1	1										
            He												
            Li	1											
            Be	0	1	2									
            B	-5	-1	0	1	2	3						
            C	-4	-3	-2	-1	0	1	2	3	4			
            N	-3	-2	-1	1	2	3	4	5				
            O	-2	-1	0	1	2							
            F	-1	0										
            Ne												
            Na	-1	1										
            Mg	0	1	2									
            Al	-2	-1	1	2	3							
            Si	-4	-3	-2	-1	0	1	2	3	4			
            P	-3	-2	-1	0	1	2	3	4	5			
            S	-2	-1	0	1	2	3	4	5	6			
            Cl	-1	1	2	3	4	5	6	7				
            Ar	0											
            K	-1	1										
            Ca	0	1	2									
            Sc	0	1	2	3								
            Ti	-2	-1	0	1	2	3	4					
            V	-3	1	0	1	2	3	4	5				
            Cr	-4	-2	-1	0	1	2	3	4	5	6		
            Mn	-3	-2	-1	0	1	2	3	4	5	6	7	
            Fe	-4	-2	-1	0	1	2	3	4	5	6	7	
            Co	-3	-1	0	1	2	3	4	5				
            Ni	-2	-1	0	1	2	3	4					
            Cu	-2	0	1	2	3	4						
            Zn	-2	0	1	2								
            Ga	-5	-4	-3	-2	-1	1	2	3				
            Ge	-3	-2	-1	0	1	2	3	4				
            As	-3	-2	-1	0	1	2	3	4	5			
            Se	-2	-1	1	2	3	4	5	6				
            Br	-1	1	3	4	5	7						
            Kr	0	1	2									
            Rb	-1	1										
            Sr	0	1	2									
            Y	0	1	2	3								
            Zr	-2	0	1	2	3	4						
            Nb	-3	-1	0	1	2	3	4	5				
            Mo	-4	-2	-1	0	1	2	3	4	5	6		
            Tc	-3	-1	0	1	2	3	4	5	6	7		
            Ru	-4	-2	0	1	2	3	4	5	6	7	8	
            Rh	-3	-1	0	1	2	3	4	5	6			
            Pd	0	1	2	3	4							
            Ag	-2	-1	1	2	3							
            Cd	-2	1	2									
            In	-5	-2	-2	1	2	3						
            Sn	-4	-3	-2	-1	0	1	2	3	4			
            Sb	-3	-2	-1	0	1	2	3	4	5			
            Te	-2	-1	1	3	4	5	6	7				
            I	-1	1	3	4	5	6	7					
            Xe	0	2	4	6	8							
            Cs	-1	1										
            Ba	0	1	2									
            La	0	1	2	3								
            Ce	0	2	3	4								
            Pr	0	1	2	3	4	5						
            Nd	0	2	3	4								
            Pm	2	3										
            Sm	0	2	3									
            Eu	0	2	3									
            Gd	0	1	2	3								
            Tb	0	1	2	3	4							
            Dy	0	2	3	4								
            Ho	0	2	3									
            Er	0	2	3									
            Tm	0	2	3									
            Yb	0	2	3									
            Lu	0	2	3									
            Hf	-2	0	1	2	3	4						
            Ta	-3	-1	0	1	2	3	4	5				
            W	-4	-2	-1	0	1	2	3	4	5	6		
            Re	-3	-1	0	1	2	3	4	5	6	7		
            Os	-4	-2	-1	0	1	2	3	4	5	6	7	8
            Ir	-3	-1	0	1	2	3	4	5	6	7	8	9
            Pt	-3	-2	-1	0	1	2	3	4	5	6		
            Au	-3	-2	-1	0	1	2	3	5				
            Hg	-2	1	2									
            Tl	-5	-2	-1	1	2	3						
            Pb	-4	-2	-1	1	2	3	4					
            Bi	-3	-2	-1	1	2	3	4	5				
            Po	-2	2	4	5	6							
            At	-1	1	3	5	7							
            Rn	2	6										
            Fr	1											
            Ra	2											
            Ac	2	3										
            Th	1	2	3	4								
            Pa	3	4	5									
            U	1	2	3	4	5	6						
            Np	2	3	4	5	6	7						
            Pu	2	3	4	5	6	7	8					
            Am	2	3	4	5	6	7						
            Cm	3	4	5	6								
            Bk	2	3	4	5								
            Cf	2	3	4	5								
            Es	2	3	4									
            Fm	2	3										
            Md	2	3										
            No	2	3										
            Lr	3											
            Rf	4											
            Db	5											
            Sg	0	6										
            Bh	7											
            Hs	8											
            Mt												
            Ds												
            Rg												
            Cn	2											
            Nh												
            Fl												
            Mc												
            Lv												
            Ts												
            Og						
        """.strip()

    HKONs: dict[str, int] = {}
    for line in data_str.splitlines():
        parts = line.split()
        elem = parts[0]
        # ignore entries with no numeric states
        if len(parts) > 1:
            # convert all remaining tokens to ints, take the max
            states = [int(x) for x in parts[1:]]
            HKONs[elem] = max(states)
        else:
            # if no oxidation state listed, set to 0 or remove entirely
            HKONs[elem] = 0

    return HKONs


def ONprobabilities() -> dict[str, list[float]]:
    """
    Reads in the probability of each oxidation state for all metal elements.
    Approximate probabilities are assessed by their relative frequency of
    occurence in the CSD metadata.

    Returns:
        ONP (dict[str, list[float]]): dictionary with metal element symbols as keys
                        and a list of the probability at the relevant oxidation
                        states as values.
    """
    data_str = """
                    Li	0	0.993162393	0.003418803	0.001709402	0.001709402	0	0	0	0	0	0
                    Be	0.023529412	0	0.976470588	0	0	0	0	0	0	0	0
                    Na	0	0.998743719	0	0.001256281	0	0	0	0	0	0	0
                    Mg	0.002141328	0.016416845	0.981441827	0	0	0	0	0	0	0	0
                    Al	0	0.020121951	0.016463415	0.962804878	0.000609756	0	0	0	0	0	0
                    K	0	0.998137803	0	0.001862197	0	0	0	0	0	0	0
                    Ca	0	0.001176471	0.998823529	0	0	0	0	0	0	0	0
                    Sc	0	0.003656307	0.010968921	0.985374771	0	0	0	0	0	0	0
                    Ti	0.004095004	0.00020475	0.028460278	0.128992629	0.838247338	0	0	0	0	0	0
                    V	0.002768166	0.007612457	0.037197232	0.133564014	0.391868512	0.426816609	0.00017301	0	0	0	0
                    Cr	0.160139589	0.017836371	0.138037999	0.612252811	0.014928267	0.018805739	0.037999225	0	0	0	0
                    Mn	0.001919631	0.044599437	0.64864346	0.260621961	0.039928334	0.003519324	0.000319939	0.000447914	0	0	0
                    Fe	0.016919249	0.011083465	0.545216014	0.42053834	0.0059715	0.000226193	4.52E-05	0	0	0	0
                    Co	0.002508241	0.021678372	0.655403469	0.31990827	0.00035832	0.000143328	0	0	0	0	0
                    Ni	0.020168417	0.010846588	0.935890772	0.031222927	0.001871296	0	0	0	0	0	0
                    Cu	6.78E-05	0.236206984	0.75977341	0.003951765	0	0	0	0	0	0	0
                    Zn	7.78E-05	0.000778362	0.999104884	3.89E-05	0	0	0	0	0	0	0
                    Ga	0.003392706	0.100932994	0.066157761	0.828668363	0.000848176	0	0	0	0	0	0
                    Rb	0	0.9875	0.0125	0	0	0	0	0	0	0	0
                    Sr	0	0	1	0	0	0	0	0	0	0	0
                    Y	0.000608273	0.000608273	0.001824818	0.996350365	0.000608273	0	0	0	0	0	0
                    Zr	0.001476378	0	0.030019685	0.039370079	0.928149606	0	0.000984252	0	0	0	0
                    Nb	0.006417112	0.03315508	0.02459893	0.100534759	0.19144385	0.643850267	0	0	0	0	0
                    Mo	0.067368705	0.007830431	0.106115836	0.052922911	0.132847307	0.184150128	0.448764682	0	0	0	0
                    Tc	0	0.1488	0.0592	0.2256	0.056	0.4288	0.04	0.0416	0	0	0
                    Ru	0.01519372	0.011901747	0.819870009	0.119355111	0.027348696	0.000759686	0.005571031	0	0	0	0
                    Rh	0.002670673	0.455850442	0.141545652	0.397095643	0.00166917	0.001168419	0	0	0	0	0
                    Pd	0.018729923	0.009834445	0.961255251	0.004101804	0.005979738	9.88E-05	0	0	0	0	0
                    Ag	0.000163733	0.988784282	0.006631191	0.004420794	0	0	0	0	0	0	0
                    Cd	0	0.000421977	0.998818466	0.000759558	0	0	0	0	0	0	0
                    In	0	0.056062581	0.024771838	0.91916558	0	0	0	0	0	0	0
                    Sn	0.00222187	0.00095223	0.194730995	0.001110935	0.800666561	0.000158705	0.000158705	0	0	0	0
                    Cs	0	0.986486486	0.013513514	0	0	0	0	0	0	0	0
                    Ba	0	0	0.998236332	0.001763668	0	0	0	0	0	0	0
                    Hf	0	0	0.015473888	0.025145068	0.959381044	0	0	0	0	0	0
                    Ta	0.003726708	0.016149068	0.01242236	0.068322981	0.116770186	0.782608696	0	0	0	0	0
                    W	0.128280743	0.003538779	0.131229726	0.027720436	0.146564435	0.130934827	0.431731053	0	0	0	0
                    Re	0.001812324	0.316149819	0.044905356	0.124043496	0.065243657	0.364679823	0.017921869	0.065243657	0	0	0
                    Os	0.009817672	0.005610098	0.460729313	0.140953717	0.204067321	0.016830295	0.152173913	0.000701262	0.00911641	0	0
                    Ir	0.001902045	0.201854494	0.027104137	0.752258678	0.011650024	0.005230623	0	0	0	0	0
                    Pt	0.017449918	0.005606652	0.844021671	0.01870984	0.114085927	0	0.000125992	0	0	0	0
                    Au	0.001029204	0.731892448	0.015824006	0.251125691	0.00012865	0	0	0	0	0	0
                    Hg	0.00251004	0.013805221	0.983433735	0.000251004	0	0	0	0	0	0	0
                    Tl	0.001342282	0.616107383	0.013422819	0.369127517	0	0	0	0	0	0	0
                    Pb	0	0.000545852	0.964792576	0.000272926	0.034388646	0	0	0	0	0	0
                    Bi	0	0	0	1	0	0	0	0	0	0	0
                    Po	0	0	0	0	0	0	0	0	0	0	0
                    Fr	0	0	0	0	0	0	0	0	0	0	0
                    Ra	0	0	0	0	0	0	0	0	0	0	0
                    Rf	0	0	0	0	0	0	0	0	0	0	0
                    Db	0	0	0	0	0	0	0	0	0	0	0
                    Sg	0	0	0	0	0	0	0	0	0	0	0
                    Bh	0	0	0	0	0	0	0	0	0	0	0
                    Hs	0	0	0	0	0	0	0	0	0	0	0
                    Mt	0	0	0	0	0	0	0	0	0	0	0
                    Ds	0	0	0	0	0	0	0	0	0	0	0
                    Rg	0	0	0	0	0	0	0	0	0	0	0
                    Cn	0	0	0	0	0	0	0	0	0	0	0
                    Nh	0	0	0	0	0	0	0	0	0	0	0
                    Fl	0	0	0	0	0	0	0	0	0	0	0
                    Mc	0	0	0	0	0	0	0	0	0	0	0
                    Lv	0	0	0	0	0	0	0	0	0	0	0
                    La	0	0	0.002631579	0.996842105	0.000526316	0	0	0	0	0	0
                    Ce	0	0.000666223	0.005996003	0.742838108	0.250499667	0	0	0	0	0	0
                    Pr	0	0	0.002636204	0.993848858	0.003514938	0	0	0	0	0	0
                    Nd	0	0	0.005083884	0.994916116	0	0	0	0	0	0	0
                    Pm	0	0	0	1	0	0	0	0	0	0	0
                    Sm	0	0	0.153056235	0.846943765	0	0	0	0	0	0	0
                    Eu	0	0	0.081716037	0.918283963	0	0	0	0	0	0	0
                    Gd	0.000430663	0	0.005167959	0.994401378	0	0	0	0	0	0	0
                    Tb	0	0	0.132275132	0.866717057	0.001007811	0	0	0	0	0	0
                    Dy	0	0	0.007676768	0.992323232	0	0	0	0	0	0	0
                    Ho	0.001218027	0	0.007308161	0.990255786	0.001218027	0	0	0	0	0	0
                    Er	0	0	0.099644793	0.899233502	0.001121705	0	0	0	0	0	0
                    Tm	0	0	0.112903226	0.887096774	0	0	0	0	0	0	0
                    Yb	0	0	0.268343816	0.731656184	0	0	0	0	0	0	0
                    Lu	0	0	0.003663004	0.996336996	0	0	0	0	0	0	0
                    Ac	0	0	0	0	0.333333333	0.666666667	0	0	0	0	0
                    Th	0	0	0.008896797	0.019572954	0.96975089	0	0.001779359	0	0	0	0
                    Pa	0	0	0	0	0.333333333	0.666666667	0	0	0	0	0
                    U	0.000308452	0	0.003701419	0.098704503	0.333127699	0.064466379	0.499691548	0	0	0	0
                    Np	0	0	0.004048583	0.04048583	0.246963563	0.453441296	0.251012146	0.004048583	0	0	0
                    Pu	0	0	0	0.231292517	0.510204082	0.040816327	0.217687075	0	0	0	0
                    Am	0	0	0	0.857142857	0	0.057142857	0.085714286	0	0	0	0
                    Cm	0	0	0	1	0	0	0	0	0	0	0
                    Bk	0	0	0	1	0	0	0	0	0	0	0
                    Cf	0	0	0	1	0	0	0	0	0	0	0
                    Es	0	0	0	0	0	0	0	0	0	0	0
                    Fm	0	0	0	0	0	0	0	0	0	0	0
                    Md	0	0	0	0	0	0	0	0	0	0	0
                    No	0	0	0	0	0	0	0	0	0	0	0
                    Lr	0	0	0	0	0	0	0	0	0	0	0
                    Ge	0.00867052	0.012524085	0.512524085	0.002890173	0.462427746	0	0.000963391	0	0	0	0
                    As	0	0.021459227	0	0.811158798	0	0.167381974	0	0	0	0	0
                    Se	0	0.023809524	0.619047619	0	0.325396825	0	0.031746032	0	0	0	0
                    Sb	0	0.002742732	0	0.541963796	0.000548546	0.45419638	0.000548546	0	0	0	0
                    Te	0	0.001194743	0.338112306	0.001194743	0.628434886	0	0.031063321	0	0	0	0
                    """.strip()

    ONP: dict[str, list[float]] = {}
    for line in data_str.splitlines():
        parts = line.split()  # tab/space 都能分
        element = parts[0]
        probs = [float(x) for x in parts[1:]]
        ONP[element] = probs

    return ONP


def ONmostprob(iONP: dict[str, list[float]]) -> dict[str, int]:
    """
    Determines the highest probability oxidation state for each metal element.
    These values are utilized during charge distribution routines.

    Parameters:
        iONP (dict[str, list[float]]): dictionary with metal element symbols as
                        keys and a list of the probability at the relevant oxidation
                        states as values.

    Returns:
        MPOS (dict[str, int]) : dictionary with metal element symbols as keys and
                        their oxidation state with the highest probability as values.
    """
    MPOS = {}
    for metal in iONP:
        hprob = 0
        for index, prob in enumerate(iONP[metal]):
            if prob >= hprob:
                hprob = prob
                MPOS[metal] = index
    return MPOS


def getCN(lsites: dict[Atom, list[Atom]]) -> dict[Molecule, int]:
    """
    Determines the highest probability oxidation state for each metal element.
    These values are utilized during charge distribution routines.

    Parameters:
        lsites (dict[Atom, list[Atom]]): dictionary with metal Atom object as keys
                        and the list of ligand atoms which bind them as values.

    Returns:
        CNdict (dict[Molecule, int]): dictionary with as metal Atom objects as
                        keys and effective coordination number as values.
    """
    CNdict = {}
    for metal in lsites:
        CNdict[metal] = 0
        for ligand in lsites[metal]:
            if hapticity(ligand, metal):
                CNdict[metal] += 0.5
            else:
                CNdict[metal] += 1
        for neighbour in metal.neighbours:
            if neighbour.is_metal:
                CNdict[metal] += 0.5
    return CNdict


KnownON  = KnownONs()
KnownIE  = IonizationEnergies()
HighestKnownON = HighestKnownONs()
ONProb   = ONprobabilities()
HighestProbON = ONmostprob(ONProb)


def check(file_path: str) -> dict[str, list[str]]:
    """
    Process a single .cif or .mol2 file and return a dict mapping each metal site
    label to the list of flag-names that are not 'GOOD'.
    """
    # Read structure
    try:
        if file_path.endswith(".cif"):
            cif = readentry(file_path)
            mol = cif.molecule
            asymmol = cif.asymmetric_unit_molecule
        else:
            mol = readSBU(file_path)
            asymmol = mol
    except RuntimeError:
        return {}

    # Identify sites
    uniquesites = get_unique_sites(mol, asymmol)
    metalsites = get_metal_sites(uniquesites)
    if not metalsites:
        return {}

    # Bonding and oxidation contributions
    dVBO = delocalisedLBO(mol)
    rVBO = ringVBOs(mol)
    AON = iVBS_Oxidation_Contrib(uniquesites, rVBO, dVBO)
    rAON = redundantAON(AON, mol)
    ligand_sites = get_ligand_sites(metalsites, uniquesites)
    binding_sites = get_binding_sites(metalsites, uniquesites)
    binding_sphere = binding_domain(binding_sites, rAON, mol, {u.label: u for u in uniquesites})
    bindingAON = binding_contrib(binding_sphere, binding_sites, rAON)
    connected = get_metal_networks(ligand_sites, binding_sphere, bindingAON)
    mCN = getCN(ligand_sites)

    # ONE-center inner
    ONEC_inner = {}
    for metal in ligand_sites:
        oxidation = 0
        val = valence_e(metal)
        for lig in ligand_sites[metal]:
            LBO = bindingAON[lig]
            Nbr = bridging(lig)
            Ox = LBO / Nbr
            oxidation += Ox
            if Ox >= 2:
                mCN[metal] += 1
            if Ox >= 3:
                mCN[metal] += 1
        ONEC_inner[metal] = [oxidation, val + 2*mCN[metal] - oxidation]

    # Network distribution
    noint_balance = distribute_ONEC(ONEC_inner, connected, KnownIE, ONProb, HighestKnownON, mCN, HighestProbON)

    # Outer-sphere
    OSC = outer_sphere_contrib(outer_sphere_domain(uniquesites, binding_sphere), rAON)
    noint_outer = distribute_OuterSphere(noint_balance, OSC, KnownIE, ONProb, HighestKnownON, mCN)

    # Prepare flag dictionaries
    flags = {name: {} for name in [
        "impossible", "unknown", "zero_valence", "noint_flag",
        "low_prob_1", "low_prob_2", "low_prob_3", "low_prob_multi"
    ]}

    # Evaluate flags
    for metal in noint_outer:
        lbl = metal.label
        flags["impossible"][lbl] = "GOOD"
        flags["unknown"][lbl]    = "GOOD"
        flags["zero_valence"][lbl]       = "GOOD"
        flags["noint_flag"][lbl]         = "GOOD"
        flags["low_prob_1"][lbl]         = "GOOD"
        flags["low_prob_2"][lbl]         = "GOOD"
        flags["low_prob_3"][lbl]         = "GOOD"
        flags["low_prob_multi"][lbl]     = "GOOD"
        
        val = valence_e(metal)
        os_val = noint_outer[metal][0]

        if os_val > val:
            flags["impossible"][lbl] = "BAD"
        if os_val == 0:
            flags["zero_valence"][lbl] = "BAD"
        if not any(math.isclose(os_val, i, abs_tol=0.5) for i in KnownON[metal.atomic_symbol]):
            flags["unknown"][lbl] = "BAD"
        # low-prob checks
        if math.isclose(os_val % 1, 0, abs_tol=1e-4):
            prob = ONProb[metal.atomic_symbol][round(os_val)] if round(os_val) < len(ONProb[metal.atomic_symbol]) else 0
            if prob < 0.01:
                flags["low_prob_1"][lbl] = "LIKELY_BAD"
            if prob < 0.001:
                flags["low_prob_2"][lbl] = "BAD"
            if prob < 0.0001:
                flags["low_prob_3"][lbl] = "BAD"
        
    return flags

result = check(file_path = "UDOJOQ_ION_pacman.cif")
print(result)