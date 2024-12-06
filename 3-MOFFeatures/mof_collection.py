import os
import glob
import json
import time
import pickle
import shutil
import random
import warnings
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from multiprocessing import Process, cpu_count, Array
from sys import exit
pd.options.display.max_rows = 1000
import copy
from pymatgen.core import Structure
import sys
import itertools
import hashlib
import datetime
import math
from ase.io import read
from pymatgen.io.cif import CifParser
from pymatgen.io.ase import AseAtomsAdaptor


class Atom:
    """A class to hold atomic information, and check bonds."""

    def __init__(self, element):
        """Create an Atom object given an element.

        :param element: Element of the Atom object. This determines it's
        properties.
        """
        # Covalent radii taken from DOI: Covalent radii revisited
        # Beatriz Cordero,   Verónica Gómez,   Ana E. Platero-Prats,
        # Marc Revés,   Jorge Echeverría,  Eduard Cremades, Flavia Barragána
        # and Santiago Alvarez Dalton Trans., 2008, 2832-2838
        # DOI: 10.1039/B801115J
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
        """Get the maximum possible distance between the Atom object and
        another atom of type ele2.

        Returns the some of covelant radii for Atom and Atom(ele2) plus their
        covalant radii.

        :param bond_tol: Bold tolerance to use, if not set then the default will
        be used.
        :param ele2: Element of the atom that the max_bond corresponds to.
        :return:
        """
        if bond_tol is None:
            bond_tol = self.bond_tolerance(ele2)
        return Atom(ele2).co + self.co + bond_tol


class MofStructure(Structure):
    """Extend the pymatgen Structure class to add MOF specific features"""

    def __init__(self, lattice, species, coords, charge=None,
                 validate_proximity=False, to_unit_cell=False,
                 coords_are_cartesian=False, site_properties=None, name="N/A"):

        """Create a MOf structure. The arguments are the same as in the
        pymatgen Structure class with the addition of the name argument.
        The super constructor is called and additional MOF specific properties
        are initialized.

        :param name: MOF name, used to identify the structure.
        """
        super().__init__(lattice, species, coords,
                         charge=charge,
                         validate_proximity=validate_proximity,
                         to_unit_cell=to_unit_cell,
                         coords_are_cartesian=coords_are_cartesian,
                         site_properties=site_properties)

        self._all_coord_spheres_indices = None
        self._all_distances = None
        self._metal_coord_spheres = []
        self._name = name
        self.metal = None
        self.metal_indices = []
        self.organic = None
        self.species_str = [str(s) for s in self.species]

        metal_set = set([s for s in self.species_str if Atom(s).is_metal])
        non_metal_set = set([s for s in self.species_str
                             if not Atom(s).is_metal])
        todays_date = datetime.datetime.now().isoformat()
        self.summary = {'cif_okay': 'N/A',
                        'problematic': 'N/A',
                        'has_oms': 'N/A',
                        'metal_sites': [],
                        'oms_density': 'N/A',
                        'checksum': 'N/A',
                        'metal_species': list(metal_set),
                        'non_metal_species': list(non_metal_set),
                        'name': name,
                        'uc_volume': self.volume,
                        'density': self.density,
                        'date_created': str(todays_date)}
        
        self._tolerance = None
        self._split_structure_to_organic_and_metal()

    @classmethod
    def from_file(cls, filename, primitive=False, sort=False, merge_tol=0.0):
        """Create a MofStructure from a CIF file.

        This makes use of the from_file function of the Structure class and
        catches the exception in case a CIF file cannot be read.
        If the CIF is read successfully then the MofStructure is marked as okay,
        and the file checksum is added to the summary. If the CIF file cannot be
        read then it is marked as not okay and all the other properties are
        set to None and because there cannot be an empty Structure a carbon atom
        is added as placeholder at 0,0,0.

        :param filename: (str) The filename to read from.
        :param primitive: (bool) Whether to convert to a primitive cell
        Only available for cifs. Defaults to False.
        :param sort: (bool) Whether to sort sites. Default to False.
        :param merge_tol: (float) If this is some positive number, sites that
        are within merge_tol from each other will be merged. Usually 0.01
        should be enough to deal with common numerical issues.
        :return: Return the created MofStructure
        """
        mof_name = os.path.splitext(os.path.basename(filename))[0]
        try:
            try:
                atoms = read(filename)
                s = AseAtomsAdaptor.get_structure(atoms)
            except:
                try:
                   s = Structure.from_file(filename, primitive=primitive, sort=sort,
                                        merge_tol=merge_tol)
                except:
                    s = CifParser(filename, occupancy_tolerance=10)
                    s.get_structures()
            try:
                s_mof = cls(s.lattice, s.species, s.frac_coords, name=mof_name)
            except:
                s_mof = cls(s.lattice, s.specie, s.frac_coords, name=mof_name)
            s_mof.summary['cif_okay'] = True
            s_mof.summary['checksum'] = Helper.get_checksum(filename)
        except Exception as e:
            print('\nAn Exception occurred: {}'.format(e))
            print('Cannot load {}\n'.format(filename))
            # Make a placeholder MOF object, set all its summary entries
            # to None and set cif_okay to False
            s_mof = cls([[10, 0, 0], [0, 10, 0], [0, 0, 10]],
                        ["C"], [[0, 0, 0]], name=mof_name)
            s_mof._mark_failed_to_read()
            s_mof.summary['cif_okay'] = False

        return s_mof

    def analyze_metals(self, output_folder, verbose='normal'):
        """Run analysis to detect all open metal sites in a MofStructure. In
        addition the metal sites are marked as unique.

        :param output_folder: Folder where OMS analysis results will be stored.
        :param verbose: Verbosity level for the output of the analysis.
        """

        Helper.make_folder(output_folder)
        running_indicator = output_folder + "/analysis_running"
        open(running_indicator, 'w').close()

        self.summary['problematic'] = False

        ms_cs_list = {True: [], False: []}
        for m, omc in enumerate(self.metal_coord_spheres):
            m_index = self.metal_indices[m]
            omc.check_if_open()
            if not self.summary['problematic']:
                self.summary['problematic'] = omc.is_problematic

            cs = self._find_coordination_sequence(m_index)
            cs = [self.species_str[m_index]] + cs
            omc.is_unique = self._check_if_new_site(ms_cs_list[omc.is_open], cs)
            if omc.is_unique:
                ms_cs_list[omc.is_open].append(cs)

            self.summary['metal_sites'].append(omc.metal_summary)

        unique_sites = [s['unique'] for s in self.summary['metal_sites']]
        open_sites = [s['is_open'] for s in self.summary['metal_sites']]

        self.summary['oms_density'] = sum(unique_sites) / self.volume
        self.summary['has_oms'] = any(open_sites)

        self.write_results(output_folder, verbose)
        os.remove(running_indicator)

    def write_results(self, output_folder, verbose='normal'):
        """Store summary dictionary holding all MOF and OMS information to a
        JSON file, store CIF files for the metal and non-metal parts of the MOF
        as well as all the identified coordination spheres.

        :param output_folder: Location to be used to store
        :param verbose: Verbosity level (default: 'normal')
        """
        Helper.make_folder(output_folder)
        for index, mcs in enumerate(self.metal_coord_spheres):
            mcs.write_cif_file(output_folder, index)
        if self.metal:
            output_fname = "{}/{}_metal.cif".format(output_folder,
                                                    self.summary['name'])
            self.metal.to(filename=output_fname)
        output_fname = "{}/{}_organic.cif".format(output_folder,
                                                  self.summary['name'])
        self.organic.to(filename=output_fname)

        json_file_out = "{}/{}.json".format(output_folder, self.summary['name'])
        summary = copy.deepcopy(self.summary)
        if verbose == 'normal':
            for ms in summary["metal_sites"]:
                ms.pop('all_dihedrals', None)
                ms.pop('min_dihedral', None)
        with open(json_file_out, 'w') as outfile:
            json.dump(summary, outfile, indent=3)

    @property
    def tolerance(self):
        """Tolerance values for dihedral checks. If not set, defaults are given.
        """
        if self._tolerance is None:
            self._tolerance = {'on_plane': 15}
        return self._tolerance

    @property
    def name(self):
        """Name of the MofStructure."""
        return self._name

    @name.setter
    def name(self, name):
        """Setter for the name of the MofStructure."""
        self._name = name
        self.summary['name'] = name

    @property
    def all_distances(self):
        """Distances between all atoms in the MofStructure"""
        if self._all_distances is None:
            self._all_distances = self.lattice.get_all_distances(
                self.frac_coords, self.frac_coords)
        return self._all_distances

    @property
    def all_coord_spheres_indices(self):
        """Compute the indices of the atoms in the first coordination shell
        for all atoms in the MofStructure
        """
        if self._all_coord_spheres_indices:
            return self._all_coord_spheres_indices

        self._all_coord_spheres_indices = [self._find_cs_indices(i)
                                           for i in range(len(self))]
        return self._all_coord_spheres_indices

    @property
    def metal_coord_spheres(self):
        """For all metal atoms in a MofStructure compute the first coordination
        sphere as a MetalSite object.
        """
        if not self._metal_coord_spheres:
            self._metal_coord_spheres = [self._find_metal_coord_sphere(c)
                                         for c in self.metal_indices]
        return self._metal_coord_spheres

    def _mark_failed_to_read(self):
        """If a CIF cannot be read set certain properties to None"""
        self.summary['metal_species'] = None
        self.summary['non_metal_species'] = None
        self.summary['uc_volume'] = None
        self.summary['density'] = None

    def _split_structure_to_organic_and_metal(self):
        """Split a MOF to two pymatgen Structures, one containing only metal
         atoms and one containing only non-metal atoms."""
        self.metal = Structure(self.lattice, [], [])
        self.organic = Structure(self.lattice, [], [])
        i = 0
        for s, fc in zip(self.species, self.frac_coords):
            if Atom(str(s)).is_metal:
                self.metal.append(s, fc)
                self.metal_indices.append(i)
            else:
                self.organic.append(s, fc)
            i += 1

    def _find_cs_indices(self, center):
        """Find the indices of the atoms in the coordination sphere.

        :param center: Central atom of coordination sphere.
        :return: c_sphere_indices: Return in the coordination sphere of center.
        """
        dist = list(self.all_distances[center])
        if dist[center] > 0.0000001:
            sys.exit('The self distance appears to be non-zero')

        a = Atom(self.species_str[center])
        c_sphere_indices = [i for i, dis in enumerate(dist)
                            if i != center
                            and a.check_bond(self.species_str[i], dis)]
        c_sphere_indices.insert(0, center)
        return c_sphere_indices

    def _find_metal_coord_sphere(self, center):
        """Identify the atoms in the first coordination sphere of a metal atom.

        Obtain all atoms connecting to the metal using the
        all_coord_spheres_indices values and keeping only valid bonds as well as
        center the atoms around the metal center for visualization purposes.

        :param center:
        :return:
        """
        dist = self.all_distances[center]
        if dist[center] > 0.0000001:
            sys.exit('The self distance appears to be non-zero')

        c_sphere = MetalSite(self.lattice, [self.species[center]],
                             [self.frac_coords[center]],
                             tolerance=self.tolerance)

        cs_i = self.all_coord_spheres_indices[center]
        for i in cs_i[1:]:
            c_sphere.append(self.species_str[i], self.frac_coords[i])
        c_sphere.keep_valid_bonds()
        c_sphere.center_around_metal()
        return c_sphere

    @staticmethod
    def _check_if_new_site(cs_list, cs):
        """Check if a given site is unique based on its coordination sequence"""
        for cs_i in cs_list:
            if Helper.compare_lists(cs_i, cs):
                return False
        return True  # len(cs_list),

    def _find_coordination_sequence(self, center):
        """Compute the coordination sequence up to the 6th coordination shell.

        :param center: Atom to compute coordination sequence for
        :return cs: Coordination sequence for center
        """

        shell_list = {(center, (0, 0, 0))}
        shell_list_prev = set([])
        all_shells = set(shell_list)
        n_shells = 6
        cs = []
        count_total = 0
        for n in range(0, n_shells):
            c_set = set([])
            for a_uc in shell_list:
                a = a_uc[0]
                lattice = a_uc[1]
                coord_sphere = self.all_coord_spheres_indices[a]
                count_total += 1
                coord_sphere_with_uc = []
                for c in coord_sphere:
                    diff = self.frac_coords[a] - self.frac_coords[c]
                    new_lat_i = [round(d, 0) for d in diff]
                    uc = tuple(l-nl for l, nl in zip(lattice, new_lat_i))
                    coord_sphere_with_uc.append((c, uc))
                coord_sphere_with_uc = tuple(coord_sphere_with_uc)
                c_set = c_set.union(set(coord_sphere_with_uc))
            for a in shell_list_prev:
                c_set.discard(a)
            for a in shell_list:
                c_set.discard(a)

            cs.append(len(c_set))
            all_shells = all_shells.union(c_set)
            shell_list_prev = shell_list
            shell_list = c_set

        return cs


class MetalSite(MofStructure):

    def __init__(self, lattice, species, coords, validate_proximity=False,
                 to_unit_cell=False, coords_are_cartesian=False,
                 site_properties=None, tolerance=None, name='N/A'):
        super().__init__(lattice, species, coords,
                         validate_proximity=validate_proximity,
                         to_unit_cell=to_unit_cell,
                         coords_are_cartesian=coords_are_cartesian,
                         site_properties=site_properties,
                         name=name)

        self._metal_type = "unknown"
        self._tolerance = tolerance
        self._is_open = None
        self._is_unique = None
        self._is_problematic = None
        self._t_factor = None
        self._min_dihedral = None
        self._all_dihedrals = {}

    @property
    def tolerance(self):
        """Tolerance values for dihedral checks. If not set, defaults are given.
        """
        if self._tolerance is None:
            self._tolerance = {'on_plane': 15}
        return self._tolerance

    @property
    def num_linkers(self):
        """Number of linkers in coordination sphere of MetalSite."""
        return self.num_sites - 1

    @property
    def is_open(self):
        """Whether the MetalSite is open or not."""
        return self._is_open

    @property
    def is_problematic(self):
        """Whether the MetalSite is problematic or not."""
        return self._is_problematic

    @property
    def is_unique(self):
        """Whether the MetalSite is unique or not."""
        return self._is_unique

    @is_unique.setter
    def is_unique(self, value):
        if not isinstance(value, bool):
            sys.exit('is_unique can only be boolean.')
        self._is_unique = value

    @property
    def metal_type(self):
        """The type of the metal center."""
        return self._metal_type

    @property
    def metal_summary(self):
        """Whether the MetalSite is problematic or not."""

        _summary = {"metal": str(self.species[0]),
                    "type": self.metal_type,
                    "is_open": self.is_open,
                    "unique": self.is_unique,
                    "problematic": self.is_problematic,
                    "number_of_linkers": self.num_linkers,
                    "min_dihedral": 0.0,
                    "all_dihedrals": 0.0,
                    't_factor': self._t_factor}

        return _summary

    def keep_valid_bonds(self):
        """Loop over atoms in the coordination sphere and remove any extraneous
        sites.
        """
        if len(self) == 0:
            return
        all_dists = self.lattice.get_all_distances(self.frac_coords,
                                                   self.frac_coords)
        for i, j in itertools.combinations(range(1, len(self)), 2):
            assert all_dists[i][j] == all_dists[j][i]
            dis = all_dists[i][j]
            if not self._valid_pair(i, j, dis):
                dist_ij_c = [all_dists[i][0], all_dists[j][0]]
                if len(set(dist_ij_c)) == 1:
                    index_to_remove = i
                else:
                    index_to_remove = [i, j][dist_ij_c.index(max(dist_ij_c))]
                self.remove_sites([index_to_remove])
                return self.keep_valid_bonds()

    def center_around_metal(self):
        """Shift atoms across periodic boundary conditions to have the
        coordination appear centered around the metal atom for visualisation
        purposes
        """
        gc = self.lattice.get_cartesian_coords
        center = self.frac_coords[0]
        center_cart_coords = gc(center)
        for i in range(1, self.num_sites):
            c_i = self.frac_coords[i]
            dist_vector = center - c_i
            dist_vector_r = []
            for j in range(0, 3):
                dist_vector_r.append(round(dist_vector[j]))
            dist_before = np.linalg.norm(center_cart_coords - gc(c_i))
            c_i_centered = c_i + dist_vector_r
            dist_after = np.linalg.norm(center_cart_coords - gc(c_i_centered))
            if dist_after > dist_before:
                for j in range(0, 3):
                    dist_vector_r[j] = np.rint(dist_vector[j])
                c_i_centered = c_i + dist_vector_r
                if dist_after > dist_before:
                    c_i_centered = c_i
            self.replace(i, self.species[i], c_i_centered)

    def check_if_open(self):
        """Get t-factor, check if problematic based on number of linkers and
         if necessary call to check the dihedrals to determine if the metal site
         is open.
         """

        self.get_t_factor()

        if Atom(str(self.species[0])).is_lanthanide_or_actinide:
            self._is_problematic = self.num_linkers < 5
        else:
            self._is_problematic = self.num_linkers < 3

        self._is_open = False
        self._metal_type = "Closed"
        if self.num_linkers <= 3:
            self._mark_oms(oms_type='3_or_less')
            return
        else:
            # 0 should always correspond to the
            self._check_planes(0)

    def _mark_oms(self, oms_type):
        self._metal_type = oms_type
        self._is_open = True

    def get_t_factor(self):
        """Compute t-factors, only meaningful for 4-,5-, and 6-coordinated
        metals, if not the value of -1 is assigned.
        """
        nl = self.num_sites - 1
        index_range = range(1, self.num_sites)
        all_angles = []
        for i in itertools.combinations(index_range, 2):
            angle = self.get_angle(i[0], 0, i[1])
            all_angles.append([angle, i[0], i[1]])

        all_angles.sort(key=lambda x: x[0])
        if nl == 5 or nl == 4:
            # beta is the largest angle and alpha is the second largest angle
            # in the coordination sphere; using the same convention
            # as Yang et al. DOI: 10.1039/b617136b
            beta = all_angles[-1][0]
            alpha = all_angles[-2][0]
            if nl == 4:
                tau = self.get_t4_factor(alpha, beta)
            else:
                tau = self.get_t5_factor(alpha, beta)
        elif nl == 6:
            max_indices_all = all_angles[-1][1:3]
            l3_l4_angles = [x for x in all_angles if
                            x[1] not in max_indices_all and
                            x[2] not in max_indices_all]
            max_indices_all_3_4 = max(l3_l4_angles, key=lambda x: x[0])[1:3]
            l5_l6_angles = [x for x in l3_l4_angles
                            if x[1] not in max_indices_all_3_4 and
                            x[2] not in max_indices_all_3_4]
            gamma = max(l5_l6_angles, key=lambda x: x[0])[0]
            tau = self.get_t6_factor(gamma)
        else:
            tau = -1
        self._t_factor = tau

    @staticmethod
    def get_t4_factor(a, b):
        return (360 - (a + b)) / 141.0

    @staticmethod
    def get_t5_factor(a, b):
        return (b - a) / 60.0

    @staticmethod
    def get_t6_factor(c):
        return c / 180.0

    def write_cif_file(self, output_folder, index):
        """Write MofSite to specified output_folder as a CIF file and use index
        to name it.
        """
        Helper.make_folder(output_folder)
        output_fname = output_folder
        output_fname += '/first_coordination_sphere'+str(index)+'.cif'
        self.to(filename=output_fname)

    def _valid_pair(self, i, j, dis):
        """Determine whether two atoms in the coordination sphere form a valid
        pair.

        A pair is not valid if it forms a bond unless both atoms are metals of
        the same kind as the center or both atoms are carbon atoms (e.g. in the
        case of a ferocene type coordination sphere).

        :param i:
        :param j:
        :param dis:
        :return:
        """
        s_one = str(self.species[i])
        s_two = str(self.species[j])
        a_one = Atom(s_one)
        a_two = Atom(s_two)

        bond = a_one.check_bond(s_two, dis, a_one.bond_tolerance(s_two))

        same_atoms = s_one == s_two == str(self.species[0])
        two_same_metals = same_atoms and a_one.is_metal and a_two.is_metal

        carbon_atoms = s_one == s_two == 'C'

        return (not bond) or two_same_metals or carbon_atoms

    def _check_planes(self, site):
        """Determine whether a site is open using the dihedral angles
        between the atoms in the coordination sphere.
        :param site: Index of site to be checked.
        """

        for i, j, k in itertools.combinations(range(self.num_sites), 3):
            plane = self._compute_plane(i, j, k)
            if all([abs(p-0.0) < 1e-5 for p in plane]):
                continue
            sides = self._sides([i, j, k], plane)
            # Side of the site in question.
            s_site = sides[site]
            # All sites that are not on the plane and are not the site in
            # question.
            s_o = [s for i, s in enumerate(sides) if i != site and s != 0]
            # Keep only the unique sides
            s_o_unique = list(set(s_o))
            # Number of unique sides for other sites (sites not on plane and
            # not the site in question)
            ls = len(s_o_unique)
            # ls = 0 : all other sites are on the plane
            # ls = 1 : all other sites on one side of plane
            # ls = 2 : other sites are on both sides of plane
            if ls == 0 or (ls == 1 and s_site != s_o_unique[0]):
                # Site is open if:
                # a) If all other sites fall on the plane. (ls == 0)
                # b) The metal site falls on the plane and all other sites
                # fall on one side of the plane. (ls == 1, s_site == 0 and
                # s_site != s_o_unique[0])
                # c) The metal site falls on one side of the plane and all
                # other sites fall on the oposite side of the plane.  (ls == 1,
                # s_site == 1,-1 and s_site != s_o_unique[0])
                place = {0: "over", 1: "on"}[abs(s_site)]
                msg = "{}_{}L_{}_open_plane".format(self[site].specie,
                                                    self.num_linkers, place)
                self._mark_oms(msg)
                break
            assert self.is_open is False

    def _sides(self, p_i, plane):
        """Given a plane p defined by 3 of the atoms in the MetalSite determine
        on which side of the plane all the atoms in the MetalSite fall (-1 or 1)
        or if it falls on the plane (0).

        :param p_i: Indices of the 3 atoms that define the plane
        :param plane: Plane constants
        :return: List of side value for all atoms in the MetalSite, possible
        values can -1, 0, and 1.
        """
        atoms_on_plane = [True if i in p_i
                          else self._is_point_on_plane(self[i].coords, p_i,
                                                       plane)
                          for i in range(len(self))]

        dists = [self._get_distance_from_plane(s.coords, plane) for s in self]
        sides = [0 if a or d == 0.0
                 else int(d/abs(d))
                 for d, a in zip(dists, atoms_on_plane)]
        return sides

    def _is_point_on_plane(self, point, p_i, p):
        """Given a point and plane determine if the point falls on the plane,
        using the angle between the projection of the point, each atom on the
        plane and the actual position of the point with a specified tolerance
        value.
        :param point: Cartesian coordinates of point to check.
        :param p: plane in the form of a list with the 4 constants defining
        a plane.
        :return: True if the point falls on plane and False otherwise.

        """
        tol = self.tolerance['on_plane']
        point_on_plane = self._project_point_onto_plane(point, p)
        angles = [self._get_angle_c(point_on_plane, self[ii].coords, point)
                  for ii in p_i]
        return all([a < tol for a in angles])

    def _get_angle_c(self, c1, c2, c3):
        """
        Calculates the angle between three points in degrees.

        :param c1: Coordinates of first point.
        :param c2: Coordinates of second point.
        :param c3: Coordinates of third point.
        :return: Angle between them in degrees.
        """
        v1 = c1 - c2
        v2 = c3 - c2
        return self._get_angle_v(v1, v2)

    @staticmethod
    def _get_angle_v(v1, v2):
        """
        Calculates the angle between two vectors in degrees.

        :param v1: First vector.
        :param v2: Second vector.
        :return: Angle between them in degrees.
        """
        if np.dot(v1, v2) == 0.0:
            return 0.0
        d = np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)
        d = min(d, 1.0)
        d = max(d, -1.0)
        angle = math.acos(d)
        return math.degrees(angle)

    @staticmethod
    def _get_distance_from_plane(point, plane):
        """Given a point and a plane compute the distance between the point and
        the projection of the point on the plane."""
        plane_xyz = plane[0:3]
        distance = np.inner(plane_xyz, point) - plane[3]
        return distance / np.linalg.norm(plane_xyz)

    def _compute_plane(self, i, j, k):
        """Given three atom indices, compute the plane that passes through them.
        """
        c1 = self[i].coords
        c2 = self[j].coords
        c3 = self[k].coords
        return self._compute_plane_c(c1, c2, c3)

    @staticmethod
    def _compute_plane_c(c1, c2, c3):
        """Given three atom coordinates, compute the plane that passes
        through them.
        """
        ij = c1 - c2
        kj = c3 - c2
        p_vector = np.cross(ij, kj)
        c = np.dot(c1, p_vector)
        plane = list(p_vector) + [c]
        return plane

    @staticmethod
    def _project_point_onto_plane(point, plane):
        """Given a point and plane compute the projection of the point onto the
        plane.
        """
        vector = plane[0:3]
        constant = plane[3]
        nom = np.inner(vector, point) - constant
        denom = np.inner(vector, vector)
        const = nom / denom
        return np.array([po - v * const for po, v in zip(point, vector)])


class Helper:

    @classmethod
    def make_folder(cls, output_folder):
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

    @classmethod
    def compare_lists(cls, l1, l2):
        if len(l1) != len(l2):
            return False
        return all([i == j for i, j in zip(l1, l2)])

    @classmethod
    def copy_folder(cls, dest, src):
        if not os.path.exists(dest):
            os.makedirs(dest)
        s = src.split('/')[-1]
        d = os.path.join(dest, s)
        if not os.path.exists(d):
            shutil.copytree(src, d)

    @classmethod
    def get_checksum(cls, filename):
        with open(filename, 'rb') as f:
            file = f.read()
        return hashlib.sha256(file).hexdigest()


class MofCollection:
    """A collection to hold and analyse MOF structures from CIF files"""

    separator = "".join(['-'] * 50)

    def __init__(self, path_list, analysis_folder='analysis_folder'):
        """Create a MofCollection from a list of path names.

        :param path_list: List of paths to MOF CIF files to be added to the
        collection.
        :param analysis_folder: Path to the folder where the results will
        be stored. (default: 'analysis_folder')
        """
        self._analysis_folder = analysis_folder
        self.path_list = path_list
        self.mof_coll = []
        self.batches = []
        self._metal_site_df = None
        self._mof_oms_df = None
        self._properties = {}
        self.load_balance_index = {}
        self.analysis_limit = None

        self.filter_functions = {
            "density": self._apply_filter_range,
            "oms_density": self._apply_filter_range,
            "uc_volume": self._apply_filter_range,
            "metal_species": self._apply_filter_in_value,
            "non_metal_species": self._apply_filter_in_value,
            "cif_okay": self._apply_filter_value,
            "has_oms": self._apply_filter_value,
            "mof_name": self._apply_value_in_filter
        }

        self._load_mofs()

    def __len__(self):
        return len(self.mof_coll)

    def __repr__(self):
        print_str = self.separator
        print_str += "\nThis collection holds information for "
        print_str += "{} MOFs.\n".format(len(self))
        if self.analysis_folder is None:
            print_str += "Analysis folder is not set.\n"
        else:
            f = os.path.abspath(self.analysis_folder)
            print_str += "Analysis folder is: {}\n\n".format(f)
        print_str += "List of cif files in collection:\n\n"
        for mc in self.mof_coll:
            print_str += "{}\n".format(mc['mof_file'])
        print_str += self.separator

        return print_str

    @property
    def analysis_folder(self):
        """Get value of the analysis folder."""
        Helper.make_folder(self._analysis_folder)
        return self._analysis_folder

    @analysis_folder.setter
    def analysis_folder(self, analysis_folder):
        """Set value of the analysis folder."""
        self._analysis_folder = analysis_folder

    @property
    def oms_results_folder(self):
        """Get value of the OMS results folder."""
        orf = self.analysis_folder + '/oms_results'
        Helper.make_folder(orf)
        return orf

    @property
    def summary_folder(self):
        """Get value of the summary folder."""
        sf = self.analysis_folder + '/summary'
        Helper.make_folder(sf)
        return sf

    @property
    def _properties_filename(self):
        """Get value of the properties pickle file."""
        return self.analysis_folder + '/properties.pickle'

    @property
    def properties(self):
        """Get value for the MOF properties. If the property variable is not
        None and the pickle file exists, then load the file and return it."""
        if not self._properties and os.path.isfile(self._properties_filename):
            with open(self._properties_filename, 'rb') as properties_file:
                self._properties = pickle.load(properties_file)
        return self._properties

    @property
    def mof_oms_df(self):
        """Get a pandas DataFrame that lists for each MOF whether it has an OMS
        or not and if it has an OMS what metal types it is.
        """
        if self._mof_oms_df is not None:
            return self._mof_oms_df
        if not self._validate_properties(['has_oms'])[1]:
            print('OMS analysis not finished for all MOFs in collection.')
            return False
        mof_info = {}
        for mi in self.mof_coll:
            mp = self.properties[mi['checksum']]
            if 'metal_sites' not in mp:
                continue
            metal_sites = mp['metal_sites']
            if len(metal_sites) == 0:
                print('No Metal Found in {}'.format(mp['name']))
            oms_types = [ms["metal"] for ms in metal_sites
                         if ms["is_open"] and ms["unique"]]
            oms_types = list(set(oms_types))
            if oms_types:
                oms_types = ",".join(oms_types)
            else:
                oms_types = "N/A"
            if mp['has_oms']:
                has_oms = 'Yes'
            else:
                has_oms = 'No'
            try:
                all_metal_species = ",".join(set(mp['metal_species']))
            except:
                all_metal_species="unknown"
            mof_info[mp['name']] = {'Metal Types': all_metal_species,
                                    'Has OMS': has_oms,
                                    'OMS Types': oms_types}
        self._metal_site_df = pd.DataFrame.from_dict(mof_info,
                                                     orient='index')
        return self._metal_site_df

    @property
    def metal_site_df(self):
        """Get a pandas DataFrame that lists the OMS results for each metal
        type.
        """
        if self._metal_site_df is not None:
            return self._metal_site_df
        if not self._validate_properties(['has_oms'])[1]:
            print('OMS analysis not finished for all MOFs in collection.')
            return False
        site_info = {}
        for mi in self.mof_coll:
            mp = self.properties[mi['checksum']]
            if 'metal_sites' not in mp:
                continue
            metal_sites = mp['metal_sites']
            if len(metal_sites) == 0:
                print('No Metal Found in {}'.format(mp['name']))
            for i, ms in enumerate(metal_sites):
                key = mp['name'] + '_' + str(i)
                site_info[key] = ms
                if 'all_dihedrals' in ms:
                    del site_info[key]['all_dihedrals']
                if 'min_dihedral' in ms:
                    del site_info[key]['min_dihedral']
                site_info[key]['mof_name'] = mp['name']
        self._metal_site_df = pd.DataFrame.from_dict(site_info, orient='index')
        return self._metal_site_df

    @classmethod
    def from_folder(cls, collection_folder, analysis_folder='analysis_folder',
                    name_list=None):
        """Create a MofCollection from a the CIF files in a folder.

        :param collection_folder: Path to the folder containing the CIF files to
        be added to the collection.
        :param analysis_folder: Path to the folder where the results will
        be stored. (default: 'analysis_folder')
        :param name_list: List of MOF names to include in the collection. If
        set, all the other CIF files in the folder will be excluded.
        (default: None)
        :return: A MofCollection object holding the specified MOF structures.
        """

        if name_list:
            print(cls.separator)
            print('Using only MOFs in the name list.')
            print(cls.separator)
            d = collection_folder
            path_list = [d+'/'+name for name in name_list]
        else:
            path_list = glob.glob(collection_folder + "/*.cif")
        return cls(path_list, analysis_folder)

    def analyse_mofs(self, overwrite=False, num_batches=1, analysis_limit=None):
        """Run OMS analysis for the MOFs in the collection.

        :param overwrite: Controls if the results will be overwritten or not
        (default: False)
        :param num_batches: Sets the number of batches the structures will be
        split in and analyzed on a separate process. (default: 1)
        :param analysis_limit: Analyze only up to the number of MOFs set by
        analysis_limit, if set to None all MOFs will be analyzed (default: None)
        """
        print(self.separator)
        print("Running OMS Analysis...")
        self.analysis_limit = analysis_limit

        t0 = time.time()

        self._make_batches(num_batches, overwrite)

        status = Array('i', [0 for i in range(num_batches)])
        for i, batch in enumerate(self.batches):
            p = Process(target=self._run_batch,
                        args=(i, batch, overwrite,status))
            p.start()

        lbs = [len(batch)/100.0 for batch in self.batches]
        wait_time = 0.0
        status_prev = [0 for i in range(num_batches)]
        while True:
            # Create a list from the shared array to make sure it doesnt change
            # during the iteration
            status_ = list(status)
            if all([sp == s for sp, s in zip(status_prev, status_)]):
                wait_time = min(25, 0.1+wait_time)
                time.sleep(wait_time)
            status_prev = status_

            sout = ["Batch {} Finished.".format(b + 1)
                    if len(self.batches[b]) == 0 or s < 0 else
                    "Batch {} {:.2f} % : Analysing {:}"
                    "".format(b+1, (s+1)/lbs[b], self.batches[b][s]['mof_name'])
                    for b, s in enumerate(status_)]
            print("|**| ".join(sout) + 100 * " ", end='\r', flush=True)

            if all([s < 0 for s in status_]):
                break

        if overwrite:
            for mi in self.mof_coll:
                self._update_property_from_oms_result(mi)
        self._validate_properties(['has_oms'])

        t1 = time.time()
        print('\nAnalysis Finished. Time required:{:.2f} sec'.format(t1 - t0))
        print(self.separator)

    def check_structures(self):
        """Iterate over all the MOFs in the collection and validate that they
        can be read and a MofStructure can be created.
        """
        self._validate_properties(['cif_okay'])
        not_read = [mi for mi in self.mof_coll
                    if not self.properties[mi['checksum']]['cif_okay']]
        read_len = len(self.mof_coll) - len(not_read)
        print('\nChecked {} structures.'.format(len(self.mof_coll)))
        msg1 = {0: '\r',
                1: '{} was read.'.format(read_len),
                2: '{} were read.'.format(read_len)}
        msg2 = {0: '\r',
                1: '{} was NOT read.'.format(len(not_read)),
                2: '{} were NOT read.'.format(len(not_read))}
        print(msg1[min(2, read_len)])
        print(msg2[min(2, len(not_read))])

        msg = {0: "\r", 1: "\nThe following structures could not be read:"}
        print(msg[min(1, len(not_read))])
        for i, mi in enumerate(not_read):
            print("{}".format(mi['mof_name']))

        mofs_no_metal = [mi for mi in self.mof_coll
                         if self.properties[mi['checksum']]['cif_okay']
                         and not
                         self.properties[mi['checksum']]['metal_species']]
        msg = {0: "\r", 1: "The following structures contain no metal:"}
        print(msg[min(1, len(mofs_no_metal))])
        for mi in mofs_no_metal:
            p = self.properties[mi['checksum']]
            print("{}.cif {}".format(p['name'],
                                     p['metal_species']+p['non_metal_species']))

        print('\nFinished checking structures.')

    def check_analysis_status(self):
        """Iterate over all the MOFs in the collection and check if the results
        from the OMS analysis exist.
        """
        print(self.separator)
        not_done = [mi['mof_file'] for mi in self.mof_coll
                    if not self._check_if_results_exist(mi['mof_name'])]
        done = len(self.mof_coll) - len(not_done)
        msg1 = {0: '\nAnalysis for no structures has been completed.',
                1: '\nAnalysis for {} out of {} structures have been completed.'
                   .format(done, len(self.mof_coll))}
        msg2 = {0: "\r", 1: "\nThe following structures are missing:"}

        print(msg1[min(1, done)])
        print(msg2[min(1, len(not_done))])
        for nd in not_done:
            print(nd)
        print(self.separator)

    def sample_collection(self, sample_size=50):
        """Randomly select a sample of MOFs in the collection and
        return a new collection with the MOFs in the sample.

        :param sample_size: Number of MOFs to be selected. Default value is 50.

        """
        ll = len(self.mof_coll)
        if sample_size > ll:
            sample_size = ll
            print(f"Can only sample up to the number of MOFs "
                  f"in the collection ({ll}).")
        mof_list = [mi['mof_file'] for mi in self.mof_coll]
        sampled_list = random.sample(mof_list, sample_size)
        return MofCollection(sampled_list, analysis_folder=self.analysis_folder)

    def filter_collection(self, using_filter=None,
                          new_collection_folder=None,
                          new_analysis_folder=None):
        """Filter a collection given a number of filters.

        Calling this method of a MofCollection applies the filter and creates a
        new collection for the MOFs that match the filter. The cif files that
        match the filter are  copied to the new_collection_folder.
        The filters can be one or more of the following:

        'density': [min, max] (range of values)
        'oms_density': [min, max] (range of values)
        'uc_volume':  [min, max] (range of values)
        'metal_species': ["Cu", "Zn", ...] (list of metal species)
        'non_metal_species': ["C", "N", ...] (list of non metal species)
        'cif_okay': True (boolean value)
        'has_oms': True (boolean value)
        'mof_name':  [mof_name1, mof_name2] (string values)

        :param using_filter: Filter used to identify MOFs with certain
        characteristics. Has to be a python dictionary (default: None)
        :param new_collection_folder: Path to the folder where the CIF files of
        the filtered collection will be stored. If set to None the CIF files
        will not be copied. (default: None)
        :param new_analysis_folder: Path to the folder where the OMS result
        files of the filtered collection will be stored. If set to None the
        result files will not be copied. (default: None)
        :return: A MofCollection with only the filtered MOFs. If
        new_collection_folder or new_analysis_folder is not set then the
        collection will point to the original location of these files.
        """
        print(self.separator)
        if any([f not in self.filter_functions for f in using_filter]):
            print('Unknown filter. Try again using one of the following '
                  'filters:\n\"{}\"'.format(", ".join(self.filter_functions)))
            print(self.separator)
            return

        validation_level, cf = self._validate_properties(using_filter)
        if validation_level == 1 and not cf:
            print('Properties from CIF files could not be validated.'
                  'Check that all CIF files can be read')
            return
        elif validation_level == 2 and not cf:
            print('Requested a filter that needs OMS information but the '
                  'OMS analysis does not appear to be complete.\n'
                  'Run it first and try again.')
            return

        print(self.separator)
        print('Filtering collection.')
        filtered_list = []
        for i, mi in enumerate(self.mof_coll):
            mp = self.properties[mi['checksum']]
            fun = self._apply_filter
            if all([fun(f, mp[f], using_filter[f]) for f in using_filter]):
                filtered_list.append(mi['mof_file'])

        found_s = {0: "No", 1: len(filtered_list)}[min(1, len(filtered_list))]
        print('\n{} MOFs were matched using the provided'
              ' filter.'.format(found_s))
        if len(filtered_list) == 0:
            print('No collection returned.')
            return None
        print('Returning a new collection using the matched MOFs.')
        sub_collection = MofCollection(filtered_list,
                                       analysis_folder=self.analysis_folder)
        print(self.separator)

        sub_collection.copy_cifs(new_collection_folder)
        sub_collection.copy_results(new_analysis_folder)

        return sub_collection

    def read_cif_files(self):
        """Iterate over all MOF files in the collection, load each CIF and
        store MOF properties such as density, unit cell volume etc.
        """
        print(self.separator)
        print('Reading CIF files and updating properties...')
        self._loop_over_collection(self._update_property_from_cif_file)
        self._store_properties()
        print('Done')
        print(self.separator)

    def read_oms_results(self):
        """Iterate over all MOF files in the collection, load each OMS result
        file and store OMS information to the MOF properties.
        """
        print(self.separator)
        print('Adding results to properties.')
        self._loop_over_collection(self._update_property_from_oms_result)
        print('Done')
        self._store_properties()
        print(self.separator)

    def copy_cifs(self, target_folder):
        """Copy cif files from their existing location to the specified
        target_folder.

        :param target_folder: Path of folder to copy collection CIF files to.
        """
        if target_folder is None:
            return
        tf_abspath = os.path.abspath(target_folder)
        Helper.make_folder(tf_abspath)
        print(self.separator)
        print('The cif files for this collection will be copied to'
              ' the specified folder:\n\"{}\"'.format(tf_abspath))
        print('The cif paths will be updated.')

        for i, mi in enumerate(list(self.mof_coll)):
            destination_path = "{}/{}.cif".format(tf_abspath, mi['mof_name'])
            self.mof_coll[i] = {"mof_name": mi['mof_name'],
                                "mof_file": destination_path,
                                "checksum": mi['checksum']}
            if not os.path.isfile(destination_path):
                shutil.copyfile(mi['mof_file'], destination_path)
        print(self.separator)

    def copy_results(self, target_folder):
        """Copy OMS result files from their existing location to the specified
        target_folder.

        :param target_folder: Path of folder to copy collection OMS result
        files to.
        """
        if target_folder is None:
            return

        print(self.separator)
        tf_abspath = os.path.abspath(target_folder)
        destination_path = tf_abspath + '/oms_results'

        print('The result files for this collection will be copied to the '
              'specified folder:\n{}\nThe analysis folder will be updated.'
              ''.format(tf_abspath))

        Helper.make_folder(tf_abspath)
        Helper.make_folder(destination_path)

        for i, mi in enumerate(self.mof_coll):
            mof_name = mi['mof_name']
            if self._check_if_results_exist(mof_name):
                source_path = "{}/{}".format(self.oms_results_folder, mof_name)
                Helper.copy_folder(destination_path, source_path)
        self.analysis_folder = tf_abspath
        self._validate_properties(['has_oms'])
        print(self.separator)

    def summarize_results(self, max_atomic_number=None):
        """Create a summary table for the OMS results of the collection, group
        results by metal type.

        :param max_atomic_number: Maximum atomic number to be included in
        summary table. If not defined all metal atoms will be considered
        (default: None)
        """
        df = self.metal_site_df.copy()
        site_df_u = df.loc[df['unique']]
        site_df_o = site_df_u.loc[site_df_u['is_open']]

        all_sites = self._group_and_summarize(site_df_u, ['MOFs',
                                                          'Metal Sites'])
        open_sites = self._group_and_summarize(site_df_o, ['MOFs_with_OMS',
                                                           'OMS'])

        s_df = pd.concat([all_sites, open_sites], axis=1)
        s_df.fillna(0.0, inplace=True)
        s_df = s_df.astype(int)

        s_df['MOFs_with_OMS(%)'] = 100.0 * s_df['MOFs_with_OMS']/s_df['MOFs']
        s_df['OMS (%)'] = 100.0 * s_df['OMS'] / s_df['Metal Sites']
        cols = ['MOFs', 'MOFs_with_OMS', 'Metal Sites', 'OMS',
                'MOFs_with_OMS(%)', 'OMS (%)']
        s_df = s_df[cols]

        s_df['MOFs_with_OMS(%)'] = s_df['MOFs_with_OMS(%)'].apply('{:.2f} %'
                                                                  ''.format)
        s_df['OMS (%)'] = s_df['OMS (%)'].apply('{:.2f} %'.format)
        s_df.sort_values("MOFs", inplace=True, ascending=False)

        num_mofs = df['mof_name'].nunique()
        num_oms_mofs = df[df['is_open']]['mof_name'].nunique()
        num_sites = len(site_df_u)
        num_oms_sites = len(site_df_u[site_df_u['is_open']])

        print(self.separator)
        print('Number of total MOFs: {}'.format(num_mofs))
        print('Number of total MOFs with open metal sites: {}'
              ''.format(num_oms_mofs))
        print('Number of total unique sites: {}'.format(num_sites))
        print('Number of total unique open metal sites: {}'
              ''.format(num_oms_sites))
        print(self.separator)

        msg = "Summary Table\n"
        fname = "{0}/stats.out".format(self.summary_folder, max_atomic_number)
        if max_atomic_number:
            subset = pd.Series(s_df.index).apply(
                lambda x: Atom(x).atomic_number <= max_atomic_number)
            s_df = s_df.loc[subset.values]
            fname = "{0}/stats_less_{1}.out".format(self.summary_folder,
                                                    max_atomic_number)
            msg = "Summary Table for metal atoms with atomic number smaller " \
                  "than {}.\n".format(max_atomic_number)
        print(msg)
        print(s_df)
        s_df.to_csv(fname, sep=' ')

    def summarize_tfactors(self):
        """Summarize the t-factor information and make histograms for all the
        MOFs in the collection.
        """
        tfac_analysis_folder = self.summary_folder + '/tfac_analysis'
        Helper.make_folder(self.summary_folder)
        Helper.make_folder(tfac_analysis_folder)

        df = self.metal_site_df.copy()
        sites_u = df[df['unique']]

        for n in range(4, 7):
            self._write_t_factors(sites_u, n, tfac_analysis_folder)

    def _load_mofs(self):
        """Add MOfs to collection, use CIF file checksum as an identifier."""
        print('Loading CIF files...')
        li = max(int(len(self.path_list) / 1000), 1)
        lm = len(self.path_list) / 100.0
        for i, mof_file in enumerate(self.path_list):
            if i % li == 0:
                print("{:4.1f} %".format((i+1) / lm), end="\r", flush=True)
            checksum = Helper.get_checksum(mof_file)
            mof_name = os.path.splitext(os.path.basename(mof_file))[0]
            mof_info = {"mof_name": mof_name,
                        "mof_file": mof_file,
                        "checksum": checksum}
            self.mof_coll.append(mof_info)
            if checksum not in self.properties:
                self.properties[checksum] = {"mof_name": mof_name}
            else:
                if self.properties[checksum]["mof_name"] != mof_name:
                    print("Warning: MOF name and CIF checksum mismatch for {}.cif "
                        "{}.cif. Either the CIF files has already been "
                        "processed with a different name, or the CIF file "
                        "has changed since it was processed. Continuing with processing."
                        "".format(mof_name,
                                    self.properties[checksum]['mof_name']))

            # else:
            #     if self.properties[checksum]["mof_name"] != mof_name:
            #         exit("MOF name and CIF checksum mismatch for {}.cif "
            #              "{}.cif. Either the CIF files has already been "
            #              "processed with a different name, or the CIF file "
            #              "has changed since it was processed."
            #              "".format(mof_name,
            #                        self.properties[checksum]['mof_name']))
            if self._check_if_results_exist(mof_name):
                self._compare_checksums(mof_file, mof_name, checksum)
        print("\nAll Done.")
        self._store_properties()

    def _compare_checksums(self, mof_file, mof_name, checksum):
        """If OMS results exist for one of the CIF names in the collection then
        ensure that the CIF checksum matches the one in the result file.
        """
        mof_folder = "{0}/{1}/".format(self.oms_results_folder,
                                       mof_name)
        results_file = "{0}/{1}.json".format(mof_folder, mof_name)
        with open(results_file, 'r') as f:
            results_dict = json.load(f)
        if results_dict['checksum'] != checksum:
            print("Results for a MOF named {0} appear to already exist"
                  " in the analysis folder \n\"{1}\".\nHowever the "
                  "file checksum in the result file does not match the "
                  "checksum of \n\"{2}\".\n\nHave the CIF files in the "
                  "collection changed since the results were computed?"
                  "\nClear results and try again.".format(mof_name,
                                                          mof_folder,
                                                          mof_file))
            exit(1)

    def _run_batch(self, b, batch, overwrite, status):
        """Run OMS analysis for each of the batches."""
        for i, mi in enumerate(batch):
            status[b] = i
            self._analyse(mi, overwrite)
        status[b] = -1

    def _analyse(self, mi, overwrite):
        """For a given CIF file, create MofStructure object and run OMS
        analysis. If overwrite is false check if results already exist first.
        """
        mof_folder = "{}/{}".format(self.oms_results_folder, mi['mof_name'])
        results_exist = self._check_if_results_exist(mi['mof_name'])
        if not overwrite and results_exist:
            print("Skipping {}. Results already exist and overwrite is set "
                  "to False.".format(mi['mof_name']))
            return
        mof = self._create_mof_from_cif_file(mi['mof_file'])
        if mof.summary['cif_okay']:
            mof.analyze_metals(output_folder=mof_folder)

    def _make_batches(self, num_batches=1, overwrite=False):
        """Split collection into number of batches

        :param num_batches: Number of batches (default: 1)
        :param overwrite: Controls if the results will be overwritten or not
        (default: False)
        """
        print(self.separator)
        if cpu_count() < num_batches:
            warnings.warn('You requested {} batches but there are only {}'
                          ' CPUs available.'.format(num_batches, cpu_count()))
        b_s = {1: 'batch', 2: 'batches'}[min(num_batches, 2)]
        print('{} {} requested. '.format(num_batches, b_s))
        print('Overwrite is set to {}. '.format(overwrite))
        print('Storing results in {}. '.format(self.oms_results_folder))
        print(self.separator)
        self._validate_properties(['load_balancing_index'])
        print(self.separator)
        lbi = {}
        for mi in self.mof_coll:
            mp = self.properties[mi['checksum']]
            lbi[mi['mof_name']] = mp['load_balancing_index']
        # Remove any structures not in load balancing index.
        subset = [mc for mc in self.mof_coll if mc['mof_name'] in lbi]

        # If there is no balancing info for a MOF at this point it means
        # that it could not be read.
        if len(self.mof_coll) != len(subset):
            print('\nSkipping {} structures that could not be read.'
                  ' '.format(len(self.mof_coll)-len(subset)))

        # Remove any structures already completed
        if not overwrite:
            print('Checking if results for any of the MOFs exist...')
            all_ = len(subset)
            subset = [mc for mc in subset if not
                      self._check_if_results_exist(mc['mof_name'])]
            msg = {0: "Will not skip any MOFs",
                   1: "Skipping {} MOFs because results were found. "
                      "".format(all_ - len(subset))}
            print(msg[min(1, all_ - len(subset))])

        # Sort mof list using the load balancing index
        subset.sort(key=lambda x: lbi[x['mof_name']])

        sum_load_balance = sum(lbi[mi["mof_name"]] for mi in subset)
        lb_per_batch = sum_load_balance / num_batches

        # Select only up to analysis_limit to work with
        if self.analysis_limit and len(subset) > self.analysis_limit:
            subset = subset[0:self.analysis_limit]

        self.batches = [[] for b in range(num_batches)]
        for i, mi in enumerate(subset):
            sum_lb = sum([lbi[mi["mof_name"]] for mi in subset[0:i]])
            batch = int(sum_lb / lb_per_batch)
            self.batches[batch].append(mi)
        print(self.separator)
        for i, batch in enumerate(self.batches):
            print("Batch {0} has {1} MOFs".format(i+1, len(batch)))
        print(self.separator)

    def _check_if_results_exist(self, mof_name):
        """Check if OMS results already exist for a MOF"""
        mof_folder = "{}/{}".format(self.oms_results_folder, mof_name)
        if os.path.isfile(mof_folder+'/'+mof_name+'.json'):
            if not os.path.isfile(mof_folder + '/' + 'analysis_running'):
                return True
        return False

    def _loop_over_collection(self, func):
        """Iterate over all the MOFs in the collection and run the specified
        function.
        :param func: Function to use.
        """
        li = max(int(len(self.mof_coll) / 1000), 1)
        lm = len(self.mof_coll) / 100
        for i, mi in enumerate(self.mof_coll):
            if i % li == 0:
                print("{:4.1f} % {} {:100}".format((i+1)/lm, mi['mof_name'],
                                                   " "), end="\r", flush=True)
            func(mi)
        print()

    def _apply_filter(self, filter_, v, f):
        """Apply the proper filter_function for the given filter"""
        return self.filter_functions[filter_](v, f)

    @staticmethod
    def _apply_filter_value(v, f):
        """Filter function to match a value. Returns false if values is None"""
        if not v:
            return False
        return v == f

    @staticmethod
    def _apply_filter_in_value(v, f):
        """Filter function to match all values of a list"""
        if not v:
            return False
        return all([f_ in v for f_ in f])

    @staticmethod
    def _apply_value_in_filter(v, f):
        """Filter function to match any of the values of a list"""
        if not v:
            return False
        return v in f

    @staticmethod
    def _apply_filter_range(v, f):
        """Filter function to match a range of values"""
        if not v:
            return False
        return min(f) <= v <= max(f)

    def _validate_properties(self, keys):
        """Check if a given property can be found in the properties dictionary.
        If not try to read the CIF file and check again. If the check fails
        again try to read the OMS results and check again. If the check fails
        a third time return False, the property cannot be validated."""
        msg = {1: "Validating property", 2: "Validating properties"}
        print('\n{} : '.format(msg[min(2, len(keys))]), end='')
        print("\"{}\"".format(", ".join([k for k in keys])))
        validation_level = 0
        li = max(int(len(self.mof_coll)/1000), 1)
        lm = len(self.mof_coll) / 100
        for i, mi in enumerate(self.mof_coll):
            if i % li == 0:
                print("{:4.1f} % {} {:100}".format((i+1) / lm, mi['mof_name'],
                                                   " "), end="\r", flush=True)
            mp = self.properties[mi['checksum']]
            if not self._validate_property(mp, keys):
                self._update_property_from_cif_file(mi)
                validation_level = 1
            if not self._validate_property(mp, keys):
                self._update_property_from_oms_result(mi)
                validation_level = 2
            if not self._validate_property(mp, keys):
                self._store_properties()
                print('\nProperty Missing\n{}'.format(self.separator))
                return validation_level, False
        self._store_properties()
        print("Validated 100 % "+100*" ", end="\r")
        print()
        return validation_level, True

    @staticmethod
    def _validate_property(mp, keys):
        """Check if property exists."""
        test1 = all([f in mp for f in keys])
        if test1 and all([mp[f] != 'N/A' for f in keys]):
            return True
        if test1 and not mp['cif_okay']:
            return True
        return False

    def _update_property_from_cif_file(self, mi):
        """Update properties dictionary from a CIF file."""
        mp = self.properties[mi['checksum']]
        mof = self._create_mof_from_cif_file(mi['mof_file'])
        if mof:
            mp.update(mof.summary)
            self.load_balance_index[mi['mof_name']] = len(mof) * len(mof)
            mp['load_balancing_index'] = self.load_balance_index[mi['mof_name']]

    def _update_property_from_oms_result(self, mi):
        """Update properties dictionary from an OMS result file."""
        mp = self.properties[mi['checksum']]
        mof_name = mp["mof_name"]
        mof_folder = "{0}/{1}/".format(self.oms_results_folder, mof_name)
        results_file = "{0}/{1}.json".format(mof_folder, mof_name)
        results_dict = None
        if os.path.isfile(results_file):
            results_dict = json.load(open(results_file))
        if isinstance(results_dict, dict):
            results_dict['source_name'] = mof_folder
            mp.update(results_dict)

    def _store_properties(self):
        """Store properties dictionary as a python pickle file."""
        with open(self._properties_filename, 'wb') as properties_file:
            pickle.dump(self._properties, properties_file)

    @staticmethod
    def _create_mof_from_cif_file(path_to_mof):
        """Create and return a MofStructure object from a path to a CIF file."""
        mof = MofStructure.from_file(path_to_mof, primitive=False)
        return mof

    def _write_t_factors(self, sites, n, target):
        """Summarize the findings in table form and histograms for a give
        t-factor.
        """
        s_n = sites.loc[sites['number_of_linkers'] == n].copy()
        s_n['is_open_yn'] = np.where(s_n['is_open'], 'yes', 'no')
        s_n = s_n[['mof_name', 'is_open_yn', 't_factor']]
        for flag in ['yes', 'no']:
            outpath = "{}/{}_{}.out".format(target, flag, str(n))
            s = s_n[s_n['is_open_yn'] == flag]
            s.to_csv(outpath, index=False)
            fout = "{}/{}_{}_hist.out".format(target, flag, n)
            self._write_histogram(s['t_factor'], True, fout)
            fout = "{}/{}_{}_hist_abs.out".format(target, flag, n)
            self._write_histogram(s['t_factor'], False, fout)

        fig = plt.figure(figsize=(10, 5))
        plt.title('t-{} factor'.format(n))
        s_yes = s_n[s_n['is_open_yn'] == 'yes']
        s_yes['t_factor'].hist(bins=50, range=(0, 1), density=False)
        s_no = s_n[s_n['is_open_yn'] == 'no']
        s_no['t_factor'].hist(bins=50, range=(0, 1), density=False)
        plt.show()

    @staticmethod
    def _write_histogram(sites, dens, target):
        """Generate histograms to be used for summarizing the t-factor
        results.
        """
        hist, edges = np.histogram(sites, bins=50, range=(0, 1), density=dens)
        with open(target, 'w') as hist_file:
            w = (edges[1] - edges[0]) / 2
            for e, h in zip(edges, hist):
                print(e + w, h, file=hist_file)

    @staticmethod
    def _group_and_summarize(df, names=None):
        """Group the DataFrame holding the OMS results by metal type and rename
        its columns.
        """
        rename = {"mof_name": names[0], "is_open": names[1]}
        agg_dict = {"mof_name": pd.Series.nunique, "is_open": "count"}
        return df.groupby('metal').agg(agg_dict).rename(columns=rename)
