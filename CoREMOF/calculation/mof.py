"""Classes used for open metal site analysis.
"""

import json
import copy
from pymatgen.core.structure import Structure
from CoREMOF.calculation.atomic_parameters import Atom
import numpy as np
import sys
import itertools
import os
import shutil
import hashlib
import datetime
import math


class MofStructure(Structure):

    """Extend the pymatgen Structure class to add MOF specific features. Create a MOf structure. The arguments are the same as in the pymatgen Structure class with the addition of the name argument.
    The super constructor is called and additional MOF specific properties are initialized.

    Args:
        Structure: MOF name, used to identify the structure."""

    def __init__(self, lattice, species, coords, charge=None,
                 validate_proximity=False, to_unit_cell=False,
                 coords_are_cartesian=False, site_properties=None, name="N/A"):

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
        """Create a MofStructure from a CIF file. This makes use of the from_file function of the Structure class and catches the exception in case a CIF file cannot be read. If the CIF is read successfully then the MofStructure is marked as okay, and the file checksum is added to the summary. If the CIF file cannot be read then it is marked as not okay and all the other properties are set to None and because there cannot be an empty Structure a carbon atom is added as placeholder at 0,0,0.

        Args:
            filename (str): The filename to read from.
            primitive (bool): Whether to convert to a primitive cell Only available for cifs. Defaults to False.
            sort (bool): Whether to sort sites. Default to False.
            merge_tol (float): If this is some positive number, sites that are within merge_tol from each other will be merged. Usually 0.01 should be enough to deal with common numerical issues.

        Returns:
            -   The created MofStructure
        """

        mof_name = os.path.splitext(os.path.basename(filename))[0]
        try:
            s = Structure.from_file(filename, primitive=primitive, sort=sort,
                                    merge_tol=merge_tol)
            s_mof = cls(s.lattice, s.species, s.frac_coords, name=mof_name)
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
        """Run analysis to detect all open metal sites in a MofStructure. In addition the metal sites are marked as unique.

        Args:
            output_folder (str): Folder where OMS analysis results will be stored.
            verbose: Verbosity level for the output of the analysis.
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

        """Store summary dictionary holding all MOF and OMS information to a JSON file, store CIF files for the metal and non-metal parts of the MOF as well as all the identified coordination spheres.

        Args:
            output_folder: Location to be used to store.
            verbose: Verbosity level.
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
        """Compute the indices of the atoms in the first coordination shell for all atoms in the MofStructure
        """
        if self._all_coord_spheres_indices:
            return self._all_coord_spheres_indices

        self._all_coord_spheres_indices = [self._find_cs_indices(i)
                                           for i in range(len(self))]
        return self._all_coord_spheres_indices

    @property
    def metal_coord_spheres(self):
        """For all metal atoms in a MofStructure compute the first coordination sphere as a MetalSite object.
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
        """Split a MOF to two pymatgen Structures, one containing only metal atoms and one containing only non-metal atoms."""
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

        Args:
            center: Central atom of coordination sphere.

        Returns:
            c_sphere_indices: Return in the coordination sphere of center.
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

        """Identify the atoms in the first coordination sphere of a metal atom. Obtain all atoms connecting to the metal using the all_coord_spheres_indices values and keeping only valid bonds as well as center the atoms around the metal center for visualization purposes.

        Args:
            center:
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

        Args:
            center: Atom to compute coordination sequence

        Returns:
            cs: Coordination sequence for center
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
        """Loop over atoms in the coordination sphere and remove any extraneous sites.
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
        """Shift atoms across periodic boundary conditions to have the coordination appear centered around the metal atom for visualisation purposes
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
        """Get t-factor, check if problematic based on number of linkers and if necessary call to check the dihedrals to determine if the metal site is open.
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
        """Compute t-factors, only meaningful for 4-,5-, and 6-coordinated metals, if not the value of -1 is assigned.
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
        """Write MofSite to specified output_folder as a CIF file and use index to name it.
        """
        Helper.make_folder(output_folder)
        output_fname = output_folder
        output_fname += '/first_coordination_sphere'+str(index)+'.cif'
        self.to(filename=output_fname)

    def _valid_pair(self, i, j, dis):

        """Determine whether two atoms in the coordination sphere form a valid pair. A pair is not valid if it forms a bond unless both atoms are metals of the same kind as the center or both atoms are carbon atoms (e.g. in the case of a ferocene type coordination sphere).

        Args:
            i:
            j:
            dis:
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

        """Determine whether a site is open using the dihedral angles between the atoms in the coordination sphere.

        Args:
            site: Index of site to be checked.
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

        """Given a plane p defined by 3 of the atoms in the MetalSite determine on which side of the plane all the atoms in the MetalSite fall (-1 or 1) or if it falls on the plane (0).

        Args:
            p_i: Indices of the 3 atoms that define the plane
            plane: Plane constants

        Returns:
            List of side value for all atoms in the MetalSite, possible values can -1, 0, and 1.
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

        """Given a point and plane determine if the point falls on the plane, using the angle between the projection of the point, each atom on the plane and the actual position of the point with a specified tolerance value.

        Args:
            point: Cartesian coordinates of point to check.
            p: plane in the form of a list with the 4 constants defining a plane.

        Returns:
            True if the point falls on plane and False otherwise.
        """

        tol = self.tolerance['on_plane']
        point_on_plane = self._project_point_onto_plane(point, p)
        angles = [self._get_angle_c(point_on_plane, self[ii].coords, point)
                  for ii in p_i]
        return all([a < tol for a in angles])

    def _get_angle_c(self, c1, c2, c3):

        """
        Calculates the angle between three points in degrees.

        Args:
            c1: Coordinates of first point.
            c2: Coordinates of second point.
            c3: Coordinates of third point.

        Returns:
            Angle between them in degrees.
        """

        v1 = c1 - c2
        v2 = c3 - c2
        return self._get_angle_v(v1, v2)

    @staticmethod
    def _get_angle_v(v1, v2):

        """
        Calculates the angle between two vectors in degrees.

        Args:
            v1: First vector.
            v2: Second vector.

        Returns:
            Angle between them in degrees.
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
        """Given a point and a plane compute the distance between the point and the projection of the point on the plane."""
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
        """Given three atom coordinates, compute the plane that passes through them.
        """
        ij = c1 - c2
        kj = c3 - c2
        p_vector = np.cross(ij, kj)
        c = np.dot(c1, p_vector)
        plane = list(p_vector) + [c]
        return plane

    @staticmethod
    def _project_point_onto_plane(point, plane):
        """Given a point and plane compute the projection of the point onto the plane.
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
