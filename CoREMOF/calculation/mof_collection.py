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
from CoREMOF.calculation.mof import Helper
from CoREMOF.calculation.mof import MofStructure
from CoREMOF.calculation.atomic_parameters import Atom
from sys import exit
pd.options.display.max_rows = 1000


class MofCollection:

    """A collection to hold and analyse MOF structures from CIF files. Create a MofCollection from a list of path names.

    Args:
        path_list: List of paths to MOF CIF files to be added to the collection.
        analysis_folder: Path to the folder where the results will be stored."""

    separator = "".join(['-'] * 50)

    def __init__(self, path_list, analysis_folder='analysis_folder'):
        """
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
        """Get value for the MOF properties. If the property variable is not None and the pickle file exists, then load the file and return it."""
        if not self._properties and os.path.isfile(self._properties_filename):
            with open(self._properties_filename, 'rb') as properties_file:
                self._properties = pickle.load(properties_file)
        return self._properties

    @property
    def mof_oms_df(self):
        """Get a pandas DataFrame that lists for each MOF whether it has an OMS or not and if it has an OMS what metal types it is.
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
            all_metal_species = ",".join(set(mp['metal_species']))
            mof_info[mp['name']] = {'Metal Types': all_metal_species,
                                    'Has OMS': has_oms,
                                    'OMS Types': oms_types}
        self._metal_site_df = pd.DataFrame.from_dict(mof_info,
                                                     orient='index')
        return self._metal_site_df

    @property
    def metal_site_df(self):
        """Get a pandas DataFrame that lists the OMS results for each metal type.
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
    def from_folder(cls, collection_folder, analysis_folder='analysis_folder', name_list=None):
        
        """Create a MofCollection from a the CIF files in a folder.

        Args:
            collection_folder: Path to the folder containing the CIF files to be added to the collection.
            analysis_folder: Path to the folder where the results will be stored.
            name_list: List of MOF names to include in the collection. If set, all the other CIF files in the folder will be excluded.

        Returns:
            -   A MofCollection object holding the specified MOF structures.
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

        Args:
            overwrite: Controls if the results will be overwritten or not.
            num_batches: Sets the number of batches the structures will be split in and analyzed on a separate process.
            analysis_limit: Analyze only up to the number of MOFs set by analysis_limit, if set to None all MOFs will be analyzed.
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
        """Iterate over all the MOFs in the collection and validate that they can be read and a MofStructure can be created.
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
        """Iterate over all the MOFs in the collection and check if the results from the OMS analysis exist.
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

        """Randomly select a sample of MOFs in the collection and return a new collection with the MOFs in the sample.

        Args:
            sample_size: Number of MOFs to be selected. Default value is 50.
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
        
        """Filter a collection given a number of filters. Calling this method of a MofCollection applies the filter and creates a new collection for the MOFs that match the filter. The cif files that match the filter are  copied to the new_collection_folder. The filters can be one or more of the following: 'density': [min, max] (range of values), 'oms_density': [min, max] (range of values), 'uc_volume':  [min, max] (range of values), 'metal_species': ["Cu", "Zn", ...] (list of metal species), 'non_metal_species': ["C", "N", ...] (list of non metal species), 'cif_okay': True (boolean value), 'has_oms': True (boolean value), 'mof_name':  [mof_name1, mof_name2] (string values)

        Args:
            using_filter: Filter used to identify MOFs with certain characteristics. Has to be a python dictionary.
            new_collection_folder: Path to the folder where the CIF files of the filtered collection will be stored. If set to None the CIF files will not be copied.
            new_analysis_folder: Path to the folder where the OMS result files of the filtered collection will be stored. If set to None the result files will not be copied.

        Returns:
            -   A MofCollection with only the filtered MOFs. If new_collection_folder or new_analysis_folder is not set then the collection will point to the original location of these files.
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
        """Iterate over all MOF files in the collection, load each CIF and store MOF properties such as density, unit cell volume etc.
        """
        print(self.separator)
        print('Reading CIF files and updating properties...')
        self._loop_over_collection(self._update_property_from_cif_file)
        self._store_properties()
        print('Done')
        print(self.separator)

    def read_oms_results(self):
        """Iterate over all MOF files in the collection, load each OMS result file and store OMS information to the MOF properties.
        """
        print(self.separator)
        print('Adding results to properties.')
        self._loop_over_collection(self._update_property_from_oms_result)
        print('Done')
        self._store_properties()
        print(self.separator)

    def copy_cifs(self, target_folder):

        """Copy cif files from their existing location to the specified target_folder.

        Args:
            target_folder: Path of folder to copy collection CIF files to.
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

        """Copy OMS result files from their existing location to the specified target_folder.

        Args:
            target_folder: Path of folder to copy collection OMS result files to.
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

        """Create a summary table for the OMS results of the collection, group results by metal type.

        Args:
            max_atomic_number: Maximum atomic number to be included in summary table. If not defined all metal atoms will be considered.
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
        """Summarize the t-factor information and make histograms for all the MOFs in the collection.
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
                    exit("MOF name and CIF checksum mismatch for {}.cif "
                         "{}.cif. Either the CIF files has already been "
                         "processed with a different name, or the CIF file "
                         "has changed since it was processed."
                         "".format(mof_name,
                                   self.properties[checksum]['mof_name']))
            if self._check_if_results_exist(mof_name):
                self._compare_checksums(mof_file, mof_name, checksum)
        print("\nAll Done.")
        self._store_properties()

    def _compare_checksums(self, mof_file, mof_name, checksum):
        """If OMS results exist for one of the CIF names in the collection then ensure that the CIF checksum matches the one in the result file.
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
        """For a given CIF file, create MofStructure object and run OMS analysis. If overwrite is false check if results already exist first.
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

        Args:
            num_batches: Number of batches (default: 1)
            overwrite: Controls if the results will be overwritten or not (default: False)
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

        """Iterate over all the MOFs in the collection and run the specified function.
        Args:
            func: Function to use.
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
        """Check if a given property can be found in the properties dictionary. If not try to read the CIF file and check again. If the check fails again try to read the OMS results and check again. If the check fails a third time return False, the property cannot be validated."""
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
        """Summarize the findings in table form and histograms for a give t-factor.
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
        s_yes['t_factor'].hist(bins=50, range=(0, 1), normed=False)
        s_no = s_n[s_n['is_open_yn'] == 'no']
        s_no['t_factor'].hist(bins=50, range=(0, 1), normed=False)
        plt.show()

    @staticmethod
    def _write_histogram(sites, dens, target):
        """Generate histograms to be used for summarizing the t-factor results.
        """
        hist, edges = np.histogram(sites, bins=50, range=(0, 1), density=dens)
        with open(target, 'w') as hist_file:
            w = (edges[1] - edges[0]) / 2
            for e, h in zip(edges, hist):
                print(e + w, h, file=hist_file)

    @staticmethod
    def _group_and_summarize(df, names=None):
        """Group the DataFrame holding the OMS results by metal type and rename its columns.
        """
        rename = {"mof_name": names[0], "is_open": names[1]}
        agg_dict = {"mof_name": pd.Series.nunique, "is_open": "count"}
        return df.groupby('metal').agg(agg_dict).rename(columns=rename)