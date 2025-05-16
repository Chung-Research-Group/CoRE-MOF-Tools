# -*- coding: utf-8 -*-

"""The methods to featurize a dataset of porous materials for heat capacity."""

import string
from pathlib import Path
import numpy as np
import pandas as pd
from typing import List
from pymatgen.io.cif import CifParser
from matminer.featurizers.site import GaussianSymmFunc, SiteElementalProperty,AGNIFingerprints 
from matminer.featurizers.base import BaseFeaturizer
from matminer.utils.data import MagpieData
from pymatgen.analysis.local_env import VoronoiNN


def featurize_dataset(cifs: list, verbos=False, saveto: str="features.csv")-> pd.DataFrame:

    """Featurize crystal structures using elemetal, geometric, and chemical descriptors for local environments.

    Args:
        cifs: list of paths to crystal structure in cif format.
        verbos: printing the steps.
        saveto: filename to save the generated features.
        
    """

    features={}
    for cif in cifs:
        structure = CifParser(cif).get_structures()[0]
        structure_name = Path(cif).name
        features[structure_name]={}
        features[structure_name]["structure"]=structure
        features[structure_name]["structure_name"]=Path(cif).name
        features[structure_name]["structure_path"]=str(Path(cif).parent)


    data=pd.DataFrame.from_dict(features).T
    features_dict={}
    ## 1. initialize the dictionary for each site
    for index,row in data.iterrows():
        structure=row["structure"]
        for atomidx in range(structure.num_sites):
            site_name="%s_%i"%(index,atomidx)
            features_dict[site_name]={}
            features_dict[site_name]={"structure_name": row["structure_name"]}
            features_dict[site_name].update({"structure_path": row["structure_path"]})

    ## 2. Site Elemental Property
    print("site elemental properties")
    property_list=("Number","AtomicWeight","Row","Column","Electronegativity","CovalentRadius")
    SEP = SiteElementalProperty(properties=property_list)
    colnames=SEP._generate_column_labels(multiindex=False,return_errors=False)
    for index,row in data.iterrows():
        structure=row["structure"]
        if verbos:
            print(index)
        for atomidx in range(structure.num_sites):
            feat=SEP.featurize(structure,idx=atomidx)
            site_name="%s_%i"%(index,atomidx)
            features_dict[site_name].update(dict(zip(colnames, feat)))

    ## 3. AGNI
    print("AGNI")
    property_list=("Number","AtomicWeight","Row","Column","Electronegativity","CovalentRadius")
    AGNI = AGNIFingerprints(cutoff=5,directions=[None])
    colnames=AGNI._generate_column_labels(multiindex=False,return_errors=False)
    for index,row in data.iterrows():
        structure=row["structure"]
        if verbos:
            print(index)
        for atomidx in range(structure.num_sites):
            feat=AGNI.featurize(structure,idx=atomidx)
            site_name="%s_%i"%(index,atomidx)
            features_dict[site_name].update(dict(zip(colnames, feat)))


    ## 4. Gaussian Symmetry Functions 
    print("GSF")
    GSF = GaussianSymmFunc(cutoff=5)
    colnames=GSF._generate_column_labels(multiindex=False,return_errors=False)
    for index,row in data.iterrows():
        structure=row["structure"]
        if verbos:
            print(index)
        for atomidx in range(structure.num_sites):
            feat=GSF.featurize(structure,idx=atomidx)
            site_name="%s_%i"%(index,atomidx)
            features_dict[site_name].update(dict(zip(colnames, feat)))

    ## 5. site difference stats 
    print("LPD")
    LPD = LocalPropertyStatsNew(properties=property_list)
    colnames=LPD._generate_column_labels(multiindex=False,return_errors=False)
    for index,row in data.iterrows():
        structure=row["structure"]
        if verbos:
            print(index)
        for atomidx in range(structure.num_sites):
            feat=LPD.featurize(structure,idx=atomidx)
            site_name="%s_%i"%(index,atomidx)
            features_dict[site_name].update(dict(zip(colnames, feat)))


    df_features=pd.DataFrame.from_dict(features_dict).T
    if saveto:
        df_features.to_csv(saveto)

    return df_features


def featurize_structure(cif: string, verbos=False, saveto: str="features.csv")-> pd.DataFrame:
    """Featurize a crystal structure using elemetal, geometric, and chemical descriptors for local environments.

    Args:
        cifs: list of paths to crystal structure in cif format.
        verbos: printing the steps.
        saveto: filename to save the generated features.
    """

    structure = CifParser(cif).get_structures()[0]
    structure_name = Path(cif).name
    features = {structure_name:{}}
    features[structure_name]["structure"]=structure
    features[structure_name]["structure_name"]=Path(cif).name
    features[structure_name]["structure_path"]=str(Path(cif).parent)
    
    data=pd.DataFrame.from_dict(features).T
    features_dict={}
    ## 1. initialize the dictionary for each site
    for index,row in data.iterrows():
        structure=row["structure"]
        for atomidx in range(structure.num_sites):
            site_name="%s_%i"%(index,atomidx)
            features_dict[site_name]={"structure_name": row["structure_name"]}
            features_dict[site_name].update({"structure_path": row["structure_path"]})

    ## 2. Site Elemental Property
    property_list=("Number","AtomicWeight","Row","Column","Electronegativity","CovalentRadius")
    SEP = SiteElementalProperty(properties=property_list)
    colnames=SEP._generate_column_labels(multiindex=False,return_errors=False)
    for index,row in data.iterrows():
        structure=row["structure"]
        if verbos:
            print(index)
        for atomidx in range(structure.num_sites):
            feat=SEP.featurize(structure,idx=atomidx)
            site_name="%s_%i"%(index,atomidx)
            features_dict[site_name].update(dict(zip(colnames, feat)))

    ## 3. AGNI
    property_list=("Number","AtomicWeight","Row","Column","Electronegativity","CovalentRadius")
    AGNI = AGNIFingerprints(cutoff=5,directions=[None])
    colnames=AGNI._generate_column_labels(multiindex=False,return_errors=False)
    for index,row in data.iterrows():
        structure=row["structure"]
        if verbos:
            print(index)
        for atomidx in range(structure.num_sites):
            feat=AGNI.featurize(structure,idx=atomidx)
            site_name="%s_%i"%(index,atomidx)
            features_dict[site_name].update(dict(zip(colnames, feat)))


    ## 4. Gaussian Symmetry Functions 
    GSF = GaussianSymmFunc(cutoff=5)
    colnames=GSF._generate_column_labels(multiindex=False,return_errors=False)
    for index,row in data.iterrows():
        structure=row["structure"]
        if verbos:
            print(index)
        for atomidx in range(structure.num_sites):
            feat=GSF.featurize(structure,idx=atomidx)
            site_name="%s_%i"%(index,atomidx)
            features_dict[site_name].update(dict(zip(colnames, feat)))

    ## 5. site difference stats 
    LPD = LocalPropertyStatsNew(properties=property_list)
    colnames=LPD._generate_column_labels(multiindex=False,return_errors=False)
    for index,row in data.iterrows():
        structure=row["structure"]
        if verbos:
            print(index)
        for atomidx in range(structure.num_sites):
            feat=LPD.featurize(structure,idx=atomidx)
            site_name="%s_%i"%(index,atomidx)
            features_dict[site_name].update(dict(zip(colnames, feat)))


    df_features=pd.DataFrame.from_dict(features_dict).T
    if saveto:
        df_features.to_csv(saveto)

    return df_features


class LocalPropertyStatsNew(BaseFeaturizer):

    """
    Differences, minima and maxima in elemental properties between site and its neighboring sites.
    Uses the Voronoi tessellation of the structure to determine the
    neighbors of the site, and assigns each neighbor (:math:`n`) a
    weight (:math:`A_n`) that corresponds to the area of the facet
    on the tessellation corresponding to that neighbor.
    The local property difference is then computed by
    :math:`\\frac{\sum_n {A_n |p_n - p_0|}}{\sum_n {A_n}}`
    where :math:`p_n` is the property (e.g., atomic number) of a neighbor
    and :math:`p_0` is the property of a site. If signed parameter is assigned
    True, signed difference of the properties is returned instead of absolute
    difference. Taken from oximachine featurizer: https://github.com/kjappelbaum/oximachine_featurizer                          
                                             
    Features:
        -   "local property stat in [property]"                               
    References:
        -   `Ward et al. _PRB_ 2017 <http://link.aps.org/doi/10.1103/PhysRevB.96.024104>`_
    """

    def __init__(
        self,
        data_source=MagpieData(),
        weight: str = "area",
        properties: List[str] = ("Electronegativity",),
        cutoff: List[str] = 5,
    ):
        """Initialize the featurizer.

        Args:
            data_source (AbstractData): Class from which to retrieve elemental properties.
            weight (str): What aspect of each voronoi facet to use to weigh each neighbor (see VoronoiNN).
            properties (List[str]): List of properties to use.
            cutoff (float): cut off value.
        """
        self.data_source = data_source
        self.properties = properties
        self.weight = weight
        self.cutoff = cutoff

    @staticmethod
    def from_preset(preset: str, cutoff: float = 13):
        """
        Create a new LocalPropertyStats class according to a preset.

        Args:
            preset (str): Name of preset
            cutoff (float): Cutoff for the nearest neighbor search
        """

        if preset == "interpretable":
            return LocalPropertyStatsNew(
                data_source=MagpieData(),
                properties=[
                    "MendeleevNumber",
                    "Column",
                    "Row",
                    "Electronegativity",
                    "NsValence",
                    "NpValence",
                    "NdValence",
                    "NfValence",
                    "NValence",
                    "NsUnfilled",
                    "NpUnfilled",
                    "NdUnfilled",
                    "NfUnfilled",
                    "NUnfilled",
                    "GSbandgap",
                ],
                cutoff=cutoff,
            )
        else:
            raise ValueError("Unrecognized preset: " + preset)

    def featurize(self, strc, idx: int):
        # Get the targeted site
        my_site = strc[idx]

        # Get the tessellation of a site
        nn = VoronoiNN(
            weight=self.weight,
            tol=0.0,
            cutoff=self.cutoff,
            compute_adj_neighbors=False,
        ).get_nn_info(strc, idx)

        # Get the element and weight of each site
        elems = [n["site"].specie for n in nn]
        weights = [n["weight"] for n in nn]

        # Compute the difference for each property
        output = np.zeros((len(self.properties),))
        output_signed = np.zeros((len(self.properties),))
        output_max = np.zeros((len(self.properties),))
        output_min = np.zeros((len(self.properties),))

        total_weight = np.sum(weights)
        for i, p in enumerate(self.properties):
            my_prop = self.data_source.get_elemental_property(my_site.specie, p)
            n_props = self.data_source.get_elemental_properties(elems, p)
            output[i] = (
                np.dot(weights, np.abs(np.subtract(n_props, my_prop))) / total_weight
            )
            output_signed[i] = (
                np.dot(weights, np.subtract(n_props, my_prop)) / total_weight
            )
            output_max[i] = np.max(np.subtract(n_props, my_prop))
            output_min[i] = np.min(np.subtract(n_props, my_prop))
        return np.hstack([output, output_signed, output_max, output_min])

    def feature_labels(self):

        return (
            ["local difference in " + p for p in self.properties]
            + ["local signed difference in " + p for p in self.properties]
            + ["maximum local difference in " + p for p in self.properties]
            + ["minimum local difference in " + p for p in self.properties]
        )

    def citations(self):
        return [
            "@article{Ward2017,"
            "author = {Ward, Logan and Liu, Ruoqian "
            "and Krishna, Amar and Hegde, Vinay I. "
            "and Agrawal, Ankit and Choudhary, Alok "
            "and Wolverton, Chris},"
            "doi = {10.1103/PhysRevB.96.024104},"
            "journal = {Physical Review B},"
            "pages = {024104},"
            "title = {{Including crystal structure attributes "
            "in machine learning models of formation energies "
            "via Voronoi tessellations}},"
            "url = {http://link.aps.org/doi/10.1103/PhysRevB.96.014107},"
            "volume = {96},year = {2017}}",
            "@article{jong_chen_notestine_persson_ceder_jain_asta_gamst_2016,"
            "title={A Statistical Learning Framework for Materials Science: "
            "Application to Elastic Moduli of k-nary Inorganic Polycrystalline Compounds}, "
            "volume={6}, DOI={10.1038/srep34256}, number={1}, journal={Scientific Reports}, "
            "author={Jong, Maarten De and Chen, Wei and Notestine, Randy and Persson, "
            "Kristin and Ceder, Gerbrand and Jain, Anubhav and Asta, Mark and Gamst, Anthony}, "
            "year={2016}, month={Mar}}",
        ]

    def implementors(self):
        return ["Logan Ward", "Aik Rui Tan"]


def add_type_feature(mydict,atomtype,name,feat):
    if atomtype in mydict:
        mydict[atomtype][name]=feat
    else:
        mydict[atomtype]={name:feat}
    return mydict
