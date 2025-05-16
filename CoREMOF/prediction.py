"""ML-predicted features.
"""

from PACMANCharge import pmcharge

import os, warnings, shutil, joblib, requests, cloudpickle
import keras
import keras.backend as K
import pandas as pd

warnings.filterwarnings('ignore')

import numpy as np
import pickle as pkl

from CoREMOF.calculation import Zeopp
from CoREMOF.calculation.mof_features import RACs, Volume

package_directory = os.path.abspath(__file__).replace("prediction.py","")

def get_files_from_github(repo, path):

    """query files from github due to limit of uploading size by PyPi.

    Args:
        repo (str): github repository.
        path (str): the path of models.

    Returns:
        Dictionary:
            -   response of downloading.
    """
        
    url = f"https://api.github.com/repos/{repo}/contents/{path}"
    headers = {'Accept': 'application/vnd.github.v3+json'}
    response = requests.get(url, headers=headers)
    response.raise_for_status() 
    return response.json()

def download_file(url, save_path):

    """download models from github due to limit of uploading size by PyPi.

    Args:
        url (str): link of downloading file.
        save_path (str): the path to save files.

    """    

    if not os.path.exists(save_path):
        response = requests.get(url)
        response.raise_for_status()
        os.makedirs(os.path.dirname(save_path), exist_ok=True) 
        with open(save_path, 'wb') as file:
            file.write(response.content)
        print(f"Downloaded {save_path}")
    else:
        pass

repo = 'sxm13/CoREMOF_tools'
github_paths = [
    'models/cp_app/ensemble_models_smallML_120_100/300',
    'models/cp_app/ensemble_models_smallML_120_100/350',
    'models/cp_app/ensemble_models_smallML_120_100/400',
    'models/stability'
]
local_directories = [
    '/models/cp_app/ensemble_models_smallML_120_100/300',
    '/models/cp_app/ensemble_models_smallML_120_100/350',
    '/models/cp_app/ensemble_models_smallML_120_100/400',
    '/models/stability'
]

path_map = dict(zip(github_paths, local_directories))

for github_path in github_paths:
    files_info = get_files_from_github(repo, github_path)
    for file_info in files_info:
        if file_info['type'] == 'file':
            raw_url = file_info['download_url']
            local_path_suffix = path_map[github_path].lstrip('/')
            local_file_path = os.path.join(package_directory, local_path_suffix, os.path.basename(file_info['path']))
            download_file(raw_url, local_file_path)

from CoREMOF.models.cp_app.descriptors import cv_features
from CoREMOF.models.cp_app.featurizer import featurize_structure
from CoREMOF.models.cp_app.predictions import predict_Cv_ensemble_structure_multitemperatures


def pacman(structure, output_folder="result_pacman", charge_type="DDEC6", digits=10, atom_type=True, neutral=True, keep_connect=False):

    """predict partial atom charge by PACMAN Charge: https://doi.org/10.1021/acs.jctc.4c00434.

    Args:
        structure (str): path to your structure.
        output_folder (str): the path to save CIF with predicted charges.
        charge_type (str): models of DDEC6, Bader, CM5 or REPEAT.
        digits (int): number of decimal places to print for partial atomic charges. ML models were trained on a 6-digit dataset.
        atom_type (bool): keep the same partial atomic charge for the same atom types (based on the similarity of partial atomic charges up to 2 decimal places).
        neutral (bool): keep the net charge is zero. We use "mean" method to neuralize the system where the excess charges are equally distributed across all atoms.
        keep_connect (bool): retain the atomic and connection information (such as _atom_site_adp_type, bond) for the structure.

    Returns:
        Dictionary & cif:
            -   predicted PBE energy and bandgap of your structure.
            -   CIF with predicted charges.
    """  
        
    results_eb = {}
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder, exist_ok=True)
    
    name = os.path.basename(structure).replace(".cif", "")
    results_eb["Name"] = name

    try:
        pmcharge.predict(
            cif_file=structure,
            charge_type=charge_type,
            digits=digits,
            atom_type=atom_type,
            neutral=neutral,
            keep_connect=keep_connect
        )
        
        pbe, bandgap = pmcharge.Energy(cif_file=structure)
        
        results_eb["PBE Energy"] = pbe
        results_eb["Bandgap"] = bandgap
      
        shutil.move(structure.replace(".cif","_pacman.cif"),
                    os.path.join(output_folder,
                                structure.split("/")[-1].replace(".cif","_pacman.cif")))

        return results_eb
    
    except Exception as e:
        print(e)

def cp(structure, T=[300, 350, 400]):
    
    """predict heat capacity by GBR models at different temperatures: https://doi.org/10.1038/s41563-022-01374-3.

    Args:
        structure (str): path to your structure.
        T (list): the temperatures of your system.

    Returns:
        Dictionary:
            -   unit by ["unit"], always "J/g/K", "J/mol/K".
            -   predicted heat capacity of your structure.    
    """  
        
    name = os.path.basename(structure).replace(".cif", "")
    featurize_structure(structure, verbos=False, saveto="features.csv")
    
    predict_Cv_ensemble_structure_multitemperatures(
                                                    path_to_models=package_directory+"models/cp_app/ensemble_models_smallML_120_100",
                                                    structure_name=name + ".cif",
                                                    features_file="features.csv", 
                                                    FEATURES=cv_features,
                                                    temperatures=T,
                                                    save_to="cp.csv"
                                                    )
    result_ = pd.read_csv("cp.csv")
    result_cp = {}
    result_cp["unit"] = "J/g/K", "J/mol/K"

    for t in T:
        result_cp[str(t)+"_mean"] = [result_["Cv_gravimetric_"+str(t)+"_mean"].iloc[0],
                                        result_["Cv_molar_"+str(t)+"_mean"].iloc[0]]
        result_cp[str(t)+"_std"] = [result_["Cv_gravimetric_"+str(t)+"_std"].iloc[0],
                                        result_["Cv_molar_"+str(t)+"_std"].iloc[0]]

    os.remove("features.csv")
    os.remove("cp.csv")

    return result_cp


def precision(y_true, y_pred):

    """
    Computes the precision metric for binary classification.

    Precision is defined as the ratio of correctly predicted positive observations 
    to the total predicted positive observations. It is given by:

        precision = TP / (TP + FP)

    where:
        - TP (True Positives)  : Correctly predicted positive cases
        - FP (False Positives)  : Incorrectly predicted positive cases

    Args:
        y_true (tensor): Ground truth binary labels (0 or 1).
        y_pred (tensor): Predicted probabilities or binary predictions.

    Returns:
        tensor: Precision score (between 0 and 1).
    """
        
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
    precision = true_positives / (predicted_positives + K.epsilon())
    return precision


def recall(y_true, y_pred):

    """
    Computes the recall metric for binary classification.

    Recall (also known as Sensitivity or True Positive Rate) measures the proportion 
    of actual positives that are correctly identified by the model. It is given by:

        recall = TP / (TP + FN)

    where:
        - TP (True Positives)  : Correctly predicted positive cases
        - FN (False Negatives)  : Actual positive cases that were predicted as negative

    Args:
        y_true (tensor): Ground truth binary labels (0 or 1).
        y_pred (tensor): Predicted probabilities or binary predictions.

    Returns:
        tensor: Recall score (between 0 and 1).
    """

    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    recall = true_positives / (possible_positives + K.epsilon())
    return recall


def f1(y_true, y_pred):

    """
    Computes the F1-score for binary classification.

    The F1-score is the harmonic mean of precision and recall, providing a balanced 
    metric when the dataset is imbalanced. It is given by:

        F1 = 2 * (precision * recall) / (precision + recall)

    where:
        - Precision = TP / (TP + FP)
        - Recall    = TP / (TP + FN)

    Args:
        y_true (tensor): Ground truth binary labels (0 or 1).
        y_pred (tensor): Predicted probabilities or binary predictions.

    Returns:
        tensor: F1-score (between 0 and 1).
    """
        
    p = precision(y_true, y_pred)
    r = recall(y_true, y_pred)
    return 2 * ((p * r) / (p + r + K.epsilon()))


import tensorflow as tf
# tf.config.set_visible_devices([], 'GPU')

# os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
def stability(structure):
    
    """predict stability of MOFs: https://doi.org/10.1021/jacs.1c07217, https://doi.org/10.1021/jacs.4c05879.

    Args:
        structure (str): path to your structure.

    Returns:
        Dictionary:
            -   unit by ["unit"], always "nan, °C, nan".
            -   predicted thermal, solvent and water stabilities.     
    """  
        
    solvent_feature_names = [
                            'f-chi-0-all', 'f-chi-1-all', 'f-chi-2-all', 'f-chi-3-all', 'f-Z-0-all',
                            'f-Z-1-all', 'f-Z-2-all', 'f-Z-3-all', 'f-I-0-all', 'f-I-1-all', 'f-I-2-all',
                            'f-I-3-all', 'f-T-0-all', 'f-T-1-all', 'f-T-2-all', 'f-T-3-all', 'f-S-0-all',
                            'f-S-1-all', 'f-S-2-all', 'f-S-3-all', 'mc-chi-0-all', 'mc-chi-1-all',
                            'mc-chi-2-all', 'mc-chi-3-all', 'mc-Z-0-all', 'mc-Z-1-all', 'mc-Z-2-all',
                            'mc-Z-3-all', 'mc-I-1-all', 'mc-I-2-all', 'mc-I-3-all', 'mc-T-0-all',
                            'mc-T-1-all', 'mc-T-2-all', 'mc-T-3-all', 'mc-S-0-all', 'mc-S-1-all',
                            'mc-S-2-all', 'mc-S-3-all', 'D_mc-chi-1-all', 'D_mc-chi-2-all',
                            'D_mc-chi-3-all', 'D_mc-Z-1-all', 'D_mc-Z-2-all', 'D_mc-Z-3-all',
                            'D_mc-T-1-all', 'D_mc-T-2-all', 'D_mc-T-3-all', 'D_mc-S-1-all',
                            'D_mc-S-2-all', 'D_mc-S-3-all', 'f-lig-chi-0', 'f-lig-chi-1',
                            'f-lig-chi-2', 'f-lig-chi-3', 'f-lig-Z-0', 'f-lig-Z-1', 'f-lig-Z-2',
                            'f-lig-Z-3', 'f-lig-I-0', 'f-lig-I-1', 'f-lig-I-2', 'f-lig-I-3',
                            'f-lig-T-0', 'f-lig-T-1', 'f-lig-T-2', 'f-lig-T-3', 'f-lig-S-0',
                            'f-lig-S-1', 'f-lig-S-2', 'f-lig-S-3', 'lc-chi-0-all', 'lc-chi-1-all',
                            'lc-chi-2-all', 'lc-chi-3-all', 'lc-Z-0-all', 'lc-Z-1-all', 'lc-Z-2-all',
                            'lc-Z-3-all', 'lc-I-1-all', 'lc-I-2-all', 'lc-I-3-all', 'lc-T-0-all',
                            'lc-T-1-all', 'lc-T-2-all', 'lc-T-3-all', 'lc-S-0-all', 'lc-S-1-all',
                            'lc-S-2-all', 'lc-S-3-all', 'D_lc-chi-1-all', 'D_lc-chi-2-all',
                            'D_lc-chi-3-all', 'D_lc-Z-1-all', 'D_lc-Z-2-all', 'D_lc-Z-3-all',
                            'D_lc-T-1-all', 'D_lc-T-2-all', 'D_lc-T-3-all', 'D_lc-S-1-all',
                            'D_lc-S-2-all', 'D_lc-S-3-all', 'func-chi-0-all', 'func-chi-1-all',
                            'func-chi-2-all', 'func-chi-3-all', 'func-Z-0-all', 'func-Z-1-all',
                            'func-Z-2-all', 'func-Z-3-all', 'func-I-0-all', 'func-I-1-all',
                            'func-I-2-all', 'func-I-3-all', 'func-T-0-all', 'func-T-1-all',
                            'func-T-2-all', 'func-T-3-all', 'func-S-0-all', 'func-S-1-all',
                            'func-S-2-all', 'func-S-3-all', 'D_func-chi-1-all', 'D_func-chi-2-all',
                            'D_func-chi-3-all', 'D_func-Z-1-all', 'D_func-Z-2-all', 'D_func-Z-3-all',
                            'D_func-T-1-all', 'D_func-T-2-all', 'D_func-T-3-all', 'D_func-S-1-all',
                            'D_func-S-2-all', 'D_func-S-3-all']

    thermal_feature_names = ['f-chi-0-all', 'f-chi-1-all', 'f-chi-2-all', 'f-chi-3-all',
                            'f-Z-0-all', 'f-Z-1-all', 'f-Z-2-all', 'f-Z-3-all', 'f-I-0-all',
                            'f-I-1-all', 'f-I-2-all', 'f-I-3-all', 'f-T-0-all', 'f-T-1-all',
                            'f-T-2-all', 'f-T-3-all', 'f-S-0-all', 'f-S-1-all', 'f-S-2-all',
                            'f-S-3-all', 'mc-chi-0-all', 'mc-chi-1-all', 'mc-chi-2-all',
                            'mc-chi-3-all', 'mc-Z-0-all', 'mc-Z-1-all', 'mc-Z-2-all',
                            'mc-Z-3-all', 'mc-I-1-all', 'mc-I-2-all', 'mc-I-3-all',
                            'mc-T-0-all', 'mc-T-1-all', 'mc-T-2-all', 'mc-T-3-all',
                            'mc-S-0-all', 'mc-S-1-all', 'mc-S-2-all', 'mc-S-3-all',
                            'D_mc-chi-1-all', 'D_mc-chi-2-all', 'D_mc-chi-3-all',
                            'D_mc-Z-1-all', 'D_mc-Z-2-all', 'D_mc-Z-3-all', 'D_mc-T-1-all',
                            'D_mc-T-2-all', 'D_mc-T-3-all', 'D_mc-S-1-all', 'D_mc-S-2-all',
                            'D_mc-S-3-all', 'f-lig-chi-0', 'f-lig-chi-1', 'f-lig-chi-2',
                            'f-lig-chi-3', 'f-lig-Z-0', 'f-lig-Z-1', 'f-lig-Z-2', 'f-lig-Z-3',
                            'f-lig-I-0', 'f-lig-I-1', 'f-lig-I-2', 'f-lig-I-3', 'f-lig-T-0',
                            'f-lig-T-1', 'f-lig-T-2', 'f-lig-T-3', 'f-lig-S-0', 'f-lig-S-1',
                            'f-lig-S-2', 'f-lig-S-3', 'lc-chi-0-all', 'lc-chi-1-all', 'lc-chi-2-all',
                            'lc-chi-3-all', 'lc-Z-0-all', 'lc-Z-1-all', 'lc-Z-2-all', 'lc-Z-3-all',
                            'lc-I-1-all', 'lc-I-2-all', 'lc-I-3-all', 'lc-T-0-all', 'lc-T-1-all',
                            'lc-T-2-all', 'lc-T-3-all', 'lc-S-0-all', 'lc-S-1-all', 'lc-S-2-all',
                            'lc-S-3-all', 'D_lc-chi-1-all', 'D_lc-chi-2-all', 'D_lc-chi-3-all',
                            'D_lc-Z-1-all', 'D_lc-Z-2-all', 'D_lc-Z-3-all', 'D_lc-T-1-all',
                            'D_lc-T-2-all', 'D_lc-T-3-all', 'D_lc-S-1-all', 'D_lc-S-2-all',
                            'D_lc-S-3-all', 'func-chi-0-all', 'func-chi-1-all', 'func-chi-2-all',
                            'func-chi-3-all', 'func-Z-0-all', 'func-Z-1-all', 'func-Z-2-all',
                            'func-Z-3-all', 'func-I-0-all', 'func-I-1-all', 'func-I-2-all',
                            'func-I-3-all', 'func-T-0-all', 'func-T-1-all', 'func-T-2-all',
                            'func-T-3-all', 'func-S-0-all', 'func-S-1-all', 'func-S-2-all',
                            'func-S-3-all', 'D_func-chi-1-all', 'D_func-chi-2-all', 'D_func-chi-3-all',
                            'D_func-Z-1-all', 'D_func-Z-2-all', 'D_func-Z-3-all', 'D_func-T-1-all',
                            'D_func-T-2-all', 'D_func-T-3-all', 'D_func-S-1-all', 'D_func-S-2-all',
                            'D_func-S-3-all']

    water_feature_names = ['mc-Z-3-all', 'D_mc-Z-3-all', 'D_mc-Z-2-all',
                            'D_mc-Z-1-all', 'mc-chi-3-all', 'mc-Z-1-all',
                            'mc-Z-0-all', 'D_mc-chi-2-all', 'f-lig-Z-2',
                            "ASA",
                            'f-lig-I-0', 'func-S-1-all']

    result_stability = {}
    result_stability["unit"] = "nan, °C, nan"

    dependencies = {'precision':precision,'recall':recall,'f1':f1}
    solvent_model = keras.models.load_model(package_directory+'/models/stability/final_model_flag_few_epochs.h5', custom_objects=dependencies)
    thermal_model = keras.models.load_model(package_directory+'/models/stability/final_model_T_few_epochs.h5', custom_objects=dependencies)
    with open(package_directory+'/models/stability/solvent_scaler.pkl', 'rb') as f:
        solvent_scaler = pkl.load(f)
    with open(package_directory+'/models/stability/thermal_x_scaler.pkl', 'rb') as f:
        thermal_x_scaler = pkl.load(f)
    with open(package_directory+'/models/stability/thermal_y_scaler.pkl', 'rb') as f:
        thermal_y_scaler = pkl.load(f)
    with open(package_directory+'/models/stability/water_model.pkl', 'rb') as f:
        water_model = cloudpickle.load(f)
    with open(package_directory+'/models/stability/water_scaler.pkl', 'rb') as f:
        water_scaler = pkl.load(f)

    results_pd = Zeopp.PoreDiameter(structure)
    results_pv = Zeopp.PoreVolume(structure, 1.86, 1.86, 10000, True)
    results_sa = Zeopp.SurfaceArea(structure, 1.86, 1.86, 10000, True)
    result_v = Volume(structure)
    results_sa_1_4 = Zeopp.SurfaceArea(structure, 1.4, 1.4, 10000, True)
    result_RACs = RACs(structure)
    
    '''
        zeopp4sol in the start
        LCD,
        PLD,
        LFPD,
        GPOAV_1_86,
        GPONAV_1_86,
        GPOV_1_86,
        GSA_1_86,
        POAV_1_86,
        POAV_vol_frac_1_86,
        PONAV_1_86,
        PONAV_vol_frac_1_86,
        VPOV_1_86,
        VSA_1_86,
        cell_volume
    '''

    X_solvent = [results_pd["LCD"],
                    results_pd["PLD"],
                    results_pd["LFPD"],
                    results_pv["PV"][1],
                    results_pv["NPV"][1],
                    results_pv["PV"][1] + results_pv["NPV"][1],
                    results_sa["ASA"][2],
                    results_pv["PV"][0],
                    results_pv["VF"],
                    results_pv["NPV"][0],
                    results_pv["NVF"],
                    results_pv["VF"] + results_pv["NVF"],
                    results_sa["ASA"][1],
                    result_v["total_volume"]
                ]
    
    for fn_sol in solvent_feature_names:
        try:
            X_solvent.append(result_RACs["Metal"][fn_sol])
        except:
            try:
                X_solvent.append(result_RACs["Linker"][fn_sol])
            except:
                X_solvent.append(result_RACs["Function-group"][fn_sol])

    # X_solvent = np.array(X_solvent).reshape(-1, 1).flatten().reshape(1, -1).tolist()
    
    X_solvent = solvent_scaler.transform([X_solvent])
    solvent_model_prob = solvent_model.predict(X_solvent)
    solvent_model_prob = solvent_model_prob.flatten()
    result_stability["solvent removal probability"] = float(solvent_model_prob[0])

    '''
        zeopp4therm in the end
        LCD,
        PLD,
        LFPD,
        GPOAV_1_86,
        GPONAV_1_86,
        GPOV_1_86,
        GSA_1_86,
        POAV_1_86,
        POAV_vol_frac_1_86,
        PONAV_1_86,
        PONAV_vol_frac_1_86,
        VPOV_1_86,
        VSA_1_86,
        cell_volume
    '''

    X_thermal = []
    for fn_them in thermal_feature_names:
        try:
            X_thermal.append(result_RACs["Metal"][fn_them])
        except:
            try:
                X_thermal.append(result_RACs["Linker"][fn_them])
            except:
                X_thermal.append(result_RACs["Function-group"][fn_them])

    for zeo_them in [results_pd["LCD"],
                    results_pd["PLD"],
                    results_pd["LFPD"],
                    results_pv["PV"][1],
                    results_pv["NPV"][1],
                    results_pv["PV"][1] + results_pv["NPV"][1],
                    results_sa["ASA"][2],
                    results_pv["PV"][0],
                    results_pv["VF"],
                    results_pv["NPV"][0],
                    results_pv["NVF"],
                    results_pv["VF"] + results_pv["NVF"],
                    results_sa["ASA"][1],
                    result_v["total_volume"]
                ]:
        
        X_thermal.append(zeo_them)
    X_thermal = thermal_x_scaler.transform([X_thermal])
    thermal_model_pred = thermal_y_scaler.inverse_transform(thermal_model.predict(X_thermal))
    thermal_model_pred = np.round(thermal_model_pred, 1)
    thermal_model_pred = thermal_model_pred.flatten()
    result_stability["thermal stability"] = float(thermal_model_pred[0])

    X_water = []
    for fn_water in water_feature_names:
        try:
            X_water.append(result_RACs["Metal"][fn_water])
        except:
            try:
                X_water.append(result_RACs["Linker"][fn_water])
            except:
                try:
                    X_water.append(result_RACs["Function-group"][fn_water])
                except:
                    X_water.append(results_sa_1_4[fn_water][2])

    X_water = water_scaler.transform([X_water])

    water_model_prob = water_model.predict_proba(X_water)[:,1]
    # water_model_label = water_model.predict(X_water)
    result_stability["water probability"] = float(water_model_prob[0])

    return result_stability