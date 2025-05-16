import numpy as np
import pandas as pd
from .descriptors import cv_features
import joblib
import glob
import copy

FEATURES = cv_features

def predict_Cv_ensemble_structure(models: list, FEATURES: list, df_features: pd.DataFrame, temperature: float) -> list:
     
    """Predict heat capacity using an ensemble of ML models for one structure.

    Args:
        models: list (ensemble) of ML models.
        FEATURES: features for ML model.
        df_features: pandas dataframe containing the features.
        temperature: target temperature.

    Returns:
        a list containing the gravimetric and molar heat capacity together with the uncertainty of the models.
    """
     
    if len(df_features["structure_name"].unique())>2:
        raise ValueError("More than one structure in the features file...")
        
    df_site_structure = copy.deepcopy(df_features)
    structure_name = df_site_structure["structure_name"].unique()[0]
    print(structure_name)
    for model_idx,model in enumerate(models):
        df_site_structure["pCv_{}_predicted_{}".format(temperature, model_idx)]=model.predict(df_site_structure[FEATURES])
    results=[]
    predicted_mol=[]
    predicted_gr=[]
    for model_idx in range(len(models)):
        sites=df_site_structure.loc[df_site_structure["structure_name"]==structure_name]
        predicted_mol.append(np.sum(sites["pCv_{}_predicted_{}".format(temperature,model_idx)])/len(sites))
        predicted_gr.append(np.sum(sites["pCv_{}_predicted_{}".format(temperature,model_idx)])/np.sum(sites["site AtomicWeight"]))
    results.append({
        "name": structure_name,
        "Cv_gravimetric_{}_mean".format(temperature): np.mean(predicted_gr),
        "Cv_gravimetric_{}_std".format(temperature): np.std(predicted_gr),
        "Cv_molar_{}_mean".format(temperature): np.mean(predicted_mol),
        "Cv_molar_{}_std".format(temperature): np.std(predicted_mol),
    })
    return results


def predict_Cv_ensemble_structure_multitemperatures(path_to_models: str, structure_name: str, features_file: str="features.csv", 
                FEATURES: list=cv_features, temperatures: list=[300.00], save_to: str="cv_predicted.csv") -> pd.DataFrame:

    """Predict heat capacity using an ensemble of ML models for a dataset.

    Args:
        models: ensemble of ML models.
        FEATURES: features for ML model.
        df_features: pandas dataframe containing the features.
        temperature: target temperature.

    Returns:
        a list containing the gravimetric and molar heat capacity together with the uncertainty of the models.
    """
     
    df_features = pd.read_csv(features_file)
    df_features["structure_name"]=["_".join(n.split("_")[:-1]) for n in df_features["Unnamed: 0"]]
    df_features=df_features.loc[df_features["structure_name"]==structure_name]
    for i,temperature in enumerate(temperatures):
        models=[]
        print("loading models for:", temperature)
        modelnames = glob.glob("{}/{:.0f}/*".format(path_to_models, temperature))
        models = [joblib.load(n) for n in modelnames]
        print("{} models loaded, predicting...".format(len(models)))
        if i==0:
            res= pd.DataFrame(predict_Cv_ensemble_structure(models, FEATURES, df_features, temperature))
            all_results=res
        else:
            res= pd.DataFrame(predict_Cv_ensemble_structure(models, FEATURES, df_features, temperature))
            all_results=all_results.merge(res, how="inner",on="name")

    if save_to:
        all_results.to_csv(save_to)
    return all_results