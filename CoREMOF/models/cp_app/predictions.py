import numpy as np
import pandas as pd
from .descriptors import cv_features
import joblib
import glob
import copy
from CoREMOF._heat_capacity import normalize_temperatures

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
     
    if not models:
        raise ValueError(f"No heat-capacity models were supplied for {temperature} K")
    if df_features.empty:
        raise ValueError("No atomic features were found for the requested structure")
    if "structure_name" not in df_features:
        raise ValueError("The feature table must contain a 'structure_name' column")
    if df_features.columns.has_duplicates:
        raise ValueError("The feature table contains duplicate column names")
    if df_features["structure_name"].isna().any():
        raise ValueError("Every atomic feature row must have a structure_name")
    structure_names = df_features["structure_name"].unique()
    if len(structure_names) != 1:
        raise ValueError(
            "Expected features for exactly one structure; found "
            f"{len(structure_names)}: {structure_names.tolist()}"
        )
    missing_features = [name for name in FEATURES if name not in df_features.columns]
    if missing_features:
        preview = ", ".join(missing_features[:5])
        suffix = " ..." if len(missing_features) > 5 else ""
        raise ValueError(f"Feature table is missing required columns: {preview}{suffix}")
    if "site AtomicWeight" not in df_features:
        raise ValueError("Feature table is missing required column 'site AtomicWeight'")
    if not FEATURES or len(set(FEATURES)) != len(FEATURES):
        raise ValueError("Feature names must be nonempty and distinct")
    try:
        values = df_features[FEATURES].to_numpy(dtype=float)
        weights = df_features["site AtomicWeight"].to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("Atomic features and weights must be numeric") from exc
    if not np.isfinite(values).all():
        raise ValueError("Atomic features must all be finite; no imputation is performed")
    if not np.isfinite(weights).all() or (weights <= 0).any():
        raise ValueError("Atomic weights must all be finite and positive")
    if not np.isfinite(np.sum(weights)):
        raise ValueError("Total atomic weight must be finite")
        
    df_site_structure = copy.deepcopy(df_features)
    structure_name = structure_names[0]
    for model_idx,model in enumerate(models):
        # Preserve the model's floating-point dtype and historical summation
        # arithmetic. Promoting XGBoost float32 outputs changes valid values.
        predicted = np.asarray(model.predict(df_site_structure[FEATURES]))
        if (predicted.shape != (len(df_site_structure),)
                or predicted.dtype.kind not in "biuf"
                or not np.isfinite(predicted).all()):
            raise ValueError(
                f"Heat-capacity model {model_idx} must return one finite value per atom"
            )
        if not np.isfinite(np.sum(predicted)):
            raise ValueError(f"Heat-capacity model {model_idx} has a non-finite total")
        df_site_structure["pCv_{}_predicted_{}".format(temperature, model_idx)] = predicted
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
    if not all(np.isfinite(value) for key, value in results[0].items() if key != "name"):
        raise ValueError("Heat-capacity ensemble aggregation produced non-finite values")
    return results


def predict_Cv_ensemble_structure_multitemperatures(
    path_to_models: str,
    structure_name: str,
    features_file: str = "features.csv",
    FEATURES: list = cv_features,
    temperatures=None,
    save_to: str = "cv_predicted.csv",
    verbose: bool = False,
) -> pd.DataFrame:

    """Predict heat capacity using an ensemble of ML models for a dataset.

    Args:
        path_to_models: directory containing one subdirectory per temperature.
        structure_name: CIF filename stored in the feature table.
        features_file: CSV produced by :func:`featurize_structure`.
        FEATURES: feature columns expected by each model.
        temperatures: temperatures in kelvin. Defaults to ``[300.0]``.
        save_to: optional output CSV path. Use ``None`` to skip writing.
        verbose: print model-loading progress.

    Returns:
        One-row dataframe containing gravimetric and molar heat capacities and
        ensemble standard deviations.
    """
     
    if temperatures is None:
        temperatures = [300.0]
    temperatures = normalize_temperatures(temperatures)

    all_features = pd.read_csv(features_file)
    if "structure_name" not in all_features.columns:
        if "Unnamed: 0" not in all_features.columns:
            raise ValueError(
                "Feature CSV must contain 'structure_name' or the legacy 'Unnamed: 0' column"
            )
        all_features["structure_name"] = [
            "_".join(str(name).rsplit("_", 1)[:-1])
            for name in all_features["Unnamed: 0"]
        ]
    df_features = all_features.loc[all_features["structure_name"] == structure_name]
    if df_features.empty:
        available = sorted(all_features["structure_name"].dropna().unique())
        preview = ", ".join(map(str, available[:5])) or "none"
        raise ValueError(
            f"Structure {structure_name!r} was not found in {features_file}. "
            f"Available structures: {preview}"
        )

    for i,temperature in enumerate(temperatures):
        if verbose:
            print("loading models for:", temperature)
        modelnames = sorted(glob.glob("{}/{:.0f}/*".format(path_to_models, temperature)))
        if not modelnames:
            raise FileNotFoundError(
                f"No ensemble models found for {temperature:g} K under {path_to_models}"
            )
        models = [joblib.load(n) for n in modelnames]
        if verbose:
            print("{} models loaded, predicting...".format(len(models)))
        if i==0:
            res= pd.DataFrame(predict_Cv_ensemble_structure(models, FEATURES, df_features, temperature))
            all_results=res
        else:
            res= pd.DataFrame(predict_Cv_ensemble_structure(models, FEATURES, df_features, temperature))
            all_results=all_results.merge(res, how="inner",on="name")

    if save_to:
        all_results.to_csv(save_to, index=False)
    return all_results
