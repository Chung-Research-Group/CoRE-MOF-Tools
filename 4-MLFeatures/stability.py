import pandas as pd
import numpy as np
import os,keras
import keras.backend as K
import pickle as pkl
import joblib
import cloudpickle

os.environ["CUDA_VISIBLE_DEVICES"] = "-1" 

def precision(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
    precision = true_positives / (predicted_positives + K.epsilon())
    return precision
def recall(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    recall = true_positives / (possible_positives + K.epsilon())
    return recall
def f1(y_true, y_pred):
    p = precision(y_true, y_pred)
    r = recall(y_true, y_pred)
    return 2 * ((p * r) / (p + r + K.epsilon()))

solvent_feature_names = ['Df', 'Di', 'Dif', 'GPOAV_1_86', 'GPONAV_1_86', 'GPOV_1_86',
                            'GSA_1_86', 'POAV_1_86', 'POAV_vol_frac_1_86', 'PONAV_1_86', 'PONAV_vol_frac_1_86',
                            'VPOV_1_86', 'VSA_1_86', 'unit_cell_volume_1_86', 'f-chi-0-all', 'f-chi-1-all',
                            'f-chi-2-all', 'f-chi-3-all', 'f-Z-0-all', 'f-Z-1-all', 'f-Z-2-all', 'f-Z-3-all',
                            'f-I-0-all', 'f-I-1-all', 'f-I-2-all', 'f-I-3-all', 'f-T-0-all', 'f-T-1-all',
                            'f-T-2-all', 'f-T-3-all', 'f-S-0-all', 'f-S-1-all', 'f-S-2-all', 'f-S-3-all',
                            'mc-chi-0-all', 'mc-chi-1-all', 'mc-chi-2-all', 'mc-chi-3-all', 'mc-Z-0-all',
                            'mc-Z-1-all', 'mc-Z-2-all', 'mc-Z-3-all', 'mc-I-1-all', 'mc-I-2-all', 'mc-I-3-all',
                            'mc-T-0-all', 'mc-T-1-all', 'mc-T-2-all', 'mc-T-3-all', 'mc-S-0-all', 'mc-S-1-all',
                            'mc-S-2-all', 'mc-S-3-all', 'D_mc-chi-1-all', 'D_mc-chi-2-all', 'D_mc-chi-3-all', 
                            'D_mc-Z-1-all', 'D_mc-Z-2-all', 'D_mc-Z-3-all', 'D_mc-T-1-all', 'D_mc-T-2-all',
                            'D_mc-T-3-all', 'D_mc-S-1-all', 'D_mc-S-2-all', 'D_mc-S-3-all', 'f-lig-chi-0',
                            'f-lig-chi-1', 'f-lig-chi-2', 'f-lig-chi-3', 'f-lig-Z-0', 'f-lig-Z-1', 'f-lig-Z-2',
                            'f-lig-Z-3', 'f-lig-I-0', 'f-lig-I-1', 'f-lig-I-2', 'f-lig-I-3', 'f-lig-T-0', 'f-lig-T-1',
                            'f-lig-T-2', 'f-lig-T-3', 'f-lig-S-0', 'f-lig-S-1', 'f-lig-S-2', 'f-lig-S-3', 'lc-chi-0-all',
                            'lc-chi-1-all', 'lc-chi-2-all', 'lc-chi-3-all', 'lc-Z-0-all', 'lc-Z-1-all', 'lc-Z-2-all', 'lc-Z-3-all',
                            'lc-I-1-all', 'lc-I-2-all', 'lc-I-3-all', 'lc-T-0-all', 'lc-T-1-all', 'lc-T-2-all', 'lc-T-3-all',
                            'lc-S-0-all', 'lc-S-1-all', 'lc-S-2-all', 'lc-S-3-all', 'D_lc-chi-1-all', 'D_lc-chi-2-all',
                            'D_lc-chi-3-all', 'D_lc-Z-1-all', 'D_lc-Z-2-all', 'D_lc-Z-3-all', 'D_lc-T-1-all', 'D_lc-T-2-all',
                            'D_lc-T-3-all', 'D_lc-S-1-all', 'D_lc-S-2-all', 'D_lc-S-3-all', 'func-chi-0-all', 'func-chi-1-all',
                            'func-chi-2-all', 'func-chi-3-all', 'func-Z-0-all', 'func-Z-1-all', 'func-Z-2-all', 'func-Z-3-all',
                            'func-I-0-all', 'func-I-1-all', 'func-I-2-all', 'func-I-3-all', 'func-T-0-all', 'func-T-1-all',
                            'func-T-2-all', 'func-T-3-all', 'func-S-0-all', 'func-S-1-all', 'func-S-2-all', 'func-S-3-all',
                            'D_func-chi-1-all', 'D_func-chi-2-all', 'D_func-chi-3-all', 'D_func-Z-1-all', 'D_func-Z-2-all',
                            'D_func-Z-3-all', 'D_func-T-1-all', 'D_func-T-2-all', 'D_func-T-3-all', 'D_func-S-1-all',
                            'D_func-S-2-all', 'D_func-S-3-all']

thermal_feature_names = ['f-chi-0-all', 'f-chi-1-all', 'f-chi-2-all', 'f-chi-3-all',
                        'f-Z-0-all', 'f-Z-1-all', 'f-Z-2-all', 'f-Z-3-all', 'f-I-0-all', 'f-I-1-all', 'f-I-2-all',
                        'f-I-3-all', 'f-T-0-all', 'f-T-1-all', 'f-T-2-all', 'f-T-3-all', 'f-S-0-all', 'f-S-1-all',
                        'f-S-2-all', 'f-S-3-all', 'mc-chi-0-all', 'mc-chi-1-all', 'mc-chi-2-all', 'mc-chi-3-all',
                        'mc-Z-0-all', 'mc-Z-1-all', 'mc-Z-2-all', 'mc-Z-3-all', 'mc-I-1-all', 'mc-I-2-all',
                        'mc-I-3-all', 'mc-T-0-all', 'mc-T-1-all', 'mc-T-2-all', 'mc-T-3-all', 'mc-S-0-all',
                        'mc-S-1-all', 'mc-S-2-all', 'mc-S-3-all', 'D_mc-chi-1-all', 'D_mc-chi-2-all',
                        'D_mc-chi-3-all', 'D_mc-Z-1-all', 'D_mc-Z-2-all', 'D_mc-Z-3-all', 'D_mc-T-1-all',
                        'D_mc-T-2-all', 'D_mc-T-3-all', 'D_mc-S-1-all', 'D_mc-S-2-all', 'D_mc-S-3-all', 'f-lig-chi-0',
                        'f-lig-chi-1', 'f-lig-chi-2', 'f-lig-chi-3', 'f-lig-Z-0', 'f-lig-Z-1', 'f-lig-Z-2', 'f-lig-Z-3',
                        'f-lig-I-0', 'f-lig-I-1', 'f-lig-I-2', 'f-lig-I-3', 'f-lig-T-0', 'f-lig-T-1', 'f-lig-T-2', 'f-lig-T-3',
                        'f-lig-S-0', 'f-lig-S-1', 'f-lig-S-2', 'f-lig-S-3', 'lc-chi-0-all', 'lc-chi-1-all', 'lc-chi-2-all',
                        'lc-chi-3-all', 'lc-Z-0-all', 'lc-Z-1-all', 'lc-Z-2-all', 'lc-Z-3-all', 'lc-I-1-all', 'lc-I-2-all',
                        'lc-I-3-all', 'lc-T-0-all', 'lc-T-1-all', 'lc-T-2-all', 'lc-T-3-all', 'lc-S-0-all', 'lc-S-1-all', 
                        'lc-S-2-all', 'lc-S-3-all', 'D_lc-chi-1-all', 'D_lc-chi-2-all', 'D_lc-chi-3-all', 'D_lc-Z-1-all',
                        'D_lc-Z-2-all', 'D_lc-Z-3-all', 'D_lc-T-1-all', 'D_lc-T-2-all', 'D_lc-T-3-all', 'D_lc-S-1-all',
                        'D_lc-S-2-all', 'D_lc-S-3-all', 'func-chi-0-all', 'func-chi-1-all', 'func-chi-2-all', 'func-chi-3-all',
                        'func-Z-0-all', 'func-Z-1-all', 'func-Z-2-all', 'func-Z-3-all', 'func-I-0-all', 'func-I-1-all', 'func-I-2-all',
                        'func-I-3-all', 'func-T-0-all', 'func-T-1-all', 'func-T-2-all', 'func-T-3-all', 'func-S-0-all', 'func-S-1-all',
                        'func-S-2-all', 'func-S-3-all', 'D_func-chi-1-all', 'D_func-chi-2-all', 'D_func-chi-3-all', 'D_func-Z-1-all',
                        'D_func-Z-2-all', 'D_func-Z-3-all', 'D_func-T-1-all', 'D_func-T-2-all', 'D_func-T-3-all', 'D_func-S-1-all',
                        'D_func-S-2-all', 'D_func-S-3-all', 'Df', 'Di', 'Dif', 'GPOAV_1_86', 'GPONAV_1_86', 'GPOV_1_86', 'GSA_1_86',
                        'POAV_1_86', 'POAV_vol_frac_1_86', 'PONAV_1_86', 'PONAV_vol_frac_1_86', 'VPOV_1_86', 'VSA_1_86',
                        'unit_cell_volume_1_86']

rfa_2_class_water_features = ['mc-Z-3-all', 'D_mc-Z-3-all', 'D_mc-Z-2-all',
                                'D_mc-Z-1-all', 'mc-chi-3-all', 'mc-Z-1-all',
                                'mc-Z-0-all', 'D_mc-chi-2-all', 'f-lig-Z-2',
                                'GSA_1_4', 'f-lig-I-0', 'func-S-1-all']


class run:
    def __init__(self, zeo_input, rac_input, output_folder):
        self.zeo_input = zeo_input
        self.rac_input = rac_input
        self.output = output_folder
        self.csv_path = os.path.join(output_folder, "stability_log.csv")
        self.process()

    def process(self):
        current_path = os.getcwd()
        dependencies = {'precision':precision,'recall':recall,'f1':f1}
        solvent_model = keras.models.load_model(current_path+'/features/ML/models/final_model_flag_few_epochs.h5', custom_objects=dependencies)
        thermal_model = keras.models.load_model(current_path+'/features/ML/models/final_model_T_few_epochs.h5', custom_objects=dependencies)
        with open(current_path+'/features/ML/models/solvent_scaler.pkl', 'rb') as f:
            solvent_scaler = pkl.load(f)
        with open(current_path+'/features/ML/models/thermal_x_scaler.pkl', 'rb') as f:
            thermal_x_scaler = pkl.load(f)
        with open(current_path+'/features/ML/models/thermal_y_scaler.pkl', 'rb') as f:
            thermal_y_scaler = pkl.load(f)
        with open(current_path+'/features/ML/models/water_model.pkl', 'rb') as f:
            water_model = cloudpickle.load(f)
        with open(current_path+'/features/ML/models/water_scaler.pkl', 'rb') as f:
            water_scaler = pkl.load(f)

        os.makedirs(self.output, exist_ok=True)

        df_zeo = pd.read_csv(self.zeo_input)
        df_rac = pd.read_csv(self.rac_input)
        df =  pd.merge(df_zeo, df_rac, on="name", how="inner") 
        MOF_names = df['name'].tolist()

        X_solvent = df[solvent_feature_names].to_numpy()
        X_solvent = solvent_scaler.transform(X_solvent)

        solvent_model_prob = solvent_model.predict(X_solvent)
        solvent_model_prob = solvent_model_prob.flatten()
        solvent_model_label = solvent_model_prob > 0.5
        solvent_model_label = solvent_model_label.astype('int')

        df2 = pd.DataFrame({'name': MOF_names, 
                            'solvent_pred': solvent_model_prob,
                            'solvent_label': solvent_model_label,
                            })
        df2.to_csv(self.output +'/activation_stability.csv', index=False)

        X_thermal = df[thermal_feature_names].to_numpy()
        X_thermal = thermal_x_scaler.transform(X_thermal)

        thermal_model_pred = thermal_y_scaler.inverse_transform(thermal_model.predict(X_thermal))
        thermal_model_pred = np.round(thermal_model_pred, 1)
        thermal_model_pred = thermal_model_pred.flatten()

        df3 = pd.DataFrame({'name': MOF_names, 
                            'thermal_pred': thermal_model_pred,
                            })
        df3.to_csv(self.output + '/thermal_stability.csv', index=False)

        X_water = df[rfa_2_class_water_features].to_numpy()
        X_water = water_scaler.transform(X_water)

        water_model_prob = water_model.predict_proba(X_water)[:,1]
        water_model_label = water_model.predict(X_water)

        df4 = pd.DataFrame({'name': MOF_names, 
            'water_pred': water_model_prob,
            'water_label': water_model_label,
            })
        df4.to_csv(self.output + '/water_stability.csv', index=False)