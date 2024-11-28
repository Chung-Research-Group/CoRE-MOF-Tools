#!/usr/bin/env python
# -*- coding: utf-8 -*-
import flask
import tensorflow as tf
import pandas as pd
import numpy as np
import time
import json
import os
import shutil
import subprocess
import stat
import keras
import keras.backend as K
import sklearn
import smtplib
import glob
import uuid
import pickle as pkl
import joblib 
from sklearn.metrics import pairwise_distances
from bokeh.plotting import figure
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.palettes import Inferno256
from flask import request, session
# import flask_login
# from flask_login import LoginManager, UserMixin, login_required, current_user
from molSimplify.Informatics.MOF.MOF_descriptors import get_primitive, get_MOF_descriptors
from flask_cors import CORS
from datetime import datetime
from pymongo import MongoClient
from werkzeug.utils import secure_filename
# from list_content.list_content import my_linkers, my_sbus, my_nets, my_MOFs

cmap_bokeh = Inferno256

MOFSIMPLIFY_PATH = os.path.abspath('.') # the main directory
MOFSIMPLIFY_PATH += '/'
USE_SPLASH_PAGE = False

app = flask.Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 1024 * 1024 * 20 # Upload max 20 megabytes
app.config['UPLOAD_EXTENSIONS'] = ['.jpg', '.jpeg', '.png', '.pdf', '.tiff', '.tif', '.eps'] # acceptable file types
    # Note: these are a superset of the extensions indicated on the form and allowed on the front end, so some of these extensions don't actually apply
# app.secret_key = str(json.load(open('secret_key.json','r'))['key']) # secret key
# cors = CORS(app)

operation_counter = 0 # This variable keeps track of server traffic to alert the user if they should wait until later to make a request. 
    # I only change operation_counter for the more time-consuming operations (the two predictions, and component analysis).
    # operation_counter is periodically set back to zero, since if a user closes their browser in the middle of an operation operation_counter is not properly -=1'd
last_operation_counter_clear = time.time() # the current time when server is started

MAX_OPERATIONS = 4 # This variable dictates the maximum number of concurrent operations, to prevent server overload.

# The following three functions are needed for the ANN models.
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

# global variables for the ANN models, loaded once when server first started
# solvent_ANN_path = MOFSIMPLIFY_PATH + 'model/solvent/ANN/' 
thermal_ANN_path = MOFSIMPLIFY_PATH + 'model/thermal/ANN/'
dependencies = {'precision':precision,'recall':recall,'f1':f1}
# loading the ANN models to save time later
tf_session = tf.Session()
from tensorflow import keras as tf_keras
tf_keras.backend.set_session(tf_session)
##### Note: the h5 models for the solvent removal stability prediction ANN and the thermal stability prediction ANN should be based on the same version of TensorFlow (here, 1.14). #####
# solvent_model = keras.models.load_model(solvent_ANN_path + 'final_model_flag_few_epochs.h5',custom_objects=dependencies)
thermal_model = keras.models.load_model(thermal_ANN_path + 'final_model_T_few_epochs.h5',custom_objects=dependencies)

# # These two models use scikit-learn 0.23.2
# water_model = joblib.load('model/water_and_acid/models/water_model.joblib')
# acid_model = joblib.load('model/water_and_acid/models/acid_model.joblib')
# water_scaler = joblib.load('model/water_and_acid/models/water_scaler.joblib')
# acid_scaler = joblib.load('model/water_and_acid/models/acid_scaler.joblib')

# ethane_ANN_path = MOFSIMPLIFY_PATH + 'model/c2/ethane/'
# ethylene_ANN_path = MOFSIMPLIFY_PATH + 'model/c2/ethylene/'

# ethane_model = keras.models.load_model(ethane_ANN_path + 'ethane_model.h5')
# ethylene_model = keras.models.load_model(ethylene_ANN_path + 'ethylene_model.h5')

# global variable dictionary for stable building block aliases. See https://zenodo.org/record/7091192
bb_mapping = {'XUDNAN_clean_linker_2': 'E0', 'MIFKUJ_clean_linker_0': 'E1', 'AJUNOK_clean_linker_0': 'E2', 'OJICUG_clean_linker_0': 'E3',
 'FUNLAD_clean_linker_1': 'E4', 'CIFDUS_clean_linker_0': 'E5', 'EYOCIG_clean_linker_4': 'E6', 'TEMPAE_clean_linker_10': 'E7',
 'KIFKEQ_clean_linker_1': 'E8', 'FUNLAD_clean_linker_0': 'E9', 'EZIPEK_clean_linker_0': 'E10', 'UKIBIB_clean_linker_0': 'E11',
 'UKIBUN_clean_linker_0': 'E12', 'UVAHIK_clean_linker_0': 'E13', 'XIVTED_clean_linker_0': 'E14', 'ZUJKAS_clean_linker_0': 'E15',
 'APAZEY_clean_linker_0': 'orgN0', 'AZAVOO_clean_linker_0': 'orgN1', 'BELTOD_clean_linker_0': 'orgN2', 'BETDAH_clean_linker_0': 'orgN3',
 'BETDEL_clean_linker_0': 'orgN4', 'BETDIP_clean_linker_0': 'orgN5', 'CEKHIL_clean_linker_0': 'orgN6', 'CUKXOW_clean_linker_0': 'orgN7',
 'CUQRIR_clean_linker_1': 'orgN8', 'ESEHIV_clean_linker_0': 'orgN9', 'GIZVER_clean_linker_0': 'orgN10', 'INOVEN_clean_linker_0': 'orgN11',
 'IYUCEM_clean_linker_1': 'orgN12', 'IZENUY_clean_linker_0': 'orgN13', 'KAKCAD_clean_linker_0': 'orgN14', 'KUFVIS_clean_linker_0': 'orgN15',
 'MAKGUD_clean_linker_0': 'orgN16', 'MIDCAF_clean_linker_0': 'orgN17', 'MIFKUJ_clean_linker_1': 'orgN18', 'MUDLON_clean_linker_0': 'orgN19',
 'NAHDIM_clean_linker_0': 'orgN20', 'NAWXER_clean_linker_0': 'orgN21', 'PAMTOO_clean_linker_0': 'orgN22', 'QEWDON_clean_linker_0': 'orgN23',
 'QUQFOY_clean_linker_2': 'orgN24', 'TATPOV_clean_linker_0': 'orgN25', 'VETTIZ_charged_linker_0': 'orgN26', 'WOSHET_clean_linker_1': 'orgN27',
 'BEQFEK_clean_linker_11': 'orgN28', 'LEJCEK_clean_linker_0': 'orgN29', 'UFIREI_clean_linker_0': 'orgN30', 'UWAGAB01_clean_linker_0': 'orgN31',
 'BEQFEK_clean_sbu_0': 'N0', 'BOTCEU_clean_sbu_0': 'N1', 'ENOWUB_clean_sbu_0': 'N2', 'FANWIC_clean_sbu_0': 'N3',
 'ICAMEG_clean_sbu_1': 'N4', 'LIZSOE_clean_sbu_0': 'N5', 'UKALOJ_clean_sbu_0': 'N6', 'UKIBUN_clean_sbu_0': 'N7',
 'ZALLEG_clean_sbu_0': 'N8', 'AJUNOK_clean_sbu_0': 'N9', 'AZAVOO_clean_sbu_0': 'N10', 'BELTOD_clean_sbu_0': 'N11',
 'BETDAH_clean_sbu_0': 'N12', 'BETDEL_clean_sbu_0': 'N13', 'BETDIP_clean_sbu_0': 'N14', 'BETFAJ_clean_sbu_0': 'N15',
 'BETFEN_clean_sbu_0': 'N16', 'BETGAK_clean_sbu_0': 'N17', 'CEKHIL_clean_sbu_0': 'N18', 'CIFDUS_clean_sbu_1': 'N19',
 'CUKXOW_clean_sbu_0': 'N20', 'CUWYAW_clean_sbu_0': 'N21', 'EBIMEJ_clean_sbu_0': 'N22', 'EQERAU_clean_sbu_0': 'N23',
 'ESEHIV_clean_sbu_0': 'N24', 'EYACOX_clean_sbu_0': 'N25', 'EYACOX_clean_sbu_1': 'N26', 'EZIPEK_clean_sbu_0': 'N27',
 'FUNLAD_clean_sbu_0': 'N28', 'GALJAG_clean_sbu_0': 'N29', 'GEDQOX_clean_sbu_0': 'N30', 'GIZVER_clean_sbu_0': 'N31',
 'HICVOG_clean_sbu_0': 'N32', 'HISSIN_clean_sbu_0': 'N33', 'HISSIN_clean_sbu_1': 'N34', 'ICIZOL_clean_sbu_0': 'N35',
 'INOVEN_clean_sbu_0': 'N36', 'IYUCEM_clean_sbu_0': 'N37', 'IZENUY_clean_sbu_0': 'N38', 'JUFBIX_clean_sbu_0': 'N39',
 'KAKCAD_clean_sbu_0': 'N40', 'KOZSID_clean_sbu_0': 'N41', 'KUFVIS_clean_sbu_0': 'N42', 'KUMBOL_clean_sbu_0': 'N43',
 'KUMBOL_clean_sbu_1': 'N44', 'KUMBUR_clean_sbu_0': 'N45', 'KUMBUR_clean_sbu_1': 'N46', 'KUMJIN_clean_sbu_0': 'N47',
 'KUMJIN_clean_sbu_1': 'N48', 'LEVNOQ01_clean_sbu_1': 'N49', 'MAKGOX_clean_sbu_0': 'N50', 'MAKGUD_clean_sbu_1': 'N51',
 'MIDCAF_clean_sbu_0': 'N52', 'MIFKUJ_clean_sbu_0': 'N53', 'MUDLON_clean_sbu_0': 'N54', 'NAHDIM_clean_sbu_0': 'N55',
 'NAHDIM_clean_sbu_2': 'N56', 'NAWXER_clean_sbu_0': 'N57', 'NUHQIS_clean_sbu_0': 'N58', 'NUHQUE_clean_sbu_2': 'N59',
 'NUHRAL_clean_sbu_3': 'N60', 'OJICUG_clean_sbu_0': 'N61', 'OLANAS_clean_sbu_1': 'N62', 'OLANAS_clean_sbu_2': 'N63',
 'OLANEW_clean_sbu_0': 'N64', 'OLANEW_clean_sbu_1': 'N65', 'OLANEW_clean_sbu_2': 'N66', 'OLANEW_clean_sbu_3': 'N67',
 'OLANEW_clean_sbu_4': 'N68', 'PAMTOO_clean_sbu_0': 'N69', 'PAMTUU_clean_sbu_0': 'N70', 'PORLAL_clean_sbu_0': 'N71',
 'QEWDON_clean_sbu_0': 'N72', 'QUQFOY_clean_sbu_1': 'N73', 'QUQFOY_clean_sbu_3': 'N74', 'SARMOO_clean_sbu_0': 'N75',
 'TAGTUT_clean_sbu_1': 'N76', 'TAGTUT_clean_sbu_3': 'N77', 'TATPOV_clean_sbu_0': 'N78', 'UXUYUI_clean_sbu_0': 'N79',
 'VETTIZ_charged_sbu_0': 'N80', 'VOLQOD_clean_sbu_1': 'N81', 'WAQDOJ_charged_sbu_0': 'N82', 'WUSLED_clean_sbu_0': 'N83',
 'XADDAJ01_clean_sbu_3': 'N84', 'XOMCOT_clean_sbu_0': 'N85', 'XOMCOT_clean_sbu_1': 'N86', 'XUDNAN_clean_sbu_1': 'N87'}

# # Loading the names of MOFs from our hypothetical database. Constructed with ultrastable building blocks.
# with open('pickle_files/MOF_names/1inorganic_1edge_cifs_list.pkl', 'rb') as f:
#     list_1inorganic_1edge_MOFs = pkl.load(f)
# with open('pickle_files/MOF_names/1inorganic_1organic_1edge_cifs_list.pkl', 'rb') as f:
#     list_1inorganic_1organic_1edge_MOFs = pkl.load(f)
# with open('pickle_files/MOF_names/2inorganic_1edge_cifs_list.pkl', 'rb') as f:
#     list_2inorganic_1edge_MOFs = pkl.load(f)

# # Loading the names of the ultrastable MOFs from our hypothetical database.
# with open('pickle_files/ultrastable_MOFs/1inor_1edge_cifs_list_ultrastable.pkl', 'rb') as f:
#     list_1inorganic_1edge_MOFs_ultrastable = pkl.load(f)
# with open('pickle_files/ultrastable_MOFs/1inor_1org_1edge_cifs_list_ultrastable.pkl', 'rb') as f:
#     list_1inorganic_1organic_1edge_MOFs_ultrastable = pkl.load(f)
# with open('pickle_files/ultrastable_MOFs/2inor_1edge_cifs_list_ultrastable.pkl', 'rb') as f:
#     list_2inorganic_1edge_MOFs_ultrastable = pkl.load(f)

# ultrastable_MOFs = list_1inorganic_1edge_MOFs_ultrastable + list_1inorganic_1organic_1edge_MOFs_ultrastable + list_2inorganic_1edge_MOFs_ultrastable

def conditional_diminish(counter):
    """
    conditional_diminish decreases the input by one and returns it, unless the input is already zero.
    The input is intended to be operation_counter.
    operation_counter might be zero when conditional_diminish is called b/c of the periodic zero-ing of operation_counter. It is unlikely though.

    :param counter: An int, which should be operation_counter.
    :return: The new value of operation_counter.
    """ 
    if counter != 0:
        counter -= 1

    return counter

def operation_counter_periodic_clear():
    global operation_counter # global variable
    global last_operation_counter_clear

    if time.time() - last_operation_counter_clear > 300: # 5 minutes or more since last time operation_counter was zero'd
        last_operation_counter_clear = time.time()
        operation_counter = 0

# app.route takes jquery ($) requests from index.html and executes the associated function in app.py.
# Output can then be returned to index.html.

def set_ID():
    """
    set_ID sets the session user ID. 
    This is also used to generate unique folders, so that multiple users can use the website at a time. 
        The user's folder is temp_file_creation_[ID]
    Specifically, the function copies the temp_file_creation folder for the current user so that the rest of MOFSimplify functionality can be executed for the current user in their own folder.
        This lets multiple operations be run for different users concurrently.
    This function also deletes temp_file_creation copies from other users that have not been used for a while, in order to reduce clutter.

    :return: The session ID for this user.
    """ 

    session['ID'] = uuid.uuid4() # a unique ID for this session
    session['permission'] = True # keeps track of if user gave us permission to store the MOFs they predict on; defaults to Yes

    # make a version of the temp_file_creation folder for this user
    new_folder = MOFSIMPLIFY_PATH + '/temp_file_creation_' + str(session['ID'])
    shutil.copytree('temp_file_creation', new_folder)
    os.remove(new_folder + '/temp_cif.cif') # remove this, for sole purpose of updating time stamp on the new folder (copytree doesn't)

    # delete all temp_file_creation clone folders that haven't been used for a while, to prevent folder accumulation
    for root, dirs, files in os.walk(MOFSIMPLIFY_PATH):
        for dir in dirs:
            target_str = 'temp_file_creation'
            if len(dir) > len(target_str) and target_str in dir and file_age_in_seconds(dir) > 7200: # 7200s is two hours
                # target_str in dir since the names of all copies start with the sequence "temp_file_creation"
                # len(dir) > len(target_str) to prevent deleting the original temp_file_creation folder
                shutil.rmtree(dir)

    return str(session['ID']) # return a string

@app.route('/get_ID', methods=['GET'])
def get_ID():
    """
    get_ID gets the session user ID. 
    This is used for getting building block generated MOFs.

    :return: string, The session ID for this user.
    """ 
    return str(session['ID']) # return a string

@app.route('/permission', methods=['POST'])
def change_permission():
    """
    change_permission adjusts whether or not MOFSimplify stores information on the MOFs the user predicts on.

    If the user clicks "Yes" or "No" before get_lists() finishes running upon the website's startup in the browser (~4 seconds), their "Yes" or "No" input will not register for the session and will register for another session (the most recently started one) instead.

    :return: string, The boolean sent from the front end. We return this because we have to return something, but nothing is done with the returned value on the front end.
    """

    # Grab data
    permission = json.loads(flask.request.get_data())
    session['permission'] = permission
    print('Permission check')
    print(permission)

    return str(permission)

@app.route('/list_getter', methods=['GET'])
def get_lists():
    """
    get_lists gets the dropdown lists. 

    :return: dictionary. The dropdown lists
    """ 

    # Initializes the new user, since this function is called when the browser is first opened up.
    set_ID()

    return {'my_linkers':my_linkers, 'my_sbus':my_sbus, 'my_nets':my_nets, 'my_MOFs':my_MOFs}

# The send_from_directory functions that follow provide images from the MOFSimplify server to the website. The images are in a folder called images.
@app.route('/TGA_graphic.png')
def serve_TGA_graphic():
    return flask.send_from_directory('images', 'TGA_graphic.png')

@app.route('/banner_light')
def serve_banner_light():
    return flask.send_from_directory('images', 'MOF_light.webp') # Google's webp format. It is optimized for websites and loads quickly.

@app.route('/banner_dark')
def serve_banner_dark():
    return flask.send_from_directory('images', 'MOF_dark.webp') # Google's webp format. It is optimized for websites and loads quickly.

@app.route('/MOF_logo.png')
def serve_MOFSimplify_logo():
    return flask.send_from_directory('images', 'MOF_logo.png')

## Handle feedback
@app.route('/process_feedback', methods=['POST'])
def process_feedback():
    """
    process_feedback inserts MOFSimplify form feedback into the MongoDB feedback database. 
    If an uploaded file has an incorrect extension (i.e. is a disallowed file format), the user is directed to an error page.
    """ 
    client = MongoClient('18.18.63.68',27017) # connect to mongodb
        # The first argument is the IP address. The second argument is the port.
    db = client.feedback
    collection = db.MOFSimplify # The MOFSimplify collection in the feedback database.
    fields = ['feedback_form_name', 'rating', 'email', 'reason', 'comments', 'cif_file_name', 'structure', 'solvent']
    #$meta_fields = ['IP', 'datetime', 'cif_file', 'MOF_name']
    final_dict = {}
    for field in fields:
        final_dict[field] = request.form.get(field)

    # Populate special fields
    uploaded_file = request.files['file']
    if uploaded_file.filename == '' and request.form.get('feedback_form_name') != 'upload_form':
        # User did not upload the optional TGA trace
        print('No TGA trace')
    # if final_dict['file']==b'':
    #     file_ext = ''
    else:
        final_dict['filetype'] = uploaded_file.content_type
        filename = secure_filename(uploaded_file.filename)
        final_dict['filename'] = filename
        final_dict['file'] = uploaded_file.read()
        file_ext = os.path.splitext(filename)[1].lower()
        if file_ext not in app.config['UPLOAD_EXTENSIONS']:
            return ('', 204) # 204 no content response
            # return flask.send_from_directory('./splash_page/', 'error.html')

    # Special tasks if the form is upload_form
    if request.form.get('feedback_form_name') == 'upload_form':
        uploaded_cif = request.files['cif_file']
        cif_filename = secure_filename(uploaded_cif.filename)
        file_ext = os.path.splitext(cif_filename)[1].lower()
        if file_ext != '.cif':
            return ('', 204) # 204 no content response
            # return flask.send_from_directory('./splash_page/', 'error.html')
        final_dict['cif_file_name'] = cif_filename
        final_dict['structure'] = uploaded_cif.read()

    final_dict['ip'] = request.remote_addr
    final_dict['timestamp'] = datetime.now().isoformat()
    
    print(final_dict)
    collection.insert(final_dict) # insert the dictionary into the mongodb collection
    return ('', 204) # 204 no content response
    # return flask.send_from_directory('./splash_page/', 'success.html')

## Handle feedback, for the 2-class water and acid models
@app.route('/process_feedback_water', methods=['POST'])
def process_feedback_water():
    """
    process_feedback_water inserts MOFSimplify form feedback into the MongoDB feedback database. 
    """ 
    client = MongoClient('18.18.63.68',27017) # connect to mongodb
        # The first argument is the IP address. The second argument is the port.
    db = client.feedback
    collection = db.MOFSimplify_water # The MOFSimplify_water collection in the feedback database.
    fields = ['feedback_form_name', 'rating', 'email', 'reason', 'comments', 'cif_file_name', 'structure']
    final_dict = {}
    for field in fields:
        final_dict[field] = request.form.get(field)

    final_dict['ip'] = request.remote_addr
    final_dict['timestamp'] = datetime.now().isoformat()
    
    print(final_dict)
    collection.insert(final_dict) # insert the dictionary into the mongodb collection
    return ('', 204) # 204 no content response

## Handle removal request
@app.route('/process_removal', methods=['POST'])
def process_removal():
    """
    process_removal emails mofsimplify@mit.edu when the removal form is filled out.
    """ 

    email = request.form.get('email')
    comments = request.form.get('comments')
    ip = request.remote_addr
    timestamp = datetime.now().isoformat()

    # grabbing environment variables (should be in .bashrc_conda or .zshrc)
    EMAIL_ADDRESS = os.environ.get('EMAIL_USER')
    EMAIL_PASSWORD = os.environ.get('EMAIL_PASS')

    # https://www.youtube.com/watch?v=JRCJ6RtE3xU
    with smtplib.SMTP('smtp.gmail.com', 587) as smtp:
        smtp.ehlo()
        smtp.starttls()
        smtp.ehlo()

        # logging in to mail server
        smtp.login(EMAIL_ADDRESS, EMAIL_PASSWORD)

        subject = 'MOFSimplify removal request'
        body = f'email: {email}\ncomments: {comments}\nip: {ip}\ntimestamp: {timestamp}'

        msg = f'Subject: {subject}\n\n{body}'

        smtp.sendmail(EMAIL_ADDRESS, 'mofsimplify@mit.edu', msg) # sending to mofsimplify@mit.edu

    return ('', 204) # 204 no content response

## Splash page management. Splash page is currently disabled.
@app.route('/', methods=['GET', 'POST'])
@app.route('/<path:path>', methods=['GET', 'POST'])
def index(path='index.html'):
  if request.method == 'POST':
    username = 'user1'
    # user = User()
    # user.id = username
    # flask_login.login_user(user)
    # This is the section where the password is checked.
    # if request.form.get('password') == users[username]['password']:
    #   user = User()
    #   user.id = 'user1'
    #   flask_login.login_user(user)
  # print('is user authenticated?')
  # print(current_user.is_authenticated)
  # print('input check')
  # print(request.form.get('password'))
  # if current_user.is_authenticated:
  #   return flask.send_from_directory('.', 'index.html')
  # elif request.form.get('password') == None:
  #   return flask.send_from_directory('./splash_page/', path)
  # else:
  #   return flask.send_from_directory('./splash_page/', 'index_wrong_password.html')
  return flask.send_from_directory('.', 'index.html')

@app.route('/mof_examples/<path:path>') # needed for fetch
def serve_example(path):
    """
    serve_example returns a file to MOFSimplify.
    The file is intended to be a cif file of an example MOF. 
    So, this function serves the example MOF.

    :param path: The path to the desired example MOF in the mof_examples folder. For example, HKUST1.cif.
    :return: The MOF specified in the path input.
    """ 
    return flask.send_from_directory('mof_examples', path)

@app.route('/how_to_cite.html')
def serve_cite():
    """
    serve_cite serves the how to cite page.
    So the user is redirected to the how to cite page.

    :return: The how to cite page.
    """ 
    return flask.send_from_directory('.', 'how_to_cite.html')

@app.route('/water_stability_prediction.html')
def serve_water_stability_page():
    """
    serve_water_stability_page serves the water stability page.
    So the user is redirected to the water stability page.

    :return: The water stability page.
    """     
    return flask.send_from_directory('.', 'water_stability_prediction.html')

@app.route('/C2_uptake_prediction.html')
def serve_C2_uptake_page():
    """
    serve_C2_uptake_page serves the C2 uptake page.
    So the user is redirected to the C2 uptake page.

    :return: The C2 uptake page.
    """     
    return flask.send_from_directory('.', 'C2_uptake_prediction.html')

@app.route('/stable_MOFs.html')
def serve_stable_bb_page():
    """
    serve_stable_bb_page serves the stable MOFs page.
    So the user is redirected to the stable MOFs page.

    :return: The stable MOFs page.
    """ 
    return flask.send_from_directory('.', 'stable_MOFs.html')

@app.route('/libraries/<path:path>')
def serve_library_files(path):
    """
    serve_library_files returns a file to MOFSimplify.
    The file is intended to be a library file, either .js or .css.

    :param path: The path to the desired library in the libraries folder. For example, jquery-3.4.1.min.js.
    :return: The library specified in the path input.
    """ 
    return flask.send_from_directory('libraries', path)

@app.route('/list_content/<path:path>')
def serve_list_files(path):
    """
    serve_list_files returns a file to MOFSimplify.
    The file is intended to be a list file, which contains information for dropdowns.

    :param path: The path to the desired library in the libraries folder.
    :return: The list file specified in the path input.
    """ 
    return flask.send_from_directory('list_content', path)

@app.route('/bbcif/<path:path>')
def serve_bbcif(path):
    """
    serve_bbcif returns a file to MOFSimplify.
    The file is intended to be a cif file for a MOF that was constructed using MOFSimplify's building block functionality.
    So, this function serves the building block generated MOF.

    :param path: The path to the desired MOF in the user's building block folder.
    :return: The cif file for the building block generated MOF.
    """ 

    path_parts = path.split('~')
    cif_name = path_parts[0]
    user_ID = path_parts[1]
    return flask.send_from_directory('temp_file_creation_' + user_ID + '/tobacco_3.0/output_cifs', cif_name);

@app.route('/CoRE2019/<path:path>') # needed for fetch
def serve_CoRE_MOF(path):
    """
    serve_CoRE_MOF returns a file to MOFSimplify.
    The file is intended to be a cif file. It should be a CoRE MOF.

    :param path: The path to the desired MOF in the CoRE2019 folder.
    :return: The cif file for the neighbor MOF.
    """ 
    return flask.send_from_directory('CoRE2019', path)

@app.route('/unoptimized_geo/<path:path>')
def serve_unoptimized_MOF(path):
    """
    serve_unoptimized_MOF returns a file to MOFSimplify.
    The file is intended to be a cif file. It should be a MOF.

    :param path: The path to the desired MOF in the stable_building_blocks folder.
    :return: The cif file for the non geometry optimized MOF.
    """ 

    MOF_type = type_determination(path)

    return flask.send_from_directory(f'stable_building_blocks/initial_structures/{MOF_type}', path)

@app.route('/optimized_geo/<path:path>')
def serve_optimized_MOF(path):
    """
    serve_optimized_MOF returns a file to MOFSimplify.
    The file is intended to be a cif file. It should be a MOF.

    :param path: The path to the desired MOF in the stable_building_blocks folder.
    :return: The cif file for the geometry optimized MOF.
    """ 

    MOF_type = type_determination(path)

    return flask.send_from_directory(f'stable_building_blocks/optimized_structures/{MOF_type}', 'optimized_' + path)

@app.route('/ris_files/MOFSimplify_citation.ris')
def serve_ris():
    """
    serve_ris returns a file to MOFSimplify.
    The file is intended the citation file for the MOFSimplify paper.

    :return: The rif citation file.
    """ 

    return flask.send_from_directory('ris_files', 'MOFSimplify_citation.ris')

def listdir_nohidden(path): # used for bb_generate. Ignores hidden files
    """
    listdir_nohidden returns files in the current directory that are not hidden. 
    It is used as a helper function in the bb_generate function.

    :param path: The path to be examined.
    :return: The non-hidden files in the current path.
    """ 
    myList = os.listdir(path);
    for i in myList:
        if i.startswith('.'):
            myList.remove(i)
    return myList

def file_age_in_seconds(pathname): 
    """
    file_age_in_seconds returns the age of the file/folder specified in pathname since the last modification.
    It is used as a helper function in the set_ID function.

    :return: The age of the file/folder specified in pathname since the last modification, in seconds.
    """ 

    return time.time() - os.stat(pathname)[stat.ST_MTIME] # time since last modification

def extract_info(my_data):
    structure = my_data['structure']
    name = my_data['name']
    if name[-4:] == '.cif':
        name = name[:-4] # remove the .cif part of the name
    return structure, name

def extract_info_C2(my_data):
    structure, name = extract_info(my_data)
    temperature = my_data['T']
    pressure = my_data['P']

    return structure, name, temperature, pressure

@app.route('/curr_users', methods=['GET'])
def curr_num_users():
    """
    curr_num_users returns the current number of users on MOFSimplify.
    This is determined by looking at user specific folders. User specific folders that have not been used for a while are deleted (see the set_ID function). 

    :return: The number of extant user-specific folders on MOFSimplify.
    """ 

    sum = 0

    for root, dirs, files in os.walk(MOFSIMPLIFY_PATH):
        for dir in dirs:
            target_str = 'temp_file_creation'
            if len(dir) > len(target_str) and target_str in dir: 
                # target_str in dir since all copies start with temp_file_creation
                # len(dir) > len(target_str) to prevent counting the original temp_file_creation folder
                sum += 1

    return str(sum+1)

@app.route('/get_bb_generated_MOF', methods=['POST']) 
def bb_generate():
    """
    bb_generated generates a MOF using the building blocks and net specified by the user. 
    The function uses ToBaCCo code, version 3.0. 
    It returns the constructed MOF's name to the front end

    :return: The name of the building block MOF.
    """

    tobacco_folder = MOFSIMPLIFY_PATH + "temp_file_creation_" + str(session['ID']) + "/tobacco_3.0/"

    # Grab data
    my_data = json.loads(flask.request.get_data())

    linker = my_data['linker']
    sbu = my_data['sbu']
    net = my_data['net']

    # clear the edges, nodes, templates, and output cifs folders to start fresh
        # when running python tobacco.py, it looks in these folders

    shutil.rmtree(tobacco_folder + 'edges')
    os.mkdir(tobacco_folder + 'edges')
    shutil.rmtree(tobacco_folder + 'nodes')
    os.mkdir(tobacco_folder + 'nodes')
    shutil.rmtree(tobacco_folder + 'templates')
    os.mkdir(tobacco_folder + 'templates')
    shutil.rmtree(tobacco_folder + 'output_cifs')
    os.mkdir(tobacco_folder + 'output_cifs')

    # copy over the linker, sbu, and net specified by the user in the edges, nodes, and templates folders

    shutil.copy(tobacco_folder + 'edges_database/' + linker + '.cif', tobacco_folder + 'edges/' + linker + '.cif')
    shutil.copy(tobacco_folder + 'nodes_database/' + sbu + '.cif', tobacco_folder + 'nodes/' + sbu + '.cif')
    shutil.copy(tobacco_folder + 'template_database/' + net + '.cif', tobacco_folder + 'templates/' + net + '.cif')


    # run the command to construct the MOF
    os.chdir(tobacco_folder) 
    # note: os.chdir here could cause issues if multiple users are using the website and try to make a building block generated MOF at the same time, since MOFSimplify server might chdir when someone else is in the middle of an operation
    # luckily, it is a quick operation, so this is unlikely
    subprocess.run(['python', 'tobacco.py']) 
    os.chdir(MOFSIMPLIFY_PATH)

    # if successful, there will be an output cif in the folder output_cifs
    if listdir_nohidden(tobacco_folder + 'output_cifs') == []: # no files in folder
        print('Construction failed.')
        return 'FAILED'
    
    constructed_MOF = listdir_nohidden(tobacco_folder + 'output_cifs')
    constructed_MOF = constructed_MOF[0] # getting the first, and only, element out of the list

    dictionary = {};

    dictionary['mof_name'] = constructed_MOF

    # getting the primitive cell using molSimplify
    get_primitive(tobacco_folder + 'output_cifs/' + constructed_MOF, tobacco_folder + 'output_cifs/primitive_' + constructed_MOF);

    json_object = json.dumps(dictionary, indent = 4);

    return json_object


### New section. The code for thermal and solvent stability. ### 


def normalize_data_solvent(df_train, df_newMOF, fnames, lname, debug=False):
    """
    normalize_data_solvent takes in two DataFrames df_train and df_newMOF, one for the training data (many rows) and one for the new MOF (one row) for which a prediction is to be generated.
    This function also takes in fnames (the feature names) and lname (the target property name).
    This function normalizes the X values from the pandas DataFrames and returns them as X_train and X_newMOF.
    It also standardizes y_train, which are the solvent removal stability flags in the training data DataFrame, and returns x_scaler (which scaled X_train).
        By standardizes, I mean that it makes the values of y_train either 0 or 1

    :param df_train: A pandas DataFrame of the training data.
    :param df_newMOF: A pandas DataFrame of the new MOF being analyzed.
    :param fnames: An array of column names of the descriptors.
    :param lname: An array of the column name of the target.
    :param debug: A boolean that determines whether extra information is printed.
    :return: numpy.ndarray X_train, the descriptors of the training data. Its number of rows is the number of MOFs in the training data. Its number of columns is the number of descriptors.
    :return: numpy.ndarray X_newMOF, the descriptors of the new MOF being analyzed by MOFSimplify. It contains only one row.
    :return: numpy.ndarray y_train, the solvent removal stabilities of the training data. 
    :return: sklearn.preprocessing._data.StandardScaler x_scaler, the scaler used to normalize the descriptor data to unit mean and a variance of 1.
    """ 
    _df_train = df_train.copy().dropna(subset=fnames+lname)
    _df_newMOF = df_newMOF.copy().dropna(subset=fnames) 
    X_train, X_newMOF = _df_train[fnames].values, _df_newMOF[fnames].values # takes care of ensuring ordering is same for both X
    y_train = _df_train[lname].values
    if debug:
        print("training data reduced from %d -> %d because of nan." % (len(df_train), y_train.shape[0]))
    x_scaler = sklearn.preprocessing.StandardScaler()
    x_scaler.fit(X_train)
    X_train = x_scaler.transform(X_train)
    X_newMOF = x_scaler.transform(X_newMOF)
    y_train = np.array([1 if x == 1 else 0 for x in y_train.reshape(-1, )])
    return X_train, X_newMOF, y_train, x_scaler

def standard_labels(df, key="flag"):
    """
    standard_labels makes the solvent removal stability either 1 (stable upon solvent removal) or 0 (unstable upon solvent removal)
    "flag" is the column under which solvent removal stability is reported in the DataFrame

    :param df: A pandas DataFrame to modify.
    :param key: The column in the pandas DataFrame to look at.
    :return: The modified pandas DataFrame.
    """ 
    flags = [1 if row[key] == 1 else 0 for _, row in df.iterrows()] # Look through all rows of the DataFrame df.
    df[key] = flags
    return df

def run_solvent_ANN(user_id, path, MOF_name, solvent_ANN):
    """
    run_solvent_ANN runs the solvent removal stability ANN with the desired MOF as input.
    It returns a prediction between zero and one. This prediction corresponds to an assessment of MOF stability upon solvent removal.
    The further the prediction is from 0.5, the more sure the ANN is.

    :param user_id: str, the session ID of the user
    :param path: str, the server's path to the MOFSimplify folder on the server
    :param MOF_name: str, the name of the MOF for which a prediction is being generated
    :param solvent_ANN: keras.engine.training.Model, the ANN itself
    :return: str str(new_MOF_pred[0][0]), the model solvent removal stability prediction 
    :return: list neighbor_names, the latent space nearest neighbor MOFs in the solvent removal stability ANN
    :return: list neighbor_distances, the latent space distances of the latent space nearest neighbor MOFs in neighbor_names 
    """ 

    RACs = ['D_func-I-0-all','D_func-I-1-all','D_func-I-2-all','D_func-I-3-all',
     'D_func-S-0-all', 'D_func-S-1-all', 'D_func-S-2-all', 'D_func-S-3-all',
     'D_func-T-0-all', 'D_func-T-1-all', 'D_func-T-2-all', 'D_func-T-3-all',
     'D_func-Z-0-all', 'D_func-Z-1-all', 'D_func-Z-2-all', 'D_func-Z-3-all',
     'D_func-chi-0-all', 'D_func-chi-1-all', 'D_func-chi-2-all',
     'D_func-chi-3-all', 'D_lc-I-0-all', 'D_lc-I-1-all', 'D_lc-I-2-all',
     'D_lc-I-3-all', 'D_lc-S-0-all', 'D_lc-S-1-all', 'D_lc-S-2-all',
     'D_lc-S-3-all', 'D_lc-T-0-all', 'D_lc-T-1-all', 'D_lc-T-2-all',
     'D_lc-T-3-all', 'D_lc-Z-0-all', 'D_lc-Z-1-all', 'D_lc-Z-2-all',
     'D_lc-Z-3-all', 'D_lc-chi-0-all', 'D_lc-chi-1-all', 'D_lc-chi-2-all',
     'D_lc-chi-3-all', 'D_mc-I-0-all', 'D_mc-I-1-all', 'D_mc-I-2-all',
     'D_mc-I-3-all', 'D_mc-S-0-all', 'D_mc-S-1-all', 'D_mc-S-2-all',
     'D_mc-S-3-all', 'D_mc-T-0-all', 'D_mc-T-1-all', 'D_mc-T-2-all',
     'D_mc-T-3-all', 'D_mc-Z-0-all', 'D_mc-Z-1-all', 'D_mc-Z-2-all',
     'D_mc-Z-3-all', 'D_mc-chi-0-all', 'D_mc-chi-1-all', 'D_mc-chi-2-all',
     'D_mc-chi-3-all', 'f-I-0-all', 'f-I-1-all', 'f-I-2-all', 'f-I-3-all',
     'f-S-0-all', 'f-S-1-all', 'f-S-2-all', 'f-S-3-all', 'f-T-0-all', 'f-T-1-all',
     'f-T-2-all', 'f-T-3-all', 'f-Z-0-all', 'f-Z-1-all', 'f-Z-2-all', 'f-Z-3-all',
     'f-chi-0-all', 'f-chi-1-all', 'f-chi-2-all', 'f-chi-3-all', 'f-lig-I-0',
     'f-lig-I-1', 'f-lig-I-2', 'f-lig-I-3', 'f-lig-S-0', 'f-lig-S-1', 'f-lig-S-2',
     'f-lig-S-3', 'f-lig-T-0', 'f-lig-T-1', 'f-lig-T-2', 'f-lig-T-3', 'f-lig-Z-0',
     'f-lig-Z-1', 'f-lig-Z-2', 'f-lig-Z-3', 'f-lig-chi-0', 'f-lig-chi-1',
     'f-lig-chi-2', 'f-lig-chi-3', 'func-I-0-all', 'func-I-1-all',
     'func-I-2-all', 'func-I-3-all', 'func-S-0-all', 'func-S-1-all',
     'func-S-2-all', 'func-S-3-all', 'func-T-0-all', 'func-T-1-all',
     'func-T-2-all', 'func-T-3-all', 'func-Z-0-all', 'func-Z-1-all',
     'func-Z-2-all', 'func-Z-3-all', 'func-chi-0-all', 'func-chi-1-all',
     'func-chi-2-all', 'func-chi-3-all', 'lc-I-0-all', 'lc-I-1-all', 'lc-I-2-all',
     'lc-I-3-all', 'lc-S-0-all', 'lc-S-1-all', 'lc-S-2-all', 'lc-S-3-all',
     'lc-T-0-all', 'lc-T-1-all', 'lc-T-2-all', 'lc-T-3-all', 'lc-Z-0-all',
     'lc-Z-1-all', 'lc-Z-2-all', 'lc-Z-3-all', 'lc-chi-0-all', 'lc-chi-1-all',
     'lc-chi-2-all', 'lc-chi-3-all', 'mc-I-0-all', 'mc-I-1-all', 'mc-I-2-all',
     'mc-I-3-all', 'mc-S-0-all', 'mc-S-1-all', 'mc-S-2-all', 'mc-S-3-all',
     'mc-T-0-all', 'mc-T-1-all', 'mc-T-2-all', 'mc-T-3-all', 'mc-Z-0-all',
     'mc-Z-1-all', 'mc-Z-2-all', 'mc-Z-3-all', 'mc-chi-0-all', 'mc-chi-1-all',
     'mc-chi-2-all', 'mc-chi-3-all']
    geo = ['Df','Di', 'Dif','GPOAV','GPONAV','GPOV','GSA','POAV','POAV_vol_frac',
      'PONAV','PONAV_vol_frac','VPOV','VSA','cell_v']
     
    ANN_path = path + 'model/solvent/ANN/'
    temp_file_path = path + 'temp_file_creation_' + user_id + '/'
    df_train = pd.read_csv(ANN_path+'dropped_connectivity_dupes/train.csv')
    df_train = df_train.loc[:, (df_train != df_train.iloc[0]).any()]
    df_newMOF = pd.read_csv(temp_file_path + 'merged_descriptors/' + MOF_name + '_descriptors.csv') # assumes that temp_file_creation/ is in parent folder
    features = [val for val in df_train.columns.values if val in RACs+geo]

    df_train = standard_labels(df_train, key="flag")

    # The normalize_data_solvent function is expecting a DataFrame with each MOF in a separate row, and features in columns

    ### Utilize the function below to normalize the RACs + geo features of the new MOF
    # newMOF refers to the MOF that has been uploaded to MOFSimplify, for which a prediction will be generated
    X_train, X_newMOF, y_train, x_scaler = normalize_data_solvent(df_train, df_newMOF, features, ["flag"], debug=False)
    # Order of values in X_newMOF matters, but this is taken care of in normalize_data_solvent.
    X_train.shape, y_train.reshape(-1, ).shape
    model = solvent_ANN

    from tensorflow.python.keras.backend import set_session
    with tf_session.as_default(): # session stuff is needed because the model was loaded from h5 a while ago
        with tf_session.graph.as_default():
            ### new_MOF_pred will be a decimal value between 0 and 1, below 0.5 is unstable, above 0.5 is stable
            new_MOF_pred = np.round(model.predict(X_newMOF),2) # round to 2 decimals

            # Define the function for the latent space. This will depend on the model. We want the layer before the last, in this case this was the 12th one.
            get_latent = K.function([model.layers[0].input],
                                    [model.layers[12].output]) # Last layer before dense-last

            # Get the latent vectors for the training data first, then the latent vectors for the test data.
            training_latent = get_latent([X_train, 0])[0]
            design_latent = get_latent([X_newMOF, 0])[0]

    # Compute the pairwise distances between the test latent vectors and the train latent vectors to get latent distances
    d1 = pairwise_distances(design_latent,training_latent,n_jobs=30)
    df1 = pd.DataFrame(data=d1, columns=df_train['CoRE_name'].tolist())
    df1.to_csv(temp_file_path + 'solvent_test_latent_dists.csv')

    # Want to find the closest points (let's say the closest 5 points); so, smallest values in df1
    neighbors = 5 # number of closest points

    # will make arrays of length neighbors, where each entry is the next closest neighbor (will do this for both names and distances)
    neighbors_names = []
    neighbors_distances = []

    df_reformat = df1.min(axis='index')

    for i in range(neighbors):
        name = df_reformat.idxmin() # name of next closest complex in the training data
        distance = df_reformat.min() # distance of the next closest complex in the training data to the new MOF
        df_reformat = df_reformat.drop(name) # dropping the next closest complex, in order to find the next-next closest complex

        neighbors_names.append(name)
        neighbors_distances.append(str(distance))

    return str(new_MOF_pred[0][0]), neighbors_names, neighbors_distances


def normalize_data_thermal(df_train, df_newMOF, fnames, lname, debug=False): # Function assumes it gets pandas DataFrames with MOFs as rows and features as columns
    """
    normalize_data_thermal takes in two DataFrames df_train and df_newMOF, one for the training data (many rows) and one for the new MOF (one row) for which a prediction is to be generated.
    This function also takes in fnames (the feature names) and lname (the target property name).
    This function normalizes the X values from the pandas DataFrames and returns them as X_train and X_newMOF.
    It also normalizes y_train, which are the thermal breakdown temperatures in the training data DataFrame, and returns x_scaler (which scaled X_train) and y_scaler (which scaled y_train).

    :param df_train: A pandas DataFrame of the training data.
    :param df_newMOF: A pandas DataFrame of the new MOF being analyzed.
    :param fnames: An array of column names of the descriptors.
    :param lname: An array of the column name of the target.
    :param debug: A boolean that determines whether extra information is printed.
    :return: numpy.ndarray X_train, the descriptors of the training data. Its number of rows is the number of MOFs in the training data. Its number of columns is the number of descriptors.
    :return: numpy.ndarray X_newMOF, the descriptors of the new MOF being analyzed by MOFSimplify. It contains only one row.
    :return: numpy.ndarray y_train, the thermal stabilities of the training data. 
    :return: sklearn.preprocessing._data.StandardScaler x_scaler, the scaler used to normalize the descriptor data to unit mean and a variance of 1. 
    :return: sklearn.preprocessing._data.StandardScaler y_scaler, the scaler used to normalize the target data to unit mean and a variance of 1.
    """ 
    _df_train = df_train.copy().dropna(subset=fnames+lname)
    _df_newMOF = df_newMOF.copy().dropna(subset=fnames) 
    X_train, X_newMOF = _df_train[fnames].values, _df_newMOF[fnames].values # takes care of ensuring ordering is same for both X
    y_train = _df_train[lname].values
    if debug:
        print("training data reduced from %d -> %d because of nan." % (len(df_train), y_train.shape[0]))
    x_scaler = sklearn.preprocessing.StandardScaler()
    x_scaler.fit(X_train)
    X_train = x_scaler.transform(X_train)
    X_newMOF = x_scaler.transform(X_newMOF)
    y_scaler = sklearn.preprocessing.StandardScaler()
    y_scaler.fit(y_train)
    y_train = y_scaler.transform(y_train)
    return X_train, X_newMOF, y_train, x_scaler, y_scaler

def run_thermal_ANN(user_id, path, MOF_name, thermal_ANN):
    """
    run_thermal_ANN runs the thermal stability ANN with the desired MOF as input.
    It returns a prediction for the thermal breakdown temperature of the chosen MOF.

    :param user_id: str, the session ID of the user
    :param path: str, the server's path to the MOFSimplify folder on the server
    :param MOF_name: str, the name of the MOF for which a prediction is being generated
    :param thermal_ANN: keras.engine.training.Model, the ANN itself
    :return: str new_MOF_pred, the model thermal stability prediction 
    :return: list neighbor_names, the latent space nearest neighbor MOFs in the thermal stability ANN
    :return: list neighbor_distances, the latent space distances of the latent space nearest neighbor MOFs in neighbor_names 
    """ 

    RACs = ['D_func-I-0-all','D_func-I-1-all','D_func-I-2-all','D_func-I-3-all',
     'D_func-S-0-all', 'D_func-S-1-all', 'D_func-S-2-all', 'D_func-S-3-all',
     'D_func-T-0-all', 'D_func-T-1-all', 'D_func-T-2-all', 'D_func-T-3-all',
     'D_func-Z-0-all', 'D_func-Z-1-all', 'D_func-Z-2-all', 'D_func-Z-3-all',
     'D_func-chi-0-all', 'D_func-chi-1-all', 'D_func-chi-2-all',
     'D_func-chi-3-all', 'D_lc-I-0-all', 'D_lc-I-1-all', 'D_lc-I-2-all',
     'D_lc-I-3-all', 'D_lc-S-0-all', 'D_lc-S-1-all', 'D_lc-S-2-all',
     'D_lc-S-3-all', 'D_lc-T-0-all', 'D_lc-T-1-all', 'D_lc-T-2-all',
     'D_lc-T-3-all', 'D_lc-Z-0-all', 'D_lc-Z-1-all', 'D_lc-Z-2-all',
     'D_lc-Z-3-all', 'D_lc-chi-0-all', 'D_lc-chi-1-all', 'D_lc-chi-2-all',
     'D_lc-chi-3-all', 'D_mc-I-0-all', 'D_mc-I-1-all', 'D_mc-I-2-all',
     'D_mc-I-3-all', 'D_mc-S-0-all', 'D_mc-S-1-all', 'D_mc-S-2-all',
     'D_mc-S-3-all', 'D_mc-T-0-all', 'D_mc-T-1-all', 'D_mc-T-2-all',
     'D_mc-T-3-all', 'D_mc-Z-0-all', 'D_mc-Z-1-all', 'D_mc-Z-2-all',
     'D_mc-Z-3-all', 'D_mc-chi-0-all', 'D_mc-chi-1-all', 'D_mc-chi-2-all',
     'D_mc-chi-3-all', 'f-I-0-all', 'f-I-1-all', 'f-I-2-all', 'f-I-3-all',
     'f-S-0-all', 'f-S-1-all', 'f-S-2-all', 'f-S-3-all', 'f-T-0-all', 'f-T-1-all',
     'f-T-2-all', 'f-T-3-all', 'f-Z-0-all', 'f-Z-1-all', 'f-Z-2-all', 'f-Z-3-all',
     'f-chi-0-all', 'f-chi-1-all', 'f-chi-2-all', 'f-chi-3-all', 'f-lig-I-0',
     'f-lig-I-1', 'f-lig-I-2', 'f-lig-I-3', 'f-lig-S-0', 'f-lig-S-1', 'f-lig-S-2',
     'f-lig-S-3', 'f-lig-T-0', 'f-lig-T-1', 'f-lig-T-2', 'f-lig-T-3', 'f-lig-Z-0',
     'f-lig-Z-1', 'f-lig-Z-2', 'f-lig-Z-3', 'f-lig-chi-0', 'f-lig-chi-1',
     'f-lig-chi-2', 'f-lig-chi-3', 'func-I-0-all', 'func-I-1-all',
     'func-I-2-all', 'func-I-3-all', 'func-S-0-all', 'func-S-1-all',
     'func-S-2-all', 'func-S-3-all', 'func-T-0-all', 'func-T-1-all',
     'func-T-2-all', 'func-T-3-all', 'func-Z-0-all', 'func-Z-1-all',
     'func-Z-2-all', 'func-Z-3-all', 'func-chi-0-all', 'func-chi-1-all',
     'func-chi-2-all', 'func-chi-3-all', 'lc-I-0-all', 'lc-I-1-all', 'lc-I-2-all',
     'lc-I-3-all', 'lc-S-0-all', 'lc-S-1-all', 'lc-S-2-all', 'lc-S-3-all',
     'lc-T-0-all', 'lc-T-1-all', 'lc-T-2-all', 'lc-T-3-all', 'lc-Z-0-all',
     'lc-Z-1-all', 'lc-Z-2-all', 'lc-Z-3-all', 'lc-chi-0-all', 'lc-chi-1-all',
     'lc-chi-2-all', 'lc-chi-3-all', 'mc-I-0-all', 'mc-I-1-all', 'mc-I-2-all',
     'mc-I-3-all', 'mc-S-0-all', 'mc-S-1-all', 'mc-S-2-all', 'mc-S-3-all',
     'mc-T-0-all', 'mc-T-1-all', 'mc-T-2-all', 'mc-T-3-all', 'mc-Z-0-all',
     'mc-Z-1-all', 'mc-Z-2-all', 'mc-Z-3-all', 'mc-chi-0-all', 'mc-chi-1-all',
     'mc-chi-2-all', 'mc-chi-3-all']
    geo = ['Df','Di', 'Dif','GPOAV','GPONAV','GPOV','GSA','POAV','POAV_vol_frac',
      'PONAV','PONAV_vol_frac','VPOV','VSA','cell_v']
     
    ANN_path = path + 'model/thermal/ANN/'
    temp_file_path = path + 'temp_file_creation_' + user_id + '/'
    df_train_all = pd.read_csv(ANN_path+"train.csv").append(pd.read_csv(ANN_path+"val.csv"))
    df_train = pd.read_csv(ANN_path+"train.csv")
    df_train = df_train.loc[:, (df_train != df_train.iloc[0]).any()]
    df_newMOF = pd.read_csv(temp_file_path + 'merged_descriptors/' + MOF_name + '_descriptors.csv') # Assume temp_file_creation/ in parent directory
    features = [val for val in df_train.columns.values if val in RACs+geo]

    X_train, X_newMOF, y_train, x_scaler, y_scaler = normalize_data_thermal(df_train, df_newMOF, features, ["T"], debug=False)
    X_train.shape, y_train.reshape(-1, ).shape 

    model = thermal_ANN

    from tensorflow.python.keras.backend import set_session
    with tf_session.as_default():
        with tf_session.graph.as_default():
            new_MOF_pred = y_scaler.inverse_transform(model.predict(X_newMOF))
            new_MOF_pred = np.round(new_MOF_pred,1) # round to 1 decimal

            # isolating just the prediction, since the model spits out the prediction like [[PREDICTION]], as in, in hard brackets
            new_MOF_pred = new_MOF_pred[0][0]
            new_MOF_pred = str(new_MOF_pred)

            # adding units
            degree_sign= u'\N{DEGREE SIGN}'
            new_MOF_pred = new_MOF_pred + degree_sign + 'C' # degrees Celsius

            # Define the function for the latent space. This will depend on the model. We want the layer before the last, in this case this was the 8th one.
            get_latent = K.function([model.layers[0].input],
                                    [model.layers[8].output]) # Last layer before dense-last

            # Get the latent vectors for the training data first, then the latent vectors for the test data.
            training_latent = get_latent([X_train, 0])[0]
            design_latent = get_latent([X_newMOF, 0])[0]

            print(training_latent.shape,design_latent.shape)

    # Compute the pairwise distances between the test latent vectors and the train latent vectors to get latent distances
    d1 = pairwise_distances(design_latent,training_latent,n_jobs=30)
    df1 = pd.DataFrame(data=d1, columns=df_train['CoRE_name'].tolist())
    df1.to_csv(temp_file_path + 'solvent_test_latent_dists.csv')

    # Want to find the closest points (let's say the closest 5 points); so, smallest values in df1
    neighbors = 5 # number of closest points

    # will make arrays of length neighbors, where each entry is the next closest neighbor (will do this for both names and distances)
    neighbors_names = []
    neighbors_distances = []

    df_reformat = df1.min(axis='index')

    for i in range(neighbors):
        name = df_reformat.idxmin() # name of next closest complex in the training data
        distance = df_reformat.min() # distance of the next closest complex in the training data to the new MOF
        df_reformat = df_reformat.drop(name) # dropping the next closest complex, in order to find the next-next closest complex

        neighbors_names.append(name)
        neighbors_distances.append(str(distance))

    return new_MOF_pred, neighbors_names, neighbors_distances


def descriptor_generator(name, structure, prediction_type, is_entry):
    """
    # descriptor_generator is used by both ss_predict() and ts_predict() to generate RACs and Zeo++ descriptors.
    # These descriptors are subsequently used in ss_predict() and ts_predict() for the ANN models.
    # Inputs are the name of the MOF and the structure (cif file text) of the MOF for which descriptors are to be generated.
    # The third input indicates the type of prediction (solvent removal or thermal).

    :param name: str, the name of the MOF being analyzed.
    :param structure: str, the text of the cif file of the MOF being analyzed.
    :param prediction_type: str, the type of prediction being run. Can either be 'solvent' or 'thermal'.
    :param is_entry: boolean, indicates whether the descriptor CSV has already been written.
    :return: Depends, either the string 'FAILED' if descriptor generation fails, a dictionary myDict (if the MOF being analyzed is in the training data), or an array myResult (if the MOF being analyzed is not in the training data) 
    """ 

    print('TIME CHECK 2')
    timeStarted = time.time() # save start time (debugging)

    temp_file_folder = MOFSIMPLIFY_PATH + "temp_file_creation_" + str(session['ID']) + '/'
    cif_folder = temp_file_folder + 'cifs/'

    # Write the data back to a cif file.
    try:
        cif_file = open(cif_folder + name + '.cif', 'w')
    except FileNotFoundError:
        return 'FAILED'
    cif_file.write(structure)
    cif_file.close()

    # There can be a RACs folder for solvent predictions and a RACs folder for thermal predictions. Same for Zeo++.
    RACs_folder = temp_file_folder + 'feature_generation/' + prediction_type + '_RACs/'
    zeo_folder = temp_file_folder + 'feature_generation/' + prediction_type + '_zeo++/'

    # Delete the RACs folder, then remake it (to start fresh for this prediction).
    shutil.rmtree(RACs_folder)
    os.mkdir(RACs_folder)

    # Doing the same with the Zeo++ folder.
    shutil.rmtree(zeo_folder)
    os.mkdir(zeo_folder)

    if not is_entry: # have to generate the CSV

        # Next, running MOF featurization
        try:
            get_primitive(cif_folder + name + '.cif', cif_folder + name + '_primitive.cif');
        except ValueError:
            return 'FAILED'

        timeDelta = time.time() - timeStarted # get execution time
        print('Finished process in ' + str(timeDelta) + ' seconds')

        print('TIME CHECK 3')

        timeStarted = time.time() # save start time (debugging)

        # get_MOF_descriptors is used in RAC_getter.py to get RAC features.
            # The files that are generated from RAC_getter.py: lc_descriptors.csv, sbu_descriptors.csv, linker_descriptors.csv

        # cmd1, cmd2, and cmd3 are for Zeo++. cm4 is for RACs.
        cmd1 = MOFSIMPLIFY_PATH + 'zeo++-0.3/network -ha -res ' + zeo_folder + name + '_pd.txt ' + cif_folder + name + '_primitive.cif'
        cmd2 = MOFSIMPLIFY_PATH + 'zeo++-0.3/network -sa 1.86 1.86 10000 ' + zeo_folder + name + '_sa.txt ' + cif_folder + name + '_primitive.cif'
        cmd3 = MOFSIMPLIFY_PATH + 'zeo++-0.3/network -volpo 1.86 1.86 10000 ' + zeo_folder + name + '_pov.txt '+ cif_folder + name + '_primitive.cif'
        cmd4 = 'python ' + MOFSIMPLIFY_PATH + 'model/RAC_getter.py %s %s %s' %(cif_folder, name, RACs_folder)

        # four parallelized Zeo++ and RAC commands
        process1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=None, shell=True)
        process2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, stderr=None, shell=True)
        process3 = subprocess.Popen(cmd3, stdout=subprocess.PIPE, stderr=None, shell=True)
        process4 = subprocess.Popen(cmd4, stdout=subprocess.PIPE, stderr=None, shell=True)

        output1 = process1.communicate()[0]
        output2 = process2.communicate()[0]
        output3 = process3.communicate()[0]
        output4 = process4.communicate()[0]

        # Have written output of Zeo++ commands to files. Now, code below extracts information from those files.

        ''' 
        The geometric descriptors are:
        - the maximum included sphere (Di)
        - maximum free sphere (Df)
        - maximum included sphere in the free sphere path (Dif)
        - gravimetric pore volume (GPOV)
        - volumetric pore volume (VPOV)
        - gravimetric surface area (GSA)
        - volumetric surface area (VSA)
        - cell volume (cell_v)
        - gravimetric pore accessible volume (GPOAV)
        - gravimetric pore non-accessible volume (GPONAV)
        - pore-accessible volume (POAV)
        - pore non-accessible volume (PONAV)
        - pore accessible void fraction (POAVF)
        - pore nonaccessible void fraction (PONAVF)

        All Zeo++ calculations use a probe radius of 1.86 angstrom, and zeo++ is called by subprocess.
        '''

        dict_list = []
        cif_file = name + '_primitive.cif' 
        basename = cif_file.strip('.cif')
        largest_included_sphere, largest_free_sphere, largest_included_sphere_along_free_sphere_path = np.nan, np.nan, np.nan
        unit_cell_volume, crystal_density, VSA, GSA = np.nan, np.nan, np.nan, np.nan
        VPOV, GPOV = np.nan, np.nan
        POAV, PONAV, GPOAV, GPONAV, POAV_volume_fraction, PONAV_volume_fraction = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

        if (os.path.exists(zeo_folder + name + '_pd.txt') & os.path.exists(zeo_folder + name + '_sa.txt') &
            os.path.exists(zeo_folder + name + '_pov.txt')):
            with open(zeo_folder + name + '_pd.txt') as f:
                pore_diameter_data = f.readlines()
                for row in pore_diameter_data:
                    largest_included_sphere = float(row.split()[1]) # largest included sphere
                    largest_free_sphere = float(row.split()[2]) # largest free sphere
                    largest_included_sphere_along_free_sphere_path = float(row.split()[3]) # largest included sphere along free sphere path
            with open(zeo_folder + name + '_sa.txt') as f:
                surface_area_data = f.readlines()
                for i, row in enumerate(surface_area_data):
                    if i == 0:
                        unit_cell_volume = float(row.split('Unitcell_volume:')[1].split()[0]) # unit cell volume
                        crystal_density = float(row.split('Density:')[1].split()[0]) # crystal density
                        VSA = float(row.split('ASA_m^2/cm^3:')[1].split()[0]) # volumetric surface area
                        GSA = float(row.split('ASA_m^2/g:')[1].split()[0]) # gravimetric surface area
            with open(zeo_folder + name + '_pov.txt') as f:
                pore_volume_data = f.readlines()
                for i, row in enumerate(pore_volume_data):
                    if i == 0:
                        density = float(row.split('Density:')[1].split()[0])
                        POAV = float(row.split('POAV_A^3:')[1].split()[0]) # Probe accessible pore volume
                        PONAV = float(row.split('PONAV_A^3:')[1].split()[0]) # Probe non-accessible probe volume
                        GPOAV = float(row.split('POAV_cm^3/g:')[1].split()[0])
                        GPONAV = float(row.split('PONAV_cm^3/g:')[1].split()[0])
                        POAV_volume_fraction = float(row.split('POAV_Volume_fraction:')[1].split()[0]) # probe accessible volume fraction
                        PONAV_volume_fraction = float(row.split('PONAV_Volume_fraction:')[1].split()[0]) # probe non accessible volume fraction
                        VPOV = POAV_volume_fraction+PONAV_volume_fraction
                        GPOV = VPOV/density
        else:
            print('Not all 3 files exist, so at least one Zeo++ call failed!', 'sa: ',os.path.exists(zeo_folder + name + '_sa.txt'), 
                  '; pd: ',os.path.exists(zeo_folder + name + '_pd.txt'), '; pov: ', os.path.exists(zeo_folder + name + '_pov.txt'))
            return 'FAILED'
        geo_dict = {'name': basename, 'cif_file': cif_file, 'Di': largest_included_sphere, 'Df': largest_free_sphere, 'Dif': largest_included_sphere_along_free_sphere_path,
                    'cell_v': unit_cell_volume, 'VSA': VSA, 'GSA': GSA, 'VPOV': VPOV, 'GPOV': GPOV, 'POAV_vol_frac': POAV_volume_fraction, 
                    'PONAV_vol_frac': PONAV_volume_fraction, 'GPOAV': GPOAV, 'GPONAV': GPONAV, 'POAV': POAV, 'PONAV': PONAV}
        dict_list.append(geo_dict)
        geo_df = pd.DataFrame(dict_list)
        geo_df.to_csv(zeo_folder + 'geometric_parameters.csv',index=False)

        # error handling for cmd4
        with open(RACs_folder + 'RAC_getter_log.txt', 'r') as f:
            if f.readline() == 'FAILED':
                print('RAC generation failed.')
                return 'FAILED'


        timeDelta = time.time() - timeStarted # get execution time
        print('Finished process in ' + str(timeDelta) + ' seconds')

        print('TIME CHECK 4')

        timeStarted = time.time() # save start time

        # Merging geometric information with the RAC information that is in the get_MOF_descriptors-generated files (lc_descriptors.csv, sbu_descriptors.csv, linker_descriptors.csv)
        try:
            lc_df = pd.read_csv(RACs_folder + "lc_descriptors.csv") 
            sbu_df = pd.read_csv(RACs_folder + "sbu_descriptors.csv")
            linker_df = pd.read_csv(RACs_folder + "linker_descriptors.csv")
        except Exception: # csv files have been deleted
            return 'FAILED' 

        lc_df = lc_df.mean().to_frame().transpose() # averaging over all rows. Convert resulting Series into a DataFrame, then transpose
        sbu_df = sbu_df.mean().to_frame().transpose()
        linker_df = linker_df.mean().to_frame().transpose()

        merged_df = pd.concat([geo_df, lc_df, sbu_df, linker_df], axis=1)

        merged_df.to_csv(temp_file_folder + '/merged_descriptors/' + name + '_descriptors.csv',index=False) # written in /temp_file_creation_SESSIONID

    else: # CSV is already written
        merged_df = pd.read_csv(temp_file_folder + '/merged_descriptors/' + name + '_descriptors.csv')

    if prediction_type == 'solvent':

        ANN_folder = MOFSIMPLIFY_PATH + 'model/solvent/ANN/'
        train_df = pd.read_csv(ANN_folder + 'dropped_connectivity_dupes/train.csv')

    if prediction_type == 'thermal':
        ANN_folder = MOFSIMPLIFY_PATH + 'model/thermal/ANN/'
        train_df = pd.read_csv(ANN_folder + 'train.csv')

    ### Here, I do a check to see if the current MOF is in the training data. ###
    # If it is, then I return the known truth for the MOF, rather than make a prediction.

    # Will iterate through the rows of the train pandas DataFrame

    in_train = False

    for index, row in train_df.iterrows(): # iterate through rows of the training data MOFs

        row_match = True # gets set to false if any values don't match 
        matching_MOF = None # gets assigned a value if the current MOF is in the training data

        for col in merged_df.columns: # iterate through columns of the single new MOF we are predicting on (merged_df is just one row)
            if col == 'name' or col == 'cif_file' or col == 'Dif':
                continue # skip these
                # Dif was sometimes differing between new Zeo++ call and training data value, for the same MOF

            # If for any property a training MOF and the new MOF we are predicting on differ too much, we know they are not the same MOF
            # So row_match is set to false for this training MOF
            if np.absolute(row[col] - merged_df.iloc[0][col]) > 0.05 * np.absolute(merged_df.iloc[0][col]): # row[col] != merged_df.iloc[0][col] was leading to some same values being identified as different b/c of some floating 10^-15 values 
                row_match = False
                break
        
        if row_match and prediction_type == 'solvent': # all columns for the row match! Some training MOF is the same as the new MOF
            in_train = True
            match_truth = row['flag'] # the flag for the MOF that matches the current MOF
            matching_MOF = row['CoRE_name']
            print(row['CoRE_name'])
            break
        if row_match and prediction_type =='thermal': # all columns for the row match! Some training MOF is the same as the new MOF
            in_train = True
            match_truth = row['T'] # the flag for the MOF that matches the current MOF
            match_truth = np.round(match_truth,1) # round to 1 decimal
            matching_MOF = row['CoRE_name']
            # adding units
            degree_sign= u'\N{DEGREE SIGN}'
            match_truth = str(match_truth) + degree_sign + 'C' # degrees Celsius
            break

    if in_train:
        myDict = {'in_train': True, 'truth': match_truth, 'match': matching_MOF}
        return myDict

    ### End of training data check. ###

    timeDelta = time.time() - timeStarted # get execution time
    print('Finished process in ' + str(timeDelta) + ' seconds')

    print('TIME CHECK 5')

    myResult = [temp_file_folder, ANN_folder]

    return myResult 
    
@app.route('/predict_solvent_stability', methods=['POST']) 
def ss_predict():
    """
    ss_predict generates the solvent removal stability prediction for the selected MOF.
        Or it will return a ground truth if the MOF is in the solvent stability ANN training data.
        Or it will return a signal that descriptor generation failed and thus a prediction cannot be made. 
    RAC featurization and Zeo++ geometry information for the selected MOF is generated, using descriptor_generator.
    Then, Aditya's model is applied to make a prediction using run_solvent_ANN.

    If the structure is already in our history.MOFSimplify_v2 collection, we use information from the database to speed things up.

    :return: dict results, contains the prediction and information on latent space nearest neighbors.
        May instead return a dictionary output if the MOF is in the training data, containing the ground truth and the name of the matching MOF in the training data
        May instead return a string 'FAILED' if descriptor generation fails.
        May instead return a string 'OVERLOAD' if the server is being asked too much of, in order to avoid 50X errors (like 504, aka Gateway Time-out).
    """ 

    global operation_counter # global variable
    global last_operation_counter_clear

    operation_counter_periodic_clear()

    print(f'Current operation_counter: {operation_counter}')

    if operation_counter >= MAX_OPERATIONS:
        return 'OVERLOAD'

    operation_counter += 1

    print('TIME CHECK 1')
    timeStarted = time.time() # save start time (debugging)

    # Grab data and tweak the MOF name.
    my_data = json.loads(flask.request.get_data())
    structure, name = extract_info(my_data)

    timeDelta = time.time() - timeStarted # get execution time
    print('Finished process in ' + str(timeDelta) + ' seconds')

    temp_file_folder = MOFSIMPLIFY_PATH + "temp_file_creation_" + str(session['ID']) + '/'

    ### Check in MongoDB history.MOFSimplify_v2 collection to see if this structure has been predicted on before. 
    # Comment out this section if you are running MOFSimplify on your computer, and define is_entry as False. Also select "No" for the question May MOFSimplify store information on your MOFs?
    client = MongoClient('18.18.63.68',27017) # connect to mongodb. The first argument is the IP address. The second argument is the port.
    db = client.history # The history database
    collection = db.MOFSimplify_v2 # The MOFSimplify_v2 collection in the history database.
        # The collection client.history.MOFSimplify, or db.MOFSimplify, was for a prior version, where the cell_v feature was erroneously called rho
    my_documents = collection.find({'structure':structure})
    is_entry = my_documents.count() # will be zero or one, depending if structure is in the database
    print(f'is_entry is {is_entry}')
    if is_entry: # return statements in this if statement are of same form as those from the workflow if the structure isn't in the database. The frontend needs to receive information of the same form.
        entry_data = my_documents[0]
        if entry_data['failure'] == True: # MOF was not featurizable
            print('entry_data["failure"] is True!')
            if session['permission']:
                db_push_lite(structure=structure, prediction_type='solvent_stability_prediction') # add 1 to the s_times in the database
            operation_counter = conditional_diminish(operation_counter)
            return 'FAILED'
        else: 
            # Making the csv; can skip csv making in descriptor_generator with this
            csv_content = entry_data['csv_content']
            with open(temp_file_folder + 'merged_descriptors/' + name + '_descriptors.csv', 'w') as f:
                f.write(csv_content)

        if entry_data['s_intrain'] == True:
            my_dict = {'in_train':True, 'truth':entry_data['s_result'],'match':None} # don't care about the MOF that matches in CoRE here
            if session['permission']:
                db_push_lite(structure=structure, prediction_type='solvent_stability_prediction') # add 1 to the s_times in the database
            operation_counter = conditional_diminish(operation_counter)
            return my_dict
        elif entry_data['s_intrain'] == False:
            my_dict = {'prediction':entry_data['s_result'],'neighbor_names':entry_data['s_neighbornames'],'neighbor_distances':entry_data['s_neighbordists'],'in_train':False}
            if session['permission']:
                db_push_lite(structure=structure, prediction_type='solvent_stability_prediction') # add 1 to the s_times in the database
            operation_counter = conditional_diminish(operation_counter)
            return my_dict

        # if haven't returned anything by now and entered the if is_entry statement, have a featurizable MOF for which a solvent prediction has not been run, but a thermal prediction has

    ###

    output = descriptor_generator(name, structure, 'solvent', is_entry) # generate descriptors
    if output == 'FAILED': # Description generation failure
        operation_counter = conditional_diminish(operation_counter)
        if session['permission']:
            db_push(structure=structure, prediction_type='solvent_stability_prediction', in_train='', result='', neighbor_names='', neighbor_dists='', failure=True, csv_content='')
        return 'FAILED'
    elif isinstance(output, dict): # MOF was in the training data. Return the ground truth.
        operation_counter = conditional_diminish(operation_counter)

        if session['permission']:
            # Getting the contents of the CSV containing features.
            descriptors_folder = temp_file_folder + "merged_descriptors/"
            with open(descriptors_folder + name + '_descriptors.csv', 'r') as f:
                csv_contents = f.read()

            db_push(structure=structure, prediction_type='solvent_stability_prediction', in_train=True, result=output['truth'], neighbor_names='', neighbor_dists='', failure=False, csv_content=csv_contents)
        return output
    else: # grabbing some variables (that aren't used). We will make a prediction
        temp_file_folder = output[0]
        ANN_folder = output[1]

    # Applying the model next

    timeStarted = time.time() # save start time (debugging)
    prediction, neighbor_names, neighbor_distances = run_solvent_ANN(str(session['ID']), MOFSIMPLIFY_PATH, name, solvent_model)

    results = {'prediction': prediction,
        'neighbor_names': neighbor_names, # The names of the five nearest latent space nearest neighbors in the ANN's training data. The ANNs used CoRE MOFs as training data.
        'neighbor_distances': neighbor_distances, # The latent space distances of the five nearest latent space neighbors.
        'in_train': False} # A prediction was made. The requested (selected) MOF was not in the training data.

    timeDelta = time.time() - timeStarted # get execution time
    print('Finished process in ' + str(timeDelta) + ' seconds') 

    print('TIME CHECK 6')

    if session['permission']:
        # database push of MOF structure and ANN results

        # Getting the contents of the CSV containing features.
        descriptors_folder = temp_file_folder + "merged_descriptors/"
        with open(descriptors_folder + name + '_descriptors.csv', 'r') as f:
            csv_contents = f.read()

        db_push(structure=structure, prediction_type='solvent_stability_prediction', in_train=False, result=prediction, neighbor_names=neighbor_names, neighbor_dists=neighbor_distances, failure=False, csv_content=csv_contents)

    operation_counter = conditional_diminish(operation_counter)
    return results

@app.route('/predict_thermal_stability', methods=['POST']) 
def ts_predict():
    """
    ts_predict generates the thermal stability prediction for the selected MOF.
        Or it will return a ground truth if the MOF is in the thermal stability ANN training data.
        Or it will return a signal that descriptor generation failed and thus a prediction cannot be made. 
    RAC featurization and Zeo++ geometry information for the selected MOF is generated, using descriptor_generator.
    Then, Aditya's model is applied to make a prediction using run_thermal_ANN.

    If the structure is already in our history.MOFSimplify_v2 collection, we use information from the database to speed things up.

    :return: dict results, contains the prediction and information on latent space nearest neighbors.
        May instead return a dictionary output if the MOF is in the training data, containing the ground truth and the name of the matching MOF in the training data
        May instead return a string 'FAILED' if descriptor generation fails.
        May instead return a string 'OVERLOAD' if the server is being asked too much of, in order to avoid 50X errors (like 504, aka Gateway Time-out).
    """ 

    global operation_counter # global variable
    global last_operation_counter_clear

    operation_counter_periodic_clear()

    print(f'Current operation_counter: {operation_counter}')

    if operation_counter >= MAX_OPERATIONS:
        return 'OVERLOAD'

    operation_counter += 1

    # Grab data and tweak the MOF name.
    my_data = json.loads(flask.request.get_data()) # This is a dictionary.
    structure, name = extract_info(my_data)

    temp_file_folder = MOFSIMPLIFY_PATH + "temp_file_creation_" + str(session['ID']) + '/'

    output = descriptor_generator(name, structure, 'thermal', is_entry) # generate descriptors
    if output == 'FAILED': # Description generation failure
        operation_counter = conditional_diminish(operation_counter)
        if session['permission']:
            db_push(structure=structure, prediction_type='thermal_stability_prediction', in_train='', result='', neighbor_names='', neighbor_dists='', failure=True, csv_content='')
        return 'FAILED'
    elif isinstance(output, dict): # MOF was in the training data. Return the ground truth.
        operation_counter = conditional_diminish(operation_counter)

        return output
    else: # grabbing some variables (that aren't used). We will make a prediction
        temp_file_folder = output[0]
        ANN_folder = output[1]

    # Applying the model next

    timeStarted = time.time() # save start time (debugging)

    prediction, neighbor_names, neighbor_distances = run_thermal_ANN(str(session['ID']), MOFSIMPLIFY_PATH, name, thermal_model)

    timeDelta = time.time() - timeStarted # get execution time
    print('Finished process in ' + str(timeDelta) + ' seconds')

    results = {'prediction': prediction,
        'neighbor_names': neighbor_names, # The names of the five nearest latent space nearest neighbors in the ANN's training data. The ANNs used CoRE MOFs as training data.
        'neighbor_distances': neighbor_distances, # The latent space distances of the five nearest latent space neighbors.
        'in_train': False} # A prediction was made. The requested (selected) MOF was not in the training data.

    print('TIME CHECK 6 ' + str(session['ID']))

    operation_counter = conditional_diminish(operation_counter)
    return results

def db_push(structure, prediction_type, in_train, result, neighbor_names, neighbor_dists, failure, csv_content):
    """
    db_push sends the structure of the MOF being predicted on, and other information, to the MOFSimplify_v2 collection in the MongoDB history database.

    :param structure: str, the text of the cif file of the MOF being analyzed.
    :param prediction_type: str, the type of prediction being run. Can either be 'solvent_stability_prediction' or 'thermal_stability_prediction'.
    :param in_train: boolean, indicates whether the current structure is in the training data or not.
    :param result: float, the ground truth or prediction for this structure.
    :param neighbor_names: list, the latent space nearest neighbor MOFs in the thermal stability ANN.
    :param neighbor_distances: list, the latent space distances of the latent space nearest neighbor MOFs in neighbor_names.
    :param failure: boolean, whether the featurization failed.
    :param csv_content: string, the content of the csv of features.
    :return: A 204 no content response, so the front end does not display a page different than the one it is on.
    """ 

    print(f'Have entered the db_push function')

    client = MongoClient('18.18.63.68',27017) # connect to mongodb
    # The first argument is the IP address. The second argument is the port.
    db = client.history # The history database
    collection = db.MOFSimplify_v2 # The MOFSimplify_v2 collection in the history database.
    # s refers to a solvent removal stability prediction; t refers to a thermal stability prediction
    fields = ['structure', 's_intrain', 's_result', 's_neighbornames', 's_neighbordists', 't_intrain', 't_result', 't_neighbornames', 't_neighbordists', 'failure', 'csv_content', 's_times', 't_times'] 
        # s_times and t_times indicate the number of times this structure has had its solvent removal stability or thermal stability predicted
    final_dict = {}
    for field in fields:
        final_dict[field] = '' # start with everything empty

    # https://docs.mongodb.com/manual/reference/method/db.collection.update/
    my_documents = collection.find({'structure':structure})
    is_entry = my_documents.count() # will be zero or one, depending if structure is in the database

    final_dict['structure'] = structure
    final_dict['failure'] = failure
    final_dict['csv_content'] = csv_content
    if is_entry == 0: # structure has not been seen before. Make a new entry
        # Populate certain fields
        if prediction_type == 'solvent_stability_prediction':
            final_dict['s_intrain'] = in_train
            final_dict['s_result'] = result
            final_dict['s_neighbornames'] = neighbor_names
            final_dict['s_neighbordists'] = neighbor_dists
            final_dict['s_times'] = 1
            final_dict['t_times'] = 0
        elif prediction_type == 'thermal_stability_prediction':
            final_dict['t_intrain'] = in_train
            final_dict['t_result'] = result
            final_dict['t_neighbornames'] = neighbor_names
            final_dict['t_neighbordists'] = neighbor_dists
            final_dict['s_times'] = 0
            final_dict['t_times'] = 1

        # print(final_dict)
        collection.insert(final_dict) # insert the dictionary into the mongodb collection
    else:  # existing structure. Update existing document (couldn't get actual update function to work, so I delete entry and write a new one.)
        my_document = my_documents[0]
        final_dict['s_intrain'] = my_document['s_intrain']
        final_dict['s_result'] = my_document['s_result']
        final_dict['s_neighbornames'] = my_document['s_neighbornames']
        final_dict['s_neighbordists'] = my_document['s_neighbordists']
        final_dict['t_intrain'] = my_document['t_intrain']
        final_dict['t_result'] = my_document['t_result']
        final_dict['t_neighbornames'] = my_document['t_neighbornames']
        final_dict['t_neighbordists'] = my_document['t_neighbordists']
        final_dict['s_times'] = my_document['s_times']
        final_dict['t_times'] = my_document['t_times']

        if prediction_type == 'solvent_stability_prediction': # structure has had a thermal prediction but not a solvent prediction
            final_dict['s_intrain'] = in_train
            final_dict['s_result'] = result
            final_dict['s_neighbornames'] = neighbor_names
            final_dict['s_neighbordists'] = neighbor_dists
            final_dict['s_times'] = final_dict['s_times'] + 1
        elif prediction_type == 'thermal_stability_prediction': # structure has had a solvent prediction but not a thermal prediction
            final_dict['t_intrain'] = in_train
            final_dict['t_result'] = result
            final_dict['t_neighbornames'] = neighbor_names
            final_dict['t_neighbordists'] = neighbor_dists
            final_dict['t_times'] = final_dict['t_times'] + 1

        collection.remove({'structure':structure}) # delete entry
        collection.insert(final_dict) # insert the dictionary into the mongodb collection (so, write new entry)
    return ('', 204) # 204 no content response

def db_push_lite(structure, prediction_type):
    """
    db_push_lite increases either s_times or t_times in the document corresponding to structure, in the MOFSimplify_v2 collection in the MongoDB history database.
    The prediction_type determines whether s_times or t_times is increased by one.

    :param structure: str, the text of the cif file of the MOF being analyzed.
    :param prediction_type: str, the type of prediction being run. Can either be 'solvent_stability_prediction' or 'thermal_stability_prediction'.
    :return: A 204 no content response, so the front end does not display a page different than the one it is on.
    """ 

    print('Entering db_push_lite!')

    client = MongoClient('18.18.63.68',27017) # connect to mongodb
    # The first argument is the IP address. The second argument is the port.
    db = client.history # The history database
    collection = db.MOFSimplify_v2 # The MOFSimplify_v2 collection in the history database.
    final_dict = {}

    # https://docs.mongodb.com/manual/reference/method/db.collection.update/
    my_documents = collection.find({'structure':structure})

    # if db_push_lite is called, structure passed to it will already have a document in the history.MOFSimplify_v2 collection

    # Update existing document (couldn't get actual update function to work, so I delete entry and write a new one.)
    # s refers to a solvent removal stability prediction; t refers to a thermal stability prediction
    # s_times and t_times indicate the number of times this structure has had its solvent removal stability or thermal stability predicted
    my_document = my_documents[0]
    final_dict['structure'] = structure
    final_dict['s_intrain'] = my_document['s_intrain']
    final_dict['s_result'] = my_document['s_result']
    final_dict['s_neighbornames'] = my_document['s_neighbornames']
    final_dict['s_neighbordists'] = my_document['s_neighbordists']
    final_dict['t_intrain'] = my_document['t_intrain']
    final_dict['t_result'] = my_document['t_result']
    final_dict['t_neighbornames'] = my_document['t_neighbornames']
    final_dict['t_neighbordists'] = my_document['t_neighbordists']
    final_dict['s_times'] = my_document['s_times']
    final_dict['t_times'] = my_document['t_times']
    final_dict['failure'] = my_document['failure']
    final_dict['csv_content'] = my_document['csv_content']

    if prediction_type == 'solvent_stability_prediction': 
        final_dict['s_times'] += 1
    elif prediction_type == 'thermal_stability_prediction': 
        final_dict['t_times'] += 1

    collection.remove({'structure':structure}) # delete entry
    collection.insert(final_dict) # insert the dictionary into the mongodb collection (so, write new entry)
    return ('', 204) # 204 no content response

@app.route('/plot_thermal_stability', methods=['POST']) 
def plot_thermal_stability():
    """
    plot_thermal_stability returns a plot of the distribution of thermal breakdown temperatures of the MOFs our ANN was trained on.
    Additionally, it displays the position of the current MOF's predicted or ground truth thermal breakdown temperature.

    :return: HTML document, a plot of the kernel-density estimate of the breakdown temperatures in the training data 
        with a red dot indicating the position of the current MOF's predicted or ground truth thermal breakdown temperature.
    """ 

    # Grab data
    info = json.loads(flask.request.get_data()) 
    my_data = info['temperature'] # this is the current MOF's predicted thermal breakdown temperature
    if info['temperature_data_type'] == 'string':
        my_data = my_data[:-3] # getting rid of the celsius symbol, left with just the number
        my_data = float(my_data)

    # Getting the temperature data
    temps_df = pd.read_csv(MOFSIMPLIFY_PATH + "model/thermal/ANN/adjusted_TSD_df_all.csv")

    import matplotlib
    matplotlib.use('Agg') # noninteractive backend
    import matplotlib.pyplot as plt
    plt.close("all")
    plt.rcParams.update({'font.size': 16}) # large font
    import scipy.stats as stats

    # In training data, smallest T breakdown is 35, and largest T breakdown is 654.

    # Use stats.gaussian_kde to estimate the probability density function from the histogram. This is a kernel-density estimate using Gaussian kernels.
    density = stats.gaussian_kde(temps_df['T'])
    x = np.arange(30,661,1) # in training data, smallest T breakdown is 35, and largest T breakdown is 654
    fig = plt.figure(figsize=(10, 6), dpi=80) # set figure size
    ax = fig.add_subplot(1,1,1)
    plt.plot(x, density(x))
    plt.plot(my_data, density(my_data), "or") # The current MOF's predicted or ground truth thermal breakdown temperature.

    ax.set_xlabel('Breakdown temperature (°C)')
    ax.set_ylabel('Frequency in the training data')

    # The title of the plot differs slightly depending on if the selected MOF is in the training data or not.
    if info['prediction']: # MOF wasn't in training data, and its ANN predicted breakdown temperature is used.
        ax.set_title('Current MOF\'s predicted breakdown temperature relative to others')
    else: # MOF was in the training data, and its reported breakdown temperature is used.
        ax.set_title('Current MOF\'s breakdown temperature relative to others')

    import mpld3

    return mpld3.fig_to_html(fig)

@app.route('/thermal_stability_percentile', methods=['POST']) 
def thermal_stability_percentile():
    """
    thermal_stability_percentile returns what percentile the thermal breakdown temperature (prediction or ground truth) of the selected MOF lies in
        with respect to the MOFs used to train the ANN for thermal stability predictions.
    For example, a 30th percentile breakdown temperature is higher than 30% of the breakdown temperatures in the CoRE training data for the thermal stability ANN.

    :return: str our_percentile, the percentile of the thermal breakdown temperature.
    """ 

    # Grab data.
    my_data = json.loads(flask.request.get_data()) # this is the current MOF's predicted thermal breakdown temperature
    my_data = my_data[:-3] # getting rid of the Celsius symbol, left with just the number
    my_data = float(my_data)

    # Getting the temperature data.
    temps_df = pd.read_csv(MOFSIMPLIFY_PATH + "model/thermal/ANN/adjusted_TSD_df_all.csv")

    # Will find what percentile our prediction belongs to, by checking the 100 percentiles and seeing which is closest to our prediction.
    difference = np.Infinity

    breakdown_Ts = temps_df['T']
    for i in np.arange(0,100.1,1): # 0,1, ..., 99, 100
        current_percentile = np.percentile(breakdown_Ts, i) # ith percentile
        current_difference = np.absolute(my_data - current_percentile) # absolute difference
        if current_difference < difference:
            difference = current_difference
            our_percentile = i

    our_percentile = int(our_percentile) # no decimal points
    our_percentile = str(our_percentile)

    return our_percentile

### Helper functions for TGA_plot. ###

# I use https://stackoverflow.com/questions/3252194/numpy-and-line-intersections to get code for the intersection of the two TGA lines
# for use in /TGA_plot call. perp and seg_intersect are from that stackoverflow page

#
# line segment intersection using vectors
# see Computer Graphics by F.S. Hill
#
def perp( a ) :
    b = np.empty_like(a)
    b[0] = -a[1]
    b[1] = a[0]
    return b

# line segment a given by endpoints a1, a2
# line segment b given by endpoints b1, b2
# return 
def seg_intersect(a1,a2, b1,b2) :
    da = a2-a1
    db = b2-b1
    dp = a1-b1
    dap = perp(da)
    denom = np.dot( dap, db)
    num = np.dot( dap, dp )
    return (num / denom.astype(float))*db + b1

### End of helper functions for TGA_plot. ###

@app.route('/TGA_plot', methods=['POST'])
def TGA_plot():
    """
    TGA_plot makes the TGA plot for the current thermal ANN nearest neighbor.
    The TGA plot contains the four points that constitute the simplified TGA data. These points are yellow stars.
    The two points representing the flatter part of the trace are connected with a line (line A), 
        and the two points representing the steeper part of the trace are connected with a line (line B).
    The intersection of lines A and B is designated with a red dot.
    Lines A and B extend past the red dot for a bit. Lines A and B are blue dashed lines.
    The simplified TGA plot is sent to the front end (index.html) for display.

    :return: HTML document, the simplified TGA plot.
    """ 

    # Grab data
    my_data = json.loads(flask.request.get_data()); # This is the neighbor complex's name.

    # cut off these endings, in order to access the TGA file correctly:
    # _clean, _ion_b, _neutral_b, _SL, _charged, _clean_h, _manual, _auto, _charged, etc
    cut_index = my_data.find('_') # gets index of the first underscore
    my_data = my_data[:cut_index] 

    # Grab data 
    slopes_df = pd.read_csv(MOFSIMPLIFY_PATH + "TGA/raw_TGA_digitization_data/digitized_csv/" + my_data + ".csv")

    x_values = []
    y_values = []
    for i in range(4): # 0, 1, 2, 3
        x_values.append(slopes_df.iloc[[i]]['T (degrees C)'][i])
        y_values.append(slopes_df.iloc[[i]]['mass (arbitrary units)'][i])

    # Making the four points.
    p1 = np.array( [x_values[0], y_values[0]] )
    p2 = np.array( [x_values[1], y_values[1]] )

    p3 = np.array( [x_values[2], y_values[2]] )
    p4 = np.array( [x_values[3], y_values[3]] )

    intersection_point = seg_intersect(p1, p2, p3, p4)

    # Want lines to extend beyond intersection point. I will make it 20%.
    x_extension = [None] * 2 # This variable will hold the x values of the two extension points.
    y_extension = [None] * 2 # This variable will hold the y values of the two extension points.
    slope = [None] * 2 # This variable will hold the slopes of the two lines.

    x_extension[0] = intersection_point[0] + (intersection_point[0] - x_values[0])*0.2 # x_values[0] is the smallest value, x_values[3] is the largest value
    x_extension[1] = intersection_point[0] - (x_values[3] - intersection_point[0])*0.2
    slope[0] = (y_values[1] - y_values[0])/(x_values[1] - x_values[0])
    slope[1] = (y_values[3] - y_values[2])/(x_values[3] - x_values[2])
    y_extension[0] = intersection_point[1] + slope[0] * (x_extension[0] - intersection_point[0]) # y2 = y1 + m(x2-x1)
    y_extension[1] = intersection_point[1] + slope[1] * (x_extension[1] - intersection_point[0]) # y2 = y1 + m(x2-x1)

    # Instantiating the figure object. 
    graph = figure(title = "Simplified literature TGA plot of selected thermal ANN neighbor")  
         
    # The points to be plotted.
    xs = [[x_values[0], x_values[1],intersection_point[0], x_extension[0]], [x_values[2], x_values[3],intersection_point[0], x_extension[1]]] 
    ys = [[y_values[0], y_values[1],intersection_point[1], y_extension[0]], [y_values[2], y_values[3],intersection_point[1], y_extension[1]]] 
        
    # Plotting the graph.
    graph.multi_line(xs, ys, line_dash='dashed') 
    graph.circle([intersection_point[0]], [intersection_point[1]], size=20, fill_color="red", line_color='black')
    graph.star(x=x_values, y=y_values, size=20, fill_color="yellow", line_color='black') # four stars
    graph.xaxis.axis_label = 'Temperature (°C)'
    graph.yaxis.axis_label = 'Percentage mass remaining or Mass'    

    return file_html(graph,CDN,'my plot')

@app.route('/get_components', methods=['POST']) 
def get_components():
    """
    get_components uses Aditya's MOF code to get components (linkers and sbus).
    It returns a dictionary with the linker and sbu xyz files' text, along with information about the number of linkers and sbus.
    The dictionary also contains the SMILES string for each of the linkers and sbus.

    :return: str json_object, encodes a dictionary. The dictionary contains the xyz file geometry information for each component, 
        the number of each type of component, the unique linkers and sbus (kept track of by their index, see linker_num and sbu_num in the code),
        and the SMILES string for each component.

        May instead return a string 'OVERLOAD' if the server is being asked too much of, in order to avoid 50X errors (like 504, aka Gateway Time-out).
    """ 

    global operation_counter # global variable
    global last_operation_counter_clear

    operation_counter_periodic_clear()

    print(f'Current operation_counter: {operation_counter}')

    if operation_counter >= MAX_OPERATIONS:
        return 'OVERLOAD'

    operation_counter += 1

    # Grab data
    my_data = json.loads(flask.request.get_data());
    structure, name = extract_info(my_data)

    temp_file_folder = MOFSIMPLIFY_PATH + "temp_file_creation_" + str(session['ID']) + '/'
    cif_folder = temp_file_folder + 'cifs/'

    # Write the data back to a cif file.
    try:
        cif_file = open(cif_folder + name + '.cif', 'w')
    except FileNotFoundError:
        operation_counter = conditional_diminish(operation_counter)
        return 'FAILED'
    cif_file.write(structure)
    cif_file.close()

    RACs_folder = temp_file_folder +  'components_RACs/'

    # Delete the RACs folder, then remake it (to start fresh for this prediction).
    shutil.rmtree(RACs_folder)
    os.mkdir(RACs_folder)

    # Next, running MOF featurization
    try:
        get_primitive(cif_folder + name + '.cif', cif_folder + name + '_primitive.cif');
    except ValueError:
        operation_counter = conditional_diminish(operation_counter)
        return 'FAILED'

    try:
        full_names, full_descriptors = get_MOF_descriptors(cif_folder + name + '_primitive.cif',3,path= RACs_folder, xyzpath= RACs_folder + name + '.xyz');
            # makes the linkers and sbus folders, complete with linkers and sbus. We don't care about full_names and full_descriptors in this function
    except ValueError:
        operation_counter = conditional_diminish(operation_counter)
        return 'FAILED'
    except NotImplementedError:
        operation_counter = conditional_diminish(operation_counter)
        return 'FAILED'
    except AssertionError:
        operation_counter = conditional_diminish(operation_counter)
        return 'FAILED'

    if (len(full_names) <= 1) and (len(full_descriptors) <= 1): # this is a featurization check from MOF_descriptors.py
        operation_counter = conditional_diminish(operation_counter)
        return 'FAILED'

    # At this point, have the RAC featurization. 

    # will return a json object (dictionary)
    # The fields are string representations of the linkers and sbus, however many there are.

    dictionary = {};
    
    linker_num = 0;
    while True:
        if not os.path.exists(RACs_folder + 'linkers/' + name + '_primitive_linker_' + str(linker_num) + '.xyz'):
            break
        else:
            linker_file = open(RACs_folder + 'linkers/' + name + '_primitive_linker_' + str(linker_num) + '.xyz', 'r');
            linker_info = linker_file.read();
            linker_file.close();

            dictionary['linker_' + str(linker_num)] = linker_info;

            linker_num = linker_num + 1;


    sbu_num = 0;
    while True:
        if not os.path.exists(RACs_folder + 'sbus/' + name + '_primitive_sbu_' + str(sbu_num) + '.xyz'):
            break
        else:
            sbu_file = open(RACs_folder + 'sbus/' + name + '_primitive_sbu_' + str(sbu_num) + '.xyz', 'r');
            sbu_info = sbu_file.read();
            sbu_file.close();

            dictionary['sbu_' + str(sbu_num)] = sbu_info;

            sbu_num = sbu_num + 1;


    dictionary['total_linkers'] = linker_num;
    dictionary['total_sbus'] = sbu_num;


    # Identifying which linkers and sbus have different connectivities.

    # Code below uses molecular graph determinants.
    # two ligands that are the same (via connectivity) will have the same molecular graph determinant. 
    # Molecular graph determinant fails to distinguish between isomers, but so would RCM (Reverse Cuthill McKee).

    # This simple script is meant to be run within the linkers directory, and it will give a bunch of numbers. 
    # If those numbers are the same, the linker is the same, if not, the linkers are different, etc.

    from molSimplify.Classes.mol3D import mol3D
    MOF_of_interest = name + '_primitive'
    XYZs = sorted(glob.glob(RACs_folder + 'linkers/*'+MOF_of_interest+'*xyz'))
    det_list = []

    for xyz in XYZs:
        net = xyz.replace('xyz', 'net') # substring replacement; getting the appropriate .net file
        linker_mol = mol3D()
        linker_mol.readfromxyz(xyz)
        with open(net,'r') as f:
            data = f.readlines()
        for i, line in enumerate(data):
            if i==0:
                # Skip first line
                continue
            elif i == 1:
                graph = np.array([[int(val) for val in line.strip('\n').split(',')]])
            else:
                graph = np.append(graph, [np.array([int(val) for val in line.strip('\n').split(',')])],axis=0)
        linker_mol.graph = graph
        safedet = linker_mol.get_mol_graph_det(oct=False)
        det_list.append(safedet)

    #### Linkers with the same molecular graph determinant are the same.
    #### Molecular graph determinant does not catch isomers.
    linker_det_list = det_list

    unique_linker_det = set(linker_det_list) # getting the unique determinants
    unique_linker_indices = []
    for item in unique_linker_det:
        unique_linker_indices.append(linker_det_list.index(item)) # indices of the unique linkers in the list of linkers

    XYZs = sorted(glob.glob(RACs_folder + 'sbus/*'+MOF_of_interest+'*xyz'))
    det_list = []
    for xyz in XYZs:
        net = xyz.replace('xyz', 'net') # substring replacement; getting the appropriate .net file
        linker_mol = mol3D()
        linker_mol.readfromxyz(xyz)
        with open(net,'r') as f:
            data = f.readlines()
        for i, line in enumerate(data):
            if i==0:
                # Skip first line
                continue
            elif i == 1:
                graph = np.array([[int(val) for val in line.strip('\n').split(',')]])
            else:
                graph = np.append(graph, [np.array([int(val) for val in line.strip('\n').split(',')])],axis=0)
        linker_mol.graph = graph
        safedet = linker_mol.get_mol_graph_det(oct=False)
        det_list.append(safedet)

    sbu_det_list = det_list

    unique_sbu_det = set(sbu_det_list) # getting the unique determinants
    unique_sbu_indices = []
    for item in unique_sbu_det:
        unique_sbu_indices.append(sbu_det_list.index(item)) # indices of the unique sbus in the list of linkers

    # adding the unique indices to the dictionary
    dictionary['unique_linker_indices'] = unique_linker_indices
    dictionary['unique_sbu_indices'] = unique_sbu_indices


    # In this next section, getting the SMILES strings for all of the linkers and sbus using pybel.
    # Write the smiles strings to file, then read the file.
    import pybel

    linkers_folder = RACs_folder + 'linkers/'

    for i in range(dictionary['total_linkers']): # 0, 1, 2, ..., numberoflinkersminus1
        smilesFile = pybel.Outputfile('smi', linkers_folder + name + '_primitive_linker_' + str(i) + '.txt') # smi refers to SMILES
        smilesFile.write(next(pybel.readfile('xyz', linkers_folder + name + '_primitive_linker_' + str(i) + '.xyz'))) # writes SMILES string to the text file

        # Next, get the SMILES string from the text file.
        f = open(linkers_folder + name + '_primitive_linker_' + str(i) + '.txt', 'r')
        line = f.readline()
        line = line.split('\t') # split at tabs
        smiles_ID = line[0]
        dictionary['linker_' + str(i) + '_SMILES'] = smiles_ID
        f.close()

    sbus_folder = RACs_folder + 'sbus/'

    for i in range(dictionary['total_sbus']):
        smilesFile = pybel.Outputfile('smi', sbus_folder + name + '_primitive_sbu_' + str(i) + '.txt') # smi refers to SMILES
        smilesFile.write(next(pybel.readfile('xyz', sbus_folder + name + '_primitive_sbu_' + str(i) + '.xyz'))) # writes SMILES string to the text file

        # Next, get the SMILES string from the text file.
        f = open(sbus_folder + name + '_primitive_sbu_' + str(i) + '.txt', 'r')
        line = f.readline()
        line = line.split('\t') # split at tabs
        smiles_ID = line[0]
        dictionary['sbu_' + str(i) + '_SMILES'] = smiles_ID
        f.close()


    json_object = json.dumps(dictionary, indent = 4);

    operation_counter = conditional_diminish(operation_counter)
    return json_object

@app.route('/solvent_neighbor_flag', methods=['POST']) 
def is_stable():
    """
    is_stable returns the flag (i.e. the stability upon solvent removal, yes or no) and DOI of the neighbor sent over from the front end.
    The flag is returned to the front end (index.html). 

    :return: dict myDict, contains the solvent removal stability of the latent space nearest neighbor MOF, and its DOI.
    """ 

    # Grab data.
    my_data = json.loads(flask.request.get_data()); # This is the neighbor complex's name.

    print(my_data)
    
    my_folder = MOFSIMPLIFY_PATH + 'model/solvent/ANN/dropped_connectivity_dupes/'
    solvent_flags_df = pd.read_csv(my_folder + 'train.csv')

    this_neighbor = solvent_flags_df[solvent_flags_df['CoRE_name'] == my_data] # getting the row with the MOF of interest

    this_neighbor_flag = this_neighbor['flag'] # getting the flag value
    this_neighbor_flag = str(this_neighbor_flag.iloc[0]) # extract the flag value and return it as a string

    # next, getting DOI
    this_neighbor_doi = this_neighbor['doi'] # getting the doi
    this_neighbor_doi = this_neighbor_doi.iloc[0]

    myDict = {'flag': this_neighbor_flag, 'doi': this_neighbor_doi}

    return myDict

@app.route('/thermal_neighbor_T', methods=['POST']) 
def breakdown_T():
    """
    breakdown_T returns the thermal breakdown temperature and DOI of the neighbor sent over from the front end.
    The information is returned to the front end (index.html).
    
    :return: dict myDict, contains the breakdown temperature of the latent space nearest neighbor MOF, and its DOI.
    """ 

    # Grab data
    my_data = json.loads(flask.request.get_data()); # This is the neighbor complex's name.

    print(my_data)
    
    # To begin, always go to main directory 
    ANN_folder = MOFSIMPLIFY_PATH + 'model/thermal/ANN/'

    breakdown_T_df = pd.read_csv(ANN_folder + 'train.csv')
    this_neighbor = breakdown_T_df[breakdown_T_df['CoRE_name'] == my_data] # getting the row with the MOF of interest
    this_neighbor_T = this_neighbor['T'] # getting the breakdown temperature value
    this_neighbor_T = str(round(this_neighbor_T.iloc[0], 1)) # extract the breakdown temperature value and return it as a string. Want just one decimal place


    TGA_folder = MOFSIMPLIFY_PATH + 'TGA/'
    TGA_df = pd.read_excel(TGA_folder + 'TGA_info_log.xlsx')
    this_neighbor = TGA_df[TGA_df['CoRE_name'] == my_data] # getting the row with the MOF of interest
    this_neighbor_doi = this_neighbor['doi'] # getting the doi
    this_neighbor_doi = this_neighbor_doi.iloc[0]

    myDict = {'T': this_neighbor_T, 'doi': this_neighbor_doi}

    return myDict

@app.route('/neighbor_writer', methods=['POST']) 
def neighbor_writer():
    """
    neighbor_writer writes information to a txt file about the selected latent space nearest neighbor.
    It returns the information in that txt file to the front end for downloading.

    :return: dict myDict, contains the content for the text file to be downloaded, and the name of the text file to be downloaded.
    """ 

    # Grab data
    my_data = json.loads(flask.request.get_data()); # This is a dictionary with information about the neighbor and the MOF that was analyzed by an ANN

    prediction_type = my_data['prediction_type']
    current_MOF = my_data['current_MOF']
    selected_neighbor = my_data['selected_neighbor']
    latent_space_distance = my_data['latent_space_distance']

    # next, getting neighbor_truth
    if prediction_type == 'solvent':
        my_folder = MOFSIMPLIFY_PATH + 'model/solvent/ANN/dropped_connectivity_dupes/'
        solvent_flags_df = pd.read_csv(my_folder + 'train.csv')
        this_neighbor = solvent_flags_df[solvent_flags_df['CoRE_name'] == selected_neighbor] # getting the row with the MOF of interest
        this_neighbor_flag = this_neighbor['flag'] # getting the flag value
        neighbor_truth = str(this_neighbor_flag.iloc[0]) # extract the flag value and convert it into a string
        if neighbor_truth == '1':
            neighbor_truth = 'Stable upon solvent removal'
        else:
            neighbor_truth = 'Unstable upon solvent removal'

    else: # thermal stability
        my_folder = MOFSIMPLIFY_PATH + 'model/thermal/ANN/'
        breakdown_T_df = pd.read_csv(my_folder + 'train.csv')
        this_neighbor = breakdown_T_df[breakdown_T_df['CoRE_name'] == selected_neighbor] # getting the row with the MOF of interest
        this_neighbor_T = this_neighbor['T'] # getting the breakdown temperature value
        neighbor_truth = str(round(this_neighbor_T.iloc[0], 1)) # extract the flag value and convert it into a string. Want just one decimal place


    # next, getting DOI
    if prediction_type == 'solvent':
        this_neighbor = solvent_flags_df[solvent_flags_df['CoRE_name'] == selected_neighbor] # getting the row with the MOF of interest
        this_neighbor_doi = this_neighbor['doi'] # getting the doi
        this_neighbor_doi = this_neighbor_doi.iloc[0]
    else: # thermal
        TGA_df = pd.read_excel(MOFSIMPLIFY_PATH + 'TGA/TGA_info_log.xlsx')
        this_neighbor = TGA_df[TGA_df['CoRE_name'] == selected_neighbor] # getting the row with the MOF of interest
        this_neighbor_doi = this_neighbor['doi'] # getting the doi
        this_neighbor_doi = this_neighbor_doi.iloc[0]

    # Now, writing all of this information to a string

    # Delete and remake the latent_neighbor folder each time, so only one file is ever in it
    shutil.rmtree(MOFSIMPLIFY_PATH + 'temp_file_creation_' + str(session['ID']) + '/latent_neighbor')
    os.mkdir(MOFSIMPLIFY_PATH + 'temp_file_creation_' + str(session['ID']) + '/latent_neighbor')

    with open(MOFSIMPLIFY_PATH + 'temp_file_creation_' + str(session['ID']) + '/latent_neighbor/' + prediction_type + '-' + current_MOF + '-' + selected_neighbor + '.txt','w') as f:
        f.write('Prediction type: ' + prediction_type + '\n')
        f.write('Current MOF: ' + current_MOF + '\n')
        f.write('Selected CoRE nearest neighbor in latent space: ' + selected_neighbor + '\n')
        f.write('Latent space distance: ' + latent_space_distance + '\n')
        if prediction_type == 'solvent':
            f.write('Property for nearest neighbor: ' + neighbor_truth + '\n')
        else:
            degree_sign= u'\N{DEGREE SIGN}'
            f.write('Property for nearest neighbor: ' + neighbor_truth + ' ' + degree_sign + 'C\n')
        f.write('Neighbor DOI: ' + this_neighbor_doi + '\n')

    with open(MOFSIMPLIFY_PATH + 'temp_file_creation_' + str(session['ID']) + '/latent_neighbor/' + prediction_type + '-' + current_MOF + '-' + selected_neighbor + '.txt','r') as f:
        contents = f.read()

    myDict = {'contents': contents, 'file_name': prediction_type + '-' + current_MOF + '-' + selected_neighbor + '.txt'}

    return myDict

@app.route('/get_descriptors', methods=['POST']) 
def descriptor_getter():
    """
    descriptor_getter returns the contents of the csv with the descriptors of the desired MOF.
    These descriptors are RACs and Zeo++ descriptors. There are 189 total.

    :return: str contents, the contents of the csv with the descriptor data. To be written to a csv in the front end for download.
    """ 

    # Grab data
    name = json.loads(flask.request.get_data()); # This is the selected MOF
    if name[-4:] == '.cif':
        name = name[:-4] # remove the .cif part of the name

    temp_file_folder = MOFSIMPLIFY_PATH + "temp_file_creation_" + str(session['ID']) + '/'
    descriptors_folder = temp_file_folder + "merged_descriptors/"

    with open(descriptors_folder + name + '_descriptors.csv', 'r') as f:
        contents = f.read()

    return contents

@app.route('/get_TGA', methods=['POST']) 
def TGA_getter():
    """
    TGA_getter returns the contents of the csv with the simplified TGA data of the desired MOF.

    :return: str contents, the contents of the csv with the simplified TGA data. To be written to a csv in the front end for download.
    """ 

    # Grab data
    name = json.loads(flask.request.get_data()); # This is the thermal ANN latent space nearest neighbor
    
    # cut off these endings, in order to access the TGA file correctly:
    # _clean, _ion_b, _neutral_b, _SL, _charged, _clean_h, _manual, _auto, _charged, etc
    cut_index = name.find('_') # gets index of the first underscore
    name = name[:cut_index] 

    tga_folder = MOFSIMPLIFY_PATH + "TGA/raw_TGA_digitization_data/digitized_csv/"

    with open(tga_folder + name + '.csv', 'r') as f:
        contents = f.read()

    return contents

@app.route('/bb', methods=['POST'])
def stable_bb_getter():
    """
    stable_bb_getter returns the contents of the XYZ of the building block specified by the front end.

    :return: str contents, the contents of the XYZ file of the requested building block.
    """

    # Grab data
    name = json.loads(flask.request.get_data()); # This is the building block to get.

    bb_file = open(MOFSIMPLIFY_PATH + 'stable_building_blocks/building_blocks/' + name + '.xyz', 'r');
    bb_info = bb_file.read();
    bb_file.close();

    return bb_info

@app.route('/filter_MOFs', methods=['POST'])
def filter_MOFs():
    """
    filter_MOFs returns the list of stable MOFs after applying the filter criteria passed from the front end.

    :return: list contents, the stable MOFs that survive all of the filters (e.g. inorganic node filter).
    """

    # Grab data
    filter_dictionary = json.loads(flask.request.get_data()); # This dictionary contains the filters to apply.
    print(f'check, filter_dictionary is {filter_dictionary}')

    ultrastable_filter = filter_dictionary['ultrastable']
    inorg_node_filter = filter_dictionary['inorg_node']
    org_node_filter = filter_dictionary['org_node']
    edge_filter = filter_dictionary['edge']
    metal_filter = filter_dictionary['metal']
    type_filter = filter_dictionary['type']
    net_filter = filter_dictionary['net']

    # Applying the type and ultrastable filter filter.
    if ultrastable_filter == 'No':
        if type_filter == 'Any' or len(type_filter) == 0: # # No filter applied in the 'Any' case. The empty string occurs if the user clicks on the selection box, and then clicks off without selecting anything
            MOF_list = list_1inorganic_1edge_MOFs + list_1inorganic_1organic_1edge_MOFs + list_2inorganic_1edge_MOFs
        elif type_filter == '1, 0, 1':
            MOF_list = list_1inorganic_1edge_MOFs
        elif type_filter == '1, 1, 1':
            MOF_list = list_1inorganic_1organic_1edge_MOFs
        elif type_filter == '2, 0, 1':
            MOF_list = list_2inorganic_1edge_MOFs
    elif ultrastable_filter == 'Yes':
        if type_filter == 'Any' or len(type_filter) == 0: # # No filter applied in the 'Any' case. The empty string occurs if the user clicks on the selection box, and then clicks off without selecting anything
            MOF_list = list_1inorganic_1edge_MOFs_ultrastable + list_1inorganic_1organic_1edge_MOFs_ultrastable + list_2inorganic_1edge_MOFs_ultrastable
        elif type_filter == '1, 0, 1':
            MOF_list = list_1inorganic_1edge_MOFs_ultrastable
        elif type_filter == '1, 1, 1':
            MOF_list = list_1inorganic_1organic_1edge_MOFs_ultrastable
        elif type_filter == '2, 0, 1':
            MOF_list = list_2inorganic_1edge_MOFs_ultrastable
    else:
        raise Exception('Unexpected value of ultrastable_filter')

    filtered_MOFs = MOF_list.copy()

    # Getting rid of '.cif' from each string. For nicer display on the front end (MOF name without .cif).
    filtered_MOFs = [i.replace('.cif', '') for i in filtered_MOFs]

    # Dictionary that indicates which SBUs contain which metals.
    metal_dict = {'Gd': ['PORLAL_clean_sbu_0', 'MAKGOX_clean_sbu_0', 'INOVEN_clean_sbu_0', 'KUMBOL_clean_sbu_0', 'KUMBOL_clean_sbu_1'], 'Tb': ['MAKGUD_clean_sbu_1', 'KUMBUR_clean_sbu_0', 'KUMBUR_clean_sbu_1', 'NUHQIS_clean_sbu_0', 'NUHQUE_clean_sbu_2', 'ZALLEG_clean_sbu_0', 'NUHRAL_clean_sbu_3', 'CUKXOW_clean_sbu_0', 'XUDNAN_clean_sbu_1', 'BETDIP_clean_sbu_0', 'XADDAJ01_clean_sbu_3', 'PAMTUU_clean_sbu_0'], 'Mn': ['BELTOD_clean_sbu_0', 'KUFVIS_clean_sbu_0', 'CIFDUS_clean_sbu_1', 'ICIZOL_clean_sbu_0'], 'Dy': ['BETDEL_clean_sbu_0', 'PAMTOO_clean_sbu_0'], 'Co': ['CEKHIL_clean_sbu_0', 'JUFBIX_clean_sbu_0', 'GIZVER_clean_sbu_0', 'BEQFEK_clean_sbu_0', 'OLANAS_clean_sbu_2', 'EBIMEJ_clean_sbu_0', 'KOZSID_clean_sbu_0', 'XOMCOT_clean_sbu_1', 'OLANEW_clean_sbu_3', 'XOMCOT_clean_sbu_0', 'ICAMEG_clean_sbu_1', 'OLANEW_clean_sbu_0', 'LIZSOE_clean_sbu_0'], 'Nd': ['GEDQOX_clean_sbu_0'], 'U': ['WUSLED_clean_sbu_0'], 'Zn': ['QEWDON_clean_sbu_0', 'NAHDIM_clean_sbu_0', 'NAWXER_clean_sbu_0', 'NAHDIM_clean_sbu_2', 'FANWIC_clean_sbu_0', 'HICVOG_clean_sbu_0', 'TATPOV_clean_sbu_0', 'VOLQOD_clean_sbu_1', 'MIDCAF_clean_sbu_0', 'OLANAS_clean_sbu_1', 'FUNLAD_clean_sbu_0', 'UKALOJ_clean_sbu_0', 'MIFKUJ_clean_sbu_0'], 'Sm': ['BETFAJ_clean_sbu_0'], 'Eu': ['HISSIN_clean_sbu_1', 'HISSIN_clean_sbu_0', 'NUHQIS_clean_sbu_0', 'KUMJIN_clean_sbu_1', 'NUHQUE_clean_sbu_2', 'CUWYAW_clean_sbu_0', 'KUMJIN_clean_sbu_0', 'EYACOX_clean_sbu_1', 'EYACOX_clean_sbu_0', 'NUHRAL_clean_sbu_3', 'SARMOO_clean_sbu_0', 'IYUCEM_clean_sbu_0', 'EZIPEK_clean_sbu_0'], 'Mg': ['TAGTUT_clean_sbu_1', 'AZAVOO_clean_sbu_0', 'TAGTUT_clean_sbu_3', 'WAQDOJ_charged_sbu_0', 'LEVNOQ01_clean_sbu_1', 'EQERAU_clean_sbu_0', 'MUDLON_clean_sbu_0'], 'La': ['BETGAK_clean_sbu_0'], 'Fe': ['VETTIZ_charged_sbu_0', 'IZENUY_clean_sbu_0'], 'Na': ['OLANAS_clean_sbu_1'], 'Cd': ['AJUNOK_clean_sbu_0', 'OLANEW_clean_sbu_2', 'OLANEW_clean_sbu_3', 'QUQFOY_clean_sbu_3', 'QUQFOY_clean_sbu_1', 'OLANEW_clean_sbu_1', 'OLANEW_clean_sbu_0', 'OLANEW_clean_sbu_4'], 'Pr': ['BETFEN_clean_sbu_0'], 'Zr': ['ENOWUB_clean_sbu_0'], 'Hf': ['UKIBUN_clean_sbu_0'], 'Ho': ['BETDAH_clean_sbu_0'], 'Cu': ['BOTCEU_clean_sbu_0', 'ESEHIV_clean_sbu_0'], 'Sc': ['OJICUG_clean_sbu_0'], 'Sr': ['KAKCAD_clean_sbu_0'], 'Li': ['UXUYUI_clean_sbu_0'], 'In': ['GALJAG_clean_sbu_0']}

    # Applying the rest of the filters.
    if inorg_node_filter != 'Any' and len(inorg_node_filter) != 0: # No filter applied in the 'Any' case. The empty string occurs if the user clicks on the selection box, and then clicks off without selecting anything   
        filtered_MOFs = [i for i in filtered_MOFs if (bb_mapping[inorg_node_filter]+'_') in i] # The underscore helps distinguish between N4 and N46, for example
    if org_node_filter != 'Any' and len(org_node_filter) != 0:
        filtered_MOFs = [i for i in filtered_MOFs if (bb_mapping[org_node_filter]+'_') in i] # The underscore helps distinguish between N4 and N46, for example
    if edge_filter != 'Any' and len(edge_filter) != 0:
        if edge_filter == 'none':
            filtered_MOFs = [i for i in filtered_MOFs if i.endswith('none')]
        else: # In this case, we need to de-alias using the mapping.
            filtered_MOFs = [i for i in filtered_MOFs if i.endswith(bb_mapping[edge_filter])]
    if metal_filter != 'Any' and len(metal_filter) != 0:
        suitable_SBUs = metal_dict[metal_filter]
        suitable_SBU_aliases = [bb_mapping[SBU] for SBU in suitable_SBUs]

        filtered_MOFs = [i for i in filtered_MOFs if any(j in i for j in suitable_SBU_aliases)] # See if any of the SBUs with the requested metal in it are in the assessed MOF i. If so, keep.
    if net_filter != 'Any' and len(net_filter) != 0:
        filtered_MOFs = [i for i in filtered_MOFs if f'net-{net_filter}_node1' in i]

    response_dict = {'filtered_MOFs': filtered_MOFs} # List is not a valid return type, but dictionary is.

    return response_dict

def type_determination(MOF_name):
    """
    type_determination returns a string indicating the type of the MOF passed from the front end. Either 1,0,1; 1,1,1; or 2,0,1 for #inorganic nodes, #organic nodes, #edges 

    :param MOF_name: str, the name of the MOF for which the type will be ascertained.
    :return: str MOF_type, the type of the MOF in question.
    """
    if MOF_name in list_1inorganic_1edge_MOFs:
        MOF_type = '1inorganic_1edge'
    elif MOF_name in list_1inorganic_1organic_1edge_MOFs:
        MOF_type = '1inorganic_1organic_1edge'
    elif MOF_name in list_2inorganic_1edge_MOFs:
        MOF_type = '2inorganic_1edge'
    else:
        raise Exception('CIF name is unexpected. Not in any of the MOF lists.')

    return MOF_type

@app.route('/data_pull', methods=['POST'])
def grab_data():
    """
    grab_data returns a dictionary of information on the MOF whose name is sent from the front end.

    :return: dict data_dict, a dictionary with stability and uptake information on the selected MOF.
    """

    # Grab data
    name = json.loads(flask.request.get_data()); # This is the MOF about which to gather data.    
    MOF_type = type_determination(name + '.cif')
    if MOF_type == '1inorganic_1edge':
        type_abbrev = '1inor_1edge'
    elif MOF_type == '1inorganic_1organic_1edge':
        type_abbrev = '1inor_1org_1edge'
    elif MOF_type == '2inorganic_1edge':
        type_abbrev = '2inor_1edge'

    # First, getting ANN solvent removal and thermal stabilities
    df = pd.read_csv(f'stable_building_blocks/predictions/{type_abbrev}_predictions.csv')
    adjusted_name = name.replace('optimized_', '') + '.cif' # To match the format of entries in the predictions excel files.

    try:
        desired_row = df[df['filename'] == adjusted_name]
        predicted_solvent_removal_stability = desired_row['predicted_solvent_removal_stability'].iloc[0]
        predicted_thermal_stability = desired_row['predicted_thermal_stability'].iloc[0]

        # Reducing decimal places. 2 for solvent removal stability, 1 for thermal stability.
        predicted_solvent_removal_stability = predicted_solvent_removal_stability.round(2)
        predicted_thermal_stability = predicted_thermal_stability.round(1)

    except IndexError: # adjusted_name is not present in the CSV
        predicted_solvent_removal_stability = None
        predicted_thermal_stability = None

    # Next, getting mechanical moduli
    df = pd.read_csv(f'stable_building_blocks/computed_properties/stable/moduli_{type_abbrev}.csv')
    try:
        desired_row = df[df['name'] == 'optimized_'+name] # Don't need to use adjusted_name, since the moduli CSV files have a slightly different format for the name column.
        moduli_info = {} # Dictionary to be populated.
        moduli_properties = ['KR', 'KV', 'KVRH', 'GR', 'GV', 'GVRH']
        for prop in moduli_properties:
            moduli_info[prop] = desired_row[prop].iloc[0].round(1) # Round to one decimal place.
    except IndexError: # adjusted_name is not present in the CSV
        moduli_info = None

    # Next, getting uptake information
    df = pd.read_csv(f'stable_building_blocks/computed_properties/stable/uptake_{type_abbrev}.csv')
    try:
        desired_row = df[df['filename'] == adjusted_name]
        uptake_info = {} # Dictionary to be populated.
        uptake_properties = ['uptake_low_p (cm3 methane/g framework)', 'uptake_low_p (cm3 methane/cm3 framework)', 'uptake_high_p (cm3 methane/g framework)', 'uptake_high_p (cm3 methane/cm3 framework)']
        for prop in uptake_properties:
            uptake_info[prop] = desired_row[prop].iloc[0].round(2) # Round to two decimal places.
    except IndexError: # adjusted_name is not present in the CSV
        uptake_info = None

    # Lastly, assessing whether MOF is ultrastable
    is_ultrastable = adjusted_name in ultrastable_MOFs

    data_dict = {
    'solvent_removal_stability': predicted_solvent_removal_stability, 
    'thermal_stability': predicted_thermal_stability,
    'moduli_info': moduli_info,
    'uptake_info': uptake_info,
    'is_ultrastable': is_ultrastable
    }
    
    print(data_dict)

    return data_dict

@app.route('/opt_geo_existence_check', methods=['POST'])
def existence_check():
    """
    existence_check returns a string indicating whether the optimized geometry for the requested MOF exists.

    :return: str exists, a string indicating whether the optimized CIF file is available for download.
    """

    # Grab data
    name = json.loads(flask.request.get_data()); # This is the MOF about which to gather data.    
    MOF_type = type_determination(name + '.cif')

    exists = os.path.isfile(f'stable_building_blocks/optimized_structures/{MOF_type}/optimized_{name}.cif')
    if exists:
        exists = 'Yes'
    else:
        exists = 'No'

    return exists

### New section. The code for water and acid stability. ### 

def db_push_water(structure, prediction_type, result, failure, csv_content):
    """
    db_push_water sends the structure of the MOF being predicted on, and other information, to the MOFSimplify_water_v2 collection in the MongoDB history database.

    :param structure: str, the text of the cif file of the MOF being analyzed.
    :param prediction_type: str, the type of prediction being run. Can either be 'water_stability_prediction' or 'acid_stability_prediction'.
    :param result: float, the prediction for this structure.
    :param failure: boolean, whether the featurization failed.
    :param csv_content: string, the content of the csv of features.
    :return: A 204 no content response, so the front end does not display a page different than the one it is on.
    """
    client = MongoClient('18.18.63.68',27017) # connect to mongodb
    # The first argument is the IP address. The second argument is the port.
    db = client.history # The history database
    collection = db.MOFSimplify_water_v2 # The MOFSimplify_water_v2 collection in the history database.
    # w refers to a 2-class water prediction; a refers to an acid prediction
    fields = ['structure', 'w_result', 'a_result', 'failure', 'csv_content', 'w_times', 'a_times'] 
        # w_times and a_times indicate the number of times this structure has had its water stability or acid stability predicted
    final_dict = {}
    for field in fields:
        final_dict[field] = '' # start with everything empty

    # https://docs.mongodb.com/manual/reference/method/db.collection.update/
    my_documents = collection.find({'structure':structure})
    is_entry = my_documents.count() # will be zero or one, depending if structure is in the database

    final_dict['structure'] = structure
    final_dict['failure'] = failure
    final_dict['csv_content'] = csv_content
    if is_entry == 0: # structure has not been seen before. Make a new entry
        # Populate certain fields
        if prediction_type == 'water_stability_prediction':
            final_dict['w_result'] = result
            final_dict['w_times'] = 1
            final_dict['a_times'] = 0
        elif prediction_type == 'acid_stability_prediction':
            final_dict['a_result'] = result
            final_dict['w_times'] = 0
            final_dict['a_times'] = 1

        collection.insert(final_dict) # insert the dictionary into the mongodb collection
    else:  # existing structure. Update existing document (couldn't get actual update function to work, so I delete entry and write a new one.)
        my_document = my_documents[0]
        final_dict['w_result'] = my_document['w_result']
        final_dict['a_result'] = my_document['a_result']
        final_dict['w_times'] = my_document['w_times']
        final_dict['a_times'] = my_document['a_times']

        if prediction_type == 'water_stability_prediction': # structure has had an acid prediction but not a water prediction
            final_dict['w_result'] = result
            final_dict['w_times'] = final_dict['w_times'] + 1
        elif prediction_type == 'acid_stability_prediction': # structure has had a water prediction but not an acid prediction
            final_dict['a_result'] = result
            final_dict['a_times'] = final_dict['a_times'] + 1

        collection.remove({'structure':structure}) # delete entry
        collection.insert(final_dict) # insert the dictionary into the mongodb collection (so, write new entry)
    return ('', 204) # 204 no content response

def db_push_lite_water(structure, prediction_type):
    """
    # db_push_lite_water increases either w_times or a_times in the document corresponding to structure, in the MOFSimplify_water_v2 collection in the MongoDB history database.
    # The prediction_type determines whether w_times or a_times is increased by one.

    :param structure: str, the text of the cif file of the MOF being analyzed.
    :param prediction_type: str, the type of prediction being run. Can either be 'water_stability_prediction' or 'acid_stability_prediction'.
    :return: A 204 no content response, so the front end does not display a page different than the one it is on.
    """ 

    client = MongoClient('18.18.63.68',27017) # connect to mongodb
    # The first argument is the IP address. The second argument is the port.
    db = client.history # The history database
    collection = db.MOFSimplify_water_v2 # The MOFSimplify_water_v2 collection in the history database.
    final_dict = {}

    # https://docs.mongodb.com/manual/reference/method/db.collection.update/
    my_documents = collection.find({'structure':structure})

    # if db_push_lite_water is called, structure passed to it will already have a document in the history.MOFSimplify_water_v2 collection

    # Update existing document (couldn't get actual update function to work, so I delete entry and write a new one.)
    # w refers to a 2-class water stability prediction; a refers to an acid stability prediction
    # w_times and a_times indicate the number of times this structure has had its water stability or acid stability predicted
    my_document = my_documents[0]
    final_dict['structure'] = structure
    final_dict['w_result'] = my_document['w_result']
    final_dict['a_result'] = my_document['a_result']
    final_dict['w_times'] = my_document['w_times']
    final_dict['a_times'] = my_document['a_times']
    final_dict['failure'] = my_document['failure']
    final_dict['csv_content'] = my_document['csv_content']

    if prediction_type == 'water_stability_prediction': 
        final_dict['w_times'] += 1
    elif prediction_type == 'acid_stability_prediction': 
        final_dict['a_times'] += 1

    collection.remove({'structure':structure}) # delete entry
    collection.insert(final_dict) # insert the dictionary into the mongodb collection (so, write new entry)
    return ('', 204) # 204 no content response

@app.route('/predict_water_stability', methods=['POST']) 
def ws_predict():
    """
    ws_predict generates the 2-class water stability prediction for the selected MOF.
        Or it will return a signal that descriptor generation failed and thus a prediction cannot be made. 
    RAC featurization and Zeo++ geometry information for the selected MOF is generated, using descriptor_generator_water.
    Then, Gianmarco's model is applied to make a prediction using run_water_RF.

    :return: dict results, contains the prediction.
        May instead return a string 'FAILED' if descriptor generation fails.
        May instead return a string 'OVERLOAD' if the server is being asked too much of, in order to avoid 50X errors (like 504, aka Gateway Time-out).
    """ 

    global operation_counter # global variable
    global last_operation_counter_clear

    operation_counter_periodic_clear()

    if operation_counter >= MAX_OPERATIONS:
        return 'OVERLOAD'

    operation_counter += 1

    # Grab data and tweak the MOF name.
    my_data = json.loads(flask.request.get_data())
    structure, name = extract_info(my_data)

    temp_file_folder = MOFSIMPLIFY_PATH + "temp_file_creation_" + str(session['ID']) + '/'

    ### Check in MongoDB history.MOFSimplify_water_v2 collection to see if this structure has been predicted on before. 
    # Comment out this section if you are running MOFSimplify on your computer, and define is_entry as False. Also select "No" for the question May MOFSimplify store information on your MOFs?
    client = MongoClient('18.18.63.68',27017) # connect to mongodb. The first argument is the IP address. The second argument is the port.
    db = client.history # The history database
    collection = db.MOFSimplify_water_v2 # The MOFSimplify_water_v2 collection in the history database.
    my_documents = collection.find({'structure':structure})
    is_entry = my_documents.count() # will be zero or one, depending if structure is in the database
    if is_entry: # return statements in this if statement are of same form as those from the workflow if the structure isn't in the database. The frontend needs to receive information of the same form.
        entry_data = my_documents[0]
        if entry_data['failure'] == True: # MOF was not featurizable
            if session['permission']:
                db_push_lite_water(structure=structure, prediction_type='water_stability_prediction') # add 1 to the w_times in the database
            operation_counter = conditional_diminish(operation_counter)
            return 'FAILED'
        else: 
            # Making the csv; can skip csv making in descriptor_generator with this
            csv_content = entry_data['csv_content']
            with open(temp_file_folder + 'merged_descriptors/' + name + '_descriptors_water.csv', 'w') as f:
                f.write(csv_content)

        if entry_data['w_times'] > 0:
            my_dict = {'prediction':entry_data['w_result']}
            if session['permission']:
                db_push_lite_water(structure=structure, prediction_type='water_stability_prediction') # add 1 to the w_times in the database
            operation_counter = conditional_diminish(operation_counter)
            return my_dict

        # if haven't returned anything by now and entered the if is_entry statement, have a featurizable MOF for which a water prediction has not been run, but an acid prediction has

    output = descriptor_generator_water(name, structure, 'water', is_entry) # generate descriptors
    if output == 'FAILED': # Description generation failure
        operation_counter = conditional_diminish(operation_counter)
        if session['permission']:
            db_push_water(structure=structure, prediction_type='water_stability_prediction', result='', failure=True, csv_content='')
        return 'FAILED'

    # Applying the model next
    prediction = run_water_RF(str(session['ID']), MOFSIMPLIFY_PATH, name)

    results = {
    'prediction': prediction,
    }
    
    if session['permission']:
        # Database push of MOF structure and ML results.
        # Getting the contents of the CSV containing features.
        descriptors_folder = temp_file_folder + "merged_descriptors/"
        with open(descriptors_folder + name + '_descriptors_water.csv', 'r') as f:
            csv_contents = f.read()
        db_push_water(structure=structure, prediction_type='water_stability_prediction', result=prediction, failure=False, csv_content=csv_contents)

    operation_counter = conditional_diminish(operation_counter)
    return results

@app.route('/predict_acid_stability', methods=['POST']) 
def as_predict():
    """
    as_predict generates the acid stability prediction for the selected MOF.
        Or it will return a signal that descriptor generation failed and thus a prediction cannot be made. 
    RAC featurization and Zeo++ geometry information for the selected MOF is generated, using descriptor_generator_water.
    Then, Gianmarco's model is applied to make a prediction using run_acid_RF.

    :return: dict results, contains the prediction.
        May instead return a string 'FAILED' if descriptor generation fails.
        May instead return a string 'OVERLOAD' if the server is being asked too much of, in order to avoid 50X errors (like 504, aka Gateway Time-out).
    """ 

    global operation_counter # global variable
    global last_operation_counter_clear

    operation_counter_periodic_clear()

    if operation_counter >= MAX_OPERATIONS:
        return 'OVERLOAD'

    operation_counter += 1

    # Grab data and tweak the MOF name.
    my_data = json.loads(flask.request.get_data())
    structure, name = extract_info(my_data)

    temp_file_folder = MOFSIMPLIFY_PATH + "temp_file_creation_" + str(session['ID']) + '/'

    ### Check in MongoDB history.MOFSimplify_water_v2 collection to see if this structure has been predicted on before. 
    # Comment out this section if you are running MOFSimplify on your computer, and define is_entry as False. Also select "No" for the question May MOFSimplify store information on your MOFs?
    client = MongoClient('18.18.63.68',27017) # connect to mongodb. The first argument is the IP address. The second argument is the port.
    db = client.history # The history database
    collection = db.MOFSimplify_water_v2 # The MOFSimplify_water_v2 collection in the history database.
    my_documents = collection.find({'structure':structure})
    is_entry = my_documents.count() # will be zero or one, depending if structure is in the database
    if is_entry: # return statements in this if statement are of same form as those from the workflow if the structure isn't in the database. The frontend needs to receive information of the same form.
        entry_data = my_documents[0]
        if entry_data['failure'] == True: # MOF was not featurizable
            if session['permission']:
                db_push_lite_water(structure=structure, prediction_type='acid_stability_prediction') # add 1 to the a_times in the database
            operation_counter = conditional_diminish(operation_counter)
            return 'FAILED'
        else: 
            # Making the csv; can skip csv making in descriptor_generator with this
            csv_content = entry_data['csv_content']
            with open(temp_file_folder + 'merged_descriptors/' + name + '_descriptors_water.csv', 'w') as f:
                f.write(csv_content)

        if entry_data['a_times'] > 0:
            my_dict = {'prediction':entry_data['a_result']}
            if session['permission']:
                db_push_lite_water(structure=structure, prediction_type='acid_stability_prediction') # add 1 to the a_times in the database
            operation_counter = conditional_diminish(operation_counter)
            return my_dict

        # if haven't returned anything by now and entered the if is_entry statement, have a featurizable MOF for which an acid prediction has not been run, but a water prediction has

    output = descriptor_generator_water(name, structure, 'acid', is_entry) # generate descriptors
    if output == 'FAILED': # Description generation failure
        operation_counter = conditional_diminish(operation_counter)
        if session['permission']:
            db_push_water(structure=structure, prediction_type='acid_stability_prediction', result='', failure=True, csv_content='')
        return 'FAILED'

    # Applying the model next
    prediction = run_acid_RF(str(session['ID']), MOFSIMPLIFY_PATH, name)

    results = {
    'prediction': prediction,
    }
    
    if session['permission']:
        # Database push of MOF structure and ML results.
        # Getting the contents of the CSV containing features.
        descriptors_folder = temp_file_folder + "merged_descriptors/"
        with open(descriptors_folder + name + '_descriptors_water.csv', 'r') as f:
            csv_contents = f.read()
        db_push_water(structure=structure, prediction_type='acid_stability_prediction', result=prediction, failure=False, csv_content=csv_contents)

    operation_counter = conditional_diminish(operation_counter)
    return results

@app.route('/get_descriptors_water', methods=['POST']) 
def descriptor_getter_water():
    """
    descriptor_getter_water returns the contents of the csv with the descriptors of the desired MOF.
    These descriptors are RACs and Zeo++ descriptors. There are 189 total.
    The Zeo++ features are slightly different than those from descriptor_getter.

    :return: str contents, the contents of the csv with the descriptor data. To be written to a csv in the front end for download.
    """ 

    # Grab data
    name = json.loads(flask.request.get_data()); # This is the selected MOF
    if name[-4:] == '.cif':
        name = name[:-4] # remove the .cif part of the name

    temp_file_folder = MOFSIMPLIFY_PATH + "temp_file_creation_" + str(session['ID']) + '/'
    descriptors_folder = temp_file_folder + "merged_descriptors/"

    with open(descriptors_folder + name + '_descriptors_water.csv', 'r') as f:
        contents = f.read()

    return contents

def descriptor_generator_water(name, structure, prediction_type, is_entry):
    """
    # descriptor_generator_water is used by both ws_predict() and as_predict() to generate RACs and Zeo++ descriptors.
    # These descriptors are subsequently used in ws_predict() and as_predict() for the RF models.
    # Inputs are the name of the MOF and the structure (cif file text) of the MOF for which descriptors are to be generated.
    # The third input indicates the type of prediction (water or acid).
    # Unlike in descriptor_generator, no check to see if the new MOF is in the training data.

    :param name: str, the name of the MOF being analyzed.
    :param structure: str, the text of the cif file of the MOF being analyzed.
    :param prediction_type: str, the code for the type of description (water or acid)
    :param is_entry: boolean, indicates whether the descriptor CSV has already been written.
    :return: str, either the string 'FAILED' if descriptor generation fails, otherwise 'SUCCESS'
    """ 

    temp_file_folder = MOFSIMPLIFY_PATH + "temp_file_creation_" + str(session['ID']) + '/'
    cif_folder = temp_file_folder + 'cifs/'

    # Write the data back to a cif file.
    try:
        cif_file = open(cif_folder + name + '.cif', 'w')
    except FileNotFoundError:
        return 'FAILED'
    cif_file.write(structure)
    cif_file.close()

    # There can be a RACs folder for water predictions and a RACs folder for acid predictions. Same for Zeo++.
    RACs_folder = temp_file_folder + 'feature_generation/' + prediction_type + '_RACs/'
    zeo_folder = temp_file_folder + 'feature_generation/' + prediction_type + '_zeo++/'

    # Delete the RACs folder, then remake it (to start fresh for this prediction).
    shutil.rmtree(RACs_folder)
    os.mkdir(RACs_folder)

    # Doing the same with the Zeo++ folder.
    shutil.rmtree(zeo_folder)
    os.mkdir(zeo_folder)

    # If is_entry is True, do not need to generate features, since we already have them.
    if is_entry:
        return 'SUCCESS'

    # Next, running MOF featurization
    try:
        get_primitive(cif_folder + name + '.cif', cif_folder + name + '_primitive.cif');
    except ValueError:
        # return 'FAILED'

        # Will proceed even if primitive unit cell not obtained by pymatgen
        shutil.copy(cif_folder + name + '.cif', cif_folder + name + '_primitive.cif')

    # get_MOF_descriptors is used in RAC_getter.py to get RAC features.
        # The files that are generated from RAC_getter.py: lc_descriptors.csv, sbu_descriptors.csv, linker_descriptors.csv

    # cmd1, cmd2, and cmd3 are for Zeo++. cm4 is for RACs.
    cmd1 = MOFSIMPLIFY_PATH + 'zeo++-0.3/network -ha -res ' + zeo_folder + name + '_pd.txt ' + cif_folder + name + '_primitive.cif'
    cmd2 = MOFSIMPLIFY_PATH + 'zeo++-0.3/network -sa 1.4 1.4 10000 ' + zeo_folder + name + '_sa.txt ' + cif_folder + name + '_primitive.cif'
    cmd3 = MOFSIMPLIFY_PATH + 'zeo++-0.3/network -volpo 1.4 1.4 10000 ' + zeo_folder + name + '_pov.txt '+ cif_folder + name + '_primitive.cif'
    cmd4 = 'python ' + MOFSIMPLIFY_PATH + 'model/RAC_getter.py %s %s %s' %(cif_folder, name, RACs_folder)

    # four parallelized Zeo++ and RAC commands
    process1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=None, shell=True)
    process2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, stderr=None, shell=True)
    process3 = subprocess.Popen(cmd3, stdout=subprocess.PIPE, stderr=None, shell=True)
    process4 = subprocess.Popen(cmd4, stdout=subprocess.PIPE, stderr=None, shell=True)

    output1 = process1.communicate()[0]
    output2 = process2.communicate()[0]
    output3 = process3.communicate()[0]
    output4 = process4.communicate()[0]

    # Have written output of Zeo++ commands to files. Now, code below extracts information from those files.

    ''' 
    The geometric descriptors are:
    - the maximum included sphere (Di)
    - maximum free sphere (Df)
    - maximum included sphere in the free sphere path (Dif)
    - gravimetric pore volume (GPOV)
    - volumetric pore volume (VPOV)
    - gravimetric surface area (GSA)
    - volumetric surface area (VSA)
    - cell volume (cell_v)
    - gravimetric pore accessible volume (GPOAV)
    - gravimetric pore non-accessible volume (GPONAV)
    - pore-accessible volume (POAV)
    - pore non-accessible volume (PONAV)
    - pore accessible void fraction (POAVF)
    - pore nonaccessible void fraction (PONAVF)

    All Zeo++ calculations use a probe radius of 1.4 angstrom, and zeo++ is called by subprocess.
    '''

    dict_list = []
    cif_file = name + '_primitive.cif' 
    basename = cif_file.strip('.cif')
    largest_included_sphere, largest_free_sphere, largest_included_sphere_along_free_sphere_path = np.nan, np.nan, np.nan
    unit_cell_volume, crystal_density, VSA, GSA = np.nan, np.nan, np.nan, np.nan
    VPOV, GPOV = np.nan, np.nan
    POAV, PONAV, GPOAV, GPONAV, POAV_volume_fraction, PONAV_volume_fraction = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    if (os.path.exists(zeo_folder + name + '_pd.txt') & os.path.exists(zeo_folder + name + '_sa.txt') &
        os.path.exists(zeo_folder + name + '_pov.txt')):
        with open(zeo_folder + name + '_pd.txt') as f:
            pore_diameter_data = f.readlines()
            for row in pore_diameter_data:
                largest_included_sphere = float(row.split()[1]) # largest included sphere
                largest_free_sphere = float(row.split()[2]) # largest free sphere
                largest_included_sphere_along_free_sphere_path = float(row.split()[3]) # largest included sphere along free sphere path
        with open(zeo_folder + name + '_sa.txt') as f:
            surface_area_data = f.readlines()
            for i, row in enumerate(surface_area_data):
                if i == 0:
                    unit_cell_volume = float(row.split('Unitcell_volume:')[1].split()[0]) # unit cell volume
                    crystal_density = float(row.split('Density:')[1].split()[0]) # crystal density
                    VSA = float(row.split('ASA_m^2/cm^3:')[1].split()[0]) # volumetric surface area
                    GSA = float(row.split('ASA_m^2/g:')[1].split()[0]) # gravimetric surface area
        with open(zeo_folder + name + '_pov.txt') as f:
            pore_volume_data = f.readlines()
            for i, row in enumerate(pore_volume_data):
                if i == 0:
                    density = float(row.split('Density:')[1].split()[0])
                    POAV = float(row.split('POAV_A^3:')[1].split()[0]) # Probe accessible pore volume
                    PONAV = float(row.split('PONAV_A^3:')[1].split()[0]) # Probe non-accessible probe volume
                    GPOAV = float(row.split('POAV_cm^3/g:')[1].split()[0])
                    GPONAV = float(row.split('PONAV_cm^3/g:')[1].split()[0])
                    POAV_volume_fraction = float(row.split('POAV_Volume_fraction:')[1].split()[0]) # probe accessible volume fraction
                    PONAV_volume_fraction = float(row.split('PONAV_Volume_fraction:')[1].split()[0]) # probe non accessible volume fraction
                    VPOV = POAV_volume_fraction+PONAV_volume_fraction
                    GPOV = VPOV/density
    else:
        print('Not all 3 files exist, so at least one Zeo++ call failed!', 'sa: ',os.path.exists(zeo_folder + name + '_sa.txt'), 
              '; pd: ',os.path.exists(zeo_folder + name + '_pd.txt'), '; pov: ', os.path.exists(zeo_folder + name + '_pov.txt'))
        return 'FAILED'
    geo_dict = {'name': basename, 'cif_file': cif_file, 'Di': largest_included_sphere, 'Df': largest_free_sphere, 'Dif': largest_included_sphere_along_free_sphere_path,
                'cell_v': unit_cell_volume, 'VSA': VSA, 'GSA': GSA, 'VPOV': VPOV, 'GPOV': GPOV, 'POAV_vol_frac': POAV_volume_fraction, 
                'PONAV_vol_frac': PONAV_volume_fraction, 'GPOAV': GPOAV, 'GPONAV': GPONAV, 'POAV': POAV, 'PONAV': PONAV}
    dict_list.append(geo_dict)
    geo_df = pd.DataFrame(dict_list)
    geo_df.to_csv(zeo_folder + 'geometric_parameters.csv',index=False)

    # error handling for cmd4
    with open(RACs_folder + 'RAC_getter_log.txt', 'r') as f:
        if f.readline() == 'FAILED':
            print('RAC generation failed.')
            return 'FAILED'

    # Merging geometric information with the RAC information that is in the get_MOF_descriptors-generated files (lc_descriptors.csv, sbu_descriptors.csv, linker_descriptors.csv)
    try:
        lc_df = pd.read_csv(RACs_folder + "lc_descriptors.csv") 
        sbu_df = pd.read_csv(RACs_folder + "sbu_descriptors.csv")
        linker_df = pd.read_csv(RACs_folder + "linker_descriptors.csv")
    except Exception: # csv files have been deleted
        return 'FAILED' 

    lc_df = lc_df.mean().to_frame().transpose() # averaging over all rows. Convert resulting Series into a DataFrame, then transpose
    sbu_df = sbu_df.mean().to_frame().transpose()
    linker_df = linker_df.mean().to_frame().transpose()

    merged_df = pd.concat([geo_df, lc_df, sbu_df, linker_df], axis=1)

    merged_df.to_csv(temp_file_folder + '/merged_descriptors/' + name + '_descriptors_water.csv',index=False) # written in /temp_file_creation_SESSIONID

    return 'SUCCESS' 

def run_water_RF(user_id, path, MOF_name):
    """
    run_water_RF runs the water stability RF with the desired MOF as input.
    It returns a prediction between zero and one. This prediction corresponds to an assessment of MOF water stability.
    The further the prediction is from 0.5, the more sure the ANN is.

    :param user_id: str, the session ID of the user
    :param path: str, the server's path to the MOFSimplify folder on the server
    :param MOF_name: str, the name of the MOF for which a prediction is being generated
    :return: str, the model water stability prediction 
    """ 

    rfa_2_class_water_features = [
    'mc-Z-3-all', 
    'D_mc-Z-3-all', 
    'D_mc-Z-2-all', 
    'D_mc-Z-1-all', 
    'mc-chi-3-all', 
    'mc-Z-1-all', 
    'mc-Z-0-all', 
    'D_mc-chi-2-all', 
    'f-lig-Z-2', 
    'GSA', 
    'f-lig-I-0', 
    'func-S-1-all',
    ]
     
    temp_file_path = path + 'temp_file_creation_' + user_id + '/'
    df_newMOF = pd.read_csv(temp_file_path + 'merged_descriptors/' + MOF_name + '_descriptors_water.csv') # assumes that temp_file_creation/ is in parent folder
    X_newMOF = df_newMOF[rfa_2_class_water_features].to_numpy()
    X_newMOF = water_scaler.transform(X_newMOF)
    water_pred = water_model.predict_proba(X_newMOF)[:,1][0]
    water_pred = np.round(water_pred, 2) # round to 2 decimals

    return str(water_pred)

def run_acid_RF(user_id, path, MOF_name):
    """
    run_acid_RF runs the acid stability RF with the desired MOF as input.
    It returns a prediction between zero and one. This prediction corresponds to an assessment of MOF acid stability.
    The further the prediction is from 0.5, the more sure the ANN is.

    :param user_id: str, the session ID of the user
    :param path: str, the server's path to the MOFSimplify folder on the server
    :param MOF_name: str, the name of the MOF for which a prediction is being generated
    :return: str, the model acid stability prediction 
    """ 

    rfa_acid_features = [
    'mc-chi-3-all', 
    'Dif', 
    'mc-Z-2-all', 
    'Di', 
    'f-T-2-all', 
    'D_mc-chi-2-all', 
    'lc-chi-3-all', 
    'D_mc-S-1-all',
    ]
     
    temp_file_path = path + 'temp_file_creation_' + user_id + '/'
    df_newMOF = pd.read_csv(temp_file_path + 'merged_descriptors/' + MOF_name + '_descriptors_water.csv') # assumes that temp_file_creation/ is in parent folder
    X_newMOF = df_newMOF[rfa_acid_features].to_numpy()
    X_newMOF = acid_scaler.transform(X_newMOF)
    acid_pred = acid_model.predict_proba(X_newMOF)[:,1][0]
    acid_pred = np.round(acid_pred, 2) # round to 2 decimals

    return str(acid_pred)

### New section. The code for C2 uptake prediction. ### 

def db_push_C2(structure, temperature, pressure, prediction_type, result, failure, csv_content):
    """
    db_push_C2 sends the structure of the MOF being predicted on, and other information, to the MOFSimplify_C2_v2 collection in the MongoDB history database.

    :param structure: str, the text of the cif file of the MOF being analyzed.
    :param temperature: str, the temperature of the data point.
    :param pressure: str, the pressure of the data point.
    :param prediction_type: str, the type of prediction being run. Can either be 'ethane_uptake' or 'ethylene_uptake'.
    :param result: float, the prediction for this structure.
    :param failure: boolean, whether the featurization failed.
    :param csv_content: string, the content of the csv of features.
    :return: A 204 no content response, so the front end does not display a page different than the one it is on.
    """
    client = MongoClient('18.18.63.68',27017) # connect to mongodb
    # The first argument is the IP address. The second argument is the port.
    db = client.history # The history database
    collection = db.MOFSimplify_C2_v2 # The MOFSimplify_C2_v2 collection in the history database.
    # ethane_result refers to an ethane uptake prediction; ethylene_result refers to an ethylene uptake prediction
    fields = ['structure', 'temperature', 'pressure', 'ethane_result', 'ethylene_result', 'failure', 'csv_content', 'ethane_times', 'ethylene_times'] 
        # ethane_times and ethylene_times indicate the number of times this structure has had its ethane or ethylene uptake predicted, at this temperature and pressure
    final_dict = {}
    for field in fields:
        final_dict[field] = '' # start with everything empty

    # https://docs.mongodb.com/manual/reference/method/db.collection.update/
    my_documents = collection.find({'structure':structure, 'temperature':temperature, 'pressure':pressure})
    is_entry = my_documents.count() # will be zero or one, depending if structure is in the database

    final_dict['structure'] = structure
    final_dict['temperature'] = temperature
    final_dict['pressure'] = pressure
    final_dict['failure'] = failure
    final_dict['csv_content'] = csv_content

    if is_entry == 0: # structure has not been seen before. Make a new entry
        # Populate certain fields
        if prediction_type == 'ethane_uptake':
            final_dict['ethane_result'] = result
            final_dict['ethane_times'] = 1
            final_dict['ethylene_times'] = 0
        elif prediction_type == 'ethylene_uptake':
            final_dict['ethylene_result'] = result
            final_dict['ethane_times'] = 0
            final_dict['ethylene_times'] = 1

        collection.insert(final_dict) # insert the dictionary into the mongodb collection
    else:  # existing structure. Update existing document (couldn't get actual update function to work, so I delete entry and write a new one.)
        my_document = my_documents[0]
        final_dict['ethane_result'] = my_document['ethane_result']
        final_dict['ethylene_result'] = my_document['ethylene_result']
        final_dict['ethane_times'] = my_document['ethane_times']
        final_dict['ethylene_times'] = my_document['ethylene_times']

        if prediction_type == 'ethane_uptake': # structure has had an ethylene prediction but not an ethane prediction
            final_dict['ethane_result'] = result
            final_dict['ethane_times'] = final_dict['ethane_times'] + 1
        elif prediction_type == 'ethylene_uptake': # structure has had an ethane prediction but not an ethylene prediction
            final_dict['ethylene_result'] = result
            final_dict['ethylene_times'] = final_dict['ethylene_times'] + 1

        collection.remove({'structure':structure, 'temperature':temperature, 'pressure':pressure}) # delete entry
        collection.insert(final_dict) # insert the dictionary into the mongodb collection (so, write new entry)
    return ('', 204) # 204 no content response

def db_push_lite_C2(structure, temperature, pressure, prediction_type):
    """
    # db_push_lite_C2 increases either ethane_times or ethylene_times in the document corresponding to structure, in the MOFSimplify_C2_v2 collection in the MongoDB history database.
    # The prediction_type determines whether ethane_times or ethylene_times is increased by one.

    :param structure: str, the text of the cif file of the MOF being analyzed.
    :param temperature: str, the temperature of the data point.
    :param pressure: str, the pressure of the data point.
    :param prediction_type: str, the type of prediction being run. Can either be 'ethane_uptake' or 'ethylene_uptake'.
    :return: A 204 no content response, so the front end does not display a page different than the one it is on.
    """ 

    client = MongoClient('18.18.63.68',27017) # connect to mongodb
    # The first argument is the IP address. The second argument is the port.
    db = client.history # The history database
    collection = db.MOFSimplify_C2_v2 # The MOFSimplify_C2_v2 collection in the history database.
    final_dict = {}

    # https://docs.mongodb.com/manual/reference/method/db.collection.update/
    my_documents = collection.find({'structure':structure, 'temperature':temperature, 'pressure':pressure})

    # if db_push_lite_C2 is called, structure passed to it will already have a document in the history.MOFSimplify_C2_v2 collection
    # for the given temperature and pressure

    # Update existing document (couldn't get actual update function to work, so I delete entry and write a new one.)
    my_document = my_documents[0]
    final_dict['structure'] = structure
    final_dict['temperature'] = temperature
    final_dict['pressure'] = pressure
    final_dict['ethane_result'] = my_document['ethane_result']
    final_dict['ethylene_result'] = my_document['ethylene_result']
    final_dict['ethane_times'] = my_document['ethane_times']
    final_dict['ethylene_times'] = my_document['ethylene_times']
    final_dict['failure'] = my_document['failure']
    final_dict['csv_content'] = my_document['csv_content']

    if prediction_type == 'ethane_uptake': 
        final_dict['ethane_times'] += 1
    elif prediction_type == 'ethylene_uptake': 
        final_dict['ethylene_times'] += 1

    collection.remove({'structure':structure, 'temperature':temperature, 'pressure':pressure}) # delete entry
    collection.insert(final_dict) # insert the dictionary into the mongodb collection (so, write new entry)
    return ('', 204) # 204 no content response

@app.route('/get_descriptors_C2', methods=['POST']) 
def descriptor_getter_C2():
    """
    descriptor_getter_C2 returns the contents of the csv with the descriptors of the desired MOF.
    These descriptors are RACs and Zeo++ descriptors, plus temperature and pressure.

    :return: str contents, the contents of the csv with the descriptor data. To be written to a csv in the front end for download.
    """ 

    # Grab data
    name = json.loads(flask.request.get_data()); # This is the selected MOF
    if name[-4:] == '.cif':
        name = name[:-4] # remove the .cif part of the name

    temp_file_folder = MOFSIMPLIFY_PATH + "temp_file_creation_" + str(session['ID']) + '/'
    descriptors_folder = temp_file_folder + "merged_descriptors/"

    with open(descriptors_folder + name + '_descriptors_C2.csv', 'r') as f:
        contents = f.read()

    return contents

@app.route('/predict_ethane_uptake', methods=['POST']) 
def ethane_predict():
    """
    ethane_predict generates the ethane uptake prediction for the selected MOF.
        Or it will return a signal that descriptor generation failed and thus a prediction cannot be made. 
    RAC featurization and Zeo++ geometry information for the selected MOF is generated, using descriptor_generator_C2.
    Then, Matt's model is applied to make a prediction using run_C2_model.

    :return: dict results, contains the prediction.
        May instead return a string 'FAILED' if descriptor generation fails.
        May instead return a string 'OVERLOAD' if the server is being asked too much of, in order to avoid 50X errors (like 504, aka Gateway Time-out).
    """ 

    global operation_counter # global variable
    global last_operation_counter_clear

    operation_counter_periodic_clear()

    if operation_counter >= MAX_OPERATIONS:
        return 'OVERLOAD'

    operation_counter += 1

    # Grab data and tweak the MOF name.
    my_data = json.loads(flask.request.get_data())
    structure, name, temperature, pressure = extract_info_C2(my_data)

    temp_file_folder = MOFSIMPLIFY_PATH + "temp_file_creation_" + str(session['ID']) + '/'

    ### Check in MongoDB history.MOFSimplify_C2_v2 collection to see if this structure has been predicted on before. 
    # Comment out this section if you are running MOFSimplify on your computer, and define is_entry as False. Also select "No" for the question May MOFSimplify store information on your MOFs?
    client = MongoClient('18.18.63.68',27017) # connect to mongodb. The first argument is the IP address. The second argument is the port.
    db = client.history # The history database
    collection = db.MOFSimplify_C2_v2 # The MOFSimplify_C2_v2 collection in the history database.
    my_documents = collection.find({'structure':structure, 'temperature':temperature, 'pressure':pressure})
    is_entry = my_documents.count() # will be zero or one, depending if structure is in the database
    if is_entry: # return statements in this if statement are of same form as those from the workflow if the structure isn't in the database. The frontend needs to receive information of the same form.
        entry_data = my_documents[0]
        if entry_data['failure'] == True: # MOF was not featurizable
            if session['permission']:
                db_push_lite_C2(structure=structure, temperature=temperature, pressure=pressure, prediction_type='ethane_uptake') # add 1 to the ethane_times in the database
            operation_counter = conditional_diminish(operation_counter)
            return 'FAILED'
        else: 
            # Making the csv; can skip csv making in descriptor_generator with this
            csv_content = entry_data['csv_content']
            with open(temp_file_folder + 'merged_descriptors/' + name + '_descriptors_C2.csv', 'w') as f:
                f.write(csv_content)

        if entry_data['ethane_times'] > 0:
            my_dict = {'prediction':entry_data['ethane_result']}
            if session['permission']:
                db_push_lite_C2(structure=structure, temperature=temperature, pressure=pressure, prediction_type='ethane_uptake') # add 1 to the ethane_times in the database
            operation_counter = conditional_diminish(operation_counter)
            return my_dict

        # if haven't returned anything by now and entered the if is_entry statement, have a featurizable MOF for which an ethylene prediction has not been run, but an ethane prediction has

    output = descriptor_generator_C2(name, structure, temperature, pressure, 'ethane') # generate descriptors
    if output == 'FAILED': # Description generation failure
        operation_counter = conditional_diminish(operation_counter)
        if session['permission']:
            db_push_C2(structure=structure, temperature=temperature, pressure=pressure, prediction_type='ethane_uptake', result='', failure=True, csv_content='')
        return 'FAILED'

    # Applying the model next
    prediction = run_C2_model(str(session['ID']), MOFSIMPLIFY_PATH, name, 'ethane')

    results = {
    'prediction': prediction,
    }

    if session['permission']:
        # Database push of MOF structure and ML results.
        # Getting the contents of the CSV containing features.
        descriptors_folder = temp_file_folder + "merged_descriptors/"
        with open(descriptors_folder + name + '_descriptors_C2.csv', 'r') as f:
            csv_contents = f.read()
        db_push_C2(structure=structure, temperature=temperature, pressure=pressure, prediction_type='ethane_uptake', result=prediction, failure=False, csv_content=csv_contents)    
    
    operation_counter = conditional_diminish(operation_counter)
    return results

@app.route('/predict_ethylene_uptake', methods=['POST']) 
def ethylene_predict():
    """
    ethylene_predict generates the ethylene uptake prediction for the selected MOF.
        Or it will return a signal that descriptor generation failed and thus a prediction cannot be made. 
    RAC featurization and Zeo++ geometry information for the selected MOF is generated, using descriptor_generator_C2.
    Then, Matt's model is applied to make a prediction using run_C2_model.

    :return: dict results, contains the prediction.
        May instead return a string 'FAILED' if descriptor generation fails.
        May instead return a string 'OVERLOAD' if the server is being asked too much of, in order to avoid 50X errors (like 504, aka Gateway Time-out).
    """ 

    global operation_counter # global variable
    global last_operation_counter_clear

    operation_counter_periodic_clear()

    if operation_counter >= MAX_OPERATIONS:
        return 'OVERLOAD'

    operation_counter += 1

    # Grab data and tweak the MOF name.
    my_data = json.loads(flask.request.get_data())
    structure, name, temperature, pressure = extract_info_C2(my_data)

    temp_file_folder = MOFSIMPLIFY_PATH + "temp_file_creation_" + str(session['ID']) + '/'

    ### Check in MongoDB history.MOFSimplify_C2_v2 collection to see if this structure has been predicted on before. 
    # Comment out this section if you are running MOFSimplify on your computer, and define is_entry as False. Also select "No" for the question May MOFSimplify store information on your MOFs?
    client = MongoClient('18.18.63.68',27017) # connect to mongodb. The first argument is the IP address. The second argument is the port.
    db = client.history # The history database
    collection = db.MOFSimplify_C2_v2 # The MOFSimplify_C2_v2 collection in the history database.
    my_documents = collection.find({'structure':structure, 'temperature':temperature, 'pressure':pressure})
    is_entry = my_documents.count() # will be zero or one, depending if structure is in the database
    if is_entry: # return statements in this if statement are of same form as those from the workflow if the structure isn't in the database. The frontend needs to receive information of the same form.
        entry_data = my_documents[0]
        if entry_data['failure'] == True: # MOF was not featurizable
            if session['permission']:
                db_push_lite_C2(structure=structure, temperature=temperature, pressure=pressure, prediction_type='ethylene_uptake') # add 1 to the ethylene_times in the database
            operation_counter = conditional_diminish(operation_counter)
            return 'FAILED'
        else: 
            # Making the csv; can skip csv making in descriptor_generator with this
            csv_content = entry_data['csv_content']
            with open(temp_file_folder + 'merged_descriptors/' + name + '_descriptors_C2.csv', 'w') as f:
                f.write(csv_content)

        if entry_data['ethylene_times'] > 0:
            my_dict = {'prediction':entry_data['ethylene_result']}
            if session['permission']:
                db_push_lite_C2(structure=structure, temperature=temperature, pressure=pressure, prediction_type='ethylene_uptake') # add 1 to the ethylene_times in the database
            operation_counter = conditional_diminish(operation_counter)
            return my_dict

        # if haven't returned anything by now and entered the if is_entry statement, have a featurizable MOF for which an ethane prediction has not been run, but an ethylene prediction has

    output = descriptor_generator_C2(name, structure, temperature, pressure, 'ethylene') # generate descriptors
    if output == 'FAILED': # Description generation failure
        operation_counter = conditional_diminish(operation_counter)
        if session['permission']:
            db_push_C2(structure=structure, temperature=temperature, pressure=pressure, prediction_type='ethylene_uptake', result='', failure=True, csv_content='')
        return 'FAILED'

    # Applying the model next
    prediction = run_C2_model(str(session['ID']), MOFSIMPLIFY_PATH, name, 'ethylene')

    results = {
    'prediction': prediction,
    }
    
    if session['permission']:
        # Database push of MOF structure and ML results.
        # Getting the contents of the CSV containing features.
        descriptors_folder = temp_file_folder + "merged_descriptors/"
        with open(descriptors_folder + name + '_descriptors_C2.csv', 'r') as f:
            csv_contents = f.read()
        db_push_C2(structure=structure, temperature=temperature, pressure=pressure, prediction_type='ethylene_uptake', result=prediction, failure=False, csv_content=csv_contents)  

    operation_counter = conditional_diminish(operation_counter)
    return results

def channel_surface_area_extract(row):
    # The row will look like 
    # Number_of_channels: 2 Channel_surface_area_A^2: 712.983  709.065
    # or
    # Number_of_channels: 1 Channel_surface_area_A^2: 2038.46 
    # etc

    split_row = row.split()
    num_channels = int(split_row[1])

    if num_channels == 0:
        CSA_1, CSA_2 = 0, 0
    elif num_channels == 1:
        CSA_1, CSA_2 = float(split_row[3]), 0
    else: # There are 2+ channels
        CSA_1, CSA_2 = float(split_row[3]), float(split_row[4])

    return CSA_1, CSA_2

def descriptor_generator_C2(name, structure, temperature, pressure, prediction_type):
    """
    # descriptor_generator_C2 is used by both ethane_predict() and ethylene_predict() to generate RACs and Zeo++ descriptors.
    # These descriptors are subsequently used in ethane_predict() and ethylene_predict() for the ANN models.
    # First two inputs are the name of the MOF and the structure (cif file text) of the MOF for which descriptors are to be generated.
    # Next two inputs are the operating conditions.
    # The fifth input indicates the type of prediction (ethane or ethylene).
    # Unlike in descriptor_generator, no check to see if the new MOF is in the training data.

    :param name: str, the name of the MOF being analyzed.
    :param structure: str, the text of the cif file of the MOF being analyzed.
    :param temperature: str, temperature in Kelvin
    :param pressure: str, pressure in kPa
    :param prediction_type: str, indicating whether descriptors are made for ethane or ethylene
    :return: str, either the string 'FAILED' if descriptor generation fails, otherwise 'SUCCESS'
    """ 

    temp_file_folder = MOFSIMPLIFY_PATH + "temp_file_creation_" + str(session['ID']) + '/'
    cif_folder = temp_file_folder + 'cifs/'

    # Write the data back to a cif file.
    try:
        cif_file = open(cif_folder + name + '.cif', 'w')
    except FileNotFoundError:
        return 'FAILED'
    cif_file.write(structure)
    cif_file.close()

    # There can be a RACs folder for ethane predictions and a RACs folder for ethylene predictions. Same for Zeo++.
    RACs_folder = temp_file_folder + 'feature_generation/' + prediction_type + '_RACs/'
    zeo_folder = temp_file_folder + 'feature_generation/' + prediction_type + '_zeo++/'

    # Delete the RACs folder, then remake it (to start fresh for this prediction).
    shutil.rmtree(RACs_folder)
    os.mkdir(RACs_folder)

    # Doing the same with the Zeo++ folder.
    shutil.rmtree(zeo_folder)
    os.mkdir(zeo_folder)

    # # If is_entry is True, do not need to generate features, since we already have them.
    # if is_entry:
    #     return 'SUCCESS'

    # Next, running MOF featurization
    try:
        get_primitive(cif_folder + name + '.cif', cif_folder + name + '_primitive.cif');
    except ValueError:
        # return 'FAILED'

        # Will proceed even if primitive unit cell not obtained by pymatgen
        shutil.copy(cif_folder + name + '.cif', cif_folder + name + '_primitive.cif')

    # get_MOF_descriptors is used in RAC_getter.py to get RAC features.
        # The files that are generated from RAC_getter.py: lc_descriptors.csv, sbu_descriptors.csv, linker_descriptors.csv

    # cmd1, cmd2, and cmd3 are for Zeo++. cm4 is for RACs.
    cmd1 = MOFSIMPLIFY_PATH + 'zeo++-0.3/network -ha -res ' + zeo_folder + name + '_pd.txt ' + cif_folder + name + '_primitive.cif'
    cmd2 = MOFSIMPLIFY_PATH + 'zeo++-0.3/network -sa 1.86 1.86 10000 ' + zeo_folder + name + '_sa.txt ' + cif_folder + name + '_primitive.cif'
    cmd3 = MOFSIMPLIFY_PATH + 'zeo++-0.3/network -volpo 1.86 1.86 10000 ' + zeo_folder + name + '_pov.txt '+ cif_folder + name + '_primitive.cif'
    cmd4 = 'python ' + MOFSIMPLIFY_PATH + 'model/RAC_getter.py %s %s %s' %(cif_folder, name, RACs_folder)

    # four parallelized Zeo++ and RAC commands
    process1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=None, shell=True)
    process2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, stderr=None, shell=True)
    process3 = subprocess.Popen(cmd3, stdout=subprocess.PIPE, stderr=None, shell=True)
    process4 = subprocess.Popen(cmd4, stdout=subprocess.PIPE, stderr=None, shell=True)

    output1 = process1.communicate()[0]
    output2 = process2.communicate()[0]
    output3 = process3.communicate()[0]
    output4 = process4.communicate()[0]

    # Have written output of Zeo++ commands to files. Now, code below extracts information from those files.

    dict_list = []
    cif_file = name + '_primitive.cif' 
    basename = cif_file.strip('.cif')
    largest_included_sphere, largest_free_sphere, largest_included_sphere_along_free_sphere_path = np.nan, np.nan, np.nan
    unit_cell_volume, crystal_density, VSA, GSA = np.nan, np.nan, np.nan, np.nan
    VPOV, GPOV = np.nan, np.nan
    POAV, PONAV, GPOAV, GPONAV, POAV_volume_fraction, PONAV_volume_fraction = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    if (os.path.exists(zeo_folder + name + '_pd.txt') & os.path.exists(zeo_folder + name + '_sa.txt') &
        os.path.exists(zeo_folder + name + '_pov.txt')):
        with open(zeo_folder + name + '_pd.txt') as f:
            pore_diameter_data = f.readlines()
            for row in pore_diameter_data:
                largest_included_sphere = float(row.split()[1]) # largest included sphere
                largest_free_sphere = float(row.split()[2]) # largest free sphere
                largest_included_sphere_along_free_sphere_path = float(row.split()[3]) # largest included sphere along free sphere path
        with open(zeo_folder + name + '_sa.txt') as f:
            surface_area_data = f.readlines()
            for i, row in enumerate(surface_area_data):
                if i == 0:
                    unit_cell_volume = float(row.split('Unitcell_volume:')[1].split()[0]) # unit cell volume
                    crystal_density = float(row.split('Density:')[1].split()[0]) # crystal density
                    VSA = float(row.split('ASA_m^2/cm^3:')[1].split()[0]) # volumetric surface area
                    GSA = float(row.split('ASA_m^2/g:')[1].split()[0]) # gravimetric surface area
                    ASA = float(row.split('ASA_A^2:')[1].split()[0]) # accessible surface area
                elif i == 1:
                    # Extracting channel surface areas
                    CSA_1, CSA_2 = channel_surface_area_extract(row)
        with open(zeo_folder + name + '_pov.txt') as f:
            pore_volume_data = f.readlines()
            for i, row in enumerate(pore_volume_data):
                if i == 0:
                    density = float(row.split('Density:')[1].split()[0])
                    POAV = float(row.split('POAV_A^3:')[1].split()[0]) # Probe accessible pore volume
                    PONAV = float(row.split('PONAV_A^3:')[1].split()[0]) # Probe non-accessible probe volume
                    GPOAV = float(row.split('POAV_cm^3/g:')[1].split()[0])
                    GPONAV = float(row.split('PONAV_cm^3/g:')[1].split()[0])
                    POAV_volume_fraction = float(row.split('POAV_Volume_fraction:')[1].split()[0]) # probe accessible volume fraction
                    PONAV_volume_fraction = float(row.split('PONAV_Volume_fraction:')[1].split()[0]) # probe non accessible volume fraction
                    VPOV = POAV_volume_fraction+PONAV_volume_fraction
                    GPOV = VPOV/density
    else:
        print('Not all 3 files exist, so at least one Zeo++ call failed!', 'sa: ',os.path.exists(zeo_folder + name + '_sa.txt'), 
              '; pd: ',os.path.exists(zeo_folder + name + '_pd.txt'), '; pov: ', os.path.exists(zeo_folder + name + '_pov.txt'))
        return 'FAILED'
    geo_dict = {'name': basename, 'cif_file': cif_file, 'Di': largest_included_sphere, 'Df': largest_free_sphere, 'Dif': largest_included_sphere_along_free_sphere_path,
                'cell_v': unit_cell_volume, 'VSA': VSA, 'GSA': GSA, 'ASA': ASA, 'VPOV': VPOV, 'GPOV': GPOV, 'POAV_vol_frac': POAV_volume_fraction, 
                'PONAV_vol_frac': PONAV_volume_fraction, 'GPOAV': GPOAV, 'GPONAV': GPONAV, 'POAV': POAV, 'PONAV': PONAV, 'density': crystal_density,
                'CSA_1': CSA_1, 'CSA_2': CSA_2}
    dict_list.append(geo_dict)
    geo_df = pd.DataFrame(dict_list)
    geo_df.to_csv(zeo_folder + 'geometric_parameters.csv',index=False)

    # error handling for cmd4
    with open(RACs_folder + 'RAC_getter_log.txt', 'r') as f:
        if f.readline() == 'FAILED':
            print('RAC generation failed.')
            return 'FAILED'

    # Merging geometric information with the RAC information that is in the get_MOF_descriptors-generated files (lc_descriptors.csv, sbu_descriptors.csv, linker_descriptors.csv)
    try:
        lc_df = pd.read_csv(RACs_folder + "lc_descriptors.csv") 
        sbu_df = pd.read_csv(RACs_folder + "sbu_descriptors.csv")
        linker_df = pd.read_csv(RACs_folder + "linker_descriptors.csv")
    except Exception: # csv files have been deleted
        return 'FAILED' 

    lc_df = lc_df.mean().to_frame().transpose() # averaging over all rows. Convert resulting Series into a DataFrame, then transpose
    sbu_df = sbu_df.mean().to_frame().transpose()
    linker_df = linker_df.mean().to_frame().transpose()

    merged_df = pd.concat([geo_df, lc_df, sbu_df, linker_df], axis=1)

    # Adding temperature and pressure information
    merged_df['temperature'] = float(temperature)
    merged_df['pressure'] = float(pressure)

    merged_df.to_csv(temp_file_folder + '/merged_descriptors/' + name + '_descriptors_C2.csv',index=False) # written in /temp_file_creation_SESSIONID

    return 'SUCCESS' 

def normalize_data_C2(df_train, df_newMOF, fnames):
    """
    normalize_data_C2 takes in two DataFrames df_train and df_newMOF, one for the training data (many rows) and one for the new MOF (one row) for which a prediction is to be generated.
    This function also takes in fnames (the feature names).
    This function normalizes the X values from the pandas DataFrames and returns them as X_train and X_newMOF.
    The function also returns x_scaler (which scaled X_train).

    :param df_train: A pandas DataFrame of the training data.
    :param df_newMOF: A pandas DataFrame of the new MOF being analyzed.
    :param fnames: An array of column names of the descriptors.
    :return: numpy.ndarray X_train, the descriptors of the training data. Its number of rows is the number of MOFs in the training data. Its number of columns is the number of descriptors.
    :return: numpy.ndarray X_newMOF, the descriptors of the new MOF being analyzed by MOFSimplify. It contains only one row.
    :return: sklearn.preprocessing._data.StandardScaler x_scaler, the scaler used to normalize the descriptor data to unit mean and a variance of 1.
    """ 
    X_train, X_newMOF = df_train.values, df_newMOF[fnames].values # takes care of ensuring ordering is same for both X
    x_scaler = sklearn.preprocessing.StandardScaler()
    x_scaler.fit(X_train)
    X_train = x_scaler.transform(X_train)
    X_newMOF = x_scaler.transform(X_newMOF)
    return X_train, X_newMOF, x_scaler

def run_C2_model(user_id, path, MOF_name, gas='ethane'):
    """
    run_C2_ANN runs the ethane or ethylene uptake ANN with the desired MOF as input.

    :param user_id: str, the session ID of the user
    :param path: str, the server's path to the MOFSimplify folder on the server
    :param MOF_name: str, the name of the MOF for which a prediction is being generated
    :param gas: str, either "ethane" or "ethylene"
    :return: str str(new_MOF_pred[0][0]), the model ethane uptake prediction 
    """ 

    features = ['f-chi-0-all', 'f-chi-1-all', 'f-chi-2-all', 'f-chi-3-all', 'f-Z-0-all',
       'f-Z-1-all', 'f-Z-2-all', 'f-Z-3-all', 'f-I-0-all', 'f-I-1-all',
       'f-I-2-all', 'f-I-3-all', 'f-T-0-all', 'f-T-1-all', 'f-T-2-all',
       'f-T-3-all', 'f-S-0-all', 'f-S-1-all', 'f-S-2-all', 'f-S-3-all',
       'mc-chi-0-all', 'mc-chi-1-all', 'mc-chi-2-all', 'mc-chi-3-all',
       'mc-Z-0-all', 'mc-Z-1-all', 'mc-Z-2-all', 'mc-Z-3-all', 'mc-I-0-all',
       'mc-I-1-all', 'mc-I-2-all', 'mc-I-3-all', 'mc-T-0-all', 'mc-T-1-all',
       'mc-T-2-all', 'mc-T-3-all', 'mc-S-0-all', 'mc-S-1-all', 'mc-S-2-all',
       'mc-S-3-all', 'D_mc-chi-0-all', 'D_mc-chi-1-all', 'D_mc-chi-2-all',
       'D_mc-chi-3-all', 'D_mc-Z-0-all', 'D_mc-Z-1-all', 'D_mc-Z-2-all',
       'D_mc-Z-3-all', 'D_mc-I-0-all', 'D_mc-I-1-all', 'D_mc-I-2-all',
       'D_mc-I-3-all', 'D_mc-T-0-all', 'D_mc-T-1-all', 'D_mc-T-2-all',
       'D_mc-T-3-all', 'D_mc-S-0-all', 'D_mc-S-1-all', 'D_mc-S-2-all',
       'D_mc-S-3-all', 'f-lig-chi-0', 'f-lig-chi-1', 'f-lig-chi-2',
       'f-lig-chi-3', 'f-lig-Z-0', 'f-lig-Z-1', 'f-lig-Z-2', 'f-lig-Z-3',
       'f-lig-I-0', 'f-lig-I-1', 'f-lig-I-2', 'f-lig-I-3', 'f-lig-T-0',
       'f-lig-T-1', 'f-lig-T-2', 'f-lig-T-3', 'f-lig-S-0', 'f-lig-S-1',
       'f-lig-S-2', 'f-lig-S-3', 'ASA', 'GSA',
       'CSA_1', 'CSA_2', 'density',
       'Df', 'Di',
       'Dif', 'POAV',
       'POAV_vol_frac', 'GPOAV', 'cell_v', 'temperature',
       'pressure']
       # CSA: channel surface area
       # Only taking the first two channels

    if gas == 'ethane':   
        ANN_path = path + 'model/c2/ethane/'
        model = ethane_model
    elif gas == 'ethylene':
        ANN_path = path + 'model/c2/ethylene/'
        model = ethylene_model
    else:
        raise Exception('Invalid value of gas parameter.')
    temp_file_path = path + 'temp_file_creation_' + user_id + '/'
    df_train = pd.read_csv(ANN_path+'unnormalized_features/x_train.csv')
    df_newMOF = pd.read_csv(temp_file_path + 'merged_descriptors/' + MOF_name + '_descriptors_C2.csv')

    # The normalize_data_C2 function is expecting a DataFrame with each MOF in a separate row, and features in columns

    ### Utilize the function below to normalize the features of the new MOF
    # newMOF refers to the MOF that has been uploaded to MOFSimplify, for which a prediction will be generated
    X_train, X_newMOF, x_scaler = normalize_data_C2(df_train, df_newMOF, features)
    # Order of values in X_newMOF matters, but this is taken care of in normalize_data_C2.

    from tensorflow.python.keras.backend import set_session
    with tf_session.as_default(): # session stuff is needed because the model was loaded from h5 a while ago
        with tf_session.graph.as_default():
            new_MOF_pred = np.round(model.predict(X_newMOF),2) # round to 2 decimals

    return str(new_MOF_pred[0][0])


######

# ts_predict()

# Run with MOFSimplify conda env


RACs = ['D_func-I-0-all','D_func-I-1-all','D_func-I-2-all','D_func-I-3-all',
 'D_func-S-0-all', 'D_func-S-1-all', 'D_func-S-2-all', 'D_func-S-3-all',
 'D_func-T-0-all', 'D_func-T-1-all', 'D_func-T-2-all', 'D_func-T-3-all',
 'D_func-Z-0-all', 'D_func-Z-1-all', 'D_func-Z-2-all', 'D_func-Z-3-all',
 'D_func-chi-0-all', 'D_func-chi-1-all', 'D_func-chi-2-all',
 'D_func-chi-3-all', 'D_lc-I-0-all', 'D_lc-I-1-all', 'D_lc-I-2-all',
 'D_lc-I-3-all', 'D_lc-S-0-all', 'D_lc-S-1-all', 'D_lc-S-2-all',
 'D_lc-S-3-all', 'D_lc-T-0-all', 'D_lc-T-1-all', 'D_lc-T-2-all',
 'D_lc-T-3-all', 'D_lc-Z-0-all', 'D_lc-Z-1-all', 'D_lc-Z-2-all',
 'D_lc-Z-3-all', 'D_lc-chi-0-all', 'D_lc-chi-1-all', 'D_lc-chi-2-all',
 'D_lc-chi-3-all', 'D_mc-I-0-all', 'D_mc-I-1-all', 'D_mc-I-2-all',
 'D_mc-I-3-all', 'D_mc-S-0-all', 'D_mc-S-1-all', 'D_mc-S-2-all',
 'D_mc-S-3-all', 'D_mc-T-0-all', 'D_mc-T-1-all', 'D_mc-T-2-all',
 'D_mc-T-3-all', 'D_mc-Z-0-all', 'D_mc-Z-1-all', 'D_mc-Z-2-all',
 'D_mc-Z-3-all', 'D_mc-chi-0-all', 'D_mc-chi-1-all', 'D_mc-chi-2-all',
 'D_mc-chi-3-all', 'f-I-0-all', 'f-I-1-all', 'f-I-2-all', 'f-I-3-all',
 'f-S-0-all', 'f-S-1-all', 'f-S-2-all', 'f-S-3-all', 'f-T-0-all', 'f-T-1-all',
 'f-T-2-all', 'f-T-3-all', 'f-Z-0-all', 'f-Z-1-all', 'f-Z-2-all', 'f-Z-3-all',
 'f-chi-0-all', 'f-chi-1-all', 'f-chi-2-all', 'f-chi-3-all', 'f-lig-I-0',
 'f-lig-I-1', 'f-lig-I-2', 'f-lig-I-3', 'f-lig-S-0', 'f-lig-S-1', 'f-lig-S-2',
 'f-lig-S-3', 'f-lig-T-0', 'f-lig-T-1', 'f-lig-T-2', 'f-lig-T-3', 'f-lig-Z-0',
 'f-lig-Z-1', 'f-lig-Z-2', 'f-lig-Z-3', 'f-lig-chi-0', 'f-lig-chi-1',
 'f-lig-chi-2', 'f-lig-chi-3', 'func-I-0-all', 'func-I-1-all',
 'func-I-2-all', 'func-I-3-all', 'func-S-0-all', 'func-S-1-all',
 'func-S-2-all', 'func-S-3-all', 'func-T-0-all', 'func-T-1-all',
 'func-T-2-all', 'func-T-3-all', 'func-Z-0-all', 'func-Z-1-all',
 'func-Z-2-all', 'func-Z-3-all', 'func-chi-0-all', 'func-chi-1-all',
 'func-chi-2-all', 'func-chi-3-all', 'lc-I-0-all', 'lc-I-1-all', 'lc-I-2-all',
 'lc-I-3-all', 'lc-S-0-all', 'lc-S-1-all', 'lc-S-2-all', 'lc-S-3-all',
 'lc-T-0-all', 'lc-T-1-all', 'lc-T-2-all', 'lc-T-3-all', 'lc-Z-0-all',
 'lc-Z-1-all', 'lc-Z-2-all', 'lc-Z-3-all', 'lc-chi-0-all', 'lc-chi-1-all',
 'lc-chi-2-all', 'lc-chi-3-all', 'mc-I-0-all', 'mc-I-1-all', 'mc-I-2-all',
 'mc-I-3-all', 'mc-S-0-all', 'mc-S-1-all', 'mc-S-2-all', 'mc-S-3-all',
 'mc-T-0-all', 'mc-T-1-all', 'mc-T-2-all', 'mc-T-3-all', 'mc-Z-0-all',
 'mc-Z-1-all', 'mc-Z-2-all', 'mc-Z-3-all', 'mc-chi-0-all', 'mc-chi-1-all',
 'mc-chi-2-all', 'mc-chi-3-all']
geo = ['Df','Di', 'Dif','GPOAV','GPONAV','GPOV','GSA','POAV','POAV_vol_frac',
  'PONAV','PONAV_vol_frac','VPOV','VSA','cell_v']
 
# ANN_path = path + 'model/thermal/ANN/'
ANN_path = 'model/thermal/ANN/'
df_train_all = pd.read_csv(ANN_path+"train.csv").append(pd.read_csv(ANN_path+"val.csv"))
df_train = pd.read_csv(ANN_path+"train.csv")
df_train = df_train.loc[:, (df_train != df_train.iloc[0]).any()]
features = [val for val in df_train.columns.values if val in RACs+geo]

print(f'features is {features}')
with open('features.txt', 'w') as f:
    f.write(', '.join(features))

# Using a stand-in for df_newMOF
df_newMOF = df_train

X_train, X_newMOF, y_train, x_scaler, y_scaler = normalize_data_thermal(df_train, df_newMOF, features, ["T"], debug=False)
X_train.shape, y_train.reshape(-1, ).shape 

# Saving the x and y scalers as pkl files
with open('thermal_x_scaler.pkl', 'wb') as f:
    pkl.dump(x_scaler, f, protocol=pkl.HIGHEST_PROTOCOL)

with open('thermal_y_scaler.pkl', 'wb') as f:
    pkl.dump(y_scaler, f, protocol=pkl.HIGHEST_PROTOCOL)

# model = thermal_ANN

# new_MOF_pred = y_scaler.inverse_transform(model.predict(X_newMOF))
# new_MOF_pred = np.round(new_MOF_pred,1) # round to 1 decimal

# # isolating just the prediction, since the model spits out the prediction like [[PREDICTION]], as in, in hard brackets
# new_MOF_pred = new_MOF_pred[0][0]
# new_MOF_pred = str(new_MOF_pred)