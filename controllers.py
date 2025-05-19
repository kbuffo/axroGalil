"""
Contains information about the controllers avialable. Setup assumes that all controllers are the same model.
"""

import os
import gclib
import numpy as np

N_ao = 8 # number of analog outputs per controller
N_ai = 8 # number of analog inputs per controller
N_di = 16 # number of digital inputs per controller

ao_pins = np.arange(N_ao) # array of AO pins
ai_pins = np.arange(N_ai) # array of AI pins
di_pins = np.arange(N_di) # array of DI pins

ao_min_v = -10. # minimum voltage for ao pins
ao_max_v = 10. # maximum voltage for ao pins
ai_min_v = -10. # minimum voltage for ai pins
ai_max_v = 10. # maximum voltage for ai pins
# configurable voltage output and input ranges and the associated GCommand code used to set the
# controller voltage range
vranges = {'AO': {(ao_min_v, ao_max_v): 4,
                  (0., ao_max_v): 2}, 
           'AI': {(ai_min_v, ai_max_v): 2,
                  (0., ai_max_v): 4}}
di_min_v = 0. # minimum voltage for di pins
di_max_v = 24. # maximum voltage for di pins
model = 'RIO-47142' # controller model


# avialable controllers
controllers = {1: {'ip_address': '192.168.42.100', # IP address used to connect to controller
                   'gclib': gclib.py(),            # unique instance of gclib used to communicate with controller
                   'ai_cal_fit': None,             # The calibration fit for the AI pins of the controller can be manually loaded here               
                   'ai_cal_data':None,             # The calibration data for the AI pins of the controller can be manually loaded here
                   'ao_cal_fit': None,             # The calibration fit for the AO pins of the controller can be manually loaded here
                   'ao_cal_data':None,},           # The calibration data for the AO pins of the controller can be manually loaded here
               2: {'ip_address': '192.168.42.200', 
                   'gclib': gclib.py(),
                   'ai_cal_fit': None,
                   'ai_cal_data': None,
                   'ao_cal_fit': None, 
                   'ao_cal_data': None},
               3: {'ip_address': '192.168.42.30', 
                   'gclib': gclib.py(),
                   'ai_cal_fit': None,
                   'ai_cal_data': None,
                   'ao_cal_fit': None, 
                   'ao_cal_data': None},
               4: {'ip_address': '192.168.42.40', 
                   'gclib': gclib.py(),
                   'ai_cal_fit': None,
                   'ai_cal_data': None,
                   'ao_cal_fit': None, 
                   'ao_cal_data': None}
               }

N_controllers = len(controllers) # number of available controllers

"""
The following parameters and file signatures can be altered to automatically load the AI/AO calibration data for each controller.
"""

# directories for storing the AI and AO calibration fits and data for writing and reading amplified voltages
repo_dir = os.path.dirname(os.path.abspath(__file__)) + '\\'
default_ai_cal_fit_dir = repo_dir + 'calibration_data\\AI_calibration\\'
default_ai_cal_data_dir = repo_dir + 'calibration_data\\AI_calibration\\'
default_ao_cal_fit_dir = repo_dir + 'calibration_data\\AO_calibration\\'
default_ao_cal_data_dir = repo_dir + 'calibration_data\\AO_calibration\\'

N_cal_points = 6000 # The number of test points used during AI and AO calibration
ai_cal_write_vrange = [0.2, 310.0] # The voltage range of the PSU during AI calibration
ai_cal_read_vrange = [0.0, 10.0] # The voltage range of that the AI pins read during AI calibration
ao_cal_write_vrange = [0.0, 10.0] # The voltage range of the AO pins during AO calibration
ao_cal_read_vrange = [0.0, 310.0] # The voltage range of the that the AI pins read during AO calibration
# File signature for AI calibration fit to automatically load for each controller
ai_cal_fit_file_sig = 'C{}_'+'N_{}_{}-{}V_to_{}-{}V_AI_cal_fit.npy'.format(N_cal_points, *ai_cal_write_vrange, *ai_cal_read_vrange)
# File signature for AI calibration data to automatically load for each controller
ai_cal_data_file_sig = 'C{}_'+'N_{}_{}-{}V_to_{}-{}V_AI_cal_vdata.npy'.format(N_cal_points, *ai_cal_write_vrange, *ai_cal_read_vrange)
# File signature for AO calibration data to automatically load for each controller
ao_cal_fit_file_sig = 'C{}_'+'N_{}_{}-{}V_to_{}-{}V_AO_cal_fit.npy'.format(N_cal_points, *ao_cal_write_vrange, *ao_cal_read_vrange)
# File signature for AO calibration data to automatically load for each controller
ao_cal_data_file_sig = 'C{}_'+'N_{}_{}-{}V_to_{}-{}V_AO_cal_vdata.npy'.format(N_cal_points, *ao_cal_write_vrange, *ao_cal_read_vrange)

def autoload_controller_cal_data(cids, ai_cal_fit_file_sig, ai_cal_data_file_sig, ao_cal_fit_file_sig, ao_cal_data_file_sig):
    """
    Automatically loads the AI and AO calibration data for each avialable controller assuming that each controller has been calibrated under
    the same voltage constraints.
    cids: list of ints -> The controller IDs for the controllers whose calibration data will be loaded.
    ai_cal_fit_file_sig: str -> The path and file signature for the the AI calibration fit
    ai_cal_data_file_sig: str -> The path and file signature for the AI calibration data
    ao_cal_fit_file_sig: str -> The path and file signature for the AO calibration fit
    ao_cal_data_file_sig: str -> The path and file signature for the AO calibration data
    If any of these file signatures are None, that data will not be loaded.
    """
    for cid in cids:
        for controller_prop, file_sig in zip(['ai_cal_fit', 'ai_cal_data', 'ao_cal_fit', 'ao_cal_data'], 
                                             [ai_cal_fit_file_sig, ai_cal_data_file_sig, ao_cal_fit_file_sig, ao_cal_data_file_sig]):
            if file_sig is not None:
                controllers[cid][controller_prop] = np.load(file_sig.format(cid))

# This will auto load the calibration data when this script is imported
#autoload_controller_cal_data(list(controllers.keys()), 
#                             default_ai_cal_fit_dir+ai_cal_fit_file_sig, 
#                             default_ai_cal_data_dir+ai_cal_data_file_sig, 
#                             None, 
#                             None)