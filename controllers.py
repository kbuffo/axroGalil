"""
Contains information about the controllers avialable. Setup assumes that all controllers are the same model.
"""

import gclib
import numpy as np

N_ao = 8 # number of analog outputs per controller
N_ai = 8 # number of analog inputs per controller
N_di = 16 # number of digital inputs per controller

ao_pins = np.arange(N_ao) # array of AO pins
ai_pins = np.arange(N_ai) # array of AI pins
di_pins = np.arange(N_di) # array of DI pins

controllers = {1: {'ip_address': '192.168.42.100', # avialable controllers
                   'gclib': gclib.py()}, 
               2: {'ip_address': '192.168.42.200', 
                   'gclib': gclib.py()}}

N_controllers = len(controllers) # number of available controllers

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
