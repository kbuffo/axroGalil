import numpy as np
import time
import os
from .controllers import controllers as ctrls
from .controllers import default_ai_cal_data_dir
from .controllers import default_ai_cal_fit_dir
from .controllers import default_ao_cal_data_dir
from .controllers import default_ao_cal_fit_dir
from .controllers import ao_min_v, ao_max_v
from .controllers import ai_min_v, ai_max_v
from .controllers import ao_pins as full_ao_pins
from .controllers import ai_pins as full_ai_pins
from .high_voltage_power_supply import psu_min_voltage as min_allowed_voltage
from .high_voltage_power_supply import psu_max_voltage as max_allowed_voltage
from .high_voltage_power_supply import psu_min_current as min_allowed_current
from .high_voltage_power_supply import psu_max_current as max_allowed_current
from .high_voltage_power_supply import psu_visa_address
from . import electronics as el
from . import psu_control as pcon
from . import amp_electronics as ampel

def calibrate_RIO_AO(cids, psu_ref_voltage, test_vrange=[0., ao_max_v], N_points=1000, ao_range=[0., ao_max_v], ao_pins=None, ai_range=[0., ai_max_v], ai_pins=None, 
                     iter_delay=2., samples_to_mean=20, psu_visa_address=psu_visa_address, ai_cal_fit_dir=default_ai_cal_fit_dir, ao_cal_data_dir=default_ao_cal_data_dir, 
                     ao_cal_fit_dir=default_ao_cal_fit_dir, printout=True):
    """
    Calibrates the analog outputs of a set of RIO controllers to read voltages from the controllers that have been amplified using a set of amplifiers and high voltage PSU.
    INPUTS:
    cids: list of ints -> specifies the controller IDs whose AO pins will be calibrated. See controllers.py for list of available controllers.
    psu_ref_voltage: float -> the reference voltage of PSU (in volts) used to amplify the AO voltage.
    test_vrange:list of two floats -> specifies the voltage range [min, max] (in volts) that will be set during calibration.
    N_points: int -> specifies the number of data points on psu_vrange to test.
    ao_range: list of two floats -> specifies the voltage range of the analog output pins in volts. For a list of the available voltage ranges the controllers
        can be set to,  see controllers.vranges
    ao_pins: (optional) list 1D arrays -> specifies thee analog output pins from each controller to calibrate. For example if cids=[1, 2] and 
        ao_pins=[[0, 1, 2, 3, 4, 5, 6, 7], [1, 4, 6]] then all pins of controller will be calibrated and pins [1, 4, 6] will be calibrated for controller 2. If None, all pins
        are tested.
    ai_range: list of two floats -> specifies the voltage range of the analog input pins in volts. For a list of the available voltage ranges the controllers
        can be set to, see controllers.vranges.
    ai_pins: (optional) list of 1D arrays -> specifies the analog input pins from each controller to calibrate. These pins are assumed to be given in the order such that 
        correspond to their AO pin connection given in ao_pins. If None, all pins are tested and it assumed that there is a one-to-one matching connection between AO and AI pins.
    iter_delay: (optional) float -> the number of seconds between when the voltage is set and the AIs read the voltage for data point tested. If None, there
        is no delay between setting and reading the voltages.
    samples_to_mean: int -> The number of AI voltage read samples to average per voltage tested.
    psu_visa_address: str -> the VISA address used to establish a connection session with the PSU. For the default VISA address used, see 
        high_voltage_power_supply.psu_visa_address
    ai_cal_fit_dir: str -> the directory where the AI pin calibration data is stored. This the same directory that cal_fit was saved to when running calibrate_RIO_AI()
    ao_cal_data_dir: (optional) str -> the directory in which to save cal_vdata. File names are automatically generated with the naming convention 
        "C{ID}_N_{N}_{pmin}-{pmax}V_to_{amin}-{amax}V_AO_cal_vdata.npy" where {ID} is the ID of controller that was calibrated, {N} is the value of N_points, 
        {pmin} and {pmax} are the values in test_vrange, and {amin} = 0 and {amax} is the voltage set by psu_ref_voltage. User will be prompted before altering existing 
        files in cal_data_dir with the new calibration data.
    ao_cal_fit_dir: (optional) str -> the directory in which to save cal_fit. File names are automatically generated with the naming convention 
        "C{ID}_N_{N}_{pmin}-{pmax}V_to_{amin}-{amax}V_AO_cal_fit.npy" where {ID} is the ID of controller that was calibrated, {N} is the value of N_points, 
        {pmin} and {pmax} are the values in test_vrange, and {amin} = 0 and {amax} is the voltage set by psu_ref_voltage. User will be prompted before altering existing 
        files in cal_fit_dir with the new calibration fit.
    printout: bool -> if True, prints the result of each test iteration.
    RETURNS:
    cal_vdata: List of 2D arrays -> The calibration write and read voltages. cal_vdata[i] specifies the voltage array for cids[i], and in 
        the corresponding array:
            * Has shape (P+1, N_points), P is the total number of AI pins that the controller has (specified by full_ai_pins).
            * The first row is the voltages set by the controller.
            * The last P rows are the voltages read by full_ai_pins. Any pins not tested are set to NaN unless altering an existing file.
    cal_fit: List of 2D arrays -> The calibration linear fits. cal_fit[i] specifies the fit for cids[i], and in the corresponding array:
            * Has shape (P, 2), P is the total number of AI pins that the controller has (specified by full_ai_pins).
            * The first P rows index the fits by full_ai_pins. Any pins not tested are set to NaN unless altering an existing file.
            * The first column is the slope of the fit, the second column is the intercept. The fit computes
                (AO pin voltage) = (Amplified voltage)*(AO pin fit slope) + (AO pin fit intercept)
    """
    if ao_pins is None:
        ao_pins = [full_ao_pins for i in range(len(cids))]
    if ai_pins is None:
        ai_pins = [full_ai_pins for i in range(len(cids))]
    for i in range(len(cids)): # check that ao_pins and ai_pins match in length for each controller
        if len(ao_pins[i]) != len(ai_pins[i]):
            raise CalibrationError(" Controller ({}) ao_pins ({}) do not match length of ai_pins ({}).".format(cids[i], ao_pins[i], ai_pins[i]))
    # initialize data products and filenames
    ao_cal_vdata, ao_cal_fit, ai_cal_fit, ao_cal_data_fnames, ao_cal_fit_fnames = get_ao_cal_props(cids, psu_ref_voltage, test_vrange, N_points, 
                                                                                                   ao_range, ao_cal_data_dir, ao_cal_fit_dir, 
                                                                                                   ai_cal_fit_dir)
    if ao_cal_vdata is None or ao_cal_fit is None:
        print('Calibration cancelled.')
        return None, None
    test_volts = np.linspace(*test_vrange, num=N_points)
    # initialize PSU
    for cid in cids: # initialize controllers
        el.open_controller(ctrls[cid]['gclib'], ctrls[cid]['ip_address'])
        for ao_pin in full_ao_pins:
            el.set_controller_ao_vrange(ctrls[cid]['gclib'], ao_pin, vrange=ao_range)
        for ai_pin in full_ai_pins:
            el.set_controller_ai_vrange(ctrls[cid]['gclib'], ai_pin, vrange=ai_range)
        print('Open connection to controller: {}, AO voltage range: {} V, AI voltage range: {} V'.format(cid, ao_range, ai_range))
    print('Calibrating...')
    if printout:
        print('='*80)
    for i, tvolt in enumerate(test_volts): # test each voltage
        if printout:
            print('-'*10+'Testing: {:.2f} V, (test {} of {}, {:.1f}% complete)'.format(tvolt, i+1, len(test_volts), (i/len(test_volts))*100.)+'-'*10)
        for j, cid in enumerate(cids): # set the test voltage on the pins of the controllers
            for k, ao_pin in enumerate(full_ao_pins):
                if ao_pin in ao_pins[j]:
                    el.set_ao(ctrls[cid]['gclib'], int(ao_pin), tvolt)
        if iter_delay is not None:
            time.sleep(iter_delay)
        # get the voltage of each AI pin of each controller
        for j, cid in enumerate(cids):
            for k, ai_pin in enumerate(full_ai_pins):
                if ai_pin in ai_pins[j]:
                    ai_voltage = np.mean([ampel.read_amp_ai(cid, ai_pin, slope=ai_cal_fit[j][k][0], intercept=ai_cal_fit[j][k][1]) 
                                          for n in range(samples_to_mean)])
                    l = np.where(ai_pins[j]==ai_pin)[0][0] # index of ai pin
                    m = np.where(full_ao_pins==ao_pins[j][l])[0][0] # index of ao pin
                    matching_ao_pin = full_ao_pins[m]
                    ao_cal_vdata[j][m+1][i] = ai_voltage
                    if printout:
                        print('Controller {}, AO {}, AI {} reads {} V'.format(cid, matching_ao_pin, ai_pin, ai_voltage))
    if printout:
        print('='*80)
    for cid in cids: # close controllers
        for ao_pin in full_ao_pins:
            el.set_ao(ctrls[cid]['gclib'], int(ao_pin), 0.)
        el.close_controller(ctrls[cid]['gclib'])
        print('Closed connection to controller: {}'.format(cid))
    print('Calibration testing complete.')
    # write the PSU voltages to cal_vdata
    for vdata_arr in ao_cal_vdata:
        vdata_arr[0] = test_volts
    # calcluate the linear fits of AI pin
    ao_cal_fit = calc_ao_fit(cids, ao_pins, ao_cal_vdata, ao_cal_fit)
    # save data
    if all(var is not None for var in [ao_cal_vdata, ao_cal_fit, ao_cal_data_dir]):
        for i in range(len(cids)):
            np.save(ao_cal_data_dir+ao_cal_data_fnames[i], ao_cal_vdata[i])
        print('Calibration voltage data saved.')
    if all(var is not None for var in [ao_cal_vdata, ao_cal_fit, ao_cal_fit_dir]):
        for i in range(len(cids)):
            np.save(ao_cal_fit_dir+ao_cal_fit_fnames[i], ao_cal_fit[i])
        print('Calibration fits data saved.')
    return ao_cal_vdata, ao_cal_fit

def calibrate_RIO_AI(cids, psu_vrange=[.2, 310.], N_points=1000, ai_range=[0., ai_max_v], ai_pins=None, iter_delay=10., 
                     samples_to_mean=20, psu_visa_address=psu_visa_address, cal_data_dir=default_ai_cal_data_dir, cal_fit_dir=default_ai_cal_fit_dir,
                     readback_psu_voltage=True, printout=True):
    """
    Calibrates the analog inputs of a set of RIO controllers to read voltages from the PSU that have been stepped down using voltage dividers.
    INPUTS:
    cids: list of ints -> specifies the controller IDs whose AI pins will be calibrated. See controllers.py for list of available controllers.
    psu_vrange:list of two floats -> specifies the voltage range [min, max] in volts of the PSU that will be set during calibration.
    N_points: int -> specifies the number of data points on psu_vrange to test.
    ai_range: list of two floats -> specifies the voltage range of the analog input pins in volts. For a list of the available voltage ranges the controllers
        can be set to, see controllers.vranges.
    ai_pins: (optional) list of 1D arrays -> specifies the analog input pins from each controller to calibrate. For example if cids=[1, 2] and 
        ai_pins=[[0, 1, 2, 3, 4, 5, 6, 7], [1, 4, 6]] then all pins of controller 1 will be calibrated and pins [1, 4, 6] will be calibrated for controller 2.
        ai_pins are always sorted into numerical order. If None, all pins for cids are calibrated.
    iter_delay: (optional) float -> the number of seconds between when the PSU voltage is set and the AIs read the voltage for data point tested. If None, there
        is no delay between setting and reading the voltages.
    samples_to_mean: int -> The number of AI voltage read samples to average per voltage tested.
    psu_visa_address: str -> the VISA address used to establish a connection session with the PSU. For the default VISA address used, see 
        high_voltage_power_supply.psu_visa_address
    cal_data_dir: (optional) str -> the directory in which to save cal_vdata. File names are automatically generated with the naming convention 
        "C{ID}_N_{N}_{pmin}-{pmax}V_to_{amin}-{amax}V_AI_cal_vdata.npy" where {ID} is the ID of controller that was calibrated, {N} is the value of N_points, 
        {pmin} and {pmax} are the values in psu_vrange, and {amin} and {amax} are the values in ai_range. User will be prompted before altering existing 
        files in cal_data_dir with the new calibration data.
    cal_fit_dir: (optional) str -> the directory in which to save cal_fit. File names are automatically generated with the naming convention 
        "C{ID}_N_{N}_{pmin}-{pmax}V_to_{amin}-{amax}V_AI_cal_fit.npy" where {ID} is the ID of controller that was calibrated, {N} is the value of N_points, 
        {pmin} and {pmax} are the values in psu_vrange, and {amin} and {amax} are the values in ai_range. User will be prompted before altering existing 
        files in cal_fit_dir with the new calibration fit.
    readback_psu_voltage: bool -> if True, when testing voltage in psu_vrange, the PSU will read back the voltage that was set and use that in the 
        calibration rather than assuming that the test value that was sent to it is the true value of the PSU. NOTE: The Keysight N5771A can't internally 
        read voltages below 0.15 V, hence if this is set to True, the recommended minimum PSU voltage to test is 0.2. 
    printout: bool -> if True, prints the result of each test iteration.
    RETURNS:
    cal_vdata: List of 2D arrays -> The calibration write and read voltages. cal_vdata[i] specifies the voltage array for cids[i], and in 
        the corresponding array:
            * Has shape (P+1, N_points), P is the total number of AI pins that the controller has (specified by full_ai_pins).
            * The first row is the voltages set by the PSU.
            * The last P rows are the voltages read by full_ai_pins. Any pins not tested are set to NaN unless altering an existing file.
    cal_fit: List of 2D arrays -> The calibration linear fits. cal_fit[i] specifies the fit for cids[i], and in the corresponding array:
            * Has shape (P, 2), P is the total number of AI pins that the controller has (specified by full_ai_pins).
            * The first P rows index the fits by full_ai_pins. Any pins not tested are set to NaN unless altering an existing file.
            * The first column is the slope of the fit, the second column is the intercept. The fit computes 
                (PSU_voltage) = (AI pin voltage)*(AI pin fit slope) + (AI pin fit intercept)
    """
    cids = sorted(cids)
    if ai_pins is None:
        ai_pins = [full_ai_pins for i in range(len(cids))]
    else:
        ai_pins = [np.sort(ai_pins[i]) for i in range(len(ai_pins))]
    # initialize data products and filenames
    cal_vdata, cal_fit, cal_data_fnames, cal_fit_fnames = get_ai_cal_props(cids, psu_vrange, N_points, ai_range, cal_data_dir, cal_fit_dir)
    if cal_vdata is None or cal_fit is None:
        print('Calibration cancelled.')
        return None, None
    test_volts = np.linspace(psu_vrange[0], psu_vrange[1], num=N_points)
    psu_volts = np.zeros(test_volts.shape)
    # initialize PSU
    psu, rm = pcon.open_psu_connection(visa_address=psu_visa_address, init_psu=True, printout=True)
    pcon.enable_psu_output(psu, printout=True)
    for cid in cids: # initialize controllers
        el.open_controller(ctrls[cid]['gclib'], ctrls[cid]['ip_address'])
        for ai_pin in full_ai_pins:
            el.set_controller_ai_vrange(ctrls[cid]['gclib'], ai_pin, vrange=ai_range)
        print('Open connection to controller: {}, AI voltage range: {} V'.format(cid, ai_range))
    print('Calibrating...')
    if printout:
        print('='*80)
    for i, tvolt in enumerate(test_volts): # test each voltage
        pcon.set_psu_voltage(psu, tvolt, readback=False)
        if iter_delay is not None:
            time.sleep(iter_delay)
        # get the PSU voltage
        if readback_psu_voltage:
            psu_volts[i] = pcon.read_psu_voltage(psu) #np.mean([pcon.read_psu_voltage(psu) for n in range(samples_to_mean)])
        else:
            psu_volts[i] = tvolt
        if printout:
            print('-'*10+'Testing: {:.2f} V, PSU: {:.2f} V, (test {} of {}, {:.0f}% complete)'.format(tvolt, psu_volts[i], i+1, len(test_volts), (i/len(test_volts))*100.)+'-'*10)
        # get the voltage of each AI pin of each controller
        for j, cid in enumerate(cids):
            for k, ai_pin in enumerate(full_ai_pins):
                if ai_pin in ai_pins[j]:
                    ai_voltage = np.mean([el.read_ai(ctrls[cid]['gclib'], int(ai_pin)) for n in range(samples_to_mean)])
                    cal_vdata[j][k+1][i] = ai_voltage
                    if printout:
                        print('Controller {}, AI {} reads {} V'.format(cid, ai_pin, ai_voltage))
    if printout:
        print('='*80)
    pcon.reset_psu(psu, printout=True)
    pcon.get_psu_output_state(psu)
    for cid in cids: # close controllers
        el.close_controller(ctrls[cid]['gclib'])
        print('Closed connection to controller: {}'.format(cid))
    print('Calibration testing complete.')
    # write the PSU voltages to cal_vdata
    for vdata_arr in cal_vdata:
        vdata_arr[0] = psu_volts
    # calcluate the linear fits of AI pin
    cal_fit = calc_ai_fit(cids, ai_pins, cal_vdata, cal_fit)
    # save data
    if all(var is not None for var in [cal_vdata, cal_fit, cal_data_dir]):
        for i in range(len(cids)):
            np.save(cal_data_dir+cal_data_fnames[i], cal_vdata[i])
        print('Calibration voltage data saved.')
    if all(var is not None for var in [cal_vdata, cal_fit, cal_fit_dir]):
        for i in range(len(cids)):
            np.save(cal_fit_dir+cal_fit_fnames[i], cal_fit[i])
        print('Calibration fits data saved.')
    return cal_vdata, cal_fit

def calc_ao_fit(cids, ao_pins, cal_vdata, cal_fit, fit_deg=1):
    """
    Calculates the linear fits of the AO pins tested. cal_vdata and cal_fit are defined in calibrate_RIO_AO().
    Returns cal_fit with the linear fit values (slope, intercept) appended, or None in the event fitting is 
    cancelled.
    """
    nan_datasets = []
    # verify integrity of the AO and AI data
    for i, cid in enumerate(cids):
        print('C{} cal_vdata:\n{}'.format(cid, cal_vdata[i]))
        if np.any(np.isnan(cal_vdata[i][0])): # Bad controller AO data
            nan_datasets.append('C{}: AO VOLTAGES'.format(cid))
        for j, ao_pin in enumerate(full_ao_pins):
            if ao_pin in ao_pins[i] and np.any(np.isnan(cal_vdata[i][j+1])): # Bad AI data
                print('AO {} is in {}'.format(ao_pin, ao_pins[i]))
                nan_datasets.append('C{}: AI{} VOLTAGES'.format(cid, ao_pin))
    if nan_datasets:
        prompt = 'WARNING the following datasets contain NaN values after testing:\n--- '+'\n--- '.join(nan_datasets)+'\n Proceed with data fitting? (y/n): '
        error_prompt = "Error: enter 'y' to proceed with fitting or 'n' to cancel."
        if not yesno_prompt(prompt, error_prompt): # Cancel fitting
            print('Fitting cancelled.')
            return None
    for i in range(len(cids)): # compute linear fits for AO pins tested
        for j, ao_pin in enumerate(full_ao_pins):
            if ao_pin in ao_pins[i]:
                cal_fit[i][j] = np.polyfit(cal_vdata[i][j+1], cal_vdata[i][0], fit_deg)
    return cal_fit

def calc_ai_fit(cids, ai_pins, cal_vdata, cal_fit, fit_deg=1):
    """
    Calculates the linear fits of the AI pins tested. cal_vdata and cal_fit are defined in calibrate_RIO_AI().
    Returns cal_fit with the linear fit values (slope, intercept) appended, or None in the event fitting is cancelled.
    """
    nan_datasets = []
    # verify the integrity of the PSU and AI data
    for i, cid in enumerate(cids):
        if np.any(np.isnan(cal_vdata[i][0])): # Bad controller PSU data
            nan_datasets.append('C{}: PSU VOLTAGES'.format(cid))
        for j, ai_pin in enumerate(full_ai_pins):
            if ai_pin in ai_pins[i] and np.any(np.isnan(cal_vdata[i][j+1])): # Bad AI data
                nan_datasets.append('C{}: AI{} VOLTAGES'.format(cid, ai_pin))
    if nan_datasets:
        prompt = 'WARNING the following datasets contain NaN values after testing:\n--- '+'\n--- '.join(nan_datasets)+'\n Proceed with data fitting? (y/n): '
        error_prompt = "Error: enter 'y' to proceed with fitting or 'n' to cancel."
        if not yesno_prompt(prompt, error_prompt): # Cancel fitting
            print('Fitting cancelled.')
            return None
    for i in range(len(cids)): # compute linear fits for AI pins tested.
        for j, ai_pin in enumerate(full_ai_pins):
            if ai_pin in ai_pins[i]:
                cal_fit[i][j] = np.polyfit(cal_vdata[i][j+1], cal_vdata[i][0], fit_deg)
    return cal_fit

def get_ao_cal_props(cids, psu_ref_voltage, test_vrange, N_points, ao_range, ao_cal_data_dir, ao_cal_fit_dir, ai_cal_fit_dir):
    """
    Generates filenames and data arrays for calibrating controller AO pins.
    """
    # Check for existing calibration files.
    existing_files = []
    if ao_cal_data_dir is not None: # saving cal data files
        # cal data files that will be saved
        ao_cal_data_fnames = ['C{}_N_{}_{}-{}V_to_{}-{}V_AO_cal_vdata'.format(cid, N_points, *test_vrange, 0., psu_ref_voltage) 
                              for cid in cids]
        # cal data files that exist and (might) be altered
        existing_ao_cal_data_files = [ao_cal_data_fname for ao_cal_data_fname in ao_cal_data_fnames 
                                      if os.path.isfile(ao_cal_data_dir+ao_cal_data_fname+'.npy')]
        if existing_ao_cal_data_files: # add to list of existing files
            [existing_files.append(fname) for fname in existing_ao_cal_data_files]
    else: # Not saving cal data files
        ao_cal_data_fnames = None
    if ao_cal_fit_dir is not None: # saving cal fit files
        # cal fit files that will be saved
        ao_cal_fit_fnames = ['C{}_N_{}_{}-{}V_to_{}-{}V_AO_cal_fit'.format(cid, N_points, *test_vrange, 0., psu_ref_voltage) 
                             for cid in cids]
        # cal fit files that exist and (might) be altered
        existing_ao_cal_fit_files = [ao_cal_fit_fname for ao_cal_fit_fname in ao_cal_fit_fnames 
                                     if os.path.isfile(ao_cal_fit_dir+ao_cal_fit_fname+'.npy')]
        if existing_ao_cal_fit_files: # add to list of existing files
            [existing_files.append(fname) for fname in existing_ao_cal_fit_files]
    else:
        ao_cal_fit_fnames = None

    if existing_files: # If there are existing files
        existing_files_print_str = '.npy\n---'.join([fname for fname in existing_files]) + '.npy\n'
        prompt = 'WARNING: The following calibration files will be altered:\n---'+existing_files_print_str+'\nProceed? (y/n): '
        error_prompt = "Error: enter 'y' to proceed with calibration or 'n' to cancel."
        if yesno_prompt(prompt, error_prompt): # Load existing files to alter them.
            cal_vdata = []
            for fname in ao_cal_data_fnames:
                if fname in existing_ao_cal_data_files:
                    cal_vdata_array = np.load(ao_cal_data_dir+fname+'.npy')
                else:
                    cal_vdata_array = np.full((len(full_ai_pins)+1, N_points), np.nan)
                cal_vdata.append(cal_vdata_array)
            cal_fit = []
            for fname in ao_cal_fit_fnames:
                if fname in existing_ao_cal_fit_files:
                    cal_fit_array = np.load(ao_cal_fit_dir+fname+'.npy')
                else:
                    cal_fit_array = np.fulll((len(full_ai_pins), 2), np.nan)
                cal_fit.append(cal_fit_array)
        else: # Cancel calibration
            cal_vdata = None
            cal_fit = None
    else: # If there are no existing files
        cal_vdata = [np.full((len(full_ai_pins)+1, N_points), np.nan) for i in range(len(cids))]
        cal_fit = [np.full((len(full_ai_pins), 2), np.nan) for i in range(len(cids))]
    
    # Find controller AI calibration files
    matching_ai_cal_fit_files = [] # load the AI calibration fits
    ai_cal_files = os.listdir(ai_cal_fit_dir)
    #print('ai_cal_fit_files:\n{}'.format(ai_cal_fit_files))
    for cid in cids:
        found_match = False
        ai_cal_fit_files = [fname for fname in ai_cal_files if fname.startswith('C{}'.format(cid)) 
                        and fname.endswith('cal_fit.npy')]
        for fname in ai_cal_fit_files:
            print('current fname to check: {}'.format(fname))
            if fname.startswith('C{}'.format(cid)) and fname.endswith('cal_fit.npy'):
                if not found_match:
                    matching_ai_cal_fit_files.append(fname)
                else:
                    raise CalibrationError('In AI calibration directory:\n{},\nMultiple AI calibration fit files found for (controller {}). Ensure only one calibration fit file per controller in ai_cal_fit_dir.'.format(ai_cal_fit_dir, cid))
                found_match = True
            else:
                raise CalibrationError('In AI calibration directory:\n{},\nController ({}) calibration fit file not found. Calibration fit file should be named with the signature C{}_N_{}_{}-{}V_to_{}-{}V_AI_cal_fit.npy'.format(ai_cal_fit_dir, cid, cid, *(['{}']*5)))
    cal_ai_fit = [np.load(ai_cal_fit_dir+fname) for fname in matching_ai_cal_fit_files]
    return cal_vdata, cal_fit, cal_ai_fit, ao_cal_data_fnames, ao_cal_fit_fnames



def get_ai_cal_props(cids, psu_vrange, N_points, ai_range, cal_data_dir, cal_fit_dir):
    """
    Generates filenames and data arrays for calibrating controller AI pins.
    """

    # Check for existing calibration files.
    existing_files = []
    if cal_data_dir is not None: # saving cal data files
        # cal data files that will be saved
        cal_data_fnames = ['C{}_N_{}_{}-{}V_to_{}-{}V_AI_cal_vdata'.format(cid, N_points, psu_vrange[0], psu_vrange[1], ai_range[0], ai_range[1]) 
                           for cid in cids]
        # cal data files that exist and (might) be altered
        existing_cal_data_files = [cal_data_fname for cal_data_fname in cal_data_fnames if os.path.isfile(cal_data_dir+cal_data_fname+'.npy')]
        if existing_cal_data_files: # add to list of existing files
            [existing_files.append(fname) for fname in existing_cal_data_files]
    else: # Not saving cal data files
        cal_data_fnames = None
    if cal_fit_dir is not None: # saving cal fit files
        # cal fit files that will be saved
        cal_fit_fnames = ['C{}_N_{}_{}-{}V_to_{}-{}V_AI_cal_fit'.format(cid, N_points, psu_vrange[0], psu_vrange[1], ai_range[0], ai_range[1]) 
                          for cid in cids]
        # cal fit files that exist and (might) be altered
        existing_cal_fit_files = [cal_fit_fname for cal_fit_fname in cal_fit_fnames if os.path.isfile(cal_fit_dir+cal_fit_fname+'.npy')]
        if existing_cal_fit_files: # add to list of existing files
            [existing_files.append(fname) for fname in existing_cal_fit_files]
    else: # Not saving cal fit files
        cal_fit_fnames = None

    if existing_files: # If there are existing files
        existing_files_print_str = '.npy\n--- '.join([fname for fname in existing_files]) + '.npy\n'
        prompt = 'WARNING: The following calibration files will be altered:\n--- '+existing_files_print_str+'\nProceed? (y/n): '
        error_prompt = "Error: enter 'y' to proceed with calibration or 'n' to cancel."
        if yesno_prompt(prompt, error_prompt): # Load existing files to alter them.
            cal_vdata = []
            for fname in cal_data_fnames:
                if fname in existing_cal_data_files:
                    cal_vdata_array = np.load(cal_data_dir+fname+'.npy')
                else:
                    cal_vdata_array = np.full((len(full_ai_pins)+1, N_points), np.nan)
                cal_vdata.append(cal_vdata_array)
            cal_fit = []
            for fname in cal_fit_fnames:
                if fname in existing_cal_fit_files:
                    cal_fit_array = np.load(cal_fit_dir+fname+'.npy')
                else:
                    cal_fit_array = np.full((len(full_ai_pins), 2), np.nan)
                cal_fit.append(cal_fit_array)
        else: # Cancel calibration
            cal_vdata = None
            cal_fit = None
    else: # If there are no existing files
        cal_vdata = [np.full((len(full_ai_pins)+1, N_points), np.nan) for i in range(len(cids))]
        cal_fit = [np.full((len(full_ai_pins), 2), np.nan) for i in range(len(cids))]
    return cal_vdata, cal_fit, cal_data_fnames, cal_fit_fnames

def yesno_prompt(prompt, error_prompt):
    """
    Function that writes prompt to user. Returns True if user inputs "y" and False if user inputs "n".
    """
    if not error_prompt:
        error_prompt = "Error: please enter 'y' or 'n': "
    while True:
        user_input = input(prompt).strip().lower()
        if user_input == 'y':
            return True
        elif user_input == 'n':
            return False
        else:
            print(error_prompt)

class CalibrationError(Exception):
    """
    Custom exception related to calibrating the amp system
    """
    pass