"""
Allows setting and reading amplified voltages through combined use of the controllers and high voltage PSU
"""

from .high_voltage_power_supply import psu_min_voltage
from .high_voltage_power_supply import psu_max_voltage
from .high_voltage_power_supply import psu_min_current
from .high_voltage_power_supply import psu_max_current
from .high_voltage_power_supply import psu_visa_address
from .controllers import controllers as ctrl
from .controllers import ao_min_v, ao_max_v
from .controllers import ai_min_v, ai_max_v
from .controllers import ao_pins as full_ao_pins
from .controllers import ai_pins as full_ai_pins
from .controllers import vranges as hardware_vranges
import axroGalil.psu_control as pcon
import axroGalil.electronics as el
import numpy as np
import time

write_read_delay = 2 # seconds, the time to wait after setting an amplified voltage before reading it back

####################### CELL CONTROL FUNCTIONS #######################

def open_cell_amp_system(cells, psu_visa_address=psu_visa_address, init_psu=True, 
                         set_psu_current=psu_max_current, psu_scale_voltage=300., 
                         enable_psu_output=False, ai_unscaled_vrange=[0., ai_max_v], 
                         ao_unscaled_vrange=[0., ao_max_v], ground_ao=True, printout=True):
    """
    Wrapper function for open_amp_system() that takes in a list of cell objects and configures the
    PSU and appropriate controllers.
    """
    cids = el.get_controller_ids(cells)
    psu, rm = open_amp_system(cids, psu_visa_address=psu_visa_address, init_psu=init_psu, 
                              set_psu_current=set_psu_current, psu_scale_voltage=psu_scale_voltage, 
                              enable_psu_output=enable_psu_output, ai_unscaled_vrange=ai_unscaled_vrange, 
                              ao_unscaled_vrange=ao_unscaled_vrange, ground_ao=ground_ao, printout=printout)
    return psu, rm

def close_cell_amp_system(cells, psu_session, resource_manager, reset_psu=True, ground_ao=True,
                          printout=True):
    """
    Wrapper function for close_amp_system() that takes in a list of cell objects and closes the
    connection to the PSU and controllers.
    """
    cids = el.get_controller_ids(cells)
    close_amp_system(cids, psu_session, resource_manager, reset_psu=reset_psu, ground_ao=ground_ao,
                     printout=True)

def setAmpCell(cell_no, cells, voltage, readback=True, delay=write_read_delay, printout=True):
    """
    Sets the amplified voltage on a cell. If readback is True, readCell() is run on the same cell and 
    cell.voltage and cell.volt_hist are update with the voltage that is read back. If readback is True,
    delay specifies how long to wait (in seconds) after setting the voltage to read it back.
    cell_no: int -> cell number to control
    cells: list -> list of cell objects to read from
    voltage: float -> voltage to set (in volts)
    """
    cell = el.filter_cells(cells, 'no', cell_no)[0]
    set_amp_ao(cell.cid, cell.ao_pin, voltage)
    if readback:
        read_voltage = readAmpCell(cell_no, cells, delay=delay)
        cell.voltage = read_voltage
        cell.volt_hist.append(read_voltage)
        if printout:
            print('Cell {} is now at {} V'.format(cell.no, cell.voltage))
            
def readAmpCell(cell_no, cells, delay=None):
    """
    Returns the amplified voltage for a cell number cell_no in a given list of cells.
    cell_no: int ->  the cell number to read.
    cells: list -> the list of cell objects to read from.
    delay: (optional) float -> the time (in seconds) to wait before reading the voltage. 
    """
    if delay is not None:
        time.sleep(delay)
    cell = el.filter_cells(cells, 'no', cell_no)[0]
    voltage = read_amp_ai(cell.cid, cell.ai_pin)
    return voltage

####################### DIRECT SYSTEM CONTROL FUNCTIONS #######################

def open_amp_system(cids, psu_visa_address=psu_visa_address, init_psu=True, 
                    set_psu_current=psu_max_current, psu_scale_voltage=300., 
                    enable_psu_output=False, ai_unscaled_vrange=[0., ai_max_v], 
                    ao_unscaled_vrange=[0., ao_max_v], ground_ao=True, printout=True):
    """
    Opens connection to PSU and controllers.
    cids: list of ints -> list of controller IDs whose controllers will be connected to.
    psu_visa_address: str -> the VISA address of the PSU to connect to.
    init_psu: bool -> If True, the PSU is set to its reset state.
    set_psu_current: (optional) float -> If not None and init_psu=True, then the PSU current
        will be set to the set_psu_current value (in Amps).
    psu_scale_voltage: (optional) float -> If not None, the PSU will be set to this voltage as the 
        reference voltage for the amplifiers (in Volts).
    enable_psu_output: bool -> If true, the output of the PSU will be enabled.
    ai_unscaled_vrange: (optional) list of floats -> If not None, all AI pins of the controllers dictated by 
        cids will be set to this voltage range.
    ao_unscaled_vrange: (optional) list of floats -> If not None, all AO pins of the controllers dictated by
        cids will be set to this voltage range.
    ground_ao: bool -> If True, then all AO pins of all cids will be set to 0 V.
    Returns:
    psu: an instance of the PSU connection session
    rm: an instance of the resource manager which tracks all connected instruments
        through the pyvisa API.
    """
    psu, rm = pcon.open_psu_connection(visa_address=psu_visa_address, init_psu=init_psu, 
                                       set_current=set_psu_current, printout=printout)
    if psu_scale_voltage is not None:
        pcon.set_psu_voltage(psu, psu_scale_voltage, readback=False)
    for cid in cids:
        el.open_controller(ctrl[cid]['gclib'], ctrl[cid]['ip_address'])
        printout_str = 'Opened connection to controller: {},'.format(cid)
        if ai_unscaled_vrange is not None:
            for pin in full_ai_pins:
                el.set_controller_ai_vrange(ctrl[cid]['gclib'], pin, vrange=ai_unscaled_vrange)
            printout_str += ' AI unscaled range: {} V,'.format(ai_unscaled_vrange)
        if ao_unscaled_vrange is not None:
            for pin in full_ao_pins:
                el.set_controller_ao_vrange(ctrl[cid]['gclib'], pin, vrange=ao_unscaled_vrange)
            printout_str += ' AO unscaled range: {} V'.format(ao_unscaled_vrange)
        if ground_ao:
            for pin in full_ao_pins:
                set_amp_ao(cid, pin, 0.0)
        if printout:
            print(printout_str)
    if ground_ao and printout:
        time.sleep(0.5)
        mean_ai_read_volt = np.mean([[read_amp_ai(cid, pin) for pin in full_ai_pins] for cid in cids])
        print('All AO of all controllers grounded. Mean voltage across all AI is: {:.6f} V'.format(mean_ai_read_volt))
    if enable_psu_output:
        pcon.enable_psu_output(psu, printout=printout)
        if printout:
            time.sleep(0.5)
            psu_output_state = pcon.get_psu_output_state(psu, return_output_state=True)
            psu_read_voltage = pcon.read_psu_voltage(psu)
            print(psu_output_state+'. Internal PSU voltmeter reads: {} V.'.format(psu_read_voltage))
    return psu, rm

def close_amp_system(cids, psu_session, resource_manager, reset_psu=True, ground_ao=True, printout=True):
    """
    Closes connection to PSU and controllers.
    cids: list of ints -> list of controller IDs whose controllers will be connected to.
    psu_session: pyvisa.resources.usb.USBInstrument -> the instance of the PSU session that was used to
        connect to the PSU.
    resource_manager: pyvisa.highlevel.ResourceManager -> the instance of the resource manager that was
        used to connect to the PSU.
    reset_psu: bool -> If True (recommended), then the PSU will enter the reset state.
    ground_ao: bool -> If True (recommended), then all AO pins of all cids will be set to 0 V.
    """
    if reset_psu:
        pcon.reset_psu(psu_session, printout=printout)
    if ground_ao:
        [[set_amp_ao(cid, pin, 0.0) for pin in full_ai_pins] for cid in cids]
    if printout:
        time.sleep(0.5)
        mean_ai_read_volt = np.mean([[read_amp_ai(cid, pin) for pin in full_ai_pins] for cid in cids])
        print('All AO of all controllers grounded. Mean voltage across all AI is: {:.6f} V'.format(mean_ai_read_volt))
    for cid in cids:
        el.close_controller(ctrl[cid]['gclib'])
        if printout:
            print('Closed connection to controller: {}'.format(cid))

def set_amp_ao(cid, ao_pin, amp_voltage, slope=None, intercept=None):
    """
    Sets an amplified voltage on an AO pin of a controller. Default is to obtain slope and intercept from
    controllers[cid]['ao_cal_fit'][ao_pin].
    cid: int -> controller ID
    ao_pin: int -> AO pin to set
    amp_voltage: float -> Amplified voltage (volts)
    slope: (optional) float -> slope to convert amplified voltage to AO voltage.
    intercept: (optional) float -> intercept to convert amplified voltage to AO voltage.
    """
    if amp_voltage < psu_min_voltage or amp_voltage > psu_max_voltage:
        raise AmpSystemError('Voltage ({}) is outside the voltage range that the PSU can amplify ({}-{}).'.format(amp_voltage, psu_min_voltage, psu_max_voltage))
    if ao_pin not in full_ao_pins:
        raise AmpSystemError('AO pin {} not in list of avialable pins: {}'.format(ao_pin, full_ao_pins))
    else:
        if slope is None:
            slope = ctrl[cid]['ao_cal_fit'][ao_pin][0]
        if intercept is None:
            intercept = ctrl[cid]['ao_cal_fit'][ao_pin][1]
        ao_voltage = slope*amp_voltage + intercept
        el.set_ao(ctrl[cid]['gclib'], ao_pin, ao_voltage)

def read_amp_ai(cid, ai_pin, slope=None, intercept=None):
    """
    Returns an amplified voltage reading of an AI pin of a controller (in vollts). 
    Default is to obtain slope and intercept from controllers[cid]['ai_cal_fit'][ai_pin].
    cid: int -> controller ID
    ai_pin: int -> AI pin to read
    slope: (optional) float -> slope to convert AI pin voltage to amplified voltage.
    intercept: (optional) float -> intercept to convert AI pin voltage to amplified voltage.
    """
    ai_voltage = el.read_ai(ctrl[cid]['gclib'], ai_pin)
    if slope is None:
        slope = ctrl[cid]['ai_cal_fit'][ai_pin][0]
    if intercept is None:
        intercept = ctrl[cid]['ai_cal_fit'][ai_pin][1]
    amp_voltage = slope*ai_voltage + intercept
    #print('AI{}: slope: {}, int: {}'.format(ai_pin, slope, intercept))
    return amp_voltage

class AmpSystemError(Exception):
    """
    Custom exception related to using the amp system
    """
    pass
