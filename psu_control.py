"""
Controls high voltage power supply unit (PSU) through USB and pyvisa API
"""

import pyvisa
from .high_voltage_power_supply import psu_min_voltage as min_allowed_voltage
from .high_voltage_power_supply import psu_max_voltage as max_allowed_voltage
from .high_voltage_power_supply import psu_min_current as min_allowed_current
from .high_voltage_power_supply import psu_max_current as max_allowed_current
from .high_voltage_power_supply import psu_visa_address

def open_psu_connection(visa_address=psu_visa_address, init_psu=True, 
                        set_current=max_allowed_current, printout=True):
    """
    Opens the connection to the PSU. If init_psu is True, the PSU enters the
    reset state. If reset_psu is True and set_current is not None, the PSU current
    is set to this value (in Amps) to keep the PSU operating in CV mode.
    Returns:
    psu: an instance of the PSU connection session
    rm: an instance of the resource manager which tracks all connected instruments
        through the pyvisa API.
    """
    # create a connection (session) to the instrument
    rm = pyvisa.ResourceManager()
    psu = rm.open_resource(visa_address)
    # For Serial and TCP/IP socket connections enable the read Termination Character, 
    # or read's will timeout
    if psu.resource_name.startswith('ASRL') or psu.resource_name.endswith('SOCKET'):
        psu.read_termination = '\n'
    if printout: # identify connected instrument
        psu.write('*IDN?')
        response = psu.read().strip()
        print('Opened connection to: {}'.format(response))
    if init_psu:
        reset_psu(psu, set_current_limit=set_current, printout=printout)
    return psu, rm

def close_psu_connection(psu_session, resource_manager, printout=True):
    psu_session.close()
    resource_manager.close()
    if printout:
        print('Connection to PSU closed.')
        
def reset_psu(psu_session, set_current_limit=max_allowed_current,
              printout=True):
    """
    Resets the power supply to the factory-defined state:
    CAL:STAT = OFF, INIT:CONT = OFF, OUTP = OFF
    [SOUR:]CURR = 0, [SOUR:]CURR:TRIG = 0, [SOUR:]CURR:PROT:STAT = OFF
    [SOUR:]VOLT = 0, [SOUR:]VOLT:LIM = 0, [SOUR:]VOLT:TRIG = 0, [SOUR:]VOLT:PROT = MAXimum
    If set_current_limit is not None, the current is set to this value (in Amps). 
    This is necessary to keep the PSU operating in CV mode when the load is connected.
    """
    psu_session.write('*RST')
    if printout:
        print('PSU RESET')
    if set_current_limit is not None:
        set_psu_current(psu_session, set_current_limit)
        if printout:
            print('PSU current limit set to: {} A'.format(set_current_limit))
        
def set_psu_power_on_state(psu_session, printout=True):
    """
    Sets the PSU power-on state to be the same as the state
    described by reset_psu() the next time the PSU is powered on.
    """
    psu_session.write('OUTPUT:PON:STATE RST')
    if printout:
        psu_session.write('OUTPUT:PON:STATE?')
        response = psu_session.read().strip()
        if response == 'RST':
            print('PSU power-on state: RESET STATE')
        else:
            raise PSU_Error('Cannot determine power-on state')
        
def enable_psu_output(psu_session, printout=True):
    """
    Sets the output of the PSU to be on. The voltage 
    and current are set to the values that were last passed to set_psu_voltage()
    and/or set_psu_current(). Running this after reset_psu() sets the voltage and
    current to 0.
    """
    psu_session.write('OUTPUT 1')
    if printout:
        get_psu_output_state(psu_session)

def disable_psu_output(psu_session, printout=True):
    """
    Sets the output of the PSU to be off. The voltage and 
    current outputted are 0 V and 0 A. For setting the PSU to 0V, 
    this method is preferred since the voltage cannot be programmed to
    lower than about 5% above the UVL setting.
    """
    psu_session.write('OUTPUT 0')
    if printout:
        get_psu_output_state(psu_session)

def get_psu_output_state(psu_session, return_output_state=False):
    """
    Prints or returns the output state of the PSU.
    """
    psu_session.write('OUTPUT?')
    output_state = psu_session.read().strip()
    if output_state == '1':
        psu_state = 'PSU output: ON'
    elif output_state == '0':
        psu_state = 'PSU output: OFF'
    else:
        raise PSU_Error('Cannot determine PSU output state.')
    if not return_output_state:
        print(psu_state)
    else:
        return psu_state
        
def set_psu_voltage(psu_session, voltage, readback=False):
    """
    Sets the output voltage of the PSU. voltage is in volts. 
    enable_psu_output() must be run prior to setting the voltage
    for the output to take effect!
    """
    if voltage < min_allowed_voltage or voltage > max_allowed_voltage:
        raise PSU_Error('Error: {} V is outside the range of allowed voltages ({} V).'.format(voltage, [min_allowed_voltage, max_allowed_voltage]))
    else:
        psu_session.write('VOLTAGE {}'.format(float(voltage)))
    if readback:
        read_voltage = read_psu_voltage(psu_session)
        print('PSU voltage output is now at: {} V'.format(read_voltage))
        
def set_psu_current(psu_session, current, readback=False):
    """
    Sets the ouput current (or rather the current limit to keep the PSU in CV mode)
    for the PSU. current in Amps.
    """
    if current < min_allowed_current or current > max_allowed_current:
        raise PSU_Error('Error: {} A is outside the range of allowed current ({} A).'.format([min_allowed_current, max_allowed_current]))
    else:
        psu_session.write('CURRENT {}'.format(float(current)))
    if readback:
        read_current = read_psu_current(psu_session)
        print('PSU current is now at: {} A'.format(read_current))
        
def read_psu_voltage(psu_session):
    """
    Reads the voltage setting of the PSU (in volts).
    The minimum voltage readable is 0.15 V
    """
    psu_session.write('MEASURE:VOLTAGE?')
    response = psu_session.read()
    voltage = float(response.strip())
    return voltage

def read_psu_current(psu_session):
    """
    Reads the current setting of the PSU (in Amps).
    The minimum current readable is 0.001 A.
    """
    psu_session.write('MEASURE:CURRENT?')
    response = psu_session.read()
    current = float(response.strip())
    return current

class PSU_Error(Exception):
    """
    Custom exception related to using the PSU
    """
    pass