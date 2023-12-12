from axroGalil.controllers import controllers as ctrl
from axroGalil.controllers import ao_min_v, ao_max_v, ao_pins
from axroGalil.controllers import ai_min_v, ai_max_v, ai_pins
from axroGalil.controllers import model
from axroGalil.controllers import vranges as hardware_vranges
import gclib
import numpy as np
import time


min_v_allowed = ao_min_v # minimum voltage that can be outputted by an ao pin
max_v_allowed = ao_max_v # maximum voltage that can be outputted by an ao pin
write_read_delay = 0.5 # seconds, the time it takes after setting a voltage to accurately read it back

################ CELL CONTROL FUNCTIONS ################

def open_cell_controllers(cells, setVrange=True, printout=True):
    """
    Opens the connection for the controllers used by a list of cell objects
    If setVrange is True, the controllers will be set to 0 - 10 V analog output and input
    """
    cids = get_controller_ids(cells)
    for cid in cids:
        open_controller(ctrl[cid]['gclib'], ctrl[cid]['ip_address'])
        if printout:
            print('Opened connection to controller: {}'.format(cid))
    if setVrange:
        set_cells_ao_vrange(cells, printout=printout)
        set_cells_ai_vrange(cells, printout=printout)

def close_cell_controllers(cells, printout=True):
    """
    Closes the connection to the controllers associated with a list of cells
    """
    cids = get_controller_ids(cells)
    for cid in cids:
        close_controller(ctrl[cid]['gclib'])
        if printout:
            print('Closed controller: {}'.format(cid))

def set_cells_ao_vrange(cells, vrange=[0., ao_max_v], printout=True):
    """
    Sets the allowed analog output voltage range for a list of cells' associated controllers
    """
    for cell in cells:
        set_controller_ao_vrange(ctrl[cell.cid]['gclib'], cell.ao_pin, vrange)
        if printout:
            print('Controller {} AO pin {} range set to {}-{} V'.format(cell.cid, cell.ao_pin, 
                                                                        vrange[0], vrange[1]))

def set_cells_ai_vrange(cells, vrange=[0., ai_max_v], printout=True):
    """
    Sets the allowed analog input voltage range for a list of cells' associated controllers
    """
    for cell in cells:
        set_controller_ai_vrange(ctrl[cell.cid]['gclib'], cell.ai_pin, vrange)
        if printout:
            print('Controller {} AI pin {} range set to {}-{} V'.format(cell.cid, cell.ai_pin, 
                                                                        vrange[0], vrange[1]))

def setCell(cell_no, cells, voltage, readback=True, delay=write_read_delay, printout=True):
    """
    Sets the voltage on a cell. If readback is True, readCell() is run on the same cell and cell.voltage
    and cell.volt_hist are updated with the voltage that is read back. If readback is True, delay 
    specifies how long to wait after setting the voltage to read it back.
    """
    cell = filter_cells(cells, 'no', cell_no)[0]
    set_ao(ctrl[cell.cid]['gclib'], cell.ao_pin, voltage)
    if readback:
        time.sleep(delay)
        read_voltage = readCell(cell_no, cells)
        cell.voltage = read_voltage
        cell.volt_hist.append(read_voltage)
        if printout:
            print('Cell {} is now at {} V'.format(cell.no, cell.voltage))

def readCell(cell_no, cells, delay=None):
    """
    Returns the voltage for a cell with number cell_no given a list of cells
    """
    if delay is not None:
        time.sleep(delay)
    cell = filter_cells(cells, 'no', cell_no)[0]
    voltage = read_ai(ctrl[cell.cid]['gclib'], cell.ai_pin)
    return voltage

def setAllCells(cells, voltage, setAllPins=False, readback=True, delay=write_read_delay, printout=False):
    """
    Sets all AO pins associated with a list of cells to a voltage. If setAllPins=True, then all AO pins
    of the controllers associated with cells are set to voltage, not just the ones connected to the cells.
    If readback is True, then cells' voltages are read and updated. If readback is True, then the time 
    between writing and reading the voltage is set by delay.
    """
    cids = get_controller_ids(cells)
    for cid in cids:
        for pin in ao_pins:
            cells_matching_controller = filter_cells(cells, 'cid', cid)
            cell_matching_pin = filter_cells(cells_matching_controller, 'ao_pin', pin)
            if cell_matching_pin:
                cell = cell_matching_pin[0]
                setCell(cell.no, cells, voltage, readback=readback, delay=delay, printout=printout)
            else:
                if setAllPins:
                    set_ao(ctrl[cid]['gclib'], pin, voltage)
                    if printout:
                        print('Controller {} AO pin {} was set to {} V'.format(cid, pin, voltage))
                else:
                    continue

def filter_cells(cells, attribute_name, attribute_value):
    """
    Takes in a list of cells and returns a list of cells whose attribute_name is equal to attribute_value
    """
    filtered_cells = [cell for cell in cells if getattr(cell, attribute_name)==attribute_value]
    return filtered_cells

################ DIRECT CONTROLLER FUNCTIONS ################
    
def open_controller(g, ip_address):
    """
    Takes in an instance of gclib and opens a direct connection to controller with a known IP address
    You can create an instance of gclib by using g = gclib.py()
    """
    g.GOpen('{} --direct'.format(ip_address))

def close_controller(g):
    """
    Takes in an instance of gclib and closes the connection to the controller
    """
    g.GClose()
         
def set_controller_ao_vrange(g, ao_pin, vrange=[ao_min_v, ao_max_v]):
    """
    Takes in an instance of gclib and sets the voltage range on AO pin on that controller
    """
    global min_v_allowed
    global max_v_allowed
    vrange = tuple(sorted(vrange))
    if vrange in hardware_vranges['AO'].keys():
        command_code = hardware_vranges['AO'][vrange]
        if ao_pin in ao_pins:
            g.GCommand('DQ {},{}'.format(ao_pin, command_code))
            min_v_allowed = vrange[0]
            max_v_allowed = vrange[1]
        else:
            raise ControllerError('AO pin {} not in list of available pins: {}'.format(ao_pin, ao_pins))
    else:
        raise ControllerError('AO voltage range {} is outside the hardware limited voltage ranges. Options are {}'.format(vrange, hardware_vranges['AO'].keys()))

def set_controller_ai_vrange(g, ai_pin, vrange=[ai_min_v, ai_max_v]):
    """
    Takes in an instance of gclib and sets the voltage range on AI pin on that controller
    """
    vrange = tuple(sorted(vrange))
    if vrange in hardware_vranges['AI'].keys():
        command_code = hardware_vranges['AI'][vrange]
        if ai_pin in ai_pins:
            g.GCommand('AQ {},{}'.format(ai_pin, command_code))
        else:
            raise ControllerError('AI pin {} not in list of available pins: {}'.format(ai_pin, ai_pins))
    else:
        raise ControllerError('AI voltage range {} is outside the hardware limited voltage ranges. Options are {}'.format(vrange, hardware_vranges['AI'].keys()))
    
def set_ao(g, ao_pin, voltage):
    """
    Takes in an instance of gclib and sets the voltage on an analog output pin
    """
    if voltage<min_v_allowed or voltage>max_v_allowed:
        raise ControllerError('Voltage ({}) is outside voltage range set for controller ({}-{}). You can alter the controller voltage range using set_controller_ao_vrange()'
                              .format(voltage, min_v_allowed, max_v_allowed))
    elif ao_pin not in ao_pins:
        raise ControllerError('AO pin {} not in list of available pins: {}'.format(ao_pin, ao_pins))
    else:
        g.GCommand('AO {}={}'.format(ao_pin, voltage))
    
def read_ai(g, ai_pin):
    """
    Takes in an instance of gclib and reads the voltage on an analog input pin
    """
    if ai_pin not in ai_pins:
        raise ControllerError('AI pin {} not in list of available pins: {}'.format(ai_pin, ai_pins))
    else:
        voltage = float(g.GCommand('MG @AN[{}]'.format(ai_pin)))
        return voltage

def get_controller_ids(cells):
    """
    Returns the list of unique controller IDs used to control a list of cells
    """
    return sorted(list(set([cell.cid for cell in cells])))

class ControllerError(Exception):
    """
    Custom exception related to using the controllers
    """
    pass