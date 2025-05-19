from . import controllers as ctrl
import numpy as np
import math
import copy
import pickle

def construct_cells(cell_num_array, cell_trans_array=None, controllers=None, ao_pins=None, ai_pins=None):
    """
    Returns a list of cell objects
    cell_num_array: a 2D array whose elements represent the cell numbers in the locations you want to assign. Cell numbers 
    start at 1.NaNs are allowed to in cell_num_array to form rugged or uneven cell arrangements.
    
    cell_trans_array (optional): represents is an array that is the same shape as cell num array that represents how the controllers
    are physically connected to the controllers. For an optic with N cells, the convention is:
    controller 1: ao pins: [0, 1, 2, ..., 7] -> cells: [1, 2, 3, ..., 8]
    controller 2: ao pins: [0, 1, 2, ..., 7] -> cells: [9, 10, 11, ..., 16]
    ...
    controller J: ao pins: [0, 1, 2, ..., 7] -> cells: [N-7, N-6, N-5, ..., N]
    If cell_trans_array is None, then it is assumed there is a one-to-one mapping between the controllers and 
    cell_num_array and thus cell_trans_array=cell_num_array.
    
    controllers (optional): A list of ints that specify which controllers should be used. If controllers=None, controllers are gathered
    from the controllers.controllers starting at 1.
    
    ao_pins (optional): A list of 1D arrays. ao_pins[i] specifies the array of analog output pins to use for controllers[i]. The total
    number of analog output pins must equal the number of cells. If ao_pins is None, pins 0-7 are used for each controller starting
    from controllers[0] until the number of pins is equal to the number of cells.
    
    ai_pins (optional): A list of 1D arrays. ai_pins[i] specifies the array of analog input pins to use for controllers[i]. ai_pins[i][j]
    is the analog input to read the voltage of ao_pins[i][j]. The total number of analog input pins must equal the number of cells. 
    If ai_pins is None, pins 0-7 are used for each controller starting from controllers[0] until the number of pins is equal to the 
    number of cells.
    """
    # check that controllers is formatted based on number of cells
    cell_nos = order_nanarray(cell_num_array)
    N_cells = len(cell_nos)
    controllers = validate_controllers(controllers, N_cells)
    # check that cell_num_array and cell_trans_array are formatted correctly
    cell_trans_array = match_nan_locations(cell_num_array, cell_trans_array)
    # check that ao_pins is formatted correctly
    ao_pins = validate_ao_pins(ao_pins, controllers, N_cells)
    flat_ao_pins = np.concatenate(ao_pins, axis=0)
    # check that ai_pins is formatted correctly
    ai_pins = validate_ai_pins(ai_pins, ao_pins, controllers, N_cells)
    flat_ai_pins = np.concatenate(ai_pins, axis=0)
    """
    iterate through cell_nos
    from the current cell no you get the grid coords and the trans cell no
    from the trans cell no you get the controller, ao, and di pins
    """
    cells = []
    for i in range(N_cells):
        cell = Cell()
        cell.no = int(cell_nos[i])
        cell.coords = np.nonzero(cell_num_array==cell.no) # np.argwhere(cell_num_array==cell.no)[0]
        cell.trans_cell_no = int(cell_trans_array[np.nonzero(cell_num_array==cell.no)][0])
        cell.idx = int(cell.trans_cell_no - 1)
        cell.ao_pin = int(flat_ao_pins[cell.idx])
        cell.ai_pin = int(flat_ai_pins[cell.idx])
        cell.cid = int(find_controller(controllers, cell.idx, ao_pins))
        cell.ifunc = np.array([])
        cell.rms = -1.
        cell.pv = -1.
        cell.pvq = -1.
        cell.pvr = -1.
        cell.maxInd = np.array([])
        cell.boxCoords = np.array([])
        cell.deadCell = False
        cell.short_cell_nos = np.array([])
        cell.short_volts = np.array([])
        cell.voltage = None
        cell.volt_hist = []
        cell.gnd_voltMap = np.array([])
        cell.high_voltMap = np.array([])
        cell.voltMap = np.array([])
        cell.gnd_figMap = np.array([])
        cell.high_figMap = np.array([])
        cell.figMap = np.array([])
        cell.badIF = False
        cell.no_array = cell_num_array
        cell.trans_no_array = cell_trans_array
        cells.append(cell)
    return cells

def save_cells(cells, pickle_file):
    """
    Save a list of cell objects using pickling.
    """
    attribs = []
    for cell in cells:
        cell_copy = copy.deepcopy(cell)
        keys = [key for key in cell_copy.__dict__.keys()]
        attribs.append([keys, [cell_copy.__dict__[key] for key in keys]])
        f = open(pickle_file, 'wb')
        pickle.dump(attribs, f)
        f.close()

def load_cells(pickle_file):
    """
    Load a list of cell objects from a pickle file.
    """
    f = open(pickle_file,'rb')
    attribs = pickle.load(f)
    f.close()
    N_cells = len(attribs)
    N_values = len(attribs[0][0])
    keys = attribs[0][0]
    blank_cells = [Cell() for i in range(N_cells)]
    for i, blank_cell in enumerate(blank_cells):
        for j in range(N_values):
            blank_cell.__dict__[keys[j]] = attribs[i][1][j]
    return blank_cells

def print_cells_info(cells):
    """
    Print all the attributes of a list of cell objects.
    """
    dash_num = 15
    for cell in cells:
        print('='*dash_num+'CELL #: {}'.format(cell.no)+'='*dash_num)
        print('-'*dash_num+'Location'+'-'*dash_num)
        print('Coords: {}, Index: {}'.format(cell.coords, cell.idx))
        print('Translated cell #: {}'.format(cell.trans_cell_no))
        print('maxInd:\n {}\nDead Cell: {}, Shorted Cell #: {}'.format(cell.maxInd, cell.deadCell, cell.short_cell_nos))
        print('-'*dash_num+'Control'+'-'*dash_num)
        print('CTRL: {}, AO: {}, AI: {}'.format(cell.cid, cell.ao_pin, cell.ai_pin))
        print('-'*dash_num+'IF'+'-'*dash_num)
        print('IF Shape: {}, RMS: {:.2f} um, PV: {:.2f} um'.format(cell.ifunc.shape, cell.rms, cell.pv))
        print('='*((dash_num*2)+11)+'\n')
        
class Cell:
    """
    Cell object that has attributes that correspond to a single actuator on an adjustable optic
    """

    def __init__(self):
        self.idx = None # int, the index of the cell in the list of cells generated
        self.no = None # int, cell number
        self.coords = None # tuple, the coordinates to access self.no from cell_num_array
        self.trans_cell_no = None # int, the translated cell number
        self.idx = None # int, the index of the translated cell number, self.idx = self.trans_cell_no - 1
        self.ao_pin = None # int, the analog output pin assigned to the cell
        self.ai_pin = None # int, the analog input pin assinged to the cell
        self.cid = None # int, the controller ID number that corresponds to the cell
        self.ifunc = None # np array, the 2D array representing the influence function of the cell,
                          # self.ifunc = self.high_figMap - self.gnd_figMap
        self.rms = None # float, the root mean square of self.ifunc
        self.pv = None # float, the peak-to-valley of self.ifunc
        self.maxInd = None # np array, [y_coord, x_coord], the coordinates to the centroid maximum of self.ifunc
        self.boxCoords = None # np array, the coordinates of the corners that form a box around the centroid of self.ifunc
        self.deadCell = False # bool, If the cell cannot maintain a voltage when assigned, badCell is True
        self.short_cell_nos = None # np array, the array of cell numbers that the cell is electrically shorted to
        self.short_volts = None # np array, the voltages read for short_cell_nos
        self.voltage = None # float, the current voltage of the cell
        self.volt_hist = None # list, list of floats that traces the history of the cell's voltage
        self.gnd_voltMap = None # np array, 2D array that represents the voltages of all cells
                                # when all cells are grounded prior to the cell being energized
        self.high_voltMap = None # np array, 2D array that represents the voltages of all cells
                                # when the cell was energized
        self.voltMap = None # np array, 2D array that rperesents the voltages of all cells for an arbitrary voltage measurement
        self.gnd_figMap = None # np array, 2D array that represents the interferometer measurement of the cell 
                                # when all cells are grounded prior to the cell being energized
        self.high_figMap = None # np array, 2D array that represents the interferometer measurement of the cell 
                                # when cell was energized
        self.figMap = None # np array, 2D array that represents an arbitrary interferometer measurement
        self.dx = None # float, the pixel spacing in mm/pix of self.gnd_figMap, self.high_figMap, and self.figMap
        self.badIF = None # bool, determines whether the IF that was measured is observable
        self.no_array = None # np array, the cell_num_array that was passed when construct_cells() was used to create the cells
        self.trans_no_array = None # np array, the cell_trans_array that was passed when construct_cells() was used to create the cells

    def add_if(self, ifunc):
        self.ifunc = ifunc
        self.rms = alsis.rms(ifunc)
        self.pv = alsis.ptov(ifunc)

    def add_maxInd(self, maxInd, short_cell_no):
        self.maxInd = maxInd
        if short_cell_no is not None:
            self.short_cell_no = short_cell_no

def order_nanarray(a):
    """
    Returns the sorted non-nan elements of an array as a 1D array
    """
    return np.sort(a[~np.isnan(a)].flatten())

def validate_controllers(controllers, N_cells):
    """
    Checks that the choice of controllers is avialable to support the number of cells
    """
    if controllers is None:
        N_ctrl_req = int(math.ceil(N_cells/ctrl.N_ao))
        print('Number of cells:', N_cells)
        print('controllers required: {}'.format(N_ctrl_req))
        if N_ctrl_req > ctrl.N_controllers:
            raise CellError('Not enough outputs available ({}) to support {} cells.'.format(ctrl.N_controllers*ctrl.N_ao, 
                                                                                            N_cells))
        else:
            return [i for i in range(1, N_ctrl_req+1)]
    else:
        controller_IDs = list(ctrl.controllers.keys())
        for ID in controllers:
            if ID not in controller_IDs:
                raise CellError('Controller: {} not in list of known controllers: {}'.format(ID, controller_IDs))
        return controllers

def match_nan_locations(a1, a2):
    """
    Checks that cell_num_array and cell_trans_array have the same shape and have NaNs (if any) in the same lcoations
    """
    if a2 is None:
        a2 = np.copy(a1)
    if a1.shape != a2.shape:
        raise CellError('Unequal shapes. cell_num_array has shape {} but cell_trans_array is of shape {}.'.format(a1.shape, a2.shape))
    if np.isnan(a1).any() and not np.isnan(a2).any():
        raise CellError('cell_num_array contains NaNs but cell_trans_array does not.')
    elif np.isnan(a2).any() and not np.isnan(a1).any():
        raise CellError('cell_trans_array contains NaNs but cell_num_array does not.')
    elif np.isnan(a1).any() and np.isnan(a2).any():
        if np.all(np.isnan(a1) == np.isnan(a2)):
            return a2
        else:
            raise CellError('NaN locations in cell_num_array and cell_trans_array do not match')
    else:
        return a2

def validate_ao_pins(ao_pins, controllers, N_cells):
    """
    Check that the choice of ao_pins can support the number of cells.
    """
    if ao_pins is None:
        ao_pins = [ctrl.ao_pins] * (N_cells//len(ctrl.ao_pins))
        if N_cells % len(ctrl.ao_pins) != 0:
            ao_pins.append(ctrl.ao_pins[:N_cells%len(ctrl.ao_pins)])    
    else:
        N_ao_pins = sum(len(pins) for pins in ao_pins)
        if N_ao_pins != N_cells:
            raise CellError('Number of AO pins ({}) is unequal to the number of cells () provided'.format(N_ao_pins, N_cells))
        for i, controller_pins in enumerate(ao_pins):
            for pin in controller_pins:
                if pin not in ctrl.ao_pins:
                    raise CellError('Controller {}, AO pin {} not in list of supported AO pins: {}'.format(controllers[i], pin, 
                                                                                                           ctrl.ao_pins))
                else:
                    continue
    return ao_pins

def validate_ai_pins(ai_pins, ao_pins, controllers, N_cells):
    """
    Check that the choice of ai_pins can support the number of cells and choice of ao_pins
    """
    if ai_pins is None:
        ai_pins = []
        for controller_pins in ao_pins:
            ai_pins.append(ctrl.ai_pins[:len(controller_pins)])
    else:
        for i, controller_ai_pins in enumerate(ai_pins):
            if len(controller_ai_pins) != len(ao_pins[i]):
                raise CellError('Controller {}, number of AI pins ({}) does not match number of AO pins ({})'.format(controllers[i], 
                                                                                                                 controller_ai_pins, 
                                                                                                                 ao_pins[i]))
            for pin in controller_ai_pins:
                if pin not in ctrl.ai_pins:
                    raise CellError('Controller {}, AI pin {} not in list of supported AI pins: {}'.format(controllers[i], pin, 
                                                                                                           ctrl.ai_pins))
                else:
                    continue
    return ai_pins

def find_controller(controllers, idx, ao_pins):
    """
    Return the corresponding controller based on the current index
    """
    pins_taken = len(ao_pins[0])
    for controller_idx, controller_pins in enumerate(ao_pins):
        if idx < pins_taken:
            return controllers[controller_idx]
        else:
            pins_taken += len(ao_pins[controller_idx+1])

class CellError(Exception):
    """
    Custom exception related to constructing cells.
    """
    pass
