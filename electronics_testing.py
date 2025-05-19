import axroGalil.electronics as el
from axroGalil.webIF import WebIF
from axroGalil.construct_connections import order_nanarray
import utilities.metrology as met
from utilities.imaging.man import newGridSize
import numpy as np
import copy
import time
import os

def cells_voltage_test(cells_input, tvolts=[0., 1., 5., 10.], ao_range=[0., 10.], ai_range=[0., 10.], 
                       wr_delay=0.125, thresh=0.5, printout=False):
        """
        Tests setting all cells to each voltage in tvolts. If a cell's voltage is < tvolt-thresh or > tvolt+thresh,
        cell.deadCell will be set to True for each cell.
        """
        cells = copy.deepcopy(cells_input)
        el.open_cell_controllers(cells, setVrange=False, printout=True) # initialize controllers, ground cells
        el.set_cells_ao_vrange(cells, vrange=ao_range, printout=False)
        el.set_cells_ai_vrange(cells, vrange=ai_range, printout=False)
        start_time = time.time()
        el.setAllCells(cells, 0., setAllPins=True, readback=False, delay=wr_delay, printout=False)
        for cell in cells: # clear cells' voltage history
            cell.volt_hist = []
        for i in range(len(tvolts)): # set and read all cells at each test voltage
            el.setAllCells(cells, tvolts[i], setAllPins=False, readback=False, delay=wr_delay)
            time.sleep(wr_delay)
            for cell in cells: # read back the voltage that was set
                read_voltage = el.readCell(cell.no, cells, delay=None)
                cell.volt_hist.append(read_voltage)
                if (read_voltage<tvolts[i]-thresh) or (read_voltage>tvolts[i]+thresh): # check if cell is bad
                    cell.deadCell = True
                if printout:
                    print('Cell {}: {} V'.format(cell.no, read_voltage))
            print('Tested {} V'.format(tvolts[i]))
        el.setAllCells(cells, 0., setAllPins=True, readback=False, delay=wr_delay, printout=False)
        end_time = time.time() - start_time
        el.close_cell_controllers(cells, printout=True)
        print('-'*20+'\nTest complete. Time elapsed: {:.2f} s = {:.2f} min = {:.2f}hr'
              .format(end_time, end_time/60, end_time/3600))
        return cells

def cells_short_test(cells_input, tvolt=10., ao_range=[0., 10.], ai_range=[0., 10.], wr_delay=0.125, thresh=0.5):
    """
    Tests for which cells are shorted. All cells are grounded. Then a single cell is set to tvolt. If any other cells
    are measured to be greater than thresh, then those cells are marked as being shorted to the cell that was enegized
    via cell.short_cell_nos.
    """
    cells = copy.deepcopy(cells_input)
    el.open_cell_controllers(cells, setVrange=False, printout=True) # initialize controllers, ground cells
    el.set_cells_ao_vrange(cells, vrange=ao_range, printout=False)
    el.set_cells_ai_vrange(cells, vrange=ai_range, printout=False)
    start_time = time.time()
    el.setAllCells(cells, 0., setAllPins=True, readback=False, delay=wr_delay, printout=False)
    for cell in cells: # clear voltage history of all cells
        cell.volt_hist = []
        cell.short_cell_nos = np.array([])
        cell.short_volts = np.array([])
    for i, cell in enumerate(cells): # set cell to tvolt
        el.setCell(cell.no, cells, tvolt, readback=True, delay=wr_delay, printout=False)
        for j, check_cell in enumerate(cells): # read the voltage of another cell
            if i == j: # skip measuring the cell that was just energized
                continue
            read_voltage = el.readCell(check_cell.no, cells, delay=None)
            check_cell.voltage = read_voltage
            check_cell.volt_hist.append(read_voltage)
            if read_voltage > thresh: # check_cell that was read is shorted
                if cell.no not in cell.short_cell_nos: # add the energized cell's own number to its short list
                    cell.short_cell_nos = np.append(cell.short_cell_nos, cell.no)
                    cell.short_volts = np.append(cell.short_volts, cell.voltage)
                if check_cell.no not in check_cell.short_cell_nos: # add the checked cell's own number to its short list
                    check_cell.short_cell_nos = np.append(check_cell.short_cell_nos, check_cell.no)
                    check_cell.short_volts = np.append(check_cell.short_volts, read_voltage)
                cell.short_cell_nos = np.append(cell.short_cell_nos, check_cell.no) # update energized cell
                cell.short_volts = np.append(cell.short_volts, read_voltage)
                check_cell.short_cell_nos = np.append(check_cell.short_cell_nos, cell.no) # update checked cell
                check_cell.short_volts = np.append(check_cell.short_volts, cell.voltage)
        el.setCell(cell.no, cells, 0., readback=False, delay=wr_delay, printout=False)
    for cell in cells: # order each short_cell_nos and short_volts array
        sort_indices = np.argsort(cell.short_cell_nos)
        cell.short_volts = cell.short_volts[sort_indices].astype(int)
        cell.short_cell_nos = np.sort(cell.short_cell_nos).astype(int)
    for i, cell in enumerate(cells): # check that all short_cell_nos arrays agree appropriately
        for j, check_cell in enumerate(cells):
            if i == j:
                continue
            if cell.short_cell_nos.size>0 and np.all(np.in1d(cell.short_cell_nos, check_cell.short_cell_nos)):
                cell.short_cell_nos = np.copy(check_cell.short_cell_nos)
                cell.short_volts = np.copy(check_cell.short_volts)
    end_time = time.time() - start_time
    el.close_cell_controllers(cells, printout=True)
    print('-'*20+'\nTest complete. Time elapsed: {:.2f} s = {:.2f} min = {:.2f}hr'
            .format(end_time, end_time/60, end_time/3600))
    return cells

#def check_cells_for_shorts(cells_input, thresh):
#    """
#    Checks cells' volt_hist for shorts between cells. Updates cell.short_cell_nos for each cell. thresh is the voltage
#    above 0 V to determine if a cell is voltage.
#    """
#    cells = copy.deepcopy(cells_input)
#    for cell in cells:
#        cell.short_cell_nos = np.array([])
#    for i in range(len(cells)): # indexes the cell that was energized
#        print('-'*20)
#        print('Cell {} was energized'.format(i+1))
#        for j, cell in enumerate(cells): # indexes the cell that is being checked for shorts
#            if i == j: # skip measuring the cell that was energized
#                continue
#            else:
#                voltage = cell.volt_hist[i]
#                print('Cell {} was at {:.2f} V'.format(cell.no, voltage))
#                if voltage>thresh:
#                    print('Cell {} was above the threshold'.format(cell.no))
#                    #if cell.short_cell_nos is None:
#                    #    cell.short_cell_nos = np.array([])
#                    print("I am adding cell {} to cell {}'s short_cell_nos".format(cells[i].no, cell.no))
#                    cell.short_cell_nos = np.sort(np.append(cell.short_cell_nos, np.array(cells[i].no))).astype(int)
#                    cells[i].short_cell_nos = np.sort(np.append(cells[i].short_cell_nos, np.array(cell.no))).astype(int)
#    for cell in cells:
#        if cell.short_cell_nos.size == 0:
#            continue
#        else:
#            for short_cell_no in cell.short_cell_nos:
#                shorted_cell = el.filter_cells(cells, 'no', short_cell_no)[0]
#                cell.short_cell_nos = np.sort(np.unique(np.append(cell.short_cell_nos, 
#                                                                  shorted_cell.short_cell_nos))).astype(int)

#    return cells

def cell_vpulse_test(cells_input, pulseRange=[0., 10.], ao_range=[0., 10.], ai_range=[0., 10.], 
                     N_steps=150, wr_delay=0.125):
    """
    Tests setting and reading voltages for a list of cells over a reflected pulseRange over a number measurements set by N_steps
    """
    cells = copy.deepcopy(cells_input)
    rising_volts = np.linspace(pulseRange[0], pulseRange[1], int(N_steps/2)) # construct test volt array
    test_volts = np.concatenate((rising_volts, np.flip(rising_volts)))
    el.open_cell_controllers(cells, setVrange=False, printout=True) # initialize controllers, ground cells
    el.set_cells_ao_vrange(cells, vrange=ao_range, printout=False)
    el.set_cells_ai_vrange(cells, vrange=ai_range, printout=False)
    start_time = time.time()
    el.setAllCells(cells, 0., setAllPins=True, readback=False, delay=wr_delay, printout=False)
    for cell in cells: # clear voltage history
        cell.volt_hist = []
    for i in range(len(test_volts)): # set and read all cells at each test voltage
        el.setAllCells(cells, test_volts[i], setAllPins=False, readback=False, delay=wr_delay)
        time.sleep(wr_delay)
        for cell in cells:
            cell.volt_hist.append(el.readCell(cell.no, cells, delay=None))
        print('Tested {:.3f} V'.format(test_volts[i]))
    el.setAllCells(cells, 0., setAllPins=True, readback=False, delay=wr_delay, printout=False)
    end_time = time.time() - start_time
    el.close_cell_controllers(cells, printout=True)
    print('-'*20+'\nTest complete. Time elapsed: {:.2f} s = {:.2f} min = {:.2f} hr'
          .format(end_time, end_time/60., end_time/3600))
    return cells, test_volts

def get_vpulse_test_voltMaps(cells):
    """
    Get the 3D array of voltMaps from a list of cells that were returned from cell_vpulse_test().
    This is used so that the voltMaps can be passed to anime.displayVoltMaps()
    """
    cell_num_array = cells[0].no_array
    N_steps = len(cells[0].volt_hist)
    voltMaps = np.full((N_steps, cell_num_array.shape[0], cell_num_array.shape[1]), 
                       np.nan)
    for i in range(voltMaps.shape[0]):
        for cell in cells:
            voltMaps[i][cell.coords] = cell.volt_hist[i]
    return voltMaps

def measureIFs(input_cells, IFmeas_dir, IF_voltage=10, selectCellNos=None, frames=32, addIFs_to_cells=False, metGeo=met.readFlat4D_h5,
               figMapNewShape=None, cells_fname=None, ao_range=[0., 10.], ai_range=[0., 10.], alternate_file_handle=None, alternate_IF_dir=None, 
               alternate_meas_num=0, volt_label_ground_measurements=False, overwrite_existing_files=True):
    """
    Measures a set of influence functions for a list of cells. The grounded figMap and high figMap associated with each IF
    must be saved to IFmeas_dir before they can be used to construct the IFs.
    input_cells: list, list of the FULL cell objects returned from construct_cells()
    IFmeas_dir: str, directory to save raw inteferometer measurements
    IF_voltage: float, the voltage each IF is measured at
    selectCellNos (optional): list, list of integers that specifies which cells should be measured. If None, then all cells in input_cells are measured
    frames: int, the number of frames the inteferometer gathers to construct a single measurement
    addIFs_to_cells: bool, if True, the IFs will be automatically constructed by loading the inteferometer data from IFmeas_dir.
        The ground figmaps, high figmaps, and influence functions will added as attributes to each corresponding cell.
    metGeo: If addIFs_to_cells is True, metGeo specifies the function used to load the H5 files to construct the IFs.
    figMapNewShape: If addIFs_to_cells is True, the ground figmaps, high figmaps, and influence functions will be reshaped to figMapNewShape
        before being added as attributes to the cells
    cells_fname (optional): str, if specified, saves the list of cells with updated voltMaps (and optionally the figmaps if
        addIfs_to_cells is True)
    """
    cells = copy.deepcopy(input_cells)
    if selectCellNos is None:
        selectCellNos = [cell.no for cell in cells]
    if alternate_IF_dir is None:
        alternate_IF_dir = IFmeas_dir
    w = WebIF() # create instance of 4D web service
    w.setTimeout(150000) # set the timeout to 2.5 minutes

    el.open_cell_controllers(cells, setVrange=False, printout=True) # initialize controllers
    el.set_cells_ao_vrange(cells, vrange=ao_range, printout=False)
    el.set_cells_ai_vrange(cells, vrange=ai_range, printout=False)
    el.setAllCells(cells, 0., setAllPins=True, readback=True, printout=False)
    meas_start_time = time.time()
    for i, cell in enumerate(cells):
        if cell.no not in selectCellNos:
            print('cell {} is not in selectCellNos: {}. Skipping...'.format(cell.no, selectCellNos))
            continue
        print('-'*12+'CELL # {}'.format(cell.no)+'-'*12)
        cell.gnd_voltMap = get_voltMap(cells) # get the voltage map with all cells grounded
        w.averageMeasure(frames) # take an average measurement with all cells grounded

        if not volt_label_ground_measurements:
            ground_meas_file_name = 'cell_{}_grounded'.format(cell.no)
        else:
            ground_meas_file_name = 'cell_{}_{}V_grounded'.format(cell.no, str(IF_voltage))
            print('overwrite existing files:', overwrite_existing_files)
            print('does {} already exist?\n {}'.format(IFmeas_dir+ground_meas_file_name, os.path.isfile(IFmeas_dir+ground_meas_file_name+'.h5')))
        if os.path.isfile(IFmeas_dir+ground_meas_file_name+'.h5') and not overwrite_existing_files:
            dupl_meas_num = 1
            while os.path.isfile(IFmeas_dir+ground_meas_file_name+'.h5'):
                if dupl_meas_num == 1:
                    ground_meas_file_name = ground_meas_file_name + '_{}'.format(dupl_meas_num)
                else:
                    char_list = list(ground_meas_file_name)
                    char_list[-1] = str(dupl_meas_num)
                    ground_meas_file_name = ''.join(char_list)
                dupl_meas_num += 1

        w.saveMeasurement(IFmeas_dir+ground_meas_file_name) # save the grounded measurement
        if alternate_file_handle:
            alternate_meas_num_str = get_alternate_meas_num_str(alternate_meas_num)
            w.saveMeasurement(alternate_IF_dir+alternate_file_handle+'_4DRAW_{}_{}'.format(alternate_meas_num_str, '000'))
            alternate_meas_num += 1
        el.setCell(cell.no, cells, IF_voltage, readback=True, printout=True) # energize the cell
        cell.high_voltMap = get_voltMap(cells) # get the voltage map with the cell energized
        w.averageMeasure(frames) # take an average measurement with the cell energized
        if alternate_file_handle:
            alternate_meas_num_str = get_alternate_meas_num_str(alternate_meas_num)
            w.saveMeasurement(alternate_IF_dir+alternate_file_handle+'_4DRAW_{}_{}'.format(alternate_meas_num_str, '001'))
            alternate_meas_num += 1

        high_meas_file_name = 'cell_{}_{}V'.format(cell.no, str(IF_voltage))
        if os.path.isfile(IFmeas_dir+high_meas_file_name+'.h5') and not overwrite_existing_files:
            dupl_meas_num = 1
            while os.path.isfile(IFmeas_dir+high_meas_file_name+'.h5'):
                if dupl_meas_num == 1:
                    high_meas_file_name = high_meas_file_name + '_{}'.format(dupl_meas_num)
                else:
                    char_list = list(high_meas_file_name)
                    char_list[-1] = str(dupl_meas_num)
                    high_meas_file_name = ''.join(char_list)
                dupl_meas_num += 1

        w.saveMeasurement(IFmeas_dir+high_meas_file_name) # save the energized measurement
        el.setAllCells(cells, 0., setAllPins=False, readback=True, printout=False) # ground all cells
        print('Cells grounded.')
        time.sleep(1)
    el.close_cell_controllers(cells, printout=True) # close controllers
    meas_dt = time.time() - meas_start_time
    print('-'*30) # print the time elapsed
    print('Measurements complete. Cells are now grounded.')
    print('Total time to measure {} cells: {:.2f} min = {:.2f} hrs.\nIF measure rate is {:.2f} cells/hr.'
          .format(len(selectCellNos), meas_dt/60, meas_dt/3600, len(selectCellNos)/(meas_dt/3600)))
    if addIFs_to_cells: # construct the gnd figmaps, high figmaps, and IFs and add them to the cell objects
        load_start_time = time.time()
        load_cells = [cell for cell in cells if cell.no in selectCellNos]
        cells = construct_IFs_from_meas_dir(load_cells, IFmeas_dir, IF_voltage, metGeo=metGeo, newShape=figMapNewShape, 
                                            alternate_file_handle=alternate_file_handle, volt_label_ground_measurements=volt_label_ground_measurements,
                                            printout=True)
        load_dt = time.time() - load_start_time
        print('Total time to construct {} IFs: {:.2f} min = {:.2f} hrs.\nIF construction rate is {:.2f} IFs/hr.'
          .format(len(load_cells), load_dt/60, load_dt/3600, len(load_cells)/(load_dt/3600)))
    if cells_fname is not None: # save the list of cell objects
        cc.save_cells(cells, cells_fname)
    return cells

def get_alternate_meas_num_str(alternate_meas_num):
    "Returns an integer as a 3 digit string"
    meas_num_str = str(alternate_meas_num)
    while len(meas_num_str) < 3:
        meas_num_str = '0' + meas_num_str
    return meas_num_str

def get_voltMap(cells, readCells=False):
    """
    Takes in a list of cells and generates a single voltMap. The voltMap is an array that is the same shape of cell.no_array whose
    elements are the voltages of cells at the instance the function was called. The voltages are arranged in the same locations
    corresponding to each cell number in cell.no_array.
    If readCells is True, the voltage of each cell is explicitly measured before constructing the voltage map. If False, the 
    voltage is taken from cell.voltage.
    """
    voltMap = np.full(cells[0].no_array.shape, np.nan) # array of all NaNs
    cell_nos = [int(no) for no in order_nanarray(cells[0].no_array)]
    for cell in cells:
        if readCells:
            read_volt = el.readCell(cell.no, cells)
        else:
            read_volt = cell.voltage
        voltMap[cell.coords] = read_volt
    return voltMap

def construct_IFs_from_meas_dir(cells_input, IFmeas_dir, IF_voltage, metGeo=met.readFlat4D_h5, newShape=None, printout=True, 
                                volt_label_ground_measurements=False, alternate_file_handle=None):
    """
    Loads the ground figmaps and high figmaps from IFmeas_dir and uses them to construct an influence function for each cell.
    The ground figmaps, high figmaps, and IFs are then attached to each corresponding cell.
    IF_voltage: float or int, the voltage the IFs were taken at.
    metgeo is the function used to load the h5 files and might need to be changed for different optic geometries.
    newShape (optional): tuple, if specified, the ground figmaps, high figmaps, and IFs are resized to newShape.
    """
    cells = copy.deepcopy(cells_input)
    if printout:
        print('-'*20+'Constructing IFs'+'-'*20)
    for cell in cells:
        if not volt_label_ground_measurements:
            ground_meas_file_name = 'cell_{}_grounded'.format(cell.no)
        else:
            ground_meas_file_name = 'cell_{}_{}V_grounded'.format(cell.no, str(IF_voltage))
        d_gnd, dx_gnd = metGeo(IFmeas_dir+ground_meas_file_name+'.h5') # load the figmaps
        d_high, dx_high = metGeo(IFmeas_dir+'cell_{}_{}V.h5'.format(cell.no, str(IF_voltage)))
        dx = np.mean([dx_gnd, dx_high])
        optic_length = dx * d_gnd.shape[0]
        if newShape is not None: # resize the figmap
            d_gnd = newGridSize(d_gnd, newShape)
            d_high = newGridSize(d_high, newShape)
            dx = optic_length/d_gnd.shape[0]
        cell.gnd_figMap = d_gnd # attach figmaps to cells
        cell.high_figMap = d_high
        cell.dx = dx
        cell.ifunc = d_high - d_gnd
        if printout:
            print('Cell {} IF constructed.'.format(cell.no))
    return cells

def measure_randomVolt_figChanges(input_cells, N_measurements, measDir, minVoltage=0., maxVoltage=10., refVoltage=0., meas_frames=50, start_meas_idx=0, returnMeasurements=True, metGeo=met.readFlat4D_h5, ao_range=[0., 10.], ai_range=[0., 10.]):
    cells = copy.deepcopy(input_cells)
    w = WebIF() # create instance of 4D web service
    w.setTimeout(150000) # set the timeout to 2.5 minutes
    el.open_cell_controllers(cells, setVrange=False, printout=True)
    el.set_cells_ao_vrange(cells, vrange=ao_range, printout=False)
    el.set_cells_ai_vrange(cells, vrange=ai_range, printout=False)
    el.setAllCells(cells, refVoltage, setAllPins=True, readback=True, printout=False)
    gnd_voltMaps = []
    high_voltMaps = []
    start_time = time.time()
    if not os.path.isdir(measDir+'gnd_voltMaps\\'):
        os.makedirs(measDir+'gnd_voltMaps\\')
    if not os.path.isdir(measDir+'high_voltMaps\\'):
        os.makedirs(measDir+'high_voltMaps\\')
    for i in range(N_measurements):
        gnd_voltMap = get_voltMap(cells)
        gnd_voltMaps.append(gnd_voltMap)
        w.averageMeasure(meas_frames)
        w.saveMeasurement(measDir+'meas_ref_{}'.format(start_meas_idx+i))
        np.save(measDir+'gnd_voltMaps\\gnd_voltMap_{}'.format(start_meas_idx+i), gnd_voltMap)
        rand_volts_array = np.random.uniform(low=minVoltage, high=maxVoltage, size=len(cells))
        for cell, rand_volt in zip(cells, rand_volts_array):
            if not cell.badIF:
                el.setCell(cell.no, cells, rand_volt, readback=True, printout=False)
        high_voltMap = get_voltMap(cells)
        high_voltMaps.append(high_voltMap)
        w.averageMeasure(meas_frames)
        w.saveMeasurement(measDir+'meas_energized_{}'.format(start_meas_idx+i))
        np.save(measDir+'high_voltMaps\\high_voltMap_{}'.format(start_meas_idx+i), high_voltMap)
        print('Figure change measurement {} complete. Elapesd time: {:.2f} hrs.'.format(start_meas_idx+i, (time.time()-start_time)/3600.))
        el.setAllCells(cells, refVoltage, setAllPins=False, readback=True, printout=False)
    el.setAllCells(cells, 0., setAllPins=True, readback=True, printout=False)
    print('Cells grounded.')
    time.sleep(1)
    el.close_cell_controllers(cells, printout=True)
    gnd_voltMaps = np.stack(gnd_voltMaps, axis=0)
    high_voltMaps = np.stack(high_voltMaps, axis=0)
    if returnMeasurements:
        figChange_maps, dx = load_random_figChange_measurements(measDir, N_measurements, start_meas_idx, metGeo=metGeo)
    else:
        figChange_maps, dx = None, None
    return figChange_maps, dx, gnd_voltMaps, high_voltMaps

def load_random_figChange_measurements(measDir, N_measurements, start_meas_idx, metGeo=met.readFlat4D_h5, newShape=None, printout=True, ref_meas_file_name='meas_ref', energ_meas_file_name='meas_energized'):
    figChange_maps = []
    dxs = []
    for i in range(N_measurements):
        d_gnd, dx_gnd = metGeo(measDir+ref_meas_file_name+'_{}'.format(start_meas_idx+i)+'.h5') # load the figmaps
        d_high, dx_high = metGeo(measDir+energ_meas_file_name+'_{}'.format(start_meas_idx+i)+'.h5') # load the figmaps
        dx = np.nanmean([dx_gnd, dx_high])
        dxs.append(dx)
        optic_length = dx * d_gnd.shape[0]
        if newShape is not None: # resize the figmap
            d_gnd = newGridSize(d_gnd, newShape)
            d_high = newGridSize(d_high, newShape)
            dx = optic_length/d_gnd.shape[0]
        figChange_map = d_high - d_gnd
        figChange_maps.append(figChange_map)
        if printout:
            print('Figure change map {} constructed.'.format(start_meas_idx+i))
    figChange_maps = np.stack(figChange_maps, axis=0)
    dx = np.nanmean(dxs)
    return figChange_maps, dx

def load_random_voltageMaps(measDir, N_measurements, start_meas_idx):
    gnd_voltMaps, high_voltMaps = [], []
    for i in range(start_meas_idx, N_measurements):
        gnd_voltMaps.append(np.load(measDir+'gnd_voltMaps\\gnd_voltMap_{}.npy'.format(i)))
        high_voltMaps.append(np.load(measDir+'high_voltMaps\\high_voltMap_{}.npy'.format(i)))
    gnd_voltMaps = np.stack(gnd_voltMaps, axis=0)
    high_voltMaps = np.stack(high_voltMaps, axis=0)
    return gnd_voltMaps, high_voltMaps

def hex_ad5535(v_des, v_ref=200, n=14, cal_slope=1.7008880289855073, cal_int=0.01753800000001159):
    """
    Converts a desired voltage to a hexadecimal code for the Analog Devices Eval-AD5535BSDZ.
    v_des: float, desired voltage in the range (0, v_ref).
    v_ref: the reference voltage that is set by the high voltage power supply that is plugged into VPP
    on the eval board.
    n: int, the resolution of the board in bits
    cal_slope: float, the slope for the linear calibration that converts v_des to v_adj for the eval
        board
    cal_int: float, the y intercept for the linear calibration that converts v_des to v_adj for the eval
        board
    """
    v_adj = (v_des-cal_int)/cal_slope
    dac_code = int((2**n - 1) * (v_adj/v_ref))
    hex_code = hex(dac_code)[2:].upper()
    return hex_code

def make_energized_measurement(input_cells, voltages, frames=50, filename='energized_figMap', subtract_gndFigMap=False, gnd_figMap_name=None, figMapNewShape=None, metGeo=met.readFlat4D_h5, ao_range=[0., 10.], ai_range=[0., 10.]):
    """
    Makes an energized measurement of a mirror.
    INPUTS:
    input_cells: list, list of the FULL cell objects returned from construct_cells()
    voltages: 1d array, array of voltages the same length as input cells that represent the voltage to apply to each cell for the measurement.
    frames: int, the number of frames captured by the interferometer to average over for the measurment
    filename: str, the path/filename to save the raw .h5 interferometer measurement to
    figMapNewShape: tuple of ints, optional shape to resize the measurement to after it is read using metGeo
    metGeo: func, the metrology function used to reload and process the raw interferometer measurement into a numpy array
    ao_range: list of floats, the analog output range to set the controllers to
    ai_range: list of floats, the analog input range to set the controllers to
    RETURNS:
    energ_figMap: 2d array, the processed interferometer measurement
    dx: float, the pixel spacing in mm/pix of the measurement
    meas_volts: 1d array, the measured voltages of the cells at the time of the measurement
    """

    cells = copy.deepcopy(input_cells)
    meas_volts = np.zeros(voltages.shape)
    w = WebIF()
    w.setTimeout(150000)
    if gnd_figMap_name is None:
        dirname = os.path.dirname(filename)
        gnd_figMap_name = dirname+'\\ground_figMap'

    el.open_cell_controllers(cells, setVrange=False, printout=True) # initialize controllers
    el.set_cells_ao_vrange(cells, vrange=ao_range, printout=False)
    el.set_cells_ai_vrange(cells, vrange=ai_range, printout=False)
    el.setAllCells(cells, 0., setAllPins=True, readback=True, printout=False)
    time.sleep(2)
    if subtract_gndFigMap:
        dirname = os.path.dirname(filename)
        print('Measuring ground figure map...')
        w.averageMeasure(frames)
        w.saveMeasurement(gnd_figMap_name)
        print('...finished')
    for cell, volt in zip(cells, voltages):
        el.setCell(cell.no, cells, volt, readback=True, printout=True) # energize the cell
    time.sleep(2)
    for i, cell in enumerate(cells):
        meas_volts[i] = el.readCell(cell.no, cells, delay=None) # read the voltage of each cell

    print('Measurement in progress...')
    w.averageMeasure(frames)
    w.saveMeasurement(filename)
    print('...finished')
    el.setAllCells(cells, 0., setAllPins=True, readback=True, printout=False)
    print('Cells grounded.')
    el.close_cell_controllers(cells, printout=True)
    
    energ_figMap, dx = metGeo(filename+'.h5')
    if subtract_gndFigMap:
        gnd_figMap, _ = metGeo(gnd_figMap_name+'.h5')
        energ_figMap = energ_figMap - gnd_figMap
    if figMapNewShape:
        optic_length = dx * energ_figMap.shape[0]
        energ_figMap = newGridSize(energ_figMap, figMapNewShape)
        dx = optic_length/energ_figMap.shape[0]
    return energ_figMap, dx, meas_volts