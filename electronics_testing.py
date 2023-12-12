import axroGalil.electronics as el
import numpy as np
import copy
import time

def cell_vpulse_test(cells_input, pulseRange=[0., 10.], ao_range=[0., 10.], ai_range=[0., 10.], 
                     N_steps=150, wr_delay=0.5):
    """
    Tests setting and reading voltages for a list of cells over a reflected pulseRange over a number measurements set by N_steps
    """
    cells = copy.deepcopy(cells_input)
    rising_volts = np.linspace(pulseRange[0], pulseRange[1], int(N_steps/2))
    test_volts = np.concatenate((rising_volts, np.flip(rising_volts)))
    el.open_cell_controllers(cells, setVrange=False, printout=True)
    el.set_cells_ao_vrange(cells, vrange=ao_range, printout=False)
    el.set_cells_ai_vrange(cells, vrange=ai_range, printout=False)
    start_time = time.time()
    el.setAllCells(cells, 0., setAllPins=True, readback=False, delay=wr_delay, printout=False)
    for cell in cells:
        cell.volt_hist = []
    for i in range(len(test_volts)):
        el.setAllCells(cells, test_volts[i], setAllPins=False, readback=True, delay=wr_delay)
        print('Tested {:.3f} V'.format(test_volts[i]))
    el.setAllCells(cells, 0., setAllPins=True, readback=False, delay=wr_delay, printout=False)
    end_time = time.time() - start_time
    el.close_cell_controllers(cells, printout=True)
    print('-'*20+'\nTest complete. Time elapsed: {:.2f} s = {:.2f} min = {:.2f} hr'
          .format(end_time, end_time/60., end_time/3600))
    return cells