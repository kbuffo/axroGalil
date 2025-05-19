from axroGalil.electronics import filter_cells
import axroGalil.plotting.plot_utility_functions as puf
import math
import numpy as np
import matplotlib.pyplot as plt
import copy

def plot_cell_voltageTests(cells, tvolts, date=None, tdelays=None,
                            title='Cell Test Voltages Response',
                            ylim=None, figsize=None, fontsize=12, 
                            cells_per_subplot=16,
                            legendloc='upper center', legendcoords=(0.5, 1), 
                            N_legendcols=1, title_vertGap=1, 
                            pinNumAxesLabels=False):
    """
    Plots the results of running electronics_testing.cells_voltage_test()
    """
    N_subplots = math.ceil(len(cells)/cells_per_subplot)
    if date is not None:
        if '_' in date: 
            date = date.replace('_', '')
    if figsize is None: 
        figsize = (7, 7)
    fig, axs = plt.subplots(N_subplots, 1, figsize=figsize)
    keyword, keyunit = 'voltage', 'V'
    if tdelays is not None: 
        keyword, keyunit = 'delay', 's'
    for i in range(N_subplots):
        ax = np.ravel(axs)[i]
        subplot_cells = [cell for cell in cells[cells_per_subplot*i:cells_per_subplot*(i+1)]]
        cell_nos = [cell.no for cell in subplot_cells]
        ax.set_ylabel('Voltage (V)', fontsize=fontsize)
        for j in range(len(tvolts)):
            volts = [cell.volt_hist[j] for cell in subplot_cells]
            legend_val = tvolts[j]
            if tdelays is not None: 
                legend_val = tdelays[j]
            mean_volt = np.mean([cell.volt_hist[j] for cell in cells])
            ax.plot(cell_nos, volts, marker='.',
                    label='Applied {}: {:.2f} {}, Mean voltage: {:.4f} V'\
                    .format(keyword, legend_val, keyunit, mean_volt))
        if i == 0:
            ax.legend(loc=legendloc, bbox_to_anchor=legendcoords, ncols=N_legendcols,
            fancybox=True, shadow=True, fontsize=fontsize)
        ax.set_xticks(cell_nos)
        if pinNumAxesLabels:
            ticklabels = ['C{},AO{}'.format(cell.cid, cell.ao_pin) for cell in subplot_cells]
            ax.set_xticklabels(ticklabels, fontsize=fontsize-2, rotation=45)
            ax.set_xlabel('Controller Number (C), Analog Output Channel (AO)')
        else:
            ax.set_xticklabels(cell_nos, fontsize=fontsize, rotation=45)
            ax.set_xlabel('Cell #', fontsize=fontsize)
        ax.grid(axis='x')
        all_volts = []
        for cell in cells:
            all_volts += cell.volt_hist
        if ylim:#is None: ylim = [np.min(all_volts)-0.5, np.max(all_volts)+0.5]
            ax.set_ylim(ylim)
        ax.set_xlim([np.min(cell_nos)-1, np.max(cell_nos)+1])
    if date is None:
        combined_title = title
    else:
        combined_title = title + '\nDate: ' + date
    fig.suptitle(combined_title, fontsize=fontsize+2)
    fig.tight_layout(rect=[0, 0, 1, title_vertGap])
    return fig

def plot_cell_shortTest(cells_input, tvolt, shortThresh, figsize=None, ax_fontsize=12, title_fontsize=14, 
                        title='Shorted Cells', pinNumAxesLabels=False):
    """
    Plots the results of running electronics_testing.cells_short_test().
    """
    cells = copy.deepcopy(cells_input)
    short_groups = puf.get_cells_short_groups(cells)
    short_groups_voltages = puf.get_cell_short_groups_voltages(cells)
    if short_groups is None:
        print("Didn't find any short groups.")
        return None
    N_subplots = len(short_groups)
    if figsize is None:
        figsize = (7, 8)
    fig, axs = plt.subplots(N_subplots, 1, figsize=figsize)
    for i in range(N_subplots):
        if N_subplots > 1:
            ax = axs.flatten()[i]
        else:
            ax = axs
        short_group = short_groups[i]
        short_group_voltages = short_groups_voltages[i]
        for j in range(len(short_group)):
            ax.plot(j, short_group_voltages[j], color='black', linestyle='dashed', marker='.', markersize=8)
        if pinNumAxesLabels:
            short_group_cells = [filter_cells(cells, 'no', cell_no)[0] for cell_no in short_group]
            xticklabels = ['C{},AO{}'.format(cell.cid, cell.ao_pin) for cell in short_group_cells]
            xlabel = 'Controller Number (C), Analog Output Channel (AO)'
        else:
            xticklabels = short_group
            xlabel = 'Cell #'
        xvals = [i for i in range(len(short_group))]
        ax.set_xticks(xvals)
        ax.set_xticklabels(xticklabels, fontsize=ax_fontsize-2)
        ax.set_xlabel(xlabel, fontsize=ax_fontsize)
        ax.set_ylabel('Voltage (V)', fontsize=ax_fontsize)
        ax.grid()
        fig.suptitle(title+'\nApplied Voltage: {} V\nVoltage Threshold to Detect Short: {} V'.format(tvolt, shortThresh), 
                     fontsize=title_fontsize)
        fig.subplots_adjust(hspace=0.5)
        fig.tight_layout(rect=[0, 0, 1, 0.98])
    return fig

def plot_cell_status_convex(cells_input, tvolt=None, deadCell_thresh=None, shortCell_thresh=None, 
                            figsize=None, title_fontsize=16, ax_fontsize=14,
                            title='Cell Status Viewed\nFrom Convex (Non-reflective) Side', 
                            cbarlabel='Controller ID / Cell Status',
                            deadCell_label=None, shortCell_label=None, badIF_label=None, 
                            controllerColors=None, deadCell_boxColor='black', shortCell_boxColor='white', 
                            badIF_markerColor='darkgrey', goodCell_numColor='black', deadCell_numColor='white', 
                            shortCell_numColors=None, badIF_numColor='black', cell_borderColor='black', 
                            badIFtextloc='lower right', badIFtextcoords=(1.56, -0.1), badIFmarkersize=40, 
                            badIFlegendmarkerscale=0.4):
    """
    Plots the voltage cell status for an adjustable optic. cells_voltage_test() should have been run at tvolt on cells
    and cells_short_test() at tvolt should have been run on cells before being passed to this function. The cells with
    unobservable IFs should be set via cell.badIF = True before passing cells to this function.
    
    INPUTS:
    cells_input: list, list of cells that have had cells_voltage_test() and cells_short_test() run on them
    tvolt: float, the test voltage that was used during cells_voltage_test() and cells_short_test()
    deadCell_thresh: the threshold that was used during cells_voltage_test() to determine if the cell is dead
    shortCell_thresh: the threshold used during cells_short_test() to determine if the cell is shorted
    """
    cells = copy.deepcopy(cells_input)
    if figsize is None:
        figsize = (8, 8)
    if deadCell_label is None: # set colorbar labels if None are provided
        if deadCell_thresh is not None:
            deadCell_label = 'Dead Cell: < {} V'.format(tvolt-deadCell_thresh)
        else:
            deadCell_label = 'Dead Cell'
    if shortCell_label is None:
        if shortCell_thresh is not None:
            shortCell_label = 'Shorted Cell: > {} V'.format(shortCell_thresh)
        else:
            shortCell_thresh = 'Shorted Cell'
    if badIF_label is None:
        badIF_label = 'Unobservable IF'
    # get colorbar properties
    cbar_labels, cmap, cbar_label_vals, cbar_tick_bounds, cbar_labelpad = puf.cellStatus_cbar_properties(cells,
                                                                                                         deadCell_label,
                                                                                                         shortCell_label,
                                                                                                         controllerColors, 
                                                                                                         deadCell_boxColor, 
                                                                                                         shortCell_boxColor)
    fig, ax = plt.subplots(figsize=figsize)
    cell_value_array = np.where(~np.isnan(cells[0].no_array), 0, cells[0].no_array)
    image = puf.paint_cellStatus(ax, cells, cell_value_array, cbar_labels, cmap, goodCell_numColor, deadCell_numColor, 
                             shortCell_numColors, badIF_numColor, deadCell_label, shortCell_label, ax_fontsize) # paint cell_value_array
    puf.cellStatus_cellBorders(ax, cells, cell_borderColor, ax_fontsize, title, title_fontsize) # paint cell borders
    puf.cellStatus_cbar(ax, fig, image, cbar_tick_bounds, cbar_label_vals, cbarlabel, # construct colorbar
                        ax_fontsize, cbar_labelpad, cbar_labels)
    puf.mark_badIF_cells(ax, cells, badIF_markerColor, badIF_label, badIFtextloc, badIFtextcoords, 
                         badIFmarkersize, badIFlegendmarkerscale, ax_fontsize)
    return fig