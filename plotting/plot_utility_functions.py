from axroGalil.electronics import get_controller_ids
from axroGalil.electronics import filter_cells
import numpy as np
import random
import copy
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable

def cellStatus_cbar_properties(cells, deadCell_label, shortCell_label, controllerColors, deadCell_boxColor, 
                               shortCell_boxColor):
    """
    Helper function for plot_cell_status_convex(). Returns the colorbar labels, the colormap, the colorbar ticks, and
    the colorbar title padding to construct the colorbar.
    """
    cids = get_controller_ids(cells)
    cbar_labels = [shortCell_label, deadCell_label] + cids
    if controllerColors is None: # get colors for controllers
        if len(cids) <= len(colors_36): # if there's less or equal to 36 controllers being used
            controllerColors = copy.deepcopy(colors_36[:len(cids)])
        elif (len(cids) > len(colors_36)) and (len(cids) <= len(rand_nonblack_colors)): # if there's less or equal to 148 controllers being used
            controllerColors = copy.deepcopy(rand_nonblack_colors[len(cids)])
        else:
            raise PlottingError("The number of controllers used ({}) exceeds the number of standard colors available ({})"
                                           .format(len(cids), len(rand_nonblack_colors)))
    cbar_colors = [shortCell_boxColor, deadCell_boxColor] + controllerColors
    shorts = any(cell.short_cell_nos.size > 0 for cell in cells) # bool for if any cells are shorted
    nonShort_deadCells = any(cell.deadCell==True and cell.short_cell_nos.size==0 for cell in cells) # bool for if any cells are dead but not shorted
    if shorts and nonShort_deadCells:
        print('There are shorts and exclusively dead cells.')
        labelpad = -40
    elif shorts and not nonShort_deadCells:
        print('There are shorts but no exclusivley dead cells.')
        cbar_labels.pop(1)
        cbar_colors.pop(1)
        labelpad = -10
    elif not shorts and nonShort_deadCells:
        print('There are no shorts but there are exclusivley dead cells.')
        cbar_labels.pop(0)
        cbar_colors.pop(0)
        labelpad = -40
    else:
        print('There are no shorts and no exclusively dead cells.')
        cbar_labels.pop(0)
        cbar_colors.pop(0)
        cbar_labels.pop(0)
        cbar_colors.pop(0)
        labelpad = 10
    cmap = ListedColormap(cbar_colors)
    cbar_label_vals = np.arange(len(cbar_labels))
    cbar_tick_bounds = np.arange(cbar_label_vals[0]-0.5, cbar_label_vals[-1]+1, 1)
    return cbar_labels, cmap, cbar_label_vals, cbar_tick_bounds, labelpad

def get_rand_nonblack_colors():
    """
    Returns a randomized list of CSS colors excluding colors in the black to white range.
    """
    colors = mcolors.CSS4_COLORS
    sorted_colors = sorted(colors, key=lambda c: tuple(mcolors.rgb_to_hsv(mcolors.to_rgb(c))))
    nonblack_colors = sorted_colors[14:]
    rand_nonblack_colors = random.sample(nonblack_colors, len(nonblack_colors))
    return rand_nonblack_colors

colors_36 = ['red', 'dodgerblue', 'limegreen', 'mediumorchid', 'yellow', 'chocolate', 'pink', 'mediumturquoise', 
             'coral', 'skyblue', 'forestgreen', 'mediumslateblue', 'khaki', 'burlywood', 'mediumvioletred', 'aqua', 
             'rosybrown', 'steelblue', 'palegreen', 'plum', 'gold', 'darkorange', 'magenta', 'teal', 
             'firebrick', 'blue', 'darkseagreen', 'thistle', 'darkgoldenrod', 'olive', 'thistle', 'deepskyblue', 
             'greenyellow', 'darkblue', 'saddlebrown', 'darkorange']

rand_nonblack_colors = get_rand_nonblack_colors()

def paint_cellStatus(ax, cells, cell_value_array, cbar_labels, cmap, goodCell_numColor, deadCell_numColor, 
                     shortCell_numColors, badIF_numColor, deadCell_label, shortCell_label, ax_fontsize):
    """
    Helper function for plot_cell_status_convex(). Paints the cell_value_array the correct colors.
    """
    if shortCell_numColors is None:
        shortCell_numColors = list(copy.deepcopy(mcolors.TABLEAU_COLORS).keys()) + \
                              list(copy.deepcopy(mcolors.BASE_COLORS).keys())
    short_groups = get_cells_short_groups(cells)
    for cell in cells:
        if cell.deadCell == True and cell.short_cell_nos.size == 0: # cell is dead but not shorted
            cell_value_array[cell.coords] = cbar_labels.index(deadCell_label)
            ax.text(cell.coords[1][0], cell.coords[0][0], int(cell.no), color=deadCell_numColor, 
                    ha='center', va='center', fontsize=ax_fontsize-6, zorder=5.)
        elif cell.short_cell_nos.size > 0: # cell is shorted
            cell_value_array[cell.coords] = cbar_labels.index(shortCell_label)
            for i, short_group in enumerate(short_groups):
                if cell.no in short_group:
                    ax.text(cell.coords[1][0], cell.coords[0][0], int(cell.no), color=shortCell_numColors[i],
                            ha='center', va='center', fontsize=ax_fontsize-6)
                    ax.add_patch(Rectangle((cell.coords[1][0]-0.5, cell.coords[0][0]-0.5), 1, 1,
                                edgecolor=shortCell_numColors[i], facecolor="None", lw=1.5,
                                zorder=5.0))
                    break
        else: # cell has no shorts and is not dead
            cell_value_array[cell.coords] = cbar_labels.index(cell.cid)
            if cell.badIF:
                numColor = badIF_numColor
            else:
                numColor = goodCell_numColor
            ax.text(cell.coords[1][0], cell.coords[0][0], int(cell.no), color=numColor, 
                    ha='center', va='center', fontsize=ax_fontsize-6, zorder=5.)
    image = ax.imshow(cell_value_array, aspect='equal', cmap=cmap)
    return image

def get_cells_short_groups(cells):
    """
    Returns a list of lists where each sub-list contains the cell numbers that forms a short group.
    """
    short_groups = [list(cell.short_cell_nos) for cell in cells if cell.short_cell_nos.size>0]
    if not short_groups:
        unique_short_groups = None
    else:
        seen = set()
        unique_short_groups = [tuple(short_group) for short_group in short_groups if tuple(short_group) not 
                               in seen and not seen.add(tuple(short_group))] # preserve unique short groups
        unique_short_groups = [list(short_group) for short_group in unique_short_groups]
    return unique_short_groups

def get_cell_short_groups_voltages(cells):
    """
    Returns the voltages that the short groups were at when the short was measured (list of lists of voltages).
    """
    short_groups = get_cells_short_groups(cells)
    if short_groups is None:
        short_groups_voltages = None
    else:
        short_groups_voltages = []
        for i, short_group in enumerate(short_groups):
            short_group_cell = filter_cells(cells, 'no', short_group[0])[0]
            short_groups_voltages.append(list(short_group_cell.short_volts))
    return short_groups_voltages

def cellStatus_cellBorders(ax, cells, cell_borderColor, ax_fontsize, title, title_fontsize):
    """
    Helper function for plot_cell_status_convex. Draws the borders between cells and labels the x and y axes.
    """
    row_nums = np.array([i for i in range(cells[0].no_array.shape[0])])
    col_nums = np.array([i for i in range(cells[0].no_array.shape[1])])
    hline_locs = (row_nums[1:]+row_nums[:-1]) / 2
    vline_locs = (col_nums[1:]+col_nums[:-1]) / 2
    for loc in hline_locs: 
        ax.axhline(loc, color=cell_borderColor)
    for loc in vline_locs: 
        ax.axvline(loc, color=cell_borderColor)
    ax.set_yticks(row_nums)
    ax.set_xticks(col_nums)
    ax.set_yticklabels(row_nums+1)
    ax.set_xticklabels(col_nums+1)
    ax.set_ylabel('Row', fontsize=ax_fontsize)
    ax.set_xlabel('Column', fontsize=ax_fontsize)
    ax.set_title(title, fontsize=title_fontsize)
    
def cellStatus_cbar(ax, fig, image, cbar_tick_bounds, cbar_label_vals, cbar_title, ax_fontsize, labelpad, 
                    cbar_labels):
    """
    Constructs the colorbar for plot_cell_status_convex().
    """
    div = make_axes_locatable(ax)
    cax = div.append_axes('right', size='5%', pad=0.50)
    cbar = fig.colorbar(image, cax=cax, boundaries=cbar_tick_bounds, ticks=cbar_label_vals)
    cbar.ax.set_yticklabels(cbar_labels, fontsize=ax_fontsize)
    cbar.set_label(cbar_title, fontsize=ax_fontsize, labelpad=10)
    cbar.ax.yaxis.set_label_position('left')
    fig.subplots_adjust(top=0.85, hspace=0.5, wspace=0.8)

def mark_badIF_cells(ax, cells, badIF_markerColor, badIF_label, legendloc, legendcoords, markersize, 
                     legendmarkerscale, ax_fontsize):
    """
    Helper function for plot_cell_status_convex(). Marks the cells whose IFs are unobservable.
    """
    # get the cells that have bad IFs, but are not dead and not shorted
    badIF_cells = [cell for cell in cells if cell.badIF==True and cell.deadCell==False and cell.short_cell_nos.size==0]
    if badIF_cells:
        y_coords = [cell.coords[0][0] for cell in badIF_cells]
        x_coords = [cell.coords[1][0] for cell in badIF_cells]
        ax.plot(x_coords, y_coords, color=badIF_markerColor, marker="X", label=badIF_label, 
                linestyle="None", markersize=markersize, zorder=1.)
        legend = ax.legend(loc=legendloc, bbox_to_anchor=(legendcoords[0], legendcoords[1]), 
                           fancybox=True, fontsize=ax_fontsize, markerscale=legendmarkerscale)



def get_imbounds(imbounds, nos, array_ls):
    """
    Formats the user provided imbounds to get the appropriate images to display
    from the array stacks.
    """
    displaySingle = False
    if isinstance(imbounds, int):
        imbounds = [imbounds, imbounds]
    elif isinstance(imbounds, list) and (len(imbounds) == 1):
        imbounds += imbounds
    if imbounds is not None: # the user provided imbounds
        if imbounds[0] == imbounds[1]:
            displaySingle = True # show only a single frame
        if type(nos) != type(None): # match the imbounds to the cell numbers
            try:
                lowerBound = int(np.where(nos == imbounds[0])[0])
                upperBound = int(np.where(nos == imbounds[1])[0] + 1)
            except:
                raise PlottingError("One or both imbounds provided {} is not in the given number labels.".format(imbounds))
        else: # explicitly use the imbounds as indexes when nos are not present
            lowerBound, upperBound = imbounds[0], imbounds[1]+1
    else: # show all images supplied
        lowerBound, upperBound = 0, array_ls[0].shape[0]
        if array_ls[0].shape[0] == 1: 
            displaySingle = True
    return lowerBound, upperBound, displaySingle

class PlottingError(Exception):
    """
    Custom exception related to plotting data.
    """
    pass