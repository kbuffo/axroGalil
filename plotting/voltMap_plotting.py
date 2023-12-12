import plot_utility_functions as puf
import numpy as np
from scipy import ndimage as nd
from scipy.interpolate import griddata
from operator import itemgetter
from itertools import chain
import copy
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.animation as animation
plt.rcParams['savefig.facecolor']='white'

def displayVoltMaps(voltMaps_ls_input,
                    imbounds=None, 
                    vbounds=None,
                    cell_num_arrays=None,
                    suppress_thresh=None,
                    map_nos=None,
                    colormap='inferno',
                    plot_titles=None,
                    global_title=None,
                    global_title_y_pos=0.94,
                    voltMaplabel='Map:',
                    ylabel='Row',
                    xlabel='Column',
                    cbarlabel='Voltage (V)',
                    global_title_fontsize=14,
                    plot_title_fontsize=14,
                    ax_fontsize=12,
                    includeMeanVolt=False,
                    date=None,
                    cell_border_color='white',
                    cell_num_lowV_color='white',
                    cell_num_highV_color='black',
                    figsize=None,  
                    frame_time=1000, 
                    repeat_ani=False):
    """
    Displays a list of voltage map stacks.
    voltMaps_ls_input: list, list of numpy arrays that represent the voltage values of cells. Each element
        of the list will generate an independent subplot.
        The elements of voltMaps_ls_input should all be either 2D arrays or 3D arrays. If they are all
        3D arrays, they should all have same length along the last axis (which is the number of frames
        to render).
    imbounds (optional): int or list, specifies which frame(s) of voltMaps_ls_input to render if the 
        elements of voltMaps_ls_input are 3D arrays. If a int j is provided, only a single frame is 
        rendered of voltMaps_ls_input[:][j]. If a list [j,k] is provided, a range of frames is rendered 
        from voltMaps_ls_input[:][j:k]. If map_nos is None, imbounds accesses the indexes of 
        voltMaps_ls_input. If map_nos is provided, imbounds then refers to map_nos to get the correct 
        frames to render.
    vbounds (optional): None or list or list of lists. Sets the data range to display on the colorbar.
        If vbounds is None, separate vbounds are generated for each subplot. If vbounds is a list, all
        subplots will have the same vbounds. If vbounds is a list of lists, each subplot will have its own
        subplot.
    cell_num_arrays (optional): list, list of 1D and/or 2D arrays that is the same length of 
        voltMaps_ls_input that specifies the cell number locations for each subplot.
    suppress_thresh (optional): float, if specified, all voltage values below the suppress_thresh will not
        be displayed.
    map_nos (optional): np array, 1D array that describes the number labels to use when displaying frames
        of the voltMaps.
    colormap: str, specifies the colormap to display the voltMaps.
    plot_titles (optional): list, specifies the subplot titles for voltMaps_ls_input
    global_tile (optional): str, sets the overall title for the figure.
    globab_title_y_pos: float, sets the y position of the global title on the figure.
    voltMaplabel (optional): str or None, specifies the label for map_nos. If None, the voltMap number will
        not be animated.
    ylabel (optional): str or list, sets the y axis label. Each subplot y axis label can be set separately
        if a list is passed.
    xlabel (optional): str or list, sets the x axis label. Each subplot x axis label can be set separately
        if a list is passed.
    cbarlabel (optional): str or list, sets the colorbar label. Each subplot colorbar label can be set
        separately if a list is passed.
    global_title_fontsize: int or float, sets the fontsize of the global title of the figure.
    plot_title_fontsize: int or float, sets the fontsize of each supblot's title.
    ax_fontsize: int or float, sets the fontsize of ylabel, xlabel, and cbarlabel for each subplot.
    includeMeanVolt (optional): bool or list, adds text that specifies the mean voltage value of each 
        subplot
    date (optional): str, if specified, adds the date to the figure
    """
    voltMaps = [np.copy(voltMap) for voltMap in voltMaps_ls_input]
    N_plots = len(voltMaps)

    pass