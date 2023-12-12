import numpy as np

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