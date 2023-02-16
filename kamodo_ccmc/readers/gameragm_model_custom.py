# Collaboration Example
# ---------------------
# This script outlines the items needed from the modeling team to collaborate
# on getting the model outputs accessible through Kamodo. This includes the
# model_varnames dictionary, the custom interpolation routine, the custom
# coordinate conversion routine (if needed), and a coordinate grid in both the
# standard and custom coordinate systems for plotting purposes.


# Below is an example of the model_varnames metadata dictionary. The keys are
# the names of the variables as given in the files. The list begins with the
# LaTeX representation of the variable and is followed by a description of the
# variable. These are followed by an integer, coordinate information, a list of
# the coordinate names, and the units. See the ConcertCoord function in
# kamodo_ccmc/flythrough/utils.py for more details on the coordinate systems.
# See https://ensemblegovservices.github.io/kamodo-core/notebooks/Syntax/
# for details on naming syntax.
# In the example below, '_ilev1' is added to LaTeX representation of the
# variable that depend on pressure level (the model specific coordinate system)
# to distinguish from the version of the variable that depends on the standard
# coordinate system.
model_varnames = {
    'Bx': [
        'B_x',
        'X-component of magnetic field',
        0,
        'SM',
        'cart',
        ['X', 'Y', 'Z'],
        'nT'
    ],
    'By': [
        'B_y',
        'Y-component of magnetic field',
        1,
        'SM',
        'cart',
        ['X', 'Y', 'Z'],
        'nT'
    ],
    'Bz': [
        'B_x',
        'Z-component of magnetic field',
        2,
        'SM',
        'cart',
        ['X', 'Y', 'Z'],
        'nT'
    ],
    'D': [
        'rho',
        'Number density',
        3,
        'SM',
        'cart',
        ['X', 'Y', 'Z'],
        '#/cm**3'
    ],
    'Jx': [
        'J_x',
        'X-component of current density',
        4,
        'SM',
        'cart',
        ['X', 'Y', 'Z'],
        'nA/m**2'
    ],
    'Jy': [
        'J_y',
        'Y-component of current density',
        5,
        'SM',
        'cart',
        ['X', 'Y', 'Z'],
        'nA/m**2'
    ],
    'Jz': [
        'J_z',
        'Z-component of current density',
        6,
        'SM',
        'cart',
        ['X', 'Y', 'Z'],
        'nA/m**2'
    ],
    'P': [
        'P',
        'Pressure',
        7,
        'SM',
        'cart',
        ['X', 'Y', 'Z'],
        'nPa'
    ],
    'Pb': [
        'Pb',
        'Magnetic Pressure',
        8,
        'SM',
        'cart',
        ['X', 'Y', 'Z'],
        'nPa'
    ],
    'SrcD': [
        'SrcD',
        '???',
        9,
        'SM',
        'cart',
        ['X', 'Y', 'Z'],
        '#/cm**3'
    ],
    'SrcDT': [
        'SrcDT',
        'Internal diagnostic variable',
        10,
        'SM',
        'cart',
        ['X', 'Y', 'Z'],
        's'
    ],
    'SrcP': [
        'SrcP',
        'Internal diagnostic variable',
        11,
        'SM',
        'cart',
        ['X', 'Y', 'Z'],
        'nPa'
    ],
    'SrcX1': [
        'SrcX1',
        'Internal diagnostic variable',
        12,
        'SM',
        'cart',
        ['X', 'Y', 'Z'],
        'deg'
    ],
    'SrcX2': [
        'SrcX2',
        'Internal diagnostic variable',
        13,
        'SM',
        'cart',
        ['X', 'Y', 'Z'],
        'deg'
    ],
    'Vx': [
        'V_x',
        'X-component of velocity',
        14,
        'SM',
        'cart',
        ['X', 'Y', 'Z'],
        'km/s'
    ],
    'Vy': [
        'V_y',
        'Y-component of velocity',
        15,
        'SM',
        'cart',
        ['X', 'Y', 'Z'],
        'km/s'
    ],
    'Vz': [
        'V_z',
        'Z-component of velocity',
        16,
        'SM',
        'cart',
        ['X', 'Y', 'Z'],
        'km/s'
    ],
}

# If the file format is binary, compressed, or otherwise difficult to access,
# we need a function that returns a list of the variable names that are
# available in the chosen file(s).
# The beginning timestamp of each file should be calculable from the filename.
# If more than one timestamp is present in each file, and the file format is
# difficult to access, then we need a funtion that returns the first and last
# timestamp for the file or the interpolating range set in the custom
# interpolator, and the time resolution of the data.

# For each type of coordinate, we will need the associated unit given in a
# separate dictionary as in the example below.
coord_dict = {
    "X": "R_E",
    "Y": "R_E",
    "Z": "R_E",
}


# The 'interp_library' function below (line 72) should return the interpolating
# function for the given variable. The returned interpolating function should
# accept inputs of the form of the standard_track variable in the example
# above, either in numpy array or list form. The return value of the
# interpolation function should be the interpolated value in the MODEL-SPECIFIC
# coordinate system. E.g. value = interp(custom_track)

# Example track
import numpy as np
t = np.linspace(0., 24., 25)
c1 = np.linspace(-5, 25, 15)
c2 = np.linspace(-25, 2, 20)
c3 = np.linspace(-15, 15, 11)
tt, cc1, cc2, cc3 = np.meshgrid(t, c1, c2, c3)
t, c1, c2, c3 = np.ravel(tt), np.ravel(cc1), np.ravel(cc2), np.ravel(cc3)
standard_track = np.array([t, c1, c2, c3]).T
standard_track.shape  # (82500, 4)
# each row of the array/list should be a coordinate set
# e.g. [t_val, c1_val, c2_val, c3_val] for each row


# The names of the functions and script can be chosen by the developing team
# and should reflect the purpose of the  function/script. The name of the
# script should also begin with the name of the model.
def CustomInterp(variable_name):
    '''Custom Interpolator for given variable name in the model-specific
    coordinate system.'''
    # import custom interpolator script, likely via f2py interface
    from model_special import interp_library
    # call interpolator function 
    interp = interp_library(variable_name)
    return interp

# Include an example of how to use the interpolating function
T_interp = CustomInterp('T')
T_values = T_interp(custom_track)  # should return a 1D array or list of values
# custom_track should be an array similar to standard_track but with
# model-specific coordinate values of shape
# (number of positions, number of coordinate values per position)
# The interpolator should cover the entirety of the model grid with the
# exception of time, which can be sectioned off into single or multi-hour
# segments (preferably 24-hour segments if it will fit into 15 GB of memory).
# If the time grid is sectioned off, then the interpolator should
# properly handle interpolations between the time sections.


def CustomCoordConv():
    '''Custom function converting from a satellite track in a standard
    coordinate system to a track in the model-specific coordinate system.'''
    # import custom coordinate conversion script, likely via f2py interface
    from model_special import coord_conv
    return coord_conv

# Include an example of how to use the function, both with and without a 
# variable-specific interpolator
coord_conv = CustomCoordConv()
custom_track = coord_conv(standard_track)
# should return the values of the model-specific
# coordinates for the given coordinates in a standard system.
# def test(t, standard_coord1, standard_coord2, standard_coord3):
#   ---convert coordinates here----
#   return np.array([t, custom_coord1, custom_coord2, custom_coord3, ...])

# The two functions should be able to be used together, as below.
T_values = T_interp(coord_conv(standard_track))


# Finally, I need a set of 1D arrays with default values for the model-specific
# coordinate system. The max and min of the arrays should correspond to the max
# and min of the grid in standard coordinates. These will be used as the
# default coordinate grids for plotting. If these are different for each run,
# then a function should be created to return them, such as the one below.
def Coord_Ranges():
    # import custom script, likely via f2py
    from model_special import coord_ranges
    t, c1, c2, c3 = coord_ranges()
    return t, c1, c2, c3  # can add additional coordinate dimensions as needed


# We will work together on getting the gridded versions of these functions
# working once I get the reader together. Plotting in kamodo is problematic
# at best until the gridded versions of these functions are working.