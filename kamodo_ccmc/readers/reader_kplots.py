# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 15:58:53 2021

@author: rringuet
"""
from numpy import float32, float64, array
from kamodo import kamodofy, partial, Kamodo, gridify


def convert_to_array(value):
    '''Check for floats and integers. Convert to arrays if found.

    Input:
        value: float, integer, np.float32, no.float64, list or numpy array
    Output:
        numpy array containing the same value(s)
    '''

    type_check = [isinstance(value, float), isinstance(value, int),
                  isinstance(value, float32), isinstance(value, float64),
                  isinstance(value, list)]
    if sum(type_check) > 0:
        return array([value])
    else:
        return value


def plot2D(kamodo_object, varname, plottype, t, lon, lat, h=-1):
    '''Use Kamodo's native plotting to generate 2D plot.

    Inputs:
    -------
    kamodo_object: A model reader instance containing at least one function.
    varname: A string indicating the standardized name of the functionalized
        variable in the given kamodo_object. Variable names ending in '_ijk'
        are assumed to be gridded functions.
    t: A value, list, or array of time values in hours.
    lon: A value, list, or array of longitude or X values in deg or R_E.
    lat: A value, list, or array of latitude or Y values in deg or R_E.
    h: A value, list, or array of vertical, pressure level, or Z values in deg,
        unitless, or R_E.
    Note: lon, lat, and h also double as x, y, and z for cartesian inputs.
    plottype: A string indicating the plot type. Possible plot types are
        LonLat, LatH, LonH, TimeLat, TimeLon, and TimeH for spherical
        coordinates; and TimeX, TimeY, TimeZ, XY, XZ, and YZ for
        cartesian coordinates.

    If the variable depends on 4 dimensions, h must be given.
    If a LonLat plot is requested, then the function expects a single value
        (integer, float, float32, or float64) for t and h (if h is given).
        In this case, lon and lat should be 1D arrays or flat lists. Similar
        data formatting is required for coordinates not plotted for all plot
        types.

    Output: Input object to an iplot function
        (from plotly.offline import iplot)
    '''

    # initialize new kamodo object
    plot_kamodo = Kamodo()

    # determine if kamodo function is griddified or not, and function units
    gridified = (varname[-3:] == 'ijk')
    units = kamodo_object.variables[varname]['units']
    xvec = kamodo_object.variables[varname]['xvec']

    # determine vertical dependency of variable
    coord_list = list(xvec.keys())
    if len(coord_list) == 4:
        vert = coord_list[-1]  # height, ilev, ilev1, or milev (always last)
    else:
        vert = 'none'
        if 'H' in plottype or 'Z' in plottype:
            raise AttributeError(
                f'Cannot produce {plottype} plot for a variable ' +
                f'that does not depend on height.\n{varname}: {xvec}\n')

    # convert inputs to arrays
    t = convert_to_array(t)
    lon = convert_to_array(lon)  # doubles as x
    lat = convert_to_array(lat)  # doubles as y
    h = convert_to_array(h)  # doubles as z and pressure level

    # save data in a dictionary with keys named with coordinate dependencies
    data = {coord_list[0]: t, coord_list[1]: lon, coord_list[2]: lat}
    if len(coord_list) == 4:
        data[coord_list[3]] = h

    # create printing message for heading of plot
    if t.shape[0] == 1:
        t_message = f'Time slice at {t[0]:.3f} hrs. '
    else:
        t_message = ''
    if lon.shape[0] == 1:
        if 'z' in vert:
            lon_message = f'X slice at {lon[0]:.3f} R_E. '
        else:
            lon_message = f'Longitude slice at {lon[0]:.3f} deg. '
    else:
        lon_message = ''
    if lat.shape[0] == 1:
        if 'z' in vert:
            lat_message = f'Y slice at {lat[0]:.3f} R_E. '
        else:
            lat_message = f'Latitude slice at {lat[0]:.3f} deg. '
    else:
        lat_message = ''
    if vert == 'none':
        h_message = ''
    elif h.shape[0] > 1:
        h_message = ''
    else:
        if vert in ['ilev', 'ilev1', 'milev']:
            h_message = f'Pressure level slice at {h[0]}.'
        elif vert == 'height':
            h_message = f'Height slice at {h[0]:.3f} km.'
        elif vert == 'radius':
            h_message = f'Radius slice at {h[0]:.7f} R_E.'
        elif 'z' in vert:
            h_message = f'Z slice at {h[0]:.7f} R_E.'
    print(t_message+lon_message+lat_message+h_message)

    # create 2D kamodo function for plotting desired plottype
    # plotting_dict is a data dictionary with the numpy arrays for the grids
    # partial_dict is a data dictionary with the values of two coordinates
    if plottype in ['LonLat', 'XY']:
        plotting_dict = {coord_list[1]: data[coord_list[1]],
                         coord_list[2]: data[coord_list[2]]}
        if vert != 'none':
            partial_dict = {coord_list[0]: data[coord_list[0]],
                            coord_list[3]: data[coord_list[3]]}
        else:
            partial_dict = {coord_list[0]: data[coord_list[0]]}
    elif plottype in ['TimeLon', 'TimeX']:
        plotting_dict = {coord_list[0]: data[coord_list[0]],
                         coord_list[1]: data[coord_list[1]]}
        if vert != 'none':
            partial_dict = {coord_list[2]: data[coord_list[2]],
                            coord_list[3]: data[coord_list[3]]}
        else:
            partial_dict = {coord_list[2]: data[coord_list[2]]}
    elif plottype in ['TimeLat', 'TimeY']:
        plotting_dict = {coord_list[0]: data[coord_list[0]],
                         coord_list[2]: data[coord_list[2]]}
        if vert != 'none':
            partial_dict = {coord_list[1]: data[coord_list[1]],
                            coord_list[3]: data[coord_list[3]]}
        else:
            partial_dict = {coord_list[1]: data[coord_list[1]]}
    elif plottype in ['TimeH', 'TimeZ']:
        plotting_dict = {coord_list[0]: data[coord_list[0]],
                         coord_list[3]: data[coord_list[3]]}
        partial_dict = {coord_list[1]: data[coord_list[1]],
                        coord_list[2]: data[coord_list[2]]}
    elif plottype in ['LonH', 'XZ']:
        plotting_dict = {coord_list[1]: data[coord_list[1]],
                         coord_list[3]: data[coord_list[3]]}
        partial_dict = {coord_list[0]: data[coord_list[0]],
                        coord_list[2]: data[coord_list[2]]}
    elif plottype in ['LatH', 'YZ']:
        plotting_dict = {coord_list[2]: data[coord_list[2]],
                         coord_list[3]: data[coord_list[3]]}
        partial_dict = {coord_list[0]: data[coord_list[0]],
                        coord_list[1]: data[coord_list[1]]}
    if not gridified:
        plot_kamodo[varname] = kamodofy(partial(gridify(
            kamodo_object[varname], **data), **partial_dict),
            units=units, arg_units=xvec)
    else:  # use non-gridified version instead
        plot_kamodo[varname] = kamodofy(partial(gridify(
            kamodo_object[varname[:-4]], **data), **partial_dict),
            units=units, arg_units=xvec)
    return plot_kamodo.plot(**{varname: plotting_dict})
