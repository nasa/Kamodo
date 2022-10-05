# -*- coding: utf-8 -*-
"""
Created on Thu May 13 18:28:18 2021
@author: rringuet
"""
from kamodo import kamodofy, gridify
from numpy import NaN, vectorize, append, array, ravel
from scipy.interpolate import RegularGridInterpolator as rgiND
from scipy.interpolate import interp1d as rgi1D
import forge


def register_interpolator(kamodo_object, varname, interpolator, coord_units):
    '''Register interpolators for each variable.

    Inputs:
        - kamodo_object: A kamodo object produced by the Kamodo core package.
        - varname: A string indicating the standardized variable name
            associated with the given interpolator.
        - interpolator: A kamodofied function produced by one of the functions
            above (gridded or non-gridded).
        - coord_units: A dictionary of key, value pairs indicating the
            argument units. In this case, coord_units =
            {'time':'hr','lon':'deg',...} or similar.
    Output: The same kamodo object given, except with the new function
        included.
    '''

    kamodo_object[varname] = interpolator
    kamodo_object.variables[varname]['xvec'] = coord_units
    kamodo_object._registered += 1
    return kamodo_object


def define_griddedinterp(data_dict, coord_units, coord_data, interp):
    '''Define a gridded interpolator.
    Inputs:
        data_dict: a dictionary containing the data information.
            {'units': 'data_units', 'data': data_array}
            data_array should have the same shape as (c1, c2, c3, ..., cN) 
        coord_units: a dictionary containing the coordinate information.
            {'name_of_coord1': coord1_units', 'name_of_coord2': 'coord2_units',
             etc...}. All units should be strings.
        coord_data: a dictionary containing the coordinate data.
            {'name_of_coord1': coord1_data', 'name_of_coord2': 'coord2_data',
             etc...}. All arrays should be 1D arrays.
        rgi: an interpolator
    Returns: a gridded kamodo interpolator
    '''
    interpolator_grid = kamodofy(gridify(interp, **coord_data),
                                 units=data_dict['units'],
                                 data=data_dict['data'], arg_units=coord_units)
    return interpolator_grid


def create_interp(coord_data, data_dict):
    '''Create an interpolator depending on the dimensions of the input.
    Inputs:
        coord_data: a dictionary containing the coordinate data.
            {'name_of_coord1': coord1_data', 'name_of_coord2': 'coord2_data',
             etc...}. All arrays should be 1D arrays.
        data_dict: a dictionary containing the data information.
            {'units': 'data_units', 'data': data_array}
            data_array should have the same shape as (c1, c2, c3, ..., cN)
        Note: The dataset must also depend upon ALL of the coordinate arrays
            given.
    Output: an interpolator created with the given dataset and coordinates.
    '''
    
    # determine the number of coordinates, create the interpolator
    coord_list = [value for key, value in coord_data.items()]  # list of arrays
    n_coords = len(coord_data.keys())
    if n_coords == 1:
        rgi = rgi1D(*coord_list, data_dict['data'],
                          bounds_error=False, fill_value=NaN)
    else:
        rgi = rgiND(coord_list, data_dict['data'],
                                      bounds_error=False, fill_value=NaN)

    # wrap in a function and return the function
    def interp(xvec):
        return rgi(xvec)
    return interp


def create_funcsig(coord_data, coord_str):
    '''Create a custom function signature based on the dimensions and the
    coordinate string given (e.g. "SMcar").
    Inputs:
        coord_data: a dictionary containing the coordinate data.
            {'name_of_coord1': coord1_data', 'name_of_coord2': 'coord2_data',
             etc...}. All arrays should be 1D arrays.
        coord_str: a string indicating the coordinate system of the data
            (e.g. "SMcar" or "GEOsph").
    Outputs: A forge parameter object.
    '''
    
    # determine the number of coordinates
    n_coords = len(coord_data.keys())
    if n_coords == 1:
        # prepare the custom function signature
        if coord_str != '':
            coord_name = list(coord_data.keys())[0]+'_'+coord_str
        else:
            coord_name = list(coord_data.keys())[0]
    else:
        # prepare the custom function signature
        if 'sph' in coord_str:
            coord_name = 'rvec_'+coord_str+str(n_coords)+'D'
        else:
            coord_name = 'xvec_'+coord_str+str(n_coords)+'D'
    # e.g. 'xvec_SMcar4D' or 'rvecGEOsph3D'
    param_xvec = forge.FParameter(
        name=coord_name, interface_name='xvec',
        kind=forge.FParameter.POSITIONAL_OR_KEYWORD)
    return param_xvec


def Functionalize_Dataset(kamodo_object, coord_dict, variable_name,
                          data_dict, gridded_int, coord_str):
    '''Determine and call the correct functionalize routine.
    Inputs:
        kamodo_object: the previously created kamodo object.
        coord_dict: a dictionary containing the coordinate information.
            {'name_of_coord1': {'units': 'coord1_units', 'data': coord1_data},
             'name_of_coord2': {'units': 'coord2_units', 'data': coord2_data},
             etc...}
            coordX_data should be a 1D array. All others should be strings.
        variable_name: a string giving the LaTeX representation of the variable
        data_dict: a dictionary containing the data information.
            {'units': 'data_units', 'data': data_array}
            data_array should have the same shape as (c1, c2, c3, ..., cN)
        Note: The dataset must also depend upon ALL of the coordinate arrays
            given.
        
        gridded_int: True to create a gridded interpolator (necessary for
            plotting higher dimensions and slicing). False otherwise.
        coord_str: a string indicating the coordinate system of the data
            (e.g. "SMcar" or "GEOsph").

    Output: A kamodo object with the functionalized dataset added.
    '''

    # split the coord_dict into data and units
    coord_data = {key: value['data'] for key, value in coord_dict.items()}
    coord_units = {key: value['units'] for key, value in coord_dict.items()}

    # create interpolator function and function signature
    interp = create_interp(coord_data, data_dict)
    param_xvec = create_funcsig(coord_data, coord_str)

    # Functionalize the dataset
    new_interp = forge.replace('xvec', param_xvec)(interp)
    interp = kamodofy(units=data_dict['units'], data=data_dict['data'],
             arg_units=coord_units)(new_interp)
    kamodo_object = register_interpolator(kamodo_object, variable_name, interp,
                                          coord_units)

    # Add gridded version if requested, even for 1D functions
    if gridded_int:
        interp_grid = define_griddedinterp(data_dict, coord_units, coord_data,
                                           interp)
        kamodo_object.variables[variable_name+'_ijk'] = data_dict
        kamodo_object = register_interpolator(kamodo_object,
                                              variable_name+'_ijk',
                                              interp_grid, coord_units)
    return kamodo_object


def time_interp(kamodo_object, coord_dict, variable_name, data_dict,
                  gridded_int, coord_str):
    '''Create a functionalized interpolator by splitting into timesteps.
    Inputs:
        kamodo_object: the previously created kamodo object.
        coord_dict: a dictionary containing the coordinate information.
            {'name_of_coord1': {'units': 'coord1_units', 'data': coord1_data},
             'name_of_coord2': {'units': 'coord2_units', 'data': coord2_data},
             etc...}
            coordX_data should be a 1D array. All others should be strings.
        variable_name: a string giving the LaTeX representation of the variable
        data_dict: a dictionary containing the data information.
            {'units': 'data_units', 'data': data_array}
            data_array should have the same shape as (c1, c2, c3, ..., cN)
        Note: The dataset must also depend upon ALL of the coordinate arrays
            given.
        
        gridded_int: True to create a gridded interpolator (necessary for
            plotting higher dimensions and slicing). False otherwise.
        coord_str: a string indicating the coordinate system of the data
            (e.g. "SMcar" or "GEOsph").

    Output: A kamodo object with the functionalized dataset added.    
    '''
    
    # create a list of interpolators per timestep, not storing the data, and 
    # interpolate from there. such as done in superdarnea_interp.py
    # Assumes that time is the first dimension
    # Assumes that data_dict['data'] is a cdf_data.variable object
    # does this also work for h5 files? ******************************************************

    # split the coord_dict into data and units
    coord_data = {key: value['data'] for key, value in coord_dict.items()}
    coord_units = {key: value['units'] for key, value in coord_dict.items()}
    coord_list = [value for key, value in coord_data.items()]  # list of arrays
    
    # Create ND interpolators for each time grid value
    if len(coord_list) > 2:
        time_interps = [rgiND(coord_list[1:], array(data_dict['data'][i]),
                              bounds_error=False, fill_value=NaN) for i in
                        range(len(coord_list[0])-1)]
    elif len(coord_list) <= 2:  # call normal routine and return
        data_dict['data'] = array(data_dict['data'])  # ensure it is an array
        return Functionalize_Dataset(kamodo_object, coord_dict, variable_name,
                                  data_dict, gridded_int, coord_str)

    @vectorize
    def interp_i(*args):

        position = [*args]  # time coordinate value must be first
        # Choose indices of time grid values surrounding desired time
        idx = sorted(append(coord_list[0], position[0])).index(position[0])  # get location
        if idx > 1 and idx < len(coord_list[0])-1:  # middle somewhere
            idx_list = [idx-1, idx, idx+1]
        elif idx <= 1:  # beginning
            idx_list = [0, 1, 2]
        elif idx > 1 and idx >= len(coord_list[0])-1:  # at end
            idx_list = [idx-2, idx-1, idx]

        # Interpolate values for chosen latitude grid values
        interp_values = ravel([time_interps[i]([*position[1:]]) for i in idx_list])
        time_interp = rgi1D(coord_list[0][idx_list], interp_values,
                              bounds_error=False, fill_value=NaN)
        return time_interp(position[0])

    # retrieve the updated function signature
    param_xvec = create_funcsig(coord_data, coord_str)

    @forge.replace('xvec', param_xvec)
    @kamodofy(units=data_dict['units'], data=data_dict['data'],
              arg_units=coord_units)
    def total_interp(xvec):
        return interp_i(*array(xvec).T)

    # register the interpolator in the kamodo object
    kamodo_object = register_interpolator(kamodo_object, variable_name,
                                          total_interp, coord_units)

    # Add gridded version if requested, even for 1D functions
    if gridded_int:
        interp = define_griddedinterp(data_dict, coord_units, coord_data,
                                      total_interp)
        kamodo_object.variables[variable_name+'_ijk'] = data_dict
        kamodo_object = register_interpolator(kamodo_object,
                                              variable_name+'_ijk', interp,
                                              coord_units)
    return kamodo_object