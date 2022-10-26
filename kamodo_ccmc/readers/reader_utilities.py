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
                          data_dict, gridded_int, coord_str, interp_flag=0,
                          func=None, start_idx=None):
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
        interp_flag: the method of interpolation required. 
            Options: 
                0: (default option) assumes the given data_dict['data'] is an
                    N-dimensional numpy array and creates a standard scipy
                    interpolator to functionalize the data. 
                1: Assumes the given data_dict['data'] is a list of cdf or h5
                    file objects with data for one time slice per object.
                    Lazy interpolation is applied on top of the standard scipy
                    interpolators (one per time slice).
                2: Assumes the given data_dict['data'] is a list of cdf or h5
                    file objects with data for a range of times per object.
                    Lazy chunked interpolation is applied on top of the
                    standard scipy interpolators (one per time chunk).
                    Interpolation between time chunks occurs automatically.
        func: a function defining the logic to be executed on a given time
            slice or chunk (e.g. converting to an array and then transposing
            it). 
            *** Only needed for interp_flag greater than zero. ***
            - For interp_flag=1, the default (None) is to convert to a numpy
            array.
            - For interp_flag=2, the default (None) is to do nothing. The
            chunks will be converted to arrays during the preparation to
            interpolate between time chunks.
            - This option is useful when a dataset's shape order does not match
            the required order (e.g. t, lat, lon -> t, lon, lat) and allows
            such operations to be done on the fly.
        start_idx: a list of indices indicating the position of the start times
            for each data chunk in the time grid (coord_dict['time']['data']).
            *** Only needed for interpolation over time chunks. ***
            The length should be one longer than the number of chunks, with the
            last value equal to the length of the time grid.

    Output: A kamodo object with the functionalized dataset added.
    '''

    # split the coord_dict into data and units
    coord_data = {key: value['data'] for key, value in coord_dict.items()}
    coord_units = {key: value['units'] for key, value in coord_dict.items()}

    # create a functionalized interpolator and modified function signature
    param_xvec = create_funcsig(coord_data, coord_str)
    if interp_flag == 0:  # standard logic
        interp = create_interp(coord_data, data_dict)
    elif interp_flag == 1:
        interp = time_interp(coord_dict, data_dict, func)
    elif interp_flag == 2:
        interp = multitime_interp(coord_dict, data_dict, start_idx, func)
    new_interp = forge.replace('xvec', param_xvec)(interp)
    interp = kamodofy(units=data_dict['units'], data=data_dict['data'],
             arg_units=coord_units)(new_interp)

    # Register and add gridded version if requested, even for 1D functions
    kamodo_object = register_interpolator(kamodo_object, variable_name, interp,
                                          coord_units)
    if gridded_int:
        interp_grid = define_griddedinterp(data_dict, coord_units, coord_data,
                                           interp)
        kamodo_object.variables[variable_name+'_ijk'] = data_dict
        kamodo_object = register_interpolator(kamodo_object,
                                              variable_name+'_ijk',
                                              interp_grid, coord_units)
    return kamodo_object


def time_interp(coord_dict, data_dict, func=None):
    '''Create a functionalized interpolator by splitting into timesteps.
    Inputs:
        coord_dict: a dictionary containing the coordinate information.
            {'name_of_coord1': {'units': 'coord1_units', 'data': coord1_data},
             'name_of_coord2': {'units': 'coord2_units', 'data': coord2_data},
             etc...}
            coordX_data should be a 1D array. All others should be strings.
        data_dict: a dictionary containing the data information.
            {'units': 'data_units', 'data': data_array}
            data_array should have the same shape as (c1, c2, c3, ..., cN)
            Note: The dataset must also depend upon ALL of the coordinate
                arrays given.
        func: a function defining the logic to be executed on a given time
            slice (e.g. converting to an array and then transposing it).
            Default is to only convert to an array (e.g. func=None).

    Output: A lazy time interpolator. 
    '''

    # create a list of interpolators per timestep, not storing the data, and 
    # interpolate from there. such as done in superdarnea_interp.py
    # Assumes that time is the first dimension
    # Assumes that data_dict['data'] is a cdf_data.variable object
    # does this also work for h5 files? ******************************************************

    # split the coord_dict into data and units
    coord_list = [value['data'] for key, value in coord_dict.items()]  # list of arrays
    idx_map, time_interps = [], []  # initialize variables
    if func == None:  # define the default operation
        def func(cdf_data_object):
            return array(cdf_data_object)

    # define method for how to add time slices to memory
    def add_timeinterp(i):
        idx_map.append(i)
        if len(coord_list) > 1:
            time_interps.append(rgiND(coord_list[1:],
                                      func(data_dict['data'][i]),
                                      bounds_error=False, fill_value=NaN))
        else:  # has to be different for time series data
            ''''***TEST THIS. DON'T NEED AN INTERPOLATOR FOR A SINGLE POINT!***'''
            def dummy_1D(spatial_position):
                return func(data_dict['data'][i])
            time_interps.append(dummy_1D)
        return time_interps, idx_map

    @vectorize
    def interp_i(*args, time_interps=time_interps,
                 idx_map=idx_map):

        position = [*args]  # time coordinate value must be first
        # Choose indices of time grid values surrounding desired time
        idx = sorted(append(coord_list[0], position[0])).index(position[0])  # get location
        if idx > 1 and idx < len(coord_list[0])-1:  # middle somewhere
            idx_list = [idx-1, idx]
        elif idx <= 1:  # beginning
            idx_list = [0, 1]
        elif idx > 1 and idx >= len(coord_list[0])-1:  # at end
            idx_list = [len(coord_list[0])-2, len(coord_list[0])-1]

        # Interpolate values for chosen time grid values
        for i in idx_list:
            if i not in idx_map:
                len_map = len(idx_map)  # length before adding another
                print(idx, idx_list, idx_map)
                print('Adding idx ', i)
                try:  # allow for memory errors
                    time_interps, idx_map = add_timeinterp(i)
                except MemoryError:  # remove one time slice first
                    print('Avoiding memory error...')
                    if len_map == 0:
                        print('Not enough memory to load one time slice. ' +
                              'Please close some applications.')
                    elif len_map < 2:  # tried to add a second slice but failed
                        print('Not enough memory to load two time slices. ' +
                              'Please close some applications.')
                    elif abs(i-idx_map[0]) > 2:
                        del idx_map[0], time_interps[0]
                    else:
                        del idx_map[-1], time_interps[-1]
                    time_interps, idx_map = add_timeinterp(i)                   
        interp_location = [idx_map.index(val) for val in idx_list]
        interp_values = ravel(array([time_interps[i]([*position[1:]])
                                     for i in interp_location]))
        time_interp = rgi1D(coord_list[0][idx_list], interp_values,
                              bounds_error=False, fill_value=NaN)
        return time_interp(position[0])

    def interp(xvec):
        return interp_i(*array(xvec).T)
    return interp


def multitime_interp(coord_dict, data_dict, start_times, func):
    '''Create a functionalized interpolator by splitting into timesteps.
    Inputs:
        coord_dict: a dictionary containing the coordinate information.
            {'name_of_coord1': {'units': 'coord1_units', 'data': coord1_data},
             'name_of_coord2': {'units': 'coord2_units', 'data': coord2_data},
             etc...}  # time not given
            coordX_data should be a 1D array. All others should be strings.
        data_dict: a dictionary containing the data information.
            {'units': 'data_units', 'data': empty}
            data_array should have the same shape as (c1, c2, c3, ..., cN)
            Note: The dataset must also depend upon ALL of the coordinate
                arrays given.
        start_times: the start times of each file in the file_dir.
            Interpolation between time chunks should be taken care of by appending
            the first time slice of the next time chunk to the end of the
            current time chunk (see func description).
        func: a function defining the logic to be executed on a given time
            chunk, including retrieving the data from the file.

    Output: A time-chunked lazy interpolator.
    '''
    
    # create a list of interpolators per time chunk, not storing the data, and 
    # interpolate from there. such as done in superdarnea_interp.py
    # Assumes that time is the first dimension.
    # Assumes that data_dict['data'][i] is a cdf_data.variable object with 
    # multiple time slices.
    # does this also work for h5 files???????????????????????????????????????????????????

    # split the coord_dict into data and units, initialize variables/func
    coord_list = [value['data'] for key, value in coord_dict.items()]  # list of arrays
    idx_map, time_interps = [], []  # initialize variables

    # define method for how to add time chunks to memory
    def add_timeinterp(i):
        idx_map.append(i)
        # set up to interpolate between time chunks
        if i < len(start_times)-1:  # append 1st time slice from next chunk
            data, time = func([i, i+1])
        else:
            data, time = func([i])
        if len(coord_list) > 1:
            coord_list_i = [time] + coord_list
            time_interps.append(rgiND(coord_list_i, data,
                                      bounds_error=False, fill_value=NaN))
        else:  # has to be different for time series data
            time_interps.append(rgi1D(time, data, bounds_error=False,
                                      fill_value=NaN))
        return time_interps, idx_map
    
    @vectorize
    def interp_i(*args, time_interps=time_interps, idx_map=idx_map):

        position = [*args]  # time coordinate value must be first
        # get location of interpolator
        idx = sorted(append(start_times, position[0])).index(position[0])
        if idx == 0:
            i = 0  # beg of first file
        else:
            i = idx - 1  # somewhere else

        # Interpolate values for chosen time chunk
        if i not in idx_map:
            print(idx, i, idx_map)
            print('Adding idx ', i)
            len_map = len(idx_map)
            try:  # allow for memory errors
                time_interps, idx_map = add_timeinterp(i)
            except MemoryError:  # remove two(?) items first
                print('Avoiding memory error...')
                if len_map == 0:  # tried but failed
                    print('Not enough memory to load a time chunk. ' +
                          'Please close some applications.')
                elif i != idx_map[0]:
                    del idx_map[0], time_interps[0]
                else:
                    del idx_map[-1], time_interps[-1]
                time_interps, idx_map = add_timeinterp(i)             
        interp_location = idx_map.index(i)
        #print(time_interps[interp_location])
        #print(*position)
        #print(time_interps[interp_location](position))
        return time_interps[interp_location](position)  # single time chunk

    def interp(xvec):
        return interp_i(*array(xvec).T)

    return interp
