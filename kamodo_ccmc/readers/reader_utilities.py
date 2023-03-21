# -*- coding: utf-8 -*-
"""
Created on Thu May 13 18:28:18 2021
@author: rringuet
"""
from kamodo import kamodofy, gridify
from numpy import NaN, vectorize, append, array, meshgrid, ravel, diff
from numpy import unique, zeros, ndarray, floor, float32, all, insert
from scipy.interpolate import RegularGridInterpolator as rgiND
from scipy.interpolate import interp1d as rgi1D
import forge
from datetime import datetime, timezone, timedelta
from os.path import basename


def create_interp(coord_data, data_dict, func=None, func_default='data'):
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
    if n_coords == 1 and func_default == 'data':
        rgi = rgi1D(*coord_list, data_dict['data'], bounds_error=False,
                    fill_value=NaN)
        # wrap in a function and return the function
        def interp(xvec):
            return rgi(xvec)
        return interp
    elif n_coords > 1 and func_default == 'data':
        rgi = rgiND(coord_list, data_dict['data'], bounds_error=False,
                    fill_value=NaN)
        # wrap in a function and return the function
        def interp(xvec):
            return rgi(xvec)
        return interp
    elif func_default == 'custom':
        return func  # returns the custom interpolator


def create_funcsig(coord_data, coord_str, bounds):
    '''Create a custom function signature based on the dimensions and the
    coordinate string given (e.g. "SMcar").
    Inputs:
        coord_data: a dictionary containing the coordinate data.
            {'name_of_coord1': coord1_data', 'name_of_coord2': 'coord2_data',
             etc...}. All arrays should be 1D arrays.
        coord_str: a string indicating the coordinate system of the data
            (e.g. "SMcar" or "GEOsph").
        bounds: an array of positions for the max and min values of each
            coordinate grid transposed.
            aa, bb, cc, dd = np.meshgrid(trange, lon_range, lat_range, h_range)
            a, b, c, d = np.ravel(aa), np.ravel(bb), np.ravel(cc), np.ravel(dd)
            bounds = array([a, b, c, d]).T
            where trange, lon_range, ... = [t_min, t_max], [lon_min, lon_max]
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
        kind=forge.FParameter.POSITIONAL_OR_KEYWORD,
        default=bounds)
    return param_xvec


def Functionalize_Dataset(kamodo_object, coord_dict, variable_name,
                          data_dict, gridded_int, coord_str, interp_flag=0,
                          func=None, times_dict=None, func_default='data'):
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
            If interp_flag=0 is used, data_array should be a numpy array of
            shape (c1, c2, c3, ..., cN). Otherwise, it could be a string or
            other thing that feeds into the given func.
        Note: The dataset must depend upon ALL of the coordinate arrays given.

        gridded_int: True to create a gridded interpolator (necessary for
            plotting and slicing). False otherwise (e.g. for the flythrough).
        coord_str: a string indicating the coordinate system of the data
            (e.g. "SMcar" or "GEOsph").
        interp_flag: the chosen method of interpolation.
            Options:
                0: (default option) assumes the given data_dict['data'] is an
                    N-dimensional numpy array and creates a standard scipy
                    interpolator to functionalize the data. This option should
                    be used when the entire dataset is contained in a single
                    array AND can easily fit into memory (e.g. 1D time series)
                1: Lazy interpolation is applied on top of the standard scipy
                    interpolators (one per time slice). This option is designed
                    to be used when the data is stored in one time step per
                    file and the time step can easily fit into memory.
                2: Lazy chunked interpolation is applied on top of the
                    standard scipy interpolators (one per time chunk).
                    Interpolation between time chunks occurs in the given
                    func. This option is designed to be used when the data is
                    stored with more than one time step per file AND the time
                    chunks are small enough to easily fit into memory.
                3. Combination of 1 and 2, used when the time chunks are too
                    large for the typical computer memory. Searches for the
                    correct time chunk, then the correct time slice in the
                    time chunk. This option is designed to be used when the
                    files contain more than one time step and are too large to
                    easily read into memory.
        func: a function defining the logic to be executed on a given time
            slice or chunk (e.g. converting to an array and then transposing
            it).
            *** Only needed for interp_flag greater than zero. ***
            - For interp_flag=1, the function should return only the data for
                the time slice. Required syntax is data = func(i), where i is
                the index of the slice on the full time grid.
            - For interp_flag=2, the function should return the data for the
                time chunk. Required syntax is data = func(i), where i is the
                file number in a list of files of the same naming pattern.
            - For interp_flag=3, the function should return only the data for
                the time slice. Required syntax is data = func(i, fi), where i
                is the file number in a list of files of the same naming
                pattern, and fi is the index of the time slice in that file.
            - ALL functions must include the method to retrieve the data from
                the file.
            - This option is useful when a dataset's shape order does not match
            the required order (e.g. t, lat, lon -> t, lon, lat) and allows
            such operations to be done on the fly.
        start_times: a list of the start times for each data chunk.
            *** Only needed for interpolation over time chunks. ***
            (interp_flag options 2 and 3).
        func_default: a string indicating the type of object returned by func.
            The default is 'data', indicating that the returned object is a
            numpy array. Set this to 'custom' to indicate that func returns
            a custom interpolator.

    Output: A kamodo object with the functionalized dataset added.
    '''

    # split the coord_dict into data and units
    coord_data = {key: value['data'] for key, value in coord_dict.items()}
    coord_units = {key: value['units'] for key, value in coord_dict.items()}

    # create the bounds array (see docstring of create_funcsig)
    # positions of corners of coordinate box
    data_list = [array([value['data'].min(), value['data'].max()]) for
                 key, value in coord_dict.items()]
    mesh_list = meshgrid(*data_list)
    bounds = array([ravel(item) for item in mesh_list], dtype=float).T

    # create a functionalized interpolator and modified function signature
    param_xvec = create_funcsig(coord_data, coord_str, bounds)
    if interp_flag == 0:  # standard logic
        interp = create_interp(coord_data, data_dict, func,
                               func_default=func_default)
    elif interp_flag == 1:
        interp = time_interp(coord_dict, data_dict, func,
                             func_default=func_default)
    elif interp_flag == 2:
        interp = multitime_interp(coord_dict, data_dict, times_dict, func,
                                  func_default=func_default)
    elif interp_flag == 3:
        interp = multitime_biginterp(coord_dict, data_dict, times_dict, func,
                                     func_default=func_default)
    new_interp = forge.replace('xvec', param_xvec)(interp)
    interp = kamodofy(units=data_dict['units'], data=data_dict['data'],
                      arg_units=coord_units)(new_interp)

    # add gridded version if requested, even for 1D functions
    kamodo_object[variable_name] = interp
    if gridded_int:
        interp_grid = kamodofy(gridify(interp, **coord_data),
                               units=data_dict['units'],
                               data=data_dict['data'], arg_units=coord_units)
        kamodo_object[variable_name+'_ijk'] = interp_grid
    return kamodo_object


def time_interp(coord_dict, data_dict, func, func_default='data'):
    '''Create a functionalized interpolator by splitting into timesteps.
    interp_flag=1 uses this function.
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
            slice i, which must include data retrieval from the file. The
            should return the data for the time slice. Required syntax is
            data = func(i).
        func_default: a string indicating the type of object returned by func.
            The default is 'data', indicating that the returned object is a
            numpy array. Set this to 'custom' to indicate that func returns
            a custom interpolator.
    Output: A lazy time interpolator.
    '''

    # create a list of interpolators per timestep, not storing the data, and
    # interpolate from there. such as done in superdarnea_interp.py
    # Assumes that time is the first dimension

    # split the coord_dict into data and units
    coord_list = [value['data'] for key, value in coord_dict.items()]  # arrays
    idx_map, time_interps = [], []  # initialize variables

    # define method for how to add time slices to memory
    def add_timeinterp(i, time_interps):
        if len(coord_list) > 1 and func_default == 'data':
            time_interps.append(rgiND(coord_list[1:], func(i),
                                      bounds_error=False, fill_value=NaN))
        elif len(coord_list) == 1:  # has to be different for time series data
            def dummy_1D(spatial_position):
                return func(i)  # data only
            time_interps.append(dummy_1D[0])
        elif func_default == 'custom':
            time_interps.append(func(i))
        return time_interps

    def interp_i(*args, time_interps=time_interps, idx_map=idx_map):

        position = array([*args])  # time coordinate value must be first
        times = position[0]
        sposition = position[1:].T  # spatial coordinates (len(times), ncoords)
        if len(position.shape) == 1:
            out_vals = zeros(1) * NaN
        else:
            out_vals = zeros(times.shape[0]) * NaN
        # loop through time grid instead of input time array
        for i in range(len(coord_list[0])):
            # figure out what times given are in this file, if any
            if i < len(coord_list[0]) - 1:
                end_time, idx_list = coord_list[0][i+1], [i, i+1]
            else:
                end_time, idx_list = coord_list[0][i], [i]
            if isinstance(times, ndarray):
                st_idx = [j for j, time_val in enumerate(times) if
                          time_val >= coord_list[0][i] and
                          time_val <= end_time]
            else:
                if times >= coord_list[0][i] and times <= end_time:
                    st_idx = [i]
                else:
                    continue
            if len(st_idx) == 0:  # skip if none
                continue

            # Interpolate values for chosen time grid values
            for i in idx_list:
                if i not in idx_map:
                    len_map = len(idx_map)  # length before adding another
                    try:  # allow for memory errors
                        time_interps = add_timeinterp(i, time_interps)
                    except MemoryError:  # remove one time slice first
                        print('Avoiding memory error...')
                        if len_map == 0:
                            print('Not enough memory to load one time slice.' +
                                  ' Please close some applications.')
                        elif len_map < 2:  # failed to add a second slice
                            print('Not enough memory to load two time slices' +
                                  '. Please close some applications.')
                        elif abs(i-idx_map[0]) > 2:
                            del idx_map[0], time_interps[0]
                        else:
                            del idx_map[-1], time_interps[-1]
                        time_interps = add_timeinterp(i, time_interps)
                    print(f'Time slice index {i} added from file.')
                    idx_map.append(i)
            if len(idx_list) > 1:
                interp_locations = [idx_map.index(val) for val in idx_list]
                if isinstance(times, ndarray):
                    interp_values = array([time_interps[i](sposition[st_idx])
                                           for i in interp_locations]).T
                    for j, vals in enumerate(interp_values):  # loop positions
                        time_int = rgi1D(coord_list[0][idx_list], vals,
                                         bounds_error=False,
                                         fill_value=NaN)
                        out_vals[st_idx[j]] = time_int(times[st_idx[j]])
                else:
                    interp_values = array([time_interps[i](sposition)
                                           for i in interp_locations]).T
                    time_int = rgi1D(coord_list[0][idx_list],
                                     interp_values,
                                     bounds_error=False, fill_value=NaN)
                    out_vals = time_int(times)
            else:
                interp_location = idx_map.index(idx_list[0])
                if isinstance(times, ndarray):
                    out_vals[st_idx] = time_interps[interp_location](
                        sposition[st_idx])
                else:
                    out_vals = time_interps[interp_location](sposition)
        return out_vals

    def interp(xvec):
        return interp_i(*array(xvec).T)
    return interp


def multitime_interp(coord_dict, data_dict, times_dict, func,
                     func_default='data'):
    '''Create a functionalized interpolator by splitting into time chunks.
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
            Interpolation between time chunks should be taken care of by
            appending the first time slice of the next time chunk to the end of
            the current time chunk (see func description).
        func: a function defining the logic to be executed on a given time
            chunk, including retrieving the data from the file. The function
            must return the data for the time chunk. Interpolation between time
            chunks can be easily achieved by appending a time slice from the
            next chunk to the current one in the logic of the function. The
            required syntax is data = func(i), where i is the file number.
        func_default: a string indicating the type of object returned by func.
            The default is 'data', indicating that the returned object is a
            numpy array. Set this to 'custom' to indicate that func returns
            a custom interpolator.
    Output: A time-chunked lazy interpolator.
    '''

    # create a list of interpolators per time chunk, not storing the data, and
    # interpolate from there. such as done in superdarnea_interp.py
    # Assumes that time is the first dimension.
    # Assumes that data_dict['data'][i] is a string used to find the right file

    # split the coord_dict into data and units, initialize variables/func
    coord_list = [value['data'] for key, value in coord_dict.items()]  # arrays
    idx_map, time_interps = [], []  # initialize variables

    # define method for how to add time chunks to memory
    def add_timeinterp(i, time_interps):  # i is the file number
        idx_map.append(i)
        # determine time grid for file
        start_idx = list(coord_list[0]).index(times_dict['start'][i])
        end_idx = list(coord_list[0]).index(times_dict['end'][i]) + 1
        if i < len(times_dict['end'])-1:  # interpolating btwn files
            end_idx += 1
        time = coord_list[0][start_idx:end_idx]
        # get data for chunk
        data = func(i)
        if len(coord_list) > 1 and func_default == 'data':
            coord_list_i = [time] + coord_list[1:]
            time_interps.append(rgiND(coord_list_i, data,
                                      bounds_error=False, fill_value=NaN))
        elif len(coord_list) == 1:  # has to be different for time series data
            time_interps.append(rgi1D(time, data, bounds_error=False,
                                      fill_value=NaN))
        elif func_default == 'custom':  # when func returns an interpolator
            time_interps.append(func(i))
        return time_interps, idx_map

    def interp_i(*args, time_interps=time_interps, idx_map=idx_map):
        position = array([*args])  # time coordinate value must be first
        if len(position.shape) == 1:
            out_vals = zeros(1) * NaN
        else:
            out_vals = zeros(position.shape[1]) * NaN  # (num_cgrids, num_pos)
        # loop through start times instead
        for i in range(len(times_dict['start'])):
            # figure out what times given are in this file, if any
            if i < len(times_dict['start']) - 1:
                end_time = times_dict['start'][i+1]
            else:
                end_time = times_dict['end'][i]
            if isinstance(position[0], ndarray):
                st_idx = [j for j, time_val in enumerate(position[0]) if
                          time_val >= times_dict['start'][i] and
                          time_val <= end_time]
            else:
                if position[0] >= times_dict['start'][i] and \
                        position[0] <= end_time:
                    st_idx = [i]
                else:
                    continue
            if len(st_idx) == 0:  # skip if none
                continue

            # Interpolate values for chosen time chunk
            if i not in idx_map:
                len_map = len(idx_map)
                try:  # allow for memory errors
                    time_interps, idx_map = add_timeinterp(i, time_interps)
                except MemoryError:  # remove two(?) items first
                    print('Avoiding memory error...')
                    if len_map == 0:  # tried but failed
                        print('Not enough memory to load a time chunk. ' +
                              'Please close some applications.')
                    elif i != idx_map[0]:
                        del idx_map[0], time_interps[0]
                    else:
                        del idx_map[-1], time_interps[-1]
                    time_interps, idx_map = add_timeinterp(i, time_interps)
                    print(f'Time chunk added from file {i+1}.')
            # interpolate all the positions with times in this time chunk
            if isinstance(position[0], ndarray):
                out_vals[st_idx] = time_interps[idx_map.index(i)](
                    position[:, st_idx].T)
            else:
                out_vals = time_interps[idx_map.index(i)](position)
        return out_vals

    def interp(xvec):
        return interp_i(*array(xvec).T)

    return interp


def multitime_biginterp(coord_dict, data_dict, times_dict, func,
                        func_default='data'):
    '''Create a functionalized interpolator by splitting into time chunks as
    stored in the files and then into time steps. Use this option (3) when the
    time chunks are large compared to the typically available computer memory.
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
            Interpolation between time chunks should be taken care of by
            appending the first time slice of the next time chunk to the end of
            the current time chunk (see func description).
        func: a function defining the logic to be executed on a given time
            chunk, including retrieving the correct time slice of data from the
            file. The function must return the data for the time SLICE.
            Required syntax is data = func(f, i), where f is the file number
            of the time chunk and i is the time slice number in that file.
            Counting starts at zero.
        func_default: a string indicating the type of object returned by func.
            The default is 'data', indicating that the returned object is a
            numpy array. Set this to 'custom' to indicate that func returns
            a custom interpolator.
    Output: A time-chunked lazy interpolator that loads two slices each time.
    '''

    # create a list of interpolators per time chunk, not storing the data, and
    # interpolate from there. such as done in superdarnea_interp.py
    # Assumes that time is the first dimension.
    # Assumes that data_dict['data'][i] is a string used to find the right file

    # split the coord_dict into data and units
    coord_list = [value['data'] for key, value in coord_dict.items()]  # arrays
    idx_map, time_interps = [], []  # initialize variables

    # define method for how to add time slices to memory
    # i is the file#, fi is the time slice # in file
    def add_timeinterp(i, fi, time_interps):
        if len(coord_list) > 1 and func_default == 'data':
            time_interps.append(rgiND(coord_list[1:], func(i, fi),
                                      bounds_error=False, fill_value=NaN))
        # has to be different for time series data
        elif len(coord_list) == 1 and func_default == 'data':
            def dummy_1D(spatial_position):
                return func(i, fi)  # data only
            time_interps.append(dummy_1D[0])
        elif func_default == 'custom':
            time_interps.append(func(i, fi))
        return time_interps

    def interp_i(*args, time_interps=time_interps, idx_map=idx_map):

        position = array([*args])  # time coordinate value must be first
        times = position[0]
        sposition = position[1:].T  # spatial coordinates (len(times), ncoords)
        if len(position.shape) == 1:
            out_vals = zeros(1) * NaN
        else:
            out_vals = zeros(times.shape[0]) * NaN
        # loop through time grid instead of input time array
        for ti in range(len(coord_list[0])):
            # figure out what times given are in this slice, if any
            if ti < len(coord_list[0]) - 1:
                end_time, idx_list = coord_list[0][ti+1], [ti, ti+1]
            else:
                end_time, idx_list = coord_list[0][ti], [ti]
            if isinstance(times, ndarray):
                st_idx = [j for j, time_val in enumerate(times) if
                          time_val >= coord_list[0][ti] and
                          time_val <= end_time]
            else:
                if times >= coord_list[0][ti] and times <= end_time:
                    st_idx = [ti]
                else:
                    continue
            if len(st_idx) == 0:  # skip if none
                continue

            # Interpolate values for chosen time grid values
            for idx in idx_list:
                if idx not in idx_map:
                    # Get index of file containing the desired time
                    i = get_file_index(coord_list[0][idx], times_dict['start'])
                    # determine which slice in the file corresponds to time val
                    # aka time position relative to beginning of file
                    start_idx = list(coord_list[0]).index(
                        times_dict['start'][i])
                    fi = list(coord_list[0][start_idx:]).index(
                        coord_list[0][idx])

                    # add time slice
                    len_map = len(idx_map)  # length before adding another
                    try:  # allow for memory errors
                        time_interps = add_timeinterp(i, fi, time_interps)
                    except MemoryError:  # remove one time slice first
                        print('Avoiding memory error...')
                        if len_map == 0:
                            print('Not enough memory to load one time slice.' +
                                  ' Please close some applications.')
                        elif len_map < 2:  # failed to add a second slice
                            print('Not enough memory to load two time slices' +
                                  '. Please close some applications.')
                        elif abs(idx-idx_map[0]) > 2:
                            del idx_map[0], time_interps[0]
                            del idx_map[0], time_interps[0]
                        else:
                            del idx_map[-1], time_interps[-1]
                            del idx_map[-1], time_interps[-1]
                        time_interps = add_timeinterp(i, fi, time_interps)
                    idx_map.append(idx)
                    print(f'Time slice index {idx} (file time {fi}) added ' +
                          f'from file {i+1}.')
            if len(idx_list) > 1:
                interp_locations = [idx_map.index(val) for val in idx_list]
                if isinstance(times, ndarray):
                    interp_values = array([time_interps[ii](sposition[st_idx])
                                           for ii in interp_locations]).T
                    for j, vals in enumerate(interp_values):  # loop positions
                        time_int = rgi1D(coord_list[0][idx_list], vals,
                                         bounds_error=False,
                                         fill_value=NaN)
                        out_vals[st_idx[j]] = time_int(times[st_idx[j]])
                else:  # one position
                    interp_values = array([time_interps[ii](sposition)
                                           for ii in interp_locations]).T
                    time_int = rgi1D(coord_list[0][idx_list],
                                     interp_values, bounds_error=False,
                                     fill_value=NaN)
                    out_vals = time_int(times)
            else:
                interp_location = idx_map.index(idx_list[0])
                if isinstance(times, ndarray):
                    out_vals[st_idx] = time_interps[interp_location](
                        sposition[st_idx])
                else:
                    out_vals = time_interps[interp_location](sposition)
        return out_vals

    def interp(xvec):
        return interp_i(*array(xvec).T)
    return interp


def get_slice_idx(time_val, time_array):
    '''Figures out where the time_val is in time_array.
    Inputs:
        time_val: float value of time
        time_array: array of float time values
    Returns: A list of indices for one time value on either side OR a list
        containing the single index that corresponds exactly to the given
        time_val.
    '''
    if time_val not in time_array:
        idx = sorted(append(time_array, time_val)).index(time_val)
        if idx > 1 and idx < len(time_array)-1:  # middle somewhere
            idx_list = [idx-1, idx]
        elif idx <= 1:  # beginning
            idx_list = [0, 1]
        elif idx > 1 and idx >= len(time_array)-1:  # at end
            idx_list = [len(time_array)-2, len(time_array)-1]
    else:
        idx = list(time_array).index(time_val)
        idx_list = [idx]
    return idx_list


def get_file_index(time_val, time_grid):
    '''Figures out which file start time is directly before the given time
    value.
    Inputs:
        time_val: float value of desired coordinate
        time_grid: an array of coordinate grid values of the same coordinate
    Returns: index of the start time directly preceding the given time value.
    '''

    if time_val not in time_grid:
        idx = sorted(append(time_grid, time_val)
                     ).index(time_val)  # get file #
        if idx == 0:
            i = 0  # beg of first file
        else:
            i = idx - 1  # somewhere else
    else:
        i = list(time_grid).index(time_val)
    return i


# ############# BEGIN PRESSURE LEVEL INVERSION ROUTINE #######################

# -*- coding: utf-8 -*-
# The pressure level inversion routine used by many ITM models


def PLevelInterp(h_func, time, longitude, latitude, ilev, units, km_grid,
                 km_min_max):
    '''Custom interpolator functionto convert from altitude in km to pressure
    level.
    Parameters:
        h_func - the H_ilev gridded function/interpolator in km
        time - a 1D array with the grid values for time
        longitude - a 1D array with the grid values for longitude
        latitude - a 1D array with the grid values for latitude
        ilev - a 1D array with the grid values for pressure level
        units - a string indicating the units of the pressure level grid
        km_grid - a 1D array of the median altitude values for each pressure
            level
        km_min_max - a 2-element list of the absolute min and max of the
            altitude in km

    Output: Two interpolators. The first interpolator accepts time, lon, lat,
        and pressure level and returns time, lon, lat, and height. The second
        interpolator accepts the same arguments and returns only the height.
        Both interpolators are 'kamodofied'.

    This code changes a function from being dependent on pressure level to
    instead depending upon altitude. The model reader calling this function
    needs to first create a gridded height variable using typical methods and
    use Kamodo-core's unit conversion capability to convert the gridded height
    function into km. This routine creates a function that uses the slicing
    technique to find the height values at a given time, longitude and latitude
    for the default grid values of pressure level. These height values and
    pressure level grid values are then used as the X and Y input values to a
    1D interpolator (height is X, pressure level is Y) to invert the function.
    This interpolator is then evaluated at the requested height (or heights)
    associated with the same time, longitude and latitude and returns the
    pressure level(s). The code repeats this process for each unique
    (t, lon, lat) trio. We use the interp1d interpolator for the 1D inversion
    (from SciPy). Regular and gridded interpolators based on this logic are
    returned.
    '''

    def km_to_ilev(t, lon, lat, km):
        '''Inputs t, lon, lat, and km are arrays;
        the interpolated pressure levels are returned.'''
        # sort lats and lons per time value
        input_arr = array([t, lon, lat]).T  # array of all positions given
        if not isinstance(t, ndarray):  # only one time/position
            km_vals = h_func(**{'time': input_arr[0], 'lon': input_arr[1],
                                'lat': input_arr[2]})
            km_interp = rgi1D(km_vals, ilev, bounds_error=False,
                              fill_value=NaN)
            return km_interp(km)

        # remaining logic is for if there is more than one position given
        pos_arr = unique(input_arr, axis=0)  # only create interp for unique
        out_ilev = zeros(len(t)) * NaN  # (t, lon, lat) positions to save time
        for pos in pos_arr:
            pos_idx = [i for i, p in enumerate(input_arr) if all(p == pos)]
            km_vals = h_func(**{'time': pos[0], 'lon': pos[1], 'lat': pos[2]})
            km_interp = rgi1D(km_vals, ilev, bounds_error=False,
                              fill_value=NaN)
            out_ilev[pos_idx] = km_interp(km[pos_idx])  # interp for all km
        return out_ilev

    # Convert f(ilev) to f(km)
    def plevconvert(xvec):
        '''Interpolator to convert from height to pressure level.
        Input xvec is a tuple or list of arrays/values.
        Returns the 4D position as a tuple of four 1D arrarys to feed into
        the original function.'''
        print('Inverting the pressure level grid. Please wait...', end="")
        t, lon, lat, km = array(xvec).T
        out_ilev = km_to_ilev(t, lon, lat, km)
        data = array([t, lon, lat, out_ilev]).T
        print('done.')
        return data

    def plevconvert_ijk(xvec):
        '''Interpolator to convert from height to pressure level.
        Input xvec is a tuple or list of arrays/values.
        Returns the 4D position as a tuple of four 1D arrarys to feed into
        the original function.'''
        print('Inverting the pressure level grid. Please wait...', end="")
        t, lon, lat, km = array(xvec).T
        data = km_to_ilev(t, lon, lat, km)
        print('done.')
        return data

    # split the coord_dict into data and units
    km_grid_altered = km_grid.copy()
    if diff(km_grid_altered).min() < 1:  # pressure grid inverted
        km_grid_altered[-1], km_grid_altered[0] = km_min_max
    else:
        km_grid_altered[0], km_grid_altered[-1] = km_min_max
    coord_data = {'time': time, 'lon': longitude, 'lat': latitude,
                  'height': km_grid_altered}
    coord_units = {'time': 'hr', 'lon': 'deg', 'lat': 'deg', 'height': 'km'}

    # create the bounds array (see docstring of create_funcsig)
    # positions of corners of coordinate box
    data_list = [array([value.min(), value.max()]) for key, value in
                 coord_data.items()]
    mesh_list = meshgrid(*data_list)
    bounds = array([ravel(item) for item in mesh_list], dtype=float).T

    # create the functionalized interpolators and modified function signatures
    fake_data = zeros((2, 2, 2, 2)) * NaN  # avoiding computation
    param_xvec = create_funcsig(coord_data, 'GDZsphkm', bounds)
    new_interp = forge.replace('xvec', param_xvec)(plevconvert)
    interp = kamodofy(units=units, data=fake_data, arg_units=coord_units
                      )(new_interp)
    new_interp = forge.replace('xvec', param_xvec)(plevconvert_ijk)
    interp_ijk = kamodofy(units=units, data=fake_data, arg_units=coord_units
                          )(new_interp)

    return interp, interp_ijk


def register_griddedPlev(kamodo_object, new_varname, units, interp_ijk,
                         coord_dict, kms, km_min_max):
    '''Properly register the gridded interpolators defined either from the
    pressure level inversion or by function composition. Need height coord
    grid to be median values in gridded creation but have absolute max/min in
    non-gridded creation.'''

    new_coord_units = {'time': 'hr', 'lon': 'deg',
                       'lat': 'deg', 'height': 'km'}
    fake_data = zeros((2, 2, 2, 2)) * NaN  # avoiding computation
    coord_data = {key: value['data'] for key, value in
                  coord_dict.items() if key in
                  new_coord_units.keys()}  # exclude ilev
    coord_data['height'] = insert(insert(kms, len(kms), km_min_max[-1]),
                                  0, km_min_max[0])  # add min and max
    kamodo_object.variables[new_varname+'_ijk'] = {'data': fake_data,
                                                   'units': units}
    kamodo_object[new_varname+'_ijk'] = kamodofy(
        gridify(interp_ijk, **coord_data),
        units=kamodo_object.variables[new_varname+'_ijk']['units'],
        data=kamodo_object.variables[new_varname+'_ijk']['data'],
        arg_units=new_coord_units)
    return kamodo_object


# ################ BEGIN TIME LIST FUNCTIONS #########################


@vectorize
def hrs_to_str(hrs, filedate):
    '''Convert hrs since midnight of first day to a string for the file list
    of format "Date: YYYY-MM-DD  Time: HH:MM:SS".'''
    h, m = floor(hrs), round(hrs % 1 * 60., 4)
    sec = round(((hrs - h) * 60. - m) * 60.)
    return datetime.strftime(filedate + timedelta(hours=int(h), minutes=int(m),
                                                  seconds=int(sec)),
                             '  Date: %Y-%m-%d  Time: %H:%M:%S')


@vectorize
def str_to_hrs(dt_str, filedate, format_string='%Y-%m-%d %H:%M:%S'):
    '''Convert datetime string of format "YYYY-MM-DD HH:MM:SS" to hrs since
    midnight of filedate.'''
    tmp = datetime.strptime(dt_str, format_string).replace(
        tzinfo=timezone.utc)
    return float32((tmp - filedate).total_seconds()/3600.)


@vectorize
def hrs_to_tstr(hrs, ms_timing=False):
    '''Convert number of hours into HH:MM:SS format. If ms_timing is True,
    str is in format HH:MM:SS.mms'''
    h, m = floor(hrs), round(hrs % 1 * 60., 4)
    if not ms_timing:
        sec = round(((hrs - h) * 60. - m) * 60.)
        return str(f'\n{int(h):02d}:{int(m):02d}:{int(sec):02d}')
    else:
        sec = ((hrs - h) * 60. - int(m)) * 60.
        sec_str = f'{sec:2.3f}'.zfill(6)
        return str(f'\n{int(h):02d}:{int(m):02d}:'+sec_str)


@vectorize
def tstr_to_hrs(time_str, ms_timing=False):
    '''Convert str from HH:MM:SS format to float 32. If ms_timing is True,
    str is in format HH:MM:SS.mms'''
    hh, mm, ss = time_str.split(':')
    t = float32(hh) + float32(mm)/60. + float32(ss)/3600.
    return float32(t)


def create_timelist(list_file, time_file, modelname, times, pattern_files,
                    filedate, ms_timing=False):
    '''Used by all the readers to create the time_file and list_file files.
    Inputs:
        list_file - the name of the file, including the file directory, to
            write the list of all of the files, start dates, start times,
            end dates, and end times.
        time_file - the name of the file, including the file directory, to
            write the time grid for each file pattern.
        modelname - a string with the name of the model
        times - a dictionary with keys indicating the file patterns.
            ['all'] - the complete time grid across the entire file_dir for
                the indicated file pattern.
            ['start'] - the starting time for each of the files of the given
                file pattern.
            ['end'] - the end time for each of the files of the given file
                pattern.
            All times are stored as the number of hours since midnight of the
                first file of all the files in the given directory.
        pattern_files - a dictionary with keys indicating the file patterns.
            Each key has as its value the list of files in the current file
            directory that match the file pattern.
        filedate - a datetime object indicating the date of the first file of
            all the files in the file_dir.
        ms_timing - a boolean indicating whether millisecond timing is needed.
            Default is False.
    Returns nothing.
    '''

    # create time list file if DNE
    print('Creating the time files...', end="")
    list_out = open(list_file, 'w')
    list_out.write(f'{modelname} file list start and end ' +
                   'dates and times')
    time_out = open(time_file, 'w')
    time_out.write(f'{modelname} time grid per pattern')
    for p in pattern_files.keys():
        # print out time grid to time file
        time_out.write('\nPattern: '+p)
        str_out = hrs_to_tstr(times[p]['all'], ms_timing)
        time_out.write(''.join(str_out))
        # print start and end dates and times to list file
        start_time_str = hrs_to_str(times[p]['start'], filedate)
        end_time_str = hrs_to_str(times[p]['end'], filedate)
        files = pattern_files[p]
        for i in range(len(files)):
            list_out.write('\n' + files[i].replace('\\', '/') +
                           start_time_str[i] + '  ' + end_time_str[i])
    time_out.close()
    list_out.close()
    print('done.')
    return


def read_timelist(time_file, list_file, ms_timing=False):
    '''Used by all readers to read in the time grid from the time_file and
    the list of files with start and end times from the list_file.
    Returns:
        times - a dictionary with keys indicating the file patterns.
            ['all'] - the complete time grid across the entire file_dir for
                the indicated file pattern.
            ['start'] - the starting time for each of the files of the given
                file pattern.
            ['end'] - the end time for each of the files of the given file
                pattern.
            All times are stored as the number of hours since midnight of the
                first file of all the files in the given directory.
        pattern_files - a dictionary with keys indicating the file patterns.
            Each key has as its value the list of files in the current file
            directory that match the file pattern.
        filedate - a datetime object indicating the date of the first file of
            all the files in the file_dir.
        filename - a string containing the names of all the files in the file
            directory, separated by commas.
        ms_timing - a boolean indicating whether millisecond timing is needed.
            Default is False.
    '''

    # get time grids and initialize self.times structure
    times, pattern_files = {}, {}
    time_obj = open(time_file)
    data = time_obj.readlines()
    for line in data[1:]:
        if 'Pattern' in line:
            p = line.strip()[9:]
            times[p] = {'start': [], 'end': [], 'all': []}
        else:
            times[p]['all'].append(line.strip())
    time_obj.close()
    for p in times.keys():
        times[p]['all'] = tstr_to_hrs(times[p]['all'], ms_timing)

    # get filenames, dates and times from list file
    files, start_date_times, end_date_times = [], [], []
    list_obj = open(list_file)
    data = list_obj.readlines()
    for line in data[1:]:
        file, tmp, date_start, tmp, start_time, tmp, date_end, \
            tmp, end_time = line.strip().split()
        files.append(file)
        start_date_times.append(date_start+' '+start_time)
        end_date_times.append(date_end+' '+end_time)
    list_obj.close()
    filedate = datetime.strptime(start_date_times[0][:10], '%Y-%m-%d').replace(
        tzinfo=timezone.utc)
    for p in times.keys():
        pattern_files[p] = [f for f in files if p in basename(f)]
        start_times = [t for f, t in zip(files, start_date_times)
                       if p in f]
        times[p]['start'] = str_to_hrs(start_times, filedate)
        end_times = [t for f, t in zip(files, end_date_times)
                     if p in f]
        times[p]['end'] = str_to_hrs(end_times, filedate)
    filename = ''.join([f+',' for f in files])[:-1]

    return times, pattern_files, filedate, filename
