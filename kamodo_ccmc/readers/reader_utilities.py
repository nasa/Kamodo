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
from time import perf_counter


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
                          func=None, start_times=None):
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
                1: Lazy interpolation is applied on top of the standard scipy
                    interpolators (one per time slice).
                2: Lazy chunked interpolation is applied on top of the
                    standard scipy interpolators (one per time chunk).
                    Interpolation between time chunks occurs in the given
                    func.
                3. Combination of 2 and 3, used when the time chunks are too
                    large for the typical computer memory. Searches for the
                    correct time chunk, then the correct time slice in the
                    time chunk.
        func: a function defining the logic to be executed on a given time
            slice or chunk (e.g. converting to an array and then transposing
            it). 
            *** Only needed for interp_flag greater than zero. ***
            - For interp_flag=1, the function should return only the data for
                the time slice. Required syntax is data = func(i), where i is
                the index of the slice on the full time grid.
            - For interp_flag=2, the function should return the data for the
                time chunk AND the time grid for that time chunk. Required
                syntax is data, time = func(i), where i is the file number in
                a list of files of the same naming pattern.
            - For interp_flag=3, the function should return only the data for
                the time slice. Required syntax is data = func(i), where i is
                a string of the combined indices, such as '1_3' indicating
                the second time chunk and the fourth time slice in that time
                chunk (counting starts at 0 for both).
            - ALL functions must include the method to retrieve the data from
                the file.
            - This option is useful when a dataset's shape order does not match
            the required order (e.g. t, lat, lon -> t, lon, lat) and allows
            such operations to be done on the fly.
        start_times: a list of the start times for each data chunk.
            *** Only needed for interpolation over time chunks. ***

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
        interp = multitime_interp(coord_dict, data_dict, start_times, func)
    elif interp_flag == 3:
        interp = multitime_biginterp(coord_dict, data_dict, start_times, func)
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


def time_interp(coord_dict, data_dict, func):
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
            slice i, which must include data retrieval from the file. The
            should return the data for the time slice. Required syntax is
            data = func(i).
    Output: A lazy time interpolator. 
    '''

    # create a list of interpolators per timestep, not storing the data, and 
    # interpolate from there. such as done in superdarnea_interp.py
    # Assumes that time is the first dimension

    # split the coord_dict into data and units
    coord_list = [value['data'] for key, value in coord_dict.items()]  # list of arrays
    idx_map, time_interps = [], []  # initialize variables

    # define method for how to add time slices to memory
    def add_timeinterp(i, time_interps):
        if len(coord_list) > 1:
            time_interps.append(rgiND(coord_list[1:], func(i), 
                                      bounds_error=False, fill_value=NaN))
        else:  # has to be different for time series data
            def dummy_1D(spatial_position):
                return func(i)  # data only
            time_interps.append(dummy_1D[0])
        return time_interps

    @vectorize
    def interp_i(*args, time_interps=time_interps, idx_map=idx_map):

        position = [*args]  # time coordinate value must be first
        # Choose indices of time grid values surrounding desired time
        idx_list = get_slice_idx(position[0], coord_list[0])

        # Interpolate values for chosen time grid values
        for i in idx_list:
            if i not in idx_map:
                len_map = len(idx_map)  # length before adding another
                try:  # allow for memory errors
                    time_interps = add_timeinterp(i, time_interps)
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
                    time_interps = add_timeinterp(i, time_interps)
                print(f'Time slice index {i} added from file.')
                idx_map.append(i)
        if len(idx_list) > 1:
            interp_locations = [idx_map.index(val) for val in idx_list]
            interp_values = ravel(array([time_interps[i]([*position[1:]])
                                         for i in interp_locations]))
            time_interps = rgi1D(coord_list[0][idx_list], interp_values,
                                bounds_error=False, fill_value=NaN)
            return time_interps(position[0])
        else: 
            interp_location = idx_map.index(idx_list[0])
            return time_interps[interp_location](position[1:])  # single time

    def interp(xvec):
        return interp_i(*array(xvec).T)
    return interp


def multitime_interp(coord_dict, data_dict, start_times, func):
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
            Interpolation between time chunks should be taken care of by appending
            the first time slice of the next time chunk to the end of the
            current time chunk (see func description).
        func: a function defining the logic to be executed on a given time
            chunk, including retrieving the data from the file. The function
            must return the data for the time chunk AND the time grid for
            that time chunk. Interpolation between time chunks can be easily
            achieved by appending a time slice from the next chunk to the
            current one in the logic of the function. Required syntax is
            data, time = func(i).
    Output: A time-chunked lazy interpolator.
    '''
    
    # create a list of interpolators per time chunk, not storing the data, and 
    # interpolate from there. such as done in superdarnea_interp.py
    # Assumes that time is the first dimension.
    # Assumes that data_dict['data'][i] is a string used to find the right file

    # split the coord_dict into data and units, initialize variables/func
    coord_list = [value['data'] for key, value in coord_dict.items()]  # list of arrays
    idx_map, time_interps = [], []  # initialize variables

    # define method for how to add time chunks to memory
    def add_timeinterp(i):  # i is the file number
        idx_map.append(i)
        # determine time grid for file
        start_idx = list(coord_list[0]).index(start_times[i])
        if i == len(start_times)-1:
            time = coord_list[0][start_idx:]
        else:  # determine the first index of the next file
            end_idx = list(coord_list[0]).index(start_times[i+1])
            time = coord_list[0][start_idx:end_idx]
        # get data for chunk
        data = func(i)
        if len(coord_list) > 1:
            coord_list_i = [time] + coord_list[1:]
            time_interps.append(rgiND(coord_list_i, data,
                                      bounds_error=False, fill_value=NaN))
        else:  # has to be different for time series data
            time_interps.append(rgi1D(time, data, bounds_error=False,
                                      fill_value=NaN))
        return time_interps, idx_map
    
    @vectorize
    def interp_i(*args, time_interps=time_interps, idx_map=idx_map):

        position = [*args]  # time coordinate value must be first
        # get location of file
        i = get_file_index(position[0], start_times)

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
        print(f'Time chunk index {i} added from file.')
        interp_location = idx_map.index(i)
        return time_interps[interp_location](position)  # single time chunk

    def interp(xvec):
        return interp_i(*array(xvec).T)

    return interp


def multitime_biginterp(coord_dict, data_dict, start_times, func):
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
            Interpolation between time chunks should be taken care of by appending
            the first time slice of the next time chunk to the end of the
            current time chunk (see func description).
        func: a function defining the logic to be executed on a given time
            chunk, including retrieving the correct time slice of data from the
            file. The function must return the data for the time SLICE.
            Required syntax is data = func(f, i), where f is the file number
            of the time chunk and i is the time slice number in that file.
            Counting starts at zero.
    Output: A time-chunked lazy interpolator that loads two slices each time.
    '''
    
    # create a list of interpolators per time chunk, not storing the data, and 
    # interpolate from there. such as done in superdarnea_interp.py
    # Assumes that time is the first dimension.
    # Assumes that data_dict['data'][i] is a string used to find the right file

    # split the coord_dict into data and units
    coord_list = [value['data'] for key, value in coord_dict.items()]  # list of arrays
    idx_map, time_interps = [], []  # initialize variables

    # define method for how to add time slices to memory
    def add_timeinterp(i, fi, time_interps):  # i is the file#, fi is the time slice # in file
        if len(coord_list) > 1:
            time_interps.append(rgiND(coord_list[1:], func(i, fi), 
                                      bounds_error=False, fill_value=NaN))
        else:  # has to be different for time series data
            def dummy_1D(spatial_position):
                return func(i, fi)  # data only
            time_interps.append(dummy_1D[0])
        return time_interps

    @vectorize
    def interp_i(*args, time_interps=time_interps, idx_map=idx_map):

        position = [*args]  # time coordinate value must be first
        # Choose indices of time grid values surrounding or equal to time
        idx_list = get_slice_idx(position[0], coord_list[0])

        # add time slices as needed
        for idx in idx_list:
            if idx not in idx_map:
                # Get index of file containing the desired time
                i = get_file_index(coord_list[0][idx], start_times)
    
                # determine which slice in the file corresponds to the time val
                # aka time position relative to beginning of file
                start_idx = list(coord_list[0]).index(start_times[i])
                fi = list(coord_list[0][start_idx:]
                                ).index(coord_list[0][idx])

                # add index
                len_map = len(idx_map)  # length before adding another
                try:  # allow for memory errors
                    time_interps = add_timeinterp(i, fi, time_interps)
                except MemoryError:  # remove one time slice first
                    print('Avoiding memory error...')
                    if len_map == 0:
                        print('Not enough memory to load one time slice. ' +
                              'Please close some applications.')
                    elif len_map < 2:  # tried to add a second slice but failed
                        print('Not enough memory to load two time slices. ' +
                              'Please close some applications.')
                    elif abs(idx-idx_map[0]) > 2:
                        del idx_map[0], time_interps[0]
                        del idx_map[0], time_interps[0]
                    else:
                        del idx_map[-1], time_interps[-1]
                        del idx_map[-1], time_interps[-1]
                    time_interps = add_timeinterp(i, fi, time_interps)
                idx_map.append(idx)
                print(f'Time slice index {idx} added from file.')
        if len(idx_list) > 1:
            interp_locations = [idx_map.index(val) for val in idx_list]
            interp_values = ravel(array([time_interps[i]([*position[1:]])
                                         for i in interp_locations]))
            time_interps = rgi1D(coord_list[0][idx_list], interp_values,
                                bounds_error=False, fill_value=NaN)
            return time_interps(position[0])
        else: 
            interp_location = idx_map.index(idx_list[0])
            return time_interps[interp_location](position[1:])  # single time

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


def get_file_index(time_val, start_times):
    '''Figures out which file start time is directly before the given time
    value.
    Inputs:
        time_val: float value of time
        start_times: an array of time values of the start of each chunk/file
    Returns: index of the start time directly preceding the given time value.
    '''
    
    if time_val not in start_times:
        idx = sorted(append(start_times, time_val)
                     ).index(time_val)  # get file #
        if idx == 0:
            i = 0  # beg of first file
        else:
            i = idx - 1  # somewhere else
    else:
        i = list(start_times).index(time_val)
    return i


################# BEGIN NON-INTERPOLATOR FUNCTIONS #########################
from datetime import datetime, timezone, timedelta


@vectorize
def hrs_to_str(hrs, filedate):
    '''Convert hrs since midnight of first day to a string for the file list
    of format "Date: YYYY-MM-DD  Time: HH:MM:SS".'''
    return datetime.strftime(filedate + timedelta(hours=hrs),
                             '  Date: %Y-%m-%d  Time: %H:%M:%S')


@vectorize
def str_to_hrs(dt_str, filedate):
    '''Convert datetime string of format "YYYY-MM-DD HH:MM:SS" to hrs since
    midnight of filedate.'''
    tmp = datetime.strptime(dt_str, '%Y-%m-%d %H:%M:%S').replace(
        tzinfo=timezone.utc)
    return (tmp - filedate).total_seconds()/3600.


def create_timelist(list_file, time_file, modelname, times, pattern_files,
                    filedate):
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
    Returns nothing.
    '''

    # create time list file if DNE
    list_out = open(list_file, 'w')
    list_out.write(f'{modelname} file list start and end ' +
                   'dates and times')
    time_out = open(time_file, 'w')
    time_out.write(f'{modelname} time grid per pattern')
    for p in pattern_files.keys():
        # print out time grid to time file
        time_out.write('\nPattern: '+p)
        for t in times[p]['all']:
            time_out.write('\n'+str(t))
        # print start and end dates and times to list file
        start_time_str = hrs_to_str(times[p]['start'], filedate)
        end_time_str = hrs_to_str(times[p]['end'], filedate)
        files = pattern_files[p]
        for i in range(len(files)):
            list_out.write('\n'+files[i].replace('\\', '/')+
                           start_time_str[i]+'  '+end_time_str[i])
    time_out.close()
    list_out.close()
    return

def read_timelist(time_file, list_file):
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
            directory, separated by commas.'''
    
    # get time grids and initialize self.times structure
    times, pattern_files = {}, {}
    time_obj = open(time_file)
    data = time_obj.readlines()
    for line in data[1:]:
        if 'Pattern' in line:
            p = line.strip()[9:]
            times[p] = {'start': [], 'end': [], 'all': []}
        else:
            times[p]['all'].append(float(line.strip()))
    time_obj.close()
    for p in times.keys():
        times[p]['all'] = array(times[p]['all'])
    
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
    filedate = datetime.strptime(start_date_times[0][:10],
                                      '%Y-%m-%d').replace(
                                          tzinfo=timezone.utc)
    for p in times.keys():
        pattern_files[p] = [f for f in files if p in f]
        start_times = [t for f, t in zip(files, start_date_times)
                       if p in f]
        times[p]['start'] = str_to_hrs(start_times, filedate)
        end_times = [t for f, t in zip(files, end_date_times)
                     if p.replace('_', '.') in f]
        times[p]['end'] = str_to_hrs(end_times, filedate)
    filename = ''.join([f+',' for f in files])[:-1]

    return times, pattern_files, filedate, filename

