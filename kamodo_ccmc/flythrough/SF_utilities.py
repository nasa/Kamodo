# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 17:09:14 2021
@author: rringuet
Most needed functions to support the SatelliteFlythrough software.
"""
from numpy import vectorize, array, diff, where, float64, unique, concatenate
from time import perf_counter
from os.path import basename, isfile
from datetime import datetime, timedelta, timezone
import kamodo_ccmc.flythrough.model_wrapper as MW
from kamodo_ccmc.flythrough.utils import ConvertCoord


@vectorize
def ts_to_hrs(time_val, filedate):
    '''Convert array of timestamps to hours since midnight of filedate string
    '''
    file_datetime = datetime.strptime(filedate+' 00:00:00',
                                      '%Y-%m-%d %H:%M:%S')
    return (datetime.utcfromtimestamp(time_val) -
            file_datetime).total_seconds()/3600.


@vectorize
def hrs_to_ts(time_val, filedate):
    '''Convert array of hours since midnight of filedate string to timestamps
    '''
    file_datetime = datetime.strptime(filedate + ' 00:00:00',
                                      '%Y-%m-%d %H:%M:%S').replace(
                                          tzinfo=timezone.utc)
    return datetime.timestamp(file_datetime + timedelta(hours=time_val))


def ts_to_ISOstring(utc_ts):
    '''Convert timestamp to string of format 2017-05-28T00:00:00'''
    return datetime.utcfromtimestamp(utc_ts).isoformat()


def check_plot_dir(plot_dir):
    '''If plot_dir does not exist, create it'''
    from os import path, mkdir
    if not path.isdir(plot_dir):
        mkdir(plot_dir)
    return


def write_timescsv(csv_filename, times):
    '''writes times dict from day_files to csv for faster execution time next
    time.'''
    data_out = open(csv_filename, 'w')
    data_out.write('# '+csv_filename)
    data_out.write('\n#file_date, filename, datetimes[0], datetimes[1], ' +
                   'filetimes[0], filetimes[1], dt')
    for key in times.keys():
        data_out.write('\n'+key+','+''.join([f'{value},' for value in
                                             times[key]]).strip(','))
    data_out.close()
    return


def read_timescsv(csv_filename):
    '''reads times dict from csv_filename for faster execution time.'''
    times = {}
    data_in = open(csv_filename, 'r')
    lines = data_in.readlines()
    data_in.close()
    for line in lines[2:]:
        data_line = line.strip('\n').split(',')
        times[data_line[0]] = data_line[1:4]
        times[data_line[0]].extend([float64(data_line[4]),
                                    float64(data_line[5]),
                                    float(data_line[6])])
    return times


def day_files(file_pattern, model, call_type):
    '''Retrieve file times. Convert files if necessary.'''
    # file_pattern could be a list if more than one pattern in file_dir exists
    if not isinstance(file_pattern, str):  # if a list/array of strings
        files, times = file_pattern, {}  # run reader with file_prefixes given
    else:
        from glob import glob
        files, times = sorted(glob(file_pattern)), {}

    # collect only time information from files for full time range
    reader = MW.Model_Reader(model)
    for f in files:
        k = reader(f, variables_requested=[], filetime=True, fulltime=True,
                   printfiles=False)
        if hasattr(k, 'conversion_test'):
            if not k.conversion_test:
                continue  # if file conversion errors, skip file_pattern
        if call_type == 'normal':
            file_date = k.datetimes[0][0:10]  # 'YYYY-MM-DD'
        elif call_type == 'single':
            file_date = k.datetimes[0][0:13].replace(' ', '_')  # YYYY-MM-DD_HH
        if file_date not in times.keys():  # prevent overwriting
            times[file_date] = [f, k.datetimes[0], k.datetimes[1],
                                k.filetimes[0], k.filetimes[1], k.dt]
        else:
            file_date += '_' + k.datetimes[0][11:13]  # 'YYYY-MM-DD_HH'
            times[file_date] = [f, k.datetimes[0], k.datetimes[1],
                                k.filetimes[0], k.filetimes[1], k.dt]
    return times


def check_timescsv(file_pattern, model, call_type='normal'):
    '''check for times csv file, write if not found in file_dir or if outdated
    '''
    # file_pattern could be a list if more than one pattern in file_dir exists
    if not isinstance(file_pattern, str):  # if a list/array of strings
        sample_pattern = file_pattern[0]
    else:
        sample_pattern = file_pattern

    # determine csv filename
    sample_prefix = basename(sample_pattern)
    file_dir = sample_pattern.split(sample_prefix)[0]
    if call_type == 'normal':
        csv_filename = file_dir + model + '_times.csv'
    elif call_type == 'single':
        csv_filename = file_dir + model + '_singletimes.csv'

    # if file DNE or outdated, write and return, else read and return
    if not isfile(csv_filename):
        times = day_files(file_pattern, model, call_type)
        write_timescsv(csv_filename, times)
    else:
        times = read_timescsv(csv_filename)

        # compare file contents to data in dir
        # file_patterns in times.csv file
        oldfile_pattern = [value[0] for key, value in times.items()]
        # if same length -> compare contents
        if len(oldfile_pattern) == len(file_pattern):
            compare = sum(oldfile_pattern == file_pattern)
            if compare == len(file_pattern):  # if match, then return
                return times
            else:  # not matching -> delete old file and recreate
                times = day_files(file_pattern, model, call_type)
                write_timescsv(csv_filename, times)
        else:  # different lengths -> delete old file and recreate
            times = day_files(file_pattern, model, call_type)
            write_timescsv(csv_filename, times)
    return times


def save_times(file_patterns, sat_time, model, verbose=False):
    '''Adjust times between files to filetime within half of dt in seconds
    (7.5min by default). Works for models with one day of data per file.'''
    times = check_timescsv(file_patterns, model)
    # look for sat_times not in files
    l_idx, file_dates = 0, list(times.keys())
    for i in range(len(file_dates)):  # filter out times not in files
        idx = where((sat_time >= times[file_dates[i]][3]) &
                    (sat_time <= times[file_dates[i]][4] +
                     times[file_dates[i]][5]))[0]  # end_time+dt
        times[file_dates[i]].append(idx)
        # remove indices from previous idx list if it occurs in this one
        if i > 0:
            # prefer time to be after beg time and not in dt section after end
            tmp_idx = [ival for ival in times[file_dates[i-1]][6] if ival not
                       in idx]
            times[file_dates[i-1]][6] = tmp_idx

    # collect indices into one array for plotting
    net_idx = array(concatenate(tuple([times[file_date][6] for file_date in
                                       times.keys()])), dtype=int)
    l_idx = len(net_idx)
    test_idx = unique(net_idx)
    if len(test_idx) != l_idx:
        print(l_idx, len(test_idx))
        raise AttributeError("net_idx has repeating values. Idx filtering " +
                             "didn't work correctly.")

    # print errors for any remaining 'bad' times
    nbad_times = len(sat_time) - l_idx
    if nbad_times == 1:
        print(f'{nbad_times} time is not in model output files and is ' +
              'excluded from the flythrough.')
    elif nbad_times > 0:
        print(f'{nbad_times} times are not in model output files and are ' +
              'excluded from the flythrough.')
    return sat_time, times, net_idx


def sat_tracks(sat_time, c1, c2, c3, z_dependencies, verbose=False):
    '''Calculate satellite tracks for interpolation'''

    # Create satellite tracks with appropriate inputs
    sat_track = {}  # initialize list of satellite tracks
    if '3D' in z_dependencies.keys():
        if verbose:
            print('Building height-independent satellite track.')
        sat_track['3D'] = [[t, c1_val, c2_val] for t, c1_val, c2_val in
                           zip(sat_time, c1, c2)]
    if '4D' in z_dependencies.keys():
        if verbose:
            print('Building height-dependent satellite track.')
        sat_track['4D'] = [[t, c1_val, c2_val, c3_val] for t, c1_val, c2_val,
                           c3_val in zip(sat_time, c1, c2, c3)]
    return sat_track


def Model_FlyAway(reader, filename, variable_list, sat_time, c1, c2, c3,
                  z_dependencies, verbose=False):
    '''Perform flythrough for one day of data and one coordinate system.'''

    # create kamodo object, initialize some variables
    var_list = variable_list.copy()  # save copy before it gets altered
    kamodo_object = reader(filename, variables_requested=variable_list,
                           gridded_int=False)

    # remove requested variables not found in data from var_list and z_dependen
    newvar_list = [var for var in var_list if var in
                   kamodo_object.variables.keys()]
    z_dependencies = {key: [var for var in value if var in newvar_list]
                      for key, value in z_dependencies.items()}
    # e.g. {'3D': a list of variables, '4D': a list of variables}
    # used as a mapping of coordinate types to variables

    # create satellite tracks of types needed based on vertical dependencies
    sat_track = sat_tracks(sat_time, c1, c2, c3, z_dependencies,
                           verbose=verbose)

    # retrieve interpolator and interpolate data for each variable, using track
    #   type appropriate for each variable.
    results = {var: kamodo_object[var](sat_track[[key for key, value in
                                                  z_dependencies.items() if var
                                                  in value][0]]) for var in
               newvar_list}
    del kamodo_object   # save memory
    return results


def coordinate_systems(model, sat_time, c1, c2, c3, variable_list, coord_type,
                       coord_grid):
    '''Detect what coordinate system is needed per variable, convert and return
    per type.
    Returns a dictionary. Keys are the coordinate systems. Values are
    0: a list of the variable names needing those coordinates
    1: x/lon 1D array
    2: y/lat 1D array
    3: z/rad/height 1D array
    4: z_dependencies dictionary
    '''

    # determine coordinate types needed
    # convert to alternative coordinates if necessary
    var_dict = MW.Model_Variables(model, return_dict=True)
    # {varname: [desc, int, coord_name, grid_type, coord_list, unit]}
    var_coord_strs = unique([value[2] + ',' + value[3] for
                             key, value in var_dict.items() if key in
                             variable_list])
    # 'SPH,sph','MAG,sph','GDZ,sph', etc
    if len(var_coord_strs) != 1 or var_coord_strs[0] != (coord_type + ',' +
                                                         coord_grid):
        # then coordinate conversion needed
        new_coords = {coord_name: [[key for key, value in var_dict.items()
                                    if (value[2].split('_')[0]+','+value[3] ==
                                        coord_name) and key in variable_list]]
                      for coord_name in var_coord_strs}  # key is 'name,type'
        # first value is a list of the variable names needing those coordinates
        for key in new_coords.keys():   # e.g. key= 'GDZ,sph'
            # convert to needed coordinates
            alt_c1, alt_c2, alt_c3, units_out = \
                    ConvertCoord(sat_time, c1, c2, c3, coord_type, coord_grid,
                                 *key.split(','))
            new_coords[key].extend([alt_c1, alt_c2, alt_c3])  # elements 1,2,3
        # second, third and fourth values are 1D position arrays
    else:
        new_coords = {coord_type+','+coord_grid: [variable_list, c1, c2, c3]}

    # determine z_dependency of relevant variables for each coordinate system
    for key in new_coords.keys():
        z_dependencies = {}
        coord_lengths = [len(value[4]) for keyv, value in var_dict.items() if
                         keyv in new_coords[key][0]]
        if 3 in coord_lengths:  # determine which variables are 3D
            z_dependencies['3D'] = [keyv for keyv, value in var_dict.items() if
                                    len(value[4]) == 3 and keyv in
                                    new_coords[key][0]]
        if 4 in coord_lengths:  # determine which variables are 4D
            z_dependencies['4D'] = [keyv for keyv, value in var_dict.items() if
                                    len(value[4]) == 4 and keyv in
                                    new_coords[key][0]]
        new_coords[key].append(z_dependencies)  # element 5
        # e.g. {'3D': a list of variables, '4D': a list of variables}
        # used as a mapping of coordinate types to variables
    return new_coords


def Model_SatelliteFlythrough(model, file_dir, variable_list, sat_time, c1, c2,
                              c3, coord_type, coord_grid, verbose=False):
    '''
    Execute flythrough for model data. Returns results_dict.
    results_dict is a dictionary of the interpolated data for the entire data
        set sorted by variable name.
    file_dir is a string indicating where the data files are located.
    variable_list is a list of strings of the desired variable names.
    sat_time is an array of timestamp values.
    c1, c2, c3 = x, y, z or lon, lat, height (or radius)
    if x, y, z, then must be in R_E units
    if radius -> R_E. if height -> km
    '''

    # Check that sat data is all the same length, will error if not
    if max(diff(array([len(sat_time), len(c3), len(c2), len(c1)]))) > 0:
        raise AttributeError('Satellite arrays or lists must all be the ' +
                             'same length. Current array lengths are ' +
                             f'{len(sat_time)}, {len(c1)}, {len(c2)}, and ' +
                             f'{len(c3)}')

    # reader prefers converted filename, even if it does not exist.
    # will create if no wrapped data found.
    file_patterns = MW.FileSearch(model, file_dir)
    reader = MW.Model_Reader(model)  # Kamodo gets imported here

    # match trajectory times to model data output files
    sat_time, times, net_idx = save_times(file_patterns, sat_time, model,
                                          verbose=verbose)

    # initialize results dictionary with given trajectory
    results_dict = {'utc_time': sat_time[net_idx], 'c1': c1[net_idx],
                    'c2': c2[net_idx], 'c3': c3[net_idx],
                    'net_idx': net_idx}
    # net_idx for comparison with other data from real satellite

    # perform coordinate conversions and sort variables by coordinate systems
    coord_dict = coordinate_systems(model, sat_time, c1, c2, c3, variable_list,
                                    coord_type, coord_grid)

    # perform flythroughs
    if verbose:
        print('Interpolating through model data...', end="")
    interp_time = perf_counter()
    for key in coord_dict.keys():
        # interpolate requested data for each day.
        # reader, file_name, variable_list, sat_time in hrs, c1, c2, c3,
        #   z_dependencies
        list_results = [Model_FlyAway(reader, times[file_date][0],
                                      coord_dict[key][0],
                                      ts_to_hrs(sat_time[times[file_date][6]],
                                                file_date.split('_')[0]),
                                      coord_dict[key][1][times[file_date][6]],
                                      coord_dict[key][2][times[file_date][6]],
                                      coord_dict[key][3][times[file_date][6]],
                                      coord_dict[key][4], verbose=verbose)
                        for file_date in times.keys() if
                        len(sat_time[times[file_date][6]]) > 0]

        # get new variable list from results dictionaries
        newvar_list = []
        [newvar_list.extend(list(results.keys())) for results in list_results]

        # collect interpolated data into the same dictionary
        for var in newvar_list:  # sort and combine arrays for the same var
            results_dict[var] = concatenate(tuple([results[var] for results
                                                   in list_results]))
    if verbose:
        print(f'done in {perf_counter()-interp_time:.5f} s.')
    return results_dict


def Prepare_Files(model, file_dir, call_type='normal'):
    '''Return a list of the required height input for each variable. Create
    wrapped files if needed.'''

    # Determine possible file patterns. Create wrapped files if needed.
    file_patterns = MW.FileSearch(model, file_dir, call_type=call_type)
    times = check_timescsv(file_patterns, model, call_type=call_type)
    # reader creates converted files if DNE
    return
