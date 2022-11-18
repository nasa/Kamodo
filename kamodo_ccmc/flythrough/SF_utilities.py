# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 17:09:14 2021
@author: rringuet
Most needed functions to support the SatelliteFlythrough software.
"""
from numpy import vectorize, array, diff, where, unique
from time import perf_counter
from datetime import datetime, timedelta
from os.path import isfile
import kamodo_ccmc.flythrough.model_wrapper as MW
from kamodo_ccmc.readers.reader_utilities import read_timelist
from kamodo_ccmc.flythrough.utils import ConvertCoord


@vectorize
def ts_to_hrs(utcts, filedate):
    '''Convert array of timestamps to hours since midnight of filedate.'''
    return (utcts - filedate.timestamp())/3600.


def ts_to_ISOstring(utc_ts):
    '''Convert timestamp to string of format 2017-05-28T00:00:00'''
    return datetime.utcfromtimestamp(utc_ts).isoformat()


def sat_tracks(sat_time, c1, c2, c3, z_dependencies, verbose=False):
    '''Calculate satellite tracks for interpolation'''
    # Create satellite tracks with appropriate inputs
    sat_track = {}  # initialize list of satellite tracks
    if '3D' in z_dependencies.keys():
        sat_track['3D'] = array([sat_time, c1, c2]).T
    if '4D' in z_dependencies.keys():
        sat_track['4D'] = array([sat_time, c1, c2, c3]).T
    return sat_track


def Model_FlyAway(reader, file_dir, variable_list, sat_time, c1, c2, c3,
                  z_dependencies, verbose=False):
    '''Perform flythrough for one day of data and one coordinate system.'''

    # create kamodo object, initialize some variables
    var_list = variable_list.copy()  # save copy before it gets altered
    kamodo_object = reader(file_dir, variables_requested=variable_list,
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
    1: x/lon 1D array, converted to the coordinate system key
    2: y/lat 1D array, converted to the coordinate system key
    3: z/rad/height 1D array, converted to the coordinate system key
    4: z_dependencies dictionary {3D: ([t, c1, c2]).T, 4D: ([t, c1, c2, c3]).T}
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

    # initialize model reader and times dictionary
    reader = MW.Model_Reader(model)  # Kamodo gets imported here
    start_utcts, end_utcts, filedate = File_UTCTimes(model, file_dir)
    
    # cut off trajectory times not found in data in file_dir
    idx = where((sat_time >= start_utcts) & (sat_time <= end_utcts))[0]
    results_dict = {'utc_time': sat_time[idx], 'c1': c1[idx],
                    'c2': c2[idx], 'c3': c3[idx], 'net_idx': idx}
    # net_idx for comparison with other data from real satellite
    if len(idx) < 1:
        print('None of the given times are in the model output files.')
        return results_dict
    if len(idx) < len(sat_time):
        print(f'{len(sat_time)-len(idx)} times are not in model output files' +
              ' and are excluded from the flythrough.')
    
    # perform coordinate conversions and sort variables by coordinate systems
    coord_dict = coordinate_systems(model, sat_time[idx], c1[idx], c2[idx],
                                    c3[idx], variable_list, coord_type,
                                    coord_grid)

    # perform flythroughs
    if verbose:
        print('Interpolating through model data...', end="")
    interp_time = perf_counter()
    for key in coord_dict.keys():  # keys = coordinate system names
        # interpolate requested data for each coordinate system
        # reader, file_dir, variable_list, sat_time in hrs, c1, c2, c3,
        #   z_dependencies
        
        results = Model_FlyAway(reader, file_dir, coord_dict[key][0],
                                ts_to_hrs(sat_time[idx], filedate),
                                coord_dict[key][1], coord_dict[key][2],
                                coord_dict[key][3], coord_dict[key][4],
                                verbose=verbose)

        # collect interpolated data into the same dictionary
        for var in results.keys():
            results_dict[var] = results[var]
    if verbose:
        print(f'done in {perf_counter()-interp_time:.5f} s.')
    return results_dict


def File_UTCTimes(model, file_dir):
    '''Returns the min start utc timestamp, the max end utc timestamp, and the
    datetime object containing midnight of the start date of the data. Creates
    preprocessed files if needed.'''

    # get times dictionary and datetime filedate object from files
    list_file = file_dir + model +'_list.txt'
    time_file = file_dir + model +'_times.txt'
    if isfile(list_file):
        times, tmp, filedate, tmp = read_timelist(time_file, list_file)
    else:
        reader = MW.Model_Reader(model)
        ko = reader(file_dir, filetime=True)  # creates any preprocessed files
        times, filedate = ko.times, ko.filedate
        del ko
        
    # calculate minimum start time
    start, end = [], []
    for p in times.keys():
        start.append(times[p]['start'][0])
        end.append(times[p]['end'][-1])
    start_utcts = (filedate + timedelta(hours=float(min(start)))).timestamp()
    end_utcts = (filedate + timedelta(hours=float(max(end)))).timestamp()
    return start_utcts, end_utcts, filedate
        
