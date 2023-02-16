# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 14:57:55 2021
author: rringuet

Satellite Flythrough code to call for any model.
All desired model data for the desired model output run should be in a single
directory.

Example calls from the command line:
    Get list of models and function calls possible:
            python ./SatelliteFlythrough_new.py
    Get call syntax for one of the functions:
            python ./SatelliteFlythrough_new.py MyFlight
    Get information on what variables are available for a given model:
        python ./SatelliteFlythrough_new.py CTIPe
    For funtion calls, all of the parameters in the function call should be
        given after the name of the function in the example above, even if the
        default value is desired
See ConvertCoord for details on the coordinate systems.
"""
import numpy as np
from os.path import isfile
import kamodo_ccmc.flythrough.SF_output as O
import kamodo_ccmc.flythrough.SF_utilities as U


def SatelliteTrajectory(dataset, start_ts, stop_ts, coord_type='GEO',
                        verbose=False):
    '''Retrieve and return satellite trajectory from HAPI/CDAWeb
    Parameters:
    ----------
    dataset: name of the satellite data set to pull trajectory from
    start_ts: utc timestamp for start of desired time interval
    stop_ts: utc timestamp for end of desired time interval
    coord_type: Pick from GEO, GSM, GSE, or SM
    verbose: Set to true to be overwhelmed with information.

    Returns a dictionary with keys: sat_time, c1, c2, and c3.
        sat_time is an array in UTC seconds since 1970-01-01.
        (c1,c2,c3) = (lon, lat, alt) in (deg,deg,km) in the 'GDZ', 'sph'
        coordinate system in SpacePy. See
        kamodo_ccmc.flythrough.utils.ConvertCoord for more info on the
        coordinate systems.

    Coordinates are retrieved on a cartesian grid.
    See kamodo_ccmc.flythrough.utils.ConvertCoord for info on the coordinate
    systems.
    '''
    from kamodo_ccmc.readers.hapi import HAPI

    # convert from utc timestamps to isoformat
    start = U.ts_to_ISOstring(start_ts)
    stop = U.ts_to_ISOstring(stop_ts)

    # check input coord_type
    if coord_type not in ['GEO', 'GSM', 'GSE', 'SM']:
        raise AttributeError(f'Coordinate type {coord_type} not available. ' +
                             'Pick from GEO, GSM, GSE, or SM.')
    else:
        parameters = 'X_'+coord_type+',Y_'+coord_type+',Z_'+coord_type

    # retrieve satellite trajectory
    server = 'http://hapi-server.org/servers/SSCWeb/hapi'  # for coord data
    hapi = HAPI(server, dataset, parameters, start, stop, verbose=verbose)
    satellite_dict = {'sat_time': hapi.tsarray,  # utc timestamps
                      'c1': hapi.variables[parameters.split(',')[0]]['data'],
                      'c2': hapi.variables[parameters.split(',')[1]]['data'],
                      'c3': hapi.variables[parameters.split(',')[2]]['data']}

    print(f'Attribute/Key names of return dictionary: {satellite_dict.keys()}')

    return satellite_dict, coord_type+'-car'


def SampleTrajectory(start_time, stop_time, max_lat=65., min_lat=-65.,
                     lon_perorbit=363., max_height=450., min_height=400.,
                     p=0.01, n=2.):
    '''
    Given start and stop times in timestamp form, return a test satellite
    trajectory.
    Parameters:
    ----------
        start_time: utc timestamp in seconds for start
        stop_time: utc timestamp in seconds for stop
        max_lat: maximum latitude for sample trajectory, in degrees
            (default=65.)
        min_lat: minimum latitude for sample trajectory, in degrees
            (default=-65.)
        lon_perorbit: the degrees of longitude per about 90 minute orbit
            (set less than 360 for precession forward in longitude, set less
            than 360 for precession backwards) (default=363.)
        max_height: maximum starting height of orbit in km (default=450.)
        min_height: minimum starting height of orbit in km (default=400.)
        p: a rough precession variable, applied as an overall height decrease
            as a percentage of the min_height value: p =  (default=0.01).
        n: the time cadence of the sample trajectory generated
            (default = 2 seconds)
    Returns a dictionary with keys: sat_time, c1, c2, and c3.
        sat_time is an array in UTC seconds since 1970-01-01.
        (c1,c2,c3) = (lon, lat, alt) in (deg,deg,km) in the 'GDZ', 'sph'
        coordinate system in SpacePy. See
        kamodo_ccmc.flythrough.utils.ConvertCoord for more info on the
        coordinate systems.
    '''

    # determine basic parameters
    # determine number of samples per 90min orbit
    orbit_seconds = int(90. * 60. / float(n))
    n_orbits = (stop_time - start_time) / float(orbit_seconds*n)
    h_scale, h_offset = (max_height - min_height)/2., np.mean([max_height,
                                                               min_height])
    lat_scale, lat_offset = (max_lat-min_lat)/2., np.mean([max_lat, min_lat])
    time_left = (stop_time-start_time) / float(n) -\
        int(n_orbits) * orbit_seconds

    # create orbital tracks
    pi_arr = np.linspace(0., 2.*np.pi, orbit_seconds)
    lat = np.tile(np.cos(pi_arr), int(n_orbits))
    height = np.tile(np.sin(pi_arr), int(n_orbits))
    if time_left > 0:  # add partial last orbit
        lat = np.append(lat, np.cos(pi_arr[0:int(time_left)]))
        height = np.append(height, np.sin(pi_arr[0:int(time_left)]))
    lon = np.linspace(0., float(lon_perorbit)*n_orbits,
                      int((stop_time-start_time)/float(n)))
    while max(lon) > 360.:
        lon[np.where(lon > 360.)[0]] -= 360.
    while max(lon) < 0.:
        lon[np.where(lon < 0.)[0]] += 360.
    height = height * h_scale + h_offset -\
        np.linspace(0., p, int((stop_time-start_time)/float(n))) * min_height

    # store results in dictionary to return
    sample_dict = {'sat_time': np.linspace(start_time, stop_time,
                                           int((stop_time-start_time) /
                                               float(n))),
                   'c1': lon-180., 'c2': lat*lat_scale+lat_offset,
                   'c3': height}

    print(f'Attribute/Key names of return dictionary: {sample_dict.keys()}')
    print('(c1,c2,c3) = (lon, lat, alt) in (deg,deg,km) in the GDZ, sph ' +
          'coordinate system.\nsat_time contains the utc timestamps.')
    return sample_dict, 'GDZ-sph'


def TLETrajectory(tle_file, start_utcts, stop_utcts, time_cadence,
                  method='forward', verbose=False):
    '''Use sgp4 to calculate a satellite trajectory given TLEs.
    Parameters:
        tle_file: The file name, including complete file path, of a file
            containing two-line elements. It is assumed that the file has no
            header and no other content.
        start_utcts: The UTC timestamp corresponding to the desired start time.
            Should be an integer.
        stop_utcts: The UTC timestamp corresponding to the desired stop time.
            Should be an integer.
        time_cadence: The number of seconds desired between trajectory
            positions. Should be an integer.
        method: 'forward' or 'nearest'. This keyword changes the propagation
            method for timestamps between TLEs, not for timestamps before the
            first TLE or after the last TLE. The 'forward' (default) method
            uses the previous TLE to propagate forward for all timestamps
            between the selected TLE and the next one, while the 'nearest'
            method finds the TLE nearest to the timestamp and propagates either
            forward or backward for the timestamp.

    If the time cadence does not evenly divide into the range of timestamps
        given, then the ending time value will be extended so that the entire
        requested range will be covered.

    Returns a dictionary with keys: sat_time, c1, c2, and c3.
        sat_time is an array in UTC seconds since 1970-01-01.
        (c1,c2,c3) = (x, y, z) in (km,km,km) in the 'teme', 'car'
        coordinate system in AstroPy. See
        kamodo_ccmc.flythrough.utils.ConvertCoord for more info on the
        coordinate systems.
    '''
    from sgp4.api import days2mdhms, Satrec
    from pandas import to_datetime
    from datetime import datetime, timezone
    from astropy.constants import R_earth

    if verbose:
        print('Using the sgp4 propagator...')

    # Create an array of timestamps with the requested spacing.
    # Initialize return structure.
    n = int((stop_utcts-start_utcts)/time_cadence)
    if (start_utcts+time_cadence*n < stop_utcts):
        n += 1
        stop_utcts = start_utcts+time_cadence * n
    UTCtimestamps = np.linspace(start_utcts, stop_utcts, n+1)
    results = {'sat_time': UTCtimestamps,
               'c1': np.zeros(UTCtimestamps.shape, dtype=float),
               'c2': np.zeros(UTCtimestamps.shape, dtype=float),
               'c3': np.zeros(UTCtimestamps.shape, dtype=float)}
    # Keep independent or the arrays will be linked together in memory!

    # Convert timestamp array to decimal Julian dates using pandas.
    @np.vectorize
    def UTCtstoJulianDate(UTC_ts):
        '''Converts a UTC timestamp to a Julian date using pandas.
        Input: UTC timestamp integer value.
        Output: Julian date calculated by pandas.
        '''

        dt = datetime.utcfromtimestamp(UTC_ts).replace(tzinfo=timezone.utc)
        pd_dt = to_datetime(dt)
        return pd_dt.to_julian_date()
    julian_dates = UTCtstoJulianDate(UTCtimestamps)

    # Read in TLEs
    tle_data = {'TLE_line1': [], 'TLE_line2': [],
                'Element Set Epoch (UTC)': []}
    with open(tle_file) as file:
        for line in file:
            tle_data['TLE_line1'].append(line.strip())
            tle_data['Element Set Epoch (UTC)'].append(line[18:32].strip())
            line = next(file)
            tle_data['TLE_line2'].append(line.strip())

    # Calculate UTC timestamps for TLEs.
    tmp = []
    for item in tle_data['Element Set Epoch (UTC)']:
        yr = int(item[0:2])
        if yr >= 57:
            year = 1900 + yr
        else:
            year = 2000 + yr
        month, day, hr, minute, sec = days2mdhms(year, float(item[2:]))
        ts = datetime(year, month, day, hr, minute, round(sec),
                      tzinfo=timezone.utc).timestamp()
        tmp.append(ts)
    tle_data['UTC_timestamps'] = np.array(tmp, dtype=int)
    if verbose:
        print('TLE data:')
        for key in tle_data.keys():
            print(key, tle_data[key])

    # If nearest propagation method chosen, prepare timestamp assignments
    if method == 'nearest':
        ts_nearest = np.zeros((len(tle_data['UTC_timestamps']),
                               len(UTCtimestamps)), dtype=int)
        for i in range(len(tle_data['UTC_timestamps'])):
            ts_nearest[i, :] = abs(UTCtimestamps -
                                   tle_data['UTC_timestamps'][i])
        ts_assignments = np.argmin(ts_nearest, axis=0)

    # Assign timestamps to TLE indices
    for i in range(len(tle_data['UTC_timestamps'])+1):
        # Determine relevant timestamps in user array
        if method == 'forward':  # default method typically used in forecasting
            if i == 0:  # all timestamps before first TLE, propagate backwards
                idx = np.where(UTCtimestamps <
                               tle_data['UTC_timestamps'][i])[0]
                if len(idx) > 0:
                    tle_idx = 0
                    test_diff = abs(UTCtimestamps[idx] -
                                    tle_data['UTC_timestamps'][i])
                    diff_idx = np.where(test_diff > 86400)[0]
                    if len(diff_idx) > 0:
                        print(f'Caution: {len(diff_idx)} times are more than' +
                              ' 24 hrs before the first TLE.')
                else:
                    if verbose:
                        print('Skipped TLE index:', i-1)
                    continue
            elif i < len(tle_data['UTC_timestamps']):  # propagate forwards
                idx = np.where((UTCtimestamps >=
                                tle_data['UTC_timestamps'][i-1]) &
                               (UTCtimestamps <
                                tle_data['UTC_timestamps'][i]))[0]
                tle_idx = i-1
                if len(idx) < 1:
                    if verbose:
                        print('Skipped TLE index:', i-1)
                    continue
            elif i == len(tle_data['UTC_timestamps']):  # propagate forwards
                idx = np.where(UTCtimestamps >=
                               tle_data['UTC_timestamps'][i-1])[0]
                tle_idx = i-1  # same as last loop
                if len(idx) > 0:
                    test_diff = abs(UTCtimestamps[idx] -
                                    tle_data['UTC_timestamps'][i-1])
                    diff_idx = np.where(test_diff > 86400)[0]
                    if len(diff_idx) > 0:
                        print(f'Caution: {len(diff_idx)} times are more than' +
                              ' 24 hrs after the last TLE.')
                else:
                    if verbose:
                        print('Skipped TLE index:', i-1)
                    continue
        # more accurate method when future TLEs are known
        elif method == 'nearest':
            idx = np.where(ts_assignments == i)[0]
            if len(idx) > 0:
                tle_idx = i
                if (i == 0) or (i == len(tle_data['UTC_timestamps'])-1):
                    test_diff = abs(UTCtimestamps[idx] -
                                    tle_data['UTC_timestamps'][i])
                    diff_idx = np.where(test_diff > 86400)[0]
                    if len(diff_idx) > 0:
                        if i == 0:
                            print(f'Caution: {len(diff_idx)} times are more' +
                                  ' than 24 hrs before the first TLE.')
                        if i == len(tle_data['UTC_timestamps'])-1:
                            print(f'Caution: {len(diff_idx)} times are more' +
                                  ' than 24 hrs after the last TLE.')
            elif len(idx) < 1:
                continue

        # Initialize TLE information with chosen tle
        if verbose:
            print('TLE assignment:', tle_idx, min(idx), max(idx))
        satellite = Satrec.twoline2rv(tle_data['TLE_line1'][tle_idx],
                                      tle_data['TLE_line2'][tle_idx])

        # Calculate position vector for each relevant timestamp
        @np.vectorize
        def calc_position(jd, satellite):
            '''Use sgp4 to calculate the position vector for a given julian
                date and initiliazed satellite instance.
            Returns r, the position in km in the idiosyncratic True Equator
            Mean Equinox (TEME) coordinate frame.'''
            e, r, v = satellite.sgp4(jd, 0.0)
            return r

        # Convert to earth radii and store
        R_earth_km = R_earth.value/1000.
        results['c1'][idx], results['c2'][idx], results['c3'][idx] = \
            calc_position(julian_dates[idx], satellite) / R_earth_km

    return results, 'teme-car'


# want to enable call of this from C++ for flexibility, so return only one item
# keep so users can call this if they have their own satellite trajectory data
def ModelFlythrough(model, file_dir, variable_list, sat_time, c1, c2, c3,
                    coord_sys, output_name='', plot_coord='GEO',
                    verbose=False, _print_units=True):
    '''Call satellite flythrough wrapper specific to the model chosen.
    Parameters:
    ------------
    model: 'CTIPe','IRI', ...
    file_dir: complete path to where model data files are stored
    variable_list: List of standardized variable names. See model variable
        output for details. Dimensions must be at least time + 2D spatial.
    sat_time: a numpy array of the utc timestamps
    c1, c2, c3: numpy arrays of the positions correlating to the utc timestamps
        (c1, c2, c3) should be (x,y,z) in R_E for cartesian coordinates, and
        (lon, lat, radius (R_E) or altitude (km)) for spherical coordinates.
    coord_sys: one of 'GDZ', 'GEO', 'GSM', 'GSE', 'SM', 'GEI', 'MAG', 'SPH',
        'RLL', etc, combined with '-sph' or '-car'. E.g. 'SM-car' or 'GDZ-sph'.
    output_name: complete path with filename (with the extension) for the file
        to write the results to. Plotting filenames are determined by
        output_name - extension + variable names + plot type (1D or 3D).
        Extensions must be one of 'nc' for netCDF4 files, 'csv' for comma
        separated data, or 'txt' for tab separated data files.
    plot_coord: one of 'GDZ', 'GEO', 'GSM', 'GSE', 'SM', 'GEI', 'MAG'.
        Indicates the coordinate system the plot will be generated in. Only
        plots in cartesian coordinates systems are supported, so 'SPH' and
        'RLL' are not accepted. Default is 'GEO'.
    verbose: Set to true to be overwhelmed with information.

    Returns a dictionary with keys: 'utc_time', 'c1', 'c2', 'c3', and 'net_idx'
        - utc_time is an array in UTC seconds since 1970-01-01 of the given
            timestamps with any occuring outside of the model data removed.
        - 'c1', 'c2', and 'c3' are arrays of the given coordinate values for
            each surviving timestamp
        - 'net_idx' is the original index value of the surviving timestamps.
            This is kept for easier comparison with the original dataset.
        - additional keys are included for each variable and label an array of
            the values of the indicated variable for each time+spatial
            coordinate given.
        - The units of each array in the returned dictionary are printed to
            the screen.

    See kamodo_ccmc.flythrough.utils.ConvertCoord for info on the coordinate
    systems.'''

    # if input types are lists, correct to be numpy arrays
    # (important for calling from C++)
    if isinstance(sat_time, list):
        sat_time = np.array(sat_time)
    if isinstance(c1, list):
        c1 = np.array(c1)
    if isinstance(c2, list):
        c2 = np.array(c2)
    if isinstance(c3, list):
        c3 = np.array(c3)

    # give error if unknown output type given BEFORE running flythrough
    output_type = output_name.split('.')[-1]
    output_typelist = ['nc', 'csv', 'txt', '']
    if output_type not in output_typelist:
        raise AttributeError('Output extension not recognized. Must be one' +
                             ' of ', *output_typelist)

    # Check if output file already exists. Break if it does.
    if isfile(output_name):
        raise AttributeError(f'{output_name} already exists. ' +
                             'Please either remove the file or change the ' +
                             'output_name and then rerun the command.')

    # correct variables with 'ilev' in the name 
    new_list = [''.join([i+'_' for i in var.split('_')[:-1]])[:-1]
                if 'ilev' in var else var for var in variable_list]
    coord_type, coord_grid = coord_sys.split('-')

    # get interpolated results
    # coord_type should be one of SpacePy's or AstroPy's coordinates
    # coord_grid is either 'sph' or 'car'
    results = U.Model_SatelliteFlythrough(model, file_dir, new_list,
                                          sat_time, c1, c2, c3,
                                          coord_type, coord_grid,
                                          verbose=verbose)

    # remove requested variables not found in the data
    var_list = [key for key in results.keys() if key not in
                ['utc_time', 'c1', 'c2', 'c3', 'net_idx']]
    if len(var_list) == 0:
        return results

    # retrieve coordinate and results units
    coord_units = U.MW.coord_units(coord_type, coord_grid)
    results_units = U.MW.Var_units(model, var_list)
    for key in coord_units:
        results_units[key] = coord_units[key]
    if _print_units:
        print(results_units)

    if verbose:
        print(f'Units from the {model} model by variable name:' +
              f'\n{results_units}')
        print(f'Dictionary key names in results:\n{results.keys()}')
        print('The units of the trajectory variables are unchanged from the ' +
              'inputs.')

    if output_name != '':
        # retrieve list of files used in the execution
        filenames = U.MW.File_List(model, file_dir)
        # perform output type desired
        output_filename = O.SF_write(output_name, filenames, model, results,
                                     results_units, coord_sys)
        # no access to model filenames
        print(f"Output saved in {output_filename}.")
        print('Generating interactive plots...')

        # check plot_coord variable
        if plot_coord in ['SPH', 'RLL']:
            raise AttributeError('Plots can only be requested in coordinate ' +
                                 'grids where cartesian coordinates are ' +
                                 "supported. The 'SPH' and 'RLL' coordinate " +
                                 'systems do not support cartesian grids so ' +
                                 'are not allowed.')

        # generate and save plots without displaying
        from kamodo_ccmc.flythrough.plots import SatPlot4D
        # presentation options: all, day, hour, minute, N, orbitE, orbitM
        plot_output = output_name.split('.'+output_type)[0]  # remove extension
        for var in var_list:
            SatPlot4D(var, results['utc_time'], results['c1'], results['c2'],
                      results['c3'], results[var], results_units[var],
                      coord_type, coord_grid, plot_coord, 'all', model,
                      body='black', divfile=plot_output+'_'+var+'_3D.html',
                      displayplot=False)
            SatPlot4D(var, results['utc_time'], results['c1'], results['c2'],
                      results['c3'], results[var], results_units[var],
                      coord_type, coord_grid, plot_coord, 'all', model,
                      type='1D', divfile=plot_output+'_'+var+'_1D.html',
                      displayplot=False)
    return results


def FakeFlight(start_time, stop_time, model, file_dir, variable_list,
               max_lat=65., min_lat=-65., lon_perorbit=363., max_height=450.,
               min_height=400., p=0.01, n=2., verbose=False,
               output_name='', plot_coord='GEO'):
    '''Generates a sample trajectory and then flies that trajectory through the
    model data chosen.

    Parameters:
        start_time: utc timestamp in seconds for start
        stop_time: utc timestamp in seconds for stop
        model: CTIPe, IRI, ....
        file_dir: complete path to where model data is stored
        variable_list: list of standardized variable names desired.
            Integers allowed. Dimensions must be at least time + 2D spatial.
        max_lat: maximum latitude for sample trajectory, in degrees
            (default=65.)
        min_lat: minimum latitude for sample trajectory, in degrees
            (default=-65.)
        lon_perorbit: the degrees of longitude per about 90 minute orbit
            (set less than 360 for precession forward in longitude, set less
            than 360 for precession backwards) (default=363.)
        max_height: maximum starting height of orbit in km (default=450.)
        min_height: minimum starting height of orbit in km (default=400.)
        p: a rough precession variable, applied as an overall height decrease
            as a percentage of the min_height value: p =  (default=0.01).
        n: the time cadence of the sample trajectory generated
            (default = 2 seconds)
        output_name: complete path with filename (with the extension) for the
            file to write the results to. Plotting filenames are determined by
            output_name - extension + variable names + plot type (1D or 3D).
            Extensions must be one of 'nc' for netCDF4 files, 'csv' for comma
            separated data, or 'txt' for tab separated data files.
        plot_coord: one of 'GDZ', 'GEO', 'GSM', 'GSE', 'SM', 'GEI', 'MAG', etc.
            Indicates the coordinate system the plot will be generated in.
            Only plots in cartesian coordinates systems are supported, so 'SPH'
            and 'RLL' are not accepted. Default is 'GEO'. Astropy coordinates
            also allowed.

    Returns a dictionary with keys: 'utc_time', 'c1', 'c2', 'c3', and 'net_idx'
    - utc_time is an array in UTC seconds since 1970-01-01 of the generated
        timestamps with any occuring outside of the model data removed.
    - 'c1', 'c2', and 'c3' are arrays of the given coordinate values for each
        surviving timestamp
    - 'net_idx' is the original index value of the surviving timestamps.
        This is kept for easier comparison with the original dataset.
    - additional keys are included for each variable and label an array of the
        values of the indicated variable for each time+spatial coordinate given
    - The units of each array in the returned dictionary are printed to
        the screen.

    See kamodo_ccmc.flythrough.utils.ConvertCoord for info on the coordinate
    systems.'''

    # generate a sample satellite trajectory
    sat_dict, coord_sys = SampleTrajectory(start_time, stop_time,
                                           max_lat=max_lat, min_lat=min_lat,
                                           lon_perorbit=lon_perorbit,
                                           max_height=max_height,
                                           min_height=min_height,
                                           p=p, n=n)

    # call satellite flythrough code
    results = ModelFlythrough(model, file_dir, variable_list,
                              sat_dict['sat_time'], sat_dict['c1'],
                              sat_dict['c2'], sat_dict['c3'],
                              coord_sys, verbose=verbose,
                              output_name=output_name,
                              plot_coord=plot_coord)
    return results


def RealFlight(dataset, start, stop, model, file_dir, variable_list,
               coord_type='GEO', output_name='', plot_coord='GEO',
               verbose=False):
    '''
    Retrieves the trajectory for the satellite requested and then flies that
    trajectory through the model data requested.

    dataset: name of the satellite data set to pull trajectory from
    start: utc timestamp for start of desired time interval
    stop: utc timestamp for end of desired time interval
    model: 'CTIPe','IRI', ...
    file_dir: complete path to where model data files are stored
    variable_list: List of standardized variable names. See model variable
        output for details. Dimensions must be at least time + 2D spatial.
    coord_type: Pick from GEO, GSM, GSE, or SM for the satellite trajectory.
    output_name: complete path with filename (with the extension) for the file
        to write the results to. Plotting filenames are determined by
        output_name - extension + variable names + plot type (1D or 3D).
        Extensions must be one of 'nc' for netCDF4 files, 'csv' for comma
        separated data, or 'txt' for tab separated data files.
    plot_coord: one of 'GDZ', 'GEO', 'GSM', 'GSE', 'SM', 'GEI', 'MAG', etc.
        Indicates the coordinate system the plot will be generated in.
        Only plots in cartesian coordinates systems are supported, so 'SPH'
        and 'RLL' are not accepted. Default is 'GEO'. Astropy coordinates
        also allowed.
    verbose: Set to true to be overwhelmed with information.

    Returns a dictionary with keys: 'utc_time', 'c1', 'c2', 'c3', and 'net_idx'
        - utc_time is an array in UTC seconds since 1970-01-01 of the satellite
            timestamps with any occuring outside of the model data removed.
        - 'c1', 'c2', and 'c3' are arrays of the given coordinate values for
            each surviving timestamp.
        - 'net_idx' is the original index value of the surviving timestamps.
            This is kept for easier comparison with the original dataset.
        - additional keys are included for each variable and label an array of
            the values of the indicated variable for each time+spatial
            coordinate given.
        - The units of each array in the returned dictionary are printed to
            the screen.

    See kamodo_ccmc.flythrough.utils.ConvertCoord for info on the coordinate
    systems.'''

    # retrieve satellite trajectory from HAPI/CDAWeb
    sat_dict, coord_sys = SatelliteTrajectory(dataset, start, stop,
                                              coord_type=coord_type,
                                              verbose=verbose)
    # call satellite flythrough code
    results = ModelFlythrough(model, file_dir, variable_list,
                              sat_dict['sat_time'], sat_dict['c1'],
                              sat_dict['c2'], sat_dict['c3'],
                              coord_sys, output_name=output_name,
                              plot_coord=plot_coord,
                              verbose=verbose)
    return results


def TLEFlight(tle_file, start, stop, time_cadence, model, file_dir,
              variable_list, output_name='', plot_coord='GEO',
              method='forward', verbose=False):
    '''
    Use sgp4 to calculate a satellite trajectory given TLEs, then fly the
    trajectory through the chosen model data. If the time cadence does not
    evenly divide into the range of timestamps given, then the ending time
    value will be extended so that the entire requested range will be covered.
    Parameters:
        tle_file: The file name, including complete file path, of a file
            containing two-line elements. It is assumed that the file has no
            header and no other content.
        start_utcts: The UTC timestamp corresponding to the desired start time.
            Should be an integer.
        stop_utcts: The UTC timestamp corresponding to the desired stop time.
            Should be an integer.
        time_cadence: The number of seconds desired between trajectory
            positions. Should be an integer.
        model: 'CTIPe','IRI', ...
        file_dir: complete path to where model data files are stored.
        variable_list: List of standardized variable names. See model variable
            output for details. Dimensions must be at least time + 2D spatial.
        coord_type: Pick from GEO, GSM, GSE, or SM for the satellite
            trajectory.
        output_name: complete path with filename (with the extension) for the
            file to write the results to. Plotting filenames are determined by
            output_name - extension + variable names + plot type (1D or 3D).
            Extensions must be one of 'nc' for netCDF4 files, 'csv' for comma
            separated data, or 'txt' for tab separated data files.
        plot_coord: Any of the coordinate systems in SpacePy or AstroPy.
            See the ConvertCoord function for more information on the choices
            of string values.
        method: 'forward' or 'nearest'. This keyword changes the propagation
            method for timestamps between TLEs, not for timestamps before the
            first TLE or after the last TLE. The 'forward' (default) method
            uses the previous TLE to propagate forward for all timestamps
            between the selected TLE and the next one, while the 'nearest'
            method finds the TLE nearest to the timestamp and propagates either
            forward or backward for the timestamp.
        verbose: Set to true to be overwhelmed with information.

    Returns a dictionary with keys: 'utc_time', 'c1', 'c2', 'c3', and 'net_idx'
        - utc_time is an array in UTC seconds since 1970-01-01 of the satellite
            timestamps with any occuring outside of the model data removed.
        - 'c1', 'c2', and 'c3' are arrays of the given coordinate values for
            each surviving timestamp
        - 'net_idx' is the original index value of the surviving timestamps.
            This is kept for easier comparison with the original dataset.
        - additional keys are included for each variable and label an array of
            the values of the indicated variable for each time+spatial
            coordinate given.
        - The units of each array in the returned dictionary are printed to
            the screen.

    See kamodo_ccmc.flythrough.utils.ConvertCoord for info on the coordinate
    systems.'''

    # retrieve satellite trajectory from HAPI/CDAWeb
    sat_dict, coord_sys = TLETrajectory(tle_file, start, stop, time_cadence,
                                        method=method, verbose=verbose)
    # call satellite flythrough code
    results = ModelFlythrough(model, file_dir, variable_list,
                              sat_dict['sat_time'], sat_dict['c1'],
                              sat_dict['c2'], sat_dict['c3'],
                              coord_sys, output_name=output_name,
                              plot_coord=plot_coord, verbose=verbose)
    return results


def MyFlight(traj_file, model, file_dir, variable_list, output_name='',
             plot_coord='GEO', verbose=False):
    '''Read in a trajectory from a file, then fly through the model data
    selected.

    traj_file: complete path and filename for file containing trajectory data.
    file_type: one of 'cdf4' for netCDF4 files, 'csv' for comma-separated
        files, or 'txt' for a tab-separated text file. Indicates the format of
        the input trajectory file.
    coord_sys: one of 'GDZ', 'GEO', 'GSM', 'GSE', 'SM', 'GEI', 'MAG', 'SPH',
        'RLL' combined with '-sph' or '-car'. E.g. 'SM-car' or 'GDZ-sph'.
        Astropy coordinate systems also supported.
    model: 'CTIPe', 'IRI', ...
    file_dir: complete path to model data files
    variable_list: List of standardized variable names. See model variable
        output for details. Dimensions must be at least time + 2D spatial.
    output_name: complete path with filename (with the extension) for the file
        to write the results to. Plotting filenames are determined by
        output_name - extension + variable names + plot type (1D or 3D).
        Extensions must be one of 'nc' for netCDF4 files, 'csv' for comma
        separated data, or 'txt' for tab separated data files.
    plot_coord: one of 'GDZ', 'GEO', 'GSM', 'GSE', 'SM', 'GEI', 'MAG', etc.
        Indicates the coordinate system the plot will be generated in.
        Only plots in cartesian coordinates systems are supported, so 'SPH'
        and 'RLL' are not accepted. Default is 'GEO'. Astropy coordinates
        also allowed.
    verbose: Set to true to be overwhelmed with information.

    Returns a dictionary with keys: 'utc_time', 'c1', 'c2', 'c3', and 'net_idx'
        - utc_time is an array in UTC seconds since 1970-01-01 of the given
            timestamps with any occuring outside of the model data removed.
        - 'c1', 'c2', and 'c3' are arrays of the given coordinate values for
            each surviving timestamp.
        - 'net_idx' is the original index value of the surviving timestamps.
            This is kept for easier comparison with the original dataset.
        - additional keys are included for each variable and label an array of
            the values of the indicated variable for each time+spatial
            coordinate given.
        - The units of each array in the returned dictionary are printed to
            the screen.

    See kamodo_ccmc.flythrough.utils.ConvertCoord for info on the coordinate
    systems.'''

    # read in trajectory from file into dictionary, including metadata
    traj_data = O.SF_read(traj_file)

    # figure out key for time data
    for key in traj_data:
        if 'time' in key:
            time_key = key
            break

    # call satellite flythrough code
    coord_sys = traj_data['metadata']['coord_type'] + '-' + \
        traj_data['metadata']['coord_grid']
    results = ModelFlythrough(model, file_dir, variable_list,
                              traj_data[time_key]['data'],
                              traj_data['c1']['data'], traj_data['c2']['data'],
                              traj_data['c3']['data'], coord_sys,
                              output_name=output_name, plot_coord=plot_coord,
                              verbose=verbose)

    return results


# allow calls from the command line
# these calls require all variables be given
if __name__ == '__main__':
    from sys import argv

    # print info if called without arguments
    if len(argv) == 2:
        if argv[1] == 'FakeFlight':
            help(FakeFlight)
        elif argv[1] == 'RealFlight':
            help(RealFlight)
        elif argv[1] == 'TLEFlight':
            help(TLEFlight)
        elif argv[1] == 'MyFlight':
            help(MyFlight)
        else:
            U.MW.Model_Variables(argv[1])
    elif len(argv) == 4:
        search_str = argv[1]
        model = argv[2]
        file_dir = argv[3]
        U.MW.Variable_Search(search_string=search_str, model=model,
                             file_dir=file_dir)
    elif len(argv) > 7:
        if argv[1] == 'FakeFlight':  # gather variables and call FakeFlight
            start_time = int(argv[2])
            stop_time = int(argv[3])
            if len(argv[4]) == 1:
                model = int(argv[4])
            else:
                model = argv[4]
            file_dir = argv[5]
            temp_str = argv[6][1:-1].replace("'",
                                             "").replace(' ',
                                                         '').replace('"', '')
            variable_list = temp_str.split(',')   # ['rho','N_n']
            max_lat = float(argv[7])
            min_lat = float(argv[8])
            lon_perorbit = float(argv[9])
            max_height = float(argv[10])
            min_height = float(argv[11])
            p = float(argv[12])
            n = float(argv[13])
            output_name = argv[14]
            plot_coord = argv[15]

            # check input
            print(f'\nstart_time: {start_time}, \nstop_time: {stop_time}, ' +
                  f'\nmodel: {model}, \nfile_dir: {file_dir}, ' +
                  f'\nvariable_list: {variable_list}, \nmax_lat: {max_lat}, ' +
                  f'\nmin_lat: {min_lat}, \nlon_perorbit: {lon_perorbit}, ' +
                  f'\nmax_height: {max_height}, \nmin_height {min_height}, ' +
                  f'\np: {p}, n: {n}, \noutput_name: {output_name}, ' +
                  f'\nplot_coord: {plot_coord}')

            results = FakeFlight(start_time, stop_time, model, file_dir,
                                 variable_list, max_lat=max_lat,
                                 min_lat=min_lat, lon_perorbit=lon_perorbit,
                                 max_height=max_height, min_height=min_height,
                                 p=p, n=n, verbose=False,
                                 output_name=output_name,
                                 plot_coord=plot_coord)

        elif argv[1] == 'RealFlight':  # gather variables and call RealFlight
            dataset = argv[2]
            start = int(argv[3])
            stop = int(argv[4])
            if len(argv[5]) == 1:
                model = int(argv[5])
            else:
                model = argv[5]
            file_dir = argv[6]
            temp_str = argv[7][1:-1].replace("'",
                                             "").replace(' ',
                                                         '').replace('"', '')
            variable_list = temp_str.split(',')   # ['rho','N_n']
            coord_type = argv[8]
            output_name = argv[9]
            plot_coord = argv[10]

            # check input
            print(f'\ndataset: {dataset}, \nstart: {start}, \nstop: {stop},' +
                  f'\nmodel: {model}, \nfile_dir: {file_dir},' +
                  f'\nvariable_list: {variable_list}, ' +
                  f'\ncoord_type: {coord_type},' +
                  f'\noutput_name: {output_name}, ' +
                  f'\nplot_coord: {plot_coord}')

            results = RealFlight(dataset, start, stop, model, file_dir,
                                 variable_list, coord_type=coord_type,
                                 output_name=output_name,
                                 plot_coord=plot_coord)

        elif argv[1] == 'TLEFlight':  # gather variables and call TLEFlight
            tle_file = argv[2]
            start = int(argv[3])
            stop = int(argv[4])
            time_cadence = int(argv[5])
            if len(argv[6]) == 1:
                model = int(argv[6])
            else:
                model = argv[6]
            file_dir = argv[7]
            temp_str = argv[8][1:-1].replace("'",
                                             "").replace(' ',
                                                         '').replace('"', '')
            variable_list = temp_str.split(',')   # ['rho','N_n']
            output_name = argv[9]
            plot_coord = argv[10]
            method = argv[11]

            # check input
            print(f'\ntle_file: {tle_file}, \nstart: {start}, ' +
                  f'\nstop: {stop},' +
                  f'\ntime_cadence: {time_cadence},' +
                  f'\nmodel: {model}, \nfile_dir: {file_dir},' +
                  f'\nvariable_list: {variable_list}, ' +
                  f'\noutput_name: {output_name}, ' +
                  f'\nplot_coord: {plot_coord}, ' +
                  f'\nmethod: {method}\n')

            results = TLEFlight(tle_file, start, stop, time_cadence,
                                model, file_dir, variable_list,
                                output_name=output_name,
                                plot_coord=plot_coord, method=method)

        elif argv[1] == 'MyFlight':  # gather variables and call MyFlight
            traj_file = argv[2]
            if len(argv[3]) == 1:
                model = int(argv[3])
            else:
                model = argv[3]
            file_dir = argv[4]
            temp_str = argv[5][1:-1].replace("'",
                                             "").replace(' ',
                                                         '').replace('"', '')
            variable_list = temp_str.split(',')   # ['rho','N_n']
            output_name = argv[6]
            plot_coord = argv[7]

            # check inputs
            print(f'\ntraj_file: {traj_file},' +
                  f'\nmodel: {model}, \nfile_dir: {file_dir},' +
                  f'\nvariable_list: {variable_list}, ' +
                  f'\noutput_name: {output_name}, ' +
                  f'\nplot_coord: {plot_coord}\n')

            results = MyFlight(traj_file, model, file_dir, variable_list,
                               output_name=output_name, plot_coord=plot_coord)
        else:
            print('Call signature not recognized.')
    else:
        print('\nPossible call types (first argument): FakeFlight, ' +
              'RealFlight, TLEFlight, MyFlight')
        print('Use the call type as the first input to get call syntax.\n')
        U.MW.Choose_Model('')  # asking for a list of possible models
