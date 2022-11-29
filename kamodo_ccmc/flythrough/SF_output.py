# -*- coding: utf-8 -*-
"""
Routines to output satellite flythrough time series into various file formats.
Oct 12, 2021: adding coordinate system to metadata of each file (in and out)
Mar 8, 2022: Adding function to convert outputs into Kamodo functions
"""
import numpy as np
from kamodo import Kamodo, kamodofy
from scipy.interpolate import interp1d
from datetime import datetime, timezone


@np.vectorize
def ts_to_datetime(time_val):
    '''Convert from utc timestamp to datetime object.'''
    return datetime.utcfromtimestamp(int(time_val)).replace(
        tzinfo=timezone.utc)


@np.vectorize
def datetime_to_utcts(date_time):
    '''Convert from datetime UTC object to utc timestamp.'''
    return datetime.timestamp(date_time)


def Functionalize_TimeSeries(utc_time, variable_name, variable_units,
                             variable_data, kamodo_object=None):
    '''Funtionalize the given time series data, return in a kamodo object.

    utc_time = a 1D array of utc timestamps.
    variable_name = a string representation of the name of the variable.
    variable_units = a string representation of the units of the variable.
    variable_data = a 1D array of the variable data of the same length as the
        utc_time array.
    Each command returns a kamodo object with the data functionalized. If you
        want additional functions added to
    a pre-existing kamodo object, then input the desired kamodo object through
        the kamodo_object keyword.
    '''

    # Perform the conversion on the array of UTC timestamps
    time_date = ts_to_datetime(utc_time)

    # Define an interpolator for the given time series array.
    interp = interp1d(utc_time, variable_data, bounds_error=False,
                      fill_value=np.NaN)

    # Functionalize the time series array
    @kamodofy(units=variable_units, data=variable_data)  # units
    # function depends on array of datetime objects
    def timeseries_func(time=time_date):
        # convert the given datetime objects to UTC timestamps
        utc_ts = datetime_to_utcts(time)
        return interp(utc_ts)     # return the scipy interpolator

    # Define a Kamodo object to include the desired functions.
    if kamodo_object is None:
        kamodo_object = Kamodo()
    kamodo_object[variable_name] = timeseries_func
    return kamodo_object


def Functionalize_SFResults(model, results, kamodo_object=None):
    '''Functionalizes all non-coordinate variables in the results dictionary.

    model = name of model (string) from outut of MW.Choose_Model('').
    results = dictionary returned from any flythrough function.
    kamodo_object (default=None) = a predefined Kamodo object to add
       functionalized data to. If none is provided, a new one will be created.
    '''
    import kamodo_ccmc.flythrough.model_wrapper as MW

    if kamodo_object is None:
        kamodo_object = Kamodo()
    variable_list = [item for item in list(results.keys()) if item not in
                     ['utc_time', 'c1', 'c2', 'c3', 'net_idx', 'metadata']]
    for varname in variable_list:
        if isinstance(results['utc_time'], (dict)):   # file version
            kamodo_object = Functionalize_TimeSeries(
                results['utc_time']['data'], varname,
                MW.Var_units(model, varname)[varname],
                results[varname]['data'],
                kamodo_object=kamodo_object)
        else:  # version of results directly from flythrough
            kamodo_object = \
                Functionalize_TimeSeries(results['utc_time'], varname,
                                         MW.Var_units(model, varname)[varname],
                                         results[varname],
                                         kamodo_object=kamodo_object)
    return kamodo_object


def SFcdf_reader(filename):
    '''Loads the data from a cdf file that was written by the SFdata_tocdf
    routine below into a nested dictionary.'''
    from netCDF4 import Dataset

    cdf_data = Dataset(filename, 'r')
    cdf_dict = {key: {'units': cdf_data.variables[key].units,
                      'data': np.array(cdf_data.variables[key])}
                for key in cdf_data.variables.keys()}
    cdf_dict['metadata'] = {'model_files': cdf_data.modelfile,
                            'model_used': cdf_data.model,
                            'coord_type': cdf_data.coord_type,
                            'coord_grid': cdf_data.coord_grid}  # add metadata
    cdf_data.close()
    return cdf_dict


def SFdata_tocdf(filename, model_filename, model_name, results_dict,
                 results_units, coord_type, coord_grid):
    '''Write satellite flythrough time series data to a netCDF4 file.'''

    from netCDF4 import Dataset

    # start new output object
    data_out = Dataset(filename, 'w', format='NETCDF4')
    data_out.modelfile = model_filename
    data_out.model = model_name
    data_out.coord_type = coord_type
    data_out.coord_grid = coord_grid

    # store time dimension
    for key in results_dict:
        if 'time' in key:
            time_key = key
            break
    new_dim = data_out.createDimension(time_key, results_dict[time_key].size)
    # create dimension

    # copy over variables to file
    for key in results_dict.keys():
        new_var = data_out.createVariable(key, np.float64, (time_key))
        new_var[:] = results_dict[key]  # store data in variable
        new_var.units = results_units[key]  # store units for variable

    # close file
    data_out.close()
    return filename


def SFcsv_reader(filename, delimiter=','):
    '''Loads the data from a csv file that was written by the SFdata_tocsv
    routine below into a nested dict.'''

    # open file
    from csv import reader
    read_obj = open(filename, 'r')
    csv_reader = reader(read_obj, delimiter=delimiter)

    # sort out header
    model_files = next(csv_reader)
    model_used = next(csv_reader)
    coord_info = next(csv_reader)
    variable_keys = next(csv_reader)
    variable_keys[0] = variable_keys[0][1:]  # cut out # from first name
    variable_units = next(csv_reader)
    variable_units[0] = variable_units[0][1:]  # cut out # from first unit
    # trim off empty variable keys with empty variable units
    while variable_keys[-1] == '':
        variable_keys = variable_keys[:-1]
        variable_units = variable_units[:-1]
    # trim [ ] off unit strings
    trimmed_units = [string[1:-1] for string in variable_units]

    # create dictionary to store data in
    data_dict = {variable_keys[i]: {'units': trimmed_units[i], 'data': []}
                 for i in range(len(variable_keys))}

    # store data into dictionary
    for row in csv_reader:
        for i in range(len(variable_keys)):  # skip empty block(s) at the end
            data_dict[variable_keys[i]]['data'].append(row[i])

    # convert to numpy float arrays, except for net_idx
    for key in data_dict.keys():
        if key == 'net_idx':
            data_dict[key]['data'] = np.array(data_dict[key]['data'],
                                              dtype=int)
        else:
            data_dict[key]['data'] = np.array(data_dict[key]['data'],
                                              dtype=float)

    # add metadata
    data_dict['metadata'] = {'model_files': model_files[1].strip(),
                             'model_used': model_used[1].strip(),
                             'coord_type': coord_info[1].strip(),
                             'coord_grid': coord_info[2].strip()}
    read_obj.close()
    return data_dict


def SFdata_tocsv(filename, model_filename, model_name, results_dict,
                 results_units, coord_type, coord_grid):
    '''Write satellite flythrough time series data to a csv file'''

    # get key name for time information
    for key in results_dict:
        if 'time' in key:
            time_key = key
            break

    data_out = open(filename, 'w')
    if not isinstance(model_filename, list):
        data_out.write(f'#Model files used:, {model_filename}')
    else:
        data_out.write('#Model files used:,'+''.join([f+',' for f in
                                                      model_filename]))
    data_out.write(f'\n#Model used:, {model_name}')
    data_out.write(f'\n#Coordinates:, {coord_type}, {coord_grid}')
    data_out.write('\n#'+''.join([key+',' for key in results_dict.keys()]))
    data_out.write('\n#'+''.join(['['+results_units[key]+'],' for key in
                                  results_dict.keys()]))
    for i in range(len(results_dict[time_key])):
        data_out.write('\n'+''.join([f'{values[i]},' for key, values in
                                     results_dict.items()]))
    data_out.close()
    return filename


def SFascii_reader(filename):
    '''Loads the data from a csv file that was written by the SFdata_toascii
    routine below into a nested dict'''

    return SFcsv_reader(filename, delimiter='\t')


def SFdata_toascii(filename, model_filename, model_name, results_dict,
                   results_units, coord_type, coord_grid):
    '''Write satellite flythrough time series data to a txt file'''

    # get key name for time information
    for key in results_dict:
        if 'time' in key:
            time_key = key
            break

    data_out = open(filename, 'w')
    if not isinstance(model_filename, list):
        data_out.write(f'#Model files used:\t {model_filename}')
    else:
        data_out.write('#Model files used:\t'+''.join([f+'\t' for f in
                                                       model_filename]))
    data_out.write(f'\n#Model used:\t {model_name}')
    data_out.write(f'\n#Coordinates:\t {coord_type}\t {coord_grid}')
    data_out.write('\n#'+''.join([key+'\t' for key in results_dict.keys()]))
    data_out.write('\n#'+''.join(['['+results_units[key]+']\t' for key in
                                  results_dict.keys()]))
    for i in range(len(results_dict[time_key])):
        data_out.write('\n'+''.join([f'{values[i]}\t' for key, values in
                                     results_dict.items()]))
    data_out.close()
    return filename


def SF_read(filename):
    '''Collect input function calls into one function.

    filename = string with complete filepath. The file extension must be one
        of 'nc'for a netCDF4 file, 'csv' for a comma separated file, or 'txt'
        for a tab separated file.

    Output: a nested dictionary containing the metadata, data, and units.'''

    file_type = filename.split('.')[-1]
    if file_type == 'nc':
        traj_data = SFcdf_reader(filename)
    elif file_type == 'csv':
        traj_data = SFcsv_reader(filename)
    elif file_type == 'txt':
        traj_data = SFascii_reader(filename)
    else:
        raise AttributeError('File type not recognized. Must be one of' +
                             ' cdf4, csv, or txt.')
    return traj_data


def SF_write(filename, model_filename, model_name, results_dict, results_units,
             coord_sys):
    '''Collect output function calls into one function.

    Inputs:
        filename = string with complete filepath. The file extension must be
            one of 'nc' for a netCDF4 file, 'csv' for a comma separated file,
            or 'txt' for a tab separated file.
        model_filename = A list of the model data filenames or prefixes used
            to generate the data. Filenames should include the full file path.
        model_name = A string indicating the model name.
        results_dict = A dictionary with variable names as keys (strings) and
            the time series data as the values (one array per key).
        results_units = A dictionary with variable names as keys (strings) and
            the units as the values (one value per key).
        coord_sys = one of 'GDZ', 'GEO', 'GSM', 'GSE', 'SM', 'GEI', 'MAG',
            'SPH', or 'RLL' combined with '-sph' or '-car'. E.g. 'SM-car' or
            'GDZ-sph'. Astropy coordinate systems supported. See ConvertCoord
            for details.
    '''
    coord_type, coord_grid = coord_sys.split('-')
    output_type = filename.split('.')[-1]
    if output_type == 'csv':
        output_filename = SFdata_tocsv(filename, model_filename, model_name,
                                       results_dict, results_units, coord_type,
                                       coord_grid)
    elif output_type == 'nc':
        output_filename = SFdata_tocdf(filename, model_filename, model_name,
                                       results_dict, results_units, coord_type,
                                       coord_grid)
    elif output_type == 'txt':
        output_filename = SFdata_toascii(filename, model_filename, model_name,
                                         results_dict, results_units,
                                         coord_type, coord_grid)
    elif output_type == 'CAMEL':
        # will put a call to the coming CAMEL output routine here (or diff place?)
        pass
    
    elif output_type == 'coupling':
        import kamodo_ccmc.flythrough.model_wrapper as MW
        output_routine = MW.Model_Coupling(model_name)
        output_filename = output_routine()
        
        
    return output_filename
