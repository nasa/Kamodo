# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 18:17:56 2021
@author: rringuet
Convert data in CTIPe output files to wrapped versions
Oct 3: Add timestep from earlier file
Oct 5: calculate median height in km for each pressure level and store
"""
from numpy import transpose, zeros, array, append, where, insert, median
from time import perf_counter
from netCDF4 import Dataset
from os.path import isfile, basename
from datetime import datetime, timedelta

file_varnames = {'density': ['rho', 0, ['time', 'lon_d', 'lat_d', 'lev'],
                             'kg/m**3'],
                 'temperature': ['T', 1, ['time', 'lon_d', 'lat_d', 'lev'],
                                 'K'],
                 'electron_temperature': ['T_e', 2, ['time', 'lon_h', 'lat_h',
                                                     'height'], 'K'],
                 'ion_temperature': ['T_i', 3, ['time', 'lon_h', 'lat_h',
                                                'height'], 'K'],
                 'height_d': ['H_lev', 4, ['time', 'lon_d', 'lat_d', 'lev'],
                              'm'],
                 'height_n': ['H_ilev', 5, ['time', 'lon_n', 'lat_n', 'ilev'],
                              'm'],
                 'height': ['H_ilev', 5, ['time', 'lon_n', 'lat_n', 'ilev'],
                              'm'],
                 'meridional_neutral_wind': ['Vn_lat', 6, ['time', 'lon_n',
                                                           'lat_n', 'ilev'],
                                             'm/s'],
                 'zonal_neutral_wind': ['Vn_lon', 7, ['time', 'lon_n', 'lat_n',
                                                      'ilev'], 'm/s'],
                 'vertical_neutral_wind': ['Vn_H', 8, ['time', 'lon_n',
                                                       'lat_n', 'ilev'],
                                           'm/s'],
                 'neutral_temperature': ['T_n', 9, ['time', 'lon_n', 'lat_n',
                                                    'ilev'], 'K'],
                 'mean_molecular_mass': ['Rmt', 10, ['time', 'lon_n', 'lat_n',
                                                     'ilev'], 'amu'],
                 'electron_density': ['N_e', 11, ['time', 'lon_h', 'lat_h',
                                                  'height'], '1/m**3'],
                 # 'neutral_density': ['N_n',11,['time','lon_n','lat_n',
                 # 'ilev'],'1/m**3'],
                 'solar_heating': ['Q_Solar', 12, ['time', 'lon_n', 'lat_n',
                                                   'ilev'], 'J/kg/s'],
                 'joule_heating': ['Q_Joule', 13, ['time', 'lon_n', 'lat_n',
                                                   'ilev'], 'J/kg/s'],
                 'radiation_heat_cool': ['Q_radiation', 14, ['time', 'lon_n',
                                                             'lat_n', 'ilev'],
                                         'J/kg/s'],
                 'atomic_oxygen_density': ['N_O', 15, ['time', 'lon_n',
                                                       'lat_n', 'ilev'],
                                           '1/m**3'],
                 'molecular_oxygen_density': ['N_O2', 16, ['time', 'lon_n',
                                                           'lat_n', 'ilev'],
                                              '1/m**3'],
                 'molecular_nitrogen_density': ['N_N2', 17, ['time', 'lon_n',
                                                             'lat_n', 'ilev'],
                                                '1/m**3'],
                 'nitric_oxide_density': ['N_NO', 18, ['time', 'lon_n',
                                                       'lat_n', 'ilev'],
                                          '1/m**3'],
                 'nitric_oxide_ion_density': ['N_NOplus', 19,
                                              ['time', 'lon_n', 'lat_n',
                                               'ilev'], '1/m**3'],
                 'molecular_nitrogen_ion_density': ['N_N2plus', 20,
                                                    ['time', 'lon_n', 'lat_n',
                                                     'ilev'], '1/m**3'],
                 'molecular_oxygen_ion_density': ['N_O2plus', 21,
                                                  ['time', 'lon_n', 'lat_n',
                                                   'ilev'], '1/m**3'],
                 'atomic_nitrogen_ion_density': ['N_Nplus', 22,
                                                 ['time', 'lon_n', 'lat_n',
                                                  'ilev'], '1/m**3'],
                 'atomic_oxygen_ion_density': ['N_Oplus', 23,
                                               ['time', 'lon_h', 'lat_h',
                                                'height'], '1/m**3'],
                 'atomic_hydrogen_ion_density': ['N_Hplus', 24,
                                                 ['time', 'lon_h', 'lat_h',
                                                  'height'], '1/m**3'],
                 'pedersen_conductivity': ['Sigma_P', 25,
                                           ['time', 'lon_n', 'lat_n', 'ilev'],
                                           'S/m'],
                 'hall_conductivity': ['Sigma_H', 26,
                                       ['time', 'lon_n', 'lat_n', 'ilev'],
                                       'S/m'],
                 'zonal_ion_velocity': ['Vi_lon', 27,
                                        ['time', 'lon_n', 'lat_n', 'ilev'],
                                        'm/s'],
                 'meridional_ion_velocity': ['Vi_lat', 28,
                                             ['time', 'lon_n', 'lat_n',
                                              'ilev'], 'm/s'],

                 # start 3D variables
                 'height_integrated_joule_heating': ['W_Joule', 29,
                                                     ['time', 'lon_n',
                                                      'lat_n'], 'W/m**2'],
                 'energy_influx': ['Eflux_precip', 30,
                                   ['time', 'lon_n', 'lat_n'], 'W/m**2'],
                 'mean_energy': ['Eavg_precip', 31,
                                 ['time', 'lon_n', 'lat_n'], 'keV'],
                 # '10**16/m**2' for TEC
                 'total_electron_content': ['TEC', 32,
                                            ['time', 'lon_n', 'lat_n'],
                                            '1/m**2'],
                 'theta_electric_field_at_140km': ['E_theta140km', 33,
                                                   ['time', 'Elon', 'Elat'],
                                                   'V/m'],
                 'lambda_electric_field_at_140km': ['E_lambda140km', 34,
                                                    ['time', 'Elon', 'Elat'],
                                                    'V/m'],
                 'theta_electric_field_at_300km': ['E_theta300km', 35,
                                                   ['time', 'Elon', 'Elat'],
                                                   'V/m'],
                 'lambda_electric_field_at_300km': ['E_lambda300km', 36,
                                                    ['time', 'Elon', 'Elat'],
                                                    'V/m']}


def convert_all(file_prefix):
    '''Get all data from last timestep of previous file. Finds first file
    and converts them all in order.'''

    # find first file
    print(file_prefix)
    file_date = basename(file_prefix)  # YYYY-MM-DD
    file_dir = file_prefix.split(file_date)[0]
    earlier_date = datetime.strptime(file_date, '%Y-%m-%d') - timedelta(days=1)
    new_file_date = datetime.strftime(earlier_date, '%Y-%m-%d')
    filetype_list = ['-plot-density.nc', '-plot-height.nc', '-plot-neutral.nc']
    files = [file_dir+new_file_date+filetype for filetype in filetype_list]
    while sum([isfile(file) for file in files]) > 1:
        good_files = files
        earlier_date -= timedelta(days=1)
        new_file_date = datetime.strftime(earlier_date, '%Y-%m-%d')
        files = [file_dir+new_file_date+filetype for filetype in filetype_list]
    else:
        good_files = [file_dir+file_date+filetype for filetype in
                      filetype_list]

    # convert first file with an empty first time step
    new_file, in_ddict, in_hdict, in_ndict, last_time = \
        ctipe_combine_files(good_files, verbose=False)

    # convert remaining file groups while adding a timestep from the previous
    first_file_prefix = good_files[0].split('-plot-density.nc')[0]
    file_date = basename(first_file_prefix)
    later_date = datetime.strptime(file_date, '%Y-%m-%d') + timedelta(days=1)
    new_file_date = datetime.strftime(later_date, '%Y-%m-%d')
    files = [file_dir+new_file_date+filetype for filetype in filetype_list]
    while sum([isfile(file) for file in files]) > 1:
        new_file, in_ddict, in_hdict, in_ndict, last_time = \
            ctipe_combine_files(files, in_ddict=in_ddict, in_hdict=in_hdict,
                                in_ndict=in_ndict, earlier_time=last_time,
                                verbose=False)
        later_date += timedelta(days=1)
        new_file_date = datetime.strftime(later_date, '%Y-%m-%d')
        files = [file_dir+new_file_date+filetype for filetype in filetype_list]
    return file_prefix + '.nc'


def ctipe_wrap_variables(var_dict, variable_name, lon):
    '''wrap variables in longitude and transpose as needed for ctipe model
    output'''

    lon_le180 = where(lon <= 180)[0]
    lon_ge180 = where(lon >= 180)[0]  # repeat 180 for -180 values
    if 'electric_field' in variable_name:
        # CTIPe efield variables from neutral file do not need to be wrapped
        # but need to be transposed from (time,lon,lat) to (time,lat,lon)
        # new_variable = np.transpose(variable,[0,2,1])
        pass  # want in (time, lon, lat)
    # 3D variable, wrap and shift in longitude then transpose
    elif len(var_dict['data'].shape) == 3:
        shape_list = list(var_dict['data'].shape)  # time, lat, lon
        shape_list[2] += 1  # need one more place in longitude
        tmp_arr = zeros(shape_list)  # array to set-up wrapped data in
        tmp_arr[:, :, :len(lon_ge180)] = var_dict['data'][:, :, lon_ge180]
        tmp_arr[:, :, len(lon_ge180):] = var_dict['data'][:, :, lon_le180]
        # (t,lat,lon) -> (t,lon,lat)
        var_dict['data'] = transpose(tmp_arr, (0, 2, 1))
        var_dict['size'] = (shape_list[0], shape_list[2], shape_list[1])
    elif len(var_dict['data'].shape) == 4:  # 4D variables
        shape_list = list(var_dict['data'].shape)  # time, lat, lon, height
        shape_list[3] += 1  # need one more place in longitude
        tmp_arr = zeros(shape_list)  # array to set-up wrapped data in
        tmp_arr[:, :, :, :len(lon_ge180)] = var_dict['data'][:, :, :,
                                                             lon_ge180]
        tmp_arr[:, :, :, len(lon_ge180):] = var_dict['data'][:, :, :,
                                                             lon_le180]
        # (t,h,lat,lon) -> (t,lon,lat,h)
        var_dict['data'] = transpose(tmp_arr, (0, 3, 2, 1))
        var_dict['size'] = (shape_list[0], shape_list[3], shape_list[2],
                            shape_list[1])
    return var_dict


def ctipe_combine_files(files, in_ddict=None, in_hdict=None, in_ndict=None,
                        earlier_time=None, verbose=False):
    '''Combine data from 3 files, wrapping in longitude and transposing as
    necessary.'''

    tic = perf_counter()
    dim_dict = {}
    # determine file names of group
    filename_density, filename_height, filename_neutral = files
    file_prefix = filename_density.split('-plot-density.nc')[0]

    # open data files, retrieve data, and close. wrap longitude
    if isfile(filename_density):
        ctipe_density = Dataset(filename_density)
        d_dict = {key: {'data': array(ctipe_density.variables[key]),
                        'datatype': ctipe_density.variables[key].datatype,
                        'size': ctipe_density.variables[key].size}
                  for key in ctipe_density.variables.keys()}
        ctipe_density.close()

        # wrap and shift in longitude
        for key in d_dict.keys():
            if key not in ['time', 'lat', 'lon', 'elat', 'elon', 'ht', 'plev']:
                d_dict[key] = ctipe_wrap_variables(d_dict[key], key,
                                                   d_dict['lon']['data'])
        d_dict['lon']['data'] = append(d_dict['lon']['data'], 360.) - 180.
        d_dict['lon']['size'] += 1

        # collect dimensions data, assuming time and ilev are all the same
        # assuming lon and lat in diff files might be different
        dim_dict['lat_d'] = d_dict['lat']
        dim_dict['lon_d'] = d_dict['lon']
        dim_dict['lev'] = d_dict['plev']
        dim_dict['time'] = d_dict['time']

        # adjust 'height' variable names to better distinguish
        d_dict['height_d'] = d_dict['height']
        del d_dict['height']

        # remove keys from file dictionaries for coordinates
        for key in ['time', 'lat', 'lon', 'elat', 'elon', 'ht', 'plev']:
            if key in d_dict.keys():
                del d_dict[key]

        # add time to beginning of arrays if earlier file existed
        if in_ddict != None:
            for key in in_ddict.keys():
                d_dict[key]['data'] = insert(d_dict[key]['data'], 0,
                                             in_ddict[key]['data'], axis=0)
                d_dict[key]['size'] = d_dict[key]['data'].shape

        # Calculate median height for each pressure level and store
        dim_dict['km_lev'] = {'size': len(dim_dict['lev']['data']),
                              'datatype': float,
                              'data': median(d_dict['height_d']['data'],
                                             axis=[0,1,2])/1000.}

    if isfile(filename_height):
        ctipe_height = Dataset(filename_height)  # in meters
        h_dict = {key: {'data': array(ctipe_height.variables[key]),
                        'datatype': ctipe_height.variables[key].datatype,
                        'size': ctipe_height.variables[key].size}
                  for key in ctipe_height.variables.keys()}
        ctipe_height.close()

        # wrap and shift in longitude
        for key in h_dict.keys():
            if key not in ['time', 'lat', 'lon', 'elat', 'elon', 'ht', 'plev']:
                h_dict[key] = ctipe_wrap_variables(h_dict[key], key,
                                                   h_dict['lon']['data'])
        h_dict['lon']['data'] = append(h_dict['lon']['data'], 360.) - 180.
        h_dict['lon']['size'] += 1

        # collect dimensions data, assuming time and ilev are all the same
        # assuming lon and lat in diff files might be different
        dim_dict['lat_h'] = h_dict['lat']
        dim_dict['lon_h'] = h_dict['lon']
        dim_dict['height'] = h_dict['ht']  # in km

        # collect time if prev file DNE
        if not isfile(filename_density):
            dim_dict['time'] = h_dict['time']

        # remove keys from file dictionaries for coordinates
        for key in ['time', 'lat', 'lon', 'elat', 'elon', 'ht', 'plev']:
            if key in h_dict.keys():
                del h_dict[key]

        # add time to beginning of arrays if earlier file existed
        if in_hdict != None:
            for key in in_hdict.keys():
                h_dict[key]['data'] = insert(h_dict[key]['data'], 0,
                                             in_hdict[key]['data'], axis=0)
                h_dict[key]['size'] = h_dict[key]['data'].shape

    if isfile(filename_neutral):
        ctipe_neutral = Dataset(filename_neutral)
        n_dict = {key: {'data': array(ctipe_neutral.variables[key]),
                        'datatype': ctipe_neutral.variables[key].datatype,
                        'size': ctipe_neutral.variables[key].size}
                  for key in ctipe_neutral.variables.keys()}
        ctipe_neutral.close()

        # wrap and shift in longitude
        for key in n_dict.keys():
            if key not in ['time', 'lat', 'lon', 'elat', 'elon', 'ht', 'plev']:
                if 'electric_field' in key:
                    tmp_lon = n_dict['elon']['data']
                else:
                    tmp_lon = n_dict['lon']['data']
                n_dict[key] = ctipe_wrap_variables(n_dict[key], key, tmp_lon)
        n_dict['lon']['data'] = append(n_dict['lon']['data'], 360.) - 180.
        n_dict['lon']['size'] += 1
        n_dict['elon']['data'] -= 180.

        # collect dimensions data, assuming time and ilev are all the same
        # assuming lon and lat in diff files might be different
        dim_dict['lat_n'] = n_dict['lat']
        dim_dict['lon_n'] = n_dict['lon']
        dim_dict['ilev'] = n_dict['plev']
        dim_dict['Elat'] = n_dict['elat']
        dim_dict['Elon'] = n_dict['elon']

        # collect time if prev files DNE
        if not isfile(filename_density) and not isfile(filename_height):
            dim_dict['time'] = n_dict['time']

        # remove keys from file dictionaries for coordinates
        for key in ['time', 'lat', 'lon', 'elat', 'elon', 'ht', 'plev']:
            if key in n_dict.keys():
                del n_dict[key]

        # add time to beginning of arrays if earlier file existed
        if in_ndict != None:
            for key in in_ndict.keys():
                n_dict[key]['data'] = insert(n_dict[key]['data'], 0,
                                             in_ndict[key]['data'], axis=0)
                n_dict[key]['size'] = n_dict[key]['data'].shape

        # Calculate median height for each pressure level and store
        # neutral height data sometimes has different names
        if 'height' in n_dict.keys():
            hkey = 'height'
        elif 'height_n' in n_dict.keys():
            hkey = 'height_n'
        dim_dict['km_ilev'] = {'size': len(dim_dict['ilev']['data']),
                               'datatype': float, 'data':
                                   median(n_dict[hkey]['data'],
                                                     axis=[0,1,2])/1000.}

    # add time value to time dimension if previous file existed
    if in_ddict != None or in_hdict != None or in_ndict != None:
        dim_dict['time']['size'] += 1
        dim_dict['time']['data'] = insert(dim_dict['time']['data'], 0,
                                          earlier_time)
        print('Earlier time added')

    # initialize single output file
    data_out = Dataset(file_prefix+'.nc', 'w', format='NETCDF4')
    data_out.model = 'CTIPe'
    # csv list of files
    data_out.file = ''.join([f+',' for f in [filename_density, filename_height,
                                             filename_neutral]
                             if isfile(f)]).strip(',')

    # store dimensions
    for dim in dim_dict.keys():
        if verbose:
            print(dim)
        new_dim = data_out.createDimension(dim, dim_dict[dim]['size'])
        new_var = data_out.createVariable(dim, dim_dict[dim]['datatype'],
                                          tuple((dim, )))
        new_var[:] = dim_dict[dim]['data']
    last_time = dim_dict['time']['data'][-1]
    if verbose:
        print('Dimensions complete.\n')

    # add variable data from density file to output file
    if isfile(filename_density):
        key_list = list(d_dict.keys())
        for key in key_list:
            if key in ['ZMAG', 'mean_molecular_mass'] or key not in\
                    file_varnames.keys():
                del d_dict[key]
                continue  # using Rmt in neutral file
            if verbose:
                print(key, file_varnames[key][2])
            new_var = data_out.createVariable(file_varnames[key][0],
                                              d_dict[key]['datatype'],
                                              tuple(file_varnames[key][2]))
            new_var[:] = d_dict[key]['data']
            tmp = d_dict[key]['data'][-1]
            d_dict[key]['data'] = tmp  # save last time step for next file
        if verbose:
            print('Density file complete.\n')
    else:
        d_dict = None

    # add variable data from height file to output file
    if isfile(filename_height):
        key_list = list(h_dict.keys())
        for key in key_list:
            if key == 'ZMAG' or key not in file_varnames.keys():
                del h_dict[key]
                continue  # no other variables depend on ZMAG, so ignore
            if verbose:
                print(key, file_varnames[key][2])
            new_var = data_out.createVariable(file_varnames[key][0],
                                              h_dict[key]['datatype'],
                                              tuple(file_varnames[key][2]))
            new_var[:] = h_dict[key]['data']
            tmp = h_dict[key]['data'][-1]
            h_dict[key]['data'] = tmp  # save last time step for next file
        if verbose:
            print('Height file complete.\n')
    else:
        h_dict = None

    # add variable data from neutral file to output file
    if isfile(filename_neutral):
        key_list = list(n_dict.keys())
        for key in key_list:
            if key in ['ZMAG', 'electron_density', 'atomic_oxygen_ion_density',
                       'atomic_hydrogen_ion_density'] or key not in\
                    file_varnames.keys():
                if key == 'electron_density' and not isfile(filename_height):
                    tmp = ['time', 'lon_n', 'lat_n', 'ilev']
                    new_var = data_out.createVariable(file_varnames[key][0],
                                                      n_dict[key]['datatype'],
                                                      tuple(tmp))
                    new_var[:] = n_dict[key]['data']
                    tmp = n_dict[key]['data'][-1]
                    n_dict[key]['data'] = tmp  # save last time step for next file
                else:
                    del n_dict[key]
                continue
            else:
                new_var = data_out.createVariable(file_varnames[key][0],
                                                  n_dict[key]['datatype'],
                                                  tuple(file_varnames[key][2]))
                new_var[:] = n_dict[key]['data']
                tmp = n_dict[key]['data'][-1]
                n_dict[key]['data'] = tmp  # save last time step for next file
            if verbose:
                print(key, file_varnames[key][2])
        if verbose:
            print('Neutral file complete.\n')
    else:
        n_dict = None

    # close file and return data for last time step
    print(f"Data for {file_prefix} converted in {perf_counter()-tic:.6f}s.")
    data_out.close()
    return file_prefix+'.nc', d_dict, h_dict, n_dict, last_time


if __name__ == '__main__':
    # define file names (input and output)
    file_dir = 'C:/Users/rringuet/Kamodo_Data/CTIPe/Data/'
    file_prefix = file_dir+'2015-03-18'
    new_filename = ctipe_combine_files(file_prefix)
