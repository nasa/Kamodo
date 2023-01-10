# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 11:57:13 2022

@author: rringuet
"""
from csv import reader
from datetime import datetime, timezone
from netCDF4 import Dataset
import numpy as np
from glob import glob
from os.path import basename, isfile


def convert_all(file_dir):
    '''Converting all files to netCDF4.'''

    files = sorted(glob(file_dir + '*.txt'))  # want YYYYMMDD
    dates = np.unique([basename(f).split('.')[0][-8:] for f in files])
    # print(dates)
    file_patterns = np.unique([basename(f).split('.')[0][:-9] for f in files])
    # file(s) could be from either hemisphere, but only one file naming pattern
    new_pattern = np.unique([''.join(*[basename(p).split('_N')]) if '_N' in p
                             else ''.join(*[basename(p).split('_S')])
                             for p in file_patterns])[0]  # only one pattern
    for d in dates:
        new_file = file_dir + new_pattern + '_' + d + '.nc'
        if not isfile(new_file):  # only convert if converted file not found
            to_CDF(new_file[:-3])  # slice off '.nc'
    return True


def adelphi_data_hemisphere(file_name):
    '''Read in data from given file and return in a dictionary of arrays.'''

    # read in data from file
    file_obj = open(file_name, 'r')
    csv_reader = reader(file_obj, delimiter='\n')  # splits file by \n char

    # pull in date and time
    line = next(csv_reader)
    data = {'Time': [float(line[0][9:].strip().split(' ')[0])]}  # in hours
    year, month, day = int(line[0][:4]), int(line[0][4:6]), int(line[0][6:8])
    dt_str = datetime(year, month, day, 0, 0, 0, tzinfo=timezone.utc
                      ).strftime('%Y-%m-%d %H:%M:%S')

    # set up data structure with next two lines
    line = next(csv_reader)
    data['MLT'] = [float(line[0].split('=')[1].split('Hour')[0].strip())]
    line = next(csv_reader)
    keys = [key for key in line[0].split(' ') if key != '']  # MLAT, PED, etc
    for key in keys:
        data[key] = [[[]]]

    # read in data
    time_idx, lon_idx = 0, 0
    for line in csv_reader:
        if line[0][:4] == dt_str[:4]:  # retrieve time value
            data['Time'].append(float(line[0][9:].strip().split(' ')[0]))
            time_idx += 1
            lon_idx = -1
            for i in range(len(keys)):
                data[keys[i]].append([])
        elif line[0][:3] == 'MLT':  # retrive coordinate value
            data['MLT'].append(float(line[0].split('=')[1].split('Hour')[0].
                                     strip()))
            lon_idx += 1
            for i in range(len(keys)):
                data[keys[i]][time_idx].append([])
        elif line[0][:4] == 'MLAT':  # skip repeat of keys
            continue
        else:  # read in data values
            # check for numbers running against each other, insert space if yes
            items = line[0][1:].split('-')  # avoid preceding - for S hemispher
            if len(items) > 1:
                # print(items)
                for i in range(len(items[1:])):
                    if items[i][-1] != ' ':
                        items[i] += ' '
                line[0] = line[0][0] + ''.join(['-'+item for item in items
                                                ])[1:]

            # put data in correct list
            values = [float(value) for value in line[0].split(' ') if
                      value != '']
            for i in range(len(keys)):
                data[keys[i]][time_idx][lon_idx].append(values[i])
    file_obj.close()

    # convert to arrays an reshape variable data
    data['Time'] = np.array(data['Time'])  # time in hours
    data['MLT'] = np.unique(np.array(data['MLT'])) * 15.  # magnetic longitude
    data['MLAT'] = np.array(data['MLAT'][0][0])  # magnetic latitude
    # in the data MLAT is fastest changing, then MLT, then Time
    for key in keys[1:]:
        data[key] = np.array(data[key])
    return data, dt_str


def vector_average3D(top, shape_list, varname, latval, lon):
    '''find vector average at pole for array with shape (time, lon,
    height)'''

    # CAN NOT TEST BECAUSE NUMBERS ARE ALL ZEROS!!!
    # find net x and y components, final array shapes are (t, lon, ht)
    # reverse and .T?
    lon_arr = np.broadcast_to(lon[1:], (shape_list[0], shape_list[1]-1))
    # need to put 'old' shape at end in broadcast_to call ^
    xval = np.sum(top*np.cos((lon_arr+180.)*np.pi/180.), axis=1)  # same shape
    yval = np.sum(top*np.sin((lon_arr+180.)*np.pi/180.), axis=1)  # as t and z
    xarr = np.broadcast_to(xval, (shape_list[1]-1, shape_list[0])).T
    yarr = np.broadcast_to(yval, (shape_list[1]-1, shape_list[0])).T

    # convert to proper unit vector (see wiki on spherical coordinates)
    if 'EAST' in varname:
        # Zonal/east components -> convert to psi_hat vector (lon)
        # -xsin(psi)+ycos(psi), psi = longitude (0 to 360)
        new_top = -xarr*np.sin((lon_arr+180.)*np.pi/180.) +\
            yarr*np.cos((lon_arr+180.)*np.pi/180.)
    elif 'NORTH' in varname:
        # meridional/north -> convert to theta_hat vector (latitude)
        # xcos(psi)cos(theta)+ysin(psi)cos(theta)
        # sin(theta) is always zero at the poles
        # theta = latitude (0 to 180), psi = longitude (0 to 360)
        new_top = xarr * np.cos((lon_arr+180.) * np.pi/180.) *\
            np.cos((90.-latval)*np.pi/180.) + yarr *\
            np.sin((lon_arr+180.) * np.pi/180.) *\
            np.cos((90.-latval) * np.pi/180.)
    return new_top


def combine_hemispheres(north_file, south_file):
    '''Given the filenames for the northern and southern hemispheres, combine
    the data into a single dictionary and return.'''

    data_N, date_str = adelphi_data_hemisphere(north_file)
    data_S, date_str = adelphi_data_hemisphere(south_file)
    key_list = [key for key in data_N.keys() if key not in ['Time', 'MLT',
                                                            'MLAT']]
    # for key in key_list:
    #   print(key, np.max(data_S[key][:, :, 0]), np.max(data_S[key][:, :, -1]),
    #          np.max(data_N[key][:, :, 0]), np.max(data_N[key][:, :, -1]))

    # Combine coordinate grids, carefully piecing together the latitudes
    # ********* CHECK WITH SCIENTISTS ABOUT THE 180 DEGREE SHIFT **************
    len_lat = len(data_N['MLAT'])
    data = {'Time': data_N['Time'], 'Lon': np.append(data_N['MLT'], 360.)-180.}
    data['Lat'] = np.zeros(len_lat*2+6)  # need poles and NaN rows
    data['Lat'][0], data['Lat'][-1] = -90., 90.  # add pole latitude values
    data['Lat'][1:len_lat+1] = np.flip(data_S['MLAT'])
    # add buffer rows for S hemisphere
    data['Lat'][len_lat+1] = data['Lat'][len_lat] +\
        abs(np.diff(data_S['MLAT'])[0])/10.
    data['Lat'][len_lat+2] = data['Lat'][len_lat] +\
        abs(np.diff(data_S['MLAT'])[0])
    data['Lat'][len_lat+5:-1] = data_N['MLAT']  # add in N hemi lats
    data['Lat'][len_lat+4] = data['Lat'][len_lat+5] -\
        np.diff(data_N['MLAT'])[0]/10.  # add buffer rows for N hemisphere
    data['Lat'][len_lat+3] = data['Lat'][len_lat+5] -\
        np.diff(data_N['MLAT'])[0]
    # print('Latitude:', data['Lat'])

    # combine data from different hemispheres
    new_shape = (len(data['Time']), len(data['Lon']), len(data['Lat']))
    for key in key_list:
        # pull in data and set buffer rows on equator side for both hemispheres
        data[key] = np.zeros(new_shape)
        data[key][:, 1:, 1:len_lat+1] = np.flip(data_S[key], axis=2)  # S data
        data[key][:, 1:, len_lat+1] = data[key][:, 1:, len_lat]  # value buffer
        data[key][:, 1:, len_lat+2:len_lat+4] = \
            np.tile(np.NaN, (new_shape[0], new_shape[1]-1, 2))  # NaN buffer
        data[key][:, 1:, len_lat+4] = data[key][:, 1:, len_lat+5]  # value buff
        data[key][:, 1:, len_lat+5:-1] = data_N[key]  # N hemisphere data

        # wrap in longitude
        data[key][:, 0, :] = data[key][:, -1, :]

        ''' The # of rows cut off is different for each variable, so avoiding
        # slice off zeros at each pole before averaging and functionalizing
        # find first latitude row that is nonzero, only near poles
        # this purposefully ignores the buffer rows near the equator
        # this is not in the file converter due to the dynamic nature of
        #   the latitude grid.  ????
        # south pole at beginning
        SP_idx, slice_idx = 1, []
        zero_check = np.count_nonzero(data[key][:, :, SP_idx])
        while zero_check == 0:
            slice_idx.append(SP_idx)
            SP_idx += 1
            zero_check = np.count_nonzero(data[key][:, :, SP_idx])
        # north pole at the end
        NP_idx = -2
        zero_check = np.count_nonzero(data[key][:, :, NP_idx])
        while zero_check == 0:
            slice_idx.append(NP_idx)
            NP_idx -= 1
            zero_check = np.count_nonzero(data[key][:, :, NP_idx])
        # remove 'extra' latitude values from coordinate grid and from data
        #data['Lat'] = np.delete(data['Lat'], slice_idx)
        print(key, len(slice_idx), slice_idx)
        data[key] = np.delete(data[key], slice_idx, axis=2)
        '''

        # perform scalar averaging at the poles
        # new_shape =(data['Time'].shape, data['Lon'].shape, data['Lat'].shape)
        if key not in ['E_east', 'E_north', 'J_east', 'J_north']:
            # south pole
            top = np.mean(data[key][:, 1:, 1], axis=1)  # same shape as time
            data[key][:, 1:, 0] = np.broadcast_to(
                top, (new_shape[1]-1, new_shape[0])).T
            # north pole
            top = np.mean(data[key][:, 1:, -2], axis=1)  # same shape as time
            data[key][:, 1:, -1] = np.broadcast_to(
                top, (new_shape[1]-1, new_shape[0])).T
        else:  # perform vector averaging at the poles
            # calculate net vector magnitude for south pole
            data[key][:, 1:, 0] = \
                vector_average3D(data[key][:, 1:, 1], list(data[key].shape),
                                 key, data['Lat'][0], data['Lon'])
            # repeat for north pole
            data[key][:, 1:, -1] = \
                vector_average3D(data[key][:, 1:, -2], list(data[key].shape),
                                 key, data['Lat'][-1], data['Lon'])

    # data wrangling complete, return data dictionary
    return data, date_str


def to_CDF(file_prefix):
    '''Given the file_prefix of form ADELPHI_2D_MAG_YYYMMDD, find and combine
    data from both hemispheres into a netCDF4 file.'''
    from time import perf_counter

    print('Converted data file not found. Converting files with ' +
          f'{file_prefix} prefix to a netCDF4 file.')  # , end="")
    north_file = file_prefix[:-13]+'_N'+file_prefix[-13:]+'.txt'
    south_file = file_prefix[:-13]+'_S'+file_prefix[-13:]+'.txt'
    time0 = perf_counter()
    data, date_str = combine_hemispheres(north_file, south_file)  # ~ 18-19 s
    time1 = perf_counter()
    print('Time to read and combine the data sets:', time1-time0)

    # file_prefix = file_dir + 'ADELPHI_2D_MAG_20130317'
    cdf_filename = file_prefix + '.nc'
    data_out = Dataset(cdf_filename, 'w', format='NETCDF4')
    data_out.file = north_file + ',' + south_file
    data_out.model = 'ADELPHI'
    data_out.filedate = date_str

    # establish coordinates (time, lon, lat)
    # time: create dimension and variable
    new_dim = data_out.createDimension('time', data['Time'].size)
    new_var = data_out.createVariable('time', np.float32, 'time')
    new_var[:] = data['Time']
    new_var.units = 'hr'
    # lon: create dimension and variable
    new_dim = data_out.createDimension('lon', data['Lon'].size)
    new_var = data_out.createVariable('lon', np.float32, 'lon')
    new_var[:] = data['Lon']  # store data for dimension in variable
    new_var.units = 'deg'
    # lat: create dimension and variable
    new_dim = data_out.createDimension('lat', data['Lat'].size)
    new_var = data_out.createVariable('lat', np.float32, 'lat')
    new_var[:] = data['Lat']  # store data for dimension in variable
    new_var.units = 'deg'

    # copy over variables to file
    var_list = [key for key in data.keys() if key not in
                ['Time', 'Lon', 'Lat']]
    for variable_name in var_list:
        new_var = data_out.createVariable(variable_name, np.float32,
                                          ('time', 'lon', 'lat'))
        new_var[:] = data[variable_name]  # store data in variable
        new_var.units = 'V'

    # close file
    data_out.close()
    print('Time to create netCDF4 file:', perf_counter()-time0)
    return cdf_filename
