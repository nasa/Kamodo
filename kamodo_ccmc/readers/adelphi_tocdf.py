# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 11:57:13 2022

@author: rringuet
"""
from csv import reader
from datetime import datetime, timezone
from netCDF4 import Dataset
import numpy as np
from os.path import basename


def adelphi_data_hemisphere(file_name):
    '''Read in data from given file and return in a dictionary of arrays.'''

    # read in data from file
    file_obj = open(file_name, 'r')
    csv_reader = reader(file_obj, delimiter='\n') # splits file by \n character

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
    keys = [key for key in line[0].split(' ') if key != '']
    for key in keys:
        data[key] = []

    # read in data
    for line in csv_reader:
        if line[0][:4] == dt_str[:4]:  # retrieve time value
            data['Time'].append(float(line[0][9:].strip().split(' ')[0]))
        elif line[0][:3] == 'MLT':  # retrive coordinate value
            data['MLT'].append(float(line[0].split('=')[1].split('Hour')[0].
                                     strip()))
        elif line[0][:4] == 'MLAT':  # skip repeat of keys
            continue
        else:  # read in data values
            # check for numbers running against each other, insert space if yes
            items = line[0].split('-')
            for i in range(len(items[1:])):
                if items[i][-1] != ' ':
                    items[i] += ' '
            line[0] = ''.join(['-'+item for item in items])[1:]
    
            # put data in correct list
            values = [float(value) for value in line[0].split(' ') if
                      value != '']
            for i in range(len(keys)):
                data[keys[i]].append(values[i])
    file_obj.close()

    # convert to arrays an reshape variable data
    data['Time'] = np.array(data['Time'])  # time in hours
    data['MLT'] = np.unique(data['MLT']) * 15.  # magnetic longitude
    data['MLAT'] = np.unique(data['MLAT'])  # magnetic latitude
    array_shape = (len(data['Time']), len(data['MLT']), len(data['MLAT']))
    for key in keys[1:]:
        data[key] = np.reshape(data[key], array_shape)

    return data, dt_str


def combine_hemispheres(north_file, south_file):
    '''Given the filenames for the northern and southern hemispheres, combine
    the data into a single dictionary and return.'''

    data_N, date_str = adelphi_data_hemisphere(north_file)
    data_S, date_str = adelphi_data_hemisphere(south_file)
    
    # Test hemispheres for identical grids, REMOVE LATER
    for key in ['Time', 'MLT', 'MLAT']:
        testbyitem = sum(data_N[key] == data_S[key])
        test = (testbyitem == len(data_S[key])) & (testbyitem ==
                                                   len(data_N[key]))
        if not test:
            print(key, testbyitem, len(data_N[key]), len(data_S[key]))
    
    # Combine coordinate grids, carefully piecing together the latitudes
    # ********* CHECK WITH SCIENTISTS ABOUT THE 180 DEGREE SHIFT **************
    len_lon = len(data_N['MLAT'])
    data = {'Time': data_N['Time'], 'Lon': np.append(data_N['MLT'], 360.)-180.}
    data['Lat'] = np.zeros(len_lon*2+6)  # need poles and NaN rows
    data['Lat'][0], data['Lat'][-1] = -90., 90.  # add pole latitude values
    data['Lat'][1:len_lon+1] = np.flip(data_S['MLAT']) * -1  # make S lats neg
    # add buffer rows for S hemisphere
    data['Lat'][len_lon+1] = data['Lat'][len_lon] +\
        np.diff(data_S['MLAT'])[0]/10.
    data['Lat'][len_lon+2] = data['Lat'][len_lon] + np.diff(data_S['MLAT'])[0]
    data['Lat'][len_lon+5:-1] = data_N['MLAT']  # add in N hemi lats
    data['Lat'][len_lon+4] = data['Lat'][len_lon+5] -\
        np.diff(data_N['MLAT'])[0]/10. # add buffer rows for N hemisphere
    data['Lat'][len_lon+3] = data['Lat'][len_lon+5] -\
        np.diff(data_N['MLAT'])[0]
        
    # combine data from different hemispheres
    new_shape = (len(data['Time']), len(data['Lon']), len(data['Lat']))
    key_list = [key for key in data_N.keys() if key not in ['Time', 'MLT',
                                                            'MLAT']]
    for key in key_list:
        # pull in data and set buffer rows
        data[key] = np.zeros(new_shape)
        data[key][:, 1:, 1:len_lon+1] = np.flip(data_S[key], axis=2)  # S data
        data[key][:, 1:, len_lon+1] = data[key][:, 1:, len_lon]  # value buffer
        data[key][:, 1:, len_lon+2:len_lon+4] = \
            np.tile(np.NaN, (new_shape[0], new_shape[1]-1, 2))  # NaN buffer
        data[key][:, 1:, len_lon+5:-1] = data_N[key]  # N hemisphere data
        data[key][:, 1:, len_lon+4] = data[key][:, 1:, len_lon+5]  # value buff
    
        # perform scalar averaging at the poles, SP at beginning then NP at end
        top = np.mean(data[key][:, 1:, 1], axis=1)  # same shape as time axis
        data[key][:, 1:, 0] = np.broadcast_to(top, (new_shape[1]-1,
                                                    new_shape[0])).T
        top = np.mean(data[key][:, 1:, -2], axis=1)  # same shape as time axis
        data[key][:, 1:, -1] = np.broadcast_to(top, (new_shape[1]-1,
                                                     new_shape[0])).T
    
        # wrap in longitude
        data[key][:, 0, :] = data[key][:, -1, :]
    
    # data wrangling complete, return data dictionary
    return data, date_str

# files
#file_dir = 'C:/Users/rringuet/Kamodo_Data/ADELPHI/Robert_Robinson_20201023_IT_1/'
#north_file = file_dir + 'ADELPHI_2D_N_MAG_20130317.txt'
#south_file = file_dir + 'ADELPHI_2D_S_MAG_20130317.txt'
#file_prefix = file_dir + 'ADELPHI_2D_MAG_20130317'
def to_CDF(file_prefix):
    '''Given the file_prefix of form ADELPHI_2D_MAG_YYYMMDD, find and combine
    data from both hemispheres into a netCDF4 file.'''
    from time import perf_counter

    print('Converted data file not found. Converting files with ' +
          f'{file_prefix} prefix to a netCDF4 file.') # , end="")
    north_file = file_prefix[:-13]+'_N'+file_prefix[-13:]+'.txt'
    south_file = file_prefix[:-13]+'_S'+file_prefix[-13:]+'.txt'
    time0 = perf_counter()
    data, date_str = combine_hemispheres(north_file, south_file)  # takes about 18-19 seconds
    time1 = perf_counter()
    print('Time to read and combine the data sets:', time1-time0)
    
    #file_prefix = file_dir + 'ADELPHI_2D_MAG_20130317'
    cdf_filename = file_prefix + '.nc'
    data_out = Dataset(cdf_filename, 'w', format='NETCDF4')
    data_out.file = north_file + ',' + south_file
    data_out.model = 'ADELPHI'
    data_out.filedate = date_str
    
    # establish coordinates (time, lon, lat)
    # time: create dimension and variable
    new_dim = data_out.createDimension('time', data['Time'].size)  # create dimension
    new_var = data_out.createVariable('time', np.float32, 'time')  # variable
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
    var_list = [key for key in data.keys() if key not in ['Time', 'Lon', 'Lat']]
    for variable_name in var_list:
        new_var = data_out.createVariable(variable_name, np.float32,
                                          ('time', 'lon', 'lat'))
        new_var[:] = data[variable_name]  # store data in variable
        #units = [value[-1] for key, value in model_varnames.items() if value[0] == variable_name][0]
        new_var.units = 'V'
    
    # close file
    data_out.close()
    print('Time to create netCDF4 file:', perf_counter()-time0)
    return cdf_filename