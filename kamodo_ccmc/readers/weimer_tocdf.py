# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 11:57:13 2022

@author: rringuet
"""
from csv import reader
from time import perf_counter
from datetime import datetime, timezone
from netCDF4 import Dataset
from numpy import array, float32, unique, append, reshape
from os.path import getsize
from statistics import mode
from kamodo_ccmc.readers.reader_utilities import str_to_hrs


def convert_all(files):
    '''Converting all files into one netCDF4. Skips files of different sizes
    to avoid file errors. Typically, this is the first file and any list files.
    '''

    sizes = [getsize(f)/1024. for f in files]  # file sizes in KB
    good_files = [f for i, f in enumerate(files) if sizes[i] == mode(sizes)]
    bad_files = [f for f in files if f not in good_files]
    if len(bad_files) > 0:
        print('Skipping the following files: ', bad_files)
    cdf_file = good_files[0][:-18] + '.nc'
    to_CDF(cdf_file, good_files)  # only perform analysis on good data
    return True


def read_weimerfile(file_name):
    '''Read in data from given file and return in a dictionary of arrays.
    Weimer data starts at 12 MLT and ends at 12 MLT. Longitude is already
    wrapped in the data.'''

    # read in data from file
    file_obj = open(file_name, 'r')
    csv_reader = reader(file_obj, delimiter='\n')  # splits file by \n char

    # skip header
    for i in range(9):
        line = next(csv_reader)
    # initialize based on what variables are present
    data = {key: [] for key in line[0][1:].split()}
    var_list = list(data.keys())
    line = next(csv_reader)  # skip units line
    # read in data
    for line in csv_reader:
        for i, val in enumerate(line[0].split()):
            data[var_list[i]].append(val)
    for key in data.keys():
        data[key] = array(data[key], dtype=float32)
    file_obj.close()
    r = data['R'][0]  # constant value
    del data['R']
    return data, r


def to_CDF(cdf_filename, files):
    '''Given the filename of the new cdf file and the list of files, put all
    the data into a single netcdf4 file.'''

    print('Converted data file not found. Converting files with ' +
          f'{cdf_filename.split(".")[0]} naming pattern to a netCDF4 file.')
    time0 = perf_counter()
    # choosing netCDF format that works for large files
    data_out = Dataset(cdf_filename, 'w', format='NETCDF3_64BIT_OFFSET')
    data_out.file = ''.join([f+',' for f in files])[:-1]
    data_out.model = 'Weimer'

    # figure out date and time information from filenames
    filedate = datetime.strptime(
        files[0][-17:-9]+' 00:00:00', '%Y%m%d %H:%M:%S').replace(
            tzinfo=timezone.utc)
    data_out.filedate = filedate.strftime('%Y-%m-%d %H:%M:%S')
    date_strs = [f[-17:-4] for f in files]
    time = str_to_hrs(date_strs, filedate, '%Y%m%d_%H%M')

    # establish time coordinate (time, lon, lat)
    # time: create dimension and variable
    new_dim = data_out.createDimension('time', len(files))
    new_var = data_out.createVariable('time', float32, 'time')
    new_var[:] = time
    new_var.units = 'hr'

    # get lon and lat information from first file
    data, r = read_weimerfile(files[0])
    var_list = [var for var in data.keys() if var not in ['LAT', 'MLT']]
    data_out.radius = float32(r)
    tmp = unique(data['MLT']) * 15.  # files starts at -180 and ends at 180
    lon = append(tmp - 180., 180.)  # lon already wrapped in file
    lat = unique(data['LAT'])

    # lon: create dimension and variable
    new_dim = data_out.createDimension('lon', lon.size)
    new_var = data_out.createVariable('lon', float32, 'lon')
    new_var[:] = lon  # store data for dimension in variable
    new_var.units = 'deg'
    # lat: create dimension and variable
    new_dim = data_out.createDimension('lat', lat.size)
    new_var = data_out.createVariable('lat', float32, 'lat')
    new_var[:] = lat  # store data for dimension in variable
    new_var.units = 'deg'

    # copy over variable to file (only first time in data['PHI'])
    for var in var_list:
        new_var = data_out.createVariable(var, float32,
                                          ('time', 'lon', 'lat'))
        new_var[0] = reshape(data[var], (lon.shape[0], lat.shape[0]))

    # loop through remaining files
    for i, f in enumerate(files[1:]):
        data, r = read_weimerfile(f)
        for var in var_list:
            data_out[var][i+1] = reshape(data[var],
                                         (lon.shape[0], lat.shape[0]))

    # close file
    data_out.close()
    print('Time to create netCDF4 file:', perf_counter()-time0)
    return None
