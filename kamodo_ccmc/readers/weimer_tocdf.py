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
from kamodo_ccmc.readers.reader_utilities import str_to_hrs


def convert_all(file_dir, txt_files):
    '''Converting all files into one netCDF4.'''

    pattern = txt_files[0][:-18]
    cdf_file = pattern + '.nc'
    to_CDF(cdf_file, txt_files)
    return True


def read_weimerfile(file_name):
    '''Read in data from given file and return in a dictionary of arrays.
    Weimer data starts at 12 MLT and ends at 12 MLT. Longitude is already
    wrapped in the data.'''

    # read in data from file
    file_obj = open(file_name, 'r')
    csv_reader = reader(file_obj, delimiter='\n')  # splits file by \n char
    data = {'LAT': [], 'MLT': [], 'PHI': []}  # R is constant

    # skip header
    for i in range(10):
        line = next(csv_reader)
    # read in data
    for line in csv_reader:
        lat, mlt, r, phi = line[0].split()
        data['LAT'].append(lat)
        data['MLT'].append(mlt)
        data['PHI'].append(phi)
    for key in data.keys():
        data[key] = array(data[key], dtype=float32)
    file_obj.close()
    return data, r


def to_CDF(cdf_filename, files):
    '''Given the filename of the new cdf file and the list of files, put all
    the data into a single netcdf4 file.'''

    print('Converted data file not found. Converting files with ' +
          f'{cdf_filename.split(".")[0]} naming pattern to a netCDF4 file.')  # , end="")
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
    new_var = data_out.createVariable('PHI', float32,
                                      ('time', 'lon', 'lat'))
    new_var[0] = reshape(data['PHI'], (lon.shape[0], lat.shape[0]))
    new_var.units = 'kV'

    # loop through remaining files
    for i, f in enumerate(files[1:]):
        data, r = read_weimerfile(f)
        new_var[i+1] = reshape(data['PHI'], (lon.shape[0], lat.shape[0]))

    # close file
    data_out.close()
    print('Time to create netCDF4 file:', perf_counter()-time0)
    return None
