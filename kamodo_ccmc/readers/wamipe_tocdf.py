# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 15:05:37 2022
@author: rringuet
Data is already in netCDF4, so just doing a few special things here.
The pressure level variables will need to be inverted, so using height to calc
    the median km value for each pressure level across the range of times,
    latitudes and longitudes.
"""
from numpy import array, median, unique
from glob import glob
from os.path import isfile, basename
from time import perf_counter
from netCDF4 import Dataset


def convert_all(file_dir):
    '''Prepare one gsm10H file per run.'''

    # find first file and associated files
    t0 = perf_counter()
    files = sorted(glob(file_dir+'*.nc'))
    patterns = unique([basename(f)[:-19] for f in files])
    # figure out which pattern has the 'height' variable, if any
    height_pattern = ''
    for p in patterns:
        files = sorted(glob(file_dir+p+'*.nc'))
        cdf_data = Dataset(files[0])
        if 'height' not in cdf_data.variables.keys():
            cdf_data.close()
            continue
        height_pattern = p
    if height_pattern == '':
        print('No files found with a height(ilev) to invert.')
        return None

    # prepare files
    files = sorted(glob(file_dir+height_pattern+'*.nc'))
    if isfile(files[0][:-18] + 'h0.nc'):
        return None  # calculations already done
    prepare_h0file(files)  # prepare file
    print(f'Completed in {perf_counter()-t0}s.')
    return None


def prepare_h0file(dstr_files):
    '''Loop through files and perform data wranging. Split into one file per
    timestep.'''

    # set name of output file to have different pattern than original files
    h0_file = dstr_files[0][:-18] + 'h0.nc'
    print('\nPreparing', h0_file)
    t_file = perf_counter()

    # make output file
    data_out = Dataset(h0_file, 'w', format='NETCDF3_64BIT_OFFSET')
    for i, file in enumerate(dstr_files):
        cdf_data = Dataset(file)
        # no time grid in files, and assume no pressure level
        # figure out length of vertical grid
        tmp = cdf_data.variables['height']
        if i == 0:  # initialize output dimensions and variables
            new_dim = data_out.createDimension('time', None)  # unlimited
            new_dim = data_out.createDimension('ilev', tmp.shape[0])
            new_var_t = data_out.createVariable(
                'km_ilev_t', tmp.datatype, tuple(['time', 'ilev']))
            new_var_max = data_out.createVariable(
                'km_ilev_max', tmp.datatype, tuple(['time', 'ilev']))
            new_var_min = data_out.createVariable(
                'km_ilev_min', tmp.datatype, tuple(['time', 'ilev']))
        data = array(tmp)
        new_var_t[i] = median(data, axis=[1, 2])/1000.
        new_var_max[i] = data.max()/1000.
        new_var_min[i] = data.min()/1000.
        cdf_data.close()  # close per time step
    # calculate median over time and save to file
    km_arr = array(data_out.variables['km_ilev_t'])
    new_var = data_out.createVariable(
        'km_ilev', data_out.variables['km_ilev_t'].datatype,
        tuple(['ilev']))
    new_var[:] = median(km_arr, axis=[0])  # median over time
    data_out.km_max = array(data_out.variables['km_ilev_max']).max()
    data_out.km_min = array(data_out.variables['km_ilev_min']).min()
    # close and return
    data_out.close()
    print(f'{h0_file} created in {perf_counter()-t_file}s.')

    # all done. return.
    return None
