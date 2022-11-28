# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 15:05:37 2022
@author: rringuet
Data is already in netCDF4, so just doing a few special things here.
The pressure level variables will need to be inverted, so using height in both
    the neutral and density files to calculate
    the median km value for each pressure level across the range of times,
    latitudes and longitudes.
Adapted from the waccmx post processing script.
"""
from numpy import array, median, float32
from time import perf_counter
from netCDF4 import Dataset


def convert_all(file_dir, pattern_files, times):
    '''Find files with height and prepare km post-processing file.'''

    t0 = perf_counter()
    pattern = list(pattern_files.keys())[0]  # only one pattern for TIEGCM data
    out_file = file_dir + 'TIEGCM_km.nc'
    data_out = Dataset(out_file, 'w')
    new_dim = data_out.createDimension('time', len(times[pattern]['start']))
    new_var = data_out.createVariable('time', float32, ('time'))
    new_var[:] = times[pattern]['start']
    # loop through neutral/density files and add median/max/min vals to outfile
    for i in range(len(pattern_files[pattern])):
        data = Dataset(pattern_files[pattern][i])
        ilev = data.variables['lev']  # grid is constant in time
        if 'ZGMID' in data.variables.keys():  # H_ilev logic
            h_ilev = data.variables['ZGMID']
            if i == 0:
                new_dim = data_out.createDimension('ilev', ilev.shape[0])
                new_var_t = data_out.createVariable(
                    'km_ilev_t', h_ilev.datatype, tuple(['time', 'ilev']))
                new_var_max = data_out.createVariable(
                    'km_ilev_max', h_ilev.datatype, ('time'))
                new_var_min = data_out.createVariable(
                    'km_ilev_min', h_ilev.datatype, ('time'))
            height = array(h_ilev)  # in cm
            new_var_t[i] = median(height, axis=[0, 2, 3])/100000.  # ilev = 1
            new_var_max[i] = height.max()/100000.
            new_var_min[i] = height.min()/100000.
        ilev1 = data.variables['ilev']  # grid is constant in time
        if 'ZG' in data.variables.keys():  # H_ilev1 logic
            h_ilev1 = data.variables['ZG']
            if i == 0:
                new_dim = data_out.createDimension('ilev1', ilev1.shape[0])
                new_var1_t = data_out.createVariable(
                    'km_ilev1_t', h_ilev1.datatype, tuple(['time', 'ilev1']))
                new_var1_max = data_out.createVariable(
                    'km_ilev1_max', h_ilev1.datatype, ('time'))
                new_var1_min = data_out.createVariable(
                    'km_ilev1_min', h_ilev1.datatype, ('time'))
            height = array(h_ilev1)  # in cm
            new_var1_t[i] = median(height, axis=[0, 2, 3])/100000.  # plev = 1
            new_var1_max[i] = height.max()/100000.
            new_var1_min[i] = height.min()/100000.
        data.close()

    # calculate median/max/min over time
    if 'km_ilev_t' in data_out.variables.keys():
        new_var = data_out.createVariable('km_ilev', new_var_t.datatype,
                                          ('ilev'))
        new_var[:] = median(array(new_var_t), axis=0)
        data_out.km_ilev_max = array(new_var_max).max()
        data_out.km_ilev_min = array(new_var_min).min()
    if 'km_ilev1_t' in data_out.variables.keys():
        new_var = data_out.createVariable('km_ilev1', new_var1_t.datatype,
                                          ('ilev1'))
        new_var[:] = median(array(new_var1_t), axis=0)
        data_out.km_ilev1_max = array(new_var1_max).max()
        data_out.km_ilev1_min = array(new_var1_min).min()

    # all done. close output file and return.
    data_out.close()
    print(f'{out_file} created in {perf_counter()-t0}s.')
    return None
