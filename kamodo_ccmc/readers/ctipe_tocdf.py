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
from numpy import array, median
from time import perf_counter
from netCDF4 import Dataset


def convert_all(file_dir, pattern_files, times):
    '''Find files with height and prepare km post-processing file.'''

    t0 = perf_counter()
    patterns = list(pattern_files.keys())  # e.g. 'height', 'density', 'neutral
    out_file = file_dir + 'CTIPe_km.nc'
    data_out = Dataset(out_file, 'w')
    # loop through neutral/density files and add median/max/min vals to outfile
    for i in range(len(pattern_files[patterns[0]])):
        if 'density' in patterns:  # H_ilev1 logic
            d_data = Dataset(pattern_files['density'][i])
            tmp = d_data.variables['plev']  # plev grid is constant in time
            d_ht = d_data.variables['height']
            if i == 0:
                new_dim = data_out.createDimension(
                    'time1', len(times['density']['start']))
                new_var = data_out.createVariable(
                    'time1', d_ht.datatype, ('time1'))
                new_var[:] = times['density']['start']
                new_dim = data_out.createDimension('ilev1', tmp.shape[0])
                new_var1_t = data_out.createVariable(
                    'km_ilev1_t', d_ht.datatype, tuple(['time1', 'ilev1']))
                new_var1_max = data_out.createVariable(
                    'km_ilev1_max', d_ht.datatype, tuple(['time1']))
                new_var1_min = data_out.createVariable(
                    'km_ilev1_min', d_ht.datatype, tuple(['time1']))
            d_height = array(d_ht)
            new_var1_t[i] = median(d_height, axis=[0, 2, 3])/1000.  # plev = 1
            new_var1_max[i] = d_height.max()/1000.
            new_var1_min[i] = d_height.min()/1000.
            d_data.close()
        if 'neutral' in patterns:  # H_ilev logic
            n_data = Dataset(pattern_files['neutral'][i])
            tmp = n_data.variables['plev']  # plev grid is constant in time
            n_ht = n_data.variables['height']
            if i == 0:
                new_dim = data_out.createDimension(
                    'time', len(times['neutral']['start']))
                new_var = data_out.createVariable(
                    'time', n_ht.datatype, ('time'))
                new_var[:] = times['neutral']['start']
                new_dim = data_out.createDimension('ilev', tmp.shape[0])
                new_var_t = data_out.createVariable(
                    'km_ilev_t', n_ht.datatype, tuple(['time', 'ilev']))
                new_var_max = data_out.createVariable(
                    'km_ilev_max', n_ht.datatype, tuple(['time']))
                new_var_min = data_out.createVariable(
                    'km_ilev_min', n_ht.datatype, tuple(['time']))
            n_height = array(n_ht)
            new_var_t[i] = median(n_height, axis=[0, 2, 3])/1000.  # plev = 1
            new_var_max[i] = n_height.max()/1000.
            new_var_min[i] = n_height.min()/1000.
            n_data.close()
    # calculate median/max/min over time
    if 'density' in patterns:
        new_var = data_out.createVariable('km_ilev1', new_var1_t.datatype,
                                          tuple(['ilev1']))
        new_var[:] = median(array(new_var1_t), axis=0)
        data_out.km_ilev1_max = array(new_var1_max).max()
        data_out.km_ilev1_min = array(new_var1_min).min()
    if 'neutral' in patterns:
        new_var = data_out.createVariable('km_ilev', new_var_t.datatype,
                                          tuple(['ilev']))
        new_var[:] = median(array(new_var_t), axis=0)
        data_out.km_ilev_max = array(new_var_max).max()
        data_out.km_ilev_min = array(new_var_min).min()
    # all done. close output file and return.
    data_out.close()
    print(f'{out_file} created in {perf_counter()-t0}s.')
    return None
