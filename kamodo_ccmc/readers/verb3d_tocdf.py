"""
@author: xandrd
"""
from netCDF4 import Dataset
import numpy as np
import os
import pyverbplt
import kamodo_ccmc.readers.reader_utilities as RU
from datetime import datetime

def convert_all(file_dir):
    '''Converting all plt files to netCDF4.'''

    # TODO: add check for files existence

    perp_grid_filename = os.path.join(file_dir, 'Output', 'perp_grid.plt')

    # Work with grid
    grid = pyverbplt.load_plt(perp_grid_filename, squeeze=True)
    data = {'L': grid[0]['arr'],
            'E': grid[1]['arr'],
            'alpha': grid[2]['arr'],
            'pc': grid[3]['arr']}

    cdf_filename = os.path.join(file_dir, 'Output', 'rad_grid.nc')
    var_shape = grid[0]['arr'].shape
    with Dataset(cdf_filename, 'w', format='NETCDF4') as ncfile:
        # Create dimensions
        ncfile.createDimension('L', var_shape[0])
        ncfile.createDimension('E', var_shape[1])
        ncfile.createDimension('alpha', var_shape[2])

        for key, var in data.items():
            data_var = ncfile.createVariable(key, np.float32, ('L', 'E', 'alpha'))
            data_var[:] = var

    # Work with PSD
    psd_filename = os.path.join(file_dir, 'Output', 'OutPSD.dat')
    psd = pyverbplt.load_plt(psd_filename)

    time = np.array([np.float32(t) for t in psd['zone']])
    psd_size = psd['arr'].shape

    nc_files = []
    for t in range(psd_size[0]):
        cdf_filename = os.path.join(file_dir, 'Output', f'OutPSD{t}.nc')
        nc_files.append(cdf_filename)
        with Dataset(cdf_filename, 'w', format='NETCDF4') as ncfile:
            # Create dimensions
            ncfile.createDimension('time', 1)
            ncfile.createDimension('L', psd_size[1])
            ncfile.createDimension('E', psd_size[2])
            ncfile.createDimension('alpha', psd_size[3])

            psd_var = ncfile.createVariable('PSD', np.float32, ('L', 'E', 'alpha'))
            psd_var[:] = psd['arr'][t, :, :, :]

            time_var = ncfile.createVariable('time', np.float32, ('time'))
            time_var[:] = time[t]



    modelname = 'VERB-3D'
    # List of files for Kamodo reader
    list_file = file_dir + modelname + '_list.txt'
    time_file = file_dir + modelname + '_times.txt'

    pattern_files = {'OutPSD': nc_files}
    times_list = list(time)
    times_end = times_list
    if len(times_list) > 1:
        times_end = times_list[1:]
        times_end.append(times_end[-1] + times_list[-1] - times_list[-2])

    # TODO: Change to
    # All times are stored as the number of hours since midnight of the
    # first file of all the files in the given directory.
    times = {'OutPSD': {'start': times_list, 'end': times_end, 'all': times_list}}

    RU.create_timelist(list_file, time_file, modelname,
                       times, pattern_files,
                       datetime(1, 1, 1))

    return True
