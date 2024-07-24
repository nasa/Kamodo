"""
@author: xandrd
"""
from netCDF4 import Dataset
import numpy as np
import os
import pyverbplt
import kamodo_ccmc.readers.reader_utilities as RU
from datetime import datetime
import re
import rbamlib

def get_start_date(file_dir):
    '''Return the start date based on the information from DatabaseInfo1'''

    # Default date_start
    date_start = datetime(1970, 1, 1)

    # Determine if there is a file that contains userinput
    database_filename = os.path.join(file_dir, '..', 'DatabaseInfo1')
    if RU._isfile(database_filename):

        # Define the regex pattern for the date and time
        pattern = r'(\d{4}/\d{2}/\d{2} \d{2}:\d{2})\s+# start_time'

        with open(database_filename, 'r') as file:
            for line in file:
                match = re.search(pattern, line)
                if match:
                    date_str = match.group(1)  # Return only the date part
                    date_start = datetime.strptime(date_str, '%Y/%m/%d %H:%M')

    return date_start


def convert_all(file_dir):
    '''Converting all plt files to netCDF4.'''

    # perp_grid_filename = os.path.join(file_dir, 'Output', 'perp_grid.plt')
    # psd_filename = os.path.join(file_dir, 'Output', 'OutPSD.dat')
    perp_grid_filename = os.path.join(file_dir, 'perp_grid.plt')
    psd_filename = os.path.join(file_dir, 'OutPSD.dat')

    # check for files existence
    if not RU._isfile(perp_grid_filename) or not RU._isfile(psd_filename):
        return False

    # Get the start date
    start_date = get_start_date(file_dir)

    # Work with grid
    grid = pyverbplt.load_plt(perp_grid_filename, squeeze=True)

    data = {'L': grid[0]['arr'],
            'E': grid[1]['arr'],
            'Alpha': grid[2]['arr'],
            'pc': grid[3]['arr']}

    # Kamodo data design rely on concept that all variables that are available by the model should be stored in the files
    # Below we calculate additional variables and store them in the corresponding nc files
    # Create virtual variables
    mu = rbamlib.conv.en2mu(data['E'], data['L'], np.deg2rad(data['Alpha']))
    K = rbamlib.conv.Lal2K(data['L'], np.deg2rad(data['Alpha']))
    data['Mu'] = mu
    data['K'] = K

    # cdf_filename = os.path.join(file_dir, 'Output', 'perp_grid.nc')
    cdf_filename = os.path.join(file_dir, 'perp_grid.nc')
    var_shape = grid[0]['arr'].shape
    grid_file = [cdf_filename]
    with Dataset(cdf_filename, 'w', format='NETCDF4') as ncfile:
        # Create dimensions
        ncfile.createDimension('L', var_shape[0])
        ncfile.createDimension('E', var_shape[1])
        ncfile.createDimension('Alpha', var_shape[2])

        for key, var in data.items():
            data_var = ncfile.createVariable(key, np.float32, ('L', 'E', 'Alpha'))
            data_var[:] = var

    # Work with PSD, PSD on LMK grid and flux
    psd = pyverbplt.load_plt(psd_filename)

    time = np.array([np.float32(t) for t in psd['zone']])
    psd_size = psd['arr'].shape

    nc_psd_files = []
    nc_psdlmk_files = []
    reshaped_pc = np.expand_dims(grid[3]['arr'], axis=0)
    units_constant = 1/3e7 # (c/MeV/cm)^3 after that

    for t in range(psd_size[0]):
        # PSD
        # cdf_filename = os.path.join(file_dir, 'Output', f'OutPSD{t}.nc')
        cdf_filename = os.path.join(file_dir, f'OutPSD_Flux{t}.nc')
        nc_psd_files.append(cdf_filename)
        with Dataset(cdf_filename, 'w', format='NETCDF4') as ncfile:
            # Create dimensions
            ncfile.createDimension('time', 1)
            ncfile.createDimension('L', psd_size[1])
            ncfile.createDimension('E', psd_size[2])
            ncfile.createDimension('Alpha', psd_size[3])

            psd_var = ncfile.createVariable('PSD', np.float32, ('L', 'E', 'Alpha'))
            psd_var[:] = psd['arr'][t, :, :, :] * units_constant

            # Add flux
            flux_var = ncfile.createVariable('Flux', np.float32, ('L', 'E', 'Alpha'))
            flux_var[:] = psd['arr'][t, :, :, :] * reshaped_pc**2

            # Time is number of days from zero - start of the simulation, directly from zone
            time_var = ncfile.createVariable('time', np.float32, ('time'))
            time_var[:] = time[t]

        # PSD_LMK
        cdf_filename = os.path.join(file_dir, f'OutPSD_lmk{t}.nc')
        nc_psdlmk_files.append(cdf_filename)
        with Dataset(cdf_filename, 'w', format='NETCDF4') as ncfile:
            # Create dimensions
            ncfile.createDimension('time', 1)
            ncfile.createDimension('L', psd_size[1])
            ncfile.createDimension('Mu', psd_size[2])
            ncfile.createDimension('K', psd_size[3])

            psd_var = ncfile.createVariable('PSD_2', np.float32, ('L', 'Mu', 'K'))
            psd_var[:] = psd['arr'][t, :, :, :] * units_constant

            # Time is number of days from zero - start of the simulation, directly from zone
            time_var = ncfile.createVariable('time', np.float32, ('time'))
            time_var[:] = time[t]

    modelname = 'VERB-3D'
    # List of files for Kamodo reader
    list_file = file_dir + modelname + '_list.txt'
    time_file = file_dir + modelname + '_times.txt'

    pattern_files = {'OutPSD_Flux': nc_psd_files, 'OutPSD_lmk': nc_psdlmk_files, 'perp_grid': grid_file}
    times_list = list(time * 24)  # Convert to number of hours
    times_end = times_list

    # if len(times_list) > 1:
    #     times_end = times_list[1:]
    #     times_end.append(times_end[-1] + times_list[-1] - times_list[-2])

    # All times are stored as the number of hours since midnight of the
    # first file of all the files in the given directory.
    times = {'OutPSD_Flux': {'start': times_list, 'end': times_end, 'all': times_list},
             'OutPSD_lmk': {'start': times_list, 'end': times_end, 'all': times_list},
             'perp_grid': {'start': [times_list[0]], 'end': [times_end[-1]], 'all': [times_list[0]]}}

    RU.create_timelist(list_file, time_file, modelname,
                       times, pattern_files,
                       start_date)

    return True
