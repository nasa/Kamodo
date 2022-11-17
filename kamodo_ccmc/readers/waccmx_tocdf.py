# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 15:05:37 2022
@author: rringuet
Data is already in netCDF4, so just doing a few special things here.
The pressure level variables will need to be inverted, so using Z3 to calculate
    the median km value for each pressure level across the range of times,
    latitudes and longitudes. Also interpolating RHO_CLUBB (the air density)
    from the secondary pressure level grid to the primary one.
NOTE: This assumes that the Z3 and RHO_CLUBB variables are found in the same
    file.
"""
from numpy import array, median, unique, transpose
from glob import glob
from time import perf_counter
from netCDF4 import Dataset
import kamodo_ccmc.readers.reader_utilities as RU
from kamodo import Kamodo

# variable names to keep
var_list = ['Z3', 'RHO_CLUBB', 'datesec']


def convert_all(file_dir):
    '''Find files with RHO_CLUBB and Z3 and prepare h0 files per file.'''

    # find first file and associated files
    print(f'Preparing h0 files for new files found in {file_dir}...', end="")
    t0 = perf_counter()
    files = sorted(glob(file_dir+'*.h?.*.nc'))
    patterns = unique([file[-19:-9] for file in files])  # date strings

    # prepare files
    H_key = False
    for dstr in patterns:
        files = sorted(glob(file_dir+'*'+dstr+'*.nc'))
        h0_test = sum(['h0' in f for f in files])
        if h0_test > 0:
            continue  # skip date string pattern if h0 file already exists
        prepare_h0file(files)  # prepare file
        H_key = True
    if not H_key:
        print('none found.')
    else:
        print(f'Completed in {perf_counter()-t0}s.')
    return True


def prepare_h0file(dstr_files):
    '''Loop through files and perform data wranging. Split into one file per
    timestep.'''

    # set name of output file to have different pattern than original files
    h_type = dstr_files[0][-23:-19]
    h0_file = dstr_files[0].replace(h_type, '.h0.')
    print('\nPreparing', h0_file)
    t_file = perf_counter()

    # set up output file
    data_out = Dataset(h0_file, 'w', format='NETCDF3_64BIT_OFFSET')
    # save coords for ilev interpolation later
    coord_dict = {'ilev': {'data': [], 'units': 'm/m'},
                  'lat': {'data': [], 'units': 'deg'},
                  'lon': {'data': [], 'units': 'deg'}}  # order of keys matters

    # get coord_dict and file dimensions from RHO_CLUBB file
    for file in dstr_files:
        cdf_data = Dataset(file)
        if not any([True for key in cdf_data.variables.keys() if key in
                    ['Z3', 'RHO_CLUBB']]):  # assume datesec is in every file
            cdf_data.close()
            continue
        # loop through dimensions and save
        lev = array(cdf_data.variables['lev'])
        for dim in cdf_data.dimensions.keys():
            if dim in coord_dict.keys():
                coord_dict[dim]['data'] = array(cdf_data.variables[dim])
            if dim in ['time', 'lon', 'lat', 'lev'] and (
                    dim not in data_out.variables.keys()):
                # save dimensions and values to file
                tmp = array(cdf_data.variables[dim])
                if dim == 'time':
                    new_dim = data_out.createDimension(dim, None)
                    len_time = len(tmp)
                else:
                    new_dim = data_out.createDimension(dim, len(tmp))
                new_var = data_out.createVariable(
                    dim, cdf_data.variables[dim].datatype, (dim, ))
                new_var[:] = tmp

        # save datesec variable
        if 'datesec' in cdf_data.variables.keys() and (
                'datesec' not in data_out.variables.keys()):  # copy datesec
            tmp = cdf_data.variables['datesec']
            new_var = data_out.createVariable(
                'datesec', tmp.datatype, tuple(tmp.dimensions))
            new_var[:] = array(tmp)

        # calculate km_ilev grid and store as km_ilev
        if 'Z3' in cdf_data.variables.keys() and (  # (time,) ilev, lat, lon
                'km_ilev' not in data_out.variables.keys()):
            # calculate the median height in km per ilev level
            print('Calculating km grid for pressure level inversion...',
                  end="")
            tmp = cdf_data.variables['Z3']
            data = array(tmp)
            km = median(data, axis=[0, 2, 3])/1000.
            new_dim = data_out.createDimension('km_ilev', len(km))
            new_var = data_out.createVariable('km_ilev', tmp.datatype,
                                              tuple(['km_ilev']))
            new_var[:] = km
            data_out.km_ilev_max = data.max()/1000.
            data_out.km_ilev_min = data.min()/1000.
            print('done.')

        # perform ilev -> lev interpolation for air density
        if 'RHO_CLUBB' in cdf_data.variables.keys() and (  # (t,) ilev,lat,lon
                'RHO_CLUBB_lev' not in data_out.variables.keys()):
            # initialize new variable in output file
            print('Interpolating density to primary pressure level grid...')
            tmp = cdf_data.variables['RHO_CLUBB']
            dim_list = tmp.dimensions
            new_dims = tuple([dim_list[0], 'lev', dim_list[2], dim_list[3]])
            new_var = data_out.createVariable('RHO_CLUBB_lev', tmp.datatype,
                                              new_dims)
            # initialize object for a gridded interpolator
            ko = Kamodo()
            ko.variables, ko._registered = {}, 0
            ko.variables['rho_ilev'] = {'units': 'kg/m**3'}
            # perform interpolation per time slice
            for t in range(len_time):
                ko.variables['rho_ilev']['data'] = array(tmp[t])
                ko = RU.Functionalize_Dataset(ko, coord_dict, 'rho_ilev',
                                              ko.variables['rho_ilev'], True,
                                              'GDZsph')
                # interpolate onto lev grid and save
                new_data = transpose(ko.rho_ilev_ijk(ilev=lev), axes=(1, 0, 2))
                new_var[t] = new_data
                print(f'Completed for time {t+1} out of {len_time}.')
        cdf_data.close()  # close per time step
    data_out.close()
    print(f'{h0_file} created in {perf_counter()-t_file}s.')

    # all done. return.
    return None
