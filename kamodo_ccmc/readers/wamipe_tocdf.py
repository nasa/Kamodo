# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 15:05:37 2022
@author: rringuet
Data is already in netCDF4, so just doing a few special things here.
The pressure level variables will need to be inverted, so using height to calc
    the median km value for each pressure level across the range of times,
    latitudes and longitudes.
"""
from numpy import array, median, unique, transpose
from glob import glob
from os.path import isfile
from time import perf_counter
from netCDF4 import Dataset
import kamodo_ccmc.readers.reader_utilities as RU
from kamodo import Kamodo


def convert_all(file_dir):
    '''Prepare gsm10H files per day.'''

    # find first file and associated files
    print(f'Preparing gsm10H files for new files found in {file_dir}...',
          end="")
    t0 = perf_counter()
    files = sorted(glob(file_dir+'*.gsm10.*.nc'))
    patterns = unique([file[-18:-10] for file in files])  # YYYYMMDD
    
    # prepare files
    H_key = False
    for dstr in patterns:  # loop through days
        if isfile(files[0][:-18].replace('.gsm10.', '.gsm10H.')+dstr+'.nc'):
            continue
        files = sorted(glob(file_dir+'*.gsm10.'+dstr+'*.nc'))  # files for day
        #try:
        prepare_h0file(files)  # prepare file
        H_key = True
        #except:
        #    print('h0 file preparation failed for', files)
        #    return False
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
        if 'RHO_CLUBB' not in cdf_data.variables.keys():
            cdf_data.close()
            continue
        # loop through dimensions and save
        lev = array(cdf_data.variables['lev'])
        for dim in cdf_data.dimensions.keys():
            if dim in coord_dict.keys():
                coord_dict[dim]['data'] = array(cdf_data.variables[dim])
            if dim in ['time', 'lon', 'lat', 'lev']:
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
        if 'datesec' in cdf_data.variables.keys():  # copy datesec to file
            tmp = cdf_data.variables['datesec']
            new_var = data_out.createVariable(
                'datesec', tmp.datatype, tuple(tmp.dimensions))
            new_var[:] = array(tmp)

        # calculate km_ilev grid and store as km_ilev
        if 'Z3' in cdf_data.variables.keys():  # (time,) ilev, lat, lon
            # calculate the median height in km per ilev level
            print('Calculating km grid for pressure level inversion...',
                  end="")
            tmp = cdf_data.variables['Z3']
            km = median(array(tmp), axis=[0, 2, 3])/1000.
            new_dim = data_out.createDimension('km_ilev', len(km))
            new_var = data_out.createVariable('km_ilev', tmp.datatype,
                                              tuple(['km_ilev']))
            new_var[:] = km
            print('done.')

        # perform ilev -> lev interpolation for air density
        if 'RHO_CLUBB' in cdf_data.variables.keys():  # (time,) ilev, lat, lon
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
