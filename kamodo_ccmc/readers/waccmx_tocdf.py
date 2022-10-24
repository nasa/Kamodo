# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 15:05:37 2022
@author: rringuet
Data is already in netCDF4, but is too large to handle transposing and adding
    the neighboring time step to the arrays. Doing here.
"""
from numpy import array, append, where, insert, median, transpose
from glob import glob
from time import perf_counter
from netCDF4 import Dataset
from os.path import isfile
import kamodo_ccmc.readers.reader_utilities as RU
from kamodo import Kamodo

# variable names to keep
var_list = ['co2vmr', 'ch4vmr', 'n2ovmr', 'f11vmr', 'f12vmr', 'sol_tsi',
            'colat_crit1', 'colat_crit2', 'ED1', 'ED2', 'EDYN_ZIGM11_PED',
            'EDYN_ZIGM2_HAL', 'ElecColDens', 'H', 'O', 'O2', 'OMEGA', 'PHIM2D',
            'PS', 'T', 'TElec', 'TIon', 'U', 'UI', 'V', 'VI', 'WI', 'Z3', 'e',
            'RHO_CLUBB', 'Z3GM', 'EDens', 'HMF2', 'NMF2', 'CO2', 'N', 'NO', 
            'OpDens', 'QCO2', 'QHC2S', 'QJOULE', 'QNO', 'QO3', 'QO3P',
            'QRS_TOT', 'SolIonRate_Tot', 'TTGW', 'UTGW_TOTAL', 'OMEGA_08_COS',
            'OMEGA_08_SIN', 'OMEGA_12_COS', 'OMEGA_12_SIN', 'OMEGA_24_COS',
            'OMEGA_24_SIN', 'T_08_COS', 'T_08_SIN', 'T_12_COS', 'T_12_SIN',
            'T_24_COS', 'T_24_SIN', 'U_08_COS', 'U_08_SIN', 'U_12_COS',
            'U_12_SIN', 'U_24_COS', 'U_24_SIN', 'V_08_COS', 'V_08_SIN',
            'V_12_COS', 'V_12_SIN', 'V_24_COS', 'V_24_SIN', 'OPLUS']

def convert_all(file_dir):
    '''Get all data from last timestep of previous file. Finds first file
    and converts them all in order. The CCMC files have an open time at the
    beginning, but the 1979 files have an open time at the end. Converts into
    one file per time step with all coords in first file only.'''

    # find first file and associated files
    print(f'Converting new files found in {file_dir}. ', end="")
    t0 = perf_counter()
    files = sorted(glob(file_dir+'*.h1.*.nc'))
    
    # collect file types/lists of each type
    file_dict = {}
    for n in [1, 2, 3, 4]:
        if isfile(files[0].replace('h1', 'h'+str(n))):
            file_dict['h'+str(n)] = sorted(glob(file_dir+'*.h'+str(n)+'.*.nc'))
            file_dict['h'+str(n)+'v2'] = []

    # convert files while adding a timestep from the previous, one type at a time
    for key in file_dict.keys():
        first_file = True  # write coords to file for first file only
        for file in file_dict[key]:
            # skip if already converted
            file_pattern = file[:-8].replace(key, key+'v2')
            files = sorted(glob(file_pattern+'*.nc'))
            if len(files) > 0:
                file_dict[key+'v2'].extend(files)
                first_file = False
                continue
            #try:
            convert_files(file, first_file)
            first_file = False  # save space and don't write coords to file
            new_files = sorted(glob(file_pattern+'*.nc'))
            file_dict[key+'v2'].extend(new_files)
            #except:
            #    print('File conversion failed for', file)
            #    return False, file_dict
        if first_file:  # no new files found of pattern 'key'
            print('No new files found of pattern', key)
    for n in [1, 2, 3, 4]:
        del file_dict['h'+str(n)]  # return only lists of hnv2 files
    print(f'Completed in {perf_counter()-t0}s.')

    return True, file_dict


def convert_files(file, first_file):
    '''Loop through files and perform data wranging. Split into one file per
    timestep. The coordinates are only written to the first file to save space.
    '''

    print('Converting', file)
    t_file = perf_counter()

    # set name of output file
    for n in [1, 2, 3, 4]:
        if '.h'+str(n)+'.' in file:
            new_file = file.replace('.h'+str(n)+'.', '.h'+str(n)+'v2.')

    # too memory intensive to wrangle as a whole
    # try having both in and out files open
    cdf_data = Dataset(file)
    # get variable list set up
    cdf_list = [key for key in cdf_data.variables.keys() if key in var_list]
    print(cdf_list)

    # perform data wrangling and save to one file per timestep
    # save coords for ilev interpolation later
    time = array(cdf_data.variables['time'])  # hrs since Jan 1, 1979 at 12am.
    datesec = array(cdf_data.variables['datesec'])  # seconds since midnight
    coord_dict = {'lon': {'data': [], 'units': 'deg'},
                  'lat': {'data': [], 'units': 'deg'},
                  'ilev': {'data': [], 'units': 'm/m'}}
    lev = array(cdf_data.variables['lev'])
    for t in range(len(time)):
        # format needed to allow for large files
        t_counter = perf_counter()
        time_filename = new_file[:-8] + str(f'{datesec[t]:05}') + '.nc'
        data_out = Dataset(time_filename, 'w', format='NETCDF3_64BIT_OFFSET')
        
        # loop through dimensions
        for dim in cdf_data.dimensions.keys():
            if dim == 'nbnd' or dim == 'chars':
                continue  # skip these dimensions
            tmp = array(cdf_data.variables[dim])
            # prep for longitude wrapping
            if (dim == 'lon'):
                lon_le180 = list(where(tmp <= 180)[0])
                lon_ge180 = list(where(tmp >= 180)[0])  # repeat 180 for -180 values
                out = append(tmp, 360.) - 180.
            # prep for magnetic longitude wrapping
            elif dim == 'mlon':
                out = append(tmp, 180.)
            else: # loop through other dimensions without changing anything
                out = tmp
            if dim in coord_dict.keys():
                coord_dict[dim]['data'] = out
            if t == 0:
                print(dim, len(out), ', ')
            if dim == 'time':  # use one time per file
                new_dim = data_out.createDimension(dim, 1)
                new_var = data_out.createVariable(
                    dim, cdf_data.variables[dim].datatype, (dim, ))
                new_var[:] = time[t]
            else:
                new_dim = data_out.createDimension(dim, len(out))
                if first_file:  # only write coords for first file.
                    new_var = data_out.createVariable(
                        dim, cdf_data.variables[dim].datatype, (dim, ))
                    new_var[:] = out
    
        # loop through variables
        for var in cdf_list:
            var_start = perf_counter()
            dim_list = cdf_data.variables[var].dimensions
            datatype = cdf_data.variables[var].datatype
            # (ilev,) lat/mlat, lon/mlon -> lon/mlon, lat/mlat, (ilev)
            tmp = array(cdf_data.variables[var][t]).T  # also fine for 1D
            if len(dim_list) < 4:
                if len(dim_list) == 3:
                    # wrap and shift in longitude, different for mlon!
                    if 'mlon' in dim_list:  # only wrap and don't shift
                        variable = insert(tmp, tmp.shape[0],
                                     tmp[0], axis=0)
                    else:  # use normal longitude shifting and wrapping
                        variable = tmp[lon_ge180+lon_le180]
                    new_dims = tuple([dim_list[0], dim_list[2], dim_list[1]])
                else:  # nothing special needed for time series data, just copy
                    new_dims = tuple(dim_list)
                    variable = tmp
            elif len(dim_list) == 4:
                # (time,) ilev, lat, lon -> (time,) lon, lat, ilev
                new_dims = tuple([dim_list[0], dim_list[3], dim_list[2],
                                 dim_list[1]])
                # wrap and shift in longitude (all 4D variables)
                variable = tmp[lon_ge180+lon_le180]
            
            # write data for time step to file
            new_var = data_out.createVariable(var, datatype, new_dims)
            new_var[:] = variable

            # calculate km_ilev grid and store as km_ilev
            if var == 'Z3':
                # retrieve variable and calculate the median height
                km = median(variable, axis=[0, 1])/1000.  # lon, lat, ilev
                new_dim = data_out.createDimension('km_ilev', len(km))
                new_var = data_out.createVariable('km_ilev', datatype,
                                                  tuple(['km_ilev']))
                new_var[:] = km

            # perform ilev -> lev interpolation for air density
            if var == 'RHO_CLUBB':  # convert from ilev grid to lev grid
                # initialize new variable in output file
                new_dims = tuple([dim_list[0], dim_list[3], dim_list[2],
                                  'lev'])
                new_var = data_out.createVariable(var+'_lev', datatype,
                                                  new_dims)

                # initialize object for a gridded interpolator and create
                ko = Kamodo()
                ko.variables, ko._registered = {}, 0
                ko.variables['rho_ilev'] = {'units': 'kg/m**3'}
                ko = RU.Functionalize_Dataset(ko, coord_dict, 'rho_ilev',
                                              {'data': variable, 'units':
                                               'kg/m**3'}, True, 'GDZsph')
                # interpolate onto lev grid and save
                new_data = transpose(ko.rho_ilev_ijk(ilev=lev), axes=(1, 0, 2))
                new_var[:] = new_data
                del ko  # clean memory before continuing
        data_out.close()  # close per time step
        print(f'Time {t} out of {len(time)} times done in ' +
              f'{perf_counter()-t_counter}s.')
    cdf_data.close()
    print(f'{file} converted in {perf_counter()-t_file}s.\n')
    
    # all done. return.
    return None
