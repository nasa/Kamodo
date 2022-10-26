# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 15:05:37 2022
@author: rringuet
Data is already in netCDF4, but is too large to handle transposing and adding
    the neighboring time step to the arrays. Doing here.
Converting chunk of data into a list of time slices directly from the file
    takes about 30-40 seconds per file PER VARIABLE! This is too long to
    execute each time the model reader is called. Better in the long run to
    convert the files and leave the data in the cdf file (MUCH faster).
"""
from numpy import array, append, where, insert, median, transpose, unique
from glob import glob
from os.path import isfile
from time import perf_counter
from netCDF4 import Dataset
from datetime import datetime, timedelta
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
            'V_12_COS', 'V_12_SIN', 'V_24_COS', 'V_24_SIN', 'OPLUS', 'datesec']

def convert_all(file_dir):
    '''Get all data from last timestep of previous file. Converts into
    one file per time step because files are too large.'''

    # find first file and associated files
    print(f'Converting new files found in {file_dir}. ', end="")
    t0 = perf_counter()
    files = sorted(glob(file_dir+'*.h?.*.nc'))
    patterns = unique([file[-22:-20] for file in files])  # h1, h2, etc
    file_dict = {}
    for p in patterns:
        file_dict[p] = sorted(glob(file_dir+'*.'+p+'.*.nc'))
    
    # convert files while adding a timestep from the previous, one type at a time
    key_list = list(file_dict.keys())
    H_key = False
    for key in key_list:
        first_file = True
        file_dict[key+'v2'] = []
        for file in file_dict[key]:
            # skip if already converted
            file_pattern = file[:-8].replace('.'+key+'.', '.'+key+'v2.')
            files = sorted(glob(file_pattern+'*.nc'))
            if len(files) > 1:  # one file could be from previous file's data
                file_dict[key+'v2'].extend(files)
                continue
            try:
                convert_files(file, key)
                first_file = False
                new_files = sorted(glob(file_pattern+'*.nc'))
                file_dict[key+'v2'].extend(new_files)
            except:
                print('File conversion failed for', file)
                return False, file_dict
        if first_file:  # no new files found of pattern 'key'
            print('No new files found of pattern', key)
        del file_dict[key]

        # find which file pattern has Z3 (need to collect all into one file)
        cdf_data = Dataset(file_dict[key+'v2'][0])
        if 'Z3' in cdf_data.variables.keys() and not H_key:
            H_key = key
        cdf_data.close()

    # checking for H_ilev chunked files (h0v2 files)
    print('Checking for pressure level files...', end="")
    Z3chunk_files = sorted(glob(file_dir+'*.'+H_key+'.*.nc'))  # original files
    Z3_files = sorted(glob(file_dir+'*.'+H_key+'v2.*.nc'))  # all slice files
    H_files = [f.replace(H_key, 'h0v2') for f in Z3chunk_files]  # new files
    if sum([isfile(f) for f in H_files]) < len(H_files):
        print('some found missing. Preparing files.')
        chunk_files = []
        for f in H_files:
            # collect time slices from files and append to h0v2 chunk files
            data_out = Dataset(f, 'w', format='NETCDF3_64BIT_OFFSET')
            slice_files = sorted(glob(f[:-8].replace('h0v2', H_key+'v2')+'*.nc'))
            chunk_files.append(slice_files)  # catch missing one at end
            # get grids from first file, except time
            cdf_data = Dataset(slice_files[0])
            for dim in cdf_data.variables['Z3'].dimensions:
                out = array(cdf_data.variables[dim])
                if dim == 'time':
                    new_dim = data_out.createDimension(dim, None)
                    time_var = data_out.createVariable(
                        dim, cdf_data.variables[dim].datatype, (dim, ))
                else:
                    new_dim = data_out.createDimension(dim, len(out))
                    new_var = data_out.createVariable(
                        dim, cdf_data.variables[dim].datatype, (dim, ))
                    new_var[:] = out
            Z3_var = data_out.createVariable(
                'Z3', cdf_data.variables['Z3'].datatype,
                cdf_data.variables['Z3'].dimensions)
            datesec_var = data_out.createVariable(
                'datesec', cdf_data.variables['datesec'].datatype,
                cdf_data.variables['datesec'].dimensions)
            cdf_data.close()
            
            # loop through files for times and Z3 slices
            for s in slice_files:
                cdf_data = Dataset(s)
                time_var[s] = array(cdf_data.variables['time'])[0]
                datesec_var[s] = array(cdf_data.variables['datesec'])[0]
                Z3_var[s] = array(cdf_data.variables['Z3'])[0]
                cdf_data.close()
            
            # catch extra file at end
            if f == H_files[-1]:
                missing_file = [f for f in Z3_files if f not in chunk_files]
                for file in missing_file:
                    # append remaining time slices to data_out
                    cdf_data = Dataset(file)
                    time_var[s] = array(cdf_data.variables['time'])[0]
                    datesec_var[s] = array(cdf_data.variables['datesec'])[0]
                    Z3_var[s] = array(cdf_data.variables['Z3'])[0]
                    cdf_data.close()
            data_out.close()
            print(f, 'prepared.')
    else:
        print('found.')
    print(f'Completed in {perf_counter()-t0}s.')

    return True, file_dict


def convert_files(file, h_type):
    '''Loop through files and perform data wranging. Split into one file per
    timestep.'''

    print('Converting', file)
    t_file = perf_counter()

    # set name of output file to have different pattern that original files
    new_file = file.replace('.'+h_type+'.', '.'+h_type+'v2.')

    # too memory intensive to wrangle as a whole
    # try having both in and out files open
    cdf_data = Dataset(file)
    # get variable list set up
    cdf_list = [key for key in cdf_data.variables.keys() if key in var_list]

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

        # time at end of files is 0 to indicate midnight on NEXT day
        if int(datesec[t]) == 0 and t == len(time)-1:
            corrected_date = datetime.strftime(datetime.strptime(
                new_file[-19:-9], '%Y-%m-%d') + timedelta(days=1), '%Y-%m-%d')
            time_filename = new_file[:-19] + corrected_date +\
                str(f'-{datesec[t]:05}') + '.nc'
        else:
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
            if dim == 'time':  # use one time per file
                out = [time[t]]
            # save dimensions and values to file
            new_dim = data_out.createDimension(dim, len(out))
            new_var = data_out.createVariable(
                dim, cdf_data.variables[dim].datatype, (dim, ))
            new_var[:] = out
            del new_var
    
        # loop through variables
        for var in cdf_list:
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
            del new_var

            # calculate km_ilev grid and store as km_ilev
            if var == 'Z3':
                # retrieve variable and calculate the median height
                km = median(variable, axis=[0, 1])/1000.  # lon, lat, ilev
                new_dim = data_out.createDimension('km_ilev', len(km))
                new_var = data_out.createVariable('km_ilev', datatype,
                                                  tuple(['km_ilev']))
                new_var[:] = km
                del new_var

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
                del ko, new_var  # clean memory before continuing
        data_out.close()  # close per time step
        print(f'Time {t+1} out of {len(time)} times done in ' +
              f'{perf_counter()-t_counter}s.')
    cdf_data.close()
    print(f'{file} converted in {perf_counter()-t_file}s.\n')
    
    # all done. return.
    return None
