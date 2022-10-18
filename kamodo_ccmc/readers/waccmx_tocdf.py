# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 15:05:37 2022
@author: rringuet
Data is already in netCDF4, but is too large to handle transposing and adding
    the neighboring time step to the arrays. Doing here.
"""
from numpy import transpose, array, append, where, insert, median
from glob import glob
from time import perf_counter
from netCDF4 import Dataset
from os.path import isfile

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
    beginning, but the 1979 files have an open time at the end.'''

    # find first file and associated files
    print(f'Converting files found in {file_dir}. ', end="")
    t0 = perf_counter()
    files = sorted(glob(file_dir+'*.h1.*.nc'))
    
    # check for time convention: open at beginning or end?
    cdf_data = Dataset(files[0])
    check_first = cdf_data['datesec'][0]
    check_last = cdf_data['datesec'][-1]
    print(array(cdf_data['datesec']))
    cdf_data.close()
    #cdf_data = Dataset(files[-1])
    #print(array(cdf_data['datesec']))
    #check_last = cdf_data['datesec'].shape[0]
    #cdf_data.close()

    # collect file types/lists of each type
    file_dict = {}
    for n in [1, 2, 3, 4]:
        if isfile(files[0].replace('h1', 'h'+str(n))):
            file_dict['h'+str(n)] = sorted(glob(file_dir+'*.h'+str(n)+'.*.nc'))

    # put in the right cycling order based on time convention
    if check_first == check_last:
        print('help!')
        stop
    elif check_last == 0:  # open time at beginning, go forwards (already good)
        check_first = True  # cycle forwards
    elif check_first == 0:  # open time at the end, go backwards
        for key in file_dict.keys():
            file_dict[key].reverse()
        check_first = False  # cycle backwards
    print('check_first:', check_first, check_last)
    print('File order:', file_dict)

    # convert files while adding a timestep from the previous, one type at a time
    for key in file_dict.keys():
        out_dict, time_out = None, None
        for file in file_dict[key]:
            # convert next file with a time step from the previous file
            #try:
            out_dict, time_out = convert_files(
                    file, check_first, in_dict=out_dict, time_in=time_out)
            #except:
            #    return False
            print('convert_all', time_out, out_dict.keys())

    print(f'Completed in {perf_counter()-t0}s.')
    return True


def convert_files(file, check_first, in_dict=None, time_in=None,
                  verbose=False):
    '''Loop through files and perform data wranging. This includes adding a
    neighboring time step for each variable and transposing into the expected
    coordinate order (time, lon, lat, ilev). This is too avoid performing
    memory intensive calculations each time the reader is executed.
    check_first = True if open time at beg, False if open at end.'''

    print('Converting', file)
    t_file = perf_counter()

    # set name of output file
    for n in [1, 2, 3, 4]:
        if '.h'+str(n)+'.' in file:
            new_file = file.replace('.h'+str(n)+'.', '.h'+str(n)+'v2.')
            h0_file = file.replace('.h'+str(n)+'.', '.h0v2.')

    # too memory intensive to wrangle as a whole
    # try having both in and out files open
    cdf_data = Dataset(file)
    # get variable list set up
    cdf_list = ['datesec'] + [key for key in cdf_data.variables.keys()
                              if key in var_list]
    if 'Z3' in cdf_list:
        cdf_list.remove('Z3')
        cdf_list = ['Z3'] + cdf_list  # do Z3 before others if memory fails

        # perform km averaging and write to h0 file
        # retrieve variable and calculate the median height
        km_time = perf_counter()
        datatype = cdf_data.variables['Z3'].datatype
        height = array(cdf_data.variables['Z3'])
        km = median(array(height), axis=[0, 2, 3])/1000.  # time, ilev, lat, lon
    
        # save in a netCDF4 file with the name h0
        data_out = Dataset(h0_file, 'w', format='NETCDF3_64BIT_OFFSET')
        data_out.model = 'WACCM-X'
        new_dim = data_out.createDimension('km_ilev', len(km))
        new_var = data_out.createVariable('km_ilev', datatype,
                                          tuple(['km_ilev']))
        new_var[:] = km
        data_out.close()
        print(f'Height inversion grid calculated in {perf_counter()-km_time}s.')
    if 'RHO_CLUBB' in cdf_list:
        cdf_list.remove('RHO_CLUBB')
        cdf_list.insert(0, 'RHO_CLUBB')  # do first
    print(cdf_list)

    # send to next file if only one timestep
    if cdf_data.variables['time'].shape[0] == 1:  # just send to next file
        time_out = array(cdf_data.variables['time'])
        print('One time value detected.')

        # loop through variables
        out_dict = {var: array(cdf_data.variables[var]) for var in
                    cdf_list}  # output single time as input time for next
        if 'datesec' in out_dict.keys():
            print('datesec before', out_dict['datesec'])
        #if 'datesec' in out_dict.keys() and check_first:  don't need to change b/c already zero!
            #out_dict['datesec'] -= 86400
        if 'datesec' in out_dict.keys() and not check_first:
            out_dict['datesec'] += 86400
        if 'datesec' in out_dict.keys():
            print('datesec after', out_dict['datesec'])

    # perform data wrangling and save to a new file
    else:
        # format needed to allow for large files
        data_out = Dataset(new_file, 'w', format='NETCDF3_64BIT_OFFSET')

        # loop through dimensions
        for dim in cdf_data.dimensions.keys():
            if dim == 'nbnd' or dim == 'chars':
                continue  # skip these dimensions
            tmp = array(cdf_data.variables[dim])
            # add time at beginning, save last time for output
            if dim == 'time':
                if check_first:
                    time_out = tmp[-1]
                else:
                    time_out = tmp[0]
                print('time_out:', time_out)
                print('time_in:', time_in)
                if time_in != None and check_first:  # add time at beg
                    out = insert(tmp, 0, time_in)
                    print('new time')
                elif time_in != None and not check_first:  # at end
                    out = insert(tmp, len(tmp), time_in)
                    print('new time')
                else:
                    out = tmp
            # prep for longitude wrapping
            elif (dim == 'lon'):
                lon_le180 = list(where(tmp <= 180)[0])
                lon_ge180 = list(where(tmp >= 180)[0])  # repeat 180 for -180 values
                out = append(tmp, 360.) - 180.
            # prep for magnetic longitude wrapping
            elif dim == 'mlon':
                out = append(tmp, 180.)
            else: # loop through other dimensions without changing anything
                out = tmp
            print(dim, len(out), ', ')
            if dim == 'time':  # use unlimited length for time to save memory
                new_dim = data_out.createDimension(dim, None)
                len_time = len(out)
            else:
                new_dim = data_out.createDimension(dim, len(out))
            new_var = data_out.createVariable(
                dim, cdf_data.variables[dim].datatype, (dim, ))
            new_var[:] = out
        time_dim = perf_counter()
        print(f'Dimensions complete in {time_dim-t_file}s.')
            
        # loop through variables
        out_dict = {}
        for var in cdf_list:
            var_start = perf_counter()
            print(f'Converting {var}', end="")
            dim_list = cdf_data.variables[var].dimensions
            print(f'{dim_list}...', end="")
            if len(dim_list) < 4:
                variable = array(cdf_data.variables[var])
                # set output aside BEFORE data wrangling
                if check_first:
                    # deal with time grid (seconds since midnight)
                    if var == 'datesec':  # open at beg
                        print('before', variable[-1], end="")
                        out_dict[var] = variable[-1]  # sec since 12am, already zero!
                        variable[-1] += 86400  # need to be at end of day, not zero
                        print('after', out_dict[var], end="")
                    else:
                        out_dict[var] = variable[-1]  # save last time for output
                else:
                    if var == 'datesec':  # open at end
                        print('before', variable[0], end="")
                        out_dict[var] = variable[0] + 86400  # sec since 12am
                        print('after', out_dict[var], end="")
                    else:
                        out_dict[var] = variable[0]  # save first time for output
    
                # add data from earlier/later file
                if in_dict != None:  
                    if var == 'datesec':
                        print('before', variable, in_dict[var])
                    if check_first:
                        variable = insert(variable, 0, in_dict[var], axis=0)
                        print('earlier time added...', end="")
                    else:
                        variable = insert(variable, variable.shape[0],
                                          in_dict[var], axis=0)
                        print('later time added...', end="")
                    if var == 'datesec':
                        print('after', variable)

                # time, lat/mlat, lon/mlon -> time, lon/mlon, lat/mlat
                if len(dim_list) == 3:
                    # wrap and shift in longitude, different for mlon!
                    if 'mlon' in dim_list:  # only wrap and don't shift
                        tmp = insert(variable, variable.shape[-1],
                                     variable[:, :, 0], axis=3)
                    else:  # use normal longitude shifting and wrapping
                        tmp = variable[:, :, lon_ge180+lon_le180]
                    variable = transpose(tmp, (0, 2, 1))
                    new_dims = tuple([dim_list[0], dim_list[2], dim_list[1]])
                else:  # nothing special needed for time series data, just copy
                    new_dims = tuple(dim_list)

            # loop through time for 4D to save memory
            elif len(dim_list) == 4:
                # set output aside BEFORE data wrangling
                cdf_var = cdf_data.variables[var]  # do NOT convert to array yet
                if check_first:
                    out_dict[var] = array(cdf_var[-1])  # save last time
                else:
                    out_dict[var] = array(cdf_var[0])  # save first time

                # initialize file output variable before time looping
                # (time,) ilev, lat, lon -> (time,) lon, lat, ilev
                new_dims = tuple([dim_list[0], dim_list[3], dim_list[2],
                                 dim_list[1]])
                new_var = data_out.createVariable(
                    var, cdf_data.variables[var].datatype, new_dims)

                # move copied time data to file first, set iterators
                if in_dict != None:
                    # wrap and shift in longitude (all 4D variables)
                    tmp = in_dict[var][:, :, lon_ge180+lon_le180]
                    variable = transpose(tmp, (2, 1, 0))
                    if check_first:
                        new_var[0] = variable
                        t_range = range(1, len_time)
                    else:
                        new_var[len_time-1] = variable
                        t_range = range(0, len_time-1)
                    cdf_range = range(0, len_time-1)  # added time not in file
                else:
                    cdf_range, t_range = range(0, len_time), range(0, len_time)

                # loop through time dimension and copy over to new file
                for i, j in zip(cdf_range, t_range):
                    # wrap and shift in longitude (all 4D variables)
                    tmp = array(cdf_var[i])[:, :, lon_ge180+lon_le180]
                    variable = transpose(tmp, (2, 1, 0))
                    new_var[j] = variable

            # perform ilev -> lev interpolation for air density
            if var == 'RHO_CLUBB':  # convert from ilev grid to lev grid
                # prepare variables for the gridify decorator, one per time
                lev = array(data_out.variables['lev'])
                coord_dict = {'lon': {'data': array(data_out.variables['lon']),
                                      'units': 'deg'},
                              'lat': {'data': array(data_out.variables['lat']),
                                      'units': 'deg'},
                              'ilev': {
                                  'data': array(data_out.variables['ilev']),
                                  'units': 'm/m'}}

                # initialize object for a gridded interpolator
                from kamodo import Kamodo
                ko = Kamodo()
                ko.variables, ko._registered = {}, 0
                ko.variables['rho_ilev'] = {'units': 'kg/m**3'}
                import kamodo_ccmc.readers.reader_utilities as RU

                # initialize new variable in output file
                rho_dims = tuple([dim_list[0], dim_list[3], dim_list[2],
                                  'lev'])
                new_var_lev = data_out.createVariable(
                    var+'_lev', cdf_data.variables[var].datatype, rho_dims)
                
                for i in range(len_time):
                    # copy over array (+ added time) into dict
                    rho_dict = {'data': array(new_var[i]), 'units': 'kg/m**3'}
                    ko = RU.Functionalize_Dataset(ko, coord_dict, 'rho_ilev',
                                                  rho_dict, True, 'GDZsph')
                    # interpolate onto lev grid and save
                    print(f'\nInterpolating for time number {i} out of ' +
                          f'{len_time}')
                    new_var_lev[i] = ko.rho_ilev_ijk(ilev=lev)
                    del ko['rho_ilev_ijk'], ko['rho_ilev']
                del rho_dict, coord_dict  # keep memory usage clean
            var_end = perf_counter()
            print(f'done in {var_end-var_start}s.')
        cdf_data.close()
        data_out.close()
        del cdf_data
        del data_out
        print(f'{file} converted in {perf_counter()-t_file}s.\n')
    
    # all done. return.
    return out_dict, time_out
