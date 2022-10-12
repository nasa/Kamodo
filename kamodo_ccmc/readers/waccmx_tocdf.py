# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 15:05:37 2022
@author: rringuet
Data is already in netCDF4, but is too large to handle transposing and adding
    the neighboring time step to the arrays. Doing here.
"""
from numpy import transpose, array, append, where, insert
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
    
    # check for time convenion: open at beginning or end?
    cdf_data = Dataset(files[0])
    check_first = cdf_data['datesec'][0]
    check_last = cdf_data['datesec'][-1]
    print(array(cdf_data['datesec']))
    cdf_data.close()
    #cdf_data = Dataset(files[-1])
    #print(array(cdf_data['datesec']))
    #check_last = cdf_data['datesec'].shape[0]
    #cdf_data.close()
    if check_first == check_last:
        print('help!')
        stop
    elif check_last == 0:  # open time at beginning
        good_files = [files[0]]
        files.remove(files[0])
        check_first = True  # cycle forwards
    elif check_first == 0:  # open time at the end
        good_files = [files[-1]]
        files.remove(files[-1])
        files.reverse()
        check_first = False  # cycle backwards
    print('check_first:', check_first, check_last)

    # collect other files for the same day
    for n in [2, 3, 4]:
        if isfile(good_files[0].replace('h1', 'h'+str(n))):
            good_files.append(good_files[0].replace('h1', 'h'+str(n)))
    print('File order:', files)
    print('good_files:', good_files)

    # convert first file with an empty first time step
    #try:
    in_h1dict, in_h2dict, in_h3dict, in_h4dict, time_out = \
            convert_files(good_files, check_first)  # *.h1v2.*.nc etc
    #except:
    #    return False

    # convert remaining file groups while adding a timestep from the previous
    for file in files:
        good_files = [file]
        for n in [2, 3, 4]:
            if isfile(good_files[0].replace('h1', 'h'+str(n))):
                good_files.append(good_files[0].replace('h1', 'h'+str(n)))

        # convert next file with an  first time step
        try:
            in_h1dict, in_h2dict, in_h3dict, in_h4dict, time_out = \
                convert_files(good_files, check_first, in_h1dict=in_h1dict,
                                     in_h2dict=in_h2dict, in_h3dict=in_h3dict,
                                     in_h4dict=in_h4dict, time_in=time_out)
        except:
            return False
        print('convert_all', time_out)

    print(f'Completed in {perf_counter()-t0}s.')
    return True


def convert_files(good_files, check_first, in_h1dict=None, in_h2dict=None,
                  in_h3dict=None, in_h4dict=None, time_in=None,
                  verbose=False):
    '''Loop through files and perform data wranging. This includes adding a
    neighboring time step for each variable and transposing into the expected
    coordinate order (time, lon, lat, ilev). This is too avoid performing
    memory intensive calculations each time the reader is executed.
    check_first = True if open time at beg, False if open at end.'''

    # initialize output variables
    out_h1dict, out_h2dict, out_h3dict, out_h4dict = None, None, None, None

    for file in good_files:
        print('Converting', file)
        t_file = perf_counter()
        # which file is it, and is there info from a previous file?
        if '.h1.' in file and in_h1dict != None:
            t0_dict = in_h1dict
            new_file = file.replace('.h1.', '.h1v2.')
        elif '.h2.' in file and in_h2dict != None:
            t0_dict = in_h2dict
        elif '.h3.' in file and in_h3dict != None:
            t0_dict = in_h3dict
        elif '.h4.' in file and in_h4dict != None:
            t0_dict = in_h4dict
        else:
            t0_dict = None  # if so, then earlier_time should also be None

        # set name of output file
        for n in [1, 2, 3, 4]:
            if '.h'+str(n)+'.' in file:
                new_file = file.replace('.h'+str(n)+'.', '.h'+str(n)+'v2.')

        # too memory intensive to wrangle as a whole
        # try having both in and out files open
        cdf_data = Dataset(file)
        # get variable list set up
        cdf_list = ['datesec'] + [key for key in cdf_data.variables.keys()
                                  if key in var_list]
        if 'Z3' in cdf_list:
            cdf_list.remove('Z3')
            cdf_list = ['Z3'] + cdf_list  # do Z3 first b/c memory
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
                if t0_dict != None:  
                    if var == 'datesec':
                        print('before', variable, t0_dict[var])
                    if check_first:
                        variable = insert(variable, 0, t0_dict[var], axis=0)
                        print('earlier time added...', end="")
                    else:
                        variable = insert(variable, variable.shape[0],
                                          t0_dict[var], axis=0)
                        print('later time added...', end="")
                    if var == 'datesec':
                        print('after', variable)
    
                # perform data wrangling based on dimensions
                # time, lev/ilev, lat, lon -> time, lon, lat, lev/ilev
                if len(dim_list) == 4:
                    # wrap and shift in longitude (all 4D variables)
                    tmp = variable[:, :, :, lon_ge180+lon_le180]
                    variable = transpose(tmp, (0, 3, 2, 1))
                    new_dims = tuple([dim_list[0], dim_list[3], dim_list[2],
                                     dim_list[1]])
                # time, lat/mlat, lon/mlon -> time, lon/mlon, lat/mlat
                elif len(dim_list) == 3:
                    # wrap and shift in longitude, different for mlon!
                    if 'mlon' in dim_list:  # only wrap and don't shift
                        tmp = insert(variable, variable.shape[-1],
                                     variable[:, :, 0], axis=2)
                        #tmp = append(variable, variable[:, :, 0], axis=2)
                    else:  # use normal longitude shifting and wrapping
                        tmp = variable[:, :, lon_ge180+lon_le180]
                    variable = transpose(tmp, (0, 2, 1))
                    new_dims = tuple([dim_list[0], dim_list[2], dim_list[1]])
                else:  # nothing special needed for time series data, just copy
                    new_dims = tuple(dim_list)

                # store in file
                print(variable.shape, new_dims, end="")
                new_var = data_out.createVariable(
                    var, cdf_data.variables[var].datatype, new_dims)
                new_var[:] = variable
                var_end = perf_counter()


                print(f'done in {var_end-var_start}s.')
            cdf_data.close()
            data_out.close()
            del cdf_data
            del data_out

        # assign file-specific values to outputs based on filename
        if '.h1' in file:
            out_h1dict = out_dict
        elif '.h2.' in file:
            out_h2dict = out_dict
        elif '.h3.' in file:
            out_h3dict = out_dict
        elif '.h4.' in file:
            out_h4dict = out_dict
        print(f'{file} converted in {perf_counter()-t_file}s.\n')
    
    # all done. return.
    return out_h1dict, out_h2dict, out_h3dict, out_h4dict, time_out
