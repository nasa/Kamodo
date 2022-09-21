'''
Written by Rebecca Ringuette, 2022
'''

from datetime import datetime, timezone
from numpy import vectorize

# 'C:/Users/rringuet/Kamodo_WinDev1/AmGEO/'
model_varnames = {'E_ph': ['E_east', 'Electric Field (eastward)',
                           0, 'SM', 'sph', ['time', 'lon', 'lat'], 'V/m'],
                  'E_th': ['E_north', 'Electric Field (equatorward)',
                           1, 'SM', 'sph', ['time', 'lon', 'lat'], 'V/m'],
                  'cond_hall': ['Sigma_H', 'Ovation Pyme Hall Conductance',
                                2, 'SM', 'sph', ['time', 'lon', 'lat'], 'S'],
                  # units=mho=S
                  'cond_ped': ['Sigma_P', 'Ovation Pyme Pedersen Conductance',
                               3, 'SM', 'sph', ['time', 'lon', 'lat'], 'S'],
                  # units=mho=S
                  'epot': ['V', 'Electric Potential',
                           4, 'SM', 'sph', ['time', 'lon', 'lat'], 'V'],
                  'int_joule_heat_n': ['W_JouleN', 'Northern Hemisphere ' +
                                       'Integrated Joule Heating',
                                       5, 'SM', 'sph', ['time'], 'GW'],
                  'int_joule_heat_s': ['W_JouleS', 'Southern Hemisphere ' +
                                       'Integrated Joule Heating',
                                       6, 'SM', 'sph', ['time'], 'GW'],
                  # 1D time series only, but different for each hemisphere
                  'jfac': ['j_fac', 'Field Aligned Current',
                           7, 'SM', 'sph', ['time', 'lon', 'lat'], 'muA/m**2'],
                  'joule_heat': ['Q_Joule',
                                 'Joule Heating (E-field^2*Pedersen)',
                                 8, 'SM', 'sph', ['time', 'lon', 'lat'],
                                 'mW/m**2'],
                  'mpot': ['psi', 'Magnetic Potential',
                           9, 'SM', 'sph', ['time', 'lon', 'lat'], 'cT/m'],
                  'sdB_ph': ['dB_east', 'Spacecraft-Observed Magnetic ' +
                             'Perturbations (eastward)',
                             10, 'SM', 'sph', ['time', 'lon', 'lat'], 'nT'],
                  'sdB_th': ['dB_north', 'Spacecraft-Observed Magnetic ' +
                             'Perturbations (equatorward)',
                             11, 'SM', 'sph', ['time', 'lon', 'lat'], 'nT'],
                  'v_ph': ['v_ieast', 'Ion Drift Velocity (eastward)',
                           12, 'SM', 'sph', ['time', 'lon', 'lat'], 'm/s'],
                  'v_th': ['v_inorth', 'Ion Drift Velocity (equatorward)',
                           13, 'SM', 'sph', ['time', 'lon', 'lat'], 'm/s'],
                  # in attrs of time dataset
                  'imf_By': ['B_y', 'measured y component of IMF magnetic ' +
                             'field from OMNI',
                             14, 'SM', 'sph', ['time'], 'nT'],
                  # in attrs of time dataset
                  'imf_Bz': ['B_z', 'measured z component of IMF magnetic ' +
                             'field from OMNI',
                             15, 'SM', 'sph', ['time'], 'nT'],
                  # in attrs of time dataset
                  'solar_wind_speed': ['v_sw', 'measured solar wind speed ' +
                                       'from OMNI',
                                       16, 'SM', 'sph', ['time'], 'km/s']
                  }


@vectorize
def timestr_hrs(time_str, filedate):
    '''Converts time string from data file into hours since midnight UTC.'''
    dt = datetime.strptime(time_str,
                           '%Y%m%d_%H%M%S').replace(tzinfo=timezone.utc)
    return (dt-filedate).total_seconds()/3600.


def timestr_datetime(time_str):
    '''Converts time string into a datetime object at midnight.'''
    return datetime.strptime(time_str[:8],
                             '%Y%m%d').replace(tzinfo=timezone.utc)


def timestr_utcts(time_str):
    '''Converts time string into utc timestamp.'''
    return datetime.strptime(time_str, '%Y%m%d_%H%M%S').replace(
        tzinfo=timezone.utc).timestamp()


def timestr_datetimestr(time_str):
    '''Converts time string into standard format (‘YYYY-MM-DD HH:MM:SS’).'''
    dt = datetime.strptime(time_str,
                           '%Y%m%d_%H%M%S').replace(tzinfo=timezone.utc)
    return dt.isoformat(sep=' ')[:19]


def ts_to_hrs(time_val, filedate):
    '''Convert utc timestamp to hours since midnight on filedate.'''
    return (datetime.utcfromtimestamp(time_val).replace(tzinfo=timezone.utc) -
            filedate).total_seconds()/3600.


def MODEL():
    from kamodo import Kamodo
    from glob import glob
    import h5py
    from os.path import basename
    from numpy import array, NaN, unique, append, zeros, abs, diff, sin, cos
    from numpy import where, flip, concatenate, insert, mean, broadcast_to
    from numpy import transpose
    from numpy import pi as nppi
    from numpy import sum as npsum
    from time import perf_counter
    from astropy.constants import R_earth
    from kamodo_ccmc.readers.reader_utilities import regdef_1D_interpolators
    from kamodo_ccmc.readers.reader_utilities import regdef_3D_interpolators

    class MODEL(Kamodo):
        '''AmGEO model data reader.

        Inputs:
            full_filenameN: a string representing the file pattern of the model
                output data.
                Note: This reader takes a file pattern of the format
                file_dir+*YYYYMMDDN.h5, where file_dir is the complete file
                path to the data files, and YYYYMMDD is the four digit year,
                two digit month, and two digit day in the desired output file
                names (e.g. *20150622N.h5).
            variables_requested = a list of variable name strings chosen from
                the model_varnames dictionary in this script, specifically the
                first item in the list associated with a given key.
                - If empty, the reader functionalizes all possible variables
                    (default)
                - If 'all', the reader returns the model_varnames dictionary
                    above for only the variables present in the given files.
                    Note: the fulltime keyword must be False to acheive this
                    behavior.
            filetime = boolean (default = False)
                - if False, the script fully executes.
                - If True, the script only executes far enough to determine the
                    time values associated with the chosen data.
                Note: The behavior of the script is determined jointly by the
                    filetime and fulltime keyword values.
            printfiles = boolean (default = False)
                - If False, the filenames associated with the data retrieved
                    ARE NOT printed.
                - If True, the filenames associated with the data retrieved ARE
                    printed.
            gridded_int = boolean (default = True)
                - If True, the variables chosen are functionalized in both the
                    standard method and a gridded method.
                - If False, the variables chosen are functionalized in only the
                    standard method.
            fulltime = boolean (default = True)
                - If True, linear interpolation in time between files is
                    included in the returned interpolator functions.
                - If False, no linear interpolation in time between files is
                    included.
            verbose = boolean (False)
                - If False, script execution and the underlying Kamodo
                    execution is quiet except for specified messages.
                - If True, be prepared for a plethora of messages.
        All inputs are described in further detail in
            KamodoOnboardingInstructions.pdf.

        Returns: a kamodo object (see Kamodo core documentation) containing all
            requested variables in functionalized form.
        '''
        def __init__(self, file_pattern, variables_requested=[],
                     printfiles=False, filetime=False, gridded_int=True,
                     fulltime=True, verbose=False, hemi='?', **kwargs):
            super(MODEL, self).__init__(**kwargs)
            self.modelname = 'AMGeO'
            t0 = perf_counter()

            # collect filenames
            files = sorted(glob(file_pattern+hemi+'.h5'))
            if 'N.h5' in files[0]:
                full_filenameN = files[0]
                if len(files) == 2:  # S file should be next
                    full_filenameS = files[1]
                    hemi = '*'  # data for both hemispheres available
                else:
                    full_filenameS = ''  # file DNE
                    hemi = 'N'  # only N hemisphere data found
                filename = basename(full_filenameN)
                file_dir = full_filenameN.split(filename)[0]
            elif 'S.h5' in files[0]:
                full_filenameS = files[0]
                full_filenameN = ''  # file DNE, N file should have been first
                hemi = 'S'  # only S hemisphere data found
                filename = basename(full_filenameS)
                file_dir = full_filenameS.split(filename)[0]
            self.filename = full_filenameN+', '+full_filenameS

            # establish time attributes first
            if full_filenameN != '':
                f_data1 = h5py.File(full_filenameN, 'r')
            else:
                f_data1 = h5py.File(full_filenameS, 'r')
            # remove 'N' or 'S' at end of time
            time_list = [key[:-1] for key in f_data1.keys() if key not in
                         ['lats', 'lons']]
            # datetime object for midnight on date
            self.filedate = timestr_datetime(time_list[0])
            # convert to hours since midnight of file
            time = timestr_hrs(time_list, self.filedate)
            self.datetimes = [timestr_datetimestr(time_list[0]),
                              timestr_datetimestr(time_list[-1])]
            self.filetimes = [timestr_utcts(time_list[0]),
                              timestr_utcts(time_list[-1])]
            self.dt = diff(time).max()*3600.  # convert time resolution to s

            if filetime and not fulltime:
                return  # return times as is to prevent recursion

            # if variables are given as integers, convert to standard names
            if len(variables_requested) > 0:
                if isinstance(variables_requested[0], int):
                    tmp_var = [value[0] for key, value in
                               model_varnames.items() if value[2] in
                               variables_requested]
                    variables_requested = tmp_var

            if fulltime:  # add boundary time (default value)
                # find other files with same pattern
                if hemi != '*':
                    file_pattern = file_dir+'*'+hemi+'.h5'  # returns a string
                else:
                    file_pattern = file_dir+'*.h5'
                files = sorted(glob(file_pattern))  # method may change for AWS
                filenames = unique([basename(f) for f in files])

                # find closest file by utc timestamp
                # files are sorted by YYMMDD, so next file is next in the list
                current_idx = where(filenames == filename)[0]
                if current_idx+1 == len(files):
                    if verbose:
                        print('No later file available.')
                    filecheck = False
                    if filetime:
                        return
                else:
                    min_file = file_dir+'*'+filenames[current_idx+1][0][-12:-4]
                    kamodo_test = MODEL(min_file, filetime=True,
                                        fulltime=False, hemi=hemi)
                    time_test = abs(kamodo_test.filetimes[0]-self.filetimes[1])
                    # if nearest file time at least within one timestep (hrs)
                    if time_test <= self.dt:
                        filecheck = True
                        self.datetimes[1] = kamodo_test.datetimes[0]
                        self.filetimes[1] = kamodo_test.filetimes[0]

                        # time only version if returning time for searching
                        if filetime:
                            return  # return object with additional time

                        # get kamodo object with same requested variables
                        if verbose:
                            print(f'Took {perf_counter()-t0:.3f}s to find ' +
                                  'closest file.')
                        kamodo_neighbor = MODEL(
                            min_file, variables_requested=variables_requested,
                            fulltime=False, hemi=hemi)
                        short_data = kamodo_neighbor.short_data
                        if verbose:
                            print(f'Took {perf_counter()-t0:.3f}s to get ' +
                                  'data from closest file.')
                    else:
                        if verbose:
                            print('No later file found within ' +
                                  f'{diff(time).max()*3600.:.1f}s.')
                        filecheck = False
                        if filetime:
                            return

            # perform initial check on variables_requested list
            if len(variables_requested) > 0 and fulltime and \
                    variables_requested != 'all':
                test_list = [value[0] for key, value in model_varnames.items()]
                err_list = [item for item in variables_requested if item
                            not in test_list]
                if len(err_list) > 0:
                    print('Variable name(s) not recognized:', err_list)

            # collect variable list (both in attributes and in datasets)
            if hemi in ['*', 'N']:
                key_list = list(f_data1[time_list[0]+'N'].attrs.keys()) + \
                    list(f_data1[time_list[0]+'N'].keys())
            elif hemi == 'S':
                key_list = list(f_data1[time_list[0]+'S'].attrs.keys()) + \
                    list(f_data1[time_list[0]+'S'].keys())
            if 'int_joule_heat' in key_list:  # replace with separate names
                key_list.remove('int_joule_heat')
                if hemi in ['*', 'N']:
                    key_list.append('int_joule_heat_n')
                if hemi in ['*', 'S']:
                    key_list.append('int_joule_heat_s')
            if len(variables_requested) > 0 and variables_requested != 'all':
                gvar_list = [key for key, value in model_varnames.items()
                             if value[0] in variables_requested and
                             key in key_list]  # file variable names

                # check for variables requested but not available
                if len(gvar_list) != len(variables_requested):
                    err_list = [value[0] for key, value in
                                model_varnames.items() if value[0] in
                                variables_requested and key not in gvar_list]
                    if len(err_list) > 0:
                        print('Some requested variables are not available:',
                              err_list)
            else:  # only input variables on the avoid_list if requested
                gvar_list = [key for key in key_list if key in
                             model_varnames.keys()]
                # return list of variables included in data files
                if not fulltime and variables_requested == 'all':
                    self.var_dict = {value[0]: value[1:] for key, value in
                                     model_varnames.items()
                                     if key in gvar_list}
                    return

            # store data for each variable desired
            if hemi == '*':
                f_south = h5py.File(full_filenameS, 'r')
                s = 'N'
            else:
                s = hemi
            variables = {model_varnames[var][0]:
                         {'units': model_varnames[var][-1],
                          'data': 0.} for var in gvar_list}
            # list of variables to be retrieved from attrs
            attrs_var = ['imf_By', 'imf_Bz', 'solar_wind_speed']
            for var in gvar_list:
                if var in attrs_var:  # save time series data from attributes
                    variables[model_varnames[var][0]]['data'] = \
                        array([f_data1[time+s].attrs[var] for time in
                               time_list], dtype=float)
                elif var == 'int_joule_heat_n' and hemi in ['*', 'N']:
                    # only time series variable in north dataset
                    variables[model_varnames[var][0]]['data'] = \
                        array([array(f_data1[time+'N']['int_joule_heat'])[0]
                               for time in time_list], dtype=float)
                elif var == 'int_joule_heat_s' and hemi in ['*', 'S']:
                    # only time series variable in south dataset
                    if hemi == '*':
                        variables[model_varnames[var][0]]['data'] = \
                            array([array(f_south[time+'S'][
                                'int_joule_heat'])[0] for time in time_list],
                                dtype=float)
                    elif hemi == 'S':
                        variables[model_varnames[var][0]]['data'] = \
                            array([array(f_data1[time+'S'][
                                'int_joule_heat'])[0] for time in time_list],
                                dtype=float)
                else:  # pull from datasets
                    if hemi == '*':
                        tmp = array([array(f_data1[time+'N'][var]) for
                                     time in time_list], dtype=float,
                                    order='F')
                        tmp1 = transpose(tmp, (0, 2, 1))  # now t, lon, lat
                        data1 = flip(tmp1, axis=2)
                        tmp = array([flip(array(f_south[time+'S'][var]),
                                          axis=0) for time in time_list],
                                    dtype=float, order='F')
                        tmp1 = transpose(tmp, (0, 2, 1))
                        south_data = flip(tmp1, axis=2)
                    elif hemi in ['N', 'S']:
                        tmp = array([array(f_data1[time+hemi][var]) for
                                     time in time_list], dtype=float,
                                    order='F')
                        tmp1 = transpose(tmp, (0, 2, 1))
                        data1 = flip(tmp1, axis=2)

                    # need to reverse order along lat axis for southern hemi
                    if '_th' in var and hemi != 'S':
                        # change from equatorward to northward in N hemi data
                        data1 *= -1
                    if hemi == '*':
                        total_data = concatenate((south_data, data1), axis=2)
                    else:
                        total_data = data1
                    variables[model_varnames[var][0]]['data'] = total_data

            # prepare and return data
            if not fulltime:
                f_data1.close()
                if hemi == '*':
                    f_south.close()
                variables['time'] = self.filetimes[0]
                self.short_data = variables
                return

            # Store coordinate data as class attributes
            if filecheck:  # add new time in hours since midnight
                new_time = ts_to_hrs(short_data['time'], self.filedate)
                self._time = append(time, new_time)
            else:
                self._time = time

            # collect and arrange lat grid to be increasing (from neg to pos)
            if hemi in ['*', 'N']:
                lat_N = unique(f_data1['lats'])  # 24 values
                if hemi == '*':
                    # reverse order and make negative
                    lat_S = flip(unique(f_south['lats']))*(-1)
                    # -88.?? to +88.?? (south pole first, north pole last)
                    lat = append(lat_S, lat_N)
                    lat = insert(lat, 0, -90.)  # insert south pole value
                    self._lat = append(lat, 90.)  # append north pole value
                else:
                    self._lat = append(lat_N, 90.)  # append north pole value
            elif hemi == 'S':  # N hemisphere file DNE
                # reverse order and make negative
                lat_S = flip(unique(f_data1['lats']))*(-1)
                self._lat = insert(lat_S, 0, -90.)  # insert south pole value
            self._lon = unique(f_data1['lons'])-180.  # shift noon to 0 deg lon

            # convert height in km to radius in R_E to conform to MAG coord sys
            # 110 km altitude
            self._radius = (110.+R_earth.value/1000.)/(R_earth.value/1000.)
            f_data1.close()   # close files
            if hemi == '*':
                f_south.close()

            # store a few items in object
            self.missing_value = NaN
            self._registered = 0
            if verbose:
                print(f'Took {perf_counter()-t0:.6f}s to read in data')
            if printfiles:
                print(self.filename)

            # need latitude wrapping (scalars and vectors), then
            #  register interpolators for each requested variable
            t_reg = perf_counter()
            varname_list = [key for key in variables.keys()]
            self.variables = {}
            for varname in varname_list:
                if len(variables[varname]['data'].shape) == 1:
                    if filecheck:  # if neighbor found
                        # append data for last time stamp, transpose
                        data_shape = list(variables[varname]['data'].shape)
                        data_shape[0] += 1  # add space for time
                        new_data = zeros(data_shape)
                        # put in current data
                        new_data[:-1] = variables[varname]['data']
                        # add in data for additional time
                        new_data[-1] = short_data[varname]['data'][0]
                        variable = new_data  # no transposing needed
                    else:
                        variable = variables[varname]['data']
                    self.variables[varname] = dict(
                        units=variables[varname]['units'], data=variable)
                    self.register_1D_variable(self.variables[varname]['units'],
                                              self.variables[varname]['data'],
                                              varname, gridded_int)
                elif len(variables[varname]['data'].shape) == 3:
                    if filecheck:
                        # append data for last time stamp, transpose
                        data_shape = list(variables[varname]['data'].shape)
                        data_shape[0] += 1  # add space for time
                        new_data = zeros(data_shape)
                        # put in current data
                        new_data[:-1, :, :] = variables[varname]['data']
                        # add in data for additional time
                        new_data[-1, :, :] = \
                            short_data[varname]['data'][0, :, :]
                        variable = new_data  # no transposing needed
                    else:
                        variable = variables[varname]['data']
                    self.variables[varname] = dict(
                        units=variables[varname]['units'], data=variable)
                    self.register_3D_variable(self.variables[varname]['units'],
                                              self.variables[varname]['data'],
                                              varname, gridded_int, hemi)
            if verbose:
                print(f'Took {perf_counter()-t_reg:.5f}s to register ' +
                      f'{len(varname_list)} variables.')
            if verbose:
                print(f'Took a total of {perf_counter()-t0:.5f}s to kamodofy' +
                      f' {len(varname_list)} variables.')

        def wrap3Dlat(self, varname, variable, hemi):
            '''Wrap latitude values for scalars and vectors.'''

            shape_list = list(variable.shape)  # time, lon, lat
            if hemi == '*':
                shape_list[2] += 2  # need two more places in latitude
            else:
                shape_list[2] += 1
            tmp_arr = zeros(shape_list)  # array to set-up wrapped data in
            if hemi == '*':
                tmp_arr[:, :, 1:-1] = variable  # copy data into grid
            elif hemi == 'N':
                tmp_arr[:, :, :-1] = variable  # N hemi lat appended
            elif hemi == 'S':
                tmp_arr[:, :, 1:] = variable  # S hemi lat inserted

            # wrapping in latitude for scalar variables
            if 'east' not in varname and 'north' not in varname:
                if hemi in ['*', 'S']:
                    # put in top values
                    top = mean(tmp_arr[:, :-1, 1], axis=1)  # same shape as t
                    tmp_arr[:, :-1, 0] = broadcast_to(top, (shape_list[1]-1,
                                                            shape_list[0])).T
                if hemi in ['*', 'N']:
                    # same for bottom, reusing variable names
                    top = mean(tmp_arr[:, :-1, -2], axis=1)  # same shape as t
                    tmp_arr[:, :-1, -1] = broadcast_to(top, (shape_list[1]-1,
                                                             shape_list[0])).T

            # wrapping in latitude for relevant vector variables
            elif 'east' in varname or 'north' in varname:
                if hemi in ['*', 'S']:
                    # calculate net vector magnitude for top
                    tmp_arr[:, :-1, 0] = self.vector_average3D(
                        tmp_arr[:, :-1, 1], shape_list, varname, self._lat[0])
                if hemi in ['*', 'N']:
                    # repeat for bottom
                    tmp_arr[:, :-1, -1] = self.vector_average3D(
                        tmp_arr[:, :-1, -2], shape_list, varname,
                        self._lat[-1])
            tmp_arr[:, -1, :] = tmp_arr[:, 0, :]  # wrap value in lon after
            self.variables[varname]['data'] = tmp_arr  # store result
            return tmp_arr

        def vector_average3D(self, top, shape_list, varname, latval):
            '''find vector average at pole for array with shape (time, lon)'''

            # find net x and y components, final array shapes are (time, lon)
            lon_arr = broadcast_to(self._lon[:-1], (shape_list[0],
                                                    shape_list[1]-1))
            # sum over lon axis, same shape as time
            xval = npsum(top*cos((lon_arr+180.)*nppi/180.), axis=1)
            yval = npsum(top*sin((lon_arr+180.)*nppi/180.), axis=1)
            # xval.shape must be last in broadcast_to call
            xarr = broadcast_to(xval, (shape_list[1]-1, shape_list[0])).T
            yarr = broadcast_to(yval, (shape_list[1]-1, shape_list[0])).T

            # convert to proper unit vector (see wiki on spherical coordinates)
            if 'east' in varname:
                # Zonal / east components -> convert to psi_hat vector (lon)
                # -xsin(psi)+ycos(psi), psi = longitude (0 to 360)
                new_top = -xarr*sin((lon_arr+180.)*nppi/180.) +\
                    yarr*cos((lon_arr+180.)*nppi/180.)
            elif 'north' in varname:
                # meridional/north -> convert to theta_hat vector (latitude)
                # xcos(psi)cos(theta)+ysin(psi)cos(theta)
                # sin(theta) is always zero at the poles
                # theta = latitude (0 to 180), psi = longitude (0 to 360)
                new_top = xarr * cos((lon_arr+180.)*nppi/180.) *\
                    cos((90.-latval)*nppi/180.) +\
                    yarr * sin((lon_arr+180.)*nppi/180.) *\
                    cos((90.-latval)*nppi/180.)
            return new_top

        # define and register a 1D variable
        def register_1D_variable(self, units, variable, varname, gridded_int):
            """Registers a 1d interpolator with 1d signature"""

            # define and register the interpolators
            xvec_dependencies = {'time': 'hr'}
            coord_str = [value[3]+value[4] for key, value in
                         model_varnames.items() if value[0] == varname][0]+'1D'
            self = regdef_1D_interpolators(self, units, variable, self._time,
                                           varname, xvec_dependencies,
                                           gridded_int, coord_str)
            return

        # define and register a 3D variable
        def register_3D_variable(self, units, variable, varname, gridded_int,
                                 hemi):
            """Registers a 3d interpolator with 3d signature"""

            # define and register the fast interpolator
            xvec_dependencies = {'time': 'hr', 'lon': 'deg', 'lat': 'deg'}
            variable = self.wrap3Dlat(varname, variable, hemi)
            coord_str = [value[3]+value[4] for key, value in
                         model_varnames.items() if value[0] == varname][0]+'3D'
            self = regdef_3D_interpolators(self, units, variable, self._time,
                                           self._lon, self._lat, varname,
                                           xvec_dependencies, gridded_int,
                                           coord_str)
            return
    return MODEL
