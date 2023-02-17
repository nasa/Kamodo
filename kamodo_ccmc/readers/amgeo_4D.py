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
                  'epot': ['phi', 'Electric Potential',
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


# list of variables found in h5_data.attrs.keys() instead of h5_data.keys()
attrs_var = ['imf_By', 'imf_Bz', 'solar_wind_speed']


def MODEL():
    from kamodo import Kamodo
    from glob import glob
    import h5py
    from os.path import isfile
    from numpy import array, NaN, unique, append, zeros, diff, sin, cos, insert
    from numpy import flip, concatenate, mean, broadcast_to, ones
    from numpy import transpose
    from numpy import pi as nppi
    from numpy import sum as npsum
    from time import perf_counter
    from astropy.constants import R_earth
    import kamodo_ccmc.readers.reader_utilities as RU

    class MODEL(Kamodo):
        '''AmGEO model data reader.

        Inputs:
            file_dir: a string representing the file directory of the
                model output data.
                Note: This reader 'walks' the entire dataset in the directory.
            variables_requested = a list of variable name strings chosen from
                the model_varnames dictionary in this script, specifically the
                first item in the list associated with a given key.
                - If empty, the reader functionalizes all possible variables
                    (default)
                - If 'all', the reader returns the model_varnames dictionary
                    above for only the variables present in the given files.
            filetime = boolean (default = False)
                - If False, the script fully executes.
                - If True, the script only executes far enough to determine the
                    time values associated with the chosen data.
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
            verbose = boolean (False)
                - If False, script execution and the underlying Kamodo
                    execution is quiet except for specified messages.
                - If True, be prepared for a plethora of messages.
        All inputs are described in further detail in
            KamodoOnboardingInstructions.pdf.

        Returns: a kamodo object (see Kamodo core documentation) containing all
            requested variables in functionalized form.

        Notes:
            - AMGeO model outputs are given in h5 files with one file per N/S
              hemisphere per day. The data from both hemispheres are assembled
              at each interpolation call.
            - The outputs do not provide values at the poles, so scalar and
              vector averaging are used as appropriate to determine these
              values.
            - The files are small and contain multiple time steps per file, so
              interpolation method 2 is chosen. The standard SciPy interpolator
              is used.
        '''
        def __init__(self, file_dir, variables_requested=[],
                     printfiles=False, filetime=False, gridded_int=True,
                     verbose=False, **kwargs):
            super(MODEL, self).__init__(**kwargs)
            self.modelname = 'AMGeO'
            t0 = perf_counter()

            # first, check for file list, create if DNE
            list_file = file_dir + self.modelname + '_list.txt'
            time_file = file_dir + self.modelname + '_times.txt'
            self.times, self.pattern_files = {}, {}
            if not isfile(list_file) or not isfile(time_file):
                # collect filenames
                files = sorted(glob(file_dir+'*.h5'))
                patterns = unique([f[-4] for f in files])  # N or S
                self.filename = ''.join([f+',' for f in files])[:-1]
                self.filedate = datetime.strptime(
                    files[0][-12:-4], '%Y%m%d').replace(tzinfo=timezone.utc)

                # establish time attributes
                for p in patterns:
                    # get list of files to loop through later
                    pattern_files = sorted(glob(file_dir+'*'+p+'.h5'))
                    self.pattern_files[p] = pattern_files
                    self.times[p] = {'start': [], 'end': [], 'all': []}

                    # loop through to get times
                    for f in range(len(pattern_files)):
                        h5_data = h5py.File(pattern_files[f])
                        tmp_var = [key[:-1] for key in h5_data.keys() if key
                                   not in ['lats', 'lons']]
                        h5_data.close()
                        tmp = RU.str_to_hrs(tmp_var, self.filedate,
                                            '%Y%m%d_%H%M%S')
                        self.times[p]['start'].append(tmp[0])
                        self.times[p]['end'].append(tmp[-1])
                        self.times[p]['all'].extend(tmp)
                    self.times[p]['start'] = array(self.times[p]['start'])
                    self.times[p]['end'] = array(self.times[p]['end'])
                    self.times[p]['all'] = array(self.times[p]['all'])

                # create time list file if DNE
                RU.create_timelist(list_file, time_file, self.modelname,
                                   self.times, self.pattern_files,
                                   self.filedate)
            else:  # read in data and time grids from file list
                self.times, self.pattern_files, self.filedate, self.filename =\
                    RU.read_timelist(time_file, list_file)
            if filetime:
                return  # return times only

            # perform initial check on variables_requested list
            if len(variables_requested) > 0 and variables_requested != 'all':
                test_list = [value[0] for key, value in model_varnames.items()]
                err_list = [item for item in variables_requested if item
                            not in test_list]
                if len(err_list) > 0:
                    print('Variable name(s) not recognized:', err_list)
                for item in err_list:
                    variables_requested.remove(item)
                if len(variables_requested) == 0:
                    return

            # collect variable list (in attributes of datasets)
            patterns = list(self.pattern_files.keys())
            h5_data = h5py.File(self.pattern_files[patterns[0]][0])
            key_list = list(h5_data[list(h5_data.keys())[0]].attrs.keys())
            key_list.extend(list(h5_data[list(h5_data.keys())[0]].keys()))
            h5_data.close()
            if 'int_joule_heat' in key_list:  # replace with separate names
                key_list.remove('int_joule_heat')
                if 'N' in patterns:
                    key_list.append('int_joule_heat_n')
                if 'S' in patterns:
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
                if variables_requested == 'all':
                    self.var_dict = {value[0]: value[1:] for key, value in
                                     model_varnames.items()
                                     if key in gvar_list}
                    return

            # initialize data mapping for each variable desired
            self.variables = {model_varnames[var][0]:
                              {'units': model_varnames[var][-1],
                               'data': var} for var in gvar_list}

            # collect and arrange lat grid to be increasing (from neg to pos)
            # add buffer rows for NaN over equator region
            if 'N' in patterns:
                h5_data = h5py.File(self.pattern_files['N'][0])
                lat_N = unique(h5_data['lats'])  # 24 values
                ldiff = diff(lat_N).min()
                lat_N = array([lat_N[0]-ldiff, lat_N[0]-ldiff/10.] +
                              list(lat_N) + [90.])
                if 'S' in patterns:
                    h5_dataS = h5py.File(self.pattern_files['S'][0])
                    # reverse order and make negative
                    lat_S = flip(unique(h5_dataS['lats']))*(-1)
                    lat_S = array([-90.] + list(lat_S) + [lat_S[-1]+ldiff/10.,
                                                          lat_S[-1]+ldiff])
                    h5_dataS.close()
                    # -88.?? to +88.?? (south pole first, north pole last)
                    self._lat = append(lat_S, lat_N)  # append north pole value
                else:
                    self._lat = lat_N  # append north pole value
            elif 'S' in patterns:  # N hemisphere file DNE
                h5_data = h5py.File(self.pattern_files['S'][0])
                # reverse order and make negative
                lat_S = flip(unique(h5_data['lats']))*(-1)
                ldiff = diff(lat_N).min()
                lat_S = array([-90.] + list(lat_S) + [lat_S[-1]+ldiff/10.,
                                                      lat_S[-1]+ldiff])
                self._lat = lat_S  # insert south pole value
            self._lon = unique(h5_data['lons'])-180.  # shift noon to 0 deg lon
            h5_data.close()

            # convert height in km to radius in R_E to conform to MAG coord sys
            # 110 km altitude
            self._radius = (110.+R_earth.value/1000.)/(R_earth.value/1000.)

            # store a few items in object
            if verbose:
                print(f'Took {perf_counter()-t0:.6f}s to read in data')
            if printfiles:
                print(self.filename)

            # need latitude wrapping (scalars and vectors), then
            #  register interpolators for each requested variable
            t_reg = perf_counter()
            varname_list = [key for key in self.variables.keys()]
            for varname in varname_list:
                self.register_variable(varname, gridded_int)
            if verbose:
                print(f'Took {perf_counter()-t_reg:.5f}s to register ' +
                      f'{len(varname_list)} variables.')
            if verbose:
                print(f'Took a total of {perf_counter()-t0:.5f}s to kamodofy' +
                      f' {len(varname_list)} variables.')

        def wrap3Dlat(self, varname, variable, ps):
            '''Wrap latitude values for scalars and vectors.'''

            shape_list = list(variable.shape)  # time, lon, lat
            if len(ps) == 2:
                shape_list[2] += 2  # need two more places in latitude
            else:
                shape_list[2] += 1  # pole row
            tmp_arr = zeros(shape_list)  # array to set-up wrapped data in
            if len(ps) == 2:
                tmp_arr[:, :, 1:-1] = variable  # copy data into grid
            elif ps[0] == 'N':
                tmp_arr[:, :, :-1] = variable  # N hemi lat appended
            elif ps[0] == 'S':
                tmp_arr[:, :, 1:] = variable  # S hemi lat inserted

            # wrapping in latitude for scalar variables
            if 'east' not in varname and 'north' not in varname:
                if 'S' in ps:
                    # put in top values
                    top = mean(tmp_arr[:, :-1, 1], axis=1)  # same shape as t
                    tmp_arr[:, :-1, 0] = broadcast_to(top, (shape_list[1]-1,
                                                            shape_list[0])).T
                if 'N' in ps:
                    # same for bottom, reusing variable names
                    top = mean(tmp_arr[:, :-1, -2], axis=1)  # same shape as t
                    tmp_arr[:, :-1, -1] = broadcast_to(top, (shape_list[1]-1,
                                                             shape_list[0])).T

            # wrapping in latitude for relevant vector variables
            elif 'east' in varname or 'north' in varname:
                if 'S' in ps:
                    # calculate net vector magnitude for top
                    tmp_arr[:, :-1, 0] = self.vector_average3D(
                        tmp_arr[:, :-1, 1], shape_list, varname, self._lat[0])
                if 'N' in ps:
                    # repeat for bottom
                    tmp_arr[:, :-1, -1] = self.vector_average3D(
                        tmp_arr[:, :-1, -2], shape_list, varname,
                        self._lat[-1])
            tmp_arr[:, -1, :] = tmp_arr[:, 0, :]  # wrap value in lon after
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
        def register_variable(self, varname, gridded_int):
            """Creates and registers the interpolator"""

            # determine coordinate variables and xvec by coord list
            gvar = self.variables[varname]['data']
            ps = list(self.pattern_files.keys())
            coord_list = [value[5] for key, value in model_varnames.items()
                          if value[0] == varname][0]
            # time grids are the same in both hemispheres
            coord_dict = {'time': {'units': 'hr',
                                   'data': self.times[ps[0]]['all']}}
            if 'lon' in coord_list:
                coord_dict['lon'] = {'units': 'deg', 'data': self._lon}
            if 'lat' in coord_list:
                coord_dict['lat'] = {'units': 'deg', 'data': self._lat}
            coord_str = [value[3]+value[4] for key, value in
                         model_varnames.items() if value[0] == varname][0]

            # functionalize time series data, converting to an array first
            if len(coord_list) == 1:
                if gvar in attrs_var:
                    data = []
                    for f in self.pattern_files[ps[0]]:  # N, or S if N DNE
                        h5_data = h5py.File(f)
                        data.extend([h5_data[tkey].attrs[gvar] for tkey in
                                     h5_data.keys() if tkey not in
                                     ['lats', 'lons']])  # one float per time
                        h5_data.close()
                elif gvar[:-2] == 'int_joule_heat' and gvar[-1].upper() in ps:
                    data, p = [], gvar[-1].upper()
                    for f in self.pattern_files[p]:  # N, or S if N DNE
                        h5_data = h5py.File(f)
                        data.extend([h5_data[tkey][gvar[:-2]][0] for tkey in
                                     h5_data.keys() if tkey not in
                                     ['lats', 'lons']])  # one float per time
                        h5_data.close()
                self.variables[varname]['data'] = array(data)
                self = RU.Functionalize_Dataset(self, coord_dict, varname,
                                                self.variables[varname],
                                                gridded_int, coord_str,
                                                interp_flag=0)  # single array
                return

            # remaining logic is for 3D data, define func
            def func(i):  # i is the file number
                # get first hemisphere data
                if 'N' in ps:  # first one is N hemisphere
                    h5_data = h5py.File(self.pattern_files['N'][i])
                    tkeys = [key for key in h5_data.keys() if key not in
                             ['lats', 'lons']]
                    tmp = array([array(h5_data[tkey][gvar]) for tkey in
                                 tkeys], dtype=float, order='F')
                    # fill_value = h5_data[tkeys[0]][gvar].fillvalue
                    h5_data.close()
                    data = flip(transpose(tmp, (0, 2, 1)), axis=2)
                    NaN_row = ones(data[:, :, 0].shape) * NaN
                    data = insert(data, 0, data[:, :, 0], axis=2)  # add buffer
                    data = insert(data, 0, NaN_row, axis=2)
                    # change from equatorward to northward in N hemi data
                    if '_th' in varname:
                        data *= -1
                    # add S hemisphere data
                    if len(ps) == 2:
                        h5_data = h5py.File(self.pattern_files['S'][i])
                        tmp = array([array(h5_data[tkey][gvar]) for tkey in
                                     h5_data.keys() if tkey not in
                                     ['lats', 'lons']],
                                    dtype=float, order='F')
                        h5_data.close()
                        south_data = transpose(tmp, (0, 2, 1))
                        # add buffer row
                        south_data = insert(south_data, south_data.shape[2],
                                            south_data[:, :, -1], axis=2)
                        south_data = insert(south_data, south_data.shape[2],
                                            NaN_row, axis=2)
                        data = concatenate((south_data, data), axis=2)
                elif 'S' in ps:  # only S hemisphere data
                    h5_data = h5py.File(self.pattern_files['S'][i])
                    tkeys = [key for key in h5_data.keys() if key not in
                             ['lats', 'lons']]
                    tmp = array([array(h5_data[tkey][gvar]) for tkey in
                                 tkeys], dtype=float, order='F')
                    # fill_value = h5_data[tkeys[0]][gvar].fillvalue
                    h5_data.close()
                    data = transpose(tmp, (0, 2, 1))
                    NaN_row = ones(data[:, :, 0].shape) * NaN
                    # add buffer row
                    data = insert(data, data.shape[2], data[:, :, -1], axis=2)
                    data = insert(data, data.shape[2], NaN_row, axis=2)

                # add one time step from the next file if not the last file
                if i != len(self.pattern_files[ps[0]])-1:
                    if 'N' in ps:  # N hemisphere data is first
                        h5_data = h5py.File(self.pattern_files['N'][i+1])
                        time_list = list(h5_data.keys())[:-2]
                        tmp = array(h5_data[time_list[-1]][gvar], dtype=float,
                                    order='F')
                        h5_data.close()
                        data1 = flip(tmp.T, axis=1)  # shape = (lon, lat)
                        NaN_row1 = ones(data1[:, 0].shape) * NaN
                        # add buffer row
                        data1 = insert(data1, 0, data1[:, 0], axis=1)
                        data1 = insert(data1, 0, NaN_row1, axis=1)
                        # change from equatorward to northward in N hemi data
                        if '_th' in varname:
                            data1 *= -1
                        # add S hemisphere data
                        if len(ps) == 2:
                            h5_data = h5py.File(self.pattern_files[ps[1]][i+1])
                            time_list = list(h5_data.keys())[:-2]
                            tmp = array(h5_data[time_list[-1]][gvar],
                                        dtype=float, order='F')
                            h5_data.close()
                            south_data1 = tmp.T
                            south_data1 = insert(
                                south_data1, south_data1.shape[1],
                                south_data1[:, 0], axis=1)  # add buffer row
                            south_data1 = insert(
                                south_data1, south_data1.shape[1], NaN_row1,
                                axis=1)
                            data1 = concatenate((south_data1, data1), axis=1)
                    elif 'S' in ps:  # only S hemisphere data
                        h5_data = h5py.File(self.pattern_files['S'][i])
                        tmp = array(h5_data[time_list[-1]][gvar], dtype=float,
                                    order='F')
                        h5_data.close()
                        data1 = tmp.T
                        NaN_row1 = ones(data1[:, 0].shape) * NaN
                        data1 = insert(data1, data1.shape[1], data1[:, 0],
                                       axis=1)  # add buffer row
                        data1 = insert(data1, data1.shape[1], NaN_row1, axis=1)
                    data = append(data, [data1], axis=0)
                # No fill values were defined in sample dataset.
                # Not sure how this will affect in wrap3Dlat function.
                # if fill_value != 0.:# if fill value defined, replace with NaN
                #    data = where(data != fill_value, data, NaN)# none in test
                final_data = self.wrap3Dlat(varname, data, ps)
                return final_data

            self = RU.Functionalize_Dataset(
                self, coord_dict, varname, self.variables[varname],
                gridded_int, coord_str, interp_flag=2, func=func,
                times_dict=self.times[ps[0]])
            return
    return MODEL
