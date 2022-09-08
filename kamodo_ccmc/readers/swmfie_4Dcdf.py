# -*- coding: utf-8 -*-
"""
Created on Mon May 17 18:53:35 2021

@author: rringuet
"""
from datetime import datetime, timezone
from numpy import vectorize
from os.path import isfile, basename

model_varnames = {"Sigma_H": ['Sigma_H', '3D Hall conductivity',
                              0, 'SM', 'sph', ['time', 'lon', 'lat'], "S"],
                  "Sigma_P": ['Sigma_P', '3D Pederson conductivity',
                              1, 'SM', 'sph', ['time', 'lon', 'lat'], "S"],
                  "Phi_E": ['Phi_E', 'energy flux',
                            2, 'SM', 'sph', ['time', 'lon', 'lat'], "W/m**2"],
                  "AveE_avgE": ['E_avg', 'average energy',
                                3, 'SM', 'sph', ['time', 'lon', 'lat'], 'eV'],
                  "j_R": ["j_R", 'radial current density',
                          4, 'SM', 'sph', ['time', 'lon', 'lat'], "muA/m**2"],
                  "Phi": ["phi", 'electric potential',
                          5, 'SM', 'sph', ['time', 'lon', 'lat'], "kV"],
                  "E_x": ["E_x", 'electric field, x component',
                          6, 'SM', 'sph', ['time', 'lon', 'lat'], "mV/m"],
                  "E_y": ["E_y", 'electric field, y component',
                          7, 'SM', 'sph', ['time', 'lon', 'lat'], "mV/m"],
                  "E_z": ["E_z", 'electric field, z component',
                          8, 'SM', 'sph', ['time', 'lon', 'lat'], "mV/m"],
                  "j_x": ["j_x", 'current density, x component',
                          9, 'SM', 'sph', ['time', 'lon', 'lat'], "muA/m**2"],
                  "j_y": ["j_y", 'current density, y component',
                          10, 'SM', 'sph', ['time', 'lon', 'lat'], "muA/m**2"],
                  "j_z": ["j_z", 'current density, z component',
                          11, 'SM', 'sph', ['time', 'lon', 'lat'], "muA/m**2"],
                  "v_x": ['v_x', 'total velocity, x component',
                          12, 'SM', 'sph', ['time', 'lon', 'lat'], "km/s"],
                  "v_y": ['v_y', 'total velocity, y component',
                          13, 'SM', 'sph', ['time', 'lon', 'lat'], "km/s"],
                  "v_z": ['v_z', 'total velocity, z component',
                          14, 'SM', 'sph', ['time', 'lon', 'lat'], "km/s"],
                  "Q_Joule": ['W_JouleH', 'height integrated joule heating',
                              15, 'SM', 'sph', ['time', 'lon', 'lat'],
                              "mW/m**2"],
                  "Phi_nion": ['Phi_Nion', 'flux of ions in number density',
                               16, 'SM', 'sph', ['time', 'lon', 'lat'],
                               "1/cm**2/s"],
                  "Binv_RT": ['Binv_RT', 'inverse magnetic field (RT) ' +
                              '(ray tracing integrated along field line)',
                              17, 'SM', 'sph', ['time', 'lon', 'lat'], "1/T"],
                  "rho_RT": ['rho_RTamu', 'molecular mass density (RT) ' +
                             '(ray tracing integrated along field line)',
                             18, 'SM', 'sph', ['time', 'lon', 'lat'],
                             "amu/cm**3"],
                  "P_RT": ['P_RT', 'pressure (RT) ' +
                           '(ray tracing integrated along field line)',
                           19, 'SM', 'sph', ['time', 'lon', 'lat'], "Pa"],
                  "dLat_star": ['lat_star', 'conjugate latitude',
                                20, 'SM', 'sph', ['time', 'lon', 'lat'],
                                "deg"],
                  "dlon_star": ['lon_star', 'conjugate longitude',
                                21, 'SM', 'sph', ['time', 'lon', 'lat'],
                                "deg"]
                  }
'''
Documentation variables:
(ignored, given in R_E on a unit sphere)
"X":['x', 'km'], "Y":['y', 'km'], "Z":['z', 'km'],
"Theta":['theta',"deg"],"Psi":['psi',"deg"],     (used as coordinates)
"Btilt_theta":['theta_Btilt',"deg"], "Btilt_psi":['psi_Btilt',"deg"]
(added directly to object for documentation purposes)
'''


@vectorize
def dts_to_hrs(datetime_string, filedate):
    '''Get hours since midnight from datetime string'''
    return (datetime.strptime(datetime_string, '%Y-%m-%d %H:%M:%S').replace(
        tzinfo=timezone.utc)-filedate).total_seconds()/3600.


@vectorize
def filename_to_dts(filename, string_date):
    '''Get datetime string in format "YYYY-MM-SS HH:mm:SS" from filename'''
    mmhhss = basename(filename)[12:18]
    return string_date+' '+mmhhss[:2]+':'+mmhhss[2:4]+':'+mmhhss[4:]


def dts_to_ts(file_dts):
    '''Get datetime timestamp in UTC from datetime string'''
    return datetime.timestamp(datetime.strptime(file_dts, '%Y-%m-%d %H:%M:%S'
                                                ).replace(tzinfo=timezone.utc))


def ts_to_hrs(time_val, filedate):
    '''Convert utc timestamp to hours since midnight on filedate.'''
    return (datetime.utcfromtimestamp(time_val).replace(tzinfo=timezone.utc) -
            filedate).total_seconds()/3600.


def MODEL():
    from numpy import array, NaN, abs, unique, append, zeros, diff, where
    from time import perf_counter
    from netCDF4 import Dataset
    from kamodo import Kamodo
    from kamodo_ccmc.readers.reader_utilities import regdef_3D_interpolators

    class MODEL(Kamodo):
        '''SWMF_IE model data reader.

        Inputs:
            full_file_prefix: a string representing the file pattern of the
                model output data.
                Note: This reader takes a file pattern of the format
                file_dir+'i_eYYYYMMDD', where YYYY is the four digit year,
                MM is the two digit month, and DD is the two digit day
                (e.g. 20080502 for May 2, 2008).
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
        def __init__(self, full_file_prefix, variables_requested=[],
                     filetime=False, verbose=False, gridded_int=True,
                     printfiles=False, fulltime=True, **kwargs):
            super(MODEL, self).__init__()
            self.modelname = 'SWMF_IE'

            # check if given .nc file exists. If not, convert files to netCDF
            file_prefix = basename(full_file_prefix)
            file_dir = full_file_prefix.split(file_prefix)[0]
            if not isfile(full_file_prefix+'.nc'):
                from kamodo_ccmc.readers.swmfie_tocdf import \
                    convertSWMFIE_toCDF
                test = convertSWMFIE_toCDF(full_file_prefix)
                if not test:
                    self.conversion_test = test
                    return    # if file conversion fails, return
                else:
                    self.conversion_test = True
            else:
                self.conversion_test = True
            t0 = perf_counter()

            # determine type of prefix: for a day or for a hour
            if '-' in file_prefix:
                day_flag = False
            else:
                day_flag = True

            # establish time attributes first for file searching
            file_datestr = file_prefix[3:11]
            # string_date = 'YYYY-MM-DD'
            string_date = file_datestr[:4] + '-' + file_datestr[4:6] + '-' +\
                file_datestr[6:8]
            self.filedate = datetime.strptime(string_date+' 00:00:00',
                                              '%Y-%m-%d %H:%M:%S'
                                              ).replace(tzinfo=timezone.utc)

            # establish beginning and end time of file list
            cdf_data = Dataset(full_file_prefix+'.nc', 'r')
            files = cdf_data.file.split(',')
            # strings in format = YYYY-MM-DD HH:MM:SS
            self.datetimes = list(filename_to_dts([files[0], files[-1]],
                                                  string_date))
            # timestamps in UTC
            self.filetimes = [dts_to_ts(file_dts) for file_dts in
                              self.datetimes]
            t = array(cdf_data.variables['time'])  # hours since midnight
            if len(t) > 1:
                self.dt = diff(t).max()*3600.
            else:
                self.dt = 0

            if filetime and not fulltime:
                return  # return times as is to prevent recursion

            # if variables are given as integers, convert to standard names
            if len(variables_requested) > 0:
                if isinstance(variables_requested[0], int):
                    tmp_var = [value[0] for key, value in
                               model_varnames.items()
                               if value[2] in variables_requested]
                    variables_requested = tmp_var

            if fulltime:  # add boundary time (default value)
                # find other files with same pattern
                from glob import glob

                files = sorted(glob(file_dir+'i_e*'))
                if day_flag:
                    file_prefixes = unique([basename(f)[:11] for f in files
                                            if '.nc' not in basename(f)])
                else:  # give prefix for hourly files
                    file_prefixes = unique([basename(f)[:14] for f in files
                                            if '.nc' not in basename(f)])

                # find closest file by utc timestamp
                # swmf_ie has an open time at the end
                current_idx = where(file_prefixes == file_prefix)[0]
                if current_idx+1 == len(file_prefixes):
                    if verbose:
                        print('No later file available.')
                    filecheck = False
                    if filetime:
                        return
                else:
                    # +1 for adding an end time
                    min_file_prefix = file_dir+file_prefixes[current_idx+1][0]
                    kamodo_test = MODEL(min_file_prefix, filetime=True,
                                        fulltime=False)
                    if not kamodo_test.conversion_test:
                        if verbose:
                            print('No later file available.')
                        filecheck = False
                        if filetime:
                            return
                    else:
                        time_test = abs(kamodo_test.filetimes[0] -
                                        self.filetimes[1])
                        if time_test <= self.dt:
                            filecheck = True
                            self.datetimes[1] = kamodo_test.datetimes[0]
                            self.filetimes[1] = kamodo_test.filetimes[0]

                            # time only version if returning time for searching
                            if filetime:
                                return

                            # get kamodo object with same requested variables
                            if verbose:
                                print(f'Took {perf_counter()-t0:.3f}s to ' +
                                      'find closest file.')
                            kamodo_neighbor = MODEL(
                                min_file_prefix,
                                variables_requested=variables_requested,
                                fulltime=False)
                            short_data = kamodo_neighbor.short_data
                            if verbose:
                                print(f'Took {perf_counter()-t0:.3f}s to get' +
                                      ' data from closest file.')
                        else:
                            if verbose:
                                print('No later file found within ' +
                                      f'{self.dt:.1f}s.')
                            filecheck = False
                            if filetime:
                                return

            # perform initial check on variables_requested list
            if len(variables_requested) > 0 and fulltime and\
                    variables_requested != 'all':
                test_list = [value[0] for key, value in model_varnames.items()]
                err_list = [item for item in variables_requested if item not in
                            test_list]
                if len(err_list) > 0:
                    print('Variable name(s) not recognized:', err_list)

            # get list of variables possible in these files using first file
            if len(variables_requested) > 0 and variables_requested != 'all':
                gvar_list = [key for key, value in model_varnames.items()
                             if value[0] in variables_requested and
                             key in cdf_data.variables.keys()]

                # check for variables requested but not available
                if len(gvar_list) != len(variables_requested):
                    err_list = [value[0] for key, value in
                                model_varnames.items()
                                if value[0] in variables_requested and
                                key not in cdf_data.variables.keys()]
                    if len(err_list) > 0:
                        print('Some requested variables are not available:',
                              err_list)
            else:
                avoid_list = ['theta_Btilt', 'psi_Btilt']
                gvar_list = [key for key in cdf_data.variables.keys()
                             if key in model_varnames.keys() and
                             key not in avoid_list]
                if not fulltime and variables_requested == 'all':
                    self.var_dict = {value[0]: value[1:] for key, value in
                                     model_varnames.items() if key in
                                     gvar_list}
                    return

            # Store variable's data and units, transposing the 2D+time array.
            variables = {model_varnames[key][0]: {
                'units': model_varnames[key][-1],
                'data': array(cdf_data.variables[key])}
                for key in gvar_list}
            self.theta_Btilt = array(cdf_data.variables['theta_Btilt'])
            self.psi_Btilt = array(cdf_data.variables['psi_Btilt'])

            # prepare and return data only
            if not fulltime:
                cdf_data.close()
                variables['time'] = self.filetimes[0]
                variables['theta_Btilt'] = self.theta_Btilt[0]
                variables['psi_Btilt'] = self.psi_Btilt[0]
                self.short_data = variables
                return

            # return if only one file found - interpolator code will break
            if len(files) < 2:
                print('Not enough files found with given file prefix.')
                return

            # store variables
            self.filename = files
            self.missing_value = NaN
            self.modelname = 'SWMF_IE'
            self._registered = 0
            if printfiles:
                print('Files:')
                for file in self.filename:
                    print(file)

            # Store coordinate data as class attributes
            if filecheck:
                # new time in hours since midnight
                new_time = ts_to_hrs(short_data['time'], self.filedate)
                self._time = append(t, new_time)
            else:
                self._time = t

            # store coordinate data
            self._lat = array(cdf_data.variables['lat'])
            self._lon = array(cdf_data.variables['lon'])
            cdf_data.close()
            if verbose:
                print(f'Took {perf_counter()-t0:.6f}s to read in data')

            # register interpolators for each variable
            # store original list b/c gridded interpolators
            varname_list = [key for key in variables.keys()]
            self.variables = {}
            t_reg = perf_counter()
            for varname in varname_list:  # all are 3D variables
                if filecheck:  # if neighbor found
                    # append data for first time stamp, transpose
                    data_shape = list(variables[varname]['data'].shape)
                    data_shape[0] += 1  # add space for time
                    new_data = zeros(data_shape)
                    # put in current data
                    new_data[:-1, :, :] = variables[varname]['data']
                    # add in data for additional time
                    new_data[-1, :, :] = short_data[varname]['data'][0, :, :]
                else:
                    new_data = variables[varname]['data']
                self.variables[varname] = dict(
                    units=variables[varname]['units'], data=new_data)
                self.register_3D_variable(self.variables[varname]['units'],
                                          self.variables[varname]['data'],
                                          varname, gridded_int)
            if verbose:
                print(f'Took {perf_counter()-t_reg:.5f}s to register ' +
                      f'{len(varname_list)} variables.')
            if verbose:
                print(f'Took a total of {perf_counter()-t0:.5f}s to kamodofy' +
                      f' {len(gvar_list)} variables.')

        # define and register a 3D variable------------------------------------
        def register_3D_variable(self, units, variable, varname, gridded_int):
            """Registers a 3d interpolator with 3d signature"""

            # define and register the interpolators
            xvec_dependencies = {'time': 'hr', 'lon': 'deg', 'lat': 'deg'}
            self = regdef_3D_interpolators(self, units, variable, self._time,
                                           self._lon, self._lat, varname,
                                           xvec_dependencies, gridded_int)
            return
    return MODEL
