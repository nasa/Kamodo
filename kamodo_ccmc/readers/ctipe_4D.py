'''
Spatial-only interpolation form:  Lutz Raestatter
Modified from spatial-only interpolation into this form:  Rebecca Ringuette
Kamodofication of the CTIPe model output

To check for pep8 standards:
    - pip install pep8  (or pycodestyle)
    - pep8 --first ctipe_4D.py   #to only show first instance of error type
    - pep8 --show-source ctipe_4D.py  #to show full text per error
https://pypi.org/project/pep8/ for more info.
'''

from datetime import datetime,  timezone
from numpy import vectorize


# constants and dictionaries
model_varnames = {'rho': ['rho', 'total mass density',
                          0, 'GDZ_plev', 'sph', ['time', 'lon', 'lat', 'ilev1'], 'kg/m**3'],
                  'T': ['T', 'temperature',
                        1, 'GDZ_plev', 'sph', ['time', 'lon', 'lat', 'ilev1'], 'K'],
                  'T_e': ['T_e', 'electron temperature',
                          2, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'], 'K'],
                  'T_i': ['T_i', 'ion temperature',
                          3, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'], 'K'],
                  'H_ilev': ['H_ilev', 'height dependent on primary pressure level',
                             4, 'GDZ_plev', 'sph', ['time', 'lon', 'lat', 'ilev'], 'm'],
                  'H_lev': ['H_ilev1', 'height dependent on secondary pressure level',
                            5, 'GDZ_plev', 'sph', ['time', 'lon', 'lat', 'ilev1'], 'm'],
                  'Vn_lat': ['v_nnorth', 'meridional neutral wind velocity (north)',
                             6, 'GDZ_plev', 'sph', ['time', 'lon', 'lat', 'ilev'], 'm/s'],
                  'Vn_lon': ['v_neast', 'zonal neutral wind velocity (east)',
                             7, 'GDZ_plev', 'sph', ['time', 'lon', 'lat', 'ilev'], 'm/s'],
                  'Vn_H': ['v_nup', 'vertical neutral wind velocity (up)',
                           8, 'GDZ_plev', 'sph', ['time', 'lon', 'lat', 'ilev'], 'm/s'],
                  'T_n': ['T_n', 'neutral temperature',
                          9, 'GDZ_plev', 'sph', ['time', 'lon', 'lat', 'ilev'], 'K'],
                  'Rmt': ['m_avgmol', 'mean molecular mass',
                          10, 'GDZ_plev', 'sph', ['time', 'lon', 'lat', 'ilev1'], 'amu'],
                  'N_e': ['N_e', 'electron number density',
                          11, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'], '1/m**3'],
                  # 'N_n': ['N_n', 'variable description',
                  # 12, 'GDZ_plev', 'sph', ['time', 'lon', 'lat', 'ilev'], '1/m**3'],
                  'Q_Solar': ['Q_Solar', 'solar heating',
                              13, 'GDZ_plev', 'sph', ['time', 'lon', 'lat', 'ilev'], 'J/kg/s'],
                  'Q_Joule': ['Q_Joule', 'joule heating',
                              14, 'GDZ_plev', 'sph', ['time', 'lon', 'lat', 'ilev'], 'J/kg/s'],
                  'Q_radiation': ['Q_rad', 'radiative heating or cooling',
                                  15, 'GDZ_plev', 'sph', ['time', 'lon', 'lat', 'ilev'], 'J/kg/s'],
                  'N_O': ['N_O', 'number density of atomic oxygen',
                          16, 'GDZ_plev', 'sph', ['time', 'lon', 'lat', 'ilev'], '1/m**3'],
                  'N_O2': ['N_O2', 'number density of molecular oxygen',
                           17, 'GDZ_plev', 'sph', ['time', 'lon', 'lat', 'ilev'], '1/m**3'],
                  'N_N2': ['N_N2', 'number density of molecular nitrogen',
                           18, 'GDZ_plev', 'sph', ['time', 'lon', 'lat', 'ilev'], '1/m**3'],
                  'N_NO': ['N_NO', 'number density of molecular nitric oxide',
                           19, 'GDZ_plev', 'sph', ['time', 'lon', 'lat', 'ilev'], '1/m**3'],
                  'N_NOplus': ['N_NOplus', 'number density of nitric oxide ion',
                               20, 'GDZ_plev', 'sph', ['time', 'lon', 'lat', 'ilev'], '1/m**3'],
                  'N_N2plus': ['N_N2plus', 'number density of molecular nitrogen ion',
                               21, 'GDZ_plev', 'sph', ['time', 'lon', 'lat', 'ilev'], '1/m**3'],
                  'N_O2plus': ['N_O2plus', 'number density of molecular oxygen ion',
                               22, 'GDZ_plev', 'sph', ['time', 'lon', 'lat', 'ilev'], '1/m**3'],
                  'N_Nplus': ['N_Nplus', 'number density of atomic nitrogen ion',
                              23, 'GDZ_plev', 'sph', ['time', 'lon', 'lat', 'ilev'], '1/m**3'],
                  'N_Oplus': ['N_Oplus', 'number density of atomic oxygen ion',
                              24, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'], '1/m**3'],
                  'N_Hplus': ['N_Hplus', 'number density of atomic hydrogen ion',
                              25, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'], '1/m**3'],
                  'Sigma_P': ['sigma_P', 'Pedersen conductivity',
                              26, 'GDZ_plev', 'sph', ['time', 'lon', 'lat', 'ilev'], 'S/m'],
                  'Sigma_H': ['sigma_H', 'Hall conductivity',
                              27, 'GDZ_plev', 'sph', ['time', 'lon', 'lat', 'ilev'], 'S/m'],
                  'Vi_lon': ['v_inorth', 'meridional ion wind velocity (north)',
                             28, 'GDZ_plev', 'sph', ['time', 'lon', 'lat', 'ilev'], 'm/s'],
                  'Vi_lat': ['v_ieast', 'zonal ion wind velocity (east)',
                             29, 'GDZ_plev', 'sph', ['time', 'lon', 'lat', 'ilev'], 'm/s'],
                  'W_Joule': ['W_JouleH', 'height integrated joule heating',
                              30, 'GDZ', 'sph', ['time', 'lon', 'lat'], 'W/m**2'],
                  'Eflux_precip': ['Phi_E', 'energy flux',
                                   31, 'GDZ', 'sph', ['time', 'lon', 'lat'], 'mW/m**2'],
                  'Eavg_precip': ['E_avg', 'average energy',
                                  32, 'GDZ', 'sph', ['time', 'lon', 'lat'], 'keV'],
                  'TEC': ['TEC', 'vertical total electron content (height integrated from bottom to top boundary)',
                          33, 'GDZ', 'sph', ['time', 'lon', 'lat'], '1/m**2'],  # '10**16/m**2'
                  'E_theta140km': ['E_theta140km', 'Electric field at 140 km,  theta component',
                                   34, 'GDZ_E', 'sph', ['time', 'Elon', 'Elat'], 'V/m'],
                  'E_lambda140km': ['E_lambda140km', 'Electric field at 140 km, lambda component',
                                    35, 'GDZ_E', 'sph', ['time', 'Elon', 'Elat'], 'V/m'],
                  'E_theta300km': ['E_theta300km', 'Electric field at 300 km, theta component',
                                   36, 'GDZ_E', 'sph', ['time', 'Elon', 'Elat'], 'V/m'],
                  'E_lambda300km': ['E_lambda300km', 'Electric field at 300 km, lambda component',
                                    37, 'GDZ_E', 'sph', ['time', 'Elon', 'Elat'], 'V/m']}


@vectorize
def ts_to_hrs(time_val, filedate):
    '''Convert utc timestamp to hours since midnight on filedate.'''

    return (datetime.utcfromtimestamp(time_val).replace(tzinfo=timezone.utc) -
            filedate).total_seconds()/3600.


def MODEL():
    from numpy import array, zeros, abs, NaN, unique, insert, diff, where
    from time import perf_counter
    from os.path import isfile, basename
    from kamodo import Kamodo
    from netCDF4 import Dataset
    from kamodo_ccmc.readers.reader_utilities import regdef_4D_interpolators
    from kamodo_ccmc.readers.reader_utilities import regdef_3D_interpolators

    # main class
    class MODEL(Kamodo):
        '''CTIPe model data reader.

        Inputs:
            full_file_prefix:  a string representing the file pattern of the model
                output data.
                Note:  This reader takes a file pattern of the format
                file_dir+YYYY-MM-DD*,  where file_dir is the complete file path to the data
                files, and YYYY-MM-DD is the four digit year, two digit month, and two
                digit day in the desired output file names (e.g. 2015-03-15 for
                March 15, 2015).
            variables_requested = a list of variable name strings chosen from the
                model_varnames dictionary in this script, specifically the first
                item in the list associated with a given key.
                - If empty, the reader functionalizes all possible variables (default)
                - If 'all', the reader returns the model_varnames dictionary above
                    for only the variables present in the given files. Note: the
                    fulltime keyword must be False to acheive this behavior.
            filetime = boolean (default = False)
                - if False, the script fully executes.
                - If True, the script only executes far enough to determine the
                    time values associated with the chosen data.
                Note: The behavior of the script is determined jointly by the
                    filetime and fulltime keyword values.
            printfiles = boolean (default = False)
                - If False, the filenames associated with the data retrieved ARE
                    NOT printed.
                - If True, the filenames associated with the data retrieved ARE
                    printed.
            gridded_int = boolean (default = True)
                - If True, the variables chosen are functionalized in both the
                    standard method and a gridded method.
                - If False, the variables chosen are functionalized in only the
                    standard method.
            fulltime = boolean (default = True)
                - If True, linear interpolation in time between files is included
                    in the returned interpolator functions.
                - If False, no linear interpolation in time between files is included.
            verbose = boolean (False)
                - If False, script execution and the underlying Kamodo execution
                    is quiet except for specified messages.
                - If True, be prepared for a plethora of messages.
        All inputs are described in further detail in KamodoOnboardingInstructions.pdf.

        Returns: a kamodo object (see Kamodo core documentation) containing all
            requested variables in functionalized form.
        '''
        def __init__(self, full_file_prefix, variables_requested=[],
                     filetime=False, printfiles=False, gridded_int=True,
                     fulltime=True, verbose=False, **kwargs):

            # only the density, height and neutral files are combined
            super(MODEL, self).__init__()
            self.modelname = 'CTIPe'

            # check for prepared file of given prefix
            t0 = perf_counter()
            file_prefix = basename(full_file_prefix)  # YYYY-MM-DD
            file_dir = full_file_prefix.split(file_prefix)[0]
            if isfile(full_file_prefix+'.nc'):   # file already prepared!
                cdf_file = full_file_prefix+'.nc'  # input file name
                self.conversion_test = True
            else:   # file not prepared,  prepare it
                from kamodo_ccmc.readers.ctipe_tocdf import ctipe_combine_files
                cdf_file = ctipe_combine_files(full_file_prefix)
                self.conversion_test = True

            # establish time attributes first
            cdf_data = Dataset(cdf_file, 'r')
            self.filedate = datetime.strptime(file_prefix+' 00:00:00',
                                              '%Y-%m-%d %H:%M:%S').replace(
                                                  tzinfo=timezone.utc)
            t = array(cdf_data.variables['time'])
            self.datetimes = [datetime.utcfromtimestamp(t[0]).isoformat(sep=' '),
                              datetime.utcfromtimestamp(t[-1]).isoformat(sep=' ')]
            self.filetimes = [t[0], t[-1]]   # timestamps in hours
            self.dt = diff(t).max()  # t is in seconds

            # execute logic for finding nearest time in neighboring file
            # (used when searching for neighboring files below)
            if filetime and not fulltime:
                cdf_data.close()
                return  # return times as is to prevent recursion

            # if variables are given as integers,  convert to standard names
            if len(variables_requested) > 0:
                if isinstance(variables_requested[0], int):
                    print('Integers detected. Converting.')
                    tmp_var = [value[0] for key, value in model_varnames.items()
                               if value[2] in variables_requested]
                    variables_requested = tmp_var
                    print('Converted: ', variables_requested)

            if fulltime:   # add boundary time (default value)
                # find other files with same pattern
                from glob import glob

                file_pattern = file_dir + '*.nc'  # returns a string
                files = sorted(glob(file_pattern))
                prefix_list = unique([basename(f)[:10] for f in files
                                      if 'CTIPe' not in basename(f)])

                # find closest file by utc timestamp
                # ctipe has an open time at the beginning,
                # so need an earlier time from the closest file
                # files are automatically sorted by YYMMDD,
                # so next file is next in the list
                current_idx = where(prefix_list == file_prefix)[0]
                if current_idx == 0:
                    if verbose:
                        print('No earlier file available.')
                    filecheck = False
                    if filetime:
                        cdf_data.close()
                        return
                else:
                    # -1 for adding an earlier time
                    min_file_prefix = prefix_list[current_idx-1][0]
                    kamodo_test = MODEL(file_dir+min_file_prefix, filetime=True,
                                        fulltime=False)
                    if not kamodo_test.conversion_test:
                        if verbose:
                            print('No earlier file available.')
                        filecheck = False
                        if filetime:
                            cdf_data.close()
                            return
                    else:
                        time_test = abs(kamodo_test.filetimes[1] - self.filetimes[0])
                        # if nearest file time at least within one timestep (s)
                        if time_test <= self.dt:
                            filecheck = True
                            self.datetimes[0] = kamodo_test.datetimes[1]
                            self.filetimes[0] = kamodo_test.filetimes[1]

                            # time only version if returning time for searching
                            if filetime:
                                cdf_data.close()
                                return
                            # return object with additional time (for SF code)

                            # get kamodo object with same requested variables
                            # to add to each array below
                            if verbose:
                                print(f'Took {perf_counter()-t0:.3f}s to find closest file.')
                            kamodo_neighbor = MODEL(file_dir+min_file_prefix,
                                                    variables_requested=variables_requested,
                                                    fulltime=False)
                            short_data = kamodo_neighbor.short_data
                            if verbose:
                                print(f'Took {perf_counter()-t0:.3f}s to get ' +
                                      'data from closest file.')
                        else:
                            if verbose:
                                print(f'{file_prefix} No earlier file found ' +
                                      'within {diff(t).max():.1f}s.')
                            filecheck = False
                            if filetime:
                                cdf_data.close()
                                return

            # collect variable dependency lists
            self.ilev1_list = [value[0] for key, value in model_varnames.items()
                               if value[5][-1] == 'ilev1']
            self.ilev_list = [value[0] for key, value in model_varnames.items()
                              if value[5][-1] == 'ilev']
            # self.Elat_list = [value[0] for key,  value in model_varnames.items()
            #                  if value[5][-1]=='Elat']

            # perform initial check on variables_requested list
            if len(variables_requested) > 0 and fulltime and variables_requested != 'all':
                test_list = [value[0] for key, value in model_varnames.items()]
                err_list = [item for item in variables_requested if item not in
                            test_list]
                if len(err_list) > 0:
                    print('Variable name(s) not recognized: ', err_list)

            # collect variable list
            if len(variables_requested) > 0 and variables_requested != 'all':
                gvar_list = [key for key, value in model_varnames.items()
                             if value[0] in variables_requested and key in
                             cdf_data.variables.keys()]  # file variable names

                # check for variables requested but not available
                if len(gvar_list) != len(variables_requested):
                    err_list = [value[0] for key, value in model_varnames.items()
                                if value[0] in variables_requested and
                                key not in cdf_data.variables.keys()]
                    if len(err_list) > 0:
                        print('Some requested variables are not available: ',
                              err_list)

                # check that the appropriate height variable is added if needed
                check_list = [key for key, value in model_varnames.items()
                              if value[0] in self.ilev_list and key in gvar_list]
                # force addition of H for conversion of ilev to H and back
                if 'H_ilev' not in gvar_list and len(check_list) > 0:
                    gvar_list.append('H_ilev')
                check_list = [key for key, value in model_varnames.items()
                              if value[0] in self.ilev1_list and key in gvar_list]
                if 'H_ilev1' not in gvar_list and len(check_list) > 0:
                    gvar_list.append('H_lev')
                # check_list = [key for key,  value in model_varnames.items()\
                #              if value[0] in self.Elat_list and key in gvar_list]
                # if 'ZMAG' not in gvar_list and len(check_list)>0:
                #    gvar_list.append('ZMAG')
            else:  # only input variables on the avoid_list if specifically requested
                avoid_list = []
                gvar_list = [key for key in cdf_data.variables.keys()
                             if key in model_varnames.keys() and key not in
                             avoid_list]
                if not fulltime and variables_requested == 'all':
                    self.var_dict = {value[0]: value[1:] for key, value in
                                     model_varnames.items() if key in gvar_list}
                    cdf_data.close()
                    return

            # Store the requested variables into a dictionary
            variables = {model_varnames[key][0]: {'units': model_varnames[key][-1],
                                                  'data': array(cdf_data.variables[key])}
                         for key in gvar_list}  # store with key = standardized name

            # prepare and return data only for last timestamp
            if not fulltime:
                cdf_data.close()
                variables['time'] = self.filetimes[1]
                self.short_data = variables
                return

            # print files to screen if option requested
            self.filename = cdf_data.file.replace(',', '\n')
            if printfiles:
                print('Files: \n',  self.filename)
            self.missing_value = NaN
            self._registered = 0

            # Store coordinate data as class attributes
            if filecheck:   # don't need time correction b/c in utctimestamps
                # new time in hours since midnight
                new_time = ts_to_hrs(short_data['time'], self.filedate)
                # convert to hours since midnight
                self._time = insert(ts_to_hrs(t, self.filedate), 0, new_time)
            else:
                self._time = ts_to_hrs(t, self.filedate)

            # store dimensions if requested variable(s) require it
            if 'ilev' in cdf_data.variables.keys():  # if neutral file existed
                self._ilev = array(cdf_data.variables['ilev'])  # neutral file
                self._lon0 = array(cdf_data.variables['lon_n'])
                self._lat0 = array(cdf_data.variables['lat_n'])
                self._Elat = array(cdf_data.variables['Elat'])
                self._Elon = array(cdf_data.variables['Elon'])
            if 'lev' in cdf_data.variables.keys():  # if density file existed
                self._ilev1 = array(cdf_data.variables['lev'])      # density file
                self._lon1 = array(cdf_data.variables['lon_d'])
                self._lat1 = array(cdf_data.variables['lat_d'])
            if 'height' in cdf_data.variables.keys():  # if height file existed
                self._height = array(cdf_data.variables['height'])  # height file
                self._lon = array(cdf_data.variables['lon_h'])
                self._lat = array(cdf_data.variables['lat_h'])
            cdf_data.close()

            # register interpolators for each requested variable
            # store original list b/c gridded interpolators
            varname_list, self.variables = [key for key in variables.keys()], {}
            t_reg = perf_counter()
            for varname in varname_list:
                if len(variables[varname]['data'].shape) == 3:
                    if filecheck:   # if neighbor found
                        # append data for first time stamp, transpose and register
                        data_shape = list(variables[varname]['data'].shape)
                        data_shape[0] += 1  # add space for time
                        new_data = zeros(data_shape)
                        # put in current data
                        new_data[1:, :, :] = variables[varname]['data']
                        # add in data for additional time
                        new_data[0, :, :] = short_data[varname]['data'][-1, :, :]
                    else:
                        new_data = variables[varname]['data']
                    self.variables[varname] = dict(units=variables[varname]['units'],
                                                   data=new_data)
                    self.register_3D_variable(self.variables[varname]['units'],
                                              self.variables[varname]['data'],
                                              varname, gridded_int)
                elif len(variables[varname]['data'].shape) == 4:
                    if filecheck:
                        # append data for first time stamp,  transpose and register
                        data_shape = list(variables[varname]['data'].shape)
                        data_shape[0] += 1  # add space for time
                        new_data = zeros(data_shape)
                        # put in current data
                        new_data[1:, :, :, :] = variables[varname]['data']
                        # add in data for additional time
                        new_data[0, :, :, :] = short_data[varname]['data'][-1, :, :, :]
                    else:
                        new_data = variables[varname]['data']
                    self.variables[varname] = dict(units=variables[varname]['units'],
                                                   data=new_data)
                    self.register_4D_variable(self.variables[varname]['units'],
                                              self.variables[varname]['data'],
                                              varname, gridded_int)
            if verbose:
                print(f'Took {perf_counter() - t_reg:.5f}s to register ' +
                      f'{len(varname_list)} variables.')
            if verbose:
                print(f'Took a total of {perf_counter() - t0:.5f}s to kamodofy ' +
                      f'{len(varname_list)} variables.')
            return

        # define and register a 3D variable
        def register_3D_variable(self, units, variable, varname, gridded_int):
            """Registers a 3d interpolator with 3d signature"""

            # determine coordinate variables and xvec by coord list
            coord_list = [value[5][-1] for key, value in model_varnames.items()
                          if value[0] == varname][0]
            if 'lat' in coord_list:   # 3D variables only come from the neutral file
                lon,  lat = self._lon0, self._lat0
                xvec_dependencies = {'time': 'hr', 'lon': 'deg', 'lat': 'deg'}
            if 'Elat' in coord_list:
                lon,  lat = self._Elon, self._Elat
                xvec_dependencies = {'time': 'hr', 'Elon': 'deg', 'Elat': 'deg'}

            # define and register the interpolators
            self = regdef_3D_interpolators(self, units, variable, self._time,
                                           lon, lat, varname, xvec_dependencies,
                                           gridded_int)
            return

        # define and register a 4D variable
        def register_4D_variable(self, units, variable, varname, gridded_int):
            """Registers a 4d interpolator with 4d signature"""

            # determine coordinate variables by coord list
            coord_list = [value[5][-1] for key, value in model_varnames.items()
                          if value[0] == varname][0]
            if 'ilev1' in coord_list and hasattr(self, '_ilev1'):
                lon, lat, z = self._lon1, self._lat1, self._ilev1
                xvec_dependencies = {'time': 'hr', 'lon': 'deg', 'lat': 'deg',
                                     'ilev1': 'm/m'}
            if 'ilev' in coord_list and hasattr(self, '_ilev'):
                lon, lat, z = self._lon0, self._lat0, self._ilev
                xvec_dependencies = {'time': 'hr', 'lon': 'deg', 'lat': 'deg',
                                     'ilev': 'm/m'}
            if varname == 'N_e':  # special case b/c in more than one file
                if not hasattr(self, '_height'):  # height file DNE
                    lon, lat, z = self._lon0, self._lat0, self._ilev
                    xvec_dependencies = {'time': 'hr', 'lon': 'deg', 'lat': 'deg',
                                         'ilev': 'm/m'}
                else:
                    lon, lat, z = self._lon, self._lat, self._height
                    xvec_dependencies = {'time': 'hr', 'lon': 'deg', 'lat': 'deg',
                                         'height': 'km'}
            if 'height' in coord_list and varname != 'N_e' and hasattr(self, '_height'):
                lon, lat, z = self._lon, self._lat, self._height
                xvec_dependencies = {'time': 'hr', 'lon': 'deg', 'lat': 'deg',
                                 'height': 'km'}

            # define and register the interpolators
            self = regdef_4D_interpolators(self, units, variable, self._time,
                                           lon, lat, z, varname,
                                           xvec_dependencies, gridded_int)
            return

    return MODEL
