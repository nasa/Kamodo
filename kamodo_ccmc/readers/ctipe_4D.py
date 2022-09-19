'''
Spatial-only interpolation form:  Lutz Raestatter
Modified from spatial-only interpolation into this form:  Rebecca Ringuette
Kamodofication of the CTIPe model output

To check for pep8 standards:
    - pip install pycodestyle
    - pycodestyle ctipe_4D.py  # to show all errors
    - pycodestyle --first ctipe_4D.py  # to show first instances of error type
    - pycodestyle --show-source ctipe_4D.py  # to show full text per error
https://pypi.org/project/pycodestyle/ for more info.
'''
from datetime import datetime,  timezone
from numpy import vectorize

# constants and dictionaries
model_varnames = {'rho': ['rho_ilev1', 'total mass density', 0, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'ilev1'], 'kg/m**3'],
                  'rho_2': ['rho', 'total mass density', 0, 'GDZ', 'sph',
                            ['time', 'lon', 'lat', 'height'], 'kg/m**3'],
                  'T': ['T_ilev1', 'temperature', 1, 'GDZ', 'sph',
                        ['time', 'lon', 'lat', 'ilev1'], 'K'],
                  'T_2': ['T', 'temperature', 1, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'height'], 'K'],
                  'T_e': ['T_e', 'electron temperature', 2, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'height'], 'K'],
                  'T_i': ['T_i', 'ion temperature', 3, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'height'], 'K'],
                  'H_ilev': ['H_ilev',
                             'height dependent on primary pressure level',
                             4, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                             'm'],
                  'H_lev': ['H_ilev1',
                            'height dependent on secondary pressure level',
                            5, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev1'],
                            'm'],
                  'Vn_lat': ['v_nnorth_ilev',
                             'meridional neutral wind velocity (north)',
                             6, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                             'm/s'],
                  'Vn_lat_2': ['v_nnorth',
                               'meridional neutral wind velocity (north)',
                               6, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                 'height'], 'm/s'],
                  'Vn_lon': ['v_neast_ilev',
                             'zonal neutral wind velocity (east)',
                             7, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                             'm/s'],
                  'Vn_lon_2': ['v_neast', 'zonal neutral wind velocity (east)',
                               7, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                 'height'], 'm/s'],
                  'Vn_H': ['v_nup_ilev', 'vertical neutral wind velocity (up)',
                           8, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                           'm/s'],
                  'Vn_H_2': ['v_nup', 'vertical neutral wind velocity (up)',
                             8, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                             'm/s'],
                  'T_n': ['T_n_ilev', 'neutral temperature', 9, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'ilev'], 'K'],
                  'T_n_2': ['T_n', 'neutral temperature', 9, 'GDZ', 'sph',
                            ['time', 'lon', 'lat', 'height'], 'K'],
                  'Rmt': ['m_avgmol_ilev1', 'mean molecular mass', 10, 'GDZ',
                          'sph', ['time', 'lon', 'lat', 'ilev1'], 'amu'],
                  'Rmt_2': ['m_avgmol', 'mean molecular mass', 10, 'GDZ',
                            'sph', ['time', 'lon', 'lat', 'height'], 'amu'],
                  'N_e': ['N_e', 'electron number density', 11, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'height'], '1/m**3'],
                  # 'N_n': ['N_n', 'variable description',
                  # 12, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],'1/m**3'],
                  'Q_Solar': ['Q_Solar_ilev', 'solar heating', 13, 'GDZ',
                              'sph', ['time', 'lon', 'lat', 'ilev'], 'J/kg/s'],
                  'Q_Solar_2': ['Q_Solar', 'solar heating', 13, 'GDZ', 'sph',
                                ['time', 'lon', 'lat', 'height'], 'J/kg/s'],
                  'Q_Joule': ['Q_Joule_ilev', 'joule heating', 14, 'GDZ',
                              'sph', ['time', 'lon', 'lat', 'ilev'], 'J/kg/s'],
                  'Q_Joule_2': ['Q_Joule', 'joule heating', 14, 'GDZ', 'sph',
                                ['time', 'lon', 'lat', 'height'], 'J/kg/s'],
                  'Q_radiation': ['Q_rad_ilev', 'radiative heating or cooling',
                                  15, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                     'ilev'], 'J/kg/s'],
                  'Q_radiation_2': ['Q_rad', 'radiative heating or cooling',
                                    15, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                       'height'], 'J/kg/s'],
                  'N_O': ['N_O_ilev', 'number density of atomic oxygen', 16,
                          'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                          '1/m**3'],
                  'N_O_2': ['N_O', 'number density of atomic oxygen', 16,
                            'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                            '1/m**3'],
                  'N_O2': ['N_O2_ilev', 'number density of molecular oxygen',
                           17, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                           '1/m**3'],
                  'N_O2_2': ['N_O2', 'number density of molecular oxygen',
                             17, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                'height'], '1/m**3'],
                  'N_N2': ['N_N2_ilev', 'number density of molecular nitrogen',
                           18, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                           '1/m**3'],
                  'N_N2_2': ['N_N2', 'number density of molecular nitrogen',
                             18, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                'height'], '1/m**3'],
                  'N_NO': ['N_NO_ilev',
                           'number density of molecular nitric oxide', 19,
                           'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                           '1/m**3'],
                  'N_NO_2': ['N_NO',
                             'number density of molecular nitric oxide', 19,
                             'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                             '1/m**3'],
                  'N_NOplus': ['N_NOplus_ilev',
                               'number density of nitric oxide ion', 20, 'GDZ',
                               'sph', ['time', 'lon', 'lat', 'ilev'],
                               '1/m**3'],
                  'N_NOplus_2': ['N_NOplus',
                                 'number density of nitric oxide ion', 20,
                                 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                'height'], '1/m**3'],
                  'N_N2plus': ['N_N2plus_ilev',
                               'number density of molecular nitrogen ion', 21,
                               'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                               '1/m**3'],
                  'N_N2plus_2': ['N_N2plus',
                                 'number density of molecular nitrogen ion',
                                 21, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                    'height'], '1/m**3'],
                  'N_O2plus': ['N_O2plus_ilev',
                               'number density of molecular oxygen ion', 22,
                               'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                               '1/m**3'],
                  'N_O2plus_2': ['N_O2plus',
                                 'number density of molecular oxygen ion', 22,
                                 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                'height'], '1/m**3'],
                  'N_Nplus': ['N_Nplus_ilev',
                              'number density of atomic nitrogen ion', 23,
                              'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                              '1/m**3'],
                  'N_Nplus_2': ['N_Nplus',
                                'number density of atomic nitrogen ion', 23,
                                'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                                '1/m**3'],
                  'N_Oplus': ['N_Oplus', 'number density of atomic oxygen ion',
                              24, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                 'height'], '1/m**3'],
                  'N_Hplus': ['N_Hplus',
                              'number density of atomic hydrogen ion', 25,
                              'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                              '1/m**3'],
                  'Sigma_P': ['sigma_P_ilev', 'Pedersen conductivity', 26,
                              'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                              'S/m'],
                  'Sigma_P_2': ['sigma_P', 'Pedersen conductivity', 26, 'GDZ',
                                'sph', ['time', 'lon', 'lat', 'height'],
                                'S/m'],
                  'Sigma_H': ['sigma_H_ilev', 'Hall conductivity', 27, 'GDZ',
                              'sph', ['time', 'lon', 'lat', 'ilev'], 'S/m'],
                  'Sigma_H_2': ['sigma_H', 'Hall conductivity', 27, 'GDZ',
                                'sph', ['time', 'lon', 'lat', 'height'],
                                'S/m'],
                  'Vi_lon': ['v_inorth_ilev',
                             'meridional ion wind velocity (north)', 28, 'GDZ',
                             'sph', ['time', 'lon', 'lat', 'ilev'], 'm/s'],
                  'Vi_lon_2': ['v_inorth',
                               'meridional ion wind velocity (north)', 28,
                               'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                               'm/s'],
                  'Vi_lat': ['v_ieast_ilev', 'zonal ion wind velocity (east)',
                             29, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                             'm/s'],
                  'Vi_lat_2': ['v_ieast', 'zonal ion wind velocity (east)',
                               29, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                  'height'], 'm/s'],
                  'W_Joule': ['W_JouleH', 'height integrated joule heating',
                              30, 'GDZ', 'sph', ['time', 'lon', 'lat'],
                              'W/m**2'],
                  'Eflux_precip': ['Phi_E', 'energy flux', 31, 'GDZ', 'sph',
                                   ['time', 'lon', 'lat'], 'mW/m**2'],
                  'Eavg_precip': ['E_avg', 'average energy', 32, 'GDZ', 'sph',
                                  ['time', 'lon', 'lat'], 'keV'],
                  'TEC': ['TEC', 'vertical total electron content (height ' +
                          'integrated from bottom to top boundary)', 33, 'GDZ',
                          'sph', ['time', 'lon', 'lat'],
                          '1/m**2'],  # '10**16/m**2'
                  'E_theta140km': ['E_theta140km', 'Electric field at 140 km' +
                                   ', theta component', 34, 'GDZ', 'sph',
                                   ['time', 'Elon', 'Elat'], 'V/m'],
                  'E_lambda140km': ['E_lambda140km', 'Electric field at 140 ' +
                                    'km, lambda component', 35, 'GDZ', 'sph',
                                    ['time', 'Elon', 'Elat'], 'V/m'],
                  'E_theta300km': ['E_theta300km', 'Electric field at 300 ' +
                                   'km, theta component', 36, 'GDZ', 'sph',
                                   ['time', 'Elon', 'Elat'], 'V/m'],
                  'E_lambda300km': ['E_lambda300km', 'Electric field at 300 ' +
                                    'km, lambda component', 37, 'GDZ', 'sph',
                                    ['time', 'Elon', 'Elat'], 'V/m']}


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
    import kamodo_ccmc.readers.reader_utilities as RU

    # main class
    class MODEL(Kamodo):
        '''CTIPe model data reader.

        Inputs:
            full_file_prefix:  a string representing the file pattern of the
                model output data.
                Note:  This reader takes a file pattern of the format
                file_dir+YYYY-MM-DD*,  where file_dir is the complete file path
                to the data files, and YYYY-MM-DD is the four digit year,
                two digit month, and two digit day in the desired output file
                names (e.g. 2015-03-15 for March 15, 2015).
            variables_requested = a list of variable name strings chosen from
                the model_varnames dictionary in this script, specifically the
                firstcitem in the list associated with a given key.
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
            self.datetimes = \
                [datetime.utcfromtimestamp(t[0]).isoformat(sep=' '),
                 datetime.utcfromtimestamp(t[-1]).isoformat(sep=' ')]
            self.filetimes = [t[0], t[-1]]   # timestamps in hours
            self.dt = diff(t).max()  # t is in seconds

            # execute logic for finding nearest time in neighboring file
            # (used when searching for neighboring files below)
            if filetime and not fulltime:
                cdf_data.close()
                return  # return times as is to prevent recursion

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
                    kamodo_test = MODEL(file_dir+min_file_prefix,
                                        filetime=True, fulltime=False)
                    if not kamodo_test.conversion_test:
                        if verbose:
                            print('No earlier file available.')
                        filecheck = False
                        if filetime:
                            cdf_data.close()
                            return
                    else:
                        time_test = abs(kamodo_test.filetimes[1] -
                                        self.filetimes[0])
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
                                print(f'Took {perf_counter()-t0:.3f}s to ' +
                                      'find the closest file.')
                            kamodo_neighbor = \
                                MODEL(file_dir+min_file_prefix,
                                      variables_requested=variables_requested,
                                      fulltime=False)
                            short_data = kamodo_neighbor.short_data
                            if verbose:
                                print(f'Took {perf_counter()-t0:.3f}s to get' +
                                      ' data from closest file.')
                        else:
                            if verbose:
                                print(f'{file_prefix} No earlier file found ' +
                                      'within {diff(t).max():.1f}s.')
                            filecheck = False
                            if filetime:
                                cdf_data.close()
                                return

            # collect variable dependency lists
            ilev1_list = [value[0] for key, value in model_varnames.items()
                          if value[5][-1] == 'ilev1']
            ilev1_replace = [item.split('_ilev1')[0] for item in ilev1_list if
                             item != 'H_ilev1']
            ilev_list = [value[0] for key, value in model_varnames.items() if
                         value[5][-1] == 'ilev']
            ilev_replace = [item.split('_ilev')[0] for item in ilev_list if
                            item != 'H_ilev']
            self.total_ilev = [item for item in ilev_list + ilev1_list if item
                               not in ['H_ilev', 'H_ilev1']]  # total ilev list
            self.total_replace = ilev_replace + ilev1_replace
            # dictionary mapping to navigate related variable names
            self.ilev_map = {item1: item2 for item1, item2 in
                             zip(self.total_replace, self.total_ilev)}

            # perform initial check on variables_requested list
            if (len(variables_requested) > 0) and (variables_requested !=
                                                   'all') and fulltime:
                test_list = [value[0] for key, value in
                             model_varnames.items()]
                err_list = [item for item in variables_requested if item
                            not in test_list]
                if len(err_list) > 0:
                    print('Variable name(s) not recognized: ', err_list)

            # collect variable list
            if len(variables_requested) > 0 and variables_requested != 'all':
                # add ilev version of variables to the list, adding H_ilev(1)
                add_ilev = [var+'_ilev' for var in variables_requested if var
                            in ilev_replace]
                add_ilev1 = [var+'_ilev1' for var in variables_requested if var
                             in ilev1_replace]
                if len(add_ilev) > 0:
                    add_ilev += ['H_ilev']
                if len(add_ilev1) > 0:
                    add_ilev1 += ['H_ilev1']
                new_var = variables_requested + add_ilev + add_ilev1
                short_var = [item for item in new_var if item not
                             in self.ilev_map.keys()]  # remove replaced items
                gvar_list = [key for key, value in model_varnames.items()
                             if value[0] in short_var and key in
                             cdf_data.variables.keys()]  # file variable names

                # check for variables requested but not available
                if len(gvar_list) != len(short_var):
                    err_list = [value[0] for key, value in
                                model_varnames.items() if value[0] in
                                short_var and key not in
                                cdf_data.variables.keys()]
                    if len(err_list) > 0:
                        print('Some requested variables are not available: ',
                              err_list)
                    if 'H_ilev' in err_list or 'H_ilev1' in err_list:
                        other_name = ['H_ilev', 'H_ilev1']
                        err_name = [name for name in other_name if name in
                                    err_list][0]
                        other_name.remove(err_name)  # see first element
                        gvar_name = [key for key, value in
                                     model_varnames.items() if value[0] ==
                                     other_name[0]]
                        gvar_list += gvar_name  # add other height name to list
                        print(f'Retrieving the {other_name[0]} variable' +
                              ' instead.')

            else:  # return full possible variable list
                gvar_list = [key for key in cdf_data.variables.keys()
                             if key in model_varnames.keys()]
                if not fulltime and variables_requested == 'all':
                    self.var_dict = {value[0]: value[1:] for key, value in
                                     model_varnames.items() if key in
                                     gvar_list}
                    # add non-ilev versions of the variables in the files
                    key_list = list(self.var_dict.keys())
                    for var_key in key_list:
                        if var_key in self.total_ilev:
                            # retrieve equivalent non-ilev variable name
                            new_key = [key for key, value in
                                       self.ilev_map.items() if value ==
                                       var_key][0]
                            # retrieve key for model_varnames mapping
                            model_key = [key for key, value in
                                         model_varnames.items() if
                                         value[0] == new_key][0]
                            # add to returned dictionary
                            self.var_dict[model_varnames[model_key][0]] =\
                                model_varnames[model_key][1:]

                    cdf_data.close()
                    return

            # Store the requested variables into a dictionary
            variables = {model_varnames[key][0]:
                         {'units': model_varnames[key][-1],
                          'data': array(cdf_data.variables[key])}
                         for key in gvar_list}  # key = standardized name

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
                self._ilev1 = array(cdf_data.variables['lev'])      # density
                self._lon1 = array(cdf_data.variables['lon_d'])
                self._lat1 = array(cdf_data.variables['lat_d'])
            if 'height' in cdf_data.variables.keys():  # if height file existed
                self._height = array(cdf_data.variables['height'])  # height
                self._lon = array(cdf_data.variables['lon_h'])
                self._lat = array(cdf_data.variables['lat_h'])
            cdf_data.close()

            # Check for presence of necessary height variables in varname_list.
            varname_list = [key for key in variables.keys()]
            ilev1_check = unique([True for item in varname_list if 'ilev1' ==
                                  item[-5:]])
            ilev_check = unique([True for item in varname_list if 'ilev' ==
                                 item[-4:]])
            if ilev1_check and 'H_ilev1' not in varname_list:
                self.ilev_sub = 'H_ilev1'
            elif ilev_check and 'H_ilev' not in varname_list:
                self.ilev_sub = 'H_ilev'
            else:
                self.ilev_sub = False

            # register interpolators for each requested variable
            # rearrange to deal with H_ilev and H_ilev1 first if there
            if 'H_ilev' in varname_list:
                varname_list.remove('H_ilev')
                varname_list = ['H_ilev'] + varname_list
            if 'H_ilev1' in varname_list:
                varname_list.remove('H_ilev1')
                varname_list = ['H_ilev1'] + varname_list
            self.variables = {}
            t_reg = perf_counter()
            for varname in varname_list:
                if len(variables[varname]['data'].shape) == 3:
                    if filecheck:   # if neighbor found
                        # append data for first time stamp and transpose
                        data_shape = list(variables[varname]['data'].shape)
                        data_shape[0] += 1  # add space for time
                        new_data = zeros(data_shape)
                        # put in current data
                        new_data[1:, :, :] = variables[varname]['data']
                        # add in data for additional time
                        new_data[0, :, :] =\
                            short_data[varname]['data'][-1, :, :]
                    else:
                        new_data = variables[varname]['data']
                    self.variables[varname] =\
                        dict(units=variables[varname]['units'], data=new_data)
                    self.register_3D_variable(self.variables[varname]['units'],
                                              self.variables[varname]['data'],
                                              varname, gridded_int)
                elif len(variables[varname]['data'].shape) == 4:
                    if filecheck:
                        # append data for first time stamp, transpose
                        data_shape = list(variables[varname]['data'].shape)
                        data_shape[0] += 1  # add space for time
                        new_data = zeros(data_shape)
                        # put in current data
                        new_data[1:, :, :, :] = variables[varname]['data']
                        # add in data for additional time
                        new_data[0, :, :, :] = \
                            short_data[varname]['data'][-1, :, :, :]
                    else:
                        new_data = variables[varname]['data']
                    self.variables[varname] = \
                        dict(units=variables[varname]['units'],
                             data=new_data)
                    self.register_4D_variable(self.variables[varname]['units'],
                                              self.variables[varname]['data'],
                                              varname, gridded_int)
            if verbose:
                print(f'Took {perf_counter() - t_reg:.5f}s to register ' +
                      f'{len(varname_list)} variables.')
            if verbose:
                print(f'Took a total of {perf_counter() - t0:.5f}s to ' +
                      f'kamodofy {len(varname_list)} variables.')
            return

        # define and register a 3D variable
        def register_3D_variable(self, units, variable, varname, gridded_int):
            """Registers a 3d interpolator with 3d signature"""

            # determine coordinate variables and xvec by coord list
            coord_list = [value[5] for key, value in model_varnames.items()
                          if value[0] == varname][0]
            if 'lat' in coord_list:   # 3D variables come from neutral file
                lon,  lat = self._lon0, self._lat0
                xvec_dependencies = {'time': 'hr', 'lon': 'deg', 'lat': 'deg'}
            if 'Elat' in coord_list:
                lon,  lat = self._Elon, self._Elat
                xvec_dependencies = {'time': 'hr', 'Elon': 'deg',
                                     'Elat': 'deg'}

            # define and register the interpolators
            coord_str = [value[3]+value[4] for key, value in
                         model_varnames.items() if value[0] == varname][0]+'3D'
            self = RU.regdef_3D_interpolators(self, units, variable,
                                              self._time, lon, lat, varname,
                                              xvec_dependencies, gridded_int,
                                              coord_str)
            return

        # define and register a 4D variable
        def register_4D_variable(self, units, variable, varname, gridded_int):
            """Registers a 4d interpolator with 4d signature"""

            # determine coordinate variables by coord list
            coord_list = [value[5] for key, value in model_varnames.items()
                          if value[0] == varname][0]
            # both pressure level grids are present, just not both H functions
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
                    xvec_dependencies = {'time': 'hr', 'lon': 'deg',
                                         'lat': 'deg', 'ilev': 'm/m'}
                else:
                    lon, lat, z = self._lon, self._lat, self._height
                    xvec_dependencies = {'time': 'hr', 'lon': 'deg',
                                         'lat': 'deg', 'height': 'km'}
            if 'height' in coord_list and hasattr(self, '_height') and (
                    varname != 'N_e'):
                lon, lat, z = self._lon, self._lat, self._height
                xvec_dependencies = {'time': 'hr', 'lon': 'deg',
                                     'lat': 'deg', 'height': 'km'}

            # define and register the interpolators
            coord_str = [value[3]+value[4] for key, value in
                         model_varnames.items() if value[0] == varname][0]+'4D'
            # need H functions to be gridded regardless of gridded_int value
            h_grid = True if varname in ['H_ilev', 'H_ilev1'] else gridded_int
            self = RU.regdef_4D_interpolators(self, units, variable,
                                              self._time, lon, lat, z,
                                              varname, xvec_dependencies,
                                              h_grid, coord_str)

            # perform substitution if needed
            if isinstance(self.ilev_sub, str) and varname == self.ilev_sub:
                other_name = ['H_ilev', 'H_ilev1']
                other_name.remove(varname)  # first element is the other name
                print(f'{other_name[0]} missing in data and is needed to ' +
                      'convert the requested variables to depend on height.' +
                      f' Using {varname} instead.')
                xvec_dependencies = {'time': 'hr', 'lon': 'deg', 'lat': 'deg',
                                     other_name[0][2:]: 'm/m'}
                self.variables[other_name[0]] = dict(units=units,
                                                     data=variable)
                # register the other variable by linking to this one
                self = RU.register_interpolator(self, other_name[0], varname,
                                                xvec_dependencies)

            # create pressure level -> km function once per ilev type
            if varname in ['H_ilev', 'H_ilev1'] or varname in self.total_ilev:
                if varname in ['H_ilev', 'H_ilev1']:  # create custom interp
                    new_varname = 'P'+coord_list[-1][1:]
                    # Import and call custom interpolator
                    from ctipe_ilevinterp import PLevelInterp
                    interpolator, interp_ijk, kms = PLevelInterp(
                        self, self._time, lon, lat, z, 'H_'+coord_list[-1])
                    setattr(self, '_kms_'+coord_list[-1], kms)
                    units = 'm/m'
                    # kms is a 1D array of the median height values in km
                else:  # define by function composition
                    new_varname = varname.split('_ilev')[0]
                    # substitute kms array if height was also substituted
                    if isinstance(self.ilev_sub, str):
                        other_name = ['ilev', 'ilev1']
                        other_name.remove(self.ilev_sub[2:])
                        kms = getattr(self, '_kms_'+other_name[0])
                        interpolator = varname+'(P'+other_name[0][1:]+')'
                    else:
                        kms = getattr(self, '_kms_'+coord_list[-1])
                        interpolator = varname+'(P'+coord_list[-1][1:]+')'

                # Register in kamodo object
                new_xvec_dependencies = {'time': 'hr', 'lon': 'deg',
                                         'lat': 'deg', 'height': 'km'}
                self.variables[new_varname] = dict(units=units)
                self = RU.register_interpolator(self, new_varname,
                                                interpolator,
                                                new_xvec_dependencies)
                if varname in self.total_ilev:  # different if H vs not
                    interp_ijk = self[new_varname]

                # Create 'gridified' interpolators in the kamodo_object
                if gridded_int:
                    fake_data = zeros((len(self._time), len(lon), len(lat),
                                       len(kms)))  # saves execution time
                    self.variables[new_varname+'_ijk'] = dict(units=units)
                    gridded_interpolator = RU.define_4d_gridded_interpolator(
                        units, fake_data, self._time, lon, lat, kms,
                        new_xvec_dependencies, interp_ijk)
                    self = RU.register_interpolator(
                        self, new_varname+'_ijk', gridded_interpolator,
                        new_xvec_dependencies)
            return

    return MODEL
