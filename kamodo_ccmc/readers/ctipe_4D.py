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

filling on, default _FillValue of 9.969209968386869e+36 used
'''

# constants and dictionaries
model_varnames = {'density': ['rho_ilev1', 'total mass density', 0, 'GDZ',
                              'sph', ['time', 'lon', 'lat', 'ilev1'],
                              'kg/m**3'],
                  'density_2': ['rho', 'total mass density', 0, 'GDZ', 'sph',
                                ['time', 'lon', 'lat', 'height'], 'kg/m**3'],
                  'temperature': ['T_ilev1', 'temperature', 1, 'GDZ', 'sph',
                                  ['time', 'lon', 'lat', 'ilev1'], 'K'],
                  'temperature_2': ['T', 'temperature', 1, 'GDZ', 'sph',
                                    ['time', 'lon', 'lat', 'height'], 'K'],
                  'electron_temperature': ['T_e', 'electron temperature', 2,
                                           'GDZ', 'sph', ['time', 'lon', 'lat',
                                                          'height'], 'K'],
                  'ion_temperature': ['T_i', 'ion temperature', 3, 'GDZ',
                                      'sph', ['time', 'lon', 'lat', 'height'],
                                      'K'],
                  'height_n': ['H_ilev', 'height dependent on primary ' +
                               'pressure level', 4, 'GDZ', 'sph',
                               ['time', 'lon', 'lat', 'ilev'], 'm'],
                  'height_d': ['H_ilev1', 'height dependent on secondary ' +
                               'pressure level', 5, 'GDZ', 'sph',
                               ['time', 'lon', 'lat', 'ilev1'], 'm'],
                  'meridional_neutral_wind': ['v_nnorth_ilev', 'meridional ' +
                                              'neutral wind velocity (north)',
                                              6, 'GDZ', 'sph',
                                              ['time', 'lon', 'lat', 'ilev'],
                                              'm/s'],
                  'meridional_neutral_wind_2': ['v_nnorth', 'meridional ' +
                                                'neutral wind velocity ' +
                                                '(north)', 6, 'GDZ', 'sph',
                                                ['time', 'lon', 'lat',
                                                 'height'], 'm/s'],
                  'zonal_neutral_wind': ['v_neast_ilev', 'zonal neutral wind' +
                                         ' velocity (east)', 7, 'GDZ', 'sph',
                                         ['time', 'lon', 'lat', 'ilev'],
                                         'm/s'],
                  'zonal_neutral_wind_2': ['v_neast', 'zonal neutral wind ' +
                                           'velocity (east)', 7, 'GDZ', 'sph',
                                           ['time', 'lon', 'lat', 'height'],
                                           'm/s'],
                  'vertical_neutral_wind': ['v_nup_ilev', 'vertical neutral ' +
                                            'wind velocity (up)', 8, 'GDZ',
                                            'sph', ['time', 'lon', 'lat',
                                                    'ilev'], 'm/s'],
                  'vertical_neutral_wind_2': ['v_nup', 'vertical neutral ' +
                                              'wind velocity (up)', 8, 'GDZ',
                                              'sph', ['time', 'lon', 'lat',
                                                      'height'], 'm/s'],
                  'neutral_temperature': ['T_n_ilev', 'neutral temperature', 9,
                                          'GDZ', 'sph', ['time', 'lon', 'lat',
                                                         'ilev'], 'K'],
                  'neutral_temperature_2': ['T_n', 'neutral temperature', 9,
                                            'GDZ', 'sph',
                                            ['time', 'lon', 'lat', 'height'],
                                            'K'],
                  'mean_molecular_mass': ['m_avgmol_ilev1', 'mean molecular ' +
                                          'mass', 10, 'GDZ', 'sph',
                                          ['time', 'lon', 'lat', 'ilev1'],
                                          'amu'],
                  'mean_molecular_mass_2': ['m_avgmol', 'mean molecular mass',
                                            10, 'GDZ', 'sph',
                                            ['time', 'lon', 'lat', 'height'],
                                            'amu'],
                  'electron_density': ['N_e', 'electron number density', 11,
                                       'GDZ', 'sph',
                                       ['time', 'lon', 'lat', 'height'],
                                       '1/m**3'],
                  # 'N_n': ['N_n', 'variable description',
                  # 12, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],'1/m**3'],
                  'solar_heating': ['Q_Solar_ilev', 'solar heating', 13, 'GDZ',
                                    'sph', ['time', 'lon', 'lat', 'ilev'],
                                    'J/kg/s'],
                  'solar_heating_2': ['Q_Solar', 'solar heating', 13, 'GDZ',
                                      'sph', ['time', 'lon', 'lat', 'height'],
                                      'J/kg/s'],
                  'joule_heating': ['Q_Joule_ilev', 'joule heating', 14, 'GDZ',
                                    'sph', ['time', 'lon', 'lat', 'ilev'],
                                    'J/kg/s'],
                  'joule_heating_2': ['Q_Joule', 'joule heating', 14, 'GDZ',
                                      'sph', ['time', 'lon', 'lat', 'height'],
                                      'J/kg/s'],
                  'radiation_heat_cool': ['Q_rad_ilev', 'radiative heating ' +
                                          'or cooling', 15, 'GDZ', 'sph',
                                          ['time', 'lon', 'lat', 'ilev'],
                                          'J/kg/s'],
                  'radiation_heat_cool_2': ['Q_rad', 'radiative heating or ' +
                                            'cooling', 15, 'GDZ', 'sph',
                                            ['time', 'lon', 'lat', 'height'],
                                            'J/kg/s'],
                  'atomic_oxygen_density': ['N_O_ilev', 'number density of ' +
                                            'atomic oxygen', 16, 'GDZ', 'sph',
                                            ['time', 'lon', 'lat', 'ilev'],
                                            '1/m**3'],
                  'atomic_oxygen_density_2': ['N_O', 'number density of ' +
                                              'atomic oxygen', 16, 'GDZ',
                                              'sph', ['time', 'lon', 'lat',
                                                      'height'], '1/m**3'],
                  'molecular_oxygen_density': ['N_O2_ilev', 'number density ' +
                                               'of molecular oxygen', 17,
                                               'GDZ', 'sph',
                                               ['time', 'lon', 'lat', 'ilev'],
                                               '1/m**3'],
                  'molecular_oxygen_density_2': ['N_O2', 'number density of ' +
                                                 'molecular oxygen', 17, 'GDZ',
                                                 'sph', ['time', 'lon', 'lat',
                                                         'height'], '1/m**3'],
                  'molecular_nitrogen_density': ['N_N2_ilev', 'number ' +
                                                 'density of molecular ' +
                                                 'nitrogen', 18, 'GDZ', 'sph',
                                                 ['time', 'lon', 'lat',
                                                  'ilev'], '1/m**3'],
                  'molecular_nitrogen_density_2': ['N_N2', 'number density ' +
                                                   'of molecular nitrogen',
                                                   18, 'GDZ', 'sph',
                                                   ['time', 'lon', 'lat',
                                                    'height'], '1/m**3'],
                  'nitric_oxide_density': ['N_NO_ilev', 'number density of ' +
                                           'molecular nitric oxide', 19, 'GDZ',
                                           'sph', ['time', 'lon', 'lat',
                                                   'ilev'], '1/m**3'],
                  'nitric_oxide_density_2': ['N_NO', 'number density of ' +
                                             'molecular nitric oxide', 19,
                                             'GDZ', 'sph', ['time', 'lon',
                                                            'lat', 'height'],
                                             '1/m**3'],
                  'nitric_oxide_ion_density': ['N_NOplus_ilev', 'number ' +
                                               'density of nitric oxide ion',
                                               20, 'GDZ', 'sph',
                                               ['time', 'lon', 'lat', 'ilev'],
                                               '1/m**3'],
                  'nitric_oxide_ion_density_2': ['N_NOplus', 'number density' +
                                                 ' of nitric oxide ion', 20,
                                                 'GDZ', 'sph',
                                                 ['time', 'lon', 'lat',
                                                  'height'], '1/m**3'],
                  'molecular_nitrogen_ion_density': ['N_N2plus_ilev',
                                                     'number density of ' +
                                                     'molecular nitrogen ion',
                                                     21, 'GDZ', 'sph',
                                                     ['time', 'lon', 'lat',
                                                      'ilev'], '1/m**3'],
                  'molecular_nitrogen_ion_density_2': ['N_N2plus',
                                                       'number density of ' +
                                                       'molecular nitrogen ' +
                                                       'ion', 21, 'GDZ', 'sph',
                                                       ['time', 'lon', 'lat',
                                                        'height'], '1/m**3'],
                  'molecular_oxygen_ion_density': ['N_O2plus_ilev', 'number ' +
                                                   'density of molecular ' +
                                                   'oxygen ion', 22, 'GDZ',
                                                   'sph', ['time', 'lon',
                                                           'lat', 'ilev'],
                                                   '1/m**3'],
                  'molecular_oxygen_ion_density_2': ['N_O2plus', 'number ' +
                                                     'density of molecular ' +
                                                     'oxygen ion', 22, 'GDZ',
                                                     'sph', ['time', 'lon',
                                                             'lat', 'height'],
                                                     '1/m**3'],
                  'atomic_nitrogen_ion_density': ['N_Nplus_ilev', 'number ' +
                                                  'density of atomic ' +
                                                  'nitrogen ion', 23, 'GDZ',
                                                  'sph', ['time', 'lon', 'lat',
                                                          'ilev'], '1/m**3'],
                  'atomic_nitrogen_ion_density_2': ['N_Nplus', 'number ' +
                                                    'density of atomic ' +
                                                    'nitrogen ion', 23, 'GDZ',
                                                    'sph', ['time', 'lon',
                                                            'lat', 'height'],
                                                    '1/m**3'],
                  'atomic_oxygen_ion_density': ['N_Oplus', 'number density ' +
                                                'of atomic oxygen ion', 24,
                                                'GDZ', 'sph',
                                                ['time', 'lon', 'lat',
                                                 'height'], '1/m**3'],
                  'atomic_hydrogen_ion_density': ['N_Hplus', 'number density' +
                                                  ' of atomic hydrogen ion',
                                                  25, 'GDZ', 'sph',
                                                  ['time', 'lon', 'lat',
                                                   'height'], '1/m**3'],
                  'pedersen_conductivity': ['sigma_P_ilev', 'Pedersen ' +
                                            'conductivity', 26, 'GDZ', 'sph',
                                            ['time', 'lon', 'lat', 'ilev'],
                                            'S/m'],
                  'pedersen_conductivity_2': ['sigma_P', 'Pedersen ' +
                                              'conductivity', 26, 'GDZ', 'sph',
                                              ['time', 'lon', 'lat', 'height'],
                                              'S/m'],
                  'hall_conductivity': ['sigma_H_ilev', 'Hall conductivity',
                                        27, 'GDZ', 'sph', ['time', 'lon',
                                                           'lat', 'ilev'],
                                        'S/m'],
                  'hall_conductivity_2': ['sigma_H', 'Hall conductivity', 27,
                                          'GDZ', 'sph', ['time', 'lon', 'lat',
                                                         'height'], 'S/m'],
                  'zonal_ion_velocity': ['v_inorth_ilev', 'meridional ion ' +
                                         'wind velocity (north)', 28, 'GDZ',
                                         'sph', ['time', 'lon', 'lat', 'ilev'],
                                         'm/s'],
                  'zonal_ion_velocity_2': ['v_inorth', 'meridional ion wind ' +
                                           'velocity (north)', 28, 'GDZ',
                                           'sph', ['time', 'lon', 'lat',
                                                   'height'], 'm/s'],
                  'meridional_ion_velocity': ['v_ieast_ilev', 'zonal ion ' +
                                              'wind velocity (east)', 29,
                                              'GDZ', 'sph',
                                              ['time', 'lon', 'lat', 'ilev'],
                                              'm/s'],
                  'meridional_ion_velocity_2': ['v_ieast', 'zonal ion wind ' +
                                                'velocity (east)', 29, 'GDZ',
                                                'sph', ['time', 'lon', 'lat',
                                                        'height'], 'm/s'],
                  'height_integrated_joule_heating': ['W_JouleH', 'height ' +
                                                      'integrated joule ' +
                                                      'heating', 30, 'GDZ',
                                                      'sph', ['time', 'lon',
                                                              'lat'],
                                                      'W/m**2'],
                  'energy_influx': ['Phi_E', 'energy flux', 31, 'GDZ', 'sph',
                                    ['time', 'lon', 'lat'], 'mW/m**2'],
                  'mean_energy': ['E_avg', 'average energy', 32, 'GDZ', 'sph',
                                  ['time', 'lon', 'lat'], 'keV'],
                  'total_electron_content': ['TEC', 'vertical total electron' +
                                             ' content (height integrated ' +
                                             'from bottom to top boundary)',
                                             33, 'GDZ', 'sph',
                                             ['time', 'lon', 'lat'], '1/m**2'],
                  # '10**16/m**2'
                  'theta_electric_field_at_140km': ['E_theta140km',
                                                    'Electric field at 140 ' +
                                                    'km, theta component', 34,
                                                    'GDZ', 'sph',
                                                    ['time', 'Elon', 'Elat'],
                                                    'V/m'],
                  'lambda_electric_field_at_140km': ['E_lambda140km',
                                                     'Electric field at 140 ' +
                                                     'km, lambda component',
                                                     35, 'GDZ', 'sph',
                                                     ['time', 'Elon', 'Elat'],
                                                     'V/m'],
                  'theta_electric_field_at_300km': ['E_theta300km',
                                                    'Electric field at 300 ' +
                                                    'km, theta component', 36,
                                                    'GDZ', 'sph',
                                                    ['time', 'Elon', 'Elat'],
                                                    'V/m'],
                  'lambda_electric_field_at_300km': ['E_lambda300km',
                                                     'Electric field at 300 ' +
                                                     'km, lambda component',
                                                     37, 'GDZ', 'sph',
                                                     ['time', 'Elon', 'Elat'],
                                                     'V/m']}


def MODEL():
    from numpy import array, NaN, unique, append, where, transpose
    from time import perf_counter
    from glob import glob
    from os.path import basename, isfile
    from datetime import datetime, timezone
    from kamodo import Kamodo
    from netCDF4 import Dataset
    import kamodo_ccmc.readers.reader_utilities as RU

    # main class
    class MODEL(Kamodo):
        '''CTIPe model data reader.

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
            - CTIPe output files are given in multiple netCDF files per day. No
              attempt is made to combine these files, but some pre-processing
              is needed to calculate the median height across the entire
              dataset (time, lon, lat) for a given pressure level value.
            - CTIPe data is given in several coordinate systems, two depending
              on pressure level - a unitless representation of height
              corresponding to a given atmospheric pressure. The preset values
              of pressure level in the two coordinate systems are not
              guaranteed to be identical, so they are inverted independently of
              the other unless only one is provided. In that case, the are
              assumed to be identical.
            - Pressure level inversion is performed in the reader_utilities
              script, specifically the PLevelInterp function. Two versions of
              all variables that depend on pressure level are created: the
              original and one dependent on height, which is created through
              Kamodo's function composition feature.
            - The files are small and contain multiple time steps per file, so
              interpolation method 2 is chosen. The standard SciPy interpolator
              is used.
        '''
        def __init__(self, file_dir, variables_requested=[],
                     filetime=False, printfiles=False, gridded_int=True,
                     verbose=False, **kwargs):

            # only the density, height and neutral files are combined
            super(MODEL, self).__init__()
            self.modelname = 'CTIPe'
            t0 = perf_counter()

            # first, check for file list, create if DNE
            list_file = file_dir + self.modelname + '_list.txt'
            time_file = file_dir + self.modelname + '_times.txt'
            self.times, self.pattern_files = {}, {}
            if not isfile(list_file) or not isfile(time_file):
                # collect filenames
                all_files = sorted(glob(file_dir+'*-plot-*.nc'))
                files = [f for f in all_files if 'plasma' not in f]
                patterns = unique([basename(f)[16:-3] for f in files])
                # e.g. 'density', 'neutral', 'height'
                self.filename = ''.join([f+',' for f in files])[:-1]
                # datetime object for midnight on date
                date = basename(files[0])[:10]  # 'YYYY-MM-DD'
                self.filedate = datetime.strptime(date+' 00:00:00',
                                                  '%Y-%m-%d %H:%M:%S').replace(
                                                      tzinfo=timezone.utc)
                filedate_ts = self.filedate.timestamp()

                # establish time attributes, using 3D time as default
                for p in patterns:
                    # get list of files to loop through later
                    pattern_files = sorted(glob(file_dir+'*'+p+'*.nc'))
                    self.pattern_files[p] = pattern_files
                    self.times[p] = {'start': [], 'end': [], 'all': []}

                    # loop through to get times
                    for f in range(len(pattern_files)):
                        cdf_data = Dataset(pattern_files[f])
                        tmp = array(cdf_data.variables['time'])  # utc tstamps
                        self.times[p]['start'].append(tmp[0])
                        self.times[p]['end'].append(tmp[-1])
                        self.times[p]['all'].extend(tmp)
                        cdf_data.close()
                    # convert timestamps to be hrs since midnight of 1st file
                    self.times[p]['start'] = (array(self.times[p]['start']) -
                                              filedate_ts)/3600.
                    self.times[p]['end'] = (array(self.times[p]['end']) -
                                            filedate_ts)/3600.
                    self.times[p]['all'] = (array(self.times[p]['all']) -
                                            filedate_ts)/3600.

                # perform post-processing to speed up pressure level conversion
                from kamodo_ccmc.readers.ctipe_tocdf import convert_all
                convert_all(file_dir, self.pattern_files, self.times)

                # create time list file if DNE
                RU.create_timelist(list_file, time_file, self.modelname,
                                   self.times, self.pattern_files,
                                   self.filedate)
            else:  # read in data and time grids from file list
                self.times, self.pattern_files, self.filedate, self.filename =\
                    RU.read_timelist(time_file, list_file)
            # return time info
            if filetime:
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
                                                   'all'):
                test_list = [value[0] for key, value in
                             model_varnames.items()]
                err_list = [item for item in variables_requested if item
                            not in test_list]
                if len(err_list) > 0:
                    print('Variable name(s) not recognized: ', err_list)
                for item in err_list:
                    variables_requested.remove(item)
                if len(variables_requested) == 0:
                    return

                # determine variable mapping
                # add ilev version of variables to the list, adding H_ilev(1)
                add_ilev = [var+'_ilev' for var in variables_requested if var
                            in ilev_replace]
                add_ilev1 = [var+'_ilev1' for var in variables_requested if var
                             in ilev1_replace]
                if len(add_ilev) > 0 or len(add_ilev1) > 0:  # add both
                    add_ilev += ['H_ilev']  # might need to replace one
                    add_ilev1 += ['H_ilev1']  # with the other
                new_var = variables_requested + add_ilev + add_ilev1
                short_var = [item for item in new_var if item not
                             in self.ilev_map.keys()]  # remove replaced items

            # collect variables per pattern type
            self.varfiles, self.gvarfiles = {}, {}
            self.err_list, self.var_dict = [], {}
            for p in self.pattern_files.keys():
                # check var_list for variables not possible in this file set
                cdf_data = Dataset(self.pattern_files[p][0], 'r')
                if len(variables_requested) > 0 and\
                        variables_requested != 'all':
                    gvar_list = [key for key, value in model_varnames.items()
                                 if value[0] in short_var and key in
                                 cdf_data.variables.keys()]  # file var names
                    # check for variables requested but not available in p file
                    if len(gvar_list) != len(short_var):
                        err_list = [value[0] for key, value in
                                    model_varnames.items() if value[0] in
                                    short_var and key not in
                                    cdf_data.variables.keys()]
                        # deal with height variable in two files (same name!)
                        if p == 'density':
                            if 'height' in cdf_data.variables.keys() and (
                                    'H_ilev1' in err_list):
                                gvar_list.append('height')
                                err_list.remove('H_ilev1')
                            if 'H_ilev' in err_list:  # remove neutral H
                                err_list.remove('H_ilev')
                        if p == 'neutral':
                            if 'height' in cdf_data.variables.keys() and (
                                    'H_ilev' in err_list):
                                gvar_list.append('height')
                                err_list.remove('H_ilev')
                            if 'H_ilev1' in err_list:  # remove density H
                                err_list.remove('H_ilev1')
                        if p == 'height':  # remove H from err list
                            if 'H_ilev' in err_list:
                                err_list.remove('H_ilev')
                            if 'H_ilev1' in err_list:
                                err_list.remove('H_ilev1')
                        if len(err_list) > 0:
                            self.err_list.extend(err_list)

                else:  # return full possible variable list
                    gvar_list = [key for key in cdf_data.variables.keys()
                                 if key in model_varnames.keys()]
                    # deal with height variable in two files (same name!)
                    if p in ['density', 'neutral'] and (
                            'height' in cdf_data.variables.keys()):
                        gvar_list.append('height')
                    if variables_requested == 'all':
                        var_dict = {value[0]: value[1:] for key, value in
                                    model_varnames.items() if key in gvar_list}
                        # deal with height again
                        if p == 'density' and 'height' in gvar_list:
                            var_dict['H_ilev1'] = \
                                model_varnames['height_d'][1:]
                        if p == 'neutral' and 'height' in gvar_list:
                            var_dict['H_ilev'] = \
                                model_varnames['height_n'][1:]
                        # add non-ilev versions of the variables in the files
                        key_list = list(var_dict.keys())
                        for var_key in key_list:
                            self.var_dict[var_key] = var_dict[var_key]
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
                # store which file these variables came from
                self.varfiles[p] = [model_varnames[key][0] for key in gvar_list
                                    if key in model_varnames.keys()]
                if 'height' in gvar_list and p == 'density':
                    self.varfiles[p].append('H_ilev1')
                elif 'height' in gvar_list and p == 'neutral':
                    self.varfiles[p].append('H_ilev')
                self.gvarfiles[p] = gvar_list

            # return var_dict (already created)
            if variables_requested == 'all':
                return
            # clean up error list and then take action
            var_list = []
            for p in self.varfiles.keys():
                var_list.extend(self.varfiles[p])
            err_list = [var for var in self.err_list if var not in var_list]
            if len(err_list) > 0:
                print('Some requested variables are not available: ',
                      err_list)
            # perform height substitution if missing
            if 'H_ilev1' in self.err_list:  # density H missing
                self.gvarfiles['neutral'].append('height')
                print('Retrieving the H_ilev variable instead.')
            if 'H_ilev' in self.err_list:  # neutral H missing
                self.gvarfiles['density'].append('height')
                print('Retrieving the H_ilev1 variable instead.')

            # get coordinates from first file of each type
            self.variables = {}
            for p in self.pattern_files.keys():
                cdf_data = Dataset(self.pattern_files[p][0], 'r')
                lon = array(cdf_data.variables['lon'])
                lon_le180 = list(where(lon <= 180)[0])  # 0 to 180
                lon_ge180 = list(where(lon >= 180)[0])
                tmp = append(lon, 360.) - 180.  # now -180. to +180.
                setattr(self, '_lon_'+p, tmp)
                setattr(self, '_lon_idx_'+p, lon_ge180+lon_le180)
                setattr(self, '_lat_'+p,
                        array(cdf_data.variables['lat']))
                if p == 'density':
                    setattr(self, '_ilev1',
                            array(cdf_data.variables['plev']))
                    if isfile(file_dir+'CTIPe_km.nc'):  # km_ilev1 from file
                        km_data = Dataset(file_dir+'CTIPe_km.nc')
                        setattr(self, '_km_ilev1',
                                array(km_data.variables['km_ilev1']))
                        setattr(self, '_km_ilev1_max', km_data.km_ilev1_max)
                        setattr(self, '_km_ilev1_min', km_data.km_ilev1_min)
                        km_data.close()
                elif p == 'neutral':
                    setattr(self, '_ilev',
                            array(cdf_data.variables['plev']))
                    if isfile(file_dir+'CTIPe_km.nc'):  # km_ilev from file
                        km_data = Dataset(file_dir+'CTIPe_km.nc')
                        setattr(self, '_km_ilev',
                                array(km_data.variables['km_ilev']))
                        setattr(self, '_km_ilev_max', km_data.km_ilev_max)
                        setattr(self, '_km_ilev_min', km_data.km_ilev_min)
                        km_data.close()
                    lon = array(cdf_data.variables['elon'])  # 0 to 360
                    lon_le180 = list(where(lon <= 180)[0])  # 0 to 180
                    lon_ge180 = list(where((lon >= 180) & (lon < 360.))[0])
                    setattr(self, '_Elon_'+p, lon-180.)
                    setattr(self, '_Elon_idx_'+p, lon_ge180+lon_le180)
                    setattr(self, '_Elat_'+p,
                            array(cdf_data.variables['elat']))
                elif p == 'height':
                    setattr(self, '_height',
                            array(cdf_data.variables['ht']))
                cdf_data.close()

                # initialize variable dictionaries
                for var in self.gvarfiles[p]:
                    if var != 'height':
                        self.variables[model_varnames[var][0]] = {
                            'units': model_varnames[var][-1], 'data': p}
                    elif var == 'height' and p == 'density':
                        self.variables['H_ilev1'] = {
                            'units': model_varnames['height_d'][-1],
                            'data': p}
                    elif var == 'height' and p == 'neutral':
                        self.variables['H_ilev'] = {
                            'units': model_varnames['height_n'][-1],
                            'data': p}

            # print files to screen if option requested
            if printfiles:
                print('Files: \n',  self.filename)
            self.missing_value = NaN
            self._registered = 0

            # Check for presence of necessary height variables in varname_list.
            varname_list = [key for key in self.variables.keys()]
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
            t_reg = perf_counter()
            for varname in varname_list:
                self.register_variable(varname, gridded_int)
            if verbose:
                print(f'Took {perf_counter() - t_reg:.5f}s to register ' +
                      f'{len(varname_list)} variables.')
            if verbose:
                print(f'Took a total of {perf_counter() - t0:.5f}s to ' +
                      f'kamodofy {len(varname_list)} variables.')
            return

        # define and register a 3D variable
        def register_variable(self, varname, gridded_int):
            """Registers an interpolator with proper signature"""

            # determine coordinate variables and xvec by coord list
            key = self.variables[varname]['data']
            coord_list = [value[5] for key, value in model_varnames.items()
                          if value[0] == varname][0]
            coord_dict = {'time': {'units': 'hr',
                                   'data': self.times[key]['all']}}
            # special cases b/c in more than one file
            if (varname == 'N_e' and key != 'height') or \
                    (varname == 'm_avgmol_ilev1' and key != 'density'):
                key = 'neutral'
                coord_list = ['time', 'lon', 'lat', 'ilev']
            # retrieve coordinate grids
            if 'lat' in coord_list:
                coord_dict['lon'] = {'units': 'deg', 'data':
                                     getattr(self, '_lon_'+key)}
                coord_dict['lat'] = {'units': 'deg', 'data':
                                     getattr(self, '_lat_'+key)}
                lon_idx = getattr(self, '_lon_idx_'+key)
            if 'Elat' in coord_list:
                coord_dict['Elon'] = {'units': 'deg', 'data':
                                      getattr(self, '_Elon_'+key)}
                coord_dict['Elat'] = {'units': 'deg', 'data':
                                      getattr(self, '_Elat_'+key)}
                lon_idx = getattr(self, '_Elon_idx_'+key)
            if 'height' in coord_list:
                coord_dict['height'] = {'units': 'km', 'data':
                                        getattr(self, '_height')}
            if 'ilev' in coord_list:
                coord_dict['ilev'] = {'units': 'm/m', 'data':
                                      getattr(self, '_ilev')}
            if 'ilev1' in coord_list:
                coord_dict['ilev1'] = {'units': 'm/m', 'data':
                                       getattr(self, '_ilev1')}
            if varname not in ['H_ilev', 'H_ilev1']:
                gvar = [key for key, value in model_varnames.items() if
                        value[0] == varname][0]
            else:
                gvar = 'height'  # same in both files
            coord_str = [value[-4]+value[-3] for key, value in
                         model_varnames.items() if value[0] == varname][0]

            # define operations for each variable when given the key
            def func(i):
                '''key is the file pattern, start_idxs is a list of one or two
                indices matching the file start times in self.start_times[key].
                '''
                # get data from file
                file = self.pattern_files[key][i]
                cdf_data = Dataset(file)
                data = array(cdf_data.variables[gvar])
                if hasattr(cdf_data.variables[gvar][0], 'fill_value'):
                    fill_value = cdf_data.variables[gvar][0].fill_value
                else:
                    fill_value = None
                cdf_data.close()
                # if not the last file, tack on first time from next file
                if file != self.pattern_files[key][-1]:  # interp btwn files
                    next_file = self.pattern_files[key][i+1]
                    cdf_data = Dataset(next_file)
                    data_slice = array(cdf_data.variables[gvar][0])
                    cdf_data.close()
                    data = append(data, [data_slice], axis=0)
                # data wrangling
                if fill_value is not None:  # if defined, replace with NaN
                    data = where(data != fill_value, data, NaN)
                if len(data.shape) == 3 and 'Elon' not in coord_dict.keys():
                    variable = transpose(data, (0, 2, 1))
                elif 'Elon' in coord_dict.keys():
                    variable = data  # already correct dimension order
                elif len(data.shape) == 4:
                    variable = transpose(data, (0, 3, 2, 1))
                return variable[:, lon_idx]

            # define and register the interpolators, pull 3D data into arrays
            # need H functions to be gridded regardless of gridded_int value
            if varname in ['H_ilev', 'H_ilev1']:
                gridded_int = True
            self = RU.Functionalize_Dataset(
                self, coord_dict, varname, self.variables[varname],
                gridded_int, coord_str, interp_flag=2, func=func,
                times_dict=self.times[key])

            # perform substitution if needed
            if isinstance(self.ilev_sub, str) and varname == self.ilev_sub:
                other_name = ['H_ilev', 'H_ilev1']
                other_name.remove(varname)  # first element is the other name
                print(f'{other_name[0]} missing in data and is needed to ' +
                      'convert the requested variables to depend on height.' +
                      f' Using {varname} instead.')
                self.variables[other_name[0]] = self.variables[varname]
                # register the other variable by linking to this one
                self[varname] = other_name[0]

            # create pressure level -> km function once per ilev type
            if varname in ['H_ilev', 'H_ilev1'] or varname in self.total_ilev:
                if varname in ['H_ilev', 'H_ilev1']:  # create custom interp
                    new_varname = 'P'+coord_list[-1][1:]
                    kms = getattr(self, '_km_'+coord_list[-1])
                    kms_max = getattr(self, '_km_'+coord_list[-1]+'_max')
                    kms_min = getattr(self, '_km_'+coord_list[-1]+'_min')
                    # perform unit conversion if needed
                    if self.variables[varname]['units'] != 'km':
                        self[varname+'km_ijk[km]'] = varname + '_ijk'
                        km_interp = self[varname+'km_ijk']
                    else:
                        km_interp = self[varname+'_ijk']
                    # Import and call custom interpolator
                    units = 'm/m'
                    self[new_varname], interp_ijk = RU.PLevelInterp(
                        km_interp, coord_dict['time']['data'],
                        coord_dict['lon']['data'], coord_dict['lat']['data'],
                        coord_dict[coord_list[-1]]['data'], units, kms,
                        [kms_min, kms_max])
                    # kms is a 1D array of the median height values in km
                else:  # define by function composition
                    new_varname = varname.split('_ilev')[0]
                    units = self.variables[varname]['units']
                    # substitute kms array if height was also substituted
                    if isinstance(self.ilev_sub, str):
                        other_name = ['ilev', 'ilev1']
                        other_name.remove(self.ilev_sub[2:])
                        kms = getattr(self, '_km_'+other_name[0])
                        kms_max = getattr(self, '_km_'+other_name[0]+'_max')
                        kms_min = getattr(self, '_km_'+other_name[0]+'_min')
                        self[new_varname] = varname+'(P'+other_name[0][1:]+')'
                        interp_ijk = self[new_varname]
                        self[new_varname].meta['arg_units'] = \
                            self['P'+other_name[0][1:]].meta['arg_units']
                    else:
                        kms = getattr(self, '_km_'+coord_list[-1])
                        kms_max = getattr(self, '_km_'+coord_list[-1]+'_max')
                        kms_min = getattr(self, '_km_'+coord_list[-1]+'_min')
                        self[new_varname] = varname+'(P'+coord_list[-1][1:]+')'
                        interp_ijk = self[new_varname]
                        self[new_varname].meta['arg_units'] = \
                            self['P'+coord_list[-1][1:]].meta['arg_units']
                self.variables[new_varname] = {'units': units, 'data': key}

                # Register in kamodo object with a dedicated logic
                if gridded_int:
                    self = RU.register_griddedPlev(
                        self, new_varname, units, interp_ijk, coord_dict, kms,
                        [kms_min, kms_max])
            return

    return MODEL
