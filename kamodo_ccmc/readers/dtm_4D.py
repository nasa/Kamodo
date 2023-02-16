# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 12:59:08 2022

@author: rringuet
DTM model reader
"""

model_varnames = {'Temp_exo': ['T_exo', 'Exospheric temperature', 0, 'GDZ',
                               'sph', ['time', 'lon', 'lat'], 'K'],
                  'Temp': ['T', 'temperature', 0, 'GDZ', 'sph',
                           ['time', 'lon', 'lat', 'height'], 'K'],
                  'DEN': ['rho', 'Total mass density', 0, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'height'], 'g/cm**3'],
                  'MU': ['m_avgmol', 'Mean molecular mass', 0, 'GDZ', 'sph',
                         ['time', 'lon', 'lat', 'height'], 'g'],
                  'H': ['N_H', 'Atomic hydrogen partial density', 0, 'GDZ',
                        'sph', ['time', 'lon', 'lat', 'height'], 'g/cm**3'],
                  'He': ['N_He', 'Atomic helium partial density', 0, 'GDZ',
                         'sph', ['time', 'lon', 'lat', 'height'], 'g/cm**3'],
                  'O': ['N_O', 'Atomic oxygen partial density', 0, 'GDZ',
                        'sph', ['time', 'lon', 'lat', 'height'], 'g/cm**3'],
                  'N2': ['N_N2', 'Molecular nitrogen partial density', 0,
                         'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                         'g/cm**3'],
                  'O2': ['N_O2', 'Molecular oxygen partial density', 0, 'GDZ',
                         'sph', ['time', 'lon', 'lat', 'height'], 'g/cm**3']
                  }


def MODEL():
    from time import perf_counter
    from glob import glob
    from os.path import basename, isfile
    from numpy import array, unique, NaN, append, transpose, where
    from datetime import datetime, timezone
    from netCDF4 import Dataset
    from kamodo import Kamodo
    import kamodo_ccmc.readers.reader_utilities as RU

    class MODEL(Kamodo):
        '''DTM model data reader.

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
            - This model reader is the most basic example of what is required
              in a model reader. The only 'data wrangling' needed is a simple
              longitude wrapping and numpy array transposition to get the
              coordinate order correct.
            - DTM outputs files are given in one netCDF file per day.
            - The files are small and contain multiple time steps per file, so
              interpolation method 2 is chosen. The standard SciPy interpolator
              is used.
        '''

        def __init__(self, file_dir, variables_requested=[],
                     filetime=False, verbose=False, gridded_int=True,
                     printfiles=False, **kwargs):
            super(MODEL, self).__init__()
            self.modelname = 'DTM'

            # first, check for file list, create if DNE
            list_file = file_dir + self.modelname + '_list.txt'
            time_file = file_dir + self.modelname + '_times.txt'
            self.times, self.pattern_files = {}, {}
            if not isfile(list_file) or not isfile(time_file):
                t0 = perf_counter()  # begin timer
                # figure out types of files present (2DTEC, 3DALL, 3DLST, etc)
                files = sorted(glob(file_dir+'*.nc'))
                patterns = sorted(unique([basename(f)[:-11] for f in
                                          files]))  # cut off date
                self.filename = ''.join([f+',' for f in files])[:-1]
                self.filedate = datetime.strptime(
                    basename(files[0])[-10:-3]+' 00:00:00', '%Y%j %H:%M:%S'
                    ).replace(tzinfo=timezone.utc)

                # establish time attributes
                for p in patterns:  # only one pattern
                    # get list of files to loop through later
                    pattern_files = sorted(glob(file_dir+p+'*.nc'))
                    self.pattern_files[p] = pattern_files
                    self.times[p] = {'start': [], 'end': [], 'all': []}

                    # loop through to get times, one day per file
                    for f in range(len(pattern_files)):
                        cdf_data = Dataset(pattern_files[f])
                        # minutes since 12am EACH file -> hrs since 12am 1st f
                        tmp = array(cdf_data.variables['time'])/60. + f*24.
                        self.times[p]['start'].append(tmp[0])
                        self.times[p]['end'].append(tmp[-1])
                        self.times[p]['all'].extend(tmp)
                        cdf_data.close()
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
                return  # return times as is to prevent infinite recursion

            # store variables
            self.missing_value = NaN
            self.varfiles = {}  # store which variable came from which file
            self.gvarfiles = {}  # store file variable name similarly
            self.err_list = []
            self.variables = {}

            # perform initial check on variables_requested list
            if len(variables_requested) > 0 and variables_requested != 'all':
                test_list = [value[0] for key, value in model_varnames.items()]
                err_list = [item for item in variables_requested if item not in
                            test_list]
                if len(err_list) > 0:
                    print('Variable name(s) not recognized:', err_list)
                for item in err_list:
                    variables_requested.remove(item)
                if len(variables_requested) == 0:
                    return

            # there is only one pattern for DTM, so just save the one grid
            p = list(self.pattern_files.keys())[0]
            pattern_files = self.pattern_files[p]

            # get coordinate grids from first file
            cdf_data = Dataset(pattern_files[0], 'r')
            self._lat = array(cdf_data.variables['lat'])  # -90 to 90
            lon = array(cdf_data.variables['lon'])  # 0 to 360
            lon_le180 = list(where(lon <= 180)[0])  # 0 to 180
            lon_ge180 = list(where((lon >= 180) & (lon < 360.))[0])
            self._lon_idx = lon_ge180 + lon_le180
            self._lon = lon - 180.
            self._height = array(cdf_data.variables['ht'])  # km

            # check var_list for variables not possible in this file set
            if len(variables_requested) > 0 and\
                    variables_requested != 'all':
                gvar_list = [key for key in model_varnames.keys()
                             if key in cdf_data.variables.keys() and
                             model_varnames[key][0] in variables_requested]
                if len(gvar_list) != len(variables_requested):
                    err_list = [value[0] for key, value in
                                model_varnames.items()
                                if key not in cdf_data.variables.keys() and
                                value[0] in variables_requested]
                    self.err_list.extend(err_list)  # add to master list
            else:
                gvar_list = [key for key in model_varnames.keys()
                             if key in cdf_data.variables.keys()]

            # store which file these variables came from
            self.varfiles[p] = [model_varnames[key][0] for
                                key in gvar_list]
            self.gvarfiles[p] = gvar_list
            cdf_data.close()

            # print message if variables not found
            if len(self.err_list) > 0:
                print('Some requested variables are not available: ',
                      self.err_list)

            # collect all possible variables in set of files and return
            if variables_requested == 'all':
                var_list = self.varfiles[p]
                self.var_dict = {value[0]: value[1:] for key, value in
                                 model_varnames.items() if value[0] in
                                 var_list}
                return

            # initialize storage structure
            self.variables = {model_varnames[gvar][0]: {
                'units': model_varnames[gvar][-1], 'data': p} for gvar in
                self.gvarfiles[p]}

            # option to print files
            if printfiles:
                print(f'{len(self.filename)} Files:')
                files = self.filename.split(',')
                for f in files:
                    print(f)

            # register interpolators for each variable
            t_reg = perf_counter()
            # store original list b/c gridded interpolators change keys list
            varname_list = [key for key in self.variables.keys()]
            for varname in varname_list:
                self.register_variable(varname, gridded_int)

            if verbose:
                print(f'Took {perf_counter()-t_reg:.5f}s to register ' +
                      f'{len(varname_list)} variables.')
            if verbose:
                print(f'Took a total of {perf_counter()-t0:.5f}s to kamodofy' +
                      f' {len(varname_list)} variables.')

        def register_variable(self, varname, gridded_int):
            """Registers an interpolator with proper signature"""

            # determine which file the variable came from, retrieve the coords
            key = self.variables[varname]['data']
            gvar = [key for key, value in model_varnames.items() if
                    value[0] == varname][0]  # variable name in file
            coord_list = [value[-2] for key, value in
                          model_varnames.items() if value[0] == varname][0]
            coord_dict = {'time': {'units': 'hr',
                                   'data': self.times[key]['all']}}
            # get the correct coordinates
            coord_dict['lon'] = {'data': self._lon, 'units': 'deg'}
            coord_dict['lat'] = {'data': self._lat, 'units': 'deg'}
            if len(coord_list) == 4:
                coord_dict['height'] = {'data': self._height, 'units': 'km'}
            coord_str = [value[3]+value[4] for key, value in
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
                if len(data.shape) == 3:
                    variable = transpose(data, (0, 2, 1))
                elif len(data.shape) == 4:
                    variable = transpose(data, (0, 3, 2, 1))
                return variable[:, self._lon_idx]

            self = RU.Functionalize_Dataset(
                self, coord_dict, varname, self.variables[varname],
                gridded_int, coord_str, interp_flag=2, func=func,
                times_dict=self.times[key])

    return MODEL
