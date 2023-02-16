'''
Written by Rebecca Ringuette, 2021
'''

# variable name in file: [standardized variable name, descriptive term, units]
model_varnames = {'Ne': ['N_e', 'electron number density',
                         0, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                         '1/m**3'],
                  'Te': ['T_e', 'electron temperature',
                         1, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                         'K'],
                  'Ti': ['T_i', 'ion temperature',
                         2, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                         'K'],
                  'Tn': ['T_n', 'neutral temperature',
                         3, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                         'K'],
                  'O+': ['N_Oplus', 'number density of atomic oxygen ion',
                         4, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                         '1/m**3'],
                  'H+': ['N_Hplus', 'number density of atomic hydrogen ion',
                         5, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                         '1/m**3'],
                  'He+': ['N_Heplus', 'number density of atomic helium ion',
                          6, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                          '1/m**3'],
                  'O2+': ['N_O2plus', 'number density of molecular oxygen ion',
                          7, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                          '1/m**3'],
                  'NO+': ['N_NOplus', 'number density of molecular nitric ' +
                          'oxide', 8, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                     'height'], '1/m**3'],
                  'N+': ['N_Nplus', 'number density of atomic nitrogen ion',
                         9, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                         '1/m**3'],
                  'TEC': ['TEC', 'vertical total electron content (height ' +
                          'integrated from bottom to top boundary)',
                          10, 'GDZ', 'sph', ['time', 'lon', 'lat'],
                          '10**16/m**2'],
                  'NmF2': ['NmF2', 'maximum electron number density in F2 ' +
                           'layer', 11, 'GDZ', 'sph', ['time', 'lon', 'lat'],
                           '1/m**3'],
                  'HmF2': ['HmF2', 'height of maximum electron number ' +
                           'density in F2 layer', 12, 'GDZ', 'sph',
                           ['time', 'lon', 'lat'], 'km']
                  }


# times from file converted to seconds since midnight of filedate
# plotting input times will be datetime strings of format 'YYYY-MM-DD HH:mm:ss'
# filedate is self.filedate from iri object
# converts to hours since midnight of filedate for plotting
def MODEL():

    from kamodo import Kamodo
    from netCDF4 import Dataset
    from glob import glob
    from os.path import basename, isfile
    from numpy import array, transpose, NaN, unique
    from numpy import where, append
    from time import perf_counter
    from datetime import datetime, timedelta, timezone
    import kamodo_ccmc.readers.reader_utilities as RU

    class MODEL(Kamodo):
        '''IRI model data reader.

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
            - IRI output files are given in multiple netCDF files per day.
            - The files are small and contain multiple time steps per file, so
              interpolation method 2 is chosen. The standard SciPy interpolator
              is used.
        '''
        def __init__(self, file_dir, variables_requested=[],
                     printfiles=False, filetime=False, gridded_int=True,
                     verbose=False, **kwargs):
            super(MODEL, self).__init__(**kwargs)
            self.modelname = 'IRI'
            t0 = perf_counter()

            # first, check for file list, create if DNE
            list_file = file_dir + self.modelname + '_list.txt'
            time_file = file_dir + self.modelname + '_times.txt'
            self.times, self.pattern_files = {}, {}
            if not isfile(list_file) or not isfile(time_file):
                # collect filenames
                files = sorted(glob(file_dir+'*.nc'))
                patterns = unique([basename(f)[:-10] for f in files])  # 2D, 3D
                self.filename = ''.join([f+',' for f in files])[:-1]

                # establish time attributes
                for p in patterns:
                    # get list of files to loop through later
                    pattern_files = sorted(glob(file_dir+p+'*.nc'))
                    self.pattern_files[p] = pattern_files
                    self.times[p] = {'start': [], 'end': [], 'all': []}

                    # loop through to get times
                    for f in range(len(pattern_files)):
                        cdf_data = Dataset(pattern_files[f])
                        tmp = array(cdf_data.variables['time'])/60. + \
                            float(f)*24.  # hrs since midnite 1st file
                        self.times[p]['start'].append(tmp[0])
                        self.times[p]['end'].append(tmp[-1])
                        self.times[p]['all'].extend(tmp)
                        cdf_data.close()
                    self.times[p]['start'] = array(self.times[p]['start'])
                    self.times[p]['end'] = array(self.times[p]['end'])
                    self.times[p]['all'] = array(self.times[p]['all'])

                # datetime object for midnight on date
                self.filedate = datetime(int(files[0][-10:-6]), 1, 1, 0, 0, 0
                                         ).replace(tzinfo=timezone.utc) + \
                    timedelta(days=int(files[0][-6:-3]) - 1)

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
                err_list = [item for item in variables_requested if item not in
                            test_list]
                if len(err_list) > 0:
                    print('Variable name(s) not recognized:', err_list)
                for item in err_list:
                    variables_requested.remove(item)
                if len(variables_requested) == 0:
                    return

            # loop through file patterns for var mapping
            self.gvarfiles, self.varfiles, self.err_list = {}, {}, []
            for p in self.pattern_files.keys():
                # check var_list for variables not possible in this file set
                cdf_data = Dataset(self.pattern_files[p][0], 'r')
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
                cdf_data.close()
                # store which file these variables came from
                self.varfiles[p] = [model_varnames[key][0] for
                                    key in gvar_list]
                self.gvarfiles[p] = gvar_list

            # clean up error list and then take action
            var_list = []
            for p in self.varfiles.keys():
                var_list.extend(self.varfiles[p])
            err_list = [var for var in self.err_list if var not in var_list]
            if len(err_list) > 0:
                print('Some requested variables are not available: ',
                      err_list)

            # collect all possible variables in set of files and return
            if variables_requested == 'all':
                self.var_dict = {value[0]: value[1:] for key, value in
                                 model_varnames.items() if value[0] in
                                 var_list}
                return

            # get coordinates and data from files
            self.variables = {}
            for p in self.pattern_files.keys():
                # get coordinates from first file of each type
                cdf_data = Dataset(self.pattern_files[p][0], 'r')
                lon = array(cdf_data.variables['lon'])
                lon_le180 = list(where(lon <= 180)[0])  # 0 to 180
                lon_ge180 = list(where((lon >= 180) & (lon < 360.))[0])
                tmp = lon - 180.  # now -180. to +180.
                setattr(self, '_lon_'+p, tmp)
                setattr(self, '_lon_idx_'+p, lon_ge180+lon_le180)
                setattr(self, '_lat_'+p,
                        array(cdf_data.variables['lat']))
                if '3D' in p:
                    setattr(self, '_height_'+p,
                            array(cdf_data.variables['ht']))
                cdf_data.close()
                # initialize variable dictionaries
                for var in self.gvarfiles[p]:
                    self.variables[model_varnames[var][0]] = {
                        'units': model_varnames[var][-1], 'data': p}

            # store a few items
            self.missing_value = NaN
            self._registered = 0
            if verbose:
                print(f'Took {perf_counter()-t0:.6f}s to read in data')
            if printfiles:
                print(self.filename)

            # register interpolators for each requested variable
            t_reg = perf_counter()
            # store original list b/c gridded interpolators change key list
            varname_list = list(self.variables.keys())
            for varname in varname_list:
                # determine which time grid applies
                self.register_variable(varname, gridded_int)
            if verbose:
                print(f'Took {perf_counter()-t_reg:.5f}s to register ' +
                      f'{len(varname_list)} variables.')
            if verbose:
                print(f'Took a total of {perf_counter()-t0:.5f}s to kamodofy' +
                      f' {len(varname_list)} variables.')

        # define and register a variable
        def register_variable(self, varname, gridded_int):
            """Registers an interpolator with proper signature"""

            # define and register the interpolators
            key = self.variables[varname]['data']
            coord_list = [value[5] for key, value in model_varnames.items()
                          if value[0] == varname][0]
            coord_dict = {'time': {'units': 'hr',
                                   'data': self.times[key]['all']}}
            if 'lat' in coord_list:   # 3D variables come from neutral file
                coord_dict['lon'] = {'units': 'deg', 'data':
                                     getattr(self, '_lon_'+key)}
                coord_dict['lat'] = {'units': 'deg', 'data':
                                     getattr(self, '_lat_'+key)}
            if 'height' in coord_list:
                coord_dict['height'] = {'units': 'km', 'data':
                                        getattr(self, '_height_'+key)}
            lon_idx = getattr(self, '_lon_idx_'+key)
            gvar = [key for key, value in model_varnames.items() if
                    value[0] == varname][0]  # variable name in file
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
                    # time, lat, lon -> time, lon, lat
                    variable = transpose(data, (0, 2, 1))
                elif len(data.shape) == 4:
                    # time, height, lat, lon -> time, lon, lat, height
                    variable = transpose(data, (0, 3, 2, 1))
                return variable[:, lon_idx]

            # functionalize the variable data (time chunked interpolation)
            self = RU.Functionalize_Dataset(
                self, coord_dict, varname, self.variables[varname],
                gridded_int, coord_str, interp_flag=2, func=func,
                times_dict=self.times[key])
            return
    return MODEL
