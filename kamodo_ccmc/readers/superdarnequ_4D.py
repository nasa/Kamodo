'''
Written by Rebecca Ringuette, 2021
'''
from datetime import datetime, timezone

# variable name in file: [standardized variable name, descriptive term, units]
model_varnames = {"V": ['V', 'Electric potential', 0, 'SM', 'sph',
                        ['time', 'lon', 'lat'], 'kV'],
                  "theta_v": ['theta_v', 'Azimuthal angle of convection ' +
                              'velocity', 1, 'SM', 'sph', ['time', 'lon',
                                                           'lat'], 'deg'],
                  "v": ['v', 'magnitude of convection velocity', 2, 'SM',
                        'sph', ['time', 'lon', 'lat'], 'm/s'],
                  # remaining variables are time series
                  "theta_Btilt": ['theta_Btilt', 'Dipole tilt', 3, 'SM',
                                  'sph', ['time'], 'deg'],
                  'E_sw': ['E_sw', 'Electric field of the solar wind', 4,
                           'SM', 'sph', ['time'], 'mV/m'],
                  'theta_B': ['theta_B', 'IMF clock angle', 5, 'SM', 'sph',
                              ['time'], 'deg']}


# times from file converted to seconds since midnight of filedate
# plotting input times will be datetime strings of format 'YYYY-MM-DD HH:mm:ss'
# filedate is self.filedate from superdarn object
# converts to hours since midnight of filedate for plotting
def MODEL():

    from kamodo import Kamodo
    from netCDF4 import Dataset
    from glob import glob
    from os.path import basename, isfile
    from numpy import array, NaN, unique, zeros, linspace
    from time import perf_counter
    import kamodo_ccmc.readers.reader_utilities as RU

    class MODEL(Kamodo):
        '''SuperDARN model data reader for the equal area coordinate grid.

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
            - The SuperDARN equal are grid model output is given in one ascii
              file per timestep per N/S hemisphere. The file converter combines
              these files into one netCDF4 file per timestep.
            - The data is only given within 40-60 degrees of the poles, so a
              buffer row of the same values is added in the data to avoid
              losing the most equatorward ring of data in the interpolation.
            - The outputs do not provide values at the poles, so scalar
              averaging is used to determine these values. (The three
              spatially dependent variables are all scalars.)
            - The converted files are small and are created with one time step
              per file, so interpolation method 1 is chosen for spatially
              dependent variables. Interpolation method 0 is chosen for the
              spatially independent variables, meaning the entire time series
              is read into memory and the SciPy interp1d interpolator is used.
            - For the spatially dependent variables, a custom interpolator is
              required to properly handle the longitude grid's changes with
              changing latitudes and the changing range of the latitude grid
              with time. The logic of the interpolator is almost identical to
              that used in the lazy interpolation logic included in the
              reader_utilities.py script, except applied to latitude instead of
              time.
        '''
        def __init__(self, file_dir, variables_requested=[],
                     filetime=False, verbose=False, gridded_int=True,
                     printfiles=False, **kwargs):
            super(MODEL, self).__init__(**kwargs)
            self.modelname = 'SuperDARN_equ'
            t0 = perf_counter()

            # first, check for file list, create if DNE
            list_file = file_dir + self.modelname + '_list.txt'
            time_file = file_dir + self.modelname + '_times.txt'
            self.times, self.pattern_files = {}, {}
            if not isfile(list_file) or not isfile(time_file):
                # make sure all files are converted
                from kamodo_ccmc.readers.superdarn_tocdf import \
                    convert_all
                tmp = convert_all(file_dir, '*equ.txt')

                # continue with converted files
                files = sorted(glob(file_dir+'*_equ.nc'))
                self.filename = ''.join([f+',' for f in files])[:-1]
                patterns = unique([basename(f)[:-20] for f in files])
                self.filedate = datetime.strptime(
                    basename(files[0])[-20:-12], '%Y%m%d').replace(
                        tzinfo=timezone.utc)

                # establish time attributes from filenames
                for p in patterns:
                    # get list of files to loop through later
                    pattern_files = sorted(glob(file_dir+p+'*_equ.nc'))
                    self.pattern_files[p] = pattern_files
                    self.times[p] = {'start': [], 'end': [], 'all': []}
                    if pattern_files[0][-12] == '-':
                        format_string = '%Y%m%d-%H%M'
                    else:
                        format_string = '%Y%m%d_%H%M'

                    # loop through to get times, one file per time step
                    self.times[p]['start'] = array([
                        RU.str_to_hrs(f[-20:-7], self.filedate,
                                      format_string=format_string) for f in
                        pattern_files])
                    self.times[p]['end'] = self.times[p]['start'].copy()
                    self.times[p]['all'] = self.times[p]['start'].copy()

                # create time list file if DNE
                RU.create_timelist(list_file, time_file, self.modelname,
                                   self.times, self.pattern_files,
                                   self.filedate)
            else:  # read in data and time grids from file list
                self.times, self.pattern_files, self.filedate, self.filename =\
                    RU.read_timelist(time_file, list_file)
            if filetime:
                return  # return times as is to prevent infinite recursion
            # only one pattern, so simplifying code
            p = list(self.pattern_files.keys())[0]
            cdf_data = Dataset(self.pattern_files[p][0])

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

            # collect variable list
            sample_groupname = list(cdf_data.groups.keys())[0]
            file_variables = list(cdf_data.variables.keys()) +\
                list(cdf_data.groups[sample_groupname].variables.keys())
            if len(variables_requested) > 0 and variables_requested != 'all':
                gvar_list = [key for key, value in model_varnames.items()
                             if value[0] in variables_requested and
                             key in file_variables]

                # check for variables requested but not available
                if len(gvar_list) != len(variables_requested):
                    err_list = [value[0] for key, value in
                                model_varnames.items() if value[0] in
                                variables_requested and key not in gvar_list]
                    if len(err_list) > 0:
                        print('Some requested variables are not available:',
                              err_list)
            else:  # only input variables on the avoid_list if requested
                gvar_list = [key for key in file_variables
                             if key in model_varnames.keys()]
                # returns list of variables included in data files
                if variables_requested == 'all':
                    self.var_dict = {value[0]: value[1:] for key, value in
                                     model_varnames.items() if key in
                                     gvar_list}
                    cdf_data.close()
                    return

            # retrieve latitude grid
            self._lat = array(cdf_data.variables['lat'])

            # store kay for each variable desired
            self.variables = {model_varnames[var][0]: {
                'units': model_varnames[var][-1],
                'data': p} for var in gvar_list}

            # retrieve dictionary of longitude grids. Keys are lat values.
            # e.g. 'n75_5' for -75.5
            self.lat_keys = [str(lat_val).replace('.', '_').replace('-', 'n')
                             if lat_val < 0 else 'p'+str(lat_val).replace(
                                     '.', '_') for lat_val in self._lat]
            self._lon_dict = {lat_val: array(cdf_data.groups[lat_key]['lon'])
                              for lat_key, lat_val in zip(self.lat_keys,
                                                          self._lat)}
            cdf_data.close()

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
            varname_list = [key for key in self.variables.keys()]
            for varname in varname_list:
                self.register_variables(varname, gridded_int)
            if verbose:
                print(f'Took {perf_counter()-t_reg:.5f}s to register ' +
                      f'{len(varname_list)} variables.')
            if verbose:
                print(f'Took a total of {perf_counter()-t0:.5f}s to kamodofy' +
                      f' {len(varname_list)} variables.')

        # define and register the variable
        def register_variables(self, varname, gridded_int):
            '''Functionalizes the indicated dataset.'''

            # determine coordinate variables and xvec by coord list
            key = self.variables[varname]['data']
            coord_list = [value[5] for key, value in model_varnames.items()
                          if value[0] == varname][0]
            coord_dict = {'time': {'units': 'hr',
                                   'data': self.times[key]['all']}}
            gvar = [key for key, value in model_varnames.items() if
                    value[0] == varname][0]
            coord_str = [value[3]+value[4] for key, value in
                         model_varnames.items() if value[0] == varname][0]

            # retrieve coordinate grids
            if len(coord_list) == 1:  # read in entire time series
                data = []
                for f in self.pattern_files[key]:
                    cdf_data = Dataset(f)
                    tmp = array(cdf_data.variables[gvar])
                    cdf_data.close()
                    data.extend(tmp)
                self.variables[varname]['data'] = array(data)
                self = RU.Functionalize_Dataset(self, coord_dict, varname,
                                                self.variables[varname],
                                                gridded_int, coord_str,
                                                interp_flag=0)  # single array
                return

            # rest of data is time + 2D spatial
            # sample lon grid, real lon grid is a dictionary
            fake_lon = linspace(-180, 180, len(self._lon_dict.keys()))
            coord_dict['lon'] = {'units': 'deg', 'data': fake_lon}
            coord_dict['lat'] = {'units': 'deg', 'data': self._lat}

            # variable is a dictionary = {lat_val: array}
            # need an interpolator for each lat_val
            # then a 1D interpolator for the output values vs latitude
            from kamodo_ccmc.readers.superdarnequ_interp import custom_interp

            def func(i):  # i is the file/time number
                # get data from file(s)
                cdf_data = Dataset(self.pattern_files[key][i])
                lat = array(cdf_data.variables['lat'])
                lat_keys = [str(lat_val).replace('.', '_').replace('-', 'n')
                            if lat_val < 0 else 'p'+str(lat_val).replace(
                                    '.', '_') for lat_val in lat]
                lon_dict = {lat_val: array(cdf_data.groups[lat_key]['lon'])
                            for lat_key, lat_val in zip(lat_keys, lat)}
                data = {lat_val: array(cdf_data.groups[lat_key][gvar])
                        for lat_key, lat_val in zip(lat_keys, lat)}
                cdf_data.close()

                # assign custom interpolator
                interp = custom_interp(lon_dict, lat, data)
                return interp

            # functionalize the variable dataset
            tmp = self.variables[varname]
            tmp['data'] = zeros((2, 2, 2))  # saves execution time
            self = RU.Functionalize_Dataset(
                self, coord_dict, varname, tmp, gridded_int, coord_str,
                interp_flag=1, func=func, func_default='custom')
            return
    return MODEL
