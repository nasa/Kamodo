'''
Written by Rebecca Ringuette, 2021
'''
from datetime import datetime, timedelta, timezone

# variable name in file: [standardized variable name, descriptive term, units]
model_varnames = {"V": ['V', 'Electric potential', 0, 'SM', 'sph',
                        ['time', 'lon', 'lat'], 'kV'],
                  "theta_v": ['theta_v', 'Azimuthal angle of convection ' +
                              'velocity', 1, 'SM', 'sph', ['time', 'lon',
                                                            'lat'], 'deg'],
                  "v": ['v', 'SMnetiude of convection velocity', 2, 'SM',
                        'sph', ['time', 'lon', 'lat'], 'm/s'],
                  # remaining variables are time series
                  "theta_Btilt": ['theta_Btilt', 'Dipole tilt', 3, 'SM',
                                  'sph', ['time'], 'deg'],
                  'E_sw': ['E_sw', 'Electric field of the solar wind', 4,
                           'SM', 'sph', ['time'], 'mV/m'],
                  'theta_B': ['theta_B', 'IMF clock angle', 5, 'SM', 'sph',
                              ['time'], 'deg']}


def ts_to_hrs(time_val, filedate):
    '''Convert utc timestamp to hours since midnight on filedate.'''
    return (datetime.utcfromtimestamp(time_val).replace(tzinfo=timezone.utc) -
            filedate).total_seconds()/3600.


# times from file converted to seconds since midnight of filedate
# plotting input times will be datetime strings of format 'YYYY-MM-DD HH:mm:ss'
# filedate is self.filedate from superdarn object
# converts to hours since midnight of filedate for plotting
def MODEL():

    from kamodo import Kamodo
    from netCDF4 import Dataset
    from os.path import basename, isfile
    from numpy import array, where, NaN, unique, append, zeros, abs, diff
    from numpy import linspace, ndarray
    from time import perf_counter
    import kamodo_ccmc.readers.reader_utilities as RU

    class MODEL(Kamodo):
        '''SuperDARN model data reader for the default coordinate grid.

        Inputs:
            full_filenamedf: a string representing the file pattern of the
                model output data for the default coordinate grid format.
                Note: This reader takes the file prefix, typically of
                the naming convention
                file_dir+'modelYYYYMMDD',
                where YYYY is the four digit year, MM is the two digit month,
                and DD is the two digit day (e.g. 20200505 for May 5, 2020).
                All files with the default grid should be in the chosen
                directory, with all files from the equal-area grid
                (if generated) in a separate directory.
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
        def __init__(self, file_prefixea, variables_requested=[],
                     printfiles=False, filetime=False, gridded_int=True,
                     fulltime=True, verbose=False, **kwargs):
            super(MODEL, self).__init__(**kwargs)
            self.modelname = 'SuperDARN_ea'
            t0 = perf_counter()

            # check for prepared file of given prefix
            t0 = perf_counter()
            if isfile(file_prefixea + '_ea.nc'):   # file already prepared!
                cdf_file = file_prefixea + '_ea.nc'  # input file name
                self.conversion_test = True
            else:   # file not prepared,  prepare it
                from kamodo_ccmc.readers.superdarn_tocdf import \
                    convert_files
                cdf_file = convert_files(file_prefixea)
                self.conversion_test = True
            filename = basename(cdf_file)
            file_dir = cdf_file.split(filename)[0]
            self.filename = cdf_file

            # establish time attributes first
            cdf_data = Dataset(cdf_file, 'r')
            # hours since midnight of file
            time = array(cdf_data.variables['time'])
            # datetime object for midnight on date
            self.filedate = datetime(
                int(file_prefixea[-8:-4]), int(file_prefixea[-4:-2]),
                int(file_prefixea[-2:]), 0, 0, 0
                ).replace(tzinfo=timezone.utc)
            # strings with timezone info chopped off (UTC anyway).
            # Format: ‘YYYY-MM-DD HH:MM:SS’
            self.datetimes = [
                (self.filedate + timedelta(seconds=int(time[0]*3600.))
                 ).isoformat(sep=' ')[:19],
                (self.filedate + timedelta(seconds=int(time[-1]*3600.))
                 ).isoformat(sep=' ')[:19]]
            self.filetimes = [datetime.timestamp(datetime.strptime(
                dt, '%Y-%m-%d %H:%M:%S').replace(tzinfo=timezone.utc)) for dt
                in self.datetimes]   # utc timestamp
            self.dt = diff(time).max()*3600.  # convert time resolution to sec

            if filetime:
                cdf_data.close()
                return  # return times as is to prevent recursion

            # if variables are given as integers, convert to standard names
            if len(variables_requested) > 0:
                if isinstance(variables_requested[0], int):
                    tmp_var = [value[0] for key, value in
                               model_varnames.items()
                               if value[2] in variables_requested]
                    variables_requested = tmp_var

            # perform initial check on variables_requested list
            if len(variables_requested) > 0 and fulltime and \
                    variables_requested != 'all':
                test_list = [value[0] for key, value in model_varnames.items()]
                err_list = [item for item in variables_requested if item not in
                            test_list]
                if len(err_list) > 0:
                    print('Variable name(s) not recognized:', err_list)

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
                if not fulltime and variables_requested == 'all':
                    self.var_dict = {value[0]: value[1:] for key, value in
                                     model_varnames.items() if key in
                                     gvar_list}
                    cdf_data.close()
                    return

            # retrieve latitude grid
            self._lat = array(cdf_data.variables['lat'])

            # store data for each 1D variable desired
            variables = {model_varnames[var][0]: {
                'units': model_varnames[var][-1],
                'data': array(cdf_data.variables[var])} for var in gvar_list
                if len(model_varnames[var][5]) == 1}
            # gives an empty dictionary if no 1D variables are requested

            # collect 3D data
            var3D_list = [model_varnames[var][0] for var in gvar_list
                          if len(model_varnames[var][5]) == 3]
            lat_keys = [str(lat_val).replace('.', '_').replace('-', 'n')
                        if lat_val < 0 else 'p'+str(lat_val).replace('.', '_')
                        for lat_val in self._lat]  # e.g. 'n75_5' for -75.5
            for var in var3D_list:
                variables[var] = {'units': model_varnames[var][-1]}
                variables[var]['data'] = {
                    lat_val: array(cdf_data.groups[lat_key][var])
                    for lat_key, lat_val in zip(lat_keys, self._lat)}

            # retrieve dictionary of longitude grids. Keys are lat values.
            self._lon = {lat_val: array(cdf_data.groups[lat_key]['lon']) for
                         lat_key, lat_val in zip(lat_keys, self._lat)}

            # Store coordinate data as class attributes
            self._time = time
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
            varname_list = [key for key in variables.keys()]
            self.variables = {}
            for varname in varname_list:
                self.variables[varname] = dict(
                    units=variables[varname]['units'],
                    data=variables[varname]['data'])
                if isinstance(variables[varname]['data'], dict):  # 3D var
                    self.register_3D_variable(varname, gridded_int)
                elif isinstance(variables[varname]['data'], ndarray):  # 1D var
                    self.register_1D_variable(varname, gridded_int)
            if verbose:
                print(f'Took {perf_counter()-t_reg:.5f}s to register ' +
                      f'{len(varname_list)} variables.')
            if verbose:
                print(f'Took a total of {perf_counter()-t0:.5f}s to kamodofy' +
                      f' {len(varname_list)} variables.')

        # define and register a 3D variable
        def register_3D_variable(self, varname, gridded_int):
            """Registers a 3d interpolator with 3d signature"""

            # define the coordinate dictionary
            coord_dict = {'time': {'units': 'hr', 'data': self._time},
                          'lon': {'units': 'deg', 'data': self._lon},
                          'lat': {'units': 'deg', 'data': self._lat}}
            coord_data = [value['data'] for key, value in coord_dict.items()]
            coord_units = {key: value['units'] for key, value in
                           coord_dict.items()}

            # variable is a dictionary = {lat_val: array}
            # need an interpolator for each lat_val
            # then a 1D interpolator for the output values vs latitude
            from kamodo_ccmc.readers.superdarnea_interp import custom_interp
            interp = custom_interp(*coord_dict,
                                   self.variables[varname]['data'],
                                   self.variables[varname]['units'])

            # Register in kamodo object, with gridded version if desired
            self = RU.register_interpolator(self, varname, interp,
                                            coord_units)
            if gridded_int:
                fake_lon = linspace(-180, 180, 20)  # sample lon grid
                coord_data[-2] = fake_lon  # replace lon dict with sample grid
                fake_data = zeros((2, 2, 2))  # saves execution time
                self.variables[varname+'_ijk'] = {
                    'units': self.variables[varname]['units'],
                    'data': fake_data}
                gridded_interpolator = RU.define_griddedinterp(
                    self.variables[varname+'_ijk'], coord_units,
                    coord_data, interp)
                self = RU.register_interpolator(
                    self, varname+'_ijk', gridded_interpolator,
                    coord_units)
            return

        # define and register a 1D variable
        def register_1D_variable(self, varname, gridded_int):
            """Registers a 4d interpolator with 4d signature"""

            # define and register the fast interpolator
            coord_dict = {'time': {'units': 'hr', 'data': self._time}}
            coord_str = [value[3]+value[4] for key, value in
                         model_varnames.items() if value[0] == varname][0]
            self = RU.Functionalize_Dataset(self, coord_dict, varname,
                                         self.variables[varname],
                                         gridded_int, coord_str)
            return
    return MODEL
