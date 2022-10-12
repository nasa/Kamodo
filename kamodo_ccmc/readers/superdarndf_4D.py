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
    from numpy import array, NaN, diff
    from time import perf_counter
    from kamodo_ccmc.readers.reader_utilities import Functionalize_Dataset

    class MODEL(Kamodo):
        '''SuperDARN model data reader for the default coordinate grid.

        Inputs:
            file_prefixdf: a string representing the file pattern of the
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
        def __init__(self, file_prefixdf, variables_requested=[],
                     printfiles=False, filetime=False, gridded_int=True,
                     fulltime=True, verbose=False, **kwargs):
            super(MODEL, self).__init__(**kwargs)
            self.modelname = 'SuperDARN_df'
            t0 = perf_counter()

            # check for prepared file of given prefix
            t0 = perf_counter()
            if isfile(file_prefixdf + '_df.nc'):   # file already prepared!
                cdf_file = file_prefixdf + '_df.nc'  # input file name
                self.conversion_test = True
            else:   # file not prepared,  prepare it
                from kamodo_ccmc.readers.superdarn_tocdf import \
                    convert_files
                cdf_file = convert_files(file_prefixdf)
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
                int(file_prefixdf[-8:-4]), int(file_prefixdf[-4:-2]),
                int(file_prefixdf[-2:]), 0, 0, 0
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
            if len(variables_requested) > 0 and variables_requested != 'all':
                gvar_list = [key for key, value in model_varnames.items()
                             if value[0] in variables_requested and
                             key in cdf_data.variables.keys()]

                # check for variables requested but not available
                if len(gvar_list) != len(variables_requested):
                    err_list = [value[0] for key, value in
                                model_varnames.items() if value[0] in
                                variables_requested and key not in gvar_list]
                    if len(err_list) > 0:
                        print('Some requested variables are not available:',
                              err_list)
            else:  # only input variables on the avoid_list if requested
                gvar_list = [key for key in cdf_data.variables.keys()
                             if key in model_varnames.keys()]
                # returns list of variables included in data files
                if not fulltime and variables_requested == 'all':
                    self.var_dict = {value[0]: value[1:] for key, value in
                                     model_varnames.items() if key in
                                     gvar_list}
                    cdf_data.close()
                    return

            # store data for each variable desired
            variables = {model_varnames[var][0]: {
                'units': model_varnames[var][-1],
                'data': array(cdf_data.variables[var])} for var in gvar_list}

            # prepare and return data
            if not fulltime:
                cdf_data.close()
                variables['time'] = self.filetimes[0]
                self.short_data = variables
                return

            # Store coordinate data as class attributes
            self._time = time
            self._lon = array(cdf_data.variables['lon'])  # -180 to 180
            self._lat = array(cdf_data.variables['lat'])
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
                coord_str = [value[3]+value[4] for key, value in
                             model_varnames.items() if value[0] == varname][0]
                # set coordinate dictionary
                if len(variables[varname]['data'].shape) == 3:
                    coord_dict = {'time': {'units': 'hr', 'data': self._time},
                                  'lon': {'units': 'deg', 'data': self._lon},
                                  'lat': {'units': 'deg', 'data': self._lat}}
                elif len(variables[varname]['data'].shape) == 1:
                    coord_dict = {'time': {'units': 'hr', 'data': self._time}}
                # functionalize the variable dataset
                self = Functionalize_Dataset(self, coord_dict, varname,
                                             self.variables[varname],
                                             gridded_int, coord_str)
            if verbose:
                print(f'Took {perf_counter()-t_reg:.5f}s to register ' +
                      f'{len(varname_list)} variables.')
            if verbose:
                print(f'Took a total of {perf_counter()-t0:.5f}s to kamodofy' +
                      f' {len(varname_list)} variables.')
            return
    return MODEL
