'''
Written by Rebecca Ringuette, 2021
'''
from datetime import datetime, timedelta, timezone

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


def ts_to_hrs(time_val, filedate):
    '''Convert utc timestamp to hours since midnight on filedate.'''
    return (datetime.utcfromtimestamp(time_val).replace(tzinfo=timezone.utc) -
            filedate).total_seconds()/3600.


# times from file converted to seconds since midnight of filedate
# plotting input times will be datetime strings of format 'YYYY-MM-DD HH:mm:ss'
# filedate is self.filedate from iri object
# converts to hours since midnight of filedate for plotting
def MODEL():

    from kamodo import Kamodo
    from netCDF4 import Dataset
    from os.path import basename
    from numpy import array, transpose, NaN, unique, append, zeros, abs, diff
    from numpy import where
    from time import perf_counter
    from astropy.constants import R_earth
    from kamodo_ccmc.readers.reader_utilities import regdef_4D_interpolators
    from kamodo_ccmc.readers.reader_utilities import regdef_3D_interpolators

    class MODEL(Kamodo):
        '''IRI model data reader.

        Inputs:
            full_filename3d: a string representing the file pattern of the
                model output data.
                Note: This reader takes the full filename of the 3D output
                file, typically of the naming convention
                file_dir+'IRI.3D.YYYYDDD.nc',
                where YYYY is the four digit year and DDD is the three digit
                day of year (e.g. 2017148 for May 28, 2017).
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
        def __init__(self, full_filename3d, variables_requested=[],
                     printfiles=False, filetime=False, gridded_int=True,
                     fulltime=True, verbose=False, **kwargs):
            super(MODEL, self).__init__(**kwargs)
            self.modelname = 'IRI'
            t0 = perf_counter()

            # collect filenames
            filename = basename(full_filename3d)
            file_dir = full_filename3d.split(filename)[0]
            if '.2D.' in filename:  # require that input filename be 3D file
                full_filename2d = full_filename3d
                # can't replace in place
                f = full_filename3d.replace('.2D.', '.3D.')
                full_filename3d = f
            else:
                full_filename2d = full_filename3d.replace('.3D.', '.2D.')
            self.filename = full_filename3d + ',' + full_filename2d

            # establish time attributes first
            iri3D = Dataset(full_filename3d, 'r')
            # convert to hours since midnight of file
            time = array(iri3D.variables['time']) / 60.
            # datetime object for midnight on date
            self.filedate = datetime(int(filename[-10:-6]), 1, 1, 0, 0, 0
                                     ).replace(tzinfo=timezone.utc) + \
                timedelta(days=int(filename[-6:-3]) - 1)
            # strings with timezone info chopped off (UTC anyway).
            # Format: ‘YYYY-MM-DD HH:MM:SS’
            self.datetimes = [
                (self.filedate+timedelta(hours=time[0])).isoformat(
                    sep=' ')[:19],
                (self.filedate+timedelta(hours=time[-1])).isoformat(
                    sep=' ')[:19]]
            self.filetimes = [datetime.timestamp(datetime.strptime(
                dt, '%Y-%m-%d %H:%M:%S').replace(tzinfo=timezone.utc)) for dt
                in self.datetimes]   # utc timestamp
            self.dt = diff(time).max()*3600.  # convert time resolution to sec

            if filetime and not fulltime:
                iri3D.close()
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

                file_pattern = file_dir + 'IRI.3D.*.nc'  # returns a string
                files = sorted(glob(file_pattern))  # method may change for AWS
                filenames = unique([basename(f) for f in files])

                # find closest file by utc timestamp
                # iri has an open time at the end
                # need a beginning time from the closest file
                # files are automatically sorted by YYMMDD
                # next file is next in the list
                current_idx = where(filenames == filename)[0]
                if current_idx+1 == len(files):
                    if verbose:
                        print('No later file available.')
                    filecheck = False
                    if filetime:
                        iri3D.close()
                        return
                else:
                    # +1 for adding an end time
                    min_file = file_dir + filenames[current_idx+1][0]
                    kamodo_test = MODEL(min_file, filetime=True,
                                        fulltime=False)
                    time_test = abs(kamodo_test.filetimes[0] -
                                    self.filetimes[1])
                    # if nearest file time at least within one timestep
                    if time_test <= self.dt:
                        filecheck = True
                        self.datetimes[1] = kamodo_test.datetimes[0]
                        self.filetimes[1] = kamodo_test.filetimes[0]

                        # time only version if returning time for searching
                        if filetime:
                            iri3D.close()
                            return  # return object with additional time

                        # get kamodo object with same requested variables
                        if verbose:
                            print(f'Took {perf_counter()-t0:.3f}s to find ' +
                                  'closest file.')
                        kamodo_neighbor = MODEL(
                            min_file, variables_requested=variables_requested,
                            fulltime=False)
                        short_data = kamodo_neighbor.short_data
                        if verbose:
                            print(f'Took {perf_counter()-t0:.3f}s to get ' +
                                  'data from closest file.')
                    else:
                        if verbose:
                            print('No later file found within ' +
                                  f'{diff(time).max()*3600.:.1f}s.')
                        filecheck = False
                        if filetime:
                            iri3D.close()
                            return

            # perform initial check on variables_requested list
            if len(variables_requested) > 0 and fulltime and \
                    variables_requested != 'all':
                test_list = [value[0] for key, value in model_varnames.items()]
                err_list = [item for item in variables_requested if item not in
                            test_list]
                if len(err_list) > 0:
                    print('Variable name(s) not recognized:', err_list)

            # collect variable list
            iri2D = Dataset(full_filename2d, 'r')
            if len(variables_requested) > 0 and variables_requested != 'all':
                gvar_list_2d = [key for key, value in model_varnames.items()
                                if value[0] in variables_requested and
                                key in iri2D.variables.keys()]
                gvar_list_3d = [key for key, value in model_varnames.items()
                                if value[0] in variables_requested and
                                key in iri3D.variables.keys()]

                # check for variables requested but not available
                if len(gvar_list_2d)+len(gvar_list_3d) != \
                        len(variables_requested):
                    err_list = [value[0] for key, value in
                                model_varnames.items() if value[0] in
                                variables_requested and key not in gvar_list_2d
                                and key not in gvar_list_3d]
                    if len(err_list) > 0:
                        print('Some requested variables are not available:',
                              err_list)
            else:  # only input variables on the avoid_list if requested
                gvar_list_2d = [key for key in iri2D.variables.keys()
                                if key in model_varnames.keys()]
                gvar_list_3d = [key for key in iri3D.variables.keys()
                                if key in model_varnames.keys()]
                # returns list of variables included in data files
                if not fulltime and variables_requested == 'all':
                    self.var_dict = {value[0]: value[1:] for key, value in
                                     model_varnames.items() if key in
                                     gvar_list_2d+gvar_list_3d}
                    iri3D.close()
                    iri2D.close()
                    return

            # store data for each variable desired
            variables_2d = {model_varnames[var][0]: {
                'units': model_varnames[var][-1],
                'data': array(iri2D.variables[var])} for var in gvar_list_2d}
            variables_3d = {model_varnames[var][0]: {
                'units': model_varnames[var][-1],
                'data': array(iri3D.variables[var])} for var in gvar_list_3d}
            variables = variables_3d
            for key in variables_2d:
                variables[key] = variables_2d[key]

            # prepare and return data
            if not fulltime:
                iri3D.close()
                iri2D.close()
                variables['time'] = self.filetimes[0]
                self.short_data = variables
                return

            # Store coordinate data as class attributes
            if filecheck:
                # new time in hours since midnight
                new_time = ts_to_hrs(short_data['time'], self.filedate)
                self._time = append(time, new_time)
            else:
                self._time = time

            # collect data and make dimensional grid from 3D file
            lon = array(iri3D.variables['lon'])  # 0 to 360
            lon_le180 = where(lon <= 180)[0]
            lon_ge180 = where(lon >= 180)[0]  # repeat 180 for -180 values
            self._lon = lon - 180.
            self._lat = array(iri3D.variables['lat'])
            self._height = array(iri3D.variables['ht'])
            iri3D.close()   # close netCDF4 files
            iri2D.close()

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
                if len(variables[varname]['data'].shape) == 3:
                    if filecheck:  # if neighbor found
                        # append data for last time stamp
                        data_shape = list(variables[varname]['data'].shape)
                        data_shape[0] += 1  # add space for time
                        new_data = zeros(data_shape)
                        # put in current data
                        new_data[:-1, :, :] = variables[varname]['data']
                        # add in data for additional time
                        new_data[-1, :, :] =\
                            short_data[varname]['data'][0, :, :]
                        variables[varname]['data'] = new_data  # save

                    # shift longitude
                    data_shape = list(variables[varname]['data'].shape)
                    new_data = zeros(data_shape)
                    new_data[:, :, :len(lon_ge180)] = \
                        variables[varname]['data'][:, :, lon_ge180]
                    new_data[:, :, len(lon_ge180)-1:] = \
                        variables[varname]['data'][:, :, lon_le180]
                    # (t, lat, lon) -> (t, lon, lat)
                    variable = transpose(new_data, (0, 2, 1))
                    self.variables[varname] = dict(
                        units=variables[varname]['units'], data=variable)
                    self.register_3D_variable(self.variables[varname]['units'],
                                              self.variables[varname]['data'],
                                              varname, gridded_int)
                elif len(variables[varname]['data'].shape) == 4:
                    if filecheck:
                        # append data for last time stamp, transpose
                        data_shape = list(variables[varname]['data'].shape)
                        data_shape[0] += 1  # add space for time
                        new_data = zeros(data_shape)
                        # put in current data
                        new_data[:-1, :, :, :] = variables[varname]['data']
                        # add in data for additional time
                        new_data[-1, :, :, :] =\
                            short_data[varname]['data'][0, :, :, :]
                        variables[varname]['data'] = new_data  # save

                    # shift longitude
                    data_shape = list(variables[varname]['data'].shape)
                    new_data = zeros(data_shape)
                    new_data[:, :, :, :len(lon_ge180)] = \
                        variables[varname]['data'][:, :, :, lon_ge180]
                    new_data[:, :, :, len(lon_ge180)-1:] = \
                        variables[varname]['data'][:, :, :, lon_le180]
                    # (t, h, lat, lon) -> (t, lon, lat, h)
                    variable = transpose(new_data, (0, 3, 2, 1))
                    self.variables[varname] = dict(
                        units=variables[varname]['units'], data=variable)
                    self.register_4D_variable(self.variables[varname]['units'],
                                              self.variables[varname]['data'],
                                              varname, gridded_int)
            if verbose:
                print(f'Took {perf_counter()-t_reg:.5f}s to register ' +
                      f'{len(varname_list)} variables.')
            if verbose:
                print(f'Took a total of {perf_counter()-t0:.5f}s to kamodofy' +
                      f' {len(varname_list)} variables.')

        # define and register a 3D variable
        def register_3D_variable(self, units, variable, varname, gridded_int):
            """Registers a 3d interpolator with 3d signature"""

            # define and register the interpolators
            xvec_dependencies = {'time': 'hr', 'lon': 'deg', 'lat': 'deg'}
            self = regdef_3D_interpolators(self, units, variable, self._time,
                                           self._lon, self._lat, varname,
                                           xvec_dependencies, gridded_int)
            return

        # define and register a 4D variable
        def register_4D_variable(self, units, variable, varname, gridded_int):
            """Registers a 4d interpolator with 4d signature"""

            # define and register the fast interpolator
            xvec_dependencies = {'time': 'hr', 'lon': 'deg', 'lat': 'deg',
                                 'height': 'km'}
            self = regdef_4D_interpolators(self, units, variable, self._time,
                                           self._lon, self._lat, self._height,
                                           varname, xvec_dependencies,
                                           gridded_int)
            return
    return MODEL
