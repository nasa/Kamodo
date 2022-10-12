'''
Written by Rebecca Ringuette, 2021
'''
from datetime import datetime, timedelta, timezone

# variable name in file: [standardized variable name, descriptive term, units]
model_varnames = {'PED': ['Sigma_P', 'Pedersen conductance ', 0, 'MAG', 'sph',
                          ['time', 'lon', 'lat'], 'S'],
                  'HALL': ['Sigma_H', 'hall Conductance',
                           1, 'MAG', 'sph', ['time', 'lon', 'lat'],
                           'S'],
                  'PHI': ['phi', 'Electric potential',
                          2, 'MAG', 'sph', ['time', 'lon', 'lat'],
                          'kV'],
                  'EEAST': ['E_east', 'electric field in eastern direction ' +
                            '(increasing longitude) ',
                            3, 'MAG', 'sph', ['time', 'lon', 'lat'],
                         'mV/m'],
                  'ENORTH': ['E_north', 'electric field in north direction ' +
                             '(increasing latitude) ',
                         4, 'MAG', 'sph', ['time', 'lon', 'lat'],
                         'mV/m'],
                  'JEAST': ['J_east', 'electric current in eastern direction' +
                            '(increasing longitude) (height integrated ' +
                            'current density)',
                         5, 'MAG', 'sph', ['time', 'lon', 'lat'],
                         'A/m'],
                  'JNORTH': ['J_north', 'electric current in north direction' +
                             ' (increasing latitude) (height integrated ' +
                             'current density)',
                          6, 'MAG', 'sph', ['time', 'lon', 'lat'],
                          'A/m'],
                  'EFLUX': ['Phi', 'energy flux',
                          7, 'MAG', 'sph', ['time', 'lon', 'lat'],
                          'mW/m**2'],
                  'JHEAT': ['J_heat', 'Joule hating rate ', 8, 'MAG', 'sph',
                            ['time', 'lon', 'lat'], 'mW/m**2'],
                  'JRIN': ['J_Rin', 'Input radial current',
                         9, 'MAG', 'sph', ['time', 'lon', 'lat'],
                         'mA/m**2'],
                  'JROUT': ['J_Rout', 'model-generated radial current ' +
                            'would be identical to J_Rin if the model were ' +
                            'perfect',
                          10, 'MAG', 'sph', ['time', 'lon', 'lat'],
                          'mA/m**2'],
                  }


def ts_to_hrs(time_val, filedate):
    '''Convert utc timestamp to hours since midnight on filedate.'''
    return (datetime.utcfromtimestamp(time_val).replace(tzinfo=timezone.utc) -
            filedate).total_seconds()/3600.


# times from file converted to seconds since midnight of filedate
# plotting input times will be datetime strings of format 'YYYY-MM-DD HH:mm:ss'
# filedate is self.filedate from adelphi object
# converts to hours since midnight of filedate for plotting
def MODEL():

    from kamodo import Kamodo
    from netCDF4 import Dataset
    from os.path import basename, isfile
    from numpy import array, NaN, diff, count_nonzero, broadcast_to, mean
    from numpy import delete, sum, cos, sin, pi
    from time import perf_counter
    import kamodo_ccmc.readers.reader_utilities as RU

    class MODEL(Kamodo):
        '''ADELPHI model data reader.

        Inputs:
            file_prefix: a string representing the file pattern of the
                model output data.
                Note: This reader takes the file prefix of the output
                file, typically of the naming convention
                file_dir+'ADELPHI_2D_MAG_YYYYMMDD',
                where YYYY is the four digit year, MM is the two digit month,
                and DD is the two digit day. (e.g. 20170528 for May 28, 2017).
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
        def __init__(self, file_prefix, variables_requested=[],
                     printfiles=False, filetime=False, gridded_int=True,
                     fulltime=True, verbose=False, **kwargs):
            super(MODEL, self).__init__(**kwargs)
            self.modelname = 'ADELPHI'
            t0 = perf_counter()

            # check for prepared file of given prefix
            t0 = perf_counter()
            cdf_file = file_prefix + '.nc'  # input file name
            if isfile(file_prefix + '.nc'):   # file already prepared!
                self.conversion_test = True
            else:   # file not prepared,  prepare it
                filename = basename(file_prefix)
                file_dir = file_prefix.split(filename)[0]
                from kamodo_ccmc.readers.adelphi_tocdf import convert_all
                self.conversion_test = convert_all(file_dir)
            self.filename = cdf_file

            # establish time attributes first
            cdf_data = Dataset(cdf_file, 'r')
            # convert to hours since midnight of file
            time = array(cdf_data.variables['time'])  # in hours
            # datetime object for midnight on date
            self.filedate = datetime.strptime(cdf_data.filedate,
                                              '%Y-%m-%d %H:%M:%S').replace(
                                                  tzinfo=timezone.utc)
            # strings with timezone info chopped off (UTC anyway).
            # Format: ‘YYYY-MM-DD HH:MM:SS’
            self.datetimes = [
                (self.filedate+timedelta(seconds=int(time[0]*3600.))).isoformat(
                    sep=' ')[:19],
                (self.filedate+timedelta(seconds=int(time[-1]*3600.))).isoformat(
                    sep=' ')[:19]]
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

            # Store coordinate data as class attributes
            self._time = time
            self._lon = array(cdf_data.variables['lon'])  # 0 to 360
            self._lat = array(cdf_data.variables['lat'])  # -180 to 180
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
                if len(variables[varname]['data'].shape) == 3:
                    self.variables[varname] = dict(
                        units=variables[varname]['units'],
                        data=variables[varname]['data'])
                    self.register_3D_variable(varname, gridded_int)
            if verbose:
                print(f'Took {perf_counter()-t_reg:.5f}s to register ' +
                      f'{len(varname_list)} variables.')
            if verbose:
                print(f'Took a total of {perf_counter()-t0:.5f}s to kamodofy' +
                      f' {len(varname_list)} variables.')

        def vector_average3D(top, shape_list, varname, latval, lon):
            '''find vector average at pole for array with shape (time, lon,
            height)'''
        
            # CAN NOT TEST BECAUSE NUMBERS ARE ALL ZEROS!!!
            # find net x and y components, final array shapes are (t, lon, ht)
            lon_arr = broadcast_to(lon[1:], (shape_list[0], shape_list[1]-1))  # reverse and .T?
            # need to put 'old' shape at end in broadcast_to call ^
            xval = sum(top*cos((lon_arr+180.)*pi/180.), axis=1)  # same shape
            yval = sum(top*sin((lon_arr+180.)*pi/180.), axis=1)  # as t and z
            xarr = broadcast_to(xval, (shape_list[1]-1, shape_list[0])).T
            yarr = broadcast_to(yval, (shape_list[1]-1, shape_list[0])).T
        
            # convert to proper unit vector (see wiki on spherical coordinates)
            if 'EAST' in varname:
                # Zonal/east components -> convert to psi_hat vector (lon)
                # -xsin(psi)+ycos(psi), psi = longitude (0 to 360)
                new_top = -xarr*sin((lon_arr+180.)*pi/180.) +\
                    yarr*cos((lon_arr+180.)*pi/180.)
            elif 'NORTH' in varname:
                # meridional/north -> convert to theta_hat vector (latitude)
                # xcos(psi)cos(theta)+ysin(psi)cos(theta)
                # sin(theta) is always zero at the poles
                # theta = latitude (0 to 180), psi = longitude (0 to 360)
                new_top = xarr * cos((lon_arr+180.) * pi/180.) *\
                    cos((90.-latval)*pi/180.) + yarr *\
                    sin((lon_arr+180.) * pi/180.) * cos((90.-latval) * pi/180.)
            return new_top

        # define and register a 3D variable
        def register_3D_variable(self, varname, gridded_int):
            """Registers a 3d interpolator with 3d signature"""

            # define and register the interpolators
            coord_dict = {'time': {'units': 'hr', 'data': self._time},
                          'lon': {'units': 'deg', 'data': self._lon},
                          'lat': {'units': 'deg', 'data': self._lat}}
            coord_str = [value[3]+value[4] for key, value in
                         model_varnames.items() if value[0] == varname][0]

            # slice off zeros at each pole before averaging and functionalizing
            # find first latitude row that is nonzero, only near poles
            # this purposefully ignores the buffer rows near the equator
            # this is not in the file converter due to the dynamic nature of 
            #   the latitude grid.
            # south pole at beginning
            SP_idx, slice_idx = 1, []
            zero_check = count_nonzero(
                self.variable[varname]['data'][:, :, SP_idx])
            while zero_check == 0:
                slice_idx.append(SP_idx)
                SP_idx += 1
                zero_check = count_nonzero(
                    self.variable[varname]['data'][:, :, SP_idx])
            # north pole at the end
            NP_idx = -2
            zero_check = count_nonzero(
                self.variable[varname]['data'][:, :, NP_idx])
            while zero_check == 0:
                slice_idx.append(NP_idx)
                NP_idx -= 1
                zero_check = count_nonzero(
                    self.variable[varname]['data'][:, :, NP_idx])
            # remove 'extra' latitude values from coordinate grid and from data
            coord_dict['lat']['data'] = delete(coord_dict['lat']['data'],
                                               slice_idx)
            self.variable[varname]['data'] = delete(
                self.variable[varname]['data'], slice_idx, axis=2)

            # perform scalar averaging at the poles
            new_shape = (len(coord_dict['time']['data']),
                         len(coord_dict['lon']['data']),
                         len(coord_dict['lat']['data']))
            if varname not in ['E_east', 'E_north', 'J_east', 'J_north']:
                # south pole
                top = mean(self.variable[varname]['data'][:, 1:, 1],
                           axis=1)  # same shape as time axis
                self.variable[varname]['data'][:, 1:, 0] = broadcast_to(
                    top, (new_shape[1]-1, new_shape[0])).T
                # north pole
                top = mean(self.variable[varname]['data'][:, 1:, -2],
                           axis=1)  # same shape as time axis
                self.variable[varname]['data'][:, 1:, -1] = broadcast_to(
                    top, (new_shape[1]-1, new_shape[0])).T
            else:  # perform vector averaging at the poles
                # calculate net vector magnitude for south pole
                self.variable[varname]['data'][:, 1:, 0] = \
                    self.vector_average3D(
                        self.variable[varname]['data'][:, 1:, 1],
                        list(self.variable[varname]['data'].shape), varname,
                        coord_dict['lat']['data'][0],
                        coord_dict['lon']['data'])
                # repeat for north pole
                self.variable[varname]['data'][:, 1:, -1] = \
                    self.vector_average3D(
                        self.variable[varname]['data'][:, 1:, -2],
                        list(self.variable[varname]['data'].shape), varname,
                        coord_dict['lat']['data'][-1],
                        coord_dict['lon']['data'])

            # now functionalize the cleaned data
            self = RU.Functionalize_Dataset(self, coord_dict, varname,
                                            self.variables[varname],
                                            gridded_int, coord_str)
            return

    return MODEL
