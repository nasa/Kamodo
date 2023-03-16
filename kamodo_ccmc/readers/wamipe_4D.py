from datetime import datetime, timezone


# varnames in cdf files are standardized (value[0])
model_varnames = {'den400': ['rho_400km', 'Density at 400 km.',
                             0, 'GDZ', 'sph', ['time', 'lon', 'lat'],
                             'kg/m**3'],
                  'ON2': ['ON2', 'mmr or ratio?',
                          0, 'GDZ', 'sph', ['time', 'lon', 'lat'], 'm/m'],
                  'height': ['H_ilev', 'Height dependent upon pressure level',
                             0, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                             'm'],
                  'temp_neutral': ['T_n_ilev', 'neutral temperature', 0, 'GDZ',
                                   'sph', ['time', 'lon', 'lat', 'ilev'],
                                   'K'],
                  'temp_neutral_2': ['T_n', 'neutral temperature', 0, 'GDZ',
                                     'sph', ['time', 'lon', 'lat', 'height'],
                                     'K'],
                  'O_Density': ['N_O_ilev', 'Oxygen number density', 0, 'GDZ',
                                'sph', ['time', 'lon', 'lat', 'ilev'],
                                '1/m**3'],
                  'O_Density_2': ['N_O', 'Oxygen number density', 0, 'GDZ',
                                  'sph', ['time', 'lon', 'lat', 'height'],
                                  '1/m**3'],
                  'O2_Density': ['N_O2_ilev', 'Molecular oxygen number ' +
                                 'density', 0, 'GDZ', 'sph',
                                 ['time', 'lon', 'lat', 'ilev'], '1/m**3'],
                  'O2_Density_2': ['N_O2', 'Molecular oxygen number density',
                                   0, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                     'height'], '1/m**3'],
                  'N2_Density': ['N_N2_ilev', 'Molecular nitrogen number ' +
                                 'density', 0, 'GDZ', 'sph',
                                 ['time', 'lon', 'lat', 'ilev'], '1/m**3'],
                  'N2_Density_2': ['N_N2', 'Molecular nitrogen number density',
                                   0, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                     'height'], '1/m**3'],
                  'u_neutral': ['v_neast_ilev', 'Eastward component of the ' +
                                'neutral wind velocity', 0, 'GDZ', 'sph',
                                ['time', 'lon', 'lat', 'ilev'], 'm/s'],
                  'u_neutral_2': ['v_neast', 'Eastward component of the ' +
                                  'neutral wind velocity', 0, 'GDZ', 'sph',
                                  ['time', 'lon', 'lat', 'height'], 'm/s'],
                  'v_neutral': ['v_nnorth_ilev', 'Northward component of the' +
                                ' neutral wind velocity', 0, 'GDZ', 'sph',
                                ['time', 'lon', 'lat', 'ilev'], 'm/s'],
                  'v_neutral_2': ['v_nnorth', 'Northward component of the ' +
                                  'neutral wind velocity', 0, 'GDZ', 'sph',
                                  ['time', 'lon', 'lat', 'height'], 'm/s'],
                  'w_neutral': ['v_nup_ilev', 'Upwards component of the ' +
                                'neutral wind velocity', 0, 'GDZ', 'sph',
                                ['time', 'lon', 'lat', 'ilev'], 'm/s'],
                  'w_neutral_2': ['v_nup', 'Upwards component of the ' +
                                  'neutral wind velocity', 0, 'GDZ', 'sph',
                                  ['time', 'lon', 'lat', 'height'], 'm/s'],
                  'tec': ['TEC', 'Total electron content', 0, 'GDZ', 'sph',
                          ['time', 'lon', 'lat'], '10**16/m**2'],
                  'NmF2': ['NmF2', 'Maximum number density in the F2 layer',
                           0, 'GDZ', 'sph', ['time', 'lon', 'lat'], '1/m**3'],
                  'HmF2': ['HmF2', 'Height of the maximum number density in ' +
                           'the F2 layer', 0, 'GDZ', 'sph',
                           ['time', 'lon', 'lat'], 'km'],
                  'O_plus_density': ['N_Oplus', 'O+ number density', 0, 'GDZ',
                                     'sph', ['time', 'lon', 'lat', 'height'],
                                     '1/m**3'],
                  'H_plus_density': ['N_Hplus', 'H+ number density', 0, 'GDZ',
                                     'sph', ['time', 'lon', 'lat', 'height'],
                                     '1/m**3'],
                  'He_plus_density': ['N_Heplus', 'He+ number density', 0,
                                      'GDZ', 'sph', ['time', 'lon', 'lat',
                                                     'height'], '1/m**3'],
                  'N_plus_density': ['N_Nplus', 'N+ number density', 0, 'GDZ',
                                     'sph', ['time', 'lon', 'lat', 'height'],
                                     '1/m**3'],
                  'NO_plus_density': ['N_NOplus', 'NO+ number density', 0,
                                      'GDZ', 'sph', ['time', 'lon', 'lat',
                                                     'height'], '1/m**3'],
                  'O2_plus_density': ['N_O2plus', 'O2+ number density', 0,
                                      'GDZ', 'sph', ['time', 'lon', 'lat',
                                                     'height'], '1/m**3'],
                  'N2_plus_density': ['N_O2plus', 'N2+ number density', 0,
                                      'GDZ', 'sph', ['time', 'lon', 'lat',
                                                     'height'], '1/m**3'],
                  'ion_temperature': ['T_i', 'Ion temperature', 0, 'GDZ',
                                      'sph', ['time', 'lon', 'lat', 'height'],
                                      'K'],
                  'electron_temperature': ['T_e', 'Electron temperature', 0,
                                           'GDZ', 'sph', ['time', 'lon', 'lat',
                                                          'height'], 'K'],
                  'eastward_exb_velocity': ['v_eastExB', 'Eastward component' +
                                            ' of the ExB velocity', 0, 'GDZ',
                                            'sph', ['time', 'lon', 'lat',
                                                    'height'], 'm/s'],
                  'northward_exb_velocity': ['v_northExB', 'Northward ' +
                                             'component of the ExB velocity',
                                             0, 'GDZ', 'sph',
                                             ['time', 'lon', 'lat', 'height'],
                                             'm/s'],
                  'upward_exb_velocity': ['v_upExB', 'Upward component' +
                                          ' of the ExB velocity', 0, 'GDZ',
                                          'sph', ['time', 'lon', 'lat',
                                                  'height'], 'm/s']
                  }


def MODEL():
    from time import perf_counter
    from glob import glob
    from os.path import basename, isfile
    from numpy import array, unique, NaN, append, linspace, where, median
    from numpy import zeros
    from netCDF4 import Dataset
    from kamodo import Kamodo
    import kamodo_ccmc.readers.reader_utilities as RU

    class MODEL(Kamodo):
        '''WAM-IPE model data reader.

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
            - WAM-IPE output files are given in multiple netCDF files per day.
              No attempt is made to combine these files, but some
              pre-processing is needed to calculate the median height across
              the entire dataset (time, lon, lat) for a given pressure level
              value. If a tarball of files is detected in the absense of the
              netCDF files, then the tarball will be decompressed and the
              netCDF files extracted automatically.
            - WAM-IPE data is given in several coordinate systems, one
              depending on pressure level - a unitless representation of height
              corresponding to a given atmospheric pressure. The pressure
              level values are not given in the data, so an arbitrary integer
              representation is chosen. 
            - Pressure level inversion is performed in the reader_utilities
              script, specifically the PLevelInterp function. Two versions of
              all variables that depend on pressure level are created: the
              original and one dependent on height, which is created through
              Kamodo's function composition feature.
            - The files are small and contain one time step per file, so
              interpolation method 1 is chosen. The standard SciPy interpolator
              is used.
        '''
        def __init__(self, file_dir, variables_requested=[],
                     filetime=False, verbose=False, gridded_int=True,
                     printfiles=False, **kwargs):
            super(MODEL, self).__init__()
            self.modelname = 'WAMIPE'
            t0 = perf_counter()

            # first, check for file list, create if DNE
            list_file = file_dir + self.modelname + '_list.txt'
            time_file = file_dir + self.modelname + '_times.txt'
            self.times, self.pattern_files = {}, {}
            if not isfile(list_file) or not isfile(time_file):
                # collect filenames
                files = sorted(glob(file_dir+'*.nc'))
                if len(files) == 0:  # find tar files and untar them
                    print('Decompressing files...this may take a moment.')
                    import tarfile
                    tar_files = glob(file_dir+'*.tar')
                    for file in tar_files:
                        tar = tarfile.open(file)
                        tar.extractall(file_dir)
                        tar.close()
                    files = sorted(glob(file_dir+'*.nc'))

                # create h0 file containing km_ilev
                from kamodo_ccmc.readers.wamipe_tocdf import convert_all
                tmp = convert_all(file_dir)

                # continue
                files = sorted(glob(file_dir+'*.nc'))
                h0_file = [f for f in files if 'h0' in f]
                if len(h0_file) > 0:
                    files.remove(h0_file[0])
                patterns = unique([basename(f)[:-19] for f in files])
                self.filename = ''.join([f+',' for f in files])[:-1]
                self.filedate = datetime.strptime(
                    basename(files[0])[-18:-10]+' 00:00:00', '%Y%m%d %H:%M:%S'
                    ).replace(tzinfo=timezone.utc)

                # establish time attributes from filenames
                for p in patterns:
                    # get list of files to loop through later
                    pattern_files = sorted(glob(file_dir+p+'*.nc'))
                    h0_file = [f for f in pattern_files if 'h0' in f]
                    if len(h0_file) > 0:
                        pattern_files.remove(h0_file[0])  # should be only one
                    self.pattern_files[p] = pattern_files
                    self.times[p] = {'start': [], 'end': [], 'all': []}

                    # loop through to get times, one file per time
                    time_str = [basename(f)[-18:-3] for f in pattern_files]
                    times = RU.str_to_hrs(time_str, self.filedate,
                                          format_string='%Y%m%d_%H%M%S')
                    self.times[p]['start'] = array(times)
                    self.times[p]['end'] = array(times)
                    self.times[p]['all'] = array(times)

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

            if printfiles:
                print(f'{len(self.filename)} Files:')
                for file in self.filename:
                    print(file)

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

            # collect variable dependency lists
            self.total_ilev = [value[0] for key, value in
                               model_varnames.items() if value[-2][-1] ==
                               'ilev' and value[0] != 'H_ilev']
            self.total_replace = [item.split('_ilev')[0] for item in
                                  self.total_ilev if item != 'H_ilev']
            # dictionary mapping to navigate related variable names
            self.ilev_map = {item1: item2 for item1, item2 in
                             zip(self.total_replace, self.total_ilev)}
            # determine variable mapping
            if len(variables_requested) > 0 and variables_requested != 'all':
                # add ilev version of variables to the list, adding H_ilev
                add_ilev = [var+'_ilev' for var in variables_requested if var
                            in self.ilev_map.keys()]
                if len(add_ilev) > 0:  # add to list
                    add_ilev += ['H_ilev']  # add for plev inversion
                new_var = variables_requested + add_ilev
                short_var = [item for item in new_var if item not
                             in self.ilev_map.keys()]  # remove replaced items

            # collect variable names and coordinate grids per pattern type
            self.varfiles, self.gvarfiles = {}, {}
            self.err_list, self.var_dict = [], {}
            for p in self.pattern_files.keys():
                # check var_list for variables not possible in this file set
                cdf_data = Dataset(self.pattern_files[p][0], 'r')
                if len(variables_requested) > 0 and \
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
                        self.err_list.extend(err_list)  # add to master list
                else:
                    gvar_list = [key for key in model_varnames.keys()
                                 if key in cdf_data.variables.keys()]
                # store which file these variables came from
                self.varfiles[p] = [model_varnames[key][0] for key in
                                    gvar_list]
                self.gvarfiles[p] = gvar_list

                # collect all possible variables in set of files and return
                if variables_requested == 'all':
                    var_dict = {value[0]: value[1:] for key, value in
                                model_varnames.items() if key in gvar_list}
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

                # retrieve coordinate grids, assuming constant in time
                lon = array(cdf_data.variables['lon'])
                lon_le180 = list(where(lon <= 180)[0])  # 0 to 180
                lon_ge180 = list(where(lon >= 180)[0])
                tmp = append(lon, 360.) - 180.  # now -180. to +180.
                setattr(self, '_lon_'+p, tmp)
                setattr(self, '_lon_idx_'+p, lon_ge180+lon_le180)
                setattr(self, '_lat_'+p,
                        array(cdf_data.variables['lat']))
                if 'alt' in cdf_data.variables.keys():
                    setattr(self, '_height_'+p,
                            array(cdf_data.variables['alt']))  # km
                    setattr(self, '_heightunits_'+p, 'km')
                # determine if pressure leve is needed (not included)
                ilev_check = [True if 'ilev' in var else False for var in
                              self.varfiles[p]]
                if sum(ilev_check) > 0:
                    # create a pressure level coordinate grid
                    idx = ilev_check.index(True)  # get first index of 3D var
                    varname = self.gvarfiles[p][idx]
                    grid_length = cdf_data.variables[varname].shape[0]
                    setattr(self, '_ilev_'+p,
                            linspace(1, grid_length, grid_length))

                    # get median km grid from h0 file
                    if 'height' in cdf_data.variables.keys():
                        h0_file = self.pattern_files[p][0][:-18] + 'h0.nc'
                        if isfile(h0_file):
                            cdf_h = Dataset(h0_file)
                            self._km_ilev = array(cdf_h.variables['km_ilev'])
                            self._km_ilev_max = cdf_h.km_max
                            self._km_ilev_min = cdf_h.km_min
                            cdf_h.close()
                cdf_data.close()
            # return var_dict (already created)
            if variables_requested == 'all':
                return

            # loop through patterns to initialize variables dictionary
            for p in self.pattern_files.keys():
                gvar_list = self.gvarfiles[p]
                variables = {model_varnames[key][0]: {
                    'units': model_varnames[key][-1], 'data': p}
                    for key in gvar_list}
                for key in variables.keys():
                    self.variables[key] = variables[key]

            # remove successful variables from err_list
            self.err_list = list(unique(self.err_list))
            self.err_list = [item for item in self.err_list if item not in
                             self.variables.keys()]
            if len(self.err_list) > 0:
                print('Some requested variables are not available in the ' +
                      'files found:\n', self.err_list)

            # register interpolators for each requested variable
            # rearrange to deal with H_ilev and H_ilev1 first if there
            varname_list = list(self.variables.keys())
            if 'H_ilev' in varname_list:
                varname_list.remove('H_ilev')
                varname_list = ['H_ilev'] + varname_list
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
            gvar = [key for key, value in model_varnames.items() if
                    value[0] == varname][0]
            # retrieve coordinate grids
            if 'lat' in coord_list:
                coord_dict['lon'] = {'units': 'deg', 'data':
                                     getattr(self, '_lon_'+key)}
                coord_dict['lat'] = {'units': 'deg', 'data':
                                     getattr(self, '_lat_'+key)}
                lon_idx = getattr(self, '_lon_idx_'+key)
            if 'ilev' in coord_list:
                try:
                    coord_dict['ilev'] = {'units': 'm/m', 'data':
                                          getattr(self, '_ilev_'+key)}
                    if hasattr(self, '_ilevunits_'+key):
                        coord_dict['ilev']['units'] = \
                            getattr(self, '_ilevunits_'+key)
                except:
                    coord_list.remove('ilev')
                    coord_list.append('height')
                    new_varname = varname.split('_ilev')[0]
                    self.variables[new_varname] = self.variables[varname]
                    varname = new_varname
                    if gvar[-2:] == '_2':
                        gvar = gvar[:-2]
            if 'height' in coord_list:
                coord_dict['height'] = {'units': 'km', 'data':
                                        getattr(self, '_height_'+key)}
                if hasattr(self, '_heightunits_'+key):
                    coord_dict['height']['units'] = \
                        getattr(self, '_heightunits_'+key)
            coord_str = [value[3]+value[4] for key, value in
                         model_varnames.items() if value[0] == varname][0]

            # define operations for each variable when given the key
            def func(i):
                '''i is the time slice. WAM-IPE has one time slice per file.'''
                # get data from file
                file = self.pattern_files[key][i]
                cdf_data = Dataset(file)
                data = array(cdf_data.variables[gvar])
                # data wrangling
                if hasattr(cdf_data.variables[gvar], '_FillValue'):
                    fill_value = cdf_data.variables[gvar]._FillValue
                    if fill_value in data:
                        data = where(data == fill_value, NaN, data)
                cdf_data.close()
                return data.T[lon_idx]

            # define and register the interpolators
            # need H functions to be gridded regardless of gridded_int value
            if varname == 'H_ilev':
                gridded_int = True
            self = RU.Functionalize_Dataset(
                self, coord_dict, varname, self.variables[varname],
                gridded_int, coord_str, interp_flag=1, func=func)

            # create pressure level -> km function once per ilev type
            if (varname == 'H_ilev' or 'ilev' in coord_list) and (
                    hasattr(self, '_km_ilev')):  # only if inversion is possibl
                if varname == 'H_ilev':  # create custom interp
                    new_varname = 'Plev'
                    # perform unit conversion if needed
                    if self.variables[varname]['units'] != 'km':
                        self[varname+'km_ijk[km]'] = varname + '_ijk'
                        km_interp = self[varname+'km_ijk']
                    else:
                        km_interp = self[varname+'_ijk']
                    # Import and call custom interpolator
                    units = coord_dict['ilev']['units']
                    self[new_varname], interp_ijk = RU.PLevelInterp(
                        km_interp, coord_dict['time']['data'],
                        coord_dict['lon']['data'], coord_dict['lat']['data'],
                        coord_dict['ilev']['data'], units, self._km_ilev,
                        [self._km_ilev_min, self._km_ilev_max])
                    # kms is a 1D array of the median height values in km
                else:  # define by function composition
                    new_varname = varname.split('_ilev')[0]
                    units = self.variables[varname]['units']
                    self[new_varname] = varname+'(Plev)'
                    interp_ijk = self[new_varname]
                    self[new_varname].meta['arg_units'] = \
                        self['Plev'].meta['arg_units']
                self.variables[new_varname] = {'units': units, 'data': key}

                # create gridded interpolator if requested
                if gridded_int:
                    self = RU.register_griddedPlev(
                        self, new_varname, units, interp_ijk, coord_dict,
                        self._km_ilev, [self._km_ilev_min, self._km_ilev_max])
            return
    return MODEL
