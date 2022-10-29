from datetime import datetime, timezone, timedelta
from numpy import vectorize


# varnames in cdf files are standardized (value[0])
model_varnames = {'den400': ['rho_400km', 'Density at 400 km.',
                           0, 'GDZ', 'sph', ['time', 'lon', 'lat'],
                           'kg/m**3'],
                  'ON2': ['ON2', 'mmr or ratio?',
                          0, 'GDZ', 'sph', ['time', 'lon', 'lat'], 'm/m'],
                  'temp_neutral': ['T_n', 'neutral temperature', 0, 'GDZ',
                                   'sph', ['time', 'lon', 'lat', 'height'],
                                   'K'],
                  'O_Density': ['N_O', 'Oxygen number density', 0, 'GDZ',
                                'sph', ['time', 'lon', 'lat', 'height'],
                                '1/m**3'],
                  'O2_Density': ['N_O2', 'Molecular oxygen number density', 0,
                                 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                'height'], '1/m**3'],
                  'N2_Density': ['N_N2', 'Molecular nitrogen number density',
                                 0, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                   'height'], '1/m**3'],
                  'u_neutral': ['v_neast', 'Eastward component of the ' +
                                'neutral wind velocity', 0, 'GDZ', 'sph',
                                ['time', 'lon', 'lat', 'height'], 'm/s'],
                  'v_neutral': ['v_nnorth', 'Northward component of the ' +
                                'neutral wind velocity', 0, 'GDZ', 'sph',
                                ['time', 'lon', 'lat', 'height'], 'm/s'],
                  'w_neutral': ['v_nup', 'Upwards component of the ' +
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
                                                  'height'], 'm/s']}


def dts_to_hrs(datetime_string, filedate):
    '''Get hours since midnight from datetime string'''
    return (datetime.strptime(datetime_string, '%Y-%m-%d %H:%M:%S'
                              ).replace(tzinfo=timezone.utc) -
            filedate).total_seconds()/3600.


def filename_to_dts(file_dts):
    '''Get datetime string in format "YYYY-MM-DD HH:mm:SS" from filename'''
    str_date = file_dts[:4] + '-' + file_dts[4:6] + '-' + file_dts[6:8]
    str_time = file_dts[9:11] + ':' + file_dts[11:13] + ':' + file_dts[13:]
    return str_date + ' ' + str_time


@vectorize
def filename_to_hrs(file_dts, filedate):
    '''Get hrs since midnight from filename.'''
    dts = filename_to_dts(file_dts)  # datetime string
    dt = datetime.strptime(dts, '%Y-%m-%d %H:%M:%S'
                                                ).replace(tzinfo=timezone.utc)
    return (dt-filedate).total_seconds()/3600.


def dts_to_ts(file_dts):
    '''Get datetime timestamp in UTC from datetime string'''
    return datetime.timestamp(datetime.strptime(file_dts, '%Y-%m-%d %H:%M:%S'
                                                ).replace(tzinfo=timezone.utc))


def hrs_to_ts(hrs, filedate):
    '''Add hours to filedate and return utc timestamp.'''
    return datetime.timestamp(filedate+timedelta(hours=hrs))


def ts_to_hrs(time_val, filedate):
    '''Convert utc timestamp to hours since midnight on filedate.'''
    return (datetime.utcfromtimestamp(time_val).replace(tzinfo=timezone.utc) -
            filedate).total_seconds()/3600.


def MODEL():
    from time import perf_counter
    from glob import glob
    from os.path import basename, isfile
    from numpy import array, unique, NaN, diff, append, linspace, where
    from netCDF4 import Dataset
    from kamodo import Kamodo
    import kamodo_ccmc.readers.reader_utilities as RU

    class MODEL(Kamodo):
        '''WAM-IPE model data reader.

        Inputs:
            file_dir: a string representing the file pattern of the************************
                model output data.
                Note: This reader takes a file pattern of the format
                file_dir+*YYMMDD to load data for an entire day or *YYMMDD_HH
                to load data for one hour, where file_dir is the complete file
                path to the data files, and YYMMDD is the two digit year, two
                digit month, and two digit day in the desired output file names
                (e.g. 150317 for March 17 of 2015). For the hourly option, HH
                is the two digit hour assuming a 24 hour convention.
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
        def __init__(self, file_dir, variables_requested=[],
                     filetime=False, verbose=False, gridded_int=True,
                     printfiles=False, fulltime=True, **kwargs):
            super(MODEL, self).__init__()
            self.modelname = 'WAM-IPE'
            t0 = perf_counter()

            # first, check for file list, create if DNE
            list_file = file_dir + self.modelname +'_list.txt'
            time_file = file_dir + self.modelname +'_times.txt'
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
                patterns = unique([basename(f)[:-19] for f in files]) 
                self.filename = ''.join([f+',' for f in files])[:-1]
                datetime_str = filename_to_dts(basename(files[0])[-18:-3])
                # strings in format = YYYY-MM-DD HH:MM:SS
                self.filedate = datetime.strptime(
                    datetime_str[:10]+' 00:00:00', '%Y-%m-%d %H:%M:%S'
                    ).replace(tzinfo=timezone.utc)

                # establish time attributes from filenames
                for p in patterns:
                    # get list of files to loop through later
                    pattern_files = sorted(glob(file_dir+p+'*.nc'))
                    self.pattern_files[p] = pattern_files
                    self.times[p] = {'start': [], 'end': [], 'all': []}
                    
                    # loop through to get times, one file per time
                    time_str = [basename(f)[-18:-3] for f in pattern_files]
                    times = filename_to_hrs(time_str, self.filedate)
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

            # if variables are given as integers, convert to standard names
            if len(variables_requested) > 0:
                if isinstance(variables_requested[0], int):
                    tmp_var = [value[0] for key, value in
                               model_varnames.items() if value[2] in
                               variables_requested]
                    variables_requested = tmp_var

            # store variables
            self.missing_value = NaN
            self._registered = 0
            self.varfiles = {}  # store which variable came from which file
            self.gvarfiles = {}  # store file variable name similarly
            self.err_list = []
            self.variables = {}

            if printfiles:
                print(f'{len(self.filename)} Files:')
                for file in self.filename:
                    print(file)

            # perform initial check on variables_requested list
            if len(variables_requested) > 0 and fulltime and\
                    variables_requested != 'all':
                test_list = [value[0] for key, value in model_varnames.items()]
                err_list = [item for item in variables_requested if item not in
                            test_list]
                if len(err_list) > 0:
                    print('Variable name(s) not recognized:', err_list)

            # loop through file patterns for coordinate grids + var mapping
            for p in self.pattern_files.keys():
                # use first file for coordinates and variable mapping
                file = self.pattern_files[p][0]
                cdf_data = Dataset(file, 'r')

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
                
                # retrieve coordinate grids, assuming constant in time
                for key in cdf_data.variables.keys():
                    print(key, cdf_data.variables[key].dimensions)
                setattr(self, '_lat_'+p,
                        array(cdf_data.variables['lat']))
                setattr(self, '_lon_'+p,
                        append(array(cdf_data.variables['lon']) - 180., 180.))  # CHECK THIS WITH PLOTS!!!!
                if 'alt' in cdf_data.variables.keys():
                    setattr(self, '_height_'+p,
                            array(cdf_data.variables['alt']))  # km
                    print(p, getattr(self, '_height_'+p).shape)
                elif 'gsm10' in p:  # add missing pressure level for gsm10
                    var_3D = [key for key, value in model_varnames.items()
                              if key in gvar_list and len(value[-2])==4][0]
                    ilev_name = cdf_data.variables[var_3D].dimensions[0]
                    if ilev_name not in cdf_data.variables.keys():
                        len_ilev = cdf_data.variables[var_3D].shape[0]
                        setattr(self, '_ilev_'+p,
                                linspace(0, len_ilev-1, len_ilev))
                    else:
                        setattr(self, '_ilev_'+p,
                                array(cdf_data.variables[ilev_name]))
                
                cdf_data.close()
            print('self.varfiles', self.varfiles)
            print('self.gvarfiles', self.gvarfiles)

            # collect all possible variables in set of files and return
            if not fulltime and variables_requested == 'all':
                self.var_dict, gvar_list = {}, []
                # loop through gvar_list stored for files to make a master list
                for i in range(len(self.patterns)):
                    gvar_list += self.gvarfiles[self.patterns[i]]
                self.var_dict = {value[0]: value[1:] for key, value in
                                 model_varnames.items() if key in gvar_list}
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
                      'files found:\n',
                      self.patterns, self.err_list)

            # register interpolators for each variable
            t_reg = perf_counter()
            # store original list b/c gridded interpolators change keys list
            varname_list = [key for key in self.variables.keys()]
            print(varname_list)
            for varname in varname_list:
                print(varname)
                # determine which file the variable came from, retrieve the coords
                coord_list = [value[-2] for key, value in
                             model_varnames.items() if value[0] == varname][0]
                coord_dict = {}
                for i in range(len(self.patterns)):
                    if varname in self.varfiles[self.patterns[i]]:
                        # get the correct coordinates
                        coord_dict['time'] = {
                            'data': getattr(self, '_time_'+self.patterns[i]),
                            'units': 'hr'}
                        coord_dict['lon'] = {
                            'data': getattr(self, '_lon_'+self.patterns[i]),
                            'units': 'deg'}
                        lon = coord_dict['lon']['data']
                        # need to do fancy lon shift, ignoring appended value
                        lon_le180 = list(where(lon <= 0)[0])  # after shift
                        lon_ge180 = list(where((lon >= 0) & (lon < 180.))[0])
                        coord_dict['lat'] = {
                            'data': getattr(self, '_lat_'+self.patterns[i]),
                            'units': 'deg'}
                        if len(coord_list) == 4:
                            coord_dict['height'] = {'units': 'km'}
                            if hasattr(self, '_height_'+self.patterns[i]):
                                coord_dict['height']['data'] = \
                                    getattr(self, '_height_'+self.patterns[i])
                            else:  # find a height grid of matching dimensions
                                height_size = \
                                    self.variables[varname]['data'][0].shape[0]
                                print('Finding an alternative alt grid of ' +
                                      f'size {height_size}...')
                                for pattern in self.patterns:  # search for it
                                    if hasattr(self, '_height_'+pattern):
                                        h = getattr(self, '_height_'+pattern
                                                    ).size
                                        print(pattern, h)
                                        if h == height_size:
                                            coord_dict['height']['data'] = h
                                if 'data' not in coord_dict['height'].keys():
                                    print(f'Cannot functionalize {varname} ' +
                                          'due to missing vertical coordinat' +
                                          'e grid. Skipping.')
                        break  # move on when coordinate grids are found

                if 'height' in coord_dict.keys() and 'data' not in \
                    coord_dict['height'].keys():
                        continue  # skip variable

                # define the operation to occur when the time slice is loaded
                def func(cdf_data_object):
                    '''define the logic to prepare the data per time slice.'''
                    tmp = array(cdf_data_object).T  # transpose (2D and 3D)
                    return tmp[lon_ge180 + lon_le180]

                # define and register the interpolators
                coord_str = [value[3]+value[4] for key, value in
                             model_varnames.items() if value[0] == varname][0]
                self = RU.time_interp(self, coord_dict, varname,
                                      self.variables[varname], gridded_int,
                                      coord_str, func)
            if verbose:
                print(f'Took {perf_counter()-t_reg:.5f}s to register ' +
                      f'{len(varname_list)} variables.')
            if verbose:
                print(f'Took a total of {perf_counter()-t0:.5f}s to kamodofy' +
                      f' {len(varname_list)} variables.')
    return MODEL
