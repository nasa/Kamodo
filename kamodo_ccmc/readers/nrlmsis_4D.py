'''
NRLMSIS model reader adapted from the IRI reader by Darren De Zeeuw, 2026
'''

from numpy import nan

# variable name in file: [standardized variable name, descriptive term, 
#                         index, coord_system, type, dimensions, units]
model_varnames = {
    'O': ['N_O', 'number density of atomic oxygen', 
          0, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'], '1/cm**3'],
    'N2': ['N_N2', 'number density of molecular nitrogen', 
           1, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'], '1/cm**3'],
    'O2': ['N_O2', 'number density of molecular oxygen', 
           2, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'], '1/cm**3'],
    'MASS': ['rho', 'mass density', 
             3, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'], 'g/cm**3'],
    'NT': ['T_n', 'neutral temperature', 
           4, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'], 'K'],
    'ET': ['T_exo', 'exospheric temperature', 
           5, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'], 'K'],
    'He': ['N_He', 'number density of atomic helium', 
           6, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'], '1/cm**3'],
    'AR': ['N_Ar', 'number density of argon', 
           7, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'], '1/cm**3'],
    'AO': ['N_AO', 'number density of anomalous oxygen', 
           8, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'], '1/cm**3'],
    'H': ['N_H', 'number density of atomic hydrogen', 
          9, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'], '1/cm**3'],
    'N': ['N_N', 'number density of atomic nitrogen', 
          10, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'], '1/cm**3']
}

def MODEL():

    from kamodo import Kamodo
    from os.path import basename
    from numpy import array, transpose, nan, unique
    from numpy import where, append
    from time import perf_counter
    from datetime import datetime, timedelta, timezone
    import kamodo_ccmc.readers.reader_utilities as RU

    from scipy.interpolate import RegularGridInterpolator as rgiND
    from numpy import log, exp, delete, reshape
    
    class MODEL(Kamodo):
        '''NRLMSIS model data reader.'''
        
        def __init__(self, file_dir, variables_requested=[],
                     printfiles=False, filetime=False, gridded_int=True,
                     verbose=False, logD=True, **kwargs):
            super(MODEL, self).__init__(**kwargs)
            self.modelname = 'NRLMSIS'
            self.doublemidnight = False
            self.missinglon = False
            self.logD = logD
            t0 = perf_counter()

            # first, check for file list, create if DNE
            list_file = file_dir + self.modelname + '_list.txt'
            time_file = file_dir + self.modelname + '_times.txt'
            self.times, self.pattern_files = {}, {}
            if not RU._isfile(list_file) or not RU._isfile(time_file):
                # collect filenames
                files = sorted(RU.glob(file_dir+'*.nc'))
                # For NRLMSIS00.3D.2025260.nc, strip the last 11 chars ('.2025260.nc')
                patterns = unique([basename(f)[:-11] for f in files])
                self.filename = ''.join([f+',' for f in files])[:-1]

                # establish time attributes
                for p in patterns:
                    pattern_files = sorted(RU.glob(file_dir+p+'*.nc'))
                    self.pattern_files[p] = pattern_files
                    self.times[p] = {'start': [], 'end': [], 'all': []}

                    # loop through to get times
                    for f in range(len(pattern_files)):
                        cdf_data = RU.Dataset(pattern_files[f], filetype='netCDF3')
                        tmp = array(cdf_data.variables['time'])/60. + \
                            float(f)*24.  # hrs since midnite 1st file
                        
                        self.times[p]['start'].append(tmp[0])
                        self.times[p]['end'].append(tmp[-1])
                        self.times[p]['all'].extend(tmp)
                        cdf_data.close()
                    self.times[p]['start'] = array(self.times[p]['start'])
                    self.times[p]['end'] = array(self.times[p]['end'])
                    self.times[p]['all'] = array(self.times[p]['all'])

                # datetime object for midnight on date (YYYYDDD)
                # files[0][-10:-6] is YYYY, files[0][-6:-3] is DDD
                self.filedate = datetime(int(files[0][-10:-6]), 1, 1, 0, 0, 0
                                         ).replace(tzinfo=timezone.utc) + \
                    timedelta(days=int(files[0][-6:-3]) - 1)

                # create time list file if DNE
                RU.create_timelist(list_file, time_file, self.modelname,
                                   self.times, self.pattern_files,
                                   self.filedate)
            else:  
                self.times, self.pattern_files, self.filedate, self.filename =\
                    RU.read_timelist(time_file, list_file)
            
            for p in self.pattern_files.keys():
                for i in range(len(self.times[p]['all'])-1):
                    if self.times[p]['all'][i] == self.times[p]['all'][i+1]:
                        self.doublemidnight = True
                        self.times[p]['all'][i] += -.0001
                        for j in range(len(self.times[p]['end'])):
                            if self.times[p]['end'][j] == self.times[p]['all'][i+1]:
                                self.times[p]['end'][j] += -.0001
            if filetime:
                return  

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
                cdf_data = RU.Dataset(self.pattern_files[p][0], filetype='netCDF3')
                if len(variables_requested) > 0 and variables_requested != 'all':
                    gvar_list = [key for key in model_varnames.keys()
                                 if key in cdf_data.variables.keys() and
                                 model_varnames[key][0] in variables_requested]
                    if len(gvar_list) != len(variables_requested):
                        err_list = [value[0] for key, value in
                                    model_varnames.items()
                                    if key not in cdf_data.variables.keys() and
                                    value[0] in variables_requested]
                        self.err_list.extend(err_list)  
                else:
                    gvar_list = [key for key in model_varnames.keys()
                                 if key in cdf_data.variables.keys()]
                cdf_data.close()
                
                self.varfiles[p] = [model_varnames[key][0] for key in gvar_list]
                self.gvarfiles[p] = gvar_list

            var_list = []
            for p in self.varfiles.keys():
                var_list.extend(self.varfiles[p])
            err_list = [var for var in self.err_list if var not in var_list]
            if len(err_list) > 0:
                print('Some requested variables are not available: ', err_list)

            if variables_requested == 'all':
                self.var_dict = {value[0]: value[1:] for key, value in
                                 model_varnames.items() if value[0] in
                                 var_list}
                return

            self.variables = {}
            for p in self.pattern_files.keys():
                cdf_data = RU.Dataset(self.pattern_files[p][0], filetype='netCDF3')
                lon = array(cdf_data.variables['lon'][:], copy=True)
                lon_le180 = list(where(lon <= 180)[0])  
                lon_ge180 = list(where((lon >= 180) & (lon < 360.))[0])
                tmp = lon - 180.  
                setattr(self, '_lon_'+p, tmp)
                setattr(self, '_lon_idx_'+p, lon_ge180+lon_le180)
                setattr(self, '_lat_'+p,
                        array(cdf_data.variables['lat'][:], copy=True))
                # Map NRLMSIS 'ht' to kamodo 'height'
                if '3D' in p:
                    setattr(self, '_height_'+p,
                            array(cdf_data.variables['ht'][:], copy=True))
                cdf_data.close()
                for var in self.gvarfiles[p]:
                    self.variables[model_varnames[var][0]] = {
                        'units': model_varnames[var][-1], 'data': p}

            self.missing_value = nan
            self._registered = 0
            if verbose:
                print(f'Took {perf_counter()-t0:.6f}s to read in data')
            if printfiles:
                print(self.filename)

            t_reg = perf_counter()
            varname_list = list(self.variables.keys())
            for varname in varname_list:
                self.register_variable(varname, gridded_int)
            if verbose:
                print(f'Took {perf_counter()-t_reg:.5f}s to register {len(varname_list)} variables.')
                print(f'Took a total of {perf_counter()-t0:.5f}s to kamodofy {len(varname_list)} variables.')

        def register_variable(self, varname, gridded_int):
            key = self.variables[varname]['data']
            coord_list = [value[5] for key, value in model_varnames.items()
                          if value[0] == varname][0]
            coord_dict = {'time': {'units': 'hr',
                                   'data': self.times[key]['all']}}
            if 'lat' in coord_list:   
                coord_dict['lon'] = {'units': 'deg', 'data':
                                     getattr(self, '_lon_'+key)}
                step = coord_dict['lon']['data'][1] - coord_dict['lon']['data'][0]
                if coord_dict['lon']['data'][-1] < 179.5 and \
                   coord_dict['lon']['data'][0] == -180. and \
                   coord_dict['lon']['data'][-1]+step == 180.:
                    self.missinglon = True
                    coord_dict['lon']['data'] = append(coord_dict['lon']['data'], [180.], axis=0)
                coord_dict['lat'] = {'units': 'deg', 'data':
                                     getattr(self, '_lat_'+key)}
            if 'height' in coord_list:
                coord_dict['height'] = {'units': 'km', 'data':
                                        getattr(self, '_height_'+key)}
            lon_idx = getattr(self, '_lon_idx_'+key)
            gvar = [key for key, value in model_varnames.items() if
                    value[0] == varname][0]  
            coord_str = [value[3]+value[4] for key, value in
                         model_varnames.items() if value[0] == varname][0]

            def func_custom(i):
                file = self.pattern_files[key][i]
                cdf_data = RU.Dataset(file, filetype='netCDF3')
                data = array(cdf_data.variables[gvar])
                if hasattr(cdf_data.variables[gvar][0], 'fill_value'):
                    fill_value = cdf_data.variables[gvar][0].fill_value
                else:
                    fill_value = None
                cdf_data.close()
                time_index_start = self.times[key]['start_index'][i]
                time_index_end = self.times[key]['start_index'][i+1]+1
                if self.doublemidnight:
                    time_index_end += 1
                    if i > 0:
                        time_index_start += 1
                if file != self.pattern_files[key][-1]:  
                    next_file = self.pattern_files[key][i+1]
                    cdf_data = RU.Dataset(next_file, filetype='netCDF3')
                    data_slice = array(cdf_data.variables[gvar][0])
                    cdf_data.close()
                    data = append(data, [data_slice], axis=0)
                coord_dict_data = [ coord_dict[key]['data'] for key in coord_dict ]
                times_file = self.times[key]['all'][time_index_start:time_index_end]
                coord_dict_data[0] = array(times_file)

                if fill_value is not None:  
                    data = where(data != fill_value, data, nan)
                if len(data.shape) == 3:
                    log_variable = log(transpose(data, (0, 2, 1)))
                    if self.missinglon:
                        lvs = log_variable[:,0,:]
                        lvs2 = reshape(lvs, (lvs.shape[0], 1, lvs.shape[1]))
                        log_variable = append(log_variable, lvs2, axis=1)
                elif len(data.shape) == 4:
                    # NRLMSIS order: time, ht, lat, lon -> Kamodo expected: time, lon, lat, height
                    log_variable = log(transpose(data, (0, 3, 2, 1)))
                    if self.missinglon:
                        lvs = log_variable[:,0,:,:]
                        lvs2 = reshape(lvs, (lvs.shape[0], 1, lvs.shape[1], lvs.shape[2]))
                        log_variable = append(log_variable, lvs2, axis=1)
                rgi = rgiND(coord_dict_data, log_variable[:, lon_idx],
                            bounds_error=False,fill_value=nan)
                def interp4d_custom(xvec):
                    return exp(rgi(xvec))
                return interp4d_custom
            
            def func(i):
                file = self.pattern_files[key][i]
                cdf_data = RU.Dataset(file, filetype='netCDF3')
                data = array(cdf_data.variables[gvar][:], copy=True)
                if hasattr(cdf_data.variables[gvar][0], 'fill_value'):
                    fill_value = cdf_data.variables[gvar][0].fill_value
                else:
                    fill_value = None
                cdf_data.close()
                if file != self.pattern_files[key][-1]:  
                    next_file = self.pattern_files[key][i+1]
                    cdf_data = RU.Dataset(next_file, filetype='netCDF3')
                    data_slice = array(cdf_data.variables[gvar][0])
                    cdf_data.close()
                    data = append(data, [data_slice], axis=0)
                if fill_value is not None:  
                    data = where(data != fill_value, data, nan)
                if len(data.shape) == 3:
                    variable = transpose(data, (0, 2, 1))
                elif len(data.shape) == 4:
                    # NRLMSIS order: time, ht, lat, lon -> time, lon, lat, height
                    variable = transpose(data, (0, 3, 2, 1))
                return variable[:, lon_idx]
            
            if (varname[0:3] == "rho" or varname[0:2] == "N_") and self.logD:
                self = RU.Functionalize_Dataset(
                    self, coord_dict, varname, self.variables[varname],
                    gridded_int, coord_str, interp_flag=2, func=func_custom,
                    func_default='custom',
                    times_dict=self.times[key])
            else:
                self = RU.Functionalize_Dataset(
                    self, coord_dict, varname, self.variables[varname],
                    gridded_int, coord_str, interp_flag=2, func=func,
                    times_dict=self.times[key])                
            return
            
    return MODEL

