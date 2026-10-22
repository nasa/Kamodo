'''
HASDM model reader for Kamodo
'''

from numpy import nan

# variable name in file: [standardized variable name, descriptive term, 
#                         index, coord_system, type, dimensions, units]
model_varnames = {
    'rho': ['rho', 'mass density', 
             0, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'], 'kg/m**3']
}

def MODEL():

    from kamodo import Kamodo
    from os.path import basename
    from numpy import array, transpose, nan, unique
    from numpy import where, append
    from time import perf_counter
    from datetime import datetime, timezone
    import kamodo_ccmc.readers.reader_utilities as RU

    from scipy.interpolate import RegularGridInterpolator as rgiND
    from numpy import log, exp
    
    class MODEL(Kamodo):
        '''HASDM model data reader.'''
        
        def __init__(self, file_dir, variables_requested=[],
                     printfiles=False, filetime=False, gridded_int=True,
                     verbose=False, logD=True, **kwargs):
            super(MODEL, self).__init__(**kwargs)
            self.modelname = 'HASDM'
            self.doublemidnight = False
            self.logD = logD
            t0 = perf_counter()

            list_file = file_dir + self.modelname + '_list.txt'
            time_file = file_dir + self.modelname + '_times.txt'
            self.times, self.pattern_files = {}, {}
            
            if not RU._isfile(list_file) or not RU._isfile(time_file):
                files = sorted(RU.glob(file_dir+'*.nc'))
                # Filename format: HASDM_YYYY.nc
                # Strip '_YYYY.nc' (8 characters) to group all files under the 'HASDM' pattern
                patterns = unique([basename(f)[:-8] for f in files])
                self.filename = ''.join([f+',' for f in files])[:-1]

                for p in patterns:
                    pattern_files = sorted(RU.glob(file_dir+p+'*.nc'))
                    self.pattern_files[p] = pattern_files
                    self.times[p] = {'start': [], 'end': [], 'all': []}

                    # Extract year from the first filename (e.g., HASDM_2019.nc -> 2019)
                    year = int(basename(pattern_files[0])[-7:-3])
                    if not hasattr(self, 'filedate'):
                        self.filedate = datetime(year, 1, 1, 0, 0, 0, tzinfo=timezone.utc)
                    
                    filedate_timestamp = self.filedate.timestamp()

                    for f in range(len(pattern_files)):
                        cdf_data = RU.Dataset(pattern_files[f], filetype='netCDF4')
                        # Cast Unix timestamps to float64 immediately to preserve precision
                        raw_times = array(cdf_data.variables['time'][:], dtype='float64')
                        # Convert to hours since self.filedate
                        tmp = (raw_times - filedate_timestamp) / 3600.0
                        
                        self.times[p]['start'].append(tmp[0])
                        self.times[p]['end'].append(tmp[-1])
                        self.times[p]['all'].extend(tmp)
                        cdf_data.close()
                        
                    self.times[p]['start'] = array(self.times[p]['start'])
                    self.times[p]['end'] = array(self.times[p]['end'])
                    self.times[p]['all'] = array(self.times[p]['all'])

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
                err_list = [item for item in variables_requested if item not in test_list]
                if len(err_list) > 0:
                    print('Variable name(s) not recognized:', err_list)
                for item in err_list:
                    variables_requested.remove(item)
                if len(variables_requested) == 0:
                    return

            self.gvarfiles, self.varfiles, self.err_list = {}, {}, []
            for p in self.pattern_files.keys():
                cdf_data = RU.Dataset(self.pattern_files[p][0], filetype='netCDF4')
                if len(variables_requested) > 0 and variables_requested != 'all':
                    gvar_list = [key for key in model_varnames.keys()
                                 if key in cdf_data.variables.keys() and
                                 model_varnames[key][0] in variables_requested]
                    if len(gvar_list) != len(variables_requested):
                        err_list = [value[0] for key, value in model_varnames.items()
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
                                 model_varnames.items() if value[0] in var_list}
                return

            self.variables = {}
            for p in self.pattern_files.keys():
                cdf_data = RU.Dataset(self.pattern_files[p][0], filetype='netCDF4')
                
                # HASDM Lons are [0, 15, ..., 345]
                lon = array(cdf_data.variables['lons'][:], copy=True)
                lon_ge180 = list(where(lon >= 180)[0])
                lon_lt180 = list(where(lon < 180)[0])
                
                # Shift longitude values to [-180, ..., 165]
                tmp = append(lon[lon_ge180] - 360., lon[lon_lt180])
                
                # Pad to +180 to guarantee seamless wrapping
                if tmp[-1] < 180.:
                    tmp = append(tmp, [180.], axis=0)
                    # The +180 data is identical to -180 data, which is at the first ge180 index
                    lon_idx = lon_ge180 + lon_lt180 + [lon_ge180[0]]
                else:
                    lon_idx = lon_ge180 + lon_lt180
                
                setattr(self, '_lon_'+p, tmp)
                setattr(self, '_lon_idx_'+p, lon_idx)
                
                setattr(self, '_lat_'+p, array(cdf_data.variables['lats'][:], copy=True))
                setattr(self, '_height_'+p, array(cdf_data.variables['alts'][:], copy=True))
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

        def register_variable(self, varname, gridded_int):
            key = self.variables[varname]['data']
            coord_dict = {'time': {'units': 'hr', 'data': self.times[key]['all']},
                          'lon': {'units': 'deg', 'data': getattr(self, '_lon_'+key)},
                          'lat': {'units': 'deg', 'data': getattr(self, '_lat_'+key)},
                          'height': {'units': 'km', 'data': getattr(self, '_height_'+key)}}
            
            lon_idx = getattr(self, '_lon_idx_'+key)
            gvar = [key for key, value in model_varnames.items() if value[0] == varname][0]  
            coord_str = [value[3]+value[4] for key, value in model_varnames.items() if value[0] == varname][0]

            def func_custom(i):
                file = self.pattern_files[key][i]
                cdf_data = RU.Dataset(file, filetype='netCDF4')
                data = array(cdf_data.variables[gvar])
                fill_value = getattr(cdf_data.variables[gvar], 'fill_value', None)
                cdf_data.close()
                
                time_index_start = self.times[key]['start_index'][i]
                time_index_end = self.times[key]['start_index'][i+1]+1
                if self.doublemidnight:
                    time_index_end += 1
                    if i > 0:
                        time_index_start += 1
                        
                if file != self.pattern_files[key][-1]:  
                    next_file = self.pattern_files[key][i+1]
                    cdf_data = RU.Dataset(next_file, filetype='netCDF4')
                    data = append(data, [array(cdf_data.variables[gvar][0])], axis=0)
                    cdf_data.close()
                    
                coord_dict_data = [coord_dict[k]['data'] for k in coord_dict]
                coord_dict_data[0] = array(self.times[key]['all'][time_index_start:time_index_end])

                if fill_value is not None:  
                    data = where(data != fill_value, data, nan)
                
                # HASDM dims: (time, alt, lat, lon) -> Kamodo dims: (time, lon, lat, alt)
                log_variable = log(transpose(data, (0, 3, 2, 1)))
                
                rgi = rgiND(coord_dict_data, log_variable[:, lon_idx],
                            bounds_error=False, fill_value=nan)
                
                def interp4d_custom(xvec):
                    return exp(rgi(xvec))
                return interp4d_custom
            
            def func(i):
                file = self.pattern_files[key][i]
                cdf_data = RU.Dataset(file, filetype='netCDF4')
                data = array(cdf_data.variables[gvar][:], copy=True)
                fill_value = getattr(cdf_data.variables[gvar], 'fill_value', None)
                cdf_data.close()
                
                if file != self.pattern_files[key][-1]:  
                    next_file = self.pattern_files[key][i+1]
                    cdf_data = RU.Dataset(next_file, filetype='netCDF4')
                    data = append(data, [array(cdf_data.variables[gvar][0])], axis=0)
                    cdf_data.close()
                    
                if fill_value is not None:  
                    data = where(data != fill_value, data, nan)
                
                variable = transpose(data, (0, 3, 2, 1))
                return variable[:, lon_idx]
            
            if (varname[0:3] == "rho" or varname[0:2] == "N_") and self.logD:
                self = RU.Functionalize_Dataset(
                    self, coord_dict, varname, self.variables[varname],
                    gridded_int, coord_str, interp_flag=2, func=func_custom,
                    func_default='custom', times_dict=self.times[key])
            else:
                self = RU.Functionalize_Dataset(
                    self, coord_dict, varname, self.variables[varname],
                    gridded_int, coord_str, interp_flag=2, func=func,
                    times_dict=self.times[key])                
            return
            
    return MODEL

