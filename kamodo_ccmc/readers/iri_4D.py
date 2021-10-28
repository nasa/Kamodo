'''
Written by Rebecca Ringuette, 2021
'''
from datetime import datetime, timedelta, timezone


#variable name in file: [standardized variable name, descriptive term, units]
model_varnames = {'Ne':['N_e','variable description',0,'SPH','sph',['time','lon','lat','radius'],'1/m**3'], 
                'Te':['T_e','variable description',1,'SPH','sph',['time','lon','lat','radius'],'K'],
                'Ti':['T_i','variable description',2,'SPH','sph',['time','lon','lat','radius'],'K'], 
                'Tn':['T_n','variable description',3,'SPH','sph',['time','lon','lat','radius'],'K'],
                'O+':['N_Oplus','variable description',4,'SPH','sph',['time','lon','lat','radius'],'1/m**3'],
                'H+':['N_Hplus','variable description',5,'SPH','sph',['time','lon','lat','radius'],'1/m**3'],
                'He+':['N_Heplus','variable description',6,'SPH','sph',['time','lon','lat','radius'],'1/m**3'],
                'O2+':['N_O2plus','variable description',7,'SPH','sph',['time','lon','lat','radius'],'1/m**3'],
                'NO+':['N_NOplus','variable description',8,'SPH','sph',['time','lon','lat','radius'],'1/m**3'],
                'N+':['N_Nplus','variable description',9,'SPH','sph',['time','lon','lat','radius'],'1/m**3'],
                'TEC':['TEC','variable description',10,'SPH','sph',['time','lon','lat'],'10**16/m**2'],
                'NmF2':['NmF2','variable description',11,'SPH','sph',['time','lon','lat'],'1/m**3'],
                'HmF2':['HmF2','variable description',12,'SPH','sph',['time','lon','lat'],'km']}

def ts_to_hrs(time_val, filedate):
    '''Convert utc timestamp to hours since midnight on filedate.'''
    
    return (datetime.utcfromtimestamp(time_val).replace(tzinfo=timezone.utc)-filedate).total_seconds()/3600.

#times from file converted to seconds since midnight of filedate
#plotting input times will be datetime strings of format 'YYYY-MM-DD HH:mm:ss'
#filedate is self.filedate from iri object
#converts to hours since midnight of filedate for plotting
def MODEL():
    from kamodo import Kamodo
    #print('KAMODO IMPORTED!')
    from netCDF4 import Dataset
    from os.path import basename
    from numpy import array, transpose, NaN, unique, append, zeros, abs, diff, where
    from time import perf_counter
    from astropy.constants import R_earth
    from kamodo_ccmc.readers.reader_utilities import regdef_4D_interpolators, regdef_3D_interpolators

    class MODEL(Kamodo):
        '''IRI model data reader.'''
        def __init__(self, full_filename3d, variables_requested = [], runname = "noname",
                     printfiles=False, filetime=False, gridded_int=True, fulltime=True,
                     verbose=False,**kwargs): #                 time_index=None, time_seconds=None,
            # Prepare model for function registration for the input argument
            super(MODEL, self).__init__(**kwargs)
            t0 = perf_counter()
    
            #collect filenames
            filename = basename(full_filename3d)
            file_dir = full_filename3d.split(filename)[0]            
            if '.2D.' in filename:  #require that input filename be for 3D file
                full_filename2d = full_filename3d
                f = full_filename3d.replace('.2D.', '.3D.')  #can't replace in place
                full_filename3d = f
            else: 
                full_filename2d = full_filename3d.replace('.3D.','.2D.')
            self.filename = full_filename3d+','+full_filename2d
            
            #establish time attributes first
            iri3D = Dataset(full_filename3d, 'r')
            time = array(iri3D.variables['time'])/60.  #convert to hours since midnight of file        
            self.filedate = datetime(int(filename[-10:-6]),1,1,0,0,0).replace(tzinfo=timezone.utc)+\
                timedelta(days=int(filename[-6:-3])-1)
            #strings with timezone info chopped off (UTC anyway)
            self.datetimes=[(self.filedate+timedelta(hours=time[0])).isoformat(sep=' ')[:19], 
                            (self.filedate+timedelta(hours=time[-1])).isoformat(sep=' ')[:19]]  #strings
            self.filetimes=[datetime.timestamp(datetime.strptime(dt, '%Y-%m-%d %H:%M:%S').replace(\
                tzinfo=timezone.utc)) for dt in self.datetimes]   #timestamp in seconds, for value matching in wrapper
            self.dt = diff(time).max()*3600.  #time is in hours since midnight
                
            if filetime and not fulltime: #(used when searching for neighboring files below)
                return  #return times as is to prevent recursion
            
            #if variables are given as integers, convert to standard names
            if len(variables_requested)>0:
                if isinstance(variables_requested[0], int):
                    tmp_var = [value[0] for key, value in model_varnames.items()\
                                           if value[2] in variables_requested]
                    variables_requested = tmp_var        
            
            if fulltime:  #add boundary time (default value)
                #find other files with same pattern
                from glob import glob
                
                file_pattern = file_dir+'IRI.3D.*.nc' #returns a string for iri
                files = sorted(glob(file_pattern))
                filenames = unique([basename(f) for f in files])
                
                #find closest file by utc timestamp
                #iri has an open time at the end, so need a beginning time from the closest file
                #files are automatically sorted by YYMMDD, so next file is next in the list
                current_idx = where(filenames==filename)[0]
                if current_idx+1==len(files):
                    print('No later file available.')
                    filecheck = False  
                    if filetime:
                        return   
                else:
                    min_file = file_dir+filenames[current_idx+1][0]  #+1 for adding an end time
                    kamodo_test = MODEL(min_file, filetime=True, fulltime=False)               
                    time_test = abs(kamodo_test.filetimes[0]-self.filetimes[1])  #never empty b/c includes itself 
                    if time_test<=self.dt:  #if nearest file time at least within one timestep (hrs)                
                        filecheck = True
                    
                        #time only version if returning time for searching
                        if filetime:
                            kamodo_neighbor = MODEL(min_file, fulltime=False, filetime=True)
                            self.datetimes[1] = kamodo_neighbor.datetimes[0]
                            self.filetimes[1] = kamodo_neighbor.filetimes[0]
                            return  #return object with additional time (for SF code) 
                        
                        #get kamodo object with same requested variables to add to each array below
                        if verbose: print(f'Took {perf_counter()-t0:.3f}s to find closest file.')
                        kamodo_neighbor = MODEL(min_file, 
                                                variables_requested=variables_requested, 
                                                fulltime=False)
                        self.datetimes[1] = kamodo_neighbor.datetimes[0]
                        self.filetimes[1] = kamodo_neighbor.filetimes[0]
                        short_data = kamodo_neighbor.short_data                          
                        if verbose: print(f'Took {perf_counter()-t0:.3f}s to get data from closest file.')
                    else:
                        print(f'No later file found within {diff(time).max()*3600.:.1f}s.')
                        filecheck = False        
                        if filetime:
                            return                    

            #perform initial check on variables_requested list
            if len(variables_requested)>0 and fulltime:
                test_list = [value[0] for key, value in model_varnames.items()]
                err_list = [item for item in variables_requested if item not in test_list]
                if len(err_list)>0: print('Variable name(s) not recognized:', err_list)
                
            #collect variable list
            iri2D = Dataset(full_filename2d, 'r')
            if len(variables_requested)>0:
                gvar_list_2d = [key for key, value in model_varnames.items() \
                                 if value[0] in variables_requested and \
                                     key in iri2D.variables.keys()]  # file variable names
                gvar_list_3d = [key for key, value in model_varnames.items() \
                                 if value[0] in variables_requested and \
                                     key in iri3D.variables.keys()]  # file variable names                
                    
                #check for variables requested but not available
                if len(gvar_list_2d)+len(gvar_list_3d)!=len(variables_requested):
                    err_list = [value[0] for key, value in model_varnames.items() \
                                 if value[0] in variables_requested and \
                                     key not in gvar_list_2d and key not in gvar_list_3d]
                    if len(err_list)>0: print('Some requested variables are not available:', err_list)
            else:  #only input variables on the avoid_list if specifically requested
                gvar_list_2d = [key for key in iri2D.variables.keys() if key in model_varnames.keys()]
                gvar_list_3d = [key for key in iri3D.variables.keys() if key in model_varnames.keys()]

            #store data for each variable desired
            variables_2d = {model_varnames[var][0]: {'units': model_varnames[var][-1], 
                                                       'data': array(iri2D.variables[var])} \
                                                       for var in gvar_list_2d}
            variables_3d = {model_varnames[var][0]: {'units': model_varnames[var][-1], 
                                                       'data': array(iri3D.variables[var])} \
                                                       for var in gvar_list_3d}            
            variables = variables_3d
            for key in variables_2d: variables[key] = variables_2d[key]
                
            #prepare and return data only for first timestamp
            if not fulltime:  
                iri3D.close()
                iri2D.close()
                variables['time'] = self.filetimes[0]
                self.short_data = variables
                return        
    
            #### Store coordinate data as class attributes   
            if filecheck:
                new_time = ts_to_hrs(short_data['time'], self.filedate)  #new time in hours since midnight
                self._time = append(time, new_time) #append new time in hours since midnight
            else: 
                self._time = time           
                
            #collect data and make dimensional grid from 3D file
            self._lon = array(iri3D.variables['lon'])
            self._lat = array(iri3D.variables['lat'])
            #convert height in km to radius in R_E to conform to SPH coord sys in spacepy
            self._radius = (array(iri3D.variables['ht'])+R_earth.value/1000.)/(R_earth.value/1000.)
            iri3D.close()   #close netCDF4 files
            iri2D.close()
            
            #store a few items in iri object
            self.missing_value = NaN
            self._registered = 0
            self.runname=runname
            self.modelname = 'IRI'
            if verbose: print(f'Took {perf_counter()-t0:.6f}s to read in data')
            if printfiles: print(self.filename)
            
            # register interpolators for each requested variable
            t_reg = perf_counter()
            varname_list, self.variables = [key for key in variables.keys()], {}  #store original list b/c gridded interpolators        
            for varname in varname_list:
                if len(variables[varname]['data'].shape)==3:
                    if filecheck:  #if neighbor found
                        #append data for last time stamp, transpose and register
                        data_shape = list(variables[varname]['data'].shape)
                        data_shape[0]+=1  #add space for time
                        new_data = zeros(data_shape)                
                        new_data[:-1,:,:] = variables[varname]['data']  #put in current data
                        new_data[-1,:,:] = short_data[varname]['data'][0,:,:]  #add in data for additional time
                        variable = transpose(new_data, (0,2,1)) #(t,lat,lon) -> (t,lon,lat)
                    else:
                        variable = transpose(variables[varname]['data'], (0,2,1)) #(t,lat,lon) -> (t,lon,lat)
                    self.variables[varname] = dict(units = variables[varname]['units'], data = variable)
                    self.register_3D_variable(self.variables[varname]['units'], 
                                          self.variables[varname]['data'], varname,
                                          gridded_int)
                elif len(variables[varname]['data'].shape)==4:
                    if filecheck:
                        #append data for last time stamp, transpose and register
                        data_shape = list(variables[varname]['data'].shape)
                        data_shape[0]+=1  #add space for time
                        new_data = zeros(data_shape)
                        new_data[:-1,:,:,:] = variables[varname]['data']  #put in current data
                        new_data[-1,:,:,:] = short_data[varname]['data'][0,:,:,:]   #add in data for additional time                
                        variable = transpose(new_data, (0,3,2,1)) #(t,h,lat,lon) -> (t,lon,lat,h)
                    else:
                        variable = transpose(variables[varname]['data'], (0,3,2,1)) #(t,h,lat,lon) -> (t,lon,lat,h)
                    self.variables[varname] = dict(units = variables[varname]['units'], data = variable)
                    self.register_4D_variable(self.variables[varname]['units'], 
                                          self.variables[varname]['data'], varname,
                                          gridded_int)
            if verbose: print(f'Took {perf_counter()-t_reg:.5f}s to register '+\
                              f'{len(varname_list)} variables.')
            if verbose: print(f'Took a total of {perf_counter()-t0:.5f}s to kamodofy '+\
                              f'{len(varname_list)} variables.')            
                
    
        #define and register a 3D variable
        def register_3D_variable(self, units, variable, varname, gridded_int):
            """Registers a 3d interpolator with 3d signature"""
            
            #define and register the interpolators
            xvec_dependencies = {'time':'hr','lon':'deg','lat':'deg'}
            self = regdef_3D_interpolators(self, units, variable, self._time, 
                                           self._lon, self._lat, varname, 
                                           xvec_dependencies, gridded_int)
            return 
        
        #define and register a 4D variable
        def register_4D_variable(self, units, variable, varname, gridded_int):
            """Registers a 4d interpolator with 4d signature"""
            
            #define and register the fast interpolator
            xvec_dependencies = {'time':'hr','lon':'deg','lat':'deg','radius':'R_E'}
            self = regdef_4D_interpolators(self, units, variable, self._time,
                                              self._lon, self._lat, self._radius,
                                              varname, xvec_dependencies, gridded_int)
            return
    return MODEL