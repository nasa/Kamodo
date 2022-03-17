'''
Written by Rebecca Ringuette, 2022

file_dir = 'C:/Users/rringuet/Kamodo_WinDev1/AmGEO/'
file name: '*_YYYYMMDDx.h5' where x is either N for northern hemisphere or S for southern hemisphere
data structure: file object is a nested dictionary, keys are format 'YYYYMMDD_HHMMSSx' and 'lats' 'lons'
    - subkeys are dataset objects all of the same shape, one per timestep
    - 'lats' has same shape as timestep datasets, all [:,x] are same list of values, 
        all [x,:] are a single repeated value -> constant latitude grid for northern hemisphere
        (24 gridpoints in latitude, ranges from about 50 to 88.3 degrees in steps of 1.67 deg)
        ----NEED LATITUDE WRAPPING----
    - 'lons' has same shape as timestep datasets, all [:,x] are a single repeated value,
        all [x,:] are the same list of values -> constant longitude grid for northern hemisphere
        (37 gridpoints in longitude, from 0 to 360 degrees in steps of 10 deg)
        ----NEED LONGITUDE SWAPPING----
    - timestep datasets have shape (lat,lon), except int_joule_heat is a time series
    - data type is float

- the values reported in the southern hemisphere should be taken with the given sign 
    and the values reported in the northern hemisphere should be multiplied by negative one.
    (variables affected: ‘Electric Field (equatorward)’, ‘Spacecraft-Observed Magnetic Perturbations (equatorward)’ 
     and ‘Ion Drift Velocity (equatorward)’)
- solar wind speed, IMF By and IMF Bz variables are in attributes. Should be included in output.
- time res is 5 min, beg/end times are 2.5 min from midnight, so choosing to add time at end
'''

from datetime import datetime, timezone
from numpy import vectorize

#'C:/Users/rringuet/Kamodo_WinDev1/AmGEO/'
model_varnames = {'E_ph':['E_east','Electric Field (eastward)',0,'MAG','sph',['time','lon','lat'],'V/m'], 
                  'E_th':['E_north','Electric Field (equatorward)',1,'MAG','sph',['time','lon','lat'],'V/m'], 
                  'cond_hall':['Sigma_H','Ovation Pyme Hall Conductance',2,'MAG','sph',['time','lon','lat'],'S'], #units=mho=S
                  'cond_ped':['Sigma_P','Ovation Pyme Pedersen Conductance',3,'MAG','sph',['time','lon','lat'],'S'], #units=mho=S
                  'epot':['V','Electric Potential',4,'MAG','sph',['time','lon','lat'],'V'], 
                  'int_joule_heat_n':['W_JouleN','Northern Hemisphere Integrated Joule Heating',5,'MAG','sph',['time'],'GW'],   #varname???
                  'int_joule_heat_s':['W_JouleS','Southern Hemisphere Integrated Joule Heating',6,'MAG','sph',['time'],'GW'],   #varname???
                  # 1D time series only, but different for each hemisphere
                  'jfac':['j_fac','Field Aligned Current',7,'MAG','sph',['time','lon','lat'],'muA/m**2'], 
                  'joule_heat':['Q_Joule','Joule Heating (E-field^2*Pedersen)',8,'MAG','sph',['time','lon','lat'],'mW/m**2'], 
                  'mpot':['psi','Magnetic Potential',9,'MAG','sph',['time','lon','lat'],'cT/m'], 
                  'sdB_ph':['dB_east','Spacecraft-Observed Magnetic Perturbations (eastward)',10,'MAG','sph',['time','lon','lat'],'nT'], 
                  'sdB_th':['dB_north','Spacecraft-Observed Magnetic Perturbations (equatorward)',11,'MAG','sph',['time','lon','lat'],'nT'], 
                  'v_ph':['v_ieast','Ion Drift Velocity (eastward)',12,'MAG','sph',['time','lon','lat'],'m/s'],
                  'v_th':['v_inorth','Ion Drift Velocity (equatorward)',13,'MAG','sph',['time','lon','lat'],'m/s'],
                  'imf_By':['B_y','measured y component of IMF magnetic field from OMNI',14,'MAG','sph',['time'],'nT'],  #in attrs of time dataset
                  'imf_Bz':['B_z','measured z component of IMF magnetic field from OMNI',15,'MAG','sph',['time'],'nT'], #in attrs of time dataset
                  'solar_wind_speed':['v_sw','measured solar wind speed from OMNI',16,'MAG','sph',['time'],'km/s']  #in attrs of time dataset
                  }

@vectorize
def timestr_hrs(time_str, filedate):
    '''Converts time string from data file into hours since midnight UTC.'''
    
    dt = datetime.strptime(time_str, '%Y%m%d_%H%M%S').replace(tzinfo=timezone.utc)
    return (dt-filedate).total_seconds()/3600.

def timestr_datetime(time_str):
    '''Converts time string into a datetime object at midnight.'''
    
    return datetime.strptime(time_str[:8], '%Y%m%d').replace(tzinfo=timezone.utc)

def timestr_utcts(time_str):
    '''Converts time string into utc timestamp.'''
    
    return datetime.strptime(time_str, '%Y%m%d_%H%M%S').replace(tzinfo=timezone.utc).timestamp()

def timestr_datetimestr(time_str):
    '''Converts time string into standard format (‘YYYY-MM-DD HH:MM:SS’).'''
    
    dt = datetime.strptime(time_str, '%Y%m%d_%H%M%S').replace(tzinfo=timezone.utc)
    return dt.isoformat(sep=' ')[:19]

def ts_to_hrs(time_val, filedate):
    '''Convert utc timestamp to hours since midnight on filedate.'''
    
    return (datetime.utcfromtimestamp(time_val).replace(tzinfo=timezone.utc)-filedate).total_seconds()/3600.

def MODEL():
    from kamodo import Kamodo
    import h5py
    from os.path import basename
    from numpy import array, NaN, unique, append, zeros, abs, diff, sin, cos
    from numpy import where, flip, concatenate, insert, mean, broadcast_to
    from numpy import pi as nppi
    from time import perf_counter
    from astropy.constants import R_earth
    from kamodo_ccmc.readers.reader_utilities import regdef_1D_interpolators, regdef_3D_interpolators

    class MODEL(Kamodo):
        '''AmGEO model data reader.'''
        def __init__(self, full_filenameN, variables_requested = [], 
                     printfiles=False, filetime=False, gridded_int=True, fulltime=True,
                     verbose=False,**kwargs):                
            super(MODEL, self).__init__(**kwargs)
            self.modelname = 'AmGEO'
            t0 = perf_counter()
    
            #collect filenames
            filename = basename(full_filenameN)
            file_dir = full_filenameN.split(filename)[0]            
            if 'S.h5' in filename:  #require that input filename be for N file
                full_filenameS = full_filenameN
                f = full_filenameN.replace('S.h5', 'N.h5')  #can't replace in place
                full_filenameN = f
            else: 
                full_filenameS = full_filenameN.replace('N.h5','S.h5')
            self.filename = full_filenameN+','+full_filenameS
            
            #establish time attributes first
            f_north = h5py.File(full_filenameN, 'r')
            time_list = [key[:-1] for key in f_north.keys() if key not in ['lats','lons']]  #remove 'N' or 'S' at end of time
            self.filedate = timestr_datetime(time_list[0])  #datetime object for midnight on date
            time = timestr_hrs(time_list, self.filedate)  #convert to hours since midnight of file  
            self.datetimes=[timestr_datetimestr(time_list[0]),timestr_datetimestr(time_list[-1])]
            self.filetimes=[timestr_utcts(time_list[0]),timestr_utcts(time_list[-1])]
            self.dt = diff(time).max()*3600.  #convert time resolution to seconds
            
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
                
                file_pattern = file_dir+'*N.h5' #returns a string for amgeo
                files = sorted(glob(file_pattern))  #method may change for AWS
                filenames = unique([basename(f) for f in files])
                
                #find closest file by utc timestamp
                #iri has an open time at the end, so need a beginning time from the closest file
                #files are automatically sorted by YYMMDD, so next file is next in the list
                current_idx = where(filenames==filename)[0]
                if current_idx+1==len(files):
                    if verbose: print('No later file available.')
                    filecheck = False  
                    if filetime:
                        return   
                else:
                    min_file = file_dir+filenames[current_idx+1][0]  #+1 for adding an end time
                    kamodo_test = MODEL(min_file, filetime=True, fulltime=False)               
                    time_test = abs(kamodo_test.filetimes[0]-self.filetimes[1])  #never empty b/c includes itself 
                    if time_test<=self.dt:  #if nearest file time at least within one timestep (hrs)                
                        filecheck = True
                        self.datetimes[1] = kamodo_test.datetimes[0]
                        self.filetimes[1] = kamodo_test.filetimes[0]
                            
                        #time only version if returning time for searching
                        if filetime:
                            return  #return object with additional time (for SF code) 
                        
                        #get kamodo object with same requested variables to add to each array below
                        if verbose: print(f'Took {perf_counter()-t0:.3f}s to find closest file.')
                        kamodo_neighbor = MODEL(min_file, 
                                                variables_requested=variables_requested, 
                                                fulltime=False)
                        short_data = kamodo_neighbor.short_data                          
                        if verbose: print(f'Took {perf_counter()-t0:.3f}s to get data from closest file.')
                    else:
                        if verbose: print(f'No later file found within {diff(time).max()*3600.:.1f}s.')
                        filecheck = False        
                        if filetime:
                            return                    

            #perform initial check on variables_requested list
            if len(variables_requested)>0 and fulltime and variables_requested!='all':
                test_list = [value[0] for key, value in model_varnames.items()]
                err_list = [item for item in variables_requested if item not in test_list]
                if len(err_list)>0: print('Variable name(s) not recognized:', err_list)
                
            #collect variable list (both in attributes and in datasets)
            key_list = list(f_north[time_list[0]+'N'].attrs.keys())+list(f_north[time_list[0]+'N'].keys())
            if 'int_joule_heat' in key_list: #replace with separate names for each hemisphere
                key_list.remove('int_joule_heat')
                key_list.append('int_joule_heat_n')
                key_list.append('int_joule_heat_s')
            if len(variables_requested)>0 and variables_requested!='all':
                gvar_list = [key for key, value in model_varnames.items() \
                                 if value[0] in variables_requested and \
                                     key in key_list]  # file variable names
                    
                #check for variables requested but not available
                if len(gvar_list)!=len(variables_requested):
                    err_list = [value[0] for key, value in model_varnames.items() \
                                 if value[0] in variables_requested and \
                                     key not in gvar_list]
                    if len(err_list)>0: print('Some requested variables are not available:', err_list)
            else:  #only input variables on the avoid_list if specifically requested
                gvar_list = [key for key in key_list if key in model_varnames.keys()]
                if not fulltime and variables_requested=='all':  #returns list of variables included in data files
                    self.var_dict = {value[0]: value[1:] for key, value in model_varnames.items() \
                            if key in gvar_list}
                    return 

            #store data for each variable desired
            f_south = h5py.File(full_filenameS, 'r')
            variables = {model_varnames[var][0]: {'units': model_varnames[var][-1], 
                                                       'data': 0.} for var in gvar_list}
            attrs_var = ['imf_By','imf_Bz','solar_wind_speed'] #list of variables to be retrieved from attrs
            lon = unique(f_north['lons'])  #NEED lon here for longitude swapping logic
            lon_le180 = where(lon<=180)[0]
            lon_ge180 = where(lon>=180)[0][:-1]  #ignore 360 values, repeat 180 for -180 values 
            for var in gvar_list:  
                if var in attrs_var:  #save time series data from attributes
                    variables[model_varnames[var][0]]['data'] = \
                        array([f_north[time+'N'].attrs[var] for time in time_list], dtype=float)
                elif var=='int_joule_heat_n':  #only time series variable in north dataset
                    variables[model_varnames[var][0]]['data'] = \
                        array([array(f_north[time+'N']['int_joule_heat'])[0] for time in time_list], dtype=float)
                elif var=='int_joule_heat_s':  #only time series variable in south dataset
                    variables[model_varnames[var][0]]['data'] = \
                        array([array(f_south[time+'S']['int_joule_heat'])[0] for time in time_list], dtype=float)
                else:  #pull from datasets
                    north_data = array([array(f_north[time+'N'][var]).T for time in time_list], dtype=float)
                    south_data = array([flip(array(f_south[time+'S'][var]),axis=0).T for time in time_list], dtype=float)
                    #need to reverse order along latitude axis for southern hemisphere array
                    if '_th' in var: north_data*= -1 #change from equatorward to northward in N hemisphere data
                    total_data = concatenate((south_data,north_data),axis=2)
                    
                    #swap longitudes, repeat 180 values for -180 position
                    variables[model_varnames[var][0]]['data'] = zeros(total_data.shape, dtype=float)
                    variables[model_varnames[var][0]]['data'][:,:len(lon_ge180),:] = total_data[:,lon_ge180,:]
                    variables[model_varnames[var][0]]['data'][:,len(lon_ge180):,:] = total_data[:,lon_le180,:]
                
            #prepare and return data 
            if not fulltime:  
                f_north.close()
                f_south.close()
                variables['time'] = self.filetimes[0]
                self.short_data = variables
                return        
    
            #### Store coordinate data as class attributes   
            if filecheck:
                new_time = ts_to_hrs(short_data['time'], self.filedate)  #new time in hours since midnight
                self._time = append(time, new_time) #append new time in hours since midnight
            else: 
                self._time = time  
                
            #collect and arrange lat grid to be increasing (from neg to pos)
            lat_N = unique(f_north['lats'])  #24 values
            lat_S = flip(unique(f_south['lats']))*(-1)  #reverse order and make negative
            lat = append(lat_S,lat_N)  #-88.?? to +88.?? (south pole first, north pole last)
            lat = insert(lat, 0, -90.)  #insert south pole value
            self._lat = append(lat, 90.)  #insert north pole value
            
            #rearrange lon grid to be from -180 to 180 to match data
            self._lon = zeros(lon.shape)
            self._lon[:len(lon_ge180)] = lon[lon_ge180]-360.
            self._lon[len(lon_ge180):] = lon[lon_le180]
            
            #convert height in km to radius in R_E to conform to MAG coord sys in spacepy
            self._radius = (110.+R_earth.value/1000.)/(R_earth.value/1000.)  #110 km altitude
            f_north.close()   #close files
            f_south.close()
            
            #store a few items in object
            self.missing_value = NaN
            self._registered = 0
            if verbose: print(f'Took {perf_counter()-t0:.6f}s to read in data')
            if printfiles: print(self.filename)
            
            #need latitude wrapping (scalars and vectors): TIEGCM reader, wrap_3Dlatlon and vector_average4D functions 
            # register interpolators for each requested variable
            t_reg = perf_counter()
            varname_list, self.variables = [key for key in variables.keys()], {}  #store original list b/c gridded interpolators        
            for varname in varname_list:
                if len(variables[varname]['data'].shape)==1:
                    if filecheck:  #if neighbor found
                        #append data for last time stamp, transpose and register
                        data_shape = list(variables[varname]['data'].shape)
                        data_shape[0]+=1  #add space for time
                        new_data = zeros(data_shape)                
                        new_data[:-1] = variables[varname]['data']  #put in current data
                        new_data[-1] = short_data[varname]['data'][0]  #add in data for additional time
                        variable = new_data  #no transposing needed
                    else:
                        variable = variables[varname]['data']
                    self.variables[varname] = dict(units = variables[varname]['units'], data = variable)
                    self.register_1D_variable(self.variables[varname]['units'], 
                                          self.variables[varname]['data'], varname,
                                          gridded_int)
                elif len(variables[varname]['data'].shape)==3:
                    if filecheck:
                        #append data for last time stamp, transpose and register
                        data_shape = list(variables[varname]['data'].shape)
                        data_shape[0]+=1  #add space for time
                        new_data = zeros(data_shape)
                        new_data[:-1,:,:] = variables[varname]['data']  #put in current data
                        new_data[-1,:,:] = short_data[varname]['data'][0,:,:]   #add in data for additional time                
                        variable = new_data  #no transposing needed
                    else:
                        variable = variables[varname]['data']
                    self.variables[varname] = dict(units = variables[varname]['units'], data = variable)
                    self.register_3D_variable(self.variables[varname]['units'], 
                                          self.variables[varname]['data'], varname,
                                          gridded_int)
            if verbose: print(f'Took {perf_counter()-t_reg:.5f}s to register '+\
                              f'{len(varname_list)} variables.')
            if verbose: print(f'Took a total of {perf_counter()-t0:.5f}s to kamodofy '+\
                              f'{len(varname_list)} variables.')            
                
        def wrap3Dlat(self, varname, variable):
            '''Wrap latitude values for scalars and vectors.'''
            
            shape_list = list(variable.shape)  #time, lon, lat
            shape_list[2]+=2  #need two more places in latitude
            tmp_arr = zeros(shape_list)  #array to set-up wrapped data in
            tmp_arr[:,:,1:-1] = variable  #copy data into grid
    
            #wrapping in latitude for scalar variables
            if '_th' not in varname and '_ph' not in varname:
                # put in top values
                top = mean(tmp_arr[:,:-1,1],axis=1)  #same shape as time axis
                tmp_arr[:,:-1,0] = broadcast_to(top, (shape_list[1]-1,shape_list[0])).T
                #same for bottom, reusing variable names
                top = mean(tmp_arr[:,:-1,-2],axis=1)  #same shape as time axis
                tmp_arr[:,:-1,-1] = broadcast_to(top, (shape_list[1]-1,shape_list[0])).T
                
            #wrapping in latitude for relevant vector variables
            elif '_th' in varname or '_ph' in varname:
                #calculate net vector magnitude for top
                tmp_arr[:,:-1,0] = self.vector_average3D(tmp_arr[:,:-1,1], 
                                                      shape_list, varname, self._lat[0])
                #repeat for bottom
                tmp_arr[:,:-1,-1] = self.vector_average3D(tmp_arr[:,:-1,-2],
                                                      shape_list, varname, self._lat[-1])
            tmp_arr[:,-1,:] = tmp_arr[:,0,:]  #wrap value in longitude after to prevent double-counting           
            self.variables[varname]['data'] = tmp_arr  #store result
            return tmp_arr    
    
        def vector_average3D(self, top, shape_list, varname, latval):
            '''find vector average at pole for array with shape (time, lon)'''
        
            #find net x and y components, final array shapes are (time, lon)
            lon_arr = broadcast_to(self._lon[:-1], (shape_list[0],shape_list[1]-1))
            xval = sum(top*cos((lon_arr+180.)*nppi/180.), axis=1)  #sum over lon axis
            yval = sum(top*sin((lon_arr+180.)*nppi/180.), axis=1)  #same shape as time
            xarr = broadcast_to(xval, (shape_list[1]-1,shape_list[0])).T  #xval.shape must be last in broadcast_to call
            yarr = broadcast_to(yval, (shape_list[1]-1,shape_list[0])).T

            #convert to proper unit vector (see wiki on spherical coordinates)
            if '_ph' in varname:  #Zonal / east components -> convert to psi_hat vector (longitude)
                # -xsin(psi)+ycos(psi), psi = longitude (0 to 360)
                new_top = -xarr*sin((lon_arr+180.)*nppi/180.)+yarr*cos((lon_arr+180.)*nppi/180.)
            elif '_th' in varname:  #meridional/north -> convert to theta_hat vector (latitude)
                # xcos(psi)cos(theta)+ysin(psi)cos(theta),  sin(theta) is always zero at the poles
                # theta = latitude (0 to 180), psi = longitude (0 to 360)
                new_top = xarr*cos((lon_arr+180.)*nppi/180.)*cos((90.-latval)*nppi/180.)+\
                    yarr*sin((lon_arr+180.)*nppi/180.)*cos((90.-latval)*nppi/180.)
            
            return new_top    
    
        #define and register a 1D variable
        def register_1D_variable(self, units, variable, varname, gridded_int):
            """Registers a 1d interpolator with 1d signature"""
            
            #define and register the interpolators
            xvec_dependencies = {'time':'hr'}
            self = regdef_1D_interpolators(self, units, variable, self._time, 
                                           varname, xvec_dependencies, gridded_int)
            return 
        
        #define and register a 3D variable
        def register_3D_variable(self, units, variable, varname, gridded_int):
            """Registers a 3d interpolator with 3d signature"""
            
            #define and register the fast interpolator
            xvec_dependencies = {'time':'hr','lon':'deg','lat':'deg'}
            variable = self.wrap3Dlat(varname, variable)
            self = regdef_3D_interpolators(self, units, variable, self._time,
                                              self._lon, self._lat,
                                              varname, xvec_dependencies, gridded_int)
            return
    return MODEL
