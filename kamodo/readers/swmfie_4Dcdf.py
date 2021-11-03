# -*- coding: utf-8 -*-
"""
Created on Mon May 17 18:53:35 2021

@author: rringuet
"""
from datetime import datetime, timezone
from numpy import vectorize
from os.path import isfile, basename
#read 1 day of data from cdf instead of from multiple .tec files


model_varnames={"Sigma_H":['Sigma_H','variable description',0,'SM','sph',['time','lon','lat'],"S"],
                "Sigma_P":['Sigma_P','variable description',1,'SM','sph',['time','lon','lat'],"S"],
                 "Phi_E":['Phi_E','variable description',2,'SM','sph',['time','lon','lat'],"W/m**2"], 
                 "AveE_avgE":['E_avg','variable description',3,'SM','sph',['time','lon','lat'],'eV'],
                 "j_R":["j_R",'variable description',4,'SM','sph',['time','lon','lat'],"muA/m**2"],
                 "Phi":["Phi",'variable description',5,'SM','sph',['time','lon','lat'],"kV"],
                 "E_x":["E_x",'variable description',6,'SM','sph',['time','lon','lat'],"mV/m"],
                 "E_y":["E_y",'variable description',7,'SM','sph',['time','lon','lat'],"mV/m"],
                 "E_z":["E_z",'variable description',8,'SM','sph',['time','lon','lat'],"mV/m"],
                 "j_x":["j_x",'variable description',9,'SM','sph',['time','lon','lat'],"muA/m**2"],
                 "j_y":["j_y",'variable description',10,'SM','sph',['time','lon','lat'],"muA/m**2"],
                 "j_z":["j_z",'variable description',11,'SM','sph',['time','lon','lat'],"muA/m**2"],
                 "v_x":['v_x','variable description',12,'SM','sph',['time','lon','lat'],"km/s"],
                 "v_y":['v_y','variable description',13,'SM','sph',['time','lon','lat'],"km/s"],
                 "v_z":['v_z','variable description',14,'SM','sph',['time','lon','lat'],"km/s"],
                 "Q_Joule":['Q_Joule','variable description',15,'SM','sph',['time','lon','lat'],"mW/m**2"], 
                 "Phi_nion":['Phi_nion','variable description',16,'SM','sph',['time','lon','lat'],"1/cm**2/s"],
                 "Binv_RT":['Binv_RT','variable description',17,'SM','sph',['time','lon','lat'],"1/T"],
                 "rho_RT":['rho_RT','variable description',18,'SM','sph',['time','lon','lat'],"amu/cm**3"],
                 "P_RT":['P_RT','variable description',19,'SM','sph',['time','lon','lat'],"Pa"],
                 "dLat_star":['dLat_star','variable description',20,'SM','sph',['time','lon','lat'],"deg"],
                 "dlon_star":['dlon_star','variable description',21,'SM','sph',['time','lon','lat'],"deg"]}
                 
 
'''                
#documentation variables:
"X":['x','km'],"Y":['y','km'],"Z":['z','km'],    (ignored, given in R_E on a unit sphere)
"Theta":['theta',"deg"],"Psi":['psi',"deg"],     (used as coordinates)
"Btilt_theta":['theta_Btilt',"deg"], "Btilt_psi":['psi_Btilt',"deg"] 
(added directly to object for documentation purposes)
'''
@vectorize
def dts_to_hrs(datetime_string, filedate):
    '''Get hours since midnight from datetime string'''
    
    return (datetime.strptime(datetime_string, '%Y-%m-%d %H:%M:%S').replace(tzinfo=timezone.utc)\
            -filedate).total_seconds()/3600.

@vectorize
def filename_to_dts(filename, string_date):
    '''Get datetime string in format "YYYY-MM-SS HH:mm:SS" from filename'''
    
    mmhhss = basename(filename)[12:18]
    return string_date+' '+mmhhss[:2]+':'+mmhhss[2:4]+':'+mmhhss[4:] 

def dts_to_ts(file_dts):
    '''Get datetime timestamp in UTC from datetime string'''    
    
    return datetime.timestamp(datetime.strptime(file_dts, '%Y-%m-%d %H:%M:%S'
                                                ).replace(tzinfo=timezone.utc))
 
def ts_to_hrs(time_val, filedate):
    '''Convert utc timestamp to hours since midnight on filedate.'''
    
    return (datetime.utcfromtimestamp(time_val).replace(tzinfo=timezone.utc)-filedate).total_seconds()/3600.
   
#base for main class object----------------------------------------
#main class object
#file_dir = 'C:/Users/rringuet/Kamodo_WinDev1/SWMF_IE/'
#file_prefix = 'C:/Users/rringuet/Kamodo_WinDev1/SWMF_IE/i_e20180801'  #example
#files = glob.glob(file_dir+'i_e*.tec')  #for wrapper, this and the next line
#file_patterns = unique([file_dir+f.split('/')[-1].split('\\')[-1][:11] for f in files])
def MODEL():
    from numpy import array, NaN, abs, unique, append, zeros, diff, where, insert, flip
    from time import perf_counter
    from netCDF4 import Dataset
    from kamodo import Kamodo
    #print('KAMODO IMPORTED!')
    from kamodo.readers.reader_utilities import regdef_3D_interpolators    
       
    class MODEL(Kamodo): 
        def __init__(self, full_file_prefix, variables_requested=[], runname="noname",
                     filetime=False, verbose=False, gridded_int=True, printfiles=False,
                     fulltime=True, **kwargs): 
            '''file_prefix must be of form "3D***_tYYMMDD" to load all files for one day
             and include a complete path to the files'''
            super(MODEL, self).__init__()
            
            #check if given .nc file exists. If not, convert files with same prefix to netCDF
            file_prefix = basename(full_file_prefix)
            file_dir = full_file_prefix.split(file_prefix)[0]   
            if not isfile(full_file_prefix+'.nc'):
                from kamodo.readers.swmfie_tocdf import convertSWMFIE_toCDF
                test = convertSWMFIE_toCDF(full_file_prefix)  
                if not test: 
                    self.conversion_test = test  #only needed for 1 file/time cases
                    return    #if file conversion fails, return 
                else: self.conversion_test = True
            else: self.conversion_test = True
            t0 = perf_counter()
            
            #determine type of prefix: for a day or for a hour
            if '-' in file_prefix: day_flag=False
            else: day_flag=True
            
            #establish time attributes first for file searching
            file_datestr = file_prefix[3:11]
            string_date = file_datestr[:4]+'-'+file_datestr[4:6]+'-'+file_datestr[6:8]  #'YYYY-MM-DD'
            self.filedate = datetime.strptime(string_date+' 00:00:00', \
                                              '%Y-%m-%d %H:%M:%S').replace(tzinfo=timezone.utc) #dt object
            
            #establish beginning and end time of file list
            cdf_data = Dataset(full_file_prefix+'.nc', 'r')
            files = cdf_data.file.split(',')
            self.datetimes = list(filename_to_dts([files[0], files[-1]], string_date))  #strings in format = YYYY-MM-DD HH:MM:SS 
            self.filetimes=[dts_to_ts(file_dts) for file_dts in self.datetimes]   #timestamps in UTC  
            t = array(cdf_data.variables['time'])  #hours since midnight
            if len(t)>1: self.dt = diff(t).max()*3600.  #t is in hours since midnight
            else: self.dt = 0
            
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
                
                files = sorted(glob(file_dir+'i_e*'))
                if day_flag: 
                    file_prefixes = unique([basename(f)[:11] for f in files\
                                            if '.nc' not in basename(f)])
                else:  #give prefix for hourly files
                    file_prefixes = unique([basename(f)[:14] for f in files\
                                            if '.nc' not in basename(f)])
                
                #find closest file by utc timestamp
                #swmf_ie has an open time at the end, so need a beginning time from the next file
                #files are automatically sorted by YYMMDD, so next file is next in the list
                current_idx = where(file_prefixes==file_prefix)[0]
                if current_idx+1==len(file_prefixes):
                    print('No later file available.')
                    filecheck = False  
                    if filetime:
                        return   
                else:
                    min_file_prefix = file_dir+file_prefixes[current_idx+1][0]  #+1 for adding an end time
                    kamodo_test = MODEL(min_file_prefix, filetime=True, fulltime=False)
                    if not kamodo_test.conversion_test: 
                        print('No later file available.')
                        filecheck = False  
                        if filetime:
                            return        
                    else:
                        time_test = abs(kamodo_test.filetimes[0]-self.filetimes[1])  
                        if time_test<=self.dt:  #if nearest file time at least within one timestep (hrs)
                            filecheck = True
                        
                            #time only version if returning time for searching
                            if filetime:
                                kamodo_neighbor = MODEL(min_file_prefix, fulltime=False, filetime=True)
                                self.datetimes[1] = kamodo_neighbor.datetimes[0]
                                self.filetimes[1] = kamodo_neighbor.filetimes[0]
                                return  #return object with additional time (for SF code) 
                            
                            #get kamodo object with same requested variables to add to each array below
                            if verbose: print(f'Took {perf_counter()-t0:.3f}s to find closest file.')
                            kamodo_neighbor = MODEL(min_file_prefix, 
                                                    variables_requested=variables_requested, 
                                                    fulltime=False)
                            self.datetimes[1] = kamodo_neighbor.datetimes[0]
                            self.filetimes[1] = kamodo_neighbor.filetimes[0]
                            short_data = kamodo_neighbor.short_data                                
                            if verbose: print(f'Took {perf_counter()-t0:.3f}s to get data from closest file.')
                        else:
                            print(f'No later file found within {self.dt:.1f}s.')
                            filecheck = False 
                            if filetime:
                                return                    

            #perform initial check on variables_requested list
            if len(variables_requested)>0 and fulltime:
                test_list = [value[0] for key, value in model_varnames.items()]
                err_list = [item for item in variables_requested if item not in test_list]
                if len(err_list)>0: print('Variable name(s) not recognized:', err_list)
                
            #get list of variables possible in these files using first file
            if len(variables_requested)>0:
                gvar_list = [key for key, value in model_varnames.items() \
                                 if value[0] in variables_requested and \
                                     key in cdf_data.variables.keys()]  # file variable names
                #check for variables requested but not available
                if len(gvar_list)!=len(variables_requested):
                    err_list = [value[0] for key, value in model_varnames.items() \
                                 if value[0] in variables_requested and \
                                     key not in cdf_data.variables.keys()]
                    if len(err_list)>0: print('Some requested variables are not available:', err_list)
            else:
                avoid_list = ['theta_Btilt', 'psi_Btilt']
                gvar_list = [key for key in cdf_data.variables.keys() \
                             if key in model_varnames.keys() and \
                                 key not in avoid_list]
                                   
            # Store variable's data and units, transposing the 2D+time array.
            variables = {model_varnames[key][0]:{'units':model_varnames[key][-1],
                                   'data':array(cdf_data.variables[key])}\
                              for key in gvar_list} 
            self.theta_Btilt = array(cdf_data.variables['theta_Btilt'])
            self.psi_Btilt = array(cdf_data.variables['psi_Btilt'])
                
            #prepare and return data only for first timestamp
            if not fulltime:  
                cdf_data.close()
                variables['time'] = self.filetimes[0]
                variables['theta_Btilt'] = self.theta_Btilt[0]
                variables['psi_Btilt'] = self.psi_Btilt[0]
                self.short_data = variables
                return            
    
            #return if only one file found because interpolator code will break
            if len(files)<2:
                print('Not enough files found with given file prefix.')
                return 
    
            #store variables
            self.filename = files
            self.runname = runname
            self.missing_value = NaN
            self.modelname = 'SWMF_IE'
            self._registered = 0
            if printfiles: 
                print('Files:')
                for file in self.filename: print(file)
    
            #### Store coordinate data as class attributes   
            if filecheck:
                new_time = ts_to_hrs(short_data['time'], self.filedate)  #new time in hours since midnight
                self._time = append(t, new_time) 
            else: 
                self._time = t
                
            #store coordinate data
            #self._radius = array(cdf_data.variables['radius'])
            self._lat = array(cdf_data.variables['lat'])  
            self._lon = array(cdf_data.variables['lon'])
            cdf_data.close()
            if verbose: print(f'Took {perf_counter()-t0:.6f}s to read in data')
    
            #register interpolators for each variable
            varname_list, self.variables = [key for key in variables.keys()], {}  #store original list b/c gridded interpolators
            t_reg = perf_counter()
            for varname in varname_list:  #all are 3D variables
                if filecheck:  #if neighbor found
                    #append data for first time stamp, transpose and register
                    data_shape = list(variables[varname]['data'].shape)
                    data_shape[0]+=1  #add space for time
                    new_data = zeros(data_shape)                
                    new_data[:-1,:,:] = variables[varname]['data']  #put in current data
                    new_data[-1,:,:] = short_data[varname]['data'][0,:,:]  #add in data for additional time
                else:
                    new_data = variables[varname]['data']
                self.variables[varname] = dict(units = variables[varname]['units'], data = new_data)           
                self.register_3D_variable(self.variables[varname]['units'], 
                                          self.variables[varname]['data'], varname,
                                          gridded_int)
            if verbose: print(f'Took {perf_counter()-t_reg:.5f}s to register '+\
                              f'{len(varname_list)} variables.')
            if verbose: print(f'Took a total of {perf_counter()-t0:.5f}s to kamodofy '+\
                              f'{len(gvar_list)} variables.')
        
        #define and register a 3D variable-----------------------------------------
        def register_3D_variable(self, units, variable, varname, gridded_int):
            """Registers a 3d interpolator with 3d signature"""
            
            #define and register the interpolators
            xvec_dependencies = {'time':'hr','lon':'deg','lat':'deg'}
            self = regdef_3D_interpolators(self, units, variable, self._time, 
                                           self._lon, self._lat, varname, 
                                           xvec_dependencies, gridded_int)       
            return 
    return MODEL