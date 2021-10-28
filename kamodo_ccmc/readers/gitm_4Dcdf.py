from datetime import datetime, timezone, timedelta


#file_prefix = 'C:/Users/rringuet/Kamodo_WinDev1/GITM/3DLST_t150317'
#keep for easy access of variable names from satellite flythrough software
#line for MODEL wrapper if data scale isn't important:
# return {key:[value[0],value[2]] for key, value in model_varnames.items()}
#varnames in cdf files are standardized (value[0])
model_varnames = {'r_Ar': ['r_Ar','variable description',0,'SPH','sph',['time','lon','lat','radius'], ''],
                 'rho_Ar': ['rho_Ar','variable description',1,'SPH','sph',['time','lon','lat','radius'], 'kg/m**3'],
                 'r_CH4': ['r_CH4','variable description',2,'SPH','sph',['time','lon','lat','radius'], ''],
                 'k': ['k','variable description',3,'SPH','sph',['time','lon','lat','radius'], 'W/m/K'],
                 'Q_EUV': ['Q_EUV','variable description',4,'SPH','sph',['time','lon','lat','radius'], 'K per timestep'],
                 'rho_H': ['rho_H','variable description',5,'SPH','sph',['time','lon','lat','radius'], 'kg/m**3'],
                 'rho_Hplus': ['rho_Hplus','variable description',6,'SPH','sph',['time','lon','lat','radius'], 'kg/m**3'],
                 'r_H2': ['r_H2','variable description',7,'SPH','sph',['time','lon','lat','radius'], ''],
                 'r_HCN': ['r_HCN','variable description',8,'SPH','sph',['time','lon','lat','radius'], ''],
                 'rho_He': ['rho_He','variable description',9,'SPH','sph',['time','lon','lat','radius'], 'kg/m**3'],
                 'rho_Heplus': ['rho_Heplus','variable description',10,'SPH','sph',['time','lon','lat','radius'], 'kg/m**3'],
                 'HeatingEfficiency': ['HeatingEfficiency','variable description',11,'SPH','sph',['time','lon','lat','radius'], ''],
                 'HeatBalanceTotall': ['HeatBalanceTotal','variable description',12,'SPH','sph',['time','lon','lat','radius'], ''],
                 'rho_N2': ['rho_N2','variable description',13,'SPH','sph',['time','lon','lat','radius'], 'kg/m**3'],
                 'rho_N2plus': ['rho_N2plus','variable description',14,'SPH','sph',['time','lon','lat','radius'], 'kg/m**3'],
                 'rho_Nplus': ['rho_Nplus','variable description',15,'SPH','sph',['time','lon','lat','radius'], 'kg/m**3'],
                 'rho_N2D': ['rho_N2D','variable description',16,'SPH','sph',['time','lon','lat','radius'], 'kg/m**3'],
                 'rho_N2P': ['rho_N2P','variable description',17,'SPH','sph',['time','lon','lat','radius'], 'kg/m**3'],
                 'rho_N4S': ['rho_N4S','variable description',18,'SPH','sph',['time','lon','lat','radius'], 'kg/m**3'],
                 'r_N2': ['r_N2','variable description',19,'SPH','sph',['time','lon','lat','radius'], ''],
                 'rho_NO': ['rho_NO','variable description',20,'SPH','sph',['time','lon','lat','radius'], 'kg/m**3'],
                 'rho_NOplus': ['rho_NOplus','variable description',21,'SPH','sph',['time','lon','lat','radius'], 'kg/m**3'],
                 'rho_O2': ['rho_O2','variable description',22,'SPH','sph',['time','lon','lat','radius'], 'kg/m**3'],
                 'rho_O1D': ['rho_O1D', 'variable description',23,'SPH','sph',['time','lon','lat','radius'],'kg/m**3'],
                 'rho_O2plus': ['rho_O2plus','variable description',24,'SPH','sph',['time','lon','lat','radius'], 'kg/m**3'],
                 'rho_O2D': ['rho_O2D','variable description',25,'SPH','sph',['time','lon','lat','radius'], '1/m**3'],
                 'rho_Oplus2P': ['rho_Oplus2P','variable description',26,'SPH','sph',['time','lon','lat','radius'], 'kg/m**3'],
                 'rho_O3P': ['rho_O3P','variable description',27,'SPH','sph',['time','lon','lat','radius'], 'kg/m**3'],
                 'rho_Oplus4SP': ['rho_Oplus4SP','variable description',28,'SPH','sph',['time','lon','lat','radius'], 'kg/m**3'],
                 'L_Rad': ['L_Rad', 'variable description',29,'SPH','sph',['time','lon','lat','radius'],''],
                 'rho': ['rho_n','variable description',30,'SPH','sph',['time','lon','lat','radius'], 'kg/m**3'],
                 'T_n': ['T_n','variable description',31,'SPH','sph',['time','lon','lat','radius'], 'K'],
                 'vi_east': ['vi_east','variable description',32,'SPH','sph',['time','lon','lat','radius'], 'm/s'],
                 'vi_north': ['vi_north','variable description',33,'SPH','sph',['time','lon','lat','radius'], 'm/s'],
                 'vi_up': ['vi_up','variable description',34,'SPH','sph',['time','lon','lat','radius'], 'm/s'],
                 'vn_east': ['vn_east','variable description',35,'SPH','sph',['time','lon','lat','radius'], 'm/s'],
                 'vn_north': ['vn_north','variable description',36,'SPH','sph',['time','lon','lat','radius'], 'm/s'],
                 'vn_up': ['vn_up','variable description',37,'SPH','sph',['time','lon','lat','radius'], 'm/s'],
                 'v_N2_up': ['v_N2_up','variable description',38,'SPH','sph',['time','lon','lat','radius'], 'm/s'],
                 'v_N4S_up': ['v_N4S_up','variable description',39,'SPH','sph',['time','lon','lat','radius'], 'm/s'],
                 'v_N_up': ['v_N_up','variable description',40,'SPH','sph',['time','lon','lat','radius'], 'm/s'],
                 'v_O2_up': ['v_O2_up','variable description',41,'SPH','sph',['time','lon','lat','radius'], 'm/s'],
                 'v_O3P_up': ['v_O3P_up','variable description',42,'SPH','sph',['time','lon','lat','radius'], 'm/s'],
                 'v_He_up': ['v_He_up','variable description',43,'SPH','sph',['time','lon','lat','radius'], 'm/s'],
                 'N_e': ['N_e','variable description',44,'SPH','sph',['time','lon','lat','radius'], '1/m**3'],
                 'ElectronAverageEnergy': ['ElectronAverageEnergy','variable description',45,'SPH','sph',['time','lon','lat','radius'], 'J'],
                 'T_e': ['T_e','variable description',46,'SPH','sph',['time','lon','lat','radius'], 'K'],
                 'T_i': ['T_i','variable description',47,'SPH','sph',['time','lon','lat','radius'], 'K'],
                 'SolarZenithAngle': ['SolarZenithAngle','variable description',48,'SPH','sph',['time','lon','lat'], 'radians'],
                 'rho_CO2': ['rho_CO2','variable description',49,'SPH','sph',['time','lon','lat','radius'], 'kg/m**3'],
                 'DivJu FL': ['DivJuFL','variable description',50,'SPH','sph',['time','lon','lat','radius'], ''],
                 'DivJuAlt': ['DivJuAlt','variable description',51,'SPH','sph',['time','lon','lat','radius'], ''],
                 'ElectronEnergyFlux': ['ElectronEnergyFlux','variable description',52,'SPH','sph',['time','lon','lat','radius'], 'J/m**2'],
                 'Field Line Length': ['Field Line Length','variable description',53,'SPH','sph',['time','lon','lat','radius'], 'm'],
                 'sigma_P': ['sigma_P','variable description',54,'SPH','sph',['time','lon','lat','radius'], 'S/m'],
                 'V': ['V','variable description',57,'SPH','sph',['time','lon','lat','radius'], 'V'],
                 'sigma_H': ['sigma_H','variable description',58,'SPH','sph',['time','lon','lat','radius'], 'S/m'],
                 'I_R2': ['I_R2','variable description',59,'SPH','sph',['time','lon','lat','radius'], 'A/m**2'],
                 'I_R1': ['I_R1','variable description',60,'SPH','sph',['time','lon','lat','radius'], 'A/m**2'],
                 'Ed1': ['Ed1','variable description',61,'SPH','sph',['time','lon','lat','radius'], ''],
                 'Ed2': ['Ed2','variable description',62,'SPH','sph',['time','lon','lat','radius'], ''],
                 'SolarLocalTime': ['SolarLocalTime','variable description',63,'SPH','sph',['time','lon','lat'], 'hr'],
                 'E_up': ['E_up','variable description',64,'SPH','sph',['time','lon','lat','radius'], 'V/m'],
                 'E_east': ['E_east','variable description',65,'SPH','sph',['time','lon','lat','radius'], 'V/m'],
                 'E_north': ['E_north','variable description',66,'SPH','sph',['time','lon','lat','radius'], 'V/m'],
                 'E_mag': ['E_mag','variable description',67,'SPH','sph',['time','lon','lat','radius'], 'V/m'],
                 'B_up': ['B_up','variable description',68,'SPH','sph',['time','lon','lat','radius'], 'nT'],
                 'B_east': ['B_east','variable description',69,'SPH','sph',['time','lon','lat','radius'], 'nT'],
                 'B_north': ['B_north','variable description',70,'SPH','sph',['time','lon','lat','radius'], 'nT'],
                 'B_mag': ['B_mag','variable description',71,'SPH','sph',['time','lon','lat','radius'], 'nT'],
                 'MagLat': ['MagLat','variable description',72,'SPH','sph',['time','lon','lat','radius'], 'deg'],
                 'MagLon': ['MagLon','variable description',73,'SPH','sph',['time','lon','lat','radius'], 'deg'],
                 'g': ['g','variable description',74,'SPH','sph',['time','lon','lat','radius'], 'm/s**2'],
                 'GradP_east': ['GradP_east','variable description',75,'SPH','sph',['time','lon','lat','radius'], 'Pa/m'],
                 'GradP_north': ['GradP_north','variable description',76,'SPH','sph',['time','lon','lat','radius'], 'Pa/m'],
                 'GradP_up': ['GradP_up','variable description',77,'SPH','sph',['time','lon','lat','radius'], 'Pa/m'],
                 'nu_in': ['nu_in','variable description',78,'SPH','sph',['time','lon','lat','radius'], '1/s'],
                 'ChemicalHeatingRate': ['ChemicalHeatingRate','variable description',79,'SPH','sph',['time','lon','lat','radius'], ''],
                 'TotalAbsoluteEUV': ['TotalAbsoluteEUV','variable description',80,'SPH','sph',['time','lon','lat','radius'], 'K per timestep'],
                 'Q_Ocool': ['Q_Ocool','variable description',81,'SPH','sph',['time','lon','lat','radius'], 'K per timestep'],
                 'Q_Joule': ['Q_Joule','variable description',82,'SPH','sph',['time','lon','lat','radius'], 'K per timestep'],
                 'Q_Auroral': ['Q_Auroral','variable description',83,'SPH','sph',['time','lon','lat','radius'], 'K per timestep'],
                 'Q_PhotoE': ['Q_PhotoE','variable description',84,'SPH','sph',['time','lon','lat','radius'], 'K per timestep'],
                 'k_eddy': ['k_eddy','variable description',85,'SPH','sph',['time','lon','lat','radius'], ''],
                 'k_adiabaticeddy': ['k_adiabaticeddy','variable description',86,'SPH','sph',['time','lon','lat','radius'], ''],
                 'Q_NOcool': ['Q_NOcool','variable description',87,'SPH','sph',['time','lon','lat','radius'], 'K per timestep'],
                 'k_molecular': ['k_molecular','variable description',88,'SPH','sph',['time','lon','lat','radius'], ''],
                 'NmF2': ['NmF2','variable description',89,'SPH','sph',['time','lon','lat'], ''],
                 'hmF2': ['hmF2','variable description',90,'SPH','sph',['time','lon','lat'], 'km'],
                 'TEC': ['TEC','variable description',91,'SPH','sph',['time','lon','lat'], '10**16/m**2'],
                 'Phi_Joule':['Phi_Joule','variable description',92,'SPH','sph',['time','lon','lat'],'W/m**2'], 
                 'Phi_Q':['Phi_Q','variable description',93,'SPH','sph',['time','lon','lat'],'W/m**2'], 
                 'Phi_EUV':['Phi_EUV','variable description',94,'SPH','sph',['time','lon','lat'],'W/m**2'], 
                 'Phi_NOCooling':['Phi_NOCooling','variable description',95,'SPH','sph',['time','lon','lat'],'W/m**2']}

#times from file converted to seconds since midnight of filedate
#plotting input times will be datetime strings of format 'YYYY-MM-DD HH:mm:ss'
#filedate is self.filedate from MODEL object
#converts to hours since midnight of filedate for plotting
#def dtsns_to_hrs(datetime_string, filedate):
#    '''Get hours since midnight from datetime string including nanoseconds'''
#    
#    return (datetime.strptime(datetime_string, '%Y-%m-%d %H:%M:%S.%f').replace(tzinfo=timezone.utc)\
#            -filedate).total_seconds()/3600.

def dts_to_hrs(datetime_string, filedate):
    '''Get hours since midnight from datetime string'''
    
    return (datetime.strptime(datetime_string, '%Y-%m-%d %H:%M:%S').replace(tzinfo=timezone.utc)\
            -filedate).total_seconds()/3600.

def filename_to_dts(filename, string_date):
    '''Get datetime string in format "YYYY-MM-SS HH:mm:SS" from filename'''
    
    mmhhss = filename.split('_t')[-1].split('_')[1].split('.bin')[0]
    return string_date+' '+mmhhss[:2]+':'+mmhhss[2:4]+':'+mmhhss[4:]  

def dts_to_ts(file_dts):
    '''Get datetime timestamp in UTC from datetime string'''    
    
    return datetime.timestamp(datetime.strptime(file_dts, '%Y-%m-%d %H:%M:%S'
                                                ).replace(tzinfo=timezone.utc))

def hrs_to_ts(hrs, filedate):
    '''Add hours to filedate and return utc timestamp.'''
    
    return datetime.timestamp(filedate+timedelta(hours=hrs))

def ts_to_hrs(time_val, filedate):
    '''Convert utc timestamp to hours since midnight on filedate.'''
    
    return (datetime.utcfromtimestamp(time_val).replace(tzinfo=timezone.utc)-filedate).total_seconds()/3600.

def MODEL():
    from time import perf_counter
    from glob import glob
    from os.path import basename
    from numpy import array, unique, NaN, diff, abs, append, zeros, where
    from netCDF4 import Dataset
    from kamodo import Kamodo
    #print('KAMODO IMPORTED!')
    from kamodo_ccmc.readers.reader_utilities import regdef_4D_interpolators, regdef_3D_interpolators    

    #main class object
    class MODEL(Kamodo): 
        '''GITM model data reader.'''
        def __init__(self, full_file_prefix, variables_requested=[], runname="noname",
                     filetime=False, verbose=False, gridded_int=True, printfiles=False,
                     fulltime=True,**kwargs): 
            '''filename must be of form "***_tYYMMDD" to load all files for one day
            or of form "***_tYYMMDD_HH" to load only files for one hour
             and must include a complete path to the files with the files stored in the 
             "Data" directory.'''
            super(MODEL, self).__init__()
            
            #check for prepared .nc files
            file_prefix = basename(full_file_prefix)
            file_dir = full_file_prefix.split(file_prefix)[0]
            total_files = sorted(glob(file_dir+file_prefix+'*'))
            if len(file_prefix.split('_')[-1])==2:  #asking for hourly files (looping)
                self.patterns = unique([basename(f)[:16] for f in total_files \
                                        if basename(f)[13:16]!='.nc' and \
                                            'GITM' not in basename(f)])
                day_tag=False
            else:  #asking for daily files (satellite flythrough)
                self.patterns = unique([basename(f)[:13] for f in total_files\
                                        if 'GITM' not in basename(f)])
                day_tag=True
            n_ncfiles = len(glob(file_dir+file_prefix+'.nc'))
            
            #if .nc files not there, create then
            self.conversion_test=True  #default value
            if len(self.patterns)!=n_ncfiles:  #then need to create .nc files for each pattern   
                from kamodo_ccmc.readers.gitm_tocdf import GITMbin_toCDF as toCDF
                test = [toCDF(file_dir+p) for p in self.patterns]  #get/convert files with given prefix
                if sum(test)!=len(self.patterns): 
                    self.conversion_test = False
                    return   #if file conversion fails, return
            
            t0 = perf_counter()
            #establish time attributes first for file searching, preferring either 3DLST or 3DALL
            for name in self.patterns:
                if '3DLST' in name or '3DALL' in name: 
                    cdf_file = name  #need a 3D file to get reference radius correct
            cdf_data = Dataset(file_dir+cdf_file+'.nc', 'r')
            string_date = cdf_data.filedate       
            self.filedate = datetime.strptime(string_date+' 00:00:00', \
                                              '%Y-%m-%d %H:%M:%S').replace(tzinfo=timezone.utc) #dt object
    
            #establish beginning and end times of files found
            files = cdf_data.file.split(',')
            self.datetimes = [filename_to_dts(filename, string_date) for filename \
                              in [files[0], files[-1]]]  #strings in format = YYYY-MM-DD HH:MM:SS 
            self.filetimes=[dts_to_ts(file_dts) for file_dts in self.datetimes]   #timestamps in UTC   
            t = array(cdf_data.variables['time'])  #accurate to the sec
            if len(t)>1: self.dt = diff(t).max()*3600.  #t is in hours since midnight
            else: self.dt = 0
            
            if filetime and not fulltime: #(used when searching for neighboring files below)
                return  #return times as is to prevent infinite recursion
            
            #if variables are given as integers, convert to standard names
            if len(variables_requested)>0:
                if isinstance(variables_requested[0], int):
                    tmp_var = [value[0] for key, value in model_varnames.items()\
                                           if value[2] in variables_requested]
                    variables_requested = tmp_var        
            
            if fulltime:  #add boundary time (default value)
                t_search = perf_counter()
                #find other files with same pattern
                files = sorted(glob(file_dir+'*'))
                if day_tag:  #perform task with day files
                    file_prefixes = unique(['*'+basename(f)[7:13] for f in files\
                                            if 'GITM' not in basename(f) and \
                                                '.nc' not in basename(f)])
                else:   #perform task with hourly files
                    file_prefixes = unique(['*'+basename(f)[7:16] for f in files \
                                            if '.nc' not in basename(f) and\
                                                'GITM' not in basename(f)])
    
                #find closest file by utc timestamp
                #gitm has an open time at the end, so need a beginning time from the closest file
                #files are automatically sorted by YYMMDD, so next file is next in the list
                current_idx = where(file_prefixes==file_prefix)[0]
                if current_idx+1==len(file_prefixes):
                    print('No later file available.')
                    filecheck = False  
                    if filetime:
                        return   
                else:
                    min_file_prefix = file_prefixes[current_idx+1][0]  #+1 for adding an end time
                    kamodo_test = MODEL(file_dir+min_file_prefix, filetime=True, fulltime=False)
                    if not kamodo_test.conversion_test: 
                        print('No later file available.')
                        filecheck = False  
                        if filetime:
                            return 
                    else:
                        time_test = abs(kamodo_test.filetimes[0]-self.filetimes[1])  #never empty b/c includes itself 
                        if verbose: print(f'Took {perf_counter()-t_search:.5f}s to find the best file.')
                        if time_test<=self.dt:  #if nearest file time at least within one timestep (hrs)
                            filecheck = True
                        
                            #time only version if returning time for searching
                            if filetime:
                                kamodo_neighbor = MODEL(file_dir+min_file_prefix, fulltime=False, filetime=True)
                                self.datetimes[1] = kamodo_neighbor.datetimes[0]
                                self.filetimes[1] = kamodo_neighbor.filetimes[0]
                                return  #return object with additional time (for SF code) 
                            
                            #get kamodo object with same requested variables to add to each array below
                            if verbose: print(f'Took {perf_counter()-t0:.3f}s to find closest file.')
                            kamodo_neighbor = MODEL(file_dir+min_file_prefix, 
                                                    variables_requested=variables_requested, 
                                                    fulltime=False)
                            self.datetimes[1] = kamodo_neighbor.datetimes[0]
                            self.filetimes[1] = kamodo_neighbor.filetimes[0]
                            short_data = kamodo_neighbor.short_data                            
                            if verbose: print(f'Took {perf_counter()-t0:.3f}s to get data from closest file.')
                        else:
                            if verbose: print(f'No later file found within {self.dt:.1f}s.')
                            filecheck = False  
                            if filetime:
                                return                    
            
            #store coordinates for reference, time is NOT the same for all
            self._lat = array(cdf_data.variables['lat'])  #-90 to 90
            self._lon = array(cdf_data.variables['lon'])  #0 to 360
            self._radius = array(cdf_data.variables['radius'])  #need a 3D file to get this right
            cdf_data.close()
            
            #store variables
            self.filename = []
            self.runname = runname
            self.missing_value = NaN
            self._registered = 0
            self.varfiles = {}  #store which variable came from which file for easier association with coords
            self.err_list = []
            self.modelname = 'GITM'
            
            #perform initial check on variables_requested list
            if len(variables_requested)>0 and fulltime:
                test_list = [value[0] for key, value in model_varnames.items()]
                err_list = [item for item in variables_requested if item not in test_list]
                if len(err_list)>0: print('Variable name(s) not recognized:', err_list)
            
            #loop through files
            for i in range(len(self.patterns)):
                cdf_data = Dataset(file_dir+self.patterns[i]+'.nc', 'r')
                files = cdf_data.file.split(',')
                self.filename.extend(files)
                
                #store coordinates, time is always the same between files
                setattr(self, '_time'+str(i), array(cdf_data.variables['time'])) #in deg
                setattr(self, '_lat'+str(i), array(cdf_data.variables['lat'])) #in deg
                setattr(self, '_lon'+str(i), array(cdf_data.variables['lon'])) #in deg
                setattr(self, '_radius'+str(i), array(cdf_data.variables['radius'])) #in km
            
                #check var_list for variables not possible in this file set
                if len(variables_requested)>0:
                    gvar_list = [key for key in model_varnames.keys() \
                                 if key in cdf_data.variables.keys() and \
                                     model_varnames[key][0] in variables_requested]
                    if len(gvar_list)!=len(variables_requested):
                        err_list = [value[0] for key, value in model_varnames.items() \
                                    if key not in cdf_data.variables.keys() \
                                        and value[0] in variables_requested] 
                        self.err_list.extend(err_list)  #add to master list
                else:
                    gvar_list = [key for key in model_varnames.keys() \
                                 if key in cdf_data.variables.keys()]
                self.varfiles[str(i)] = [model_varnames[key][0] for key in gvar_list]  #store which file these variables came from
            
                # Store variable data, units from netCDF file.
                variables = {model_varnames[key][0]:{'units':cdf_data.variables[key].units,
                                            'data':array(cdf_data.variables[key])} \
                                  for key in gvar_list}
                cdf_data.close()
                #NOTE: this overwrites the TEC data from the 2D files with the calculated TEC data from the 3D files
    
            #prepare and return data only for last timestamp
            if not fulltime:  
                variables['time'] = self.filetimes[0]
                for i in range(len(self.patterns)):
                    variables['time'+str(i)] = hrs_to_ts(getattr(self, '_time'+str(i))[0], self.filedate)
                self.short_data = variables
                return   
    
            #remove successful variables from err_list
            self.err_list = list(unique(self.err_list))
            for item in self.err_list:
                if item in variables.keys(): self.err_list.remove(item)
            if len(self.err_list)>0:
                print('Some requested variables are not available in the files found:\n', \
                      self.patterns, self.err_list)
    
            if printfiles: 
                print(f'{len(self.filename)} Files:')
                for file in self.filename: print(file)
    
            #### Add new time and store time coordinate data as class attribute
            if filecheck:
                new_time = ts_to_hrs(short_data['time'], self.filedate)  #new time in hours since midnight
                self._time = append(t, new_time) 
                for i in range(len(self.patterns)):
                    new_time = ts_to_hrs(short_data['time'], self.filedate)
                    setattr(self, '_time'+str(i), append(getattr(self, '_time'+str(i)), new_time))
            else: 
                self._time = t
                
            #return if only one file found because interpolator code will break
            if len(self._time)<2:
                print('Not enough times found with given file prefix.')
                return   
                
            #register interpolators for each variable
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
                    else:
                        new_data = variables[varname]['data']
                    self.variables[varname] = dict(units = variables[varname]['units'], data = new_data)
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
                    else:
                        new_data = variables[varname]['data']
                    self.variables[varname] = dict(units = variables[varname]['units'], data = new_data)
                    self.register_4D_variable(self.variables[varname]['units'], 
                                          self.variables[varname]['data'], varname,
                                          gridded_int)
            if verbose: print(f'Took {perf_counter()-t_reg:.5f}s to register '+\
                              f'{len(varname_list)} variables.')
            if verbose: print(f'Took a total of {perf_counter()-t0:.5f}s to kamodofy '+\
                              f'{len(varname_list)} variables.') 
            
        #define and register a 3D variable-----------------------------------------
        def register_3D_variable(self, units, variable, varname, gridded_int):
            """Registers a 3d interpolator with 3d signature"""
            
            #determine which file the variable came from
            for i in range(len(self.patterns)-1,-1,-1):  #go in reverse to account for overwritten variables
                if varname in self.varfiles[str(i)]:
                    time = getattr(self, '_time'+str(i))
                    lat = getattr(self, '_lat'+str(i))  #get the correct coordinates
                    lon = getattr(self, '_lon'+str(i))
                    #print(varname, self.patterns[i])
                    break
            
            #define and register the interpolators
            xvec_dependencies = {'time':'hr','lon':'deg','lat':'deg'}
            self = regdef_3D_interpolators(self, units, variable, time, 
                                           lon, lat, varname, 
                                           xvec_dependencies, gridded_int)       
            return 
        
        #define and register a 4D variable -----------------------------------------
        def register_4D_variable(self, units, variable, varname, gridded_int):
            """Registers a 4d interpolator with 4d signature"""
    
            #determine which file the variable came from
            for i in range(len(self.patterns)-1,-1,-1):  #go in reverse to account for overwritten variables
                if varname in self.varfiles[str(i)]:
                    time = getattr(self, '_time'+str(i))
                    lat = getattr(self, '_lat'+str(i))  #get the correct coordinates
                    lon = getattr(self, '_lon'+str(i))
                    radius = getattr(self, '_radius'+str(i))
                    #print(varname, self.patterns[i])
                    break
                else: print(varname, self.patterns[i])
            
            #define and register the interpolators
            xvec_dependencies = {'time':'hr','lon':'deg','lat':'deg','radius':'R_E'}
            self = regdef_4D_interpolators(self, units, variable, time,
                                              lon, lat, radius, varname, 
                                              xvec_dependencies, gridded_int)
            return
    return MODEL