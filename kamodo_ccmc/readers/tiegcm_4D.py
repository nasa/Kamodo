'''
TIEGCM Kamodo reader, adapted to new structure for satellite flythrough software
Initial version - Asher Pembroke (?)
Initial version of model_varnames contributed by Zachary Waldron
New code: Rebecca Ringuette (June 2021 and on)

NOTE: The current logic for variables that depend on imlev slices off self._imlev coordinate
    This only works because there is one variable that depends on imlev: H_imlev
    The logic on lines 311-313 will have to be reworked a bit if other variables depend on imlev later.

Remaining tasks:
    - check variable dictionary for proper naming, and use of units in Kamodo
    
    
'''
from numpy import vectorize
from datetime import datetime, timezone, timedelta


### Make a dict of the possible variable names in TIEGCM
###       the intended format is:   "Output VarName":['Latex_Varname', 'Latex_Unit' ]

###  The constituent species are output in units of mass mixing ratio (mmr).
###  Denoted by \psi_i, mmr is the fractional contribution of
###    a species i to the total mass density \rho_{total}, and 
###    is calculated as \psi_i = \rho_i / \rho_{total}, where 
###    \rho_i is the mass density of a constituent species.
model_varnames={
                 ### 4D Variables, vertical coordinate on midpoint levels (lev)
                 "ZGMID"    : ["H_ilev",'height dependent on primary pressure level',0,'GDZ','sph',['time','lon','lat','ilev'],"cm"],    # geometric height- interpolated to the mid points
                 "TN"       : ["T_n",'neutral temperature',1,'GDZ','sph',['time','lon','lat','ilev'],"K"],          # neutral temperature    
                 "O2"       : ["mmr_O2",'mass mixing ratio of molecular oxygen',2,'GDZ','sph',['time','lon','lat','ilev'],""],        # molecular oxygen,   mmr
                 "O1"       : ["mmr_O",'mass mixing ratio of atomic oxygen',3,'GDZ','sph',['time','lon','lat','ilev'],""],        # atomic oxygen ,   mmr
                 "N2"       : ["mmr_N2",'mass mixing ratio of molecular nitrogen',4,'GDZ','sph',['time','lon','lat','ilev'],""],        # molecular nitrogen,mmr
                 "HE"       : ["mmr_He",'mass mixing ratio of atomic helium',5,'GDZ','sph',['time','lon','lat','ilev'],""],         # helium  ,   mmr
                 "NO"       : ["mmr_NO",'mass mixing ratio of molecular nitric oxide',6,'GDZ','sph',['time','lon','lat','ilev'],""],         # nitric oxide , mmr
                 "N4S"      : ["mmr_Nstate4S",'mass mixing ratio of atomic nitrogen (4S state)',7,'GDZ','sph',['time','lon','lat','ilev'],""],        #  N4S ?,mmr
                 "N2D"      : ["mmr_Nstate2D",'mass mixing ratio of atomic nitrogen (2D state)',8,'GDZ','sph',['time','lon','lat','ilev'],""],     # N(2D) mmr
                 "TE"  : ["T_e",'electron temperature',9,'GDZ','sph',['time','lon','lat','ilev'],"K"],         #  ELECTRON TEMPERATURE,
                 "TI"  : ["T_i",'ion temperature',10,'GDZ','sph',['time','lon','lat','ilev'],"K"],         #  ION TEMPERATURE
                 "O2P" : ["N_O2plus",'number density of molecular oxygen ion',11,'GDZ','sph',['time','lon','lat','ilev'],"1/cm**3"],  #  O2+ ION
                 "OP"  : ["N_Oplus",'number density of atomic oxygen ion',12,'GDZ','sph',['time','lon','lat','ilev'],"1/cm**3"],    #   O+ ION
                 "N2N"      : ["N_N2",'number density of molecular nitrogen',13,'GDZ','sph',['time','lon','lat','ilev'],"1/cm**3"],    # molecular nitrogen (maybe number density),mmr
                 "CO2_COOL" : ["Q_CO2cool",'cooling rate of carbon dioxide',14,'GDZ','sph',['time','lon','lat','ilev'],"erg/g/s"],  #  CO2 cooling rates
                 "NO_COOL"  : ["Q_NOcool",'cooling rate of nitric oxide',15,'GDZ','sph',['time','lon','lat','ilev'],"erg/g/s"],      #  NO cooling rates
                 "UN"  : ["v_neast",'zonal neutral wind velocity (east)',16,'GDZ','sph',['time','lon','lat','ilev'],"cm/s"],            #  neutral ZONAL wind (+EAST)
                 "VN"  : ["v_nnorth",'meridional neutral wind velocity (north)',17,'GDZ','sph',['time','lon','lat','ilev'],"cm/s"],            #  neutral MERIDIONAL wind (+NORTH)
                 #"O2P_ELD"  : ['O2P_ELD','???',18,'GDZ','sph',['time','lon','lat','ilev'],''],     #NO DESCRIPTION GIVEN
                 #"N2P_ELD"  :['N2P_ELD','???',19,'GDZ','sph',['time','lon','lat','ilev'],''],      #NO DESCRIPTION GIVEN
                 "NPLUS"    :['N_Nplus','number density of atomic nitrogen ion',20,'GDZ','sph',['time','lon','lat','ilev'],'1/cm**3'],  #GUESS ONLY based on other number densities
                 #"NOP_ELD"  :['NOP_ELD','???',21,'GDZ','sph',['time','lon','lat','ilev'],''],    #NO DESCRIPTION GIVEN
                 "SIGMA_PED" :['sigma_P','Pedersen conductivity',22,'GDZ','sph',['time','lon','lat','ilev'],'S/m'],  #Pedersen Conductivity
                 "SIGMA_HAL" :['sigma_H','Hall conductivity',23,'GDZ','sph',['time','lon','lat','ilev'],'S/m'], #Hall Conductivity
                 "QJOULE"   :['Q_Joule','joule heating',24,'GDZ','sph',['time','lon','lat','ilev'],'erg/g/s'], #Joule Heating
                 "O_N2"    :['OtoN2','Oxygen/molecular nitrogen ratio',25,'GDZ','sph',['time','lon','lat','ilev'],''],  #O/N2 RATIO
                 #"N2D_ELD"  :['N2D_ELD','???',26,'GDZ','sph',['time','lon','lat','ilev'],''],  #NO DESCRIPTION GIVEN
                 "O2N"    :['N_O2','number density of molecular oxygen',27,'GDZ','sph',['time','lon','lat','ilev'],'1/cm**3'],  #GUESS ONLY 
                 #
                ### 4D Variables, vertical coordinate on interface levels (ilev)
                 "DEN"      :["rho",'total mass density',28,'GDZ','sph',['time','lon','lat','ilev1'],"g/cm**3"],     # total neutral mass density  
                 "ZG"       :["H_ilev1",'height dependent on secondary pressure level',29,'GDZ','sph',['time','lon','lat','ilev1'],"cm"],          # geometric height  
                 "Z"        :["H_geopot",'geopotential height',30,'GDZ','sph',['time','lon','lat','ilev1'],"cm"],            # geopotential height (cm)  
                 "NE"       : ["N_e",'electron number density',31,'GDZ','sph',['time','lon','lat','ilev1'],"1/cm**3"],    #  ELECTRON DENSITY
                 "OMEGA"    : ["omega",'Vertical motion frequency',32,'GDZ','sph',['time','lon','lat','ilev1'],"1/s"],      #  VERTICAL MOTION
                 "POTEN"    : ["V",'electric potential',33,'GDZ','sph',['time','lon','lat','ilev1'],"V"],        #  ELECTRIC POTENTIAL
                 "UI_ExB"   : ["v_iExBeast",'zonal ExB ion velocity (east)',34,'GDZ','sph',['time','lon','lat','ilev1'],'cm/s'],  #Zonal ExB Velocity
                 "VI_ExB"   :["v_iExBnorth",'meridional ExB ion velocity (north)',35,'GDZ','sph',['time','lon','lat','ilev1'],'cm/s'],  #Meridional ExB Velocity
                 "WI_ExB"   :["v_iExBup",'vertical ExB ion velocity (up)',36,'GDZ','sph',['time','lon','lat','ilev1'], 'cm/s'], #Vertical ExB Velocity
                ### 4D Variables, vertical coordinate on interface mag levels (imlev)
                 "ZMAG"  : ["H_milev",'height dependent on geomagnetic pressure level',37,'MAG','sph',['time','mlon','mlat','milev'],"km"],     #  Geopotential Height on Geomagnetic Grid
                #
                ### 3D Variables,    (time, lat, lon)
                 "TEC"  : ["TEC",'vertical total electron content (height integrated from bottom to top boundary)',38,'GDZ','sph',['time','lon','lat'],"1/cm**2"],     #  Total Electron Content
                 #"TLBC"  : ["T_nLBC",'Lower boundary condition for T_n',39,'GDZ','sph',['time','lon','lat'],"K"],       #  Lower boundary condition for TN
                 #"ULBC"  : ["v_neastLBC",'Lower boundary condition for v_n east component',40,'GDZ','sph',['time','lon','lat'],"cm/s"],    #  Lower boundary condition for UN
                 #"VLBC"  : ["v_nnorthLBC",'Lower boundary condition for v_n north component',41,'GDZ','sph',['time','lon','lat'],"cm/s"],    #  Lower boundary condition for VN
                 #"TLBC_NM"  : ["T_nLBCNminus1",'Lower boundary condition for T_n at t=N-1',42,'GDZ','sph',['time','lon','lat'],"K"],  #  Lower boundary condition for TN (TIME N-1)
                 #"ULBC_NM"  : ["v_neastLBCNminus1",'Lower boundary condition for v_n north component at t=N-1',43,'GDZ','sph',['time','lon','lat'],"cm/s"],  #  Lower boundary condition for UN (TIME N-1)
                 #"VLBC_NM"  : ["v_nnorthLBCNminus1",'Lower boundary condition for v_n east component at t=N-1',44,'GDZ','sph',['time','lon','lat'],"cm/s"],  #  Lower boundary condition for VN (TIME N-1)
                 "QJOULE_INTEG":["W_JouleH",'height integrated joule heating',45,'GDZ','sph',['time','lon','lat'],'erg/cm**2/s'],  #Height-integrated Joule Heating
                 "EFLUX"  :['Phi_E','energy flux',46,'GDZ','sph',['time','lon','lat'],'erg/cm**2/s'],  #Aurora Energy Flux
                 "HMF2" :['HmF2','height of maximum electron number density in F2 layer',47,'GDZ','sph',['time','lon','lat'],'km'], #  Height of the F2 Layer
                 "NMF2"  :['NmF2','maximum electron number density in F2 layer',48,'GDZ','sph',['time','lon','lat'],'1/cm**3'], #Peak Density of the F2 Layer
                  }

#####--------------------------------------------------------------------------------------
##### Define some helpful functions for dealing with time systems

def dts_to_ts(file_dts):
    '''Get datetime timestamp in UTC from datetime string'''    
    
    return datetime.timestamp(datetime.strptime(file_dts, '%Y-%m-%d %H:%M:%S'
                                                ).replace(tzinfo=timezone.utc))

def year_mtime_todt0(year, mtime):  #self.filedate
    '''Convert year and day to datetime object in UTC at midnight'''
    
    day, hour, minute = mtime  #unpack mtime values
    return datetime(int(year),1,1).replace(tzinfo=timezone.utc)+\
        timedelta(days=int(day-1))
        
def year_mtime_todt(year, mtime):
    '''Convert year and [day,hour,minute] to datetime object in UTC'''
    
    day, hour, minute = mtime  #unpack mtime values
    return datetime(int(year),1,1).replace(tzinfo=timezone.utc)+\
        timedelta(days=int(day-1),hours=int(hour),minutes=int(minute))        
        
def year_mtime_todts(year, mtime):
    '''Convert year and mtime to a datetime string'''
    
    return datetime.strftime(year_mtime_todt(year, mtime), '%Y-%m-%d %H:%M:%S')

def year_mtime_todate(year, mtime):
    '''Use year and mtime to determine the date in the file. Returns a datetime object.'''
    
    date_string = datetime.strftime(year_mtime_todt(year, mtime), '%Y-%m-%d')  #'YYYY-MM-DD'
    return datetime.strptime(date_string, '%Y-%m-%d').replace(tzinfo=timezone.utc)

@vectorize
def year_mtime_tohrs(year, day, hour, minute, filedate):
    '''Convert year and mtime to hours since midnight using predetermined datetime object.'''
    
    mtime = [day, hour, minute]
    return (year_mtime_todt(year, mtime)-filedate).total_seconds()/3600.

def ts_to_hrs(time_val, filedate):
    '''Convert utc timestamp to hours since midnight on filedate.'''
    
    return (datetime.utcfromtimestamp(time_val).replace(tzinfo=timezone.utc)-filedate).total_seconds()/3600.


def MODEL():
    from time import perf_counter
    from os.path import basename
    from numpy import zeros, transpose, array, append, insert, where, unique 
    from numpy import NaN, diff, abs, mean, broadcast_to, cos, sin, sum
    from numpy import pi as nppi
    from netCDF4 import Dataset
    from kamodo import Kamodo
    #print('KAMODO IMPORTED!')
    from kamodo_ccmc.readers.reader_utilities import regdef_3D_interpolators, regdef_4D_interpolators    
    
    class MODEL(Kamodo): 
        '''TIEGCM model data reader.'''
        def __init__(self, full_filename, variables_requested=[], runname="noname",
                     filetime=False, verbose=False, gridded_int=True, printfiles=False,
                     fulltime=True, **kwargs):  #filename should include the full path
            
            #### Use a super init so that your class inherits any methods from Kamodo
            super(MODEL, self).__init__()
    
            #store time information for satellite flythrough layer to choose the right file
            t0 = perf_counter()
            filename = basename(full_filename)
            file_dir = full_filename.split(filename)[0]            
            cdf_data = Dataset(full_filename, 'r')
            
            #calculate time information
            year = array(cdf_data.variables['year'])
            mtime = array(cdf_data.variables['mtime'])
            day, hour, minute = mtime.T #only matters for the vectorized function
            self.filedate = year_mtime_todt0(year[0], mtime[0])  #datetime object for file date at midnight UTC
            self.datetimes = [year_mtime_todts(y, m) for y, m \
                              in zip([year[0], year[-1]],[mtime[0],mtime[-1]])]  #strings in format = YYYY-MM-DD HH:MM:SS
            self.filetimes=[dts_to_ts(file_dts) for file_dts in self.datetimes]   #timestamps in UTC 
            time = year_mtime_tohrs(year, day, hour, minute, self.filedate)
            #time = array([year_mtime_tohrs(y, m, self.filedate) for y, m in \
            #              zip(year, mtime)])  #hours since midnight of self.filedate
            self.dt = diff(time).max()*3600.  #time is in hours
            
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
                
                file_pattern = file_dir+'s*.nc' #returns a string for tiegcm
                files = sorted(glob(file_pattern))
                filenames = unique([basename(f) for f in files])
                
                #find closest file by utc timestamp
                #tiegcm has an open time at the beginning, so need an end time from the previous file
                #files are automatically sorted by YYMMDD, so previous file is previous in the list
                current_idx = where(filenames==filename)[0]
                if current_idx==0:
                    print('No earlier file available.')
                    filecheck = False  
                    if filetime:
                        return   
                else:
                    min_filename = file_dir+filenames[current_idx-1][0]  #-1 for adding a beginning time
                    kamodo_test = MODEL(min_filename, filetime=True, fulltime=False)   
                    time_test = abs(kamodo_test.filetimes[1]-self.filetimes[0])  
                    if time_test<=self.dt:  #if nearest file time at least within one timestep (hrs)
                        filecheck = True                
                    
                        #time only version if returning time for searching
                        if filetime:
                            kamodo_neighbor = MODEL(min_filename, fulltime=False, filetime=True)
                            self.datetimes[0] = kamodo_neighbor.datetimes[1]
                            self.filetimes[0] = kamodo_neighbor.filetimes[1]
                            return  #return object with additional time (for SF code) 
                        
                        #get kamodo object with same requested variables to add to each array below
                        if verbose: print(f'Took {perf_counter()-t0:.3f}s to find closest file.')
                        kamodo_neighbor = MODEL(min_filename, variables_requested=variables_requested, 
                                               fulltime=False)
                        self.datetimes[0] = kamodo_neighbor.datetimes[1]
                        self.filetimes[0] = kamodo_neighbor.filetimes[1]
                        short_data = kamodo_neighbor.short_data
                        if verbose: print(f'Took {perf_counter()-t0:.3f}s to get data from previous file.')
                    else:
                        print(f'No earlier file found within {self.dt:.1f}s')
                        filecheck = False
                        if filetime:
                            return                    
    
            #These lists need to be the standardized variable name to match that above,
            #not the names from the data file.
            self.ilev1_list  = [value[0] for key, value in model_varnames.items() if value[5][-1]=='ilev1']
            self.ilev_list = [value[0] for key, value in model_varnames.items() if value[5][-1]=='ilev']
            self.milev_list = [value[0] for key, value in model_varnames.items() if value[5][-1]=='milev']
            #don't need an internal coord dict because there is only one lat/lon (other than magnetic)
            
            #perform initial check on variables_requested list
            if len(variables_requested)>0 and fulltime and variables_requested!='all':
                test_list = [value[0] for key, value in model_varnames.items()]
                err_list = [item for item in variables_requested if item not in test_list]
                if len(err_list)>0: print('Variable name(s) not recognized:', err_list)
                
            #translate from standardized variables to names in file
            #remove variables requested that are not in the file
            if len(variables_requested)>0 and variables_requested!='all':
                gvar_list = [key for key, value in model_varnames.items() \
                                 if value[0] in variables_requested and \
                                     key in cdf_data.variables.keys()]  # file variable names
                
                #check for variables requested but not available
                if len(gvar_list)!=len(variables_requested):
                    err_list = [value[0] for key, value in model_varnames.items() \
                                 if value[0] in variables_requested and \
                                     key not in cdf_data.variables.keys()]
                    if len(err_list)>0: print('Some requested variables are not available:', err_list)
                
                #check that the appropriate height variable is added for the variables requested
                check_list = [key for key, value in model_varnames.items()\
                              if value[0] in self.ilev1_list and key in gvar_list]
                if 'ZG' not in gvar_list and len(check_list)>0: 
                    gvar_list.append('ZG')  #force addition of H for conversion of ilev to H and back
                check_list = [key for key, value in model_varnames.items()\
                              if value[0] in self.ilev_list and key in gvar_list]
                if 'ZGMID' not in gvar_list and len(check_list)>0: 
                    if 'ZGMID' in cdf_data.variables.keys():
                        gvar_list.append('ZGMID')  #no ZGMID in files from CCMC!
                    else: 
                        print('ERROR! H_ilev function not in data!')
                check_list = [key for key, value in model_varnames.items()\
                              if value[0] in self.milev_list and key in gvar_list]    
                if 'ZMAG' not in gvar_list and len(check_list)>0: 
                    gvar_list.append('ZMAG')
            else:  #only input variables on the avoid_list if specifically requested
                avoid_list = ['TLBC','ULBC','VLBC','TLBC_NM','ULBC_NM','VLBC_NM',
                              'NOP_ELD','O2P_ELD','N2P_ELD','N2D_ELD'] 
                gvar_list = [key for key in cdf_data.variables.keys() \
                             if key in model_varnames.keys() and \
                                 key not in avoid_list]
                if not fulltime and variables_requested=='all':
                    self.var_dict = {value[0]: value[1:] for key, value in model_varnames.items() \
                            if key in gvar_list}
                    return                     

            # Store the requested variables into a dictionary 
            variables = {model_varnames[key][0]:{'units':model_varnames[key][-1], 
                               'data':array(cdf_data.variables[key])}\
                          for key in gvar_list}  #store with key = standardized name
            
            #prepare and return data only for last timestamp
            if not fulltime:  
                cdf_data.close()
                variables['time'] = self.filetimes[1]  #utc timestamp
                self.short_data = variables
                return
    
            #### Store our inputs as class attributes to the class
            self.filename      = full_filename
            self.runname       = runname
            self.missing_value = NaN
            self._registered   = 0
            self.variables     = dict()
            self.modelname = 'TIEGCM'
            if printfiles: 
                print('Files:', self.filename)
            
            #### Store coordinate data as class attributes   
            if filecheck:  #new_time iis a utc timestamp
                new_time = ts_to_hrs(short_data['time'], self.filedate)  #new time in hours since midnight
                self._time = insert(time, 0, new_time)  #insert new value
            else: 
                self._time = time
                
            #store coordinates
            lat = array(cdf_data.variables['lat'])  #NOT FULL RANGE IN LATITIUDE!!!
            lat = insert(lat, 0, -90)  #insert a grid point at beginning (before -87.5)
            self._lat = append(lat, 90.)   #and at the end (after 87.5)
            lon = array(cdf_data.variables['lon'])  #NOT WRAPPED IN LONGITUDE!!!!!
            self._lon = append(lon, 180.)  #add 180. to end of array            
            self._ilev = array(cdf_data.variables['lev'])
            self._ilev1 = array(cdf_data.variables['ilev'])
            self._milev = array(cdf_data.variables['imlev'])  #'imlev' isn't used by any of the variables except H_imlev
            self._mlat = array(cdf_data.variables['mlat'])
            self._mlon = array(cdf_data.variables['mlon'])  #-180 to 180  
    
            #close file
            cdf_data.close()
            if verbose: print(f'Took {perf_counter()-t0:.6f}s to read in data')
    
            # register interpolators for each requested variable
            varname_list, self.variables = [key for key in variables.keys()], {}  #store original list b/c gridded interpolators
            t_reg = perf_counter()
            for varname in varname_list:
                #print(varname, len(variables[varname]['data'].shape))
                if len(variables[varname]['data'].shape)==3:
                    if filecheck:  #if neighbor found
                        #append data for first time stamp, transpose and register
                        data_shape = list(variables[varname]['data'].shape)
                        data_shape[0]+=1  #add space for time
                        new_data = zeros(data_shape)                
                        new_data[1:,:,:] = variables[varname]['data']  #put in current data
                        new_data[0,:,:] = short_data[varname]['data'][-1,:,:]  #add in data for additional time
                        variable = transpose(new_data, (0,2,1)) #(t,lat,lon) -> (t,lon,lat)
                    else:
                        variable = transpose(variables[varname]['data'], (0,2,1)) #(t,lat,lon) -> (t,lon,lat)
                    self.variables[varname] = dict(units = variables[varname]['units'], data = variable)
                    self.register_3D_variable(self.variables[varname]['units'], 
                                          self.variables[varname]['data'], varname,
                                          gridded_int)
                elif len(variables[varname]['data'].shape)==4:
                    if filecheck:
                        #append data for first time stamp, transpose and register
                        data_shape = list(variables[varname]['data'].shape)
                        data_shape[0]+=1  #add space for time
                        new_data = zeros(data_shape)
                        new_data[1:,:,:,:] = variables[varname]['data']  #put in current data
                        new_data[0,:,:,:] = short_data[varname]['data'][-1,:,:,:]   #add in data for additional time                
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
  
        def wrap_3Dlatlon(self, varname, variable):
            '''Wraps the data array in longitude (-180=180), and latitude'''
        
            shape_list = list(variable.shape)  #e.g. time, lat, lon -> time, lon, lat!!!
            shape_list[2]+=2  #need two more places in latitude
            shape_list[1]+=1  #need one more place in longitude
            tmp_arr = zeros(shape_list)  #array to set-up wrapped data in
            tmp_arr[0:,:-1,1:-1]=variable  #copy data into grid
            
            #wrapping in latitude for scalar variables
            # put in top values
            top = mean(tmp_arr[:,:-1,1],axis=1)  #same shape as time axis
            tmp_arr[:,:-1,0] = broadcast_to(top, (shape_list[1]-1,shape_list[0])).T
            #same for bottom, reusing variable names
            top = mean(tmp_arr[:,:-1,-2],axis=1)  #same shape as time axis
            tmp_arr[:,:-1,-1] = broadcast_to(top, (shape_list[1]-1,shape_list[0])).T
            
            #wrap in longitude after to prevent double counting in average
            tmp_arr[:,-1,:] = tmp_arr[:,0,:]  
            self.variables[varname]['data'] = tmp_arr  #store result
            return tmp_arr        
    
        def wrap_4Dlatlon(self, varname, variable):
            '''Wraps the data array in longitude (-180=180), and latitude (0=-2, -1=1)'''
        
            shape_list = list(variable.shape)  #time, lon, lat, ilev
            shape_list[2]+=2  #need two more places in latitude
            shape_list[1]+=1  #need one more place in longitude
            tmp_arr = zeros(shape_list)  #array to set-up wrapped data in
            tmp_arr[:,:-1,1:-1,:] = variable  #copy data into grid
            
            #wrapping in latitude for scalar and cartesian/radial variables
            if varname not in ['u_n','v_n','u_iExB','v_iExB']:
                #print('Calculating scalar sum for', varname)
                # put in top values
                top = mean(tmp_arr[:,:-1,1,:],axis=1)  #average over longitudes
                tmp_arr[:,:-1,0,:] = transpose(broadcast_to(top, (shape_list[1]-1,shape_list[0],
                                                       shape_list[3])), (1,0,2))
                #same for bottom, reusing variable names
                top = mean(tmp_arr[:,:-1,-2,:],axis=1)  #average over longitudes
                tmp_arr[:,:-1,-1,:] = transpose(broadcast_to(top, (shape_list[1]-1,shape_list[0],
                                                      shape_list[3])), (1,0,2))
            #wrapping in latitude for relevant vector variables
            elif varname in ['u_n','v_n','u_iExB','v_iExB']:
                #print('Calculating vector sum for', varname)
                #calculate net vector magnitude for top
                tmp_arr[:,:-1,0,:] = self.vector_average4D(tmp_arr[:,:-1,1,:], 
                                                      shape_list, varname, self._lat[0])
                #repeat for bottom
                tmp_arr[:,:-1,-1,:] = self.vector_average4D(tmp_arr[:,:-1,-2,:],
                                                      shape_list, varname, self._lat[-1])
            tmp_arr[:,-1,:,:] = tmp_arr[:,0,:,:]  #wrap value in longitude                
            self.variables[varname]['data'] = tmp_arr  #store result
            return tmp_arr

        def vector_average4D(self, top, shape_list, varname, latval):
            '''find vector average at pole for array with shape (time, lon, height)'''
        
            #find net x and y components, final array shapes are time, lon, and height/ilev
            lon_arr = transpose(broadcast_to(self._lon[:-1], (shape_list[0],
                                        shape_list[3], shape_list[1]-1)), (0,2,1))      
            #need to put 'old' shape at end in broadcast_to call ^
            xval = sum(top*cos((lon_arr+180.)*nppi/180.), axis=1)  #same shape as
            yval = sum(top*sin((lon_arr+180.)*nppi/180.), axis=1)  #time and vertical 
            xarr = transpose(broadcast_to(xval, (shape_list[1]-1,shape_list[0], 
                                                 shape_list[3])), (1,0,2))
            yarr = transpose(broadcast_to(yval, (shape_list[1]-1,shape_list[0], 
                                                 shape_list[3])), (1,0,2))

            #convert to proper unit vector (see wiki on spherical coordinates)
            if 'u' in varname:  #Zonal / east components -> convert to psi_hat vector (longitude)
                # -xsin(psi)+ycos(psi), psi = longitude (0 to 360)
                new_top = -xarr*sin((lon_arr+180.)*nppi/180.)+yarr*cos((lon_arr+180.)*nppi/180.)
            elif 'v' in varname:  #meridional/north -> convert to theta_hat vector (latitude)
                # xcos(psi)cos(theta)+ysin(psi)cos(theta),  sin(theta) is always zero at the poles
                # theta = latitude (0 to 180), psi = longitude (0 to 360)
                new_top = xarr*cos((lon_arr+180.)*nppi/180.)*cos((90.-latval)*nppi/180.)+\
                    yarr*sin((lon_arr+180.)*nppi/180.)*cos((90.-latval)*nppi/180.)
            
            #flip around so values match longitude location (to keep zero reference point true)
            zero_idx = min(where(self._lon>=0.)[0])  #find splitting index
            top = zeros((shape_list[0],shape_list[1]-1,shape_list[3]))  #empty array for destination
            div = (shape_list[1]-1)%2  #deal with even or odd number of non-wrapped longitude elements
            #print(shape_list[1]-1, zero_idx, div)
            top[:,:zero_idx-div,:] = new_top[:,zero_idx:,:]  #move last half to first half
            top[:,zero_idx-div:,:] = new_top[:,:zero_idx,:]  #and vice versa
            return top
        
        ##### Define and register a 3D variable -----------------------------------------
        def register_3D_variable(self, units, variable, varname, gridded_int):
            """Registers a 3d interpolator with 3d signature"""
            
            #define and register the interpolators
            xvec_dependencies = {'time':'hr','lon':'deg','lat':'deg'}
            wrapped_data = self.wrap_3Dlatlon(varname, variable)
            self = regdef_3D_interpolators(self, units, wrapped_data, self._time, 
                                           self._lon, self._lat, varname, 
                                           xvec_dependencies, gridded_int)       
            return 
        
        
        #### Define and register a 4D variable -----------------------------------------
        def register_4D_variable(self, units, variable, varname, gridded_int):
            """Registers a 4d interpolator with 4d signature"""
                 
            ####  Get the correct coordinates
            if varname in self.ilev1_list:
                h = self._ilev1
                coord_lat, coord_lon = self._lat, self._lon
                xvec_dependencies = {'time':'hr','lon':'deg','lat':'deg','ilev1':'m/m'}
            elif varname in self.ilev_list:
                h = self._ilev
                coord_lat, coord_lon = self._lat, self._lon
                xvec_dependencies = {'time':'hr','lon':'deg','lat':'deg','ilev':'m/m'}
            elif varname in self.milev_list:
                h = self._milev
                coord_lat, coord_lon = self._mlat, self._mlon
                xvec_dependencies = {'time':'hr','mlon':'deg','mlat':'deg','milev':'m/m'}
            else:
                print(varname, 'error')
    
            #### define and register the interpolators
            if 'lat' in xvec_dependencies.keys():
                wrapped_data = self.wrap_4Dlatlon(varname, variable)
            else: 
                top_shape = list(variable[:,:,:,-1].shape)
                top_size = top_shape[0]*top_shape[1]*top_shape[2]  #3D array
                idx_top = where(variable[:,:,:,-1]>1e+35)[0]
                tmp_data = variable
                while top_size==len(idx_top):   #then top row is undefined! Remove it.
                    print(f'All values at max milev are 1e+36 for {varname}. Slicing off top array.')
                    if self._milev.shape[0]==len(tmp_data[0,0,0,:]):
                        self._milev = self._milev[0:-1]
                    tmp_data = tmp_data[:,:,:,0:-1]
                    top_shape = list(tmp_data[:,:,:,-1].shape)
                    top_size = top_shape[0]*top_shape[1]*top_shape[2]  #3D array
                    idx_top = where(tmp_data[:,:,:,-1]>1e+35)[0]
                wrapped_data = tmp_data
                h = self._milev
            self = regdef_4D_interpolators(self, units, wrapped_data, 
                                              self._time, coord_lon, coord_lat, h,
                                              varname, xvec_dependencies, gridded_int)
            return
    return MODEL
