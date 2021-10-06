'''
Original version: Lutz Raestatter Oct 1(?), 2021
Modify to work with flythrough: Oct 5, 2021 (Rebecca Ringuette)
'''
from datetime import datetime,timedelta,timezone

#standard model dictionary for reference
model_varnames={
    ### 3D variables - to be aggregated into 4D ###
    'bx':['B_x','x component of magnetic field',0,'GSE','car',['time','x','y','z'],'nT'],
    'by':['B_y','y component of magnetic field',1,'GSE','car',['time','x','y','z'],'nT'],
    'bz':['B_z','z component of magnetic field',2,'GSE','car',['time','x','y','z'],'nT'],
    'bx1':['B1_x','x component of magnetic field (on grid cell faces)',3,'GSE','car',['time','x','x','x'],'nT'],
    'by1':['B1_y','y component of magnetic field (on grid cell faces)',4,'GSE','car',['time','y','y','y'],'nT'],
    'bz1':['B1_z','z component of magnetic field (on grid cell faces)',5,'GSE','car',['time','z','z','z'],'nT'],
    'ex':['E_x','x component of electric field (on grid cell edges)',6,'GSE','car',['time','x','x','x'],'mV/m'],
    'ey':['E_y','y component of electric field (on grid cell edges)',7,'GSE','car',['time','y','y','y'],'mV/m'],
    'ez':['E_z','z component of electric field (on grid cell edges)',8,'GSE','car',['time','z','z','z'],'mV/m'],
    'vx':['V_x','x component of plasma velocity',9,'GSE','car',['time','x','y','z'],'km/s'],
    'vy':['V_y','y component of plasma velocity',10,'GSE','car',['time','x','y','z'],'km/s'],
    'vz':['V_z','z component of plasma velocity',11,'GSE','car',['time','x','y','z'],'km/s'],
    'rr':['N_plasma','plasma number denstity (hydrogen equivalent)',12,'GSE','car',['time','x','y','z'],'1/cm**3'],
    'resis':['eta','resistivity',13,'GSE','car',['time','x','y','z'],'m**2/s'],
    'pp':['P_plasma','plasma pressure',14,'GSE','car',['time','x','y','z'],'pPa'],
}

# variable linkage to grid position vectors are established during variable registration
# these are gx_bx, gy_bx, ... gz_ez affecting magnetic field (b1x,b1y,b1z) and electric field (ex,ey,ez)

#convert an array of timestamps to an array of hrs since midnight
def ts_to_hrs(time_val, filedate):
    '''Convert utc timestamp to hours since midnight on filedate.'''
    
    return (datetime.utcfromtimestamp(time_val).replace(tzinfo=timezone.utc)-filedate).total_seconds()/3600.

def hrs_to_ts(hrs, filedate):
    '''Add hours to filedate and return utc timestamp.'''
    
    return datetime.timestamp(filedate+timedelta(hours=float(hrs)))

#main function/class definition
def MODEL():
    from numpy import array, zeros, abs, NaN, unique, diff, where, append, sqrt, sum
    from numpy import float32
    from time import perf_counter
    from netCDF4 import Dataset
    from os.path import isfile, basename
    from glob import glob
    from kamodo import Kamodo, kamodofy
    #from kamodo.readers.reader_utilities import regdef_4D_interpolators
    from kamodo.readers.reader_utilities import register_interpolator, define_4d_gridded_interpolator
    from scipy.interpolate import RegularGridInterpolator
    #print('KAMODO IMPORTED!')
    
    class MODEL(Kamodo):
        '''OpenGGCM_GM magnetosphere reader'''
        def __init__(self,full_file_prefix, variables_requested=[], runname = "noname",
                     filetime=False, verbose=False, gridded_int=True, printfiles=False, 
                     fulltime=True, missing_value=NaN, **kwargs):
            super(MODEL, self).__init__()
            t0=perf_counter() # profiling time stamp
            
            #separate file directory from file name
            file_prefix = basename(full_file_prefix)  # runname.3df.ssss
            file_dir = full_file_prefix.split(file_prefix)[0]
            
            #convert files to netcdf4 if needed
            nc_file = full_file_prefix+'.nc'  # input file name: file_dir/YYYY-MM-DD_HH.nc
            if isfile(nc_file):  #file already prepared!
                self.conversion_test = True  #default value
            else:  #file not prepared, prepare it
                try:  #I don't have the file converter, so leave in try/except for now
                    from openggcm_to_cdf import openggcm_combine_magnetosphere_files as gmconv
                    self.conversion_test = gmconv(full_file_prefix)
                    #should return a boolean (True is successful, False if not)
                except:
                    self.conversion_test = False
            day_flag=False  #files are too big to do more than one hour (~40 GB for one hour)

            #establish time attributes first
            cdf_data = Dataset(nc_file, 'r')
            self.filedate = datetime.strptime(cdf_data.filedate+' 00:00:00', 
                                              '%Y-%m-%d %H:%M:%S').replace(tzinfo=timezone.utc)
            t = array(cdf_data.variables['_time'])  #hours since midnight
            if len(t)>1: self.dt = diff(t).max()*3600.  #t is in hours since midnight
            else: self.dt = 0
            #print("dt:",self.dt)  # this dt may not be constant            
            self.datetimes=[datetime.utcfromtimestamp(hrs_to_ts(t[0], self.filedate)).isoformat(sep=' '),
                            datetime.utcfromtimestamp(hrs_to_ts(t[-1], self.filedate)).isoformat(sep=' ')]
            self.filetimes=[hrs_to_ts(t[0], self.filedate), hrs_to_ts(t[-1], self.filedate)]   #timestamps for matching in wrapper
            #print(self.filedate, self.datetimes, self.filetimes, full_file_prefix)

            #execute logic for finding nearest time in neighboring file if requested
            if filetime and not fulltime: #(used when searching for neighboring files below)
                return  #return times as is to prevent recursion
            
            #if variables are given as integers, convert to standard names
            if len(variables_requested)>0:
                if isinstance(variables_requested[0], int):
                    print('Integers detected. Converting...', end="")
                    tmp_var = [value[0] for key, value in model_varnames.items()\
                                           if value[2] in variables_requested]
                    variables_requested = tmp_var
                    print('Converted:', variables_requested)

            if fulltime:  #add boundary time for interp btwm files (default value of fulltime)
                t_search = perf_counter()
                #find other files with same pattern
                if day_flag:  #if entire day in one file (data too big for now)
                    files = sorted(glob(file_dir+file_prefix[:-10]+'*.nc')) #chop date off
                    file_prefixes = unique([basename(f).split('3df_')[0]+'3df_'+f.split('3df_')[1][:10]\
                                            for f in files])  #full day files only, no .nc
                else: #if data split into hourly files
                    files = sorted(glob(file_dir+file_prefix[:-13]+'*.nc'))
                    file_prefixes = unique([basename(f).split('3df_')[0]+'3df_'+f.split('3df_')[1][:13]\
                                            for f in files])  #hourly files only
            
                #find closest file by utc timestamp
                #openggcm_gm has an open time at the end, so need a beginning time from the closest file
                #files are automatically sorted by YYMMDD, so next file is next in the list
                current_idx = where(file_prefixes==file_prefix)[0]
                if current_idx+1==len(file_prefixes):
                    print('No later file available.')
                    filecheck = False  
                    if filetime:
                        return   
                else:
                    min_file_prefix = file_prefixes[current_idx+1][0]  #+1 for adding an end time
                    #print(min_file_prefix)
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
                                self.datetimes[1] = kamodo_neighbor.datetimes[0]  #add first time from next file
                                self.filetimes[1] = kamodo_neighbor.filetimes[0]
                                return  #return object with additional time (for SF code) 
                            
                            #get kamodo object with same requested variables to add to each array below
                            kamodo_neighbor = MODEL(file_dir+min_file_prefix, 
                                                    variables_requested=variables_requested, 
                                                    fulltime=False)
                            self.datetimes[1] = kamodo_neighbor.datetimes[0]
                            self.filetimes[1] = kamodo_neighbor.filetimes[0]
                            short_data = kamodo_neighbor.short_data     
                            del kamodo_neighbor  #to avoid too much data in memory
                            if verbose: print(f'Took {perf_counter()-t0:.3f}s to get data from closest file.')
                        else:
                            if verbose: print(f'No later file found within {self.dt:.1f}s.')
                            filecheck = False  
                            if filetime:
                                return           

            #perform initial check on variables_requested list
            if len(variables_requested)>0 and fulltime:
                test_list = [value[0] for key, value in model_varnames.items()]
                err_list = [item for item in variables_requested if item not in test_list]
                if len(err_list)>0: print('Variable name(s) not recognized:', err_list)
                
            #collect variable list            
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
            else:  #only input variables on the avoid_list if specifically requested
                avoid_list = []   #empty for now
                gvar_list = [key for key in cdf_data.variables.keys() \
                             if key in model_varnames.keys() and \
                                 key not in avoid_list]            

            # Store variable's data and units, transposing the 2D+time array.
            variables = {model_varnames[key][0]:{'units':model_varnames[key][-1],
                                   'data':array(cdf_data.variables[key])}\
                              for key in gvar_list} 
            if verbose: print('Done reading in variable data.', full_file_prefix)
                
            #prepare and return data with first utc timestamp
            if not fulltime:  
                cdf_data.close()
                #print(self.filetimes[0], full_file_prefix)
                variables['time'] = self.filetimes[0]
                self.short_data = variables
                return            
    
            #store variables
            self.near_Earth_boundary_radius = cdf_data.near_Earth_boundary_radius
            self.near_Earth_boundary_radius_unit = cdf_data.near_Earth_boundary_radius_units
            #print('Inner boundary radius:', self.near_Earth_boundary_radius, 
            #      self.near_Earth_boundary_radius_unit)  #should be 2.5 R_E
            self.missing_value = NaN
            self.verbose = verbose
            self.filename = cdf_data.file.split(',')
            self.modelname = cdf_data.model
            self.runname = runname
            self.modelname = 'OpenGGCM_GM'
            self._registered = 0
            if printfiles: 
                print('Files:')
                for file in self.filename: print(file)
    
            #### Store coordinate time data as class attributes   
            if filecheck:
                new_time = ts_to_hrs(short_data['time'], self.filedate)  #new time in hours since midnight
                self._time = append(t, new_time) 
            else: 
                self._time = t
            #print(self._time)

            #add coordinate grids as needed
            #grid_list = list of coordinate names in cdf file
            grid_list = ['_x','_y','_z','_x_bx','_y_bx','_z_bx','_x_by','_y_by',
                         '_z_by','_x_bz','_y_bz','_z_bz','_x_ex','_y_ex','_z_ex',
                         '_x_ey','_y_ey','_z_ey','_x_ez','_y_ez','_z_ez']
            #trim down to only save coordinate grids needed for variables requested
            #  and available in file
            if 'B1_x' not in variables.keys(): 
                grid_list.remove('_x_bx')
                grid_list.remove('_y_bx')
                grid_list.remove('_z_bx')
            if 'B1_y' not in variables.keys(): 
                grid_list.remove('_x_by')
                grid_list.remove('_y_by')
                grid_list.remove('_z_by')
            if 'B1_z' not in variables.keys(): 
                grid_list.remove('_x_bz')
                grid_list.remove('_y_bz')
                grid_list.remove('_z_bz')
            if 'E1_x' not in variables.keys(): 
                grid_list.remove('_x_ex')
                grid_list.remove('_y_ex')
                grid_list.remove('_z_ex')
            if 'E1_y' not in variables.keys(): 
                grid_list.remove('_x_ey')
                grid_list.remove('_y_ey')
                grid_list.remove('_z_ey')
            if 'E1_z' not in variables.keys(): 
                grid_list.remove('_x_ez')
                grid_list.remove('_y_ez')
                grid_list.remove('_z_ez')     
                
            #separate identical coordinate grid values by a small amount to 
            #  prevent interpolator from breaking
            for grid in grid_list:
                tol, data = 0.000001, array(cdf_data.variables[grid])
                if diff(data).min()==0.0: 
                    idx = where(diff(data)==0.0)[0]
                    while len(idx)>0:  #while identical #s remain
                        data[idx]-=tol
                        idx = where(diff(data)==0.0)[0]
                        tol+=tol  #slightly increase value subtracted from coordinate
                else:
                    data = array(cdf_data.variables[grid])
                if verbose: print(grid, tol)
                setattr(self, grid, data)  #store coordinate data
                
            cdf_data.close()
            if verbose: print(f'Took {perf_counter()-t0:.6f}s to read in data')
    
            #register interpolators for each variable
            varname_list, self.variables = [key for key in variables.keys()], {}  #store original list b/c gridded interpolators
            t_reg = perf_counter()
            for varname in varname_list:  #all are 3D variables
                if filecheck:  #if neighbor found
                    #append data for first time stamp and register
                    data_shape = list(variables[varname]['data'].shape)
                    data_shape[0]+=1  #add space for time
                    new_data = zeros(data_shape, dtype=float32)                
                    new_data[:-1,:,:,:] = variables[varname]['data']  #put in current data
                    new_data[-1,:,:,:] = short_data[varname]['data'][0,:,:,:]  #add in data for additional time
                else:
                    new_data = variables[varname]['data']
                self.variables[varname] = dict(units = variables[varname]['units'], 
                                               data = {})     #not saving data to decrease memory demand      
                self.register_variable(self.variables[varname]['units'], 
                                       new_data, varname, gridded_int)
            if verbose: print(f'Took {perf_counter()-t_reg:.5f}s to register '+\
                              f'{len(varname_list)} variables.')
            if verbose: print(f'Took a total of {perf_counter()-t0:.5f}s to kamodofy '+\
                              f'{len(gvar_list)} variables.')
        
        #define and register a 4D variable-----------------------------------------
        def register_variable(self, units, variable, varname, gridded_int):
            x_, y_, z_ = self.get_grid(varname) # variable may have different grid positions in the staggered grid of the model
            xvec_dependencies = {'time':'hr','x':'R_E','y':'R_E','z':'R_E'}

            self = self.custom_interp(units, variable, self._time, x_, y_, z_ , 
                                      varname, xvec_dependencies, gridded_int)   #regdef_4D_interpolators     
            return
        
        def get_grid(self, varname):
            """fetch the grid positon for this variable"""

            if varname == 'B1_x':
                return self._x_bx, self._y_bx, self._z_bx
            elif varname == 'B1_y':
                return self._x_by, self._y_by, self._z_by
            elif varname == 'B1_z':
                return self._x_bz, self._y_bz, self._z_bz
            elif varname == 'E1_x':
                return self._x_ex, self._y_ex, self._z_ex
            elif varname == 'E1_y':
                return self._x_ey, self._y_ey, self._z_ey
            elif varname == 'E1_z':
                return self._x_ez, self._y_ez, self._z_ez
            else: # (default) positions on plasma grid
                return self._x, self._y, self._z     
            
        def custom_interp(self, units, variable, t, x, y, z, varname, 
                          xvec_dependencies, gridded_int):
            '''define interpolator same as normal one EXCEPT need inner boundary limit.
            Decrease memory demand by not saving data arrays.'''

            rgi = RegularGridInterpolator((t, x, y, z), 
                                          variable, bounds_error = False, fill_value=NaN)
            
            @kamodofy(units=units, data={})
            def interpolator(xvec):  #xvec = [[t1,x1,y1,z1],[t2,x2,y2,z2],...]
                """Interpolates 4d variable without a grid"""
                
                xvec = array(xvec)
                r = sqrt(sum(xvec[:,1:]**2, axis=1))  #calculate radius
                return where(r>self.near_Earth_boundary_radius, rgi(xvec), NaN)
                #return interpolator if true, NaN if false
            
            self = register_interpolator(self, varname, interpolator, 
                                                     xvec_dependencies)
            
            #define and register the gridded interpolator if desired
            if gridded_int:
                self.variables[varname+'_ijk'] = dict(units = units, data = {}) 
                gridded_interpolator = define_4d_gridded_interpolator(units,{},t,x,
                                                                      y,z,xvec_dependencies,
                                                                      interpolator)
                self = register_interpolator(self, varname+'_ijk', 
                                                         gridded_interpolator, 
                                                         xvec_dependencies)
            return            
                    
        
    return MODEL
        
        
