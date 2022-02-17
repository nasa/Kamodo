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
    #'bx1':['B1_x','x component of magnetic field (on grid cell faces)',3,'GSE','car',['time','x','x','x'],'nT'],
    #'by1':['B1_y','y component of magnetic field (on grid cell faces)',4,'GSE','car',['time','y','y','y'],'nT'],
    #'bz1':['B1_z','z component of magnetic field (on grid cell faces)',5,'GSE','car',['time','z','z','z'],'nT'],
    'ex':['E_x','x component of electric field',6,'GSE','car',['time','x','x','x'],'mV/m'],
    'ey':['E_y','y component of electric field',7,'GSE','car',['time','y','y','y'],'mV/m'],
    'ez':['E_z','z component of electric field',8,'GSE','car',['time','z','z','z'],'mV/m'],
    'vx':['v_plasmax','x component of plasma velocity',9,'GSE','car',['time','x','y','z'],'km/s'],
    'vy':['v_plasmay','y component of plasma velocity',10,'GSE','car',['time','x','y','z'],'km/s'],
    'vz':['v_plasmaz','z component of plasma velocity',11,'GSE','car',['time','x','y','z'],'km/s'],
    'rr':['N_plasma','number density of plasma (hydrogen equivalent)',12,'GSE','car',['time','x','y','z'],'1/cm**3'],
    'resis':['eta','resistivity',13,'GSE','car',['time','x','y','z'],'m**2/s'],
    'pp':['P_plasma','plasma pressure',14,'GSE','car',['time','x','y','z'],'pPa'],
    'xjx':['j_x','current density, x component',15,'GSE','car',['time','x','y','z'],'muA/m**2'],
    'xjy':['j_y','current density, y component',16,'GSE','car',['time','x','y','z'],'muA/m**2'],
    'xjz':['j_z','current density, z component',17,'GSE','car',['time','x','y','z'],'muA/m**2'],
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

#sample file name: 'D:/OpenGGCM_GM/Data/Yihua_1/Yihua_Zheng_090721_1.3df_2015-10-16_12.nc'
#main function/class definition
def MODEL():
    from numpy import array, NaN, diff, sqrt, sum, expand_dims, zeros, where
    from time import perf_counter
    import xarray, dask
    from os.path import isfile
    from kamodo import Kamodo, kamodofy
    from kamodo_ccmc.readers.reader_utilities import register_interpolator, define_4d_gridded_interpolator
    
    class MODEL(Kamodo):
        '''OpenGGCM_GM magnetosphere reader'''
        def __init__(self,full_file_prefix, variables_requested=[], runname = "noname",
                     filetime=False, verbose=False, gridded_int=True, printfiles=False, 
                     fulltime=True, missing_value=NaN, **kwargs):
            super(MODEL, self).__init__()
            t0=perf_counter() # profiling time stamp
            
            #convert files to netcdf4 if needed
            nc_file = full_file_prefix+'.nc'  # input file name: file_dir/YYYY-MM-DD_HH.nc
            if isfile(nc_file):  #file already prepared!
                self.conversion_test = True  #default value
            else:  #file not prepared, prepare it
                try:  #I don't have the file converter, so leave in try/except for now
                    from kamodo_ccmc.readers.openggcm_to_cdf import openggcm_combine_magnetosphere_files as gmconv
                    self.conversion_test = gmconv(full_file_prefix)
                    #should return a boolean (True is successful, False if not)
                except:
                    self.conversion_test = False

            #data are time-wrapped in files. This logic prevents the flythrough from breaking from this decision.
            #cdf_data.added_time_at_beginning and cdf_data.added_time_at_end = 0 if not, 1 if yes
            cdf_data = xarray.open_dataset(nc_file, chunks={'time':100,'x':100,'y':100,'z':100})
            if not fulltime and cdf_data.added_time_at_end:  #use unwrapped time values
                t = cdf_data.variables['_time'].values[:-1]  #skip added time at end
            else: #need wrapped time array for interpolation
                t = cdf_data.variables['_time'].values
            self._time = t

            #establish time attributes first
            self.filedate = datetime.strptime(cdf_data.filedate+' 00:00:00', 
                                              '%Y-%m-%d %H:%M:%S').replace(tzinfo=timezone.utc)
            if len(t)>1: self.dt = diff(t).max()*3600.  #t is in hours since midnight
            else: self.dt = 0
            self.datetimes=[datetime.utcfromtimestamp(hrs_to_ts(t[0], self.filedate)).isoformat(sep=' '),
                            datetime.utcfromtimestamp(hrs_to_ts(t[-1], self.filedate)).isoformat(sep=' ')]
            self.filetimes=[hrs_to_ts(t[0], self.filedate), hrs_to_ts(t[-1], self.filedate)]   #timestamps for matching in wrapper

            #return time information only for flythrough
            if filetime: 
                return  
            
            #if variables are given as integers, convert to standard names
            if len(variables_requested)>0:
                if isinstance(variables_requested[0], int):
                    print('Integers detected. Converting...', end="")
                    tmp_var = [value[0] for key, value in model_varnames.items()\
                                           if value[2] in variables_requested]
                    variables_requested = tmp_var
                    print('Converted:', variables_requested)

            #perform initial check on variables_requested list
            if len(variables_requested)>0 and fulltime and variables_requested!='all':
                test_list = [value[0] for key, value in model_varnames.items()]
                err_list = [item for item in variables_requested if item not in test_list]
                if len(err_list)>0: print('Variable name(s) not recognized:', err_list)
                
            #collect variable list            
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
            else:  #only input variables on the avoid_list if specifically requested
                avoid_list = []   #empty for now
                gvar_list = [key for key in cdf_data.variables.keys() \
                             if key in model_varnames.keys() and \
                                 key not in avoid_list]    
                if not fulltime and variables_requested=='all':
                    self.var_dict = {value[0]: value[1:] for key, value in model_varnames.items() \
                            if key in gvar_list}
                    return                     

            # Store variable's units and reference to Dataset object in memory
            variables = {model_varnames[key][0]:{'units':model_varnames[key][-1],
                                   'data':getattr(cdf_data, key)}\
                              for key in gvar_list} 
            if verbose: print('Done reading in variable data.', full_file_prefix)
    
            #store variables
            self.near_Earth_boundary_radius = cdf_data.near_Earth_boundary_radius
            self.near_Earth_boundary_radius_unit = cdf_data.near_Earth_boundary_radius_units
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
            for grid in grid_list:
                setattr(self, grid, getattr(cdf_data, grid).values)  #store coordinate data
            cdf_data.close()  #done with file
    
            #register interpolators for each variable
            varname_list, self.variables = [key for key in variables.keys()], {}  #store original list b/c gridded interpolators
            t_reg = perf_counter()
            for varname in varname_list:  #all are 3D variables
                #make dimension names uniform for easier interpolation later
                dim_names, dims_dict = list(variables[varname]['data'].dims), {}
                dims_dict[dim_names[0]], dims_dict[dim_names[1]] = 'time', 'x'
                dims_dict[dim_names[2]], dims_dict[dim_names[3]] = 'y', 'z'
                variables[varname]['data'] = variables[varname]['data'].rename(
                    **dims_dict)
                
                #store and register data                
                self.variables[varname] = dict(units = variables[varname]['units'], 
                                               data = variables[varname]['data'])     #not saving data to decrease memory demand      
                self.register_variable(self.variables[varname]['units'], 
                                       self.variables[varname]['data'], varname, 
                                       gridded_int)
            #cdf_data.close()
            if verbose: print(f'Took {perf_counter()-t_reg:.5f}s to register '+\
                              f'{len(varname_list)} variables.')
            if verbose: print(f'Took a total of {perf_counter()-t0:.5f}s to kamodofy '+\
                              f'{len(gvar_list)} variables.')
        
        #define and register a 4D variable-----------------------------------------
        def register_variable(self, units, variable, varname, gridded_int):
            x_, y_, z_ = self.get_grid(varname) # variable may have different grid positions in the staggered grid of the model
            xvec_dependencies = {'time':'hr','x':'R_E','y':'R_E','z':'R_E'}
            
            #variable is a DataArray object with a built-in interpolator, add coordinates
            variable = variable.assign_coords({'time':self._time,'x':x_,'y':y_,'z':z_})

            self = self.custom_interp(units, variable, varname, xvec_dependencies, 
                                      gridded_int)   #regdef_4D_interpolators     
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
            
        def custom_interp(self, units, variable, varname, 
                          xvec_dependencies, gridded_int):
            '''define interpolator based on xarray's except need inner boundary limit.
            Decrease memory demand by not saving data arrays.'''

            @kamodofy(units=units, data=variable)
            def interpolator(xvec):  #xvec = [[t1,x1,y1,z1],[t2,x2,y2,z2],...]
                """Interpolates 4d variable without a grid"""
                
                t0 = perf_counter()
                xvec_arr = array(xvec)
                if len(xvec_arr.shape)==1: xvec_arr = expand_dims(xvec_arr, 0)
                r = sqrt(sum(xvec_arr[:,1:]**2, axis=1))  #calculate radius 1D array
                result = zeros(len(r))  #option 2, not much time diff from option 1
                for i in range(len(r)):
                    if r[i]>self.near_Earth_boundary_radius:
                        result[i] = float(variable.interp(time=xvec_arr[i][0], x=xvec_arr[i][1], 
                                                          y=xvec_arr[i][2], z=xvec_arr[i][3]).values)
                    else: 
                        result[i] = NaN
                print(f'Took {perf_counter()-t0:.3f} s for {len(xvec)} positions.')
                return result
            
            self = register_interpolator(self, varname, interpolator, xvec_dependencies)
            
            #define and register the gridded interpolator if desired
            if gridded_int:
                self.variables[varname+'_ijk'] = dict(units = units, data = variable) 
                
                @kamodofy(units=units, data={}, arg_units=xvec_dependencies)
                def gridded_interpolator(time, x, y, z):  
                    """Interpolates 4d variable on a grid"""
                    
                    t0 = perf_counter()
                    time, x, y, z = array(time), array(x), array(y), array(z)
                    r = sqrt(x**2+y**2+z**2)
                    return where(r>self.near_Earth_boundary_radius, 
                                 variable.interp(time=time, x=x, y=y, z=z).values, NaN)

                self = register_interpolator(self, varname+'_ijk', 
                                                         gridded_interpolator, 
                                                         xvec_dependencies)
            return            
        
    return MODEL
