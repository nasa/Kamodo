# -*- coding: utf-8 -*-
"""
Convert files with uniform grid to netcdf4
@author: rringuet, 2022

Date:  2020-05-05 00:00
Model: TS18
Bin:   Esw    1.1 mV/m, Bang   178 deg., tilt  13.1 deg.
Grid:  Uniform (lat_step: 1.00, lon_step: 2.00 [deg])

MLAT [deg]   MLT [hr]   Pot [kV] Vazm [deg] Vmag [m/s]
---------- ---------- ---------- ---------- ----------

"""
from glob import glob
import numpy as np
from time import perf_counter
from datetime import datetime, timezone
#from astropy.constants import R_earth
from netCDF4 import Dataset
import re

model_varnames={"Pot":['V','kV'],"Vazm":['theta_v','deg'],"Vmag":['v','m/s'],
                 #remaining variables are time series
                 "tilt":['theta_Btilt',"deg"], 'Esw':['E_sw','mV/m'],
                 'Bang':['theta_B','deg'],
                 #these are the coordinate variables
                 'MLAT':['MLAT','deg'], 'MLT':['MLT','hr']}
 
def dts_to_hrs(datetime_string, filedate):
    '''Get hours since midnight from datetime string'''
    
    return (datetime.strptime(datetime_string, '%Y-%m-%d %H:%M').replace(tzinfo=timezone.utc)\
            -filedate).total_seconds()/3600.
        
def grid_type(filename):
    '''Determine grid type of data file. True if uniform, False if equal-area.'''
    
    read_obj = open(filename, 'r')
    line = read_obj.readline().strip()  #Date:  2020-05-05 00:00
    line = read_obj.readline().strip()  #Model: TS18
    line = read_obj.readline().strip()  #Bin:   Esw    1.1 mV/m, Bang   178 deg., tilt  13.1 deg.
    line = read_obj.readline().strip()  #Grid:  Uniform (lat_step: 1.00, lon_step: 2.00 [deg])  
    read_obj.close()
    return 'Uniform' in line

#filename='C:/Users/rringuet/Kamodo_WinDev1/SuperDARN/fullday/model20200505-0000.txt'
def ascii_reader(filename):
    '''Loads the data from a superdarn txt file into a nested dict'''
    
    #open file
    read_obj = open(filename, 'r')
    
    #extract header
    date_string = read_obj.readline().strip()  #Date:  2020-05-05 00:00
    model_string = read_obj.readline().strip()  #Model: TS18
    bin_string = read_obj.readline().strip()  #Bin:   Esw    1.1 mV/m, Bang   178 deg., tilt  13.1 deg.
    grid_string = read_obj.readline().strip()  #Grid:  Uniform (lat_step: 1.00, lon_step: 2.00 [deg])
    trash = read_obj.readline().strip()  #empty line
    variable_keys = read_obj.readline().strip() #MLAT [deg]   MLT [hr]   Pot [kV] Vazm [deg] Vmag [m/s]
    trash = read_obj.readline().strip()  #---------- ---------- ---------- ---------- ----------
    
    #extract info from header strings
    time_str = date_string[6:].strip()  #'2020-05-05 00:00'  date_string[2]+' '+date_string[3]
    filedate = datetime.strptime(time_str[:10], '%Y-%m-%d').replace(tzinfo=timezone.utc)
    hrs = dts_to_hrs(time_str, filedate)
    bin_list = bin_string[4:].split(',')
    Esw = float(bin_list[0].strip()[3:].strip().split(' ')[0])
    Bang = float(bin_list[1].strip()[4:].strip().split(' ')[0])
    tilt = float(bin_list[2].strip()[4:].strip().split(' ')[0])
    var_list = re.split(' +',variable_keys)
    header_keys = ['tilt','Esw','Bang']
    variable_keys = [item for item in var_list if '[' not in item]
    
    #create dictionary to store data in
    variables = {model_varnames[var][0]: {'units': model_varnames[var][-1], 
                                          'data': []} for var in variable_keys+header_keys}
    #store time series values
    variables[model_varnames['tilt'][0]]['data'] = tilt
    variables[model_varnames['Esw'][0]]['data'] = Esw
    variables[model_varnames['Bang'][0]]['data'] = Bang
    variables['time'] = {'units':'hr', 'data': hrs}
        
    #store array data into dictionary
    for line in read_obj:
        vals = re.split(' +', line.strip())
        for i in range(len(variable_keys)): #skip empty block(s) at the end
            variables[model_varnames[variable_keys[i]][0]]['data'].append(vals[i])
    
    #convert to numpy float arrays
    for key in variables.keys():
        if isinstance(variables[key]['data'],(list)):
            variables[key]['data'] = np.array(variables[key]['data'], dtype=float)
    
    #add metadata
    variables['metadata'] = {'grid': grid_string[0][5:].strip(), 
                             'model': model_string[0][6:].strip(),
                             'filedate': time_str[:10]}    
    return variables

def _toCDF(files, file_prefix):
    '''Reads in data from all files, writes to a netcdf4 file. Used for uniform grid
    data files for faster data access.'''
        
    #get data from first file, set lat/lon arrays
    file_data = ascii_reader(files[0])
    lat = np.unique(file_data['MLAT']['data'])
    if sum(lat>0)==lat.size: 
        cdf_filename = file_prefix+'_default_N.nc'  #northern hemisphere data
        lat = np.append(lat, 90.)
    else: 
        cdf_filename = file_prefix+'_default_S.nc'  #southern hemisphere data
        lat = np.insert(lat, 0, -90.)
    lon = np.unique(file_data['MLT']['data'])*15.
    
    #set up net variables dictionary and time coordinate lists
    time = [file_data['time']['data']]
    var1D_list = ['theta_Btilt', 'E_sw', 'theta_B']
    var3D_list = ['V', 'theta_v', 'v']
    variables = {var: [np.reshape(file_data[var]['data'], (lat.size-1, lon.size)).T] \
                 for var in var3D_list}  #reshape into lon/lat array
    for var in var1D_list: variables[var] = [file_data[var]['data']]
    
    #loop through files and add data to variables dict
    for file in files[1:]:
        data = ascii_reader(file)
        time.append(data['time']['data'])
        for var in var1D_list: variables[var].append(file_data[var]['data'])
        for var in var3D_list: variables[var].append(np.reshape(file_data[var]['data'], 
                                                                (lat.size-1, lon.size)).T)
    for var in var1D_list+var3D_list: variables[var] = np.array(variables[var])
    
    #perform longitude wrapping in coordinate grid
    lon_le180 = np.where(lon<=180)[0]  
    lon_ge180 = np.where(lon>=180)[0]  #repeat 180 for -180 values 
    if not 180. in lon:  #add a cushion value for proper interpolation range (-180 to 180)
        lon_le180 = np.append(lon_le180, lon_le180.max()+1)
        lon_ge180 = np.insert(lon_ge180, 0, lon_ge180.min()-1)
    lon_size = len(lon_le180)+len(lon_ge180)
    tmp = np.zeros(lon_size)
    tmp[:len(lon_ge180)] = lon[lon_ge180]-360.
    tmp[len(lon_ge180):] = lon[lon_le180]
    lon = tmp    
    
    #perform lon and lat wrapping in variable data
    for var in var3D_list:
        #perform scalar averaging for pole values (latitude wrapping)
        data_shape = variables[var].shape
        total_shape = (data_shape[0],data_shape[1],data_shape[2]+1)
        tmp = np.zeros(total_shape, dtype=float)        
        if '_N.nc' in cdf_filename:  #north pole at end of array
            tmp[:,:,:-1] = variables[var]  #copy data into grid
            top = np.mean(tmp[:,:,-2],axis=1)  #same shape as time axis
            tmp[:,:,-1] = np.broadcast_to(top, (total_shape[1],total_shape[0])).T            
        elif '_S.nc' in cdf_filename:  #south pole at beginning of array
            tmp[:,:,1:] = variables[var]  #copy data into grid
            top = np.mean(tmp[:,:,1],axis=1)  #same shape as time axis
            tmp[:,:,0] = np.broadcast_to(top, (total_shape[1],total_shape[0])).T
        variables[var] = tmp 
            
        #swap longitudes, repeat 180 values for -180 position
        data_shape = variables[var].shape
        total_shape = (data_shape[0],lon_size,data_shape[2])
        tmp = np.zeros(total_shape, dtype=float)
        tmp[:,:len(lon_ge180),:] = variables[var][:,lon_ge180,:]
        tmp[:,len(lon_ge180):,:] = variables[var][:,lon_le180,:]   
        variables[var] = tmp
    
    #Data wrangling complete. Start new output file
    data_out = Dataset(cdf_filename, 'w', format='NETCDF4')
    data_out.file = ''.join([f+',' for f in files]).strip(',')
    data_out.model = 'SuperDARN'
    data_out.filedate = file_data['metadata']['filedate']
    data_out.grid = file_data['metadata']['grid']
    data_out.internal_model = file_data['metadata']['model']
    
    #establish coordinates (lon, lat, then time open)
    #lon
    new_dim = data_out.createDimension('lon', lon.size)  #create dimension
    new_var = data_out.createVariable('lon', np.float32, 'lon')  #create variable
    new_var[:] = lon  #store data for dimension in variable
    new_var.units = 'deg'
    #lat
    new_dim = data_out.createDimension('lat', lat.size)  #create dimension
    new_var = data_out.createVariable('lat', np.float32, 'lat')  #create variable
    new_var[:] = lat  #store data for dimension in variable
    new_var.units = 'deg'    
    #time
    new_dim = data_out.createDimension('time', len(files))  #create dimension
    new_var = data_out.createVariable('time', np.float32, 'time')  #create variable
    new_var[:] = np.array(time)     
    new_var.units = 'hr'
    
    #copy over variables to file
    for variable_name in variables.keys(): 
        if variable_name in var3D_list:
            new_var = data_out.createVariable(variable_name, np.float32, ('time','lon','lat'))
            new_data = variables[variable_name]
        elif variable_name in var1D_list:
            new_var = data_out.createVariable(variable_name, np.float32, ('time'))
            new_data = variables[variable_name]
        else: 
            continue
        new_var[:] = new_data  #store data in variable
        units = [value[-1] for key, value in model_varnames.items() if value[0]==variable_name][0]
        new_var.units = units
        
    #close file
    data_out.close()     
    
    return cdf_filename

def _toCDFGroup(files, file_prefix):
    '''Reads in data from all files, writes to h5 files. Used for equal-area
    data files so that lon grids from different lat vals can be stored in groups.'''

    #get data from first file
    file_data = ascii_reader(files[0])
    lat = np.unique(file_data['MLAT']['data'])
    if sum(lat>0)==lat.size: 
        cdf_filename = file_prefix+'_equalarea_N.nc'  #northern hemisphere data
        lat = np.append(lat, 90.)
    else: 
        cdf_filename = file_prefix+'_equalarea_S.nc'  #southern hemisphere data
        lat = np.insert(lat, 0, -90.)
    
    #set up net variables dictionary and time coordinate lists
    time = [file_data['time']['data']]
    var1D_list = ['theta_Btilt', 'E_sw', 'theta_B']
    var3D_list = ['V', 'theta_v', 'v']
    variables_1D = {var: [file_data[var]['data']] for var in var1D_list}
    variables_3D = {var: {latval: [] for latval in lat} for var in var3D_list}
    
    #store longitude locations and grids for each latitude value
    lonval_dict, lonidx_dict = {}, {}
    for latval in lat:
        idx = np.where(file_data['MLAT']['data']==latval)[0]
        lonidx_dict[latval] = idx  #store indices
        if len(idx)==0: continue  #skip latval=90.
        lon = file_data['MLT']['data'][idx]*15.  #convert from MLT to degrees
        if lon.max()>360.: lon[np.argmax(lon)]-=360.  #grid at pole has issues
        lonval_dict[latval] = np.unique(lon)
        if len(lonval_dict[latval])!=len(lon):  #last lon value repeated near pole
            lonidx_dict[latval] = idx[:-1]  #remove repeated value  
        
    #Store variable data in nested fashion
    for file in files:
        data = ascii_reader(file)
        for var in var3D_list:
            for latval in lat:  #1D array of vals for given time and lat 
                variables_3D[var][latval].append(data[var]['data'][lonidx_dict[latval]])   

    #perform latitude wrapping first to avoid repeated values in average
    for var in var3D_list:
        #store which latval key is closest to pole
        if '_N.nc' in cdf_filename: 
            pole_lat = 90. 
            latval = lat[-2]  #north pole at end of array
        elif '_S.nc' in cdf_filename: 
            pole_lat = -90.
            latval = lat[1] #south pole at beginning of array
            
        #perform scalar averaging and store
        #variables_3D[var][latval] has shape (time,lon), average over lon values
        variables_3D[var][pole_lat] = np.mean(np.array(variables_3D[var][latval]),axis=1)  #same shape as time
                
    #perform longitude swapping and lon wrapping per lat value
    for latval in lat:  
        if latval!=pole_lat:  #if not at the poles
            #wrap longitude coordinate values and store
            lon = lonval_dict[latval]
            lon_le180 = np.where(lon<=180)[0]  
            lon_ge180 = np.where(lon>=180)[0]  #repeat 180 for -180 values 
            if not 180. in lon:  #add a cushion value for proper interpolation range (-180 to 180)
                lon_le180 = np.append(lon_le180, lon_le180.max()+1)
                lon_ge180 = np.insert(lon_ge180, 0, lon_ge180.min()-1)
            lon_size = len(lon_le180)+len(lon_ge180)
            tmp = np.zeros(lon_size)
            tmp[:len(lon_ge180)] = lon[lon_ge180]-360.
            tmp[len(lon_ge180):] = lon[lon_le180]
            lonval_dict[latval] = tmp           
            
            #swap longitude dimension of variables, each of shape (time,lon)
            for var in var3D_list:
                variables_3D[var][latval] = np.array(variables_3D[var][latval])
                #swap longitudes, repeat 180 values for -180 position
                data_shape = variables[var].shape
                tmp = np.zeros((data_shape[0],lon_size), dtype=float)
                tmp[:,:len(lon_ge180)] = variables[var][:,lon_ge180]
                tmp[:,len(lon_ge180):] = variables[var][:,lon_le180]   
                variables[var] = tmp            

    #Data wrangling complete. Start new output file
    data_out = Dataset(cdf_filename, 'w', format='NETCDF4')
    data_out.file = ''.join([f+',' for f in files]).strip(',')
    data_out.model = 'SuperDARN'
    data_out.filedate = file_data['metadata']['filedate']
    data_out.grid = file_data['metadata']['grid']
    data_out.internal_model = file_data['metadata']['model']
    
    #establish coordinates (lat, then time open)
    #lat
    new_dim = data_out.createDimension('lat', lat.size)  #create dimension
    new_var = data_out.createVariable('lat', np.float32, 'lat')  #create variable
    new_var[:] = lat  #store data for dimension in variable
    new_var.units = 'deg'    
    #time
    new_dim = data_out.createDimension('time', len(files))  #create dimension
    new_var = data_out.createVariable('time', np.float32, 'time')  #create variable
    new_var[:] = np.array(time)     
    new_var.units = 'hr'
    
    
    #CHANGE TO NESTED/GROUP FORMAT!!!!
    #copy over variables to file
    for variable_name in variables.keys(): 
        if variable_name in var3D_list:
            new_var = data_out.createVariable(variable_name, np.float32, ('time','lon','lat'))
            new_data = variables[variable_name]
        elif variable_name in var1D_list:
            new_var = data_out.createVariable(variable_name, np.float32, ('time'))
            new_data = variables[variable_name]
        else: 
            continue
        new_var[:] = new_data  #store data in variable
        units = [value[-1] for key, value in model_varnames.items() if value[0]==variable_name][0]
        new_var.units = units
        
    #close file
    data_out.close()   
        







    return cdf_filename

def convert_files(file_prefix):
    '''Convert files of given pattern into one netCDF4 or h5 file'''
    #convert N and S hemisphere files separately or combine?
    
    print(f'Converted data file not found. Converting files with {file_prefix} prefix.')
    ftic = perf_counter()
    files = sorted(glob(file_prefix+'*.txt'))
    if grid_type(files[0]):  #If grid is uniform, writte to netcdf4
        print('Uniform grid detected. Converting to a netcdf4 file.')
        out_file = _toCDF(files, file_prefix)
    else:
        print('Equal-area grid detected. Converting to a grouped netcdf4 file.')
        out_file = _toCDFGroup(files, file_prefix)
    print(out_file)
    
    print(f'{len(files)} files with prefix {file_prefix} now combined into {out_file} '+\
          f'in {perf_counter()-ftic:.6f}s.')
    return True