# -*- coding: utf-8 -*-
"""
Created on Mon May 17 18:53:35 2021

@author: rringuet
"""
from glob import glob
import numpy as np
from time import perf_counter
from datetime import datetime, timezone
from os.path import basename
#from astropy.constants import R_earth
from netCDF4 import Dataset

swmfie_varnames={"X":['x','km'],"Y":['y','km'],"Z":['z','km'],    #(ignored, given in R_E on a unit sphere)
                 "Theta":['theta',"deg"],"Psi":['psi',"deg"],     #(used as coordinates)
                 "Btilt_theta":['theta_Btilt',"deg"], "Btilt_psi":['psi_Btilt',"deg"], 
                 #(added directly to object for documentation purposes)
                 "SigmaH":['Sigma_H',"S"],"SigmaP":['Sigma_P',"S"],
                 "E-Flux":['Phi_E',"W/m**2"], "Ave-E":['E_avg','eV'],
                 "JR":["j_R","mA/m**2"],"PHI":["Phi","kV"],
                 "Ex":["E_x","mV/m"],"Ey":["E_y","mV/m"],"Ez":["E_z","mV/m"],
                 "Jx":["j_x","mA/m**2"],"Jy":["j_y","mA/m**2"],"Jz":["j_z","mA/m**2"],
                 "Ux":['v_x',"km/s"],"Uy":['v_y',"km/s"],"Uz":['v_z',"km/s"],
                 "JouleHeat":['Q_Joule',"mW/m**2"], "IonNumFlux":['Phi_nion',"1/cm**2/s"],
                 "RT 1/B":['Binv_RT',"1/T"],"RT Rho":['rho_RT',"amu/cm**3"],"RT P":['P_RT',"Pa"],
                 "conjugate dLat":['dLat_star',"deg"],"conjugate dLon":['dlon_star',"deg"]}
                 
 
'''                
#documentation variables:
"X":['x','km'],"Y":['y','km'],"Z":['z','km'],    (ignored, given in R_E on a unit sphere)
"Theta":['theta',"deg"],"Psi":['psi',"deg"],     (used as coordinates)
"Btilt_theta":['theta_Btilt',"deg"], "Btilt_psi":['psi_Btilt',"deg"] 
(added directly to object for documentation purposes)
'''
@np.vectorize
def dts_to_hrs(datetime_string, filedate):
    '''Get hours since midnight from datetime string'''
    
    return (datetime.strptime(datetime_string, '%Y-%m-%d %H:%M:%S').replace(tzinfo=timezone.utc)\
            -filedate).total_seconds()/3600.

@np.vectorize
def filename_to_dts(filename, string_date):
    '''Get datetime string in format "YYYY-MM-SS HH:mm:SS" from filename'''
    
    mmhhss = basename(filename)[12:18]
    return string_date+' '+mmhhss[:2]+':'+mmhhss[2:4]+':'+mmhhss[4:] 

def dts_to_ts(file_dts):
    '''Get datetime timestamp in UTC from datetime string'''    
    
    return datetime.timestamp(datetime.strptime(file_dts, '%Y-%m-%d %H:%M:%S'
                                                ).replace(tzinfo=timezone.utc))

def read_SWMFIE_header(filename):
    '''parse header for a representative file'''
    
    #get file data
    file_object = open(filename, 'r')
    file_contents = file_object.read()
    file_object.close()
    
    #sort info with TITLE keyword
    title_lines = file_contents.split('TITLE=')[1].split('VARIABLES=')[0].\
        replace('\n',',').replace('"','').strip().strip(',')
    title, time, BtiltDeg_str = title_lines.split(',')
    #get time into 'YYYY-MM-DD HH:mm:SS' format, clean up BtiltDeg
    dts = time[:10]+' '+''.join([time[11:].split('-')[i]+':' \
                                 for i in [0,1]])+time[11:].split('-')[2]
    BtiltDeg = BtiltDeg_str.split('=')[1].split()

    #sort info with VARIABLE keyword
    variable_lines = file_contents.split('VARIABLES=')[1].split('ZONE')[0].\
        replace('\n',',').replace('"','').strip().strip(',')
    variable_list = [var.strip().split(']')[0].split(' [')[0] for var in \
                 variable_lines.split(',')]  #units in dict above
    
    #sort info with ZONE keyword 'IonN N=0000241 T=0000:20:00'
    zone_lines = file_contents.split('ZONE T="')[1].split('"')[0].strip()
    run_type1, N1, t = zone_lines.split(' ')
    N1 = int(N1.split('=')[1].strip())
    zone_lines = file_contents.split('ZONE T="')[2].split('"')[0].strip()
    run_type2, N2, t = zone_lines.split(' ')
    N2 = int(N2.split('=')[1].strip())  
    
    #sort coordinate info   'I=           91  J=          181  F=POINT'
    coord_info = file_contents.split(zone_lines)[1].split('\n')[1].strip()
    text, i, text, j, f = coord_info.split()
    i, j, f = int(i), int(j), f.split('=')[1]
    #if sometimes different, need another set for this
    
    #determine number of lines to skip for data sections
    line_test = [f in c for c in file_contents.split('\n')]
    skip1, skip2 = np.where(np.array(line_test)==True)[0]
    
    #determine how many lines constitute one data list
    data_lines = file_contents.split('\n')[skip1+1:]
    num_test = [len(line) for line in data_lines]  #list of lengths of each line
    ndata_lines = np.diff(np.where(np.array(num_test)==num_test[0])[0])[0]
    
    #cleanup and return
    del file_contents
    header = {'title':title, 'time':dts, 'BtiltDeg':BtiltDeg, 
              'variable_list':variable_list, 'run_type1':run_type1, 'N1':N1,
              'run_type2':run_type2, 'N2':N2, 'i':i, 'j':j, 'f':f, 'skip1':skip1,
              'skip2':skip2, 'ndata_lines':ndata_lines}
    return header

def read_SWMFIE_data(filename, header):
    '''read only data from file and return in labeled dictionary'''
    
    #get data from file using skip as number of lines to skip
    file_object = open(filename, 'r')
    file_contents = file_object.readlines()
    file_object.close()

    #Get theta_Btilt value from file
    title_line = file_contents[0].replace('\n',',').replace('"','').strip().strip(',')
    theta_Btilt = float(title_line.split(',')[-1].strip().split()[1])

    #collect data into an array, one for each run_type
    file_data1 = file_contents[header['skip1']+1:header['skip2']-1]
    data1 = np.array([''.join(file_data1[i:i+header['ndata_lines']]).strip('\n').split() \
                      for i in range(0, len(file_data1), header['ndata_lines'])], dtype=float)
    file_data2 = file_contents[header['skip2']+1:]
    data2 = np.array([''.join(file_data2[i:i+header['ndata_lines']]).strip('\n').split() \
                      for i in range(0, len(file_data2), header['ndata_lines'])], dtype=float)
        
    #determine correct dimension sizes from theta and psi (3 and 4), combine data and return in dict
    nLatA, nLon = np.unique(data1[:,3]).size, np.unique(data1[:,4]).size    
    data = np.concatenate((np.reshape(data1,(nLon,nLatA,len(header['variable_list']))),
                           np.reshape(data2,(nLon,nLatA,len(header['variable_list'])))),
                          axis=1)[:,::-1,:]  #leave since grids are identical in N/S hemispheres
    
    return {swmfie_varnames[header['variable_list'][i]][0]:np.transpose(data[:,:,i],[1,0]) for i in \
                     range(len(header['variable_list']))}, theta_Btilt      
    

#file_prefix = 'C:/Users/rringuet/Kamodo_WinDev1/SWMF_IE/i_e20180801'  #example
def _read_SWMFIE(file_prefix, verbose=False): 
    '''file_prefix must be of form "3D***_tYYMMDD" to load all files for one day
     and include a complete path to the files'''
    t0 = perf_counter()
    
    #establish time attributes first for file searching
    files = glob(file_prefix+'*.tec')   #take one day of data
    file_datestr = basename(file_prefix)[3:11]
    string_date = file_datestr[:4]+'-'+file_datestr[4:6]+'-'+file_datestr[6:8]  #'YYYY-MM-DD' 
    print('CONVERTER:', file_prefix, file_datestr, string_date)
    filedate = datetime.strptime(string_date+' 00:00:00', \
                                      '%Y-%m-%d %H:%M:%S').replace(tzinfo=timezone.utc) #dt object

    #get list of variables possible in these files using first file
    #t_header = perf_counter()
    if len(files)>1:
        file_header = read_SWMFIE_header(files[1])  #first file has a different header format
    else:
        file_header = read_SWMFIE_header(files[0])
    
    #if verbose: print(f'Took {perf_counter()-t_header:.6f}s to parse header.')
    file_varlist = file_header['variable_list']
    gvar_list = [swmfie_varnames[key][0] for key in file_varlist if key not in \
                     ['X','Y','Z','Theta','Psi','Btilt_theta','Btilt_psi']]  
        #avoid returning coordinates stored elsewhere (or ignored)
    #print(gvar_list)
    
    #set time dimension from list of files
    time = dts_to_hrs(filename_to_dts(files, string_date), filedate)
    
    #collect coordinates and intialize data storage from first file
    data, Btilt = read_SWMFIE_data(files[0], file_header)
    lon, lat0 = np.unique(data['psi']), np.unique(data['theta'])-90.
    lon -= 180 #shifting longitude to be in range -180 to 180 for SM coordinates
    lat0 = np.insert(lat0, np.where(lat0==0.)[0], 0.)
    #interpolator requires unique ascending values, offset two zero values by 0.0001
    lat0[np.where(lat0==0.)[0]] = -0.0001, 0.0001
    lat = np.insert(lat0,0,-90.)
    lat = np.append(lat,90.)  #add values for poles for wrapping
    #height = np.array([R_earth.value/1000.])
    radius = np.array([1.])  #placeholder value of one earth radius
    #data is 2D (lon, lat) with time, so no height
    theta_Btilt = [Btilt]   #to append more values later
    psi_Btilt = np.repeat(0., len(files))   #this value is always zero

    # Store variable data and units for first file.
    var_units = {value[0]:value[-1] for key, value in swmfie_varnames.items() \
                 if value[0] in gvar_list}
    variables = {key:[[data[key]]] for key in gvar_list}

    #determine where to split arrays in longitude
    lon_idx = min(np.where(lon>0)[0])
    #print(lon.size, lon_idx, lon.size%2)    

    #add data for other files
    if verbose: print(f'Reading {len(files)} files...')
    for f in files[1:]: 
        data, Btilt = read_SWMFIE_data(f, file_header)
        theta_Btilt.append(Btilt)
        for var_key in variables: 
            variables[var_key].append([data[var_key]])
            
    #convert to arrays
    theta_Btilt = np.array(theta_Btilt, dtype=float)
    for var_key in variables:
        variables[var_key] = np.concatenate(tuple(variables[var_key]))
        new_data = np.transpose(variables[var_key], (0,2,1))   #(t,lat,lon) -> (t,lon,lat)
        #print(var_key, new_data.shape)
        
        #original longitude range is from 0 to 360, moving data to be from -180 to 180 instead
        #without changing the longitude reference point
        tmp = np.zeros(new_data.shape)
        tmp[:,:lon_idx-lon.size%2,:] = new_data[:,lon_idx:,:]
        tmp[:,lon_idx-lon.size%2:,:] = new_data[:,:lon_idx,:]
        
        #wrap latitude for scalar variables (all to be treated as scalars)
        new_shape = list(tmp.shape)
        new_shape[2]+=2  #one value on each end
        tmp2 = np.zeros(new_shape)
        tmp2[:,:,1:-1] = tmp
        # put in top values
        top = np.mean(tmp2[:,:,1],axis=1)  #same shape as time axis
        new_top = np.broadcast_to(top, (new_shape[1],new_shape[0])).T
        tmp2[:,:,0] = new_top
        #same for bottom, reusing variable names
        top = np.mean(tmp2[:,:,-2],axis=1)  #same shape as time axis
        new_top = np.broadcast_to(top, (new_shape[1],new_shape[0])).T
        tmp2[:,:,-1] = new_top        
        variables[var_key] = tmp2  #store result
        
    coords = {'time':time, 'radius':radius, 'lat':lat, 'lon':lon}
    variables['theta_Btilt'], variables['psi_Btilt'] = theta_Btilt, psi_Btilt
    var_units['theta_Btilt'], var_units['psi_Btilt'] = 'deg','deg'
    if verbose: print(f'Took {perf_counter()-t0:.6f}s to read and '+\
                      f'assemble data from {len(files)} files.')
    return files, coords, variables, var_units, filedate
        
def _toCDF(filename, files, coords, variables, var_units, filedate):
        
    #start new wrapped output object
    cdf_filename = filename+'.nc'
    data_out = Dataset(cdf_filename, 'w', format='NETCDF4')
    data_out.file = ''.join([f+',' for f in files]).strip(',')
    data_out.model = 'SWMF_IE'
    data_out.filedate = filedate.strftime('%Y-%m-%d %H:%M:%S')
    for dim in coords.keys():  #check that the datatype works correctly when opening the file******
        new_dim = data_out.createDimension(dim, coords[dim].size)  #create dimension
        new_var = data_out.createVariable(dim, np.float64, dim)  #create variable
        new_var[:] = coords[dim]  #store data for dimension in variable
        if dim=='radius': units = 'R_E'
        elif dim=='time': units = 'hr'
        else: units='deg'
        new_var.units = units
        
    #copy over variables to file
    var_1D = ['theta_Btilt','psi_Btilt']
    for variable_name in variables.keys(): 
        if len(variables[variable_name].shape)==3:
            new_var = data_out.createVariable(variable_name, np.float64, ('time','lon','lat'))
            new_data = variables[variable_name]
        elif variable_name in var_1D:
            new_var = data_out.createVariable(variable_name, np.float64, ('time'))
            new_data = variables[variable_name]
        new_var[:] = new_data  #store data in variable
        new_var.units = var_units[variable_name]        
        
    #close file
    data_out.close()     
    
    return cdf_filename

def convertSWMFIE_toCDF(file_prefix):
    '''Convert files of given pattern into one netCDF4 file'''
    
    print(f'NetCDF version of data not found. Converting files with {file_prefix} prefix to netCDF.')
    ftic = perf_counter()
    files, coords, variables, var_units, filedate = _read_SWMFIE(file_prefix, verbose=True)
    if coords==0: return False #if only one file, perform a clean return
    cdf_filename = _toCDF(file_prefix, files, coords, variables, var_units, filedate)
    print(f'{len(files)} files with prefix {file_prefix} now combined into {cdf_filename} '+\
          f'in {perf_counter()-ftic:.6f}s.')
    return True