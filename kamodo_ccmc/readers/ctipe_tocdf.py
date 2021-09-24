# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 18:17:56 2021
@author: rringuet
Convert data in CTIPe output files to wrapped versions
"""
#import numpy as np
from numpy import transpose, zeros, array, append
from time import perf_counter
from netCDF4 import Dataset
from astropy.constants import R_earth

file_varnames = {'density':['rho',0,['time','lon_d','lat_d','lev'],'kg/m**3'],
                  'temperature':['T',1,['time','lon_d','lat_d','lev'],'K'],
                  'electron_temperature':['T_e',2,['time','lon_h','lat_h','radius'],'K'],
                  'ion_temperature':['T_i',3,['time','lon_h','lat_h','radius'],'K'],
                  'height_d':['H_lev',4,['time','lon_d','lat_d','lev'],'m'],   
                  'height_n':['H_ilev',4,['time','lon_n','lat_n','ilev'],'m'],                    
                  'meridional_neutral_wind':['Vn_lat',5,['time','lon_n','lat_n','ilev'],'m/s'],
                  'zonal_neutral_wind':['Vn_lon',6,['time','lon_n','lat_n','ilev'],'m/s'],
                  'vertical_neutral_wind':['Vn_H',7,['time','lon_n','lat_n','ilev'],'m/s'],
                  'neutral_temperature':['T_n',8,['time','lon_n','lat_n','ilev'],'K'],
                  'mean_molecular_mass':['Rmt',9,['time','lon_d','lat_d','lev'],'amu'],
                  'electron_density':['N_e',10,['time','lon_h','lat_h','radius'],'1/m**3'],
                  'neutral_density':['N_n',11,['time','lon_n','lat_n','ilev'],'1/m**3'],
                  'solar_heating':['Q_Solar',12,['time','lon_n','lat_n','ilev'],'J/kg/s'],
                  'joule_heating':['Q_Joule',13,['time','lon_n','lat_n','ilev'],'J/kg/s'],
                  'radiation_heat_cool':['Q_radiation',14,['time','lon_n','lat_n','ilev'],'J/kg/s'],
                  'atomic_oxygen_density':['N_O',15,['time','lon_n','lat_n','ilev'],'1/m**3'],
                  'molecular_oxygen_density':['N_O2',16,['time','lon_n','lat_n','ilev'],'1/m**3'],
                  'molecular_nitrogen_density':['N_N2',17,['time','lon_n','lat_n','ilev'],'1/m**3'],
                  'nitric_oxide_density':['N_NO',18,['time','lon_n','lat_n','ilev'],'1/m**3'],
                  'nitric_oxide_ion_density':['N_NOplus',19,['time','lon_n','lat_n','ilev'],'1/m**3'],
                  'molecular_nitrogen_ion_density':['N_N2plus',20,['time','lon_n','lat_n','ilev'],'1/m**3'],  
                  'molecular_oxygen_ion_density':['N_O2plus',21,['time','lon_n','lat_n','ilev'],'1/m**3'],
                  'atomic_nitrogen_ion_density':['N_Nplus',22,['time','lon_n','lat_n','ilev'],'1/m**3'],
                  'atomic_oxygen_ion_density':['N_Oplus',23,['time','lon_h','lat_h','radius'],'1/m**3'],
                  'atomic_hydrogen_ion_density':['N_Hplus',24,['time','lon_h','lat_h','radius'],'1/m**3'],
                  'pedersen_conductivity':['Sigma_P',25,['time','lon_n','lat_n','ilev'],'S/m'],
                  'hall_conductivity':['Sigma_H',26,['time','lon_n','lat_n','ilev'],'S/m'],
                  'zonal_ion_velocity':['Vi_lon',27,['time','lon_n','lat_n','ilev'],'m/s'],
                  'meridional_ion_velocity':['Vi_lat',28,['time','lon_n','lat_n','ilev'],'m/s'],
                  #start 3D variables
                  'height_integrated_joule_heating':['W_Joule',29,['time','lon_n','lat_n'],'W/m**2'],
                  'energy_influx':['Eflux_precip',30,['time','lon_n','lat_n'],'W/m**2'],
                  'mean_energy':['Eavg_precip',31,['time','lon_n','lat_n'],'keV'],
                  'total_electron_content':['TEC',32,['time','lon_n','lat_n'],'1/m**2'], #'10**16/m**2'
                  'theta_electric_field_at_140km':['E_theta140km',33,['time','Elon','Elat'],'V/m'],
                  'lambda_electric_field_at_140km':['E_lambda140km',34,['time','Elon','Elat'],'V/m'],
                  'theta_electric_field_at_300km':['E_theta300km',35,['time','Elon','Elat'],'V/m'],
                  'lambda_electric_field_at_300km':['E_lambda300km',36,['time','Elon','Elat'],'V/m']}


def ctipe_wrap_variables(var_dict, variable_name):
    '''wrap variables in longitude and transpose as needed for ctipe model output'''
    
    if 'electric_field' in variable_name:
        # CTIPe efield variables from neutral file do not need to be wrapped but
        # need to be transposed from (time,lon,lat) to (time,lat,lon)
        #new_variable = np.transpose(variable,[0,2,1])
        pass  #want in (time, lon, lat)
    elif len(var_dict['data'].shape) == 3: # 3D variable, wrap in longitude then transpose
        shape_list = list(var_dict['data'].shape)  # time, lat, lon 
        shape_list[2]+=1  #need one more place in longitude
        tmp_arr = zeros(shape_list)  #array to set-up wrapped data in
        tmp_arr[:,:,:-1]=var_dict['data']  #copy data into grid
        tmp_arr[:,:,-1] = var_dict['data'][:,:,0]  #wrap in longitude 
        var_dict['data'] = transpose(tmp_arr, (0,2,1))  #(t,lat,lon) -> (t,lon,lat)
        var_dict['size'] = (shape_list[0],shape_list[2],shape_list[1])
    elif len(var_dict['data'].shape) == 4: # 4D variable
        shape_list = list(var_dict['data'].shape)  # time, lat, lon, height 
        shape_list[3]+=1  #need one more place in longitude
        tmp_arr = zeros(shape_list)  #array to set-up wrapped data in
        tmp_arr[:,:,:,:-1]=var_dict['data']  #copy data into grid
        tmp_arr[:,:,:,-1] = var_dict['data'][:,:,:,0]  #wrap in longitude 
        var_dict['data'] = transpose(tmp_arr, (0,3,2,1))   #(t,h,lat,lon) -> (t,lon,lat,h)
        var_dict['size'] = (shape_list[0],shape_list[3],shape_list[2],shape_list[1])
    return var_dict
 
'''
    elif (variable_name == 'lon') and (variable.max() < 360.):
        new_variable = np.append(variable,360.)
'''    
def ctipe_combine_files(file_prefix, verbose=False):
    '''Combine data from 3 files, wrapping in longitude and transposing as necessary.'''
    
    tic=perf_counter() 
    #determine file names of group
    filetype_list = ['-plot-density.nc','-plot-height.nc','-plot-neutral.nc']
    filename_density, filename_height, filename_neutral = [
        file_prefix+filetype for filetype in filetype_list]
    
    #open data files
    ctipe_density = Dataset(filename_density)
    ctipe_height = Dataset(filename_height)  #in meters
    ctipe_neutral = Dataset(filename_neutral)
    
    #retrieve data and key properties from each file
    d_dict={key:{'data':array(ctipe_density.variables[key]),
                      'datatype':ctipe_density.variables[key].datatype,
                      'size':ctipe_density.variables[key].size} \
              for key in ctipe_density.variables.keys()}
    h_dict={key:{'data':array(ctipe_height.variables[key]),
                  'datatype':ctipe_height.variables[key].datatype,
                  'size':ctipe_height.variables[key].size} \
          for key in ctipe_height.variables.keys()}
    n_dict={key:{'data':array(ctipe_neutral.variables[key]),
                      'datatype':ctipe_neutral.variables[key].datatype,
                      'size':ctipe_neutral.variables[key].size} \
              for key in ctipe_neutral.variables.keys()}
        
    #close files
    ctipe_density.close()
    ctipe_height.close()
    ctipe_neutral.close()
    
    #wrap longitude dimensions
    d_dict['lon']['data'] = append(d_dict['lon']['data'], 360.)    
    d_dict['lon']['size']+=1
    h_dict['lon']['data'] = append(h_dict['lon']['data'], 360.)    
    h_dict['lon']['size']+=1        
    n_dict['lon']['data'] = append(n_dict['lon']['data'], 360.)    
    n_dict['lon']['size']+=1      

    #collect dimensions data, assuming time and ilev are all the same
    #assuming lon and lat in diff files might be different
    dim_dict={}
    dim_dict['time'] = d_dict['time']
    dim_dict['lat_d'] = d_dict['lat']
    dim_dict['lat_h'] = h_dict['lat']
    dim_dict['lat_n'] = n_dict['lat']
    dim_dict['lon_d'] = d_dict['lon']
    dim_dict['lon_h'] = h_dict['lon']
    dim_dict['lon_n'] = n_dict['lon']
    dim_dict['lev'] = d_dict['plev']
    dim_dict['ilev'] = n_dict['plev']
    dim_dict['Elat'] = n_dict['elat']
    dim_dict['Elon'] = n_dict['elon']

    #convert height in km to radius in R_E to align with SPH coord sys (since long is 0 to 360)
    height_dict = h_dict['ht']
    height_dict['data'] = (height_dict['data']+R_earth.value/1000.)/(R_earth.value/1000.)
    dim_dict['radius'] = height_dict
    
    #remove keys from file dictionaries for coordinates
    for key in ['time','lat','lon','elat','elon','ht','plev']:
        if key in d_dict.keys(): del d_dict[key]
        if key in h_dict.keys(): del h_dict[key]
        if key in n_dict.keys(): del n_dict[key]
        
    #adjust 'height' variable names to better distinguish
    d_dict['height_d'] = d_dict['height']
    del d_dict['height']
    n_dict['height_n'] = n_dict['height']
    del n_dict['height']
    
    #initialize single output file
    data_out = Dataset(file_prefix+'.nc', 'w', format='NETCDF4')
    data_out.model = 'CTIPe'
    data_out.file = ''.join([f+',' for f in [filename_density, filename_height, 
                                             filename_neutral]]).strip(',')  #csv list of files

    # store dimensions
    for dim in dim_dict.keys():
        if verbose: print(dim)
        new_dim = data_out.createDimension(dim, dim_dict[dim]['size'])
        new_var = data_out.createVariable(dim, dim_dict[dim]['datatype'],tuple((dim,)))
        new_var[:] = dim_dict[dim]['data']
    if verbose: print('Dimensions complete.\n')

    #add variable data from density file to output file
    for key in d_dict.keys():
        if key in ['ZMAG','mean_molecular_mass'] or key not in file_varnames.keys(): 
            continue  #using Rmt in neutral file
        if verbose: print(key,file_varnames[key][2])
        d_dict[key] = ctipe_wrap_variables(d_dict[key], key) 
        new_var = data_out.createVariable(file_varnames[key][0], d_dict[key]['datatype'],
                                          tuple(file_varnames[key][2]))
        new_var[:] = d_dict[key]['data'] 
    if verbose: print('Density file complete.\n')

    #add variable data from height file to output file
    for key in h_dict.keys():
        if key == 'ZMAG' or key not in file_varnames.keys(): 
            continue  #no other variables depend on ZMAG, so ignore
        if verbose: print(key,file_varnames[key][2])
        h_dict[key] = ctipe_wrap_variables(h_dict[key], key) 
        new_var = data_out.createVariable(file_varnames[key][0], h_dict[key]['datatype'],
                                          tuple(file_varnames[key][2]))
        new_var[:] = h_dict[key]['data'] 
    if verbose: print('Height file complete.\n')
        
    #add variable data from neutral file to output file
    for key in n_dict.keys():
        if key in ['ZMAG','electron_density','atomic_oxygen_ion_density',
                   'atomic_hydrogen_ion_density'] or key not in file_varnames.keys(): 
            continue  #N_e is in height file
        if verbose: print(key,file_varnames[key][2])
        n_dict[key] = ctipe_wrap_variables(n_dict[key], key) 
        new_var = data_out.createVariable(file_varnames[key][0], n_dict[key]['datatype'],
                                          tuple(file_varnames[key][2]))
        new_var[:] = n_dict[key]['data'] 
    if verbose: print('Neutral file complete.\n')
        
    #close file
    print(f"Data for {file_prefix} converted in {perf_counter()-tic:.6f}s.")
    data_out.close()  
    return file_prefix+'.nc'
    

if __name__=='__main__':
    #define file names (input and output)
    file_dir = 'C:/Users/rringuet/Kamodo_WinDev1/CTIPe/Data/'
    file_prefix = file_dir+'2015-03-18'
    new_filename = ctipe_combine_files(file_prefix)