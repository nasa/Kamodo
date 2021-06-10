# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 18:17:56 2021
@author: rringuet
Convert data in CTIPe output files to wrapped versions
"""
import numpy as np
import time as ti
from netCDF4 import Dataset


def ctipe_wrap_variables(variable, variable_name):
    '''Copied from Lutz. wrap variables as needed for ctipe model ouput'''
    
    if 'electric_field' in variable_name:
        # CTIPe efield variables from neutral file do not need to be wrapped but
        # need to be transposed from (time,lon,lat) to (time,lat,lon)
        new_variable = np.transpose(variable,[0,2,1])
    elif len(variable.shape) == 3: # 3D variable
        new_variable = np.concatenate((variable, variable[:,:,0:1]), axis = 2)
    elif len(variable.shape) == 4: # 4D variable
   	    new_variable = np.concatenate((variable, variable[:,:,:,0:1]), axis = 3)
    elif (variable_name == 'lon') and (variable.max() < 360.):
        new_variable = np.append(variable,360.)
    else: 
        new_variable = variable
    return new_variable.__array__()
     
def ctipe_wrap_files(filename):
    filetype_list = ['plot-density','plot-height','plot-neutral','plot-plasma']
    file_beg, file_end = [filename.split(filetype) for filetype in filetype_list 
                          if len(filename.split(filetype))>1][0]
    filename_density, filename_height, filename_neutral = [
        file_beg+filetype+file_end for filetype in filetype_list[:-1]]
    
    #looping through datasets
    for filename in [filename_density, filename_height, filename_neutral]:
        data = Dataset(filename)  #get data from file
        
        #start new wrapped output object
        data_out = Dataset(filename[:-3]+'-wrapped.nc', 'w', format='NETCDF4')
        for dim in data.dimensions.keys():  
            if dim == 'lon':
                new_dim = data_out.createDimension(dim, data.dimensions[dim].size+1)
            else: 
                new_dim = data_out.createDimension(dim, data.dimensions[dim].size)
            
        #copy over variables to file, wrapping for longitude 
        for variable_name in data.variables.keys():  #wrap longitude and dependent variables
            var = data.variables[variable_name]
            name, datatype, dimensions = var.name, var.datatype, var.dimensions
            if 'electric_field' in variable_name: #transpose dimensions
                dimensions = (var.dimensions[0], var.dimensions[2], var.dimensions[1])
            new_var = data_out.createVariable(name, datatype, dimensions)
            var_data = ctipe_wrap_variables(var.__array__(), variable_name)
            new_var[:] = var_data
            
        #close file
        print(f'{filename} converted.')
        data_out.close()  
    return filename_density[:-3]+'-wrapped.nc'
    

if __name__=='__main__':
    tic=ti.perf_counter()                 
    #define file names (input and output)
    file_dir = 'C:/Users/rringuet/Kamodo_WinDev1/CTIPe/'
    filename = file_dir+'2015-03-18-plot-density.nc'
    new_filename = ctipe_wrap_files(filename)
    print(f'{ti.perf_counter()-tic:.6f}s')