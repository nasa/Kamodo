# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 11:55:58 2021

@author: rringuet

Instead of writing individual wrappers for each model, this file is used to
do anything model-specific, such as importing the right readers.


Help on coordinate systems:
    The resource information on SpacePy's coordinate conversion function is 
        sparse at best, so the below information has been collected via other
        resources and our own testing of the function. The data concerning the
        spherical coordinate systems are collected into a table format for 
        easier perusal.
    For cartesian coordinates, all of the input values should be in earth radii (R_E)
        in order (x, y, z) to work properly.
    For spherical coordinates, all of the input values should be in order 
        (longitude, latitude, altitude or radius). The longitude and latitude
        values should be in degrees, altitude values in kilometers, and radius
        values in earth radii (R_E) from the Earth's center. All latitude values 
        should fall between -90 and 90 degrees. The longitude range differs 
        between the coordinate systems and is given for each in the table below.
        
SpacePy 
Abbrev.   Full Name                       Lon. range     vertical variable
--------------------------------------------------------------------------
GDZ    Geodetic (WGS 84)                  (-180, 180)    Altitude (km)
GEO    Geographic                         (-180, 180)    Radius (R_E)
GSM    Geocentric Solar Magnetospheric    (-180, 180)    Radius (R_E)
GSE    Geocentric Solar Ecliptic          (-180, 180)    Radius (R_E)
SM     Solar Magnetic                     (-180, 180)    Radius (R_E)
GEI    Geocentric Equatorial Inertial     (-180, 180)    Radius (R_E)
      (also ECI = Earth-Centered Inertial)
MAG    Geomagnetic                        (-180, 180)    Radius (R_E)
SPH    Spherical                            (0, 360)     Radius (R_E)
RLL    Radius, Latitude, Longitude        (-180, 180)    Radius (R_E)

For descriptions of most of the coordinate systems, see 
https://sscweb.gsfc.nasa.gov/users_guide/Appendix_C.shtml and it's reference,
"Geophysical Coordinate Transformations", C.T. Russell, Cosmic Electrodynamics, Vol. 2, pp. 184 - 196, 1971.


"""
from glob import glob
from numpy import unique
from os.path import basename
import numpy as np


model_dict = {0:'CTIPe', 1:'GITM', 2:'IRI', 3:'SWMF_IE', 4:'TIEGCM', 5:'OpenGGCM_GM'}

def convert_model_string(model_int):
    '''Converts numerical model reference to string.'''
    
    if isinstance(model_int,int): 
        return model_dict[model_int]
    else:
        return model_int


def Choose_Model(model):
    '''Returns module specific to the model requested.'''
    #UPDATE THIS AS MORE MODELS ARE ADDED

    if model == '':  #Give a list of possible values
        print(f"Possible models are: {model_dict}")
        print('Integers or strings allowed.')
        return
    
    model = convert_model_string(model)  #convert to string
    
    if model=='CTIPe':
        import kamodo_ccmc.readers.ctipe_4D as module
        return module
    
    elif model=='IRI':
        import kamodo_ccmc.readers.iri_4D as module
        return module
    
    elif model=='GITM':
        import kamodo_ccmc.readers.gitm_4Dcdf as module
        return module
    
    elif model=='SWMF_IE':
        import kamodo_ccmc.readers.swmfie_4Dcdf as module
        return module
    
    elif model=='TIEGCM':
        import kamodo_ccmc.readers.tiegcm_4D as module
        return module

    elif model=='OpenGGCM_GM':
        import kamodo_ccmc.readers.openggcm_gm_4Dcdf_xarray as module
        return module
    
    else:
        raise AttributeError('Model not yet added.')    


def FileSearch(model, file_dir, call_type='normal'):
    '''Returns list of model data files for each model based on the name pattern.
    If only one file per day, or reader knows of the different filenames, then return string.
    Else, return an array of filename patterns.'''
    #UPDATE THIS AS NEW MODELS ARE ADDED
    
    if isinstance(model,int): model = model_dict[model]  #convert to string
    
    if model=='CTIPe':
        files = glob(file_dir+'*.nc')  #look for wrapped and original data
        file_patterns = unique([file_dir+basename(f)[:10] for f in files \
                                            if 'CTIPe' not in basename(f)])
        return file_patterns  
    
    elif model=='IRI':
        return file_dir+'IRI.3D.*.nc'
    
    elif model=='GITM':  #whole day version of filesearch
        files = glob(file_dir+'*')  #next line returns list of prefixes: e.g. 3DALL_t20150315
        if call_type=='normal': #give prefix for full day files
            file_patterns = unique([file_dir+'*'+basename(f)[7:13] for f in files\
                                    if 'GITM' not in basename(f) and '.nc' not in basename(f)])
        else:  #give prefix for hourly files
            file_patterns = unique([file_dir+'*'+basename(f)[7:16] for f in files \
                                    if '.nc' not in basename(f) and 'GITM' not in basename(f)])
        return file_patterns     
    
    elif model=='SWMF_IE':
        files = glob(file_dir+'i_e*')  #next line returns list of prefixes: e.g. i_e20150315
        if call_type=='normal': #give prefix for full day files
            file_patterns = unique([file_dir+basename(f)[:11] for f in files\
                                            if '.nc' not in basename(f)])
        else:  #give prefix for hourly files
            file_patterns = unique([file_dir+basename(f)[:14] for f in files\
                                            if '.nc' not in basename(f)])
        return file_patterns        

    elif model=='TIEGCM':
        print('Please remove all pxxx.nc files if present.')
        return file_dir+'*.nc'
    
    elif model=='OpenGGCM_GM':
        files = sorted(glob(file_dir+'*.nc'))
        file_patterns = unique([file_dir+basename(f).split('3df_')[0]+'3df_'+f.split('3df_')[1][:13]\
                                            for f in files])  #hourly files only  
        return file_patterns
    
    else:
        raise AttributeError('Model not yet added.')
        

def Model_Reader(model):
    '''Returns model reader for requested model. Model agnostic.'''
    
    module = Choose_Model(model)
    return module.MODEL()  #imports Kamodo
        

def Model_Variables(model, return_dict=False):
    '''Returns model variables for requested model. Model agnostic.'''
    
    #choose the model-specific function to retrieve the variables
    module = Choose_Model(model)
    variable_dict = module.model_varnames
    var_dict = {value[0]:value[1:] for key, value in variable_dict.items()}
        
    #retrieve and print model specific and standardized variable names
    if return_dict: 
        return var_dict
    else:
        print('\nThe model accepts the standardized variable names listed below.')
        #print('Units for the chosen variables are printed during the satellite flythrough if available.')
        print('-----------------------------------------------------------------------------------')
        for key, value in sorted(var_dict.items()): print(f"{key} : '{value}'")
        print()
        return    
    
def File_Variables(model, file_dir, return_dict=False):
    '''Print list of variables in model data output stored in file_dir.'''

    reader = Model_Reader(model)
    file_patterns = FileSearch(model, file_dir)
    file_variables={}
    #collect file variables in a nested dictionary
    if isinstance(file_patterns, list) or isinstance(file_patterns,np.ndarray):
        #print(model, file_patterns)
        for file_pattern in file_patterns:
            kamodo_object = reader(file_pattern, fulltime=False, variables_requested='all')
            file_variables[file_pattern] = {key:value for key, value \
                                            in sorted(kamodo_object.var_dict.items())}

    else:  #reader requires full filenames, not a file pattern
        files = glob(file_patterns)  #find full filenames for given pattern
        #print(model, files)
        for file in files:
            kamodo_object = reader(file, fulltime=False, variables_requested='all')
            file_variables[file] = {key:value for key, value \
                                            in sorted(kamodo_object.var_dict.items())}
        
    #either return or print nested_dictionary
    if return_dict: 
        return file_variables
    else:
        #print file pattern standardized variable names
        for file_key in file_variables.keys():
            print(f'\nThe file {file_key} contains the following standardized variable names:')
            #print('Units for the chosen variables are printed during the satellite flythrough if available.')
            print('-----------------------------------------------------------------------------------')
            for key, value in file_variables[file_key].items(): print(f"{key} : '{value}'")
            print()
        return
    
def File_Times(model, file_dir):
    '''Return/print time ranges available in the data in the given dir. Also
    performs file conversions of new data if needed.'''
    
    #get time ranges from data
    from kamodo_ccmc.flythrough.SF_utilities import check_timescsv
    file_patterns = FileSearch(model, file_dir)
    times_dict = check_timescsv(file_patterns, model)
    
    #print time ranges for given file groups: file_pattern, beg time, end time, etc
    print('File pattern: UTC time ranges')
    print('------------------------------------------')
    for key in times_dict.keys():
        print(f"{times_dict[key][0]} : {times_dict[key][1:]}")
    
    return times_dict

def Var_3D(model):
    '''Return list of model variables that are three-dimensional. Model agnostic.'''
    
    #choose the model-specific function to retrieve the 3D variable list
    variable_dict = Model_Variables(model, return_dict=True)    
    return [value[0] for key, value in variable_dict.items() if len(value[4])==3]

def Var_ilev(model):
    '''Return list of possible ilev coordinate names for model given.'''

    variable_dict = Model_Variables(model, return_dict=True)   
    ilev_list = list(unique([value[4][-1] for key, value in variable_dict.items() \
                                if len(value[4])==4 and 'ilev' in value[4][-1]]))
    return ilev_list

def Var_units(model, variable_list):
    '''Returns dictionary of key, value = varname, units.'''

    variable_dict = Model_Variables(model, return_dict=True)  
    return {key:value[-1] for key, value in variable_dict.items() if key in variable_list}    

def convert_variablenames(model, variable_list):
    '''Given list of integers for variable names, convert to names for given model.'''

    if isinstance(variable_list[0], int):
        variable_dict = Model_Variables(model, return_dict=True)          
        tmp_var = [value[0] for key, value in variable_dict.items()\
                               if value[2] in variable_list]
        variable_list = tmp_var
    return variable_list
    
def convert_coordnames(coord_type, coord_grid):
    '''convert integers to strings for coordinate names and grid types.'''
    
    if isinstance(coord_type, int) and isinstance(coord_grid, int):
        spacepy_dict = {0:'GDZ',1:'GEO',2:'GSM',3:'GSE',4:'SM',5:'GEI',6:'MAG',
                        7:'SPH',8:'RLL'}
        coord_type = spacepy_dict[coord_type]
        if coord_grid==0:
            coord_grid = 'car'
        elif coord_grid==1:
            coord_grid = 'sph'
    return coord_type, coord_grid
    
def coord_units(coord_type, coord_grid):
    '''return proper units given coordinate system'''
    
    if coord_grid=='car':
        if coord_type in ['GDZ','SPH','RLL']:
            print(f'There is no cartesian version in the {coord_type} coordinate system.')
            return
        return {'utc_time':'s','net_idx':'','c1':'R_E','c2':'R_E','c3':'R_E'}
    elif coord_grid=='sph':
        if coord_type=='GDZ':
            return {'utc_time':'s','net_idx':'','c1':'deg','c2':'deg','c3':'km'}
        else:
            return {'utc_time':'s','net_idx':'','c1':'deg','c2':'deg','c3':'R_E'}
        
def coord_names(coord_type, coord_grid):
    '''given type and grid of coordinates, return dicitonary of coordinate names.'''
    
    if coord_grid=='car':
        return {'c1':'X_'+coord_type, 'c2':'Y_'+coord_type, 'c3':'Z_'+coord_type}
    elif coord_grid=='sph' and coord_type=='GDZ':
        return {'c1':'Longitude', 'c2':'Latitude', 'c3':'Altitude'}
    elif coord_grid=='sph' and coord_type!='GDZ':
        return {'c1':'Longitude', 'c2':'Latitude', 'c3':'Radius'}
