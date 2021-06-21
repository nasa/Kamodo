# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 11:55:58 2021

@author: rringuet

Instead of writing individual wrappers for each model, this file is used to
do anything model-specific, such as importing the right readers.
"""
import glob
import numpy as np


def Choose_Model(model):
    '''Returns module specific to the model requested.'''
    #UPDATE THIS AS MORE MODELS ARE ADDED

    if model == '':  #Give a list of possible values
        print("Possible models are: 'CTIPe','IRI', 'GITM', 'SWMF_IE', and 'TIEGCM'")
        return
    
    if model=='CTIPe':
        import kamodo.readers.ctipe_4D as module
        return module
    
    elif model=='IRI':
        import kamodo.readers.iri_4D as module
        return module
    
    elif model=='GITM':
        import kamodo.readers.gitm_4Dcdf as module
        return module
    
    elif model=='SWMF_IE':
        import kamodo.readers.swmfie_4Dcdf as module
        return module
    
    elif model=='TIEGCM':
        import kamodo.readers.tiegcm_4D as module
        return module
    
    else:
        raise AttributeError('Model not yet added.')    


def FileSearch(model, file_dir):
    '''Returns list of model data files for each model based on the name pattern.
    If only one file per day, or reader knows of the different filenames, then return string.
    Else, return an array of filename patterns.'''
    #UPDATE THIS AS NEW MODELS ARE ADDED
    
    if model=='CTIPe':
        files = glob.glob(file_dir+'Data/*-plot-density*.nc')  #look for wrapped and original data
        file_patterns = np.unique([file_dir+'Data/'+f.split('/')[-1].split('\\')[-1][:23]+\
                                   '-wrapped.nc' for f in files]) 
        return file_patterns  
    
    elif model=='IRI':
        return file_dir+'Data/IRI.3D.*.nc'
    
    elif model=='GITM':
        files = glob.glob(file_dir+'Data/*')  #next line returns list of prefixes: e.g. 3DALL_t20150315
        file_patterns = np.unique([file_dir+'Data/*'+f.split('/')[-1].split('\\')[-1][5:13] for f in files])
        return file_patterns     
    
    elif model=='SWMF_IE':
        files = glob.glob(file_dir+'Data/*')  #next line returns list of prefixes: e.g. 3DALL_t20150315
        file_patterns = np.unique([file_dir+'Data/'+f.split('/')[-1].split('\\')[-1][:11] for f in files])
        return file_patterns        

    elif model=='TIEGCM':
        return file_dir+'Data/*.nc'
    
    else:
        raise AttributeError('Model not yet added.')
        

def Model_Reader(model):
    '''Returns model reader for requested model. Model agnostic.'''
    
    module = Choose_Model(model)
    return module.MODEL
        

def Model_Variables(model, return_dict=False):
    '''Returns model variables for requested model. Model agnostic.'''
    
    #choose the model-specific function to retrieve the variables
    module = Choose_Model(model)
    variable_dict = module.model_varnames
        
    #retrieve and print model specific and standardized variable names
    if return_dict: 
        return variable_dict
    else:
        print('\nThe model accepts the standardized variable names listed below.')
        print('Units for the chosen variables are printed during the satellite flythrough if available.')
        print('-----------------------------------------------------------------------------------')
        for key, value in variable_dict.items(): print(f"{key} : '{value}'")
        print()
        return    
    
    
#saving list of 3D variables for now. Need to move into variable dicts instead.
def Var_3D(model):
    '''Return list of model variables that are three-dimensional. Model agnostic.'''
    
    #choose the model-specific function to retrieve the 3D variable list
    module = Choose_Model(model)
    variable_dict = module.model_varnames    
    return [value[0] for key, value in variable_dict.items() if value[-2]=='3D']
