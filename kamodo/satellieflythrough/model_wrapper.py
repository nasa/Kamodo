# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 11:55:58 2021

@author: rringuet

Instead of writing individual wrappers for each model, this file is used to
do anything model-specific, such as importing the right readers.
"""
import glob
import numpy as np


def FileSearch(model, file_dir):
    '''Returns list of model data files for each model based on the name pattern.
    If only one file per day, or reader knows of the different filenames, then return string.
    Else, return an array of filename patterns.'''
    
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
        return glob.glob(file_dir+'Data/*.nc')
    
    else:
        raise AttributeError('Model not yet added.')
        

def Model_Reader(model):
    '''Returns model reader for requested model.'''
    
    if model=='CTIPe':
        from kamodo.readers.ctipe_4D import CTIPe
        return CTIPe
    
    elif model=='IRI':
        from kamodo.readers.iri_4D import IRI
        return IRI
    
    elif model=='GITM':
        from kamodo.readers.gitm_4Dcdf import GITM
        return GITM
    
    elif model=='SWMF_IE':
        from kamodo.readers.swmfie_4Dcdf import SWMF_IE
        return SWMF_IE
    
    elif model=='TIEGCM':
        from kamodo.readers.tiegcm_4D import TIEGCM
        return TIEGCM
    
    else:
        raise AttributeError('Model not yet added.')
        

def Model_Variables(model, return_dict=False):
    '''Returns model variables for requested model.'''
    
    if model == '':  #Give a list of possible values
        print("Possible models are: 'CTIPe','IRI', 'GITM', 'SWMF_IE', and 'TIEGCM'")
        return
    
    #choose the model-specific function to retrieve the variables
    if model=='CTIPe':
        from kamodo.readers.ctipe_4D import ctipe_varnames as variable_dict
    
    elif model=='IRI':
        from kamodo.readers.iri_4D import iri_varnames as variable_dict
    
    elif model=='GITM':
        from kamodo.readers.gitm_4Dcdf import gitm_varnames as variable_dict
    
    elif model=='SWMF_IE':
        from kamodo.readers.swmfie_4Dcdf import swmfie_varnames as variable_dict
    
    elif model=='TIEGCM':
        from kamodo.readers.tiegcm_4D import tiegcm_varnames as variable_dict
        
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
    '''Return list of model variables that are three-dimensional.'''
    
    if model=='CTIPe':
        return ['W_Joule', 'Eflux_precip', 'Eavg_precip', 'TEC', 'E_theta140km',
       'E_lambda140km', 'E_theta300km', 'E_lambda300km']
    
    elif model=='IRI':
        return ['TEC', 'NmF2', 'HmF2']
    
    elif model=='GITM':
        return ['TEC', 'NmF2', 'hmF2','SolarLocalTime','SolarZenithAngle',
              'phi_qJoule','phi_q','phi_qEUV','phi_qNOCooling']

    elif model=='SWMF_IE':
        Var = Model_Variables(model, return_dict=True)
        return [value[0] for key, value in Var.items() if value[0] not in \
                ['x','y','z','theta','psi','theta_Btilt', 'psi_Btilt']]

    elif model=='TIEGCM':
        return ['T_nLBC','u_nLBC','v_nLBC','T_nLBCNM','u_nLBCNM','v_nLBCNM','TEC']