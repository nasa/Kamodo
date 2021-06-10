# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 15:21:50 2021

@author: rringuet
"""
#coding imports
import glob
import numpy as np
from kamodo.readers.swmfie_4Dcdf import SWMF_IE, swmfie_varnames
import kamodo.satelliteflythrough.wrapper_utilities as U
 
#need input times to be timestamps since 1970 to pull data from the correct files
   
def SWMFIEVariables():
    '''return a list of all the possible variables in SWMFIE'''
    
    return swmfie_varnames

def FlyAway(file_prefix, variable_list, sat_time, sat_height, sat_lat, sat_lon, 
                  plot_sampling=4000, plot_file='', verbose=False):
    '''fly satellite through SWMFIE model data, per day'''

    #Check that sat data is all the same length, will error if not
    U.sat_data_check(sat_time, sat_height, sat_lat, sat_lon)
    
    #create SWMFIE kamodo object
    swmfie = SWMF_IE(file_prefix, variables_requested=variable_list, printfiles=True, 
                gridded_int=False)
    
    #create satellite tracks and interpolate data for each variable
    results = U.Generic_FlyAway(swmfie, variable_list, sat_time, sat_height, sat_lat,
                    sat_lon, plot_file, plot_sampling, verbose=False)
        
    if verbose: print(f'Done for {file_prefix}\n')
    return results

def SWMFIE_SatelliteFlythrough(file_dir, variable_list, sat_time, sat_height, sat_lat, 
                            sat_lon, dt=450., plots=False, daily_plots=False, plot_close=True, 
                            plot_sampling=4000, verbose=False):
    '''
    Execute flythrough for SWMFIE model data. Returns results_dict, results_units.
    results_dict is a dictionary of the interpolated data for the entire data set
        sorted by variable name.
    results_units is a dictionary of the units from the model for each variable.
    file_dir is a string indicating where the data files are located.
    variable_list is a list of strings of the desired variables. 
    sat_time is an array of timestamp values.
    sat_height is an array of heights in meters.
    sat_lat and sat_lon ar arrays of latitudes and longitudes in degrees.
    Set plots=True to get plots of the entire data set for each variable.
    Set daily_plots=True to get plots of the data set for each variable and each file.
    Set plot_close=False to keep plots of the entire data set open.
    '''
    
    files = glob.glob(file_dir+'Data/*')  #next line returns list of prefixes: e.g. 3DALL_t20150315
    file_patterns = np.unique([file_dir+'Data/'+f.split('/')[-1].split('\\')[-1][:11] for f in files])
    varlist_3d = [value[0] for key, value in swmfie_varnames.items() \
                         if value[0] not in ['x','y','z','theta','psi','theta_Btilt', 'psi_Btilt']]   #add 3D items to list as discovered ----------
    results_dict, results_units = U.generic_SatelliteFlythrough(file_dir, file_patterns, 
                                sat_time, SWMF_IE, FlyAway, variable_list, sat_height, 
                                sat_lat, sat_lon, daily_plots, plot_sampling, 
                                swmfie_varnames, varlist_3d, plots, plot_close, 
                                dt=dt, verbose=verbose)
    return results_dict, results_units

if __name__=='__main__':
    ''' Begin program '''
    #initialize input parameters (filename, variable_name)
    file_dir = 'C:/Users/rringuet/Kamodo_WinDev1/SWMF_IE/'
    variable_list = ['Q_Joule']  
    
    #generate a fake satellite track    
    from kamodo.satelliteflythrough.SatelliteFlythrough import SampleTrajectory as ST
    traj_dict = ST(1533081600.0-30000., 1533167760.0+1000., n=1)
    print(f'{len(traj_dict["sat_time"])} satellite locations.')
    
    results_dict, results_units = SWMFIE_SatelliteFlythrough(file_dir, variable_list, 
                                    traj_dict['sat_time'], traj_dict['sat_height'], 
                                    traj_dict['sat_lat'], traj_dict['sat_lon'], dt=450.,
                                    plots=True, daily_plots=True, plot_close=True,
                                    verbose=False)
