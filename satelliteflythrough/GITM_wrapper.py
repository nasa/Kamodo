# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 15:21:50 2021

@author: rringuet
"""
#coding imports
import glob
import numpy as np
import time as ti
from kamodo.readers.gitm_4Dcdf import GITM, gitm_varnames
import kamodo.satelliteflythrough.wrapper_utilities as U
 
#need input times to be timestamps since 1970 to pull data from the correct files
   
def GITMVariables():
    '''return a list of all the possible variables in GITM'''
    
    return gitm_varnames

def FlyAway(file_prefix, variable_list, sat_time, sat_height, sat_lat, sat_lon, 
                  plot_sampling=4000, plot_file='', verbose=False):
    '''fly satellite through GITM model data, per day'''

    #Check that sat data is all the same length, will error if not
    U.sat_data_check(sat_time, sat_height, sat_lat, sat_lon)
    
    #create gitm kamodo object
    gitm = GITM(file_prefix, variables_requested=variable_list, printfiles=True, 
                gridded_int=False)
    
    #create satellite tracks and interpolate data for each variable
    results = U.Generic_FlyAway(gitm, variable_list, sat_time, sat_height, sat_lat,
                    sat_lon, plot_file, plot_sampling, verbose=False)
        
    if verbose: print(f'Done for {file_prefix}\n')
    return results

def GITM_SatelliteFlythrough(file_dir, variable_list, sat_time, sat_height, sat_lat, 
                            sat_lon, dt=450., plots=False, daily_plots=False, plot_close=True, 
                            plot_sampling=4000, verbose=False):
    '''
    Execute flythrough for GITM model data. Returns results_dict, results_units.
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
    file_patterns = np.unique([file_dir+'Data/*'+f.split('/')[-1].split('\\')[-1][5:13] for f in files])
    print(file_patterns)
    varlist_3d = ['TEC', 'NmF2', 'hmF2','SolarLocalTime','SolarZenithAngle',
              'phi_qJoule','phi_q','phi_qEUV','phi_qNOCooling']  #add 3D items to list as discovered ----------
    results_dict, results_units = U.generic_SatelliteFlythrough(file_dir, file_patterns, 
                                sat_time, GITM, FlyAway, variable_list, sat_height, 
                                sat_lat, sat_lon, daily_plots, plot_sampling, 
                                gitm_varnames, varlist_3d, plots, plot_close, 
                                dt=dt, verbose=verbose)
    return results_dict, results_units

if __name__=='__main__':
    ''' Begin program '''
    #initialize input parameters (filename, variable_name)
    file_dir = 'C:/Users/rringuet/Kamodo_WinDev1/GITM/'
    variable_list = ['T_n','phi_qJoule']#,'TEC']  #Test ilev with N_n, without ilev with T_e, 3D with TEC
    t0 = ti.perf_counter()
    
    #generate a fake satellite track    
    from kamodo.satelliteflythrough.SatelliteFlythrough import SampleTrajectory as ST
    traj_dict = ST(1165968000.0-30000., 1166053801.0+1000., n=1)
    print(f'{len(traj_dict["sat_time"])} satellite locations.')
    
    results_dict, results_units = GITM_SatelliteFlythrough(file_dir, variable_list, 
                                    traj_dict['sat_time'], traj_dict['sat_height'], 
                                    traj_dict['sat_lat'], traj_dict['sat_lon'], dt=450.,
                                    plots=True, daily_plots=True, plot_close=True,
                                    verbose=False)
    print(f'{ti.perf_counter()-t0:.5f}s')

    from kamodo.satelliteflythrough import FlythroughPlots as FPlot
    for var in variable_list:
            FPlot.Plot4Dcar(var,results_dict['sat_time'],results_dict['sat_lat'],
                results_dict['sat_lon'],results_dict['sat_height'], results_dict[var], 
                results_units[var], 'C:/Users/rringuet/Kamodo_WinDev1/GITM/Plots/', plot_close=False, 
                plot_sampling=10000)        