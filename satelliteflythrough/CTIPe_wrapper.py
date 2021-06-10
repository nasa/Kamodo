# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 15:21:50 2021

@author: rringuet
"""
#coding imports
import glob
import numpy as np
from kamodo.readers.ctipe_4D import CTIPe, ctipe_varnames
import kamodo.satelliteflythrough.wrapper_utilities as U
 
#need input times to be timestamps since 1970 to pull data from the correct files
   
def CTIPeVariables():
    '''return a list of all the possible variables in CTIPe'''
    
    return ctipe_varnames

def FlyAway(filename, variable_list, sat_time, sat_height, sat_lat, sat_lon, 
                  plot_sampling=4000, plot_file='', verbose=False):
    '''fly satellite through CTIPe model data, per file'''
    
    #Check that sat data is all the same length, will error if not
    U.sat_data_check(sat_time, sat_height, sat_lat, sat_lon)
    
    #create ctipe kamodo object, initialize some variables
    ctipe = CTIPe(filename, variables_requested=variable_list, printfiles=True, 
                  gridded_int=False)
    if 'H' in variable_list: variable_list.remove('H')  #H only needed if other functions require ilev
    
    #create satellite tracks and interpolate data for each variable
    results = U.Generic_FlyAway(ctipe, variable_list, sat_time, sat_height, sat_lat,
                                sat_lon, plot_file, plot_sampling, verbose=False)
        
    if verbose: print(f'Done for {filename}\n')
    return results

def CTIPe_SatelliteFlythrough(file_dir, variable_list, sat_time, sat_height, sat_lat, 
                              sat_lon, dt=450., plots=False, daily_plots=False, plot_close=True, 
                              plot_sampling=4000, verbose=False):
    '''
    Execute flythrough for CTIPe model data. Returns results_dict, results_units.
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
    
    files = glob.glob(file_dir+'Data/*-plot-density*.nc')  #look for wrapped and original data
    #reader prefers wrapped filename, even if it does not exist. will create if no wrapped data found.
    file_patterns = np.unique([file_dir+'Data/'+f.split('/')[-1].split('\\')[-1][:23]+\
                               '-wrapped.nc' for f in files]) 
    varlist_3d = ['W_Joule', 'Eflux_precip', 'Eavg_precip', 'TEC', 'E_theta140km',
       'E_lambda140km', 'E_theta300km', 'E_lambda300km'] #add 3D items to list as discovered ----------
    results_dict, results_units = U.generic_SatelliteFlythrough(file_dir, file_patterns, 
                                sat_time, CTIPe, FlyAway, variable_list, sat_height, 
                                sat_lat, sat_lon, daily_plots, plot_sampling, 
                                ctipe_varnames, varlist_3d, plots, plot_close, 
                                dt=dt, verbose=verbose)
    return results_dict, results_units

if __name__=='__main__':
    ''' Begin program '''
    #initialize input parameters (filename, variable_name)
    file_dir = 'C:/Users/rringuet/Kamodo_WinDev1/CTIPe/'
    variable_list = ['T_e', 'TEC']  #Test ilev with N_n, without ilev with T_e, 3D with TEC
    #variable_list = ['rho', 'T', 'T_e', 'T_i', 'H', 'Vn_lat', 'Vn_lon', 'Vn_H', 
    #                 'T_n', 'Rmt', 'N_e', 'N_n', 'Q_Solar', 'Q_Joule', 'Q_radiation', 
    #                 'N_O', 'N_O2', 'N_N2', 'N_NO', 'N_NOplus', 'N_N2plus', 'N_O2plus', 
    #                 'N_Nplus', 'N_Oplus', 'N_Hplus', 'sigma_P', 'sigma_H', 'Vi_lon', 
    #                 'Vi_lat', 'W_Joule', 'Eflux_precip', 'Eavg_precip', 'TEC', 
    #                 'E_theta140km', 'E_lambda140km', 'E_theta300km', 'E_lambda300km']
    
    #generate a fake satellite track    
    from kamodo.satelliteflythrough.SatelliteFlythrough import SampleTrajectory as ST
    traj_dict = ST(1426660000.0-30000., 1426791280.0)
    
    results_dict, results_units = CTIPe_SatelliteFlythrough(file_dir, variable_list, 
                                    traj_dict['sat_time'], traj_dict['sat_height'], 
                                    traj_dict['sat_lat'], traj_dict['sat_lon'], dt=450.,
                                    plots=True, daily_plots=False, plot_close=False,
                                    verbose=False)