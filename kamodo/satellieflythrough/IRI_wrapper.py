# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 15:21:50 2021

@author: rringuet
"""
#coding imports
from kamodo.readers.iri_4D import IRI, iri_varnames
import kamodo.satelliteflythrough.wrapper_utilities as U
 
#need input times to be timestamps since 1970 to pull data from the correct files
   
def IRIVariables():
    '''return a list of all the possible variables in IRI'''
    
    return iri_varnames

def FlyAway(filename, variable_list, sat_time, sat_height, sat_lat, sat_lon, 
                  plot_sampling=4000, plot_file='', verbose=False):
    '''fly satellite through IRI model data, per file'''
    
    #Check that sat data is all the same length
    #Check that sat data is all the same length, will error if not
    U.sat_data_check(sat_time, sat_height, sat_lat, sat_lon)
    
    #create iri kamodo object
    iri = IRI(filename, variables_requested=variable_list, printfiles=True, 
              gridded_int=False)  #speed up calculation by not creating the gridifed interpolator
    
    #create satellite tracks and interpolate data for each variable
    results = U.Generic_FlyAway(iri, variable_list, sat_time, sat_height, sat_lat,
                    sat_lon, plot_file, plot_sampling, verbose=False)
        
    if verbose: print(f'Done for {filename}\n')
    return results

def IRI_SatelliteFlythrough(file_dir, variable_list, sat_time, sat_height, sat_lat, 
                            sat_lon, dt=450., plots=False, daily_plots=False, plot_close=True, 
                            plot_sampling=4000, verbose=False):
    '''
    Execute flythrough for IRI model data. Returns results_dict, results_units.
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
    
    pattern = file_dir+'Data/IRI.3D.*.nc'
    print(pattern)
    varlist_3d = ['TEC', 'NmF2', 'HmF2']  #add 3D items to list as discovered -----------
    results_dict, results_units = U.generic_SatelliteFlythrough(file_dir, pattern, 
                                sat_time, IRI, FlyAway, variable_list, sat_height, 
                                sat_lat, sat_lon, daily_plots, plot_sampling, 
                                iri_varnames, varlist_3d, plots, plot_close, 
                                dt=dt, verbose=verbose)
    return results_dict, results_units

if __name__=='__main__':
    ''' Begin program '''
    #initialize input parameters (filename, variable_name)
    file_dir = 'C:/Users/rringuet/Kamodo_WinDev1/IRI/'
    variable_list = ['T_e','TEC']  #Test ilev with N_n, without ilev with T_e, 3D with TEC
    
    #generate a fake satellite track    
    from kamodo.satelliteflythrough.SatelliteFlythrough import SampleTrajectory as ST
    traj_dict = ST(1495945560.0-30000., 1496014100.0, n=1)
    print(f'{len(traj_dict["sat_time"])} satellite locations.')
    
    results_dict, results_units = IRI_SatelliteFlythrough(file_dir, variable_list, 
                                    traj_dict['sat_time'], traj_dict['sat_height'], 
                                    traj_dict['sat_lat'], traj_dict['sat_lon'], dt=450.,
                                    plots=True, daily_plots=True, plot_close=False,
                                    verbose=False)
