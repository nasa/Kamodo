# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 14:57:55 2021
author: rringuet

Code to be called from other languages for a single trajectory point.
All desired model data should be in a single directory.


#Syntax of possible command line calls 
-------------------------------------
(1) python ./Kamodo-master/kamodo/flythrough/SingleSatelliteFlythrough.py
(prints possible models)

(2) python ./Kamodo-master/kamodo/flythrough/SingleSatelliteFlythrough.py CTIPe
(prints variable dictionary for given model)

(3) python ./Kamodo-master/kamodo/flythrough/SingleSatelliteFlythrough.py CTIPe 
    C:/Users/rringuet/Kamodo_WinDev1/CTIPe/Data/ 
(initializes files)

(4) python ./Kamodo-master/kamodo/flythrough/SingleSatelliteFlythrough.py CTIPe 
    C:/Users/rringuet/Kamodo_WinDev1/CTIPe/Data/ ['rho','TEC'] [1,0] 
    1426660000.0 -25. 90. 400. GDZ sph 0.02 
(calculates value at given location and vertical derivatives if requested)
syntax: python program_path/program_name.py model file_dir variable_list
    z_dependence dz sat_time sat_height sat_lat sat_lon model_dt high_res
    
    model: the model string chosen from the output of the first syntax choice.
        An integer representing the model name is also allowed. See result from
        first call version for more informaiton.
    file_dir: the directory above the ./Data/ directory where the model data is stored
    variable_list: the list of desired variables, chosen from the second syntax choice.
        Users may choose either the string or integer representation of the variable
        for the desired model. See the result of the second call type for more 
        information.
    dz: a list of 0s and 1s indicating which variables the vertical derivative 
        is desired (1) and which is not (0). The code will error if you request
        a vertical derivative for a variable that has no vertical dependence.
    sat_time: the UTC time in seconds
    c1: either the x (R_E) or longitude (deg) value
    c2: either the y (R_E) or latitude (deg) value
    c3: either the z (R_E), altitude (km), or radius (R_E) value
    coord_type: one of GDZ or 0, GEO or 1, GSM or 2, GSE or 3, SM or 4, GEI or 5, 
                MAG or 6, SPH or 7, RLL or 8 to describe your coordinates. Strings
                or integers allowed.
    coord_grid: either car (0) or sph (1) to describe your coordinates. Strings
                or integers allowed.
    high_res: the accuracy of the conversion of height to pressure level (meters).
        If not given, 20 meters will be used. This is optional since some models
        do not have variables that depend on pressure level.

Example times and variables for testing
--------------------------    
CTIPE: 1426660000.0, 'rho'
GITM: 1166051800.0, 'rho_n'
IRI: 1495844200.0, 'N_e'
SWMF_IE: 1533081600.0, 'Q_Joule'
TIEGCM: 974420200.0, 'rho'

"""  
#import numpy as np
from time import perf_counter
t0=perf_counter()
from numpy import array
from kamodo.flythrough.SF_utilities import Prepare_Files, Single_FlyAway, MW


def Prerun(model, file_dir, verbose=True):
    '''Determine vertical dependency and model time resolution of given model'''
    
    Prepare_Files(model, file_dir, call_type='single')
    if verbose:
        print('Full call syntax:')
        #print('python programpath/SingleSatelliteFlythrough.py model file_dir var_list z_list dz_list UTC_time height lat lon dt')
        print('python programpath/SingleSatelliteFlythrough.py model file_dir var_list dz_list UTC_time c1 c2 c3')
        print('If cartesian coordinates, then (c1,c2,c3)=(x,y,z). If spherical, then (c1,c2,c3)=(lon,lat,radius or altitude).')
    return #z_dependence, dt


#want to enable call of this from C++ for flexibility, so return only one value
#keep so users can call this if they have their own satellite trajectory data
def SingleModelFlythrough(model, file_dir, variable_list, dz, sat_time, c1, c2, 
                          c3, coord_type, coord_grid, high_res):  
    '''Call satellite flythrough wrapper specific to the model chosen.
    Parameters:    
        Name of model: model (Options: 'CTIPe', ...)
        Absolute path to where model data is stored: file_dir
        List of desired standardized variable names: variable_list
        Single float of satellite trajectory timestamps: sat_time
            (in number of seconds since 1970-01-01)
        Single float of satellite trajectory heights in meters: sat_height
        Single float of satellite trajectory latitudes in degrees: sat_lat
        Single float of satellite trajectory longitudes in degrees: sat_lon   
        List of vertical dependence for each variable in variable_list in same
            order as variables: z_dependence. Options for each variable: 'none', 'height', 'ilev'
            for no dependence on height, dependence on height, dependence on pressure
            altitude (ilev)
        List of booleans or 0's and 1's: dz. Indicates whether partial derivative
            of each variable with height in variable_list is requested. If
            z_dependence for that variable is 'none', then this term is ignored.
    Returns a dictionary with keys: sat_time, sat_height, sat_lat, sat_lon, 
    and keys naming the requested variables.
        sat_time is a single float in seconds since 1970-01-01.
        sat_height is a single float in meters.
        sat_lat and sat_lon are single floats in degrees.
        model variable keys are returned in the units printed out.
    ''' 
    
    var_list = variable_list.copy()  #avoid addition of H from CTIPe and similar
    var_units = MW.Var_units(model, var_list)
    results = Single_FlyAway(model, file_dir, variable_list, dz, sat_time, c1, 
                             c2, c3, coord_type, coord_grid, high_res) 
    
    #print results to a csv file
    file = open(file_dir+model+'_results.txt', 'w')
    for i in range(len(var_list)):
        if dz[i]==1 and var_list[i]+'_dz' in results.keys():  #if dz was determined for this variable, output both
            file.write(f"{var_list[i]}, {results[var_list[i]]:.20f}, {results[var_list[i]+'_dz']:.20f}, {var_units[var_list[i]]}")
        else:  #otherwise, just output the value
            file.write(f"{var_list[i]}, {results[var_list[i]]:.20f}, {var_units[var_list[i]]}")
        if i<len(var_list):
            file.write('\n' )
    file.close()
    print(results)
    
    return results  #not sure that than C++ can take more than one return variable

if __name__=='__main__':
    '''Code for translating cmd line call into action. Disregard any return values.'''
    t1=perf_counter()  
    from sys import argv  #for accepting command line inputs
    

    if len(argv)==2:  #first term is a string of a model name
        #python programpath/SingleSatelliteFlythrough CTIPe
        model = argv[1]
        MW.Model_Variables(model)  #asking possible model variables
    elif len(argv)==3:  #asking to prepare data
        #python programpath/SingleSatelliteFlythrough CTIPe C:/Users/rringuet/Kamodo_WinDev1/CTIPe/
        model = argv[1]
        file_dir = argv[2]   
        Prerun(model, file_dir)
    elif len(argv)>3:  #asking to run flythrough
        #python ./Kamodo-master/kamodo/satelliteflythrough/SingleSatelliteFlythrough_testing.py 
        #   CTIPe C:/Users/rringuet/Kamodo_WinDev1/CTIPe/ ['rho','TEC'] 
        #   [1,0] 1426637500.0 90. -25. 400. 'GDZ','sph'
        model = argv[1]   #CTIPe
        if len(model)==1: model = int(model)
        file_dir = argv[2]   #'C:/Users/rringuet/Kamodo_WinDev1/CTIPe/'
        temp_str = argv[3][1:-1].replace("'","").replace(' ','').replace('"','')
        variable_list = temp_str.split(',')   #['rho','N_n']
        if array([len(item) for item in variable_list]).max()==1:
            variable_list = array(variable_list, dtype=int)
        temp_str = argv[4][1:-1].replace("'","").replace(' ','').replace('"','')
        dz = list(array(temp_str.split(','),dtype=int))  #[1,0], later treated as a boolean
        sat_time = float(argv[5]) #1426637500.0
        c1 = float(argv[6]) #x[R_E] or lon[deg]
        c2 = float(argv[7]) #y[R_E] or lat[deg]
        c3 = float(argv[8]) #z[R_E] or radius[R_E] or altitude[km]
        coord_type = argv[9]  #'SPH', 'GDZ', etc
        #print(coord_type, len(coord_type), type(coord_type))
        if len(coord_type)==1: coord_type=int(coord_type)
        coord_grid = argv[10]  #'car' or 'sph'
        if len(coord_grid)==1: coord_type=int(coord_grid)
        if len(argv)>11: 
            high_res = float(argv[11])
        else:
            high_res = 20.      
        ta=perf_counter()    
        results = SingleModelFlythrough(model, file_dir, variable_list, 
                                        dz, sat_time, c1, c2, c3, coord_type,
                                        coord_grid, high_res)
        print(f'SingleModelFlythrough Time: {perf_counter()-ta:.5f}')   #0.305475 s
    else:
        MW.Choose_Model('')  #asking for a list of possible models
    print(f'Total Time: {perf_counter()-t1:.5f}')   #0.30689 s

if __name__!='__main__': print('import time:', perf_counter()-t0)