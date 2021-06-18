# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 14:57:55 2021
author: rringuet

Code to be called from other languages for a single trajectory point.
All desired model data should be in a single directory.


#Syntax of possible command line calls 
-------------------------------------
(1) python ./Kamodo-master/kamodo/satelliteflythrough/SingleSatelliteFlythrough.py
(prints possible models)

(2) python ./Kamodo-master/kamodo/satelliteflythrough/SingleSatelliteFlythrough.py CTIPe
(prints variable dictionary for given model)

(3) python ./Kamodo-master/kamodo/satelliteflythrough/SingleSatelliteFlythrough.py CTIPe 
    C:/Users/rringuet/Kamodo_WinDev1/CTIPe/ ['rho',TEC']
(prints vertical dependence per variable and model time resolution and initializes files)

(4) python ./Kamodo-master/kamodo/satelliteflythrough/SingleSatelliteFlythrough.py CTIPe 
    C:/Users/rringuet/Kamodo_WinDev1/CTIPe/ ['rho','TEC'] ['ilev','none'] [1,0] 
    1426637500.0 400. -25. 90. 450. 0.02
(calculates value at given location and vertical derivatives if requested)
syntax: python program_path/program_name.py model file_dir variable_list
    z_dependence dz sat_time sat_height sat_lat sat_lon model_dt high_res
    
    model: the model string chosen from the output of the first syntax choice
    file_dir: the directory above the ./Data/ directory where the model data is stored
    variable_list: the list of desired variables, chosen from the second syntax choice
    z_dependence: the list of vertical dependence strings given by the third syntax choice
    dz: a list of 0s and 1s indicating which variables the vertical derivative 
        is desired (1) and which is not (0). The code will error if you request
        a vertical derivative for a variable that has no vertical dependence.
    sat_time: the UTC time in seconds
    sat_height: the height of the satellite in kilometers
    sat_lat: the latitude of the satellite in degrees, assuming a range of -90 to +90 degrees
    sat_lon: the longitude of the satellite in degrees, assuming a range of 0 to 360 degrees
    model_dt: half of the time resolution of the model data, given by the output 
        of the third syntax choice. This value is used as a maximum allowed time 
        shift for times requested that are between the model data files.
    high_res: the accuracy of the conversion of height to pressure level (meters).
        If not given, 20 meters will be used. This is optional since some models
        do not have variables that depend on pressure level.

Example times for testing
--------------------------    
CTIPE: 1426660000.0
IRI: 1495945560.0
GITM: 1165968000.0
SWMF_IE: 1533081600.0
TIEGCM: 974264400.0

"""  
import numpy as np
import kamodo.satelliteflythrough.model_wrapper as MW
import kamodo.satelliteflythrough.wrapper_utilities as U


def ModelPrerun(model, file_dir, variable_list, verbose=True):
    '''Determine vertical dependency and model time resolution of given model'''
    
    z_dependence, dt = U.Single_Prerun(model, file_dir, variable_list)
    if verbose:
        print('Full call syntax:')
        print('python programpath/SingleSatelliteFlythrough.py model file_dir var_list z_list dz_list UTC_time height lat lon dt')
    return z_dependence, dt


#want to enable call of this from C++ for flexibility, so return only one value
#keep so users can call this if they have their own satellite trajectory data
def SingleModelFlythrough(model, file_dir, variable_list, z_dependence, dz, 
                          sat_time, sat_height, sat_lat, sat_lon, dt, high_res):  
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
    Var = MW.Model_Variables(model, return_dict=True)
    var_units = [value[-1] for key, value in Var.items() if value[0] in var_list]
    results = U.Single_FlyAway(model, file_dir, variable_list, z_dependence, dz, 
                                     sat_time, sat_height, sat_lat, sat_lon, dt,
                                     high_res) 
    
    #print results to a csv file
    file = open(file_dir+model+'_results.txt', 'w')
    for i in range(len(var_list)):
        if dz[i]==1:  #if dz was determined for this variable, output both
            file.write(f"{var_list[i]}, {results[var_list[i]]:.20f}, {results[var_list[i]+'_dz']:.20f}, {var_units[i]}")
        else:  #otherwise, just output the value
            file.write(f"{var_list[i]}, {results[var_list[i]]:.20f}, {var_units[i]}")
        if i<len(var_list):
            file.write('\n' )
    file.close()
    print(results)
    
    return results  #not sure that than C++ can take more than one return variable

if __name__=='__main__':
    '''Code for translating cmd line call into action. Disregard any return values.'''
    import sys  #for accepting command line inputs
    

    if len(sys.argv)==2:  #first term is a string of a model name
        #python programpath/SingleSatelliteFlythrough CTIPe
        model = sys.argv[1]
        MW.Model_Variables(model)  #asking possible model variables
    elif len(sys.argv)==4:  #asking for prerun
        #python programpath/SingleSatelliteFlythrough CTIPe C:/Users/rringuet/Kamodo_WinDev1/CTIPe/
        #    ['rho','N_n']
        model = sys.argv[1]
        file_dir = sys.argv[2]   
        temp_str = sys.argv[3][1:-1].replace("'","").replace(' ','').replace('"','')
        variable_list = temp_str.split(',')  
        z_dependence, dt = ModelPrerun(model, file_dir, variable_list)
    elif len(sys.argv)>4:  #asking to run flythrough
        #python ./Kamodo-master/kamodo/satelliteflythrough/SingleSatelliteFlythrough_testing.py 
        #   CTIPe C:/Users/rringuet/Kamodo_WinDev1/CTIPe/ ['rho','TEC'] ['ilev','none'] 
        #   [1,0] 1426637500.0 400. -25. 90. 450.
        model = sys.argv[1]   #CTIPe
        file_dir = sys.argv[2]   #'C:/Users/rringuet/Kamodo_WinDev1/CTIPe/'
        temp_str = sys.argv[3][1:-1].replace("'","").replace(' ','').replace('"','')
        variable_list = temp_str.split(',')   #['rho','N_n']
        temp_str = sys.argv[4][1:-1].replace("'","").replace(' ','').replace('"','')
        z_dependence = temp_str.split(',')  #['ilev','ilev']
        temp_str = sys.argv[5][1:-1].replace("'","").replace(' ','').replace('"','')
        dz = list(np.array(temp_str.split(','),dtype=int))  #[1,0], later treated as a boolean
        sat_time = float(sys.argv[6]) #1426637500.0
        sat_height = float(sys.argv[7]) #400. (in km)
        sat_lat = float(sys.argv[8]) #-25.
        sat_lon = float(sys.argv[9]) #90.
        model_dt = float(sys.argv[10])  #450.
        if len(sys.argv)>11: 
            high_res = float(sys.argv[11])
        else:
            high_res = 20.      
        results = SingleModelFlythrough(model, file_dir, variable_list, 
                                        z_dependence, dz, sat_time, sat_height, 
                                        sat_lat, sat_lon, model_dt, high_res)
    else:
        MW.Model_Variables('')  #asking for a list of possible models
  