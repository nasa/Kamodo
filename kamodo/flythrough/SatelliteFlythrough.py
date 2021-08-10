# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 14:57:55 2021
author: rringuet

Satellite Flythrough code to call for any model.
All desired model data for the desired model should be in a single directory.

Example calls from the command line:
    Get list of models and function calls possible:  python ./SatelliteFlythrough_new.py
    Get call syntax for one of the functions: python ./SatelliteFlythrough_new.py MyFlight
    Get information on what variables are available for a given model:
        python ./SatelliteFlythrough_new.py CTIPe
    For funtion calls, all of the parameters in the function call should be given
        after the name of the function in the example above, even if the default
        value is desired


"""
import numpy as np
from os.path import basename
import kamodo.flythrough.wrapper_output as WO
import kamodo.flythrough.SF_utilities as U


def SatelliteTrajectory(dataset, start_ts, stop_ts, coord_type='GEO',
                        verbose=False):
    '''Retrieve and return satellite trajectory from HAPI/CDAWeb
    Parameters:
    ----------
    dataset: name of the satellite data set to pull trajectory from
    start_ts: utc timestamp for start of desired time interval
    stop_ts: utc timestamp for end of desired time interval    
    coord_type: Pick from GEO, GSM, GSE, or SM
    verbose: Set to true to be overwhelmed with information.
    
    Coordinates are retrieved on a cartesian grid.
    '''
    from kamodo.readers.hapi import HAPI


    #convert from utc timestamps to isoformt
    start = U.ts_to_ISOstring(start_ts)
    stop = U.ts_to_ISOstring(stop_ts)
    
    #convert from integer input of coord_type to string
    coord_type, coord_grid = U.MW.convert_coordnames(coord_type, 'car')

    #check input coord_type
    if coord_type not in ['GEO','GSM','GSE','SM']: 
        raise AttributeError(f'Coordinate type {coord_type} not available. '+
              'Pick from GEO, GSM, GSE, or SM.')
    else:
        parameters = 'X_'+coord_type+',Y_'+coord_type+',Z_'+coord_type

    #retrieve satellite trajectory
    server = 'http://hapi-server.org/servers/SSCWeb/hapi'  #for coordinate data
    hapi = HAPI(server, dataset, parameters, start, stop, verbose=verbose)
    satellite_dict = {'sat_time': hapi.tsarray,  #utc timestamps
                      'c1': hapi.variables[parameters.split(',')[0]]['data'],  #x coord
                      'c2': hapi.variables[parameters.split(',')[1]]['data'],  #y coord
                      'c3': hapi.variables[parameters.split(',')[2]]['data']}  #z coord

    print(f'Attribute/Key names of return dictionary: {satellite_dict.keys()}')
        
    return satellite_dict, coord_type, 'car'

def SampleTrajectory(start_time, stop_time, max_lat=65., min_lat=-65.,
                     lon_perorbit=363., max_height=450., min_height=400., 
                     p=0.01, n=2.):
    '''
    Given start and stop times in timestamp form, return a test satellite trajectory.
    Parameters:
    ----------
        start_time: utc timestamp in seconds for start
        stop_time: utc timestamp in seconds for stop
        max_lat: maximum latitude for sample trajectory, in degrees (default=65.)
        min_lat: minimum latitude for sample trajectory, in degrees (default=-65.)
        lon_perorbit: the degrees of longitude per about 90 minute orbit 
            (set less than 360 for precession forward in longitude, set less 
            than 360 for precession backwards) (default=363.)
        max_height: maximum starting height of orbit in km (default=450.)
        min_height: minimum starting height of orbit in km (default=400.)
        p: a rough precession variable, applied as an overall height decrease 
            as a percentage of the min_height value: p =  (default=0.01).  
        n: the time cadence of the sample trajectory generated (default = 2 seconds)
    Returns a dictionary with keys: sat_time, c1, c2, and c3.
        sat_time is an array in seconds since 1970-01-01.
        (c1,c2,c3) = (lon, lat, alt) in (deg,deg,km) in the 'GDZ', 'sph' 
        coordinate system in SpacePy.
    '''
    
    #determine basic parameters
    orbit_seconds = int(90.*60./float(n))  #determine number of samples per 90min orbit
    n_orbits = (stop_time-start_time)/float(orbit_seconds*n) #orbits of 90 min each
    h_scale, h_offset = (max_height-min_height)/2., np.mean([max_height,min_height])
    lat_scale, lat_offset = (max_lat-min_lat)/2., np.mean([max_lat,min_lat])
    time_left = (stop_time-start_time)/float(n)-int(n_orbits)*orbit_seconds
    
    #create orbital tracks 
    pi_arr = np.linspace(0.,2.*np.pi,orbit_seconds)  
    lat, height = np.tile(np.cos(pi_arr), int(n_orbits)), np.tile(np.sin(pi_arr), int(n_orbits))
    if time_left>0:
        lat = np.append(lat, np.cos(pi_arr[0:int(time_left)]))  #add partial last orbit
        height = np.append(height, np.sin(pi_arr[0:int(time_left)]))
    lon = np.linspace(0.,float(lon_perorbit)*n_orbits,int((stop_time-start_time)/float(n)))
    while max(lon)>360.: lon[np.where(lon>360.)[0]]-=360.
    while max(lon)<0.: lon[np.where(lon<0.)[0]]+=360.
    height = height*h_scale+h_offset-np.linspace(0.,p,int((stop_time-start_time)/float(n)))*min_height
    
    #store results in dictionary to return
    sample_dict={'sat_time': np.linspace(start_time,stop_time,int((stop_time-start_time)/float(n))),
                 'c1': lon-180., 'c2': lat*lat_scale+lat_offset, 'c3': height}   
    
    print(f'Attribute/Key names of return dictionary: {sample_dict.keys()}')
    print('(c1,c2,c3) = (lon, lat, alt) in (deg,deg,km) in the GDZ, sph coordinate system.'+\
          'sat_time contains the utc timestamps.')    
    return sample_dict, 'GDZ', 'sph'


#want to enable call of this from C++ for flexibility, so return only one value
#keep so users can call this if they have their own satellite trajectory data
def ModelFlythrough(model, file_dir, variable_list, sat_time, c1, c2, c3, 
                    coord_type, coord_grid, high_res=20., verbose=False, 
                    csv_output='', plot_output=''):  
    '''Call satellite flythrough wrapper specific to the model chosen.
    Parameters:   
    ------------
    model: 'CTIPe','IRI', ...
    file_dir: complete path to where model data files are stored
    variable_list: List of standardized variable names. Corresponding integers 
        are allowed. See model variable output for details.
    sat_time: a numpy array of the utc timestamps
    c1, c2, c3: numpy arrays of the positions correlating to the utc timestamps
        (c1, c2, c3) should be (x,y,z) in R_E for cartesian coordinates, and (lon, lat, 
        radius (R_E) or altitude (km)) for spherical coordinates. 
    coord_type: one of 'GDZ', 'GEO', 'GSM', 'GSE', 'SM', 'GEI', 'MAG', 'SPH', 'RLL'
        integers also allowed with 'GDZ'=0 and so on
    coord_grid: either 'car' or 'sph' (0 or 1). Note that not all combinations 
        make sense (e.g. 'SPH' and 'car') and are not allowed.    
    high_res: the accuracy of the conversion from radius or altitude to pressure
        level. Ignore if no conversion is needed for the variable(s) selected.        
    csv_output: complete path pluts filename (without the .csv) for the file to
        write the results to.
    plot_output: complete path pluts file naming convention (without the .html)
        for the file to write the plots to.
    verbose: Set to true to be overwhelmed with information.        
    ''' 

    #if input types are lists, correct to be numpy arrays (important for calling from C++)
    if isinstance(sat_time, list): sat_time = np.array(sat_time)
    if isinstance(c1, list): c1 = np.array(c1)
    if isinstance(c2, list): c2 = np.array(c2)
    if isinstance(c3, list): c3 = np.array(c3)
        
    #if model is given as an integer, then convert to a string
    model = U.MW.convert_model_string(model)
    
    #if variable_list is a list of integers, convert to standard names for given model
    variable_list = U.MW.convert_variablenames(model, variable_list)
    
    #convert integer coordinate names or grids to strings
    coord_type, coord_grid = U.MW.convert_coordnames(coord_type, coord_grid)
    
    #retrieve coordinate and results units
    coord_units = U.MW.coord_units(coord_type, coord_grid)
    results_units = U.MW.Var_units(model, variable_list)    
    for key in coord_units: results_units[key] = coord_units[key]
    print(results_units)
    
    #prepare files for run
    U.Prepare_Files(model, file_dir)
    
    #get interpolated results
    #coord_type should be one of SpacePy's coordinates, coord_grid is either 'sph' or 'car'
    results = U.Model_SatelliteFlythrough(model, file_dir, variable_list, 
                                sat_time, c1, c2, c3, 
                                coord_type, coord_grid, high_res,
                                verbose=verbose)  
    
    print('Done.\n')
    if verbose: 
        print(f'Units from the {model} model by variable name:\n{results_units}')
        print(f'Dictionary key names in results:\n{results.keys()}')
        print('The units of the trajectory variables are unchanged from the inputs.')
        
    if csv_output!='':
        #correct input filename
        if model not in basename(csv_output):
            file_dir = csv_output.split(basename(csv_output))[0]
            csv_output = file_dir+model+'_'+basename(csv_output)
        csv_filename = WO.SFdata_tocsv(csv_output, '', model, results, results_units)
        print(f"Output saved in {csv_filename}.")  #no access to model filenames
        
    if plot_output!='':
        print('Generating interactive plots...')
        
        #correct input filename and split into useful pieces
        file_prefix = basename(plot_output)
        file_dir = plot_output.split(file_prefix)[0]
        if model not in file_prefix:
            file_dir+=model+'_'     
        file_prefix+='_'            
        
        #generate and save plots without displaying
        from kamodo.flythrough.plots import SatPlot4D
        #presentation options: all, day, hour, minute, N, orbitE, orbitM
        for var in variable_list:
            SatPlot4D(var,results['utc_time'],results['c1'],results['c2'],results['c3'],
                      results[var],results_units[var],
                      coord_type, coord_grid, 'GEO','all',model,body='black', 
                      divfile=file_dir+file_prefix+var+'_3D.html', displayplot=False)
            SatPlot4D(var,results['utc_time'],results['c1'],results['c2'],
                      results['c3'],results[var],results_units[var],
                  coord_type, coord_grid, 'GEO','all',model,type='1D',
                  divfile=file_dir+file_prefix+var+'_1D.html', displayplot=False)
    
    return results  #not sure that than C++ can take more than one return variable 

def FakeFlight(start_time, stop_time, model, file_dir, variable_list, max_lat=65., 
               min_lat=-65., lon_perorbit=363., max_height=450., min_height=400., 
               p=0.01, n=2., high_res=20., verbose=False, csv_output='',
               plot_output=''):
    '''Generates a sample trajectory and then flies that trajectory through the 
    model data chosen.
     
    Parameters: 
        start_time: utc timestamp in seconds for start
        stop_time: utc timestamp in seconds for stop
        model: CTIPe, IRI, .... (integers allowed)
        file_dir: complete path to where model data is stored
        variable_list: list of standardized variable names desired. Integers allowed.
        max_lat: maximum latitude for sample trajectory, in degrees (default=65.)
        min_lat: minimum latitude for sample trajectory, in degrees (default=-65.)
        lon_perorbit: the degrees of longitude per about 90 minute orbit 
            (set less than 360 for precession forward in longitude, set less 
            than 360 for precession backwards) (default=363.)
        max_height: maximum starting height of orbit in km (default=450.)
        min_height: minimum starting height of orbit in km (default=400.)
        p: a rough precession variable, applied as an overall height decrease 
            as a percentage of the min_height value: p =  (default=0.01).  
        n: the time cadence of the sample trajectory generated (default = 2 seconds)
        high_res: the resolution of the height conversion to pressure level
            in units of km
        csv_output: complete path pluts filename (without the .csv) for the file to
            write the results to.
        plot_output: complete path pluts file naming convention (without the .html)
            for the file to write the plots to.

    ''' 
    
    #generate a sample satellite trajectory
    sat_dict, coord_type, coord_grid = SampleTrajectory(start_time, stop_time,
                                max_lat=max_lat, min_lat=min_lat, lon_perorbit=lon_perorbit, 
                                max_height=max_height, min_height=min_height, p=p, n=n)
    
    #call satellite flythrough code
    results = ModelFlythrough(model, file_dir, variable_list, sat_dict['sat_time'], 
                              sat_dict['c1'], sat_dict['c2'], sat_dict['c3'],
                              coord_type, coord_grid, high_res=high_res,
                              verbose=verbose, csv_output=csv_output,
                              plot_output=plot_output)    
    return results    

def RealFlight(dataset, start, stop, model, file_dir, variable_list, coord_type='GEO',
               csv_output='', plot_output='', high_res=20., verbose=False):
    '''
    Retrieves the trajectory for the satellite requested and then flies that
    trajectory through the model data requested.

    dataset: name of the satellite data set to pull trajectory from
    start: utc timestamp for start of desired time interval
    stop: utc timestamp for end of desired time interval
    model: 'CTIPe','IRI', ...
    file_dir: complete path to where model data files are stored
    variable_list: List of standardized variable names. Corresponding integers 
        are allowed. See model variable output for details.
    csv_output: complete path pluts filename (without the .csv) for the file to
        write the results to.
    plot_output: complete path pluts file naming convention (without the .html)
        for the file to write the plots to.
    high_res: the accuracy of the conversion from radius or altitude to pressure
        level. Ignore if no conversion is needed for the variable(s) selected.
    verbose: Set to true to be overwhelmed with information.
    '''
    #retrieve satellite trajectory from HAPI/CDAWeb
    sat_dict, coord_type, coord_grid = SatelliteTrajectory(dataset, start, stop, 
                                   coord_type, verbose=verbose)
    #call satellite flythrough code
    results = ModelFlythrough(model, file_dir, variable_list, sat_dict['sat_time'], 
                              sat_dict['c1'], sat_dict['c2'], sat_dict['c3'],
                              coord_type, coord_grid, csv_output=csv_output,
                              plot_output=plot_output, high_res=high_res, 
                              verbose=verbose)  
    return results    


def MyFlight(traj_file, file_type, coord_type, coord_grid, model, file_dir, variable_list,
               csv_output='', plot_output='', high_res=20., verbose=False):
    '''Read in a trajectory from a file, then fly through the model data selected.
    
    traj_file: complete path and filename for file containing trajectory data.
    file_type: one of 'cdf', 'csv', or 'ascii' of the format required
    coord_type: one of 'GDZ', 'GEO', 'GSM', 'GSE', 'SM', 'GEI', 'MAG', 'SPH', 'RLL'
        integers also allowed with 'GDZ'=0 and so on
    coord_grid: either 'car' or 'sph' (0 or 1). Note that not all combinations 
        make sense (e.g. 'SPH' and 'car') and are not allowed.
    model: 'CTIPe', 'IRI', ...  
    file_dir: complete path to model data files
    variable_list: List of standardized variable names. Corresponding integers 
        are allowed. See model variable output for details.
    csv_output: complete path pluts filename (without the .csv) for the file to
        write the results to.
    plot_output: complete path pluts file naming convention (without the .html)
        for the file to write the plots to.
    high_res: the accuracy of the conversion from radius or altitude to pressure
        level. Ignore if no conversion is needed for the variable(s) selected.
    verbose: Set to true to be overwhelmed with information.
    '''

    #read in trajectory from file into dictionary
    if file_type=='cdf':
        traj_data = WO.SFcdf_reader(traj_file)
    elif file_type=='csv':
        traj_data=WO.SFcsv_reader(traj_file)
    elif file_type=='ascii':
        traj_data=WO.SFascii_reader(traj_file)
        
    #figure out key for time data
    for key in traj_data:
        if 'time' in key: 
            time_key = key
            break
    
    #call satellite flythrough code
    results = ModelFlythrough(model, file_dir, variable_list, traj_data[time_key]['data'], 
                              traj_data['c1']['data'], traj_data['c2']['data'], 
                              traj_data['c3']['data'], coord_type, coord_grid, 
                              csv_output=csv_output, plot_output=plot_output, 
                              high_res=high_res, verbose=verbose)      
        
    return results


#allow calls from the command line
#these calls require all variables be given
if __name__=='__main__':
    from sys import argv 
    
    
    #print info if called without arguments
    if len(argv)==2:
        if argv[1]=='FakeFlight':
            help(FakeFlight)
        elif argv[1]=='RealFlight':
            help(RealFlight)
        elif argv[1]=='MyFlight':
            help(MyFlight)
        else:
            U.MW.Model_Variables(argv[1])
    elif len(argv)>2:
        if argv[1]=='FakeFlight':  #gather variables and call FakeFlight
            start_time = int(argv[2])
            stop_time = int(argv[3])
            if len(argv[4])==1:
                model = int(argv[4])
            else:
                model = argv[4]
            file_dir = argv[5]
            temp_str = argv[6][1:-1].replace("'","").replace(' ','').replace('"','')
            variable_list = temp_str.split(',')   #['rho','N_n']            
            max_lat = float(argv[7])
            min_lat = float(argv[8])
            lon_perorbit = float(argv[9])
            max_height = float(argv[10])
            min_height = float(argv[11])
            p = float(argv[12])
            n = float(argv[13])
            high_res = float(argv[14])
            csv_output = argv[15]
            plot_output = argv[16]
            
            #check input
            print(f'\nstart_time: {start_time}, \nstop_time: {stop_time},'+\
                  f'\nmodel: {model}, \nfile_dir: {file_dir},'+\
                  f'\nvariable_list: {variable_list}, \nmax_lat: {max_lat}, \nmin_lat: {min_lat},'+\
                  f'\nlon_perorbit: {lon_perorbit}, \nmax_height: {max_height}, \nmin_height {min_height},'+\
                  f'\np: {p}, n: {n}, \ncsv_output: {csv_output}, \nplot_output: {plot_output}\n')

            results = FakeFlight(start_time, stop_time, model, file_dir, variable_list, 
                                 max_lat=max_lat, min_lat=min_lat, lon_perorbit=lon_perorbit, 
                                 max_height=max_height, min_height=min_height,
                                 p=p, n=n, high_res=high_res, verbose=False, 
                                 csv_output=csv_output, plot_output=plot_output)

        elif argv[1]=='RealFlight':  #gather variables and call RealFlight
            dataset = argv[2]
            start = int(argv[3])
            stop = int(argv[4])
            if len(argv[5])==1:
                model = int(argv[5])
            else:
                model = argv[5]
            file_dir = argv[6]
            temp_str = argv[7][1:-1].replace("'","").replace(' ','').replace('"','')
            variable_list = temp_str.split(',')   #['rho','N_n']              
            if len(argv[8])==1:
                coord_type = int(argv[8])
            else:
                coord_type = argv[8]
            csv_output = argv[9]
            plot_output = argv[10]
            high_res = float(argv[11])            
            
            #check input
            print(f'\ndataset: {dataset}, \nstart: {start}, \nstop: {stop},'+\
                  f'\nmodel: {model}, \nfile_dir: {file_dir},'+\
                  f'\nvariable_list: {variable_list}, \ncoord_type: {coord_type},'+\
                  f'\ncsv_output: {csv_output}, \nplot_output: {plot_output},'+\
                  f'\nhigh_res: {high_res}\n')

            results = RealFlight(dataset, start, stop, model, file_dir, variable_list, 
                                 coord_type=coord_type, csv_output=csv_output, 
                                 plot_output=plot_output, high_res=high_res)

        elif argv[1]=='MyFlight':  #gather variables and call MyFlight
            traj_file = argv[2]
            file_type = argv[3]
            if len(argv[4])==1:
                coord_type = int(argv[4])
            else:
                coord_type = argv[4]
            if len(argv[5])==1:
                coord_grid = int(argv[5])
            else:
                coord_grid = argv[5]                
            if len(argv[6])==1:
                model = int(argv[6])
            else:
                model = argv[6]
            file_dir = argv[7]
            temp_str = argv[8][1:-1].replace("'","").replace(' ','').replace('"','')
            variable_list = temp_str.split(',')   #['rho','N_n']              
            csv_output = argv[9]
            plot_output = argv[10]
            high_res = float(argv[11])
            
            #check inputs
            print(f'\ntraj_file: {traj_file}, \nfile_type: {file_type},'+\
                  f'\ncoord_type: {coord_type}, \ncoord_grid: {coord_grid},'+\
                  f'\nmodel: {model}, \nfile_dir: {file_dir},'+\
                  f'\nvariable_list: {variable_list}, \ncsv_output: {csv_output},'+\
                  f'\nplot_output: {plot_output}, \nhigh_res: {high_res}\n')
            
            results = MyFlight(traj_file, file_type, coord_type, coord_grid, 
                               model, file_dir, variable_list, csv_output=csv_output, 
                               plot_output=plot_output, high_res=high_res)
        else:
            print('Call signature not recognized.')
    else:
        print('\nPossible call types (first argument): FakeFlight, RealFlight, MyFlight')
        print('Use the call type as the first input to get call syntax.\n')
        U.MW.Choose_Model('')  #asking for a list of possible models
        
        
