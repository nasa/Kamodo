# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 14:57:55 2021
author: rringuet

Satellite Flythrough code to call for any model.
All desired model data for the desired model output run should be in a single directory.

Example calls from the command line:
    Get list of models and function calls possible:  python ./SatelliteFlythrough_new.py
    Get call syntax for one of the functions: python ./SatelliteFlythrough_new.py MyFlight
    Get information on what variables are available for a given model:
        python ./SatelliteFlythrough_new.py CTIPe
    For funtion calls, all of the parameters in the function call should be given
        after the name of the function in the example above, even if the default
        value is desired

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
The current links to SpacePy's coordinate documentation and wrapped conversion functions are:
    https://spacepy.github.io/autosummary/spacepy.coordinates.Coords.html
    http://svn.code.sf.net/p/irbem/code/trunk/manual/user_guide.html
"""
import numpy as np
from os.path import basename
import kamodo_ccmc.flythrough.wrapper_output as WO
import kamodo_ccmc.flythrough.SF_utilities as U


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
    See kamodo_ccmc.flythrough.utils.ConvertCoord for info on the coordinate systems.
    '''
    from kamodo_ccmc.readers.hapi import HAPI


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
        sat_time is an array in UTC seconds since 1970-01-01.
        (c1,c2,c3) = (lon, lat, alt) in (deg,deg,km) in the 'GDZ', 'sph' 
        coordinate system in SpacePy. See kamodo_ccmc.flythrough.utils.ConvertCoord 
        for more info on the coordinate systems.
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
          '\nsat_time contains the utc timestamps.')    
    return sample_dict, 'GDZ', 'sph'


#want to enable call of this from C++ for flexibility, so return only one value
#keep so users can call this if they have their own satellite trajectory data
def ModelFlythrough(model, file_dir, variable_list, sat_time, c1, c2, c3, 
                    coord_type, coord_grid, high_res=20., verbose=False, 
                    output_type='', output_name='', plot_output='', 
                    plot_coord='GEO', _print_units=True):  
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
        level. Ignore if no conversion is needed for the variable(s) selected. Default is 20.
    output_type: One of 'csv' for comma separated output, 'cdf4' for a netCDF4 
        output file, or 'txt' for a tab-separated text file.
    output_name: complete path with filename (without the extension) for the file to
        write the results to.
    plot_output: complete path pluts file naming convention (without the .html)
        for the file to write the plots to.
    plot_coord: one of 'GDZ', 'GEO', 'GSM', 'GSE', 'SM', 'GEI', 'MAG'
        integers also allowed with 'GDZ'=0 and so on. Indicates the coordinate
        system the plot will be generated in. Only plots in cartesian coordinates
        systems are supported, so 'SPH' and 'RLL' are not accepted. Default is 
        'GEO'.
    verbose: Set to true to be overwhelmed with information.     
    
    Returns a dictionary with keys: 'utc_time', 'c1', 'c2', 'c3', and 'net_idx'
        - utc_time is an array in UTC seconds since 1970-01-01 of the given
            timestamps with any occuring outside of the model data removed.
        - 'c1', 'c2', and 'c3' are arrays of the given coordinate values for each 
            surviving timestamp
        - 'net_idx' is the original index value of the surviving timestamps. 
            This is kept for easier comparison with the original dataset.
        - additional keys are included for each variable and label an array of the
            values of the indicated variable for each time+spatial coordinate given
        - The units of each array in the returned dictionary are printed to 
            the screen.
            
    See kamodo_ccmc.flythrough.utils.ConvertCoord for info on the coordinate systems.
    ''' 

    #if input types are lists, correct to be numpy arrays (important for calling from C++)
    if isinstance(sat_time, list): sat_time = np.array(sat_time)
    if isinstance(c1, list): c1 = np.array(c1)
    if isinstance(c2, list): c2 = np.array(c2)
    if isinstance(c3, list): c3 = np.array(c3)
    
    #give error if unknown output type given BEFORE running flythrough
    if output_type not in ['cdf4','csv','txt','', ' ']:  #allow empty strings
        raise AttributeError('Output file type not recognized. Must be one of'+\
                            ' cdf4, csv, or txt.')    
        
    #if model is given as an integer, then convert to a string
    model = U.MW.convert_model_string(model)
    
    #if variable_list is a list of integers, convert to standard names for given model
    variable_list = U.MW.convert_variablenames(model, variable_list)
    
    #convert integer coordinate names or grids to strings
    coord_type, coord_grid = U.MW.convert_coordnames(coord_type, coord_grid)
    
    #prepare files for run
    U.Prepare_Files(model, file_dir)
    
    #get interpolated results
    #coord_type should be one of SpacePy's coordinates, coord_grid is either 'sph' or 'car'
    results = U.Model_SatelliteFlythrough(model, file_dir, variable_list, 
                                sat_time, c1, c2, c3, 
                                coord_type, coord_grid, high_res,
                                verbose=verbose)  
    
    
    #remove requested variables not found in the data
    var_list = [key for key in results.keys() if key not in ['utc_time','c1','c2','c3','net_idx']]
    
    #retrieve coordinate and results units
    coord_units = U.MW.coord_units(coord_type, coord_grid)
    results_units = U.MW.Var_units(model, var_list)    
    for key in coord_units: results_units[key] = coord_units[key]
    if _print_units: print(results_units)    
    
    if verbose: 
        print(f'Units from the {model} model by variable name:\n{results_units}')
        print(f'Dictionary key names in results:\n{results.keys()}')
        print('The units of the trajectory variables are unchanged from the inputs.')
        
    if output_type!='':
        #correct input filename
        if model not in basename(output_name):
            output_file_dir = output_name.split(basename(output_name))[0]
            output_name = output_file_dir+model+'_'+basename(output_name)
            
        #retrieve file names/patterns for output
        file_times = U.MW.File_Times(model, file_dir, print_output=False)
        filenames=[]
        for key in file_times.keys(): filenames.append(file_times[key][0])
        
        #perform output type desired
        if output_type=='csv':
            output_filename = WO.SFdata_tocsv(output_name, filenames, model, results, 
                                              results_units, coord_type, coord_grid)
        elif output_type=='cdf4':
            output_filename= WO.SFdata_tocdf(output_name, filenames, model, results, 
                                             results_units, coord_type, coord_grid)
        elif output_type=='txt':
            output_filename = WO.SFdata_toascii(output_name, filenames, model, results, 
                                                results_units, coord_type, coord_grid)
        print(f"Output saved in {output_filename}.")  #no access to model filenames
        
    if plot_output!='':
        print('Generating interactive plots...')
        
        #correct input filename and split into useful pieces
        plot_file_prefix = basename(plot_output)
        plot_file_dir = plot_output.split(plot_file_prefix)[0]
        if model not in file_prefix:
            plot_file_dir+=model+'_'     
        plot_file_prefix+='_'            
        
        #check plot_coord variable, convert from integer and prevent plotting errors
        plot_coord, plot_grid = U.MW.convert_coordnames(plot_coord, 'car')
        if plot_coord in ['SPH','RLL']:
            raise AttributeError('Plots can only be requested in coordinate '+\
                                 'grids where cartesian coordinates are supported.'+\
                                 " The 'SPH' and 'RLL' coordinate systems "+
                                 ' do not support cartesian grids so are not allowed.')
        
        #generate and save plots without displaying
        from kamodo_ccmc.flythrough.plots import SatPlot4D
        #presentation options: all, day, hour, minute, N, orbitE, orbitM
        for var in var_list:
            SatPlot4D(var,results['utc_time'],results['c1'],results['c2'],results['c3'],
                      results[var],results_units[var],
                      coord_type, coord_grid, plot_coord,'all',model,body='black', 
                      divfile=plot_file_dir+plot_file_prefix+var+'_3D.html', displayplot=False)
            SatPlot4D(var,results['utc_time'],results['c1'],results['c2'],
                      results['c3'],results[var],results_units[var],
                  coord_type, coord_grid, plot_coord,'all',model,type='1D',
                  divfile=plot_file_dir+plot_file_prefix+var+'_1D.html', displayplot=False)
    
    return results  #not sure that than C++ can take more than one return variable 

def FakeFlight(start_time, stop_time, model, file_dir, variable_list, max_lat=65., 
               min_lat=-65., lon_perorbit=363., max_height=450., min_height=400., 
               p=0.01, n=2., high_res=20., verbose=False, output_type='',
               output_name='', plot_output='', plot_coord='GEO'):
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
            in units of km. Default is 20.
        output_type: One of 'csv' for comma separated output, 'cdf4' for a netCDF4 
            output file, or 'txt' for a tab-separated text file.
        output_name: complete path with filename (without the extension) for the file to
            write the results to.
        plot_output: complete path pluts file naming convention (without the .html)
            for the file to write the plots to.
        plot_coord: one of 'GDZ', 'GEO', 'GSM', 'GSE', 'SM', 'GEI', 'MAG'
            integers also allowed with 'GDZ'=0 and so on. Indicates the coordinate
            system the plot will be generated in. Only plots in cartesian coordinates
            systems are supported, so 'SPH' and 'RLL' are not accepted. Default is 
            'GEO'.
            
    Returns a dictionary with keys: 'utc_time', 'c1', 'c2', 'c3', and 'net_idx'
    - utc_time is an array in UTC seconds since 1970-01-01 of the generated
        timestamps with any occuring outside of the model data removed.
    - 'c1', 'c2', and 'c3' are arrays of the given coordinate values for each 
        surviving timestamp
    - 'net_idx' is the original index value of the surviving timestamps. 
        This is kept for easier comparison with the original dataset.
    - additional keys are included for each variable and label an array of the
        values of the indicated variable for each time+spatial coordinate given
    - The units of each array in the returned dictionary are printed to 
        the screen.

    See kamodo_ccmc.flythrough.utils.ConvertCoord for info on the coordinate systems.
    ''' 
    
    #generate a sample satellite trajectory
    sat_dict, coord_type, coord_grid = SampleTrajectory(start_time, stop_time,
                                max_lat=max_lat, min_lat=min_lat, lon_perorbit=lon_perorbit, 
                                max_height=max_height, min_height=min_height, p=p, n=n)
    
    #call satellite flythrough code
    results = ModelFlythrough(model, file_dir, variable_list, sat_dict['sat_time'], 
                              sat_dict['c1'], sat_dict['c2'], sat_dict['c3'],
                              coord_type, coord_grid, high_res=high_res,
                              verbose=verbose, output_type=output_type, 
                              output_name=output_name, plot_output=plot_output,
                              plot_coord=plot_coord)    
    return results    

def RealFlight(dataset, start, stop, model, file_dir, variable_list, coord_type='GEO',
               output_type='', output_name='', plot_output='', plot_coord='GEO',
               high_res=20., verbose=False):
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
    output_type: One of 'csv' for comma separated output, 'cdf4' for a netCDF4 
        output file, or 'txt' for a tab-separated text file.
    output_name: complete path with filename (without the extension) for the file to
        write the results to.
    plot_output: complete path pluts file naming convention (without the .html)
        for the file to write the plots to.
    high_res: the accuracy of the conversion from radius or altitude to pressure
        level. Ignore if no conversion is needed for the variable(s) selected. Default is 20.
    plot_coord: one of 'GDZ', 'GEO', 'GSM', 'GSE', 'SM', 'GEI', 'MAG'
        integers also allowed with 'GDZ'=0 and so on. Indicates the coordinate
        system the plot will be generated in. Only plots in cartesian coordinates
        systems are supported, so 'SPH' and 'RLL' are not accepted. Default is 
        'GEO'.        
    verbose: Set to true to be overwhelmed with information.
    
    Returns a dictionary with keys: 'utc_time', 'c1', 'c2', 'c3', and 'net_idx'
        - utc_time is an array in UTC seconds since 1970-01-01 of the satellite
            timestamps with any occuring outside of the model data removed.
        - 'c1', 'c2', and 'c3' are arrays of the given coordinate values for each 
            surviving timestamp
        - 'net_idx' is the original index value of the surviving timestamps. 
            This is kept for easier comparison with the original dataset.
        - additional keys are included for each variable and label an array of the
            values of the indicated variable for each time+spatial coordinate given
        - The units of each array in the returned dictionary are printed to 
            the screen.
            
    See kamodo_ccmc.flythrough.utils.ConvertCoord for info on the coordinate systems.
    '''
    
    #retrieve satellite trajectory from HAPI/CDAWeb
    sat_dict, coord_type, coord_grid = SatelliteTrajectory(dataset, start, stop, 
                                   coord_type=coord_type, verbose=verbose)
    #call satellite flythrough code
    results = ModelFlythrough(model, file_dir, variable_list, sat_dict['sat_time'], 
                              sat_dict['c1'], sat_dict['c2'], sat_dict['c3'],
                              coord_type, coord_grid, output_type=output_type,
                              output_name=output_name, plot_output=plot_output, 
                              plot_coord=plot_coord, high_res=high_res, verbose=verbose)  
    return results    


def MyFlight(traj_file, file_type, model, file_dir, 
             variable_list, output_type='', output_name='', plot_output='', 
             plot_coord='GEO', high_res=20., verbose=False):
    '''Read in a trajectory from a file, then fly through the model data selected.
    
    traj_file: complete path and filename for file containing trajectory data.
    file_type: one of 'cdf4' for netCDF4 files, 'csv' for comma-separated files, 
        or 'txt' for a tab-separated text file. Indicates the format of the input
        trajectory file. 
    coord_type: one of 'GDZ', 'GEO', 'GSM', 'GSE', 'SM', 'GEI', 'MAG', 'SPH', 'RLL'
        integers also allowed with 'GDZ'=0 and so on
    coord_grid: either 'car' or 'sph' (0 or 1). Note that not all combinations 
        make sense (e.g. 'SPH' with 'car') and are not allowed.
    model: 'CTIPe', 'IRI', ...  
    file_dir: complete path to model data files
    variable_list: List of standardized variable names. Corresponding integers 
        are allowed. See model variable output for details.
    output_type: One of 'csv' for comma separated output, 'cdf4' for a netCDF4 
        output file, or 'txt' for a tab-separated text file.
    output_name: complete path with filename (without the extension) for the file to
        write the results to.
    plot_output: complete path pluts file naming convention (without the .html)
        for the file to write the plots to.
    high_res: the accuracy of the conversion from radius or altitude to pressure
        level. Ignore if no conversion is needed for the variable(s) selected. Default is 20.
    plot_coord: one of 'GDZ', 'GEO', 'GSM', 'GSE', 'SM', 'GEI', 'MAG'
        integers also allowed with 'GDZ'=0 and so on. Indicates the coordinate
        system the plot will be generated in. Only plots in cartesian coordinates
        systems are supported, so 'SPH' and 'RLL' are not accepted. Default is 
        'GEO'.        
    verbose: Set to true to be overwhelmed with information.
    
    Returns a dictionary with keys: 'utc_time', 'c1', 'c2', 'c3', and 'net_idx'
        - utc_time is an array in UTC seconds since 1970-01-01 of the given
            timestamps with any occuring outside of the model data removed.
        - 'c1', 'c2', and 'c3' are arrays of the given coordinate values for each 
            surviving timestamp
        - 'net_idx' is the original index value of the surviving timestamps. 
            This is kept for easier comparison with the original dataset.
        - additional keys are included for each variable and label an array of the
            values of the indicated variable for each time+spatial coordinate given
        - The units of each array in the returned dictionary are printed to 
            the screen.
            
    See kamodo_ccmc.flythrough.utils.ConvertCoord for info on the coordinate systems.
    '''

    #read in trajectory from file into dictionary, including metadata
    if file_type=='cdf4':
        traj_data = WO.SFcdf_reader(traj_file)
    elif file_type=='csv':
        traj_data=WO.SFcsv_reader(traj_file)
    elif file_type=='txt':
        traj_data=WO.SFascii_reader(traj_file)
    else:
        raise AttributeError('File type not recognized. Must be one of'+\
                            ' cdf4, csv, or txt.')
        
    #figure out key for time data
    for key in traj_data:
        if 'time' in key: 
            time_key = key
            break
    
    #call satellite flythrough code
    results = ModelFlythrough(model, file_dir, variable_list, traj_data[time_key]['data'], 
                              traj_data['c1']['data'], traj_data['c2']['data'], 
                              traj_data['c3']['data'], traj_data['metadata']['coord_type'], 
                              traj_data['metadata']['coord_grid'], 
                              output_type=output_type, output_name=output_name, 
                              plot_output=plot_output, plot_coord=plot_coord, 
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
            output_type = argv[15]
            output_name = argv[16]
            plot_output = argv[17]
            plot_coord = argv[18]
            
            #check input
            print(f'\nstart_time: {start_time}, \nstop_time: {stop_time},'+\
                  f'\nmodel: {model}, \nfile_dir: {file_dir},'+\
                  f'\nvariable_list: {variable_list}, \nmax_lat: {max_lat}, '+\
                  f'\nmin_lat: {min_lat}, \nlon_perorbit: {lon_perorbit}, '+\
                  f'\nmax_height: {max_height}, \nmin_height {min_height},'+\
                  f'\np: {p}, n: {n}, \noutput_type: {output_type}, '+\
                  f'\noutput_name: {output_name}, \nplot_output: {plot_output}\n'+\
                  f'\nplot_coord: {plot_coord}')

            results = FakeFlight(start_time, stop_time, model, file_dir, variable_list, 
                                 max_lat=max_lat, min_lat=min_lat, lon_perorbit=lon_perorbit, 
                                 max_height=max_height, min_height=min_height,
                                 p=p, n=n, high_res=high_res, verbose=False, 
                                 output_type=output_type, output_name=output_name,
                                 plot_output=plot_output, plot_coord=plot_coord)

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
            output_type = argv[9]
            output_name = argv[10]
            plot_output = argv[11]
            plot_coord = argv[12]
            high_res = float(argv[13])            
            
            #check input
            print(f'\ndataset: {dataset}, \nstart: {start}, \nstop: {stop},'+\
                  f'\nmodel: {model}, \nfile_dir: {file_dir},'+\
                  f'\nvariable_list: {variable_list}, \ncoord_type: {coord_type},'+\
                  f'\noutput_type: {output_type}, \noutput_name: {output_name}, '+\
                  f'\nplot_output: {plot_output}, \nplot_coord: {plot_coord}, '+\
                  f'\nhigh_res: {high_res}\n')

            results = RealFlight(dataset, start, stop, model, file_dir, variable_list, 
                                 coord_type=coord_type, output_type=output_type, 
                                 output_name=output_name, plot_output=plot_output, 
                                 plot_coord=plot_coord, high_res=high_res)

        elif argv[1]=='MyFlight':  #gather variables and call MyFlight
            traj_file = argv[2]
            file_type = argv[3]            
            if len(argv[4])==1:
                model = int(argv[4])
            else:
                model = argv[4]
            file_dir = argv[5]
            temp_str = argv[6][1:-1].replace("'","").replace(' ','').replace('"','')
            variable_list = temp_str.split(',')   #['rho','N_n']              
            output_type = argv[7]
            output_name= argv[8]
            plot_output = argv[9]
            plot_coord = argv[10]
            high_res = float(argv[11])
            
            #check inputs
            print(f'\ntraj_file: {traj_file}, \nfile_type: {file_type},'+\
                  f'\nmodel: {model}, \nfile_dir: {file_dir},'+\
                  f'\nvariable_list: {variable_list}, \noutput_type: {output_type},'+\
                  f'\noutput_name: {output_name}, \nplot_output: {plot_output}, '+\
                  f'\nplot_coord: {plot_coord}, \nhigh_res: {high_res}\n')
            
            results = MyFlight(traj_file, file_type,  
                               model, file_dir, variable_list, output_type=output_type, 
                               output_name=output_name, plot_output=plot_output, 
                               plot_coord=plot_coord, high_res=high_res)
        else:
            print('Call signature not recognized.')
    else:
        print('\nPossible call types (first argument): FakeFlight, RealFlight, MyFlight')
        print('Use the call type as the first input to get call syntax.\n')
        U.MW.Choose_Model('')  #asking for a list of possible models
        
        
