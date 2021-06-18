# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 17:09:14 2021

@author: rringuet

Most needed functions to support the SatelliteFlythrough and SingleSatelliteFlythrough
softwares. Need to add names of additional pressure level variations to list below
as encountered. The corresponding height function to be inverted by CalcIlev
will need to be labeled H_ilev for ilev, H_lev for the lev variation, etc.

When Zach's new method of height extrapolation is ready, it is to be called in 
CalcIlev (around line 300) where the note is.
"""
import glob, os
import numpy as np
import time as ti
from datetime import datetime, timedelta, timezone
import kamodo.satelliteflythrough.model_wrapper as MW
import kamodo.satelliteflythrough.FlythroughPlots as FPlot
from pprint import pprint


ilev_list = ['ilev','lev','imlev']

@np.vectorize
def ts_to_hrs(time_val, filedate):
    '''Convert array of timestamps to hours since midnight of filedate string'''
    
    file_datetime = datetime.strptime(filedate+' 00:00:00', '%Y-%m-%d %H:%M:%S')
    return (datetime.utcfromtimestamp(time_val)-file_datetime).total_seconds()/3600.

@np.vectorize
def hrs_to_ts(time_val, filedate):
    '''Convert array of hours since midnight of filedate string to timestamps'''
    
    file_datetime = datetime.strptime(filedate+' 00:00:00', '%Y-%m-%d %H:%M:%S').replace(tzinfo=timezone.utc)    
    return datetime.timestamp(file_datetime+timedelta(hours=time_val))

def check_plot_dir(plot_file):
    '''If plot_dir does not exist, create it'''
    
    if '-' in plot_file.split('/')[-1].split('\\')[-1]: 
        plot_dir = plot_file[0:-10]  #cut off date 
    else:
        plot_dir = plot_file
    if not os.path.isdir(plot_dir): os.mkdir(plot_dir)
    return

def day_files(file_pattern, reader):
    '''Retrieve file times. Convert files if necessary.'''

    #file_pattern could be a list if more than one pattern in file_dir exists   
    if not isinstance(file_pattern, str):  #if a list/array of strings (GITM/similar)
        files, times, ts_list = file_pattern, {}, []  #run reader with file_prefixes given
    else: 
        files, times, ts_list = glob.glob(file_pattern), {}, []
    for f in files:
        k = reader(f, variables_requested=[], filetimes=True, printfiles=False)
        if hasattr(k, 'conversion_test'): 
            continue  #if reader requires multiple files and only one found, skip file_pattern
        file_date = k.datetimes[0][0:10]
        times[file_date] = [f,k.timerange['min'],k.timerange['max'],
                            k.filetimes[0], k.filetimes[1]]
        ts_list.extend([k.filetimes[0], k.filetimes[1]])
    return times, ts_list

def save_times(file_patterns, sat_time, reader, dt=450., verbose=False):
    '''Adjust times between files to filetime within half of dt in seconds (7.5min by default).
    Works for models with one day of data per file.'''

    times, ts_list = day_files(file_patterns, reader)  
    #look for sat_times not in files
    sat_time_mask = np.zeros(len(sat_time), bool)  #mask of indices to not change
    idx_arr = np.linspace(0,len(sat_time)-1,len(sat_time), dtype=int)  #list of all indices
    for file_date in times.keys():  #filter out times in files
        idx = np.where((sat_time>=times[file_date][3]) & (sat_time<=times[file_date][4]))[0]
        sat_time_mask[idx] = True  #True to not change times in files
    change_idx = np.delete(idx_arr, sat_time_mask)  #indices of times not in files
        
    #correct times not in files to best time within dt/2, else leave as is
    new_times = np.array([ts_list[abs(np.array(ts_list)-s).argmin()] if \
                          (abs(np.array(ts_list)-s).min() < dt) \
                 else -1. for s in sat_time[change_idx]])
    new_idx = np.where(new_times > -1.)[0]
    print(f'\nAdjusted {len(new_idx)} times within {dt/60.} minutes of a file to be in the nearest file.')
    if verbose:
        for i in new_idx: 
            print(f'Given: {sat_time[change_idx[i]]:.3f}, Nearest: {new_times[i]:.3f}, '+\
                  f'Time diff (minutes): {(sat_time[change_idx[i]]-new_times[i])/60.:.3f}')
    sat_time[change_idx[new_idx]] = new_times[new_idx]
    #for i in new_idx: print(f'{sat_time[change_idx[i]]:.3f}')  #shows that replacement actually happened
    bad_idx = change_idx[np.where(new_times == -1.)[0]]
    net_idx = np.delete(idx_arr, bad_idx)
    
    #print errors for any remaining 'bad' times
    if len(bad_idx)>0:
        print(f'{len(bad_idx)} data points are farther than {dt/60.:.2f} minutes from model times and are excluded.')
        if verbose: print(sat_time[bad_idx])  #print 'bad' sat_times
        print('\nCheck that data fits in file time ranges:')
        print(f'Data time range (UTC): {min(sat_time)} {max(sat_time)}')
        print('Filename, Min DateTime, Max DateTime, Min Time (UTC), Max Time (UTC)')
        for file_date in times.keys(): 
            print (times[file_date][0].split('/')[-1].split('\\')[-1], times[file_date][3:5])
            
    #distribute indices of sat_time times for each file
    for file_date in times.keys():  #filter out times in files
        idx = np.where((sat_time>=times[file_date][3]) & (sat_time<=times[file_date][4]))[0]
        times[file_date].append(idx)
        
    return sat_time, times, net_idx

def make_net_plots(var, results_dict, results_units, varlist_4d, varlist_3d, file_dir,
              plot_close, plot_sampling):
    '''Make net plots for a flythrough of one variable'''
    
    check_plot_dir(file_dir+'Plots/')  #create proper directory if DNE
    if var in varlist_4d:
        FPlot.Plot4D(var, results_dict['sat_time'], results_dict['sat_lat'], 
                     results_dict['sat_lon'], results_dict['sat_height'], 
                     results_dict[var], results_units[var], 
                     sat_zunits='km', plot_file=file_dir+'Plots/', 
                     plot_close=plot_close, plot_sampling=plot_sampling)
    elif var in varlist_3d:
        FPlot.Plot3D(var, results_dict['sat_time'], results_dict['sat_lat'], 
                     results_dict['sat_lon'], results_dict[var], results_units[var], 
                     plot_file=file_dir+'Plots/', plot_close=plot_close, 
                     plot_sampling=plot_sampling)    
        
def make_daily_plots(var, idx, sat_time, sat_lat, sat_lon, results, results_units, 
                     plot_file, plot_sampling, sat_height=0, sat_ilev=0):
    '''Make daily plots for a flythrough of one variable'''
    
    #Choose correct plotting routing  
    if idx=='0':  #independent of height
        FPlot.Plot3D(var, sat_time, sat_lat, sat_lon, results, 
                     results_units, plot_file, plot_sampling=plot_sampling)
    elif idx=='1': #ALL functions that require height have input height in km
        FPlot.Plot4D(var, sat_time, sat_lat, sat_lon, sat_height, 
                     results, results_units, plot_file, 
                     sat_zunits='km', plot_sampling=plot_sampling) 
    elif idx=='2':  #all others require ilev
        FPlot.Plot4D_ilev(var, sat_time, sat_lat, sat_lon, sat_height, 
                          sat_ilev, results, results_units, plot_file,
                          sat_zunits='km', sat_ilevunits='', 
                          sat_ilevname='Pressure\nLevel', plot_sampling=plot_sampling)

def call_FlyAway(model, file_dir, times, variable_list, sat_time, 
                 sat_height, sat_lat, sat_lon, plot_sampling, high_res, daily_plots, 
                 verbose=False):
    '''interpolate requested data for each day'''
    
    if verbose: print('\nInterpolating data for each file.')
    if daily_plots:  #interpolate data using idx list from before and make daily plots
        check_plot_dir(file_dir+'Plots/')  #create proper directory if DNE
        list_results = [Model_FlyAway(model, times[file_date][0], variable_list, 
                                      ts_to_hrs(sat_time[times[file_date][5]], file_date),
                                      sat_height[times[file_date][5]], sat_lat[times[file_date][5]], 
                                      sat_lon[times[file_date][5]], high_res, plot_sampling=plot_sampling,
                                      plot_file=file_dir+'Plots/'+file_date+'-',
                                      verbose=verbose) \
                        for file_date in times.keys() if len(sat_time[times[file_date][5]])>0]
    else:  #interpolate data using idx list from before without making daily plots
        list_results = [Model_FlyAway(model, times[file_date][0], variable_list, 
                                      ts_to_hrs(sat_time[times[file_date][5]], file_date),
                                      sat_height[times[file_date][5]], sat_lat[times[file_date][5]], 
                                      sat_lon[times[file_date][5]], high_res, plot_file='',
                                      verbose=verbose) \
                        for file_date in times.keys() if len(sat_time[times[file_date][5]])>0]
    return list_results

def collect_results(sat_time, sat_height, sat_lat, sat_lon, net_idx, varnames, 
                    varlist_3d, variable_list, list_results, plots, file_dir, 
                    plot_close, plot_sampling):
    '''Collect results from different days into one dictionary'''

    #combine idx lists for plotting, collect filtered trajectory
    results_dict = {'sat_time': sat_time[net_idx], 'sat_height': sat_height[net_idx],
                    'sat_lat': sat_lat[net_idx], 'sat_lon': sat_lon[net_idx],
                    'net_idx': net_idx}  #index for comparison with other data from real satellite
    
    #determine units for all variables
    results_units = {value[0]:value[-1] for key, value in varnames.items() \
                     if value[0] in variable_list}
    results_units['sat_height'], results_units['net_idx'] = 'km', ''
    results_units['sat_lat'], results_units['sat_lon'] = 'deg', 'deg'
    results_units['sat_time'] = 's'
        
    #collect interpolated data into the same dictionary
    varlist_4d = [value[0] for key, value in varnames.items()]
    for item in varlist_3d: varlist_4d.remove(item)  #remove 3D items from 4D list
    for var in variable_list:  #sort and combine arrays for the same variable
        results_dict[var] = np.concatenate(tuple([results[var] for results in list_results]))
        if plots:  #make a net plot per variable if desired, based function dependencies (height or no)
            make_net_plots(var, results_dict, results_units, varlist_4d, varlist_3d, file_dir,
              plot_close, plot_sampling)
            
    return results_dict, results_units

''' Using CTIPe object for these tests:
#speed test for CalcIlev section, for two step process, one satellite position at a time
#outputs n_stepd, bounds of height as delta_h, bounds of ilev for delta_ilev, ilev, 
#    and processing time
t, height, lat, lon = sat_track[0]
sample_ilev = np.linspace(1,15,15,dtype=float)   #global variable
for i in [10,50,100,250,500]:
    start = ti.time()
    rough_height = np.array([ctipe.H([t, ilev, lat, lon])[0] for ilev in sample_ilev])
    ilev_range = np.sort(sample_ilev[np.argsort(abs(height-rough_height))[0:2]])
    test_ilev = np.linspace(ilev_range[0],ilev_range[1],i,dtype=float)
    finer_height = np.array([ctipe.H([t, ilev, lat, lon])[0] for ilev in test_ilev])
    ilev_idx = np.argmin(abs(height-finer_height))
    print(i, finer_height[ilev_idx+1]-finer_height[ilev_idx-1], 
          test_ilev[ilev_idx+1]-test_ilev[ilev_idx-1], test_ilev[ilev_idx],ti.time()-start)

OUTPUT: (i, delta_h, delta_ilev, ilev, calc_time)
10 10061.973765432136 0.2222222222222232 12.777777777777779 0.01804208755493164
50 1848.1176303854445 0.040816326530611846 12.775510204081632 0.04254770278930664
100 914.7248877665261 0.020202020202019 12.787878787878787 0.07263469696044922
250 363.68579875055 0.008032128514056325 12.783132530120483 0.20418334007263184
500 181.47848474729108 0.004008016032063466 12.783567134268537 0.3631880283355713
-> Too slow to execute repeatedly
   

#speed test with one step process, same outputs
t, height, lat, lon = sat_track[0]
for i in [10,50,100,500,1000,5000,10000]:
    test_ilev = np.linspace(1,15,n,dtype=float)  #global variable each time
    start = ti.time()
    test_track = np.array([[t, ilev, lat, lon] for ilev in test_ilev])
    finer_height = ctipe.H(test_track)
    ilev_idx = np.argmin(abs(height-finer_height))
    print(i, finer_height[ilev_idx+1]-finer_height[ilev_idx-1], 
          test_ilev[ilev_idx+1]-test_ilev[ilev_idx-1], test_ilev[ilev_idx],ti.time()-start)
    
OUTPUT: (i, delta_h, delta_ilev, ilev, calc_time)
10 152055.1091820988 3.1111111111111107 13.444444444444445 0.0
50 25873.64682539692 0.571428571428573 12.714285714285714 0.0025177001953125
100 12806.148428731773 0.28282828282828376 12.737373737373737 0.0010514259338378906
500 2540.6987864619005 0.056112224448899184 12.783567134268537 0.0
1000 1269.0777722166968 0.02802802802802873 12.785785785785786 0.006573677062988281
5000 253.6124613811844 0.005601120224044465 12.784756951390278 0.002038717269897461
10000 126.79354879941093 0.0028002800280031437 12.783578357835783 0.010225057601928711  *chosen method
-> Not accurate to within 200 km until i=10000, but better speed than 1st option


#faster 2-step process?
t, height, lat, lon = sat_track[0]
sample_ilev = np.linspace(1,15,15,dtype=float)   #global variable
for i in [10,50,100,250,500,1000,5000,10000]:
    start = ti.time()
    rough_track = np.array([[t, ilev, lat, lon] for ilev in sample_ilev])
    rough_height = ctipe.H(rough_track)
    ilev_range = np.sort(sample_ilev[np.argsort(abs(height-rough_height))[0:2]])
    test_ilev = np.linspace(ilev_range[0],ilev_range[1],i,dtype=float)
    finer_track = np.array([[t, ilev, lat, lon] for ilev in test_ilev])
    finer_height = ctipe.H(finer_track)
    ilev_idx = np.argmin(abs(height-finer_height))
    print(i, finer_height[ilev_idx+1]-finer_height[ilev_idx-1], 
          test_ilev[ilev_idx+1]-test_ilev[ilev_idx-1], test_ilev[ilev_idx],ti.time()-start)
    
OUTPUT: (i, delta_h, delta_ilev, ilev, calc_time)
10 10061.973765432136 0.2222222222222232 12.777777777777779 0.0
50 1848.1176303854445 0.040816326530611846 12.775510204081632 0.0
100 914.7248877665261 0.020202020202019 12.787878787878787 0.0
250 363.68579875055 0.008032128514056325 12.783132530120483 0.015614032745361328
500 181.47848474729108 0.004008016032063466 12.783567134268537 0.0   
1000 90.64841230114689 0.0020020020020012907 12.783783783783784 0.0
5000 18.115175812970847 0.0004000800160035567 12.783956791358271 0.0
10000 9.056682057096623 0.00020002000199959014 12.783978397839784 0.015622377395629883

-> Accurate to within 200 km for i=500, and too fast to time it. Best option out of first 3.
Freezes when including i values above 500. Not sure why. Gives answer if I hit enter a few times.

-> Comparing method 2 (i=10000) and 3 (i=500) in execution showed a large time difference.
Method 2 took about 150 seconds while method 3 took about 16 seconds. Choosing method 3.

#Method 2
sample_ilev = np.linspace(1,15,10000, dtype=float)
def CalcIlev2(H, t, height, lat, lon):
    
    finer_height = H(np.array([[t, ilev, lat, lon] for ilev in sample_ilev]))
    ilev_idx = np.argmin(abs(height-finer_height))
    return sample_ilev[ilev_idx] 
'''  
def CalcIlev(H, Hunit, t, height, lat, lon, min_ilev, max_ilev, high_res, verbose=True):
    '''Approximate ilev by inverting the chosen height function for one sat point.
    high_res default is 20 meters. (high_res), input height is in meters.'''
    
    #input height is in meters, time is in hours since midnight.
    #convert for common function output units. Add as more are discovered.
    if Hunit=='m':
        Hconv=1
    elif Hunit=='cm':
        Hconv=1./100.  #convert from cm to m
    elif Hunit=='km':
        Hconv=1000.   #convert from km to m
    
    #get height function output for the ilev range allowed in model
    sample_ilev = np.linspace(min_ilev,max_ilev,75,dtype=float)  #typical range is 15
    rough_height = H(np.array([[t, ilev, lat, lon] for ilev in sample_ilev]))*Hconv
    max_height, model_ilev_max = rough_height.max(), max_ilev  #save for output
    
    #allow extrapolation for heights ABOVE height function range for this time/location
    #when/if Zach writes a more physical way to extrapolate, can call from here.
    while rough_height.max()<height:
        max_ilev+=1
        sample_ilev = np.linspace(max_ilev-1,max_ilev,5,dtype=float)
        rough_height = H(np.array([[t, ilev, lat, lon] for ilev in sample_ilev]))*Hconv
        
    #allow extrapolation for heights BELOW height function range for this time/location
    #when/if Zach writes a more physical way to extrapolate, can call from here.
    while rough_height.min()>height:
        min_ilev-=1
        sample_ilev = np.linspace(min_ilev,min_ilev+1,5,dtype=float)
        rough_height = H(np.array([[t, ilev, lat, lon] for ilev in sample_ilev]))*Hconv    
    
    #continue with numerical inversion
    ilev_range = np.sort(sample_ilev[np.argsort(abs(height-rough_height))[0:2]])  
    test_ilev = np.linspace(ilev_range[0],ilev_range[1],100,dtype=float) 
    finer_height = H(np.array([[t, ilev, lat, lon] for ilev in test_ilev]))*Hconv
    
    #loop through process until requested resolution is acheived
    loop_count=0
    while abs(height-finer_height).min()>high_res: 
        ilev_range = np.sort(test_ilev[np.argsort(abs(height-finer_height))[0:2]])
        test_ilev = np.linspace(ilev_range[0],ilev_range[1],100,dtype=float)
        finer_height = H(np.array([[t, ilev, lat, lon] for ilev in test_ilev]))*Hconv
        #limit iterations to prevent infinite loops
        loop_count+=1
        if loop_count>10: break
    
    #output info for inspection
    if verbose: 
        print(f'\nRequested height: {height:.3f} m')
        print(f'Maximum allowed height: {max_height:.3f} m')
        print(f'Maximum presure level allowed for model: {model_ilev_max}')
        if max_ilev>model_ilev_max:
            print(f'Adjusted maximum pressure level: {max_ilev}')
        print(f'Number of loops required to acheive resolution: {loop_count}')
        print(f'Calculated equivalent pressure level: {test_ilev[np.argmin(abs(height-finer_height))]}')
        print(f'Height resolution achieved: {abs(height-finer_height).min():.5f} m\n')
    return [test_ilev[np.argmin(abs(height-finer_height))], abs(height-finer_height).min()]

def call_CalcIlev(ilev_string, kamodo_object, sat_track, sat_track0, model_sat_time, 
                  sat_lat, sat_lon, high_res, call_type):
    '''ilev agnostic method to gather parameters and call CalcIlev with less code.'''
    
    #retrieve parameters for CalcIlev call
    min_ilev = getattr(kamodo_object, '_'+ilev_string).min()
    max_ilev = getattr(kamodo_object, '_'+ilev_string).max()
    Hfunc = getattr(kamodo_object, 'H_'+ilev_string) 
    Hunit = kamodo_object.variables['H_'+ilev_string]['units']
    
    #for each call version, call CalcIlev and build satellite track with pressure level
    if call_type=='single':
        sat_ilev, height_res = CalcIlev(Hfunc, Hunit, *sat_track0, min_ilev, max_ilev, 
                            high_res, verbose=True)  
        sat_track[ilev_string]=[model_sat_time,sat_ilev,sat_lat,sat_lon]
        #User feedback for single call programmed into CalcIlev function
    elif call_type=='multiple':
        sat_ilev, height_res = np.array([CalcIlev(Hfunc, Hunit, *sat_position, min_ilev, max_ilev, 
                           high_res, verbose=False) for sat_position in sat_track0]).T
        sat_track[ilev_string]=[[t, ilev, lat, lon] for t, ilev, lat, lon in zip(
                                model_sat_time,sat_ilev,sat_lat,sat_lon)]
        #Give user feedback about range of height resolution for height to ilev conversion
        print(f'\nBest height resolution achieved: {height_res.min():.5f} m')
        print(f'Worst height resolution achieved: {height_res.max():.5f} m\n')
    return sat_track, sat_ilev

def sat_tracks(variable_list, kamodo_object, sat_time, sat_height, sat_lat, sat_lon,
               high_res, verbose=False):
    '''Calculate satellite tracks for interpolation'''

    #determine which kind of satellite track(s) are needed
    var_test = {'none':[], 'height':[]}
    for ilev_string in ilev_list:  #add ilev names from list at top of file
        var_test[ilev_string] = []
    for var in variable_list:  #determine which variables require ilev(2), height(1), neither(0)
        var_list0 = list(kamodo_object.variables[var]['xvec'].keys())[1]
        if var_list0 in var_test.keys():
            var_test[var_list0].append(var)
        else:
            var_test['none'].append(var)
    
    #Create satellite tracks with appropriate inputs
    sat_track = {}  #initialize list of satellite tracks
    if len(var_test['none'])>0:
        if verbose: print('Building height-independent satellite track.')
        sat_track['none']=[[t, sat_lat, sat_lon] for t, sat_lat, sat_lon in zip(
            sat_time,sat_lat,sat_lon)]
        sat_ilev=0
    if len(var_test['height'])>0:  #if function requires height (in km)
        if verbose: print('Building height-dependent satellite track (km).')
        sat_track['height']=[[t, h, lat, lon] for t, h, lat, lon in zip(
            sat_time,sat_height,sat_lat,sat_lon)]
        sat_ilev=0
    if sum([len(var_test[i]) for i in ilev_list]):  
        #if ilev, lev, or imlev is required for at least one variable
        if verbose: print('Converting height to pressure level and building ilev-dependent satellite track.')
        start = ti.time()  #Input h in meters
        sat_track0 = [[t, h*1000., lat, lon] for t, h, lat, lon in zip(
            sat_time,sat_height,sat_lat,sat_lon)]
        for ilev_string in ilev_list:
            if len(var_test[ilev_string])>0:
                #add new track type to dictionary
                sat_track, sat_ilev = call_CalcIlev(ilev_string, kamodo_object, sat_track, 
                                          sat_track0, sat_time, sat_lat, 
                                          sat_lon, high_res, 'multiple')    
                if verbose: print(f'{ilev_string} track added')              
        if verbose: print(f'Conversion took {ti.time()-start} s for {len(sat_time)} positions.')
    
    return var_test, sat_track, sat_ilev

def interpolate_data(variable_list, var_test, kamodo_object, sat_track, plot_file,
                     sat_time, filedate, sat_lat, sat_lon, plot_sampling,
                     sat_height, sat_ilev):
    '''retrieve interpolator and interpolate data for each variable.'''
    
    results = {}
    for var in variable_list:
        #choose correct satellite track and interpolate data
        idx = [key for key, value in var_test.items() if var in value][0]
        intp = kamodo_object[var]   #gridded interpolator way too slow (the one with _ijk in name)
        results[var] = intp(sat_track[idx])
        results_units = intp.meta['units']  #give entire track of correct type 
    
        #make correct set of prints depending on the dimensionality of the function
        if plot_file != '':       
            make_daily_plots(var, idx, hrs_to_ts(sat_time, filedate), 
                               sat_lat, sat_lon, results[var], results_units, 
                               plot_file, plot_sampling, sat_height=sat_height, 
                               sat_ilev=sat_ilev)
    return results


def Model_FlyAway(model, filename, variable_list, sat_time, sat_height, sat_lat, sat_lon, 
                  high_res, plot_sampling=4000, plot_file='', verbose=False):
    '''Functions that are generic to each wrapper in daily processing'''

    #Check that sat data is all the same length, will error if not
    sat_data_check(sat_time, sat_height, sat_lat, sat_lon)
    
    #create kamodo object, initialize some variables
    reader = MW.Model_Reader(model)
    kamodo_object = reader(filename, variables_requested=variable_list, printfiles=False, 
                  gridded_int=False)
    for item in ilev_list:  #check for each H variation
        if 'H_'+item in variable_list: variable_list.remove('H_'+item)  #H added by default
    
    #create satellite tracks of types needed based on vertical dependencies
    filedate = kamodo_object.timerange['max'][0:10]
    var_test, sat_track, sat_ilev = sat_tracks(variable_list, kamodo_object, 
                                               sat_time, sat_height, sat_lat, sat_lon, 
                                               high_res, verbose=verbose)

    #retrieve interpolator and interpolate data for each variable. 
    results = interpolate_data(variable_list, var_test, kamodo_object, sat_track, 
                               plot_file, sat_time, filedate, sat_lat, sat_lon, 
                               plot_sampling, sat_height, sat_ilev)
    
    if verbose: print(f'Done for {filename}\n')
    return results

def sat_data_check(sat_time, sat_height, sat_lat, sat_lon):
    '''Check that satellite trajectory arrays are all the same length.'''
    
    if max(np.diff(np.array([len(sat_time),len(sat_height),len(sat_lat),len(sat_lon)])))>0:
        raise AttributeError (f'Satellite time, height, latitude, and longitude\
                              arrays or lists must all be the same length.\
                              Current array lengths are {len(sat_time)}, {len(sat_height)},\
                               {len(sat_lat)}, and {len(sat_lon)}')
    return

def generic_SatelliteFlythrough(file_dir, pattern, sat_time, kamodo_object, 
                                model, variable_list, sat_height, sat_lat, 
                                sat_lon, daily_plots, plot_sampling, high_res, varnames, 
                                varlist_3d, plots, plot_close, dt=450., verbose=False):
    '''Functions that are generic to each wrapper in net processing'''
    
    sat_time, times, net_idx = save_times(pattern, sat_time, kamodo_object, 
                                            dt=dt, verbose=verbose)
        
    #interpolate requested data for each day. FlyAway is specific to each wrapper
    list_results = call_FlyAway(model, file_dir, times, variable_list, sat_time, 
                                  sat_height, sat_lat, sat_lon, plot_sampling, high_res, 
                                  daily_plots, verbose=verbose)
    
    #collect results from different days into one dictionary
    results_dict, results_units = collect_results(sat_time, sat_height, sat_lat, sat_lon, 
                                     net_idx, varnames, varlist_3d, 
                                     variable_list, list_results, plots, file_dir, 
                                     plot_close, plot_sampling)
    print('Results_units: ')
    pprint(results_units)
    return results_dict, results_units 

def Model_SatelliteFlythrough(model, file_dir, variable_list, sat_time, sat_height, sat_lat, 
                              sat_lon, high_res, dt=450., plots=False, daily_plots=False, plot_close=True, 
                              plot_sampling=4000, verbose=False):
    '''
    Execute flythrough for model data. Returns results_dict, results_units.
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
    
    #reader prefers wrapped filename, even if it does not exist. will create if no wrapped data found.
    file_patterns = MW.FileSearch(model, file_dir)
    reader = MW.Model_Reader(model)
    Var = MW.Model_Variables(model, return_dict=True)
    varlist_3d = MW.Var_3D(model)
    results_dict, results_units = generic_SatelliteFlythrough(file_dir, file_patterns, 
                                sat_time, reader, model, variable_list, sat_height, 
                                sat_lat, sat_lon, daily_plots, plot_sampling, high_res, 
                                Var, varlist_3d, plots, plot_close, 
                                dt=dt, verbose=verbose)
    return results_dict, results_units

#----------Code below is for Single Satellite Flythrough software--------------------------------
def generic_single_prerun(file_patterns, reader, variable_list):
    '''Converts data files as needed, returns vertical dependency of variables given.'''
    
    times, ts_list = day_files(file_patterns, reader)  #reader creates wrapped files if DNE
            
    #create ctipe object, return vertical dependencies for variables requested
    files = [times[key][0] for key in times.keys()]
    var_list = variable_list.copy()  #keep original copy of variable_list
    kamodo_object = reader(files[0], variables_requested=variable_list, printfiles=False)
    var_test = []
    for var in var_list:  #determine which variables require ilev(2), height(1), neither(0)
        input_var_z = list(kamodo_object.variables[var]['xvec'])[1]
        if input_var_z in ilev_list or input_var_z=='height': 
            var_test.append(input_var_z)
        else: 
            var_test.append('none')
    dt = kamodo_object.dt  #determine dt from model grid
    return var_test, dt

def find_singlefiletime(file_patterns, sat_time, reader, dt=450.):
    '''Find file containing given single time. Adjust if within dt seconds.'''
    
    times, ts_list = day_files(file_patterns, reader)  #Retrieve file times
    
    #Check if time is not in a file. Correct if within dt. Break if not.
    time_check = False
    for file_date in times.keys():  #filter out times in files
        if (sat_time>=times[file_date][3]) & (sat_time<=times[file_date][4]):
            time_check = True
    if not time_check:  #correct if within dt
        new_time = [ts_list[abs(np.array(ts_list)-sat_time).argmin()] if \
                      (abs(np.array(ts_list)-sat_time).min() < dt) else -1.][0]
        if new_time==-1: 
            print(f'Time not within {dt} seconds of possible files.')
            print(f'Given time is farther than {dt/60.:.2} minutes from model times.')
            print('\nCheck that time fits in file time ranges:')
            print(f'Data time (UTC): {sat_time}')
            print('Filename, Min DateTime, Max DateTime, Min Time (UTC), Max Time (UTC)')
            for file_date in times.keys(): 
                print (times[file_date][0].split('/')[-1].split('\\')[-1], times[file_date][3:5])
            raise AttributeError('No file found with the time given.')
        else: 
            print('Time given not in possible files:', sat_time)
            print(f'Nearest time is {new_time} seconds UTC')
            print(f'This is a difference of {abs(sat_time-new_time):.3f} seconds.')
            print('Proceeding with shifted time.')
            sat_time = new_time
    
    #choose file time is in
    for file_date in times.keys():
        if ((sat_time>=times[file_date][3]) & (sat_time<=times[file_date][4])):
            filename = times[file_date][0]
            break

    return filename, sat_time    

def sat_tracks_single(variable_list, z_dependence, model_sat_time, sat_height, 
                      sat_lat, sat_lon, kamodo_object, high_res):
    '''Calculate satellite tracks for interpolation'''
    
    #Create satellite tracks with appropriate inputs
    sat_track = {}  #initialize list of satellite tracks
    if z_dependence.count('none')>0:
        #if verbose: print('Building height-independent satellite track.')
        sat_track['none']=[model_sat_time,sat_lat,sat_lon]
    if z_dependence.count('height')>0:  #if function requires height (in km)
        #if verbose: print('Building height-dependent satellite track (km).')
        sat_track['height']=[model_sat_time,sat_height,sat_lat,sat_lon]
    if sum([z_dependence.count(i) for i in ilev_list])>0:  
        #if ilev is required for at least one variable
        sat_track0 = [model_sat_time,sat_height*1000.,sat_lat,sat_lon] #set up track
        for ilev_string in ilev_list:
            if z_dependence.count(ilev_string)>0:
                #add new track type to dictionary
                sat_track, sat_ilev = call_CalcIlev(ilev_string, kamodo_object, sat_track, 
                                          sat_track0, model_sat_time, sat_lat, 
                                          sat_lon, high_res, 'single')
        #if verbose: print(f'Conversion took {ti.time()-start} s for given position.')
    #print(sat_track)
    return sat_track #sat_ilev not used here because not plotting anything

def generic_Single_FlyAway(kamodo_object, variable_list, z_dependence, dz,
                           sat_time, sat_height, sat_lat, sat_lon, high_res):
    '''fly satellite through model data, per position
    sat_time, sat_height, sat_lat, and sat_lon are all floats, not arrays
    z_dependence = ["none","height","ilev"] if variables depend on all three options
       dependence must be in same order as variable_list to match with variable names'''
    

    #assume given sat_height is in km, convert time to hrs since midnight
    str_date = datetime.strftime(kamodo_object.filedate, format='%Y-%m-%d')
    model_sat_time = ts_to_hrs(sat_time, str_date)
    sat_track = sat_tracks_single(variable_list, z_dependence, model_sat_time, 
                                  sat_height, sat_lat, sat_lon, kamodo_object, 
                                  high_res)

    #retrieve interpolator and interpolate data for each variable.
    results = {variable_list[i]: kamodo_object[variable_list[i]](sat_track[z_dependence[i]])[0] 
               for i in range(len(variable_list))}
    
    #determine vertical derivatives for each variable if requested
    if len(dz)>0:  #should be a list of booleans or 1s and 0s 
        for i in range(len(variable_list)):
            if dz[i]==1 and z_dependence[i]!='none':  #if dz requested and variable has a vertical dependence
                #generate tracks with slightly larger and smaller heights
                if z_dependence[i]=='height':
                    #stay within height (km) boundaries
                    sat_height_low = sat_height-100
                    if sat_height_low <= kamodo_object.htrange0['min']: 
                        sat_height_low = kamodo_object.htrange0['min']
                    sat_height_high = sat_height+100.
                    if sat_height_high>=kamodo_object.htrange0['max']: 
                        sat_height_high = kamodo_object.htrange0['max']
                    dz_track = [[model_sat_time,sat_height_low,sat_lat,sat_lon],
                                [model_sat_time,sat_height_high,sat_lat,sat_lon]]
                #same logic, but for pressure level.
                if z_dependence[i] in ilev_list:
                    #stay within ilev boundaries
                    sat_ilev_low = sat_track[z_dependence[i]][1]-1.
                    if sat_ilev_low<=kamodo_object.iprange0['min']: #iprange0=largest range of lev variables
                        sat_ilev_low = kamodo_object.iprange0['min']
                    sat_ilev_high = sat_track[z_dependence[i]][1]+1.
                    if sat_ilev_high >= kamodo_object.iprange0['max']: 
                        sat_ilev_high = kamodo_object.iprange0['max']
                    dz_track =  [[model_sat_time,sat_ilev_low,sat_lat,sat_lon],
                                [model_sat_time,sat_ilev_high,sat_lat,sat_lon]]
                #compute requested values for offset heights, then calculate dz
                dz_result = kamodo_object[variable_list[i]](dz_track)  #returns two values
                results[variable_list[i]+'_dz'] = dz_result[1]-dz_result[0]
    return results

def Single_Prerun(model, file_dir, variable_list):
    '''Return a list of the required height input for each variable. Create wrapped files if needed.'''
    
    #Determine possible file patterns. Create wrapped files if needed.
    file_patterns = MW.FileSearch(model, file_dir)
    reader = MW.Model_Reader(model)
    z_dependence, dt = generic_single_prerun(file_patterns, reader, variable_list) 
    print(f'Vertical dependence determined for each given variable: \n{variable_list}\n{z_dependence}')
    print(f'Time resolution is {dt*3600.:.2}s. Half this and use as model_dt, or 450 seconds will be assumed.')
    return z_dependence, dt

def Single_FlyAway(model, file_dir, variable_list, z_dependence, dz, 
                   sat_time, sat_height, sat_lat, sat_lon, dt, high_res):
    '''Fly given satellite position through model code.'''

    #find file nearest sat_time (within 450 seconds), correct sat_time if needed
    file_patterns = MW.FileSearch(model, file_dir)
    reader = MW.Model_Reader(model)
    filename, sat_time = find_singlefiletime(file_patterns, sat_time, reader, dt=dt)
    
    #get ctipe object for requested variable (+H too)
    var_list = variable_list.copy()  #variable_list is changed every time the CTIPe reader is created
    kamodo_object = reader(filename, variables_requested=variable_list, printfiles=False)
            
    #get results requested for single position given
    results_dict = generic_Single_FlyAway(kamodo_object, var_list, z_dependence, dz,
                                    sat_time, sat_height, sat_lat, sat_lon, high_res)
    return results_dict



#------ Code below here is for possible link directly into fortran------------
def test_validobject(kamodo_object, sat_time):   
    ''' Determine if a new ctipe object is needed bsed on the time given'''
    
    if isinstance(kamodo_object, list):  
        return True  #if a list, then ctipe DNE, get a new one
    elif (sat_time>=kamodo_object.filetimes[0]) and (sat_time<=kamodo_object.filetimes[1]):
        return False  #sat_time is within known time range of file, use current ctipe
    else:  #sat_time is not within known time range of file, get a new ctipe object
        return True
