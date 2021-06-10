# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 17:09:14 2021

@author: rringuet
"""
import glob, os
import numpy as np
import time as ti
from datetime import datetime, timedelta, timezone
#import kamodo.satelliteflythrough.FlythroughPlots as FPlot
from pprint import pprint

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
    '''For models with all output for one day in file, retrieve file times'''

    #file_pattern could be a list if more than one pattern in file_dir exists   
    if not isinstance(file_pattern, str):  #if a list/array of strings (GITM/similar)
        files, times, ts_list = file_pattern, {}, []  #run reader with file_prefixes given
    else: 
        files, times, ts_list = glob.glob(file_pattern), {}, []
    for f in files:
        print(f)
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
        
    #correct times not in files to best time within 7.5 minutes, else leave as is
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
        print(f'{len(bad_idx)} data points are farther than 7.5 minutes from model times and are excluded.')
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

def call_FlyAway(model_FlyAway, file_dir, times, variable_list, sat_time, 
                 sat_height, sat_lat, sat_lon, plot_sampling, daily_plots, 
                 verbose=False):
    '''interpolate requested data for each day'''
    
    if verbose: print('\nInterpolating data for each file.')
    if daily_plots:  #interpolate data using idx list from before and make daily plots
        check_plot_dir(file_dir+'Plots/')  #create proper directory if DNE
        list_results = [model_FlyAway(times[file_date][0], variable_list, 
                                      ts_to_hrs(sat_time[times[file_date][5]], file_date),
                                      sat_height[times[file_date][5]], sat_lat[times[file_date][5]], 
                                      sat_lon[times[file_date][5]], plot_sampling=plot_sampling,
                                      plot_file=file_dir+'Plots/'+file_date+'-',
                                      verbose=verbose) \
                        for file_date in times.keys() if len(sat_time[times[file_date][5]])>0]
    else:  #interpolate data using idx list from before without making daily plots
        list_results = [model_FlyAway(times[file_date][0], variable_list, 
                                      ts_to_hrs(sat_time[times[file_date][5]], file_date),
                                      sat_height[times[file_date][5]], sat_lat[times[file_date][5]], 
                                      sat_lon[times[file_date][5]], plot_file='',
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
        '''if plots:  #make a net plot per variable if desired, based function dependencies (height or no)
            make_net_plots(var, results_dict, results_units, varlist_4d, varlist_3d, file_dir,
              plot_close, plot_sampling)
        '''    
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
sample_ilev = np.linspace(1,15,75,dtype=float)   #global variable
def CalcIlev(H, t, height, lat, lon):
    '''Approximate ilev by inverting the gridded height function CTIPe.H for one sat point'''
    
    rough_height = H(np.array([[t, ilev, lat, lon] for ilev in sample_ilev]))
    ilev_range = np.sort(sample_ilev[np.argsort(abs(height-rough_height))[0:2]])
    test_ilev = np.linspace(ilev_range[0],ilev_range[1],100,dtype=float)
    finer_height = H(np.array([[t, ilev, lat, lon] for ilev in test_ilev]))
    return test_ilev[np.argmin(abs(height-finer_height))]

def sat_tracks(variable_list, variables, sat_time, sat_height, sat_lat, sat_lon,
               Hfunc=None, verbose=False):
    '''Calculate satellite tracks for interpolation'''
    
    #determine which kind of satellite track(s) are needed
    var_test = {'0':[], '1':[], '2':[]}
    for var in variable_list:  #determine which variables require ilev(2), height(1), neither(0)
        input_var_list = variables[var]['xvec']
        if 'ilev' in input_var_list.keys(): var_test['2'].append(var)
        elif 'height' in input_var_list.keys(): var_test['1'].append(var)
        else: var_test['0'].append(var)
    
    #Create satellite tracks with appropriate inputs
    sat_track = {}  #initialize list of satellite tracks
    if len(var_test['0'])>0:
        if verbose: print('Building height-independent satellite track.')
        sat_track['0']=[[t, sat_lat, sat_lon] for t, sat_lat, sat_lon in zip(sat_time,sat_lat,sat_lon)]
        sat_ilev=0
    if len(var_test['1'])>0:  #if function requires height (in km)
        if verbose: print('Building height-dependent satellite track (km).')
        sat_track['1']=[[t, h, lat, lon] for t, h, lat, lon in zip(sat_time,sat_height,
                                                                     sat_lat,sat_lon)]
        sat_ilev=0
    if len(var_test['2'])>0:  #if ilev is required for at least one variable
        if verbose: print('Converting height to ilev and building ilev-dependent satellite track.')
        start = ti.time()  #The H function outputs in m, so input h in meters
        sat_track0 = [[t, h*1000., lat, lon] for t, h, lat, lon in zip(sat_time,sat_height,
                                                                     sat_lat,sat_lon)]
        sat_ilev = np.array([CalcIlev(Hfunc, *sat_position) for 
                               sat_position in sat_track0])  #try adding a parallel idea here?
        sat_track['2']=[[t, ilev, lat, lon] for t, ilev, lat, lon in zip(sat_time,
                                                    sat_ilev,sat_lat,sat_lon)]
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
    
        '''#make correct set of prints depending on the dimensionality of the function
        if plot_file != '':       
            make_daily_plots(var, idx, hrs_to_ts(sat_time, filedate), 
                               sat_lat, sat_lon, results[var], results_units, 
                               plot_file, plot_sampling, sat_height=sat_height, 
                               sat_ilev=sat_ilev)
        '''
    return results


def Generic_FlyAway(kamodo_object, variable_list, sat_time, sat_height, sat_lat,
                    sat_lon, plot_file, plot_sampling, verbose=False):
    '''Functions that are generic to each wrapper in daily processing'''
    
    #create satellite tracks of types needed based on vertical dependencies
    if hasattr(kamodo_object, 'H'): Hfunc=kamodo_object.H
    else: Hfunc=None
    filedate = kamodo_object.timerange['max'][0:10]
    var_test, sat_track, sat_ilev = sat_tracks(variable_list, kamodo_object.variables, 
                                               sat_time, sat_height, sat_lat, sat_lon, 
                                               Hfunc=Hfunc, verbose=verbose)

    #retrieve interpolator and interpolate data for each variable. 
    results = interpolate_data(variable_list, var_test, kamodo_object, sat_track, 
                               plot_file, sat_time, filedate, sat_lat, sat_lon, 
                               plot_sampling, sat_height, sat_ilev)
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
                                FlyAway, variable_list, sat_height, sat_lat, 
                                sat_lon, daily_plots, plot_sampling, varnames, 
                                varlist_3d, plots, plot_close, dt=450., verbose=False):
    '''Functions that are generic to each wrapper in net processing'''
    
    sat_time, times, net_idx = save_times(pattern, sat_time, kamodo_object, 
                                            dt=dt, verbose=verbose)
        
    #interpolate requested data for each day. FlyAway is specific to each wrapper
    list_results = call_FlyAway(FlyAway, file_dir, times, variable_list, sat_time, 
                                  sat_height, sat_lat, sat_lon, plot_sampling, 
                                  daily_plots, verbose=verbose)
    
    #collect results from different days into one dictionary
    results_dict, results_units = collect_results(sat_time, sat_height, sat_lat, sat_lon, 
                                     net_idx, varnames, varlist_3d, 
                                     variable_list, list_results, plots, file_dir, 
                                     plot_close, plot_sampling)
    print('Results_units: ')
    pprint(results_units)
    return results_dict, results_units 