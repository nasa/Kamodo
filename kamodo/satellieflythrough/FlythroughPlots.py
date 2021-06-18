# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 17:43:34 2021
@author: rringuet

Separate out plotting routines for easier accessibility
"""
from datetime import datetime
import os
import numpy as np
import matplotlib as mpl
import matplotlib.dates as mdates
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16,'figure.autolayout': True})


def Plot4D_ilev(variable_name, sat_time, sat_lat, sat_lon, sat_z, sat_ilev, 
                result, result_units, plot_file, sat_zunits='km', sat_ilevunits='',
                sat_ilevname='Pressure\nLevel', plot_close=True, plot_sampling=4000):
    '''make a 3D heat map of result, and a group of time series plots, with ilev'''
    
    #turn off interactive commands based on backend in use
    if mpl.get_backend() in ['agg','pdf','ps','svg','pgf','cairo']:
        int_off = True
        plt.ioff()  
        #print('non-interactive mode ', mpl.get_backend())
    else:
        int_off = False
        #print('interactive mode ', mpl.get_backend())

    #Reduce number of data points for 4D plot only
    i, n = 0, len(sat_time)  #still can't handle size representation, so no s
    while n > plot_sampling: 
        i += 1
        n = len(sat_time)/i
    if i==0: i=1
    
    #make 4D plot, result=color
    fig = plt.figure(figsize=(10,8))
    fig.tight_layout()
    ax = fig.add_subplot(projection='3d')  #add size = time
    ax.xaxis.set_major_locator(MaxNLocator(5)) 
    ax.yaxis.set_major_locator(MaxNLocator(5)) 
    ax.zaxis.set_major_locator(MaxNLocator(5)) 
    p = ax.scatter(sat_lon[::i], sat_lat[::i], sat_z[::i], result[::i], 
                   c=result[::i], cmap='viridis')  #works!
    ax.set_xlabel('Longitude [deg]', labelpad=15)
    ax.set_ylabel('Latitude [deg]', labelpad=15)
    ax.set_zlabel('Height ['+sat_zunits+']', labelpad=15)
    cbar = fig.colorbar(p, pad=0.1)  #add label and units to colorbar
    cbar.set_label(variable_name+' ['+result_units+']', rotation=90, labelpad=25)
    try:  #remove pre-existing file if it exists
        os.remove(plot_file+variable_name+'_3Dheat.png')
    except:
        pass
    plt.savefig(plot_file+variable_name+'_3Dheat.png')
    if not int_off:  #only include if in interactive mode
        if plot_close: plt.close()
        
    #plot result, lon, lat, ilev as 'time' series
    fig, axs = plt.subplots(5, figsize=(10,7))
    axs[0].plot(sat_time,result, 'black')
    axs[0].set_ylabel('Interpolated \n'+variable_name+' \n['+result_units+']', labelpad=10)
    axs[0].set_xticks([])
    axs[1].plot(sat_time,sat_lon, 'black')
    axs[1].set_ylabel('Sat Lon \n[deg]', labelpad=10)
    axs[1].set_xticks([])
    axs[2].plot(sat_time,sat_lat, 'black')
    axs[2].set_ylabel('Sat Lat \n[deg]', labelpad=10)
    axs[2].set_xticks([])
    axs[3].plot(sat_time, sat_ilev, 'black')
    axs[3].set_ylabel(sat_ilevname+'\n['+sat_ilevunits+']', labelpad=10)
    axs[3].set_xticks([])
    axs[4].xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y\n%H:%M:%S'))
    axs[4].xaxis.set_major_locator(MaxNLocator(5))     
    axs[4].plot([datetime.utcfromtimestamp(t) for t in sat_time], sat_z, 'black')
    axs[4].set_ylabel('Height \n['+sat_zunits+']')
    axs[4].set_xlabel('Satellite Time [UTC]', labelpad=10)   #add datetime for x axis
    try:  #remove pre-existing file if it exists
        os.remove(plot_file+variable_name+'_1Dlines.png')
    except:
        pass
    plt.savefig(plot_file+variable_name+'_1Dlines.png')
    if not int_off:  #only include if in interactive mode
        if plot_close: plt.close()
    return

def Plot4D(variable_name, sat_time, sat_lat, sat_lon, sat_z, result, result_units, 
           plot_file, sat_zunits='km', plot_close=True, plot_sampling=4000):
    '''make a 3D heat map of result, and a group of time series plots'''

    #turn off interactive commands based on backend in use
    if mpl.get_backend() in ['agg','pdf','ps','svg','pgf','cairo']:
        int_off = True
        plt.ioff()  
        #print('non-interactive mode ', mpl.get_backend())
    else:
        int_off = False
        #print('interactive mode ', mpl.get_backend())
    
    #Reduce number of data points for 4D plot only
    i, n = 0, len(sat_time)  #still can't handle size representation, so no s
    while n > plot_sampling: 
        i += 1
        n = len(sat_time)/i
    if i==0: i=1
    #s = 2.5**((sat_time-np.min(sat_time))/(np.max(sat_time)-np.min(sat_time))*4+1) #doesn't show
    
    #make 4D plot, result=color
    fig = plt.figure(figsize=(10,8))
    fig.tight_layout()
    ax = fig.add_subplot(projection='3d')  #add size = time
    ax.xaxis.set_major_locator(MaxNLocator(5)) 
    ax.yaxis.set_major_locator(MaxNLocator(5)) 
    ax.zaxis.set_major_locator(MaxNLocator(5)) 
    p = ax.scatter(sat_lon[::i], sat_lat[::i], sat_z[::i], result[::i], 
                   c=result[::i], cmap='viridis')  #works!
    ax.set_xlabel('Longitude [deg]', labelpad=15)
    ax.set_ylabel('Latitude [deg]', labelpad=15)
    ax.set_zlabel('Height ['+sat_zunits+']', labelpad=15)
    cbar = fig.colorbar(p, pad=0.1)  #add label and units to colorbar
    cbar.set_label(variable_name+' ['+result_units+']', rotation=90, labelpad=25)
    try:  #remove pre-existing file if it exists
        os.remove(plot_file+variable_name+'_3Dheat.png')
    except:
        pass
    plt.savefig(plot_file+variable_name+'_3Dheat.png')
    if not int_off:  #only include if in interactive mode
        if plot_close: plt.close()
        
    #plot result, lon, lat, ilev as 'time' series
    fig, axs = plt.subplots(4, figsize=(10,7))
    axs[0].plot(sat_time,result, 'black')
    axs[0].set_ylabel('Interpolated \n'+variable_name+' \n['+result_units+']', labelpad=10)
    axs[0].set_xticks([])
    axs[1].plot(sat_time,sat_lon, 'black')
    axs[1].set_ylabel('Sat Lon \n[deg]', labelpad=10)
    axs[1].set_xticks([])
    axs[2].plot(sat_time,sat_lat, 'black')
    axs[2].set_ylabel('Sat Lat \n[deg]', labelpad=10)
    axs[2].set_xticks([])
    axs[3].xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y\n%H:%M:%S'))
    axs[3].xaxis.set_major_locator(MaxNLocator(5))     
    axs[3].plot([datetime.utcfromtimestamp(t) for t in sat_time], sat_z, 'black')
    axs[3].set_ylabel('Height \n['+sat_zunits+']')
    axs[3].set_xlabel('Satellite Time [UTC]', labelpad=10)   #add datetime for x axis
    try:  #remove pre-existing file if it exists
        os.remove(plot_file+variable_name+'_1Dlines.png')
    except:
        pass
    plt.savefig(plot_file+variable_name+'_1Dlines.png')
    if not int_off:  #only include if in interactive mode
        if plot_close: plt.close()
    return
    
def Plot3D(variable_name, sat_time, sat_lat, sat_lon, result, result_units, 
                 plot_file, plot_close=True, plot_sampling=4000):
    '''Make a 2D heat map of the result, and a group of time series plots'''

    #turn off interactive commands based on backend in use
    if mpl.get_backend() in ['agg','pdf','ps','svg','pgf','cairo']:
        int_off = True
        plt.ioff()  
        #print('non-interactive mode ', mpl.get_backend())
    else: 
        int_off = False
        #print('interactive mode ', mpl.get_backend())
    
    #Reduce number of data points for 3D plot only
    i, n = 0, len(sat_time)  
    while n > plot_sampling: 
        i += 1
        n = len(sat_time)/i
    if i==0: i=1
    
    #make 3D plot, result=color
    fig, ax = plt.subplots(1, figsize=(10,7))
    #fig.tight_layout()
    ax.xaxis.set_major_locator(MaxNLocator(5)) 
    ax.yaxis.set_major_locator(MaxNLocator(5)) 
    p = ax.scatter(sat_lon[::i], sat_lat[::i], c=result[::i], cmap='viridis')  #works!
    ax.set_xlabel('Longitude [deg]', labelpad=15)
    ax.set_ylabel('Latitude [deg]', labelpad=15)
    cbar = fig.colorbar(p, pad=0.1)  #add label and units to colorbar
    cbar.set_label(variable_name+' ['+result_units+']', rotation=90, labelpad=25)
    try:  #remove pre-existing file if it exists
        os.remove(plot_file+variable_name+'_2Dheat.png')
    except:
        pass    
    plt.savefig(plot_file+variable_name+'_2Dheat.png')
    if not int_off:  #only include if in interactive mode
        if plot_close: plt.close()
        
    #plot result, lon, lat as 'time' series
    fig, axs = plt.subplots(3, figsize=(10,7))
    axs[0].plot(sat_time,result, 'black')
    axs[0].set_ylabel('Interpolated \n'+variable_name+' \n['+result_units+']', labelpad=10)
    axs[0].set_xticks([])
    axs[1].plot(sat_time,sat_lon, 'black')
    axs[1].set_ylabel('Sat Lon \n[deg]', labelpad=10)
    axs[1].set_xticks([])
    axs[2].xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y\n%H:%M:%S'))
    axs[2].xaxis.set_major_locator(MaxNLocator(5))     
    axs[2].plot([datetime.utcfromtimestamp(t) for t in sat_time], sat_lat, 'black')
    axs[2].set_ylabel('Sat Lat \n[deg]', labelpad=10)
    axs[2].set_xlabel('Satellite Time [UTC]', labelpad=10)   #add datetime for x axis
    try:  #remove pre-existing file if it exists
        os.remove(plot_file+variable_name+'_1Dlines.png')
    except:
        pass
    plt.savefig(plot_file+variable_name+'_1Dlines.png')
    if not int_off:  #only include if in interactive mode
        if plot_close: plt.close()
    return

#sample call using real flight data:
# FPlot.Plot4Dcar('T_e',sat_dict['sat_time'][sat_dict['net_idx']],sat_dict['sat_lat'][sat_dict['net_idx']],
# sat_dict['sat_lon'][sat_dict['net_idx']],sat_dict['sat_height'][sat_dict['net_idx']], sat_dict['T_e'], 
#'K', 'C:/Users/rringuet/Kamodo_WinDev1/CTIPe/Plots/', plot_close=False)
#NOTE: the net_idx is needed for real data, not sample trajectories?!
def Plot4Dcar(variable_name, sat_time, sat_lat, sat_lon, sat_z, result, result_units, 
           plot_file, sat_zunits='km', plot_close=True, plot_sampling=4000):
    '''make a 3D heat map of result in cartesian coords, with same time-series plots'''
    
    from spacepy import coordinates as coord
    from spacepy.time import Ticktock
    from astropy.constants import R_earth
    

    #turn off interactive commands based on backend in use
    if mpl.get_backend() in ['agg','pdf','ps','svg','pgf','cairo']:
        int_off = True
        plt.ioff()  
        #print('non-interactive mode ', mpl.get_backend())
    else:
        int_off = False
        #print('interactive mode ', mpl.get_backend())
    
    #Reduce number of data points for 4D plot only
    i, n = 0, len(sat_time)  #still can't handle size representation, so no s
    while n > plot_sampling: 
        i += 1
        n = len(sat_time)/i
    if i==0: i=1
    #s = 2.5**((sat_time-np.min(sat_time))/(np.max(sat_time)-np.min(sat_time))*4+1) #doesn't show
    
    #convert coordinates to cartesian, height => r in R_E, utc_ts->strings
    sat_track = [[(h*1000.+R_earth.value)/R_earth.value, lat, lon] \
                 for h, lat, lon in zip(sat_z,sat_lat,sat_lon)]
    @np.vectorize
    def ts_to_spacepystr(ts):
        return datetime.strftime(datetime.utcfromtimestamp(ts), '%Y-%m-%dT%H:%M:%S')
    coord_ticks = ts_to_spacepystr(sat_time)
    
    #perform coordinate conversion
    cvals = coord.Coords(sat_track, 'GEO', 'sph', units=['Re', 'deg', 'deg']) 
    cvals.ticks = Ticktock(coord_ticks, 'ISO')
    x, y, z = cvals.convert('GEO','car').data.T*R_earth.value/1000.  #converts to r, lat, lon 
    
    #make 4D plot, result=color
    fig = plt.figure(figsize=(10,8))
    fig.tight_layout()
    ax = fig.add_subplot(projection='3d')  #add size = time
    ax.xaxis.set_major_locator(MaxNLocator(5)) 
    ax.yaxis.set_major_locator(MaxNLocator(5)) 
    ax.zaxis.set_major_locator(MaxNLocator(5)) 
    p = ax.scatter(x[::i], y[::i], z[::i], result[::i], 
                   c=result[::i], cmap='viridis')  #works!
    ax.set_xlabel('x_GEO [km]', labelpad=15)
    ax.set_ylabel('y_GEO [km]', labelpad=15)
    ax.set_zlabel('z_GEO [km]', labelpad=15)
    cbar = fig.colorbar(p, pad=0.1)  #add label and units to colorbar
    cbar.set_label(variable_name+' ['+result_units+']', rotation=90, labelpad=25)
    try:  #remove pre-existing file if it exists
        os.remove(plot_file+variable_name+'_3Dheatcar.png')
    except:
        pass
    plt.savefig(plot_file+variable_name+'_3Dheatcar.png')
    if not int_off:  #only include if in interactive mode
        if plot_close: plt.close()
        
    #plot result, lon, lat, ilev as 'time' series
    fig, axs = plt.subplots(4, figsize=(10,7))
    axs[0].plot(sat_time,result, 'black')
    axs[0].set_ylabel('Interpolated \n'+variable_name+' \n['+result_units+']', labelpad=10)
    axs[0].set_xticks([])
    axs[1].plot(sat_time,x, 'black')
    axs[1].set_ylabel('x_GEO [km]', labelpad=10)
    axs[1].set_xticks([])
    axs[2].plot(sat_time,y, 'black')
    axs[2].set_ylabel('y_GEO [km]', labelpad=10)
    axs[2].set_xticks([])
    axs[3].xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y\n%H:%M:%S'))
    axs[3].xaxis.set_major_locator(MaxNLocator(5))     
    axs[3].plot([datetime.utcfromtimestamp(t) for t in sat_time], z, 'black')
    axs[3].set_ylabel('z_GEO [km]')
    axs[3].set_xlabel('Satellite Time [UTC]', labelpad=10)   #add datetime for x axis
    try:  #remove pre-existing file if it exists
        os.remove(plot_file+variable_name+'_1Dlinescar.png')
    except:
        pass
    plt.savefig(plot_file+variable_name+'_1Dlinescar.png')
    if not int_off:  #only include if in interactive mode
        if plot_close: plt.close()
    return