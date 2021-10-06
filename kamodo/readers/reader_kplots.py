# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 15:58:53 2021

@author: rringuet
"""
#import numpy as np
from numpy import meshgrid, float32, float64, ravel, array, reshape
from kamodo import kamodofy, partial, Kamodo


def convert_to_array(value):
    '''check for floats and integers. convert to arrays if found.'''
    
    type_check = [isinstance(value, float),isinstance(value, int),
                  isinstance(value, float32),isinstance(value, float64),
                  isinstance(value, list)]
    if sum(type_check)>0:
        return array([value])
    else:
        return value

def grid4D(kamodo_object, varname, time, c1, c2, c3):
    '''return data from interpolated function'''
    
    tt, xx, yy, zz = meshgrid(time, c1, c2, c3, indexing = 'xy')
    traj = array([ravel(tt), ravel(xx), ravel(yy), ravel(zz)]).T
    return getattr(kamodo_object, varname)(traj)


def grid3D(kamodo_object, varname, time, c1, c2):
    '''return data from interpolated function'''
    
    tt, xx, yy = meshgrid(time, c1, c2, indexing = 'xy')
    traj = array([ravel(tt), ravel(xx), ravel(yy)]).T
    return getattr(kamodo_object, varname)(traj)    

def plot2D(kamodo_object, varname, plottype, t, lon, lat, h=-1):
    '''Use Kamodo's native plotting to generate 2D plot.
    t, lon, lat, and h also double as t, x, y, and z for cartesian inputs.
    Possible plot types are LonLat, LatH, LonH, TimeLat, TimeLon, and TimeH for
        spherical coordinates; and TimeX, TimeY, TimeX, XY, XZ, and YZ for
        cartesian coordinates.
    If the variable depends on 4 dimensions, h should be given.
    If a LonLat plot is requested, then the function expects a single value
        (integer, float, float32, or float64) for t and h (if h is given).
        In this case, lon and lat should be 1D arrays or flat lists. Similar 
        data formatting is required for coordinates not plotted for all plot types.
    If the variable depends on height, then a value or array should be given for h.
    '''
    
    #initialize new kamodo object
    plot_kamodo=Kamodo()
    
    #first, determine if kamodo function is griddified or not, and function units
    gridified = (varname[-3:]=='ijk')
    units = kamodo_object.variables[varname]['units']
    xvec = kamodo_object.variables[varname]['xvec']
    
    #next, determine vertical dependency of variable
    coord_list = list(xvec.keys())
    if len(coord_list)==4: 
        vert = coord_list[-1]  #height, ilev, ilev1, or milev (always last)
    else:
        vert='none'
        if 'H' in plottype:
            raise AttributeError(f'Cannot produce {plottype} plot for a variable '+\
                                 f'that does not depend on height.\n{varname}: {xvec}\n')
        
    #convert inputs to arrays
    t = convert_to_array(t)
    lon = convert_to_array(lon)  #doubles as x
    lat = convert_to_array(lat)  #doubles as y
    h = convert_to_array(h)  #doubles as z
    
    #create printing message for heading of plot    
    #print(varname, plottype, units, gridified, vert)
    if t.shape[0]==1: t_message = f'Time slice at {t[0]:.3f} hrs. '
    else: t_message=''
    if lon.shape[0]==1: 
        if 'z' in vert: lon_message = f'X slice at {lon[0]:.3f} R_E. '
        else: lon_message = f'Longitude slice at {lon[0]:.3f} deg. '
    else: lon_message=''
    if lat.shape[0]==1: 
        if 'z' in vert: lat_message = f'Y slice at {lat[0]:.3f} R_E. '
        else: lat_message = f'Latitude slice at {lat[0]:.3f} deg. '
    else: lat_message=''
    if vert=='none': 
        h_message = ''
    elif h.shape[0]>1: 
        h_message = ''
    else: 
        if vert in ['ilev','ilev1','milev']: h_message = f'Pressure level slice at {h[0]}.'
        elif vert=='height': h_message = f'Height slice at {h[0]:.3f} km.'
        elif vert=='radius': h_message = f'Radius slice at {h[0]:.7f} R_E.'
        elif 'z' in vert: h_message = f'Z slice at {h[0]:.7f} R_E.'
    print(t_message+lon_message+lat_message+h_message)    
    
    #create 2D kamodo function for plotting desired plottype with given function
    if gridified:  #logic for plotting with gridified functions
        #LonLat plots
        if plottype=='LonLat':
            arg_units = {coord_list[1]:xvec[coord_list[1]], 
                         coord_list[2]:xvec[coord_list[2]]}  #e.g. {'lon':'deg','lat':'deg'}
            #accounting for differing vertical dependencies
            if vert=='height':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,height=h)
                def pfunc(time, lon, lat, height):
                    return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat,height=height)
            if vert=='radius':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,radius=h)
                def pfunc(time, lon, lat, radius):
                    return  getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat,radius=radius)
            elif vert=='ilev':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,ilev=h)
                def pfunc(time, lon, lat, ilev):
                    return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat,ilev=ilev)
            elif vert=='ilev1':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,ilev1=h)
                def pfunc(time, lon, lat, ilev1):
                    return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat,ilev1=ilev1)
            elif vert=='milev':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,milev=h)
                def pfunc(time, mlon, mlat, milev):
                    return getattr(kamodo_object, varname)(time=time,mlon=mlon,mlat=mlat,milev=milev)
                plot_kamodo['LonLat'] = pfunc  
                return plot_kamodo.plot(LonLat=dict(mlat=lat,mlon=lon))                
            elif vert=='none':
                if coord_list[1]=='Elon':  #for ctipe 3D variables that depend on Elon and Elat
                    @kamodofy(units=units, arg_units=arg_units)
                    @partial(time=t)
                    def pfunc(time, Elon, Elat): 
                        return getattr(kamodo_object, varname)(time=time,Elon=Elon,Elat=Elat)  
                    plot_kamodo['LonLat'] = pfunc  
                    return plot_kamodo.plot(LonLat=dict(Elat=lat,Elon=lon))                    
                else:
                    @kamodofy(units=units, arg_units=arg_units)
                    @partial(time=t)
                    def pfunc(time, lon, lat):
                        return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat)
            plot_kamodo['LonLat'] = pfunc  
            return plot_kamodo.plot(LonLat=dict(lat=lat,lon=lon))    
 
        #TimeLon plots
        elif plottype=='TimeLon':
            arg_units={'time':'hr',coord_list[1]:xvec[coord_list[1]]}  #'lon':'deg'
            #accounting for differing vertical dependencies
            if vert=='height':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(lat=lat,height=h)
                def pfunc(time, lon, lat, height):
                    return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat,height=height)
            elif vert=='radius':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(lat=lat,radius=h)
                def pfunc(time, lon, lat, radius):
                    return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat,radius=radius)
            elif vert=='ilev':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(lat=lat,ilev=h)
                def pfunc(time, lon, lat, ilev):
                    return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat,ilev=ilev)
            elif vert=='ilev1':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(lat=lat,ilev1=h)
                def pfunc(time, lon, lat, ilev1):
                    return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat,ilev1=ilev1)
            elif vert=='milev':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(mlat=lat,milev=h)
                def pfunc(time, mlon, mlat, milev):
                    return getattr(kamodo_object, varname)(time=time,mlon=mlon,mlat=mlat,milev=milev)
                plot_kamodo['TimeLon'] = pfunc  
                return plot_kamodo.plot(TimeLon=dict(time=t,mlon=lon))                
            elif vert=='none':
                if coord_list[1]=='Elon':  #for ctipe 3D variables that depend on Elon and Elat
                    @kamodofy(units=units, arg_units=arg_units)
                    @partial(Elat=lat)
                    def pfunc(time, Elon, Elat): 
                        return getattr(kamodo_object, varname)(time=time,Elon=Elon,Elat=Elat)  
                    plot_kamodo['TimeLon'] = pfunc  
                    return plot_kamodo.plot(TimeLon=dict(time=t,Elon=lon))                     
                else:
                    @kamodofy(units=units, arg_units=arg_units)
                    @partial(lat=lat)
                    def pfunc(time, lon, lat):
                        return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat)
            plot_kamodo['TimeLon'] = pfunc  
            return plot_kamodo.plot(TimeLon=dict(time=t,lon=lon))              
            
        #TimeLat plots
        elif plottype=='TimeLat':
            arg_units={'time':'hr',coord_list[2]:xvec[coord_list[2]]} #'lat':'deg'
            #accounting for differing vertical dependencies
            if vert=='height':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(lon=lon,height=h)
                def pfunc(time, lon, lat, height):
                    return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat,height=height)
            elif vert=='radius':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(lon=lon,radius=h)
                def pfunc(time, lon, lat, radius):
                    return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat,radius=radius)
            elif vert=='ilev':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(lon=lon,ilev=h)
                def pfunc(time, lon, lat, ilev):
                    return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat,ilev=ilev)
            elif vert=='ilev1':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(lon=lon,ilev1=h)
                def pfunc(time, lon, lat, ilev1):
                    return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat,ilev1=ilev1)
            elif vert=='milev':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(mlon=lon,milev=h)
                def pfunc(time, mlon, mlat, milev):
                    return getattr(kamodo_object, varname)(time=time,mlon=mlon,mlat=mlat,milev=milev)
                plot_kamodo['TimeLat'] = pfunc  
                return plot_kamodo.plot(TimeLat=dict(time=t,mlat=lat))
            elif vert=='none':
                if coord_list[1]=='Elon':  #for ctipe 3D variables that depend on Elon and Elat
                    @kamodofy(units=units, arg_units=arg_units)
                    @partial(Elon=lon)
                    def pfunc(time, Elon, Elat): 
                        return getattr(kamodo_object, varname)(time=time,Elon=Elon,Elat=Elat)  
                    plot_kamodo['TimeLat'] = pfunc  
                    return plot_kamodo.plot(TimeLat=dict(time=t,Elat=lat))                     
                else:
                    @kamodofy(units=units, arg_units=arg_units)
                    @partial(lon=lon)
                    def pfunc(time, lon, lat):
                        return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat)
            plot_kamodo['TimeLat'] = pfunc  
            return plot_kamodo.plot(TimeLat=dict(time=t,lat=lat))             
        
        #TimeH plots
        elif plottype=='TimeH':
            #accounting for differing vertical dependencies
            arg_units = {'time':'hr',coord_list[-1]:xvec[coord_list[-1]]}  #'time':'hr', 'height':'km'
            if vert=='height':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(lon=lon,lat=lat)
                def pfunc(time, lon, lat, height):
                    return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat,height=height)
                plot_kamodo['TimeH'] = pfunc  
                return plot_kamodo.plot(TimeH=dict(time=t,height=h))  
            elif vert=='radius':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(lon=lon,lat=lat)
                def pfunc(time, lon, lat, radius):
                    return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat,radius=radius)
                plot_kamodo['TimeH'] = pfunc  
                return plot_kamodo.plot(TimeH=dict(time=t,radius=h))
            elif vert=='ilev':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(lon=lon,lat=lat)
                def pfunc(time, lon, lat, ilev):
                    return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat,ilev=ilev)
                plot_kamodo['TimeH'] = pfunc  
                return plot_kamodo.plot(TimeH=dict(time=t,ilev=h))  
            elif vert=='ilev1':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(lon=lon,lat=lat)
                def pfunc(time, lon, lat, ilev1):
                    return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat,ilev1=ilev1)
                plot_kamodo['TimeH'] = pfunc  
                return plot_kamodo.plot(TimeH=dict(time=t,ilev1=h))  
            elif vert=='milev':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(mlon=lon,mlat=lat)
                def pfunc(time, mlon, mlat, milev):
                    return getattr(kamodo_object, varname)(time=time,mlon=mlon,mlat=mlat,milev=milev)
                plot_kamodo['TimeH'] = pfunc  
                return plot_kamodo.plot(TimeH=dict(time=t,milev=h))  
            elif vert=='none':
                raise AttributeError('Variable does not depend on height.')          
    
        #LonH plots
        elif plottype=='LonH':
            #accounting for differing vertical dependencies
            arg_units = {coord_list[1]:xvec[coord_list[1]],
                         coord_list[-1]:xvec[coord_list[-1]]}  #'lon':'deg', 'height':'km'
            if vert=='height':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,lat=lat)
                def pfunc(time, lon, lat, height):
                    return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat,height=height)
                plot_kamodo['LonH'] = pfunc  
                return plot_kamodo.plot(LonH=dict(lon=lon,height=h)) 
            elif vert=='radius':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,lat=lat)
                def pfunc(time, lon, lat, radius):
                    return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat,radius=radius)
                plot_kamodo['LonH'] = pfunc  
                return plot_kamodo.plot(LonH=dict(lon=lon,radius=h))
            elif vert=='ilev':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,lat=lat)
                def pfunc(time, lon, lat, ilev):
                    return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat,ilev=ilev)
                plot_kamodo['LonH'] = pfunc  
                return plot_kamodo.plot(LonH=dict(lon=lon,ilev=h)) 
            elif vert=='ilev1':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,lat=lat)
                def pfunc(time, lon, lat, ilev1):
                    return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat,ilev1=ilev1)
                plot_kamodo['LonH'] = pfunc  
                return plot_kamodo.plot(LonH=dict(lon=lon,ilev1=h)) 
            elif vert=='milev':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,mlat=lat)
                def pfunc(time, mlon, mlat, milev):
                    return getattr(kamodo_object, varname)(time=time,mlon=mlon,mlat=mlat,milev=milev)
                plot_kamodo['LonH'] = pfunc  
                return plot_kamodo.plot(LonH=dict(mlon=lon,milev=h)) 
            elif vert=='none':
                raise AttributeError('Variable does not depend on height.')     
        
        #LatH plots
        elif plottype=='LatH':
            #accounting for differing vertical dependencies
            arg_units = {coord_list[2]:xvec[coord_list[2]],
                         coord_list[-1]:xvec[coord_list[-1]]}  #'lat':'deg', 'height':'km'
            if vert=='height':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,lon=lon)
                def pfunc(time, lon, lat, height):
                    return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat,height=height)
                plot_kamodo['LatH'] = pfunc  
                return plot_kamodo.plot(LatH=dict(lat=lat,height=h)) 
            elif vert=='radius':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,lon=lon)
                def pfunc(time, lon, lat, radius):
                    return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat,radius=radius)
                plot_kamodo['LatH'] = pfunc  
                return plot_kamodo.plot(LatH=dict(lat=lat,radius=h))
            elif vert=='ilev':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,lon=lon)
                def pfunc(time, lon, lat, ilev):
                    return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat,ilev=ilev)
                plot_kamodo['LatH'] = pfunc  
                return plot_kamodo.plot(LatH=dict(lat=lat,ilev=h)) 
            elif vert=='ilev1':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,lon=lon)
                def pfunc(time, lon, lat, ilev1):
                    return getattr(kamodo_object, varname)(time=time,lon=lon,lat=lat,ilev1=ilev1)
                plot_kamodo['LatH'] = pfunc  
                return plot_kamodo.plot(LatH=dict(lat=lat,ilev1=h)) 
            elif vert=='milev':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,mlon=lon)
                def pfunc(time, mlon, mlat, milev):
                    return getattr(kamodo_object, varname)(time=time,mlon=mlon,mlat=mlat,milev=milev)
                plot_kamodo['LatH'] = pfunc  
                return plot_kamodo.plot(LatH=dict(mlat=lat,milev=h)) 
            elif vert=='none':
                raise AttributeError('Variable does not depend on height.')

        #cartesian plots
        elif plottype=='TimeX':
            arg_units={'time':'hr',coord_list[1]:xvec[coord_list[1]]}  #'X':'R_E'
            @kamodofy(units=units, arg_units=arg_units)
            @partial(y=lat,z=h)
            def pfunc(time, x, y, z):
                return getattr(kamodo_object, varname)(time=time,x=x,y=y,z=z)
            plot_kamodo['TimeX'] = pfunc  
            return plot_kamodo.plot(TimeX=dict(time=t,x=lon))              
        elif plottype=='TimeY':
            arg_units={'time':'hr',coord_list[2]:xvec[coord_list[2]]}  #'X':'R_E'
            @kamodofy(units=units, arg_units=arg_units)
            @partial(x=lon,z=h)
            def pfunc(time, x, y, z):
                return getattr(kamodo_object, varname)(time=time,x=x,y=y,z=z)
            plot_kamodo['TimeY'] = pfunc  
            return plot_kamodo.plot(TimeY=dict(time=t,y=lat))       
        elif plottype=='TimeZ':
            arg_units={'time':'hr',coord_list[3]:xvec[coord_list[3]]}  #'X':'R_E'
            @kamodofy(units=units, arg_units=arg_units)
            @partial(x=lon,y=lat)
            def pfunc(time, x, y, z):
                return getattr(kamodo_object, varname)(time=time,x=x,y=y,z=z)
            plot_kamodo['TimeZ'] = pfunc  
            return plot_kamodo.plot(TimeZ=dict(time=t,z=h))       
        elif plottype=='XY':
            arg_units = {coord_list[1]:xvec[coord_list[1]], 
                         coord_list[2]:xvec[coord_list[2]]}  #e.g. {'x':'R_E','y':'R_E'}
            @kamodofy(units=units, arg_units=arg_units)
            @partial(time=t,z=h)
            def pfunc(time, x, y, z):
                return getattr(kamodo_object, varname)(time=time,x=x,y=y,z=z)
            plot_kamodo['XY'] = pfunc  
            return plot_kamodo.plot(XY=dict(x=lon,y=lat))       
        elif plottype=='XZ':
            arg_units = {coord_list[1]:xvec[coord_list[1]], 
                         coord_list[3]:xvec[coord_list[3]]}  #e.g. {'x':'R_E','z':'R_E'}
            @kamodofy(units=units, arg_units=arg_units)
            @partial(time=t,y=lat)
            def pfunc(time, x, y, z):
                return getattr(kamodo_object, varname)(time=time,x=x,y=y,z=z)
            plot_kamodo['XZ'] = pfunc  
            return plot_kamodo.plot(XZ=dict(x=lon,z=h))       
        elif plottype=='YZ':
            arg_units = {coord_list[2]:xvec[coord_list[2]], 
                         coord_list[3]:xvec[coord_list[3]]}  #e.g. {'y':'R_E','z':'R_E'}
            @kamodofy(units=units, arg_units=arg_units)
            @partial(time=t,x=lon)
            def pfunc(time, x, y, z):
                return getattr(kamodo_object, varname)(time=time,x=x,y=y,z=z)
            plot_kamodo['YZ'] = pfunc  
            return plot_kamodo.plot(YZ=dict(y=lat,z=h))   
        
    else:  #logic for plotting with not gridified function-----------------------------------------------------
        #LonLat plots
        if plottype=='LonLat':
            arg_units = {coord_list[1]:xvec[coord_list[1]], 
                         coord_list[2]:xvec[coord_list[2]]}  #{'lon':'deg','lat':'deg'}
            #accounting for differing vertical dependencies
            if vert=='height':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,height=h)
                def pfunc(time, lon, lat, height):
                    data = grid4D(kamodo_object, varname, time, lon, lat, height)
                    return reshape(data,(lon.shape[0],lat.shape[0]))
            elif vert=='radius':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,radius=h)
                def pfunc(time, lon, lat, radius):
                    data = grid4D(kamodo_object, varname, time, lon, lat, radius)
                    return reshape(data,(lon.shape[0],lat.shape[0]))
            elif vert=='ilev':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,ilev=h)
                def pfunc(time, lon, lat, ilev):
                    data = grid4D(kamodo_object, varname, time, lon, lat, ilev)
                    return reshape(data,(lon.shape[0],lat.shape[0]))                
            elif vert=='ilev1':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,ilev1=h)
                def pfunc(time, lon, lat, ilev1):
                    data = grid4D(kamodo_object, varname, time, lon, lat, ilev1)
                    return reshape(data,(lon.shape[0],lat.shape[0]))                
            elif vert=='milev':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,milev=h)
                def pfunc(time, mlon, mlat, milev):
                    data = grid4D(kamodo_object, varname, time, mlon, mlat, milev)
                    return reshape(data,(mlon.shape[0],mlat.shape[0])) 
                plot_kamodo['LonLat'] = pfunc  
                return plot_kamodo.plot(LonLat=dict(mlat=lat,mlon=lon))                 
            elif vert=='none':
                if coord_list[1]=='Elon':
                    @kamodofy(units=units, arg_units=arg_units)
                    @partial(time=t)
                    def pfunc(time, Elon, Elat):
                        data = grid3D(kamodo_object, varname, time, Elon, Elat)
                        return reshape(data,(Elon.shape[0],Elat.shape[0]))     
                    plot_kamodo['LonLat'] = pfunc  
                    return plot_kamodo.plot(LonLat=dict(Elat=lat,Elon=lon))                     
                else:
                    @kamodofy(units=units, arg_units=arg_units)
                    @partial(time=t)
                    def pfunc(time, lon, lat):
                        data = grid3D(kamodo_object, varname, time, lon, lat)
                        return reshape(data,(lon.shape[0],lat.shape[0]))                       
            plot_kamodo['LonLat'] = pfunc  
            return plot_kamodo.plot(LonLat=dict(lat=lat,lon=lon))    
 
        #TimeLon plots   #####reshape command has reverse order b/c plots looked wrong for CTIPe
        elif plottype=='TimeLon':
            arg_units={'time':'hr',coord_list[1]:xvec[coord_list[1]]}  #'lon':'deg'
            #accounting for differing vertical dependencies
            if vert=='height':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(lat=lat,height=h)
                def pfunc(time, lon, lat, height):
                    data = grid4D(kamodo_object, varname, time, lon, lat, height)
                    return reshape(data,(lon.shape[0],time.shape[0])).T
            elif vert=='radius':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(lat=lat,radius=h)
                def pfunc(time, lon, lat, radius):
                    data = grid4D(kamodo_object, varname, time, lon, lat, radius)
                    return reshape(data,(lon.shape[0],time.shape[0])).T
            elif vert=='ilev':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(lat=lat,ilev=h)
                def pfunc(time, lon, lat, ilev):
                    data = grid4D(kamodo_object, varname, time, lon, lat, ilev)
                    return reshape(data,(lon.shape[0],time.shape[0])).T            
            elif vert=='ilev1':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(lat=lat,ilev1=h)
                def pfunc(time, lon, lat, ilev1):
                    data = grid4D(kamodo_object, varname, time, lon, lat, ilev1)
                    return reshape(data,(lon.shape[0],time.shape[0])).T               
            elif vert=='milev':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(mlat=lat,milev=h)
                def pfunc(time, mlon, mlat, milev):
                    data = grid4D(kamodo_object, varname, time, mlon, mlat, milev)
                    return reshape(data,(mlon.shape[0],time.shape[0])).T        
                plot_kamodo['TimeLon'] = pfunc  
                return plot_kamodo.plot(TimeLon=dict(time=t,mlon=lon))                 
            elif vert=='none':
                if coord_list[1]=='Elon':
                    @kamodofy(units=units, arg_units=arg_units)
                    @partial(Elat=lat)
                    def pfunc(time, Elon, Elat):
                        data = grid3D(kamodo_object, varname, time, Elon, Elat)
                        return reshape(data,(Elon.shape[0],time.shape[0])) 
                    plot_kamodo['TimeLon'] = pfunc  
                    return plot_kamodo.plot(TimeLon=dict(time=t,Elon=lon))                     
                else:
                    @kamodofy(units=units, arg_units=arg_units)
                    @partial(lat=lat)
                    def pfunc(time, lon, lat):
                        data = grid3D(kamodo_object, varname, time, lon, lat)
                        return reshape(data,(lon.shape[0],time.shape[0]))                     
            plot_kamodo['TimeLon'] = pfunc  
            return plot_kamodo.plot(TimeLon=dict(time=t,lon=lon))              
            
        #TimeLat plots
        elif plottype=='TimeLat':
            arg_units={'time':'hr',coord_list[2]:xvec[coord_list[2]]}  #'lat':'deg'
            #accounting for differing vertical dependencies
            if vert=='height':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(lon=lon,height=h)
                def pfunc(time, lon, lat, height):
                    data = grid4D(kamodo_object, varname, time, lon, lat, height)
                    return reshape(data,(time.shape[0],lat.shape[0]))
            elif vert=='radius':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(lon=lon,radius=h)
                def pfunc(time, lon, lat, radius):
                    data = grid4D(kamodo_object, varname, time, lon, lat, radius)
                    return reshape(data,(time.shape[0],lat.shape[0]))
            elif vert=='ilev':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(lon=lon,ilev=h)
                def pfunc(time, lon, lat, ilev):
                    data = grid4D(kamodo_object, varname, time, lon, lat, ilev)
                    return reshape(data,(time.shape[0],lat.shape[0]))               
            elif vert=='ilev1':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(lon=lon,ilev1=h)
                def pfunc(time, lon, lat, ilev1):
                    data = grid4D(kamodo_object, varname, time, lon, lat, ilev1)
                    return reshape(data,(time.shape[0],lat.shape[0]))                
            elif vert=='milev':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(mlon=lon,milev=h)
                def pfunc(time, mlon, mlat, milev):
                    data = grid4D(kamodo_object, varname, time, mlon, mlat, milev)
                    return reshape(data,(time.shape[0],mlat.shape[0]))      
                plot_kamodo['TimeLat'] = pfunc  
                return plot_kamodo.plot(TimeLat=dict(time=t,mlat=lat))                 
            elif vert=='none':
                if coord_list[1]=='Elon':
                    @kamodofy(units=units, arg_units=arg_units)
                    @partial(Elon=lon)
                    def pfunc(time, Elon, Elat):
                        data = grid3D(kamodo_object, varname, time, Elon, Elat)
                        return reshape(data,(time.shape[0],Elat.shape[0]))  
                    plot_kamodo['TimeLat'] = pfunc  
                    return plot_kamodo.plot(TimeLat=dict(time=t,Elat=lat))                      
                else:
                    @kamodofy(units=units, arg_units=arg_units)
                    @partial(lon=lon)
                    def pfunc(time, lon, lat):
                        data = grid3D(kamodo_object, varname, time, lon, lat)
                        return reshape(data,(time.shape[0],lat.shape[0]))                      
            plot_kamodo['TimeLat'] = pfunc  
            return plot_kamodo.plot(TimeLat=dict(time=t,lat=lat))             
        
        #TimeH plots
        elif plottype=='TimeH':
            #accounting for differing vertical dependencies
            arg_units = {'time':'hr',coord_list[-1]:xvec[coord_list[-1]]}  #'time':'hr', 'height':'km'
            if vert=='height':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(lon=lon,lat=lat)
                def pfunc(time, lon, lat, height):
                    data = grid4D(kamodo_object, varname, time, lon, lat, height)
                    return reshape(data,(time.shape[0],height.shape[0]))
                plot_kamodo['TimeH'] = pfunc  
                return plot_kamodo.plot(TimeH=dict(time=t,height=h))    
            elif vert=='radius':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(lon=lon,lat=lat)
                def pfunc(time, lon, lat, radius):
                    data = grid4D(kamodo_object, varname, time, lon, lat, radius)
                    return reshape(data,(time.shape[0],radius.shape[0]))
                plot_kamodo['TimeH'] = pfunc  
                return plot_kamodo.plot(TimeH=dict(time=t,radius=h))  
            elif vert=='ilev':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(lon=lon,lat=lat)
                def pfunc(time, lon, lat, ilev):
                    data = grid4D(kamodo_object, varname, time, lon, lat, ilev)
                    return reshape(data,(time.shape[0],ilev.shape[0]))    
                plot_kamodo['TimeH'] = pfunc  
                return plot_kamodo.plot(TimeH=dict(time=t,ilev=h))    
            elif vert=='ilev1':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(lon=lon,lat=lat)
                def pfunc(time, lon, lat, ilev1):
                    data = grid4D(kamodo_object, varname, time, lon, lat, ilev1)
                    return reshape(data,(time.shape[0],ilev1.shape[0]))  
                plot_kamodo['TimeH'] = pfunc  
                return plot_kamodo.plot(TimeH=dict(time=t,ilev1=h))    
            elif vert=='milev':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(mlon=lon,mlat=lat)
                def pfunc(time, mlon, mlat, milev):
                    data = grid4D(kamodo_object, varname, time, mlon, mlat, milev)
                    return reshape(data,(time.shape[0],milev.shape[0])) 
                plot_kamodo['TimeH'] = pfunc  
                return plot_kamodo.plot(TimeH=dict(time=t,milev=h))    
            elif vert=='none':
                raise AttributeError('Variable does not depend on height.')        
    
        #LonH plots
        elif plottype=='LonH':
            #accounting for differing vertical dependencies
            arg_units = {coord_list[1]:xvec[coord_list[1]],
                         coord_list[-1]:xvec[coord_list[-1]]}  #'lon':'deg', 'height':'km'
            if vert=='height':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,lat=lat)
                def pfunc(time, lon, lat, height):
                    data = grid4D(kamodo_object, varname, time, lon, lat, height)
                    return reshape(data,(lon.shape[0],height.shape[0]))
                plot_kamodo['LonH'] = pfunc  
                return plot_kamodo.plot(LonH=dict(lon=lon,height=h))  
            elif vert=='radius':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,lat=lat)
                def pfunc(time, lon, lat, radius):
                    data = grid4D(kamodo_object, varname, time, lon, lat, radius)
                    return reshape(data,(lon.shape[0],radius.shape[0]))
                plot_kamodo['LonH'] = pfunc  
                return plot_kamodo.plot(LonH=dict(lon=lon,radius=h))  
            elif vert=='ilev':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,lat=lat)
                def pfunc(time, lon, lat, ilev):
                    data = grid4D(kamodo_object, varname, time, lon, lat, ilev)
                    return reshape(data,(lon.shape[0],ilev.shape[0]))     
                plot_kamodo['LonH'] = pfunc  
                return plot_kamodo.plot(LonH=dict(lon=lon,ilev=h))  
            elif vert=='ilev1':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,lat=lat)
                def pfunc(time, lon, lat, ilev1):
                    data = grid4D(kamodo_object, varname, time, lon, lat, ilev1)
                    return reshape(data,(lon.shape[0],ilev1.shape[0]))    
                plot_kamodo['LonH'] = pfunc  
                return plot_kamodo.plot(LonH=dict(lon=lon,ilev1=h))
            elif vert=='milev':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,mlat=lat)
                def pfunc(time, mlon, mlat, milev):
                    data = grid4D(kamodo_object, varname, time, mlon, mlat, milev)
                    return reshape(data,(mlon.shape[0],milev.shape[0])) 
                plot_kamodo['LonH'] = pfunc  
                return plot_kamodo.plot(LonH=dict(mlon=lon,milev=h))
            elif vert=='none':
                raise AttributeError('Variable does not depend on height.')
        
        #LatH plots
        elif plottype=='LatH':
            #accounting for differing vertical dependencies
            arg_units = {coord_list[2]:xvec[coord_list[2]],
                         coord_list[-1]:xvec[coord_list[-1]]}  #'lat':'deg', 'height':'km'
            if vert=='height':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,lon=lon)
                def pfunc(time, lon, lat, height):
                    data = grid4D(kamodo_object, varname, time, lon, lat, height)
                    return reshape(data,(lat.shape[0],height.shape[0]))
                plot_kamodo['LatH'] = pfunc  
                return plot_kamodo.plot(LatH=dict(lat=lat,height=h)) 
            elif vert=='radius':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,lon=lon)
                def pfunc(time, lon, lat, radius):
                    data = grid4D(kamodo_object, varname, time, lon, lat, radius)
                    return reshape(data,(lat.shape[0],radius.shape[0]))
                plot_kamodo['LatH'] = pfunc  
                return plot_kamodo.plot(LatH=dict(lat=lat,radius=h)) 
            elif vert=='ilev':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,lon=lon)
                def pfunc(time, lon, lat, ilev):
                    data = grid4D(kamodo_object, varname, time, lon, lat, ilev)
                    return reshape(data,(lat.shape[0],ilev.shape[0]))   
                plot_kamodo['LatH'] = pfunc  
                return plot_kamodo.plot(LatH=dict(lat=lat,ilev=h)) 
            elif vert=='ilev1':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,lon=lon)
                def pfunc(time, lon, lat, ilev1):
                    data = grid4D(kamodo_object, varname, time, lon, lat, ilev1)
                    return reshape(data,(lat.shape[0],ilev1.shape[0]))   
                plot_kamodo['LatH'] = pfunc  
                return plot_kamodo.plot(LatH=dict(lat=lat,ilev1=h)) 
            elif vert=='milev':
                @kamodofy(units=units, arg_units=arg_units)
                @partial(time=t,mlon=lon)
                def pfunc(time, mlon, mlat, milev):
                    data = grid4D(kamodo_object, varname, time, mlon, mlat, milev)
                    return reshape(data,(mlat.shape[0],milev.shape[0])) 
                plot_kamodo['LatH'] = pfunc  
                return plot_kamodo.plot(LatH=dict(mlat=lat,milev=h)) 
            elif vert=='none':
                raise AttributeError('Variable does not depend on height.')

        #cartesian plots
        elif plottype=='TimeX':
            arg_units={'time':'hr',coord_list[1]:xvec[coord_list[1]]}  #'X':'R_E'
            @kamodofy(units=units, arg_units=arg_units)
            @partial(y=lat,z=h)
            def pfunc(time, x, y, z):
                data = grid4D(kamodo_object, varname, time, x, y, z)
                return reshape(data,(x.shape[0],time.shape[0])).T                 
            plot_kamodo['TimeX'] = pfunc  
            return plot_kamodo.plot(TimeX=dict(time=t,x=lon))              
        elif plottype=='TimeY':
            arg_units={'time':'hr',coord_list[2]:xvec[coord_list[2]]}  #'X':'R_E'
            @kamodofy(units=units, arg_units=arg_units)
            @partial(x=lon,z=h)
            def pfunc(time, x, y, z):
                data = grid4D(kamodo_object, varname, time, x, y, z)
                return reshape(data,(time.shape[0],y.shape[0]))  
            plot_kamodo['TimeY'] = pfunc  
            return plot_kamodo.plot(TimeY=dict(time=t,y=lat))       
        elif plottype=='TimeZ':
            arg_units={'time':'hr',coord_list[3]:xvec[coord_list[3]]}  #'X':'R_E'
            @kamodofy(units=units, arg_units=arg_units)
            @partial(x=lon,y=lat)
            def pfunc(time, x, y, z):
                data = grid4D(kamodo_object, varname, time, x, y, z)
                return reshape(data,(time.shape[0],z.shape[0]))  
            plot_kamodo['TimeZ'] = pfunc  
            return plot_kamodo.plot(TimeZ=dict(time=t,z=h))       
        elif plottype=='XY':
            arg_units = {coord_list[1]:xvec[coord_list[1]], 
                         coord_list[2]:xvec[coord_list[2]]}  #e.g. {'x':'R_E','y':'R_E'}
            @kamodofy(units=units, arg_units=arg_units)
            @partial(time=t,z=h)
            def pfunc(time, x, y, z):
                data = grid4D(kamodo_object, varname, time, x, y, z)
                return reshape(data,(x.shape[0],y.shape[0]))  
            plot_kamodo['XY'] = pfunc  
            return plot_kamodo.plot(XY=dict(x=lon,y=lat))       
        elif plottype=='XZ':
            arg_units = {coord_list[1]:xvec[coord_list[1]], 
                         coord_list[3]:xvec[coord_list[3]]}  #e.g. {'x':'R_E','z':'R_E'}
            @kamodofy(units=units, arg_units=arg_units)
            @partial(time=t,y=lat)
            def pfunc(time, x, y, z):
                data = grid4D(kamodo_object, varname, time, x, y, z)
                return reshape(data,(x.shape[0],z.shape[0]))  
            plot_kamodo['XZ'] = pfunc  
            return plot_kamodo.plot(XZ=dict(x=lon,z=h))       
        elif plottype=='YZ':
            arg_units = {coord_list[2]:xvec[coord_list[2]], 
                         coord_list[3]:xvec[coord_list[3]]}  #e.g. {'y':'R_E','z':'R_E'}
            @kamodofy(units=units, arg_units=arg_units)
            @partial(time=t,x=lon)
            def pfunc(time, x, y, z):
                data = grid4D(kamodo_object, varname, time, x, y, z)
                return reshape(data,(y.shape[0],z.shape[0]))  
            plot_kamodo['YZ'] = pfunc  
            return plot_kamodo.plot(YZ=dict(y=lat,z=h))     
