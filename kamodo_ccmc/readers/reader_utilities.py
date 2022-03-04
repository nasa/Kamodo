# -*- coding: utf-8 -*-
"""
Created on Thu May 13 18:28:18 2021

@author: rringuet
"""
from kamodo import kamodofy, gridify
from scipy.interpolate import RegularGridInterpolator
from numpy import NaN
from datetime import datetime, timezone

def dts_to_hrs(datetime_string, filedate):
    '''Convert datetime string to hours since midnight in filedate datetime object.'''
    
    return (datetime.strptime(datetime_string, '%Y-%m-%d %H:%M:%S').replace(tzinfo=timezone.utc)\
            -filedate).total_seconds()/3600.

def define_3d_interpolator(units,variable,t,lon,lat):
    '''Define interpolators for 3D variables uniformly'''
    
    rgi = RegularGridInterpolator((t, lon, lat),
                                  variable, bounds_error = False, fill_value=NaN)
    @kamodofy(units=units, data=variable)
    def interpolator(xvec):
        """Interpolates 3d variable without a grid"""
        return rgi(xvec)
    return interpolator

def define_3d_gridded_interpolator(units,variable,t,lon,lat,xvec_dependencies, 
                                   function):
    '''Define interpolators for 3D variables uniformly, gridded version'''

    if 'Elat' in xvec_dependencies.keys():
        interpolator_grid = kamodofy(gridify(function, time = t, Elon = lon, Elat = lat), 
                                     units=units, data=variable, arg_units=xvec_dependencies)  
    elif 'x' in xvec_dependencies.keys() or 'X' in xvec_dependencies.keys():
        interpolator_grid = kamodofy(gridify(function, time = t, x = lon, y = lat), 
                                     units=units, data=variable, arg_units=xvec_dependencies)         
    else:
        interpolator_grid = kamodofy(gridify(function, time = t, lon = lon, lat = lat), 
                                     units=units, data=variable, arg_units=xvec_dependencies)          
    return interpolator_grid 

def define_4d_interpolator(units,variable,t,lon,lat,ht):
    
    rgi = RegularGridInterpolator((t, lon, lat, ht), 
                                  variable, bounds_error = False, fill_value=NaN)
    
    @kamodofy(units=units, data=variable)
    def interpolator(xvec):
        """Interpolates 4d variable without a grid"""
        return rgi(xvec)
    return interpolator

def define_4d_gridded_interpolator(units,variable,t,lon,lat,ht,xvec_dependencies,
                                   function):

    #rgi = RegularGridInterpolator((t, lon, lat, ht), 
    #                              variable, bounds_error = False, fill_value=NaN) 
    if 'ilev' in xvec_dependencies.keys():
        interpolator_grid = kamodofy(gridify(function, time=t, lon=lon, lat=lat, ilev=ht), 
                                     units=units, data=variable, arg_units=xvec_dependencies)        
    elif 'ilev1' in xvec_dependencies.keys():
        interpolator_grid = kamodofy(gridify(function, time=t, lon=lon, lat=lat, ilev1=ht), 
                                     units=units, data=variable, arg_units=xvec_dependencies)         
    elif 'milev' in xvec_dependencies.keys():
        if 'mlat' in xvec_dependencies.keys() and 'mlon' in xvec_dependencies.keys():
            interpolator_grid = kamodofy(gridify(function, time = t, mlon=lon, mlat = lat, milev = ht), 
                                         units=units, data=variable, arg_units=xvec_dependencies)             
        else:   #not used yet by any model *******************************
            interpolator_grid = kamodofy(gridify(function, time = t, lon=lon, lat = lat, milev = ht), 
                                         units=units, data=variable, arg_units=xvec_dependencies)           
    elif 'radius' in xvec_dependencies.keys():
        interpolator_grid = kamodofy(gridify(function, time = t, lon=lon, lat = lat,  radius = ht), 
                                     units=units, data=variable, arg_units=xvec_dependencies)   
    elif 'x' in xvec_dependencies.keys() or 'X' in xvec_dependencies.keys():
        interpolator_grid = kamodofy(gridify(function, time = t, x=lon, y = lat, z = ht), 
                                     units=units, data=variable, arg_units=xvec_dependencies)      
    else:
        interpolator_grid = kamodofy(gridify(function, time = t, lon=lon, lat = lat,  height = ht), 
                                     units=units, data=variable, arg_units=xvec_dependencies)          
    return interpolator_grid 

def register_interpolator(kamodo_object, varname, interpolator, xvec_dependencies):
    '''Register interpolators for each variable'''
    
    kamodo_object[varname] = interpolator
    kamodo_object.variables[varname]['xvec'] = xvec_dependencies
    kamodo_object._registered += 1  
    return kamodo_object

def regdef_3D_interpolators(kamodo_object, units, variable, t, lon, lat, varname, 
                            xvec_dependencies, gridded_int):
    '''Register and define 3D interpolators'''
    
    interpolator = define_3d_interpolator(units,variable,t,lon,lat)
    kamodo_object = register_interpolator(kamodo_object, varname, interpolator, 
                                             xvec_dependencies)
    
    #define and register the gridded interpolator if desired
    if gridded_int:
        kamodo_object.variables[varname+'_ijk'] = dict(units = units, data = variable) 
        gridded_interpolator = define_3d_gridded_interpolator(units,variable,t,
                                                              lon,lat,xvec_dependencies,
                                                              interpolator)
        kamodo_object = register_interpolator(kamodo_object, varname+'_ijk', 
                                                 gridded_interpolator, xvec_dependencies)
    return kamodo_object

def regdef_4D_interpolators(kamodo_object, units, variable, t, lon, lat, z, varname, 
                            xvec_dependencies, gridded_int):
    '''Register and define 4D interpolators'''
    
    interpolator = define_4d_interpolator(units,variable,t,lon,lat,z)
    kamodo_object = register_interpolator(kamodo_object, varname, interpolator, 
                                             xvec_dependencies)
    
    #define and register the gridded interpolator if desired
    if gridded_int:
        kamodo_object.variables[varname+'_ijk'] = dict(units = units, data = variable) 
        gridded_interpolator = define_4d_gridded_interpolator(units,variable,t,lon,
                                                              lat,z,xvec_dependencies,
                                                              interpolator)
        kamodo_object = register_interpolator(kamodo_object, varname+'_ijk', 
                                                 gridded_interpolator, xvec_dependencies)
    return kamodo_object
