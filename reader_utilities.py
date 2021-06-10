# -*- coding: utf-8 -*-
"""
Created on Thu May 13 18:28:18 2021

@author: rringuet
"""
from kamodo import kamodofy, gridify
from scipy.interpolate import RegularGridInterpolator
import numpy as np
from datetime import datetime, timezone

def dts_to_hrs(datetime_string, filedate):
    '''Convert datetime string to hours since midnight in filedate datetime object.'''
    
    return (datetime.strptime(datetime_string, '%Y-%m-%d %H:%M:%S').replace(tzinfo=timezone.utc)\
            -filedate).total_seconds()/3600.

def define_3d_interpolator(units,variable,t,lat,lon):
    '''Define interpolators for 3D variables uniformly'''
    
    rgi = RegularGridInterpolator((t, lat, lon),
                                  variable, bounds_error = False)
    @kamodofy(units = units, data = variable)
    def interpolator(xvec):
        """Interpolates 3d variable without a grid"""
        return rgi(xvec)
    return interpolator

def define_3d_gridded_interpolator(units,variable,t,lat,lon,xvec_dependencies):
    '''Define interpolators for 3D variables uniformly, gridded version'''

    rgi = RegularGridInterpolator((t, lat, lon), 
                                  variable, bounds_error = False)  
    if 'elat' in xvec_dependencies.keys():
        @kamodofy(units = units, data = variable, arg_units=xvec_dependencies)
        @gridify(time = t, elat = lat, elon = lon)
        def interpolator_grid(xvec):
            """Interpolates 3d variable into a grid"""
            return rgi(xvec)
    else:
        @kamodofy(units = units, data = variable, arg_units=xvec_dependencies)
        @gridify(time = t, lat = lat, lon = lon)
        def interpolator_grid(xvec):
            """Interpolates 3d variable into a grid"""
            return rgi(xvec)
    return interpolator_grid 

def define_4d_interpolator(units,variable,t,ht,lat,lon):
    
    rgi = RegularGridInterpolator((t, ht, lat, lon), 
                                  variable, bounds_error = False)
    
    @kamodofy(units = units, data = variable)
    def interpolator(xvec):
        """Interpolates 4d variable without a grid"""
        return rgi(xvec)
    return interpolator

def define_4d_gridded_interpolator(units,variable,t,ht,lat,lon,xvec_dependencies):

    rgi = RegularGridInterpolator((t, ht, lat, lon), 
                                  variable, bounds_error = False) 
    if 'ilev' in xvec_dependencies.keys():
        @kamodofy(units = units, data = variable, arg_units=xvec_dependencies)
        @gridify(time = t, ilev = ht, lat = lat, lon = lon)
        def interpolator_grid(xvec):
            """Interpolates 3d variable into a grid"""
            return rgi(xvec)
    else:
        @kamodofy(units = units, data = variable, arg_units=xvec_dependencies)
        @gridify(time = t, height = ht, lat = lat, lon = lon)
        def interpolator_grid(xvec):
            """Interpolates 3d variable into a grid"""
            return rgi(xvec)        
    return interpolator_grid 

def register_interpolator(kamodo_object, varname, interpolator, xvec_dependencies):
    '''Register interpolators for each variable'''
    
    kamodo_object[varname] = interpolator
    kamodo_object.variables[varname]['xvec'] = xvec_dependencies
    kamodo_object._registered += 1  
    return kamodo_object

def regdef_3D_interpolators(kamodo_object, units, variable, t, lat, lon, varname, 
                            xvec_dependencies, gridded_int):
    '''Register and define 3D interpolators'''
    
    interpolator = define_3d_interpolator(units,variable,t,lat,lon)
    kamodo_object = register_interpolator(kamodo_object, varname, interpolator, 
                                             xvec_dependencies)
    
    #define and register the gridded interpolator if desired
    if gridded_int:
        kamodo_object.variables[varname+'_ijk'] = dict(units = units, data = variable) 
        gridded_interpolator = define_3d_gridded_interpolator(units,variable,t,
                                                              lat,lon,xvec_dependencies)
        kamodo_object = register_interpolator(kamodo_object, varname+'_ijk', 
                                                 gridded_interpolator, xvec_dependencies)
    return kamodo_object

def regdef_4D_interpolators(kamodo_object, units, variable, t, z, lat, lon, varname, 
                            xvec_dependencies, gridded_int):
    '''Register and define 4D interpolators'''
    
    interpolator = define_4d_interpolator(units,variable,t,z,lat,lon)
    kamodo_object = register_interpolator(kamodo_object, varname, interpolator, 
                                             xvec_dependencies)
    
    #define and register the gridded interpolator if desired
    if gridded_int:
        kamodo_object.variables[varname+'_ijk'] = dict(units = units, data = variable) 
        gridded_interpolator = define_4d_gridded_interpolator(units,variable,t,z,
                                                              lat,lon,xvec_dependencies)
        kamodo_object = register_interpolator(kamodo_object, varname+'_ijk', 
                                                 gridded_interpolator, xvec_dependencies)
    return kamodo_object
    
def setup_interpolating_grids(kamodo_object, var):
    '''setup the grid for interpolation with or without height dependence'''
    
    if len(kamodo_object.variables[var]['xvec'].keys())==4: #with height dependence
        tt, xx, yy, zz = np.meshgrid(kamodo_object.newt,kamodo_object.newx,kamodo_object.newy,kamodo_object.newz, indexing = 'xy')
        kamodo_object.newgrid = np.array([np.ravel(tt), np.ravel(xx), np.ravel(yy), np.ravel(zz)]).T
    if len(kamodo_object.variables[var]['xvec'].keys())==3: #without height dependence    
        tt, yy, zz = np.meshgrid(kamodo_object.newt,kamodo_object.newy,kamodo_object.newz, indexing = 'xy')
        kamodo_object.newgrid = np.array([np.ravel(tt), np.ravel(yy), np.ravel(zz)]).T 
    #print('Gridshape', kamodo_object.newgrid.shape)
    kamodo_object.plots[kamodo_object.plottype]['newgrid'] = kamodo_object.newgrid
    return kamodo_object