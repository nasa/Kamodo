# -*- coding: utf-8 -*-
"""
Created on Thu May 13 18:28:18 2021
@author: rringuet
"""
from kamodo import kamodofy, gridify
from scipy.interpolate import RegularGridInterpolator, interp1d
from numpy import NaN
import forge
from datetime import datetime, timezone


def dts_to_hrs(datetime_string, filedate):
    '''Convert datetime string to hours since midnight in filedate datetime
    object.

    Inputs:
        datetime_string: A string of format 'YY-MM-DD HH:mm:SS', where YY is
            the two digit year, MM is the two digit month, DD is the two digit
            day, HH is the two digit hour assuming A 24 hour convention, mm is
            the two digit minute, and SS is the two digit second. The string
            should correspond to the UTC date-time to be converted into hours.
        filedate: A datetime object in UTC corresponding to midnight on the
            desired date.
    Output: The number of hours since midnight (float).
    '''

    return (datetime.strptime(datetime_string, '%Y-%m-%d %H:%M:%S'
                              ).replace(tzinfo=timezone.utc) -
            filedate).total_seconds()/3600.


def define_1d_interpolator(units, variable, t, coord_str):
    '''Define interpolators for 1D variables uniformly.

    Inputs:
        units: A string representing the units of the data in the variable
            array.
        variable: A 1D array of data values.
        t: A 1D array of utc timestamp values.
    Output: A kamodofied function to interpolate over the data given. See
        Kamodo core documentation for more details.
    '''

    rgi = interp1d(t, variable, bounds_error=False, fill_value=NaN)
    param_xvec = forge.FParameter(name='xvec_'+coord_str,
                                  interface_name='xvec',
                                  kind=forge.FParameter.POSITIONAL_OR_KEYWORD,
                                  default=t)

    @forge.replace('xvec', param_xvec)
    @kamodofy(units=units, data=variable)
    def interpolator(xvec):
        """Interpolates 1d variable without A grid"""
        return rgi(xvec)
    return interpolator


def define_2d_interpolator(units, variable, lon, lat, coord_str):
    '''Define interpolators for 2D variables uniformly.

    Inputs:
        units: A string representing the units of the data in the variable
            array.
        variable: A 3D array of data values.
        lon: A 1D array of longitude or X values. This could also be time.
        lat: A 1D array of latitude or Y values.
    Output: A kamodofied function to interpolate over the data given. See
        Kamodo core documentation for more details.
    '''

    rgi = RegularGridInterpolator((lon, lat),
                                  variable, bounds_error=False, fill_value=NaN)
    track = [[loni, lati] for loni, lati in zip(lon, lat)]
    param_xvec = forge.FParameter(name='xvec_'+coord_str,
                                  interface_name='xvec',
                                  kind=forge.FParameter.POSITIONAL_OR_KEYWORD,
                                  default=track)

    @forge.replace('xvec', param_xvec)
    @kamodofy(units=units, data=variable)
    def interpolator(xvec):
        """Interpolates 3d variable without A grid"""
        return rgi(xvec)
    return interpolator


def define_2d_gridded_interpolator(units, variable, lon, lat,
                                   xvec_dependencies, function):
    '''Define interpolators for 2D variables uniformly, gridded version.

    Inputs:
        units: A string representing the units of the data in the variable
            array.
        variable: A 3D array of data values.
        t: A 1D array of utc timestamp values.
        lon: A 1D array of longitude or X values.
        lat: A 1D array of latitude or Y values.
        xvec_dependencies: A dictionary of key, value pairs indicating the
            argument units. In this case, xvec_dependencies =
            {'time':'hr','lon':'deg',...} or similar.
        function: An ungridded kamodofied function, typically produced by the
            define_3d_interpolator function above.
    Output: A kamodofied gridded function to interpolate over the data given.
        See Kamodo core documentation for more details.
    '''

    if 'x' in xvec_dependencies.keys() or 'X' in xvec_dependencies.keys():
        interpolator_grid = kamodofy(gridify(function, x=lon, y=lat),
                                     units=units, data=variable,
                                     arg_units=xvec_dependencies)
    else:
        interpolator_grid = kamodofy(gridify(function, lon=lon, lat=lat),
                                     units=units, data=variable,
                                     arg_units=xvec_dependencies)
    return interpolator_grid


def define_3d_interpolator(units, variable, t, lon, lat, coord_str):
    '''Define interpolators for 3D variables uniformly.

    Inputs:
        units: A string representing the units of the data in the variable
            array.
        variable: A 3D array of data values.
        t: A 1D array of utc timestamp values.
        lon: A 1D array of longitude or X values.
        lat: A 1D array of latitude or Y values.
    Output: A kamodofied function to interpolate over the data given. See
        Kamodo core documentation for more details.
    '''

    rgi = RegularGridInterpolator((t, lon, lat),
                                  variable, bounds_error=False, fill_value=NaN)
    track = [[ti, loni, lati] for ti, loni, lati in zip(t, lon, lat)]
    param_xvec = forge.FParameter(name='xvec_'+coord_str,
                                  interface_name='xvec',
                                  kind=forge.FParameter.POSITIONAL_OR_KEYWORD,
                                  default=track)

    @forge.replace('xvec', param_xvec)
    @kamodofy(units=units, data=variable)
    def interpolator(xvec):
        """Interpolates 3d variable without A grid"""
        return rgi(xvec)
    return interpolator


def define_3d_gridded_interpolator(units, variable, t, lon, lat,
                                   xvec_dependencies, function):
    '''Define interpolators for 3D variables uniformly, gridded version.

    Inputs:
        units: A string representing the units of the data in the variable
            array.
        variable: A 3D array of data values.
        t: A 1D array of utc timestamp values.
        lon: A 1D array of longitude or X values.
        lat: A 1D array of latitude or Y values.
        xvec_dependencies: A dictionary of key, value pairs indicating the
            argument units. In this case, xvec_dependencies =
            {'time':'hr','lon':'deg',...} or similar.
        function: An ungridded kamodofied function, typically produced by the
            define_3d_interpolator function above.
    Output: A kamodofied gridded function to interpolate over the data given.
        See Kamodo core documentation for more details.
    '''

    if 'Elat' in xvec_dependencies.keys():
        interpolator_grid = kamodofy(gridify(function, time=t, Elon=lon,
                                             Elat=lat),
                                     units=units, data=variable,
                                     arg_units=xvec_dependencies)
    elif 'x' in xvec_dependencies.keys() or 'X' in xvec_dependencies.keys():
        interpolator_grid = kamodofy(gridify(function, time=t, x=lon, y=lat),
                                     units=units, data=variable,
                                     arg_units=xvec_dependencies)
    else:
        interpolator_grid = kamodofy(gridify(function, time=t, lon=lon,
                                             lat=lat),
                                     units=units, data=variable,
                                     arg_units=xvec_dependencies)
    return interpolator_grid


def define_4d_interpolator(units, variable, t, lon, lat, ht, coord_str):
    '''Define interpolators for 4D variables uniformly.

    Inputs:
        units: A string representing the units of the data in the variable
            array.
        variable: A 4D array of data values.
        t: A 1D array of utc timestamp values.
        lon: A 1D array of longitude or X values.
        lat: A 1D array of latitude or Y values.
        ht: A 1D array of height, radius, pressure level, or Z values.
    Output: A kamodofied function to interpolate over the data given. See
        Kamodo core documentation for more details.
    '''

    rgi = RegularGridInterpolator((t, lon, lat, ht),
                                  variable, bounds_error=False, fill_value=NaN)
    track = [[ti, loni, lati, hti] for ti, loni, lati, hti in
             zip(t, lon, lat, ht)]
    param_xvec = forge.FParameter(name='xvec_'+coord_str,
                                  interface_name='xvec',
                                  kind=forge.FParameter.POSITIONAL_OR_KEYWORD,
                                  default=track)

    @forge.replace('xvec', param_xvec)
    @kamodofy(units=units, data=variable)
    def interpolator(xvec):
        """Interpolates 4d variable without a grid"""
        return rgi(xvec)
    return interpolator


def define_4d_gridded_interpolator(units, variable, t, lon, lat, ht,
                                   xvec_dependencies, function):
    '''Define interpolators for 4D variables uniformly, gridded version.

    Inputs:
        units: A string representing the units of the data in the variable
            array.
        variable: A 4D array of data values.
        t: A 1D array of utc timestamp values.
        lon: A 1D array of longitude or X values.
        lat: A 1D array of latitude or Y values.
        ht: A 1D array of height, radius, pressure level, or Z values.
        xvec_dependencies: A dictionary of key, value pairs indicating the
            argument units. In this case, xvec_dependencies =
            {'time':'hr','lon':'deg',...} or similar.
        function: An ungridded kamodofied function, typically produced by the
            define_4d_interpolator function above.
    Output: A kamodofied gridded function to interpolate over the data given.
        See Kamodo core documentation for more details.
    '''

    if 'ilev' in xvec_dependencies.keys():
        interpolator_grid = kamodofy(gridify(function, time=t, lon=lon,
                                             lat=lat, ilev=ht), units=units,
                                     data=variable,
                                     arg_units=xvec_dependencies)
    elif 'ilev1' in xvec_dependencies.keys():
        interpolator_grid = kamodofy(gridify(function, time=t, lon=lon,
                                             lat=lat, ilev1=ht), units=units,
                                     data=variable,
                                     arg_units=xvec_dependencies)
    elif 'milev' in xvec_dependencies.keys():
        if 'mlat' in xvec_dependencies.keys() and 'mlon' in\
                xvec_dependencies.keys():
            interpolator_grid = kamodofy(gridify(function, time=t, mlon=lon,
                                                 mlat=lat, milev=ht),
                                         units=units, data=variable,
                                         arg_units=xvec_dependencies)
        else:   # not used yet by any model *******************************
            interpolator_grid = kamodofy(gridify(function, time=t, lon=lon,
                                                 lat=lat, milev=ht),
                                         units=units, data=variable,
                                         arg_units=xvec_dependencies)
    elif 'radius' in xvec_dependencies.keys():
        interpolator_grid = kamodofy(gridify(function, time=t, lon=lon,
                                             lat=lat, radius=ht), units=units,
                                     data=variable,
                                     arg_units=xvec_dependencies)
    elif 'x' in xvec_dependencies.keys() or 'X' in xvec_dependencies.keys():
        interpolator_grid = kamodofy(gridify(function, time=t, x=lon, y=lat,
                                             z=ht), units=units, data=variable,
                                     arg_units=xvec_dependencies)
    else:
        interpolator_grid = kamodofy(gridify(function, time=t, lon=lon,
                                             lat=lat, height=ht), units=units,
                                     data=variable,
                                     arg_units=xvec_dependencies)
    return interpolator_grid


def register_interpolator(kamodo_object, varname, interpolator,
                          xvec_dependencies):
    '''Register interpolators for each variable.

    Inputs:
        - kamodo_object: A kamodo object produced by the Kamodo core package.
        - varname: A string indicating the standardized variable name
            associated with the given interpolator.
        - interpolator: A kamodofied function produced by one of the functions
            above (gridded or non-gridded).
        - xvec_dependencies: A dictionary of key, value pairs indicating the
            argument units. In this case, xvec_dependencies =
            {'time':'hr','lon':'deg',...} or similar.
    Output: The same kamodo object given, except with the new function
        included.
    '''

    kamodo_object[varname] = interpolator
    kamodo_object.variables[varname]['xvec'] = xvec_dependencies
    kamodo_object._registered += 1
    return kamodo_object


def regdef_1D_interpolators(kamodo_object, units, variable, t, varname,
                            xvec_dependencies, gridded_int, coord_str):
    '''Calls all necessary functions to register and define 1D interpolators.

    Inputs:
        kamodo_object: A kamodo object produced by the Kamodo core package.
        units: A string representing the units of the data in the variable
            array.
        variable: A 1D array of data values.
        t: A 1D array of utc timestamp values.
        varname: A string indicating the standardized variable name associated
            with the given interpolator.
        xvec_dependencies: A dictionary of key, value pairs indicating the
            argument units. In this case, xvec_dependencies = {'time':'hr'}
            typically.
        gridded_int: A boolean. If True, A gridded version of the standard
            interpolator is created and registered. If False, only the standard
            interpolator is created and registered.
    Output: The same kamodo object given, except with the new function(s)
        included.
    '''

    interpolator = define_1d_interpolator(units, variable, t, coord_str)
    kamodo_object = register_interpolator(kamodo_object, varname, interpolator,
                                          xvec_dependencies)

    # do not bother with a gridded interpolator for 1D functions
    return kamodo_object


def regdef_3D_interpolators(kamodo_object, units, variable, t, lon, lat,
                            varname, xvec_dependencies, gridded_int,
                            coord_str):
    '''Calls all necessary functions to register and define 3D interpolators.

    Inputs:
        kamodo_object: A kamodo object produced by the Kamodo core package.
        units: A string representing the units of the data in the variable
            array.
        variable: A 1D array of data values.
        t: A 1D array of utc timestamp values.
        lon: A 1D array of longitude or X values.
        lat: A 1D array of latitude or Y values.
        varname: A string indicating the standardized variable name associated
            with the given interpolator.
        xvec_dependencies: A dictionary of key, value pairs indicating the
            argument units. In this case, xvec_dependencies = {'time':'hr'}
            typically.
        gridded_int: A boolean. If True, A gridded version of the standard
            interpolator is created and registered. If False, only the standard
            interpolator is created and registered.
    Output: The same kamodo object given, except with the new function(s)
        included.
    '''

    interpolator = define_3d_interpolator(units, variable, t, lon, lat,
                                          coord_str)
    kamodo_object = register_interpolator(kamodo_object, varname, interpolator,
                                          xvec_dependencies)

    # define and register the gridded interpolator if desired
    if gridded_int:
        kamodo_object.variables[varname+'_ijk'] = dict(units=units,
                                                       data=variable)
        gridded_interpolator = define_3d_gridded_interpolator(
            units, variable, t, lon, lat, xvec_dependencies, interpolator)
        kamodo_object = register_interpolator(kamodo_object, varname+'_ijk',
                                              gridded_interpolator,
                                              xvec_dependencies)
    return kamodo_object


def regdef_4D_interpolators(kamodo_object, units, variable, t, lon, lat, z,
                            varname, xvec_dependencies, gridded_int,
                            coord_str):
    '''Calls all necessary functions to register and define 4D interpolators.

    Inputs:
        kamodo_object: A kamodo object produced by the Kamodo core package.
        units: A string representing the units of the data in the variable
            array.
        variable: A 1D array of data values.
        t: A 1D array of utc timestamp values.
        lon: A 1D array of longitude or X values.
        lat: A 1D array of latitude or Y values.
        z: A 1D array of height, radius, pressure level, or Z values.
        varname: A string indicating the standardized variable name associated
            with the given interpolator.
        xvec_dependencies: A dictionary of key, value pairs indicating the
            argument units. In this case, xvec_dependencies = {'time':'hr'}
            typically.
        gridded_int: A boolean. If True, A gridded version of the standard
            interpolator is created and registered. If False, only the standard
            interpolator is created and registered.
    Output: The same kamodo object given, except with the new function(s)
        included.
    '''

    interpolator = define_4d_interpolator(units, variable, t, lon, lat, z,
                                          coord_str)
    kamodo_object = register_interpolator(kamodo_object, varname, interpolator,
                                          xvec_dependencies)

    # define and register the gridded interpolator if desired
    if gridded_int:
        kamodo_object.variables[varname+'_ijk'] = dict(units=units,
                                                       data=variable)
        gridded_interpolator = define_4d_gridded_interpolator(
            units, variable, t, lon, lat, z, xvec_dependencies, interpolator)
        kamodo_object = register_interpolator(kamodo_object, varname+'_ijk',
                                              gridded_interpolator,
                                              xvec_dependencies)
    return kamodo_object
