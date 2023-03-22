# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 13:25:49 2023

@author: rringuet
"""
import numpy as np
from datetime import datetime
from kamodo import Kamodo
from kamodo_ccmc.tools.functionalize import Functionalize_Dataset


def varlist(meta):
    '''Creates a list of variables names as found in the meta object returned
    by HAPI.
    Example usage: var_list = varlist(meta)'''
    var_list = [meta['parameters'][i+1]['name'] for i in range(len(
        meta['parameters'])-1)]
    return var_list


@np.vectorize
def fromiso_totimestamp(iso):
    '''Converts an iso string into an UTC timestamp.
    Example usage: timestamp_arr = fromiso_totimestamp(iso_arr)'''
    return datetime.fromisoformat(iso).timestamp()


def time_data(data):
    '''Retrieves the ISO time strings from HAPI's meta object; returns UTC
    timestamps for Kamodo.
    Example usage: time_grid = time_data(data)'''
    return fromiso_totimestamp([str(val[0])[2:].split('Z')[0].strip("'") for
                                val in data])


def kamodo_units(unit_str):
    '''Tests that unit is supported in Kamodo. Returns the corrected unit
    string if yes, or an empty string if no.
    Example usage: units = test_units(unit_str)'''
    x = np.linspace(0, 1, 10)  # sample data
    u = unit_str.replace(' ', '*')  # can't replace in place must copy
    try:
        ko = Kamodo('T['+u+']=x')
        return u
    except:
        return ''


def functionalize_variable(coord_dict, var_dict, meta_par, kamodo_object=None,
                      custom_interp=None, coord_str=''):
    '''Functionalize the given variable with an interpolator.
    Inputs:
        coord_dict: a dictionary containing the coordinate information.
            {'name_of_coord1': {'units': 'coord1_units', 'data': coord1_data},
             'name_of_coord2': {'units': 'coord2_units', 'data': coord2_data},
             etc...}
            coordX_data should be a 1D array. All others should be strings.
        var_dict: a dictionary containing the data information.
            {'variable_name1': {'units': 'data1_units', 'data': data1_array},
             'variable_name2': {'units': 'data2_units', 'data': data2_array},
             etc...}
            dataX_array should have the same shape as
                (coord1, coord2, coord3, ..., coordN)
        meta_par: meta['parameters'][i+1] where meta is the meta object
            returned from HAPI
        kamodo_object: a object created by Kamodo
        custom_interp: a custom interpolator with the same execution syntax as
            SciPy's interp1d and RegularGridInterpolator functions. The command
            ```py
            interp = custom_interp(coord_dict, var_dict)
            ```
            must initialize the interpolator with the given coordinate
            dictionary and the variable dictionary as described above.
            Default is none to functionalize with the standard interpolator.
        coord_str: a string indicating the coordinate system
            (e.g. 'GEOsph' for spherical geographical coordinates).
            See Kamodo's ConvertCoord function for more information.
            Default is an empty string.
    Returns a kamodo object with the variable added to it.
    Called by kamodofy_hapi function below.
    '''
    var_name = list(var_dict.keys())[0]
    if len(var_dict[var_name]['data'].shape) == 1:  # coord_dict fine as is
        if custom_interp is None:
            kamodo_object = Functionalize_Dataset(
                coord_dict, var_dict, kamodo_object, coord_str=coord_str)
        else:
            # initialize custom interpolation method with the variable data
            interp = custom_interp(coord_dict, var_dict)
            # functionalize
            kamodo_object = Functionalize_Dataset(
                coord_dict, var_dict, kamodo_object, coord_str=coord_str,
                func=interp, func_default='custom')
    elif len(var_dict[var_name]['data'].shape) > 1:  # add more coordinates
        coord_names = [meta_par['bins'][j]['name'] for j in
                       range(len(meta_par['bins']))]
        for j, c in enumerate(coord_names):
            coord_u = kamodo_units(meta_par['bins'][j]['units'])
            coord_dict[c] = {'data': np.array(meta_par['bins'][j]['centers']),
                             'units': coord_u}
            if len(coord_dict[c]['data'].shape) > 1:
                print('Coordinate name {c}, coordinate dimensions ' +
                      f'{coord_dict[c]["data"].shape}')
                print('Coordinate arrays of dimension higher than 1D require' +
                      ' a custom interpolator from the user. ' +
                      'See documentation for how this can be done.')
                return kamodo_object
        if custom_interp is None:
            kamodo_object = Functionalize_Dataset(
                coord_dict, var_dict, kamodo_object, coord_str=coord_str)
        else:
            # initialize custom interpolation method with the variable data
            interp = custom_interp(coord_dict, var_dict)
            # functionalize
            kamodo_object = Functionalize_Dataset(
                coord_dict, var_dict, kamodo_object, coord_str=coord_str,
                func=interp, func_default='custom')
    return kamodo_object


def functionalize_hapi(data, meta, custom_interp=None, coord_str=''):
    '''Functionalizes the data found in the data object given.
    Inputs:
        data: data object returned by HAPI
        meta: meta object returned by HAPI
        custom_interp: a custom interpolator with the same execution syntax as
            SciPy's interp1d and RegularGridInterpolator functions. If a list
            is given, the interpolators should be given in the same order as
            in the meta object. Default is None.
        custom_interp: a custom interpolator with the same execution syntax as
            SciPy's interp1d and RegularGridInterpolator functions. If a list
            is given, the interpolators should be given in the same order as
            in the meta object. The command
            ```py
            interp = custom_interp(coord_dict, var_dict)
            ```
            must initialize the interpolator with the given coordinate
            dictionary and the variable dictionary.
            Default is none to functionalize with the standard interpolator.
            Required inputs:
                coord_dict: a dictionary containing the coordinate information.
                    {'c1_name': {'units': 'c1_units', 'data': c1_data},
                     'c2_name': {'units': 'c2_units', 'data': c2_data},
                     etc...}
                    cX_data should be a 1D array. cX_name should be the desired
                    LaTeX representation of the coordinate name. cX_units
                    should be a string representing the unit to be
                    functionalized in Kamodo.
                var_dict: a dictionary containing the data information.
                    {'var_name1': {'units': 'var1_units', 'data': var1_data},
                     'var_name2': {'units': 'var2_units', 'data': var2_data},
                     etc...}
                    varX_data arrays should have the same shape as the
                    coordinates given in coord_dict
                    (e.g. (c1, c2, c3, ..., cN)).
                These objects are typically created with the following code:
                    ```py
                    time_grid = time_data(data)
                    var_list = varlist(meta)
                    for i, key in enumerate(var_list):
                        u = test_units(meta['parameters'][i+1]['units'])
                        var_dict = {key: {'data': np.array([val[i+1] for val in
                                                            data]),
                                          'units': u}}
                        coord_dict = {'UTC_time': {'units': 's',
                                                   'data': time_grid},
                                      'c1': {'units': 'R_E', 'data': .....},
                                      etc...}
                    ```
                where time_data and varlist are functions in the
                functionalize_hapi script. See the functionalize_hapi function
                logic in the same script for more details.
        coord_str: a string indicating the coordinate system
            (e.g. 'GEOsph' for spherical geographical coordinates).
            See Kamodo's ConvertCoord function for more information.
            Default is an empty string.
    Returns:
        A kamodo object with all possible variables functionalized.
    '''

    time_grid = time_data(data)  # convert iso strings to UTC timestamps
    var_list = varlist(meta)  # also used later
    kamodo_object = None  # initialize for the loop
    for i, key in enumerate(var_list):
        u = kamodo_units(meta['parameters'][i+1]['units'])  # check units Kamodo
        var_dict = {key: {'data': np.array([val[i+1] for val in data]),
                          'units': u}}
        coord_dict = {'UTC_time': {'units': 's', 'data': time_grid}}
        if isinstance(custom_interp, list):
            interp = custom_interp[i]
        else:
            interp = custom_interp
        kamodo_object = functionalize_variable(
            coord_dict, var_dict, meta['parameters'][i+1], kamodo_object,
            custom_interp=interp, coord_str=coord_str)
    return kamodo_object
