# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 11:51:15 2022

@author: rringuet
"""
import kamodo_ccmc.readers.reader_utilities as RU
from kamodo import Kamodo, kamodofy, gridify
import forge
from numpy import array, meshgrid, ravel


def Functionalize_Dataset(coord_dict, data_dict, kamodo_object=None,
                          coord_str='', func=None, func_default='data'):
    '''Determine and call the correct functionalize routine.
    Inputs:
        coord_dict: a dictionary containing the coordinate information.
            {'name_of_coord1': {'units': 'coord1_units', 'data': coord1_data},
             'name_of_coord2': {'units': 'coord2_units', 'data': coord2_data},
             etc...}
            coordX_data should be a 1D array. All others should be strings.
        data_dict: a dictionary containing the data information.
            {'variable_name1': {'units': 'data1_units', 'data': data1_array},
             'variable_name2': {'units': 'data2_units', 'data': data2_array},
             etc...}
            dataX_array should have the same shape as
                (coord1, coord2, coord3, ..., coordN)
        Note:The datasets given in the data_dict dictionary should all have the
            same dimensions. Datasets with different dimensions can be
            functionalized by simply calling the function again with the other
            dataset and the associated coordinate arrays. The datasets must
            also EACH depend upon ALL of the coordinate arrays given.
        coord_str: a string indicating the coordinate system of the data
            (e.g. "SMcar" or "GEOsph").
        kamodo_object: the previously created kamodo object. If one is not
            given, then one will be created.
        func: the function to be used for interpolation through the given
            datasets. The function must accept values for interpolation in an
            identical call structure as SciPy's RegularGridInterpolator or
            interp1D. See SciPy's documentation for more information.
        func_default: a string indicating whether a custom interpolation
            method is dersired. The default is 'data', indicating that the
            standard interpolation method will be used. Set this to 'custom' to
            indicate that func is a custom interpolator.

    Output: A kamodo object with the functionalized dataset.

    This is similar to RU.Functionalize_Dataset, except only the gridded
        interpolator is registered.
    '''
    # initialize kamodo object if None
    if kamodo_object is None:
        kamodo_object = Kamodo()

    # split the coord_dict into data and units
    coord_data = {key: value['data'] for key, value in coord_dict.items()}
    coord_units = {key: value['units'] for key, value in coord_dict.items()}

    # create the bounds array (see docstring of create_funcsig)
    # positions of corners of coordinate box
    data_list = [array([value['data'].min(), value['data'].max()]) for
                 key, value in coord_dict.items()]
    mesh_list = meshgrid(*data_list)
    bounds = array([ravel(item) for item in mesh_list], dtype=float).T

    for key in data_dict.keys():
        # create interpolator function and function signature
        interp = RU.create_interp(coord_data, data_dict[key], func=func,
                                  func_default=func_default)
        param_xvec = RU.create_funcsig(coord_data, coord_str, bounds)

        # Functionalize the dataset
        new_interp = forge.replace('xvec', param_xvec)(interp)
        interp = kamodofy(units=data_dict[key]['units'],
                          data=data_dict[key]['data'],
                          arg_units=coord_units)(new_interp)

        # Convert to gridded version (even for 1D functions) and register
        interp_grid = kamodofy(gridify(interp, **coord_data),
                               units=data_dict[key]['units'],
                               data=data_dict[key]['data'], arg_units=coord_units)
        kamodo_object[key] = interp_grid
        
    return kamodo_object
