# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 11:51:15 2022

@author: rringuet
"""
from numpy import NaN
from kamodo import Kamodo, kamodofy, gridify
from scipy.interpolate import RegularGridInterpolator, interp1d


def Functionalize_Dataset(coord_dict, data_dict, kamodo_object=None):
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
        kamodo_object: the previously created kamodo object. If one is not
            given, then one will be created.

    Output: A kamodo object with the functionalized dataset.
    '''

    # split the coord_dict into data and units
    coord_list = [value['data'] for key, value in coord_dict.items()]
    coord_data = {key: value['data'] for key, value in coord_dict.items()}
    coord_units = {key: value['units'] for key, value in coord_dict.items()}

    # determine the number of coordinates
    n_coords = len(coord_dict.keys())
    for key in data_dict.keys():  # repeat for each dataset
        if n_coords == 1:
            rgi = interp1d(coord_list[0], data_dict[key]['data'],
                              bounds_error=False, fill_value=NaN)
        else:
            rgi = RegularGridInterpolator(tuple(coord_list),
                                          data_dict[key]['data'],
                                          bounds_error=False, fill_value=NaN)
        # Functionalize the dataset
        @kamodofy(units=data_dict[key]['units'], data=data_dict[key]['data'],
                  arg_units=coord_units)
        @gridify(**coord_data)
        def interp(xvec):
            return rgi(xvec)     # return the scipy interpolator

        # add to a kamodo object and continue
        if kamodo_object is None:
            kamodo_object = Kamodo()
        kamodo_object[key] = interp
    return kamodo_object
