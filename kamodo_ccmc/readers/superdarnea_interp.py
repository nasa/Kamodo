# -*- coding: utf-8 -*-
# Custom interpolator script for equal area grid output of the SuperDARN model.
# given data structure is a dictionary = {lat_val: array(time, lon)}
# need an interpolator for each lat_val, then input into an interpolator for
# the outputs vs latitude.

from scipy.interpolate import interp1d
from numpy import NaN, array, unique, zeros, ravel, ndarray
from kamodo_ccmc.readers.reader_utilities import get_slice_idx


def custom_interp(lon_dict, lat_grid, data_dict):
    '''Create a custom interpolator for a given variable for the superdarn
    equal area gridded output.
    Inputs:
        - time_grid: a 1D array of time values in hours.
        - longitude: a dictionary of the form {lat_val: array of lons}
        - lat_grid: a 1D array of latitude values in degrees
        - data_dict: a dictionary of the form {lat_val: array of values}
        Both input dictionaries have arrays with dimensions (time, lon).
    Outputs:
        - An kamodofied interpolator for the given data that accepts
        xvec = [[t1, lon1, lat1], [t2, lon2, lat2], ...].'''

    # initialize lists
    lat_interps, latidx_map = [], []

    def interp_i(lon, lat):
        '''lon, and lat are either 1D arrays or float values.'''
        if isinstance(lat, ndarray):
            out_vals = zeros(lon.shape)
        ulats = unique(lat)
        # retrieve all lon values for given lat_val
        for lat_val in ulats:
            if isinstance(lat, ndarray):
                uidx = [i for i, llat in enumerate(lat) if llat == lat_val]
                position = lon[uidx]
            else:
                position = lon
            # add lat interpolators as needed per lat value
            lat_idxs = get_slice_idx(lat_val, lat_grid)
            for i in lat_idxs:
                if i not in latidx_map:  # add interpolator over longitude
                    latidx_map.append(i)
                    lat_int = interp1d(
                        lon_dict[lat_grid[i]], data_dict[lat_grid[i]],
                        bounds_error=False, fill_value=NaN)  # returns data_val
                    lat_interps.append(lat_int)
            if len(lat_idxs) > 1:
                interp_locations = [latidx_map.index(val) for val in lat_idxs]
                # get values for each latitude values in the index list
                # once transposed, same shape as (len(uidx), len(lat_idx))
                interp_values = array([lat_interps[i](position)
                                       for i in interp_locations]).T
                if isinstance(lat, ndarray):
                    for j, vals in enumerate(interp_values):  # loop lons
                        interp_lat = interp1d(lat_grid[lat_idxs], vals,
                                              bounds_error=False,
                                              fill_value=NaN)
                        out_vals[uidx[j]] = interp_lat(lat_val)
                else:  # only one lon, so interp_values.shape =? lon.shape
                    interp_lat = interp1d(lat_grid[lat_idxs],
                                          ravel(interp_values),
                                          bounds_error=False, fill_value=NaN)
                    out_vals = interp_lat(lat_val)
            else:
                if isinstance(lat, ndarray):
                    interp_location = latidx_map.index(lat_idxs[0])
                    out_vals[uidx] = lat_interps[interp_location](position)
                else:
                    interp_location = latidx_map.index(lat_idxs[0])
                    out_vals = lat_interps[interp_location](position)
        return out_vals

    def total_interp(xvec):
        lon, lat = array(xvec).T
        return interp_i(lon, lat)
    return total_interp
