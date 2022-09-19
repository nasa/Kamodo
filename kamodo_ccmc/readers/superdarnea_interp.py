# -*- coding: utf-8 -*-
# Custom interpolator script for equal area grid output of the SuperDARN model.
# given data structure is a dictionary = {lat_val: array(time, lon)}
# need an interpolator for each lat_val, then input into an interpolator for
# the outputs vs latitude.

from scipy.interpolate import RegularGridInterpolator, interp1d
from numpy import NaN, vectorize, array, append, ravel
from kamodo import kamodofy


def custom_interp(time_grid, lon_dict, lat_grid, data_dict, units):
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

    # Create 2D interpolators for each latitude grid value
    lat_interps = [RegularGridInterpolator((time_grid, lon_dict[lat_val]),
                                           data_dict[lat_val],
                                           bounds_error=False, fill_value=NaN)
                   for lat_val in lat_grid]

    @vectorize
    def interp_i(t, lon, lat):
        '''t, lon, and lat are floats.'''

        # Choose indices of latitude grid values surrounding desired latitude
        lat_idx = sorted(append(lat_grid, lat)).index(lat)  # get location
        if lat_idx > 1 and lat_idx < len(lat_grid)-1:  # middle somewhere
            lat_idxs = [lat_idx-2, lat_idx-1, lat_idx, lat_idx+1]
        elif lat_idx <= 1:  # beginning
            lat_idxs = [0, 1, 2, 3]
        elif lat_idx > 1 and lat_idx >= len(lat_grid)-1:  # at end
            lat_idxs = [lat_idx-2, lat_idx-1, lat_idx]

        # Interpolate values for chosen latitude grid values
        interp_values = ravel([lat_interps[idx]([t, lon]) for idx in lat_idxs])
        lat_interp = interp1d(lat_grid[lat_idxs], interp_values,
                              bounds_error=False, fill_value=NaN)
        return lat_interp(lat)

    @kamodofy(units=units)
    def total_interp(xvec_SMsph3D):
        t, lon, lat = array(xvec_SMsph3D).T
        return interp_i(t, lon, lat)

    return total_interp
