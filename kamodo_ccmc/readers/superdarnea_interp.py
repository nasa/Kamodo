# -*- coding: utf-8 -*-
# Custom interpolator script for equal area grid output of the SuperDARN model.
# given data structure is a dictionary = {lat_val: array(time, lon)}
# need an interpolator for each lat_val, then input into an interpolator for
# the outputs vs latitude.

from scipy.interpolate import RegularGridInterpolator, interp1d
from numpy import NaN, vectorize, array, append, ravel, unique, zeros
from kamodo import kamodofy
from kamodo_ccmc.readers.reader_utilities import get_slice_idx


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

    # initialize lists
    lat_interps, latidx_map = [], []

    def interp_i(t, lon, lat):
        '''t, lon, and lat are floats.'''

        print(t.shape, lon.shape, lat.shape)
        ulats = unique(lat)
        out_vals = zeros(len(t))
        for lat_val in ulats:
            # retrieve all (t, lon) pairs for given lat_val
            uidx = [i for i, l in enumerate(lat) if l == lat_val]
            position = array([t[uidx], lon[uidx]]).T
        
            # add lat interpolators as needed per lat value
            lat_idxs = get_slice_idx(lat_val, lat_grid)
            for i in lat_idxs:
                if i not in latidx_map:
                    latidx_map.append(i)
                    lat_int = RegularGridInterpolator(  # (t, lon) interpolator
                        (time_grid, lon_dict[lat_grid[i]]),
                        data_dict[lat_grid[i]],
                        bounds_error=False, fill_value=NaN)  # returns data_val
                    lat_interps.append(lat_int)
        
            if len(lat_idxs) > 1:
                interp_locations = [latidx_map.index(val) for val in lat_idxs]
                # get values for each latitude values in the index list
                # once transposed, same shape as (len(uidx), len(lat_idx))
                interp_values = array([lat_interps[i](position)
                                       for i in interp_locations]).T
                for j, vals in enumerate(interp_values):  # loop positions
                    interp_lat = interp1d(lat_grid[lat_idxs], vals,
                                          bounds_error=False, fill_value=NaN)
                    out_vals[uidx[j]] = interp_lat(lat_val)
            else: 
                interp_location = latidx_map.index(lat_idxs[0])
                out_vals[uidx] = lat_interps[interp_location](position)  # single lat_val
        
        return out_vals


    @kamodofy(units=units)
    def total_interp(xvec_SMsph3D):
        t, lon, lat = array(xvec_SMsph3D).T
        return interp_i(t, lon, lat)

    return total_interp
