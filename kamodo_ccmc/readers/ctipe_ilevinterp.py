# -*- coding: utf-8 -*-
# The custom interpolator routine for CTIPe
from scipy.interpolate import interp1d
from kamodo import kamodofy
from numpy import NaN, vectorize, array


def PLevelInterp(kamodo_object, c3, plev_name):
    '''Custom interpolator functionto convert from altitude in km to pressure
    level.
    Parameters:
        kamodo_object - the CTIPe kamodo object
        var_name - the name of the variable as stored in the kamodo_object
        c3 - a 1D array with the grid values for pressure level
        plev_name - the name of the height function

    Output: an interpolator that accepts time, c1, c2, and pressure level and
        returns the interpolated value for the variable given
    '''
    # Retrieve correct height function
    h_func = getattr(kamodo_object, plev_name)

    @vectorize
    def km_to_ilev(t, lon, lat, km):
        '''Inputs t, lon, lat, and km are floats;
        the interpolated pressure level is returned.'''
        track = [[t, lon, lat, c3_val] for c3_val in c3]
        km_vals = h_func(track)/1000.  # interpolate for all ilev values
        km_interp = interp1d(km_vals, c3, bounds_error=False, fill_value=NaN)
        return km_interp(km)

    # Convert f(ilev) to f(km)
    @kamodofy(units='m/m')
    def plevconvert(xvec):
        '''Interpolator to convert from height to pressure level.
        Input xvec is a tuple or list of arrays/values.
        Returns the 4D position as a tuple of four 1D arrarys to feed into
        the original function.'''
        t, lon, lat, km = array(xvec).T
        out_ilev = km_to_ilev(t, lon, lat, km)
        return t, lon, lat, array(out_ilev)

    return plevconvert
