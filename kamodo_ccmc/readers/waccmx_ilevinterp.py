# -*- coding: utf-8 -*-
# The custom interpolator routine for WACCMX
from scipy.interpolate import interp1d
from kamodo import kamodofy #, get_defaults
from numpy import NaN, vectorize, array


def PLevelInterp(h_func, time, longitude, latitude, ilev):
    '''Custom interpolator functionto convert from altitude in km to pressure
    level.
    Parameters:
        h_func - the H_geopot_ilev function/interpolator
        time - a 1D array with the grid values for time
        longitude - a 1D array with the grid values for longitude
        latitude - a 1D array with the grid values for latitude
        ilev - a 1D array with the grid values for pressure level

    Output: Two interpolators and a 1D array of the median height values for
        each pressure level on the grid (in ilev). The first interpolator
        accepts time, lon, lat, and pressure level and returns time, lon, lat,
        and height. The second interpolator accepts the same arguments and
        returns only the height. Both interpolators are 'kamodofied'.
    '''

    @vectorize
    def km_to_ilev(t, lon, lat, km):
        '''Inputs t, lon, lat, and km are floats;
        the interpolated pressure level is returned.'''
        # interpolate for all ilev values
        km_vals = h_func(**{'time': t, 'lon': lon, 'lat': lat})/1000.  # km
        km_interp = interp1d(km_vals, ilev, bounds_error=False, fill_value=NaN)
        return km_interp(km)

    # Convert f(ilev) to f(km)
    @kamodofy(units='hPa')
    def plevconvert(xvec_GDZsph4Dkm):
        '''Interpolator to convert from height to pressure level.
        Input xvec is a tuple or list of arrays/values.
        Returns the 4D position as a tuple of four 1D arrarys to feed into
        the original function.'''
        t, lon, lat, km = array(xvec_GDZsph4Dkm).T
        out_ilev = km_to_ilev(t, lon, lat, km)
        return array([t, lon, lat, out_ilev]).T

    @kamodofy(units='hPa')
    def plevconvert_ijk(xvec_GDZsph4Dkm):
        '''Interpolator to convert from height to pressure level.
        Input xvec is a tuple or list of arrays/values.
        Returns the 4D position as a tuple of four 1D arrarys to feed into
        the original function.'''
        t, lon, lat, km = array(xvec_GDZsph4Dkm).T
        return km_to_ilev(t, lon, lat, km)

    return plevconvert, plevconvert_ijk
