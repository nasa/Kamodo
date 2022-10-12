# -*- coding: utf-8 -*-
# The custom interpolator routine for TIEGCM
from scipy.interpolate import interp1d
from kamodo import kamodofy
from numpy import NaN, vectorize, array, zeros, median


def PLevelInterp(kamodo_object, time, longitude, latitude, ilev, plev_name):
    '''Custom interpolator functionto convert from altitude in km to pressure
    level.
    Parameters:
        kamodo_object - the TIEGCM kamodo object
        time - a 1D array with the grid values for time
        longitude - a 1D array with the grid values for longitude
        latitude - a 1D array with the grid values for latitude
        ilev - a 1D array with the grid values for pressure level
        plev_name - the name of the height function (e.g. 'H_ilev')

    Output: Two interpolators and a 1D array of the median height values for
        each pressure level on the grid (in ilev). The first interpolator
        accepts time, lon, lat, and pressure level and returns time, lon, lat,
        and height. The second interpolator accepts the same arguments and
        returns only the height. Both interpolators are 'kamodofied'.
    '''
    # Retrieve correct height function
    h_func = getattr(kamodo_object, plev_name+'_ijk')
    if plev_name in ['H_ilev', 'H_ilev1', 'H_geopot']:
        conv_factor = 100000.   # convert cm to km: 10^5
    elif plev_name == 'H_milev':
        conv_factor = 1.  # already in km

    # determine the median altitude for each pressure level
    avg_kms = zeros(len(ilev))
    for i in range(len(ilev)):
        avg_kms[i] = median(h_func(**{plev_name[2:]: ilev[i]})/conv_factor)

    @vectorize
    def km_to_ilev(t, lon, lat, km):
        '''Inputs t, lon, lat, and km are floats;
        the interpolated pressure level is returned.'''
        # interpolate for all ilev values at t, lon, lat location
        km_vals = h_func(**{'time': t, 'lon': lon, 'lat': lat})/conv_factor
        km_interp = interp1d(km_vals, ilev, bounds_error=False, fill_value=NaN)
        return km_interp(km)

    # Convert f(ilev) to f(km)
    @kamodofy(units='m/m')
    def plevconvert(xvec_GDZsph4Dkm):
        '''Interpolator to convert from height to pressure level.
        Input xvec is a tuple or list of arrays/values.
        Returns the 4D position as a tuple of four 1D arrarys to feed into
        the original function.'''
        t, lon, lat, km = array(xvec_GDZsph4Dkm).T
        out_ilev = km_to_ilev(t, lon, lat, km)
        return array([t, lon, lat, out_ilev]).T

    @kamodofy(units='m/m')
    def plevconvert_ijk(xvec_GDZsph4Dkm):
        '''Interpolator to convert from height to pressure level.
        Input xvec is a tuple or list of arrays/values.
        Returns the 4D position as a tuple of four 1D arrarys to feed into
        the original function.'''
        t, lon, lat, km = array(xvec_GDZsph4Dkm).T
        return km_to_ilev(t, lon, lat, km)

    return plevconvert, plevconvert_ijk, avg_kms
