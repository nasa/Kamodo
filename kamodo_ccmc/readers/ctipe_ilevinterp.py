# -*- coding: utf-8 -*-
# The custom interpolator routine for CTIPe
from scipy.interpolate import interp1d
from kamodo import kamodofy
from numpy import NaN, array, unique, zeros


def PLevelInterp(kamodo_object, time, longitude, latitude, ilev, plev_name):
    '''Custom interpolator functionto convert from altitude in km to pressure
    level.
    Parameters:
        kamodo_object - the CTIPe kamodo object
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
    # Retrieve correct height function (gridded version)
    h_func = getattr(kamodo_object, plev_name+'_ijk')

    def km_to_ilev(t, lon, lat, km):
        '''Inputs t, lon, lat, and km are arrays;
        the interpolated pressure levels are returned.'''
        # sort lats and lons per time value
        input_arr = array([t, lon, lat]).T  # array of all positions given
        pos_arr = unique(input_arr, axis=0)  # only create interp for unique
        out_ilev = zeros(len(t))  # (t, lon, lat) positions to save time
        for pos in pos_arr:
            km_vals = h_func(**{'time': pos[0], 'lon': pos[1], 'lat': pos[2]})/1000.
            km_interp = interp1d(km_vals, ilev, bounds_error=False,
                                 fill_value=NaN)
            pos_idx = [i for i, p in enumerate(input_arr) if all(p == pos)]
            out_ilev[pos_idx] = km_interp(km[pos_idx])  # interp for all km
        return out_ilev
    

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

    return plevconvert, plevconvert_ijk
