# -*- coding: utf-8 -*-
# The custom interpolator routine for WACCMX
from scipy.interpolate import interp1d
from kamodo import kamodofy #, get_defaults
from numpy import NaN, vectorize, array, median #, zeros, median
from time import perf_counter


def calc_km(hx_files):
    '''Calculate the median km for each pressure level, save to a h0 file.
    Input: hx_files = [h0_file, h1_file, h2_file, ...]
    Output: boolean for execution success/failure.'''

    # Z3 = H_geopot_ilev variable is in either the h1 or h2 file
    # Find the right file
    km_time = perf_counter()
    from netCDF4 import Dataset
    cdf_data = Dataset(hx_files[1])
    if 'Z3' not in cdf_data.variables.keys():
        cdf_data.close()
        cdf_data = Dataset(hx_files[2])
        if 'Z3' not in cdf_data.variables.keys():  # we have a problem
            cdf_data.close()
            print('Height function for pressure level inversion is not in ' +
                  'the h1 or h2 file. Inversion cannot be performed.')
            return False

    # retrieve variable and calculate the median height
    datatype = cdf_data.variables['Z3'].datatype
    height = array(cdf_data.variables['Z3'])
    cdf_data.close()
    km = median(array(height), axis=[0, 1, 2])/1000.  # time, lon, lat, ilev

    # save in a netCDF4 file with the name h0
    data_out = Dataset(hx_files[0], 'w', format='NETCDF3_64BIT_OFFSET')
    data_out.model = 'WACCM-X'
    new_dim = data_out.createDimension('km_ilev', len(km))
    new_var = data_out.createVariable('km_ilev', datatype, tuple(['km_ilev']))
    new_var[:] = km
    data_out.close()
    print(f'Height inversion grid calculated in {perf_counter()-km_time}s.')
    return True


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
