# %%
'''Function List:
conv_to_array(value): Convert object into an array.
ts_to_spacepystr(ts):  Convert from utc timestamp to datetime string.
ts_to_spacepydt(ts): Convert from utc timestamp to a datetime object.
create_spacepy(inTime, c1, c2, c3, inCoord, inType): Create a SpacePy Coords
    object and return.
create_astropy(inTime, c1, c2, c3, inCoord, inType): Create as AstroPy SkyCoord
    object and return.
spacepy_spacepy(spacepy_coord, outCoord, outType): Converts given SpacePy
    object to the requested SpacePy coordinate system.
astropy_astropy(astropy_coord, outCoord, outType): Converts given Astropy
    object to the requested AstroPy coordinate system.
spacepy_astropy(spacepy_coord, outCoord, outType): Converts given SpacePy
    object to an AstroPy object in the requested coordinate system.
astropy_spacepy(astropy_coord, outCoord, outType): Converts given Astropy
    object to a SpacePy object in the requested coordinate system.
extract_spacepy(inTime, c1, c2, c3, inCoord, inType, outCoord, outType, newC,
    verbose=False): Extract x, y, z, and units from a SpacePy Coords object and
    return.
extract_astropy(astropy_coord, outType): Extract x, y, z, and units from an
    AstroPy SkyCoord object and return.
ConvertCoord(inTime,c1,c2,c3,inCoord,inType,outCoord,outType,verbose=False):
    Main function to perform coordinate conversions.
ComputeLshell(inTime,c1,c2,c3,inCoordName,inCoordType): Computes L Shell
ConvertCoord_f(inTime,c1,c2,c3,iC1,iT1,iC2,iT2): Wrapper function for
    ConvertCoord to be compatible with integer inputs for coordinate systems
    and types.
'''

import time
import numpy as np
from datetime import datetime as dt
from datetime import timezone
import spacepy.coordinates as spc
from spacepy.coordinates import Coords
from spacepy.time import Ticktock
import astropy.units as u
from astropy.coordinates import SkyCoord


# Set option in spacepy to use ctrans instead of irbempy
spc.DEFAULTS.set_values(use_irbem=False, itol=5)

# Define coordinate lists for AstroPy and SpacePy
astropy_coordlist = ['teme', 'icrs', 'fk5', 'fk4', 'itrs', 'galactic', 'cirs',
                     'tete', 'galactocentric', 'precessedgeocentric',
                     'geocentricmeanecliptic', 'geocentrictrueecliptic',
                     'hcrs', 'barycentricmeanecliptic',
                     'heliocentricmeanecliptic', 'barycentrictrueecliptic',
                     'heliocentrictrueecliptic', 'heliocentriceclipticiau76',
                     'custombarycentricecliptic', 'lsr', 'lsrk', 'lsrd',
                     'galacticlsr', 'fk4noeterms', 'supergalactic']
spacepy_coordlist = ['GDZ', 'GEO', 'GSM', 'GSE', 'SM', 'GEI', 'MAG', 'SPH',
                     'RLL']


def convert_to_array(value):
    '''Convert given object to array.'''

    if isinstance(value, (np.ndarray)):
        return value
    else:
        if isinstance(value, (list)):
            return np.array(value)
        else:
            return np.array([value])


@np.vectorize
def ts_to_spacepystr(ts):
    '''Convert from utc timestamp to datetime string.'''
    return dt.strftime(dt.utcfromtimestamp(int(ts)), '%Y-%m-%dT%H:%M:%S')


# Steve Morley recommended converting to something other than a string to
# improve speed.
@np.vectorize
def ts_to_spacepydt(ts):
    '''Convert from utc timestamp to a datetime object.'''
    return dt.fromtimestamp(ts, tz=timezone.utc)

    '''Testing description and results:
    Given nine spherical coordinate positions...
        sat_time = np.repeat(1068865200, 9)  #UTC timestamps
        sat_lon = np.array([-180.,-180.,-180.,0.,0.,0.,180.,180.,180.])  #deg
        sat_lat = np.array([-90.,0.,90.,-90.,0.,90.,-90.,0.,90.])  #in degrees
        sat_radius = np.repeat(1.063871355909377,9)  #in R_earth
    and converting between every possible coordinate pair,
        coord_list = ['GDZ','GEO','GSM','GSE','SM','GEI','MAG','SPH','RLL']
    the ratio of completion times were compared between the original string
    method, the datetime object method and the direct timestamp method, each
    individually to the original method and also to each other (see table
    below). The ratios were chosen so the slower method was on top with the
    faster method on the bottom to indicate how many times faster the 'new'
    method would be as compared to the 'old' one.

    Avg   Std.Dev. Min   Max
    ---------------------------------
    2.05  0.67     1.00  6.2   ratio = str/dt completion times
    1.46  0.52     0.43  4.4   ratio = str/ts completion times
    1.52  0.69     0.73  4.8   ratio = ts/dt  completion times

    Result: The datetime method is the fastest of the three tested, although
    within uncertainties of the direct timestamp method. In addition, the
    datetime method is never slower than the string method unlike the timestamp
    method. Also, the datetime method was 0.55 (+/- 0.32) ms faster on average
    than the string method, and 0.24 (+/- 0.34) ms faster on average than the
    direct timestamp method. So, the datetime object method is preferred.
    '''


def create_spacepy(inTime, c1, c2, c3, inCoord, inType):
    '''Create a SpacePy Coords object and return.'''

    # pack position arrays in desired format
    if inType == "sph":
        # Need to reorder Kamodo arrays to be consistent with spacepy
        sat_track = [[c3[i], c2[i], c1[i]] for i in range(len(c1))]
    else:  # cartesian
        sat_track = [[c1[i], c2[i], c3[i]] for i in range(len(c1))]

    # create SpacePy coordinate object
    cvals = Coords(sat_track, inCoord, inType)
    tt = Ticktock(ts_to_spacepydt(inTime), 'UTC')  # datetime object method
    cvals.ticks = tt
    return cvals


def create_astropy(inTime, c1, c2, c3, inCoord, inType):
    '''Create as AstroPy SkyCoord object and return.'''

    t = ts_to_spacepydt(inTime)
    if inType == 'car':
        astropy_coord = SkyCoord(c1*u.Rearth, c2*u.Rearth, c3*u.Rearth,
                                 obstime=t, frame=inCoord,
                                 representation_type='cartesian')
    elif inType == 'sph':
        astropy_coord = SkyCoord(c1*u.deg, c2*u.deg, c3*u.Rearth, obstime=t,
                                 frame=inCoord,
                                 representation_type='spherical')
    return astropy_coord


def spacepy_spacepy(spacepy_coord, outCoord, outType):
    '''Converts given SpacePy object to the final SpacePy coordinate system.'''

    new_spacepy_coord = spacepy_coord.convert(outCoord, outType)
    assert spacepy_coord.data.shape == new_spacepy_coord.data.shape
    return new_spacepy_coord


def astropy_astropy(astropy_coord, outCoord, outType):
    '''Converts given Astropy object to the final AstroPy coordinate system.'''

    new_astropy_coord = astropy_coord.transform_to(outCoord)
    if outType == 'car':
        new_astropy_coord.representation_type = 'cartesian'
    elif outType == 'sph':
        new_astropy_coord.representation_type = 'spherical'
    return new_astropy_coord


def spacepy_astropy(spacepy_coord, outCoord, outType):
    '''Converts given SpacePy object to an AstroPy object in the final
    coordinate system.'''

    itrs = Coords.to_skycoord(spacepy_coord)  # always to ITRS cartesian
    astropy_coord = astropy_astropy(itrs, outCoord, outType)
    return astropy_coord


def astropy_spacepy(astropy_coord, outCoord, outType):
    '''Converts given Astropy object to a SpacePy object in the final
    coordinate system.'''

    cvals = Coords.from_skycoord(astropy_coord)  # always to GEO car
    spacepy_coord = spacepy_spacepy(cvals, outCoord, outType)
    return spacepy_coord


def extract_spacepy(inTime, c1, c2, c3, inCoord, inType, outCoord, outType,
                    newC, verbose=False):
    '''Extract x, y, z, and units from a SpacePy Coords object and return.'''

    if outType == "sph":
        # Need to reorder Kamodo arrays to be consistent with spacepy
        zz, yy, xx = newC.data.T
    else:
        xx, yy, zz = newC.data.T

    # check for and correct negative radii or altitudes
    # observed in some conversions to GDZ sph
    # also correct for negative longitudes in SPH sph coordinates
    # (should be 0 to 360 longitude)
    if outType == 'sph':
        idx = np.where(zz < 0.)[0]
        # neg radii error only observed to occur when lat=-90 deg in GDZ
        # also occurs at other times. ignoring those...
        close_to_pole = all(np.round(yy[idx], 1) == -90.) or \
            all(np.round(yy[idx], 1) == 90.)
        if len(idx) > 0 and not close_to_pole:
            print('Negative radii/altitudes detected away from poles.')
        if close_to_pole:
            while len(idx) > 0:
                if verbose:
                    print(f'Shifting {len(idx)} latitudes to avoid negative ' +
                          'radii or altitudes.')
                c2[idx] += 0.000000001  # slightly offset lat from -90 deg
                idx2 = np.where(c2 > 90.)[0]  # check for lat>90 after offset
                if len(idx2) > 0:
                    c2[idx2] -= 0.000000002  # correct in other direction
                sat_track_idx = [[c3[idx[i]], c2[idx[i]], c1[idx[i]]] for i in
                                 range(len(idx))]  # make new sat_track
                cvals_idx = Coords(sat_track_idx, inCoord, inType)
                tt_idx = Ticktock(ts_to_spacepydt(inTime[idx]), 'UTC')
                cvals_idx.ticks = tt_idx
                newC_idx = cvals_idx.convert(outCoord, outType)
                zz[idx], yy[idx], xx[idx] = newC_idx.data.T
                idx = np.where(zz < 0.)[0]
        if outCoord == 'SPH':  # SPH sph lon range should be 0 to 360
            idx = np.where(xx < 0)[0]  # select negative longitudes
            xx[idx] += 360  # fix
    outUnits = newC.units
    if outType == "sph":
        # Need to reorder Kamodo arrays to be consistent with spacepy
        outUnits.reverse()

    # fix spacepy bug in returned units
    if outType == 'sph':
        # check altitude for spherical return units
        if max(zz) > 300.:
            if outUnits[2] != 'km':
                outUnits[2] = 'km'
        elif min(zz) < 100.:
            if outUnits[2] != 'Re':
                outUnits[2] = 'Re'
        # else leave as it is
    else:
        # check cartesian position return units
        rr = np.sqrt(xx**2 + yy**2 + zz**2)
        if max(rr) > 300.:
            if outUnits[0] != 'km':
                outUnits = ['km', 'km', 'km']
        elif min(rr) < 100.:
            if outUnits[0] != 'Re':
                outUnits = ['Re', 'Re', 'Re']
        # else leave as it is

    # change unit Re to R_E for recognition in Kamodo
    for i in range(len(outUnits)):
        if outUnits[i] == 'Re':
            outUnits[i] = 'R_E'

    return xx, yy, zz, outUnits


def extract_astropy(astropy_coord, outType):
    '''Extract x, y, z, units from an AstroPy SkyCoord object and return.'''

    # retrieve converted values by the attribute names, converting to the
    #   standard unit based on car vs sph
    coord_names = list(astropy_coord.representation_component_names.keys())
    if outType == 'car':
        units = ['R_E', 'R_E', 'R_E']
        x = getattr(astropy_coord, coord_names[0]).to(u.Rearth).value
        y = getattr(astropy_coord, coord_names[1]).to(u.Rearth).value
    elif outType == 'sph':
        units = ['deg', 'deg', 'R_E']
        x = getattr(astropy_coord, coord_names[0]).to(u.deg).value
        y = getattr(astropy_coord, coord_names[1]).to(u.deg).value
    z = getattr(astropy_coord, coord_names[2]).to(u.Rearth).value
    return x, y, z, units


def ConvertCoord(inTime, c1, c2, c3, inCoord, inType, outCoord, outType,
                 verbose=False):
    """
    This function uses spacepy and astropy to convert time and position arrays
    from one coordinate system to another. It will correct obvious errors in
    the return units, but may not catch all incorrect values.

    INPUTS:
    inTime array:  time in UTC timestamp
    c1 array:  x (in R_earth)*, lon (in deg)
    c2 array:  y (in R_earth)*, lat (in deg)
    c3 array:  z (in R_earth)*, alt (in km), radius (in R_earth)
    inCoord string: case-sensitive string from list:
        'GDZ', 'GEO', 'GSM', 'GSE', 'SM', 'GEI', 'MAG', 'SPH', 'RLL'
        (SpacePy coordinates)
        'teme', 'icrs', 'fk5', 'fk4', 'itrs', 'galactic', 'galactocentric',
        'cirs', 'tete', 'precessedgeocentric', 'geocentricmeanecliptic',
        'geocentrictrueecliptic', 'hcrs', 'barycentricmeanecliptic',
        'heliocentricmeanecliptic', 'barycentrictrueecliptic',
        'heliocentrictrueecliptic', 'heliocentriceclipticiau76',
        'custombarycentricecliptic', 'lsr', 'lsrk', 'lsrd', 'supergalactic',
        'galacticlsr', 'fk4noeterms' (AstroPy coordinates)
        Note: Not compatible with AstroPy's HADec and AltAz coordinate systems.
        Note: Conversions using the galactocentric coordinate system are not
            conserved (a conversion to and then from this coordinate system
                       does not return the same beginning values).
    inType string: car, sph
    outCoord string: (same list as for inCoord string)
    outType string: car, sph

    OUTPUT:
    c1 array:  x (in R_earth)*, lon (in deg)
    c2 array:  y (in R_earth)*, lat (in deg)
    c3 array:  z (in R_earth)*, alt (in km), radius (in R_earth)
    units array:  [unit_c1, unit_c2, unit_c3]  (for example ['deg','deg','km']
                                                    or ['R_E','R_E','R_E'])
    *SpacePy's GDZ car coordinate system requires and produces (x, y, z) in km.

    The resource information on SpacePy's coordinate conversion function is
        sparse at best, so the below information has been collected via other
        resources and our own testing of the function. The data concerning the
        spherical coordinate systems are collected into a table format for
        easier perusal. Some details concerning the AstroPy coordinate systems
        are also below.
    For cartesian coordinates, all of the input values after time should be in
        earth radii (R_E) in order (x, y, z) to work properly, except for GDZ.
    For spherical coordinates, all of the input values after time should be in
        order (longitude, latitude, altitude or radius). The longitude and
        latitude values should be in degrees, altitude values in kilometers,
        and radius values in earth radii (R_E) from the Earth's center. All
        latitude values should fall between -90 and 90 degrees. The longitude
        range differs between the coordinate systems and is given for each in
        the table below.
    The intepretations of the input coordinates vary between the AstroPy
        coordinate systems (see below). We leave it to the user to determine
        the proper input values accordingly.
    The longitude values returned for a given coordinate converted to an
        AstroPy coordinate system are always positive (0 to 360 degrees).

    SpacePy
    Abbrev.   Full Name                       Lon. range     vertical variable
    --------------------------------------------------------------------------
    GDZ    Geodetic (WGS 84)                  (-180, 180)    Altitude (km)
    GEO    Geographic                         (-180, 180)    Radius (R_E)
    GSM    Geocentric Solar Magnetospheric    (-180, 180)    Radius (R_E)
    GSE    Geocentric Solar Ecliptic          (-180, 180)    Radius (R_E)
    SM     Solar Magnetic                     (-180, 180)    Radius (R_E)
    GEI    Geocentric Equatorial Inertial     (-180, 180)    Radius (R_E)
          (also ECI = Earth-Centered Inertial)
    MAG    Geomagnetic                        (-180, 180)    Radius (R_E)
    SPH    Spherical                            (0, 360)     Radius (R_E)
    RLL    Radius, Latitude, Longitude        (-180, 180)    Radius (R_E)

    For descriptions of most of the coordinate systems, see
    https://sscweb.gsfc.nasa.gov/users_guide/Appendix_C.shtml and
    "Geophysical Coordinate Transformations", C.T. Russell, Cosmic
        Electrodynamics, Vol. 2, pp. 184 - 196, 1971.
    The current links to SpacePy's coordinate documentation and wrapped
    conversion functions are:
        https://spacepy.github.io/autosummary/spacepy.coordinates.Coords.html
        http://svn.code.sf.net/p/irbem/code/trunk/manual/user_guide.html

    AstroPy coordinate systems:
        https://docs.astropy.org/en/stable/coordinates/skycoord.html
    Spherical coordinates in the AstroPy coordinate systems are interpreted as
    described below with the units [x,y,z] = [deg, deg, R_E].
    (All z values must be in R_E, not km.)
        teme: ['lon', 'lat', 'distance']
        icrs: ['ra', 'dec', 'distance']
        fk5: ['ra', 'dec', 'distance']
        fk4: ['ra', 'dec', 'distance']
        itrs: ['lon', 'lat', 'distance']
        galactic: ['l', 'b', 'distance']
        galactocentric: ['lon', 'lat', 'distance']
        cirs: ['ra', 'dec', 'distance']
        tete: ['ra', 'dec', 'distance']
        precessedgeocentric: ['ra', 'dec', 'distance']
        geocentricmeanecliptic: ['lon', 'lat', 'distance']
        geocentrictrueecliptic: ['lon', 'lat', 'distance']
        hcrs: ['ra', 'dec', 'distance']
        barycentricmeanecliptic: ['lon', 'lat', 'distance']
        heliocentricmeanecliptic: ['lon', 'lat', 'distance']
        barycentrictrueecliptic: ['lon', 'lat', 'distance']
        heliocentrictrueecliptic: ['lon', 'lat', 'distance']
        heliocentriceclipticiau76: ['lon', 'lat', 'distance']
        custombarycentricecliptic: ['lon', 'lat', 'distance']
        lsr: ['ra', 'dec', 'distance']
        lsrk: ['ra', 'dec', 'distance']
        lsrd: ['ra', 'dec', 'distance']
        galacticlsr: ['l', 'b', 'distance']
        fk4noeterms: ['ra', 'dec', 'distance']
        supergalactic: ['sgl', 'sgb', 'distance']
    Cartesian coordinates in the AstroPy coordinate systems are interpreted as
    (x,y,z) with a few exceptions (below) and always in R_E (earth radii).
        galactic: ['u', 'v', 'w']
        supergalactic: ['sgx', 'sgy', 'sgz']
    """

    # Start timer, define coordinate system lists in each package
    tic = time.perf_counter()
    if (inCoord not in astropy_coordlist) and (inCoord not in
                                               spacepy_coordlist):
        raise ValueError(f'The {inCoord} coordinate system is not recogized.' +
                         f'Must be one of: {spacepy_coordlist} ' +
                         f'{astropy_coordlist}.')

    # convert input values to arrays with a defined length to avoid errors.
    inTime = convert_to_array(inTime)
    c1 = convert_to_array(c1)
    c2 = convert_to_array(c2)
    c3 = convert_to_array(c3)

    # skip conversions to same coordinate system and simply return
    if inCoord == outCoord and inType == outType:
        # figure out units to return
        if verbose:
            print(f'No conversion necessary from {inCoord} {inType} to ' +
                  f'{outCoord} {outType}. Returning. ', end="")
        if inType == 'sph':
            units = ['deg', 'deg', 'R_E']
        if inType == 'sph' and inCoord == 'GDZ':
            units[2] = 'km'
        if inType == 'car':
            units = ['R_E', 'R_E', 'R_E']
        toc = time.perf_counter()
        if verbose:
            print(f'Elapsed time: {toc-tic:.4f} seconds.')
        return c1, c2, c3, units
    if len(c1) > 10000 or verbose:
        print(f'Converting {len(c1)} positions into {outCoord+outType} ' +
              'coordinates...', end="")
    # Create coordinate object
    if inCoord in spacepy_coordlist:
        coord_obj = create_spacepy(inTime, c1, c2, c3, inCoord, inType)
    elif inCoord in astropy_coordlist:
        coord_obj = create_astropy(inTime, c1, c2, c3, inCoord, inType)

    # Perform coordinate conversion
    if inCoord in spacepy_coordlist and outCoord in spacepy_coordlist:
        new_coord = spacepy_spacepy(coord_obj, outCoord, outType)
    elif inCoord in astropy_coordlist and outCoord in astropy_coordlist:
        new_coord = astropy_astropy(coord_obj, outCoord, outType)
    elif inCoord in spacepy_coordlist and outCoord in astropy_coordlist:
        new_coord = spacepy_astropy(coord_obj, outCoord, outType)
    elif inCoord in astropy_coordlist and outCoord in spacepy_coordlist:
        new_coord = astropy_spacepy(coord_obj, outCoord, outType)

    # Extract x, y, z, and units from converted coordinate object
    if outCoord in spacepy_coordlist:
        x, y, z, units = extract_spacepy(inTime, c1, c2, c3, inCoord, inType,
                                         outCoord, outType, new_coord,
                                         verbose=verbose)
    elif outCoord in astropy_coordlist:
        x, y, z, units = extract_astropy(new_coord, outType)

    # Ending messages and final return.
    toc = time.perf_counter()
    if len(c1) > 10000 or verbose:
        print(f'done in {toc-tic:0.4f} seconds.')
    return x, y, z, units


def ComputeLshell(inTime, c1, c2, c3, inCoord, inType):
    """
    This function uses spacepy to compute L shell and MLT from given time and
    position values. It uses the IRBEM library to do the calculations.
    (Not currently used in Kamodo.)

    INPUTS:
    inTime:       array:  time in UTC timestamp
    c1:           array:  x(R_E)*, lon(deg)
    c2:           array:  y(R_E)*, lat(deg)
    c3:           array:  z(R_E)*, alt(km), radius(R_E)
    *SpacePy's GDZ car coordinate system requires (x, y, z) in km.
    inCoordName:  string: case-sensitive string from list in ConvertCood. See
        that function's documentation for more details.

    OUTPUT:
    c1:       array:  Lm
    c2:       array:  Lstar
    c3:       array:  MLT
    """

    from spacepy.irbempy import get_Lstar

    # Start timer
    tic = time.perf_counter()

    # convert input values to arrays with a defined length to avoid errors
    inTime = convert_to_array(inTime)
    c1 = convert_to_array(c1)
    c2 = convert_to_array(c2)
    c3 = convert_to_array(c3)

    # create Spacepy Coord, converting from astropy coord sys if needed
    if inCoord in spacepy_coordlist:
        spacepy_coord = create_spacepy(inTime, c1, c2, c3, inCoord, inType)
        cvals = spacepy_spacepy(spacepy_coord, 'GSM', 'sph')
    elif inCoord in astropy_coordlist:
        astropy_coord = create_astropy(inTime, c1, c2, c3, inCoord, inType)
        cvals = astropy_spacepy(astropy_coord, 'GSM', 'sph')

    # compute Lstar
    # Need to test what spacepy coord sys this works fastest in
    # (GSM sph as above?)
    qq = get_Lstar(cvals.ticks, cvals)
    Lm = abs(np.reshape(qq['Lm'], -1))
    Lstar = abs(np.reshape(qq['Lstar'], -1))
    MLT = qq['MLT']

    # Ending statements and return
    toc = time.perf_counter()
    print('Computed Lm, lstar, MLT in ', "{:0.4f}".format(toc-tic),
          ' seconds for', len(inTime), ' positions.')
    return Lm, Lstar, MLT


# Not sure we need this. If so, will need to add integers for astropy
# coordinates.
def ConvertCoord_f(inTime, c1, c2, c3, iC1, iT1, iC2, iT2):
    """
    Wrapper function to ConvertCoord for Fortran calls.

    iC1,iC2 mapping:
        0 = GDZ   3 = GSE   6 = MAG
        1 = GEO   4 = SM    7 = SPH
        2 = GSM   5 = GEI   8 = RLL

    iT1,iT2 mapping:
        0 = car   1 = sph
    """

    Coords = ['GDZ', 'GEO', 'GSM', 'GSE', 'SM', 'GEI', 'MAG', 'SPH', 'RLL']
    Types = ['car', 'sph']

    inCoord = Coords[iC1]
    outCoord = Coords[iC2]
    inType = Types[iT1]
    outType = Types[iT2]

    x, y, z, units = ConvertCoord(inTime, c1, c2, c3, inCoord, inType,
                                  outCoord, outType)
    return x, y, z
