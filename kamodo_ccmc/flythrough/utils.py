# %%
import time
import numpy as np
from datetime import datetime as dt
from datetime import timezone
from spacepy.coordinates import Coords
from spacepy.time import Ticktock
from spacepy.irbempy import get_Lstar

def ConvertCoord_f(inTime,c1,c2,c3,iC1,iT1,iC2,iT2):
    """
    Wrapper function to ConvertCoord for Fortran calls.
    
    iC1,iC2 mapping: 
        0 = GDZ   3 = GSE   6 = MAG
        1 = GEO   4 = SM    7 = SPH
        2 = GSM   5 = GEI   8 = RLL

    iT1,iT2 mapping: 
        0 = car   1 = sph
    """
    
    Coords=['GDZ','GEO','GSM','GSE','SM','GEI','MAG','SPH','RLL']
    Types=['car','sph']

    inCoord  = Coords[iC1]
    outCoord = Coords[iC2]
    inType  = Types[iT1]
    outType = Types[iT2]

    x,y,z,units = ConvertCoord(inTime,c1,c2,c3,inCoord,inType,outCoord,outType)
    return x,y,z

def ConvertCoord(inTime,c1,c2,c3,inCoord,inType,outCoord,outType,verbose=False):
    """
    This function uses spacepy to convert time and position arrays from one coordinate system to another.
    It will correct obvious errors in the return units, but may not catch all incorrect values.
    
    INPUTS:
    inTime:   array:  time in UTC timestamp
    c1:       array:  x, lon
    c2:       array:  y, lat
    c3:       array:  z, alt
    inCoord:  string: GDZ, GEO, GSM, GSE, SM, GEI, MAG, SPH, RLL
    inType:   string: car, sph
    outCoord: string: GDZ, GEO, GSM, GSE, SM, GEI, MAG, SPH, RLL
    outType:  string: car, sph
    
    OUTPUT:
    c1:       array:  x, alt
    c2:       array:  y, lat
    c3:       array:  z, lon
    units:    array:  [unit_c1, unit_c2, unit_c3]  (for example ['km','deg','deg'] 
                                                    or ['R_E','R_E','R_E'])

    The resource information on SpacePy's coordinate conversion function is 
        sparse at best, so the below information has been collected via other
        resources and our own testing of the function. The data concerning the
        spherical coordinate systems are collected into a table format for 
        easier perusal.
    For cartesian coordinates, all of the input values after time should be in 
        earth radii (R_E) in order (x, y, z) to work properly.
    For spherical coordinates, all of the input values after time should be in order 
        (longitude, latitude, altitude or radius). The longitude and latitude
        values should be in degrees, altitude values in kilometers, and radius
        values in earth radii (R_E) from the Earth's center. All latitude values 
        should fall between -90 and 90 degrees. The longitude range differs 
        between the coordinate systems and is given for each in the table below.
        
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
    https://sscweb.gsfc.nasa.gov/users_guide/Appendix_C.shtml and it's reference,
    "Geophysical Coordinate Transformations", C.T. Russell, Cosmic Electrodynamics, Vol. 2, pp. 184 - 196, 1971.
    The current links to SpacePy's coordinate documentation and wrapped conversion functions are:
        https://spacepy.github.io/autosummary/spacepy.coordinates.Coords.html
        http://svn.code.sf.net/p/irbem/code/trunk/manual/user_guide.html

    """

    # Start timer
    tic = time.perf_counter()
    #print("Starting ...")

    #skip conversions to same coordinate system and simply return
    if inCoord==outCoord and inType==outType:
        #figure out units to return
        if verbose: print(f'No conversion necessary from {inCoord} {inType} to {outCoord} {outType}. Returning. ', end="")
        if inType=='sph': units=['deg','deg','R_E']
        if inType=='sph' and inCoord=='GDZ': units[2]='km'
        if inType=='car': units=['R_E','R_E','R_E']
        toc = time.perf_counter()
        if verbose: print(f'Elapsed time: {toc-tic:.4f} seconds.')
        return c1, c2, c3, units

    #pack position arrays in desired format
    if inType == "sph":
        # Need to reorder Kamodo arrays to be consistent with spacepy
        sat_track = [[c3[i], c2[i], c1[i] ] for i in range(len(c1))]
    else:
        sat_track = [[c1[i], c2[i], c3[i] ] for i in range(len(c1))]

    #convert utc_ts to spacepy time strings
    @np.vectorize
    def ts_to_spacepystr(ts):
        return dt.strftime(dt.utcfromtimestamp(int(ts)), '%Y-%m-%dT%H:%M:%S')
    
    #Steve Morley recommended converting to something other that a string to improve speed
    @np.vectorize
    def ts_to_spacepydt(ts):
        return dt.fromtimestamp(ts, tz=timezone.utc)

    #perform coordinate conversion
    cvals = Coords(sat_track, inCoord, inType)
    #tt = Ticktock(ts_to_spacepystr(inTime), 'ISO')  #original string method
    tt = Ticktock(ts_to_spacepydt(inTime), 'UTC')  #datetime object method
    #tt = Ticktock(inTime, 'UNX')  #direct timestamp method
    
    '''Testing description and results:
    Given nine spherical coordinate positions...
        sat_time = np.repeat(1068865200, 9)  #UTC timestamps
        sat_lon = np.array([-180.,-180.,-180.,0.,0.,0.,180.,180.,180.])  #in degrees
        sat_lat = np.array([-90.,0.,90.,-90.,0.,90.,-90.,0.,90.])  #in degrees
        sat_radius = np.repeat(1.063871355909377,9)  #in R_earth        
    and converting between every possible coordinate pair,
        coord_list = ['GDZ','GEO','GSM','GSE','SM','GEI','MAG','SPH','RLL']
    the ratio of completion times were compared between the original string method,
    the datetime object method and the direct timestamp method, each individually
    to the original method and also to each other (see table below). The ratios
    were chosen so the slower method was on top with the faster method on the 
    bottom to indicate how many times faster the 'new' method would be as compared
    to the 'old' one.
    
    Avg   Std.Dev. Min   Max
    ---------------------------------
    2.05  0.67     1.00  6.2   ratio = str/dt completion times
    1.46  0.52     0.43  4.4   ratio = str/ts completion times
    1.52  0.69     0.73  4.8   ratio = ts/dt  completion times
    
    Result: The datetime method is the fastest of the three tested, although 
    within uncertainties of the direct timestamp method. In addition, the datetime
    method is never slower than the string method unlike the timestamp method.
    Also, the datetime method was 0.55 (+/- 0.32) ms faster on average than the 
    string method, and 0.24 (+/- 0.34) ms faster on average than the direct timestamp 
    method. So, the datetime object method is preferred.
    
    '''
    
    cvals.ticks = tt
    newC = cvals.convert(outCoord,outType)
    #print(newC)
    if outType == "sph":
        # Need to reorder Kamodo arrays to be consistent with spacepy
        zz, yy, xx = newC.data.T
    else:
        xx, yy, zz = newC.data.T

    #check for and correct negative radii or altitudes, observed in some conversions to GDZ sph
    #also correct for negative longitudes in SPH sph coordinates (should be 0 to 360 longitude)
    if outType == 'sph':
        idx = np.where(zz<0.)[0]
        while len(idx)>0:  #neg radii error only observed to occur when lat=-90 deg in GDZ
            if verbose: print(f'Shifting {len(idx)} latitudes to avoid negative radii or altitudes.')
            c2[idx] += 0.000000001  #slightly offset lat from -90 deg
            idx2 = np.where(c2>90.)[0]  #safety check for lat>90 after offset
            if len(idx2)>0: c2[idx2] -= 0.000000002  #if so, correct in other direction
            sat_track_idx = [[c3[idx[i]],c2[idx[i]],c1[idx[i]] ] for i in range(len(idx))]  #make new sat_track
            cvals_idx = Coords(sat_track_idx, inCoord, inType)
            tt_idx = Ticktock(ts_to_spacepydt(inTime[idx]), 'UTC') 
            cvals_idx.ticks = tt_idx
            newC_idx = cvals_idx.convert(outCoord,outType)
            zz[idx], yy[idx], xx[idx] = newC_idx.data.T
            #print(zz[idx])
            idx = np.where(zz<0.)[0]
        if outCoord=='SPH': #SPH sph lon range should be 0 to 360
            idx = np.where(xx<0)[0]  #select negative longitudes
            xx[idx] += 360  #fix

    outUnits = newC.units
    if outType == "sph":
        # Need to reorder Kamodo arrays to be consistent with spacepy
        outUnits.reverse()

    #fix spacepy bug in returned units
    if outType == 'sph':
        #check altitude for spherical return units
        if max(zz)>300.:
            if outUnits[2] != 'km':
                #print(' -changing return alt units from',outUnits[2],'to km')
                outUnits[2]='km'
        elif min(zz)<100.:
            if outUnits[2] != 'Re':
                #print(' -changing return alt units from',outUnits[2],'to Re')
                outUnits[2]='Re'
        #else leave as it is
    else:
        #check cartesian position return units
        rr=np.sqrt(xx**2 + yy**2 + zz**2)
        if max(rr)>300.:
            if outUnits[0] != 'km':
                #print(' -changing spacepy return position units from',outUnits[0],'to km')
                outUnits[0]='km'
                outUnits[1]='km'
                outUnits[2]='km'
        elif min(rr)<100.:
            if outUnits[0] != 'Re':
                #print(' -changing spacepy return position units from',outUnits[0],'to Re')
                outUnits[0]='Re'
                outUnits[1]='Re'
                outUnits[2]='Re'
        #else leave as it is

    #change unit Re to R_E for recognition in Kamodo
    for i in range(len(outUnits)):
        if outUnits[i]=='Re':
            outUnits[i]='R_E'

    toc = time.perf_counter()
    if verbose: print('Converted from ',inCoord,inType,'to:',outCoord,outType,outUnits,'in',"{:0.4f}".format(toc-tic),'seconds.')

    return xx,yy,zz,outUnits

def ComputeLshell(inTime,c1,c2,c3,inCoordName,inCoordType):
    """
    This function uses spacepy to compute L shell and MLT from given time and position values.
    It uses the IRBEM library to do the calculations.
    
    INPUTS:
    inTime:       array:  time in UTC timestamp
    c1:           array:  x, alt
    c2:           array:  y, lat
    c3:           array:  z, lon
    inCoordName:  string: GDZ, GEO, GSM, GSE, SM, GEI, MAG, SPH, RLL
    inCoordType:  string: car, sph
    
    OUTPUT:
    c1:       array:  Lm
    c2:       array:  Lstar
    c3:       array:  MLT
    """

    # Start timer
    tic = time.perf_counter()
    #print("Starting ...")

    #pack position arrays in desired format
    sat_track = [[c1[i], c2[i], c3[i] ] for i in range(len(c1))]

    #convert utc_ts to spacepy time strings
    @np.vectorize
    def ts_to_spacepystr(ts):
        return dt.strftime(dt.utcfromtimestamp(int(ts)), '%Y-%m-%dT%H:%M:%S')
    
    #Steve Morley recommended converting to something other that a string to improve speed
    @np.vectorize
    def ts_to_spacepydt(ts):
        return dt.fromtimestamp(ts, tz=timezone.utc)

    #perform coordinate conversion
    cvals = Coords(sat_track, inCoordName, inCoordType)
    tt = Ticktock(ts_to_spacepydt(inTime), 'UTC')  #changed to UTC datetime object

    qq=get_Lstar(tt,cvals)
    Lm=abs(np.reshape(qq['Lm'],-1))
    Lstar=abs(np.reshape(qq['Lstar'],-1))
    MLT=qq['MLT']

    toc = time.perf_counter()
    print('Computed Lm,lstar,MLT in',"{:0.4f}".format(toc-tic),'seconds for',len(inTime),' positions.')

    return Lm,Lstar,MLT

