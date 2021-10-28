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

def ConvertCoord(inTime,c1,c2,c3,inCoord,inType,outCoord,outType):
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
    units:    array:  [unit_c1, unit_c2, unit_c3]  (for example ['km','deg','deg'] or ['R_E','R_E','R_E'])
    """

    # Start timer
    tic = time.perf_counter()
    #print("Starting ...")

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

    #perform coordinate conversion
    cvals = Coords(sat_track, inCoord, inType)
    tt = Ticktock(ts_to_spacepystr(inTime), 'ISO')
    cvals.ticks = tt
    newC = cvals.convert(outCoord,outType)
    if outType == "sph":
        # Need to reorder Kamodo arrays to be consistent with spacepy
        zz, yy, xx = newC.data.T
    else:
        xx, yy, zz = newC.data.T

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
    print('Converted from ',inCoord,inType,'to:',outCoord,outType,outUnits,'in',"{:0.4f}".format(toc-tic),'seconds.')

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

    #perform coordinate conversion
    cvals = Coords(sat_track, inCoordName, inCoordType)
    tt = Ticktock(ts_to_spacepystr(inTime), 'ISO')

    qq=get_Lstar(tt,cvals)
    Lm=abs(np.reshape(qq['Lm'],-1))
    Lstar=abs(np.reshape(qq['Lstar'],-1))
    MLT=qq['MLT']

    toc = time.perf_counter()
    print('Computed Lm,lstar,MLT in',"{:0.4f}".format(toc-tic),'seconds for',len(inTime),' positions.')

    return Lm,Lstar,MLT

