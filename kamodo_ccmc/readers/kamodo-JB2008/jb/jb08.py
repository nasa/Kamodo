# +
import sys
sys.path.append('/DESTOPy')

from python_utils.astro_function_jit import generateObservationsMEE, generateROMdensityModel, getDensityJB2008llajd, ep2pv_jit, UKF, readDTCFILE, from_jd, getTLEsForEstimation, propagateState_MeeBcRom
from python_utils.astro_function_jit import inputEOP_Celestrak_Full, readSOLFSMY, readDTCFILE
from python_utils.astro_function_jit import computeJB2000SWinputs, get_julian_datetime
from python_utils.JB2008_subfunc import JB2008
import spiceypy as spice

import numpy as np
import datetime
from kamodo import Kamodo, kamodofy


# -

class KJB08(Kamodo): 
    """Kamodofied JB2008 model"""
    def __init__(self, **kwargs):
        ## Load space weather data and Earth orientation parameters needed for JB2008 density model from file
        # Read Earth orientation parameters
        # need to specify the full path to these files
        spice.furnsh("/DESTOPy/Data/kernel.txt")
        
        self._eopdata = inputEOP_Celestrak_Full('/DESTOPy/Data/EOP-All.txt')

        # Read space weather data: solar activity indices
        self._SOLdata = readSOLFSMY('/DESTOPy/Data/SOLFSMY.txt')

        # Read geomagnetic storm DTC values
        self._DTCdata = readDTCFILE('/DESTOPy/Data/DTCFILE.txt')

    
        super(KJB08, self).__init__(**kwargs)
        
        @kamodofy(units='kg/m^3')
        def jb2008_density(lon, lat, alt, t):
            eop = self._eopdata  
            sol = self._SOLdata 
            dtc = self._DTCdata
            
            jd = get_julian_datetime(t)
            
            TEMP1, TEMP2, dens = self.getDensityJB2008llajd(lon, lat, alt, jd)
                        
            return dens
            
        self['rho'] = jb2008_density

    def getDensityJB2008llajd(self, lon, lat, alt, jdate):
        """
        Compute density of a particular point in atmosphere based on JB2008 density model and space weather input

        """
        eopdata = self._eopdata
        SOLdata = self._SOLdata
        DTCdata = self._DTCdata
        
        # Date and time
        dt_jdf = from_jd(np.ceil(round(jdate * 1e7,7))/ 1e7);  # End date of TLE collection window 
        yyUTC = dt_jdf.year
        mmUTC = dt_jdf.month
        ddUTC = dt_jdf.day
        hhUTC = dt_jdf.hour
        mnmnUTC = dt_jdf.minute
        ssUTC = dt_jdf.second

        doyUTC = dt_jdf.timetuple().tm_yday

        # Get JB2008 space weather data
        MJD,GWRAS,SUN,F10,F10B,S10,S10B,XM10,XM10B,Y10,Y10B,DSTDTC = computeJB2000SWinputs(yyUTC,doyUTC,hhUTC,mnmnUTC,ssUTC,SOLdata,DTCdata,eopdata,spice);

        SAT = np.zeros((3,1))
        XLON = np.deg2rad(lon); # Lon
        SAT[0] = np.mod(GWRAS + XLON, 2*np.pi);
        SAT[1] = np.deg2rad(lat); # Lat
        SAT[2] = alt;

        YRDAY = doyUTC + ((hhUTC*3600 + mnmnUTC*60 + ssUTC) / 86400)
        (TEMP1, TEMP2), rho = JB2008(MJD,YRDAY,SUN.flatten(),SAT.flatten(),F10,F10B,S10,S10B,XM10,XM10B,Y10,Y10B,DSTDTC)

        return TEMP1, TEMP2, rho





