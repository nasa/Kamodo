#!/usr/bin/env python
#-----------------------------------------------------------------------------
# $Id: gitm.py,v 1.16 2014/04/09 10:44:28 agburr Exp $
#
# GITM.py, Dan Welling, UMich
#
# Comments: Defines a class for GITM binary output files and modifies the input
#           to improve the data analysis and plotting experience
#
# Contains: class GitmBin    - The class for the GITM binary, which will read
#                              a single GITM output binary
#           def calc_magdi   - Reads a single GITM ion or mag output binary and
#                              computes the magnetic inclination and declination
#           def calc_magvel  - Reads a single GITM ion output binary and
#                              uses data from the standard GITM output (3DAll)
#                              to compute the ion characteristics in magnetic
#                              coordinates
#           def calc_deg     - Computes and appends latitude and longitude
#                              in degrees from values in radians
#           def calc_lt      - Computes and appends local time in hours from
#                              the universal time for the file and longitude
#           def append_units - Appends unit, descriptive name, and scale
#                              attributes to each known data type
#           def append_data  - Appends a list of data variables to a GitmBin
#                              data structure, where only a limited number of
#                              data variables from that file have been read
#                              in before
#           def calc_tec     - Calculate the VTEC
#           def calc_2dion   - Calculate the 2D ionospheric parameters (VTEC,
#                              hmF2, NmF2)
#
# Updates: Angeline Burrell (AGB) - 1/7/13: Added calc_lt, append_units, and
#                                           calc_magvel
#          AGB - 11/7/13: Improved calc_2dion, added Aaron Ridley's calc_tec
#          AGB - 12/6/13: Added Inclination/Declination calculation
#          Darren De Zeeuw (DDZ) - 06/24/19: Updated code to python3,
#                 Aaron Ridley approved reader for open source use in Kamodo
#------------------------------------------------------------------------------

'''
PyBats submodule for handling input/output for the Global 
Ionosphere-Thermosphere Model (GITM), one of the choices for the UA module
in the SWMF.
'''

# Global imports:
import numpy as np
import datetime as dt
from spacepy.pybats import PbData
from spacepy.datamodel import dmarray

class GitmBin(PbData):
    '''
    Object to open, manipulate and visualize 1-3 dimensional GITM output
    stored in binary format.  Object inherits from spacepy.pybats.PbData; see
    that documentation for general information on how these objects work.

    GITM index ordering is [lon, lat, altitude]; data arrays read from file
    will always be of the same shape and size.

    kwargs may be specified for:
    magfile = 3DION or 3DMAG file, allows computation of velocities in magnetic
              coordinates
    varlist = list of variable keys.  Will limit the variables read in to those
              listed.  Time and position will always be read in.  If the list
              is empty, all variables will be read in.
    '''

    def __init__(self, filename, *args, **kwargs):
        # Remove any known kwargs, as we don't want them to be included in
        # the GITM data keys!

        if 'varlist' not in kwargs:
            varlist = list()
        else:
            varlist = kwargs.pop('varlist')

        if 'magfile' not in kwargs:
            magfile = None
        else:
            magfile = kwargs.pop('magfile')

        # Initialize the GITM data structure
        super(GitmBin, self).__init__(*args, **kwargs) # Init as PbData.
        self.attrs['file']=filename

        # Load the GITM data
        self._read(varlist)
        self.calc_deg()
        self.calc_lt()

        self.append_units()

        if magfile or filename.find("3DION") >= 0:
            self.attrs['magfile']=magfile
            self.calc_magdi()
            self.calc_magvel()

    def __repr__(self):
        return 'GITM binary output file %s' % (self.attrs['file'])

    def _read(self, varlist, newfile=True):
        '''
        Read binary file.
        '''

        from re import sub
        from struct import unpack
        import sys
        
        # Read data and header info
        f=open(self.attrs['file'], 'rb')

        # Using the first FORTRAN header, determine endian.
        # Default is little.
        self.attrs['endian']='little'
        endChar='>'
        rawRecLen=f.read(4)
        if len(rawRecLen) < 4:
            print("GitmBin ERROR: empty file [", self.attrs['file'], "]")
            sys.exit(1)
        recLen=(unpack(endChar+'l',rawRecLen))[0]
        if (recLen>10000)or(recLen<0):
            # Ridiculous record length implies wrong endian.
            self.attrs['endian']='big'
            endChar='<'
            recLen=(unpack(endChar+'l',rawRecLen))[0]

        # Read version; read fortran footer+header.
        self.attrs['version']=unpack(endChar+'d',f.read(recLen))[0]
        (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

        # Read grid size information.
        (self.attrs['nLon'],self.attrs['nLat'],self.attrs['nAlt'])=\
            unpack(endChar+'lll',f.read(recLen))
        (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

        # Read number of variables.
        self.attrs['nVars']=unpack(endChar+'l',f.read(recLen))[0]
        (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

        # Collect variable names.
        var=[]
        for i in range(self.attrs['nVars']):
            var.append(unpack(endChar+'%is'%(recLen),f.read(recLen))[0])
            (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

        # Extract time. 
        (yy,mm,dd,hh,mn,ss,ms)=unpack(endChar+'lllllll',f.read(recLen))
        #print(yy,mm,dd,hh,mn,ss,ms)
        #self['time']=dt.datetime(yy,mm,dd,hh,mn,ss,ms/1000)
        self['time']=dt.datetime(yy,mm,dd,hh,mn,ss,ms)
        (oldLen)=unpack(endChar+'l',f.read(4))


        # Read the rest of the data.
        nTotal=self.attrs['nLon']*self.attrs['nLat']*self.attrs['nAlt']
        for val in var:
            # Trim variable names.
            v=sub('\[|\]', '', val.decode('utf-8')).strip()
            s=unpack(endChar+'l',f.read(4))[0]
            # Test to see if this variable is desired
            gvar=True
            if len(varlist) > 0:
                try:
                    varlist.index(v)
                except ValueError:
                    if((v.find('Altitude') < 0 and v.find('Longitude') < 0
                        and v.find('Latitude') < 0) or not newfile):
                        gvar=False
            # Unpack the data and save, if desired
            temp=unpack(endChar+'%id'%(nTotal),f.read(s))
            if gvar:
                self[v]=dmarray(np.array(temp))
                # Reshape arrays, note that ordering in file is Fortran-like.
                self[v]=self[v].reshape( 
                    (self.attrs['nLon'],self.attrs['nLat'],self.attrs['nAlt']),
                    order='F')
                
            f.read(4)


    def calc_deg(self):
        '''
        Gitm defaults to radians for lat and lon, which is sometimes difficult
        to use.  This routine leaves the existing latitude and longitude
        arrays intact and creates *dLat* and *dLon*, which contain the lat and
        lon in degrees.
        '''
        from numpy import pi
        import string

        self['dLat'] = dmarray(self['Latitude']*180.0/pi, 
                               attrs={'units':'degrees', 'scale':'linear',
                                      'name':'Latitude'})
        self['dLon'] = dmarray(self['Longitude']*180.0/pi, 
                               attrs={'units':'degrees', 'scale':'linear',
                                      'name':'Longitude'})

        # Do not correct for over/under limits because these are duplicates
        # that allow contour plots (and more) to display correctly
        #
        #for i in range(self.attrs['nLon']):
        #    for j in range(self.attrs['nLat']):
        #        if self['dLon'][i][j][0] < 0.0:
        #            self['dLon'][i][j] += 360.0
        #        elif self['dLon'][i][j][0] >= 360.0:
        #            self['dLon'][i][j] -= 360.0

    def calc_lt(self):
        '''
        Gitm defaults to universal time.  Compute local time from date and
        longitude.
        '''

        from numpy import pi
        import math

        ut = (self['time'].hour * 3600 + self['time'].minute * 60
              + self['time'].second + self['time'].microsecond * 1e-6) / 3600.0
        self['LT'] = dmarray(ut + self['Longitude']*12.0/pi,
                             attrs={'units':'hours', 'scale':'linear',
                                    'name':'Local Time'})

        # Because datetime won't do lists or arrays
        if dmarray.max(self['LT']) >= 24.0:
            for i in range(self.attrs['nLon']):
                # All local times will be the same for each lat/alt
                # in the same longitude index
                ltmax = dmarray.max(self['LT'][i])
                if ltmax >= 24.0:
                    self['LT'][i] -= 24.0 * math.floor(ltmax / 24.0)
    def calc_magdi(self):
        '''
        GITM 3DION and 3DMAG files contain the magnetic field in
        North-East-Vertical coordinates.  This routine computes the magnetic
        inclination and declination.
        '''
        import math

        mag = None
        # Test to determine the type of file we have
        if(self.attrs['file'].find("ION") > 0
           or self.attrs['file'].find("MAG") > 0):
            mag = self
        else:
            if 'magfile' not in self.attrs:
                print("No 3D MAG/ION file associated with this GITM Binary")
            elif(self.attrs['magfile'].find("ION") <= 0 and
                 self.attrs['magfile'].find("MAG") <= 0):
                print("No 3D MAG/ION file associated with this GITM Binary")
            else:
                mag = GitmBin(self.attrs['magfile'])

        if mag is not None:
            dec_frac = mag['B.F. East'] / mag['B.F. North']
            inc_sign = -1.0 * mag['B.F. Vertical'] / abs(mag['B.F. Vertical'])
            inc_pmag = mag['B.F. North']**2 + mag['B.F. East']**2

            for ilon in range(self.attrs['nLon']):
                for ilat in range(self.attrs['nLat']):
                    for ialt,df in enumerate(dec_frac[ilon,ilat]):
                        dec_frac[ilon,ilat,ialt] = math.atan(df) * 180.0 / np.pi
                        i = (inc_sign[ilon,ilat,ialt]
                             * math.sqrt(inc_pmag[ilon,ilat,ialt])
                             / mag['B.F. Magnitude'][ilon,ilat,ialt])
                        inc_pmag[ilon,ilat,ialt] = math.acos(i) * 180.0 / np.pi
                        if inc_pmag[ilon,ilat,ialt] > 90.0:
                            inc_pmag[ilon,ilat,ialt] -= 180.0

            self['Declination'] = dmarray.copy(dec_frac)
            self['Declination'].attrs = {"name":"Declination", "scale":"linear",
                                         "units":"degrees"}
            self['Inclination'] = dmarray.copy(inc_pmag)
            self['Inclination'].attrs = {"name":"Inclination", "scale":"linear",
                                         "units":"degrees"}
            del dec_frac, inc_sign, inc_pmag

            if 'B.F. East' not in self:
                self['B.F. East'] = dmarray.copy(mag['B.F. East'])
                self['B.F. North'] = dmarray.copy(mag['B.F. North'])
                self['B.F. Vertical'] = dmarray.copy(mag['B.F. Vertical'])
                self['B.F. Magnitude'] = dmarray.copy(mag['B.F. Magnitude'])
            if 'Magnetic Latitude' not in self:
                self['Magnetic Latitude']=dmarray.copy(mag['Magnetic Latitude'])
                self['Magnetic Longitude']=dmarray.copy(mag['Magnetic Longitude'])

    def calc_magvel(self):
        '''
        Gitm 3DIon files contain the magnetic coordinates that allow the
        field-aligned and field-perpendicular velocities to be computed.
        The 3DIon file must be from the same run as the 3DAll file so that
        the geographic grid is the same.  If a 3D Ion file was not produced
        in the original run, don't fret!  You can get one by using the same
        UAM.in file.  Unless the magnetic field is varying secularly, you
        can get away with producing just one 3DIon file.  It is better to have
        a matching 3D Ion file for every 3D All file, however.
        '''
        import math
        import string
        import sys

        ion = None
        # If this is a 3D ION file we don't need a mag file
        if self.attrs['file'].find("ION") > 0:
            ion = self
        else:
            if 'magfile' not in self.attrs:
                print("No 3D MAG/ION file associated with this GITM Binary")
            elif(self.attrs['magfile'].find("ION") <= 0 and
                 self.attrs['magfile'].find("MAG") <= 0):
                print("No 3D MAG/ION file associated with this GITM Binary")
            else:
                ion = GitmBin(self.attrs['magfile'])

        # Compute the field-aligned unit vector in East, North,
        # and Vertical coordinates
        if ion:
            bhat_e = ion['B.F. East'] / ion['B.F. Magnitude']
            bhat_n = ion['B.F. North'] / ion['B.F. Magnitude']
            bhat_v = ion['B.F. Vertical'] / ion['B.F. Magnitude']

            # Compute the zonal unit vector in East, North, Vertical coord

            mag = np.sqrt(np.square(ion['B.F. East'])
                          + np.square(ion['B.F. North']))

            zhat_e = ion['B.F. North'] / mag
            zhat_n = ion['B.F. East'] / mag
            # zhat_v is identically zero

            # Compute the meridional unit vector in East, North, Vertical coord

            mhat_e = (-ion['B.F. East']*ion['B.F. Vertical']
                       / (mag * ion['B.F. Magnitude']))
            mhat_n = (-ion['B.F. North']*ion['B.F. Vertical']
                       / (mag * ion['B.F. Magnitude']))
            mhat_v = mag / ion['B.F. Magnitude']

            # Compute the magnetic coordinate velocities for each overlapping
            # latitude, longitude, and altitude.  Also include the mag coord.

            for item in list(self.keys()):
                if(string.find(item, "V!") >= 0 or
                   string.find(item, "Gravity") >= 0 or
                   string.find(item, "PressGrad") >= 0):
                    sp = string.split(item)
                    units = self[item].attrs['units']
                    scale = self[item].attrs['scale']
                    # extract the non-directional part of the name
                    if sp[0] not in self:
                        east = string.join([sp[0], "(east)"], " ")
                        north = string.join([sp[0], "(north)"], " ")
                        up = string.join([sp[0], "(up)"], " ")
                        par = string.join([sp[0], "(par)"], " ")
                        zon = string.join([sp[0], "(zon)"], " ")
                        mer = string.join([sp[0], "(mer)"], " ")
                        name = self[east].attrs['name']

                        vp=bhat_e*self[east]+bhat_n*self[north]+bhat_v*self[up]
                        vz=zhat_e*self[east]+zhat_n*self[north]
                        vm=mhat_e*self[east]+mhat_n*self[north]+mhat_v*self[up]
                    elif sp[0].find("Gravity") >= 0:
                        par = string.join([sp[0], "(par)"], " ")
                        zon = string.join([sp[0], "(zon)"], " ")
                        mer = string.join([sp[0], "(mer)"], " ")
                        name = "%s$_{east}$" % (self[item].attrs['name'])

                        vp=bhat_v*self[item]
                        vz=0.0*self[item]
                        vm=mhat_v*self[item]

                    self[par] = dmarray.copy(vp)
                    self[par].attrs = {'units':'%s{\mathdefault{,\,positive\,mag\,north}}'%(units), 'scale':scale, 'name':name.replace("east", "\parallel")}
                    self[zon] = dmarray.copy(vz)
                    self[zon].attrs = {'units':'%s{\mathdefault{,\,positive\,east}}'%(units), 'scale':scale, 'name':name.replace("east", "zon")}
                    self[mer] = dmarray.copy(vm)
                    self[mer].attrs = {'units':'%s{\mathdefault{,\,positive\,up}}'%(units), 'scale':scale, 'name':name.replace("east", "mer")}
        
            if 'B.F. East' not in self:
                self['B.F. East'] = dmarray.copy(ion['B.F. East'])
                self['B.F. North'] = dmarray.copy(ion['B.F. North'])
                self['B.F. Vertical'] = dmarray.copy(ion['B.F. Vertical'])
                self['B.F. Magnitude'] = dmarray.copy(ion['B.F. Magnitude'])
            if 'Magnetic Latitude' not in self:
                self['Magnetic Latitude'] = dmarray.copy(ion['Magnetic Latitude'])
                self['Magnetic Longitude'] = dmarray.copy(ion['Magnetic Longitude'])

    def append_units(self):
        '''
        Append units, descriptive names, and plot scaling (e.g. linear,
        exponential) to the attributes of known data types
        '''
        unit_dict = {"Altitude":"m", "Ar Mixing Ratio":"", "Ar":"kg \, m^{-3}",
                     "CH4 Mixing Ratio":"", "Conduction":"W m$^{-1}$ K$^{-1}$",
                     "EuvHeating":"K per timestep", "H":"kg \, m^{-3}",
                     "H!U+!N":"kg \, m^{-3}", "H2 Mixing Ratio":"",
                     "HCN Mixing Ratio":"", "He":"kg \, m^{-3}",
                     "He!U+!N":"kg \, m^{-3}", "Heating Efficiency":"",
                     "Heat Balance Total":"", "Latitude":"radians",
                     "Longitude":"radians", 
                     "N!D2!N":"kg \, m^{-3}",
                     "N!D2!U+!N":"kg \, m^{-3}",
                     "N!D2!U+!N            (/m3)":"m^{-3}",
                     "N!U+!N":"kg \, m^{-3}",
                     "N(!U2!ND)":"kg \, m^{-3}",
                     "N(!U2!NP)":"kg \, m^{-3}",
                     "N(!U4!NS)":"kg \, m^{-3}",
                     "N2 Mixing Ratio":"",
                     "NO":"kg \, m^{-3}",
                     "NO!U+!N":"kg \, m^{-3}",
                     "O!D2!N":"kg \, m^{-3}", "O(!U1!ND)":"kg \, m^{-3}",
                     "O!D2!U+!N":"kg \, m^{-3}", "O(!U2!ND)!":"m^{-3}",
                     "O(!U2!ND)!U+!N":"m^{-3}", 
                     "O(!U2!NP)!U+!N":"kg \, m^{-3}",
                     "O(!U2!NP)!U+!N":"kg \, m^{-3}",
                     "O(!U3!NP)":"kg \, m^{-3}",
                     "O_4SP_!U+!N":"kg \, m^{-3}",
                     "RadCooling":"", 
                     "Rho":"kg \, m^{-3}", 
                     "Temperature":"K",
                     "V!Di!N (east)":"m s^{-1}",
                     "Vi (east) (m/s)":"m s^{-1}",
                     "V!Di!N (north)":"m s^{-1}",
                     "Vi (north) (m/s)":"m s^{-1}",
                     "V!Di!N (up)":"m s^{-1}",
                     "Vi (up) (m/s)": "m s^{-1}",
                     "V!Dn!N (east)":"m s^{-1}",
                     "Vn (east) (m/s)":"m s^{-1}",
                     "V!Dn!N (north)":"m s^{-1}",
                     "Vn (north) (m/s)":"m s^{-1}", 
                     "V!Dn!N (up)":"m s^{-1}",
                     "Vn (up) (m/s)":"m s^{-1}",
                     "V!Dn!N (up,N!D2!N              )":"m s^{-1}",
                     "V!Dn!N (up,N(!U4!NS)           )":"m s^{-1}",
                     "V!Dn!N (up,NO                  )":"m s^{-1}",
                     "V!Dn!N (up,O!D2!N              )":"m s^{-1}",
                     "V!Dn!N (up,O(!U3!NP)           )":"m s^{-1}",
                     "V!Dn!N (up,He                  )":"m s^{-1}",
                     "e-":"m^{-3}", 
                     "Electron_Average_Energy":"J",
                     "eTemperature":"K", "iTemperature":"K", "LT":"h",
                     "Solar Zenith Angle":"radians", 
                     "Vertical TEC":"TECU",
                     "CO!D2!N":"kg \, m^{-3}",  
                     "DivJu FL":"", "DivJuAlt":"",
                     "Electron_Energy_Flux":"J m$^{-2}$", "FL Length":"m",
                     "Pedersen FL Conductance":"S m^{-1}", "dLon":"degrees",
                     "Pedersen Conductance":"S m^{-1}", "dLat":"degrees",
                     "Hall FL Conductance":"S m^{-1}", "Potential":"V",
                     "Hall Conductance":"S m^{-1}", "Je2":"A m^{-2}",
                     "Je1":"A m^{-2}", "Magnetic Longitude":"degrees",
                     "E.F. Vertical":"V m^{-1}", "E.F. East":"V m^{-1}",
                     "E.F. North":"V m^{-1}", "E.F. Magnitude":"V m^{-1}",
                     "B.F. Vertical":"nT", "B.F. East":"nT", "B.F. North":"nT",
                     "B.F. Magnitude":"nT", "Magnetic Latitude":"degrees",
                     "Ed1":"", "Ed2":"", "Gravity":"m s^{-2}",
                     "PressGrad (east)":"Pa\;m^{-1}", "Joule Heating":"K per timestep",
                     "PressGrad (north)":"Pa\;m^{-1}", "O Cooling":"K per timestep",
                     "PressGrad (up)":"Pa\;m^{-1}", "Total Abs EUV":"K per timestep",
                     "IN Collision Freq":"s^{-1}", "Chemical Heating":"",
                     "Auroral Heating":"K per timestep", "Photoelectron Heating":"K per timestep",
                     "Eddy Conduction":"", "Eddy Adiabatic Conduction":"",
                     "NO Cooling":"K per timestep", "Molecular Conduction":""}

        scale_dict = {"Altitude":"linear", "Ar Mixing Ratio":"linear",
                      "Ar":"exponential", "CH4 Mixing Ratio":"linear",
                      "Conduction":"linear", "EuvHeating":"linear",
                      "H":"exponential", "H!U+!N":"exponential",
                      "H2 Mixing Ratio":"linear", "HCN Mixing Ratio":"linear",
                      "He":"exponential", "He!U+!N":"exponential",
                      "Heating Efficiency":"linear", "DivJuAlt":"linear",
                      "Heat Balance Total":"linear", "Latitude":"linear",
                      "Longitude":"linear", 
                      "N!D2!N":"exponential",
                      "N!D2!U+!N":"exponential",
                      "N!U+!N":"exponential",
                      "N(!U2!ND)":"exponential",
                      "N(!U2!NP)":"exponential",
                      "N(!U4!NS)":"exponential",
                      "N2 Mixing Ratio":"linear",
                      "NO":"exponential", 
                      "NO!U+!N":"exponential",
                      "O!D2!N":"exponential", 
                      "O!D2!U+!N":"exponential",
                      "O(!U2!ND)!":"exponential", "O(!U1!ND)":"exponential",
                      "O(!U2!ND)!U+!N":"exponential", 
                      "CO!D2!N":"exponential",
                      "O(!U2!NP)!U+!N":"exponential", "DivJu FL":"",
                      "O(!U2!NP)!U+!N":"exponential", 
                      "O(!U3!NP)":"exponential",
                      "O_4SP_!U+!N":"exponential",
                      "RadCooling":"linear",
                      "Rho":"exponential", "Temperature":"linear",
                      "V!Di!N (east)":"linear", 
                      "Vi (east)":"linear", 
                      "V!Di!N (north)":"linear",
                      "Vi (north) (m/s)":"linear",
                      "V!Di!N (up)":"linear", 
                      "Vi (up) (m/s)":"linear",
                      "V!Dn!N (east)":"linear",
                      "Vn (east) (m/s)":"linear",
                      "V!Dn!N (north)":"linear",
                      "Vn (north) (m/s)":"linear",
                      "V!Dn!N (up)":"linear",
                      "Vn (up) (m/s)":"linear",
                      "V!Dn!N (up,N!D2!N              )":"linear",
                      "V!Dn!N (up,N(!U4!NS)           )":"linear",
                      "V!Dn!N (up,NO                  )":"linear",
                      "V!Dn!N (up,O!D2!N              )":"linear",
                      "V!Dn!N (up,O(!U3!NP)           )":"linear",
                      "V!Dn!N (up,He                  )":"linear",
                      "e-":"linear", 
                      # "e-                   (/m3)":"linear",
                      "Electron_Average_Energy":"linear",
                      "eTemperature":"linear", "iTemperature":"linear",
                      "Solar Zenith Angle":"linear", "Vertical TEC":"linear",
                      "Electron_Energy_Flux":"exponential",
                      "FL Length":"linear", "Pedersen FL Conductance":"linear",
                      "Hall Conductance":"linear", "Potential":"linear",
                      "Hall FL Conductance":"linear", "dLon":"linear",
                      "Pedersen Conductance":"linear", "Je2":"linear",
                      "Je1":"linear", "Ed1":"linear", "Ed2":"linear",
                      "E.F. Vertical":"linear", "E.F. East":"linear",
                      "E.F. North":"linear", "E.F. Magnitude":"linear",
                      "B.F. Vertical":"linear", "B.F. East":"linear",
                      "B.F. North":"linear", "B.F. Magnitude":"linear",
                      "Magnetic Latitude":"linear", "LT":"linear",
                      "Magnetic Longitude":"linear", "dLat":"linear",
                      "Gravity":"linear", "PressGrad (east)":"linear",
                      "PressGrad (north)":"linear", "PressGrad (up)":"linear",
                      "IN Collision Freq":"linear", "Chemical Heating":"linear",
                      "Total Abs EUV":"linear", "O Cooling":"linear",
                      "Joule Heating":"linear", "Auroral Heating":"linear",
                      "Photoelectron Heating":"linear", "NO Cooling":"linear",
                      "Eddy Conduction":"linear",
                      "Molecular Conduction":"linear",
                      "Eddy Adiabatic Conduction":"linear"}

        name_dict = {"Altitude":"Altitude",
                     "Ar Mixing Ratio":"Argon Mixing Ratio",
                     "Ar":"Ar Mass Density",
                     "CH4 Mixing Ratio":"Methane Mixing Ratio",
                     "Conduction":"Conduction", "EuvHeating":"EUV Heating",
                     "H":"H Mass Density", "H!U+!N":"H$^+$ Mass Density",
                     "H2 Mixing Ratio":"H$_2$ Mixing Ratio",
                     "HCN Mixing Ratio":"Hydrogen Cyanide Mixing Ratio",
                     "He":"He Mass Density", "He!U+!N":"He$^+$ Mass Density",
                     "Heating Efficiency":"Heating Efficiency",
                     "Heat Balance Total":"Heat Balance Total",
                     "Latitude":"Latitude", "Longitude":"Longitude",
                     "N!D2!N":"N$_2$ Mass Density",
                     "N!D2!U+!N":"N$_2$$^+$ Mass Density",
                     "N!U+!N":"N$^+$ Mass Density",
                     "N(!U2!ND)":"N($^2$D) Mass Density",
                     "N(!U2!NP)":"N($^2$P) Mass Density",
                     "N(!U4!NS)":"N($^4$S) Mass Density",
                     "N2 Mixing Ratio":"N$_2$ Mixing Ratio",
                     "NO":"NO Mass Density", 
                     "NO!U+!N":"NO$^+$ Mass Density",
                     "O!D2!N":"O$_2$ Mass Density",
                     "O(!U1!ND)":"O($^1$D) Mass Density",
                     "O!D2!U+!N":"O$_2$$^+$ Mass Density",
                     "O(!U2!ND)!":"O($^2$D) Mass Density",
                     "O(!U2!ND)!U+!N":"O($^2$D) Mass Density",
                     "O(!U2!NP)!U+!N":"O($^2$P)$^+$ Mass Density",
                     "O(!U2!NP)!U+!N":"O($^2$P) Mass Density",
                     "O(!U3!NP)":"O($^3$P) Mass Density",
                     "O_4SP_!U+!N":"O($^4$SP)$^+$ Mass Density",
                     "RadCooling":"Radiative Cooling", 
                     "Rho":"Neutral Density",
                     "Temperature":"T$_n$", "V!Di!N (east)":"v$_{east}$",
                     "V!Di!N (north)":"v$_{north}$", 
                     "Vi (north) (m/s)":"v$_{north}$", 
                     "V!Di!N (up)":"v$_{up}$",
                     "Vi (up) (m/s)": "v$_{up}$",
                     "V!Dn!N (east)":"u$_{east}$",
                     "Vn (east) (m/s)":"u$_{east}$",
                     "V!Dn!N (north)":"u$_{north}$",
                     "Vn (north) (m/s)":"u$_{north}$", 
                     "V!Dn!N (up)":"u$_{up}$",
                     "Vn (up) (m/s)":"u$_{up}$",
                     "V!Dn!N (up,N!D2!N              )":"u$_{Up, N_2}$",
                     "V!Dn!N (up,N(!U4!NS)           )":"u$_{Up, N(^4S)}$",
                     "V!Dn!N (up,NO                  )":"u$_{Up, NO}$",
                     "V!Dn!N (up,O!D2!N              )":"u$_{Up, O_2}$",
                     "V!Dn!N (up,O(!U3!NP)           )":"u$_{Up, O(^3P)}$",
                     "V!Dn!N (up,He                  )":"u$_{Up, He}$",
                     "e-":"[e-]",
                     "Electron_Average_Energy":"Electron Average Energy",
                     "eTemperature":"T$_e$", "iTemperature":"T$_i$",
                     "Solar Zenith Angle":"Solar Zenith Angle",
                     "Vertical TEC":"VTEC", "CO!D2!N":"CO$_2$ Mass Density",
                     "DivJu FL":"DivJu FL", "DivJuAlt":"DivJuAlt",
                     "Electron_Energy_Flux":"Electron Energy Flux",
                     "FL Length":"Field Line Length",
                     "Pedersen FL Conductance":"$\sigma_P$",
                     "Pedersen Conductance":"$\Sigma_P$",
                     "Hall FL Conductance":"$\sigma_H$",
                     "Potential":"Potential", "Hall Conductance":"$\Sigma_H$",
                     "Je2":"Region 2 Current", "Je1":"Region 1 Current",
                     "Ed1":"Ed1", "Ed2":"Ed2", "LT":"Solar Local Time",
                     "E.F. Vertical":"Vertical Electric Field",
                     "E.F. East":"Eastward Electric Field",
                     "E.F. North":"Northward Electric Field",
                     "E.F. Magnitude":"Electric Field Magnitude",
                     "B.F. Vertical":"Vertical Magnetic Field",
                     "B.F. East":"Eastward Magnetic Field",
                     "B.F. North":"Northward Magnetic Field",
                     "B.F. Magnitude":"Magnetic Field Magnitude",
                     "Magnetic Latitude":"Magnetic Latitude",
                     "Magnetic Longitude":"Magnetic Longitude",
                     "dLat":"Latitude", "dLon":"Longitude", "Gravity":"g",
                     "PressGrad (east)":r"$\nabla_{east}$ (P$_i$ + P$_e$)",
                     "PressGrad (north)":r"$\nabla_{north}$ (P$_i$ + P$_e$)",
                     "PressGrad (up)":r"$\nabla_{up}$ (P$_i$ + P$_e$)",
                     "IN Collision Freq":r"$\nu_{in}$",
                     "Chemical Heating":"Chemical Heating Rate",
                     "Total Abs EUV":"Total Absolute EUV",
                     "O Cooling":"O Cooling", "Joule Heating":"Joule Heating",
                     "Auroral Heating":"Auroral Heating",
                     "Photoelectron Heating":"Photoelectron Heating",
                     "Eddy Conduction":"Eddy Conduction",
                     "Eddy Adiabatic Conduction":"Adiabatic Eddy Conduction",
                     "NO Cooling":"NO Cooling",
                     "Molecular Conduction":"Molecular Conduction"}

        for k in list(self.keys()):
            if type(self[k]) is dmarray:
                nk = k
                # Temporary fix for misspelled key (9/30/13)
                if nk.find("Heaing Efficiency") >= 0:
                    nk = "Heating Efficiency"
                elif nk.find("EUV Heating") >= 0:
                    nk = "EuvHeating"
                elif nk.find("Rho (kg/m3)") >= 0:
                    nk = "Rho"
                elif nk.find("Neutral Temperature (K)") >= 0:
                    nk = 'Temperature'
                elif nk.find("Vn (up) (m/s)") >= 0:
                    nk = 'V!Dn!N (up)'
                elif nk.find("Vi (east) (m/s)") >= 0:
                    nk = "V!Di!N (east)"

                try:
                    #print ('registering',k)
                    self.register_name(k, nk, unit_dict, scale_dict, name_dict)
                except:
                    if k.split()[0] in name_dict:
                        nk = k.split()[0]
                        if "(/m3)" in k:
                            # print 'found (/m3) in key'
                            unit_dict[k] = "m^{-3}"
                            name_dict[k] = name_dict[nk]
                            scale_dict[k] = scale_dict[nk]
                            self.register_name(k, k, unit_dict, scale_dict, name_dict)
                        else:
                            raise
                    else:
                        raise
            
    def register_name(self, k, nk, unit_dict, scale_dict, name_dict):
        # Different versions of GITM differ in header capitalization
        if nk not in name_dict:
            # Try to capitalize or lowercase the key
            if nk == nk.capitalize():
                nk = k.lower()
            else:
                nk = k.capitalize()
        if 'units' not in self[k].attrs:
            self[k].attrs['units'] = unit_dict[nk]
        if 'scale' not in self[k].attrs:
            self[k].attrs['scale'] = scale_dict[nk]
        if 'name' not in self[k].attrs:
            self[k].attrs['name'] = name_dict[nk]

    def append_data(self, varlist):
        '''
        A routine to append more variables to an existing GITM data structure.
        New variables can only be added from the same file (specified in
        the 'file' attribute).
        '''

        temp = GitmBin(self.attrs['file'], varlist, False)

        for nkey in list(temp.keys()):
            if nkey not in self:
                self[nkey] = dmarray.copy(temp[nkey])

    def calc_tec(self):
        '''
        A routine to calculate the 2D VTEC.
        To perform these calculations, electron density ("e-") must be one of
        the available data types.
        '''
        import scipy.integrate as integ
        from scipy.interpolate import interp1d

        if 'e-' in self:
            self['VTEC'] = dmarray(self['e-'] * 1.0e-16,
                                   attrs={'units':'TECU', 'scale':'linear',
                                          'name':'Vertical TEC'})

            for ilon in range(self.attrs['nLon']):
                for ilat in range(self.attrs['nLat']):
                    # Integrate electron density over altitude, not including
                    # ghost cells
                    vtec = integ.simps(self['VTEC'][ilon,ilat,2:-2],
                                       self['Altitude'][ilon,ilat,2:-2], "avg")
                    self['VTEC'][ilon,ilat,:] = vtec

    def calc_2dion(self):
        '''
        A routine to calculate the 2D ionospheric parameters: VTEC, hmF2, NmF2.
        To perform these calculations, electron density ("e-") must be one of
        the available data types.
        '''
        import scipy.integrate as integ
        from scipy.interpolate import interp1d
        from scipy.signal import argrelextrema

        calc_tec = False
        if 'e-' in self:
            if ('VTEC' in self) is False:
                calc_tec = True
                self['VTEC'] = dmarray(self['e-'] * 1.0e-16,
                                       attrs={'units':'TECU', 'scale':'linear',
                                              'name':'Vertical TEC'})
            self['NmF2'] = dmarray.copy(self['e-'])
            self['NmF2'].attrs['name'] = 'N$_m$F$_2$'
            self['hmF2'] = dmarray(self['Altitude'] / 1000.0,
                                   attrs={'units':'km', 'scale':'linear',
                                          'name':'h$_m$F$_2$'})
            alt = np.linspace(min(self['hmF2'][0,0,2:-2]),
                              max(self['hmF2'][0,0,2:-2]),1000)

            for ilon in range(self.attrs['nLon']):
                for ilat in range(self.attrs['nLat']):
                    if calc_tec is True:
                        # Integrate electron density over altitude
                        vtec = integ.simps(self['VTEC'][ilon,ilat,2:-2],
                                           self['Altitude'][ilon,ilat,2:-2],
                                           "avg")
                        self['VTEC'][ilon,ilat,:] = vtec

                    # Interpolate over the electron density altitude profile
                    eprof = interp1d(self['hmF2'][ilon,ilat,2:-2],
                                     self['NmF2'][ilon,ilat,2:-2], kind="cubic")
                    try:
                        edens = eprof(alt)
                    except:
                        alt = np.linspace(min(self['hmF2'][ilon,ilat,2:-2]),
                                          max(self['hmF2'][ilon,ilat,2:-2]),
                                          1000)
                        edens = eprof(alt)

                    # Find the local maxima of the electron density profile
                    emax = argrelextrema(edens, np.greater)
                    emax_list = list(emax[0])
                    saddle = False
                    if len(emax_list) == 0:
                        # No local maxima were found, try identifying
                        # saddle points
                        saddle = True
                    elif len(emax_list) == 1:
                        if max(edens) == edens[emax_list[0]]:
                            # Only one maxima exists and is realistic. Sometimes
                            # a profile will have inflection points instead of
                            # local maxima and this can confuse the routine
                            self['NmF2'][ilon,ilat,:] = edens[emax_list[0]]
                            self['hmF2'][ilon,ilat,:] = alt[emax_list[0]]
                        else:
                            saddle = True
                    elif alt[emax_list[-1]] < 120.0:
                        saddle = True
                    else:
                        # More than one maxima exists.  Seperate hmF2 from hmF1
                        # and spurious local maxima
                        NmF2 = list(edens[emax_list])
                        HmF2 = list(alt[emax_list])

                        # If the global maximum is over 200 km,
                        # this is the F2 peak
                        eindex = NmF2.index(max(NmF2))
                        if HmF2[eindex] <= 200.0 and HmF2[eindex] == min(HmF2):
                            # The global max may be the F1 peak, see if the
                            # secondary (or lessor) maxima is at the upper
                            # limit of the model.  If so, remove this point
                            # from consideration
                            if max(HmF2) > self['hmF2'][ilon,ilat,-5]:
                                eindex = HmF2.index(max(HmF2))
                                emax_list.pop(eindex)
                                NmF2.pop(eindex)
                                HmF2.pop(eindex)
                                eindex = NmF2.index(max(NmF2))

                            if len(emax_list) > 1:
                                # If there are multiple maxima after the upper
                                # boundary has been removed, choose the largest
                                # maxima above 200 km since the hmF1 is often
                                # larger than the hmF2
                                emax_list.pop(eindex)
                                NmF2.pop(eindex)
                                eindex = NmF2.index(max(NmF2))

                        # Set the hmF2 and NmF2
                        self['NmF2'][ilon,ilat,:] = edens[emax_list[eindex]]
                        self['hmF2'][ilon,ilat,:] = alt[emax_list[eindex]]

                    if saddle:
                        # It is difficult to find saddle points.  Examine
                        # the rate of change of density
                        delta_alt = alt[1] - alt[0] # Equally spaced
                        edens_dot = np.diff(edens) / delta_alt

                        # Identify inflection points by looking for
                        # minima in the derivative of electron density.
                        edot_min = argrelextrema(edens_dot, np.less)
                        emin_list = list(edot_min[0])
                        edens_min = list(edens[emin_list])
                        emax = np.max(edens_min)
                        eindex = emin_list[edens_min.index(emax)]

                        # Assign the inflection with the largest
                        # electron density to the ion peak
                        self['NmF2'][ilon,ilat,:] = edens[eindex]
                        self['hmF2'][ilon,ilat,:] = alt[eindex]

    def lon_lt_ticks(self, x, pos):
        '''
        Define ticks to include local time in addition to longitude.  Assumes
        that x is longitude in degrees.
        '''
        import math
        from . import gitm_plot_rout as gpr

        # Calculate the local time in hours
        lth = gpr.glon_to_localtime(self['time'], x)
        ltm = int((lth - math.floor(lth)) * 60.0)

        # Build the format string
        fmtstring = "{:g}$^\circ$\n {:02d}:{:02d}".format(x, int(lth), ltm)

        return(fmtstring)

# END

