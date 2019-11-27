#!/usr/bin/env python
#-----------------------------------------------------------------------------
# $Id: gitm_plot_rout.py,v 1.18 2014/02/17 16:41:48 agburr Exp $
# gitm_plot_rout
#
# Author: Angeline G. Burrell, UMichigan, Jan 2013
#
# Comments: Common routine used to make GITM plots.
#
# AGB: 10/17/13: Added several routines and improved colorbar formatting
# Darren De Zeeuw (DDZ) - 06/24/19: Updated code to python3, 
#                 Aaron Ridley approved reader for open source use in Kamodo
#
# Includes: choose_contour_map         - choose a color map depending on certain
#                                        specified plot characteristics
#           add_colorbar               - add a colorbar to a contour plot
#           find_order_of_magnitude    - find the order of magnitude
#           center_polar_cap           - center radial coordinates for a
#                                        polar plot
#           find_data_limits           - find the upper and lower limits
#                                        for a list of GITM data arrays
#           find_data_limits_irange    - find the upper and lower limits
#                                        for a list of GITM data arrays and
#                                        a range of lon/lat/alt indexes
#           find_data_limits_ivalues   - find the upper and lower limits
#                                        for a list of GITM data arrays and
#                                        a range or a specific lon/lat/alt index
#           glon_to_localtime          - convert longitude to local time
#           localtime_to_glon          - convert local time to longitude
#           find_lon_lat_index         - find the appropriate index for
#                                        a specified location
#           retrieve_key_from_web_name - a routine to retrieve the data key
#                                        from a website-friendly data name
#           find_alt_index             - A routine to find the appropriate index
#                                        for a specified altitude
#           match_cindi_key            - a routine to retrieve the CINDI data
#                                        key from a GITM key
#           add_subsolar_point         - calc and add subsolar point to a plot
#           add_geomagnetic_equator    - a routine to add the geomagnetic
#                                        equator to an existing plot
#           add_solar_terminator       - a routine to calculate the lat/lon
#                                        of the solar terminator at a given
#                                        altitude.  Can add this to a plot.
#           create_contour_input_array - Create a contour input 2D numpy array
#                                        at a specified geo/mag lat/lon/alt
#           create_linear_input_array  - interpolate linear data at a location
#           get_meq_offset             - Calculate the latitude offset between
#                                        the geographic and geomagnetic equators
#           add_dipole_fieldline       - Calculate and add dipole field lines
#                                        with a specified apex height to a plot
#           
#----------------------------------------------------------------------------

'''
Plot data from a 2D GITM file (or 3D GITM file at a single altitude) for
different geographic configurations
'''

# Import modules
import numpy as np

def choose_contour_map(color, center):
    '''
    Choose a standard contour color map based on whether the output image
    will be black and white or color, centered around zero or not.

    Input: color  = True for color, False for black and white
           center = True for centered about zero, False if not
    '''

    if color:
        if center:
            return("seismic_r")
        else:
            return("Spectral_r")
    else:
        if center:
            return("binary")
        else:
            return("Greys")

def add_colorbar(contour_handle, zmin, zmax, zinc, orient, scale, name, units):
    '''
    Add a colorbar

    Input: contour_handle = handle to contour plot
           zmin           = minimum z value
           zmax           = maximum z value
           zinc           = z tick incriment (recommend 6)
           orient         = orientation of the colorbar (horizontal or vertical)
           scale          = linear or exponential?
           name           = z variable name
           units          = z variable units
    '''
    import matplotlib.pyplot as plt
    import math
    from matplotlib.ticker import FormatStrFormatter, FuncFormatter

    if scale.find("exponential") >= 0:
        magmin = find_order_of_magnitude(zmin)
        magmax = find_order_of_magnitude(zmax)
        expinc = abs(magmin - magmax) + 1
        w = np.logspace(magmin, magmax, expinc)
    else:
        w = np.linspace(zmin, zmax, zinc, endpoint=True)

    cb = plt.colorbar(contour_handle, ticks=w, pad=.15, orientation=orient,
                      fraction=.07)
    zscale = max(abs(zmin), abs(zmax))

    if scale.find("exponential") >= 0:
        def exp_ticks(x, pos):
            '''
            Define ticks so that they use scientific notation
            '''
            omag = find_order_of_magnitude(x)
            lexp = True
            rord = int(expinc / zinc)

            # If there are too many ticks to label, label half of them
            if expinc >= 2.0 * zinc and (omag - magmin) % rord != 0:
                lexp = False

            if x / math.pow(10.0, omag) == 1.0 and lexp is True:
                tckstr = "10$^{{{:.0f}}}$".format(omag)
            else:
                tckstr = ""
            return tckstr
        cb.formatter = FuncFormatter(exp_ticks)
        cb.set_label(r'{:s} (${:s}$)'.format(name, units))
    else:
        if zscale > 1.0e3 or zscale < 1.0e-3:
            omag = find_order_of_magnitude(zmax)

            def scaled_ticks(x, pos):
                '''
                Define ticks so that they are scaled by the order of magnitude.
                The two arguements are the value (x) and the tick position (pos)
                and are required for this function to be used by FuncFormatter
                '''
                tckstr = "{:.1f}".format(x / math.pow(10.0, omag))
                return tckstr

            # Use the previously defined function to scale the ticks
            cb.formatter = FuncFormatter(scaled_ticks)
            # Set the label
            cb.set_label(r'{:s} (${:s} \times 10^{{{:.0f}}}$)'.format(name,
                                                                      units,
                                                                      omag))
        else:
            if zscale < 1.0e1:
                cb.formatter=FormatStrFormatter('%.2f')
            else:
                cb.formatter=FormatStrFormatter('%.0f')
            # Set the label
            cb.set_label(r'{:s} (${:s}$)'.format(name, units))

    # Update the ticks to reflect formatting
    cb.update_ticks()
    return cb

def find_order_of_magnitude(value):
    '''
    Find the order of magnitude of a number.  Returns the exponent.
    Ex: -4000.0 = -4 x 10^3 will return 3
    '''
    import math

    return math.floor(math.log10(abs(value)))

def center_polar_cap(rcenter, redge, r):
    '''
    Adjust the radial axis in a polar plot so that it is centered about
    the northern or southern pole
    '''

    if(rcenter > redge):
        return rcenter - r
    else:
        return r

def find_data_limits(gDataList, xkey, lat_index=-1, lon_index=-1, alt_index=-2,
                     inc=6, raw=False, *args, **kwargs):
    '''
    Establish the appropriate axis limits for a list of GitmBin files at
    a particular latitude/longitude index

    Input: gDataList = A list of GitmBin data structures
           xkey      = key for the desired values
           lat_index = latitude index (default -1 for no index)
           lon_index = longitude index (default -1 for no index)
           alt_index = altitude index (default -2 for no index, -1 for 2D)
           inc       = number of tick incriments (default is 6)
           raw       = Keep limits (True) or round of data limits to a more
                       palatable value (False)? (default=False)
    '''
    import math

    hold_min = []
    hold_max = []

    for gData in gDataList:
        if(lat_index < 0 and lon_index < 0):
            if(alt_index > -2):
                flat = gData[xkey][:,:,alt_index].reshape(-1)
            else:
                flat = gData[xkey][:,:,:].reshape(-1)
        elif(lat_index < 0):
            if(alt_index > -2):
                flat = gData[xkey][lon_index,:,alt_index].reshape(-1)
            else:
                flat = gData[xkey][lon_index,:,:].reshape(-1)
        elif(lon_index < 0):
            if(alt_index > -2):
                flat = gData[xkey][:,lat_index,alt_index].reshape(-1)
            else:
                flat = gData[xkey][:,lat_index,:].reshape(-1)
        else:
            if(alt_index > -1):
                flat = gData[xkey][lon_index,lat_index,alt_index].reshape(-1)
            else:
                flat = gData[xkey][lon_index,lat_index,:].reshape(-1)

        hold_min.append(min(flat))
        hold_max.append(max(flat))

    xmin = min(hold_min)
    xmax = max(hold_max)
    xran = round((xmax-xmin)/inc)

    if(xran != 0.0 and raw is False):
        xmin = math.floor(float("%.14f" % (xmin / xran))) * xran
        xmax = math.ceil(float("%.14f" % (xmax / xran))) * xran

    # Consider physical limits for Latitude and Longitude keys.
    if(xkey == "dLat"):
        if(xmin < -90.0):
            xmin = -90.0
        if(xmax > 90.0):
            xmax = 90.0
    elif(xkey == "Latitude"):
        if(xmin < -np.pi / 2.0):
            xmin = -np.pi / 2.0
        if(xmax > np.pi / 2.0):
            xmax = np.pi / 2.0
    elif(xkey == "dLon" or xkey == "Longitude"):
        if(xmin < 0.0):
            xmin = 0.0
        if(xkey == "dLon" and xmax > 360.0):
            xmax = 360
        elif(xkey == "Longitude" and xmax > np.pi):
            xmax = np.pi

    return xmin, xmax

def find_data_limits_irange(gDataList, xkey, min_ilat=-1, max_ilat = -1,
                            min_ilon=-1, max_ilon=-1, min_ialt=-2, max_ialt=-2,
                            inc=6, *args, **kwargs):
    '''
    Establish the appropriate axis limits for a list of GitmBin files at
    a specified range of indexes.  If you only want one index, the maximum
    index should have a value on integer higher than the minimum index, which
    contains the desired index integer.

    Input: gDataList = A list of GitmBin data structures
           xkey      = key for the desired values
           min_ilat  = minimum latitude index (default -1 for no index)
           max_ilat  = maximum latitude index (default -1 for no index)
           min_ilon  = minimum longitude index (default -1 for no index)
           max_ilon  = maximum llongitude index (default -1 for no index)
           min_ialt  = minimum altitude index (default -2 for no index, -1
                       for 2D)
           max_ialt  = maximum laltitude index (default -2 for no index, -1
                       for 2D)
           inc       = number of tick incriments (default is 6)
    '''
    import math

    hold_min = []
    hold_max = []

    for gData in gDataList:
        if(min_ilat < 0 and min_ilon < 0):
            if(min_ialt > -2):
                flat = gData[xkey][:,:,min_ialt:max_ialt].reshape(-1)
            else:
                flat = gData[xkey][:,:,:].reshape(-1)
        elif(min_ilat < 0):
            if(min_ialt > -2):
                flat = gData[xkey][min_ilon:max_ilon,:,
                                   min_ialt:max_ialt].reshape(-1)
            else:
                flat = gData[xkey][min_ilon:max_ilon,:,:].reshape(-1)
        elif(min_ilon < 0):
            if(min_ialt > -2):
                flat = gData[xkey][:,min_ilat:max_ilat,
                                     min_ialt:max_ialt].reshape(-1)
            else:
                flat = gData[xkey][:,min_ilat:max_ilat,:].reshape(-1)
        else:
            if(min_ialt > -1):
                flat = gData[xkey][min_ilon:max_ilon,min_ilat:max_ilat,
                                   min_ialt:max_ialt].reshape(-1)
            else:
                flat = gData[xkey][min_ilon:max_ilon,min_ilat:max_ilat,
                                   :].reshape(-1)

        hold_min.append(min(flat))
        hold_max.append(max(flat))

    xmin = min(hold_min)
    xmax = max(hold_max)
    xran = round((xmax-xmin)/inc)

    if(xran != 0.0):
        xmin = math.floor(float("%.14f" % (xmin / xran))) * xran
        xmax = math.ceil(float("%.14f" % (xmax / xran))) * xran

    # Consider physical limits for Latitude and Longitude keys.

    if(xkey == "dLat"):
        if(xmin < -90.0):
            xmin = -90.0
        if(xmax > 90.0):
            xmax = 90.0
    elif(xkey == "Latitude"):
        if(xmin < -np.pi / 2.0):
            xmin = -np.pi / 2.0
        if(xmax > np.pi / 2.0):
            xmax = np.pi / 2.0
    elif(xkey == "dLon" or xkey == "Longitude"):
        if(xmin < 0.0):
            xmin = 0.0
        if(xkey == "dLon" and xmax > 360.0):
            xmax = 360
        elif(xkey == "Longitude" and xmax > np.pi):
            xmax = np.pi

    return xmin, xmax

def find_data_limits_ivalues(gDataList, xkey, lat_indices, lon_indices,
                             alt_indices, lat_range=False, lon_range=False, 
                             alt_range=False, inc=6, rvals=True, *args,
                             **kwargs):
    '''
    Establish the appropriate axis limits for a list of GitmBin files at
    a specified list of indexes.  If you want to include the entire range of
    latitudes, longitudes, or altitudes instead of a particular index, set
    the range flag to True.  If this is set to True, the first value in the
    xxx_indices list will be used as the lower index and the second value as
    the upper index.  If you are using indexes for more than one coordinate,
    keep in mind that the indices will be paired unless they do not have the
    same length.

    Example:

    xkey = "e-"
    lat_indices = [1, 3, 5] with lat_range = False
    lon_indices = [0, 2] with lon_range = False
    alt_indices = [0, gdata.attrs['nAlt']] with alt_range = True

    Finds the maximum and minimum electron density over all altitudes at
    (lon index, lat index) locations: (0,1), (2,3), and (2,5)

    Input: gDataList   = A list of GitmBin data structures
           xkey        = key for the desired values
           lat_indices = list of latitude indices (or lower index, upper index)
           lon_indices = list of longitude indices (or lower index, upper index)
           alt_indices = list of altitude indices (or lower index, upper index)
           lat_range   = use range for lat instead of indexes (default is False)
           lat_range   = use range for lon instead of indexes (default is False)
           lat_range   = use range for alt instead of indexes (default is False)
           inc         = number of tick incriments (default is 6)
           rvals       = Round the min and max to a sensible number based
                         on the desired number of incriments?  Will never
                         decrease the maximum or increase the minimum.  Will
                         also assess the max/min for latitude and longitudes
                         to ensure they fall within physically sensible limits.
                         (default is True)

    Output: xmin = miminum data value
            xmax = maximum data value
    '''
    import math
    import sys
    import copy

    # Initialize the variables

    hold_min = []
    hold_max = []
    ilat = -1
    ilon = -1
    ialt = -1
    latlist = copy.deepcopy(lat_indices)
    lonlist = copy.deepcopy(lon_indices)
    altlist = copy.deepcopy(alt_indices)

    # Set the maximum and minimum values for the range coordinates

    if lat_range:
        ilatmax = latlist.pop()
        ilatmin = latlist.pop()

    if lon_range:
        ilonmax = lonlist.pop()
        ilonmin = lonlist.pop()

    if alt_range:
        ialtmax = altlist.pop()
        ialtmin = altlist.pop()

    # Cycle over the index coordinates
    while len(altlist) > 0 or len(lonlist) > 0 or len(latlist) > 0:
        # Initialize the index coordiates
        if len(latlist) > 0:
            ilat = latlist.pop()

        if len(lonlist) > 0:
            ilon = lonlist.pop()

        if len(altlist) > 0:
            ialt = altlist.pop()

        if((not lat_range and ilat == -1) or (not lon_range and ilon == -1)
           or (not alt_range and ialt == -1)):
            print("INPUT ERROR in find_data_limit_ivalues")
            sys.exit(1)

        # Cycle over the GITM data structures
        for gData in gDataList:
            # Make the appropriate data selection
            if lat_range:
                if lon_range:
                    if alt_range:
                        # This is silly to do, but sometimes people are silly
                        flat = gData[xkey][ilonmin:ilonmax,ilatmin:ilatmax,
                                           ialtmin:ialtmax].reshape(-1)
                    else:
                        flat = gData[xkey][ilonmin:ilonmax,ilatmin:ilatmax,
                                           ialt].reshape(-1)
                else:
                    if alt_range:
                        flat = gData[xkey][ilon,ilatmin:ilatmax,
                                           ialtmin:ialtmax].reshape(-1)
                    else:
                        flat = gData[xkey][ilon,ilatmin:ilatmax,
                                           ialt].reshape(-1)
            else:
                if lon_range:
                    if alt_range:
                        flat = gData[xkey][ilonmin:ilonmax,ilat,
                                           ialtmin:ialtmax].reshape(-1)
                    else:
                        flat = gData[xkey][ilonmin:ilonmax,ilat,
                                           ialt].reshape(-1)
                else:
                    if alt_range:
                        flat = gData[xkey][ilon,ilat,
                                           ialtmin:ialtmax].reshape(-1)
                    else:
                        # This is silly to do, but sometimes people are silly
                        flat = gData[xkey][ilon,ilat,ialt].reshape(-1)
        
            # Maintain the maximum and mimimum values for this data selection
            hold_min.append(min(flat))
            hold_max.append(max(flat))

    # Find the data max and min from the individual selection max and mins
    xmin = min(hold_min)
    xmax = max(hold_max)
    xran = round((xmax-xmin)/inc)

    if rvals:
        if(xran != 0.0):
            xmin = math.floor(float("%.14f" % (xmin / xran))) * xran
            xmax = math.ceil(float("%.14f" % (xmax / xran))) * xran

        # Consider physical limits for Latitude and Longitude keys.

        if(xkey == "dLat"):
            if(xmin < -90.0):
                xmin = -90.0
            if(xmax > 90.0):
                xmax = 90.0
        elif(xkey == "Latitude"):
            if(xmin < -np.pi / 2.0):
                xmin = -np.pi / 2.0
            if(xmax > np.pi / 2.0):
                xmax = np.pi / 2.0
        elif(xkey == "dLon" or xkey == "Longitude"):
            if(xmin < 0.0):
                xmin = 0.0
            if(xkey == "dLon" and xmax > 360.0):
                xmax = 360
            elif(xkey == "Longitude" and xmax > np.pi):
                xmax = np.pi

    return xmin, xmax

def glon_to_localtime(ut_datetime, glon, units="degrees"):
    '''
    Routine to compute the local time where the longitude is at a specified
    value for a specified universal time

    ut_datetime = Universal time as a datetime object
    glon        = Longitude
    units       = units of longitude (default degrees)

    Output: lt = local time in hours
    '''

    scale = 1.0 / 15.0 # (hours / degree)
    if units.find("rad") >= 0:
        scale *= (180.0 / np.pi) # Convert to hours / radians
    uth = ut_datetime.hour+(ut_datetime.minute/60.0)+(ut_datetime.second/3600.0)
    lt = glon * scale + uth

    # Ensure that the local time falls between 00:00 and 23:59
    if lt < 0.0:
        lt += 24.0

    if lt >= 24.0:
        lt -= 24.0

    return lt

def localtime_to_glon(ut_datetime, localtime):
    '''
    Routine to compute the longitude where the local time is at a specified
    value for a specified universal time

    ut_datetime = Universal time as a datetime object
    localtime   = Local time in hours
    '''

    uth = ut_datetime.hour+(ut_datetime.minute/60.0)+(ut_datetime.second/3600.0)
    lon = (localtime - uth) * 15.0 # 15 = 360 degrees / 24 hours

    return lon

def find_lon_lat_index(gData, glon, glat, units="degrees"):
    '''
    Routine to locate the appropriate longitude and latitude indexes for a
    given location.  The location may be specified in degrees (default) or
    radians.

    Input:
    gData = GitmBin data structure
    glon  = longitude
    glat  = latitude
    units = units of longitude and latitude (default = degrees)
    '''

    import string

    # Set the keys to look for the location in the appropriate units
    if string.lower(units) == "degrees":
        latkey = "dLat"
        lonkey = "dLon"
    else:
        latkey = "Latitude"
        lonkey = "Longitude"

    # First identify the appropriate Longitude.  All longitude values are
    # the same for any latitude and altitude index.


    for (lonindex,clon) in enumerate(gData[lonkey][:,0,0]):
        if clon >= glon:
            if (clon - glon) > (glon - gData[lonkey][lonindex-1,0,0]):
                lonindex = lonindex - 1
            break
        
    # Next identify the appropriate Latitude at the specified longitude

    for (latindex,clat) in enumerate(gData[latkey][lonindex,:,0]):
        if clat >= glat:
            if (clat - glat) > (glat - gData[latkey][lonindex,latindex-1,0]):
                latindex = latindex - 1
            break

    return(lonindex, latindex)

def retrieve_key_from_web_name(name):
    '''
    A routine to retrieve a GITM key corresponding to a descriptive name.
    '''
    key_dict = {"Altitude":"Altitude", "Argon Mixing Ratio":"Ar Mixing Ratio",
                "[Ar]":"Ar","Methane Mixing Ratio":"CH4 Mixing Ratio",
                "Conduction":"Conduction", "EUV Heating":"EuvHeating",
                "[H]":"H", "[H+]":"H!U+!N", "[He]":"He",
                "H2 Mixing Ratio":"H2 Mixing Ratio", "[He+]":"He!U+!N",
                "Hydrogen Cyanide Mixing Ratio":"HCN Mixing Ratio",
                "Heating Efficiency":"Heating Efficiency",
                "Heat Balance Total":"Heat Balance Total",
                "Latitude (rad)":"Latitude", "Longitude (rad)":"Longitude",
                "[N2]":"N!D2!N", "[N2+]":"N!D2!U+!N",
                "[N+]":"N!U+!N", "[N(2D)]":"N(!U2!ND)",
                "[N(2P)]":"N(!U2!NP)", "[N(4S)]":"N(!U4!NS)",
                "N2 Mixing Ratio":"N2 Mixing Ratio", "[NO]":"NO",
                "[NO+]":"NO!U+!N", "[O(4SP)+]":"O_4SP_!U+!N",
                "[O(1D)]":"O(!U1!ND)", "[O2+]":"O!D2!U+!N",
                "[O(2D)]":"O(!U2!ND)!", "[O(2D)]":"O(!U2!ND)!U+!N",
                "[O(2P)+]":"O(!U2!NP)!U+!N", "[O2]":"O!D2!N",
                "[O(2P)]":"O(!U2!NP)!U+!N", "[O(3P)]":"O(!U3!NP)",
                "Radiative Cooling":"RadCooling", "Neutral Density":"Rho",
                "Tn":"Temperature", "v(East)":"V!Di!N (east)",
                "v(North)":"V!Di!N (north)", "v(Up)":"V!Di!N (up)",
                "u(East)":"V!Dn!N (east)","u(North)":"V!Dn!N (north)",
                "u(Up)":"V!Dn!N (up)", "[e-]":"e-",
                "u(Up, N_2)":"V!Dn!N (up,N!D2!N              )",
                "u(Up, N(4S))":"V!Dn!N (up,N(!U4!NS)           )",
                "u(Up, NO)":"V!Dn!N (up,NO                  )",
                "u(Up, O_2)":"V!Dn!N (up,O!D2!N              )",
                "u(Up, O(3P))":"V!Dn!N (up,O(!U3!NP)           )",
                "Electron Average Energy":"Electron_Average_Energy",
                "Te":"eTemperature", "Ti":"iTemperature",
                "Solar Zenith Angle":"Solar Zenith Angle", "[CO2]":"CO!D2!N",
                "Vertical TEC":"Vertical TEC", "DivJu FL":"DivJu FL",
                "DivJuAlt":"DivJuAlt", "Field Line Length":"FL Length",
                "Electron Energy Flux":"Electron_Energy_Flux",
                "sigmaP":"Pedersen FL Conductance", "Potential":"Potential",
                "SigmaP":"Pedersen Conductance", "Region 2 Current":"Je2",
                "sigmaH":"Hall FL Conductance", "Region 1 Current":"Je1",
                "Hall Conductance":"SigmaH", "Ed1":"Ed1", "Ed2":"Ed2",
                "Vertical Electric Field":"E.F. Vertical",
                "Eastward Electric Field":"E.F. East", "dLat":"Latitude (deg)",
                "Northward Electric Field":"E.F. North",
                "Electric Field Magnitude":"E.F. Magnitude",
                "Magnetic Latitude":"Magnetic Latitude",
                "Magnetic Longitude":"Magnetic Longitude",
                "dLon":"Longitude (deg)", "LT":"Solar Local Time"}

    if name in key_dict:
        return key_dict[name]
    else:
        print("ERROR: unknown data type [", name, "], known names are: ")
        print(list(key_dict.keys()))

def find_alt_index(gData, ilon, ilat, alt, units="km"):
    '''
    Routine to locate the appropriate altitude index in a given array.
    The altitude may be specified in km (default) or m.
    '''
    import string

    if string.lower(units) == "km":
        alt *= 1000.0

    # The GITM arrays are sorted, so we can use searchsorted to find
    # a close index

    ialt = np.searchsorted(gData['Altitude'][ilon,ilat,:], alt)

    # If the search index is zero, it is at or below the minimum altitude
    # and doesn't need to be changed.  If the distance between the desired
    # altitude is closer to the search index than the the previous index,
    # it does not need to be changed either

    if(ialt >= gData.attrs['nAlt'] or
       (ialt> 0 and abs(alt - gData['Altitude'][ilon,ilat,ialt]) >
        abs(alt + gData['Altitude'][ilon,ilat,ialt-1]))):
        # If this location is above the maximum height, return maximum index 
        # by reducing the count by one, or if this location is closer to the
        # previous altitude than the one at this index, reduce the count by one
        ialt -= 1
        
    return(ialt)
#End find_alt_index

def match_cindi_key(in_key, out_type="CINDI"):
    '''
    A routine to retrieve a CINDI/GITM data structure key from a GITM/CINDI key.

    in_key   = Input key (GITM or CINDI) that you want to pair.  May specify ALL
               to retrieve a list of all possible keys.
    out_type = Should the output key(s) be the CINDI (default) or GITM keys?
    '''
    key_dict = {"Altitude":"Alt.", "H!U+!N":"FrH", "He!U+!N":"FrHe",
                "NO!U+!N":"FrNO", "O_4SP_!U+!N":"FracO",
                "O(!U2!NP)!U+!N":"FracO", "e-":"Ni(cm^-3)", "iTemperature":"Ti",
                "dLat":"GLAT", "Magnetic Latitude":"MLAT", "dLon":"GLONG",
                "LT":"SLT", "V!Di!N (zon)":"Vzonal", "V!Di!N (par)":"Vpara",
                "V!Di!N (mer)":"Vmerid"}

    if in_key == "ALL":
        # Return a list of all keys of the desired type
        if out_type == "CINDI":
            return(list(key_dict.values()))
        else:
            return(list(key_dict.keys()))
    else:
        if out_type == "CINDI" and in_key in key_dict:
            return key_dict[in_key]
        elif out_type == "GITM":
            out_list = [k for k, v in key_dict.items() if v == 'Alt.']
            if len(out_list) == 1:
                return(out_list[0])
            else:
                print("WARNING: unknown CINDI data type [",in_key,"]")
                return None
        elif out_type == "CINDI":
            print("WARNING: unknown GITM data type [",in_key,"], known names are:")
            print(list(key_dict.keys()))
            return None
        else:
            print("WARNING: unknown data source [", out_type, "]")
            return None
# END match_cindi_key

def add_geomagnetic_equator(ax, color="k", linestyle="-", markerstyle=None):
    '''
    A routine to add a line along the geomagnetic equator to a plot.  Color,
    linestyle, and markerstyle can be specified using matplotlib symbols.  The
    location of the equator was determined using IGRF-10

    Input: ax          = axis handle
           color       = Output color for equator line (default is black)
           linestyle   = Linestyle for equator line (default is solid)
           markerstyle = Marker style (default is None)
    '''

    meq_lon = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0,
               55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0,
               105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0,
               150.0, 155.0, 160.0, 165.0, 170.0, 175.0, 180.0, 185.0, 190.0,
               195.0, 200.0, 205.0, 210.0, 215.0, 220.0, 225.0, 230.0, 235.0,
               240.0, 245.0, 250.0, 255.0, 260.0, 265.0, 270.0, 275.0, 280.0,
               285.0, 290.0, 295.0, 300.0, 305.0, 310.0, 315.0, 320.0, 325.0,
               330.0, 335.0, 340.0, 345.0, 350.0, 355.0, 360.,]
    meq_lat = [10.7, 10.6, 10.5, 10.3, 9.9, 9.5, 9.0, 8.4, 7.9, 7.5, 7.2, 7.0,
               7.1, 7.2, 7.3, 7.5, 7.7, 7.9, 8.0, 8.1, 8.0, 8.0, 7.9, 7.8, 7.8,
               7.8, 7.9, 8.0, 7.9, 7.8, 7.5, 7.1, 6.5, 5.7, 4.9, 4.0, 3.1, 2.3,
               1.6, 0.9, 0.3, -0.4, -1.0, -1.6, -2.3, -3.0, -3.7, -4.4, -5.1,
               -5.9, -6.7, -7.6, -8.6, -9.6, -10.6, -11.4, -11.9, -12.0, -11.6,
               -10.5, -8.8, -6.5, -3.8, -1.0, 1.9, 4.4, 6.4, 8.0, 9.2, 9.9,
               10.4, 10.7, 10.7,]

    stylestring = "{:s}".format(color)
    if markerstyle:
        stylestring = "{:s}{:s}".format(stylestring, markerstyle)
    if linestyle:
        stylestring = "{:s}{:s}".format(stylestring, linestyle)

    ax.plot(meq_lon, meq_lat, stylestring)

    return
# END add_geomagnetic_equator

try:
    from . import solar_rout as sr

    def add_subsolar_point(ut, ax=None, color='none', edge="k", style="o",
                           size=5):
        '''
        A routine to calculate and add a marker at the subsolar point to a plot.
        Fill color, edge color, and markerstyle can be specified using
        matplotlib symbols.  The location can be returned without adding the
        marker to a plot by setting ax=None.

        Input: ut    = datetime holding the current UT
               ax    = axis handle (default=None)
               color = Marker face (fill) color (default='none' (empty))
               edge  = Marker edge color (default='k' (black))
               style = Marker style (default='o' (circle))
               size  = Marker size (default=5)

        Output: sslon = Subsolar Longitude (degrees)
                sslat = Subsolar Latitude (degrees)
        '''
        sslon, sslat = sr.subsolar_point(ut)

        if ax:
            ax.plot([sslon], [sslat], marker=style, ms=size,
                    markerfacecolor=color, markeredgecolor=edge)

        return(sslon, sslat)
    # END add_subsolar_point

    def add_solar_terminator(ut, alt=0.0, ax=None, color='k', style="-",
                             width=1):
        '''
        A routine to calculate and add a line marking the solar terminator to 
        a plot. Line color, style, and width can be specified using matplotlib
        symbols.  Lists containing the coordinates of this line can be returned
        without adding the marker to a plot by setting ax=None.

        Input: ut    = datetime holding the current UT
               alt   = elevation altitude (default=0.0 m)
               ax    = axis handle (default=None)
               color = Line color (default='k' (black))
               style = Line style (default='-' (solid))
               width = Line width/weight (default=1)

        Output: sslon = Subsolar Longitude (degrees)
                sslat = Subsolar Latitude (degrees)
        '''
        # Get the sunrise and sunset coordinates
        term_lat, rise_lon, set_lon = sr.get_solar_terminator_lat_lon(ut,
                                                                      alt=alt,
                                                                      nlocs=500)

        # Align the sunrise and sunset data into a single array
        tlat = np.array(term_lat)
        if set_lon.index(min(set_lon)) < rise_lon.index(min(rise_lon)):
            set_lon = [l if l < 360.0 else l - 360.0 for l in set_lon]
            tlon = np.array(set_lon)
            set_lon = rise_lon
        else:
            rise_lon = [l if l < 360.0 else l - 360.0 for l in rise_lon]
            tlon = np.array(rise_lon)

        tlat = np.append(term_lat, term_lat[::-1])
        tlon = np.append(tlon, set_lon[::-1])

        # Sort the arrays by longitude
        sindex = np.argsort(tlon)
        tlon = tlon[sindex]
        tlat = tlat[sindex]

        # Pad arrays on each side to ensure linear continuity
        tlon = np.append(tlon - 360.0, [tlon, tlon + 360.0])
        tlat = np.append(tlat, [tlat, tlat])

        if ax:
            ax.plot(tlon, tlat, color=color, linestyle=style, linewidth=width)

        return(tlon, tlat)
# END add_solar_terminator
except:
    print("solar_rout.py unavailable: can't load add_solar_terminator, add_subsolar_point")

try:
    import Pysolar as solar

    def find_sunside_twilight_sza(ut, lat_top, lat_bot, lon_top, lon_bot,
                                  mlon_top,  mlon_bot):
        '''
        A routine to find the maximum angluar distance between the solar
        terminator and the sunlight side with conjugate flux tube feet in
        darkness.  The solar zenith angle corresponding to the sunlight
        boundary is returned.  The geographic coordinates of the upper and
        lower boundaries of interest (ie +/- 15 degrees magnetic latitude) must
        be provided along with the magnetic longitude along these boundaries
        and a universal time day.  The entire day is searched to identify the
        sunlight boundary

        Input: ut       = datetime object containing UT year, month and day
               lat_top  = numpy array containing geog. lat at top boundary
               lat_bot  = numpy array containing geog. lat at bottom boundary
               lon_top  = numpy array containing geog. lon at top boundary
               lon_bot  = numpy array containing geog. lon at bottom boundary
               mlon_top = numpy array containing geog. mlon at top boundary
               mlon_bot = numpy array containing geog. mlon at bottom boundary

        Output: sza = Solar zenith angle in degrees of inner boundary
        '''
        import datetime as dt
        import gitm_loc_rout as glr

        # Initialize output
        sza = 90.0
        bad_time = 0

        for sec in np.arange(0.0, 86400.0, 1800.0):
            # Incriment time
            t = ut + dt.timedelta(0.0, sec)

            # Find the solar zenith angle along the boundaries, decrease the
            # angle so that values that occur after noon are negative

            sza_top = np.array([90.0 - solar.GetAltitude(lat_top[i], l, t)
                                if glon_to_localtime(t, l) < 12.0
                                else -90.0 + solar.GetAltitude(lat_top[i], l, t)
                                for i,l in enumerate(lon_top)])
            sza_bot = np.array([90.0 - solar.GetAltitude(lat_bot[i], l, t)
                                if glon_to_localtime(t, l) < 12.0
                                else -90.0 + solar.GetAltitude(lat_bot[i], l, t)
                                for i,l in enumerate(lon_bot)])

            # Find the magnetic longitudes where SZA is +/- 90.0 by identifying
            # the value closest to the desired SZA
            ti = list()
            bi = list()
            delta, i = glr.find_nearest_value(sza_top, 90.0)
            if(delta < 1.0):
                ti.append(i)
            delta, i = glr.find_nearest_value(sza_top, -90.0)
            if(delta < 1.0):
                ti.append(i)
            delta, i = glr.find_nearest_value(sza_bot, 90.0)
            if(delta < 1.0):
                bi.append(i)
            delta, i = glr.find_nearest_value(sza_bot, -90.0)
            if(delta < 1.0):
                bi.append(i)

            # Find the SZA at the opposite borders by locating the correspoding
            # magnetic longitude at the opposite border for the indices at
            # the solar terminator
            new_sza = list()
            for i in ti:
                delta, j = glr.find_nearest_value(mlon_bot, mlon_top[i])
                if sza_bot[j] >= 0.0:
                    new_sza.append(sza_bot[j])
                else:
                    new_sza.append(-sza_bot[j])

            for i in bi:
                delta, j = glr.find_nearest_value(mlon_top, mlon_bot[i])
                if sza_bot[j] >= 0.0:
                    new_sza.append(sza_top[j])
                else:
                    new_sza.append(-sza_top[j])

            # Test to see if this value is higher than the last
            try:
                smin = np.nanmin(new_sza)
                if sza > smin:
                    sza = smin
            except:
                bad_time += 1

        return(sza, bad_time)
except:
    print("PySolar not installed, cannot load find_sunside_twilight_sza")

# Create a contour input 2D numpy array at a specified geo/mag lat/lon/alt

def create_contour_input_array(hold_key, hold_value, xkey, ykey, zkeys, gdata,
                               xarray, yarray, lonmin=0, lonmax=None, latmin=0,
                               latmax=None, altmin=0, altmax=None,
                               *argv, **kwargs):
    '''
    create_contour_input_array: A routine to create contour input at a 
                                specific location.  For example, plotting 
                                Te at a specific latitude or longitude 
                                between two grid points.
    Input: hold_key   = key of coordinate in which desired location is specified
           hold_value = value of desired, specific coordinate location
           xkey       = key of x-coordinate
           ykey       = key of y-coordinate
           zkeys      = list of keys containing data to be interpolated at
                        the specified location as a function of the x and y
                        coordinates
           gdata      = GitmBin structure
           xarray     = numpy array with the desired x-coordinate positions
           yarray     = numpy array with the desired y-coordinate positions
           lonmin     = Minimum longitude index to include in interpolation grid
                        (default=0).  Set as -1 if the coordinate specified
                        by hold_key is analigous to this coordinate.  For
                        example, if data is desired at declination=10 deg,
                        then a limited range of longitues are desired.
           lonmax     = Maximum longitude index to include in interpolation grid
                        (default=None, will use gdata to find maximum)
           latmin     = Minimum latitude index to include in interpolation grid
                        (default=0)  Set as -1 if the coordinate specified
                        by hold_key is analigous to this coordinate.  For
                        example, if data is desired at mlat=0 deg,
                        then a limited range of latitues are desired.
           latmax     = Maximum latitude index to include in interpolation grid
                        (default=None, will use gdata to find maximum)
           altmin     = Minimum altitude index to include in interpolation grid
                        (default=0)  Set as -1 if the coordinate specified
                        by hold_key is analigous to this coordinate.  For
                        example, if data is desired at P=1e10-14 Pa,
                        then a limited range of atlitues are desired.
           altmax     = Maximum atlitude index to include in interpolation grid
                        (default=None, will use gdata to find maximum)
    Output: out = a dictionary containing xarray and the interpolated data
                  for the keys specified in ykeys
    '''
    from scipy import interpolate
    module_name = "create_contour_input_array"

    # Establish the lon/lat/alt indices to navigate the GitmBin structure
    if not lonmax:
        lonmax = gdata.attrs['nLon']

    if not latmax:
        latmax = gdata.attrs['nLat']

    if not altmax:
        altmax = gdata.attrs['nAlt']

    if lonmin < 0:
        iflag = 0
        smax = gdata.attrs['nLon']
        itermax = [[np.searchsorted(gdata[hold_key][:,ilat,ialt], hold_value)
                    for ilat in range(gdata.attrs['nLat'])]
                   for ialt in range(gdata.attrs['nAlt'])]
    elif latmin < 0:
        iflag = 1
        smax = gdata.attrs['nLat']
        itermax = [[np.searchsorted(gdata[hold_key][ilon,:,ialt], hold_value)
                    for ilon in range(gdata.attrs['nLon'])]
                   for ialt in range(gdata.attrs['nAlt'])]
    else:
        iflag = 2
        smax = gdata.attrs['nAlt']
        itermax = [[np.searchsorted(gdata[hold_key][ilon,ilat,:], hold_value)
                    for ilon in range(gdata.attrs['nLon'])]
                   for ilat in range(gdata.attrs['nLat'])]

    # Set up the locations where the data will be interpolated
    out = dict()
    out[xkey], out[ykey] = np.meshgrid(xarray, yarray)

    # Iterate through one data type
    for i,ilist in enumerate(itermax):
        # Ensure the minimum and maximum are not out of bounds
        imin = np.min(ilist)
        while imin <= 0:
            ilist.pop(ilist.index(imin))
            imin = np.min(ilist)

        imax = np.max(ilist)
        while imax >= smax:
            ilist.pop(ilist.index(imax))
            imax = np.max(ilist)

        if iflag == 2:
            altmax = imax + 1
            altmin = imin - 1
            latmin = i
            latmax = i + 1
        else:
            altmin = i
            altmax = i + 1
            if iflag == 1:
                latmax = imax + 1
                latmin = imin - 1
            else:
                lonmax = imax + 1
                lonmin = imin - 1

        # Set up a grid to use to hold the data to be interpolated
        hdata = gdata[hold_key][lonmin:lonmax,latmin:latmax,
                                altmin:altmax].flatten()
        xdata = gdata[xkey][lonmin:lonmax,latmin:latmax,altmin:altmax].flatten()
        N = len(hdata)
        points = np.ndarray(shape=(N, 2), dtype=float,
                            buffer=np.array([[h, xdata[j]]
                                             for j,h in enumerate(hdata)]))

        xinc = out[xkey][i,:].flatten()
        xi = np.ndarray(shape=(len(xinc), 2), dtype=float,
                        buffer=np.array([[hold_value, x] for x in xinc]))

        # Cycle through the z keys to find the desired data values at the
        # desired locations
        for zk in zkeys:
            values = np.ndarray(shape=N, dtype=float, buffer=np.array(gdata[zk][lonmin:lonmax,latmin:latmax,altmin:altmax].flatten()))
            v = interpolate.griddata(points, values, xi, method="linear")
            if zk in out:
                out[zk] = np.append(out[zk], v)
            else:
                out[zk] = np.array(v)

    for zk in zkeys:
        out[zk] = out[zk].reshape(out[xkey].shape)

    return out
# End create_contour_input_array

def create_linear_input_array(hold_key, hold_value, xkey, ykeys, gdata, xarray,
                              lonmin=0, lonmax=1, latmin=0, latmax=1, altmin=0,
                              altmax=1, *argv, **kwargs):
    '''
    create_linear_input_array: A routine to create linear input at a 
                                specific location.  For example, plotting 
                                VTEC at a specific latitude or longitude 
                                between two grid points.
    Input: hold_key   = key of coordinate in which desired location is specified
           hold_value = value of desired, specific coordinate location
           xkey       = key of x-coordinate
           ykeys      = list of keys containing data to be interpolated at
                        the specified location as a function of the x-coordiante
           gdata      = GitmBin structure
           xarray     = numpy array with the desired x-coordinate positions
           lonmin     = Minimum longitude index to include in interpolation grid
                        (default=0)
           lonmax     = Maximum longitude index to include in interpolation grid
                        (default=1)
           latmin     = Minimum latitude index to include in interpolation grid
                        (default=0)
           latmax     = Maximum latitude index to include in interpolation grid
                        (default=1)
           altmin     = Minimum altitude index to include in interpolation grid
                        (default=0)
           altmax     = Maximum atlitude index to include in interpolation grid
                        (default=1)
    Output: out = a dictionary containing xarray and the interpolated data
                  for the keys specified in ykeys
    '''
    from scipy import interpolate
    module_name = "create_linear_input_array"

    # Set up the locations where the data will be interpolated
    out = dict()
    out[xkey] = xarray

    # Set up a grid to use to hold the data to be interpolated
    hdata = gdata[hold_key][lonmin:lonmax,latmin:latmax,altmin:altmax].flatten()
    xdata = gdata[xkey][lonmin:lonmax,latmin:latmax,altmin:altmax].flatten()
    N = len(hdata)
    points = np.ndarray(shape=(N, 2), dtype=float,
                        buffer=np.array([[h, xdata[j]]
                                         for j,h in enumerate(hdata)]))

    xi = np.ndarray(shape=(len(xarray), 2), dtype=float,
                    buffer=np.array([[hold_value, x] for x in out[xkey]]))

    # Cycle through the y keys to find the desired data values at the
    # desired locations
    for yk in ykeys:
        values = np.ndarray(shape=N, dtype=float, buffer=np.array(gdata[yk][lonmin:lonmax,latmin:latmax,altmin:altmax].flatten()))
        v = interpolate.griddata(points, values, xi, method="linear")
        if yk in out:
            out[yk] = np.append(out[yk], v)
        else:
            out[yk] = np.array(v)

    return out
# End create_linear_input_array

def add_dipole_fieldline(ax, apex_alt, lat_array, alt_units="km",
                         lat_units="degrees", lat_type="geographic",
                         longitude=0.0, lon_units="degrees", declination=0.0,
                         dec_units="degrees", color="k", linestyle="-",
                         markerstyle=None):
    '''
    A routine to add a dipole field line to a plot.  Color, linestyle, and
    markerstyle can be specified using matplotlib symbols.  If the latitude
    is geographic, the longitude and declination must be specified so that the
    equatorial offset can be calculated.  Otherwise these are not used.

    Input: ax          = axis handle
           apex_alt    = Field line apex altitude
           lat_array   = list of latitudes to calculate the dipole altitude at
           alt_units   = Altitude units (default is km)
           lat_units   = Latitude units (default is degrees)
           lat_type    = Latitude type: geographic (default), magnetic,
                         inclination
           longitude   = Longitude at the geomagnetic equator (default is 0.0)
           lon_units   = Longitude units (default is degrees)
           longitude   = Longitude at the geomagnetic equator (default is 0.0)
           lon_units   = Longitude units (default is degrees)
           color       = Output color for equator line (default is black)
           linestyle   = Linestyle for equator line (default is solid)
           markerstyle = Marker style (default is None)
    '''
    import math

    # Set the latitude unit conversion
    ldeg = 1.0
    lrad = np.pi / 180.0
    if lat_units.find("rad") >= 0:
        ldeg = 180.0 / np.pi
        lrad = 1.0

    if lat_type.find("geog") >= 0: 
        # Determine the latitude of the geomagnetic equator
        if lon_units.find("rad") >= 0:
            longitude = math.degrees(longitude)
        meq_offset = math.radians(get_meq_offset(longitude))

        # Calculate the cosine of the declination
        if dec_units.find("deg") >= 0:
            declination *= (np.pi / 180.0)
        cos_dec = math.cos(declination)

    # Calculate the height ratio: (apex_alt + earth_rad) / earth_rad
    Re = 6378.1370 # Equatorial radius in km
    if apex_alt == "m":
        Re *= 1000.0
    elif apex_alt == "L":
        Re = 1.0
    
    hsum = (apex_alt + Re)

    # Calculate the dipole field line heights at the desired locations
    dipole_alt = list()
    for lat in lat_array:
        if lat_type.find("geog") >= 0:
            angle = (lat * lrad - meq_offset) / cos_dec
        elif lat_type.find("mag") >= 0:
            angle = lat * lrad
        elif lat_type.find("inc") >= 0:
            angle = math.atan(0.5 * math.tan(lat * lrad))
        else:
            print(func_name, "ERROR: unknown latitude type")
            return

        dipole_alt.append(hsum * math.pow(math.cos(angle), 2.0) - Re)

    # Format the plotting options
    stylestring = "{:s}".format(color)
    if markerstyle:
        stylestring = "{:s}{:s}".format(stylestring, markerstyle)
    if linestyle:
        stylestring = "{:s}{:s}".format(stylestring, linestyle)

    ax.plot(lat_array, dipole_alt, stylestring)

    return(dipole_alt)
# END add_dipole_fieldline


def get_meq_offset(longitude):
    '''
    A routine to determine the geographic latitude in degrees at the
    geomagnetic equator for a specified longitude (provided in degrees).
    Longitude may be a single value, a list, or a numpy array.  The location of
    the equator was determined using IGRF-10.

    Input: longitude = single or multiple longitudes in degrees
    Output: offset = numpy.ndarray containing the geographic latitude in 
                     degrees at the geomagnetic equator
    '''
    from scipy import interpolate

    meq_lon = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0,
               55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0,
               105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0,
               150.0, 155.0, 160.0, 165.0, 170.0, 175.0, 180.0, 185.0, 190.0,
               195.0, 200.0, 205.0, 210.0, 215.0, 220.0, 225.0, 230.0, 235.0,
               240.0, 245.0, 250.0, 255.0, 260.0, 265.0, 270.0, 275.0, 280.0,
               285.0, 290.0, 295.0, 300.0, 305.0, 310.0, 315.0, 320.0, 325.0,
               330.0, 335.0, 340.0, 345.0, 350.0, 355.0, 360.,]
    meq_lat = [10.7, 10.6, 10.5, 10.3, 9.9, 9.5, 9.0, 8.4, 7.9, 7.5, 7.2, 7.0,
               7.1, 7.2, 7.3, 7.5, 7.7, 7.9, 8.0, 8.1, 8.0, 8.0, 7.9, 7.8, 7.8,
               7.8, 7.9, 8.0, 7.9, 7.8, 7.5, 7.1, 6.5, 5.7, 4.9, 4.0, 3.1, 2.3,
               1.6, 0.9, 0.3, -0.4, -1.0, -1.6, -2.3, -3.0, -3.7, -4.4, -5.1,
               -5.9, -6.7, -7.6, -8.6, -9.6, -10.6, -11.4, -11.9, -12.0, -11.6,
               -10.5, -8.8, -6.5, -3.8, -1.0, 1.9, 4.4, 6.4, 8.0, 9.2, 9.9,
               10.4, 10.7, 10.7,]

    # Set up the interpolation curve
    tck = interpolate.splrep(meq_lon, meq_lat, s=0)

    # Return the data at the desired point(s)
    return(interpolate.splev(longitude, tck, der=0))
# END add_meq_offset
