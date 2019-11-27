#!/usr/bin/env python
#-----------------------------------------------------------------------------
# $Id: solar_rout.py,v 1.1 2013/12/23 15:53:07 agburr Exp $
# solar_rout
#
# Author: Angeline G. Burrell, UMichigan, Dec 2013
#
# Comments: Routines to calculate the solar terminator, subsolar point
#-----------------------------------------------------------------------------

import Pysolar as solar
import numpy as np
import math

def subsolar_point(ut):
    '''
    A routine to calculate the subsolar point. 

    Input: ut = datetime holding the current UT

    Output: sslon = Subsolar Longitude (degrees)
            sslat = Subsolar Latitude (degrees)
    '''
    # subsolar latitude is given by the solar declination
    doy = solar.GetDayOfYear(ut)
    sslat = solar.GetDeclination(doy)

    # subsolar longitude occurs at noon SLT
    hour = ut.hour + ut.minute / 60.0 + ut.second / 3600.0
    sslon = 15.0 * (12.0 - hour)

    if sslon > 360.0:
        sslon -= 360.0
    if sslon < 0.0:
        sslon += 360.0

    return(sslon, sslat)

def lat_lon2spherical_xyz(lat, lon, alt):
    '''
    Transformation from angular coordinates assuming a spherical earth to
    spherical cartesian coordinates.

    Input: lat = latitude in degrees
           lon = longitude in degrees
           alt = altitude in m
    '''
    r = solar.constants.earth_radius + alt
    xyz = list()
    xyz.append(r * math.cos(math.radians(lat)) * math.cos(math.radians(lon)))
    xyz.append(r * math.cos(math.radians(lat)) * math.sin(math.radians(lon)))
    xyz.append(r * math.sin(math.radians(lat)))

    return(xyz)

def spherical_xyz2lat_lon(x, y, z):
    '''
    Transformation from spherical cartesian coordinates to angular coordinates
    assuming a spherical earth.
    '''
    alt = math.sqrt(x**2 + y**2 + z**2);
    lon = math.degrees(math.atan2(y, x))
    lat = math.degrees(math.asin(z / alt))
    alt -= solar.constants.earth_radius # in meters

    if(lon < 0.0):
        lon += 360.0

    return(lon, lat, alt)

def get_solar_terminator_xyz_matrix(ut, alt=0.0):
    '''
    Input: ut  = datetime object containing universial time
           alt = altitude in m above the surface of the earth (default=0.0)

    Output: t_mag  = Magnitude of T vector
            t_unit = T unit vector (3-element list)
            s_mag  = Magnitude of S vector
            s_unit = S unit vector (3-element list)
            z_mag  = Magnitude of S vector
            z_unit = S unit vector (3-element list)
    '''
    # Initialize the time-dependent variables
    jd = solar.julian.GetJulianDay(ut)
    jde = solar.julian.GetJulianEphemerisDay(jd, 65)
    jce = solar.julian.GetJulianEphemerisCentury(jde)
    jme = solar.julian.GetJulianEphemerisMillenium(jce)
    #doy = solar.GetDayOfYear(ut)

    # Get the position the sun at the current time
    (sslon, sslat) = subsolar_point(ut)

    #------------------------------------------------------------------------
    # Find the TSZ matrix, the transformation matrix that converts the
    # Cartesian XYZ earth centric coordinate system to the cartesian TSZ
    #  earth centric coordinate system.
    #
    # T: Terminal Plane Normal - the axis of the TSZ system that runs
    #    normal to the plane that slices the earth along the terminal
    #    sunlight line.  The positive T vector points towards noon.
    #
    # S: SunRise/sunSet - the axis of the TSZ system that runs along the
    #    terminal plane and orthogonal to the other axes.  The positive
    #    axis points towards sunset.
    #
    # Z: the vertical axis that bisects the earth at the north and south
    #    points where the sunrise time equals the sunset time.  This axis
    #    also runs along the terminal plane. Positive is towards the north
    #------------------------------------------------------------------------
    # Find the T axis using the projected position of the sun as the noon
    # point. 
    #------------------------------------------------------------------------

    t_vect = lat_lon2spherical_xyz(sslat, sslon, alt)
    # Get the unit vector
    t_mag = math.sqrt(t_vect[0]**2 + t_vect[1]**2 + t_vect[2]**2)
    t_unit = [t / t_mag for t in t_vect]

    # Find the Z axis using 90 degrees plus the position of the subsolar point
    # at the same longitude as the point where sunrise and sunset times are
    # equal.  If the northern sunrise/sunset equality point is in the opposite 
    # hemisphere as the subsolar point (as it is during the June solstice) then
    # we'll use the southern sunrise/sunset point for our calc. 
    plat = sslat + 90.0
    sign = 1.0

    if plat > 90.0:
        plat -= 180.0
        sign = -1.0

    z_vect = lat_lon2spherical_xyz(plat, sslon, alt)
    # Get the unit vector
    z_mag = math.sqrt(z_vect[0]**2 + z_vect[1]**2 + z_vect[2]**2)
    z_unit = [sign * z / z_mag for z in z_vect]

    # Take the cross product of T and Z to find S axis
    s_unit = np.cross(z_unit, t_unit)

    # Set the s mag to the equatorial radius
    s_mag = solar.constants.earth_radius

    return(t_mag, t_unit, s_mag, s_unit, z_mag, z_unit)

def get_terminator_lat_lon_coordinates(s_vect, z_vect, num_latlon=10):
    '''
    Input: s_vect     = S vector or unit vector of TSZ coordinate system
           z_vect     = Z vector  or unit vector of TSZ coordinate system
           num_latlon = number of sunrise/set locations to output (default=10)

    Output: term_lat = Terminator latitude (degrees)
            rise_lon = Sunrise longitude (degrees)
            set_lon  = Sunset longiutde (degrees)
    '''
    zinc = np.arange(-1.0, 1.0, 2.0 / num_latlon, dtype=float)
    term_lat = list()
    rise_lon = list()
    set_lon = list()

    for z_scale in zinc:
        # When finding terminator coordinates, t component is zero
        s_scale = math.sqrt(1.0 - z_scale**2)
	xyz = [z_scale * z_vect[i] - s_scale * s for i,s in enumerate(s_vect)]

        # Find the sunrise location
        (gdlon, gdlat, gdalt) = spherical_xyz2lat_lon(xyz[0], xyz[1], xyz[2])

        term_lat.append(float("{:.3f}".format(np.round(gdlat, 3))))
        rise_lon.append(float("{:.3f}".format(np.round(gdlon, 3))))

        # Find the sunset location
	xyz = [s_scale * s + z_scale * z_vect[i] for i,s in enumerate(s_vect)]
        (gdlon, gdlat, gdalt) = spherical_xyz2lat_lon(xyz[0], xyz[1], xyz[2])

        gdlat = float("{:.3f}".format(np.round(gdlat, 3)))

        if term_lat[-1] == gdlat:
            set_lon.append(float("{:.3f}".format(np.round(gdlon, 3))))
        else:
            print("ERROR: can't to find terminator lat at", gdlat, term_lat[-1])
            return(list(), list(), list())

    return(term_lat, rise_lon, set_lon)

def get_solar_terminator_lat_lon(ut, alt=0.0, nlocs=10):
    '''
    This routine computes the sunrise and sunset times closest to the given
    time and location provided.  It follows the procedure outlined on RADEX
    Research Notebooks 3 (p 1-7) and 6 (p 107-109)

    Input: ut    = datetime object containing desired UT
           alt   = altitude above the surface of the earth (default=0.0 m)
           nlocs = number of locs to find the sunrise/set location (default=10)
    '''
    # Get the solar terminator centric xyz coordinates for the current time
    (tmag,tunit,smag,sunit,zmag,zunit)=get_solar_terminator_xyz_matrix(ut,alt)

    sunit = smag * np.array(sunit)
    zunit = zmag * np.array(zunit)

    # Return the desired solar terminator in geographic coordinates
    return(get_terminator_lat_lon_coordinates(sunit, zunit, nlocs))
