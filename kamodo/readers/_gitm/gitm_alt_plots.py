#!/usr/bin/env python
#-----------------------------------------------------------------------------
# gitm_alt_plots
#
# Author: Angeline G. Burrell, UMichigan, Feb 2013
#
# Comments: Routine to make altitude plots of ionospheric and thermospheric
#           characteristics at a specified location from GITM.
#
# AGB Oct 2013, Adapted to use more general plotting subroutines
#
# Includes: gitm_single_alt_image - plots a single linear or location slice as
#                                   a function of altitude
#           gitm_mult_alt_images  - plot multiple locations of linear or 3D
#                                   altitude slices
#           gitm_alt_slices       - plot a single 3D altitude contour with
#                                   several linear slices
#----------------------------------------------------------------------------

'''
Plot data from a 3D GITM file for different spatiotemporal coordinates
'''

# Import modules
import math
import numpy as np
from spacepy.pybats import gitm
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from . import gitm_plot_rout as gpr
from . import plot_alt_profiles as pap
from scipy import interpolate

def gitm_single_alt_image(plot_type, zkey, gData, lat_index=-1, lon_index=-1,
                          title=None, figname=None, draw=True, xkey="dLat",
                          lon_lt=False, hold_key=None, geo_key=None,
                          hold_value=None, color="b", marker="o",  line=":",
                          zmax=None, zmin=None, amax=None, amin=None, xmax=None,
                          xmin=None, zcenter=False, add_hmf2=False, hcolor="k",
                          hline="--", add_fieldlines=False, apex_alts=list(),
                          acolor="k", aline="-", *args, **kwargs):
    '''
    Creates a rectangular or polar map projection plot for a specified latitude
    range.
    Input: plot_type  = key to determine plot type (linear, contour, scatter)
           zkey       = key for z variable (ie 'Vertical TEC')
           gData      = gitm bin structure
           lat_index  = index of constant latitude (default -1, none)
           lon_index  = index of constant longitude (default -1, none)
           title      = plot title
           figname    = file name to save figure as (default is none)
           xkey       = for contour plots specify an x key (default dLat)
                        (options dLat/dLon/Latitude/Longitude)
           lon_lt     = Add local time to the x axis (default=False).  Will
                        cause problems unless xkey is dLon.
           hold_key   = Key to coordinate in which desired location is specified
           geo_key    = Geographic coordinate corresponding to hold_key coord.
           hold_value = value of desired, specific coordinate location
           color      = linear color (default blue) for b&w contour, enter 'k'
           marker     = linear marker type (default circles)
           line       = linear line type (default dotted)
           zmax       = z key maximum (None to compute automatically, default)
           zmin       = z key minimum (None to compute automatically, default)
           amax       = a key maximum (None to compute automatically, default)
           amin       = a key minimum (None to compute automatically, default)
           xmax       = x key maximum for countour plots (or None to compute
                        automatically, default)
           xmin       = x key minimum for contour plots (or None to compute
                        automatically, default)
           zcenter    = Center z axis about zero? (default is False)
           add_hmf2   = Add line showing hmF2? (default=False)
           hcolor     = Line color for hmF2 line (default=black)
           hline      = Line style for hmF2 line (default=dashes)
           add_fieldlines = Add dipole field lines (default=False)
           apex_alts      = List of apex altitudes to plot field lines at
                            (default=list(), will pick a range of apex alts)
           acolor         = Line color for field lines (default="k")
           aline          = Line style for field lines (default="-")
    '''
    # Initialize the altitude and x variable limits
    if(amin == None or amax == None):
        tmin, tmax = gpr.find_data_limits([gData], "Altitude", lat_index,
                                          lon_index, -2, 30)
        if amin == None:
            amin = tmin
            amin = math.ceil(amin / 10000.0) * 10.0
        if amax == None:
            amax = tmax
            amax = math.floor(amax / 10000.0) * 10.0

    if(plot_type.find("linear") < 0 and (xmin == None or xmax == None)):
        if gData[xkey].attrs['scale'].find("exp") >= 0:
            raw = True
        else:
            raw = False

        tmin, tmax = gpr.find_data_limits([gData], xkey, lat_index, lon_index,
                                          -2, 6, raw=raw)
        if xmin == None:
            xmin = tmin
        if xmax == None:
            xmax = tmax

    # Initialize the input data indices
    datadim = [gData.attrs['nAlt']]
    if lat_index >= 0:
        latmin = lat_index
        latmax = latmin + 1
    else:
        latmin = 0
        latmax = gData.attrs['nLat']
        datadim.insert(0, gData.attrs['nLat'])

    if lon_index >= 0:
        lonmin = lon_index
        lonmax = lonmin + 1
    else:
        lonmin = 0
        lonmax = gData.attrs['nLon']
        datadim.insert(0, gData.attrs['nLon'])

    # Initialize the data
    alt_data = np.array(gData['Altitude'][lonmin:lonmax,latmin:latmax,:])
    alt_data = alt_data.reshape(datadim)

    if plot_type.find("linear") >= 0:
        # Perform interpolation to plot data at a specific location instead of a
        # grid point
        if(hold_key in gData and hold_value is not None):
            if geo_key.find("Lon") >= 0:
                lonmin = 0
                lonmax = gData.attrs['nLon']
            elif geo_key.find("Lat") >= 0:
                latmin = 0
                latmax = gData.attrs['nLat']
            i_data = gpr.create_linear_input_array(hold_key, hold_value,
                                                   "Altitude", [zkey], gData,
                                                   alt_data, lonmin=lonmin,
                                                   lonmax=lonmax, latmin=latmin,
                                                   latmax=latmax,
                                                   altmax=gData.attrs['nAlt'])
            x_data = np.array(i_data[zkey])
            if zmin == None:
                xmin = np.nanmin(x_data)
            if zmax == None:
                xmax = np.nanmax(x_data)
        else:
            if(zmin == None or zmax == None):
                if gData[zkey].attrs['scale'].find("exp") >= 0:
                    raw = True
                else:
                    raw = False
                tmin, tmax = gpr.find_data_limits([gData],zkey,lat_index,
                                                  lon_index,-2,6,raw=raw)
                if zmin == None:
                    xmin = tmin
                if zmax == None:
                    xmax = tmax
            x_data = np.array(gData[zkey][lonmin:lonmax,latmin:latmax,:])
            x_data = x_data.reshape(datadim)
        x_name = gData[zkey].attrs['name']
        x_scale = gData[zkey].attrs['scale']
        x_units = gData[zkey].attrs['units']
        z_data = []
        z_name = ""
        z_scale = ""
        z_units = ""
    else:
        if color.find('k') == 0:
            color = False
        else:
            color = True
        marker=True
        line=zcenter

        # Perform interpolation to plot data at a specific location instead of a
        # grid point
        if(hold_key in gData and hold_value is not None):
            ilon = 0
            ilat = 0
            x_data = np.array(gData[xkey][lonmin:lonmax,latmin:latmax,
                                          0]).flatten()
            if geo_key.find("Lon") >= 0:
                ilon = -1
            elif geo_key.find("Lat") >= 0:
                ilat = -1

            i_data = gpr.create_contour_input_array(hold_key, hold_value, xkey,
                                                    "Altitude", [zkey], gData,
                                                    x_data, alt_data[0,:],
                                                    lonmin=ilon, latmin=ilat)
            x_data = np.array(i_data[xkey])
            z_data = np.array(i_data[zkey])
            alt_data = np.array(i_data['Altitude'])

            if zmin == None:
                zmin = np.nanmin(z_data)
            if zmax == None:
                zmax = np.nanmax(z_data)
        else:
            if(zmin == None or zmax == None):
                if gData[zkey].attrs['scale'].find("exp") >= 0:
                    raw = True
                else:
                    raw = False
                tmin, tmax = gpr.find_data_limits([gData],zkey,lat_index,
                                                  lon_index, -2, 6, raw=raw)
                if zmin == None:
                    zmin = tmin
                if zmax == None:
                    zmax = tmax

            x_data = np.array(gData[xkey][lonmin:lonmax,latmin:latmax,:])
            x_data = x_data.reshape(datadim)
            z_data = np.array(gData[zkey][lonmin:lonmax,latmin:latmax,:])
            z_data = z_data.reshape(datadim)
        x_name = gData[xkey].attrs['name']
        x_scale = gData[xkey].attrs['scale']
        x_units = gData[xkey].attrs['units']
        z_name = gData[zkey].attrs['name']
        z_scale = gData[zkey].attrs['scale']
        z_units = gData[zkey].attrs['units']

    # Initialize the new figure
    f, ax = pap.plot_single_alt_image(plot_type, x_data, alt_data/1000.0,
                                      z_data, x_name, x_scale, x_units, "km",
                                      z_name=z_name, z_scale=z_scale,
                                      z_units=z_units, xmin=xmin, xmax=xmax,
                                      amin=amin, amax=amax, zmin=zmin,
                                      zmax=zmax, title=title, draw=False,
                                      color1=color, color2=marker, color3=line)

    # Add local time to longitude axis, if desired
    if lon_lt:
        xfmt = FuncFormatter(gData.lon_lt_ticks)
        ax.xaxis.set_major_formatter(xfmt)
        ax.set_xlabel("Longitude \ Local Time")
        plt.subplots_adjust(bottom=.13)

    # Add hmF2 if desired
    if add_hmf2 and "hmF2" not in gData:
        gData.calc_2dion()
        if "hmF2" not in gData:
            print(module_name, "WARNING: hmF2 data is not available")
            add_hmF2 = False

    if add_hmf2:
        if plot_type.find("linear") >= 0:
            x = np.array([xmin, xmax])
            if(hold_key in gData and hold_value is not None):
                if geo_key.find("Lon") >= 0:
                    x_data = gData[hold_key][:,latmin:latmax, int(gData.attrs['nAlt']/2)].flatten()
                    y_data = gData['hmF2'][:,latmin:latmax,0].flatten()
                elif geo_key.find("Lat") >= 0:
                    x_data = gData[hold_key][lonmin:lonmax,:,int(gData.attrs['nAlt']/2)].flatten()
                    y_data = gData['hmF2'][20:21,:,0].flatten()
                hold = interpolate.interp1d(x_data, y_data)
                hmf2 = hold(hold_value)
            else:
                hmf2 = gData['hmF2'][lonmin:lonmax,latmin:latmax,0]
            y = np.array([hmf2, hmf2])
        else:
            if(hold_key in gData and hold_value is not None):
                if geo_key.find("Lon") >= 0:
                    lonmin = 0
                    lonmax = gData.attrs['nLon']
                elif geo_key.find("Lat") >= 0:
                    latmin = 0
                    latmax = gData.attrs['nLat']
                x = x_data[0,:]
                i_data = gpr.create_linear_input_array(hold_key, hold_value,
                                                       xkey, ['hmF2'], gData,
                                                       x, lonmin=lonmin,
                                                       lonmax=lonmax,
                                                       latmin=latmin,
                                                       latmax=latmax)
                y = np.array(i_data['hmF2'])
            else:
                x = x_data
                y = np.array(gData['hmF2'][lonmin:lonmax,latmin:latmax,0])
                y = y.reshape(datadim[0:-1])
        ax.plot(x, y, color=hcolor, linestyle=hline, linewidth=2)

    # Add field lines, if desired
    if add_fieldlines and plot_type.find("contour")>=0:
        lon = None
        dec = None
        lon_units = "radians"
        lat_units = "degrees"
        lat_type = "geographic"

        if hold_key == None:
            hold_key = "Longitude"
            hold_value = float(gData[hold_key][lonmin:lonmax,0,0])

        if xkey == "Latitude" or xkey == "dLat":
            from scipy import interpolate
            x = x_data[0]

            if xkey == "Latitude":
                lat_units = "radians"
            # If the geographic equator is on the x-axis, find the longitude
            # and declination at the geomagnetic equator
            ikeys = ["Declination"]
            if hold_key != "Longitude" and hold_key != "dLon":
                ikeys.append("Longitude")
                i_data = gpr.create_linear_input_array(hold_key, hold_value,
                                                       "Inclination", ikeys,
                                                       gData, x_data[0],
                                                       lonmax=lonmax,
                                                       lonmin=lonmin,
                                                       latmax=latmax,
                                                       latmin=latmin, altmax=gData.attrs['nAlt'])
            else:
                x = x_data[:,0]
                i_data = dict()
                i_data['Inclination'] = gData['Inclination'][lonmin:lonmax,
                                                             latmin:latmax,0]
                i_data['Inclination'] = i_data['Inclination'].reshape(latmax)
                i_data['Declination'] = gData['Declination'][lonmin:lonmax,
                                                             latmin:latmax,0]
                i_data['Declination'] = i_data['Declination'].reshape(latmax)

            # Remove any NaN
            ikeys.append("Inclination")
            test = list(range(len(i_data[ikeys[0]])))
            for i in ikeys:
                test *= i_data[i]

            good = [i for i,t in enumerate(test) if not np.isnan(t)]

            if len(good) > 2:
                for i in ikeys:
                    i_data[i] = i_data[i].take(good)

                # Interpolate good data
                tckdec = interpolate.splrep(i_data['Inclination'],
                                            i_data['Declination'], s=0)
                dec = interpolate.splev(0.0, tckdec, der=0)

                if "Longitude" in i_data:
                    tcklon = interpolate.splrep(i_data['Inclination'],
                                                i_data['Longitude'], s=0)
                    lon = interpolate.splev(0.0, tcklon, der=0)
                else:
                    lon = hold_value

                if np.isnan(dec) or np.isnan(lon):
                    print("WARNING: unable to interpolate lon and dec at meq")
                    lat_type = None
                else:
                    # Add local time to the title
                    lt = gpr.glon_to_localtime(gData['time'], lon, lon_units)
                    h = int(lt)
                    m = int((lt - h) * 60.0)
                    ax.set_title("{:}, {:02d}:{:02d} SLT".format(title, h, m),
                                 size="medium")
            else:
                print("WARNING: unable to find longitude and declination at meq")
                lat_type = None

        elif xkey == "Inclination":
            lat_type = "inclination"
        elif xkey == "Magnetic Latitude":
            lat_type = "magnetic"
        else:
            print("WARNING: can't output field lines when xaxis is", xkey)
            lat_type = None

        if lat_type is not None:
            if len(apex_alts) < 1:
                # Set default apex altitudes
                apex_alts = [amin+float(i+1)*200.0 for i in
                             range(int(math.ceil((amax-amin)/200.0)))]

            for alt in apex_alts:
                gpr.add_dipole_fieldline(ax, alt, x, lat_units=lat_units,
                                         lat_type=lat_type, longitude=lon,
                                         lon_units=lon_units, declination=dec,
                                         color=acolor, linestyle=aline)
    if draw:
        # Draw to screen.
        if plt.isinteractive():
            plt.draw() #In interactive mode, you just "draw".
        else:
            # W/o interactive mode, "show" stops the user from typing more 
            # at the terminal until plots are drawn.
            plt.show()

    # Save output file
    if figname is not None:
        plt.savefig(figname)

    return f, ax  

# End plot_single_alt_image

def gitm_mult_alt_images(plot_type, zkey, gData, lat_index, lon_index,
                         title=None, figname=None, draw=True, xkey="dLat",
                         color="b", marker="o", line=":", zmax=None, zmin=None,
                         amax=None, amin=None, xmax=None, xmin=None,
                         zcenter=False, add_hmf2=False, hcolor="k", hline="--",
                         *args, **kwargs):
    '''
    Creates a linear or contour altitude map for a specified altitude range.
    A list of latitude and longitude indexes should be specified.  They may
    be of equal length (for paired values), or for a constant value in one
    coordinate, a length of list one can be specified.
    
    Input: plot_type = key to determine plot type (linear, contour)
           zkey      = key for z variable (ie 'Vertical TEC')
           gData     = gitm bin structure
           lat_index = list of latitude indices (empty list for all)
           lon_index = list of longitude indices (empty list for all)
           title     = plot title
           figname   = file name to save figure as (default is none)
           draw      = Draw figure to screen (default is True)
           xkey      = x coordinate for contour plots (default dLat)
           color     = line color for linear plots (default blue)
           marker    = marker type for linear plots (default circle)
           line      = line type for linear plots (default dotted line)
           zmax      = z key maximum (or None to compute automatically, default)
           zmin      = z key minimum (or None to compute automatically, default)
           amax      = a key maximum (or None to compute automatically, default)
           amin      = a key minimum (or None to compute automatically, default)
           xmax      = x key maximum for countour plots (or None to compute
                       automatically, default)
           xmin      = x key minimum for contour plots (or None to compute
                       automatically, default)
           zcenter   = Should the z range be centered about zero (default is
                       False, for uncentered)
           add_hmf2  = Add line showing hmF2? (default=False)
           hcolor    = Line color for hmF2 line (default=black)
           hline     = Line style for hmF2 line (default=dashes)
    '''
    module_name = "gitm_mult_alt_images"

    # Process the index lists
    lat_len = len(lat_index)
    lon_len = len(lon_index)

    if lat_len != lon_len and lat_len > 1 and lon_len > 1:
        print(module_name, "ERROR: improperly paired lat/lon indices")
        return

    if lat_len <= 1:
        y_label = ["{:.1f}$^\circ$ Lon".format(gData['dLon'][l,1,1])
                   for l in lon_index]
    elif lon_len <= 1:
        y_label = ["{:.1f}$^\circ$ Lat".format(gData['dLat'][1,l,1])
                   for l in lat_index]
    else:
        y_label = ["{:.1f}$^\circ$ Lat, {:.1f}$^\circ$ Lon".format(gData['dLat'][1,l,1], gData['dLon'][lon_index[i],1,1]) for i,l in enumerate(lat_index)]

    # Initialize the input data indices
    subindices = [lon_index, lat_index]

    # Initialize the x,y,z variable limits if desired
    alt_index = [0, gData.attrs['nAlt']]
    alt_range = True
    lat_range = False
    lon_range = False
    if lat_len == 0:
        lat_index = [0, gData.attrs['nLat']]
        lat_range = True
    if lon_len == 0:
        lon_index = [0, gData.attrs['nLon']]
        lon_range = True

    if amin == None or amax == None:
        tmin, tmax = gpr.find_data_limits_ivalues([gData],"Altitude",lat_index,
                                                  lon_index,alt_index,
                                                  lat_range=lat_range,
                                                  lon_range=lon_range,
                                                  alt_range=alt_range,
                                                  rvals=False)
        if amin == None:
            amin = math.ceil(tmin / 10000.0) * 10.0
        if amax == None:
            amax = math.floor(tmax / 10000.0) * 10.0

    if zmin == None or zmax == None:
        if gData[zkey].attrs['scale'].find("exp") >= 0:
            rvals = False
        else:
            rvals = True
        tmin, tmax = gpr.find_data_limits_ivalues([gData], zkey, lat_index,
                                                  lon_index, alt_index, 
                                                  lat_range=lat_range,
                                                  lon_range=lon_range,
                                                  alt_range=alt_range,
                                                  rvals=rvals)
        if zmin == None:
            zmin = tmin
        if zmax == None:
            zmax = tmax

    if (xmin == None or xmax == None) and plot_type.find("linear") < 0:
        tmin, tmax = gpr.find_data_limits_ivalues([gData], xkey, lat_index,
                                                  lon_index, alt_index,  
                                                  lat_range=lat_range,
                                                  lon_range=lon_range,
                                                  alt_range=alt_range,
                                                  rvals=True)
        if xmin == None:
            xmin = tmin
        if xmax == None:
            xmax = tmax

    # Initialize the x and z data
    if plot_type.find("linear") >= 0:
        x_data = np.array(gData[zkey])
        x_name = gData[zkey].attrs['name']
        x_scale = gData[zkey].attrs['scale']
        x_units = gData[zkey].attrs['units']
        xmin = zmin
        xmax = zmax
        z_data = []
        z_name = ""
        z_scale = ""
        z_units = ""
    else:
        x_data = np.array(gData[xkey])
        x_name = gData[xkey].attrs['name']
        x_scale = gData[xkey].attrs['scale']
        x_units = gData[xkey].attrs['units']
        z_data = np.array(gData[zkey])
        z_name = gData[zkey].attrs['name']
        z_scale = gData[zkey].attrs['scale']
        z_units = gData[zkey].attrs['units']

        if color.find('k') == 0:
            color = False
        else:
            color = True

        marker=True
        line=zcenter

    # Initialize the new figure
    alt_data = np.array(gData['Altitude'] / 1000.0)
    f, ax = pap.plot_mult_alt_images(plot_type, subindices, x_data, alt_data,
                                     z_data, x_name, x_scale, x_units, "km",
                                     y_label=y_label, z_name=z_name,
                                     z_scale=z_scale, z_units=z_units,
                                     xmin=xmin, xmax=xmax, amin=amin, amax=amax,
                                     zmin=zmin, zmax=zmax, title=title,
                                     figname=None, draw=False, color1=color,
                                     color2=marker, color3=line)

    # Add the hmF2 lines, if desired
    if add_hmf2 and "hmF2" not in gData:
        gData.calc_2dion()
        if "hmF2" not in gData:
            print(module_name, "WARNING: hmF2 data is not available")
            add_hmF2 = False

    if add_hmf2:
        # Initialize lat/lon indexes that won't change
        datadim = list()
        if lon_len == 0:
            datadim.append(gData.attrs['nLon'])
            lonmin = 0
            lonmax = gData.attrs['nLon']
        elif lon_len == 1:
            lonmin = lon_index[0]
            lonmax = lonmin + 1
        
        if lat_len == 0:
            datadim.append(gData.attrs['nLat'])
            latmin = 0
            latmax = gData.attrs['nLat']
        elif lat_len == 1:
            latmin = lat_index[0]
            latmax = latmin + 1

        # Iterate over all subplot axes
        for i,iax in enumerate(ax):
            # Initialize changing lat/lon indexes
            if lon_len > 1:
                lonmin = lon_index[i]
                lonmax = lonmin + 1
            if lat_len > 1:
                latmin = lat_index[i]
                latmax = latmin + 1

            # Initialize the hmF2 data by plot type
            if plot_type.find("linear") >= 0:
                x = np.array([xmin, xmax])
                y = np.array([gData['hmF2'][lonmin:lonmax,latmin:latmax,0],
                              gData['hmF2'][lonmin:lonmax,latmin:latmax,0]])
                y = y.reshape(2)
            else:
                x = x_data[lonmin:lonmax,latmin:latmax,0]
                x = x.reshape(datadim)
                y = np.array(gData['hmF2'][lonmin:lonmax,latmin:latmax,0])
                y = y.reshape(datadim)

            # Plot the hmF2 line
            iax.plot(x, y, color=hcolor, linestyle=hline, linewidth=2)

    if draw:
        # Draw to screen.
        if plt.isinteractive():
            plt.draw() #In interactive mode, you just "draw".
        else:
            # W/o interactive mode, "show" stops the user from typing more 
            # at the terminal until plots are drawn.
            plt.show()

    # Save output file
    if figname is not None:
        plt.savefig(figname)

    return f
# End gitm_mult_alt_images

def gitm_alt_slices(zkey, gData, lat_index, lon_index, title=None, figname=None,
                    draw=True, degrees=True, color="k", marker="o", line=":",
                    zmax=None, zmin=None, amax=None, amin=None, xmax=None,
                    xmin=None, zcolor=True, zcenter=False, add_hmf2=False,
                    hcolor="k", hline="--", *args, **kwargs):
    '''
    Creates a contour altitude map with several linear slices as a function of
    altitude for a specified GITM variable.  A list of latitude and longitude
    indexes should be specified.  One list should consist of a single value,
    the other will be the x variable in the contour plot.  The degree flag
    determines whether x will be ploted in radians or degrees.

    Input: zkey      = key for z variable (ie 'e-')
           gData     = gitm bin structure
           lat_index = list of latitude indices
           lon_index = list of longitude indices
           title     = plot title
           figname   = file name to save figure as (default is none)
           draw      = Draw the figure to screen (default=True)
           degrees   = plot x label in radians (False) or degrees (default True)
           color     = line color for linear plots (default black)
           marker    = marker type for linear plots (default circle)
           line      = line type for linear plots (default dotted line)
           zmax      = z key maximum (or None to compute automatically, default)
           zmin      = z key minimum (or None to compute automatically, default)
           amax      = a key maximum (or None to compute automatically, default)
           amin      = a key minimum (or None to compute automatically, default)
           xmax      = x key maximum for countour plots (or None to compute
                       automatically, default)
           xmin      = x key minimum for contour plots (or None to compute
                       automatically, default)
           zcolor    = Color plot or B&W (default is True for color)
           zcenter   = Should the z range be centered about zero (default is
                       False, for uncentered)
           add_hmf2  = Add line showing hmF2? (default=False)
           hcolor    = Line color for hmF2 line (default=black)
           hline     = Line style for hmF2 line (default=dashes)
    '''
    module_name = "gitm_alt_slices"

    # Process the index lists
    lat_len = len(lat_index)
    lon_len = len(lon_index)
    pnum = max([lat_len, lon_len])

    if(pnum < 1):
        print(module_name, "ERROR: no altitude slices specified")
        return

    if(lat_len > 1 and lon_len > 1):
        print(module_name, "ERROR: one geographic variable must be constant")
        return

    lat_range = False
    lon_range = False
    if lat_len == 1:
        xkey = "Longitude"
        x_indices = lon_index
        lon_index = [0, gData.attrs['nLon']]
        lon_range = True

        if title:
            title = "{:s} at {:5.2f}$^\circ$ N".format(title, gData['dLat'][1,ilat_index[0],1])
        else:
            title = " {:5.2f}$^\circ$ N".format(gData['dLat'][1,lat_index[0],1])

        if degrees:
            xkey = "dLon"
        x_data = np.array(gData[xkey][:,lat_index[0],:])
        alt_data = np.array(gData['Altitude'][:,lat_index[0],:] / 1000.0)
        z_data = np.array(gData[zkey][:,lat_index[0],:])
    else:
        xkey = "Latitude"
        x_indices = lat_index
        lat_index = [0, gData.attrs['nLat']]
        lat_range = True

        if title:
            title = "{:s} at {:5.2f}$^\circ$ E".format(title, gData['dLon'][lon_index[0],1,1])
        else:
            title = "{:5.2f}$^\circ$ E".format(gData['dLon'][lon_index[0],1,1])

        if degrees:
            xkey = "dLat"
        x_data = np.array(gData[xkey][lon_index[0],:,:])
        alt_data = np.array(gData['Altitude'][lon_index[0],:,:] / 1000.0)
        z_data = np.array(gData[zkey][lon_index[0],:,:])

    # Initialize the x,y,z variable limits
    alt_index = [0, gData.attrs['nAlt']]
    alt_range = True

    if amin == None or amax == None:
        tmin, tmax = gpr.find_data_limits_ivalues([gData],"Altitude",lat_index,
                                                  lon_index,alt_index,
                                                  lat_range=lat_range,
                                                  lon_range=lon_range,
                                                  alt_range=alt_range,
                                                  rvals=False)
        if amin == None:
            amin = math.ceil(tmin / 10000.0) * 10.0
        if amax == None:
            amax = math.floor(tmax / 10000.0) * 10.0

    if zmin == None or zmax == None:
        tmin, tmax = gpr.find_data_limits_ivalues([gData], zkey, lat_index,
                                                  lon_index, alt_index,
                                                  lat_range=lat_range,
                                                  lon_range=lon_range,
                                                  alt_range=alt_range,
                                                  rvals=True)
        if zmin == None:
            zmin = tmin
        if zmax == None:
            zmax = tmax

    if xmin == None or xmax == None:
        tmin, tmax = gpr.find_data_limits_ivalues([gData], xkey, lat_index,
                                                  lon_index, alt_index,
                                                  lat_range=lat_range,
                                                  lon_range=lon_range,
                                                  alt_range=alt_range,
                                                  rvals=True)
        if xmin == None:
            xmin = tmin
        if xmax == None:
            xmax = tmax

    # Initialize the new figure
    f, axc, axl = pap.plot_alt_slices(x_data, alt_data, z_data, 0, x_indices,
                                      gData[xkey].attrs['name'],
                                      gData[xkey].attrs['scale'],
                                      gData[xkey].attrs['units'], "km",
                                      gData[zkey].attrs['name'],
                                      gData[zkey].attrs['scale'],
                                      gData[zkey].attrs['units'], xmin=xmin,
                                      xmax=xmax, amin=amin, amax=amax,
                                      zmin=zmin, zmax=zmax, title=title,
                                      draw=False, color=color, marker=marker,
                                      line=line, zcolor=zcolor, zcenter=zcenter)
    # Add hmF2 lines, if desired
    if add_hmf2:
        # Add hmF2 lines to the linear plots
        ilon = lon_index[0]
        ilat = lat_index[0]
        xlen = len(x_indices) - 1
        for i,iax in enumerate(axl):
            if lon_len > 1:
                ilon = x_indices[xlen-i] 
            else:
                ilat = x_indices[xlen-i]
            x = np.array([zmin,zmax])
            y = np.array([gData['hmF2'][ilon,ilat,0],
                          gData['hmF2'][ilon,ilat,0]])
            y = y.reshape(2)

            iax.plot(x, y, color=hcolor, linestyle=hline, linewidth=2)

        # Add hmF2 lines to the contour plot
        datadim = list()
        if lat_len == 1:
            lonmin = 0
            lonmax = gData.attrs['nLon']
            latmin = lat_index[0]
            latmax = latmin + 1
            datadim.append(lonmax)
        else:
            lonmin = lon_index[0]
            lonmax = lonmin + 1
            latmin = 0
            latmax = gData.attrs['nLat']
            datadim.append(latmax)

        x = np.array(x_data[:,0])
        x = x.reshape(datadim)
        y = np.array(gData['hmF2'][lonmin:lonmax,latmin:latmax,0])
        y = y.reshape(datadim)
        axc.plot(x, y, color=hcolor, linestyle=hline, linewidth=2)

    if draw:
        # Draw to screen.
        if plt.isinteractive():
            plt.draw() #In interactive mode, you just "draw".
        else:
            # W/o interactive mode, "show" stops the user from typing more 
            # at the terminal until plots are drawn.
            plt.show()

    # Save output file
    if figname is not None:
        plt.savefig(figname)

    return f, axc, axl
# End gitm_alt_slices
