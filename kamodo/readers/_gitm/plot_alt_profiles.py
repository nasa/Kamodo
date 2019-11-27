#!/usr/bin/env python
#-----------------------------------------------------------------------------
# $Id: plot_alt_profiles.py,v 1.2 2014/02/17 16:57:29 agburr Exp $
# plot_alt_profiles
#
# Author: Angeline G. Burrell, UMichigan, Oct 2013
#
# Comments: Routines to make linear and contour altitude plots.
#
# Includes: plot_single_alt_image - plots a single linear or location slice as
#                                   a function of altitude
#           plot_mult_alt_images  - plot multiple locations of linear or 3D
#                                   altitude slices
#           plot_alt_slices       - plot a single 3D altitude contour with
#                                   several linear slices
#           -----------------------------------------------------------------
#           plot_linear_alt       - plot the linear altitude dependence of a
#                                   quantity
#           plot_3D_alt           - plot the altitude dependence of a quantity
#                                   as the function of another spatiotemporal
#                                   coordinate
#----------------------------------------------------------------------------

'''
Plot data from a 3D GITM file for different spatiotemporal coordinates
'''

# Import modules
import sys
import string 
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.cm import get_cmap
from matplotlib.colors import LogNorm
from matplotlib.ticker import FormatStrFormatter,MultipleLocator
from . import gitm_plot_rout as gpr

def plot_single_alt_image(plot_type, x_data, alt_data, z_data, x_name, x_scale,
                          x_units, alt_units, z_name="", z_scale="", z_units="",
                          xmin=None, xmax=None, amin=None, amax=None, zmin=None,
                          zmax=None, xinc=6, ainc=6, zinc=6, title=None,
                          tloc="t", figname=None, draw=True, color1="b",
                          color2="o", color3=":", xl=True, xt=True, yl=True,
                          yt=True, *args, **kwargs):
    '''
    Creates a rectangular or polar map projection plot for a specified latitude
    range.
    Input: plot_type = key to determine plot type (linear, contour, scatter)
           x_data    = 2D numpy array containing x-axis data
           alt_data  = 2D numpy array containing y-axis altitude data
           z_data    = 1D or 2D numpy array containing data to plot using a
                       color scale or an empty list for a linear plot
           x_name    = Name of x-axis data
           x_scale   = Plot x-axis data using a linear or exponetial scale?
           x_units   = x-axis data units
           alt_units = y-axis altitude units (m or km)
           z_name    = Name of z-axis data (default="")
           z_scale   = Plot z-axis data using a linear or exponetial scale?
                       (default="")
           z_units   = z-axis data units (default="")
           xmin      = minimum value for x variable (default=None)
           xmax      = maximum value for x variable (default=None)
           amin      = minimum value for altitude (default=None)
           amax      = maximum value for altitude (default=None)
           zmin      = minimum value for z variable (default=None)
           zmax      = maximum value for z variable (default=None)
           xinc      = number of tick incriments for x variable (default 6)
           ainc      = number of tick incriments for altitude (default 6)
           zinc      = number of tick incriments for z variable (default 6)
           title     = plot title
           tloc      = title location (t=top, r=right, l=left, b=bottom,
                       default is top)
           figname   = file name to save figure as (default is none)
           draw      = draw to screen? (default is True)
           xkey      = for contour plots specify an x key (default dLat)
                       (options dLat/dLon/Latitude/Longitude)
           color1    = linear color (default blue) or True/False for color/B&W
                       for the contour colorscale
           color2    = linear marker type (default circles) or True/False for
                       including a colorbar
           color3    = linear line type (default dotted) or True/False if the
                       colorscale using in the contour plot should be centered
                       about zero.
           xl        = Include x label (default is True)
           xt        = Include x ticks (default is True)
           yl        = Include y (altitude) label (default is True)
           yt        = Include y ticks (default is True)
    '''
    # Initialize the new figure
    f = plt.figure()
    ax = f.add_subplot(111)

    if(string.lower(plot_type)=="linear"):
        con = plot_linear_alt(ax, x_data, alt_data, x_name, x_scale, x_units,
                              alt_units, xmin=xmin, xmax=xmax, amin=amin,
                              amax=amax, xinc=xinc, ainc=ainc, title=title,
                              tloc=tloc, xl=xl, xt=xt, yl=yl, yt=yt,
                              color=color1, marker=color2, line=color3)
    else:
        con = plot_3D_alt(ax, x_data, alt_data, z_data, x_name, x_scale,
                          x_units, alt_units, z_name, z_scale, z_units,
                          xmin=xmin, xmax=xmax, amin=amin, amax=amax, zmin=zmin,
                          zmax=zmax, xinc=xinc, ainc=ainc, zinc=zinc, cb=color2,
                          color=color1, zcenter=color3, title=title, tloc=tloc,
                          xl=xl, xt=xt, yl=yl, yt=yt, plot_type=plot_type)
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

def plot_mult_alt_images(plot_type, subindices, x_data, alt_data, z_data,
                         x_name, x_scale, x_units, alt_units, y_label=False,
                         z_name="", z_scale="", z_units="", xmin=None,
                         xmax=None, amin=None, amax=None, zmin=None, zmax=None,
                         xinc=6, ainc=6, zinc=6, title=None, tloc="t",
                         figname=None, draw=True, color1="b", color2="o",
                         color3=":", *args, **kwargs):
    '''
    Creates a linear or contour altitude map for a specified altitude range.
    A list of latitude and longitude indexes should be specified.  They may
    be of equal length (for paired values), or for a constant value in one
    coordinate, a length of list one can be specified.
    
    Input: plot_type  = key to determine plot type (linear, contour, scatter)
           subindices = 2D or 3D list of lists containing the index or indices
                        to include in subplots.  How it works:
                        subindices = [[1], [], [2, 3, 4]]
                                     First index of data only uses index=1,
                                     the third index of data will be used to
                                     select data for linear plots with the
                                     third data index= 2, 3, and 4.  The second
                                     data index includes all of the data.
                        subindices = [[1,2], [40, 50]]
                                     Two subplots will be made, with the data
                                     arrays using data[1,40], data[2,50]
                                     for the two subplots, and the third index
                                     will implicitly include all the data
           x_data     = 2D or 3D numpy array containing x-axis data
           alt_data   = 2D or 3D  numpy array containing y-axis altitude data
           z_data     = 1D or 2D numpy array containing data to plot using a
                        color scale or an empty list for a linear plot
           x_name     = Name of x-axis data
           x_scale    = Plot x-axis data using a linear or exponetial scale?
           x_units    = x-axis data units
           alt_units  = y-axis altitude units (m or km)
           y_label    = List of right-side y-axis labels (labeling subplots)
                        or False to provide no labels (default=False)
           z_name     = Name of z-axis data (default="")
           z_scale    = Plot z-axis data using a linear or exponetial scale?
                        (default="")
           z_units    = z-axis data units (default="")
           xmin       = minimum value for x variable (default=None)
           xmax       = maximum value for x variable (default=None)
           amin       = minimum value for altitude (default=None)
           amax       = maximum value for altitude (default=None)
           zmin       = minimum value for z variable (default=None)
           zmax       = maximum value for z variable (default=None)
           xinc       = number of tick incriments for x variable (default 6)
           ainc       = number of tick incriments for altitude (default 6)
           zinc       = number of tick incriments for z variable (default 6)
           title      = plot title
           tloc       = title location (t=top, r=right, l=left, b=bottom,
                        default is top)
           figname    = file name to save figure as (default is none)
           draw       = draw to screen? (default is True)
           xkey       = for contour plots specify an x key (default dLat)
                        (options dLat/dLon/Latitude/Longitude)
           color1     = linear color (default blue) or True/False for color/B&W
                        for the contour colorscale
           color2     = linear marker type (default circles) or True/False for
                        including a colorbar
           color3     = linear line type (default dotted) or True/False if the
                        colorscale using in the contour plot should be centered
                        about zero.
    '''
    module_name = "plot_mult_alt_images"

    # Test the subindices input
    slen = [len(i) for i in subindices]
    pnum = max(slen)

    if(pnum < 1):
        print(module_name, "ERROR: no subplot regions specified")
        return

    if(pnum == 1):
        print(module_name, "WARNING: only one region, better to use plot_single_alt_image")

    for sl in slen:
        if sl > 1 and sl != pnum:
            print(module_name, "ERROR: subindices input format is incorrect")
            return

    # Initialize the x,y,z variable limits if desired
        if len(slen) == 1 or (len(slen) == 2 and slen[1] == 0):
            if xmin == None:
                xmin = np.nanmin(x_data[subindices[0]])
            if xmax == None:
                xmax = np.nanmax(x_data[subindices[0]])
            if amin == None:
                amin = np.nanmin(alt_data[subindices[0]])
            if amax == None:
                amax = np.nanmax(alt_data[subindices[0]])
            if len(z_data) > 0:
                if zmin == None:
                    zmin = np.nanmin(z_data[subindices[0]])
                if zmax == None:
                    zmax = np.nanmax(z_data[subindices[0]])
        elif len(slen) == 2 or (len(slen) == 3 and slen[2] == 0):
            if slen[0] == 0:
                if xmin == None:
                    xmin = np.nanmin(x_data[:,subindices[1]])
                if xmax == None:
                    xmax = np.nanmax(x_data[:,subindices[1]])
                if amin == None:
                    amin = np.nanmin(alt_data[:,subindices[1]])
                if amax == None:
                    amax = np.nanmax(alt_data[:,subindices[1]])
                if len(z_data) > 0:
                    if zmin == None:
                        zmin = np.nanmin(z_data[:,subindices[1]])
                    if zmax == None:
                        zmax = np.nanmax(z_data[:,subindices[1]])
            else:
                if slen[0] == pnum:
                    s0 = subindices[0]
                else:
                    s0 = subindices[0][0]
                if slen[1] == pnum:
                    s1 = subindices[1]
                else:
                    s1 = subindices[1][0]
                if xmin == None:
                    xmin = np.nanmin(x_data[s0,s1])
                if xmax == None:
                    xmax = np.nanmax(x_data[s0,s1])
                if amin == None:
                    amin = np.nanmin(alt_data[s0,s1])
                if amax == None:
                    amax = np.nanmax(alt_data[s0,s1])
                if len(z_data) > 0:
                    if zmin == None:
                        zmin = np.nanmin(z_data[s0,s1])
                    if zmax == None:
                        zmax = np.nanmax(z_data[s0,s1])
        elif len(slen) == 3:
            if slen[0] == 0 and slen[1] == 0:
                if xmin == None:
                    xmin = np.nanmin(x_data[:,:,subindices[2]])
                if xmax == None:
                    xmax = np.nanmax(x_data[:,:,subindices[2]])
                if amin == None:
                    amin = np.nanmin(alt_data[:,:,subindices[2]])
                if amax == None:
                    amax = np.nanmax(alt_data[:,:,subindices[2]])
                if len(z_data) > 0:
                    if zmin == None:
                        zmin = np.nanmin(z_data[:,:,subindices[2]])
                    if zmax == None:
                        zmax = np.nanmax(z_data[:,:,subindices[2]])
            elif slen[0] == 0:
                if slen[1] == pnum:
                    s1 = subindices[1]
                else:
                    s1 = subindices[1][0]
                if slen[2] == pnum:
                    s2 = subindices[2]
                else:
                    s2 = subindices[2][0]
                if xmin == None:
                    xmin = np.nanmin(x_data[:,s1,s2])
                if xmax == None:
                    xmax = np.nanmax(x_data[:,s1,s2])
                if amin == None:
                    amin = np.nanmin(alt_data[:,s1,s2])
                if amax == None:
                    amax = np.nanmax(alt_data[:,s1,s2])
                if len(z_data) > 0:
                    if zmin == None:
                        zmin = np.nanmin(z_data[:,s1,s2])
                    if zmax == None:
                        zmax = np.nanmax(z_data[:,s1,s2])
            elif slen[1] == 0:
                if slen[0] == pnum:
                    s0 = subindices[0]
                else:
                    s0 = subindices[0][0]
                if slen[2] == pnum:
                    s2 = subindices[2]
                else:
                    s2 = subindices[2][0]
                if xmin == None:
                    xmin = np.nanmin(x_data[s0,:,s2])
                if xmax == None:
                    xmax = np.nanmax(x_data[s0,:,s2])
                if amin == None:
                    amin = np.nanmin(alt_data[s0,:,s2])
                if amax == None:
                    amax = np.nanmax(alt_data[s0,:,s2])
                if len(z_data) > 0:
                    if zmin == None:
                        zmin = np.nanmin(z_data[s0,:,s2])
                    if zmax == None:
                        zmax = np.nanmax(z_data[s0,:,s2])
            else:
                print(module_name, "WARNING: Including restrictions in 3D will cause this program to crash unless the input data is 4D (for linear plots) or 5D (for contour/scatter plots).")
                if slen[0] == pnum:
                    s0 = subindices[0]
                else:
                    s0 = subindices[0][0]
                if slen[1] == pnum:
                    s1 = subindices[1]
                else:
                    s1 = subindices[1][0]
                if slen[2] == pnum:
                    s2 = subindices[2]
                else:
                    s2 = subindices[2][0]
                if xmin == None:
                    xmin = np.nanmin(x_data[s0,s1,s2])
                if xmax == None:
                    xmax = np.nanmax(x_data[s0,s1,s2])
                if amin == None:
                    amin = np.nanmin(alt_data[s0,s1,s2])
                if amax == None:
                    amax = np.nanmax(alt_data[s0,s1,s2])
                if len(z_data) > 0:
                    if zmin == None:
                        zmin = np.nanmin(z_data[s0,s1,s2])
                    if zmax == None:
                        zmax = np.nanmax(z_data[s0,s1,s2])
        else:
            print(module_name, "ERROR: subplot index out of range")
            return

    # Initialize the new figure
    f = plt.figure()
    ax = list()
    tl = " "
    f.text(0.01,0.55,"Altitude (${:s}$)".format(alt_units),rotation="vertical")

    if title:
        f.suptitle(title, size="medium")

    # Adjust the figure height to accomadate the number of subplots
    if(pnum > 2):
        fheight = f.get_figheight()
        f.set_figheight(fheight * 0.5 * pnum)

    for snum in reversed(list(range(0, pnum))):
        cl = False
        xl = False
        fnum = (pnum * 100) + 11 + snum
        ax.append(f.add_subplot(fnum))

        try:
            yl = y_label[snum]
        except:
            yl = False

        if(pnum == snum + 1):
            xl = True

        if(snum == 0):
            cl = True

        if len(slen) == 1 or (len(slen) == 2 and slen[1] == 0):
            xdat = np.array(x_data[subindices[0][snum]])
            altdat = np.array(alt_data[subindices[0][snum]])
            if len(z_data) > 0:
                zdat = np.array(z_data[subindices[0][snum]])
        elif len(slen) == 2 or (len(slen) == 3 and slen[2] == 0):
            if slen[0] == 0:
                xdat = np.array(x_data[:,subindices[1][snum]])
                altdat = np.array(alt_data[:,subindices[1][snum]])
                if len(z_data) > 0:
                    zdat = np.array(z_data[:,subindices[1][snum]])
            else:
                if slen[0] == pnum:
                    s0 = subindices[0][snum]
                else:
                    s0 = subindices[0][0]
                if slen[1] == pnum:
                    s1 = subindices[1][snum]
                else:
                    s1 = subindices[1][0]
                xdat = np.array(x_data[s0,s1])
                altdat = np.array(alt_data[s0,s1])
                if len(z_data) > 0:
                    zdat = np.array(z_data[s0,s1])
        elif len(slen) == 3:
            if slen[0] == 0 and slen[1] == 0:
                xdat = np.array(x_data[:,:,subindices[2][snum]])
                altdat = np.array(alt_data[:,:,subindices[2][snum]])
                if len(z_data) > 0:
                    zdat = np.array(z_data[:,:,subindices[2][snum]])
            elif slen[0] == 0:
                if slen[1] == pnum:
                    s1 = subindices[1][snum]
                else:
                    s1 = subindices[1][0]
                if slen[2] == pnum:
                    s2 = subindices[2][snum]
                else:
                    s2 = subindices[2][0]
                xdat = np.array(x_data[:,s1,s2])
                altdat = np.array(alt_data[:,s1,s2])
                if len(z_data) > 0:
                    zdat = np.array(z_data[:,s1,s2])
            elif slen[1] == 0:
                if slen[0] == pnum:
                    s0 = subindices[0][snum]
                else:
                    s0 = subindices[0][0]
                if slen[2] == pnum:
                    s2 = subindices[2][snum]
                else:
                    s2 = subindices[2][0]
                xdat = np.array(x_data[s0,:,s2])
                altdat = np.array(alt_data[s0,:,s2])
                if len(z_data) > 0:
                    zdat = np.array(z_data[s0,:,s2])
            else:
                print(module_name, "WARNING: Including restrictions in 3D will cause this program to crash unless the input data is 4D (for linear plots) or 5D (for contour/scatter plots)")
                if slen[0] == pnum:
                    s0 = subindices[0][snum]
                else:
                    s0 = subindices[0][0]
                if slen[1] == pnum:
                    s1 = subindices[1][snum]
                else:
                    s1 = subindices[1][0]
                if slen[2] == pnum:
                    s2 = subindices[2][snum]
                else:
                    s2 = subindices[2][0]
                xdat = np.array(x_data[s0,s1,s2])
                altdat = np.array(alt_data[s0,s1,s2])
                if len(z_data) > 0:
                    zdat = np.array(z_data[s0,s1,s2])
        else:
            print(module_name, "ERROR: subplot index out of range")
            return

        if(string.lower(plot_type)=="linear"):
            con = plot_linear_alt(ax[-1], xdat, altdat, x_name, x_scale,
                                  x_units, alt_units, xmin=xmin, xmax=xmax,
                                  amin=amin, amax=amax, xinc=xinc, ainc=ainc,
                                  xl=xl, yl=yl, color=color1, marker=color2,
                                  line=color3)
        else:
            con = plot_3D_alt(ax[-1], xdat, altdat, zdat, x_name, x_scale,
                              x_units, alt_units, z_name, z_scale, z_units,
                              xmin=xmin, xmax=xmax, amin=amin, amax=amax,
                              zmin=zmin, zmax=zmax, xinc=xinc, ainc=ainc,
                              zinc=zinc, cb=False, color=color1,
                              zcenter=color3, title=False, tloc=tloc,
                              xl=xl, yl=yl, plot_type=plot_type)


            if plot_type.find("scatter") >= 0:
                cax = con.axes
            else:
                cax = con.ax

            cpr = list(cax.get_position().bounds)
            if cl is True:
                cpr[2] = new_width
                cax.set_position(cpr)
            else:
                # Add and adjust colorbar
                cbar = gpr.add_colorbar(con, zmin, zmax, zinc, "vertical",
                                        z_scale, z_name, z_units)

                bp = list(cbar.ax.get_position().bounds)
                cp = list(cax.get_position().bounds)

                new_width = cp[2] + 0.075
                cp[2] = new_width
                bp[1] = bp[1] + bp[3] * (float(pnum - 1) / 2.0)
                bp[0] = bp[0] + 0.015

                cbar.ax.set_position(bp)
                cax.set_position(cp)

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
# End plot_mult_alt_images

def plot_alt_slices(x_data, alt_data, z_data, xdim, x_indices, x_name, x_scale,
                    x_units, alt_units, z_name, z_scale, z_units, xmin=None,
                    xmax=None, amin=None, amax=None, zmin=None, zmax=None,
                    xinc=6, ainc=6, zinc=6, title=None, figname=None, draw=True,
                    color="k", marker="o", line=":", zcolor=True, zcenter=False,
                    plot_type="contour", *args, **kwargs):
    '''
    Creates a contour altitude map with several linear slices as a function of
    altitude for a specified variable.  A list of latitude and longitude
    indexes should be specified.  One list should consist of a single value,
    the other will be the x variable in the contour plot.  The degree flag
    determines whether x will be ploted in radians or degrees.

    Input: x_data    = 3D numpy array containing x data for the contour plot
           alt_data  = 3D numpy array containing altitude data
           z_data    = 3D numpy array containing data to plot using a
                       color scale and as the x variable in the linear plots
           xdim      = index (0-2) of the dimension to hold the x data constant
           x_indices = List of index values to plot alt and z at a constant x
           x_name    = Name of x data
           x_scale   = Plot x data using a linear or exponetial scale?
           x_units   = x data units
           alt_units = altitude units (m or km)
           z_name    = Name of z data 
           z_scale   = Plot z data using a linear or exponetial scale?
           z_units   = z data units
           xmin      = minimum value for x variable (default=None)
           xmax      = maximum value for x variable (default=None)
           amin      = minimum value for altitude (default=None)
           amax      = maximum value for altitude (default=None)
           zmin      = minimum value for z variable (default=None)
           zmax      = maximum value for z variable (default=None)
           xinc      = number of tick incriments for x variable (default 6)
           ainc      = number of tick incriments for altitude (default 6)
           zinc      = number of tick incriments for z variable (default 6)
           title     = plot title
           figname   = file name to save figure as (default is none)
           color     = line color for linear plots (default black)
           marker    = marker type for linear plots (default circle)
           line      = line type for linear plots (default dotted line)
           zcolor    = Color plot or B&W (default is True for color)
           zcenter   = Should the z range be centered about zero (default is
                        False, for uncentered)
           plot_type = Make a contour or scatter plot
    '''
    module_name = "plot_alt_slices"

    # Process the index lists
    pnum = len(x_indices)

    if(pnum < 1):
        print(module_name, "ERROR: no subplot slices specified")
        return

    # Initialize the x,y,z variable limits
    if(xmin is None):
        xmin = np.nanmin(x_data)
    if(xmax is None):
        xmax = np.nanmax(x_data)

    if(zmin is None):
        zmin = np.nanmin(z_data)
    if(zmax is None):
        zmax = np.nanmax(z_data)

    if zcenter and abs(zmin) != zmax:
        zmax = max(abs(zmin), zmax)
        zmin = -1.0 * zmax

    if(amin is None):
        amin = np.nanmin(alt_data)
    if(amax is None):
        amax = np.nanmax(alt_data)

    # Initialize the new figure
    f = plt.figure()

    if title:
        f.suptitle(title, size="medium")

    # Adjust the figure size to accomadate the number of subplots
    if(pnum > 2):
        fwidth = f.get_figwidth()
        f.set_figwidth(fwidth * 0.5 * pnum)

    # Display the 3D contour plot on top
    gs = gridspec.GridSpec(2,pnum)
    gs.update(hspace=.4)

    axc = plt.subplot(gs[0,:])
    con = plot_3D_alt(axc, x_data, alt_data, z_data, x_name, x_scale, x_units,
                      alt_units, z_name, z_scale, z_units, xmin=xmin, xmax=xmax,
                      amin=amin, amax=amax, zmin=zmin, zmax=zmax, xinc=xinc,
                      ainc=ainc, zinc=zinc, cb=True, cloc="t", color=zcolor,
                      zcenter=zcenter, plot_type=plot_type)
  
    # Determine how many x tics to include in each linear slice
    if(pnum > 4):
        zinc *= 0.5

    # Add the linear slices
    y = [amin, amax]
    axl = list()

    for snum in reversed(list(range(0, pnum))):
        xl = False
        yl = False
        axl.append(plt.subplot(gs[1,snum]))

        if(snum == 0):
            yl = True

        if(math.floor(pnum * 0.5) == snum):
            xl = True

        if xdim == 0:
            xloc = x_data[x_indices[snum]]
            for i in range(len(x_data[x_indices[snum]].shape)):
                xloc = xloc[0]
            tl = "{:.2f} ${:s}$ {:s}".format(xloc, x_units, x_name)
            x = [x_data[x_indices[snum]], x_data[x_indices[snum]]]
            altdat = alt_data[x_indices[snum]]
            zdat = z_data[x_indices[snum]]
        elif xdim == 1:
            xloc = x_data[0,x_indices[snum]]
            for i in range(len(x_data[x_indices[snum]].shape)):
                xloc = xloc[0]
            tl = "{:.2f} ${:s}$ {:s}".format(xloc, x_units, x_name)
            x = [x_data[0,x_indices[snum]], x_data[0,x_indices[snum]]]
            altdat = alt_data[:,x_indices[snum]]
            zdat = z_data[:,x_indices[snum]]
        else:
            tl = "{:.2f} ${:s}$ {:s}".format(x_data[0,0,x_indices[snum]],
                                             x_units, x_name)
            x = [x_data[0,0,x_indices[snum]], x_data[0,0,x_indices[snum]]]
            altdat = alt_data[:,:,x_indices[snum]]
            zdat = z_data[:,:,x_indices[snum]]

        # Add a line to the contour plot
        axc.plot(x, y, color=color, linestyle=line, linewidth=2)

        # Add a linear slice below the contour plot
        plot_linear_alt(axl[-1], zdat, altdat, z_name, z_scale, z_units,
                        alt_units, xmin=zmin, xmax=zmax, amin=amin, amax=amax,
                        xinc=zinc, ainc=ainc, title=tl, tloc="t", xl=xl, yl=yl,
                        color=color, marker=marker, line=line)

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
# End plot_alt_slices

def plot_linear_alt(ax, x_data, alt_data, x_name, x_scale, x_units, alt_units,
                    xmin=None, xmax=None, amin=None, amax=None, xinc=6, ainc=6,
                    title=None, tloc="t", xl=True, xt=True, yl=True, yt=True,
                    color="b", marker="+", line="-", *args, **kwargs):
    '''
    Creates a rectangular map projection plot for a specified latitude range.
    Input: ax        = axis handle
           x_data    = 2D numpy array containing x-axis data
           alt_data  = 2D numpy array containing y-axis altitude data
           x_name    = Name of x-axis data
           x_scale   = Plot x-axis data using a linear or exponetial scale?
           x_units   = x-axis data units
           alt_units = y-axis altitude units (m or km)
           xmin      = minimum value for x variable
           xmax      = maximum value for x variable
           amin      = minimum altitude
           amax      = maximum altitude
           xinc      = number of x variable tick incriments (default 6)
           ainc      = number of alt variable tick incriments (default 6)
           title     = plot title (default is None)
           tloc      = Specify the title location (t=top, r=right, l=left,
                       b=bottom, default is top)
           xl        = Include x (z variable) label (default is True)
           xt        = Include x ticks (default is True)
           yl        = Include y label.  By default this label is placed on 
                       the left to label the altitude.  If a non-Boolian value
                       is provided, this will be assumed to be a string to be
                       placed as a label on the right. (default is True)
           yt        = Include y ticks (default is True)
           color     = line/marker color (default b [blue])
           marker    = marker type (default +)
           line      = line type (default - [solid line])
    '''

    # Set the x, a, and z ranges
    if(xmin is None):
        xmin = np.nanmin(x_data)
    if(xmax is None):
        xmax = np.nanmax(x_data)
    arange = xmax - xmin
    xwidth = arange / xinc

    if(amin is None):
        amin = np.nanmin(alt_data)
    if(amax is None):
        amax = np.nanmax(alt_data)
    arange = amax - amin
    awidth = arange / ainc

    # Plot the values
    con = ax.plot(x_data, alt_data, color=color, marker=marker, linestyle=line)

    # Configure axis
    if yt:
        ytics = MultipleLocator(awidth)
        ax.yaxis.set_major_locator(ytics)
    else:
        ax.yaxis.set_major_formatter(FormatStrFormatter(""))

    if yl is True:
        ax.set_ylabel('Altitude ($km$)')
    elif yl is not False:
        ax.set_ylabel(yl)
        ax.yaxis.set_label_position("right")
    plt.ylim(amin, amax)

    if x_scale.find("exponential"):
        ax.set_xscale('log')
    elif xt:
        xtics = MultipleLocator(xwidth)
        ax.xaxis.set_major_locator(xtics)
    else:
        ax.xaxis.set_major_formatter(FormatStrFormatter(""))

    if xl:
        ax.set_xlabel(r'%s ($%s$)' % (x_name, x_units))
    plt.xlim(xmin, xmax)

    # Set the title
    if title:
        rot  = 'horizontal'
        yloc = 1.05
        xloc = .5

        if(tloc == "l" or tloc == "r"):
            xloc = -.1
            yloc = .5
            rot  = 'vertical'

            if(tloc == "r"):
                xloc = 1.1

        if(tloc == "b"):
            yloc = -.1
            
        ax.set_title(title,size='medium', rotation=rot, y=yloc, x=xloc)

    return con
#End plot_linear_alt

def plot_3D_alt(ax, x_data, alt_data, z_data, x_name, x_scale, x_units,
                alt_units, z_name, z_scale, z_units, xmin=None, xmax=None,
                amin=None, amax=None, zmin=None, zmax=None, xinc=6, ainc=6,
                zinc=6, cb=True, cloc="r", color=True, zcenter=False,
                title=None, tloc="t", xl=True, xt=True, yl=True, yt=True,
                plot_type="contour", *args, **kwargs):
    '''
    Creates a single polar projection, with the latitude center and range
    determined by the input.
    Input: ax        = axis handle
           x_data    = 2D numpy array containing x-axis data
           alt_data  = 2D numpy array containing y-axis altitude data
           z_data    = 1D or 2D numpy array containing data to plot using a
                       color scale
           x_name    = Name of x-axis data
           x_scale   = Plot x-axis data using a linear or exponetial scale?
           x_units   = x-axis data units
           alt_units = y-axis altitude units (m or km)
           z_name    = Name of z-axis data
           z_scale   = Plot z-axis data using a linear or exponetial scale?
           z_units   = z-axis data units
           xmin      = minimum value for x variable (default=None)
           xmax      = maximum value for x variable (default=None)
           amin      = minimum value for altitude (default=None)
           amax      = maximum value for altitude (default=None)
           zmin      = minimum value for z variable (default=None)
           zmax      = maximum value for z variable (default=None)
           xinc      = number of tick incriments for x variable (default 6)
           ainc      = number of tick incriments for altitude (default 6)
           zinc      = number of tick incriments for z variable (default 6)
           cb        = Add a colorbar (default is True)
           cloc      = Colorbar location (t=top, r=right, l=left, b=bottom, 
                       default is right)
           color     = Color plot or B&W (default is True for color)
           zcenter   = Should the z range be centered about zero (default is
                       False, for uncentered)
           title     = plot title (default is none)
           tloc      = title location (t=top, r=right, l=left, b=bottom,
                       default is top)
           xl        = Include x label (default is True)
           xt        = Include x ticks (default is True)
           yl        = Include y label.  This defaults to placing an altitude
                       label on the left axis.  If a non-Boolian value is
                       provided, it is assumed to be a string that will be
                       used as a right axis label.  (default is True)
           yt        = Include y ticks (default is True)
           plot_type = Make a scatter or contour plot? (default=contour)
    '''
    # Set the x, a, and z ranges
    if(xmin is None):
        xmin = np.nanmin(x_data)
    if(xmax is None):
        xmax = np.nanmax(x_data)
    arange = xmax - xmin
    xwidth = arange / xinc

    if(zmin is None):
        zmin = np.nanmin(z_data)
    if(zmax is None):
        zmax = np.nanmax(z_data)

    if zcenter and abs(zmin) != zmax:
        arange = max(abs(zmin), zmax)
        zmax = arange
        zmin = -1.0 * arange

    arange = zmax - zmin
    zwidth = arange / zinc

    if(amin is None):
        amin = np.nanmin(alt_data)
    if(amax is None):
        amax = np.nanmax(alt_data)
    arange = amax - amin
    awidth = arange / ainc

    # Determine the z scale
    if z_scale.find("exp") >= 0:
        v = np.logspace(math.log10(zmin), math.log10(zmax), zinc*10,
                        endpoint=True)
        norm = LogNorm(vmin=zmin, vmax=zmax)
    else:
        norm = None
        v = np.linspace(zmin, zmax, zinc*10, endpoint=True)

    # Plot the data
    col = gpr.choose_contour_map(color, zcenter)
    if plot_type.find("scatter") >= 0:
        con = ax.scatter(x_data, alt_data, c=z_data, cmap=get_cmap(col),
                         norm=norm, vmin=zmin, vmax=zmax, edgecolors="none",
                         s=10)
        cax = con.axes
    else:
        con = ax.contourf(x_data, alt_data, z_data, v, cmap=get_cmap(col),
                          norm=norm, vmin=zmin, vmax=zmax)
        cax = con.ax

    # Configure axis
    if yt:
        ytics = MultipleLocator(awidth)
        ax.yaxis.set_major_locator(ytics)
    else:
        ax.yaxis.set_major_formatter(FormatStrFormatter(""))

    if yl is True:
        ax.set_ylabel('Altitude ($km$)')
    elif yl is not False:
        ax.set_ylabel(yl)
        ax.yaxis.set_label_position("right")
    plt.ylim(amin, amax)

    if x_scale.find("exponential") >= 0:
        ax.set_xscale('log')
    elif xt:
        xtics = MultipleLocator(xwidth)
        ax.xaxis.set_major_locator(xtics)
    else:
        ax.xaxis.set_major_formatter(FormatStrFormatter(""))

    if xl:
        ax.set_xlabel(r'%s ($%s$)' % (x_name, x_units))
    plt.xlim(xmin, xmax)
           
    # Set the title
    if title:
        rot  = 'horizontal'
        yloc = 1.05
        xloc = 0.5

        if tloc == "b":
            yloc = -.1
        elif tloc != "t":
            rot  = 'vertical'
            yloc = 0.5
            xloc = -.2

            if tloc == "r":
                xloc = 1.1

        title = ax.set_title(title,y=yloc,size='medium',x=xloc,rotation=rot)
 
    # Change the background color
    ax.patch.set_facecolor('#747679')

    # Add a colorbar
    if cb:
        orient = 'vertical'

        if(cloc == 't' or cloc == 'b'):
            orient = 'horizontal'

        cbar = gpr.add_colorbar(con, zmin, zmax, zinc, orient, z_scale, z_name,
                                z_units)

        if(cloc == 'l' or cloc == 't'):
            bp = list(cbar.ax.get_position().bounds)
            cp = list(cax.get_position().bounds)

            if(cloc == 't'):
                cp[1] = bp[1]
                bp[1] = cp[1] + cp[3] + 0.085
            else:
                bp[0] = 0.125
                cp[0] = bp[0] + 0.1 + bp[2]

            cax.set_position(cp)
            cbar.ax.set_position(bp)

    return con
#End plot_3D_alt

#End plot_alt_profiles

