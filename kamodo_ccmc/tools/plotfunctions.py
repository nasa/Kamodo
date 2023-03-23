def figMods(fig, log10=False, lockAR=False, ncont=-1, colorscale='',
            cutInside=-1., returnGrid=False, enhanceHover=False, newTitle=''):
    '''
    Function to modify a plotly figure object in multiple ways.

    Arguments:
      fig         A plotly figure object
      log10       Logical, if true take the log of the contour values
      lockAR      Logical, if true lock the X/Y axis aspect ratio to 1
      ncont       number of contours to include in the plot
      colorscale  Set the desired colorscale for the plot
                  NOTE: appending '_r' to colorscale reverses it, ie RdBu_r
      cutInside   Make values inside a radius of cutInside NaNs
      returnGrid  Take the plot and return a new grid only plotly object
    '''

    import re
    import numpy as np
    import plotly.graph_objects as go

    if returnGrid:
        xx = fig.data[0]['x']
        yy = fig.data[0]['y']
        zz = fig.data[0]['z']
        xx_mg, yy_mg = np.meshgrid(xx, yy)
        xm = np.reshape(xx_mg, -1)
        ym = np.reshape(yy_mg, -1)
        fig2 = go.Figure(data=go.Scattergl(x=xm, y=ym, mode='markers'))
        return fig2

    if log10:
        # grab values from old plot
        val = fig.data[0]['z']
        # set negative and zero values to NaN
        # compute log10 of values, NaN will stay NaN
        val[val <= 0.] = np.nan
        val = np.log10(val)
        # assign back to plot object
        fig.data[0]['z'] = val
        # add log10 to colorbar label
        newlabel = 'log10<br>'+fig.data[0]['colorbar']['title']['text']
        fig.data[0]['colorbar']['title']['text'] = newlabel

    if lockAR:
        fig.update_xaxes(scaleanchor='y')

    if ncont > 0:
        fig.update_traces(ncontours=ncont,
                          contours=dict(coloring="fill", showlines=False))

    if colorscale != '':
        if colorscale == "BlueRed":
            fig.update_traces(colorscale="RdBu_r")
        elif colorscale == "Rainbow":
            fig.update_traces(
                colorscale=[[0.00, 'rgb(0,0,255)'],
                            [0.25, 'rgb(0,255,255)'],
                            [0.50, 'rgb(0,255,0)'],
                            [0.75, 'rgb(255,255,0)'],
                            [1.00, 'rgb(255,0,0)']])
        else:
            fig.update_traces(colorscale=colorscale)

    if cutInside > 0.:
        xx = fig.data[0]['x']
        yy = fig.data[0]['y']
        zz = fig.data[0]['z']
        xx_mg, yy_mg = np.meshgrid(xx, yy)  # Now the same dimensions as zz
        rr_mg = np.sqrt(xx_mg*xx_mg + yy_mg*yy_mg)
        maskR = rr_mg < cutInside
        zz[maskR] = np.nan
        fig.data[0]['z'] = zz

    if enhanceHover:
        xvar = re.sub(' +', ' ',
                      fig.layout.xaxis.title.text.strip('$').strip())
        yvar = re.sub(' +', ' ',
                      fig.layout.yaxis.title.text.strip('$').strip())
        cvar = re.sub(' +', ' ',
                      fig.data[0]['colorbar'].title.text.strip('$').strip())
        fig.update_traces(
            hovertemplate=xvar + ": %{x:.3f} <br>" +
            yvar + ": %{y:.3f} <br>" +
            cvar + ": %{z:.4g} <br>" +
            "<extra></extra>"
        )

    if newTitle != '':
        fig.layout.title.text = newTitle

    return fig


def toColor(fig, colorscale='Viridis'):
    '''
    Placeholder function for compatibility with old notebooks.
    Use figMods instead.

    Arguments:
      fig         A plotly figure object
      colorscale  Set the desired colorscale for the plot

    Returns a modified plotly figure object
    '''

    return figMods(fig, colorscale=colorscale, ncont=200)


def toLog10(fig):
    """
    Placeholder function for compatibility with old notebooks.
    Use figMods instead.

    Arguments:
      fig         A plotly figure object

    Returns a modified plotly figure object
    """

    return figMods(fig, log10=True)


def XYC(Xlabel, X, Ylabel, Y, Clabel, C, title='Plot Title',
        colorscale='Viridis', crange='',):
    """
    Simple 2D plot. Send in array of X points, Y points, and 2d array of
      values with labels for each and a basic plot is created.

    Arguments:
      Xlabel     Text label of X points
      X          1D array of points in X direction
      Ylabel     Text label of Y points
      Y          1D array of points in Y direction
      Clabel     Text label of plotted value
      C           2D array of values to plot, dimensioned[len(X),len(Y)]
      title       Optional plot title
      colorscale  Text string of colorscale to use for contours
      crange      Two position array with min/max contour values, [cmin,cmax]

    Returns a plotly figure object.
    """

    import numpy as np
    from kamodo import Kamodo

    # Compute/set contour range
    cmin = np.min(C)
    cmax = np.max(C)
    if crange != '':
        cmin = float(crange[0])
        cmax = float(crange[1])

    # Make basic figure to start from

    def plot_2D(X=X, Y=Y):
        return C

    plot1 = Kamodo(plot_2D=plot_2D)
    fig = plot1.plot(plot_2D=dict())

    # Set colorscale
    fig = figMods(fig, colorscale=colorscale, ncont=200)

    # Set plot options
    fig.update_traces(
        zmin=cmin, zmax=cmax, colorbar=dict(title=Clabel, tickformat=".3g"),
        hovertemplate=Xlabel + ": %{x:.4g}<br>" + Ylabel + ": %{y:.4g}<br>" +
        Clabel + ": %{z:.6g}<br>" + "<extra></extra>"
    )
    fig.update_layout(
        scene_aspectmode='data',
        title=dict(text=title,
                   yref="container", yanchor="top", x=0.01, y=0.97,
                   font=dict(size=16, family="sans serif", color="#000000")),
        margin=dict(l=0), xaxis_title=Xlabel, yaxis_title=Ylabel
    )

    return fig


def ReplotLL3D(figIn, model, altkm, plotts, plotCoord='GEO',
               title='Plot Title', colorscale='Viridis', crange='',
               opacity=0.70, axis=True, debug=0):
    """
    Takes a gridified 2D lon/lat figure and creates new plots
    in 3D and for chosen coordinate systems.
    """

    import numpy as np
    import datetime as dt
    import pytz
    from kamodo import Kamodo
    import kamodo_ccmc.flythrough.model_wrapper as MW
    from kamodo_ccmc.flythrough.utils import ConvertCoord
    from kamodo_ccmc.tools.shoreline import shoreline

    # Get current coordinate system and type. Create string for labels
    tmp = list(MW.Model_Variables(model, return_dict=True).values())[0][2:4]
    co = tmp[0]
    cot = tmp[1]

    # Pull out lon, lat, values and min/max from passed in figure
    lon = figIn.data[0]['x']
    lat = figIn.data[0]['y']
    val = figIn.data[0]['z']
    val2 = np.reshape(val, (len(lat), len(lon)))
    varn = figIn.data[0]['colorbar']['title']['text']
    cmin = np.min(val)
    cmax = np.max(val)
    if crange != '':
        cmin = float(crange[0])
        cmax = float(crange[1])

    # Prepare variables for use later
    rscale = (altkm + 6.3781E3)/6.3781E3
    lon_mg, lat_mg = np.meshgrid(np.array(lon), np.array(lat))
    x = lon_mg
    y = lat_mg
    if co == 'GDZ':
        z = np.full(x.shape, altkm)
    else:
        z = np.full(x.shape, rscale)
    t = np.full(x.shape, plotts)
    full_x = np.reshape(x, -1)
    full_y = np.reshape(y, -1)
    full_z = np.reshape(z, -1)
    full_t = np.reshape(t, -1)
    ilon = np.linspace(-180, 180, 181)
    # ilat = np.linspace(-90, 90, 91)

    # Convert incoming coordinates into plot coordinages (cartesian)
    if co == 'GDZ':
        # GDZ does not convert into other coordinates well, first go to GEO-car
        xx, yy, zz, units = ConvertCoord(full_t, full_x, full_y, full_z,
                                         co, cot, 'GEO', 'car')
        full_x, full_y, full_z = xx, yy, zz
        if plotCoord != 'GEO':
            xx, yy, zz, units = ConvertCoord(full_t, full_x, full_y, full_z,
                                             'GEO', 'car', plotCoord, 'car')
    else:
        xx, yy, zz, units = ConvertCoord(full_t, full_x, full_y, full_z,
                                         co, cot, plotCoord, 'car')
    x = np.reshape(xx, (len(lat), len(lon)))
    y = np.reshape(yy, (len(lat), len(lon)))
    z = np.reshape(zz, (len(lat), len(lon)))
    # NaNs in original plot don't hide properly, so set positions to NaN
    mask2 = np.isnan(val2)
    x[mask2] = np.nan
    y[mask2] = np.nan
    z[mask2] = np.nan

    # Generate initial figure to build upon

    def plot3d_var(X=x, Y=y, Z=z):
        return val2

    kobject = Kamodo(plot_var=plot3d_var)
    fig = kobject.plot(plot_var=dict())

    # Set colorscale
    if colorscale == "BlueRed":
        fig.update_traces(colorscale="RdBu_r")
    elif colorscale == "Rainbow":
        fig.update_traces(
            colorscale=[[0.00, 'rgb(0,0,255)'],
                        [0.25, 'rgb(0,255,255)'],
                        [0.50, 'rgb(0,255,0)'],
                        [0.75, 'rgb(255,255,0)'],
                        [1.00, 'rgb(255,0,0)']])
    else:
        fig.update_traces(colorscale=colorscale)

    # Set plot options
    fig.update_traces(
        cmin=cmin, cmax=cmax, colorbar=dict(title=varn, tickformat=".3g"),
        opacity=opacity,
        customdata=np.dstack((
            np.transpose(val2, axes=[1, 0]),
            np.transpose(lon_mg, axes=[1, 0]),
            np.transpose(lat_mg, axes=[1, 0]))),
        hovertemplate="<b>" + model + "</b><br>" +
        "X: %{x:.2f} " + units[0] + " " + plotCoord + "<br>"
        "Y: %{y:.2f} " + units[1] + " " + plotCoord + "<br>"
        "Z: %{z:.2f} " + units[2] + " " + plotCoord + "<br>"
        "Lon: %{customdata[1]:.2f} " + co + "<br>" +
        "Lat: %{customdata[2]:.2f} " + co + "<br>" +
        "Alt: " + "{:.2f}".format(altkm) + " km<br>" +
        varn + ": %{customdata[0]:.4g} <br>" +
        "<extra></extra>"
    )
    if not axis:
        fig.update_scenes(xaxis=dict(visible=False),
                          yaxis=dict(visible=False),
                          zaxis=dict(visible=False))
    plotdt = dt.datetime.fromtimestamp(plotts, tz=pytz.utc)
    timestr = plotdt.strftime("%Y/%m/%d %H:%M:%S")
    lltxt = model + ',  ' + plotCoord + ' Coordinates,  '
    lltxt = lltxt + str(altkm) + ' km Altitude,  ' + timestr
    fig.update_layout(
        scene_aspectmode='data',
        title=dict(text=title,
                   yref="container", yanchor="top", x=0.01, y=0.95,
                   font=dict(size=16, family="sans serif", color="#000000")),
        annotations=[
            dict(text=lltxt, x=0.0, y=0.0, ax=0, ay=0, xanchor="left",
                 xshift=0, yshift=-20, xref="paper", yref="paper",
                 font=dict(size=16, family="sans serif", color="#000000"))
        ],
        margin=dict(l=0),
        width=750,
        height=450,
    )

    # Add poles to plot
    xt = [0., 0.]
    yt = [0., 0.]
    zt = [-1.2, 1.2]
    fig.add_scatter3d(mode='lines', x=xt, y=yt, z=zt,
                      line=dict(width=4, color='black'),
                      showlegend=False,
                      hovertemplate=plotCoord+' Pole<extra></extra>')
    if plotCoord != 'GEO':
        # Add pole for GEO as well
        xt = np.array([0., 0.])
        yt = np.array([0., 0.])
        zt = np.array([-1., 1.])
        tt = np.full(xt.shape, plotts)
        xt, yt, zt, un = ConvertCoord(tt, xt, yt, zt,
                                      'GEO', 'car', plotCoord, 'car')
        xt = 1.2*xt
        yt = 1.2*yt
        zt = 1.2*zt
        fig.add_scatter3d(mode='lines', x=xt, y=yt, z=zt,
                          line=dict(width=4, color='rgb(79,79,79)'),
                          showlegend=False,
                          hovertemplate='GEO Pole<extra></extra>')

    # Zero latitude
    xt = rscale*(np.cos(ilon*np.pi/180.))
    yt = rscale*(np.sin(ilon*np.pi/180.))
    zt = np.zeros(xt.shape)
    fig.add_scatter3d(mode='lines', x=xt, y=yt, z=zt,
                      line=dict(width=2, color='black'),
                      showlegend=False,
                      hovertemplate=plotCoord+' Latitude=0<extra></extra>')

    # Create blank earth surface
    elon = np.linspace(-180, 180, 181)
    elat = np.linspace(-90, 90, 91)
    elon_mg, elat_mg = np.meshgrid(np.array(elon), np.array(elat))
    ex = -(np.cos(elat_mg*np.pi/180.)*np.cos(elon_mg*np.pi/180.))
    ey = -(np.cos(elat_mg*np.pi/180.)*np.sin(elon_mg*np.pi/180.))
    ez = (np.sin(elat_mg*np.pi/180.))
    colorse = np.zeros(shape=ex.shape)
    # colorse[ex<0.]=1  # Option to make day side a lighter color
    colorscalee = ['rgb(99,99,99)', 'rgb(0,0,0)']
    fig.add_surface(x=ex, y=ey, z=ez, surfacecolor=colorse,
                    cmin=0, cmax=1, colorscale=colorscalee,
                    showlegend=False, showscale=False, hoverinfo='skip')

    # Shoreline (land/water boundaries)
    pos = shoreline(rscale=1.001, coord=plotCoord, utcts=plotts)
    fig.add_scatter3d(mode='lines', x=pos[:, 0], y=pos[:, 1], z=pos[:, 2],
                      line=dict(width=2, color='white'),
                      showlegend=False, hoverinfo='skip')

    return fig


def GDZSlice4D(interp, varname, model, date, plotType, plotCoord='GEO',
               fixed_time='', fixed_lon='', fixed_lat='', fixed_alt='',
               title='Plot Title', shoreline=False, colorscale='Viridis',
               crange='', logscale=False):
    """
    Function takes 4D (Time/Lon/Lat/Alt) Kamodo interpolator in GDZ coordinates
    and creates a 2D plot by fixing two dimensions to create a slice in the
    other two dimensions.
    This can be displayed in any supported coordinate system.

    Arguments:
      interp      Kamodo interpolator for variable to plot
      varname     String name of variable to plot
      model       String name of model
      date        Datetime object for date of plot
      plotType    String with two variables to plot, ie. 'Lon-Lat'
                    choose from Time, Lon, Lat, Alt

    Optional Arguments:
      plotCoord   Coordinate system to display the plot in
      fixed_time  Fixed time value if not a plot axis [hours]
      fixed_lon   Fixed longitude value if not a plot axis [degrees]
      fixed_lat   Fixed latgitude value if not a plot axis [degrees]
      fixed_alt   Fixed altitude value if not a plot axis [km]
      title       String used for plot title
      shoreline   Logical to show shoreline on map if Lon-Lat plot
      colorscale  String name of colorscale to use
      crange      Contour range to override min,max, ie. [999,2001]
      logscale    Logical to set the colorscale range to log scale

    Returs a plotly figure object.
    """
    import numpy as np
    import time
    import datetime as dt
    from kamodo import Kamodo
    import kamodo_ccmc.flythrough.model_wrapper as MW
    from kamodo_ccmc.flythrough.utils import ConvertCoord
    from kamodo_ccmc.tools.shoreline import shoreline

    tic = time.perf_counter()

    # Check passed in plotType and parse for use later
    Types = ['Lon', 'Lat', 'Alt', 'Time']
    axis = plotType.split('-')
    if len(axis) != 2:
        print("ERROR, invalid plotType passed.")
        return
    if axis[0] not in Types or axis[1] not in Types:
        print("ERROR, plotType elements invalid.")
        return

    # Set base timestamp to start of day passed in.
    basets = date.replace(hour=0, minute=0, second=0,
                          microsecond=0).timestamp()
    if 'Time' not in axis:
        slicets = basets + fixed_time*3600.

    # Get current coordinate system and type from model. Create label strings
    tmp = list(MW.Model_Variables(model, return_dict=True).values())[0][2:4]
    co = tmp[0]
    cot = tmp[1]
    costr = co + '(' + cot + ')'
    if co != 'GDZ' or cot != 'sph':
        print('Expecting coordinates to be GDZ(sph) but got ', costr)

    # Set itim array
    if 'Time' in axis:
        itim = np.linspace(0, 24, 97)
        if axis[0] == 'Time':
            i1 = itim
        if axis[1] == 'Time':
            i2 = itim
        fixed_time = ''
    else:
        if fixed_time == '':
            print("ERROR, fixed_time not passed.")
            return
        itim = np.array([fixed_time])

    # Set ilon array
    if 'Lon' in axis:
        ilon = np.linspace(-180, 180, 91)
        if axis[0] == 'Lon':
            i1 = ilon
        if axis[1] == 'Lon':
            i2 = ilon
        fixed_lon = ''
    else:
        if fixed_lon == '':
            print("ERROR, fixed_lon not passed.")
            return
        ilon = np.array([fixed_lon])

    # Set ilat array
    if 'Lat' in axis:
        ilat = np.linspace(-90, 90, 73)
        if axis[0] == 'Lat':
            i1 = ilat
        if axis[1] == 'Lat':
            i2 = ilat
        fixed_lat = ''
    else:
        if fixed_lat == '':
            print("ERROR, fixed_lat not passed.")
            return
        ilat = np.array([fixed_lat])

    # Set ialt array
    if 'Alt' in axis:
        # Using quick bisection solve to find alt range.
        dh, h = 50., 50.
        v = interp([12, 0, 0, h])
        for _ in range(100):
            if np.isnan(v):
                if dh < 0.:
                    dh = -.5 * dh
            else:
                if abs(dh) < 1e-6:
                    break
                if dh > 0.:
                    dh = -.5 * dh
            h += dh
            v = interp([12, 0, 0, h])
        h1 = h
        dh = 200.
        for _ in range(100):
            if np.isnan(v):
                if dh > 0.:
                    dh = -.5 * dh
            else:
                if abs(dh) < 1e-6:
                    break
                if dh < 0.:
                    dh = -.5 * dh
            h += dh
            v = interp([12, 0, 0, h])
        ialt = np.linspace(h1, h, 101)
        if axis[0] == 'Alt':
            i1 = ialt
        if axis[1] == 'Alt':
            i2 = ialt
        fixed_alt = ''
    else:
        if fixed_alt == '':
            print("ERROR, fixed_alt not passed.")
            return
        ialt = np.array([fixed_alt])

    # Set grid to interpolate to
    i1mg, i2mg = np.meshgrid(np.array(i1), np.array(i2))
    w = np.reshape(i1mg, -1)  # Dummy to store shape for next steps

    if axis[0] == 'Lon':
        x = np.reshape(i1mg, -1)
    elif axis[1] == 'Lon':
        x = np.reshape(i2mg, -1)
    else:
        x = np.full(w.shape, fixed_lon)

    if axis[0] == 'Lat':
        y = np.reshape(i1mg, -1)
    elif axis[1] == 'Lat':
        y = np.reshape(i2mg, -1)
    else:
        y = np.full(w.shape, fixed_lat)

    z = np.ones(w.shape)  # Either 1km or 1Re are fine for coordinate transform

    # The t array needs to be timestamps
    if axis[0] == 'Time':
        t = np.reshape(i1mg, -1)
        t = basets + t*3600.
    elif axis[1] == 'Time':
        t = np.reshape(i2mg, -1)
        t = basets + t*3600.
    else:
        t = np.full(w.shape, slicets)

    # Convert from plotCoord coordinates into model coordinates to interpolate
    if plotCoord == 'GDZ':
        # No change to coordinates, so don't transform
        xx, yy, zz = x, y, z
    elif plotCoord == 'GEO':
        # Only altitude changes, and that is reset later
        xx, yy, zz = x, y, z
    else:
        # GDZ only converts to/from GEO, so  have to convert twice
        xx, yy, zz, units = ConvertCoord(t, x, y, z, plotCoord,
                                         'sph', 'GEO', 'sph')
        x, y, z = xx, yy, zz
        xx, yy, zz, units = ConvertCoord(t, x, y, z,
                                         'GEO', 'sph', 'GDZ', 'sph')

    # Don't use transformed zz, reset to desired altitude (km) to interpolate
    if axis[0] == 'Alt':
        zz = np.reshape(i1mg, -1)
    elif axis[1] == 'Alt':
        zz = np.reshape(i2mg, -1)
    else:
        zz = np.full(w.shape, fixed_alt)

    # The tt array needs to be in hours
    if axis[0] == 'Time':
        tt = np.reshape(i1mg, -1)
    elif axis[1] == 'Time':
        tt = np.reshape(i2mg, -1)
    else:
        tt = np.full(w.shape, fixed_time)

    positions = np.transpose([tt, xx, yy, zz])
    val = interp(positions)

    if logscale:
        val[val <= 0.] = np.nan
        val = np.log10(val)
        varname = 'log10<br>'+varname

    val2 = np.reshape(val, i1mg.shape)

    cmin = np.min(val)
    cmax = np.max(val)
    if crange != '':
        cmin = float(crange[0])
        cmax = float(crange[1])

    def plot_2D(v1=i1, v2=i2):
        return val2
    plot1 = Kamodo(plot_2D=plot_2D)
    fig = plot1.plot(plot_2D=dict())

    # Set colorscale
    fig = figMods(fig, colorscale=colorscale, ncont=200)

    # Set plot options
    fig.update_traces(
        zmin=cmin, zmax=cmax, colorbar=dict(title=varname, tickformat=".3g"),
        hovertemplate=str(axis[0]) + ": %{x:.2f}<br>" +
        str(axis[1]) + ": %{y:.2f}<br>" +
        varname + ": %{z:.6g}<br>" + "<extra></extra>"
    )
    fig.update_xaxes(title=axis[0])
    fig.update_yaxes(title=axis[1])
    if axis[0] == 'Lon' or axis[0] == 'Lat':
        fig.update_xaxes(tick0=0., dtick=45.)
    if axis[0] == 'Time':
        fig.update_xaxes(tick0=0., dtick=3.)
    if axis[1] == 'Lon' or axis[1] == 'Lat':
        fig.update_yaxes(tick0=0., dtick=45.)
    if axis[1] == 'Time':
        fig.update_yaxes(tick0=0., dtick=3.)
    timestr = date.strftime("%Y/%m/%d")
    # Check if fixed_time is >= 24 hours, if so adjust timestr and fixed_time
    if fixed_time != "":
        date2 = date
        while fixed_time >= 24.:
            fixed_time += -24.
            date2 = date2 + dt.timedelta(days=1)
            timestr = date2.strftime("%Y/%m/%d")
    subt = model + ',  ' + plotCoord + ' Coordinates,  '
    subt += timestr + ',  Slice at:'
    if fixed_time != "":
        subt = subt + '  ' + str(fixed_time) + ' hrs'
    if fixed_alt != "":
        subt = subt + '  ' + str(fixed_alt) + ' km Altitude'
    if fixed_lon != "":
        subt = subt + '  ' + str(fixed_lon) + ' deg Longitude'
    if fixed_lat != "":
        subt = subt + '  ' + str(fixed_lat) + ' deg Latitude'
    fig.update_layout(
        plot_bgcolor="white", scene_aspectmode='data',
        title=dict(text=title + '<br>' + subt,
                   yref="container", yanchor="top", x=0.01, y=0.97,
                   font=dict(size=16, family="sans serif", color="#000000")),
        margin=dict(l=0), width=750, height=375)

    if plotType == 'Lon-Lat':
        # Shoreline (land/water boundaries)
        if plotCoord == 'GDZ':
            # Convert into GEO as we only use lon/lat
            pos = shoreline(coord='GEO', coordT='sph', utcts=slicets)
        else:
            pos = shoreline(coord=plotCoord, coordT='sph', utcts=slicets)
        fig.add_scattergl(mode='lines', x=pos[:, 0], y=pos[:, 1],
                          hoverinfo='skip', showlegend=False,
                          line=dict(width=1, color='#EEEEEE'))

    toc = time.perf_counter()
    print(f"Time: {toc - tic:0.4f} seconds")
    return fig


def swmfgm3D(ko, var, time=0., title='',
             xcut=-999., ycut=-999., zcut=-999.,
             zgrid=False, znodes=False):
    '''
    Function to combine default X,Y,Z constant slices into one 3D
      plot. Returns a full 3D plotly figure.

    ko:      Kamodo object
    var:     string variable name
    time:    floating time in hours from start of first day of data
    title:   string text for plot title top left
    xcut:    if set, value for X slice through data
    ycut:    if set, value for Y slice through data
    zcut:    if set, value for Z slice through data
    zgrid:   logical to return plot of grid dz size
    znodes:  logical to show dots at grid locations
    '''
    import re
    import numpy as np
    import pandas as pd
    from kamodo import Kamodo
    import plotly.graph_objs as go

    useX = False
    useY = False
    useZ = False
    if xcut > -999.:
        useX = True
    if ycut > -999.:
        useY = True
    if zcut > -999.:
        useZ = True

    cmin = 1.e+99
    cmax = -1.e+99
    opacity = 1.

    if useZ:
        print('Extracting z slice at', zcut)
        figZ = ko.plot(var, plot_partial={var: {'time': time, 'Z': zcut}})
        x0 = figZ.data[0]['x']
        y0 = figZ.data[0]['y']
        z0 = np.array([zcut])
        val0 = figZ.data[0]['z']
        cmin = min(cmin, np.min(val0))
        cmax = max(cmax, np.max(val0))

    if useY:
        print('Extracting y slice at', ycut)
        figY = ko.plot(var, plot_partial={var: {'time': time, 'Y': ycut}})
        x1 = figY.data[0]['x']
        z1 = figY.data[0]['y']
        y1 = np.array([ycut])
        val1 = figY.data[0]['z']
        cmin = min(cmin, np.min(val1))
        cmax = max(cmax, np.max(val1))

    if useX:
        print('Extracting x slice at', xcut)
        figX = ko.plot(var, plot_partial={var: {'time': time, 'X': xcut}})
        y2 = figX.data[0]['x']
        z2 = figX.data[0]['y']
        x2 = np.array([xcut])
        val2 = figX.data[0]['z']
        cmin = min(cmin, np.min(val2))
        cmax = max(cmax, np.max(val2))

    print('starting plot')
    fig = go.Figure()
    if useZ:
        if zgrid:
            x_mg, y_mg, z_mg = np.meshgrid(x0, y0, z0)
            x_mg_1d = x_mg.reshape(-1)
            y_mg_1d = y_mg.reshape(-1)
            z_mg_1d = z_mg.reshape(-1)
            iv = []
            jv = []
            kv = []
            dx = []
            for iy in range(len(y0)-1):
                for ix in range(len(x0)-1):
                    # For each cell, create two triangles
                    iv.append(ix + iy*len(x0))
                    jv.append(ix+1 + iy*len(x0))
                    kv.append(ix + (iy+1)*len(x0))
                    dx.append(max(x_mg_1d[iv[-1]],
                                  x_mg_1d[jv[-1]],
                                  x_mg_1d[kv[-1]]) -
                              min(x_mg_1d[iv[-1]],
                                  x_mg_1d[jv[-1]],
                                  x_mg_1d[kv[-1]]))
                    iv.append(ix + (iy+1)*len(x0))
                    jv.append(ix+1 + (iy+1)*len(x0))
                    kv.append(ix+1 + iy*len(x0))
                    dx.append(max(x_mg_1d[iv[-1]],
                                  x_mg_1d[jv[-1]],
                                  x_mg_1d[kv[-1]]) -
                              min(x_mg_1d[iv[-1]],
                                  x_mg_1d[jv[-1]],
                                  x_mg_1d[kv[-1]]))
            dxlog2 = np.log2(dx)
            figG = go.Figure(data=[
                go.Mesh3d(
                    x=x_mg_1d,
                    y=y_mg_1d,
                    z=z_mg_1d,
                    colorbar_title='log2(dx)',
                    colorscale='Viridis',
                    intensity=dxlog2,
                    intensitymode='cell',
                    # intensity=x_mg_1d,  # If using, fix hover and labeling
                    # intensitymode='vertex',
                    i=iv,
                    j=jv,
                    k=kv,
                    name='y',
                    showscale=True
                )
            ])
            figG.update_traces(
                flatshading=True,
                customdata=dx,
                hovertemplate="<b>Grid Size</b><br>" +
                # "i: %{i:9i}<br>" +
                # "j: %{j:9i}<br>" +
                # "k: %{k:9i}<br>" +
                "log2(dx): %{intensity:.4g}<br>" +
                "dx: %{customdata:.4g}<br>" +
                "<extra></extra>"
            )
            if znodes:
                figG.add_scatter3d(x=x_mg_1d, y=y_mg_1d, z=z_mg_1d,
                                   mode='markers',
                                   marker=dict(size=1, color='black'),
                                   line=dict(width=1))

            return figG

        val0b = np.reshape(val0, (len(y0), len(x0), len(z0)))
        xvar = re.sub(' +', ' ',
                      figZ.layout.xaxis.title.text.strip('$').strip())
        yvar = re.sub(' +', ' ',
                      figZ.layout.yaxis.title.text.strip('$').strip())
        zvar = xvar.replace('X', 'Z')
        cvar = re.sub(' +', ' ',
                      figZ.data[0]['colorbar'].title.text.strip('$').strip())

        def plot3d_var0(X=x0, Y=y0, Z=z0):
            return val0b

        ko0 = Kamodo(plot_var=plot3d_var0)
        figZ2 = ko0.plot(plot_var=dict())
        figZ2.update_traces(colorscale='Viridis')
        figZ2.update_traces(
            cmin=cmin, cmax=cmax,
            colorbar=dict(title=cvar, tickformat=".3g"),
            opacity=opacity,
            customdata=np.transpose(val0b, axes=[1, 0, 2]),
            hovertemplate="<b>Position Values</b><br>" +
            xvar + ": %{x:.3f}<br>" +
            yvar + ": %{y:.3f}<br>" +
            zvar + ": %{z:.3f}<br>" +
            cvar + ": %{customdata:.4g}<br>" +
            "<extra></extra>"
        )
        fig.add_trace(figZ2.data[0])

    if useY:
        xm, ym, zm = np.meshgrid(x1, y1, z1)
        val1tp = np.transpose(val1, axes=[1, 0])
        val1b = np.reshape(val1tp, (len(y1), len(x1), len(z1)))
        xvar = re.sub(' +', ' ',
                      figY.layout.xaxis.title.text.strip('$').strip())
        zvar = re.sub(' +', ' ',
                      figY.layout.yaxis.title.text.strip('$').strip())
        yvar = xvar.replace('X', 'Y')
        cvar = re.sub(' +', ' ',
                      figY.data[0]['colorbar'].title.text.strip('$').strip())

        def plot3d_var1(X=x1, Y=y1, Z=z1):
            return val1b

        ko1 = Kamodo(plot_var=plot3d_var1)
        figY2 = ko1.plot(plot_var=dict())
        figY2.update_traces(colorscale='Viridis')
        figY2.update_traces(
            cmin=cmin, cmax=cmax,
            colorbar=dict(title=cvar, tickformat=".3g"),
            opacity=opacity,
            customdata=np.transpose(val1b, axes=[2, 1, 0]),
            hovertemplate="<b>Position Values</b><br>" +
            xvar + ": %{x:.3f} <br>" +
            yvar + ": %{y:.3f} <br>" +
            zvar+": %{z:.3f} <br>" +
            cvar+": %{customdata:.4g} <br>" +
            "<extra></extra>"
        )
        fig.add_trace(figY2.data[0])

    if useX:
        xm, ym, zm = np.meshgrid(x2, y2, z2)
        val2tp = np.transpose(val2, axes=[1, 0])
        val2b = np.reshape(val2tp, (len(y2), len(x2), len(z2)))
        yvar = re.sub(' +', ' ',
                      figX.layout.xaxis.title.text.strip('$').strip())
        zvar = re.sub(' +', ' ',
                      figX.layout.yaxis.title.text.strip('$').strip())
        xvar = yvar.replace('Y', 'Z')
        cvar = re.sub(' +', ' ',
                      figX.data[0]['colorbar'].title.text.strip('$').strip())

        def plot3d_var2(X=x2, Y=y2, Z=z2):
            return val2b

        ko2 = Kamodo(plot_var=plot3d_var2)
        figX2 = ko2.plot(plot_var=dict())
        figX2.update_traces(colorscale='Viridis')
        figX2.update_traces(
            cmin=cmin, cmax=cmax,
            colorbar=dict(title=cvar, tickformat=".3g"),
            opacity=opacity,
            customdata=np.transpose(val2b, axes=[2, 0, 1]),
            hovertemplate="<b>Position Values</b><br>" +
            xvar + ": %{x:.3f} <br>" + yvar + ": %{y:.3f} <br>" +
            zvar + ": %{z:.3f} <br>" + cvar + ": %{customdata:.4g} <br>" +
            "<extra></extra>"
        )
        fig.add_trace(figX2.data[0])

    # Create blank inner boundary
    elon = np.linspace(-180, 180, 181)
    elat = np.linspace(-90, 90, 91)
    elon_mg, elat_mg = np.meshgrid(np.array(elon), np.array(elat))
    ex = -2.5*(np.cos(elat_mg*np.pi/180.)*np.cos(elon_mg*np.pi/180.))
    ey = -2.5*(np.cos(elat_mg*np.pi/180.)*np.sin(elon_mg*np.pi/180.))
    ez = 2.5*(np.sin(elat_mg*np.pi/180.))
    colorse = np.zeros(shape=ex.shape)
    # colorse[ex<0.]=1  # Option to make day side a lighter color
    colorscalee = ['rgb(99,99,99)', 'rgb(0,0,0)']
    fig.add_surface(x=ex, y=ey, z=ez, surfacecolor=colorse,
                    cmin=0, cmax=1, colorscale=colorscalee,
                    showlegend=False, showscale=False, hoverinfo='skip')

    camera = dict(
        center=dict(x=0.35, y=0., z=0.),
    )
    fig.update_layout(
        scene_camera=camera,
        scene_aspectmode='data',
        title=dict(text=title,
                   yref="container", yanchor="top", x=0.01, y=0.97,
                   font=dict(size=16, family="sans serif", color="#000000"))
    )

    return fig


def swmfgm3Darb(ko, var, time=0., pos=[0, 0, 0], normal=[0, 1, 0],
                title='', lowerlabel='', showgrid=False, showibs=False,
                log10=False):
    '''
    Function to create a slice for any point and normal and interpolate
      from Kamodo object onto that grid. Returns a full 3D plotly figure.

    ko:    Kamodo object
    var:   string variable name
    time:  floating time in hours from start of first day of data
    pos, normal:   position and normal vector for slice plane
    title:         string text for plot title top left
    lowerlabel:    string text to label plot lower left
    showgrid:      logical to show dots at grid locations
    showibs:       logical to show inner boundary sphere
    log10:         logical to take log10 of display value
    '''
    import numpy as np
    import plotly.graph_objs as go
    import kamodo_ccmc.flythrough.model_wrapper as MW

    # Set interpolator and variable labels
    interp = getattr(ko, var)
    varu = ko.variables[var]['units']
    varlabel = var+" ["+varu+"]"

    # Compute values from pos, normal values
    uvec = normal/np.linalg.norm(normal)  # unit normal vector
    odist = np.dot(uvec, pos)  # closest distance to slice from origin
    opos = odist*uvec  # vector from origin to closest point

    # Compute base grid
    dg = np.linspace(-180, 180, 145)  # degree grid
    rv = 0.  # radius value
    if abs(odist) < 2.5:
        rv = 2.5*np.sin(((2.5-odist)/2.5)*np.pi/2.)
    rg = []  # radius grid
    rg.append(rv)
    rtrans = 2.5  # radius transition point from fixed dr to variable
    for _ in range(250):
        if rv < rtrans:
            rv += .11
        else:
            rv += .11*rv/rtrans
        if rv > 300.:
            break
        rg.append(rv)
    rg = np.array(rg)
    dg_mg, rg_mg = np.meshgrid(dg, rg)
    gx = rg_mg*np.cos(dg_mg*np.pi/180.)
    gy = rg_mg*np.sin(dg_mg*np.pi/180.)
    gx_1d = gx.reshape(-1)
    gy_1d = gy.reshape(-1)
    gz_1d = np.zeros([len(gx_1d)])
    time_1d = np.full((len(gx_1d)), time)
    grid0 = np.stack((gx_1d, gy_1d, gz_1d), axis=-1)  # nx3 position grid

    # Transform base grid to pos/normal and trim if needed
    # Rotate
    if abs(uvec[2]) < 1.:  # No rotation for +/- Z normal
        new_xaxis = np.cross([0, 0, 1], uvec)
        new_yaxis = np.cross(new_xaxis, uvec)
        transform = np.array([new_xaxis, new_yaxis, uvec]).T
        grid0 = np.inner(grid0, transform)
    # Shift
    grid0[:, 0] += opos[0]
    grid0[:, 1] += opos[1]
    grid0[:, 2] += opos[2]
    # back to 2D grid
    gx2 = grid0[:, 0].reshape(len(rg), -1)
    gy2 = grid0[:, 1].reshape(len(rg), -1)
    gz2 = grid0[:, 2].reshape(len(rg), -1)
    # Trim points beyond data, first extracting range
    cr = MW.Coord_Range(ko, [var], return_dict=True, print_output=False)
    x1, x2 = cr[var]['X'][0], cr[var]['X'][1]
    y1, y2 = cr[var]['Y'][0], cr[var]['Y'][1]
    z1, z2 = cr[var]['Z'][0], cr[var]['Z'][1]
    for j in range(len(dg)):
        for i in range(len(rg)):
            found = False
            frac = 0.
            if gx2[i, j] < x1:
                found = True
                frac = max(frac, (gx2[i, j] - x1)/(gx2[i, j] - gx2[i-1, j]))
            if gx2[i, j] > x2:
                found = True
                frac = max(frac, (gx2[i, j] - x2)/(gx2[i, j] - gx2[i-1, j]))
            if gy2[i, j] < y1:
                found = True
                frac = max(frac, (gy2[i, j] - y1)/(gy2[i, j] - gy2[i-1, j]))
            if gy2[i, j] > y2:
                found = True
                frac = max(frac, (gy2[i, j] - y2)/(gy2[i, j] - gy2[i-1, j]))
            if gz2[i, j] < z1:
                found = True
                frac = max(frac, (gz2[i, j] - z1)/(gz2[i, j] - gz2[i-1, j]))
            if gz2[i, j] > z2:
                found = True
                frac = max(frac, (gz2[i, j] - z2)/(gz2[i, j] - gz2[i-1, j]))
            if found:
                gx2[i, j] -= frac*(gx2[i, j] - gx2[i-1, j])
                gy2[i, j] -= frac*(gy2[i, j] - gy2[i-1, j])
                gz2[i, j] -= frac*(gz2[i, j] - gz2[i-1, j])
                gx2[i+1:, j] = gx2[i, j]
                gy2[i+1:, j] = gy2[i, j]
                gz2[i+1:, j] = gz2[i, j]
                break

    # Create 4D (nx4) grid and interpolate plot values
    grid = np.ndarray(shape=(len(gx_1d), 4), dtype=np.float32)
    grid[:, 0] = time_1d
    grid[:, 1] = gx2.reshape(-1)
    grid[:, 2] = gy2.reshape(-1)
    grid[:, 3] = gz2.reshape(-1)
    value = interp(grid)
    if log10:
        value[value <= 0.] = np.nan
        value = np.log10(value)
        varlabel = "log10("+varlabel+")"

    # Build connectivity grid cell by cell looping over positions
    iv = []
    jv = []
    kv = []
    for iy in range(len(rg)-1):
        for ix in range(len(dg)-1):
            # For each cell, create two triangular connectivity entries
            iv.append(ix+iy*len(dg))
            jv.append(ix+1+iy*len(dg))
            kv.append(ix+(iy+1)*len(dg))

            iv.append(ix+(iy+1)*len(dg))
            jv.append(ix+1+(iy+1)*len(dg))
            kv.append(ix+1+iy*len(dg))

    # Build resulting plot
    fig = go.Figure(data=[
        go.Mesh3d(
            x=grid[:, 1], y=grid[:, 2], z=grid[:, 3], i=iv, j=jv, k=kv,
            colorbar_title=varlabel, colorscale='Viridis',
            intensity=value, intensitymode='vertex',
            name='cont', showscale=True
        )
    ])
    fig.update_traces(
        flatshading=True,
        hovertemplate="X: %{x:.4g}<br>" + "Y: %{y:.4g}<br>" +
        "Z: %{z:.4g}<br>" + varlabel + ": %{intensity:.4g}<br>" +
        "<extra></extra>"
    )
    # Add grid points
    if showgrid:
        fig.add_scatter3d(
            x=grid[:, 1], y=grid[:, 2], z=grid[:, 3], mode='markers',
            marker=dict(size=1, color='white'), line=dict(width=1),
            name='grid'
        )

    # Create blank inner boundary
    if showibs:
        elon = np.linspace(-180, 180, 181)
        elat = np.linspace(-90, 90, 91)
        elon_mg, elat_mg = np.meshgrid(np.array(elon), np.array(elat))
        ex = -2.5*(np.cos(elat_mg*np.pi/180.)*np.cos(elon_mg*np.pi/180.))
        ey = -2.5*(np.cos(elat_mg*np.pi/180.)*np.sin(elon_mg*np.pi/180.))
        ez = +2.5*(np.sin(elat_mg*np.pi/180.))
        colorse = np.zeros(shape=ex.shape)
        # colorse[ex<0.]=1  # Option to make day side a lighter color
        colorscalee = ['rgb(99,99,99)', 'rgb(0,0,0)']
        fig.add_surface(x=ex, y=ey, z=ez, surfacecolor=colorse,
                        cmin=0, cmax=1, colorscale=colorscalee,
                        showlegend=False, showscale=False, hoverinfo='skip')

    fig.update_layout(
        scene_aspectmode='data',
        title=dict(text=title,
                   yref="container", yanchor="top", x=0.01, y=0.97,
                   font=dict(size=16, family="sans serif", color="#000000")),
        scene=dict(
            xaxis=dict(showbackground=False, showgrid=False),
            yaxis=dict(showbackground=False, showgrid=False),
            zaxis=dict(showbackground=False, showgrid=False)),
        annotations=[
            dict(text=lowerlabel, x=0.0, y=0.0, ax=0, ay=0, xanchor="left",
                 xshift=0, yshift=-20, xref="paper", yref="paper",
                 font=dict(size=16, family="sans serif", color="#000000"))
        ],
        margin=dict(l=0, t=35),
    )

    return fig
