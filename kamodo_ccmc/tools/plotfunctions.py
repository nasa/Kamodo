'''
Plot Functions:
'''
colorset = ['#00B7FF', '#004DFF', '#00FFFF', '#826400', '#580041', '#FF00FF', 
            '#00FF00', '#C500FF', '#B4FFD7', '#FFCA00', '#969600', '#B4A2FF',
            '#C20078', '#0000C1', '#FF8B00', '#FFC8FF', '#666666', '#FF0000',
            '#CCCCCC', '#009E8F', '#D7A870', '#8200FF', '#960000', '#006F00']

### ====================================================================================== ###
def rotatePoints(x, y, degA):
    '''
    Rotate arrays of x,y values around 0,0 by degA given in degrees counter-clockwise
    Return new x,y values
    '''
    import numpy as np

    cosA = np.cos(degA*np.pi/180.)
    sinA = np.sin(degA*np.pi/180.)

    newx = cosA * x - sinA * y
    newy = sinA * x + cosA * y
    
    return newx, newy


### ====================================================================================== ###
def figMods(fig, log10=False, lockAR=False, ncont=-1, colorscale='',
            cutInside=-1., returnGrid=False, enhanceHover=False,
            enhanceHover1D=False, newTitle='', cText='',
            llText='', llText2='', coText='', resetAnnotations=True,
            zcrange=False, crange='', xtic='', ytic=''):
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
      returnGrid  Take the plot and return a new grid only plotly object (beta)
      enhanceHover Logical, if true update the hover information
      newTitle    String to use for new plot title
      cText       String to use for colorbar label
      llText      String to add to lower left of plot
      llText2     String to add just above llText
      coText      String to add coordinate system to lower right of plot
      resetAnnotations Logical, if true will clear annotations if
                  llText, llText2, or coText are not empty
      crange      Two position array with min/max contour values, [cmin,cmax]
      xtic        X axis tick spacing
      ytic        Y axis tick spacing
    '''

    import re
    import numpy as np
    import plotly.graph_objects as go

    # Create list of the plot types in the figure object
    #   Some figure mods can only apply to some traces in the fig
    typelist = []
    for i in range(len(fig.data)):
        typelist.append(fig.data[i].type)
    typelist = np.unique(typelist)

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
                          contours=dict(coloring="fill", showlines=False),
                          selector=dict(type='contour'))

    if colorscale != '':
        usedcolorscale = colorscale
        if colorscale == "BlueRed":
            usedcolorscale = "RdBu_r"
        elif colorscale == "Rainbow":
            usedcolorscale=[[0.00, 'rgb(0,0,255)'],
                            [0.25, 'rgb(0,255,255)'],
                            [0.50, 'rgb(0,255,0)'],
                            [0.75, 'rgb(255,255,0)'],
                            [1.00, 'rgb(255,0,0)']]
        fig.update_traces(colorscale=usedcolorscale, selector=dict(type='contour'))
        fig.update_traces(colorscale=usedcolorscale, selector=dict(type='surface'))
        fig.update_traces(colorscale=usedcolorscale, selector=dict(type='mesh3d'))

    if crange != '':
        cmin = float(crange[0])
        cmax = float(crange[1])
        for val in typelist:
            if val in ['contour']:
                fig.update_traces(zmin=cmin, zmax=cmax, selector={'type': val})
            if val in ['surface', 'mesh3d']:
                fig.update_traces(cmin=cmin, cmax=cmax, selector={'type': val})

    if zcrange:  # ONLY for data block 0
        val = fig.data[0].type
        if val in ['contour']:
            if fig.data[0].zmin and fig.data[0].zmax:
                zmin = fig.data[0].zmin
                zmax = fig.data[0].zmax
            else:
                zmin = np.min(fig.data[0].z)
                zmax = np.max(fig.data[0].z)
            zmm = max(abs(zmin),abs(zmax))
            fig.update_traces(zmin=-zmm, zmax=zmm, selector={'type': val})
        if val in ['surface', 'mesh3d']:
            if fig.data[0].cmin and fig.data[0].cmax:
                cmin = fig.data[0].cmin
                cmax = fig.data[0].cmax
            else:
                if val in ['surface']:
                    cmin = np.min(fig.data[0].surfacecolor)
                    cmax = np.max(fig.data[0].surfacecolor)
                else:
                    cmin = np.min(fig.data[0].intensity)
                    cmax = np.max(fig.data[0].intensity)
            cmm = max(abs(cmin),abs(cmax))
            fig.update_traces(cmin=-cmm, cmax=cmm, selector={'type': val})

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

    if enhanceHover1D:
        xvar = re.sub(' +', ' ',
                      fig.layout.xaxis.title.text.strip('$').strip())
        yvar = re.sub(' +', ' ',
                      fig.layout.yaxis.title.text.strip('$').strip())
        fig.update_traces(
            hovertemplate=xvar + ": %{x:.3f} <br>" +
            yvar + ": %{y:.4g} <br>" +
            "<extra></extra>"
        )

    if newTitle != '':
        fig.layout.title.text = newTitle

    if cText != '':
        fig.update_traces(colorbar=dict(title=cText, tickformat=".3g"))

    if llText != '' or llText2 != '' or coText != '':
        if resetAnnotations:
            fig.layout.annotations = []

        # Default text offsets unless specific plot types
        xs1, xs2 = -70,  -6
        ys1, ys2 = -44, -28
        if 'surface' in typelist or 'mesh3d' in typelist:
            xs1, xs2 =   6, -6
            ys1, ys2 = -23, -8
        if fig.data[0].name:
            if '2Dpolar' in fig['data'][0]['name']:
                xs1, xs2 = -15,  40
                ys1, ys2 = -36, -20
        if llText != '':
            fig.add_annotation(text=llText, x=0.0, y=0.0, ax=0, ay=0, xanchor="left",
                     xshift=xs1, yshift=ys1, xref="paper", yref="paper",
                     font=dict(size=12, family="sans serif", color="#000000"))
        if llText2 != '':
            fig.add_annotation(text=llText2, x=0.0, y=0.0, ax=0, ay=0, xanchor="left",
                     xshift=xs1, yshift=ys2, xref="paper", yref="paper",
                     font=dict(size=12, family="sans serif", color="#000000"))
        if coText != '':
            fig.add_annotation(text=coText, x=1.0, y=0.0, ax=0, ay=0, xanchor="right",
                     xshift=xs2, yshift=ys1, xref="paper", yref="paper",
                     font=dict(size=12, family="sans serif", color="#000000"))

    if xtic != '':
        fig.update_layout(xaxis = dict(tick0 = 0., dtick = xtic))

    if ytic != '':
        fig.update_layout(yaxis = dict(tick0 = 0., dtick = ytic))

    return fig


### ====================================================================================== ###
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


### ====================================================================================== ###
def toLog10(fig):
    """
    Placeholder function for compatibility with old notebooks.
    Use figMods instead.

    Arguments:
      fig         A plotly figure object

    Returns a modified plotly figure object
    """

    return figMods(fig, log10=True)


### ====================================================================================== ###
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


### ====================================================================================== ###
def ReplotLL3D(figIn, model, altkm, plotts, plotCoord='GEO',
               title='Plot Title', colorscale='Viridis', crange='',
               opacity=0.70, axis=True, debug=0, zcrange=False,
               useCo='', useCot='', plotType='surface',
               showPole=True, showEarth=True,
               showshore=True, shorewidth=2, shorecolor='white'):
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
    import plotly.graph_objs as go
    from scipy.interpolate import griddata
    from plotly.subplots import make_subplots

    if plotCoord == 'GDZ':
        plotCoord = 'GEO';
    if useCo != '' and useCot != '':
        co = useCo
        cot = useCot
    else:
        # Get current coordinate system and type. Create string for labels
        tmp = list(MW.Model_Variables(model, return_dict=True).values())[0][2:4]
        co = tmp[0]
        cot = tmp[1]
    if cot != 'sph':
        print('ERROR, coordinate is not spherical! Returning ...')
        return

    # Pull out lon, lat, values and min/max from passed in figure
    lon = figIn.data[0]['x']
    lat = figIn.data[0]['y']
    val = figIn.data[0]['z']
    if len(lon) == len(lat):
        # Pad duplicate values to end of arrays to avoid array sizes bug
        lat = np.append(lat, lat[-1])
        aaa = np.array([val[-1,:]])
        val = np.append(val, aaa, axis=0)
    val2 = np.reshape(val, (len(lat), len(lon)))
    val1 = np.reshape(val2, -1)
    varn = figIn.data[0]['colorbar']['title']['text']
    varn2 = varn
    pos = varn.find('_ijk')
    if pos != -1:
        varn2 = varn[0:pos]+varn[pos+4:]
    cmin = np.min(val)
    cmax = np.max(val)
    if crange != '':
        cmin = float(crange[0])
        cmax = float(crange[1])
    if zcrange:  # Set zero value to center of contour range
        cmm = max(abs(cmin),abs(cmax))
        cmin = -cmm
        cmax =  cmm

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

    # Convert incoming coordinates into plot coordinates (cartesian)
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
    if co != plotCoord:
        lon4, lat4, alt4, units = ConvertCoord(full_t, full_x, full_y, full_z, co, cot, plotCoord, 'sph')
    x = np.reshape(xx, (len(lat), len(lon)))
    y = np.reshape(yy, (len(lat), len(lon)))
    z = np.reshape(zz, (len(lat), len(lon)))
    # NaNs in original plot don't hide properly, so set positions to NaN
    mask2 = np.isnan(val2)
    x[mask2] = np.nan
    y[mask2] = np.nan
    z[mask2] = np.nan

    # \\\ Start of 2Dpolar
    #   -this plot will return and other plot types will continue
    if plotType == '2Dpolar':
        # verified rotations for GEO, GSM, GSE, SM, GEI
        degRotate = 90.
        if plotCoord == 'GEO':
            tmpHr = pytz.utc.localize(dt.datetime.utcfromtimestamp(plotts)).hour
            degRotate = (360.*tmpHr/24.) - 90.

        # Constants for lat/lon lines on plot
        c80=rscale*np.sin(((90.-80.)/90.)*np.pi/2.)
        c70=rscale*np.sin(((90.-70.)/90.)*np.pi/2.)
        c60=rscale*np.sin(((90.-60.)/90.)*np.pi/2.)
        c50=rscale*np.sin(((90.-50.)/90.)*np.pi/2.)
        c50d=c50/np.sqrt(2.)

        # Original coordinates lon,lat
        full_lon = np.reshape(lon_mg, -1)
        full_lat = np.reshape(lat_mg, -1)

        # Rotated points
        xxR, yyR = rotatePoints(xx, yy, degRotate)

        # Polar plot base grid
        x3 = np.linspace(-.65, .65, 131)
        y3 = np.linspace(-.65, .65, 131)
        xx3, yy3 = np.meshgrid(x3, y3)

        # North
        mask = zz > 0.
        xxN = xxR[mask]
        yyN = yyR[mask]
        vvN = val1[mask]
        lonN = full_lon[mask]
        latN = full_lat[mask]
        z3N = griddata((xxN, yyN), vvN, (xx3, yy3), method='linear')
        lon3N = griddata((xxN, yyN), lonN, (xx3, yy3), method='linear')
        lat3N = griddata((xxN, yyN), latN, (xx3, yy3), method='linear')
        figN = go.Figure(go.Contour(x=x3, y=y3, z=z3N, showscale=False, name='2DpolarN'))
        if co == plotCoord:
            figN.update_traces(zmin=cmin, zmax=cmax,
                colorbar=dict(title=varn2, tickformat=".3g"),
                customdata=np.dstack((z3N, lon3N, lat3N)),
                hovertemplate="<b>" + model + " <> NORTH</b><br>" +
                    "Model Lon: %{customdata[1]:.2f} " + co + "<br>" +
                    "Model Lat: %{customdata[2]:.2f} " + co + "<br>" +
                    "Alt: " + "{:.2f}".format(altkm) + " km<br>" +
                    varn2 + ": %{customdata[0]:.4g} <br>" + 
                    "<extra></extra>"
                              )
        else:
            lonNplot = lon4[mask]
            latNplot = lat4[mask]
            plotlon3N = griddata((xxN, yyN), lonNplot, (xx3, yy3), method='linear')
            plotlat3N = griddata((xxN, yyN), latNplot, (xx3, yy3), method='linear')
            figN.update_traces(zmin=cmin, zmax=cmax,
                colorbar=dict(title=varn2, tickformat=".3g"),
                customdata=np.dstack((z3N, lon3N, lat3N, plotlon3N, plotlat3N)),
                hovertemplate="<b>" + model + " <> NORTH</b><br>" +
                    "Model Lon: %{customdata[1]:.2f} " + co + "<br>" +
                    "Model Lat: %{customdata[2]:.2f} " + co + "<br>" +
                    "Alt: " + "{:.2f}".format(altkm) + " km<br>" +
                    "Plot Lon: %{customdata[3]:.2f} " + plotCoord + "<br>" +
                    "Plot Lat: %{customdata[4]:.2f} " + plotCoord + "<br>" +
                    varn2 + ": %{customdata[0]:.4g} <br>" + 
                    "<extra></extra>"
                              )

        # South
        mask = zz < 0.
        xxS = xxR[mask]
        yyS = yyR[mask]
        vvS = val1[mask]
        lonS = full_lon[mask]
        latS = full_lat[mask]
        z3S = griddata((xxS, yyS), vvS, (xx3, yy3), method='linear')
        lon3S = griddata((xxS, yyS), lonS, (xx3, yy3), method='linear')
        lat3S = griddata((xxS, yyS), latS, (xx3, yy3), method='linear')
        figS = go.Figure(go.Contour(x=x3, y=y3, z=z3S, showscale=True, name='2DpolarS'))
        if co == plotCoord:
            figS.update_traces(zmin=cmin, zmax=cmax,
                colorbar=dict(title=varn2, tickformat=".3g"),
                customdata=np.dstack((z3S, lon3S, lat3S)),
                hovertemplate="<b>" + model + " <> SOUTH</b><br>" +
                    "Model Lon: %{customdata[1]:.2f} " + co + "<br>" +
                    "Model Lat: %{customdata[2]:.2f} " + co + "<br>" +
                    "Alt: " + "{:.2f}".format(altkm) + " km<br>" +
                    varn2 + ": %{customdata[0]:.4g} <br>" + 
                    "<extra></extra>"
                              )
        else:
            lonSplot = lon4[mask]
            latSplot = lat4[mask]
            plotlon3S = griddata((xxS, yyS), lonSplot, (xx3, yy3), method='linear')
            plotlat3S = griddata((xxS, yyS), latSplot, (xx3, yy3), method='linear')
            figS.update_traces(zmin=cmin, zmax=cmax,
                colorbar=dict(title=varn2, tickformat=".3g"),
                customdata=np.dstack((z3S, lon3S, lat3S, plotlon3N, plotlat3N)),
                hovertemplate="<b>" + model + " <> SOUTH</b><br>" +
                    "Model Lon: %{customdata[1]:.2f} " + co + "<br>" +
                    "Model Lat: %{customdata[2]:.2f} " + co + "<br>" +
                    "Alt: " + "{:.2f}".format(altkm) + " km<br>" +
                    "Plot Lon: %{customdata[3]:.2f} " + plotCoord + "<br>" +
                    "Plot Lat: %{customdata[4]:.2f} " + plotCoord + "<br>" +
                    varn2 + ": %{customdata[0]:.4g} <br>" + 
                    "<extra></extra>"
                              )

        # Shoreline
        if showshore:
            # Get points, rotate
            posN = shoreline(rscale=rscale, coord=plotCoord, utcts=plotts)
            posN[:,0], posN[:,1] = rotatePoints(posN[:,0], posN[:,1], degRotate)
            # Mask to trim to plot size
            mask = (posN[:,0] > max(x3)) | (posN[:,0] < min(x3)) | (posN[:,1] > max(y3)) | (posN[:,1] < min(y3))
            posN[mask] = np.nan
            # Duplicate for N and S copies
            posS = posN.copy()
            # Mask to blank out unused hemisphere (the ~ reverses logical)
            mask = (posN[:,2] < 0.)
            posN[mask] = np.nan
            posS[~mask] = np.nan
            figN.add_scattergl(mode='lines', x=posN[:,0], y=posN[:,1], 
                              hoverinfo='skip', showlegend=False, name='N Shoreline',
                              line=dict(width=1, color='#EEEEEE'))
            figS.add_scattergl(mode='lines', x=posS[:,0], y=posS[:,1], 
                              hoverinfo='skip', showlegend=False, name='S Shoreline',
                              line=dict(width=1, color='#EEEEEE'))

        fig = make_subplots(rows=1, cols=2)
        for i in range(len(figN.data)):
            fig.add_trace(figN.data[i], row=1, col=1)
        for i in range(len(figS.data)):
            fig.add_trace(figS.data[i], row=1, col=2)

        fig.update_layout(
            xaxis  = dict(scaleanchor="y", scaleratio=1, range=[-.6, .6], visible=False,
                domain = [0.0, 0.47]),
            xaxis2 = dict(scaleanchor="y", scaleratio=1, range=[-.6, .6], visible=False,
                domain = [0.53, 1.0], autorange="reversed"),
            yaxis  = dict(scaleanchor="x", scaleratio=1, range=[-.6, .6], visible=False),
            yaxis2 = dict(scaleanchor="x", scaleratio=1, range=[-.6, .6], visible=False),
            title_text="2D Polar Plots, North and South",
            width=800, height=400,
            margin=dict(t=50, b=60, l=25, r=120),
        )
        fig.add_annotation(xref='paper', yref='paper', x=.485, y=.50, textangle=90,
                           yanchor='middle', showarrow=False, text='Dawn')
        fig.add_annotation(xref='paper', yref='paper', x=.515, y=.50, textangle=-90,
                           yanchor='middle', showarrow=False, text='Dawn')
        fig.add_annotation(xref='paper', yref='paper', x=1.03, y=.50, textangle=90,
                           yanchor='middle', showarrow=False, text='Dusk')
        fig.add_annotation(xref='paper', yref='paper', x=-.03, y=.50, textangle=-90,
                           yanchor='middle', showarrow=False, text='Dusk')
        fig.add_annotation(xref='paper', yref='paper', x=.235, y=1.06,
                           xanchor='center', showarrow=False, text='Noon')
        fig.add_annotation(xref='paper', yref='paper', x=.765, y=1.06,
                           xanchor='center', showarrow=False, text='Noon')
        fig.add_annotation(xref='paper', yref='paper', x=.235, y=-.06,
                           xanchor='center', showarrow=False, text='Midnight')
        fig.add_annotation(xref='paper', yref='paper', x=.765, y=-.06,
                           xanchor='center', showarrow=False, text='Midnight')
        fig.add_annotation(xref='paper', yref='paper', x=.4, y=1.06,
                           showarrow=False, text='<b>NORTH</b>')
        fig.add_annotation(xref='paper', yref='paper', x=.6, y=1.06,
                           showarrow=False, text='<b>SOUTH</b>')

        for i in range(2):
            fig.add_shape(row=1, col=i+1, type="circle", line=dict(color="black", width=1),
                xref="x", yref="y", x0=-c80, y0=-c80, x1=c80, y1=c80 )
            fig.add_shape(row=1, col=i+1, type="circle", line=dict(color="black", width=1),
                xref="x", yref="y", x0=-c70, y0=-c70, x1=c70, y1=c70 )
            fig.add_shape(row=1, col=i+1, type="circle", line=dict(color="black", width=1),
                xref="x", yref="y", x0=-c60, y0=-c60, x1=c60, y1=c60 )
            fig.add_shape(row=1, col=i+1, type="circle", line=dict(color="black", width=1),
                xref="x", yref="y", x0=-c50, y0=-c50, x1=c50, y1=c50 )
            fig.add_shape(row=1, col=i+1, type="line", line=dict(color="black", width=1),
                xref="x", yref="y", x0=-c50, y0=0., x1=c50, y1=.0 )
            fig.add_shape(row=1, col=i+1, type="line", line=dict(color="black", width=1),
                xref="x", yref="y", x0=0., y0=-c50, x1=0., y1=c50 )
            fig.add_shape(row=1, col=i+1, type="line", line=dict(color="black", width=1),
                xref="x", yref="y", x0=-c50d, y0=-c50d, x1=c50d, y1=c50d )
            fig.add_shape(row=1, col=i+1, type="line", line=dict(color="black", width=1),
                xref="x", yref="y", x0=c50d, y0=-c50d, x1=-c50d, y1=c50d )

        # Set colorscale
        fig = figMods(fig, colorscale=colorscale, ncont=200, resetAnnotations=False) 

        return fig
    # /// End of 2Dpolar

    # Generate initial figure to build upon
    if plotType == 'surface':
        def plot3d_var(X=x, Y=y, Z=z):
            return val2

        kobject = Kamodo(plot_var=plot3d_var)
        fig = kobject.plot(plot_var=dict())

        fig.update_traces(
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
            varn2 + ": %{customdata[0]:.4g} <br>" +
            "<extra></extra>"
        )

    elif plotType == 'mesh3d':
        # Build connectivity grid cell by cell looping over positions
        iv = []
        jv = [] 
        kv = [] 
        for iy in range(len(lat)-1):
            for ix in range(len(lon)-1):
                # For each cell, create two triangular connectivity entries
                iv.append(ix+iy*len(lon))
                jv.append(ix+1+iy*len(lon))
                kv.append(ix+(iy+1)*len(lon))

                iv.append(ix+(iy+1)*len(lon))
                jv.append(ix+1+(iy+1)*len(lon))
                kv.append(ix+1+iy*len(lon))

        # Build resulting plot
        fig = go.Figure(data=[
            go.Mesh3d(
                x=xx, y=yy, z=zz, i=iv, j=jv, k=kv,
                colorbar_title=varn2, 
                intensity=val1, intensitymode='vertex',
                name='cont', showscale=True
            )
        ])
        fig.update_traces(
            flatshading=True,
            hovertemplate="X: %{x:.4g}<br>" + "Y: %{y:.4g}<br>" +
            "Z: %{z:.4g}<br>" + varn2 + ": %{intensity:.4g}<br>" +
            "<extra></extra>"
        )

    else:
        print('ERROR, unknown plotType: ',plotType)
        return

    # Set colorscale
    fig = figMods(fig, colorscale=colorscale, ncont=200)

    # Set plot options
    fig.update_traces(
        cmin=cmin, cmax=cmax,
        colorbar=dict(title=varn2, tickformat=".3g"),
        opacity=opacity,
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
    if showPole:
        xt = [0., 0.]
        yt = [0., 0.]
        zt = [-1.2, 1.2]
        fig.add_scatter3d(mode='lines', x=xt, y=yt, z=zt,
                          line=dict(width=4, color='black'),
                          showlegend=False, name=plotCoord+' Pole',
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
                              showlegend=False, name='GEO Pole',
                              hovertemplate='GEO Pole<extra></extra>')

        # Zero latitude (in plotCoord)
        xt = rscale*(np.cos(ilon*np.pi/180.))
        yt = rscale*(np.sin(ilon*np.pi/180.))
        zt = np.zeros(xt.shape)
        fig.add_scatter3d(mode='lines', x=xt, y=yt, z=zt,
                          line=dict(width=2, color='black'),
                          showlegend=False, name=plotCoord+' Lat=0',
                          hovertemplate=plotCoord+' Latitude=0<extra></extra>')

    # Create blank earth surface
    if showEarth:
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
                        showlegend=False, showscale=False, hoverinfo='skip',
                        name='Earth Sphere')

    # Shoreline (land/water boundaries)
    if showshore:
        pos = shoreline(rscale=1.001, coord=plotCoord, utcts=plotts)
        fig.add_scatter3d(mode='lines', x=pos[:, 0], y=pos[:, 1], z=pos[:, 2],
                          line=dict(width=shorewidth, color=shorecolor),
                          showlegend=False, hoverinfo='skip', name='Shoreline')

    return fig


### ====================================================================================== ###
def LLChangeCoord(figIn, model, altkm, plotts, inCoord, plotCoord):
    """
    Takes a gridified 2D lon/lat figure and changes the coordinate system.
    """

    import numpy as np
    from kamodo_ccmc.flythrough.utils import ConvertCoord
    from kamodo_ccmc.tools.shoreline import shoreline
    import plotly.graph_objs as go
    from scipy.interpolate import griddata

    if inCoord == plotCoord:
        print('WARNING, no change in coordinates, returning existing plot unchanged.')
        return figIn

    if figIn.data[0].type != 'contour':
        print('ERROR, plotly plot type is not contour. Returning ...')
        return

    # Pull out lon, lat, values from passed in figure
    lon = figIn.data[0]['x']
    lat = figIn.data[0]['y']
    val = figIn.data[0]['z']
    val1d = np.reshape(val, -1)

    # Prepare variables for use later
    rscale = (altkm + 6.3781E3)/6.3781E3
    lon2d, lat2d = np.meshgrid(np.array(lon), np.array(lat))
    lon1d = np.reshape(lon2d, -1)
    lat1d = np.reshape(lat2d, -1)
    alt1d = np.full(lon1d.shape, rscale)
    t = np.full(lon1d.shape, plotts)

    # Convert from plot coordinates to model coordinates
    newlon1d, newlat1d, newalt1d, units = ConvertCoord(t, lon1d, lat1d, alt1d,
        plotCoord, 'sph', inCoord, 'sph')
    newlon2d = np.reshape(newlon1d, val.shape)
    newlat2d = np.reshape(newlat1d, val.shape)

    # Interpolate on converted coordinates
    newval = griddata((lon1d, lat1d), val1d, (newlon2d, newlat2d), method='linear')
    newval2d = np.reshape(newval, val.shape)

    # Make new figure, updating values
    figOut = go.Figure(data=figIn.data[0],layout=figIn.layout)
    figOut.data[0]['z'] = newval2d

    return figOut


### ====================================================================================== ###
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


### ====================================================================================== ###
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


### ====================================================================================== ###
def swmfgm3Darb(ko, var, time=0., pos=[0, 0, 0], normal=[0, 1, 0],
                title='', lowerlabel='', showgrid=False, showibs=False,
                showE=False, log10=False, vectorMag=False):
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
    if vectorMag:
        interpX = getattr(ko, var+'_x')
        interpY = getattr(ko, var+'_y')
        interpZ = getattr(ko, var+'_z')
        varu = ko.variables[var+'_x']['units']
        varlabel = var+" ["+varu+"]"
        var = var+'_x'
    else:
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
    if vectorMag:
        valueX = interpX(grid)
        valueY = interpY(grid)
        valueZ = interpZ(grid)
        value = np.sqrt(valueX**2 + valueY**2 + valueZ**2)
    else:
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

    # Create Earth sphere
    if showE:
        elon = np.linspace(-180, 180, 181)
        elat = np.linspace(-90, 90, 91)
        elon_mg, elat_mg = np.meshgrid(np.array(elon), np.array(elat))
        ex = -1.*(np.cos(elat_mg*np.pi/180.)*np.cos(elon_mg*np.pi/180.))
        ey = -1.*(np.cos(elat_mg*np.pi/180.)*np.sin(elon_mg*np.pi/180.))
        ez = +1.*(np.sin(elat_mg*np.pi/180.))
        colorse = np.zeros(shape=ex.shape)
        colorse[ex<0.]=1  # Option to make day side a lighter color
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

### ====================================================================================== ###
def SatPosFig(satid, plotDT, coord='GSM', padHR=6, nPts=200,
              color='#d6d622', showE=True):
    '''
    Function to create a simple plotly figure from a satellite ID via SSCWeb HAPI
    
    Arguments:
      satid   ID of the satellite from
              https://hapi-server.org/servers/SSCWeb/hapi/catalog
      plotDT  datetime of desired plot time
      coord   coordinate system of positions (TOD, J2K, GEO, GM, GSE, GSM, SM)
      padHR   number of hours to show positions before plot time
      nPts    number of points to spread over the time period before and after
      color   string color of satellite positions
      showE   logical for showing a sphere to represent Earth
    '''
    import numpy as np
    import plotly.graph_objects as go
    from kamodo_ccmc.readers.hapi import HAPI
    from kamodo_ccmc.readers.hapi import hapi_get_dataset_title
    from datetime import datetime, timedelta

    # Set values to get data from the HAPI server
    server = 'https://hapi-server.org/servers/SSCWeb/hapi'
    dataset = satid
    parameters = 'X_'+coord+',Y_'+coord+',Z_'+coord
    start = (plotDT + timedelta(hours=-(1+padHR))).strftime("%Y-%m-%dT%H:%M:%S.%fZ")
    stop  = (plotDT + timedelta(hours=+(1+padHR))).strftime("%Y-%m-%dT%H:%M:%S.%fZ")
    ko_sat = HAPI(server, dataset, parameters, start, stop, register_components=True)
    satname = hapi_get_dataset_title(server, satid)
    tmp = satname.split("(")
    satname = tmp[0].strip()
    sat_vars = [*ko_sat.variables]

    # array of times, +/- padHR hours with 2*nPts+1 total points
    nPts2 = 1 + 2*nPts
    delt = np.linspace(-3600.*padHR, 3600.*padHR, nPts2)

    # set array of marker size
    deg = np.linspace(-np.pi/2., np.pi/2., nPts2)
    msize = 2.+4*(np.cos(deg)**2)
    msize[nPts+1] = 15.

    # timestamps to interpolate positions
    tss=(plotDT.timestamp() + delt)
    xx = ko_sat.variables[sat_vars[0]]['interpolator'](tss)
    yy = ko_sat.variables[sat_vars[1]]['interpolator'](tss)
    zz = ko_sat.variables[sat_vars[2]]['interpolator'](tss)

    # datetime strings for hover labels
    timestrings = [datetime.utcfromtimestamp(t).strftime("%Y-%m-%d %H:%M:%S") for t in tss]

    fig = go.Figure()
    fig.add_trace(go.Scatter3d(x=xx, y=yy, z=zz, mode='lines+markers',
        line=dict(color=color, width=2),
        marker=dict(color=color, size=msize, line=dict(color=color, width=1)),
        name=satname, showlegend=True))
    fig.update_traces(
        customdata = timestrings,
        hovertemplate=coord + " Coordinates:<br>" +
        "X [R_E]: %{x:.5g}<br>" +
        "Y [R_E]: %{y:.5g}<br>" +
        "Z [R_E]: %{z:.5g}<br>" + "%{customdata}" +
        "<extra>" + satname + "</extra>",
    )

    # Create Earth sphere
    if showE:
        elon = np.linspace(-180, 180, 181)
        elat = np.linspace(-90, 90, 91)
        elon_mg, elat_mg = np.meshgrid(np.array(elon), np.array(elat))
        ex = -1.*(np.cos(elat_mg*np.pi/180.)*np.cos(elon_mg*np.pi/180.))
        ey = -1.*(np.cos(elat_mg*np.pi/180.)*np.sin(elon_mg*np.pi/180.))
        ez = +1.*(np.sin(elat_mg*np.pi/180.))
        colorse = np.zeros(shape=ex.shape)
        colorse[ex<0.]=1  # Option to make day side a lighter color
        colorscalee = ['rgb(199,199,199)', 'rgb(0,0,0)']
        fig.add_surface(x=ex, y=ey, z=ez, surfacecolor=colorse,
                        cmin=0, cmax=1, colorscale=colorscalee,
                        showlegend=False, showscale=False, hoverinfo='skip')

    fig.update_layout(scene_aspectmode='data',
        title='Plot of '+satname+' at '+timestrings[0]+'<br>(with '+str(padHR)+
              ' hours before and after in '+coord+' coordinates)')
    
    return fig

### ====================================================================================== ###
def SatOrbitPlane(satid, refDT, coord='GSM', color='#d6d622'):
    '''
    A function to look at a satellite at a given time and it's current orbit.
      The orbit is broken when it crosses midnight

    Input arguements:
      satid   ID for the satellite (using SSCWeb HAPI naming)
      refDT   a datetime object representing the satellite current position
      coord   the coordinate system of the satellite
      color   the color for the satellite trajectory in the plot

    Output:
      Odict   a dictionary with orbit plane centroid/normal and other metadata
      Oko     a Kamodo object for the trajectory
      Ofig    a plotly figure of the orbit
                -trace 0 is orbit path, trace 1 is satellite position, trace 2 is Earth
    '''
    import numpy as np
    from kamodo_ccmc.readers.hapi import HAPI
    from kamodo_ccmc.readers.hapi import hapi_get_dataset_title
    from datetime import datetime, timedelta
    import plotly.graph_objects as go

    # NOTE, max +/- 48 hours on orbit for now
    # NOTE, orbit cutoff is negative X and change of sign in Y (crossing midnight)

    parameters = 'X_'+coord+',Y_'+coord+',Z_'+coord
    server = 'https://hapi-server.org/servers/SSCWeb/hapi'

    # Start looking back in time from reference point
    start = (refDT + timedelta(hours=-48)).strftime("%Y-%m-%dT%H:%M:%S.%fZ")
    stop  = refDT.strftime("%Y-%m-%dT%H:%M:%S.%fZ")
    Oko = HAPI(server, satid, parameters, start, stop, register_components=True)
    vars = [*Oko.variables]
    pts = len(Oko[vars[0]].data)
    found1DT = Oko.dtarray[0]
    found1TS = Oko.tsarray[0]
    for i in range(pts-1):
        j = (pts-1) - i
        if Oko[vars[0]].data[j] < 0. and (Oko[vars[1]].data[j] * Oko[vars[1]].data[j-1]) <= 0.:
            found1DT = Oko.dtarray[j]
            found1TS = Oko.tsarray[j]
            break

    # Then forward in time from reference point
    start  = refDT.strftime("%Y-%m-%dT%H:%M:%S.%fZ")
    stop  = (refDT + timedelta(hours=+48)).strftime("%Y-%m-%dT%H:%M:%S.%fZ")
    Oko = HAPI(server, satid, parameters, start, stop, register_components=True)
    pts = len(Oko[vars[0]].data)
    found2DT = Oko.dtarray[-1]
    found2TS = Oko.tsarray[-1]
    for i in range(pts-1):
        if Oko[vars[0]].data[i] < 0. and (Oko[vars[1]].data[i] * Oko[vars[1]].data[i+1]) <= 0.:
            found2DT = Oko.dtarray[i]
            found2TS = Oko.tsarray[i]
            break

    # Get actual range to check
    start = (found1DT + timedelta(seconds=-1)).strftime("%Y-%m-%dT%H:%M:%S.%fZ")
    stop  = (found2DT + timedelta(seconds=+1)).strftime("%Y-%m-%dT%H:%M:%S.%fZ")
    Oko = HAPI(server, satid, parameters, start, stop, register_components=True)
    aveX = np.mean(Oko[vars[0]].data)
    aveY = np.mean(Oko[vars[1]].data)
    aveZ = np.mean(Oko[vars[2]].data)
    pts = len(Oko[vars[0]].data)
    vecNormal = [0, 0, 0]
    for i in range(pts-1):
        vec1 = [Oko[vars[0]].data[i]-Oko[vars[0]].data[i+1],
                Oko[vars[1]].data[i]-Oko[vars[1]].data[i+1],
                Oko[vars[2]].data[i]-Oko[vars[2]].data[i+1]]
        vec2 = [Oko[vars[0]].data[i]-aveX,
                Oko[vars[1]].data[i]-aveY,
                Oko[vars[2]].data[i]-aveZ]
        vec3 = np.cross(vec1, vec2)
        vec3 = vec3/np.linalg.norm(vec3)
        vecNormal += vec3
    vecNormal = vecNormal/np.linalg.norm(vecNormal)

    # Other values
    refTS = refDT.timestamp()
    refSTR = refDT.strftime("%Y-%m-%d %H:%M:%S UT")
    satname = hapi_get_dataset_title(server, satid)
    tmp = satname.split("(") 
    satname = tmp[0].strip() 

    Odict = dict(satid=satid, satname=satname, satvars=vars, coord=coord,
                 refDT=refDT, refTS=refTS, refSTR=refSTR,
                 startDT=found1DT, startTS=found1TS,
                 stopDT=found2DT, stopTS=found2TS,
                 orbitPoint=[aveX, aveY, aveZ], orbitNormal=vecNormal)

    timestrings = [datetime.utcfromtimestamp(t).strftime("%Y-%m-%d %H:%M:%S") for t in Oko.tsarray]
    Ofig = go.Figure()
    Ofig.add_trace(go.Scatter3d(x=Oko.variables[vars[0]]['data'], 
                                y=Oko.variables[vars[1]]['data'], 
                                z=Oko.variables[vars[2]]['data'], mode='lines',
                                line=dict(color=color, width=3),
                                name=satname+' orbit', showlegend=True))
    Ofig.update_traces(customdata = timestrings,
        hovertemplate=satname + " orbit<br>" + coord + " Coordinates:<br>" +
        "X [R_E]: %{x:.5g}<br>" +
        "Y [R_E]: %{y:.5g}<br>" +
        "Z [R_E]: %{z:.5g}<br>" + "%{customdata}" +
        "<extra></extra>",
        selector=dict(name=satname+' orbit')
    )

    sx = Oko.variables[vars[0]]['interpolator'](refTS)
    sy = Oko.variables[vars[1]]['interpolator'](refTS)
    sz = Oko.variables[vars[2]]['interpolator'](refTS)
    Ofig.add_trace(go.Scatter3d(x=[sx], y=[sy], z=[sz], mode='markers',
                                marker=dict(color='black', size=5, line=dict(color='black', width=1)),
                                name=satname+' position', showlegend=True))
    Ofig.update_traces(
        hovertemplate=satname + " position<br>" + coord + " Coordinates:<br>" +
        "X [R_E]: %{x:.5g}<br>" +
        "Y [R_E]: %{y:.5g}<br>" +
        "Z [R_E]: %{z:.5g}<br>" + refSTR +
        "<extra></extra>",
        selector=dict(name=satname+' position')
    )

    elon = np.linspace(-180, 180, 181)
    elat = np.linspace(-90, 90, 91)
    elon_mg, elat_mg = np.meshgrid(np.array(elon), np.array(elat))
    ex = -1.*(np.cos(elat_mg*np.pi/180.)*np.cos(elon_mg*np.pi/180.))
    ey = -1.*(np.cos(elat_mg*np.pi/180.)*np.sin(elon_mg*np.pi/180.))
    ez = +1.*(np.sin(elat_mg*np.pi/180.))
    colorse = np.zeros(shape=ex.shape)
    colorse[ex<0.]=1  # Option to make day side a lighter color
    colorscalee = ['rgb(199,199,199)', 'rgb(0,0,0)']
    Ofig.add_surface(x=ex, y=ey, z=ez, surfacecolor=colorse,
                     cmin=0, cmax=1, colorscale=colorscalee,
                     name='Earth', showlegend=False, showscale=False, hoverinfo='skip')

    Ofig.update_layout(scene_aspectmode='data',
        title='Plot of '+satname+' orbit and position at '+refSTR+' in '+coord+' coordinates')

    return Odict, Oko, Ofig

### ====================================================================================== ###
def fixFigOrigin(fig, setX='', setY='', setZ=''):
    '''
    Simple function that takes a passed in 3D plotly figure object
    and computes a new camera center near the origin (0,0,0) by examining
    the x, y, and z values of the data blocks, returning an updated figure.
    '''
    minX, maxX, minY, maxY, minZ, maxZ = 0., 0., 0., 0., 0., 0.
    for i in range(len(fig.data)):
        if len(fig.data[i]['x'].shape) == 1:
            minX = min(minX, min(v for v in fig.data[i]['x'] if v is not None))
            minY = min(minY, min(v for v in fig.data[i]['y'] if v is not None))
            minZ = min(minZ, min(v for v in fig.data[i]['z'] if v is not None))
            maxX = max(maxX, max(v for v in fig.data[i]['x'] if v is not None))
            maxY = max(maxY, max(v for v in fig.data[i]['y'] if v is not None))
            maxZ = max(maxZ, max(v for v in fig.data[i]['z'] if v is not None))
    aveX = (minX + maxX)/2.
    aveY = (minY + maxY)/2.
    aveZ = (minZ + maxZ)/2.
    newX = 0.7*(0. - aveX)/(maxX - aveX)
    newY = 0.7*(0. - aveY)/(maxY - aveY)
    newZ = 0.7*(0. - aveZ)/(maxZ - aveZ)
    if setX != '': newX = float(setX)
    if setY != '': newY = float(setY)
    if setZ != '': newZ = float(setZ)
    camera = dict(center=dict(x=newX, y=newY, z=newZ))
    fig.update_layout(scene_camera=camera)

    return fig

### ====================================================================================== ###
def B3Dfig(fullfile, showE = True):
    '''
    Function to take a specially generated file of last closed fieldline traces
      and generate a plotly figure of it. This can be displayed as is or added
      to other figures.
    
    Arguments:
      fullfile:  The full file path to the file with fieldline extractions
      showE:     A logical for showing an Earth sphere in the figure
    '''
    import numpy as np
    import pandas as pd
    import plotly.graph_objects as go

    file = fullfile.split('/')[-1]
    try:
        tmpD=pd.read_csv(fullfile, delimiter=r"\s+", header=None,
			 comment='#', skiprows=0).to_numpy()
    except FileNotFoundError:
        print('ERROR, no file found.')
        return None

    # The precomputed fieldlines file has 29 lines with 41 points each, FIXED size
    tmpD3=tmpD.reshape(29,41,3)
    # Placing a None between lines keeps them separated so we can use
    #   one trace to plot them all.
    vnone = np.full((29,1,3), None)
    newD = np.hstack((tmpD3,vnone)).reshape(29*42,3)
    # Make the figure
    fig = go.Figure()
    fig.add_trace(go.Scatter3d(x=newD[:,0], y=newD[:,1], z=newD[:,2], mode='lines',
        line=dict(color="#a1a1a1", width=2), hoverinfo='skip', name='Last Closed B'))

    # Create Earth sphere
    if showE:
        elon = np.linspace(-180, 180, 181)
        elat = np.linspace(-90, 90, 91)
        elon_mg, elat_mg = np.meshgrid(np.array(elon), np.array(elat))
        ex = -1.*(np.cos(elat_mg*np.pi/180.)*np.cos(elon_mg*np.pi/180.))
        ey = -1.*(np.cos(elat_mg*np.pi/180.)*np.sin(elon_mg*np.pi/180.))
        ez = +1.*(np.sin(elat_mg*np.pi/180.))
        colorse = np.zeros(shape=ex.shape)
        colorse[ex<0.]=1  # Option to make day side a lighter color
        colorscalee = ['rgb(199,199,199)', 'rgb(0,0,0)']
        fig.add_surface(x=ex, y=ey, z=ez, surfacecolor=colorse,
                        cmin=0, cmax=1, colorscale=colorscalee,
                        showlegend=False, showscale=False, hoverinfo='skip')

    # Add hidden plot ranges (may need to adjust later)
    xx = np.full((3), None)
    yy = np.full((3), None)
    zz = np.full((3), None)
    xx[0], xx[2] = -75.,  15.
    yy[0], yy[2] =  25., -25.
    zz[0], zz[2] = -20.,  20.
    fig.add_trace(go.Scatter3d(x=xx, y=yy, z=zz, mode='markers',
        marker=dict(color="black", size=1), showlegend=False, hoverinfo='skip' ))

    # Set 3D view
    camera = dict(
        up=dict(x=0, y=0, z=1),
        eye=dict(x=1.1, y=0.85, z=0.25)
    )
    fig.update_layout(title=file, scene_camera=camera)
    fig.update_layout(scene_aspectmode='data')

    return fig

### ====================================================================================== ###
def gm3DSlicePlus(ko, var, timeHrs=0., pos=[0, 0, 0], normal=[0, 0, 1],
                  lowerlabel='', colorscale='RdBu', addTraceTime=False,
                  showgrid=False, gdeg=2., showE=False, log10=False, 
                  showMP=True, showBS=True, wireframe=False, crange='', 
                  xrange=[-999., 999.], yrange=[-999., 999.], zrange=[-999., 999.]):
    '''
    Function to create a slice for any point and normal and interpolate
      from Kamodo object onto that grid. Returns a full 3D plotly figure.
      Options to show grid (dots for vertices), magnetopause boundary, and
      bow shock boundary. Many options to customize the figure are available.

    ko:           Kamodo object
    var:          string variable name
    timeHrs:      floating time in hours from start of first day of data
    pos, normal:  position and normal vector for slice plane
    lowerlabel:   string text to label plot lower left
    colorscale:   string name of Python colorscale, ie. Viridis, Cividis, RdBu
    addTraceTime: logical to add time string to trace label
    showgrid:     logical to show dots at grid locations
    gdeg:         slice grid degree resolution, default is 2.
    showE:        logical to show a sphere at R=1
    log10:        logical to take log10 of display value
    showMP:       logical to show magnetopause boundary (requires 'status' variable)
    showBS:       logical to show bow shock (requires 'v_x' variable)
    wireframe:    logical to show MP and BS as wireframe
    crange:       2 value array for min/max contour range
    xrange:       2 value array for min/max extent of X values in slice
    yrange:       2 value array for min/max extent of Y values in slice
    zrange:       2 value array for min/max extent of Z values in slice
      -note: xyz range values default as large values that will trim to
             actual model ranges. User can set small ranges.
    '''
    import numpy as np
    import plotly.graph_objs as go
    import kamodo_ccmc.flythrough.model_wrapper as MW
    from datetime import datetime

    # Error checking  =============================================== Validate
    if xrange[0] > xrange[1]:
        print('Error in xrange: ',xrange)
        return
    if yrange[0] > yrange[1]:
        print('Error in yrange: ',yrange)
        return
    if zrange[0] > zrange[1]:
        print('Error in zrange: ',zrange)
        return
    if "_ijk" in var:
        #print('Removing _ijk from variable name')
        var = var.replace('_ijk','')
    vMag = False
    if var not in ko:
        if var+'_x' in ko:
            vMag = True
            interpX = getattr(ko, var+'_x')
            interpY = getattr(ko, var+'_y')
            interpZ = getattr(ko, var+'_z')
            vunits = ko.variables[var+'_x']['units']
            var2 = var+'_x_ijk'
            cr = MW.Coord_Range(ko, [var2], return_dict=True, print_output=False)
            xunits = cr[var2]['X'][2]
            coord = MW.Variable_Search('', model=ko.modelname, return_dict=True)[var+'_x'][2]
        else:
            print('Error, variable not in Kamodo object')
            return
    else:
        interp = getattr(ko, var)
        vunits = ko.variables[var]['units']
        var2 = var+'_ijk'
        cr = MW.Coord_Range(ko, [var2], return_dict=True, print_output=False)
        xunits = cr[var2]['X'][2]
        coord = MW.Variable_Search('', model=ko.modelname, return_dict=True)[var][2]
    varlabel = var+" ["+vunits+"]"

    # Make sure range is not larger than actual data range
    x1, x2 = cr[var2]['X'][0], cr[var2]['X'][1]
    y1, y2 = cr[var2]['Y'][0], cr[var2]['Y'][1]
    z1, z2 = cr[var2]['Z'][0], cr[var2]['Z'][1]
    xrange[0] = max(x1, xrange[0])
    xrange[1] = min(x2, xrange[1])
    yrange[0] = max(y1, yrange[0])
    yrange[1] = min(y2, yrange[1])
    zrange[0] = max(z1, zrange[0])
    zrange[1] = min(z2, zrange[1])

    plotTS = ko.filedate.timestamp() + timeHrs*3600.
    plotDT = datetime.utcfromtimestamp(plotTS)
    plotDateStr = plotDT.strftime("%Y-%m-%d %H:%M:%S UT")

    sliceDir = ''
    if normal[0] == 1: sliceDir = 'X '
    if normal[1] == 1: sliceDir = 'Y '
    if normal[2] == 1: sliceDir = 'Z ' 

    # Compute base 1D grid values  ====================================== Grid
    # Compute max radius of range box
    xm = max(abs(xrange[0]),abs(xrange[1]))
    ym = max(abs(yrange[0]),abs(yrange[1]))
    zm = max(abs(zrange[0]),abs(zrange[1]))
    rmax = np.sqrt(xm*xm + ym*ym + zm*zm)

    # Compute values from pos, normal values
    uvec = normal/np.linalg.norm(normal)  # unit normal vector
    odist = np.dot(uvec, pos)  # closest distance to slice from origin
    opos = odist*uvec  # vector from origin to closest point

    # Compute base grid
    ndeg = 1+int(360./gdeg)
    dg = np.linspace(-180, 180, ndeg)  # degree grid
    deldg = dg[1]-dg[0]
    rv = 0.  # radius value
    if abs(odist) < 2.75:
        rv = 2.75*np.sin(((2.75-odist)/2.75)*np.pi/2.)
    rg = []  # radius grid
    rg.append(rv)
    rtrans = 2.5  # radius transition point from fixed dr to variable
    rvfactor = 0.11/(2.5/deldg)
    for _ in range(999):
        if rv < rtrans:
            rv += rvfactor
        else:
            rv += rvfactor*rv/rtrans
        rg.append(rv)
        if rv > rmax:
            break
    rg = np.array(rg)
    dg_mg, rg_mg = np.meshgrid(dg, rg)
    gx = rg_mg*np.cos(dg_mg*np.pi/180.)
    gy = rg_mg*np.sin(dg_mg*np.pi/180.)
    gx_1d = gx.reshape(-1)
    gy_1d = gy.reshape(-1)
    gz_1d = np.zeros([len(gx_1d)])
    time_1d = np.full((len(gx_1d)), timeHrs)
    grid0 = np.stack((gx_1d, gy_1d, gz_1d), axis=-1)  # nx3 position grid

    # Transform base grid to pos/normal and trim if needed
    # Rotate
    if abs(uvec[2]) < 1.:  # No rotation for +/- Z normal
        new_xaxis = np.cross([0, 0, 1], uvec)
        new_xaxis = new_xaxis/np.linalg.norm(new_xaxis)  # unit vector
        new_yaxis = np.cross(new_xaxis, uvec)
        new_yaxis = new_yaxis/np.linalg.norm(new_yaxis)  # unit vector
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
    x1, x2 = xrange[0], xrange[1]
    y1, y2 = yrange[0], yrange[1]
    z1, z2 = zrange[0], zrange[1]
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

    # Interpolate  ======================================================= Interp
    # Create 4D (nx4) grid and interpolate plot values
    grid = np.ndarray(shape=(len(gx_1d), 4), dtype=np.float32)
    grid[:, 0] = time_1d
    grid[:, 1] = gx2.reshape(-1)
    grid[:, 2] = gy2.reshape(-1)
    grid[:, 3] = gz2.reshape(-1)
    if vMag:
        valueX = interpX(grid)
        valueY = interpY(grid)
        valueZ = interpZ(grid)
        value = np.sqrt(valueX**2 + valueY**2 + valueZ**2)
    else:
        value = interp(grid)
    if log10:
        value[value <= 0.] = np.nan
        value = np.log10(value)
        varlabel = "log10("+varlabel+")"
    vmin = np.nanmin(value)
    vmax = np.nanmax(value)

    # Build connectivity grid cell by cell looping over positions
    iv, jv, kv = [], [], []
    for iy in range(len(rg)-1):
        for ix in range(len(dg)-1):
            # For each cell, create two triangular connectivity entries
            iv.append(ix  + iy   *len(dg))
            jv.append(ix+1+ iy   *len(dg))
            kv.append(ix  +(iy+1)*len(dg))
            iv.append(ix  +(iy+1)*len(dg))
            jv.append(ix+1+(iy+1)*len(dg))
            kv.append(ix+1+ iy   *len(dg))

    # Build resulting plot
    traceTime = ''
    if addTraceTime:
        traceTime = ' at '+plotDateStr
    fig = go.Figure(data=[
        go.Mesh3d(
            x=grid[:, 1], y=grid[:, 2], z=grid[:, 3], i=iv, j=jv, k=kv,
            colorbar_title=varlabel, colorscale=colorscale,
            intensity=value, intensitymode='vertex',
            name=sliceDir+'Slice'+traceTime, showscale=True, showlegend=True
        )
    ])
    if crange != '':
        cmin = float(crange[0])
        cmax = float(crange[1])
        fig.update_traces(cmin=cmin, cmax=cmax)
    fig.update_traces(colorbar=dict(lenmode='fraction', len=0.4, y=0.2))
    xlabel = 'X [' + xunits + ']'
    ylabel = 'Y [' + xunits + ']'
    zlabel = 'Z [' + xunits + ']'
    fig.update_traces(
        flatshading=True,
        hovertemplate=xlabel + ": %{x:.4g}<br>" + 
                      ylabel + ": %{y:.4g}<br>" +
                      zlabel + ": %{z:.4g}<br>" + 
                      varlabel + ": %{intensity:.4g}<br>" +
        "<extra></extra>"
    )
    # add invisible trace to keep view consistent
    ix = [xrange[0], None, xrange[1]]
    iy = [yrange[0], None, yrange[1]]
    iz = [zrange[0], None, zrange[1]]
    sz = [0.1, 0.1, 0.1]
    fig.add_scatter3d(x=np.array(ix), y=np.array(iy), z=np.array(iz), mode='markers',
        marker=dict(size=1, color='red', opacity=0.10),
        showlegend=False, hoverinfo='skip', name='FixedRange' )
    # Add grid points
    if showgrid:
        fig.add_scatter3d(name='grid',
            x=grid[:, 1], y=grid[:, 2], z=grid[:, 3], mode='markers',
            marker=dict(size=1, color='grey'), line=dict(width=1) )

    # Create Earth sphere
    if showE:
        elon = np.linspace(-180, 180, 181)
        elat = np.linspace(-90, 90, 91)
        elon_mg, elat_mg = np.meshgrid(np.array(elon), np.array(elat))
        ex = -1.*(np.cos(elat_mg*np.pi/180.)*np.cos(elon_mg*np.pi/180.))
        ey = -1.*(np.cos(elat_mg*np.pi/180.)*np.sin(elon_mg*np.pi/180.))
        ez = +1.*(np.sin(elat_mg*np.pi/180.))
        ecolors = np.zeros(shape=ex.shape)
        ecolors[ex<0.]=1  # Option to make day side a lighter color
        ecolorscale = ['rgb(199,199,199)', 'rgb(0,0,0)']
        fig.add_surface(x=ex, y=ey, z=ez, surfacecolor=ecolors,
            cmin=0, cmax=1, colorscale=ecolorscale,
            showlegend=False, showscale=False, hoverinfo='skip')

    # Magnetopause surface
    if showMP:
        fig2,success = gmGetSurfacePlot(ko=ko,timeHrs=timeHrs,wireframe=wireframe,
            what='MP',traceTime=False)
        if success:
            fig.add_trace(fig2.data[0])

    # Bow shock surface
    if showBS:
        fig2,success = gmGetSurfacePlot(ko=ko,timeHrs=timeHrs,wireframe=wireframe,
            what='BS',traceTime=False)
        if success:
            fig.add_trace(fig2.data[0])

    # Final figure modifications
    xs2 = 18
    ys1 = -40
    ys2 = -24
    ys3 = -8
    lrText1 = coord+' Coordinates'
    lrText2 = 'Min = '+"{:.4e}".format(vmin)
    lrText3 = 'Max = '+"{:.4e}".format(vmax)
    fig.update_layout(
        scene_aspectmode='data',
        title=dict(text=plotDateStr,
                   yref="container", yanchor="top", x=0.01, y=0.97,
                   font=dict(size=16, family="sans serif", color="#000000")),
        scene=dict(
            xaxis=dict(showbackground=False, showgrid=False),
            yaxis=dict(showbackground=False, showgrid=False),
            zaxis=dict(showbackground=False, showgrid=False)),
        annotations=[
            dict(text=lowerlabel, x=0.0, y=0.0, ax=0, ay=0, xanchor="left",
                 xshift=0, yshift=-20, xref="paper", yref="paper",
                 font=dict(size=16, family="sans serif", color="#000000")),
            dict(text=lrText1, x=1.0, y=0.0, ax=0, ay=0, xanchor="left",
                xshift=xs2, yshift=ys1, xref="paper", yref="paper",
                font=dict(size=12, family="sans serif", color="#000000")),
            dict(text=lrText2, x=1.0, y=0.0, ax=0, ay=0, xanchor="left",
                xshift=xs2, yshift=ys2, xref="paper", yref="paper",
                font=dict(size=12, family="sans serif", color="#000000")),
            dict(text=lrText3, x=1.0, y=0.0, ax=0, ay=0, xanchor="left",
                xshift=xs2, yshift=ys3, xref="paper", yref="paper",
                font=dict(size=12, family="sans serif", color="#000000")),
         ],
        margin=dict(l=10, t=35),
    )

    return fig

### ====================================================================================== ###
def gmGetSurfacePlot(ko='', timeHrs='', wireframe=False, Gridsize=21, what='BS',
                     traceTime=True, sfile=''):
    '''
    Function to get bow shock location in a plotly figure.

    ko:           Kamodo object
    timeHrs:      time to interpolate values (hrs from midnight of 1st day of data)
    wireframe:     logical to set wireframe or opaque surface for returned plot
    Gridsize:     surface grid size, default is 21x21 grid of points
    what:         string with what to retrieve from: BS, MP
    traceTime:    logical to add time to trace label
    sfile:        string of relative directory path to read surface file
    '''
    import os
    import numpy as np
    import plotly.graph_objs as go
    from datetime import datetime

    surfaces = ['MP', 'BS', 'CS']
    if what in surfaces:
        color, name = '#666666', 'ERROR'
        if sfile == '':
            x,y,z = gmComputeSurface(ko,timeHrs,Gridsize=Gridsize,what=what)
            DT = datetime.utcfromtimestamp( ko.filedate.timestamp() + timeHrs*3600. )
            DateStr = DT.strftime("%Y/%m/%d %H:%M:%S UT")
            titleStr = what+' Surface for '+DateStr
        else:
            success,x,y,z = gmLoadSurface(sfile)
            if not success:
                return None,False
            bname = os.path.basename(sfile)
            bname = bname[:-5]
            p = bname.split("_")
            DateStr = p[1][0:4]+'/'+p[1][4:6]+'/'+p[1][6:8]+\
                ' '+p[2][0:2]+':'+p[2][2:4]+':'+p[2][4:6]+' UT'
            titleStr = p[0]+' Surface for '+DateStr
            what = p[0]
        if what == 'MP':
            color, name = '#0000c1', 'Magnetopause'
        elif what == 'BS':
            color, name = '#00ff00', 'Bow Shock'
        elif what == 'CS':
            color, name = '#223344', 'Tail Current Sheet'
            return None,False
        if traceTime:
            name2 = name+'<br>'+DateStr
        else:
            name2 = name
    else:
        print('ERROR, unknown surface in gmGetSurfacePlot.')
        return None,False

    Tpts = x.shape[0]
    if Tpts < 3:
        # Grid size must be at least 3x3
        return None,False
        
    if wireframe:
        # Create wire mesh view of BS
        xx, yy, zz = [], [], []
        for i in range(Tpts):
            xx = np.concatenate((xx, [None], x[:,i]))
            yy = np.concatenate((yy, [None], y[:,i]))
            zz = np.concatenate((zz, [None], z[:,i]))
        for i in range(Tpts):
            xx = np.concatenate((xx, [None], x[i,:]))
            yy = np.concatenate((yy, [None], y[i,:]))
            zz = np.concatenate((zz, [None], z[i,:]))
        fig2 = go.Figure()
        fig2.add_trace(go.Scatter3d(x=xx, y=yy, z=zz, mode='lines', name=name2,
            line=dict(color=color, width=2), hoverinfo='skip'))
    else:
        # Create transparent shell view of BS
        iv, jv, kv = [], [], []
        for iy in range(Tpts-1):
            for ix in range(Tpts-1):
                # For each cell, create two triangular connectivity entries
                iv.append(ix   +  iy   *Tpts)
                jv.append(ix+1 +  iy   *Tpts)
                kv.append(ix   + (iy+1)*Tpts)
                iv.append(ix   + (iy+1)*Tpts)
                jv.append(ix+1 + (iy+1)*Tpts)
                kv.append(ix+1 +  iy   *Tpts)
        fig2 = go.Figure(data=[go.Mesh3d(name=name2, flatshading=True, 
            x=x.reshape(-1), y=y.reshape(-1), z=z.reshape(-1), i=iv, j=jv, k=kv, 
            color=color, opacity=0.40, hoverinfo='skip', showlegend=True, showscale=False,
        )])

    # put an Earth sphere on plot
    elon = np.linspace(-180, 180, 181)
    elat = np.linspace(-90, 90, 91)
    elon_mg, elat_mg = np.meshgrid(np.array(elon), np.array(elat))
    ex = -1.*(np.cos(elat_mg*np.pi/180.)*np.cos(elon_mg*np.pi/180.))
    ey = -1.*(np.cos(elat_mg*np.pi/180.)*np.sin(elon_mg*np.pi/180.))
    ez = +1.*(np.sin(elat_mg*np.pi/180.))
    ecolors = np.zeros(shape=ex.shape)
    ecolors[ex<0.]=1  # Option to make day side a lighter color
    ecolorscale = ['rgb(199,199,199)', 'rgb(0,0,0)']
    fig2.add_surface(x=ex, y=ey, z=ez, surfacecolor=ecolors,
        cmin=0, cmax=1, colorscale=ecolorscale,
        showlegend=False, showscale=False, hoverinfo='skip')

    # add invisible trace to keep view consistent
    ix = [-15., None, 35.]
    iy = [-25., None, 25.]
    iz = [-25., None, 25.]
    sz = [0.1, 0.1, 0.1]
    fig2.add_scatter3d(x=ix, y=iy, z=iz, mode='markers',
        marker=dict(size=1, color='red', opacity=0.10),
        showlegend=False, hoverinfo='skip' )

    fig2.update_layout(
        scene_aspectmode='data',
        title=dict(text=titleStr,
                   yref="container", yanchor="top", x=0.01, y=0.97,
                   font=dict(size=16, family="sans serif", color="#000000")) )

    return fig2,True

### ====================================================================================== ###
def gmComputeSurface(ko, timeHrs, Gridsize=21, what='MP'):
    '''
    Function to compute grid points defining a surface.

    ko:           Kamodo object
    timeHrs:      time to interpolate values (hrs from midnight of 1st day of data)
    Gridsize:     surface grid size, default is 21x21 grid of points
    what:         string with what to retrieve from: BS, MP
    '''
    import numpy as np
    import kamodo_ccmc.flythrough.model_wrapper as MW

    # Error array to return if needed.
    e = np.full((1,1), None)

    if what == 'MP':
        # Make sure the 'status' variable is in the Kamodo object and create interpolator.
        if 'status' in ko:
            interpMP = getattr(ko, 'status')
        else:
            showBS = False
            print('Warning, no status variable in Kamodo object. No magnetopause computed.')
            return e,e,e

        # Get min/max of X range
        cr = MW.Coord_Range(ko, ['status'], return_dict=True, print_output=False)
        x1, x2 = cr['status']['X'][0], cr['status']['X'][1]
        y1, y2 = cr['status']['Y'][0], cr['status']['Y'][1]
        z1, z2 = cr['status']['Z'][0], cr['status']['Z'][1]

        # Look at status=0.5
        x = np.full((Gridsize,Gridsize), None)
        y = np.full((Gridsize,Gridsize), None)
        z = np.full((Gridsize,Gridsize), None)
        a1 = np.linspace(-np.pi/2., np.pi/2., Gridsize)
        a2 = np.linspace(-np.pi/(12./7.), np.pi/(12./7.), Gridsize) # +/- 7 hrs from noon
        for i, a in enumerate(a1):
            for j, b in enumerate(a2):
                xval = np.cos(b)
                yval = np.sin(b)*np.sin(a)
                zval = np.sin(b)*np.cos(a)
                dr, r = 2., 5.
                newx, newy, newz = r*xval, r*yval, r*zval
                v = interpMP([timeHrs, newx, newy, newz])
                for _ in range(300):
                    if v[0] < 0.5:
                        if dr > 0.: dr = -0.5 * dr
                    else:
                        if dr < 0.: dr = -0.5 * dr
                    r += dr
                    newx, newy, newz = r*xval, r*yval, r*zval
                    if newx < x1 or newx > x2: break
                    if newy < y1 or newy > y2: break
                    if newz < z1 or newz > z2: break
                    #\\ TEST NEW STUFF
                    if newx > 15.: break
                    newr = np.sqrt(newx*newx + newy*newy + newz*newz)
                    if newr > 22.: break
                    #//
                    v = interpMP([timeHrs, newx, newy, newz])
                    if abs(dr) < 1.e-3:
                        x[i, j] = newx
                        y[i, j] = newy
                        z[i, j] = newz
                        break
        return x,y,z
    elif what == 'BS':
        # Make sure the 'v_x' variable is in the Kamodo object and create interpolator.
        if 'v_x' in ko:
            interpBS = getattr(ko, 'v_x')
        else:
            showBS = False
            print('Warning, no v_x variable in Kamodo object. No bow shock computed.')
            return e,e,e

        # Get min/max of X range
        cr = MW.Coord_Range(ko, ['v_x'], return_dict=True, print_output=False)
        x1, x2 = cr['v_x']['X'][0], cr['v_x']['X'][1]
        y1, y2 = cr['v_x']['Y'][0], cr['v_x']['Y'][1]
        z1, z2 = cr['v_x']['Z'][0], cr['v_x']['Z'][1]

        # Shock looking at v_x where it drops to 85% of value upstream (note v_x negative in SW)
        x = np.full((Gridsize,Gridsize), None)
        y = np.full((Gridsize,Gridsize), None)
        z = np.full((Gridsize,Gridsize), None)
        a1 = np.linspace(-np.pi/2., np.pi/2., Gridsize)
        a2 = np.linspace(-np.pi/3., np.pi/3., Gridsize)
        Npts = 200
        rbs = np.linspace(x2-1., -10., Npts)
        Rfactor = 25./np.sin(a2[0]) # R=25 from X axis for surface
        ii = 0
        for i, a in enumerate(a1):
            for j, b in enumerate(a2):
                xval = x2-1.
                yval = Rfactor*np.sin(b)*np.sin(a)
                zval = Rfactor*np.sin(b)*np.cos(a)
                v = interpBS([timeHrs, xval, yval, zval])
                basevalue = v[0]
                dx = -1.
                xval += dx
                v = interpBS([timeHrs, xval, yval, zval])
                for _ in range(300):
                    if v[0] > (0.85*basevalue):
                        if dx < 0.: dx = -0.5 * dx
                    else:
                        if dx > 0.: dx = -0.5 * dx
                    xval += dx
                    if xval < x1 or xval > x2: break
                    if yval < y1 or yval > y2: break
                    if zval < z1 or zval > z2: break
                    if xval < -10.: break
                    v = interpBS([timeHrs, xval, yval, zval])
                    if abs(dx) < 1.e-3:
                        x[i, j] = xval
                        y[i, j] = yval
                        z[i, j] = zval
                        break
        return x,y,z

    # Default if it somehow falls through to the end
    return e,e,e

### ====================================================================================== ###
def gmSaveSurface(ko, timeHrs, Gridsize=21, what='MP', where='.', runname='unknown'):
    '''
    Function to save computed grid points defining a surface.

    ko:           Kamodo object
    timeHrs:      time to interpolate values (hrs from midnight of 1st day of data)
    Gridsize:     surface grid size, default is 21x21 grid of points
    what:         string with what to save from: BS, MP
    where:        string of directory path to save output
    '''
    from datetime import datetime
    import kamodo_ccmc.flythrough.model_wrapper as MW
    import json

    # Error checks
    surfaces = ['MP', 'BS', 'CS']
    if what in surfaces:
        # Extract surface
        x,y,z = gmComputeSurface(ko,timeHrs,Gridsize=Gridsize,what=what)
    else:
        print('ERROR, unknown surface in gmSaveSurface.')
        return
    Tpts = x.shape[0]
    if Tpts < 3:
        # An error occured
        print('ERROR, ',what,' surface could not be computed.')
        return False

    # Set some metadata
    title = 'ERROR'
    if what == 'MP':
        title = 'Magnetopause Surface Extraction'
    elif what == 'BS':
        title = 'Bow Shock Surface Extraction'
    elif what == 'CS':
        title = 'Tail Current Sheet Surface Extraction'
    model = ko.modelname
    coord = MW.Variable_Search('', model=model, return_dict=True)['v_x'][2]
    DT = datetime.utcfromtimestamp( ko.filedate.timestamp() + timeHrs*3600. )
    DateStr = DT.strftime("%Y-%m-%dT%H:%M:%SZ")
    FileStr = DT.strftime("_%Y%m%d_%H%M%S")
    DTnow = datetime.now()
    NowStr = DTnow.strftime("%Y-%m-%dT%H:%M:%SZ")
    # Save file
    filename = where+'/'+what+FileStr+'.json'
    data = {}
    data['Title'] = title
    data['ExtractDate'] = NowStr
    data['RunName'] = runname
    data['Model'] = model
    data['CoordinateSystem'] = coord
    data['ExtractedTime'] = DateStr
    data['ExtractedHours'] = timeHrs
    data['Gridsize'] = Gridsize
    data['DataArrays'] = ['x','y','z']
    # Positions rounded to 4 digits after decimal. Must make sure
    #   type is not object before around(). Convert to list
    data['x'] = np.around(x.astype(np.double),4).tolist()
    data['y'] = np.around(y.astype(np.double),4).tolist()
    data['z'] = np.around(z.astype(np.double),4).tolist()
    with open(filename, "w") as f:
        json.dump(data, f)

    return

### ====================================================================================== ###
def gmLoadSurface(sfile):
    '''
    Function to load computed grid points defining a surface.

    sfile:  string of relative directory path and filename to load surface
    '''
    import numpy as np
    from os.path import isfile,isdir
    import json

    # Error checks
    e = np.full((1,1), None)
    if not isfile(sfile):
        print('Error, file does not exist. ',sfile)
        return False,e,e,e

    # Load json file
    with open(sfile, 'r') as f:
        contents = json.load(f)
    x = np.array(contents['x'])
    y = np.array(contents['y'])
    z = np.array(contents['z'])

    return True,x,y,z

### ====================================================================================== ###
def gmLoadBtraces(sfile):
    '''
    Function to load computed B field traces.

    sfile:  string of relative directory path and filename to load surface
    '''
    import numpy as np
    from os.path import isfile,isdir
    import json

    # Error checks
    e = np.full((1,1), None)
    if not isfile(sfile):
        print('Error, file does not exist. ',sfile)
        return False,e,e,e

    # Load json file
    with open(sfile, 'r') as f:
        contents = json.load(f)
    x = np.array(contents['x'])
    y = np.array(contents['y'])
    z = np.array(contents['z'])
    #print(contents.keys())

    return True,x,y,z

### ====================================================================================== ###
def BlinesFig(fullfile, showE = True): 
    '''
    Function to take a specially generated file of last closed fieldline traces
      and generate a plotly figure of it. This can be displayed as is or added
      to other figures.
      
    Arguments:    
      fullfile:  The full file path to the file with fieldline extractions
      showE:     A logical for showing an Earth sphere in the figure
    '''
    import numpy as np
    import pandas as pd
    import plotly.graph_objects as go
    from os.path import isfile,isdir
    import json

    # Get filename for labels
    file = fullfile.split('/')[-1]

    # Make empty figure
    fig = go.Figure()

    # Load json file
    with open(fullfile, 'r') as f:
        contents = json.load(f)
    x = np.array(contents['x']).astype(float)
    y = np.array(contents['y']).astype(float)
    z = np.array(contents['z']).astype(float)
    nlines = contents['nlines']
    typ = contents['type']
    typset = set(typ)
    colors = ["#a1a1a1", "#e1a1a1", "#a1e1a1", "#a1a1e1"]
    k = 23
    for j in typset:
        label = 'B '+j
        x3 = np.array([None])
        y3 = np.array([None])
        z3 = np.array([None])
        for i in range(nlines):
            if typ[i] == j:
                x1 = np.array(x[i,:])
                y1 = np.array(y[i,:])
                z1 = np.array(z[i,:])
                x2 = x1[np.logical_not(np.isnan(x1))]
                y2 = y1[np.logical_not(np.isnan(y1))]
                z2 = z1[np.logical_not(np.isnan(z1))]
                x3 = np.concatenate([x3, x2, [None]])
                y3 = np.concatenate([y3, y2, [None]])
                z3 = np.concatenate([z3, z2, [None]])
        fig.add_trace(go.Scatter3d(x=x3, y=y3, z=z3, mode='lines',
                                   line=dict(color=colorset[k], width=3), 
                                   hoverinfo='skip', name=label, showlegend=True))
        k -= 1

    # Create Earth sphere
    if showE:
        elon = np.linspace(-180, 180, 181)
        elat = np.linspace(-90, 90, 91)
        elon_mg, elat_mg = np.meshgrid(np.array(elon), np.array(elat))
        ex = -1.*(np.cos(elat_mg*np.pi/180.)*np.cos(elon_mg*np.pi/180.))
        ey = -1.*(np.cos(elat_mg*np.pi/180.)*np.sin(elon_mg*np.pi/180.))
        ez = +1.*(np.sin(elat_mg*np.pi/180.))
        colorse = np.zeros(shape=ex.shape)
        colorse[ex<0.]=1  # Option to make day side a lighter color
        colorscalee = ['rgb(199,199,199)', 'rgb(0,0,0)']
        fig.add_surface(x=ex, y=ey, z=ez, surfacecolor=colorse,
                        cmin=0, cmax=1, colorscale=colorscalee, name='Earth',
                        showlegend=False, showscale=False, hoverinfo='skip')

    # Add hidden plot ranges (may need to adjust later)
    xx = np.full((3), None)
    yy = np.full((3), None)
    zz = np.full((3), None)
    xx[0], xx[2] = -75.,  15.
    yy[0], yy[2] =  25., -25.
    zz[0], zz[2] = -20.,  20.
    fig.add_trace(go.Scatter3d(x=xx, y=yy, z=zz, mode='markers', name='bounds',
        marker=dict(color="black", size=1, opacity=0.1), showlegend=False, hoverinfo='skip' ))

    # Set 3D view
    camera = dict(
        up=dict(x=0, y=0, z=1),
        eye=dict(x=1.1, y=0.85, z=0.25)
    )
    fig.update_layout(title=file, scene_camera=camera)
    fig.update_layout(scene_aspectmode='data')

    return fig

### ====================================================================================== ###
def BlinesMovie(where='.', showSW=False, showDate=''):
    '''
    Function to take a directory of extracted magnetic fieldlines and make a movie from them
      Limited to a max of 999 files at a time

    Arguments:
      where     Relative directory to get json files from
      showSW    Logical to hide or show solar wind field lines
      showDate  If specified, will only show values from this date
    '''
    import os
    import numpy as np
    import plotly.graph_objects as go
    import glob
    import time

    tic = time.perf_counter()

    if showDate == '':
        sfiles = glob.glob(where+'/*.json')
    else:
        sfiles = glob.glob(where+'/*'+showDate+'*.json')
    sfiles.sort()
    nF = len(sfiles)
    print(' ... found ',nF,' files ...')
    if nF > 999:
        print('ERROR, too many files to process, exiting.')
        return

    # make figure
    fig_dict = {
        "data": [],
        "layout": {},
        "frames": []
    }

    # bounds initialization
    minX, maxX, minY, maxY, minZ, maxZ = 0., 0., 0., 0., 0., 0.

    # fill in most of layout
    fig_dict["layout"]["scene_aspectmode"] = "data"
    fig_dict["layout"]["title"] = {"text": "Magnetic Fieldline Animation"}
    fig_dict["layout"]["margin"] = {"l": 0, "t": 25, "b": 5}
    fig_dict["layout"]["updatemenus"] = [
        {
            "buttons": [
                {
                    "args": [None, {"frame": {"duration": 0, "redraw": True},
                                    "fromcurrent": True,
                                    "transition": {"duration": 0}}],
                    "label": "Play",
                    "method": "animate"
                },
                {
                    "args": [[None], {"frame": {"duration": 0, "redraw": True},
                                     "mode": "immediate",
                                     "transition": {"duration": 0}}],
                    "label": "Pause",
                    "method": "animate"
                }
            ],
            "direction": "left",
            "pad": {"r": 10, "t": 25},
            "showactive": False,
            "type": "buttons",
            "x": 0.16,
            "xanchor": "right",
            "y": 0,
            "yanchor": "top"
        }
    ]

    sliders_dict = {
        "active": 0,
        "yanchor": "top", "xanchor": "left",
        "currentvalue": {
            "prefix": "Currently showing step: ",
            "visible": True,
            "xanchor": "left" },
        "transition": {"duration": 0},
        "pad": {"b": 10, "t": 5},
        "len": 0.83,
        "x": 0.17,
        "y": 0,
        "steps": []
    }

    # make frames
    ii = 0
    for sfile in sfiles:
        bname = os.path.basename(sfile)
        bname = bname[:-5]
        p1 = bname.split("_B_")
        p = p1[1]
        DateStr = p[0:4]+'-'+p[4:6]+'-'+p[6:8]+' '+p[9:11]+':'+p[11:13]+':'+p[13:15]+' UT'
        frame = {"data": [], "name": sfile}
        fig1 = pf.BlinesFig(sfile, showE=False)
        for i in range(len(fig1.data)):
            if len(fig1.data[i]['x'].shape) == 1:
                minX = min(minX, min(v for v in fig1.data[i]['x'] if v is not None))
                minY = min(minY, min(v for v in fig1.data[i]['y'] if v is not None))
                minZ = min(minZ, min(v for v in fig1.data[i]['z'] if v is not None))
                maxX = max(maxX, max(v for v in fig1.data[i]['x'] if v is not None))
                maxY = max(maxY, max(v for v in fig1.data[i]['y'] if v is not None))
                maxZ = max(maxZ, max(v for v in fig1.data[i]['z'] if v is not None))
        fig1.add_trace(go.Scatter3d(x=[0,0], y=[0,0], z=[0,0], mode='markers', name=DateStr,
            marker=dict(color="black", size=1, opacity=0.1), showlegend=True, hoverinfo='skip' ))
        frame["data"].append(fig1.data[len(fig1.data)-1])
        for i in range(len(fig1.data)):
            if 'B ' in fig1.data[i]['name']:
                if 'solar wind' not in fig1.data[i]['name']:
                    frame["data"].append(fig1.data[i])
                else:
                    if showSW: frame["data"].append(fig1.data[i])
        fig_dict["frames"].append(frame)
        slider_step = {"args": [[sfile],
            {"frame": {"duration": 0, "redraw": True},
             "mode": "immediate",
             "transition": {"duration": 0}}],
             "label": ii,
             "method": "animate"}
        sliders_dict["steps"].append(slider_step)
        ii += 1

    fig_dict["layout"]["sliders"] = [sliders_dict]

    # make data
    xx = np.full((3), None)
    yy = np.full((3), None)
    zz = np.full((3), None)
    xx[0], xx[2] = minX, maxX
    yy[0], yy[2] = minY, maxY
    zz[0], zz[2] = minZ, maxZ
    fig1 = pf.BlinesFig(sfiles[0], showE=True)
    fig1.add_trace(go.Scatter3d(x=xx, y=yy, z=zz, mode='markers', name='newbounds',
        marker=dict(color="black", size=1, opacity=0.1), showlegend=False, hoverinfo='skip' ))

    fig_dict["data"].append(fig1.data[len(fig1.data)-1])  # First one gets replaced with animation (make small)
    for i in range(len(fig1.data)):
        if 'Earth' in fig1.data[i]['name']: fig_dict["data"].append(fig1.data[i])  # Earth
    fig_dict["data"].append(fig1.data[len(fig1.data)-1])  # Bounds
    fig_dict["data"].append(fig1.data[len(fig1.data)-1])  # Bounds
    fig_dict["data"].append(fig1.data[len(fig1.data)-1])  # Bounds

    # make final figure
    fig = go.Figure(fig_dict)

    # print out time to make movie
    toc = time.perf_counter()
    print(f"Creation time: {toc - tic:0.4f} seconds")

    return fig

### ====================================================================================== ###
def gmSurfaceMovie(where='.', where2='', wireframe=False):
    import os
    import numpy as np
    import plotly.graph_objects as go
    import glob

    sfiles = glob.glob(where+'/*.json')
    sfiles.sort()
    if where2 != '':
        sfiles2 = glob.glob(where2+'/*.json')
        sfiles2.sort()

    # make figure
    fig_dict = {
        "data": [],
        "layout": {},
        "frames": []
    }

    # make data
    fig1,success = gmGetSurfacePlot(sfile=sfiles[0], wireframe=wireframe)
    fig_dict["data"].append(fig1.data[0])  # this one gets replaced with animation
    fig_dict["data"].append(fig1.data[1])
    fig_dict["data"].append(fig1.data[2])

    # fill in most of layout
    fig_dict["layout"]["scene_aspectmode"] = "data"
    fig_dict["layout"]["title"] = {"text": "Surface Extraction Animation"}
    fig_dict["layout"]["margin"] = {"l": 0, "t": 25, "b": 5}
    fig_dict["layout"]["updatemenus"] = [
        {
            "buttons": [
                {
                    "args": [None, {"frame": {"duration": 0, "redraw": True},
                                    "fromcurrent": True, 
                                    "transition": {"duration": 0}}],
                    "label": "Play",
                    "method": "animate"
                },
                {
                    "args": [[None], {"frame": {"duration": 0, "redraw": True},
                                     "mode": "immediate",
                                     "transition": {"duration": 0}}],
                    "label": "Pause",
                    "method": "animate"
                }
            ],
            "direction": "left",
            "pad": {"r": 10, "t": 25},
            "showactive": False,
            "type": "buttons",
            "x": 0.16,
            "xanchor": "right",
            "y": 0,
            "yanchor": "top"
        }
    ]

    sliders_dict = { 
        "active": 0,
        "yanchor": "top", "xanchor": "left",
        "currentvalue": {
            "prefix": "Currently showing step: ",
            "visible": True,
            "xanchor": "left" },
        "transition": {"duration": 0},
        "pad": {"b": 10, "t": 5},
        "len": 0.83,
        "x": 0.17,
        "y": 0,
        "steps": []
    }

    # make frames
    ii = 0
    for sfile in sfiles:
        bname = os.path.basename(sfile)
        bname = bname[:-5]
        p = bname.split("_")
        DateStr = p[1][0:4]+'/'+p[1][4:6]+'/'+p[1][6:8]+\
            ' '+p[2][0:2]+':'+p[2][2:4]+':'+p[2][4:6]+' UT'
        frame = {"data": [], "name": sfile}
        fig1,success = gmGetSurfacePlot(sfile=sfile, wireframe=wireframe)
        if where2 != '':
            fig2,success = gmGetSurfacePlot(sfile=sfiles2[ii], wireframe=wireframe)
        frame["data"].append(fig1.data[0])
        frame["data"].append(fig2.data[0])
        fig_dict["frames"].append(frame)
        slider_step = {"args": [[sfile],  
            {"frame": {"duration": 0, "redraw": True},
             "mode": "immediate",
             "transition": {"duration": 0}}],
             "label": ii,
             "method": "animate"}
        sliders_dict["steps"].append(slider_step)
        ii += 1

    fig_dict["layout"]["sliders"] = [sliders_dict]

    fig = go.Figure(fig_dict)

    return fig

### ====================================================================================== ###
def gm2DSliceFig(ko, timeHrs=1., var='P', pco='GSM', slicedir='Z', sliceval=0., 
                runname='UNKNOWN', resolution=0.25, showBS=True, showMP=True,
                xrange=[-95., 25.], yrange=[-30., 30.], crange='',
                xtic=10., ytic=10., log10=False, logCustom=True, 
                colorscale='Viridis'):
    '''
    Function to create a 2D slice for a Y or Z constant value. It will interpolate
      from a Kamodo object onto that grid in given coordinate system.
    Returns a full 2D plotly figure.
    Many options to customize the figure are available.
    
    ko:           Kamodo object
    timeHrs:      floating time in hours from start of first day of data
    var:          string variable name
    pco:          plot coordinate system
    slicedir:     direction of slice, Y or Z
    sliceval:     value of slice position
    runname:      string that describes run name
    resolution:   grid dx,dy in R_E of plot positions
    showMP:       logical to show magnetopause boundary (requires 'status' variable)
    showBS:       logical to show bow shock (requires 'v_x' variable)
    crange:       2 value array for min/max contour range
    xrange:       2 value array for min/max extent of X values in slice
    yrange:       2 value array for min/max extent of Y values in slice
      -note: xy range values may be ajusted to keep overall plot aspect ratio
    xtic:         tic spacing on x axis
    ytic:         tic spacing on y axis
    log10:        logical to use log10 of variable values
    logCustom:    logical to adjust label on colorscale to more human interpretable values
    colorscale:   string value of colorscale range from official plotly list
    '''

    import time
    import numpy as np
    from datetime import datetime
    from kamodo_ccmc.tools.plotfunctions import XYC
    from kamodo_ccmc.flythrough.utils import ConvertCoord
    import kamodo_ccmc.flythrough.model_wrapper as MW

    tic = time.perf_counter()

    # Set some constants
    pcot = 'car'
    slicedir = slicedir.upper()
    if slicedir != 'Y' and slicedir != 'Z':
        print('ERROR: Invalid slice passed.')
        return

    # Pad range to make plot nice for 800x400 size
    xmin, xmax = xrange[0], xrange[1]
    ymin, ymax = yrange[0], yrange[1]
    ratio = (xmax - xmin)/(ymax - ymin)
    if ratio > 2.05:
        # extend Y
        pad = 0.5*(((xmax - xmin)/2.05)-(ymax - ymin))
        ymin += -pad
        ymax += pad
    elif ratio < 2.05:
        # extend X
        pad = 0.5*(-(2.05*(ymax-ymin)-xmax)-xmin)
        xmin += -pad
        xmax += pad
    xmin = int(xmin-6.)
    xmax = int(xmax+6.)
    ymin = int(ymin-3.)
    ymax = int(ymax+3.)
    xpts = 1+int((xmax-xmin)/resolution)
    ypts = 1+int((ymax-ymin)/resolution)

    # Get current coordinate system and type from Kamodo object
    tmp = list(MW.Model_Variables(model=ko.modelname, return_dict=True).values())[0][2:4]
    co, cot = tmp[0], tmp[1]
    print('Coords:',co,cot)

    # Setup interpolator and some labels
    vMag = False  
    if var not in ko:
        if var+'_x' in ko:
            vMag = True
            interpX = getattr(ko, var+'_x')
            interpY = getattr(ko, var+'_y')
            interpZ = getattr(ko, var+'_z')
            vunits = ko.variables[var+'_x']['units']
            var2 = var+'_x_ijk'
            cr = MW.Coord_Range(ko, [var2], return_dict=True, print_output=False)
            xunits = cr[var2]['X'][2]
            coord = MW.Variable_Search('', model=ko.modelname, return_dict=True)[var+'_x'][2]
        else:
            print('Error, variable not in Kamodo object')
            return
    else:
        interp = getattr(ko, var)
        vunits = ko.variables[var]['units']
        var2 = var+'_ijk' 
        cr = MW.Coord_Range(ko, [var2], return_dict=True, print_output=False)
        xunits = cr[var2]['X'][2] 
        coord = MW.Variable_Search('', model=ko.modelname, return_dict=True)[var][2]
    varlabel = var+" ["+vunits+"]"
    if 'v_x' in ko:
        interpBS = getattr(ko, 'v_x') 
    else:
        showBS = False
    if 'status' in ko:
        interpMP = getattr(ko, 'status') 
    else:
        showMP = False

    # Make base grid
    x = np.linspace(xmin, xmax, xpts)
    y = np.linspace(ymin, ymax, ypts)
    xx, yy = np.array(np.meshgrid(x, y, indexing='ij'))
    c = xx
    x1 = np.reshape(xx, -1)
    xlabel = 'X ['+xunits+'] '+pco
    if slicedir == 'Z':
        y1 = np.reshape(yy, -1)
        z1 = np.full((len(x1)), sliceval)
        ylabel = 'Y ['+xunits+'] '+pco
    else:
        z1 = np.reshape(yy, -1)
        y1 = np.full((len(x1)), sliceval)
        ylabel = 'Z ['+xunits+'] '+pco
    t1 = np.full((len(x1)), timeHrs)

    # Take grid in preferred plot coordinates and transform into model coordinates
    if co == 'GDZ':
        print('Fix this edge case later')
        return
    elif co == pco and cot == pcot:
        # same
        g2 = np.stack((t1, x1, y1, z1), axis=-1)  # nx4 grid
    else:
        # transform
        x2, y2, z2, units = ConvertCoord(t1, x1, y1, z1, pco, pcot, co, cot)
        g2 = np.stack((t1, x2, y2, z2), axis=-1)  # nx4 grid


    # Interpolate
    if vMag:
        valueX = interpX(g2)
        valueY = interpY(g2)
        valueZ = interpZ(g2)
        value = np.sqrt(valueX**2 + valueY**2 + valueZ**2)
    else:
        value = interp(g2)
    if log10:
        value[value <= 0.] = np.nan
        value = np.log10(value)
        if logCustom:
            varlabel = "log scale<br>"+varlabel
        else:
            varlabel = "log10("+varlabel+")"

    # Set values inside R=2.25 to None
    r1 = np.sqrt(x1**2 + y1**2 + z1**2)
    mask = r1 < 2.25
    value[mask] = None
    vmin = np.nanmin(value)
    vmax = np.nanmax(value)
    c = np.reshape(value, (len(x),len(y)))

    # Bow shock
    if showBS:
        valueVX = interpBS(g2)
        vx = np.reshape(valueVX, (len(x),len(y)))
        bsx = np.full((len(y)), None)
        bsy = np.full((len(y)), None)
        for j in range(len(y)):
            bsy[j] = y[j]
            refVX = vx[-1,j]
            for ii in range(len(x)):
                i = len(x)-1-ii
                if x[i] < 1.:
                    break
                if vx[i,j] > (0.85*refVX):
                    bsx[j] = x[i]
                    break

    # Magnetopause
    if showMP:
        valueMP = interpMP(g2)
        status = np.reshape(valueMP, (len(x),len(y)))
        mpx = np.full((len(y)), None)
        mpy = np.full((len(y)), None)
        for j in range(len(y)):
            mpy[j] = y[j]
            for ii in range(len(x)):
                i = len(x)-1-ii
                if x[i] < 1.:
                    break
                if status[i,j] > 2.5:
                    mpx[j] = x[i]
                    break

    # Plot values
    pTS = ko.filedate.timestamp() + timeHrs*3600.
    pDT = datetime.utcfromtimestamp(pTS)
    pDateStr = pDT.strftime("%Y-%m-%d %H:%M:%S UT")
    label1 = slicedir+' = '+str(sliceval)+' slice at  '+pDateStr
    label2 = 'Model: '+ko.modelname
    ptitle = 'CCMC Run: '+runname

    fig = XYC(xlabel,x,ylabel,y,varlabel,c,title=ptitle,crange=crange)
    fig.update_traces(colorscale=colorscale)
    if log10 and logCustom:
        fig.update_traces(colorbar=dict(
            tickvals=[-5,-4,-3,-2,-1,0,1,2],
            ticktext=['.00001','.0001','.001','.01','.1','1.','10.','100.'],
        ))

    if showBS:
        fig.add_scatter(x=bsx, y=bsy, showlegend=False, hoverinfo='skip',
                       line=dict(color="grey", width=2))
    if showMP:
        fig.add_scatter(x=mpx, y=mpy, showlegend=False, hoverinfo='skip',
                       line=dict(color="grey", width=2))
    degs = np.linspace(-180, 180, 181)
    ex = np.cos(degs*np.pi/180.)
    ey = np.sin(degs*np.pi/180.)
    if abs(sliceval) < 2.5:
        r = np.sqrt(2.5**2 - sliceval**2)
        fig.add_scatter(x=r*ex, y=r*ey, fill='toself', showlegend=False, hoverinfo='skip',
                       fillcolor="black", line=dict(color="black", width=1))
    if sliceval == 0.0:
        # Only show when slicing through 0
        r = np.sqrt(1. - sliceval**2)
        fig.add_scatter(x=r*ex, y=r*ey, showlegend=False, hoverinfo='skip',
                       fillcolor="black", line=dict(color="white", width=1))
        degs = np.linspace(-90, 90, 91)
        ex = np.cos(degs*np.pi/180.)
        ey = np.sin(degs*np.pi/180.)
        fig.add_scatter(x=r*ex, y=r*ey, fill='toself', showlegend=False, hoverinfo='skip',
                       fillcolor="white", line=dict(color="white", width=1))

    fig.update_layout(
        xaxis=dict(scaleanchor="y", scaleratio=1, tick0 = 0., dtick = xtic,
                   range=[xrange[0], xrange[1]]),
        yaxis=dict(scaleanchor="x", scaleratio=1, tick0 = 0., dtick = ytic,
                   range=[yrange[0], yrange[1]]),
        annotations=[
            dict(text=label1, x=1.18, y=1.07, ax=0, ay=0, xanchor="right",
                 xshift=0, yshift=0, xref="paper", yref="paper",
                 font=dict(size=16, family="sans serif", color="#000000")),
            dict(text=label2, x=-0.07, y=-0.1, ax=0, ay=0, xanchor="left",
                 xshift=0, yshift=0, xref="paper", yref="paper",
                 font=dict(size=16, family="sans serif", color="#000000")),
        ],
        width=800, height=400,
    )
    fig.update_layout(margin=dict(r=120))
    fig.data[0]['colorbar']['xpad'] = 5
    fig.data[0]['colorbar']['y'] = -.03
    fig.data[0]['colorbar']['yanchor'] = 'bottom'

    toc = time.perf_counter()
    print(f"  Time: {toc - tic:0.4f} seconds")

    return fig

