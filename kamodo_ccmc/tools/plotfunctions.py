def toLog10(fig):
    """
    Quick function to take a 2D contour figure and make the contour log10 scale.
    Pass in a plotly figure and it will return an updated plotly figure.
    """
    
    import numpy as np
    
    # grab values from old plot
    val = fig.data[0]['z']
    # set negative and zero values to NaN, compute log10 of values, NaN will stay NaN
    val[val <= 0.] = np.nan
    val = np.log10(val)
    # assign back to plot object
    fig.data[0]['z'] = val
    # add log10 to colorbar label
    newlabel = 'log10<br>'+fig.data[0]['colorbar']['title']['text']
    fig.data[0]['colorbar']['title']['text'] = newlabel
    return fig


def XYC(Xlabel, X, Ylabel, Y, Clabel, C, title='Plot Title',
        colorscale='Viridis', crange='',):
    """
    Simple 2D plot. Send in array of X points, Y points, and 2d array of values
      with labels for each and a basic plot is created.
    
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
    def plot_2D(X = X, Y = Y):
        return C
    plot1 = Kamodo(plot_2D = plot_2D)
    fig = plot1.plot(plot_2D = dict())

    # Set colorscale
    if colorscale == "BlueRed":
        fig.update_traces(colorscale="RdBu", reversescale=True)
    elif colorscale == "Rainbow":
        fig.update_traces(
            colorscale=[[0.00, 'rgb(0,0,255)'],
                        [0.25, 'rgb(0,255,255)'],
                        [0.50, 'rgb(0,255,0)'],
                        [0.75, 'rgb(255,255,0)'],
                        [1.00, 'rgb(255,0,0)']]
        )
    else:
        fig.update_traces(colorscale=colorscale)

    # Set plot options
    fig.update_traces(
        zmin=cmin, zmax=cmax, colorbar=dict(title=Clabel, tickformat=".3g"),
        ncontours=201, contours=dict(coloring="fill",showlines=False),
        hovertemplate=Xlabel+": %{x:.4g}<br>"+Ylabel+": %{y:.4g}<br>"+
            Clabel+": %{z:.6g}<br>"+"<extra></extra>"
    )
    fig.update_layout(
        scene_aspectmode='data',
        title=dict(text=title, 
                   yref="container", yanchor="top", x=0.01, y=0.97,
                   font=dict(size=16, family="sans serif", color="#000000")),
        margin=dict(l=0), xaxis_title=Xlabel, yaxis_title=Ylabel
    )

    return fig


def ReplotLL3D(figIn,model,altkm,plotTime,plotCoord='GEO',
                 title='Plot Title',colorscale='Viridis',crange='',
                 opacity=0.70,axis=True,debug=0):
    """
    Takes a gridified 2D lon/lat figure and creates new plots
    in 3D various coordinate systems.
    """

    import numpy as np
    import pandas as pd
    import datetime as dt
    import pytz
    from kamodo import Kamodo
    import kamodo_ccmc.flythrough.model_wrapper as MW
    from kamodo_ccmc.flythrough.utils import ConvertCoord
    from kamodo_ccmc.tools.shoreline import shoreline

    # Get current coordinate system and type. Create string for labels
    tmp = list(MW.Model_Variables(model,return_dict=True).values())[0][2:4]
    co = tmp[0]
    cot = tmp[1]
    costr = co+'('+cot+')'

    # Pull out lon, lat, values and min/max from passed in figure
    lon = figIn.data[0]['x']
    lat = figIn.data[0]['y']
    val = figIn.data[0]['z']
    val2 = np.reshape(val,(len(lat),len(lon)))
    varn = figIn.data[0]['colorbar']['title']['text']
    cmin = np.min(val)
    cmax = np.max(val)
    if crange != '':
        cmin = float(crange[0])
        cmax = float(crange[1])

    # Prepare variables for use later
    nlon = len(lon)
    nlat = len(lat)
    npts = (nlon*nlat)
    lon_mg, lat_mg = np.meshgrid(np.array(lon), np.array(lat))
    x = lon_mg
    y = lat_mg
    z = np.full(x.shape, altkm)
    t = np.full(x.shape, plotTime)
    full_x = np.reshape(x,-1)
    full_y = np.reshape(y,-1)
    full_z = np.reshape(z,-1)
    full_t = np.reshape(t,-1)
    rscale = (altkm + 6.3781E3)/6.3781E3
    ilon = np.linspace(-180, 180, 181)
    ilat = np.linspace(-90, 90, 91)

    #Convert incoming coordinates into plot coordinages (cartesian)
    if co == 'GDZ':
        # GDZ does not convert into other coordinates well, first go to GEO-car
        xx,yy,zz,units = ConvertCoord(full_t,full_x,full_y,full_z,co,cot,'GEO','car')
        full_x = xx
        full_y = yy
        full_z = zz
    if plotCoord != 'GEO':
        xx,yy,zz,units = ConvertCoord(full_t,full_x,full_y,full_z,'GEO','car',plotCoord,'car')
    x = np.reshape(xx,(len(lat),len(lon)))
    y = np.reshape(yy,(len(lat),len(lon)))
    z = np.reshape(zz,(len(lat),len(lon)))

    # Generate initial figure to build upon
    def plot3d_var(X = x, Y = y, Z = z):
        return val2
    kobject = Kamodo(plot_var = plot3d_var)
    fig = kobject.plot(plot_var = dict())

    # Set colorscale
    if colorscale == "BlueRed":
        fig.update_traces(colorscale="RdBu", reversescale=True)
    elif colorscale == "Rainbow":
        fig.update_traces(
            colorscale=[[0.00, 'rgb(0,0,255)'],
                        [0.25, 'rgb(0,255,255)'],
                        [0.50, 'rgb(0,255,0)'],
                        [0.75, 'rgb(255,255,0)'],
                        [1.00, 'rgb(255,0,0)']]
        )
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
        hovertemplate="<b>"+model+"</b><br>"+
            "X: %{x:.2f} "+units[0]+" "+plotCoord+"<br>"
            "Y: %{y:.2f} "+units[1]+" "+plotCoord+"<br>"
            "Z: %{z:.2f} "+units[2]+" "+plotCoord+"<br>"
            "Lon: %{customdata[1]:.2f} "+co+"<br>"+
            "Lat: %{customdata[2]:.2f} "+co+"<br>"+
            "Alt: "+"{:.2f}".format(altkm)+" km<br>"+
            varn+": %{customdata[0]:.4g} <br>"+
            "<extra></extra>"
    )
    if not axis:
        fig.update_scenes(xaxis=dict(visible=False),
                          yaxis=dict(visible=False),
                          zaxis=dict(visible=False))
    timestr = dt.datetime.fromtimestamp(plotTime, tz = pytz.utc).strftime("%Y/%m/%d %H:%M:%S")
    lltxt = model+',  '+plotCoord+' Coordinates,  '+str(altkm)+' km Altitude,  '+timestr
    camera = dict(up=dict(x=-1, y=0, z=0.),
                  eye=dict(x=0., y=0., z=1.25))
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
        #scene_camera=camera,
        width=750,
        height=450,
    )

    # Add poles to plot
    xt = [0., 0.]
    yt = [0., 0.]
    zt = [-1.2, 1.2]
    fig.add_scatter3d(mode='lines',x=xt,y=yt,z=zt,line=dict(width=4,color='black'),
                      showlegend=False,hovertemplate=plotCoord+' Pole<extra></extra>')
    if plotCoord != 'GEO':
        # Add pole for GEO as well
        xt = np.array([ 0., 0.])
        yt = np.array([ 0., 0.])
        zt = np.array([-1., 1.])
        tt = plotTime + 0.*xt
        xt,yt,zt,un = ConvertCoord(tt,xt,yt,zt,'GEO','car',plotCoord,'car')
        xt = 1.2*xt
        yt = 1.2*yt
        zt = 1.2*zt
        fig.add_scatter3d(mode='lines',x=xt,y=yt,z=zt,line=dict(width=4,color='rgb(79,79,79)'),
                          showlegend=False,hovertemplate='GEO Pole<extra></extra>')

    # Zero longitude
    #xt = rscale*(np.cos(ilat*np.pi/180.))
    #yt = 0.*ilat
    #zt = rscale*(np.sin(ilat*np.pi/180.))
    #fig.add_scatter3d(mode='lines',x=xt,y=yt,z=zt,line=dict(width=2,color='black'),
    #                  showlegend=False,hovertemplate=plotCoord+' Longitude=0<extra></extra>')

    # Zero latitude
    xt = rscale*(np.cos(ilon*np.pi/180.))
    yt = rscale*(np.sin(ilon*np.pi/180.))
    zt = 0.*ilon
    fig.add_scatter3d(mode='lines',x=xt,y=yt,z=zt,line=dict(width=2,color='black'),
                      showlegend=False,hovertemplate=plotCoord+' Latitude=0<extra></extra>')

    # Create blank earth surface
    elon = np.linspace(-180, 180, 181)
    elat = np.linspace(-90, 90, 91)
    elon_mg, elat_mg = np.meshgrid(np.array(elon), np.array(elat))
    ex=-(np.cos(elat_mg*np.pi/180.)*np.cos(elon_mg*np.pi/180.))
    ey=-(np.cos(elat_mg*np.pi/180.)*np.sin(elon_mg*np.pi/180.))
    ez= (np.sin(elat_mg*np.pi/180.))
    colorse = np.zeros(shape=ex.shape)
    #colorse[ex<0.]=1  # Option to make day side a lighter color
    colorscalee = [ 'rgb(99,99,99)', 'rgb(0,0,0)']
    fig.add_surface(x=ex, y=ey, z=ez, surfacecolor=colorse, 
                    cmin=0, cmax=1, colorscale=colorscalee, 
                    showlegend=False, showscale=False, hoverinfo='skip')

    # Shoreline (land/water boundaries)
    pos = shoreline(rscale=1.001,coord=plotCoord,utcts=plotTime)
    fig.add_scatter3d(mode='lines', x=pos[:,0], y=pos[:,1], z=pos[:,2], 
                      line=dict(width=2,color='white'),
                      showlegend=False,hoverinfo='skip')

    return fig


def TransformedLL(interp,varname,model,plotTime,fhrs,altkm,
                  plotCoord='GEO',title='Plot Title',shoreline=False,
                  colorscale='Viridis',crange='',
                  debug=0):
    """Assumes spherical coords"""

    import numpy as np
    import pytz
    import datetime as dt
    from kamodo import Kamodo
    import kamodo_ccmc.flythrough.model_wrapper as MW
    from kamodo_ccmc.flythrough.utils import ConvertCoord
    from kamodo_ccmc.tools.shoreline import shoreline

    # Get current coordinate system and type from model. Create string for labels
    tmp = list(MW.Model_Variables(model,return_dict=True).values())[0][2:4]
    co = tmp[0]
    cot = tmp[1]
    costr = co+'('+cot+')'
    if co != 'GDZ' or cot != 'sph':
        print('Expecting coordinates to be GDZ(sph) but got ',costr)

    # Set grid to interpolate to
    ilon = np.linspace(-180, 180, 91)
    ilat = np.linspace(-90, 90, 73)
    ilon_mg, ilat_mg = np.meshgrid(np.array(ilon), np.array(ilat))
    x = np.reshape(ilon_mg,-1)
    y = np.reshape(ilat_mg,-1)
    # The z value 1 [Re] so that we can compute km altitude shift (GDZ uses 0)
    if plotCoord == 'GDZ':
        z = np.zeros(x.shape)
    else:
        z = np.ones(x.shape)
    t = np.full(x.shape, plotTime)
 
    # Convert from plotCoord coordinates into model coordinates to interpolate
    if co == 'GDZ':
        if plotCoord == 'GDZ':
            # No change to coordinates, so don't transform
            xx = x; yy = y; zz = z
        else:
            # GDZ only converts to/from GEO, so have to convert twice
            xx,yy,zz,units = ConvertCoord(t,x,y,z,plotCoord,'sph','GEO','sph')
            x = xx; y = yy; z = zz
            xx,yy,zz,units = ConvertCoord(t,x,y,z,'GEO','sph',co,cot)
    else:
        xx,yy,zz,units = ConvertCoord(t,x,y,z,plotCoord,'sph',co,cot)
    if debug > 0:
        # Plot altitude variable
        print('DEBUG:',debug)
        val = zz
        varname = 'Altitude[km]'
    else:
        zz = np.full(xx.shape, altkm)
        val = np.ndarray(shape=(len(xx)), dtype=np.float32)
        for i in range(len(xx)):
            val[i] = interp([fhrs, xx[i], yy[i], zz[i]])  # time, lon, lat, height
    val2 = np.reshape(val,(len(ilat),len(ilon)))
    cmin = np.min(val)
    cmax = np.max(val)
    if crange != '':
        cmin = float(crange[0])
        cmax = float(crange[1])

    def plot_2D(Lon = ilon, Lat = ilat):
        return val2
    plot1 = Kamodo(plot_2D = plot_2D)
    fig = plot1.plot(plot_2D = dict())

    # Set colorscale
    if colorscale == "BlueRed":
        fig.update_traces(colorscale="RdBu", reversescale=True)
    elif colorscale == "Rainbow":
        fig.update_traces(colorscale=[[0.00, 'rgb(0,0,255)'],
                                      [0.25, 'rgb(0,255,255)'],
                                      [0.50, 'rgb(0,255,0)'],
                                      [0.75, 'rgb(255,255,0)'],
                                      [1.00, 'rgb(255,0,0)']]
        )
    else:
        fig.update_traces(colorscale=colorscale)

    # Set plot options
    fig.update_traces(
        zmin=cmin, zmax=cmax, colorbar=dict(title=varname, tickformat=".3g"),
        ncontours=201, contours=dict(coloring="fill",showlines=False),
        hovertemplate="Lon: %{x:.2f}<br>"+"Lat: %{y:.2f}<br>"+
            varname+": %{z:.6g}<br>"+"<extra></extra>"
    )
    fig.update_yaxes(scaleanchor = "x", scaleratio = 1)
    fig.update_xaxes(tick0=0.,dtick=45.)
    fig.update_yaxes(tick0=0.,dtick=45.)
    timestr = dt.datetime.fromtimestamp(plotTime, tz = pytz.utc).strftime("%Y/%m/%d %H:%M:%S")
    lltxt = model+',  '+plotCoord+' Coordinates,  '+str(altkm)+' km Altitude,  '+timestr
    fig.update_layout(
        plot_bgcolor="white",
        scene_aspectmode='data',
        title=dict(text=title+'<br>'+lltxt, 
                   yref="container", yanchor="top", x=0.01, y=0.97,
                   font=dict(size=16, family="sans serif", color="#000000")),
        margin=dict(l=0),
        width=750,
        height=375,
    )

    # Shoreline (land/water boundaries)
    if plotCoord == 'GDZ':
        # Can't convert into GDZ, but will convert to GEO as we only use lon/lat
        pos = shoreline(coord='GEO',coordT='sph',utcts=plotTime)
    else:
        pos = shoreline(coord=plotCoord,coordT='sph',utcts=plotTime)
    fig.add_scattergl(mode='lines', x=pos[:,0], y=pos[:,1], hoverinfo='skip', 
                      line=dict(width=1,color='white'), showlegend=False)

    return fig

