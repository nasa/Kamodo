def toLog10(fig):
    """
    Function to take a 2D contour figure and make the contour log10 scale.
    Pass in a plotly figure and it will return an updated plotly figure.
    
    Arguments:
        fig   A plotly figure object
    
    Returns a plotly figure object.
    """
    
    import numpy as np
    
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
    return fig


def toColor(fig,colorscale='Viridis'):
    """
    Function to take a 2D contour figure from Kamodo and change the colorscale
    and set the number of contours to a larger number.
    
    Arguments:
        fig   A plotly figure object
        colorscale  Optional string name of colorscale. Takes custom 
                    colorscales RdBu, Rainbow, or standard python values
    
    Returns a plotly figure object.
    """

    # Set colorscale
    if colorscale == "BlueRed":
        fig.update_traces(colorscale="RdBu", reversescale=True)
    elif colorscale == "Rainbow":
        fig.update_traces(
            colorscale=[[0.00, 'rgb(0,0,255)'],
                        [0.25, 'rgb(0,255,255)'],
                        [0.50, 'rgb(0,255,0)'],
                        [0.75, 'rgb(255,255,0)'],
                        [1.00, 'rgb(255,0,0)']])
    else:
        fig.update_traces(colorscale=colorscale)
    fig.update_traces(ncontours=201, 
                      contours=dict(coloring="fill", showlines=False))
    
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
    def plot_2D(X = X, Y = Y):
        return C
    plot1 = Kamodo(plot_2D = plot_2D)
    fig = plot1.plot(plot_2D = dict())

    # Set colorscale
    fig = toColor(fig,colorscale=colorscale)

    # Set plot options
    fig.update_traces(
        zmin=cmin, zmax=cmax, colorbar=dict(title=Clabel, tickformat=".3g"),
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


def ReplotLL3D(figIn,model,altkm,plotts,plotCoord='GEO',
                 title='Plot Title',colorscale='Viridis',crange='',
                 opacity=0.70,axis=True,debug=0):
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
    tmp = list(MW.Model_Variables(model,return_dict=True).values())[0][2:4]
    co = tmp[0]
    cot = tmp[1]

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
    lon_mg, lat_mg = np.meshgrid(np.array(lon), np.array(lat))
    x = lon_mg
    y = lat_mg
    z = np.full(x.shape, altkm)
    t = np.full(x.shape, plotts)
    full_x = np.reshape(x,-1)
    full_y = np.reshape(y,-1)
    full_z = np.reshape(z,-1)
    full_t = np.reshape(t,-1)
    rscale = (altkm + 6.3781E3)/6.3781E3
    ilon = np.linspace(-180, 180, 181)
    #ilat = np.linspace(-90, 90, 91)

    #Convert incoming coordinates into plot coordinages (cartesian)
    if co == 'GDZ':
        # GDZ does not convert into other coordinates well, first go to GEO-car
        xx,yy,zz,units = ConvertCoord(full_t,full_x,full_y,full_z,
                                      co,cot,'GEO','car')
        full_x, full_y, full_z = xx,yy,zz
    if plotCoord != 'GEO':
        xx,yy,zz,units = ConvertCoord(full_t,full_x,full_y,full_z,
                                      'GEO','car',plotCoord,'car')
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
    plotdt = dt.datetime.fromtimestamp(plotts, tz = pytz.utc)
    timestr = plotdt.strftime("%Y/%m/%d %H:%M:%S")
    lltxt = model+',  '+plotCoord+' Coordinates,  '
    lltxt = lltxt+str(altkm)+' km Altitude,  '+timestr
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
    fig.add_scatter3d(mode='lines',x=xt,y=yt,z=zt,
                      line=dict(width=4,color='black'),
                      showlegend=False,
                      hovertemplate=plotCoord+' Pole<extra></extra>')
    if plotCoord != 'GEO':
        # Add pole for GEO as well
        xt = np.array([ 0., 0.])
        yt = np.array([ 0., 0.])
        zt = np.array([-1., 1.])
        tt = np.full(xt.shape, plotts)
        xt,yt,zt,un = ConvertCoord(tt,xt,yt,zt,'GEO','car',plotCoord,'car')
        xt = 1.2*xt
        yt = 1.2*yt
        zt = 1.2*zt
        fig.add_scatter3d(mode='lines',x=xt,y=yt,z=zt,
                          line=dict(width=4,color='rgb(79,79,79)'),
                          showlegend=False,
                          hovertemplate='GEO Pole<extra></extra>')

    # Zero latitude
    xt = rscale*(np.cos(ilon*np.pi/180.))
    yt = rscale*(np.sin(ilon*np.pi/180.))
    zt = np.zeros(xt.shape)
    fig.add_scatter3d(mode='lines',x=xt,y=yt,z=zt,
                      line=dict(width=2,color='black'),
                      showlegend=False,
                      hovertemplate=plotCoord+' Latitude=0<extra></extra>')

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
    pos = shoreline(rscale=1.001,coord=plotCoord,utcts=plotts)
    fig.add_scatter3d(mode='lines', x=pos[:,0], y=pos[:,1], z=pos[:,2], 
                      line=dict(width=2,color='white'),
                      showlegend=False,hoverinfo='skip')

    return fig


def GDZSlice4D(interp,varname,model,date,plotType,plotCoord='GEO',
               fixed_time='',fixed_lon='',fixed_lat='',fixed_alt='', 
               title='Plot Title',shoreline=False,colorscale='Viridis',
               crange='',logscale=False):
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
    from kamodo import Kamodo
    import kamodo_ccmc.flythrough.model_wrapper as MW
    from kamodo_ccmc.flythrough.utils import ConvertCoord
    from kamodo_ccmc.tools.shoreline import shoreline

    tic = time.perf_counter()

    # Check passed in plotType and parse for use later
    Types = ['Lon', 'Lat', 'Alt', 'Time']
    axis = plotType.split('-')
    if len(axis) != 2:
        print("ERROR, invalid plotType passed."); return
    if axis[0] not in Types or axis[1] not in Types:
        print("ERROR, plotType elements invalid."); return

    # Set base timestamp to start of day passed in.
    basets = date.replace(hour=0,minute=0,second=0,microsecond=0).timestamp()
    if 'Time' not in axis:
        slicets = basets + fixed_time*3600.

    # Get current coordinate system and type from model. Create label strings
    tmp = list(MW.Model_Variables(model,return_dict=True).values())[0][2:4]
    co = tmp[0]
    cot = tmp[1]
    costr = co+'('+cot+')'
    if co != 'GDZ' or cot != 'sph':
        print('Expecting coordinates to be GDZ(sph) but got ',costr)

    # Set itim array 
    if 'Time' in axis:
        itim = np.linspace(0, 24, 97)
        if axis[0] == 'Time': i1 = itim
        if axis[1] == 'Time': i2 = itim
        fixed_time = ''
    else:
        if fixed_time == '':  print("ERROR, fixed_time not passed."); return
        itim = np.array([fixed_time])

    # Set ilon array 
    if 'Lon' in axis:
        ilon = np.linspace(-180, 180, 91)
        if axis[0] == 'Lon': i1 = ilon
        if axis[1] == 'Lon': i2 = ilon
        fixed_lon = ''
    else:
        if fixed_lon == '':  print("ERROR, fixed_lon not passed."); return
        ilon = np.array([fixed_lon])

    # Set ilat array 
    if 'Lat' in axis:
        ilat = np.linspace(-90, 90, 73)
        if axis[0] == 'Lat': i1 = ilat
        if axis[1] == 'Lat': i2 = ilat
        fixed_lat = ''
    else:
        if fixed_lat == '':  print("ERROR, fixed_lat not passed."); return
        ilat = np.array([fixed_lat])

    # Set ialt array 
    if 'Alt' in axis:
        # Using quick bisection solve to find alt range.
        dh, h = 50., 50.
        v = interp([12, 0, 0, h])
        for _ in range(100):
            if np.isnan(v):
                if dh < 0.: dh = -.5 * dh
            else:
                if abs(dh) < 1e-6: break
                if dh > 0.: dh = -.5 * dh
            h += dh
            v = interp([12, 0, 0, h])
        h1 = h
        dh = 200.
        for _ in range(100):
            if np.isnan(v):
                if dh > 0.: dh = -.5 * dh
            else:
                if abs(dh) < 1e-6: break
                if dh < 0.: dh = -.5 * dh
            h += dh
            v = interp([12, 0, 0, h])
        ialt = np.linspace(h1, h, 101)
        if axis[0] == 'Alt': i1 = ialt
        if axis[1] == 'Alt': i2 = ialt
        fixed_alt = ''
    else:
        if fixed_alt == '':  print("ERROR, fixed_alt not passed."); return
        ialt = np.array([fixed_alt])

    # Set grid to interpolate to
    i1mg, i2mg = np.meshgrid(np.array(i1), np.array(i2))
    w = np.reshape(i1mg,-1)  # Dummy to store shape for next steps

    if axis[0] == 'Lon':
        x = np.reshape(i1mg,-1)
    elif axis[1] == 'Lon':
        x = np.reshape(i2mg,-1)
    else:
        x = np.full(w.shape, fixed_lon)

    if axis[0] == 'Lat':
        y = np.reshape(i1mg,-1)
    elif axis[1] == 'Lat':
        y = np.reshape(i2mg,-1)
    else:
        y = np.full(w.shape, fixed_lat)
    
    z = np.ones(w.shape)  # Either 1km or 1Re are fine for coordinate transform

    # The t array needs to be timestamps
    if axis[0] == 'Time':
        t = np.reshape(i1mg,-1)
        t = basets + t*3600.
    elif axis[1] == 'Time':
        t = np.reshape(i2mg,-1)
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
        xx,yy,zz,units = ConvertCoord(t,x,y,z,plotCoord,'sph','GEO','sph')
        x, y, z = xx, yy, zz
        xx,yy,zz,units = ConvertCoord(t,x,y,z,'GEO','sph','GDZ','sph')

    # Don't use transformed zz, reset to desired altitude (km) to interpolate
    if axis[0] == 'Alt':
        zz = np.reshape(i1mg,-1)
    elif axis[1] == 'Alt':
        zz = np.reshape(i2mg,-1)
    else:
        zz = np.full(w.shape, fixed_alt)

    # The tt array needs to be in hours
    if axis[0] == 'Time':
        tt = np.reshape(i1mg,-1)
    elif axis[1] == 'Time':
        tt = np.reshape(i2mg,-1)
    else:
        tt = np.full(w.shape, fixed_time)

    positions = np.transpose([tt,xx,yy,zz])
    val = interp(positions)

    if logscale:
        val[val <= 0.] = np.nan
        val = np.log10(val)        
        varname = 'log10<br>'+varname

    val2 = np.reshape(val,i1mg.shape)

    cmin = np.min(val)
    cmax = np.max(val)
    if crange != '':
        cmin = float(crange[0])
        cmax = float(crange[1])

    def plot_2D(v1 = i1, v2 = i2):
        return val2
    plot1 = Kamodo(plot_2D = plot_2D)
    fig = plot1.plot(plot_2D = dict())

    # Set colorscale
    fig = toColor(fig,colorscale=colorscale)

    # Set plot options
    fig.update_traces(
        zmin=cmin, zmax=cmax, colorbar=dict(title=varname, tickformat=".3g"),
        hovertemplate=str(axis[0])+": %{x:.2f}<br>"+
            str(axis[1])+": %{y:.2f}<br>"+
            varname+": %{z:.6g}<br>"+"<extra></extra>"
    )
    fig.update_xaxes(title=axis[0])
    fig.update_yaxes(title=axis[1])
    if axis[0] == 'Lon' or axis[0] == 'Lat':
        fig.update_xaxes(tick0=0.,dtick=45.)
    if axis[0] == 'Time':
        fig.update_xaxes(tick0=0.,dtick=3.)
    if axis[1] == 'Lon' or axis[1] == 'Lat':
        fig.update_yaxes(tick0=0.,dtick=45.)
    if axis[1] == 'Time':
        fig.update_yaxes(tick0=0.,dtick=3.)
    timestr = date.strftime("%Y/%m/%d")
    subt = model+',  '+plotCoord+' Coordinates,  '+timestr+',  Slice at:'
    if fixed_time != "": subt = subt+'  '+str(fixed_time)+' hrs'
    if fixed_alt != "":  subt = subt+'  '+str(fixed_alt) +' km Altitude'
    if fixed_lon != "":  subt = subt+'  '+str(fixed_lon) +' deg Longitude'
    if fixed_lat != "":  subt = subt+'  '+str(fixed_lat) +' deg Latitude'
    fig.update_layout(
        plot_bgcolor="white",scene_aspectmode='data',
        title=dict(text=title+'<br>'+subt, 
                   yref="container", yanchor="top", x=0.01, y=0.97,
                   font=dict(size=16, family="sans serif", color="#000000")),
        margin=dict(l=0),width=750,height=375)

    if plotType == 'Lon-Lat':
        # Shoreline (land/water boundaries)
        if plotCoord == 'GDZ':
            # Convert into GEO as we only use lon/lat
            pos = shoreline(coord='GEO',coordT='sph',utcts=slicets)
        else:
            pos = shoreline(coord=plotCoord,coordT='sph',utcts=slicets)
        fig.add_scattergl(mode='lines', x=pos[:,0], y=pos[:,1], 
                          hoverinfo='skip', showlegend=False,
                          line=dict(width=1,color='#EEEEEE'))

    toc = time.perf_counter()
    print(f"Time: {toc - tic:0.4f} seconds")
    return fig

