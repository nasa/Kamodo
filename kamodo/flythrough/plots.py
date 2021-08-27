# %%
import pandas as pd
import numpy as np
import time
import datetime
from datetime import datetime as dt
from datetime import timezone
from spacepy import coordinates as coord
from spacepy.time import Ticktock
from astropy.constants import R_earth
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.offline import iplot
import kamodo
from flythrough.utils import ConvertCoord

def SatPlot4D(var,time,c1,c2,c3,vard,varu,inCoordName,inCoordType,plotCoord,groupby,model,
              displayplot=True,type='3D',body='black',divfile='',htmlfile=''):
    """New 4D plotting for satellite trajectories using plotly by Darren De Zeeuw
    
    __Required variables__
    
    var: string of variable name
    time: time formatted as a timestamp in UTC
    c1: latitude  or X
    c2: longitude or Y
    c3: altitude  or Z
    vard: data of variable var, same size array as positions
    varu: string of variable var units
    inCoordName: string for incoming coordinate system.  GDZ, GEO, GSM, GSE, SM, GEI, MAG, RLL
    inCoordType: string for incoming coordinate type.  car, sph
    plotCoord: string for coordinate system used in 3D plot. Assumes cartesian type.
    groupby: grouping of data for animation, values include
        all, day, hour, minute, N, orbitE, orbitM
    model: string of name of model the data was extracted from
    
    __Optional variables__
    
    displayplot: logical to show/hide displayed plot (may want false when saving htmlfile)
    type: string for choice of plot type, values: 3D, 1D, 2D, 2DLT
    body: string for choice of 3D inner body, values: black, earth (only GEO), none
    divfile: string with filename to save a html div file of the plot
    htmlfile: string with filename to save a full html file of the plot
    """

    REkm = (R_earth.value/1000.)

    import kamodo.flythrough.model_wrapper as MW
    coord_names=MW.coord_names(inCoordName,inCoordType)
    coord_units=MW.coord_units(inCoordName,inCoordType)
    if inCoordType=='sph':
        for key in coord_names.keys(): coord_names[key]=coord_names[key][:3]

    if type == "3D":
        #Convert incoming coordinates into plot coordinages (cartesian)
        xx,yy,zz,units = ConvertCoord(time,c1,c2,c3,inCoordName,inCoordType,plotCoord,'car')

        # Create dictionary block to pass to plotting with selected options
        plot_dict=dict(
            title = 'Satellite extraction from model: '+model+"<br>"+plotCoord+" coordinates",  # Displayed title for plot, can use <br> for new lines
            sats = ["Sat1"],            # Array of satellites to include in plot
            Sat1 = dict(
                display_name = "",
                time = dict(format='timestamp', data=time),  # possible formats: datetime, timestamp (assumes UTC)
                vars = dict(
                    x = dict(units=units[0], data=xx),
                    y = dict(units=units[1], data=yy),
                    z = dict(units=units[2], data=zz)
                ),
                position_variables = ["x", "y", "z"],     # three variables to use for position
            ),
            options = dict(
                position_units = "R_E", # possible values: R_E, km, ""
                var = var,            # variable to use for colorscale
                hover_vars = [],     # other information for hoverinfo display
                quiver = False,         # logical value to display or hide quivers
                quiver_scale = 0.1,     # length scale of quivers
                quiver_skip = 0,        # points to skip between displaying quivers
                groupby = groupby,        # possible values: all, day, hour, minute, N (integer, show N values at a time)
                                    #   orbitE (break at equator crossing S->N), orbitM (break at prime meridian crossing)
                body = body,         # possible values: black, earth, and any other value is no body
                colorscale = "Viridis", # named colorscale
                REkm = REkm,  # Earth radius in km
                coord = coord,  # Coordinate system of plot
            ),
        )
        # Fixed position variables already included, now add passed in variable to dictionary
        plot_dict['Sat1']['vars'][var]=dict(units=varu, data=vard)
        plot_dict['Sat1']['vars'][coord_names['c1']]=dict(units=coord_units['c1'], data=c1)
        plot_dict['Sat1']['vars'][coord_names['c2']]=dict(units=coord_units['c2'], data=c2)
        plot_dict['Sat1']['vars'][coord_names['c3']]=dict(units=coord_units['c3'], data=c3)
        plot_dict['options']['hover_vars']=[coord_names['c1'],coord_names['c2'],coord_names['c3']]

        # Execute creation and display of figure
        fig=custom3Dsat(plot_dict,vbose=0)
        if divfile != '':
            print('-saving html div file: ',divfile)
            fig.write_html(divfile,full_html=False)
        if htmlfile != '':
            print('-saving full html file: ',htmlfile)
            fig.write_html(htmlfile,full_html=True)
        if displayplot:
            iplot(fig)

    if type == "1D" or type == "2D" or type == "2DLT":
        # Create dictionary block to pass to plotting with selected options
        plot_dict=dict(
            title = 'Satellite extraction from model: '+model,  # Displayed title for plot, can use <br> for new lines
            sats = ["Sat1"],            # Array of satellites to include in plot
            Sat1 = dict(
                display_name = "",
                time = dict(format='timestamp', data=time),  # possible formats: datetime, timestamp (assumes UTC)
                vars = dict(),
                position_variables = [],     # three variables to use for position
            ),
            options = dict(
                position_units = "", # possible values: R_E, km, ""
                var = var,            # variable to use for colorscale
                hover_vars = [],     # other information for hoverinfo display
                groupby = groupby,        # possible values: all, day, hour, minute, N (integer, show N values at a time)
                                    #   orbitE (break at equator crossing S->N), orbitM (break at prime meridian crossing)
            ),
        )
        # Fixed position variables already included, now add passed in variable to dictionary
        plot_dict['Sat1']['vars'][var]=dict(units=varu, data=vard)
        plot_dict['Sat1']['vars'][coord_names['c1']]=dict(units=coord_units['c1'], data=c1)
        plot_dict['Sat1']['vars'][coord_names['c2']]=dict(units=coord_units['c2'], data=c2)
        plot_dict['Sat1']['vars'][coord_names['c3']]=dict(units=coord_units['c3'], data=c3)
        plot_dict['Sat1']['position_variables']=[coord_names['c1'],coord_names['c2'],coord_names['c3']]
        plot_dict['options']['hover_vars']=[coord_names['c1'],coord_names['c2'],coord_names['c3']]

        # Execute creation and display of figure
        if type == "1D":
            fig=custom1Dsat(plot_dict,vbose=0)
        elif type == "2D":
            fig=custom2Dsat(plot_dict,vbose=0)
        elif type == "2DLT":
            fig=custom2Dsat(plot_dict,useLT=True,vbose=0)

        if divfile != '':
            print('-saving html div file: ',divfile)
            fig.write_html(divfile,full_html=False)
        if htmlfile != '':
            print('-saving full html file: ',htmlfile)
            fig.write_html(htmlfile,full_html=True)
        if displayplot:
            iplot(fig)


# ===============================================================================================
# ===============================================================================================
def custom3Dsat(datad, vbose=1):
    '''
    This function creates a custom 3D satellite plot, returning a plotly figure object.
    
    Parameters
    ----------
    datad: This is a data dictionary with the data used to create the plot
    vbose: An optional verbosity value, 0 will only print out warnings and errors. Default is 1.
    
    Returns
    -------
    fig: A plotly figure object that can then be visualized.
    
    Other
    -----
    The python code block below will setup and display a working demo.

import numpy as np
import datetime
from datetime import timezone
from plotly.offline import iplot

# Build a datetime array to use in the dictionary
base = datetime.datetime(2000, 1, 1).replace(tzinfo=timezone.utc)
arr = np.array([base + datetime.timedelta(minutes=30*i) for i in range(8)])

sample=dict(
    title = 'Plot Title Here',  # Displayed title for plot, can use <br> for new lines
    sats = ["Sat1"],            # Array of satellites to include in plot
    Sat1 = dict(
        display_name = "Fake Satellite",
        time = dict(format='datetime', data=arr),  # possible formats: datetime, timestamp (assumes UTC)
        vars = dict(
            x = dict(units='R_E', data=np.array(np.arange(1.,9.))),
            y = dict(units='R_E', data=np.array(np.arange(1.,9.))),
            z = dict(units='R_E', data=np.array(np.arange(1.,9.))),
            p = dict(units='nP', data=np.array(np.arange(11.,19.))),
            U_x = dict(units='km/s', data=np.array(-1.*np.arange(11.,19.))),
            U_y = dict(units='km/s', data=np.array(np.arange(21.,29.))),
            U_z = dict(units='km/s', data=np.array(np.arange(31.,39.))),
        ),
        position_variables = ["x", "y", "z"],     # three variables to use for position
        vector_variables = ["U_x", "U_y", "U_z"], # three variables to use for quiver if quiver is True
    ),
    options = dict(
        position_units = "R_E", # possible values: R_E, km, ""
        var = "p",              # variable to use for colorscale
        hover_vars = ["U_x"],   # other information for hoverinfo display
        quiver = True,          # logical value to display or hide quivers
        quiver_scale = 0.1,     # length scale of quivers
        quiver_skip = 0,        # points to skip between displaying quivers
        groupby = "orbitM",     # possible values: all, day, hour, minute, N (integer, show N values at a time)
                                #   orbitE (break at equator crossing S->N), orbitM (break at prime meridian crossing)
        body = "black",         # possible values: black, earth, and any other value is no body
        colorscale = "Viridis", # named colorscale
        REkm = 6.3781E3,        # Earth radius in km
    ),
)

fig=custom3Dsat(sample)
iplot(fig)
    '''

    # ===============================================================================================  
    # Start timer
    tic = time.perf_counter()

    # Start with error checking ...
    if 'title' not in datad:
        print("Warning, no title given for plot.")
        txttop = "No Title"
    else:
        txttop = datad['title']
    if 'var' not in datad['options']:
        print("ERROR, no variable selected to plot, returning.")
        return None
    var=datad['options']['var']
    if var == "time":
        varu=""
    else:
        varu=datad[datad['sats'][0]]['vars'][var]['units']
    if 'REkm' in datad['options']:
        REkm = datad['options']['REkm']
    else:
        REkm=6.3781E3
    scale=datad['options']['position_units']
    if 'groupby' in datad['options']:
        groupby = datad['options']['groupby']
    else:
        groupby = "all"
    if 'quiver' in datad['options']:
        quiver=datad['options']['quiver']
    else:
        quiver=False
    if quiver:
        quiverscale=datad['options']['quiver_scale']
        if scale == "km":
            quiverscale=quiverscale*REkm
        quiverskip=int(datad['options']['quiver_skip'])
    if 'body' in datad['options']:
        body=datad['options']['body']
    else:
        body = "none"
    if 'colorscale' in datad['options']:
        colorscale = datad['options']['colorscale']
    else:
        colorscale = "Viridis"
    if 'coord' in datad['options']:
        coord = datad['options']['coord']
    else:
        coord = ""

    # set initial values used later, including loop over all sats
    xmin=0.
    xmax=0.
    ymin=0.
    ymax=0.
    zmin=0.
    zmax=0.
    cmin= 1.e99
    cmax=-1.e99
    localts=dict()
    localtimestring=dict()
    agroup=dict()
    ugroup=()
    for sat in datad['sats']:
        sPts=len(datad[sat]['vars'][datad[sat]['position_variables'][0]]['data'])
        # Set localtimestring values
        notime=False
        if 'time' in datad[sat]:
            if datad[sat]['time']['format'] == "datetime":
                localts[sat]=np.array([d.timestamp() for d in datad[sat]['time']['data']])
            elif datad[sat]['time']['format'] == "timestamp":
                localts[sat]=datad[sat]['time']['data'].copy()
            else:
                print("ERROR, Unknown time format.")
                return None
            localtimestring[sat] = np.array([dt.fromtimestamp(int(d),tz=timezone.utc).strftime\
                ("%Y-%m-%d %H:%M:%S") for d in localts[sat]])
        else:
            notime=True
            if var == "time":
                print("ERROR, no time given and plot var selected is time")
                return None
            localtimestring[sat]=np.array(["point "+str(i+1) for i in range(sPts)])

        # Find global contour min/max
        if var == "time":
            c=localts[sat]
        else:
            c=datad[sat]['vars'][var]['data']
        cmin=min(cmin,min(c))
        cmax=max(cmax,max(c))

        # Create array of possible 'groupby' value
        if groupby == "day":
            if notime:
                print("ERROR, no time given and groupby value is",groupby)
                return None
            agroup[sat] = np.array([dt.fromtimestamp(int(d),tz=timezone.utc).strftime\
                ("%Y-%m-%d") for d in localts[sat]])
        elif groupby == "hour":
            if notime:
                print("ERROR, no time given and groupby value is",groupby)
                return None
            agroup[sat] = np.array([dt.fromtimestamp(int(d),tz=timezone.utc).strftime\
                ("%Y-%m-%d %H") for d in localts[sat]])
        elif groupby == "minute":
            if notime:
                print("ERROR, no time given and groupby value is",groupby)
                return None
            agroup[sat] = np.array([dt.fromtimestamp(int(d),tz=timezone.utc).strftime\
                ("%Y-%m-%d %H:%M") for d in localts[sat]])
        elif groupby == "orbitM":
            # Satellite path crosses prime meridian
            x=datad[sat]['vars'][datad[sat]['position_variables'][0]]['data']
            y=datad[sat]['vars'][datad[sat]['position_variables'][1]]['data']
            bgroup = ['orbit'] * len(x)
            j=1
            for i in range(sPts):
                if i != 0:
                    if x[i] > 0. and (y[i]*y[i-1]) < 0.:
                        j+=1
                bgroup[i] = "orbit "+str(j)
            agroup[sat]=np.array(bgroup)
        elif groupby == "orbitE":
            # Satellite path crosses equator going North
            z=datad[sat]['vars'][datad[sat]['position_variables'][2]]['data']
            bgroup = ['orbit'] * len(z)
            j=1
            for i in range(sPts):
                if i != 0:
                    if (z[i]>0. and z[i-1]<0.):
                        j+=1
                bgroup[i] = "orbit "+str(j)
            agroup[sat]=np.array(bgroup)
        elif groupby.isdigit():
            gb=int(groupby)
            agroup[sat] = np.array(["points "+str(int(i/gb)*gb+1)+" - "+str(int(i/gb)*gb+gb) for i in range(sPts)])
        else:
            agroup[sat] = np.array(["all" for i in range(sPts)])

        # Use pandas unique function rather than numpy. Its faster and does not sort the results.
        ugroup=pd.unique(np.append(ugroup, pd.unique(agroup[sat])))

    ngroup = len(ugroup)

    # Build DUMMY data block to insert as needed.
    data_dict_dummy = {
        "type": "scatter3d",
        "name": "dummy", "x": [0.], "y": [0.], "z": [0.],
        "mode": "lines", "line": {"width": 1},
        "hoverinfo": "none",
    }

    # ===============================================================================================  AAA
    # make figure dictionary pieces
    fig_dict = {"data": [], "layout": {}, "frames": []}
    fig_data_saved = {"data": []}
    sliders_dict = {
        "active": 0,
        "yanchor": "top", "xanchor": "left",
        "currentvalue": {
            "prefix": "Currently showing: ",
            "visible": True,
            "xanchor": "left"
        },
        "transition": {"duration": 0},
        "pad": {"b": 10, "t": 10},
        "len": 0.9,
        "x": 0.1,
        "y": 0,
        "steps": []
    }

    # Actual plot creation loop
    for date in ugroup:
        frame = {"data": [], "name": date}
        for sat in datad['sats']:
            x=datad[sat]['vars'][datad[sat]['position_variables'][0]]['data']
            y=datad[sat]['vars'][datad[sat]['position_variables'][1]]['data']
            z=datad[sat]['vars'][datad[sat]['position_variables'][2]]['data']
            sc=1.
            if scale == "km" and datad[sat]['vars'][datad[sat]['position_variables'][0]]['units'] == "R_E":
                sc=REkm
            elif scale == "R_E" and datad[sat]['vars'][datad[sat]['position_variables'][0]]['units'] == "km":
                sc=1./REkm
            if var == "time":
                c=localts[sat]
                varline=""
            else:
                c=datad[sat]['vars'][var]['data']
                varline=var+": %{marker.color:.4g} "+varu+"<br>"
            if quiver:
                qxvar=datad[sat]['vector_variables'][0]
                qyvar=datad[sat]['vector_variables'][1]
                qzvar=datad[sat]['vector_variables'][2]
                qx=datad[sat]['vars'][qxvar]['data']
                qy=datad[sat]['vars'][qyvar]['data']
                qz=datad[sat]['vars'][qzvar]['data']

            # Update position min/max values
            if date == ugroup[0]:
                xmin=min(xmin,min(x*sc))
                xmax=max(xmax,max(x*sc))
                ymin=min(ymin,min(y*sc))
                ymax=max(ymax,max(y*sc))
                zmin=min(zmin,min(z*sc))
                zmax=max(zmax,max(z*sc))

            # Compute mask to restrict all data in trace
            mask = date == agroup[sat]
            
            # Create hover information, including extras passed in. Quiver shows additional variables.
            Nhv = len(datad['options']['hover_vars'])
            cd=[]
            cd.append(localtimestring[sat][mask])
            qline=""
            Ndv=1
            if quiver:
                cd.append(qx[mask])
                cd.append(qy[mask])
                cd.append(qz[mask])
                qline+=qxvar+": %{customdata[1]:.2f}<br>"
                qline+=qyvar+": %{customdata[2]:.2f}<br>"
                qline+=qzvar+": %{customdata[3]:.2f}<br>"
                Ndv+=3
            for i in range(Nhv):
                cd.append(datad[sat]['vars'][datad['options']['hover_vars'][i]]['data'][mask])
                qline+=datad['options']['hover_vars'][i]+": %{customdata["+str(Ndv)+"]:.2f} "+\
                    datad[sat]['vars'][datad['options']['hover_vars'][i]]['units']+"<br>"
                Ndv+=1
            cd=np.asarray(cd).T
            dateline="%{customdata[0]}<br>"
            
            # Build data block with mask
            data_dict = {
                "type": "scatter3d",
                "name": date,
                "x": list(x[mask]*sc), "y": list(y[mask]*sc), "z": list(z[mask]*sc),
                "mode": "markers+lines",
                "marker": {
                    "size": 4, "cmin": cmin, "cmax": cmax, "color": list(c[mask]),
                    "showscale": True, "colorscale": colorscale,
                    "colorbar": { "title": "<b>"+var+"</b><br>["+varu+"]", "tickformat": ".3g" }
                },
                "line": {"width": 3, "color": "rgba(22,22,22,0.2)"},
                "customdata": cd,
                "hovertemplate": "<b>"+datad[sat]['display_name']+"</b>"+
                    "<br>X: %{x:.4f} "+scale+"<br>Y: %{y:.4f} "+scale+"<br>Z: %{z:.4f} "+scale+"<br>"+
                    qline+varline+dateline+"<extra></extra>",
            }

            # If time is colorbar variable, hide labels by selecting ticks out of range
            if var == "time":
                data_dict.marker.colorbar['tickvals']=(0,1)

            # Put each part of sequence in frame data block
            frame["data"].append(data_dict)

            # First in sequence, put dummy in main data block
            if date == ugroup[0]:
                fig_dict["data"].append(data_dict_dummy)
            
            # Start quiver
            if quiver:
                # Compute values to put in quiver trace
                # Make array max possible size (triple len(x)), fill, and trim as needed
                xx=np.concatenate([x,x,x])
                yy=np.concatenate([y,y,y])
                zz=np.concatenate([z,z,z])
                qxx=qx
                qyy=qy
                qzz=qz
                # Build new position array, element by element
                j=0
                for i in range(len(mask)):
                    if mask[i]:
                        if i%(quiverskip+1) == 0:
                            xx[j]=x[i]*sc
                            yy[j]=y[i]*sc
                            zz[j]=z[i]*sc
                            xx[j+1]=x[i]*sc+quiverscale*qx[i]
                            yy[j+1]=y[i]*sc+quiverscale*qy[i]
                            zz[j+1]=z[i]*sc+quiverscale*qz[i]
                            xx[j+2]=None
                            yy[j+2]=None
                            zz[j+2]=None
                            j+=3
                xx=np.array(xx[0:j], dtype=np.float64)
                yy=np.array(yy[0:j], dtype=np.float64)
                zz=np.array(zz[0:j], dtype=np.float64)

                # Update position min/max values
                xmin=min(xmin,min(xx))
                xmax=max(xmax,max(xx))
                ymin=min(ymin,min(yy))
                ymax=max(ymax,max(yy))
                zmin=min(zmin,min(zz))
                zmax=max(zmax,max(zz))
                
                # Build data block
                data_dict = {
                    "type": "scatter3d",
                    "name": "positions", "x": list(xx), "y": list(yy), "z": list(zz),
                    "mode": "lines", "line": {"width": 3, "color": "rgba(22,22,22,0.2)"},
                    "hoverinfo": "none",
                }
                # Put each part of sequence in frame data block
                frame["data"].append(data_dict)
                
                # First in sequence, put in main data block
                if date == ugroup[0]:
                    fig_dict["data"].append(data_dict_dummy)
                
        fig_dict["frames"].append(frame)
        slider_step = {"args": [
            [date],
            {"frame": {"duration": 300, "redraw": True},
             "mode": "immediate",
             "transition": {"duration": 0}}
        ],
            "label": date,
            "method": "animate"}
        sliders_dict["steps"].append(slider_step)
        
    # Assemble frame and slider pieces
    fig_dict["layout"]["sliders"] = [sliders_dict]

    # ===============================================================================================  BBB
    if ngroup > 1:
        for sat in datad['sats']:
            # Add trace if more than one group.
            # This shows the whole trajectory when a subsection of data is showing.
            x=datad[sat]['vars'][datad[sat]['position_variables'][0]]['data']
            y=datad[sat]['vars'][datad[sat]['position_variables'][1]]['data']
            z=datad[sat]['vars'][datad[sat]['position_variables'][2]]['data']
            sc=1.
            if scale == "km" and datad[sat]['vars'][datad[sat]['position_variables'][0]]['units'] == "R_E":
                sc=REkm
            elif scale == "R_E" and datad[sat]['vars'][datad[sat]['position_variables'][0]]['units'] == "km":
                sc=1./REkm
            
            # Build data block
            data_dict = {
                "type": "scatter3d",
                "name": "positions", "x": list(x*sc), "y": list(y*sc), "z": list(z*sc),
                "mode": "lines", "line": {"width": 3, "color": "rgba(22,22,22,0.2)"},
                "hoverinfo": "none",
            }

            # Put into main data block
            fig_dict["data"].append(data_dict)

    # ===============================================================================================  CCC
    ticE = time.perf_counter()
    # Load points and add 1 RE sphere, padded to cover all data positions
    if body == "black":
        dataXYZ = pd.read_csv('https://ccmc.gsfc.nasa.gov/Kamodo/demo/sphereXYZ.csv')
        dataIJK = pd.read_csv('https://ccmc.gsfc.nasa.gov/Kamodo/demo/sphereIJK.csv')
        if scale == "km":
            dataXYZ *= REkm
        # Build data block
        data_dict = {
            "type": "mesh3d",
            "name": '1 R_E sphere',
            "x": list(np.append(dataXYZ['x'],(xmin,xmax))),
            "y": list(np.append(dataXYZ['y'],(ymin,ymax))),
            "z": list(np.append(dataXYZ['z'],(zmin,zmax))),
            "i": list(dataIJK['i']),
            "j": list(dataIJK['j']),
            "k": list(dataIJK['k']),
            "facecolor": list(dataIJK['c']),
            "flatshading": True,
            "hovertemplate": "Earth<extra></extra>",
        }
        # Put in main data block
        fig_dict["data"].append(data_dict)
    elif body == "earth" and coord == "GEO":
        dataXYZ = pd.read_csv('https://ccmc.gsfc.nasa.gov/ungrouped/GM_IM/EarthXYZ.csv')
        dataIJK = pd.read_csv('https://ccmc.gsfc.nasa.gov/ungrouped/GM_IM/EarthIJKRGB.csv')
        if scale == "km":
            dataXYZ *= REkm
        color=np.array(["rgb("+str(dataIJK['r'][i])+","+str(dataIJK['g'][i])+","+str(dataIJK['b'][i])+")" \
                        for i in range(len(dataIJK['r']))])
        # Need to reverse x,y from file to be in proper GEO coords (180 degree rotation)
        xe=-dataXYZ['x']
        ye=-dataXYZ['y']
        ze= dataXYZ['z']
        # Build data block
        data_dict = {
            "type": "mesh3d",
            "name": '1 R_E sphere',
            "x": np.append(xe,(xmin,xmax)),
            "y": np.append(ye,(ymin,ymax)),
            "z": np.append(ze,(zmin,zmax)),
            "i": dataIJK['i'],
            "j": dataIJK['j'],
            "k": dataIJK['k'],
            "facecolor": color,
            "hovertemplate": "Earth<extra></extra>",
        }
        # Put in main data block
        fig_dict["data"].append(data_dict)
    tocE = time.perf_counter()
    if vbose > 0:
        print(f"  -time loading Earth: {tocE - ticE:0.4f} seconds")

    # ===============================================================================================  DDD
    # Set layout values
    fig_dict["layout"]["height"] = 700
    fig_dict["layout"]["width"] = 800
    fig_dict["layout"]["scene_aspectmode"] = "data"
    fig_dict["layout"]["scene"] = dict(xaxis=dict(title=dict(text="X ["+scale+"]")),
                                       yaxis=dict(title=dict(text="Y ["+scale+"]")),
                                       zaxis=dict(title=dict(text="Z ["+scale+"]")))
    fig_dict["layout"]["title_text"] = txttop
    fig_dict["layout"]["showlegend"] = False
    fig_dict["layout"]["scene_camera"] = dict(center=dict(x=0, y=0, z=0))
    fig_dict["layout"]["hoverlabel_align"] = 'right'
    if ngroup > 1:
        fig_dict["layout"]["updatemenus"] = [
            {
                "buttons": [
                    {
                        "args": [None, {"frame": {"duration": 500, "redraw": True},
                                        "fromcurrent": True, "transition": {"duration": 0}}],
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
                "pad": {"r": 10, "t": 35},
                "showactive": False,
                "type": "buttons",
                "x": 0.1,
                "xanchor": "right",
                "y": 0,
                "yanchor": "top"
            }
        ]


    # end timer                                                                               
    toc = time.perf_counter()
    if vbose > 0:
        print(f"Total time creating figure object: {toc - tic:0.4f} seconds")

    fig3 = go.Figure(fig_dict)

    return fig3

# ===============================================================================================
# ===============================================================================================
def custom1Dsat(datad, vbose=1):
    """
    This function creates a custom 1D satellite plot, returning a plotly figure object.
   
    Parameters
    ----------
    datad: This is a data dictionary with the data used to create the plot
    vbose: Optional verbosity value, 0 will only print out warnings and errors. Default is 1.
    
    Returns
    -------
    fig: A plotly figure object that can then be visualized.
    
    """

    # Start 1D fig
    var=datad['options']['var']
    varu=datad[datad['sats'][0]]['vars'][var]['units']
    Nhv = len(datad['options']['hover_vars'])
    
    localts=dict()
    localtimestring=dict()
    localdt=dict()
    txttop = datad['title']

    # For now this only makes the plot for the last sat in the dictionary.
    for sat in datad['sats']:
        sPts=len(datad[sat]['time']['data'])
        # Set localtimestring values
        notime=False
        if 'time' in datad[sat]:
            if datad[sat]['time']['format'] == "datetime":
                localts[sat]=np.array([d.timestamp() for d in datad[sat]['time']['data']])
            elif datad[sat]['time']['format'] == "timestamp":
                localts[sat]=datad[sat]['time']['data'].copy()
            else:
                print("ERROR, Unknown time format.")
                return None
            localtimestring[sat] = np.array([datetime.datetime.fromtimestamp(int(d),tz=timezone.utc).strftime\
                ("%Y-%m-%d %H:%M:%S") for d in localts[sat]])
            localdt[sat] = np.array([datetime.datetime.fromtimestamp(int(d),tz=timezone.utc) for d in localts[sat]])
        else:
            notime=True
            if var == "time":
                print("ERROR, no time given and plot var selected is time")
                return None
            localtimestring[sat]=np.array(["point "+str(i+1) for i in range(sPts)])

        # Find global contour min/max
        if var == "time":
            c=localts[sat]
        else:
            c=datad[sat]['vars'][var]['data']

        fig1 = make_subplots(rows=(Nhv+1), cols=1, shared_xaxes=True, vertical_spacing=0.04)

        fig1.add_trace(go.Scatter(x=localdt[sat], y=c, name=var,
                                  mode='lines', line= dict(shape='linear', color='black'),
                                  hovertemplate=var+': %{y:.4g}<br>%{x}<extra></extra>',
                                ),
                       row=1, col=1)
        fig1.update_yaxes(title_text='<b>'+var+'</b><br>['+varu+']', exponentformat='e', row=1, col=1)
        fig1.update_layout(yaxis=dict(title=dict(font=dict(size=12))))

        for i in range(Nhv):
            tmpv=datad['options']['hover_vars'][i]
            tmpu=datad[sat]['vars'][tmpv]['units']
            fig1.add_trace(go.Scatter(x=localdt[sat], y=datad[sat]['vars'][tmpv]['data'], name=tmpv,
                                      mode='lines', line= dict(shape='linear', color='black'),
                                      hovertemplate=tmpv+': %{y:.4g}<br>%{x}<extra></extra>',
                                     ),
                           row=(i+2), col=1)
#            if tmpv == "Alt":
#                tmpu=" [km]"
            if tmpv == "Lon":
#                tmpu=" [deg]"
                fig1.update_yaxes(tick0=0., dtick=90., row=(i+2), col=1)
            if tmpv == "Lat":
#                tmpu=" [deg]"
                fig1.update_yaxes(tick0=0., dtick=30., row=(i+2), col=1)
            ya='yaxis'+str(i+2)
            ys="dict(text='<b>"+tmpv+"</b> ["+tmpu+"]',font=dict(size=12))"
            fig1['layout'][ya]['title']=eval(ys)

        fig1.update_layout(height=600, width=800, title_text=txttop, showlegend = False,)

    return fig1

# ===============================================================================================
# ===============================================================================================
def custom2Dsat(datad, useLT=False, vbose=1):
    '''
    This function creates a custom 2D satellite plot in lat/lon or lat/LT, returning a plotly figure object.
    
    Parameters
    ----------
    datad: This is a data dictionary with the data used to create the plot
    useLT: Optional logical to modify plot to use local time instead of longitude on the X axis.
    vbose: Optional verbosity value, 0 will only print out warnings and errors. Default is 1.
    
    Returns
    -------
    fig: A plotly figure object that can then be visualized.
    
    '''

    # ===============================================================================================  
    # Start timer
    tic = time.perf_counter()

    # Start with error checking ...
    if 'title' not in datad:
        print("Warning, no title given for plot.")
        txttop = "No Title"
    else:
        txttop = datad['title']
    if 'var' not in datad['options']:
        print("ERROR, no variable selected to plot, returning.")
        return None
    var=datad['options']['var']
    if var == "time":
        varu=""
    else:
        varu=datad[datad['sats'][0]]['vars'][var]['units']
    if 'REkm' in datad['options']:
        REkm = datad['options']['REkm']
    else:
        REkm=6.3781E3
    scale=datad['options']['position_units']
    if 'groupby' in datad['options']:
        groupby = datad['options']['groupby']
    else:
        groupby = "all"
    if 'body' in datad['options']:
        body=datad['options']['body']
    else:
        body = "none"
    if 'colorscale' in datad['options']:
        colorscale = datad['options']['colorscale']
    else:
        colorscale = "Viridis"

    # set initial values used later, including loop over all sats
    cmin= 1.e99
    cmax=-1.e99
    localts=dict()
    localtimestring=dict()
    agroup=dict()
    ugroup=()
    for sat in datad['sats']:
        sPts=len(datad[sat]['vars'][datad[sat]['position_variables'][0]]['data'])
        # Set localtimestring values
        notime=False
        if 'time' in datad[sat]:
            if datad[sat]['time']['format'] == "datetime":
                localts[sat]=np.array([d.timestamp() for d in datad[sat]['time']['data']])
            elif datad[sat]['time']['format'] == "timestamp":
                localts[sat]=datad[sat]['time']['data'].copy()
            else:
                print("ERROR, Unknown time format.")
                return None
            localtimestring[sat] = np.array([dt.fromtimestamp(int(d),tz=timezone.utc).strftime\
                ("%Y-%m-%d %H:%M:%S") for d in localts[sat]])
        else:
            notime=True
            if var == "time" or useLT:
                print("ERROR, no time given and plot var selected is time")
                return None
            localtimestring[sat]=np.array(["point "+str(i+1) for i in range(sPts)])

        # Find global contour min/max
        if var == "time":
            c=localts[sat]
        else:
            c=datad[sat]['vars'][var]['data']
        cmin=min(cmin,min(c))
        cmax=max(cmax,max(c))

        # Create array of possible 'groupby' value
        if groupby == "day":
            if notime:
                print("ERROR, no time given and groupby value is",groupby)
                return None
            agroup[sat] = np.array([dt.fromtimestamp(int(d),tz=timezone.utc).strftime\
                ("%Y-%m-%d") for d in localts[sat]])
        elif groupby == "hour":
            if notime:
                print("ERROR, no time given and groupby value is",groupby)
                return None
            agroup[sat] = np.array([dt.fromtimestamp(int(d),tz=timezone.utc).strftime\
                ("%Y-%m-%d %H") for d in localts[sat]])
        elif groupby == "minute":
            if notime:
                print("ERROR, no time given and groupby value is",groupby)
                return None
            agroup[sat] = np.array([dt.fromtimestamp(int(d),tz=timezone.utc).strftime\
                ("%Y-%m-%d %H:%M") for d in localts[sat]])
        elif groupby == "orbitM":
            # Satellite path crosses prime meridian
            lon=datad[sat]['vars'][datad[sat]['position_variables'][0]]['data']
            bgroup = ['orbit'] * len(lon)
            j=1
            for i in range(sPts):
                if i != 0:
                    if abs(lon[i]-lon[i-1]) > 180.:
                        j+=1
                bgroup[i] = "orbit "+str(j)
            agroup[sat]=np.array(bgroup)
        elif groupby == "orbitE":
            # Satellite path crosses equator going North
            lat=datad[sat]['vars'][datad[sat]['position_variables'][1]]['data']
            bgroup = ['orbit'] * len(lat)
            j=1
            for i in range(sPts):
                if i != 0:
                    if (lat[i]>0. and lat[i-1]<0.):
                        j+=1
                bgroup[i] = "orbit "+str(j)
            agroup[sat]=np.array(bgroup)
        elif groupby.isdigit():
            gb=int(groupby)
            agroup[sat] = np.array(["points "+str(int(i/gb)*gb+1)+" - "+str(int(i/gb)*gb+gb) for i in range(sPts)])
        else:
            agroup[sat] = np.array(["all" for i in range(sPts)])

        # Use pandas unique function rather than numpy. Its faster and does not sort the results.
        ugroup=pd.unique(np.append(ugroup, pd.unique(agroup[sat])))

    ngroup = len(ugroup)

    # Build DUMMY data block to insert as needed.
    data_dict_dummy = {
        "type": "scatter",
        "name": "dummy", "x": [0.], "y": [0.],
        "mode": "lines", "line": {"width": 1},
        "hoverinfo": "none",
    }

    # ===============================================================================================  AAA
    # make figure dictionary pieces
    fig_dict = {"data": [], "layout": {}, "frames": []}
    fig_data_saved = {"data": []}
    sliders_dict = {
        "active": 0,
        "yanchor": "top", "xanchor": "left",
        "currentvalue": {
            "prefix": "Currently showing: ",
            "visible": True,
            "xanchor": "left"
        },
        "transition": {"duration": 0},
        "pad": {"b": 10, "t": 40},
        "len": 0.9,
        "x": 0.1,
        "y": 0,
        "steps": []
    }

    # Actual plot creation loop
    for date in ugroup:
        frame = {"data": [], "name": date}
        for sat in datad['sats']:
            lon=datad[sat]['vars'][datad[sat]['position_variables'][0]]['data']
            lat=datad[sat]['vars'][datad[sat]['position_variables'][1]]['data']

            h = np.array([dt.fromtimestamp(int(d),tz=timezone.utc).strftime("%H") for d in localts[sat]])
            m = np.array([dt.fromtimestamp(int(d),tz=timezone.utc).strftime("%M") for d in localts[sat]])
            s = np.array([dt.fromtimestamp(int(d),tz=timezone.utc).strftime("%S.%f") for d in localts[sat]])
            lt=np.array(h, dtype=np.float64)
            for i in range(len(h)):
                lt[i]=(float(h[i]))+(float(m[i])/60.)+(float(s[i])/3600.)+lon[i]*24./360.
            lt[lt > 24.]-=24.

            if var == "time":
                c=localts[sat]
                varline=""
            else:
                c=datad[sat]['vars'][var]['data']
                varline=var+": %{marker.color:.4g} "+varu+"<br>"

            # Compute mask to restrict all data in trace
            mask = date == agroup[sat]
            
            # Create hover information, including extras passed in.
            Nhv = len(datad['options']['hover_vars'])
            cd=[]
            cd.append(localtimestring[sat][mask])
            qline=""
            Ndv=1
            for i in range(Nhv):
                hovv=datad['options']['hover_vars'][i]
                if hovv != "Lat" and hovv != "Lon":
                    cd.append(datad[sat]['vars'][hovv]['data'][mask])
                    qline+=hovv+": %{customdata["+str(Ndv)+"]:.2f} "+\
                        datad[sat]['vars'][hovv]['units']+"<br>"
                    Ndv+=1
            if useLT:
                hovv='Lon'
                cd.append(datad[sat]['vars'][hovv]['data'][mask])
                qline+=hovv+": %{customdata["+str(Ndv)+"]:.2f} "+\
                    datad[sat]['vars'][hovv]['units']+"<br>"
                Ndv+=1
            else:
                cd.append(lt[mask])
                qline+="LT: %{customdata["+str(Ndv)+"]:.2f}<br>"
                Ndv+=1
            cd=np.asarray(cd).T
            dateline="%{customdata[0]}<br>"
            
            # Build data block with mask
            if useLT:
                data_dict = {
                    "type": "scatter",
                    "name": date,
                    "x": list(lt[mask]), "y": list(lat[mask]),
                    "mode": "markers",
                    "marker": {
                        "size": 9, "cmin": cmin, "cmax": cmax, "color": list(c[mask]),
                        "showscale": True, "colorscale": colorscale,
                        "colorbar": { "title": "<b>"+var+"</b><br>["+varu+"]", "tickformat": ".3g" }
                    },
                    "customdata": cd,
                    "hovertemplate": "<b>"+datad[sat]['display_name']+"</b>"+
                        "<br>LT: %{x:.4f}<br>Lat: %{y:.4f}<br>"+
                        qline+varline+dateline+"<extra></extra>",
                }
            else:
                data_dict = {
                    "type": "scatter",
                    "name": date,
                    "x": list(lon[mask]), "y": list(lat[mask]),
                    "mode": "markers",
                    "marker": {
                        "size": 9, "cmin": cmin, "cmax": cmax, "color": list(c[mask]),
                        "showscale": True, "colorscale": colorscale,
                        "colorbar": { "title": "<b>"+var+"</b><br>["+varu+"]", "tickformat": ".3g" }
                    },
                    "customdata": cd,
                    "hovertemplate": "<b>"+datad[sat]['display_name']+"</b>"+
                        "<br>Lon: %{x:.4f}<br>Lat: %{y:.4f}<br>"+
                        qline+varline+dateline+"<extra></extra>",
                }

            # If time is colorbar variable, hide labels by selecting ticks out of range
            if var == "time":
                data_dict.marker.colorbar['tickvals']=(0,1)

            if ngroup == 1:
                # Put full sequence in main data block
                fig_dict["data"].append(data_dict)
                
            else:
                # Put each part of sequence in frame data block
                frame["data"].append(data_dict)
                # First in sequence, put dummy in main data block
                if date == ugroup[0]:
                    fig_dict["data"].append(data_dict_dummy)
                
        if ngroup > 1:
            fig_dict["frames"].append(frame)
        slider_step = {"args": [
            [date],
            {"frame": {"duration": 300, "redraw": True},
             "mode": "immediate",
             "transition": {"duration": 0}}
        ],
            "label": date,
            "method": "animate"}
        sliders_dict["steps"].append(slider_step)
        
    if ngroup > 1:
        # Assemble frame and slider pieces
        fig_dict["layout"]["sliders"] = [sliders_dict]

    # ===============================================================================================  BBB
    if ngroup > 1:
        for sat in datad['sats']:
            # Add trace if more than one group.
            # This shows the whole trajectory when a subsection of data is showing.
            lon=datad[sat]['vars'][datad[sat]['position_variables'][0]]['data']
            lat=datad[sat]['vars'][datad[sat]['position_variables'][1]]['data']
            
            # Make array larger than needed (double len(lat)), fill, and trim as needed
            xx=np.concatenate([lon,lon])
            yy=np.concatenate([lat,lat])
            # Build new position array, element by element
            j=0
            for i in range(len(lat)):
                xx[j]=lon[i]
                yy[j]=lat[i]
                j+=1
                if i < (len(lat)-1) and abs(lon[i]-lon[i+1]) > 180.:
                    L=min(abs(lon[i]),abs(360-lon[i]))
                    R=min(abs(lon[i+1]),abs(360-lon[i+1]))
                    Llat=lat[i]+(lat[i+1]-lat[i])*L/(L+R)
                    Rlat=lat[i]+(lat[i+1]-lat[i])*R/(L+R)
                    xx[j]=360.
                    yy[j]=Llat
                    xx[j+1]=None
                    yy[j+1]=None
                    xx[j+2]=0.
                    yy[j+2]=Rlat
                    j+=3
            xx=np.array(xx[0:j], dtype=np.float64)
            yy=np.array(yy[0:j], dtype=np.float64)
            
            # Build data block
            data_dict = {
                "type": "scatter",
                "name": "positions", "x": list(xx), "y": list(yy),
                "mode": "lines", "line": {"width": 1, "color": "rgba(22,22,22,0.2)"},
                "hoverinfo": "none",
            }

            # Put into main data block
            if not useLT:
                fig_dict["data"].append(data_dict)

    # ===============================================================================================  DDD
    # Set layout values
    if useLT:
        fig_dict["layout"]["xaxis"] = {'dtick': 3.0, 'range': [0.0, 24.0], 'tick0': 0.0,
                                       'title': {'text': '<b>LT</b> [hrs]'}}
    else:
        fig_dict["layout"]["xaxis"] = {'dtick': 30.0, 'range': [0.0, 360.0], 'tick0': 0.0,
                                       'title': {'text': '<b>Lon</b> [deg]'}}
    fig_dict["layout"]["yaxis"] = {'dtick': 30.0, 'range': [-90.0, 90.0], 'tick0': 0.0,
                                   'title': {'text': '<b>Lat</b> [deg]'}}
    fig_dict["layout"]["title_text"] = txttop
    fig_dict["layout"]["showlegend"] = False
    fig_dict["layout"]["hoverlabel_align"] = 'right'
    if ngroup > 1:
        fig_dict["layout"]["updatemenus"] = [
            {
                "buttons": [
                    {
                        "args": [None, {"frame": {"duration": 500, "redraw": True},
                                        "fromcurrent": False, "transition": {"duration": 0}}],
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
                "pad": {"r": 10, "t": 65},
                "showactive": False,
                "type": "buttons",
                "x": 0.1,
                "xanchor": "right",
                "y": 0,
                "yanchor": "top"
            }
        ]

    # end timer                                                                               
    toc = time.perf_counter()
    if vbose > 0:
        print(f"Total time creating figure object: {toc - tic:0.4f} seconds")

    fig2 = go.Figure(fig_dict)
    return fig2

# %%
