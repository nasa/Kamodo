# -*- coding: utf-8 -*-
"""
Created on Wed May 12 19:01:01 2021

@author: rringuet

move generalizable functions here for 4D readers
"""
import numpy as np
from kamodo import Kamodo
from datetime import datetime, timezone

def get_cut(data_range, cut_name, new_type, data_unit):
    '''Give text of cut_type and value_type for appropriate message'''
    
    if new_type=='min':
        cut = data_range['min']
        print(f'{cut_name} outside of data range ({data_range["min"]} {data_unit}, '+\
              f'{data_range["max"]} {data_unit}). Using minimum: '+\
                  f'{cut} {data_unit}.')
    elif new_type=='average':  #Don't try to use for datetime strings or similar
        cut=(data_range["min"]+data_range["max"])/2.
        print(f'{cut_name} outside of data range ({data_range["min"]} {data_unit}, '+\
              f'{data_range["max"]} {data_unit}). Using average: '+\
                  f'{cut} {data_unit}.')     
    elif new_type=='max':
        cut = data_range['max']
        print(f'{cut_name} outside of data range ({data_range["min"]} {data_unit}, '+\
                      f'{data_range["max"]} {data_unit}). Using maximum: '+\
                          f'{cut} {data_unit}.')        
    return cut

def initialize_4D_plot(kamodo_object):
    '''Initialize plot related variables'''
    
    #set gridSize for each coordinate
    if len(kamodo_object._time)>1:
        kamodo_object.dt = min(np.diff(kamodo_object._time))
    else: 
        kamodo_object.dt = 0.
    kamodo_object.dlat = min(np.diff(kamodo_object._lat)) #same for _lat_height and _lat_neutral
    kamodo_object.dlon = min(np.diff(kamodo_object._lon)) #same for _lon_height and _lon_neutral
    if len(kamodo_object._height)>1:
        kamodo_object.dheight = min(np.diff(kamodo_object._height))
    else:
        kamodo_object.dheight = 0.
    if hasattr(kamodo_object,'_ilev'):
        kamodo_object.dilev = min(np.diff(kamodo_object._ilev)) #same for _ilev_neutral
    if hasattr(kamodo_object,'_lev'):
        kamodo_object.dlev = min(np.diff(kamodo_object._lev)) #same for _ilev_neutral
    if hasattr(kamodo_object,'_imlev'):
        kamodo_object.dimlev = min(np.diff(kamodo_object._imlev)) #same for _ilev_neutral
    if hasattr(kamodo_object,'_elon'):
        kamodo_object.delon = min(np.diff(kamodo_object._elon))
    if hasattr(kamodo_object,'_elat'):    
        kamodo_object.delat = min(np.diff(kamodo_object._elat))
    if hasattr(kamodo_object,'_mlon'):
        kamodo_object.dmlon = min(np.diff(kamodo_object._mlon))
    if hasattr(kamodo_object,'_mlat'):    
        kamodo_object.dmlat = min(np.diff(kamodo_object._mlat))
    #Don't set net gridsize here since it is determined by the different coordinate dependencies
    
    #set default domains of plotting coordinate variables
    kamodo_object.lonrange0={'min':0.,'max':360.,'n':len(kamodo_object._lon)}
    kamodo_object.latrange0={'min':-90.,'max':90.,'n':len(kamodo_object._lat)}
    kamodo_object.htrange0={'min':kamodo_object._height.min(),'max':kamodo_object._height.max(),
                      'n':len(kamodo_object._height)}  #height in km  
    if hasattr(kamodo_object,'_ilev'):
        kamodo_object.iprange0={'min':kamodo_object._ilev.min(),'max':kamodo_object._ilev.max(),
                      'n':len(kamodo_object._ilev)}
    if hasattr(kamodo_object,'_lev'):
        kamodo_object.prange0={'min':kamodo_object._lev.min(),'max':kamodo_object._lev.max(),
                      'n':len(kamodo_object._lev)}
    if hasattr(kamodo_object,'_imlev'):
        kamodo_object.imprange0={'min':kamodo_object._imlev.min(),'max':kamodo_object._imlev.max(),
                      'n':len(kamodo_object._imlev)}
    kamodo_object.lonrange, kamodo_object.latrange, kamodo_object.htrange = kamodo_object.lonrange0, \
        kamodo_object.latrange0, kamodo_object.htrange0  #keep originals separate
    
    # Initialize values for plotting
    kamodo_object.plots = dict()
    kamodo_object.plottype = "" # empty default to force reinterpolation and regeneration
    kamodo_object.cutVtext, kamodo_object.cutV = 'Height [km]', 400.   #typical height cut
    kamodo_object.cutLtext, kamodo_object.cutL = 'Lat [deg]', 0.
    kamodo_object.nt, kamodo_object.nX, kamodo_object.nY, kamodo_object.nZ = 1, 1, 37, 73   #t, height, lat, lon
    kamodo_object.newx = kamodo_object.cutV
    kamodo_object.newy = np.linspace(-90., 90., kamodo_object.nY)
    kamodo_object.newz = np.linspace(0., 360., kamodo_object.nX)
    kamodo_object.newt = kamodo_object._time[0]
    kamodo_object.tunit, kamodo_object.xunit ='s', 'km'
    kamodo_object.yunit, kamodo_object.zunit = 'deg', 'deg'
    kamodo_object.tol = 1.1
    kamodo_object.nDim = 4
    kamodo_object.plots[kamodo_object.plottype] = \
        dict(cutVtext=kamodo_object.cutVtext, cutV=kamodo_object.cutV, 
             cutLtext=kamodo_object.cutLtext, cutL=kamodo_object.cutL,
             tol=kamodo_object.tol, nt = kamodo_object.nt, nX=kamodo_object.nX, 
             nY=kamodo_object.nY, nZ=kamodo_object.nZ, newt=kamodo_object.newt, 
             newx=kamodo_object.newx, newy=kamodo_object.newy, newz=kamodo_object.newz, 
             nDim=kamodo_object.nDim) 
    return kamodo_object

def dts_to_hrs(datetime_string, filedate):
    '''Get hours since midnight from datetime string, round to avoid boundary issues'''
    return round((datetime.strptime(datetime_string, '%Y-%m-%d %H:%M:%S').replace(tzinfo=timezone.utc)\
            -filedate).total_seconds()/3600.,ndigits=5)

def ellipse_arc(x_center=0, y_center=0, a=1, b=1, start_angle=-np.pi/2, end_angle=np.pi/2, N=100, closed= False):
    '''small function to overlay a partial ellipse on the plot for sun-lit half of earth effect'''\
        
    t = np.linspace(start_angle, end_angle, N)
    x = x_center + a*np.cos(t)
    y = y_center + b*np.sin(t)
    path = f'M {x[0]}, {y[0]}'
    for k in range(1, len(t)):
        path += f'L{x[k]}, {y[k]}'
    if closed:
        path += ' Z'
    return path

#begin generic code for 2D heat plots
def heatplot2D(xint, yint, plot_flip, nx, ny, result, datascale, kamodo_plot, 
               xlabel, ylabel, colorscale, txtbar, xformat, yformat, zformat, 
               xunit, yunit, txttop, txtbot, ellipse=False):
    '''Generate a 2D heat plot for the given data'''    
    
    xrange={'min':xint.min(),'max':xint.max(),'n':len(xint)}
    yrange={'min':yint.min(),'max':yint.max(),'n':len(yint)}
    if plot_flip:
        result2=np.reshape(result,(nx,ny))  # Reshape interpolated values into 2D
    else:
        result2=np.reshape(result,(ny,nx))   #some plots need this
    if nx==ny: result2 = result2.T  #fix weird behavior when nx=ny
    if datascale == "linear":
        def plot_XY(xint = xint, yint = yint):
            return result2
    if datascale == "log":
        def plot_XY(xint = xint, yint = yint):
            return np.log(result2)        
    plotXY = Kamodo(plot_XY = plot_XY)  #func_name = func
    fig = plotXY.plot(plot_XY = dict())        
    if kamodo_plot=='XY':
        fig.update_xaxes(title_text='', scaleanchor='y') # label included in annotations
    if kamodo_plot=='XZ' or kamodo_plot=='YZ' or kamodo_plot=='TY' or kamodo_plot=='TZ':
        fig.update_xaxes(title_text='', scaleanchor=None)
    fig.update_yaxes(title_text=ylabel) 

    # Choose colorscale
    if colorscale == "Rainbow":
        fig.update_traces(
            colorscale=[[0.00, 'rgb(0,0,255)'],
                        [0.25, 'rgb(0,255,255)'],
                        [0.50, 'rgb(0,255,0)'],
                        [0.75, 'rgb(255,255,0)'],
                        [1.00, 'rgb(255,0,0)']]
        )
    elif colorscale == "BlueRed":
        fig.update_traces(colorscale="RdBu", reversescale=True)
    elif colorscale == "Cividis":
        fig.update_traces(colorscale="Cividis")
    else:
        fig.update_traces(colorscale="Viridis")
    fig.update_traces(ncontours=201,
        colorbar=dict(title=txtbar,tickformat=zformat),
            hovertemplate="X: %{x:"+xformat+"}<br>Y: %{y:"+yformat+\
                "}<br><b> %{z:"+zformat+"}</b><extra></extra>",
        contours=dict(coloring="fill", showlines=False)
    )
        
    #set aspect ratio
    #if xunit == yunit:    #real
    if xunit == 'deg':
        if xrange['max']-xrange['min'] > yrange['max']-yrange['min']:
            aspectratio=dict(x=np.asarray([2.,(xrange['max']-xrange['min'])/np.asarray(
                [1.e-5,(yrange['max']-yrange['min'])]).max()]).min(),y=1)
        else:
            aspectratio=dict(x=1,y=np.asarray([1.5,(yrange['max']-yrange['min'])/np.asarray(
                [1.e-5,(xrange['max']-xrange['min'])]).max()]).min())
            if yunit == 'km': aspectratio=dict(x=2, y=1)
        aspectmode='manual'
    else:
        aspectratio=dict(x=1.5,y=1)
        aspectmode='cube'
    
    width=180+300*aspectratio['x']
    height=90+300*aspectratio['y']
    fig.update_layout(
        autosize=False,
        height=height,
        width=width,
        scene=dict(aspectmode=aspectmode),
        margin=dict(t=45,b=45,l=40,r=140),
        title=dict(text=txttop),
        annotations=[
            dict(text=xlabel, x=0.5, y=-0.12, showarrow=False, xref="paper", 
                 yref="paper", font=dict(size=14)),
            dict(text=txtbot,
                 font=dict(size=16, family="sans serif", color="#000000"),
                 x=0.1, y=0.0, ax=0, ay=0, xanchor="left", xshift=-65, yshift=-42, 
                 xref="paper", yref="paper") 
        ])
    if 'Lat' in xlabel or 'Lon' in xlabel:
        fig.update_layout(xaxis = dict(tickmode = 'linear',tick0 = 0,dtick = 30))
    if 'Lat' in ylabel or 'Lon' in ylabel:
        fig.update_layout(yaxis = dict(tickmode = 'linear',tick0 = 0,dtick = 30))      
    if ellipse == True:
        fig.update_layout(
            shapes=[
                dict(type="circle", xref="x", yref="y", x0=-3, y0=-3, x1=3, y1=3, fillcolor="black", line_color="black"),
                dict(type="circle", xref="x", yref="y", x0=-1, y0=-1, x1=1, y1=1, fillcolor="black", line_color="white"),
                dict(type="path", path=ellipse_arc(N=30), fillcolor="white", line_color="white")
            ])
    return fig

def set_plotvar(kamodo_object, datascale, var):
    '''set plot variables specific to plot type'''
    
    if kamodo_object.plottype in ["LonLat_H","LonLat_IP"]:
        txttop=f'{kamodo_object.cutVtext}={kamodo_object.plots[kamodo_object.plottype]["cutV"]:.2f}'+\
            f" slice,  Time = {kamodo_object.timerange['min']}UTC"
        xint, xunit, nx, xlabel = kamodo_object.newz, kamodo_object.zunit, kamodo_object.nZ, 'Lon [deg]'
        yint, yunit, ny, ylabel = kamodo_object.newy, kamodo_object.yunit, kamodo_object.nY, 'Lat [deg]'
        xformat, yformat, zformat, kamodo_plot = '.2f', '.2f', '.4g', 'XY'
    elif kamodo_object.plottype == "LonH":    
        txttop=f'{kamodo_object.cutLtext}={kamodo_object.plots[kamodo_object.plottype]["cutL"]:.2f}'+\
            f" slice,  Time = {kamodo_object.timerange['min']}UTC"
        xint, xunit, nx, xlabel = kamodo_object.newz, kamodo_object.zunit, kamodo_object.nZ, 'Lon [deg]'
        yint, yunit, ny, ylabel = kamodo_object.newx, kamodo_object.xunit, kamodo_object.nX, 'Height [km]'    
        xformat, yformat, zformat, kamodo_plot = '.2f', '.4g', '.4g', 'XZ'   
    elif kamodo_object.plottype == "LonIP":
        txttop=f'{kamodo_object.cutLtext}={kamodo_object.plots[kamodo_object.plottype]["cutL"]:.2f}'+\
            f" slice,  Time = {kamodo_object.timerange['min']}UTC"
        xint, xunit, nx, xlabel = kamodo_object.newz, kamodo_object.zunit, kamodo_object.nZ, 'Lon [deg]'
        yint, yunit, ny, ylabel = kamodo_object.newx, kamodo_object.xunit, kamodo_object.nX, 'IP '    
        xformat, yformat, zformat, kamodo_plot = '.2f', '.2f', '.4g', 'XZ'
    elif kamodo_object.plottype == 'LatH':
        txttop=f'{kamodo_object.cutLtext}={kamodo_object.plots[kamodo_object.plottype]["cutL"]:.2f}'+\
            f" slice,  Time = {kamodo_object.timerange['min']}UTC"
        xint, xunit, nx, xlabel = kamodo_object.newy, kamodo_object.yunit, kamodo_object.nY, 'Lat [deg]'
        yint, yunit, ny, ylabel = kamodo_object.newx, kamodo_object.xunit, kamodo_object.nX, 'Height [km]'
        xformat, yformat, zformat, kamodo_plot = '.2f', '.4g', '.4g', 'YZ'
    elif kamodo_object.plottype == 'LatIP':
        txttop=f'{kamodo_object.cutLtext}={kamodo_object.plots[kamodo_object.plottype]["cutL"]:.2f}'+\
            f" slice,  Time = {kamodo_object.timerange['min']}UTC"
        xint, xunit, nx, xlabel = kamodo_object.newy, kamodo_object.yunit, kamodo_object.nY, 'Lat [deg]'
        yint, yunit, ny, ylabel = kamodo_object.newx, kamodo_object.xunit, kamodo_object.nX, 'IP '
        xformat, yformat, zformat, kamodo_plot = '.2f', '.2f', '.4g', 'YZ'
    elif kamodo_object.plottype in ["TimeLon_H","TimeLon_IP"]:
        txttop=f'{kamodo_object.cutVtext}={kamodo_object.plots[kamodo_object.plottype]["cutV"]:.2f}'+\
            f" slice, {kamodo_object.cutLtext}={kamodo_object.plots[kamodo_object.plottype]['cutL']:.2f} slice"
        xint, xunit, nx = kamodo_object.newt, kamodo_object.tunit, kamodo_object.nt
        xlabel = 'Hours on '+kamodo_object.filedate.isoformat(sep=' ')[0:10]+' [UTC]'
        yint, yunit, ny, ylabel = kamodo_object.newz, kamodo_object.zunit, kamodo_object.nZ, 'Lon [deg]'
        xformat, yformat, zformat, kamodo_plot = '.2f', '.2f', '.4g', 'TY'
    elif kamodo_object.plottype in ["TimeLat_H","TimeLat_IP"]:
        txttop=f'{kamodo_object.cutVtext}={kamodo_object.plots[kamodo_object.plottype]["cutV"]:.2f}'+\
            f" slice, {kamodo_object.cutLtext}={kamodo_object.plots[kamodo_object.plottype]['cutL']:.2f} slice"
        xint, xunit, nx = kamodo_object.newt, kamodo_object.tunit, kamodo_object.nt
        xlabel = 'Hours on '+kamodo_object.filedate.isoformat(sep=' ')[0:10]+' [UTC]'
        yint, yunit, ny, ylabel = kamodo_object.newy, kamodo_object.yunit, kamodo_object.nY, 'Lat [deg]'
        xformat, yformat, zformat, kamodo_plot = '.2f', '.2f', '.4g', 'TY'
    elif kamodo_object.plottype == 'TimeH':
        txttop=f'{kamodo_object.cutVtext}={kamodo_object.plots[kamodo_object.plottype]["cutV"]:.2f}'+\
            f" slice, {kamodo_object.cutLtext}={kamodo_object.plots[kamodo_object.plottype]['cutL']:.2f} slice"
        xint, xunit, nx = kamodo_object.newt, kamodo_object.tunit, kamodo_object.nt
        xlabel = 'Hours on '+kamodo_object.filedate.isoformat(sep=' ')[0:10]+' [UTC]'
        yint, yunit, ny, ylabel = kamodo_object.newx, kamodo_object.xunit, kamodo_object.nX, 'Height [km]'
        xformat, yformat, zformat, kamodo_plot = '.2f', '.4g', '.4g', 'YZ'
    elif kamodo_object.plottype == 'TimeIP':
        txttop=f'{kamodo_object.cutVtext}={kamodo_object.plots[kamodo_object.plottype]["cutV"]:.2f}'+\
            f" slice, {kamodo_object.cutLtext}={kamodo_object.plots[kamodo_object.plottype]['cutL']:.2f} slice"
        xint, xunit, nx = kamodo_object.newt, kamodo_object.tunit, kamodo_object.nt
        xlabel = 'Hours on '+kamodo_object.filedate.isoformat(sep=' ')[0:10]+' [UTC]'
        yint, yunit, ny, ylabel = kamodo_object.newx, kamodo_object.xunit, kamodo_object.nX, 'IP '
        xformat, yformat, zformat, kamodo_plot = '.2f', '.2f', '.4g', 'YZ'
    else:
        print(f'Plot type {kamodo_object.plottype} not yet working.')
        return
    if datascale == "linear":  txtbar=var + " [" + kamodo_object.variables[var]['units'] + "]"
    if datascale == "log":  txtbar="log("+var + ") [" + kamodo_object.variables[var]['units'] + "]"
    result=kamodo_object[var](kamodo_object.newgrid)  # Get values from interpolation already computed
        
    return xint, yint, nx, ny, kamodo_plot, xlabel, ylabel, xformat, yformat, \
        zformat, xunit, yunit, txttop, txtbar, result

def check_plot_dict(given_dict, default_dict, dict_name, dict_unit):
    '''Compare and correct given_dict to be within default_dict max/min'''

    if 'min' not in given_dict.keys(): 
        given_dict['min']=default_dict['min']
    elif given_dict['min']<default_dict['min'] or given_dict['min']>default_dict['max']: 
        given_dict['min']=get_cut(default_dict, dict_name,'min',dict_unit)
    if 'max' not in given_dict.keys(): 
        given_dict['max']=default_dict['max']
    elif given_dict['max']<default_dict['min'] or given_dict['max']>default_dict['max']: 
        given_dict['max']=get_cut(default_dict, dict_name,'max',dict_unit)
    if 'n' not in given_dict.keys(): 
        given_dict['n']=default_dict['n']
    return given_dict

def check_plot_defaults(kamodo_object, var, timerange, lonrange, latrange, htrange):
    '''Check and correct plot defaults to be within allowed ranges for data'''
    
    #set defaults if value(s) not given, esp for timerange
    if 'min' not in timerange.keys():
        timerange['min']=kamodo_object.timerange0['min']   #keep as strings here
    else:  #check that given time is in range
        testmin = dts_to_hrs(timerange['min'], kamodo_object.filedate)
        if testmin<0. or testmin>kamodo_object._time[-1]:  #hrs since midnight
            timerange['min']=get_cut(kamodo_object.timerange0, 'Time','min','')
    if 'max' not in timerange.keys():
        timerange['max']=kamodo_object.timerange0['max']   #convert to timestamps later
    else:  #check that given time is in range
        testmax = dts_to_hrs(timerange['max'], kamodo_object.filedate)
        if testmax<0. or testmax>kamodo_object._time[-1]:  #timestamp values
            timerange['max']=get_cut(kamodo_object.timerange0, 'Time','max','')
    if 'n' not in timerange.keys(): 
        timerange['n']=kamodo_object.timerange0['n']
    lonrange = check_plot_dict(lonrange, kamodo_object.lonrange0, 'Longitude', 'deg')
    latrange = check_plot_dict(latrange, kamodo_object.latrange0, 'Latitude', 'deg') 
    if 'ilev' in kamodo_object.variables[var]['xvec'].keys():
        htrange = check_plot_dict(htrange, kamodo_object.iprange0, 'Pressure Level', '')
    else:
        htrange = check_plot_dict(htrange, kamodo_object.htrange0, 'Height', 'km')
    return timerange, lonrange, latrange, htrange
    
def set_2Dplot_data(kamodo_object, plottype, cutV, cutL, timerange, lonrange, latrange, 
                  htrange):
    '''Set plot data for each plot type'''

    if plottype in ["LonLat_IP",'LonLat_H']:
        kamodo_object.plottype = plottype
        if plottype=='LonLat_IP': 
            kamodo_object.cutVtext, kamodo_object.dz, kamodo_object.dzunit, kamodo_object.numZ = \
                'IP', kamodo_object.dilev, '', len(kamodo_object._ilev)
            #check that cutV is in data range, adjust if needed
            if cutV<htrange['min']: 
                cutV=get_cut(kamodo_object.iprange0, 'cutV','average','')
            elif cutV>htrange['max']:
                cutV=get_cut(kamodo_object.iprange0, 'cutV','max','')           
        if plottype=='LonLat_H': 
            kamodo_object.cutVtext, kamodo_object.dz, kamodo_object.dzunit, kamodo_object.numZ = \
                'Height [km]', kamodo_object.dheight, 'km', len(kamodo_object._height)
            #check that cutV is in data range
            if cutV<htrange['min']: 
                cutV=get_cut(kamodo_object.htrange0, 'cutV','average','km')
            elif cutV>htrange['max']:
                cutV=get_cut(kamodo_object.htrange0, 'cutV','max','km')  
        kamodo_object.cutV, kamodo_object.newx = cutV, cutV
        kamodo_object.nt, kamodo_object.nX, kamodo_object.nY, kamodo_object.nZ = 1, 1, latrange['n'], lonrange['n']
        kamodo_object.newt = dts_to_hrs(timerange['min'], kamodo_object.filedate)
        kamodo_object.newy = np.linspace(latrange['min'],latrange['max'],latrange['n'])
        kamodo_object.newz = np.linspace(lonrange['min'],lonrange['max'],lonrange['n'])
    elif plottype == "LonIP":
        kamodo_object.plottype = plottype
        #check that cutL is in data range
        if cutL<latrange['min']: 
            cutL=get_cut(kamodo_object.latrange0, 'cutL','min','deg')
        elif cutL>latrange['max']:
            cutL=get_cut(kamodo_object.latrange0, 'cutL','max','deg') 
        kamodo_object.cutLtext, kamodo_object.cutL, kamodo_object.newy = 'Lat [deg]', cutL, cutL
        kamodo_object.dz, kamodo_object.dzunit, kamodo_object.numZ = kamodo_object.dilev, '', len(kamodo_object._ilev)
        kamodo_object.nt, kamodo_object.nX, kamodo_object.nY, kamodo_object.nZ = 1, htrange['n'], 1, lonrange['n']
        kamodo_object.newt = dts_to_hrs(timerange['min'], kamodo_object.filedate)
        kamodo_object.newx=np.linspace(htrange['min'],htrange['max'],htrange['n'])
        kamodo_object.newz=np.linspace(lonrange['min'],lonrange['max'],lonrange['n'])
    elif plottype == "LonH":
        kamodo_object.plottype = plottype
        #check that cutL is in data range
        if cutL<latrange['min']: 
            cutL=get_cut(kamodo_object.latrange0, 'cutL','min','deg')
        elif cutL>latrange['max']:
            cutL=get_cut(kamodo_object.latrange0, 'cutL','max','deg') 
        kamodo_object.cutLtext, kamodo_object.cutL, kamodo_object.newy = 'Lat [deg]', cutL, cutL
        kamodo_object.dz, kamodo_object.dzunit, kamodo_object.numZ = kamodo_object.dheight, 'km', len(kamodo_object._height)
        kamodo_object.nt, kamodo_object.nX, kamodo_object.nY, kamodo_object.nZ = 1, htrange['n'], 1, lonrange['n']
        kamodo_object.newt = dts_to_hrs(timerange['min'], kamodo_object.filedate)
        kamodo_object.newx=np.linspace(htrange['min'],htrange['max'],htrange['n'])
        kamodo_object.newz=np.linspace(lonrange['min'],lonrange['max'],lonrange['n'])
    elif plottype == "LatIP":
        kamodo_object.plottype = plottype
        #check that cutL is in data range
        if cutL<lonrange['min']: 
            cutL=get_cut(kamodo_object.lonrange0, 'cutL', 'min', 'deg')
        elif cutL>lonrange['max']:
            cutL=get_cut(kamodo_object.lonrange0, 'cutL', 'max', 'deg') 
        kamodo_object.cutLtext, kamodo_object.cutL, kamodo_object.newz = 'Lon [deg]', cutL, cutL
        kamodo_object.dz, kamodo_object.dzunit, kamodo_object.numZ = kamodo_object.dilev, '', len(kamodo_object._ilev)
        kamodo_object.nt, kamodo_object.nX, kamodo_object.nY, kamodo_object.nZ = 1, htrange['n'], latrange['n'], 1
        kamodo_object.newt = dts_to_hrs(timerange['min'], kamodo_object.filedate)
        kamodo_object.newx = np.linspace(htrange['min'],htrange['max'],htrange['n'])
        kamodo_object.newy = np.linspace(latrange['min'],latrange['max'],latrange['n'])
    elif plottype == "LatH":
        kamodo_object.plottype = plottype
        #check that cutL is in data range
        if cutL<lonrange['min']: 
            cutL = get_cut(kamodo_object.lonrange0, 'cutL', 'min', 'deg')
        elif cutL>lonrange['max']:
            cutL = get_cut(kamodo_object.lonrange0, 'cutL', 'max', 'deg')
        kamodo_object.cutLtext, kamodo_object.cutL, kamodo_object.newz = 'Lon [deg]', cutL, cutL
        kamodo_object.dz, kamodo_object.dzunit, kamodo_object.numZ = kamodo_object.dheight, 'km', len(kamodo_object._height)
        kamodo_object.nt, kamodo_object.nX, kamodo_object.nY, kamodo_object.nZ = 1, htrange['n'], latrange['n'], 1
        kamodo_object.newt = dts_to_hrs(timerange['min'], kamodo_object.filedate)
        kamodo_object.newx = np.linspace(htrange['min'],htrange['max'],htrange['n'])
        kamodo_object.newy = np.linspace(latrange['min'],latrange['max'],latrange['n'])
    elif plottype in ["TimeLon_IP",'TimeLon_H']:   #new
        kamodo_object.plottype = plottype
        if plottype=='TimeLon_IP': 
            kamodo_object.cutVtext, kamodo_object.dz, kamodo_object.dzunit, kamodo_object.numZ = \
                'IP', kamodo_object.dilev, '', len(kamodo_object._ilev)
            #check that cutV is in data range
            if cutV<htrange['min']: 
                cutV=get_cut(kamodo_object.iprange0, 'cutV','average','')
            elif cutV>htrange['max']:
                cutV=get_cut(kamodo_object.iprange0, 'cutV','max','')     
        if plottype=='TimeLon_H': 
            kamodo_object.cutVtext, kamodo_object.dz, kamodo_object.dzunit, kamodo_object.numZ = \
                'Height [km]', kamodo_object.dheight, 'km', len(kamodo_object._height)
            #check that cutV is in data range
            if cutV<htrange['min']: 
                cutV=get_cut(kamodo_object.htrange0, 'cutV','average','km')
            elif cutV>htrange['max']:
                cutV=get_cut(kamodo_object.htrange0, 'cutV','max','km')    
        #check that cutL is in data range
        if cutL<latrange['min']: 
            cutL=get_cut(kamodo_object.latrange0, 'cutL','min','deg')
        elif cutL>latrange['max']:
            cutL=get_cut(kamodo_object.latrange0, 'cutL','max','deg') 
        kamodo_object.cutLtext, kamodo_object.cutL, kamodo_object.newy = 'Lat [deg]', cutL, cutL
        kamodo_object.cutV, kamodo_object.newx = cutV, cutV
        kamodo_object.nt, kamodo_object.nX, kamodo_object.nY, kamodo_object.nZ = timerange['n'], 1, 1, lonrange['n']
        kamodo_object.newt = np.linspace(dts_to_hrs(timerange['min'], kamodo_object.filedate), 
                                dts_to_hrs(timerange['max'], kamodo_object.filedate),
                                timerange['n'], dtype=float)
        kamodo_object.newz = np.linspace(lonrange['min'],lonrange['max'],lonrange['n'])       
    elif plottype in ["TimeLat_IP",'TimeLat_H']:   #new
        kamodo_object.plottype = plottype
        if plottype=='TimeLat_IP': 
            kamodo_object.cutVtext, kamodo_object.dz, kamodo_object.dzunit, kamodo_object.numZ = \
                'IP', kamodo_object.dilev, '', len(kamodo_object._ilev)
            #check that cutV is in data range
            if cutV<htrange['min']: 
                cutV=get_cut(kamodo_object.iprange0, 'cutV','average','')
            elif cutV>htrange['max']:
                cutV=get_cut(kamodo_object.iprange0, 'cutV','max','')   
        if plottype == 'TimeLat_H': 
            kamodo_object.cutVtext, kamodo_object.dz, kamodo_object.dzunit, kamodo_object.numZ = \
                'Height [km]', kamodo_object.dheight, 'km', len(kamodo_object._height)
            #check that cutV is in data range
            if cutV<htrange['min']: 
                cutV=get_cut(kamodo_object.htrange0, 'cutV','average','km')
            elif cutV>htrange['max']:
                cutV=get_cut(kamodo_object.htrange0, 'cutV','max','km')  
        #check that cutL is in data range
        if cutL<lonrange['min']: 
            cutL = get_cut(kamodo_object.lonrange0, 'cutL', 'min', 'deg')
        elif cutL>lonrange['max']:
            cutL = get_cut(kamodo_object.lonrange0, 'cutL', 'max', 'deg') 
        kamodo_object.cutLtext, kamodo_object.cutL, kamodo_object.newz = 'Lon [deg]', cutL, cutL
        kamodo_object.cutV, kamodo_object.newx = cutV, cutV
        kamodo_object.nt, kamodo_object.nX, kamodo_object.nY, kamodo_object.nZ = timerange['n'], 1, latrange['n'], 1
        kamodo_object.newt = np.linspace(dts_to_hrs(timerange['min'], kamodo_object.filedate), 
                                dts_to_hrs(timerange['max'], kamodo_object.filedate),
                                timerange['n'], dtype=float)
        kamodo_object.newy = np.linspace(latrange['min'],latrange['max'],latrange['n'])  
    elif plottype == "TimeIP":   #new
        kamodo_object.plottype = plottype
        #check that cutV is in data range
        if cutL<latrange['min']: 
            cutL=get_cut(kamodo_object.latrange0, 'cutL','min','deg')
        elif cutL>latrange['max']:
            cutL=get_cut(kamodo_object.latrange0, 'cutL','max','deg') 
        kamodo_object.cutVtext, kamodo_object.cutV, kamodo_object.newy = 'Lat [deg]', cutV, cutV
        #check that cutL is in data range
        if cutL<lonrange['min']: 
            cutL = get_cut(kamodo_object.lonrange0, 'cutL', 'min', 'deg')
        elif cutL>lonrange['max']:
            cutL = get_cut(kamodo_object.lonrange0, 'cutL', 'max', 'deg')
        kamodo_object.cutLtext, kamodo_object.cutL, kamodo_object.newz = 'Lon [deg]', cutL, cutL
        kamodo_object.dz, kamodo_object.dzunit, kamodo_object.numZ = kamodo_object.dilev, '', len(kamodo_object._ilev)
        kamodo_object.nt, kamodo_object.nX, kamodo_object.nY, kamodo_object.nZ = timerange['n'], htrange['n'], 1, 1
        kamodo_object.newt = np.linspace(dts_to_hrs(timerange['min'], kamodo_object.filedate), 
                                dts_to_hrs(timerange['max'], kamodo_object.filedate),
                                timerange['n'], dtype=float)
        kamodo_object.newx = np.linspace(htrange['min'],htrange['max'],htrange['n'])   
    elif plottype == "TimeH":   #new
        kamodo_object.plottype = plottype
        #check that cutV is in data range
        if cutV<latrange['min']: 
            cutV = get_cut(kamodo_object.latrange0, 'cutV', 'min', 'deg') 
        elif cutV>latrange['max']:
            cutV = get_cut(kamodo_object.latrange0, 'cutV', 'max', 'deg') 
        kamodo_object.cutVtext, kamodo_object.cutV, kamodo_object.newy = 'Lat [deg]', cutV, cutV
        #check that cutL is in data range
        if cutL<lonrange['min']: 
            cutL = get_cut(kamodo_object.lonrange0, 'cutL', 'min', 'deg')
        elif cutL>lonrange['max']:
            cutL = get_cut(kamodo_object.lonrange0, 'cutL', 'max', 'deg')
        kamodo_object.cutLtext, kamodo_object.cutL, kamodo_object.newz = 'Lon [deg]', cutL, cutL
        kamodo_object.dz, kamodo_object.dzunit, kamodo_object.numZ = kamodo_object.dheight, 'km', len(kamodo_object._height)
        kamodo_object.nt, kamodo_object.nX, kamodo_object.nY, kamodo_object.nZ = timerange['n'], htrange['n'], 1, 1
        kamodo_object.newt = np.linspace(dts_to_hrs(timerange['min'], kamodo_object.filedate), 
                                dts_to_hrs(timerange['max'], kamodo_object.filedate),
                                timerange['n'], dtype=float)
        kamodo_object.newx = np.linspace(htrange['min'],htrange['max'],htrange['n'])  
    else:
        print('Error, unknown plottype. ',plottype)
        return   
     
    #store values
    kamodo_object.plots[plottype] = dict(cutVtext=kamodo_object.cutVtext, cutV=kamodo_object.cutV, 
                                cutLtext=kamodo_object.cutLtext, cutL=kamodo_object.cutL,
                                tol=kamodo_object.tol, nt=kamodo_object.nt, nX=kamodo_object.nX, nY=kamodo_object.nY, 
                                nZ=kamodo_object.nZ, newt=kamodo_object.newt, newx=kamodo_object.newx, 
                                newy=kamodo_object.newy, newz=kamodo_object.newz)
    for varname in kamodo_object.variables: #retrieve interpolators
        kamodo_object.plots[plottype][varname]=kamodo_object[varname] #kamodo_object.variables[varname]['interpolator']

    return kamodo_object

def if_new_plot(kamodo_object, var, plottype, cutV, cutL, timerange, lonrange, 
                latrange, htrange):
    '''Determine is new plot values are needed. If so, generate and store them.'''

    #compare plot data to defaults
    timerange, lonrange, latrange, htrange = check_plot_defaults(
        kamodo_object, var, timerange, lonrange, latrange, htrange)
    
    
    #correct plottype according to height dependency
    new_plottype = correct_2Dplot_type(plottype, kamodo_object.variables[var]['xvec'])
    if new_plottype==0: return 1  #if height-dependent plot requested for height independent variable
    else: plottype = new_plottype  #otherwise, adjust plottype for height vs ilev
    
    #compare previous values to new values to see if plot needs to be redone
    plotpar_listA = [kamodo_object.plottype, kamodo_object.cutV, kamodo_object.cutL, \
                    kamodo_object.lonrange, kamodo_object.latrange, kamodo_object.timerange, \
                    kamodo_object.htrange, kamodo_object.nDim]
    nDim = len(kamodo_object.variables[var]['xvec'].keys())  #check that num of dimensions in var is the same
    plotpar_listB = [plottype, cutV, cutL, lonrange, latrange, timerange,\
                     htrange, nDim]        
    plotpar_comp = [0 if aa==bb else 1 for aa, bb in zip(plotpar_listA, plotpar_listB)]
    if sum(plotpar_comp)==0:
        print('Plottype and variables are unchanged, returning.')
        return 0
    
    #prepare values/data for new plot
    kamodo_object.timerange,kamodo_object.lonrange,kamodo_object.latrange,\
        kamodo_object.htrange = timerange,lonrange,latrange,htrange  #set new values if changed
    kamodo_object.plottype, kamodo_object.cutV, kamodo_object.cutL, kamodo_object.nDim \
        = plottype, cutV, cutL, nDim
    kamodo_object = set_2Dplot_data(kamodo_object, plottype, cutV, cutL, timerange, lonrange, 
                            latrange, htrange) #prepare plot data
    return kamodo_object

def correct_2Dplot_type(plottype, xvec_dependencies):
    '''Choose appropriate plot variation based on height dependency'''
    
    #determine height dependency
    if 'ilev' in xvec_dependencies.keys(): 
        h='IP'
        print('Variable depends on pressure level.')
    elif 'height' in xvec_dependencies.keys(): 
        h='H'
        print('Variable depends on height in km.')
    else: 
        h='none'  
        print('Variable has no height dependence.')
    
    #possible plottypes: ['LonLat','LonH','LatH','TimeLon','TimeLat','TimeH']
    if 'H' in plottype and h=='none': 
        print('Choose a different plot type.')
        return 0
    elif 'H' in plottype and h=='IP':
        new_plottype = plottype.replace('H','IP')
    elif 'H' in plottype and h!='IP':
        new_plottype = plottype
    elif 'H' not in plottype and h!='IP':
        new_plottype = plottype+'_H'  #LonLat -> LonLat_H
    elif 'H' not in plottype and h=='IP':
        new_plottype = plottype+'_IP'  #LonLat -> LonLat_IP
    #print(plottype, h, new_plottype)
    
    return new_plottype
    