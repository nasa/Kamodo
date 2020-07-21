from plotly.offline import init_notebook_mode, iplot, plot
import plotly.graph_objs as go
#init_notebook_mode(connected = True)
import plotly
import spacepy
import os
import time
import numpy as np
from scipy.interpolate import RegularGridInterpolator

from kamodo.kamodo import Kamodo, kamodofy
from kamodo import readers

from kamodo.readers.gitm import gitm as gold
from kamodo.readers.gitm import gitm_alt_plots as gap
from kamodo.readers.gitm import gitm_plot_rout as gpr

class GITM(Kamodo): 
    def __init__(self, 
                 filename, 
                 runpath = "./", 
                 runname = "noname", 
                 **kwargs):
        # Start timer
        tic = time.perf_counter()
        
        # these three variables are needed to prevent errors in Kamodo
        self.verbose=False;
        self.symbol_registry=dict();
        self.signatures=dict();

        # perform any necessary I/O
        print('opening {}'.format(filename))
        self.filename = filename
        self.runpath = runpath
        self.runname = runname
        self.missing_value = np.NAN
        self.variables=dict()
        gData = gold.GitmBin(filename)
        self.gData = gData
        
        # Get time
        self.time = gData['time']
        print('... simulation time = ',self.time)
        print('... raw data array size = ',gData['Altitude'].shape)
        
        # Pull from 3D arrays into 1D array of unique values
        self.lonkey = "dLon"
        self.latkey = "dLat"
        self.altkey = "Altitude"
        self.lon = np.unique(gData[self.lonkey])
        self.lat = np.unique(gData[self.latkey])
        self.alt = np.unique(gData[self.altkey])
        self.altmin = np.min(self.alt)
        self.altmax = np.max(self.alt)
        print('... range of altitudes is ',self.altmin,' to ',self.altmax,' meters.')
        
        # Just one variable for now
        self.variables['rho'] = dict(units = 'kg/m^3', data = gData['Rho (kg/m3)'], dtype=np.float32)
        self.variables['N_O3P'] = dict(units = '1/m^3', data = gData['O(!U3!NP)            (/m3)'], dtype=np.float32)
        self.variables['N_O2'] = dict(units = '1/m^3', data = gData['O!D2!N               (/m3)'], dtype=np.float32)
        self.variables['N_N2'] = dict(units = '1/m^3', data = gData['N!D2!N               (/m3)'], dtype=np.float32)
        self.variables['N_NO'] = dict(units = '1/m^3', data = gData['NO                   (/m3)'], dtype=np.float32)
        self.variables['Vn_east'] = dict(units = 'm/s', data = gData['Vn (east) (m/s)'], dtype=np.float32)
        self.variables['Vn_north'] = dict(units = 'm/s', data = gData['Vn (north) (m/s)'], dtype=np.float32)
        self.variables['Vn_up'] = dict(units = 'm/s', data = gData['Vn (up) (m/s)'], dtype=np.float32)
        self.variables['N_O4SPplus'] = dict(units = '1/m^3', data = gData['O_4SP_!U+!N          (/m3)'], dtype=np.float32)
        self.variables['N_O2plus'] = dict(units = '1/m^3', data = gData['O!D2!U+!N            (/m3)'], dtype=np.float32)
        self.variables['N_N2plus'] = dict(units = '1/m^3', data = gData['N!D2!U+!N            (/m3)'], dtype=np.float32)
        self.variables['N_NOplus'] = dict(units = '1/m^3', data = gData['NO!U+!N              (/m3)'], dtype=np.float32)
        self.variables['N_eminus'] = dict(units = '1/m^3', data = gData['e-                   (/m3)'], dtype=np.float32)
        self.variables['Vi_east'] = dict(units = 'm/s', data = gData['Vi (east) (m/s)'], dtype=np.float32)
        self.variables['Vi_north'] = dict(units = 'm/s', data = gData['Vi (north) (m/s)'], dtype=np.float32)
        self.variables['Vi_up'] = dict(units = 'm/s', data = gData['Vi (up) (m/s)'], dtype=np.float32)
        self.variables['Tn'] = dict(units = 'K', data = gData['Neutral Temperature (K)'], dtype=np.float32)

        for varname in self.variables:
            units = self.variables[varname]['units']
            self.register_variable(varname, units)
        
        # end timer
        toc = time.perf_counter()
        print(f"Time loading file and kamodifying results: {toc - tic:0.4f} seconds")
        
    def register_variable(self, varname, units):
        """register variables into Kamodo for this model, GITM"""

        ilon = self.lon
        ilat = self.lat
        ialt = self.alt
        interpolator = self.get_grid_interpolator(ilon, ilat, ialt, varname)

        # store the interpolator
        self.variables[varname]['interpolator'] = interpolator

        def interpolate(xvec):  
            return self.variables[varname]['interpolator'](xvec)

        # update docstring for this variable
        interpolate.__doc__ = "A function that returns {} in [{}].".format(varname,units)

        self[varname] = kamodofy(interpolate, 
                                 units = units, 
                                 citation = "De Zeeuw 2020",
                                 data = self.variables[varname]['data'])
        
    def kamodofy_names(name, name_maps):
        """replaces all substrings in name with those given by name_maps"""
        for old, new in name_maps:
            name = name.replace(old, new)
        return name
    
    def get_grid_interpolator(self, lon, lat, alt, varname):
        """create a regular grid interpolator for this variable"""
        data =  self.variables[varname]['data']
        interpolator = RegularGridInterpolator((lon, lat, alt),
                                               data, 
                                               bounds_error = False,
                                               fill_value = self.missing_value)
        return interpolator

    def get_plot(self, var, value, plottype, colorscale="BlueRed", sym="F", log="F"):
        '''
        Return a plotly figure object for the plottype requested.
        var, value, and plottype are required variables. value is position of slice
        colorscale = Viridis, Cividis, Rainbow, or BlueRed
        sym = F [default] for symetric colorscale around 0
        log = F [default] for log10() of plot value
        '''

        if plottype == "2D-alt":
            # Check if altitude entered is valid
            if value < self.altmin or value > self.altmax:
                print('Altitude is out of range: alt=',value,\
                      ' min/max=',self.altmin,'/',self.altmax)
                return

            # Set grid for plot
            altkm=value/1000.
            ilon = np.linspace(0, 360, 361)
            ilat = np.linspace(-90, 90, 181)
            ialt = np.array([value])
            xx, yy = np.meshgrid(np.array(ilon), np.array(ilat))
            grid = np.ndarray(shape=(np.size(np.reshape(xx,-1)),3), dtype=np.float32)
            grid[:,0] = np.reshape(xx,-1)
            grid[:,1] = np.reshape(yy,-1)
            grid[:,2] = value
            units=self.variables[var]['units']
            test = self.variables[var]['interpolator'](grid)
            result = np.reshape(test,(ilat.shape[0],ilon.shape[0]))
            if log == "T":
                result = np.log10(result)
            if sym == "T":
                cmax = np.max(np.absolute(result))
                cmin = -cmax
            else:
                cmax = np.max(result)
                cmin = np.min(result)

            time=self.time.strftime("%Y/%m/%d %H:%M:%S UT")

            def plot_var(lon = ilon, lat = ilat):
                return result
            plotvar = Kamodo(plot_var = plot_var)

            fig = plotvar.plot(plot_var = dict())
            #fig.update_xaxes(nticks=7,title_text="",scaleanchor='y')
            fig.update_xaxes(tick0=0.,dtick=45.,title_text="")
            fig.update_yaxes(tick0=0.,dtick=45,title_text="")
            txtbar = var + " [" + units + "]"
            if log == "T":
                txtbar = "log<br>"+txtbar
            txtbot = "Model: GITM,  Run: " + self.runname
            if colorscale == "BlueRed":
                fig.update_traces(
                    colorscale="RdBu",
                    reversescale=True,
                )
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
            fig.update_traces(
                zmin=cmin, zmax=cmax,
                ncontours=201,
                colorbar=dict(title=txtbar,    tickformat=".3g"),
                contours=dict(coloring="fill", showlines=False)
            )
            if log == "T":
                fig.update_traces(
                    hovertemplate="Lon: %{x:.0f}<br>Lat: %{y:.0f}<br><b>"+"log("+var+"): %{z:.4g}</b><extra></extra>"
                )
            else:
                fig.update_traces(
                    hovertemplate="Lon: %{x:.0f}<br>Lat: %{y:.0f}<br><b>"+var+": %{z:.4g}</b><extra></extra>"
                )                
            fig.update_layout(
                title=dict(text="Altitude="+"{:.0f}".format(altkm)+" km,  Time = " + time,
                           yref="container", yanchor="top", y=0.95),
                title_font_size=16,
                annotations=[
                    dict(text="Lon [degrees]", x=0.5, y=-0.13, showarrow=False,
                         xref="paper", yref="paper", font=dict(size=12)),
                    dict(text="Lat [degrees]", x=-0.1, y=0.5, showarrow=False,
                         xref="paper", yref="paper", font=dict(size=12), textangle=-90),
                    dict(text=txtbot, x=0.0, y=0.0, ax=0, ay=0, xanchor="left",
                         xshift=-65, yshift=-42, xref="paper", yref="paper",
                         font=dict(size=16, family="sans serif", color="#000000"))
                ],
                height=340
            )
            return fig

        if plottype == "2D-lat":
            # Check if latitude entered is valid
            if value < -90. or value > 90.:
                print('Latitude is out of range: lat=',value,' min/max= -90./90.')
                return

            # Set grid for plot
            ilon = np.linspace(0, 360, 361)
            ilat = np.array([value])
            ialt = np.linspace(98000, 736000, 320)
            xx, yy = np.meshgrid(np.array(ilon), np.array(ialt))
            grid = np.ndarray(shape=(np.size(np.reshape(xx,-1)),3), dtype=np.float32)
            grid[:,0] = np.reshape(xx,-1)
            grid[:,1] = value
            grid[:,2] = np.reshape(yy,-1)
            units=self.variables[var]['units']
            test = self.variables[var]['interpolator'](grid)
            result = np.reshape(test,(ialt.shape[0],ilon.shape[0]))
            if log == "T":
                result = np.log10(result)
            if sym == "T":
                cmax = np.max(np.absolute(result))
                cmin = -cmax
            else:
                cmax = np.max(result)
                cmin = np.min(result)

            time=self.time.strftime("%Y/%m/%d %H:%M:%S UT")

            ialt = ialt/1000.
            def plot_var(lon = ilon, alt = ialt):
                return result
            plotvar = Kamodo(plot_var = plot_var)

            fig = plotvar.plot(plot_var = dict())
            #fig.update_xaxes(nticks=7,title_text="",scaleanchor='y')
            fig.update_xaxes(tick0=0.,dtick=45.,title_text="")
            fig.update_yaxes(title_text="")
            txtbar = var + " [" + units + "]"
            if log == "T":
                txtbar = "log<br>"+txtbar
            txtbot = "Model: GITM,  Run: " + self.runname
            if colorscale == "BlueRed":
                fig.update_traces(
                    colorscale="RdBu",
                    reversescale=True,
                )
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
            fig.update_traces(
                zmin=cmin, zmax=cmax,
                ncontours=201,
                colorbar=dict(title=txtbar,    tickformat=".3g"),
                contours=dict(coloring="fill", showlines=False)
            )
            if log == "T":
                fig.update_traces(
                    hovertemplate="Lon: %{x:.0f}<br>Alt: %{y:.0f}<br><b>"+"log("+var+"): %{z:.4g}</b><extra></extra>"
                )
            else:
                fig.update_traces(
                    hovertemplate="Lon: %{x:.0f}<br>Alt: %{y:.0f}<br><b>"+var+": %{z:.4g}</b><extra></extra>"
                )                
            fig.update_layout(
                title=dict(text="Latitude="+"{:.1f}".format(value)+" degrees,  Time = " + time,
                           yref="container", yanchor="top", y=0.95),
                title_font_size=16,
                annotations=[
                    dict(text="Lon [degrees]", x=0.5, y=-0.13, showarrow=False,
                         xref="paper", yref="paper", font=dict(size=12)),
                    dict(text="Altitude [km]", x=-0.1, y=0.5, showarrow=False,
                         xref="paper", yref="paper", font=dict(size=12), textangle=-90),
                    dict(text=txtbot, x=0.0, y=0.0, ax=0, ay=0, xanchor="left",
                         xshift=-65, yshift=-42, xref="paper", yref="paper",
                         font=dict(size=16, family="sans serif", color="#000000"))
                ],
                height=340
            )
            return fig

        if plottype == "2D-lon":
            # Check if longitude entered is valid
            if value < 0. or value > 360.:
                print('Latitude is out of range: lat=',value,' min/max= 0./360.')
                return

            # Set grid for plot
            ilon = np.array([value])
            ilat = np.linspace(-90, 90, 181)
            ialt = np.linspace(98000, 736000, 320)
            xx, yy = np.meshgrid(np.array(ilat), np.array(ialt))
            grid = np.ndarray(shape=(np.size(np.reshape(xx,-1)),3), dtype=np.float32)
            grid[:,0] = value
            grid[:,1] = np.reshape(xx,-1)
            grid[:,2] = np.reshape(yy,-1)
            units=self.variables[var]['units']
            test = self.variables[var]['interpolator'](grid)
            result = np.reshape(test,(ialt.shape[0],ilat.shape[0]))
            if log == "T":
                result = np.log10(result)
            if sym == "T":
                cmax = np.max(np.absolute(result))
                cmin = -cmax
            else:
                cmax = np.max(result)
                cmin = np.min(result)

            time=self.time.strftime("%Y/%m/%d %H:%M:%S UT")

            ialt = ialt/1000.
            def plot_var(lat = ilat, alt = ialt):
                return result
            plotvar = Kamodo(plot_var = plot_var)

            fig = plotvar.plot(plot_var = dict())
            #fig.update_xaxes(nticks=7,title_text="",scaleanchor='y')
            fig.update_xaxes(tick0=0.,dtick=30.,title_text="")
            fig.update_yaxes(title_text="")
            txtbar = var + " [" + units + "]"
            if log == "T":
                txtbar = "log<br>"+txtbar
            txtbot = "Model: GITM,  Run: " + self.runname
            if colorscale == "BlueRed":
                fig.update_traces(
                    colorscale="RdBu",
                    reversescale=True,
                )
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
            fig.update_traces(
                zmin=cmin, zmax=cmax,
                ncontours=201,
                colorbar=dict(title=txtbar,    tickformat=".3g"),
                contours=dict(coloring="fill", showlines=False)
            )
            if log == "T":
                fig.update_traces(
                    hovertemplate="Lat: %{x:.0f}<br>Alt: %{y:.0f}<br><b>"+"log("+var+"): %{z:.4g}</b><extra></extra>"
                )
            else:
                fig.update_traces(
                    hovertemplate="Lat: %{x:.0f}<br>Alt: %{y:.0f}<br><b>"+var+": %{z:.4g}</b><extra></extra>"
                )                
            fig.update_layout(
                title=dict(text="Longitude="+"{:.1f}".format(value)+" degrees,  Time = " + time,
                           yref="container", yanchor="top", y=0.95),
                title_font_size=16,
                annotations=[
                    dict(text="Lat [degrees]", x=0.5, y=-0.13, showarrow=False,
                         xref="paper", yref="paper", font=dict(size=12)),
                    dict(text="Altitude [km]", x=-0.1, y=0.5, showarrow=False,
                         xref="paper", yref="paper", font=dict(size=12), textangle=-90),
                    dict(text=txtbot, x=0.0, y=0.0, ax=0, ay=0, xanchor="left",
                         xshift=-65, yshift=-42, xref="paper", yref="paper",
                         font=dict(size=16, family="sans serif", color="#000000"))
                ],
                height=340
            )
            return fig

        if plottype == "3D-alt":
            # Check if altitude entered is valid
            if value < self.altmin or value > self.altmax:
                print('Altitude is out of range: alt=',value,\
                      ' min/max=',self.altmin,'/',self.altmax)
                return

            # Set grid for plot
            altkm=value/1000.
            ilon = np.linspace(0, 360, 361)
            ilat = np.linspace(-90, 90, 181)
            ialt = np.array([value])
            xx, yy = np.meshgrid(np.array(ilon), np.array(ilat))
            grid = np.ndarray(shape=(np.size(np.reshape(xx,-1)),3), dtype=np.float32)
            grid[:,0] = np.reshape(xx,-1)
            grid[:,1] = np.reshape(yy,-1)
            grid[:,2] = value
            units=self.variables[var]['units']
            test = self.variables[var]['interpolator'](grid)
            result = np.reshape(test,(ilat.shape[0],ilon.shape[0]))
            if log == "T":
                result = np.log10(result)
            r = value + 6.3781E6
            x=(r*np.cos(yy*np.pi/180.)*np.cos(xx*np.pi/180.))/6.3781E6
            y=(r*np.cos(yy*np.pi/180.)*np.sin(xx*np.pi/180.))/6.3781E6
            z=(r*np.sin(yy*np.pi/180.))/6.3781E6
            if sym == "T":
                cmax = np.max(np.absolute(result))
                cmin = -cmax
            else:
                cmax = np.max(result)
                cmin = np.min(result)

            time=self.time.strftime("%Y/%m/%d %H:%M:%S UT")

            def plot_var(x = x, y = y, z = z):
                return result
            plotvar = Kamodo(plot_var = plot_var)

            fig = plotvar.plot(plot_var = dict())
            fig.update_scenes(xaxis=dict(title=dict(text="X [Re]")),
                              yaxis=dict(title=dict(text="Y [Re]")),
                              zaxis=dict(title=dict(text="Z [Re]")))
            txtbar = var + " [" + units + "]"
            if log == "T":
                txtbar = "log<br>"+txtbar
            txtbot = "Model: GITM,  Run: " + self.runname
            if colorscale == "BlueRed":
                fig.update_traces(
                    colorscale="RdBu",
                    reversescale=True,
                )
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
            fig.update_traces(
                cmin=cmin, cmax=cmax,
                colorbar=dict(title=txtbar, tickformat=".3g")
            )
            fig.update_traces(
                hovertemplate="X [Re]: %{x:.3f}<br>Y [Re]: %{y:.3f}<br>Z [Re]: %{z:.3f}<extra></extra>"
            )
            fig.update_layout(
                title=dict(text="Altitude="+"{:.0f}".format(altkm)+" km,  Time = " + time,
                           yref="container", yanchor="top", x=0.01, y=0.95),
                title_font_size=16,
                annotations=[
                    dict(text=txtbot, x=0.0, y=0.0, ax=0, ay=0, xanchor="left",
                         xshift=0, yshift=-20, xref="paper", yref="paper",
                         font=dict(size=16, family="sans serif", color="#000000"))
                ],
                margin=dict(l=0),
                width=600
            )
            x1=[0., 0.]
            y1=[0., 0.]
            z1=[-1.2, 1.2]
            fig.add_scatter3d(mode='lines',x=x1,y=y1,z=z1,line=dict(width=4,color='black'),
                              showlegend=False,hovertemplate='Polar Axis<extra></extra>')
            r = value + 10000. + 6.3781E6
            x2=(r*np.cos(ilon*np.pi/180.))/6.3781E6
            y2=(r*np.sin(ilon*np.pi/180.))/6.3781E6
            z2=0.*ilon
            fig.add_scatter3d(mode='lines',x=x2,y=y2,z=z2,line=dict(width=2,color='black'),
                              showlegend=False,hovertemplate='Equator<extra></extra>')
            x3=(r*np.cos(ilat*np.pi/180.))/6.3781E6
            y3=0.*ilat
            z3=(r*np.sin(ilat*np.pi/180.))/6.3781E6
            fig.add_scatter3d(mode='lines',x=x3,y=y3,z=z3,line=dict(width=2,color='black'),
                              showlegend=False,hovertemplate='Midnight<extra></extra>')
            return fig

        print('Unknown plottype (',plottype,') returning.')
        return
    
