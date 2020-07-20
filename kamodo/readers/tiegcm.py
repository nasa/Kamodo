import numpy as np

from kamodo import Kamodo, kamodofy, gridify
import time
from scipy.interpolate import RegularGridInterpolator, interp1d
from datetime import datetime,timedelta

# pip install pytiegcm
from tiegcm.tiegcm import TIEGCM

def parse_units(varname, variable):
    """Parses units based on variable input"""
    try:
        units = variable.units
        units = units.replace('-','**-')
        units = units.replace('mmr','1').replace('S', 'ohm**-1')
        units = units.replace('cm3','cm**3').replace('erg/g/s', 'erg/(g*s)')
        units = units.replace('none','1')
    except:
        units = ''
        if varname == 'Z':
            units = 'cm'
    return units

def sort_variables(self, variables):
    variables_3d = []
    variables_4d = []
    for varname in variables:
        try:
            varshape = self._tiegcm.rootgrp.variables[varname].shape
        except:
            continue
        if len(varshape) == 4:
            variables_4d.append(varname)
        elif len(varshape) == 3:
            variables_3d.append(varname)
    return variables_4d + variables_3d

class TIEGCM_Kamodo(Kamodo):
    def __init__(self, filename,
                 variables_requested = None,
                 date = None,
                 time = None,
                 debug= False,
                 runpath = "./", 
                 runname = "noname",
                 **kwargs):
        self._tiegcm = TIEGCM(filename)
        
        self._ilev = np.array(self._tiegcm.rootgrp.variables['ilev'])
        self._time = np.array(self._tiegcm.rootgrp.variables['time'])
        self._year=np.array(self._tiegcm.rootgrp.variables['year'])
        self._day=np.array(self._tiegcm.rootgrp.variables['day'])
        
        self._lat = self._tiegcm.lat
        self._lon = self._tiegcm.lon
        self._registered = 0
        
        super(TIEGCM_Kamodo, self).__init__() 
        print('opening {}'.format(filename))

# needed by CCMC online visualization based on Kamodo
        self.filename = filename
        self.runpath = runpath
        self.runname = runname
        self.missing_value = np.NAN
        self.variables=dict()

# Interpolation values
        self.plots = dict()
        self.plottype = "LonLat" 
        self.cut = 'IP'
        self.cutV = 0. # within all the coordinate ranges
        self.nT = 1
        self.nX = 1
        self.nY = 37
        self.nZ = 73
        self.newt = 1440.
        self.newx = self.cutV
        self.newy = np.linspace(-90., 90., self.nY)
        self.newz = np.linspace(0., 360., self.nZ)
        self.xname="IP"
        self.xunit=''
        self.yname="Lat"
        self.yunit='deg'
        self.zname="Lon"
        self.zunit='deg'
        self.lonrange=dict(min=self._lon.min(),max=self._lon.max(),n=len(self._lon))
        self.latrange=dict(min=self._lat.min(),max=self._lat.max(),n=len(self._lat))
        self.hrange=dict(min=80,max=450,n=34)
        self._H=np.linspace(self.hrange['min'],self.hrange['max'],self.hrange['n'])
        self.tol = 1.1
        self.plots[self.plottype] = dict(cut=self.cut, cutV=self.cutV, tol=self.tol,
                                         nT=self.nT,nX=self.nX, nY=self.nY, nZ=self.nZ,
                                         newt=self.newt,
                                         newx=self.newx,
                                         newy=self.newy,
                                         newz=self.newz)  
        
        self.gridSize=len(self._ilev)*len(self._lat)*len(self._lon)
        self.gridMinDx=np.asarray([self._lon[1]-self._lon[0],self._lat[1]-self._lat[0]]).min()
        # Get grids ready to use
        self.setup_interpolating_grids()
        
        
        if variables_requested is None:
            variables = self._tiegcm.rootgrp.variables.keys()
        else:
            variables = variables_requested

        variables = sort_variables(self, variables)

        for varname in variables:            
            units = parse_units(varname, self._tiegcm.rootgrp.variables[varname])
            try:
                self._tiegcm.set_variable_boundary_condition(varname)
            except:
#                 print('can not set boundary condition for {}, skipping..'.format(varname))
                continue
            self._tiegcm.wrap_variable(varname)
            variable = self._tiegcm.rootgrp.variables[varname]

            if len(variable.shape) not in [3,4]:
                # skip this variable
                continue
                
            if varname == 'ZMAG':
                continue

            elif len(variable.shape) == 4:
                self.register_4d_variable(units, variable, varname)

            elif len(variable.shape) == 3:
                self.register_3d_variable(units, variable, varname)

            self.variables[varname]=dict(data=variable,units=units)

                
        print('registered {} variables'.format(self._registered))
        
        # register user's input variables, assuming kamodo-compatible
#        for varname, variable in kwargs.items():
#            self[varname] = variable

#        def ilev(points):
#            return self.lev(*points)
#
#        self['ilev'] = ilev

    # vertical interpolation as is done in CCMC-Vis (3DView) IDL code
    @np.vectorize
# can't do gridify here unless 'self' is somwhow known
#    @gridify(t = self._time, z=self._H, lat = self._lat, lon = self._lon)
    def vert_interp(self,t, z, lat, lon):
    # 1) meshgrid and squeeze the shapes of time, lat, lon if they are not the same
    # 2) calculate z_levels for all points
    # 3) construct rgi for (time, z_level, lat, lon)
    # 4) interpolate over z
        varname=self.varname
        z_levels=np.squeeze(self['Z'](t=t,lat=lat,lon=lon))
        data_levels=np.squeeze(self[varname](t=t,lat=lat,lon=lon))
        interp = interp1d(z_levels, data_levels, bounds_error=False)
        return interp(z)

    def register_3d_variable(self, units, variable, varname):
        """Registers a 3d interpolator with 3d signature"""

        rgi = RegularGridInterpolator((self._time, self._lat, self._lon),
                                      variable, bounds_error = False)
        @kamodofy(units = units, data = variable)
        @gridify(t = self._time, lat = self._lat, lon = self._lon)
        def interpolator(xvec):
            """Interpolates 3d variable"""
            return rgi(xvec)

        self[varname] = interpolator
        self._registered += 1

    def register_4d_variable(self, units, variable, varname):
        """Registers a 4d interpolator with 4d signature"""

        rgi = RegularGridInterpolator((self._time, self._ilev, self._lat, self._lon), 
                                      variable, bounds_error = False)
        @kamodofy(units = units, data = variable)
        @gridify(t = self._time, ilev = self._ilev, lat = self._lat, lon = self._lon) 
#        @np.vectorize
        def interpolator(xvec):
            """Interpolates 4d variable"""
            return rgi(xvec)
        self.variables['varname']=dict(data=variable,units=units,interpolator=rgi)
        self[varname] = interpolator
        self._registered += 1

    @np.vectorize
    def lev(self, t, z, lat, lon):
        """Finds ilev for a given height"""
        # 2) calculate z_levels for all points
        # 3) construct rgi for (time, z_level, lat, lon)
        # 4) interpolate over z
        
        z_levels = np.squeeze(self.Z(t = t, lat = lat, lon = lon)) # ilev has a default
        level = interp1d(z_levels, self._ilev, bounds_error=False)
        return level(z)

    def set_plot(self,
                 plottype = "XY",
                 time_of_day="12:00:00",
                 date=None,
                 cutV = 10,
                 lonrange=dict(min=0,max=360,n=73), # reference to self._lon_density not possible?
                 latrange=dict(min=-90,max=90,n=37),
                 hrange=dict(min=1,max=15,n=15) ):
        '''Set plotting variables for available preset plot types: XY, YZ, YZ, XYZ'''
        if 'min' not in lonrange.keys():
            lonrange['min']=self.lonrange['min']
        if 'max' not in lonrange.keys():
            lonrange['max']=self.lonrange['max']
        if 'n' not in lonrange.keys():
            lonrange['n']=self.lonrange['n']
        if 'min' not in latrange.keys():
            latrange['min']=self.latrange['min']
        if 'max' not in latrange.keys():
            latrange['max']=self.latrange['max']
        if 'n' not in latrange.keys():
            latrange['n']=self.latrange['n']
        if 'min' not in hrange.keys():
            hrange['min']=self.hrange['min']
        if 'max' not in hrange.keys():
            hrange['max']=self.hrange['max']
        if 'n' not in hrange.keys():
            hrange['n']=self.hrange['n']
            
        tic = time.perf_counter()
        if plottype == self.plottype:
            if cutV == self.cutV:
               if lonrange['min'] == self.lonrange['min'] and lonrange['max'] == self.lonrange['max'] and lonrange['n'] == self.lonrange['n']:
                       if latrange['min'] == self.latrange['min'] and latrange['max'] == self.latrange['max'] and latrange['n'] == self.latrange['n']:
                               if hrange['min'] == self.hrange['min'] and hrange['max'] == self.hrange['max'] and hrange['n'] == self.hrange['n']:
                       
                                       print('Plottype (',plottype,') and cut value (',cutV,') are unchanged, returning.')
                                       return
                               
        self.lonrange=lonrange
        self.latrange=latrange
        self.hrange=hrange
        ip_max=70
        h_km_max=7000
        h_m_max=7e6
                           
        if date is None or time_of_day is None:
            datetime0=datetime(int(self._year[0]),1,1)
            datetime1=datetime0+timedelta(days=int(self._day[0]))
            self.datetime=datetime1+timedelta(minutes=int(self._time[0]))
            self.time_of_day=(datetime1-datetime0).total_seconds()/60.
        else:
            self.plotdate=date
            self.plottime=time_of_day            
            datetime0=datetime.strptime(date,"%Y/%m/%d")
            datetime1=datetime.strptime(date+" "+time_of_day,"%Y/%m/%d %H:%M:%S")
            self.datetime=datetime1
            self.time_of_day=(datetime1-datetime0).total_seconds()/60.


        self.filetime=self.datetime.strftime("%Y/%m/%d %H:%M:%S")
        self.date=self.datetime.strftime("%Y/%m/%d")
        self.plottime=self.datetime.strftime("%H:%M:%S")
            
        self.filetime=self.datetime.strftime("%Y/%m/%d %H:%M:%S")
# time selection
        self.nT = 1
#        self.newt = self.plottime
        self.newt = self.time_of_day
# X=IP or Height
# Y=Lat
# Z=Lon
        if plottype == "LonLat" or plottype == "LonLatH":
# auto-detect whether we are using IP or H (or H in cm or km)
            if cutV < ip_max:
                self.cut="IP"
                self.cutunit=" []"
                self.cutfactor=1.
            else:
                self.cut="H"
                if cutV < h_km_max:
                    self.cutunit=" [km]"
                    self.cutfactor=1e5
                elif cutV < h_m_max:
                    self.cutunit=" [m]"
                    self.cutfactor=100.
                else:
                    self.cutunit=" [cm]"
                    self.cutfactor=1.
            self.plottype = plottype
            self.cutV = cutV
            self.nX = 1
            self.nY = latrange['n']
            self.nZ = lonrange['n']
            self.newx = cutV*self.cutfactor
            self.newy = np.linspace(latrange['min'],latrange['max'],latrange['n'])
            self.newz = np.linspace(lonrange['min'],lonrange['max'],lonrange['n'])
        else:
            self.cutunit="[deg]" #both for Lon=const (LatIP, LatH) and Lat=const (LonIP, LonH)
# auto-detect whether we are using IP or H (or H in cm or km)
            if hrange['max'] < ip_max:
                self.xname="IP"
                self.xunit="[]"
                self.xfactor=1.
                self.xformat=".4f"
            else:
                self.xname="H"                
                self.xformat=".4g"
                if hrange['max'] < h_km_max:
                    self.xunit=" [km]"
                    self.xfactor=1e5
                elif hrange['max'] < h_m_max:
                    self.xunit=" [m]"
                    self.xfactor=100.
                else:
                    self.xunit="[cm]"
                    self.xfactor=1.

            self.plottype = plottype
            if plottype == "LonIP" or plottype == "LonH":
                self.cut = 'Lat'
                self.cutV = cutV
                self.nX = hrange['n']
                self.nY = 1
                self.nZ = lonrange['n']
                self.newx=np.linspace(hrange['min'],hrange['max'],hrange['n'])*self.xfactor
                self.newy = cutV
                self.newz=np.linspace(lonrange['min'],lonrange['max'],lonrange['n'])
            elif plottype == "LatIP" or plottype == "LatH":
                self.plottype = plottype
                self.cut = 'Lon'
                self.cutV = cutV
                self.nX=hrange['n']
                self.nY=latrange['n']
                self.nZ=1
                self.newx = np.linspace(hrange['min'],hrange['max'],hrange['n'])*self.xfactor
                self.newy = np.linspace(latrange['min'],latrange['max'],latrange['n'])
                self.newz = cutV
            else:
                print('Error, unknown plottype. ',plottype)
                return

        self.plots[plottype] = dict(cut=self.cut, cutV=self.cutV, tol=self.tol,
                                    nT=self.nT,nX=self.nX, nY=self.nY, nZ=self.nZ,
                                    newt=self.newt,newx=self.newx, newy=self.newy, newz=self.newz)
        self.setup_interpolating_grids()
#        self.reinterpolate_values()

#        for varname in self.variables:
#            self.plots[plottype][varname]=self.variables[varname]['interpolator']

        toc = time.perf_counter()
        print(f"Time resetting plot and precomputing interpolations: {toc - tic:0.4f} seconds")
        return

    def setup_interpolating_grids(self):
        '''setup the grid to interpolate to, trim to necessary size of source grid, and compute interpolation weights'''
        tt, xx, yy, zz = np.meshgrid(self.newt,self.newx,self.newy,self.newz, indexing = 'xy')
        self.newgrid = np.ndarray(shape=(np.size(np.reshape(xx,-1)),4), dtype=np.float32)
        self.newgrid[:,0] = np.reshape(tt,-1)
        self.newgrid[:,1] = np.reshape(xx,-1)
        self.newgrid[:,2] = np.reshape(yy,-1)
        self.newgrid[:,3] = np.reshape(zz,-1)
        print("newgrid shape: ",self.newgrid.shape)

    def get_plot(self, var, colorscale="Viridis",style="linear"):
        '''
        Return a plotly figure object for the available plot types set in set_plot()..
        colorscale = Viridis [default], Cividis, BlueRed or Rainbow
        '''
        #Set some text strings
        txtbot="Model: TIEGCM,  Run: " + str(self.runname) + ",  " + str(self.gridSize) + " cells,  minimum dx=" + str(self.gridMinDx)
        if style == "linear":
            txtbar=var + " [" + self.variables[var]['units'] + "]"
        if style == "log":
            txtbar="log("+var + ") [" + self.variables[var]['units'] + "]"
        print("before interpolation")
        # Get values from interpolation already computed
        tic = time.perf_counter()
        if (self.newgrid[:,1]).max() < 50: # IP or height is column 1 (CTIPe_singletime: 0)
#            result=self[var](self.newgrid[:,0].ravel(),self.newgrid[:,1].ravel(),self.newgrid[:,2].ravel(),self.newgrid[:,3].ravel())
            result=self[var](self.newt,self.newx,self.newy,self.newz)
        else:
            self.varname=var
# once we have a gridify'd version we can use this
#            result=vert_interp(self.newt,self.newx,self.newy,self.newz)
# have to use np.vectorize'd version
            result=self.vert_interp(self,self.newgrid[:,0].ravel(),self.newgrid[:,1].ravel(),self.newgrid[:,2].ravel(),self.newgrid[:,3].ravel())
        toc = time.perf_counter()
        print(f"Time for computing interpolations: {toc - tic:0.4f} seconds")
        print("result.shape after interpolation: ",result.shape)
        
        if self.plottype == "LonLat" or self.plottype=="LonLatH":
            txttop=self.plots[self.plottype]['cut']+" "+self.cutunit+"=" + str(self.plots[self.plottype]['cutV']) + " slice,  Time = " + self.filetime
            tint=self.newt
            xint = self.newz
            yint = self.newy
            zint = self.newx
#            print("tint:",tint)
#            print("xint:",xint)
#            print("yint:",yint)
#            print("zint:",zint)

            xunit=self.zunit
            yunit=self.yunit
            xrange=dict(min=xint.min(),max=xint.max(),n=len(xint))
            yrange=dict(min=yint.min(),max=yint.max(),n=len(yint))
            # Reshape interpolated values into 2D
            print ("result shape: ",result.shape)
            result2=np.reshape(result,(self.nY,self.nZ))
            print ("result2 shape: ",result2.shape)
#            print(result2)
            if style == "linear":
                def plot_XY(xint = xint, yint = yint):
                    return result2
            if style == "log":
                def plot_XY(xint = xint, yint = yint):
                    return np.log(result2)
            plotXY = Kamodo(plot_XY = plot_XY)
            fig = plotXY.plot(plot_XY = dict())
            fig.update_xaxes(title_text="",scaleanchor='y')
            fig.update_yaxes(title_text="Lat [deg]")
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
                fig.update_traces(colorscale="RdBu",
                                  reversescale=True,
                )
            elif colorscale == "Cividis":
                fig.update_traces(colorscale="Cividis")
            else:
                fig.update_traces(colorscale="Viridis")
            fig.update_traces(
#                zmin=cmin, zmax=cmax,
                ncontours=201,
                colorbar=dict(title=txtbar,tickformat='.4g'),
                    hovertemplate="X: %{x:.2f}<br>Y: %{y:.2f}<br><b> %{z:.4g}</b><extra></extra>",
                contours=dict(coloring="fill", showlines=False)
            )
            if xunit == yunit:
            # real aspect ratio
                if xrange['max']-xrange['min'] > yrange['max']-yrange['min']:
                    aspectratio=dict(x=np.asarray([4.,(xrange['max']-xrange['min'])/np.asarray([1.e-5,(yrange['max']-yrange['min'])]).max()]).min(),y=1)
                else:
                    aspectratio=dict(x=1,y=np.asarray([4,(yrange['max']-yrange['min'])/np.asarray([1.e-5,(xrange['max']-xrange['min'])]).max()]).min())
            else:
                aspectratio=dict(x=1.25,y=1)
            
            width=180+300*aspectratio['x']
            height=90+300*aspectratio['y']
            fig.update_layout(
                autosize=False,
                height=height,
                width=width,
                margin=dict(t=45,b=45,l=40,r=140),
                title=dict(text=txttop),
                annotations=[
                    dict(text="Lon [deg]", x=0.5, y=-0.10, showarrow=False, xref="paper", yref="paper", font=dict(size=14)),
                    dict(text=txtbot,
                         font=dict(size=16, family="sans serif", color="#000000"),
                         x=0.0, y=0.0, ax=0, ay=0, xanchor="left", xshift=-65, yshift=-42, xref="paper", yref="paper"
                    )
                ],
                xaxis = dict(
                    tickmode = 'linear',
                    tick0 = 0,
                    dtick = 30
                ),
                yaxis = dict(
                    tickmode = 'linear',
                    tick0 = 0,
                    dtick = 30
                )
            )
            return fig


        if self.plottype == "LonIP" or self.plottype == "LonH":
            xint = self.newz
            yint = self.newx/self.xfactor
            xunit=self.zunit
            yunit=self.xunit

            xrange=dict(min=xint.min(),max=xint.max(),n=len(xint))
            yrange=dict(min=yint.min(),max=yint.max(),n=len(yint))

            # Reshape interpolated values into 2D
            result2=np.reshape(result,(self.nX,self.nZ))
#            def plot_XZ(xint = xint, yint = yint):
#                return result2
            print(style)
            if style=="linear":
                def plot_XZ(xint = xint, yint = yint):
                    return result2
            if style=="log":
                def plot_XZ(xint = xint, yint = yint):
                    return np.log(result2)
            plotXZ = Kamodo(plot_XZ = plot_XZ)
            fig = plotXZ.plot(plot_XZ = dict())
#            fig.update_xaxes(title_text="",scaleanchor='y')
            fig.update_xaxes(title_text="")
            vert_coord=self.xname
            vert_unit=self.xunit
            vert_format=self.xformat
            txttop=self.cut+" "+self.cutunit+"=" + str(self.plots[self.plottype]['cutV']) + " slice,  Time = " + self.filetime
            fig.update_yaxes(title_text=vert_coord+" "+vert_unit,tickformat=".4g")
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
                fig.update_traces(colorscale="RdBu",
                                  reversescale=True,
                )
            elif colorscale == "Cividis":
                fig.update_traces(colorscale="Cividis")
            else:
                fig.update_traces(colorscale="Viridis")
            fig.update_traces(
#                zmin=cmin, zmax=cmax,
                ncontours=201,
                colorbar=dict(title=txtbar,tickformat='.4g'),
                hovertemplate="Lon: %{x:.2f}<br>"+vert_coord+": %{y:"+vert_format+"}<br><b> %{z:.4g}</b><extra></extra>",
#                hovertemplate="Lon: %{x:.2f}<br>"+vert_coord+": %{y:.2f}<br><b> %{z:.4g}</b><extra></extra>",
                contours=dict(coloring="fill", showlines=False)
            )
            if xunit == yunit:
            # real aspect ratio
                if xrange['max']-xrange['min'] > yrange['max']-yrange['min']:
                    aspectratio=dict(x=np.asarray([4.,(xrange['max']-xrange['min'])/np.asarray([1.e-5,(yrange['max']-yrange['min'])]).max()]).min(),y=1)
                else:
                    aspectratio=dict(x=1,y=np.asarray([4,(yrange['max']-yrange['min'])/np.asarray([1.e-5,(xrange['max']-xrange['min'])]).max()]).min())
            else:
                aspectratio=dict(x=1.25,y=1)
            
            width=180+300*aspectratio['x']
            height=90+300*aspectratio['y']
            fig.update_layout(
                autosize=False,
                height=height,
                width=width,
                margin=dict(t=45,b=45,l=40,r=140),
                title_text=txttop,
                annotations=[
                    dict(text="Lon [deg]", x=0.5, y=-0.10, showarrow=False, xref="paper", yref="paper", font=dict(size=14)),
                    dict(text=txtbot,
                         font=dict(size=16, family="sans serif", color="#000000"),
                         x=0.0, y=0.0, ax=0, ay=0, xanchor="left", xshift=-65, yshift=-42, xref="paper", yref="paper"
                    )
                ]
            )
            return fig;
            
        if self.plottype == "LatIP" or self.plottype == "LatH":
            xint = self.newy
            yint = self.newx/self.xfactor
            txttop=self.cut+" "+self.cutunit+"=" + str(self.plots[self.plottype]['cutV']) + " slice,  Time = " + self.filetime

            xunit=self.yunit
            yunit=self.xunit
            vert_coord=self.xname
            vert_format=self.xformat
            xrange=dict(min=xint.min(),max=xint.max(),n=len(xint))
            yrange=dict(min=yint.min(),max=yint.max(),n=len(yint))
            # Reshape interpolated values into 2D
#            result2=np.transpose(np.reshape(result,(self.nY,self.nX)))
            result2=np.reshape(result,(self.nX,self.nY))
            if style=="linear":
                def plot_XY(xint = xint, yint = yint):
                    return result2
            if style=="log":
                def plot_XY(xint = xint, yint = yint):
                    return np.log(result2)
            plotXY = Kamodo(plot_XY = plot_XY)
            fig = plotXY.plot(plot_XY = dict())
#            fig.update_xaxes(title_text="",scaleanchor='y')
            fig.update_xaxes(title_text="")
            fig.update_yaxes(title_text=vert_coord+" "+yunit,tickformat=".4g")
            # Choose colorscale
            if colorscale == "Rainbow":
                fig.update_traces(
                    colorscale=[[0.00, 'rgb(0,0,255)'],
                                [0.25, 'rgb(0,255,255)'],
                                [0.50, 'rgb(0,255,0)'],
                                [0.75, 'rgb(255,255,0)'],
                                [1.00, 'rgb(255,0,0)']]
                )
            elif colorscale == "Cividis":
                fig.update_traces(colorscale="Cividis")
            elif colorscale == "BlueRed":
                fig.update_traces(colorscale="RdBu",
                                  reversescale=True,
                )
            else:
                fig.update_traces(colorscale="Viridis")
            fig.update_traces(
#                zmin=cmin, zmax=cmax,
                ncontours=201,
                colorbar=dict(title=txtbar,tickformat='.4g'),
                hovertemplate="Lat: %{x:.2f}<br>"+vert_coord+": %{y:"+vert_format+"}<br><b> %{z:.4g}</b><extra></extra>",
#                hovertemplate="Lat: %{x:.2f}<br>"+vert_coord+": %{y:.2f}<br><b> %{z:.4g}</b><extra></extra>",
                contours=dict(coloring="fill", showlines=False)
            )
            if xunit == yunit:
            # real aspect ratio
                if xrange['max']-xrange['min'] > yrange['max']-yrange['min']:
                    aspectratio=dict(x=np.asarray([4.,(xrange['max']-xrange['min'])/np.asarray([1.e-5,(yrange['max']-yrange['min'])]).max()]).min(),y=1)
                    if aspectratio['x'] > 2:
                        aspectratio_x=aspectratio['x']
                        aspectratio=dict(x=2., y=aspectratio['y']*2./aspectratio_x)
                else:
                    aspectratio=dict(x=1,y=np.asarray([4,(yrange['max']-yrange['min'])/np.asarray([1.e-5,(xrange['max']-xrange['min'])]).max()]).min())
            else:
                aspectratio=dict(x=1.25,y=1)

            width=180+300*aspectratio['x']
            height=90+300*aspectratio['y']
            fig.update_layout(
                autosize=False,
                height=height,
                width=width,
                margin=dict(t=45,b=45,l=40,r=140),
                title_text=txttop,
                annotations=[
                    dict(text="Lat [deg]", x=0.5, y=0.0, showarrow=False, xref="paper", yref="paper", yshift=-32,font=dict(size=14)),
                    dict(text=txtbot,
                         font=dict(size=16, family="sans serif", color="#000000"),
                         x=0.0, y=0.0, ax=0, ay=0, xanchor="left", xshift=-45, yshift=-42, xref="paper", yref="paper"
                    )
                ]
            )
            return fig


        print("ERROR, unknown plottype =",plottype,", exiting without definition of figure.")
        return
        
    def ellipse_arc(self, x_center=0, y_center=0, a=1, b=1, start_angle=-np.pi/2, end_angle=np.pi/2, N=100, closed= False):
        '''small function to overlay a partial ellipse on the plot for sun-lit half of earth effect'''
        t = np.linspace(start_angle, end_angle, N)
        x = x_center + a*np.cos(t)
        y = y_center + b*np.sin(t)
        path = f'M {x[0]}, {y[0]}'
        for k in range(1, len(t)):
            path += f'L{x[k]}, {y[k]}'
        if closed:
            path += ' Z'
        return path
    




