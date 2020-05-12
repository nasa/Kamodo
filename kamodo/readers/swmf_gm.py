import spacepy
import spacepy.pybats

import plotly
import plotly.graph_objs as go
from plotly.offline import init_notebook_mode, iplot, plot
init_notebook_mode(connected = True)

import numpy as np

import scipy.interpolate as spint
import scipy.spatial.qhull as qhull

from kamodo import Kamodo, kamodofy, gridify
import time
import glob

from os import path
from datetime import datetime, timedelta


class SWMF_GM(Kamodo): 
    def __init__(self, 
                 filename, 
                 runpath = "./", 
                 runname = "noname", 
                 start_time = "1900/01/01 00:00", 
                 **kwargs):
        # Start timer
        tic = time.perf_counter()

        # perform any necessary I/O
        print('opening {}'.format(filename))
        self.filename = filename
        self.runpath = runpath
        self.runname = runname
        self.start_time = start_time
        self.missing_value = np.NAN
        self.variables=dict()
        mhd = spacepy.pybats.IdlFile(filename)

        # Pull out variables and units
        vars = list(mhd.keys())
        vars = vars[4:]
        self.nvars = len(vars)
        self.vars = vars
        units = mhd.meta['header'].split("_", 1)[0].strip().split(" ")
        units = units[3:]
        # var = vars[0]
        # unit = units[vars.index(var)]
        # self.nvars = mhd.meta['nvar']
        print('... found',self.nvars,'variables')
        self.var = vars[0]
        for var in vars:
            self.variables[var] = dict(units = units[vars.index(var)], data = np.array(mhd[var], dtype=np.float32) )
        
        # Pull start_time from the DatabaseInfo file (if it exists), then add simulated time from file metadata
        import os
        from os import path
        from datetime import datetime, timedelta
        if path.isfile(runpath+runname+'/DatabaseInfo'):
            for line in open(runpath+runname+'/DatabaseInfo'):
                if "start_time" in line:
                    self.start_time = line.split("#", 1)[0].strip()
                    print('... found run start_time = ',self.start_time)
                    break

        dt = datetime.strptime(self.start_time, '%Y/%m/%d %H:%M')
        dt = dt + timedelta(seconds=int(mhd.meta['runtime']))
        self.filetime = dt.strftime("%Y/%m/%d %H:%M:%S UT")
        print('... setting filetime = ',self.filetime)
        
        # Pull grid to interpolate from into stored arrays
        self.grid = np.empty([mhd['x'].size, 3], dtype=np.float32)
        self.grid[:,0] = np.array(mhd['x'])
        self.grid[:,1] = np.array(mhd['y'])
        self.grid[:,2] = np.array(mhd['z'])
        dd = np.absolute(np.concatenate(( np.array(mhd['x']), np.array([999.], dtype=np.float32) )) - 
                         np.concatenate(( np.array([999.], dtype=np.float32), np.array(mhd['x']) )) )
        self.gridSize = mhd['x'].size
        self.gridMinDx = np.amin(dd[dd[:]>0.])
        print('... file grid has',self.gridSize,'cells with minimum dx =',self.gridMinDx)
        
        # Prepare model for function registration
        super(SWMF_GM, self).__init__(**kwargs) 
        
        # Interpolation values
        self.plottype = "XY"
        self.cutV = 0.
        self.nX = 141
        self.nY = 81
        self.nZ = 1
        self.newx = np.linspace(-50., 20., self.nX)
        self.newy = np.linspace(-20., 20., self.nY)
        self.newz = self.cutV
        self.tol = 1.1
        
        # Get grids ready to use
        self.setup_interpolating_grids()
        
        for varname in self.variables:
            units = self.variables[varname]['units']
            self.register_variable(varname, units)

        # end timer
        toc = time.perf_counter()
        print(f"Time loading file and precomputing default interpolations: {toc - tic:0.4f} seconds")

    def setup_interpolating_grids(self):
        '''setup the grid to interpolate to, trim to necessary size of source grid, and compute interpolation weights'''
        # Setup grid to interpolate to
        xx, yy, zz = np.meshgrid(self.newx,self.newy,self.newz, indexing = 'xy')
        self.newgrid = np.ndarray(shape=(np.size(np.reshape(xx,-1)),3), dtype=np.float32)
        self.newgrid[:,0] = np.reshape(xx,-1)
        self.newgrid[:,1] = np.reshape(yy,-1)
        self.newgrid[:,2] = np.reshape(zz,-1)
        
        # Reduce size of read data block to speed up the interpolation
        self.grid2 = self.grid[(self.grid[:,0] > (np.amin(self.newx)-self.tol)) & 
                               (self.grid[:,0] < (np.amax(self.newx)+self.tol)) &
                               (self.grid[:,1] > (np.amin(self.newy)-self.tol)) & 
                               (self.grid[:,1] < (np.amax(self.newy)+self.tol)) & 
                               (self.grid[:,2] > (np.amin(self.newz)-self.tol)) & 
                               (self.grid[:,2] < (np.amax(self.newz)+self.tol))]
        
        # Precompute weights for interpolation to optimize
        self.vtx, self.wts = self.GMcompute_weights(self.grid2, self.newgrid)

        return

    def register_variable(self, varname, units):
        """register variables into Kamodo for this model, SWMF-GM"""
        interpolator = self.get_grid_interpolator(varname)
        
        # store the interpolator
        self.variables[varname]['interpolator'] = interpolator

        def interpolate(xvec):  
            # Reduce size of read data block to speed up the interpolation
            tmpgrid2 = self.grid[(self.grid[:,0] > (np.amin(xvec[:,0])-self.tol)) & 
                                 (self.grid[:,0] < (np.amax(xvec[:,0])+self.tol)) &
                                 (self.grid[:,1] > (np.amin(xvec[:,1])-self.tol)) & 
                                 (self.grid[:,1] < (np.amax(xvec[:,1])+self.tol)) & 
                                 (self.grid[:,2] > (np.amin(xvec[:,2])-self.tol)) & 
                                 (self.grid[:,2] < (np.amax(xvec[:,2])+self.tol))]
            data =  self.variables[varname]['data']
            tmpdata = data[(self.grid[:,0] > (np.amin(xvec[:,0])-self.tol)) & 
                           (self.grid[:,0] < (np.amax(xvec[:,0])+self.tol)) &
                           (self.grid[:,1] > (np.amin(xvec[:,1])-self.tol)) & 
                           (self.grid[:,1] < (np.amax(xvec[:,1])+self.tol)) & 
                           (self.grid[:,2] > (np.amin(xvec[:,2])-self.tol)) & 
                           (self.grid[:,2] < (np.amax(xvec[:,2])+self.tol))]
            # Precompute weights for interpolation to optimize
            tmpvtx, tmpwts = self.GMcompute_weights(tmpgrid2, xvec)
            return self.GMinterpolate(tmpdata, tmpvtx, tmpwts)

        # update docstring for this variable
        interpolate.__doc__ = "A function that returns {} in [{}].".format(varname,units)

        self[varname] = kamodofy(interpolate, 
                                 units = units, 
                                 citation = "De Zeeuw 2020",
                                 data = None)
        #self[varname + '_ijk'] = kamodofy(gridify(self[varname], 
        #                                          x_i = self.newgrid[:,0], 
        #                                          y_j = self.newgrid[:,1], 
        #                                          z_k = self.newgrid[:,2]),
        #                                  units = units,
        #                                  citation = "De Zeeuw 2020",
        #                                  data = self.variables[varname]['data'])

    def get_grid_interpolator(self, varname):
        """create a grid interpolator for passed variable based on Delaunay triangulation of points"""
        data =  self.variables[varname]['data']
        data2 = data[(self.grid[:,0] > (np.amin(self.newx)-self.tol)) & 
                     (self.grid[:,0] < (np.amax(self.newx)+self.tol)) &
                     (self.grid[:,1] > (np.amin(self.newy)-self.tol)) & 
                     (self.grid[:,1] < (np.amax(self.newy)+self.tol)) & 
                     (self.grid[:,2] > (np.amin(self.newz)-self.tol)) & 
                     (self.grid[:,2] < (np.amax(self.newz)+self.tol))]
        interpolator = self.GMinterpolate(data2, self.vtx, self.wts)
        # interpolator = spint.griddata(self.grid2, data2, self.newgrid, method='linear')

        return interpolator
    
    def reinterpolate_values(self):
        '''recompute interpolation from saved grid and weights and store values for plotting'''
        for varname in self.variables:
            data =  self.variables[varname]['data']
            data2 = data[(self.grid[:,0] > (np.amin(self.newx)-self.tol)) & 
                         (self.grid[:,0] < (np.amax(self.newx)+self.tol)) &
                         (self.grid[:,1] > (np.amin(self.newy)-self.tol)) & 
                         (self.grid[:,1] < (np.amax(self.newy)+self.tol)) & 
                         (self.grid[:,2] > (np.amin(self.newz)-self.tol)) & 
                         (self.grid[:,2] < (np.amax(self.newz)+self.tol))]
            self.variables[varname]['interpolator'] = self.GMinterpolate(data2, self.vtx, self.wts)

        return

    def GMcompute_weights(self, xyz, uvw):
        """compute the weights to optimize griddata interpolation using a constructed Delaunay grid"""
        d = 3
        tri = qhull.Delaunay(xyz)
        simplex = tri.find_simplex(uvw)
        vertices = np.take(tri.simplices, simplex, axis=0)
        temp = np.take(tri.transform, simplex, axis=0)
        delta = uvw - temp[:, d]
        bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
        return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))

    def GMinterpolate(self, values, vtx, wts):
        """
        Use precomputed weights to complete the interpolation.
        Note: Use of a Delaunay triangulation with a cartesian grid will result in any given point
              having the possibility of 8 different point becoming the 4 used points for the tetrahedal.
              Thus, the interpolation can randomly vary the interpolated value with any changes to the
              number of points in the original grid.
        """
        return np.einsum('nj,nj->n', np.take(values, vtx), wts)
    
    def chk_plot(self):
        '''Check current values of plotting variables and print them out.'''
        print('Values set for plotting:')
        print(' plottype = ',self.plottype)
        print(' cutV = ',self.cutV)
        print(' nX',self.nX)
        print(' nY',self.nY)
        print(' nZ',self.nZ)
        print(' newx',self.newx)
        print(' newy',self.newy)
        print(' newz',self.newz)
        print(' tol',self.tol)
        print(' full grid shape: ',self.grid.shape)
        print(' reduced size grid shape: ',self.grid2.shape)
        print(' interpolation grid shape: ',self.newgrid.shape)
        return
    
    def set_plot(self, plottype = "XY", cutV = 0.):
        '''Set plotting variables for available preset plot types: XY, YZ, YZ, XYZ'''
        tic = time.perf_counter()
        if plottype == self.plottype:
            if cutV == self.cutV:
                print('Plottype is unchanged, returning. ',plottype)
                return
        
        if plottype == "XY":
            self.plottype = plottype
            self.cutV = cutV
            self.nX = 141
            self.nY = 81
            self.nZ = 1
            self.newx = np.linspace(-50., 20., self.nX)
            self.newy = np.linspace(-20., 20., self.nY)
            self.newz = cutV
            self.tol = 1.1
        elif plottype == "XZ":
            self.plottype = plottype
            self.cutV = cutV
            self.nX = 141
            self.nY = 1
            self.nZ = 81
            self.newx = np.linspace(-50., 20., self.nX)
            self.newy = cutV
            self.newz = np.linspace(-20., 20., self.nZ)
            self.tol = 1.1
        elif plottype == "YZ":
            self.plottype = plottype
            self.cutV = cutV
            self.nX = 1
            self.nY = 141
            self.nZ = 81
            self.newx = cutV
            self.newy = np.linspace(-35., 35., self.nY)
            self.newz = np.linspace(-20., 20., self.nZ)
            self.tol = 1.1
        elif plottype == "XYZ":
            self.plottype = plottype
            return
        else:
            print('Error, unknown plottype. ',plottype)
            return

        self.setup_interpolating_grids()
        self.reinterpolate_values()

        toc = time.perf_counter()
        print(f"Time resetting plot and precomputing interpolations: {toc - tic:0.4f} seconds")
        return
    
    def get_plot(self, var, colorscale="Viridis"):
        '''
        Return a plotly figure object for the available plot types set in set_plot()..
        colorscale = Viridis [default], Cividis, or Rainbow
        '''
        #Set some text strings
        txtbot="Model: BATSRUS,  Run: " + str(self.runname) + ",  " + str(self.gridSize) + " cells,  minimum dx=" + str(self.gridMinDx)
        txtbar=var + " [" + self.variables[var]['units'] + "]"
        
        # Get values from interpolation already computed
        result=self.variables[var]['interpolator']
        r = np.sqrt(np.square(self.newgrid[:,0]) + np.square(self.newgrid[:,1]) + np.square(self.newgrid[:,2]))
        cmin=np.amin(result[(r[:] > 2.999)])
        cmax=np.amax(result[(r[:] > 2.999)])
        
        if self.plottype == "XY":
            txttop="Z=" + str(self.newz) + " slice,  Time = " + self.filetime
            xint = self.newx
            yint = self.newy
            # Reshape interpolated values into 2D
            result2=np.reshape(result,(self.nY,self.nX))
            def plot_XY(xint = xint, yint = yint):
                return result2
            plotXY = Kamodo(plot_XY = plot_XY)
            fig = plotXY.plot(plot_XY = dict())
            fig.update_xaxes(title_text="",scaleanchor='y')
            fig.update_yaxes(title_text="Y [RE]")
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
            else:
                fig.update_traces(colorscale="Viridis")
            fig.update_traces(
                zmin=cmin, zmax=cmax,
                ncontours=201,
                colorbar=dict(title=txtbar),
                hovertemplate="X: %{x:.2f}<br>Y: %{y:.2f}<br><b> %{z:.4f}</b><extra></extra>",
                contours=dict(coloring="fill", showlines=False)
            )
            fig.update_layout(
                title_text=txttop,
                annotations=[
                    dict(text="X [RE]", x=0.5, y=-0.10, showarrow=False, xref="paper", yref="paper", font=dict(size=14)),
                    dict(text=txtbot,
                         font=dict(size=16, family="sans serif", color="#000000"),
                         x=0.0, y=0.0, ax=0, ay=0, xanchor="left", xshift=-65, yshift=-45, xref="paper", yref="paper"
                    )
                ]
            )
            if self.newz == 0.:
                fig.update_layout(
                    shapes=[
                        dict(type="circle", xref="x", yref="y", x0=-3, y0=-3, x1=3, y1=3, fillcolor="black", line_color="black"),
                        dict(type="circle", xref="x", yref="y", x0=-1, y0=-1, x1=1, y1=1, fillcolor="black", line_color="white"),
                        dict(type="path", path= self.ellipse_arc(N=30), fillcolor="white", line_color="white")
                    ]
                )
            return fig
        
        if self.plottype == "XZ":
            txttop="Y=" + str(self.newy) + " slice,  Time = " + self.filetime
            xint = self.newx
            zint = self.newz
            # Reshape interpolated values into 2D
            result2=np.reshape(result,(self.nX,self.nZ))
            def plot_XZ(xint = xint, zint = zint):
                return result2
            plotXZ = Kamodo(plot_XZ = plot_XZ)
            fig = plotXZ.plot(plot_XZ = dict())
            fig.update_xaxes(title_text="",scaleanchor='y')
            fig.update_yaxes(title_text="z [RE]")
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
            else:
                fig.update_traces(colorscale="Viridis")
            fig.update_traces(
                zmin=cmin, zmax=cmax,
                ncontours=201,
                colorbar=dict(title=txtbar),
                hovertemplate="X: %{x:.2f}<br>Z: %{y:.2f}<br><b> %{z:.4f}</b><extra></extra>",
                contours=dict(coloring="fill", showlines=False)
            )
            fig.update_layout(
                title_text=txttop,
                annotations=[
                    dict(text="X [RE]", x=0.5, y=-0.10, showarrow=False, xref="paper", yref="paper", font=dict(size=14)),
                    dict(text=txtbot,
                         font=dict(size=16, family="sans serif", color="#000000"),
                         x=0.0, y=0.0, ax=0, ay=0, xanchor="left", xshift=-65, yshift=-45, xref="paper", yref="paper"
                    )
                ]
            )
            if self.newy == 0.:
                fig.update_layout(
                    shapes=[
                        dict(type="circle", xref="x", yref="y", x0=-3, y0=-3, x1=3, y1=3, fillcolor="black", line_color="black"),
                        dict(type="circle", xref="x", yref="y", x0=-1, y0=-1, x1=1, y1=1, fillcolor="black", line_color="white"),
                        dict(type="path", path= self.ellipse_arc(N=30), fillcolor="white", line_color="white")
                    ]
            )
            return fig
        
        if self.plottype == "YZ":
            txttop="X=" + str(self.newx) + " slice,  Time = " + self.filetime
            yint = self.newy
            zint = self.newz
            # Reshape interpolated values into 2D
            result2=np.reshape(result,(self.nY,self.nZ))
            def plot_YZ(yint = yint, zint = zint):
                return result2
            plotYZ = Kamodo(plot_YZ = plot_YZ)
            fig = plotYZ.plot(plot_YZ = dict())
            fig.update_xaxes(title_text="",scaleanchor='y')
            fig.update_yaxes(title_text="Z [RE]")
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
            else:
                fig.update_traces(colorscale="Viridis")
            fig.update_traces(
                zmin=cmin, zmax=cmax,
                ncontours=201,
                colorbar=dict(title=txtbar),
                hovertemplate="Y: %{x:.2f}<br>Z: %{y:.2f}<br><b> %{z:.4f}</b><extra></extra>",
                contours=dict(coloring="fill", showlines=False)
            )
            fig.update_layout(
                title_text=txttop,
                annotations=[
                    dict(text="Y [RE]", x=0.5, y=-0.10, showarrow=False, xref="paper", yref="paper", font=dict(size=14)),
                    dict(text=txtbot,
                         font=dict(size=16, family="sans serif", color="#000000"),
                         x=0.0, y=0.0, ax=0, ay=0, xanchor="left", xshift=-65, yshift=-45, xref="paper", yref="paper"
                    )
                ]
            )
            if self.newx == 0.:
                fig.update_layout(
                    shapes=[
                        dict(type="circle", xref="x", yref="y", x0=-3, y0=-3, x1=3, y1=3, fillcolor="black", line_color="black"),
                        dict(type="circle", xref="x", yref="y", x0=-1, y0=-1, x1=1, y1=1, fillcolor="white", line_color="white"),
                        #dict(type="path", path= self.ellipse_arc(N=30), fillcolor="white", line_color="white")
                    ]
            )
            return fig
        
        if self.plottype == "XYZ":
            txttop="3D Plot,  Time = " + self.filetime
            xint = self.newx
            yint = self.newy
            zint = self.newz
            # Reshape interpolated values into 3D
            result2=np.reshape(result,(self.nY,self.nX,self.nZ))
            def plot_XYZ(xint = xint, yint = yint, zint = zint):
                return result2
            plotXYZ = Kamodo(plot_XYZ = plot_XYZ)
            fig = plotXYZ.plot(plot_XYZ = dict())
            # Choose colorscale
            if colorscale == "Rainbow":
                fig.update_traces(
                    colorscale=[[0.00, 'rgb(0,0,255)'],
                                [0.25, 'rgb(0,255,255)'],
                                [0.50, 'rgb(0,255,0)'],
                                [0.75, 'rgb(255,255,0)'],
                                [1.00, 'rgb(255,0,0)']]
                )
            elif colorscale == "RdBu":
                fig.update_traces(colorscale="RdBu")
            elif colorscale == "Cividis":
                fig.update_traces(colorscale="Cividis")
            else:
                fig.update_traces(colorscale="Viridis")
            fig.update_traces(
                colorbar=dict(title=txtbar),
                cmin=cmin, cmax=cmax
            )
            fig.update_layout(
                scene_aspectmode='data',
                title_text=txttop,
                scene=dict(xaxis=dict(title="X [RE]"),yaxis=dict(title="Y [RE]"),zaxis=dict(title="Z [RE]")),
                annotations=[
                    dict(text=txtbot,
                         font=dict(size=16, family="sans serif", color="#000000"),
                         x=-0.08, y=-0.06, ax=0, ay=0, xanchor="left", xref="paper", yref="paper"
                    )
                ],
                margin=dict(l=0, b=30, t=30)
            )
            return fig
        
        print("ERROR, unknown plottype =",plottype,", exiting without definition of figure.")
        return
    
    # function to draw a partial ellipse (for sun-lit half of earth)
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
    
#
# Functions outside of the module that are useful.
#
def show_GM_files(runpath = "./", 
                  runname = "noname", 
                  start_time = "1000/01/01 00:00", 
                  filepath = "/IO2/3d*.out"):
    '''
    This function will show a list of available SWMF-GM plot files with simulated date/time for that file.
    The resulting files will be read from: runpath + runname + filepath (ie. ./noname/IO2/3d*.out)
    '''
    
    print('Routine prints simulation data/time for each simulation output and returns last file.',
         ' If year is \n1000, then no start_time for the run was found and arguement should be added to call.')
    
    # Pull start_time from the DatabaseInfo file (if it exists), then add simulated time from filenames
    if path.isfile(runpath+runname+'/DatabaseInfo'):
        for line in open(runpath+runname+'/DatabaseInfo'):
            if "start_time" in line:
                start_time = line.split("#", 1)[0].strip()
                print('... found run start_time = ',start_time)
                break
    print('')
    files = glob.glob(runpath+runname+filepath)
    files.sort()
    for file in files:
        runtime = file.split("_t")[-1].strip()
        runtime = str(runtime.split("_")[0].strip())
        hour=int(runtime[0:4])
        minute=int(runtime[4:6])
        second=int(runtime[6:])
        dt = datetime.strptime(start_time, '%Y/%m/%d %H:%M')
        dt = dt + timedelta(hours=hour, minutes=minute, seconds=second)
        filetime = dt.strftime("%Y/%m/%d %H:%M:%S UT")
        print(filetime,'   ',file)
    print('')
    return file

