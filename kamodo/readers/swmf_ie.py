from kamodo import Kamodo, kamodofy
#, gridify
from netCDF4 import Dataset
import numpy as np
import os
import scipy
from scipy import interpolate
from scipy.interpolate import RegularGridInterpolator
import math
import scipy.constants as constants

from plotly.offline import init_notebook_mode, iplot, plot
from plotly.subplots import make_subplots

import time
#from util import time_in_interval
from util import *
#from util import boundary_conditions, fill_masked

class SWMF_IE(Kamodo):
    def __init__(self,filename,**kwargs):
        print('opening SWMF Ionosphere Electrodynamics file %s' % filename)
        self.filename = filename
        self.missing_value=np.NAN
        data=self.read_swmf_ie(filename)
        self.variables=dict()
        self.coords=['theta','phi','altitude']
        self.coord_units=['deg','deg','km']
        self.x=data['Theta']
        self.y=data['Psi']
        self.z=[110.]
#        self.coordinate_system='SM'
# these three variables are needed to prevent errors in Kamodo      
        self.verbose=False;
        self.symbol_registry=dict();
        self.signatures=dict();
       
        #        self.z_, self.y_, self.x_ = scipy.meshgrid(self.z, self.y, self.x, indexing = 'ij')
        self.y_, self.x_ = scipy.meshgrid(self.y, self.x, indexing = 'ij')

#        print(data['orig_units'])
#        print(data['units'])
        for ivar in range(len(data['variables'])):
            varname=data['variables'][ivar]
            var_data=data['data'][:,:,ivar].squeeze().T
            unit=data['units'][ivar]
#            print("variable:",varname,"unit:",unit,"data:",var_data)

            self.variables[varname]=dict(units=unit, data=var_data)
            self.register_variable(varname,unit)
            
    def read_swmf_ie(self,filename):
        import re
        import datetime
        arrays = []
        orig_units=[]
        units=[]
        variables=[]
        iline=0
        with open(filename, 'r') as a:
            for line in a.readlines():
                A = re.match(r'TITLE=(.*$)', line, re.M | re.I)
                B = re.match(r'(\s*)VARIABLES=(.*$)', line, re.M | re.I)
                C = re.match(r'(\s*)(\")(.*$)', line,re.M | re.I)
                D = re.match(r'(\s*)ZONE (.*$)', line, re.M | re.I)
                E = re.match(r'(\s*)I=(.*$)', line, re.M | re.I)
                if A or B or C or D or E:
                    if A:
                        (title,date,tilt)=line.split(',')
                        date.replace(" ","")
                        (year,month,day,hour,minute,second,ms)=date.split('-')
                        self.Date=datetime.datetime(int(year),int(month),int(day),hour=int(hour),minute=int(minute),second=int(second),microsecond=int(ms)*1000,tzinfo=datetime.timezone.utc)
                    if B or C:
                        for s in (line.split('"'))[1:]:
                            if s != "," and s != '':
# any strings containing variable name and unit
                                if s != "\n":
                                    if s.find('[') > 0:
                                        (var,unit) = s.split("[")
                                        unit="["+unit
                                    else:
                                        var=s
                                        unit="[]"
# cannot have " " or "-" sign in variable names
                                    if var.find('conjugate') == 0:
                                        var=var[10:]+'_conj'
                                    var=var.replace(" ","")
                                    var=var.replace("/","over")
                                    var=var.replace("-","")
                                    variables.append(var)
                            # map unit names to something SymPy may understand
                                    if len(unit) >= 2:
# strip off "[" and "]"
                                        unit=unit[1:-1]
                                        orig_unit=unit
# special rules: R -> R_E `m -> mu /(...) -> 1/(...), m2 or m^2 -> m**2
                                    if unit == 'R':
                                        unit='R_E'
                                    unit=re.sub(r"^`m","mu",unit)
                                    unit=re.sub(r"^/","1/",unit)
                                    unit=unit.replace('m2','m**2')
                                    unit=unit.replace('m^2','m**2')
                                    orig_units.append(orig_unit)
                                    units.append(unit)
                    if E:
                        index_i=line.index("I=")+2
                        index_j=line.index("J=")+2
                        NI=int(line[index_i:].split()[0])
                        NJ=int(line[index_j:].split()[0])
                        continue
                else:
                    for s in line.split():
                        arrays.append(float(s))
        print(variables)
        print(units)
                        
        nvar=len(variables)
        nelements=len(arrays)
        npos=int(nelements/nvar)           
        arrays = np.array(arrays)
    
        arrays=arrays.reshape((npos,nvar))
        arrays_N=arrays[0:int(npos/2),:].reshape((NJ,NI,nvar))
        arrays_S=arrays[int(npos/2):,:].reshape((NJ,NI,nvar))
        data=np.concatenate((arrays_N,arrays_S[:,1:,:]),axis=1)
        data[:,0,3]=0.
        data[:,-1,3]=180.
        df={'data':data,
            'variables':variables,
            'orig_units':orig_units,
            'units':units,
            'axis1':'Theta',
            'axis2':'Psi',
            'Theta': data[0,:,3].flatten(),
            'Psi': data[:,0,4].flatten()
        }
    
        return df

    def get_grid_interpolator(self, varname):
        """create a regular grid interpolator for this variable"""
        data =  self.variables[varname]['data']
        interpolator = RegularGridInterpolator((self.x, self.y),
                                               data, 
                                               bounds_error = False,
                                               fill_value = self.missing_value)
        return interpolator
        
    def register_variable(self, varname, units):
        interpolator = self.get_grid_interpolator(varname)
        
        # store the interpolator
        self.variables[varname]['interpolator'] = interpolator

        def interpolate(xvec):  
            return self.variables[varname]['interpolator'](xvec)

        # update docstring for this variable
        interpolate.__doc__ = "A function that returns {} in [{}].".format(varname,units)

        self[varname] = kamodofy(interpolate, 
                                 units = units, 
                                 citation = "Pembroke et al 2019, Rastaetter 2020",
                                 data = None)
        self[varname + '_ij'] = kamodofy(gridify(self[varname], 
                                                  x_i = self.x, 
                                                  y_j = self.y),
                                          units = units,
                                          citation = "Pembroke et al 2019, Rastaetter 2020",
                                          data = self.variables[varname]['data'])
        return

    def get_plot(self, var, plottype, runname, colorscale="BlueRed", sym="T"):
        '''
        Return a plotly figure object for the plottype requested)..
        var, plottype, and runname are required variables. 
        colorscale = BlueRed [default], Viridis, Cividis, or Rainbow
        sym = T [default] for symetric colorscale around 0
        '''

        if plottype.count('2D-') > 0:
            x = np.linspace(-.65, .65, 130)
            y = np.linspace(-.65, .65, 130)
            xx, yy = np.meshgrid(np.array(x), np.array(y))
            if plottype == "2D-N":
                loctxt = "Northern"
                llat = 180.-np.arcsin(np.sqrt(xx*xx + yy*yy))*180./np.pi
            if plottype == "2D-S":
                loctxt = "Southern"
                llat = np.arcsin(np.sqrt(xx*xx + yy*yy))*180./np.pi
            llon = 180.-(90.+np.arctan(xx/yy)*180./np.pi)
            llon[yy < 0.] += 180.
            grid = np.ndarray(shape=(np.size(np.reshape(xx,-1)),2), dtype=np.float32)
            grid[:,0] = np.reshape(llat,-1)
            grid[:,1] = np.reshape(llon,-1)
            units=self.variables[var]['units']
            test = self.variables[var]['interpolator'](grid)
            result = np.reshape(test,(y.shape[0],x.shape[0]))
            if sym == "T":
                cmax = np.max(np.absolute(self.variables[var]['data']))
                #cmax = np.max(np.absolute(result))
                cmin = -cmax
            else:
                cmax = np.max(self.variables[var]['data'])
                #cmax = np.max(result)
                cmin = np.min(self.variables[var]['data'])
            
            time=self.Date.strftime("%Y/%m/%d %H:%M:%S UT")
            
            def plot_var(y = y, x = x):
                return result
            plotvar = Kamodo(plot_var = plot_var)
            
            fig = plotvar.plot(plot_var = dict())
            fig.update_xaxes(nticks=7,title_text="",scaleanchor='y',autorange="reversed")
            fig.update_yaxes(nticks=7,title_text="")
            txtbar = var + " [" + units + "]"
            txtbot = "Model: SWMF-IE,  Run: " + runname
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
                colorbar=dict(title=txtbar),
                hovertemplate=loctxt + " Hemisphere<br>X: %{y:.2f}<br>Y: %{x:.2f}<br><b>"+var+": %{z:.4g}</b><extra></extra>",
                contours=dict(coloring="fill", showlines=False)
            )
            c80=np.sin(((90.-80.)/90.)*np.pi/2.)
            c70=np.sin(((90.-70.)/90.)*np.pi/2.)
            c60=np.sin(((90.-60.)/90.)*np.pi/2.)
            c50=np.sin(((90.-50.)/90.)*np.pi/2.)
            c50d=c50/np.sqrt(2.)
            fig.update_layout(
                title=dict(text=loctxt + " Hemisphere,  Time = " + time,
                           yref="container", yanchor="top", y=0.95),
                title_font_size=16,
                shapes=[
                    dict(type="circle",
                         xref="x", yref="y", x0=-c80, y0=-c80, x1=c80, y1=c80,
                         line=dict(color="black", width=1, dash="dash")),
                    dict(type="circle",
                         xref="x", yref="y", x0=-c70, y0=-c70, x1=c70, y1=c70,
                         line=dict(color="black", width=1)),
                    dict(type="circle",
                         xref="x", yref="y", x0=-c60, y0=-c60, x1=c60, y1=c60,
                         line=dict(color="black", width=1, dash="dash")),
                    dict(type="circle",
                         xref="x", yref="y", x0=-c50, y0=-c50, x1=c50, y1=c50,
                         line=dict(color="black", width=1)),
                    dict(type="line",
                         xref="x", yref="y", x0=-c50, y0=0., x1=c50, y1=0.,
                         line=dict(color="black", width=1)),
                    dict(type="line",
                         xref="x", yref="y", x0=0., y0=-c50, x1=0., y1=c50,
                         line=dict(color="black", width=1)),
                    dict(type="line",
                         xref="x", yref="y", x0=-c50d, y0=-c50d, x1=c50d, y1=c50d,
                         line=dict(color="black", width=1)),
                    dict(type="line",
                         xref="x", yref="y", x0=-c50d, y0=c50d, x1=c50d, y1=-c50d,
                         line=dict(color="black", width=1))
                ],
                annotations=[
                    dict(text="Y [RE]", x=0.5, y=-0.11, showarrow=False, 
                         xref="paper", yref="paper", font=dict(size=12)),
                    dict(text="X [RE]", x=-0.19, y=0.5, showarrow=False, 
                         xref="paper", yref="paper", font=dict(size=12), textangle=-90),
                    dict(text="midnight", x=0.5, y=0.0, showarrow=False, 
                         xref="paper", yref="paper", font=dict(size=10)),
                    dict(text="noon", x=0.5, y=1.0, showarrow=False, 
                         xref="paper", yref="paper", font=dict(size=10)),
                    dict(text="dawn", x=1.0, y=0.5, showarrow=False, textangle=90,
                         xref="paper", yref="paper", font=dict(size=10)),
                    dict(text="dusk", x=0.0, y=0.5, showarrow=False, textangle=-90,
                         xref="paper", yref="paper", font=dict(size=10)),
                    dict(text=txtbot, x=0.0, y=0.0, ax=0, ay=0, xanchor="left",
                         xshift=-65, yshift=-42, xref="paper", yref="paper",
                         font=dict(size=16, family="sans serif", color="#000000")
                    )
                ],
                height=375,
                width=500,
                margin=dict(t=45,r=140)
            )
            return fig


        if self.plottype == "3D":
            return

        print('Error, no valid plottype was given.')
        return

    def list_variables(self):
        '''
        Return an array of the variables that can be interpolated/plotted.
        '''
        vars=[k for k in self.variables.keys()]
        return vars
