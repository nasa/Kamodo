from kamodo import Kamodo, kamodofy
import numpy as np
import scipy
import datetime
from datetime import timezone
import urllib
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
from pandas import DatetimeIndex
from collections.abc import Iterable


def ror_get_extraction(server, runID, coord, satellite):
    '''Query for file contents from server'''
    query = '{}/{}/{}/{}_{}.txt'.format(server, runID, satellite, coord, satellite)
    response = urllib.request.urlopen(query)
    file = response.read()
    return file

class SATEXTRACT(Kamodo):
    def __init__(self, runID, coord, satellite, **kwargs):
        super(SATEXTRACT, self).__init__(**kwargs)
        self.verbose=False # overrides kwarg
        self.symbol_registry=dict() # wipes any user-defined symbols
        self.signatures=dict() # wipes any user-defined signatures
        self.RE=6.3781E3
        self.server = "https://ccmc.gsfc.nasa.gov/RoR_WWW/VMR/" # should be set by keyword
        self.runID = runID
        self.coordinates = coord
        self.satellite = satellite
        print(' -server: ',self.server)
        print(' -runID: ',runID)
        print(' -coordinate system: ',coord)
        print(' -satellite: ',satellite)
        self.variables=dict()
        self.file = ror_get_extraction(self.server, runID, coord, satellite).decode('ascii')
        self.parse_file()
        ts=self.tsarray[0]
        self.start = datetime.datetime.fromtimestamp(ts,tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        ts=self.tsarray[-1]
        self.stop  = datetime.datetime.fromtimestamp(ts,tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        print(" ")
        print(" -date start: ",self.start)
        print("         end: ",self.stop)

        for varname in self.variables:
            if varname == "N":
                continue
            units = self.variables[varname]['units']
            print('... registering ',varname,units)
            self.register_variable(varname, units)
        
        # classification of position into coordinates to assist visualizion
        self.possible_coords=('TOD','J2K','GEO','GM','GSM','GSE','SM')
        self.possible_directions=('x','y','z')
        self.coords=dict()
        for varname in self.variables:
            size = self.variables[varname]['size']
            if size == 1:
                # Look for position values
                direction = varname.lower()
                key = self.coordinates
                if key in self.possible_coords and direction in self.possible_directions:
                    if key not in self.coords:
                        self.coords[key] = dict(coord=key)
                    self.coords[key]['size'] = size
                    self.coords[key][direction] = varname

        # Change 'fill' values in data to NaNs
        self.fill2nan()
        
    def parse_file(self):
        import re
        vars=[]
        units=[]
        times=[]
        arrays = []
        print("===> Printing File Header ...")
        for line in self.file.splitlines(False):
            A = re.match('^# ', line)
            B = re.match('# Run', line)
            C = re.match('# Coordinate', line)
            D = re.match('# Satellite', line)
            E = re.match('# Year', line)
            F = re.match('# \[year\]', line)
            if A or B or C or D or E or F:
                if A:
                    print("-> ",line)
                if B:
                    # Find runname and fill value
                    parts=re.sub(' +', ' ', line).split(' ')
                    self.runname = parts[3]
                    self.fillvalue = parts[6]
                if C:
                    # Check that coordinate system matches
                    parts=re.sub(' +', ' ', line).split(' ')
                    if self.coordinates != parts[3]:
                        print("ERROR: Coordinate system does not match.",self.coordinates,parts[3])
                if D:
                    # Check that satellite name matches
                    parts=re.sub(' +', ' ', line).split(' ')
                    if self.satellite != parts[3]:
                        print("ERROR: Satellite does not match.",self.satellite,parts[3])
                if E:
                    # Variable names, remove . and change N and B_1
                    parts=re.sub(' +', ' ', line).strip().split(' ')
                    for p in parts[7:]:
                        p=re.sub("\.","",p)
                        p=re.sub("B_1","B1",p)
                        p=re.sub("^N$","rho",p)
                        vars.append(p)
                if F:
                    # Variable units, remove [] and fix exponents
                    parts=re.sub(' +', ' ', line).strip().split(' ')
                    for p in parts[7:]:
                        p=re.sub("cm\^-3","1/cm^3",p)
                        p=re.sub("m2","m^2",p)
                        p=re.sub("m3","m^3",p)
                        p=re.sub("\[","",p)
                        p=re.sub("\]","",p)
                        units.append(p)
            else:
                parts=re.sub(' +', ' ', line).strip().split(' ')
                year=parts[0]
                month=parts[1]
                day=parts[2]
                hour=parts[3]
                minute=parts[4]
                second=parts[5]
                ms=0
                if '.' in second:
                    (second,ms)=second.split('.')
                dd=datetime.datetime(int(year),int(month),int(day),
                                     hour=int(hour),minute=int(minute),second=int(second),
                                     microsecond=int(ms)*1000,tzinfo=datetime.timezone.utc)
                times.append(dd)
                for s in parts[6:]:
                    arrays.append(float(s))
        # self.dtarray=np.array([dd for dd in times])
        self.dtarray = pd.to_datetime(times)
        self.tsarray = np.array([d.timestamp() for d in self.dtarray])
        
        nvar=len(vars)
        nval=len(arrays)
        npos=int(nval/nvar)           
        arrays=np.array(arrays)
        arrays=arrays.reshape((npos,nvar))

        i=0
        for var in vars:
            self.variables[var] = dict(units=units[i],
                                       data=arrays[:,i],
                                       size=1,
                                       fill=self.fillvalue)
            i+=1
        
        return
                
    def register_variable(self, varname, units):
        """register variables into Kamodo for this service, CCMC ROR satellite extractions"""

        data =  self.variables[varname]['data']
        times = self.dtarray

        ser = pd.Series(data, index=pd.DatetimeIndex(times))


        @kamodofy(units = units, 
                  citation = "De Zeeuw 2020",
                  data = None)
        def interpolate(t=times):
            ts = t
            isiterable = isinstance(t, Iterable)

            if isinstance(ts, DatetimeIndex):
                pass
            elif isinstance(ts, float):
                ts = pd.to_datetime([ts], utc=True, unit='s')
            elif isiterable:
                if isinstance(ts[0], float):
                    ts = pd.to_datetime(ts, utc=True, unit='s')
                ts = DatetimeIndex(ts)
            else:
                raise NotImplementedError(ts)
            
            ser_ = ser.reindex(ser.index.union(ts))
            ser_interpolated = ser_.interpolate(method='time')
            result = ser_interpolated.reindex(ts)
            if isiterable:
                return result.values
            else:
                return result.values[0]

        # store the interpolator
        self.variables[varname]['interpolator'] = interpolate

        # update docstring for this variable
        interpolate.__doc__ = "A function that returns {} in [{}].".format(varname,units)

        var_reg = '{}__{}'.format(varname, self.coordinates)
        self[var_reg] = interpolate

    def fill2nan(self):
        '''
        Replaces fill value in data with NaN.
        Not Yet called by default. Call as needed.
        '''
        for varname in self.variables:
            data = self.variables[varname]['data']
            fill = self.variables[varname]['fill']
            if fill is not None:
                mask = data==float(fill)
                nbad = np.count_nonzero(mask)
                if nbad > 0:
                    print("Found",nbad,"fill values, replacing with NaN for variable",varname,
                          "of size",data.size)
                data[mask]=np.nan
                self.variables[varname]['data'] = data
        
    def get_plot(self, type="1Dpos", scale="R_E", var=""):
        '''
        Return a plotly figure object.
        type = 1Dvar => 1D plot of variable value vs Time
               1Dpos (default) => 1D location x,y,z vs Time
               3Dpos => 3D location colored by altitude
        scale = km, R_E (default)
        var = variable name for variable value plots
        '''

        coord=self.coordinates
        
        # Set plot title for plots
        txttop=self.satellite + " position extracted from run " + self.runname + "<br>"\
            + self.start + " to " + self.stop + "<br>" + coord
        
        if type == '1Dvar':
            if var == "":
                print("No plot variable passed in.")
                return
            fig=go.Figure()
            if self.variables[var]['size'] == 1:
                x=self.variables[var]['data']
                fig.add_trace(go.Scatter(x=self.dtarray, y=x, mode='lines+markers', name=var))
            elif self.variables[var]['size'] == 3:
                x=self.variables[var]['data'][:,0]
                y=self.variables[var]['data'][:,1]
                z=self.variables[var]['data'][:,2]
                fig.add_trace(go.Scatter(x=self.dtarray, y=x, mode='lines+markers', name=var))
                fig.add_trace(go.Scatter(x=self.dtarray, y=y, mode='lines+markers', name=var))
                fig.add_trace(go.Scatter(x=self.dtarray, y=z, mode='lines+markers', name=var))
            ytitle=var+" ["+self.variables[var]['units']+"]"
            fig.update_xaxes(title_text="Time")
            fig.update_yaxes(title_text=ytitle)
            fig.update_layout(hovermode="x")
            fig.update_layout(title_text=txttop)
 
            return fig
        
        if type == '1Dpos':
            fig=go.Figure()
            xvarname = self.coords[coord]['x']
            if self.coords[coord]['size'] == 1:
                x=self.variables[self.coords[coord]['x']]['data']
                y=self.variables[self.coords[coord]['y']]['data']
                z=self.variables[self.coords[coord]['z']]['data']
            elif self.coords[coord]['size'] == 3:
                x=self.variables[self.coords[coord]['x']]['data'][:,0]
                y=self.variables[self.coords[coord]['y']]['data'][:,1]
                z=self.variables[self.coords[coord]['z']]['data'][:,2]
            if scale == "km":
                if self.variables[xvarname]['units'] == "R_E":
                    x=x*self.RE
                    y=y*self.RE
                    z=z*self.RE
                ytitle="Position [km]"
            else:
                if self.variables[xvarname]['units'] == "km":
                    x=x/self.RE
                    y=y/self.RE
                    z=z/self.RE
                ytitle="Position [R_E]"
            fig.add_trace(go.Scatter(x=self.dtarray, y=x,
                                     mode='lines+markers', name=self.coords[coord]['x']))
            fig.add_trace(go.Scatter(x=self.dtarray, y=y,
                                     mode='lines+markers', name=self.coords[coord]['y']))
            fig.add_trace(go.Scatter(x=self.dtarray, y=z,
                                     mode='lines+markers', name=self.coords[coord]['z']))
            fig.update_xaxes(title_text="Time")
            fig.update_yaxes(title_text=ytitle)
            fig.update_layout(hovermode="x")
            fig.update_layout(title_text=txttop)
 
            return fig
        
        if type == "3Dpos":
            xvarname = self.coords[coord]['x']
            if self.coords[coord]['size'] == 1:
                x=self.variables[self.coords[coord]['x']]['data']
                y=self.variables[self.coords[coord]['y']]['data']
                z=self.variables[self.coords[coord]['z']]['data']
            elif self.coords[coord]['size'] == 3:
                x=self.variables[self.coords[coord]['x']]['data'][:,0]
                y=self.variables[self.coords[coord]['y']]['data'][:,1]
                z=self.variables[self.coords[coord]['z']]['data'][:,2]
            if scale == "km":
                if self.variables[xvarname]['units'] == "R_E":
                    x=x*self.RE
                    y=y*self.RE
                    z=z*self.RE
                r=(np.sqrt(x**2 + y**2 + z**2))-self.RE
                ytitle="Position [km]"
            else:
                if self.variables[xvarname]['units'] == "km":
                    x=x/self.RE
                    y=y/self.RE
                    z=z/self.RE
                r=(np.sqrt(x**2 + y**2 + z**2))-1.
                ytitle="Position [R_E]"
            fig=px.scatter_3d(
                x=x,
                y=y,
                z=z,
                color=r)
            bartitle = "Altitude [" + scale + "]"
            fig.update_layout(coloraxis=dict(colorbar=dict(title=bartitle)))
            fig.update_layout(scene=dict(xaxis=dict(title=dict(text="X ["+scale+"]")),
                                         yaxis=dict(title=dict(text="Y ["+scale+"]")),
                                         zaxis=dict(title=dict(text="Z ["+scale+"]"))))
            fig.update_layout(title_text=txttop)
            return fig
        

        print('ERROR, reached end of get_plot without any action taken.')
        return

