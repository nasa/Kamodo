from kamodo import Kamodo, kamodofy
import pandas as pd
import numpy as np
import scipy
import time
import datetime
from datetime import timezone
import urllib, json
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
from pandas import DatetimeIndex
from collections.abc import Iterable


def ror_get_info(runID):
    '''Query run information for given runID'''
    server = "https://ccmc.gsfc.nasa.gov/RoR_WWW/VMR"
    query = '{}/files.php?id={}'.format(server, runID)
    response = urllib.request.urlopen(query)
    data = json.loads(response.read())
    return data

def ror_show_info(runID):
    '''Display run information for a given runID.'''
    result=ror_get_info(runID)
    print("Run information for runID =",runID,"from the CCMC RoR system.")
    for item in result['info']:
        for k in item:
            print(k.rjust(25),':',item[k])
    sats = []
    for sat in result['satellites']:
        satname = sat['name']
        sats.append(satname)
    k='Satellite extractions'
    print(k.rjust(25),':',sats)

def ror_return_satellites(runID):
    '''Display list of satellites as python array for given runID.'''
    result=ror_get_info(runID)
    sats = []
    for sat in result['satellites']:
        satname = sat['name']
        sats.append(satname)
    return sats

def ror_get_extraction(runID, coord, satellite):
    '''Query for file contents from server'''
    server = "https://ccmc.gsfc.nasa.gov/RoR_WWW/VMR"
    query = '{}/{}/{}/{}_{}.txt'.format(server, runID, satellite, coord, satellite)
    # print(query)
    response = urllib.request.urlopen(query)
    file = response.read()
    return file

class SATEXTRACT(Kamodo):
    def __init__(self,
            runID, # ccmc runs-on-request run id
            coord, # coordinate system
            satellite, # satellite
            debug=0,
            server="https://ccmc.gsfc.nasa.gov/RoR_WWW/VMR",
            **kwargs):
        super(SATEXTRACT, self).__init__(**kwargs)
        self.verbose=False # overrides kwarg
        # self.symbol_registry=dict() # wipes any user-defined symbols
        # self.signatures=dict() # wipes any user-defined signatures
        self.RE=6.3781E3
        self.server = server # should be set by keyword
        self.runID = runID
        self.coordinates = coord
        self.satellite = satellite
        self.debug = debug
        if self.debug > 0:
            print(' -server: CCMC RoR')
            print(' -runID: ',runID)
            print(' -coordinate system: ',coord)
            print(' -satellite: ',satellite)
        self.variables=dict()
        self.file = ror_get_extraction(runID, coord, satellite).decode('ascii')
        self.parse_file()
        ts=self.tsarray[0]
        self.start = datetime.datetime.fromtimestamp(ts,tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        ts=self.tsarray[-1]
        self.stop  = datetime.datetime.fromtimestamp(ts,tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        if self.debug > 0:
            print(" ")
            print(" -date start: ",self.start)
            print("         end: ",self.stop)

        for varname in self.variables:
            if varname == "N":
                continue
            units = self.variables[varname]['units']
            if self.debug > 0:
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
        if self.debug > 0:
            print("===> Printing File Header ...")
        for line in self.file.splitlines(False):
            A = re.match('^# ', line)
            B = re.match('# Run', line)
            C = re.match('# Coordinate', line)
            D = re.match('# Satellite', line)
            E = re.match('# Year', line)
            F = re.match('# \[year\]', line)
            G = re.match('# Data type', line)
            if A or B or C or D or E or F or G:
                if A:
                    if self.debug > 0:
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
                    if self.debug > 1:
                        print(len(vars), vars)
                if F:
                    # Variable units, remove [] and fix exponents
                    parts=re.sub(' +', ' ', line).strip().split(' ')
                    if self.modelname == "GUMICS":
                        # missing x,y,z  --put before other units
                        units.append('R_E')
                        units.append('R_E')
                        units.append('R_E')
                    for p in parts[7:]:
                        if self.modelname == "BATSRUS":
                            # need a unit applied for status variable, currently missing in some
                            if len(units) != len(vars):
                                if vars[len(units)] == "status":
                                    units.append('')
                        p=re.sub("cm\^-3","1/cm^3",p)
                        p=re.sub("m2","m^2",p)
                        p=re.sub("m3","m^3",p)
                        p=re.sub("\[","",p)
                        p=re.sub("\]","",p)
                        p=re.sub("Vm/A","V/A",p) # This is wrong but use it for now
                        p=re.sub("nJ/m","J/m",p) # This is wrong but use it for now
                        units.append(p)
                        if self.modelname == "LFM" and p == "nPa":
                            # missing X,Y,Z,Vol  --unknown Vol units
                            units.append('R_E')
                            units.append('R_E')
                            units.append('R_E')
                            units.append('')
                        if self.modelname == "OpenGGCM" and p == "V/A":
                            # missing eflx,efly,eflz  --units are unknown
                            units.append('')
                            units.append('')
                            units.append('')
                    if self.debug > 1:
                        print(len(units), units)
                if G:
                    # Pull out the model name
                    parts=re.sub(' +', ' ', line).split(' ')
                    self.modelname = parts[3]
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
        self.dtarrayclean = np.array([datetime.datetime.fromtimestamp(d,tz=timezone.utc).strftime("%Y-%m-%d %H:%M:%S") for d in self.tsarray])
        
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
                    if self.debug > 0:
                        print("Found",nbad,"fill values, replacing with NaN for variable",
                              varname,"of size",data.size)
                data[mask]=np.nan
                self.variables[varname]['data'] = data

    def info(self):
        '''
        Show information stored in this Kamodo object.
        '''
        print('Kamodo object information:')
        print('  server         CCMC RoR')
        print('  runID         ',self.runID)
        print('  coordinates   ',self.coordinates)
        print('  satellite     ',self.satellite)
        print('  run start     ',self.start)
        print('  run end       ',self.stop)
        print('  variables')
        for varname in self.variables:
            units = self.variables[varname]['units']
            print('     ',varname,' [',units,']')

    def get_plot(self, type="1Dpos", scale="R_E", var="", groupby="all",
                 quiver=False, quiverscale="5.", quiverskip="0"):
        '''
        Return a plotly figure object.
        type = 1Dvar => 1D plot of variable value vs Time (also all variables option)
               1Dpos (default) => 1D location x,y,z vs Time
               3Dpos => 3D location colored by altitude
               3Dvar => View variable on 3D position (also quiver and groupby options)
        scale = km, R_E (default)
        var = variable name for variable value plots
        groupby = day, hour, all (default) => groupings for 3Dvar plots
        quiver = True, False (default) => if var is a vector value and 3Dvar plot, then
               turn on quivers and color by vector magnitude
        quiverscale = 5. (default) => length of quivers in units of RE
        quiverskip =  (default) => now many quivers to skip displaying
        '''
        quiverscale=float(quiverscale)
        if scale == "km":
            quiverscale=quiverscale*self.RE
        quiverskip=int(quiverskip)
        coord=self.coordinates
        
        # Set plot title for plots
        txttop=self.satellite + " position extracted from run " + self.runname + "<br>"\
            + self.start + " to " + self.stop + "<br>" + coord
        
        if type == '1Dvar':
            if var == "":
                print("No plot variable passed in.")
                return
            fig=go.Figure()
            if var == "all":
                # Create menu pulldown for each variable
                steps = []
                i=0
                for varname in self.variables:
                    ytitle=varname+" ["+self.variables[varname]['units']+"]"
                    step = dict(
                        label=ytitle,
                        method="update",
                        args=[{"visible": [False] * len(self.variables)}],
                    )
                    step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
                    steps.append(step)
                    if self.variables[varname]['size'] == 1:
                        x=self.variables[varname]['data']
                        fig.add_trace(go.Scatter(x=self.dtarray, y=x, mode='lines+markers',
                                                 name=varname, visible=False))
                    elif self.variables[varname]['size'] == 3:
                        x=self.variables[varname]['data'][:,0]
                        y=self.variables[varname]['data'][:,1]
                        z=self.variables[varname]['data'][:,2]
                        fig.add_trace(go.Scatter(x=self.dtarray, y=x, mode='lines+markers',
                                                 name=varname, visible=False))
                        fig.add_trace(go.Scatter(x=self.dtarray, y=y, mode='lines+markers',
                                                 name=varname, visible=False))
                        fig.add_trace(go.Scatter(x=self.dtarray, y=z, mode='lines+markers',
                                                 name=varname, visible=False))
                    i+=1
                fig.data[0].visible=True
                fig.update_layout(updatemenus=list([dict(buttons=steps)]))
                fig.update_xaxes(title_text="Time")
                fig.update_layout(hovermode="x")
                fig.update_layout(title_text=txttop)
                
            else:
                # Standard single/vector variable plot
                varn=var
                if var in self.variables:
                    size=self.variables[var]['size']
                else:
                    if var+'_x' in self.variables:
                        size=0
                    else:
                        print("Invalid plot variable passed in.")
                        return
                if size == 0:
                    var=varn+'_x'
                    x=self.variables[var]['data']
                    fig.add_trace(go.Scatter(x=self.dtarray, y=x, mode='lines+markers', name=var))
                    var=varn+'_y'
                    x=self.variables[var]['data']
                    fig.add_trace(go.Scatter(x=self.dtarray, y=x, mode='lines+markers', name=var))
                    var=varn+'_z'
                    x=self.variables[var]['data']
                    fig.add_trace(go.Scatter(x=self.dtarray, y=x, mode='lines+markers', name=var))
                elif size == 1:
                    x=self.variables[var]['data']
                    fig.add_trace(go.Scatter(x=self.dtarray, y=x, mode='lines+markers', name=var))
                elif size == 3:
                    x=self.variables[var]['data'][:,0]
                    y=self.variables[var]['data'][:,1]
                    z=self.variables[var]['data'][:,2]
                    fig.add_trace(go.Scatter(x=self.dtarray, y=x, mode='lines+markers', name=var))
                    fig.add_trace(go.Scatter(x=self.dtarray, y=y, mode='lines+markers', name=var))
                    fig.add_trace(go.Scatter(x=self.dtarray, y=z, mode='lines+markers', name=var))
                ytitle=varn+" ["+self.variables[var]['units']+"]"
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
            else:
                if self.variables[xvarname]['units'] == "km":
                    x=x/self.RE
                    y=y/self.RE
                    z=z/self.RE
                r=(np.sqrt(x**2 + y**2 + z**2))-1.
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
        
        if type == "3Dvar":
            if var == "":
                print("No plot variable passed in.")
                return
            xvarname = self.coords[coord]['x']
            vard=self.variables[var]['data']
            varu=self.variables[var]['units']

            if quiver:
                if "_x" in var or "_y" in var or "_z" in var:
                    var=var.split('_')[0]
                    qxvar=var+"_x"
                    qyvar=var+"_y"
                    qzvar=var+"_z"
                    qxvard=self.variables[qxvar]['data']
                    qyvard=self.variables[qyvar]['data']
                    qzvard=self.variables[qzvar]['data']
                    vard = np.sqrt(np.square(qxvard) +\
                                   np.square(qyvard) +\
                                   np.square(qzvard))
                else:
                    print("A vector variable was not passed, turning quiver off.")
                    quiver=False

            cmin=np.amin(vard)
            cmax=np.amax(vard)
            if self.coords[coord]['size'] == 1:
                x=self.variables[self.coords[coord]['x']]['data']
                y=self.variables[self.coords[coord]['y']]['data']
                z=self.variables[self.coords[coord]['z']]['data']
            elif self.coords[coord]['size'] == 3:
                x=self.variables[self.coords[coord]['x']]['data'][:,0]
                y=self.variables[self.coords[coord]['y']]['data'][:,1]
                z=self.variables[self.coords[coord]['z']]['data'][:,2]
            bodyscale=1.    
            if scale == "km":
                if self.variables[xvarname]['units'] == "R_E":
                    x=x*self.RE
                    y=y*self.RE
                    z=z*self.RE
                    bodyscale=self.RE    
                r=(np.sqrt(x**2 + y**2 + z**2))-self.RE
            else:
                if self.variables[xvarname]['units'] == "km":
                    x=x/self.RE
                    y=y/self.RE
                    z=z/self.RE
                r=(np.sqrt(x**2 + y**2 + z**2))-1.
            xmin=np.amin(x)
            xmax=np.amax(x)
            ymin=np.amin(y)
            ymax=np.amax(y)
            zmin=np.amin(z)
            zmax=np.amax(z)
            if quiver:
                tmpx=x+(qxvard*quiverscale/cmax)
                tmpy=y+(qyvard*quiverscale/cmax)
                tmpz=z+(qzvard*quiverscale/cmax)
                xmin=min(xmin,np.amin(tmpx))
                xmax=max(xmax,np.amax(tmpx))
                ymin=min(ymin,np.amin(tmpy))
                ymax=max(ymax,np.amax(tmpy))
                zmin=min(zmin,np.amin(tmpz))
                zmax=max(zmax,np.amax(tmpz))
            
            # Create empty figure to build pieces into
            fig=go.Figure()
            Nplot = 0

            # Start 1 RE sphere, padded to cover all data positions
            dataXYZ = pd.read_csv('https://ccmc.gsfc.nasa.gov/Kamodo/demo/sphereXYZ.csv')
            dataIJK = pd.read_csv('https://ccmc.gsfc.nasa.gov/Kamodo/demo/sphereIJK.csv')
            if scale == "km":
                dataXYZ=dataXYZ*self.RE
            fig.add_mesh3d()
            fig.data[Nplot].x = np.append(dataXYZ['x'],(xmin,xmax))
            fig.data[Nplot].y = np.append(dataXYZ['y'],(ymin,ymax))
            fig.data[Nplot].z = np.append(dataXYZ['z'],(zmin,zmax))
            fig.data[Nplot].i = dataIJK['i']
            fig.data[Nplot].j = dataIJK['j']
            fig.data[Nplot].k = dataIJK['k']
            fig.data[Nplot].facecolor = dataIJK['c']
            fig.data[Nplot].flatshading = True
            fig.data[Nplot].name = '1 R_E sphere'
            fig.data[Nplot].hovertemplate="Earth<extra></extra>"
            fig.data[Nplot].visible=True
            Nplot += 1

            # Array of 'groupby' values and unique values
            if groupby == "day":
                datearray = np.array([datetime.datetime.fromtimestamp(d,tz=timezone.utc).strftime\
                                      ("%Y-%m-%d") for d in self.tsarray])
            elif groupby == "hour":
                datearray = np.array([datetime.datetime.fromtimestamp(d,tz=timezone.utc).strftime\
                                      ("%Y-%m-%d %H") for d in self.tsarray])
            else:
                datearray = np.array([datetime.datetime.fromtimestamp(d,tz=timezone.utc).strftime\
                                      ("all") for d in self.tsarray])
            udates = np.unique(datearray)
            nsteps = len(udates)

            for date in udates:
                # Compute mast to restrict all data in trace
                mask = date == datearray
                # Create figure with data
                if quiver:
                    fig.add_trace(go.Scatter3d(
                        name=date,  x=x[mask], y=y[mask], z=z[mask],
                        mode='markers',
                        marker=dict(
                            size=4,
                            cmin=cmin, cmax=cmax,
                            color=vard[mask],
                            showscale=True,
                            colorscale='Viridis',
                            colorbar=dict(title=var+" ["+varu+"]"),
                        ),
                        customdata=np.vstack((self.dtarrayclean[mask],\
                                              qxvard[mask],\
                                              qyvard[mask],\
                                              qzvard[mask])).T,
                        hovertemplate="<b>"+self.satellite+" Position</b>"+
                        "<br>X: %{x:.4f}<br>Y: %{y:.4f}<br>Z: %{z:.4f}<br>"+
                        qxvar+": %{customdata[1]:.2f}<br>"+
                        qyvar+": %{customdata[2]:.2f}<br>"+
                        qzvar+": %{customdata[3]:.2f}<br>"+
                        " "+var+": %{marker.color:.2f}"+varu+"<br>"+
                        "%{customdata[0]}<br>"+"<extra></extra>",
                        visible=False,
                    ))
                    Nplot += 1
                else:
                    fig.add_trace(go.Scatter3d(
                        name=date,  x=x[mask], y=y[mask], z=z[mask],
                        mode='markers',
                        marker=dict(
                            size=4,
                            cmin=cmin, cmax=cmax,
                            color=vard[mask],
                            showscale=True,
                            colorscale='Viridis',
                            colorbar=dict(title=var+" ["+varu+"]"),
                        ),
                        customdata=self.dtarrayclean[mask],
                        hovertemplate="<b>"+self.satellite+" Position</b>"+
                        "<br>X: %{x:.4f}<br>Y: %{y:.4f}<br>Z: %{z:.4f}<br>"+
                        " "+var+": %{marker.color:.2f}"+varu+"<br>"+
                        "%{customdata}<br>"+"<extra></extra>",
                        visible=False,
                    ))
                    Nplot += 1
            fig.data[1].visible=True

            if quiver:
                for date in udates:
                    # Compute mask to restrict all data in trace
                    mask = date == datearray
                    # Precompute needed values
                    xm=x[mask]
                    ym=y[mask]
                    zm=z[mask]
                    vm=vard[mask]
                    qxm=qxvard[mask]
                    qym=qyvard[mask]
                    qzm=qzvard[mask]
                    xm=xm.reshape(len(xm),1)
                    ym=ym.reshape(len(ym),1)
                    zm=zm.reshape(len(zm),1)
                    vm=vm.reshape(len(vm),1)
                    qxm=qxm.reshape(len(qxm),1)
                    qym=qym.reshape(len(qym),1)
                    qzm=qzm.reshape(len(qzm),1)
                    if quiverskip > 0:
                        for i in range(len(qxm)):
                            if i%(quiverskip+1) > 0:
                                qxm[i]=0.
                                qym[i]=0.
                                qzm[i]=0.
                    xx=np.concatenate((xm,xm+qxm*quiverscale/cmax,xm),axis=1).reshape(3*len(xm))
                    yy=np.concatenate((ym,ym+qym*quiverscale/cmax,ym),axis=1).reshape(3*len(ym))
                    zz=np.concatenate((zm,zm+qzm*quiverscale/cmax,zm),axis=1).reshape(3*len(zm))
                    # Create figure with data
                    fig.add_trace(go.Scatter3d(
                        name=date,
                        x = xx,
                        y = yy,
                        z = zz,
                        mode='lines',
                        line=dict(
                            width=2,
                            color='rgba(22,22,22,0.2)',
                        ),
                        hoverinfo='skip',
                        visible=False,
                    ))
                    Nplot += 1
                fig.data[1+nsteps].visible=True

            # Start selection slider build
            steps = []
            for i in range(nsteps):
                step = dict(
                    method="update",
                    args=[{"visible": [False] * len(fig.data)}],
                    label=udates[i]
                )
                step["args"][0]["visible"][0] = True  # Set first trace to "visible"
                step["args"][0]["visible"][i+1] = True  # Toggle i'th trace to "visible"
                if quiver:
                    step["args"][0]["visible"][i+1+nsteps] = True  # Toggle i'th trace to "visible"
                steps.append(step)

            sliders = [dict(
                active=0,
                currentvalue={"prefix": "Currently showing: "},
                pad={"t": 50},
                steps=steps
            )]
            fig.update_layout(sliders=sliders)

            # Update axis labels and plot title
            fig.update_layout(scene=dict(xaxis=dict(title=dict(text="X ["+scale+"]")),
                                         yaxis=dict(title=dict(text="Y ["+scale+"]")),
                                         zaxis=dict(title=dict(text="Z ["+scale+"]"))))
            fig.update_layout(title_text=txttop)
            fig.update_layout(showlegend=False)
            
            return fig
        

        print('ERROR, reached end of get_plot without any action taken.')
        return

#======
# New class to collect all satellites
#
class SATEXTRACTALL(Kamodo):
    def __init__(self, runID, coord, debug=0, sats=[], **kwargs):
        super(SATEXTRACTALL, self).__init__(**kwargs)
        # Start timer
        tic = time.perf_counter()
        self.verbose=False
        self.symbol_registry=dict()
        self.signatures=dict()
        self.RE=6.3781E3
        self.runID = runID
        self.coordinates = coord
        self.debug = debug
        if self.debug > 0:
            print(' -server: CCMC RoR')
            print(' -runID: ',runID)
            print(' -coordinate system: ',coord)
        self.satellites = dict()
        result=ror_get_info(runID)
        for sat in result['satellites']:
            satname = sat['name']
            if len(sats) == 0 or satname in sats:
                ror = SATEXTRACT(runID, coord, satname, debug=debug)
                self.satellites[satname] = ror
        self.start=ror.start
        self.stop=ror.stop
        self.runname=ror.runname
        self.modelname=ror.modelname
        if self.debug > 0:
            print(' -data extracted from model: ',self.modelname)
        # end timer
        toc = time.perf_counter()
        print(f"Time loading files and registering satellites: {toc - tic:0.4f} seconds")

    def info(self):
        '''
        Show information stored in this Kamodo object.
        '''
        print('Kamodo object information:')
        print('  server         CCMC RoR')
        print('  runID         ',self.runID)
        print('  coordinates   ',self.coordinates)
        print('  satellites')
        for sat in self.satellites:
            print('     ',sat)
        print('  run start     ',self.start)
        print('  run end       ',self.stop)
        print('  variables')
        for sat in self.satellites:
            ror = self.satellites[sat]
            for varname in ror.variables:
                units = ror.variables[varname]['units']
                print('     ',varname,' [',units,']')
            break

    def get_plot(self, type="3Dvar", scale="R_E", var="", groupby="all",
                 quiver=False, quiverscale="5.", quiverskip="0"):
        '''
        Return a plotly figure object.
        type = 3Dvar => View variable on 3D position (also quiver and groupby options)
        scale = km, R_E (default)
        var = variable name for variable value plots
        groupby = day, hour, all (default) => groupings for 3Dvar plots
        quiver = True, False (default) => if var is a vector value and 3Dvar plot, then
               turn on quivers and color by vector magnitude
        quiverscale = 5. (default) => length of quivers in units of RE
        quiverskip =  (default) => now many quivers to skip displaying
        '''
        quiverscale=float(quiverscale)
        if scale == "km":
            quiverscale=quiverscale*self.RE
        quiverskip=int(quiverskip)
        coord=self.coordinates

        # Set plot title for plots
        txttop="Satellite positions extracted from run " + self.runname + "<br>"\
            + self.start + " to " + self.stop + "<br>" + coord

        # set initial values used later
        xmin=0.
        xmax=0.
        ymin=0.
        ymax=0.
        zmin=0.
        zmax=0.
        cmin= 1.e99
        cmax=-1.e99

        if type == "3Dvar":
            if var == "":
                print("No plot variable passed in.")
                return

            # Create empty figure to build pieces into
            fig=go.Figure()
            Nplot = 0

            # Pre-loop to find cmin,cmax to use in all plot pieces
            for sat in self.satellites:
                ror = self.satellites[sat]
                vard=ror.variables[var]['data']

                if quiver:
                    if "_x" in var or "_y" in var or"_z" in var:
                        var2=var.split('_')[0]
                        qxvar=var2+"_x"
                        qyvar=var2+"_y"
                        qzvar=var2+"_z"
                        qxvard=ror.variables[qxvar]['data']
                        qyvard=ror.variables[qyvar]['data']
                        qzvard=ror.variables[qzvar]['data']
                        vard = np.sqrt(np.square(qxvard) +\
                                       np.square(qyvard) +\
                                       np.square(qzvard))
                    else:
                        print("A vector variable was not passed, turning quiver off.")
                        quiver=False

                cmin=min(cmin,np.amin(vard))
                cmax=max(cmax,np.amax(vard))
            
            # Actual plot creation loop
            for sat in self.satellites:
                ror = self.satellites[sat]
                vard=ror.variables[var]['data']
                varu=ror.variables[var]['units']

                if quiver:
                    if "_x" in var or "_y" in var or"_z" in var:
                        var2=var.split('_')[0]
                        qxvar=var2+"_x"
                        qyvar=var2+"_y"
                        qzvar=var2+"_z"
                        qxvard=ror.variables[qxvar]['data']
                        qyvard=ror.variables[qyvar]['data']
                        qzvard=ror.variables[qzvar]['data']
                        vard = np.sqrt(np.square(qxvard) +\
                                       np.square(qyvard) +\
                                       np.square(qzvard))

                if ror.coords[coord]['size'] == 1:
                    x=ror.variables[ror.coords[coord]['x']]['data']
                    y=ror.variables[ror.coords[coord]['y']]['data']
                    z=ror.variables[ror.coords[coord]['z']]['data']
                elif ror.coords[coord]['size'] == 3:
                    x=ror.variables[ror.coords[coord]['x']]['data'][:,0]
                    y=ror.variables[ror.coords[coord]['y']]['data'][:,1]
                    z=ror.variables[ror.coords[coord]['z']]['data'][:,2]
                bodyscale=1.
                if scale == "km":
                    if ror.variables[ror.coords[coord]['x']]['units'] == "R_E":
                        x=x*ror.RE
                        y=y*ror.RE
                        z=z*ror.RE
                        bodyscale=ror.RE    
                    r=(np.sqrt(x**2 + y**2 + z**2))-ror.RE
                else:
                    if ror.variables[ror.coords[coord]['x']]['units'] == "km":
                        x=x/ror.RE
                        y=y/ror.RE
                        z=z/ror.RE
                    r=(np.sqrt(x**2 + y**2 + z**2))-1.
            
                # Array of 'groupby' values and unique values
                if groupby == "day":
                    datearray = np.array([datetime.datetime.fromtimestamp(d,tz=timezone.utc).strftime\
                                          ("%Y-%m-%d") for d in ror.tsarray])
                elif groupby == "hour":
                    datearray = np.array([datetime.datetime.fromtimestamp(d,tz=timezone.utc).strftime\
                                          ("%Y-%m-%d %H") for d in ror.tsarray])
                else:
                    datearray = np.array([datetime.datetime.fromtimestamp(d,tz=timezone.utc).strftime\
                                          ("all") for d in ror.tsarray])
                udates = np.unique(datearray)
                nsteps = len(udates)

                for date in udates:
                    # Compute mast to restrict all data in trace
                    mask = date == datearray
                    # Create figure with data
                    if quiver:
                        fig.add_trace(go.Scatter3d(
                            name=date,  x=x[mask], y=y[mask], z=z[mask],
                            mode='markers',
                            marker=dict(
                                size=4,
                                cmin=cmin, cmax=cmax,
                                color=vard[mask],
                                showscale=True,
                                colorscale='Viridis',
                                colorbar=dict(title=var2+" ["+varu+"]"),
                            ),
                            customdata=np.vstack((ror.dtarrayclean[mask],\
                                                  qxvard[mask],\
                                                  qyvard[mask],\
                                                  qzvard[mask])).T,
                            hovertemplate="<b>"+ror.satellite+" Position</b>"+
                            "<br>X: %{x:.4f}<br>Y: %{y:.4f}<br>Z: %{z:.4f}<br>"+
                            qxvar+": %{customdata[1]:.2f}<br>"+
                            qyvar+": %{customdata[2]:.2f}<br>"+
                            qzvar+": %{customdata[3]:.2f}<br>"+
                            " "+var2+": %{marker.color:.2f}"+varu+"<br>"+
                            "%{customdata[0]}<br>"+"<extra></extra>",
                            visible=False,
                        ))
                        Nplot += 1
                    else:
                        fig.add_trace(go.Scatter3d(
                            name=date,  x=x[mask], y=y[mask], z=z[mask],
                            mode='markers',
                            marker=dict(
                                size=4,
                                cmin=cmin, cmax=cmax,
                                color=vard[mask],
                                showscale=True,
                                colorscale='Viridis',
                                colorbar=dict(title=var+" ["+varu+"]"),
                            ),
                            customdata=ror.dtarrayclean[mask],
                            hovertemplate="<b>"+ror.satellite+" Position</b>"+
                            "<br>X: %{x:.4f}<br>Y: %{y:.4f}<br>Z: %{z:.4f}<br>"+
                            " "+var+": %{marker.color:.2f}"+varu+"<br>"+
                            "%{customdata}<br>"+"<extra></extra>",
                            visible=False,
                        ))
                        Nplot += 1

                if quiver:
                    for date in udates:
                        # Compute mask to restrict all data in trace
                        mask = date == datearray
                        # Precompute needed values
                        xm=x[mask]
                        ym=y[mask]
                        zm=z[mask]
                        vm=vard[mask]
                        qxm=qxvard[mask]
                        qym=qyvard[mask]
                        qzm=qzvard[mask]
                        xm=xm.reshape(len(xm),1)
                        ym=ym.reshape(len(ym),1)
                        zm=zm.reshape(len(zm),1)
                        vm=vm.reshape(len(vm),1)
                        qxm=qxm.reshape(len(qxm),1)
                        qym=qym.reshape(len(qym),1)
                        qzm=qzm.reshape(len(qzm),1)
                        if quiverskip > 0:
                            for i in range(len(qxm)):
                                if i%(quiverskip+1) > 0:
                                    qxm[i]=0.
                                    qym[i]=0.
                                    qzm[i]=0.
                        xx=np.concatenate((xm,xm+qxm*quiverscale/cmax,xm),axis=1).reshape(3*len(xm))
                        yy=np.concatenate((ym,ym+qym*quiverscale/cmax,ym),axis=1).reshape(3*len(ym))
                        zz=np.concatenate((zm,zm+qzm*quiverscale/cmax,zm),axis=1).reshape(3*len(zm))
                        # Create figure with data
                        fig.add_trace(go.Scatter3d(
                            name=date,
                            x = xx,
                            y = yy,
                            z = zz,
                            mode='lines',
                            line=dict(
                                width=2,
                                color='rgba(22,22,22,0.2)',
                            ),
                            hoverinfo='skip',
                            visible=False,
                        ))
                        Nplot += 1

            # find cmin,cmax for traces and set to global min,max
            for i in range(len(fig.data)):
                if fig.data[i].marker.cmin is not None:
                    cmin=min(cmin,float(fig.data[i].marker.cmin))
                    cmax=max(cmax,float(fig.data[i].marker.cmax))
            for i in range(len(fig.data)):
                if fig.data[i].marker.cmin is not None:
                    fig.data[i].marker.cmin=cmin
                    fig.data[i].marker.cmax=cmax

            # find min,max locations to pad sphere points
            for i in range(len(fig.data)):
                xmin=min(xmin,np.amin(fig.data[i].x))
                xmax=max(xmax,np.amax(fig.data[i].x))
                ymin=min(ymin,np.amin(fig.data[i].y))
                ymax=max(ymax,np.amax(fig.data[i].y))
                zmin=min(zmin,np.amin(fig.data[i].z))
                zmax=max(zmax,np.amax(fig.data[i].z))
                
            # Add 1 RE sphere, padded to cover all data positions
            dataXYZ = pd.read_csv('https://ccmc.gsfc.nasa.gov/Kamodo/demo/sphereXYZ.csv')
            dataIJK = pd.read_csv('https://ccmc.gsfc.nasa.gov/Kamodo/demo/sphereIJK.csv')
            if scale == "km":
                dataXYZ=dataXYZ*self.RE
            fig.add_mesh3d()
            fig.data[Nplot].x = np.append(dataXYZ['x'],(xmin,xmax))
            fig.data[Nplot].y = np.append(dataXYZ['y'],(ymin,ymax))
            fig.data[Nplot].z = np.append(dataXYZ['z'],(zmin,zmax))
            fig.data[Nplot].i = dataIJK['i']
            fig.data[Nplot].j = dataIJK['j']
            fig.data[Nplot].k = dataIJK['k']
            fig.data[Nplot].facecolor = dataIJK['c']
            fig.data[Nplot].flatshading = True
            fig.data[Nplot].name = '1 R_E sphere'
            fig.data[Nplot].hovertemplate="Earth<extra></extra>"
            fig.data[Nplot].visible=True
            Nplot += 1

            # Start selection slider build
            steps = []
            for i in range(nsteps):
                step = dict(
                    method="update",
                    args=[{"visible": [False] * len(fig.data)}],
                    label=udates[i]
                )
                step["args"][0]["visible"][len(fig.data)-1] = True  # Set last trace to "visible"
                for j in range(len(self.satellites)):
                    if quiver:
                        step["args"][0]["visible"][i+(2*j*nsteps)] = True  # Toggle i'th trace
                        step["args"][0]["visible"][i+((2*j+1)*nsteps)] = True  # Toggle quiver trace
                    else:
                        step["args"][0]["visible"][i+(j*nsteps)] = True  # Toggle i'th trace
                steps.append(step)

            # Repeat first iteration of loop to set who is visible from start
            i=0
            fig.data[len(fig.data)-1].visible=True  # Set last trace to "visible"
            for j in range(len(self.satellites)):
                if quiver:
                    fig.data[i+(2*j*nsteps)].visible = True  # Toggle i'th trace
                    fig.data[i+((2*j+1)*nsteps)].visible = True  # Toggle quiver trace
                else:
                    fig.data[i+(j*nsteps)].visible = True  # Toggle i'th trace
            
            sliders = [dict(
                active=0,
                currentvalue={"prefix": "Currently showing: "},
                pad={"t": 50},
                steps=steps
            )]
            if nsteps > 1:
                fig.update_layout(sliders=sliders)

            # Update axis labels and plot title
            fig.update_layout(scene_aspectmode='data')
            fig.update_layout(scene=dict(xaxis=dict(title=dict(text="X ["+scale+"]")),
                                         yaxis=dict(title=dict(text="Y ["+scale+"]")),
                                         zaxis=dict(title=dict(text="Z ["+scale+"]"))))
            fig.update_layout(title_text=txttop)
            fig.update_layout(showlegend=False)
            fig.update_layout(scene_camera=dict(center=dict(x=0, y=0, z=0)))

            return fig
