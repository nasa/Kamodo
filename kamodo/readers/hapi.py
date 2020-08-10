from kamodo import Kamodo, kamodofy
import numpy as np
import scipy
import datetime
import urllib, json
# Renamed hapi as PYhapi to avoid name conflict
from hapiclient import hapi as PYhapi
from hapiclient import hapitime2datetime
import plotly.express as px
import plotly.graph_objects as go

def hapi_get_info(server, dataset):
    '''Query for info json return from HAPI server'''
    query = '{}/info?id={}'.format(server, dataset)
    response = urllib.request.urlopen(query)
    data = json.loads(response.read())
    return data

def hapi_get_date_range(server, dataset):
    '''Query for start and end date for dataset'''
    info = hapi_get_info(server, dataset)
    start_date = info['startDate']
    end_date = info['stopDate']
    return start_date, end_date

class HAPI(Kamodo):
    def __init__(self, server, dataset, parameters = None, start = None, stop = None, **kwargs):
        self.verbose=False
        self.symbol_registry=dict()
        self.signatures=dict()

        self.server = server
        self.dataset = dataset
        print(' -server: ',server)
        print(' -dataset: ',dataset)
        opts       = {'logging': False, 'usecache': False}
        if start is None or stop is None:
            start, stop = hapi_get_date_range(server, dataset)
            ts = datetime.datetime.strptime(stop, '%Y-%m-%dT%H:%M:%S.%fZ').timestamp()
            ts = ts - 3600.
            start = datetime.datetime.fromtimestamp(ts).strftime("%Y-%m-%dT%H:%M:%S.%fZ")
            print(" -dates not selected, using last hour available.")
            print("       start: ",start)
            print("         end: ",stop)
        else:
            print(" -date start: ",start)
            print("         end: ",stop)
        self.start=start
        self.stop=stop
        # Get data, meta from the HAPI python client
        data, meta = PYhapi(server, dataset, parameters, start, stop, **opts)
        self.data = data
        self.meta = meta
        self.startDate = meta['startDate']
        self.stopDate = meta['stopDate']
        # Convert binary time array into datetime
        self.dtarray = hapitime2datetime(data['Time'])
        self.tsarray = np.array([d.timestamp() for d in self.dtarray])
        self.variables=dict()
        for d in self.meta['parameters']:
            var = d['name']
            if var != "Time":
                self.variables[var] = dict(units=d['units'], data=self.data[var])

        for varname in self.variables:
            units = self.variables[varname]['units']
            print('... registering ',varname,units)
            self.register_variable(varname, units)
        
        self.possible_coords=('TOD','J2K','GEO','GM','GSM','GSE','SM')
        self.coords=dict()
        for varname in self.variables:
            tmp = varname.split("_")
            dir = tmp[0].lower()
            key = tmp[1]
            if key in self.possible_coords:
                if key not in self.coords:
                    self.coords[key] = dict(coord=key)
                self.coords[key][dir] = varname
 
    def register_variable(self, varname, units):
        """register variables into Kamodo for this service, HAPI"""

        def interpolate(timestamp):  
            data =  self.variables[varname]['data']
            return np.interp(timestamp,self.tsarray,data)

        # store the interpolator
        self.variables[varname]['interpolator'] = interpolate

        # update docstring for this variable
        interpolate.__doc__ = "A function that returns {} in [{}].".format(varname,units)

        self[varname] = kamodofy(interpolate, 
                                 units = units, 
                                 citation = "De Zeeuw 2020",
                                 data = None)

    def get_plot(self, coord, type="3D", alt="R_E"):
        '''
        Return a plotly figure object for selected coordinate system.
        type = 3D, time
        color = km, R_E
        '''
        txttop="HAPI data for " + self.dataset + " in "+ coord + " coordinates.<br>"\
            + self.start + " to " + self.stop
        
        if coord not in self.coords:
            print('ERROR, sent in coords not in dataset.')
            return
        
        if type == "3D":
            x=self.variables[self.coords[coord]['x']]['data']
            y=self.variables[self.coords[coord]['y']]['data']
            z=self.variables[self.coords[coord]['z']]['data']
            r=(np.sqrt(x**2 + y**2 + z**2))-1.
            if alt == "km":
                x=x*6.3781E3
                y=y*6.3781E3
                z=z*6.3781E3
                r=(np.sqrt(x**2 + y**2 + z**2))-6.3781E3
            fig=px.scatter_3d(
                x=x,
                y=y,
                z=z,
                color=r)
            bartitle = "Altitude [" + alt + "]"
            fig.update_layout(coloraxis=dict(colorbar=dict(title=bartitle)))
            fig.update_layout(scene=dict(xaxis=dict(title=dict(text="X ["+alt+"]")),
                                         yaxis=dict(title=dict(text="Y ["+alt+"]")),
                                         zaxis=dict(title=dict(text="Z ["+alt+"]"))))
            fig.update_layout(title_text=txttop)
            return fig
        
        if type == 'time':
            x=self.variables[self.coords[coord]['x']]['data']
            y=self.variables[self.coords[coord]['y']]['data']
            z=self.variables[self.coords[coord]['z']]['data']
            ytitle="Position [R_E]"
            if alt == "km":
                x=x*6.3781E3
                y=y*6.3781E3
                z=z*6.3781E3
                ytitle="Position [km]"
            fig=go.Figure()
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
        

        print('ERROR, reached end of get_plot without any action taken.')
        return

