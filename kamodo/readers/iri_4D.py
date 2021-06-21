from kamodo import Kamodo
from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta, timezone
import time as ti
import kamodo.readers.reader_plotutilities as RPlot
import kamodo.readers.reader_utilities as RU

#variable name in file: [standardized variable name, descriptive term, units]
model_varnames = {'Ne':['N_e','electron_density','4D','1/m**3'], 
                'Te':['T_e','electron_temperature','4D','K'],
                'Ti':['T_i','ion_temperature','4D','K'], 
                'Tn':['T_n','neutral_temperature','4D','K'],
                'O+':['N_Oplus','atomic_oxygen_ion_density','4D','1/m**3'],
                'H+':['N_Hplus','atomic_hydrogen_ion_density','4D','1/m**3'],
                'He+':['N_Heplus','atomic_helium_ion_density','4D','1/m**3'],
                'O2+':['N_O2plus','molecular_oxygen_ion_density','4D','1/m**3'],
                'NO+':['N_NOplus','nitric_oxide_ion_density','4D','1/m**3'],
                'N+':['N_Nplus','atomic_nitrogen_ion_density','4D','1/m**3'],
                'TEC':['TEC','total_electron_content','3D','10**16/m**2'],
                'NmF2':['NmF2','max_electron_density','3D','1/m**3'],
                'HmF2':['HmF2','max_electron_density_height','3D','km']}



#times from file converted to seconds since midnight of filedate
#plotting input times will be datetime strings of format 'YYYY-MM-DD HH:mm:ss'
#filedate is self.filedate from iri object
#converts to hours since midnight of filedate for plotting

class MODEL(Kamodo):
    def __init__(self, filename, variables_requested = None, runname = "noname",
                 printfiles=True, filetimes=False, gridded_int=True, **kwargs): #                 time_index=None, time_seconds=None,
        # Prepare model for function registration for the input argument
        super(MODEL, self).__init__(**kwargs)

        #collect filenames
        if '.2D.' in filename:  #require that input filename be for 3D file
            filename2d = filename
            f = filename.replace('.3D.', '.2D.')  #can't replace in place
            filename = f
        else: filename2d=filename.replace('.3D.','.2D.')
        self.filename = filename
        self.filename2d=filename2d
        if printfiles: print(filename,filename2d)
        
        #establish time attributes first
        self._iri3D = Dataset(filename, 'r')
        self._time = np.array(self._iri3D.variables['time'])/60.  #convert to hours since midnight of file        
        self.filedate = datetime(int(filename[-10:-6]),1,1,0,0,0).replace(tzinfo=timezone.utc)+\
            timedelta(days=int(filename[-6:-3])-1)
        #strings with timezone info chopped off (UTC anyway)
        self.datetimes=[(self.filedate+timedelta(hours=self._time[0])).isoformat(sep=' ')[:19], 
                        (self.filedate+timedelta(hours=self._time[-1])).isoformat(sep=' ')[:19]]  #strings
        self.filetimes=[datetime.timestamp(datetime.strptime(dt, '%Y-%m-%d %H:%M:%S').replace(\
            tzinfo=timezone.utc)) for dt in self.datetimes]   #timestamp in seconds, for value matching in wrapper?
        self.timerange0={'min':self.datetimes[0], 'max':self.datetimes[1],
                            'n':len(self._time)}     #strings in format = YYYY-MM-DD HH:MM:SS 
        self.timerange = self.timerange0
        if filetimes: 
            return
                
        #collect data and make dimensional grid from 3D file
        self._iri2D = Dataset(filename2d, 'r')
        self._lon = np.array(self._iri3D.variables['lon'])
        self._lat = np.array(self._iri3D.variables['lat'])
        self._height = np.array(self._iri3D.variables['ht'])
        
        #store a few items in iri object
        self.missing_value = np.NAN
        self._registered = 0
        self.variables={}
        self.runname=runname
        self.modelname = 'MODEL'
        
        #if variables_requested not given, collect all values from dict above as a list
        if variables_requested is None:
            variables_requested = [value[0] for key,value in model_varnames.items()]
            
        #collect list of iri variable name equivalents
        var_names = [key for key, value in model_varnames.items() if value[0] in variables_requested]
        extra_variables = [var for var in variables_requested if var not in 
                     [value[0] for key, value in model_varnames.items()]]
        if len(extra_variables)>0:   #pull out variables not allowed and error if not empty
            print('Some requested variables are not available:', extra_variables)        
        
        #register each variable desired 
        for varname in var_names:
            #determine source file type for variable
            file_type=''
            if varname in self._iri3D.variables.keys(): file_type = '3D'
            elif varname in self._iri2D.variables.keys(): file_type = '2D'
            else:
                raise AttributeError(f"{varname} not found in the files' metadata.")           
                
            #set variables, units
            variable = np.array(getattr(self, '_iri'+file_type).variables[varname])  #set data        
            if (len(variable.shape) not in [3,4]): continue  #skip anything not 3D or 4D
            units = model_varnames[varname][-1]  #units stored as last item in list per varname
            kamodo_varname = model_varnames[varname][0]
            
            #register allowed 3D and 4D variables
            self.variables[kamodo_varname] = dict(units = units, data = variable)  #register in object
            if len(variable.shape) == 4:  #define and register interpolators for each
                self.register_4D_variable(units, variable, 
                                          kamodo_varname, gridded_int)  #len(var.shape) instead of file_type
            elif len(variable.shape) == 3:
                self.register_3D_variable(units, variable, 
                                          kamodo_varname, gridded_int)
        
        #close netCDF4 files, initialize plotting variables
        self._iri3D.close()
        self._iri2D.close()
        self = RPlot.initialize_4D_plot(self)  #initialize 4D plotting variables  

    #define and register a 3D variable
    def register_3D_variable(self, units, variable, varname, gridded_int):
        """Registers a 3d interpolator with 3d signature"""
        
        #define and register the interpolators
        xvec_dependencies = {'time':'hr','lat':'deg','lon':'deg'}
        self = RU.regdef_3D_interpolators(self, units, variable, self._time, 
                                       self._lat, self._lon, varname, 
                                       xvec_dependencies, gridded_int)
        return 
    
    #define and register a 4D variable
    def register_4D_variable(self, units, variable, varname, gridded_int):
        """Registers a 4d interpolator with 4d signature"""
        
        #define and register the fast interpolator
        xvec_dependencies = {'time':'hr','height':'km','lat':'deg','lon':'deg'}
        self = RU.regdef_4D_interpolators(self, units, variable, self._time,
                                          self._height, self._lat, self._lon,
                                          varname, xvec_dependencies, gridded_int)
        return

    def set_plot(self, var, plottype, cutV=400., cutL=0, 
                 timerange={}, lonrange={}, latrange={}, htrange={}):
        '''Set plotting variables for available preset plot types.'''
        
        tic = ti.perf_counter()  #start timer
        test = RPlot.if_new_plot(self, var, plottype, cutV, cutL, timerange, 
                                 lonrange, latrange, htrange)
        if test==0: return
        else: self=test
        self = RU.setup_interpolating_grids(self, var) 
        toc = ti.perf_counter()  #end timer
        print(f"Time resetting plot and precomputing interpolations: {toc - tic:0.4f} seconds")
        return
    
    def get_plot(self, var, colorscale="Viridis",datascale="linear", ellipse=False):
        '''
        Return a plotly figure object for the available plot types set in set_plot()..
        colorscale = Viridis [default], Cividis, BlueRed or Rainbow
        '''
        
        self.gridSize=len(self._height)*len(self._lat)*len(self._lon)
        print(f'set_plot::colorscale={colorscale}, datascale={datascale}')
        print(f'Run: {self.runname}')
        #Set some text strings
        txtbot=f"Model: IRI, dt={self.dt} hrs, dlat={self.dlat} deg, "+\
            f"dlon={self.dlon} deg, dz={self.dz} km." #{self.gridSize} volume cells, 
        
        #set plot variables
        xint, yint, nx, ny, kamodo_plot, xlabel, ylabel, xformat, yformat, \
            zformat, xunit, yunit, txttop, txtbar, result = RPlot.set_plotvar(self, datascale, var)
        if ('TimeLon' in self.plottype or 'TimeLat' in self.plottype) and self.nDim==4: 
            plot_flip=True
        else:
            plot_flip=False
        
        #get plot with chosen settings
        fig = RPlot.heatplot2D(xint, yint, plot_flip, nx, ny, result, datascale, kamodo_plot, 
               xlabel, ylabel, colorscale, txtbar, xformat, yformat, zformat, 
               xunit, yunit, txttop, txtbot, ellipse=ellipse)

        return fig
    
    def make_plot(self, var, plottype, cutV=0, cutL=0, timerange={},
                 lonrange={}, latrange={}, htrange={}, 
                 colorscale="Viridis", datascale="linear", ellipse=False):
        '''Simplified call for plots. Execute with iplot(self.make_plot(....))
        Possible plottypes: LonLat, LonH, LatH, TimeLon, TimeLat, TimeH'''
        
        test = self.set_plot(var, plottype, cutV=cutV, cutL=cutL, timerange=timerange,
                 lonrange=lonrange, latrange=latrange, htrange=htrange)
        if test==1: return {} #if plottype requested invalid for variable, do nothing
        fig = self.get_plot(var, colorscale=colorscale, datascale=datascale, ellipse=ellipse)
        return fig