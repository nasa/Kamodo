'''
Kamodofication of the CTIPe model output
'''

import numpy as np
import time as ti
import glob, os
from kamodo import Kamodo
from netCDF4 import Dataset
from datetime import datetime, timezone
#import kamodo.readers.reader_plotutilities as RPlot
import kamodo.readers.reader_utilities as RU


# constants and dictionaries
ctipe_varnames = {'density':['rho','kg/m**3'],
                  'temperature':['T','K'],
                  'electron_temperature':['T_e','K'],
                  'ion_temperature':['T_i','K'],
                  'height':['H_ilev','m'],                      
                  'meridional_neutral_wind':['Vn_lat','m/s'],
                  'zonal_neutral_wind':['Vn_lon','m/s'],
                  'vertical_neutral_wind':['Vn_H','m/s'],
                  'neutral_temperature':['T_n','K'],
                  'mean_molecular_mass':['Rmt','amu'],
                  'electron_density':['N_e','1/m**3'],
                  'neutral_density':['N_n','1/m**3'],
                  'solar_heating':['Q_Solar','J/kg/s'],
                  'joule_heating':['Q_Joule','J/kg/s'],
                  'radiation_heat_cool':['Q_radiation','J/kg/s'],
                  'atomic_oxygen_density':['N_O','1/m**3'],
                  'molecular_oxygen_density':['N_O2','1/m**3'],
                  'molecular_nitrogen_density':['N_N2','1/m**3'],
                  'nitric_oxide_density':['N_NO','1/m**3'],
                  'nitric_oxide_ion_density':['N_NOplus','1/m**3'],
                  'molecular_nitrogen_ion_density':['N_N2plus','1/m**3'],  
                  'molecular_oxygen_ion_density':['N_O2plus','1/m**3'],
                  'atomic_nitrogen_ion_density':['N_Nplus','1/m**3'],
                  'atomic_oxygen_ion_density':['N_Oplus','1/m**3'],
                  'atomic_hydrogen_ion_density':['N_Hplus','1/m**3'],
                  'pedersen_conductivity':['Sigma_P','S/m'],
                  'hall_conductivity':['Sigma_H','S/m'],
                  'zonal_ion_velocity':['Vi_lon','m/s'],
                  'meridional_ion_velocity':['Vi_lat','m/s'],
                  'height_integrated_joule_heating':['W_Joule','W/m**2'],
                  'energy_influx':['Eflux_precip','W/m**2'],
                  'mean_energy':['Eavg_precip','keV'],
                  'total_electron_content':['TEC','10**16/m**2'],
                  'theta_electric_field_at_140km':['E_theta140km','V/m'],
                  'lambda_electric_field_at_140km':['E_lambda140km','V/m'],
                  'theta_electric_field_at_300km':['E_theta300km','V/m'],
                  'lambda_electric_field_at_300km':['E_lambda300km','V/m']}


#convert an array of timestamps to an array of hrs since midnight
@np.vectorize
def ts_to_hrs(time_val, filedate):
    return (datetime.utcfromtimestamp(time_val).replace(tzinfo=timezone.utc)-filedate).total_seconds()/3600.

def CTIPe_filesearch(filename):
    '''check for wrapped data output in dir, create if necessary'''
    
    if 'wrapped' not in filename:
        print('Wrapped file not given, searching for file...')
        file_pattern = filename.split('plot')[0]  #.../CTIPe/2015-03-18- to choose this group of files
        files = glob.glob(file_pattern+'plot-*-wrapped.nc')
        if len(files)==3: 
            print('Found file.')
            filename = files[0]  #all three files exist, use them
        else: 
            print('Files not found. Generating...')
            from kamodo.readers.ctipe_data_wrapper import ctipe_wrap_files as wrap
            filename = wrap(filename)  #produce the wrapped files, returns the new filename
        return filename
    elif not os.path.isfile(filename): #check to make sure the file exists
        print('Files not found. Generating...')
        from kamodo.readers.ctipe_data_wrapper import ctipe_wrap_files as wrap
        filename = wrap(filename.split('-wrapped')[0]+'.nc')  #produce the wrapped files, returns the new filename   
        return filename
    elif os.path.isfile(filename):  #if file exists, return filename
        return filename
    
#main class
class CTIPe(Kamodo):
    def __init__(self, filename, variables_requested = None, filetimes=False,
                 runname = "noname", printfiles=True, gridded_int=True, **kwargs):  
                #date = None, date is in filename, so exclude? (self.date ....)
                #time=None, runpath = "./", not used

        # input file name can be one of the 4 files for each day of model outputs
        # YYYYMMDD-plot-[density|height|neutral|plasma].nc files
        # only the density, height and neutral files have data and are read
        super(CTIPe, self).__init__()   #what does this line do??
        
        filename = CTIPe_filesearch(filename)  #get/convert filename
        filetype_list = ['plot-density-wrapped','plot-height-wrapped',
                         'plot-neutral-wrapped','plot-plasma-wrapped']
        file_beg, file_end = [filename.split(filetype) for filetype in filetype_list 
                              if len(filename.split(filetype))>1][0]
        filename_density, filename_height, filename_neutral = [
            file_beg+filetype+file_end for filetype in filetype_list[:-1]]
        self.filename = [filename_density, filename_height, filename_neutral]
        if printfiles: 
            print('Files:\n'+filename_density+'\n'+filename_height+'\n'+filename_neutral)    
            
        #establish time attributes first
        self._ctipe_density = Dataset(filename_density)
        self.filedate = datetime.strptime(file_beg[-11:-1]+' 00:00:00', 
                                          '%Y-%m-%d %H:%M:%S').replace(tzinfo=timezone.utc)
        t = np.array(self._ctipe_density.variables['time'])
        self.datetimes=[datetime.utcfromtimestamp(t[0]).isoformat(sep=' '), 
                        datetime.utcfromtimestamp(t[-1]).isoformat(sep=' ')]  #strings
        self.filetimes=[t[0], t[-1]]   #timestamps in hours for matching in wrapper
        self.timerange0={'min':self.datetimes[0], 'max':self.datetimes[1],
                            'n':len(t)}     #strings in format = YYYY-MM-DD HH:MM:SS 
        self.timerange = self.timerange0
        if filetimes: 
            return

        #pull in remaining datasets from files into kamodo object
        self.modelname = 'CTIPe'
        self._ctipe_height = Dataset(filename_height)  #in meters
        self._ctipe_neutral = Dataset(filename_neutral)
        
        #pull in coordinate variables into kamodo object (density file)
        self._time = np.array(ts_to_hrs(t, self.filedate))  #convert timestamps to hrs since midnight 
        self._ilev = np.array(self._ctipe_density.variables['plev'])  #_ilev to match others
        self._lat = np.array(self._ctipe_density.variables['lat'])
        self._lon = np.array(self._ctipe_density.variables['lon'])

        #pull in coordinate variables into kamodo object (height file, in km)
        self._height = np.array(self._ctipe_height.variables['ht'])  #_height to match others
        self._lat_height = np.array(self._ctipe_height.variables['lat'])
        self._lon_height = np.array(self._ctipe_height.variables['lon'])

        ##pull in coordinate variables into kamodo object (neutral file)
        self._ilev_neutral = np.array(self._ctipe_neutral.variables['plev'])
        self._lat_neutral = np.array(self._ctipe_neutral.variables['lat'])
        self._lon_neutral = np.array(self._ctipe_neutral.variables['lon'])
        self._elat = np.array(self._ctipe_neutral.variables['elat'])
        self._elon = np.array(self._ctipe_neutral.variables['elon'])

        #initialize variables
        self._registered = 0  
        super(CTIPe, self).__init__()   #what does this line do???
        self.filename = filename
        self.runname = runname
        self.missing_value = np.NAN
        self.variables=dict()
        
        #if variables_requested not given, collect all values from dict above as a list
        if variables_requested is None:
            variables_requested = [value[0] for key,value in ctipe_varnames.items()]
        
        # add height variable needed to height (not IP-level) interpolatioms
        if 'H_ilev' not in variables_requested: variables_requested.append('H_ilev')
        #print(f'Requested {len(variables_requested)} variables: {variables_requested} \n')
        
        #collect list of ctipe variable name equivalents
        var_names = [key for key, value in ctipe_varnames.items() if value[0] in variables_requested]
        extra_variables = [var for var in variables_requested if var not in 
                     [value[0] for key, value in ctipe_varnames.items()]]
        if len(extra_variables)>0:   #pull out variables not allowed and error if not empty
            raise AttributeError("No such variable(s):{}".format(extra_variables))

        #cycle through all variables in one loop
        bad_varnames={'density': ['ZMAG','Rmt','H'], 'height': ['ZMAG'],
                  'neutral': ['ZMAG','electron_density']}   #initialize bad_var dictionary per file_type
        for varname in var_names:
            #determine source file type for variable
            file_type=''
            if varname in self._ctipe_density.variables.keys(): file_type = 'density'
            elif varname in self._ctipe_height.variables.keys(): file_type = 'height'
            elif varname in self._ctipe_neutral.variables.keys(): file_type = 'neutral'
            else:
                raise AttributeError(f"{varname} not found in the files' metadata.")  
            
            #set units, initialize variables
            variable = np.array(getattr(self, '_ctipe_'+file_type).variables[varname])  #set variable
            units = ctipe_varnames[varname][-1]
            if (len(variable.shape) not in [3,4]) or (varname in bad_varnames[file_type]):
                continue  #if not 3D or 4D or not allowed, skip to next variable
            
            #register allowed 3D and 4D variables
            kamodo_varname=ctipe_varnames[varname][0]  #retreive standardized name
            self.variables[kamodo_varname] = dict(units = units, data = variable)  #register in object
            if len(variable.shape) == 4:  #define and register interpolators for each
                self.register_4D_variable(units, variable, kamodo_varname, 
                                          file_type, gridded_int)
            elif len(variable.shape) == 3:
                self.register_3D_variable(units, variable, kamodo_varname, 
                                          file_type, gridded_int)
        
        #close netCDF4 files, initialize plotting variables
        self._ctipe_density.close()
        self._ctipe_height.close()
        self._ctipe_neutral.close()
        self = RPlot.initialize_4D_plot(self)  #initialize 4D plotting variables         

    #define and register a 3D variable
    def register_3D_variable(self, units, variable, varname, file_type, gridded_int):
        """Registers a 3d interpolator with 3d signature"""
        
        #determine coordinate variables by file_type
        if file_type=='density': 
            t, lat, lon = self._time, self._lat, self._lon
            xvec_dependencies = {'time':'hr','lat':'deg','lon':'deg'}
        if file_type=='height': 
            t, lat, lon = self._time, self._lat_height, self._lon_height
            xvec_dependencies = {'time':'hr','lat':'deg','lon':'deg'}
        if file_type=='neutral':
            if variable.shape[1] == self._elat.shape[0]:
                t, lat, lon = self._time,self._elat, self._elon  
                xvec_dependencies = {'time':'hr','elat':'deg','elon':'deg'}
            else:
                t, lat, lon = self._time, self._lat_neutral, self._lon_neutral  
                xvec_dependencies = {'time':'hr','lat':'deg','lon':'deg'}        
        
        #define and register the interpolators
        self = RU.regdef_3D_interpolators(self, units, variable, t, lat, lon, 
                                          varname, xvec_dependencies, gridded_int)
        return 
    
    #define and register a 4D variable
    def register_4D_variable(self, units, variable, varname, file_type, gridded_int):
        """Registers a 4d interpolator with 4d signature"""
        
        #determine coordinate variables by file_type
        if file_type=='density': 
            t, z, lat, lon = self._time, self._ilev, self._lat, self._lon
            xvec_dependencies = {'time':'hr','ilev':'m/m','lat':'deg','lon':'deg'}
        if file_type=='height':
            t, z, lat, lon = self._time, self._height, self._lat_height, self._lon_height
            xvec_dependencies = {'time':'hr','height':'km','lat':'deg','lon':'deg'}
        if file_type=='neutral':
            t, z, lat, lon= self._time, self._ilev_neutral, self._lat_neutral, self._lon_neutral
            xvec_dependencies = {'time':'hr','ilev':'m/m','lat':'deg','lon':'deg'}
        
        #define and register the interpolators
        self = RU.regdef_4D_interpolators(self, units, variable, t, z, lat, lon,
                                          varname, xvec_dependencies, gridded_int)
        return
    
    """----------------------- Plotting code below here --------------------

    def set_plot(self, var, plottype, cutV=10, cutL=0, timerange={},
                 lonrange={}, latrange={}, htrange={}):
        '''Set plotting variables for available preset plot types.'''
        
        #compare plot data to defaults
        tic = ti.perf_counter()  #start timer
        test = RPlot.if_new_plot(self, var, plottype, cutV, cutL, timerange, 
                                 lonrange, latrange, htrange)
        if test in [0,1]: return test  #do nothing if plot variables are unchanged or type invalid
        else: self=test
        self = RU.setup_interpolating_grids(self, var) 
        toc = ti.perf_counter()  #end timer
        print(f"Time resetting plot and precomputing interpolations: {toc - tic:0.4f} seconds")
        return 2
    
    def get_plot(self, var, colorscale="Viridis",datascale="linear", ellipse=False):
        '''
        Return a plotly figure object for the available plot types set in set_plot()..
        colorscale = Viridis [default], Cividis, BlueRed or Rainbow
        '''
        
        #need to determine elat/elon vs lat/lon dependence here with dependencies
        arg_dict = self.variables[var]['xvec']
        if 'elon' in arg_dict.keys(): 
            self.dLON, self.nLON = self.delon, len(self._elon)
        else: 
            self.dLON, self.nLON = self.dlon, len(self._lon)
        if 'elat' in arg_dict.keys(): 
            self.dLAT, self.nLAT = self.delat, len(self._elat)  #electric field latitude
        else: 
            self.dLAT, self.nLAT = self.dlat, len(self._lat)
        self.gridSize = self.nLAT*self.nLON*self.numZ
        
        print(f'set_plot::colorscale={colorscale}, datascale={datascale}')
        print(f'Runname: {self.runname}')
        #Set some text strings  
        txtbot=f"Model=CTIPe, dt={self.dt} hrs, dlat={self.dLAT} deg, "+\
            f"dlon={self.dLON} deg, dz={self.dz} {self.dzunit}" #{self.gridSize} volume cells,
        
        #set plot variables
        xint, yint, nx, ny, kamodo_plot, xlabel, ylabel, xformat, yformat, \
            zformat, xunit, yunit, txttop, txtbar, result = RPlot.set_plotvar(self, datascale, var)
        if 'TimeLon' in self.plottype: 
            plot_flip = True
        elif 'TimeLat' in self.plottype and self.nDim==4:
            plot_flip = True
        else:
            plot_flip = False
            
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
    """
