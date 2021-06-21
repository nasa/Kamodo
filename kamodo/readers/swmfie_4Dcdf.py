# -*- coding: utf-8 -*-
"""
Created on Mon May 17 18:53:35 2021

@author: rringuet
"""
from os import path
import numpy as np
import time as ti
from datetime import datetime, timezone
from netCDF4 import Dataset
from kamodo import Kamodo
import kamodo.readers.reader_plotutilities as RPlot
import kamodo.readers.reader_utilities as RU
#read 1 day of data from cdf instead of from multiple .tec files


model_varnames={"X":['x','1D','km'],"Y":['y','1D','km'],"Z":['z','1D','km'],    #(ignored, given in R_E on a unit sphere)
                 "Theta":['theta','1D',"deg"],"Psi":['psi','1D',"deg"],     #(used as coordinates)
                 "Btilt_theta":['theta_Btilt','1D',"deg"], "Btilt_psi":['psi_Btilt','1D',"deg"], 
                 #(added directly to object for documentation purposes)
                 "SigmaH":['Sigma_H','3D',"S"],"SigmaP":['Sigma_P','3D',"S"],
                 "E-Flux":['Phi_E','3D',"W/m**2"], "Ave-E":['E_avg','3D','eV'],
                 "JR":["j_R",'3D',"mA/m**2"],"PHI":["Phi",'3D',"kV"],
                 "Ex":["E_x",'3D',"mV/m"],"Ey":["E_y",'3D',"mV/m"],"Ez":["E_z",'3D',"mV/m"],
                 "Jx":["j_x",'3D',"mA/m**2"],"Jy":["j_y",'3D',"mA/m**2"],"Jz":["j_z",'3D',"mA/m**2"],
                 "Ux":['v_x','3D',"km/s"],"Uy":['v_y','3D',"km/s"],"Uz":['v_z','3D',"km/s"],
                 "JouleHeat":['Q_Joule','3D',"mW/m**2"], "IonNumFlux":['Phi_nion','3D',"1/cm**2/s"],
                 "RT 1/B":['Binv_RT','3D',"1/T"],"RT Rho":['rho_RT','3D',"amu/cm**3"],"RT P":['P_RT','3D',"Pa"],
                 "conjugate dLat":['dLat_star','3D',"deg"],"conjugate dLon":['dlon_star','3D',"deg"]}
                 
 
'''                
#documentation variables:
"X":['x','km'],"Y":['y','km'],"Z":['z','km'],    (ignored, given in R_E on a unit sphere)
"Theta":['theta',"deg"],"Psi":['psi',"deg"],     (used as coordinates)
"Btilt_theta":['theta_Btilt',"deg"], "Btilt_psi":['psi_Btilt',"deg"] 
(added directly to object for documentation purposes)
'''
@np.vectorize
def dts_to_hrs(datetime_string, filedate):
    '''Get hours since midnight from datetime string'''
    
    return (datetime.strptime(datetime_string, '%Y-%m-%d %H:%M:%S').replace(tzinfo=timezone.utc)\
            -filedate).total_seconds()/3600.

@np.vectorize
def filename_to_dts(filename, string_date):
    '''Get datetime string in format "YYYY-MM-SS HH:mm:SS" from filename'''
    
    mmhhss = filename.split('/')[-1].split('\\')[-1][12:18]
    return string_date+' '+mmhhss[:2]+':'+mmhhss[2:4]+':'+mmhhss[4:] 

def dts_to_ts(file_dts):
    '''Get datetime timestamp in UTC from datetime string'''    
    
    return datetime.timestamp(datetime.strptime(file_dts, '%Y-%m-%d %H:%M:%S'
                                                ).replace(tzinfo=timezone.utc))
    
#base for main class object----------------------------------------
#main class object
#file_dir = 'C:/Users/rringuet/Kamodo_WinDev1/SWMF_IE/'
#file_prefix = 'C:/Users/rringuet/Kamodo_WinDev1/SWMF_IE/i_e20180801'  #example
#files = glob.glob(file_dir+'i_e*.tec')  #for wrapper, this and the next line
#file_patterns = np.unique([file_dir+f.split('/')[-1].split('\\')[-1][:11] for f in files])

       
class MODEL(Kamodo): 
    def __init__(self, file_prefix, variables_requested=[], runname="noname",
                 filetimes=False, verbose=False, gridded_int=True, printfiles=True,
                 **kwargs): 
        '''file_prefix must be of form "3D***_tYYMMDD" to load all files for one day
         and include a complete path to the files'''
        super(MODEL, self).__init__()
        
        #check if given .nc file exists. If not, convert files with same prefix to netCDF
        if not path.isfile(file_prefix+'.nc'):
            from kamodo.readers.swmfie_tocdf import convertSWMFIE_toCDF
            test = convertSWMFIE_toCDF(file_prefix)  
            if not test: 
                self.conversion_test = test  #only needed for 1 file/time cases
                return    #if file conversion fails, return 0
        t0 = ti.perf_counter()
        
        #establish time attributes first for file searching
        file_datestr = file_prefix.split('/')[-1].split('\\')[-1][3:11]
        string_date = file_datestr[:4]+'-'+file_datestr[4:6]+'-'+file_datestr[6:8]  #'YYYY-MM-DD'
        self.filedate = datetime.strptime(string_date+' 00:00:00', \
                                          '%Y-%m-%d %H:%M:%S').replace(tzinfo=timezone.utc) #dt object
        
        #establish beginning and end time of file list
        cdf_data = Dataset(file_prefix+'.nc', 'r')
        files = cdf_data.file.split(',')
        self.datetimes = list(filename_to_dts([files[0], files[-1]], string_date))  #strings in format = YYYY-MM-DD HH:MM:SS 
        self.timerange0={'min':self.datetimes[0], 'max':self.datetimes[-1], 'n':len(files)}
        self.timerange = self.timerange0
        self.filetimes=[dts_to_ts(file_dts) for file_dts in self.datetimes]   #timestamps in UTC     
        if filetimes: 
            return           
    
        #return if only one file found because interpolator code will break
        if len(files)<2:
            print('Not enough files found with given file prefix.')
            return
        
        #store variables
        self.filename = files
        self.runname = runname
        self.missing_value = np.NAN
        self._registered = 0
        self.variables=dict()
        if printfiles: 
            print('Files:')
            for file in self.filename: print(file)

        #get list of variables possible in these files using first file
        if len(variables_requested)>0:
            gvar_list = [key for key in cdf_data.variables.keys() if key \
                         in variables_requested]
            if len(gvar_list)!=len(variables_requested):
                err_list = [item for item in variables_requested if item not in \
                            cdf_data.variables.keys()]
                print('Some requested variables are not available:', err_list)
        else:
            gvar_list = [key for key in cdf_data.variables.keys() \
                         if key not in cdf_data.dimensions.keys() and key not in \
                             ['theta_Btilt', 'psi_Btilt']]  
            #avoid returning coordinates stored elsewhere (or ignored)
        #print(gvar_list)
        
        #store coordinate data and Btilt (for documentation)
        self._time = np.array(cdf_data.variables['time'])  #hours since midnight
        self._height = np.array(cdf_data.variables['height'])
        self._lat = np.array(cdf_data.variables['lat'])
        self._lon = np.array(cdf_data.variables['lon'])
        self.theta_Btilt = np.array(cdf_data.variables['theta_Btilt'])
        self.psi_Btilt = np.array(cdf_data.variables['psi_Btilt'])

        # Store variable's data, units, and datatypes.
        self.variables = {key:{'units':cdf_data.variables[key].units, 'dtype':np.float32,
                               'data':np.array(cdf_data.variables[key])}\
                          for key in gvar_list}
        cdf_data.close()
        if verbose: print(f'Took {ti.perf_counter()-t0:.6f}s to read in data')
            
        #register interpolators for each variable
        varname_list = [key for key in self.variables.keys()]  #store original list b/c gridded interpolators
        t_reg = ti.perf_counter()
        for varname in varname_list:  #all are 3D variables
            self.register_3D_variable(self.variables[varname]['units'], 
                                      self.variables[varname]['data'], varname,
                                      gridded_int)
        if verbose: print(f'Took {ti.perf_counter()-t_reg:.5f}s to register '+\
                          f'{len(varname_list)} variables.')
        self = RPlot.initialize_4D_plot(self)  #initialize   
        if verbose: print(f'Took a total of {ti.perf_counter()-t0:.5f}s to kamodofy '+\
                          f'{len(gvar_list)} variables.')

    #define and register a 3D variable-----------------------------------------
    def register_3D_variable(self, units, variable, varname, gridded_int):
        """Registers a 3d interpolator with 3d signature"""
        
        #define and register the interpolators
        xvec_dependencies = {'time':'hr','lat':'deg','lon':'deg'}
        self = RU.regdef_3D_interpolators(self, units, variable, self._time, 
                                       self._lat, self._lon, varname, 
                                       xvec_dependencies, gridded_int)       
        return 

#begin plotting code -----------------------------------
    def set_plot(self, var, plottype, cutV=400., cutL=0, 
                 timerange={}, lonrange={}, latrange={}, htrange={}):
        '''Set plotting variables for available preset plot types.'''
        
        tic = ti.perf_counter()  #start timer
        test = RPlot.if_new_plot(self, var, plottype, cutV, cutL, timerange, 
                                 lonrange, latrange, htrange)
        if test==0: return
        else: self=test
        self = RU.setup_interpolating_grids(self,var)
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
        txtbot=f"Model: SWMF_IE, dt={self.dt:.4f} hrs, "+\
            f"dlat={self.dlat:.1f} deg, dlon={self.dlon:.1f} deg, dz={self.dz:.1f} km."
        #{self.gridSize} volume cells,
            
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
