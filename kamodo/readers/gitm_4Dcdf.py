from os import path
import time as ti
import glob
import numpy as np
from netCDF4 import Dataset
from datetime import datetime, timezone
from kamodo import Kamodo
import kamodo.readers.reader_plotutilities as RPlot
import kamodo.readers.reader_utilities as RU


#file_prefix = 'C:/Users/rringuet/Kamodo_WinDev1/GITM/3DLST_t150317'
#keep for easy access of variable names from satellite flythrough software
#line for gitm wrapper if data scale isn't important:
# return {key:[value[0],value[2]] for key, value in gitm_varnames.items()}
#varnames in cdf files are standardized (value[0])
gitm_varnames = {'Argon Mixing Ratio': ['r_Ar', 'linear', ''],
                 'Ar Mass Density': ['rho_Ar', 'exponential', 'kg/m**3'],
                 'Methane Mixing Ratio': ['r_CH4', 'linear', ''],
                 'Conduction': ['k', 'linear', 'W/m/K'],
                 'EUV Heating': ['Q_EUV', 'linear', 'K per timestep'],
                 'H Mass Density': ['rho_H', 'exponential', 'kg/m**3'],
                 'H+ Mass Density': ['rho_Hplus', 'exponential', 'kg/m**3'],
                 'H2 Mixing Ratio': ['r_H2', 'linear', ''],
                 'Hydrogen Cyanide Mixing Ratio': ['r_HCN', 'linear', ''],
                 'He Mass Density': ['rho_He', 'exponential', 'kg/m**3'],
                 'He+ Mass Density': ['rho_Heplus', 'exponential', 'kg/m**3'],
                 'Heating Efficiency': ['HeatingEfficiency', 'linear', ''],
                 'Heat Balance Total': ['HeatBalanceTotal', 'linear', ''],
                 'N2 Mass Density': ['rho_N2', 'exponential', 'kg/m**3'],
                 'N2+ Mass Density': ['rho_N2plus', 'exponential', 'kg/m**3'],
                 'N+ Mass Density': ['rho_Nplus', 'exponential', 'kg/m**3'],
                 'N(2D) Mass Density': ['rho_N2D', 'exponential', 'kg/m**3'],
                 'N(2P) Mass Density': ['rho_N2P', 'exponential', 'kg/m**3'],
                 'N(4S) Mass Density': ['rho_N4S', 'exponential', 'kg/m**3'],
                 'N2 Mixing Ratio': ['r_N2', 'linear', ''],
                 'NO Mass Density': ['rho_NO', 'exponential', 'kg/m**3'],
                 'NO+ Mass Density': ['rho_NOplus', 'exponential', 'kg/m**3'],
                 'O2 Mass Density': ['rho_O2', 'exponential', 'kg/m**3'],
                 'O(1D) Mass Density': ['rho_O1D', 'exponential', 'kg/m**3'],
                 'O2+ Mass Density': ['rho_O2plus', 'exponential', 'kg/m**3'],
                 'O(2D) Mass Density': ['rho_O2D', 'exponential', '1/m**3'],
                 'O+(2P) Mass Density': ['rho_Oplus2P', 'exponential', 'kg/m**3'],
                 'O(3P) Mass Density': ['rho_O3P', 'exponential', 'kg/m**3'],
                 'O+(4SP) Mass Density': ['rho_Oplus4SP', 'exponential', 'kg/m**3'],
                 'Radiative Cooling': ['L_Rad', 'linear', ''],
                 'Neutral Density': ['rho', 'exponential', 'kg/m**3'],
                 'T_n': ['T_n', 'linear', 'K'],
                 'vi_east': ['vi_east', 'linear', 'm/s'],
                 'vi_north': ['vi_north', 'linear', 'm/s'],
                 'vi_up': ['vi_up', 'linear', 'm/s'],
                 'vn_east': ['vn_east', 'linear', 'm/s'],
                 'vn_north': ['vn_north', 'linear', 'm/s'],
                 'vn_up': ['vn_up', 'linear', 'm/s'],
                 'v_N2_up': ['v_N2_up', 'linear', 'm/s'],
                 'v_N(4S)_up': ['v_N4S_up', 'linear', 'm/s'],
                 'v_N_up': ['v_N_up', 'linear', 'm/s'],
                 'v_O2_up': ['v_O2_up', 'linear', 'm/s'],
                 'v_O(3P)_up': ['v_O3P_up', 'linear', 'm/s'],
                 'v_He_up': ['v_He_up', 'linear', 'm/s'],
                 '[e-]': ['rho_e', 'linear', '1/m**3'],
                 'Electron Average Energy': ['ElectronAverageEnergy', 'linear', 'J'],
                 'T_e': ['T_e', 'linear', 'K'],
                 'T_i': ['T_i', 'linear', 'K'],
                 'Solar Zenith Angle': ['SolarZenithAngle', 'linear', 'radians'],
                 'CO2 Mass Density': ['rho_CO2', 'exponential', 'kg/m**3'],
                 'DivJu FL': ['DivJuFL', '', ''],
                 'DivJuAlt': ['DivJuAlt', 'linear', ''],
                 'Electron Energy Flux': ['ElectronEnergyFlux', 'exponential', 'J/m**2'],
                 'Field Line Length': ['Field Line Length', 'linear', 'm'],
                 'sigma_P': ['sigma_P', 'linear', 'S/m'],
                 'Sigma_P': ['Sigma_P', 'linear', 'S/m'],
                 'sigma_H': ['sigma_H', 'linear', 'S/m'],
                 'Potential': ['V', 'linear', 'V'],
                 'Sigma_H': ['Sigma_H', 'linear', 'S/m'],
                 'Region 2 Current': ['I_R2', 'linear', 'A/m**2'],
                 'Region 1 Current': ['I_R1', 'linear', 'A/m**2'],
                 'Ed1': ['Ed1', 'linear', ''],
                 'Ed2': ['Ed2', 'linear', ''],
                 'Solar Local Time': ['SolarLocalTime', 'linear', 'h'],
                 'Vertical Electric Field': ['E_up', 'linear', 'V/m'],
                 'Eastward Electric Field': ['E_east', 'linear', 'V/m'],
                 'Northward Electric Field': ['E_north', 'linear', 'V/m'],
                 'Electric Field Magnitude': ['E_mag', 'linear', 'V/m'],
                 'Vertical Magnetic Field': ['B_up', 'linear', 'nT'],
                 'Eastward Magnetic Field': ['B_east', 'linear', 'nT'],
                 'Northward Magnetic Field': ['B_north', 'linear', 'nT'],
                 'Magnetic Field Magnitude': ['B_mag', 'linear', 'nT'],
                 'Magnetic Latitude': ['MagLat', 'linear', 'deg'],
                 'Magnetic Longitude': ['MagLon', 'linear', 'deg'],
                 'g': ['g', 'linear', 'm/s**2'],
                 'GradP_east (P_i + P_e)': ['GradP_east', 'linear', 'Pa/m'],
                 'GradP_north (P_i + P_e)': ['GradP_north', 'linear', 'Pa/m'],
                 'GradP_up (P_i + P_e)': ['GradP_up', 'linear', 'Pa/m'],
                 'nu_in': ['nu_in', 'linear', '1/s'],
                 'Chemical Heating Rate': ['ChemicalHeatingRate', 'linear', ''],
                 'Total Absolute EUV': ['TotalAbsoluteEUV', 'linear', 'K per timestep'],
                 'O Cooling': ['L_O', 'linear', 'K per timestep'],
                 'Joule Heating': ['Q_Joule', 'linear', 'K per timestep'],
                 'Auroral Heating': ['Q_Auroral', 'linear', 'K per timestep'],
                 'Photoelectron Heating': ['Q_PhotoE', 'linear', 'K per timestep'],
                 'Eddy Conduction': ['k_eddy', 'linear', ''],
                 'Adiabatic Eddy Conduction': ['k_adiabaticeddy', 'linear', ''],
                 'NO Cooling': ['L_NO', 'linear', 'K per timestep'],
                 'Molecular Conduction': ['k_molecular', 'linear', ''],
                 'max_electron_density': ['NmF2', 'linear', ''],
                 'max_electron_density_height': ['hmF2', 'linear', 'km'],
                 'Vertical TEC': ['TEC', 'linear', '10**16/m**2'],
                 'HeatFlux_Joule':['phi_qJoule','linear','W/m**2'], 
                 'HeatFlux':['phi_q','linear','W/m**2'], 
                 'HeatFlux_EUV':['phi_qEUV','linear','W/m**2'], 
                 'NO CoolingFlux':['phi_qNOCooling','linear','W/m**2']}

#times from file converted to seconds since midnight of filedate
#plotting input times will be datetime strings of format 'YYYY-MM-DD HH:mm:ss'
#filedate is self.filedate from gitm object
#converts to hours since midnight of filedate for plotting
#def dtsns_to_hrs(datetime_string, filedate):
#    '''Get hours since midnight from datetime string including nanoseconds'''
#    
#    return (datetime.strptime(datetime_string, '%Y-%m-%d %H:%M:%S.%f').replace(tzinfo=timezone.utc)\
#            -filedate).total_seconds()/3600.

def dts_to_hrs(datetime_string, filedate):
    '''Get hours since midnight from datetime string'''
    
    return (datetime.strptime(datetime_string, '%Y-%m-%d %H:%M:%S').replace(tzinfo=timezone.utc)\
            -filedate).total_seconds()/3600.

def filename_to_dts(filename, string_date):
    '''Get datetime string in format "YYYY-MM-SS HH:mm:SS" from filename'''
    
    mmhhss = filename.split('_t')[-1].split('_')[1].split('.bin')[0]
    return string_date+' '+mmhhss[:2]+':'+mmhhss[2:4]+':'+mmhhss[4:]  

def dts_to_ts(file_dts):
    '''Get datetime timestamp in UTC from datetime string'''    
    
    return datetime.timestamp(datetime.strptime(file_dts, '%Y-%m-%d %H:%M:%S'
                                                ).replace(tzinfo=timezone.utc))

#main class object
class GITM(Kamodo): 
    def __init__(self, file_prefix, variables_requested=[], runname="noname",
                 filetimes=False, verbose=False, gridded_int=True, printfiles=True,
                 **kwargs): 
        '''filename must be of form "3D***_tYYMMDD" to load all files for one day,
         and must include a complete path to the files'''
        super(GITM, self).__init__()
        
        #check for .nc files
        total_files = glob.glob(file_prefix+'*')
        self.patterns = np.unique([f.split('/')[-1].split('\\')[-1][:13] \
                                    for f in total_files])
        n_ncfiles = len(glob.glob(file_prefix+'.nc'))
        
        #separate out file directory
        slash_location = file_prefix.rfind('\\')
        if slash_location==-1: slash_location = file_prefix.rfind('/')
        file_dir = file_prefix[:slash_location+1]
        
        #if .nc files not there, create then
        if len(self.patterns)!=n_ncfiles:  #then need to create .nc files for each pattern            
            from kamodo.readers.gitm_tocdf import GITMbin_toCDF as toCDF
            test = [toCDF(file_dir+p) for p in self.patterns]  #get/convert files with given prefix
            if sum(test)!=len(self.patterns): 
                self.conversion_test = False
                return   #if file conversion fails, return
        
        #establish time attributes first for file searching, preferring 3D file
        t0 = ti.perf_counter()
        cdf_data = Dataset(file_dir+self.patterns[-1]+'.nc', 'r')
        string_date = cdf_data.filedate       
        self.filedate = datetime.strptime(string_date+' 00:00:00', \
                                          '%Y-%m-%d %H:%M:%S').replace(tzinfo=timezone.utc) #dt object

        #establish beginning and end times of files found
        files = cdf_data.file.split(',')
        self.datetimes = [filename_to_dts(filename, string_date) for filename \
                          in [files[0], files[-1]]]  #strings in format = YYYY-MM-DD HH:MM:SS 
        self.timerange0={'min':self.datetimes[0], 'max':self.datetimes[-1], 'n':len(files)}
        self.timerange = self.timerange0
        self.filetimes=[dts_to_ts(file_dts) for file_dts in self.datetimes]   #timestamps in UTC     
        if filetimes: 
            return   
        
        #return if only one file found because interpolator code will break
        if len(files)<2:
            print('Not enough files found with given file prefix.')
            return        
        
        #store coordinates of last file for reference, time is the same for all
        self._time = np.array(cdf_data.variables['time'])  #accurate to the sec
        self._lat = np.array(cdf_data.variables['lat'])
        self._lon = np.array(cdf_data.variables['lon'])
        self._height = np.array(cdf_data.variables['height'])
        cdf_data.close()
        
        #store variables
        self.filename = []
        self.runname = runname
        self.missing_value = np.NAN
        self._registered = 0
        self.variables = {}
        self.varfiles = {}  #store which variable came from which file for easier association with coords
        
        #loop through files
        for i in range(len(self.patterns)):
            cdf_data = Dataset(file_dir+self.patterns[i]+'.nc', 'r')
            files = cdf_data.file.split(',')
            self.filename.extend(files)
            
            #store coordinates, time is always the same between files
            setattr(self, '_lat'+str(i), np.array(cdf_data.variables['lat'])) #in deg
            setattr(self, '_lon'+str(i), np.array(cdf_data.variables['lon'])) #in deg
            setattr(self, '_height'+str(i), np.array(cdf_data.variables['height'])) #in km
        
            #check var_list for variables not possible in this file set
            if len(variables_requested)>0:
                gvar_list = [key for key in cdf_data.variables.keys() if key \
                             in variables_requested]
                if len(gvar_list)!=len(variables_requested):
                    err_list = [item for item in variables_requested if item not in \
                                cdf_data.variables.keys()]
                    print(f'Some requested variables are not available in {self.patterns[i]}:', err_list)
                    print('These files have:', cdf_data.variables.keys())
            else:
                gvar_list = [key for key in cdf_data.variables.keys() if key not in \
                             ['time','lat','lon','height']]
            self.varfiles[str(i)] = gvar_list  #store which file these variables came from
            #print('Variables: ', gvar_list)
        
            # Store variable data, units, etc from netCDF file.
            variables = {key:{'units':cdf_data.variables[key].units, 'dtype':np.float32,
                                        'data':np.array(cdf_data.variables[key]), 
                                        'scale': cdf_data.variables[key].datascale} \
                              for key in gvar_list}
            cdf_data.close()
            for key in variables.keys(): self.variables[key] = variables[key] #save to class object
            #this overwrites the TEC data from the 2D files with the calculated TEC data from the 3D files

        if printfiles: 
            print(f'{len(self.filename)} Files:')
            for file in self.filename: print(file)
            
        #register interpolators for each variable
        varname_list = [key for key in self.variables.keys()]  #store original list b/c gridded interpolators
        for varname in varname_list:
            if len(self.variables[varname]['data'].shape)==3:
                #print('3D', varname, self.variables[varname]['data'].shape)
                self.register_3D_variable(self.variables[varname]['units'], 
                                      self.variables[varname]['data'], varname,
                                      gridded_int)
            elif len(self.variables[varname]['data'].shape)==4:
                #print('4D', varname, self.variables[varname]['data'].shape)
                self.register_4D_variable(self.variables[varname]['units'], 
                                      self.variables[varname]['data'], varname,
                                      gridded_int)
        self = RPlot.initialize_4D_plot(self)  #initialize 4D plotting variables 
        if verbose: print(f'{len(varname_list)} variables kamodofied in {ti.perf_counter()-t0:.5f}s.')
        
    #define and register a 3D variable-----------------------------------------
    def register_3D_variable(self, units, variable, varname, gridded_int):
        """Registers a 3d interpolator with 3d signature"""
        
        #determine which file the variable came from
        for i in range(len(self.patterns)-1,-1,-1):  #go in reverse to account for overwritten variables
            if varname in self.varfiles[str(i)]:
                lat = getattr(self, '_lat'+str(i))  #get the correct coordinates
                lon = getattr(self, '_lon'+str(i))
                #print(varname, self.patterns[i])
                break
        
        #define and register the interpolators
        xvec_dependencies = {'time':'hr','lat':'deg','lon':'deg'}
        self = RU.regdef_3D_interpolators(self, units, variable, self._time, 
                                       lat, lon, varname, 
                                       xvec_dependencies, gridded_int)       
        return 
    
    #define and register a 4D variable -----------------------------------------
    def register_4D_variable(self, units, variable, varname, gridded_int):
        """Registers a 4d interpolator with 4d signature"""

        #determine which file the variable came from
        for i in range(len(self.patterns)-1,-1,-1):  #go in reverse to account for overwritten variables
            if varname in self.varfiles[str(i)]:
                lat = getattr(self, '_lat'+str(i))  #get the correct coordinates
                lon = getattr(self, '_lon'+str(i))
                height = getattr(self, '_height'+str(i))
                #print(varname, self.patterns[i])
                break
        
        #define and register the interpolators
        xvec_dependencies = {'time':'hr','height':'km','lat':'deg','lon':'deg'}
        self = RU.regdef_4D_interpolators(self, units, variable, self._time,
                                          height, lat, lon,
                                          varname, xvec_dependencies, gridded_int)
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
        txtbot=f"Model: GITM, dt={self.dt:.4f} hrs, "+\
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