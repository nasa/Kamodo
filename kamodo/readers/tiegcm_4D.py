'''
TIEGCM Kamodo reader, adapted to new structure for satellite flythrough software
Initial version - Asher Pembroke (?)
Adapted by Zachary Waldron and Rebecca Ringuette (June-July 2021)

Remaining tasks:
    - (Zach) check variable dictionary for proper naming, and use of units in Kamodo
    - (Rebecca) deal with 3-4 different versions of ilev, esp. in plot initialization.
'''


import time as ti
import numpy as np
from netCDF4 import Dataset
from datetime import datetime, timezone, timedelta
from kamodo import Kamodo
#import kamodo.readers.reader_plotutilities as RPlot
import kamodo.readers.reader_utilities as RU




### Make a dict of the possible variable names in TIEGCM
###       the intended format is:   "Output VarName":['Latex_Varname', 'Latex_Unit' ]

###  The constituent species are output in units of mass mixing ratio (mmr).
###  Denoted by \psi_i, mmr is the fractional contribution of
###    a species i to the total mass density \rho_{total}, and 
###    is calculated as \psi_i = \rho_i / \rho_{total}, where 
###    \rho_i is the mass density of a constituent species.
tiegcm_varnames={
                 ### 4D Variables, vertical coordinate on midpoint levels (lev)
                 "ZGMID"    : ["H_lev","cm"],    # geometric height- interpolated to the mid points
                 "TN"       : ["T_n","K"],          # neutral temperature    
                 "O2"       : ["psi_O2",""],        # molecular oxygen,   mmr
                 "O1"       : ["psi_O",""],        # atomic oxygen ,   mmr
                 "N2"       : ["psi_N2",""],        # molecular nitrogen,mmr
                 "HE"       : ["psi_He",""],         # helium  ,   mmr
                 "NO"       : ["psi_NO",""],         # nitric oxide , mmr
                 "N4S"      : ["psi_N4S",""],        #  N4S ?,mmr
                 "TE"  : ["T_e","K"],         #  ELECTRON TEMPERATURE,
                 "TI"  : ["T_i","K"],         #  ION TEMPERATURE
                 "O2P" : ["N_O2plus","1/cm**3"],  #  O2+ ION
                 "OP"  : ["N_Oplus","1/cm**3"],    #   O+ ION
                 "N2N"      : ["psi_N2",""],    # molecular nitrogen (maybe number density),mmr
                 "CO2_COOL" : ["Q_CO2cool","erg/g/s"],  #  CO2 cooling rates
                 "NO_COOL"  : ["Q_NOcool","erg/g/s"],      #  NO cooling rates
                 "UN"  : ["u_n","cm/s"],            #  neutral ZONAL wind (+EAST)
                 "VN"  : ["v_n","cm/s"],            #  neutral MERIDIONAL wind (+NORTH)
                # 
                ### 4D Variables, vertical coordinate on interface levels (ilev)
                 "DEN"      :["rho","g/cm**3"],     # total neutral mass density  
                 "ZG"       :["H_ilev","cm"],          # geometric height  
                 "Z"        :["H_geopot","cm"],            # geopotential height (cm)  
                 "NE"       : ["N_e","1/cm**3"],    #  ELECTRON DENSITY
                 "OMEGA"    : ["omega","1/s"],      #  VERTICAL MOTION
                 "POTEN"    : ["V","V"],        #  ELECTRIC POTENTIAL
                ### 4D Variables, vertical coordinate on interface mag levels (imlev)
                 "ZMAG"  : ["H_imlev","cm"],     #  Geopotential Height on Geomagnetic Grid
                #
                ### 3D Variables,    (time, lat, lon)
                 "TEC"  : ["TEC","1/cm**2"],     #  Total Electron Content
                 "TLBC"  : ["T_nLBC","K"],       #  Lower boundary condition for TN
                 "ULBC"  : ["u_nLBC","cm/s"],    #  Lower boundary condition for UN
                 "VLBC"  : ["v_nLBC","cm/s"],    #  Lower boundary condition for VN
                 "TLBC_NM"  : ["T_nLBCNM",""],  #  Lower boundary condition for TN (TIME N-1)
                 "ULBC_NM"  : ["u_nLBCNM",""],  #  Lower boundary condition for UN (TIME N-1)
                 "VLBC_NM"  : ["v_nLBCNM",""],  #  Lower boundary condition for VN (TIME N-1)
                  }

#####--------------------------------------------------------------------------------------
##### Define some helpful functions for dealing with time systems

def dts_to_ts(file_dts):
    '''Get datetime timestamp in UTC from datetime string'''    
    
    return datetime.timestamp(datetime.strptime(file_dts, '%Y-%m-%d %H:%M:%S'
                                                ).replace(tzinfo=timezone.utc))

def min_to_dt(time_minutes):  #called by functions below
    '''Convert minutes since 2000-03-20 00:00:00 to datetime object in UTC'''
    
    return datetime(2000,3,20).replace(tzinfo=timezone.utc)+timedelta(minutes=time_minutes)
        
def min_to_dts(time_minutes):
    '''Convert minutes since 2000-03-20 00:00:00 to a datetime string'''
    
    return datetime.strftime(min_to_dt(time_minutes), '%Y-%m-%d %H:%M:%S')

def min_to_date(time_minutes):
    '''Use minutes since 2000-03-20 00:00:00 to determine the date in the file. Returns a string and a datetime object.'''
    
    date_string = datetime.strftime(min_to_dt(time_minutes), '%Y-%m-%d')  #'YYYY-MM-DD'
    return datetime.strptime(date_string, '%Y-%m-%d').replace(tzinfo=timezone.utc)

@np.vectorize
def min_to_hrs(time_minutes, date_dt):
    '''Convert  minutes since 2000-03-20 00:00:00 to hours since midnight using predetermined datetime object.'''
    
    return (min_to_dt(time_minutes)-date_dt).total_seconds()/3600.

##### --------------------------------------------------------------------------------------
###   Construct a TIEGCM class that inherits Kamodo       
class TIEGCM(Kamodo): 
    def __init__(self, filename, variables_requested=[], runname="noname",
                 filetimes=False, verbose=False, gridded_int=True, printfiles=True,
                 **kwargs):  #filename should include the full path
        
        #### Use a super init so that your class inherits any methods from Kamodo
        super(TIEGCM, self).__init__()

        #store time information for satellite flythrough layer to choose the right file
        t0 = ti.perf_counter()
        cdf_data = Dataset(filename, 'r')
        time = np.array(cdf_data.variables['time'])  #in minutes since 2000-03-20 00:00:00 UTC
        self.filedate = min_to_date(time[0])  #datetime object for file date at midnight UTC
        self.datetimes = [min_to_dts(time_minutes) for time_minutes \
                          in [time[0], time[-1]]]  #strings in format = YYYY-MM-DD HH:MM:SS
        self.timerange0={'min':self.datetimes[0], 'max':self.datetimes[-1], 'n':len(time)}
        self.timerange = self.timerange0
        self.filetimes=[dts_to_ts(file_dts) for file_dts in self.datetimes]   #timestamps in UTC     
        if filetimes:
            return      
                
        #### Store our inputs as class attributes to the class
        self.filename      = filename
        self.runname       = runname
        self.missing_value = np.NAN
        self._registered   = 0
        self.variables     = dict()
        if printfiles: 
            print('Files:', self.filename)
        
        #translate from standardized variables to names in file
        #remove variables requested that are not in the file
        if len(variables_requested)>0:
            gvar_list = [key for key, value in tiegcm_varnames.items() \
                             if value[0] in variables_requested]  # file variable names
            if len(gvar_list)!=len(variables_requested):
                err_list = [tiegcm_varnames[key][0] for key in gvar_list \
                            if key not in list(cdf_data.variables.keys())]
                print('Some requested variables are not available:', err_list)
        else:
            gvar_list = [key for key in cdf_data.variables.keys() \
                         if key in tiegcm_varnames.keys() ] 
        if 'ZG' not in gvar_list: gvar_list.append('ZG')  #force addition of H for conversion of ilev to H and back
        
        #### Store coordinate data as class attributes    
        self._time = min_to_hrs(time, self.filedate)  #convert to hours since midnight
        self._ilev = np.array(cdf_data.variables['ilev'])
        lat = np.array(cdf_data.variables['lat'])  #NOT FULL RANGE IN LATITIUDE!!!
        lat = np.insert(lat, 0, lat[0]-np.diff(lat).min())  #insert a grid point at beginning (before -87.5)
        self._lat = np.append(lat, lat[-1]+np.diff(lat).min())   #and at the end (after 87.5)
        lon = np.array(cdf_data.variables['lon'])+180.  #NOT WRAPPED IN LONGITUDE!!!!!
        self._lon = np.append(lon, 360.)  #add 360. to end of array
        self._lev = np.array(cdf_data.variables['lev'])
        self._imlev = np.array(cdf_data.variables['imlev'])  #'mlev' isn't used by any of the variables
        self._mlat = np.array(cdf_data.variables['mlat'])
        self._mlon = np.array(cdf_data.variables['mlon'])+180.  #shifting to 0-360 range
        self._height = np.array([0.])  #for compatibility only
        
        #### Store the requested variables into a dictionary 
        ####     as the variables attributes
        #### This will contain units, dtype, and the data
        #print(gvar_list)
        self.variables = {tiegcm_varnames[key][0]:{'units':tiegcm_varnames[key][-1], 
                                                   'dtype':np.float32,
                               'data':np.array(cdf_data.variables[key])}\
                          for key in gvar_list}  #store with key = standardized name
        cdf_data.close()
        if verbose: print(f'Took {ti.perf_counter()-t0:.6f}s to read in data')
            
        #### register interpolators for each requested variable
        varname_list = [key for key in self.variables.keys()]  #store original list b/c gridded interpolators
        t_reg = ti.perf_counter()
        for varname in varname_list:
            if len(self.variables[varname]['data'].shape)==3:
                self.register_3D_variable(self.variables[varname]['units'], 
                                      self.variables[varname]['data'], varname,
                                      gridded_int)
            elif len(self.variables[varname]['data'].shape)==4:
                self.register_4D_variable(self.variables[varname]['units'], 
                                      self.variables[varname]['data'], varname,
                                      gridded_int)
        if verbose: print(f'Took {ti.perf_counter()-t_reg:.5f}s to register '+\
                          f'{len(varname_list)} variables.')
        self = RPlot.initialize_4D_plot(self)  #initialize plots  
        if verbose: print(f'Took a total of {ti.perf_counter()-t0:.5f}s to kamodofy '+\
                          f'{len(varname_list)} variables.')

    def wrap_3Dlatlon(self, varname, variable):
        '''Wraps the data array in longitude (0=360), and latitude (0=-2, -1=1)'''
    
        shape_list = list(variable.shape)  #e.g. time, lat, lon
        shape_list[1]+=2  #need two more places in latitude
        shape_list[2]+=1  #need one more place in longitude
        tmp_arr = np.zeros(shape_list)  #array to set-up wrapped data in
        tmp_arr[0:,1:-1,:-1]=variable  #copy data into grid
        tmp_arr[:,1:-1,-1] = variable[:,:,0]  #wrap in longitude first
        tmp_arr[:,0,:] = np.flip(tmp_arr[:,1,:],axis=1)  #wrap in latitude...
        tmp_arr[:,-1,:] = np.flip(tmp_arr[:,-2,:],axis=1)  #reverse past poles
        self.variables[varname]['data'] = tmp_arr
        return tmp_arr        

    def wrap_4Dlatlon(self, varname, variable):
        '''Wraps the data array in longitude (0=360), and latitude (0=-2, -1=1)'''
    
        shape_list = list(variable.shape)  #e.g. time, ilev, lat, lon
        shape_list[2]+=2  #need two more places in latitude
        shape_list[3]+=1  #need one more place in longitude
        tmp_arr = np.zeros(shape_list)  #array to set-up wrapped data in
        tmp_arr[:,:,1:-1,:-1]=variable  #copy data into grid
        tmp_arr[:,:,1:-1,-1] = variable[:,:,:,0]  #wrap in longitude first
        tmp_arr[:,:,0,:] = np.flip(tmp_arr[:,:,1,:],axis=2)  #wrap in latitude...
        tmp_arr[:,:,-1,:] = np.flip(tmp_arr[:,:,-2,:],axis=2)  #reverse past poles
        self.variables[varname]['data'] = tmp_arr
        return tmp_arr
                
    ##### Define and register a 3D variable -----------------------------------------
    def register_3D_variable(self, units, variable, varname, gridded_int):
        """Registers a 3d interpolator with 3d signature"""
        
        #define and register the interpolators
        xvec_dependencies = {'time':'hr','lat':'deg','lon':'deg'}
        wrapped_data = self.wrap_3Dlatlon(varname, variable)
        self = RU.regdef_3D_interpolators(self, units, wrapped_data, self._time, 
                                       self._lat, self._lon, varname, 
                                       xvec_dependencies, gridded_int)       
        return 
    
    
    #### Define and register a 4D variable -----------------------------------------
    def register_4D_variable(self, units, variable, varname, gridded_int):
        """Registers a 4d interpolator with 4d signature"""
        
        #These lists need to be the standardized variable name to match that above,
        #not the names from the data file.
        ilev_list = ["rho","H_ilev","H_geopot","N_e","omega","V"]                                                     # index with ilev
        lev_list  = ['T_n','u_n','v_n','psi_O2','psi_O','psi_N2','psi_He','H_lev',
                     'psi_NO','psi_N4S','T_e','T_i','N_O2plus','N_Oplus','Q_CO2cool',
                     'Q_NOcool'] # index with lev
        imlev_list = ['H_imlev']
             
        ####  Get the correct coordinates
        if varname in ilev_list:
            vert_coord = self._ilev
            lat, lon = self._lat, self._lon
            xvec_dependencies = {'time':'hr', 'ilev':'m/m','lat':'deg','lon':'deg'}
        elif varname in lev_list:
            vert_coord = self._lev
            lat, lon = self._lat, self._lon
            xvec_dependencies = {'time':'hr', 'lev':'m/m','lat':'deg','lon':'deg'}
        elif varname in imlev_list:
            vert_coord = self._imlev
            lat, lon = self._mlat, self._mlon
            xvec_dependencies = {'time':'hr', 'imlev':'m/m','mlat':'deg','mlon':'deg'}

        #### define and register the interpolators
        if 'lat' in xvec_dependencies.keys():
            wrapped_data = self.wrap_4Dlatlon(varname, variable)
        else: 
            wrapped_data = variable  #no wrapping needed for mlat and mlon
        self = RU.regdef_4D_interpolators(self, units, wrapped_data, 
                                          self._time, vert_coord, lat, lon,
                                          varname, xvec_dependencies, gridded_int)
        return

"""#begin plotting code -----------------------------------
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

        #need to determine mlat/mlon vs lat/lon and ilev/lev/mlev dependence here with dependencies
        arg_dict = self.variables[var]['xvec']
        if 'mlat' in arg_dict.keys(): 
            self.dLON, self.nLON = self.dmlon, len(self._mlon)
            self.dLAT, self.nLAT = self.dmlat, len(self._mlat)
            self.dH, self.nH = self.dimlev, len(self._milev)
        else: 
            self.dLON, self.nLON = self.dlon, len(self._lon)
            self.dLAT, self.nLAT = self.dlat, len(self._lat)
        #will eventually need logic here and elsewhere to deal with multiple ilev types******************
        self.gridSize = self.nLAT*self.nLON*self.numZ

        print(f'set_plot::colorscale={colorscale}, datascale={datascale}')
        print(f'Run: {self.runname}')
        #Set some text strings
        txtbot=f"Model: TIEGCM, dt={self.dt:.4f} hrs, "+\
            f"dlat={self.dlat:.1f} deg, dlon={self.dlon:.1f} deg, dz={self.dz:.1f} {self.dzunit}."
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
"""
