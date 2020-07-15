import numpy as np

from kamodo import Kamodo, kamodofy, gridify

from scipy.interpolate import RegularGridInterpolator, interp1d
from netCDF4 import Dataset
import time
from datetime import datetime,timedelta,timezone

import numpy.ma as ma

# constants and dictionaries
ctipe_kamodo_variable_names = dict( density='rho',
                                    temperature='T',
                                    electron_temperature='T_e',
                                    ion_temperature='T_i',
                                    height='H',
                                    meridional_neutral_wind='Vn_lat',
                                    zonal_neutral_wind='Vn_lon',
                                    vertical_neutral_wind='Vn_H',
                                    neutral_temperature='T_n',
                                    mean_molecular_mass='Rmt',
                                    electron_density='N_e',
                                    neutral_density='N_n',
                                    solar_heating='Q_Solar',
                                    joule_heating='Q_Joule',
                                    radiation_heat_cool='Q_radiation',
                                    atomic_oxygen_density='N_O',
                                    molecular_oxygen_density='N_O2',
                                    molecular_nitrogen_density='N_N2',
                                    nitric_oxide_density='N_NO',
                                    nitric_oxide_ion_density='N_NOplus',
                                    molecular_nitrogen_ion_density='N_N2plus', # the "+" might not work
                                    molecular_oxygen_ion_density='N_Oplus',
                                    atomic_nitrogen_ion_density='N_Nplus',
                                    atomic_oxygen_ion_density='N_Oplus',
                                    atomic_hydrogen_ion_density='N_Hplus',
                                    pedersen_conductivity ='sigma_P',
                                    hall_conductivity='sigma_H',
                                    meridional_ion_velocity='Vi_lat',
                                    zonal_ion_velocity='Vi_lon',
                                    height_integrated_joule_heating='W_Joule',
                                    energy_influx='Eflux_precip',
                                    mean_energy='Eavg_precip',
                                    total_electron_content='TEC',
                                    theta_electric_field_at_140km='E_theta140km',
                                    lambda_electric_field_at_140km='E_lambda140km',
                                    theta_electric_field_at_300km='E_theta300km',
                                    lambda_electric_field_at_300km='E_lambda300km')

# reverse lookup in above dictionary
def ctipe_varnames(kamodo_varnames):
    keys=[]
    for varname in kamodo_varnames:
        for key,value in ctipe_kamodo_variable_names.items():
            if varname == value:
                keys.append(key)
    return(keys)


# some functions from tiegcm.util
def to_range(x, start, limit):
        """wraps x into range [start, limit]"""
        return start + (x - start) % (limit - start)

def parse_units(varname, variable):
    """Parses units based on variable input"""
    try:
        units = variable.units
        units = units.replace('-','**-')
        units = units.replace('mmr','1').replace('S', 'ohm**-1')
        units = units.replace('cm3','cm**3').replace('erg/g/s', 'erg/(g*s)')
        units = units.replace('none','1')
    except:
        units = ''
        if varname == 'Z':
            units = 'cm'
    return units

def sort_variables(self, variables):
    variables_3d = []
    variables_4d = []
    for varname in variables:
        try:
            varshape = self._ctipe.variables[varname].shape
        except:
            continue
        if len(varshape) == 4:
            variables_4d.append(varname)
        elif len(varshape) == 3:
            variables_3d.append(varname)
    return variables_4d + variables_3d

def ctipe_kamodo_variable_name(varname):

    if varname in ctipe_kamodo_variable_names:
        return(ctipe_kamodo_variable_names[varname])
    else:
        return(varname)

def totalseconds_to_datetime(seconds):
    datetime0=datetime(1970,1,1,0,0,0,tzinfo=timezone.utc)
    delta=timedelta(seconds=seconds)
    datetime1=datetime0+delta
    return(datetime1)

def seconds_since_1970(time):
    delta_t=time-datetime(1970,1,1)
    return(delta_t.total_seconds())

class CTIPe(Kamodo):
    def __init__(self, filename,
                 date = None,
                 time = None,
                 variables_requested = None,
                 debug= False,
                 runpath = "./", 
                 runname = "noname",
                 **kwargs):
        import re
# input file name can be one of the 4 files for each day of model outputs
#  YYYYMMDD-plot-[density|height|neutral|plasma].nc files
        if filename.find('plot-density') > 0:
            filename_density=filename
            filename_height=filename.replace('plot-density','plot-height')
            filename_neutral=filename.replace('plot-density','plot-neutral')
            filename_plasma=filename.replace('plot-density','plot-plasma')

        if filename.find('plot-height') > 0:
            filename_density=filename.replace('plot-height','plot-density')
            filename_height=filename
            filename_neutral=filename.replace('plot-height','plot-neutral')
            filename_plasma=filename.replace('plot-height','plot-plasma')

        if filename.find('plot-neutral') > 0:
            filename_density=filename.replace('plot-neutral','plot-density')
            filename_height=filename.replace('plot-neutral','plot-height')
            filename_neutral=filename
            filename_plasma=filename.replace('plot-neutral','plot-plasma')

        if filename.find('plot-plasma') > 0:
            filename_density=filename.replace('plot-plasma','plot-density')
            filename_height=filename.replace('plot-plasma','plot-height')
            filename_neutral=filename.replace('plot-plasma','plot-neutral')
            filename_plasma=filename

# only the density, height and neutral files have data and are read
        print(filename_density)
        print(filename_height)
        print(filename_neutral)
        print(filename_plasma)
        if filename_density == filename_height or filename_density == filename_neutral or filename_height == filename_neutral:
            raise AttributeError("Error: filenames are ill-defined")
        self.debug=debug

        self._ctipe_density = Dataset(filename_density)
        self._ctipe_height = Dataset(filename_height)
        self._ctipe_neutral = Dataset(filename_neutral)


# keep track whether periodic boundary in longiutde has been applied
#        self.wrapped=[]
        self.wrapped_height=[]
        self.wrapped_density=[]
        self.wrapped_neutral=[]

        
        self._time = np.array(self._ctipe_density.variables['time'])
        
        # implement time picker based on nearest neighbor
        if date is None or time is None:
            time_index=0
            norm_distance=0
            self.datetime=totalseconds_to_datetime(self._time[time_index])
            self.filetime=self.datetime.strftime("%Y/%m/%d %H:%M:%S")
            self.time_index=time_index
            self.time_index0=time_index
            self.time_index1=time_index+1
            self.norm_distance=norm_distance 
        else:
            date_time=datetime.strptime(date+' '+time,'%Y/%m/%d %H:%M:%S')
            seconds=seconds_since_1970(date_time)
            print(self._time)
            print(seconds)
            time_index=np.searchsorted(self._time,seconds)-1
            print(time_index)
            if time_index < 0:
                time_index = 0
            if time_index > self._time.size - 2:
                time_index = self._time.size - 2
            norm_distance=(seconds - self._time[time_index]) / (self._time[time_index + 1] - self._time[time_index])
# in case we implement time interpolation
            self.time_index0=time_index
            self.time_index1=time_index+1
            self.norm_distance=norm_distance 
# nearest-neighbor time pick
            print('Norm_distance: ',norm_distance)
            if norm_distance >=0.5:
                time_index=time_index+1
            self.time_index=time_index
            self.datetime=totalseconds_to_datetime(self._time[time_index])
            self.filetime=self.datetime.strftime("%Y/%m/%d %H:%M:%S")
            
        
        self._ilev_density = np.array(self._ctipe_density.variables['plev'])
        self._lat_density = np.array(self._ctipe_density.variables['lat'])
        self._lon_density = np.array(self._ctipe_density['lon'])
        self.lonrange=dict(min=0,max=360,n=73)
        self.latrange=dict(min=-90,max=90,n=37)
        self.iprange=dict(min=self._ilev_density.min(),max=self._ilev_density.max(),n=len(self._ilev_density))
        self.hrange=self.iprange
        self.heightrange=dict(min=90000.,max=500000.,n=len(self._ilev_density))
        
        self.gridSize=len(self._ilev_density)*len(self._lat_density)*len(self._lon_density)
        self.gridMinDx=np.asarray([self._lon_density[1]-self._lon_density[0],self._lat_density[1]-self._lat_density[0]]).min()
        
        print(self._ctipe_height.variables.keys())
        if "ht" in self._ctipe_height.variables.keys():
            self._ht_height = np.array(self._ctipe_height.variables['ht'])
        if "height" in self._ctipe_height.variables.keys():
            self._ht_height = np.array(self._ctipe_height.variables['height'])
        self._lat_height = np.array(self._ctipe_height.variables['lat'])
        self._lon_height = np.array(self._ctipe_height['lon'])

#        print(self._ctipe_neutral.variables)
        self._ilev_neutral = np.array(self._ctipe_neutral.variables['plev'])
        self._lat_neutral = np.array(self._ctipe_neutral.variables['lat'])
        self._lon_neutral = np.array(self._ctipe_neutral['lon'])
        self._elat_neutral = np.array(self._ctipe_neutral.variables['elat'])
        self._elon_neutral = np.array(self._ctipe_neutral['elon'])

#        return
    
# longitude-wrap coordinate arrays
        if self._lon_density.max() < 360:
            self._lon_density=np.append(self._lon_density,[360])
        if self._lon_height.max() < 360:
            self._lon_height=np.append(self._lon_height,[360])
        if self._lon_neutral.max() < 360:
            self._lon_neutral=np.append(self._lon_neutral,[360])
        if self._elon_neutral.max() < 360:
            self._elon_neutral=np.append(self._elon_neutral,[360])

        self._registered = 0
        
        super(CTIPe, self).__init__() 
        print('opening {}'.format(filename))
        self.filename = filename
        self.runpath = runpath
        self.runname = runname
#        self.start_time = start_time
        self.missing_value = np.NAN
        self.variables=dict()
        
# Interpolation values
        self.plots = dict()
        self.plottype = "" # empty default to force reinterpolation and regeneration
        self.cut = 'IP'
        self.cutV = 10.
        self.nX = 1
        self.nY = 37
        self.nZ = 73
        self.newx = self.cutV
        self.newy = np.linspace(-90., 90., self.nY)
        self.newz = np.linspace(0., 360., self.nX)
        self.xunit=''
        self.yunit='deg'
        self.zunit='deg'
        self.tol = 1.1
        self.plots[self.plottype] = dict(cut=self.cut, cutV=self.cutV, tol=self.tol,
                                         nX=self.nX, nY=self.nY, nZ=self.nZ,
                                         newx=self.newx, newy=self.newy, newz=self.newz)  
        
        # Get grids ready to use
        self.setup_interpolating_grids()
        
        if variables_requested is None:
            density_variables = self._ctipe_density.variables.keys()
            height_variables  = self._ctipe_height.variables.keys()
            neutral_variables = self._ctipe_neutral.variables.keys()
        else:
            density_variables = []
            height_variables  = []
            neutral_variables = []

 # add height variable needed to height (not IP-level) interpolatioms
            if 'H' not in variables_requested:
                variables_requested.append('H')
                    
            for variable in ctipe_varnames(variables_requested):
                variable_found=False
                if variable in self._ctipe_density.variables.keys():
                    density_variables.append(variable)
                    variable_found=True
                if variable in self._ctipe_height.variables.keys():
                    height_variables.append(variable)
                    variable_found=True
                if variable in self._ctipe_neutral.variables.keys():
                    neutral_variables.append(variable)
                    variable_found=True
                if variable_found is False:
                    raise AttributeError("No such variable:{}".format(variable))
                
        print(self._ctipe_density.variables.keys())
        print(self._ctipe_height.variables.keys())
        print(self._ctipe_neutral.variables.keys())
        print([density_variables,height_variables,neutral_variables])           #        if len(density_variables) == 0 and len(height_variables) == 0 and len(neutral_variables) == 0:
#            print("Error: requested variable not recognized")
#            raise
            
#        variables = sort_variables(self, variables)

        for varname in density_variables:
            print(varname)
            units = parse_units(varname, self._ctipe_density.variables[varname])
# add layer at lon=360
            self.wrap_density_variable(varname)
            variable = self._ctipe_density.variables[varname]
            if len(variable.shape) not in [3,4]:
                # skip this variable
                continue
                
            if varname == 'density':
                units='kg/m**3'
                
            if self.debug:
                print(varname)
                print(variable)
                print(units)

            if varname == 'ZMAG':
                continue

            elif len(variable.shape) == 4:
                kamodo_varname=ctipe_kamodo_variable_names[varname]
                data=np.squeeze(variable[time_index,:,:,:])
                self.variables[kamodo_varname] = dict(units = units, data = data)
                self.register_4d_density_variable(units, data, varname)

            elif len(variable.shape) == 3:
                kamodo_varname=ctipe_kamodo_variable_names[varname]
                data=np.squeeze(variable[time_index,:,:])
                self.variables[kamodo_varname] = dict(units = units, data = data)
                self.register_3d_density_variable(units, data, varname)


        for varname in height_variables:
            print(varname)
#            kamodo_varname=ctipe_kamodo_variable_names[varname]
            units = parse_units(varname, self._ctipe_height.variables[varname])
# add layer at lon=360
            self.wrap_height_variable(varname)
            variable = self._ctipe_height.variables[varname]
            if len(variable.shape) not in [3,4]:
                # skip this variable
                continue
                
            if varname == 'density':
                units='kg/m**3'
                
            if self.debug:
                print(varname)
                print(variable)
                print(units)

            if varname == 'ZMAG':
                continue

            elif len(variable.shape) == 4:
                kamodo_varname=ctipe_kamodo_variable_names[varname]
                data=np.squeeze(variable[time_index,:,:,:])
                self.variables[kamodo_varname] = dict(units = units, data = data)
                self.register_4d_height_variable(units, np.squeeze(variable[time_index,:,:,:]), varname)

            elif len(variable.shape) == 3:
                kamodo_varname=ctipe_kamodo_variable_names[varname]
                data=np.squeeze(variable[time_index,:,:])
                self.variables[kamodo_varname] = dict(units = units, data = data)
                self.register_3d_height_variable(units, np.squeeze(variable[time_index,:,:]), varname)

                
        for varname in neutral_variables:
            print(varname)
            units = parse_units(varname, self._ctipe_neutral.variables[varname])
# add layer at lon=360
            self.wrap_neutral_variable(varname)
            variable = self._ctipe_neutral.variables[varname]
            if len(variable.shape) not in [3,4]:
                # skip this variable
                continue
                
            if varname == 'density':
                units='kg/m**3'
                
            if self.debug:
                print(varname)
                print(variable.shape)
                print(units)
            
            if varname == 'ZMAG':
                continue

            elif len(variable.shape) == 4:
                kamodo_varname=ctipe_kamodo_variable_names[varname]            
                data=np.squeeze(variable[time_index,:,:,:])                    
                self.variables[kamodo_varname] = dict(units = units, data = data)
                self.register_4d_neutral_variable(units, np.squeeze(variable[time_index,:,:,:]), varname)

            elif len(variable.shape) == 3:
                kamodo_varname=ctipe_kamodo_variable_names[varname]            
                data=np.squeeze(variable[time_index,:,:])
                self.variables[kamodo_varname] = dict(units = units, data = data)
                self.register_3d_neutral_variable(units, np.squeeze(variable[time_index,:,:]), varname)

                
        print('registered {} variables'.format(self._registered))
        
        # register user's input variables, assuming kamodo-compatible
#        for varname, variable in kwargs.items():
#            self[varname] = variable

#        def ilev(points):
#            return self.lev(*points)

#        self['ilev'] = ilev

    def setup_interpolating_grids(self):
        '''setup the grid to interpolate to, trim to necessary size of source grid, and compute interpolation weights'''
        # Setup grid to interpolate to
        xx, yy, zz = np.meshgrid(self.newx,self.newy,self.newz, indexing = 'xy')
        self.newgrid = np.ndarray(shape=(np.size(np.reshape(xx,-1)),3), dtype=np.float32)
        self.newgrid[:,0] = np.reshape(xx,-1)
        self.newgrid[:,1] = np.reshape(yy,-1)
        self.newgrid[:,2] = np.reshape(zz,-1)
        self.plots[self.plottype]['newgrid'] = self.newgrid
        
        # Reduce size of read data block to speed up the interpolation
#        self.grid2 = self.grid[(self.grid[:,0] > (np.amin(self.newx)-self.tol)) & 
#                               (self.grid[:,0] < (np.amax(self.newx)+self.tol)) &
#                               (self.grid[:,1] > (np.amin(self.newy)-self.tol)) & 
#                               (self.grid[:,1] < (np.amax(self.newy)+self.tol)) & 
#                               (self.grid[:,2] > (np.amin(self.newz)-self.tol)) & 
#                               (self.grid[:,2] < (np.amax(self.newz)+self.tol))]


        # Precompute weights for interpolation to optimize
#        self.vtx, self.wts = self.GMcompute_weights(self.grid2, self.newgrid)

        return


    def register_3d_density_variable(self, units, variable, varname):
        """Registers a 3d interpolator with 3d signature"""

        rgi = RegularGridInterpolator((self._lat_density, self._lon_density),
                                      variable, bounds_error = False)
#        @kamodofy(units = units, data = variable,citation="Rastaetter 2020")
#        @gridify(lat = self._lat_density, lon = self._lon_density)
        def interpolator(xvec):
            """Interpolates 3d density variable"""
            return rgi(xvec)

        self.variables[varname]['interpolator'] = interpolator
#        self[varname] = interpolator
        self[varname]= kamodofy(interpolator, 
                                 units = units, 
                                 citation = "Rastaetter 2020",
                                 data = variable)
        self._registered += 1

    def register_4d_density_variable(self, units, variable, varname):
        """Registers a 4d interpolator with 4d signature"""
#        print(varname,self._time.shape,self._ilev_density.shape,
#              self._lat_density.shape, self._lon_density.shape, 
#              variable.shape)
            
#        print(self._lon_density)


        rgi = RegularGridInterpolator((self._ilev_density, self._lat_density, self._lon_density), 
                                      variable, bounds_error = False)

#        @kamodofy(units = units, data = variable)
#        @gridify(ilev = self._ilev_density, lat = self._lat_density, lon = self._lon_density) 
        def interpolator(xvec):
            """Interpolates 4d density variable"""
            return rgi(xvec)

        try:
            self[ctipe_kamodo_variable_name(varname)] = interpolator
            self.variables[ctipe_kamodo_variable_name(varname)]['interpolator']=interpolator
        except:
            self[varname]= kamodofy(interpolator, 
                                 units = units, 
                                 citation = "Rastaetter 2020",
                                 data = variable)
        #    self[varname] = interpolator
            self.variables[varname]['interpolator']=dict(interpolator=interpolator,data=data,units=units)

        self._registered += 1

    def register_3d_height_variable(self, units, variable, varname):
        """Registers a 3d interpolator with 3d signature"""

        rgi = RegularGridInterpolator((self._lat_height, self._lon_height),
                                      variable, bounds_error = False)
#        @kamodofy(units = units, data = variable)
#        @gridify(lat = self._lat_height, lon = self._lon_height)
        def interpolator(xvec):
            """Interpolates 3d height variable"""
            return rgi(xvec)

        try:
            self[ctipe_kamodo_variable_name(varname)] = kamodofy(interpolator, 
                                 units = units, 
                                 citation = "Rastaetter 2020",
                                 data = variable)
            self.variables[ctipe_kamodo_variable_name(varname)]['interpolator']=interpolator
#            self[ctipe_kamodo_variable_name(varname)] = interpolator
        except:
#            self[varname] = interpolator            
            self[varname]= kamodofy(interpolator, 
                                 units = units, 
                                 citation = "Rastaetter 2020",
                                 data = variable)
            self.variables[varname]=dict(interpolator=interpolator,data=variable,units=units)
        self._registered += 1

    def register_4d_height_variable(self, units, variable, varname):
        """Registers a 4d interpolator with 4d signature"""

        rgi = RegularGridInterpolator((self._ht_height, self._lat_height, self._lon_height), 
                                      variable, bounds_error = False)
#        @kamodofy(units = units, data = variable)
#        @gridify(ilev = self._ht_height, lat = self._lat_height, lon = self._lon_height) 
        def interpolator(xvec):
            """Interpolates 4d  height variable"""
            return rgi(xvec)

        try:
#            self[ctipe_kamodo_variable_name(varname)] = interpolator
            self[ctipe_kamodo_variable_name(varname)] = kamodofy(interpolator, 
                                 units = units, 
                                 citation = "Rastaetter 2020",
                                 data = variable)
            self.variables[ctipe_kamodo_variable_name(varname)]['interpolator']=interpolator
        except:
#            self[varname] = interpolator            
            self[varname]= kamodofy(interpolator, 
                                 units = units, 
                                 citation = "Rastaetter 2020",
                                 data = variable)
            self.variables[varname]=dict(interpolator=interpolator,data=variable,units=units)

        self._registered += 1

    def register_3d_neutral_variable(self, units, variable, varname):
        """Registers a 3d interpolator with 3d signature"""
#  check index 0, not index 1 since time has been stripped out
        if (variable.shape)[0] == self._elat_neutral.shape[0]:
            rgi = RegularGridInterpolator((self._elat_neutral, self._elon_neutral),
                                      variable, bounds_error = False)
#            @kamodofy(units = units, data = variable )
# to use the density grid one may have to run interpolation to provide data for the above kamodofy decorator
#           variable_on_lat_lon_neutral=rgi(self._time,???,???,????)
#            @kamodofy(units = units, data = variable_on_lat_lon_neutral )
#            @gridify(t = self._time, lat = self._lat_neutral, lon = self._lon_neutral)    
#            @gridify(lat = self._elat_neutral, lon = self._elon_neutral)    
            def interpolator(xvec):
                """Interpolates 3d neutral variable"""
                return rgi(xvec)
        else:
            rgi = RegularGridInterpolator((self._lat_neutral, self._lon_neutral),
                                      variable, bounds_error = False)
#            @kamodofy(units = units, data = variable)
#            @gridify(lat = self._lat_neutral, lon = self._lon_neutral)    
            def interpolator(xvec):
                """Interpolates 3d neutral variable"""
                return rgi(xvec)

        try:
#            self[ctipe_kamodo_variable_name(varname)] = interpolator
            self[ctipe_kamodo_variable_name(varname)] = kamodofy(interpolator, 
                                 units = units, 
                                 citation = "Rastaetter 2020",
                                 data = variable)
            self.variables[ctipe_kamodo_variable_name(varname)]['interpolator']=interpolator
        except:
#            self[varname] = interpolator
            self[varname] = kamodofy(interpolator, 
                                 units = units, 
                                 citation = "Rastaetter 2020",
                                 data = variable)
            self.variables[varname]['interpolator']=interpolator
        self._registered += 1

    def register_4d_neutral_variable(self, units, variable, varname):
        """Registers a 4d interpolator with 4d signature"""
#        print(varname,self._time.shape,self._ilev_density.shape,
#              self._lat_density.shape, self._lon_density.shape, 
#              variable.shape)
        if variable.shape[-1] < self._lon_density.shape[0]:
            self.wrap_neutral_variable(varname)

        rgi = RegularGridInterpolator((self._ilev_neutral, self._lat_neutral, self._lon_neutral), 
                                      variable, bounds_error = False)
#        @kamodofy(units = units, data = variable)
#        @gridify(ilev = self._ilev_neutral, lat = self._lat_neutral, lon = self._lon_neutral) 
        def interpolator(xvec):
            """Interpolates 4d neutral variable"""
            return rgi(xvec)

        try:
#            self[ctipe_kamodo_variable_name(varname)] = interpolator
            self[ctipe_kamodo_variable_name(varname)] = kamodofy(interpolator, 
                                 units = units, 
                                 citation = "Rastaetter 2020",
                                 data = variable)
            self.variables[ctipe_kamodo_variable_name(varname)]['interpolator']=interpolator
        except:
# this alternative is invoked when the Kamodo variable name causes problems
#            self[varname] = interpolator
            self[varname] = kamodofy(interpolator, 
                                 units = units, 
                                 citation = "Rastaetter 2020",
                                 data = variable)
            self.variables[varname]['interpolator']=interpolator
            self._registered += 1


    @np.vectorize
    def lev_density(self, z, lat, lon):
        """Finds ilev for a given height"""
        # 2) calculate z_levels for all points
        # 3) construct rgi for (time, z_level, lat, lon)
        # 4) interpolate over z
        
        z_levels = np.squeeze(self.H(lat = lat, lon = lon)) # ilev has a default
#        print(z_levels)
        level = interp1d(z_levels, self._ilev_density, bounds_error=False)
        return level(z)

    def lev_density_log(self, z, lat, lon):
        """Finds ilev for a given height using log spacing to improve interpolation accuracy"""
        # 2) calculate z_levels for all points spaced logarithmiclly
        # 3) construct rgi for (time, z_level, lat, lon)
        # 4) interpolate over log(z)
        
        z_levels = np.squeeze(np.log(self.H(lat = lat, lon = lon))) # ilev has a default
#        print(z_levels)
        level = interp1d(z_levels, self._ilev_density, bounds_error=False)
        return level(np.log(z))

    @np.vectorize
    def lev_height(self, z, lat, lon):
# return input as height data is already using ht 
#        """Finds ilev for a given height"""
#        # 2) calculate z_levels for all points
#        # 3) construct rgi for (time, z_level, lat, lon)
#        # 4) interpolate over z
        return z
#        z_levels = np.squeeze(self.Z(t = t, lat = lat, lon = lon)) # ilev has a default
#        level = interp1d(z_levels, self._ilev_height, bounds_error=False)
#        return level(z)

    @np.vectorize
    def lev_neutral(self, z, lat, lon):
        """Finds ilev for a given height"""
        # 2) calculate z_levels for all points
        # 3) construct rgi for (time, z_level, lat, lon)
        # 4) interpolate over z
        
        z_levels = np.squeeze(self.Z(lat = lat, lon = lon)) # ilev has a default
        level = interp1d(z_levels, self._ilev_neutral, bounds_error=False)
        return level(z)

    def vert_interp(self,varname,z,lat,lon, interp4=False,log_height=True):
# interp4 scheme - generate four 1D height interpolations in lon, lat and combine results to returned value
# interp4=False - generate one 1D height interpolations using regular grid interpoator
        data=self[varname].data
        lons=self._lon_density
        lats=self._lat_density
        if log_height:
            # use log scaling in altitude
            H=np.log(self['H'].data)
            zlog=np.log(z)
        else:
            # use linear scaling in altitude
            H=self['H'].data
            zlog=z
            
        ilon=np.array(np.where(lons[0:-1] <= lon)).max()
        if ilon == len(lons):
            ilon=ilon-1
        w0_lon=(lons[ilon+1]-lon)/(lons[ilon+1]-lons[ilon])
        print('ILON=%d W0=%f' %( ilon, w0_lon))
            
        ilat=np.array(np.where(lats[0:-1] <= lat)).max()
        if ilat == len(lats):
            ilat=ilat-1
        w0_lat=(lats[ilat+1]-lat)/(lats[ilat+1]-lats[ilat])
        print('ILAT=%d W0=%f' %( ilat, w0_lat) )

        print('Lons: ',lons.shape)
        print('Lats: ',lats.shape)
        print('Height: ',H.shape)
        print('Data:' ,data.shape)

        if interp4:
            ilon1=ilon+1
            ilat1=ilat+1
            w1_lon=1.-w0_lon
            w1_lat=1.-w0_lat
        
            h_000=np.squeeze(H[:,ilat,ilon])               
            d_000=np.squeeze(data[:,ilat,ilon])
            int00=interp1d(h_00,d_00,bounds_error=False)
        
            h_01=np.squeeze(H[:,ilat,ilon1])               
            d_01=np.squeeze(data[:,ilat,ilon+1])
            int01=interp1d(h_01,d_01,bounds_error=False)
        
            h_10=np.squeeze(H[:,ilat1,ilon])               
            d_10=np.squeeze(data[:,ilat+1,ilon])
            int10=interp1d(h_10,d_10,bounds_error=False)
        
            h_11=np.squeeze(H[:,ilat1,ilon1])               
            d_11=np.squeeze(data[:,ilat+1,ilon+1])
            int11=interp1d(h_11,d_11,bounds_error=False)
        
            return(+int00(zlog)*w0_lon*w0_lat
                   +int01(zlog)*w1_lon*w0_lat
                   +int10(zlog)*w0_lon*w1_lat
                   +int11(zlog)*w1_lon*w1_lat
            )
        else:
# interp4=False
# use interpolater as implemented in IDL for each time step
            rgiH = RegularGridInterpolator((self._ilev_density,self._lat_density, self._lon_density), 
                                      H[:,:,:], bounds_error = False)
            rgiD = RegularGridInterpolator((self._ilev_density, self._lat_density, self._lon_density), 
                                      data[:,:,:], bounds_error = False)
            lon_rep=np.repeat(lon,len(self._ilev_density))
            lat_rep=np.repeat(lat,len(self._ilev_density))
            xvec=np.vstack([self._ilev_density,lat_rep,lon_rep]).T
#            print(xvec)

            h_0=np.squeeze(rgiH(xvec))
            d_0=np.squeeze(rgiD(xvec))
            int0=interp1d(h_0,d_0,bounds_error=False)
            
            return(int0(zlog)*w0_time)

#
# add a layer for lon=360 (copy of lon=0) to each 3D and 4D variable
#
    def wrap_density_variable(self, variable_name):
        if variable_name not in self.wrapped_density:
            variable = self._ctipe_density.variables[variable_name].__array__()
#	    variable = self.rootgrp.variables[variable_name].__array__()
            if len(variable.shape) == 3: # 3D variable
                self._ctipe_density.variables[variable_name] = np.concatenate((variable, variable[:,:,0:1]), axis = 2)
       	        self.wrapped_density.append(variable_name)

            if len(variable.shape) == 4: # 4D variable
       	        self._ctipe_density.variables[variable_name] = np.concatenate((variable, variable[:,:,:,0:1]), axis = 3)
                self.wrapped_density.append(variable_name)


    def wrap_height_variable(self, variable_name):
        if variable_name not in self.wrapped_height:
            variable = self._ctipe_height.variables[variable_name].__array__()
            if len(variable.shape) == 3: # 3D variable
                self._ctipe_height.variables[variable_name] = np.concatenate((variable, variable[:,:,0:1]), axis = 2)
                self.wrapped_height.append(variable_name)

            if len(variable.shape) == 4: # 4D variable
                self._ctipe_height.variables[variable_name] = np.concatenate((variable, variable[:,:,:,0:1]), axis = 3)
                self.wrapped_height.append(variable_name)


    def wrap_neutral_variable(self, variable_name):
        if variable_name not in self.wrapped_neutral:
            variable = self._ctipe_neutral.variables[variable_name].__array__()
            if len(variable.shape) == 3: # 3D variable
                if variable.shape[2] == (self._lon_neutral.shape[0]-1):
                    self._ctipe_neutral.variables[variable_name] = np.concatenate((variable, variable[:,:,0:1]), axis = 2)
                else:
# CTIPe efield variables no not need to be wrapped but
# need to be transposed from (time,lon,lat) to (time,lat,lon)
                    self._ctipe_neutral.variables[variable_name] = np.transpose(variable,[0,2,1])

                    self.wrapped_neutral.append(variable_name)

            if len(variable.shape) == 4: # 4D variable
                self._ctipe_neutral.variables[variable_name] = np.concatenate((variable, variable[:,:,:,0:1]), axis = 3)
                self.wrapped_neutral.append(variable_name)


    def wrap_density_3d_variables(self):
        """Wrap all 3D variables by longitude"""
        for variable_name, variable in list(self.rootgrp.variables.items()):
            if len(variable.shape) == 3: # 3D variable
                self.wrap_density_variable(variable_name)

    def wrap_neutral_3d_variables(self):
        """Wrap all 3D variables by longitude"""
        for variable_name, variable in list(self.rootgrp.variables.items()):
            if len(variable.shape) == 3: # 3D variable
                self.wrap_neutral_variable(variable_name)

    def wrap_height_3d_variables(self):
        """Wrap all 3D variables by longitude"""
        for variable_name, variable in list(self.rootgrp.variables.items()):
            if len(variable.shape) == 3: # 3D variable
                self.wrap_height_variable(variable_name)

    def wrap_longitude(self):
        """puts longitude from [min, max] to [min, wrapped(min)]"""
        self.longitude_min = self.lon.min()
        self.longitude_max = to_range(self.longitude_min, 0, 360)
        self.lon = np.hstack((self.lon, self.longitude_max))
        self.lonrange=dict(min=0,max=360,n=73)
        
    def wrap_e_longitude(self):
        """puts efield longitude from [min, max] to [min, wrapped(min)]"""
        self.elongitude_min = self.elon.min()
        self.elongitude_max = to_range(self.elongitude_min, 0, 360)
        self.elon = np.hstack((self.elon, self.elongitude_max))

    def set_plot(self,
                 plottype = "XY",
                 cutV = 10,
                 lonrange=dict(min=0,max=360,n=73), # reference to self._lon_density not possible?
                 latrange=dict(min=-90,max=90,n=37),
                 hrange=dict(min=1,max=15,n=15) ):
        '''Set plotting variables for available preset plot types: XY, YZ, YZ, XYZ'''
        if 'min' not in lonrange.keys():
            lonrange['min']=self.lonrange['min']
        if 'max' not in lonrange.keys():
            lonrange['max']=self.lonrange['max']
        if 'n' not in lonrange.keys():
            lonrange['n']=self.lonrange['n']
        if 'min' not in latrange.keys():
            latrange['min']=self.latrange['min']
        if 'max' not in latrange.keys():
            latrange['max']=self.latrange['max']
        if 'n' not in latrange.keys():
            latrange['n']=self.latrange['n']
        if 'min' not in hrange.keys():
            hrange['min']=self.hrange['min']
        if 'max' not in hrange.keys():
            hrange['max']=self.hrange['max']
        if 'n' not in hrange.keys():
            hrange['n']=self.hrange['n']
            
        tic = time.perf_counter()
        if plottype == self.plottype:
            if cutV == self.cutV:
               if lonrange['min'] == self.lonrange['min'] and lonrange['max'] == self.lonrange['max'] and lonrange['n'] == self.lonrange['n']:
                       if latrange['min'] == self.latrange['min'] and latrange['max'] == self.latrange['max'] and latrange['n'] == self.latrange['n']:
                               if hrange['min'] == self.hrange['min'] and hrange['max'] == self.hrange['max'] and hrange['n'] == self.hrange['n']:
                       
                                       print('Plottype (',plottype,') and cut value (',cutV,') are unchanged, returning.')
                                       return
                               
        self.lonrange=lonrange
        self.latrange=latrange
        self.hrange=hrange

# X=IP or Height
# Y=Lat
# Z=Lon
        if plottype == "LonLat":
            self.plottype = plottype
            self.cut = 'IP'
            self.cutV = cutV
            self.nX = 1
            self.nY = latrange['n']
            self.nZ = lonrange['n']
            self.newx = cutV
            self.newy = np.linspace(latrange['min'],latrange['max'],latrange['n'])
            self.newz = np.linspace(lonrange['min'],lonrange['max'],lonrange['n'])
            self.scene_aspectratio=dict(x=1,y=(latrange['max']-latrange['min'])/(lonrange['max']-lonrange['min']))
#            self.newy = np.linspace(-35., 35., self.nY)
#            self.newz = np.linspace(-20., 20., self.nZ)
        elif plottype == "LonIP":
            self.plottype = plottype
            self.cut = 'Lat'
            self.cutV = cutV
#            self.nX = 141
            self.nX = hrange['n']
            self.nY = 1
            self.nZ = lonrange['n']
#            self.nZ = 81
            self.newx=np.linspace(hrange['min'],hrange['max'],hrange['n'])
#            self.newx = np.linspace(-50., 20., self.nX)
            self.newy = cutV
            self.newz=np.linspace(lonrange['min'],lonrange['max'],lonrange['n'])
#            self.newz = np.linspace(-20., 20., self.nZ)
            self.scene_aspectratio=dict(x=1,y=(hrange['max']-hrange['min'])/(lonrange['max']-lonrange['min']))
        elif plottype == "LatIP":
            self.plottype = plottype
            self.cut = 'Lon'
            self.cutV = cutV
            self.nX=hrange['n']
            self.nY=latrange['n']
            self.nZ=1
            self.newx = np.linspace(hrange['min'],hrange['max'],hrange['n'])
            self.newy = np.linspace(latrange['min'],latrange['max'],latrange['n'])
            self.newz = cutV
            self.scene_aspectratio=dict(x=1,y=(hrange['max']-hrange['min'])/(lonrange['max']-lonrange['min']))
        elif plottype == "XYZ":
            self.plottype = plottype
            return
        else:
            print('Error, unknown plottype. ',plottype)
            return

        self.plots[plottype] = dict(cut=self.cut, cutV=self.cutV, tol=self.tol,
                                    nX=self.nX, nY=self.nY, nZ=self.nZ,
                                    newx=self.newx, newy=self.newy, newz=self.newz)
        self.setup_interpolating_grids()
#        self.reinterpolate_values()

        for varname in self.variables:
            self.plots[plottype][varname]=self.variables[varname]['interpolator']

        toc = time.perf_counter()
        print(f"Time resetting plot and precomputing interpolations: {toc - tic:0.4f} seconds")
        return
    
    def get_plot(self, var, colorscale="Viridis"):
        '''
        Return a plotly figure object for the available plot types set in set_plot()..
        colorscale = Viridis [default], Cividis, or Rainbow
        '''
        #Set some text strings
        txtbot="Model: CTIPe,  Run: " + str(self.runname) + ",  " + str(self.gridSize) + " cells,  minimum dx=" + str(self.gridMinDx)
        txtbar=var + " [" + self.variables[var]['units'] + "]"
        
        # Get values from interpolation already computed
        result=self.variables[var]['interpolator'](self.newgrid)
#        r = np.sqrt(np.square(self.newgrid[:,0]) + np.square(self.newgrid[:,1]) + np.square(self.newgrid[:,2]))
#        cmin=np.amin(result[(r[:] > 2.999)])
#        cmax=np.amax(result[(r[:] > 2.999)])
        
        if self.plottype == "LonLat":
            txttop=self.cut+"=" + str(self.plots[self.plottype]['cutV']) + " slice,  Time = " + self.filetime
            xint = self.newz
            yint = self.newy
            xunit=self.zunit
            yunit=self.yunit
            xrange=dict(min=xint.min(),max=xint.max(),n=len(xint))
            yrange=dict(min=yint.min(),max=yint.max(),n=len(yint))
            # Reshape interpolated values into 2D
            result2=np.reshape(result,(self.nY,self.nZ))
            def plot_XY(xint = xint, yint = yint):
                return result2
            plotXY = Kamodo(plot_XY = plot_XY)
            fig = plotXY.plot(plot_XY = dict( xsize=600,ysize=300))
            fig.update_xaxes(title_text="",scaleanchor='y')
            fig.update_yaxes(title_text="Lat [deg]")
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
#                zmin=cmin, zmax=cmax,
                ncontours=201,
                colorbar=dict(title=txtbar,tickformat='.4g'),
                    hovertemplate="X: %{x:.2f}<br>Y: %{y:.2f}<br><b> %{z:.4g}</b><extra></extra>",
                contours=dict(coloring="fill", showlines=False)
            )
            if xunit == yunit:
            # real aspect ratio
                if xrange['max']-xrange['min'] > yrange['max']-yrange['min']:
                    aspectratio=dict(x=np.asarray([4.,(xrange['max']-xrange['min'])/np.asarray([1.e-5,(yrange['max']-yrange['min'])]).max()]).min(),y=1)
                else:
                    aspectratio=dict(x=1,y=np.asarray([4,(yrange['max']-yrange['min'])/np.asarray([1.e-5,(xrange['max']-xrange['min'])]).max()]).min())
            else:
                aspectratio=dict(x=1,y=0.25)
#            print("aspect ratio: x=",aspectratio['x'],' y=',aspectratio['y'])
            
            width=180+400*aspectratio['x']
            height=90+400*aspectratio['y']
            fig.update_layout(
                autosize=False,
                height=height,
                width=width,
                margin=dict(t=45,b=45,l=40,r=140),
                title=dict(text=txttop),
                annotations=[
                    dict(text="Lon [deg]", x=0.5, y=-0.10, showarrow=False, xref="paper", yref="paper", font=dict(size=14)),
                    dict(text=txtbot,
                         font=dict(size=16, family="sans serif", color="#000000"),
                         x=0.0, y=0.0, ax=0, ay=0, xanchor="left", xshift=-65, yshift=-42, xref="paper", yref="paper"
                    )
                ],
                xaxis = dict(
                    tickmode = 'linear',
                    tick0 = 0,
                    dtick = 30
                ),
                yaxis = dict(
                    tickmode = 'linear',
                    tick0 = 0,
                    dtick = 30
                )
            )
            if self.plots[self.plottype]['cutV'] == 0.:
                fig.update_layout(
                    shapes=[
                        dict(type="circle", xref="x", yref="y", x0=-3, y0=-3, x1=3, y1=3, fillcolor="black", line_color="black"),
                        dict(type="circle", xref="x", yref="y", x0=-1, y0=-1, x1=1, y1=1, fillcolor="black", line_color="white"),
                        dict(type="path", path= self.ellipse_arc(N=30), fillcolor="white", line_color="white")
                    ]
                )
            return fig


        if self.plottype == "LonIP":
            txttop=self.cut+"=" + str(self.plots[self.plottype]['cutV']) + " slice,  Time = " + self.filetime
            xint = self.newz
            yint = self.newx
            xunit=self.zunit
            yunit=self.xunit
            xrange=dict(min=xint.min(),max=xint.max(),n=len(xint))
            yrange=dict(min=yint.min(),max=yint.max(),n=len(yint))

            # Reshape interpolated values into 2D
            result2=np.reshape(result,(self.nX,self.nZ))
            def plot_XZ(xint = xint, yint = yint):
                return result2
            plotXZ = Kamodo(plot_XZ = plot_XZ)
            fig = plotXZ.plot(plot_XZ = dict())
            fig.update_xaxes(title_text="",scaleanchor='y')
            fig.update_yaxes(title_text="IP []")
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
#                zmin=cmin, zmax=cmax,
                ncontours=201,
                colorbar=dict(title=txtbar,tickformat='.4g'),
                hovertemplate="Lon: %{x:.2f}<br>IP: %{y:.2f}<br><b> %{z:.4g}</b><extra></extra>",
                contours=dict(coloring="fill", showlines=False)
            )
            if xunit == yunit:
            # real aspect ratio
                if xrange['max']-xrange['min'] > yrange['max']-yrange['min']:
                    aspectratio=dict(x=np.asarray([4.,(xrange['max']-xrange['min'])/np.asarray([1.e-5,(yrange['max']-yrange['min'])]).max()]).min(),y=1)
                else:
                    aspectratio=dict(x=1,y=np.asarray([4,(yrange['max']-yrange['min'])/np.asarray([1.e-5,(xrange['max']-xrange['min'])]).max()]).min())
            else:
                aspectratio=dict(x=1,y=0.25)
#            print("aspect ratio: x=",aspectratio['x'],' y=',aspectratio['y'])
            
            width=180+400*aspectratio['x']
            height=90+400*aspectratio['y']
            fig.update_layout(
                autosize=False,
                height=height,
                width=width,
                margin=dict(t=45,b=45,l=40,r=140),
                title_text=txttop,
                annotations=[
                    dict(text="Lon [deg]", x=0.5, y=-0.10, showarrow=False, xref="paper", yref="paper", font=dict(size=14)),
                    dict(text=txtbot,
                         font=dict(size=16, family="sans serif", color="#000000"),
                         x=0.0, y=0.0, ax=0, ay=0, xanchor="left", xshift=-65, yshift=-42, xref="paper", yref="paper"
                    )
                ]
            )
            return fig;
            
        if self.plottype == "LatIP":
            txttop=self.cut+"=" + str(self.plots[self.plottype]['cutV']) + " slice,  Time = " + self.filetime
            xint = self.newy
            yint = self.newx
            xunit=self.yunit
            yunit=self.xunit
            xrange=dict(min=xint.min(),max=xint.max(),n=len(xint))
            yrange=dict(min=yint.min(),max=yint.max(),n=len(yint))
            # Reshape interpolated values into 2D
            result2=np.transpose(np.reshape(result,(self.nY,self.nX)))
            def plot_XY(xint = xint, yint = yint):
                return result2
            plotXY = Kamodo(plot_XY = plot_XY)
            fig = plotXY.plot(plot_XY = dict())
            fig.update_xaxes(title_text="",scaleanchor='y')
            fig.update_yaxes(title_text="IP []")
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
#                zmin=cmin, zmax=cmax,
                ncontours=201,
                colorbar=dict(title=txtbar,tickformat='.4g'),
                hovertemplate="Lat: %{x:.2f}<br>IP: %{y:.2f}<br><b> %{z:.4g}</b><extra></extra>",
                contours=dict(coloring="fill", showlines=False)
            )
            if xunit == yunit or yunit == '':
            # real aspect ratio
                if xrange['max']-xrange['min'] > yrange['max']-yrange['min']:
                    aspectratio=dict(x=np.asarray([4.,(xrange['max']-xrange['min'])/np.asarray([1.e-5,(yrange['max']-yrange['min'])]).max()]).min(),y=1)
                    if aspectratio['x'] > 2:
                        aspectratio_x=aspectratio['x']
                        aspectratio=dict(x=2., y=aspectratio['y']*2./aspectratio_x)
                else:
                    aspectratio=dict(x=1,y=np.asarray([4,(yrange['max']-yrange['min'])/np.asarray([1.e-5,(xrange['max']-xrange['min'])]).max()]).min())
            else:
                aspectratio=dict(x=1,y=0.25)
            print('xrange: ',xrange)
            print('yrange: ',yrange)
            print("aspect ratio: x=",aspectratio['x'],' y=',aspectratio['y'])

            width=180+400*aspectratio['x']
            height=90+400*aspectratio['y']
            fig.update_layout(
                autosize=False,
                height=height,
                width=width,
                margin=dict(t=45,b=45,l=40,r=140),
                title_text=txttop,
                annotations=[
                    dict(text="Lat [deg]", x=0.5, y=0.0, showarrow=False, xref="paper", yref="paper", yshift=-32,font=dict(size=14)),
                    dict(text=txtbot,
                         font=dict(size=16, family="sans serif", color="#000000"),
                         x=0.0, y=0.0, ax=0, ay=0, xanchor="left", xshift=-45, yshift=-42, xref="paper", yref="paper"
                    )
                ]
            )
            return fig


        print("ERROR, unknown plottype =",plottype,", exiting without definition of figure.")
        return
        
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
    


