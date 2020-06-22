import numpy as np

from kamodo import Kamodo, kamodofy, gridify

from scipy.interpolate import RegularGridInterpolator, interp1d
from netCDF4 import Dataset


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

class CTIPe_Kamodo(Kamodo):
    def __init__(self, filename, variables = None, debug= False, **kwargs):
        import re
# input file name can be one of the 4 files for each day of model outputs
#  YYYYMMDD-plot-[density|height|neutral|plasma].nc files
        if filename.find('density') > 0:
            filename_density=filename
            filename_height=filename.replace('density','height')
#            filename_height.replace('density','height')
            filename_neutral=filename.replace('density','neutral')
#            filename_neutral.replace('density','neutral')
            filename_plasma=filename.replace('density','plasma')
#            filename_plasma.replace('density','plasma')

        if filename.find('height') > 0:
            filename_density=filename
            filename_density.replace('height','density')
            filename_height=filename
            filename_neutral=filename
            filename_neutral.replace('height','neutral')
            filename_plasma=filename
            filename_plasma.replace('height','plasma')

        if filename.find('neutral') > 0:
            filename_density=filename
            filename_density.replace('neutral','density')
            filename_height=filename
            filename_neutral.replace('neutral','height')
            filename_neutral=filename
            filename_plasma=filename
            filename_plasma.replace('neutral','plasma')

        if filename.find('plasma') > 0:
            filename_density=filename
            filename_density.replace('plasma','density')
            filename_height=filename
            filename_height.replace('plasma','height')
            filename_neutral=filename
            filename_neutral.replace('plasma','neutral')
            filename_plasma=filename
# only the density, height and neutral files have data and are read
        print(filename_density)
        print(filename_height)
        print(filename_neutral)
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
        
        self._ilev_density = np.array(self._ctipe_density.variables['plev'])
        self._lat_density = np.array(self._ctipe_density.variables['lat'])
        self._lon_density = np.array(self._ctipe_density['lon'])

#        print(self._ctipe_height.variables)

        self._ht_height = np.array(self._ctipe_height.variables['ht'])
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
        
        super(CTIPe_Kamodo, self).__init__() 
        
        if variables is None:
            density_variables = self._ctipe_density.variables.keys()
            height_variables  = self._ctipe_height.variables.keys()
            neutral_variables = self._ctipe_neutral.variables.keys()
        else:
            density_variables = []
            height_variables  = []
            neutral_variables = []
            for variable in variables:
                if variable in self._ctipe_density.variables.keys():
                    density_variables.append(variable)
                if variable in self._ctipe_height.variables.keys():
                    height_variables.append(variable)
                if variable in self._ctipe_neutral.variables.keys():
                    neutral_variables.append(variable)
                    
            
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
                self.register_4d_density_variable(units, variable, varname)

            elif len(variable.shape) == 3:
                self.register_3d_density_variable(units, variable, varname)


        for varname in height_variables:
            print(varname)
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
                self.register_4d_height_variable(units, variable, varname)

            elif len(variable.shape) == 3:
                self.register_3d_height_variable(units, variable, varname)

                
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
                print(variable)
                print(units)
            
            if varname == 'ZMAG':
                continue

            elif len(variable.shape) == 4:
                self.register_4d_neutral_variable(units, variable, varname)

            elif len(variable.shape) == 3:
                self.register_3d_neutral_variable(units, variable, varname)

                
        print('registered {} variables'.format(self._registered))
        
        # register user's input variables, assuming kamodo-compatible
        for varname, variable in kwargs.items():
            self[varname] = variable

#        def ilev(points):
#            return self.lev(*points)

#        self['ilev'] = ilev

    def register_3d_density_variable(self, units, variable, varname):
        """Registers a 3d interpolator with 3d signature"""

        rgi = RegularGridInterpolator((self._time, self._lat_density, self._lon_density),
                                      variable, bounds_error = False)
        @kamodofy(units = units, data = variable)
        @gridify(t = self._time, lat = self._lat_density, lon = self._lon_density)
        def interpolator(xvec):
            """Interpolates 3d density variable"""
            return rgi(xvec)

        self[varname] = interpolator
        self._registered += 1

    def register_4d_density_variable(self, units, variable, varname):
        """Registers a 4d interpolator with 4d signature"""
#        print(varname,self._time.shape,self._ilev_density.shape,
#              self._lat_density.shape, self._lon_density.shape, 
#              variable.shape)
            
#        print(self._lon_density)

        rgi = RegularGridInterpolator((self._time, self._ilev_density, self._lat_density, self._lon_density), 
                                      variable, bounds_error = False)
        @kamodofy(units = units, data = variable)
        @gridify(t = self._time, ilev = self._ilev_density, lat = self._lat_density, lon = self._lon_density) 
        def interpolator(xvec):
            """Interpolates 4d density variable"""
            return rgi(xvec)

        try:
            self[ctipe_kamodo_variable_name(varname)] = interpolator
        except:
            self[varname] = interpolator
            
        self._registered += 1

    def register_3d_height_variable(self, units, variable, varname):
        """Registers a 3d interpolator with 3d signature"""

        rgi = RegularGridInterpolator((self._time, self._lat_height, self._lon_height),
                                      variable, bounds_error = False)
        @kamodofy(units = units, data = variable)
        @gridify(t = self._time, lat = self._lat_height, lon = self._lon_height)
        def interpolator(xvec):
            """Interpolates 3d height variable"""
            return rgi(xvec)


        try:
            self[ctipe_kamodo_variable_name(varname)] = interpolator
        except:
            self[varname] = interpolator            

        self._registered += 1

    def register_4d_height_variable(self, units, variable, varname):
        """Registers a 4d interpolator with 4d signature"""

        rgi = RegularGridInterpolator((self._time, self._ht_height, self._lat_height, self._lon_height), 
                                      variable, bounds_error = False)
        @kamodofy(units = units, data = variable)
        @gridify(t = self._time, ilev = self._ht_height, lat = self._lat_height, lon = self._lon_height) 
        def interpolator(xvec):
            """Interpolates 4d  height variable"""
            return rgi(xvec)

        try:
            self[ctipe_kamodo_variable_name(varname)] = interpolator
        except:
            self[varname] = interpolator            

        self._registered += 1

    def register_3d_neutral_variable(self, units, variable, varname):
        """Registers a 3d interpolator with 3d signature"""
        if (variable.shape)[1] == self._elat_neutral.shape[0]:
            rgi = RegularGridInterpolator((self._time, self._elat_neutral, self._elon_neutral),
                                      variable, bounds_error = False)
            @kamodofy(units = units, data = variable )
# to use the density grid one may have to run interpolation to provide data for the above kamodofy decorator
#           variable_on_lat_lon_neutral=rgi(self._time,???,???,????)
#            @kamodofy(units = units, data = variable_on_lat_lon_neutral )
#            @gridify(t = self._time, lat = self._lat_neutral, lon = self._lon_neutral)    
            @gridify(t = self._time, lat = self._elat_neutral, lon = self._elon_neutral)    
            def interpolator(xvec):
                """Interpolates 3d neutral variable"""
                return rgi(xvec)
        else:
            rgi = RegularGridInterpolator((self._time, self._lat_neutral, self._lon_neutral),
                                      variable, bounds_error = False)
            @kamodofy(units = units, data = variable)
            @gridify(t = self._time, lat = self._lat_neutral, lon = self._lon_neutral)    
            def interpolator(xvec):
                """Interpolates 3d neutral variable"""
                return rgi(xvec)

        try:
            self[ctipe_kamodo_variable_name(varname)] = interpolator
        except:
            self[varname] = interpolator
        self._registered += 1

    def register_4d_neutral_variable(self, units, variable, varname):
        """Registers a 4d interpolator with 4d signature"""
#        print(varname,self._time.shape,self._ilev_density.shape,
#              self._lat_density.shape, self._lon_density.shape, 
#              variable.shape)
        if variable.shape[-1] < self._lon_density.shape[0]:
            self.wrap_neutral_variable(varname)

        rgi = RegularGridInterpolator((self._time, self._ilev_neutral, self._lat_neutral, self._lon_neutral), 
                                      variable, bounds_error = False)
        @kamodofy(units = units, data = variable)
        @gridify(t = self._time, ilev = self._ilev_neutral, lat = self._lat_neutral, lon = self._lon_neutral) 
        def interpolator(xvec):
            """Interpolates 4d neutral variable"""
            return rgi(xvec)

        try:
            self[ctipe_kamodo_variable_name(varname)] = interpolator
        except:
# this alternative is invoked when the Kamodo variable name causes problems
            self[varname] = interpolator

        self._registered += 1


    @np.vectorize
    def lev_density(self, t, z, lat, lon):
        """Finds ilev for a given height"""
        # 2) calculate z_levels for all points
        # 3) construct rgi for (time, z_level, lat, lon)
        # 4) interpolate over z
        
        z_levels = np.squeeze(self.H(t = t, lat = lat, lon = lon)) # ilev has a default
#        print(z_levels)
        level = interp1d(z_levels, self._ilev_density, bounds_error=False)
        return level(z)

    @np.vectorize
    def lev_height(self, t, z, lat, lon):
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
    def lev_neutral(self, t, z, lat, lon):
        """Finds ilev for a given height"""
        # 2) calculate z_levels for all points
        # 3) construct rgi for (time, z_level, lat, lon)
        # 4) interpolate over z
        
        z_levels = np.squeeze(self.Z(t = t, lat = lat, lon = lon)) # ilev has a default
        level = interp1d(z_levels, self._ilev_neutral, bounds_error=False)
        return level(z)

#
# add a layaer for lon=360 (copy of lon=0) to each 3D and 4D variable
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

    def wrap_e_longitude(self):
        """puts efield longitude from [min, max] to [min, wrapped(min)]"""
        self.elongitude_min = self.elon.min()
        self.elongitude_max = to_range(self.elongitude_min, 0, 360)
        self.elon = np.hstack((self.elon, self.elongitude_max))






