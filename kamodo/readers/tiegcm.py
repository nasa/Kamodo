import numpy as np

from kamodo import Kamodo, kamodofy, gridify

from scipy.interpolate import RegularGridInterpolator, interp1d

# pip install pytiegcm
from tiegcm.tiegcm import TIEGCM

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
            varshape = self._tiegcm.rootgrp.variables[varname].shape
        except:
            continue
        if len(varshape) == 4:
            variables_4d.append(varname)
        elif len(varshape) == 3:
            variables_3d.append(varname)
    return variables_4d + variables_3d

class TIEGCM_Kamodo(Kamodo):
    def __init__(self, filename, variables = None, **kwargs):
        self._tiegcm = TIEGCM(filename)
        
        self._ilev = np.array(self._tiegcm.rootgrp.variables['ilev'])
        self._time = np.array(self._tiegcm.rootgrp.variables['time'])

        self._lat = self._tiegcm.lat
        self._lon = self._tiegcm.lon
        self._registered = 0
        
        super(TIEGCM_Kamodo, self).__init__() 
        
        if variables is None:
            variables = self._tiegcm.rootgrp.variables.keys()
        
        variables = sort_variables(self, variables)

        for varname in variables:
            
            units = parse_units(varname, self._tiegcm.rootgrp.variables[varname])
            
            try:
                self._tiegcm.set_variable_boundary_condition(varname)
            except:
#                 print('can not set boundary condition for {}, skipping..'.format(varname))
                continue
            self._tiegcm.wrap_variable(varname)
            variable = self._tiegcm.rootgrp.variables[varname]

            if len(variable.shape) not in [3,4]:
                # skip this variable
                continue
                
            if varname == 'ZMAG':
                continue

            elif len(variable.shape) == 4:
                self.register_4d_variable(units, variable, varname)

            elif len(variable.shape) == 3:
                self.register_3d_variable(units, variable, varname)

                
        print('registered {} variables'.format(self._registered))
        
        # register user's input variables, assuming kamodo-compatible
        for varname, variable in kwargs.items():
            self[varname] = variable


        def ilev(points):
            return self.lev(*points)

        self['ilev'] = ilev


    def register_3d_variable(self, units, variable, varname):
        """Registers a 3d interpolator with 3d signature"""

        rgi = RegularGridInterpolator((self._time, self._lat, self._lon),
                                      variable, bounds_error = False)
        @kamodofy(units = units, data = variable)
        @gridify(t = self._time, lat = self._lat, lon = self._lon)
        def interpolator(xvec):
            """Interpolates 3d variable"""
            return rgi(xvec)

        self[varname] = interpolator
        self._registered += 1

    def register_4d_variable(self, units, variable, varname):
        """Registers a 4d interpolator with 4d signature"""

        rgi = RegularGridInterpolator((self._time, self._ilev, self._lat, self._lon), 
                                      variable, bounds_error = False)
        @kamodofy(units = units, data = variable)
        @gridify(t = self._time, ilev = self._ilev, lat = self._lat, lon = self._lon) 
        def interpolator(xvec):
            """Interpolates 4d variable"""
            return rgi(xvec)
        self[varname] = interpolator
        self._registered += 1

    @np.vectorize
    def lev(self, t, z, lat, lon):
        """Finds ilev for a given height"""
        # 2) calculate z_levels for all points
        # 3) construct rgi for (time, z_level, lat, lon)
        # 4) interpolate over z
        
        z_levels = np.squeeze(self.Z(t = t, lat = lat, lon = lon)) # ilev has a default
        level = interp1d(z_levels, self._ilev, bounds_error=False)
        return level(z)






