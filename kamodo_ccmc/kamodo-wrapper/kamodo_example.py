from kamodo import Kamodo, kamodofy, gridify

import numpy as np
from scipy.interpolate import RegularGridInterpolator


class MyModel(Kamodo):
    def __init__(self, filename, **kwargs):
        # perform any necessary I/O
        self.filename = filename
        self.missing_value = np.nan

        # store any data needed for interpolation
        self.x = np.linspace(1, 4, 4)
        self.y = np.linspace(4, 7, 4)
        self.z = np.linspace(7, 9, 3)
        xx, yy, zz = np.meshgrid(self.x, self.y, self.z, indexing='ij', sparse=True)
        density_data = 2 * xx**3 + 3 * yy**2 - zz


        self.variables = dict(rho = dict(units = 'kg/m**3', data = density_data))

        # Prepare model for function registration
        super(MyModel, self).__init__(**kwargs)

        for varname in self.variables:
            units = self.variables[varname]['units']
            self.register_variable(varname, units)

    def register_variable(self, varname, units):
        interpolator = self.get_grid_interpolator(varname)

        # store the interpolator
        self.variables[varname]['interpolator'] = interpolator

        def interpolate(xvec):
            return self.variables[varname]['interpolator'](xvec)

        # update docstring for this variable
        interpolate.__doc__ = "A function that returns {} in [{}].".format(varname,units)

        self[varname] = kamodofy(interpolate,
                           units = units,
                           citation = "Pembroke et al 2019",
                          data = None)
        self[varname + '_ijk'] = kamodofy(gridify(self[varname],
                                                  x_i = self.x,
                                                  y_j = self.y,
                                                  z_k = self.z, squeeze=False),
                            units = units,
                            citation = "Pembroke et al 2019",
                            data = self.variables[varname]['data'])


    def get_grid_interpolator(self, varname):
        """create a regulard grid interpolator for this variable"""
        data =  self.variables[varname]['data']

        interpolator = RegularGridInterpolator((self.x, self.y, self.z), data,
                                                bounds_error = False,
                                               fill_value = self.missing_value)
        return interpolator



