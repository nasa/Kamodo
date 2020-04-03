from kamodo import Kamodo, kamodofy
#, gridify
from netCDF4 import Dataset
import numpy as np
import os
import scipy
from scipy import interpolate
from scipy.interpolate import RegularGridInterpolator
import math
import scipy.constants as constants

import time
#from util import time_in_interval
from util import *
#from util import boundary_conditions, fill_masked

class SWMF_IE(Kamodo):
    def __init__(self,filename,**kwargs):
        print('opening SWMF Ionosphere Electrodynamics file %s' % filename)
        self.filename = filename
        self.missing_value=np.NAN
        data=self.read_swmf_ie(filename)
        self.variables=dict()
        self.coords=['theta','phi','altitude']
        self.coord_units=['deg','deg','km']
        self.x=data['Theta']
        self.y=data['Psi']
        self.z=[110.]
#        self.coordinate_system='SM'
# these three variables are needed to prevent errors in Kamodo      
        self.verbose=False;
        self.symbol_registry=dict();
        self.signatures=dict();
       
        #        self.z_, self.y_, self.x_ = scipy.meshgrid(self.z, self.y, self.x, indexing = 'ij')
        self.y_, self.x_ = scipy.meshgrid(self.y, self.x, indexing = 'ij')

#        print(data['orig_units'])
#        print(data['units'])
        for ivar in range(len(data['variables'])):
            varname=data['variables'][ivar]
            var_data=data['data'][:,:,ivar].squeeze().T
            unit=data['units'][ivar]
#            print("variable:",varname,"unit:",unit,"data:",var_data)

            self.variables[varname]=dict(units=unit, data=var_data)
            self.register_variable(varname,unit)
            
    def read_swmf_ie(self,filename):
        import re
        arrays = []
        orig_units=[]
        units=[]
        variables=[]
        iline=0
        with open(filename, 'r') as a:
            for line in a.readlines():
                A = re.match(r'TITLE=(.*$)', line, re.M | re.I)
                B = re.match(r'(\s*)VARIABLES=(.*$)', line, re.M | re.I)
                C = re.match(r'(\s*)(\")(.*$)', line,re.M | re.I)
                D = re.match(r'(\s*)ZONE (.*$)', line, re.M | re.I)
                E = re.match(r'(\s*)I=(.*$)', line, re.M | re.I)
                if A or B or C or D or E:
                    if B or C:
                        for s in (line.split('"'))[1:]:
                            if s != "," and s != '':
                                if len(s) >=3:
                                    (var,unit) = s.split()
# cannot have "-" sign in veriable names
                                    var=var.replace("-","")
                                    variables.append(var)
                            # map unit names to something SymPy may understand
                                    unit=unit[1:-1]
                                    orig_unit=unit
                                    if unit == 'R':
#                                        unit='6371200 m'
                                        unit='R_E'
                                    if unit[0:2] == '`m':
                                        unit="mu%s" % unit[2:]
                                    if unit[0] == '/':
                                        unit="1%s" % unit
                                    unit=unit.replace('m2','m**2')
                                    unit=unit.replace('m^2','m**2')
                                    orig_units.append(orig_unit)
                                    units.append(unit)
                    if E:
                        index_i=line.index("I=")+2
                        index_j=line.index("J=")+2
                        NI=int(line[index_i:].split()[0])
                        NJ=int(line[index_j:].split()[0])
                        continue
                else:
                    for s in line.split():
                        arrays.append(float(s))

        nvar=len(variables)
        nelements=len(arrays)
        npos=int(nelements/nvar)           
        arrays = np.array(arrays)
    
        arrays=arrays.reshape((npos,nvar))
        arrays_N=arrays[0:int(npos/2),:].reshape((NJ,NI,nvar))
        arrays_S=arrays[int(npos/2):,:].reshape((NJ,NI,nvar))
        data=np.concatenate((arrays_N,arrays_S[:,1:,:]),axis=1)
        data[:,0,3]=0.
        data[:,-1,3]=180.
        df={'data':data,
            'variables':variables,
            'orig_units':orig_units,
            'units':units,
            'axis1':'Theta',
            'axis2':'Psi',
            'Theta': data[0,:,3].flatten(),
            'Psi': data[:,0,4].flatten()
        }
    
        return df

    def get_grid_interpolator(self, varname):
        """create a regular grid interpolator for this variable"""
        data =  self.variables[varname]['data']
        interpolator = RegularGridInterpolator((self.x, self.y),
                                               data, 
                                               bounds_error = False,
                                               fill_value = self.missing_value)
        return interpolator
        
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
                                 citation = "Pembroke et al 2019, Rastaetter 2020",
                                 data = None)
        self[varname + '_ij'] = kamodofy(gridify(self[varname], 
                                                  x_i = self.x, 
                                                  y_j = self.y),
                                          units = units,
                                          citation = "Pembroke et al 2019, Rastaetter 2020",
                                          data = self.variables[varname]['data'])
        
