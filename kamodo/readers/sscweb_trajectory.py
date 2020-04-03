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
from datetime import datetime,timezone,timedelta
import time
from util import *
import sscweb_rest as sscweb

class SSCWEB_Trajectory(Kamodo):
    def __init__(self,satellite_name,start_date,end_date,**kwargs):
        print('Obtaining trajectory from SSCWeb for %s' % satellite_name)
        self.satellite_name=satellite_name
        self.missing_value=np.NAN
        self.verbose=False
        self.symbol_registry=dict();
        self.signatures=dict();
        import sscweb_rest as sscweb
        (start_year,start_month,start_day,start_hour,start_minute,start_second)=start_date
        (end_year,end_month,end_day,end_hour,end_minute,end_second)=end_date
        
        StartTime=datetime(start_year,start_month,start_day,start_hour,start_minute,start_second,tzinfo=timezone.utc)
        EndTime=datetime(end_year,end_month,end_day,end_hour,end_minute,end_second,tzinfo=timezone.utc)

#        import pickle # json does not work with DateTime objects
#        observatory_file_dir='/Users/lrastaet/Kamodo_data/SSCWEB/'
#        if not os.path.isdir(observatory_file_dir):
#            os.makedirs(observatory_file_dir,exist_ok=True)
            
#        observatory_file="/Users/lrastaet/Kamodo_data/SSCWEB/SSCWEB_observatories.pyc"
#        if os.path.isfile(observatory_file):
#            f=open(observatory_file,'rb')
#            observatories=pickle.load(f)
#            f.close()
#        else:
#            f=open(observatory_file,'wb')
#            observatories=sscweb.get_observatories()
#            pickle.dump(observatories,f)
#            f.close()
            
#        for observatory in observatories:
#            if observatory['Id'] == satellite_name or observatory['Name'] == satellite_name:
#                self.observatory=observatory
#        if observatory:
        orbit_data=sscweb.sample_orbit(satellite_name,StartTime,EndTime)
        self.orbit_data=orbit_data
        
        self.X=np.array(orbit_data['Data'][0]['Coordinates']['X']) #
        self.Y=np.array(orbit_data['Data'][0]['Coordinates']['Y']) # positions im km in GSE
        self.Z=np.array(orbit_data['Data'][0]['Coordinates']['Z']) #
        RE=6371.200 # R_E in km
        self.Time=np.array(orbit_data['Data'][0]['Time']) 
        self.time = list(self.seconds_from_20000101(self.Time))
        # need to convert X,Y,Z to vector stack and report current coordinate system and provide transformation functions to generate other coordinatesd and also long, lat altitude if needed
        
        # time interpolation function  Xvec(Time)
        self.variables=dict()
        self.variables['X']=dict(data=self.X,units='km',coordinate_system='GSE')
        self.variables['Y']=dict(data=self.Y,units='km',coordinate_system='GSE')
        self.variables['Z']=dict(data=self.Z,units='km',coordinate_system='GSE')
        self.variables['Xvec']=dict(data=np.vstack([self.X,self.Y,self.Z]),units='km',coordinate_system='GSE')

        for varname in self.variables:
            unit=self.variables[varname]['units']
            data=self.variables[varname]['data']
            print("variable:",varname,"unit:",unit,"data.size:",data.size)
            self.register_variable(varname,unit)

    def seconds_from_20000101(self,t):
        return list(map(lambda x: (x-datetime(2000,1,1,0,0,0,tzinfo=timezone.utc)).total_seconds(),t))
    
#    def time_interpolation(self,df):
#   # combine original time index with new times
#        df_ = df.reindex(df.index.union(t))
#   # apply pandas' interpolation method to fill in new data
#        df_interpolated = df_.interpolate()
#   # return only values at selected times
#        return df_interpolated.reindex

    def Xvec_interpolator(self,t):
        X_=self['X'](t)
        Y_=self['Y'](t)
        Z_=self['Z'](t)
        return np.vstack([X_,Y_,Z_]).T
    
    def get_grid_interpolator(self, varname):
        """create a regular grid interpolator for this variable"""
        if varname in ['Xvec','xvec','XVEC']:
            interpolator=self.Xvec_interpolator
        else:
            data =  self.variables[varname]['data']
            interpolator = RegularGridInterpolator(([self.time]),
                                                   data, 
                                                   bounds_error = False,
                                                   fill_value = self.missing_value)
            
        return interpolator

      
    def register_variable(self, varname, units):
        interpolator = self.get_grid_interpolator(varname)
#        var_data=pd.DataFrame(data=self.variables[varname]['data'],index=self.seconds_from_20000101(self.Time))
#        interpolator=self.time_interpolation(var_data)

# store the interpolator
        self.variables[varname]['interpolator'] = interpolator

        def interpolate(t):  
            return self.variables[varname]['interpolator'](t)

        # update docstring for this variable
        interpolate.__doc__ = "A function that returns {} in [{}].".format(varname,units)

        self[varname] = kamodofy(interpolate, 
                                 units = units, 
                                 citation = "Pembroke et al 2019, Rastaetter 2020 SSCWeb in Kamodo",
                                 data = None)
        
