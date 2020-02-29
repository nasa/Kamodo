#import t89
import numpy as np
from geopack import  geopack
from geopack import t89,t96,t01,t04
import os
import datetime
from kamodo import Kamodo, kamodofy,gridify,get_defaults
import scipy

#
# Initialization of Tsyganenko models:
# - date and time information to initialize geopack (recalc) and to obtain data:
#  T98: Kp value from hourly OMNI
#  T96: Dst (SYM_H), Solar wind Pdyn, IMF By, IMF Bz from 1-minite MNI
#  T01: Dst (SYM_H), Solar wind Pdyn, IMF By, IMF Bz from 1-minute OMNI
#  T04: Dst (SYM_H), Solar wind Pdyn, IMF By, IMF Bz, G1,G2,G3, W1-W6 from OMNI, Qin-Denton (daily files since 1995 prepared at CCMC, avilable via iSWA)
# 
#
class T89(Kamodo):
#
# using Sheng Tian's geopack module (https://github.com/tsssss/geopack)
#
    def __init__(self,year,month,day,hour,minute,use_igrf,*args,**kwargs):
        from geopack import t89
        super(T89, self).__init__(*args, **kwargs)
# time since Jan. 1 1970 00:00 UT as datetime.timedelta
        dt=datetime.datetime(year,month,day,hour,minute)-datetime.datetime(1970,1,1)
# seconds from 1970/1/1 00:00 UT
        self.dt_seconds=dt.total_seconds()
        
        self.ps = geopack.recalc(self.dt_seconds)
        self.use_igrf=use_igrf

        from geospacepy import omnireader
        sTimeIMF = datetime.datetime(year,month,day,hour)
        eTimeIMF = datetime.datetime(year,month,day,hour)+datetime.timedelta(0,0,0,1,0,0)
        omniInt = omnireader.omni_interval(sTimeIMF,eTimeIMF,'hourly')
        t = omniInt['Epoch'] #datetime timestamps
        By,Bz = omniInt['BY_GSM'],omniInt['BZ_GSM']

        kp=omniInt['KP']
        self.iopt=int(kp[0])+1
        if self.iopt > 7: self.iopt=7
        
        bounds_error = kwargs.get('bounds_error', False)
        fill_value = kwargs.get('missing_value', np.nan)        

        self.citation='Kamodo.T89 by Lutz Rastaetter (2020), Geopack/Tsyganenko by Sheng Tian (2019) and geospacepy-lite by Liam Kilkommons (2019)'
        self.unit='nT'

        self.x=np.linspace(-30.,10.,20)
        self.y=np.linspace(-10.,10.,10)
        self.z=np.linspace(-10.,10.,10)

        self.variables = dict(b_x = dict(units = 'nT', data = None),
                              b_y = dict(units = 'nT', data = None),
                              b_z = dict(units = 'nT', data = None),
                              bvec = dict(units = 'nT', data = None) )

        for varname in self.variables:
            units = self.variables[varname]['units']
            self.register_variable(varname, units)        

    def register_variable(self,varname,units):
        if varname == 'b_x':
            interpolator=self.bx
        if varname == 'b_y':
            interpolator=self.by
        if varname == 'b_z':
            interpolator=self.bz
        if varname == 'bvec':
            interpolator=self.b
        if varname == 'PF':
            interpolator=self.pressure_function
            
        self.variables[varname]['interpolator']= interpolator

        def interpolate(xvec):
            return self.variables[varname]['interpolator'](xvec)

        interpolate.__doc__ = "A function that returns {} in [{}].".format(varname,units)

        self[varname] = kamodofy(interpolate, 
                                 units = units, 
                                 citation = self.citation,
                                 data = None)
        self[varname + '_ijk'] = kamodofy(gridify(self[varname], 
                                                  x_i = self.x, 
                                                  y_j = self.y, 
                                                  z_k = self.z),
                                          units = units,
                                          citation = self.citation,
                                          data = None) 
        
    def bx(self,xvec):
        bx_,by_,bz_=np.hsplit(self.b(xvec),3)
        return(bx_)

    def by(self,xvec):
        bx_,by_,bz_=np.hsplit(self.b(xvec),3)
        return(by_)

    def bz(self,xvec):
        bx_,by_,bz_=np.hsplit(self.b(xvec),3)
        return(bz_)

    def b(self,xvec):
# x,y,z can be an array or list        
        try:
            x,y,z = xvec
        except: # assume nd array
            x,y,z = xvec.T
# we need to call recalc since common block is shared between instances
# of geopack_2008 and T89,T96,T01,T04 
        self.ps = geopack.recalc(self.dt_seconds)

        x=np.array([x])
        y=np.array([y])
        z=np.array([z])
        x=x.flatten()
        y=y.flatten()
        z=z.flatten() 
        nx=len(x)
        ny=len(y)
        nz=len(z)
        nn=min([nx,ny,nz])
        bx_out=np.zeros(nn,dtype=float)
        by_out=np.zeros(nn,dtype=float)
        bz_out=np.zeros(nn,dtype=float)
        for ix in range(nn):
            rr=np.sqrt(x[ix]**2+y[ix]**2+z[ix]**2)
            rr
            if (rr > 0.000001):
                bx_,by_,bz_=geopack.t89.t89(self.iopt,self.ps,x[ix],y[ix],z[ix])
                if self.use_igrf: bx0,by0,bz0=geopack.igrf_gsm(x[ix],y[ix],z[ix])
                else: bx0,by0,bz0=geopack.dip(x[ix],y[ix],z[ix])
                bx_out[ix]=bx_+bx0
                by_out[ix]=by_+by0
                bz_out[ix]=bz_+bz0
            else:
                bx_out[ix]=np.nan
                by_out[ix]=np.nan
                bz_out[ix]=np.nan
            
        return(np.column_stack((bx_out,by_out,bz_out)))            

# old code accepting only scalar x,y,z
#        bx_,by_,bz_=self.t89(self.iopt,self.ps,x,y,z)
#        if self.use_igrf: bx0,by0,bz0=geopack.igrf_gsm(x,y,z)
#        else: bx0,by0,bz0=geopack.dip(x,y,z)
#        return(bx_+bx0,by_+by0,bz_+bz0)


        
#    def t89(self,iopt,ps,x,y,z):
#        return geopack.t89.t89(iopt,ps,x,y,z)
        
    def trace(self,x,y,z,rlim=10.,r0=1.,dir=-1,maxloop=1000):
# returns the last x,y,z and arrays xx,yy,zz along trace
# x,y,z have to be scalar

# we need to call recalc since common block is shared between instances
# of geopack_2008 and T89,T96,T01,T04T04
        self.ps = geopack.recalc(self.dt_seconds)

        parmod=self.iopt
 #       geopack.trace(xi,yi,zi,dir,rlim=10,r0=1,parmod=2,exname='t89',inname='igrf',maxloop=1000)
        if self.use_igrf: return geopack.trace(x,y,z,dir,rlim,r0,parmod,'t89','igrf',maxloop=maxloop)
        else: return geopack.trace(x,y,z,dir,rlim,r0,parmod,'t89','dip',maxloop=maxloop)
        

#
# Initialization requires date and time information to initialize geopack (recalc) and obtain BY_GSM,BZ_GSM,Pressure and SYM_H values from 1-minute OMNI
#
class T96(Kamodo):
#
# using Sheng Tian's geopack module (https://github.com/tsssss/geopack)
#
    def __init__(self,year,month,day,hour,minute,use_igrf,*args,**kwargs):
        from geopack import t96
        super(T96, self).__init__(*args, **kwargs)
# epoch time since Jan. 1 1970 00:00 UT1
# datetime.timedelta
        dt=datetime.datetime(year,month,day,hour,minute)-datetime.datetime(1970,1,1)
# seconds from 1970/1/1 00:00 UT
        self.dt_seconds=dt.total_seconds()
        
        self.ps = geopack.recalc(self.dt_seconds)
        self.use_igrf=use_igrf

        from geospacepy import omnireader
        sTimeIMF = datetime.datetime(year,month,day,hour,minute)
        eTimeIMF = datetime.datetime(year,month,day,hour,minute)+datetime.timedelta(0,0,0,0,1,0)
        omniInt = omnireader.omni_interval(sTimeIMF,eTimeIMF,'1min')
        t = omniInt['Epoch'] #datetime timestamps
        By = omniInt['BY_GSM']
        Bz = omniInt['BZ_GSM']
        Pdyn = omniInt['Pressure']
        SYM_H = omniInt['SYM_H']

        self.parmod=np.array([Pdyn,SYM_H,By,Bz,0.,0.,0.,0.,0.,0.],dtype=float)

        bounds_error = kwargs.get('bounds_error', False)
        fill_value = kwargs.get('missing_value', np.nan)        

        units='nT'
        self.citation='Kamodo.T96 by Lutz Rastaetter (2020), Geopack/Tsyganenko by Sheng Tian (2019) and geospacepy-lite by Liam Kilkommons (2019)'

        self.x=np.linspace(-30.,10.,40) # make sure to avod (0,0,0)
        self.y=np.linspace(-10.,10.,20)
        self.z=np.linspace(-10.,10.,20)

        self.variables = dict(b_x = dict(units = 'nT', data = None),
                              b_y = dict(units = 'nT', data = None),
                              b_z = dict(units = 'nT', data = None),
                              bvec = dict(units = 'nT', data = None) )

        for varname in self.variables:
            units = self.variables[varname]['units']
            self.register_variable(varname, units)        

    def register_variable(self,varname,units):
        interpolator=None
        if varname == 'b_x':
            interpolator=self.bx
        if varname == 'b_y':
            interpolator=self.by
        if varname == 'b_z':
            interpolator=self.bz
        if varname == 'bvec':
            interpolator=self.b
            
        self.variables[varname]['interpolator']= interpolator

        def interpolate(xvec):
            return self.variables[varname]['interpolator'](xvec)

        interpolate.__doc__ = "A function that returns {} in [{}].".format(varname,units)

        self[varname] = kamodofy(interpolate, 
                                 units = units, 
                                 citation = self.citation,
                                 data = None)
        self[varname + '_ijk'] = kamodofy(gridify(self[varname], 
                                                  x_i = self.x, 
                                                  y_j = self.y, 
                                                  z_k = self.z),
                                          units = units,
                                          citation = self.citation,
                                          data = None) 
        

    def trace(self,x,y,z,rlim=10.,r0=1.,dir=-1,maxloop=1000):
# returns the last x,y,z and arrays xx,yy,zz along trace
        if self.use_igrf: return geopack.trace(x,y,z,dir,rlim,r0,self.parmod,'t96','igrf',maxloop=maxloop)
        else: return geopack.trace(x,y,z,dir,rlim,r0,self.parmod,'t96','dip',maxloop=maxloop)

    def bx(self,xvec):
        bx_,by_,bz_=self.b(xvec)
        return(bx_)

    def by(self,xvec):
        bx_,by_,bz_=self.b(xvec)
        return(by_)

    def bz(self,xvec):
        bx_,by_,bz_=self.b(xvec)
        return(bz_)

    def b(self,xvec):
# x,y,z can be an array or list        
        try:
            x,y,z=xvec
        except:
            x,y,z=xvec.T
# we need to call recalc since common block is shared between instances
# of geopack_2008 and T89,T96,T01,T04 
        self.ps = geopack.recalc(self.dt_seconds)
        x=np.array((x))
        y=np.array((y))
        z=np.array((z))
        x=x.flatten()
        y=y.flatten()
        z=z.flatten()
        nx=len(x)
        ny=len(y)
        nz=len(z)
        nn=min([nx,ny,nz])
        bx_out=np.zeros(nn,dtype=float)
        by_out=np.zeros(nn,dtype=float)
        bz_out=np.zeros(nn,dtype=float)
        for ix in range(nn):    
            rr=np.sqrt(x[ix]**2+y[ix]**2+z[ix]**2)
            if (rr > 0.000001):
                bx_,by_,bz_=geopack.t96.t96(self.parmod,self.ps,x[ix],y[ix],z[ix])
                if self.use_igrf: bx0,by0,bz0=geopack.igrf_gsm(x[ix],y[ix],z[ix])
                else: bx0,by0,bz0=geopack.dip(x[ix],y[ix],z[ix])
                bx_out[ix]=bx_+bx0
                by_out[ix]=by_+by0
                bz_out[ix]=bz_+bz0
            else:
                bx_out[ix]=np.nan
                by_out[ix]=np.nan
                bz_out[ix]=np.nan

        return(np.column_stack((bx_out,by_out,bz_out)))
    
#
# Initialization requires date and time information to initialize geopack (recalc) and obtain BY_GSM,BZ_GSM,Pressure and SYM_H values from 1-minute OMNI
#
class T01(Kamodo):
#
# using Sheng Tian's geopack module (https://github.com/tsssss/geopack)
#
    def __init__(self,year,month,day,hour,minute,use_igrf,*args,**kwargs):
        from geopack import t01
# epoch time since Jan. 1 1970 00:00 UT
# datetime.timedelta
        dt=datetime.datetime(year,month,day,hour,minute)-datetime.datetime(1970,1,1)
# seconds from 1970/1/1 00:00 UT
        self.dt_seconds=dt.total_seconds()
        
        self.ps = geopack.recalc(self.dt_seconds)
        self.use_igrf=use_igrf

        from geospacepy import omnireader
        sTimeIMF = datetime.datetime(year,month,day,hour,minute)
        eTimeIMF = datetime.datetime(year,month,day,hour,minute)+datetime.timedelta(0,0,0,0,1,0)
        omniInt = omnireader.omni_interval(sTimeIMF,eTimeIMF,'1min')
        t = omniInt['Epoch'] #datetime timestamps
        By = omniInt['BY_GSM']
        Bz = omniInt['BZ_GSM']
        Pdyn = omniInt['Pressure']
        SYM_H = omniInt['SYM_H']

        self.parmod=np.array([Pdyn,SYM_H,By,Bz,0.,0.,0.,0.,0.,0.],dtype=float)

        super(T01, self).__init__(*args, **kwargs)
#        parmod=np.zeros(10,dtype=float)
#        t89.tsyganenko.init_t89(year,month,day,hour,minute,use_igrf,0,parmod)
#        t89.tsyganenko.init_t89(int(year),int(month),int(day),int(hour),int(minute),int(use_igrf))
        bounds_error = kwargs.get('bounds_error', False)
        fill_value = kwargs.get('missing_value', np.nan)        

        self.citation='Kamodo.T01 by Lutz Rastaetter (2020), Geopack/Tsyganenko by Sheng Tian (2019) and geospacepy-lite by Liam Kilkommons (2019)'

        self.x=np.linspace(-30.,10.,40) # make sure to avoid (0,0,0)
        self.y=np.linspace(-10.,10.,20)
        self.z=np.linspace(-10.,10.,20)

        self.variables = dict(b_x = dict(units = 'nT', data = None),
                              b_y = dict(units = 'nT', data = None),
                              b_z = dict(units = 'nT', data = None),
                              bvec = dict(units = 'nT', data = None) )

        for varname in self.variables:
            units = self.variables[varname]['units']
            self.register_variable(varname, units)        

    def register_variable(self,varname,units):
        if varname == 'b_x':
            interpolator=self.bx
        if varname == 'b_y':
            interpolator=self.by
        if varname == 'b_z':
            interpolator=self.bz
        if varname == 'bvec':
            interpolator=self.b
            
        self.variables[varname]['interpolator']= interpolator

        def interpolate(xvec):
            return self[varname]['interpolator'](xvec)

        interpolate.__doc__ = "A function that returns {} in [{}].".format(varname,units)

        self[varname] = kamodofy(interpolate, 
                                 units = units, 
                                 citation = self.citation,
                                 data = None)
        self[varname + '_ijk'] = kamodofy(gridify(self[varname], 
                                                  x_i = self.x, 
                                                  y_j = self.y, 
                                                  z_k = self.z),
                                          units = units,
                                          citation = self.citation,
                                          data = None) 
        
    def trace(self,x,y,z,rlim=10.,r0=1.,dir=-1,maxloop=1000):
# returns the last x,y,z and arrays xx,yy,zz along trace
        if self.use_igrf: return geopack.trace(x,y,z,dir,rlim,r0,self.parmod,'t01','igrf',maxloop=maxloop)
        else: return geopack.trace(x,y,z,dir,rlim,r0,self.parmod,'t01','dip',maxloop=maxloop)
        
    def bx(self,xvec):
        bx_,by_,bz_=self.b(xvec)
        return(bx_)

    def by(self,xvec):
        bx_,by_,bz_=self.b(xvec)
        return(by_)

    def bz(self,xvec):
        bx_,by_,bz_=self.b(xvec)
        return(bz_)

    def b(self,xvec):
        try:
            x,y,z=xvec
        except:
            x,y,z=xvec.T
# x,y,z can be an array or list        
# we need to call recalc since common block may be shared between instances of T89
        self.ps = geopack.recalc(self.dt_seconds)

        x=np.array([x])
        y=np.array([y])
        z=np.array([z])
        x=x.flatten()
        y=y.flatten()
        z=z.flatten()
        nx=len(x)
        ny=len(y)
        nz=len(z)
        nn=min([nx,ny,nz])
        bx_out=np.zeros(nn,dtype=float)
        by_out=np.zeros(nn,dtype=float)
        bz_out=np.zeros(nn,dtype=float)

        for ix in range(nn): 
            rr=sqrt(x[ix]**2+y[ix]**2+z[ix]**2)
            if (rr > 0.000001):
                bx_,by_,bz_=geopack.t01.t01(self.parmod,self.ps,x[ix],y[ix],z[ix])
                if self.use_igrf: bx0,by0,bz0=geopack.igrf_gsm(x[ix],y[ix],z[ix])
                else: bx0,by0,bz0=geopack.dip(x[ix],y[ix],z[ix])
                bx_out[ix]=bx_+bx0
                by_out[ix]=by_+by0
                bz_out[ix]=bz_+bz0
            else:
                bx_out[ix]=np.nan
                by_out[ix]=np.nan
                bz_out[ix]=np.nan

        return(np.column_stack((bx_out,by_out,bz_out)))
    
#
# Initialization requires date and time information to initialize geopack (recalc) and obtain BY_GSM,BZ_GSM,Pressure and SYM_H values from 1-minute OMNI
#
class T04(Kamodo):
#
# using Sheng Tian's geopack module (https://github.com/tsssss/geopack)
#
    def __init__(self,year,month,day,hour,minute,use_igrf,QD_data_path='%s/Kamodo_data/Qin-Denton/' % os.environ['HOME'],*args,**kwargs):
        from geopack import t04
# epoch time since Jan. 1 1970 00:00 UT
# datetime.timedelta
        dt=datetime.datetime(year,month,day,hour,minute)-datetime.datetime(1970,1,1)
# seconds from 1970/1/1 00:00 UT
        self.dt_seconds=dt.total_seconds()
        self.use_igrf=use_igrf
        qin_denton_url_message="files for years starting with 1995 are available - CCMC are truing to update files monthly as OMNI 1-minute, 5-minute and hourly data bcome available"
# now powered by iSWA at CCMC
        qin_denton_url='https://iswa.gsfc.nasa.gov/iswa_data_tree/composite/magnetosphere/Qin-Denton/1min/%d/%.02d/' %(year,month)
#        qin_denton_url='https://rbsp-ect.newmexicoconsortium.org/data_pub/QinDenton/%d/' % (year)

        qin_denton_file='QinDenton_%d%.02d%.02d_1min.txt' % (year,month,day)
        if (QD_data_path):
            qin_denton_local_path='%s/Qin_denton_daily/%d/' % (QD_data_path,year)

# Qin-Denton file name 
        qin_denton_local_file=qin_denton_local_path+'%s' % (qin_denton_file)
        qin_denton_file_url=qin_denton_url+'/'+qin_denton_file
#        print(qin_denton_local_file)
# check local Qin-Denton file and download as needed
        import pandas as pd
        import requests
        import dateutil.parser
# create local data directory and file if necessary
        if (not os.path.isdir(qin_denton_local_path)):
            os.makedirs(qin_denton_local_path,exist_ok=True)

# were are not only testing for the presence of a file but also its size when the file does exist
# a download may have failed before and may leave an incomplete file or a file containing a HTML META tag with a redirect pointing to the iSWA splash page (https://ccmc.gsfc.nasa.gov/iswa/)
        download_qd_file=False
        if (not os.path.isfile(qin_denton_local_file)):
            download_qd_file=True
        else:
            if (os.path.getsize(qin_denton_local_file) < 347621):
                download_qd_file=True
            else:
                local_file_modified=(os.stat(qin_denton_local_file))[8]
# check last modified time of remote file against local file
                response=requests.head(qin_denton_file_url)
                last_modified=response.headers.get('Last-Modified')
                if last_modified:
                    last_modified = dateutil.parser.parse(last_modified).timestamp()
                    if last_modified > local_file_modified:
                        download_qd_file=True
                        
        if (download_qd_file):
            # download file
            print('Obtaining Qin Denton file from %s' % qin_denton_file_url)

            r=requests.get(qin_denton_file_url)
            myfiledata=r.text
            if len(myfiledata) > 1:
                f = open(qin_denton_local_file, "w")
                f.write(myfiledata)
                f.close()
            else:
                raise ValueError("no data found - the model cannot run for this day")
                    
# read daily file into frame
        print('reading Qin Denton file %s' % qin_denton_local_file)
        qindenton_frame=pd.read_table(qin_denton_local_file,skiprows=192,sep='\s+',names=['date','year','month','day','hour','minute','second','ByIMF','BzIMF','Vsw','Den_P','Pdyn','G1','G2','G3','ByIMF_status','BzIMF_status','Vsw_status','Den_P_status','Pdyn_status','G1_ststus','G2_status','G3_status','Kp','akp3','Dst','Bz','W1','W2','W3','W4','W5','W6','W1_status', 'W2_status', 'W3_status', 'W4_status', 'W5_status', 'W6_status'])


#        from geospacepy import omnireader
        sTimeIMF = datetime.datetime(year,month,day,hour,minute)
        eTimeIMF = datetime.datetime(year,month,day,hour,minute)+datetime.timedelta(0,0,0,0,1,0)
#        omniInt = omnireader.omni_interval(sTimeIMF,eTimeIMF,'1min')
        time_day = qindenton_frame['date'] #datetime timestamps
        By_day = qindenton_frame['ByIMF']
        Bz_day = qindenton_frame['BzIMF']
        Pdyn_day = qindenton_frame['Pdyn']
        Dst_day = qindenton_frame['Dst']
        w1_day=qindenton_frame['W1']
        w2_day=qindenton_frame['W2']
        w3_day=qindenton_frame['W3']
        w4_day=qindenton_frame['W4']
        w5_day=qindenton_frame['W5']
        w6_day=qindenton_frame['W6']
        time_index=hour*60+minute
# need Qin-Denton parameters
        w1=w1_day[time_index]
        w2=w2_day[time_index]
        w3=w3_day[time_index]
        w4=w4_day[time_index]
        w5=w5_day[time_index]
        w6=w6_day[time_index]
        Pdyn=Pdyn_day[time_index]
        Dst=Dst_day[time_index]
        By=By_day[time_index]
        Bz=Bz_day[time_index]
# end Qin-Dention acquisition        

        self.parmod=np.array([Pdyn,Dst,By,Bz,w1,w2,w3,w4,w5,w6],dtype=float)

        super(T04, self).__init__(*args, **kwargs)
#        parmod=np.zeros(10,dtype=float)
#        t89.tsyganenko.init_t89(year,month,day,hour,minute,use_igrf,0,parmod)
#        t89.tsyganenko.init_t89(int(year),int(month),int(day),int(hour),int(minute),int(use_igrf))
        bounds_error = kwargs.get('bounds_error', False)
        fill_value = kwargs.get('missing_value', np.nan)        

        self.units='nT'
        self.citation='Kamodo.T89 by Lutz Rastaetter (2020), Geopack/Tsyganenko by Sheng Tian (2019) and geospacepy-lite by Liam Kilkommons (2019)'

        self.x=np.linspace(-30.,10.,40) # make sure to avoid (0,0,0)
        self.y=np.linspace(-10.,10.,20)
        self.z=np.linspace(-10.,10.,20)

        self.variables = dict(b_x = dict(units = 'nT', data = None),
                              b_y = dict(units = 'nT', data = None),
                              b_z = dict(units = 'nT', data = None),
                              bvec = dict(units = 'nT', data = None))
                     
        for varname in self.variables:
            units = self.variables[varname]['units']
            self.register_variable(varname, units)        
               
                     
    def register_variable(self,varname,units):
        if varname == 'b_x':
            interpolator=self.bx
        if varname == 'b_y':
            interpolator=self.by
        if varname == 'b_z':
            interpolator=self.bz
        if varname == 'bvec':
            interpolator=self.b
            
        self.variables[varname]['interpolator']= interpolator

        def interpolate(xvec):
            return self[varname]['interpolator'](xvec)

        interpolate.__doc__ = "A function that returns {} in [{}].".format(varname,units)

        self[varname] = kamodofy(interpolate, 
                                 units = units, 
                                 citation = self.citation,
                                 data = None)
        self[varname + '_ijk'] = kamodofy(gridify(self[varname], 
                                                  x_i = self.x, 
                                                  y_j = self.y, 
                                                  z_k = self.z),
                                          units = units,
                                          citation = self.citation,
                                          data = None) 

    def register_variable_2(self,varname,units):
        interpolator=None;
        if varname == 'b_x':
            interpolator=self.bx
        if varname == 'b_y':
            interpolator=self.by
        if varname == 'b_z':
            interpolator=self.bz
        if varname == 'bvec':
            interpolator=self.b
            
        self.variables[varname]['interpolator']= interpolator

        def interpolate(xvec):
            return self[varname]['interpolator'](xvec)

        
        interpolate.__doc__ = "A function that returns {} in [{}].".format(varname,units)

        self[varname] = kamodofy(interpolate, 
                                 units = units, 
                                 citation = self.citation,
                                 data = None)
        self[varname + '_ijk'] = kamodofy(gridify(self[varname], 
                                                  x_i = self.x, 
                                                  y_j = self.y, 
                                                  z_k = self.z),
                                          units = units,
                                          citation = self.citation,
                                          data = None)            

    def trace(self,x,y,z,rlim=10.,r0=1.,dir=-1,maxloop=1000):
# returns the last x,y,z and arrays xx,yy,zz along trace
        if self.use_igrf: return geopack.trace(x,y,z,dir,rlim,r0,self.parmod,'t04','igrf',maxloop=maxloop)
        else: return geopack.trace(x,y,z,dir,rlim,r0,self.parmod,'t01','dip',maxloop=maxloop)

    def bx(self,xvec):
        bx_,by_,bz_=self.b(xvec)
        return(bx_)

    def by(self,xvec):
        bx_,by_,bz_=self.b(xvec)
        return(by_)

    def bz(self,xvec):
        bx_,by_,bz_=self.b(xvec)
        return(bz_)

    def b(self,xvec):
# x,y,z can be an array or list        
        try:
            x,y,z=xvec
        except:
            x,y,z=xvec.T
# we need to call recalc since common block may be shared between instances
# of geopack-2008 and T89,T96,T01,T04
        self.ps = geopack.recalc(self.dt_seconds)

        x=np.array([x])
        y=np.array([y])
        z=np.array([z])
        x=x.flatten()
        y=y.flatten()
        z=z.flatten()
        nx=len(x)
        ny=len(y)
        nz=len(z)
        nn=min([nx,ny,nz])
        bx_out=np.zeros(nn,dtype=float)
        by_out=np.zeros(nn,dtype=float)
        bz_out=np.zeros(nn,dtype=float)

        for ix in range(nn):
            rr=np.sqrt(x[ix]**2+y[ix]**2+z[ix]**2)
            if (rr > 0.000001):
                bx_,by_,bz_=geopack.t04.t04(self.parmod,self.ps,x[ix],y[ix],z[ix])
                if self.use_igrf: bx0,by0,bz0=geopack.igrf_gsm(x[ix],y[ix],z[ix])
                else: bx0,by0,bz0=geopack.dip(x[ix],y[ix],z[ix])
                bx_out[ix]=bx_+bx0
                by_out[ix]=by_+by0
                bz_out[ix]=bz_+bz0
            else:
                bx_out[ix]=np.nan
                by_out[ix]=np.nan
                bz_out[ix]=np.nan
                         
        return(np.column_stack((bx_out,by_out,bz_out)))



