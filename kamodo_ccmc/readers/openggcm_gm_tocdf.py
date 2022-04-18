#import numpy as np
from numpy import transpose,zeros,array,append,float32,NAN,vectorize
from time import perf_counter
from netCDF4 import Dataset
#from astropy.constants import R_earth

from datetime import datetime,timezone,timedelta
from numpy import zeros, flip,array,float32,NaN

from glob import glob
#import kamodo_ccmc.readers.OpenGGCM.read_b_grids as rbg
#import kamodo_ccmc.readers.OpenGGCM.readmagfile3d as rmhd
import kamodo_ccmc.readers.OpenGGCM.readOpenGGCM as ropgm

from os.path import sep,isfile,isdir,dirname,exists
from os import remove

# this function takes timestamps with fields separated by ":" and seconds being a float including milliseconds
#@vectorize
def dts_to_hrs(datetime_string, filedate):
    '''Get hours since midnight from datetime string'''
    
    return (datetime.strptime(datetime_string, '%Y:%m:%d:%H:%M:%S.%f').replace(tzinfo=timezone.utc)\
            -filedate).total_seconds()/3600.
           
openggcm_gm_variable_to_GSE_factors={"bx":-1,"by":-1,"bz":1,
                             "bx1":-1,"by1":-1,"bz1":1,
                             "vx":-1,"vy":-1,"vz":1,
                             "ex":-1,"ey":-1,"ez":1,
                             "xjx":-1,"xjy":-1,"xjz":1,
                             "rr":1,"pp":1,"resis":1}
openggcm_gm_variable_units={"B_x":'nT',"B_y":'nT',"B_z":'nT',
                    "B1_x":'nT',"B1_y":'nT',"B1_z":'nT',
                    "E_xedges":'mV/m',"E_yedges":'mV/m',"E_zedges":'mV/m',
                    "V_x":'km/s',"V_y":'km/s',"V_z":'km/s',
                    "J_x":'muA/m^2',"J_y":'muA/m^2',"J_z":'muA/m^2',
                    "rho":'kg/m^3',"P":'pPa',"eta":'m^2/s'}
       
openggcm_gm_varnames={
    ### 3D spatial variables - to be aggregated into 4D ###
# B field on cell centers - used in online visualization
    'bx':['B_x','x component of magnetic field',0,'GSE','car',['time','x','y','z'],'nT'],
    'by':['B_y','y component of magnetic field',1,'GSE','car',['time','x','y','z'],'nT'],
    'bz':['B_z','z component of magnetic field',2,'GSE','car',['time','x','y','z'],'nT'],
    # these are B1_x, B1_y, B1_z used for file line tracing in online visualization/Kameleon
    'bx1':['B_xfaces','x component of magnetic field (on grid cell faces)',0,'GSE','car',['time','x_bxfacces','y_bxfacces','z_bxfacces'],'nT'],
    'by1':['B_yfaces','y component of magnetic field (on grid cell faces)',1,'GSE','car',['time','x_byfacces','y_byfacces','z_byfacces'],'nT'],
    'bz1':['B_zfaces','z component of magnetic field (on grid cell faces)',2,'GSE','car',['time','x_bzfacces','y_bzfacces','z_bzfacces'],'nT'],
    'ex':['E_xedges','x component of electric field (on grid cell edges)',3,'GSE','car',['time','x_exedges','y_exedges','z_exedges'],'mV/m'],
    'ey':['E_yedges','y component of electric field (on grid cell edges)',4,'GSE','car',['time','x_eyedges','y_eyedges','z_eyedges'],'mV/m'],
    'ez':['E_zedges','z component of electric field (on grid cell edges)',5,'GSE','car',['time','x_ezedges','y_ezedges','z_ezedges'],'mV/m'],
    'vx':['V_x','x component of plasma velocity',6,'GSE','car',['time','x','y','z'],'km/s'],
    'vy':['V_y','y component of plasma velocity',7,'GSE','car',['time','x','y','z'],'km/s'],
    'vz':['V_z','z component of plasma velocity',8,'GSE','car',['time','x','y','z'],'km/s'],
    'xjx':['J_x','x component of current density',9,'GSE','car',['time','x','y','z'],'muA/m**2'],
    'xjy':['J_y','y component of current density',10,'GSE','car',['time','x','y','z'],'muA/m**2'],
    'xjz':['J_z','z component of current density',11,'GSE','car',['time','x','y','z'],'muA/m**2'],
    'rr':['N','plasma number denstity (hydrogen equivalent)',12,'GSE','car',['time','x','y','z'],'1/cm**3'],
    'resis':['eta','resistivity',13,'GSE','car',['time','x','y','z'],'m**2/s'],
    'pp':['P','plasma pressure',14,'GSE','car',['time','x','y','z'],'pPa'],
    # 2D spatial variables aggregated into 3D
}


def openggcm_combine_magnetosphere_files(full_file_prefix,cadence=None,requested_variables=None,verbose=False):
    # file prefix includes everything including '3df' for magnetosphere and intended year,month,day (and hour)
    # read matching dates, times from 3df_list to generate list of raw files.
    # coadence: Default: None -- use all available times, otherwise use cadence as input starting with first time in simualtion.
    # requested_variables default: None -- convert all available 3D fields, other vise use only the requeted variables to generate smaller NetCDF files.
    
    tic=perf_counter()
    file_prefix,file_datetime = full_file_prefix.split('.3df')
    file_path=dirname(file_prefix)+sep

    
    if verbose:
        print('file_prefix: ',file_prefix)
        print('file_datetime: ',file_datetime)

    file_prefix=file_prefix+'.3df'
    
    file_date_time=file_datetime.split('_')[1:]
    file_year,file_mon,file_day=file_date_time[0].split("-")
    file_date_fmt=file_year+'/'+file_mon+'/'+file_day
    if verbose:
        print(file_path)

    # look for runme input file in the current or the parent directory
    runme_file=file_path+'/runme'
    if not exists(runme_file):
        runme_file=file_path+'/../runme'
        if not exists(runme_file):
            print("openggcm_to_cdf: required 'runme' file not found!")
            return False

    runme_file_object=open(runme_file,'r')
    runme_inputs=runme_file_object.readlines()
    runme_file_object.close()
        
    list_file=file_prefix+'_list'
    if not isfile(list_file):
        print("Required list file not found: ",list_file)
        return False
    
    f=open(list_file,'r')
    list_data=f.readlines()
    f.close()
    if verbose:
        print('file_date_time',file_date_time)
        print('len(file_date_time): ',len(file_date_time))

    added_time_at_beginning=0
    added_time_at_end=0

    
    reduced_list_data=[list_data[1]]
    if cadence is not None:
        file_seconds0=int(((list_data[1].split())[0])[-6:])
        print(file_seconds0,cadence)
        for file_str in list_data[2:]:
            file_seconds=int(((file_str.split())[0])[-6:])
            if (file_seconds-file_seconds0) >= cadence:
                #file_seconds0=file_seconds
                file_seconds0=file_seconds0+cadence
                reduced_list_data.append(file_str)
                print('appending '+file_str)
    else:
        reduced_list_data=list_data[1:]
    
    if len(file_date_time) == 1:
        file_str_list=[s for s in reduced_list_data if file_date_fmt in s]
        if verbose:
            print(file_str_list)
        if cadence is not None:
            file_seconds0=int(file_str_list[0][-6:])
            print(file_seconds0)
            file_str_list2=[file_str_list[0]]
            for file_str in file_str_list:
                file_seconds=int(file_str[-6:])
                if file_seconds >= cadence:
                    file_seconds0=file_seconds
#                    file_seconds0=file_seconds0+cadence # 
                    file_str_list2.append(file_str)
            file_str_list=file_str_list2
                
        # add first file for next day
        this_datetime=datetime.strptime(file_date_time[0]+' 00:00:00', 
                                              '%Y-%m-%d %H:%M:%S').replace(tzinfo=timezone.utc)
        next_datetime=this_datetime.add+timedelta(days=1)
        next_date_fmt=next_datetime.strftime('%Y/%m/%d')
        next_file_str_list=[s for s in reduced_list_data if next_date_fmt in s]
        if next_file_str_list is not []:
            file_str_list=file_str_list.append(next_file_str_list[0])
            added_time_at_end=1

    if len(file_date_time) > 1:
        file_hour=file_date_time[1]
        file_str_list=[s for s in reduced_list_data if file_date_fmt in s and " "+file_hour+":" in s]

        this_datetime=datetime.strptime(file_date_time[0]+' '+file_hour+':00:00',
                                        '%Y-%m-%d %H:%M:%S').replace(tzinfo=timezone.utc)
#
# wrap in time: add first file from next hour or day
#
        next_datetime=this_datetime+timedelta(hours=1)
        next_date_fmt=next_datetime.strftime('%Y/%m/%d')
        next_hour=next_datetime.strftime('%H')
        next_file_str_list=[s for s in reduced_list_data if next_date_fmt in s and " "+next_hour+":" in s]
        if len(next_file_str_list) > 0:
            file_str_list.append(next_file_str_list[0])
            added_time_at_end=1
            

    file_names_raw=[file_str.split(' ')[0] for file_str in file_str_list]
    filedate=datetime.strptime(file_date_fmt,'%Y/%m/%d').replace(tzinfo=timezone.utc)
    
    file_names=[]
    for ifile in range(len(file_names_raw)):
        filename=file_names_raw[ifile]
        if exists(file_path+filename):
            file_names.append(filename)

    
    # find grid files
    filename=file_names[0]
    file_seconds=filename[-6:]

    # look for two possible grid files in the current or the parent directory
    grid_file_name=filename[:-11]+".grid"
    grid_file_b_name=filename[:-11]+".grid"
    grid_file_b=None
    grid_file=file_path+'/'+grid_file_name
    xgrid_file=file_path+sep+'gridx.txt'
    ygrid_file=file_path+sep+'gridy.txt'
    zgrid_file=file_path+sep+'gridz.txt'
    if exists(xgrid_file) and exists(ygrid_file) and exists(zgrid_file):
        grid_file=[xgrid_file,ygrid_file,zgrid_file]
    else:            
        if not exists(grid_file):
            grid_file=file_path+'/../'+grid_file_name
            if not exists(grid_file):
                grid_file_=None
            grid_file_b=file_path+'/'+grid_file_b_name
        if grid_file is None:
            if not exists(grid_file_b):
                grid_file_b=file_path+'/../'+grid_file_b_name
                if exists(grid_file_b):
                    grid_file=grid_file_b

    if verbose:
        print('grid file with grid for b and e: ',grid_file_b)
    openggcm_missing=1.e30
            
    if grid_file is not None:
        if len(grid_file) > 1:
            # separate text files containing (x,dx,x_bx), (y,dy,y_by) or (z,dz,z_bz)
            # (x,y,z) in the model are (-x_GSE,-y_GSE,z_GSE)
            # using flip() to make arrays with increasing values in x and y
            gridx_file_object=open(xgrid_file)
            gridx_data=gridx_file_object.readlines()
            gridx_file_object.close()

            gx=-flip(array([float((line.split())[0]) for line in gridx_data[1:]]))
            nx=len(gx)

            gx_bx=-flip(array([float((line.split())[2]) for line in gridx_data[1:]]))
            gx_by=gx
            gx_bz=gx

            gx_ex=gx
            gx_ey=gx_bx
            gx_ez=gx_bx

            gridy_file_object=open(ygrid_file)
            gridy_data=gridy_file_object.readlines()
            gridy_file_object.close()
            
            gy=-flip(array([float((line.split())[0]) for line in gridy_data[1:]]))
            ny=len(gy)

            gy_bx=gy
            gy_by=-flip(array([float((line.split())[2]) for line in gridy_data[1:]]))
            gy_bz=gy

            gy_ex=gy_by
            gy_ey=gy
            gy_ez=gy_by


            gridz_file_object=open(zgrid_file)
            gridz_data=gridz_file_object.readlines()
            gridz_file_object.close()

            gz=array([float((line.split())[0]) for line in gridz_data[1:]])
            nz=len(gz)

            gz_bx=gz
            gz_by=gz
            gz_bz=array([float((line.split())[2]) for line in gridz_data[1:]])

            gz_ex=gz_bz
            gz_ey=gz_bz
            gz_ez=gz
            
        else:
                
            gx=zeros(5000)
            gy=zeros(5000)
            gz=zeros(5000)
            nx,ny,nz,gx,gy,gz = ropgm.read_grid_for_vector(grid_file,' ',gx,gy,gz)
            if verbose:
                print('NX: %i NY: %i N: %i' % (nx,ny,nz) )
            
                gx=-flip(gx[0:nx])
                gy=-flip(gy[0:ny])
                gz=gz[0:nz]
                gx_bx=zeros(nx)
                gy_bx=zeros(ny)
                gz_bx=zeros(nz)
                gx_by=zeros(nx)
                gy_by=zeros(ny)
                gz_by=zeros(nz)
                gx_bz=zeros(nx)
                gy_bz=zeros(ny)
                gz_bz=zeros(nz)
                gx_ex=zeros(nx)
                gy_ex=zeros(ny)
                gz_ex=zeros(nz)
                gx_ey=zeros(nx)
                gy_ey=zeros(ny)
                gz_ey=zeros(nz)
                gx_ez=zeros(nx)
                gy_ez=zeros(ny)
                gz_ez=zeros(nz)
                
                nx1,ny1,nz1,gx_bx,gy_bx,gz_bx = ropgm.read_grid_for_vector(grid_file,'bx',gx_bx,gy_bx,gz_bx)
                if nx1 <= 0 or ny1 <=0 or nz1 <=0:
                    raise IOError("bx grid not found")
                gx_bx=-flip(gx_bx[0:nx1])
                gy_bx=-flip(gy_bx[0:ny1])
                gz_bx=gz_bx[0:nz1]

                nx1,ny1,nz1,gx_by,gy_by,gz_by = ropgm.read_grid_for_vector(10,grid_file,'by',gx_by,gy_by,gz_by)
                if nx1 <= 0 or ny1 <=0 or nz1 <=0:
                    raise IOError("by grid not found")
                gx_by=-flip(gx_by[0:nx1])
                gy_by=-flip(gy_by[0:ny1])
                gz_by=gz_by[0:nz1]

                nx1,ny1,nz1,gx_bz,gy_bz,gz_bz = ropgm.read_grid_for_vector(10,grid_file,'bz',gx_bz,gy_bz,gz_bz)
                if nx1 <= 0 or ny1 <=0 or nz1 <=0:
                    raise IOError("bz grid not found")

                gx_bz=-flip(gx_bz[0:nx1])
                gy_bz=-flip(gy_bz[0:ny1])
                gz_bz=gz_bz[0:nz1]

                nx1,ny1,nz1,gx_ex,gy_ex,gz_ex = ropgm.read_grid_for_vector(grid_file,'ex',gx_ex,gy_ex,gz_ex)
                if nx1 <= 0 or ny1 <=0 or nz1 <=0:
                    raise IOError("ex grid not found")
                gx_ex=-flip(gx_ex[0:nx1])
                gy_ex=-flip(gy_ex[0:ny1])
                gz_ex=gz_ex[0:nz1]

                nx1,ny1,nz1,gx_ey,gy_ey,gz_ey = ropgm.read_grid_for_vector(10,grid_file,'ey',gx_ey,gy_ey,gz_ey)
                if nx1 <= 0 or ny1 <=0 or nz1 <=0:
                    raise IOError("ey grid not found")
                gx_ey=-flip(gx_ey[0:nx1])
                gy_ey=-flip(gy_ey[0:ny1])
                gz_ey=gz_ey[0:nz1]

                nx1,ny1,nz1,gx_ez,gy_ez,gz_ez = ropgm.read_grid_for_vector(10,grid_file,'ez',gx_ez,gy_ez,gz_ez)
                if nx1 <= 0 or ny1 <=0 or nz1 <=0:
                    raise IOError("ez grid not found")
                gx_ez=-flip(gx_ez[0:nx1])
                gy_ez=-flip(gy_ez[0:ny1])
                gz_ez=gz_ez[0:nz1]
            else:
                raise IOError("grid file not found")


    #initialize single output file
    nc_file=full_file_prefix+".nc"
    if exists(nc_file):
        remove(nc_file)
    else:
        print("The file %s does not exist" % nc_file) 
    
    # NetCDF root group
    data_out = Dataset(nc_file, 'w', format='NETCDF4')
    data_out.model = 'OpenGGCM_GM'
    data_out.file=''.join([f+',' for f in file_names]).strip(',')  #csv list of files
#    x_dim=data_out.create_dimension("x",len(gx))
#    y_dim=data_out.create_dimension("y",len(gy))
#    z_dim=data_out.create_dimension("z",len(gz))
#    time_dim=data_out.create_dimension("time",ntime)
    
    # store dimensions
    dim_dict={}
#    dim_dict['time']=array(file_times,dtype=float32)  # {'data':file_times,'datatype':float32,'size':len(file_times)}

    dim_dict['x']=array(gx,dtype=float32) # {'data':gx,'datatype':float32,'size':len(gx)}
    dim_dict['y']=array(gy,dtype=float32) # {'data':gy,'datatype':float32,'size':len(gy)}
    dim_dict['z']=array(gz,dtype=float32) # {'data':gz,'datatype':float32,'size':len(gz)}
    dim_dict['x_bxfaces']=gx_bx # {'data':gx_bx,'datatype':float32,'size':len(gx_bx)}
    dim_dict['y_bxfaces']=gy_bx # {'data':gy_bx,'datatype':float32,'size':len(gy_bx)}
    dim_dict['z_bxfaces']=gz_bx # {'data':gz_bx,'datatype':float32,'size':len(gz_bx)}
    dim_dict['x_byfaces']=gx_by # {'data':gx_by,'datatype':float32,'size':len(gx_by)}
    dim_dict['y_byfaces']=gy_by # {'data':gy_by,'datatype':float32,'size':len(gy_by)}
    dim_dict['z_byfaces']=gz_by # {'data':gz_by,'datatype':float32,'size':len(gz_by)}
    dim_dict['x_bzfaces']=gx_bz # {'data':gx_bz,'datatype':float32,'size':len(gx_bz)}
    dim_dict['y_bzfaces']=gy_bz # {'data':gy_bz,'datatype':float32,'size':len(gy_bz)}
    dim_dict['z_bzfaces']=gz_bz # {'data':gz_bz,'datatype':float32,'size':len(gz_bz)}
    dim_dict['x_exedges']=gx_ex # {'data':gx_ex,'datatype':float32,'size':len(gx_ex)}
    dim_dict['y_exedges']=gy_ey # {'data':gy_ex,'datatype':float32,'size':len(gy_ex)}
    dim_dict['z_exedges']=gz_ez # {'data':gz_ex,'datatype':float32,'size':len(gz_ex)}
    dim_dict['x_eyedges']=gx_ex # {'data':gx_ey,'datatype':float32,'size':len(gx_ey)}
    dim_dict['y_eyedges']=gy_ey # {'data':gy_ey,'datatype':float32,'size':len(gy_ey)}
    dim_dict['z_eyedges']=gz_ez # {'data':gz_ey,'datatype':float32,'size':len(gz_ey)}
    dim_dict['x_ezedges']=gx_ex # {'data':gx_ez,'datatype':float32,'size':len(gx_ez)}
    dim_dict['y_ezedges']=gy_ey # {'data':gy_ez,'datatype':float32,'size':len(gy_ez)}
    dim_dict['z_ezedges']=gz_ez # {'data':gz_ez,'datatype':float32,'size':len(gz_ez)}

    
#    return(nc_file,dim_dict)

    # the dimension variables do not include time (in hours) yet
    for dim in dim_dict.keys():
        new_dim = data_out.createDimension(dim, dim_dict[dim].size)
        new_var = data_out.createVariable('_'+dim, float32, dim) # dim_dict[dim]['datatype'],dim)
        new_var[:] = dim_dict[dim]
        if dim=='time': units = 'hr'
        else: units = 'R_E'
        new_var.units = units
        
    if verbose: print('Dimensions complete.\n')
#    data_out.close() 
#    return(nc_file,dim_dict)
   
    ntime=len(file_names)
    # variable loop
    #    print(nx,ny,nz,ntime)

    data=zeros(shape=(ntime,nx,ny,nz),dtype=float32,order='C')
#    fielddata=zeros(shape=(nx,ny,nz),dtype=float32,order='F')
    missing_value=NaN
    time=zeros(ntime)
    if requested_variables is None:
        convert_varnames=openggcm_gm_varnames
    else:
        convert_varnames=[]
        for varname in openggcm_gm_varnames:
            print(varname)
            if varname in requested_variables or openggcm_gm_varnames[varname][0] in requested_variables:
                convert_varnames.append(varname)                
        
#    for varname in ['bx']: #openggcm_gm_varnames:  #model_variables:
#    for varname in openggcm_gm_varnames:  #model_variables:
    for varname in convert_varnames:  #model_variables:
        variable_found=True
        data[:]=0.
        for ifile in range(ntime):
            filename=file_names[ifile]
            fielddata=zeros(shape=(nx,ny,nz),dtype=float32,order='F')
            fieldarray,nx_field,ny_field,nz_field,asciitime = ropgm.read_3d_field(file_path+filename,fielddata,varname);
            asciitime_str=asciitime.decode('utf-8')
            if verbose:
                print('asciitime_str: ',asciitime_str)
                print('time: ',time)
            try:
                timestr,time_seconds,time_epoch,time_str_colon,datatype = asciitime_str.split()
                if verbose:
                    print(time_str_colon)
                time[ifile]=dts_to_hrs(time_str_colon, filedate)
            except:
                if verbose:
                    print('empty asciitime string')
                variable_found=False
                break
            
            factor=openggcm_gm_variable_to_GSE_factors[varname]
#            fieldarray[(fieldarray == openggcm_missing)]=missing_value
            data[ifile,:,:,:]=factor*flip(flip(fieldarray,axis=0),axis=1)
        print(varname,data.min(),data.max(),variable_found)

        if variable_found:
        # add time in hours to dimensions
            dim='time'
            try:
                new_dim = data_out.createDimension(dim, len(time)) # should error if already there
                new_var = data_out.createVariable('_'+dim, float32, dim) # dim_dict[dim]['datatype'],dim)
                new_var[:] = time
                new_var.units = 'hr'
            except:
                if verbose:
                    print("time already in NetCDF file dimensions")

            dims = openggcm_gm_varnames[varname][-2]

            new_var = data_out.createVariable(varname, float32, dims)
            new_var[:] = data 
            new_var.units=openggcm_gm_varnames[varname][-1]

    runme_inputs=array(runme_inputs)
    runme_has_isphere=[line.find('ISPHERE') == 0 for line in runme_inputs]
    runme_isphere_line=(runme_inputs[runme_has_isphere])[0]
    isphere_name_value=runme_isphere_line.split()

    data_out.near_Earth_boundary_radius=float(isphere_name_value[1])
    data_out.near_Earth_boundary_radius_units='R_E'

    data_out.added_time_at_beginning=added_time_at_beginning
    data_out.added_time_at_end=added_time_at_end

    data_out.filedate=filedate.strftime("%Y-%m-%d")

    data_out.runme=runme_inputs
    print(f"Data for {full_file_prefix} converted in {perf_counter()-tic:.6f}s.")
    data_out.close() 

    toc=perf_counter()
    print('converted files in ',toc-tic,' seconds')                       
    return True

