#import numpy as np
from numpy import transpose,reshape,array,zeros,append,float32,NAN,vectorize,arange
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
           
# standard model dictionary for reference
openggcm_ie_varnames = {'sigh': ['Sigma_H', 'Hall conductance',
                         0, 'SM', 'sph', ['time', 'lon', 'lat'], '1/Ohm'],
                  'sigp': ['Sigma_P', 'Pedersen conductance',
                         1, 'SM', 'sph', ['time', 'lon', 'lat'], '1/Ohm'],
                  'pot': ['Phi', 'electrostatic potential',
                         2, 'SM', 'sph', ['time', 'lon', 'lat'], 'kV'],
                  'delbr': ['DeltaB_r', 'r component of magnetic perturbaion on the ground',
                         3, 'SM', 'sph', ['time', 'lon', 'lat'], 'nT'],
                  'delbt': ['DeltaB_theta', 'theta (latitude pointing south) component of magnetic perturbaion on the ground',
                         3, 'SM', 'sph', ['time', 'lon', 'lat'], 'nT'],
                  'delbp': ['DeltaB_phi', 'phi (longitude) component of magnetic perturbaion on the ground',
                         3, 'SM', 'sph', ['time', 'lon', 'lat'], 'nT'],
                  'pacurr': ['JR', 'radial electric current',
                         3, 'SM', 'sph', ['time', 'lon', 'lat'], 'A/m**2'],
                  'xjh': ['Q_Joule', 'Joule heating rate',
                         3, 'SM', 'sph', ['time', 'lon', 'lat'], 'W/m**2'],
                  'ppio': ['P_mag', 'mapped plasma pressure',
                         3, 'SM', 'sph', ['time', 'lon', 'lat'], 'pPa'],
                  'rrio': ['N_mag', 'mapped plasma number density',
                         3, 'SM', 'sph', ['time', 'lon', 'lat'], '1/cm**3'],
                  'ttio': ['T_mag', 'mapped Plasma Temperature',
                         3, 'SM', 'sph', ['time', 'lon', 'lat'], 'K'],
                  'vdown': ['V_r', 'radial velocity',
                         3, 'SM', 'sph', ['time', 'lon', 'lat'], 'm/s'],
                  'xjh': ['JouleHeat', 'Joule heating rate',
                         3, 'SM', 'sph', ['time', 'lon', 'lat'], 'mW/cm**2'],
                  'prec_e_fe_1': ['Phi_eN01', 'electron number flux precipitation 1',
                                  12, 'GSE', 'car', ['time', 'lon', 'lat'], '1/cm**2/s'],
                  'prec_e_e0_1': ['Phi_eE01', 'electron energy flux precipitation 1',
                                  12, 'GSE', 'car', ['time', 'lon', 'lat'], 'mW/m**2'],
                  'prec_e_fe_2': ['Phi_eN02', 'electron number flux precipitation 2 ',
                                  12, 'GSE', 'car', ['time', 'lon', 'lat'], '1/cm**2/s'],
                  'prec_e_e0_2': ['Phi_eE02', 'electron energy flux precipitation 2 ',
                                  12, 'GSE', 'car', ['time', 'lon', 'lat'], 'mW/**2'],
                  }

def convert_all(file_dir,verbose=False,requested_variables=None):
    ''' convert all files independent of the others'''
    # file prefix includes everything including '3df' for magnetosphere and 'iof' for ionospehre electrodynamics and iof for ionosphere CTIM outputs
    # dates and times can be extracted from the files themselves, parsing teh ASCII date string written with each 2D and 3D field of data
    # Lutz Rastaetter - rewritten to convert single files only and obtain asciitime to generate date and time for NetCDF file name
    #
    import kamodo_ccmc.readers.reader_utilities as RU
    tic0 = perf_counter()

    file_path = dirname(file_dir)+sep

    file_names = sorted(RU.glob(file_path+".iof*"))
    
    if verbose:
        print('file_path: ',file_path)
        
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

    added_time_at_beginning=0
    added_time_at_end=0

    if requested_variables is None:
        convert_varnames=openggcm_ie_varnames
    else:
        convert_varnames=[]
        for varname in openggcm_ie_varnames:
            print(varname)
            if varname in requested_variables or openggcm_ie_varnames[varname][0] in requested_variables:
                convert_varnames.append(varname)                

    if verbose:
        print(convert_varnames)
    file_names=sorted(RU.glob(file_path+"*.iof*"))
    if verbose:
        print(file_names)
    nbuffer = 1000000
    for file in file_names:
        tic = perf_counter()
        varnames = [varname for varname in convert_varnames]
        fielddata=zeros(shape=(nbuffer),dtype=float32,order='F')
        results = ropgm.read_2d_field(file,fielddata,varnames[0],nbuffer);

        fieldarray,ny_field,nx_field,asciitime = results
        asciitime_str = asciitime.decode('utf-8')
        ifile = 0
        if verbose:
            print('asciitime_str: ',asciitime_str)
        timestr,time_seconds,time_epoch,time_str_colon,datatype = asciitime_str.split()
        if verbose:
            print("time_str_colon,datatype: ",time_str_colon,datatype)
            print(time_str_colon[0:10]+":00:00:00.000")
        filedate = datetime.strptime(time_str_colon[0:10]+":00:00:00.000", '%Y:%m:%d:%H:%M:%S.%f').replace(tzinfo=timezone.utc)
        
        nc_file = file[:-7]+"_"+time_str_colon[0:4]+"-"+time_str_colon[5:7]+"-"+time_str_colon[8:10]+"_"+time_str_colon[11:13]+"_"+time_str_colon[14:16]+"_"+time_str_colon[17:19]+".nc"
        print(nc_file)
            
        if exists(nc_file):
            remove(nc_file)
        else:
            print("The file %s does not exist" % nc_file) 
    
        # NetCDF root group
        data_out = Dataset(nc_file, 'w', format='NETCDF4')
        data_out.model = 'OpenGGCM_IE'
        data_out.file=file
        # data_out.file = ''.join([f+',' for f in file_names]).strip(',')  #csv list of files

        # store dimensions
        dim_dict={}
        gx = -90. + 180. * arange(nx_field)/(nx_field - 1) # degrees latitudes 
        gy = 360. * arange(ny_field)/(ny_field - 1.);      # degrees from noon local time
        
        dim_dict['lat']=array(gx,dtype=float32)
        dim_dict['lon']=array(gy,dtype=float32)
            
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
   
        ntime = 1 # len(file_names)
        ifile = 0
        nbuffer = 1000000
        data=zeros(shape=(ntime,ny_field,nx_field),dtype=float32,order='C')
        #    fielddata=zeros(shape=(nx,ny,nz),dtype=float32,order='F')
        missing_value=NaN
        time=zeros(ntime)
        #print(varnames)
        for varname in varnames:
            variable_found=True
            data[:]=0.
            fielddata=zeros(shape=(nbuffer),dtype=float32,order='F')
            fieldarray,ny_field,nx_field,asciitime = ropgm.read_2d_field(file,fielddata,varname,nbuffer);
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
            
            if variable_found:
                #            fieldarray[(fieldarray == openggcm_missing)]=missing_value
                #print(varname,fieldarray.shape,data.shape)
                fielddata=reshape(fieldarray[0:nx_field*ny_field],(nx_field,ny_field))
                data[ifile,:,:] = transpose(flip(fielddata,axis=0))
                #print(varname,data.min(),data.max(),variable_found)

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

                dims = openggcm_ie_varnames[varname][-2]
                #print(varname,dims,data.shape)
                new_var = data_out.createVariable(varname, float32, dims)
                new_var[:] = data 
                new_var.units=openggcm_ie_varnames[varname][-1]

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
        print(f"Data for {file} converted in {perf_counter()-tic:.6f}s.")
        data_out.close() 
        toc=perf_counter()
        print('converted file in ',toc-tic,' seconds')                       

        
    toc=perf_counter()
    print('converted files in ',toc-tic0,' seconds')                       

    return True



