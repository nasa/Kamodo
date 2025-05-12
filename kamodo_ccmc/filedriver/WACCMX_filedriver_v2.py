"""
@author: jmpettit

WACCMX filebuilder that takes in model data from Model A and builds a 
WACCMX ready AMIE file.

"""

def coordinate_grid_builder(model_A, input_filedir):
    
    # Build the coordinate grid using the Model B as input (GITM) 
    # Output the coordinate grid into a dictionary 

    import kamodo
    import kamodo_ccmc.flythrough.model_wrapper as MW
        
    if model_A == 'SWMF_IE':
            
        variable_list = ['phi', 'Phi_E', 'E_avg']  # SWMF_IE
        coord_grids = MW.Model_Variables(model_A, return_dict=True)
        # Sigma_H is always present so coord_grids will always return
        coord_grids = coord_grids['Sigma_H'][2:4]  
        
    elif model_A == 'AMGeO':
        
        # AMGeO data will be converted to precipitation
        variable_list = ['V']
        coord_grids = MW.Model_Variables(model_A, return_dict=True)
        coord_grids = coord_grids['Sigma_H'][2:4]
        
    elif model_A == 'Ovation-Prime':
        
        variable_list = ['']
        coord_grids = MW.Model_Variables(model_A, return_dict=True)
        coord_grids == coord_grids['Sigma_H'][2:4]
        
    elif model_A == 'SuperDARN':
        
        variable_list = ['kV']
        coord_grids = MW.Model_variables(model_A, return_dict=True)
        coord_grids == coord_grids['Sigma_H'][2:4]
    
    else:
        raise AttributeError('Model A not yet added.')      
    
    return variable_list, coord_grids


def variable_list(model_A):
    
    # map the variable list from model A to model B using dictionary, return variable list

    if model_A == 'SWMF_IE':
               
        model_B_Dictionary= {'V': ['V'], 'Phi_eE': ['W/m**2'], 'E_eavg': ['keV']}
        model_A_Dictionary = {'phi': ['kV'], 'Phi_E': ['W/m**2'], 'E_avg': ['keV']}
        
    elif model_A == 'AMGeO':
        
        model_B_Dictionary= {'V': ['V']}
        model_A_Dictionary = {'V': ['V']}
        
    elif model_A == 'SuperDARN':
         
        model_B_Dictionary= {'V': ['V']}
        model_A_Dictionary = {'kV': ['kV']}
        
    else:
        raise AttributeError('Model A not yet added.')

    print ('Variable list built')


    return model_A_Dictionary, model_B_Dictionary


def oneway_forcing(model_A, end_vars, lon_grid, lat_grid, time_grid, input_filedir, output_filedir):
    
    # Output model A variables mapped to model B grid into format Model B can read
    
    import numpy as np
    import kamodo_ccmc.flythrough.SF_utilities as SF
    from netCDF4 import Dataset
    import pandas as pd
    

    if model_A == 'SWMF_IE':

        # Inflate the arrays back to original shape using numpy reshape
        # Change variable names to Model B names instead of Model A names
        
        # Build the netcdf files       ..        
               
        output_filename = str(SF.File_UTCTimes(model_A, input_filedir)[2]) # Returning midnight
                
        year_temp = output_filename[0:4]
        month_temp = output_filename[5:7]
        day_temp = output_filename[8:10]
        output_filename = output_filename[0:10]+'_'+output_filename[11:-1]
        output_filename = output_filename.replace(':', '_')
        
        output_filename_record = output_filename[0:16].replace('-', ' ')
        output_filename_record = output_filename_record.replace('_', ' ')
        
        filedate = output_filename_record[0:10]
        filedate = filedate.replace(" ", "")
        
        NH_lat_temp = np.where(lat_grid > 0.)
        NH_lat = lat_grid[NH_lat_temp]
        
        SH_lat_temp = np.where(lat_grid < 0.)
        SH_lat = lat_grid[SH_lat_temp]
        
        Phi = np.reshape(end_vars[:, 0, :], [len(time_grid), len(lat_grid), len(lon_grid)])    
        Phi_E = np.reshape(end_vars[:, 1, :], [len(time_grid), len(lat_grid), len(lon_grid)])
        E_avg = np.reshape(end_vars[:, 2, :], [len(time_grid), len(lat_grid), len(lon_grid)])

       # Create 5-min time array 
        
        WACCMX_Time = np.arange(0, 23.92, 0.083333333333)
        Phi_final = np.empty([len(WACCMX_Time), len(lat_grid), len(lon_grid)])
        Phi_E_final = np.empty([len(WACCMX_Time), len(lat_grid), len(lon_grid)])
        E_avg_final = np.empty([len(WACCMX_Time), len(lat_grid), len(lon_grid)])            
        
        
        for itime in range(len(WACCMX_Time)-1):
            
            close_time = abs(WACCMX_Time[itime] - time_grid)
            closest_time = np.argmin(close_time)
            Phi_final[itime, :, :] = Phi[closest_time, :, :]
            Phi_E_final[itime, :, :] = Phi_E[closest_time, :, :]
            E_avg_final[itime, :, :] = E_avg[closest_time, :, :]
                
           
        
        Phi_SH = np.squeeze(Phi_final[:, SH_lat_temp, :])
        Phi_NH = np.squeeze(Phi_final[:, NH_lat_temp, :])
        
        Phi_E_SH = np.squeeze(Phi_E_final[:, SH_lat_temp, :])
        Phi_E_NH = np.squeeze(Phi_E_final[:, NH_lat_temp, :])
        
        E_avg_SH = np.squeeze(E_avg_final[:, SH_lat_temp, :])
        E_avg_NH = np.squeeze(E_avg_final[:, NH_lat_temp, :])
        
        # create variable dimensions       
        
        # Southern Hemisphere Output
        
        WACCMX_netcdf_file = Dataset(output_filedir+str(model_A)+'_to_WACCMX_'+output_filename[0:10]+'_sh.nc', 'w', format='NETCDF4')
        WACCMX_netcdf_file.description = 'Example simulation data'
                                
        
        #Time-Index = np.where(UT_grid[3] eq WACCMX_Time)
        
        WACCMX_netcdf_file.createDimension('lon', len(lon_grid))
        WACCMX_netcdf_file.createDimension('time', None)
        WACCMX_netcdf_file.createDimension('lat', len(SH_lat[0:26]))
        
        # create variables
        
        ut = WACCMX_netcdf_file.createVariable('ut', 'f8', ('time'))
        ut.long_name = "universal time (in decimal hour)"
        ut.units = "hours"
        year = WACCMX_netcdf_file.createVariable('year', 'i4', ('time'))
        year.long_name = "year"
        month = WACCMX_netcdf_file.createVariable('month', 'i4', ('time'))
        month.long_name = "month"
        day = WACCMX_netcdf_file.createVariable('day', 'i4', ('time'))
        day.long_name = "day"
        jday = WACCMX_netcdf_file.createVariable('jday', 'i4', ('time'))
        jday.long_name = "Julian day"
        lat = WACCMX_netcdf_file.createVariable('lat', 'f8', ('lat',))
        lat.long_name = "APEX latitude (south or north)"
        lat.units = "degrees"
        lon = WACCMX_netcdf_file.createVariable('lon', 'f8', ('lon',))
        lon.long_name = "APEX longtiude (0-360)"
        lon.units = "degrees"
        Phi_E_def = WACCMX_netcdf_file.createVariable('efx', 'f8', ('time', 'lat', 'lon' ))
        Phi_E_def.long_name = "SWMF_IE Energy Flux"
        Phi_E_def.units = "mW m-2"
        Phi_def = WACCMX_netcdf_file.createVariable('pot', 'f8', ('time', 'lat', 'lon' ))
        Phi_def.long_name = "SWMF_IE Electric Potential"
        Phi_def.units = "Volt"
        E_avg_def = WACCMX_netcdf_file.createVariable('ekv', 'f8', ('time', 'lat', 'lon' ))
        E_avg_def.long_name = "SWMF_IE Mean Energy"
        E_avg_def.units = "keV"
        hpi_def = WACCMX_netcdf_file.createVariable('hpi', 'f8', ('time')) 
        hpi_def.long_name = "hemispheric integrated power"
        hpi_def.units = "GW"
        pcp_def = WACCMX_netcdf_file.createVariable('pcp', 'f8', ('time')) 
        pcp_def.long_name = "cross-polar-cap potential drop"
        pcp_def.units = "kV"
        cuspmlt_def = WACCMX_netcdf_file.createVariable('cuspmlt', 'f8', ('time')) 
        cuspmlt_def.long_name = "cusp MLT"
        cuspmlt_def.units = "degree"
        cusplat_def = WACCMX_netcdf_file.createVariable('cusplat', 'f8', ('time')) 
        cusplat_def.long_name = "cusp latitude"
        cusplat_def.units = "degree"
               
                
        # fill variables

        ut[:] = WACCMX_Time
        month[:] = month_temp
        day[:] = day_temp
        year[:] = year_temp
        
        # Hemispheric Power Index (Setting Kp index to 3 for testing)
        
        KpData = np.load(input_filedir+'KP_Daily.npz')
        kp_date = KpData['date_daily']
        KpIndex = KpData['kp_daily']
        KP_Date = ["" for x in range(len(kp_date))]
        
        for i in range(0, len(kp_date)-1):
            
            KP_Date[i] = kp_date[i, 0]+kp_date[i, 1]+kp_date[i, 2]
            
            if str(KP_Date[i]) == filedate:
                
                print ("KP data successfully found!", filedate, KP_Date[i])
                
                index = i
                
                HPI = 16.82*2.7182**(KpIndex[index] - 4.86)
                print ("HPI:", HPI)
                
            else: 
                
                # Assume a Kp-index of 3
                
                HPI = 16.82*2.7182**(3 - 4.86)
        

        HPI = np.repeat(HPI, len(ut))
        hpi_def[:] = HPI        
                
        Phi_E_def[:,:,:] = Phi_E_SH[:, 0:26, :]
        Phi_def[:,:,:] = Phi_SH[:, 0:26, :]
        E_avg_def[:,:,:] = E_avg_SH[:, 0:26, :]
        lat[:] = (SH_lat[0:26])
        lon[:] = lon_grid+180.
        
        # Compute cross-polar cap potential
        
        for i in range(0, len(WACCMX_Time)):
            pcp_def[i] = abs(np.max(Phi_SH[i,0:26,:]))+abs(np.min(Phi_SH[i,0:26,:]))
            jday[i] = pd.Period(year_temp+"-"+month_temp+"-"+day_temp).day_of_year

        # Assume a constant cusp MLT and cusp Latitude   
        
        cusplat_temp = np.repeat(65., len(ut))
        cusplat_def[:] = cusplat_temp
        
        cuspmlt_temp = np.repeat(12, len(ut))
        cuspmlt_def[:] = cuspmlt_temp
                
        WACCMX_netcdf_file.close    

        
        # Northern Hemisphere Output        
        
        WACCMX_netcdf_file = Dataset(output_filedir+str(model_A)+'_to_WACCMX_'+output_filename[0:10]+'_nh.nc', 'w', format='NETCDF4')
        WACCMX_netcdf_file.description = 'Example simulation data'
        
        # create variable dimensions
        
        WACCMX_netcdf_file.createDimension('lon', len(lon_grid))
        WACCMX_netcdf_file.createDimension('time', None)
        WACCMX_netcdf_file.createDimension('lat', len(NH_lat[19:45]))
        
        # create variables

        ut = WACCMX_netcdf_file.createVariable('ut', 'f8', ('time'))
        ut.long_name = "universal time (in decimal hour)"
        ut.units = "hours"
        year = WACCMX_netcdf_file.createVariable('year', 'i4', ('time'))
        year.long_name = "year"
        month = WACCMX_netcdf_file.createVariable('month', 'i4', ('time'))
        month.long_name = "month"
        day = WACCMX_netcdf_file.createVariable('day', 'i4', ('time'))
        day.long_name = "day"
        jday = WACCMX_netcdf_file.createVariable('jday', 'i4', ('time'))
        jday.long_name = "Julian day"
        lat = WACCMX_netcdf_file.createVariable('lat', 'f8', ('lat',))
        lat.long_name = "APEX latitude (south or north)"
        lat.units = "degrees"
        lon = WACCMX_netcdf_file.createVariable('lon', 'f8', ('lon',))
        lon.long_name = "APEX longtiude (0-360)"
        lon.units = "degrees"
        Phi_E_def = WACCMX_netcdf_file.createVariable('efx', 'f8', ('time', 'lat', 'lon' ))
        Phi_E_def.long_name = "SWMF_IE Energy Flux"
        Phi_E_def.units = "mW m-2"
        Phi_def = WACCMX_netcdf_file.createVariable('pot', 'f8', ('time', 'lat', 'lon' ))
        Phi_def.long_name = "SWMF_IE Electric Potential"
        Phi_def.units = "Volt"
        E_avg_def = WACCMX_netcdf_file.createVariable('ekv', 'f8', ('time', 'lat', 'lon' ))
        E_avg_def.long_name = "SWMF_IE Mean Energy"
        E_avg_def.units = "keV"
        hpi_def = WACCMX_netcdf_file.createVariable('hpi', 'f8', ('time')) 
        hpi_def.long_name = "hemispheric integrated power"
        hpi_def.units = "GW"
        pcp_def = WACCMX_netcdf_file.createVariable('pcp', 'f8', ('time')) 
        pcp_def.long_name = "cross-polar-cap potential drop"
        pcp_def.units = "kV"
        cuspmlt_def = WACCMX_netcdf_file.createVariable('cuspmlt', 'f8', ('time')) 
        cuspmlt_def.long_name = "cusp MLT"
        cuspmlt_def.units = "degree"
        cusplat_def = WACCMX_netcdf_file.createVariable('cusplat', 'f8', ('time')) 
        cusplat_def.long_name = "cusp latitude"
        cusplat_def.units = "degree"
        
        
        # fill variables
        

        Phi_E_def[:,:,:] = Phi_E_NH[:, 19:45, :]
        Phi_def[:,:,:] = Phi_NH[:, 19:45, :]
        E_avg_def[:,:,:] = E_avg_NH[:, 19:45, :]
        lat[:] = (NH_lat[19:45])
        lon[:] = lon_grid+180.
               
        ut[:] = WACCMX_Time
        month[:] = month_temp
        day[:] = day_temp
        year[:] = year_temp
            
        # Hemispheric Power Index 
        # Since we are parameterizing it based on Kp, it will be the same for both hemispheres
        
        hpi_def[:] = HPI
        
        # Compute cross-polar cap potential
        
        for i in range(0, len(WACCMX_Time)):
            pcp_def[i] = abs(np.max(Phi_NH[i,19:45,:]))+abs(np.min(Phi_NH[i,19:45,:]))
            jday[i] = pd.Period(year_temp+"-"+month_temp+"-"+day_temp).day_of_year
        
        # Assume a constant cusp MLT and cusp Latitude
        
        cusplat_temp = np.repeat(65., len(ut))
        cusplat_def[:] = cusplat_temp
        
        cuspmlt_temp = np.repeat(12, len(ut))
        cuspmlt_def[:] = cuspmlt_temp
        
        WACCMX_netcdf_file.close
        
                                      
    else: 
        raise AttributeError('Model A not yet added.')

    return
