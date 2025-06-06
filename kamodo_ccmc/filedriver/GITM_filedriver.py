
"""
@author: jmpettit
"""

#GITM_File_Driven_Forcing


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
        
        model_B_Dictionary= {'V': ['V'], 'Phi_eE': ['erg/(s*cm**2)'], 'E_eavg': ['J']}
        model_A_Dictionary = {'phi': ['kV'], 'Phi_E': ['W/m**2'], 'E_avg': ['eV']}
        
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
    

    if model_A == 'SWMF_IE':

        # Inflate the arrays back to original shape using numpy reshape
        # Change variable names to Model B names instead of Model A names
        
        #import pdb
        #pdb.set_trace()
        
        Phi = np.reshape(end_vars[:, 0, :], [len(time_grid), len(lat_grid), len(lon_grid)])    
        Phi_E = np.reshape(end_vars[:, 1, :], [len(time_grid), len(lat_grid), len(lon_grid)])
        E_avg = np.reshape(end_vars[:, 2, :], [len(time_grid), len(lat_grid), len(lon_grid)])

        
        # Build the binary files
        
        from scipy.io import FortranFile

        MLTs = (lon_grid/15+12) % 24
        
        middleindex = middle_index(lat_grid)
        
        
        # Build Southern Hemisphere File
        
        lats = lat_grid[0:middleindex]
        latn = 90 - lat_grid[middleindex+1:-1]  
        lats = np.flip(lats)
        latn = np.flip(latn)

               
        output_filename = str(SF.File_UTCTimes(model_A, input_filedir)[2]) # Returning midnight
        output_filename = output_filename[0:10]+'_'+output_filename[11:-1]
        output_filename = output_filename.replace(':', '_')
        
        output_filename_record = output_filename[0:16].replace('-', ' ')
        output_filename_record = output_filename_record.replace('_', ' ')
        
        filedate = output_filename_record[0:10]
        filedate = filedate.replace(" ", "")
                
        f = FortranFile(output_filedir+'/'+str(model_A)+'_to_GITM_'+output_filename+'s.swmf', 'w')
        
        # Hemispheric Power Index (Setting Kp index to 3 for testing)

        KpData = np.load('KP_Daily.npz')
        kp_date = KpData['date_daily']
        KpIndex = KpData['kp_daily']
        KP_Date = ["" for x in range(len(kp_date))]
        
        for i in range(0, len(kp_date)-1):
            
            KP_Date[i] = kp_date[i, 0]+kp_date[i, 1]+kp_date[i, 2]
            
            if str(KP_Date[i]) == filedate:
                
                print ("KP data successfully found!", filedate, KP_Date[i])
                
                index = i
                
                HPI = 16.82*2.7182**KpIndex[index] - 4.86
                
            else: 
                
                # Assume a Kp-index of 3
                
                HPI = 16.82*2.7182**3 - 4.86
        

        HPI = np.repeat(HPI, len(time_grid))


        # Record record is putting in integers, not double/floats
        

        f.write_record([len(latn), len(MLTs), len(time_grid)-2])
        f.write_record(np.array(lats).reshape((len(lats))))
        f.write_record(np.array(MLTs).reshape((len(MLTs))))  
        f.write_record([len(end_vars)])
        
        
        potential_string = np.array('Potential (V)                 ', dtype='|S30')
        f.write_record(potential_string)
        Energy_flux_string = np.array('Energy Flux (ergs/cm2)        ', dtype='|S30')
        f.write_record(Energy_flux_string)
        Mean_Energy_string = np.array('Mean Energy (eV)              ', dtype='|S30')
        f.write_record(Mean_Energy_string)
        
        
        for itime in range(0, len(time_grid)-3):
            #print (itime)

            #f.write_record(itime)
            Phi_temp = Phi[:, itime, 0:middleindex]
            Phi_E_temp = Phi_E[:, itime, 0:middleindex]
            E_avg_temp = E_avg[:, itime, 0:middleindex]

            
            f.write_record(output_filename_record)
            f.write_record(np.array(Phi_temp, dtype='float32'))
            f.write_record(np.array(Phi_E_temp, dtype='float32'))
            f.write_record(np.array(E_avg_temp, dtype='float32'))
            f.write_record(1)

        f.close()
        
        f = FortranFile(output_filedir+'/'+str(model_A)+'_to_GITM_'+output_filename+'n.swmf', 'w')
        
        
        # Build Northern Hemisphere file

        f.write_record([len(latn), len(MLTs), len(time_grid)-2])
        f.write_record(latn)
        f.write_record(MLTs)   
        f.write_record(45)
        
        potential_string = np.array('Potential (V)                 ', dtype='|S30')
        f.write_record(potential_string)
        Energy_flux_string = np.array('Energy Flux (ergs/cm2)        ', dtype='|S30')
        f.write_record(Energy_flux_string)
        Mean_Energy_string = np.array('Mean Energy (eV)              ', dtype='|S30')
        f.write_record(Mean_Energy_string)
        
        
        for itime in range(0, len(time_grid)-3):

            #f.write_record(itime)
            Phi_temp = Phi[:, itime, middleindex:-1]
            Phi_E_temp = Phi_E[:, itime, middleindex:-1]
            E_avg_temp = E_avg[:, itime, middleindex:-1]

            f.write_record(output_filename_record)
            f.write_record(np.array(Phi_temp, dtype='float32'))
            f.write_record(np.array(Phi_E_temp, dtype='float32'))
            f.write_record(np.array(E_avg_temp, dtype='float32'))
            f.write_record(1)
        
        f.close()
        
                                      
    else: 
        raise AttributeError('Model A not yet added.')

    return

def middle_index(lat_grid):
    
  if (len(lat_grid) % 2 == 0):
    even_index = int(len(lat_grid)/2)-1
   
    return even_index
  
  elif (len(lat_grid) % 2 != 0):
    odd_index = int(len(lat_grid)/2)
    
    return odd_index
