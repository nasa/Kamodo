"""
@author: jmpettit
"""

import kamodo
import kamodo_ccmc.flythrough.model_wrapper as MW
import numpy as np
import kamodo_ccmc.flythrough.SatelliteFlythrough as SF
from astropy import units as u
from kamodo_ccmc.flythrough.utils import ConvertCoord

def fileforcing_oneway(model_A, model_B, input_filedir, output_filedir):

    if model_B == 'GITM':
        
        from kamodo_ccmc.filedriver import GITM_filedriver as driver 
        
    elif model_B == 'WACCMX':
        
        from kamodo_ccmc.filedriver import WACCMX_filedriver_v2 as driver
        
    elif model_B == 'TIEGCM':
        
        from kamodo_ccmc.filedriver import TIEGCM_filedriver as driver
        
    elif model_B == 'CTIPe':
        
        from kamodo.kamodo_ccmc.coupling import CTIPe_filedriver as driver
    
    else:
        
        raise AttributeError('Model not yet added.')

    # ------------------------------------------------------------------

    variable_list, coord_grids = driver.coordinate_grid_builder(model_A, input_filedir)
             
    # The model grid of model B will be needed (needs to be tested on other models)
    reader = MW.Model_Reader(model_A)
    ko_reader = reader(input_filedir, variables_requested = variable_list)        
    func = getattr(ko_reader, variable_list[0]+'_ijk')
            
    # Ror_www/output_files/kamodo_demo/
    from kamodo import get_defaults
    defaults = get_defaults(func)
    
    coord_list = list(defaults.keys()) 
    if len(coord_list) < 4:
        defaults['c3'] = 1.   # WACCMX, TIEGCM are on pressure grid (maybe choose middle point?)
        coord_list += ['c3']
    coord_data = [defaults[key] for key in coord_list[1:]]

    # Rewrite for generic grid
    UT_grid = defaults['time']
    UT_grid = UT_grid[4::5]

    time_grid = defaults['time'] * 3600. + ko_reader.filedate.timestamp()
    #time_grid = np.float64(defaults['time'])*3600.+ko_reader.filedate.timestamp()
    #time_grid = UT_grid[4::5] * 3600. + ko_reader.filedate.timestamp()

       
  # Make base grid
    
    x = np.linspace(-180, 180., 91)
    y = np.linspace(-90., 90., 91)
    xx, yy = np.array(np.meshgrid(x, y, indexing='xy'))
    x1 = np.reshape(xx, -1)
    y1 = np.reshape(yy, -1)
    z1 = np.full((len(x1)), 350.)
    values = np.zeros([len(UT_grid),len(variable_list), len(x1)])

    # ko is a Kamodo object, var is a string of the variable name
    for itime in range(len(UT_grid)):
        print (itime)
        for i in range(len(variable_list)):
            if variable_list[i] in ko_reader:
                interp = getattr(ko_reader, variable_list[i])
             
            # transform
            
            t1 = np.full(len(x1), UT_grid[itime])
            x2, y2, z2, units = ConvertCoord(t1, x1, y1, z1, 'SM', 'sph', 'GEO', 'sph')
            g2 = np.stack((t1, x2, y2), axis=-1)  # nx4 grid
             
            # Do the interpolation and reshape back into 2D array
            values[itime, i, :] = interp(g2)
  

    # ------------------------------------------------------------------        
          
    model_A_Dictionary, model_B_Dictionary = driver.variable_list(model_A)
    
    final_vars = values*0.
    
    A_list = list(model_A_Dictionary.keys())
    B_list = list(model_B_Dictionary.keys())
    A_units = [val[0] for key, val in model_A_Dictionary.items()]
    B_units = [val[0] for key, val in model_B_Dictionary.items()]

    for i in range(len(A_list)):
        initial_vars = values[:, i, :]
        unit_B = u.Unit(B_units[i])
        unit_A = u.Unit(A_units[i])
        
        try:
            results = (initial_vars * unit_A).to(unit_B)  # results.value to prevent astropy quantity
            final_vars[:, i, :] = results
            
        except:
            conv_factor = model_A_Dictionary[A_list[i]][1]
            #perform conversion
            results = initial_vars * conv_factor * unit_B # results.value to prevent astropy quantity
            final_vars[:, i, :] = results# Change to end_vars for compatibility (initialize outside)

    
   # ------------------------------------------------------------------        
 
    driver.oneway_forcing(model_A, final_vars, x, y, UT_grid, input_filedir, output_filedir)

   # ------------------------------------------------------------------  
      
    #return driver
    return final_vars