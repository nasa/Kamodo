# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 15:37:59 2023

@author: rringuet,lrastaet
"""
from numpy import vectorize
from datetime import datetime, timezone

model_varnames = {'Bx': ['B_x', 'X-component of magnetic field',
                         0, 'SM', 'car', ['time', 'X', 'Y', 'Z'], 'nT'],
                  'By': ['B_y', 'Y-component of magnetic field',
                         0, 'SM', 'car', ['time', 'X', 'Y', 'Z'], 'nT'],
                  'Bz': ['B_z', 'Z-component of magnetic field',
                         0, 'SM', 'car', ['time', 'X', 'Y', 'Z'], 'nT'],
                  'Cs': ['c_s', 'Sound speed of plasma', 0, 'SM', 'car',
                         ['time', 'X', 'Y', 'Z'], 'km/s'],
                  'D': ['N_plasma', 'Plasma number density (M/mp)', 0, 'SM',
                        'car', ['time', 'X', 'Y', 'Z'], '1/cm**3'],
                  'Jx': ['J_x', 'X-component of current density',
                         0, 'SM', 'car', ['time', 'X', 'Y', 'Z'], 'nA/m**2'],
                  'Jy': ['J_y', 'Y-component of current density',
                         0, 'SM', 'car', ['time', 'X', 'Y', 'Z'], 'nA/m**2'],
                  'Jz': ['J_z', 'Z-component of current density',
                         0, 'SM', 'car', ['time', 'X', 'Y', 'Z'], 'nA/m**2'],
                  'P': ['P', 'Pressure', 0, 'SM', 'car',
                        ['time', 'X', 'Y', 'Z'], 'nPa'],
                  'Pb': ['P_mag', 'Magnetic pressure', 0, 'SM', 'car',
                         ['time', 'X', 'Y', 'Z'], 'nPa'],
                  'Vx': ['v_x', 'X-component of velocity',
                         0, 'SM', 'car', ['time', 'X', 'Y', 'Z'], 'km/s'],
                  'Vy': ['v_y', 'Y-component of velocity',
                         0, 'SM', 'car', ['time', 'X', 'Y', 'Z'], 'km/s'],
                  'Vz': ['v_z', 'Z-component of velocity',
                         0, 'SM', 'car', ['time', 'X', 'Y', 'Z'], 'km/s'],
                  # CONSTANTS below this line.
                  'dV': ['dV', 'Simulation cell volume', 0, 'SM', 'car',
                         ['X', 'Y', 'Z'], 'R_E**3'],
                  'X': ['X', 'X grid position', 0, 'SM', 'car',
                         ['X', 'Y', 'Z'], 'R_E'],
                  'Y': ['Y', 'Y grid position', 0, 'SM', 'car',
                         ['X', 'Y', 'Z'], 'R_E'],
                  'Z': ['Z', 'Z grid position', 0, 'SM', 'car',
                         ['X', 'Y', 'Z'], 'R_E'],
                  }

#constants = ['X', 'Y', 'Z', 'dV']
constants = []

@vectorize
def timestr_timestamp(time_str):
    '''Converts time string from data file into an UTC timestamp. Cuts off
    millisecond portion of times to prevent errors later.'''
    dt = datetime.fromisoformat(time_str).replace(tzinfo=timezone.utc)
    return dt.timestamp()


def MODEL():
    from kamodo import Kamodo
    from os.path import basename
    from astropy.time import Time
    from datetime import datetime, timezone
    from numpy import array, NaN, ndarray, zeros, linspace, repeat, float32, float64, int32, sqrt, arctan2, arccos, mod, pi
    from time import perf_counter
    import kamodo_ccmc.readers.reader_utilities as RU
    import kamodo_ccmc.readers.gameragm_grids as G

    from readers.Tri2D._interpolate_tri2d import ffi as tri2d_ffi
    from readers.Tri2D._interpolate_tri2d import lib as tri2d_lib
    
    class MODEL(Kamodo):
        '''GAMERA GM model data reader.

        Inputs:
            file_dir: a string representing the file directory of the
                model output data.
                Note: This reader 'walks' the entire dataset in the directory.
            variables_requested = a list of variable name strings chosen from
                the model_varnames dictionary in this script, specifically the
                first item in the list associated with a given key.
                - If empty, the reader functionalizes all possible variables
                    (default)
                - If 'all', the reader returns the model_varnames dictionary
                    above for only the variables present in the given files.
            filetime = boolean (default = False)
                - If False, the script fully executes.
                - If True, the script only executes far enough to determine the
                    time values associated with the chosen data.
            printfiles = boolean (default = False)
                - If False, the filenames associated with the data retrieved
                    ARE NOT printed.
                - If True, the filenames associated with the data retrieved ARE
                    printed.
            gridded_int = boolean (default = True)
                - If True, the variables chosen are functionalized in both the
                    standard method and a gridded method.
                - If False, the variables chosen are functionalized in only the
                    standard method.
            verbose = boolean (False)
                - If False, script execution and the underlying Kamodo
                    execution is quiet except for specified messages.
                - If True, be prepared for a plethora of messages.
        All inputs are described in further detail in
            KamodoOnboardingInstructions.pdf.

        Returns: a kamodo object (see Kamodo core documentation) containing all
            requested variables in functionalized form.

        Notes and instructions:
        - The GAMERA global magnetosphere outputs are given in one or more h5
          files, each containing all of the time steps for the entire run (or several run segments).
          If the model is run in serial mode, only one file will be produced per segment.
          If the model is run in MPI mode, then multiple files will be produced
          with the grid of the simulation sectioned off into one piece per file.
          No file conversion is attempted, but some pre-processing is
          performed on the coordinate grids.
        - The files are typically larger than 16 GB, so interpolation method 3
          is chosen for the time-varying variables. Interpolation method 0 is
          chosen for the constants (e.g. dV). The coordinate grid is not
          uniform in any manner, so a custom interpolator is required for both
          interpolation methods.
        - Adding this interpolator to the reader is a work in progress via a
          collaboration with the GAMERA modeling team.
        '''
        def __init__(self, file_dir, variables_requested=[],
                     printfiles=False, filetime=False, gridded_int=True,
                     verbose=False, **kwargs):
            super(MODEL, self).__init__(**kwargs)
            self.modelname = 'GAMERA_GM'
            t0 = perf_counter()

            # first, check for file list, create if DNE
            list_file = file_dir + self.modelname + '_list.txt'
            time_file = file_dir + self.modelname + '_times.txt'
            self.times, self.pattern_files = {}, {}
            if not RU._isfile(list_file) or not RU._isfile(time_file):
                # collect filenames
                files = [f for f in sorted(RU.glob(file_dir+'*.h5')) if 'gam' in f]  # # if 'Res' not in f]
                if len(files) == 0:
                    files=[f for f in sorted(RU.glob(file_dir+'output_?/*.h5'))  if 'gam' in f]  # no more than 10 output_? directories, please
                    files.append([f for f in sorted(RU.glob(file_dir+'output_??/*.h5')) if 'gam' in f]) # no more than 90 more output_?? directories, please
                data_files = [f for f in files if 'Res' not in f and len(f) > 0]

                self.filename = ''.join([f+',' for f in data_files[:-1]])
                # one pattern per run: abcd_00ni_00nj_00nk
                # abcd part might not always be four characters
                p = ''.join([b+'_' for b in basename(
                    data_files[0]).split('.')[0].split('_')[:-3]])[:-1] # e.g. 'msphere_0008_0008_0001' (MPI runs)
                if p == '':
                    p = ''.join([b+'_' for b in basename(
                        data_files[0]).split('.')[0]])[:-1]  # e.g., 'msphere' (serial runs)
                    self.num_files_per_set = 1
                    self.ni_set = 1
                    self.nj_set = 1
                    self.nk_set = 1
                else:
                    self.ni_set,self.nj_set,self.nk_set = [int(j) for j in (p.split('_'))[-3:] ]  # e.g., 8, 8, 1
                    self.num_files_per_set = self.ni_set*self.nj_set*self.nk_set
                
                # establish time attributes
                # get list of files to loop through later
                self.pattern_files[p] = data_files
                self.times[p] = {'start': [], 'end': [], 'all': []}

                # all times are in each file, so just use first file of each set
                for ifile in range(0,len(data_files),self.num_files_per_set):
                    h5_data = RU.h5py(data_files[ifile])
                    timestep_keys = [key for key in h5_data.keys()
                                     if 'Step' in key]
                    mjd_list = [h5_data[key].attrs['MJD'] for key in
                                timestep_keys]
                    h5_data.close()

                    # convert from modified julian date to UTC time
                    # output format is 'YYYY-MM-DD HH:MM:SS.mms
                    t = sorted(Time(mjd_list, format='mjd').iso)
                    if ifile == 0:   # filedate -- date of start time from first file
                        self.filedate = datetime.fromisoformat(t[0][:10]).replace(
                            tzinfo=timezone.utc)  # date only
                    # convert to hrs since midnight
                    utc_ts = timestr_timestamp(t)
                    hrs = (utc_ts - self.filedate.timestamp()) / 3600.
                    # store in self.times dictionary
                    # avoid overlap when listing time steps in multiple sets of files
                    if ifile == 0:
                        i_start = 0
                    else:
                        i_start = [i_ for i_ in range(len(hrs)) if hrs[i_] == self.times['end']]
                        if len(i_start) > 0:
                            i_start = i_start[0]
                        else:
                            i_start = 0
                                                                 
                    self.times[p]['start'].append(hrs[0])
                    self.times[p]['end'].append(hrs[-1])
                    self.times[p]['all'].extend(list(hrs[i_start:])) # avoid overlap
                    
                    
                # need to figure out if millisecond accuracy is needed
                test = (len(utc_ts) - 1)/(utc_ts[-1] - utc_ts[0])  # #steps/s
                if test >= 1.5:  # need ms timing if >= 1.5 timesteps/second
                    ms_timing = True
                else:
                    ms_timing = False

                # create time list file if DNE
                RU.create_timelist(list_file, time_file, self.modelname,
                                   self.times, self.pattern_files,
                                   self.filedate, ms_timing=ms_timing,num_files_per_set=self.num_files_per_set)
            else:  # read in data and time grids from file list
                self.times, self.pattern_files, self.filedate, self.filename =\
                    RU.read_timelist(time_file, list_file)


            if filetime:
                return  # return times only

            # perform initial check on variables_requested list
            if len(variables_requested) > 0 and variables_requested != 'all':
                test_list = [value[0] for key, value in model_varnames.items()]
                err_list = [item for item in variables_requested if item
                            not in test_list]
                if len(err_list) > 0:
                    print('Variable name(s) not recognized:', err_list)
                for item in err_list:
                    variables_requested.remove(item)
                if len(variables_requested) == 0:
                    return

            # collect variable list (in attributes of datasets)
            p = list(self.pattern_files.keys())[0]  # only one pattern
            h5_data = RU.h5py(self.pattern_files[p][0])
            key_list = [key for key in h5_data.keys() if 'Step' not in key]  
                        # key not in ['X', 'Y', 'Z']]  # skip coordinates
            step_list = list(h5_data['Step#0'].keys())
            var_list = key_list + step_list
            h5_data.close()
            data_files = self.pattern_files[p]
            if len(p.split('_')) < 4:
                # p is 'msphere' or 'geo_mpi'
                self.num_files_per_set = 1
                self.ni_set = 1
                self.nj_set = 1
                self.nk_set = 1
            else:
                # p is 'mpshere_IIII_JJJJ_KKKK' or 'geo_mpi_iiii_jjjj_kkkk'
                self.ni_set,self.nj_set,self.nk_set = [int(j) for j in (p.split('_'))[-3:] ]  # e.g., 8, 8, 1
                self.num_files_per_set = self.ni_set*self.nj_set*self.nk_set
                
         
            if len(variables_requested) > 0 and variables_requested != 'all':
                gvar_list = [key for key, value in model_varnames.items()
                             if value[0] in variables_requested and
                             key in var_list]  # file variable names

                # check for variables requested but not available
                if len(gvar_list) != len(variables_requested):
                    err_list = [value[0] for key, value in
                                model_varnames.items() if value[0] in
                                variables_requested and key not in gvar_list]
                    if len(err_list) > 0:
                        print('Some requested variables are not available:',
                              err_list)
            else:  # only input variables on the avoid_list if requested
                gvar_list = [key for key in var_list if key in
                             model_varnames.keys()]
                # return list of variables included in data files
                if variables_requested == 'all':
                    self.var_dict = {value[0]: value[1:] for key, value in
                                     model_varnames.items()
                                     if key in gvar_list}
                    return

            # initialize data mapping for each variable desired
            self.variables = {model_varnames[var][0]:
                              {'units': model_varnames[var][-1],
                               'data': p} for var in gvar_list}

            # read in coordinate grid ranges and make sample grids
            # X, Y, and Z are the same shape in the data because they are the
            # x, y, and z coordinates of the same set of points.
            #
            # prepare (or read) cylindrical grid setup and triangulation in cylindrical grid
            #
            data_files = (self.pattern_files[p])[0]

            data_files = RU.glob('_'.join(data_files.split('_')[0:-3])+'*.h5')

            net_indices = G.GridMapping(data_files[0:self.num_files_per_set],True,'X')
            X = G.AssembleGrid(data_files[0:self.num_files_per_set],net_indices,True,'X')
            Y = G.AssembleGrid(data_files[0:self.num_files_per_set],net_indices,True,'Y')
            Z = G.AssembleGrid(data_files[0:self.num_files_per_set],net_indices,True,'Z')

            G.GridRanges(data_files[0:self.num_files_per_set])

            X_min, X_max, Y_min, Y_max, Z_min, Z_max =\
                G.Read_GridRanges(file_dir)
            if len(self.pattern_files[p]) > 1:
                sample_var = [gvar for gvar in model_varnames.keys()
                              if gvar in step_list][0]
                net_indices = G.GridMapping(self.pattern_files[p], False,
                                            sample_var)
                shape = net_indices['net'][1]  # shape of total array
            else:
                h5_data = RU.h5py(self.pattern_files[p][0])
                shape = list(h5_data['X'].shape)  # shape of total array
                h5_data.close()
            # grid for making slices
            # sizes of _x,_y,_Z reflect the aspect ratio of the box
            # 3 times coarser near the Earth but decent overall
            self._X = linspace(X_min, X_max, endpoint=True, num=6*shape[0]) 
            self._Y = linspace(Y_min, Y_max, endpoint=True, num=4*shape[0]) 
            self._Z = linspace(Z_min, Z_max, endpoint=True, num=4*shape[0]) 

            # store a few items
            self.missing_value = NaN
            if verbose:
                print(f'Took {perf_counter()-t0:.6f}s to read in data')
            if printfiles:
                print(self.filename)
            #
            # prepare (or read) cylindrical grid setup and triangulation in cylindrical grid
            #

            net_indices = G.GridMapping(data_files[0:self.num_files_per_set],True,'X')
            X = G.AssembleGrid(data_files[0:self.num_files_per_set],net_indices,True,'X')
            Y = G.AssembleGrid(data_files[0:self.num_files_per_set],net_indices,True,'Y')
            Z = G.AssembleGrid(data_files[0:self.num_files_per_set],net_indices,True,'Z')
            X2D,R2D,PHI = G.PrepareCylindricalGrid(X,Y,Z)
            RAD2D = sqrt(X2D*X2D+R2D*R2D)
            THETA2D = arccos(X2D/RAD2D) # THETA2D = 0 is along the negative X axis!
            tri_vertices,tri_neighbors,tri_ij = G.PrepareTriangularInterpolation(X2D.shape[1],X2D.shape[0])

            # need to add overlay grid (*_sg_*)
            xmin_sg,dx_sg,nx_sg,rmin_sg,dr_sg,nr_sg,start_index_sg = G.PrepareOverlaySphGrid2D(X2D,R2D,tri_vertices)

            nth_grid,nr_grid = X2D.shape
            nphi_grid = len(PHI)

            self._gridcorners = {'X':X,
                                 'Y':Y,
                                 'Z':Z,
                                 'net_indices':net_indices,
                                 'X2D':X2D,
                                 'R2D':R2D,
                                 'RAD2D':RAD2D,
                                 'THETA2D':THETA2D,
                                 'PHI':PHI,
                                 'n1':nr_grid,
                                 'n2':nth_grid,
                                 'nphi':nphi_grid,
                                 'ntri':(nth_grid-1)*(nr_grid-1)*2,
                                 'tri_vertices':tri_vertices,
                                 'tri_neighbors':tri_neighbors,
                                 'tri_ij':tri_ij,
                                 'xmin_sg':xmin_sg,
                                 'dx_sg':dx_sg,
                                 'nx_sg':nx_sg,
                                 'rmin_sg':rmin_sg,
                                 'dr_sg':dr_sg,
                                 'nr_sg':nr_sg,
                                 'start_index_sg':start_index_sg,
                                 # pointers
                                 'X2D_p':tri2d_ffi.new("float[]",list(X2D.flatten())),
                                 'R2D_p':tri2d_ffi.new("float[]",list(R2D.flatten())),
                                 'PHI_p':tri2d_ffi.new("float[]",list(PHI)),
                                 'RAD2D_p':tri2d_ffi.new("float[]",list(RAD2D.flatten())),
                                 'THETA2D_p':tri2d_ffi.new("float[]",list(THETA2D.flatten())),
                                 'n1_p':tri2d_ffi.new("int[]",[nr_grid]),
                                 'n2_p':tri2d_ffi.new("int[]",[nth_grid]),
                                 'nphi_p':tri2d_ffi.new("int[]",[nphi_grid]),
                                 'ntri_p':tri2d_ffi.new("int[]",[(nth_grid-1)*(nr_grid-1)*2]),
                                 'tri_vertices_p':tri2d_ffi.new("int[]",list(tri_vertices.flatten())),
                                 'tri_neighbors_p':tri2d_ffi.new("int[]",list(tri_neighbors.flatten())),
                                 'tri_ij_p':tri2d_ffi.new("int[]",list(tri_ij.flatten())),
                                 'xmin_sg_p':tri2d_ffi.new("float[]",[xmin_sg]),
                                 'dx_sg_p':tri2d_ffi.new("float[]",[dx_sg]),
                                 'nx_sg_p':tri2d_ffi.new("int[]",[nx_sg]),
                                 'rmin_sg_p':tri2d_ffi.new("float[]",[rmin_sg]),
                                 'dr_sg_p':tri2d_ffi.new("float[]",[dr_sg]),
                                 'nr_sg_p':tri2d_ffi.new("int[]",[nr_sg]),
                                 'start_index_sg_p':tri2d_ffi.new("int[]",list(start_index_sg.flatten())),
                                 }
            

            net_indices_c = G.GridMapping(data_files[0:self.num_files_per_set],False,'D')
            X2D_c,R2D_c,PHI_c = G.PrepareCylindricalCellCenteredGrid(X2D,R2D,PHI)
            RAD2D_c = sqrt(X2D_c*X2D_c+R2D_c*R2D_c)
            THETA2D_c = arccos(X2D_c/RAD2D_c)
            nth_grid_c,nr_grid_c = X2D_c.shape
            nphi_grid_c = len(PHI_c)

            X3D_c = zeros((nphi_grid_c,nth_grid_c,nr_grid_c),dtype=float32)

            for iphi in range(nphi_grid_c):
                X3D_c[iphi,:,:] = X2D_c
                      
            tri_vertices_c,tri_neighbors_c,tri_ij_c = G.PrepareTriangularInterpolation(X2D_c.shape[1],X2D_c.shape[0])
            # 2024/01/16 - removing overlay grid
            #xmin_sg_c,dx_sg_c,nx_sg_c,rmin_sg_c,dr_sg_c,nr_sg_c,start_index_sg_c = G.PrepareOverlayGrid2D(X2D_c,R2D_c,tri_vertices_c)
            xmin_sg_c,dx_sg_c,nx_sg_c,rmin_sg_c,dr_sg_c,nr_sg_c,start_index_sg_c = G.PrepareOverlaySphGrid2D(X2D_c,R2D_c,tri_vertices_c)
            self._gridcenters = {'X':X3D_c, # data that matches shape of grid in 3D
                                 'net_indices':net_indices_c,
                                 'X2D':X2D_c,
                                 'R2D':R2D_c,
                                 'RAD2D':RAD2D_c, # radial distances of grid vertices  from orign (Earth center) 
                                 'PHI':PHI_c,
                                 'THETA2D': THETA2D_c,
                                 'n1':nr_grid_c,
                                 'n2':nth_grid_c,
                                 'nphi':nphi_grid_c,
                                 'ntri':(nth_grid_c-1)*(nr_grid_c-1)*2,
                                 'tri_vertices':tri_vertices_c,
                                 'tri_neighbors':tri_neighbors_c,
                                 'tri_ij':tri_ij_c,
                                 # overlay grid has been removed
                                 'xmin_sg':xmin_sg_c,
                                 'dx_sg':dx_sg_c,
                                 'nx_sg':nx_sg_c,
                                 'rmin_sg':rmin_sg_c,
                                 'dr_sg':dr_sg_c,
                                 'nr_sg':nr_sg_c,
                                 'start_index_sg':start_index_sg_c,
                                 # pointers
                                 'X2D_p':tri2d_ffi.new("float[]",list(X2D_c.flatten())),
                                 'R2D_p':tri2d_ffi.new("float[]",list(R2D_c.flatten())),
                                 'PHI_p':tri2d_ffi.new("float[]",list(PHI_c.flatten())),
                                 'RAD2D_p':tri2d_ffi.new("float[]",list(RAD2D_c.flatten())),
                                 'THETA2D_p':tri2d_ffi.new("float[]",list(THETA2D_c.flatten())),
                                 'n1_p':tri2d_ffi.new("int[]",[nr_grid_c]),
                                 'n2_p':tri2d_ffi.new("int[]",[nth_grid_c]),
                                 'nphi_p':tri2d_ffi.new("int[]",[nphi_grid_c]),
                                 'ntri_p':tri2d_ffi.new("int[]",[(nth_grid_c-1)*(nr_grid_c-1)*2]),
                                 'tri_vertices_p':tri2d_ffi.new("int[]",list(tri_vertices_c.flatten())),
                                 'tri_neighbors_p':tri2d_ffi.new("int[]",list(tri_neighbors_c.flatten())),
                                 'tri_ij_p':tri2d_ffi.new("int[]",list(tri_ij_c.flatten())),
                                 'xmin_sg_p':tri2d_ffi.new("float[]",[xmin_sg_c]),
                                 'dx_sg_p':tri2d_ffi.new("float[]",[dx_sg_c]),
                                 'nx_sg_p':tri2d_ffi.new("int[]",[nx_sg_c]),
                                 'rmin_sg_p':tri2d_ffi.new("float[]",[rmin_sg_c]),
                                 'dr_sg_p':tri2d_ffi.new("float[]",[dr_sg_c]),
                                 'nr_sg_p':tri2d_ffi.new("int[]",[nr_sg_c]),
                                 'start_index_sg_p':tri2d_ffi.new("int[]",list(start_index_sg_c.flatten())),
                                 }
                
            # GAMERA data has all the timestamps in all the files in each set)
            # This will confuse the lazy interpolation, so stripping down the
            # times dict to one file entry (only start and end fields needed).
            # For ms_timing, 'all' has ms resolution, start and end do not.
            self.singletimes = {'start': [self.times[p]['all'][0]],
                                'end': [self.times[p]['all'][-1]]}
            
            # register interpolators for each requested variable
            t_reg = perf_counter()
            # store original list b/c gridded interpolators change key list
            varname_list = [key for key in self.variables.keys()]
            for varname in varname_list:
                self.register_variables(varname, gridded_int)
            if verbose:
                print(f'Took {perf_counter()-t_reg:.5f}s to register ' +
                      f'{len(varname_list)} variables.')
            if verbose:
                print(f'Took a total of {perf_counter()-t0:.5f}s to kamodofy' +
                      f' {len(varname_list)} variables.')

        # define and register the variable
        def register_variables(self, varname, gridded_int):
            '''Functionalizes the indicated dataset.'''

            # determine coordinate variables and xvec by coord list
            key = self.variables[varname]['data']
            coord_list = [value[5] for key, value in model_varnames.items()
                          if value[0] == varname][0]
            gvar = [key for key, value in model_varnames.items() if
                    value[0] == varname][0]
            coord_str = [value[3]+value[4] for key, value in
                         model_varnames.items() if value[0] == varname][0]
            file_dir = self.pattern_files[key][0].split(
                basename(self.pattern_files[key][0]))[0]

            if 'time' in coord_list:
                coord_dict = {'time': {'units': 'hr',
                                       'data': self.times[key]['all']},
                              'X': {'units': 'R_E', 'data': self._X},
                              'Y': {'units': 'R_E', 'data': self._Y},
                              'Z': {'units': 'R_E', 'data': self._Z}}
                #tri2d_lib.setup_tri_pointers(self._gridcenters['X2D'],
                #                               self._gridcenters['R2D'],
                #                               self._gridcenters['nx'],
                #                               self._gridcenters['nr'],
                #                               self._gridcenters['ntri'],
                #                               self._gridcenters['tri_vertices'],
                #                               self._gridcenters['tri_ij'],
                #                               self._gridcenters['tri_neighbors'],
                #                               self._gridcenters['xmin_sg'],
                #                               self._gridcenters['rmin_sg'],
                #                               self._gridcenters['dx_sg'],
                #                               self._gridcenters['dr_sg'],
                #                               self._gridcenters['nx_sg'],
                #                               self._gridcenters['nr_sg'],
                #                               self._gridcenters['start_index_sg'],
                #                               self._gridcenters['PHI'],
                #                               self._gridcenters['nphi']
                #                               )

                tri2d_lib.setup_tri_pointers(self._gridcorners['X2D_p'],
                                             self._gridcorners['R2D_p'],
                                             self._gridcorners['n1_p'],
                                             self._gridcorners['n2_p'],
                                             self._gridcorners['ntri_p'],
                                             self._gridcorners['tri_vertices_p'],
                                             self._gridcorners['tri_ij_p'],
                                             self._gridcorners['tri_neighbors_p'],
                                             self._gridcorners['xmin_sg_p'],
                                             self._gridcorners['rmin_sg_p'],
                                             self._gridcorners['dx_sg_p'],
                                             self._gridcorners['dr_sg_p'],
                                             self._gridcorners['nx_sg_p'],
                                             self._gridcorners['nr_sg_p'],
                                             self._gridcorners['start_index_sg_p'],
                                             self._gridcorners['PHI_p'],
                                             self._gridcorners['nphi_p']
                                             )

                def func(i, fi):  # i = file# (always 0), fi = slice# (= Step#)
                    # interpolate cell-centered data
                    # the interpolator will find the triangle and its location in
                    # i, j, and k then return the center value for that cell (i,j,k)
                    data_files = (self.pattern_files[key])[0]
                    data_files = RU.glob('_'.join(data_files.split('_')[0:-3])+'*.h5')
                    data = G.AssembleGrid(data_files,self._gridcenters['net_indices'],False,gvar,step=fi)
#                    data = G.AssembleGrid(self.pattern_files[key],self._gridcenters['net_indices'],False,gvar,step=fi)
                    # define custom interpolator here  *******************************
                    def interp(xvec):
                        tic = perf_counter()
                        X, Y, Z = xvec.T  # xvec can be used like this
                        if not isinstance(X, ndarray):
                            X = array([X],dtype=float32)
                            Y = array([Y],dtype=float32)
                            Z = array([Z],dtype=float32)
                        else:
                            X = array(list(X),dtype=float32)
                            Y = array(list(Y),dtype=float32)
                            Z = array(list(Z),dtype=float32)

                        # call custom interpolator here
                        X_ptr = tri2d_ffi.new("float[]",list(X))
                        R = sqrt(Y*Y+Z*Z)
                        R_ptr = tri2d_ffi.new("float[]",list(R))
                        PHI = mod(arctan2(Z,Y)+2*pi,2*pi) # have to test whether we cover all angles
                        PHI_ptr = tri2d_ffi.new("float[]",list(PHI))

                        return_data = zeros(len(X),dtype=float32)
                        return_data[:] = NaN
                        #return_data_ffi = tri2d_ffi.new("float[]",list(return_data))
                        return_data_ffi = tri2d_ffi.cast("float *",return_data.ctypes.data)
                        data_pointer = tri2d_ffi.cast("float *",data.ctypes.data)
                        is_cell_centered = 1
                        # zero-order interpolation returns value at cell center
                        failure = tri2d_lib.interpolate_tri2d_plus_1d_multipos(
                            X_ptr,R_ptr,PHI_ptr,len(X),data_pointer,return_data_ffi,is_cell_centered) 
                        
                        toc = perf_counter()
                        # print('interpolation done in ',toc-tic,' seconds')

                        return return_data

                    return interp

                # functionalize the variable dataset
                tmp = self.variables[varname]
                tmp['data'] = zeros((2, 2, 2, 2))  # saves execution time
                self = RU.Functionalize_Dataset(
                    self, coord_dict, varname, tmp, gridded_int, coord_str,
                    interp_flag=3, func=func, times_dict=self.singletimes,
                    func_default='custom')
            else:
                # functionalize the constants
                coord_dict = {'X': {'units': 'R_E', 'data': self._X},
                              'Y': {'units': 'R_E', 'data': self._Y},
                              'Z': {'units': 'R_E', 'data': self._Z}}
                tri2d_lib.setup_tri_pointers(self._gridcorners['X2D_p'],
                                             self._gridcorners['R2D_p'],
                                             self._gridcorners['n1_p'],
                                             self._gridcorners['n2_p'],
                                             self._gridcorners['ntri_p'],
                                             self._gridcorners['tri_vertices_p'],
                                             self._gridcorners['tri_ij_p'],
                                             self._gridcorners['tri_neighbors_p'],
                                             self._gridcorners['xmin_sg_p'],
                                             self._gridcorners['rmin_sg_p'],
                                             self._gridcorners['dx_sg_p'],
                                             self._gridcorners['dr_sg_p'],
                                             self._gridcorners['nx_sg_p'],
                                             self._gridcorners['nr_sg_p'],
                                             self._gridcorners['start_index_sg_p'],
                                             self._gridcorners['PHI_p'],
                                             self._gridcorners['nphi_p']
                                             )

                def func_const():
                    # interpolate data on cell corners.
                    # here X,Y,Z only 
                    data_files = (self.pattern_files[key])[0]
                    data_files = RU.glob('_'.join(data_files.split('_')[0:-3])+'*.h5')

                    if gvar == 'dV':  # dV values are cell-centered!
                        data = G.AssembleGrid(data_files,self._gridcenters['net_indices'],True,gvar)
                        is_cell_centered = 1
                    else:
                        data = G.AssembleGrid(data_files,self._gridcorners['net_indices'],True,gvar)
                        is_cell_centered = 0
                    
                    # define custom interpolator here  ***********************************
                    def interp(xvec):
                        tic = perf_counter()
                        X, Y, Z = array(xvec).T  # xvec can be used like this
                        if not isinstance(X, ndarray):
                            X = array([X],dtype=float32)
                            Y = array([Y],dtype=float32)
                            Z = array([Z],dtype=float32)
                        else:
                            X = array(list(X),dtype=float32)
                            Y = array(list(Y),dtype=float32)
                            Z = array(list(Z),dtype=float32)


                        X_ptr = tri2d_ffi.new("float[]",list(X))
                        R = sqrt(Y*Y+Z*Z)
                        R_ptr = tri2d_ffi.new("float[]",list(R))
                        PHI = mod(arctan2(Z,Y)+2*pi,2*pi) 
                        PHI_ptr = tri2d_ffi.new("float[]",list(PHI))
                        
                        return_data = zeros(len(X),dtype=float32)
                        return_data[:] = NaN
                        return_data_ffi = tri2d_ffi.cast("float *",return_data.ctypes.data)
                        data_pointer = tri2d_ffi.cast("float *",data.ctypes.data)

                        # first-order interpolation on triangles
                        failure = tri2d_lib.interpolate_tri2d_plus_1d_multipos(
                            X_ptr,R_ptr,PHI_ptr,len(X),data_pointer,return_data_ffi,is_cell_centered) 
                        toc = perf_counter()

                        #print('interpolation done in ',toc-tic,' seconds')

                        return return_data

                    return interp

                # functionalize the variable dataset
                tmp = self.variables[varname]
                tmp['data'] = zeros((2, 2, 2))  # saves execution time

                self = RU.Functionalize_Dataset(
                    self, coord_dict, varname, tmp, gridded_int, coord_str,
                    interp_flag=0, func=func_const, func_default='custom')
            return
    return MODEL
