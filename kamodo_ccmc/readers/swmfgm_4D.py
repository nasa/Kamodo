'''
Written by Rebecca Ringuette, 2021
'''
from datetime import datetime, timezone

# variable name in file: [standardized variable name, descriptive term, units]
model_varnames = {'b1x': ['B_1x', 'x component of the deviation from the ' +
                          'Earth dipole field', 0, 'GSM', 'car',
                          ['time', 'X', 'Y', 'Z'], 'nT'],
                  'b1y': ['B_1y', 'x component of the deviation from the ' +
                          'Earth dipole field', 0, 'GSM', 'car',
                          ['time', 'X', 'Y', 'Z'], 'nT'],
                  'b1z': ['B_1z', 'x component of the deviation from the ' +
                          'Earth dipole field', 0, 'GSM', 'car',
                          ['time', 'X', 'Y', 'Z'], 'nT'],
                  'bx': ['B_x', 'x component of the magnetic field', 0, 'GSM',
                         'car', ['time', 'X', 'Y', 'Z'], 'nT'],
                  'by': ['B_y', 'y component of the magnetic field', 0, 'GSM',
                         'car', ['time', 'X', 'Y', 'Z'], 'nT'],
                  'bz': ['B_z', 'z component of the magnetic field', 0, 'GSM',
                         'car', ['time', 'X', 'Y', 'Z'], 'nT'],
                  'e': ['u', 'Energy density', 0, 'GSM', 'car',
                        ['time', 'X', 'Y', 'Z'], 'J/m**3'],
                  'jx': ['J_x', 'x component of the current density', 0, 'GSM',
                         'car', ['time', 'X', 'Y', 'Z'], 'muA/m**2'],
                  'jy': ['J_y', 'y component of the current density', 0, 'GSM',
                         'car', ['time', 'X', 'Y', 'Z'], 'muA/m**2'],
                  'jz': ['J_z', 'z component of the current density', 0, 'GSM',
                         'car', ['time', 'X', 'Y', 'Z'], 'muA/m**2'],
                  'p': ['P', 'Pressure', 0, 'GSM', 'car',
                        ['time', 'X', 'Y', 'Z'], 'nPa'],
                  'rho': ['N_p', 'proton number density', 0, 'GSM', 'car',
                          ['time', 'X', 'Y', 'Z'], '10**6/cm**3'],  # Mp/cc
                  'ux': ['v_x', 'x component of velocity', 0, 'GSM', 'car',
                         ['time', 'X', 'Y', 'Z'], 'km/s'],
                  'uy': ['v_y', 'y component of velocity', 0, 'GSM', 'car',
                         ['time', 'X', 'Y', 'Z'], 'km/s'],
                  'uz': ['v_z', 'z component of velocity', 0, 'GSM', 'car',
                         ['time', 'X', 'Y', 'Z'], 'km/s'],
                  'th': ['theta_Btilt', 'Dipole tilt angle: positive for the' +
                         ' northern hemisphere dipole axis tilting towards ' +
                         'the Sun in northern hemisphere summer', 0, 'GSM',
                         'car', ['time'], 'radians']}


# coordinates stored as x, y, z in earth radii.
def MODEL():

    from kamodo import Kamodo
    import spacepy.pybats as sp
    from glob import glob
    import time
    from os.path import basename, isfile
    from numpy import array, NaN, unique, zeros, ravel, ndarray
    from numpy import min, max, floor, log
    from time import perf_counter
    import kamodo_ccmc.readers.reader_utilities as RU

    from readers.OCTREE_BLOCK_GRID._interpolate_amrdata import ffi
    from readers.OCTREE_BLOCK_GRID._interpolate_amrdata import lib

    class MODEL(Kamodo):
        '''SWMF model data reader for magnetosphere outputs.

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
            The run constants c, g, and R are saved in the kamodo.constants
            dictionary. R is the radius of the minimum inner boundary in R_E,
            c is the reduced speed of light, and g is the ratio of specific
            heats. Consult model documentation for more information.

        Notes and instructions:
        - The SWMF global magnetosphere outputs are given in one or more files
          per timestep in *.out files. The structure of these files are
          advantageous for efficient interpolation, so no file conversion is
          attempted.
        - The files are small and contain one time step per file, so
          interpolation method 1 is chosen. A custom interpolator is required.
        - The reader requires code in readers/OCTREE_BLOCK_GRID/:
          Compile the library using the command
            python interpolate_amrdata_extension_build.py
          inside the anaconda/miniconda3/miniforge3 environment for Kamodo.
          The shared library file name contains the python version and the name
          of the operation system, e.g.,
            _interpolate_amrdata.cpython-39-darwin.so
        - There is a pull request pending to implement sort_unstructured_data
          in the SpacePy software package. Until that PR is merged into
          SpacePy, you will need to manually modify the package. The
          lib/site-packages/spacepy/pybats/__init__.py script needs to be
          modified in your installation to not sort unstructured data in
          _read_Idl_bin() by default. See detailed instructions on the 2022 AGU
          poster https://doi.org/10.22541/essoar.167214301.16153548/v1 , at the
          bottom of the second to last column.
        '''
        def __init__(self, file_dir, variables_requested=[],
                     filetime=False, verbose=False, gridded_int=True,
                     printfiles=False, **kwargs):
            super(MODEL, self).__init__(**kwargs)
            self.modelname = 'SWMF_GM'
            t0 = perf_counter()

            # first, check for file list, create if DNE
            list_file = file_dir + self.modelname + '_list.txt'
            time_file = file_dir + self.modelname + '_times.txt'
            self.times, self.pattern_files = {}, {}
            if not isfile(list_file) or not isfile(time_file):
                files = sorted(glob(file_dir+'*.out')) +\
                    sorted(glob(file_dir+'*.outs'))  # .out or .outs
                self.filename = ''.join([f+',' for f in files])[:-1]
                patterns = unique([basename(f)[:10] for f in files])
                # get time grid from files
                dt = sp.IdlFile(files[0]).attrs['time']
                if dt is not None:  # filedate given not always at midnight
                    self.filedate = datetime.strptime(
                        dt.isoformat()[:10], '%Y-%m-%d').replace(
                        tzinfo=timezone.utc)
                else:
                    self.filedate = datetime(1900, 1, 1, tzinfo=timezone.utc)

                # establish time attributes from filenames
                for p in patterns:
                    # get list of files to loop through later
                    pattern_files = sorted(glob(file_dir+p+'*.out')) +\
                        sorted(glob(file_dir+p+'*.outs'))
                    self.pattern_files[p] = pattern_files
                    self.times[p] = {'start': [], 'end': [], 'all': []}

                    # loop through to get times, one file per time step
                    self.times[p]['start'] = array([RU.tstr_to_hrs(
                        sp.IdlFile(f).attrs['strtime'][:-1].replace(
                            'h', ':').replace('m', ':'))
                        for f in pattern_files])
                    self.times[p]['end'] = self.times[p]['start'].copy()
                    self.times[p]['all'] = self.times[p]['start'].copy()

                # create time list file if DNE
                RU.create_timelist(list_file, time_file, self.modelname,
                                   self.times, self.pattern_files,
                                   self.filedate)
            else:  # read in data and time grids from file list
                self.times, self.pattern_files, self.filedate, self.filename =\
                    RU.read_timelist(time_file, list_file)
            if filetime:
                return  # return times as is to prevent infinite recursion
            # only one pattern, so simplifying code
            p = list(self.pattern_files.keys())[0]

            # perform initial check on variables_requested list
            if len(variables_requested) > 0 and variables_requested != 'all':
                test_list = [value[0] for key, value in model_varnames.items()]
                err_list = [item for item in variables_requested if item not in
                            test_list]
                if len(err_list) > 0:
                    print('Variable name(s) not recognized:', err_list)
                for item in err_list:
                    variables_requested.remove(item)
                if len(variables_requested) == 0:
                    return

            # collect variable list
            mhd = sp.IdlFile(self.pattern_files[p][0])
            # mhd = sp.IdlFile(self.pattern_files[p][0],
            #                  sort_unstructured_data=False)
            file_variables = [key for key in mhd.keys() if key not in
                              ['grid', 'status', 'x', 'y', 'z']] + ['th']
            if len(variables_requested) > 0 and variables_requested != 'all':
                gvar_list = [key for key, value in model_varnames.items()
                             if value[0] in variables_requested and
                             key in file_variables]

                # check for variables requested but not available
                if len(gvar_list) != len(variables_requested):
                    err_list = [value[0] for key, value in
                                model_varnames.items() if value[0] in
                                variables_requested and key not in gvar_list]
                    if len(err_list) > 0:
                        print('Some requested variables are not available:',
                              err_list)
            else:  # only input variables on the avoid_list if requested
                gvar_list = [key for key in file_variables
                             if key in model_varnames.keys()]
                # returns list of variables included in data files
                if variables_requested == 'all':
                    self.var_dict = {value[0]: value[1:] for key, value in
                                     model_varnames.items() if key in
                                     gvar_list}
                    return

            # retrieve coordinate grids and run constants
            self._X = array(unique(mhd['x']))
            self._Y = array(unique(mhd['y']))
            self._Z = array(unique(mhd['z']))
            self.constants = {'c': mhd.meta['c'], 'g': mhd.meta['g']}
            # older runs may not have these
            if 'R' in mhd.meta.keys():
                self.constants['R'] = mhd.meta['R']
            if 'NX' in mhd.meta.keys():
                self.constants['NX'] = mhd.meta['NX']
            if 'NY' in mhd.meta.keys():
                self.constants['NY'] = mhd.meta['NY']
            if 'NZ' in mhd.meta.keys():
                self.constants['NZ'] = mhd.meta['NZ']

            # store kay for each variable desired
            self.variables = {model_varnames[var][0]: {
                'units': model_varnames[var][-1],
                'data': p} for var in gvar_list}

            # store a few items
            self.missing_value = NaN
            if verbose:
                print(f'Took {perf_counter()-t0:.6f}s to read in data')
            if printfiles:
                print(self.filename)
            self.octree = {p: {} for p in self.pattern_files.keys()}

            # register interpolators for each requested variable
            t_reg = perf_counter()
            # store original list b/c gridded interpolators change key list
            varname_list = [key for key in self.variables.keys()]
            for varname in varname_list:
                self.register_variables(varname, gridded_int, verbose=verbose)
            if verbose:
                print(f'Took {perf_counter()-t_reg:.5f}s to register ' +
                      f'{len(varname_list)} variables.')
            if verbose:
                print(f'Took a total of {perf_counter()-t0:.5f}s to kamodofy' +
                      f' {len(varname_list)} variables.')

        def setup_octree(self, x, y, z, verbose=False):
            '''
            This function requires _interpolate_amrdata*.so in
            readers/OCTREE_BLOCK_GRID/
            '''
            if verbose:
                tic = time.time()

            Ncell = int(len(x))

            if verbose:
                print("setup_octree: Ncell=", Ncell)

            if 'NX' in self.constants.keys():
                NX = int(self.constants['NX'])
            else:
                NX = 4  # minimum block size in BATSRUS is 4x4xc4 cells
            if 'NZ' in self.constants.keys():
                NY = int(self.constants['NY'])
            else:
                NY = 4
            if 'NZ' in self.constants.keys():
                NZ = int(self.constants['NZ'])
            else:
                NZ = 4

            N_blks = int(len(x)/(NX*NY*NZ))

            if 'R' in self.constants.keys():
                R = float(self.constants['R'])
            else:
                R = 3.0

            # initialize cffi arrays with Python lists may be slightly faster
            octree = {}
            x = ffi.new("float[]", list(x))
            y = ffi.new("float[]", list(y))
            z = ffi.new("float[]", list(z))
            # 2 elements per block, smallest blocks have 2x2x2=8 positions
            # at cell corners (GUMICS), so maximum array size is N/4
            # BATSRUS blocks have at least 4x4x4 cell centers
            octree['xr_blk'] = ffi.new("float[]", 2*N_blks)
            octree['yr_blk'] = ffi.new("float[]", 2*N_blks)
            octree['zr_blk'] = ffi.new("float[]", 2*N_blks)
            octree['x_blk'] = ffi.new("float[]", NX*N_blks)
            octree['y_blk'] = ffi.new("float[]", NY*N_blks)
            octree['z_blk'] = ffi.new("float[]", NZ*N_blks)
            octree['box_range'] = ffi.new("float[]", 7)
            octree['NX'] = ffi.new("int[]", [NX])
            octree['NY'] = ffi.new("int[]", [NY])
            octree['NZ'] = ffi.new("int[]", [NZ])

            N_blks = lib.xyz_ranges(Ncell, x, y, z, octree['xr_blk'],
                                    octree['yr_blk'], octree['zr_blk'],
                                    octree['x_blk'], octree['y_blk'],
                                    octree['z_blk'], octree['box_range'],
                                    octree['NX'], octree['NY'],
                                    octree['NZ'], 1)

            if N_blks < 0:
                print('NX:', list(octree['NX']))
                print('NY:', list(octree['NY']))
                print('NZ:', list(octree['NZ']))
                print('Y:', list(y[0:16]))
                print('Z:', list(z[0:16]))
                print('X_blk:', list(octree['x_blk'][0:16]))
                print('Y_blk:', list(octree['y_blk'][0:16]))
                print('Z_blk:', list(octree['z_blk'][0:16]))
                raise IOError("Block shape and size was not determined.")

            octree['N_blks'] = ffi.new("int[]", [N_blks])

            N_octree = int(N_blks*8/7)
            octree['octree_blocklist'] = ffi.new("octree_block[]", N_octree)

            dx_blk = zeros(N_blks, dtype=float)
            dy_blk = zeros(N_blks, dtype=float)
            dz_blk = zeros(N_blks, dtype=float)

            for i in range(0, N_blks):
                dx_blk[i] = octree['xr_blk'][2*i+1]-octree['xr_blk'][2*i]
                dy_blk[i] = octree['yr_blk'][2*i+1]-octree['yr_blk'][2*i]
                dz_blk[i] = octree['zr_blk'][2*i+1]-octree['zr_blk'][2*i]

            dxmin_blk = min(dx_blk)
            dxmax_blk = max(dx_blk)
            dymin_blk = min(dy_blk)
            dymax_blk = max(dy_blk)
            dzmin_blk = min(dz_blk)
            dzmax_blk = max(dz_blk)
            XMIN = octree['box_range'][0]
            XMAX = octree['box_range'][1]
            YMIN = octree['box_range'][2]
            YMAX = octree['box_range'][3]
            ZMIN = octree['box_range'][4]
            ZMAX = octree['box_range'][5]
            p1 = int(floor((XMAX-XMIN)/dxmax_blk+0.5))
            p2 = int(floor((YMAX-YMIN)/dymax_blk+0.5))
            p3 = int(floor((ZMAX-ZMIN)/dzmax_blk+0.5))
            while ((int(p1/2)*2 == p1)
                   & (int(p2/2)*2 == p2)
                   & (int(p3/2)*2 == p3)
                   ):
                p1 = int(p1/2)
                p2 = int(p2/2)
                p3 = int(p3/2)

            if verbose:
                print("XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX: ",
                      XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX)
                print("p1, p2, p3: ", p1, p2, p3, "dxmin,dymin,dzmin: ",
                      dxmin_blk, dymin_blk, dzmin_blk)

            MAX_AMRLEVEL = int(log((XMAX-XMIN)/(p1*dxmin_blk))/log(2.)+0.5)+1
            octree['MAX_AMRLEVEL'] = MAX_AMRLEVEL

            if verbose:
                print("MAX_AMRLEVEL:", MAX_AMRLEVEL)

            octree['numparents_at_AMRlevel'] = ffi.new("int[]",
                                                       N_blks*(MAX_AMRLEVEL+1))
            octree['block_at_AMRlevel'] = ffi.new("int[]",
                                                  N_blks*(MAX_AMRLEVEL+1))
            success = int(-1)
            success = lib.setup_octree(N_blks,
                                       octree['xr_blk'], octree['yr_blk'],
                                       octree['zr_blk'], MAX_AMRLEVEL,
                                       octree['box_range'],
                                       octree['octree_blocklist'], N_octree,
                                       octree['numparents_at_AMRlevel'],
                                       octree['block_at_AMRlevel'])

            if verbose:
                toc = time.time()
                print('setup_octree: ', toc-tic, 'seconds elapsed')

            return (octree)

        def register_variables(self, varname, gridded_int, verbose=False):
            '''Functionalizes the indicated dataset.'''

            # determine coordinate variables and xvec by coord list
            key = self.variables[varname]['data']
            coord_list = [value[5] for key, value in model_varnames.items()
                          if value[0] == varname][0]
            coord_dict = {'time': {'units': 'hr',
                                   'data': self.times[key]['all']}}
            gvar = [key for key, value in model_varnames.items() if
                    value[0] == varname][0]
            coord_str = [value[3]+value[4] for key, value in
                         model_varnames.items() if value[0] == varname][0]

            if len(coord_list) == 1:  # read all the values in
                self.variables[varname]['data'] = ravel([
                    sp.IdlFile(f).meta[gvar] for f in self.pattern_files[key]])
                self = RU.Functionalize_Dataset(
                    self, coord_dict, varname, self.variables[varname],
                    gridded_int, coord_str, interp_flag=0)
                return

            # remainder of logic is for time + 3D variables
            coord_dict['X'] = {'units': 'R_E', 'data': self._X}
            coord_dict['Y'] = {'units': 'R_E', 'data': self._Y}
            coord_dict['Z'] = {'units': 'R_E', 'data': self._Z}

            # all of data is time + 3D spatial
            def func(i):  # i is the file/time number
                # initialize octree instance for file if not already done
                if i not in self.octree[key].keys():
                    # files is closed after read in spacepy.pybats.IdlFile
                    mhd = sp.IdlFile(self.pattern_files[key][i])
                    # pull request pending to implement sort_unstructured_data
                    # lib/site-packages/spacepy/pybats/__init__.py needs to be
                    # modified in your installation to not sort unstructured
                    # data in _read_Idl_bin() by default
                    #  def _read_idl_bin(pbdat, header='units', start_loc=0,
                    #                     keep_case=True, headeronly=False,
                    #                     sort_unstructured_data=False):
                    # ...
                    #
                    # # Unstructured data can be in any order, so let's sort it
                    # if gtyp == 'Unstructured' and sort_unstructured_data:
                    #     gridtotal = np.zeros(npts)
                    #
                    # Once pull request is implemented in spacepy use keyword:
                    # mhd = sp.IdlFile(self.pattern_files[key][i],
                    #                  sort_unstructured_data=False)
                    x = array(mhd['x'])
                    y = array(mhd['y'])
                    z = array(mhd['z'])
                    # initialize octree object and store in dictionary
                    self.octree[key][i] = self.setup_octree(x, y, z,
                                                            verbose=verbose)

                # update which octree arrays the library points to
                lib.setup_octree_pointers(self.octree[key][i]['MAX_AMRLEVEL'],
                                          self.octree[key][i]['octree_blocklist'],
                                          self.octree[key][i]['numparents_at_AMRlevel'],
                                          self.octree[key][i]['block_at_AMRlevel'])
                # get data from file
                data = sp.IdlFile(self.pattern_files[key][i])[gvar]  # dmarray
                var_data = ffi.new("float[]", list(data))
                
                # assign custom interpolator: Lutz Rastaetter 2021
                def interp(xvec):
                    tic = perf_counter()
                    if not isinstance(xvec, ndarray):
                        xvec = array(xvec)

                    X, Y, Z = xvec.T  # xvec can be used like this
                    if not isinstance(X, ndarray):
                        X = array([X])
                        Y = array([Y])
                        Z = array([Z])

                    xpos = ffi.new("float[]", list(X))
                    ypos = ffi.new("float[]", list(Y))
                    zpos = ffi.new("float[]", list(Z))
                    npos = len(X)
                    return_data = list(zeros(npos, dtype=float))
                    return_data_ffi = ffi.new("float[]", return_data)

                    IS_ERROR = lib.interpolate_amrdata_multipos(
                        xpos, ypos, zpos, npos, var_data, return_data_ffi)
                    if IS_ERROR != 0:
                        print("Warning: SWMF/BATSRUS interpolation failed.")

                    return_data[:] = list(return_data_ffi)
                    toc = perf_counter()
                    if verbose:
                        print(varname, ' interp: ', toc-tic, 'seconds elapsed')

                    return return_data

                return interp

            # functionalize the variable dataset
            tmp = self.variables[varname]
            tmp['data'] = zeros((2, 2, 2, 2))  # saves execution time
            self = RU.Functionalize_Dataset(
                self, coord_dict, varname, tmp, gridded_int, coord_str,
                interp_flag=1, func=func, func_default='custom')
            return
    return MODEL
