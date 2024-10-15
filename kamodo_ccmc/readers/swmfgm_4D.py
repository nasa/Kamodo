'''
Written by Rebecca Ringuette, 2021
Added custom interpolator by Lutz Rastaetter, 2022
'''
from datetime import datetime, timezone
import numpy as np

# variable name in file: [standardized variable name, descriptive term, units]
model_varnames = {'b1x': ['B1_x', 'x component of the deviation from the ' +
                          'Earth dipole magnetic field', 0, 'GSM', 'car',
                          ['time', 'X', 'Y', 'Z'], 'nT'],
                  'b1y': ['B1_y', 'y component of the deviation from the ' +
                          'Earth dipole magnetic field', 0, 'GSM', 'car',
                          ['time', 'X', 'Y', 'Z'], 'nT'],
                  'b1z': ['B1_z', 'z component of the deviation from the ' +
                          'Earth dipole magnetic field', 0, 'GSM', 'car',
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
                          ['time', 'X', 'Y', 'Z'], '1/cm**3'],  # Mp/cc
                  'ux': ['v_x', 'x component of velocity', 0, 'GSM', 'car',
                         ['time', 'X', 'Y', 'Z'], 'km/s'],
                  'uy': ['v_y', 'y component of velocity', 0, 'GSM', 'car',
                         ['time', 'X', 'Y', 'Z'], 'km/s'],
                  'uz': ['v_z', 'z component of velocity', 0, 'GSM', 'car',
                         ['time', 'X', 'Y', 'Z'], 'km/s'],
                  'th': ['theta_Btilt', 'Dipole tilt angle: positive for the' +
                         ' northern hemisphere dipole axis tilting towards ' +
                         'the Sun in northern hemisphere summer', 0, 'GSM',
                         'car', ['time'], 'radians'],
                  'status': ['status',
                         'Classification of magnetic fieldlines', 0, 'GSM', 'car',
                         ['time', 'X', 'Y', 'Z'], ''],
                  'theta1': ['BfootNlat',
                         'Magnetic fieldline footpoint North (latitude)', 0, 'GSM', 'car',
                         ['time', 'X', 'Y', 'Z'], 'deg'],
                  'theta2': ['BfootSlat',
                         'Magnetic fieldline footpoint South (latitude)', 0, 'GSM', 'car',
                         ['time', 'X', 'Y', 'Z'], 'deg'],
                  'phi1': ['BfootNlon',
                         'Magnetic fieldline footpoint North (longitude)', 0, 'GSM', 'car',
                         ['time', 'X', 'Y', 'Z'], 'deg'],
                  'phi2': ['BfootSlon',
                         'Magnetic fieldline footpoint South (longitude)', 0, 'GSM', 'car',
                         ['time', 'X', 'Y', 'Z'], 'deg']}

# Alternate variable names in model output, dimensioned and dimensionless
alt_varnames = {'N_p': ['rho', 'Rho'],
                'v_x': ['ux', 'Ux'],
                'v_y': ['uy', 'Uy'],
                'v_z': ['uz', 'Uz'],
                'P': ['p', 'P'],
                'B1_x': ['b1x', 'B1x'],
                'B1_y': ['b1y', 'B1y'],
                'B1_z': ['b1z', 'B1z'],
                'B_x': ['bx', 'Bx'],
                'B_y': ['by', 'By'],
                'B_z': ['bz', 'Bz']}

# Dimensionalizing constants and scales to apply to variables
kbSI = 1.3807e-23       # Boltzmann constant
mpSI = 1.6726e-27       # proton mass
mu0SI = 4.*np.pi*1.e-7  # vacuum permeability
eSI = 1.602e-19         # elementary charge
ReSI = 6378.e+3         # radius of Earth
RsSI = 6.96e+8          # radius of Sun
AuSI = 1.4959787e+11    # astronomical unit
cSI = 2.9979e+8         # speed of light
xSI = ReSI
rhoSI = 1.e+6 * mpSI
nSI = 1.e+6                         # for Rho
uSI = ReSI                          # for Ux,Uy,Uz
pSI = rhoSI * uSI**2                # for P
bSI = uSI * np.sqrt(mu0SI * rhoSI)  # for Bx,By,Bz,B1x,B1y,B1z
jSI = bSI / (xSI * mu0SI)           # for Jx,Jy,Jz

scale_varnames = {'Rho': nSI/1.e+6,
                  'P':   pSI*1.e+9,
                  'Ux':  uSI/1.e+3,
                  'Uy':  uSI/1.e+3,
                  'Uz':  uSI/1.e+3,
                  'Bx':  bSI*1.e+9,
                  'By':  bSI*1.e+9,
                  'Bz':  bSI*1.e+9,
                  'B1x': bSI*1.e+9,
                  'B1y': bSI*1.e+9,
                  'B1z': bSI*1.e+9,
                  'Jx':  jSI*1.e+6,
                  'Jy':  jSI*1.e+6,
                  'Jz':  jSI*1.e+6}

# coordinates stored as x, y, z in earth radii.
def MODEL():

    from kamodo import Kamodo
    import spacepy.pybats as sp
    import time
    from os.path import basename
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

        The reader requires code in readers/OCTREE_BLOCK_GRID/:
        Compile the library using the command
          python interpolate_amrdata_extension_build.py
        inside the anaconda/miniconda3/miniforge3 environment for Kamodo.
        The shared library file name contains the python version
        and the name of the operation system, e.g.,
          _interpolate_amrdata.cpython-39-darwin.so
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
            if not RU._isfile(list_file) or not RU._isfile(time_file):
                files = sorted(RU.glob(file_dir+'*.out')) +\
                    sorted(RU.glob(file_dir+'*.outs'))  # .out or .outs
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
                    pattern_files = sorted(RU.glob(file_dir+p+'*.out')) +\
                        sorted(RU.glob(file_dir+p+'*.outs'))
                    self.pattern_files[p] = pattern_files
                    self.times[p] = {'start': [], 'end': [], 'all': []}

                    # loop through to get times, one file per time step
                    self.times[p]['start'] = array([(
                        sp.IdlFile(f).attrs[
                        'time'].replace(tzinfo=timezone.utc) - self.filedate
                        ).total_seconds()/3600. for f in pattern_files])
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
            mhd = sp.IdlFile(self.pattern_files[p][0])  # This fails on s3
            file_variables = [key for key in mhd.keys() if key not in
                              ['grid', 'x', 'y', 'z']] + ['th']
            if len(variables_requested) > 0 and variables_requested != 'all':
                self.check_alt_varnames(file_variables)
                gvar_list = [key for key, value in model_varnames.items()
                             if value[0] in variables_requested and
                             key in file_variables]

                # check for variables requested but not available
                if len(gvar_list) != len(variables_requested):
                    err_list = [value[0] for key, value in
                                model_varnames.items() if value[0] in
                                variables_requested and key not in gvar_list]
                    # remove value from list if another renamed value matches
                    for errval in err_list:
                        err_not = [key for key, value in model_varnames.items()
                             if value[0] == errval and key in file_variables]
                        if len(err_not) > 0: err_list.remove(errval)
                    if len(err_list) > 0:
                        print('Some requested variables are not available:',
                              err_list)
            else:  # only input variables on the avoid_list if requested
                self.check_alt_varnames(file_variables)
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

            # store key for each variable desired
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

        def check_alt_varnames(self, file_variables):
            '''
            Some variables can be saved as either dimensioned or dimensionless.
            Lower case (default) is dimensioned. Upper case is dimensionless.
            This includes B1x, B1y, B1z, Bx, By, Bz, Ux, Uy, Uz, P, Rho, Jx, Jy, Jz.
            This function looks to see if the variables in the file match the current
              model_varnames dictionary and if not changes the dictionary to match.
            It then will dimensionalize the variables if needed upon reading.
            '''
            for key, values in alt_varnames.items():
                for value in values:
                    if value in model_varnames: oldvalue = value
                for value in values:
                    if value in file_variables and value not in model_varnames:
                        model_varnames[value] = model_varnames.pop(oldvalue)

            return

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
                    sp.IdlFile(f).meta[gvar]
                    for f in self.pattern_files[key]])
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
                mhd = sp.IdlFile(self.pattern_files[key][i])
                # initialize octree instance for file if not already done
                if i not in self.octree[key].keys():
                    # files is closed after read in spacepy.pybats.IdlFile
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
                    #                  sort_unstructured=False)
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
                # get data from file and initialize
                data = mhd[gvar]

                if 'normalized variables' in mhd.meta['header']:
                    if gvar in scale_varnames:
                        data = data*scale_varnames[gvar]
                if gvar in ['theta1', 'theta2', 'phi1', 'phi2']:
                    dmask = data < -99.
                    data[dmask] = None
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
