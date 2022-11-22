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
    from os.path import basename, isfile
    from numpy import array, NaN, unique, zeros, ravel, ndarray
    from time import perf_counter
    import kamodo_ccmc.readers.reader_utilities as RU
    # from readers.OCTREE_BLOCK_GRID._interpolate_amrdata import ffi
    # from readers.OCTREE_BLOCK_GRID._interpolate_amrdata import lib

    class MODEL(Kamodo):
        '''SWMF model data reader for the magnetoshpere outputs.

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
                files = sorted(glob(file_dir+'*.out*'))  # .out or .outs
                self.filename = ''.join([f+',' for f in files])[:-1]
                patterns = unique([basename(f)[:10] for f in files])
                # get time grid from files
                dt = sp.IdlFile(files[0]).attrs['time']
                if dt is not None:  # filedate given not always at midnight
                    self.filedate = datetime.strptime(
                        dt.isoformat()[:10], '%Y-%m-%d').replace(
                        tzinfo=timezone.utc)
                    delta_t = RU.str_to_hrs(dt.isoformat(), self.filedate,
                                            '%Y-%m-%dT%H:%M:%S')
                else:
                    self.filedate = datetime(1900, 1, 1, tzinfo=timezone.utc)
                    delta_t = 0.

                # establish time attributes from filenames
                for p in patterns:
                    # get list of files to loop through later
                    pattern_files = sorted(glob(file_dir+p+'*.out*'))
                    self.pattern_files[p] = pattern_files
                    self.times[p] = {'start': [], 'end': [], 'all': []}

                    # loop through to get times, one file per time step
                    self.times[p]['start'] = array([RU.tstr_to_hrs(
                        sp.IdlFile(f).attrs['strtime'][:-1].replace(
                            'h', ':').replace('m', ':')) + delta_t
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
            self.constants = {'c': mhd.meta['c'], 'g': mhd.meta['g'],
                              'R': mhd.meta['R']}

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
                # get data from file
                data = sp.IdlFile(self.pattern_files[key][i])[gvar]  # dmarray

                # assign custom interpolator: Lutz Rastaetter 2021
                def interp(xvec):
                    tic = perf_counter()
                    X, Y, Z = xvec.T
                    '''
                    var_data = ffi.new("float[]", list(data))
                    xpos = ffi.new("float[]", list(X))
                    ypos = ffi.new("float[]", list(Y))
                    zpos = ffi.new("float[]", list(Z))
                    npos = len(X)
                    return_data = list(zeros(npos, dtype=float))
                    return_data_ffi = ffi.new("float[]", return_data)
                    # THIS NEEDS TO RETURN THE return_data_ffi VARIABLE!!!!!!!!!!!
                    IS_OK = lib.interpolate_amrdata_multipos(
                      xpos, ypos, zpos, npos, var_data, return_data_ffi)
                    return_data[:] = list(return_data_ffi)  # STILL ZEROS!
                    toc = perf_counter()
                    print(varname,' interp: ', toc-tic, 'seconds elapsed')
                    return return_data
                    '''
                    if not isinstance(X, ndarray):
                        return NaN
                    else:
                        return zeros(len(X)) * NaN
                return interp

            # functionalize the variable dataset
            tmp = self.variables[varname]
            tmp['data'] = zeros((2, 2, 2, 2))  # saves execution time
            self = RU.Functionalize_Dataset(
                self, coord_dict, varname, tmp, gridded_int, coord_str,
                interp_flag=1, func=func, func_default='custom')
            return
    return MODEL
