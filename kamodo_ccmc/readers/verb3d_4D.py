"""
@author: xandrd
"""
import numpy as np

# model_varnames = {'Name in file': ['LaTeX representation', 'Description', integer, 'Coordinate system',
#                            'Coordinate grid', [Coordinate list], 'Units'],

# It is easier to use pc, Mu, K as grid variables than calculate it within Kamodo.
model_varnames = {'PSD': ['PSD_lea', 'Phase Space Density in (L, E_e, alpha_e)', 0, 'LEA',
                          'rb', ['time', 'L', 'E_e', 'alpha_e'], '(c/MeV/cm)**3'],
                  'Flux': ['flux_lea', 'Electron flux (L, E_e, alpha_e)', 0, 'LEA',
                           'rb', ['time', 'L', 'E_e', 'alpha_e'], '1/(s*cm**2*keV*sr)'],
                  'L': ['L', 'L-shell', 1, 'LEA',
                        'rb', ['L', 'E_e', 'alpha_e'], ''],
                  'E': ['E_e', 'Electron energy', 1, 'LEA',
                        'rb', ['L', 'E_e', 'alpha_e'], 'MeV'],
                  'Alpha': ['alpha_e', 'Equatorial pitch angle', 1, 'LEA',
                            'car', ['L', 'E_e', 'alpha_e'], 'deg'],
                  'pc': ['pc', 'Momentum times speed of light', 1, 'LEA',
                         'rb', ['L', 'E_e', 'alpha_e'], 'MeV'],
                  'PSD_2': ['PSD_lmk', 'Phase Space Density in (L, mu, K)', 0, 'LMK',
                            'rb', ['time', 'L', 'mu', 'K'], '(c/MeV/cm)**3'],
                  'Mu': ['mu', '1st adiabatic invariant mu', 1, 'LMK',
                         'rb', ['L', 'mu', 'K'], 'MeV/G'],  # Kamodo should also not handle units of G...
                  'K': ['K', '2dn adiabatic invariant K', 1, 'LMK',
                        'rb', ['L', 'mu', 'K'], '10**-4*T*(km/6371*km)**1/2']
                  # Kamodo cannot handle units of G * R_E^1/2
                  }

def MODEL():
    from kamodo import Kamodo
    import numpy as np
    from numpy import array, unique, NaN, append, transpose, where
    import os
    import kamodo_ccmc.readers.reader_utilities as RU
    from time import perf_counter
    from dataclasses import dataclass
    from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
    from scipy.spatial import cKDTree

    # main class
    class MODEL(Kamodo):
        '''VERB-3D model data reader.

        Inputs:
            file_dir: a string representing the file directory of the
                model output data.
            variables_requested = a list of variable name strings (LaTeX representation) chosen from
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

        Notes:
            -
        '''

        # Static properties of this mode class
        modelname = 'VERB-3D'

        # Prefix for grid variables
        _grid_prefix = '_grid'

        # Forcing to recompute all nc files
        _force_convert_all = False

        # Forcing start date of the dataset
        _start_date = None

        # Interpolation method for PSD and flux functionalization
        #_interpolation_method = 'nearest'
        _interpolation_method = 'interp'

        def __init__(self, file_dir, variables_requested=[],
                     filetime=False, verbose=False, gridded_int=True,
                     printfiles=False, **kwargs):
            super(MODEL, self).__init__()

            t0 = perf_counter()

            # Kamodo does not handle windows path correctly. Changing to paths to replacing \\ with /
            file_dir = file_dir.replace('\\', '/')

            # Check time and list files, and generate them if needed
            time_file, list_file = self.time_list_files(file_dir)

            # read in data and time grids from file list
            self.times, self.pattern_files, self.filedate, self.filename = \
                RU.read_timelist(time_file, list_file)

            # List of patterns
            self.patterns = list(self.pattern_files.keys())

            # return time info. This behavior is used by the File_Times.
            if filetime:
                return

            # perform initial check on variables_requested list
            if variables_requested != 'all' and len(variables_requested) > 0:
                # Check if the variables_requested is a sting and convert if needed
                if isinstance(variables_requested, str):
                    variables_requested = [variables_requested]

                self._variable_requested_check(variables_requested)
                # Incorrect variables_requested could have been removed
                if len(variables_requested) == 0:
                    return

            # store variables
            # missing_value = NaN
            varfiles = {}  # store which variable came from which file {patterns: [LaTeX requested variables], ...}
            self.gvarfiles = {}  # store file variable name similarly {patterns: [File requested variables], ...}
            self.gvarfiles_full = {}  # All {patterns: [File variables], ...} from cdf files
            # err_list = []

            self.var_dict = {}  # This could be unused....

            # TODO: This needs to be wrapped in the function
            for p in self.patterns:
                # Get the file for the current pattern
                _pattern_files = self.pattern_files[p]

                # Kamodo RU.create_timelist discards file path for no reason
                # _pattern_files = [os.path.join(file_dir, 'Output', file) for file in _pattern_files]

                with RU.Dataset(_pattern_files[0], 'r') as cdf_data:
                    # check var_list for variables not possible in this file set
                    self.gvarfiles_full.update({p: [k for k in cdf_data.variables.keys()]})

                    # TODO: Check this condition. How to get here??
                    if variables_requested != 'all' and len(variables_requested) > 0:
                        # _gvar_list = [key for key in model_varnames.keys()
                        #              if key in cdf_data.variables.keys() and
                        #              model_varnames[key][0] in variables_requested]

                        # list of file variables in cdf that correspond to variables_requested
                        _gvar_list = _model_vars.file_variables(variables_requested=variables_requested,
                                                                cdf_keys=cdf_data.variables.keys())
                    else:
                        # list of file variables in cdf
                        _gvar_list = _model_vars.file_variables(cdf_keys=cdf_data.variables.keys())

                # store which file these variables came from
                # varfiles[p] = [model_varnames[key][0] for key in _gvar_list]
                varfiles[p] = [_model_vars.vars[v].var for v in _gvar_list]
                self.gvarfiles[p] = _gvar_list

            # Identify and print the errors
            self.print_err_list(variables_requested, varfiles)

            if variables_requested == 'all':
                # collect all possible variables in set of files and return
                self._update_var_dict(varfiles)
                return

            self._print_files(printfiles)

            self._create_variables()

            # Load required grid variables
            self.load_grid(file_dir)

            # Load static variables like pc, that do not need to be loaded from nc files
            self.load_static(file_dir)

            # register interpolators for each variable
            t_reg = perf_counter()
            # store original list b/c gridded interpolators change keys list
            varname_list = [key for key in self.variables.keys()]
            for varname in varname_list:
                # Register PSD variables
                if varname in ['PSD_lmk', 'PSD_lea', 'flux_lea']:
                    self.register_ncfile_variable(varname, gridded_int, file_dir)
                # If this is not a PSD variable, it is a grid variable.
                # Grid and static variables handled diffrently
                else:
                    self.register_grid_variable(varname, gridded_int, file_dir)

            if verbose:
                print(f'Took {perf_counter() - t_reg:.5f}s to register ' +
                      f'{len(varname_list)} variables.')
                print(f'Took a total of {perf_counter() - t0:.5f}s to kamodofy' +
                      f' {len(varname_list)} variables.')

        def print_err_list(self, variables_requested, varfiles):
            """
            Print variables_requested that are not in varfiles
            """
            if variables_requested != 'all' and len(variables_requested) > 0:
                err_list = [v for v in variables_requested if v not in
                            [i for sublist in varfiles.values() for i in sublist]]

                # if variables not found
                if len(err_list) > 0:
                    print(f'requested variables are not available: {err_list}')

        def register_ncfile_variable(self, varname, gridded_int, file_dir):
            """Registers an interpolator with proper signature"""

            # coord_dict = {'time': {'units': 'hr', 'data': time},
            #               'L': {'units': '', 'data': L},
            #               'mu': {'units': 'MeV/G', 'data': mu},
            #               'K': {'units': 'R_E/G**0.5', 'data': K}}
            # var_dict = {'PSD': {'units': '(c/cm/MeV)**3', 'data': PSD}}

            coord_dict = {}

            # File pattern of this varname
            pattern_key = self.variables[varname]['data']

            # Name of the variable in the file
            gvar = _model_vars.keys[varname]

            # List of coordinates
            grid_variables = _model_vars.vars[gvar].coord

            # Handle time
            if 'time' in grid_variables:
                coord_dict = {'time': {'units': 'hr',
                                       'data': self.times[pattern_key]['all']}}

            # Find corresponding name for the coordinates based on Latex name
            gvar_grid = [_model_vars.keys[v] for v in grid_variables if v in _model_vars.keys]

            grid = 'LEA'

            if all(var in grid_variables for var in ['L', 'E_e', 'alpha_e']):
                npoints = np.max(self._gridL.shape) * 2
                unique_en = np.sort(np.unique(self._gridE.flatten()))
                unique_A = np.sort(np.unique(self._gridAlpha.flatten()))
                idx_en = np.linspace(0, len(unique_en) - 1, min([len(unique_en), npoints])).astype(int)
                idx_A = np.linspace(0, len(unique_A) - 1, min([len(unique_A), npoints])).astype(int)
                coord_dict.update({'L': {'units': _model_vars.vars['L'].units,
                                         'data': self._gridL[:, 0, 0]},  # L is the same
                                   'E_e': {'units': _model_vars.vars['E'].units,
                                           'data': unique_en[idx_en]},
                                   'alpha_e': {'units': _model_vars.vars['Alpha'].units,
                                               'data': unique_A[idx_A]},
                                   })

            elif all(var in grid_variables for var in ['L', 'mu', 'K']):
                npoints = np.max(self._gridL.shape) * 2
                unique_mu = np.sort(np.unique(self._gridMu.flatten()))
                unique_K = np.sort(np.unique(self._gridK.flatten()))
                idx_mu = np.linspace(0, len(unique_mu) - 1, min([len(unique_mu), npoints])).astype(int)
                idx_K = np.linspace(0, len(unique_K) - 1, min([len(unique_K), npoints])).astype(int)

                coord_dict.update({'L': {'units': _model_vars.vars['L'].units,
                                         'data': self._gridL[:, 0, 0]},  # L is the same
                                   'mu': {'units': _model_vars.vars['Mu'].units,
                                          'data': unique_mu[idx_mu]},
                                   'K': {'units': _model_vars.vars['K'].units,
                                         'data': unique_K[idx_K]},
                                   })
                grid = 'LMK'

            coord_str = ''

            def func_near(i):  # i = file#, fi = slice#
                pattern_files = self.pattern_files[pattern_key]

                # Kamodo RU.create_timelist discards file path for no reason
                # pattern_files = [os.path.join(file_dir, 'Output', file) for file in pattern_files]

                with RU.Dataset(self.pattern_files[pattern_key][i]) as cdf_data:
                    data = array(cdf_data.variables[gvar])

                # Using nearest interpolator
                def nearest_xvec(xvec):
                    # TODO: Change to just return no data_point
                    if grid == 'LEA':
                        data_point = self._nearest_xvec_data_LEA(xvec, data)
                    elif grid == 'LMK':
                        data_point = self._nearest_xvec_data_LMK(xvec, data)
                    return data_point

                # Using interpolator
                def interp_xvec(xvec):
                    if grid == 'LEA':
                        data_point = self._interp_xvec_data_LEA(xvec, data)
                    elif grid == 'LMK':
                        data_point = self._interp_xvec_data_LMK(xvec, data)
                    return data_point

                if self._interpolation_method == 'nearest':
                    return nearest_xvec
                else:
                    return interp_xvec

            # functionalize the 3D or 4D dataset, series of time slices
            self = RU.Functionalize_Dataset(
                self, coord_dict, varname, self.variables[varname],
                gridded_int, coord_str, interp_flag=1, func=func_near,
                times_dict=self.times[pattern_key], func_default='custom')
            return

        def register_grid_variable(self, varname, gridded_int, file_dir):
            """Registers an interpolator with proper signature"""

            coord_dict = {}

            # File pattern of this varname
            pattern_key = self.variables[varname]['data']

            # Name of the variable in the file
            gvar = _model_vars.keys[varname]

            # List of coordinates
            grid_variables = _model_vars.vars[gvar].coord

            if all(var in grid_variables for var in ['L', 'E_e', 'alpha_e']):
                grid = 'LEA'
            elif all(var in grid_variables for var in ['L', 'mu', 'K']):
                grid = 'LMK'

            # Find corresponding name for the coordinates based on Latex name
            gvar_grid = [_model_vars.keys[v] for v in grid_variables if v in _model_vars.keys]

            for coord in gvar_grid:
                coord_dict.update({_model_vars.vars[coord].var:
                                       {'units': _model_vars.vars[coord].units,
                                        'data': getattr(self, self._grid_prefix + coord)}})

            # variable_name = model_varnames[gvar][0]
            coord_str = ''

            def func_near():  # xvec - interpolator input ?

                # Get index of a grid
                data = getattr(self, self._grid_prefix + gvar)

                # Using nearest interpolator
                def nearest_xvec(xvec):
                    # TODO: Change to just return no data_point
                    if grid == 'LEA':
                        data_point = self._nearest_xvec_data_LEA(xvec, data)
                    elif grid == 'LMK':
                        data_point = self._nearest_xvec_data_LMK(xvec, data)
                    return data_point

                return nearest_xvec

            # functionalize the 3D or 4D dataset, series of time slices
            self = RU.Functionalize_Dataset(
                self, coord_dict, varname, self.variables[varname],
                gridded_int=False, coord_str=coord_str, interp_flag=0, func=func_near(),
                times_dict=None, func_default='custom')
            return

        def load_grid(self, file_dir):
            """Load grid data according to the requested variables"""

            # Determine the variables
            variables = [v for gvars in self.gvarfiles.values() for v in gvars]

            # Determine the grid variables
            grid_variables = set(item for v in variables for item in _model_vars.vars[v].coord)

            # Exclude time
            if 'time' in grid_variables:
                grid_variables.remove('time')

            # Find corresponding name for the coordinates based on Latex name
            gvar_variabels = [_model_vars.keys[v] for v in grid_variables if v in _model_vars.keys]

            # Find what patterns to load
            patterns_to_load = set()
            for key, value in self.gvarfiles_full.items():
                if any(item in gvar_variabels for item in value):
                    patterns_to_load.add(key)

            # Load grids
            for p in patterns_to_load:
                pattern_files = self.pattern_files[p]

                # Kamodo RU.create_timelist discards file path for no reason
                # pattern_files = [os.path.join(file_dir, 'Output', file) for file in pattern_files]

                with RU.Dataset(pattern_files[0], 'r') as cdf_data:
                    for var in self.gvarfiles_full[p]:
                        if var in gvar_variabels:
                            setattr(self, self._grid_prefix + var, array(cdf_data.variables[var]))

        def load_static(self, file_dir):
            """Load static variables data according to the requested variables"""

            # Determine the requested variables
            variables = [v for gvars in self.gvarfiles.values() for v in gvars]

            # Find corresponding static variable name based on Latex name
            static_variables = ['pc']
            gvar_variabels = [_model_vars.keys[v] for v in variables if v in static_variables]

            # Find what patterns to load
            patterns_to_load = set()
            for key, value in self.gvarfiles_full.items():
                if any(item in gvar_variabels for item in value):
                    patterns_to_load.add(key)

            # Load static variables
            for p in patterns_to_load:
                pattern_files = self.pattern_files[p]

                with RU.Dataset(pattern_files[0], 'r') as cdf_data:
                    for var in self.gvarfiles_full[p]:
                        if var in gvar_variabels:
                            setattr(self, self._grid_prefix + var, array(cdf_data.variables[var]))

        def _update_var_dict(self, varfiles):
            for p in self.patterns:
                var_list = varfiles[p]

                self.var_dict.update({var: _model_vars.model_var_by_name(var).to_list()[1:]
                                      for var in _model_vars.latex_variables()
                                      if var in var_list})

        def _create_variables(self):
            # initialize storage structure
            self.variables = {}
            for p in self.patterns:
                var = {_model_vars.vars[key].var: {
                    'units': _model_vars.vars[key].units,
                    'data': p} for key in self.gvarfiles[p]}
                self.variables.update(var)

        def _print_files(self, printfiles):
            # option to print files
            if printfiles:
                print(f'{len(self.filename)} Files:')
                files = self.filename.split(',')
                for f in files:
                    print(f)

        def _variable_requested_check(self, variables_requested):
            """
            Initial check of the requested variable.
            Compare requested variables by user with the list of the variables that are provided by the model.
            Prints the variables that were not recognized and removes them from the requested variables list.
            """

            err_list = [item for item in variables_requested
                        if item not in _model_vars.latex_variables()]
            if len(err_list) > 0:
                print('Variable name(s) not recognized:', err_list)
            for item in err_list:
                variables_requested.remove(item)

        def time_list_files(self, file_dir):
            """
            Check if '_list.txt' and '_time.txt' files exist.
            If not, then we don't have generated files for this dataset (run) specified by file_dir.
            Creates dataset files.
            """

            # first, check for file list, create if DNE
            list_file = file_dir + self.modelname + '_list.txt'
            time_file = file_dir + self.modelname + '_times.txt'
            self.times, self.pattern_files = {}, {}
            if self._force_convert_all or not RU._isfile(list_file) or not RU._isfile(time_file):
                from kamodo_ccmc.readers.verb3d_tocdf import convert_all
                # TODO: Do we need self.conversion_test?
                self.conversion_test = convert_all(file_dir, start_date=self._start_date)
            return time_file, list_file

        def _nearest_xvec_data_LEA(self, xvec, data):
            # TODO: Understand how gridded intrepolant takes xvec

            if not isinstance(xvec, np.ndarray):
                xvec = np.array(xvec)  # convert to numpy array

            xvec = np.atleast_2d(xvec)  # Make sure we can iterate over the array

            d1 = []

            for vec in xvec:
                L, E, A = vec

                L_idx = _nearest_index(self._gridL[:, 0, 0], L)
                A_idx = _nearest_index(self._gridAlpha[L_idx, 0, :], A)
                E_idx = _nearest_index(self._gridE[L_idx, :, A_idx], E)

                d1.append(data[L_idx, E_idx, A_idx])

            # TODO: Add verification that indexes are within the range of the grid
            return d1

        def _nearest_xvec_data_LMK(self, xvec, data):
            # TODO: Understand how gridded intrepolant takes xvec

            if not isinstance(xvec, np.ndarray):
                xvec = np.array(xvec)  # convert to numpy array

            xvec = np.atleast_2d(xvec)  # Make sure we can iterate over the array

            d1 = []

            for vec in xvec:
                L, Mu, K = vec

                L_idx = _nearest_index(self._gridL[:, 0, 0], L)
                K_idx = _nearest_index(self._gridK[L_idx, 0, :], K)
                Mu_idx = _nearest_index(self._gridMu[L_idx, :, K_idx], Mu)

                d1.append(data[L_idx, Mu_idx, K_idx])

            # TODO: Add verification that indexes are within the range of the grid
            return d1

        def _intrep_xvec_data_prepare(self, xvec):
            if not isinstance(xvec, np.ndarray):
                xvec = np.array(xvec)

            xvec = np.atleast_2d(xvec)  # Ensure we can iterate over the array
            L, Y, Z = xvec[:, 0], xvec[:, 1], xvec[:, 2]

            # Extract grid data
            L_grid_sorted = np.sort(np.unique(self._gridL))
            min_L_grid, max_L_grid = np.min(L_grid_sorted), np.max(L_grid_sorted)
            num_L_grid = self._gridL.shape[0]

            nobs = L.size

            # Check for NaN values
            nan_check = np.isnan(L) | np.isnan(Y) | np.isnan(Z)
            if np.all(nan_check):
                return L, Y, Z, np.full(nobs, np.NaN), None, None

            # Handle NaN points in the input data
            L[nan_check] = np.nan
            Y[nan_check] = np.nan
            Z[nan_check] = np.nan

            # Determine interpolation limits in the L grid
            Lmin, Lmax = np.nanmin(L), np.nanmax(L)
            limax = np.argmin(np.abs(L_grid_sorted - Lmax))
            limin = np.argmin(np.abs(L_grid_sorted - Lmin))

            limin = max(limin - 1, 0)
            limax = min(limax + 1, num_L_grid - 1)

            if (Lmax > max_L_grid and Lmin > max_L_grid) or (Lmax < min_L_grid and Lmin < min_L_grid):
                return L, Y, Z, np.full(nobs, np.NaN), limin, limax

            return L, Y, Z, None, limin, limax

        def _intrep_stack_hstack(self, Xit, Xi, Yit, Yi, Zit, Zi, PSDit, PSDi):
            # Accumulate grid points and Psi values using np.hstack
            Xit = np.hstack((Xit, Xi.ravel()))
            Yit = np.hstack((Yit, Yi.ravel()))
            Zit = np.hstack((Zit, Zi.ravel()))
            PSDit = np.hstack((PSDit, PSDi.ravel()))

            return Xit, Yit, Zit, PSDit

        def _intrep_stack_list(self, Xit_list, Xi, Yit_list, Yi, Zit_list, Zi, PSDit_list, PSDi):
            # Collect grid points and Psi values in lists
            Xit_list.append(Xi.ravel())
            Yit_list.append(Yi.ravel())
            Zit_list.append(Zi.ravel())
            PSDit_list.append(PSDi.ravel())
            return Xit_list, Yit_list, Zit_list, PSDit_list

        def _finalize_stack(Xit_list, Yit_list, Zit_list, PSDit_list):
            # Convert lists to numpy arrays after all data is collected
            Xit = np.hstack(Xit_list)
            Yit = np.hstack(Yit_list)
            Zit = np.hstack(Zit_list)
            PSDit = np.hstack(PSDit_list)
            return Xit, Yit, Zit, PSDit

        def _interp_xvec_data_LEA(self, xvec, data):

            L, E, A, output, limin, limax = self._intrep_xvec_data_prepare(xvec)
            if output is not None:
                return output

            # Initialize arrays for storing grid points and corresponding PSD values
            Xit, Yit, Zit, PSDit = np.array([]), np.array([]), np.array([]), np.array([])

            # Determine the range of the input coordinates
            A_min, A_max = np.nanmin(A), np.nanmax(A)
            E_min, E_max = np.nanmin(E), np.nanmax(E)


            # Loop over L grid indices and gather relevant grid points for interpolation
            for l in range(limin, limax + 1):
                alpha_values = np.sort(np.unique(self._gridAlpha[l, :, :]))
                energy_values = np.sort(np.unique(self._gridE[l, :, :]))

                num_alpha = alpha_values.size
                num_energy = energy_values.size

                aimax = np.argmin(np.abs(alpha_values - A_max))
                aimin = np.argmin(np.abs(alpha_values - A_min))
                eimax = np.argmin(np.abs(energy_values - E_max))
                eimin = np.argmin(np.abs(energy_values - E_min))

                # Boundary limiting
                eimin1 = max(eimin - 1, 0)
                eimax1 = min(eimax + 1, num_energy - 1)
                aimin1 = max(aimin - 1, 0)
                aimax1 = min(aimax + 1, num_alpha - 1)

                PSDi = np.log10(data[l, eimin1:eimax1 + 1, aimin1:aimax1 + 1])
                Xi = self._gridL[l, eimin1:eimax1 + 1, aimin1:aimax1 + 1]
                Yi = np.log10(self._gridE[l, eimin1:eimax1 + 1, aimin1:aimax1 + 1])
                Zi = self._gridAlpha[l, eimin1:eimax1 + 1, aimin1:aimax1 + 1]

                Xit, Yit, Zit, PSDit =  self._intrep_stack_hstack(Xit, Xi, Yit, Yi, Zit, Zi, PSDit, PSDi)

            # Perform interpolation using the common interpolation logic
            return self._interpolate(L, np.log10(E), A, Xit, Yit, Zit, PSDit)

        def _interp_xvec_data_LMK(self, xvec, data):

            L, M, K, output, limin, limax = self._intrep_xvec_data_prepare(xvec)
            if output is not None:
                return output

            # Initialize arrays for storing grid points and corresponding PSD values
            Xit, Yit, Zit, PSDit = np.array([]), np.array([]), np.array([]), np.array([])

            # Determine the range of the input coordinates
            M_min, M_max = np.nanmin(M), np.nanmax(M)
            K_min, K_max = np.nanmin(K), np.nanmax(K)

            # Loop over L grid indices and gather relevant grid points for interpolation
            for l in range(limin, limax + 1):
                K_arr = self._gridK[l, 0, :]  # K arr must be unique and monotonic for this to work
                # K_values = np.sort(np.unique(self._gridK[l, :, :]))[::-1]
                num_K = K_arr.size
                Kimax = np.argmin(np.abs(K_arr - K_max))  # lowers index
                Kimin = np.argmin(np.abs(K_arr - K_min))  # highest index

                M_values_min = self._gridMu[l, :, Kimax]  # Lowest Mu at highest K
                M_values_max = self._gridMu[l, :, Kimin]  # Highest Mu at lowest K
                Mimax = np.argmin(np.abs(M_values_max - M_max))
                Mimin = np.argmin(np.abs(M_values_min - M_min))
                num_M = M_values_min.size


                # Boundary limiting
                Kimin1 = min(Kimin + 1, num_K - 1)
                Kimax1 = max(Kimax - 1, 0)

                Mimin1 = max(Mimin - 1, 0)
                Mimax1 = min(Mimax + 1, num_M - 1)

                # K indexing from max to min, because grid is monotonically decreasing
                PSDi = np.log10(data[l, Mimin1:Mimax1+1, Kimax1:Kimin1+1])
                Xi = self._gridL[l, Mimin1:Mimax1+1, Kimax1:Kimin1+1]
                Yi = np.log10(self._gridMu[l, Mimin1:Mimax1+1, Kimax1:Kimin1+1])
                Zi = self._gridK[l, Mimin1:Mimax1+1, Kimax1:Kimin1+1]

                # Accumulate grid points and PSDit values
                Xit, Yit, Zit, PSDit = self._intrep_stack_hstack(Xit, Xi, Yit, Yi, Zit, Zi, PSDit, PSDi)

            # Perform interpolation using the common interpolation logic
            return self._interpolate(L, np.log10(M), K, Xit, Yit, Zit, PSDit)

        def _interpolate(self, X, Y, Z, x, y, z, Psi_tt):
            """Common interpolation logic for both LEA and LMK methods."""

            # Perform linear interpolation
            linear_interp = LinearNDInterpolator(list(zip(x, y, z)), Psi_tt)
            PSD_interp = linear_interp(X, Y, Z)

            # Handle NaNs using nearest-neighbor interpolation if needed
            nan_mask = np.isnan(PSD_interp) & ~(np.isnan(X) | np.isnan(Y) | np.isnan(Z))
            if np.any(nan_mask):
                # Normalize the coordinates only for proximity search
                x_range = np.ptp(x)  # Peak-to-peak (max - min) range
                y_range = np.ptp(y)
                z_range = np.ptp(z)

                # Normalize the grid and input coordinates
                x_normalized = x / x_range
                y_normalized = y / y_range
                z_normalized = z / z_range

                X_normalized = X / x_range
                Y_normalized = Y / y_range
                Z_normalized = Z / z_range

                # Create a KDTree for efficient proximity search
                kdtree = cKDTree(list(zip(x_normalized, y_normalized, z_normalized)))

                # Find the nearest neighbors for the NaN points
                distances, nearest_indices = kdtree.query(
                    list(zip(X_normalized[nan_mask], Y_normalized[nan_mask], Z_normalized[nan_mask])))

                # Define a tolerance for "close" proximity
                tolerance = 1e-1  # Adjust based on your precision requirements

                # Apply nearest-neighbor interpolation for points within the tolerance
                close_enough = distances < tolerance
                if np.any(close_enough):
                    nearest_interp = NearestNDInterpolator(list(zip(x, y, z)), Psi_tt)
                    close_indices = np.where(nan_mask)[0][close_enough]
                    PSD_interp[close_indices] = nearest_interp(X[close_indices], Y[close_indices], Z[close_indices])

            # Return the interpolated values in their original scale
            return 10 ** PSD_interp

    # get nearest index
    def _nearest_index(arr, val):
        idx = np.abs(arr - val).argmin()
        return idx

    @dataclass
    class ModelVariable:
        """ Class of the Model Variables. This cals replaces the list of kamodo variable"""

        var: str  # LaTeX representation
        desc: str  # Description
        i: int  # integer
        sys: str  # Coordinate system
        grid: str  # Coordinate grid
        coord: list  # [Coordinate list]
        units: str  # Units

        def to_list(self):
            return [self.var, self.desc, self.i, self.sys, self.grid, self.coord, self.units]

    # TODO: change the class so different LaTeX variables can correspond to the same key in file
    class ModelVariables:
        """ Class which service the model variables"""

        def __init__(self, model_varnames):
            self.vars = {}  # ModelVariable
            self.keys = {}  # Name in file according to the variable name

            for key, v in model_varnames.items():
                self.vars[key] = ModelVariable(v[0], v[1], v[2], v[3], v[4], v[5], v[6])
                self.keys[v[0]] = key

        def latex_variables(self, variables_requested=None):
            """
            List of latex formatted variables
            If requested_variable is specified, return only matching variables
            """
            if variables_requested is None:
                return [v.var for v in self.vars.values()]
            else:
                return [v.var for v in self.vars.values() if variables_requested in v.var]

        def file_variables(self, variables_requested=None, cdf_keys=None):
            """
            List of files variables
            If requested_variable is specified, return only matching variables
            If cdf_key is specified, return matching variables that correspond to cdf_keys
            If only cdf_key is specified, return all variables that correspond to cdf_keys
            """

            if cdf_keys and variables_requested:
                return [self.keys[v] for v in variables_requested
                        if v in self.keys and self.keys[v] in cdf_keys]
            elif cdf_keys is None:
                return [self.keys[v] for v in variables_requested if v in self.keys]
            elif variables_requested is None:
                return [v for v in self.keys.values() if v in cdf_keys]
            else:
                return list(self.keys.values())

        def model_var_by_name(self, name):
            """ Get model variable by its name """
            key = self.keys[name]
            return self.vars[key]

    _model_vars = ModelVariables(model_varnames)

    return MODEL
