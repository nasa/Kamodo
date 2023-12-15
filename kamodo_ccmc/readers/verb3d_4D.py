"""
@author: xandrd
"""

# model_varnames = {'Name in file': ['LaTeX representation', 'Description', integer, 'Coordinate system',
#                            'Coordinate grid', [Coordinate list], 'Units'],

model_varnames = {'PSD': ['PSD_lea', 'Phase Space Density in (L, E, A)', 0, 'LEA',
                          'car', ['time', 'L', 'Energy', 'Alpha'], '1/(s*cm**2*keV*sr*MeV**2)'],
                  'PSD_2': ['PSD_lmk', 'Phase Space Density in (L, \\mu, K)', 0, 'LEA',
                            'car', ['time', 'L', 'Mu', 'K'], '1/(s*cm**2*keV*sr*MeV**2)'],
                  'L': ['L', 'L-shell', 1, 'LEA',
                        'car', ['L', 'Energy', 'Alpha'], ''],
                  'L_2': ['L_lmk', 'L-shell', 1, 'LEA',
                          'car', ['L', 'Mu', 'K'], ''],
                  'E': ['Energy', 'Energy', 1, 'LEA',
                        'car', ['L', 'Energy', 'Alpha'], 'MeV'],
                  'Alpha': ['Alpha', 'Pitch angle', 1, 'LEA',
                            'car', ['L', 'Energy', 'Alpha'], 'deg'],
                  'Mu': ['Mu', '1st adiabatic invariant \\mu', 1, 'LMK',
                         'car', ['L', 'Mu', 'K'], 'MeV/G'],
                  'K': ['K', '2dn adiabatic invariant K', 1, 'LMK',
                        'car', ['L', 'Mu', 'K'], 'GR_E**1/2']
                  }


def MODEL():
    from kamodo import Kamodo
    import numpy as np
    from numpy import array, unique, NaN, append, transpose, where
    import os
    import kamodo_ccmc.readers.reader_utilities as RU
    from time import perf_counter

    # main class
    class MODEL(Kamodo):
        '''VERB-3D model data reader.

        Inputs:
            file_dir: a string representing the file directory of the
                model output data.
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
            # gridded_int = boolean (default = True)
            #     - If True, the variables chosen are functionalized in both the
            #         standard method and a gridded method.
            #     - If False, the variables chosen are functionalized in only the
            #         standard method.
            # verbose = boolean (False)
            #     - If False, script execution and the underlying Kamodo
            #         execution is quiet except for specified messages.
            #     - If True, be prepared for a plethora of messages.
        All inputs are described in further detail in
            KamodoOnboardingInstructions.pdf.

        Returns: a kamodo object (see Kamodo core documentation) containing all
            requested variables in functionalized form.

        Notes:
            -
        '''

        def __init__(self, file_dir, variables_requested=[],
                     filetime=False, verbose=False, gridded_int=True,
                     printfiles=False, **kwargs):
            super(MODEL, self).__init__()
            self.modelname = 'VERB-3D'
            t0 = perf_counter()

            # Kamodo does handle windows path correctly. Changing to paths to replacing \\ with /
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
                self._variable_requested_check(variables_requested)
                if len(variables_requested) == 0:
                    return

            # store variables
            # missing_value = NaN
            varfiles = {}  # store which variable came from which file
            self.gvarfiles = {}  # store file variable name similarly
            self.gvarfiles_full = {}  # All {patterns: [variables], ...}
            err_list = []

            self.var_dict = {}  # This could be unused....

            for p in self.patterns:
                _pattern_files = self.pattern_files[p]

                # Kamodo RU.create_timelist discards file path for no reason
                # _pattern_files = [os.path.join(file_dir, 'Output', file) for file in _pattern_files]

                with RU.Dataset(_pattern_files[0], 'r') as cdf_data:
                    # check var_list for variables not possible in this file set
                    self.gvarfiles_full.update({p: [k for k in cdf_data.variables.keys()]})

                    # TODO: Check this condition
                    if len(variables_requested) > 0 and \
                            variables_requested != 'all':
                        _gvar_list = [key for key in model_varnames.keys()
                                     if key in cdf_data.variables.keys() and
                                     model_varnames[key][0] in variables_requested]
                        if len(_gvar_list) != len(variables_requested):
                            _err_list = [value[0] for key, value in
                                        model_varnames.items()
                                        if key not in cdf_data.variables.keys() and
                                        value[0] in variables_requested]
                            err_list.extend(_err_list)  # add to master list
                    else:
                        _gvar_list = [key for key in model_varnames.keys()
                                     if key in cdf_data.variables.keys()]

                # store which file these variables came from
                varfiles[p] = [model_varnames[key][0] for key in _gvar_list]
                self.gvarfiles[p] = _gvar_list

                self.print_err_list(err_list)

            if variables_requested == 'all':
                # collect all possible variables in set of files and return
                self._update_var_dict(varfiles)
                return

            self._print_files(printfiles)

            self._create_variables()

            # Load required grid variables
            self.load_grid(file_dir)

            # register interpolators for each variable
            t_reg = perf_counter()
            # store original list b/c gridded interpolators change keys list
            varname_list = [key for key in self.variables.keys()]
            for varname in varname_list:
                self.register_variable(varname, gridded_int, file_dir)

            if verbose:
                print(f'Took {perf_counter() - t_reg:.5f}s to register ' +
                      f'{len(varname_list)} variables.')
                print(f'Took a total of {perf_counter() - t0:.5f}s to kamodofy' +
                      f' {len(varname_list)} variables.')

        def print_err_list(self, err_list):
            # if variables not found
            if len(err_list) > 0:
                print(f'requested variables are not available: {err_list}')

        def register_variable(self, varname, gridded_int, file_dir):
            """Registers an interpolator with proper signature"""

            # coord_dict = {'time': {'units': 'hr', 'data': time},
            #               'L': {'units': '', 'data': L},
            #               'mu': {'units': 'MeV/G', 'data': mu},
            #               'K': {'units': 'R_E/G**0.5', 'data': K}}
            # var_dict = {'PSD': {'units': '(c/cm/MeV)**3', 'data': PSD}}

            gvar = [key for key, value in model_varnames.items() if
                    value[0] == varname][0]  # variable name in file

            grid_variabels = model_varnames[gvar][5]
            coord_dict = {}

            # Handle time
            key = self.variables[varname]['data']
            if 'time' in grid_variabels:
                coord_dict = {'time': {'units': 'hr',
                                       'data': self.times[key]['all']}}

            # Find corresponding name for the coordinates based on Latex name
            # TODO: this can create ambiguity if latex variables are the same
            gvar_variabels = [gvar for gvar, items in model_varnames.items() if items[0] in grid_variabels]

            for coord in gvar_variabels:
                coord_dict.update({coord: {'units': model_varnames[coord][6], 'data': getattr(self, "_" + coord)}})

            # variable_name = model_varnames[gvar][0]
            coord_str = ''

            def func(i, fi):  # i = file#, fi = slice#
                pattern_files = self.pattern_files[key]

                # Kamodo RU.create_timelist discards file path for no reason
                # pattern_files = [os.path.join(file_dir, 'Output', file) for file in pattern_files]

                with RU.Dataset(self.pattern_files[key][i]) as cdf_data:
                    data = array(cdf_data.variables[gvar][fi])
                return data

            # functionalize the 3D or 4D dataset, series of time slices
            self = RU.Functionalize_Dataset(
                self, coord_dict, varname, self.variables[varname],
                gridded_int, coord_str, interp_flag=3, func=func,
                times_dict=self.times[key])
            return

        def load_grid(self, file_dir):
            """Load grid data according to the requested variables"""

            # Determine the variables
            variables = [v for gvars in self.gvarfiles.values() for v in gvars]

            # Determine the grid variables
            grid_variables = set(item for v in variables for item in model_varnames[v][5])

            # Exclude time
            if 'time' in grid_variables:
                grid_variables.remove('time')

            # Find corresponding name for the coordinates based on Latex name
            # TODO: this can create ambiguity if latex variables are the same
            gvar_variabels = [gvar for gvar, items in model_varnames.items() if items[0] in grid_variables]

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
                            # TODO: consider adding leading _ to avoid naming problems
                            setattr(self, "_" + var, array(cdf_data.variables[var]))

        def _update_var_dict(self, varfiles):
            for p in self.patterns:
                var_list = varfiles[p]

                self.var_dict.update({value[0]: value[0:] for key, value in
                                      model_varnames.items() if value[0] in
                                      var_list})

        def _create_variables(self):
            # initialize storage structure
            self.variables = {}
            for p in self.patterns:
                variables = {model_varnames[gvar][0]: {
                    'units': model_varnames[gvar][-1],
                    'data': p} for gvar in self.gvarfiles[p]}
                self.variables.update(variables)

        def _print_files(self, printfiles):
            # option to print files
            if printfiles:
                print(f'{len(self.filename)} Files:')
                files = self.filename.split(',')
                for f in files:
                    print(f)

        def _variable_requested_check(self, variables_requested):
            # TODO: Write with dataclass
            test_list = [value[0] for key, value in model_varnames.items()]
            err_list = [item for item in variables_requested if item not in test_list]
            if len(err_list) > 0:
                print('Variable name(s) not recognized:', err_list)
            for item in err_list:
                variables_requested.remove(item)

        def time_list_files(self, file_dir):
            # For debug purpose
            FORCE_LOAD = False
            # first, check for file list, create if DNE
            list_file = file_dir + self.modelname + '_list.txt'
            time_file = file_dir + self.modelname + '_times.txt'
            self.times, self.pattern_files = {}, {}
            if FORCE_LOAD or not RU._isfile(list_file) or not RU._isfile(time_file):
                from kamodo_ccmc.readers.verb3d_tocdf import convert_all
                self.conversion_test = convert_all(file_dir)
            return time_file, list_file

    return MODEL


from dataclasses import dataclass


@dataclass
class ModelVariable:
    var: str
    desc: str
    i: int
    sys: str
    grid: str
    coord: list
    units: str


class ModelVariables():
    def __init__(self, model_varnames):
        self.vars = {}
        self.keys = {}

        for key, v in model_varnames.items():
            self.vars[key] = ModelVariable(v[0], v[1], v[2], v[3], v[4], v[5], v[6])
            self.keys[v[0]] = key


_model_vars = ModelVariables(model_varnames)
