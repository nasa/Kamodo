"""
@author: xandrd
"""

# model_varnames = {'Name in file': ['LaTeX representation', 'Description', integer, 'Coordinate system',
#                            'Coordinate grid', [Coordinate list], 'Units'],

model_varnames = {'PSD': ['PSD_{lea}', 'Phase Space Density in (L, E, A)', 0, 'LEA',
                          'car', ['time', 'L', 'Energy', 'Alpha'], '1/(s*cm**2*keV*sr*MeV**2)'],
                  'PSD_2': ['PSD_{lmk}', 'Phase Space Density in (L, \\mu, K)', 0, 'LEA',
                            'car', ['time', 'L', 'Mu', 'K'], '1/(s*cm**2*keV*sr*MeV**2)'],
                  'L': ['L', 'L-shell', 1, 'LEA',
                        'car', ['L', 'Energy', 'Alpha'], ''],
                  'L_2': ['L', 'L-shell', 1, 'LEA',
                          'car', ['L', 'Mu', 'K'], ''],
                  'E': ['E', 'Energy', 1, 'LEA',
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
            # filetime = boolean (default = False)
            #     - If False, the script fully executes.
            #     - If True, the script only executes far enough to determine the
            #         time values associated with the chosen data.
            # printfiles = boolean (default = False)
            #     - If False, the filenames associated with the data retrieved
            #         ARE NOT printed.
            #     - If True, the filenames associated with the data retrieved ARE
            #         printed.
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

            # For debug purpose
            FORCE_LOAD = False

            # first, check for file list, create if DNE
            list_file = file_dir + self.modelname + '_list.txt'
            time_file = file_dir + self.modelname + '_times.txt'
            self.times, self.pattern_files = {}, {}
            if FORCE_LOAD or not RU._isfile(list_file) or not RU._isfile(time_file):
                from kamodo_ccmc.readers.verb3d_tocdf import convert_all
                self.conversion_test = convert_all(file_dir)

            # read in data and time grids from file list
            self.times, self.pattern_files, self.filedate, self.filename = \
                RU.read_timelist(time_file, list_file)

            # return time info
            if filetime:
                return

            # perform initial check on variables_requested list
            if len(variables_requested) > 0 and variables_requested != 'all':
                test_list = [value[0] for key, value in model_varnames.items()]
                err_list = [item for item in variables_requested if item not in test_list]
                if len(err_list) > 0:
                    print('Variable name(s) not recognized:', err_list)
                for item in err_list:
                    variables_requested.remove(item)
                if len(variables_requested) == 0:
                    return

            # store variables
            self.missing_value = NaN
            self.varfiles = {}  # store which variable came from which file
            self.gvarfiles = {}  # store file variable name similarly
            self.gvarfiles_full = {}  # All {patterns: [variables], ...}
            self.err_list = []
            self.var_dict = {}

            self.patterns = list(self.pattern_files.keys())

            for p in self.patterns:
                pattern_files = self.pattern_files[p]

                # Kamodo RU.create_timelist discards file path for no reason
                pattern_files = [os.path.join(file_dir, 'Output', file) for file in pattern_files]

                with RU.Dataset(pattern_files[0], 'r') as cdf_data:
                    # check var_list for variables not possible in this file set
                    self.gvarfiles_full.update({p: [k for k in cdf_data.variables.keys()]})

                    # TODO: Check this condition
                    if len(variables_requested) > 0 and \
                            variables_requested != 'all':
                        gvar_list = [key for key in model_varnames.keys()
                                     if key in cdf_data.variables.keys() and
                                     model_varnames[key][0] in variables_requested]
                        if len(gvar_list) != len(variables_requested):
                            err_list = [value[0] for key, value in
                                        model_varnames.items()
                                        if key not in cdf_data.variables.keys() and
                                        value[0] in variables_requested]
                            self.err_list.extend(err_list)  # add to master list
                    else:
                        gvar_list = [key for key in model_varnames.keys()
                                     if key in cdf_data.variables.keys()]

                # store which file these variables came from
                self.varfiles[p] = [model_varnames[key][0] for key in gvar_list]
                self.gvarfiles[p] = gvar_list

                # print message if variables not found
                if len(self.err_list) > 0:
                    print('Some requested variables are not available: ',
                          self.err_list)

            # collect all possible variables in set of files and return
            if variables_requested == 'all':
                for p in self.patterns:
                    var_list = self.varfiles[p]

                    self.var_dict.update({value[0]: value[0:] for key, value in
                                          model_varnames.items() if value[0] in
                                          var_list})
                return

            # option to print files
            if printfiles:
                print(f'{len(self.filename)} Files:')
                files = self.filename.split(',')
                for f in files:
                    print(f)

            # initialize storage structure
            self.variables = {}
            for p in self.patterns:
                variables = {model_varnames[gvar][0]: {
                    'units': model_varnames[gvar][-1],
                    'data': p} for gvar in self.gvarfiles[p]}
                self.variables.update(variables)

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
            if verbose:
                print(f'Took a total of {perf_counter() - t0:.5f}s to kamodofy' +
                      f' {len(varname_list)} variables.')

        def register_variable(self, varname, gridded_int, file_dir):
            """Registers an interpolator with proper signature"""

            # coord_dict = {'time': {'units': 'hr', 'data': time},
            #               'L': {'units': '', 'data': L},
            #               'mu': {'units': 'MeV/G', 'data': mu},
            #               'K': {'units': 'R_E/G**0.5', 'data': K}}
            # var_dict = {'PSD': {'units': '(c/cm/MeV)**3', 'data': PSD}}

            grid_variabels = model_varnames[varname][5]
            coord_dict = {}

            # Handle time
            key = self.variables[varname]['data']
            if 'time' in grid_variabels:
                coord_dict = {'time': {'units': 'hr',
                                       'data': self.times[key]['all']}}

            for coord in grid_variabels:
                if coord != 'time':
                    coord_dict.update({coord: {'units': model_varnames[coord][6], 'data': getattr(self, coord)}})

            variable_name = model_varnames[varname][0]
            coord_str = ''

            def func(i, fi):  # i = file#, fi = slice#
                pattern_files = self.pattern_files[key]

                # Kamodo RU.create_timelist discards file path for no reason
                pattern_files = [os.path.join(file_dir, 'Output', file) for file in pattern_files]

                with RU.Dataset(self.pattern_files[key][i]) as cdf_data:
                    data = array(cdf_data.variables[varname][fi])
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

            # Find what patterns to load
            patterns_to_load = set()
            for key, value in self.gvarfiles_full.items():
                if any(item in grid_variables for item in value):
                    patterns_to_load.add(key)

            # Load grids
            for p in patterns_to_load:
                pattern_files = self.pattern_files[p]

                # Kamodo RU.create_timelist discards file path for no reason
                pattern_files = [os.path.join(file_dir, 'Output', file) for file in pattern_files]

                with RU.Dataset(pattern_files[0], 'r') as cdf_data:
                    for var in self.gvarfiles_full[p]:
                        if var in grid_variables:
                            # TODO: consider adding leading _ to avoid naming problems
                            setattr(self, var, array(cdf_data.variables[var]))

    return MODEL
