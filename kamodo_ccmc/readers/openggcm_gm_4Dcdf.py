'''
Original version: Lutz Raestatter Oct 1(?), 2021
Modify to work with flythrough: Oct 5, 2021 (Rebecca Ringuette)
'''
from datetime import datetime, timezone

# standard model dictionary for reference
model_varnames = {'bx': ['B_x', 'x component of magnetic field',
                         0, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'nT'],
                  'by': ['B_y', 'y component of magnetic field',
                         1, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'nT'],
                  'bz': ['B_z', 'z component of magnetic field',
                         2, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'nT'],
                  'ex': ['E_x', 'x component of electric field',
                         6, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'mV/m'],
                  'ey': ['E_y', 'y component of electric field',
                         7, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'mV/m'],
                  'ez': ['E_z', 'z component of electric field',
                         8, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'mV/m'],
                  'vx': ['v_plasmax', 'x component of plasma velocity',
                         9, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'km/s'],
                  'vy': ['v_plasmay', 'y component of plasma velocity',
                         10, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'km/s'],
                  'vz': ['v_plasmaz', 'z component of plasma velocity',
                         11, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'km/s'],
                  'rr': ['N_plasma', 'number density of plasma ' +
                         '(hydrogen equivalent)',
                         12, 'GSE', 'car', ['time', 'x', 'y', 'z'], '1/cm**3'],
                  'resis': ['eta', 'resistivity',
                            13, 'GSE', 'car', ['time', 'x', 'y', 'z'],
                            'm**2/s'],
                  'pp': ['P_plasma', 'plasma pressure',
                         14, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'pPa'],
                  'xjx': ['j_x', 'current density, x component',
                          15, 'GSE', 'car', ['time', 'x', 'y', 'z'],
                          'muA/m**2'],
                  'xjy': ['j_y', 'current density, y component',
                          16, 'GSE', 'car', ['time', 'x', 'y', 'z'],
                          'muA/m**2'],
                  'xjz': ['j_z', 'current density, z component',
                          17, 'GSE', 'car', ['time', 'x', 'y', 'z'],
                          'muA/m**2']}


def MODEL():
    from numpy import array, unique, squeeze
    from time import perf_counter
    from os.path import isfile, basename
    from glob import glob
    from netCDF4 import Dataset
    from kamodo import Kamodo
    import kamodo_ccmc.readers.reader_utilities as RU

    class MODEL(Kamodo):
        '''OpenGGCM_GM magnetosphere reader.

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

        Notes and special features:
            - The file converter for the OpenGGCM global magnetosphere outputs
              (compressed binary files) currently only runs on CCMC machines.
              Please contact CCMC for the desired run to be converted to
              netCDF4 (Lutz Rastaetter).
            - This model reader has two special properties called
              kamodo_object.near_Earth_boundary_radius and
              kamodo_object.near_Earth_boundary_radius_unit that give the
              inner boundaries of the radial domain for the given run. The
              inner boundary will also be readily apparent when viewing any
              plot including the coordinate origin (X, Y, Z) = (0, 0, 0). The
              unit of the inner boundary is typically earth radii (R_E).
            - The model outputs are produced with one time step per file, so
              interpolation method 1 is chosen. The standard SciPy interpolator
              is used.
                
        '''
        def __init__(self, file_dir, variables_requested=[],
                     filetime=False, verbose=False, gridded_int=True,
                     printfiles=False, **kwargs):
            super(MODEL, self).__init__()
            self.modelname = 'OpenGGCM_GM'
            t0 = perf_counter()  # profiling time stamp

            # first, check for file list, create if DNE
            list_file = file_dir + self.modelname + '_list.txt'
            time_file = file_dir + self.modelname + '_times.txt'
            self.times, self.pattern_files = {}, {}
            if not isfile(list_file) or not isfile(time_file):
                # collect filenames
                files = sorted(glob(file_dir+'*.nc'))
                if len(files) == 0:
                    try:
                        from kamodo_ccmc.readers.openggcm_to_cdf import \
                            openggcm_combine_magnetosphere_files as gmconv
                        self.conversion_test = gmconv(file_dir)
                    except:
                        print('The file converter for the OpenGGCM global ' +
                              'magnetosphere outputs currently only runs on ' +
                              'CCMC machines. Please contact CCMC for the ' +
                              'files for the desired run converted to ' +
                              'netCDF4 from Lutz Rastaetter.')
                        return
                patterns = unique([basename(file[:-19]) for file in files])
                self.filename = ''.join([f+',' for f in files])[:-1]
                self.filedate = datetime.strptime(
                    files[0][-19:-9]+' 00:00:00', '%Y-%m-%d %H:%M:%S'
                    ).replace(tzinfo=timezone.utc)

                # establish time attributes
                for p in patterns:
                    # get list of files to loop through later
                    pattern_files = sorted(glob(file_dir+p+'*.nc'))
                    self.pattern_files[p] = pattern_files
                    self.times[p] = {'start': [], 'end': [], 'all': []}

                    # loop through to get times
                    self.times[p]['start'] = array([RU.str_to_hrs(
                        f[-19:-3]+'_00', self.filedate,
                        format_string='%Y-%m-%d_%H_%M_%S') for f in
                        pattern_files])
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
                return  # return times only

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

            # there is only one pattern for OpenGGCM, so just save the one grid
            p = list(self.pattern_files.keys())[0]
            pattern_files = self.pattern_files[p]
            cdf_data = Dataset(pattern_files[0], 'r')

            # check var_list for variables not possible in this file set
            self.err_list = []
            if len(variables_requested) > 0 and\
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
            self.varfiles = [model_varnames[key][0] for
                             key in gvar_list]
            self.gvarfiles = gvar_list
            # initialize storage structure
            self.variables = {model_varnames[gvar][0]: {
                'units': model_varnames[gvar][-1], 'data': p} for gvar in
                self.gvarfiles}

            # get coordinate grids
            self.near_Earth_boundary_radius = \
                cdf_data.near_Earth_boundary_radius
            self.near_Earth_boundary_radius_unit = \
                cdf_data.near_Earth_boundary_radius_units
            for grid in ['_x', '_y', '_z']:
                setattr(self, grid, array(cdf_data.variables[grid]))
            cdf_data.close()

            # print message if variables not found
            if len(self.err_list) > 0:
                print('Some requested variables are not available: ',
                      self.err_list)

            # collect all possible variables in set of files and return
            if variables_requested == 'all':
                self.var_dict = {value[0]: value[1:] for key, value in
                                 model_varnames.items() if value[0] in
                                 self.varfiles}
                return

            # option to print files
            if printfiles:
                print(f'{len(self.filename)} Files:')
                files = self.filename.split(',')
                for f in files:
                    print(f)

            # register interpolators for each variable
            t_reg = perf_counter()
            # store original list b/c gridded interpolators change keys list
            varname_list = [key for key in self.variables.keys()]
            for varname in varname_list:
                self.register_variable(varname, gridded_int)

            if verbose:
                print(f'Took {perf_counter()-t_reg:.5f}s to register ' +
                      f'{len(varname_list)} variables.')
            if verbose:
                print(f'Took a total of {perf_counter()-t0:.5f}s to kamodofy' +
                      f' {len(varname_list)} variables.')

        # define and register a 4D variable (all are 4D)
        def register_variable(self, varname, gridded_int):
            '''Register and functionalize the variable data.'''

            # determine which file the variable came from, retrieve the coords
            key = self.variables[varname]['data']
            gvar = [key for key, value in model_varnames.items() if
                    value[0] == varname][0]  # variable name in file
            coord_dict = {'time': {'units': 'hr',
                                   'data': self.times[key]['all']}}
            coord_dict['X'] = {'data': self._x, 'units': 'R_E'}
            coord_dict['Y'] = {'data': self._y, 'units': 'R_E'}
            coord_dict['Z'] = {'data': self._z, 'units': 'R_E'}
            coord_str = [value[3]+value[4] for key, value in
                         model_varnames.items() if value[0] == varname][0]

            def func(i):
                '''i is the file/time number. OpenGGCM-mag is one file per
                timestep.'''
                # get data from file
                file = self.pattern_files[key][i]
                cdf_data = Dataset(file)
                data = array(cdf_data.variables[gvar])
                cdf_data.close()
                return squeeze(data)

            # define and register the interpolators
            self = RU.Functionalize_Dataset(self, coord_dict, varname,
                                            self.variables[varname],
                                            gridded_int, coord_str,
                                            interp_flag=1, func=func)
            return
    return MODEL
