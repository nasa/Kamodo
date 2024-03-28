'''
Original version: Lutz Raestatter Oct 1(?), 2021
Modify to work with flythrough: Oct 5, 2021 (Rebecca Ringuette)
'''
from datetime import datetime, timezone, timedelta

# standard model dictionary for reference
model_varnames = {'sigh': ['Sigma_H', 'Hall conductance',
                         0, 'SM', 'sph', ['time', 'lon', 'lat'], 'S'],
                  'sigp': ['Sigma_P', 'Pedersen conductance',
                         1, 'SM', 'sph', ['time', 'lon', 'lat'], 'S'],
                  'pot': ['Phi', 'electrostatic potential',
                         2, 'SM', 'sph', ['time', 'lon', 'lat'], 'kV'],
                  'delbr': ['DeltaB_r', 'r component of magnetic perturbaion on the ground',
                         3, 'SM', 'sph', ['time', 'lon', 'lat'], 'nT'],
                  'delbt': ['DeltaB_theta', 'theta (latitude pointing south) component of magnetic perturbaion on the ground',
                         3, 'SM', 'sph', ['time', 'lon', 'lat'], 'nT'],
                  'delbp': ['DeltaB_phi', 'phi (longitude) component of magnetic perturbaion on the ground',
                         3, 'SM', 'sph', ['time', 'lon', 'lat'], 'nT'],
                  'pacurr': ['JR', 'radial electric current',
                         3, 'SM', 'sph', ['time', 'lon', 'lat'], 'A/m**2'],
                  'xjh': ['Q_Joule', 'Joule heating rate',
                         3, 'SM', 'sph', ['time', 'lon', 'lat'], 'W/m**2'],
                  'ppio': ['P_mag', 'mapped plasma pressure',
                         3, 'SM', 'sph', ['time', 'lon', 'lat'], 'pPa'],
                  'rrio': ['N_mag', 'mapped plasma number density',
                         3, 'SM', 'sph', ['time', 'lon', 'lat'], '1/cm**3'],
                  'ttio': ['T_mag', 'mapped Plasma Temperature',
                         3, 'SM', 'sph', ['time', 'lon', 'lat'], 'K'],
                  'vdown': ['V_r', 'radial velocity',
                         3, 'SM', 'sph', ['time', 'lon', 'lat'], 'm/s'],
                  'xjh': ['JouleHeat', 'Joule heating rate',
                         3, 'SM', 'sph', ['time', 'lon', 'lat'], 'mW/cm**2'],
                  'prec_e_e0_1': ['Phi_eE01', 'electron energy flux precipitation 1',
                                  12, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'mW/m**2'],
                  'prec_e_e0_2': ['Phi_eE02', 'electron energy flux precipitation 2 ',
                                  12, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'mW/**2'],
                  'prec_e_fe_1': ['Phi_eE01', 'electron number flux precipitation 1',
                                  12, 'GSE', 'car', ['time', 'x', 'y', 'z'], '1/cm**2/s'],
                  'prec_e_fe_2': ['Phi_eE02', 'electron number flux precipitation 2 ',
                                  12, 'GSE', 'car', ['time', 'x', 'y', 'z'], '1/cm**2/s'],
                  }

def MODEL():
    from numpy import array, unique, squeeze
    from time import perf_counter
    from os.path import basename
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
            if not RU._isfile(list_file) or not RU._isfile(time_file):
                # collect filenames
                files = sorted(RU.glob(file_dir+'*.nc'))
                if len(files) == 0:
                    try:
                        from kamodo_ccmc.readers.openggcm_ie_tocdf import \
                            convert_all as ieconv
                        self.conversion_test = ieconv(file_dir)
                    except:
                        print('The file converter for the OpenGGCM global ' +
                              'magnetosphere outputs currently only runs on ' +
                              'CCMC machines. Please contact CCMC for the ' +
                              'files for the desired run converted to ' +
                              'netCDF4 from Lutz Rastaetter.')
                        return
                # *.3df_YYYY-MM-DD_hh_mm_ss.nc
                patterns = unique([basename(file[:-22]) for file in files])
                self.filename = ''.join([f+',' for f in files])[:-1]
                model_warmup = 7200. # seconds of preconditioning precending STARTTIME
                runme_file=file_dir+'/runme'
                if not RU._isfile(runme_file):
                    runme_file=file_dir+'/../runme'
                    if not RU._isfile(runme_file):
                        print("openggcm_ie_4D: required 'runme' file not found!")
                        return False
                runme_file_object=open(runme_file,'r')
                runme_inputs=runme_file_object.readlines()
                runme_file_object.close()
                run_starttime_line = [s for s in  runme_inputs if s.find("STARTTIME") == 0 ]
                run_starttime = ((run_starttime_line[-1].split())[1]).split(":")
                self.run_start_datetime = datetime.strptime(
                    ((run_starttime_line[-1].split())[1])[0:19],
                    '%Y:%m:%d:%H:%M:%S'
                    ).replace(tzinfo=timezone.utc)-timedelta(hours=2)

                #self.filedate = datetime.strptime(
                #    ((run_starttime_line[-1].split())[1])[0:10]+' 00:00:00',
                #    '%Y:%m:%d %H:%M:%S'
                #    ).replace(tzinfo=timezone.utc)
                self.filedate = datetime.strptime(
                    files[0][-22:-12]+' 00:00:00', '%Y-%m-%d %H:%M:%S'
                    ).replace(tzinfo=timezone.utc)

                # establish time attributes
                for p in patterns:
                    # get list of files to loop through later
                    pattern_files = sorted(RU.glob(file_dir+p+'*.nc'))
                    self.pattern_files[p] = pattern_files
                    self.times[p] = {'start': [], 'end': [], 'all': []}

                    # loop through to get times
                    self.times[p]['start'] = array([RU.str_to_hrs(
                        f[-22:-3], self.filedate,
                        format_string='%Y-%m-%d_%H_%M_%S') for f in
                        pattern_files])
                    #self.times[p]['start'] = array(
                    #    [( (run_start_datetime-filedate).total_seconds()
                    #       + int(f[-6:]) - model_warmup
                    #      )/3600. for f in pattern_files ] )
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
            cdf_data = RU.Dataset(pattern_files[0], 'r')

            # check var_list for variables not possible in this file set
            self.err_list = []
            if len(variables_requested) > 0 and\
                    variables_requested != 'all':
                gvar_list = [key for key in model_varnames.keys()
                             if key in cdf_data.keys() and
                             model_varnames[key][0] in variables_requested]
                if len(gvar_list) != len(variables_requested):
                    err_list = [value[0] for key, value in
                                model_varnames.items()
                                if key not in cdf_data.keys() and
                                value[0] in variables_requested]
                    self.err_list.extend(err_list)  # add to master list
            else:
                gvar_list = [key for key in model_varnames.keys()
                             if key in cdf_data.keys()]

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
            for grid in ['_lon', '_lat']:
                setattr(self, grid, array(cdf_data[grid]))
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
            coord_dict['lon'] = {'data': self._lon, 'units': 'R_E'}
            coord_dict['lat'] = {'data': self._lat, 'units': 'R_E'}

            coord_str = [value[3]+value[4] for key, value in
                         model_varnames.items() if value[0] == varname][0]

            def func(i):
                '''i is the file/time number. OpenGGCM-mag is one file per
                timestep.'''
                # get data from file
                file = self.pattern_files[key][i]
                cdf_data = RU.Dataset(file)
                data = array(cdf_data[gvar])
                cdf_data.close()
                return squeeze(data)

            # define and register the interpolators
            self = RU.Functionalize_Dataset(self, coord_dict, varname,
                                            self.variables[varname],
                                            gridded_int, coord_str,
                                            interp_flag=1, func=func)
            return
    return MODEL
