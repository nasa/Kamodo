# -*- coding: utf-8 -*-
"""
Created on Mon May 17 18:53:35 2021

@author: rringuet
"""
from datetime import datetime, timezone
from os.path import isfile, basename

model_varnames = {"Sigma_H": ['Sigma_H', '3D Hall conductivity',
                              0, 'SM', 'sph', ['time', 'lon', 'lat'], "S"],
                  "Sigma_P": ['Sigma_P', '3D Pederson conductivity',
                              1, 'SM', 'sph', ['time', 'lon', 'lat'], "S"],
                  "Phi_E": ['Phi_E', 'energy flux',
                            2, 'SM', 'sph', ['time', 'lon', 'lat'], "W/m**2"],
                  "AveE_avgE": ['E_avg', 'average energy',
                                3, 'SM', 'sph', ['time', 'lon', 'lat'], 'eV'],
                  "j_R": ["j_R", 'radial current density',
                          4, 'SM', 'sph', ['time', 'lon', 'lat'], "muA/m**2"],
                  "Phi": ["phi", 'electric potential',
                          5, 'SM', 'sph', ['time', 'lon', 'lat'], "kV"],
                  "E_x": ["E_x", 'electric field, x component',
                          6, 'SM', 'sph', ['time', 'lon', 'lat'], "mV/m"],
                  "E_y": ["E_y", 'electric field, y component',
                          7, 'SM', 'sph', ['time', 'lon', 'lat'], "mV/m"],
                  "E_z": ["E_z", 'electric field, z component',
                          8, 'SM', 'sph', ['time', 'lon', 'lat'], "mV/m"],
                  "j_x": ["j_x", 'current density, x component',
                          9, 'SM', 'sph', ['time', 'lon', 'lat'], "muA/m**2"],
                  "j_y": ["j_y", 'current density, y component',
                          10, 'SM', 'sph', ['time', 'lon', 'lat'], "muA/m**2"],
                  "j_z": ["j_z", 'current density, z component',
                          11, 'SM', 'sph', ['time', 'lon', 'lat'], "muA/m**2"],
                  "v_x": ['v_x', 'total velocity, x component',
                          12, 'SM', 'sph', ['time', 'lon', 'lat'], "km/s"],
                  "v_y": ['v_y', 'total velocity, y component',
                          13, 'SM', 'sph', ['time', 'lon', 'lat'], "km/s"],
                  "v_z": ['v_z', 'total velocity, z component',
                          14, 'SM', 'sph', ['time', 'lon', 'lat'], "km/s"],
                  "Q_Joule": ['W_JouleH', 'height integrated joule heating',
                              15, 'SM', 'sph', ['time', 'lon', 'lat'],
                              "mW/m**2"],
                  "Phi_nion": ['Phi_Nion', 'flux of ions in number density',
                               16, 'SM', 'sph', ['time', 'lon', 'lat'],
                               "1/cm**2/s"],
                  "Binv_RT": ['Binv_RT', 'inverse magnetic field (RT) ' +
                              '(ray tracing integrated along field line)',
                              17, 'SM', 'sph', ['time', 'lon', 'lat'], "1/T"],
                  "rho_RT": ['rho_RTamu', 'molecular mass density (RT) ' +
                             '(ray tracing integrated along field line)',
                             18, 'SM', 'sph', ['time', 'lon', 'lat'],
                             "amu/cm**3"],
                  "P_RT": ['P_RT', 'pressure (RT) ' +
                           '(ray tracing integrated along field line)',
                           19, 'SM', 'sph', ['time', 'lon', 'lat'], "Pa"],
                  "dLat_star": ['lat_star', 'conjugate latitude',
                                20, 'SM', 'sph', ['time', 'lon', 'lat'],
                                "deg"],
                  "dlon_star": ['lon_star', 'conjugate longitude',
                                21, 'SM', 'sph', ['time', 'lon', 'lat'],
                                "deg"]
                  }
'''
Documentation variables:
(ignored, given in R_E on a unit sphere)
"X":['x', 'km'], "Y":['y', 'km'], "Z":['z', 'km'],
"Theta":['theta',"deg"],"Psi":['psi',"deg"],     (used as coordinates)
"Btilt_theta":['theta_Btilt',"deg"], "Btilt_psi":['psi_Btilt',"deg"]
left in files, not in kamodo object
'''


def MODEL():
    from numpy import array, unique
    from glob import glob
    from time import perf_counter
    from netCDF4 import Dataset
    from kamodo import Kamodo
    import kamodo_ccmc.readers.reader_utilities as RU

    class MODEL(Kamodo):
        '''SWMF_IE model data reader.

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

        Notes:
            - SWMF Ionosphere Electrodynamics outputs are given in .tec files
              with one timestep per file. Each file is converted into a netCDF4
              file.
            - The outputs do not provide values at the poles, so scalar
              averaging is used to determine these values.
            - The converted files are small and are created with one time step
              per file, so interpolation method 1 is chosen. The standard SciPy
              interpolator is used.
        '''
        def __init__(self, file_dir, variables_requested=[],
                     filetime=False, verbose=False, gridded_int=True,
                     printfiles=False, **kwargs):
            super(MODEL, self).__init__(**kwargs)
            self.modelname = 'SWMF_IE'
            t0 = perf_counter()

            # first, check for file list, create if DNE
            list_file = file_dir + self.modelname + '_list.txt'
            time_file = file_dir + self.modelname + '_times.txt'
            self.times, self.pattern_files = {}, {}
            if not isfile(list_file) or not isfile(time_file):
                # find unconverted files and convert them
                nc_files = sorted(glob(file_dir+'*.nc'))
                tec_files = sorted(glob(file_dir+'*.tec'))
                if len(nc_files) != len(tec_files):
                    from kamodo_ccmc.readers.swmfie_tocdf import \
                        convert_all
                    convert_all(file_dir)
                else:
                    print('All files already converted.')

                # continue
                files = sorted(glob(file_dir+'*.nc'))
                patterns = unique([basename(f)[:-22] for f in files])
                self.filename = ''.join([f+',' for f in files])[:-1]
                self.filedate = datetime.strptime(
                    basename(files[0])[-22:-7], '%Y%m%d-%H%M%S').replace(
                        tzinfo=timezone.utc)

                # establish time attributes from filenames
                for p in patterns:
                    # get list of files to loop through later
                    pattern_files = sorted(glob(file_dir+p+'*.nc'))
                    self.pattern_files[p] = pattern_files
                    self.times[p] = {'start': [], 'end': [], 'all': []}

                    # loop through to get times, one time per file
                    for f in pattern_files:
                        time = RU.str_to_hrs(f[-22:-7], self.filedate,
                                             format_string='%Y%m%d-%H%M%S')
                        self.times[p]['start'].append(time)
                        self.times[p]['end'].append(time)
                        self.times[p]['all'].append(time)  # one time per file
                    self.times[p]['start'] = array(self.times[p]['start'])
                    self.times[p]['end'] = array(self.times[p]['end'])
                    self.times[p]['all'] = array(self.times[p]['all'])

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
            cdf_data = Dataset(self.pattern_files[p][0])

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

            # get list of variables possible in these files using first file
            if len(variables_requested) > 0 and variables_requested != 'all':
                gvar_list = [key for key, value in model_varnames.items()
                             if value[0] in variables_requested and
                             key in cdf_data.variables.keys()]

                # check for variables requested but not available
                if len(gvar_list) != len(variables_requested):
                    err_list = [value[0] for key, value in
                                model_varnames.items()
                                if value[0] in variables_requested and
                                key not in cdf_data.variables.keys()]
                    if len(err_list) > 0:
                        print('Some requested variables are not available:',
                              err_list)
            else:
                gvar_list = [key for key in cdf_data.variables.keys()
                             if key in model_varnames.keys()]
                if variables_requested == 'all':
                    self.var_dict = {value[0]: value[1:] for key, value in
                                     model_varnames.items() if key in
                                     gvar_list}
                    return

            # Store variable's data and units, transposing the 2D+time array.
            self.variables = {model_varnames[key][0]: {
                'units': model_varnames[key][-1], 'data': p} for key in
                gvar_list}
            # self.theta_Btilt = cdf_data.theta_Btilt   # leaving out
            # self.psi_Btilt = cdf_data.psi_Btilt

            # initialize
            self._registered = 0
            if printfiles:
                print('Files:')
                for file in self.filename:
                    print(file)

            # Store coordinate data as class attributes
            self._lat = array(cdf_data.variables['lat'])
            self._lon = array(cdf_data.variables['lon'])
            cdf_data.close()
            if verbose:
                print(f'Took {perf_counter()-t0:.6f}s to read in data')

            # register interpolators for each variable
            # store original list b/c gridded interpolators
            varname_list = [key for key in self.variables.keys()]
            t_reg = perf_counter()
            for varname in varname_list:  # all are 3D variables
                self.register_variable(varname, gridded_int)
            if verbose:
                print(f'Took {perf_counter()-t_reg:.5f}s to register ' +
                      f'{len(varname_list)} variables.')
            if verbose:
                print(f'Took a total of {perf_counter()-t0:.5f}s to kamodofy' +
                      f' {len(gvar_list)} variables.')

        # define and register a 3D variable------------------------------------
        def register_variable(self, varname, gridded_int):
            """Registers an interpolator. SWMF_IE only has 3D data."""

            # determine coordinate variables and xvec by coord list
            key = self.variables[varname]['data']
            coord_dict = {'time': {'units': 'hr',
                                   'data': self.times[key]['all']}}
            gvar = [key for key, value in model_varnames.items() if
                    value[0] == varname][0]
            coord_str = [value[3]+value[4] for key, value in
                         model_varnames.items() if value[0] == varname][0]
            coord_dict['lon'] = {'units': 'deg', 'data': self._lon}
            coord_dict['lat'] = {'units': 'deg', 'data': self._lat}

            # define operations for each variable when given the key
            def func(i):
                '''i is the file number.'''
                # get data from file
                file = self.pattern_files[key][i]
                cdf_data = Dataset(file)
                data = array(cdf_data.variables[gvar])
                cdf_data.close()
                # data wrangling all done in the file conversion step
                return data

            # functionalize the variable dataset
            self = RU.Functionalize_Dataset(
                self, coord_dict, varname, self.variables[varname],
                gridded_int, coord_str, interp_flag=1, func=func,
                times_dict=self.times[key])
            return
    return MODEL
