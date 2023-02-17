'''
Written by Rebecca Ringuette, 2022
'''

# variable name in file: [standardized variable name, descriptive term, units]
model_varnames = {'PHI': ['Phi', 'Electric potential', 0, 'SM', 'sph',
                          ['time', 'lon', 'lat'], 'kV'],
                  'Psi': ['Psi', 'Magnetic potential', 0, 'SM', 'sph',
                          ['time', 'lon', 'lat'], 'cT*m'],
                  'FAC': ['j_FAC', 'Field-aligned current density', 0, 'SM',
                          'sph', ['time', 'lon', 'lat'], 'muA/m**2'],
                  'JH': ['JH', 'Joule heating rate', 0, 'SM', 'sph',
                         ['time', 'lon', 'lat'], 'mW/m**2']}


def MODEL():

    from kamodo import Kamodo
    from netCDF4 import Dataset
    from glob import glob
    from os.path import basename, isfile
    from numpy import array
    from time import perf_counter
    from datetime import datetime, timezone
    import kamodo_ccmc.readers.reader_utilities as RU

    class MODEL(Kamodo):
        '''Weimer model data reader.

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
        Note: If you want to add a file to a directory already converted to nc,
            you need to delete the converted file and the _list.txt and
            _times.txt files to trigger a new file conversion on the next
            execution.

        Returns: a kamodo object (see Kamodo core documentation) containing all
            requested variables in functionalized form.

        Notes:
            - Weimer model outputs are given in ascii files with one timestep
              per file. Since the files are extremely small, all of the data
              for the run is collected into one netCDF4 file.
            - The converted file is assumed to be larger than 16 GB, so
              interpolation method 3 is chosen. The standard SciPy
              interpolator is used.
        '''
        def __init__(self, file_dir, variables_requested=[],
                     printfiles=False, filetime=False, gridded_int=True,
                     verbose=False, **kwargs):
            super(MODEL, self).__init__(**kwargs)
            self.modelname = 'Weimer'
            t0 = perf_counter()

            # first, check for file list, create if DNE
            list_file = file_dir + self.modelname + '_list.txt'
            time_file = file_dir + self.modelname + '_times.txt'
            self.times, self.pattern_files = {}, {}
            if not isfile(list_file) or not isfile(time_file):
                # check for nc files
                nc_files = sorted(glob(file_dir+'*.nc'))
                txt_files = sorted(glob(file_dir+'*.txt'))
                if len(nc_files) == 0:  # perform file conversion if none
                    from kamodo_ccmc.readers.weimer_tocdf import convert_all
                    convert_all(txt_files)
                    nc_files = sorted(glob(file_dir+'*.nc'))
                self.filename = ''.join([f+',' for f in txt_files])[:-1]
                p = basename(nc_files[0]).split('.')[0]  # only one nc file

                # datetime object for midnight on date from first text file
                date = basename(txt_files[0])[-17:-9]  # YYYYMMDD
                self.filedate = datetime.strptime(
                    date+' 00:00:00', '%Y%m%d %H:%M:%S').replace(
                        tzinfo=timezone.utc)

                # establish time attributes
                # store list of files to loop through later
                self.pattern_files[p] = nc_files
                self.times[p] = {'start': [], 'end': [], 'all': []}

                # get times from single nc file
                cdf_data = Dataset(nc_files[0])
                tmp = array(cdf_data.variables['time'])  # hrs since midnit
                self.times[p]['start'] = array([tmp[0]])
                self.times[p]['end'] = array([tmp[-1]])
                self.times[p]['all'] = tmp
                cdf_data.close()

                # create time list file if DNE
                RU.create_timelist(list_file, time_file, self.modelname,
                                   self.times, self.pattern_files,
                                   self.filedate)
            else:  # read in data and time grids from file list
                self.times, self.pattern_files, self.filedate, self.filename =\
                    RU.read_timelist(time_file, list_file)
            # return time info
            if filetime:
                return

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
            cdf_data = Dataset(self.pattern_files[p][0])

            # get coordinates from first file
            self._lat = array(cdf_data.variables['lat'])
            self._lon = array(cdf_data.variables['lon'])

            # check var_list for variables not possible in this file set
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
                    if len(err_list) > 0:
                        print('Some requested variables are not available: ',
                              err_list)
            else:
                gvar_list = [key for key in model_varnames.keys()
                             if key in cdf_data.variables.keys()]
            # store which file these variables came from
            self.varfiles = [model_varnames[key][0] for key in gvar_list]
            self.gvarfiles = gvar_list
            cdf_data.close()

            # collect all possible variables in set of files and return
            if variables_requested == 'all':
                self.var_dict = {value[0]: value[1:] for key, value in
                                 model_varnames.items() if value[0] in
                                 self.varfiles}
                return

            # store mapping for each variable desired
            self.variables = {model_varnames[gvar][0]: {
                'units': model_varnames[gvar][-1],
                'data': p} for gvar in gvar_list}

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

        # define and register the variable
        def register_variable(self, varname, gridded_int):
            """Registers a 3d interpolator with 3d signature"""

            # define and register the interpolators
            key = self.variables[varname]['data']
            coord_dict = {'time': {'units': 'hr', 'data':
                                   self.times[key]['all']},
                          'lon': {'units': 'deg', 'data': self._lon},
                          'lat': {'units': 'deg', 'data': self._lat}}
            coord_str = [value[3]+value[4] for key, value in
                         model_varnames.items() if value[0] == varname][0]
            gvar = [key for key, value in model_varnames.items() if
                    value[0] == varname][0]  # variable name in file

            def func(i, fi):  # i = file#, fi = slice#
                cdf_data = Dataset(self.pattern_files[key][i])
                data = array(cdf_data.variables[gvar][fi])
                cdf_data.close()
                return data

            # functionalize the 3D or 4D dataset, series of time slices
            self = RU.Functionalize_Dataset(
                self, coord_dict, varname, self.variables[varname],
                gridded_int, coord_str, interp_flag=3, func=func,
                times_dict=self.times[key])
            return

    return MODEL
