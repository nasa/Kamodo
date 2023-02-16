'''
Written by Rebecca Ringuette, 2021
'''

# variable name in file: [standardized variable name, descriptive term, units]
model_varnames = {'PED': ['Sigma_P', 'Pedersen conductance ', 0, 'SM', 'sph',
                          ['time', 'lon', 'lat'], 'S'],
                  'HALL': ['Sigma_H', 'hall Conductance', 1, 'SM', 'sph',
                           ['time', 'lon', 'lat'], 'S'],
                  'PHI': ['phi', 'Electric potential', 2, 'SM', 'sph',
                          ['time', 'lon', 'lat'], 'kV'],
                  'EEAST': ['E_east', 'electric field in eastern direction ' +
                            '(increasing longitude) ', 3, 'SM', 'sph',
                            ['time', 'lon', 'lat'], 'mV/m'],
                  'ENORTH': ['E_north', 'electric field in north direction ' +
                             '(increasing latitude) ', 4, 'SM', 'sph',
                             ['time', 'lon', 'lat'], 'mV/m'],
                  'JEAST': ['J_east', 'electric current in eastern direction' +
                            '(increasing longitude) (height integrated ' +
                            'current density)', 5, 'SM', 'sph',
                            ['time', 'lon', 'lat'], 'A/m'],
                  'JNORTH': ['J_north', 'electric current in north direction' +
                             ' (increasing latitude) (height integrated ' +
                             'current density)', 6, 'SM', 'sph',
                             ['time', 'lon', 'lat'], 'A/m'],
                  'EFLUX': ['Phi', 'energy flux', 7, 'SM', 'sph',
                            ['time', 'lon', 'lat'], 'mW/m**2'],
                  'JHEAT': ['J_heat', 'Joule heating rate ', 8, 'SM', 'sph',
                            ['time', 'lon', 'lat'], 'mW/m**2'],
                  'JRIN': ['J_Rin', 'Input radial current',
                           9, 'SM', 'sph', ['time', 'lon', 'lat'], 'mA/m**2'],
                  'JROUT': ['J_Rout', 'model-generated radial current ' +
                            'would be identical to J_Rin if the model were ' +
                            'perfect', 10, 'SM', 'sph',
                            ['time', 'lon', 'lat'], 'mA/m**2'],
                  }


# times from file converted to seconds since midnight of filedate
# plotting input times will be datetime strings of format 'YYYY-MM-DD HH:mm:ss'
# filedate is self.filedate from adelphi object
# converts to hours since midnight of filedate for plotting
def MODEL():

    from kamodo import Kamodo
    from netCDF4 import Dataset
    from glob import glob
    from os.path import basename, isfile
    from numpy import array, unique, append
    from time import perf_counter
    from datetime import datetime, timezone
    import kamodo_ccmc.readers.reader_utilities as RU

    class MODEL(Kamodo):
        '''ADELPHI model data reader.

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
        Note: If you want to add a hemisphere file to a file date already
            present in the converted file list, you need to delete that
            converted file and the _list.txt and _times.txt files to trigger
            a file conversion on the next execution.

        Returns: a kamodo object (see Kamodo core documentation) containing all
            requested variables in functionalized form.

        Notes:
            - ADELPHI model outputs are produced in ascii form with one file
              per each N/S hemisphere per day. The file converter combines the
              data from both hemispheres into one netCDF4 file per day.
            - The data is only given within 40-60 degrees of the poles, so a
              buffer row of the same values is added in the data to avoid
              losing the most equatorward ring of data in the interpolation.
            - The converted files are small and contain multiple time steps per
              file, so interpolation method 2 is chosen. The standard SciPy
              interpolator is used.
        '''
        def __init__(self, file_dir, variables_requested=[],
                     printfiles=False, filetime=False, gridded_int=True,
                     verbose=False, **kwargs):
            super(MODEL, self).__init__(**kwargs)
            self.modelname = 'ADELPHI'
            t0 = perf_counter()

            # first, check for file list, create if DNE
            list_file = file_dir + self.modelname + '_list.txt'
            time_file = file_dir + self.modelname + '_times.txt'
            self.times, self.pattern_files = {}, {}
            if not isfile(list_file) or not isfile(time_file):
                # check for nc files
                nc_files = sorted(glob(file_dir+'*.nc'))
                if len(nc_files) == 0:  # perform file conversion if none
                    from kamodo_ccmc.readers.adelphi_tocdf import convert_all
                    self.conversion_test = convert_all(file_dir)
                    nc_files = sorted(glob(file_dir+'*.nc'))
                # check for new files needing conversion by comparing dates
                nc_dates = unique([basename(f).split('.')[0][-8:] for f in
                                   nc_files])
                txt_files = sorted(glob(file_dir+'*.txt'))
                txt_dates = unique([basename(f).split('.')[0][-8:] for f in
                                    txt_files])
                if not all(txt_dates == nc_dates):
                    from kamodo_ccmc.readers.adelphi_tocdf import convert_all
                    self.conversion_test = convert_all(file_dir)
                    nc_files = sorted(glob(file_dir+'*.nc'))
                self.filename = ''.join([f+',' for f in nc_files])[:-1]

                # prepare time metadata
                patterns = unique([basename(f).split('.')[0][:-8] for f in
                                   nc_files])
                p = patterns[0]  # only one file pattern, so simplify code

                # datetime object for midnight on date
                date = basename(nc_files[0]).split('.')[0][-8:]  # 'YYYY-MM-DD'
                self.filedate = datetime.strptime(date+' 00:00:00',
                                                  '%Y%m%d %H:%M:%S').replace(
                                                      tzinfo=timezone.utc)

                # establish time attributes
                # store list of files to loop through later
                self.pattern_files[p] = nc_files
                self.times[p] = {'start': [], 'end': [], 'all': []}

                # loop through files to get times
                for i, f in enumerate(nc_files):
                    cdf_data = Dataset(f)
                    tmp = array(cdf_data.variables['time'])  # hrs since midnit
                    print(i, f, tmp[0], tmp[-1], len(tmp))
                    self.times[p]['start'].append(tmp[0] + 24.*i)
                    self.times[p]['end'].append(tmp[-1] + 24.*i)
                    self.times[p]['all'].extend(tmp + 24.*i)
                    cdf_data.close()
                # convert timestamps to be hrs since midnight of 1st file
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

        # define and register a 3D variable
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

            # function to read in data, time chunked files
            def func(i):
                '''key is the file pattern, start_idxs is a list of one or two
                indices matching the file start times in self.start_times[key].
                '''
                # get data from file
                file = self.pattern_files[key][i]
                cdf_data = Dataset(file)
                data = array(cdf_data.variables[gvar])
                cdf_data.close()
                # if not the last file, tack on first time from next file
                if file != self.pattern_files[key][-1]:  # interp btwn files
                    next_file = self.pattern_files[key][i+1]
                    cdf_data = Dataset(next_file)
                    data_slice = array(cdf_data.variables[gvar][0])
                    cdf_data.close()
                    data = append(data, [data_slice], axis=0)
                # data wrangling done in file converter
                return data

            # functionalize the variable data (time chunked interpolation)
            self = RU.Functionalize_Dataset(
                self, coord_dict, varname, self.variables[varname],
                gridded_int, coord_str, interp_flag=2, func=func,
                times_dict=self.times[key])
            return

    return MODEL
