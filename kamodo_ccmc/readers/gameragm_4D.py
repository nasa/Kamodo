# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 15:37:59 2023

@author: rringuet
"""
from numpy import vectorize
from datetime import datetime, timezone

model_varnames = {'Bx': ['B_x', 'X-component of magnetic field',
                         0, 'SM', 'car', ['time', 'X', 'Y', 'Z'], 'nT'],
                  'By': ['B_y', 'Y-component of magnetic field',
                         0, 'SM', 'car', ['time', 'X', 'Y', 'Z'], 'nT'],
                  'Bz': ['B_z', 'Z-component of magnetic field',
                         0, 'SM', 'car', ['time', 'X', 'Y', 'Z'], 'nT'],
                  'Cs': ['c_s', 'Sound speed of plasma', 0, 'SM', 'car',
                         ['time', 'X', 'Y', 'Z'], 'km/s'],
                  'D': ['N_plasma', 'Plasma number density (M/mp)', 0, 'SM',
                        'car', ['time', 'X', 'Y', 'Z'], '1/cm**3'],
                  'Jx': ['J_x', 'X-component of current density',
                         0, 'SM', 'car', ['time', 'X', 'Y', 'Z'], 'nA/m**2'],
                  'Jy': ['J_y', 'Y-component of current density',
                         0, 'SM', 'car', ['time', 'X', 'Y', 'Z'], 'nA/m**2'],
                  'Jz': ['J_z', 'Z-component of current density',
                         0, 'SM', 'car', ['time', 'X', 'Y', 'Z'], 'nA/m**2'],
                  'P': ['P', 'Pressure', 0, 'SM', 'car',
                        ['time', 'X', 'Y', 'Z'], 'nPa'],
                  'Pb': ['P_mag', 'Magnetic pressure', 0, 'SM', 'car',
                         ['time', 'X', 'Y', 'Z'], 'nPa'],
                  'Vx': ['v_x', 'X-component of velocity',
                         0, 'SM', 'car', ['time', 'X', 'Y', 'Z'], 'km/s'],
                  'Vy': ['v_y', 'Y-component of velocity',
                         0, 'SM', 'car', ['time', 'X', 'Y', 'Z'], 'km/s'],
                  'Vz': ['v_z', 'Z-component of velocity',
                         0, 'SM', 'car', ['time', 'X', 'Y', 'Z'], 'km/s'],
                  # CONSTANTS below this line.
                  'dV': ['dV', 'Simulation cell volume', 0, 'SM', 'car',
                         ['X', 'Y', 'Z'], 'R_E**3']
                  }

constants = ['X', 'Y', 'Z', 'dV']


@vectorize
def timestr_timestamp(time_str):
    '''Converts time string from data file into an UTC timestamp. Cuts off
    millisecond portion of times to prevent errors later.'''
    dt = datetime.fromisoformat(time_str).replace(tzinfo=timezone.utc)
    return dt.timestamp()


def MODEL():
    from kamodo import Kamodo
    from glob import glob
    import h5py
    from netCDF4 import Dataset
    from os.path import isfile, basename
    from astropy.time import Time
    from datetime import datetime, timezone
    from numpy import array, NaN, ndarray, zeros, linspace, repeat
    from time import perf_counter
    import kamodo_ccmc.readers.reader_utilities as RU
    import kamodo_ccmc.readers.gameragm_grids as G

    class MODEL(Kamodo):
        '''GAMERA GM model data reader.

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

        Notes and instructions:
        - The GAMERA global magnetosphere outputs are given in one or more h5
          files, each containing all of the time steps for the entire run. If
          the model is run in serial model, only one file will be produced. If
          the model is run in MPI mode, then multiple files will be produced
          with the grid of the simulation sectioned off into one piece per
          file. No file conversion is attempted, but some pre-processing is
          performed on the coordinate grids.
        - The files are typically larger than 16 GB, so interpolation method 3
          is chosen for the time-varying variables. Interpolation method 0 is
          chosen for the constants (e.g. dV). The coordinate grid is not
          uniform in any manner, so a custom interpolator is required for both
          interpolation methods.
        - Adding this interpolator to the reader is a work in process via a
          collaboration with the GAMERA modeling team.
        '''
        def __init__(self, file_dir, variables_requested=[],
                     printfiles=False, filetime=False, gridded_int=True,
                     verbose=False, **kwargs):
            super(MODEL, self).__init__(**kwargs)
            self.modelname = 'GAMERA_GM'
            t0 = perf_counter()

            # first, check for file list, create if DNE
            list_file = file_dir + self.modelname + '_list.txt'
            time_file = file_dir + self.modelname + '_times.txt'
            self.times, self.pattern_files = {}, {}
            if not isfile(list_file) or not isfile(time_file):
                # collect filenames
                files = sorted(glob(file_dir+'*.h5'))
                data_files = [f for f in files if 'Res' not in f]
                self.filename = ''.join([f+',' for f in data_files])[:-1]
                # one pattern per run: abcd_00nx_00ny_00nz
                # abcd part might not always be four characters
                p = ''.join([b+'_' for b in basename(
                    data_files[0]).split('.')[0].split('_')[:-3]])[:-1]

                # establish time attributes
                # get list of files to loop through later
                self.pattern_files[p] = data_files
                self.times[p] = {'start': [], 'end': [], 'all': []}

                # all times are in each file, so just use first file
                h5_data = h5py.File(data_files[0])
                timestep_keys = [key for key in h5_data.keys()
                                 if 'Step' in key]
                mjd_list = [h5_data[key].attrs['MJD'] for key in
                            timestep_keys]
                h5_data.close()

                # convert from modified julian date to UTC time
                # output format is 'YYYY-MM-DD HH:MM:SS.mms
                t = sorted(Time(mjd_list, format='mjd').iso)
                self.filedate = datetime.fromisoformat(t[0][:10]).replace(
                    tzinfo=timezone.utc)  # date only
                # convert to hrs since midnight
                utc_ts = timestr_timestamp(t)
                hrs = (utc_ts - self.filedate.timestamp()) / 3600.
                # store in self.times dictionary
                self.times[p]['start'] = repeat([hrs[0]], len(data_files))
                self.times[p]['end'] = repeat([hrs[-1]], len(data_files))
                self.times[p]['all'] = array(hrs)

                # need to figure out if millisecond accuracy is needed
                test = (len(utc_ts) - 1)/(utc_ts[-1] - utc_ts[0])  # #steps/s
                if test >= 1.5:  # need ms timing if >= 1.5 timesteps/second
                    ms_timing = True
                else:
                    ms_timing = False

                # Pre-calculate the ranges of the coordinate grids.
                # Needed to generate a default uniform grid for plotting.
                print('Computing grid ranges...', end="")
                G.GridRanges(data_files)
                print('done.')

                # Pre-calculate the centers of the coordinate grids for each
                # file separately. Needed for interpolation later.
                print('Computing cell centers...', end="")
                G.ComputeCellCentersFile(data_files)
                print('done.')

                # create time list file if DNE
                RU.create_timelist(list_file, time_file, self.modelname,
                                   self.times, self.pattern_files,
                                   self.filedate, ms_timing=ms_timing)
            else:  # read in data and time grids from file list
                self.times, self.pattern_files, self.filedate, self.filename =\
                    RU.read_timelist(time_file, list_file)
            if filetime:
                return  # return times only

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
            h5_data = h5py.File(self.pattern_files[p][0])
            key_list = [key for key in h5_data.keys() if 'Step' not in key and
                        key not in ['X', 'Y', 'Z']]  # skip coordinates
            step_list = list(h5_data['Step#0'].keys())
            var_list = key_list + step_list
            h5_data.close()

            if len(variables_requested) > 0 and variables_requested != 'all':
                gvar_list = [key for key, value in model_varnames.items()
                             if value[0] in variables_requested and
                             key in var_list]  # file variable names

                # check for variables requested but not available
                if len(gvar_list) != len(variables_requested):
                    err_list = [value[0] for key, value in
                                model_varnames.items() if value[0] in
                                variables_requested and key not in gvar_list]
                    if len(err_list) > 0:
                        print('Some requested variables are not available:',
                              err_list)
            else:  # only input variables on the avoid_list if requested
                gvar_list = [key for key in var_list if key in
                             model_varnames.keys()]
                # return list of variables included in data files
                if variables_requested == 'all':
                    self.var_dict = {value[0]: value[1:] for key, value in
                                     model_varnames.items()
                                     if key in gvar_list}
                    return

            # initialize data mapping for each variable desired
            self.variables = {model_varnames[var][0]:
                              {'units': model_varnames[var][-1],
                               'data': p} for var in gvar_list}

            # read in coordinate grid ranges and make sample grids
            # X, Y, and Z are the same shape in the data because they are the
            # x, y, and z coordinates of the same set of points.
            X_min, X_max, Y_min, Y_max, Z_min, Z_max =\
                G.Read_GridRanges(file_dir)
            if len(self.pattern_files[p]) > 1:
                sample_var = [gvar for gvar in model_varnames.keys()
                              if gvar in step_list][0]
                net_indices = G.GridMapping(self.pattern_files[p], False,
                                            sample_var)
                shape = net_indices['net'][1]  # shape of total array
            else:
                h5_data = h5py.File(self.pattern_files[p][0])
                shape = list(h5_data['X'].shape)  # shape of total array
                h5_data.close()
            self._X = linspace(X_min, X_max, endpoint=True, num=shape[0])
            self._Y = linspace(Y_min, Y_max, endpoint=True, num=shape[1])
            self._Z = linspace(Z_min, Z_max, endpoint=True, num=shape[2])

            # store a few items
            self.missing_value = NaN
            if verbose:
                print(f'Took {perf_counter()-t0:.6f}s to read in data')
            if printfiles:
                print(self.filename)

            # GAMERA data has all the timestamps in all the files.
            # This will confuse the lazy interpolation, so stripping down the
            # times dict to one file entry (only start and end fields needed).
            # For ms_timing, 'all' has ms resolution, start and end do not.
            self.singletimes = {'start': [self.times[p]['all'][0]],
                                'end': [self.times[p]['all'][-1]]}

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
            gvar = [key for key, value in model_varnames.items() if
                    value[0] == varname][0]
            coord_str = [value[3]+value[4] for key, value in
                         model_varnames.items() if value[0] == varname][0]
            file_dir = self.pattern_files[key][0].split(
                basename(self.pattern_files[key][0]))[0]
            if 'time' in coord_list:
                coord_dict = {'time': {'units': 'hr',
                                       'data': self.times[key]['all']},
                              'X': {'units': 'R_E', 'data': self._X},
                              'Y': {'units': 'R_E', 'data': self._Y},
                              'Z': {'units': 'R_E', 'data': self._Z}}

                def func(i, fi):  # i = file# (always 0), fi = slice# (= Step#)
                    # Leave data in blocks. Need interpolator to determine
                    # which block is needed. Setting default of block 0 for now
                    block = 0  # default, but need this value from the interp
                    # read in data from block
                    file = self.pattern_files[key][block]
                    h5_data = h5py.File(file)
                    data = array(h5_data['Step#'+str(fi)][gvar])
                    h5_data.close()
                    # read in cell centers for requested block
                    center_file = file_dir + basename(file).split('.')[0] +\
                        '_gridcenters.nc'
                    cdf_data = Dataset(center_file)
                    X_c = array(cdf_data['X_c'])
                    Y_c = array(cdf_data['Y_c'])
                    Z_c = array(cdf_data['Z_c'])
                    cdf_data.close()

                    # define custom interpolator here  ***********************************
                    def interp(xvec):
                        tic = perf_counter()
                        X, Y, Z = xvec.T  # xvec can be used like this
                        # call custom interpolator here
                        print('variable', i, fi, data.shape, X_c.shape)
                        # dummy code for now
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
                    interp_flag=3, func=func, times_dict=self.singletimes,
                    func_default='custom')
                # choosing large chunked format, same as WACCMX.
            else:  # functionalize the constants
                coord_dict = {'X': {'units': 'R_E', 'data': self._X},
                              'Y': {'units': 'R_E', 'data': self._Y},
                              'Z': {'units': 'R_E', 'data': self._Z}}

                def func():
                    # Leave data in blocks. Need interpolator to determine
                    # which block is needed. Setting default of block 0 for now
                    block = 0  # default, but need this value from the interp
                    # read in data from block
                    file = self.pattern_files[key][block]
                    h5_data = h5py.File(file)
                    data = array(h5_data[gvar])
                    h5_data.close()
                    # read in cell centers for requested block
                    center_file = file_dir + basename(file).split('.')[0] +\
                        '_gridcenters.nc'
                    cdf_data = Dataset(center_file)
                    X_c = array(cdf_data['X_c'])
                    Y_c = array(cdf_data['Y_c'])
                    Z_c = array(cdf_data['Z_c'])
                    cdf_data.close()

                    # define custom interpolator here  ***********************************
                    def interp(xvec):
                        tic = perf_counter()
                        X, Y, Z = array(xvec).T  # xvec can be used like this
                        # call custom interpolator here
                        print('constant', data.shape, X_c.shape)
                        # dummy code for now
                        if not isinstance(X, ndarray):
                            return NaN
                        else:
                            return zeros(len(X)) * NaN
                    return interp

                # functionalize the variable dataset
                tmp = self.variables[varname]
                tmp['data'] = zeros((2, 2, 2))  # saves execution time

                self = RU.Functionalize_Dataset(
                    self, coord_dict, varname, tmp, gridded_int, coord_str,
                    interp_flag=0, func=func, func_default='custom')
            return
    return MODEL
