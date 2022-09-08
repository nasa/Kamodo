'''
TIEGCM Kamodo reader, adapted to new structure for satellite flythrough
    software
Initial version - Asher Pembroke (?)
Initial version of model_varnames contributed by Zachary Waldron
New code: Rebecca Ringuette (June 2021 and on)

NOTE:
    The current logic for variables that depend on imlev slices off e36 values
        in self._imlev coordinate array. This only works because there is one
        variable that depends on imlev: H_imlev. The logic on lines 311-313
        will have to be reworked a bit if other variables depend on imlev
        later.
'''
from numpy import vectorize
from datetime import datetime, timezone, timedelta

model_varnames = {  # 4D Variables, vert coordinate on midpoint levels (ilev)
                  # geometric height- interpolated to the mid points
                  "ZGMID": ["H_ilev",
                            'height dependent on primary pressure level',
                            0, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                            "cm"],
                  "TN": ["T_n_ilev", 'neutral temperature',
                         1, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'], "K"],
                  "TN2": ["T_n", 'neutral temperature', 1, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'height'], "K"],
                  "O2": ["mmr_O2_ilev",
                         'mass mixing ratio of molecular oxygen', 2, 'GDZ',
                         'sph', ['time', 'lon', 'lat', 'ilev'], ""],
                  "O22": ["mmr_O2", 'mass mixing ratio of molecular oxygen',
                          2, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                          ""],
                  "O1": ["mmr_O_ilev", 'mass mixing ratio of atomic oxygen',
                         3, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'], ""],
                  "O12": ["mmr_O", 'mass mixing ratio of atomic oxygen', 3,
                          'GDZ', 'sph', ['time', 'lon', 'lat', 'height'], ""],
                  "N2": ["mmr_N2_ilev",
                         'mass mixing ratio of molecular nitrogen', 4, 'GDZ',
                         'sph', ['time', 'lon', 'lat', 'ilev'], ""],
                  "N22": ["mmr_N2", 'mass mixing ratio of molecular nitrogen',
                          4, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                          ""],
                  "HE": ["mmr_He_ilev", 'mass mixing ratio of atomic helium',
                         5, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'], ""],
                  "HE2": ["mmr_He", 'mass mixing ratio of atomic helium', 5,
                          'GDZ', 'sph', ['time', 'lon', 'lat', 'height'], ""],
                  "NO": ["mmr_NO_ilev",
                         'mass mixing ratio of molecular nitric oxide', 6,
                         'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'], ""],
                  "NO2": ["mmr_NO",
                          'mass mixing ratio of molecular nitric oxide', 6,
                          'GDZ', 'sph', ['time', 'lon', 'lat', 'height'], ""],
                  "N4S": ["mmr_Nstate4S_ilev",
                          'mass mixing ratio of atomic nitrogen (4S state)',
                          7, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'], ""],
                  "N4S2": ["mmr_Nstate4S",
                           'mass mixing ratio of atomic nitrogen (4S state)',
                           7, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                           ""],
                  "N2D": ["mmr_Nstate2D_ilev",
                          'mass mixing ratio of atomic nitrogen (2D state)',
                          8, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'], ""],
                  "N2D2": ["mmr_Nstate2D",
                           'mass mixing ratio of atomic nitrogen (2D state)',
                           8, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                           ""],
                  "TE": ["T_e_ilev", 'electron temperature', 9, 'GDZ', 'sph',
                         ['time', 'lon', 'lat', 'ilev'], "K"],
                  "TE2": ["T_e", 'electron temperature', 9, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'height'], "K"],
                  "TI": ["T_i_ilev", 'ion temperature', 10, 'GDZ', 'sph',
                         ['time', 'lon', 'lat', 'ilev'], "K"],
                  "TI2": ["T_i", 'ion temperature', 10, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'height'], "K"],
                  "O2P": ["N_O2plus_ilev",
                          'number density of molecular oxygen ion', 11, 'GDZ',
                          'sph', ['time', 'lon', 'lat', 'ilev'], "1/cm**3"],
                  "O2P2": ["N_O2plus",
                           'number density of molecular oxygen ion', 11, 'GDZ',
                           'sph', ['time', 'lon', 'lat', 'height'], "1/cm**3"],
                  "OP": ["N_Oplus_ilev", 'number density of atomic oxygen ion',
                         12, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                         "1/cm**3"],
                  "OP2": ["N_Oplus", 'number density of atomic oxygen ion',
                          12, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                          "1/cm**3"],
                  "N2N": ["N_N2_ilev", 'number density of molecular nitrogen',
                          13, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                          "1/cm**3"],
                  "N2N2": ["N_N2", 'number density of molecular nitrogen',
                           13, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                           "1/cm**3"],
                  "CO2_COOL": ["Q_CO2cool_ilev",
                               'cooling rate of carbon dioxide', 14, 'GDZ',
                               'sph', ['time', 'lon', 'lat', 'ilev'],
                               "erg/g/s"],
                  "CO2_COOL2": ["Q_CO2cool", 'cooling rate of carbon dioxide',
                                14, 'GDZ', 'sph',
                                ['time', 'lon', 'lat', 'height'], "erg/g/s"],
                  "NO_COOL": ["Q_NOcool_ilev", 'cooling rate of nitric oxide',
                              15, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                              "erg/g/s"],
                  "NO_COOL2": ["Q_NOcool", 'cooling rate of nitric oxide',
                               15, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                  'height'], "erg/g/s"],
                  "UN": ["v_neast_ilev", 'zonal neutral wind velocity (east)',
                         16, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                         "cm/s"],
                  "UN2": ["v_neast", 'zonal neutral wind velocity (east)',
                          16, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                          "cm/s"],
                  "VN": ["v_nnorth_ilev",
                         'meridional neutral wind velocity (north)', 17, 'GDZ',
                         'sph', ['time', 'lon', 'lat', 'ilev'], "cm/s"],
                  "VN2": ["v_nnorth",
                          'meridional neutral wind velocity (north)', 17,
                          'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                          "cm/s"],
                  # "O2P_ELD": ['O2P_ELD', '???',18,'GDZ', 'sph',
                  # ['time', 'lon', 'lat', 'ilev'], ''], #NO DESCRIPTION GIVEN
                  # "N2P_ELD": ['N2P_ELD', '???',19,'GDZ', 'sph',
                  # ['time', 'lon', 'lat', 'ilev'], ''], #NO DESCRIPTION GIVEN
                  "NPLUS": ['N_Nplus_ilev',
                            'number density of atomic nitrogen ion', 20, 'GDZ',
                            'sph', ['time', 'lon', 'lat', 'ilev'], '1/cm**3'],
                  "NPLUS2": ['N_Nplus',
                             'number density of atomic nitrogen ion', 20,
                             'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                             '1/cm**3'],
                  # "NOP_ELD": ['NOP_ELD', '???',21,'GDZ', 'sph',
                  # ['time', 'lon', 'lat', 'ilev'], ''], #NO DESCRIPTION GIVEN
                  "SIGMA_PED": ['sigma_P_ilev', 'Pedersen conductivity', 22,
                                'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                                'S/m'],
                  "SIGMA_PED2": ['sigma_P', 'Pedersen conductivity', 22, 'GDZ',
                                 'sph', ['time', 'lon', 'lat', 'height'],
                                 'S/m'],
                  "SIGMA_HAL": ['sigma_H_ilev', 'Hall conductivity', 23, 'GDZ',
                                'sph', ['time', 'lon', 'lat', 'ilev'], 'S/m'],
                  "SIGMA_HAL2": ['sigma_H', 'Hall conductivity', 23, 'GDZ',
                                 'sph', ['time', 'lon', 'lat', 'height'],
                                 'S/m'],
                  "QJOULE": ['Q_Joule_ilev', 'joule heating', 24, 'GDZ', 'sph',
                             ['time', 'lon', 'lat', 'ilev'], 'erg/g/s'],
                  "QJOULE2": ['Q_Joule', 'joule heating', 24, 'GDZ', 'sph',
                              ['time', 'lon', 'lat', 'height'], 'erg/g/s'],
                  "O_N2": ['OtoN2_ilev', 'Oxygen/molecular nitrogen ratio', 25,
                           'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'], ''],
                  "O_N22": ['OtoN2', 'Oxygen/molecular nitrogen ratio', 25,
                            'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                            ''],
                  # "N2D_ELD": ['N2D_ELD', '???',26,'GDZ', 'sph',
                  # ['time', 'lon', 'lat', 'ilev'], ''],  #NO DESCRIPTION GIVEN
                  "O2N": ['N_O2_ilev', 'number density of molecular oxygen',
                          27, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                          '1/cm**3'],
                  "O2N2": ['N_O2', 'number density of molecular oxygen', 27,
                           'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                           '1/cm**3'],

                  # 4D Variables, vert coordinate on interface levels (ilev1)
                  "DEN": ["rho_ilev1", 'total mass density', 28, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'ilev1'], "g/cm**3"],
                  "DEN2": ["rho", 'total mass density', 28, 'GDZ', 'sph',
                           ['time', 'lon', 'lat', 'height'], "g/cm**3"],
                  "ZG": ["H_ilev1",
                         'height dependent on secondary pressure level', 29,
                         'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev1'], "cm"],
                  "Z": ["H_geopot", 'geopotential height', 30, 'GDZ', 'sph',
                        ['time', 'lon', 'lat', 'ilev1'], "cm"],
                  "NE": ["N_e_ilev", 'electron number density', 31, 'GDZ',
                         'sph', ['time', 'lon', 'lat', 'ilev1'], "1/cm**3"],
                  "NE2": ["N_e", 'electron number density', 31, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'height'], "1/cm**3"],
                  "OMEGA": ["omega_ilev", 'Vertical motion frequency', 32,
                            'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev1'],
                            "1/s"],
                  "OMEGA2": ["omega", 'Vertical motion frequency', 32, 'GDZ',
                             'sph', ['time', 'lon', 'lat', 'height'], "1/s"],
                  "POTEN": ["V_ilev", 'electric potential', 33, 'GDZ', 'sph',
                            ['time', 'lon', 'lat', 'ilev1'], "V"],
                  "POTEN2": ["V", 'electric potential', 33, 'GDZ', 'sph',
                             ['time', 'lon', 'lat', 'height'], "V"],
                  "UI_ExB": ["v_iExBeast_ilev",
                             'zonal ExB ion velocity (east)', 34, 'GDZ', 'sph',
                             ['time', 'lon', 'lat', 'ilev1'], 'cm/s'],
                  "UI_ExB2": ["v_iExBeast", 'zonal ExB ion velocity (east)',
                              34, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                 'height'], 'cm/s'],
                  "VI_ExB": ["v_iExBnorth_ilev",
                             'meridional ExB ion velocity (north)', 35, 'GDZ',
                             'sph', ['time', 'lon', 'lat', 'ilev1'], 'cm/s'],
                  "VI_ExB2": ["v_iExBnorth",
                              'meridional ExB ion velocity (north)', 35, 'GDZ',
                              'sph', ['time', 'lon', 'lat', 'height'], 'cm/s'],
                  "WI_ExB": ["v_iExBup_ilev", 'vertical ExB ion velocity (up)',
                             36, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev1'],
                             'cm/s'],
                  "WI_ExB2": ["v_iExBup", 'vertical ExB ion velocity (up)',
                              36, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                 'height'], 'cm/s'],

                  # 4D Variables, vert coord on interface mag levels (imlev)
                  "ZMAG": ["H_milev",
                           'height dependent on geomagnetic pressure level',
                           37, 'MAG', 'sph', ['time', 'mlon', 'mlat', 'milev'],
                           "km"],

                  # 3D Variables,    (time, lat, lon)
                  "TEC": ["TEC", 'vertical total electron content (height ' +
                          'integrated from bottom to top boundary)',
                          38, 'GDZ', 'sph', ['time', 'lon', 'lat'], "1/cm**2"],
                  # Lower boundary condition for TN
                  # "TLBC": ["T_nLBC", 'Lower boundary condition for T_n',39,
                  # 'GDZ', 'sph', ['time', 'lon', 'lat'], "K"],
                  # Lower boundary condition for UN
                  # "ULBC": ["v_neastLBC", 'Lower boundary condition for v_n
                  # east component',40,'GDZ', 'sph', ['time', 'lon', 'lat'],
                  # "cm/s"],
                  # Lower boundary condition for VN
                  # "VLBC": ["v_nnorthLBC", 'Lower boundary condition for v_n
                  # north component',41,'GDZ', 'sph', ['time', 'lon', 'lat'],
                  # "cm/s"],
                  # Lower boundary condition for TN (TIME N-1)
                  # "TLBC_NM": ["T_nLBCNminus1", 'Lower boundary condition for
                  # T_n at t=N-1',42,'GDZ', 'sph', ['time', 'lon', 'lat'],
                  # "K"],
                  # Lower boundary condition for UN (TIME N-1)
                  # "ULBC_NM": ["v_neastLBCNminus1", 'Lower boundary condition
                  # for v_n north component at t=N-1',43,'GDZ', 'sph', ['time',
                  # 'lon', 'lat'], "cm/s"],
                  # Lower boundary condition for VN (TIME N-1)
                  # "VLBC_NM": ["v_nnorthLBCNminus1", 'Lower boundary condition
                  # for v_n east component at t=N-1',44,'GDZ', 'sph', ['time',
                  # 'lon', 'lat'], "cm/s"],
                  "QJOULE_INTEG": ["W_JouleH",
                                   'height integrated joule heating', 45,
                                   'GDZ', 'sph', ['time', 'lon', 'lat'],
                                   'erg/cm**2/s'],
                  "EFLUX": ['Phi_E', 'energy flux', 46, 'GDZ', 'sph',
                            ['time', 'lon', 'lat'], 'erg/cm**2/s'],
                  "HMF2": ['HmF2', 'height of maximum electron number ' +
                           'density in F2 layer', 47, 'GDZ', 'sph',
                           ['time', 'lon', 'lat'], 'km'],
                  "NMF2": ['NmF2', 'maximum electron number density in F2 ' +
                           'layer', 48, 'GDZ', 'sph', ['time', 'lon', 'lat'],
                           '1/cm**3'],
                  }


def dts_to_ts(file_dts):
    '''Get datetime timestamp in UTC from datetime string'''

    return datetime.timestamp(datetime.strptime(file_dts, '%Y-%m-%d %H:%M:%S'
                                                ).replace(tzinfo=timezone.utc))


def year_mtime_todt0(year, mtime):  # self.filedate
    '''Convert year and day to datetime object in UTC at midnight'''

    day, hour, minute = mtime  # unpack mtime values
    return datetime(int(year), 1, 1).replace(tzinfo=timezone.utc) +\
        timedelta(days=int(day-1))


def year_mtime_todt(year, mtime):
    '''Convert year and [day,hour,minute] to datetime object in UTC'''

    day, hour, minute = mtime  # unpack mtime values
    return datetime(int(year), 1, 1).replace(tzinfo=timezone.utc) +\
        timedelta(days=int(day-1), hours=int(hour), minutes=int(minute))


def year_mtime_todts(year, mtime):
    '''Convert year and mtime to a datetime string'''

    return datetime.strftime(year_mtime_todt(year, mtime), '%Y-%m-%d %H:%M:%S')


def year_mtime_todate(year, mtime):
    '''Use year and mtime to determine the date in the file.
    Returns a datetime object.'''

    date_string = datetime.strftime(year_mtime_todt(year, mtime), '%Y-%m-%d')
    return datetime.strptime(date_string, '%Y-%m-%d').replace(
        tzinfo=timezone.utc)


@vectorize
def year_mtime_tohrs(year, day, hour, minute, filedate):
    '''Convert year and mtime to hours since midnight using predetermined
    datetime object.'''

    mtime = [day, hour, minute]
    return (year_mtime_todt(year, mtime) - filedate).total_seconds()/3600.


def ts_to_hrs(time_val, filedate):
    '''Convert utc timestamp to hours since midnight on filedate.'''

    return (datetime.utcfromtimestamp(time_val).replace(tzinfo=timezone.utc) -
            filedate).total_seconds()/3600.


def MODEL():
    from time import perf_counter
    from os.path import basename
    from numpy import zeros, transpose, array, append, insert, where, unique
    from numpy import NaN, diff, abs, mean, broadcast_to, cos, sin, sum
    from numpy import pi as nppi
    from netCDF4 import Dataset
    from kamodo import Kamodo
    import kamodo_ccmc.readers.reader_utilities as RU

    class MODEL(Kamodo):
        '''TIEGCM model data reader.

        Inputs:
            full_filename: a string representing the file name of the model
                output data.
                Note: This reader takes a file name of the format
                file_dir+'*.nc'.
            variables_requested = a list of variable name strings chosen from
                the model_varnames dictionary in this script, specifically the
                first item in the list associated with a given key.
                - If empty, the reader functionalizes all possible variables
                (default)
                - If 'all', the reader returns the model_varnames dictionary
                    above for only the variables present in the given files.
                    Note: the fulltime keyword must be False to acheive this
                    behavior.
            filetime = boolean (default = False)
                - if False, the script fully executes.
                - If True, the script only executes far enough to determine the
                    time values associated with the chosen data.
                Note: The behavior of the script is determined jointly by the
                    filetime and fulltime keyword values.
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
            fulltime = boolean (default = True)
                - If True, linear interpolation in time between files is
                    included in the returned interpolator functions.
                - If False, no linear interpolation in time between files is
                    included.
            verbose = boolean (False)
                - If False, script execution and the underlying Kamodo
                    execution is quiet except for specified messages.
                - If True, be prepared for a plethora of messages.
        All inputs are described in further detail in
            KamodoOnboardingInstructions.pdf.

        Returns: a kamodo object (see Kamodo core documentation) containing all
            requested variables in functionalized form.
        '''
        def __init__(self, full_filename, variables_requested=[],
                     filetime=False, verbose=False, gridded_int=True,
                     printfiles=False, fulltime=True, **kwargs):
            super(MODEL, self).__init__()
            self.modelname = 'TIEGCM'

            # store time information
            t0 = perf_counter()
            filename = basename(full_filename)
            file_dir = full_filename.split(filename)[0]
            cdf_data = Dataset(full_filename, 'r')
            year = array(cdf_data.variables['year'])
            mtime = array(cdf_data.variables['mtime'])
            day, hour, minute = mtime.T
            # datetime object for file date at midnight UTC
            self.filedate = year_mtime_todt0(year[0], mtime[0])
            # strings in format = YYYY-MM-DD HH:MM:SS
            self.datetimes = [year_mtime_todts(y, m) for y, m
                              in zip([year[0], year[-1]], [mtime[0],
                                                           mtime[-1]])]
            # timestamps in UTC
            self.filetimes = [dts_to_ts(file_dts) for file_dts in
                              self.datetimes]
            time = year_mtime_tohrs(year, day, hour, minute, self.filedate)
            if len(time) > 1:
                self.dt = diff(time).max()*3600.  # time is in hours
            else:
                self.dt = 0.

            if filetime and not fulltime:
                return  # return times as is to prevent recursion

            # if variables are given as integers, convert to standard names
            if len(variables_requested) > 0:
                if isinstance(variables_requested[0], int):
                    tmp_var = [value[0] for key, value in
                               model_varnames.items() if value[2] in
                               variables_requested]
                    variables_requested = tmp_var

            if fulltime:  # add boundary time (default value)
                # find other files with same pattern
                from glob import glob

                file_pattern = file_dir+'*.nc'  # returns a string for tiegcm
                files = sorted(glob(file_pattern))
                filenames = unique([basename(f) for f in files])

                # find closest file by utc timestamp
                # tiegcm has an open time at the beginning
                current_idx = where(filenames == filename)[0]
                if current_idx == 0:
                    filecheck = False
                    if verbose:
                        print('No earlier file available.')
                    if filetime:
                        return
                else:
                    min_filename = file_dir+filenames[current_idx-1][0]
                    kamodo_test = MODEL(min_filename, filetime=True,
                                        fulltime=False)
                    time_test = abs(kamodo_test.filetimes[1]-self.filetimes[0])
                    if time_test <= self.dt or (self.dt == 0. and
                                                time_test <= 3600.*6.):
                        # if nearest file time at least within one timestep
                        filecheck = True
                        # add time at beginning
                        self.datetimes[0] = kamodo_test.datetimes[1]
                        self.filetimes[0] = kamodo_test.filetimes[1]

                        # time only version if returning time
                        if filetime:
                            return

                        # get kamodo object with same requested variables
                        if verbose:
                            print(f'Took {perf_counter()-t0:.3f}s to find ' +
                                  'closest file.')
                        kamodo_neighbor = MODEL(
                            min_filename,
                            variables_requested=variables_requested,
                            fulltime=False)
                        short_data = kamodo_neighbor.short_data
                        if verbose:
                            print(f'Took {perf_counter()-t0:.3f}s to get ' +
                                  'data from previous file.')
                    else:
                        if verbose:
                            print('No earlier file found within ' +
                                  f'{self.dt:.1f}s')
                        filecheck = False
                        if filetime:
                            return

            # These lists are the standardized variable name.
            # The only milev variable is H_milev.
            ilev1_list = [value[0] for key, value in model_varnames.items()
                          if value[5][-1] == 'ilev1']
            ilev1_replace = [item.split('_ilev1')[0] for item in ilev1_list if
                             item != 'H_ilev1']
            ilev_list = [value[0] for key, value in model_varnames.items()
                         if value[5][-1] == 'ilev']
            ilev_replace = [item.split('_ilev')[0] for item in ilev_list if
                            item != 'H_ilev']
            # milev_list = [value[0] for key, value in model_varnames.items()
            #               if value[5][-1] == 'milev']
            self.total_ilev = [item for item in ilev_list + ilev1_list if item
                               not in ['H_ilev', 'H_ilev1']]
            self.total_replace = ilev_replace + ilev1_replace
            # dictionary mapping to navigate related variable names
            self.ilev_map = {item1: item2 for item1, item2 in
                             zip(self.total_replace, self.total_ilev)}

            # perform initial check on variables_requested list
            if len(variables_requested) > 0 and fulltime and (
                    variables_requested != 'all'):
                test_list = [value[0] for key, value in model_varnames.items()]
                err_list = [item for item in variables_requested if item not in
                            test_list]
                if len(err_list) > 0:
                    print('Variable name(s) not recognized:', err_list)

            # translate from standardized variables to names in file
            # remove variables requested that are not in the file
            if len(variables_requested) > 0 and variables_requested != 'all':
                # add ilev version of variables to the list, adding H_ilev(1)
                var_replace = [var for var in variables_requested if var in
                               ilev_replace or var in ilev1_replace]
                add_ilev = [var+'_ilev' for var in variables_requested if var
                            in ilev_replace]
                add_ilev1 = [var+'_ilev1' for var in variables_requested if var
                             in ilev1_replace]
                if len(add_ilev1 + add_ilev) > 0:
                    print('The following functions will be converted from ' +
                          'pressure level to height and will not have ' +
                          'gridded versions:', var_replace)
                if len(add_ilev) > 0:
                    add_ilev += ['H_ilev']
                if len(add_ilev1) > 0:
                    add_ilev1 += ['H_ilev1']
                new_var = variables_requested + add_ilev + add_ilev1
                short_var = [item for item in new_var if item not
                             in self.ilev_map.keys()]  # remove replaced items
                gvar_list = [key for key, value in model_varnames.items()
                             if value[0] in short_var and key in
                             cdf_data.variables.keys()]  # file variable names

                # check for variables requested but not available
                if len(gvar_list) != len(short_var):
                    err_list = [value[0] for key, value in
                                model_varnames.items() if value[0] in
                                short_var and key not in
                                cdf_data.variables.keys()]
                    if len(err_list) > 0:
                        print('Some requested variables are not available: ',
                              err_list)
                    if 'H_ilev' in err_list or 'H_ilev1' in err_list:
                        other_name = ['H_ilev', 'H_ilev1']
                        err_name = [name for name in other_name if name in
                                    err_list][0]
                        other_name.remove(err_name)  # see first element
                        gvar_name = [key for key, value in
                                     model_varnames.items() if value[0] ==
                                     other_name[0]]
                        gvar_list += gvar_name  # add other height name to list
                        print(f'Retrieving the {other_name[0]} variable' +
                              ' instead.')

            else:  # only input variables on the avoid_list if requested
                avoid_list = ['TLBC', 'ULBC', 'VLBC', 'TLBC_NM', 'ULBC_NM',
                              'VLBC_NM', 'NOP_ELD', 'O2P_ELD', 'N2P_ELD',
                              'N2D_ELD']
                gvar_list = [key for key in cdf_data.variables.keys()
                             if key in model_varnames.keys() and
                             key not in avoid_list]
                print('The following functions will be converted from ' +
                      'pressure level to height and will not have ' +
                      'gridded versions:', self.total_replace)
                if not fulltime and variables_requested == 'all':
                    self.var_dict = {value[0]: value[1:] for key, value in
                                     model_varnames.items() if key in
                                     gvar_list}
                    # add non-ilev versions of the variables in the files
                    key_list = list(self.var_dict.keys())
                    for var_key in key_list:
                        if var_key in self.total_ilev:
                            # retrieve equivalent non-ilev variable name
                            new_key = [key for key, value in
                                       self.ilev_map.items() if value ==
                                       var_key][0]
                            # retrieve key for model_varnames mapping
                            model_key = [key for key, value in
                                         model_varnames.items() if
                                         value[0] == new_key][0]
                            # add to returned dictionary
                            self.var_dict[model_varnames[model_key][0]] =\
                                model_varnames[model_key][1:]

                    cdf_data.close()
                    return

            # Store the requested variables into a dictionary
            variables = {model_varnames[key][0]: {
                'units': model_varnames[key][-1],
                'data': array(cdf_data.variables[key])} for key in gvar_list}

            # prepare and return data only for last timestamp
            if not fulltime:
                cdf_data.close()
                variables['time'] = self.filetimes[1]  # utc timestamp
                self.short_data = variables
                return

            # Store inputs as class attributes
            self.filename = full_filename
            self.missing_value = NaN
            self._registered = 0
            self.variables = dict()
            self.modelname = 'TIEGCM'
            if printfiles:
                print('Files:', self.filename)

            # Store new time if neighboring file found.
            if filecheck:  # new_time is a utc timestamp
                new_time = ts_to_hrs(short_data['time'], self.filedate)
                self._time = insert(time, 0, new_time)
            else:
                self._time = time

            # store coordinates
            lat = array(cdf_data.variables['lat'])  # NOT FULL RANGE IN LAT
            lat = insert(lat, 0, -90)  # insert a grid point before -87.5
            self._lat = append(lat, 90.)   # and at the end (after 87.5)
            lon = array(cdf_data.variables['lon'])  # NOT WRAPPED IN LONGITUDE
            self._lon = append(lon, 180.)  # add 180. to end of array
            self._ilev = array(cdf_data.variables['lev'])
            self._ilev1 = array(cdf_data.variables['ilev'])
            self._milev = array(cdf_data.variables['imlev'])
            self._mlat = array(cdf_data.variables['mlat'])
            self._mlon = array(cdf_data.variables['mlon'])  # -180 to 180
            cdf_data.close()
            if verbose:
                print(f'Took {perf_counter()-t0:.6f}s to read in data')

            # Check for presence of necessary height variables in varname_list.
            varname_list = [key for key in variables.keys()]
            ilev1_check = unique([True for item in varname_list if 'ilev1' ==
                                  item[-5:]])
            ilev_check = unique([True for item in varname_list if 'ilev' ==
                                 item[-4:]])
            if ilev1_check and 'H_ilev1' not in varname_list:
                self.ilev_sub = True
            elif ilev_check and 'H_ilev' not in varname_list:
                self.ilev_sub = True
            else:
                self.ilev_sub = False

            # register interpolators for each requested variable
            # rearrange to deal with H_ilev and H_ilev1 first if there
            if 'H_ilev' in varname_list:
                varname_list.remove('H_ilev')
                varname_list = ['H_ilev'] + varname_list
            if 'H_ilev1' in varname_list:
                varname_list.remove('H_ilev1')
                varname_list = ['H_ilev1'] + varname_list
            self.variables = {}
            t_reg = perf_counter()
            for varname in varname_list:
                if len(variables[varname]['data'].shape) == 3:
                    if filecheck:  # if neighbor found
                        # append data for first time stamp, transpose
                        data_shape = list(variables[varname]['data'].shape)
                        data_shape[0] += 1  # add space for time
                        new_data = zeros(data_shape)
                        # put in current data
                        new_data[1:, :, :] = variables[varname]['data']
                        # add in data for additional time
                        new_data[0, :, :] = \
                            short_data[varname]['data'][-1, :, :]
                        # (t,lat,lon) -> (t,lon,lat)
                        variable = transpose(new_data, (0, 2, 1))
                    else:
                        # (t,lat,lon) -> (t,lon,lat)
                        variable = transpose(variables[varname]['data'],
                                             (0, 2, 1))
                    self.variables[varname] = dict(
                        units=variables[varname]['units'], data=variable)
                    self.register_3D_variable(self.variables[varname]['units'],
                                              self.variables[varname]['data'],
                                              varname, gridded_int)
                elif len(variables[varname]['data'].shape) == 4:
                    if filecheck:
                        # append data for first time stamp, transpose
                        data_shape = list(variables[varname]['data'].shape)
                        data_shape[0] += 1  # add space for time
                        new_data = zeros(data_shape)
                        # put in current data
                        new_data[1:, :, :, :] = variables[varname]['data']
                        # add in data for additional time
                        new_data[0, :, :, :] = \
                            short_data[varname]['data'][-1, :, :, :]
                        # (t,h,lat,lon) -> (t,lon,lat,h)
                        variable = transpose(new_data, (0, 3, 2, 1))
                    else:
                        # (t,h,lat,lon) -> (t,lon,lat,h)
                        variable = transpose(variables[varname]['data'],
                                             (0, 3, 2, 1))
                    self.variables[varname] = dict(
                        units=variables[varname]['units'], data=variable)
                    self.register_4D_variable(self.variables[varname]['units'],
                                              self.variables[varname]['data'],
                                              varname, gridded_int,
                                              verbose=verbose)
            if verbose:
                print(f'Took {perf_counter()-t_reg:.5f}s to register ' +
                      f'{len(varname_list)} variables.')
            if verbose:
                print(f'Took a total of {perf_counter()-t0:.5f}s to kamodofy' +
                      f' {len(varname_list)} variables.')

        def wrap_3Dlatlon(self, varname, variable):
            '''Wraps the data array in longitude (-180=180), and latitude'''

            shape_list = list(variable.shape)
            shape_list[2] += 2  # need two more places in latitude
            shape_list[1] += 1  # need one more place in longitude
            tmp_arr = zeros(shape_list)  # array to set-up wrapped data in
            tmp_arr[0:, :-1, 1:-1] = variable  # copy data into grid

            # wrapping in latitude for scalar variables
            # put in top values
            top = mean(tmp_arr[:, :-1, 1], axis=1)  # same shape as time axis
            tmp_arr[:, :-1, 0] = broadcast_to(top, (shape_list[1]-1,
                                                    shape_list[0])).T
            # same for bottom, reusing variable names
            top = mean(tmp_arr[:, :-1, -2], axis=1)  # same shape as time axis
            tmp_arr[:, :-1, -1] = broadcast_to(top, (shape_list[1]-1,
                                                     shape_list[0])).T

            # wrap in longitude after to prevent double counting in average
            tmp_arr[:, -1, :] = tmp_arr[:, 0, :]
            self.variables[varname]['data'] = tmp_arr  # store result
            return tmp_arr

        def wrap_4Dlatlon(self, varname, variable):
            '''Wraps the data array in longitude (-180=180), and latitude
            (0=-2, -1=1)'''

            shape_list = list(variable.shape)  # time, lon, lat, ilev
            shape_list[2] += 2  # need two more places in latitude
            shape_list[1] += 1  # need one more place in longitude
            tmp_arr = zeros(shape_list)  # array to set-up wrapped data in
            tmp_arr[:, :-1, 1:-1, :] = variable  # copy data into grid

            # wrapping in latitude for scalar and cartesian/radial variables
            if varname not in ['u_n', 'v_n', 'u_iExB', 'v_iExB']:
                # put in top values
                top = mean(tmp_arr[:, :-1, 1, :], axis=1)  # average over lon
                tmp_arr[:, :-1, 0, :] = transpose(broadcast_to(
                    top, (shape_list[1]-1, shape_list[0], shape_list[3])),
                    (1, 0, 2))
                # same for bottom, reusing variable names
                top = mean(tmp_arr[:, :-1, -2, :], axis=1)  # average over lon
                tmp_arr[:, :-1, -1, :] = transpose(broadcast_to(
                    top, (shape_list[1]-1, shape_list[0], shape_list[3])),
                    (1, 0, 2))
            # wrapping in latitude for relevant vector variables
            elif varname in ['u_n', 'v_n', 'u_iExB', 'v_iExB']:
                # calculate net vector magnitude for top
                tmp_arr[:, :-1, 0, :] = self.vector_average4D(
                    tmp_arr[:, :-1, 1, :], shape_list, varname, self._lat[0])
                # repeat for bottom
                tmp_arr[:, :-1, -1, :] = self.vector_average4D(
                    tmp_arr[:, :-1, -2, :], shape_list, varname, self._lat[-1])
            tmp_arr[:, -1, :, :] = tmp_arr[:, 0, :, :]  # wrap value in lon
            self.variables[varname]['data'] = tmp_arr  # store result
            return tmp_arr

        def vector_average4D(self, top, shape_list, varname, latval):
            '''find vector average at pole for array with shape (time, lon,
            height)'''

            # find net x and y components, final array shapes are (t, lon, ht)
            lon_arr = transpose(broadcast_to(self._lon[:-1],
                                             (shape_list[0], shape_list[3],
                                              shape_list[1]-1)), (0, 2, 1))
            # need to put 'old' shape at end in broadcast_to call ^
            xval = sum(top*cos((lon_arr+180.)*nppi/180.), axis=1)  # same shape
            yval = sum(top*sin((lon_arr+180.)*nppi/180.), axis=1)  # as t and z
            xarr = transpose(broadcast_to(xval, (shape_list[1]-1,
                                                 shape_list[0],
                                                 shape_list[3])), (1, 0, 2))
            yarr = transpose(broadcast_to(yval, (shape_list[1]-1,
                                                 shape_list[0],
                                                 shape_list[3])), (1, 0, 2))

            # convert to proper unit vector (see wiki on spherical coordinates)
            if 'u' in varname:
                # Zonal/east components -> convert to psi_hat vector (lon)
                # -xsin(psi)+ycos(psi), psi = longitude (0 to 360)
                new_top = -xarr*sin((lon_arr+180.)*nppi/180.) +\
                    yarr*cos((lon_arr+180.)*nppi/180.)
            elif 'v' in varname:
                # meridional/north -> convert to theta_hat vector (latitude)
                # xcos(psi)cos(theta)+ysin(psi)cos(theta)
                # sin(theta) is always zero at the poles
                # theta = latitude (0 to 180), psi = longitude (0 to 360)
                new_top = xarr * cos((lon_arr+180.) * nppi/180.) *\
                    cos((90.-latval)*nppi/180.) + yarr *\
                    sin((lon_arr+180.) * nppi/180.) * cos((90.-latval) *
                                                          nppi/180.)

            # flip lon around so values match lon location (keep zero ref true)
            zero_idx = min(where(self._lon >= 0.)[0])  # find splitting index
            # create empty array for destination
            top = zeros((shape_list[0], shape_list[1]-1, shape_list[3]))
            div = (shape_list[1]-1) % 2  # deal with even or odd number of lon
            # move last half to first half and vice versa
            top[:, :zero_idx-div, :] = new_top[:, zero_idx:, :]
            top[:, zero_idx-div:, :] = new_top[:, :zero_idx, :]
            return top

        # Define and register a 3D variable --------------------------
        def register_3D_variable(self, units, variable, varname, gridded_int):
            """Registers a 3d interpolator with 3d signature"""

            # define and register the interpolators
            xvec_dependencies = {'time': 'hr', 'lon': 'deg', 'lat': 'deg'}
            wrapped_data = self.wrap_3Dlatlon(varname, variable)
            self = RU.regdef_3D_interpolators(self, units, wrapped_data,
                                              self._time, self._lon, self._lat,
                                              varname, xvec_dependencies,
                                              gridded_int)
            return

        # Define and register a 4D variable -------------------------------
        def register_4D_variable(self, units, variable, varname, gridded_int,
                                 verbose=False):
            """Registers a 4d interpolator with 4d signature"""

            # Get the correct coordinates
            coord_list = [value[5] for key, value in model_varnames.items()
                          if value[0] == varname][0]
            if 'ilev1' in coord_list and hasattr(self, '_ilev1'):
                h = self._ilev1
                coord_lat, coord_lon = self._lat, self._lon
                xvec_dependencies = {'time': 'hr', 'lon': 'deg', 'lat': 'deg',
                                     'ilev1': 'm/m'}
            elif 'ilev' in coord_list and hasattr(self, '_ilev'):
                h = self._ilev
                coord_lat, coord_lon = self._lat, self._lon
                xvec_dependencies = {'time': 'hr', 'lon': 'deg', 'lat': 'deg',
                                     'ilev': 'm/m'}
            elif 'milev' in coord_list and hasattr(self, '_milev'):
                h = self._milev
                coord_lat, coord_lon = self._mlat, self._mlon
                xvec_dependencies = {'time': 'hr', 'mlon': 'deg',
                                     'mlat': 'deg', 'milev': 'm/m'}

            # define and register the interpolators
            if 'lat' in xvec_dependencies.keys():
                wrapped_data = self.wrap_4Dlatlon(varname, variable)
            else:
                top_shape = list(variable[:, :, :, -1].shape)
                top_size = top_shape[0] * top_shape[1] * top_shape[2]  # 3D arr
                idx_top = where(variable[:, :, :, -1] > 1e+35)[0]
                tmp_data = variable
                while top_size == len(idx_top):   # remove undefined top row(s)
                    if verbose:
                        print('All values at max milev are 1e+36 for ' +
                              f'{varname}. Slicing off top array.')
                    if self._milev.shape[0] == len(tmp_data[0, 0, 0, :]):
                        self._milev = self._milev[0:-1]
                    tmp_data = tmp_data[:, :, :, 0:-1]
                    top_shape = list(tmp_data[:, :, :, -1].shape)
                    top_size = top_shape[0]*top_shape[1]*top_shape[2]  # 3D
                    idx_top = where(tmp_data[:, :, :, -1] > 1e+35)[0]
                wrapped_data = tmp_data
                h = self._milev
            self = RU.regdef_4D_interpolators(self, units, wrapped_data,
                                              self._time, coord_lon, coord_lat,
                                              h, varname, xvec_dependencies,
                                              gridded_int)

            # perform substitution if needed
            if varname in ['H_ilev', 'H_ilev1'] and self.ilev_sub:
                other_name = ['H_ilev', 'H_ilev1']
                other_name.remove(varname)  # first element is the other name
                print(f'{other_name[0]} missing in data and is needed to ' +
                      'convert the requested variables to depend on height.' +
                      f' Using {varname} instead.')
                self.variables[other_name[0]] = dict(units=units,
                                                     data=variable)
                # register the other variable by linking to this one
                self = RU.register_interpolator(self, other_name[0], varname,
                                                xvec_dependencies)

            # implement custom interpolation method here for f(ilev)
            if varname in self.total_ilev:
                # set inputs
                new_xvec_dependencies = {'time': 'hr', 'lon': 'deg',
                                         'lat': 'deg', 'height': 'km'}
                new_varname = varname.split('_ilev')[0]
                # Import and call custom interpolator
                from tiegcm_ilevinterp import PLevelInterp
                interpolator = PLevelInterp(self, h, 'H_'+coord_list[-1])
                # Register in kamodo object
                self.variables['P'+coord_list[-1][1:]] = dict(units='m/m')
                self = RU.register_interpolator(self,
                                                'P'+coord_list[-1][1:],
                                                interpolator,
                                                xvec_dependencies)

                # Use function composition to convert variable function
                comp_string = varname+'(P'+coord_list[-1][1:]+')'
                self[new_varname] = comp_string
                self.variables[new_varname] = dict(units=units)
                self.variables[new_varname]['xvec'] = new_xvec_dependencies
                self._registered += 1
            return

    return MODEL
