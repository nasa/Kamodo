'''
TIEGCM Kamodo reader, adapted to new structure for satellite flythrough
    software
Initial version - Asher Pembroke (?)
Initial version of model_varnames contributed by Zachary Waldron
New code: Rebecca Ringuette (June 2021 and on)

The high-altitude version of TIEGCM has one timestep per file, while the normal
version of TIEGCM has multiple timesteps per file. This difference requires
different interp_flag values and different functions. Adding the logic for
this at the end. - Nov 2, 2022

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
                  "NE": ["N_e_ilev1", 'electron number density', 31, 'GDZ',
                         'sph', ['time', 'lon', 'lat', 'ilev1'], "1/cm**3"],
                  "NE2": ["N_e", 'electron number density', 31, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'height'], "1/cm**3"],
                  "OMEGA": ["omega_ilev1", 'Vertical motion frequency', 32,
                            'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev1'],
                            "1/s"],
                  "OMEGA2": ["omega", 'Vertical motion frequency', 32, 'GDZ',
                             'sph', ['time', 'lon', 'lat', 'height'], "1/s"],
                  "POTEN": ["V_ilev1", 'electric potential', 33, 'GDZ', 'sph',
                            ['time', 'lon', 'lat', 'ilev1'], "V"],
                  "POTEN2": ["V", 'electric potential', 33, 'GDZ', 'sph',
                             ['time', 'lon', 'lat', 'height'], "V"],
                  "UI_ExB": ["v_iExBeast_ilev1",
                             'zonal ExB ion velocity (east)', 34, 'GDZ', 'sph',
                             ['time', 'lon', 'lat', 'ilev1'], 'cm/s'],
                  "UI_ExB2": ["v_iExBeast", 'zonal ExB ion velocity (east)',
                              34, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                 'height'], 'cm/s'],
                  "VI_ExB": ["v_iExBnorth_ilev1",
                             'meridional ExB ion velocity (north)', 35, 'GDZ',
                             'sph', ['time', 'lon', 'lat', 'ilev1'], 'cm/s'],
                  "VI_ExB2": ["v_iExBnorth",
                              'meridional ExB ion velocity (north)', 35, 'GDZ',
                              'sph', ['time', 'lon', 'lat', 'height'], 'cm/s'],
                  "WI_ExB": ["v_iExBup_ilev1", 'vertical ExB ion velocity (up)',
                             36, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev1'],
                             'cm/s'],
                  "WI_ExB2": ["v_iExBup", 'vertical ExB ion velocity (up)',
                              36, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                 'height'], 'cm/s'],
                  "EX": ["E_east_ilev1", 'zonal component of electric field',
                             71, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev1'],
                             'V/m'],
                  "EX2": ["E_east", 'zonal component of electric field',
                              72, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                 'height'], 'V/m'],
                  "EY": ["E_north_ilev1", 'meridional component of electric field',
                             73, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev1'],
                             'V/m'],
                  "EY2": ["E_north", 'meridional component of electric field',
                              74, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                 'height'], 'V/m'],
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
                  "ALFA": ["E_Char", 'Aurora Characteristic Energy', 70, 'GDZ', 'sph',
                            ['time', 'lon', 'lat'], 'keV'],
                  }


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


@vectorize
def year_mtime_tohrs(year, day, hour, minute, filedate):
    '''Convert year and mtime to hours since midnight using predetermined
    datetime object.'''

    mtime = [day, hour, minute]
    return (year_mtime_todt(year, mtime) - filedate).total_seconds()/3600.


def MODEL():
    from time import perf_counter
    from os.path import basename
    from numpy import zeros, transpose, array, append, insert, where, unique
    from numpy import NaN, mean, broadcast_to, cos, sin, sum, squeeze
    from numpy import pi as nppi
    from kamodo import Kamodo
    import kamodo_ccmc.readers.reader_utilities as RU

    from scipy.interpolate import RegularGridInterpolator as rgiND
    from numpy import log, exp
    
    class MODEL(Kamodo):
        '''TIEGCM model data reader.

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
            - TIE-GCM model outputs are given in netCDF files, so no file
              conversion is needed. However, all of the variables depend on a
              pressure level coordinate, so some pre-processing is needed to
              calculate the median height across the entire dataset
              (time, lon, lat) for a given pressure level value for each type
              of pressure level coordinate (up to three).
            - TIE-GCM data is given in several coordinate systems depending
              on pressure level - a unitless representation of height
              corresponding to a given atmospheric pressure. The preset values
              of pressure level in the coordinate systems are not
              guaranteed to be identical, so they are inverted independently of
              the other unless only one is provided. In that case, the are
              assumed to be identical. The magnetic pressure level coordinate
              is always treated independently of the other two.
            - Pressure level inversion is performed in the reader_utilities
              script, specifically the PLevelInterp function. Two versions of
              all variables that depend on pressure level are created: the
              original and one dependent on height, which is created through
              Kamodo's function composition feature.
            - The outputs do not provide values at the poles, so scalar and
              vector averaging are used as appropriate to determine these
              values on the fly.
            - The files are small, but different versions require different
              interpolation choices. The high-altitude version of TIEGCM has
              one timestep per file, requiring that interpolation methods 1 be
              chosen, while the normal version of TIEGCM has multiple timesteps
              per file, requiring that interpolation method 2 be chosen. This
              difference requires custom code to automatically choose the
              appropriate interpolation method and the associated read logic.

        Developer note:
            The current logic for variables that depend on imlev slices off the
            e36 values in self._imlev coordinate array. This only works because
            there is only one variable that depends on imlev: H_imlev. The
            logic beginning on line 767 will have to be reworked a bit if any
            new variables depend on imlev to avoid repeatedly slicing off the
            coordinate values.
        '''
        def __init__(self, file_dir, variables_requested=[],
                     filetime=False, verbose=False, gridded_int=True,
                     printfiles=False, **kwargs):
            super(MODEL, self).__init__()
            self.modelname = 'TIEGCM'
            t0 = perf_counter()

            # first, check for file list, create if DNE
            list_file = file_dir + self.modelname + '_list.txt'
            time_file = file_dir + self.modelname + '_times.txt'
            self.times, self.pattern_files = {}, {}
            if not RU._isfile(list_file) or not RU._isfile(time_file):
                # collect filenames, all one pattern, so doesn't matter
                raw_files = sorted(RU.glob(file_dir+'*.nc'))
                # ignore TIEGCM_km.nc, pxxx.nc, and runname_tie_ files
                # This will not work if TIEGCM is in the runname
                files = [f for f in raw_files if self.modelname not in
                         basename(f) and 
                         (len(basename(f))==7 and basename(f)[0]=='s'
                          or '_sech_tie_' in basename(f))]
                if len(basename(files[0])) == 7:  # e.g. s001.nc
                    patterns = unique([basename(f)[0] for f in files])
                else:  # e.g. the high altitude outputs
                    patterns = unique([basename(f)[:-22] for f in files])
                self.filename = ''.join([f+',' for f in files])[:-1]

                # establish time attributes
                for p in patterns:
                    # get list of files to loop through later
                    pattern_files = sorted(RU.glob(file_dir+p+'*.nc'))
                    self.pattern_files[p] = pattern_files
                    self.times[p] = {'start': [], 'end': [], 'all': [], 'start_index': []}

                    # loop through to get times
                    for f in range(len(pattern_files)):
                        print(pattern_files[f])
                        ncdf_filetype = RU.netcdf_type(pattern_files[f])
                        cdf_data = RU.Dataset(pattern_files[f],filetype=ncdf_filetype)
                        year = array(cdf_data['year'])
                        mtime = array(cdf_data['mtime'])
                        day, hour, minute = mtime.T
                        # datetime object for file date at midnight UTC
                        if f == 0 and p == patterns[0]:
                            self.filedate = year_mtime_todt0(year[0], mtime[0])
                        time = year_mtime_tohrs(year, day, hour, minute,
                                                self.filedate)
                        cdf_data.close()
                        self.times[p]['start_index'].append(len(self.times[p]['all']))
                        self.times[p]['start'].append(time[0])
                        self.times[p]['end'].append(time[-1])
                        self.times[p]['all'].extend(time)
                    self.times[p]['start_index'].append(len(self.times[p]['all'])-1)

                    self.times[p]['start_index'] = array(self.times[p]['start_index'])
                    self.times[p]['start'] = array(self.times[p]['start'])
                    self.times[p]['end'] = array(self.times[p]['end'])
                    self.times[p]['all'] = array(self.times[p]['all'])

                # perform post-processing to speed up pressure level conversion
                from kamodo_ccmc.readers.tiegcm_tocdf import convert_all
                convert_all(file_dir, self.pattern_files, self.times)

                # create time list file if DNE
                RU.create_timelist(list_file, time_file, self.modelname,
                                   self.times, self.pattern_files,
                                   self.filedate)
            else:  # read in data and time grids from file list
                self.times, self.pattern_files, self.filedate, self.filename =\
                    RU.read_timelist(time_file, list_file)
            if filetime:
                return  # return times as is to prevent recursion

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

            # translate from standardized variables to names in file
            # remove variables requested that are not in the file
            p = list(self.pattern_files.keys())[0]
            pattern_files = self.pattern_files[p]
            ncdf_filetype = RU.netcdf_type(pattern_files[0])
            print('file type: ',ncdf_filetype)
            cdf_data = RU.Dataset(pattern_files[0],filetype=ncdf_filetype)
            if len(variables_requested) > 0 and variables_requested != 'all':
                # add ilev version of variables to the list, adding H_ilev(1)
                add_ilev = [var+'_ilev' for var in variables_requested if var
                            in ilev_replace]
                add_ilev1 = [var+'_ilev1' for var in variables_requested if var
                             in ilev1_replace]
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
                if variables_requested == 'all':
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
            self.variables = {model_varnames[key][0]: {
                'units': model_varnames[key][-1],
                'data': key} for key in gvar_list}

            # Store inputs as class attributes
            self.missing_value = NaN
            self._registered = 0
            if printfiles:
                print('Files:', self.filename)

            # store coordinates
            lat = array(cdf_data['lat'])  # NOT FULL RANGE IN LAT
            lat = insert(lat, 0, -90)  # insert a grid point before -87.5
            self._lat = append(lat, 90.)   # and at the end (after 87.5)
            lon = array(cdf_data['lon'])  # NOT WRAPPED IN LONGITUDE
            self._lon = append(lon, 180.)  # add 180. to end of array
            self._ilev = array(cdf_data['lev'])
            self._ilev1 = array(cdf_data['ilev'])
            self._milev = array(cdf_data['imlev'])
            self._mlat = array(cdf_data['mlat'])
            self._mlon = array(cdf_data['mlon'])  # -180 to 180
            if verbose:
                print(f'Took {perf_counter()-t0:.6f}s to read in data')

            # Check for presence of necessary height variables in varname_list.
            varname_list = [key for key in self.variables.keys()]
            ilev1_check = unique([True for item in varname_list if 'ilev1' ==
                                  item[-5:]])
            ilev_check = unique([True for item in varname_list if 'ilev' ==
                                 item[-4:]])
            if ilev1_check and 'H_ilev1' not in varname_list:
                self.ilev_sub = 'H_ilev1'  # name of H missing
            elif ilev_check and 'H_ilev' not in varname_list:
                self.ilev_sub = 'H_ilev'
            else:
                self.ilev_sub = False

            # register interpolators for each requested variable
            # rearrange to deal with H_ilev and H_ilev1 first if there
            # also calcalate median km grids
            if 'H_ilev' in varname_list:  # height in cm
                varname_list.remove('H_ilev')
                varname_list = ['H_ilev'] + varname_list
                if RU._isfile(file_dir+'TIEGCM_km.nc'):  # km_ilev from file
                    km_data = RU.Dataset(file_dir+'TIEGCM_km.nc') # ,filetype='netCDF3')
                    if hasattr(km_data, 'km_ilev_max'):
                        self._km_ilev = array(km_data['km_ilev'])
                        self._km_ilev_max = km_data.km_ilev_max
                        self._km_ilev_min = km_data.km_ilev_min
                    km_data.close()
            if 'H_ilev1' in varname_list:  # height in cm
                varname_list.remove('H_ilev1')
                varname_list = ['H_ilev1'] + varname_list
                if RU._isfile(file_dir+'TIEGCM_km.nc'):  # km_ilev1 from file
                    km_data = RU.Dataset(file_dir+'TIEGCM_km.nc')  #,filetype='netCDF3')
                    if hasattr(km_data, 'km_ilev1_max'):
                        self._km_ilev1 = array(km_data['km_ilev1'])
                        self._km_ilev1_max = km_data.km_ilev1_max
                        self._km_ilev1_min = km_data.km_ilev1_min
                    km_data.close()
            cdf_data.close()
            t_reg = perf_counter()
            for varname in varname_list:
                self.register_variable(varname, gridded_int, verbose)
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

        def register_variable(self, varname, gridded_int, verbose=False):
            '''register the variable data with the chunked interpolation method
            since the data files are small and chunked.'''

            # determine coordinate variables and xvec by coord list
            gvar = self.variables[varname]['data']
            p = list(self.pattern_files.keys())[0]
            coord_list = [value[5] for key, value in model_varnames.items()
                          if value[0] == varname][0]
            # time grids are the same in both hemispheres
            coord_dict = {'time': {'units': 'hr',
                                   'data': self.times[p]['all']}}
            if 'lon' in coord_list:
                coord_dict['lon'] = {'units': 'deg', 'data': self._lon}
                coord_dict['lat'] = {'units': 'deg', 'data': self._lat}
            if 'milev' in coord_list and hasattr(self, '_milev'):
                coord_dict['mlon'] = {'units': 'deg', 'data': self._mlon}
                coord_dict['mlat'] = {'units': 'deg', 'data': self._mlat}
                coord_dict['milev'] = {'units': 'm/m', 'data': self._milev}
            if 'ilev1' in coord_list and hasattr(self, '_ilev1'):
                coord_dict['ilev1'] = {'units': 'm/m', 'data': self._ilev1}
            elif 'ilev' in coord_list and hasattr(self, '_ilev'):
                coord_dict['ilev'] = {'units': 'm/m', 'data': self._ilev}
            coord_str = [value[3]+value[4] for key, value in
                         model_varnames.items() if value[0] == varname][0]

            # define data retrieval and wrangling logic
            # keep in mind that high-alt files have one time step per file
            # normal files have multiple time steps per file
            def func(i):  # i is the file number
                # get the data
                file = self.pattern_files[p][i]
                ncdf_filetype = RU.netcdf_type(file)
                cdf_data = RU.Dataset(file,filetype=ncdf_filetype)
                data = array(cdf_data[gvar])
                cdf_data.close()
                if data.shape[0] > 1 and file != self.pattern_files[p][-1]:
                    # if not the last file, tack on first time from next
                    next_file = self.pattern_files[p][i+1]
                    cdf_data = RU.Dataset(next_file,filetype=ncdf_filetype)
                    data_slice = array(cdf_data[gvar][0])
                    cdf_data.close()
                    data = append(data, [data_slice], axis=0)

                # wrangle the data
                if len(data.shape) == 3:
                    # (t,lat,lon) -> (t,lon,lat)
                    variable = transpose(data, (0, 2, 1))
                    out = self.wrap_3Dlatlon(varname, variable)
                    if out.shape[0] == 1:  # high-alt data has one time
                        return squeeze(out)
                    else:
                        return out
                # 4D specific logic from here on down
                # (t,h,lat,lon) -> (t,lon,lat,h)
                variable = transpose(data, (0, 3, 2, 1))
                if 'lat' in coord_list:
                    out = self.wrap_4Dlatlon(varname, variable)
                    if out.shape[0] == 1:  # high-alt data has one time
                        return squeeze(out)
                    else:
                        return out
                # otherwise, look for undefined top rows and remove them
                # they need to be removed to avoid interpolation problems
                # only occurs in mlon/mlat/milev dependent variables (H_imlev)
                top_idx = -1
                top_shape = list(variable[:, :, :, top_idx].shape)
                top_size = top_shape[0] * top_shape[1] * top_shape[2]  # 3D arr
                idx_top = where(variable[:, :, :, top_idx] > 1e+35)[0]
                while top_size == len(idx_top):  # replace undefined top row(s)
                    if verbose:
                        print('All values at max milev are 1e+36 for ' +
                              f'{varname}. Slicing off top array.')
                    variable[:, :, :, top_idx] = NaN
                    top_idx -= 1
                    idx_top = where(variable[:, :, :, top_idx] > 1e+35)[0]
                if variable.shape[0] == 1:  # high-alt data has one time
                    return squeeze(variable)
                else:
                    return variable

            # define data retrieval and wrangling logic
            # keep in mind that high-alt files have one time step per file
            # normal files have multiple time steps per file
            # func_custom returns interpolator i log space (for variable names starting with "rho" or "N_")
            def func_custom(i):  # i is the file number
                # get the data
                file = self.pattern_files[p][i]
                ncdf_filetype = RU.netcdf_type(file)
                cdf_data = RU.Dataset(file,filetype=ncdf_filetype)
                data = array(cdf_data[gvar])
                cdf_data.close()
                if data.shape[0] > 1 and file != self.pattern_files[p][-1]:
                    # if not the last file, tack on first time from next
                    next_file = self.pattern_files[p][i+1]
                    cdf_data = RU.Dataset(next_file,filetype=ncdf_filetype)
                    data_slice = array(cdf_data[gvar][0])
                    cdf_data.close()
                    data = append(data, [data_slice], axis=0)

                # wrangle the data
                if len(data.shape) == 3:
                    # (t,lat,lon) -> (t,lon,lat)
                    variable = transpose(data, (0, 2, 1))
                    out = self.wrap_3Dlatlon(varname, variable)
                    if out.shape[0] == 1:  # high-alt data has one time
                        return squeeze(out)
                    else:
                        return out
                # 4D specific logic from here on down
                # (t,h,lat,lon) -> (t,lon,lat,h)
                variable = transpose(data, (0, 3, 2, 1))

                coord_dict_data = [ coord_dict[key]['data'] for key in coord_dict ]
                time_index_start = self.times[p]['start_index'][i]
                time_index_end = self.times[p]['start_index'][i+1]+1
                times_file = self.times[p]['all'][time_index_start:time_index_end]
                coord_dict_data[0] = array(times_file)
        
                if 'lat' in coord_list:
                    out = log(self.wrap_4Dlatlon(varname, variable))
                    if out.shape[0] == 1:  # high-alt data has one time
                        rgi = rgiND(coord_dict_data[1:], squeeze(out), bounds_error=False,fill_value=NaN)
                        def interp3d_custom(xvec):
                            return exp(rgi(xvec))                    
                        return interp3d_custom
                        #return squeeze(out)
                    else:
                        rgi = rgiND(coord_dict_data, out, bounds_error=False,fill_value=NaN)
                        def interp4d_custom(xvec):
                            return exp(rgi(xvec))
                        return interp4d_custom
                        #return out
                # otherwise, look for undefined top rows and remove them
                # they need to be removed to avoid interpolation problems
                # only occurs in mlon/mlat/milev dependent variables (H_imlev)
                top_idx = -1
                top_shape = list(variable[:, :, :, top_idx].shape)
                top_size = top_shape[0] * top_shape[1] * top_shape[2]  # 3D arr
                idx_top = where(variable[:, :, :, top_idx] > 1e+35)[0]
                while top_size == len(idx_top):  # replace undefined top row(s)
                    if verbose:
                        print('All values at max milev are 1e+36 for ' +
                              f'{varname}. Slicing off top array.')
                    variable[:, :, :, top_idx] = NaN
                    top_idx -= 1
                    idx_top = where(variable[:, :, :, top_idx] > 1e+35)[0]

                log_data = log(variable)
            
                if variable.shape[0] == 1:  # high-alt data has one time
                    rgi = rgiND(coord_dict_data[1:], squeeze(log_data), bounds_error=False,fill_value=NaN)
                    def interp3d_custom(xvec):
                        return exp(rgi(xvec))                    
                    return interp3d_custom
                else:
                    rgi = rgiND(coord_dict_data, log_data, bounds_error=False,fill_value=NaN)
                    def interp4d_custom(xvec):
                        return exp(rgi(xvec))
                    return interp4d_custom
                
                
            # determine interpolation method based on file structuring
            if len(self.times[p]['all']) > len(self.times[p]['start']):
                # time chunking method
                interp_flag = 2
            elif len(self.times[p]['all']) == len(self.times[p]['start']):
                # time slicing method
                interp_flag = 1

            # need H functions to be gridded regardless of gridded_int value
            h_grid = True if varname in ['H_ilev', 'H_ilev1'] else gridded_int

            if varname[0:3] == "rho" or varname[0:2] == "N_":
                self = RU.Functionalize_Dataset(
                    self, coord_dict, varname, self.variables[varname],
                    h_grid, coord_str, interp_flag=interp_flag, func=func_custom, func_default='custom',
                    times_dict=self.times[p])
            else:
                self = RU.Functionalize_Dataset(
                    self, coord_dict, varname, self.variables[varname],
                    h_grid, coord_str, interp_flag=interp_flag, func=func,
                    times_dict=self.times[p])

            # perform H_ilev/H_ilev1 substitution if needed
            if isinstance(self.ilev_sub, str) and varname == self.ilev_sub:
                other_name = ['H_ilev', 'H_ilev1']
                other_name.remove(varname)  # first element is the other name
                print(f'{other_name[0]} missing in data and is needed to ' +
                      'convert the requested variables to depend on height.' +
                      f' Using {varname} instead.')
                del coord_dict[other_name[0]]
                # e.g. set ilev1 grid equal to ilev grid if H_ilev1 is missing
                coord_dict[varname[2:]] = {
                    'units': 'm/m', 'data': getattr(self, other_name[0][2:])}
                coord_units = {key: value['units'] for key, value in
                               coord_dict.items()}
                self.variables[other_name[0]] = self.variables[varname]
                # register the other variable by linking to this one
                self = RU.register_interpolator(self, other_name[0], varname,
                                                coord_units)

            # create pressure level -> km function once per ilev type
            if varname in ['H_ilev', 'H_ilev1'] or varname in self.total_ilev:
                if varname in ['H_ilev', 'H_ilev1']:  # create custom interp
                    new_varname = 'P'+coord_list[-1][1:]
                    kms = getattr(self, '_km_'+coord_list[-1])
                    kms_max = getattr(self, '_km_'+coord_list[-1]+'_max')
                    kms_min = getattr(self, '_km_'+coord_list[-1]+'_min')
                    # perform unit conversion if needed
                    if self.variables[varname]['units'] != 'km':
                        self[varname+'km_ijk[km]'] = varname + '_ijk'
                        km_interp = self[varname+'km_ijk']
                    else:
                        km_interp = self[varname+'_ijk']
                    # Import and call custom interpolator
                    units = 'm/m'
                    self[new_varname], interp_ijk = RU.PLevelInterp(
                        km_interp, coord_dict['time']['data'],
                        coord_dict['lon']['data'], coord_dict['lat']['data'],
                        coord_dict[coord_list[-1]]['data'], units, kms,
                        [kms_min, kms_max])
                    # kms is a 1D array of the median height values in km
                else:  # define by function composition
                    new_varname = varname.split('_ilev')[0]
                    units = self.variables[varname]['units']
                    # substitute kms array if height was also substituted
                    if isinstance(self.ilev_sub, str):
                        other_name = ['ilev', 'ilev1']
                        other_name.remove(self.ilev_sub[2:])
                        kms = getattr(self, '_km_'+other_name[0])
                        kms_max = getattr(self, '_km_'+other_name[0]+'_max')
                        kms_min = getattr(self, '_km_'+other_name[0]+'_min')
                        self[new_varname] = varname+'(P'+other_name[0][1:]+')'
                        interp_ijk = self[new_varname]
                        self[new_varname].meta['arg_units'] = \
                            self['P'+other_name[0][1:]].meta['arg_units']
                    else:
                        kms = getattr(self, '_km_'+coord_list[-1])
                        kms_max = getattr(self, '_km_'+coord_list[-1]+'_max')
                        kms_min = getattr(self, '_km_'+coord_list[-1]+'_min')
                        self[new_varname] = varname+'(P'+coord_list[-1][1:]+')'
                        interp_ijk = self[new_varname]
                        self[new_varname].meta['arg_units'] = \
                            self['P'+coord_list[-1][1:]].meta['arg_units']
                self.variables[new_varname] = {'units': units}

                # create gridded interpolator if requested
                if h_grid:
                    self = RU.register_griddedPlev(self, new_varname, units,
                                                   interp_ijk, coord_dict, kms,
                                                   [kms_min, kms_max])
            return
    return MODEL
