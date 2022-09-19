from datetime import datetime, timezone, timedelta


# file_prefix = 'C:/Users/rringuet/Kamodo_WinDev1/GITM/3DLST_t150317'
# varnames in cdf files are standardized (value[0])
model_varnames = {'r_Ar': ['mmr_Ar', 'mass mixing ratio of argon/neutrals',
                           0, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                           ''],
                  'rho_Ar': ['rho_Ar', 'mass density of argon',
                             1, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                             'kg/m**3'],
                  'r_CH4': ['mmr_CH4', 'mass mixing ratio of methane/neutrals',
                            2, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                            ''],
                  'k': ['k', 'total conduction',
                        3, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                        'W/m/K'],
                  'Q_EUV': ['Q_EUV', 'EUV heating',
                            4, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                            'K per timestep'],
                  'rho_H': ['rho_H', 'mass density of hydrogen',
                            5, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                            'kg/m**3'],
                  'rho_Hplus': ['rho_Hplus', 'mass density of hydrogen ion',
                                6, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                  'height'], 'kg/m**3'],
                  'r_H2': ['mmr_H2', 'mass mixing ratio of molecular ' +
                           'hydrogen/neutrals', 7, 'GDZ', 'sph',
                           ['time', 'lon', 'lat', 'height'], ''],
                  'r_HCN': ['mmr_HCN', 'mass mixing ratio of hydrogen ' +
                            'cyanide/neutrals', 8, 'GDZ', 'sph',
                            ['time', 'lon', 'lat', 'height'], ''],
                  'rho_He': ['rho_He', 'mass density of helium',
                             9, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                             'kg/m**3'],
                  'rho_Heplus': ['rho_Heplus', 'mass density of helium ion',
                                 10, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                    'height'], 'kg/m**3'],
                  'HeatingEfficiency': ['Q_eff', 'heating efficiency', 11,
                                        'GDZ', 'sph', ['time', 'lon', 'lat',
                                                       'height'], ''],
                  'HeatBalanceTotall': ['Q_bal', 'heat balance total', 12,
                                        'GDZ', 'sph', ['time', 'lon', 'lat',
                                                       'height'], ''],
                  'rho_N2': ['rho_N2', 'mass density of molecular nitrogen',
                             13, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                'height'], 'kg/m**3'],
                  'rho_N2plus': ['rho_N2plus', 'mass density of molecular ' +
                                 'nitrogen ion', 14, 'GDZ', 'sph',
                                 ['time', 'lon', 'lat', 'height'], 'kg/m**3'],
                  'rho_Nplus': ['rho_Nplus', 'mass density of atomic ' +
                                'nitrogen ion', 15, 'GDZ', 'sph',
                                ['time', 'lon', 'lat', 'height'], 'kg/m**3'],
                  'rho_N2D': ['rho_N2D', 'mass density of atomic nitrogen ' +
                              '(2D state)', 16, 'GDZ', 'sph',
                              ['time', 'lon', 'lat', 'height'], 'kg/m**3'],
                  'rho_N2P': ['rho_N2P', 'mass density of atomic nitrogen ' +
                              '(2P state)', 17, 'GDZ', 'sph',
                              ['time', 'lon', 'lat', 'height'], 'kg/m**3'],
                  'rho_N4S': ['rho_N4S', 'mass density of atomic nitrogen ' +
                              '(4S state)', 18, 'GDZ', 'sph',
                              ['time', 'lon', 'lat', 'height'], 'kg/m**3'],
                  'r_N2': ['mmr_N2', 'mass mixing ratio of molecular nitrogen',
                           19, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                           ''],
                  'rho_NO': ['rho_NO', 'mass density of nitric oxide',
                             20, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                'height'], 'kg/m**3'],
                  'rho_NOplus': ['rho_NOplus', 'mass density of nitric oxide' +
                                 ' ion', 21, 'GDZ', 'sph',
                                 ['time', 'lon', 'lat', 'height'], 'kg/m**3'],
                  'rho_O2': ['rho_O2', 'mass density of molecular oxygen', 22,
                             'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                             'kg/m**3'],
                  'rho_O1D': ['rho_O1D', 'mass density of atomic oxygen ' +
                              '(1D state)', 23, 'GDZ', 'sph',
                              ['time', 'lon', 'lat', 'height'], 'kg/m**3'],
                  'rho_O2plus': ['rho_O2plus', 'mass density of molecular ' +
                                 'oxygen ion', 24, 'GDZ', 'sph',
                                 ['time', 'lon', 'lat', 'height'], 'kg/m**3'],
                  'rho_O2D': ['rho_O2D', 'mass density of atomic oxygen ' +
                              '(2D state)', 25, 'GDZ', 'sph',
                              ['time', 'lon', 'lat', 'height'], '1/m**3'],
                  'rho_Oplus2P': ['rho_Oplus2P', 'mass density of atomic ' +
                                  'oxygen ion (2P state)', 26, 'GDZ', 'sph',
                                  ['time', 'lon', 'lat', 'height'], 'kg/m**3'],
                  'rho_O3P': ['rho_O3P', 'mass density of atomic oxygen ' +
                              '(3P state)', 27, 'GDZ', 'sph',
                              ['time', 'lon', 'lat', 'height'], 'kg/m**3'],
                  'rho_Oplus4SP': ['rho_Oplus4S4P', 'mass density of atomic ' +
                                   'oxygen ion (4S or 4P state)', 28, 'GDZ',
                                   'sph', ['time', 'lon', 'lat', 'height'],
                                   'kg/m**3'],
                  'L_Rad': ['Q_cool', 'radiative cooling rate', 29, 'GDZ',
                            'sph',  ['time', 'lon', 'lat', 'height'], ''],
                  'rho': ['rho_n', 'neutral mass density', 30, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'height'], 'kg/m**3'],
                  'T_n': ['T_n', 'neutral temperature', 31, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'height'], 'K'],
                  'vi_east': ['v_ieast', 'zonal ion wind velocity (east)',
                              32, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                 'height'], 'm/s'],
                  'vi_north': ['v_inorth', 'meridional ion wind velocity ' +
                               '(north)', 33, 'GDZ', 'sph',
                               ['time', 'lon', 'lat', 'height'], 'm/s'],
                  'vi_up': ['v_iup', 'vertical ion wind velocity (up)',
                            34, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                               'height'], 'm/s'],
                  'vn_east': ['v_neast', 'zonal neutral wind velocity (east)',
                              35, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                 'height'], 'm/s'],
                  'vn_north': ['v_nnorth', 'meridional neutral wind velocity' +
                               ' (north)', 36, 'GDZ', 'sph',
                               ['time', 'lon', 'lat', 'height'], 'm/s'],
                  'vn_up': ['v_nup', 'vertical neutral wind velocity (up)',
                            37, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                            'm/s'],
                  'v_N2_up': ['v_N2up', 'vertical velocity of molecular ' +
                              'nitrogen (up)', 38, 'GDZ', 'sph',
                              ['time', 'lon', 'lat', 'height'], 'm/s'],
                  'v_N4S_up': ['v_Nstate4Sup', 'vertical velocity of atomic ' +
                               'nitrogen (4S state) (up)', 39, 'GDZ', 'sph',
                               ['time', 'lon', 'lat', 'height'], 'm/s'],
                  'v_N_up': ['v_Nup', 'vertical velocity of atomic nitrogen ' +
                             '(up)', 40, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                        'height'], 'm/s'],
                  'v_O2_up': ['v_O2up', 'vertical velocity of molecular ' +
                              'oxygen (up)', 41, 'GDZ', 'sph',
                              ['time', 'lon', 'lat', 'height'], 'm/s'],
                  'v_O3P_up': ['v_Ostate3Pup', 'vertical velocity of atomic' +
                               ' (3P state) (up)', 42, 'GDZ', 'sph',
                               ['time', 'lon', 'lat', 'height'], 'm/s'],
                  'v_He_up': ['v_Heup', 'vertical velocity of helium (up)',
                              43, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                 'height'], 'm/s'],
                  'N_e': ['N_e', 'electron number density', 44, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'height'], '1/m**3'],
                  'ElectronAverageEnergy': ['E_eavg', 'average electron ' +
                                            'energy', 45, 'GDZ', 'sph',
                                            ['time', 'lon', 'lat', 'height'],
                                            'J'],
                  'T_e': ['T_e', 'electron temperature', 46, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'height'], 'K'],
                  'T_i': ['T_i', 'ion temperature', 47, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'height'], 'K'],
                  'SolarZenithAngle': ['SZA', 'solar zenith angle', 48, 'GDZ',
                                       'sph',  ['time', 'lon', 'lat'],
                                       'radians'],
                  'rho_CO2': ['rho_CO2', 'mass density of carbon dioxide', 49,
                              'GDZ', 'sph',  ['time', 'lon', 'lat', 'height'],
                              'kg/m**3'],
                  'DivJu FL': ['DivI_nfl', 'divergence of the neutral ' +
                               'wind-driven currents integrated along the ' +
                               'field-line', 50, 'GDZ', 'sph',
                               ['time', 'lon', 'lat', 'height'], ''],
                  'DivJuAlt': ['DivI_nalt', 'divergence of the neutral ' +
                               'wind-driven currents integrated along the ' +
                               'altitude', 51, 'GDZ', 'sph',
                               ['time', 'lon', 'lat', 'height'], ''],
                  'ElectronEnergyFlux': ['Phi_eE', 'electron energy flux', 52,
                                         'GDZ', 'sph', ['time', 'lon', 'lat',
                                                        'height'], 'J/m**2'],
                  'Field Line Length': ['s_Bfield', 'magnetic field arc ' +
                                        'line length', 53, 'GDZ', 'sph',
                                        ['time', 'lon', 'lat', 'height'], 'm'],
                  'sigma_P': ['sigma_P', 'Pedersen conductivity', 54, 'GDZ',
                              'sph',  ['time', 'lon', 'lat', 'height'], 'S/m'],
                  'V': ['V', 'electric potential', 57, 'GDZ', 'sph',
                        ['time', 'lon', 'lat', 'height'], 'V'],
                  'sigma_H': ['sigma_H', 'Hall conductivity', 58, 'GDZ', 'sph',
                              ['time', 'lon', 'lat', 'height'], 'S/m'],
                  'I_R2': ['j_R2', 'region 2 electric current density',
                           59, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                           'A/m**2'],
                  'I_R1': ['j_R1', 'region 1 electric current density',
                           60, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                           'A/m**2'],
                  'Ed1': ['E_perpeast', 'dynamo electric field in the' +
                          ' perpendicular to the magnetic field direction ' +
                          'that is "eastward"', 61, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'height'], ''],
                  'Ed2': ['E_perpnorth', 'dynamo electric field in the ' +
                          'perpendicular to the magnetic field direction ' +
                          'that is "northward"', 62, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'height'], ''],
                  'SolarLocalTime': ['SLT', 'solar local time', 63, 'GDZ',
                                     'sph', ['time', 'lon', 'lat'], 'hr'],
                  'E_up': ['E_up', 'vertical electric field velocity (up)',
                           64, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                           'V/m'],
                  'E_east': ['E_east', 'zonal electric field (east)',
                             65, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                'height'], 'V/m'],
                  'E_north': ['E_north', 'meridional electric field (north)',
                              66, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                 'height'], 'V/m'],
                  'E_mag': ['E_mag', 'magnitude of electric field',
                            67, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                            'V/m'],
                  'B_up': ['B_up', 'vertical magnetic field velocity (up)',
                           68, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                           'nT'],
                  'B_east': ['B_east', 'zonal magnetic field (east)',
                             69, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                'height'], 'nT'],
                  'B_north': ['B_north', 'meridional magnetic field (north)',
                              70, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                 'height'], 'nT'],
                  'B_mag': ['B_mag', 'magnitude of magnetic field', 71, 'GDZ',
                            'sph', ['time', 'lon', 'lat', 'height'], 'nT'],
                  'MagLat': ['lat_B', 'magnetic latitude', 72, 'GDZ', 'sph',
                             ['time', 'lon', 'lat', 'height'], 'deg'],
                  'MagLon': ['lon_B', 'magnetic longitude', 73, 'GDZ', 'sph',
                             ['time', 'lon', 'lat', 'height'], 'deg'],
                  'g': ['g', 'gravitational acceleration', 74, 'GDZ', 'sph',
                        ['time', 'lon', 'lat', 'height'], 'm/s**2'],
                  'GradP_east': ['GradP_east', 'zonal component of gradient ' +
                                 'of sum of ion and electron pressures (east)',
                                 75, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                    'height'], 'Pa/m'],
                  'GradP_north': ['GradP_north', 'meridional component of ' +
                                  'gradient of sum of ion and electron ' +
                                  'pressures (north)', 76, 'GDZ', 'sph',
                                  ['time', 'lon', 'lat', 'height'], 'Pa/m'],
                  'GradP_up': ['GradP_up', 'vertical component of gradient ' +
                               'of sum of ion and electron pressures (up)',
                               77, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                  'height'], 'Pa/m'],
                  'nu_in': ['nu_ion', 'ion neutral collision frequency', 78,
                            'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                            '1/s'],
                  'ChemicalHeatingRate': ['Q_chem', 'chemical heating rate',
                                          79, 'GDZ', 'sph', ['time', 'lon',
                                                             'lat', 'height'],
                                          ''],
                  'TotalAbsoluteEUV': ['Q_EUVabs', 'total absolute EUV', 80,
                                       'GDZ', 'sph',  ['time', 'lon', 'lat',
                                                       'height'],
                                       'K per timestep'],
                  'Q_Ocool': ['Q_Ocool', 'cooling rate of atomic oxygen', 81,
                              'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                              'K per timestep'],
                  'Q_Joule': ['Q_Joule', 'joule heating', 82, 'GDZ', 'sph',
                              ['time', 'lon', 'lat', 'height'],
                              'K per timestep'],
                  'Q_Auroral': ['Q_auroral', 'auroral heating', 83, 'GDZ',
                                'sph', ['time', 'lon', 'lat', 'height'],
                                'K per timestep'],
                  'Q_PhotoE': ['Q_photoe', 'heating due to the ' +
                               'photoelectric effect', 84, 'GDZ', 'sph',
                               ['time', 'lon', 'lat', 'height'],
                               'K per timestep'],
                  'k_eddy': ['k_ed', 'eddy conduction', 85, 'GDZ', 'sph',
                             ['time', 'lon', 'lat', 'height'], ''],
                  'k_adiabaticeddy': ['k_edadiab', 'adiabatic eddy conduction',
                                      86, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                         'height'], ''],
                  'Q_NOcool': ['Q_NOcool', 'cooling rate of nitric oxide', 87,
                               'GDZ', 'sph',  ['time', 'lon', 'lat', 'height'],
                               'K per timestep'],
                  'k_molecular': ['k_mol', 'molecular conduction', 88, 'GDZ',
                                  'sph', ['time', 'lon', 'lat', 'height'], ''],
                  'NmF2': ['NmF2', 'Maximum electron number density in ' +
                           'F2 layer', 89, 'GDZ', 'sph',
                           ['time', 'lon', 'lat'], ''],
                  'hmF2': ['hmF2', 'Height of maximum electron number ' +
                           'density in F2 layer', 90, 'GDZ', 'sph',
                           ['time', 'lon', 'lat'], 'km'],
                  'TEC': ['TEC', 'vertical total electron content ' +
                          '(height integrated from bottom to top boundary)',
                          91, 'GDZ', 'sph', ['time', 'lon', 'lat'],
                          '10**16/m**2'],
                  'Phi_Joule': ['Phi_Joule', 'joule heat flux', 92, 'GDZ',
                                'sph', ['time', 'lon', 'lat'], 'W/m**2'],
                  'Phi_Q': ['Phi_heat', 'heat flux', 93, 'GDZ', 'sph',
                            ['time', 'lon', 'lat'], 'W/m**2'],
                  'Phi_EUV': ['Phi_EUV', 'EUV heat flux', 94, 'GDZ', 'sph',
                              ['time', 'lon', 'lat'], 'W/m**2'],
                  'Phi_NOCooling': ['Phi_NOCooling', 'NO cooling flux', 95,
                                    'GDZ', 'sph', ['time', 'lon', 'lat'],
                                    'W/m**2']
                  }


def dts_to_hrs(datetime_string, filedate):
    '''Get hours since midnight from datetime string'''
    return (datetime.strptime(datetime_string, '%Y-%m-%d %H:%M:%S'
                              ).replace(tzinfo=timezone.utc) -
            filedate).total_seconds()/3600.


def filename_to_dts(filename, string_date):
    '''Get datetime string in format "YYYY-MM-SS HH:mm:SS" from filename'''
    mmhhss = filename.split('_t')[-1].split('_')[1].split('.bin')[0]
    return string_date+' '+mmhhss[:2]+':'+mmhhss[2:4]+':'+mmhhss[4:]


def dts_to_ts(file_dts):
    '''Get datetime timestamp in UTC from datetime string'''
    return datetime.timestamp(datetime.strptime(file_dts, '%Y-%m-%d %H:%M:%S'
                                                ).replace(tzinfo=timezone.utc))


def hrs_to_ts(hrs, filedate):
    '''Add hours to filedate and return utc timestamp.'''
    return datetime.timestamp(filedate+timedelta(hours=hrs))


def ts_to_hrs(time_val, filedate):
    '''Convert utc timestamp to hours since midnight on filedate.'''
    return (datetime.utcfromtimestamp(time_val).replace(tzinfo=timezone.utc) -
            filedate).total_seconds()/3600.


def MODEL():
    from time import perf_counter
    from glob import glob
    from os.path import basename
    from numpy import array, unique, NaN, diff, abs, append, zeros, where
    from netCDF4 import Dataset
    from kamodo import Kamodo
    from kamodo_ccmc.readers.reader_utilities import regdef_4D_interpolators
    from kamodo_ccmc.readers.reader_utilities import regdef_3D_interpolators

    class MODEL(Kamodo):
        '''GITM model data reader.

        Inputs:
            full_file_prefix: a string representing the file pattern of the
                model output data.
                Note: This reader takes a file pattern of the format
                file_dir+*YYMMDD to load data for an entire day or *YYMMDD_HH
                to load data for one hour, where file_dir is the complete file
                path to the data files, and YYMMDD is the two digit year, two
                digit month, and two digit day in the desired output file names
                (e.g. 150317 for March 17 of 2015). For the hourly option, HH
                is the two digit hour assuming a 24 hour convention.
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
        def __init__(self, full_file_prefix, variables_requested=[],
                     filetime=False, verbose=False, gridded_int=True,
                     printfiles=False, fulltime=True, **kwargs):
            super(MODEL, self).__init__()
            self.modelname = 'GITM'

            # check for prepared .nc files
            file_prefix = basename(full_file_prefix)
            file_dir = full_file_prefix.split(file_prefix)[0]
            total_files = sorted(glob(file_dir+file_prefix+'*'))
            if len(file_prefix.split('_')[-1]) == 2:  # asking for hourly files
                self.patterns = unique([basename(f)[:16] for f in total_files
                                        if basename(f)[13:16] != '.nc' and
                                        'GITM' not in basename(f)])
                day_tag = False
            else:  # asking for daily files (satellite flythrough)
                self.patterns = unique([basename(f)[:13] for f in total_files
                                        if 'GITM' not in basename(f)])
                day_tag = True
            n_ncfiles = len(glob(file_dir+file_prefix+'.nc'))

            # if .nc files not there, create then
            self.conversion_test = True  # default value
            # create .nc files for each pattern
            if len(self.patterns) != n_ncfiles:
                # determine of calc of 2D variables is necessary
                if sum(['2D' in p for p in self.patterns]) > 0:
                    flag_2D = True  # 2D files found, will not calc 2D vars
                else:
                    flag_2D = False  # calc will occur

                # convert files
                from kamodo_ccmc.readers.gitm_tocdf import GITMbin_toCDF as\
                    toCDF
                test = [toCDF(file_dir+p, flag_2D) for p in self.patterns]
                if sum(test) != len(self.patterns):
                    self.conversion_test = False
                    return   # if file conversion fails, return

            t0 = perf_counter()
            # establish time attributes first for file searching, preferring 3D
            for name in self.patterns:
                if '3DLST' in name or '3DALL' in name:
                    cdf_file = name  # need to get reference height correct
            cdf_data = Dataset(file_dir+cdf_file+'.nc', 'r')
            string_date = cdf_data.filedate
            self.filedate = datetime.strptime(string_date+' 00:00:00',
                                              '%Y-%m-%d %H:%M:%S'
                                              ).replace(tzinfo=timezone.utc)

            # establish beginning and end times of files found
            files = cdf_data.file.split(',')
            # strings in format = YYYY-MM-DD HH:MM:SS
            self.datetimes = [filename_to_dts(filename, string_date) for
                              filename in [files[0], files[-1]]]
            # timestamps in UTC
            self.filetimes = [dts_to_ts(file_dts) for file_dts in
                              self.datetimes]
            t = array(cdf_data.variables['time'])  # accurate to the sec
            if len(t) > 1:
                self.dt = diff(t).max()*3600.  # t is in hours since midnight
            else:
                self.dt = 0

            if filetime and not fulltime:
                cdf_data.close()
                return  # return times as is to prevent infinite recursion

            # if variables are given as integers, convert to standard names
            if len(variables_requested) > 0:
                if isinstance(variables_requested[0], int):
                    tmp_var = [value[0] for key, value in
                               model_varnames.items() if value[2] in
                               variables_requested]
                    variables_requested = tmp_var

            if fulltime:  # add boundary time (default value)
                t_search = perf_counter()
                # find other files with same pattern
                files = sorted(glob(file_dir+'*'))
                if day_tag:  # perform task with day files
                    file_prefixes = unique(['*'+basename(f)[7:13] for f in
                                            files if 'GITM' not in basename(f)
                                            and '.nc' not in basename(f)])
                else:   # perform task with hourly files
                    file_prefixes = unique(['*'+basename(f)[7:16] for f in
                                            files if '.nc' not in basename(f)
                                            and 'GITM' not in basename(f)])

                # find closest file by utc timestamp
                # gitm has an open time at the end
                # so need a beginning time from the closest file
                # files are automatically sorted by YYMMDD
                # so next file is next in the list
                current_idx = where(file_prefixes == file_prefix)[0]
                if current_idx+1 == len(file_prefixes):
                    if verbose:
                        print('No later file available.')
                    filecheck = False
                    if filetime:
                        cdf_data.close()
                        return
                else:
                    # +1 for adding an end time
                    min_file_prefix = file_prefixes[current_idx+1][0]
                    kamodo_test = MODEL(file_dir+min_file_prefix,
                                        filetime=True, fulltime=False)
                    if not kamodo_test.conversion_test:
                        if verbose:
                            print('No later file available.')
                        filecheck = False
                        if filetime:
                            cdf_data.close()
                            return
                    else:
                        time_test = abs(kamodo_test.filetimes[0] -
                                        self.filetimes[1])
                        if verbose:
                            print(f'Took {perf_counter()-t_search:.5f}s to ' +
                                  'find the best file.')
                        # if nearest file time at least within one timestep
                        if time_test <= self.dt:
                            filecheck = True
                            self.datetimes[1] = kamodo_test.datetimes[0]
                            self.filetimes[1] = kamodo_test.filetimes[0]

                            # time only version if returning time for searching
                            if filetime:
                                cdf_data.close()
                                return  # return object with additional time

                            # get kamodo object with same requested variables
                            if verbose:
                                print(f'Took {perf_counter()-t0:.3f}s to ' +
                                      'find closest file.')
                            kamodo_neighbor = MODEL(
                                file_dir+min_file_prefix,
                                variables_requested=variables_requested,
                                fulltime=False)
                            short_data = kamodo_neighbor.short_data
                            if verbose:
                                print(f'Took {perf_counter()-t0:.3f}s to ' +
                                      'get data from closest file.')
                        else:
                            if verbose:
                                print('No later file found within ' +
                                      f'{self.dt:.1f}s.')
                            filecheck = False
                            if filetime:
                                cdf_data.close()
                                return

            # store coordinates for reference, time is NOT the same for all
            self._lat = array(cdf_data.variables['lat'])  # -90 to 90
            self._lon = array(cdf_data.variables['lon'])  # -180 to 180
            self._height = array(cdf_data.variables['height'])
            cdf_data.close()

            # store variables
            self.filename = []
            self.missing_value = NaN
            self._registered = 0
            self.varfiles = {}  # store which variable came from which file
            self.gvarfiles = {}  # store file variable name similarly
            self.err_list = []
            self.variables = {}

            # perform initial check on variables_requested list
            if len(variables_requested) > 0 and fulltime and\
                    variables_requested != 'all':
                test_list = [value[0] for key, value in model_varnames.items()]
                err_list = [item for item in variables_requested if item not in
                            test_list]
                if len(err_list) > 0:
                    print('Variable name(s) not recognized:', err_list)

            # loop through files to collect variables and mapping
            for i in range(len(self.patterns)):
                cdf_data = Dataset(file_dir+self.patterns[i]+'.nc', 'r')
                files = cdf_data.file.split(',')
                self.filename.extend(files)

                # store coordinates, time is always the same between files
                setattr(self, '_time'+str(i),
                        array(cdf_data.variables['time']))
                setattr(self, '_lat'+str(i), array(cdf_data.variables['lat']))
                setattr(self, '_lon'+str(i), array(cdf_data.variables['lon']))
                if '2D' not in self.patterns[i]:
                    setattr(self, '_height'+str(i),
                            array(cdf_data.variables['height']))

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
                        self.err_list.extend(err_list)  # add to master list
                else:
                    gvar_list = [key for key in model_varnames.keys()
                                 if key in cdf_data.variables.keys()]
                # store which file these variables came from
                self.varfiles[self.patterns[i]] = [model_varnames[key][0] for
                                                   key in gvar_list]
                self.gvarfiles[self.patterns[i]] = gvar_list
                cdf_data.close()

            # collect all possible variables in set of files and return
            if not fulltime and variables_requested == 'all':
                self.var_dict, gvar_list = {}, []
                # loop through gvar_list stored for files to make a master list
                for i in range(len(self.patterns)):
                    gvar_list += self.gvarfiles[self.patterns[i]]
                self.var_dict = {value[0]: value[1:] for key, value in
                                 model_varnames.items() if key in gvar_list}
                return

            # loop through files to store variable data and units
            for i in range(len(self.patterns)):
                gvar_list = self.gvarfiles[self.patterns[i]]
                cdf_data = Dataset(file_dir+self.patterns[i]+'.nc', 'r')
                variables = {model_varnames[key][0]: {
                    'units': cdf_data.variables[key].units,
                    'data': array(cdf_data.variables[key])}
                             for key in gvar_list}
                for key in variables.keys():
                    self.variables[key] = variables[key]
                cdf_data.close()

            # prepare and return data only for last timestamp
            if not fulltime:
                self.variables['time'] = self.filetimes[0]
                for i in range(len(self.patterns)):
                    self.variables['time'+str(i)] = hrs_to_ts(
                        getattr(self, '_time'+str(i))[0], self.filedate)
                self.short_data = self.variables
                return

            # remove successful variables from err_list
            self.err_list = list(unique(self.err_list))
            self.err_list = [item for item in self.err_list if item not in
                             self.variables.keys()]
            if len(self.err_list) > 0:
                print('Some requested variables are not available in the ' +
                      'files found:\n',
                      self.patterns, self.err_list)

            if printfiles:
                print(f'{len(self.filename)} Files:')
                for file in self.filename:
                    print(file)

            # Add new time and store time coordinate data as class attribute
            if filecheck:
                # new time in hours since midnight
                new_time = ts_to_hrs(short_data['time'], self.filedate)
                self._time = append(t, new_time)
                for i in range(len(self.patterns)):
                    new_time = ts_to_hrs(short_data['time'], self.filedate)
                    setattr(self, '_time'+str(i),
                            append(getattr(self, '_time'+str(i)), new_time))
            else:
                self._time = t

            # return if only one file found - interpolator code will break
            if len(self._time) < 2:
                print('Not enough times found with given file prefix.')
                return

            # register interpolators for each variable
            t_reg = perf_counter()
            # store original list b/c gridded interpolators change keys list
            varname_list = [key for key in self.variables.keys()]
            for varname in varname_list:
                if len(self.variables[varname]['data'].shape) == 3:
                    if filecheck:  # if neighbor found
                        # append data for last time stamp, transpose
                        data_shape = list(
                            self.variables[varname]['data'].shape)
                        data_shape[0] += 1  # add space for time
                        new_data = zeros(data_shape)
                        # put in current data
                        new_data[:-1, :, :] = self.variables[varname]['data']
                        # add in data for additional time
                        new_data[-1, :, :] =\
                            short_data[varname]['data'][0, :, :]
                    else:
                        new_data = self.variables[varname]['data']
                    self.variables[varname]['data'] = new_data
                    self.register_3D_variable(self.variables[varname]['units'],
                                              self.variables[varname]['data'],
                                              varname, gridded_int)
                elif len(self.variables[varname]['data'].shape) == 4:
                    if filecheck:
                        # append data for last time stamp, transpose
                        data_shape = list(
                            self.variables[varname]['data'].shape)
                        data_shape[0] += 1  # add space for time
                        new_data = zeros(data_shape)
                        # put in current data
                        new_data[:-1, :, :, :] =\
                            self.variables[varname]['data']
                        # add in data for additional time
                        new_data[-1, :, :, :] =\
                            short_data[varname]['data'][0, :, :, :]
                    else:
                        new_data = self.variables[varname]['data']
                    self.variables[varname]['data'] = new_data
                    self.register_4D_variable(self.variables[varname]['units'],
                                              self.variables[varname]['data'],
                                              varname, gridded_int)
            if verbose:
                print(f'Took {perf_counter()-t_reg:.5f}s to register ' +
                      f'{len(varname_list)} variables.')
            if verbose:
                print(f'Took a total of {perf_counter()-t0:.5f}s to kamodofy' +
                      f' {len(varname_list)} variables.')

        # define and register a 3D variable------------------------------
        def register_3D_variable(self, units, variable, varname, gridded_int):
            """Registers a 3d interpolator with 3d signature"""

            # determine which file the variable came from
            for i in range(len(self.patterns)-1, -1, -1):
                if varname in self.varfiles[self.patterns[i]]:
                    # get the correct coordinates
                    time = getattr(self, '_time'+str(i))
                    lat = getattr(self, '_lat'+str(i))
                    lon = getattr(self, '_lon'+str(i))
                    break

            # define and register the interpolators
            xvec_dependencies = {'time': 'hr', 'lon': 'deg', 'lat': 'deg'}
            self = regdef_3D_interpolators(self, units, variable, time, lon,
                                           lat, varname, xvec_dependencies,
                                           gridded_int)
            return

        # define and register a 4D variable -----------------------------------
        def register_4D_variable(self, units, variable, varname, gridded_int):
            """Registers a 4d interpolator with 4d signature"""

            # determine which file the variable came from
            for i in range(len(self.patterns)-1, -1, -1):
                if varname in self.varfiles[self.patterns[i]]:
                    # get the correct coordinates
                    time = getattr(self, '_time'+str(i))
                    lat = getattr(self, '_lat'+str(i))
                    lon = getattr(self, '_lon'+str(i))
                    height = getattr(self, '_height'+str(i))
                    break

            # define and register the interpolators
            xvec_dependencies = {'time': 'hr', 'lon': 'deg', 'lat': 'deg',
                                 'height': 'km'}
            self = regdef_4D_interpolators(self, units, variable, time, lon,
                                           lat, height, varname,
                                           xvec_dependencies, gridded_int)
            return
    return MODEL
