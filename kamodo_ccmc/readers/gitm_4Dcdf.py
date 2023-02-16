
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


def MODEL():
    from time import perf_counter
    from glob import glob
    from os.path import basename, isfile
    from numpy import array, unique, NaN
    from datetime import datetime, timezone
    from netCDF4 import Dataset
    from kamodo import Kamodo
    import kamodo_ccmc.readers.reader_utilities as RU

    class MODEL(Kamodo):
        '''GITM model data reader.

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
            - GITM output files are given in multiple binary files per time
              step. No attempt is made to combine these files into a single
              file per time step, but each is converted into a netCDF4 file
              using code adapted from a script written by Dan Welling and
              Angeline Burrell (used with permission).
            - The converted files are small and are created with one time step
              per file, so interpolation method 1 is chosen. The standard SciPy
              interpolator is used.
        '''
        def __init__(self, file_dir, variables_requested=[],
                     filetime=False, verbose=False, gridded_int=True,
                     printfiles=False, **kwargs):
            super(MODEL, self).__init__()
            self.modelname = 'GITM'

            # first, check for file list, create if DNE
            list_file = file_dir + self.modelname + '_list.txt'
            time_file = file_dir + self.modelname + '_times.txt'
            self.times, self.pattern_files = {}, {}
            if not isfile(list_file) or not isfile(time_file):
                # determine if calc of 2D variables is necessary
                total_files = sorted(glob(file_dir+'*.nc'))
                if sum(['2D' in file for file in total_files]) > 0:
                    flag_2D = True  # 2D files found, will not calc 2D vars
                else:  # try bin files...
                    bin_files = sorted(glob(file_dir+'2D*.bin'))
                    if len(bin_files) > 0:
                        flag_2D = True
                    else:
                        flag_2D = False  # calc will occur

                # check for and convert any files not converted yet
                from kamodo_ccmc.readers.gitm_tocdf import GITMbin_toCDF as\
                    toCDF
                toCDF(file_dir, flag_2D)

                t0 = perf_counter()  # begin timer
                # figure out types of files present (2DTEC, 3DALL, 3DLST, etc)
                total_files = sorted(glob(file_dir+'*.nc'))
                patterns = sorted(unique([basename(f).split('_t')[0] for f in
                                          total_files]))
                self.filename = ''.join([f+',' for f in total_files])[:-1]

                # establish time attributes
                for p in patterns:
                    # get list of files to loop through later
                    pattern_files = sorted(glob(file_dir+p+'*.nc'))
                    self.pattern_files[p] = pattern_files
                    self.times[p] = {'start': [], 'end': [], 'all': []}

                    # loop through to get times
                    for f in range(len(pattern_files)):
                        cdf_data = Dataset(pattern_files[f])
                        # hours since midnight first file
                        tmp = array(cdf_data.variables['time'])
                        self.times[p]['start'].append(tmp[0])
                        self.times[p]['end'].append(tmp[-1])
                        self.times[p]['all'].extend(tmp)
                        cdf_data.close()
                    self.times[p]['start'] = array(self.times[p]['start'])
                    self.times[p]['end'] = array(self.times[p]['end'])
                    self.times[p]['all'] = array(self.times[p]['all'])

                # establish date of first file
                cdf_data = Dataset(total_files[0], 'r')
                string_date = cdf_data.filedate
                cdf_data.close()
                self.filedate = datetime.strptime(
                    string_date+' 00:00:00', '%Y-%m-%d %H:%M:%S'
                    ).replace(tzinfo=timezone.utc)

                # create time list file if DNE
                RU.create_timelist(list_file, time_file, self.modelname,
                                   self.times, self.pattern_files,
                                   self.filedate)
            else:  # read in data and time grids from file list
                self.times, self.pattern_files, self.filedate, self.filename =\
                    RU.read_timelist(time_file, list_file)
            if filetime:
                return  # return times as is to prevent infinite recursion

            # store variables
            self.missing_value = NaN
            self.varfiles = {}  # store which variable came from which file
            self.gvarfiles = {}  # store file variable name similarly
            self.err_list = []
            self.variables = {}

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

            # loop through file patterns for coordinate grids + var mapping
            for p in self.pattern_files.keys():
                pattern_files = self.pattern_files[p]

                # get coordinates from first file
                cdf_data = Dataset(pattern_files[0], 'r')
                setattr(self, '_lat_'+p, array(cdf_data.variables['lat']))
                setattr(self, '_lon_'+p, array(cdf_data.variables['lon']))
                if '2D' not in p:
                    setattr(self, '_height_'+p,
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
                self.varfiles[p] = [model_varnames[key][0] for
                                    key in gvar_list]
                self.gvarfiles[p] = gvar_list
                cdf_data.close()

            # collect all possible variables in set of files and return
            var_list = []
            for p in self.varfiles.keys():
                var_list.extend(self.varfiles[p])
            if variables_requested == 'all':
                self.var_dict = {value[0]: value[1:] for key, value in
                                 model_varnames.items() if value[0] in
                                 var_list}
                return

            # loop through patterns to store variable keys and units
            for p in self.pattern_files.keys():
                # initialize variable dictionaries
                for gvar in self.gvarfiles[p]:
                    self.variables[model_varnames[gvar][0]] = {
                        'units': model_varnames[gvar][-1], 'data': p}

            # remove successful variables from err_list
            self.err_list = list(unique(self.err_list))
            self.err_list = [item for item in self.err_list if item not in
                             self.variables.keys()]
            if len(self.err_list) > 0:
                print('Some requested variables are not available in the ' +
                      'files found:\n',
                      self.pattern_files.keys(), self.err_list)

            if printfiles:
                print(f'{len(self.filename)} Files:')
                files = self.filename.split(',')
                for f in files:
                    print(f)

            # return if only one time found - interpolator code will break
            if ',' not in self.filename:
                print('Not enough times found in given directory.')
                return

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

        def register_variable(self, varname, gridded_int):
            """Registers an interpolator with proper signature"""

            # determine which file the variable came from, retrieve the coords
            key = self.variables[varname]['data']
            gvar = [key for key, value in model_varnames.items() if
                    value[0] == varname][0]  # variable name in file
            coord_list = [value[-2] for key, value in
                          model_varnames.items() if value[0] == varname][0]
            coord_dict = {'time': {'units': 'hr',
                                   'data': self.times[key]['all']}}
            # get the correct coordinates
            coord_dict['lon'] = {
                'data': getattr(self, '_lon_'+key), 'units': 'deg'}
            coord_dict['lat'] = {
                'data': getattr(self, '_lat_'+key), 'units': 'deg'}
            if len(coord_list) == 4:
                coord_dict['height'] = {
                    'data': getattr(self, '_height_'+key),
                    'units': 'km'}
            coord_str = [value[3]+value[4] for key, value in
                         model_varnames.items() if value[0] == varname][0]

            def func(i):
                '''i is the file number.'''
                # get data from file
                file = self.pattern_files[key][i]
                cdf_data = Dataset(file)
                data = array(cdf_data.variables[gvar])
                cdf_data.close()
                return data

            # define and register the interpolators
            self = RU.Functionalize_Dataset(self, coord_dict, varname,
                                            self.variables[varname],
                                            gridded_int, coord_str,
                                            interp_flag=1, func=func)
    return MODEL
