# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 15:57:01 2022

@author: rringuet
files given have one hour of data, open at end (days since 1979-01-01 00:00:00)
lon and mlon need wrapping
lon is 0-360 (shape=144), mlon is -180 to 180 (shape=80)
lat (shape=96) and mlat (shape=97) include both poles
two pressure level variables:
    lev.shape = 130 at midpoints (hPa)
    ilev.shape = 131 at interfaces (hPa)
data structure order = (time, lev/ilev, lat/mlat, lon/mlon)
"""
# starts with variables from the h1 files. No h2 files in this set
model_varnames = {'co2vmr': ['mmr_CO2', 'co2 volume mixing ratio', 0,
                             'GDZ', 'sph', ['time'], ''],
                  'ch4vmr': ['mmr_CH4', 'ch4 volume mixing ratio', 0,
                             'GDZ', 'sph', ['time'], ''],
                  'n2ovmr': ['mmr_N2O', 'n2o volume mixing ratio', 0,
                             'GDZ', 'sph', ['time'], ''],
                  'f11vmr': ['mmr_F11', 'f11 volume mixing ratio', 0,
                             'GDZ', 'sph', ['time'], ''],
                  'f12vmr': ['mmr_F12', 'f12 volume mixing ratio', 0,
                             'GDZ', 'sph', ['time'], ''],
                  'sol_tsi': ['TSI', 'total solar irradiance', 0, 'GDZ', 'sph',
                              ['time'], 'W/m**2'],
                  'colat_crit1': ['Colat1', 'First co-latitude of electro-' +
                                  'potential critical angle', 0, 'GDZ', 'sph',
                                  ['time'],  'degrees'],
                  'colat_crit2': ['Colat2', 'Second co-latitude of electro-' +
                                  'potential critical angle', 0, 'GDZ', 'sph',
                                  ['time'], 'degrees'],
                  'ED1': ['E_east', 'ED1: Eastward Electric Field', 0,
                          'MAG', 'sph', ['time', 'mlon', 'mlat'], 'V/m'],
                  'ED2': ['E_equator', 'ED2: Equatorward Electric Field', 0,
                          'MAG', 'sph', ['time', 'mlon', 'mlat'], 'V/m'],
                  'EDYN_ZIGM11_PED': ['Sigma_P', 'Pedersen Conductance', 0,
                                      'MAG', 'sph', ['time', 'mlon', 'mlat'],
                                      'S'],
                  'EDYN_ZIGM2_HAL': ['Sigma_H', 'Hall Conductance', 0, 'MAG',
                                     'sph', ['time', 'mlon', 'mlat'], 'S'],
                  'ElecColDens': ['TEC', 'Electron Column Density', 0, 'GDZ',
                                  'sph', ['time', 'lon', 'lat'],
                                  '10**16/m**2'],
                  'H': ['Hyd_ilev', 'H concentration', 0, 'GDZ', 'sph',
                        ['time', 'lon', 'lat', 'ilev'], 'mol/mol'],
                  'H_2': ['Hyd', 'H concentration', 0, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'height'], 'mol/mol'],
                  'O': ['Oxy_ilev', 'O concentration', 0, 'GDZ', 'sph',
                        ['time', 'lon', 'lat', 'ilev'], 'mol/mol'],
                  'O_2': ['Oxy', 'O concentration', 0, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'height'], 'mol/mol'],
                  'O2': ['O2_ilev', 'O2 concentration', 0, 'GDZ', 'sph',
                         ['time', 'lon', 'lat', 'ilev'], 'mol/mol'],
                  'O2_2': ['O2', 'O2 concentration', 0, 'GDZ', 'sph',
                           ['time', 'lon', 'lat', 'height'], 'mol/mol'],
                  'OMEGA': ['v_up_ilev', 'Vertical velocity (pressure)', 0,
                            'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                            'Pa/s'],
                  'OMEGA_2': ['v_up', 'Vertical velocity (pressure)', 0, 'GDZ',
                              'sph', ['time', 'lon', 'lat', 'height'], 'Pa/s'],
                  'PHIM2D': ['phi', 'PHIM2D: Electric Potential', 0, 'MAG',
                             'sph', ['time', 'mlon', 'mlat'], 'V'],
                  'PS': ['P_surface', 'Surface pressure', 0, 'GDZ', 'sph',
                         ['time', 'lon', 'lat'], 'Pa'],
                  'T': ['T_ilev', 'Temperature', 0, 'GDZ', 'sph',
                        ['time', 'lon', 'lat', 'ilev'], 'K'],
                  'T_2': ['T', 'Temperature', 0, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'height'], 'K'],
                  'TElec': ['T_e_ilev', 'Electron Temperature', 0, 'GDZ',
                            'sph', ['time', 'lon', 'lat', 'ilev'], 'K'],
                  'TElec_2': ['T_e', 'Electron Temperature', 0, 'GDZ', 'sph',
                              ['time', 'lon', 'lat', 'height'], 'K'],
                  'TIon': ['T_i_ilev', 'Ion Temperature', 0, 'GDZ', 'sph',
                           ['time', 'lon', 'lat', 'ilev'], 'K'],
                  'TIon_2': ['T_i', 'Ion Temperature', 0, 'GDZ', 'sph',
                             ['time', 'lon', 'lat', 'height'], 'K'],
                  'U': ['v_east_ilev', 'Zonal wind', 0, 'GDZ', 'sph',
                        ['time', 'lon', 'lat', 'ilev'], 'm/s'],
                  'U_2': ['v_east', 'Zonal wind', 0, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'height'], 'm/s'],
                  'UI': ['vi_east_ilev', 'UI Zonal ion drift from edynamo', 0,
                         'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'], 'm/s'],
                  'UI_2': ['vi_east', 'UI Zonal ion drift from edynamo', 0,
                           'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                           'm/s'],
                  'V': ['v_north_ilev', 'Meridional wind', 0, 'GDZ', 'sph',
                        ['time', 'lon', 'lat', 'ilev'], 'm/s'],
                  'V_2': ['v_north', 'Meridional wind', 0, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'height'], 'm/s'],
                  'VI': ['vi_north_ilev', 'VI Meridional ion drift from ' +
                         'edynamo', 0, 'GDZ', 'sph',
                         ['time', 'lon', 'lat', 'ilev'], 'm/s'],
                  'VI_2': ['vi_north', 'VI Meridional ion drift from edynamo',
                           0, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                           'm/s'],
                  'WI': ['vi_up_ilev', 'WI Vertical ion drift from edynamo',
                         0, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                         'm/s'],
                  'WI_2': ['vi_up', 'WI Vertical ion drift from edynamo',
                           0, 'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                           'm/s'],
                  'Z3': ['H_geopot_ilev', 'Geopotential Height (above sea ' +
                         'level)', 0, 'GDZ', 'sph',
                         ['time', 'lon', 'lat', 'ilev'], 'm'],
                  'e': ['e_ilev', 'e concentration', 0, 'GDZ', 'sph',
                        ['time', 'lon', 'lat', 'ilev'], 'mol/mol'],
                  'e_2': ['e', 'e concentration', 0, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'height'], 'mol/mol'],
                  'RHO_CLUBB': ['rho_ilev1', 'Air density', 0, 'GDZ', 'sph',
                                ['time', 'lon', 'lat', 'ilev1'], 'kg/m**3'],
                  'RHO_CLUBB_lev': ['rho_ilev', 'Air density', 0, 'GDZ', 'sph',
                                    ['time', 'lon', 'lat', 'ilev'], 'kg/m**3'],
                  'RHO_CLUBB_lev_2': ['rho', 'Air density', 0, 'GDZ', 'sph',
                                      ['time', 'lon', 'lat', 'height'],
                                      'kg/m**3'],

                  # variables in h1 files also in h3 files
                  'Z3GM': ['H_geomet_ilev', 'Geometric height', 0, 'GDZ',
                           'sph', ['time', 'lon', 'lat', 'ilev'], 'm'],

                  # variables from the h2 files (also in h3 files!)
                  'EDens': ['N_e_ilev', 'e Number Density (sum of O2+,NO+' +
                            ',N2+,O+)', 0, 'GDZ', 'sph',
                            ['time', 'lon', 'lat', 'ilev'], 'cm**3'],
                  'EDens_2': ['N_e', 'e Number Density (sum of O2+,NO+,N2+' +
                              ',O+)', 0, 'GDZ', 'sph',
                              ['time', 'lon', 'lat', 'height'], 'cm**3'],
                  'HMF2': ['HmF2', 'Height of the F2 Layer', 0, 'GDZ', 'sph',
                           ['time', 'lon', 'lat'], 'km'],
                  'NMF2': ['NmF2', 'Peak Density of the F2 Layer', 0, 'GDZ',
                           'sph', ['time', 'lon', 'lat'], '1/cm**3'],
                  'OPLUS': ['N_Oplus_ilev', 'O+ number density', 0, 'GDZ',
                            'sph', ['time', 'lon', 'lat', 'ilev'], '1/cm**3'],
                  'OPLUS_2': ['N_Oplus', 'O+ number density', 0, 'GDZ', 'sph',
                              ['time', 'lon', 'lat', 'height'], '1/cm**3'],

                  # variables from the h3 files
                  'CO2': ['CO2_ilev', 'CO2 concentration', 0, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'ilev'], 'mol/mol'],
                  'CO2_2': ['CO2', 'CO2 concentration', 0, 'GDZ', 'sph',
                            ['time', 'lon', 'lat', 'height'], 'mol/mol'],
                  # 'EKGW': ['EKGW_ilev1', 'Effective Kzz due to diffusion by '
                  # gravity waves', 0, 'GDZ', 'sph',
                  #           ['time', 'lon', 'lat', 'ilev1'], 'm**2/s'],
                  # 'EKGW_2': ['EKGW', 'Effective Kzz due to diffusion by grav'
                  # ity waves', 0, 'GDZ', 'sph',
                  #          ['time', 'lon', 'lat', 'height'], 'm**2/s'],
                  'N': ['Nit_ilev', 'N concentration', 0, 'GDZ', 'sph',
                        ['time', 'lon', 'lat', 'ilev'], 'mol/mol'],
                  'N_2': ['Nit', 'N concentration', 0, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'height'], 'mol/mol'],
                  'NO': ['NO_ilev', 'NO concentration', 0, 'GDZ', 'sph',
                         ['time', 'lon', 'lat', 'ilev'], 'mol/mol'],
                  'NO_2': ['NO', 'NO concentration', 0, 'GDZ', 'sph',
                           ['time', 'lon', 'lat', 'height'], 'mol/mol'],
                  'OpDens': ['N_Oplus_ilev', 'O+ Number Density', 0, 'GDZ',
                             'sph', ['time', 'lon', 'lat', 'ilev'], '1/cm**3'],
                  'OpDens_2': ['N_Oplus', 'O+ Number Density', 0, 'GDZ', 'sph',
                               ['time', 'lon', 'lat', 'height'], '1/cm**3'],
                  'QCO2': ['Q_CO2_ilev', 'CO2 cooling', 0, 'GDZ', 'sph',
                           ['time', 'lon', 'lat', 'ilev'], 'K/s'],
                  'QCO2_2': ['Q_CO2', 'CO2 cooling', 0, 'GDZ', 'sph',
                             ['time', 'lon', 'lat', 'height'], 'K/s'],
                  'QHC2S': ['Q_HC2S_ilev', 'Cooling to Space', 0, 'GDZ', 'sph',
                            ['time', 'lon', 'lat', 'ilev'], 'K/s'],
                  'QHC2S_2': ['Q_HC2S', 'Cooling to Space', 0, 'GDZ', 'sph',
                              ['time', 'lon', 'lat', 'height'], 'K/s'],
                  'QJOULE': ['Q_Joule_ilev', 'Joule Heat', 0, 'GDZ', 'sph',
                             ['time', 'lon', 'lat', 'ilev'], 'K/s'],
                  'QJOULE_2': ['Q_Joule', 'Joule Heat', 0, 'GDZ', 'sph',
                               ['time', 'lon', 'lat', 'height'], 'K/s'],
                  'QNO': ['Q_NO_ilev', 'NO cooling', 0, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'ilev'], 'K/s'],
                  'QNO_2': ['Q_NO', 'NO cooling', 0, 'GDZ', 'sph',
                            ['time', 'lon', 'lat', 'height'], 'K/s'],
                  'QO3': ['Q_O3_ilev', 'O3 cooling', 0, 'GDZ', 'sph',
                          ['time', 'lon', 'lat', 'ilev'], 'K/s'],
                  'QO3_2': ['Q_O3', 'O3 cooling', 0, 'GDZ', 'sph',
                            ['time', 'lon', 'lat', 'height'], 'K/s'],
                  'QO3P': ['Q_O3P_ilev', 'O3P cooling', 0, 'GDZ', 'sph',
                           ['time', 'lon', 'lat', 'ilev'], 'K/s'],
                  'QO3P_2': ['Q_O3P', 'O3P cooling', 0, 'GDZ', 'sph',
                             ['time', 'lon', 'lat', 'height'], 'K/s'],
                  'QRS_TOT': ['Q_Total_ilev', 'Merged SW heating: QRS+QCP+' +
                              'QRS_EUV+QRS_CO2NIR+QRS_AUR+QTHERMAL', 0, 'GDZ',
                              'sph', ['time', 'lon', 'lat', 'ilev'], 'K/s'],
                  'QRS_TOT_2': ['Q_Total', 'Merged SW heating: QRS+QCP+' +
                                'QRS_EUV+QRS_CO2NIR+QRS_AUR+QTHERMAL', 0,
                                'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                                'K/s'],
                  'SolIonRate_Tot': ['SolIonRate_Tot_ilev', 'reaction rate ' +
                                     'group', 0, 'GDZ', 'sph',
                                     ['time', 'lon', 'lat', 'ilev'],
                                     '1/cm**3/s'],  # molecules/cm**3/s
                  'SolIonRate_Tot_2': ['SolIonRate_Tot', 'reaction rate group',
                                       0, 'GDZ', 'sph',
                                       ['time', 'lon', 'lat', 'height'],
                                       '1/cm**3/s'],  # molecules/cm**3/s
                  'TTGW': ['TTGW_ilev', 'T tendency - gravity wave drag', 0,
                           'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                           'K/s'],
                  'TTGW_2': ['TTGW', 'T tendency - gravity wave drag', 0,
                             'GDZ', 'sph', ['time', 'lon', 'lat', 'height'],
                             'K/s'],
                  'UTGW_TOTAL': ['UTGW_TOTAL_ilev', 'Total U tendency due to' +
                                 ' gravity wave drag', 0, 'GDZ', 'sph',
                                 ['time', 'lon', 'lat', 'ilev'], 'm/s**2'],
                  'UTGW_TOTAL_2': ['UTGW_TOTAL', 'Total U tendency due to ' +
                                   'gravity wave drag', 0, 'GDZ', 'sph',
                                   ['time', 'lon', 'lat', 'height'], 'm/s**2'],

                  # variables from the h4 files
                  'OMEGA_08_COS': ['OMEGA_08_COS_ilev', 'vertical pressure ' +
                                   'velocity  8hr. cos coeff.', 0, 'GDZ',
                                   'sph', ['time', 'lon', 'lat', 'ilev'],
                                   'Pa/s'],
                  'OMEGA_08_COS_2': ['OMEGA_08_COS', 'vertical pressure ' +
                                     'velocity  8hr. cos coeff.', 0, 'GDZ',
                                     'sph', ['time', 'lon', 'lat', 'height'],
                                     'Pa/s'],
                  'OMEGA_08_SIN': ['OMEGA_08_SIN_ilev', 'vertical pressure ' +
                                   'velocity  8hr. sin coeff.', 0, 'GDZ',
                                   'sph', ['time', 'lon', 'lat', 'ilev'],
                                   'Pa/s'],
                  'OMEGA_08_SIN_2': ['OMEGA_08_SIN', 'vertical pressure ' +
                                     'velocity  8hr. sin coeff.', 0, 'GDZ',
                                     'sph', ['time', 'lon', 'lat', 'height'],
                                     'Pa/s'],
                  'OMEGA_12_COS': ['OMEGA_12_COS_ilev', 'vertical pressure ' +
                                   'velocity 12hr. cos coeff.', 0, 'GDZ',
                                   'sph', ['time', 'lon', 'lat', 'ilev'],
                                   'Pa/s'],
                  'OMEGA_12_COS_2': ['OMEGA_12_COS', 'vertical pressure ' +
                                     'velocity 12hr. cos coeff.', 0, 'GDZ',
                                     'sph', ['time', 'lon', 'lat', 'height'],
                                     'Pa/s'],
                  'OMEGA_12_SIN': ['OMEGA_12_SIN_ilev', 'vertical pressure ' +
                                   'velocity 12hr. sin coeff.', 0, 'GDZ',
                                   'sph', ['time', 'lon', 'lat', 'ilev'],
                                   'Pa/s'],
                  'OMEGA_12_SIN_2': ['OMEGA_12_SIN', 'vertical pressure ' +
                                     'velocity 12hr. sin coeff.', 0, 'GDZ',
                                     'sph', ['time', 'lon', 'lat', 'height'],
                                     'Pa/s'],
                  'OMEGA_24_COS': ['OMEGA_24_COS_ilev', 'vertical pressure ' +
                                   'velocity 24hr. cos coeff.', 0, 'GDZ',
                                   'sph', ['time', 'lon', 'lat', 'ilev'],
                                   'Pa/s'],
                  'OMEGA_24_COS_2': ['OMEGA_24_COS', 'vertical pressure ' +
                                     'velocity 24hr. cos coeff.', 0, 'GDZ',
                                     'sph', ['time', 'lon', 'lat', 'height'],
                                     'Pa/s'],
                  'OMEGA_24_SIN': ['OMEGA_24_SIN_ilev', 'vertical pressure ' +
                                   'velocity 24hr. sin coeff.', 0, 'GDZ',
                                   'sph', ['time', 'lon', 'lat', 'ilev'],
                                   'Pa/s'],
                  'OMEGA_24_SIN_2': ['OMEGA_24_SIN', 'vertical pressure ' +
                                     'velocity 24hr. sin coeff.', 0, 'GDZ',
                                     'sph', ['time', 'lon', 'lat', 'height'],
                                     'Pa/s'],
                  'T_08_COS': ['T_08_COS_ilev', 'Temperature  8hr. cos coeff.',
                               0, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                               'K'],
                  'T_08_COS_2': ['T_08_COS', 'Temperature  8hr. cos coeff.', 0,
                                 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                'height'], 'K'],
                  'T_08_SIN': ['T_08_SIN_ilev', 'Temperature  8hr. sin coeff.',
                               0, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                               'K'],
                  'T_08_SIN_2': ['T_08_SIN', 'Temperature  8hr. sin coeff.',
                                 0, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                   'height'], 'K'],
                  'T_12_COS': ['T_12_COS_ilev', 'Temperature 12hr. cos coeff.',
                               0, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                               'K'],
                  'T_12_COS_2': ['T_12_COS', 'Temperature 12hr. cos coeff.',
                                 0, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                   'height'], 'K'],
                  'T_12_SIN': ['T_12_SIN_ilev', 'Temperature 12hr. sin coeff.',
                               0, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                               'K'],
                  'T_12_SIN_2': ['T_12_SIN', 'Temperature 12hr. sin coeff.',
                                 0, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                   'height'], 'K'],
                  'T_24_COS': ['T_24_COS_ilev', 'Temperature 24hr. cos coeff.',
                               0, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                               'K'],
                  'T_24_COS_2': ['T_24_COS', 'Temperature 24hr. cos coeff.',
                                 0, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                   'height'], 'K'],
                  'T_24_SIN': ['T_24_SIN_ilev', 'Temperature 24hr. sin coeff.',
                               0, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                               'K'],
                  'T_24_SIN_2': ['T_24_SIN', 'Temperature 24hr. sin coeff.',
                                 0, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                   'height'], 'K'],
                  'U_08_COS': ['U_08_COS_ilev', 'Zonal wind  8hr. cos coeff.',
                               0, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                               'm/s'],
                  'U_08_COS_2': ['U_08_COS', 'Zonal wind  8hr. cos coeff.',
                                 0, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                   'height'], 'm/s'],
                  'U_08_SIN': ['U_08_SIN_ilev', 'Zonal wind  8hr. sin coeff.',
                               0, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                               'm/s'],
                  'U_08_SIN_2': ['U_08_SIN', 'Zonal wind  8hr. sin coeff.',
                                 0, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                   'height'], 'm/s'],
                  'U_12_COS': ['U_12_COS_ilev', 'Zonal wind 12hr. cos coeff.',
                               0, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                               'm/s'],
                  'U_12_COS_2': ['U_12_COS', 'Zonal wind 12hr. cos coeff.',
                                 0, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                   'height'], 'm/s'],
                  'U_12_SIN': ['U_12_SIN_ilev', 'Zonal wind 12hr. sin coeff.',
                               0, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                               'm/s'],
                  'U_12_SIN_2': ['U_12_SIN', 'Zonal wind 12hr. sin coeff.',
                                 0, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                   'height'], 'm/s'],
                  'U_24_COS': ['U_24_COS_ilev', 'Zonal wind 24hr. cos coeff.',
                               0, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                               'm/s'],
                  'U_24_COS_2': ['U_24_COS', 'Zonal wind 24hr. cos coeff.',
                                 0, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                   'height'], 'm/s'],
                  'U_24_SIN': ['U_24_SIN_ilev', 'Zonal wind 24hr. sin coeff.',
                               0, 'GDZ', 'sph', ['time', 'lon', 'lat', 'ilev'],
                               'm/s'],
                  'U_24_SIN_2': ['U_24_SIN', 'Zonal wind 24hr. sin coeff.',
                                 0, 'GDZ', 'sph', ['time', 'lon', 'lat',
                                                   'height'], 'm/s'],
                  'V_08_COS': ['V_08_COS_ilev', 'Meridional wind  8hr. cos ' +
                               'coeff.', 0, 'GDZ', 'sph',
                               ['time', 'lon', 'lat', 'ilev'], 'm/s'],
                  'V_08_COS_2': ['V_08_COS', 'Meridional wind  8hr. cos ' +
                                 'coeff.', 0, 'GDZ', 'sph',
                                 ['time', 'lon', 'lat', 'height'], 'm/s'],
                  'V_08_SIN': ['V_08_SIN_ilev', 'Meridional wind  8hr. sin ' +
                               'coeff.', 0, 'GDZ', 'sph',
                               ['time', 'lon', 'lat', 'ilev'], 'm/s'],
                  'V_08_SIN_2': ['V_08_SIN', 'Meridional wind  8hr. sin ' +
                                 'coeff.', 0, 'GDZ', 'sph',
                                 ['time', 'lon', 'lat', 'height'], 'm/s'],
                  'V_12_COS': ['V_12_COS_ilev', 'Meridional wind 12hr. cos ' +
                               'coeff.', 0, 'GDZ', 'sph',
                               ['time', 'lon', 'lat', 'ilev'], 'm/s'],
                  'V_12_COS_2': ['V_12_COS', 'Meridional wind 12hr. cos ' +
                                 'coeff.', 0, 'GDZ', 'sph',
                                 ['time', 'lon', 'lat', 'height'], 'm/s'],
                  'V_12_SIN': ['V_12_SIN_ilev', 'Meridional wind 12hr. sin ' +
                               'coeff.', 0, 'GDZ', 'sph',
                               ['time', 'lon', 'lat', 'ilev'], 'm/s'],
                  'V_12_SIN_2': ['V_12_SIN', 'Meridional wind 12hr. sin ' +
                                 'coeff.', 0, 'GDZ', 'sph',
                                 ['time', 'lon', 'lat', 'height'], 'm/s'],
                  'V_24_COS': ['V_24_COS_ilev', 'Meridional wind 24hr. cos ' +
                               'coeff.', 0, 'GDZ', 'sph',
                               ['time', 'lon', 'lat', 'ilev'], 'm/s'],
                  'V_24_COS_2': ['V_24_COS', 'Meridional wind 24hr. cos ' +
                                 'coeff.', 0, 'GDZ', 'sph',
                                 ['time', 'lon', 'lat', 'height'], 'm/s'],
                  'V_24_SIN': ['V_24_SIN_ilev', 'Meridional wind 24hr. sin ' +
                               'coeff.', 0, 'GDZ', 'sph',
                               ['time', 'lon', 'lat', 'ilev'], 'm/s'],
                  'V_24_SIN_2': ['V_24_SIN', 'Meridional wind 24hr. sin ' +
                                 'coeff.', 0, 'GDZ', 'sph',
                                 ['time', 'lon', 'lat', 'height'], 'm/s']
                  }

# lists of known variables per file type
# assume h1 file is always present
gvar_keys = {'h0': ['RHO_CLUBB_lev'],  # also containes km_ilev
             'h1': ['co2vmr', 'ch4vmr', 'n2ovmr', 'f11vmr', 'f12vmr',
                    'sol_tsi', 'colat_crit1', 'colat_crit2', 'ED1', 'ED2',
                    'EDYN_ZIGM11_PED', 'EDYN_ZIGM2_HAL', 'ElecColDens', 'H',
                    'O', 'O2', 'OMEGA', 'PHIM2D', 'PS', 'T', 'TElec', 'TIon',
                    'U', 'UI', 'V', 'VI', 'WI', 'Z3', 'e', 'RHO_CLUBB',
                    'Z3GM'],
             # keys in h2 files that are not in the h1 files
             'h2': ['ED1', 'ED2', 'EDens', 'ElecColDens', 'HMF2', 'NMF2',
                    'OPLUS', 'TElec', 'TIon', 'UI', 'VI', 'WI', 'Z3'],
             # keys in h3 files that are not in the h1 files
             'h3': ['CO2', 'EDens', 'EKGW', 'N', 'NO', 'OpDens', 'QCO2',
                    'QHC2S', 'QJOULE', 'QNO', 'QO3', 'QO3P', 'QRS_TOT',
                    'SolIonRate_Tot', 'TTGW', 'UTGW_TOTAL', 'Z3GM'],
             # keys in the h4 files that are not in the h1 or h3 files
             'h4': ['OMEGA_08_COS', 'OMEGA_08_SIN', 'OMEGA_12_COS',
                    'OMEGA_12_SIN', 'OMEGA_24_COS', 'OMEGA_24_SIN',
                    'T_08_COS', 'T_08_SIN', 'T_12_COS', 'T_12_SIN',
                    'T_24_COS', 'T_24_SIN', 'U_08_COS', 'U_08_SIN',
                    'U_12_COS', 'U_12_SIN', 'U_24_COS', 'U_24_SIN',
                    'V_08_COS', 'V_08_SIN', 'V_12_COS', 'V_12_SIN',
                    'V_24_COS', 'V_24_SIN']}
# set up pressure level and height dependent equivalent variables
ilev_list = [value[0] for key, value in model_varnames.items()
             if value[5][-1] == 'ilev' if value not in
             ['H_geopot_ilev', 'H_geomet_ilev']]
ilev_replace = [item.split('_ilev')[0] for item in ilev_list if
                item not in ['H_geopot_ilev', 'H_geomet_ilev']]


# times from file converted to seconds since midnight of filedate
# plotting input times will be datetime strings of format 'YYYY-MM-DD HH:mm:ss'
# filedate is self.filedate from iri object
# converts to hours since midnight of filedate for plotting
def MODEL():

    from kamodo import Kamodo
    from glob import glob
    from os.path import isfile, getsize
    import psutil
    from netCDF4 import Dataset
    from numpy import array, NaN, zeros, where, unique, append, linspace
    from numpy import transpose, median
    from time import perf_counter
    from datetime import datetime, timezone
    import kamodo_ccmc.readers.reader_utilities as RU
    from waccmx_tocdf import convert_all

    class MODEL(Kamodo):
        '''WACCM-X model data reader.

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
            - WACCM-X model outputs are given in netCDF files, so no file
              conversion is needed. However, some of the variables depend on a
              pressure level coordinate, so some pre-processing is needed to
              calculate the median height across the entire dataset
              (time, lon, lat) for a given pressure level value for each type
              of pressure level coordinate (typically two).
            - WACCM-X data is given in several coordinate systems depending
              on pressure level - a representation of height in Pascals
              corresponding to a given atmospheric pressure. The preset values
              of pressure level in the coordinate systems are not
              guaranteed to be identical, so they are inverted independently of
              the other unless only one is provided. In that case, the are
              assumed to be identical.
            - Pressure level inversion is performed in the reader_utilities
              script, specifically the PLevelInterp function. Two versions of
              all variables that depend on pressure level are created: the
              original and one dependent on height, which is created through
              Kamodo's function composition feature.
            - Interpolation method 0 is chosen for the spatially independent
              variables, meaning the entire time series is read into memory and
              the standard SciPy interpolator is used.
            - All of the files have multiple time steps per file, but some of
              the files are small, while some are larger than 16 GB. This
              motivates some logic to automatically choose the appropriate
              interpolation method based on the file size compared to the
              available computer memory. Interpolation method 3 is chosen if
              the file size is large, and interpolation method 2 is chosen if
              the file size is manageable. In both cases, the standard SciPy
              interpolator is used.
        '''
        def __init__(self, file_dir, variables_requested=[],
                     printfiles=False, filetime=False, gridded_int=True,
                     verbose=False, **kwargs):
            super(MODEL, self).__init__(**kwargs)
            self.modelname = 'WACCMX'
            t0 = perf_counter()

            # first, check for file list, create if DNE
            list_file = file_dir + self.modelname + '_list.txt'
            time_file = file_dir + self.modelname + '_times.txt'
            self.times, self.pattern_files = {}, {}
            if not isfile(list_file) or not isfile(time_file):
                # collect filenames
                files = sorted(glob(file_dir+'*.h?.*.nc'))
                patterns = unique([file[-22:-20] for file in files])
                self.filename = ''.join([f+',' for f in files])[:-1]
                self.filedate = datetime.strptime(
                    files[0][-19:-9]+' 00:00:00', '%Y-%m-%d %H:%M:%S'
                    ).replace(tzinfo=timezone.utc)

                # check for h0 files, contains RHO_CLUBB_lev and km_ilev
                if 'h0' not in patterns:
                    self.conversion_test = convert_all(file_dir)
                    if not self.conversion_test:
                        return
                    patterns = append(patterns, 'h0')

                # establish time attributes
                for p in patterns:
                    # get list of files to loop through later
                    pattern_files = sorted(glob(file_dir+'*.'+p+'.*.nc'))
                    self.pattern_files[p] = pattern_files
                    self.times[p] = {'start': [], 'end': [], 'all': []}

                    # loop through to get times
                    for f in range(len(pattern_files)):
                        cdf_data = Dataset(pattern_files[f])
                        # hours since midnight first file
                        tmp = array(cdf_data.variables['time'])
                        pmt = array(cdf_data.variables['datesec'])
                        if f == 0:
                            int_time = int(tmp[0])  # date of first file
                        net = (tmp-int_time).astype(int)*24. + pmt/3600.  # hrs
                        self.times[p]['start'].append(net[0])
                        self.times[p]['end'].append(net[-1])
                        self.times[p]['all'].extend(net)
                        cdf_data.close()
                    # adjust all to be arrays
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

            # collect variable list per file
            self.gvar_dict, gvar_list, self.err_list = {}, [], []
            # dict[key] = [gvar_list, file] for each type of file
            self.total_ilev = [item for item in ilev_list if item
                               not in ['H_geopot_ilev', 'H_geomet_ilev']]
            self.total_replace = ilev_replace
            # dictionary mapping to navigate related variable names
            self.ilev_map = {item1: item2 for item1, item2 in
                             zip(self.total_replace, self.total_ilev)}
            if len(variables_requested) > 0 and variables_requested != 'all':
                # add variables from each files if the file exists and if
                # variables are requested from that file.
                for key in gvar_keys.keys():
                    self.gvar_dict[key] = self.somefile_variables(
                        key, gvar_keys[key], variables_requested,
                        self.pattern_files)
                    gvar_list += self.gvar_dict[key][0]
                net_err = [value[0] for key, value in model_varnames.items()
                           if key not in gvar_list and value[0] in
                           self.err_list]
                if len(net_err) > 0:
                    print('Some requested variables are not available: ',
                          net_err)
            else:  # collect defaults
                # collect lists per file type
                for key in gvar_keys.keys():
                    self.gvar_dict[key] = self.allfile_variables(
                        key, gvar_keys[key], self.pattern_files)
                    gvar_list += self.gvar_dict[key][0]

                # returns list of variables included in data files
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
                    return

            # store data for each variable desired and all coordinate grids
            # need to leave files open for variable data
            self.variables = {}
            for key in self.pattern_files.keys():
                # collect dimensions from first file
                cdf_datah = Dataset(self.pattern_files[key][0])
                for dim in cdf_datah.dimensions.keys():
                    # deal with dim renaming and skipping
                    if dim == 'ilev':
                        name = 'ilev1'
                    elif dim == 'lev':
                        name = 'ilev'
                    elif dim in ['lon', 'lat', 'mlon', 'mlat']:
                        name = dim
                    else:  # skip other extra string dimensions
                        continue  # time already calculated, don't overwrite

                    # prep for longitude wrapping
                    tmp = array(cdf_datah.variables[dim])
                    if (dim == 'lon'):
                        lon_le180 = list(where(tmp <= 180)[0])
                        lon_ge180 = list(where(tmp >= 180)[0])  # repeat 180
                        lon_idx = lon_ge180+lon_le180
                        setattr(self, '_lon_idx_'+key, lon_idx)
                        out = append(tmp, 360.) - 180.
                    elif dim == 'mlon':
                        lon_idx = append(linspace(0, len(tmp)-1, len(tmp),
                                                  dtype=int), 0)
                        out = append(tmp, 180.)
                        setattr(self, '_mlon_idx_'+key, lon_idx)
                    else:  # loop through other dimensions without changes
                        out = tmp
                    # store coordinate grids
                    setattr(self, '_'+name+'_'+key, out)
                cdf_datah.close()  # close for dimensions data since now arrays

                # loop through h0 files to calculate median km_ilev grid
                if 'h0' in key:
                    km_max, km_min, km = [], [], []
                    for f in self.pattern_files[key]:
                        cdf_datah0 = Dataset(f)
                        km.append(array(cdf_datah0.variables['km_ilev']))
                        km_max.append(cdf_datah0.km_ilev_max)
                        km_min.append(cdf_datah0.km_ilev_min)
                        cdf_datah0.close()
                    self._km_ilev = median(array(km), axis=0)
                    self._km_ilev_max = array(km_max).max()
                    self._km_ilev_min = array(km_min).min()

                # initialize variables dictionary
                for var in self.gvar_dict[key][0]:
                    self.variables[model_varnames[var][0]] = {
                        'units': model_varnames[var][-1], 'data': key}

            # store a few items
            self.missing_value = NaN
            if verbose:
                print(f'Took {perf_counter()-t0:.6f}s to read in data')
            if printfiles:
                print(self.filename)

            # register interpolators for each requested variable
            # rearrange to deal with H_ilev and H_ilev1 first if there
            varname_list = [key for key in self.variables.keys()]
            if 'H_geopot_ilev' in varname_list:
                varname_list.remove('H_geopot_ilev')
                varname_list = ['H_geopot_ilev'] + varname_list
            t_reg = perf_counter()
            for varname in varname_list:
                # register the variables
                self.register_variable(varname, gridded_int)
            if verbose:
                print(f'Took {perf_counter()-t_reg:.5f}s to register ' +
                      f'{len(varname_list)} variables.')
            if verbose:
                print(f'Took a total of {perf_counter()-t0:.5f}s to kamodofy' +
                      f' {len(varname_list)} variables.')

        # variable checks
        def variable_checks(self, variables_requested, g_varkeys_h,
                            cdf_datah_keys):
            '''Add pressure level variables if needed and perform variable
            checks for case where variables_requested is not empty.'''

            add_ilev = [var+'_ilev' for var in variables_requested if var
                        in self.total_replace]
            if len(add_ilev) > 0:
                add_ilev += ['H_geopot_ilev']
            new_var = variables_requested + add_ilev
            # remove replaced items and any items not expected to be in hx file
            short_var = [value[0] for key, value in model_varnames.items()
                         if value[0] in new_var and value[0] not in
                         self.total_replace and key in g_varkeys_h]
            gvar_listh = [key for key, value in model_varnames.items()
                          if value[0] in short_var and key
                          in cdf_datah_keys]

            # check for variables requested but not available
            if len(gvar_listh) != len(short_var):
                err_list = [value[0] for key, value in
                            model_varnames.items() if value[0] in
                            short_var and key not in cdf_datah_keys]
                for item in err_list:  # need a flat list
                    self.err_list.append(item)
                # add items to list, check full list later
            return gvar_listh

        # select desired variables from given file type
        def somefile_variables(self, h_type, g_varkeys_h, variables_requested,
                               file_dict):
            '''Add variables from desired file if file exists and if variables
            are requested.
            Inputs:
                - h_type is one of 'h2', 'h3', 'h4'.
                - g_varkeysh is a list of variables names known to exist in the
                    given file type.
                - variables_requested: a list of desired LaTeX variable names.
                - hx_files: a list of files for a given day.
            '''
            if h_type in file_dict.keys():
                cdf_datah = Dataset(file_dict[h_type][0])
                cdf_datah_keys = list(cdf_datah.variables.keys())
                cdf_datah.close()
                gvar_listh = self.variable_checks(variables_requested,
                                                  g_varkeys_h, cdf_datah_keys)
                return [gvar_listh, file_dict[h_type]]
            else:  # collect error list
                gvar_listh = self.variable_checks(variables_requested,
                                                  g_varkeys_h, [])
                return [gvar_listh, '']

        # return list of variables available in chosen file type
        def allfile_variables(self, h_type, g_varkeys_h, file_dict):
            '''Return list of variables from desired file if file exists.
            Inputs:
                - h_type is one of 'h2', 'h3', 'h4'.
                - g_varkeysh is a list of variables names known to exist in the
                    given file type.
            '''
            if h_type in file_dict.keys():
                cdf_datah = Dataset(file_dict[h_type][0])
                gvar_listh = [key for key, value in model_varnames.items()
                              if key in cdf_datah.variables.keys() and key in
                              g_varkeys_h]
                cdf_datah.close()
                return [gvar_listh, file_dict[h_type]]
            else:
                return [[], '']

        # dimension agnostic registration
        def register_variable(self, varname, gridded_int):
            '''Registers and functionalizes the dataset. Converts to a 1D
            numpy array for time series data.'''

            # create the coordinate dictionary
            key = self.variables[varname]['data']
            gvar = [key for key, value in model_varnames.items() if
                    value[0] == varname][0]  # variable name in file
            coord_list = [value[-2] for key, value in model_varnames.items()
                          if value[0] == varname][0]
            coord_str = [value[3]+value[4] for key, value in
                         model_varnames.items() if value[0] == varname][0]
            coord_dict = {'time': {'units': 'hr',
                                   'data': self.times[key]['all']}}
            if len(coord_list) == 1:  # convert 1D to array and functionalize
                data = []
                for f in self.pattern_files[key]:
                    cdf_data = Dataset(f)
                    tmp = array(cdf_data.variables[gvar])
                    cdf_data.close()
                    data.extend(tmp)
                self.variables[varname]['data'] = array(data)
                self = RU.Functionalize_Dataset(self, coord_dict, varname,
                                                self.variables[varname],
                                                gridded_int, coord_str,
                                                interp_flag=0)  # single array
                return
            for coord in coord_list[1:]:
                coord_dict[coord] = {'data': getattr(self,
                                                     '_'+coord+'_'+key)}
                if coord in ['lat', 'lon', 'mlon', 'mlat']:
                    coord_dict[coord]['units'] = 'deg'
                if coord in ['ilev', 'ilev1']:
                    coord_dict[coord]['units'] = 'hPa'
            lon_idx = getattr(self, '_'+coord_list[1]+'_idx_'+key)

            # overshooting guess at memory needed to read in all data
            # chunks/slices will be dropped if memory runs out later
            memory_needed = getsize(self.pattern_files['h1'][0]) * \
                len(self.pattern_files[key])
            memory_test = memory_needed < psutil.virtual_memory().available
            if varname == 'H_geopot_ilev':
                gridded_int = True
            if not memory_test:  # varname != 'H_geopot_ilev' and not
                # save memory by only retrieving time slices from the files
                # determine the operation to occur on each time slice
                def func(i, fi):  # i = file#, fi = slice#
                    cdf_data = Dataset(self.pattern_files[key][i])
                    tmp = array(cdf_data.variables[gvar][fi])
                    cdf_data.close()
                    # perform data wrangling on time slice
                    # (ilev,) lat/mlat, lon/mlon -> lon/mlon, lat/mlat, (ilev)
                    data = tmp.T[lon_idx]
                    return data

                # functionalize the 3D or 4D dataset, series of time slices
                self = RU.Functionalize_Dataset(
                    self, coord_dict, varname, self.variables[varname],
                    gridded_int, coord_str, interp_flag=3, func=func,
                    times_dict=self.times[key])
            else:  # only using this approach if there is enough memory
                # also using for H_ilev to speed up function composition
                # determine the operation to occur on each time chunk
                def func(f):  # file_idx
                    # convert slice number s to index in file
                    cdf_data = Dataset(self.pattern_files[key][f])
                    tmp = array(cdf_data.variables[gvar])
                    cdf_data.close()
                    # perform data wrangling on time chunk
                    # (ilev,) lat/mlat, lon/mlon -> lon/mlon, lat/mlat, (ilev)
                    if len(tmp.shape) == 4:
                        data = transpose(tmp, (0, 3, 2, 1))
                    elif len(tmp.shape) == 3:
                        data = transpose(tmp, (0, 2, 1))
                    return data[:, lon_idx]

                self = RU.Functionalize_Dataset(
                    self, coord_dict, varname, self.variables[varname],
                    gridded_int, coord_str, interp_flag=2, func=func,
                    times_dict=self.times[key])

            if len(coord_list) < 4:  # remaining logic is for 4D data.
                return

            # create pressure level -> km function once per ilev type
            if varname == 'H_geopot_ilev' or varname in self.total_ilev:
                if varname == 'H_geopot_ilev':
                    new_varname, units = 'Plev', 'hPa'
                    # perform unit conversion to km with Kamodo and invert
                    self['H_ilev_ijk[km]'] = 'H_geopot_ilev_ijk'
                    self[new_varname], interp_ijk = RU.PLevelInterp(
                        self['H_ilev_ijk'], coord_dict['time']['data'],
                        coord_dict['lon']['data'], coord_dict['lat']['data'],
                        coord_dict['ilev']['data'], units, self._km_ilev,
                        [self._km_ilev_max, self._km_ilev_min])
                    # max pressure level is at lowest height, so flip order
                elif varname != 'H_geopot_ilev':  # define by function comp
                    new_varname = varname.split('_ilev')[0]
                    units = self.variables[varname]['units']
                    self[new_varname] = varname+'(Plev)'
                    interp_ijk = self[new_varname]
                    self[new_varname].meta['arg_units'] = \
                        self['Plev'].meta['arg_units']
                self.variables[new_varname] = {'units': units, 'data': key}

                # create gridded interpolator if requested
                if gridded_int:
                    self = RU.register_griddedPlev(
                        self, new_varname, units, interp_ijk, coord_dict,
                        self._km_ilev, [self._km_ilev_max, self._km_ilev_min])
                    # max pressure level is at lowest height, so flip order
            return
    return MODEL
