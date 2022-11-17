# !/usr/bin/env python
# -----------------------------------------------------------------------------
# $Id: gitm.py,v 1.16 2014/04/09 10:44:28 agburr Exp $
#
# GITM.py, Dan Welling, UMich
#
# Comments: Defines a class for GITM binary output files and modifies the input
#           to improve the data analysis and plotting experience
#
# Contains: class GitmBin    - The class for the GITM binary, which will read
#                              a single GITM output binary
#           def calc_magdi   - Reads a single GITM ion or mag output binary and
#                              computes the magnetic inclination and
#                              declination
#           def calc_magvel  - Reads a single GITM ion output binary and
#                              uses data from the standard GITM output (3DAll)
#                              to compute the ion characteristics in magnetic
#                              coordinates
#           def calc_deg     - Computes and appends latitude and longitude
#                              in degrees from values in radians
#           def calc_lt      - Computes and appends local time in hours from
#                              the universal time for the file and longitude
#           def append_units - Appends unit, descriptive name, and scale
#                              attributes to each known data type
#           def append_data  - Appends a list of data variables to a GitmBin
#                              data structure, where only a limited number of
#                              data variables from that file have been read
#                              in before
#           def calc_tec     - Calculate the VTEC
#           def calc_2dion   - Calculate the 2D ionospheric parameters (VTEC,
#                              hmF2, NmF2)
#
# Updates:
#          Angeline Burrell (AGB) - 1/7/13: Added calc_lt, append_units, and
#                                           calc_magvel
#          AGB - 11/7/13: Improved calc_2dion, added Aaron Ridley's calc_tec
#          AGB - 12/6/13: Added Inclination/Declination calculation
#          Darren De Zeeuw (DDZ) - 06/24/19: Updated code to python3,
#                 Aaron Ridley approved reader for open source use in Kamodo
#          Rebecca Ringuette - 05/14/2021: Removed spacepy dependency and
#                 converted output to a netCDF4 file for increased speed
#                 in later layers. Also removed and adapted calc_magdi,
#                 calc_magvel,
#                 calc_lt for generalized use by other readers. Copied and
#                 adapted calc_tec and calc_2dion for use by other readers.
#                 Removed append_data since converting whole file to a cdf.
# -----------------------------------------------------------------------------


# Global imports:
from glob import glob
import numpy as np
from time import perf_counter
from datetime import datetime, timezone
from netCDF4 import Dataset
from kamodo import Kamodo
from kamodo_ccmc.readers.reader_utilities import Functionalize_Dataset


# replace weird name with standard names in cdf file
gitm_varnames = {'Argon Mixing Ratio': ['r_Ar', 'linear', ''],
                 'Ar Mass Density': ['rho_Ar', 'exponential', 'kg/m**3'],
                 'Methane Mixing Ratio': ['r_CH4', 'linear', ''],
                 'Conduction': ['k', 'linear', 'W/m/K'],
                 'EUV Heating': ['Q_EUV', 'linear', 'K per timestep'],
                 'H Mass Density': ['rho_H', 'exponential', 'kg/m**3'],
                 'H+ Mass Density': ['rho_Hplus', 'exponential', 'kg/m**3'],
                 'H2 Mixing Ratio': ['r_H2', 'linear', ''],
                 'Hydrogen Cyanide Mixing Ratio': ['r_HCN', 'linear', ''],
                 'He Mass Density': ['rho_He', 'exponential', 'kg/m**3'],
                 'He+ Mass Density': ['rho_Heplus', 'exponential', 'kg/m**3'],
                 'Heating Efficiency': ['HeatingEfficiency', 'linear', ''],
                 'Heat Balance Total': ['HeatBalanceTotal', 'linear', ''],
                 'N2 Mass Density': ['rho_N2', 'exponential', 'kg/m**3'],
                 'N2+ Mass Density': ['rho_N2plus', 'exponential', 'kg/m**3'],
                 'N+ Mass Density': ['rho_Nplus', 'exponential', 'kg/m**3'],
                 'N(2D) Mass Density': ['rho_N2D', 'exponential', 'kg/m**3'],
                 'N(2P) Mass Density': ['rho_N2P', 'exponential', 'kg/m**3'],
                 'N(4S) Mass Density': ['rho_N4S', 'exponential', 'kg/m**3'],
                 'N2 Mixing Ratio': ['r_N2', 'linear', ''],
                 'NO Mass Density': ['rho_NO', 'exponential', 'kg/m**3'],
                 'NO+ Mass Density': ['rho_NOplus', 'exponential', 'kg/m**3'],
                 'O2 Mass Density': ['rho_O2', 'exponential', 'kg/m**3'],
                 'O(1D) Mass Density': ['rho_O1D', 'exponential', 'kg/m**3'],
                 'O2+ Mass Density': ['rho_O2plus', 'exponential', 'kg/m**3'],
                 'O(2D) Mass Density': ['rho_O2D', 'exponential', '1/m**3'],
                 'O+(2P) Mass Density': ['rho_Oplus2P', 'exponential',
                                         'kg/m**3'],
                 'O(3P) Mass Density': ['rho_O3P', 'exponential', 'kg/m**3'],
                 'O+(4SP) Mass Density': ['rho_Oplus4SP', 'exponential',
                                          'kg/m**3'],
                 'Radiative Cooling': ['L_Rad', 'linear', ''],
                 'Neutral Density': ['rho', 'exponential', 'kg/m**3'],
                 'T_n': ['T_n', 'linear', 'K'],
                 'vi_east': ['vi_east', 'linear', 'm/s'],
                 'vi_north': ['vi_north', 'linear', 'm/s'],
                 'vi_up': ['vi_up', 'linear', 'm/s'],
                 'vn_east': ['vn_east', 'linear', 'm/s'],
                 'vn_north': ['vn_north', 'linear', 'm/s'],
                 'vn_up': ['vn_up', 'linear', 'm/s'],
                 'v_N2_up': ['v_N2_up', 'linear', 'm/s'],
                 'v_N(4S)_up': ['v_N4S_up', 'linear', 'm/s'],
                 'v_N_up': ['v_N_up', 'linear', 'm/s'],
                 'v_O2_up': ['v_O2_up', 'linear', 'm/s'],
                 'v_O(3P)_up': ['v_O3P_up', 'linear', 'm/s'],
                 'v_He_up': ['v_He_up', 'linear', 'm/s'],
                 '[e-]': ['rho_e', 'linear', '1/m**3'],
                 'Electron Average Energy': ['ElectronAverageEnergy', 'linear',
                                             'J'],
                 'T_e': ['T_e', 'linear', 'K'],
                 'T_i': ['T_i', 'linear', 'K'],
                 'Solar Zenith Angle': ['SolarZenithAngle', 'linear',
                                        'radians'],
                 'CO2 Mass Density': ['rho_CO2', 'exponential', 'kg/m**3'],
                 'DivJu FL': ['DivJuFL', '', ''],
                 'DivJuAlt': ['DivJuAlt', 'linear', ''],
                 'Electron Energy Flux': ['ElectronEnergyFlux', 'exponential',
                                          'J/m**2'],
                 'Field Line Length': ['Field Line Length', 'linear', 'm'],
                 'sigma_P': ['sigma_P', 'linear', 'S/m'],
                 'Sigma_P': ['Sigma_P', 'linear', 'S/m'],
                 'sigma_H': ['sigma_H', 'linear', 'S/m'],
                 'Potential': ['V', 'linear', 'V'],
                 'Sigma_H': ['Sigma_H', 'linear', 'S/m'],
                 'Region 2 Current': ['I_R2', 'linear', 'A/m**2'],
                 'Region 1 Current': ['I_R1', 'linear', 'A/m**2'],
                 'Ed1': ['Ed1', 'linear', ''],
                 'Ed2': ['Ed2', 'linear', ''],
                 'Solar Local Time': ['SolarLocalTime', 'linear', 'h'],
                 'Vertical Electric Field': ['E_up', 'linear', 'V/m'],
                 'Eastward Electric Field': ['E_east', 'linear', 'V/m'],
                 'Northward Electric Field': ['E_north', 'linear', 'V/m'],
                 'Electric Field Magnitude': ['E_mag', 'linear', 'V/m'],
                 'Vertical Magnetic Field': ['B_up', 'linear', 'nT'],
                 'Eastward Magnetic Field': ['B_east', 'linear', 'nT'],
                 'Northward Magnetic Field': ['B_north', 'linear', 'nT'],
                 'Magnetic Field Magnitude': ['B_mag', 'linear', 'nT'],
                 'Magnetic Latitude': ['MagLat', 'linear', 'deg'],
                 'Magnetic Longitude': ['MagLon', 'linear', 'deg'],
                 'g': ['g', 'linear', 'm/s**2'],
                 'GradP_east (P_i + P_e)': ['GradP_east', 'linear', 'Pa/m'],
                 'GradP_north (P_i + P_e)': ['GradP_north', 'linear', 'Pa/m'],
                 'GradP_up (P_i + P_e)': ['GradP_up', 'linear', 'Pa/m'],
                 'nu_in': ['nu_in', 'linear', '1/s'],
                 'Chemical Heating Rate': ['ChemicalHeatingRate', 'linear',
                                           ''],
                 'Total Absolute EUV': ['TotalAbsoluteEUV', 'linear',
                                        'K per timestep'],
                 'O Cooling': ['L_O', 'linear', 'K per timestep'],
                 'Joule Heating': ['Q_Joule', 'linear', 'K per timestep'],
                 'Auroral Heating': ['Q_Auroral', 'linear', 'K per timestep'],
                 'Photoelectron Heating': ['Q_PhotoE', 'linear',
                                           'K per timestep'],
                 'Eddy Conduction': ['k_eddy', 'linear', ''],
                 'Adiabatic Eddy Conduction': ['k_adiabaticeddy', 'linear',
                                               ''],
                 'NO Cooling': ['L_NO', 'linear', 'K per timestep'],
                 'Molecular Conduction': ['k_molecular', 'linear', ''],
                 'max_electron_density': ['NmF2', 'linear', ''],
                 'max_electron_density_height': ['hmF2', 'linear', 'km'],
                 'Vertical TEC': ['TEC', 'linear', '10**16/m**2'],
                 'HeatFlux_Joule': ['phi_qJoule', 'linear', 'W/m**2'],
                 'HeatFlux': ['phi_q', 'linear', 'W/m**2'],
                 'HeatFlux_EUV': ['phi_qEUV', 'linear', 'W/m**2'],
                 'NO CoolingFlux': ['phi_qNOCooling', 'linear', 'W/m**2']}

# if not always one value, deal with this in kamodo wrapper
name_dict = {"Altitude": "Altitude", "Ar Mixing Ratio": "Argon Mixing Ratio",
             "Ar": "Ar Mass Density",
             "CH4 Mixing Ratio": "Methane Mixing Ratio",
             "Conduction": "Conduction", "EuvHeating": "EUV Heating",
             "H": "H Mass Density", "H!U+!N": "H+ Mass Density",
             "H2 Mixing Ratio": "H2 Mixing Ratio",
             "HCN Mixing Ratio": "Hydrogen Cyanide Mixing Ratio",
             "He": "He Mass Density", "He!U+!N": "He+ Mass Density",
             "Heating Efficiency": "Heating Efficiency",
             "Heat Balance Total": "Heat Balance Total",
             "Latitude": "Latitude", "Longitude": "Longitude",
             "N!D2!N": "N2 Mass Density", "N!D2!U+!N": "N2+ Mass Density",
             "N!U+!N": "N+ Mass Density", "N(!U2!ND)": "N(2D) Mass Density",
             "N(!U2!NP)": "N(2P) Mass Density",
             "N(!U4!NS)": "N(4S) Mass Density",
             "N2 Mixing Ratio": "N2 Mixing Ratio", "NO": "NO Mass Density",
             "NO!U+!N": "NO+ Mass Density", "O!D2!N": "O2 Mass Density",
             "O(!U1!ND)": "O(1D) Mass Density",
             "O!D2!U+!N": "O2+ Mass Density",
             "O(!U2!ND)!": "O(2D) Mass Density",
             "O(!U2!ND)!U+!N": "O(2D) Mass Density",
             "O(!U2!NP)!U+!N": "O+(2P) Mass Density",
             "O(!U3!NP)": "O(3P) Mass Density",
             "O_4SP_!U+!N": "O+(4SP) Mass Density",
             "RadCooling": "Radiative Cooling",
             "Rho": "Neutral Density", "Temperature": "T_n",
             "V!Di!N (east)": "vi_east",
             "V!Di!N (north)": "vi_north", "Vi (north) (m/s)": "vi_north",
             "V!Di!N (up)": "vi_up", "Vi (up) (m/s)": "vi_up",
             "V!Dn!N (east)": "vn_east",
             "Vn (east) (m/s)": "vn_east", "V!Dn!N (north)": "vn_north",
             "Vn (north) (m/s)": "vn_north", "V!Dn!N (up)": "vn_up",
             "Vn (up) (m/s)": "vn_up",
             "V!Dn!N (up,N!D2!N              )": "v_N2_up",
             "V!Dn!N (up,N(!U4!NS)           )": "v_N(4S)_up",
             "V!Dn!N (up,NO                  )": "v_N_up",
             "V!Dn!N (up,O!D2!N              )": "v_O2_up",
             "V!Dn!N (up,O(!U3!NP)           )": "v_O(3P)_up",
             "V!Dn!N (up,He                  )": "v_He_up",
             "e-": "[e-]",
             "Electron_Average_Energy": "Electron Average Energy",
             "eTemperature": "T_e", "iTemperature": "T_i",
             "Solar Zenith Angle": "Solar Zenith Angle",
             "Vertical TEC": "Vertical TEC", "CO!D2!N": "CO2 Mass Density",
             "DivJu FL": "DivJu FL", "DivJuAlt": "DivJuAlt",
             "Electron_Energy_Flux": "Electron Energy Flux",
             "FL Length": "Field Line Length",
             "Pedersen FL Conductance": "sigma_P",
             "Pedersen Conductance": "Sigma_P",
             "Hall FL Conductance": "sigma_H",
             "Potential": "Potential", "Hall Conductance": "Sigma_H",
             "Je2": "Region 2 Current", "Je1": "Region 1 Current",
             "Ed1": "Ed1", "Ed2": "Ed2", "LT": "Solar Local Time",
             "Local Time": "Solar Local Time",
             "E.F. Vertical": "Vertical Electric Field",
             "E.F. East": "Eastward Electric Field",
             "E.F. North": "Northward Electric Field",
             "E.F. Magnitude": "Electric Field Magnitude",
             "B.F. Vertical": "Vertical Magnetic Field",
             "B.F. East": "Eastward Magnetic Field",
             "B.F. North": "Northward Magnetic Field",
             "B.F. Magnitude": "Magnetic Field Magnitude",
             "Magnetic Latitude": "Magnetic Latitude",
             "Magnetic Longitude": "Magnetic Longitude", "dLat": "Latitude",
             "dLon": "Longitude",
             "Gravity": "g", "PressGrad (east)": "GradP_east (P_i + P_e)",
             "PressGrad (north)": "GradP_north (P_i + P_e)",
             "PressGrad (up)": "GradP_up (P_i + P_e)",
             "IN Collision Freq": "nu_in",
             "Chemical Heating": "Chemical Heating Rate",
             "Total Abs EUV": "Total Absolute EUV",
             "O Cooling": "O Cooling", "Joule Heating": "Joule Heating",
             "Auroral Heating": "Auroral Heating",
             "Photoelectron Heating": "Photoelectron Heating",
             "Eddy Conduction": "Eddy Conduction",
             "Eddy Adiabatic Conduction": "Adiabatic Eddy Conduction",
             "NO Cooling": "NO Cooling",
             "Molecular Conduction": "Molecular Conduction",
             'NmF2': 'max_electron_density',
             'hmF2': "max_electron_density_height",
             'VTEC': "Vertical TEC",
             'AltIntJouleHeating (W/m2)': 'HeatFlux_Joule',
             'AltIntHeatingTransfer (W/m2)': 'HeatFlux',
             'AltIntEuvHeating (W/m2)': 'HeatFlux_EUV',
             'AltIntNOCooling (W/m2)': 'NO CoolingFlux'}

unit_dict = {"Altitude": "m", "Ar Mixing Ratio": "", "Ar": "kg/m**3",
             "CH4 Mixing Ratio": "", "Conduction": "W/m/K",
             "EuvHeating": "K per timestep",
             "H": "kg/m**3", "H!U+!N": "kg/m**3", "H2 Mixing Ratio": "",
             "HCN Mixing Ratio": "", "He": "kg/m**3", "He!U+!N": "kg/m**3",
             "Heating Efficiency": "", "Heat Balance Total": "",
             "Latitude": "radians",
             "Longitude": "radians", "N!D2!N": "kg/m**3",
             "N!D2!U+!N": "kg/m**3",
             "N!D2!U+!N            (/m3)": "1/m**3", "N!U+!N": "kg/m**3",
             "N(!U2!ND)": "kg/m**3", "N(!U2!NP)": "kg/m**3",
             "N(!U4!NS)": "kg/m**3",
             "N2 Mixing Ratio": "", "NO": "kg/m**3", "NO!U+!N": "kg/m**3",
             "O!D2!N": "kg/m**3", "O(!U1!ND)": "kg/m**3",
             "O!D2!U+!N": "kg/m**3",
             "O(!U2!ND)!": "1/m**3", "O(!U2!ND)!U+!N": "1/m**3",
             "O(!U2!NP)!U+!N": "kg/m**3",
             "O(!U2!NP)!U+!N": "kg/m**3", "O(!U3!NP)": "kg/m**3",
             "O_4SP_!U+!N": "kg/m**3",
             "RadCooling": "", "Rho": "kg/m**3", "Temperature": "K",
             "V!Di!N (east)": "m/s",
             "Vi (east) (m/s)": "m/s", "V!Di!N (north)": "m/s",
             "Vi (north) (m/s)": "m/s",
             "V!Di!N (up)": "m/s", "Vi (up) (m/s)": "m/s",
             "V!Dn!N (east)": "m/s",
             "Vn (east) (m/s)": "m/s", "V!Dn!N (north)": "m/s",
             "Vn (north) (m/s)": "m/s",
             "V!Dn!N (up)": "m/s", "Vn (up) (m/s)": "m/s",
             "V!Dn!N (up,N!D2!N              )": "m/s",
             "V!Dn!N (up,N(!U4!NS)           )": "m/s",
             "V!Dn!N (up,NO                  )": "m/s",
             "V!Dn!N (up,O!D2!N              )": "m/s",
             "V!Dn!N (up,O(!U3!NP)           )": "m/s",
             "V!Dn!N (up,He                  )": "m/s",
             "e-": "1/m**3", "Electron_Average_Energy": "J",
             "eTemperature": "K",
             "iTemperature": "K", "LT": "h", "Local Time": "h",
             "Solar Zenith Angle": "radians",
             "Vertical TEC": "10**16/m**2", "CO!D2!N": "kg/m**3",
             "DivJu FL": "", "DivJuAlt": "",
             "Electron_Energy_Flux": "J/m**2", "FL Length": "m",
             "Pedersen FL Conductance": "S/m",
             "dLon": "deg", "Pedersen Conductance": "S/m", "dLat": "deg",
             "Hall FL Conductance": "S/m", "Potential": "V",
             "Hall Conductance": "S/m",
             "Je2": "A/m**2", "Je1": "A/m**2", "Magnetic Longitude": "deg",
             "E.F. Vertical": "V/m", "E.F. East": "V/m", "E.F. North": "V/m",
             "E.F. Magnitude": "V/m",
             "B.F. Vertical": "nT", "B.F. East": "nT", "B.F. North": "nT",
             "B.F. Magnitude": "nT",
             "Magnetic Latitude": "deg", "Ed1": "", "Ed2": "",
             "Gravity": "m/s**2",
             "PressGrad (east)": "Pa/m", "Joule Heating": "K per timestep",
             "PressGrad (north)": "Pa/m", "O Cooling": "K per timestep",
             "PressGrad (up)": "Pa/m", "Total Abs EUV": "K per timestep",
             "IN Collision Freq": "1/s", "Chemical Heating": "",
             "Auroral Heating": "K per timestep",
             "Photoelectron Heating": "K per timestep",
             "Eddy Conduction": "", "Eddy Adiabatic Conduction": "",
             "NO Cooling": "K per timestep", "Molecular Conduction": "",
             'NmF2': "", 'hmF2': "km", 'VTEC': "10**16/m**2",
             'AltIntJouleHeating (W/m2)': 'W/m**2',
             'AltIntHeatingTransfer (W/m2)': 'W/m**2',
             'AltIntEuvHeating (W/m2)': 'W/m**2',
             'AltIntNOCooling (W/m2)': 'W/m**2'}

scale_dict = {"Altitude": "linear", "Ar Mixing Ratio": "linear",
              "Ar": "exponential", "CH4 Mixing Ratio": "linear",
              "Conduction": "linear", "EuvHeating": "linear",
              "H": "exponential", "H!U+!N": "exponential",
              "H2 Mixing Ratio": "linear", "HCN Mixing Ratio": "linear",
              "He": "exponential", "He!U+!N": "exponential",
              "Heating Efficiency": "linear", "DivJuAlt": "linear",
              "Heat Balance Total": "linear", "Latitude": "linear",
              "Longitude": "linear", "N!D2!N": "exponential",
              "N!D2!U+!N": "exponential",
              "N!U+!N": "exponential", "N(!U2!ND)": "exponential",
              "N(!U2!NP)": "exponential",
              "N(!U4!NS)": "exponential", "N2 Mixing Ratio": "linear",
              "NO": "exponential",
              "NO!U+!N": "exponential", "O!D2!N": "exponential",
              "O!D2!U+!N": "exponential",
              "O(!U2!ND)!": "exponential", "O(!U1!ND)": "exponential",
              "O(!U2!ND)!U+!N": "exponential", "CO!D2!N": "exponential",
              "O(!U2!NP)!U+!N": "exponential", "DivJu FL": "",
              "O(!U2!NP)!U+!N": "exponential", "O(!U3!NP)": "exponential",
              "O_4SP_!U+!N": "exponential", "RadCooling": "linear",
              "Rho": "exponential", "Temperature": "linear",
              "V!Di!N (east)": "linear",
              "Vi (east)": "linear", "V!Di!N (north)": "linear",
              "Vi (north) (m/s)": "linear",
              "V!Di!N (up)": "linear", "Vi (up) (m/s)": "linear",
              "V!Dn!N (east)": "linear",
              "Vn (east) (m/s)": "linear", "V!Dn!N (north)": "linear",
              "Vn (north) (m/s)": "linear",
              "V!Dn!N (up)": "linear", "Vn (up) (m/s)": "linear",
              "V!Dn!N (up,N!D2!N              )": "linear",
              "V!Dn!N (up,N(!U4!NS)           )": "linear",
              "V!Dn!N (up,NO                  )": "linear",
              "V!Dn!N (up,O!D2!N              )": "linear",
              "V!Dn!N (up,O(!U3!NP)           )": "linear",
              "V!Dn!N (up,He                  )": "linear",
              "e-": "linear",
              # "e-                   (/m3)": "linear",
              "Electron_Average_Energy": "linear",
              "eTemperature": "linear", "iTemperature": "linear",
              "Solar Zenith Angle": "linear", "Vertical TEC": "linear",
              "Electron_Energy_Flux": "exponential",
              "FL Length": "linear", "Pedersen FL Conductance": "linear",
              "Hall Conductance": "linear", "Potential": "linear",
              "Hall FL Conductance": "linear", "dLon": "linear",
              "Pedersen Conductance": "linear", "Je2": "linear",
              "Je1": "linear", "Ed1": "linear", "Ed2": "linear",
              "E.F. Vertical": "linear", "E.F. East": "linear",
              "E.F. North": "linear", "E.F. Magnitude": "linear",
              "B.F. Vertical": "linear", "B.F. East": "linear",
              "B.F. North": "linear", "B.F. Magnitude": "linear",
              "Magnetic Latitude": "linear", "LT": "linear",
              "Local Time": "linear", "Magnetic Longitude": "linear",
              "dLat": "linear",
              "Gravity": "linear", "PressGrad (east)": "linear",
              "PressGrad (north)": "linear", "PressGrad (up)": "linear",
              "IN Collision Freq": "linear", "Chemical Heating": "linear",
              "Total Abs EUV": "linear", "O Cooling": "linear",
              "Joule Heating": "linear", "Auroral Heating": "linear",
              "Photoelectron Heating": "linear", "NO Cooling": "linear",
              "Eddy Conduction": "linear", "Molecular Conduction": "linear",
              "Eddy Adiabatic Conduction": "linear",
              'NmF2': "linear", 'hmF2': "linear", 'VTEC': 'linear',
              'AltIntJouleHeating (W/m2)': 'linear',
              'AltIntHeatingTransfer (W/m2)': 'linear',
              'AltIntEuvHeating (W/m2)': 'linear',
              'AltIntNOCooling (W/m2)': 'linear'}


def _read(varlist, attrs, newfile=True):
    '''
    Read binary file.
    '''

    from re import sub
    from struct import unpack
    import sys

    # Read data and header info
    coords, variables = {}, {}
    f = open(attrs['file'], 'rb')

    # Using the first FORTRAN header, determine endian.
    # Default is little.
    attrs['endian'] = 'little'
    endChar = '>'
    rawRecLen = f.read(4)
    if len(rawRecLen) < 4:
        print("GitmBin ERROR: empty file [", attrs['file'], "]")
        sys.exit(1)
    recLen = (unpack(endChar+'l', rawRecLen))[0]
    if (recLen > 10000) or (recLen < 0):
        # Ridiculous record length implies wrong endian.
        attrs['endian'] = 'big'
        endChar = '<'
        recLen = (unpack(endChar+'l', rawRecLen))[0]

    # Read version; read fortran footer+header.
    attrs['version'] = unpack(endChar+'d', f.read(recLen))[0]
    (oldLen, recLen) = unpack(endChar+'2l', f.read(8))

    # Read grid size information.
    (attrs['nLon'], attrs['nLat'], attrs['nAlt']) = \
        unpack(endChar+'lll', f.read(recLen))
    (oldLen, recLen) = unpack(endChar+'2l', f.read(8))

    # Read number of variables.
    attrs['nVars'] = unpack(endChar+'l', f.read(recLen))[0]
    (oldLen, recLen) = unpack(endChar+'2l', f.read(8))

    # Collect variable names.
    var = []
    for i in range(attrs['nVars']):
        var.append(unpack(endChar+'%is' % (recLen), f.read(recLen))[0])
        (oldLen, recLen) = unpack(endChar+'2l', f.read(8))

    # Extract time.
    (yy, mm, dd, hh, mn, ss, ms) = unpack(endChar+'lllllll', f.read(recLen))
    # self._time=dt.datetime(yy,mm,dd,hh,mn,ss,ms/1000)
    attrs['time'] = datetime(yy, mm, dd, hh, mn, ss, ms).replace(
        tzinfo=timezone.utc).isoformat(sep=' ')
    (oldLen) = unpack(endChar+'l', f.read(4))

    # Read the rest of the data. Have to read through whole file.
    nTotal = attrs['nLon']*attrs['nLat']*attrs['nAlt']
    for val in var:  # loops through all variable data stored in file
        # Trim variable names.
        v = sub('\[|\]', '', val.decode('utf-8')).strip()
        s = unpack(endChar+'l', f.read(4))[0]
        # Test to see if this variable is desired
        gvar = True
        if len(varlist) > 0:
            try:
                varlist.index(v)
            except ValueError:  # include coordinate data or data from new file
                if ((v != 'Altitude' and v != 'Longitude' and
                     v != 'Latitude') or not newfile):
                    gvar = False
        # Unpack the data and save, if desired
        temp = unpack(endChar+'%id' % (nTotal), f.read(s))
        if gvar:
            # if a coordinate, only take unique values to speed up calculations
            if v == 'Longitude':
                coords['lon'] = np.unique(np.array(temp))*180.0/np.pi
            elif v == 'Latitude':
                coords['lat'] = np.unique(np.array(temp))*180.0/np.pi
            elif v == 'Altitude':
                # km, need for height integration
                coords['height'] = np.unique(np.array(temp))/1000.
            else:  # put in height, lat, lon order, noting fortran ordering
                variables[v] = np.transpose(np.array(temp).reshape(
                    (attrs['nLon'], attrs['nLat'], attrs['nAlt']),
                    order='F'), [2, 1, 0])
        f.read(4)
    return coords, variables, attrs


def append_units(variables):
    '''
    Append units, descriptive names, and plot scaling (e.g. linear,
    exponential) to the attributes of known data types
    '''

    var_dict = {}  # to store names, units, etc
    for k in list(variables.keys()):
        if type(variables[k]) is np.ndarray:
            nk = k
            # Temporary fix for misspelled key (9/30/13)
            if nk.find("Heaing Efficiency") >= 0:
                nk = "Heating Efficiency"
            elif nk.find("EUV Heating") >= 0:
                nk = "EuvHeating"
            elif nk.find("Rho (kg/m3)") >= 0:
                nk = "Rho"
            elif nk.find("Neutral Temperature (K)") >= 0:
                nk = 'Temperature'
            elif nk.find("Vn (up) (m/s)") >= 0:
                nk = 'V!Dn!N (up)'
            elif nk.find("Vi (east) (m/s)") >= 0:
                nk = "V!Di!N (east)"
            # DDZ added 2020-08-31
            if nk == "Magnetic latitude":
                nk = "Magnetic Latitude"

            try:
                # print ('registering',k)
                var_dict = register_name(k, nk, unit_dict, scale_dict,
                                         name_dict, var_dict)
            except:
                if k.split()[0] in name_dict:
                    nk = k.split()[0]
                    if "(/m3)" in k:
                        # print 'found (/m3) in key'
                        unit_dict[k] = "m^{-3}"
                        name_dict[k] = name_dict[nk]
                        scale_dict[k] = scale_dict[nk]
                        var_dict = register_name(k, k, unit_dict, scale_dict,
                                                 name_dict, var_dict)
                    else:
                        raise
                else:
                    raise
    return var_dict


def register_name(k, nk, unit_dict, scale_dict, name_dict, var_dict):
    # Different versions of GITM differ in header capitalization
    if nk not in name_dict:
        # Try to capitalize or lowercase the key
        if nk == nk.capitalize():
            nk = k.lower()
        else:
            nk = k.capitalize()
    var_dict[k] = [name_dict[nk], scale_dict[nk], unit_dict[nk]]
    return var_dict


def calc_tec(coords, attrs, variables, verbose=False):
    '''
    A routine to calculate the 2D VTEC.
    To perform these calculations, electron density ("e-") must be one of
    the available data types.
    '''
    import scipy.integrate as integ

    if 'e-' in variables.keys() or 'e-                   (/m3)' in\
            variables.keys():
        if 'e-' in variables.keys():
            e_key = 'e-'
        elif 'e-                   (/m3)' in variables.keys():
            e_key = 'e-                   (/m3)'
        temp = np.array(variables[e_key] * 1.0e-16)
        variables['VTEC'] = np.zeros((attrs['nLat'], attrs['nLon']),
                                     dtype=np.float64)
        for ilon in range(attrs['nLon']):
            for ilat in range(attrs['nLat']):
                # Integrate electron density over altitude, not including
                # ghost cells
                vtec = integ.simps(temp[2:-2, ilat, ilon],
                                   coords['height'][2:-2]*1000, "avg")
                variables['VTEC'][ilat, ilon] = vtec
    elif verbose:
        print('e-!', variables.keys())
    return variables


def calc_2dion(coords, attrs, variables, verbose=False):
    '''
    A routine to calculate the 2D ionospheric parameters: hmF2, NmF2.
    To perform these calculations, electron density ("e-") must be one of
    the available data types.
    '''
    from scipy.interpolate import interp1d
    from scipy.signal import argrelextrema

    if 'e-' in variables.keys() or 'e-                   (/m3)' in\
            variables.keys():
        if 'e-' in variables.keys():
            e_key = 'e-'
        elif 'e-                   (/m3)' in variables.keys():
            e_key = 'e-                   (/m3)'
        variables['NmF2'] = np.zeros((attrs['nLat'], attrs['nLon']),
                                     dtype=np.float64)
        variables['hmF2'] = np.zeros((attrs['nLat'], attrs['nLon']),
                                     dtype=np.float64)
        alt = np.linspace(min(coords['height'][2:-2]),
                          max(coords['height'][2:-2]), 1000)
        ilon, ilat = coords['lon'], coords['lat']
        for ilon in range(attrs['nLon']):
            for ilat in range(attrs['nLat']):

                # Interpolate over the electron density altitude profile
                eprof = interp1d(coords['height'][2:-2],
                                 variables[e_key][2:-2, ilat, ilon],
                                 kind="cubic")
                edens = eprof(alt)

                # Find the local maxima of the electron density profile
                emax = argrelextrema(edens, np.greater)
                emax_list = list(emax[0])
                saddle = False
                if len(emax_list) == 0:
                    # No local maxima were found, try identifying
                    # saddle points
                    saddle = True
                elif len(emax_list) == 1:   # emax_list is a number
                    if max(edens) == edens[emax_list[0]]:
                        # Only one maxima exists and is realistic. Sometimes
                        # a profile will have inflection points instead of
                        # local maxima and this can confuse the routine
                        variables['NmF2'][ilat, ilon] = edens[emax_list[0]]
                        variables['hmF2'][ilat, ilon] = alt[emax_list[0]]
                    else:
                        saddle = True
                elif alt[emax_list[-1]] < 120.0:
                    saddle = True
                else:
                    # More than one maxima exists.  Seperate hmF2 from hmF1
                    # and spurious local maxima
                    NmF2 = list(edens[emax_list])
                    HmF2 = list(alt[emax_list])

                    # If the global maximum is over 200 km,
                    # this is the F2 peak
                    eindex = NmF2.index(max(NmF2))
                    if HmF2[eindex] <= 200.0 and HmF2[eindex] == min(HmF2):
                        # The global max may be the F1 peak, see if the
                        # secondary (or lessor) maxima is at the upper
                        # limit of the model.  If so, remove this point
                        # from consideration
                        if max(HmF2) > coords['height'][-5]:
                            eindex = HmF2.index(max(HmF2))
                            emax_list.pop(eindex)
                            NmF2.pop(eindex)
                            HmF2.pop(eindex)
                            eindex = NmF2.index(max(NmF2))

                        if len(emax_list) > 1:
                            # If there are multiple maxima after the upper
                            # boundary has been removed, choose the largest
                            # maxima above 200 km since the hmF1 is often
                            # larger than the hmF2
                            emax_list.pop(eindex)
                            NmF2.pop(eindex)
                            eindex = NmF2.index(max(NmF2))

                    # Set the hmF2 and NmF2 (#emax_list is a number)
                    variables['NmF2'][ilat, ilon] = edens[emax_list[eindex]]
                    variables['hmF2'][ilat, ilon] = alt[emax_list[eindex]]

                if saddle:
                    # It is difficult to find saddle points.  Examine
                    # the rate of change of density
                    delta_alt = alt[1] - alt[0]  # Equally spaced
                    edens_dot = np.diff(edens) / delta_alt

                    # Identify inflection points by looking for
                    # minima in the derivative of electron density.
                    edot_min = argrelextrema(edens_dot, np.less)
                    emin_list = list(edot_min[0])
                    edens_min = list(edens[emin_list])
                    emax = np.max(edens_min)
                    eindex = emin_list[edens_min.index(emax)]

                    # Assign the inflection with the largest
                    # electron density to the ion peak
                    # emax_list is a number
                    variables['NmF2'][ilat, ilon] = edens[eindex]
                    variables['hmF2'][ilat, ilon] = alt[eindex]
    elif verbose:
        print('e- not found', variables.keys())
    return variables


def GitmBin(filename, flag_2D, verbose=False):
    '''
    Object to open, manipulate and visualize 1-3 dimensional GITM output
    stored in binary format.  Object inherits from spacepy.pybats.PbData; see
    that documentation for general information on how these objects work.

    GITM index ordering is [lon, lat, altitude]; data arrays read from file
    will always be of the same shape and size.

    kwargs may be specified for:
    magfile = 3DION or 3DMAG file, allows computation of velocities in magnetic
              coordinates
    varlist = list of variable keys.  Will limit the variables read in to those
              listed.  Time and position will always be read in.  If the list
              is empty, all variables will be read in.
    '''
    varlist = []
    attrs = {'file': filename}

    # Load the GITM data
    coords, variables, attrs = _read(varlist, attrs)
    if not flag_2D:
        variables = calc_tec(coords, attrs, variables, verbose=verbose)
        variables = calc_2dion(coords, attrs, variables, verbose=verbose)
    var_dict = append_units(variables)
    return filename, coords, variables, var_dict, attrs


def gitm_toCDF(bin_file, coords, variables, var_dict, attrs):
    '''Prepare and write data into netcdf files.'''

    # perform longitude shift, careful of extra values on ends
    # this logic is different because the lon grid doesn't start at zero
    #   or end at 360. In some cases it goes farther on both ends.
    lon_le180 = list(np.where((coords['lon'] <= 180.) &
                              (coords['lon'] >= 0.))[0])
    lon_ge180 = list(np.where((coords['lon'] >= 180.) &
                              (coords['lon'] < 360.))[0])
    if 180. not in coords['lon']:
        lon_ge180.insert(0, lon_le180[-1])
        lon_le180.append(lon_ge180[1])
    full_lon = np.append(coords['lon'][lon_ge180] - 360.,
                         coords['lon'][lon_le180])  # keep separate
    new_lon = np.append(coords['lon'][lon_ge180] - 360.,
                        coords['lon'][lon_le180])  # in memory
    new_lon[0], new_lon[-1] = -180., 180.  # begins at -180, ends at +180.
    coords['lon'] = new_lon  # use corrected grid for file

    # prepare for interpolation at +/- 180 deg lon, possible also latitude
    # interpolation is done after lon wrapping, so wrapped lon grid here
    ko = Kamodo()
    ko.variables = {'TMP': {}}
    ko._registered = 0  # need full lon and lat grids
    coord_dict = {'lon': {'units': 'deg', 'data': full_lon},
                  'lat': {'units': 'deg', 'data': coords['lat']}}

    # check for and prepare for latitude wrapping or interpolation
    lat_check = False
    if coords['lat'].max() < 90.:
        coords['lat'] = np.insert(coords['lat'], 0, -90.)
        coords['lat'] = np.append(coords['lat'], 90.)
        coord_dict['lat']['data'] = coords['lat']
        print('Performing scalar averaging at the poles.')
    elif coords['lat'].max() > 90.:
        old_lat = coords['lat']
        lat_idx = list(np.where((old_lat > -90.) & (old_lat < 90.))[0])
        lat_idx = [min(lat_idx) - 1] + lat_idx + [max(lat_idx) + 1]  # add ends
        coord_dict['lat']['data'] = old_lat[lat_idx]
        coords['lat'] = old_lat[lat_idx]  # corrected lat grid for files
        coords['lat'][0], coords['lat'][-1] = -90., 90.
        lat_check = True  # signal to perform interpolation at poles
        print('Interpolating data at the poles.')

    # start new output object
    cdf_filename = bin_file.replace('.bin', '.nc')
    data_out = Dataset(cdf_filename, 'w', format='NETCDF4')
    # csv list of files
    data_out.file = bin_file
    data_out.model = 'GITM'
    data_out.version = attrs['version']
    data_out.endian = attrs['endian']
    data_out.filedate = attrs['filedate']
    for dim in coords.keys():  # register dimensions
        if coords[dim].size == 1 and dim != 'time':
            continue  # skip height for 2D variables
        # create dimension
        new_dim = data_out.createDimension(dim, coords[dim].size)
        # create variable
        new_var = data_out.createVariable(dim, np.float32, dim)
        new_var[:] = coords[dim]  # store data for dimension in variable
        if dim == 'height':
            units = 'km'
        if dim == 'time':
            units = 'hr'
        else:
            units = 'deg'
        new_var.datascale, new_var.units = 'linear', units

    # copy over variables to file
    # var_3D = ['TEC', 'NmF2', 'hmF2', 'SolarLocalTime', 'SolarZenithAngle',
    #          'phi_qJoule', 'phi_q', 'phi_qEUV', 'phi_qNOCooling']
    for variable_name in variables.keys():
        new_name = gitm_varnames[var_dict[variable_name][0]][0]
        # remove height dimension for 2D data
        if 1 in variables[variable_name].shape:
            variables[variable_name] = np.squeeze(variables[variable_name])
        if len(variables[variable_name].shape) == 3:
            new_var = data_out.createVariable(new_name, np.float32,
                                              ('lon', 'lat', 'height'))
            if 'height' not in coord_dict.keys():
                coord_dict['height'] = {'units': 'km',
                                        'data': coords['height']}
        elif len(variables[variable_name].shape) == 2:
            new_var = data_out.createVariable(new_name, np.float32,
                                              ('lon', 'lat'))
            if 'height' in coord_dict.keys():
                del coord_dict['height']
        # 3D: (h,lat,lon) -> (lon,lat,h), 2D: (lat,lon) -> (lon,lat)
        tmp = variables[variable_name].T

        # check for and perform latitude wrapping or prepare for interpolation
        if lat_check:
            tmp = tmp[:, lat_idx]  # chop off extra latitude values
        elif tmp.shape[1] < len(coords['lat']):
            spole = np.mean(tmp[:, 0], axis=0)  # avg over lons near pole
            npole = np.mean(tmp[:, -1], axis=0)  # other pole
            tmp = np.insert(tmp, tmp.shape[1]-1, npole, axis=1)  # append
            tmp = np.insert(tmp, 0, spole, axis=1)  # insert on lat axis
        new_data = tmp[lon_ge180 + lon_le180]  # wrap in longitude

        # replace lon borders with proper interpolated boundary values
        ko = Functionalize_Dataset(ko, coord_dict, 'TMP',
                                   {'data': new_data, 'units': ''}, True,
                                   'GDZsph')
        new_data[0] = ko['TMP_ijk'](lon=-180.)
        new_data[-1] = new_data[0]  # wrap in longitude with new values

        # check for and perform latitude interpolation
        if lat_check:
            new_data[:, 0] = ko['TMP_ijk'](lat=-90.)
            new_data[:, -1] = ko['TMP_ijk'](lat=90.)
        del ko['TMP_ijk'], ko['TMP']

        # store new variable data in file
        new_var[:] = new_data
        new_var.datascale, new_var.units =\
            gitm_varnames[var_dict[variable_name][0]][1:]

    # close file
    data_out.close()
    return cdf_filename


@np.vectorize
def dts_to_hrs(datetime_string, filedate):
    '''Get hours since midnight from datetime string'''

    return (datetime.strptime(datetime_string, '%Y-%m-%d %H:%M:%S'
                              ).replace(tzinfo=timezone.utc) -
            filedate).total_seconds()/3600.


def GITMbin_toCDF(file_dir, flag_2D=False):
    '''Collect data from all files found with file_prefix into netCDF4 files
    '''
    # If flag_2D is False, then 2D files are not present and
    # 2D variables will be calculated
    # Takes much longer if False.

    from os.path import basename, isfile

    ftic = perf_counter()
    files = sorted(glob(file_dir+'*.bin'))  # collect files
    if len(files) == 0:
        print('No unconverted files found.')
        return
    print('Checking for files that need converting...', end="")

    # initialize date values from first file
    file_name = basename(files[0])
    file_dt_str = file_name.split('_t')[-1].split('_')[0]  # keep only YYMMDD
    if int(file_dt_str[:2]) > 60:  # later than 1960
        string_date = '19' + file_dt_str[:2] + '-' + file_dt_str[2:4] + '-' +\
            file_dt_str[4:6]
    else:
        string_date = '20' + file_dt_str[:2] + '-' + file_dt_str[2:4] + '-' +\
            file_dt_str[4:6]
    filedate = datetime.strptime(string_date+' 00:00:00', '%Y-%m-%d %H:%M:%S'
                                 ).replace(tzinfo=timezone.utc)  # dt object

    # loop through all .bin files
    file_count = 0
    for file in files:
        file_tic = perf_counter()
        cdf_file = file.replace('.bin', '.nc')
        if isfile(cdf_file):
            continue  # do not convert if conversion already completed

        # read data and attributes from binary file
        file_count += 1
        filename, coords, variables, var_dict, attrs = GitmBin(file, flag_2D)
        time = np.array([attrs['time'][:19]])  # accurate to the sec
        coords['time'] = dts_to_hrs(time, filedate)
        attrs['filedate'] = filedate.strftime('%Y-%m-%d')  # store in attrs

        # perform data wrangling and write to cdf file
        cdf_file = gitm_toCDF(file, coords, variables, var_dict, attrs)
        print(f'Converted {file} to {cdf_file} in ' +
              f'{perf_counter()-file_tic:.5f}s')
    if file_count > 0:
        print(f'Files of pattern {file_dir+"*.bin"} converted to ' +
              f'cdf in {perf_counter()-ftic:.5f}s')
    else:
        print('All files are already converted.')
    return True
