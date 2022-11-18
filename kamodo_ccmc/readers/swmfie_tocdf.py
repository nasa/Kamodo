# -*- coding: utf-8 -*-
"""
Created on Mon May 17 18:53:35 2021
@author: rringuet
"""
from glob import glob
import numpy as np
from time import perf_counter
from os.path import basename, isfile
from netCDF4 import Dataset

swmfie_varnames = {"X": ['x', 'km'], "Y": ['y', 'km'], "Z": ['z', 'km'],
                   # coordinates are ignored
                   "Theta": ['theta', "deg"], "Psi": ['psi', "deg"],
                   "Btilt_theta": ['theta_Btilt', "deg"],
                   "Btilt_psi": ['psi_Btilt', "deg"],
                   "SigmaH": ['Sigma_H', "S"], "SigmaP": ['Sigma_P', "S"],
                   "E-Flux": ['Phi_E', "W/m**2"], "Ave-E": ['E_avg', 'eV'],
                   "JR": ["j_R", "mA/m**2"], "PHI": ["Phi", "kV"],
                   "Ex": ["E_x", "mV/m"], "Ey": ["E_y", "mV/m"],
                   "Ez": ["E_z", "mV/m"],
                   "Jx": ["j_x", "mA/m**2"], "Jy": ["j_y", "mA/m**2"],
                   "Jz": ["j_z", "mA/m**2"],
                   "Ux": ['v_x', "km/s"], "Uy": ['v_y', "km/s"],
                   "Uz": ['v_z', "km/s"],
                   "JouleHeat": ['Q_Joule', "mW/m**2"],
                   "IonNumFlux": ['Phi_nion', "1/cm**2/s"],
                   "RT 1/B": ['Binv_RT', "1/T"],
                   "RT Rho": ['rho_RT', "amu/cm**3"],
                   "RT P": ['P_RT', "Pa"],
                   "conjugate dLat": ['dLat_star', "deg"],
                   "conjugate dLon": ['dlon_star', "deg"]}

'''
Documentation variables:
(ignored, given in R_E on a unit sphere)
"X":['x','km'],"Y":['y','km'],"Z":['z','km'],
"Theta":['theta',"deg"],"Psi":['psi',"deg"],     (used as coordinates)
"Btilt_theta":['theta_Btilt',"deg"], "Btilt_psi":['psi_Btilt',"deg"]
(added directly to object for documentation purposes)
'''


def convert_all(file_dir):
    '''Convert all files independent of the others.'''

    print('NetCDF version of data not found. Converting tec files in ' +
          f'{file_dir} to netCDF.')
    ftic, nfiles = perf_counter(), 0
    files = sorted(glob(file_dir + '*.tec'))  # want YYYYMMDD
    if len(files) == 0:
        print('No original files found.')
        return

    # convert and return
    for file in files:
        if not isfile(file.replace('.tec', '.nc')):
            coords, variables, var_units, theta_Btilt, psi_Btilt = \
                _read_SWMFIE(file, verbose=True)
            cdf_filename = _toCDF(file, coords, variables, var_units,
                                  theta_Btilt, psi_Btilt)
            nfiles += 1
            print(cdf_filename, ' converted.')
    print(f'{nfiles} files in {file_dir} now converted into ' +
          f'netCDF4 files in {perf_counter()-ftic:.6f}s.')
    return None


def read_SWMFIE_header(filename):
    '''parse header for a representative file'''

    # get file data
    file_object = open(filename, 'r')
    file_contents = file_object.read()
    file_object.close()

    # sort info with TITLE keyword
    title_lines = file_contents.split('TITLE=')[1].split('VARIABLES=')[0].\
        replace('\n', ',').replace('"', '').strip().strip(',')
    title, time, BtiltDeg_str = title_lines.split(',')
    # get time into 'YYYY-MM-DD HH:mm:SS' format, clean up BtiltDeg
    dts = time[:10]+' '+''.join([time[11:].split('-')[i]+':'
                                 for i in [0, 1]])+time[11:].split('-')[2]
    BtiltDeg = BtiltDeg_str.split('=')[1].split()

    # sort info with VARIABLE keyword
    variable_lines = file_contents.split('VARIABLES=')[1].split('ZONE')[0].\
        replace('\n', ',').replace('"', '').strip().strip(',')
    variable_list = [var.strip().split(']')[0].split(' [')[0] for var in
                     variable_lines.split(',')]  # units in dict above

    # sort info with ZONE keyword 'IonN N=0000241 T=0000:20:00'
    zone_lines = file_contents.split('ZONE T=')[1].split('\n')[0].strip(
        '"').strip()
    run_type1, N1, t = zone_lines.split(' ')
    N1 = int(N1.split('=')[1].strip())
    zone_lines = file_contents.split('ZONE T=')[2].split('\n')[0].strip(
        '"').strip()
    run_type2, N2, t = zone_lines.split(' ')
    N2 = int(N2.split('=')[1].strip())

    # sort coordinate info   'I=           91  J=          181  F=POINT'
    coord_info = file_contents.split(zone_lines)[1].split('\n')[1].strip()
    text, i, text, j, f = coord_info.split()
    i, j, f = int(i), int(j), f.split('=')[1]
    # if sometimes different, need another set for this

    # determine number of lines to skip for data sections
    line_test = [f in c for c in file_contents.split('\n')]
    skip1, skip2 = np.where(np.array(line_test, dtype=int) == 1)[0]

    # determine how many lines constitute one data list
    data_lines = file_contents.split('\n')[skip1+1:]
    num_test = [len(line) for line in data_lines]  # list of lengths of lines
    ndata_lines = np.diff(np.where(np.array(num_test) == num_test[0])[0])[0]

    # cleanup and return
    del file_contents
    header = {'title': title, 'time': dts, 'BtiltDeg': BtiltDeg,
              'variable_list': variable_list, 'run_type1': run_type1, 'N1': N1,
              'run_type2': run_type2, 'N2': N2, 'i': i, 'j': j, 'f': f,
              'skip1': skip1, 'skip2': skip2, 'ndata_lines': ndata_lines}
    return header


def read_SWMFIE_data(filename, header):
    '''read only data from file and return in labeled dictionary'''

    # get data from file using skip as number of lines to skip
    file_object = open(filename, 'r')
    file_contents = file_object.readlines()
    file_object.close()

    # Get theta_Btilt value from file
    title_line = file_contents[0].replace('\n',
                                          ',').replace('"',
                                                       '').strip().strip(',')
    theta_Btilt = float(title_line.split(',')[-1].strip().split()[1])

    # collect data into an array, one for each run_type
    file_data1 = file_contents[header['skip1']+1:header['skip2']-1]
    data1 = np.array([''.join(file_data1[i:i+header['ndata_lines']]).
                      strip('\n').split() for i in range(
                          0, len(file_data1), header['ndata_lines'])],
                     dtype=float)
    file_data2 = file_contents[header['skip2']+1:]
    data2 = np.array([''.join(file_data2[i:i+header['ndata_lines']]).
                      strip('\n').split() for i in range(
                          0, len(file_data2), header['ndata_lines'])],
                     dtype=float)

    # determine correct dimension sizes from theta and psi (3 and 4)
    # combine data and return in dict
    nLatA, nLon = np.unique(data1[:, 3]).size, np.unique(data1[:, 4]).size
    data = np.concatenate((np.reshape(data1, (nLon, nLatA,
                                              len(header['variable_list']))),
                           np.reshape(data2, (nLon, nLatA,
                                              len(header['variable_list'])))),
                          axis=1)[:, ::-1, :]
    # leave since grids are identical in N/S hemispheres

    return {swmfie_varnames[header['variable_list'][i]][0]:
            np.transpose(data[:, :, i], [1, 0])
            for i in range(len(header['variable_list']))}, theta_Btilt


def _read_SWMFIE(file, verbose=False):
    '''file_dir should be the complete path to the location of the files.'''
    t0 = perf_counter()

    # establish time attributes first for file searching
    file_datestr = basename(file)[3:11]
    # string_date = 'YYYY-MM-DD'
    string_date = file_datestr[:4]+'-'+file_datestr[4:6]+'-'+file_datestr[6:8]

    # get list of variables possible in these files using first file
    file_header = read_SWMFIE_header(file)
    file_varlist = file_header['variable_list']
    gvar_list = [swmfie_varnames[key][0] for key in file_varlist if key not in
                 ['X', 'Y', 'Z', 'Theta', 'Psi', 'Btilt_theta', 'Btilt_psi']]
    # avoid returning coordinates stored elsewhere (or ignored)

    # collect coordinates and intialize data storage from first file
    data, Btilt = read_SWMFIE_data(file, file_header)
    lon, lat0 = np.unique(data['psi']), np.unique(data['theta'])-90.
    lon_ge180 = list(np.where((lon >= 180.) & (lon < 360.))[0])
    lon_le180 = list(np.where(lon <= 180.)[0])
    lon_idx = lon_ge180 + lon_le180
    lon -= 180  # shifting longitude to be in range -180 to 180
    lat0 = np.insert(lat0, np.where(lat0 == 0.)[0], 0.)
    # interpolator requires unique ascending values
    # offset two zero values by 0.0001
    lat0[np.where(lat0 == 0.)[0]] = -0.0000001, 0.0000001
    lat = np.insert(lat0, 0, -90.)
    lat = np.append(lat, 90.)  # add values for poles for wrapping
    # data is 2D (lon, lat) with time, so no height
    theta_Btilt = Btilt   # to append more values later
    psi_Btilt = 0.   # this value is always zero

    # Store variable data and units for the file.
    # transpose (lat, lon) -> (lon, lat) and wrap longitude
    var_units = {value[0]: value[-1] for key, value in swmfie_varnames.items()
                 if value[0] in gvar_list}
    variables = {key: np.array(data[key]).T[lon_idx] for key in gvar_list}

    # wrap latitude for scalar variables (all to be treated as scalars)
    for var_key in variables:
        tmp = variables[var_key]
        new_shape = list(tmp.shape)
        new_shape[1] += 2  # one value on each end
        tmp2 = np.zeros(new_shape)
        tmp2[:, 1:-1] = tmp
        # put in top values
        top = np.mean(tmp2[:, 1], axis=0)  # single value
        new_top = np.broadcast_to(top, (new_shape[0]))
        tmp2[:, 0] = new_top
        # same for bottom, reusing variable names
        top = np.mean(tmp2[:, -2], axis=0)  # single value
        new_top = np.broadcast_to(top, (new_shape[0]))
        tmp2[:, -1] = new_top
        variables[var_key] = tmp2  # store result

    coords = {'lat': lat, 'lon': lon}
    return coords, variables, var_units, theta_Btilt, psi_Btilt
# removed files from return statement


def _toCDF(file, coords, variables, var_units, theta_Btilt, psi_Btilt):
    '''Write data to a netCDF4 file.'''

    # start new wrapped output object
    cdf_filename = file.replace('.tec', '.nc')
    data_out = Dataset(cdf_filename, 'w', format='NETCDF4')
    data_out.file = file
    data_out.model = 'SWMF_IE'
    for dim in coords.keys():
        units = 'deg'
        # create dimension
        new_dim = data_out.createDimension(dim, coords[dim].size)
        # create variable
        new_var = data_out.createVariable(dim, np.float64, dim)
        new_var[:] = coords[dim]  # store data for dimension in variable
        new_var.units = units

    # copy over variables to file
    for variable_name in variables.keys():
        if len(variables[variable_name].shape) == 2:
            new_var = data_out.createVariable(variable_name, np.float64,
                                              ('lon', 'lat'))
            new_data = variables[variable_name]
        else:
            continue
        new_var[:] = new_data  # store data in variable
        new_var.units = var_units[variable_name]
    data_out.close()
    return cdf_filename
