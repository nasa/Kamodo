# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 10:21:46 2023

@author: rringuet

"""
from os.path import basename
from numpy import array, zeros, float32, mean
import h5py

# Is this even necessary? Does the interpolator require the entire dataset
# for a variable to be assembled, or does it load sections lazily?
# This is needed to calculate cell centers and load the coord grids.
# Need the entire coordinate grid to understand where the data pieces are.
# The requested variables need to have the same grid size as each other.
# Cannot ask for both X and dV because their sizes are different by 1 in each
# block.
def GridMapping(data_files, constant, sample_var):
    '''Determine what grid sections are in what files. Returns index mapping
    dictionary of beginning and ending indices for each block. Note that the
    constant variables are one size smaller than the coordinates because they
    are on cell-centers while the coordinates are the cell edges. Otherwise,
    the variables all have the same size, so the same index mapping may be used
    for all of them. This is called to map the coordinate variables together
    for calculation of the cell centers in the pre-processing step.
    Inputs:
        data_files - a list of the h5 files containing the data
        constant - boolean, whether the variable is constant or not. This
            decides whether the program looks in the 'Step#0' group or the main
            group for the data.
        sample_var - the name of a variable with the same grid (string
    Example usage of returned dictionary:
        net_array = np.zeros(net_indices['all'][1])
        net_array[net_indices[key][0][0]:net_indices[key][1][0]] = data_block
        for a 1D array, where data_block is from the file corresponding to the
        key name, and key is the block number of format '00xi_00yi_00zi'.
    '''
    # first, get file grid from filename of first file
    # pattern: abcd_00nx_00ny_00nz_0nxi_0nyi_0nzi....
    nx, ny, nz = array(basename(data_files[0]).split('_')[1:4],
                       dtype=int)
    # loop through files and get size of each grid
    indices = {}
    for f in data_files:
        key = ''.join([b+'_' for b in basename(f).split('.')[0].split('_')[-3:]
                       ])[:-1]  # 0000_0000_0000 for first block, etc
        h5_data = h5py.File(f)
        if constant:
            indices[key] = h5_data[sample_var].shape
        else:
            indices[key] = h5_data['Step#0'][sample_var].shape
        h5_data.close()
        
    # then get the total shape by summing the data shapes
    # this assumes that the assembled block will have a single (x, y, z) shape
    x_shape = sum([indices[str(i).zfill(4)+'_0000_0000'][0] for i in
                   range(nx)])
    y_shape = sum([indices['0000_'+str(i).zfill(4)+'_0000'][1] for i in
                   range(ny)])
    z_shape = sum([indices['0000_0000_'+str(i).zfill(4)][2] for i in
                   range(nz)])

    # determine beginning and ending indices for each block
    net_indices = {'net': [[0, 0, 0], [x_shape, y_shape, z_shape]]}
    for key in indices.keys():
        # get beginning position mapping in array
        nxi, nyi, nzi = array(key.split('_'), dtype=int)
        xi = sum([indices[str(i).zfill(4)+'_0000_0000'][0] for i in
                  range(nxi)])
        yi = sum([indices['0000_'+str(i).zfill(4)+'_0000'][1] for i in
                  range(nyi)])
        zi = sum([indices['0000_0000_'+str(i).zfill(4)][2] for i in
                  range(nzi)])
        # get ending position+1
        xip1, yip1, zip1 = array([xi, yi, zi]) + indices[key]
        net_indices[key] = [[xi, yi, zi], [xip1, yip1, zip1]]
    print(net_indices)
    return net_indices


def AssembleGrid(data_files, net_indices, constant, var, step=0):
    '''Assemble the requested variable into a single grid from all the files.
    Used for the coordinate assembly. This assumes that there are no ghost
    cells.

    Inputs:
        data_files - a list of h5 files containing the data to be assembled
        net_indices - the output from the GridMapping function above
        constant - boolean, whether the variable is constant or not. This
            decides whether the program looks in the 'Step#?' group or the main
            group for the data.
        var - a string indicating the name of the desired variable in the files
        step - an integer indicating the time step of the data to be assembled,
            default is zero.
    Output:
        A float32 array of the assembled data.
    '''

    data = zeros(net_indices['net'][1], dtype=float32)
    for f in data_files:
        key = ''.join([b+'_' for b in basename(f).split('.')[0].split('_')[-3:]
                       ])[:-1]  # 0000_0000_0000 for first block, etc
        h5_data = h5py.File(f)
        x0, y0, z0 = net_indices[key][0]
        x1, y1, z1 = net_indices[key][1]
        if constant:
            data[x0:x1, y0:y1, z0:z1] = array(h5_data[var])
        else:
            data[x0:x1, y0:y1, z0:z1] = array(h5_data['Step#'+str(step)][var])
        h5_data.close()
    return data


def ComputeCellCenters(data_files, write_file=True):
    '''Computes the centers of the X, Y, and Z coordinate grids.
    Input:
        data_files - list of h5 files containing the data
        write_file - boolean (default=True), writes cells centers to a netCDF4
            file if True, doesn't if False.
    Outputs:
        X_c - numpy array of X locations of cell centers, one size smaller than
            the assembled X grid from the file(s)
        Y_c - numpy array of Y locations of cell centers, one size smaller than
            the assembled Y grid from the file(s)
        Z_c - numpy array of Z locations of cell centers, one size smaller than
            the assembled Z grid from the file(s)
    '''
    # assemble coordinate grid from the file(s)
    if len(data_files) > 1:
        # get index mapping for X, identical to that for Y and Z
        net_indices = GridMapping(data_files, True, 'X')
        # assemble coordinate grids from MPI files
        X = AssembleGrid(data_files, net_indices, True, 'X')
        Y = AssembleGrid(data_files, net_indices, True, 'Y')
        Z = AssembleGrid(data_files, net_indices, True, 'Z')
    else:
        # read in data from serial file
        h5_data = h5py.File(data_files[0])
        X = array(h5_data['X'])
        Y = array(h5_data['Y'])
        Z = array(h5_data['Z'])
        h5_data.close()
        
    # compute cell centers
    X_c = zeros([X.shape[0]-1, X.shape[1]-1, X.shape[2]-1])
    Y_c, Z_c = X_c.copy(), X_c.copy()
    for i in range(X_c.shape[0]):
        for j in range(X_c.shape[1]):
            for k in range(X_c.shape[2]):
                X_c[i,j,k] = mean(X[i:i+2, j:j+2, k:k+2])
                Y_c[i,j,k] = mean(Y[i:i+2, j:j+2, k:k+2])
                Z_c[i,j,k] = mean(Z[i:i+2, j:j+2, k:k+2])

    # write to file if requested
    if write_file:
        from netCDF4 import Dataset

        file_dir = data_files[0].split(basename(data_files[0]))[0]
        center_file = file_dir + 'gridcenters.nc'
        data_out = Dataset(center_file, 'w', format='NETCDF4')
        data_out.Description = 'Cell centers of the coordinate grids.'
        # create dimensions, all coord grids have the same shapes
        new_dim = data_out.createDimension('c1', X_c.shape[0])
        new_dim = data_out.createDimension('c2', X_c.shape[1])
        new_dim = data_out.createDimension('c3', X_c.shape[2])
        # save X center positions
        new_var = data_out.createVariable('X_center', float32,
                                          ('c1', 'c2', 'c3'))  
        new_var[:] = X_c
        new_var.units = 'R_E'
        new_var.Description = 'X coordinate of grid cell centers'
        # save Y center positions
        new_var = data_out.createVariable('Y_center', float32,
                                          ('c1', 'c2', 'c3'))  
        new_var[:] = Y_c
        new_var.units = 'R_E'
        new_var.Description = 'Y coordinate of grid cell centers'
        # save Z center positions
        new_var = data_out.createVariable('Z_center', float32,
                                          ('c1', 'c2', 'c3'))  
        new_var[:] = Z_c
        new_var.units = 'R_E'
        new_var.Description = 'Z coordinate of grid cell centers'
        # close file
        data_out.close()   
    
    return X_c, Y_c, Z_c
