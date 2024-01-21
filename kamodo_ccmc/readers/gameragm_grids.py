# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 10:21:46 2023

@author: rringuet

Modification History:
  2023/09/06 Lutz Rastaetter - initial code to defing 2D slices+azimuth grid
  2023/09/28 Lutz Rastaetter - set up triangular interpolation on 2D slices
  2023/11/17 Lutz Rastaetter - link to shared library for 2D+1D interpolation

"""
from os.path import basename, isfile
from numpy import array, zeros, float32, float64, mean, arctan2, sqrt, squeeze, flip, int32, where, prod, pi, reshape, all, any, NaN, floor, min, max, mod, pi, arccos
from kamodo_ccmc.readers.reader_utilities import _isfile,Dataset
import h5py                 # should use RU's S3 wrapper
#from kamodo_ccmc.readers.reader_utilities import h5py
from netCDF4 import Dataset # should use RU's S3 wrapper
#from kamodo_ccmc.readers.reader_utilities import Dataset # should use RU's S3 wrapper

def GridRanges(data_files):
    '''Finds the max and min of each coordinate grid across the entire dataset.
    Writes the result to a simple text file named 'GAMERA_GridRanges.txt' in
    the same directory as the dataset. Called by the block of logic in the
    gamera reader that creates the time files.

    Input (data_files) - a list of the data files in the dataset
    Output: None, but writes the text file described.'''

    # calculate the max and min of each coordinate across the dataset
    x_max, x_min, y_max, y_min, z_max, z_min = [], [], [], [], [], []
    coord_shape = []
    for f in data_files:
        h5_data = h5py.File(f)
        X = array(h5_data['X'])
        Y = array(h5_data['Y'])
        Z = array(h5_data['Z'])
        h5_data.close()
        x_min.append(X.min())
        x_max.append(X.max())
        y_min.append(Y.min())
        y_max.append(Y.max())
        z_min.append(Z.min())
        z_max.append(Z.max())
        coord_shape.append(list(X.shape))
    del X, Y, Z
    X_min, X_max = array(x_min).min(), array(x_max).max()
    Y_min, Y_max = array(y_min).min(), array(y_max).max()
    Z_min, Z_max = array(z_min).min(), array(z_max).max()

    # write results to a simple text file
    file_dir = data_files[0].split(basename(data_files[0]))[0]
    files_str = ''.join([basename(f)+',' for f in data_files])[:-1]
    out_file = file_dir + 'GAMERA_GridRanges.txt'
    out = open(out_file, 'w')
    out.write(f'GAMERA_GridRanges output for files: {files_str}\n')
    out.write(f'X_min = {X_min:.5f}\n')
    out.write(f'X_max = {X_max:.5f}\n')
    out.write(f'Y_min = {Y_min:.5f}\n')
    out.write(f'Y_max = {Y_max:.5f}\n')
    out.write(f'Z_min = {Z_min:.5f}\n')
    out.write(f'Z_max = {Z_max:.5f}\n')
    out.close()

    return None


def Read_GridRanges(file_dir):
    '''Reads in the previously calculated grid ranges from the file produced by
    the GridRanges function above.
    Input: file_dir (the complete file path where the data is stored)
    Outputs:
        X_min - float value indicating the minimum X coordinate value across
            the dataset in the given file directory.
        X_max - float value indicating the maximum X coordinate value across
            the dataset in the given file directory.
        Y_min - float value indicating the minimum Y coordinate value across
            the dataset in the given file directory.
        Y_max - float value indicating the maximum Y coordinate value across
            the dataset in the given file directory.
        Z_min - float value indicating the minimum Z coordinate value across
            the dataset in the given file directory.
        Z_max - float value indicating the maximum Z coordinate value across
            the dataset in the given file directory.
    '''

    range_file = file_dir + 'GAMERA_GridRanges.txt'
    range_data = open(range_file, 'r')
    for line in range_data.readlines():
        if line[0] == 'G':
            continue
        tmp = float(line.split('=')[1].strip())
        if line[:5] == 'X_min':
            X_min = tmp
        if line[:5] == 'X_max':
            X_max = tmp
        if line[:5] == 'Y_min':
            Y_min = tmp
        if line[:5] == 'Y_max':
            Y_max = tmp
        if line[:5] == 'Z_min':
            Z_min = tmp
        if line[:5] == 'Z_max':
            Z_max = tmp
    range_data.close()

    return X_min, X_max, Y_min, Y_max, Z_min, Z_max


# Is this even necessary? Does the interpolator require the entire dataset
# for a variable to be assembled, or does it load sections lazily?
# Eric Winter: The interpolator loads the sections as needed.

# This calculates the cell centers and loads the coordinate grids.
# Need the entire coordinate grid to understand where the data pieces are.
# The requested variables need to have the same grid size as each other.
# Cannot ask for both X and dV because their sizes are different by 1 in each
# block.

# actually the grid mapping is the same in terns of initial index with the cell corner (constant) data
# repeating at the borders, so the next block may overwrite the last layer of the preceding block.
# Only the tital number of positins in each direction is one larger

def GridMapping(data_files, constant, sample_var):
    '''Determine what grid sections are in what files. Returns index mapping
    dictionary of beginning and ending indices for each block. Note that the
    constant variables are one size smaller than the coordinates because they
    are on cell-centers while the coordinates are the cell edges. Otherwise,
    the variables all have the same size, so the same index mapping may be used
    for all of them. This is called to map the coordinate variables together
    for calculation of the cell centers in the pre-processing step. Avoids
    converting data to arrays, so execution time it quite fast.
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

        indices returned are [i,j,k] but array indexing in python may use [k,j,i]
        we may find that the time variable data are arranged in [k,i,j]
    '''
    # first, get file grid from filename of first file
    # pattern: abcd_00ni_00nj_00nk_0i_0j_0k....
    ni, nj, nk = array(basename(data_files[0]).split('_')[1:4],
                       dtype=int)
    di = 0
    if constant:
        di = 1
    # loop through files and get size of each grid
    indices = {}
    for f in data_files:
        key = ''.join([b+'_' for b in basename(f).split('.')[0].split('_')[-3:]
                       ])[:-1]  # 0000_0000_0000 for first block, etc
        h5_data = h5py.File(f)
        if constant:
            indices[key] = flip(h5_data[sample_var].shape)
        else:
            indices[key] = flip(h5_data['Step#0'][sample_var].shape)
        h5_data.close()
    # then get the total shape by summing the data shapes
    # this assumes that the assembled block will have a single (x, y, z) shape
    x_shape = sum([indices[str(i).zfill(4)+'_0000_0000'][0]-di for i in
                   range(ni)])+di
    y_shape = sum([indices['0000_'+str(i).zfill(4)+'_0000'][1]-di for i in
                   range(nj)])+di
    z_shape = sum([indices['0000_0000_'+str(i).zfill(4)][2]-di for i in
                   range(nk)])+di

    # determine beginning and ending indices for each block
    net_indices = {'net': [[0, 0, 0], [x_shape, y_shape, z_shape]]}
    for key in indices.keys():
        # get beginning position mapping in array
        nxi, nyi, nzi = array(key.split('_'), dtype=int)
        xi = sum([indices[str(i).zfill(4)+'_0000_0000'][0]-di for i in
                  range(nxi)])
        yi = sum([indices['0000_'+str(i).zfill(4)+'_0000'][1]-di for i in
                  range(nyi)])
        zi = sum([indices['0000_0000_'+str(i).zfill(4)][2]-di for i in
                  range(nzi)])
        # get ending position+1
        xip1, yip1, zip1 = array([xi, yi, zi]) + indices[key]
        net_indices[key] = [[xi, yi, zi], [xip1, yip1, zip1]]
    return net_indices


def AssembleGrid(data_files, net_indices, constant, var, step=0):
    '''Assemble the requested variable into a single grid from all the files.
    Used for the coordinate assembly.
    This assumes that there are no ghost cells in the GAMERA outputs.
    Constant-in-time variables (grid postions X,Y,Z) share surface positions with neighboring blocks.

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

    data = zeros(flip(net_indices['net'][1]), dtype=float32)
    for f in data_files:
        key = ''.join([b+'_' for b in basename(f).split('.')[0].split('_')[-3:]
                       ])[:-1]  # 0000_0000_0000 for first block, etc
        h5_data = h5py.File(f)
        x0, y0, z0 = net_indices[key][0]
        x1, y1, z1 = net_indices[key][1]
        if constant:
            data[z0:z1, y0:y1, x0:x1] = array(h5_data[var])
        else:
            data[z0:z1, y0:y1, x0:x1] = array(h5_data['Step#'+str(step)][var])
        h5_data.close()
    return data


def ComputeCellCentersFile(data_files):
    '''Computes the centers of the X, Y, and Z coordinate grids individually
    for each file. Avoids the complication of ghost cells during grid assembly.
    Input:
        data_files - list of h5 files containing the data
    Outputs: None, but writes a set of files ending in 'gridcenters.nc'
        containing the following variables for each file.
        X_c - numpy array of X locations of cell centers, one size smaller than
            the assembled X grid from the file(s)
        Y_c - numpy array of Y locations of cell centers, one size smaller than
            the assembled Y grid from the file(s)
        Z_c - numpy array of Z locations of cell centers, one size smaller than
            the assembled Z grid from the file(s)
    '''

    # read in data from each file
    file_dir = data_files[0].split(basename(data_files[0]))[0]
    for f in data_files:
        center_file = file_dir + basename(f).split('.')[0] + '_gridcenters.nc'
        if isfile(center_file):
            continue
        h5_data = h5py.File(f)
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
                    X_c[i, j, k] = mean(X[i:i+2, j:j+2, k:k+2])
                    Y_c[i, j, k] = mean(Y[i:i+2, j:j+2, k:k+2])
                    Z_c[i, j, k] = mean(Z[i:i+2, j:j+2, k:k+2])

        # prepare output file
        data_out = Dataset(center_file, 'w', format='NETCDF4')
        data_out.Description = 'Cell centers of the coordinate grids for ' + f

        # create dimensions, all coord grids have the same shapes,
        # but shapes could be diff in diff files
        new_dim = data_out.createDimension('c1', X_c.shape[0])
        new_dim = data_out.createDimension('c2', X_c.shape[1])
        new_dim = data_out.createDimension('c3', X_c.shape[2])
        # save X center positions
        new_var = data_out.createVariable('X_c', float32, ('c1', 'c2', 'c3'))
        new_var[:] = X_c
        new_var.units = 'R_E'
        new_var.Description = 'X coordinate of grid cell centers'
        # save Y center positions
        new_var = data_out.createVariable('Y_c', float32, ('c1', 'c2', 'c3'))
        new_var[:] = Y_c
        new_var.units = 'R_E'
        new_var.Description = 'Y coordinate of grid cell centers'
        # save Z center positions
        new_var = data_out.createVariable('Z_c', float32, ('c1', 'c2', 'c3'))
        new_var[:] = Z_c
        new_var.units = 'R_E'
        new_var.Description = 'Z coordinate of grid cell centers'
        # close file
        data_out.close()

    return None


def ComputeCellCentersTotal(data_files, write_file=True, grid=True):
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
    if len(data_files) > 1 and not grid:  # only do this if grid needs assembln
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
                X_c[i, j, k] = mean(X[i:i+2, j:j+2, k:k+2])
                Y_c[i, j, k] = mean(Y[i:i+2, j:j+2, k:k+2])
                Z_c[i, j, k] = mean(Z[i:i+2, j:j+2, k:k+2])

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


def PrepareCylindricalGrid(X,Y,Z):
    # the GAMERA (and LFM) grid is acutlly in a cylindrical arrangement
    # planes around the X-axis are formed by positions in (X,R=sqrt(Y^2+Z^2)) and azimuth=atan(Z,Y)
    # to cover the whole cylinder (ellipsoid) we need to have added a ghost cell layer in R (at the axis) and azimuth (both boundaries)

    nk,nj,ni = X.shape
    cylgrid_file= f'cylgrid_{ni:04d}_{nj:04d}_{nk:04d}.nc'
    
    if not _isfile(cylgrid_file+'foo'): # never true while developing
        X2D = squeeze(X[0,:,:])                              # along X-axis
        R2D = squeeze(sqrt(Y[0,:,:]*Y[0,:,:] + Z[0,:,:]*Z[0,:,:]))       # distance from X-axis
        R2D[0,:] = 0.
        # azimuth -- use positions farthest away from X-axis
        # any position except the first or last of second index should work
        PHI = mod(arctan2(squeeze(Z[:,int(nj/2),ni-1]),squeeze(Y[:,int(nj/2),ni-1]))+2*pi,2*pi) 
        PHI[-1] = PHI[0]+2*pi
        
        cylfile_out = Dataset(cylgrid_file,'w',format='NETCDF4')
        cylfile_out.Description = 'GAMERA grid converted to cylindrical in X, R (radial distance) and azumith PHI around X-axis'

        new_dim = cylfile_out.createDimension('c1', X2D.shape[0])
        new_dim = cylfile_out.createDimension('c2', X2D.shape[1])
        new_dim = cylfile_out.createDimension('c3', PHI.shape[0])

        X2D_out = cylfile_out.createVariable('X2D', float32, ('c1', 'c2'))
        R2D_out = cylfile_out.createVariable('R2D', float32, ('c1', 'c2'))
        PHI_out = cylfile_out.createVariable('PHI', float32, ('c3'))

        X2D_out[:] = X2D
        X2D_out.units='R_E'
        X2D_out.Description = 'grid positons: along X-axis'
        R2D_out[:] = R2D
        R2D_out.units='R_E'
        R2D_out.Description = 'grid postions: radial distance from X-axis'
        PHI_out[:] = PHI
        PHI_out.units='radian'
        PHI_out.Description = 'grid postions: azimuth angle around X-Axis'
        # close file
        cylfile_out.close()
    else:
        cylfile_ds = Dataset(cylgrid_file)
        X2D = array(cylfile_ds.variables['X2D'])
        R2D = array(cylfile_ds.variables['R2D'])
        PHI = array(cylfile_ds.variables['PHI'])
        cylfile_ds.close()
    return([X2D,R2D,PHI])


def PrepareCylindricalCellCenteredGrid(X2D,R2D,PHI):
    ''' Computes cell centered positions in X, R, PHI:
    adds a layer at the X-axis: zeros to R2D_c and copied values for X2D_C (at both ends of the second dimension) 
    a value is added to PHI at both ends to cover the full range [0,2*pi]
    '''
    n1,n2 = X2D.shape
    n3 = PHI.shape[0]
    n11 = n1+1
    n21 = n2+1
    n31 = n3+1
    X2D_c = zeros((n11,n2-1))
    R2D_c = zeros((n11,n2-1))
    PHI_c = zeros((n31))
    
    PHI_c[1:n3] = (PHI[1:n3]+PHI[0:n3-1])/2
    PHI_c[0] = PHI_c[n3-1]-2*pi
    PHI_c[n3] = PHI_c[1]+2*pi
    
    for i in range(X2D.shape[0]-1):
        for j in range(X2D.shape[1]-1):
            X2D_c[i+1, j] = mean(X2D[i:i+2, j:j+2])
            R2D_c[i+1, j] = mean(R2D[i:i+2, j:j+2])
    R2D_c[0,:] = 0.
    R2D_c[-1,:] = 0.
    X2D_c[0,:] = X2D_c[1,:]
    X2D_c[-1,:] = X2D_c[-2,:]

    return(X2D_c,R2D_c,PHI_c)

    
# position in 2D array to flat array index with -1 being out of range        
def i2(ij,n1,n2,perio1,perio2):  # position in 2D array to flat array index with -1 being out of range
    ''' Utility function that comoputes index in a (flattened) 2D array at positions  ij=(i,j) in a rectangle of n1 x n2 vertex positions.
        Used by PrepareTriangularInterpolation().
    '''
    ij_shape=ij[0].shape
    i_arr=ij[0].flatten()
    j_arr=ij[1].flatten()
    n1_1=n1-1
    n2_1=n2-1
    if perio1:
        w = where(i_arr[:] < 0)[0]
        for k in range(len(w)):
            i_arr[w[k]] = i_arr[w[k]] + n1
        w = where(i_arr[:] > n1_1)[0]
        for k in range(len(w)):
            i_arr[w[k]] = i_arr[w[k]] - n1
        
    if perio2:
        w = where(j_arr[:] < 0)[0]
        for k in range(len(w)):
            j_arr[w[k]] = j_arr[w[k]] + n2
        w = where(j_arr[:] > n2_1)[0]
        for k in range(len(w)):
            j_arr[w[k]] = j_arr[w[k]] - n2
    i2_tmp = i_arr[:] + n1 * j_arr[:]
    missing = [k for k in range(len(i_arr[:])) if (i_arr[k] < 0 or i_arr[k] > n1_1 or j_arr[k] < 0 or j_arr[k] > n2_1 ) ]
    for k in range(len(missing)): 
        i2_tmp[missing[k]] = -1
    return(reshape(i2_tmp,ij_shape))

#
# decompose logically carteian arrangement of grid cells into triangles
# this is independent of the spatial disctribution and size of triangles in the 2D space
#
def PrepareTriangularInterpolation(ni,nj):
    ''' given the dimensions ni,nj of the 2D grids in X,R (nj,ni)
    this function sets up arrays with the vertex indices, neighbor indices and i,j position for each triangle in the 2D grid.
    '''
    tri_file= f'tri_{ni:04d}_{nj:04d}.nc'    
    if not _isfile(tri_file+'_'):
        # create list of triangle vertices
        cell_tris = zeros((2,3,2),dtype=int32)
        cell_tris[0,:,0] = [0,1,1]
        cell_tris[0,:,1] = [0,0,1]
        cell_tris[1,:,0] = [1,0,0]
        cell_tris[1,:,1] = [1,1,0]
        cell_neighbors = zeros((2,3,4),dtype=int32)
        # triangle 0, itype = 0
        cell_neighbors[0,0,:] = [ 1, 0, 1, 0]
        cell_neighbors[0,1,:] = [ 0, 0, 1, 1]
        cell_neighbors[0,2,:] = [ 0,-1, 1, 2]
        # triangle 1, itype = 0
        cell_neighbors[1,0,:] = [-1, 0, 0, 0]
        cell_neighbors[1,1,:] = [ 0, 0, 0, 1]
        cell_neighbors[1,2,:] = [ 0, 1, 0, 2]

        cell_zeros = zeros((2,3),dtype=int32)
        cell_ones = cell_zeros + 1

        n1 = ni 
        n2 = nj
        
        n1_1=n1-1
        n2_1=n2-1
        
        ntri = 2*(n1-1)*(n2-1)
        tri_vertices=zeros((ntri,3),dtype=int32)
        tri_neighbors = zeros((ntri,3),dtype=int32)
        tri_ij = zeros((ntri,2),dtype=int32)

        for j in range(n2_1):
            jneigh = j + cell_neighbors[:,:,1]
            for i in range(n1_1):
                iv  = i + j * n1_1
                ineigh = i + cell_neighbors[:,:,0]                
                ivv2 = i2([ineigh,jneigh],n1_1,n2_1,False,False)
                ineighbor = ( where(ivv2 >= 0,(ivv2*2+cell_neighbors[:,:,2]),cell_zeros)
                             -where(ivv2 < 0,cell_ones,cell_zeros)
                            )
                if any(tri_neighbors[2*iv:2*iv+2,:].shape != ineighbor.shape):
                    print(i,j,iv,ineigh,jneigh,ineighbor,n1,n2,2*iv,2*iv+2)
                    print(tri_neighbors.shape)
                    print('tri_neighbors.shape: ',tri_neighbors[2*iv:2*iv+2,:].shape,
                          'ineighbor.shape: ', ineighbor.shape)
                tri_neighbors[2*iv:2*iv+2,:] = ineighbor
                tri_vertices[2*iv:2*iv+2,:] = i2([i+cell_tris[:,:,0],j+cell_tris[:,:,1]],n1,n2,False,False)
                tri_ij[2*iv:2*iv+2,:]=[[i,j],[i,j]]

        trifile_out = Dataset(tri_file,'w',format='NETCDF4')
        trifile_out.Description = 'triangulation of GAMERA grid converted to cylindrical in X, R (radial distance) and azimuth PHI around X-axis'

        new_dim = trifile_out.createDimension('c1', ntri)
        new_dim = trifile_out.createDimension('c2', 3)
        new_dim = trifile_out.createDimension('c3', 2)

        tri_verts_out = trifile_out.createVariable('tri_vertices', int32, ('c1', 'c2'))
        tri_neighbors_out = trifile_out.createVariable('tri_neighbors', int32, ('c1', 'c2'))
        tri_ij_out = trifile_out.createVariable('tri_ij', int32, ('c1','c3'))

        tri_verts_out[:] = tri_vertices
        tri_verts_out.units=''
        tri_verts_out.Description = 'vertex indices for each triangle in X-R grid'
        tri_neighbors_out[:] = tri_neighbors
        tri_neighbors_out.units=''
        tri_neighbors_out.Description = 'triangle number of neighboring triangles'
        tri_ij_out[:] = tri_ij
        tri_ij_out.units=''
        tri_ij_out.Dscription = 'i,j position of each triangle in X-R grid'
        # close file
        trifile_out.close()                
    else:
        tri_file_ds = Dataset(tri_file)
        tri_vertices = array(tri_file_ds.variables['tri_vertices'])
        tri_neighbors = array(tri_file_ds.variables['tri_neighbors'])
        tri_ij = array(tri_file_ds.variables['tri_ij']) 

    return([tri_vertices,tri_neighbors,tri_ij])

#
# Cartesian overlay grid with triangle index estimate for for LFM and GAMERA grids
#
def PrepareOverlaySphGrid2D(X2D,Y2D,tri_vertices):
    ''' Prepare overlay grid to estimate strat triangle index in a 2D grid that is logically arranged like a spherical grid
    we are given X and Y and convert into theta (angle away from positive X-axis) and r (distance from origin)
    the 2D grid is assumed to be nearly soherically arranged near the origin with a hole of radius 2 and elliptical in shape on the outs=ide
    we assume that rows of grid position are arranged at nearly constant elevation angle theta '''

    xmin = X2D.min()
    xmax = X2D.max()
    ymin = Y2D.min()
    ymax = Y2D.max()

    nj,ni = X2D.shape

    dx = 0.5 # nominal spacing of overlay grid  
    dy = 0.5 # to resolve hole near Earth of radius 2 R_E
    nx = int((xmax-xmin)/dx) # that is a fairly fine grid but ensures we should get good start triangle positons in all cases
    ny = int((ymax-ymin)/dy) # the final triangle position should be reached within 10 steps or so (dr>=0.2 for quad-resolutiion grid)
    dx = (xmax-xmin)/nx
    dy = (ymax-ymin)/ny
    
    R2D = sqrt(X2D*X2D + Y2D*Y2D)
    TH2D = arccos(X2D/R2D)

    # we use X=THETA and Y=RADIUS in spherical coordinates
    # we use THETA first and tehn extract columns of RADIUS (Y2D) to get an estimate of (j,i) grid cell to set start triangle
    start_index = zeros([ny,nx],dtype=int32)

    for ix in range(nx):
        x_ = xmin+dx*(ix+0.5)
        for iy in range(ny):
            y_ = ymin+dy*(iy+0.5)
            
            r_ = sqrt(x_ * x_ + y_ * y_)
            th_ = arccos(x_/r_)
            
            w_j = [j for j in range(nj-1) if TH2D[j,0] <= y_ and TH2D[j+1,0] > th_]
            if len(w_j) == 0:
                if th_ < TH2D[1,0]:
                    j=0
                if th_ >= TH2D[nj-2,0]:
                    j=nj-2
            else:
                j = w_j[0]

            # two radial grid position arrays that bracket the above elevation angle THETA
            R1D_0 = R2D[j,:].flatten()
            R1D_1 = R2D[j+1,:].flatten()

            w_i = [i for i in range(ni-1) if R1D_0[i] <= r_ and R1D_0[i+1] > r_]
            if len(w_i) == 0:
                # an overlay cell's corners are 0.5*dx and 0.5*dy away frm the center
                # adding and subtracting 0.6*dx and 06*dy to the center position in either direction
                # ensures that a start triangle is found when a cormer is within the model domain 
                # if that enlarged cell is still entirely outside we can assume that
                # any position in the original overlay cell is outsode the simul;ation domain
                dx_ = dx*0.6 
                dy_ = dy*0.6
                x_corners = array([x_+dx_,x_+dx_,x_+dx_,x_-dx_])
                y_corners = array([y_+dy_,y_+dy_,y_+dy_,y_-dy_])
                r_corners = sqrt(x_corners*x_corners+y_corners*y_corners)
                if r_ < R1D_0[0]:
                    if all(r_corners < R1D_0[0]) and all(r_corners < R1D_1[0]):
                        i = -1
                    else:
                        i = 0
                if r_ > R1D_0[ni-1]:
                    if all(r_corners > R1D_0[ni-1]) and all(r_corners > R1D_1[ni-1]):
                        i = -1 
                    else:
                        i = ni-2
            else:
                i = w_i[0]

            if i >= 0 and j >= 0:
                # return index of a triangle that is in or near the overlay cell
                start_index[iy,ix] = 2 * (i + (ni-1)*j) #  triangulation has (ni-1) * (nj-1) cells with 2 triangles each
            else:
                start_index[iy,ix] = -1 # overlay cells fully outside of ellipsoid with cutout sphere

    return(xmin,dx,nx,ymin,dy,ny,start_index)

#
# cartesian overlay grid with estimate of triangle index for arbitrarily shaped X2D and Y2D in space
#
def PrepareOverlayGrid2D(X2D,Y2D,tri_vertices):
    xmin = X2D.min()
    xmax = X2D.max()
    ymin = Y2D.min()
    ymax = Y2D.max()
    ncells = prod(X2D.shape)
    nj,ni = X2D.shape
    nx = ni-1
    ny = nj-1
    dx = (xmax - xmin)/nx
    dy = (ymax - ymin)/ny

    start_index = zeros([ny,nx],dtype=int32)

    average_cell_area = dx * dy / ncells
    average_cell_length = sqrt(average_cell_area)
    nx = int(dx/average_cell_length/5)
    ny = int(dy/average_cell_length/5)

    X2D_flat = X2D.flatten()
    Y2D_flat = Y2D.flatten()
    x_tri = array(X2D_flat[tri_vertices])
    y_tri = array(Y2D_flat[tri_vertices])
    ix_tri = floor((x_tri-xmin)/dx).astype(int32)
    iy_tri = floor((y_tri-ymin)/dy).astype(int32)
    min_ix = min(ix_tri,axis=1) # minimum overlay grid position among triangle vertices in X
    min_iy = min(iy_tri,axis=1) # maximum ...
    max_ix = max(ix_tri,axis=1) # minimum overlay grid position among triangle vertices in Y
    max_iy = max(iy_tri,axis=1) # maximum ...
    for ix in range(nx):
        for iy in range(ny):    
            if start_index[iy,ix] == -1:
                start_index_area = start_index[max([0,iy-1]):min([iy+1,ny-1])+1,max([0,ix-1]):min([ix+1,nx-1])+1]
                start_index[iy,ix] = max(start_index_area)
            if start_index[iy,ix] < 0:
                print('setup overlay grid : ix,iy,conter_index: ',ix,iy,start_index[iy,ix])
                
    return(xmin,dx,nx,ymin,dy,ny,start_index)


def is_inside_tri(position,i_tri,grid,debug=False):
    X2D = grid['X2D'].flatten()
    R2D = grid['R2D'].flatten()
    vertices = grid['tri_vertices']
    neighbors = grid['tri_neighbors']
    ntri = grid['ntri']
    
    if i_tri < 0 or i_tri > ntri:
        return(False,None,None)

    x_tri = X2D[vertices[i_tri,:]]
    y_tri = R2D[vertices[i_tri,:]]

    v01 = array([x_tri[1]-x_tri[0],y_tri[1]-y_tri[0]])
    v12 = array([x_tri[2]-x_tri[1],y_tri[2]-y_tri[1]])
    v20 = array([x_tri[0]-x_tri[2],y_tri[0]-y_tri[2]])

    f0 = array([1,2])
    f1 = array([2,0])
    f2 = array([0,1])
    c0 = array([sum(x_tri[f0]),sum(y_tri[f0])])/2.
    c1 = array([sum(x_tri[f1]),sum(y_tri[f1])])/2.
    c2 = array([sum(x_tri[f2]),sum(y_tri[f2])])/2.
                    
    tri_c=array([sum(x_tri),sum(y_tri)])/3. 
    p0 = position-c0
    p1 = position-c1
    p2 = position-c2

    area0= array([-v12[1],v12[0]])
    area1= array([-v20[1],v20[0]]) 
    area2= array([-v01[1],v01[0]])

    vol = -sum(v01*area0)
    if vol == 0.:
        vol=1e-5
    if vol < 0:
        sign = -1
    else:
        sign = 1
    
    vol0 = sum( p0*area0)
    vol1 = sum( p1*area1)
    vol2 = sum( p2*area2)
        
    weights = array([vol0,vol1,vol2])/vol

    w = [i for i in range(3) if weights[i] < 0]
    w_comp = [i for i in range(3) if weights[i] >= 0]
    if debug > 0:
        print('volume,sign: ',vol,sign)
        print('w: ',w)
        print('w_comp: ',w_comp)
        print("tri_vertices: ",vertices[i_tri,:])
        print('weights: ',weights)

    if len(w) == 0:
        is_inside=True
        i_neighbor = None

    if len(w) == 1:
        i_neighbor = neighbors[i_tri,w[0]]
        is_inside=False

    if len(w) == 2:
        i_neighbor = neighbors[i_tri,w[0]]
        c_to_w = array([x_tri[w_comp[0]],y_tri[w_comp[0]] ]) - tri_c
        c_to_p = position-tri_c
        # cross product of c_to_w,c_to_p
        w_cross_p = c_to_w[0]*c_to_p[1] - c_to_w[1]*c_to_p[0]
        #  i_neighbor = tri_neighbors[w[w_cross_p gt 0],i_tri]
        if w_cross_p < 0:
            i_neighbor = neighbors[i_tri,w[1]]
        else:
            i_neighbor = neighbors[i_tri,w[0]]
        is_inside=False

    if len(w) == 3:
        is_inside = False
        i_neighbor = -2
        
    return(is_inside,weights,i_neighbor)    
            
def find_triangle(position,old_tri,grid,debug=False):
    tri_dx_sg = grid['dx_sg']
    tri_dy_sg = grid['dr_sg']
    tri_xmin_sg = grid['xmin_sg'] # -300 or so
    tri_ymin_sg = grid['rmin_sg'] # 0        
    sg_start_index = grid['start_index_sg']
    tri_nx_sg = grid['nx_sg']
    tri_ny_sg = grid['nr_sg']

    n_tri = grid['ntri']

    tri_n1 = grid['n1']
    tri_n2 = grid['n2']
    tri_x = grid['X2D']
    tri_y = grid['R2D']
    tri_vertices = grid['tri_vertices']
    tri_ij = grid['tri_ij']
    
    if old_tri is None:
        old_tri = -1
    if old_tri == -1:
        i_sg = int((position[0]-tri_xmin_sg)/tri_dx_sg)
        j_sg = int((position[1]-tri_ymin_sg)/tri_dy_sg)
        if i_sg < 0 or i_sg > (tri_nx_sg-1):
            return(-1,Null)
        if j_sg < 0 or j_sg > (tri_ny_sg-1):
            return(-1,Null)
        old_tri = sg_start_index[j_sg,i_sg]
            
    i_neighbor = -1
    is_inside,weights,i_neighbor = is_inside_tri(position,old_tri,grid,debug=debug)
    if is_inside:
        return([old_tri,weights])

    n_tries = 1
    n_tries_max = min([n_tri,4*sqrt(n_tri)])
    tri_x_flat = array(list(tri_x)).flatten()
    tri_y_flat = array(list(tri_y)).flatten()
    
    if debug:
        n_tries_max = 10
        print('i_neighbor,n_tries_max: ',i_neighbor,n_tries_max)
    while is_inside is False and i_neighbor > 0 and n_tries < n_tries_max:
        if debug:
            print('n_tries: ',n_tries,' old_tri: ',old_tri,' is_inside: ',is_inside,' neighbor: ',i_neighbor)
            if old_tri > -1:
                print('position: ',position,' tri_vertices: ',
                      tri_x_flat[tri_vertices[old_tri,:]],
                      tri_y_flat[tri_vertices[old_tri,:]] )
        old_tri = i_neighbor
        n_tries = n_tries+1
        is_inside,weights,i_neighbor = is_inside_tri(position,old_tri,grid,debug=debug)
       
    if is_inside:
        return([old_tri,weights])
    else:
        return([False,None])


def Interpolate2D1D(Xvec,data3D,grid,debug=False,verbose=False):
    from time import perf_counter

    tic = perf_counter()
    
    x,y,z = array(Xvec)
    x_shape = x.shape
    x = x.flatten()
    y = y.flatten()
    z = z.flatten()
    [nx] = x.shape
    output_data = zeros(nx)
    
    r = sqrt(y*y+z*z)
    phi = arctan2(z,y)

    ij = grid['tri_ij']
    neighbors=grid['tri_neighbors']
    vertices=grid['tri_vertices']
    tri_ij = grid['tri_ij']
    ntri,two = ij.shape

    X2D = grid['X2D'] #.flatten()
    R2D = grid['R2D']   # .flatten()
    PHI = grid['PHI']     #.flatten()
    data_shape = data3D.shape
    (nk,nj,ni) = data_shape
    x2d_shape = X2D.shape

    if not any(array(data_shape[1:3])-array(x2d_shape)):
        interp_order = 1
    else:
        interp_order = 0

    old_tri = -1
    nk = len(PHI)

    nxvec=len(r)
    for ix in range(nxvec):
        output_data[ix] = NaN
        phi_ = phi[ix]
        if phi_ >= PHI[0] and phi_ <= PHI[-1]:
            k = 0
            while PHI[k+1] <= phi_ and k<nk-2:
                k = k+1
            wz0 = (PHI[k+1]-phi_)/(PHI[k+1]-PHI[k])
            wz1 = 1. - wz0

            i_tri,weights = find_triangle(array([x[ix],r[ix]]),old_tri,grid,debug=debug)
            if i_tri is not False and i_tri != -1:
                v_ = vertices[i_tri,:]
                i,j = tri_ij[i_tri,:]

                if interp_order == 1:
                    # implement interpolations in 2D triangle and the perpendicular coordinate
                    j = array(v_/ni,dtype=int32)
                    i = array(mod(v_, ni),dtype=int32)
                    output_data[ix] = (
                        + wz0 * sum(weights*[data3D[k,  j[l],i[l]] for l in range(3)])
                        + wz1 * sum(weights*[data3D[k+1,j[l],i[l]] for l in range(3)])
                    )
                else:
                    # return value at cell ceenter
                    i_ = min(i)
                    j_ = min(j)
                    output_data[ix] = data3D[k,j_,i_]                
        else:
            print('phi out of range: ',ix,phi_,PHI[0],PHI[-1])
    toc = perf_counter()
    if verbose:
        print('interpolation completed after ',toc-tic,' seconds')
    return(reshape(output_data,x_shape))
        
            
