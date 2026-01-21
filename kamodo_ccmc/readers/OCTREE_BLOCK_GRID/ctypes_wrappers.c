/*
 * ctypes_wrappers.c - Wrapper functions for ctypes-based Python interface
 *
 * These functions were originally embedded inline in the CFFI build script
 * (interpolate_amrdata_extension_build.py). They have been extracted here
 * for use with ctypes instead of CFFI.
 *
 * Author: Original CFFI code by Kamodo team, extracted for ctypes by Claude
 */

#include <stdlib.h>
#include <math.h>
#include "fl_extern.h"
#include "octree_block.h"
#include "setup_octree.h"

/*
 * xyz_ranges - Analyze block structure and compute ranges
 *
 * Parameters:
 *   N - Total number of data points
 *   x, y, z - Coordinate arrays
 *   XR_BLK, YR_BLK, ZR_BLK - Output block range arrays
 *   X_BLK, Y_BLK, Z_BLK - Output block coordinate arrays
 *   box_range - Output array for overall domain bounds
 *   NX_in, NY_in, NZ_in - Output block dimensions
 *   positions_in_cell_center - Flag indicating if positions are at cell centers
 *
 * Returns:
 *   Number of blocks on success, negative error code on failure
 */
int xyz_ranges(int N, float *x, float *y, float *z,
               float *XR_BLK, float *YR_BLK, float *ZR_BLK,
               float *X_BLK, float *Y_BLK, float *Z_BLK,
               float *box_range,
               int *NX_in, int *NY_in, int *NZ_in,
               int positions_in_cell_center)
{
    int ib, ix = 0, iy = 0, iz = 0, NV;
    float XMIN_blk, YMIN_blk, ZMIN_blk, XMAX_blk, YMAX_blk, ZMAX_blk;
    float XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX;
    float dx_cell_max = 0., dy_cell_max = 0., dz_cell_max = 0.;
    float dx_cell_min = 1e10, dy_cell_min = 1e10, dz_cell_min = 1e10;

    if (XR_BLK == NULL) { return (-1); }
    if (YR_BLK == NULL) { return (-1); }
    if (ZR_BLK == NULL) { return (-1); }

    if (x == NULL) { return -2; }
    if (y == NULL) { return -2; }
    if (z == NULL) { return -2; }

    if (X_BLK == NULL) { return (-3); }
    if (Y_BLK == NULL) { return (-3); }
    if (Z_BLK == NULL) { return (-3); }

    if (NX_in == NULL) { return (-4); }
    if (NY_in == NULL) { return (-4); }
    if (NZ_in == NULL) { return (-4); }

    /* Assign to common block pointers for later use */
    x_blk = X_BLK;
    y_blk = Y_BLK;
    z_blk = Z_BLK;

    NX = 0;
    while (ix < N && NX == 0) {
        if (fabs(x[0] - x[ix]) < TOLERANCE) {
            NX = ix;
        }
        ix++;
    }
    if (NX == 0) { return (-5); }

    NY = 0;
    ix = NX;
    while (ix < N && NY == 0) {
        if (fabs(y[0] - y[ix]) < TOLERANCE) {
            NY = ix / NX;
        }
        ix += NX;
    }
    if (NY == 0) { return (-5); }

    /* Find block size in Z */
    NZ = 1;
    long *where_z_eq_z0 = (long *)malloc(N * sizeof(long));
    long N_z_eq_z0 = 0;
    for (iz = 0; iz < N; iz++) {
        if (fabs(z[0] - z[iz]) < TOLERANCE) {
            where_z_eq_z0[N_z_eq_z0] = iz;
            N_z_eq_z0++;
        }
    }
    if (N_z_eq_z0 == 0) {
        NZ = 1;
    } else {
        NZ = N;
        for (iz = 1; iz < N_z_eq_z0; iz++) {
            if ((where_z_eq_z0[iz] - where_z_eq_z0[iz - 1]) < NZ) {
                NZ = where_z_eq_z0[iz] - where_z_eq_z0[iz - 1];
            }
        }
    }
    if (NZ == 0) { return (-5); }
    if (NZ == 1) {
        NZ = NX;
#ifdef debug
        fprintf(stderr,
                "xyz_ranges WARNING: NZ was found to be 1 and set equal to NX (NZ=%i).\nNY=%i\n",
                NZ, NY);
#endif
        if (NZ == 1) { return (-6); }
    }
    free(where_z_eq_z0);

#ifdef debug
    fprintf(stderr, "Block size NX=%i  NY=%i NZ=%i\n", NX, NY, NZ);
#endif

    NV = NX * NY;
    NVV = NX * NY * NZ;
    N_blks = N / NVV;

    /* Return block shape */
    NX_in[0] = NX;
    NY_in[0] = NY;
    NZ_in[0] = NZ;

    /* Get box range */
    XMIN = 1e30;
    XMAX = -1e30;
    YMIN = 1e30;
    YMAX = -1e30;
    ZMIN = 1e30;
    ZMAX = -1e30;

    for (ib = 0; ib < N_blks; ib++) {
        float dx_cell, dy_cell, dz_cell;

        dx_cell = x[NVV * ib + 1] - x[NVV * ib];
        XMIN_blk = x[NVV * ib] - positions_in_cell_center * 0.5 * dx_cell - TOLERANCE;
        XMAX_blk = x[NVV * ib + NX - 1] + positions_in_cell_center * 0.5 * dx_cell - TOLERANCE;
        XR_BLK[2 * ib] = XMIN_blk;
        XR_BLK[2 * ib + 1] = XMAX_blk;
        if (XMIN > XMIN_blk) { XMIN = XMIN_blk; }
        if (XMAX < XMIN_blk) { XMAX = XMAX_blk; }
        if (dx_cell_max < dx_cell) { dx_cell_max = dx_cell; }
        if (dx_cell_min > dx_cell) { dx_cell_min = dx_cell; }

        dy_cell = y[NVV * ib + NX] - y[NVV * ib];
        YMIN_blk = y[NVV * ib] - positions_in_cell_center * 0.5 * dy_cell - TOLERANCE;
        YMAX_blk = y[NVV * ib + NV - 1] + positions_in_cell_center * 0.5 * dy_cell - TOLERANCE;
        YR_BLK[2 * ib] = YMIN_blk;
        YR_BLK[2 * ib + 1] = YMAX_blk;
        if (YMIN > YMIN_blk) { YMIN = YMIN_blk; }
        if (YMAX < YMIN_blk) { YMAX = YMAX_blk; }
        if (dy_cell_max < dy_cell) { dy_cell_max = dy_cell; }
        if (dy_cell_min > dy_cell) { dy_cell_min = dy_cell; }

        dz_cell = z[NVV * ib + NV] - z[NVV * ib];
        ZMIN_blk = z[NVV * ib] - positions_in_cell_center * 0.5 * dz_cell - TOLERANCE;
        ZMAX_blk = z[NVV * ib + NVV - 1] + positions_in_cell_center * 0.5 * dz_cell - TOLERANCE;
        ZR_BLK[2 * ib] = ZMIN_blk;
        ZR_BLK[2 * ib + 1] = ZMAX_blk;
        if (ZMIN > ZMIN_blk) { ZMIN = ZMIN_blk; }
        if (ZMAX < ZMIN_blk) { ZMAX = ZMAX_blk; }
        if (dz_cell_max < dz_cell) { dz_cell_max = dz_cell; }
        if (dz_cell_min > dz_cell) { dz_cell_min = dz_cell; }
    }

    /* Populate x_blk, y_blk and z_blk arrays */
    for (ib = 0; ib < N_blks; ib++) {
        for (ix = 0; ix < NX; ix++) {
            x_blk[ib * NX + ix] = x[NVV * ib + ix];
        }
        for (iy = 0; iy < NY; iy++) {
            y_blk[ib * NY + iy] = y[NVV * ib + NX * iy];
        }
        for (iz = 0; iz < NZ; iz++) {
            z_blk[ib * NZ + iz] = z[NVV * ib + NV * iz];
        }
    }

    /* Populate box_range array */
    box_range[0] = XMIN;
    box_range[1] = XMAX;
    box_range[2] = YMIN;
    box_range[3] = YMAX;
    box_range[4] = ZMIN;
    box_range[5] = ZMAX;

#ifdef debug
    fprintf(stderr,
            "Number of blocks: %i  min(dx_cell): %f  max(dx_cell): %f\nmin(dy_cell): %f  max(dy_cell): %f\nmin(dz_cell): %f  max(dz_cell): %f\n",
            N_blks, dx_cell_min, dx_cell_max, dy_cell_min, dy_cell_max, dz_cell_min, dz_cell_max);
    fprintf(stderr, "Box range: X: %f %f   Y: %f %f  Z: %f %f\n", XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX);
#endif

    return (N_blks);
}

/*
 * interpolate_amrdata_multipos - Interpolate at multiple positions
 *
 * Wrapper function that calls interpolate_amrdata for each position.
 *
 * Parameters:
 *   xx, yy, zz - Arrays of coordinates
 *   npos - Number of positions
 *   data - Field data array
 *   output - Output array for interpolated values
 *
 * Returns:
 *   0 on success
 */
int interpolate_amrdata_multipos(float *xx, float *yy, float *zz, int npos,
                                 float *data, float *output)
{
    int ipos;
    for (ipos = 0; ipos < npos; ipos++) {
        output[ipos] = interpolate_amrdata(xx[ipos], yy[ipos], zz[ipos], data, 1);
    }
    return (0);
}
