/*
 * ctypes_wrappers.c - Wrapper functions for ctypes-based Python interface
 *
 * These functions were originally embedded inline in the CFFI build script
 * (interpolate_tri2d_extension_build.py). They have been extracted here
 * for use with ctypes instead of CFFI.
 *
 * Extracted from interpolate_tri2d_extension_build.py for standalone compilation.
 */

#include <stdlib.h>
#include "tri_extern.h"

/*
 * interpolate_tri2d_plus_1d_multipos - Interpolate at multiple positions
 *
 * Wrapper function that calls interpolate_tri2d_plus_1d for each position.
 *
 * Parameters:
 *   xx, yy, zz - Arrays of coordinates
 *   npos - Number of positions
 *   data - Field data array
 *   output - Output array for interpolated values
 *   data_in_center - Flag indicating if data values are at cell centers
 *
 * Returns:
 *   0 on success
 */
int interpolate_tri2d_plus_1d_multipos(float *xx, float *yy, float *zz, int npos,
                                       float *data, float *output, int data_in_center)
{
    int ipos;
    for (ipos = 0; ipos < npos; ipos++) {
        output[ipos] = interpolate_tri2d_plus_1d(xx[ipos], yy[ipos], zz[ipos], data, data_in_center);
    }
    return (0);
}
