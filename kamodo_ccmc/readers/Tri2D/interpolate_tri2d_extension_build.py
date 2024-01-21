import os
from cffi import FFI
ffibuilder = FFI()

ffibuilder.cdef("""
typedef int IDL_LONG;
void setup_tri_pointers(float *tri_x_in,
			float *tri_y_in,
			/* float tri_y_t_in, // transposed y2d array */
                        IDL_LONG *tri_nx_in,
			IDL_LONG *tri_ny_in,
			IDL_LONG *n_tri_in,
			IDL_LONG *tri_vertices_in,
			IDL_LONG *tri_ij_in,
			IDL_LONG *tri_neighbors_in,
			// overlay grid
			float *tri_xmin_sg_in,
			float *tri_ymin_sg_in,
			float *tri_dx_sg_in,
			float *tri_dy_sg_in,
			IDL_LONG *tri_nx_sg_in,
			IDL_LONG *tri_ny_sg_in,
			IDL_LONG *tri_start_index_sg_in,
			float *tri_z_in,
			IDL_LONG *tri_nz_in
                      );
int interpolate_tri2d_plus_1d_multipos(float *xx,float *yy,float *zz, int npos,
                          float *field, float *output, int data_in_center);
/*
int trace_fieldline_tri(float x_start,float y_start,float z_start, float r_end,
                           float *bxfield, float *byfield, float *bzfield,
                           float *flx, float *fly, float *flz, IDL_LONG *step_max,
                           float dn, float bdp, float *tilt, float spherical_deg2rad);
*/
""")


# on Windows (os.name == 'nt'), no explicit list of libraries is needed
libraries=[]    
# link with the math library (libm) in Unix-style OS (Linux, Mac-OS)
if os.name == 'posix':
    libraries=['m']

ffibuilder.set_source("_interpolate_tri2d",  # name of the output C extension
"""
// always add any public function declaration in the separate ffibuilder.cdef declaration!
    #include <stdlib.h>
    #include "tri_extern.h"

    int interpolate_tri2d_plus_1d_multipos(float *xx, float *yy, float *zz, int npos, 
                                     float *data, float *output, int data_in_center){
        int ipos;
        for (ipos=0; ipos<npos; ipos++){
            output[ipos]=interpolate_tri2d_plus_1d(xx[ipos], yy[ipos], zz[ipos], data, data_in_center);
        }
        return(0);
    }

""",
    sources=['interpolate_tri2d_plus_1d.c',
             'setup_tri.c','find_tri.c','hunt.c'],
    libraries=libraries)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True,debug=False)
    #ffibuilder.compile(verbose=True)
    #ffibuilder.compile()

