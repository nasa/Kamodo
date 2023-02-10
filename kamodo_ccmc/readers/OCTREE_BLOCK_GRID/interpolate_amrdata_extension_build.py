import os
from cffi import FFI
ffibuilder = FFI()

ffibuilder.cdef("""
typedef int IDL_LONG;
typedef struct {int refinement_level;
                int parent_ID; int child_count; int child_IDs[8];
                float XMIN, XMAX, XCenter, YMIN, YMAX, YCenter, ZMIN, ZMAX, ZCenter;
               } octree_block;
int find_octree_block(float x, float y, float z, int old_blocknumber, int max_level);
int xyz_ranges(int N, float *x, float *y, float *z, float *XR_BLK, float *YR_BLK, float *ZR_BLK,
               float *X_BLK, float *Y_BLK, float *Z_BLK, float *box_range,
               int *NX_in, int *NY_in, int *NZ_in, int positions_in_cell_center);
void setup_octree_pointers(
			IDL_LONG MAX_AMRLEVELS_in,
			octree_block *octree_blocklist_in, 
			IDL_LONG *numparents_at_AMRlevel_in,
			IDL_LONG *block_at_AMRlevel_in);
void   setup_octree(IDL_LONG N_blks_in,
		    float *xr_blk, float *yr_blk, float *zr_blk,
		    IDL_LONG MAX_AMRLEVELS, float *box_range,
		    octree_block *octree_blocklist_in, 
		    IDL_LONG N_octree_blocks,
		    IDL_LONG *numparents_at_AMRlevel_in,
		    IDL_LONG *block_at_AMRlevel_in);
float interpolate_amrdata(float xx,float yy,float zz,
                          float *field, /* fields data array */
                          IDL_LONG is_new_position);
int interpolate_amrdata_multipos(float *xx,float *yy,float *zz, int npos,
                          float *field, float *output);
int trace_fieldline_octree(float x_start,float y_start,float z_start, float r_end,
                           float *bxfield, float *byfield, float *bzfield,
                           float *flx, float *fly, float *flz, IDL_LONG *step_max,
                           float dn, float bdp, float *tilt, float spherical_deg2rad);
""")


# on Windows (os.name == 'nt'), no explicit list of libraries is needed
libraries=['']    
# link with the math library (libm) in Unix-style OS (Linux, Mac-OS)
if os.name == 'posix':
    libraries=['m']

ffibuilder.set_source("_interpolate_amrdata",  # name of the output C extension
"""
// always add any public function declaration in the separate ffibuilder.cdef declaration!
    #include <stdlib.h>
    #include "fl_extern.h"
    #include "octree_block.h"
    #include "setup_octree.h"
    int xyz_ranges(int N, float *x, float *y, float *z,
                   float *XR_BLK, float *YR_BLK, float *ZR_BLK,
                   float *X_BLK, float *Y_BLK, float *Z_BLK,float
                   *box_range,
                   int *NX_in,
                   int *NY_in,
                   int *NZ_in,
                   int positions_in_cell_center)
    {
      int ib,ix=0,iy=0,iz=0,NV;
      float XMIN_blk, YMIN_blk, ZMIN_blk, XMAX_blk, YMAX_blk, ZMAX_blk;
      float XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX;
      float dx_cell_max=0., dy_cell_max=0., dz_cell_max=0., dx_cell_min=1e10, dy_cell_min=1e10, dz_cell_min=1e10;

      if (XR_BLK == NULL){return(-1);}
      if (YR_BLK == NULL){return(-1);}
      if (ZR_BLK == NULL){return(-1);}

      if (x == NULL){return -2;}
      if (y == NULL){return -2;}
      if (z == NULL){return -2;}

      if (X_BLK == NULL){return(-3);}
      if (Y_BLK == NULL){return(-3);}
      if (Z_BLK == NULL){return(-3);}

      if (NX_in == NULL){return(-4);}
      if (NY_in == NULL){return(-4);}
      if (NZ_in == NULL){return(-4);}
// assign to common block pointers for later use
      x_blk = X_BLK;
      y_blk = Y_BLK;
      z_blk = Z_BLK;

      NX = 0;
      while (ix<N && NX==0){
        if (fabs(x[0] - x[ix]) < TOLERANCE){
	  NX = ix;
	}
        ix++;
      }
      if (NX == 0){return(-5);}
      NY = 0;
      ix = NX;
      while (ix<N && NY==0){
        if (fabs(y[0] - y[ix]) < TOLERANCE){
	  NY = ix/NX;
	}
        ix += NX;
      }
      if (NY == 0){return(-5);}

//      NY=NX; // this should work in (nearly) 
//      NZ=NX; // all cases at the CCMC

// find block size in Z
      NZ = 1;
      long* where_z_eq_z0=(long*)malloc(N*sizeof(long));
      long N_z_eq_z0 = 0;
      for (iz=0; iz<N; iz++){
	if (fabs(z[0]- z[iz])< TOLERANCE){
	  where_z_eq_z0[N_z_eq_z0]=iz;
	  N_z_eq_z0++;
	}
      }
      if (N_z_eq_z0==0){
	NZ = 1;
      } else {
	NZ = N;
	for (iz=1; iz<N_z_eq_z0; iz++){
	  if ( (where_z_eq_z0[iz]-where_z_eq_z0[iz-1]) < NZ){
	    NZ = where_z_eq_z0[iz]-where_z_eq_z0[iz-1];
	  }
	}
      }
      if (NZ == 0){return(-5);}
      if (NZ == 1){
          NZ = NX;
#ifdef debug
          fprintf(stderr,
                  "xyz_ranges WARNING: NZ was found to be 1 and set equal to NX (NZ=%i).\\nNY=%i\\n",
                  NZ, NY);
#endif
if (NZ==1){return(-6);}
      }
      free(where_z_eq_z0);
#ifdef debug
      fprintf(stderr,"Block size NX=%i  NY=%i NZ=%i\\n", NX, NY, NZ);
#endif
      NV = NX*NY;
      NVV = NX*NY*NZ;
      N_blks = N/NVV;
/* return block shape */
      NX_in[0] = NX;
      NY_in[0] = NY;
      NZ_in[0] = NZ;
/* get box range */
      XMIN = 1e30;
      XMAX =-1e30;
      YMIN = 1e30;
      YMAX =-1e30;
      ZMIN = 1e30;
      ZMAX =-1e30;

      for (ib=0; ib<N_blks; ib++){
        float dx_cell, dy_cell, dz_cell;
        dx_cell = x[NVV*ib+1]-x[NVV*ib];
        XMIN_blk = x[NVV*ib]-positions_in_cell_center*0.5*dx_cell-TOLERANCE;
        XMAX_blk = x[NVV*ib+NX-1]+positions_in_cell_center*0.5*dx_cell-TOLERANCE;
        XR_BLK[2*ib] = XMIN_blk;
        XR_BLK[2*ib+1] = XMAX_blk;
        if(XMIN > XMIN_blk){XMIN = XMIN_blk;}
        if(XMAX < XMIN_blk){XMAX = XMAX_blk;}
        if (dx_cell_max < dx_cell){dx_cell_max = dx_cell;}
        if (dx_cell_min > dx_cell){dx_cell_min = dx_cell;}

        dy_cell = y[NVV*ib+NX]-y[NVV*ib];
        YMIN_blk = y[NVV*ib]-positions_in_cell_center*0.5*dy_cell-TOLERANCE;
        YMAX_blk = y[NVV*ib+NV-1]+positions_in_cell_center*0.5*dy_cell-TOLERANCE;
        YR_BLK[2*ib] = YMIN_blk;
        YR_BLK[2*ib+1] = YMAX_blk;
        if(YMIN > YMIN_blk){YMIN = YMIN_blk;}
        if(YMAX < YMIN_blk){YMAX = YMAX_blk;}

        if (dy_cell_max < dy_cell){dy_cell_max = dy_cell;}
        if (dy_cell_min > dy_cell){dy_cell_min = dy_cell;}

        dz_cell = z[NVV*ib+NV]-z[NVV*ib];
        ZMIN_blk = z[NVV*ib]-positions_in_cell_center*0.5*dz_cell-TOLERANCE;
        ZMAX_blk = z[NVV*ib+NVV-1]+positions_in_cell_center*0.5*dz_cell-TOLERANCE;
        ZR_BLK[2*ib] = ZMIN_blk;
        ZR_BLK[2*ib+1] = ZMAX_blk;
        if(ZMIN > ZMIN_blk){ZMIN = ZMIN_blk;}
        if(ZMAX < ZMIN_blk){ZMAX = ZMAX_blk;}
        if (dz_cell_max < dz_cell){dz_cell_max = dz_cell;}
        if (dz_cell_min > dz_cell){dz_cell_min = dz_cell;}
      }
// populate x_blk, y_blk and z_blk arrays
      for (ib=0; ib<N_blks; ib++){
        for (ix=0; ix<NX; ix++){
           x_blk[ib*NX+ix] = x[NVV*ib+ix];
        }
        for (iy=0; iy<NY; iy++){
           y_blk[ib*NY+iy] = y[NVV*ib+NX*iy];
        }
        for (iz=0; iz<NZ; iz++){
           z_blk[ib*NZ+iz] = z[NVV*ib+NV*iz];
        }
      }
// populate box_range array
      box_range[0] = XMIN;
      box_range[1] = XMAX;
      box_range[2] = YMIN;
      box_range[3] = YMAX;
      box_range[4] = ZMIN;
      box_range[5] = ZMAX;

#ifdef debug
      fprintf(stderr,
              "Number of blocks: %i  min(dx_cell): %f  max(dx_cell): %f\\nmin(dy_cell): %f  max(dy_cell): %f\\nmin(dz_cell): %f  max(dz_cell): %f\\n",
              N_blks, dx_cell_min, dx_cell_max, dy_cell_min, dy_cell_max, dz_cell_min, dz_cell_max);
      fprintf(stderr,"Box range: X: %f %f   Y: %f %f  Z: %f %f\\n", XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX);
#endif
      return(N_blks);
    }

    int interpolate_amrdata_multipos(float *xx, float *yy, float *zz, int npos, 
                                     float *data, float *output){
        int ipos;
        for (ipos=0; ipos<npos; ipos++){
            output[ipos]=interpolate_amrdata(xx[ipos], yy[ipos], zz[ipos], data, 1);
        }
        return(0);
    }

""",
    sources=['interpolate_amrdata.c', 'setup_parent.c', 'setup_octree.c', 'interpolate_in_block.c',
             'find_octree_block.c','find_in_block.c','find_block.c','trace_fieldline.c'],
    libraries=libraries)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True,debug=False)
    #ffibuilder.compile(verbose=True)
    #ffibuilder.compile()

