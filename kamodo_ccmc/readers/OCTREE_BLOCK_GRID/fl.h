/************************************************/
/************************************************/
/*						*/
/*     Global variables for program fl_int      */
/*						*/
/************************************************/
#ifndef __fl_h__
#define __fl_h__
#define no_idl

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

/************************************************/
#define max(a,b) ( (a)>(b) ? (a) : (b) )
#define min(a,b) ( (a)<(b) ? (a) : (b) )

/* supported IDL versions are determined by the 'make.sh' script:
   correct versions of IDL's 'export.h' are linked to the
   generic name 'IDL_export.h' */
#include "IDL_export.h"

/**********************************************/
/* using the octree structure of BATSRUS data */
/* see br_export.c                            */
/**********************************************/
#include "octree_block.h"

octree_block *octree_blocklist;

float *fields,*x_blk,*y_blk,*z_blk,*xr_blk,*yr_blk,*zr_blk;
float MISSING;
IDL_LONG  NVV,NX,NX1,NX2,NY,NY1,NY2,NZ,NZ1,NZ2,nvar,NV_blk,NVV_blk,N_blks;

IDL_LONG  MAX_AMRLEVELS,parent_count,
    *block_at_AMRlevel,*numblocks_at_AMRlevel,*numparents_at_AMRlevel;


IDL_LONG setup_parent(IDL_LONG iblock,IDL_LONG ilev,
		      float XMIN, float YMIN, float ZMIN,
		      IDL_LONG N_blks,
		      IDL_LONG N_parents_max, 
		      IDL_LONG *N_parents);

IDL_LONG climb_octree(IDL_LONG root,float x, float y, float z,
		      IDL_LONG max_level);
IDL_LONG find_octree_block(float x,float y,float z,
			   IDL_LONG old_blocknumber,IDL_LONG max_level);
/* TOLERANCE allows for small numerical deviation from integer */
#define TOLERANCE 0.001

float p1,p2,p3, /* number of root blocks in box
                   should be within TOLERANCE of integer */
    dxmin_blk,dymin_blk,dzmin_blk,
    dxmin_cell,dymin_cell,dzmin_cell;

/* smart subsampling in AMR grid */
int grid_select(float *, int, int,
                float, float, float,
                float*, long, long);

/* interpolation in AMR block */
float interpolate_in_block(float, float, float,
                           float*,
                           float*,float*,float*,
                           IDL_LONG,IDL_LONG,IDL_LONG);
#ifndef no_idl
float interpolate_amrdata(float xx, float yy, float zz,
                          char *variable_name,
                          IDL_LONG VAR_ID, IDL_LONG NVAR,
                          float *field, IDL_LONG new_position);
#else
float interpolate_amrdata(float xx, float yy, float zz,
                          float *field, IDL_LONG new_position);
#endif
/***********************/
/* function prototypes */
/***********************/
void sami_grid_pos(double y00, double y01, double y10,
		   double z00, double z01, double z10,
		   double y_p, double z_p,
		   double *dy, double *dz);
float interpolate_sami_data(const float lon_,const float lat_,const float hh_,
			    char *variable_name,
			    IDL_LONG VAR_ID, /* offset in fields array */
			    IDL_LONG NVAR, /* stride in fields array */
			    float *field, IDL_LONG new_position,
			    float *dY, float *dZ);

float interpolate_lfm_data(float XX, float YY, float ZZ,
			   char *variable_name,
			   IDL_LONG VAR_ID, /* offset in fields array */
			   IDL_LONG NVAR, /* stride in fields array */
			   float *field, 
			   IDL_LONG new_position,
			   float *dx, 
			   float *dy);

IDL_LONG trace_fieldline_analytic_b(float x_start,float y_start,float z_start,
				    float r_end,
				    float *xr_blk, float *yr_blk, float *zr_blk,
				    float *flx, float *fly, float *flz, float *v_mag,
				    IDL_LONG *step_max,
				    float dn,
				    float mag_dipole_strength,
				    float *mag_dipole_axis,
				    float mirror_dipole_strength,
				    float *mirror_dipole_axis,
				    float *mirror_dipole_xyz,
				    float *b_sw);

IDL_LONG trace_fieldline_lfm(
		    float *flx, float *fly, float *flz, float *v_mag,
                    float *fields,
		    IDL_LONG *VAR_IDs,IDL_LONG N_VAR_IDs,IDL_LONG NVAR,
                    IDL_LONG *step_max, float r_end,
                    float dn, float bdp,float tilt);

IDL_LONG find_block(float,float,float,
                    float *,float *,float *,long,long);

void find_in_block(float, float, float, int*, int*, int*,
                   float* , float*, float*,
                   IDL_LONG, IDL_LONG, IDL_LONG);

void hunt(float *, IDL_LONG, float, int *);


IDL_LONG trace_fieldline(float x_start,float y_start,float z_start,
                         float r_end,
                         IDL_LONG NX, IDL_LONG NY, IDL_LONG NZ,
                         IDL_LONG N_blks,
                         float* x_blk, float *y_blk, float *z_blk,
                         float *xr_blk,float *yr_blk, float *zr_blk,
			 float *fields,IDL_LONG nf, IDL_LONG *v_indices,
                         float *flx, float *fly, float *flz, float *b_mag,
                         IDL_LONG *step_max,
                         float dn, float bdp,float *tilt);

IDL_LONG trace_fieldline_spherical_staggered(
                    float x_start,float y_start,float z_start,
                    float r_end,
                    IDL_LONG NX, IDL_LONG NY, IDL_LONG NZ, IDL_LONG N_blks,
                    float* x, float *y, float *z,
                    float* x_bx, float *y_by, float *z_bz,
                    float *xr,float *yr, float *zr,
		    float *fields,IDL_LONG nf, IDL_LONG *v_indices,
                    float *flx, float *fly, float *flz, IDL_LONG *step_max,
                    float dn, float deg2rad, IDL_LONG usePol);


IDL_LONG trace_fieldline_staggered(float x_start,float y_start,float z_start,
                    float r_end,
                    IDL_LONG NX, IDL_LONG NY, IDL_LONG NZ, IDL_LONG N_blks,
                    float *x, float *y, float *z,
                    float *x_bx, float *y_by, float *z_bz,
                    float *xr,float *yr, float *zr,
		    float *fields, IDL_LONG nf, IDL_LONG *v_indices,
                    float *flx, float *fly, float *flz, float *v_mag,
                    IDL_LONG *step_max,
                    float dn, float bdp,float *tilt);

IDL_LONG trace_fieldline_cdf(char *cdf_file,char *cdf_variables,
			     float *flx, float *fly, float *flz, float *v_mag,
			     float r_end,
			     IDL_LONG *step_max,
			     float dn, float bdp,float *tilt,float missing);

IDL_LONG trace_fieldline_octree(
    float x_start,float y_start,float z_start,
    float r_end,
    // IDL_LONG NX, IDL_LONG NY, IDL_LONG NZ, IDL_LONG N_blks,
    // float *fields, 
    //    IDL_LONG nf, //IDL_LONG *v_indices,
    float *bxfield,float *byfield,float *bzfield,
    float *flx, float *fly, float *flz, IDL_LONG *step_max,
    float dn, float bdp, float *tilt, float deg2rad);

IDL_LONG trace_fieldline_octree_spherical(
    float x_start,float y_start,float z_start,
    float r_end,
    // IDL_LONG NX, IDL_LONG NY, IDL_LONG NZ, IDL_LONG N_blks,
    // float *fields, 
    //    IDL_LONG nf, IDL_LONG *v_indices,
    float *bxfield, float *byfield, float *bzfield,
    float *flx, float *fly, float *flz, IDL_LONG *step_max,
    float dn, float bdp, float *tilt, float deg2rad);

IDL_LONG trace_particle_octree(
    float x_start,float y_start,float z_start,
    float r_end,
    float *b_x, float *b_y, float *b_z,
    float *v_x, float *v_y, float *v_z,
    float *flx, float *fly, float *flz,
    float *v_start,
    float *energies,
    float *part_time,
    IDL_LONG *step_max,
    float bdp, float tilt);

#endif
