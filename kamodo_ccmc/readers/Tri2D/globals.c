/*
 * globals.c - Global variable definitions for Tri2D library
 *
 * These global variables are declared as "extern" in tri_extern.h and used
 * throughout the Tri2D code. They need to be defined exactly once in
 * the shared library. In the original CFFI-based build, these were defined
 * in the generated wrapper code; for the ctypes build we define them here.
 *
 * Author: Kamodo team (extracted for ctypes build)
 */

#include <stddef.h>  /* for NULL */

typedef int IDL_LONG;

/* Grid arrays: X,Y is triangulated, z,u (v,w) additional dimensions */
float *tri_x = NULL;
float *tri_y = NULL;
float *tri_z = NULL;
float *tri_u = NULL;
float *tri_v = NULL;
float *tri_w = NULL;
float tri_missing = -1e30f;
float *tri_x1d = NULL;

IDL_LONG tri_nx = 0;
IDL_LONG tri_ny = 0;
IDL_LONG tri_nz = 0;
IDL_LONG tri_nu = 0;
IDL_LONG tri_nv = 0;
IDL_LONG tri_nw = 0;

/* Overlay grid in 2D (X,Y) */
float tri_xmin_sg = 0.0f;
float tri_ymin_sg = 0.0f;
float tri_dx_sg = 0.0f;
float tri_dy_sg = 0.0f;

IDL_LONG *tri_nelems_in_sg = NULL;
IDL_LONG *tri_vertices = NULL;
IDL_LONG *tri_start_index_sg = NULL;
IDL_LONG *start_index_sg = NULL;
IDL_LONG *end_index_sg = NULL;
IDL_LONG *tri_ij = NULL;
IDL_LONG *tri_neighbors = NULL;
IDL_LONG n_tri = 0;
IDL_LONG tri_nx_sg = 0;
IDL_LONG tri_ny_sg = 0;
IDL_LONG have_outside_sg = 0;
