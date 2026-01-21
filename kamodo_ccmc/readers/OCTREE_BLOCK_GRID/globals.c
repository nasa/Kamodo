/*
 * globals.c - Global variable definitions for OCTREE_BLOCK_GRID library
 *
 * These global variables are declared as "extern" in fl_extern.h and used
 * throughout the OCTREE code. They need to be defined exactly once in
 * the shared library. In the original CFFI-based build, these were defined
 * in the generated wrapper code; for the ctypes build we define them here.
 *
 * Author: Kamodo team (extracted for ctypes build)
 */

#include <stddef.h>  /* for NULL */
#include "IDL_export.h"
#include "octree_block.h"

/* Octree block list pointer */
octree_block *octree_blocklist = NULL;

/* Field and coordinate arrays */
float *fields = NULL;
float *x_blk = NULL;
float *y_blk = NULL;
float *z_blk = NULL;
float *xr_blk = NULL;
float *yr_blk = NULL;
float *zr_blk = NULL;

/* Missing value indicator */
float MISSING = -1e30f;

/* Grid dimensions */
IDL_LONG NVV = 0;
IDL_LONG NX = 0;
IDL_LONG NX1 = 0;
IDL_LONG NX2 = 0;
IDL_LONG NY = 0;
IDL_LONG NY1 = 0;
IDL_LONG NY2 = 0;
IDL_LONG NZ = 0;
IDL_LONG NZ1 = 0;
IDL_LONG NZ2 = 0;
IDL_LONG nvar = 0;
IDL_LONG NV_blk = 0;
IDL_LONG NVV_blk = 0;
IDL_LONG N_blks = 0;

/* AMR level information */
IDL_LONG MAX_AMRLEVELS = 0;
IDL_LONG parent_count = 0;
IDL_LONG *block_at_AMRlevel = NULL;
IDL_LONG *numblocks_at_AMRlevel = NULL;
IDL_LONG *numparents_at_AMRlevel = NULL;

/* Root block layout and cell dimensions */
float p1 = 0.0f;
float p2 = 0.0f;
float p3 = 0.0f;
float dxmin_blk = 0.0f;
float dymin_blk = 0.0f;
float dzmin_blk = 0.0f;
float dxmin_cell = 0.0f;
float dymin_cell = 0.0f;
float dzmin_cell = 0.0f;
