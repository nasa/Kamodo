#include "tri_extern.h"
#include <stdlib.h>

void setup_tri_pointers(// 2D grid decomposed into triangles
			float *tri_x_in,
			float *tri_y_in,
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
			// 3rd coordinate perpendicular to 2D triangulated grid
			float *tri_z_in,
			IDL_LONG *tri_nz_in
			){
  /* 3D grid */
  tri_x = tri_x_in;
  tri_y = tri_y_in;
  tri_z = tri_z_in;
  tri_nx = tri_nx_in[0];
  tri_ny = tri_ny_in[0];
  tri_nz = tri_nz_in[0];
  /* triangulation in 2D */
  n_tri = n_tri_in[0];
  tri_vertices = tri_vertices_in;
  tri_ij = tri_ij_in;
  tri_neighbors = tri_neighbors_in;
  /* overlay grid */
  tri_start_index_sg = tri_start_index_sg_in;
  tri_xmin_sg = tri_xmin_sg_in[0];
  tri_ymin_sg = tri_ymin_sg_in[0];
  tri_dx_sg = tri_dx_sg_in[0];
  tri_dy_sg = tri_dy_sg_in[0];
  tri_nx_sg = tri_nx_sg_in[0];
  tri_ny_sg = tri_ny_sg_in[0];
  
  if (tri_x1d != NULL){free(tri_x1d);}
  tri_x1d = (float*)malloc(tri_nx*sizeof(float));

  
}  

