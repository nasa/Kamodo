typedef int IDL_LONG;
extern void crossp(float factor, float *v1, float *v2, float *v_out);
extern float dotprod(float *v1, float *v2,int N);
/* grid arrays: X,Y is triangulated, z,u (v,w) additional dimensions */
extern float *tri_x,*tri_y,*tri_z,*tri_u,*tri_v,*tri_w,tri_missing,*tri_x1d;
extern IDL_LONG   tri_nx,tri_ny,tri_nz,tri_nu,tri_nv,tri_nw;
/* overlay grid in 2D (X,Y) */
extern float  tri_xmin_sg,tri_ymin_sg,tri_dx_sg,tri_dy_sg;
extern IDL_LONG *tri_nelems_in_sg,*tri_vertices,
  *tri_start_index_sg,*start_index_sg,*end_index_sg,
  *tri_ij,*tri_neighbors,
  n_tri,tri_nx_sg,tri_ny_sg,have_outside_sg;

extern int is_inside_tri(float x, float y, IDL_LONG i_tri,
		  float *weights, IDL_LONG *i_neighbor,
		  IDL_LONG debug);

extern long int ijk3(long int i, long int j,long int k, long int NI, long int NJ, long int NK);
extern long int ijkl(long int i, long int j, long int k, long int l, long int NI, long int NJ, long int NK, long int NL);

extern void hunt(float *xx, IDL_LONG n, float x, int *jlo);
extern IDL_LONG find_tri(float xx, float yy, 
		  IDL_LONG i_tri_old, float *weights,IDL_LONG debug);
extern int is_inside_tri(float x, float y, IDL_LONG i_tri,
		  float *weights, IDL_LONG *i_neighbor,
		  IDL_LONG debug);
extern float interpolate_tri2d_plus_1d(float xx, float yy, float zz, float *data, int data_in_center);
extern void setup_tri_pointers(// 2D grid decomposed into triangles
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
                        IDL_LONG *tri_nz_in);

extern int interpolate_tri2d_plus_1d_multipos(float *xx, float *yy, float *zz, int npos, 
					      float *data, float *output, int data_in_center);
