/*#include "fl_extern.h" */
#include "tri_extern.h"
#include "stdio.h"
#include "math.h"

/* #define DEBUG */

/* index in flat array with (NI,NJ) shape */
/* -1 if any index exceeds bounds */
long int  ij2(long int i, long int j,long int k, long int NI, long int NJ){
  long int   ij2_tmp=j*NI+i;
  if (i < 0 || i > (NI-1) ||
      j < 0 || j > (NJ-1) ) {ij2_tmp=-1;}
  return(ij2_tmp);
}

/* vectors can be of any length 
   spatial dot-product has vectors of length N=3 */
float dotprod(float *v1, float *v2,int N){
  float prod;
  int i;
  prod=0.;
  for (i=0;i<N;i++){
    prod+=v1[i]*v2[i];
  }
  return(prod);
}

int is_inside_tri(float x, float y, IDL_LONG i_tri,
		  float *weights, IDL_LONG *i_neighbor,
		  IDL_LONG debug){
  float x_tri[3], y_tri[3];
  /* directions along edges and center points */ 
  float v01[2], v12[2], v20[2], c0[2], c1[2], c2[2];
  float area, areas[3], p0[2], p1[2], p2[2], area0[2], area1[2], area2[2];
  float tri_center_x, tri_center_y;
  float vol, vols[3];
  int ivert, count=0, inside_index, outside_index;
  if (i_tri < 0 || i_tri > n_tri){
    return(0);
  }
  for (ivert = 0; ivert < 3; ivert++){
    x_tri[ivert] = tri_x[tri_vertices[ivert+3*i_tri]];
    y_tri[ivert] = tri_y[tri_vertices[ivert+3*i_tri]];
  }

  /* vectors along edges */
  v01[0] = x_tri[1] - x_tri[0];
  v01[1] = y_tri[1] - y_tri[0];

  v12[0] = x_tri[2] - x_tri[1];
  v12[1] = y_tri[2] - y_tri[1];

  v20[0] = x_tri[0] - x_tri[2];
  v20[1] = y_tri[0] - y_tri[2];

  c0[0] = (x_tri[1] + x_tri[2])/2.; 
  c0[1] = (y_tri[1] + y_tri[2])/2.; 

  c1[0] = (x_tri[2] + x_tri[0])/2.; 
  c1[1] = (y_tri[2] + y_tri[0])/2.; 

  c2[0] = (x_tri[0] + x_tri[1])/2.; 
  c2[1] = (y_tri[0] + y_tri[1])/2.; 

  tri_center_x=(x_tri[0]+x_tri[1]+x_tri[2])/3.;
  tri_center_y=(y_tri[0]+y_tri[1]+y_tri[2])/3.;

  p0[0]=x-c0[0]; 
  p0[1]=y-c0[1];

  p1[0]=x-c1[0];
  p1[1]=y-c1[1];

  p2[0]=x-c2[0];
  p2[1]=y-c2[1];

    /*
      area and volume elements
    */

  area0[0] = -v12[1];
  area0[1] =  v12[0];

  area1[0] = -v20[1];
  area1[1] =  v20[0];

  area2[0] = -v01[1];
  area2[1] =  v01[0];
  
  vol    =-dotprod(v01,area0,2); /*-total(v01,crossp(v12,v23)) */
  vols[0]= dotprod(p0,area0,2);  /* total( p0*area0) */
  vols[1]= dotprod(p1,area1,2);  /* total( p1*area1) */
  vols[2]= dotprod(p2,area2,2);  /* total( p2*area2) */

  if (debug > 2){
    fprintf(stderr,"Postion and vertices: \n");
    fprintf(stderr,"  %f %f\n",x,y);
    fprintf(stderr,"  %f %f\n",x_tri[0],y_tri[0]);
    fprintf(stderr,"  %f %f\n",x_tri[1],y_tri[1]);
    fprintf(stderr,"  %f %f\n",x_tri[2],y_tri[2]);

  }
  if (debug > 1 ){
    fprintf(stderr,"Facecenter position vectors: \n");
    fprintf(stderr,"  %f %f\n",c0[0],c0[1]);
    fprintf(stderr,"  %f %f\n",c1[0],c1[1]);
    fprintf(stderr,"  %f %f\n",c2[0],c2[1]);

    fprintf(stderr,"Facecenter to position vectors: \n");
    fprintf(stderr,"  %f %f\n",p0[0],p0[1]);
    fprintf(stderr,"  %f %f\n",p1[0],p1[1]);
    fprintf(stderr,"  %f %f\n",p2[0],p2[1]);

    fprintf(stderr,"Area normal vectors: \n");
    fprintf(stderr,"  %f %f\n",area0[0],area0[1]);
    fprintf(stderr,"  %f %f\n",area1[0],area1[1]);
    fprintf(stderr,"  %f %f\n",area2[0],area2[1]);

    fprintf(stderr,"Volumes: \n");
    fprintf(stderr,"  %f %f %f %f\n",vol,vols[0],vols[1],vols[2]);
  }

  count=0;
  outside_index = -1;
  for (ivert=0;ivert<3;ivert++){
    weights[ivert]=vols[ivert]/vol;  /* these should be > 0 and add up to 1. */
    if (weights[ivert] >= 0){      
      count=count+1;
      inside_index=ivert;
    } else {
      outside_index=ivert;
    }
    if (count == 0) {i_neighbor[0]=-1;}
    if (count == 1){
      float c_to_p0,c_to_p1,c_to_w0,c_to_w1;
      c_to_p0 = x - tri_center_x;
      c_to_p1 = y - tri_center_y;
      c_to_w0 = x_tri[inside_index] - tri_center_x;
      c_to_w1 = y_tri[inside_index] - tri_center_y;
      if (debug){fprintf(stderr,"neighbors: %li %li %li\n",
			 tri_neighbors[3*i_tri],
			 tri_neighbors[3*i_tri+1],
			 tri_neighbors[3*i_tri+2]);
      }
      float det = (c_to_w0*c_to_p1-c_to_w1*c_to_p0);
      if (inside_index == 0){
	if (det < 0) {outside_index=1;} else {outside_index=2;}
      }
      if (inside_index == 1){
	if (det < 0) {outside_index=2;} else {outside_index=0;}
      }
      if (inside_index == 2){
	if (det < 0) {outside_index=0;} else {outside_index=1;}
      }
      i_neighbor[0] = tri_neighbors[3*i_tri + outside_index];
    }
    if (count == 2){i_neighbor[0] = tri_neighbors[outside_index+3*i_tri];}
  }
#ifdef DEBUG
  fprintf(stderr,"Element: %ld Vol: %f Vols: %f %f %f\n",i_tri,vol,vols[0],vols[1],vols[2]);
  fprintf(stderr,"Weights: %f %f %f\n",weights[0],weights[1],weights[2]);
  fprintf(stderr,"Inside: %ld  i_neighbor: %ld\n",(count == 3),i_neighbor[0]);
  fprintf(stderr,"Count: %ld\n",count);
#endif
  if (count == 3) { return(1); } else { return(0); }
}

/* #define DEBUG_FIND_TRI */

IDL_LONG find_tri(float xx, float yy,
		  IDL_LONG i_tri_old, float *weights,IDL_LONG debug){
  IDL_LONG i_tri, index, index_start, index_end, i_neighbor;
  long int i_sg, j_sg;
  long int ijk, ijk_stencil[9];
  long int n_tries, n_tries_max;
  int is_inside=0;

  if (i_tri_old < 0 || i_tri_old >= n_tri){    
    if (tri_start_index_sg != NULL){
      i_sg=floor((xx-tri_xmin_sg)/tri_dx_sg);
      j_sg=floor((yy-tri_ymin_sg)/tri_dy_sg);
      if (debug){
	fprintf(stderr,"SG grid start x: %f R: %f  dx: %f  dr: %f\n",
		tri_xmin_sg,tri_ymin_sg,tri_dx_sg,tri_dy_sg);
      }
      if (i_sg < 0 || i_sg > (tri_nx_sg-1) ) {return(-1);}
      if (j_sg < 0 || j_sg > (tri_ny_sg-1) ) {return(-1);}
      i_tri_old  = tri_start_index_sg[i_sg+tri_nx_sg*j_sg];
      if (debug){
	fprintf(stderr,"Start sg index: %li %li nx_sg: %li ny_sg: %li i_tri: %li\n",i_sg,j_sg,tri_nx_sg,tri_ny_sg,i_tri_old);
      }
      if (i_tri_old < 0){return(-1);}
    } else {
      i_tri_old = tri_nx/2 + tri_nx*(tri_ny/2); /* a triangle in the center of the (i,j) grid */
      if (debug){
	fprintf(stderr,"Start sg index (default): %li %li nx_sg: %li ny_sg: %li  i_tri: %li\n",i_sg,j_sg,tri_nx_sg,tri_ny_sg,i_tri_old);
      }
    }
  }
  n_tries = 0;
  n_tries_max = 10*sqrt(n_tri);
  if (debug > 3){n_tries_max=10;} 

  if (debug){fprintf(stderr,"Start triangle index: %li\n",i_tri_old);}
  is_inside = is_inside_tri(xx,yy,i_tri_old,weights,&i_neighbor,debug);

  while (is_inside == 0 && i_neighbor >= 0 && n_tries < n_tries_max){
    i_tri_old = i_neighbor;
    n_tries++;
    if (debug){fprintf(stderr,"Neighbor triangle index: %li ntry: %li\n",i_tri_old,n_tries);}
    is_inside = is_inside_tri(xx,yy,i_tri_old,weights,&i_neighbor,debug);
  }
  if (debug){fprintf(stderr,"Final triangle index: %li  Inside: %li\n",i_tri_old,is_inside);}

  if (is_inside == 0){return(-1);} else {return(i_tri_old);}

}

