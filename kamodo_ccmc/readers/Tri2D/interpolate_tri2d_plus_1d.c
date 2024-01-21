/* interpolation for single-variable data in 2D+1D */
typedef int IDL_LONG;
#include "tri.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"

#define debug 0

float interpolate_tri2d_plus_1d(float xx, float yy, float zz, float *data, int data_in_center){

  float return_value=tri_missing;
  float w0,w1,interpolation_weights[3];
  int ix,iy,iz,iv_0,iv_1,i_tri,i_tri_old,ivert;
  
  hunt(tri_z,tri_nz,zz,&iz);
  if (zz == tri_z[tri_nz-1]){iz = tri_nz-2;}
  if (debug > 0){
    fprintf(stderr,"zz: %f nz: %li iz: %li Z[iz]: %f z: %g dz: %g\n",
	    zz,tri_nz,iz,tri_z[iz],zz, tri_z[iz+1] - tri_z[iz]);
      
      fprintf(stderr,"w0: %g\n",(zz - tri_z[iz])/(tri_z[iz+1] - tri_z[iz]));
  }
  
  if (iz >= 0 && iz <= (tri_nz-2) ){
    w0 = (tri_z[iz+1] - zz)/(tri_z[iz+1] - tri_z[iz]);
    w1 = 1. - w0;
    i_tri_old = -1; // initiate search using overlay grid
    i_tri = find_tri(xx,yy,i_tri_old,interpolation_weights,debug);

    if (debug > 0){
      fprintf(stderr,"data_in_center: %li i_tri: %li\n",data_in_center,i_tri);
    }
    
    if (i_tri >= 0){
      long int idum=0;
      return_value = 0.;
      if (data_in_center){
	long int ix,iy;
	/* zero-order interpolation - return value at cell center for triangle */
	ix = tri_ij[i_tri*2];
	iy = tri_ij[i_tri*2+1];
	iv_0 = ix+(tri_nx-1)*(iy+(tri_ny-1)*iz);
	iv_1 = iv_0;
	if (debug > 0){
	  fprintf(stderr,"I_tri: %li IX:  %li  IY: %li  IZ: %li IV: %li data: %f\n",
		  i_tri,ix,iy,iz,iv_0,data[iv_0]);
	}
	return_value  = data[iv_0];
      } else {
	/* first-order interpolation - data values at each vertex positon interpolated to position in triangle */
	for (ivert=0;ivert<3;ivert++){
	  long int ix,iy;
	  ix = (tri_vertices[ivert+3*i_tri] % tri_nx  );
	  iy = (tri_vertices[ivert+3*i_tri] % (tri_nx*tri_ny) ) / tri_nx;	
	  iv_0 = ix+tri_nx*(iy+tri_ny*iz);
	  iv_1 = ix+tri_nx*(iy+tri_ny*(iz+1));
	  if (debug > 0){
	    fprintf(stderr,"I_tri: %li IVERT: %li VERTEX_Index: %li IX:  %li  IY: %li  IZ: %li IV4: %li %lii Weight: %f\n",i_tri,ivert,tri_vertices[ivert+4*i_tri],ix,iy,iz,iv_0,iv_1,interpolation_weights[ivert]);
	    fprintf(stderr,"z weights: %g %g\n",w0,w1);
	  }
	  return_value +=
	    ( + data[iv_0]*w0
	      + data[iv_1]*w1)
	    * interpolation_weights[ivert]; 
	}
      }
    } else {
      return_value = NAN;
    }
  }

  return(return_value);
    
}
