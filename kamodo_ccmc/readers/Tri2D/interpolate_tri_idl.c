/* #include <find_tetra.h> */

long int ijkl(long int i, long int j, long int k, long int l, long int NI, long int NJ, long int NK, long int NL){
  return (i+NI*(j+NJ*(k+NK*(l))));
}

/*
#define DEBUG_MAIN
#define DEBUG
#define DEBUG2 
#define debug_main
*/
#include "fl.h"
#include "tri.h"

extern void hunt(float *xx, IDL_LONG n, float x, int *jlo);


/*****************************************************************/
/* return multiple variables with single call                    */
/* input: fields[nvar,nx,ny,nz,n_blks], N_var, N_var_ID          */
/*        -> var_array[NV,N_var]                                 */
/*****************************************************************/

IDL_LONG  interpolate_tri_multivar_IDL(int argc, char *argv[]) {
  /* interpolator on 2D triangluar grid + 1D in Z     */
  /* for use with call_+external() in IDL             */
  /* Modifcation History:                             */
  /*   2023/08/10 Lutz Rastaetter - i8nitial version  */
  /*                                                  */
  /* get multiple variables from data at given location(s) */
  float *xx,*yy,*zz,*dx,*dy,dxdum,dydum;
  float *fields,*value_arr,interpolation_weights[4];
  float z_,w0,w1;
  IDL_LONG NV,iv,iz,
    ivert,NVAR,*VAR_IDs,N_VAR_IDs,ivar,NXY,NXYZ,i_tri,i_tri_old,debug;

  if (argc != 27){
    printf("Interpolate_tri_multivar_IDL: Number of arguments %ld incorrect!\n",argc);
    return(1);
  }

/* array of X,Y,Z positions */  
  xx=(float*)argv[0];
  yy=(float*)argv[1];
  zz=(float*)argv[2];
/* array of returned values */  
  value_arr=(float*)argv[3];
/* length of above arrays */
  NV=((IDL_LONG*)argv[4])[0];
  VAR_IDs=(IDL_LONG*)argv[5];
/* value_arr must be at least NV*N_VAR_IDs long */
  N_VAR_IDs=((IDL_LONG*)argv[6])[0]; 
  NVAR=((IDL_LONG*)argv[7])[0];
  /* block and grid information needed */
  fields=(float*)argv[8];

  NX=((IDL_LONG*)argv[9])[0];
  NY=((IDL_LONG*)argv[10])[0];
  NZ=((IDL_LONG*)argv[11])[0];

  /* triangle arrays */
  tri_x=(float*)argv[12];
  tri_y=(float*)argv[13];
  /* array of NZ z values - not the same shape as x and y! */
  tri_z=(float*)argv[14]; 

  tri_vertices=(IDL_LONG*)argv[15];
  n_tri=((IDL_LONG*)argv[16])[0];
  /* unstructured search grid */
  /*  index_tri=(IDL_LONG*)argv[17]; */
  /*  start_index_sg=(IDL_LONG*)argv[18];
      end_index_sg=(IDL_LONG*)argv[19]; */
  tri_center_index_sg = (IDL_LONG*)argv[17];
	
  tri_nx_sg=((IDL_LONG*)argv[18])[0];
  tri_ny_sg=((IDL_LONG*)argv[19])[0];
  /*  tri_nz_sg=((IDL_LONG*)argv[22])[0]; */

  tri_xmin_sg = ((float*)argv[20])[0];
  tri_ymin_sg = ((float*)argv[21])[0];
  /*  tetra_zmin_sg=((float*)argv[25])[0]; */
  tri_dx_sg = ((float*)argv[22])[0];
  tri_dy_sg = ((float*)argv[23])[0];
  /*  tetra_dz_sg=((float*)argv[28])[0]; */
/* missing value to be used in return */
  tri_neighbors = (IDL_LONG*)argv[24];
  MISSING = ((float*)argv[25])[0];
  debug = ((IDL_LONG*)argv[26])[0];

  NXY = NX*NY;
  NXYZ = NXY*NZ;
/* multiple variables returned - position them into separate fields
   of length NV  */
/*  i_tetra_old=1086121; */
  i_tri_old = -1L;
  i_tri = -1L;
  iz = -1;
  for (iv=0;iv<NV;iv++){
    for (ivar=0;ivar<N_VAR_IDs;ivar++) {
	value_arr[ivar+iv*NVAR] = MISSING;
    }
    z_ = zz[iv];
    hunt(tri_z,NZ,z_,&iz);
    if (z_ == tri_z[NZ-1]){iz = NZ-2;}
    
    if (iz >= 0 && iz <= (NZ-2) ){
      w0 = (tri_z[iz+1] - z_)/(tri_z[iz+1] - tri_z[iz]);
      w1 = 1. - w0;
      i_tri_old = -1; // initiate search using overlay grid
      i_tri = find_tri(xx[iv],yy[iv],i_tri_old,interpolation_weights,debug);

      if (debug > 0){
	fprintf(stderr,"iz: %ld Z[iz]: %f z: %g dz: %g\n",iz,tri_z[iz],z_, tri_z[iz+1] - tri_z[iz]); 
	fprintf(stderr,"w0: %g\n",(z_ - tri_z[iz])/(tri_z[iz+1] - tri_z[iz]));
      }

      if (i_tri >= 0){
	long int idum=0;
	for (ivar=0;ivar<N_VAR_IDs;ivar++) {
	  value_arr[ivar+iv*NVAR] = 0;
	}
	for (ivert=0;ivert<3;ivert++){
	  long int ix,iy,iv4_0,iv4_1;
	  ix = (tri_vertices[ivert+3*i_tri] % NX  );
	  iy = (tri_vertices[ivert+3*i_tri] % NXY ) / NX;

	  iv4_0 = ijkl(0,ix,iy,iz,NVAR,NX,NY,NZ);
	  iv4_1 = ijkl(0,ix,iy,iz+1,NVAR,NX,NY,NZ);

	  if (debug > 0){
	    fprintf(stderr,"I_tri: %ld IVERT: %ld VERTEX_Index: %ld IX:  %ld  IY: %ld  IZ: %ld IV4: %ld %ld Weight: %f\n",i_tri,ivert,tri_vertices[ivert+4*i_tri],ix,iy,iz,iv4_0,iv4_1,interpolation_weights[ivert]);
	    fprintf(stderr,"z weights: %g %g\n",w0,w1);
	  }
	  for (ivar=0;ivar<N_VAR_IDs;ivar++) {
	    value_arr[ivar+iv*NVAR]+=
	      ( + fields[iv4_0+VAR_IDs[ivar]]*w0
		+ fields[iv4_1+VAR_IDs[ivar]]*w1)
	      * interpolation_weights[ivert]; 
	  }
	}
	i_tri_old=i_tri;
      }
    } 
#ifdef debug_main 
    printf("IV: %ld Value 0: %f\n",iv,value_arr[iv*NVAR]);
#endif
  }
  
  return(0);
}

