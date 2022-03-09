    #include "fl.h"
    #include "octree_block.h"
    int octree_blocklist_init(octree_block *octree_blocklist,
                                float *x, float *y, float *z, int N, int postions_centered)
//                              float *XMIN, float *XMAX, float *XCenter, 
//                              float *YMIN, float *YMAX, float *YCenter, 
//                              float *ZMIN, float *ZMAX, float *ZCenter, 
//                              int arrlen)
    {
      int ix=1,iy=1,iz=1;
      NX=-1;
      NY=-1;
      NZ=-1;
      N_blks=-1;
      while (NX == -1 && ix < N){
	if (x[0] == x[ix]){
	  NX=ix;
	}
	ix++;
      }
      if (NX == -1) return(-1);   // x are all unique
      while (NY == -1 && iy < N){
	if (y[0] == y[NX*ix]){
	  NY=iy; // found repetition 
	}
	iy++;
      }
      if (NY == -1) return(-2);   // y are all unique
      long* where_z_eq_z0=(long*)malloc(N*sizeof(long));
      long N_z_eq_z0=0;
      for (iz=0;iz<N;iz++){
	if (z[0] == z[iz]){
	  where_z_eq_z0[N_z_eq_z0]=iz;
	  N_z_eq_z0++;
	}
      }
      if (N_z_eq_z0==0){
	NZ=1;
      } else {
	NZ=N;
	for (iz=1;iz<N_z_eq_z0;iz++){
	  if ( (where_z_eq_z0[iz]-where_z_eq_z0[iz-1]) < NZ){
	    NZ=where_z_eq_z0[iz]-where_z_eq_z0[iz-1];
	  }
	}
      }
      N_blks = N/(NX*NY*NZ);
      if ((N_blks*NX*NY*NZ) != N){
	fprintf(stderr,'Error: number N not multiple of block size! N=%i NX=%i NY=%i NZ=%i',N,NX,NY,NZ);
	return(-4);
      }
      xr_blk=(float*)malloc(2*N_blks);
      yr_blk=(float*)malloc(2*N_blks);
      zr_blk=(float*)malloc(2*N_blks);
      int NVV=NX*NY*NZ;
      int NV=NX*NY;
      for (int ib=0;ib<N_blks;ib++){
	float dx_blk=x[NVV*ib+1]-x[NX*ib];
	xr_blk[2*ib]=x[NVV*ib]-positions_in_cell_center*0.5*dx_blk-TOLERANCE;
	xr_blk[2*ib+1]=x[NVV*ib+NX-1]+positions_in_cell_center*0.5*dx_blk-TOLERANCE;
      }
      for (ib=0;ib<N_blks;ib++){
	float dy_blk=y[NVV*ib+NX]-y[NVV*ib];
	yr_blk[2*ib]=y[NVV*ib]-positions_in_cell_center*0.5*dy_blk-TOLERANCE;
	yr_blk[2*ib+1]=y[NVV*ib+NV-1]+positions_in_cell_center*0.5*dx_blk-TOLERANCE;
      }
      for (int ib=0;ib<N_blks;ib++){
	float dx_blk=z[NVV*ib+NV]-z[NVV*ib];
	zr_blk[2*ib]=z[NVV*ib]-positions_in_cell_center*0.5*dx_blk-TOLERANCE;
	zr_blk[2*ib+1]=z[NVV*ib+NVV-1]+positions_in_cell_center*0.5*dx_blk-TOLERANCE;
      }
      dxmin=1e10;
      for (ix=0;ix<N_blk;ix++){
	float dx;
	dx=xr_blk[1+2*ix]-xr_blk[2*ix];
      }
    

       octree_block block;
       for (int i=0;i<arrlen;i++){
         block=octree_blocklist[i];
         block.XMIN=XMIN[i];
         block.XMAX=XMIN[i];
         block.YMIN=YMIN[i];
         block.YMAX=YMIN[i];
         block.ZMIN=ZMIN[i];
         block.ZMAX=ZMIN[i];
         block.XCenter=XCenter[i];
         block.YCenter=YCenter[i];
         block.ZCenter=ZCenter[i];
         block.child_count=0;
         block.refinement_level=-1;
         block.parent_ID=-1;
         for (int j=0;j<8;j++){block.child_IDs[j]=-1;}
       }
    }
