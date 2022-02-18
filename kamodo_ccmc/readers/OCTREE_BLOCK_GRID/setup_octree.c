/**********************************************/
/* using the octree structure of BATSRUS data */
/* see br_export.c                            */
/* fl.h includes standard headers and         */
/* IDL_export.h that is linked to the         */
/* appropriate IDL version's expot.h file.    */
/* Author: lutz Rastaetter NASA/GSFC          */
/**********************************************/
/* Modification History:                      */
/*                                            */
/* 2012/04/13 LR - addressed case when there are blocks at AMR level 0 */
/*                 that do not fully cover the simulation box. */
/**********************************************/

#include "fl.h"

//#define DEBUG
#define no_idl
#ifndef no_idl
// define setup_octree_IDL with untyped inputs
void setup_octree_IDL(int argc, char *argv[])
{
        /* arguments:
         long N_blks
         float[0:2*N_blk-1] xr_blk,yr_blk,zr_blk
        */
//  float p1,p2,p3,aspect_xy,aspect_yz,aspect_xz,
    float aspect_xy,aspect_yz,aspect_xz,
      XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,RMIN,level_factor=1.,dx_min,dy_min,dz_min;

    IDL_LONG MAX_AMRLEVELS,*success;
    float *xr_blk,*yr_blk,*zr_blk,*box_range; // ,TOLERANCE=0.001; already in fl.h
    //    octree_block *octree_blocklist;
    IDL_LONG N_octree_blocks;
    //  *block_at_AMRlevel,
    //  *numparents_at_AMRlevel,
    //  *numblocks_at_AMRlevel,
    // N_blks;
    IDL_LONG ilev,ix,iy,iz,i,ib,ib2,iblk, N_parents,N_parents_max;
    if (argc != 12) {
        printf("Setup Octree_IDL: incorrect number of arguments: %ld\n",argc );
    return;}

    N_blks=*((IDL_LONG*)argv[0]);
    xr_blk=(float*)argv[1];
    yr_blk=(float*)argv[2];
    zr_blk=(float*)argv[3];
    MAX_AMRLEVELS=((IDL_LONG*)argv[4])[0];
    box_range=(float*)argv[5];
    octree_blocklist=(octree_block*)argv[6];
    N_octree_blocks=*((IDL_LONG*)argv[7]);
    numparents_at_AMRlevel=(IDL_LONG*)argv[8];
    numblocks_at_AMRlevel=(IDL_LONG*)argv[9];
    block_at_AMRlevel=(IDL_LONG*)argv[10];
    success=(IDL_LONG*)argv[11];
    
    XMIN=box_range[0];
    XMAX=box_range[1];
    YMIN=box_range[2];
    YMAX=box_range[3];
    ZMIN=box_range[4];
    ZMAX=box_range[5];
    RMIN=box_range[6]; // body_boundary_radius, if defined
    
/*    numblocks_at_AMRlevel=(IDL_LONG*)malloc(sizeof(IDL_LONG)*(1+MAX_AMRLEVELS));
 */
#ifdef DEBUG
    fprintf(stderr,"Inputs: Box: %f %f %f %f %f %f Blocks: %ld %ld\n",
           XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,N_blks,MAX_AMRLEVELS);    
#endif
    success[0]=-1;

#else
// if defined no_idl
// define setup_octree with list of inputs and types
IDL_LONG   setup_octree(IDL_LONG N_blks_in,
		    float *xr_blk,float *yr_blk,float *zr_blk,
		    IDL_LONG MAX_AMRLEVELS, float *box_range,
		    octree_block *octree_blocklist_in, 
		    IDL_LONG N_octree_blocks,
		    IDL_LONG *numparents_at_AMRlevel_in,
		    IDL_LONG *block_at_AMRlevel_in //,
			//		    IDL_LONG *success
		    ){
    IDL_LONG success[1];
    float p1,p2,p3,aspect_xy,aspect_yz,aspect_xz,
      XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,RMIN,level_factor=1.;

    //    IDL_LONG         *numblocks_at_AMRlevel=NULL;
    
    IDL_LONG ilev,ix,iy,iz,i,ib,ib2,iblk,N_parents,N_parents_max;

    // copy inputs to common blcok variables
    
    octree_blocklist=octree_blocklist_in;
    numparents_at_AMRlevel=numparents_at_AMRlevel_in;
    block_at_AMRlevel=block_at_AMRlevel_in;
    N_blks=N_blks_in;

    if (numblocks_at_AMRlevel == NULL) {
      numblocks_at_AMRlevel=(IDL_LONG*)malloc(sizeof(IDL_LONG)*(1+MAX_AMRLEVELS));
    }
    XMIN=box_range[0];
    XMAX=box_range[1];
    YMIN=box_range[2];
    YMAX=box_range[3];
    ZMIN=box_range[4];
    ZMAX=box_range[5];
    RMIN=box_range[6];
    // end: defined no_idl
#endif
    
/* older data without "th", "p1", "p2", "p3" */
/* This alghorithm fits the lowest number of cubic "root" blocks
   into the simulation box */
/* aspect ratio of first block (all blocks have the same) */
    aspect_xz=(xr_blk[1]-xr_blk[0])/(zr_blk[1]-zr_blk[0]);
    aspect_yz=(yr_blk[1]-yr_blk[0])/(zr_blk[1]-zr_blk[0]);
    aspect_xy=(xr_blk[1]-xr_blk[0])/(yr_blk[1]-yr_blk[0]);
        
    //    p1=max((XMAX-XMIN)/(aspect_xy*(YMAX-YMIN)),
    //       (XMAX-XMIN)/(aspect_xz*(ZMAX-ZMIN)) );
    //p2=max((YMAX-YMIN)/(aspect_yz*(ZMAX-ZMIN)),1);
    //p3=max((ZMAX-ZMIN)*aspect_yz/(YMAX-YMIN),1);
    p1=100000;
    p2=100000;
    p3=100000;
    for (ib=0;ib<N_blks;ib++){
      float  p1_tmp,p2_tmp,p3_tmp;
      p1_tmp=(XMAX-XMIN)/(xr_blk[2*ib+1]-xr_blk[2*ib]);
      p2_tmp=(YMAX-YMIN)/(yr_blk[2*ib+1]-yr_blk[2*ib]);
      p3_tmp=(ZMAX-ZMIN)/(zr_blk[2*ib+1]-zr_blk[2*ib]);
      if (p1_tmp < p1) {p1=p1_tmp;}
      if (p2_tmp < p2) {p2=p2_tmp;}
      if (p3_tmp < p3) {p3=p3_tmp;}
    }
    printf("P1: %f P2: %f P3: %f DP: %f %f %f\n",
	   p1,p2,p3,
	   fabs(p1-floor(p1+0.5)),
	   fabs(p2-floor(p2+0.5)),
	   fabs(p3-floor(p3+0.5))
	   );
    
    if (fabs(p1-floor(p1+0.5)) > TOLERANCE) { 
        level_factor=floor(p1+0.5)/p1; // 1./min( ceil(p1)-p1 , p1-floor(p1) );
        p1*=level_factor;
        p2*=level_factor;
        p3*=level_factor;
    }
    if (fabs(p2-floor(p2+0.5)) > TOLERANCE) {
        level_factor=floor(p2+0.5)/p2; // 1./min( ceil(p2)-p2 , p2-floor(p2) );
        p1*=level_factor;
        p2*=level_factor;
        p3*=level_factor;
    }
    if (fabs(p3-floor(p3+0.5)) > TOLERANCE) {
        level_factor=floor(p3+0.5)/p3; // 1./min( ceil(p3)-p3 , p3-floor(p3) );
        p1*=level_factor;
        p2*=level_factor;
        p3*=level_factor;
    }
    p1=floor(p1+0.5);
    p2=floor(p2+0.5);
    p3=floor(p3+0.5);
    // reduce p1,p2,p3 by powers of 2
    while (    ( (p1/2.-floor(p1/2.+0.25)) < 0.4) 
	    && ( (p2/2.-floor(p2/2.+0.25)) < 0.4) 
	    && ( (p3/2.-floor(p3/2.+0.25)) < 0.4) ){
            p1/=2.;
            p2/=2.;
            p3/=2.;
    }

    printf("root outlay: %f %f %f Factor: %f Aspect: %f %f %f\n",
	    p1,p2,p3,level_factor,
	    aspect_xy,aspect_yz,aspect_xz);
    printf("box:: %f %f %f %f %f %f\n",XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX);
    
    // make integer numbers
    p1=floor(p1+0.1);
    p2=floor(p2+0.1);
    p3=floor(p3+0.1);

#ifdef DEBUG
    fprintf(stderr,"Root block outlay: %f %f %f\n", p1,p2,p3);
    fprintf(stderr,"Max AMR level: %i\n",MAX_AMRLEVELS);
#endif    
/*    MAX_AMRLEVELS=floor(log((XMAX-XMIN)/(p1*dxmin_blk) )/log(2.)+1.5);
    fprintf(stderr,"Log((XMAX-XMIN)/(p1*dxmin_blk)): %f %ld\n",
            log((XMAX-XMIN)/(p1*dxmin_blk) )/log(2.),MAX_AMRLEVELS);
*/
//#ifndef no_idl
    for(ilev=0;ilev<MAX_AMRLEVELS;ilev++){
        numblocks_at_AMRlevel[ilev]=0;
        numparents_at_AMRlevel[ilev]=0;
    }
//#endif
    
#ifdef DEBUG
    fprintf(stderr,"Box size: %f %f  %f %f  %f %f\n",XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX);
    fprintf(stderr,"Initial root block outlay: %f %f %f N_blks:  %ld AMR levels: %ld\n",p1,p2,p3,N_blks,MAX_AMRLEVELS);
    fprintf(stderr,"Number of leaf blocks: %i  Number of all octree blocks: %i\n",N_blks,N_octree_blocks);
#endif    
/* sort blocks into AMR refinement levels */
    for (iblk=N_blks;iblk<N_octree_blocks;iblk++){
      int ic;
      for (ic=0;ic<8;ic++) {
	octree_blocklist[iblk].child_IDs[ic]=-1;
      }
      octree_blocklist[iblk].child_count=0;
    }
#ifdef DEBUG
    fprintf(stderr,"Got here: N_blks=%i\n",N_blks);
#endif
      
    for (iblk=0;iblk<N_blks;iblk++){
        ilev=floor(log((XMAX-XMIN)/(p1*(xr_blk[iblk*2+1]-xr_blk[iblk*2])))
                   /log(2.)+0.5);
	//#ifdef DEBUG
	//	fprintf(stderr,"iblk: %i ilev: %i max level: %i numblocks_at_AMRlevel: %p\n",iblk,ilev,MAX_AMRLEVELS,numblocks_at_AMRlevel);
	//#endif
        if (ilev > (MAX_AMRLEVELS-1)){
	  printf("AMR level %i too high (>%i)\n",ilev,MAX_AMRLEVELS-1);
#ifndef no_idl
            success[0]=-1;
            return;
#else
	    return(-1);
#endif
        }

	block_at_AMRlevel[ilev*N_blks+numblocks_at_AMRlevel[ilev]]=iblk;
        numblocks_at_AMRlevel[ilev]++;
    }
#ifdef DEBUG
    for (ilev=0;ilev<MAX_AMRLEVELS;ilev++){
      fprintf(stderr,"Blocks ar AMR level %i: =%i\n",ilev,numblocks_at_AMRlevel[ilev]);
    }
#endif
    for (iblk=0;iblk<N_blks;iblk++){    
	octree_blocklist[iblk].refinement_level=ilev;
        octree_blocklist[iblk].XMIN=xr_blk[2*iblk];
        octree_blocklist[iblk].XMAX=xr_blk[2*iblk+1];
        octree_blocklist[iblk].XCenter=0.5*(xr_blk[2*iblk]+xr_blk[2*iblk+1]);
        octree_blocklist[iblk].YMIN=yr_blk[2*iblk];
        octree_blocklist[iblk].YMAX=yr_blk[2*iblk+1];
        octree_blocklist[iblk].YCenter=0.5*(yr_blk[2*iblk]+yr_blk[2*iblk+1]);
        octree_blocklist[iblk].ZMIN=zr_blk[2*iblk];
        octree_blocklist[iblk].ZMAX=zr_blk[2*iblk+1];        
        octree_blocklist[iblk].ZCenter=0.5*(zr_blk[2*iblk]+zr_blk[2*iblk+1]);
	//        octree_blocklist[iblk].ID=iblk;
        octree_blocklist[iblk].parent_ID=-1;
	//        octree_blocklist[iblk].AMR_level=ilev;
        octree_blocklist[iblk].child_IDs[0]=-1;
        octree_blocklist[iblk].child_IDs[1]=-1;
        octree_blocklist[iblk].child_IDs[2]=-1;
        octree_blocklist[iblk].child_IDs[3]=-1;
        octree_blocklist[iblk].child_IDs[4]=-1;
        octree_blocklist[iblk].child_IDs[5]=-1;
        octree_blocklist[iblk].child_IDs[6]=-1;
        octree_blocklist[iblk].child_IDs[7]=-1;
        octree_blocklist[iblk].child_count=0;
    }
#ifdef DEBUG
    fprintf(stderr,"Got here (2)\n");
#endif
    for (iblk=N_blks;iblk< N_octree_blocks;iblk++){
        octree_blocklist[iblk].parent_ID=-1;
	octree_blocklist[iblk].refinement_level=-1;
        octree_blocklist[iblk].child_IDs[0]=-1;
        octree_blocklist[iblk].child_IDs[1]=-1;
        octree_blocklist[iblk].child_IDs[2]=-1;
        octree_blocklist[iblk].child_IDs[3]=-1;
        octree_blocklist[iblk].child_IDs[4]=-1;
        octree_blocklist[iblk].child_IDs[5]=-1;
        octree_blocklist[iblk].child_IDs[6]=-1;
        octree_blocklist[iblk].child_IDs[7]=-1;
        octree_blocklist[iblk].child_count=0;        
    }
    
    for (ilev=0;ilev<=MAX_AMRLEVELS-1;ilev++){
        fprintf(stderr,"AMR level: %ld number of blocks: %ld\n",
		ilev,numblocks_at_AMRlevel[ilev]);
    }
    
    // new algorithm - set up blocks from the coarsest level and 
    // employ fast search to get to the parent level for each existing block

    N_parents=0;
    // allocated number of octree blocks (handed in from IDL)
    N_parents_max=N_octree_blocks-N_blks;
#ifdef DEBUG
    fprintf(stderr,"N_octree_blocks: %i N_blks: %i N_parents_max: %i\n", N_octree_blocks,N_blks,N_parents_max);
#endif
    // need to check coarsest level (ilev=0) for existing blocks and add new ones
    int ilev_start=0; 

    if (numblocks_at_AMRlevel[0] > 0 && numblocks_at_AMRlevel[0] < (p1*p2*p3) ){
	float dx,dy,dz; // block size at level ilev
	int N_parents_old,num_new_blocks,N_ix,N_iy,N_iz,ix,iy,iz; 
	float XC,YC,ZC;
	int ib,iblk,nblk0;
	ilev=0;
	dx=(XMAX-XMIN)/p1;
	dy=(YMAX-YMIN)/p2;
	dz=(ZMAX-ZMIN)/p3;
	iblk=N_blks;
	nblk0=numblocks_at_AMRlevel[0];
	for (iz=0;iz<p3;iz++){
	  ZC=ZMIN+(iz+0.5)*dz;
	  for (iy=0;iy<p2;iy++){
	    YC=YMIN+(iy+0.5)*dy;
	    for (ix=0;ix<p1;ix++){
	      int block_exists=0;
	      XC=XMIN+(ix+0.5)*dx;
	      // test existing blocks;
	      for (ib=0;ib<nblk0;ib++){
		int test_block;
		test_block=block_at_AMRlevel[ib];
		//		fprintf(stderr,"Test block: %i\n",test_block);
		//fprintf(stderr,"testing X: %f against %f %f\n",XC,octree_blocklist[test_block].XMIN,octree_blocklist[test_block].XMAX);
		//fprintf(stderr,"testing Y: %f against %f %f\n",YC,octree_blocklist[test_block].YMIN,octree_blocklist[test_block].YMAX);
		//fprintf(stderr,"testing Z: %f against %f %f\n",ZC,octree_blocklist[test_block].ZMIN,octree_blocklist[test_block].ZMAX);

		if  ( (XC > octree_blocklist[test_block].XMIN && XC < octree_blocklist[test_block].XMAX )
		      && (YC > octree_blocklist[test_block].YMIN && YC < octree_blocklist[test_block].YMAX )
		      && (ZC > octree_blocklist[test_block].ZMIN && ZC < octree_blocklist[test_block].ZMAX ) )
		  {
		    block_exists=1;
		  }
	      }
	      if (block_exists == 0){
		// increment local counter (level = 0)
		fprintf(stderr,"Setting up new block %i at AMR level 0: %i\n",iblk,numblocks_at_AMRlevel[0]);
		// set new block number into block_at_AMRlevel array
		block_at_AMRlevel[numblocks_at_AMRlevel[0]]=iblk;
		
		octree_blocklist[iblk].refinement_level=ilev;
		octree_blocklist[iblk].XMIN=XMIN+ix*dx;
		octree_blocklist[iblk].XMAX=XMIN+(ix+1)*dx;
		octree_blocklist[iblk].XCenter=XMIN+(ix+0.5)*dx;
		octree_blocklist[iblk].YMIN=YMIN+iy*dy;
		octree_blocklist[iblk].YMAX=YMIN+(iy+1)*dy;
		octree_blocklist[iblk].YCenter=YMIN+(iy+0.5)*dy;
		octree_blocklist[iblk].ZMIN=ZMIN+iz*dz;
		octree_blocklist[iblk].ZMAX=ZMIN+(iz+1)*dz;
		octree_blocklist[iblk].ZCenter=ZMIN+(iz+0.5)*dz;
		// increment global counters		  
		iblk++;
		N_parents++;
		numblocks_at_AMRlevel[0]=numblocks_at_AMRlevel[0]+1;
	      }
	    }
	  }
	}
	// we need to set these for any interpolation to succeed;
	numblocks_at_AMRlevel[0]=p1*p2*p3;
	numparents_at_AMRlevel[0]=p1*p2*p3; // this is needed to all lookups 
	fprintf(stderr,"Added %i blocks to coarsest refinement level\n",N_parents);
	ilev_start=1; 
    }

    for (ilev=ilev_start;ilev<=MAX_AMRLEVELS-1;ilev++){
    //    for (ilev=0;ilev<=min(3,MAX_AMRLEVELS-1);ilev++){
      level_factor=pow(2.,ilev);
      if (numblocks_at_AMRlevel[ilev] <= 0){
	float dx,dy,dz; // block size at level ilev
	int num_new_blocks,N_ix,N_iy,N_iz,ix,iy,iz; 
	dx=(XMAX-XMIN)/(level_factor*p1);
	dy=(YMAX-YMIN)/(level_factor*p2);
	dz=(ZMAX-ZMIN)/(level_factor*p3);
	// setup virtual parent blocks
	N_ix=level_factor*p1;
	N_iy=level_factor*p2;
	N_iz=level_factor*p3;

	//	fprintf(stderr,"New Algorithm: setting up %i new blocks \n with (dx,dy,dz)=(%f,%f,%f)\n and refinement level: %i\n",
	//	N_ix*N_iy*N_iz,dx,dy,dz,ilev);

	//	if (ilev == 0) { // establish root blocks only
	for (iz=0;iz<N_iz;iz++){
	  int iz_c;
	  iz_c=iz % 2;
	  for (iy=0;iy<N_iy;iy++){
	    int iy_c;
	    iy_c=iy % 2;
	    for (ix=0;ix<N_ix;ix++){
	      int ix_c,iv_c,parent_block,ixyz;
	      ix_c=ix % 2;
	      iv_c=ix_c+2*iy_c+4*iz_c;
	      ixyz=ix+N_ix*(iy+N_iy*iz);
	      iblk=N_blks+N_parents+ixyz;
	      octree_blocklist[iblk].refinement_level=ilev;
	      octree_blocklist[iblk].XMIN=XMIN+ix*dx;
	      octree_blocklist[iblk].XMAX=XMIN+(ix+1)*dx;
	      octree_blocklist[iblk].XCenter=XMIN+(ix+0.5)*dx;
	      octree_blocklist[iblk].YMIN=YMIN+iy*dy;
	      octree_blocklist[iblk].YMAX=YMIN+(iy+1)*dy;
	      octree_blocklist[iblk].YCenter=YMIN+(iy+0.5)*dy;
	      octree_blocklist[iblk].ZMIN=ZMIN+iz*dz;
	      octree_blocklist[iblk].ZMAX=ZMIN+(iz+1)*dz;
	      octree_blocklist[iblk].ZCenter=ZMIN+(iz+0.5)*dz;
	      if (ilev > 0){
		int iparent,ip;
		iparent=((ix/2)+(N_ix/2)*((iy/2)+(N_iy/2)*(iz/2)));
		// update child-parent connections
		ip=block_at_AMRlevel[(ilev-1)*N_blks+iparent];
		if (octree_blocklist[ip].refinement_level != (ilev-1)){
		  fprintf(stderr,"Setup Octree: ix %i N_ix: %i iy: %i N_iy: %i iz: %i N_iz: %i parent %i in level %i.\n",
			  ix,N_ix,iy,N_iy,iz,N_iz,iparent,ilev-1);
		  fprintf(stderr,"Setup Octree: parent block %i in wrong refinement level %i!=%i (ilev).\n",ip,octree_blocklist[ip].refinement_level,ilev);
#ifndef no_idl
		  return;
#else
		  return(-1);
#endif
		}
		if (octree_blocklist[iblk].parent_ID < 0) {
		  octree_blocklist[iblk].parent_ID=ip;
		} 
		if (octree_blocklist[ip].child_IDs[iv_c] == -1){
		  octree_blocklist[ip].child_IDs[iv_c]=iblk;
		  octree_blocklist[ip].child_count++;
		} else {
		  fprintf(stderr,"Block ip=%ld:  Child block already assigned: iv_c=%ld ix_c=%ld iy_c=%ld iz_c=%ld\n",ip,iv_c,ix_c,iy_c,iz_c);
		}
	      }
	    }
	  }
      	}
	  // add blocks to total number of blocks 
	num_new_blocks=N_ix*N_iy*N_iz;
	numparents_at_AMRlevel[ilev]=num_new_blocks;
	fprintf(stderr,"Setup_octree: adding %i blocks to block_at_amrlevel\n %i\n",num_new_blocks,ilev);
	for (i=0;i<=num_new_blocks-1;i++){
#ifdef DEBUG
	  if (ilev == 7 && i < 10)  fprintf(stderr,"ilev: %i i: %i ilev*N_blks+i: %i N_blks+N_parents+i: %i\n",ilev,i,ilev*N_blks+i,N_blks+N_parents+i);
#endif
	  block_at_AMRlevel[ilev*N_blks+i]=N_blks+N_parents+i;
	}
	N_parents+=num_new_blocks;
	numblocks_at_AMRlevel[ilev]=num_new_blocks;
	//      } // add root blocks only
      } else {
	// run through existing blocks
	for (iblk=0;iblk<numblocks_at_AMRlevel[ilev];iblk++){
	  int ib,parent_block,ichild,ix_c,iy_c,iz_c;
	  float xc,yc,zc;
	  ib=block_at_AMRlevel[ilev*N_blks+iblk];
	  xc=octree_blocklist[ib].XCenter;
	  yc=octree_blocklist[ib].YCenter;
	  zc=octree_blocklist[ib].ZCenter;
	  //#ifdef DEBUG 
	  //fprintf(stderr,"About to setup_parent()\n");
	  //#endif
#ifdef DEBUG
	  fprintf(stderr,"Existing block loop 1\n");
#endif
	  parent_block=setup_parent(ib,ilev-1,
				    //p1,p2,p3,
				    XMIN,YMIN,ZMIN,
				    //numparents_at_AMRlevel,
				    //block_at_AMRlevel,
				    //octree_blocklist,
				    N_blks,N_parents_max,&N_parents);
#ifdef DEBUG
	  fprintf(stderr,"Existing block loop 2: iblock: %i iParent: %i\n",ib,parent_block);
#endif

	  if (parent_block < 0){
	    fprintf(stderr,"Parent block missing for block %i\n",ib);
	    // error encountered
#ifndef no_idl
	    return;
#else
	    return(-1);
#endif
	  }
// establish connections
	  ix_c=xc > octree_blocklist[parent_block].XCenter;
	  iy_c=yc > octree_blocklist[parent_block].YCenter;
	  iz_c=zc > octree_blocklist[parent_block].ZCenter;
	  ichild=ix_c+2*iy_c+4*iz_c;
	  if (octree_blocklist[parent_block].child_IDs[ichild] < 0){
	    octree_blocklist[parent_block].child_IDs[ichild]=ib;
	    octree_blocklist[parent_block].child_count++;
	  } else {
	    fprintf(stderr,"Child ID was already assigned\n");
	  }
	  octree_blocklist[ib].parent_ID=parent_block;
	}
      }
      fprintf(stderr,"Setup_Octree::new_algorithm: level: %i of %i level_factor: %f N_Parents: Total: %i Max: %i\n",ilev,MAX_AMRLEVELS-1,level_factor,N_parents,N_parents_max);
    }

    for (ilev=0;ilev<=MAX_AMRLEVELS-1;ilev++){
      fprintf(stderr,"Parents at AMR level %i: %i\n",
	      ilev,numparents_at_AMRlevel[ilev]);
    }
/* add some coverage to root blocks */
    for (ib=0;ib<numparents_at_AMRlevel[0];ib++){
        float eps = 0.00025;
        iblk=block_at_AMRlevel[ib];
        octree_blocklist[iblk].XMIN=octree_blocklist[iblk].XMIN-eps;
        octree_blocklist[iblk].XMAX=octree_blocklist[iblk].XMAX+eps;
        octree_blocklist[iblk].YMIN=octree_blocklist[iblk].YMIN-eps;
        octree_blocklist[iblk].YMAX=octree_blocklist[iblk].YMAX+eps;
        octree_blocklist[iblk].ZMIN=octree_blocklist[iblk].ZMIN-eps;
        octree_blocklist[iblk].ZMAX=octree_blocklist[iblk].ZMAX+eps;
    }

    // I don't know how to eassiloy tell success for new algorithm 
    // so we return this: 
    //        success[0]=numparents_at_AMRlevel[0] > 0;

    // test setup for position in one most hightly refined block
    // we need body_boundary_radius to do this reliably for all models 
    // this routine may operate on: SWMF/GM, SWMF/SC, SWMF/IH, GUMICS

//    fprintf(stderr,"testing whether:\n    find_octree_block(1.1*%f,0.1,0.1,-1,MAX_AMRLEVELS) < %f\n",
//  RMIN,N_blks);
    ib=find_octree_block(1.1*RMIN,0.,0.,-1,MAX_AMRLEVELS);
//    fprintf(stderr,"setup_octree: test block: %ld N_blk: %ld\n",ib,N_blks);
    success[0] = ib < N_blks;

    /* debugging output */
#ifdef DEBUG
    if (success[0] == 0){
      for (ib=N_blks;ib<N_blks+N_parents;ib++){
	if (octree_blocklist[ib].child_count != 8 ){
	  fprintf(stderr,
		  "Block: %ld #children: %ld ilev: %ld\n  XMIN: %f XMAX: %f YMIN: %f YMAX: %f ZMIN: %f ZMAX: %f\n",
		  ib,
		  octree_blocklist[ib].child_count,
		  octree_blocklist[ib].refinement_level,
		  octree_blocklist[ib].XMIN,
		  octree_blocklist[ib].XMAX,
		  octree_blocklist[ib].YMIN,
		  octree_blocklist[ib].YMAX,
		  octree_blocklist[ib].ZMIN,
		  octree_blocklist[ib].ZMAX);
	}
      }
    } 
#endif

#ifndef no_idl
    return;
#else
    return(0); // success!
#endif
    
}
