#include "fl.h"
#define DEBUG 0
IDL_LONG setup_parent(IDL_LONG iblock,IDL_LONG ilev,
		      //		      int p1, int p2, int p3, 
		      float XMIN, float YMIN, float ZMIN,
		      // IDL_LONG *numparents_at_AMRlevel,
		      // IDL_LONG *block_at_AMRlevel,
		      // octree_block *octree_blocklist,
		      IDL_LONG N_blks,
		      IDL_LONG N_parents_max, 
		      IDL_LONG *N_parents){

  float xc,yc,zc,dx_root,dy_root,dz_root;
  IDL_LONG root,iroot,root_tmp,level_factor,ix,iy,iz,ixyz,jlev;
  xc=octree_blocklist[iblock].XCenter;
  yc=octree_blocklist[iblock].YCenter;
  zc=octree_blocklist[iblock].ZCenter;
  if (DEBUG == 1){
    fprintf(stderr,"Setup_parent: XMIN: %f YMIN: %f ZMIN: %f N_blks: %i N_parent_max: %i N_parents: %p numblocks_at_AMRlevel: %p\n",XMIN,YMIN,ZMIN,N_blks,N_parents_max,N_parents,numblocks_at_AMRlevel);    
  }
  // 2014/02/07 we can't do this any more since root blocks may not be sorted
  // now we need to loop over all knwon root blocks
  root=-1;
  iroot=0;
  if (DEBUG == 1){
    fprintf(stderr,"Setup Parent: checking %i root blocks\n",numblocks_at_AMRlevel[0]);
  }
  while (root < 0 && iroot < numblocks_at_AMRlevel[0]){
    root_tmp=block_at_AMRlevel[iroot];
    if (   octree_blocklist[root_tmp].XMIN < xc &&  octree_blocklist[root_tmp].XMAX > xc
	&& octree_blocklist[root_tmp].YMIN < yc &&  octree_blocklist[root_tmp].YMAX > yc
	&& octree_blocklist[root_tmp].ZMIN < zc &&  octree_blocklist[root_tmp].ZMAX > zc){
      root=root_tmp;
    }
    iroot++;
  }
  if (DEBUG == 1){
    fprintf(stderr,"Setup_parent: Root block: %i XMIN: %f XMAX: %f YMIN: %f YMAX: %f ZMIN: %f ZMAX: %f\n",root,octree_blocklist[root].XMIN,octree_blocklist[root].XMAX,octree_blocklist[root].YMIN,octree_blocklist[root].YMAX,octree_blocklist[root].ZMIN,octree_blocklist[root].ZMAX);
  }
  jlev=0;
  if (DEBUG >= 2){
    fprintf(stderr,"Setup Parent: xc %f yc %f zc %f\np1: %f p2: %f p3: %f\n",xc,yc,zc,p1,p2,p3);
    dx_root=octree_blocklist[root].XMAX-octree_blocklist[root].XMIN;
    dy_root=octree_blocklist[root].YMAX-octree_blocklist[root].YMIN;
    dz_root=octree_blocklist[root].ZMAX-octree_blocklist[root].ZMIN;
    fprintf(stderr,"Setup Parent: Root: %i dx %f dy %f dz %f\n",root,dx_root,dy_root,dz_root);
    //    fprintf(stderr,"Setup Parent: ix %i iy %i iz %i ixyz: %i\n",ix,iy,iz,ixyz);
  }
  while (jlev < ilev){
    int ichild,ix_c,iy_c,iz_c;
    jlev++;
    ichild=climb_octree(root,xc,yc,zc,1); // only go one level
    if (DEBUG >= 1){
      fprintf(stderr,"climb_octree(root=%i,xc=%f,yc=%f,zc=%f,1)\n",
	      root,xc,yc,zc);
      fprintf(stderr,"jlev: %i child: %i\n",jlev,ichild);
    }
    if ( ichild < 0 || ichild == root ){
      // create child block here and add to parent we came from
      // the child block IDs will NOT be consecutive for a given parent
      ichild=N_blks+(*N_parents);
      if (DEBUG >= 2){
	fprintf(stderr,"Setup_Parent: jlev: %i N_blks %i N_parents %i ",
		jlev,N_blks,(*N_parents));
      }      
      // test available memory allocation
      if ( (*N_parents) >= N_parents_max){
	fprintf(stderr,"setup_parent: running out of memory - returning");
	return(-2);
      }
      block_at_AMRlevel[jlev*N_blks+numparents_at_AMRlevel[jlev]]=ichild;
      numparents_at_AMRlevel[jlev]=numparents_at_AMRlevel[jlev]+1;
      N_parents[0]=N_parents[0]+1;
      ix_c=xc>octree_blocklist[root].XCenter;
      iy_c=yc>octree_blocklist[root].YCenter;
      iz_c=zc>octree_blocklist[root].ZCenter;
      octree_blocklist[root].child_IDs[ix_c+2*iy_c+4*iz_c]=ichild;
      octree_blocklist[root].child_count++;
      octree_blocklist[ichild].refinement_level=jlev;
      octree_blocklist[ichild].parent_ID=root;
      if (ix_c == 0) {
	octree_blocklist[ichild].XMIN=octree_blocklist[root].XMIN;
	octree_blocklist[ichild].XMAX=octree_blocklist[root].XCenter;
      } else {
	octree_blocklist[ichild].XMIN=octree_blocklist[root].XCenter;
	octree_blocklist[ichild].XMAX=octree_blocklist[root].XMAX;
      }
      if (iy_c == 0) {
	octree_blocklist[ichild].YMIN=octree_blocklist[root].YMIN;
	octree_blocklist[ichild].YMAX=octree_blocklist[root].YCenter;
      } else {
	octree_blocklist[ichild].YMIN=octree_blocklist[root].YCenter;
	octree_blocklist[ichild].YMAX=octree_blocklist[root].YMAX;
      }
      if (iz_c == 0) {
	octree_blocklist[ichild].ZMIN=octree_blocklist[root].ZMIN;
	octree_blocklist[ichild].ZMAX=octree_blocklist[root].ZCenter;
      } else {
	octree_blocklist[ichild].ZMIN=octree_blocklist[root].ZCenter;
	octree_blocklist[ichild].ZMAX=octree_blocklist[root].ZMAX;
      }
      octree_blocklist[ichild].XCenter=
	(octree_blocklist[ichild].XMIN+octree_blocklist[ichild].XMAX)/2.;
      octree_blocklist[ichild].YCenter=
	(octree_blocklist[ichild].YMIN+octree_blocklist[ichild].YMAX)/2.;
      octree_blocklist[ichild].ZCenter=
	(octree_blocklist[ichild].ZMIN+octree_blocklist[ichild].ZMAX)/2.;
      if (DEBUG >= 1){
	fprintf(stderr,"Setup_Parent: ilev: %i Root block: %i\nNew child block: jlev: %i ID: %i\n",
		ilev,root,jlev,ichild);
      }
    } 
    root=ichild;
  }
  if (DEBUG >= 3){
    fprintf(stderr,"Exiting Setup_Parent(): Returning root: %i\n",root);
  }
  return(root);
  //     return(-1);
}
