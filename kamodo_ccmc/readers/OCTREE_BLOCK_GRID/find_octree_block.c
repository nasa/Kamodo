#include <stdlib.h>
#include <stdio.h>
#include "fl.h"

IDL_LONG find_octree_block(float x, float y, float z, 
			   IDL_LONG old_blocknumber, IDL_LONG max_level)
{
    IDL_LONG ib,i_block,i_root;
    int debug=0;
/* Are we still in the block we have used before? */
/* If we check this block first we may dramatically increase
   performance of the field line tracing  */
/* check root blocks */
    i_block=-2; // special missing value so we can distinguish it from result of climb_octree()
    i_root=0;

    ib=old_blocknumber;
    if (old_blocknumber >= 0){ // we might get -2 or -1 here as missing values
      if (numparents_at_AMRlevel == NULL) {
        if ( (xr_blk[2*ib] <= x) && (xr_blk[2*ib+1] >= x) &&
             (yr_blk[2*ib] <= y) && (yr_blk[2*ib+1] >= y) &&
             (zr_blk[2*ib] <= z) && (zr_blk[2*ib+1] >= z) ) {
	  return(ib);
        }
      } else {
        if (octree_blocklist[ib].XMIN <= x && octree_blocklist[ib].XMAX >= x &&
	    octree_blocklist[ib].YMIN <= x && octree_blocklist[ib].YMAX >= y &&
	    octree_blocklist[ib].ZMIN <= x && octree_blocklist[ib].ZMAX >= z ) {
	  return(ib);
        }
      }
    }
    if (numparents_at_AMRlevel == NULL)
    {
/* old-style linear search through all blocks */
        for (ib=0;ib<N_blks;ib++){
            if ( (xr_blk[2*ib] <= x) && (xr_blk[2*ib+1] >= x) &&
                 (yr_blk[2*ib] <= y) && (yr_blk[2*ib+1] >= y) &&
                 (zr_blk[2*ib] <= z) && (zr_blk[2*ib+1] >= z) ) {
		return(ib);
            }
        }
/* no block found */
        return(-1);
    }

/* octree search */
/*    fprintf(stderr,"Searching block for position: X: %f Y: %f Z: %f\n",x,y,z); */
    while((i_root < numparents_at_AMRlevel[0] ) && (i_block == -2)){
      // this one segfaults in GM
      //    while((i_root < numblocks_at_AMRlevel[0] ) && (i_block == -2)){
        ib=block_at_AMRlevel[i_root]; /* block number of root block */
        if (debug && octree_blocklist != NULL) {
            fprintf(stderr,"Checking root: %ld \nXMIN: %f XMAX: %f X: %f \nYMIN: %f YMAX: %f Y: %f\nZMIN: %f ZMAX: %f Z: %f\nChildren: %i\n",
                    ib,octree_blocklist[ib].XMIN,octree_blocklist[ib].XMAX,x,
                    octree_blocklist[ib].YMIN,octree_blocklist[ib].YMAX,y,
                    octree_blocklist[ib].ZMIN,octree_blocklist[ib].ZMAX,z,
		    octree_blocklist[ib].child_count
                    );
        }
        if ( (octree_blocklist[ib].XMIN <= x) &&
             (octree_blocklist[ib].XMAX >= x) &&
             (octree_blocklist[ib].YMIN <= y) &&
             (octree_blocklist[ib].YMAX >= y) &&
             (octree_blocklist[ib].ZMIN <= z) &&
             (octree_blocklist[ib].ZMAX >= z) ) {
	  i_block=climb_octree(ib,x,y,z,max_level-1);
        } else {
                /* next root block */
            i_root++;
        }
    }
    return(i_block);
}

IDL_LONG climb_octree(IDL_LONG root, float x, float y, float z,
                      IDL_LONG max_level){
    int ix,iy,iz;
    int debug=0;
    IDL_LONG child_id;
    
    // if (octree_blocklist[root].child_count != 8) {return(root);}
    // new: provide for termination of climb at predefined refinement level 
    // (parent block search)
    // max_level<0 will allow infinite number or climbs!
    if (max_level == 0 || octree_blocklist[root].child_count == 0) {
      return(root);
    }
    
    ix=x > octree_blocklist[root].XCenter;
    iy=y > octree_blocklist[root].YCenter;
    iz=z > octree_blocklist[root].ZCenter;

    if (debug) {
        fprintf(stderr,"Climbing tree at root: %ld \nXMIN: %f XMAX: %f X: %f IX: %ld\nYMIN: %f YMAX: %f Y: %f IY: %ld\nZMIN: %f ZMAX: %f Z: %f IZ: %ld\n",
                root,
                octree_blocklist[root].XMIN,octree_blocklist[root].XMAX,x,ix,
                octree_blocklist[root].YMIN,octree_blocklist[root].YMAX,y,iy,
                octree_blocklist[root].ZMIN,octree_blocklist[root].ZMAX,z,iz
                );
    }
    // recursion?
    child_id=octree_blocklist[root].child_IDs[ix+2*iy+4*iz];
    // account for numerical inaccuracies: consider a point within epsilon=5e-5 of the center as being covered by a nearby child block
    // this improves coverage near the body refinement in GUMICS
    if (child_id < 0){
      int flip_ix=0,flip_iy=0,flip_iz=0;
      if (fabs(octree_blocklist[root].XCenter-x) <= 5e-5) {
	child_id=octree_blocklist[root].child_IDs[(1-ix)+2*iy+4*iz];
	flip_ix=1; 
      }
      if (child_id < 0){
	if (fabs(octree_blocklist[root].YCenter-y) <= 5e-5) {
	  child_id=octree_blocklist[root].child_IDs[ix+2*(1-iy)+4*iz];
	  flip_iy=1; 
	}
      }
      if (child_id < 0){
	if (fabs(octree_blocklist[root].ZCenter-z) <= 5e-5) {
	  child_id=octree_blocklist[root].child_IDs[ix+2*iy+4*(1-iz)];
	  flip_iz=1; 
	}
      }
      if (child_id < 0 && flip_ix && flip_iy && flip_iz){
	child_id=octree_blocklist[root].child_IDs[(1-ix)+2*(1-iy)+4*(1-iz)];
      }
    }
    if (child_id < 0){
      // child_id=-1: end here
      return(child_id);
    } else {
      // child_id >=0: continue climb
      return(climb_octree(child_id,x,y,z,max_level-1));
    }
}

