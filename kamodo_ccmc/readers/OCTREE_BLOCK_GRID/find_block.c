#include <stdio.h>
#include "fl.h"
IDL_LONG find_block(float x, float y, float z,
                    float *xr_blk,float *yr_blk,float *zr_blk,
                    long N_blks,
                    long old_blocknumber)
{
    long int ib=0,i_block=-1;
    int debug=0;
    if (debug) {
        fprintf(stderr,"X %f Y %f Z %f\n0th block: %f %f %f %f %f %f\n",
               x,y,z,
               xr_blk[0],xr_blk[1],
               yr_blk[0],yr_blk[1],
               zr_blk[0],zr_blk[1]);
    }
    ib=old_blocknumber;
    if ( (ib < N_blks) && (ib > -1) &&
         (xr_blk[2*ib] < x) && (xr_blk[2*ib+1] > x) &&
         (yr_blk[2*ib] < y) && (yr_blk[2*ib+1] > y) &&
         (zr_blk[2*ib] < z) && (zr_blk[2*ib+1] > z) ) {
        return(ib);
    }
    ib=0;
    while((ib < N_blks) && (i_block == -1)){
        if ( (xr_blk[2*ib] < x) && (xr_blk[2*ib+1] > x) &&
             (yr_blk[2*ib] < y) && (yr_blk[2*ib+1] > y) &&
             (zr_blk[2*ib] < z) && (zr_blk[2*ib+1] > z) ) {i_block=ib;}
        ib++;
    }
    return(i_block);
}
