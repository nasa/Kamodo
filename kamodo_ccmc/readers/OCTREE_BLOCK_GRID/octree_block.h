// This defines the octree_block structure that is used by setup_octree.c
// find_octree_block.c
// An array of these is initialized in IDL prodecures read_batsrus.pro 
// and read_gumics_tec.pro.
// Author: Lutz Rastaetter, 2007/12/12
#ifndef __octree_block__
#define __octree_block__
typedef struct {
    IDL_LONG refinement_level;
    IDL_LONG parent_ID;
    IDL_LONG child_count;
    IDL_LONG child_IDs[8];
    float XMIN,XMAX,XCenter,YMIN,YMAX,YCenter,ZMIN,ZMAX,ZCenter;
} octree_block;
#endif
