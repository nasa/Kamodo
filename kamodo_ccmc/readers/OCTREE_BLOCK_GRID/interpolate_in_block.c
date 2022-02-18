#include "fl.h"
float interpolate_in_block(float x, float y, float z, 
                           float *x_blk, float *y_blk, float *z_blk,
                           float *data,
                           int NX,int NY,int NZ)
{
    float value,dx_blk,dy_blk,dz_blk,m_x,m_y,m_z;
    int ix,iy,iz;
    long int NV_blk;

    NV_blk=NX*NY;
    find_in_block(x,y,z,&ix,&iy,&iz,
                  x_blk,y_blk,z_blk,
                  NX,NY,NZ);
    
    dx_blk = x_blk[ix+1]-x_blk[ix];
    dy_blk = y_blk[iy+1]-y_blk[iy];
    dz_blk = z_blk[iz+1]-z_blk[iz];
    m_x = (x-x_blk[ix])/dx_blk;
    m_y = (y-y_blk[iy])/dy_blk;
    m_z = (z-z_blk[iz])/dz_blk;
/*    printf("Interpolate_in_block: M: %f %f %f\n",m_x,m_y,m_z);  */
    value=(1-m_z)*(
        (1-m_y)*(
            (1-m_x)*data[ix  +iy*NX+iz*NV_blk]
            +  m_x *data[ix+1+iy*NX+iz*NV_blk]
            )
        + m_y*(
            + (1-m_x)*data[ix  +(iy+1)*NX+iz*NV_blk]
            +    m_x *data[ix+1+(iy+1)*NX+iz*NV_blk]
            )
        )
        + m_z*(
            (1-m_y)*(
                +(1-m_x)*data[ix  +iy*NX+(iz+1)*NV_blk]
                +   m_x *data[ix+1+iy*NX+(iz+1)*NV_blk]
                )
            + m_y*(
                +(1-m_x)*data[ix  +(iy+1)*NX+(iz+1)*NV_blk]
                +   m_x *data[ix+1+(iy+1)*NX+(iz+1)*NV_blk]
                )
            );
/*    printf("X: %f Y: %f Z: %f V: %f\n",x,y,z,value); */
    
    return(value);
}



