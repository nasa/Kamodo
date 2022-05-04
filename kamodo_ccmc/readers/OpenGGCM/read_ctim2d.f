!--------------------------------------------------------------------------
! Purpose:  to read 2D field from OpenGGCM *iof* file
!
! Modification history:
!     2022/03/23 Lutz Rastaetter - adjust to stricter compilation with f2py
c--------------------------------------------------------------------------
      subroutine read_ctim2d_field(filename,fieldname,nnx,nny,asciitime,
     &                             field2d)
      implicit none
      
      INTEGER     ARRAY_SIZE,CHAR_SIZE,nnx,nny,it
      PARAMETER ( CHAR_SIZE  = 120 )
      PARAMETER ( ARRAY_SIZE = 2000000 )

      character*80 fieldname,asciitime
      character*120 filename
      real*4  field2d(ARRAY_SIZE)

!f2py intent(in) :: filename
!f2py intent(in) :: fieldname
!f2py intent(out) :: nnx,nny,asciitime
!f2py intent(in,out) :: field2d

!      write (0,*) 'opening file: |',filename,'| to read |',fieldname,'|'
      
      open(unit=10,file=filename,status='old',action='read',err=222)
      call getf21(10,field2d,nnx,nny,fieldname,asciitime,it)
!      write(0,*)' got ',fieldname,' ',nnx,nny

      close(10)
      return
 222  write(0,*) "could not open file"

      end
      
c--------------------------------------------------------
      subroutine read_ctim2d(f2d,nx,ny,gx,gy,
     *   prec_e_fe_1,prec_e_e0_1,prec_e_fe_2,prec_e_e0_2,pacurr,fac_tot,
     *   pot,delphi,ppio,rrio,ttio,vdown,sigh,sigp,fac_dyn,ctiot,ctiop,
     *   tau,ctaup,ctaut,cpolp,cpolt,delbr,delbp,delbt,epio,etio,xjh)
c--------------------------------------------------------
c
c   example program to read UCLA MHD fields
c
c      common /a1/bx(1500000),by(1500000),bz(1500000),vx(1500000),vy(
c     *1500000),vz(1500000),rr(1500000),pp(1500000)
c      common /a2/gx(200),gy(200),gz(200)
C	Declare a Fortran Record type that has the same form as the 
C	IDL C struct STRING. While Fortran records are not part of 
C	F77, most compilers have this option.
C
C   	Declare the string structure

C	Declare a parameter for the length of the Fortran character strings
      implicit none
      
      INTEGER     ARRAY_SIZE,CHAR_SIZE,GRID_SIZE,nx,ny,nnx,nny,it,i
      PARAMETER ( CHAR_SIZE  = 120 )
      PARAMETER ( ARRAY_SIZE = 2000000 )
      PARAMETER ( GRID_SIZE  = 400 )

      character*120 l1,l2,f2d
      real gx(GRID_SIZE),gy(GRID_SIZE),
     *     prec_e_fe_1(ARRAY_SIZE),
     *     prec_e_e0_1(ARRAY_SIZE),
     *     prec_e_fe_2(ARRAY_SIZE),
     *     prec_e_e0_2(ARRAY_SIZE),
     *     pacurr(ARRAY_SIZE),
     *     fac_tot(ARRAY_SIZE),
     *     pot(ARRAY_SIZE),
     *     delphi(ARRAY_SIZE),
     *     ppio(ARRAY_SIZE),
     *     rrio(ARRAY_SIZE),
     *     ttio(ARRAY_SIZE),
     *     vdown(ARRAY_SIZE),
     *     sigh(ARRAY_SIZE),
     *     sigp(ARRAY_SIZE),
     *     fac_dyn(ARRAY_SIZE),
     *     ctiot(ARRAY_SIZE),
     *     ctiop(ARRAY_SIZE),
     *     tau(ARRAY_SIZE),
     *     ctaup(ARRAY_SIZE),
     *     ctaut(ARRAY_SIZE),
     *     cpolp(ARRAY_SIZE),
     *     cpolt(ARRAY_SIZE),
     *     delbr(ARRAY_SIZE),
     *     delbp(ARRAY_SIZE),
     *     delbt(ARRAY_SIZE),
     *     epio(ARRAY_SIZE),
     *     etio(ARRAY_SIZE),
     *     xjh(ARRAY_SIZE)     
!f2py intent(in) :: f2d
!f2py intent(in) :: gridName
!f2py intent(out) :: nx,ny
!f2py intent(in,out) :: gx,gy
!f2py intent(in,out) :: prec_e_fe_1,prec_e_e0_1,prec_e_fe_2,prec_e_e0_2
!f2py intent(in,out) :: pacurr,fac_tot,pot,delphi,ppio,rrio,ttio
!f2py intent(in,out) :: vdown,sigh,sigp,fac_dyn,ctiot,ctiop
!f2py intent(in,out) :: tau,ctaup,ctaut,cpolp,cpolt,delbr,delbp,delbt,epio,etio,xjh
      
c      real gx(GRID_SIZE),gy(GRID_SIZE),gz(GRID_SIZE)
c
c      l1='gridx'
cc      write(*,*) l1,gx(1:10)
c     call getf11(10,gx,nnx,l1,l2,it) c      nx=nnx
c      write(0,*)' got grid x ',l1,' ',nx,'\n'
c      write(0,*)' from ',gx(1),' to ',gx(nx),'\n'
c      l1='gridy'
c      call getf11(10,gy,nny,l1,l2,it)
c      ny=nny
c      write(0,*)' got grid y ',l1,' ',ny,'\n'
c      write(0,*)' from ',gy(1),' to ',gy(ny),'\n'
c      l1='gridz'
c      call getf11(10,gz,nnz,l1,l2,it)
c      nz=nnz
c      write(0,*)' got grid z ',l1,' ',nz,'\n'
c      write(0,*)' from ',gz(1),' to ',gz(nz),'\n'
c      close(10)
cc
c..... read fields
c      Bx,By,Bz in nT, rr (density) in cm**-3
c
c      call getarg(2,f3d)
      write (0,*) 'opening file: ',f2d
      open(10,file=f2d,status='old')
c
      l1='prec_e_fe_1'
      call getf21(10,prec_e_fe_1,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
      if ((nnx .lt. 0) .or. (nny .lt. 0)) then 
         close(10)
         return
      endif
c
c     setting up grids - there is no grid file
c
      write(0,*)' Setting up default grids in ionosphere: ', nnx,nny
      do 100 i=1,nnx
         gx(i)=180./120.*(i-1)
 100  continue
      do 200 i=1,nny
         gy(i)=i-1
 200  continue
      nx=nnx
      ny=nny
c
      l1='prec_e_e0_1'
      call getf21(10,prec_e_e0_1,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='prec_e_fe_2'
      call getf21(10,prec_e_fe_2,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='prec_e_e0_2'
      call getf21(10,prec_e_e0_2,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='pacurr'
      call getf21(10,pacurr,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='fac_tot'
      call getf21(10,fac_tot,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='pot'
      call getf21(10,pot,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='delphi'
      call getf21(10,delphi,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='ppio'
      call getf21(10,ppio,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='rrio'
      call getf21(10,rrio,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='ttio'
      call getf21(10,ttio,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='vdown'
      call getf21(10,vdown,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='sigh'
      call getf21(10,sigh,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='sigp'
      call getf21(10,sigp,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='fac_dyn'
      call getf21(10,fac_dyn,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='ctiot'
      call getf21(10,ctiot,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='ctiop'
      call getf21(10,ctiop,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='tau'
      call getf21(10,tau,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='ctaup'
      call getf21(10,ctaup,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='ctaut'
      call getf21(10,ctaut,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='cpolp'
      call getf21(10,cpolp,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='cpolt'
      call getf21(10,cpolt,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='delbr'
      call getf21(10,delbr,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='delbp'
      call getf21(10,delbp,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='delbt'
      call getf21(10,delbt,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='epio'
      call getf21(10,epio,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='etio'
      call getf21(10,etio,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      l1='xjh'
      call getf21(10,xjh,nnx,nny,l1,l2,it)
      write(0,*)' got ',l1,' ',nnx,nny
c
      close(10)
c
c..... ok, everything is in memory, do one interpolation for example
c      "iout" is returned 0 if (x,y,z) is in the box, 1 otherwise
      end


c-----------------------------------------------------------
      subroutine getf21(iu,a1,nx,ny,l1,l2,it)
c-----------------------------------------------------------
      real a1(*)
      real*8 rid
      character*80 rec,l1,l2
      m=len(l1)
      call end0(l1,m)
      isany=0
      if(l1(1:m).eq.'any') isany=1
100   continue
      read(iu,1000,end=190,err=190) rec
      write(0,*)'rec: ',m,rec
      if(rec(1:10).eq.'FIELD-2D-1') then
      read(iu,1000,end=190,err=190) rec
      write(0,*)'rec2: ',rec
      if((isany.eq.0).and.l1(1:m).ne.rec(1:m)) goto 100
      l1=rec
      read(iu,1000,end=190,err=190) l2
      write(0,*)'L2: ',l2
      read(iu,*,end=190,err=190)it,nx,ny
      call rdn2(iu,a1,nn,rec,it,rid)
      write(0,*)'NX,NY,NN: ',nx,ny,nn
      if(nn.lt.0)nx=nn
      return
      endif
      goto 100
190   continue
      nx=-1
1000  format(a)
      return
      end

