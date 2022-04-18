! Subroutine to open the FileName grid data file and fill the 
!     gx,gy,gz (cell centerd grid) and 
!      (gx_bx1,gy_bx1, ... ,gz_bz1) defining the staggered grids for B fields
!
! Modification history:
!  2020/10/20 Lutz Rastaetter - added to read encoded file with staggered grids
!  2022/03/23 Lutz Rastaetter - adjust to stricter compilation with f2py
!                             
!      
      subroutine read_grid_dir(FileName,gridName,
     *       gridarr,nn)

      implicit none
      
      INTEGER    fileNoGrid
      INTEGER NN
      character*120 FileName,gridName
      character*80 l2
      real gridarr(*)
      integer it,my_iostat
!f2py intent(in) :: FileName
!f2py intent(in) :: gridName
!f2py intent(out) :: nn
!f2py intent(in,out) :: gridarr
      fileNoGrid=10
      write(*,*)filename
      write(*,*)gridname
c      open(10,file=tmp_str1(1:fgrid_slen),status='old')
      open(fileNoGrid,file=filename,status='old',iostat=my_iostat)
c
c... read grid parameters
c    gridpoint (ix,iy,iz) is at (gx(ix),gy(iy),gz(iz)) in RE
c
      if (my_iostat .eq. 0) then
        call getf11(fileNoGrid,gridarr,nn,gridName,l2,it)
        close(fileNoGrid)
      else
        write(*,*) "file ",filename," not found: "
        nn=-1
      endif
      
      end      
      
      subroutine read_grid_for_vector(FileName,
     *  vectorCompName,nx,ny,nz,gx_v,gy_v,gz_v)

      implicit none
C	Declare a parameter for the length of the Fortran character strings

      INTEGER CHAR_SIZE,nx,ny,nz
      INTEGER filenoGrid
      INTEGER nx_g
      PARAMETER (CHAR_SIZE = 80)
c this is the only one used to allocate a temporary array xy_tmp
      PARAMETER (NX_G = 5000 )

      character*7 vectorCompName
      character*80 gx_vname,gy_vname,gz_vname
      character*120 ,fileName

      real gx_v(*),gy_v(*),gz_v(*)
      integer it,my_iostat
      character*80 l2
!f2py intent(in) :: FileName
!f2py intent(in) :: vectorCompName
!f2py intent(out) :: nx,ny,nz
!f2py intent(in,out) :: gx_v,gy_v,gz_v
      fileNoGrid=10
      my_iostat=0
c      open(10,file=tmp_str1(1:fgrid_slen),status='old')
      open(fileNoGrid,file=filename,status='old',iostat=my_iostat)
      write(*,*) 'vector component name: ',vectorCompName
      write(*,*) 'String length: ',len(vectorCompName)
      if (vectorCompName(1:1) .ne. ' ') then
         gx_vname='gx-'//vectorCompName
         gy_vname='gy-'//vectorCompName
         gz_vname='gz-'//vectorCompName
      else
         gx_vname='gridx'
         gy_vname='gridy'
         gz_vname='gridz'
      endif
      write(*,*)filename
      write(*,*) gx_vname,gy_vname,gz_vname

      if (my_iostat .eq. 0) then
c     
c... read grid parameters
c    gridpoint (ix,iy,iz) is at (gx(ix),gy(iy),gz(iz)) in RE
c
         call getf11(fileNoGrid,gx_v,nx,gx_vname,l2,it)
         call getf11(fileNoGrid,gy_v,ny,gy_vname,l2,it)
         call getf11(fileNoGrid,gz_v,nz,gz_vname,l2,it)
         close(fileNoGrid)
      else
         nx=-1
         ny=-1
         nz=-1
      endif
      end

      
c-----------------------------------------------------------
      subroutine getf11(iu,a1,nx,l1,l2,it)
c-----------------------------------------------------------
      implicit none
      real a1(*)
      character*80 rec
      character*80 l1,l2
      integer iu,it,nx,isany,m
      real*8 rid
      m=len(l1)
      call end0(l1,m)
      isany=0
      if(l1(1:m).eq.'any') isany=1
100   continue
      read(iu,1000,end=190) rec
      if(rec(1:10).eq.'FIELD-1D-1') then
      read(iu,1000,end=190) rec
      if((isany.eq.0).and.l1(1:m).ne.rec(1:m)) goto 100
      l1=rec
      read(iu,1000,end=190) l2
      read(iu,*,end=190)it,nx
      call rdn2(iu,a1,nx,rec,it,rid)
      return
      endif
      goto 100
190   continue
      nx=-1
1000  format(a)
      return
      end

c-----------------------------------------------------------
      subroutine end0(r,m)
c-----------------------------------------------------------
      character*80 r
      integer m,n
      n=min0(m,len(r))
      do 100 i=1,n
      m=i-1
      if(r(i:i).eq.' ') return
100   continue
      m=n
      return
      end
