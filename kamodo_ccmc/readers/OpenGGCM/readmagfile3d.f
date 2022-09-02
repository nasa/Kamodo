! Read 3D field threeDfield named var from path3d (*.edf* output file).
! Also obtain grid dimension xn.ny.nz and time string (asciiTime).
!
! Modification History:
!  2022/03/23 Lutz Rastaetter call_read() renamed to read_3d_field()
!                             dropped logical unit number in input of routine
!
c----------------------------------------------------------------------      
      subroutine read_3d_field(path3d,
     &              threeDfield,var,nx,ny,nz,asciiTime)
c----------------------------------------------------------------------              
!      character(len=40) :: path3d
      character(len=110) :: path3d           
      character*80 var
      character*80 asciiTime
      real threeDfield(nx,ny,nz) 
      integer fileNo2
      integer nx
      integer ny
      integer nz
      integer it
!f2py intent(in) :: path3d
!f2py intent(in,out) :: threeDfield
!f2py intent(in) :: var
!f2py intent(out) :: nx,ny,nz
!f2py intent(out) :: asciiTime
      fileNo2=11
c     write (0,*) 'opening file: |',path3d,'|'              
      open(unit=fileNo2,file=path3d,status='old',action='read',err=222)
      call read_mhd3d(fileNo2,threeDfield,nx,ny,nz,var,asciiTime,it)
      close(fileNo2)
      return
 222  write(0,*) "could not open file"
      return
      end
      
! Change log ----
!
! Burlen Loring 2008-01-03
! 
! added comments. I am looking for a place to OpenMPI-ify
! to leverage our multicore systems. But, not here or any lower level 
! subs without a major rewrite. :( See notes in lower level subs
!
!
c     ChangeLog--1/21/10: Mgilson added a rewind in the end of file
c     condition.
c     capabilities
c     ChangeLog--2/1/10: Mgilson refactored getf31.  Now it detects
c     filetype and reads accordingly.
c-----------------------------------------------------------
      subroutine read_mhd3d( fileNo,
     &                   fieldData,
     &                   nx,ny,nz,
     &                   fieldName,
     &                   asciiTime,
     &                   it )
! nx,ny,nz are dimensions of the array
! nx is set to < 0 if an error occured
! other wise the decompressed field specified by fieldName
! is returned in fieldData along with asciiTime and 'it'.
! still working on what 'it' is... may be integerTime or itteration 
! Note: do not implicit none this unless you are certain
! that you can count on fieldName and asciiTime being a 
! a specific fixed length
c-----------------------------------------------------------
      include 'DBGDEV.F'
! dbg_level set to zero or negative to get diagnostic outputs
      integer dbg_level
      parameter (dbg_level=1) 
      real fieldData(nx,ny,nz)
      character*80 thisFieldType
      character*80 thisFieldName
      character*80 fieldName
      character*80 asciiTime
      integer isany
      integer fileNo
      integer nx
      integer ny
      integer nz
      integer nn
      integer it

!      print *, 'debug_level= ', debug_level.gt.dbg_level-1
      if(debug_level.gt.dbg_level-1)then                     !(.gt. >)
         write(dbglog,'("ENTERING: read_mhd3d")')
         write(dbglog,'(3X,"nx,ny,nz=",3(1X,I8))')nx,ny,nz
      endif
!      print *, 'hello there'

c     if ascii, execute this branch      
c==============================================      
      ! string length, needed for string comparisions
      m=len(fieldName)
cc      write(*,'("m=",I3)')m
c      call end0(fieldName,m)

      ! determine if the string is "any"
      isany=0
      if(fieldName(1:m).eq.'any') isany=1
      ! read records of the file until we find what we are looking for.
      ! or we hit the end of the file.
100   continue

!      print *, 'fileNo1= ', fileNo
      ! read field type
      read( UNIT=fileNo,FMT='(A)',END=190, ERR=190 ) thisFieldType
!    print *, 'thisFieldType= ', thisFieldType

        if(debug_level.gt.dbg_level+1)then
           write(dbglog,'("thisFieldType=",A10)')thisFieldType(1:10)
        endif

      ! we are looking for a 3D field and this could be it
      if ( thisFieldType(1:10).eq.'FIELD-3D-1' ) then

        ! get this fields name
        read( UNIT=fileNo,
     &        FMT='(A)',
     &        END=190 ) thisFieldName


        if(debug_level.gt.dbg_level)then
           write(dbglog,'("thisFieldName=",A)')thisFieldName
        endif

!       write(*,*) '<'//fieldName(1:m)//'>','<'//thisFieldName(1:m)//'>'
!       write(*,*) fieldName(1:m).ne.thisFieldName(1:m)

        ! this isn't the one we want
cc        write(*,'("m=",I3,2(1X,A))')m,fieldname,thisFieldName
        if (  (isany.eq.0).and.
     &        (fieldName(1:m).ne.thisFieldName(1:m)) ) then
           goto 100
        endif

        ! this is the one we want (but why re-asign??)
        fieldName=thisFieldName

        ! read the time string
        read( UNIT=fileNo,
     &        FMT='(A)',
     &        END=190 ) asciiTime

        ! read the array dimensions and 'it'
        read( UNIT=fileNo,
     &        FMT=*,
     &        END=190 ) it,nx,ny,nz

        ! read and decompress the data
        call rdn2( fileNo,fieldData,nn,thisFieldName,it,rid )
        if(debug_level.gt.dbg_level-1)then
          print *, 'fieldData= ', shape(fieldData)
          print *, 'data sample: ', fieldData(1,1,1)
        endif
!     output_fieldData=fieldData       !--reasign for f2py function jcc 09/20/19
        ! an error occured, set nx to error code 
        if(nn.lt.0)nx=nn

        ! done
        return
      endif

      ! we didn't find what we wanted so we try again
      goto 100

      ! we hit EOF before we found what we wanted
190   continue
      nx=-1
      rewind(fileNo)                ! hit end of file, rewind so that subsequent searches make sense
      if(debug_level.gt.dbg_level-1)then
         write(dbglog,'(3X,"x: nx,ny,nz=",3(1X,I8))')nx,ny,nz
         write(dbglog,*) "EXITING: read_mhd3d"
      endif
      return
      end

! Change log---
! 
! Burlen Loring 2008-01-03
! 
! added comments, indetantion, and implicit none
! I had wanted to OpenMPI the decompression loops
! however a major rewrite would be needed since these are
! interleaved with *sequential* reads of the data file :(
c.---------------------------------------
      subroutine rdn2( fileNo,a,n,cid,it,rid )
!
!
!
c character decoding
c.---------------------------------------
      real a(*)
      character*8 cid
      character*4 did
      real a1(0:63)
      integer intRec1(0:63),intRec2(0:63),i3(0:63)
      integer q
      integer fileNo
      integer n,it
      real zmin,zmax,rid
100   continue
      read( UNIT=fileNo,
     &      FMT='(a,i8,3e14.7,i8,a)',
     &      ERR=100,
     &      END=900 )did,n,zmin,zmax,rid,it,cid
      ! The format of the header record is as follows:
      ! 1, 4 char string                                    ( did )
      ! 8 chars to 1 int                                    ( n )
      ! 3 times 14 chars to double w/ 7 digits precision    ( zmin,zmax,rid )
      ! 8 chars to 1 int                                    ( it )
      ! 1, 8 char string                                    ( cid )
!      print *, 'second routine'
      ! try read again, if we don't have WRN2 data
      if( did.ne.'WRN2' ) then
        goto 100
      endif

      ! data array is constant, initialize a and return
      if(zmin.eq.zmax) then
        do i=1,n
          a(i)=zmin
        enddo
        return
      ! compute data resolution
      else
        ! what is the significance of 4410??
        ! and why the explicit cast rather than 4410.0??
        dzi=(zmax-zmin)/float(4410)
      endif

      ! process data records, these are 64 characters long
      q=0 ! index into decompressed real*8 data
      do k=1,n,64

        ! records should be 64 ints long but last one may be shorter
        nk = min0(63,n-k)

        ! decompress a record( expands and converts ascii to hex)
        call wrndec(fileNo,intRec1,nn)

        ! error in record length bail
        ! note: if returned nn=-5 explicitly indicates an error
        if(nn.ne.nk) then
          write(0,*)'rdn2: nn .ne. nk ',nn,nk,n,k
          n=-2
          return
        endif

        ! decompress next record( expands and converts ascii to hex)
        call wrndec(fileNo,intRec2,nn)

        ! error in record length bail
        ! note: if returned nn=-5 explicitly indicates an error
        if(nn.ne.nk) then
          write(0,*)'rdn2: nn .ne. nk ',nn,nk
          n=-2
          return
        endif

        ! convert these records from hex-ascii to double
        ! not sure about the specifics here except that
        ! this is the propritary OpenGGCM format
        ! what is significance of 33,47 and 94 ??
        do i=0,nk
          intRec1(i) = intRec1(i) - 33
          intRec2(i) = intRec2(i) - 33
          sig        = 1.

          if( intRec1(i).ge.47 ) then
            sig        = -1.
            intRec1(i) = intRec1(i) - 47
          endif

          i3(i) = intRec2(i) + 94*intRec1(i)
          a1(i) = FLOAT(i3(i))
          a1(i) = dzi*a1(i) + zmin
          a1(i) = sig*exp(a1(i))

          q    = q + 1
          a(q) = a1(i)
        enddo
      enddo

      return

      ! io error occured
 900  continue
      n=-1
      return
      end

! Change log ---
! Burlen Loring 2008-01-03
! commented, renamed some variables and adjusted syntax
c.---------------------------------------
      subroutine wrndec(fileNo,expandedDatOut,n)
!
! decode compressed data sequence read in from a file
! will decode a single record (typically 64 chars ion length)
!
! the decompressed record is returned in exdpandedDataOut
! the legth of the de-compressed record is returned in n
! n=-5 indicates something bad happened
!
c.---------------------------------------
      IMPLICIT NONE
      integer fileNo
      integer n
      integer nRep,hexChar
      integer checkSum
      integer i,j,q
      ! results returned
      integer expandedDatOut(0:63)
      ! work array has expanded data + check sum
      integer expandedDat(-2:67)
      ! read buffer
      character*72 asciiDat        

      ! initialize some of the variables
      ! note: read buffer is space padded, spaces indicate end of field
      n=-5
      expandedDatOut(0)=0
      asciiDat='                                    '//
     &'                                    '

      ! read 72 characters from file(record assumed to be 64 chars)
      read( UNIT=fileNo,
     &      FMT='(A)',
     &      END=900,
     &      ERR=900 ) asciiDat

!      print *, 'third routine'
      ! expand sequences in buffer. An ascii value with 7th bit set
      ! is a count(2^7=128), this tells how many times following ascii 
      ! char is repeated. If the 7th bit is not set then the ascii char
      ! appears once.
      ! this is an infinite loop which breaks once a space is encountered
      ! or too many chars have been processed (ie record length exceeded).
      i=-2 ! index into expanded data
      j=0  ! index into read buffer
100   continue
        j = j + 1

        ! error records should not be this long, bail 
        if(j.gt.67) then
          write(0,*)'wrndec: cannot find end of encoded record'
          n=-5
          return
        endif

        ! convert j-th value from ascii to hex
        hexChar=ICHAR( asciiDat(j:j) )

        ! if its a space then we are done expanding this record, break
        if ( hexChar.eq.32 ) then
          goto 190
        endif

        ! if its not an 7 bit ascii character
        ! then it tells how many time to repeat the following ascii char
        if ( hexChar.gt.127 ) then

          nRep = hexChar-170 ! convert to count ??

          ! convert next j-th value from ascii to hex
          ! this is repeated nRep times
          j  = j + 1
          hexChar = ICHAR( asciiDat(j:j) )
          do q=1,nRep
            i = i + 1
            expandedDat(i) = hexChar
          enddo

        ! otherwise it is 7 bit ascii, and apears only once
        else
          i = i + 1
          expandedDat(i) = hexChar
        endif

      ! expand next sequence in buffer
      goto 100

190   continue

      ! error we expect field to be 64 chars if its more then this
      ! there is a problem, bail
      n=i
      if( n.gt.63 ) then
        write(0,*)'wrndec: n gt 63, n=',n
        n=-5
        write(0,'(a)')'rec:'
        write(0,'(a)')asciiDat
        return
      endif

      ! copy and compute check sum 
      checkSum=0
      do i=0,n
         ! copy
         expandedDatOut(i) = expandedDat(i)
         ! compute check sum
         checkSum = checkSum + expandedDatOut(i) 
      enddo

      ! check sum error bail
      if(expandedDat(-1).ne.33+mod(checkSum,92)) then
        write(0,*)'wrndec: checksum error '
        write(0,'(a,a)')'rec:',asciiDat
        write(0,*)expandedDat(-1),33+mod(checkSum,92)
        n=-5
        return
      endif

      return

      ! file io error 
900   continue
        write(0,*)' wrndec eof/err '
        n=-5
        return

      end 
