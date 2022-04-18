! **********************************************************************
!
! Purpose :: Provide developers with program status detail
!
! Programmer :: Douglas Larson 
! Creation Date :: Saturday, September 21, 2002
!
! Integer Variables:
!     ------------------
      ! 0-No messages, 1-Simple messages, 2-verbose
      integer debug_level
      ! Unit number for the debug output log
      integer dbglog      
      common/dbgmode/debug_level,dbglog
      SAVE /dbgmode/
! ======================================================================
