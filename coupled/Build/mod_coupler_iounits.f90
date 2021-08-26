      MODULE mod_coupler_iounits
!
!svn $Id: mod_coupler_iounits.F jcwarner $
!=======================================================================
!                                                                      !
!  stdinp      Unit number for standard input (often 5).               !
!  stdout      Unit number for standard output (often 6).              !
!  Aname       Atmosphere model stadard input file name.               !
!  IWBNDname   Input boundary data file name for InWave model          !
!  IWSWNname   Input spectral SWAN data file name for InWave model     !
!                                                                      !
!=======================================================================
!
      USE mct_coupler_params
      USE swan_iounits
      USE mod_param
      implicit none
!      integer, parameter :: IOnamesize = 160
!      character (len=IOnamesize) :: Wname
! this flag is temporary for SCRIP option.
      integer :: scrip_opt
      character (len=80) :: SCRIPname
      CONTAINS
      SUBROUTINE allocate_coupler_iounits
!=======================================================================
!                                                                      !
!  This routine initialize all the coupler io vars.                    !
!                                                                      !
!=======================================================================
      character (len=1 ), parameter :: blank = ' '
      integer :: i, io, ia, iw, ih
      RETURN
      END SUBROUTINE allocate_coupler_iounits
      END MODULE mod_coupler_iounits
