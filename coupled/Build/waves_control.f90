      MODULE waves_control_mod
!
!svn $Id: waves_control.F 814 2008-10-29 01:42:17Z jcwarner $
!================================================== John C. Warner  ====
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  SWAN model:                                                         !
!                                                                      !
!  This driver executes SWAN by controlling initialization,            !
!  time-stepping, and finalization for nested grids.                   !
!                                                                      !
!     SWAN_initialize                                                  !
!     SWAN_run                                                         !
!     SWAN_finalize                                                    !
!                                                                      !
!=======================================================================
!
      USE mct_coupler_params
      USE M_COUPLING
      implicit none
      PRIVATE
      PUBLIC  :: SWAN_driver_init
      PUBLIC  :: SWAN_driver_run
      PUBLIC  :: SWAN_driver_finalize
      CONTAINS
      SUBROUTINE SWAN_driver_init (MyCOMM)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes SWAN variables               !
!  and internal and external parameters.                               !
!                                                                      !
!=======================================================================
!
      USE swan_iounits
      USE mod_coupler_kinds
      USE M_PARALL
      USE M_PARALL2
      USE SWCOMM3
      USE TIMECOMM
      USE M_MPI
      USE WAVES_COUPLER_MOD, ONLY: INITIALIZE_WAV_ROUTERS
!
!  Imported variable declarations.
!
      integer, intent(in) :: MyCOMM
!
!  Local variable declarations.
!
      integer :: ng, i, iw, MyRank, MyError
      integer :: ngp, ngc
      character (len=160) :: tempname
!
!  Initialize the time counter.
!
      IF (.not.ALLOCATED(iics)) ALLOCATE (iics(NUM_SGRIDS))
      DO iw=1,NUM_SGRIDS
        iics(iw)=0
      END DO
!
!  Initialize restart counters.
!
      DO iw=1,NUM_SGRIDS
        dtswanrst(iw)=0.
        SwanRstFnum(iw)=0
      END DO
!
!  Initialize the grids.
!
      DO ng=1,NUM_SGRIDS
        ngp=1
        ngc=1
        CALL SWAN_INITIALIZE (ng, ngc, ngp, NUM_SGRIDS, MyCOMM, Wname(ng))
        CALL SWSYNC
      END DO
!
!  Initialize the MCT routers to ROMS. This has to be
!  done outside an 'ng' loop so that it is synchronous 
!  with ROMS.
!
      CALL INITIALIZE_WAV_ROUTERS
!
!  The call to run here does not do a time step, it fills the bc arrays,
!  fill AC2 array of bound spec data for child grids, and enters into MCT.
!
      DO ng=1,NUM_SGRIDS
        CALL SWSYNC
        ngp=1
        ngc=1
        CALL SWAN_RUN (0, ng, ngc, ngp, 0, NUM_SGRIDS)
        CALL SWAN_RST (0, ng)
      END DO
!
!     --- couple models during output computations
!
        CALL SWSYNC
        CALL SWAN_CPL (0)
      RETURN
      END SUBROUTINE SWAN_driver_init
      SUBROUTINE SWAN_driver_run
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes ROMS/TOMS state variables    !
!  and internal and external parameters.                               !
!                                                                      !
!=======================================================================
!
      USE swan_iounits
      USE mct_coupler_params
      USE M_PARALL
      USE M_PARALL2
      USE SWCOMM3
      USE TIMECOMM
      USE M_MPI
!
!  Imported variable declarations.
!
!  Local variable declarations.
!
      integer :: i, iw, ng, MyRank, MyError
      integer, allocatable :: run_grid(:)
      integer, allocatable :: run_cplr(:)
      real :: rtime, rtime_start, cff
      integer :: ngp, ngc
      CALL mpi_comm_rank (WAV_COMM_WORLD, MyRank, MyError)
!
!  Set some initial run time parameters here.
!
      IF (.not.ALLOCATED(run_grid)) ALLOCATE (run_grid(NUM_SGRIDS))
      IF (.not.ALLOCATED(run_cplr)) ALLOCATE (run_cplr(NUM_SGRIDS))
      DO iw=1,NUM_SGRIDS
        run_grid(iw)=1
        run_cplr(iw)=1
      END DO
      rtime_start=0.
      rtime=rtime_start
!
!  Main job control loop here.
!
      DO WHILE (iics(NUM_SGRIDS).lt.MTC_G(NUM_SGRIDS))
!
!  Advance grids in time that have run_grid flag == 1
!  For the first entry, all grids step individual dts.
!
        DO ng=1,NUM_SGRIDS
          IF (run_grid(ng).eq.1) THEN
            iics(ng)=iics(ng)+1
            ngp=1
            ngc=1
            CALL SWAN_RUN (iics(ng), ng, ngc, ngp, 1, NUM_SGRIDS)
            run_grid(ng)=0
!           CALL SWAN_RST (iics(ng), ng)
          END IF
        END DO
!
!  Advance the time counter by the smallest dt.
!
        rtime=rtime+DT_G(NUM_SGRIDS)
!
!  Determine what grids can be time stepped. This is determined
!  by comparing dt(each grid) to global time rtime.
!
        DO ng=1,NUM_SGRIDS
          cff=rtime-rtime_start
          IF (MOD(cff,REAL(DT_G(ng))).eq.0) THEN
            run_grid(ng)=1
          END IF
        END DO
        IF (run_grid(1).eq.1) THEN
!
!     --- receive data from other coupled models.
          CALL SWSYNC
          CALL SWAN_CPL (1)
        END IF
        DO ng=1,NUM_SGRIDS
          IF (run_grid(1).eq.1) THEN
            CALL SWAN_RST (iics(ng), ng)
          END IF
        END DO
      END DO
      IF (ALLOCATED(run_grid)) DEALLOCATE (run_grid)
      IF (ALLOCATED(run_cplr)) DEALLOCATE (run_cplr)
      IF (ALLOCATED(iics)) DEALLOCATE (iics)
      RETURN
      END SUBROUTINE SWAN_driver_run
      SUBROUTINE SWAN_driver_finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates SWAN.                                       !
!                                                                      !
!=======================================================================
!
      USE swan_iounits
!
!  Local variable declarations.
!
      integer :: ng
      DO ng=1,NUM_SGRIDS
        CALL SWAN_FINALIZE (ng, NUM_SGRIDS)
      END DO
      RETURN
      END SUBROUTINE SWAN_driver_finalize
      END MODULE waves_control_mod
