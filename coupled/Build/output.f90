      SUBROUTINE output (ng)
!
!svn $Id: output.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine manages nonlinear model output. It creates output   !
!  NetCDF files and writes out data into NetCDF files. If requested,   !
!  it can create several history and/or time-averaged files to avoid   !
!  generating too large files during a single model run.               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE distribute_mod, ONLY : mp_bcasts
      USE strings_mod,    ONLY : FoundError
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical :: Ldefine, NewFile
!
      integer :: Fcount, ifile, status, tile
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Nonlinear/output.F"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Turn on output data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 8, 67, MyFile)
!
!-----------------------------------------------------------------------
!  If appropriate, process nonlinear history NetCDF file.
!-----------------------------------------------------------------------
!
!  Set tile for local array manipulations in output routines.
!
      tile=MyRank
!
!  Turn off checking for analytical header files.
!
      IF (Lanafile) THEN
        Lanafile=.FALSE.
      END IF
!
!  Create output history NetCDF file or prepare existing file to
!  append new data to it.  Also,  notice that it is possible to
!  create several files during a single model run.
!
      IF (LdefHIS(ng)) THEN
        IF (ndefHIS(ng).gt.0) THEN
          IF (idefHIS(ng).lt.0) THEN
            idefHIS(ng)=((ntstart(ng)-1)/ndefHIS(ng))*ndefHIS(ng)
            IF (idefHIS(ng).lt.iic(ng)-1) THEN
              idefHIS(ng)=idefHIS(ng)+ndefHIS(ng)
            END IF
          END IF
          IF ((nrrec(ng).ne.0).and.(iic(ng).eq.ntstart(ng))) THEN
            IF ((iic(ng)-1).eq.idefHIS(ng)) THEN
              HIS(ng)%load=0                  ! restart, reset counter
              Ldefine=.FALSE.                 ! finished file, delay
            ELSE                              ! creation of next file
              Ldefine=.TRUE.
              NewFile=.FALSE.                 ! unfinished file, inquire
            END IF                            ! content for appending
            idefHIS(ng)=idefHIS(ng)+nHIS(ng)  ! restart offset
          ELSE IF ((iic(ng)-1).eq.idefHIS(ng)) THEN
            idefHIS(ng)=idefHIS(ng)+ndefHIS(ng)
            IF (nHIS(ng).ne.ndefHIS(ng).and.iic(ng).eq.ntstart(ng)) THEN
              idefHIS(ng)=idefHIS(ng)+nHIS(ng)  ! multiple record offset
            END IF
            Ldefine=.TRUE.
            NewFile=.TRUE.
          ELSE
            Ldefine=.FALSE.
          END IF
          IF (Ldefine) THEN                     ! create new file or
            IF (iic(ng).eq.ntstart(ng)) THEN    ! inquire existing file
              HIS(ng)%load=0                    ! reset filename counter
            END IF
            ifile=(iic(ng)-1)/ndefHIS(ng)+1     ! next filename suffix
            HIS(ng)%load=HIS(ng)%load+1
            IF (HIS(ng)%load.gt.HIS(ng)%Nfiles) THEN
              IF (Master) THEN
                WRITE (stdout,10) 'HIS(ng)%load = ', HIS(ng)%load,      &
     &                             HIS(ng)%Nfiles, TRIM(HIS(ng)%base),  &
     &                             ifile
              END IF
              exit_flag=4
              IF (FoundError(exit_flag, NoError,                        &
     &                       134, MyFile)) RETURN
            END IF
            Fcount=HIS(ng)%load
            HIS(ng)%Nrec(Fcount)=0
            IF (Master) THEN
              WRITE (HIS(ng)%name,20) TRIM(HIS(ng)%base), ifile
            END IF
            CALL mp_bcasts (ng, iNLM, HIS(ng)%name)
            HIS(ng)%files(Fcount)=TRIM(HIS(ng)%name)
            IF (HIS(ng)%ncid.ne.-1) THEN
              CALL netcdf_close (ng, iNLM, HIS(ng)%ncid)
            END IF
            CALL def_his (ng, NewFile)
            IF (FoundError(exit_flag, NoError, 149, MyFile)) RETURN
          END IF
          IF ((iic(ng).eq.ntstart(ng)).and.(nrrec(ng).ne.0)) THEN
            LwrtHIS(ng)=.FALSE.                 ! avoid writing initial
          ELSE                                  ! fields during restart
            LwrtHIS(ng)=.TRUE.
          END IF
        ELSE
          IF (iic(ng).eq.ntstart(ng)) THEN
            CALL def_his (ng, ldefout(ng))
            IF (FoundError(exit_flag, NoError, 159, MyFile)) RETURN
            LwrtHIS(ng)=.TRUE.
            LdefHIS(ng)=.FALSE.
          END IF
        END IF
      END IF
!
!  Write out data into history NetCDF file.  Avoid writing initial
!  conditions in perturbation mode computations.
!
      IF (LwrtHIS(ng)) THEN
        IF (LwrtPER(ng)) THEN
          IF ((iic(ng).gt.ntstart(ng)).and.                             &
     &        (MOD(iic(ng)-1,nHIS(ng)).eq.0)) THEN
            IF (nrrec(ng).eq.0.or.iic(ng).ne.ntstart(ng)) THEN
              CALL wrt_his (ng, tile)
            END IF
            IF (FoundError(exit_flag, NoError, 176, MyFile)) RETURN
          END IF
        ELSE
          IF (MOD(iic(ng)-1,nHIS(ng)).eq.0) THEN
            CALL wrt_his (ng, tile)
            IF (FoundError(exit_flag, NoError, 181, MyFile)) RETURN
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  If appropriate, process nonlinear quicksave NetCDF file.
!-----------------------------------------------------------------------
!
!  Create output quicksave NetCDF file or prepare existing file to
!  append new data to it.  Also,  notice that it is possible to
!  create several files during a single model run.
!
      IF (LdefQCK(ng)) THEN
        IF (ndefQCK(ng).gt.0) THEN
          IF (idefQCK(ng).lt.0) THEN
            idefQCK(ng)=((ntstart(ng)-1)/ndefQCK(ng))*ndefQCK(ng)
            IF (idefQCK(ng).lt.iic(ng)-1) THEN
              idefQCK(ng)=idefQCK(ng)+ndefQCK(ng)
            END IF
          END IF
          IF ((nrrec(ng).ne.0).and.(iic(ng).eq.ntstart(ng))) THEN
            IF ((iic(ng)-1).eq.idefQCK(ng)) THEN
              QCK(ng)%load=0                  ! restart, reset counter
              Ldefine=.FALSE.                 ! finished file, delay
            ELSE                              ! creation of next file
              Ldefine=.TRUE.
              NewFile=.FALSE.                 ! unfinished file, inquire
            END IF                            ! content for appending
            idefQCK(ng)=idefQCK(ng)+nQCK(ng)  ! restart offset
          ELSE IF ((iic(ng)-1).eq.idefQCK(ng)) THEN
            idefQCK(ng)=idefQCK(ng)+ndefQCK(ng)
            IF (nQCK(ng).ne.ndefQCK(ng).and.iic(ng).eq.ntstart(ng)) THEN
              idefQCK(ng)=idefQCK(ng)+nQCK(ng)  ! multiple record offset
            END IF
            Ldefine=.TRUE.
            NewFile=.TRUE.
          ELSE
            Ldefine=.FALSE.
          END IF
          IF (Ldefine) THEN                     ! create new file or
            IF (iic(ng).eq.ntstart(ng)) THEN    ! inquire existing file
              QCK(ng)%load=0                    ! reset filename counter
            END IF
            ifile=(iic(ng)-1)/ndefQCK(ng)+1     ! next filename suffix
            QCK(ng)%load=QCK(ng)%load+1
            IF (QCK(ng)%load.gt.QCK(ng)%Nfiles) THEN
              IF (Master) THEN
                WRITE (stdout,10) 'QCK(ng)%load = ', QCK(ng)%load,      &
     &                             QCK(ng)%Nfiles, TRIM(QCK(ng)%base),  &
     &                             ifile
              END IF
              exit_flag=4
              IF (FoundError(exit_flag, NoError,                        &
     &                       236, MyFile)) RETURN
            END IF
            Fcount=QCK(ng)%load
            QCK(ng)%Nrec(Fcount)=0
            IF (Master) THEN
              WRITE (QCK(ng)%name,20) TRIM(QCK(ng)%base), ifile
            END IF
            CALL mp_bcasts (ng, iNLM, QCK(ng)%name)
            QCK(ng)%files(Fcount)=TRIM(QCK(ng)%name)
            IF (QCK(ng)%ncid.ne.-1) THEN
              CALL netcdf_close (ng, iNLM, QCK(ng)%ncid)
            END IF
            CALL def_quick (ng, NewFile)
            IF (FoundError(exit_flag, NoError, 251, MyFile)) RETURN
          END IF
          IF ((iic(ng).eq.ntstart(ng)).and.(nrrec(ng).ne.0)) THEN
            LwrtQCK(ng)=.FALSE.                 ! avoid writing initial
          ELSE                                  ! fields during restart
            LwrtQCK(ng)=.TRUE.
          END IF
        ELSE
          IF (iic(ng).eq.ntstart(ng)) THEN
            CALL def_quick (ng, ldefout(ng))
            IF (FoundError(exit_flag, NoError, 261, MyFile)) RETURN
            LwrtQCK(ng)=.TRUE.
            LdefQCK(ng)=.FALSE.
          END IF
        END IF
      END IF
!
!  Write out data into quicksave NetCDF file.  Avoid writing initial
!  conditions in perturbation mode computations.
!
      IF (LwrtQCK(ng)) THEN
        IF (LwrtPER(ng)) THEN
          IF ((iic(ng).gt.ntstart(ng)).and.                             &
     &        (MOD(iic(ng)-1,nQCK(ng)).eq.0)) THEN
            IF (nrrec(ng).eq.0.or.iic(ng).ne.ntstart(ng)) THEN
              CALL wrt_quick (ng, tile)
            END IF
            IF (FoundError(exit_flag, NoError, 278, MyFile)) RETURN
          END IF
        ELSE
          IF (MOD(iic(ng)-1,nQCK(ng)).eq.0) THEN
            CALL wrt_quick (ng, tile)
            IF (FoundError(exit_flag, NoError, 283, MyFile)) RETURN
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  If appropriate, process restart NetCDF file.
!-----------------------------------------------------------------------
!
!  Create output restart NetCDF file or prepare existing file to
!  append new data to it.
!
      IF (LdefRST(ng)) THEN
        CALL def_rst (ng)
        IF (FoundError(exit_flag, NoError, 555, MyFile)) RETURN
        LwrtRST(ng)=.TRUE.
        LdefRST(ng)=.FALSE.
      END IF
!
!  Write out data into restart NetCDF file.
!
      IF (LwrtRST(ng)) THEN
        IF ((iic(ng).gt.ntstart(ng)).and.                               &
     &      (MOD(iic(ng)-1,nRST(ng)).eq.0)) THEN
          CALL wrt_rst (ng)
          IF (FoundError(exit_flag, NoError, 569, MyFile)) RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Turn off output data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iNLM, 8, 615, MyFile)
!
 10   FORMAT (/,' OUTPUT - multi-file counter ',a,i0,                   &
     &          ', is greater than Nfiles = ',i0,1x,'dimension',        &
     &        /,10x,'in structure when creating next file: ',           &
     &           a,'_',i4.4,'.nc',                                      &
     &        /,10x,'Incorrect OutFiles logic in ''read_phypar''.')
 20   FORMAT (a,'_',i5.5,'.nc')
!
      RETURN
      END SUBROUTINE output
