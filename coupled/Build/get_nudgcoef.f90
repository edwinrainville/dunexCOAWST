      SUBROUTINE get_nudgcoef (ng, model)
!
!svn $Id: get_nudgcoef.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine reads grid time scales for nudging to climatology   !
!  fields.                                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_clima
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE exchange_2d_mod
      USE exchange_3d_mod
      USE mp_exchange_mod, ONLY : mp_exchange2d
      USE mp_exchange_mod, ONLY : mp_exchange3d
      USE mp_exchange_mod, ONLY : mp_exchange4d
      USE nf_fread2d_mod,  ONLY : nf_fread2d
      USE nf_fread3d_mod,  ONLY : nf_fread3d
      USE strings_mod,     ONLY : FoundError, find_string
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      logical :: got_generic
      logical :: got_specific(NT(ng))
!
      integer :: tile, LBi, UBi, LBj, UBj
      integer :: gtype, i, ic, itrc
      integer :: nvatt, nvdim, status, vindex
      integer :: Vsize(4)
!
      real(dp) :: Fscl
      real(r8) :: Fmax, Fmin
!
      character (len=40 ) :: tunits
      character (len=256) :: ncname
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/get_nudgcoef.F"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Inquire about the contents of grid NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, 80, MyFile)) RETURN
      ncname=NUD(ng)%name
!
!  Open grid NetCDF file for reading.
!
      IF (NUD(ng)%ncid.eq.-1) THEN
        CALL netcdf_open (ng, model, ncname, 0, NUD(ng)%ncid)
        IF (FoundError(exit_flag, NoError, 87, MyFile)) THEN
          WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
      END IF
!
!  Check grid file dimensions for consitency
!
      CALL netcdf_check_dim (ng, model, ncname,                         &
     &                       ncid = NUD(ng)%ncid)
      IF (FoundError(exit_flag, NoError, 97, MyFile)) RETURN
!
!  Inquire about the variables.
!
      CALL netcdf_inq_var (ng, model, ncname,                           &
     &                     ncid = NUD(ng)%ncid)
      IF (FoundError(exit_flag, NoError, 103, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Check if required variables are available.
!-----------------------------------------------------------------------
!
!  Nudging coefficients for 2D momentum.
!
      IF (LnudgeM2CLM(ng)) THEN
        IF (.not.find_string(var_name,n_var,Vname(1,idM2nc),            &
     &                       NUD(ng)%Vid(idM2nc))) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idM2nc)),          &
     &                                  TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END IF
!
!  Nudging coefficients for 3D momentum.
!
      IF (LnudgeM3CLM(ng)) THEN
        IF (.not.find_string(var_name,n_var,Vname(1,idM3nc),            &
     &                       NUD(ng)%Vid(idM3nc))) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idM3nc)),          &
     &                                  TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END IF
!
!  Tracers nudging coefficients.
!
      IF (ANY(LnudgeTCLM(:,ng))) THEN
        got_generic=find_string(var_name,n_var,Vname(1,idgTnc),         &
     &                          NUD(ng)%Vid(idgTnc))
        DO itrc=1,NT(ng)
          IF (LnudgeTCLM(itrc,ng)) THEN
            vindex=idTnud(itrc)
            got_specific(itrc)=find_string(var_name,n_var,              &
     &                                     Vname(1,vindex),             &
     &                                     NUD(ng)%Vid(vindex))
            IF (.not.(got_generic.or.got_specific(itrc))) THEN
              IF (Master) WRITE (stdout,30) TRIM(Vname(1,idgTnc)),      &
     &                                      TRIM(Vname(1,vindex)),      &
     &                                      TRIM(ncname)
              exit_flag=2
              RETURN
            END IF
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Read in nudging coefficients.
!-----------------------------------------------------------------------
!
!  Set 2D arrays bounds.
!
      tile=MyRank
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!  Set Vsize to zero to deativate interpolation of input data to model
!  grid in "nf_fread2d".
!
      DO i=1,4
        Vsize(i)=0
      END DO
!
!  If appropriate, read nudging coefficients for 2D momentum  (inverse
!  time scales, 1/day).
!
      IF (Master) WRITE (stdout,'(1x)')
!
      IF (LnudgeM2CLM(ng)) THEN
        gtype=r2dvar
        vindex=idM2nc
        Fscl=1.0_dp/day2sec                    ! default units: 1/day
!
        CALL netcdf_inq_var (ng, model, ncname,                         &
     &                       ncid = NUD(ng)%ncid,                       &
     &                       MyVarName = TRIM(Vname(1,vindex)),         &
     &                       VarID = NUD(ng)%Vid(vindex),               &
     &                       nVarDim = nvdim,                           &
     &                       nVarAtt = nvatt)
        IF (FoundError(exit_flag, NoError, 197, MyFile)) RETURN
        DO i=1,nvatt
          IF (TRIM(var_Aname(i)).eq.'units') THEN
            tunits=TRIM(var_Achar(i))
            IF (tunits(1:3).eq.'day') THEN
              Fscl=1.0_dp/day2sec
            ELSE IF (tunits(1:6).eq.'second') THEN
              Fscl=1.0_dp
            END IF
          END IF
        END DO
!
        status=nf_fread2d(ng, model, ncname, NUD(ng)%ncid,              &
     &                    Vname(1,vindex), NUD(ng)%Vid(vindex),         &
     &                    0, gtype, Vsize,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Fscl, Fmin, Fmax,                             &
     &                    GRID(ng) % rmask,                             &
     &                    CLIMA(ng) % M2nudgcof)
        IF (FoundError(status, nf90_noerr, 223, MyFile)) THEN
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,vindex))//': '//             &
     &                        TRIM(Vname(1,vindex)),                    &
     &                        ng, TRIM(ncname), Fmin, Fmax
          END IF
        END IF
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            CLIMA(ng) % M2nudgcof)
        END IF
        CALL mp_exchange2d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      CLIMA(ng) % M2nudgcof)
      END IF
!
!  If appropriate, read nudging coefficients for 3D momentum (inverse
!  time scales, 1/day).
!
      IF (LnudgeM3CLM(ng)) THEN
        gtype=r3dvar
        vindex=idM3nc
        Fscl=1.0_dp/day2sec                    ! default units: 1/day
!
        CALL netcdf_inq_var (ng, model, ncname,                         &
     &                       ncid = NUD(ng)%ncid,                       &
     &                       MyVarName = TRIM(Vname(1,vindex)),         &
     &                       VarID = NUD(ng)%Vid(vindex),               &
     &                       nVarDim = nvdim,                           &
     &                       nVarAtt = nvatt)
        IF (FoundError(exit_flag, NoError, 268, MyFile)) RETURN
        DO i=1,nvatt
          IF (TRIM(var_Aname(i)).eq.'units') THEN
            tunits=TRIM(var_Achar(i))
            IF (tunits(1:3).eq.'day') THEN
              Fscl=1.0_dp/day2sec
            ELSE IF (tunits(1:6).eq.'second') THEN
              Fscl=1.0_dp
            END IF
          END IF
        END DO
!
        status=nf_fread3d(ng, model, ncname, NUD(ng)%ncid,              &
     &                    Vname(1,vindex), NUD(ng)%Vid(vindex),         &
     &                    0, gtype, Vsize,                              &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    Fscl, Fmin, Fmax,                             &
     &                    GRID(ng) % rmask,                             &
     &                    CLIMA(ng) % M3nudgcof)
        IF (FoundError(status, nf90_noerr, 294, MyFile)) THEN
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,vindex))//': '//             &
     &                        TRIM(Vname(1,vindex)),                    &
     &                        ng, TRIM(ncname), Fmin, Fmax
          END IF
        END IF
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            CLIMA(ng) % M3nudgcof)
        END IF
        CALL mp_exchange3d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      CLIMA(ng) % M3nudgcof)
      END IF
!
!  If appropriate, read nudging coefficients for tracer variables
!  (inverse time scales, 1/day). If the nudging coefficients for a
!  specific tracer are available, it will read that NetCDF variable.
!  If NOT AND the generic coefficients are available, it will process
!  those values instead.
!
      ic=0
      DO itrc=1,NT(ng)
        IF (LnudgeTCLM(itrc,ng)) THEN
          ic=ic+1
          gtype=r3dvar
          IF (got_specific(itrc)) THEN
            vindex=idTnud(itrc)                ! specific variable
            Fscl=1.0_dp/day2sec                ! default units: 1/day
          ELSE IF (got_generic) THEN
            vindex=idgTnc                      ! generic variable
            Fscl=1.0_dp/day2sec                ! default units: 1/day
          END IF
!
          CALL netcdf_inq_var (ng, model, ncname,                       &
     &                         ncid = NUD(ng)%ncid,                     &
     &                         MyVarName = TRIM(Vname(1,vindex)),       &
     &                         VarID = NUD(ng)%Vid(vindex),             &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (FoundError(exit_flag, NoError, 348, MyFile)) RETURN
          DO i=1,nvatt
            IF (TRIM(var_Aname(i)).eq.'units') THEN
              tunits=TRIM(var_Achar(i))
              IF (tunits(1:3).eq.'day') THEN
                Fscl=1.0_dp/day2sec
              ELSE IF (tunits(1:6).eq.'second') THEN
                Fscl=1.0_dp
              END IF
            END IF
          END DO
!
          status=nf_fread3d(ng, model, ncname, NUD(ng)%ncid,            &
     &                      Vname(1,vindex), NUD(ng)%Vid(vindex),       &
     &                      0, gtype, Vsize,                            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      Fscl, Fmin, Fmax,                           &
     &                      GRID(ng) % rmask,                           &
     &                      CLIMA(ng) % Tnudgcof(:,:,:,ic))
          IF (FoundError(status, nf90_noerr, 374, MyFile)) THEN
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,40) TRIM(Vname(2,vindex))//': '//           &
     &                          TRIM(Vname(1,vindex)),                  &
     &                          ng, TRIM(ncname), Fmin, Fmax
            END IF
          END IF
!
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              CLIMA(ng) % Tnudgcof(:,:,:,ic))
          END IF
        END IF
      END DO
      IF (ANY(LnudgeTCLM(:,ng))) THEN
        CALL mp_exchange4d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj, 1, N(ng), 1, NTCLM(ng), &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      CLIMA(ng) % Tnudgcof)
      END IF
!
  10  FORMAT (/,' GET_NUDGCOEF - unable to open nudging NetCDF file: ',a)
  20  FORMAT (/,' GET_NUDGCOEF - unable to find nudging variable: ',a,  &
     &        /,16x,'in NetCDF file: ',a)
  30  FORMAT (/,' GET_NUDGCOEF - unable to find nudging variable: ',a,  &
     &        /,16x,            '     or generc nudging variable: ',a,  &
     &        /,16x,'in nudging NetCDF file: ',a)
  40  FORMAT(2x,' GET_NUDGCOEF - ',a,/,19x,                             &
     &        '(Grid = ',i2.2,', File: ',a,')',/,19x,                   &
     &        '(Min = ', 1p,e15.8,0p,' Max = ',1p,e15.8,0p,')')
!
      RETURN
      END SUBROUTINE get_nudgcoef
