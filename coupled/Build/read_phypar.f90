      SUBROUTINE read_PhyPar (model, inp, out, Lwrite)
!
!svn $Id: read_phypar.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads and reports physical model input parameters.     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mct_coupler_params
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_sediment
      USE mod_stepping
      USE mod_strings
!
      USE inp_decode_mod
!
      USE dateclock_mod, ONLY : ref_clock
      USE strings_mod,   ONLY : FoundError
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: model, inp, out
!
!  Local variable declarations.
!
      logical :: got_Ngrids, got_NestLayers
      logical :: obc_data
      logical :: Lvalue(1)
      logical, allocatable :: Lswitch(:)
      logical, allocatable :: Lbottom(:,:)
      logical, allocatable :: Ltracer(:,:)
!
      integer :: Npts, Nval, i, itrc, ivar, k, lstr, ng, nl, status
      integer :: ifield, ifile, igrid, itracer, nline, max_Ffiles
      integer :: ibcfile, iclmfile
      integer :: Cdim, Clen, Rdim
      integer :: nPETs, maxPETs
      integer :: OutFiles
      integer :: Ivalue(1)
      integer, allocatable :: Nfiles(:)
      integer, allocatable :: Ncount(:,:)
      integer, allocatable :: NBCcount(:,:)
      integer, allocatable :: NCLMcount(:,:)
!
      real(dp), allocatable :: Dtracer(:,:)
      real(r8), allocatable :: Rtracer(:,:)
      real(r8), allocatable :: tracer(:,:)
      real(dp) :: Dvalue(1)
      real(r8) :: Rvalue(1)
      real(dp), dimension(nRval) :: Rval
!
      character (len=1  ), parameter :: blank = ' '
      character (len=40 ) :: KeyWord, text
      character (len=50 ) :: label
      character (len=256) :: fname, line
      character (len=256), dimension(nCval) :: Cval
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/read_phypar.F"
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
      ifile=1                            ! multiple file counter
      ibcfile=1                          ! multiple BC file counter
      iclmfile=1                         ! multiple CLM file counter
      igrid=1                            ! nested grid counter
      itracer=0                          ! LBC tracer counter
      nline=0                            ! LBC multi-line counter
      DO i=1,LEN(label)
        label(i:i)=blank
      END DO
      got_Ngrids=.FALSE.
      got_NestLayers=.FALSE.
      Cdim=SIZE(Cval,1)
      Clen=LEN(Cval(1))
      Rdim=SIZE(Rval,1)
      Nval=0
!
!-----------------------------------------------------------------------
!  Read in physical model parameters. Then, load input data into module.
!  Take into account nested grid configurations.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
            CASE ('TITLE')
              IF (Nval.eq.1) THEN
                title=TRIM(ADJUSTL(Cval(Nval)))
              ELSE
                WRITE(title,'(a,1x,a)') TRIM(ADJUSTL(title)),           &
     &                                  TRIM(ADJUSTL(Cval(Nval)))
              END IF
            CASE ('MyAppCPP')
              DO i=1,LEN(MyAppCPP)
                MyAppCPP(i:i)=blank
              END DO
              MyAppCPP=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('VARNAME')
              DO i=1,LEN(varname)
                varname(i:i)=blank
              END DO
              varname=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('Ngrids')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Ngrids=Ivalue(1)
              IF (Ngrids.le.0) THEN
                IF (Master) WRITE (out,290) 'Ngrids', Ngrids,           &
     &            'must be greater than zero.'
                exit_flag=5
                RETURN
              END IF
              got_Ngrids=.TRUE.
              CALL allocate_param       ! Start allocating variables in
              CALL allocate_parallel    ! modules that solely depend on
              CALL allocate_iounits     ! the number of nested grids
              CALL allocate_stepping
              IF (.not.allocated(Lswitch)) THEN
                allocate ( Lswitch(Ngrids) )
              END IF
              IF (.not.allocated(Lbottom)) THEN
                allocate ( Lbottom(MBOTP,Ngrids) )
              END IF
              IF (.not.allocated(Nfiles)) THEN
                allocate ( Nfiles(Ngrids) )
                Nfiles(1:Ngrids)=0
              END IF
            CASE ('NestLayers')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              NestLayers=Ivalue(1)
              IF (NestLayers.lt.1) THEN
                IF (Master) WRITE (out,290) 'NestLayers', NestLayers,   &
     &            'must be greater or equal than one.'
                exit_flag=5
                RETURN
              END IF
              IF (NestLayers.gt.1) THEN
                IF (Master) WRITE (out,290) 'NestLayers', NestLayers,   &
     &            'must be equal to one in non-nesting applications.'
                exit_flag=5
                RETURN
              END IF
              got_NestLayers=.TRUE.
              IF (.not.allocated(GridsInLayer)) THEN
                allocate ( GridsInLayer(NestLayers) )
                GridsInLayer(1:NestLayers)=1
              END IF
              IF (.not.allocated(GridNumber)) THEN
                allocate ( GridNumber(Ngrids,NestLayers) )
                GridNumber(1:Ngrids,1:NestLayers)=0        ! Important
              END IF
            CASE ('GridsInLayer')
              IF (.not.got_NestLayers) THEN
                IF (Master) WRITE (out,320) 'NestLayers',               &
     &            'Add "NestLayers" keyword before GridsInLayer.'
                exit_flag=5
                RETURN
              END IF
              Npts=load_i(Nval, Rval, NestLayers, GridsInLayer)
              ng=0
              DO nl=1,NestLayers
                DO i=1,GridsInLayer(nl)
                  ng=ng+1                  ! order of grids are very in
                  GridNumber(i,nl)=ng      ! nesting applications. See
                END DO                     ! WikiROMS for details.
              END DO
            CASE ('Lm')
              IF (.not.got_Ngrids) THEN
                IF (Master) WRITE (out,320) 'Ngrids',                   &
     &            'Add "Ngrids" keyword before grid dimension (Lm, Mm).'
                exit_flag=5
                RETURN
              END IF
              Npts=load_i(Nval, Rval, Ngrids, Lm)
              DO ng=1,Ngrids
                IF (Lm(ng).le.0) THEN
                  IF (Master) WRITE (out,300) 'Lm', ng,                 &
     &              'must be greater than zero.'
                  exit_flag=5
                  RETURN
                END IF
              END DO
            CASE ('Mm')
              Npts=load_i(Nval, Rval, Ngrids, Mm)
              DO ng=1,Ngrids
                IF (Mm(ng).le.0) THEN
                  IF (Master) WRITE (out,300) 'Mm', ng,                 &
     &              'must be greater than zero.'
                  exit_flag=5
                  RETURN
                END IF
              END DO
            CASE ('N')
              Npts=load_i(Nval, Rval, Ngrids, N)
              DO ng=1,Ngrids
                IF (N(ng).lt.0) THEN
                  IF (Master) WRITE (out,300) 'N', ng,                  &
     &              'must be greater than zero.'
                  exit_flag=5
                  RETURN
                END IF
              END DO
            CASE ('NAT')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              NAT=Ivalue(1)
              IF ((NAT.lt.1).or.(NAT.gt.2)) THEN
                IF (Master) WRITE (out,290) 'NAT = ', NAT,              &
     &            'make sure that NAT is either 1 or 2.'
                exit_flag=5
                RETURN
              END IF
            CASE ('NtileI')
              Npts=load_i(Nval, Rval, Ngrids, NtileI)
              NtileX(1:Ngrids)=NtileI(1:Ngrids)
            CASE ('NtileJ')
              Npts=load_i(Nval, Rval, Ngrids, NtileJ)
              NtileE(1:Ngrids)=NtileJ(1:Ngrids)
              CALL initialize_sediment
              CALL initialize_param    ! Continue allocating/initalizing
              CALL allocate_scalars    ! variables since the application
              CALL initialize_scalars  ! number of nested grids and
              CALL allocate_ncparam    ! domain parameters are known
              CALL initialize_ncparam
              IF (.not.allocated(Ltracer)) THEN
                allocate (Ltracer(NAT+NPT,Ngrids))
              END IF
              IF (.not.allocated(Dtracer)) THEN
                allocate (Dtracer(NAT+NPT,Ngrids))
              END IF
              IF (.not.allocated(Rtracer)) THEN
                allocate (Rtracer(NAT+NPT,Ngrids))
              END IF
              IF (.not.allocated(tracer)) THEN
                allocate (tracer(MT,Ngrids))
              END IF
            CASE ('Hadvection')
              IF (itracer.lt.(NAT+NPT)) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              itrc=itracer
              Npts=load_tadv(Nval, Cval, line, nline, itrc, igrid,      &
     &                       itracer, 1, NAT+NPT,                       &
     &                       Vname(1,idTvar(itrc)),                     &
     &                       Hadvection)
              IF (FoundError(exit_flag, NoError,                        &
     &                       455, MyFile)) RETURN
            CASE ('Vadvection')
              IF (itracer.lt.(NAT+NPT)) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              itrc=itracer
              Npts=load_tadv(Nval, Cval, line, nline, itrc, igrid,      &
     &                       itracer, 1, NAT+NPT,                       &
     &                       Vname(1,idTvar(itrc)),                     &
     &                       Vadvection)
              IF (FoundError(exit_flag, NoError,                        &
     &                       468, MyFile)) RETURN
            CASE ('LBC(isFsur)')
              Npts=load_lbc(Nval, Cval, line, nline, isFsur, igrid,     &
     &                      0, 0, Vname(1,idFsur), LBC)
            CASE ('LBC(isUbar)')
              Npts=load_lbc(Nval, Cval, line, nline, isUbar, igrid,     &
     &                      0, 0, Vname(1,idUbar), LBC)
            CASE ('LBC(isVbar)')
              Npts=load_lbc(Nval, Cval, line, nline, isVbar, igrid,     &
     &                      0, 0, Vname(1,idVbar), LBC)
            CASE ('LBC(isU2Sd)')
              Npts=load_lbc(Nval, Cval, line, nline, isU2Sd, igrid,     &
     &                      0, 0, Vname(1,idU2Sd), LBC)
            CASE ('LBC(isV2Sd)')
              Npts=load_lbc(Nval, Cval, line, nline, isV2Sd, igrid,     &
     &                      0, 0, Vname(1,idV2Sd), LBC)
            CASE ('LBC(isUvel)')
              Npts=load_lbc(Nval, Cval, line, nline, isUvel, igrid,     &
     &                      0, 0, Vname(1,idUvel), LBC)
            CASE ('LBC(isVvel)')
              Npts=load_lbc(Nval, Cval, line, nline, isVvel, igrid,     &
     &                      0, 0, Vname(1,idVvel), LBC)
            CASE ('LBC(isU3Sd)')
              Npts=load_lbc(Nval, Cval, line, nline, isU3Sd, igrid,     &
     &                      0, 0, Vname(1,idU3Sd), LBC)
            CASE ('LBC(isV3Sd)')
              Npts=load_lbc(Nval, Cval, line, nline, isV3Sd, igrid,     &
     &                      0, 0, Vname(1,idV3Sd), LBC)
            CASE ('LBC(isMtke)')
              Npts=load_lbc(Nval, Cval, line, nline, isMtke, igrid,     &
     &                      0, 0, Vname(1,idMtke), LBC)
            CASE ('LBC(isTvar)')
              IF (itracer.lt.(NAT+NPT)) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              ifield=isTvar(itracer)
              Npts=load_lbc(Nval, Cval, line, nline, ifield, igrid,     &
     &                      1, NAT+NPT, Vname(1,idTvar(itracer)), LBC)
            CASE ('VolCons(west)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              VolCons(iwest,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('VolCons(east)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              VolCons(ieast,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('VolCons(south)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              VolCons(isouth,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('VolCons(north)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              VolCons(inorth,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('NTIMES')
              Npts=load_i(Nval, Rval, Ngrids, ntimes)
            CASE ('DT')
              Npts=load_r(Nval, Rval, Ngrids, dt)
            CASE ('NDTFAST')
              Npts=load_i(Nval, Rval, Ngrids, ndtfast)
            CASE ('ERstr')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              ERstr=Ivalue(1)
            CASE ('ERend')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              ERend=Ivalue(1)
            CASE ('Nouter')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Nouter=Ivalue(1)
            CASE ('Ninner')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Ninner=Ivalue(1)
            CASE ('Nsaddle')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Nsaddle=Ivalue(1)
            CASE ('Nintervals')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Nintervals=Ivalue(1)
            CASE ('NRREC')
              Npts=load_i(Nval, Rval, Ngrids, nrrec)
              DO ng=1,Ngrids
                IF (nrrec(ng).lt.0) THEN
                  LastRec(ng)=.TRUE.
                ELSE
                  LastRec(ng)=.FALSE.
                END IF
              END DO
            CASE ('LcycleRST')
              Npts=load_l(Nval, Cval, Ngrids, LcycleRST)
            CASE ('NRST')
              Npts=load_i(Nval, Rval, Ngrids, nRST)
            CASE ('NSTA')
              Npts=load_i(Nval, Rval, Ngrids, nSTA)
            CASE ('NFLT')
              Npts=load_i(Nval, Rval, Ngrids, nFLT)
            CASE ('NINFO')
              Npts=load_i(Nval, Rval, Ngrids, ninfo)
              DO ng=1,Ngrids
                IF (ninfo(ng).le.0) THEN
                  WRITE (text,'(a,i2.2,a)') 'ninfo(', ng, ') = '
                  IF (Master) WRITE (out,260) TRIM(text), ninfo(ng),    &
     &               'must be greater than zero.'
                  exit_flag=5
                  RETURN
                END IF
              END DO
            CASE ('LDEFOUT')
              Npts=load_l(Nval, Cval, Ngrids, ldefout)
            CASE ('NHIS')
              Npts=load_i(Nval, Rval, Ngrids, nHIS)
            CASE ('NDEFHIS')
              Npts=load_i(Nval, Rval, Ngrids, ndefHIS)
            CASE ('NQCK')
              Npts=load_i(Nval, Rval, Ngrids, nQCK)
            CASE ('NDEFQCK')
              Npts=load_i(Nval, Rval, Ngrids, ndefQCK)
            CASE ('NTSAVG')
              Npts=load_i(Nval, Rval, Ngrids, ntsAVG)
            CASE ('NAVG')
              Npts=load_i(Nval, Rval, Ngrids, nAVG)
            CASE ('NDEFAVG')
              Npts=load_i(Nval, Rval, Ngrids, ndefAVG)
            CASE ('NTSDIA')
              Npts=load_i(Nval, Rval, Ngrids, ntsDIA)
            CASE ('NDIA')
              Npts=load_i(Nval, Rval, Ngrids, nDIA)
            CASE ('NDEFDIA')
              Npts=load_i(Nval, Rval, Ngrids, ndefDIA)
            CASE ('LcycleTLM')
              Npts=load_l(Nval, Cval, Ngrids, LcycleTLM)
            CASE ('NTLM')
              Npts=load_i(Nval, Rval, Ngrids, nTLM)
            CASE ('NDEFTLM')
              Npts=load_i(Nval, Rval, Ngrids, ndefTLM)
            CASE ('LcycleADJ')
              Npts=load_l(Nval, Cval, Ngrids, LcycleADJ)
            CASE ('NADJ')
              Npts=load_i(Nval, Rval, Ngrids, nADJ)
            CASE ('NDEFADJ')
              Npts=load_i(Nval, Rval, Ngrids, ndefADJ)
            CASE ('NOBC')
              Npts=load_i(Nval, Rval, Ngrids, nOBC)
            CASE ('NSFF')
              Npts=load_i(Nval, Rval, Ngrids, nSFF)
            CASE ('LmultiGST')
              Npts=load_l(Nval, Cval, 1, Lvalue)
              LmultiGST=Lvalue(1)
            CASE ('LrstGST')
              Npts=load_l(Nval, Cval, 1, Lvalue)
              LrstGST=Lvalue(1)
            CASE ('MaxIterGST')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              MaxIterGST=Ivalue(1)
            CASE ('NGST')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              nGST=Ivalue(1)
            CASE ('TNU2')
              Npts=load_r(Nval, Rval, NAT+NPT, Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  nl_tnu2(itrc,ng)=Rtracer(itrc,ng)
                END DO
              END DO
            CASE ('TNU4')
              Npts=load_r(Nval, Rval, NAT+NPT, Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  nl_tnu4(itrc,ng)=Rtracer(itrc,ng)
                END DO
              END DO
            CASE ('ad_TNU2')
              Npts=load_r(Nval, Rval, NAT+NPT, Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  ad_tnu2(itrc,ng)=Rtracer(itrc,ng)
                  tl_tnu2(itrc,ng)=Rtracer(itrc,ng)
                END DO
              END DO
            CASE ('ad_TNU4')
              Npts=load_r(Nval, Rval, NAT+NPT, Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  ad_tnu4(itrc,ng)=Rtracer(itrc,ng)
                  tl_tnu4(itrc,ng)=Rtracer(itrc,ng)
                END DO
              END DO
            CASE ('VISC2')
              Npts=load_r(Nval, Rval, Ngrids, nl_visc2)
            CASE ('VISC4')
              Npts=load_r(Nval, Rval, Ngrids, nl_visc4)
            CASE ('ad_VISC2')
              Npts=load_r(Nval, Rval, Ngrids, ad_visc2)
              DO ng=1,Ngrids
                tl_visc2(ng)=ad_visc2(ng)
              END DO
            CASE ('ad_VISC4')
              Npts=load_r(Nval, Rval, Ngrids, ad_visc4)
              DO ng=1,Ngrids
                tl_visc4(ng)=ad_visc4(ng)
              END DO
            CASE ('LuvSponge')
              Npts=load_l(Nval, Cval, Ngrids, LuvSponge)
            CASE ('LtracerSponge')
              Npts=load_l(Nval, Cval, NAT+NPT, Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  LtracerSponge(itrc,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('AKT_BAK')
              Npts=load_r(Nval, Rval, NAT+NPT, Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Akt_bak(itrc,ng)=Rtracer(itrc,ng)
                END DO
              END DO
            CASE ('AKT_LIMIT')
              Npts=load_r(Nval, Rval, NAT+NPT, Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Akt_limit(itrc,ng)=Rtracer(itrc,ng)
                END DO
              END DO
            CASE ('ad_AKT_fac')
              Npts=load_r(Nval, Rval, NAT+NPT, Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  ad_Akt_fac(itrc,ng)=Rtracer(itrc,ng)
                  tl_Akt_fac(itrc,ng)=Rtracer(itrc,ng)
                END DO
              END DO
            CASE ('AKV_BAK')
              Npts=load_r(Nval, Rval, Ngrids, Akv_bak)
            CASE ('AKV_LIMIT')
              Npts=load_r(Nval, Rval, Ngrids, Akv_limit)
            CASE ('ad_AKV_fac')
              Npts=load_r(Nval, Rval, Ngrids, ad_Akv_fac)
              DO ng=1,Ngrids
                tl_Akv_fac(ng)=ad_AKv_fac(ng)
              END DO
            CASE ('AKK_BAK')
              Npts=load_r(Nval, Rval, Ngrids, Akk_bak)
            CASE ('AKP_BAK')
              Npts=load_r(Nval, Rval, Ngrids, Akp_bak)
            CASE ('TKENU2')
              Npts=load_r(Nval, Rval, Ngrids, tkenu2)
            CASE ('TKENU4')
              Npts=load_r(Nval, Rval, Ngrids, tkenu4)
            CASE ('GLS_P')
              Npts=load_r(Nval, Rval, Ngrids, gls_p)
            CASE ('GLS_M')
              Npts=load_r(Nval, Rval, Ngrids, gls_m)
            CASE ('GLS_N')
              Npts=load_r(Nval, Rval, Ngrids, gls_n)
            CASE ('GLS_Kmin')
              Npts=load_r(Nval, Rval, Ngrids, gls_Kmin)
            CASE ('GLS_Pmin')
              Npts=load_r(Nval, Rval, Ngrids, gls_Pmin)
            CASE ('GLS_CMU0')
              Npts=load_r(Nval, Rval, Ngrids, gls_cmu0)
            CASE ('GLS_C1')
              Npts=load_r(Nval, Rval, Ngrids, gls_c1)
            CASE ('GLS_C2')
              Npts=load_r(Nval, Rval, Ngrids, gls_c2)
            CASE ('GLS_C3M')
              Npts=load_r(Nval, Rval, Ngrids, gls_c3m)
            CASE ('GLS_C3P')
              Npts=load_r(Nval, Rval, Ngrids, gls_c3p)
            CASE ('GLS_SIGK')
              Npts=load_r(Nval, Rval, Ngrids, gls_sigk)
            CASE ('GLS_SIGP')
              Npts=load_r(Nval, Rval, Ngrids, gls_sigp)
            CASE ('CHARNOK_ALPHA')
              Npts=load_r(Nval, Rval, Ngrids, charnok_alpha)
            CASE ('ZOS_HSIG_ALPHA')
              Npts=load_r(Nval, Rval, Ngrids, zos_hsig_alpha)
            CASE ('SZ_ALPHA')
              Npts=load_r(Nval, Rval, Ngrids, sz_alpha)
            CASE ('CRGBAN_CW')
              Npts=load_r(Nval, Rval, Ngrids, crgban_cw)
            CASE ('WEC_ALPHA')
              Npts=load_r(Nval, Rval, Ngrids, wec_alpha)
            CASE ('RDRG')
              Npts=load_r(Nval, Rval, Ngrids, rdrg)
            CASE ('RDRG2')
              Npts=load_r(Nval, Rval, Ngrids, rdrg2)
            CASE ('Zob')
              Npts=load_r(Nval, Rval, Ngrids, Zob)
            CASE ('Zos')
              Npts=load_r(Nval, Rval, Ngrids, Zos)
            CASE ('DCRIT')
              Npts=load_r(Nval, Rval, Ngrids, Dcrit)
            CASE ('WTYPE')
              Npts=load_i(Nval, Rval, Ngrids, lmd_Jwt)
            CASE ('LEVSFRC')
              Npts=load_i(Nval, Rval, Ngrids, levsfrc)
            CASE ('LEVBFRC')
              Npts=load_i(Nval, Rval, Ngrids, levbfrc)
            CASE ('Vtransform')
              Npts=load_i(Nval, Rval, Ngrids, Vtransform)
              DO ng=1,Ngrids
                IF ((Vtransform(ng).lt.0).or.                           &
     &              (Vtransform(ng).gt.2)) THEN
                  IF (Master) WRITE (out,260) 'Vtransform = ',          &
     &                                        Vtransform(ng),           &
     &                                        'Must be either 1 or 2'
                  exit_flag=5
                  RETURN
                END IF
              END DO
            CASE ('Vstretching')
              Npts=load_i(Nval, Rval, Ngrids, Vstretching)
              DO ng=1,Ngrids
                IF ((Vstretching(ng).lt.0).or.                          &
     &              (Vstretching(ng).gt.5)) THEN
                  IF (Master) WRITE (out,260) 'Vstretching = ',         &
     &                                        Vstretching(ng),          &
     &                                        'Must between 1 and 5'
                  exit_flag=5
                  RETURN
                END IF
              END DO
            CASE ('THETA_S')
              Npts=load_r(Nval, Rval, Ngrids, theta_s)
            CASE ('THETA_B')
              Npts=load_r(Nval, Rval, Ngrids, theta_b)
            CASE ('TCLINE')
              Npts=load_r(Nval, Rval, Ngrids, Tcline)
              DO ng=1,Ngrids
                hc(ng)=Tcline(ng)
              END DO
            CASE ('RHO0')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              rho0=Rvalue(1)
            CASE ('BVF_BAK')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              bvf_bak=Rvalue(1)
            CASE ('DSTART')
              Npts=load_r(Nval, Rval, 1, Dvalue)
              dstart=Dvalue(1)
            CASE ('TIDE_START')
              Npts=load_r(Nval, Rval, 1, Dvalue)
              tide_start=Dvalue(1)
            CASE ('TIME_REF')
              Npts=load_r(Nval, Rval, 1, Dvalue)
              time_ref=Dvalue(1)
              CALL ref_clock (time_ref)
            CASE ('TNUDG')
              Npts=load_r(Nval, Rval, NAT+NPT, Ngrids, Dtracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Tnudg(itrc,ng)=Dtracer(itrc,ng)
                END DO
              END DO
            CASE ('ZNUDG')
              Npts=load_r(Nval, Rval, Ngrids, Znudg)
            CASE ('M2NUDG')
              Npts=load_r(Nval, Rval, Ngrids, M2nudg)
            CASE ('M3NUDG')
              Npts=load_r(Nval, Rval, Ngrids, M3nudg)
            CASE ('TNUDG_SSS')
              Npts=load_r(Nval, Rval, Ngrids, Tnudg_SSS)
            CASE ('SSS_MISMATCH_THRESHOLD')
              Npts=load_r(Nval, Rval, Ngrids, sss_mismatch_threshold)
            CASE ('OBCFAC')
              Npts=load_r(Nval, Rval, Ngrids, obcfac)
            CASE ('R0')
              Npts=load_r(Nval, Rval, Ngrids, R0)
              DO ng=1,Ngrids
                IF (R0(ng).lt.100.0_r8) R0(ng)=R0(ng)+1000.0_r8
              END DO
            CASE ('T0')
              Npts=load_r(Nval, Rval, Ngrids, T0)
            CASE ('S0')
              Npts=load_r(Nval, Rval, Ngrids, S0)
            CASE ('TCOEF')
              Npts=load_r(Nval, Rval, Ngrids, Tcoef)
              DO ng=1,Ngrids
                Tcoef(ng)=ABS(Tcoef(ng))
              END DO
            CASE ('SCOEF')
              Npts=load_r(Nval, Rval, Ngrids, Scoef)
              DO ng=1,Ngrids
                Scoef(ng)=ABS(Scoef(ng))
              END DO
            CASE ('GAMMA2')
              Npts=load_r(Nval, Rval, Ngrids, gamma2)
            CASE ('LuvSrc')
              Npts=load_l(Nval, Cval, Ngrids, LuvSrc)
            CASE ('LwSrc')
              Npts=load_l(Nval, Cval, Ngrids, LwSrc)
            CASE ('LtracerSrc')
              Npts=load_l(Nval, Cval, NAT+NPT, Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  LtracerSrc(itrc,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('LsshCLM')
              Npts=load_l(Nval, Cval, Ngrids, LsshCLM)
            CASE ('Lm2CLM')
              Npts=load_l(Nval, Cval, Ngrids, Lm2CLM)
            CASE ('Lm3CLM')
              Npts=load_l(Nval, Cval, Ngrids, Lm3CLM)
            CASE ('LtracerCLM')
              Npts=load_l(Nval, Cval, NAT+NPT, Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  LtracerCLM(itrc,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('LnudgeM2CLM')
              Npts=load_l(Nval, Cval, Ngrids, LnudgeM2CLM)
            CASE ('LnudgeM3CLM')
              Npts=load_l(Nval, Cval, Ngrids, LnudgeM3CLM)
            CASE ('LnudgeTCLM')
              Npts=load_l(Nval, Cval, NAT+NPT, Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  LnudgeTCLM(itrc,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idUvel)')
              IF (idUvel.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUvel'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idUvel,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idVvel)')
              IF (idVvel.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVvel'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idVvel,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idWvel)')
              IF (idWvel.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWvel'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idWvel,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idOvel)')
              IF (idOvel.eq.0) THEN
                IF (Master) WRITE (out,280) 'idOvel'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idOvel,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idUbar)')
              IF (idUbar.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUbar'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idUbar,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idVbar)')
              IF (idVbar.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVbar'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idVbar,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idFsur)')
              IF (idFsur.eq.0) THEN
                IF (Master) WRITE (out,280) 'idFsur'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idFsur,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idu2dE)')
              IF (idu2dE.eq.0) THEN
                IF (Master) WRITE (out,280) 'idu2dE'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idu2dE,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idv2dN)')
              IF (idv2dN.eq.0) THEN
                IF (Master) WRITE (out,280) 'idv2dN'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idv2dN,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idu3dE)')
              IF (idu3dE.eq.0) THEN
                IF (Master) WRITE (out,280) 'idu3dE'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idu3dE,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idv3dN)')
              IF (idv3dN.eq.0) THEN
                IF (Master) WRITE (out,280) 'idv3dN'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idv3dN,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idTvar)')
              IF (MAXVAL(idTvar).eq.0) THEN
                IF (Master) WRITE (out,280) 'idTvar'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NAT+NPT, Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idTvar(itrc)
                  Hout(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idpthR)')
              IF (idpthR.eq.0) THEN
                IF (Master) WRITE (out,280) 'idpthR'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idpthR,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idpthU)')
              IF (idpthU.eq.0) THEN
                IF (Master) WRITE (out,280) 'idpthU'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idpthU,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idpthV)')
              IF (idpthV.eq.0) THEN
                IF (Master) WRITE (out,280) 'idpthV'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idpthV,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idpthW)')
              IF (idpthW.eq.0) THEN
                IF (Master) WRITE (out,280) 'idpthW'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idpthW,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idUsms)')
              IF (idUsms.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUsms'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idUsms,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idVsms)')
              IF (idVsms.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVsms'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idVsms,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idUbms)')
              IF (idUbms.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUbms'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idUbms,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idVbms)')
              IF (idVbms.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVbms'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idVbms,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idUbrs)')
              IF (idUbrs.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUbrs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idUbrs,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idVbrs)')
              IF (idVbrs.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVbrs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idVbrs,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idUbws)')
              IF (idUbws.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUbws'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idUbws,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idVbws)')
              IF (idVbws.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVbws'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idVbws,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idUbcs)')
              IF (idUbcs.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUbcs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idUbcs,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idVbcs)')
              IF (idVbcs.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVbcs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idVbcs,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idUVwc)')
              IF (idUVwc.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUVwc'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idUVwc,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idUbot)')
              IF (idUbot.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUbot'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idUbot,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idVbot)')
              IF (idVbot.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVbot'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idVbot,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idUbur)')
              IF (idUbur.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUbur'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idUbur,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idVbvr)')
              IF (idVbvr.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVbvr'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idVbvr,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idWztw)')
              IF (idWztw.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWztw'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idWztw,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idWqsp)')
              IF (idWqsp.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWqsp'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idWqsp,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idWbeh)')
              IF (idWbeh.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWbeh'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idWbeh,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idU2rs)')
              IF (idU2rs.eq.0) THEN
                IF (Master) WRITE (out,280) 'idU2rs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idU2rs,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idV2rs)')
              IF (idV2rs.eq.0) THEN
                IF (Master) WRITE (out,280) 'idV2rs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idV2rs,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idU3rs)')
              IF (idU3rs.eq.0) THEN
                IF (Master) WRITE (out,280) 'idU3rs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idU3rs,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idV3rs)')
              IF (idV3rs.eq.0) THEN
                IF (Master) WRITE (out,280) 'idV3rs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idV3rs,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idU2Sd)')
              IF (idU2Sd.eq.0) THEN
                IF (Master) WRITE (out,280) 'idU2Sd'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idU2Sd,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idV2Sd)')
              IF (idV2Sd.eq.0) THEN
                IF (Master) WRITE (out,280) 'idV2Sd'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idV2Sd,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idU3Sd)')
              IF (idU3Sd.eq.0) THEN
                IF (Master) WRITE (out,280) 'idU3Sd'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idU3Sd,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idV3Sd)')
              IF (idV3Sd.eq.0) THEN
                IF (Master) WRITE (out,280) 'idV3Sd'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idV3Sd,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idW3Sd)')
              IF (idW3Sd.eq.0) THEN
                IF (Master) WRITE (out,280) 'idW3Sd'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idW3Sd,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idW3St)')
              IF (idW3St.eq.0) THEN
                IF (Master) WRITE (out,280) 'idW3St'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idW3St,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idWamp)')
              IF (idWamp.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWamp'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idWamp,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idWlen)')
              IF (idWlen.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWlen'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idWlen,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idWdir)')
              IF (idWdir.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWdir'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idWdir,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idWdip)')
              IF (idWdip.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWdip'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idWdip,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idWptp)')
              IF (idWptp.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWptp'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idWptp,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idWpbt)')
              IF (idWpbt.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWpbt'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idWpbt,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idWorb)')
              IF (idWorb.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWorb'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idWorb,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idUwav)')
              IF (idUwav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUwav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idUwav,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idVwav)')
              IF (idVwav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVwav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idVwav,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idWdif)')
              IF (idWdif.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWdif'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idWdif,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idWdib)')
              IF (idWdib.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWdib'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idWdib,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idWdiw)')
              IF (idWdiw.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWdiw'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idWdiw,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idWdis)')
              IF (idWdis.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWdis'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idWdis,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idWrol)')
              IF (idWrol.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWrol'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idWrol,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idTsur)')
              IF (idTsur(itemp).eq.0) THEN
                IF (Master) WRITE (out,280) 'idTsur'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NAT+NPT, Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idTsur(itrc)
                  Hout(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idLhea)')
              IF (idLhea.eq.0) THEN
                IF (Master) WRITE (out,280) 'idLhea'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idLhea,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idShea)')
              IF (idShea.eq.0) THEN
                IF (Master) WRITE (out,280) 'idShea'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idShea,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idLrad)')
              IF (idLrad.eq.0) THEN
                IF (Master) WRITE (out,280) 'idLrad'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idLrad,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idSrad)')
              IF (idSrad.eq.0) THEN
                IF (Master) WRITE (out,280) 'idSrad'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idSrad,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idEmPf)')
              IF (idEmPf.eq.0) THEN
                IF (Master) WRITE (out,280) 'idEmPf'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idEmPf,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idevap)')
              IF (idevap.eq.0) THEN
                IF (Master) WRITE (out,280) 'idevap'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idevap,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idrain)')
              IF (idrain.eq.0) THEN
                IF (Master) WRITE (out,280) 'idrain'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idrain,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idDano)')
              IF (idDano.eq.0) THEN
                IF (Master) WRITE (out,280) 'idDano'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idDano,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idVvis)')
              IF (idVvis.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVvis'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idVvis,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idTdif)')
              IF (idTdif.eq.0) THEN
                IF (Master) WRITE (out,280) 'idTdif'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idTdif,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idHsbl)')
              IF (idHsbl.eq.0) THEN
                IF (Master) WRITE (out,280) 'idHsbl'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idHsbl,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idHbbl)')
              IF (idHbbl.eq.0) THEN
                IF (Master) WRITE (out,280) 'idHbbl'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idHbbl,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idMtke)')
              IF (idMtke.eq.0) THEN
                IF (Master) WRITE (out,280) 'idMtke'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idMtke,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idMtls)')
              IF (idMtls.eq.0) THEN
                IF (Master) WRITE (out,280) 'idMtls'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idMtls,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idUvel)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idUvel,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idVvel)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idVvel,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idWvel)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idWvel,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idOvel)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idOvel,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idUbar)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idUbar,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idVbar)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idVbar,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idFsur)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idFsur,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idu2dE)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idu2dE,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idv2dN)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idv2dN,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idu3dE)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idu3dE,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idv3dN)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idv3dN,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idTvar)')
              Npts=load_l(Nval, Cval, NAT+NPT, Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idTvar(itrc)
                  Qout(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Qout(idUsur)')
              IF (idUsur.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUsur'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idUsur,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idVsur)')
              IF (idUsur.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVsur'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idVsur,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idUsuE)')
              IF (idUsuE.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUsuE'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idUsuE,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idVsuN)')
              IF (idVsuN.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVsuN'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idVsuN,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idsurT)')
              IF (MAXVAL(idsurT).eq.0) THEN
                IF (Master) WRITE (out,280) 'idsurT'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NAT+NPT, Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idsurT(itrc)
                  Qout(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Qout(idpthR)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idpthR,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idpthU)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idpthU,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idpthV)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idpthV,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idpthW)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idpthW,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idUsms)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idUsms,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idVsms)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idVsms,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idUbms)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idUbms,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idVbms)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idVbms,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idUbrs)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idUbrs,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idVbrs)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idVbrs,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idUbws)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idUbws,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idVbws)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idVbws,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idUbcs)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idUbcs,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idVbcs)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idVbcs,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idUbot)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idUbot,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idVbot)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idVbot,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idUbur)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idUbur,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idVbvr)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idVbvr,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idW2xx)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idW2xx,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idW2xy)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idW2xy,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idW2yy)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idW2yy,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idU2rs)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idU2rs,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idV2rs)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idV2rs,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idU2Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idU2Sd,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idV2Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idV2Sd,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idW3xx)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idW3xx,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idW3xy)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idW3xy,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idW3yy)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idW3yy,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idW3zx)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idW3zx,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idW3zy)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idW3zy,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idU3rs)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idU3rs,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idV3rs)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idV3rs,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idU3Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idU3Sd,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idV3Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idV3Sd,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idWamp)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idWamp,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idWlen)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idWlen,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idWdir)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idWdir,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idWdip)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idWdip,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idWptp)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idWptp,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idWpbt)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idWpbt,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idWorb)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idWorb,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idWdis)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idWdis,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idTsur)')
              Npts=load_l(Nval, Cval, NAT+NPT, Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idTsur(itrc)
                  Qout(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Qout(idLhea)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idLhea,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idShea)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idShea,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idLrad)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idLrad,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idSrad)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idSrad,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idEmPf)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idEmPf,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idevap)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idevap,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idrain)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idrain,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idDano)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idDano,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idVvis)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idVvis,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idTdif)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idTdif,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idHsbl)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idHsbl,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idHbbl)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idHbbl,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idMtke)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idMtke,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idMtls)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idMtls,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('NUSER')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Nuser=Ivalue(1)
            CASE ('USER')
              Npts=load_r(Nval, Rval, MAX(1,Nuser), user)
            CASE ('NC_SHUFFLE')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              shuffle=Ivalue(1)
            CASE ('NC_DEFLATE')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              deflate=Ivalue(1)
            CASE ('NC_DLEVEL')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              deflate_level=Ivalue(1)
            CASE ('DAINAME')
              label='DAI - Data Assimilation Initial/Restart fields'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, DAI)
            CASE ('GSTNAME')
              label='GST - generalized stability theory analysis'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, GST)
            CASE ('RSTNAME')
              label='RST - restart fields'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, RST)
            CASE ('HISNAME')
              label='HIS - nonlinear model history fields'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, HIS)
            CASE ('QCKNAME')
              label='QCK - nonlinear model quicksave fields'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, QCK)
            CASE ('TLMNAME')
              label='TLM - tangent linear model history fields'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, TLM)
            CASE ('TLFNAME')
              label='TLF - tangent linear model impulse forcing'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, TLF)
            CASE ('ADJNAME')
              label='ADM - adjoint model history fields'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, ADM)
            CASE ('AVGNAME')
              label='AVG - time-averaged history fields'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, AVG)
            CASE ('HARNAME')
              label='HAR - least-squares detiding harmonics'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, HAR)
            CASE ('DIANAME')
              label='DIA - time-averaged diagnostics fields'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, DIA)
            CASE ('STANAME')
              label='STA - stations time-series'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, STA)
            CASE ('FLTNAME')
              label='FLT - Lagragian particles trajectories'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, FLT)
            CASE ('GRDNAME')
              label='GRD - application grid'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, GRD)
            CASE ('ININAME')
              label='INI - nonlinear model initial conditions'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, INI)
            CASE ('IRPNAME')
              label='IRP - representer model initial conditions'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, IRP)
            CASE ('ITLNAME')
              label='ITL - tangent linear model initial conditions'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, ITL)
            CASE ('IADNAME')
              label='IAD - adjoint model initial conditions'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, IAD)
            CASE ('FWDNAME')
              label='FWD - basic state forward trajectory'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, FWD)
            CASE ('ADSNAME')
              label='ADS - adjoint sensitivity functional'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, ADS)
            CASE ('NGCNAME')
              DO i=1,LEN(NGCname)
                NGCname(i:i)=blank
              END DO
              NGCname=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('NBCFILES')
              Npts=load_i(Nval, Rval, Ngrids, nBCfiles)
              DO ng=1,Ngrids
                IF (nBCfiles(ng).le.0) THEN
                  IF (Master) WRITE (out,260) 'NBCFILES', nBCfiles(ng), &
     &                            'Must be equal or greater than one.'
                  exit_flag=4
                  RETURN
                END IF
              END DO
              max_Ffiles=MAXVAL(nBCfiles)
              allocate ( BRY(max_Ffiles,Ngrids) )
              allocate ( BRYids(max_Ffiles,Ngrids) )
              allocate ( NBCcount(max_Ffiles,Ngrids) )
              BRYids(1:max_Ffiles,1:Ngrids)=-1
              NBCcount(1:max_Ffiles,1:Ngrids)=0
            CASE ('BRYNAME')
              label='BRY - lateral open boundary conditions'
              DO ng=1,Ngrids
                IF (nBCfiles(ng).lt.0) THEN
                  IF (Master) WRITE (out,290) 'nBCfiles = ',            &
     &                                        nBCfiles(ng),             &
     &              'KeyWord ''NBCFILES'' unread or misssing from '//   &
     &              'input script ''roms.in''.'
                  exit_flag=4
                  RETURN
                END IF
              END DO
              Npts=load_s2d(Nval, Cval, Cdim, line, label, ibcfile,     &
     &                      igrid, Ngrids, nBCfiles, NBCcount,          &
     &                      max_Ffiles, BRY)
            CASE ('NCLMFILES')
              Npts=load_i(Nval, Rval, Ngrids, nCLMfiles)
              DO ng=1,Ngrids
                IF (nCLMfiles(ng).le.0) THEN
                  IF (Master) WRITE (out,260) 'NCLMFILES',              &
     &               nCLMfiles(ng),                                     &
     &              'Must be equal or greater than one.'
                  exit_flag=4
                  RETURN
                END IF
              END DO
              max_Ffiles=MAXVAL(nCLMfiles)
              allocate ( CLM(max_Ffiles,Ngrids) )
              allocate ( CLMids(max_Ffiles,Ngrids) )
              allocate ( NCLMcount(max_Ffiles,Ngrids) )
              CLMids(1:max_Ffiles,1:Ngrids)=-1
              NCLMcount(1:max_Ffiles,1:Ngrids)=0
            CASE ('CLMNAME')
              label='CLM - climatology fields'
              DO ng=1,Ngrids
                IF (nCLMfiles(ng).lt.0) THEN
                  IF (Master) WRITE (out,290) 'nCLMfiles = ',           &
     &                                        nCLMfiles(ng),            &
     &              'KeyWord ''NCLMFILES'' unread or misssing from '//  &
     &              'input script ''roms.in''.'
                  exit_flag=4
                  RETURN
                END IF
              END DO
              Npts=load_s2d(Nval, Cval, Cdim, line, label, iclmfile,    &
     &                      igrid, Ngrids, nCLMfiles, NCLMcount,        &
     &                      max_Ffiles, CLM)
            CASE ('NUDNAME')
              label='NUD - nudging coefficients'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, NUD)
            CASE ('SSFNAME')
              label='SSF - Sources/Sinks forcing fields'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, SSF)
            CASE ('NFFILES')
              Npts=load_i(Nval, Rval, Ngrids, nFfiles)
              DO ng=1,Ngrids
                IF (nFfiles(ng).le.0) THEN
                  IF (Master) WRITE (out,260) 'NFFILES', nFfiles(ng),   &
     &              'Must be equal or greater than one.'
                  exit_flag=4
                  RETURN
                END IF
              END DO
              max_Ffiles=MAXVAL(nFfiles)
              allocate ( FRC(max_Ffiles,Ngrids) )
              allocate ( FRCids(max_Ffiles,Ngrids) )
              allocate ( Ncount(max_Ffiles,Ngrids) )
              FRCids(1:max_Ffiles,1:Ngrids)=-1
              Ncount(1:max_Ffiles,1:Ngrids)=0
            CASE ('FRCNAME')
              label='FRC - forcing fields'
              DO ng=1,Ngrids
                IF (nFfiles(ng).lt.0) THEN
                  IF (Master) WRITE (out,290) 'nFfiles = ',             &
     &                                        nFfiles(ng),              &
     &              'KeyWord ''NFFILES'' unread or misssing from '//    &
     &              'input script ''roms.in''.'
                  exit_flag=4
                  RETURN
                END IF
              END DO
              Npts=load_s2d(Nval, Cval, Cdim, line, label, ifile,       &
     &                      igrid, Ngrids, nFfiles, Ncount, max_Ffiles, &
     &                      FRC)
            CASE ('APARNAM')
              DO i=1,LEN(aparnam)
                aparnam(i:i)=blank
              END DO
              aparnam=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('SPOSNAM')
              DO i=1,LEN(sposnam)
                sposnam(i:i)=blank
              END DO
              sposnam=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('FPOSNAM')
              DO i=1,LEN(fposnam)
                fposnam(i:i)=blank
              END DO
              fposnam=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('BPARNAM')
              DO i=1,LEN(bparnam)
                bparnam(i:i)=blank
              END DO
              bparnam=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('SPARNAM')
              DO i=1,LEN(sparnam)
                sparnam(i:i)=blank
              END DO
              sparnam=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('VEGNAM')
              DO i=1,LEN(vegnam)
                vegnam(i:i)=blank
              END DO
              vegnam=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('USRNAME')
              DO i=1,LEN(USRname)
                USRname(i:i)=blank
              END DO
              USRname=TRIM(ADJUSTL(Cval(Nval)))
          END SELECT
          IF (FoundError(exit_flag, NoError, 4169, MyFile)) RETURN
        END IF
      END DO
  10  IF (Master) WRITE (out,50) line
      exit_flag=4
      RETURN
  20  CLOSE (inp)
!
!-----------------------------------------------------------------------
!  Process input parameters.
!-----------------------------------------------------------------------
!
!  Check if nesting parameters "NestLayers", "GridsInLayer", and
!  "GridNumber" have been assigned.  The code below is necessary
!  for compatability with old "roms.in" input scripts.
!
      IF (.not.got_NestLayers) THEN
        NestLayers=1
        IF (.not.allocated(GridsInLayer)) THEN
          allocate ( GridsInLayer(NestLayers) )
        END IF
        IF (.not.allocated(GridNumber)) THEN
          allocate ( GridNumber(Ngrids,NestLayers) )
        END IF
      END IF
      GridsInLayer=1              ! In case that users set illegal
      GridNumber=1                ! values in non-nesting applications
!
!  Check if both point sources methodologies are activated.  Only one
!  method is allowed for a particular grid.  Otherwise, the point
!  source will be applies twice.
!
      DO ng=1,Ngrids
        IF (LuvSrc(ng).and.LwSrc(ng)) THEN
          IF (Master) THEN
            WRITE (out,260) 'LuvSrc', LuvSrc(ng),                       &
     &            'Because LwSrc  is also ON; only one method is legal.'
            WRITE (out,260) 'LwSrc', LwSrc(ng),                         &
     &            'Because LuvSrc is also ON; only one method is legal.'
          END IF
          exit_flag=4
          RETURN
        END IF
      END DO
!
!  Make sure that both component switches are activated when processing
!  (Eastward,Northward) momentum components at RHO-points.
!
      DO ng=1,Ngrids
        IF (.not.Hout(idu2dE,ng).and.Hout(idv2dN,ng)) THEN
          Hout(idu2dE,ng)=.TRUE.
        END IF
        IF (.not.Hout(idv2dN,ng).and.Hout(idu2dE,ng)) THEN
          Hout(idv2dN,ng)=.TRUE.
        END IF
        IF (.not.Hout(idu3dE,ng).and.Hout(idv3dN,ng)) THEN
          Hout(idu3dE,ng)=.TRUE.
        END IF
        IF (.not.Hout(idv3dN,ng).and.Hout(idu3dE,ng)) THEN
          Hout(idv3dN,ng)=.TRUE.
        END IF
      END DO
!
!  Set various parameters.
!
      DO ng=1,Ngrids
!
!  Set switch to create history NetCDF file.
!
        IF ((nHIS(ng).gt.0).and.ANY(Hout(:,ng))) THEN
          LdefHIS(ng)=.TRUE.
        END IF
!
!  Set switch to create quicksave NetCDF file.
!
        IF ((nQCK(ng).gt.0).and.ANY(Qout(:,ng))) THEN
          LdefQCK(ng)=.TRUE.
        END IF
!
!  Set switch to process climatology file.
!
        IF (LsshCLM(ng)) CLM_FILE(ng)=.TRUE.
        IF (Lm2CLM(ng)) CLM_FILE(ng)=.TRUE.
        IF (Lm3CLM(ng)) CLM_FILE(ng)=.TRUE.
        IF (ANY(LtracerCLM(:,ng))) CLM_FILE(ng)=.TRUE.
!
!  If appropriate, deactive outpur NetCDF files switches.
!
        IF (((nrrec(ng).eq.0).and.(nAVG(ng).gt.ntimes(ng))).or.         &
     &      (nAVG(ng).eq.0)) THEN
          LdefAVG(ng)=.FALSE.
        END IF
        IF (((nrrec(ng).eq.0).and.(nDIA(ng).gt.ntimes(ng))).or.         &
     &      (nDIA(ng).eq.0)) THEN
          LdefDIA(ng)=.FALSE.
        END IF
        IF (((nrrec(ng).eq.0).and.(nFLT(ng).gt.ntimes(ng))).or.         &
     &      (nFLT(ng).eq.0)) THEN
          LdefFLT(ng)=.FALSE.
        END IF
        IF (((nrrec(ng).eq.0).and.(nHIS(ng).gt.ntimes(ng))).or.         &
     &      (nHIS(ng).eq.0)) THEN
          LdefHIS(ng)=.FALSE.
        END IF
        IF (((nrrec(ng).eq.0).and.(nQCK(ng).gt.ntimes(ng))).or.         &
     &      (nQCK(ng).eq.0)) THEN
          LdefQCK(ng)=.FALSE.
        END IF
        IF (((nrrec(ng).eq.0).and.(nRST(ng).gt.ntimes(ng))).or.         &
     &      (nRST(ng).eq.0)) THEN
          LdefRST(ng)=.FALSE.
        END  IF
        IF (((nrrec(ng).eq.0).and.(nSTA(ng).gt.ntimes(ng))).or.         &
     &      (nSTA(ng).eq.0)) THEN
          LdefSTA(ng)=.FALSE.
        END IF
!
!  Determine switch to process boundary NetCDF file.
!
        ObcData(ng)=.FALSE.
        ObcData(ng)=ObcData(ng).or.ANY(LBC(:,isUvel,ng)%acquire)        &
     &                         .or.ANY(LBC(:,isVvel,ng)%acquire)
        ObcData(ng)=ObcData(ng).or.ANY(LBC(:,isTvar(:),ng)%acquire)
      END DO
!
!  If multiple output files, edit derived type structure to store the
!  information about all multi-files.
!
      DO ng=1,Ngrids
        IF ((nHIS(ng).gt.0).and.(ndefHIS(ng).gt.0)) THEN
          OutFiles=ntimes(ng)/ndefHIS(ng)
          IF ((nHIS(ng).eq.ndefHIS(ng)).or.                             &
     &        (MOD(ntimes(ng),ndefHIS(ng)).ge.nHIS(ng))) THEN
            OutFiles=Outfiles+1
          END IF
          CALL edit_file_struct (ng, OutFiles, HIS)
        END IF
        IF ((nQCK(ng).gt.0).and.(ndefQCK(ng).gt.0)) THEN
          OutFiles=ntimes(ng)/ndefQCK(ng)
          IF ((nQCK(ng).eq.ndefQCK(ng)).or.                             &
     &        (MOD(ntimes(ng),ndefQCK(ng)).ge.nQCK(ng))) THEN
            OutFiles=Outfiles+1
          END IF
          CALL edit_file_struct (ng, OutFiles, QCK)
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Master.and.Lwrite) THEN
        lstr=INDEX(my_fflags, 'free')-2
        IF (lstr.le.0) lstr=LEN_TRIM(my_fflags)
        WRITE (out,60) TRIM(title), TRIM(my_os), TRIM(my_cpu),          &
     &                 TRIM(my_fort), TRIM(my_fc), my_fflags(1:lstr),   &
     &                 OCN_COMM_WORLD, numthreads,                      &
     &                 TRIM(Iname),                                     &
     &                 TRIM(svn_url), TRIM(svn_rev),                    &
     &                 TRIM(Rdir), TRIM(Hdir), TRIM(Hfile), TRIM(Adir)
!
        DO ng=1,Ngrids
!
!  Report grid size and domain decomposition.  Check for correct tile
!  decomposition.
!
          WRITE (out,70) ng, Lm(ng), Mm(ng), N(ng), numthreads,         &
     &                   NtileI(ng), NtileJ(ng)
          maxPETs=numthreads
          nPETs=NtileI(ng)*NtileJ(ng)
          label='NtileI * NtileJ ='
          IF (nPETs.ne.maxPETs) THEN
            WRITE (out,80) ng, TRIM(label), nPETS, maxPETs
            exit_flag=6
            RETURN
          END IF
!
!  Report physical parameters.
!
          WRITE (out,110) ng
          WRITE (out,120) ntimes(ng), 'ntimes',                         &
     &          'Number of timesteps for 3-D equations.'
          WRITE (out,140) dt(ng), 'dt',                                 &
     &          'Timestep size (s) for 3-D equations.'
          WRITE (out,130) ndtfast(ng), 'ndtfast',                       &
     &          'Number of timesteps for 2-D equations between',        &
     &          'each 3D timestep.'
          WRITE (out,120) ERstr, 'ERstr',                               &
     &          'Starting ensemble/perturbation run number.'
          WRITE (out,120) ERend, 'ERend',                               &
     &          'Ending ensemble/perturbation run number.'
          WRITE (out,120) nrrec(ng), 'nrrec',                           &
     &          'Number of restart records to read from disk.'
          WRITE (out,170) LcycleRST(ng), 'LcycleRST',                   &
     &          'Switch to recycle time-records in restart file.'
          WRITE (out,130) nRST(ng), 'nRST',                             &
     &          'Number of timesteps between the writing of data',      &
     &          'into restart fields.'
          WRITE (out,130) ninfo(ng), 'ninfo',                           &
     &          'Number of timesteps between print of information',     &
     &          'to standard output.'
          WRITE (out,170) ldefout(ng), 'ldefout',                       &
     &          'Switch to create a new output NetCDF file(s).'
          WRITE (out,130) nHIS(ng), 'nHIS',                             &
     &          'Number of timesteps between the writing fields',       &
     &          'into history file.'
          IF (ndefHIS(ng).gt.0) THEN
            WRITE (out,130) ndefHIS(ng), 'ndefHIS',                     &
     &            'Number of timesteps between creation of new',        &
     &            'history files.'
          END IF
          WRITE (out,130) nQCK(ng), 'nQCK',                             &
     &          'Number of timesteps between the writing fields',       &
     &          'into quicksave file.'
          IF (ndefQCK(ng).gt.0) THEN
            WRITE (out,130) ndefQCK(ng), 'ndefQCK',                     &
     &            'Number of timesteps between creation of new',        &
     &            'brief snpashots files.'
          END IF
          DO i=1,NAT+NPT
            itrc=i
            WRITE (out,190) nl_tnu2(itrc,ng), 'nl_tnu2', itrc,          &
     &            'NLM Horizontal, harmonic mixing coefficient',        &
     &            '(m2/s) for tracer ', itrc,                           &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
          WRITE (out,210) nl_visc2(ng), 'nl_visc2',                     &
     &          'NLM Horizontal, harmonic mixing coefficient',          &
     &          '(m2/s) for momentum.'
          IF (LuvSponge(ng)) THEN
            WRITE (out,170) LuvSponge(ng), 'LuvSponge',                 &
     &          'Turning ON  sponge on horizontal momentum.'
          ELSE
            WRITE (out,170) LuvSponge(ng), 'LuvSponge',                 &
     &          'Turning OFF sponge on horizontal momentum.'
          END IF
          DO i=1,NAT
            IF (LtracerSponge(i,ng)) THEN
              WRITE (out,185) LtracerSponge(i,ng), 'LtracerSponge', i,  &
     &            'Turning ON  sponge on tracer ', i,                   &
     &            TRIM(Vname(1,idTvar(i)))
            ELSE
              WRITE (out,185) LtracerSponge(i,ng), 'LtracerSponge', i,  &
     &            'Turning OFF sponge on tracer ', i,                   &
     &            TRIM(Vname(1,idTvar(i)))
            END IF
          END DO
          DO i=1,NAT+NPT
            itrc=i
            WRITE (out,190) Akt_bak(itrc,ng), 'Akt_bak', itrc,          &
     &            'Background vertical mixing coefficient (m2/s)',      &
     &            'for tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT
            WRITE (out,190) Akt_limit(itrc,ng), 'Akt_limit', itrc,      &
     &            'Vertical diffusion upper threshold (m2/s)',          &
     &            'for tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
          WRITE (out,210) Akv_bak(ng), 'Akv_bak',                       &
     &          'Background vertical mixing coefficient (m2/s)',        &
     &          'for momentum.'
          WRITE (out,210) Akk_bak(ng), 'Akk_bak',                       &
     &          'Background vertical mixing coefficient (m2/s)',        &
     &          'for turbulent energy.'
          WRITE (out,210) Akp_bak(ng), 'Akp_bak',                       &
     &          'Background vertical mixing coefficient (m2/s)',        &
     &          'for turbulent generic statistical field.'
          WRITE (out,140) gls_p(ng), 'gls_p',                           &
     &          'GLS stability exponent.'
          WRITE (out,140) gls_m(ng), 'gls_m',                           &
     &          'GLS turbulent kinetic energy exponent.'
          WRITE (out,140) gls_n(ng), 'gls_n',                           &
     &          'GLS turbulent length scale exponent.'
          WRITE (out,200) gls_Kmin(ng), 'gls_Kmin',                     &
     &          'GLS minimum value of turbulent kinetic energy.'
          WRITE (out,200) gls_Pmin(ng), 'gls_Pmin',                     &
     &          'GLS minimum value of dissipation.'
          WRITE (out,200) gls_cmu0(ng), 'gls_cmu0',                     &
     &          'GLS stability coefficient.'
          WRITE (out,200) gls_c1(ng), 'gls_c1',                         &
     &          'GLS shear production coefficient.'
          WRITE (out,200) gls_c2(ng), 'gls_c2',                         &
     &          'GLS dissipation coefficient.'
          WRITE (out,200) gls_c3m(ng), 'gls_c3m',                       &
     &          'GLS stable buoyancy production coefficient.'
          WRITE (out,200) gls_c3p(ng), 'gls_c3p',                       &
     &          'GLS unstable buoyancy production coefficient.'
          WRITE (out,200) gls_sigk(ng), 'gls_sigk',                     &
     &          'GLS constant Schmidt number for TKE.'
          WRITE (out,200) gls_sigp(ng), 'gls_sigp',                     &
     &          'GLS constant Schmidt number for PSI.'
          WRITE (out,140) charnok_alpha(ng), 'charnok_alpha',           &
     &          'Charnok factor for Zos calculation.'
          WRITE (out,140) zos_hsig_alpha(ng), 'zos_hsig_alpha',         &
     &          'Factor for Zos calculation using Hsig(Awave).'
          WRITE (out,140) sz_alpha(ng), 'sz_alpha',                     &
     &          'Factor for Wave dissipation surface tke flux .'
          WRITE (out,140) crgban_cw(ng), 'crgban_cw',                   &
     &          'Factor for Craig/Banner surface tke flux.'
          WRITE (out,140) wec_alpha(ng), 'wec_alpha',                   &
     &          'WEC factor for roller/breaking energy distribution.'
          WRITE (out,200) rdrg(ng), 'rdrg',                             &
     &          'Linear bottom drag coefficient (m/s).'
          WRITE (out,200) rdrg2(ng), 'rdrg2',                           &
     &          'Quadratic bottom drag coefficient.'
          WRITE (out,200) Zob(ng), 'Zob',                               &
     &          'Bottom roughness (m).'
          IF (Zob(ng).le.0.0_r8) THEN
            WRITE (out,265) 'Zob = ', Zob(ng),                          &
     &            'It must be greater than zero when BBL is activated.'
            exit_flag=5
            RETURN
          END IF
          WRITE (out,200) Zos(ng), 'Zos',                               &
     &          'Surface roughness (m).'
          WRITE (out,200) Dcrit(ng), 'Dcrit',                           &
     &          'Minimum depth for wetting and drying (m).'
          WRITE (out,120) Vtransform(ng), 'Vtransform',                 &
     &          'S-coordinate transformation equation.'
          WRITE (out,120) Vstretching(ng), 'Vstretching',               &
     &          'S-coordinate stretching function.'
          WRITE (out,200) theta_s(ng), 'theta_s',                       &
     &          'S-coordinate surface control parameter.'
          WRITE (out,200) theta_b(ng), 'theta_b',                       &
     &          'S-coordinate bottom  control parameter.'
          IF (Tcline(ng).gt.1.0E+5_r8) THEN
            WRITE (out,210) Tcline(ng), 'Tcline',                       &
     &            'S-coordinate surface/bottom layer width (m) used',   &
     &            'in vertical coordinate stretching.'
          ELSE
            WRITE (out,160) Tcline(ng), 'Tcline',                       &
     &            'S-coordinate surface/bottom layer width (m) used',   &
     &            'in vertical coordinate stretching.'
          END IF
          WRITE (out,140) rho0, 'rho0',                                 &
     &          'Mean density (kg/m3) for Boussinesq approximation.'
          WRITE (out,140) dstart, 'dstart',                             &
     &          'Time-stamp assigned to model initialization (days).'
          WRITE (out,150) time_ref, 'time_ref',                         &
     &          'Reference time for units attribute (yyyymmdd.dd)'
          DO i=1,NAT+NPT
            itrc=i
            WRITE (out,190) Tnudg(itrc,ng), 'Tnudg', itrc,              &
     &            'Nudging/relaxation time scale (days)',               &
     &            'for tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
          WRITE (out,210) Znudg(ng), 'Znudg',                           &
     &          'Nudging/relaxation time scale (days)',                 &
     &          'for free-surface.'
          WRITE (out,210) M2nudg(ng), 'M2nudg',                         &
     &          'Nudging/relaxation time scale (days)',                 &
     &          'for 2D momentum.'
          WRITE (out,210) M3nudg(ng), 'M3nudg',                         &
     &          'Nudging/relaxation time scale (days)',                 &
     &          'for 3D momentum.'
          WRITE (out,210) obcfac(ng), 'obcfac',                         &
     &          'Factor between passive and active',                    &
     &          'open boundary conditions.'
          WRITE (out,170) VolCons(1,ng), 'VolCons(1)',                  &
     &          'NLM western  edge boundary volume conservation.'
          WRITE (out,170) VolCons(2,ng), 'VolCons(2)',                  &
     &          'NLM southern edge boundary volume conservation.'
          WRITE (out,170) VolCons(3,ng), 'VolCons(3)',                  &
     &          'NLM eastern  edge boundary volume conservation.'
          WRITE (out,170) VolCons(4,ng), 'VolCons(4)',                  &
     &          'NLM northern edge boundary volume conservation.'
          WRITE (out,140) T0(ng), 'T0',                                 &
     &          'Background potential temperature (C) constant.'
          WRITE (out,140) S0(ng), 'S0',                                 &
     &          'Background salinity (PSU) constant.'
          WRITE (out,160) R0(ng), 'R0',                                 &
     &          'Background density (kg/m3) used in linear Equation',   &
     &          'of State.'
          WRITE (out,200) Tcoef(ng), 'Tcoef',                           &
     &          'Thermal expansion coefficient (1/Celsius).'
          WRITE (out,200) Scoef(ng), 'Scoef',                           &
     &          'Saline contraction coefficient (1/PSU).'
          WRITE (out,160) gamma2(ng), 'gamma2',                         &
     &          'Slipperiness variable: free-slip (1.0) or ',           &
     &          '                     no-slip (-1.0).'
          IF (LuvSrc(ng)) THEN
            WRITE (out,170) LuvSrc(ng), 'LuvSrc',                       &
     &          'Turning ON  momentum point Sources/Sinks.'
          ELSE
            WRITE (out,170) LuvSrc(ng), 'LuvSrc',                       &
     &          'Turning OFF momentum point Sources/Sinks.'
          END IF
          IF (LwSrc(ng)) THEN
            WRITE (out,170) LwSrc(ng), 'LwSrc',                         &
     &          'Turning ON  volume influx point Sources/Sinks.'
          ELSE
            WRITE (out,170) LwSrc(ng), 'LwSrc',                         &
     &          'Turning OFF volume influx point Sources/Sinks.'
          END IF
          DO itrc=1,NAT
            IF (LtracerSrc(itrc,ng)) THEN
              WRITE (out,185) LtracerSrc(itrc,ng), 'LtracerSrc', itrc,  &
     &            'Turning ON  point Sources/Sinks on tracer ', itrc,   &
     &            TRIM(Vname(1,idTvar(itrc)))
            ELSE
              WRITE (out,185) LtracerSrc(itrc,ng), 'LtracerSrc', itrc,  &
     &            'Turning OFF point Sources/Sinks on tracer ', itrc,   &
     &            TRIM(Vname(1,idTvar(itrc)))
            END IF
          END DO
          IF (LsshCLM(ng)) THEN
            WRITE (out,170) LsshCLM(ng), 'LsshCLM',                     &
     &          'Turning ON  processing of SSH climatology.'
          ELSE
            WRITE (out,170) LsshCLM(ng), 'LsshCLM',                     &
     &          'Turning OFF processing of SSH climatology.'
          END IF
          IF (Lm2CLM(ng)) THEN
            WRITE (out,170) Lm2CLM(ng), 'Lm2CLM',                       &
     &          'Turning ON  processing of 2D momentum climatology.'
          ELSE
            WRITE (out,170) Lm2CLM(ng), 'Lm2CLM',                       &
     &          'Turning OFF processing of 2D momentum climatology.'
          END IF
          IF (Lm3CLM(ng)) THEN
            WRITE (out,170) Lm3CLM(ng), 'Lm3CLM',                       &
     &          'Turning ON  processing of 3D momentum climatology.'
          ELSE
            WRITE (out,170) Lm3CLM(ng), 'Lm3CLM',                       &
     &          'Turning OFF processing of 3D momentum climatology.'
          END IF
          DO i=1,NAT
            IF (LtracerCLM(i,ng)) THEN
              WRITE (out,185) LtracerCLM(i,ng), 'LtracerCLM', i,        &
     &            'Turning ON  processing of climatology tracer ', i,   &
     &            TRIM(Vname(1,idTvar(i)))
            ELSE
              WRITE (out,185) LtracerCLM(i,ng), 'LtracerCLM', i,        &
     &            'Turning OFF processing of climatology tracer ', i,   &
     &            TRIM(Vname(1,idTvar(i)))
            END IF
          END DO
          IF (LnudgeM2CLM(ng)) THEN
            WRITE (out,170) LnudgeM2CLM(ng), 'LnudgeM2CLM',             &
     &          'Turning ON  nudging of 2D momentum climatology.'
          ELSE
            WRITE (out,170) LnudgeM2CLM(ng), 'LnudgeM2CLM',             &
     &          'Turning OFF nudging of 2D momentum climatology.'
          END IF
          IF (LnudgeM3CLM(ng)) THEN
            WRITE (out,170) LnudgeM3CLM(ng), 'LnudgeM3CLM',             &
     &          'Turning ON  nudging of 3D momentum climatology.'
          ELSE
            WRITE (out,170) LnudgeM3CLM(ng), 'LnudgeM3CLM',             &
     &          'Turning OFF nudging of 3D momentum climatology.'
          END IF
          DO i=1,NAT
            IF (LnudgeTCLM(i,ng)) THEN
              WRITE (out,185) LnudgeTCLM(i,ng), 'LnudgeTCLM', i,        &
     &            'Turning ON  nudging of climatology tracer ', i,      &
     &            TRIM(Vname(1,idTvar(i)))
            ELSE
              WRITE (out,185) LnudgeTCLM(i,ng), 'LnudgeTCLM', i,        &
     &            'Turning OFF nudging of climatology tracer ', i,      &
     &            TRIM(Vname(1,idTvar(i)))
            END IF
          END DO
          IF ((nHIS(ng).gt.0).and.ANY(Hout(:,ng))) THEN
            WRITE (out,'(1x)')
            IF (Hout(idFsur,ng)) WRITE (out,170) Hout(idFsur,ng),       &
     &         'Hout(idFsur)',                                          &
     &         'Write out free-surface.'
            IF (Hout(idUbar,ng)) WRITE (out,170) Hout(idUbar,ng),       &
     &         'Hout(idUbar)',                                          &
     &         'Write out 2D U-momentum component.'
            IF (Hout(idVbar,ng)) WRITE (out,170) Hout(idVbar,ng),       &
     &         'Hout(idVbar)',                                          &
     &         'Write out 2D V-momentum component.'
            IF (Hout(idu2dE,ng)) WRITE (out,170) Hout(idu2dE,ng),       &
     &         'Hout(idu2dE)',                                          &
     &         'Write out 2D U-eastward  at RHO-points.'
            IF (Hout(idv2dN,ng)) WRITE (out,170) Hout(idv2dN,ng),       &
     &         'Hout(idv2dN)',                                          &
     &         'Write out 2D V-northward at RHO-points.'
            IF (Hout(idUvel,ng)) WRITE (out,170) Hout(idUvel,ng),       &
     &         'Hout(idUvel)',                                          &
     &         'Write out 3D U-momentum component.'
            IF (Hout(idVvel,ng)) WRITE (out,170) Hout(idVvel,ng),       &
     &         'Hout(idVvel)',                                          &
     &         'Write out 3D V-momentum component.'
            IF (Hout(idu3dE,ng)) WRITE (out,170) Hout(idu3dE,ng),       &
     &         'Hout(idu3dE)',                                          &
     &         'Write out 3D U-wastward  component at RHO-points.'
            IF (Hout(idv3dN,ng)) WRITE (out,170) Hout(idv3dN,ng),       &
     &         'Hout(idv3dN)',                                          &
     &         'Write out 3D V-northward component at RHO-points.'
            IF (Hout(idWvel,ng)) WRITE (out,170) Hout(idWvel,ng),       &
     &         'Hout(idWvel)',                                          &
     &         'Write out W-momentum component.'
            IF (Hout(idOvel,ng)) WRITE (out,170) Hout(idOvel,ng),       &
     &         'Hout(idOvel)',                                          &
     &         'Write out omega vertical velocity.'
            DO itrc=1,NAT
              IF (Hout(idTvar(itrc),ng)) WRITE (out,180)                &
     &            Hout(idTvar(itrc),ng), 'Hout(idTvar)',                &
     &            'Write out tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
            END DO
            IF (Hout(idpthR,ng)) WRITE (out,170) Hout(idpthR,ng),       &
     &         'Hout(idpthR)',                                          &
     &         'Write out time-varying dephts of RHO-points.'
            IF (Hout(idpthU,ng)) WRITE (out,170) Hout(idpthU,ng),       &
     &         'Hout(idpthU)',                                          &
     &         'Write out time-varying dephts of U-points.'
            IF (Hout(idpthV,ng)) WRITE (out,170) Hout(idpthV,ng),       &
     &         'Hout(idpthV)',                                          &
     &         'Write out time-varying dephts of V-points.'
            IF (Hout(idpthW,ng)) WRITE (out,170) Hout(idpthW,ng),       &
     &         'Hout(idpthW)',                                          &
     &         'Write out time-varying dephts of W-points.'
            IF (Hout(idUsms,ng)) WRITE (out,170) Hout(idUsms,ng),       &
     &         'Hout(idUsms)',                                          &
     &         'Write out surface U-momentum stress.'
            IF (Hout(idVsms,ng)) WRITE (out,170) Hout(idVsms,ng),       &
     &         'Hout(idVsms)',                                          &
     &         'Write out surface V-momentum stress.'
            IF (Hout(idUbms,ng)) WRITE (out,170) Hout(idUbms,ng),       &
     &         'Hout(idUbms)',                                          &
     &         'Write out bottom U-momentum stress.'
            IF (Hout(idVbms,ng)) WRITE (out,170) Hout(idVbms,ng),       &
     &         'Hout(idVbms)',                                          &
     &         'Write out bottom V-momentum stress.'
            IF (Hout(idUbrs,ng)) WRITE (out,170) Hout(idUbrs,ng),       &
     &         'Hout(idUbrs)',                                          &
     &         'Write out bottom U-current stress.'
            IF (Hout(idVbrs,ng)) WRITE (out,170) Hout(idVbrs,ng),       &
     &         'Hout(idVbrs)',                                          &
     &         'Write out bottom V-current stress.'
            IF (Hout(idUbws,ng)) WRITE (out,170) Hout(idUbws,ng),       &
     &         'Hout(idUbws)',                                          &
     &         'Write out wind-induced, bottom U-wave stress.'
            IF (Hout(idVbws,ng)) WRITE (out,170) Hout(idVbws,ng),       &
     &         'Hout(idVbws)',                                          &
     &         'Write out wind-induced, bottom V-wave stress.'
            IF (Hout(idUbcs,ng)) WRITE (out,170) Hout(idUbcs,ng),       &
     &         'Hout(idUbcs)',                                          &
     &         'Write out max wind + current, bottom U-wave stress.'
            IF (Hout(idVbcs,ng)) WRITE (out,170) Hout(idVbcs,ng),       &
     &         'Hout(idVbcs)',                                          &
     &         'Write out max wind + current, bottom V-wave stress.'
            IF (Hout(idUVwc,ng)) WRITE (out,170) Hout(idUVwc,ng),       &
     &         'Hout(idUVwc)',                                          &
     &         'Write out max wind + current, bottom UV-wave stress.'
            IF (Hout(idUbot,ng)) WRITE (out,170) Hout(idUbot,ng),       &
     &         'Hout(idUbot)',                                          &
     &         'Write out bed wave orbital U-velocity.'
            IF (Hout(idVbot,ng)) WRITE (out,170) Hout(idVbot,ng),       &
     &         'Hout(idVbot)',                                          &
     &         'Write out bed wave orbital V-velocity.'
            IF (Hout(idUbur,ng)) WRITE (out,170) Hout(idUbur,ng),       &
     &         'Hout(idUbur)',                                          &
     &         'Write out bottom U-momentum above bed.'
            IF (Hout(idVbvr,ng)) WRITE (out,170) Hout(idVbvr,ng),       &
     &         'Hout(idVbvr)',                                          &
     &         'Write out bottom V-momentum above bed.'
            IF (Hout(idU2rs,ng)) WRITE (out,170) Hout(idU2rs,ng),       &
     &         'Hout(idU2rs)',                                          &
     &         'Write out total 2D u-radiation stress.'
            IF (Hout(idV2rs,ng)) WRITE (out,170) Hout(idV2rs,ng),       &
     &         'Hout(idV2rs)',                                          &
     &         'Write out total 2D v-radiation stress.'
            IF (Hout(idU2Sd,ng)) WRITE (out,170) Hout(idU2Sd,ng),       &
     &         'Hout(idU2Sd)',                                          &
     &         'Write out 2D u-momentum stokes velocity.'
            IF (Hout(idV2Sd,ng)) WRITE (out,170) Hout(idV2Sd,ng),       &
     &         'Hout(idV2Sd)',                                          &
     &         'Write out 2D v-momentum stokes velocity.'
            IF (Hout(idU3rs,ng)) WRITE (out,170) Hout(idU3rs,ng),       &
     &         'Hout(idU3rs)',                                          &
     &         'Write out total 3D u-radiation stress.'
            IF (Hout(idV3rs,ng)) WRITE (out,170) Hout(idV3rs,ng),       &
     &         'Hout(idV3rs)',                                          &
     &         'Write out total 3D v-radiation stress.'
            IF (Hout(idU3Sd,ng)) WRITE (out,170) Hout(idU3Sd,ng),       &
     &         'Hout(idU3Sd)',                                          &
     &         'Write out 3D u-momentum stokes velocity.'
            IF (Hout(idV3Sd,ng)) WRITE (out,170) Hout(idV3Sd,ng),       &
     &         'Hout(idV3Sd)',                                          &
     &         'Write out 3D v-momentum stokes velocity.'
            IF (Hout(idW3Sd,ng)) WRITE (out,170) Hout(idW3Sd,ng),       &
     &         'Hout(idW3Sd)',                                          &
     &         'Write out 3D omega-momentum stokes velocity.'
            IF (Hout(idW3St,ng)) WRITE (out,170) Hout(idW3St,ng),       &
     &         'Hout(idW3St)',                                          &
     &         'Write out 3D w-momentum stokes velocity.'
            IF (Hout(idWztw,ng)) WRITE (out,170) Hout(idWztw,ng),       &
     &         'Hout(idWztw)',                                          &
     &         'Write out WEC quasi-static sea level adjustment.'
            IF (Hout(idWqsp,ng)) WRITE (out,170) Hout(idWqsp,ng),       &
     &         'Hout(idWqsp)',                                          &
     &         'Write out WEC quasi-static pressure.'
            IF (Hout(idWbeh,ng)) WRITE (out,170) Hout(idWbeh,ng),       &
     &         'Hout(idWbeh)',                                          &
     &         'Write out WEC Bernoulli head.'
            IF (Hout(idWamp,ng)) WRITE (out,170) Hout(idWamp,ng),       &
     &         'Hout(idWamp)',                                          &
     &         'Write out wave height.'
            IF (Hout(idWlen,ng)) WRITE (out,170) Hout(idWlen,ng),       &
     &         'Hout(idWlen)',                                          &
     &         'Write out wavelength.'
            IF (Hout(idWdir,ng)) WRITE (out,170) Hout(idWdir,ng),       &
     &         'Hout(idWdir)',                                          &
     &         'Write out mean wave direction.'
            IF (Hout(idWdip,ng)) WRITE (out,170) Hout(idWdip,ng),       &
     &         'Hout(idWdip)',                                          &
     &         'Write out peak wave direction.'
            IF (Hout(idWptp,ng)) WRITE (out,170) Hout(idWptp,ng),       &
     &         'Hout(idWptp)',                                          &
     &         'Write out wave surface period.'
            IF (Hout(idWpbt,ng)) WRITE (out,170) Hout(idWpbt,ng),       &
     &         'Hout(idWpbt)',                                          &
     &         'Write out wave bottom period.'
            IF (Hout(idWorb,ng)) WRITE (out,170) Hout(idWorb,ng),       &
     &         'Hout(idWorb)',                                          &
     &         'Write out wave bottom orbital velocity.'
            IF (Hout(idWdif,ng)) WRITE (out,170) Hout(idWdif,ng),       &
     &         'Hout(idWdif)',                                          &
     &         'Write out wave dissipation due to bottom friction.'
            IF (Hout(idWdib,ng)) WRITE (out,170) Hout(idWdib,ng),       &
     &         'Hout(idWdib)',                                          &
     &         'Write out wave dissipation due to breaking.'
            IF (Hout(idWdiw,ng)) WRITE (out,170) Hout(idWdiw,ng),       &
     &         'Hout(idWdiw)',                                          &
     &         'Write out wave dissipation due to white capping.'
            IF (Hout(idWdis,ng)) WRITE (out,170) Hout(idWdis,ng),       &
     &         'Hout(idWdis)',                                          &
     &         'Write out wave roller dissipation.'
            IF (Hout(idWrol,ng)) WRITE (out,170) Hout(idWrol,ng),       &
     &         'Hout(idWrol)',                                          &
     &         'Write out wave roller action density.'
            IF (Hout(idUwav,ng)) WRITE (out,170) Hout(idUwav,ng),       &
     &         'Hout(idUwav)',                                          &
     &         'Wave-avg surface u-velocity.'
            IF (Hout(idVwav,ng)) WRITE (out,170) Hout(idVwav,ng),       &
     &         'Hout(idVwav)',                                          &
     &         'Wave-avg surface v-velocity.'
            IF (Hout(idTsur(itemp),ng)) WRITE (out,170)                 &
     &          Hout(idTsur(itemp),ng), 'Hout(idTsur)',                 &
     &         'Write out surface net heat flux.'
            IF (Hout(idDano,ng)) WRITE (out,170) Hout(idDano,ng),       &
     &         'Hout(idDano)',                                          &
     &         'Write out density anomaly.'
            IF (Hout(idVvis,ng)) WRITE (out,170) Hout(idVvis,ng),       &
     &         'Hout(idVvis)',                                          &
     &         'Write out vertical viscosity: AKv.'
            IF (Hout(idTdif,ng)) WRITE (out,170) Hout(idTdif,ng),       &
     &         'Hout(idTdif)',                                          &
     &         'Write out vertical diffusion: AKt(itemp).'
            IF (Hout(idMtke,ng)) WRITE (out,170) Hout(idMtke,ng),       &
     &         'Hout(idMtke)',                                          &
     &         'Write out turbulent kinetic energy.'
            IF (Hout(idMtls,ng)) WRITE (out,170) Hout(idMtls,ng),       &
     &         'Hout(idMtls)',                                          &
     &         'Write out turbulent generic length-scale.'
          END IF
          IF ((nQCK(ng).gt.0).and.ANY(Qout(:,ng))) THEN
            WRITE (out,'(1x)')
            IF (Qout(idFsur,ng)) WRITE (out,170) Qout(idFsur,ng),       &
     &         'Qout(idFsur)',                                          &
     &         'Write out free-surface.'
            IF (Qout(idUbar,ng)) WRITE (out,170) Qout(idUbar,ng),       &
     &         'Qout(idUbar)',                                          &
     &         'Write out 2D U-momentum component.'
            IF (Qout(idVbar,ng)) WRITE (out,170) Qout(idVbar,ng),       &
     &         'Qout(idVbar)',                                          &
     &         'Write out 2D V-momentum component.'
            IF (Qout(idu2dE,ng)) WRITE (out,170) Qout(idu2dE,ng),       &
     &         'Qout(idu2dE)',                                          &
     &         'Write out 2D U-eastward  at RHO-points.'
            IF (Qout(idv2dN,ng)) WRITE (out,170) Qout(idv2dN,ng),       &
     &         'Qout(idv2dN)',                                          &
     &         'Write out 2D V-northward at RHO-points.'
            IF (Qout(idUvel,ng)) WRITE (out,170) Qout(idUvel,ng),       &
     &         'Qout(idUvel)',                                          &
     &         'Write out 3D U-momentum component.'
            IF (Qout(idVvel,ng)) WRITE (out,170) Qout(idVvel,ng),       &
     &         'Qout(idVvel)',                                          &
     &         'Write out 3D V-momentum component.'
            IF (Qout(idUsur,ng)) WRITE (out,170) Qout(idUsur,ng),       &
     &         'Qout(idUsur)',                                          &
     &         'Write out surface U-momentum component.'
            IF (Qout(idVsur,ng)) WRITE (out,170) Qout(idVsur,ng),       &
     &         'Qout(idVsur)',                                          &
     &         'Write out surface V-momentum component.'
            IF (Qout(idu3dE,ng)) WRITE (out,170) Qout(idu3dE,ng),       &
     &         'Qout(idu3dE)',                                          &
     &         'Write out 3D U-wastward  component at RHO-points.'
            IF (Qout(idv3dN,ng)) WRITE (out,170) Qout(idv3dN,ng),       &
     &         'Qout(idv3dN)',                                          &
     &         'Write out 3D V-northward component at RHO-points.'
            IF (Qout(idu3dE,ng)) WRITE (out,170) Qout(idu3dE,ng),       &
     &         'Qout(idu3dE)',                                          &
     &         'Write out surface U-wastward  component at RHO-points.'
            IF (Qout(idv3dN,ng)) WRITE (out,170) Qout(idv3dN,ng),       &
     &         'Qout(idv3dN)',                                          &
     &         'Write out surface V-northward component at RHO-points.'
            IF (Qout(idWvel,ng)) WRITE (out,170) Qout(idWvel,ng),       &
     &         'Qout(idWvel)',                                          &
     &         'Write out W-momentum component.'
            IF (Qout(idOvel,ng)) WRITE (out,170) Qout(idOvel,ng),       &
     &         'Qout(idOvel)',                                          &
     &         'Write out omega vertical velocity.'
            DO itrc=1,NAT
              IF (Qout(idTvar(itrc),ng)) WRITE (out,180)                &
     &            Qout(idTvar(itrc),ng), 'Qout(idTvar)',                &
     &            'Write out tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO itrc=1,NAT
              IF (Qout(idsurT(itrc),ng)) WRITE (out,180)                &
     &            Qout(idsurT(itrc),ng), 'Qout(idsurT)',                &
     &            'Write out surface tracer ', itrc,                    &
     &            TRIM(Vname(1,idsurT(itrc)))
            END DO
            IF (Qout(idpthR,ng)) WRITE (out,170) Qout(idpthR,ng),       &
     &         'Qout(idpthR)',                                          &
     &         'Write out time-varying dephts of RHO-points.'
            IF (Qout(idpthU,ng)) WRITE (out,170) Qout(idpthU,ng),       &
     &         'Qout(idpthU)',                                          &
     &         'Write out time-varying dephts of U-points.'
            IF (Qout(idpthV,ng)) WRITE (out,170) Qout(idpthV,ng),       &
     &         'Qout(idpthV)',                                          &
     &         'Write out time-varying dephts of V-points.'
            IF (Qout(idpthW,ng)) WRITE (out,170) Qout(idpthW,ng),       &
     &         'Qout(idpthW)',                                          &
     &         'Write out time-varying dephts of W-points.'
            IF (Qout(idUsms,ng)) WRITE (out,170) Qout(idUsms,ng),       &
     &         'Qout(idUsms)',                                          &
     &         'Write out surface U-momentum stress.'
            IF (Qout(idVsms,ng)) WRITE (out,170) Qout(idVsms,ng),       &
     &         'Qout(idVsms)',                                          &
     &         'Write out surface V-momentum stress.'
            IF (Qout(idUbms,ng)) WRITE (out,170) Qout(idUbms,ng),       &
     &         'Qout(idUbms)',                                          &
     &         'Write out bottom U-momentum stress.'
            IF (Qout(idVbms,ng)) WRITE (out,170) Qout(idVbms,ng),       &
     &         'Qout(idVbms)',                                          &
     &         'Write out bottom V-momentum stress.'
            IF (Qout(idUbrs,ng)) WRITE (out,170) Qout(idUbrs,ng),       &
     &         'Qout(idUbrs)',                                          &
     &         'Write out bottom U-current stress.'
            IF (Qout(idVbrs,ng)) WRITE (out,170) Qout(idVbrs,ng),       &
     &         'Qout(idVbrs)',                                          &
     &         'Write out bottom V-current stress.'
            IF (Qout(idUbws,ng)) WRITE (out,170) Qout(idUbws,ng),       &
     &         'Qout(idUbws)',                                          &
     &         'Write out wind-induced, bottom U-wave stress.'
            IF (Qout(idVbws,ng)) WRITE (out,170) Qout(idVbws,ng),       &
     &         'Qout(idVbws)',                                          &
     &         'Write out wind-induced, bottom V-wave stress.'
            IF (Qout(idUbcs,ng)) WRITE (out,170) Qout(idUbcs,ng),       &
     &         'Qout(idUbcs)',                                          &
     &         'Write out max wind + current, bottom U-wave stress.'
            IF (Qout(idVbcs,ng)) WRITE (out,170) Qout(idVbcs,ng),       &
     &         'Qout(idVbcs)',                                          &
     &         'Write out max wind + current, bottom V-wave stress.'
            IF (Qout(idUbot,ng)) WRITE (out,170) Qout(idUbot,ng),       &
     &         'Qout(idUbot)',                                          &
     &         'Write out bed wave orbital U-velocity.'
            IF (Qout(idVbot,ng)) WRITE (out,170) Qout(idVbot,ng),       &
     &         'Qout(idVbot)',                                          &
     &         'Write out bed wave orbital V-velocity.'
            IF (Qout(idUbur,ng)) WRITE (out,170) Qout(idUbur,ng),       &
     &         'Qout(idUbur)',                                          &
     &         'Write out bottom U-momentum above bed.'
            IF (Qout(idVbvr,ng)) WRITE (out,170) Qout(idVbvr,ng),       &
     &         'Qout(idVbvr)',                                          &
     &         'Write out bottom V-momentum above bed.'
            IF (Qout(idWamp,ng)) WRITE (out,170) Qout(idWamp,ng),       &
     &         'Qout(idWamp)',                                          &
     &         'Write out wave height.'
            IF (Qout(idWlen,ng)) WRITE (out,170) Qout(idWlen,ng),       &
     &         'Qout(idWlen)',                                          &
     &         'Write out wavelength.'
            IF (Qout(idWdir,ng)) WRITE (out,170) Qout(idWdir,ng),       &
     &         'Qout(idWdir)',                                          &
     &         'Write out mean wave direction.'
            IF (Qout(idWdip,ng)) WRITE (out,170) Qout(idWdip,ng),       &
     &         'Qout(idWdip)',                                          &
     &         'Write out peak wave direction.'
            IF (Qout(idWptp,ng)) WRITE (out,170) Qout(idWptp,ng),       &
     &         'Qout(idWptp)',                                          &
     &         'Write out wave surface period.'
            IF (Qout(idWpbt,ng)) WRITE (out,170) Qout(idWpbt,ng),       &
     &         'Qout(idWpbt)',                                          &
     &         'Write out wave bottom period.'
            IF (Qout(idWorb,ng)) WRITE (out,170) Qout(idWorb,ng),       &
     &         'Qout(idWorb)',                                          &
     &         'Write out wave bottom orbital velocity.'
            IF (Qout(idWdis,ng)) WRITE (out,170) Qout(idWdis,ng),       &
     &         'Qout(idWdis)',                                          &
     &         'Write out wave dissipation.'
            IF (Qout(idTsur(itemp),ng)) WRITE (out,170)                 &
     &          Qout(idTsur(itemp),ng), 'Qout(idTsur)',                 &
     &         'Write out surface net heat flux.'
            IF (Qout(idDano,ng)) WRITE (out,170) Qout(idDano,ng),       &
     &         'Qout(idDano)',                                          &
     &         'Write out density anomaly.'
            IF (Qout(idVvis,ng)) WRITE (out,170) Qout(idVvis,ng),       &
     &         'Qout(idVvis)',                                          &
     &         'Write out vertical viscosity: AKv.'
            IF (Qout(idTdif,ng)) WRITE (out,170) Qout(idTdif,ng),       &
     &         'Qout(idTdif)',                                          &
     &         'Write out vertical diffusion: AKt(itemp).'
            IF (Qout(idMtke,ng)) WRITE (out,170) Qout(idMtke,ng),       &
     &         'Qout(idMtke)',                                          &
     &         'Write out turbulent kinetic energy.'
            IF (Qout(idMtls,ng)) WRITE (out,170) Qout(idMtls,ng),       &
     &         'Qout(idMtls)',                                          &
     &         'Write out turbulent generic length-scale.'
          END IF
!
!-----------------------------------------------------------------------
!  Report output/input files and check availability of input files.
!-----------------------------------------------------------------------
!
          WRITE (out,220)
          WRITE (out,230) '             Output Restart File:  ',        &
     &                    TRIM(RST(ng)%name)
          IF (LdefHIS(ng)) THEN
            IF (ndefHIS(ng).eq.0) THEN
              WRITE (out,230) '             Output History File:  ',    &
     &                        TRIM(HIS(ng)%name)
            ELSE
              WRITE (out,230) '        Prefix for History Files:  ',    &
     &                        TRIM(HIS(ng)%head)
            END IF
          END IF
          WRITE (out,230) '        Physical parameters File:  ',        &
     &                    TRIM(Iname)
          fname=GRD(ng)%name
          IF (.not.find_file(ng, fname, 'GRDNAME')) GO TO 30
          WRITE (out,230) '                 Input Grid File:  ',        &
     &                    TRIM(fname)
          fname=INI(ng)%name
          IF (.not.find_file(ng, fname, 'ININAME')) GO TO 30
          WRITE (out,230) '    Input Nonlinear Initial File:  ',        &
     &                    TRIM(fname)
          IF (LuvSrc(ng).or.LwSrc(ng).or.(ANY(LtracerSrc(:,ng)))) THEN
            fname=SSF(ng)%name
            IF (.not.find_file(ng, fname, 'SSFNAME')) GO TO 30
            WRITE (out,230) '        Input Sources/Sinks File:  ',      &
     &                      TRIM(fname)
          END IF
          DO i=1,nCLMfiles(ng)
            IF (CLM_FILE(ng)) THEN
              DO ifile=1,CLM(i,ng)%Nfiles
                fname=CLM(i,ng)%files(ifile)
                IF (.not.find_file(ng, fname, 'CLMNAME')) GO TO 30
                IF (ifile.eq.1) THEN
                  WRITE (out,310) '       Input Climatology File ', i,  &
     &                            ':  ', TRIM(fname)
                ELSE
                  WRITE (out,'(37x,a)') TRIM(fname)
                END IF
              END DO
            END IF
          END DO
          IF (LnudgeM2CLM(ng).or.LnudgeM3CLM(ng).or.                    &
     &        (ANY(LnudgeTCLM(:,ng)))) THEN
            fname=NUD(ng)%name
            IF (.not.find_file(ng, fname, 'NUDNAME')) GO TO 30
            WRITE (out,230) '   Input Nudge Coefficients File:  ',      &
     &                      TRIM(fname)
          END IF
!
          IF (ObcData(ng)) THEN
            DO i=1,nBCfiles(ng)
              DO ifile=1,BRY(i,ng)%Nfiles
                fname=BRY(i,ng)%files(ifile)
                IF (.not.find_file(ng, fname, 'BRYNAME')) GO TO 30
                IF (ifile.eq.1) THEN
                  WRITE (out,310) '  Input Lateral Boundary File ', i,  &
     &                            ':  ', TRIM(fname)
                ELSE
                  WRITE (out,'(37x,a)') TRIM(fname)
                END IF
              END DO
            END DO
          END IF
          fname=varname
          IF (.not.find_file(ng, fname, 'VARNAME')) GO TO 30
          WRITE (out,230) 'ROMS I/O variables Metadata File:  ',        &
     &                    TRIM(fname)
          GO TO 40
  30      IF (Master) THEN
            IF (LEN_TRIM(fname).eq.0) THEN
              WRITE (out,270) ng, 'Oops unassigned file name. '//       &
     &                            'Check standard input script...'
            ELSE
              WRITE (out,270) ng, TRIM(fname)
            END IF
          END IF
          exit_flag=4
          IF (FoundError(exit_flag, NoError, 7346, MyFile)) RETURN
  40      CONTINUE
        END DO
        IF (Nuser.gt.0) THEN
          WRITE (out,230) '          Input/Output USER File:  ',        &
     &                    TRIM(USRname)
        END IF
!
!-----------------------------------------------------------------------
!  Report generic USER parameters.
!-----------------------------------------------------------------------
!
        IF (Nuser.gt.0) THEN
          WRITE (out,240)
          DO i=1,Nuser
            WRITE (out,250) user(i), i, i
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Rescale active tracer parameters
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO i=1,NAT+NPT
          itrc=i
!
!  Take the square root of the biharmonic coefficients so it can
!  be applied to each harmonic operator.
!
          nl_tnu4(itrc,ng)=SQRT(ABS(nl_tnu4(itrc,ng)))
!
!  Compute inverse nudging coefficients (1/s) used in various tasks.
!
          IF (Tnudg(itrc,ng).gt.0.0_r8) THEN
            Tnudg(itrc,ng)=1.0_r8/(Tnudg(itrc,ng)*86400.0_r8)
          ELSE
            Tnudg(itrc,ng)=0.0_r8
          END IF
        END DO
        IF (Tnudg_SSS(ng).gt.0.0_r8) THEN
          Tnudg_SSS(ng)=1.0_r8/(Tnudg_SSS(ng)*86400.0_r8)
        END IF
      END DO
  50  FORMAT (/,' READ_PhyPar - Error while processing line: ',/,a)
  60  FORMAT (/,1x,a,/,                                                 &
     &        /,1x,'Operating system  : ',a,                            &
     &        /,1x,'CPU/hardware      : ',a,                            &
     &        /,1x,'Compiler system   : ',a,                            &
     &        /,1x,'Compiler command  : ',a,                            &
     &        /,1x,'Compiler flags    : ',a,                            &
     &        /,1x,'OCN Communicator  : ',i0,',  PET size = ',i0,/,     &
     &        /,1x,'Input Script      : ',a,/,                          &
     &        /,1x,'SVN Root URL     : ',a,                             &
     &        /,1x,'SVN Revision     : ',a,/,                           &
     &        /,1x,'Local Root       : ',a,                             &
     &        /,1x,'Header Dir       : ',a,                             &
     &        /,1x,'Header file      : ',a,                             &
     &        /,1x,'Analytical Dir   : ',a)
  70  FORMAT (/,' Resolution, Grid ',i2.2,': ',i0,'x',i0,'x',i0,        &
     &        ',',2x,'Parallel Nodes: ',i0,',',2x,'Tiling: ',i0,        &
     &        'x',i0)
  80  FORMAT (/,' ROMS/TOMS: Wrong choice of grid ',i2.2,1x,            &
     &        'partition or number of parallel nodes.',                 &
     &        /,12x,a,1x,i0,/,12x,                                      &
     &        'must be equal to the number of parallel processes = ',   &
     &        i0,/,12x,'Change -np value to mpirun or',                 &
     &        /,12x,'change domain partition in input script.')
  90  FORMAT (/,' Resolution, Grid ',i2.2,': ',i0,'x',i0,'x',i0,        &
     &        ',',2x,'Parallel Threads: ',i0,',',2x,'Tiling: ',i0,      &
     &        'x',i0)
 100  FORMAT (/,' ROMS/TOMS: Wrong choice of grid ',i2.2,1x,            &
     &        'partition or number of parallel threads.',               &
     &        /,12x,'NtileI*NtileJ must be a positive multiple of the', &
     &        ' number of threads.',                                    &
     &        /,12x,'Change number of threads (environment variable) ', &
     &        'or',/,12x,'change domain partition in input script.')
 110  FORMAT (/,/,' Physical Parameters, Grid: ',i2.2,                  &
     &        /,  ' =============================',/)
 120  FORMAT (1x,i10,2x,a,t32,a)
 130  FORMAT (1x,i10,2x,a,t32,a,/,t34,a)
 140  FORMAT (f11.3,2x,a,t32,a)
 150  FORMAT (f11.2,2x,a,t32,a)
 160  FORMAT (f11.3,2x,a,t32,a,/,t34,a)
 170  FORMAT (10x,l1,2x,a,t32,a)
 180  FORMAT (10x,l1,2x,a,t32,a,i2.2,':',1x,a)
 185  FORMAT (10x,l1,2x,a,'(',i2.2,')',t32,a,i2.2,':',1x,a)
 190  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t32,a,/,t34,a,i2.2,':',1x,a)
 195  FORMAT (1p,e11.4,2x,a,t32,a,i2.2,':',1x,a)
 200  FORMAT (1p,e11.4,2x,a,t32,a)
 210  FORMAT (1p,e11.4,2x,a,t32,a,/,t34,a)
 220  FORMAT (/,' Output/Input Files:',/)
 230  FORMAT (2x,a,a)
 240  FORMAT (/,' Generic User Parameters:',/)
 250  FORMAT (1p,e11.4,2x,'user(',i2.2,')',t32,                         &
     &        'User parameter ',i2.2,'.')
 260  FORMAT (/,' READ_PHYPAR - Invalid input parameter, ',a,           &
     &        i4,/,15x,a)
 265  FORMAT (/,' READ_PHYPAR - Invalid input parameter, ',a,           &
     &        1p,e11.4,/,15x,a)
 270  FORMAT (/,' READ_PHYPAR - Grid ',i2.2,                            &
     &        ', could not find input file:  ',a)
 280  FORMAT (/,' READ_PHYPAR - Variable index not yet loaded, ', a)
 290  FORMAT (/,' READ_PHYPAR - Invalid dimension parameter, ',a,i0,    &
     &        /,15x,a)
 300  FORMAT (/,' READ_PHYPAR - Invalid dimension parameter, ',a,'(',   &
     &        i2.2,')',/,15x,a)
 310  FORMAT (2x,a,i2.2,a,a)
 320  FORMAT (/,' READ_PHYPAR - Could not find input parameter: ', a,   &
     &        /,15x,'in ROMS standard input script.',/,15x,a)
 330  FORMAT (/,' READ_PHYPAR - Invalid input parameter, ',a,i4,/,15x,a)
 340  FORMAT (/,' READ_PHYPAR - Inconsistent time-stepping period:',    &
     &        /,15x,'Grid ',i2.2,':',f14.1,' (sec)',2x,f14.2,' (days)', &
     &        /,15x,'Grid ',i2.2,':',f14.1,' (sec)',2x,f14.2,' (days)', &
     &        /,15x,'Adjust standard input parameter NTIMES in ',       &
     &              '''roms.in''.'/)
 350  FORMAT (/,' READ_PHYPAR - Invalid input parameter, ',a,i0,        &
     &        ', for grid ',i2.2,/,15x,a,i0,', ',a,i0,/,15x,a,/,15x,a)
      RETURN
      END SUBROUTINE read_PhyPar
