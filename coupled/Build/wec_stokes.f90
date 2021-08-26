      MODULE wec_stokes_mod
!
!svn $Id: wec_stokes.F 1428 2008-03-12 13:07:21Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!                                                   Nirnimesh Kumar    !
!================================================== John C. Warner ====!
!                                                                      !
!  This routine computes the stokes transport terms                    !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: wec_stokes
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wec_stokes (ng, tile)
!***********************************************************************
!
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_param
!
      integer, intent(in) :: ng, tile
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private
!  storage arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J) directions and
!  MAX(I,J) directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
      CALL wclock_on (ng, iNLM, 21)
      CALL wec_stokes_tile (ng, tile, LBi, UBi, LBj, UBj, N(ng),        &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                            GRID(ng) % angler,                    &
     &                            GRID(ng) % h,                         &
     &                            GRID(ng) % rmask,                     &
     &                            GRID(ng) % umask,                     &
     &                            GRID(ng) % vmask,                     &
     &                            GRID(ng) % umask_wet,                 &
     &                            GRID(ng) % vmask_wet,                 &
     &                            GRID(ng) % on_u,                      &
     &                            GRID(ng) % on_v,                      &
     &                            GRID(ng) % om_u,                      &
     &                            GRID(ng) % om_v,                      &
     &                            GRID(ng) % Hz,                        &
     &                            GRID(ng) % z_r,                       &
     &                            GRID(ng) % z_w,                       &
     &                            FORCES(ng) % Hwave,                   &
     &                            FORCES(ng) % Dwave,                   &
     &                            FORCES(ng) % Lwave,                   &
     &                            FORCES(ng) % rollA,                   &
     &                            OCEAN(ng) % zeta,                     &
     &                            OCEAN(ng) % ubar_stokes,              &
     &                            OCEAN(ng) % vbar_stokes,              &
     &                            OCEAN(ng) % u_stokes,                 &
     &                            OCEAN(ng) % v_stokes,                 &
     &                            OCEAN(ng) % W_stokes)
      CALL wclock_off (ng, iNLM, 21)
      RETURN
      END SUBROUTINE wec_stokes
!
!***********************************************************************
      SUBROUTINE wec_stokes_tile (ng, tile, LBi, UBi, LBj, UBj, UBk,    &
     &                                  IminS, ImaxS, JminS, JmaxS,     &
     &                                  angler, h,                      &
     &                                  rmask, umask, vmask,            &
     &                                  umask_wet, vmask_wet,           &
     &                                  on_u, on_v, om_u, om_v,         &
     &                                  Hz, z_r, z_w,                   &
     &                                  Hwave, Dwave, Lwave,            &
     &                                  rollA,                          &
     &                                  zeta,                           &
     &                                  ubar_stokes, vbar_stokes,       &
     &                                  u_stokes, v_stokes, W_stokes)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE exchange_2d_mod
      USE exchange_3d_mod
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
      USE bc_2d_mod
      USE bc_3d_mod
      USE us2dbc_mod, ONLY : us2dbc_tile
      USE vs2dbc_mod, ONLY : vs2dbc_tile
      USE us3dbc_mod, ONLY : us3dbc_tile
      USE vs3dbc_mod, ONLY : vs3dbc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: Hwave(LBi:,LBj:)
      real(r8), intent(in) :: Dwave(LBi:,LBj:)
      real(r8), intent(in) :: Lwave(LBi:,LBj:)
      real(r8), intent(in) :: rollA(LBi:,LBj:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(inout) :: ubar_stokes(LBi:,LBj:)
      real(r8), intent(inout) :: vbar_stokes(LBi:,LBj:)
      real(r8), intent(inout) :: u_stokes(LBi:,LBj:,:)
      real(r8), intent(inout) :: v_stokes(LBi:,LBj:,:)
      real(r8), intent(inout) :: W_stokes(LBi:,LBj:,0:)
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8) :: cff, cff2, cff3, cff4, cff5, cff6, cff7
      real(r8) :: fac1, fac2, fac3, ofac3
      real(r8), parameter :: eps = 1.0E-14_r8
      real(r8), parameter :: kDmax = 200.0_r8
      real(r8), parameter :: Lwave_min = 1.0_r8
      real(r8), dimension(IminS:ImaxS) :: wrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Dstp
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: kD
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wavec
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waven
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: owaven
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wavenx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waveny
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waveE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: sigma
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Huons
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hvoms
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =BOUNDS(ng) % Istr   (tile)
      IstrB  =BOUNDS(ng) % IstrB  (tile)
      IstrM  =BOUNDS(ng) % IstrM  (tile)
      IstrP  =BOUNDS(ng) % IstrP  (tile)
      IstrR  =BOUNDS(ng) % IstrR  (tile)
      IstrT  =BOUNDS(ng) % IstrT  (tile)
      IstrU  =BOUNDS(ng) % IstrU  (tile)
      Iend   =BOUNDS(ng) % Iend   (tile)
      IendB  =BOUNDS(ng) % IendB  (tile)
      IendP  =BOUNDS(ng) % IendP  (tile)
      IendR  =BOUNDS(ng) % IendR  (tile)
      IendT  =BOUNDS(ng) % IendT  (tile)
      Jstr   =BOUNDS(ng) % Jstr   (tile)
      JstrB  =BOUNDS(ng) % JstrB  (tile)
      JstrM  =BOUNDS(ng) % JstrM  (tile)
      JstrP  =BOUNDS(ng) % JstrP  (tile)
      JstrR  =BOUNDS(ng) % JstrR  (tile)
      JstrT  =BOUNDS(ng) % JstrT  (tile)
      JstrV  =BOUNDS(ng) % JstrV  (tile)
      Jend   =BOUNDS(ng) % Jend   (tile)
      JendB  =BOUNDS(ng) % JendB  (tile)
      JendP  =BOUNDS(ng) % JendP  (tile)
      JendR  =BOUNDS(ng) % JendR  (tile)
      JendT  =BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =BOUNDS(ng) % Iendp3 (tile)            ! Iend+3
      Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
      fac1=1.0_r8/dt(ng)
      fac2=1.0_r8/g
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
!
!  Compute total depth
!
          Dstp(i,j)=z_w(i,j,N(ng))-z_w(i,j,0)
!
!  Compute wave amplitude (0.5*Hrms), wave number, intrinsic frequency.
!
          waven(i,j)=2.0_r8*pi/MAX(Lwave(i,j),Lwave_min)
          owaven(i,j)=1.0_r8/waven(i,j)
          cff=1.5_r8*pi-Dwave(i,j)-angler(i,j)
          wavenx(i,j)=waven(i,j)*COS(cff)
          waveny(i,j)=waven(i,j)*SIN(cff)
          sigma(i,j)=SQRT(MAX(g*waven(i,j)*TANH(waven(i,j)*Dstp(i,j)),  &
     &                    eps))
          waveE(i,j)=0.0625_r8*g*Hwave(i,j)*Hwave(i,j)
!
!  Compute wave celerity and kD
!
          kD(i,j)=MIN(waven(i,j)*Dstp(i,j),kDmax)
          wavec(i,j)=SQRT(MAX(g*owaven(i,j)*TANH(kD(i,j)),eps))
        END DO
      END DO
!
!---------------------------------------------------------------------------
! Stokes velocities.
!---------------------------------------------------------------------------
!
!  Compute u-stokes velocities, Eqn. 2.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff2=(waveE(i-1,j)+waveE(i,j))
          cff2=cff2+(rollA(i-1,j)*sigma(i-1,j)+                         &

     &               rollA(i,j)*sigma(i,j))
          cff3=(kD(i-1,j)+kD(i,j))
!         cff3=(waven(i-1,j)+waven(i,j))*(Dstp(i-1,j)+Dstp(i,j))
          cff4=wavenx(i-1,j)+wavenx(i,j)
          cff5=wavec(i-1,j)+wavec(i,j)
          fac3=Dstp(i-1,j)+Dstp(i,j)
          ofac3=1.0_r8/fac3
          DO k=1,N(ng)
            cff6=-1.0_r8+((z_w(i-1,j,k)+z_w(i,j,k))-                    &

     &                   (z_w(i-1,j,0)+z_w(i,j,0)))*ofac3
            cff7=-1.0_r8+((z_w(i-1,j,k-1)+z_w(i,j,k-1))-                &

     &                   (z_w(i-1,j,0)+z_w(i,j,0)))*ofac3
            u_stokes(i,j,k)=0.25_r8*cff2*cff4*cff5*fac2/                &

     &                      (fac3*(cff6-cff7))/                         &

     &                      (SINH(0.5_r8*cff3)**2.0_r8)*                &

     &                      SINH(0.5_r8*cff3*(cff6-cff7))*              &

     &                      COSH(0.5_r8*cff3*(cff6+cff7+2.0_r8))
            u_stokes(i,j,k)=u_stokes(i,j,k)*umask(i,j)
            u_stokes(i,j,k)=u_stokes(i,j,k)*umask_wet(i,j)
          END DO
        END DO
      END DO
!
!  Compute v-stokes velocity, Eqn. 2.
!
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff2=(waveE(i,j-1)+waveE(i,j))
          cff2=cff2+(rollA(i,j-1)*sigma(i,j-1)+                         &

     &               rollA(i,j)*sigma(i,j))
          cff3=(kD(i,j-1)+kD(i,j))
!         cff3=(waven(i,j-1)+waven(i,j))*(Dstp(i,j-1)+Dstp(i,j))
          cff4=waveny(i,j-1)+waveny(i,j)
          cff5=wavec(i,j-1)+wavec(i,j)
          fac3=Dstp(i,j-1)+Dstp(i,j)
          ofac3=1.0_r8/fac3
          DO k=1,N(ng)
            cff6=-1.0_r8+((z_w(i,j-1,k  )+z_w(i,j,k))-                  &
     &                    (z_w(i,j-1,0  )+z_w(i,j,0)))*ofac3
            cff7=-1.0_r8+((z_w(i,j-1,k-1)+z_w(i,j,k-1))-                &
     &                    (z_w(i,j-1,0  )+z_w(i,j,0)))*ofac3
            v_stokes(i,j,k)=0.25_r8*cff2*cff4*cff5*fac2/                &
     &                      (fac3*(cff6-cff7))/                         &
     &                      (SINH(0.5_r8*cff3)**2.0_r8)*                &
     &                      SINH(0.5_r8*cff3*(cff6-cff7))*              &
     &                      COSH(0.5_r8*cff3*(cff6+cff7+2.0_r8))
            v_stokes(i,j,k)=v_stokes(i,j,k)*vmask(i,j)
            v_stokes(i,j,k)=v_stokes(i,j,k)*vmask_wet(i,j)
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
! Set lateral boundary conditions.
!-----------------------------------------------------------------------
!
      CALL us3dbc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, N(ng),                      &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  u_stokes)
      CALL vs3dbc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, N(ng),                      &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  v_stokes)
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!
      CALL exchange_u3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        u_stokes(:,:,:))
      CALL exchange_v3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        v_stokes(:,:,:))
      END IF
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    u_stokes, v_stokes)
!
!  Compute vertical stokes velocity, Eqn. 31.
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
             Huons(i,j)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*                 &
     &                  u_stokes(i,j,k)*on_u(i,j)
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
             Hvoms(i,j)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*                 &
     &                  v_stokes(i,j,k)*om_v(i,j)
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend
             W_stokes(i,j,k)=W_stokes(i,j,k-1)-                         &
     &                      (Huons(i+1,j)-Huons(i,j)+                   &
     &                       Hvoms(i,j+1)-Hvoms(i,j))
          END DO
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=Istr,Iend
          wrk(i)=W_stokes(i,j,N(ng))/(z_w(i,j,N(ng))-z_w(i,j,0))
        END DO
        DO k=N(ng)-1,1,-1
         DO i=Istr,Iend
            W_stokes(i,j,k)=W_stokes(i,j,k)-                            &
     &                      wrk(i)*(z_w(i,j,k)-z_w(i,j,0))
          END DO
        END DO
        DO i=Istr,Iend
          W_stokes(i,j,N(ng))=0.0_r8
        END DO
      END DO
!
!  For a 3D application, compute associated 2D fields by taking the
!  vertical integral of 3D fields.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff5=0.5_r8*(Hz(i-1,j,1)+Hz(i,j,1))
          ubar_stokes(i,j)=cff5*u_stokes(i,j,1)
          DO k=2,N(ng)
            cff5=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))
            ubar_stokes(i,j)=ubar_stokes(i,j)+cff5*u_stokes(i,j,k)
          END DO
          cff4=2.0_r8/(Dstp(i-1,j)+Dstp(i,j))
          ubar_stokes(i,j)=ubar_stokes(i,j)*cff4
        END DO
      END DO
      DO i=Istr,Iend
        DO j=JstrV,Jend
          cff5=0.5_r8*(Hz(i,j-1,1)+Hz(i,j,1))
          vbar_stokes(i,j)=cff5*v_stokes(i,j,1)
          DO k=2,N(ng)
            cff5=0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))
            vbar_stokes(i,j)=vbar_stokes(i,j)+cff5*v_stokes(i,j,k)
          END DO
          cff4=2.0_r8/(Dstp(i,j-1)+Dstp(i,j))
          vbar_stokes(i,j)=vbar_stokes(i,j)*cff4
        END DO
      END DO
!
!  Apply boundary conditions.
!
        CALL us2dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    ubar_stokes)
        CALL vs2dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    vbar_stokes)
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ubar_stokes(:,:))
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        vbar_stokes(:,:))
      END IF
      CALL bc_w3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  W_stokes)
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    ubar_stokes, vbar_stokes)
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    W_stokes)
      RETURN
      END SUBROUTINE wec_stokes_tile
      END MODULE wec_stokes_mod
