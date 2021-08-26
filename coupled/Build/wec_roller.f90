      MODULE wec_roller_mod
!
!svn $Id: wec_roller.F 1428 2008-03-12 13:07:21Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!                                                   Nirnimesh Kumar    !
!================================================== John C. Warner ====!
!                                                                      !
!  This routine computes the terms corresponding to vortex forces in   !
!  momentum equations.                                                 !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Svendsen, I.A., 1984. Mass flux and undertow in a surf zone.        !
!  Coastal Engineering 8, pp. 347–365.                              !
!                                                                      !
!  Reniers, A.J.M.H., Roelvink, J.A., and Thornton, E.B., 2004.        !
!  Morphodynamic modeling of an embayed beach under wave group forcing.!
!   J. Geophys. Res., 109: C01030, doi:10.1029/2002JC001586.           !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: wec_roller
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wec_roller (ng, tile)
!***********************************************************************
!
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
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
      CALL wec_roller_tile (ng, tile,    LBi, UBi, LBj, UBj, N(ng),     &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nrhs(ng),                             &
     &                            GRID(ng) % angler,                    &
     &                            GRID(ng) % h,                         &
     &                            GRID(ng) % Hz,                        &
     &                            GRID(ng) % on_u,                      &
     &                            GRID(ng) % om_v,                      &
     &                            GRID(ng) % pm,                        &
     &                            GRID(ng) % pn,                        &
     &                            OCEAN(ng) % ubar,                     &
     &                            OCEAN(ng) % vbar,                     &
     &                            OCEAN(ng) % zeta,                     &
     &                            FORCES(ng) % Hwave,                   &
     &                            FORCES(ng) % Dwave,                   &
     &                            FORCES(ng) % Lwave,                   &
     &                            FORCES(ng) % Dissip_break,            &
     &                            FORCES(ng) % Dissip_roller,           &
     &                            FORCES(ng) % rollA)
      CALL wclock_off (ng, iNLM, 21)
      RETURN
      END SUBROUTINE wec_roller
!
!***********************************************************************
      SUBROUTINE wec_roller_tile (ng, tile, LBi, UBi, LBj, UBj, UBk,    &
     &                                  IminS, ImaxS, JminS, JmaxS,     &
     &                                  nrhs,                           &
     &                                  angler, h, Hz,                  &
     &                                  on_u, om_v, pm, pn,             &
     &                                  ubar, vbar, zeta,               &
     &                                  Hwave, Dwave, Lwave,            &
     &                                  Dissip_break,                   &
     &                                  Dissip_roller,                  &
     &                                  rollA)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE exchange_2d_mod
      USE mp_exchange_mod, ONLY : mp_exchange2d
      USE bc_2d_mod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(in) :: Hwave(LBi:,LBj:)
      real(r8), intent(in) :: Dwave(LBi:,LBj:)
      real(r8), intent(in) :: Lwave(LBi:,LBj:)
      real(r8), intent(in) :: Dissip_break(LBi:,LBj:)
      real(r8), intent(inout) :: Dissip_roller(LBi:,LBj:)
      real(r8), intent(inout) :: rollA(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, j, k, numits, it
      real(r8) :: cff, cff1, cff2, cff3
      real(r8), parameter :: sinb=0.1_r8
      real(r8), parameter :: eps = 1.0E-14_r8
      real(r8), parameter :: kDmax = 5.0_r8
      real(r8), parameter :: Lwave_min = 1.0_r8
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Dstp
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: kD
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wavec
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waven
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: owaven
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wavenx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waveny
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: sigma
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: osigma
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: gamr
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
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
      DO j=Jstr,Jend+1
        DO i=Istr,Iend+1
!
!  Compute total depth
!
          Dstp(i,j)=zeta(i,j,1)+h(i,j)
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
          osigma(i,j)=1.0_r8/sigma(i,j)
!
!  Compute wave celerity and nonlinear water depth
!
          kD(i,j)=MIN(waven(i,j)*Dstp(i,j)+eps,kDmax)
          wavec(i,j)=SQRT(MAX(g*owaven(i,j)*TANH(kD(i,j)),eps))
        END DO
      END DO
!
!
!  Solve roller evolution equation for rollA.
!
      numits=30
      DO it=1,numits
!
!  Compute roller breaking source term (Eqn 40) and 
!  roller disspation sink term (Eqn 41).
!
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            cff3=(ubar(i,j,nrhs)+wavenx(i,j)*owaven(i,j)*               &
     &           wavec(i,j))*on_u(i,j)
            cff1=MAX(cff3,0.0_r8)
            cff2=MIN(cff3,0.0_r8)
            FX(i,j)=cff1*rollA(i-1,j)+cff2*rollA(i,j)
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            cff3=(vbar(i,j,nrhs)+waveny(i,j)*owaven(i,j)*               &
     &           wavec(i,j))*om_v(i,j)
            cff1=MAX(cff3,0.0_r8)
            cff2=MIN(cff3,0.0_r8)
            FE(i,j)=cff1*rollA(i,j-1)+cff2*rollA(i,j)
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff=dt(ng)*pm(i,j)*pn(i,j)/REAL(numits,r8)
            cff1=cff*(FX(i+1,j)-FX(i,j)+FE(i,j+1)-FE(i,j))
            rollA(i,j)=rollA(i,j)-cff1
          END DO
        END DO
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            Dissip_roller(i,j)=g*sinb*rollA(i,j)*sigma(i,j)/wavec(i,j) 
          END DO
        END DO
!
!  Add roller source / sink term.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff=dt(ng)/REAL(numits,r8)
            rollA(i,j)=rollA(i,j)+cff*osigma(i,j)*                      &
     &                 (wec_alpha(ng)*Dissip_break(i,j)-                &
     &                 Dissip_roller(i,j))
          END DO
        END DO
!
!  Call bc's.
!
        CALL bc_r2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    rollA)
        CALL mp_exchange2d (ng, tile, iNLM, 2,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      rollA)
      END DO
!
!  Apply boundary conditions.
!
        CALL bc_r2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Dissip_roller)
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Dissip_roller)
      RETURN
      END SUBROUTINE wec_roller_tile
      END MODULE wec_roller_mod
