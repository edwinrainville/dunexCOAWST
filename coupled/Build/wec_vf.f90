      MODULE wec_vf_mod
!
!svn $Id: wec_vf.F 1428 2008-03-12 13:07:21Z jcwarner $
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
!  Uchiyama, Y., McWilliams, J.C., and Shchepetkin, A.F. (2010).       !
!  Wave current interacation in an oceanic circulation model with a    !
!  vortex-force formalism: Applications to surf zone, Ocean Modeling,  !
!  34, 16-35.                                                          !
!                                                                      !
!  Kumar, N., Voulgaris, G., Warner, J.C., and M., Olabarrieta (2012). !
!  Implementation of a vortex force formalism in the coupled           !
!  ocean-atmosphere-wave-sediment transport (COAWST) modeling system   !
!  for inner-shelf and surf-zone applications.                         !
!  Ocean Modeling 47, pp 65-95.                                        !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: wec_vf
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wec_vf (ng, tile)
!***********************************************************************
!
      USE mod_forces
      USE mod_grid
      USE mod_mixing
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
      CALL wec_vf_tile (ng, tile, LBi, UBi, LBj, UBj, N(ng),            &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nrhs(ng),                             &
     &                            GRID(ng) % pmask,                     &
     &                            GRID(ng) % rmask,                     &
     &                            GRID(ng) % umask,                     &
     &                            GRID(ng) % vmask,                     &
     &                            GRID(ng) % rmask_wet,                 &
     &                            GRID(ng) % umask_wet,                 &
     &                            GRID(ng) % vmask_wet,                 &
     &                            GRID(ng) % om_u,                      &
     &                            GRID(ng) % om_v,                      &
     &                            GRID(ng) % on_u,                      &
     &                            GRID(ng) % on_v,                      &
     &                            GRID(ng) % pm,                        &
     &                            GRID(ng) % pn,                        &
     &                            GRID(ng) % angler,                    &
     &                            GRID(ng) % h,                         &
     &                            GRID(ng) % Hz,                        &
     &                            GRID(ng) % z_r,                       &
     &                            GRID(ng) % z_w,                       &
     &                            FORCES(ng) % Hwave,                   &
     &                            FORCES(ng) % Dwave,                   &
     &                            FORCES(ng) % Lwave,                   &
     &                            FORCES(ng) % Dissip_break,            &
     &                            FORCES(ng) % Dissip_roller,           &
     &                            FORCES(ng) % Dissip_wcap,             &
     &                            OCEAN(ng) % zeta,                     &
     &                            OCEAN(ng) % u,                        &
     &                            OCEAN(ng) % v,                        &
     &                            OCEAN(ng) % u_stokes,                 &
     &                            OCEAN(ng) % v_stokes,                 &
     &                            OCEAN(ng) % bh,                       &
     &                            OCEAN(ng) % qsp,                      &
     &                            OCEAN(ng) % zetaw,                    &
     &                            OCEAN(ng) % zetat,                    &
     &                            MIXING(ng) % rustr3d,                 &
     &                            MIXING(ng) % rvstr3d,                 &
     &                            MIXING(ng) % rurol2d,                 &
     &                            MIXING(ng) % rvrol2d,                 &
     &                            MIXING(ng) % rubrk2d,                 &
     &                            MIXING(ng) % rvbrk2d,                 &
     &                            MIXING(ng) % rukvf2d,                 &
     &                            MIXING(ng) % rvkvf2d)
      CALL wclock_off (ng, iNLM, 21)
      RETURN
      END SUBROUTINE wec_vf
!
!***********************************************************************
      SUBROUTINE wec_vf_tile (ng, tile, LBi, UBi, LBj, UBj, UBk,        &
     &                                  IminS, ImaxS, JminS, JmaxS,     &
     &                                  nrhs,                           &
     &                                  pmask, rmask, umask, vmask,     &
     &                                  rmask_wet, umask_wet, vmask_wet,&
     &                                  om_u, om_v, on_u, on_v,         &
     &                                  pm, pn,                         &
     &                                  angler,                         &
     &                                  h,                              &
     &                                  Hz, z_r, z_w,                   &
     &                                  Hwave, Dwave, Lwave,            &
     &                                  Dissip_break,                   &
     &                                  Dissip_roller,                  &
     &                                  Dissip_wcap,                    &
     &                                  zeta,                           &
     &                                  u, v,                           &
     &                                  u_stokes, v_stokes,             &
     &                                  bh, qsp, zetaw, zetat,          &
     &                                  rustr3d,  rvstr3d,              &
     &                                  rurol2d, rvrol2d,               &
     &                                  rubrk2d, rvbrk2d,               &
     &                                  rukvf2d, rvkvf2d)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
      USE bc_2d_mod
      USE bc_3d_mod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
      real(r8), intent(in) :: pmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: Hwave(LBi:,LBj:)
      real(r8), intent(in) :: Dwave(LBi:,LBj:)
      real(r8), intent(in) :: Lwave(LBi:,LBj:)
      real(r8), intent(in) :: Dissip_break(LBi:,LBj:)
      real(r8), intent(in) :: Dissip_roller(LBi:,LBj:)
      real(r8), intent(in) :: Dissip_wcap(LBi:,LBj:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: u_stokes(LBi:,LBj:,:)
      real(r8), intent(in) :: v_stokes(LBi:,LBj:,:)
      real(r8), intent(inout) :: bh(LBi:,LBj:)
      real(r8), intent(inout) :: qsp(LBi:,LBj:)
      real(r8), intent(inout) :: zetaw(LBi:,LBj:)
      real(r8), intent(inout) :: zetat(LBi:,LBj:)
      real(r8), intent(inout) :: rustr3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: rvstr3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: rurol2d(LBi:,LBj:)
      real(r8), intent(inout) :: rvrol2d(LBi:,LBj:)
      real(r8), intent(inout) :: rubrk2d(LBi:,LBj:)
      real(r8), intent(inout) :: rvbrk2d(LBi:,LBj:)
      real(r8), intent(inout) :: rukvf2d(LBi:,LBj:)
      real(r8), intent(inout) :: rvkvf2d(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5, cff6
      real(r8) :: fac1, fac2, RB1, RB2
      real(r8) :: vi1, vi2, VF
!     real(r8), parameter :: eps = 1.0E-14_r8
      real(r8), parameter :: eps = 1.0E-8_r8
      real(r8), parameter :: kDmax = 300.0_r8                !5.0_r8
      real(r8), parameter :: Lwave_min = 1.0_r8
      real(r8), dimension(IminS:ImaxS) :: CF
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Dstp
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: oDstp
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: kD
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wavec
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waven
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: owaven
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wavenx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waveny
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waveAA
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: sigma
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: sigmat
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: osigma
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: gamr
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: odiss
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: BBXbrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: BBYbrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: BBXrol
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: BBYrol
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ku
      real(r8), dimension(IminS:ImaxS,N(ng)) :: kv
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: shear1
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: shear2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hzk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Huston
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hvston
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hustom
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hvstom
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: shearuv
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
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
!
!  Compute total depth
!
          Dstp(i,j)=z_w(i,j,N(ng))-z_w(i,j,0)
          oDstp(i,j)=1.0_r8/Dstp(i,j)
!
!  Compute wave amplitude (0.5*Hrms), wave number, intrinsic frequency.
!
          waven(i,j)=2.0_r8*pi/MAX(Lwave(i,j),Lwave_min)
          owaven(i,j)=1.0_r8/waven(i,j)
          cff=1.5_r8*pi-Dwave(i,j)-angler(i,j)
          wavenx(i,j)=waven(i,j)*COS(cff)
          waveny(i,j)=waven(i,j)*SIN(cff)
          cff=0.25_r8*SQRT(2.0_r8)*Hwave(i,j)
          waveAA(i,j)=cff*cff
          waveAA(i,j)=waveAA(i,j)*rmask(i,j)
          waveAA(i,j)=waveAA(i,j)*rmask_wet(i,j)
          sigma(i,j)=SQRT(MAX(g*waven(i,j)*TANH(waven(i,j)*Dstp(i,j)),  &
     &                    eps))
          osigma(i,j)=1.0_r8/sigma(i,j)
!
!  Compute wave celerity and nonlinear water depth
!
          kD(i,j)=MIN(waven(i,j)*Dstp(i,j), kDmax)
          wavec(i,j)=SQRT(MAX(g*owaven(i,j)*TANH(kD(i,j)),eps))
!
!  Compute metrics for vertical dissipation distribution.
!
          odiss(i,j)=0.0_r8
!         gamr(i,j)=MAX(MIN(0.707_r8*Dstp(i,j)/                         &
          gamr(i,j)=MAX(MIN(2.0_r8*0.707_r8*Dstp(i,j)/                  &
     &              (1.25_r8*Hwave(i,j)+eps),1.0_r8),0.0_r8)
          DO k=1,N(ng)
            cff2=MAX(MIN((z_r(i,j,k)-z_w(i,j,0))*oDstp(i,j)*gamr(i,j),  &
     &           1.0_r8),0.0_r8)
            odiss(i,j)=odiss(i,j)+Hz(i,j,k)*                            &
     &                   COSH(2.0_r8*pi*cff2)
          END DO
          odiss(i,j)=1.0_r8/(odiss(i,j)+eps)
!
!  Initialize depth independent arrays for summation
!
        END DO
      END DO
!
!---------------------------------------------------------------------------
!  Compute Bernoulli's Head (BH), represented as the 
!  symbol cursive K in Eqn. 5.
!---------------------------------------------------------------------------
!
      DO j=Jstr,Jend
!
! Compute k*u and k*v at rho pts.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            ku(i,k)=0.5_r8*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))*             &
                    wavenx(i,j)
            kv(i,k)=0.5_r8*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))*             &
                    waveny(i,j)
          END DO
        END DO
!
!  Compute gradients for shear d(k_x*u)/dz at cell faces.
!
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            cff=1.0_r8/(z_r(i,j,k+1)-z_r(i,j,k))
            shear1(i,k)=cff*((ku(i,k+1)-ku(i,k))+                       &
     &                       (kv(i,k+1)-kv(i,k)))
          END DO
        END DO
!
! Applying boundary conditions.
!
        DO i=Istr,Iend
          shear1(i,N(ng))=shear1(i,N(ng)-1)
          shear1(i,0)=shear1(i,1)
        END DO
!
! Compute second derivative d(k_x*u)^2/dz^2
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            cff=1.0_r8/(z_w(i,j,k)-z_w(i,j,k-1))
            shear2(i,k)=cff*(shear1(i,k)-shear1(i,k-1))
          END DO
        END DO
!
!  Perform integration on the shear2 term calcuated above.
!
        DO i=Istr,Iend
          cff=Hz(i,j,1)
          fac2=(z_r(i,j,N(ng))-z_r(i,j,1))*oDstp(i,j)
          vi1=cff*shear2(i,1)*SINH(2.0_r8*kD(i,j)*fac2)
          DO k=2,N(ng)
            cff=Hz(i,j,k)
            fac2=(z_r(i,j,N(ng))-z_r(i,j,k))*oDstp(i,j)
            vi1=vi1+cff*shear2(i,k)*SINH(2.0_r8*kD(i,j)*fac2)             
          END DO
!
!  Calculate Bernoulli's Head at the surface
!
          bh(i,j)=0.25_r8*vi1*sigma(i,j)*waveAA(i,j)/                   & 
     &                    (waven(i,j)*(SINH(kD(i,j)))**2.0_r8)
        END DO
!
!---------------------------------------------------------------------------
!  Compute quasi-static components of pressure terms (i.e. Cursive P)
!  in Eqn. 9.
!---------------------------------------------------------------------------
!
        DO i=Istr,Iend
          fac2=(z_r(i,j,1)-z_w(i,j,0))*oDstp(i,j)
          cff=Hz(i,j,1)
          vi2=cff*shear2(i,1)*COSH(2.0_r8*kD(i,j)*fac2)
          DO k=2,N(ng)
             fac2=(z_r(i,j,k)-z_w(i,j,0))*oDstp(i,j)
             cff=Hz(i,j,k)
             vi2=vi2+cff*shear2(i,k)*COSH(2.0_r8*kD(i,j)*fac2)
          END DO
          cff1=0.5_r8*waveAA(i,j)*osigma(i,j)
          cff2=TANH(kD(i,j))/SINH(2.0_r8*kD(i,j))
          qsp(i,j)=cff1*(cff2*(-shear1(i,N(ng))+                        &
     &                         COSH(2.0_r8*kD(i,j))*shear1(i,0)+        &
     &                         vi2)-                                    &
     &                   2.0_r8*waven(i,j)*TANH(kD(i,j))*               &
     &                   (ku(i,N(ng))+kv(i,N(ng))))
        END DO
      END DO
!
!---------------------------------------------------------------------------
! Evaluate geopotential function (Eqn. 33, Uchiyama et al., 2010)
!---------------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
!
!  Compute quasi static sea surface elevation, Eq(7).
!
          zetaw(i,j)=-(0.5_r8*waveAA(i,j)*waven(i,j)/                   &
     &               SINH(2.0_r8*kD(i,j)))
!
!  Compute barotropic component of geopotential function, without zeta_c.
!  This gets sent to prsgrd functions.  The zeta_c is already included in
!  those algorithms. Also send individual zetaw, qsp, and bh to step2d for
!  computation of total pressure gradient.
!
          zetat(i,j)=-g*zetaw(i,j)-(g*qsp(i,j)-bh(i,j))
        END DO
      END DO
!
!---------------------------------------------------------------------------
! Compute Vertical Vortex Force terms, denoted as K in Eqn. 5.
! This next section needs to be re-written for piped J.
!---------------------------------------------------------------------------
!
      J_LOOP : DO j=JstrV-1,Jend
! Compute avg u and v at rho pts. Notice reuse of ku and kv.
!
        DO k=1,N(ng)
          DO i=IstrU-1,Iend
            ku(i,k)=0.5_r8*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
            kv(i,k)=0.5_r8*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
          END DO
        END DO
!
!  Compute gradients for shear du/dz at cell faces. notice reuse of 
!  shear1 and shear2.
!
        DO k=1,N(ng)-1
          DO i=IstrU-1,Iend
            cff=1.0_r8/(z_r(i,j,k+1)-z_r(i,j,k))
            shear1(i,k)=cff*(ku(i,k+1)-ku(i,k))
            shear2(i,k)=cff*(kv(i,k+1)-kv(i,k))
          END DO
        END DO
        DO i=IstrU-1,Iend
          shear1(i,0)=0.0_r8
          shear1(i,N(ng))=0.0_r8
          shear2(i,0)=0.0_r8
          shear2(i,N(ng))=0.0_r8
        END DO
!
!  Compute avg of shear1 at rho points.
!
        DO k=1,N(ng)
          DO i=IstrU-1,Iend
            cff1=0.5_r8*(u_stokes(i,j,k)+u_stokes(i+1,j  ,k))
            cff2=0.5_r8*(v_stokes(i,j,k)+v_stokes(i  ,j+1,k))
            shearuv(i,j,k)=(0.5_r8*cff1*(shear1(i,k-1)+shear1(i,k))+    &
     &                     0.5_r8*cff2*(shear2(i,k-1)+shear2(i,k)))*    &
     &                     Hz(i,j,k)
          END DO
        END DO
      END DO J_LOOP
!
! Compute contribution to rhs tersm.
!
      DO j=Jstr,Jend
        DO k=1,N(ng)
          DO i=IstrU,Iend
            Hzk(i,k)=0.5_r8*(Hz(i-1,j,k)+                               &
     &                       Hz(i  ,j,k))
          END DO
        END DO
        DO i=IstrU,Iend
          cff=Hzk(i,1)*on_u(i,j)
          VF=cff*(shearuv(i,j,1)-shearuv(i-1,j,1))
          rustr3d(i,j,1)=-VF
          rukvf2d(i,j)=VF
        END DO
        DO k=2,N(ng)
          DO i=IstrU,Iend
            cff=Hzk(i,k)*on_u(i,j)
            VF=cff*(shearuv(i,j,k)-shearuv(i-1,j,k))
            rustr3d(i,j,k)=-VF
            rukvf2d(i,j)=rukvf2d(i,j)+VF
          END DO
        END DO
        IF (j.ge.JstrV) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              Hzk(i,k)=0.5_r8*(Hz(i,j  ,k)+                             &
     &                         Hz(i,j-1,k))
            END DO
          END DO
          DO i=Istr,Iend
            cff=Hzk(i,1)*om_v(i,j)
            VF=cff*(shearuv(i,j,1)-shearuv(i,j-1,1))
            rvstr3d(i,j,1)=-VF
            rvkvf2d(i,j)=VF
          END DO
          DO k=2,N(ng)
            DO i=Istr,Iend
              cff=Hzk(i,k)*om_v(i,j)
              VF=cff*(shearuv(i,j,k)-shearuv(i,j-1,k))
              rvstr3d(i,j,k)=-VF
              rvkvf2d(i,j)=rvkvf2d(i,j)+VF
            END DO
          END DO
        END IF
      END DO
!
!------------------------------------------------------------------------
! Compute non conservative wave acceleration terms here
!------------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          rubrk2d(i,j)=0.0_r8
          rvbrk2d(i,j)=0.0_r8
          rurol2d(i,j)=0.0_r8
          rvrol2d(i,j)=0.0_r8
        END DO
      END DO
      DO k=1,N(ng)
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            fac2=(z_r(i,j,k)-z_w(i,j,0))*oDstp(i,j)
            cff2=fac2*gamr(i,j)
            cff3=COSH(2.0_r8*pi*cff2)*odiss(i,j)
!           cff1=Dissip_wcap(i,j)+Dissip_break(i,j)
            cff1=Dissip_break(i,j)
            cff4=(1.0_r8-wec_alpha(ng))*cff1*osigma(i,j)
            BBXbrk(i,j)=cff3*cff4*wavenx(i,j)
            BBYbrk(i,j)=cff3*cff4*waveny(i,j)
            cff5=Dissip_roller(i,j)*osigma(i,j)
            BBXrol(i,j)=cff3*cff5*wavenx(i,j)
            BBYrol(i,j)=cff3*cff5*waveny(i,j)
          END DO
        END DO
!
! Compute contribution to U-momentum
!
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=0.5_r8*(BBXbrk(i  ,j)*Hz(i  ,j,k)+                      &
     &                  BBXbrk(i-1,j)*Hz(i-1,j,k))*                     &
     &                  om_u(i,j)*on_u(i,j)
            cff=cff*umask(i,j)
            cff=cff*umask_wet(i,j)
            rustr3d(i,j,k)=rustr3d(i,j,k)-cff
            rubrk2d(i,j)=rubrk2d(i,j)+cff
            cff=0.5_r8*(BBXrol(i  ,j)*Hz(i  ,j,k)+                      &
     &                  BBXrol(i-1,j)*Hz(i-1,j,k))*                     &
     &                  om_u(i,j)*on_u(i,j)
            cff=cff*umask(i,j)
            cff=cff*umask_wet(i,j)
            rustr3d(i,j,k)=rustr3d(i,j,k)-cff
            rurol2d(i,j)=rurol2d(i,j)+cff
          END DO
          IF (j.ge.JstrV) THEN
!
! Compute contribution to V-momentum
!
            DO i=Istr,Iend
              cff=0.5_r8*(BBYbrk(i,j  )*Hz(i,j  ,k)+                    &
     &                    BBYbrk(i,j-1)*Hz(i,j-1,k))*                   &
     &                    om_v(i,j)*on_v(i,j)
              cff=cff*vmask(i,j)
              cff=cff*vmask_wet(i,j)
              rvstr3d(i,j,k)=rvstr3d(i,j,k)-cff
              rvbrk2d(i,j)=rvbrk2d(i,j)+cff
              cff=0.5_r8*(BBYrol(i,j  )*Hz(i,j  ,k)+                    &
     &                    BBYrol(i,j-1)*Hz(i,j-1,k))*                   &
     &                    om_v(i,j)*on_v(i,j)
              cff=cff*vmask(i,j)
              cff=cff*vmask_wet(i,j)
              rvstr3d(i,j,k)=rvstr3d(i,j,k)-cff
              rvrol2d(i,j)=rvrol2d(i,j)+cff
            END DO
          END IF
        END DO
      END DO
      K_LOOP : DO k=1,N(ng)
!
!---------------------------------------------------------------------------
! Contribution of a term corresponding to product of 
! Stokes and Eulerian Velocity Eqn. 26 and 27.
! This removes terms that were unneccessarily added in flux form.
!---------------------------------------------------------------------------
!
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))
            Huston(i,j)=cff*on_u(i,j)*u_stokes(i,j,k)
            Hvston(i,j)=0.25_r8*cff*on_u(i,j)*                          &
     &                  (v_stokes(i  ,j  ,k)+                           &
     &                   v_stokes(i  ,j+1,k)+                           &
     &                   v_stokes(i-1,j  ,k)+                           &
     &                   v_stokes(i-1,j+1,k))
          END DO
          DO i=IstrU-1,Iend
            UFx(i,j)=0.5_r8*(u(i  ,j  ,k,nrhs)+                         &
                             u(i+1,j  ,k,nrhs))
            VFx(i,j)=0.5_r8*(v(i  ,j  ,k,nrhs)+                         &
     &                       v(i  ,j+1,k,nrhs))
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))
            Hustom(i,j)=cff*0.25_r8*om_v(i,j)*                          &
     &                  (u_stokes(i  ,j  ,k)+                           &
     &                   u_stokes(i+1,j  ,k)+                           &
     &                   u_stokes(i  ,j-1,k)+                           &
     &                   u_stokes(i+1,j-1,k))
            Hvstom(i,j)=cff*om_v(i,j)*v_stokes(i,j,k)
          END DO
        END DO
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
            cff=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))
            UFe(i,j)=0.5_r8*(u(i+1,j  ,k,nrhs)+                         &
     &                       u(i  ,j  ,k,nrhs))
            VFe(i,j)=0.5_r8*(v(i  ,j  ,k,nrhs)+                         &
     &                       v(i  ,j+1,k,nrhs))
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff1=UFx(i,j)-UFx(i-1,j)
            cff2=VFx(i,j)-VFx(i-1,j)
            cff3=Huston(i,j)*cff1
            cff4=Hvston(i,j)*cff2
            rustr3d(i,j,k)=rustr3d(i,j,k)-cff3-cff4
!
!  Conversion to m2/s2.
!
            cff1=0.25_r8*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
            rustr3d(i,j,k)=rustr3d(i,j,k)*cff1
          END DO
        END DO
        DO i=Istr,Iend
          DO j=JstrV,Jend 
            cff1=UFe(i,j)-UFe(i,j-1)
            cff2=VFe(i,j)-VFe(i,j-1)
            cff3=Hustom(i,j)*cff1
            cff4=Hvstom(i,j)*cff2 
            rvstr3d(i,j,k)=rvstr3d(i,j,k)-cff3-cff4
!
!  Conversion to m2/s2.
!
            cff1=0.25_r8*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
            rvstr3d(i,j,k)=rvstr3d(i,j,k)*cff1
          END DO
        END DO
      END DO K_LOOP
!
!  Apply boundary conditions.
!
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  rubrk2d)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  rvbrk2d)
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  rurol2d)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  rvrol2d)
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  rukvf2d)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  rvkvf2d)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bh)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  qsp)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  zetat)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  zetaw)
      CALL bc_u3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  rustr3d)
      CALL bc_v3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  rvstr3d)
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rubrk2d, rvbrk2d)
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rurol2d, rvrol2d)
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rukvf2d, rvkvf2d)
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bh, qsp, zetat, zetaw)
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rustr3d, rvstr3d)
      RETURN
      END SUBROUTINE wec_vf_tile
      END MODULE wec_vf_mod
