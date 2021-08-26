      MODULE bbl_mod
!
!svn $Id: bbl.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes bottom momentum stress via a bottom boundary  !
!  layer formulation.                                                  !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: bblm
      CONTAINS
      SUBROUTINE bblm (ng, tile)
!
!svn $Id: ssw_bbl.h 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group        Chris Sherwood   !
!    Licensed under a MIT/X style license               Rich Signell   !
!    See License_ROMS.txt                             John C. Warner   !
!=======================================================================
!                                                                      !
!  This routine compute bottom stresses for the case when the wave     !
!  solution in the wave boundary layer is based on a  2-layer eddy     !
!  viscosity that is linear increasing above Zo and constant above     !
!  Z1.                                                                 !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Styles, R. and S.M. glenn,  2000: Modeling stratified wave and      !
!    current bottom boundary layers in the continental shelf, JGR,     !
!    105, 24119-24139.                                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_bbl
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_sedbed
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Nonlinear/ssw_bbl.h"
!
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
!
      CALL wclock_on (ng, iNLM, 37, 48, MyFile)
      CALL bblm_tile (ng, tile,                                         &
     &                LBi, UBi, LBj, UBj,                               &
     &                IminS, ImaxS, JminS, JmaxS,                       &
     &                nrhs(ng),                                         &
     &                GRID(ng) % h,                                     &
     &                GRID(ng) % z_r,                                   &
     &                GRID(ng) % z_w,                                   &
     &                GRID(ng) % angler,                                &
     &                GRID(ng) % ZoBot,                                 &
     &                FORCES(ng) % Uwave_rms,                           &
     &                FORCES(ng) % Dwave,                               &
     &                FORCES(ng) % Pwave_bot,                           &
     &                SEDBED(ng) % bottom,                              &
     &                OCEAN(ng) % rho,                                  &
     &                OCEAN(ng) % u,                                    &
     &                OCEAN(ng) % v,                                    &
     &                OCEAN(ng) % u_stokes,                             &
     &                OCEAN(ng) % v_stokes,                             &
     &                BBL(ng) % Iconv,                                  &
     &                BBL(ng) % Ubot,                                   &
     &                BBL(ng) % Vbot,                                   &
     &                BBL(ng) % Ur,                                     &
     &                BBL(ng) % Vr,                                     &
     &                BBL(ng) % bustrc,                                 &
     &                BBL(ng) % bvstrc,                                 &
     &                BBL(ng) % bustrw,                                 &
     &                BBL(ng) % bvstrw,                                 &
     &                BBL(ng) % bustrcwmax,                             &
     &                BBL(ng) % bvstrcwmax,                             &
     &                FORCES(ng) % bustr,                               &
     &                FORCES(ng) % bvstr)
      CALL wclock_off (ng, iNLM, 37, 105, MyFile)
!
      RETURN
      END SUBROUTINE bblm
!
!***********************************************************************
      SUBROUTINE bblm_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs,                                       &
     &                      h, z_r, z_w, angler, ZoBot,                 &
     &                      Uwave_rms,                                  &
     &                      Dwave, Pwave_bot,                           &
     &                      bottom, rho,                                &
     &                      u, v,                                       &
     &                      u_stokes, v_stokes,                         &
     &                      Iconv,                                      &
     &                      Ubot, Vbot, Ur, Vr,                         &
     &                      bustrc, bvstrc,                             &
     &                      bustrw, bvstrw,                             &
     &                      bustrcwmax, bvstrcwmax,                     &
     &                      bustr, bvstr)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
      USE mod_sediment
!
      USE bc_2d_mod
      USE mp_exchange_mod, ONLY : mp_exchange2d
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
      integer, intent(inout) :: Iconv(LBi:,LBj:)
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: ZoBot(LBi:,LBj:)
      real(r8), intent(in) :: Uwave_rms(LBi:,LBj:)
      real(r8), intent(in) :: Dwave(LBi:,LBj:)
      real(r8), intent(in) :: Pwave_bot(LBi:,LBj:)
      real(r8), intent(inout) :: bottom(LBi:,LBj:,:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: u_stokes(LBi:,LBj:,:)
      real(r8), intent(in) :: v_stokes(LBi:,LBj:,:)
      real(r8), intent(out) :: Ubot(LBi:,LBj:)
      real(r8), intent(out) :: Vbot(LBi:,LBj:)
      real(r8), intent(out) :: Ur(LBi:,LBj:)
      real(r8), intent(out) :: Vr(LBi:,LBj:)
      real(r8), intent(out) :: bustrc(LBi:,LBj:)
      real(r8), intent(out) :: bvstrc(LBi:,LBj:)
      real(r8), intent(out) :: bustrw(LBi:,LBj:)
      real(r8), intent(out) :: bvstrw(LBi:,LBj:)
      real(r8), intent(out) :: bustrcwmax(LBi:,LBj:)
      real(r8), intent(out) :: bvstrcwmax(LBi:,LBj:)
      real(r8), intent(out) :: bustr(LBi:,LBj:)
      real(r8), intent(out) :: bvstr(LBi:,LBj:)
!
!  Local variable declarations.
!
      logical :: ITERATE
      integer :: Iter, i, j, k
      real(r8), parameter :: eps = 1.0E-10_r8
      real(r8) :: Kbh, Kbh2, Kdh
      real(r8) :: taucr, wsedr, tstar, coef_st
      real(r8) :: coef_b1, coef_b2, coef_b3, d0
      real(r8) :: dolam, dolam1, doeta1, doeta2, fdo_etaano
      real(r8) :: lamorb, lamanorb
      real(r8) :: m_ubr, m_wr, m_ucr, m_zr, m_phicw, m_kb
      real(r8) :: m_ustrc, m_ustrwm, m_ustrr, m_fwc, m_zoa, m_dwc
      real(r8) :: zo, Dstp
      real(r8) :: Kb, Kdelta, Ustr
      real(r8) :: anglec, anglew
      real(r8) :: cff, cff1, cff2, cff3, og, fac, fac1, fac2
      real(r8) :: sg_ab, sg_abokb, sg_a1, sg_b1, sg_chi, sg_c1, d50
      real(r8) :: sg_epsilon, ssw_eta, sg_fofa, sg_fofb, sg_fofc, sg_fwm
      real(r8) :: sg_kbs, ssw_lambda, sg_mu, sg_phicw, sg_ro, sg_row
      real(r8) :: sg_shdnrm, sg_shld, sg_shldcr, sg_scf, rhos, sg_star
      real(r8) :: sg_ub, sg_ubokur, sg_ubouc, sg_ubouwm, sg_ur
      real(r8) :: sg_ustarc, sg_ustarcw, sg_ustarwm, sg_znot, sg_znotp
      real(r8) :: sg_zr, sg_zrozn, sg_z1, sg_z1ozn, sg_z2, z1, z2
      real(r8) :: zoMIN, zoMAX
      real(r8) :: coef_fd
      real(r8), parameter :: twopi=2.0_r8*pi
!
      real(r8), parameter :: absolute_zoMIN = 5.0d-5  ! in Harris-Wiberg
      real(r8), parameter :: Cd_fd = 0.5_r8
      real(r8), parameter :: K1 = 0.6666666666_r8     ! Coefficients for
      real(r8), parameter :: K2 = 0.3555555555_r8     ! explicit
      real(r8), parameter :: K3 = 0.1608465608_r8     ! wavenumber
      real(r8), parameter :: K4 = 0.0632098765_r8     ! calculation
      real(r8), parameter :: K5 = 0.0217540484_r8     ! (Dean and
      real(r8), parameter :: K6 = 0.0065407983_r8     !  Dalrymple, 1991)
      real(r8), parameter :: coef_a1=0.095_r8         ! Coefficients for
      real(r8), parameter :: coef_a2=0.442_r8         ! ripple predictor
      real(r8), parameter :: coef_a3=2.280_r8         ! (Wiberg-Harris)
      real(r8), parameter :: ar = 0.267_r8         ! Nielsen (1992)
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Ab
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Fwave_bot
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Tauc
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Tauw
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Taucwmax
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Ur_sg
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Vr_sg
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Ub
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Ucur
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Umag
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Vcur
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Zr
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: phic
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: phicw
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rheight
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rlength
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: u100
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: znot
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: znotc
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: zoN
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: zoST
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: zoBF
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: zoDEF
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: zoBIO
!
      real(r8), dimension(1:N(ng)) :: Urz, Vrz
!
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
!
!-----------------------------------------------------------------------
!  Set currents above the bed.
!-----------------------------------------------------------------------
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
!
! Calculate bottom cell thickness
!
          Zr(i,j)=z_r(i,j,1)-z_w(i,j,0)
!
!
          Dstp=z_r(i,j,N(ng))-z_w(i,j,0)
!
!
! Use the original coded ssw_logint formulation
!
          cff1=sg_z1min
!         cff1=sg_zwbl(ng)
!
!  If using the logarithmic interpolation
!
          DO k=1,N(ng)
            Urz(k)=0.5_r8*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
            Vrz(k)=0.5_r8*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
            Urz(k)=Urz(k)+0.5_r8*(u_stokes(i,j,k)+u_stokes(i+1,j,k))
            Vrz(k)=Vrz(k)+0.5_r8*(v_stokes(i,j,k)+v_stokes(i,j+1,k))
          END DO
          CALL log_interp( N(ng), Dstp, cff1,                           &
     &                 Urz, Vrz,                                        &
     &                 z_r(i,j,:),        z_w(i,j,:),                   &
     &                 bottom(i,j,isd50), bottom(i,j,izapp),            &
     &                 Zr(i,j),                                         &
     &                 Ur_sg(i,j),        Vr_sg(i,j) )
!
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute bottom stresses.
!-----------------------------------------------------------------------
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
!
!  Set bed wave orbital velocity and excursion amplitude.  Use data
!  from wave models (SWAN) or use Dean and Dalrymple (1991) 6th-degree
!  polynomial to approximate wave number on shoaling water.
          Fwave_bot(i,j)=twopi/MAX(Pwave_bot(i,j),1.0_r8)
          Ub(i,j)=MAX(Uwave_rms(i,j),0.0_r8)+eps
          Ab(i,j)=Ub(i,j)/Fwave_bot(i,j)+eps
!
!  Compute bottom current magnitude at RHO-points.
!
          Ucur(i,j)=Ur_sg(i,j)
          Vcur(i,j)=Vr_sg(i,j)
!
          Umag(i,j)=SQRT(Ucur(i,j)*Ucur(i,j)+Vcur(i,j)*Vcur(i,j)+eps)
!
!  Compute angle between currents and waves (radians)
!
          IF (Ucur(i,j).eq.0.0_r8) THEN
            phic(i,j)=0.5_r8*pi*SIGN(1.0_r8,Vcur(i,j))
          ELSE
            phic(i,j)=ATAN2(Vcur(i,j),Ucur(i,j))
          ENDIF
          phicw(i,j)=1.5_r8*pi-Dwave(i,j)-phic(i,j)-angler(i,j)
!
        END DO
      END DO
!
!  Loop over RHO points.
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
!
!  Load sediment properties, stresses, and roughness from previous time
!  step (stresses are in m2 s-2).
!
          d50=bottom(i,j,isd50)
          rhos=bottom(i,j,idens)/(rho(i,j,1)+1000.0_r8)
          wsedr=bottom(i,j,iwsed)
          Taucr=bottom(i,j,itauc)
          Tauc(i,j)=SQRT(bustrc(i,j)**2+bvstrc(i,j)**2)
          Tauw(i,j)=SQRT(bustrw(i,j)**2+bvstrw(i,j)**2)
          Taucwmax(i,j)=SQRT( bustrcwmax(i,j)**2+bvstrcwmax(i,j)**2)
!
          rheight(i,j)=bottom(i,j,irhgt)
          rlength(i,j)=bottom(i,j,irlen)
          zoMAX=0.9_r8*Zr(i,j)
          zoMIN=MAX(absolute_zoMIN,2.5_r8*d50/30.0_r8)
!
!  Initialize arrays.
!
          zoN(i,j)=MIN(MAX(2.5_r8*d50/30.0_r8, zoMIN ),zoMAX)
          zoST(i,j)=0.0_r8
          zoBF(i,j)=0.0_r8
          zoBIO(i,j)=0.0_r8
          zoDEF(i,j)=ZoBot(i,j)
!
!  Calculate components of roughness and sum: zo = zoN + zoST + zoBF
!  Determine whether sediment is in motion. Use Shields criterion to
!  determine if sediment is mobile.
!
          tstar=Taucwmax(i,j)/(Taucr+eps)
          IF (tstar.lt.1.0_r8) THEN                         ! no motion
            zoST(i,j)=0.0_r8
            zoBF(i,j)=ar*rheight(i,j)**2/(rlength(i,j)+eps)
          ELSE
!
!  Threshold of motion exceeded - calculate new zoST and zoBF
!  Calculate saltation roughness according to Wiberg & Rubin (1989)
!  (Eqn. 11 in Harris & Wiberg, 2001)
!  (d50 is in m, but this formula needs cm)
!
             coef_st=0.0204_r8*LOG(100.0_r8*d50+eps)**2+                &
     &               0.0220_r8*LOG(100.0_r8*d50+eps)+0.0709_r8
             zoST(i,j)=0.056_r8*d50*0.68_r8*tstar/                      &
     &                 (1.0_r8+coef_st*tstar)
             IF (zoST(i,j).lt.0.0_r8) THEN
               IF (Master) THEN
                 PRINT *, ' Warning: zoST<0  tstar, d50, coef_st:'
                 PRINT *, tstar,d50,coef_st
               END IF
             END IF
!
!  Calculate ripple height and wavelength.
!  Use Malarkey & Davies (2003) explict version of Wiberg & Harris.
!
             coef_b1=1.0_r8/coef_a1
             coef_b2=0.5_r8*(1.0_r8 + coef_a2)*coef_b1
             coef_b3=coef_b2**2-coef_a3*coef_b1
             d0=2.0_r8*Ab(i,j)
             IF ((d0/d50).gt.13000.0_r8) THEN              ! sheet flow
               rheight(i,j)=0.0_r8
               rlength(i,j)=535.0_r8*d50        ! does not matter since
             ELSE                               ! rheight=0
               dolam1=d0/(535.0_r8*d50)
               doeta1=EXP(coef_b2-SQRT(coef_b3-coef_b1*LOG(dolam1)))
               lamorb=0.62_r8*d0
               lamanorb=535.0_r8*d50
               IF (doeta1.lt.20.0_r8) THEN
                 dolam=1.0_r8/0.62_r8
               ELSE IF (doeta1.gt.100.0_r8) THEN
                 dolam=dolam1
               ELSE
                 fdo_etaano=-LOG(lamorb/lamanorb)*                      &
     &                       LOG(0.01_r8*doeta1)/LOG(5.0_r8)
                 dolam=dolam1*EXP(-fdo_etaano)
               END IF
               doeta2=EXP(coef_b2-SQRT(coef_b3-coef_b1*LOG(dolam)))
               rheight(i,j)=d0/doeta2
               rlength(i,j)=d0/dolam
             END IF
!
!  Value of ar can range from 0.3 to 3 (Soulsby, 1997, p. 124)
!
             zoBF(i,j)=ar*rheight(i,j)**2/rlength(i,j)
          END IF
          zo=zoN(i,j)
!
!  Compute stresses.
!
!  Default stress calcs for pure currents
!
          zo=MIN(MAX(zo,zoMIN),zoMAX)
!
          cff1=vonKar/LOG(Zr(i,j)/zo)
          cff2=MIN(Cdb_max,MAX(Cdb_min,cff1*cff1))
          Tauc(i,j)=cff2*Umag(i,j)*Umag(i,j)
          Tauw(i,j)=0.0_r8
          Taucwmax(i,j)=Tauc(i,j)
          znot(i,j)=zo
          znotc(i,j)=zo
!
          IF ((Umag(i,j).le.eps).and.(Ub(i,j).ge.eps)) THEN
!
!  Pure waves - use wave friction factor approach from Madsen
!  (1994, eqns 32-33).
!
            sg_abokb=Ab(i,j)/(30.0_r8*zo)
            sg_fwm=0.3_r8
            IF ((sg_abokb.gt.0.2_r8).and.(sg_abokb.le.100.0_r8)) THEN
              sg_fwm=EXP(-8.82_r8+7.02_r8*sg_abokb**(-0.078_r8))
            ELSE IF (sg_abokb.gt.100.0_r8)THEN
              sg_fwm=EXP(-7.30_r8+5.61_r8*sg_abokb**(-0.109_r8))
            END IF
            Tauc(i,j)= 0.0_r8
            Tauw(i,j)= 0.5_r8*sg_fwm*Ub(i,j)*Ub(i,j)
            Taucwmax(i,j)=Tauw(i,j)
            znot(i,j)=zo
            znotc(i,j)=zo
!         ELSE IF ((Umag(i,j).gt.0.0_r8).and.(Ub(i,j).ge.eps).and.      &
          ELSE IF ((Umag(i,j).gt.eps).and.(Ub(i,j).ge.eps).and.         &
     &             ((Zr(i,j)/zo).le.1.0_r8)) THEN
!
!  Waves and currents, but zr <= zo.
!
            IF (Master) THEN
              PRINT *,' Warning: w-c calcs ignored because zr <= zo'
            END IF
!         ELSE IF ((Umag(i,j).gt.0.0_r8).and.(Ub(i,j).gt.eps).and.      &
          ELSE IF ((Umag(i,j).gt.eps).and.(Ub(i,j).ge.eps).and.         &
     &             ((Zr(i,j)/zo).gt.1.0_r8)) THEN
!
!  Waves and currents, zr > zo.
!
            m_ubr=Ub(i,j)
            m_wr=Fwave_bot(i,j)
            m_ucr=Umag(i,j)
            m_zr=Zr(i,j)
            m_phicw=phicw(i,j)
            m_kb=30.0_r8*zo
            Dstp=z_r(i,j,N(ng))-z_w(i,j,0)
            CALL madsen94 (m_ubr, m_wr, m_ucr,                          &
     &                     m_zr, m_phicw, m_kb, Dstp,                   &
     &                     m_ustrc, m_ustrwm, m_ustrr, m_fwc, m_zoa,    &
     &                     m_dwc)
            Tauc(i,j)=m_ustrc*m_ustrc
            Tauw(i,j)=m_ustrwm*m_ustrwm
            Taucwmax(i,j)=m_ustrr*m_ustrr
            znotc(i,j)=min( m_zoa, zoMAX )
            u100(i,j)=(m_ustrc/vonKar)*LOG(1.0_r8/m_zoa)
          END IF
        END DO
      END DO
!
!
!-----------------------------------------------------------------------
!  Compute kinematic bottom stress components due current and wind-
!  induced waves.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          anglec=0.5_r8*(Ur_sg(i,j)+Ur_sg(i-1,j))/(0.5_r8*(Umag(i-1,j)+Umag(i,j)))
          bustr(i,j)=0.5_r8*(Tauc(i-1,j)+Tauc(i,j))*anglec
          cff2=0.75_r8*0.5_r8*(z_w(i-1,j,1)+z_w(i,j,1)-                 &
     &                         z_w(i-1,j,0)-z_w(i,j,0))
          bustr(i,j)=SIGN(1.0_r8,bustr(i,j))*MIN(ABS(bustr(i,j)),       &
     &               ABS(u(i,j,1,nrhs))*cff2/dt(ng))
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          anglec=0.5_r8*(Vr_sg(i,j)+Vr_sg(i,j-1))/(0.5_r8*(Umag(i,j-1)+Umag(i,j)))
          bvstr(i,j)=0.5_r8*(Tauc(i,j-1)+Tauc(i,j))*anglec
          cff2=0.75_r8*0.5_r8*(z_w(i,j-1,1)+z_w(i,j,1)-                 &
     &                         z_w(i,j-1,0)-z_w(i,j,0))
          bvstr(i,j)=SIGN(1.0_r8,bvstr(i,j))*MIN(ABS(bvstr(i,j)),       &
     &               ABS(v(i,j,1,nrhs))*cff2/dt(ng))
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=Istr,Iend
          anglec=Ucur(i,j)/Umag(i,j)
          anglew=COS(1.5_r8*pi-Dwave(i,j)-angler(i,j))
          bustrc(i,j)=Tauc(i,j)*anglec
          bustrw(i,j)=Tauw(i,j)*anglew
          bustrcwmax(i,j)=Taucwmax(i,j)*anglew
          Ubot(i,j)=Ub(i,j)*anglew
          Ur(i,j)=Ucur(i,j)
!
          anglec=Vcur(i,j)/Umag(i,j)
          anglew=SIN(1.5_r8*pi-Dwave(i,j)-angler(i,j))
          bvstrc(i,j)=Tauc(i,j)*anglec
          bvstrw(i,j)=Tauw(i,j)*anglew
          bvstrcwmax(i,j)=Taucwmax(i,j)*anglew
          Vbot(i,j)=Ub(i,j)*anglew
          Vr(i,j)=Vcur(i,j)
!
          bottom(i,j,ibwav)=Ab(i,j)
          bottom(i,j,irhgt)=rheight(i,j)
          bottom(i,j,irlen)=rlength(i,j)
          bottom(i,j,izdef)=zoDEF(i,j)
          bottom(i,j,izapp)=znotc(i,j)
          bottom(i,j,izNik)=zoN(i,j)
          bottom(i,j,izbio)=zoBIO(i,j)
          bottom(i,j,izbfm)=zoBF(i,j)
          bottom(i,j,izbld)=zoST(i,j)
          bottom(i,j,izwbl)=znot(i,j)
        END DO
      END DO
!
!  Apply periodic or gradient boundary conditions for output
!  purposes only.
!
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bustr)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bvstr)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bustrc)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bvstrc)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bustrw)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bvstrw)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bustrcwmax)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bvstrcwmax)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Ubot)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Vbot)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Ur)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Vr)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,ibwav))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,irhgt))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,irlen))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,izdef))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,izapp))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,izNik))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,izbio))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,izbfm))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,izbld))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,izwbl))
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bustr, bvstr, bustrc, bvstrc)
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bustrw, bvstrw, bustrcwmax, bvstrcwmax)
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Ubot, Vbot, Ur, Vr)
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bottom(:,:,ibwav),                            &
     &                    bottom(:,:,irhgt),                            &
     &                    bottom(:,:,irlen))
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bottom(:,:,izdef),                            &
     &                    bottom(:,:,izapp))
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bottom(:,:,izNik),                            &
     &                    bottom(:,:,izbio),                            &
     &                    bottom(:,:,izbfm),                            &
     &                    bottom(:,:,izwbl))
      RETURN
      END SUBROUTINE bblm_tile
      SUBROUTINE madsen94 (ubr, wr, ucr, zr, phiwc, kN, Dstp,           &
     &                     ustrc, ustrwm, ustrr, fwc, zoa, dwc_va)
!
!=======================================================================
!                                                                      !
!  Grant-Madsen model from Madsen (1994).                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ubr     Rep. wave-orbital velocity amplitude outside WBL (m/s).  !
!     wr      Rep. angular wave frequency,  2* pi/T (rad/s).           !
!     ucr     Current velocity at height zr (m/s).                     !
!     zr      Reference height for current velocity (m).               !
!     phiwc   Angle between currents and waves at zr (radians).        !
!     kN      Bottom roughness height, like Nikuradse k, (m).          !
!     Dstp    Total water depth (m).                                   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ustrc   Current friction velocity, u*c (m/s).                    !
!     ustrwm  Wave maximum friction velocity, u*wm (m/s).              !
!     ustrr   Wave-current combined friction velocity, u*r (m/s).      !
!     fwc     Wave friction factor (nondimensional).                   !
!     zoa     Apparent bottom roughness (m).                           !
!     dwc_va  Wave-boundary layer thickness (m).                       !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      real(r8), intent(in)    :: ubr, wr, ucr, zr, phiwc, Dstp
      real(r8), intent(inout) :: kN
      real(r8), intent(out)   :: ustrc, ustrwm, ustrr, fwc, zoa
      real(r8), intent(out)   :: dwc_va
!
!  Local variable declarations.
!
      integer, parameter :: MAXIT = 20
      integer :: i, nit
      integer :: iverbose = 1
      real(r8) :: bigsqr, cosphiwc, cukw, diff, lndw, lnln, lnzr
      real(r8) :: phicwc, zo
      real(r8) :: dval = 99.99_r8
      real(r8), dimension(MAXIT) :: Cmu
      real(r8), dimension(MAXIT) :: dwc
      real(r8), dimension(MAXIT) :: fwci
      real(r8), dimension(MAXIT) :: rmu
      real(r8), dimension(MAXIT) :: ustrci
      real(r8), dimension(MAXIT) :: ustrr2
      real(r8), dimension(MAXIT) :: ustrwm2
!
!-----------------------------------------------------------------------
!  Compute bottom friction velocities and roughness.
!-----------------------------------------------------------------------
!
!  Set special default values.
!
      ustrc=dval
      ustrwm=dval
      ustrr=dval
      kN=MIN(kN,0.9_r8*zr)
      zoa=kN/30.0_r8
      phicwc=phiwc
      zo = kN/30.0_r8
      IF (ubr.le.0.01_r8) THEN
        dwc_va=kN
        zoa=dwc_va
        fwc=0.0_r8
        IF (ucr.le. 0.01_r8) THEN          ! no waves or currents
          ustrc=0.0_r8
          ustrwm=0.0_r8
          ustrr=0.0_r8
          RETURN
        END IF
        ustrc=ucr*vonKar/LOG(zr/zo)        ! no waves
        ustrwm=0.0_r8
        ustrr=ustrc
        RETURN
      END IF
!
!  Iterate to compute friction velocities, roughness, and wave friction
!  factor.  Notice that the computation of the wave friction factor
!  has been inlined for efficiency.
!
      cosphiwc=ABS(COS(phiwc))
      rmu(1)=0.0_r8
      Cmu(1)=1.0_r8
      cukw=Cmu(1)*ubr/(kN*wr)
!
! New fwc CRS calculation
!
      fwci(1)=Cmu(1)*0.3_r8
      IF ((cukw.gt.0.352_r8).and.(cukw.le.100.0_r8)) THEN       ! Eq 32/33
        fwci(1)=Cmu(1)*EXP(7.02_r8*cukw**(-0.078_r8)-8.82_r8)
      ELSE IF (cukw.gt.100.0_r8) THEN
        fwci(1)=Cmu(1)*EXP(5.61_r8*cukw**(-0.109_r8)-7.30_r8)
      END IF
!
      ustrwm2(1)=0.5_r8*fwci(1)*ubr*ubr                       ! Eq 29
      ustrr2(1)=Cmu(1)*ustrwm2(1)                             ! Eq 26
      ustrr=SQRT(ustrr2(1))
      IF (cukw.ge.8.0_r8) THEN
        dwc(1)=2.0_r8*vonKar*ustrr/wr
        dwc(1)=MIN( 0.9_r8*zr, dwc(1) )
      ELSE
        dwc(1)=kN
      END IF
      lnzr=LOG(zr/dwc(1))
      lndw=LOG(dwc(1)/zo)
      lnln=lnzr/lndw
      bigsqr=-1.0_r8+SQRT(1.0_r8+((4.0_r8*vonKar*lndw)/                 &
     &                            (lnzr*lnzr))*ucr/ustrr)
      ustrci(1)=0.5_r8*ustrr*lnln*bigsqr
!
      i=1
      diff=1.0_r8
      DO WHILE ((i.lt.MAXIT).and.(diff.gt.0.000005_r8))
        i=i+1
        rmu(i)=ustrci(i-1)*ustrci(i-1)/ustrwm2(i-1)
        Cmu(i)=SQRT(1.0_r8+                                             &
     &              2.0_r8*rmu(i)*cosphiwc+rmu(i)*rmu(i))     ! Eq 27
        cukw=Cmu(i)*ubr/(kN*wr)
!
! New fwc CRS calculation
!
        fwci(i)=Cmu(i)*0.3_r8
        IF ((cukw.gt.0.352_r8).and.(cukw.le.100.0_r8)) THEN       ! Eq 32/33
          fwci(i)=Cmu(i)*EXP(7.02_r8*cukw**(-0.078_r8)-8.82_r8)
        ELSE IF (cukw.gt.100.0_r8) THEN
          fwci(i)=Cmu(i)*EXP(5.61_r8*cukw**(-0.109_r8)-7.30_r8)
        END IF
!
        ustrwm2(i)=0.5_r8*fwci(i)*ubr*ubr                     ! Eq 29
        ustrr2(i)=Cmu(i)*ustrwm2(i)                           ! Eq 26
        ustrr=SQRT(ustrr2(i))
        IF (cukw.ge.8.0_r8) THEN
          dwc(i)=2.0_r8*vonKar*ustrr/wr                       ! Eq 36
          dwc(i)=MIN( 0.9_r8*zr, dwc(i) )
        ELSE
          dwc(i)=kN
        END IF
        lnzr=LOG(zr/dwc(i))
        lndw=LOG(dwc(i)/zo)
        lnln=lnzr/lndw
        bigsqr=-1.0_r8+SQRT(1.0_r8+((4.0_r8*vonKar*lndw)/               &
     &                              (lnzr*lnzr))*ucr/ustrr)
        ustrci(i)=0.5_r8*ustrr*lnln*bigsqr                    ! Eq 38
        diff=ABS((fwci(i)-fwci(i-1))/fwci(i))
      END DO
      ustrwm=SQRT(ustrwm2(i))
      ustrc=ustrci(i)
      ustrr=SQRT(ustrr2(i))
      phicwc=phiwc
      zoa=EXP(LOG(dwc(i))-(ustrc/ustrr)*LOG(dwc(i)/zo))       ! Eq 11
      dwc_va=dwc(i)
      fwc=fwci(i)
      RETURN
      END SUBROUTINE madsen94
!
!
      SUBROUTINE log_interp( kmax, Dstp, sg_loc, u_1d, v_1d,            &
     &                        z_r_1d, z_w_1d,                           &
     &                        d50, zapp_loc,                            &
     &                        Zr_sg,                                    &
     &                        Ur_sg, Vr_sg)
!
!=======================================================================
!   Find the near-bottom current velocity in x, y dir. that corresponds!
!   to user input elevation to get near-bottom current vel. (m)        !
!                                                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_sediment
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer,  intent(in)  :: kmax
      real(r8), intent(in)  :: Dstp, sg_loc
      real(r8), intent(in)  :: u_1d(1:kmax), v_1d(1:kmax)
      real(r8), intent(in)  :: z_r_1d(1:kmax), z_w_1d(0:kmax)
      real(r8), intent(in)  :: d50, zapp_loc
      real(r8), intent(out) :: Zr_sg
      real(r8), intent(out) :: Ur_sg, Vr_sg
!
!  Local variables.
!
      integer  :: k
      real(r8) :: z1, z2, Zr
      real(r8) :: fac, fac1, fac2
!
      Zr=z_r_1d(1)-z_w_1d(0)
!
      IF ( sg_loc.ge.Zr ) THEN
!
!  If chosen height to get near bottom-current velocity lies
!  within any vertical level, perform logarithmic interpolation.
!
        DO k=2,kmax
          z1=z_r_1d(k-1)-z_w_1d(0)
          z2=z_r_1d(k  )-z_w_1d(0)
          IF ( ( z1.le.sg_loc ).and.( sg_loc.lt.z2 )) THEN
            fac=1.0_r8/LOG(z2/z1)
            fac1=fac*LOG(z2/sg_loc)
            fac2=fac*LOG(sg_loc/z1)
!
            Ur_sg=fac1*u_1d(k-1)+fac2*u_1d(k)
            Vr_sg=fac1*v_1d(k-1)+fac2*v_1d(k)
            Zr_sg=sg_loc
          ENDIF
          IF ((k.eq.kmax).and.(sg_loc.ge.z2)) THEN
            Ur_sg=u_1d(k)
            Vr_sg=v_1d(k)
            Zr_sg=z2
          END IF
        END DO
      ELSE    !IF ( sg_loc.lt.Zr ) THEN
        z1=MAX( 2.5_r8*d50/30.0_r8, zapp_loc, 1.0e-10_r8 )
        z2=Zr
!
        IF ( sg_loc.lt.z1 ) THEN
!
! If chosen height is less than the bottom roughness
! perform linear interpolation.
!
          z1=sg_loc
          fac=z1/z2
!
          Ur_sg=fac*u_1d(1)
          Vr_sg=fac*v_1d(1)
          Zr_sg=sg_loc
!
        ELSE   !IF ( sg_loc.gt.z1 ) THEN
!
! If chosen height is less than the bottom cell thickness
! perform logarithmic interpolation with bottom roughness.
!
          fac=1.0_r8/LOG(z2/z1)
          fac2=fac*LOG(sg_loc/z1)
!
          Ur_sg=fac2*u_1d(1)
          Vr_sg=fac2*v_1d(1)
          Zr_sg=sg_loc
        END IF
      END IF
!
      RETURN
      END SUBROUTINE log_interp
      END MODULE bbl_mod
