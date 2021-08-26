      MODULE wstvelocity_mod
!
!svn $Id: wstvelocity.F 732 2008-09-07 01:55:51Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutines computes vertical stokes velocity (m/s) at W-points!
!  from the vertical mass flux (omega*hz/m*n).  This computation       !
!  is done solely for output purposes.                                 !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: wstvelocity
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wstvelocity (ng, tile, Ninp)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Ninp
!
!  Local variable declarations.
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
      CALL wstvelocity_tile (ng, tile,                                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     Ninp,                                        &
     &                     GRID(ng) % pm,                               &
     &                     GRID(ng) % pn,                               &
     &                     GRID(ng) % z_r,                              &
     &                     GRID(ng) % z_w,                              &
     &                     GRID(ng) % h,                                &
     &                     GRID(ng) % om_u,                             &
     &                     GRID(ng) % on_u,                             &
     &                     GRID(ng) % om_v,                             &
     &                     GRID(ng) % on_v,                             &
     &                     OCEAN(ng) % ubar_stokes,                     &
     &                     OCEAN(ng) % vbar_stokes,                     &
     &                     OCEAN(ng) % u_stokes,                        &
     &                     OCEAN(ng) % v_stokes,                        &
     &                     OCEAN(ng) % W_stokes,                        &
     &                     OCEAN(ng) % wstvel)
      RETURN
      END SUBROUTINE wstvelocity
!
!***********************************************************************
      SUBROUTINE wstvelocity_tile (ng, tile,                            &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           Ninp,                                  &
     &                           pm, pn, z_r, z_w,                      &
     &                           h,                                     &
     &                           om_u, on_u, om_v, on_v,                &
     &                           ubar_stokes, vbar_stokes,              &
     &                           u_stokes, v_stokes, W_stokes,          &
     &                           wstvel)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_3d_mod, ONLY : bc_w3d_tile
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Ninp
!
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: ubar_stokes(LBi:,LBj:)
      real(r8), intent(in) :: vbar_stokes(LBi:,LBj:)
      real(r8), intent(in) :: u_stokes(LBi:,LBj:,:)
      real(r8), intent(in) :: v_stokes(LBi:,LBj:,:)
      real(r8), intent(in) :: W_stokes(LBi:,LBj:,0:)
      real(r8), intent(out):: wstvel(LBi:,LBj:,0:)
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8) :: cff1, cff2, cff3, cff4, cff5, slope
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: vert
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: DUst_avg1
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: DVst_avg1
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
!  Compute "true" vertical velocity (m/s).
!-----------------------------------------------------------------------
!
!  In ROMS, the terrain-following vertical velocity, omega, is given by:
!
!         Hz * omega = w - d(z)/d(t) - div(z)
!
!  where w is the "true" vertical velocity and
!
!         div(z) = pm * u * d(z)/d(xi) + pn * v * d(z)/d(eta)
!
!  The vertical coordinate is a function of several parameter but only
!  the free-surface is time dependent. However, in sediment applications
!  with stratigraphy, the bathymetry (h) also evolves in time.
!
!  Compute contribution due to quasi-horizontal motions along
!  S-coordinate surfaces:  (Usti + VStj)*GRADs(z).
!
!  This step is not performed for the stokes velocities because
!  their compuation is on a horizontal plane, and they are not
!  rotated to s-coordinates. This rotation is assumed to be
!  on the same order or smaller as the approximation to compute
!  the stokes vels.
!
      DO k=1,N(ng)
!       DO j=Jstr,Jend
!         DO i=Istr,Iend+1
!           wrk(i,j)=u_stokes(i,j,k)*(z_r(i,j,k)-z_r(i-1,j,k))*         &
!    &                             (pm(i-1,j)+pm(i,j))
!         END DO
!         DO i=Istr,Iend
!           vert(i,j,k)=0.25_r8*(wrk(i,j)+wrk(i+1,j))
!         END DO
!       END DO
!       DO j=Jstr,Jend+1
!         DO i=Istr,Iend
!           wrk(i,j)=v_stokes(i,j,k)*(z_r(i,j,k)-z_r(i,j-1,k))*         &
!    &                             (pn(i,j-1)+pn(i,j))
!         END DO
!       END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend
!           vert(i,j,k)=vert(i,j,k)+0.25_r8*(wrk(i,j)+wrk(i,j+1))
            vert(i,j,k)=0.0_r8
          END DO
        END DO
      END DO
!
! Calculate DUst_avg1 and DVst_avg1 here:
!
      DO j=Jstr,Jend
        DO i=Istr,Iend+1
          DUst_avg1(i,j)=0.5_r8*on_u(i,j)*ubar_stokes(i,j)*             &
     &                   ((z_w(i  ,j,N(ng))-z_w(i  ,j,0))+              &
     &                    (z_w(i-1,j,N(ng))-z_w(i-1,j,0)))
        END DO
      END DO
      DO j=Jstr,Jend+1
        DO i=Istr,Iend
          DVst_avg1(i,j)=0.5_r8*om_v(i,j)*vbar_stokes(i,j)*             &
     &                   ((z_w(i,j  ,N(ng))-z_w(i,j  ,0))+              &
     &                    (z_w(i,j-1,N(ng))-z_w(i,j-1,0)))
        END DO
      END DO
!
!  Compute contribution due to time tendency of the free-surface,
!  d(zeta)/d(t), which is the vertical velocity at the free-surface
!  and it is expressed in terms of barotropic mass flux divergence.
!  Notice that it is divided by the total depth of the water column.
!  This is needed because this contribution is linearly distributed
!  throughout the water column by multiplying it by the distance from
!  the bottom to the depth at which the vertical velocity is computed.
!
      cff1=3.0_r8/8.0_r8
      cff2=3.0_r8/4.0_r8
      cff3=1.0_r8/8.0_r8
      cff4=9.0_r8/16.0_r8
      cff5=1.0_r8/16.0_r8
      J_LOOP : DO j=Jstr,Jend
        DO i=Istr,Iend
          wrk(i,j)=(DUst_avg1(i,j)-DUst_avg1(i+1,j)+                    &
     &              DVst_avg1(i,j)-DVst_avg1(i,j+1))/                   &
     &             (z_w(i,j,N(ng))-z_w(i,j,0))
        END DO
!
!  Notice that a cubic interpolation is used to shift the "vert"
!  contribution from vertical RHO- to W-points.
!
        DO i=Istr,Iend
          slope=(z_r(i,j,1)-z_w(i,j,0))/                                &
     &          (z_r(i,j,2)-z_r(i,j,1))            ! extrapolation slope
          wstvel(i,j,0)=cff1*(vert(i,j,1)-                              &
     &                      slope*(vert(i,j,2)-                         &
     &                             vert(i,j,1)))+                       &
     &                cff2*vert(i,j,1)-                                 &
     &                cff3*vert(i,j,2)
          wstvel(i,j,1)=pm(i,j)*pn(i,j)*                                &
     &                (W_stokes(i,j,1)+                                 &
     &                 wrk(i,j)*(z_w(i,j,1)-z_w(i,j,0)))+               &
     &                cff1*vert(i,j,1)+                                 &
     &                cff2*vert(i,j,2)-                                 &
     &                cff3*vert(i,j,3)
        END DO
        DO k=2,N(ng)-2
          DO i=Istr,Iend
            wstvel(i,j,k)=pm(i,j)*pn(i,j)*                              &
     &                  (W_stokes(i,j,k)+                               &
     &                   wrk(i,j)*(z_w(i,j,k)-z_w(i,j,0)))+             &
     &                  cff4*(vert(i,j,k  )+vert(i,j,k+1))-             &
     &                  cff5*(vert(i,j,k-1)+vert(i,j,k+2))
          END DO
        END DO
        DO i=Istr,Iend
          slope=(z_w(i,j,N(ng))-z_r(i,j,N(ng)  ))/                      &
     &          (z_r(i,j,N(ng))-z_r(i,j,N(ng)-1))  ! extrapolation slope
          wstvel(i,j,N(ng))=pm(i,j)*pn(i,j)*                            &
     &                    wrk(i,j)*(z_w(i,j,N(ng))-z_w(i,j,0))+         &
     &                    cff1*(vert(i,j,N(ng))+                        &
     &                          slope*(vert(i,j,N(ng)  )-               &
     &                                 vert(i,j,N(ng)-1)))+             &
     &                    cff2*vert(i,j,N(ng)  )-                       &
     &                    cff3*vert(i,j,N(ng)-1)
          wstvel(i,j,N(ng)-1)=pm(i,j)*pn(i,j)*                          &
     &                      (W_stokes(i,j,N(ng)-1)+                     &
     &                       wrk(i,j)*(z_w(i,j,N(ng)-1)-z_w(i,j,0)))+   &
     &                      cff1*vert(i,j,N(ng)  )+                     &
     &                      cff2*vert(i,j,N(ng)-1)-                     &
     &                      cff3*vert(i,j,N(ng)-2)
        END DO
      END DO J_LOOP
!
!  Set lateral boundary conditions.
!
      CALL bc_w3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  wstvel)
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    wstvel)
      RETURN
      END SUBROUTINE wstvelocity_tile
      END MODULE wstvelocity_mod
