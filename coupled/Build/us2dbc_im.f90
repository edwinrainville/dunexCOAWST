      MODULE us2dbc_mod
!
!svn $Id: u2sdbc_im.F 779 2008-10-04 23:15:47Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for vertically     !
!  integrated Ustokes-velocity.                                              !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: us2dbc, us2dbc_tile
      CONTAINS
!
!***********************************************************************
      SUBROUTINE us2dbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
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
      CALL us2dbc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  OCEAN(ng) % ubar_stokes)
      RETURN
      END SUBROUTINE us2dbc
!
!***********************************************************************
      SUBROUTINE us2dbc_tile (ng, tile,                                 &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       ubar_stokes)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
      USE mod_grid
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      real(r8), intent(inout) :: ubar_stokes(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, j, know
      integer :: Imin, Imax
      real(r8), parameter :: eps = 1.0E-20_r8
      real(r8) :: Ce, Cx, cff5
      real(r8) :: bry_pgr, bry_cor, bry_str, bry_val
      real(r8) :: cff, cff1, cff2, dUde, dUdt, dUdx, tau
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad
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
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!
!  Western edge, implicit upstream radiation condition.
!
        IF (LBC(iwest,isU2Sd,ng)%radiation) THEN
          DO j=Jstr,Jend+1
            grad(Istr  ,j)=ubar_stokes(Istr  ,j  )-                     &
     &                     ubar_stokes(Istr  ,j-1)
            grad(Istr+1,j)=ubar_stokes(Istr+1,j  )-                     &
     &                     ubar_stokes(Istr+1,j-1)
          END DO
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              dUdt=ubar_stokes(Istr+1,j)-ubar_stokes(Istr+1,j)
              dUdx=ubar_stokes(Istr+1,j)-ubar_stokes(Istr+2,j)
              IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
              IF ((dUdt*(grad(Istr+1,j)+grad(Istr+1,j+1))).gt.          &
     &          0.0_r8) THEN
                dUde=grad(Istr+1,j  )
              ELSE
                dUde=grad(Istr+1,j+1)
              END IF
              cff=MAX(dUdx*dUdx+dUde*dUde,eps)
              Cx=dUdt*dUdx
              Ce=0.0_r8
              ubar_stokes(Istr,j)=(cff*ubar_stokes(Istr  ,j)+           &
     &                            Cx *ubar_stokes(Istr+1,j)-            &
     &                            MAX(Ce,0.0_r8)*grad(Istr,j  )-        &
     &                            MIN(Ce,0.0_r8)*grad(Istr,j+1))/       &
     &                           (cff+Cx)
              ubar_stokes(Istr,j)=ubar_stokes(Istr,j)*                  &
     &                            GRID(ng)%umask(Istr,j)
              ubar_stokes(Istr,j)=ubar_stokes(Istr,j)*                  &
     &                            GRID(ng)%umask_wet(Istr,j)
            END IF
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (LBC(iwest,isU2Sd,ng)%clamped) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              ubar_stokes(Istr,j)=BOUNDARY(ng)%ubarstokes_west(j)
              ubar_stokes(Istr,j)=ubar_stokes(Istr,j)*                  &
     &                            GRID(ng)%umask(Istr,j)
              ubar_stokes(Istr,j)=ubar_stokes(Istr,j)*                  &
     &                            GRID(ng)%umask_wet(Istr,j)
            END IF
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (LBC(iwest,isU2Sd,ng)%closed) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              ubar_stokes(Istr,j)=0.0_r8
            END IF
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (LBC(iwest,isU2Sd,ng)%gradient) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              ubar_stokes(Istr,j)=ubar_stokes(Istr+1,j)
              ubar_stokes(Istr,j)=ubar_stokes(Istr,j)*                  &
     &                            GRID(ng)%umask(Istr,j)
              ubar_stokes(Istr,j)=ubar_stokes(Istr,j)*                  &
     &                            GRID(ng)%umask_wet(Istr,j)
            END IF
          END DO
!
!  Western edge, nested.
!
        ELSE IF (LBC(iwest,isU2Sd,ng)%nested) THEN
          DO j=Jstr,Jend
              ubar_stokes(Istr,j)=ubar_stokes(Istr+1,j)
              ubar_stokes(Istr,j)=ubar_stokes(Istr,j)*                  &
     &                            GRID(ng)%umask(Istr,j)
              ubar_stokes(Istr,j)=ubar_stokes(Istr,j)*                  &
     &                            GRID(ng)%umask_wet(Istr,j)
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
!
!  Eastern edge, implicit upstream radiation condition.
!
        IF (LBC(ieast,isU2Sd,ng)%radiation) THEN
          DO j=Jstr,Jend+1
            grad(Iend  ,j)=ubar_stokes(Iend  ,j  )-                     &
     &                     ubar_stokes(Iend  ,j-1)
            grad(Iend+1,j)=ubar_stokes(Iend+1,j  )-                     &
     &                     ubar_stokes(Iend+1,j-1)
          END DO
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              dUdt=ubar_stokes(Iend,j)-ubar_stokes(Iend  ,j)
              dUdx=ubar_stokes(Iend,j)-ubar_stokes(Iend-1,j)
              IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
              IF ((dUdt*(grad(Iend,j)+grad(Iend,j+1))).gt.0.0_r8) THEN
                dUde=grad(Iend,j)
              ELSE
                dUde=grad(Iend,j+1)
              END IF
              cff=MAX(dUdx*dUdx+dUde*dUde,eps)
              Cx=dUdt*dUdx
              Ce=0.0_r8
              ubar_stokes(Iend+1,j)=(cff*ubar_stokes(Iend+1,j)+         &
     &                              Cx *ubar_stokes(Iend  ,j)-          &
     &                              MAX(Ce,0.0_r8)*grad(Iend+1,j  )-    &
     &                              MIN(Ce,0.0_r8)*grad(Iend+1,j+1))/   &
     &                             (cff+Cx)
              ubar_stokes(Iend+1,j)=ubar_stokes(Iend+1,j)*              &
     &                              GRID(ng)%umask(Iend+1,j)
              ubar_stokes(Iend+1,j)=ubar_stokes(Iend+1,j)*              &
     &                              GRID(ng)%umask_wet(Iend+1,j)
            END IF
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (LBC(ieast,isU2Sd,ng)%clamped) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              ubar_stokes(Iend+1,j)=BOUNDARY(ng)%ubarstokes_east(j)
              ubar_stokes(Iend+1,j)=ubar_stokes(Iend+1,j)*              &
     &                              GRID(ng)%umask(Iend+1,j)
              ubar_stokes(Iend+1,j)=ubar_stokes(Iend+1,j)*              &
     &                              GRID(ng)%umask_wet(Iend+1,j)
            END IF
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (LBC(ieast,isU2Sd,ng)%closed) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              ubar_stokes(Iend+1,j)=0.0_r8
            END IF
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (LBC(ieast,isU2Sd,ng)%gradient) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              ubar_stokes(Iend+1,j)=ubar_stokes(Iend,j)
              ubar_stokes(Iend+1,j)=ubar_stokes(Iend+1,j)*              &
     &                              GRID(ng)%umask(Iend+1,j)
              ubar_stokes(Iend+1,j)=ubar_stokes(Iend+1,j)*              &
     &                              GRID(ng)%umask_wet(Iend+1,j)
            END IF
          END DO
!
!  Eastern edge, nested.
!
        ELSE IF (LBC(ieast,isU2Sd,ng)%nested) THEN
          DO j=Jstr,Jend
              ubar_stokes(Iend+1,j)=ubar_stokes(Iend,j)
              ubar_stokes(Iend+1,j)=ubar_stokes(Iend+1,j)*              &
     &                              GRID(ng)%umask(Iend+1,j)
              ubar_stokes(Iend+1,j)=ubar_stokes(Iend+1,j)*              &
     &                              GRID(ng)%umask_wet(Iend+1,j)
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
!
!  Southern edge, implicit upstream radiation condition.
!
        IF (LBC(isouth,isU2Sd,ng)%radiation) THEN
          DO i=IstrU-1,Iend
            grad(i,Jstr-1)=ubar_stokes(i+1,Jstr-1)-                     &
     &                     ubar_stokes(i  ,Jstr-1)
            grad(i,Jstr  )=ubar_stokes(i+1,Jstr  )-                     &
     &                     ubar_stokes(i  ,Jstr  )
          END DO
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              dUdt=ubar_stokes(i,Jstr)-ubar_stokes(i,Jstr  )
              dUde=ubar_stokes(i,Jstr)-ubar_stokes(i,Jstr+1)
              IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
              IF ((dUdt*(grad(i-1,Jstr)+grad(i,Jstr))).gt.0.0_r8) THEN
                dUdx=grad(i-1,Jstr)
              ELSE
                dUdx=grad(i  ,Jstr)
              END IF
              cff=MAX(dUdx*dUdx+dUde*dUde,eps)
              Cx=0.0_r8
              Ce=dUdt*dUde
              ubar_stokes(i,Jstr-1)=(cff*ubar_stokes(i,Jstr-1)+         &
     &                              Ce*ubar_stokes(i,Jstr  )-           &
     &                              MAX(Cx,0.0_r8)*grad(i-1,Jstr-1)-    &
     &                              MIN(Cx,0.0_r8)*grad(i  ,Jstr-1))/   &
     &                              (cff+Ce)
              ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr-1)*              &
     &                              GRID(ng)%umask(i,Jstr-1)
              ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr-1)*              &
     &                              GRID(ng)%umask_wet(i,Jstr-1)
            END IF
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (LBC(isouth,isU2Sd,ng)%clamped) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              ubar_stokes(i,Jstr-1)=BOUNDARY(ng)%ubarstokes_south(i)
              ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr-1)*              &
     &                              GRID(ng)%umask(i,Jstr-1)
              ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr-1)*              &
     &                              GRID(ng)%umask_wet(i,Jstr-1)
            END IF
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (LBC(isouth,isU2Sd,ng)%gradient) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr)
              ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr-1)*              &
     &                              GRID(ng)%umask(i,Jstr-1)
            END IF
          END DO
!
!  Southern edge, nested.
!
        ELSE IF (LBC(isouth,isU2Sd,ng)%nested) THEN
          DO i=IstrU,Iend
              ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr)
              ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr-1)*              &
     &                              GRID(ng)%umask(i,Jstr-1)
          END DO
!
!  Southern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
        ELSE IF (LBC(isouth,isU2Sd,ng)%closed) THEN
          IF (EWperiodic(ng)) THEN
            Imin=IstrU
            Imax=Iend
          ELSE
            Imin=Istr
            Imax=IendR
          END IF
          DO i=Imin,Imax
            IF (LBC_apply(ng)%south(i)) THEN
              ubar_stokes(i,Jstr-1)=gamma2(ng)*ubar_stokes(i,Jstr)
              ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr-1)*              &
     &                              GRID(ng)%umask(i,Jstr-1)
              ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr-1)*              &
     &                              GRID(ng)%umask_wet(i,Jstr-1)
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
!
!  Northern edge, implicit upstream radiation condition.
!
        IF (LBC(inorth,isU2Sd,ng)%radiation) THEN
          DO i=IstrU-1,Iend
            grad(i,Jend  )=ubar_stokes(i+1,Jend  )-                     &
     &                     ubar_stokes(i  ,Jend  )
            grad(i,Jend+1)=ubar_stokes(i+1,Jend+1)-                     &
     &                     ubar_stokes(i  ,Jend+1)
          END DO
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              dUdt=ubar_stokes(i,Jend)-ubar_stokes(i,Jend  )
              dUde=ubar_stokes(i,Jend)-ubar_stokes(i,Jend-1)
              IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
              IF ((dUdt*(grad(i-1,Jend)+grad(i,Jend))).gt.0.0_r8) THEN
                dUdx=grad(i-1,Jend)
              ELSE
                dUdx=grad(i  ,Jend)
              END IF
              cff=MAX(dUdx*dUdx+dUde*dUde,eps)
              Cx=0.0_r8
              Ce=dUdt*dUde
              ubar_stokes(i,Jend+1)=(cff*ubar_stokes(i,Jend+1)+         &
     &                              Ce *ubar_stokes(i,Jend  )-          &
     &                              MAX(Cx,0.0_r8)*grad(i-1,Jend+1)-    &
     &                              MIN(Cx,0.0_r8)*grad(i  ,Jend+1))/   &
     &                              (cff+Ce)
              ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend+1)*              &
     &                              GRID(ng)%umask(i,Jend+1)
              ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend+1)*              &
     &                              GRID(ng)%umask_wet(i,Jend+1)
            END IF
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (LBC(inorth,isU2Sd,ng)%clamped) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              ubar_stokes(i,Jend+1)=BOUNDARY(ng)%ubarstokes_north(i)
              ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend+1)*              &
     &                              GRID(ng)%umask(i,Jend+1)
              ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend+1)*              &
     &                              GRID(ng)%umask_wet(i,Jend+1)
            END IF
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (LBC(inorth,isU2Sd,ng)%gradient) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend)
              ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend+1)*              &
     &                              GRID(ng)%umask(i,Jend+1)
            END IF
          END DO
!
!  Northern edge, nested.
!
        ELSE IF (LBC(inorth,isU2Sd,ng)%nested) THEN
          DO i=IstrU,Iend
              ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend)
              ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend+1)*              &
     &                              GRID(ng)%umask(i,Jend+1)
          END DO
!
!  Northern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
        ELSE IF (LBC(inorth,isU2Sd,ng)%closed) THEN
          IF (EWperiodic(ng)) THEN
            Imin=IstrU
            Imax=Iend
          ELSE
            Imin=Istr
            Imax=IendR
          END IF
          DO i=Imin,Imax
            IF (LBC_apply(ng)%north(i)) THEN
              ubar_stokes(i,Jend+1)=gamma2(ng)*ubar_stokes(i,Jend)
              ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend+1)*              &
     &                              GRID(ng)%umask(i,Jend+1)
              ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend+1)*              &
     &                              GRID(ng)%umask_wet(i,Jend+1)
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF ((LBC_apply(ng)%south(Istr  ).and.                         &
     &        LBC_apply(ng)%west (Jstr-1)).or.                          &
     &        (LBC(iwest,isU2Sd,ng)%nested.and.                         &
     &         LBC(isouth,isU2Sd,ng)%nested)) THEN
              ubar_stokes(Istr,Jstr-1)=0.5_r8*                          &
     &                                 (ubar_stokes(Istr+1,Jstr-1)+     &
     &                                  ubar_stokes(Istr  ,Jstr  ))
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF ((LBC_apply(ng)%south(Iend+1).and.                         &
     &        LBC_apply(ng)%east (Jstr-1)).or.                          &
     &        (LBC(ieast,isU2Sd,ng)%nested.and.                         &
     &         LBC(isouth,isU2Sd,ng)%nested)) THEN
              ubar_stokes(Iend+1,Jstr-1)=0.5_r8*                        &
     &                                   (ubar_stokes(Iend  ,Jstr-1)+   &
     &                                    ubar_stokes(Iend+1,Jstr  ))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF ((LBC_apply(ng)%north(Istr  ).and.                         &
     &        LBC_apply(ng)%west (Jend+1)).or.                          &
     &        (LBC(iwest,isU2Sd,ng)%nested.and.                         &
     &         LBC(inorth,isU2Sd,ng)%nested)) THEN
              ubar_stokes(Istr,Jend+1)=0.5_r8*                          &
     &                                 (ubar_stokes(Istr  ,Jend  )+     &
     &                                  ubar_stokes(Istr+1,Jend+1))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF ((LBC_apply(ng)%north(Iend+1).and.                         &
     &        LBC_apply(ng)%east (Jend+1)).or.                          &
     &        (LBC(ieast,isU2Sd,ng)%nested.and.                         &
     &         LBC(inorth,isU2Sd,ng)%nested)) THEN
              ubar_stokes(Iend+1,Jend+1)=0.5_r8*                        &
     &                                   (ubar_stokes(Iend+1,Jend  )+   &
     &                                    ubar_stokes(Iend  ,Jend+1))
          END IF
        END IF
      END IF
      RETURN
      END SUBROUTINE us2dbc_tile
      END MODULE us2dbc_mod
