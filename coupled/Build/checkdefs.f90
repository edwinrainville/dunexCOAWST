      SUBROUTINE checkdefs
!
!svn $Id: checkdefs.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine checks activated C-preprocessing options for        !
!  consistency.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE mod_strings
!
      USE strings_mod, ONLY : uppercase
!
      implicit none
!
!  Local variable declarations.
!
      integer :: iatms = 0
      integer :: ibbl = 0
      integer :: ibiology = 0
      integer :: idriver = 0
      integer :: itrcHadv = 0
      integer :: itrcVadv = 0
      integer :: itrcHadvtl = 0
      integer :: itrcVadvtl = 0
      integer :: ivelHadv = 0
      integer :: ivelVadv = 0
      integer :: ivmix = 0
      integer :: nearshore = 0
      integer :: is, lstr, ng
!
!-----------------------------------------------------------------------
!  Report activated C-preprocessing options.
!-----------------------------------------------------------------------
!
      Coptions=' '
      IF (Master) WRITE (stdout,10)
  10  FORMAT (/,' Activated C-preprocessing Options:',/)
  20  FORMAT (1x,a,t27,a)
!
      IF (Master) THEN
        WRITE (stdout,20) TRIM(ADJUSTL(MyAppCPP)), TRIM(ADJUSTL(title))
      END IF
      is=LEN_TRIM(Coptions)+1
      lstr=LEN_TRIM(MyAppCPP)
      Coptions(is:is+lstr)=TRIM(ADJUSTL(MyAppCPP))
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is)=','
!
      IF (Master) WRITE (stdout,20) 'ANA_BSFLUX',                       &
     &   'Analytical kinematic bottom salinity flux'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_BSFLUX,'
!
      IF (Master) WRITE (stdout,20) 'ANA_BTFLUX',                       &
     &   'Analytical kinematic bottom temperature flux'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_BTFLUX,'
!
      IF (Master) WRITE (stdout,20) 'ANA_FSOBC',                        &
     &   'Analytical free-surface boundary conditions'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ANA_FSOBC,'
!
      IF (Master) WRITE (stdout,20) 'ANA_INITIAL',                      &
     &   'Analytical initial conditions'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' ANA_INITIAL,'
!
      IF (Master) WRITE (stdout,20) 'ANA_M2OBC',                        &
     &   'Analytical 2D momentum boundary conditions'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ANA_M2OBC,'
!
      IF (Master) WRITE (stdout,20) 'ANA_SMFLUX',                       &
     &   'Analytical kinematic surface momentum flux'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_SMFLUX,'
!
      IF (Master) WRITE (stdout,20) 'ANA_SSFLUX',                       &
     &   'Analytical kinematic surface salinity flux'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_SSFLUX,'
!
      IF (Master) WRITE (stdout,20) 'ANA_STFLUX',                       &
     &   'Analytical kinematic surface temperature flux'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_STFLUX,'
!
      IF (Master) WRITE (stdout,20) 'ASSUMED_SHAPE',                    &
     &   'Using assumed-shape arrays'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' ASSUMED_SHAPE,'
      IF (Master) WRITE (stdout,20) '!BOUNDARY_ALLGATHER',              &
     &   'Using mpi_allreduce in mp_boundary routine'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+21)=' !BOUNDARY_ALLGATHER,'
!
      IF (Master) WRITE (stdout,20) '!COLLECT_ALL...',                  &
     &   'Using mpi_isend/mpi_recv in mp_collect routine'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' !COLLECT_ALL...,'
!
      IF (Master) WRITE (stdout,20) 'DJ_GRADPS',                        &
     &   'Parabolic Splines density Jacobian (Shchepetkin, 2002)'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' DJ_GRADPS,'
!
      IF (Master) WRITE (stdout,20) 'DOUBLE_PRECISION',                 &
     &   'Double precision arithmetic numerical kernel.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' DOUBLE_PRECISION,'
!
      IF (Master) WRITE (stdout,20) 'GLS_MIXING',                       &
     &   'Generic Length-Scale turbulence closure'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' GLS_MIXING,'
      ivmix=ivmix+1
!
      IF (Master) WRITE (stdout,20) 'KANTHA_CLAYSON',                   &
     &   'Kantha and Clayson stability function formulation'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' KANTHA_CLAYSON,'
!
      IF (Master) WRITE (stdout,20) 'LIMIT_VDIFF',                      &
     &   'Impose an upper limit on vertical diffusion coefficient'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' LIMIT_VDIFF,'
!
      IF (Master) WRITE (stdout,20) 'MASKING',                          &
     &   'Land/Sea masking'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' MASKING,'
      IF (Master) WRITE (stdout,20) 'MCT_LIB',                          &
     &   'Using Model Coupling Toolkit library'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' MCT_LIB,'
!
      IF (Master) WRITE (stdout,20) 'MIX_S_TS',                         &
     &   'Mixing of tracers along constant S-surfaces'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' MIX_S_TS,'
!
      IF (Master) WRITE (stdout,20) 'MIX_S_UV',                         &
     &   'Mixing of momentum along constant S-surfaces'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' MIX_S_UV,'
!
      IF (Master) WRITE (stdout,20) 'MPI',                              &
     &   'MPI distributed-memory configuration'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+4)=' MPI,'
!
      IF (Master) WRITE (stdout,20) 'WEC_VF',                           &
     &   'Vortex Force wave current interaction'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' WEC_VF,'
      nearshore=nearshore+1
      IF (Master) WRITE (stdout,20) 'WDISS_WAVEMOD',                    &
     &   'Wave energy dissipation acquired from coupled wave model'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' WDISS_WAVEMOD,'
      IF (Master) WRITE (stdout,20) 'ROLLER_RENIERS',                   &
     &   'Wave energy roller based on Reniers 2004'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' ROLLER_RENIERS,'
!
      IF (Master) WRITE (stdout,20) 'NONLINEAR',                        &
     &   'Nonlinear Model'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' NONLINEAR,'
!
      IF (Master) WRITE (stdout,20) '!NONLIN_EOS',                      &
     &   'Linear Equation of State for seawater'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' !NONLIN_EOS,'
!
      IF (Master) WRITE (stdout,20) 'N2S2_HORAVG',                      &
     &   'Horizontal smoothing of buoyancy and shear'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' N2S2_HORAVG,'
!
      IF (Master) WRITE (stdout,20) 'POWER_LAW',                        &
     &   'Power-law shape time-averaging barotropic filter'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' POWER_LAW,'
!
      IF (Master) WRITE (stdout,20) 'PROFILE',                          &
     &   'Time profiling activated'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' PROFILE,'
!
      IF (Master) WRITE (stdout,20) 'K_GSCHEME',                        &
     &   'Third-order upstream advection of TKE fields'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' K_GSCHEME,'
!
      IF (Master) WRITE (stdout,20) 'REDUCE_ALLREDUCE',                 &
     &   'Using mpi_allreduce in mp_reduce routine'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' REDUCE_ALLREDUCE,'
!
      IF (Master) WRITE (stdout,20) 'RI_SPLINES',                       &
     &   'Parabolic Spline Reconstruction for Richardson Number'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' RI_SPLINES,'
!
      IF (Master) WRITE (stdout,20) '!RST_SINGLE',                      &
     &   'Double precision fields in restart NetCDF file'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' !RST_SINGLE,'
!
      IF (Master) WRITE (stdout,20) 'SOLVE3D',                          &
     &   'Solving 3D Primitive Equations'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' SOLVE3D,'
!
       IF (Master) WRITE (stdout,20) 'SSW_BBL',                         &
     &   'Styles and Glenn Bottom Boundary Layer - modified'
       is=LEN_TRIM(Coptions)+1
       Coptions(is:is+8)=' SSW_BBL,'
       ibbl=ibbl+1
!
       IF (Master) WRITE (stdout,20) 'SSW_CALC_ZNOT',                   &
     &   'Internal computation of bottom roughness'
       is=LEN_TRIM(Coptions)+1
       Coptions(is:is+14)=' SSW_CALC_ZNOT,'
!
       IF (Master) WRITE (stdout,20) 'SSW_LOGINT',                      &
     &   'Bottom currents logarithmic interpolation'
       is=LEN_TRIM(Coptions)+1
       Coptions(is:is+11)=' SSW_LOGINT,'
!
      IF (Master) WRITE (stdout,20) 'SWAN_COUPLING',                    &
     &   'SWAN model coupling'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' SWAN_COUPLING,'
!
      IF (Master) WRITE (stdout,20) 'TKE_WAVEDISS',                     &
     &   'Wave breaking surface tke flux based on wave amplitude'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' TKE_WAVEDISS,'
!
      IF (Master) WRITE (stdout,20) 'TS_DIF2',                          &
     &   'Harmonic mixing of tracers'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' TS_DIF2,'
!
      IF (Master) WRITE (stdout,20) 'UV_ADV',                           &
     &   'Advection of momentum'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' UV_ADV,'
!
      IF (Master) WRITE (stdout,20) 'UV_U3HADVECTION',                  &
     &   'Third-order upstream horizontal advection of 3D momentum'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' UV_U3HADVECTION,'
      ivelHadv=ivelHadv+1
!
      IF (Master) WRITE (stdout,20) 'UV_C4VADVECTION',                  &
     &   'Fourth-order centered vertical advection of momentum'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' UV_C4VADVECTION,'
      ivelVadv=ivelVadv+1
      IF (Master) WRITE (stdout,20) 'UV_KIRBY',                         &
     &   'Compute uwave and vwave Kirby avg velocities.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' UV_KIRBY,'
!
      IF (Master) WRITE (stdout,20) 'UV_VIS2',                          &
     &   'Harmonic mixing of momentum'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' UV_VIS2,'
!
      IF (Master) WRITE (stdout,20) 'VAR_RHO_2D',                       &
     &   'Variable density barotropic mode'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' VAR_RHO_2D,'
!
      IF (Master) WRITE (stdout,20) 'WAVES_OCEAN',                      &
     &   'Two-way wave-ocean models coupling.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' WAVES_OCEAN,'
!
      IF (Master) WRITE (stdout,20) 'WET_DRY',                          &
     &   'Wetting and drying activated'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' WET_DRY,'
!
      IF (Master) WRITE (stdout,20) 'ZOS_HSIG',                         &
     &   'Wave amplitude used for Zos calculation'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ZOS_HSIG,'
!
!-----------------------------------------------------------------------
!  Stop if unsupported C-preprocessing options or report issues with
!  particular options.
!-----------------------------------------------------------------------
!
      CALL checkadj
!
!-----------------------------------------------------------------------
!  Check C-preprocessing options.
!-----------------------------------------------------------------------
!
!  Stop if more than one vertical closure scheme is selected.
!
      IF (Master.and.(ivmix.gt.1)) THEN
        WRITE (stdout,30)
  30    FORMAT (/,' CHECKDEFS - only one vertical closure scheme',      &
     &            ' is allowed.')
        exit_flag=5
      END IF
!
!  Stop if more that one bottom stress formulation is selected.
!
      IF (Master.and.(ibbl.gt.1)) THEN
        WRITE (stdout,40)
  40    FORMAT (/,' CHECKDEFS - only one bottom stress formulation is', &
     &            ' allowed.')
        exit_flag=5
      END IF
!
!  Stop if no bottom stress formulation is selected.
!
      IF (Master.and.(ibbl.eq.0)) THEN
        WRITE (stdout,50)
  50    FORMAT (/,' CHECKDEFS - no bottom stress formulation is',       &
     &            ' selected.')
        exit_flag=5
      END IF
!
!  Stop if more than one biological module is selected.
!
      IF (Master.and.(ibiology.gt.1)) THEN
        WRITE (stdout,60)
  60    FORMAT (/,' CHECKDEFS - only one biology MODULE is allowed.')
        exit_flag=5
      END IF
!
!  Stop if more that one model driver is selected.
!
      IF (Master.and.(idriver.gt.1)) THEN
        WRITE (stdout,70)
  70    FORMAT (/,' CHECKDEFS - only one model example is allowed.')
        exit_flag=5
      END IF
!
!  Stop if more than one advection scheme has been activated.
!
      IF (Master.and.(ivelHadv.gt.1)) THEN
        WRITE (stdout,140) 'horizontal','momentum','ivelHadv =',ivelHadv
        exit_flag=5
      END IF
      IF (Master.and.(ivelVadv.gt.1)) THEN
        WRITE (stdout,140) 'vertical','momentum','ivelVadv =',ivelVadv
        exit_flag=5
      END IF
 140  FORMAT (/,' CHECKDEFS - only one ',a,' advection scheme',         &
     &        /,13x,'is allowed for ',a,', ',a,1x,i1)
!
!  Stop if more that one radiation stress formulation is activated.
!
      IF (Master.and.(nearshore.gt.1)) THEN
        WRITE (stdout,160)
 160    FORMAT (/,' CHECKDEFS - only one wave current formulation'      &
     &            ' is allowed.')
        exit_flag=5
      END IF
      RETURN
      END SUBROUTINE checkdefs
