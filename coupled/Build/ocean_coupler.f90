      MODULE ocean_coupler_mod
!
!svn $Id: ocean_coupler.F 830 2017-01-24 21:21:11Z arango $
!==================================================== John C. Warner ===
!  Copyright (c) 2002-2017 The ROMS/TOMS Group      Hernan G. Arango   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module is used to communicate and exchange data between        !
!  ROMS/TOMS and other coupled model(s)  via the Model Coupling        !
!  Toolkit (MCT), developed at the Argonne National Laboratory.        !
!                                                                      !
!=======================================================================
!
!  Component Model Registry.
!
      USE m_MCTWorld, ONLY : MCTWorld_init => init
      USE m_MCTWorld, ONLY : MCTWorld_clean => clean
!
!  Domain Decomposition Descriptor DataType and associated methods.
!
      USE m_GlobalSegMap, ONLY : GlobalSegMap
      USE m_GlobalSegMap, ONLY : GlobalSegMap_init => init
      USE m_GlobalSegMap, ONLY : GlobalSegMap_lsize => lsize
      USE m_GlobalSegMap, ONLY : GlobalSegMap_clean => clean
      USE m_GlobalSegMap, ONLY : GlobalSegMap_Ordpnts => OrderedPoints
!
!  Field Storage DataType and associated methods.
!
      USE m_AttrVect, ONLY : AttrVect
      USE m_AttrVect, ONLY : AttrVect_init => init
      USE m_AttrVect, ONLY : AttrVect_zero => zero
      USE m_AttrVect, ONLY : AttrVect_lsize => lsize
      USE m_AttrVect, ONLY : AttrVect_clean => clean
      USE m_AttrVect, ONLY : AttrVect_copy => copy
      USE m_AttrVect, ONLY : AttrVect_importRAttr => importRAttr
      USE m_AttrVect, ONLY : AttrVect_exportRAttr => exportRAttr
!
!  Intercomponent communications scheduler.
!
      USE m_Router, ONLY : Router
      USE m_Router, ONLY : Router_init => init
      USE m_Router, ONLY : Router_clean => clean
!
!  Intercomponent transfer.
!
      USE m_Transfer, ONLY: MCT_send => send
      USE m_Transfer, ONLY: MCT_recv => recv
      USE m_Transfer, ONLY: MCT_isend => isend
      USE m_Transfer, ONLY: MCT_irecv => irecv
      USE m_Transfer, ONLY: MCT_waitr => waitrecv
      USE m_Transfer, ONLY: MCT_waits => waitsend
!
      implicit none
!
      PRIVATE
      PUBLIC :: ocean_coupling
      PUBLIC :: initialize_ocn2wav_coupling
      PUBLIC :: initialize_ocn2wav_routers
      PUBLIC :: ocn2wav_coupling
      PUBLIC :: ocnfwav_coupling
      PUBLIC :: finalize_ocn2wav_coupling
!
!  Declarations.
!
      TYPE T_GlobalSegMap_G
        TYPE(GlobalSegMap) :: GSMapROMS       ! GloabalSegMap variables
      END TYPE T_GlobalSegMap_G
      TYPE (T_GlobalSegMap_G), ALLOCATABLE :: GlobalSegMap_G(:)
      TYPE T_AttrVect_G
        TYPE(AttrVect) :: wav2ocn_AV          ! AttrVect variables
        TYPE(AttrVect) :: ocn2wav_AV
      END TYPE T_AttrVect_G
      TYPE (T_AttrVect_G), ALLOCATABLE :: AttrVect_G(:)
      TYPE T_Router_W
        TYPE(Router)   :: ROMStoSWAN          ! Router variables
      END TYPE T_Router_W
      TYPE (T_Router_W), ALLOCATABLE :: Router_W(:,:)
      CONTAINS
      SUBROUTINE initialize_ocn2wav_coupling (ng, tile)
!
!=======================================================================
!                                                                      !
!  Initialize ocean and wave models coupling stream.  This is the      !
!  training phase used to constuct MCT parallel interpolators and      !
!  and stablish communication patterns.                                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mct_coupler_params
      USE mod_kinds
      USE mod_scalars
      USE mod_iounits
!
!  Imported variable definitions.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: Istr, Iend, Jstr, Jend
      integer :: IstrT, IendT, JstrT, JendT
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: Asize, Isize, Jsize, MyError
      integer :: i, ic, iw, j, jc, nprocs
      integer :: nRows, nCols, num_sparse_elems
      integer :: cid, cad
      integer, allocatable :: length(:)
      integer, allocatable :: start(:)
!      integer, dimension(2) :: src_grid_dims, dst_grid_dims
      character (len=70)    :: nc_name
      character (len=20)    :: to_add
      character (len=120)   :: wostring
      character (len=120)   :: owstring
      real(r8) :: cff
!
!-----------------------------------------------------------------------
!  Compute lower and upper bounds over a particular domain partition or
!  tile for RHO-, U-, and V-variables. Notice that "set_bounds.h" is
!  not used here because of implementation of periodicity in other
!  models.
!-----------------------------------------------------------------------
!
      Istr=BOUNDS(ng)%Istr(tile)
      Iend=BOUNDS(ng)%Iend(tile)
      Jstr=BOUNDS(ng)%Jstr(tile)
      Jend=BOUNDS(ng)%Jend(tile)
      IstrT=BOUNDS(ng)%IstrT(tile)
      IendT=BOUNDS(ng)%IendT(tile)
      JstrT=BOUNDS(ng)%JstrT(tile)
      JendT=BOUNDS(ng)%JendT(tile)
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        IstrR=BOUNDS(ng)%Istr(tile)-1
      ELSE
        IstrR=BOUNDS(ng)%Istr(tile)
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        IendR=BOUNDS(ng)%Iend(tile)+1
      ELSE
        IendR=BOUNDS(ng)%Iend(tile)
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        JstrR=BOUNDS(ng)%Jstr(tile)-1
      ELSE
        JstrR=BOUNDS(ng)%Jstr(tile)
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        JendR=BOUNDS(ng)%Jend(tile)+1
      ELSE
        JendR=BOUNDS(ng)%Jend(tile)
      END IF
!
!-----------------------------------------------------------------------
!  Establish MCT communicator.
!-----------------------------------------------------------------------
!
!  Get communicator local rank and size.
!
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
      CALL mpi_comm_size (OCN_COMM_WORLD, nprocs, MyError)
!
      IF (ng.eq.1) THEN
        ALLOCATE(GlobalSegMap_G(Nocn_grids))
        ALLOCATE(AttrVect_G(Nocn_grids))
      END IF
!
!  Initialize MCT coupled model registry.
!
      OCNid=ocnids(ng)
      IF (Nocn_grids.gt.1) THEN
        CALL MCTWorld_init (N_mctmodels, MPI_COMM_WORLD,                &
     &                      OCN_COMM_WORLD,myids=ocnids)
      ELSE
        CALL MCTWorld_init (N_mctmodels, MPI_COMM_WORLD,                &
     &                      OCN_COMM_WORLD,OCNid)
      END IF
!
!  Determine the part of the grid we are working on and develop
!  local segment of the global map.
!
      Jsize=JendR-JstrR+1
      IF (.not.allocated(start)) THEN
        allocate ( start(Jsize) )
      END IF
      IF (.not.allocated(length)) THEN
        allocate ( length(Jsize) )
      END IF
      jc=0
      DO j=JstrR,JendR
        jc=jc+1
        start (jc)=(j)*(Lm(ng)+2)+IstrR+1
        length(jc)=(IendR-IstrR+1)
      END DO
      CALL GlobalSegMap_init (GlobalSegMap_G(ng)%GSMapROMS,             &
     &                        start, length, 0, OCN_COMM_WORLD, OCNid)
!
!  Deallocate working arrays.
!
      IF (allocated(start)) THEN
        deallocate (start)
      END IF
      IF (allocated(length)) THEN
        deallocate (length)
      END IF
!
!
!  Initialize the list of fields from the wave model.
!
      cad=LEN(wostring)
      DO i=1,cad
        wostring(i:i)=''
      END DO
      cid=1
!
      to_add='DISBOT'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':DISSURF'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':DISWCAP'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':HSIGN'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':RTP'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':TMBOT'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':UBOT'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':DIRE'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
      to_add=':DIRN'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':DIREP'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
      to_add=':DIRNP'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':WLEN'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':WLENP'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':QB'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':WDSPR'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':WQP'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
!
!  Finalize and remove trailing spaces from the wostring
!  for the rlist.
!
      cad=LEN_TRIM(wostring)
      wostring=wostring(1:cad)
!
!  Initialize attribute vector holding the export data of
!  the wav model. The Asize is the number of grid point on this
!  processor.
!
      Asize=GlobalSegMap_lsize(GlobalSegMap_G(ng)%GSMapROMS,            &
     &                         OCN_COMM_WORLD)
      CALL AttrVect_init(AttrVect_G(ng)%wav2ocn_AV,                     &
     &                   rList=TRIM(wostring),lsize=Asize)
      CALL AttrVect_zero(AttrVect_G(ng)%wav2ocn_AV)
!
!  Initialize attribute vector that contain the data strings from
!  the ocean model.
!
      cad=LEN(owstring)
      DO i=1,cad
        owstring(i:i)=''
      END DO
      cid=1
!
      to_add='DEPTH'
      cad=LEN_TRIM(to_add)
      write(owstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':WLEV'
      cad=LEN_TRIM(to_add)
      write(owstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':VELX'
      cad=LEN_TRIM(to_add)
      write(owstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':VELY'
      cad=LEN_TRIM(to_add)
      write(owstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':ZO'
      cad=LEN_TRIM(to_add)
      write(owstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
!  Finalize and remove trailing spaces from the owstring
!  for the rlist.
!
      cad=LEN_TRIM(owstring)
      owstring=owstring(1:cad)
!
      CALL AttrVect_init(AttrVect_G(ng)%ocn2wav_AV,                     &
     &                   rList=TRIM(owstring),lsize=Asize)
      CALL AttrVect_zero(AttrVect_G(ng)%ocn2wav_AV)
!
      RETURN
      END SUBROUTINE initialize_ocn2wav_coupling
      SUBROUTINE initialize_ocn2wav_routers (tile)
!
!=======================================================================
!                                                                      !
!  Initialize ocean and wave models coupling stream.  This is the      !
!  training phase used to constuct MCT parallel interpolators and      !
!  and stablish communication patterns.                                !
!                                                                      !
!=======================================================================
!
      USE mod_parallel
      USE mct_coupler_params
!
!  Imported variable definitions.
!
      integer, intent(in) :: tile
!
!  Local variable declarations.
!
      integer :: MyError, nprocs
      integer :: ng, iw
!
!-----------------------------------------------------------------------
!  Establish MCT router.
!-----------------------------------------------------------------------
!
      ALLOCATE(Router_W(Nocn_grids,Nwav_grids))
!
!  Initialize routers to the wave model component.
!
      DO ng=1,Nocn_grids
        DO iw=1,Nwav_grids
          WAVid=wavids(iw)
          CALL Router_init (WAVid, GlobalSegMap_G(ng)%GSMapROMS,        &
     &                      OCN_COMM_WORLD, Router_W(ng,iw)%ROMStoSWAN)
        END DO
      END DO
      RETURN
      END SUBROUTINE initialize_ocn2wav_routers
      SUBROUTINE ocn2wav_coupling (ng, iw, tile)
!
!=======================================================================
!                                                                      !
!  This routine acquires the coupling data streams between waves       !
!  and ocean models.                                                   !
!  coded:                                                              !
!                                                                      !
!  Fields exported to SWAN model:                                      !
!                                                                      !
!     * Bathymetry, bottom elevation (m), [m]                          !
!     * Free-surface, water surface elevation (m), [m]                 !
!     * Depth integrated u-momentum (m/s), [m/s]                       !
!     * Depth integrated v-momentum (m/s), [m/s]                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, iw, tile
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
      CALL wclock_on (ng, iNLM, 48)
      CALL ocn2wav_coupling_tile (ng, iw, tile,                         &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS)
      CALL wclock_off (ng, iNLM, 48)
      RETURN
      END SUBROUTINE ocn2wav_coupling
!
!***********************************************************************
      SUBROUTINE ocn2wav_coupling_tile (ng, iw, tile,                   &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_coupler
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
      USE mod_iounits
      USE mod_sedbed
      USE mod_sediment
      USE mod_coupling
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
      USE exchange_2d_mod, ONLY : exchange_u2d_tile
      USE exchange_2d_mod, ONLY : exchange_v2d_tile
      USE mp_exchange_mod, ONLY : mp_exchange2d
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, iw, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: Asize, MyError, Tag
      integer :: gtype, i, id, ifield, ij, j, k, status
      real(r8), parameter ::  Lwave_min = 1.0_r8
      real(r8), parameter ::  Lwave_max = 500.0_r8
      real(r8) :: add_offset, scale
      real(r8) :: cff, ramp
      real(r8) :: cff1, cff2, cff3, cff4, kwn, prof, u_cff, v_cff
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ubar_rho
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vbar_rho
      real(r8), pointer :: A(:)
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
!  Modify ranges to allow full exchange of fields for periodic applications.
!
      IF (EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          IstrR=Istr-1
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          IendR=Iend+1
        END IF
      END IF
      IF (NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          JstrR=Jstr-1
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          JendR=Jend+1
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Allocate communications array.
!-----------------------------------------------------------------------
!
      Asize=GlobalSegMap_lsize (GlobalSegMap_G(ng)%GSMapROMS,           &
     &      OCN_COMM_WORLD)
      allocate ( A(Asize) )
      A=0.0_r8
!
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
!
!-----------------------------------------------------------------------
!  Export fields from ocean (ROMS) to wave (SWAN) model.
!-----------------------------------------------------------------------
!
!  Depth (bathymetry).
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          A(ij)=GRID(ng)%h(i,j)
          IF (GRID(ng)%rmask_wet(i,j).lt.0.5_r8) THEN
            A(ij)=9999.0_r8
          END IF
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(ng)%ocn2wav_AV, "DEPTH",    &
     &                           A, Asize)
!
!  Water level (free-surface).
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          A(ij)=OCEAN(ng)%zeta(i,j,knew(ng))
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(ng)%ocn2wav_AV, "WLEV",     &
     &                           A, Asize)
!
!  U-velocity at RHO-points.
!
      DO j=JstrR,JendR
        DO i=Istr,Iend
!
! Compute the coupling current according to Kirby and Chen (1989).
!
          kwn=2.0_r8*pi/FORCES(ng)%Lwave(i,j)
          prof=GRID(ng)%h(i,j)+COUPLING(ng)%Zt_avg1(i,j)
          cff1=0.0_r8
          cff2=2.0_r8*kwn*prof
          IF (cff2.lt.700.0_r8) THEN
            cff2=2.0_r8*kwn
          ELSE
            cff2=700.0_r8/prof
          ENDIF
          cff3=0.0_r8
          DO k=1,N(ng)
            u_cff=0.5_r8*(OCEAN(ng)%u(i,  j,k,nrhs(ng))+                    &
     &                    OCEAN(ng)%u(i+1,j,k,nrhs(ng)))
            cff4=cosh(cff2*(GRID(ng)%h(i,j)+GRID(ng)%z_r(i,j,k)))*      &
     &           GRID(ng)%Hz(i,j,k)
            cff1=cff1+cff4*u_cff
            cff3=cff3+cff4
          END DO
          ubar_rho(i,j)=cff1/cff3
        END DO
      END DO
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=JstrR,JendR
          ubar_rho(IstrR,j)=ubar_rho(IstrR+1,j)
        END DO
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=JstrR,JendR
          ubar_rho(IendR,j)=ubar_rho(IendR-1,j)
        END DO
      END IF
        DO j=JstrR,JendR
          DO i=IstrR,IendR
             OCEAN(ng)%uwave(i,j)=ubar_rho(i,j)
          ENDDO
        ENDDO
!
!  V-velocity at RHO-points.
!
      DO j=Jstr,Jend
        DO i=IstrR,IendR
!
! Compute the coupling current according to Kirby and Chen (1989).
!
          kwn=2.0_r8*pi/FORCES(ng)%Lwave(i,j)
          prof=GRID(ng)%h(i,j)+COUPLING(ng)%Zt_avg1(i,j)
          cff1=0.0_r8
          cff2=2.0_r8*kwn*prof
          IF (cff2.lt.700.0_r8) THEN
            cff2=2.0_r8*kwn
          ELSE
            cff2=700.0_r8/prof
          ENDIF
          cff3=0.0_r8
          DO k=1,N(ng)
             v_cff=0.5_r8*(OCEAN(ng)%v(i,  j,k,nrhs(ng))+                   &
     &                     OCEAN(ng)%v(i,j+1,k,nrhs(ng)))
             cff4=cosh(cff2*(GRID(ng)%h(i,j)+GRID(ng)%z_r(i,j,k)))*     &
     &            GRID(ng)%Hz(i,j,k)
             cff1=cff1+cff4*v_cff
             cff3=cff3+cff4
          END DO
          vbar_rho(i,j)=cff1/cff3
        END DO
      END DO
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=IstrR,IendR
          vbar_rho(i,JendR)=vbar_rho(i,JendR-1)
        END DO
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=IstrR,IendR
          vbar_rho(i,JstrR)=vbar_rho(i,JstrR+1)
        END DO
      END IF
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          OCEAN(ng)%vwave(i,j)=vbar_rho(i,j)
        ENDDO
      ENDDO
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff1=ubar_rho(i,j)
          A(ij)=cff1
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(ng)%ocn2wav_AV, "VELX",     &
     &                           A, Asize)
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff1=vbar_rho(i,j)
          A(ij)=cff1
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(ng)%ocn2wav_AV, "VELY",     &
     &                           A, Asize)
!
!  bottom roughness.
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
!         Specify this to be Madsen 0.05 minimum.
          A(ij)=MAX(0.05_r8, SEDBED(ng)%bottom(i,j,izNik)*30.0_r8)
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(ng)%ocn2wav_AV, "ZO",       &
     &                           A, Asize)
!
!  Send ocean fields to wave model.
!
      Tag=ng*100+0*10+iw
      CALL MCT_isend (AttrVect_G(ng)%ocn2wav_AV,                        &
     &                Router_W(ng,iw)%ROMStoSWAN, Tag)
      CALL MCT_waits (Router_W(ng,iw)%ROMStoSWAN)
      IF (MyError.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,20) 'wave model, MyError = ', MyError
        END IF
        exit_flag=2
        RETURN
      ELSE
        IF (Master) THEN
        WRITE (stdout,36) ' ** ROMS grid ',ng,                          &
     &                    ' sent data to SWAN grid ',iw
 36     FORMAT (a14,i2,a24,i2)
        END IF
      END IF
!
!  Deallocate communication arrays.
!
      deallocate (A)
!
 10   FORMAT (' OCN2WAV_COUPLING - error while receiving fields from ', &
     &        a, i4)
 20   FORMAT (' OCN2WAV_COUPLING - error while sending fields to: ',    &
     &        a, i4)
      RETURN
      END SUBROUTINE ocn2wav_coupling_tile
      SUBROUTINE ocnfwav_coupling (ng, iw, tile)
!
!=======================================================================
!                                                                      !
!  This routine acquires the coupling data streams between waves       !
!  and ocean models.                                                   !
!                                                                      !
!  Fields imported from SWAN model:                                    !
!                                                                      !
!     * Mean wave direction (degrees), [radians]                       !
!     * Peak wave direction (degrees), [radians]                       !
!     * Significant wave height (m), [m]                               !
!     * Average wave length (m), [m]                                   !
!     * Surface wave relative peak period (s), [s]                     !
!     * Bottom wave period (s), [s]                                    !
!     * Percent of breakig waves (nondimensional), [nondimensional]    !
!     * Wave energy dissipation (W/m2), [m3/s3]                        !
!     * Wave bottom orbital velocity (m/s), [m/s]                      !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, iw, tile
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
      CALL wclock_on (ng, iNLM, 48)
      CALL ocnfwav_coupling_tile (ng, iw, tile,                         &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS)
      CALL wclock_off (ng, iNLM, 48)
      RETURN
      END SUBROUTINE ocnfwav_coupling
!
!***********************************************************************
      SUBROUTINE ocnfwav_coupling_tile (ng, iw, tile,                       &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mct_coupler_params
      USE mod_param
      USE mod_parallel
      USE mod_coupler
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
      USE mod_iounits
      USE mod_sedbed
      USE mod_sediment
      USE mod_coupling
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
      USE exchange_2d_mod, ONLY : exchange_u2d_tile
      USE exchange_2d_mod, ONLY : exchange_v2d_tile
      USE distribute_mod,  ONLY : mp_reduce
      USE mp_exchange_mod, ONLY : mp_exchange2d
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, iw, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: Asize, MyError, Tag
      integer :: gtype, i, id, ifield, ij, j, k, status
      real(r8), parameter ::  Lwave_min = 1.0_r8
      real(r8), parameter ::  Lwave_max = 500.0_r8
!     real(r8), parameter ::  Large = 1.0E+20_r8
      real(r8) :: add_offset, scale
      real(r8) :: cff, fac, ramp
      real(r8) :: cff1, cff2, cff3, cff4, kwn, prof, u_cff, v_cff
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ubar_rho
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vbar_rho
      real(r8), dimension(2) :: range
      real(r8), pointer :: A(:)
      real(r8), pointer :: A1(:)
      character (len=3), dimension(2) :: op_handle
!
!-----------------------------------------------------------------------
!  Compute lower and upper bounds over a particular domain partition or
!  tile for RHO-, U-, and V-variables. Notice that "set_bounds.h" is
!  not used here because of implementation of periodicity in other
!  models.
!-----------------------------------------------------------------------
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
      op_handle(1)='MIN'
      op_handle(2)='MAX'
!
!  Modify ranges to allow full exchange of fields for periodic applications.
!
      IF (EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          IstrR=Istr-1
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          IendR=Iend+1
        END IF
      END IF
      IF (NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          JstrR=Jstr-1
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          JendR=Jend+1
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Allocate communications array.
!-----------------------------------------------------------------------
!
      Asize=GlobalSegMap_lsize (GlobalSegMap_G(ng)%GSMapROMS,           &
     &      OCN_COMM_WORLD)
      allocate ( A(Asize) )
      allocate ( A1(Asize) )
      A=0.0_r8
      A1=0.0_r8
!
!-----------------------------------------------------------------------
!  Import fields from wave model (SWAN) to ocean model (ROMS).
!-----------------------------------------------------------------------
!
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
      Tag=ng*100+0*10+iw
      CALL MCT_irecv (AttrVect_G(ng)%wav2ocn_AV,                        &
     &                Router_W(ng,iw)%ROMStoSWAN, Tag)
!     Wait to make sure the SWAN data has arrived.
      CALL MCT_waitr (AttrVect_G(ng)%wav2ocn_AV,                        &
     &                Router_W(ng,iw)%ROMStoSWAN)
!
      IF (MyError.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,10) 'wave model, MyError = ', MyError
        END IF
        exit_flag=2
        RETURN
      ELSE
        IF (Master) THEN
        WRITE (stdout,36) ' ** ROMS grid ',ng,                          &
     &                    ' recv data from SWAN grid ',iw
 36     FORMAT (a14,i2,a26,i2)
        END IF
      END IF
!
!  Set ramp coefficient.
!
      ramp=1.0_r8
!
!  Receive fields from wave model.
 40         FORMAT (a36,1x,2(1pe14.6))
!
!  Wave dissipation due to bottom friction.
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "DISBOT",   &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      fac=1.0_r8/rho0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=MAX(0.0_r8,A(ij)*ramp)*fac
          IF (iw.eq.1) THEN
            FORCES(ng)%Dissip_fric(i,j)=cff
          ELSE
            FORCES(ng)%Dissip_fric(i,j)=FORCES(ng)%Dissip_fric(i,j)+    &
     &                                  cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'SWANtoROMS Min/Max DISBOT  (Wm-2):  ',        &
     &                    range(1),range(2)
      END IF
!
!  Wave dissipation due to surface breaking.
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "DISSURF",  &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      fac=1.0_r8/rho0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=MAX(0.0_r8,A(ij)*ramp)*fac
          IF (iw.eq.1) THEN
            FORCES(ng)%Dissip_break(i,j)=cff
          ELSE
            FORCES(ng)%Dissip_break(i,j)=FORCES(ng)%Dissip_break(i,j)+  &
     &                                   cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'SWANtoROMS Min/Max DISSURF (Wm-2):  ',        &
     &                    range(1),range(2)
      END IF
!
!  Wave dissipation due to white capping.
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "DISWCAP",  &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      fac=1.0_r8/rho0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=MAX(0.0_r8,A(ij)*ramp)*fac
          IF (iw.eq.1) THEN
            FORCES(ng)%Dissip_wcap(i,j)=cff
          ELSE
            FORCES(ng)%Dissip_wcap(i,j)=FORCES(ng)%Dissip_wcap(i,j)+    &
     &                                  cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'SWANtoROMS Min/Max DISWCAP (Wm-2):  ',        &
     &                    range(1),range(2)
      END IF
!
!  Wave height.
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "HSIGN",    &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=MAX(0.0_r8,A(ij)*ramp)
          IF (iw.eq.1) THEN
            FORCES(ng)%Hwave(i,j)=cff
          ELSE
            FORCES(ng)%Hwave(i,j)=FORCES(ng)%Hwave(i,j)+                &
     &                            cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'SWANtoROMS Min/Max HSIGN   (m):     ',        &
     &                    range(1),range(2)
      END IF
!
!  Surface peak wave period.
!
      CALL AttrVect_exportRAttr(AttrVect_G(ng)%wav2ocn_AV, "RTP",       &
     &                          A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=MAX(0.0_r8,A(ij))
          IF (iw.eq.1) THEN
            FORCES(ng)%Pwave_top(i,j)=cff
          ELSE
            FORCES(ng)%Pwave_top(i,j)=FORCES(ng)%Pwave_top(i,j)+        &
     &                                cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'SWANtoROMS Min/Max RTP     (s):     ',        &
     &                    range(1),range(2)
      END IF
!
!  Bottom mean wave period.
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "TMBOT",    &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=MAX(0.0_r8,A(ij))
          IF (iw.eq.1) THEN
            FORCES(ng)%Pwave_bot(i,j)=cff
          ELSE
            FORCES(ng)%Pwave_bot(i,j)=FORCES(ng)%Pwave_bot(i,j)+        &
     &                                cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'SWANtoROMS Min/Max TMBOT   (s):     ',        &
     &                    range(1),range(2)
      END IF
!
!  Bottom orbital velocity (m/s).
!
      CALL AttrVect_exportRAttr(AttrVect_G(ng)%wav2ocn_AV, "UBOT",      &
     &                          A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=MAX(0.0_r8,A(ij)*ramp)
          IF (iw.eq.1) THEN
            FORCES(ng)%Uwave_rms(i,j)=cff
          ELSE
            FORCES(ng)%Uwave_rms(i,j)=FORCES(ng)%Uwave_rms(i,j)+        &
     &                                cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'SWANtoROMS Min/Max UBOT    (ms-1):  ',        &
     &                    range(1),range(2)
      END IF
!
!  Mean wave direction (radians).
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "DIRE",     &
     &                           A, Asize)
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "DIRN",     &
     &                           A1, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=ATAN2(A(ij),A1(ij))
          IF (cff.lt.0.0_r8) cff=cff+2.0_r8*pi
          IF (iw.eq.1) THEN
            FORCES(ng)%Dwave(i,j)=cff
          ELSE
            FORCES(ng)%Dwave(i,j)=FORCES(ng)%Dwave(i,j)+cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'SWANtoROMS Min/Max DIR     (rad):   ',        &
     &                    range(1),range(2)
      END IF
!
!  Peak wave direction (radians).
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "DIREP",    &
     &                           A, Asize)
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "DIRNP",    &
     &                           A1, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=ATAN2(A(ij),A1(ij))
          IF (cff.lt.0.0_r8) cff=cff+2.0_r8*pi
          IF (iw.eq.1) THEN
            FORCES(ng)%Dwavep(i,j)=cff
          ELSE
            FORCES(ng)%Dwavep(i,j)=FORCES(ng)%Dwavep(i,j)+cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'SWANtoROMS Min/Max DIRP    (rad):   ',        &
     &                    range(1),range(2)
      END IF
!
!  Wave length (m).
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "WLEN",     &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=MIN(Lwave_max,MAX(1.0_r8,A(ij)))
          IF (iw.eq.1) THEN
            FORCES(ng)%Lwave(i,j)=cff
          ELSE
            FORCES(ng)%Lwave(i,j)=FORCES(ng)%Lwave(i,j)+                &
     &                            cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'SWANtoROMS Min/Max WLEN    (m):     ',        &
     &                    range(1),range(2)
      END IF
!
!
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Dissip_break)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Dissip_fric)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Dissip_wcap)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Hwave)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Pwave_top)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Pwave_bot)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Uwave_rms)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Dwave)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Dwavep)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Lwave)
      END IF
!
!-----------------------------------------------------------------------
!  Exchange tile boundaries.
!-----------------------------------------------------------------------
!
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    FORCES(ng)%Dissip_break,                      &
     &                    FORCES(ng)%Dissip_fric,                       &
     &                    FORCES(ng)%Dissip_wcap)
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    FORCES(ng)%Hwave, FORCES(ng)%Pwave_top,       &
     &                    FORCES(ng)%Pwave_bot)
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    FORCES(ng)%Uwave_rms)
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    FORCES(ng)%Dwave, FORCES(ng)%Dwavep,          &
     &                    FORCES(ng)%Lwave)
!
!  Deallocate communication arrays.
!
      deallocate (A)
      deallocate (A1)
!
 10   FORMAT (' OCN2WAV_COUPLING - error while receiving fields from ', &
     &        a, i4)
 20   FORMAT (' OCN2WAV_COUPLING - error while sending fields to: ',    &
     &        a, i4)
      RETURN
      END SUBROUTINE ocnfwav_coupling_tile
      SUBROUTINE finalize_ocn2wav_coupling
!
!========================================================================
!                                                                       !
!  This routine finalizes ocean and wave models coupling data streams.  !
!                                                                       !
!========================================================================
      USE mod_scalars
      USE mct_coupler_params
!
!  Local variable declarations.
!
      integer :: ng, iw, MyError
!
!-----------------------------------------------------------------------
!  Deallocate MCT environment.
!-----------------------------------------------------------------------
!
      deallocate ( wavids )
      deallocate ( ocnids )
      DO ng=1,Nocn_grids
!       CALL GlobalSegMap_clean (GlobalSegMap_G(ng)%GSMapROMS,          &
!    &                           MyError)
        DO iw=1,Nwav_grids
          CALL Router_clean (Router_W(ng,iw)%ROMStoSWAN, MyError)
          CALL AttrVect_clean (AttrVect_G(ng)%ocn2wav_AV, MyError)
        END DO
      END DO
      RETURN
      END SUBROUTINE finalize_ocn2wav_coupling
      SUBROUTINE ocean_coupling (nl)
!
!=======================================================================
!                                                                      !
!  Determine which roms grids are going to exchange data to otehr      !
!  model grids and call those exchange.                                !
!                                                                      !
!=======================================================================
!
      USE mod_parallel
      USE mct_coupler_params
      USE mod_scalars
!
!  Imported variable definitions.
!
      integer, intent(in) :: nl
!
!  Local variable declarations.
!
      integer :: MyError, nprocs, tile
      integer :: ng, iw, ia, ig, ih, nlay, offset
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!  Exchange fields from ocn to wav every TI_OCN2WAV steps and
!                  from wav to ocn every TI_WAV2OCN steps.
!-----------------------------------------------------------------------
!
      IF (nl.eq.1) THEN
        DO iw=1,Nwav_grids
          DO nlay=1,NestLayers
            DO ig=1,GridsInLayer(nlay)
              ng=GridNumber(ig,nlay)
              offset=-1 !nlay-NestLayers
              IF (MOD(iic(1)+offset,nOCNFWAV(1,1)).eq.0) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL ocnfwav_coupling (ng, iw, tile)
                END DO
              END IF
            END DO
          END DO
        END DO
      END IF
!
      IF (nl.eq.1) THEN
        DO iw=1,Nwav_grids
          DO nlay=1,NestLayers
            DO ig=1,GridsInLayer(nlay)
              ng=GridNumber(ig,nlay)
              offset=-1 !nlay-NestLayers
              IF (MOD(iic(1)+offset,nOCN2WAV(1,1)).eq.0) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL ocn2wav_coupling (ng, iw, tile)
                END DO
              END IF
            END DO
          END DO
        END DO
      END IF
      RETURN
      END SUBROUTINE ocean_coupling
      END MODULE ocean_coupler_mod
