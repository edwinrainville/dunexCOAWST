      MODULE WAVES_COUPLER_MOD
!
!svn $Id: waves_coupler.F 756 2008-09-14 20:18:28Z jcwarner $
!==================================================== John C. Warner ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group      Hernan G. Arango   !
!   Licensed under a MIT/X style license                               !
!   See License_ROMS.txt                                               !
!=======================================================================
!                                                                      !
!  This module is used to communicate and exchange data between SWAN   !
!  other coupled model(s) using the Model Coupling Toolkit (MCT).      !
!                                                                      !
!=======================================================================
!
!  Componenet model registry.
!
      USE m_MCTWorld, ONLY : MCTWorld_init => init
      USE m_MCTWorld, ONLY : MCTWorld_clean => clean
!
!  Domain decompositin descriptor datatype and assocoiated methods.
!
      USE m_GlobalSegMap, ONLY : GlobalSegMap
      USE m_GlobalSegMap, ONLY : GlobalSegMap_init => init
      USE m_GlobalSegMap, ONLY : GlobalSegMap_lsize => lsize
      USE m_GlobalSegMap, ONLY : GlobalSegMap_clean => clean
      USE m_GlobalSegMap, ONLY : GlobalSegMap_Ordpnts => OrderedPoints
!
!  Field storage data types and associated methods.
!
      USE m_AttrVect, ONLY : AttrVect
      USE m_AttrVect, ONLY : AttrVect_init => init
      USE m_AttrVect, ONLY : AttrVect_zero => zero
      USE m_AttrVect, ONLY : AttrVect_clean => clean
      USE m_AttrVect, ONLY : AttrVect_indxR => indexRA
      USE m_AttrVect, ONLY : AttrVect_importRAttr => importRAttr
      USE m_AttrVect, ONLY : AttrVect_exportRAttr => exportRAttr
!
!  Intercomponent communitcations scheduler.
!
      USE m_Router, ONLY : Router
      USE m_Router, ONLY : Router_init => init
      USE m_Router, ONLY : Router_clean => clean
!
!  Intercomponent transfer.
!
      USE m_Transfer, ONLY : MCT_isend => isend
      USE m_Transfer, ONLY : MCT_irecv => irecv
      USE m_Transfer, ONLY : MCT_waitr => waitrecv
      USE m_Transfer, ONLY : MCT_waits => waitsend
!
      implicit none
!
      PRIVATE
      PUBLIC :: initialize_wav_coupling
      PUBLIC :: initialize_wav_routers
      PUBLIC :: wav2ocn_coupling
      PUBLIC :: wavfocn_coupling
      PUBLIC :: finalize_wav_coupling
!
!  Declarations.
!
      TYPE T_GlobalSegMap_G
        TYPE(GlobalSegMap) :: GSMapSWAN         ! GloabalSegMap variables
      END TYPE T_GlobalSegMap_G
      TYPE (T_GlobalSegMap_G), ALLOCATABLE :: GlobalSegMap_G(:)
      TYPE T_AttrVect_G
        TYPE(AttrVect) :: wav2ocn_AV            ! AttrVect variables
        TYPE(AttrVect) :: ocn2wav_AV
      END TYPE T_AttrVect_G
      TYPE (T_AttrVect_G), ALLOCATABLE :: AttrVect_G(:)
      TYPE T_Router_O
        type(Router)   :: SWANtoROMS            ! Router variables
      END TYPE T_Router_O
      TYPE (T_Router_O), ALLOCATABLE :: Router_O(:,:)
      CONTAINS
      SUBROUTINE INITIALIZE_WAV_COUPLING (ng)
!
!=======================================================================
!                                                                      !
!  Initialize waves and ocean models coupling stream.  This is the     !
!  training phase use to constuct  MCT  parallel interpolators and     !
!  stablish communication patterns.                                    !
!                                                                      !
!=======================================================================
!
      USE OCPCOMM4
      USE SWCOMM3
      USE SWCOMM4
      USE M_GENARR
      USE M_PARALL
      USE swan_iounits
      USE mct_coupler_params
!
      include 'mpif.h'
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: MyError, MyRank
      integer :: gsmsize, nprocs
      integer :: i, j, io, ia, Isize, Jsize, Asize
      integer :: nRows, nCols, num_sparse_elems
      integer :: cid, cad, MXCGLc
      character (len=70)  :: nc_name
      character (len=20)  :: to_add
      character (len=120) :: wostring
      character (len=120) :: owstring
      real :: cff
!      integer, dimension(2) :: src_grid_dims, dst_grid_dims
      integer, allocatable :: start(:), length(:)
!
!-----------------------------------------------------------------------
!  Begin initialization phase.
!-----------------------------------------------------------------------
!
!  Get communicator local rank and size.
!
      CALL mpi_comm_rank (WAV_COMM_WORLD, MyRank, MyError)
      CALL mpi_comm_size (WAV_COMM_WORLD, nprocs, MyError)
!
!  Initialize MCT coupled model registry.
!
      IF (ng.eq.1) THEN
        ALLOCATE(GlobalSegMap_G(Nwav_grids))
        ALLOCATE(AttrVect_G(Nwav_grids))
      END IF
!
      WAVid=wavids(ng)
      IF (Nwav_grids.gt.1) THEN
        CALL MCTWorld_init (N_mctmodels, MPI_COMM_WORLD,                &
     &                      WAV_COMM_WORLD,myids=wavids)
      ELSE
        CALL MCTWorld_init (N_mctmodels, MPI_COMM_WORLD,                &
     &                      WAV_COMM_WORLD,WAVid)
      END IF
!
!  Initialize a Global Segment Map for non-haloed transfer of data for
!  SWAN. Determine non-haloed start and length arrays for this
!  processor.
!
      IF (nprocs.eq.1) THEN
        IF (KREPTX.eq.1) THEN
          Isize=MXCGL+1
        ELSE
          Isize=MXCGL
        END IF
        Jsize=MYCGL
      ELSE
        IF (MXCGL.gt.MYCGL) THEN
          IF (KREPTX.eq.1) THEN
            Isize=MXC-IHALOX*IBLKAD(1)+1
          ELSE
            Isize=MXC-IHALOX*IBLKAD(1)
          END IF
          Jsize=MYC
        ELSE
          IF (KREPTX.eq.1) THEN
            Isize=MXC+1
          ELSE
            Isize=MXC
          END IF
          Jsize=MYC-IHALOY*IBLKAD(1)
        END IF
      END IF
!
      allocate ( start(Jsize) )
      allocate ( length(Jsize) )
!
!  Add offset of 1 for repeating x-dir grids.
!
      IF (KREPTX.eq.1) THEN
        MXCGLc=MXCGL+1
      ELSE
        MXCGLc=MXCGL
      END IF
      DO j=1,Jsize
        length(j)=Isize
        IF (MXCGL.gt.MYCGL) THEN
          IF (MyRank.eq.0) THEN
            start(j)=MXF+(j-1)*MXCGLc
          ELSE
            start(j)=MXF+(j-1)*MXCGLc+IHALOX
          END IF
        ELSE
          IF (MyRank.eq.0) THEN
            start(j)=MYF+(j-1)*MXCGLc
          ELSE
            start(j)=(MYF+IHALOY-1)*MXCGLc+1+(j-1)*MXCGLc
          END IF
        END IF
      END DO
      gsmsize=Isize*Jsize
!
      CALL GlobalSegMap_init (GlobalSegMap_G(ng)%GSMapSWAN, start,      &
     &                        length, 0, WAV_COMM_WORLD, WAVid)
      deallocate (start)
      deallocate (length)
!
!  Initialize attribute vector holding the export data code strings of
!  the wave model.
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
      CALL AttrVect_init(AttrVect_G(ng)%wav2ocn_AV,                     &
     &                   rList=TRIM(wostring),lsize=gsmsize)
      CALL AttrVect_zero(AttrVect_G(ng)%wav2ocn_AV)
!
!  Initialize attribute vector holding the export data code string of
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
      CALL AttrVect_init (AttrVect_G(ng)%ocn2wav_AV,                    &
     &                    rList=TRIM(owstring),lsize=gsmsize)
      CALL AttrVect_zero (AttrVect_G(ng)%ocn2wav_AV)
!
      RETURN
      END SUBROUTINE INITIALIZE_WAV_COUPLING
      SUBROUTINE INITIALIZE_WAV_ROUTERS
!
!=======================================================================
!                                                                      !
!  Initialize waves routers for wave model.                            !
!                                                                      !
!=======================================================================
!
      USE mct_coupler_params
      USE M_PARALL
!
!      include 'mpif.h'
!
!  Local variable declarations.
!
      integer :: MyError, MyRank
      integer :: ng, iw, ia
!
!  Initialize MCT Routers.
!
      ALLOCATE(Router_O(Nwav_grids,Nocn_grids))
!
!  Initialize a router to the ocean model component.
!
      DO ng=1,Nocn_grids
        DO iw=1,Nwav_grids
          OCNid=ocnids(ng)
          CALL Router_init (OCNid, GlobalSegMap_G(iw)%GSMapSWAN,        &
     &                      WAV_COMM_WORLD, Router_O(iw,ng)%SWANtoROMS)
        END DO
      END DO
      RETURN
      END SUBROUTINE INITIALIZE_WAV_ROUTERS
      SUBROUTINE WAV2OCN_COUPLING (MIP, NVOQP, VOQR, VOQ, IRQ,          &
     &                             IVTYPE, COMPDA, Numcouple, ng, io)
!
!=======================================================================
!                                                                      !
!  This subroutine reads and writes the coupling data streams between  !
!  ocean and wave models. Currently, the following data streams are    !
!  processed:                                                          !
!                                                                      !
!  Fields exported to the OCEAN model:                                 !
!                                                                      !
!     * Wave direction mean (degrees)                                  !
!     * Wave direction peak (degrees)                                  !
!     * Significant wave height (m)                                    !
!     * Average wave length (m)                                        !
!     * Surface wave relative peak period (s)                          !
!     * Bottom wave period (s)                                         !
!     * Percent of breakig waves (nondimensional)                      !
!     * Wave energy dissipation (W/m2)                                 !
!     * Wave bottom orbital velocity (m/s)                             !
!                                                                      !
!  Fields imported from the OCEAN Model:                               !
!                                                                      !
!     * Bathymetry, bottom elevation (m)                               !
!     * Free-surface, water surface elevation (m)                      !
!     * Depth integrated u-momentum (m/s)                              !
!     * Depth integrated v-momentum (m/s)                              !
!                                                                      !
!=======================================================================
!
      USE SWCOMM1
      USE SWCOMM3
      USE SWCOMM4
      USE OUTP_DATA
      USE M_PARALL
      USE M_GENARR
      USE M_MPI
      USE OCPCOMM4
!     USE mod_coupler_kinds
      USE mct_coupler_params
!
      implicit none
!
!  Imported variable declarations.
!
      integer :: MIP, IRQ, nvoqp, Numcouple, ng, io
      integer :: VOQR(NMOVAR), IVTYPE, IP, IPI, IX, IY
      real :: COMPDA(MCGRD,MCMVAR)
      real :: VOQ(MIP,NVOQP)
!
!  Local variable declarations.
!
      integer :: MyStatus, MyError, MySize, MyRank
      integer :: i, id, j, gsmsize, ierr, indx, Tag
      integer :: Istr, Iend, Jstr, Jend, MXCGLc, start
      integer :: Isize, Jsize, INDXG, NPROCS, OFFSET
      integer, pointer :: points(:)
      real(m8), pointer :: avdata(:)
      real(m8), pointer :: DIRE(:)
      real(m8), pointer :: DIRN(:)
      real(m8), pointer :: DIREP(:)
      real(m8), pointer :: DIRNP(:)
      character (len=40) :: code
!
!-----------------------------------------------------------------------
!  Send wave fields to ROMS.
!-----------------------------------------------------------------------
!
      CALL MPI_COMM_RANK (WAV_COMM_WORLD, MyRank, MyError)
      CALL MPI_COMM_SIZE (WAV_COMM_WORLD, NPROCS, MyError)
!
!  Get the number of grid point on this processor.
!
      gsmsize=GlobalSegMap_lsize(GlobalSegMap_G(ng)%GSMapSWAN,          &
     &                           WAV_COMM_WORLD)
!
!  Allocate attribute vector array used to export/import data.
!
      allocate ( avdata(gsmsize),stat=ierr )
      allocate ( DIRE(gsmsize),stat=ierr )
      allocate ( DIRN(gsmsize),stat=ierr )
      allocate ( DIREP(gsmsize),stat=ierr )
      allocate ( DIRNP(gsmsize),stat=ierr )
      avdata=0.0_m8
      DIRE=0.0_m8
      DIRN=0.0_m8
      DIREP=0.0_m8
      DIRNP=0.0_m8
!
!  Ask for points in this tile.
!
      CALL GlobalSegMap_Ordpnts (GlobalSegMap_G(ng)%GSMapSWAN,          &
     &                           MyRank, points)
!
!  Load SWAN exporting data into MCT storage buffers.  Since this
!  routine is called several times from main, only load field
!  according to the IVTYPE flag.  The data is exported using ROMS
!  definition for real kind m8.
!
      IF (KREPTX.eq.1) THEN
        IF (nprocs.eq.1) THEN
          Isize=MXCGL+1
          Jsize=MYCGL
        ELSE
          IF (MXCGL.gt.MYCGL) THEN
            Isize=MXC-IHALOX*IBLKAD(1)+1
            Jsize=MYC
          ELSE
            Isize=MXC+1
            Jsize=MYC-IHALOY*IBLKAD(1)
          END IF
        END IF
        MXCGLc=MXCGL
        IF (MXCGL.gt.MYCGL) THEN
          IF (MyRank.eq.0) THEN
            start=MXF+(j-1)*MXCGLc
          ELSE
            start=MXF+(j-1)*MXCGLc+IHALOX
          END IF
        ELSE
          IF (MyRank.eq.0) THEN
            start=MYF
          ELSE
            start=(MYF+IHALOY-1)*MXCGLc+1
          END IF
        END IF
!
        IP=start-1
        IPI=0
        DO j=1,Jsize
          DO i=1,Isize-1
            IP=IP+1
            IPI=IPI+1
            avdata(IPI)=REAL( VOQ(IP,VOQR(IVTYPE)),m8 )
          END DO
          IPI=IPI+1
        END DO
!
! now fill the far right
!
        IP=start-1
        IPI=0
        DO j=1,Jsize
          IP=IP+Isize-1
          IPI=IPI+Isize
          avdata(IPI)=REAL( VOQ(IP,VOQR(IVTYPE)),m8 )
        END DO
      ELSE
        DO IP=1,gsmsize
          avdata(IP)=REAL( VOQ(points(IP),VOQR(IVTYPE)),m8 )
        END DO
      END IF
!
      IF (IVTYPE.eq.54) THEN
        CALL AttrVect_importRAttr (AttrVect_G(ng)%wav2ocn_AV,           &
     &                             "DISBOT",avdata)
      ELSE IF (IVTYPE.eq.55) THEN
        CALL AttrVect_importRAttr (AttrVect_G(ng)%wav2ocn_AV,           &
     &                             "DISSURF",avdata)
      ELSE IF (IVTYPE.eq.56) THEN
        CALL AttrVect_importRAttr (AttrVect_G(ng)%wav2ocn_AV,           &
     &                             "DISWCAP",avdata)
      ELSE IF (IVTYPE.eq.10) THEN
        CALL AttrVect_importRAttr (AttrVect_G(ng)%wav2ocn_AV,           &
     &                             "HSIGN",avdata)
      ELSE IF (IVTYPE.eq.12) THEN
        CALL AttrVect_importRAttr (AttrVect_G(ng)%wav2ocn_AV,           &
     &                             "RTP",avdata)
      ELSE IF (IVTYPE.eq.50) THEN
        CALL AttrVect_importRAttr (AttrVect_G(ng)%wav2ocn_AV,           &
     &                             "TMBOT",avdata)
      ELSE IF (IVTYPE.eq.6) THEN
        CALL AttrVect_importRAttr (AttrVect_G(ng)%wav2ocn_AV,           &
     &                             "UBOT",avdata)
      ELSE IF (IVTYPE.eq.13) THEN
        DO IP=1,gsmsize
          DIRE(IP)=1.0*SIN(avdata(IP)*PI/180.0)
          DIRN(IP)=1.0*COS(avdata(IP)*PI/180.0)
        END DO
        CALL AttrVect_importRAttr (AttrVect_G(ng)%wav2ocn_AV,           &
     &                             "DIRE",DIRE)
        CALL AttrVect_importRAttr (AttrVect_G(ng)%wav2ocn_AV,           &
     &                             "DIRN",DIRN)
      ELSE IF (IVTYPE.eq.14) THEN
        DO IP=1,gsmsize
          DIREP(IP)=1.0*SIN(avdata(IP)*PI/180.0)
          DIRNP(IP)=1.0*COS(avdata(IP)*PI/180.0)
        END DO
        CALL AttrVect_importRAttr (AttrVect_G(ng)%wav2ocn_AV,           &
     &                             "DIREP",DIREP)
        CALL AttrVect_importRAttr (AttrVect_G(ng)%wav2ocn_AV,           &
     &                             "DIRNP",DIRNP)
      ELSE IF (IVTYPE.eq.17) THEN
        CALL AttrVect_importRAttr (AttrVect_G(ng)%wav2ocn_AV,           &
     &                             "WLEN",avdata)
      ELSE IF (IVTYPE.eq.71) THEN
        CALL AttrVect_importRAttr (AttrVect_G(ng)%wav2ocn_AV,           &
     &                             "WLENP",avdata)
      ELSE IF (IVTYPE.eq.8) THEN
        CALL AttrVect_importRAttr (AttrVect_G(ng)%wav2ocn_AV,           &
     &                             "QB",avdata)
      ELSE IF (IVTYPE.eq.16) THEN
        CALL AttrVect_importRAttr (AttrVect_G(ng)%wav2ocn_AV,           &
     &                             "WDSPR",avdata)
      ELSE IF (IVTYPE.eq.58) THEN
        CALL AttrVect_importRAttr (AttrVect_G(ng)%wav2ocn_AV,           &
     &                             "WQP",avdata)
!
      END IF
!
      IF (IRQ.eq.Numcouple) THEN
!
!-----------------------------------------------------------------------
!  Create a restart file.
!-----------------------------------------------------------------------
!
!jcw    CALL BACKUP (AC2, SPCSIG, SPCDIR, KGRPNT, XCGRID, YCGRID, ng)
!
!-----------------------------------------------------------------------
!  Send wave fields bundle to ocean model, ROMS.
!-----------------------------------------------------------------------
!
          Tag=io*100+0*10+ng
          CALL MCT_isend (AttrVect_G(ng)%wav2ocn_AV,                    &
     &                   Router_O(ng,io)%SWANtoROMS, Tag)
          CALL MCT_waits (Router_O(ng,io)%SWANtoROMS)
          IF (MyRank.EQ.0) THEN
            WRITE (SCREEN,36)' == SWAN grid ',ng,                        &
     &                       ' sent wave data to ROMS grid ', io
 36         FORMAT (a14,i2,a29,i2)
          END IF
          IF (MyError.ne.0) THEN
            WRITE (*,*)'coupling send fail swancplr, Error= ', MyError
            CALL FINALIZE_WAV_COUPLING(ng)
          END IF
      END IF
      deallocate (avdata, points, DIRE, DIRN, DIREP, DIRNP)
!
      RETURN
      END SUBROUTINE WAV2OCN_COUPLING
      SUBROUTINE WAVFOCN_COUPLING (COMPDA, ng, io)
!
!=======================================================================
!                                                                      !
!  This subroutine reads and writes the coupling data streams between  !
!  ocean and wave models. Currently, the following data streams are    !
!  processed:                                                          !
!                                                                      !
!  Fields exported to the OCEAN model:                                 !
!                                                                      !
!     * Mean wave direction (degrees)                                  !
!     * Peak wave direction (degrees)                                  !
!     * Significant wave height (m)                                    !
!     * Average wave length (m)                                        !
!     * Surface wave relative peak period (s)                          !
!     * Bottom wave period (s)                                         !
!     * Percent of breakig waves (nondimensional)                      !
!     * Wave energy dissipation (W/m2)                                 !
!     * Wave bottom orbital velocity (m/s)                             !
!                                                                      !
!  Fields imported from the OCEAN Model:                               !
!                                                                      !
!     * Bathymetry, bottom elevation (m)                               !
!     * Free-surface, water surface elevation (m)                      !
!     * Depth integrated u-momentum (m/s)                              !
!     * Depth integrated v-momentum (m/s)                              !
!                                                                      !
!=======================================================================
!
      USE SWCOMM1
      USE SWCOMM3
      USE SWCOMM4
      USE OUTP_DATA
      USE M_PARALL
      USE M_GENARR
      USE M_MPI
      USE OCPCOMM4
!     USE mod_coupler_kinds
      USE mct_coupler_params
!
      implicit none
!
!  Imported variable declarations.
!
      integer :: ng, io
      integer :: IP, IX, IY
      real :: COMPDA(MCGRD,MCMVAR)
!
!  Local variable declarations.
!
      integer :: MyStatus, MyError, MySize, MyRank
      integer :: i, id, j, gsmsize, ierr, indx, Tag
      integer :: Istr, Iend, Jstr, Jend
      integer :: Isize, Jsize, INDXG, NPROCS, OFFSET
      integer :: NUMTRANSFER, NNEIGH, HALOSIZE, NUMSENT, INB
      integer :: WHICHWAY, GDEST, GSRC, TAGT, TAGB, TAGR, TAGL
      integer :: TREQUEST,BREQUEST,RREQUEST,LREQUEST,MSIZE
      integer :: iddep, idwlv, idvlx, idvly, idveg
      integer :: idruf, idice
      integer, dimension(MPI_STATUS_SIZE,4) :: status
      real :: cff
      real, parameter ::  Large = 1.0E+20
      real, dimension(2) :: range
      real(m8), pointer :: avdata(:)
      real, pointer :: TEMPMCT(:,:)
      real, pointer :: GRECVT(:), GRECVB(:), GRECVR(:), GRECVL(:)
      real, pointer :: GSENDT(:), GSENDB(:), GSENDR(:), GSENDL(:)
      character (len=40) :: code
!
!-----------------------------------------------------------------------
!  Send wave fields to ROMS.
!-----------------------------------------------------------------------
!
      CALL MPI_COMM_RANK (WAV_COMM_WORLD, MyRank, MyError)
      CALL MPI_COMM_SIZE (WAV_COMM_WORLD, NPROCS, MyError)
!
!  Get the number of grid point on this processor.
!
      gsmsize=GlobalSegMap_lsize(GlobalSegMap_G(ng)%GSMapSWAN,          &
     &                           WAV_COMM_WORLD)
!
!  Allocate attribute vector array used to export/import data.
!
      allocate ( avdata(gsmsize),stat=ierr )
      avdata=0.0_m8
!
!-----------------------------------------------------------------------
!  Receive from ROMS: Depth, Water Level, VELX, and VELY.
!-----------------------------------------------------------------------
!
!  Schedule receiving field from ocean model.
!
      Tag=io*100+0*10+ng
      CALL MCT_irecv (AttrVect_G(ng)%ocn2wav_AV,                        &
     &                Router_O(ng,io)%SWANtoROMS, Tag)
!
!     Wait to make sure the OCN data has arrived.
!
      CALL MCT_waitr (AttrVect_G(ng)%ocn2wav_AV,                        &
     &                 Router_O(ng,io)%SWANtoROMS)
!
      IF (MyRank.EQ.0) THEN
        WRITE (SCREEN,35) ' == SWAN grid ',ng,                          &
     &                    ' recv data from ROMS grid ', io
      END IF
      IF (MyError.ne.0) THEN
       WRITE (*,*) 'coupling fail swancplr, MyStatus= ', MyError
        CALL FINALIZE_WAV_COUPLING(ng)
      END IF
 35   FORMAT (a14,i2,a26,i2)
!
! Pass the non-halo data from MCT into tempmct array.
!
        NNEIGH = IBLKAD(1)
        IF (NPROCS.eq.1) THEN
          Istr=1
          Iend=MXC
          Jstr=1
          Jend=MYC
        ELSE
          IF (MXCGL.GT.MYCGL) THEN
            IF (MyRank.eq.0) THEN
              Istr=1
            ELSE
              Istr=IHALOX+1
            END IF
            Isize=MXC-IHALOX*IBLKAD(1)
            Iend=Istr+Isize-1
            Jstr=1
            Jend=MYC
            HALOSIZE=IHALOX*MYC
          ELSE
            IF (MyRank.eq.0) THEN
              Jstr=1
            ELSE
              Jstr=IHALOY+1
            END IF
            Jsize=MYC-IHALOY*IBLKAD(1)
            Jend=Jstr+Jsize-1
            Istr=1
            Iend=MXC
            HALOSIZE=IHALOY*MXC
          END IF
        END IF
!
!  Determine the amount of fields and assign id numbers.
!
        IP=0
        IP=IP+1
        iddep=IP
!
        IP=IP+1
        idwlv=IP
!
        IP=IP+1
        idvlx=IP
!
        IP=IP+1
        idvly=IP
!
        IP=IP+1
        idruf=IP
!
!
        NUMTRANSFER=IP
!
 40     FORMAT (a36,1x,2(1pe14.6))
        allocate ( TEMPMCT(MXC*MYC,NUMTRANSFER),stat=ierr )
        TEMPMCT=0.0
!
!  Bottom elevation.
!
        CALL AttrVect_exportRAttr (AttrVect_G(ng)%ocn2wav_AV,           &
     &                             "DEPTH",avdata,gsmsize)
        range(1)= Large
        range(2)=-Large
        IP=0
        DO IY=Jstr,Jend
          DO IX=Istr,Iend
            IP=IP+1
            INDXG=(IY-1)*MXC+IX
            TEMPMCT(INDXG,iddep)=avdata(IP)
            range(1)=MIN(range(1),REAL(avdata(IP)))
            range(2)=MAX(range(2),REAL(avdata(IP)))
          END DO
          IF (KREPTX.EQ.1) THEN
            IP=IP+1
          END IF
        END DO
        CALL SWREDUCE(range(1),1,SWREAL,SWMIN)
        CALL SWREDUCE(range(2),1,SWREAL,SWMAX)
        IF (MyRank.eq.0) THEN
          write(SCREEN,40) 'ROMStoSWAN Min/Max DEPTH   (m):     ',      &
     &                      range(1),range(2)
        END IF
!
!  Water surface elevation.
!
        CALL AttrVect_exportRAttr (AttrVect_G(ng)%ocn2wav_AV,           &
     &                             "WLEV",avdata,gsmsize)
        range(1)= Large
        range(2)=-Large
        IP=0
        DO IY=Jstr,Jend
          DO IX=Istr,Iend
            IP=IP+1
            INDXG=(IY-1)*MXC+IX
            TEMPMCT(INDXG,idwlv)=avdata(IP)
            range(1)=MIN(range(1),REAL(avdata(IP)))
            range(2)=MAX(range(2),REAL(avdata(IP)))
          END DO
          IF (KREPTX.EQ.1) THEN
            IP=IP+1
          END IF
        END DO
        CALL SWREDUCE(range(1),1,SWREAL,SWMIN)
        CALL SWREDUCE(range(2),1,SWREAL,SWMAX)
        IF (MyRank.eq.0) THEN
          write(SCREEN,40) 'ROMStoSWAN Min/Max WLEV    (m):     ',      &
     &                      range(1),range(2)
        END IF
!
!  Depth-integrated u-velocity.
!
        CALL AttrVect_exportRAttr(AttrVect_G(ng)%ocn2wav_AV,            &

     &                             "VELX",avdata,gsmsize)
        range(1)= Large
        range(2)=-Large
        IP=0
        DO IY=Jstr,Jend
          DO IX=Istr,Iend
            IP=IP+1
            INDXG=(IY-1)*MXC+IX
            TEMPMCT(INDXG,idvlx)=avdata(IP)
            range(1)=MIN(range(1),REAL(avdata(IP)))
            range(2)=MAX(range(2),REAL(avdata(IP)))
          END DO
          IF (KREPTX.EQ.1) THEN
            IP=IP+1
          END IF
        END DO
        CALL SWREDUCE(range(1),1,SWREAL,SWMIN)
        CALL SWREDUCE(range(2),1,SWREAL,SWMAX)
        IF (MyRank.eq.0) THEN
          write(SCREEN,40) 'ROMStoSWAN Min/Max VELX    (ms-1):  ',      &
     &                      range(1),range(2)
        END IF
!
!  Depth-integrated v-velocity.
!
        CALL AttrVect_exportRAttr(AttrVect_G(ng)%ocn2wav_AV,            &
     &                            "VELY",avdata,gsmsize)
        range(1)= Large
        range(2)=-Large
        IP=0
        DO IY=Jstr,Jend
          DO IX=Istr,Iend
            IP=IP+1
            INDXG=(IY-1)*MXC+IX
            TEMPMCT(INDXG,idvly)=avdata(IP)
            range(1)=MIN(range(1),REAL(avdata(IP)))
            range(2)=MAX(range(2),REAL(avdata(IP)))
          END DO
          IF (KREPTX.EQ.1) THEN
            IP=IP+1
          END IF
        END DO
        CALL SWREDUCE(range(1),1,SWREAL,SWMIN)
        CALL SWREDUCE(range(2),1,SWREAL,SWMAX)
        IF (MyRank.eq.0) THEN
          write(SCREEN,40) 'ROMStoSWAN Min/Max VELY    (ms-1):  ',      &
     &                      range(1),range(2)
        END IF
!
!  Bottom roughness.
!
        CALL AttrVect_exportRAttr(AttrVect_G(ng)%ocn2wav_AV,            &

     &                            "ZO",avdata,gsmsize)
        range(1)= Large
        range(2)=-Large
        IP=0
        DO IY=Jstr,Jend
          DO IX=Istr,Iend
            IP=IP+1
            INDXG=(IY-1)*MXC+IX
            TEMPMCT(INDXG,idruf)=avdata(IP)
            range(1)=MIN(range(1),REAL(avdata(IP)))
            range(2)=MAX(range(2),REAL(avdata(IP)))
          END DO
          IF (KREPTX.EQ.1) THEN
            IP=IP+1
          END IF
        END DO
        CALL SWREDUCE(range(1),1,SWREAL,SWMIN)
        CALL SWREDUCE(range(2),1,SWREAL,SWMAX)
        IF (MyRank.eq.0) THEN
          write(SCREEN,40) 'ROMStoSWAN Min/Max ZO      (m):      ',     &

     &                      range(1),range(2)
        END IF
!
!  Pack and send halo regions to be exchanged with adjacent tiles.
!  IBLKAD contains the tile data.
!  WHICHWAY: [top, bot, right, left] = [1 2 3 4]
!
        IF (NPROCS.GT.1) THEN
          MSIZE=HALOSIZE*NUMTRANSFER
          IF (MXCGL.GT.MYCGL) THEN
            allocate ( GSENDR(MSIZE),stat=ierr )
            allocate ( GSENDL(MSIZE),stat=ierr )
            allocate ( GRECVR(MSIZE),stat=ierr )
            allocate ( GRECVL(MSIZE),stat=ierr )
            GSENDR=0.0
            GSENDL=0.0
            GRECVR=0.0
            GRECVL=0.0
          ELSE
            allocate ( GSENDT(MSIZE),stat=ierr )
            allocate ( GSENDB(MSIZE),stat=ierr )
            allocate ( GRECVT(MSIZE),stat=ierr )
            allocate ( GRECVB(MSIZE),stat=ierr )
            GSENDT=0.0
            GSENDB=0.0
            GRECVT=0.0
            GRECVB=0.0
          END IF
          TAGT=1
          TAGB=2
          TAGR=3
          TAGL=4
          DO INB=1,NNEIGH
            OFFSET=0
            WHICHWAY=IBLKAD(3*INB)
            DO NUMSENT=1,NUMTRANSFER
              IP=OFFSET
              IF (WHICHWAY.EQ.1) THEN
                DO IY=MYC-IHALOX-2,MYC-3
                  DO IX=1,MXC
                    IP=IP+1
                    INDXG=(IY-1)*MXC+IX
                    GSENDT(IP)=TEMPMCT(INDXG,NUMSENT)
                  END DO
                END DO
              ELSE IF (WHICHWAY.EQ.2) THEN
                DO IY=IHALOY+1,IHALOY+3
                  DO IX=1,MXC
                    IP=IP+1
                    INDXG=(IY-1)*MXC+IX
                    GSENDB(IP)=TEMPMCT(INDXG,NUMSENT)
                  END DO
                END DO
              ELSE IF (WHICHWAY.EQ.3) THEN
                DO IY=1,MYC
                  DO IX=MXC-IHALOX-2,MXC-3
                    IP=IP+1
                    INDXG=(IY-1)*MXC+IX
                    GSENDR(IP)=TEMPMCT(INDXG,NUMSENT)
                  END DO
                END DO
              ELSE IF (WHICHWAY.EQ.4) THEN
                DO IY=1,MYC
                  DO IX=IHALOX+1,IHALOX+3
                    IP=IP+1
                    INDXG=(IY-1)*MXC+IX
                    GSENDL(IP)=TEMPMCT(INDXG,NUMSENT)
                  END DO
                END DO
              END IF
              OFFSET=OFFSET+HALOSIZE
            END DO
          END DO
          DO INB=1,NNEIGH
            GSRC=IBLKAD(3*INB-1)-1
            WHICHWAY=IBLKAD(3*INB)
            IF (WHICHWAY.EQ.1) THEN
              CALL mpi_irecv (GRECVT,MSIZE,SWREAL,                      &
     &                        GSRC,TAGB,WAV_COMM_WORLD,TREQUEST,MyError)
            ELSE IF (WHICHWAY.EQ.2) THEN
              CALL mpi_irecv (GRECVB,MSIZE,SWREAL,                      &
     &                        GSRC,TAGT,WAV_COMM_WORLD,BREQUEST,MyError)
            ELSE IF (WHICHWAY.EQ.3) THEN
              CALL mpi_irecv (GRECVR,MSIZE,SWREAL,                      &
     &                        GSRC,TAGL,WAV_COMM_WORLD,RREQUEST,MyError)
            ELSE IF (WHICHWAY.EQ.4) THEN
              CALL mpi_irecv (GRECVL,MSIZE,SWREAL,                      &
     &                        GSRC,TAGR,WAV_COMM_WORLD,LREQUEST,MyError)
            END IF
          END DO
          DO INB=1,NNEIGH
            GDEST=IBLKAD(3*INB-1)-1
            WHICHWAY=IBLKAD(3*INB)
            IF (WHICHWAY.EQ.1) THEN
              CALL mpi_send (GSENDT,MSIZE,SWREAL,                       &
     &                       GDEST,TAGT,WAV_COMM_WORLD,MyError)
            ELSE IF (WHICHWAY.EQ.2) THEN
              CALL mpi_send (GSENDB,MSIZE,SWREAL,                       &
     &                       GDEST,TAGB,WAV_COMM_WORLD,MyError)
            ELSE IF (WHICHWAY.EQ.4) THEN
              CALL mpi_send (GSENDL,MSIZE,SWREAL,                       &
     &                       GDEST,TAGL,WAV_COMM_WORLD,MyError)
            ELSE IF (WHICHWAY.EQ.3) THEN
              CALL mpi_send (GSENDR,MSIZE,SWREAL,                       &
     &                       GDEST,TAGR,WAV_COMM_WORLD,MyError)
            END IF
          END DO
!
! Receive and unpack halo regions exchanged with adjacent tiles.
! [top, bot, right, left] = [1 2 3 4]
!
          DO INB=1,NNEIGH
            WHICHWAY=IBLKAD(3*INB)
            IF (WHICHWAY.EQ.1) THEN
              CALL mpi_wait (TREQUEST,status(1,1),MyError)
            ELSE IF (WHICHWAY.EQ.2) THEN
              CALL mpi_wait (BREQUEST,status(1,2),MyError)
            ELSE IF (WHICHWAY.EQ.3) THEN
              CALL mpi_wait (RREQUEST,status(1,3),MyError)
            ELSE IF (WHICHWAY.EQ.4) THEN
              CALL mpi_wait (LREQUEST,status(1,4),MyError)
            END IF
          END DO
!
          DO INB=1,NNEIGH
            OFFSET=0
            WHICHWAY=IBLKAD(3*INB)
            IF (WHICHWAY.EQ.1) THEN
              DO NUMSENT=1,NUMTRANSFER
                IP=OFFSET
                DO IY=MYC-2,MYC
                  DO IX=1,MXC
                    IP=IP+1
                    INDXG=(IY-1)*MXC+IX
                    TEMPMCT(INDXG,NUMSENT)=GRECVT(IP)
                  END DO
                END DO
                OFFSET=OFFSET+HALOSIZE
              END DO
            ELSE IF (WHICHWAY.EQ.2) THEN
              DO NUMSENT=1,NUMTRANSFER
                IP=OFFSET
                DO IY=1,IHALOY
                  DO IX=1,MXC
                    IP=IP+1
                    INDXG=(IY-1)*MXC+IX
                    TEMPMCT(INDXG,NUMSENT)=GRECVB(IP)
                  END DO
                END DO
                OFFSET=OFFSET+HALOSIZE
              END DO
            ELSE IF (WHICHWAY.EQ.3) THEN
              DO NUMSENT=1,NUMTRANSFER
                IP=OFFSET
                DO IY=1,MYC
                  DO IX=MXC-2,MXC
                    IP=IP+1
                    INDXG=(IY-1)*MXC+IX
                    TEMPMCT(INDXG,NUMSENT)=GRECVR(IP)
                  END DO
                END DO
                OFFSET=OFFSET+HALOSIZE
              END DO
            ELSE IF (WHICHWAY.EQ.4) THEN
              DO NUMSENT=1,NUMTRANSFER
                IP=OFFSET
                DO IY=1,MYC
                  DO IX=1,IHALOX
                    IP=IP+1
                    INDXG=(IY-1)*MXC+IX
                    TEMPMCT(INDXG,NUMSENT)=GRECVL(IP)
                  END DO
                END DO
                OFFSET=OFFSET+HALOSIZE
              END DO
            END IF
          END DO
          IF (MXCGL.GT.MYCGL) THEN
            deallocate (GRECVR,GRECVL,GSENDR,GSENDL)
          ELSE
            deallocate (GRECVT,GRECVB,GSENDT,GSENDB)
          END IF
        END IF
!
! Finally insert the full (MXC*MYC) TEMPMCT array into the SWAN
! array for DEPTH and computational array COMPDA. Only insert
! active (wet points) using array KGRPNT.
!
!
! Move values at 'present' time level 2 to 'old' time level 1.
! MCGRD = MXC*MYC+1-#masked cells.
! MXC = # cells x-dir in this tile including halox.
! MYC = # cells y-dir in this tile including haloy.
! COMPDA has only active wet points + 1.
!
        IF (io.eq.1) THEN
          DO INDX = 1, MCGRD
            COMPDA(INDX,JWLV1)=COMPDA(INDX,JWLV2)
            COMPDA(INDX,JVX1) =COMPDA(INDX,JVX2)
            COMPDA(INDX,JVY1) =COMPDA(INDX,JVY2)
            COMPDA(INDX,JFRC3)=COMPDA(INDX,JFRC2)
          END DO
        END IF
!
! Insert bot level, water level, velx, vely, fric, and winds 
! into SWAN arrays.
!
! If not using a BBL model, then determine the non-spatially 
! varying friction coeff to enter into the JFRC2 array.
!
        IP=0
        DO IY=1,MYC
          DO IX=1,MXC
            IP=IP+1
            INDX = KGRPNT(IX,IY)
            IF (INDX.GT.1) THEN
              IF (io.eq.1) THEN
                COMPDA(INDX,JFRC2)=MAX(REAL(TEMPMCT(IP,idruf)), 0.0001)
              ELSE
                COMPDA(INDX,JFRC2)=COMPDA(INDX,JFRC2)+                  &
     &                             MAX(REAL(TEMPMCT(IP,idruf)), 0.0001)
              END IF
!             cff=TEMPMCT(IP,iddep)+TEMPMCT(IP,idwlv)
!             IF (cff.gt.0.) THEN
                IF (io.eq.1) THEN
                  COMPDA(INDX,JWLV2)=TEMPMCT(IP,idwlv)
                  COMPDA(INDX,JVX2)=TEMPMCT(IP,idvlx)
                  COMPDA(INDX,JVY2)=TEMPMCT(IP,idvly)
                ELSE
                  COMPDA(INDX,JWLV2)=COMPDA(INDX,JWLV2)+                &
     &                               TEMPMCT(IP,idwlv)
                  COMPDA(INDX,JVX2)=COMPDA(INDX,JVX2)+TEMPMCT(IP,idvlx)
                  COMPDA(INDX,JVY2)=COMPDA(INDX,JVY2)+TEMPMCT(IP,idvly)
                END IF
!             END IF
            END IF
          END DO
        END DO
!
      deallocate (TEMPMCT)
      deallocate (avdata)
!
      RETURN
      END SUBROUTINE WAVFOCN_COUPLING
      SUBROUTINE FINALIZE_WAV_COUPLING(ng)
!
!=======================================================================
!                                                                    ===
!  This routines terminates execution during coupling error.         ===
!                                                                    ===
!=======================================================================
      USE mct_coupler_params
!
!  Local variable declarations.
!
      integer :: ng, io, ia, MyError
!
!-----------------------------------------------------------------------
!  Deallocate MCT environment.
!-----------------------------------------------------------------------
!
      DO io=1,Nocn_grids
        CALL Router_clean (Router_O(ng,io)%SWANtoROMS, MyError)
      END DO
      CALL AttrVect_clean (AttrVect_G(ng)%wav2ocn_AV, MyError)
      CALL GlobalSegMap_clean (GlobalSegMap_G(ng)%GSMapSWAN, MyError)
      END SUBROUTINE FINALIZE_WAV_COUPLING
      END MODULE WAVES_COUPLER_MOD
