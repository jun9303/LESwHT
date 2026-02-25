!==========================================
MODULE MOD_POSTAVG
!---- WRITE THE INPUTS
!-----GEOMETRY
  INTEGER*8, PARAMETER :: M1 = 641, M2 = 641, M3 = 641
  INTEGER*8, PARAMETER :: M1M = M1 - 1, M2M = M2 - 1, M3M = M3 - 1
  INTEGER*8 :: N1, N1M, N2, N2M, N3, N3M
  REAL*8    :: XL, YL, ZL, X(0:M1), Y(0:M2), Z(0:M3), XMP(0:M1), YMP(0:M2), ZMP(0:M3)
  INTEGER*8 :: IPV(M1M), JPV(M2M), KPV(M3M), IMV(M1M), JMV(M2M), KMV(M3M)
  REAL*8    :: FIXIL(M1M), FIXIU(M1M), FIXJL(M2M), FIXJU(M2M), FIXKL(M3M), FIXKU(M3M)
  REAL*8    :: SDX(0:M1),SDY(0:M2),SDZ(0:M3),VDX(M1),VDY(M2),VDZ(M3),SSDX(0:M1),SSDY(0:M2),SSDZ(0:M3),VVDX(M1),VVDY(M2),VVDZ(M3)

!------INPUT FILES OPTION!
!     : Refer to the routine READINPUTFILE and post_inst.in file for detailed explanation

  INTEGER*8    :: IBMON, ITOT, ILD2, IUVWP, IWXYZ, IALL, ISKIP, JSKIP, KSKIP, IHTRANS, IOUTFMT, IUNIGRID
  INTEGER*8    :: IPZ, IPX, NIP, NJP, NKP, NFLD, IND_FILM
  INTEGER*8    :: ISTART, IEND, JSTART, JEND, KSTART, KEND
  REAL*8       :: XIP(20), YJP(20), ZKP(20)
  CHARACTER*20 :: FLDNAME(100)
  CHARACTER*9  :: ftailijk
  CHARACTER*5  :: tfn1, tfn2
  CHARACTER*9  :: tfn3
  CHARACTER*10 :: tname1, tname2
  CHARACTER*14 :: tname3

!------VARIABLES
  INTEGER*8   :: IP(20), JP(20), KP(20), IUNI(M1M), JUNI(M2M), KUNI(M3M), INUMUNI, JNUMUNI, KNUMUNI
  REAL*8       :: XUNI(M1), YUNI(M2), ZUNI(M3), FACUNIX(M1M), FACUNIY(M2M), FACUNIZ(M3M)
  REAL*8, ALLOCATABLE :: UAVG(:, :, :), VAVG(:, :, :), WAVG(:, :, :), UIUJAVG(:, :, :, :)
  REAL*8, ALLOCATABLE :: PAVG(:, :, :), P2AVG(:, :, :), TAVG(:, :, :), T2AVG(:, :, :)
  REAL*8, ALLOCATABLE :: VORAVG(:, :, :, :), VOR2AVG(:, :, :, :), SSAVG(:, :, :)
  LOGICAL :: ARR_ALLOC = .FALSE.

!------TIME
  REAL*8      :: TIMEINIT, TIMEEND
  REAL*8      :: IHISTINIT, IHISTEND

END MODULE
!=======================================================================
!========================== MAIN PROCEDURE =============================
!=======================================================================
PROGRAM POST_2D
  USE MOD_POSTAVG
  IMPLICIT NONE
  INTEGER*8    :: N, L
  CHARACTER*20 :: DUMMY
  INTEGER*8    :: IPOINT, JPOINT, KPOINT

  IPZ = 1
  IPX = 1

  CALL READINPUTFILE                             ! Read 'post_inst.in' file
  CALL GEOM('../output/grid/grid.dat')           ! Section for grid geometry
  CALL ADJUST_POST_BOUNDS
  CALL ALLOC_POST_ARRAYS                         ! Allocate arrays using actual grid size

  CALL FINDIPJPKP                                ! Find indices of planes on which 2D fields will be printed
  IF (IUNIGRID .EQ. 1) CALL GRIDUNI

  DO L = 1, NFLD
    WRITE (*, *) ' '
    WRITE (*, 108)
    WRITE (*, *) 'WORKING ON ', FLDNAME(L)
    CALL PREFLD(FLDNAME(L))                     ! Read an input field file
    CALL DATAINIT
    CALL MAKEFHEAD                              ! Make head of an output field file

!======== Work on a plane and write a 2D output field
    DO N = 1, NIP                                 ! Work on a x-plane and write a 2D output field
      WRITE (*, *) '   WORKING ON X= ', XMP(IP(N))
      IPOINT = IP(N)
      CALL MAKEFTAIL(IPOINT, 1)                   ! Make tail of an output field file
      CALL OUTPUTI(IPOINT)                       ! Write an output field file (center velocity, pressure, vorticity, vorticity magnitude, and lambda2)
    END DO

    DO N = 1, NJP
      WRITE (*, *) '   WORKING ON Y= ', YMP(JP(N))
      JPOINT = JP(N)
      CALL MAKEFTAIL(JPOINT, 2)
      CALL OUTPUTJ(JPOINT)
    END DO

    DO N = 1, NKP
      WRITE (*, *) '   WORKING ON Z= ', ZMP(KP(N))
      KPOINT = KP(N)
      CALL MAKEFTAIL(KPOINT, 3)
      CALL OUTPUTK(KPOINT)
    END DO
!=========

!======== Work on the whole 3D field
    IF (ITOT .EQ. 1) THEN
      WRITE (*, *) '   WORKING ON 3D FIELD OUTPUT '
      CALL OUTPUT_3D
    END IF

    WRITE (*, 108)
    IND_FILM = IND_FILM + 1
  END DO
108 FORMAT('====================================================================')

  STOP
END

!=======================================================================
SUBROUTINE MAKEFTAIL(IDEX, LL)
!=======================================================================
!
!     Write a tail of output file name
!     ex) ftailijk = _i048.dat if a working plane is a x-plane with
!         I = 48 i.e. if writing a 2d output field file on I = 48
!
!=======================================================================
  USE MOD_POSTAVG
  IMPLICIT NONE
  INTEGER*8    :: IDEX, LL
  INTEGER*8    :: I1, I2, I3
  CHARACTER*1  :: NN1, NN2, NN3
  CHARACTER*2  :: IJK

  I3 = IDEX/100
  I2 = IDEX/10 - I3*10
  I1 = IDEX - I3*100 - I2*10
  NN3 = CHAR(I3 + 48)
  NN2 = CHAR(I2 + 48)
  NN1 = CHAR(I1 + 48)

  IF (LL .EQ. 1) THEN
    IJK = '_i'
  ELSEIF (LL .EQ. 2) THEN
    IJK = '_j'
  ELSEIF (LL .EQ. 3) THEN
    IJK = '_k'
  END IF

  ftailijk = IJK//NN3(1:1)//NN2(1:1)//NN1(1:1)//'.vtk'

  RETURN
END

!=======================================================================
SUBROUTINE MAKEFHEAD
!=======================================================================
!
!     Write a head of output file name
!     ex) If IND_FILM is 5 (in the post_inst.in), head names for the first
!         input fields are
!            2dfm_00005     (: 2d output for a certain plane)
!            3dfm_00005     (: 3d output)
!            2dfm_uni_00005 (: 2d output on a uniform grid)
!
!=======================================================================
  USE MOD_POSTAVG
  IMPLICIT NONE
  INTEGER*8    :: idg1, idg2, idg3, idg4, idg5

  idg1 = IND_FILM/10000
  idg2 = (IND_FILM - idg1*10000)/1000
  idg3 = (IND_FILM - idg1*10000 - idg2*1000)/100
  idg4 = (IND_FILM - idg1*10000 - idg2*1000 - idg3*100)/10
  idg5 = (IND_FILM - idg1*10000 - idg2*1000 - idg3*100 - idg4*10)

  tfn1 = '2dav_'
  tfn2 = '3dav_'
  tfn3 = '2dav_uni_'

  tname1 = tfn1//char(idg1 + 48)//char(idg2 + 48)//char(idg3 + 48)//char(idg4 + 48)//char(idg5 + 48)
  tname2 = tfn2//char(idg1 + 48)//char(idg2 + 48)//char(idg3 + 48)//char(idg4 + 48)//char(idg5 + 48)
  tname3 = tfn3//char(idg1 + 48)//char(idg2 + 48)//char(idg3 + 48)//char(idg4 + 48)//char(idg5 + 48)

  RETURN
END

!=======================================================================
SUBROUTINE READINPUTFILE
!=======================================================================
!
!     Read 'post_inst.in' file
!
!     Explanation for the 'post_inst.in': (Open the 'post_inst.in' and
!     read it together)
!========================     IBM OPTION      ==========================
! IBMON (IMMERSED BOUNDARY METHOD ON/OFF) : 0= OFF, 1= ON
!========================   2D POST OPTION    ==========================
! IUNIGRID (UNIFROM GRID DATA) : 0= OFF, 1= ON
! XYZ_position : 2D CUT FIELDS
!   1. The number of positions (0= OFF)
!   2. Position of cell where one want to print a 2D field
!    ex)
!   X_position
!   2
!   0.1
!   0.2
!     -> 2D flow fields at (the nearest points from) x = 0.1, 0.2
!=====================   3D POST OPTION    =============================
! ITOT  (3D fields)                  : 0=OFF, 1=ON
! ILD2  (Lambda2 criterion fields)   : 0=OFF, 1=ON
! IUVWP (Velocity & pressure fields) : 0=OFF, 1=ON
! IWXYZ (Vorticity fields)           : 0=OFF, 1=ON
! IALL  (Velocity, pressure, vorticity, lambda2 criterion) : 0=OFF, 1=ON
!
! ISKIP,JSKIP,KSKIP     : Grid interval in X,Y and Z directions
!  ex)
! ISKIP   JSKIP   KSKIP
! 1         3       1
!  -> Print all data in X and Z directions, and skip 2 data for every
!     3 cells in Y direction
!
! ISTART,JSTART,KSTART  : Initial grid indices in X,Y and Z directions
! IEND,JEND,KEND        : Final grid indices in X,Y and Z directions
! ANIFLD (NFLD)         : Number of fields
! INDFILM               : Initial number of field
!  ex)
! ANIFLD  IND_FILM
! 3       5
! fld010000
! fld020000
! fld022000
!  -> you load 3(ANIFLD) files 'fld010000','fld020000',and 'fld022000'
!  -> after post-processing,
!  file number of 'fld010000' becomes '2dfm_00005 ... .dat' ... etc
!  file number of 'fld020000' becomes '2dfm_00006 ... .dat' ... etc
!  file number of 'fld022000' becomes '2dfm_00007 ... .dat' ... etc
!=======================================================================
  USE MOD_POSTAVG
  IMPLICIT NONE
  INTEGER*8    :: N, L
  CHARACTER*20 :: DUMMY

  OPEN (2, FILE='post_avg.in')
  READ (2, *) DUMMY
  READ (2, *) IBMON, IHTRANS
  READ (2, *) DUMMY
  READ (2, *) IOUTFMT

  READ (2, *) DUMMY
  READ (2, *) DUMMY
  READ (2, *) IUNIGRID
  READ (2, *) DUMMY
  READ (2, *) NIP
  DO N = 1, NIP
    READ (2, *) XIP(N)
  END DO
  READ (2, *) DUMMY
  READ (2, *) NJP
  DO N = 1, NJP
    READ (2, *) YJP(N)
  END DO
  READ (2, *) DUMMY
  READ (2, *) NKP
  DO N = 1, NKP
    READ (2, *) ZKP(N)
  END DO
  READ (2, *) DUMMY
  READ (2, *) DUMMY
  READ (2, *) ITOT
  READ (2, *) DUMMY
  READ (2, *) ILD2, IUVWP, IWXYZ, IALL
  READ (2, *) DUMMY
  READ (2, *) ISKIP, JSKIP, KSKIP
  READ (2, *) DUMMY, DUMMY
  READ (2, *) ISTART, IEND
  READ (2, *) DUMMY, DUMMY
  READ (2, *) JSTART, JEND
  READ (2, *) DUMMY, DUMMY
  READ (2, *) KSTART, KEND
  READ (2, *) DUMMY, DUMMY
  READ (2, *) NFLD, IND_FILM
  DO L = 1, NFLD
    READ (2, *) FLDNAME(L)
  END DO
  CLOSE (2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  WRITE (*, 101) NIP
  WRITE (*, 102) NJP
  WRITE (*, 103) NKP
  DO N = 1, NIP
    WRITE (*, 104) N, XIP(N)
  END DO
  DO N = 1, NJP
    WRITE (*, 105) N, YJP(N)
  END DO
  DO N = 1, NKP
    WRITE (*, 106) N, ZKP(N)
  END DO
  WRITE (*, 107) ITOT
  WRITE (*, 108) ILD2, IUVWP, IWXYZ, IALL
  WRITE (*, 109) ISKIP, JSKIP, KSKIP
  WRITE (*, 110) ISTART, IEND
  WRITE (*, 111) JSTART, JEND
  WRITE (*, 112) KSTART, KEND
  WRITE (*, 113) NFLD, IND_FILM
  WRITE (*, 115) IUNIGRID
  WRITE (*, 116) IOUTFMT

101 FORMAT('# OF XPOSITION  = ', I5)
102 FORMAT('# OF YPOSITION  = ', I5)
103 FORMAT('# OF ZPOSITION  = ', I5)
104 FORMAT('X_POSITION ', I3, ' : ', F13.5)
105 FORMAT('Y_POSITION ', I3, ' : ', F13.5)
106 FORMAT('Z_POSITION ', I3, ' : ', F13.5)
107 FORMAT('ITOT   =', I5)
108 FORMAT('ILD2   =', I5, '  IUVWP =', I5, '  IWXYZ =', I5, '  IALL  =', I5)
109 FORMAT('ISKIP  =', I5, '  JSKIP =', I5, '  KSKIP =', I5)
110 FORMAT('ISTART =', I5, '  IEND =', I5)
111 FORMAT('JSTART =', I5, '  JEND =', I5)
112 FORMAT('KSTART =', I5, '  KEND =', I5)
113 FORMAT('NFLD   =', I5, '  IND_FILM =', I5)
115 FORMAT('IUNI_GRID =', I5)
116 FORMAT('IOUTFMT =', I5, '  (0:VTK, 1:TEC, 2:BOTH)')

  RETURN
END

!=======================================================================
SUBROUTINE ADJUST_POST_BOUNDS
!=======================================================================
  USE MOD_POSTAVG
  IMPLICIT NONE
  INTEGER*8 :: TMP

  IF (ISKIP .LE. 0) ISKIP = 1
  IF (JSKIP .LE. 0) JSKIP = 1
  IF (KSKIP .LE. 0) KSKIP = 1

  IF (ISTART .LE. 0) ISTART = 1
  IF (JSTART .LE. 0) JSTART = 1
  IF (KSTART .LE. 0) KSTART = 1

  IF (IEND .LE. 0) IEND = N1M
  IF (JEND .LE. 0) JEND = N2M
  IF (KEND .LE. 0) KEND = N3M

  ISTART = MAX(1_8, MIN(ISTART, N1M))
  IEND = MAX(1_8, MIN(IEND, N1M))
  JSTART = MAX(1_8, MIN(JSTART, N2M))
  JEND = MAX(1_8, MIN(JEND, N2M))
  KSTART = MAX(1_8, MIN(KSTART, N3M))
  KEND = MAX(1_8, MIN(KEND, N3M))

  IF (ISTART .GT. IEND) THEN
    TMP = ISTART
    ISTART = IEND
    IEND = TMP
  END IF
  IF (JSTART .GT. JEND) THEN
    TMP = JSTART
    JSTART = JEND
    JEND = TMP
  END IF
  IF (KSTART .GT. KEND) THEN
    TMP = KSTART
    KSTART = KEND
    KEND = TMP
  END IF

  WRITE (*, '(A,I5,A,I5)') 'AUTO ISTART =', ISTART, '  IEND =', IEND
  WRITE (*, '(A,I5,A,I5)') 'AUTO JSTART =', JSTART, '  JEND =', JEND
  WRITE (*, '(A,I5,A,I5)') 'AUTO KSTART =', KSTART, '  KEND =', KEND

  RETURN
END

!=======================================================================
SUBROUTINE ALLOC_POST_ARRAYS
!=======================================================================
  USE MOD_POSTAVG
  IMPLICIT NONE

  IF (ARR_ALLOC) RETURN

  ALLOCATE (UAVG(N1M, N2M, N3M), VAVG(N1M, N2M, N3M), WAVG(N1M, N2M, N3M))
  ALLOCATE (UIUJAVG(N1M, N2M, N3M, 6))
  ALLOCATE (PAVG(N1M, N2M, N3M), P2AVG(N1M, N2M, N3M))
  ALLOCATE (VORAVG(N1M, N2M, N3M, 3), VOR2AVG(N1M, N2M, N3M, 6), SSAVG(N1M, N2M, N3M))
  ALLOCATE (TAVG(N1M, N2M, N3M), T2AVG(N1M, N2M, N3M))

  UAVG = 0.D0
  VAVG = 0.D0
  WAVG = 0.D0
  UIUJAVG = 0.D0
  PAVG = 0.D0
  P2AVG = 0.D0
  VORAVG = 0.D0
  VOR2AVG = 0.D0
  SSAVG = 0.D0
  TAVG = 0.D0
  T2AVG = 0.D0

  ARR_ALLOC = .TRUE.

  RETURN
END

!=======================================================================
SUBROUTINE DATAINIT
!=======================================================================
!
!     Calaulate cell center velocities, vorticities and lambda 2 values
!
!=======================================================================
  USE MOD_POSTAVG
  IMPLICIT NONE
  INTEGER*8  :: I, J, K, L
  REAL*8     :: FUNCBODY

  ! Let U,V,W = 0 if a given point is inside the body
  DO K = 1, N3M
    DO J = 1, N2M
      DO I = 1, N1M
        IF (FUNCBODY(X(I), YMP(J), ZMP(K), 0) .LE. 1.E-10) THEN
          UAVG(I, J, K) = 0.
        END IF
        IF (FUNCBODY(XMP(I), Y(J), ZMP(K), 0) .LE. 1.E-10) THEN
          VAVG(I, J, K) = 0.
        END IF
        IF (FUNCBODY(XMP(I), YMP(J), Z(K), 0) .LE. 1.E-10) THEN
          WAVG(I, J, K) = 0.
        END IF
      END DO
    END DO
  END DO

  IF (IBMON .EQ. 1) THEN
    DO K = 1, N3M
      DO J = 1, N2M
        DO I = 1, N1M

          ! Modify flow variables at (XMP(I),YMP(J),ZMP(K)) if (XMP(I),YMP(J),ZMP(K)) is inside the body
          IF (FUNCBODY(XMP(I), YMP(J), ZMP(K), 0) .LE. 1.E-10) THEN
            UAVG(I, J, K) = 0.
            VAVG(I, J, K) = 0.
            WAVG(I, J, K) = 0.
            DO L = 1, 6
              IF (L .LE. 3) VORAVG(I, J, K, L) = 0.
              UIUJAVG(I, J, K, L) = 0.
              VOR2AVG(I, J, K, L) = 0.
            END DO
            PAVG(I, J, K) = 1000.
            P2AVG(I, J, K) = 1000.
            SSAVG(I, J, K) = 0.
          END IF

        END DO
      END DO
    END DO
  END IF

  RETURN
END
!=======================================================================
SUBROUTINE OUTPUTI(IPOINT)
!=======================================================================
!
!     Write a 2D flow field on a x-plane (I = IPOINT)
!
!=======================================================================
  USE MOD_POSTAVG
  IMPLICIT NONE
  INTEGER*8    :: IPOINT
  INTEGER*8    :: I, J, K, L
  INTEGER*8    :: II, JJ, KK
  REAL*8       :: VXX(M2M, M3M), WXX(M2M, M3M)
  II = IPOINT

  IF (IOUTFMT .EQ. 1) THEN
    CALL OUTPUTI_TEC(IPOINT)
    RETURN
  END IF

  OPEN (102, FILE='../output/post_avg/'//tname1//ftailijk)
  WRITE (102, '(A)') '# vtk DataFile Version 3.0'
  WRITE (102, '(A)') 'LESwHT averaged x-plane'
  WRITE (102, '(A)') 'ASCII'
  WRITE (102, '(A)') 'DATASET STRUCTURED_GRID'
  WRITE (102, *) 'DIMENSIONS ', N3M, N2M, 1
  WRITE (102, *) 'POINTS ', N3M*N2M, ' double'
  DO J = 1, N2M
    DO K = 1, N3M
      WRITE (102, *) XMP(II), YMP(J), ZMP(K)
    END DO
  END DO
  WRITE (102, *) 'POINT_DATA ', N3M*N2M
  WRITE (102, '(A)') 'VECTORS mean_velocity double'
  DO J = 1, N2M
    DO K = 1, N3M
      WRITE (102, *) UAVG(II, J, K), VAVG(II, J, K), WAVG(II, J, K)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS pavg double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO K = 1, N3M
      WRITE (102, *) PAVG(II, J, K)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS p2avg double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO K = 1, N3M
      WRITE (102, *) P2AVG(II, J, K)
    END DO
  END DO
  WRITE (102, '(A)') 'VECTORS mean_vorticity double'
  DO J = 1, N2M
    DO K = 1, N3M
      WRITE (102, *) (VORAVG(II, J, K, L), L=1, 3)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS uxux double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO K = 1, N3M
      WRITE (102, *) UIUJAVG(II, J, K, 1)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS uxuy double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO K = 1, N3M
      WRITE (102, *) UIUJAVG(II, J, K, 2)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS uxuz double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO K = 1, N3M
      WRITE (102, *) UIUJAVG(II, J, K, 3)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS uyuy double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO K = 1, N3M
      WRITE (102, *) UIUJAVG(II, J, K, 4)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS uyuz double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO K = 1, N3M
      WRITE (102, *) UIUJAVG(II, J, K, 5)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS uzuz double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO K = 1, N3M
      WRITE (102, *) UIUJAVG(II, J, K, 6)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS wxwx double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO K = 1, N3M
      WRITE (102, *) VOR2AVG(II, J, K, 1)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS wxwy double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO K = 1, N3M
      WRITE (102, *) VOR2AVG(II, J, K, 2)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS wxwz double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO K = 1, N3M
      WRITE (102, *) VOR2AVG(II, J, K, 3)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS wywy double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO K = 1, N3M
      WRITE (102, *) VOR2AVG(II, J, K, 4)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS wywz double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO K = 1, N3M
      WRITE (102, *) VOR2AVG(II, J, K, 5)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS wzwz double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO K = 1, N3M
      WRITE (102, *) VOR2AVG(II, J, K, 6)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS ssavg double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO K = 1, N3M
      WRITE (102, *) SSAVG(II, J, K)
    END DO
  END DO
  IF (IHTRANS .EQ. 1) THEN
    WRITE (102, '(A)') 'SCALARS tavg double 1'
    WRITE (102, '(A)') 'LOOKUP_TABLE default'
    DO J = 1, N2M
      DO K = 1, N3M
        WRITE (102, *) TAVG(II, J, K)
      END DO
    END DO
    WRITE (102, '(A)') 'SCALARS t2avg double 1'
    WRITE (102, '(A)') 'LOOKUP_TABLE default'
    DO J = 1, N2M
      DO K = 1, N3M
        WRITE (102, *) T2AVG(II, J, K)
      END DO
    END DO
  END IF
  CLOSE (102)

  IF (IUNIGRID .EQ. 1) THEN
    DO K = 1, N3M
      DO J = 1, N2M
        JJ = JUNI(J)
        KK = KUNI(K)

        VXX(J, K) = &
          FACUNIY(J)*FACUNIZ(K)*VAVG(II, JJ + 1, KK + 1) &
          + FACUNIY(J)*(1.-FACUNIZ(K))*VAVG(II, JJ + 1, KK) &
          + (1.-FACUNIY(J))*FACUNIZ(K)*VAVG(II, JJ, KK + 1) &
          + (1.-FACUNIY(J))*(1.-FACUNIZ(K))*VAVG(II, JJ, KK)

        WXX(J, K) = &
          FACUNIY(J)*FACUNIZ(K)*WAVG(II, JJ + 1, KK + 1) &
          + FACUNIY(J)*(1.-FACUNIZ(K))*WAVG(II, JJ + 1, KK) &
          + (1.-FACUNIY(J))*FACUNIZ(K)*WAVG(II, JJ, KK + 1) &
          + (1.-FACUNIY(J))*(1.-FACUNIZ(K))*WAVG(II, JJ, KK)
      END DO
    END DO

    OPEN (104, FILE='../output/post_avg/'//tname3//ftailijk)
    WRITE (104, '(A)') '# vtk DataFile Version 3.0'
    WRITE (104, '(A)') 'LESwHT averaged x-plane uniform'
    WRITE (104, '(A)') 'ASCII'
    WRITE (104, '(A)') 'DATASET STRUCTURED_GRID'
    WRITE (104, *) 'DIMENSIONS ', KNUMUNI, JNUMUNI, 1
    WRITE (104, *) 'POINTS ', KNUMUNI*JNUMUNI, ' double'
    DO J = 1, N2M, 2
      DO K = 1, N3M, 2
        WRITE (104, *) XMP(II), YUNI(J), ZUNI(K)
      END DO
    END DO
    WRITE (104, *) 'POINT_DATA ', KNUMUNI*JNUMUNI
    WRITE (104, '(A)') 'SCALARS vavg double 1'
    WRITE (104, '(A)') 'LOOKUP_TABLE default'
    DO J = 1, N2M, 2
      DO K = 1, N3M, 2
        WRITE (104, *) VXX(J, K)
      END DO
    END DO
    WRITE (104, '(A)') 'SCALARS wavg double 1'
    WRITE (104, '(A)') 'LOOKUP_TABLE default'
    DO J = 1, N2M, 2
      DO K = 1, N3M, 2
        WRITE (104, *) WXX(J, K)
      END DO
    END DO
    CLOSE (104)
  END IF

  IF (IOUTFMT .EQ. 2) CALL OUTPUTI_TEC(IPOINT)

  RETURN
END
!=======================================================================
SUBROUTINE OUTPUTJ(JPOINT)
!=======================================================================
!
!     Write a 2D flow field on a y-plane (J = JPOINT)
!
!=======================================================================
  USE MOD_POSTAVG
  IMPLICIT NONE
  INTEGER*8    :: JPOINT
  INTEGER*8    :: I, J, K, L
  INTEGER*8    :: II, JJ, KK, IIP, KKP
  REAL*8       :: UXX(M1M, M3M), WXX(M1M, M3M)

  JJ = JPOINT

  IF (IOUTFMT .EQ. 1) THEN
    CALL OUTPUTJ_TEC(JPOINT)
    RETURN
  END IF
  OPEN (102, FILE='../output/post_avg/'//tname1//ftailijk)
  WRITE (102, '(A)') '# vtk DataFile Version 3.0'
  WRITE (102, '(A)') 'LESwHT averaged y-plane'
  WRITE (102, '(A)') 'ASCII'
  WRITE (102, '(A)') 'DATASET STRUCTURED_GRID'
  WRITE (102, *) 'DIMENSIONS ', N1M, N3M, 1
  WRITE (102, *) 'POINTS ', N1M*N3M, ' double'
  DO K = 1, N3M
    DO I = 1, N1M
      WRITE (102, *) XMP(I), YMP(JJ), ZMP(K)
    END DO
  END DO
  WRITE (102, *) 'POINT_DATA ', N1M*N3M
  WRITE (102, '(A)') 'VECTORS mean_velocity double'
  DO K = 1, N3M
    DO I = 1, N1M
      WRITE (102, *) UAVG(I, JJ, K), VAVG(I, JJ, K), WAVG(I, JJ, K)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS pavg double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO K = 1, N3M
    DO I = 1, N1M
      WRITE (102, *) PAVG(I, JJ, K)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS p2avg double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO K = 1, N3M
    DO I = 1, N1M
      WRITE (102, *) P2AVG(I, JJ, K)
    END DO
  END DO
  WRITE (102, '(A)') 'VECTORS mean_vorticity double'
  DO K = 1, N3M
    DO I = 1, N1M
      WRITE (102, *) (VORAVG(I, JJ, K, L), L=1, 3)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS uxux double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO K = 1, N3M
    DO I = 1, N1M
      WRITE (102, *) UIUJAVG(I, JJ, K, 1)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS uxuy double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO K = 1, N3M
    DO I = 1, N1M
      WRITE (102, *) UIUJAVG(I, JJ, K, 2)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS uxuz double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO K = 1, N3M
    DO I = 1, N1M
      WRITE (102, *) UIUJAVG(I, JJ, K, 3)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS uyuy double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO K = 1, N3M
    DO I = 1, N1M
      WRITE (102, *) UIUJAVG(I, JJ, K, 4)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS uyuz double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO K = 1, N3M
    DO I = 1, N1M
      WRITE (102, *) UIUJAVG(I, JJ, K, 5)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS uzuz double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO K = 1, N3M
    DO I = 1, N1M
      WRITE (102, *) UIUJAVG(I, JJ, K, 6)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS wxwx double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO K = 1, N3M
    DO I = 1, N1M
      WRITE (102, *) VOR2AVG(I, JJ, K, 1)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS wxwy double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO K = 1, N3M
    DO I = 1, N1M
      WRITE (102, *) VOR2AVG(I, JJ, K, 2)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS wxwz double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO K = 1, N3M
    DO I = 1, N1M
      WRITE (102, *) VOR2AVG(I, JJ, K, 3)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS wywy double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO K = 1, N3M
    DO I = 1, N1M
      WRITE (102, *) VOR2AVG(I, JJ, K, 4)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS wywz double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO K = 1, N3M
    DO I = 1, N1M
      WRITE (102, *) VOR2AVG(I, JJ, K, 5)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS wzwz double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO K = 1, N3M
    DO I = 1, N1M
      WRITE (102, *) VOR2AVG(I, JJ, K, 6)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS ssavg double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO K = 1, N3M
    DO I = 1, N1M
      WRITE (102, *) SSAVG(I, JJ, K)
    END DO
  END DO
  IF (IHTRANS .EQ. 1) THEN
    WRITE (102, '(A)') 'SCALARS tavg double 1'
    WRITE (102, '(A)') 'LOOKUP_TABLE default'
    DO K = 1, N3M
      DO I = 1, N1M
        WRITE (102, *) TAVG(I, JJ, K)
      END DO
    END DO
    WRITE (102, '(A)') 'SCALARS t2avg double 1'
    WRITE (102, '(A)') 'LOOKUP_TABLE default'
    DO K = 1, N3M
      DO I = 1, N1M
        WRITE (102, *) T2AVG(I, JJ, K)
      END DO
    END DO
  END IF
  CLOSE (102)

  IF (IUNIGRID .EQ. 1) THEN
    DO K = 1, N3M
      DO I = 1, N1M
        II = IUNI(I)
        KK = KUNI(K)
        IIP = II + 1
        KKP = KK + 1

        UXX(I, K) = &
          FACUNIX(I)*FACUNIZ(K)*UAVG(IIP, JJ, KKP) &
          + FACUNIX(I)*(1.-FACUNIZ(K))*UAVG(IIP, JJ, KK) &
          + (1.-FACUNIX(I))*FACUNIZ(K)*UAVG(II, JJ, KKP) &
          + (1.-FACUNIX(I))*(1.-FACUNIZ(K))*UAVG(II, JJ, KK)

        WXX(I, K) = &
          FACUNIX(I)*FACUNIZ(K)*WAVG(IIP, JJ, KKP) &
          + FACUNIX(I)*(1.-FACUNIZ(K))*WAVG(IIP, JJ, KK) &
          + (1.-FACUNIX(I))*FACUNIZ(K)*WAVG(II, JJ, KKP) &
          + (1.-FACUNIX(I))*(1.-FACUNIZ(K))*WAVG(II, JJ, KK)
      END DO
    END DO

    OPEN (104, FILE='../output/post_avg/'//tname3//ftailijk)
    WRITE (104, '(A)') '# vtk DataFile Version 3.0'
    WRITE (104, '(A)') 'LESwHT averaged y-plane uniform'
    WRITE (104, '(A)') 'ASCII'
    WRITE (104, '(A)') 'DATASET STRUCTURED_GRID'
    WRITE (104, *) 'DIMENSIONS ', INUMUNI, KNUMUNI, 1
    WRITE (104, *) 'POINTS ', INUMUNI*KNUMUNI, ' double'
    DO K = 1, N3M, 2
      DO I = 1, N1M, 2
        WRITE (104, *) XUNI(I), YMP(JJ), ZUNI(K)
      END DO
    END DO
    WRITE (104, *) 'POINT_DATA ', INUMUNI*KNUMUNI
    WRITE (104, '(A)') 'SCALARS uavg double 1'
    WRITE (104, '(A)') 'LOOKUP_TABLE default'
    DO K = 1, N3M, 2
      DO I = 1, N1M, 2
        WRITE (104, *) UXX(I, K)
      END DO
    END DO
    WRITE (104, '(A)') 'SCALARS wavg double 1'
    WRITE (104, '(A)') 'LOOKUP_TABLE default'
    DO K = 1, N3M, 2
      DO I = 1, N1M, 2
        WRITE (104, *) WXX(I, K)
      END DO
    END DO
    CLOSE (104)
  END IF

  IF (IOUTFMT .EQ. 2) CALL OUTPUTJ_TEC(JPOINT)

  RETURN
END

!=======================================================================
SUBROUTINE OUTPUTK(KPOINT)
!=======================================================================
!
!     Write a 2D flow field on a z-plane (K = KPOINT)
!
!=======================================================================
  USE MOD_POSTAVG
  IMPLICIT NONE
  INTEGER*8    :: KPOINT
  INTEGER*8    :: I, J, K, L
  INTEGER*8    :: II, JJ, KK, IIP, JJP
  REAL*8       :: UXX(M1M, M2M), VXX(M1M, M2M)

  KK = KPOINT

  IF (IOUTFMT .EQ. 1) THEN
    CALL OUTPUTK_TEC(KPOINT)
    RETURN
  END IF

  OPEN (102, FILE='../output/post_avg/'//tname1//ftailijk)
  WRITE (102, '(A)') '# vtk DataFile Version 3.0'
  WRITE (102, '(A)') 'LESwHT averaged z-plane'
  WRITE (102, '(A)') 'ASCII'
  WRITE (102, '(A)') 'DATASET STRUCTURED_GRID'
  WRITE (102, *) 'DIMENSIONS ', N1M, N2M, 1
  WRITE (102, *) 'POINTS ', N1M*N2M, ' double'
  DO J = 1, N2M
    DO I = 1, N1M
      WRITE (102, *) XMP(I), YMP(J), ZMP(KK)
    END DO
  END DO
  WRITE (102, *) 'POINT_DATA ', N1M*N2M
  WRITE (102, '(A)') 'VECTORS mean_velocity double'
  DO J = 1, N2M
    DO I = 1, N1M
      WRITE (102, *) UAVG(I, J, KK), VAVG(I, J, KK), WAVG(I, J, KK)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS pavg double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO I = 1, N1M
      WRITE (102, *) PAVG(I, J, KK)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS p2avg double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO I = 1, N1M
      WRITE (102, *) P2AVG(I, J, KK)
    END DO
  END DO
  WRITE (102, '(A)') 'VECTORS mean_vorticity double'
  DO J = 1, N2M
    DO I = 1, N1M
      WRITE (102, *) (VORAVG(I, J, KK, L), L=1, 3)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS uxux double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO I = 1, N1M
      WRITE (102, *) UIUJAVG(I, J, KK, 1)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS uxuy double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO I = 1, N1M
      WRITE (102, *) UIUJAVG(I, J, KK, 2)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS uxuz double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO I = 1, N1M
      WRITE (102, *) UIUJAVG(I, J, KK, 3)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS uyuy double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO I = 1, N1M
      WRITE (102, *) UIUJAVG(I, J, KK, 4)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS uyuz double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO I = 1, N1M
      WRITE (102, *) UIUJAVG(I, J, KK, 5)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS uzuz double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO I = 1, N1M
      WRITE (102, *) UIUJAVG(I, J, KK, 6)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS wxwx double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO I = 1, N1M
      WRITE (102, *) VOR2AVG(I, J, KK, 1)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS wxwy double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO I = 1, N1M
      WRITE (102, *) VOR2AVG(I, J, KK, 2)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS wxwz double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO I = 1, N1M
      WRITE (102, *) VOR2AVG(I, J, KK, 3)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS wywy double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO I = 1, N1M
      WRITE (102, *) VOR2AVG(I, J, KK, 4)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS wywz double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO I = 1, N1M
      WRITE (102, *) VOR2AVG(I, J, KK, 5)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS wzwz double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO I = 1, N1M
      WRITE (102, *) VOR2AVG(I, J, KK, 6)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS ssavg double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO I = 1, N1M
      WRITE (102, *) SSAVG(I, J, KK)
    END DO
  END DO
  IF (IHTRANS .EQ. 1) THEN
    WRITE (102, '(A)') 'SCALARS tavg double 1'
    WRITE (102, '(A)') 'LOOKUP_TABLE default'
    DO J = 1, N2M
      DO I = 1, N1M
        WRITE (102, *) TAVG(I, J, KK)
      END DO
    END DO
    WRITE (102, '(A)') 'SCALARS t2avg double 1'
    WRITE (102, '(A)') 'LOOKUP_TABLE default'
    DO J = 1, N2M
      DO I = 1, N1M
        WRITE (102, *) T2AVG(I, J, KK)
      END DO
    END DO
  END IF
  CLOSE (102)

  IF (IUNIGRID .EQ. 1) THEN
    DO J = 1, N2M
      DO I = 1, N1M
        II = IUNI(I)
        JJ = JUNI(J)
        IIP = II + 1
        JJP = JJ + 1

        UXX(I, J) = &
          FACUNIX(I)*FACUNIY(J)*UAVG(IIP, JJP, KK) &
          + FACUNIX(I)*(1.-FACUNIY(J))*UAVG(IIP, JJ, KK) &
          + (1.-FACUNIX(I))*FACUNIY(J)*UAVG(II, JJP, KK) &
          + (1.-FACUNIX(I))*(1.-FACUNIY(J))*UAVG(II, JJ, KK)

        VXX(I, J) = &
          FACUNIX(I)*FACUNIY(J)*VAVG(IIP, JJP, KK) &
          + FACUNIX(I)*(1.-FACUNIY(J))*VAVG(IIP, JJ, KK) &
          + (1.-FACUNIX(I))*FACUNIY(J)*VAVG(II, JJP, KK) &
          + (1.-FACUNIX(I))*(1.-FACUNIY(J))*VAVG(II, JJ, KK)
      END DO
    END DO

    OPEN (104, FILE='../output/post_avg/'//tname3//ftailijk)
    WRITE (104, '(A)') '# vtk DataFile Version 3.0'
    WRITE (104, '(A)') 'LESwHT averaged z-plane uniform'
    WRITE (104, '(A)') 'ASCII'
    WRITE (104, '(A)') 'DATASET STRUCTURED_GRID'
    WRITE (104, *) 'DIMENSIONS ', INUMUNI, JNUMUNI, 1
    WRITE (104, *) 'POINTS ', INUMUNI*JNUMUNI, ' double'
    DO J = 1, N2M, 2
      DO I = 1, N1M, 2
        WRITE (104, *) XUNI(I), YUNI(J), ZMP(KK)
      END DO
    END DO
    WRITE (104, *) 'POINT_DATA ', INUMUNI*JNUMUNI
    WRITE (104, '(A)') 'SCALARS uavg double 1'
    WRITE (104, '(A)') 'LOOKUP_TABLE default'
    DO J = 1, N2M, 2
      DO I = 1, N1M, 2
        WRITE (104, *) UXX(I, J)
      END DO
    END DO
    WRITE (104, '(A)') 'SCALARS vavg double 1'
    WRITE (104, '(A)') 'LOOKUP_TABLE default'
    DO J = 1, N2M, 2
      DO I = 1, N1M, 2
        WRITE (104, *) VXX(I, J)
      END DO
    END DO
    CLOSE (104)
  END IF

  IF (IOUTFMT .EQ. 2) CALL OUTPUTK_TEC(KPOINT)

  RETURN
END

!=======================================================================
SUBROUTINE OUTPUT_3D
!=======================================================================
!
!     Write a 3D output field file
!
!=======================================================================
  USE MOD_POSTAVG
  IMPLICIT NONE
  INTEGER*8    :: I, J, K, L
  INTEGER*8    :: INUM, JNUM, KNUM

  INUM = 0
  JNUM = 0
  KNUM = 0
  DO I = ISTART, IEND, ISKIP
    INUM = INUM + 1
  END DO
  DO J = JSTART, JEND, JSKIP
    JNUM = JNUM + 1
  END DO
  DO K = KSTART, KEND, KSKIP
    KNUM = KNUM + 1
  END DO

  IF (IOUTFMT .EQ. 1) THEN
    CALL OUTPUT_3D_TEC(INUM, JNUM, KNUM)
    RETURN
  END IF

  IF (IOUTFMT .EQ. 2) CALL OUTPUT_3D_TEC(INUM, JNUM, KNUM)

!====== (When ILD2 = 1) Write scalar criterion-only field
  IF (ILD2 .EQ. 1) THEN
    OPEN (100, FILE='../output/post_avg/'//tname2//'_ld2.vtk')
    WRITE (100, '(A)') '# vtk DataFile Version 3.0'
    WRITE (100, '(A)') 'LESwHT averaged 3D ld2-like (from ssavg)'
    WRITE (100, '(A)') 'ASCII'
    WRITE (100, '(A)') 'DATASET STRUCTURED_GRID'
    WRITE (100, *) 'DIMENSIONS ', INUM, JNUM, KNUM
    WRITE (100, *) 'POINTS ', INUM*JNUM*KNUM, ' double'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (100, *) XMP(I), YMP(J), ZMP(K)
        END DO
      END DO
    END DO
    WRITE (100, *) 'POINT_DATA ', INUM*JNUM*KNUM
    WRITE (100, '(A)') 'SCALARS lambda2_like_ssavg double 1'
    WRITE (100, '(A)') 'LOOKUP_TABLE default'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (100, *) SSAVG(I, J, K)
        END DO
      END DO
    END DO
    CLOSE (100)
  END IF

!====== (When IUVWP = 1) Write UVWP values only
  IF (IUVWP .EQ. 1) THEN
    OPEN (101, FILE='../output/post_avg/'//tname2//'_uvwp.vtk')
    WRITE (101, '(A)') '# vtk DataFile Version 3.0'
    WRITE (101, '(A)') 'LESwHT averaged 3D uvwp'
    WRITE (101, '(A)') 'ASCII'
    WRITE (101, '(A)') 'DATASET STRUCTURED_GRID'
    WRITE (101, *) 'DIMENSIONS ', INUM, JNUM, KNUM
    WRITE (101, *) 'POINTS ', INUM*JNUM*KNUM, ' double'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (101, *) XMP(I), YMP(J), ZMP(K)
        END DO
      END DO
    END DO
    WRITE (101, *) 'POINT_DATA ', INUM*JNUM*KNUM
    WRITE (101, '(A)') 'VECTORS mean_velocity double'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (101, *) UAVG(I, J, K), VAVG(I, J, K), WAVG(I, J, K)
        END DO
      END DO
    END DO
    WRITE (101, '(A)') 'SCALARS pavg double 1'
    WRITE (101, '(A)') 'LOOKUP_TABLE default'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (101, *) PAVG(I, J, K)
        END DO
      END DO
    END DO
    WRITE (101, '(A)') 'SCALARS p2avg double 1'
    WRITE (101, '(A)') 'LOOKUP_TABLE default'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (101, *) P2AVG(I, J, K)
        END DO
      END DO
    END DO
    DO L = 1, 6
      WRITE (101, '(A,I1,A)') 'SCALARS uiuj_', L, ' double 1'
      WRITE (101, '(A)') 'LOOKUP_TABLE default'
      DO K = KSTART, KEND, KSKIP
        DO J = JSTART, JEND, JSKIP
          DO I = ISTART, IEND, ISKIP
            WRITE (101, *) UIUJAVG(I, J, K, L)
          END DO
        END DO
      END DO
    END DO
    CLOSE (101)
  END IF
!======

!====== (When IWXYZ = 1) Write WXYZ values only
  IF (IWXYZ .EQ. 1) THEN
    OPEN (102, FILE='../output/post_avg/'//tname2//'_wxyz.vtk')
    WRITE (102, '(A)') '# vtk DataFile Version 3.0'
    WRITE (102, '(A)') 'LESwHT averaged 3D wxyz'
    WRITE (102, '(A)') 'ASCII'
    WRITE (102, '(A)') 'DATASET STRUCTURED_GRID'
    WRITE (102, *) 'DIMENSIONS ', INUM, JNUM, KNUM
    WRITE (102, *) 'POINTS ', INUM*JNUM*KNUM, ' double'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (102, *) XMP(I), YMP(J), ZMP(K)
        END DO
      END DO
    END DO
    WRITE (102, *) 'POINT_DATA ', INUM*JNUM*KNUM
    WRITE (102, '(A)') 'VECTORS mean_vorticity double'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (102, *) (VORAVG(I, J, K, L), L=1, 3)
        END DO
      END DO
    END DO
    DO L = 1, 6
      WRITE (102, '(A,I1,A)') 'SCALARS vor2_', L, ' double 1'
      WRITE (102, '(A)') 'LOOKUP_TABLE default'
      DO K = KSTART, KEND, KSKIP
        DO J = JSTART, JEND, JSKIP
          DO I = ISTART, IEND, ISKIP
            WRITE (102, *) VOR2AVG(I, J, K, L)
          END DO
        END DO
      END DO
    END DO
    WRITE (102, '(A)') 'SCALARS ssavg double 1'
    WRITE (102, '(A)') 'LOOKUP_TABLE default'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (102, *) SSAVG(I, J, K)
        END DO
      END DO
    END DO
    CLOSE (102)
  END IF
!======

!====== (When IALL = 1) Write all
  IF (IALL .EQ. 1) THEN
    OPEN (103, FILE='../output/post_avg/'//tname2//'_all.vtk')
    WRITE (103, '(A)') '# vtk DataFile Version 3.0'
    WRITE (103, '(A)') 'LESwHT averaged 3D all'
    WRITE (103, '(A)') 'ASCII'
    WRITE (103, '(A)') 'DATASET STRUCTURED_GRID'
    WRITE (103, *) 'DIMENSIONS ', INUM, JNUM, KNUM
    WRITE (103, *) 'POINTS ', INUM*JNUM*KNUM, ' double'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (103, *) XMP(I), YMP(J), ZMP(K)
        END DO
      END DO
    END DO
    WRITE (103, *) 'POINT_DATA ', INUM*JNUM*KNUM
    WRITE (103, '(A)') 'VECTORS mean_velocity double'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (103, *) UAVG(I, J, K), VAVG(I, J, K), WAVG(I, J, K)
        END DO
      END DO
    END DO
    WRITE (103, '(A)') 'SCALARS pavg double 1'
    WRITE (103, '(A)') 'LOOKUP_TABLE default'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (103, *) PAVG(I, J, K)
        END DO
      END DO
    END DO
    WRITE (103, '(A)') 'SCALARS p2avg double 1'
    WRITE (103, '(A)') 'LOOKUP_TABLE default'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (103, *) P2AVG(I, J, K)
        END DO
      END DO
    END DO
    WRITE (103, '(A)') 'VECTORS mean_vorticity double'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (103, *) (VORAVG(I, J, K, L), L=1, 3)
        END DO
      END DO
    END DO
    DO L = 1, 6
      WRITE (103, '(A,I1,A)') 'SCALARS uiuj_', L, ' double 1'
      WRITE (103, '(A)') 'LOOKUP_TABLE default'
      DO K = KSTART, KEND, KSKIP
        DO J = JSTART, JEND, JSKIP
          DO I = ISTART, IEND, ISKIP
            WRITE (103, *) UIUJAVG(I, J, K, L)
          END DO
        END DO
      END DO
    END DO
    DO L = 1, 6
      WRITE (103, '(A,I1,A)') 'SCALARS vor2_', L, ' double 1'
      WRITE (103, '(A)') 'LOOKUP_TABLE default'
      DO K = KSTART, KEND, KSKIP
        DO J = JSTART, JEND, JSKIP
          DO I = ISTART, IEND, ISKIP
            WRITE (103, *) VOR2AVG(I, J, K, L)
          END DO
        END DO
      END DO
    END DO
    WRITE (103, '(A)') 'SCALARS ssavg double 1'
    WRITE (103, '(A)') 'LOOKUP_TABLE default'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (103, *) SSAVG(I, J, K)
        END DO
      END DO
    END DO
    IF (IHTRANS .EQ. 1) THEN
      WRITE (103, '(A)') 'SCALARS tavg double 1'
      WRITE (103, '(A)') 'LOOKUP_TABLE default'
      DO K = KSTART, KEND, KSKIP
        DO J = JSTART, JEND, JSKIP
          DO I = ISTART, IEND, ISKIP
            WRITE (103, *) TAVG(I, J, K)
          END DO
        END DO
      END DO
      WRITE (103, '(A)') 'SCALARS t2avg double 1'
      WRITE (103, '(A)') 'LOOKUP_TABLE default'
      DO K = KSTART, KEND, KSKIP
        DO J = JSTART, JEND, JSKIP
          DO I = ISTART, IEND, ISKIP
            WRITE (103, *) T2AVG(I, J, K)
          END DO
        END DO
      END DO
    END IF

    CLOSE (103)
  END IF
!======

  RETURN
END

!=======================================================================
SUBROUTINE OUTPUTI_TEC(IPOINT)
!=======================================================================
  USE MOD_POSTAVG
  IMPLICIT NONE
  INTEGER*8 :: IPOINT, J, K
  CHARACTER*9 :: FTAILTEC

  FTAILTEC = FTAILIJK(1:5)//'.tec'
  OPEN (202, FILE='../output/post_avg/'//TNAME1//FTAILTEC)
  IF (IHTRANS .EQ. 1) THEN
                        WRITE(202,'(A)') 'VARIABLES = "X","Y","Z","UAVG","VAVG","WAVG","PAVG","P2AVG","VORX","VORY","VORZ","UXUX","UXUY","UXUZ","UYUY","UYUZ","UZUZ","WXWX","WXWY","WXWZ","WYWY","WYWZ","WZWZ","SSAVG","TAVG","T2AVG"'
  ELSE
                        WRITE(202,'(A)') 'VARIABLES = "X","Y","Z","UAVG","VAVG","WAVG","PAVG","P2AVG","VORX","VORY","VORZ","UXUX","UXUY","UXUZ","UYUY","UYUZ","UZUZ","WXWX","WXWY","WXWZ","WYWY","WYWZ","WZWZ","SSAVG"'
  END IF
  WRITE (202, '(A,I0,A,I0,A)') 'ZONE I=', N3M, ', J=', N2M, ', F=POINT'
  DO J = 1, N2M
    DO K = 1, N3M
      IF (IHTRANS .EQ. 1) THEN
     WRITE(202,*) XMP(IPOINT),YMP(J),ZMP(K),UAVG(IPOINT,J,K),VAVG(IPOINT,J,K),WAVG(IPOINT,J,K),PAVG(IPOINT,J,K),P2AVG(IPOINT,J,K), &
 VORAVG(IPOINT,J,K,1),VORAVG(IPOINT,J,K,2),VORAVG(IPOINT,J,K,3),UIUJAVG(IPOINT,J,K,1),UIUJAVG(IPOINT,J,K,2),UIUJAVG(IPOINT,J,K,3), &
                                                                   UIUJAVG(IPOINT,J,K,4),UIUJAVG(IPOINT,J,K,5),UIUJAVG(IPOINT,J,K,6),VOR2AVG(IPOINT,J,K,1),VOR2AVG(IPOINT,J,K,2),VOR2AVG(IPOINT,J,K,3), &
          VOR2AVG(IPOINT,J,K,4),VOR2AVG(IPOINT,J,K,5),VOR2AVG(IPOINT,J,K,6),SSAVG(IPOINT,J,K),TAVG(IPOINT,J,K),T2AVG(IPOINT,J,K)
      ELSE
     WRITE(202,*) XMP(IPOINT),YMP(J),ZMP(K),UAVG(IPOINT,J,K),VAVG(IPOINT,J,K),WAVG(IPOINT,J,K),PAVG(IPOINT,J,K),P2AVG(IPOINT,J,K), &
 VORAVG(IPOINT,J,K,1),VORAVG(IPOINT,J,K,2),VORAVG(IPOINT,J,K,3),UIUJAVG(IPOINT,J,K,1),UIUJAVG(IPOINT,J,K,2),UIUJAVG(IPOINT,J,K,3), &
                                                                   UIUJAVG(IPOINT,J,K,4),UIUJAVG(IPOINT,J,K,5),UIUJAVG(IPOINT,J,K,6),VOR2AVG(IPOINT,J,K,1),VOR2AVG(IPOINT,J,K,2),VOR2AVG(IPOINT,J,K,3), &
          VOR2AVG(IPOINT, J, K, 4), VOR2AVG(IPOINT, J, K, 5), VOR2AVG(IPOINT, J, K, 6), SSAVG(IPOINT, J, K)
      END IF
    END DO
  END DO
  CLOSE (202)
  RETURN
END

!=======================================================================
SUBROUTINE OUTPUTJ_TEC(JPOINT)
!=======================================================================
  USE MOD_POSTAVG
  IMPLICIT NONE
  INTEGER*8 :: JPOINT, I, K
  CHARACTER*9 :: FTAILTEC

  FTAILTEC = FTAILIJK(1:5)//'.tec'
  OPEN (202, FILE='../output/post_avg/'//TNAME1//FTAILTEC)
  IF (IHTRANS .EQ. 1) THEN
                        WRITE(202,'(A)') 'VARIABLES = "X","Y","Z","UAVG","VAVG","WAVG","PAVG","P2AVG","VORX","VORY","VORZ","UXUX","UXUY","UXUZ","UYUY","UYUZ","UZUZ","WXWX","WXWY","WXWZ","WYWY","WYWZ","WZWZ","SSAVG","TAVG","T2AVG"'
  ELSE
                        WRITE(202,'(A)') 'VARIABLES = "X","Y","Z","UAVG","VAVG","WAVG","PAVG","P2AVG","VORX","VORY","VORZ","UXUX","UXUY","UXUZ","UYUY","UYUZ","UZUZ","WXWX","WXWY","WXWZ","WYWY","WYWZ","WZWZ","SSAVG"'
  END IF
  WRITE (202, '(A,I0,A,I0,A)') 'ZONE I=', N1M, ', J=', N3M, ', F=POINT'
  DO K = 1, N3M
    DO I = 1, N1M
      IF (IHTRANS .EQ. 1) THEN
     WRITE(202,*) XMP(I),YMP(JPOINT),ZMP(K),UAVG(I,JPOINT,K),VAVG(I,JPOINT,K),WAVG(I,JPOINT,K),PAVG(I,JPOINT,K),P2AVG(I,JPOINT,K), &
 VORAVG(I,JPOINT,K,1),VORAVG(I,JPOINT,K,2),VORAVG(I,JPOINT,K,3),UIUJAVG(I,JPOINT,K,1),UIUJAVG(I,JPOINT,K,2),UIUJAVG(I,JPOINT,K,3), &
                                                                   UIUJAVG(I,JPOINT,K,4),UIUJAVG(I,JPOINT,K,5),UIUJAVG(I,JPOINT,K,6),VOR2AVG(I,JPOINT,K,1),VOR2AVG(I,JPOINT,K,2),VOR2AVG(I,JPOINT,K,3), &
          VOR2AVG(I,JPOINT,K,4),VOR2AVG(I,JPOINT,K,5),VOR2AVG(I,JPOINT,K,6),SSAVG(I,JPOINT,K),TAVG(I,JPOINT,K),T2AVG(I,JPOINT,K)
      ELSE
     WRITE(202,*) XMP(I),YMP(JPOINT),ZMP(K),UAVG(I,JPOINT,K),VAVG(I,JPOINT,K),WAVG(I,JPOINT,K),PAVG(I,JPOINT,K),P2AVG(I,JPOINT,K), &
 VORAVG(I,JPOINT,K,1),VORAVG(I,JPOINT,K,2),VORAVG(I,JPOINT,K,3),UIUJAVG(I,JPOINT,K,1),UIUJAVG(I,JPOINT,K,2),UIUJAVG(I,JPOINT,K,3), &
                                                                   UIUJAVG(I,JPOINT,K,4),UIUJAVG(I,JPOINT,K,5),UIUJAVG(I,JPOINT,K,6),VOR2AVG(I,JPOINT,K,1),VOR2AVG(I,JPOINT,K,2),VOR2AVG(I,JPOINT,K,3), &
          VOR2AVG(I, JPOINT, K, 4), VOR2AVG(I, JPOINT, K, 5), VOR2AVG(I, JPOINT, K, 6), SSAVG(I, JPOINT, K)
      END IF
    END DO
  END DO
  CLOSE (202)
  RETURN
END

!=======================================================================
SUBROUTINE OUTPUTK_TEC(KPOINT)
!=======================================================================
  USE MOD_POSTAVG
  IMPLICIT NONE
  INTEGER*8 :: KPOINT, I, J
  CHARACTER*9 :: FTAILTEC

  FTAILTEC = FTAILIJK(1:5)//'.tec'
  OPEN (202, FILE='../output/post_avg/'//TNAME1//FTAILTEC)
  IF (IHTRANS .EQ. 1) THEN
                        WRITE(202,'(A)') 'VARIABLES = "X","Y","Z","UAVG","VAVG","WAVG","PAVG","P2AVG","VORX","VORY","VORZ","UXUX","UXUY","UXUZ","UYUY","UYUZ","UZUZ","WXWX","WXWY","WXWZ","WYWY","WYWZ","WZWZ","SSAVG","TAVG","T2AVG"'
  ELSE
                        WRITE(202,'(A)') 'VARIABLES = "X","Y","Z","UAVG","VAVG","WAVG","PAVG","P2AVG","VORX","VORY","VORZ","UXUX","UXUY","UXUZ","UYUY","UYUZ","UZUZ","WXWX","WXWY","WXWZ","WYWY","WYWZ","WZWZ","SSAVG"'
  END IF
  WRITE (202, '(A,I0,A,I0,A)') 'ZONE I=', N1M, ', J=', N2M, ', F=POINT'
  DO J = 1, N2M
    DO I = 1, N1M
      IF (IHTRANS .EQ. 1) THEN
     WRITE(202,*) XMP(I),YMP(J),ZMP(KPOINT),UAVG(I,J,KPOINT),VAVG(I,J,KPOINT),WAVG(I,J,KPOINT),PAVG(I,J,KPOINT),P2AVG(I,J,KPOINT), &
 VORAVG(I,J,KPOINT,1),VORAVG(I,J,KPOINT,2),VORAVG(I,J,KPOINT,3),UIUJAVG(I,J,KPOINT,1),UIUJAVG(I,J,KPOINT,2),UIUJAVG(I,J,KPOINT,3), &
                                                                   UIUJAVG(I,J,KPOINT,4),UIUJAVG(I,J,KPOINT,5),UIUJAVG(I,J,KPOINT,6),VOR2AVG(I,J,KPOINT,1),VOR2AVG(I,J,KPOINT,2),VOR2AVG(I,J,KPOINT,3), &
          VOR2AVG(I,J,KPOINT,4),VOR2AVG(I,J,KPOINT,5),VOR2AVG(I,J,KPOINT,6),SSAVG(I,J,KPOINT),TAVG(I,J,KPOINT),T2AVG(I,J,KPOINT)
      ELSE
     WRITE(202,*) XMP(I),YMP(J),ZMP(KPOINT),UAVG(I,J,KPOINT),VAVG(I,J,KPOINT),WAVG(I,J,KPOINT),PAVG(I,J,KPOINT),P2AVG(I,J,KPOINT), &
 VORAVG(I,J,KPOINT,1),VORAVG(I,J,KPOINT,2),VORAVG(I,J,KPOINT,3),UIUJAVG(I,J,KPOINT,1),UIUJAVG(I,J,KPOINT,2),UIUJAVG(I,J,KPOINT,3), &
                                                                   UIUJAVG(I,J,KPOINT,4),UIUJAVG(I,J,KPOINT,5),UIUJAVG(I,J,KPOINT,6),VOR2AVG(I,J,KPOINT,1),VOR2AVG(I,J,KPOINT,2),VOR2AVG(I,J,KPOINT,3), &
          VOR2AVG(I, J, KPOINT, 4), VOR2AVG(I, J, KPOINT, 5), VOR2AVG(I, J, KPOINT, 6), SSAVG(I, J, KPOINT)
      END IF
    END DO
  END DO
  CLOSE (202)
  RETURN
END

!=======================================================================
SUBROUTINE OUTPUT_3D_TEC(INUM, JNUM, KNUM)
!=======================================================================
  USE MOD_POSTAVG
  IMPLICIT NONE
  INTEGER*8 :: INUM, JNUM, KNUM, I, J, K

  IF (ILD2 .EQ. 1) THEN
    OPEN (210, FILE='../output/post_avg/'//TNAME2//'_ld2.tec')
    WRITE (210, '(A)') 'VARIABLES = "X","Y","Z","LAMBDA2_LIKE_SSAVG"'
    WRITE (210, '(A,I0,A,I0,A,I0,A)') 'ZONE I=', INUM, ', J=', JNUM, ', K=', KNUM, ', F=POINT'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (210, *) XMP(I), YMP(J), ZMP(K), SSAVG(I, J, K)
        END DO
      END DO
    END DO
    CLOSE (210)
  END IF

  IF (IUVWP .EQ. 1) THEN
    OPEN (211, FILE='../output/post_avg/'//TNAME2//'_uvwp.tec')
    IF (IHTRANS .EQ. 1) THEN
                              WRITE(211,'(A)') 'VARIABLES = "X","Y","Z","UAVG","VAVG","WAVG","PAVG","P2AVG","UXUX","UXUY","UXUZ","UYUY","UYUZ","UZUZ","TAVG","T2AVG"'
    ELSE
      WRITE (211, '(A)') 'VARIABLES = "X","Y","Z","UAVG","VAVG","WAVG","PAVG","P2AVG","UXUX","UXUY","UXUZ","UYUY","UYUZ","UZUZ"'
    END IF
    WRITE (211, '(A,I0,A,I0,A,I0,A)') 'ZONE I=', INUM, ', J=', JNUM, ', K=', KNUM, ', F=POINT'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          IF (IHTRANS .EQ. 1) THEN
            WRITE (211, *) XMP(I), YMP(J), ZMP(K), UAVG(I, J, K), VAVG(I, J, K), WAVG(I, J, K), PAVG(I, J, K), P2AVG(I, J, K), &
      UIUJAVG(I,J,K,1),UIUJAVG(I,J,K,2),UIUJAVG(I,J,K,3),UIUJAVG(I,J,K,4),UIUJAVG(I,J,K,5),UIUJAVG(I,J,K,6),TAVG(I,J,K),T2AVG(I,J,K)
          ELSE
            WRITE (211, *) XMP(I), YMP(J), ZMP(K), UAVG(I, J, K), VAVG(I, J, K), WAVG(I, J, K), PAVG(I, J, K), P2AVG(I, J, K), &
        UIUJAVG(I, J, K, 1), UIUJAVG(I, J, K, 2), UIUJAVG(I, J, K, 3), UIUJAVG(I, J, K, 4), UIUJAVG(I, J, K, 5), UIUJAVG(I, J, K, 6)
          END IF
        END DO
      END DO
    END DO
    CLOSE (211)
  END IF

  IF (IWXYZ .EQ. 1) THEN
    OPEN (212, FILE='../output/post_avg/'//TNAME2//'_wxyz.tec')
    WRITE (212, '(A)') 'VARIABLES = "X","Y","Z","VORX","VORY","VORZ","WXWX","WXWY","WXWZ","WYWY","WYWZ","WZWZ","SSAVG"'
    WRITE (212, '(A,I0,A,I0,A,I0,A)') 'ZONE I=', INUM, ', J=', JNUM, ', K=', KNUM, ', F=POINT'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (212, *) XMP(I), YMP(J), ZMP(K), VORAVG(I, J, K, 1), VORAVG(I, J, K, 2), VORAVG(I, J, K, 3), &
            VOR2AVG(I,J,K,1),VOR2AVG(I,J,K,2),VOR2AVG(I,J,K,3),VOR2AVG(I,J,K,4),VOR2AVG(I,J,K,5),VOR2AVG(I,J,K,6),SSAVG(I,J,K)
        END DO
      END DO
    END DO
    CLOSE (212)
  END IF

  IF (IALL .EQ. 1) THEN
    OPEN (213, FILE='../output/post_avg/'//TNAME2//'_all.tec')
    IF (IHTRANS .EQ. 1) THEN
                              WRITE(213,'(A)') 'VARIABLES = "X","Y","Z","UAVG","VAVG","WAVG","PAVG","P2AVG","VORX","VORY","VORZ","UXUX","UXUY","UXUZ","UYUY","UYUZ","UZUZ","WXWX","WXWY","WXWZ","WYWY","WYWZ","WZWZ","SSAVG","TAVG","T2AVG"'
    ELSE
                              WRITE(213,'(A)') 'VARIABLES = "X","Y","Z","UAVG","VAVG","WAVG","PAVG","P2AVG","VORX","VORY","VORZ","UXUX","UXUY","UXUZ","UYUY","UYUZ","UZUZ","WXWX","WXWY","WXWZ","WYWY","WYWZ","WZWZ","SSAVG"'
    END IF
    WRITE (213, '(A,I0,A,I0,A,I0,A)') 'ZONE I=', INUM, ', J=', JNUM, ', K=', KNUM, ', F=POINT'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          IF (IHTRANS .EQ. 1) THEN
            WRITE (213, *) XMP(I), YMP(J), ZMP(K), UAVG(I, J, K), VAVG(I, J, K), WAVG(I, J, K), PAVG(I, J, K), P2AVG(I, J, K), &
        VORAVG(I, J, K, 1), VORAVG(I, J, K, 2), VORAVG(I, J, K, 3), UIUJAVG(I, J, K, 1), UIUJAVG(I, J, K, 2), UIUJAVG(I, J, K, 3), &
     UIUJAVG(I, J, K, 4), UIUJAVG(I, J, K, 5), UIUJAVG(I, J, K, 6), VOR2AVG(I, J, K, 1), VOR2AVG(I, J, K, 2), VOR2AVG(I, J, K, 3), &
              VOR2AVG(I, J, K, 4), VOR2AVG(I, J, K, 5), VOR2AVG(I, J, K, 6), SSAVG(I, J, K), TAVG(I, J, K), T2AVG(I, J, K)
          ELSE
            WRITE (213, *) XMP(I), YMP(J), ZMP(K), UAVG(I, J, K), VAVG(I, J, K), WAVG(I, J, K), PAVG(I, J, K), P2AVG(I, J, K), &
        VORAVG(I, J, K, 1), VORAVG(I, J, K, 2), VORAVG(I, J, K, 3), UIUJAVG(I, J, K, 1), UIUJAVG(I, J, K, 2), UIUJAVG(I, J, K, 3), &
     UIUJAVG(I, J, K, 4), UIUJAVG(I, J, K, 5), UIUJAVG(I, J, K, 6), VOR2AVG(I, J, K, 1), VOR2AVG(I, J, K, 2), VOR2AVG(I, J, K, 3), &
              VOR2AVG(I, J, K, 4), VOR2AVG(I, J, K, 5), VOR2AVG(I, J, K, 6), SSAVG(I, J, K)
          END IF
        END DO
      END DO
    END DO
    CLOSE (213)
  END IF

  RETURN
END

!=======================================================================
SUBROUTINE PREFLD(fileprevel)
!=======================================================================
!
!     Read an input field file
!
!=======================================================================
  USE MOD_POSTAVG
  IMPLICIT NONE
  CHARACTER*20 :: fileprevel
  INTEGER*8    :: I, J, K, L
  INTEGER*8    :: NN1, NN2, NN3
  REAL*8       :: RRE
  INTEGER*8    :: IOS
  CHARACTER*256 :: FPATH

  UAVG = 0.
  VAVG = 0.
  WAVG = 0.
  PAVG = 0.
  UIUJAVG = 0.
  P2AVG = 0.
  VORAVG = 0.
  VOR2AVG = 0.
  SSAVG = 0.

  FPATH = '../output/field_avg/'//TRIM(fileprevel)
  OPEN (12, FILE=TRIM(FPATH), FORM='UNFORMATTED', STATUS='OLD', ACTION='READ', IOSTAT=IOS)
  IF (IOS .NE. 0) THEN
    WRITE (*, *) 'ERROR: CANNOT OPEN AVERAGE FIELD FILE: ', TRIM(FPATH)
    WRITE (*, *) 'IOSTAT=', IOS
    STOP
  END IF

!     dum for future use
  READ (12, IOSTAT=IOS) NN1, NN2, NN3, RRE
  IF (IOS .NE. 0) GOTO 901
  READ (12, IOSTAT=IOS) TIMEINIT, TIMEEND, IHISTINIT, IHISTEND
  IF (IOS .NE. 0) GOTO 902
  READ (12, IOSTAT=IOS) (((UAVG(I, J, K), I=1, N1M), J=1, N2M), K=1, N3M)
  IF (IOS .NE. 0) GOTO 903
  READ (12, IOSTAT=IOS) (((VAVG(I, J, K), I=1, N1M), J=1, N2M), K=1, N3M)
  IF (IOS .NE. 0) GOTO 904
  READ (12, IOSTAT=IOS) (((WAVG(I, J, K), I=1, N1M), J=1, N2M), K=1, N3M)
  IF (IOS .NE. 0) GOTO 905
  READ (12, IOSTAT=IOS) ((((UIUJAVG(I, J, K, L), I=1, N1M), J=1, N2M), K=1, N3M), L=1, 6)
  IF (IOS .NE. 0) GOTO 906
  READ (12, IOSTAT=IOS) (((PAVG(I, J, K), I=1, N1M), J=1, N2M), K=1, N3M)
  IF (IOS .NE. 0) GOTO 907
  READ (12, IOSTAT=IOS) (((P2AVG(I, J, K), I=1, N1M), J=1, N2M), K=1, N3M)
  IF (IOS .NE. 0) GOTO 908
  READ (12, IOSTAT=IOS) ((((VORAVG(I, J, K, L), I=1, N1M), J=1, N2M), K=1, N3M), L=1, 3)
  IF (IOS .NE. 0) GOTO 909
  READ (12, IOSTAT=IOS) ((((VOR2AVG(I, J, K, L), I=1, N1M), J=1, N2M), K=1, N3M), L=1, 6)
  IF (IOS .NE. 0) GOTO 910
  READ (12, IOSTAT=IOS) (((SSAVG(I, J, K), I=1, N1M), J=1, N2M), K=1, N3M)
  IF (IOS .NE. 0) GOTO 911
  IF (IHTRANS .EQ. 1) THEN
    READ (12, IOSTAT=IOS) (((TAVG(I, J, K), I=1, N1M), J=1, N2M), K=1, N3M)
    IF (IOS .NE. 0) GOTO 912
    READ (12, IOSTAT=IOS) (((T2AVG(I, J, K), I=1, N1M), J=1, N2M), K=1, N3M)
    IF (IOS .NE. 0) GOTO 913
  END IF

  CLOSE (12)

  IF ((TIMEEND - TIMEINIT) .LE. 0.D0) THEN
    WRITE (*, *) 'ERROR: INVALID AVERAGING WINDOW IN ', TRIM(FPATH)
    WRITE (*, *) 'TIMEINIT=', TIMEINIT, ' TIMEEND=', TIMEEND
    STOP
  END IF

  UAVG = UAVG/(TIMEEND - TIMEINIT)
  VAVG = VAVG/(TIMEEND - TIMEINIT)
  WAVG = WAVG/(TIMEEND - TIMEINIT)
  PAVG = PAVG/(TIMEEND - TIMEINIT)
  UIUJAVG = UIUJAVG/(TIMEEND - TIMEINIT)
  P2AVG = P2AVG/(TIMEEND - TIMEINIT)
  VORAVG = VORAVG/(TIMEEND - TIMEINIT)
  VOR2AVG = VOR2AVG/(TIMEEND - TIMEINIT)
  SSAVG = SSAVG/(TIMEEND - TIMEINIT)
  IF (IHTRANS .EQ. 1) THEN
    TAVG = TAVG/(TIMEEND - TIMEINIT)
    T2AVG = T2AVG/(TIMEEND - TIMEINIT)
  END IF

  WRITE (*, 100)
  WRITE (*, 101)
  WRITE (*, 102) fileprevel
  WRITE (*, 103) RRE
  WRITE (*, 104) NN1, NN2, NN3
  WRITE (*, 105) TIMEINIT, TIMEEND, IHISTINIT, IHISTEND
  WRITE (*, *) UAVG(1, 1, 1)

100 FORMAT('----------- INITIAL FIELD INFORMATION -----------')
101 FORMAT('INITIAL FIELD      : READING DONE')
102 FORMAT('INITIAL FIELD NAME : ', A30)
103 FORMAT('RE=', ES12.3)
104 FORMAT('N1=', I10, '  N2=', I10, '  N3=', I10)
105 FORMAT('TI=', F12.5, '  TF=', F12.5, '  HISTI= ', I10, '  HISTF=', I10)

  RETURN
901 WRITE (*, *) 'ERROR READING HEADER-1 FROM ', TRIM(FPATH), ' IOSTAT=', IOS
  STOP
902 WRITE (*, *) 'ERROR READING HEADER-2 FROM ', TRIM(FPATH), ' IOSTAT=', IOS
  STOP
903 WRITE (*, *) 'ERROR READING UAVG FROM ', TRIM(FPATH), ' IOSTAT=', IOS
  STOP
904 WRITE (*, *) 'ERROR READING VAVG FROM ', TRIM(FPATH), ' IOSTAT=', IOS
  STOP
905 WRITE (*, *) 'ERROR READING WAVG FROM ', TRIM(FPATH), ' IOSTAT=', IOS
  STOP
906 WRITE (*, *) 'ERROR READING UIUJAVG FROM ', TRIM(FPATH), ' IOSTAT=', IOS
  STOP
907 WRITE (*, *) 'ERROR READING PAVG FROM ', TRIM(FPATH), ' IOSTAT=', IOS
  STOP
908 WRITE (*, *) 'ERROR READING P2AVG FROM ', TRIM(FPATH), ' IOSTAT=', IOS
  STOP
909 WRITE (*, *) 'ERROR READING VORAVG FROM ', TRIM(FPATH), ' IOSTAT=', IOS
  STOP
910 WRITE (*, *) 'ERROR READING VOR2AVG FROM ', TRIM(FPATH), ' IOSTAT=', IOS
  STOP
911 WRITE (*, *) 'ERROR READING SSAVG FROM ', TRIM(FPATH), ' IOSTAT=', IOS
  STOP
912 WRITE (*, *) 'ERROR READING TAVG FROM ', TRIM(FPATH), ' IOSTAT=', IOS
  STOP
913 WRITE (*, *) 'ERROR READING T2AVG FROM ', TRIM(FPATH), ' IOSTAT=', IOS
  STOP
END

!=======================================================================
!     Section for Geometry
!=======================================================================
SUBROUTINE GEOM(gridfile)
!=======================================================================
!
!     Read 'grid.dat' and define various variables related to grid
!
!     GEOM
!       K_BGPZ,I_BGPX: intial index in x-,z-directions
!     READGRID
!       N1,N2,N3   : number of grid in x-,y-,z-directions
!       N1M,N2M,N3M: values obtained by subtracting 1 from N1,N2,N3
!       XL,YL,ZL   : domain size
!       X,Y,Z      : cell face positions in each direction
!     INDICES
!       IPV,JPV,KPV: I,J,K => X,Y,Z and P => +1, IPV(I)=I+1
!       IMV,JMV,KMV:                    M => -1, IMV(I)=I-1
!     INDXFIX
!       FIXIL,FIXJL,FIXKL: lower boundary indicators, which are only 1
!                          at the lower (I,J,K=1) boundary
!       FIXIU,FIXJU,FIXKU: upper boundary indicators, which are only 1
!                          at the upper (I,J,K=Max) boundary
!     MESHES
!       SDX,SDY,SDZ: distance between cell face,   SDX(I)=X(I+1)-X(I)
!       VDX,VDY,VDZ: distance between cell center, VDX(I)=(SDX(I)+SDX(I-1))/2
!       SSDX,SSDY,SSDZ: reciprocal of SDX,SDY,SDZ for spatial discretization
!       VVDX,VVDY,VVDZ: reciprocal of VDX,VDY,VDZ for spatial discretization
!     PHYSPOS
!       XMP,YMP,ZMP: cell center positions,        XMP(I)=X(I)+SDX(I)/2
!
!-----------------------------------------------------------------------
  USE MOD_POSTAVG
  IMPLICIT NONE
  INTEGER*8    :: I, K
  CHARACTER*23  :: gridfile

  OPEN (11, FILE=gridfile)
  CALL READGRID
  CLOSE (11)

  CALL INDICES
  CALL INDXFIX
  CALL MESHES
  CALL PHYSPOS

  IF (IPZ .EQ. 1) THEN
    FIXKL(1) = 0.
    FIXKU(N3M) = 0.
    SDZ(0) = SDZ(N3M)
    SDZ(N3) = SDZ(1)
    SSDZ(0) = SSDZ(N3M)
    SSDZ(N3) = SSDZ(1)
    VDZ(1) = 0.5*(SDZ(1) + SDZ(N3M))
    VDZ(N3) = 0.5*(SDZ(1) + SDZ(N3M))
    VVDZ(1) = 1./VDZ(1)
    VVDZ(N3) = 1./VDZ(N3)
    KMV(1) = N3M
    KPV(N3M) = 1
  END IF

  IF (IPX .EQ. 1) THEN
    FIXIL(1) = 0.
    FIXIU(N1M) = 0.
    SDX(0) = SDX(N1M)
    SDX(N1) = SDX(1)
    SSDX(0) = SSDX(N1M)
    SSDX(N1) = SSDX(1)
    VDX(1) = 0.5*(SDX(1) + SDX(N1M))
    VDX(N1) = 0.5*(SDX(1) + SDX(N1M))
    VVDX(1) = 1./VDX(1)
    VVDX(N1) = 1./VDX(N1)
    IMV(1) = N1M
    IPV(N1M) = 1
  END IF

  RETURN
END
!=======================================================================
SUBROUTINE READGRID
!=======================================================================
!
!     Read grid information from 'grid.dat'
!
!=======================================================================
  USE MOD_POSTAVG
  IMPLICIT NONE
  INTEGER*8    :: I, J, K

  READ (11, *) N1, N2, N3
  WRITE (*, 101) N1, N2, N3
  READ (11, *) XL, YL, ZL
  WRITE (*, 102) XL, YL, ZL
101 FORMAT('NX=', I12, '  NY=', I12, '  NZ=', I12)
102 FORMAT('XL=', F12.4, '  YL=', F12.4, '  ZL=', F12.4)

  N1M = N1 - 1
  N2M = N2 - 1
  N3M = N3 - 1

  IF ((N1 .GT. M1) .OR. (N2 .GT. M2) .OR. (N3 .GT. M3)) THEN
    PRINT *, 'ARRAY SIZE CAN NOT HANDLE THIS GRID.'
    STOP
  END IF

  READ (11, *) (X(I), I=1, N1)
  READ (11, *) (Y(J), J=1, N2)
  READ (11, *) (Z(K), K=1, N3)

  RETURN
END
!=======================================================================
SUBROUTINE INDICES
!=======================================================================
!
!     Calculate adjacent grid's indices
!
!=======================================================================
  USE MOD_POSTAVG
  IMPLICIT NONE
  INTEGER*8    :: I, J, K

!-----STREAMWISE DIRECTION
  DO I = 1, N1M
    IPV(I) = I + 1
    IMV(I) = I - 1
  END DO

!-----NORMAL DIRECTION
  DO J = 1, N2M
    JPV(J) = J + 1
    JMV(J) = J - 1
  END DO

!-----SPANWISE DIRECTION
  DO K = 1, N3M
    KPV(K) = K + 1
    KMV(K) = K - 1
  END DO

  RETURN
END

!=======================================================================
SUBROUTINE INDXFIX
!=======================================================================
!
!     Variables FIXxx(ex FIXIU) are 0 if a point is an interior point,
!     and 1 if a point is a boundary point (These variables are used
!     in the routine VORNLAMBDA2)
!
!=======================================================================
  USE MOD_POSTAVG
  IMPLICIT NONE
  INTEGER*8    :: I, J, K

  DO I = 1, N1M
    FIXIL(I) = 0.
    FIXIU(I) = 0.
  END DO
  FIXIL(1) = 1.
  FIXIU(N1M) = 1.

  DO J = 1, N2M
    FIXJL(J) = 0.
    FIXJU(J) = 0.
  END DO
  FIXJL(1) = 1.
  FIXJU(N2M) = 1.

  DO K = 1, N3M
    FIXKL(K) = 0.
    FIXKU(K) = 0.
  END DO
  FIXKL(1) = 1.
  FIXKU(N3M) = 1.

  RETURN
END

!=======================================================================
SUBROUTINE MESHES
!=======================================================================
!
!     Calculate side length (SDX) and cell center length (VDX) of each cell
!
!=======================================================================
  USE MOD_POSTAVG
  IMPLICIT NONE
  INTEGER*8    :: I, J, K

  DO I = 1, N1M
    SDX(I) = X(I + 1) - X(I)
  END DO
  DO I = 2, N1M
    VDX(I) = 0.5*(SDX(I) + SDX(I - 1))
  END DO
  VDX(1) = 0.5*SDX(1)             ! LEFT BOUNDED COND.
  VDX(N1) = 0.5*SDX(N1M)          ! RIGHT BOUNDED COND.

  DO I = 1, N1M
    SSDX(I) = 1./SDX(I)
  END DO
  DO I = 1, N1
    VVDX(I) = 1./VDX(I)
  END DO

  DO J = 1, N2M
    SDY(J) = Y(J + 1) - Y(J)
  END DO
  DO J = 2, N2M
    VDY(J) = 0.5*(SDY(J) + SDY(J - 1))
  END DO
  VDY(1) = 0.5*SDY(1)             ! LOWER WALL BOUNDED COND.
  VDY(N2) = 0.5*SDY(N2M)          ! UPPER WALL BOUNDED COND.

  DO J = 1, N2M
    SSDY(J) = 1./SDY(J)
  END DO
  DO J = 1, N2
    VVDY(J) = 1./VDY(J)
  END DO

  DO K = 1, N3M
    SDZ(K) = Z(K + 1) - Z(K)
  END DO
  DO K = 2, N3M
    VDZ(K) = 0.5*(SDZ(K) + SDZ(K - 1))
  END DO
  VDZ(1) = 0.5*SDZ(1)             ! LOWER WALL BOUNDED COND.
  VDZ(N3) = 0.5*SDZ(N3M)          ! UPPER WALL BOUNDED COND.

  DO K = 1, N3M
    SSDZ(K) = 1./SDZ(K)
  END DO
  DO K = 1, N3
    VVDZ(K) = 1./VDZ(K)
  END DO

!-----set by zero
  SDX(0) = 0.
  SDX(M1) = 0.
  SDY(0) = 0.
  SDY(M2) = 0.
  SDZ(0) = 0.
  SDZ(M3) = 0.

  RETURN
END

!=======================================================================
SUBROUTINE PHYSPOS
!=======================================================================
!
!     Calculate midpoint coordinates of each cell
!
!=======================================================================
  USE MOD_POSTAVG
  IMPLICIT NONE
  INTEGER*8    :: I, J, K

  DO I = 1, N1M
    XMP(I) = X(I) + 0.5*SDX(I)
  END DO
  XMP(N1) = X(N1)
  XMP(0) = X(1)

  DO J = 1, N2M
    YMP(J) = Y(J) + 0.5*SDY(J)
  END DO
  YMP(N2) = Y(N2)
  YMP(0) = Y(1)

  DO K = 1, N3M
    ZMP(K) = Z(K) + 0.5*SDZ(K)
  END DO
  ZMP(N3) = Z(N3)
  ZMP(0) = Z(1)

  RETURN
END

!=======================================================================
SUBROUTINE GRIDUNI
!=======================================================================
  USE MOD_POSTAVG
  IMPLICIT NONE
  INTEGER*8    :: I, J, K
  REAL*8       :: XINIT, XEND, YINIT, YEND, ZINIT, ZEND, UNIXL, UNIYL, UNIZL
  REAL*8       :: DXUNI, DYUNI, DZUNI

  INUMUNI = 0
  DO I = 1, N1M, 2
    INUMUNI = INUMUNI + 1
  END DO
  JNUMUNI = 0
  DO J = 1, N2M, 2
    JNUMUNI = JNUMUNI + 1
  END DO
  KNUMUNI = 0
  DO K = 1, N3M, 2
    KNUMUNI = KNUMUNI + 1
  END DO

  XINIT = X(1)
  XEND = X(N1)
  YINIT = Y(1)
  YEND = Y(N2)
  ZINIT = Z(1)
  ZEND = Z(N3)
  UNIXL = XEND - XINIT
  UNIYL = YEND - YINIT
  UNIZL = ZEND - ZINIT
  DXUNI = UNIXL/N1M
  DYUNI = UNIYL/N2M
  DZUNI = UNIZL/N3M

  WRITE (*, 11) DXUNI
  WRITE (*, 12) DYUNI
  WRITE (*, 13) DZUNI
11 FORMAT('DXUNI = ', F8.5)
12 FORMAT('DYUNI = ', F8.5)
13 FORMAT('DZUNI = ', F8.5)

  XUNI(1) = XINIT + 0.5*DXUNI
  YUNI(1) = YINIT + 0.5*DYUNI
  ZUNI(1) = ZINIT + 0.5*DZUNI

  DO I = 2, N1
    XUNI(I) = XUNI(I - 1) + DXUNI
  END DO
  DO J = 2, N2
    YUNI(J) = YUNI(J - 1) + DYUNI
  END DO
  DO K = 2, N3
    ZUNI(K) = ZUNI(K - 1) + DZUNI
  END DO

  DO I = 1, N1M
    IF (XUNI(I) .LE. XMP(1)) THEN
      IUNI(I) = 1
      FACUNIX(I) = 0.D0
    ELSEIF (XUNI(I) .GE. XMP(N1M)) THEN
      IUNI(I) = N1M - 1
      FACUNIX(I) = 1.D0
    ELSE
      DO J = N1M - 1, 1, -1
        IF (XUNI(I) .GE. XMP(J)) THEN
          IUNI(I) = J
          FACUNIX(I) = (XUNI(I) - XMP(J))/(XMP(J + 1) - XMP(J))
          EXIT
        END IF
      END DO
    END IF
  END DO

  DO I = 1, N2M
    IF (YUNI(I) .LE. YMP(1)) THEN
      JUNI(I) = 1
      FACUNIY(I) = 0.D0
    ELSEIF (YUNI(I) .GE. YMP(N2M)) THEN
      JUNI(I) = N2M - 1
      FACUNIY(I) = 1.D0
    ELSE
      DO J = N2M - 1, 1, -1
        IF (YUNI(I) .GE. YMP(J)) THEN
          JUNI(I) = J
          FACUNIY(I) = (YUNI(I) - YMP(J))/(YMP(J + 1) - YMP(J))
          EXIT
        END IF
      END DO
    END IF
  END DO

  DO I = 1, N3M
    IF (ZUNI(I) .LE. ZMP(1)) THEN
      KUNI(I) = 1
      FACUNIZ(I) = 0.D0
    ELSEIF (ZUNI(I) .GE. ZMP(N3M)) THEN
      KUNI(I) = N3M - 1
      FACUNIZ(I) = 1.D0
    ELSE
      DO J = N3M - 1, 1, -1
        IF (ZUNI(I) .GE. ZMP(J)) THEN
          KUNI(I) = J
          FACUNIZ(I) = (ZUNI(I) - ZMP(J))/(ZMP(J + 1) - ZMP(J))
          EXIT
        END IF
      END DO
    END IF
  END DO

  RETURN
END

!=======================================================================
SUBROUTINE FINDIPJPKP
!=======================================================================
!
!     Find indices which correspond to input points (X_position, Y_position,
!     Z_position of post_inst.in file)
!
!=======================================================================
  USE MOD_POSTAVG
  IMPLICIT NONE
  INTEGER*8    :: I, J, K

  DO I = 1, NIP
    DO J = N1M, 1, -1
      IF (XIP(I) .GE. XMP(J)) THEN
        IP(I) = J
        GOTO 111
      END IF
    END DO
111 CONTINUE
  END DO

  DO I = 1, NJP
    DO J = N2M, 1, -1
      IF (YJP(I) .GE. YMP(J)) THEN
        JP(I) = J
        GOTO 112
      END IF
    END DO
112 CONTINUE
  END DO

  DO I = 1, NKP
    DO J = N3M, 1, -1
      IF (ZKP(I) .GE. ZMP(J)) THEN
        KP(I) = J
        GOTO 113
      END IF
    END DO
113 CONTINUE
  END DO

  RETURN
END
!=======================================================================
!     End of the Section for Geometry
!=======================================================================

!=======================================================================
SUBROUTINE SORTER(N, RA)
!=======================================================================
!
!     Sort the input array RA with size N
!
!=======================================================================
  REAL RA(3), RRA
  INTEGER I, IR, J, L

  IF (N .LE. 1) GOTO 30
  L = N/2 + 1
  IR = N
10 CONTINUE
  IF (L .GT. 1) THEN
    L = L - 1
    RRA = RA(L)
  ELSE
    RRA = RA(IR)
    RA(IR) = RA(1)
    IR = IR - 1
    IF (IR .EQ. 1) THEN
      RA(1) = RRA
      RETURN
    END IF
  END IF
  I = L
  J = L + L
20 IF (J .LE. IR) THEN
    IF (J .LT. IR) THEN
      IF (RA(J) .GT. RA(J + 1)) J = J + 1
    END IF
    IF (RRA .GT. RA(J)) THEN
      RA(I) = RA(J)
      I = J
      J = J + J
    ELSE
      J = IR + 1
    END IF
    GOTO 20
  END IF
  RA(I) = RRA

  GOTO 10
30 RETURN
END

!=======================================================================
SUBROUTINE jacobi(a, n, np, d, v, nrot)
!=======================================================================
!
!     Compute eigenvalues of a matrix (from Numerical recipe in FORTRAN)
!
!=======================================================================
  INTEGER n, np, nrot, NMAX
  REAL a(np, np), d(np), v(np, np)
  PARAMETER(NMAX=500)
  INTEGER i, ip, iq, j
  REAL c, g, h, s, sm, t, tau, theta, tresh, b(NMAX), z(NMAX)
  do ip = 1, n
    do iq = 1, n
      v(ip, iq) = 0.
    end do
    v(ip, ip) = 1.
  end do
  do ip = 1, n
    b(ip) = a(ip, ip)
    d(ip) = b(ip)
    z(ip) = 0.
  end do
  nrot = 0
  do i = 1, 50
    sm = 0.
    do ip = 1, n - 1
      do iq = ip + 1, n
        sm = sm + abs(a(ip, iq))
      end do
    end do
    if (sm .eq. 0.) return
    if (i .lt. 4) then
      tresh = 0.2*sm/n**2
    else
      tresh = 0.
    end if
    do ip = 1, n - 1
      do iq = ip + 1, n
        g = 100.*abs(a(ip, iq))
        if ((i .gt. 4) .and. (abs(d(ip)) + g .eq. abs(d(ip))) .and. (abs(d(iq)) + g .eq. abs(d(iq)))) then
          a(ip, iq) = 0.
        else if (abs(a(ip, iq)) .gt. tresh) then
          h = d(iq) - d(ip)
          if (abs(h) + g .eq. abs(h)) then
            t = a(ip, iq)/h

          else
            theta = 0.5*h/a(ip, iq)
            t = 1./(abs(theta) + sqrt(1.+theta**2))
            if (theta .lt. 0.) t = -t
          end if
          c = 1./sqrt(1 + t**2)
          s = t*c
          tau = s/(1.+c)
          h = t*a(ip, iq)
          z(ip) = z(ip) - h
          z(iq) = z(iq) + h
          d(ip) = d(ip) - h
          d(iq) = d(iq) + h
          a(ip, iq) = 0.
          do j = 1, ip - 1
            g = a(j, ip)
            h = a(j, iq)

            a(j, ip) = g - s*(h + g*tau)
            a(j, iq) = h + s*(g - h*tau)
          end do
          do j = ip + 1, iq - 1
            g = a(ip, j)
            h = a(j, iq)
            a(ip, j) = g - s*(h + g*tau)
            a(j, iq) = h + s*(g - h*tau)
          end do
          do j = iq + 1, n
            g = a(ip, j)
            h = a(iq, j)
            a(ip, j) = g - s*(h + g*tau)
            a(iq, j) = h + s*(g - h*tau)
          end do
          do j = 1, n
            g = v(j, ip)

            h = v(j, iq)
            v(j, ip) = g - s*(h + g*tau)
            v(j, iq) = h + s*(g - h*tau)
          end do
          nrot = nrot + 1
        end if
      end do
    end do
    do ip = 1, n
      b(ip) = b(ip) + z(ip)
      d(ip) = b(ip)
      z(ip) = 0.
    end do
  end do
  pause 'too many iterations in jacobi'
  return
END
!=======================================================================

