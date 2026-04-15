FUNCTION FUNCBODY(X, Y, Z, T) RESULT(VAL)
    IMPLICIT NONE

    REAL(8), INTENT(IN) :: X, Y, Z, T
    REAL(8) :: VAL

    REAL(8), PARAMETER :: EPS = 1.0D-4
    REAL(8), PARAMETER :: EPS_PROJ = 1.0D-12
    REAL(8), PARAMETER :: EPS_BARY = 1.0D-10
    REAL(8), PARAMETER :: EPS_HIT = 1.0D-10
    REAL(8), PARAMETER :: INF = 1.0D99

    LOGICAL, SAVE :: GEOM_READY = .FALSE.
    LOGICAL, SAVE :: GEOM_FAILED = .FALSE.

    INTEGER(8), SAVE :: NTRI = 0_8
    INTEGER, SAVE :: NBX = 0, NBZ = 0, NBINS = 0

    REAL(8), SAVE :: GMIN(3), GMAX(3)
    REAL(8), SAVE :: XMIN, XMAX, ZMIN, ZMAX, DXBIN, DZBIN
    REAL(8), SAVE :: YPAD_BOTTOM

    REAL(8), ALLOCATABLE, SAVE :: V1(:, :), V2(:, :), V3(:, :)
    REAL(8), ALLOCATABLE, SAVE :: BMINXZ(:, :), BMAXXZ(:, :)
    REAL(8), ALLOCATABLE, SAVE :: DEN2D(:)
    LOGICAL, ALLOCATABLE, SAVE :: PROJ_OK(:)

    INTEGER(8), ALLOCATABLE, SAVE :: BIN_START(:), BIN_ITEMS(:)

    LOGICAL :: IS_INSIDE

    IF (Y .GT. 2.D0 - EPS) THEN
        VAL = -1.D0
        RETURN
    END IF
    
    IF (T < -INF) THEN
        VAL = 1.0D0
        RETURN
    END IF

    IF (.NOT. GEOM_READY .AND. .NOT. GEOM_FAILED) THEN
!$OMP CRITICAL(FUNCBODY_GEOM_INIT)
        IF (.NOT. GEOM_READY .AND. .NOT. GEOM_FAILED) THEN
            CALL LOAD_GEOMETRY()
        END IF
!$OMP END CRITICAL(FUNCBODY_GEOM_INIT)
    END IF

    IF (GEOM_FAILED .OR. NTRI <= 0_8) THEN
        VAL = 1.0D0
        RETURN
    END IF

    ! Solid padding below the dimple-depth floor inferred from STL.
    if (Y .LT. YPAD_BOTTOM + EPS) then
        VAL = -1.D0
        RETURN
    end if

    IF (X < GMIN(1) - EPS .OR. X > GMAX(1) + EPS .OR. &
        Y > GMAX(2) + EPS .OR. &
        Z < GMIN(3) - EPS .OR. Z > GMAX(3) + EPS) THEN
        VAL = 1.0D0
        RETURN
    END IF

    ! Vertical two-way rays in wall-normal direction (y):
    ! inside only if we hit the slab wall both upward and downward.
    IS_INSIDE = HAS_BOTH_Y_HITS(X, Y, Z)

    IF (IS_INSIDE) THEN
        VAL = -1.0D0
    ELSE
        VAL = 1.0D0
    END IF

    RETURN

CONTAINS

    SUBROUTINE LOAD_GEOMETRY()
        IMPLICIT NONE

        INTEGER :: IUNIT, IOS, ISTAT, LEN_OUT, ISTAT_Y, LEN_Y
        INTEGER(8) :: NVERT, IVERT, ITRI, IVLOC, I
        CHARACTER(LEN=1024) :: STL_PATH
        CHARACTER(LEN=64) :: YTOP_STR
        CHARACTER(LEN=1024) :: LINE
        CHARACTER(LEN=16) :: TAG
        REAL(8) :: VX, VY, VZ
        REAL(8) :: YTOP_TARGET, YRAW_MIN, YRAW_MAX, YSHIFT
        REAL(8) :: YTOP_SURF_MIN, YTHR
        REAL(8), ALLOCATABLE :: VTMP(:, :, :)

        GMIN = (/ INF, INF, INF /)
        GMAX = (/ -INF, -INF, -INF /)

        CALL GET_ENVIRONMENT_VARIABLE('LESWHT_GEOMETRY_STL', STL_PATH, LENGTH=LEN_OUT, STATUS=ISTAT)
        IF (ISTAT /= 0 .OR. LEN_OUT <= 0) THEN
            STL_PATH = '../output/geometry/dimple_slab_ascii.stl'
        ELSE
            STL_PATH = STL_PATH(1:LEN_OUT)
        END IF

        OPEN(NEWUNIT=IUNIT, FILE=TRIM(STL_PATH), STATUS='OLD', ACTION='READ', IOSTAT=IOS)
        IF (IOS /= 0) THEN
            GEOM_FAILED = .TRUE.
            RETURN
        END IF

        NVERT = 0_8
        DO
            READ(IUNIT, '(A)', IOSTAT=IOS) LINE
            IF (IOS /= 0) EXIT
            LINE = ADJUSTL(LINE)
            IF (INDEX(LINE, 'vertex') == 1) NVERT = NVERT + 1_8
        END DO
        CLOSE(IUNIT)

        IF (NVERT <= 0_8 .OR. MOD(NVERT, 3_8) /= 0_8) THEN
            GEOM_FAILED = .TRUE.
            RETURN
        END IF

        NTRI = NVERT / 3_8

        ALLOCATE(VTMP(3, 3, NTRI), STAT=IOS)
        IF (IOS /= 0) THEN
            GEOM_FAILED = .TRUE.
            RETURN
        END IF

        ALLOCATE(V1(3, NTRI), V2(3, NTRI), V3(3, NTRI), STAT=IOS)
        IF (IOS /= 0) THEN
            DEALLOCATE(VTMP)
            GEOM_FAILED = .TRUE.
            RETURN
        END IF

        ALLOCATE(BMINXZ(2, NTRI), BMAXXZ(2, NTRI), DEN2D(NTRI), PROJ_OK(NTRI), STAT=IOS)
        IF (IOS /= 0) THEN
            DEALLOCATE(VTMP, V1, V2, V3)
            GEOM_FAILED = .TRUE.
            RETURN
        END IF

        OPEN(NEWUNIT=IUNIT, FILE=TRIM(STL_PATH), STATUS='OLD', ACTION='READ', IOSTAT=IOS)
        IF (IOS /= 0) THEN
            CALL RELEASE_GEOMETRY()
            DEALLOCATE(VTMP)
            GEOM_FAILED = .TRUE.
            RETURN
        END IF

        IVERT = 0_8
        DO
            READ(IUNIT, '(A)', IOSTAT=IOS) LINE
            IF (IOS /= 0) EXIT
            LINE = ADJUSTL(LINE)
            IF (INDEX(LINE, 'vertex') == 1) THEN
                READ(LINE, *, IOSTAT=IOS) TAG, VX, VY, VZ
                IF (IOS /= 0) CYCLE
                IVERT = IVERT + 1_8
                ITRI = (IVERT - 1_8) / 3_8 + 1_8
                IVLOC = MOD(IVERT - 1_8, 3_8) + 1_8
                VTMP(IVLOC, 1, ITRI) = VX
                VTMP(IVLOC, 2, ITRI) = VY
                VTMP(IVLOC, 3, ITRI) = VZ
            END IF
        END DO
        CLOSE(IUNIT)

        IF (IVERT /= NVERT) THEN
            CALL RELEASE_GEOMETRY()
            DEALLOCATE(VTMP)
            GEOM_FAILED = .TRUE.
            RETURN
        END IF

        ! Add a tiny positive offset so the y=0 V-face is included as wall.
        YTOP_TARGET = 0.0D0 + EPS
        CALL GET_ENVIRONMENT_VARIABLE('LESWHT_GEOMETRY_Y_TOP', YTOP_STR, LENGTH=LEN_Y, STATUS=ISTAT_Y)
        IF (ISTAT_Y == 0 .AND. LEN_Y > 0) THEN
            READ(YTOP_STR(1:LEN_Y), *, IOSTAT=IOS) YTOP_TARGET
            IF (IOS /= 0) YTOP_TARGET = 0.0D0
        END IF

        YRAW_MIN = INF
        YRAW_MAX = -INF
        DO I = 1_8, NTRI
            YRAW_MIN = MIN(YRAW_MIN, MIN(VTMP(1, 2, I), MIN(VTMP(2, 2, I), VTMP(3, 2, I))))
            YRAW_MAX = MAX(YRAW_MAX, MAX(VTMP(1, 2, I), MAX(VTMP(2, 2, I), VTMP(3, 2, I))))
        END DO

        ! Estimate the lowest top-surface level (about -depth) from STL vertices,
        ! excluding deep slab-bottom vertices by a conservative threshold.
        YTHR = YRAW_MAX - 1.0D0
        YTOP_SURF_MIN = INF
        DO I = 1_8, NTRI
            if (VTMP(1, 2, I) .GE. YTHR) YTOP_SURF_MIN = MIN(YTOP_SURF_MIN, VTMP(1, 2, I))
            if (VTMP(2, 2, I) .GE. YTHR) YTOP_SURF_MIN = MIN(YTOP_SURF_MIN, VTMP(2, 2, I))
            if (VTMP(3, 2, I) .GE. YTHR) YTOP_SURF_MIN = MIN(YTOP_SURF_MIN, VTMP(3, 2, I))
        END DO
        if (YTOP_SURF_MIN .EQ. INF) YTOP_SURF_MIN = YRAW_MIN

        YSHIFT = YTOP_TARGET - YRAW_MAX
        YPAD_BOTTOM = YTOP_SURF_MIN + YSHIFT

        DO I = 1_8, NTRI
            V1(:, I) = VTMP(1, :, I)
            V2(:, I) = VTMP(2, :, I)
            V3(:, I) = VTMP(3, :, I)

            V1(2, I) = V1(2, I) + YSHIFT
            V2(2, I) = V2(2, I) + YSHIFT
            V3(2, I) = V3(2, I) + YSHIFT

            BMINXZ(1, I) = MIN(V1(1, I), MIN(V2(1, I), V3(1, I)))
            BMINXZ(2, I) = MIN(V1(3, I), MIN(V2(3, I), V3(3, I)))

            BMAXXZ(1, I) = MAX(V1(1, I), MAX(V2(1, I), V3(1, I)))
            BMAXXZ(2, I) = MAX(V1(3, I), MAX(V2(3, I), V3(3, I)))

            DEN2D(I) = (V2(3, I) - V3(3, I)) * (V1(1, I) - V3(1, I)) + &
                       (V3(1, I) - V2(1, I)) * (V1(3, I) - V3(3, I))

            PROJ_OK(I) = (ABS(DEN2D(I)) > EPS_PROJ)

            GMIN(1) = MIN(GMIN(1), MIN(V1(1, I), MIN(V2(1, I), V3(1, I))))
            GMIN(2) = MIN(GMIN(2), MIN(V1(2, I), MIN(V2(2, I), V3(2, I))))
            GMIN(3) = MIN(GMIN(3), MIN(V1(3, I), MIN(V2(3, I), V3(3, I))))

            GMAX(1) = MAX(GMAX(1), MAX(V1(1, I), MAX(V2(1, I), V3(1, I))))
            GMAX(2) = MAX(GMAX(2), MAX(V1(2, I), MAX(V2(2, I), V3(2, I))))
            GMAX(3) = MAX(GMAX(3), MAX(V1(3, I), MAX(V2(3, I), V3(3, I))))
        END DO

        DEALLOCATE(VTMP)

        XMIN = GMIN(1)
        XMAX = GMAX(1)
        ZMIN = GMIN(3)
        ZMAX = GMAX(3)

        CALL BUILD_XZ_BINS()
        IF (GEOM_FAILED) RETURN

        GEOM_READY = .TRUE.
    END SUBROUTINE LOAD_GEOMETRY


    SUBROUTINE BUILD_XZ_BINS()
        IMPLICIT NONE

        INTEGER(8) :: ITRI
        INTEGER(8) :: TOTAL_REFS
        INTEGER :: BASE
        INTEGER :: I0, I1, J0, J1, I, J, BIN
        INTEGER(8), ALLOCATABLE :: BIN_COUNT(:), BIN_CURSOR(:)
        REAL(8) :: XR, ZR

        IF (NTRI <= 0_8) THEN
            GEOM_FAILED = .TRUE.
            RETURN
        END IF

        XR = MAX(XMAX - XMIN, EPS)
        ZR = MAX(ZMAX - ZMIN, EPS)

        BASE = INT(SQRT(REAL(NTRI)))
        BASE = MAX(16, BASE)
        BASE = MIN(192, BASE)

        NBX = BASE
        NBZ = BASE
        NBINS = NBX * NBZ

        DXBIN = XR / REAL(NBX, 8)
        DZBIN = ZR / REAL(NBZ, 8)

        IF (DXBIN <= 0.0D0 .OR. DZBIN <= 0.0D0) THEN
            GEOM_FAILED = .TRUE.
            RETURN
        END IF

        ALLOCATE(BIN_COUNT(NBINS), STAT=I)
        IF (I /= 0) THEN
            GEOM_FAILED = .TRUE.
            RETURN
        END IF
        BIN_COUNT = 0_8

        DO ITRI = 1_8, NTRI
            IF (.NOT. PROJ_OK(ITRI)) CYCLE
            CALL TRI_BIN_RANGE(ITRI, I0, I1, J0, J1)
            DO J = J0, J1
                DO I = I0, I1
                    BIN = (J - 1) * NBX + I
                    BIN_COUNT(BIN) = BIN_COUNT(BIN) + 1_8
                END DO
            END DO
        END DO

        ALLOCATE(BIN_START(NBINS + 1), STAT=I)
        IF (I /= 0) THEN
            DEALLOCATE(BIN_COUNT)
            GEOM_FAILED = .TRUE.
            RETURN
        END IF

        BIN_START(1) = 1_8
        DO BIN = 1, NBINS
            BIN_START(BIN + 1) = BIN_START(BIN) + BIN_COUNT(BIN)
        END DO

        TOTAL_REFS = BIN_START(NBINS + 1) - 1_8
        IF (TOTAL_REFS <= 0_8) THEN
            DEALLOCATE(BIN_COUNT)
            GEOM_FAILED = .TRUE.
            RETURN
        END IF

        ALLOCATE(BIN_ITEMS(TOTAL_REFS), BIN_CURSOR(NBINS), STAT=I)
        IF (I /= 0) THEN
            DEALLOCATE(BIN_COUNT)
            GEOM_FAILED = .TRUE.
            RETURN
        END IF

        BIN_CURSOR = 0_8

        DO ITRI = 1_8, NTRI
            IF (.NOT. PROJ_OK(ITRI)) CYCLE
            CALL TRI_BIN_RANGE(ITRI, I0, I1, J0, J1)
            DO J = J0, J1
                DO I = I0, I1
                    BIN = (J - 1) * NBX + I
                    BIN_ITEMS(BIN_START(BIN) + BIN_CURSOR(BIN)) = ITRI
                    BIN_CURSOR(BIN) = BIN_CURSOR(BIN) + 1_8
                END DO
            END DO
        END DO

        DEALLOCATE(BIN_COUNT, BIN_CURSOR)
    END SUBROUTINE BUILD_XZ_BINS


    SUBROUTINE TRI_BIN_RANGE(ITRI, I0, I1, J0, J1)
        IMPLICIT NONE

        INTEGER(8), INTENT(IN) :: ITRI
        INTEGER, INTENT(OUT) :: I0, I1, J0, J1

        I0 = CLAMP_INT(INT((BMINXZ(1, ITRI) - XMIN) / DXBIN) + 1, 1, NBX)
        I1 = CLAMP_INT(INT((BMAXXZ(1, ITRI) - XMIN) / DXBIN) + 1, 1, NBX)
        J0 = CLAMP_INT(INT((BMINXZ(2, ITRI) - ZMIN) / DZBIN) + 1, 1, NBZ)
        J1 = CLAMP_INT(INT((BMAXXZ(2, ITRI) - ZMIN) / DZBIN) + 1, 1, NBZ)
    END SUBROUTINE TRI_BIN_RANGE


    LOGICAL FUNCTION HAS_BOTH_Y_HITS(XP, YP, ZP)
        IMPLICIT NONE

        REAL(8), INTENT(IN) :: XP, YP, ZP

        INTEGER :: IBIN, II, JJ
        INTEGER(8) :: POS, ITRI
        REAL(8) :: YINT
        LOGICAL :: VALID
        LOGICAL :: HIT_UP, HIT_DOWN

        HAS_BOTH_Y_HITS = .FALSE.

        IF (XP < XMIN - EPS .OR. XP > XMAX + EPS .OR. &
            ZP < ZMIN - EPS .OR. ZP > ZMAX + EPS) RETURN

        II = CLAMP_INT(INT((XP - XMIN) / DXBIN) + 1, 1, NBX)
        JJ = CLAMP_INT(INT((ZP - ZMIN) / DZBIN) + 1, 1, NBZ)
        IBIN = (JJ - 1) * NBX + II

        HIT_UP = .FALSE.
        HIT_DOWN = .FALSE.

        DO POS = BIN_START(IBIN), BIN_START(IBIN + 1) - 1_8
            ITRI = BIN_ITEMS(POS)
            CALL VERTICAL_Y_INTERSECTION(ITRI, XP, ZP, YINT, VALID)
            IF (.NOT. VALID) CYCLE

            IF (YINT >= YP + EPS_HIT) HIT_UP = .TRUE.
            IF (YINT <= YP - EPS_HIT) HIT_DOWN = .TRUE.

            IF (HIT_UP .AND. HIT_DOWN) THEN
                HAS_BOTH_Y_HITS = .TRUE.
                RETURN
            END IF
        END DO
    END FUNCTION HAS_BOTH_Y_HITS


    SUBROUTINE VERTICAL_Y_INTERSECTION(ITRI, XP, ZP, YINT, VALID)
        IMPLICIT NONE

        INTEGER(8), INTENT(IN) :: ITRI
        REAL(8), INTENT(IN) :: XP, ZP
        REAL(8), INTENT(OUT) :: YINT
        LOGICAL, INTENT(OUT) :: VALID

        REAL(8) :: W1, W2, W3

        VALID = .FALSE.
        YINT = 0.0D0

        IF (.NOT. PROJ_OK(ITRI)) RETURN

        IF (XP < BMINXZ(1, ITRI) - EPS .OR. XP > BMAXXZ(1, ITRI) + EPS) RETURN
        IF (ZP < BMINXZ(2, ITRI) - EPS .OR. ZP > BMAXXZ(2, ITRI) + EPS) RETURN

        W1 = ((V2(3, ITRI) - V3(3, ITRI)) * (XP - V3(1, ITRI)) + &
              (V3(1, ITRI) - V2(1, ITRI)) * (ZP - V3(3, ITRI))) / DEN2D(ITRI)

        W2 = ((V3(3, ITRI) - V1(3, ITRI)) * (XP - V3(1, ITRI)) + &
              (V1(1, ITRI) - V3(1, ITRI)) * (ZP - V3(3, ITRI))) / DEN2D(ITRI)

        W3 = 1.0D0 - W1 - W2

        IF (W1 < -EPS_BARY .OR. W2 < -EPS_BARY .OR. W3 < -EPS_BARY) RETURN

        YINT = W1 * V1(2, ITRI) + W2 * V2(2, ITRI) + W3 * V3(2, ITRI)
        VALID = .TRUE.
    END SUBROUTINE VERTICAL_Y_INTERSECTION


    INTEGER FUNCTION CLAMP_INT(V, LO, HI)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: V, LO, HI

        CLAMP_INT = V
        IF (CLAMP_INT < LO) CLAMP_INT = LO
        IF (CLAMP_INT > HI) CLAMP_INT = HI
    END FUNCTION CLAMP_INT


    SUBROUTINE RELEASE_GEOMETRY()
        IMPLICIT NONE

        IF (ALLOCATED(V1)) DEALLOCATE(V1)
        IF (ALLOCATED(V2)) DEALLOCATE(V2)
        IF (ALLOCATED(V3)) DEALLOCATE(V3)

        IF (ALLOCATED(BMINXZ)) DEALLOCATE(BMINXZ)
        IF (ALLOCATED(BMAXXZ)) DEALLOCATE(BMAXXZ)
        IF (ALLOCATED(DEN2D)) DEALLOCATE(DEN2D)
        IF (ALLOCATED(PROJ_OK)) DEALLOCATE(PROJ_OK)

        IF (ALLOCATED(BIN_START)) DEALLOCATE(BIN_START)
        IF (ALLOCATED(BIN_ITEMS)) DEALLOCATE(BIN_ITEMS)
    END SUBROUTINE RELEASE_GEOMETRY

END FUNCTION FUNCBODY
