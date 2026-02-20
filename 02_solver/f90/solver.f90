!=======================================================================
!
!     Codebase (LICA 2017 Version) by 
!     H. Choi / Department of Mechanical & Aerospace Engineering
!     Seoul National University
!
!=======================================================================
!
!     LESwHT (c) 2026 S. Lee (ORCID: 0000-0002-2063-6298)
!     2018.02.28. Updated (DOI:10.1016/j.ijheatmasstransfer.2019.01.019)
!     2026.02.18. Code modernization
!
!=======================================================================
PROGRAM SOLVER
    USE MOD_COMMON
    USE MOD_FLOWARRAY
    IMPLICIT NONE
    
    ! INTEGER(8) :: I, J, K, N
    REAL(8)    :: OMP_GET_WTIME

    !===================================================================
    !     INITIALIZATION PHASE
    !===================================================================
    CALL PRINT_REAL_TIME()                   ! AT MISC_INIT LIBRARY
                     
    TOTAL_TIME_B = OMP_GET_WTIME()           ! INTRINSIC SUBROUTINE
    CALL READSETTINGS()                      ! AT MISC_INIT LIBRARY
    CALL READBCS()                           ! AT MISC_INIT LIBRARY
    CALL READGEOM()                          ! AT MISC_INIT LIBRARY
    CALL ALLOINIT()                          ! AT MISC_INIT LIBRARY
    CALL ALLO_ARRAY()                        ! AT MISC_INIT LIBRARY

    ! --- LOAD PRE-PROCESSED DATA ---
    IF (IBMON .EQ. 1) THEN
        CALL IBMPREREAD()                    ! AT MISC_INIT LIBRARY
        IF (ILES .EQ. 1) THEN
            CALL NUTZEROREAD()               ! AT MISC_INIT LIBRARY
        ENDIF
        IF (ICONJG .EQ. 1) CALL CONJGREAD()  ! AT MISC_INIT LIBRARY
    ENDIF

    IF (ILES .EQ. 1) CALL SGSFILTERINIT()    ! AT SGS LIBRARY

    ! --- FIELD LOADING ---
    IF (IREAD .EQ. 1) THEN
        CALL PREFLD()                        ! AT MISC_INIT LIBRARY
    ELSE
        CALL MAKEFLD()                       ! AT MISC_INIT LIBRARY
        IF (EPS_PTR .NE. 0.0d0) THEN
            CALL ADDPERTURB() 
            ! Synchronize ghost cells so the Poisson solver doesn't blow up
            CALL PRDIC_ADJ_UVW(0) 
        ENDIF
        CALL WRITEFIELD()
    ENDIF

    IF (ICH .EQ. 1) CALL MEANPG()             ! AT SLV_MMTM LIBRARY

    ! --- MISCELLANEOUS SETTINGS ---
    CALL DTTIMEINIT()                        ! AT MISC_INIT LIBRARY
    CALL RK3COEFINIT()                       ! AT MISC_INIT LIBRARY
    CALL LHSINIT()                           ! AT SLV_MMTM LIBRARY
    CALL POISINIT()                          ! AT POISS LIBRARY
    CALL FTRFILES(1)

    CFLMAX = CFLFAC
    NV = 101
    NAV = 2001

    !===================================================================
    !     TIME-DEPENDENT CALCULATION (MAIN SOLVER)
    !===================================================================
    DO M = 1, NTST
        TIME_BEGIN = OMP_GET_WTIME()         ! INTRINSIC SUBROUTINE
        CALL STEPINIT()                      ! AT SLV_AUXI LIBRARY
        SUBDT = 0.0d0

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !     RK3 SUB-STEP LOOP
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        DO MSUB = 1, 3
            CALL SUBSTEPINIT(MSUB)           ! AT SLV_AUXI LIBRARY
            
            IF (ICH .EQ. 1) THEN
                CALL QVOLCALC(QVOL(1))       ! AT SLV_AUXI LIBRARY
            ENDIF

            SGSTIME_B(MSUB) = OMP_GET_WTIME()
            IF (ILES .EQ. 1) THEN
                CALL SGSCALC()               ! AT SGS LIBRARY
                IF (IHTRANS .EQ. 1) CALL SGSCALC_T() ! AT SGS LIBRARY
            ENDIF
            SGSTIME_E(MSUB) = OMP_GET_WTIME()

            IF (IMOVINGON .EQ. 1) CALL FINDFORCING() ! AT IBM_BODY (GLOBAL)

            RHSNLHSTIME_B(MSUB) = OMP_GET_WTIME()
            CALL RHSNLHS()                   ! AT SLV_MMTM LIBRARY
            RHSNLHSTIME_E(MSUB) = OMP_GET_WTIME()

            ALLOCATE(DIVGSUM(0:N1, 0:N2, 0:N3))
            ALLOCATE(PHI(0:N1, 0:N2, 0:N3))
            DIVGSUM = 0.0d0
            PHI = 0.0d0

            POISSTIME_B(MSUB) = OMP_GET_WTIME()
            CALL DIVGS(DIVGSUM)              ! AT SLV_CONT LIBRARY
            CALL POISSON(PHI, DIVGSUM)       ! AT POISS LIBRARY
            POISSTIME_E(MSUB) = OMP_GET_WTIME()

            IF (ICH .EQ. 1) THEN
                CALL QVOLCALC(QVOL(2))       ! AT SLV_AUXI LIBRARY
                CALL QVOLCORR()              ! AT SLV_AUXI LIBRARY
            ENDIF
            
            CALL UCALC(PHI)

            IF (ICH .EQ. 1) THEN
                CALL QVOLCALC(QVOL(0))       ! AT SLV_AUXI LIBRARY
            ENDIF

            CALL PCALC(PHI, DIVGSUM)

            IF (ICH .EQ. 1) CALL MEANPG()    ! AT SLV_AUXI LIBRARY
            IF (IBMON .EQ. 1) CALL LAGFORCE()

            IF (IHTRANS .EQ. 1) THEN
                IF (ICH .EQ. 1) CALL TVOLCALC(TVOL(1)) ! AT SLV_AUXI LIBRARY
                CALL RHSNLHS_T()                       ! AT SLV_ENGY LIBRARY
                IF (ICH .EQ. 1) CALL TVOLCALC(TVOL(2)) ! AT SLV_AUXI LIBRARY
                IF (ICH .EQ. 1) CALL TVOLCORR()        ! AT SLV_AUXI LIBRARY
                CALL TCALC()                           ! AT SLV_AUXI LIBRARY
                IF (ICH .EQ. 1) CALL TVOLCALC(TVOL(0)) ! AT SLV_AUXI LIBRARY
            ENDIF

            DEALLOCATE(DIVGSUM, PHI)
        ENDDO
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TIME = TIME + DT

        CALL CONVERGENCE_CHECK()
        IF ((IBMON .EQ. 1) .AND. (MASSON .EQ. 1)) CALL MASSCHECK()

        IF (MOD(NTIME, NPIN) .EQ. 0) THEN
            IHIST = IHIST + 1
            CALL WRITEHISTORY()
        ENDIF

        IF ((NTRACE .GT. 0) .AND. (MOD(M, NTR) .EQ. 0)) CALL TRACER()
        IF (IBMON .EQ. 1) THEN
            CALL DRAGLIFT()
        ENDIF
        IF (IHTRANS .EQ. 1) CALL CALC_BOUNDARY_HEAT_FLUX()

        IF (MOD(NTIME, NPRINT) .EQ. 0) CALL WRITEFIELD()
        IF (IAVG .EQ. 1) CALL FIELD_AVG()

        TIME_END = OMP_GET_WTIME()
        CALL WRITEFTRTIME()

    ENDDO

    !===================================================================
    !     FINISH AND CLEANUP
    !===================================================================
    CALL FTRFILES(0)

    IF (MOD(NTIME, NPRINT) .NE. 0) CALL WRITEFIELD()
    IF (IAVG .EQ. 1) CALL FIELD_AVG()

    CALL CPU_TIME(TOTAL_TIME_E)

    CALL PRINT_REAL_TIME()
    CALL DEALLO_ARRAY()

    STOP
END PROGRAM SOLVER