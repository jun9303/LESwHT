!=======================================================================
!
!     Codebase (LICA 2017 Version) by 
!     H. Choi / Department of Mechanical & Aerospace Engineering
!     Seoul National University
!
!=======================================================================
!
!     LESwHT (c) 2026 S. Lee (ORCID: 0000-0002-2063-6298)
!     2018.02.28. Modified for F2PY (Fortran-to-Python) usage
!     2026.02.18. Code modernization
!
!=======================================================================
SUBROUTINE FINDFORCING()
    USE MOD_COMMON
    USE MOD_FLOWARRAY
    IMPLICIT NONE
    
    !-------------------------------------------------------------------
    ! 1. U-Velocity Component (Face X)
    !-------------------------------------------------------------------
    CALL FIND_INOUT(N1, N2, N3, X, YM, ZM, NBODY_U, INOUT_U, TIME)
    
    ! Only proceed if interpolation is ON
    CALL FINDBDY_INTP(1_8, N1, N2, N3, NBODY_U, INOUT_U, &
                      UFIX_X, UFIX_Y, UFIX_Z, LFIX_X, LFIX_Y, LFIX_Z, &
                      NINTP_U, NINNER_U, FCP_U, INTPTYPE_U, INTPINDX_U)
                      
    CALL GEOMFAC_INTP(N1, N2, N3, XPRE, YPRE, ZPRE, NINTP_U, &
                      NBODY_U, FCP_U, INTPINDX_U, GEOMFAC_U, TIME)

    !-------------------------------------------------------------------
    ! 2. V-Velocity Component (Face Y)
    !-------------------------------------------------------------------
    CALL FIND_INOUT(N1, N2, N3, XM, Y, ZM, NBODY_V, INOUT_V, TIME)
    
    CALL FINDBDY_INTP(2_8, N1, N2, N3, NBODY_V, INOUT_V, &
                      UFIX_X, UFIX_Y, UFIX_Z, LFIX_X, LFIX_Y, LFIX_Z, &
                      NINTP_V, NINNER_V, FCP_V, INTPTYPE_V, INTPINDX_V)
                      
    CALL GEOMFAC_INTP(N1, N2, N3, XPRE, YPRE, ZPRE, NINTP_V, &
                      NBODY_V, FCP_V, INTPINDX_V, GEOMFAC_V, TIME)

    !-------------------------------------------------------------------
    ! 3. W-Velocity Component (Face Z)
    !-------------------------------------------------------------------
    CALL FIND_INOUT(N1, N2, N3, XM, YM, Z, NBODY_W, INOUT_W, TIME)
    
    CALL FINDBDY_INTP(3_8, N1, N2, N3, NBODY_W, INOUT_W, &
                      UFIX_X, UFIX_Y, UFIX_Z, LFIX_X, LFIX_Y, LFIX_Z, &
                      NINTP_W, NINNER_W, FCP_W, INTPTYPE_W, INTPINDX_W)
                      
    CALL GEOMFAC_INTP(N1, N2, N3, XPRE, YPRE, ZPRE, NINTP_W, &
                      NBODY_W, FCP_W, INTPINDX_W, GEOMFAC_W, TIME)

    !-------------------------------------------------------------------
    ! 4. Temperature (Scalar) - Optional
    !-------------------------------------------------------------------
    IF (IHTRANS .EQ. 1) THEN
        CALL FIND_INOUT(N1, N2, N3, XM, YM, ZM, NBODY_T, INOUT_T, TIME)
        
        CALL FINDBDY_INTP(4_8, N1, N2, N3, NBODY_T, INOUT_T, &
                          UFIX_X, UFIX_Y, UFIX_Z, LFIX_X, LFIX_Y, LFIX_Z, &
                          NINTP_T, NINNER_T, FCP_T, INTPTYPE_T, INTPINDX_T)
                          
        CALL GEOMFAC_INTP(N1, N2, N3, XPRE, YPRE, ZPRE, NINTP_T, &
                          NBODY_T, FCP_T, INTPINDX_T, GEOMFAC_T, TIME)
    ENDIF

END SUBROUTINE FINDFORCING