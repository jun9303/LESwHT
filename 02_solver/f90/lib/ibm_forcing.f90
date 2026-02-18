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
    
    ! Local Temporary Variables
    INTEGER(8) :: NINNER_U, NINNER_V, NINNER_W, NINNER_T
    INTEGER(8), ALLOCATABLE :: FCP_TEMP(:,:)
    
    ! Allocate temporary array for Forcing Points (FCP)
    ! Size it safely to total grid points to handle any body size
    ALLOCATE(FCP_TEMP(N1*N2*N3, 3))

    !-------------------------------------------------------------------
    ! 1. U-Velocity Component (Face X)
    !    Grid: X (Face), YMP (Center), ZMP (Center)
    !-------------------------------------------------------------------
    CALL FIND_INOUT(N1, N2, N3, X, YMP, ZMP, NBODY(1), INOUT(:,:,:,1), TIME)
    
    CALL FINDBDY_INTP(1_8, N1, N2, N3, NBODY(1), INOUT(:,:,:,1), &
                      FIXIU, FIXJU, FIXKU, FIXIL, FIXJL, FIXKL, &
                      NINTP(1), NINNER_U, FCP_TEMP, INTPTYPE(:,1), INTPINDX(:,1,:))
    
    ! Copy FCP data to Module Arrays
    IFC(1:NBODY(1), 1) = FCP_TEMP(1:NBODY(1), 1)
    JFC(1:NBODY(1), 1) = FCP_TEMP(1:NBODY(1), 2)
    KFC(1:NBODY(1), 1) = FCP_TEMP(1:NBODY(1), 3)

    CALL GEOMFAC_INTP(N1, N2, N3, X, YMP, ZMP, NINTP(1), &
                      NBODY(1), FCP_TEMP, INTPINDX(:,1,:), GEOMFAC(:,1,:,:,:), TIME)

    !-------------------------------------------------------------------
    ! 2. V-Velocity Component (Face Y)
    !    Grid: XMP (Center), Y (Face), ZMP (Center)
    !-------------------------------------------------------------------
    CALL FIND_INOUT(N1, N2, N3, XMP, Y, ZMP, NBODY(2), INOUT(:,:,:,2), TIME)
    
    CALL FINDBDY_INTP(2_8, N1, N2, N3, NBODY(2), INOUT(:,:,:,2), &
                      FIXIU, FIXJU, FIXKU, FIXIL, FIXJL, FIXKL, &
                      NINTP(2), NINNER_V, FCP_TEMP, INTPTYPE(:,2), INTPINDX(:,2,:))
    
    IFC(1:NBODY(2), 2) = FCP_TEMP(1:NBODY(2), 1)
    JFC(1:NBODY(2), 2) = FCP_TEMP(1:NBODY(2), 2)
    KFC(1:NBODY(2), 2) = FCP_TEMP(1:NBODY(2), 3)

    CALL GEOMFAC_INTP(N1, N2, N3, XMP, Y, ZMP, NINTP(2), &
                      NBODY(2), FCP_TEMP, INTPINDX(:,2,:), GEOMFAC(:,2,:,:,:), TIME)

    !-------------------------------------------------------------------
    ! 3. W-Velocity Component (Face Z)
    !    Grid: XMP (Center), YMP (Center), Z (Face)
    !-------------------------------------------------------------------
    CALL FIND_INOUT(N1, N2, N3, XMP, YMP, Z, NBODY(3), INOUT(:,:,:,3), TIME)
    
    CALL FINDBDY_INTP(3_8, N1, N2, N3, NBODY(3), INOUT(:,:,:,3), &
                      FIXIU, FIXJU, FIXKU, FIXIL, FIXJL, FIXKL, &
                      NINTP(3), NINNER_W, FCP_TEMP, INTPTYPE(:,3), INTPINDX(:,3,:))
    
    IFC(1:NBODY(3), 3) = FCP_TEMP(1:NBODY(3), 1)
    JFC(1:NBODY(3), 3) = FCP_TEMP(1:NBODY(3), 2)
    KFC(1:NBODY(3), 3) = FCP_TEMP(1:NBODY(3), 3)

    CALL GEOMFAC_INTP(N1, N2, N3, XMP, YMP, Z, NINTP(3), &
                      NBODY(3), FCP_TEMP, INTPINDX(:,3,:), GEOMFAC(:,3,:,:,:), TIME)

    !-------------------------------------------------------------------
    ! 4. Temperature (Scalar) - Optional
    !    Grid: XMP, YMP, ZMP (All Center)
    !-------------------------------------------------------------------
    IF (IHTRANS .EQ. 1) THEN
        CALL FIND_INOUT(N1, N2, N3, XMP, YMP, ZMP, NBODY(4), INOUT(:,:,:,4), TIME)
        
        CALL FINDBDY_INTP(4_8, N1, N2, N3, NBODY(4), INOUT(:,:,:,4), &
                          FIXIU, FIXJU, FIXKU, FIXIL, FIXJL, FIXKL, &
                          NINTP(4), NINNER_T, FCP_TEMP, INTPTYPE(:,4), INTPINDX(:,4,:))
                          
        IFC(1:NBODY(4), 4) = FCP_TEMP(1:NBODY(4), 1)
        JFC(1:NBODY(4), 4) = FCP_TEMP(1:NBODY(4), 2)
        KFC(1:NBODY(4), 4) = FCP_TEMP(1:NBODY(4), 3)

        CALL GEOMFAC_INTP(N1, N2, N3, XMP, YMP, ZMP, NINTP(4), &
                          NBODY(4), FCP_TEMP, INTPINDX(:,4,:), GEOMFAC(:,4,:,:,:), TIME)
    ENDIF

    DEALLOCATE(FCP_TEMP)

END SUBROUTINE FINDFORCING