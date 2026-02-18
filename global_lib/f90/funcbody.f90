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
FUNCTION FUNCBODY(X, Y, Z, T) RESULT(VAL)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: X, Y, Z, T
    REAL(8) :: VAL

    !-------------------------------------------------------------------
    ! Define your body geometry here.
    ! VAL <= 0 : Inside Body (Solid)
    ! VAL >  0 : Outside Body (Fluid)
    !-------------------------------------------------------------------

    ! Example: A simple sphere
    ! VAL = SQRT(X**2 + Y**2 + Z**2) - 0.5d0

    ! Example: A channel (-0.5 < Y < 0.5) of the slab thickness of 0.5
    VAL = 0.5d0 - ABS(Y)

    RETURN
END FUNCTION FUNCBODY