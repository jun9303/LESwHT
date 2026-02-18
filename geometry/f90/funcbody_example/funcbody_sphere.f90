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
    VAL = SQRT(X**2 + Y**2 + Z**2) - 0.5d0

    RETURN
END FUNCTION FUNCBODY