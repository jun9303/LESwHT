FUNCTION FUNCBODY(X, Y, Z, T) RESULT(VAL)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: X, Y, Z, T
    REAL(8) :: VAL
    REAL(8), PARAMETER :: EPS = 1.0d-6 ! Small offset to prevent exact node alignment

    !-------------------------------------------------------------------
    ! Define your body geometry here.
    ! VAL <= 0 : Inside Body (Solid)
    ! VAL >  0 : Outside Body (Fluid)
    !-------------------------------------------------------------------

    ! Example: A channel (-1. < Y < 1.)
    VAL = (1.d0 + EPS) - ABS(Y)

    RETURN
END FUNCTION FUNCBODY