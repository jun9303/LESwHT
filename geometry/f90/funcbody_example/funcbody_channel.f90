FUNCTION FUNCBODY(X, Y, Z, T) RESULT(VAL)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: X, Y, Z, T
    REAL(8) :: VAL

    !-------------------------------------------------------------------
    ! Define your body geometry here.
    ! VAL <= 0 : Inside Body (Solid)
    ! VAL >  0 : Outside Body (Fluid)
    !-------------------------------------------------------------------

    ! Example: A channel (-1. < Y < 1.)
    VAL = 1.d0 - ABS(Y)

    RETURN
END FUNCTION FUNCBODY