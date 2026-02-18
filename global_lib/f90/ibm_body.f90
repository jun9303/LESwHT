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
SUBROUTINE FIND_INOUT(NX, NY, NZ, XCOORD, YCOORD, ZCOORD, NBODY, INOUT, T)
    !$ USE OMP_LIB
    IMPLICIT NONE
    
    INTEGER(8), INTENT(IN)  :: NX, NY, NZ
    REAL(8),    INTENT(IN)  :: XCOORD(0:NX), YCOORD(0:NY), ZCOORD(0:NZ), T
    INTEGER(8), INTENT(OUT) :: INOUT(0:NX, 0:NY, 0:NZ)
    INTEGER(8), INTENT(OUT) :: NBODY

    INTEGER(8) :: I, J, K
    REAL(8)    :: VAL
    ! EXTERNAL FUNCTION DECLARATION
    REAL(8), EXTERNAL :: FUNCBODY

    NBODY = 0
    INOUT = 1

    !$OMP PARALLEL DO REDUCTION(+:NBODY) PRIVATE(I, J, K, VAL)
    DO K = 0, NZ
        DO J = 0, NY
            DO I = 0, NX
                VAL = FUNCBODY(XCOORD(I), YCOORD(J), ZCOORD(K), T)
                IF (VAL .LE. 1.0D-10) THEN
                    NBODY = NBODY + 1
                    INOUT(I, J, K) = 0
                ENDIF
            ENDDO
        ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    RETURN
END SUBROUTINE FIND_INOUT


!=======================================================================
SUBROUTINE FINDBDY_INTP(DIR, NX, NY, NZ, NBODY, INOUT, UFIX_X,   &
                        UFIX_Y, UFIX_Z, LFIX_X, LFIX_Y, LFIX_Z,  &
                        NINTP, NINNER, FCP, INTPTYPE, INTPINDX)
    !$ USE OMP_LIB
    IMPLICIT NONE
    
    INTEGER(8), INTENT(IN)  :: DIR, NX, NY, NZ, NBODY
    INTEGER(8), INTENT(IN)  :: INOUT(0:NX, 0:NY, 0:NZ)
    INTEGER(8), INTENT(IN)  :: UFIX_X(NX-1), UFIX_Y(NY-1), UFIX_Z(NZ-1)
    INTEGER(8), INTENT(IN)  :: LFIX_X(NX-1), LFIX_Y(NY-1), LFIX_Z(NZ-1)
    
    INTEGER(8), INTENT(OUT) :: NINTP, NINNER
    INTEGER(8), INTENT(OUT) :: FCP(NBODY, 3)
    INTEGER(8), INTENT(OUT) :: INTPTYPE(NBODY, 1), INTPINDX(NBODY, 3)

    INTEGER(8) :: I, J, K, IP, IM, JP, JM, KP, KM
    INTEGER(8) :: ISTART, JSTART, KSTART
    INTEGER(8) :: INOUT_SUM
    
    ! Temporary array for inner points to avoid race conditions or complex sorting
    INTEGER(8), ALLOCATABLE :: FCP_TEMP(:,:)

    ALLOCATE(FCP_TEMP(NBODY, 3))

    ISTART = 1; JSTART = 1; KSTART = 1
    IF (DIR .EQ. 1) ISTART = 2
    IF (DIR .EQ. 2) JSTART = 2
    IF (DIR .EQ. 3) KSTART = 2

    NINTP = 0
    NINNER = 0

    ! Serial loop to classify points (could be parallelized with offset calculation)
    DO K = KSTART, NZ-1
        KM = K-1; KP = K+1
        DO J = JSTART, NY-1
            JM = J-1; JP = J+1
            DO I = ISTART, NX-1
                IM = I-1; IP = I+1

                INOUT_SUM = (ABS(INOUT(IP,J,K)-INOUT(IM,J,K)))*(1-UFIX_X(I))*(1-LFIX_X(I)) &
                          + (ABS(INOUT(I,JP,K)-INOUT(I,JM,K)))*(1-UFIX_Y(J))*(1-LFIX_Y(J)) &
                          + (ABS(INOUT(I,J,KP)-INOUT(I,J,KM)))*(1-UFIX_Z(K))*(1-LFIX_Z(K))

                IF ((INOUT(I,J,K) .EQ. 0) .AND. (INOUT_SUM .GT. 0)) THEN
                    NINTP = NINTP + 1
                    FCP(NINTP, 1) = I
                    FCP(NINTP, 2) = J
                    FCP(NINTP, 3) = K
                    INTPTYPE(NINTP, 1) = INOUT_SUM
                    INTPINDX(NINTP, 1) = INOUT(IP,J,K) - INOUT(IM,J,K)
                    INTPINDX(NINTP, 2) = INOUT(I,JP,K) - INOUT(I,JM,K)
                    INTPINDX(NINTP, 3) = INOUT(I,J,KP) - INOUT(I,J,KM)
                ELSE IF ((INOUT(I,J,K) .EQ. 0) .AND. (INOUT_SUM .EQ. 0)) THEN
                    NINNER = NINNER + 1
                    FCP_TEMP(NINNER, 1) = I
                    FCP_TEMP(NINNER, 2) = J
                    FCP_TEMP(NINNER, 3) = K
                ENDIF
            ENDDO
        ENDDO
    ENDDO

    ! Copy inner points to the end of FCP
    !$OMP PARALLEL DO PRIVATE(I)
    DO I = 1, NINNER
        FCP(NINTP+I, 1) = FCP_TEMP(I, 1)
        FCP(NINTP+I, 2) = FCP_TEMP(I, 2)
        FCP(NINTP+I, 3) = FCP_TEMP(I, 3)
    ENDDO
    !$OMP END PARALLEL DO

    DEALLOCATE(FCP_TEMP)

END SUBROUTINE FINDBDY_INTP


!=======================================================================
SUBROUTINE GEOMFAC_PRESET(NI, ICOORD, IM, PRDIC, I_ADJ)
    !$ USE OMP_LIB
    IMPLICIT NONE
    
    INTEGER(8), INTENT(IN)  :: NI
    REAL(8),    INTENT(IN)  :: ICOORD(0:NI), IM(0:NI)
    CHARACTER(*), INTENT(IN):: PRDIC
    REAL(8),    INTENT(OUT) :: I_ADJ(-1:NI+1, 3)
    INTEGER(8) :: L

    I_ADJ = 0.0d0

    DO L = 1, NI
        I_ADJ(L, 1) = ICOORD(L)
        I_ADJ(L, 2) = IM(L)
        I_ADJ(L, 3) = IM(L)
    ENDDO
    I_ADJ(0, 2) = IM(0)
    I_ADJ(0, 3) = IM(0)

    IF (PRDIC .EQ. 'ON') THEN
        I_ADJ(0, 1)    = ICOORD(1) - (ICOORD(NI) - ICOORD(NI-1))
        I_ADJ(-1, 1)   = ICOORD(1) - (ICOORD(NI) - ICOORD(NI-2))
        I_ADJ(NI+1, 1) = ICOORD(NI) + (ICOORD(2) - ICOORD(1))

        I_ADJ(0, 2)    = ICOORD(1) - 0.5d0*(ICOORD(NI) - ICOORD(NI-1))
        I_ADJ(-1, 2)   = ICOORD(1) - (ICOORD(NI) - ICOORD(NI-1)) - 0.5d0*(ICOORD(NI-1) - ICOORD(NI-2))
        I_ADJ(NI, 2)   = ICOORD(NI) + 0.5d0*(ICOORD(2) - ICOORD(1))
        I_ADJ(NI+1, 2) = ICOORD(NI) + ICOORD(2) - ICOORD(1) + 0.5d0*(ICOORD(3) - ICOORD(2))

        I_ADJ(0, 3)    = I_ADJ(0, 2)
        I_ADJ(-1, 3)   = I_ADJ(-1, 2)
        I_ADJ(NI, 3)   = I_ADJ(NI, 2)
        I_ADJ(NI+1, 3) = I_ADJ(NI+1, 2)
    ENDIF
END SUBROUTINE GEOMFAC_PRESET


!=======================================================================
SUBROUTINE GEOMFAC_INTP(NX, NY, NZ, XPRE, YPRE, ZPRE, NINTP, &
                        NBODY, FCP, INTPINDX, GEOMFAC, T)
    !$ USE OMP_LIB
    IMPLICIT NONE
    
    INTEGER(8), INTENT(IN)  :: NX, NY, NZ, NINTP, NBODY
    REAL(8),    INTENT(IN)  :: XPRE(-1:NX+1), YPRE(-1:NY+1), ZPRE(-1:NZ+1)
    INTEGER(8), INTENT(IN)  :: FCP(NBODY, 3), INTPINDX(NINTP, 3)
    REAL(8),    INTENT(OUT) :: GEOMFAC(NINTP, 3, 3, 3)
    REAL(8),    INTENT(IN)  :: T

    INTEGER(8) :: L, M, N, I, J, K, INDIC
    REAL(8)    :: X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3
    REAL(8)    :: XTEMP1, YTEMP1, ZTEMP1, XTEMP2, YTEMP2, ZTEMP2
    REAL(8)    :: FFS, FFE, FF1, FF2
    REAL(8)    :: XX1, YY1, ZZ1, XX2, YY2, ZZ2
    REAL(8)    :: DDX, DDY, DDZ, DX1, DX2, DX3, DY1, DY2, DY3, DZ1, DZ2, DZ3
    REAL(8)    :: A0, B0, C0, A1, B1, C1
    
    REAL(8), EXTERNAL :: FUNCBODY

    !$OMP PARALLEL DO PRIVATE(L, M, N, I, J, K, INDIC, X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3) &
    !$OMP PRIVATE(XTEMP1, YTEMP1, ZTEMP1, XTEMP2, YTEMP2, ZTEMP2, FFS, FFE, FF1, FF2) &
    !$OMP PRIVATE(XX1, YY1, ZZ1, XX2, YY2, ZZ2, DDX, DDY, DDZ) &
    !$OMP PRIVATE(DX1, DX2, DX3, DY1, DY2, DY3, DZ1, DZ2, DZ3, A0, B0, C0, A1, B1, C1)
    DO L = 1, NINTP
        X1 = XPRE(FCP(L,1))
        Y1 = YPRE(FCP(L,2))
        Z1 = ZPRE(FCP(L,3))

        X2 = XPRE(FCP(L,1) + INTPINDX(L,1))
        Y2 = YPRE(FCP(L,2) + INTPINDX(L,2))
        Z2 = ZPRE(FCP(L,3) + INTPINDX(L,3))

        INDIC = 0
        X0 = X1; Y0 = Y1; Z0 = Z1 ! Initialize to avoid warning

        DO M = 1, 3
            IF (M .EQ. 1) THEN
                XTEMP1 = X1; YTEMP1 = Y1; ZTEMP1 = Z1
                XTEMP2 = X2; YTEMP2 = Y2; ZTEMP2 = Z2
                FFS = FUNCBODY(X1, Y1, Z1, T)
                FFE = FUNCBODY(X2, Y2, Z2, T)
                
                IF (FFS * FFE .GT. 0.0d0) THEN
                    ! Should not happen if surface is between 1 and 2, but fail-safe
                    GEOMFAC(L,:,:,:) = 0.0d0
                    GEOMFAC(L,1,1,1) = 1.0d0
                    GOTO 45
                ENDIF
            ELSE
                XTEMP1 = XX1; YTEMP1 = YY1; ZTEMP1 = ZZ1
                XTEMP2 = XX2; YTEMP2 = YY2; ZTEMP2 = ZZ2
            ENDIF

            ! Bisect/Linear search for surface
            DO N = 0, 19
                DDX = XTEMP2 - XTEMP1
                DDY = YTEMP2 - YTEMP1
                DDZ = ZTEMP2 - ZTEMP1

                XX1 = XTEMP1 + DDX * DBLE(N) / 20.0d0
                XX2 = XTEMP1 + DDX * DBLE(N+1) / 20.0d0
                YY1 = YTEMP1 + DDY * DBLE(N) / 20.0d0
                YY2 = YTEMP1 + DDY * DBLE(N+1) / 20.0d0
                ZZ1 = ZTEMP1 + DDZ * DBLE(N) / 20.0d0
                ZZ2 = ZTEMP1 + DDZ * DBLE(N+1) / 20.0d0
                
                FF1 = FUNCBODY(XX1, YY1, ZZ1, T)
                FF2 = FUNCBODY(XX2, YY2, ZZ2, T)

                IF (FF1 .EQ. 0.0d0) THEN
                    X0 = XX1; Y0 = YY1; Z0 = ZZ1
                    GOTO 33
                ELSE IF (FF2 .EQ. 0.0d0) THEN
                    X0 = XX2; Y0 = YY2; Z0 = ZZ2
                    GOTO 33
                ELSE IF (FF1 * FF2 .LT. 0.0d0) THEN
                    X0 = 0.5d0 * (XX1 + XX2)
                    Y0 = 0.5d0 * (YY1 + YY2)
                    Z0 = 0.5d0 * (ZZ1 + ZZ2)
                    IF (M .EQ. 3) GOTO 33
                    GOTO 22
                ENDIF
            ENDDO
 22         CONTINUE
        ENDDO
 33     CONTINUE

        X3 = XPRE(FCP(L,1) + INTPINDX(L,1)*2)
        Y3 = YPRE(FCP(L,2) + INTPINDX(L,2)*2)
        Z3 = ZPRE(FCP(L,3) + INTPINDX(L,3)*2)

        DX1 = ABS(X1-X0); DX2 = ABS(X2-X0); DX3 = ABS(X3-X0)
        DY1 = ABS(Y1-Y0); DY2 = ABS(Y2-Y0); DY3 = ABS(Y3-Y0)
        DZ1 = ABS(Z1-Z0); DZ2 = ABS(Z2-Z0); DZ3 = ABS(Z3-Z0)

        ! Weight Calculation (unchanged logic)
        IF (INTPINDX(L,1) .EQ. 0) THEN
            A0 = 1.0d0; A1 = 1.0d0
        ELSE IF (DX2 .GE. DX1) THEN
            A0 = DX2 / (DX1 + DX2)
            A1 = 1.0d0
        ELSE
            A0 = 0.5d0
            A1 = (DX3 - DX1) / (DX3 - DX2)
        ENDIF

        IF (INTPINDX(L,2) .EQ. 0) THEN
            B0 = 1.0d0; B1 = 1.0d0
        ELSE IF (DY2 .GE. DY1) THEN
            B0 = DY2 / (DY1 + DY2)
            B1 = 1.0d0
        ELSE
            B0 = 0.5d0
            B1 = (DY3 - DY1) / (DY3 - DY2)
        ENDIF

        IF (INTPINDX(L,3) .EQ. 0) THEN
            C0 = 1.0d0; C1 = 1.0d0
        ELSE IF (DZ2 .GE. DZ1) THEN
            C0 = DZ2 / (DZ1 + DZ2)
            C1 = 1.0d0
        ELSE
            C0 = 0.5d0
            C1 = (DZ3 - DZ1) / (DZ3 - DZ2)
        ENDIF

        ! Populate GEOMFAC (Direct assignments)
        GEOMFAC(L,1,1,1) =  1./(A0*B0*C0)
        GEOMFAC(L,1,1,2) = -1./(A0*B0*C0)*(A0)*(B0)*(1.-C0)*(C1)
        GEOMFAC(L,1,1,3) = -1./(A0*B0*C0)*(A0)*(B0)*(1.-C0)*(1.-C1)
        GEOMFAC(L,1,2,1) = -1./(A0*B0*C0)*(A0)*(1.-B0)*(C0)*(B1)
        GEOMFAC(L,1,2,2) = -1./(A0*B0*C0)*(A0)*(1.-B0)*(1.-C0)*(B1)*(C1)
        GEOMFAC(L,1,2,3) = -1./(A0*B0*C0)*(A0)*(1.-B0)*(1.-C0)*(B1)*(1.-C1)
        GEOMFAC(L,1,3,1) = -1./(A0*B0*C0)*(A0)*(1.-B0)*(C0)*(1.-B1)
        GEOMFAC(L,1,3,2) = -1./(A0*B0*C0)*(A0)*(1.-B0)*(1.-C0)*(1.-B1)*(C1)
        GEOMFAC(L,1,3,3) = -1./(A0*B0*C0)*(A0)*(1.-B0)*(1.-C0)*(1.-B1)*(1.-C1)

        GEOMFAC(L,2,1,1) = -1./(A0*B0*C0)*(1.-A0)*(B0)*(C0)*(A1)
        GEOMFAC(L,2,1,2) = -1./(A0*B0*C0)*(1.-A0)*(B0)*(1.-C0)*(A1)*(C1)
        GEOMFAC(L,2,1,3) = -1./(A0*B0*C0)*(1.-A0)*(B0)*(1.-C0)*(A1)*(1.-C1)
        GEOMFAC(L,2,2,1) = -1./(A0*B0*C0)*(1.-A0)*(1.-B0)*(C0)*(A1)*(B1)
        GEOMFAC(L,2,2,2) = -1./(A0*B0*C0)*(1.-A0)*(1.-B0)*(1.-C0)*(A1)*(B1)*(C1)
        GEOMFAC(L,2,2,3) = -1./(A0*B0*C0)*(1.-A0)*(1.-B0)*(1.-C0)*(A1)*(B1)*(1.-C1)
        GEOMFAC(L,2,3,1) = -1./(A0*B0*C0)*(1.-A0)*(1.-B0)*(C0)*(A1)*(1.-B1)
        GEOMFAC(L,2,3,2) = -1./(A0*B0*C0)*(1.-A0)*(1.-B0)*(1.-C0)*(A1)*(1.-B1)*(C1)
        GEOMFAC(L,2,3,3) = -1./(A0*B0*C0)*(1.-A0)*(1.-B0)*(1.-C0)*(A1)*(1.-B1)*(1.-C1)

        GEOMFAC(L,3,1,1) = -1./(A0*B0*C0)*(1.-A0)*(B0)*(C0)*(1.-A1)
        GEOMFAC(L,3,1,2) = -1./(A0*B0*C0)*(1.-A0)*(B0)*(1.-C0)*(1.-A1)*(C1)
        GEOMFAC(L,3,1,3) = -1./(A0*B0*C0)*(1.-A0)*(B0)*(1.-C0)*(1.-A1)*(1.-C1)
        GEOMFAC(L,3,2,1) = -1./(A0*B0*C0)*(1.-A0)*(1.-B0)*(C0)*(1.-A1)*(B1)
        GEOMFAC(L,3,2,2) = -1./(A0*B0*C0)*(1.-A0)*(1.-B0)*(1.-C0)*(1.-A1)*(B1)*(C1)
        GEOMFAC(L,3,2,3) = -1./(A0*B0*C0)*(1.-A0)*(1.-B0)*(1.-C0)*(1.-A1)*(B1)*(1.-C1)
        GEOMFAC(L,3,3,1) = -1./(A0*B0*C0)*(1.-A0)*(1.-B0)*(C0)*(1.-A1)*(1.-B1)
        GEOMFAC(L,3,3,2) = -1./(A0*B0*C0)*(1.-A0)*(1.-B0)*(1.-C0)*(1.-A1)*(1.-B1)*(C1)
        GEOMFAC(L,3,3,3) = -1./(A0*B0*C0)*(1.-A0)*(1.-B0)*(1.-C0)*(1.-A1)*(1.-B1)*(1.-C1)
45      CONTINUE
    ENDDO
    !$OMP END PARALLEL DO

END SUBROUTINE GEOMFAC_INTP


!=======================================================================
SUBROUTINE FINDBDY_NOINTP(DIR, NX, NY, NZ, NBODY, INOUT, NINNER, FCP)
    IMPLICIT NONE
    INTEGER(8), INTENT(IN)  :: DIR, NX, NY, NZ, NBODY
    INTEGER(8), INTENT(IN)  :: INOUT(0:NX, 0:NY, 0:NZ)
    INTEGER(8), INTENT(OUT) :: NINNER
    INTEGER(8), INTENT(OUT) :: FCP(NBODY, 3)

    INTEGER(8) :: I, J, K, ISTART, JSTART, KSTART

    ISTART = 1; JSTART = 1; KSTART = 1
    IF (DIR .EQ. 1) ISTART = 2
    IF (DIR .EQ. 2) JSTART = 2
    IF (DIR .EQ. 3) KSTART = 2

    NINNER = 0

    DO K = KSTART, NZ-1
        DO J = JSTART, NY-1
            DO I = ISTART, NX-1
                IF (INOUT(I,J,K) .EQ. 0) THEN
                    NINNER = NINNER + 1
                    FCP(NINNER, 1) = I
                    FCP(NINNER, 2) = J
                    FCP(NINNER, 3) = K
                ENDIF
            ENDDO
        ENDDO
    ENDDO

END SUBROUTINE FINDBDY_NOINTP


!=======================================================================
SUBROUTINE FIND_ZERO_NU_SGS(NX, NY, NZ, XM, YM, ZM, NZERO, INOUT, T)
    !$ USE OMP_LIB
    IMPLICIT NONE
    
    INTEGER(8), INTENT(IN)  :: NX, NY, NZ
    REAL(8),    INTENT(IN)  :: XM(0:NX), YM(0:NY), ZM(0:NZ), T
    INTEGER(8), INTENT(OUT) :: INOUT(1:NX-1, 1:NY-1, 1:NZ-1)
    INTEGER(8), INTENT(OUT) :: NZERO

    INTEGER(8) :: I, J, K
    REAL(8), EXTERNAL :: FUNCBODY

    NZERO = 0
    INOUT = 1

    !$OMP PARALLEL DO REDUCTION(+:NZERO) PRIVATE(I, J, K)
    DO K = 1, NZ-1
        DO J = 1, NY-1
            DO I = 1, NX-1
                ! Check 7-point stencil. If any point is inside body, zero viscosity.
                IF ( (FUNCBODY(XM(I),   YM(J),   ZM(K),   T) .LE. 1.0D-10) .OR. &
                     (FUNCBODY(XM(I-1), YM(J),   ZM(K),   T) .LE. 1.0D-10) .OR. &
                     (FUNCBODY(XM(I+1), YM(J),   ZM(K),   T) .LE. 1.0D-10) .OR. &
                     (FUNCBODY(XM(I),   YM(J-1), ZM(K),   T) .LE. 1.0D-10) .OR. &
                     (FUNCBODY(XM(I),   YM(J+1), ZM(K),   T) .LE. 1.0D-10) .OR. &
                     (FUNCBODY(XM(I),   YM(J),   ZM(K-1), T) .LE. 1.0D-10) .OR. &
                     (FUNCBODY(XM(I),   YM(J),   ZM(K+1), T) .LE. 1.0D-10) ) THEN
                    NZERO = NZERO + 1
                    INOUT(I, J, K) = 0
                ENDIF
            ENDDO
        ENDDO
    ENDDO
    !$OMP END PARALLEL DO

END SUBROUTINE FIND_ZERO_NU_SGS


!=======================================================================
SUBROUTINE CONJG_INTP(NX, NY, NZ, CRATIO, KRATIO, XM, YM, ZM, X, Y, Z, &
                      ISZERO, T, CSTAR, KSTAR)
    !$ USE OMP_LIB
    IMPLICIT NONE
    
    INTEGER(8), INTENT(IN)  :: NX, NY, NZ
    REAL(8),    INTENT(IN)  :: CRATIO, KRATIO
    REAL(8),    INTENT(IN)  :: XM(0:NX), YM(0:NY), ZM(0:NZ)
    REAL(8),    INTENT(IN)  :: X(0:NX), Y(0:NY), Z(0:NZ)
    INTEGER(8), INTENT(IN)  :: ISZERO(1:NX-1, 1:NY-1, 1:NZ-1)
    REAL(8),    INTENT(IN)  :: T
    
    REAL(8),    INTENT(OUT) :: CSTAR(1:NX-1, 1:NY-1, 1:NZ-1)
    REAL(8),    INTENT(OUT) :: KSTAR(1:NX-1, 1:NY-1, 1:NZ-1, 6)

    INTEGER(8) :: I, J, K, SUBC
    REAL(8)    :: FPTEMP, AA
    REAL(8), EXTERNAL :: FUNCBODY, FLUID_PORTION

    ! Subdivision level for Flood Fill
    SUBC = 6  ! Optimized value from Lee & Hwang (2019)

    !$OMP PARALLEL DO PRIVATE(I, J, K, FPTEMP, AA)
    DO K = 1, NZ-1
        DO J = 1, NY-1
            DO I = 1, NX-1
                ! AA > 0 if Center is Fluid. AA <= 0 if Center is Solid.
                AA = FUNCBODY(XM(I), YM(J), ZM(K), T)

                ! Heat Capacity C*
                FPTEMP = FLUID_PORTION(X(I), X(I+1), Y(J), Y(J+1), Z(K), Z(K+1), T, SUBC)
                CSTAR(I,J,K) = (1.0d0 - FPTEMP)*CRATIO + FPTEMP*1.0d0

                ! Thermal Conductivity K* (6 faces)
                
                ! East Face (i+1/2) - Center (XM) to (XM+1)
                ! Using 'Interim' cell logic from paper
                FPTEMP = FLUID_PORTION(XM(I), XM(I+1), Y(J), Y(J+1), Z(K), Z(K+1), T, SUBC)
                IF (AA * FUNCBODY(XM(I+1), YM(J), ZM(K), T) .GE. 0.0d0) THEN
                    ! Same phase
                    KSTAR(I,J,K,1) = (1.0d0 - FPTEMP)*KRATIO + FPTEMP*1.0d0
                ELSE
                    ! Harmonic mean for interface
                    KSTAR(I,J,K,1) = KRATIO / (KRATIO*FPTEMP + 1.0d0*(1.0d0-FPTEMP))
                ENDIF

                ! West Face (i-1/2)
                FPTEMP = FLUID_PORTION(XM(I-1), XM(I), Y(J), Y(J+1), Z(K), Z(K+1), T, SUBC)
                IF (AA * FUNCBODY(XM(I-1), YM(J), ZM(K), T) .GE. 0.0d0) THEN
                    KSTAR(I,J,K,2) = (1.0d0 - FPTEMP)*KRATIO + FPTEMP*1.0d0
                ELSE
                    KSTAR(I,J,K,2) = KRATIO / (KRATIO*FPTEMP + 1.0d0*(1.0d0-FPTEMP))
                ENDIF

                ! North Face (j+1/2)
                FPTEMP = FLUID_PORTION(X(I), X(I+1), YM(J), YM(J+1), Z(K), Z(K+1), T, SUBC)
                IF (AA * FUNCBODY(XM(I), YM(J+1), ZM(K), T) .GE. 0.0d0) THEN
                    KSTAR(I,J,K,3) = (1.0d0 - FPTEMP)*KRATIO + FPTEMP*1.0d0
                ELSE
                    KSTAR(I,J,K,3) = KRATIO / (KRATIO*FPTEMP + 1.0d0*(1.0d0-FPTEMP))
                ENDIF

                ! South Face (j-1/2)
                FPTEMP = FLUID_PORTION(X(I), X(I+1), YM(J-1), YM(J), Z(K), Z(K+1), T, SUBC)
                IF (AA * FUNCBODY(XM(I), YM(J-1), ZM(K), T) .GE. 0.0d0) THEN
                    KSTAR(I,J,K,4) = (1.0d0 - FPTEMP)*KRATIO + FPTEMP*1.0d0
                ELSE
                    KSTAR(I,J,K,4) = KRATIO / (KRATIO*FPTEMP + 1.0d0*(1.0d0-FPTEMP))
                ENDIF

                ! Top Face (k+1/2)
                FPTEMP = FLUID_PORTION(X(I), X(I+1), Y(J), Y(J+1), ZM(K), ZM(K+1), T, SUBC)
                IF (AA * FUNCBODY(XM(I), YM(J), ZM(K+1), T) .GE. 0.0d0) THEN
                    KSTAR(I,J,K,5) = (1.0d0 - FPTEMP)*KRATIO + FPTEMP*1.0d0
                ELSE
                    KSTAR(I,J,K,5) = KRATIO / (KRATIO*FPTEMP + 1.0d0*(1.0d0-FPTEMP))
                ENDIF

                ! Bottom Face (k-1/2)
                FPTEMP = FLUID_PORTION(X(I), X(I+1), Y(J), Y(J+1), ZM(K-1), ZM(K), T, SUBC)
                IF (AA * FUNCBODY(XM(I), YM(J), ZM(K-1), T) .GE. 0.0d0) THEN
                    KSTAR(I,J,K,6) = (1.0d0 - FPTEMP)*KRATIO + FPTEMP*1.0d0
                ELSE
                    KSTAR(I,J,K,6) = KRATIO / (KRATIO*FPTEMP + 1.0d0*(1.0d0-FPTEMP))
                ENDIF
            ENDDO
        ENDDO
    ENDDO
    !$OMP END PARALLEL DO

END SUBROUTINE CONJG_INTP


!=======================================================================
FUNCTION FLUID_PORTION(X1, X2, Y1, Y2, Z1, Z2, T, DIV) RESULT(VAL)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: X1, X2, Y1, Y2, Z1, Z2, T
    INTEGER(8), INTENT(IN) :: DIV
    REAL(8) :: VAL

    ! Local variables for Flood Fill
    REAL(8) :: DX, DY, DZ, XC, YC, ZC
    INTEGER(8) :: I, J, K, II, JJ, KK
    INTEGER(8) :: SOLID_COUNT, START_PHASE
    
    ! Automatic Arrays (Thread-safe on Stack)
    ! PHASE: -1 (Unknown), 0 (Solid), 1 (Fluid)
    INTEGER(8) :: PHASE(DIV, DIV, DIV)
    
    ! Stack for recursion (Iterative implementation)
    INTEGER(8) :: STACK(DIV*DIV*DIV, 3) 
    INTEGER(8) :: STACK_PTR
    
    REAL(8), EXTERNAL :: FUNCBODY

    DX = (X2 - X1) / DBLE(DIV)
    DY = (Y2 - Y1) / DBLE(DIV)
    DZ = (Z2 - Z1) / DBLE(DIV)

    ! Initialize
    PHASE = -1
    STACK_PTR = 0
    SOLID_COUNT = 0

    ! --- STEP 1: INITIATION (Check Corner 1,1,1) ---
    XC = X1 + DX * 0.5d0
    YC = Y1 + DY * 0.5d0
    ZC = Z1 + DZ * 0.5d0
    
    IF (FUNCBODY(XC, YC, ZC, T) .LE. 1.0D-10) THEN
        START_PHASE = 0 ! Solid
    ELSE
        START_PHASE = 1 ! Fluid
    ENDIF
    
    PHASE(1, 1, 1) = START_PHASE
    STACK_PTR = STACK_PTR + 1
    STACK(STACK_PTR, 1) = 1
    STACK(STACK_PTR, 2) = 1
    STACK(STACK_PTR, 3) = 1

    ! --- STEP 2: FLOOD FILL ---
    DO WHILE (STACK_PTR .GT. 0)
        ! Pop
        I = STACK(STACK_PTR, 1)
        J = STACK(STACK_PTR, 2)
        K = STACK(STACK_PTR, 3)
        STACK_PTR = STACK_PTR - 1

        ! Check 6 neighbors
        DO KK = K-1, K+1
        DO JJ = J-1, J+1
        DO II = I-1, I+1
            ! Check Manhatten distance = 1 (Neighbors only, no diagonals, no self)
            IF (ABS(II-I)+ABS(JJ-J)+ABS(KK-K) .NE. 1) CYCLE
            
            ! Boundary Check
            IF (II .LT. 1 .OR. II .GT. DIV) CYCLE
            IF (JJ .LT. 1 .OR. JJ .GT. DIV) CYCLE
            IF (KK .LT. 1 .OR. KK .GT. DIV) CYCLE
            
            ! If already visited/identified, skip
            IF (PHASE(II, JJ, KK) .NE. -1) CYCLE

            ! Identify Neighbor
            XC = X1 + (X2-X1)/DBLE(DIV) * (DBLE(2*II-1)/2.0d0)
            YC = Y1 + (Y2-Y1)/DBLE(DIV) * (DBLE(2*JJ-1)/2.0d0)
            ZC = Z1 + (Z2-Z1)/DBLE(DIV) * (DBLE(2*KK-1)/2.0d0)

            IF (FUNCBODY(XC, YC, ZC, T) .LE. 1.0D-10) THEN
                PHASE(II, JJ, KK) = 0 ! Solid
            ELSE
                PHASE(II, JJ, KK) = 1 ! Fluid
            ENDIF

            ! If Neighbor has SAME phase as Start, push to stack (continue flood)
            ! If different, it is a boundary, so we stop flooding that path.
            IF (PHASE(II, JJ, KK) .EQ. START_PHASE) THEN
                STACK_PTR = STACK_PTR + 1
                STACK(STACK_PTR, 1) = II
                STACK(STACK_PTR, 2) = JJ
                STACK(STACK_PTR, 3) = KK
            ENDIF
        ENDDO
        ENDDO
        ENDDO
    ENDDO

    ! --- STEP 3: TERMINATION & COUNTING ---
    ! Any sub-cell with PHASE == -1 is inaccessible from (1,1,1).
    ! Assuming single interface, these must be the OPPOSITE phase of START_PHASE.
    DO K = 1, DIV
        DO J = 1, DIV
            DO I = 1, DIV
                IF (PHASE(I, J, K) .EQ. -1) THEN
                    PHASE(I, J, K) = 1 - START_PHASE
                ENDIF
                
                IF (PHASE(I, J, K) .EQ. 0) THEN
                    SOLID_COUNT = SOLID_COUNT + 1
                ENDIF
            ENDDO
        ENDDO
    ENDDO

    VAL = 1.0d0 - DBLE(SOLID_COUNT) / DBLE(DIV**3)

END FUNCTION FLUID_PORTION