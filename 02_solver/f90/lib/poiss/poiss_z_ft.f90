!=======================================================================
!
!     MAIN SOLVER OF POISSON EQUATION
!
!     AX=b
!     A: coefficient of discretized poisson equation
!     X: PHI: Pseudo pressure, Output of subroutine POISSON
!     b: DIVGSUM: OUTPUT of subroutine DIVGS.
!
!     x-direction: MULTI-GRID iteration/GSOR method
!     y-direction: TDMA
!     z-direction: Fourier transform
!
!     Jun. 2017, J. Park
!     Feb. 2026, S. Lee (Several bug fixes)
!
!=======================================================================
      SUBROUTINE POISINIT
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8)   :: I, J, K, ILEV

       CALL Z_FT_ALLO

       ALLOCATE(AC(N1MD, N2M, N3MH))
       ALLOCATE(GAM(N1MD, N2M, N3MH))
       ALLOCATE(BET(N1MD, N2M, N3MH))
 
       CALL PMAT
!------MULTIGRID METHOD
       CALL COEFMG

       DO ILEV = 0, NLEV
           DO K = 1, N3MH
               DO I = IIMG(ILEV, 1), IIMG(ILEV, 2)
                   BET(I, 1, K) = 1.0_8 / AC(I, 1, K)
               END DO
               DO J = 2, N2M
                   DO I = IIMG(ILEV, 1), IIMG(ILEV, 2)
                       GAM(I, J, K) = AN(J - 1) * BET(I, J - 1, K)
                       BET(I, J, K) = 1.0_8 / (AC(I, J, K) - AS(J) * GAM(I, J, K))
                   END DO
               END DO
           END DO
       END DO

      OPEN(77, FILE='../output/ftr/mgftresiduemax.dat')

      RETURN
      END SUBROUTINE POISINIT

!=======================================================================
      SUBROUTINE PMAT
!=======================================================================
! --- CONSTRUCT MATRIX FOR POISSON EQ. ---------------------------------
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8)   :: I, J, K, IP, JP

!-----FOR TOP LEVEL
      DO I = 1, N1M
          IP = IPV(I)
          AW(I) = (1.0_8 - FIXIL(I)) * C2CXI(I) * F2FXI(I)
          AE(I) = (1.0_8 - FIXIU(I)) * C2CXI(IP) * F2FXI(I)
      END DO

      DO J = 1, N2M
          JP = JPV(J)
          AS(J) = (1.0_8 - FIXJL(J)) * C2CYI(J) * F2FYI(J)
          AN(J) = (1.0_8 - FIXJU(J)) * C2CYI(JP) * F2FYI(J)
      END DO

      DO J = 1, N2M
          DO I = 1, N1M
              AC(I, J, 1) = -1.0_8 * (AW(I) + AE(I) + AS(J) + AN(J))
          END DO
      END DO

      N3MH = N3M / 2 + 1
      IF(N3M > 1) CALL MWAVENUMBER     ! INIT. MODIFIED WAVE #.

      DO K = 2, N3MH
          DO J = 1, N2M
              DO I = 1, N1M
                  AC(I, J, K) = AC(I, J, 1) - AK3(K)
              END DO
          END DO
      END DO

      RETURN
      END SUBROUTINE PMAT

!=======================================================================
      SUBROUTINE MWAVENUMBER
!=======================================================================
! --- MODIFIED WAVE NUMBER DEFINITION
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8)   :: K
       REAL(8)      :: PI

      PI = ACOS(-1.0_8)

      DO K = 1, N3MH
          AK3(K) = 2.0_8 * (1.0_8 - COS(2.0_8 * FLOAT(K - 1) * PI / FLOAT(N3M))) * F2FZI(1) * F2FZI(1)
      END DO

      RETURN
      END SUBROUTINE MWAVENUMBER

!=======================================================================
      SUBROUTINE POISSON(PHI, DIVGSUM)
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8)   :: I, J, K, KKK
       REAL(8)      :: PHI(0:N1,0:N2,0:N3), DIVGSUM(0:N1,0:N2,0:N3)
       COMPLEX(8)   :: CCAP(N1,N2,N3MH)
       COMPLEX(8), DIMENSION(:,:), ALLOCATABLE :: XXX, CP, TMP
       COMPLEX(8), DIMENSION(:),   ALLOCATABLE :: XXXX, XXXX_B

       REAL(8)      :: TEST

      IF (N3M == 1) THEN
          ALLOCATE(CP(0:N1, 0:N2), TMP(N1, N2))
          CP = 0.0_8
!$OMP PARALLEL DO
          DO J = 1, N2M
              DO I = 1, N1M
                  TMP(I, J) = DIVGSUM(I, J, 1)
              END DO
          END DO

          CALL MG2D(CP, TMP, 1_8, TEST1, 0.0_8)

!$OMP PARALLEL DO
          DO J = 1, N2M
              DO I = 1, N1M
                  PHI(I, J, 1) = REAL(CP(I, J), 8)
              END DO
          END DO
          DEALLOCATE(CP, TMP)
          RETURN
      END IF

! --- DO THE FORWARD FFT
!$OMP PARALLEL private(XXX, XXXX, XXXX_B)
      ALLOCATE(XXX(N3M, N1M))
      ALLOCATE(XXXX(N3M), XXXX_B(N3M*2))
      CALL ZFFT1D(XXXX, N3M, 0, XXXX_B)
!$OMP DO
      DO J = 1, N2M
          DO K = 1, N3M
              DO I = 1, N1M
                  XXX(K, I) = DIVGSUM(I, J, K)
              END DO
          END DO
        
          DO I = 1, N1M
              CALL ZFFT1D(XXX(1, I), N3M, -1, XXXX_B)
          END DO
        
          DO K = 1, N3MH
              DO I = 1, N1M
                  CCAP(I, J, K) = XXX(K, I)
              END DO
          END DO
      END DO
!$OMP END DO
      DEALLOCATE(XXX, XXXX, XXXX_B)
!$OMP END PARALLEL

! --- SOLVE A SET OF POISSON EQS.
       TEST = TEST1 / FLOAT(N3MH) * FLOAT(N3M) * 0.8_8   ! convergence criteria

!$OMP PARALLEL private(CP)
      ALLOCATE(CP(0:N1, 0:N2))
!$OMP DO
      DO KKK = 1, N3MH
          CP = 0.0_8
          IF (KKK <= IMGSOR) THEN    
              CALL MG2D(CP, CCAP(1, 1, KKK), KKK, TEST, 0.0_8)
          ELSE
              CALL GSOR2D(CP, CCAP(1, 1, KKK), KKK, TEST, 0.0_8)
          END IF

          DO J = 1, N2M
              DO I = 1, N1M
                  CCAP(I, J, KKK) = CP(I, J)
              END DO
          END DO
      END DO
!$OMP END DO
      DEALLOCATE(CP)
!$OMP END PARALLEL

! --- DO THE INVERSE FFT
!$OMP PARALLEL private(XXX, XXXX, XXXX_B)
      ALLOCATE(XXX(N3M, N1M))
      ALLOCATE(XXXX(N3M), XXXX_B(N3M*2))
      CALL ZFFT1D(XXXX, N3M, 0, XXXX_B)
!$OMP DO
      DO J = 1, N2M
          DO K = 1, N3MH
              DO I = 1, N1M
                  XXX(K, I) = CCAP(I, J, K)
              END DO
          END DO

          DO K = N3MH + 1, N3M
              DO I = 1, N1M
                  XXX(K, I) = CONJG(CCAP(I, J, N3M + 2 - K))
              END DO
          END DO
        
          DO I = 1, N1M
              CALL ZFFT1D(XXX(1, I), N3M, 1, XXXX_B)
          END DO
        
          DO K = 1, N3M
              DO I = 1, N1M
                  PHI(I, J, K) = REAL(XXX(K, I), 8)
              END DO
          END DO
      END DO
!$OMP END DO
      DEALLOCATE(XXX, XXXX, XXXX_B)
!$OMP END PARALLEL

      RETURN
      END SUBROUTINE POISSON

!=======================================================================
       SUBROUTINE MG2D(PC, RHS, KV, TEST, OLDV)
!=======================================================================
!     MULTIGRID ENVIRONMENT VARIABLES
!     PC  : SOLUTION OF THE MG2D.
!     RHS : RHS OF THE POISSON EQUATION FOR EACH WAVENUMBER
!     KV  : wavenumber index
!     TEST    : CONDITION FOR CONVERGENCE
!     OLDV    : 0; INITALIZE TO ZERO, 1; USE PREVIOUS SOLUTION
!     NLEV    : TOTAL CELS = (MINROW1+MINROW2)*(2**NLEV)
!     LEVHALF : HALF OF NLEV
!     IWC     : FLAG FOR USING W-CYCLE ( 1 IF YOU WANT TO USE W-CYCLE)
!     NBLI    : NUMBER OF BASE LEVEL ITERATION
!     MGITR   : NUMBER OF MAXIMUM ITERATIONS
!     IIMG(ILEV,1)  : START POINT OF INDEX AT ILEV LEVEL
!     IIMG(ILEV,2)  : END POINT OF INDEX AT ILEV LEVEL
!     IIMG(ILEV,3)  : NUMBER OF ROW AT ILEV LEVEL
!-----------------------------------------------------------------------
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8)   :: I, J, KV, KK, ILEV
       REAL(8)      :: TEST, SUMRES, OLDV
       COMPLEX(8)   :: PC(0:N1,0:N2), RHS(N1,N2)
       COMPLEX(8)   :: RESD(N1MD,N2M), GGII(0:N1MD,0:N2)

       KK = 0
       RESD = 0.0_8
       GGII = 0.0_8

       CALL TOPLEVEL(0_8, PC, RHS, SUMRES, OLDV, TEST, KV, RESD, GGII)
       IF(SUMRES < TEST) GOTO 205

       DO KK = 1, MGITR            ! main iteration
           DO ILEV = NLEV - 1, 1, -1 
               CALL RELAX(ILEV, 0.0_8, 1_8, KV, RESD, GGII) 
               CALL GODOWN(ILEV, KV, RESD, GGII)
           END DO

           CALL RELAX(0_8, 0.0_8, NBLI, KV, RESD, GGII)

           DO ILEV = 0, NLEV - 2
               CALL GOUP(ILEV, RESD, GGII)
               CALL RELAX(ILEV + 1, 1.0_8, 1_8, KV, RESD, GGII)
           END DO

           CALL TOPLEVEL(1_8, PC, RHS, SUMRES, 1.0_8, TEST, KV, RESD, GGII)

           IF(SUMRES < TEST) GOTO 205
       END DO

       PRINT *, 'ITERATION LIMIT EXCEEDED.'

205    CONTINUE
       IF (KV <= 3) THEN
           WRITE(77, 999) KV, KK, SUMRES * DTCONST
       END IF

999    FORMAT('KV=', I4, 3X, 'KK=', I4, 3X, 'RM=', ES24.16)

      RETURN
      END SUBROUTINE MG2D

!=======================================================================
      SUBROUTINE TOPLEVEL(ID, PC, RHS, SUMRES, OLDV, TEST, KV, RESD, GGII)    ! Zebra Version
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
        INTEGER(8)  :: I, J, IC
        INTEGER(8)  :: ID, KV
        REAL(8)     :: SUMRES, TEST, OLDV
        COMPLEX(8)  :: TT
        COMPLEX(8)  :: PC(0:N1,0:N2), RHS(N1,N2)
        COMPLEX(8)  :: RESD(N1MD,N2M), GGII(0:N1MD,0:N2)

      IF (ID == 1) THEN
!         INTERPOLATE & ADD
          IC = 0
          DO I = IIMG(NLEV - 1, 1), IIMG(NLEV - 1, 2) - 1
              IC = IC + 2
              DO J = 1, N2M
                  PC(IC, J) = PC(IC, J) + COI1(I) * GGII(I, J) + COI2(I) * GGII(IPM(I), J)
              END DO
          END DO

          I = IIMG(NLEV - 1, 2)
          IC = IC + 2
          DO J = 1, N2M
              PC(IC, J) = PC(IC, J) + COI1(I) * GGII(I, J)      ! use one point COI(~)=1.
          END DO
      END IF

!  RELAX
      DO J = 1, N2M
          DO I = 1, N1M, 2
              GGII(I, J) = RHS(I, J) - OLDV * (AW(I) * PC(IMM(I), J) + AE(I) * PC(IPM(I), J))
          END DO
      END DO

      CALL TRDIAG1M(GGII, GGII, 1_8, N1M, 1_8, N2M, KV)

      IF (KV == 1) THEN
          TT = GGII(1, N2M - 1)
      ELSE
          TT = 0.0_8
      END IF

      DO J = 1, N2M
          DO I = 2, N1M, 2
              GGII(I - 1, J) = GGII(I - 1, J) - TT
              GGII(I, J) = RHS(I, J) - AW(I) * GGII(IMM(I), J) - AE(I) * (GGII(IPM(I), J) - TT)
          END DO
      END DO

      CALL TRDIAG1M(GGII, GGII, 2_8, N1M, 1_8, N2M, KV)

!  CALCULATE RESIDUAL
      SUMRES = 0.0_8

      DO J = 1, N2M
          DO I = 1, N1M
              RESD(I, J) = RHS(I, J) - AW(I) * GGII(IMM(I), J) - AE(I) * GGII(IPM(I), J)  &
                           - AS(J) * GGII(I, J - 1) - AN(J) * GGII(I, J + 1)              &
                           - AC(I, J, KV) * GGII(I, J)
              SUMRES = MAX(SUMRES, ABS(RESD(I, J)))
              PC(I, J) = GGII(I, J)
          END DO
      END DO

      IF(SUMRES < TEST .OR. ID == 2) RETURN

!  RESTRICT
      IC = -1
      DO I = IIMG(NLEV - 1, 1), IIMG(NLEV - 1, 2)
          IC = IC + 2
          DO J = 1, N2M
              RESD(I, J) = RESD(IC, J) * COR1(I) + RESD(IC + 1, J) * COR2(I)
          END DO
      END DO

      RETURN
      END SUBROUTINE TOPLEVEL

!=======================================================================
      SUBROUTINE TRDIAG1M(RR, UU, L1, L2, LL1, LL2, KV)
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
        INTEGER(8)  :: I, J, L1, L2, LL1, LL2, KV
        COMPLEX(8)  :: RR(0:N1MD,0:N2), UU(0:N1MD,0:N2)

      DO I = L1, L2, 2
          UU(I, LL1) = RR(I, LL1) * BET(I, 1, KV)
      END DO

      DO J = LL1 + 1, LL2
          DO I = L1, L2, 2
              UU(I, J) = (RR(I, J) - AS(J) * UU(I, J - 1)) * BET(I, J, KV)
          END DO
      END DO
      
      DO J = LL2 - 1, LL1, -1
          DO I = L1, L2, 2
              UU(I, J) = UU(I, J) - GAM(I, J + 1, KV) * UU(I, J + 1)
          END DO
      END DO

      RETURN
      END SUBROUTINE TRDIAG1M

!=======================================================================
      SUBROUTINE GOUP(ILEV, RESD, GGII)! INTERPOLATE RESIDUAL & ADD IT TO HIGH LEVEL
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
        INTEGER(8)  :: I, J, ILEV, IBGH
        COMPLEX(8)  :: RESD(N1MD,N2M), GGII(0:N1MD,0:N2)

      IBGH = IIMG(ILEV + 1, 1) - 1

      DO I = IIMG(ILEV, 1), IIMG(ILEV, 2) - 1
          IBGH = IBGH + 2
          DO J = 1, N2M
              GGII(IBGH, J) = GGII(IBGH, J) + COI1(I) * GGII(I, J) + COI2(I) * GGII(IPM(I), J)
          END DO
      END DO

      I = IIMG(ILEV, 2)
      IBGH = IBGH + 2
      DO J = 1, N2M
          GGII(IBGH, J) = GGII(IBGH, J) + COI1(I) * GGII(I, J)           ! use one point
      END DO

      RETURN
      END SUBROUTINE GOUP

!=======================================================================
      SUBROUTINE GODOWN(ILEV, KV, RESD, GGII)        ! COMPUTE RESIDUAL & RESTRICT IT
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
        INTEGER(8)  :: I, J, ILEV, KV, IBG
        COMPLEX(8)  :: RESD(N1MD,N2M), GGII(0:N1MD,0:N2)

      DO J = 1, N2M
          DO I = IIMG(ILEV, 1), IIMG(ILEV, 2)
              RESD(I, J) = RESD(I, J) - AW(I) * GGII(IMM(I), J) - AE(I) * GGII(IPM(I), J)   &
                           - AS(J) * GGII(I, J - 1) - AN(J) * GGII(I, J + 1) - AC(I, J, KV) * GGII(I, J)
          END DO
      END DO

      IBG = IIMG(ILEV, 1) - 2

      DO I = IIMG(ILEV - 1, 1), IIMG(ILEV - 1, 2)
          IBG = IBG + 2
          DO J = 1, N2M
              RESD(I, J) = RESD(IBG, J) * COR1(I) + RESD(IBG + 1, J) * COR2(I)
          END DO
      END DO

      RETURN
      END SUBROUTINE GODOWN

!=======================================================================
      SUBROUTINE RELAX(ILEV, OLDV, IITER, KV, RESD, GGII)   ! Zebra Version
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
        INTEGER(8)  :: I, J, ILEV, IITER, KV, II
        COMPLEX(8)  :: RESD(N1MD,N2M), GGII(0:N1MD,0:N2)
        REAL(8)     :: OLDV

      DO J = 1, N2M
          DO I = IIMG(ILEV, 1), IIMG(ILEV, 2), 2
              GGII(I, J) = RESD(I, J) - OLDV * (AW(I) * GGII(IMM(I), J) + AE(I) * GGII(IPM(I), J))
          END DO
      END DO

      CALL TRDIAG1M(GGII, GGII, IIMG(ILEV, 1), IIMG(ILEV, 2), 1_8, N2M, KV)

      DO J = 1, N2M
          DO I = IIMG(ILEV, 1) + 1, IIMG(ILEV, 2), 2
              GGII(I, J) = RESD(I, J) - AW(I) * GGII(IMM(I), J) - AE(I) * GGII(IPM(I), J)
          END DO
      END DO

      CALL TRDIAG1M(GGII, GGII, IIMG(ILEV, 1) + 1, IIMG(ILEV, 2), 1_8, N2M, KV)

      DO II = 1, IITER - 1

          DO J = 1, N2M
              DO I = IIMG(ILEV, 1), IIMG(ILEV, 2), 2
                  GGII(I, J) = RESD(I, J) - (AW(I) * GGII(IMM(I), J) + AE(I) * GGII(IPM(I), J))
              END DO
          END DO

          CALL TRDIAG1M(GGII, GGII, IIMG(ILEV, 1), IIMG(ILEV, 2), 1_8, N2M, KV)

          DO J = 1, N2M
              DO I = IIMG(ILEV, 1) + 1, IIMG(ILEV, 2), 2
                  GGII(I, J) = RESD(I, J) - AW(I) * GGII(IMM(I), J) - AE(I) * GGII(IPM(I), J)
              END DO
          END DO

          CALL TRDIAG1M(GGII, GGII, IIMG(ILEV, 1) + 1, IIMG(ILEV, 2), 1_8, N2M, KV)

      END DO
      RETURN
      END SUBROUTINE RELAX
!=======================================================================
      SUBROUTINE GSOR2D(U, RHS, KV, TEST, OLDV)    ! 1 EQ. TYPE
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8)   :: I, J, KV
       REAL(8)      :: TEST, OLDV
       COMPLEX(8)   :: U(0:N1,0:N2), RHS(N1,N2)
       COMPLEX(8)   :: GGII(0:N1MD,0:N2)
       INTEGER(8)   :: KK
       REAL(8)      :: WW, WW2, ERRMAX


      WW  = WWSOR
      WW2 = 1.0_8 - WW
      KK  = 0

!  HALF RELAX
!  ----------

      DO J = 1, N2M
          DO I = 1, N1M, 2
              GGII(I, J) = RHS(I, J) - OLDV * (AW(I) * U(IMM(I), J) + AE(I) * U(IPM(I), J))
          END DO
      END DO

      CALL TRDIAG1M(GGII, GGII, 1_8, N1M, 1_8, N2M, KV)

      DO J = 1, N2M
          DO I = 1, N1M, 2
              U(I, J) = WW * GGII(I, J) + OLDV * WW2 * U(I, J)
          END DO
      END DO

!  ANOTHER HALF
!  ------------


      DO J = 1, N2M
          DO I = 2, N1M, 2
              GGII(I, J) = RHS(I, J) - (AW(I) * U(IMM(I), J) + AE(I) * U(IPM(I), J))
          END DO
      END DO

      CALL TRDIAG1M(GGII, GGII, 2_8, N1M, 1_8, N2M, KV)

      DO J = 1, N2M
          DO I = 2, N1M, 2
              U(I, J) = WW * GGII(I, J) + OLDV * WW2 * U(I, J)
          END DO
      END DO

      CALL RESID3(U, RHS, KV, ERRMAX)
      IF(ERRMAX < TEST) GOTO 88

!  MAIN ITERATION
!  ==============
      DO KK = 1, MGITR

!  HALF RELAX
!  ----------

          DO J = 1, N2M
              DO I = 1, N1M, 2
                  GGII(I, J) = RHS(I, J) - (AW(I) * U(IMM(I), J) + AE(I) * U(IPM(I), J))
              END DO
          END DO

          CALL TRDIAG1M(GGII, GGII, 1_8, N1M, 1_8, N2M, KV)

          DO J = 1, N2M
              DO I = 1, N1M, 2
                  U(I, J) = WW * GGII(I, J) + WW2 * U(I, J)
              END DO
          END DO

!  ANOTHER HALF
!  ------------

          DO J = 1, N2M
              DO I = 2, N1M, 2
                  GGII(I, J) = RHS(I, J) - (AW(I) * U(IMM(I), J) + AE(I) * U(IPM(I), J))
              END DO
          END DO

          CALL TRDIAG1M(GGII, GGII, 2_8, N1M, 1_8, N2M, KV)

          DO J = 1, N2M
              DO I = 2, N1M, 2
                  U(I, J) = WW * GGII(I, J) + WW2 * U(I, J)
              END DO
          END DO

          CALL RESID3(U, RHS, KV, ERRMAX)
          WRITE(77, 999) KV, KK, ERRMAX * DTCONST * FLOAT(N3MH)
          IF(ERRMAX < TEST) GOTO 88

      END DO

      PRINT *,'ITERATION LIMIT EXCEEDED.'
88    CONTINUE

999   FORMAT('KV=', I3, 3X, 'KK=', I3, 3X, 'RG=', ES24.16)

      RETURN
      END SUBROUTINE GSOR2D

!=======================================================================
      SUBROUTINE RESID3(U, RHS, KV, ERRMAX)
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8)   :: I, J, KV
       REAL(8)      :: ERRMAX
       COMPLEX(8)   :: U(0:N1,0:N2), RHS(N1,N2)
       COMPLEX(8)   :: ERR

       ERRMAX = 0.0_8

      DO J = 1, N2M
          DO I = 1, N1M
              ERR = RHS(I, J) - AW(I) * U(IMM(I), J) - AE(I) * U(IPM(I), J)   &
                              - AS(J) * U(I, J - 1) - AN(J) * U(I, J + 1)     &
                              - AC(I, J, KV) * U(I, J)
              ERRMAX = MAX(ERRMAX, ABS(ERR))
          END DO
      END DO

      RETURN
      END SUBROUTINE RESID3
!=======================================================================
      SUBROUTINE COEFMG
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8)   :: I, J, K, ILEV
       INTEGER(8)   :: MINROW, IBG, IEND, ISP, IC, IBGH
       REAL(8)      :: IWZ(N1MD), IEZ(N1MD)
       REAL(8)      :: XMPM(0:N1MD)
       REAL(8)      :: VDX_IBG, VDX_IEND, SDX_IBG, SDX_IEND

       MINROW = N1M / (2**NLEV)
       XMPM = 0.0_8

       LEVHALF = NINT(NLEV / 2.0_8)
       IIMG(NLEV, 1) = 1
       IIMG(NLEV, 2) = N1M
       IIMG(NLEV, 3) = N1M

       DO I = 1, N1M
           XMPM(I) = XMP(I)
           IWZ(I) = 1 - 1/I                  ! 0 only if I=1
           IEZ(I) = 1 - I/N1M                ! 0 only if I=N1M
       END DO

!  COMPUTE FOR LOWER LEVELS
      DO ILEV = NLEV - 1, 0, -1

          IIMG(ILEV, 1) = IIMG(ILEV + 1, 1) + IIMG(ILEV + 1, 3)
          IIMG(ILEV, 3) = MINROW * (2**ILEV)
          IIMG(ILEV, 2) = IIMG(ILEV, 1) + IIMG(ILEV, 3) - 1

          IBG  = IIMG(ILEV, 1)
          IEND = IIMG(ILEV, 2)
          ISP  = 2**(NLEV - ILEV)           ! width of one cell at low level

          IC = 0
          DO I = IBG, IEND
              IC = IC + 1
              XMPM(I) = 0.5_8 * (X(IC * ISP + 1) + X((IC - 1) * ISP + 1))
              IWZ(I)  = 1 - IBG/I                ! 0 only if I=IBG
              IEZ(I)  = 1 - I/IEND               ! 0 only if I=IEND
          END DO

          IC = 0
          DO I = IBG, IEND
              IC = IC + 1
              AW(I) = IWZ(I) / ((XMPM(I) - XMPM(I - 1)) * (X(IC * ISP + 1) - X((IC - 1) * ISP + 1)))
              AE(I) = IEZ(I) / ((XMPM(I + 1) - XMPM(I)) * (X(IC * ISP + 1) - X((IC - 1) * ISP + 1)))
          END DO

          DO J = 1, N2M
              DO I = IIMG(ILEV, 1), IIMG(ILEV, 2)
                  AC(I, J, 1) = -1.0_8 * (AW(I) + AE(I) + AS(J) + AN(J))
              END DO
          END DO

          DO K = 2, N3MH
              DO J = 1, N2M
                  DO I = IIMG(ILEV, 1), IIMG(ILEV, 2)
                      AC(I, J, K) = AC(I, J, 1) - AK3(K)
                  END DO
              END DO
          END DO

      END DO

!  CALCULATE RESTRICTION COEFFS
       DO ILEV = NLEV, 1, -1
           IBGH = IIMG(ILEV, 1)
           DO I = IIMG(ILEV - 1, 1), IIMG(ILEV - 1, 2)
               COR1(I) = (XMPM(IBGH + 1) - XMPM(I)) / (XMPM(IBGH + 1) - XMPM(IBGH))
               COR2(I) = 1.0_8 - COR1(I)
               IBGH = IBGH + 2
           END DO
       END DO

!  CALCULATE INTERPOLATION COEFFS
       DO ILEV = 0, NLEV - 1
           IBGH = IIMG(ILEV + 1, 1) + 1
           DO I = IIMG(ILEV, 1), IIMG(ILEV, 2) - 1
               COI1(I) = (XMPM(I + 1) - XMPM(IBGH)) / (XMPM(I + 1) - XMPM(I))  ! * lower value
               COI2(I) = 1.0_8 - COI1(I)
               IBGH = IBGH + 2
           END DO
           I = IIMG(ILEV, 2)
           COI1(I) = 1.0_8                 ! use only one lower point at upper wall
       END DO

!===== FOR THE X PERIODICITY
!       INTRODUCE IPM & IMM
       ! Ensure safe default initialization for non-periodic boundaries as well
       DO ILEV = NLEV, 0, -1
           IBG  = IIMG(ILEV, 1)
           IEND = IIMG(ILEV, 2)
           DO I = IBG, IEND
               IPM(I) = I + 1
               IMM(I) = I - 1
           END DO
       END DO

       IF (XPRDIC == 1) THEN
           DO ILEV = NLEV, 0, -1
               IBG  = IIMG(ILEV, 1)
               IEND = IIMG(ILEV, 2)
               IPM(IEND) = IBG
               IMM(IBG)  = IEND
           END DO

           DO ILEV = NLEV - 1, 0, -1
               IBG  = IIMG(ILEV, 1)
               IEND = IIMG(ILEV, 2)
               ISP  = 2**(NLEV - ILEV)
               VDX_IBG  = XMPM(IBG) - X(1) + X(N1) - XMPM(IEND)
               VDX_IEND = XMPM(IBG) - X(1) + X(N1) - XMPM(IEND)
               SDX_IBG  = X(1 + ISP) - X(1)
               SDX_IEND = X(N1) - X(N1 - ISP)
               AW(IBG)  = 1.0_8 / (VDX_IBG * SDX_IBG)
               AE(IEND) = 1.0_8 / (VDX_IEND * SDX_IEND)
           END DO

           DO ILEV = NLEV - 1, 0, -1
               DO J = 1, N2M
                   DO I = IIMG(ILEV, 1), IIMG(ILEV, 2)
                       AC(I, J, 1) = -1.0_8 * (AW(I) + AE(I) + AS(J) + AN(J))
                   END DO
               END DO

               DO K = 2, N3MH
                   DO J = 1, N2M
                       DO I = IIMG(ILEV, 1), IIMG(ILEV, 2)
                           AC(I, J, K) = AC(I, J, 1) - AK3(K)
                       END DO
                   END DO
               END DO
           END DO

!  CALCULATE INTERPOLATION COEFFS
           DO ILEV = 0, NLEV - 1
               IBG  = IIMG(ILEV, 1)
               IEND = IIMG(ILEV, 2)
               VDX_IEND = XMPM(IBG) - X(1) + X(N1) - XMPM(IEND)
               COI2(IEND) = (XMPM(IIMG(ILEV + 1, 2)) - XMPM(IEND)) / VDX_IEND
               COI1(IEND) = 1.0_8 - COI2(IEND)
           END DO

       END IF

      DO ILEV = NLEV, 0, -1
          WRITE(*,*) 'IIMG(1', ILEV, ')=', IIMG(ILEV, 1)
          WRITE(*,*) 'IIMG(2', ILEV, ')=', IIMG(ILEV, 2)
          WRITE(*,*) 'IIMG(3', ILEV, ')=', IIMG(ILEV, 3)
      END DO

      RETURN
      END SUBROUTINE COEFMG