!=======================================================================
!
!     MAIN SOLVER OF POISSON EQUATION
!
!     AX=b
!     A: coefficient of discretized poisson equation
!     X: PHI: Pseudo pressure, Output of subroutine POISSON
!     b: DIVGSUM: OUTPUT of subroutine DIVGS.
!
!     x-direction: Fourier transform
!     y-direction: TDMA
!     z-direction: MULTI-GRID iteration/GSOR method
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

       CALL X_FT_ALLO

       ALLOCATE(AC(N3MD, N2M, N1MH))
       ALLOCATE(GAM(N3MD, N2M, N1MH))
       ALLOCATE(BET(N3MD, N2M, N1MH))
 
       CALL PMAT
!------MULTIGRID METHOD
       CALL COEFMG

       DO ILEV = 0, NLEV
           DO I = 1, N1MH
               DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
                   BET(K, 1, I) = 1.0_8 / AC(K, 1, I)
               END DO
               DO J = 2, N2M
                   DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
                       GAM(K, J, I) = AN(J - 1) * BET(K, J - 1, I)
                       BET(K, J, I) = 1.0_8 / (AC(K, J, I) - AS(J) * GAM(K, J, I))
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
       INTEGER(8)  :: I, J, K, JP, KP

!-----FOR TOP LEVEL
      DO K = 1, N3M
          KP = KPV(K)
          AB(K) = (1.0_8 - FIXKL(K)) * C2CZI(K) * F2FZI(K)
          AF(K) = (1.0_8 - FIXKU(K)) * C2CZI(KP) * F2FZI(K)
      END DO

      DO J = 1, N2M
          JP = JPV(J)
          AS(J) = (1.0_8 - FIXJL(J)) * C2CYI(J) * F2FYI(J)
          AN(J) = (1.0_8 - FIXJU(J)) * C2CYI(JP) * F2FYI(J)
      END DO

      DO J = 1, N2M
          DO K = 1, N3M
              AC(K, J, 1) = -1.0_8 * (AB(K) + AF(K) + AS(J) + AN(J))
          END DO
      END DO

      N1MH = N1M / 2 + 1

      CALL MWAVENUMBER     ! INIT. MODIFIED WAVE #.

      DO I = 2, N1MH
          DO J = 1, N2M
              DO K = 1, N3M
                  AC(K, J, I) = AC(K, J, 1) - AI3(I)
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
       INTEGER(8)  :: I
       REAL(8)     :: PI

       PI = ACOS(-1.0_8)

      DO I = 1, N1MH
          AI3(I) = 2.0_8 * (1.0_8 - COS(2.0_8 * PI * FLOAT(I - 1) / FLOAT(N1M))) * F2FXI(1) * F2FXI(1)
      END DO

      RETURN
      END SUBROUTINE MWAVENUMBER

!=======================================================================
      SUBROUTINE POISSON(PHI, DIVGSUM)
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8)  :: I, J, K, III
       REAL(8)     :: PHI(0:N1,0:N2,0:N3), DIVGSUM(0:N1,0:N2,0:N3)
       COMPLEX(8)  :: CCAP(N3,N2,N1MH)
       COMPLEX(8), DIMENSION(:,:), ALLOCATABLE :: XXX, CP
       COMPLEX(8), DIMENSION(:),   ALLOCATABLE :: XXXX, XXXX_B

       REAL(8)     :: TEST, PHIREF

! --- DO THE FORWARD FFT
!$OMP PARALLEL private(XXX, XXXX, XXXX_B)
      ALLOCATE(XXX(N1M, N3M))
      ALLOCATE(XXXX(N1M), XXXX_B(N1M*2))
      CALL ZFFT1D(XXXX, N1M, 0, XXXX_B)
!$OMP DO
      DO J = 1, N2M

          DO I = 1, N1M
              DO K = 1, N3M
                  XXX(I, K) = DIVGSUM(I, J, K)
              END DO
          END DO

          DO K = 1, N3M
              CALL ZFFT1D(XXX(1, K), N1M, -1, XXXX_B)
          END DO

          DO I = 1, N1MH
              DO K = 1, N3M
                  CCAP(K, J, I) = XXX(I, K)
              END DO
          END DO

      END DO
!$OMP END DO
      DEALLOCATE(XXX, XXXX, XXXX_B)
!$OMP END PARALLEL

! --- SOLVE A SET OF POISSON EQS.
       TEST = TEST1 / FLOAT(N1MH) * FLOAT(N1M) * 0.9_8

!$OMP PARALLEL private(CP)
      ALLOCATE(CP(0:N3, 0:N2))
!$OMP DO
      DO III = 1, N1MH
          CP = 0.0_8

          IF (III <= IMGSOR) THEN
              CALL MG2D(CP, CCAP(1, 1, III), III, TEST, 0.0_8)
          ELSE
              CALL GSOR2D(CP, CCAP(1, 1, III), III, TEST, 0.0_8)
          END IF
        
          DO J = 1, N2M
              DO K = 1, N3M
                  CCAP(K, J, III) = CP(K, J)
              END DO
          END DO
      END DO
!$OMP END DO
      DEALLOCATE(CP)
!$OMP END PARALLEL
 
! --- DO THE INVERSE FFT
!$OMP PARALLEL private(XXX, XXXX, XXXX_B)
      ALLOCATE(XXX(N1M, N3M))
      ALLOCATE(XXXX(N1M), XXXX_B(N1M*2))
      CALL ZFFT1D(XXXX, N1M, 0, XXXX_B)
!$OMP DO
      DO J = 1, N2M
          DO I = 1, N1MH
              DO K = 1, N3M
                  XXX(I, K) = CCAP(K, J, I)
              END DO
          END DO
        
          DO I = N1MH + 1, N1M
              DO K = 1, N3M  
                  XXX(I, K) = CONJG(CCAP(K, J, N1M + 2 - I))
              END DO
          END DO
        
          DO K = 1, N3M
              CALL ZFFT1D(XXX(1, K), N1M, 1, XXXX_B)
          END DO

          DO I = 1, N1M        
              DO K = 1, N3M
                  PHI(I, J, K) = REAL(XXX(I, K), 8)
              END DO
          END DO
      END DO
!$OMP END DO
      DEALLOCATE(XXX, XXXX, XXXX_B)
!$OMP END PARALLEL

      IF (ICH == 1) THEN
!     SET THE AVERAGE PHI AT THE UPPER WALL TO BE ZERO.
          PHIREF = 0.0_8
!$OMP PARALLEL DO reduction(+:PHIREF)
          DO I = 1, N1M
              PHIREF = PHIREF + PHI(I, N2M, 1) * F2FX(I)
          END DO
          PHIREF = PHIREF / XL

!$OMP PARALLEL DO
          DO K = 1, N3M
              DO J = 1, N2M
                  DO I = 1, N1M
                      PHI(I, J, K) = PHI(I, J, K) - PHIREF
                  END DO
              END DO
          END DO
      END IF

      RETURN
      END SUBROUTINE POISSON

!=======================================================================
       SUBROUTINE MG2D(PC, RHS, IV, TEST, OLDV)
!=======================================================================
!     IV      : WAVE NUMBER INDEX
!     MULTIGRID ENVIRONMENT VARIABLES
!     TEST    : CONDITION FOR CONVERGENCE
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
       INTEGER(8)  :: K, J, IV, II, ILEV
       REAL(8)     :: TEST, SUMRES, OLDV
       COMPLEX(8)  :: PC(0:N3,0:N2), RHS(N3,N2)
       COMPLEX(8)  :: RESD(N3MD,N2M), GGII(0:N3MD,0:N2)

       II = 0
       RESD = 0.0_8
       GGII = 0.0_8

       CALL TOPLEVEL(0_8, PC, RHS, SUMRES, OLDV, TEST, IV, RESD, GGII)
       IF (SUMRES < TEST) GOTO 2000

       DO II = 1, MGITR            ! main iteration
           DO ILEV = NLEV - 1, 1, -1
               CALL RELAX(ILEV, 0.0_8, 1_8, IV, RESD, GGII)
               CALL GODOWN(ILEV, IV, RESD, GGII)
           END DO

           CALL RELAX(0_8, 0.0_8, NBLI, IV, RESD, GGII)

           DO ILEV = 0, NLEV - 2
               CALL GOUP(ILEV, RESD, GGII)
               CALL RELAX(ILEV + 1, 1.0_8, 1_8, IV, RESD, GGII)
           END DO

           CALL TOPLEVEL(1_8, PC, RHS, SUMRES, 1.0_8, TEST, IV, RESD, GGII)

           IF (SUMRES < TEST) GOTO 2000
       END DO
       
       WRITE(*,*) 'ITERATION LIMIT EXCEEDED.'

 2000  CONTINUE

       IF (IV <= 2) WRITE(77, 201) ' MG, TIME = ', TIME, IV, II, SUMRES * DTCONST * FLOAT(N1MH) / FLOAT(N1M)

  201  FORMAT(A13, F13.5, ' IV=', I5, ' II=', I5, ' RG=', ES16.8)

      RETURN
      END SUBROUTINE MG2D
!=======================================================================
      SUBROUTINE TOPLEVEL(ID, PC, RHS, SUMRES, OLDV, TEST, IV, RESD, GGII)    ! Zebra Version
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
        INTEGER(8)  :: K, J, KC
        INTEGER(8)  :: ID, IV
        REAL(8)     :: SUMRES, TEST, OLDV
        COMPLEX(8)  :: TT
        COMPLEX(8)  :: PC(0:N3,0:N2), RHS(N3,N2)
        COMPLEX(8)  :: RESD(N3MD,N2M), GGII(0:N3MD,0:N2)

      IF (ID == 1) THEN
!         INTERPOLATE & ADD
          KC = 0
          DO K = KKMG(NLEV - 1, 1), KKMG(NLEV - 1, 2)
              KC = KC + 2
              DO J = 1, N2M
                  PC(KC, J) = PC(KC, J) + COI1(K) * GGII(K, J) + COI2(K) * GGII(KPM(K), J)
              END DO
          END DO
      END IF

!  RELAX
      DO J = 1, N2M
          DO K = 1, N3M, 2
              GGII(K, J) = RHS(K, J) - OLDV * (AB(K) * PC(KMM(K), J) + AF(K) * PC(KPM(K), J))
          END DO
      END DO

      CALL TRDIAG1M(GGII, GGII, 1_8, N3M, 1_8, N2M, IV)

      IF (IV == 1) THEN
          TT = GGII(1, N2M)
          DO J = 1, N2M
              DO K = 1, N3M, 2
                  GGII(K, J) = GGII(K, J) - TT
              END DO
          END DO
      END IF

      DO J = 1, N2M
          DO K = 2, N3M, 2
              GGII(K, J) = RHS(K, J) - (AB(K) * GGII(KMM(K), J) + AF(K) * (GGII(KPM(K), J)))
          END DO
      END DO

      CALL TRDIAG1M(GGII, GGII, 2_8, N3M, 1_8, N2M, IV)

!  CALCULATE RESIDUAL
      SUMRES = 0.0_8

      DO J = 1, N2M
          DO K = 1, N3M
              RESD(K, J) = RHS(K, J) - AB(K) * GGII(KMM(K), J) - AF(K) * GGII(KPM(K), J)  &
                           - AS(J) * GGII(K, J - 1) - AN(J) * GGII(K, J + 1)              &
                           - AC(K, J, IV) * GGII(K, J)
              SUMRES = MAX(SUMRES, ABS(RESD(K, J)))
              PC(K, J) = GGII(K, J)
          END DO
      END DO

      IF (SUMRES < TEST) RETURN
      IF (ID == 2) RETURN

!  RESTRICT
      KC = -1
      DO K = KKMG(NLEV - 1, 1), KKMG(NLEV - 1, 2)
          KC = KC + 2
          DO J = 1, N2M
              RESD(K, J) = RESD(KC, J) * COR1(K) + RESD(KC + 1, J) * COR2(K)
          END DO
      END DO

      RETURN
      END SUBROUTINE TOPLEVEL
!=======================================================================
      SUBROUTINE TRDIAG1M(RR, UU, L1, L2, LL1, LL2, IV)
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
        INTEGER(8)  :: K, J, L1, L2, LL1, LL2, IV
        COMPLEX(8)  :: RR(0:N3MD,0:N2), UU(0:N3MD,0:N2)

      DO K = L1, L2, 2
          UU(K, LL1) = RR(K, LL1) * BET(K, 1, IV)
      END DO

      DO J = LL1 + 1, LL2
          DO K = L1, L2, 2
              UU(K, J) = (RR(K, J) - AS(J) * UU(K, J - 1)) * BET(K, J, IV)
          END DO
      END DO

      DO J = LL2 - 1, LL1, -1
          DO K = L1, L2, 2
              UU(K, J) = UU(K, J) - GAM(K, J + 1, IV) * UU(K, J + 1)
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
        INTEGER(8)  :: K, J, ILEV, KBGH
        COMPLEX(8)  :: RESD(N3MD,N2M), GGII(0:N3MD,0:N2)

      KBGH = KKMG(ILEV + 1, 1) - 1

      DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
          KBGH = KBGH + 2
          DO J = 1, N2M
              GGII(KBGH, J) = GGII(KBGH, J) + COI1(K) * GGII(K, J) + COI2(K) * GGII(KPM(K), J)
          END DO
      END DO

      RETURN
      END SUBROUTINE GOUP

!=======================================================================
      SUBROUTINE GODOWN(ILEV, IV, RESD, GGII)        ! COMPUTE RESIDUAL & RESTRICT IT
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
        INTEGER(8)  :: K, J, ILEV, IV, KBG
        COMPLEX(8)  :: RESD(N3MD,N2M), GGII(0:N3MD,0:N2)

      DO J = 1, N2M
          DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
              RESD(K, J) = RESD(K, J) - AB(K) * GGII(KMM(K), J) - AF(K) * GGII(KPM(K), J)   &
                           - AS(J) * GGII(K, J - 1) - AN(J) * GGII(K, J + 1) - AC(K, J, IV) * GGII(K, J)
          END DO
      END DO

      KBG = KKMG(ILEV, 1) - 2

      DO K = KKMG(ILEV - 1, 1), KKMG(ILEV - 1, 2)
          KBG = KBG + 2
          DO J = 1, N2M
              RESD(K, J) = RESD(KBG, J) * COR1(K) + RESD(KBG + 1, J) * COR2(K)
          END DO
      END DO

      RETURN
      END SUBROUTINE GODOWN

!=======================================================================
      SUBROUTINE RELAX(ILEV, OLDV, IITER, IV, RESD, GGII)   ! Zebra Version
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
        INTEGER(8)  :: K, J, ILEV, IITER, IV, KK
        COMPLEX(8)  :: RESD(N3MD,N2M), GGII(0:N3MD,0:N2)
        REAL(8)     :: OLDV

      DO J = 1, N2M
          DO K = KKMG(ILEV, 1), KKMG(ILEV, 2), 2
              GGII(K, J) = RESD(K, J) - OLDV * (AB(K) * GGII(KMM(K), J) + AF(K) * GGII(KPM(K), J))
          END DO
      END DO

      CALL TRDIAG1M(GGII, GGII, KKMG(ILEV, 1), KKMG(ILEV, 2), 1_8, N2M, IV)

      DO J = 1, N2M
          DO K = KKMG(ILEV, 1) + 1, KKMG(ILEV, 2), 2
              GGII(K, J) = RESD(K, J) - AB(K) * GGII(KMM(K), J) - AF(K) * GGII(KPM(K), J)
          END DO
      END DO

      CALL TRDIAG1M(GGII, GGII, KKMG(ILEV, 1) + 1, KKMG(ILEV, 2), 1_8, N2M, IV)

      DO KK = 1, IITER - 1

          DO J = 1, N2M
              DO K = KKMG(ILEV, 1), KKMG(ILEV, 2), 2
                  GGII(K, J) = RESD(K, J) - (AB(K) * GGII(KMM(K), J) + AF(K) * GGII(KPM(K), J))
              END DO
          END DO

          CALL TRDIAG1M(GGII, GGII, KKMG(ILEV, 1), KKMG(ILEV, 2), 1_8, N2M, IV)

          DO J = 1, N2M
              DO K = KKMG(ILEV, 1) + 1, KKMG(ILEV, 2), 2
                  GGII(K, J) = RESD(K, J) - AB(K) * GGII(KMM(K), J) - AF(K) * GGII(KPM(K), J)
              END DO
          END DO

          CALL TRDIAG1M(GGII, GGII, KKMG(ILEV, 1) + 1, KKMG(ILEV, 2), 1_8, N2M, IV)

      END DO
      RETURN
      END SUBROUTINE RELAX
!=======================================================================
      SUBROUTINE GSOR2D(U, RHS, IV, TEST, OLDV)    ! 1 EQ. TYPE
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8)   :: K, J, IV
       REAL(8)      :: TEST, OLDV
       COMPLEX(8)   :: U(0:N3,0:N2), RHS(N3,N2)
       COMPLEX(8)   :: GGII(0:N3MD,0:N2)
       COMPLEX(8)   :: TT
       INTEGER(8)   :: II
       REAL(8)      :: WW, WW2, ERRMAX

      GGII = 0.0_8
      WW   = WWSOR
      WW2  = 1.0_8 - WW
      II   = 0

!  HALF RELAX
!  ----------
      DO J = 1, N2M
          DO K = 1, N3M, 2
              GGII(K, J) = RHS(K, J) - OLDV * (AB(K) * U(KMM(K), J) + AF(K) * U(KPM(K), J))
          END DO
      END DO

      CALL TRDIAG1M(GGII, GGII, 1_8, N3M, 1_8, N2M, IV)

      DO J = 1, N2M
          DO K = 1, N3M, 2
              U(K, J) = WW * GGII(K, J) + OLDV * WW2 * U(K, J)
          END DO
      END DO

!  ANOTHER HALF
!  ------------
      IF (IV == 1) THEN
          TT = U(1, N2M)
          DO J = 1, N2M
              DO K = 1, N3M, 2
                  U(K, J) = U(K, J) - TT
              END DO
          END DO
      END IF

      DO J = 1, N2M
          DO K = 2, N3M, 2
              GGII(K, J) = RHS(K, J) - (AB(K) * U(KMM(K), J) + AF(K) * U(KPM(K), J))
          END DO
      END DO

      CALL TRDIAG1M(GGII, GGII, 2_8, N3M, 1_8, N2M, IV)

      DO J = 1, N2M
          DO K = 2, N3M, 2
              U(K, J) = WW * GGII(K, J) + OLDV * WW2 * U(K, J)
          END DO
      END DO

      CALL RESID3(U, RHS, IV, ERRMAX)
      IF (ERRMAX < TEST) GOTO 1000

!  MAIN ITERATION
!  ==============
      DO II = 1, MGITR

!  HALF RELAX
!  ----------
          DO J = 1, N2M
              DO K = 1, N3M, 2
                  GGII(K, J) = RHS(K, J) - (AB(K) * U(KMM(K), J) + AF(K) * U(KPM(K), J))
              END DO
          END DO

          CALL TRDIAG1M(GGII, GGII, 1_8, N3M, 1_8, N2M, IV)

          DO J = 1, N2M
              DO K = 1, N3M, 2
                  U(K, J) = WW * GGII(K, J) + WW2 * U(K, J)
              END DO
          END DO

!  ANOTHER HALF
!  ------------

          IF (IV == 1) THEN
              TT = U(1, N2M)
              DO J = 1, N2M
                  DO K = 1, N3M, 2
                      U(K, J) = U(K, J) - TT
                  END DO
              END DO
          END IF
 
          DO J = 1, N2M
              DO K = 2, N3M, 2
                  GGII(K, J) = RHS(K, J) - (AB(K) * U(KMM(K), J) + AF(K) * U(KPM(K), J))
              END DO
          END DO

          CALL TRDIAG1M(GGII, GGII, 2_8, N3M, 1_8, N2M, IV)

          DO J = 1, N2M
              DO K = 2, N3M, 2
                  U(K, J) = WW * GGII(K, J) + WW2 * U(K, J)
              END DO
          END DO

          CALL RESID3(U, RHS, IV, ERRMAX)
          IF (ERRMAX < TEST) GOTO 1000

      END DO

      PRINT *, 'ITERATION LIMIT EXCEEDED.'
      WRITE(77, 201) 'SOR', IV, II, ERRMAX * DTCONST * FLOAT(N1MH) / FLOAT(N1M)
 1000 CONTINUE

201   FORMAT(A5, '  IV=', I5, '  II=', I5, '  RG=', ES23.15)

      RETURN
      END SUBROUTINE GSOR2D

!=======================================================================
      SUBROUTINE RESID3(U, RHS, IV, ERRMAX)
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8)   :: K, J, IV
       REAL(8)      :: ERRMAX
       COMPLEX(8)   :: U(0:N3,0:N2), RHS(N3,N2)
       COMPLEX(8)   :: ERR

       ERRMAX = 0.0_8

      DO J = 1, N2M
          DO K = 1, N3M
              ERR = RHS(K, J) - AB(K) * U(KMM(K), J) - AF(K) * U(KPM(K), J)   &
                              - AS(J) * U(K, J - 1) - AN(J) * U(K, J + 1)     &
                              - AC(K, J, IV) * U(K, J)
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
       INTEGER(8)  :: I, J, K, ILEV
       INTEGER(8)  :: MINROW, KBG, KEND, KSP, KC, KBGH
       REAL(8)     :: KBZ(N3MD), KFZ(N3MD)
       REAL(8)     :: ZMPM(0:N3MD)
       REAL(8)     :: VDZ_KBG, VDZ_KEND, SDZ_KBG, SDZ_KEND


       MINROW = N3M / (2**NLEV)
       ZMPM = 0.0_8

       LEVHALF = NINT(NLEV / 2.0_8)
       KKMG(NLEV, 1) = 1    ! START INDEX
       KKMG(NLEV, 2) = N3M  ! END INDEXT
       KKMG(NLEV, 3) = N3M  ! THE NUMBER OF POINTS AT THE ILEV

       DO K = 1, N3M
           ZMPM(K) = ZMP(K)
           KBZ(K) = 1 - 1/K                  ! 0 only if K=1
           KFZ(K) = 1 - K/N3M                ! 0 only if K=N3M
       END DO


!  COMPUTE FOR LOWER LEVELS
      DO ILEV = NLEV - 1, 0, -1

          KKMG(ILEV, 1) = KKMG(ILEV + 1, 1) + KKMG(ILEV + 1, 3)
          KKMG(ILEV, 3) = MINROW * (2**ILEV)
          KKMG(ILEV, 2) = KKMG(ILEV, 1) + KKMG(ILEV, 3) - 1

          KBG = KKMG(ILEV, 1)
          KEND = KKMG(ILEV, 2)

          KSP = 2**(NLEV - ILEV)           ! width of one cell at low level

          KC = 0
          DO K = KBG, KEND
              KC = KC + 1
              ZMPM(K) = 0.5_8 * (Z(KC * KSP + 1) + Z((KC - 1) * KSP + 1))
              KBZ(K) = 1 - KBG/K                  ! 0 onlyif K=KBG
              KFZ(K) = 1 - K/KEND                 ! 0 only if K=KEND
          END DO

          KC = 0
          DO K = KBG, KEND
              KC = KC + 1
              AB(K) = KBZ(K) / ((ZMPM(K) - ZMPM(K - 1)) * (Z(KC * KSP + 1) - Z((KC - 1) * KSP + 1)))
              AF(K) = KFZ(K) / ((ZMPM(K + 1) - ZMPM(K)) * (Z(KC * KSP + 1) - Z((KC - 1) * KSP + 1)))
          END DO

          DO J = 1, N2M
              DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
                  AC(K, J, 1) = -1.0_8 * (AB(K) + AF(K) + AS(J) + AN(J))
              END DO
          END DO

          DO I = 2, N1MH
              DO J = 1, N2M
                  DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
                      AC(K, J, I) = AC(K, J, 1) - AI3(I)
                  END DO
              END DO
          END DO

      END DO

!  CALCULATE RESTRICTION COEFFS
       DO ILEV = NLEV, 1, -1
           KBGH = KKMG(ILEV, 1)
           DO K = KKMG(ILEV - 1, 1), KKMG(ILEV - 1, 2)
               COR1(K) = (ZMPM(KBGH + 1) - ZMPM(K)) / (ZMPM(KBGH + 1) - ZMPM(KBGH))
               COR2(K) = 1.0_8 - COR1(K)
               KBGH = KBGH + 2
           END DO
       END DO

!  CALCULATE INTERPOLATION COEFFS
       DO ILEV = 0, NLEV - 1
           KBGH = KKMG(ILEV + 1, 1) + 1
           DO K = KKMG(ILEV, 1), KKMG(ILEV, 2) - 1
               COI1(K) = (ZMPM(K + 1) - ZMPM(KBGH)) / (ZMPM(K + 1) - ZMPM(K))  ! * lower value
               COI2(K) = 1.0_8 - COI1(K)
               KBGH = KBGH + 2
           END DO
           K = KKMG(ILEV, 2)
           COI1(K) = 1.0_8                 ! use only one lower point at upper wall
       END DO


!===== FOR THE Z PERIODICIRY
!       INTRODUCE KPM & KMM
       IF (ZPRDIC == 1) THEN
           DO ILEV = NLEV, 0, -1
               KBG = KKMG(ILEV, 1)
               KEND = KKMG(ILEV, 2)
               DO K = KBG, KEND
                   KPM(K) = K + 1
                   KMM(K) = K - 1
               END DO
               KPM(KEND) = KBG
               KMM(KBG) = KEND
           END DO

           DO ILEV = NLEV - 1, 0, -1
               KBG = KKMG(ILEV, 1)
               KEND = KKMG(ILEV, 2)
               KSP = 2**(NLEV - ILEV)
               VDZ_KBG  = ZMPM(KBG) - Z(1) + Z(N3) - ZMPM(KEND)
               VDZ_KEND = ZMPM(KBG) - Z(1) + Z(N3) - ZMPM(KEND)
               SDZ_KBG  = Z(1 + KSP) - Z(1)
               SDZ_KEND = Z(N3) - Z(N3 - KSP)
               AB(KBG)  = 1.0_8 / (VDZ_KBG * SDZ_KBG)
               AF(KEND) = 1.0_8 / (VDZ_KEND * SDZ_KEND)
           END DO

           DO ILEV = NLEV - 1, 0, -1
               DO J = 1, N2M
                   DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
                       AC(K, J, 1) = -1.0_8 * (AB(K) + AF(K) + AS(J) + AN(J))
                   END DO
               END DO
   
               DO I = 2, N1MH
                   DO J = 1, N2M
                       DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
                           AC(K, J, I) = AC(K, J, 1) - AI3(I)
                       END DO
                   END DO
               END DO
           END DO

!  CALCULATE INTERPOLATION COEFFS
           DO ILEV = 0, NLEV - 1
               KBG = KKMG(ILEV, 1)
               KEND = KKMG(ILEV, 2)
               VDZ_KEND = ZMPM(KBG) - Z(1) + Z(N3) - ZMPM(KEND)
               COI2(KEND) = (ZMPM(KKMG(ILEV + 1, 2)) - ZMPM(KEND)) / VDZ_KEND
               COI1(KEND) = 1.0_8 - COI2(KEND)
           END DO

       END IF


       DO ILEV = NLEV, 0, -1
           WRITE(*,*) 'IIMG(1', ILEV, ')=', KKMG(ILEV, 1)
           WRITE(*,*) 'IIMG(2', ILEV, ')=', KKMG(ILEV, 2)
           WRITE(*,*) 'IIMG(3', ILEV, ')=', KKMG(ILEV, 3)
       END DO

      RETURN
      END SUBROUTINE COEFMG