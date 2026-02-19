!=======================================================================
!
!     MAIN SOLVER OF POISSON EQUATION
!
!     AX=b
!     A: coefficient of discretized poisson equation
!     X: PHI: Pseudo pressure, Output of subroutine POISSON
!     b: DIVGSUM: OUTPUT of subroutine DIVGS.
!
!     Main algorithm: Preconditioned BICGSTAB 
!     [Biconjugate gradient stabilized method]
!     algorithm: https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
!
!     Precondition: Gauss-Seidal with Multigrid acceleration
!                   x & z direction: Gauss-Seidal + Multigrid
!                   y direction: TDMA
!                   only one V-cycle for the precondition
!
!
!     Sept. 1998, S. Kang: MG3D
!     Jun. 2017,  J. Park: BICG and f90
!     Feb. 2026, S. Lee: Periodic Domain Support and Several Bug Fixes
!
!======================================================================
      SUBROUTINE POISSON(PHI, DIVGSUM)
!======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8) :: I, J, K
       INTEGER(8) :: I_CG, IMAX, JMAX, KMAX
       REAL(8)    :: RHO_CG, RHO_OL, ALP_CG, W_CG, BET_CG, RV_TMP, ACC
       REAL(8)    :: PHIREF, RESM_CG, BTTT, UPPP
       
       REAL(8), DIMENSION(0:N1, 0:N2, 0:N3) :: PHI
       REAL(8), DIMENSION(0:N1, 0:N2, 0:N3) :: DIVGSUM
       REAL(8), DIMENSION(0:N1, 0:N2, 0:N3) :: P_CG
       REAL(8), DIMENSION(0:N1, 0:N2, 0:N3) :: S_CG
       REAL(8), DIMENSION(0:N1, 0:N2, 0:N3) :: Y_CG
       REAL(8), DIMENSION(0:N1, 0:N2, 0:N3) :: Z_CG
       REAL(8), DIMENSION(0:N1, 0:N2, 0:N3) :: V_CG
       REAL(8), DIMENSION(N1, N2, N3)       :: T_CG
       REAL(8), DIMENSION(N1, N2, N3)       :: RES_CG

       P_CG = 0.0_8
       S_CG = 0.0_8
       Y_CG = 0.0_8
       Z_CG = 0.0_8
       V_CG = 0.0_8
       T_CG = 0.0_8

      IF (IOLDV == 0) THEN
          PHI = 0.0_8
!$OMP PARALLEL DO private(I, J)
          DO K = 1, N3M
              DO J = 1, N2M
                  DO I = 1, N1M
                      RES_CG(I, J, K) = DIVGSUM(I, J, K)
                  END DO
              END DO
          END DO
      ELSE
!$OMP PARALLEL DO private(I, J, ACC)
          DO K = 1, N3M
              DO J = 1, N2M
                  DO I = 1, N1M
                      ACC = -1.0_8 * (AW(I) + AE(I) + AS(J) + AN(J) + AB(K) + AF(K))
                      RES_CG(I, J, K) = DIVGSUM(I, J, K) + (-ACC * PHI(I, J, K) &
                                    - AW(I) * PHI(IMM(I), J, K) - AE(I) * PHI(IPM(I), J, K) &
                                    - AS(J) * PHI(I, J - 1, K) - AN(J) * PHI(I, J + 1, K) &
                                    - AB(K) * PHI(I, J, KMM(K)) - AF(K) * PHI(I, J, KPM(K)))
                  END DO
              END DO
          END DO
      END IF

       RHO_OL = 1.0_8
       ALP_CG = 1.0_8
       W_CG   = 1.0_8

      DO I_CG = 1, MGITR

          RHO_CG = 0.0_8
!$OMP PARALLEL DO private(I, J) reduction(+:RHO_CG)
          DO K = 1, N3M
              DO J = 1, N2M
                  DO I = 1, N1M
                      RHO_CG = RHO_CG + DIVGSUM(I, J, K) * RES_CG(I, J, K)
                  END DO
              END DO
          END DO

          BET_CG = RHO_CG / RHO_OL * ALP_CG / W_CG

!$OMP PARALLEL DO private(I, J)
          DO K = 1, N3M
              DO J = 1, N2M
                  DO I = 1, N1M
                      P_CG(I, J, K) = RES_CG(I, J, K) + BET_CG * (P_CG(I, J, K) - W_CG * V_CG(I, J, K))
                  END DO
              END DO
          END DO

          CALL MG3D(Y_CG, P_CG, TEST1, IOLDV)

!$OMP PARALLEL DO private(I, J, ACC)
          DO K = 1, N3M
              DO J = 1, N2M
                  DO I = 1, N1M
                      ACC = -1.0_8 * (AW(I) + AE(I) + AS(J) + AN(J) + AB(K) + AF(K))
                      V_CG(I, J, K) = AW(I) * Y_CG(IMM(I), J, K) + AE(I) * Y_CG(IPM(I), J, K) &
                                    + AS(J) * Y_CG(I, J - 1, K) + AN(J) * Y_CG(I, J + 1, K) &
                                    + AB(K) * Y_CG(I, J, KMM(K)) + AF(K) * Y_CG(I, J, KPM(K)) &
                                    + ACC * Y_CG(I, J, K)
                  END DO
              END DO
          END DO

          RV_TMP = 0.0_8
!$OMP PARALLEL DO private(I, J) reduction(+:RV_TMP)
          DO K = 1, N3M
              DO J = 1, N2M
                  DO I = 1, N1M
                      RV_TMP = RV_TMP + DIVGSUM(I, J, K) * V_CG(I, J, K)
                  END DO
              END DO
          END DO

          ALP_CG = RHO_CG / RV_TMP

!$OMP PARALLEL DO private(I, J)
          DO K = 1, N3M
              DO J = 1, N2M
                  DO I = 1, N1M
                      S_CG(I, J, K) = RES_CG(I, J, K) - ALP_CG * V_CG(I, J, K)
                  END DO
              END DO
          END DO

          CALL MG3D(Z_CG, S_CG, TEST1, IOLDV)

!$OMP PARALLEL DO private(I, J, ACC)
          DO K = 1, N3M
              DO J = 1, N2M
                  DO I = 1, N1M
                      ACC = -1.0_8 * (AW(I) + AE(I) + AS(J) + AN(J) + AB(K) + AF(K))
                      T_CG(I, J, K) = AW(I) * Z_CG(IMM(I), J, K) + AE(I) * Z_CG(IPM(I), J, K) &
                                    + AS(J) * Z_CG(I, J - 1, K) + AN(J) * Z_CG(I, J + 1, K) &
                                    + AB(K) * Z_CG(I, J, KMM(K)) + AF(K) * Z_CG(I, J, KPM(K)) &
                                    + ACC * Z_CG(I, J, K)
                  END DO
              END DO
          END DO

          BTTT = 0.0_8
          UPPP = 0.0_8
!$OMP PARALLEL DO private(I, J) reduction(+:UPPP) reduction(+:BTTT)
          DO K = 1, N3M
              DO J = 1, N2M
                  DO I = 1, N1M
                      UPPP = UPPP + T_CG(I, J, K) * S_CG(I, J, K)
                      BTTT = BTTT + T_CG(I, J, K)**2.0_8
                  END DO
              END DO
          END DO
          W_CG = UPPP / BTTT

!$OMP PARALLEL DO private(I, J)
          DO K = 1, N3M
              DO J = 1, N2M
                  DO I = 1, N1M
                      PHI(I, J, K) = PHI(I, J, K) + ALP_CG * Y_CG(I, J, K) + W_CG * Z_CG(I, J, K)
                  END DO
              END DO
          END DO

          RESM_CG = 0.0_8
!$OMP PARALLEL DO private(I, J) reduction(MAX:RESM_CG)
          DO K = 1, N3M
              DO J = 1, N2M
                  DO I = 1, N1M
                      RES_CG(I, J, K) = S_CG(I, J, K) - W_CG * T_CG(I, J, K)
                      RESM_CG = MAX(RESM_CG, ABS(RES_CG(I, J, K)))
                      IF(RESM_CG == ABS(RES_CG(I, J, K))) THEN
                          IMAX = I
                          JMAX = J
                          KMAX = K
                      END IF
                  END DO
              END DO
          END DO
          RHO_OL = RHO_CG

          WRITE(77, 102) I_CG, RESM_CG * DTCONST, IMAX, JMAX, KMAX
 102      FORMAT(I4, ES15.7, 3I4)

          IF (RESM_CG < TEST1) THEN
              WRITE(*, 300) I_CG, RESM_CG * DTCONST
              GOTO 1000
          END IF

      END DO

      PRINT*, '=== MUTLGRID : NOT CONVERGED ==='
      WRITE(*, 300) I_CG, RESM_CG * DTCONST
 300  FORMAT('ICYC=', I10, '  RESMAX=', ES25.12)

 1000 CONTINUE

      PHIREF = PHI(1, N2M, 1)
!$OMP PARALLEL DO private(I, J)
      DO K = 1, N3M
          DO J = 1, N2M
              DO I = 1, N1M
                  PHI(I, J, K) = PHI(I, J, K) - PHIREF
              END DO
          END DO
      END DO

      WRITE(78, 103) TIME, DT, I_CG, RESM_CG * DTCONST
 103  FORMAT(2F13.5, I5, ES15.7)

      RETURN
      END SUBROUTINE POISSON

!======================================================================
        SUBROUTINE MG3D(SOL, RHS, EPSIL, IOLD)
!======================================================================
!       contain overall multigrid procedures
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8)  :: I, J, K
       REAL(8)     :: SOL(0:N1, 0:N2, 0:N3), RHS(0:N1, 0:N2, 0:N3), EPSIL
       REAL(8)     :: RESID(N1MD, N2M, N3MD), GI(0:N1MD, 0:N2, 0:N3MD)
       INTEGER(8)  :: IOLD, ICYC, ILEV
       REAL(8)     :: RESM

        CALL TOPLEVEL(GI, RESID, RHS, RESM, IOLD)

        IF (NLEV == 0) GOTO 777
        DO ILEV = NLEV - 1, 1, -1
           CALL RELAX(GI, RESID, ILEV, 1_8, 0_8)
           CALL RESCAL(GI, RESID, ILEV)
           CALL GODOWN(RESID, ILEV - 1)
        END DO

        CALL RELAX(GI, RESID, 0_8, NBLI, 0_8)

        DO ILEV = 1, NLEV - 1
           CALL GOUP(GI, ILEV)
           CALL RELAX(GI, RESID, ILEV, 1_8, 1_8)
        END DO

        CALL GOUP(GI, NLEV)
 777    CALL TOPLEVEL(GI, RESID, RHS, RESM, 1_8)

!$OMP PARALLEL DO private(I, J)
        DO K = 1, N3M
            DO J = 1, N2M
                DO I = 1, N1M
                    SOL(I, J, K) = GI(I, J, K)
                END DO
            END DO
        END DO

        RETURN
        END SUBROUTINE MG3D

!======================================================================
        SUBROUTINE TOPLEVEL(GI, RESID, RHS, RESM, IOLD)   ! Zebra Version
!======================================================================
!       solve the Poisson equation with Zebra GS at the top level
!       with only 1 iteration
!       IOLD = 0 (GODOWN), 1 (GOUP)
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8) :: I, J, K
       INTEGER(8) :: ISTART, IOLD
       REAL(8)    :: RHS(0:N1, 0:N2, 0:N3)
       REAL(8)    :: RESID(N1MD, N2M, N3MD), GI(0:N1MD, 0:N2, 0:N3MD)
       REAL(8)    :: ACC, RESM

!------ make RHS of Zebra GS (odd line)
!$OMP PARALLEL DO private(ISTART, I, J)
        DO K = 1, N3M
            ISTART = 1 + MOD(KPM(K), 2_8)
            DO J = 1, N2M
                DO I = ISTART, N1M, 2
                    GI(I, J, K) = RHS(I, J, K) - IOLD * &
                                  (AW(I) * GI(IMM(I), J, K) + AE(I) * GI(IPM(I), J, K) &
                                 + AB(K) * GI(I, J, KMM(K)) + AF(K) * GI(I, J, KPM(K)))
                END DO
            END DO
        END DO

!------ solve TDMA
!$OMP PARALLEL DO private(ISTART)
        DO K = 1, N3M
            ISTART = 1 + MOD(KPM(K), 2_8)
            DO I = ISTART, N1M, 2
                GI(I, 1, K) = GI(I, 1, K) * BET(I, 1, K)
            END DO
        END DO

!$OMP PARALLEL DO private(ISTART)
        DO K = 1, N3M
            ISTART = 1 + MOD(KPM(K), 2_8)
            DO J = 2, N2M
                DO I = ISTART, N1M, 2
                    GI(I, J, K) = (GI(I, J, K) - AS(J) * GI(I, J - 1, K)) * BET(I, J, K)
                END DO
            END DO
        END DO

!$OMP PARALLEL DO private(ISTART)
        DO K = 1, N3M
            ISTART = 1 + MOD(KPM(K), 2_8)
            DO J = N2M - 1, 1, -1
                DO I = ISTART, N1M, 2
                    GI(I, J, K) = GI(I, J, K) - GAM(I, J + 1, K) * GI(I, J + 1, K)
                END DO
            END DO
        END DO

!------ make RHS of Zebra GS (even line)
!$OMP PARALLEL DO private(ISTART, I, J)
        DO K = 1, N3M
            ISTART = 1 + MOD(K, 2_8)
            DO J = 1, N2M
                DO I = ISTART, N1M, 2
                    GI(I, J, K) = RHS(I, J, K) &
                                 - (AW(I) * GI(IMM(I), J, K) + AE(I) * GI(IPM(I), J, K) &
                                  + AB(K) * GI(I, J, KMM(K)) + AF(K) * GI(I, J, KPM(K)))
                END DO
            END DO
        END DO

!------ solve TDMA
!$OMP PARALLEL DO private(ISTART)
        DO K = 1, N3M
            ISTART = 1 + MOD(K, 2_8)
            DO I = ISTART, N1M, 2
                GI(I, 1, K) = GI(I, 1, K) * BET(I, 1, K)
            END DO
        END DO

!$OMP PARALLEL DO private(ISTART)
        DO K = 1, N3M
            ISTART = 1 + MOD(K, 2_8)
            DO J = 2, N2M
                DO I = ISTART, N1M, 2
                    GI(I, J, K) = (GI(I, J, K) - AS(J) * GI(I, J - 1, K)) * BET(I, J, K)
                END DO
            END DO
        END DO

!$OMP PARALLEL DO private(ISTART)
        DO K = 1, N3M
            ISTART = 1 + MOD(K, 2_8)
            DO J = N2M - 1, 1, -1
                DO I = ISTART, N1M, 2
                    GI(I, J, K) = GI(I, J, K) - GAM(I, J + 1, K) * GI(I, J + 1, K)
                END DO
            END DO
        END DO

!------ calculate residual

        RESM = 0.0_8
!$OMP PARALLEL DO private(I, J, ACC) reduction(MAX:RESM)
        DO K = 1, N3M
            DO J = 1, N2M
                DO I = 1, N1M
                    ACC = -1.0_8 * (AW(I) + AE(I) + AS(J) + AN(J) + AB(K) + AF(K))
                    RESID(I, J, K) = RHS(I, J, K) - ACC * GI(I, J, K) &
                                    - AW(I) * GI(IMM(I), J, K) - AE(I) * GI(IPM(I), J, K) &
                                    - AS(J) * GI(I, J - 1, K) - AN(J) * GI(I, J + 1, K) &
                                    - AB(K) * GI(I, J, KMM(K)) - AF(K) * GI(I, J, KPM(K))
                    RESM = MAX(RESM, ABS(RESID(I, J, K)))
                END DO
            END DO
        END DO

        IF (NLEV == 0) GOTO 97

 101    FORMAT(ES15.8, 3I6)

!------ restrict residual to a lower level
!$OMP PARALLEL DO private(I, J)
        DO K = KKMG(NLEV - 1, 1), KKMG(NLEV - 1, 2)
            DO J = 1, N2M
                DO I = IIMG(NLEV - 1, 1), IIMG(NLEV - 1, 2)
                    RESID(I, J, K) = (1.0_8 - FIDW(I)) * (1.0_8 - FKDW(K)) * RESID(IH1(I), J, KH1(K)) &
                                   + FIDW(I) * (1.0_8 - FKDW(K)) * RESID(IH2(I), J, KH1(K)) &
                                   + (1.0_8 - FIDW(I)) * FKDW(K) * RESID(IH1(I), J, KH2(K)) &
                                   + FIDW(I) * FKDW(K) * RESID(IH2(I), J, KH2(K))
                END DO
            END DO
        END DO

 97     RETURN
        END SUBROUTINE TOPLEVEL

!======================================================================
        SUBROUTINE RESCAL(GI, RESID, ILEV)
!======================================================================
!      calculate residual at each level
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8)  :: I, J, K
       INTEGER(8)  :: ILEV
       REAL(8)     :: RESID(N1MD, N2M, N3MD), GI(0:N1MD, 0:N2, 0:N3MD)
       REAL(8)     :: ACC

!$OMP PARALLEL DO private(I, J, ACC)
        DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
            DO J = 1, N2M
                DO I = IIMG(ILEV, 1), IIMG(ILEV, 2)
                    ACC = -1.0_8 * (AW(I) + AE(I) + AS(J) + AN(J) + AB(K) + AF(K))
                    RESID(I, J, K) = RESID(I, J, K) - ACC * GI(I, J, K) &
                                    - AW(I) * GI(IMM(I), J, K) - AE(I) * GI(IPM(I), J, K) &
                                    - AS(J) * GI(I, J - 1, K) - AN(J) * GI(I, J + 1, K) &
                                    - AB(K) * GI(I, J, KMM(K)) - AF(K) * GI(I, J, KPM(K))
                END DO
            END DO
        END DO

        RETURN
        END SUBROUTINE RESCAL

!======================================================================
        SUBROUTINE GODOWN(RESID, ILEV)
!======================================================================
!      restrict residual to a lower level
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8)  :: I, J, K
       INTEGER(8)  :: ILEV
       REAL(8)     :: RESID(N1MD, N2M, N3MD)

!$OMP PARALLEL DO private(I, J)
        DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
            DO J = 1, N2M
                DO I = IIMG(ILEV, 1), IIMG(ILEV, 2)
                    RESID(I, J, K) = (1.0_8 - FIDW(I)) * (1.0_8 - FKDW(K)) * RESID(IH1(I), J, KH1(K)) &
                                   + FIDW(I) * (1.0_8 - FKDW(K)) * RESID(IH2(I), J, KH1(K)) &
                                   + (1.0_8 - FIDW(I)) * FKDW(K) * RESID(IH1(I), J, KH2(K)) &
                                   + FIDW(I) * FKDW(K) * RESID(IH2(I), J, KH2(K))
                END DO
            END DO
        END DO

        RETURN
        END SUBROUTINE GODOWN

!======================================================================
        SUBROUTINE GOUP(GI, ILEV)
!======================================================================
!       interpolate residual & add it to a higher level
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8)  :: I, J, K
       INTEGER(8)  :: ILEV, ISTART
       REAL(8)     :: GI(0:N1MD, 0:N2, 0:N3MD)

!$OMP PARALLEL DO private(ISTART, I, J)
        DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
            ISTART = IIMG(ILEV, 1) + MOD(K, 2_8)
            DO J = 1, N2M
                DO I = ISTART, IIMG(ILEV, 2), 2
                    GI(I, J, K) = GI(I, J, K) &
                                + (1.0_8 - FIUP(I)) * (1.0_8 - FKUP(K)) * GI(IL1(I), J, KL1(K)) &
                                + FIUP(I) * (1.0_8 - FKUP(K)) * GI(IL2(I), J, KL1(K)) &
                                + (1.0_8 - FIUP(I)) * FKUP(K) * GI(IL1(I), J, KL2(K)) &
                                + FIUP(I) * FKUP(K) * GI(IL2(I), J, KL2(K))
                END DO
            END DO
        END DO

        RETURN
        END SUBROUTINE GOUP

!======================================================================
        SUBROUTINE RELAX(GI, RESID, ILEV, ITER, IOLD)   ! Zebra Version
!======================================================================
!       solve the Poisson equation with Zebra GS
!       IOLD = 0 (GODOWN), 1 (GOUP)
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8)  :: I, J, K
       INTEGER(8)  :: ILEV, ITER, IOLD, ISTART, II
       REAL(8)     :: RESID(N1MD, N2M, N3MD), GI(0:N1MD, 0:N2, 0:N3MD)

!====== 1st iteration
!------ make RHS of Zebra GS (odd line)
!$OMP PARALLEL DO private(ISTART, I, J)
        DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
            ISTART = IIMG(ILEV, 1) + MOD(KPM(K), 2_8)
            DO J = 1, N2M
                DO I = ISTART, IIMG(ILEV, 2), 2
                    GI(I, J, K) = RESID(I, J, K) - IOLD * &
                                  (AW(I) * GI(IMM(I), J, K) + AE(I) * GI(IPM(I), J, K) &
                                 + AB(K) * GI(I, J, KMM(K)) + AF(K) * GI(I, J, KPM(K)))
                END DO
            END DO
        END DO

!------ solve TDMA
!$OMP PARALLEL DO private(ISTART)
        DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
            ISTART = IIMG(ILEV, 1) + MOD(KPM(K), 2_8)
            DO I = ISTART, IIMG(ILEV, 2), 2
                GI(I, 1, K) = GI(I, 1, K) * BET(I, 1, K)
            END DO
        END DO

!$OMP PARALLEL DO private(ISTART)
        DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
            ISTART = IIMG(ILEV, 1) + MOD(KPM(K), 2_8)
            DO J = 2, N2M
                DO I = ISTART, IIMG(ILEV, 2), 2
                    GI(I, J, K) = (GI(I, J, K) - AS(J) * GI(I, J - 1, K)) * BET(I, J, K)
                END DO
            END DO
        END DO

!$OMP PARALLEL DO private(ISTART)
        DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
            ISTART = IIMG(ILEV, 1) + MOD(KPM(K), 2_8)
            DO J = N2M - 1, 1, -1
                DO I = ISTART, IIMG(ILEV, 2), 2
                    GI(I, J, K) = GI(I, J, K) - GAM(I, J + 1, K) * GI(I, J + 1, K)
                END DO
            END DO
        END DO

!------ make RHS of Zebra GS (even line)
!$OMP PARALLEL DO private(ISTART, I, J)
        DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
            ISTART = IIMG(ILEV, 1) + MOD(K, 2_8)
            DO J = 1, N2M
                DO I = ISTART, IIMG(ILEV, 2), 2
                    GI(I, J, K) = RESID(I, J, K) &
                                 - (AW(I) * GI(IMM(I), J, K) + AE(I) * GI(IPM(I), J, K) &
                                  + AB(K) * GI(I, J, KMM(K)) + AF(K) * GI(I, J, KPM(K)))
                END DO
            END DO
        END DO

!------ solve TDMA
!$OMP PARALLEL DO private(ISTART)
        DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
            ISTART = IIMG(ILEV, 1) + MOD(K, 2_8)
            DO I = ISTART, IIMG(ILEV, 2), 2
                GI(I, 1, K) = GI(I, 1, K) * BET(I, 1, K)
            END DO
        END DO

!$OMP PARALLEL DO private(ISTART)
        DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
            ISTART = IIMG(ILEV, 1) + MOD(K, 2_8)
            DO J = 2, N2M
                DO I = ISTART, IIMG(ILEV, 2), 2
                    GI(I, J, K) = (GI(I, J, K) - AS(J) * GI(I, J - 1, K)) * BET(I, J, K)
                END DO
            END DO
        END DO

!$OMP PARALLEL DO private(ISTART)
        DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
            ISTART = IIMG(ILEV, 1) + MOD(K, 2_8)
            DO J = N2M - 1, 1, -1
                DO I = ISTART, IIMG(ILEV, 2), 2
                    GI(I, J, K) = GI(I, J, K) - GAM(I, J + 1, K) * GI(I, J + 1, K)
                END DO
            END DO
        END DO

!====== repeat previous procedures
        DO II = 1, ITER - 1

!------ make RHS of Zebra GS (odd line)
!$OMP PARALLEL DO private(ISTART, I, J)
            DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
                ISTART = IIMG(ILEV, 1) + MOD(KPM(K), 2_8)
                DO J = 1, N2M
                    DO I = ISTART, IIMG(ILEV, 2), 2
                        GI(I, J, K) = RESID(I, J, K) &
                                     - (AW(I) * GI(IMM(I), J, K) + AE(I) * GI(IPM(I), J, K) &
                                      + AB(K) * GI(I, J, KMM(K)) + AF(K) * GI(I, J, KPM(K)))
                    END DO
                END DO
            END DO

!------ solve TDMA
!$OMP PARALLEL DO private(ISTART)
            DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
                ISTART = IIMG(ILEV, 1) + MOD(KPM(K), 2_8)
                DO I = ISTART, IIMG(ILEV, 2), 2
                    GI(I, 1, K) = GI(I, 1, K) * BET(I, 1, K)
                END DO
            END DO

!$OMP PARALLEL DO private(ISTART)
            DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
                ISTART = IIMG(ILEV, 1) + MOD(KPM(K), 2_8)
                DO J = 2, N2M
                    DO I = ISTART, IIMG(ILEV, 2), 2
                        GI(I, J, K) = (GI(I, J, K) - AS(J) * GI(I, J - 1, K)) * BET(I, J, K)
                    END DO
                END DO
            END DO

!$OMP PARALLEL DO private(ISTART)
            DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
                ISTART = IIMG(ILEV, 1) + MOD(KPM(K), 2_8)
                DO J = N2M - 1, 1, -1
                    DO I = ISTART, IIMG(ILEV, 2), 2
                        GI(I, J, K) = GI(I, J, K) - GAM(I, J + 1, K) * GI(I, J + 1, K)
                    END DO
                END DO
            END DO

!------ make RHS of Zebra GS (even line)
!$OMP PARALLEL DO private(ISTART, I, J)
            DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
                ISTART = IIMG(ILEV, 1) + MOD(K, 2_8)
                DO J = 1, N2M
                    DO I = ISTART, IIMG(ILEV, 2), 2
                        GI(I, J, K) = RESID(I, J, K) &
                                     - (AW(I) * GI(IMM(I), J, K) + AE(I) * GI(IPM(I), J, K) &
                                      + AB(K) * GI(I, J, KMM(K)) + AF(K) * GI(I, J, KPM(K)))
                    END DO
                END DO
            END DO

!------ solve TDMA
!$OMP PARALLEL DO private(ISTART)
            DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
                ISTART = IIMG(ILEV, 1) + MOD(K, 2_8)
                DO I = ISTART, IIMG(ILEV, 2), 2
                    GI(I, 1, K) = GI(I, 1, K) * BET(I, 1, K)
                END DO
            END DO

!$OMP PARALLEL DO private(ISTART)
            DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
                ISTART = IIMG(ILEV, 1) + MOD(K, 2_8)
                DO J = 2, N2M
                    DO I = ISTART, IIMG(ILEV, 2), 2
                        GI(I, J, K) = (GI(I, J, K) - AS(J) * GI(I, J - 1, K)) * BET(I, J, K)
                    END DO
                END DO
            END DO

!$OMP PARALLEL DO private(ISTART)
            DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
                ISTART = IIMG(ILEV, 1) + MOD(K, 2_8)
                DO J = N2M - 1, 1, -1
                    DO I = ISTART, IIMG(ILEV, 2), 2
                        GI(I, J, K) = GI(I, J, K) - GAM(I, J + 1, K) * GI(I, J + 1, K)
                    END DO
                END DO
            END DO

        END DO

        RETURN
        END SUBROUTINE RELAX

!======================================================================
      SUBROUTINE POISINIT  ! CALCULATE COEFS FOR POISSON EQ.
!======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8)  :: I, J, K
       INTEGER(8)  :: IP, JP, KP

       CALL MGRD_ALLO

       ALLOCATE(GAM(N1MD, N2M, N3MD))
       ALLOCATE(BET(N1MD, N2M, N3MD))

!-----COMPUTE COEFS OF POISSON EQ.
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

      DO K = 1, N3M
          KP = KPV(K)
          AB(K) = (1.0_8 - FIXKL(K)) * C2CZI(K) * F2FZI(K)
          AF(K) = (1.0_8 - FIXKU(K)) * C2CZI(KP) * F2FZI(K)
      END DO

      OPEN(77, FILE='../output/ftr/poiss_itr.dat')
      OPEN(78, FILE='../output/ftr/ftrpoittr.dat')

      CALL MGCOEF

      RETURN
      END SUBROUTINE POISINIT

!======================================================================
      SUBROUTINE MGCOEF  ! CALCULATE COEFS FOR POISSON EQ.
!======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8) :: I, J, K
       REAL(8)    :: XMPM(0:N1MD), ZMPM(0:N3MD)
       INTEGER(8) :: IWEST(N1MD), IEAST(N1MD), KBACK(N3MD), KFORW(N3MD)
       INTEGER(8) :: ILEV
       INTEGER(8) :: IC, IBG, ISP, IEND, KC, KBG, KSP, KEND, IBGH, KBGH, IBGL, KBGL, KENDL, IENDL
       REAL(8)    :: VDZ_KBG, VDZ_KEND, SDZ_KBG, SDZ_KEND, DZ, DZ1, DZ2
       REAL(8)    :: VDX_IBG, VDX_IEND, SDX_IBG, SDX_IEND, DX, DX1, DX2

!====== calculate indices at each level
!------ IIMG(),KKMG(),IKMG()
        IIMG(NLEV, 1) = 1
        IIMG(NLEV, 2) = N1M
        IIMG(NLEV, 3) = N1M
        KKMG(NLEV, 1) = 1
        KKMG(NLEV, 2) = N3M
        KKMG(NLEV, 3) = N3M

        DO ILEV = NLEV - 1, 0, -1
            IIMG(ILEV, 1) = IIMG(ILEV + 1, 1) + IIMG(ILEV + 1, 3)
            IIMG(ILEV, 3) = (N1M / (2**NLEV)) * (2**ILEV)
            IIMG(ILEV, 2) = IIMG(ILEV, 1) + IIMG(ILEV, 3) - 1
            
            KKMG(ILEV, 1) = KKMG(ILEV + 1, 1) + KKMG(ILEV + 1, 3)
            KKMG(ILEV, 3) = (N3M / (2**NLEV)) * (2**ILEV)
            KKMG(ILEV, 2) = KKMG(ILEV, 1) + KKMG(ILEV, 3) - 1
        END DO

!====== compute for the finest grid
        DO I = 1, N1M
            XMPM(I) = XMP(I)
            IWEST(I) = 1 - 1/I          ! 1 for I > 1
            IEAST(I) = 1 - I/N1M        ! 1 for I < N1M
        END DO
        
        DO K = 1, N3M
            ZMPM(K) = ZMP(K)
            KBACK(K) = 1 - 1/K          ! 1 for K > 1
            KFORW(K) = 1 - K/N3M        ! 1 for K < N3M
        END DO

!====== compute for coarse grids
        DO ILEV = NLEV - 1, 0, -1

            IBG = IIMG(ILEV, 1)
            IEND = IIMG(ILEV, 2)
            ISP = 2**(NLEV - ILEV)
            
            KBG = KKMG(ILEV, 1)
            KEND = KKMG(ILEV, 2)
            KSP = 2**(NLEV - ILEV)

            IC = 0
            DO I = IBG, IEND
                IC = IC + 1
                XMPM(I) = 0.5_8 * (X((IC - 1) * ISP + 1) + X(IC * ISP + 1))
                IWEST(I) = 1 - IBG/I        ! 1 for I > IBG
                IEAST(I) = 1 - I/IEND       ! 1 for I < IEND
            END DO

            KC = 0
            DO K = KBG, KEND
                KC = KC + 1
                ZMPM(K) = 0.5_8 * (Z((KC - 1) * KSP + 1) + Z(KC * KSP + 1))
                KBACK(K) = 1 - KBG/K        ! 1 for K > KBG
                KFORW(K) = 1 - K/KEND       ! 1 for K < KEND
            END DO

!------ calculate Poisson coefficients for coarse grids
            IC = 0
            DO I = IBG, IEND
                IC = IC + 1
                AW(I) = IWEST(I) * 1.0_8 / ((XMPM(I) - XMPM(I - 1)) * (X(IC * ISP + 1) - X((IC - 1) * ISP + 1)))
                AE(I) = IEAST(I) * 1.0_8 / ((XMPM(I + 1) - XMPM(I)) * (X(IC * ISP + 1) - X((IC - 1) * ISP + 1)))
            END DO

            KC = 0
            DO K = KBG, KEND
                KC = KC + 1
                AB(K) = KBACK(K) * 1.0_8 / ((ZMPM(K) - ZMPM(K - 1)) * (Z(KC * KSP + 1) - Z((KC - 1) * KSP + 1)))
                AF(K) = KFORW(K) * 1.0_8 / ((ZMPM(K + 1) - ZMPM(K)) * (Z(KC * KSP + 1) - Z((KC - 1) * KSP + 1)))
            END DO

        END DO

!====== calculate restriction coefficients
        DO ILEV = 0, NLEV - 1

            IBGH = IIMG(ILEV + 1, 1)             ! at higher level
            DO I = IIMG(ILEV, 1), IIMG(ILEV, 2)
                FIDW(I) = (XMPM(I) - XMPM(IBGH)) / (XMPM(IBGH + 1) - XMPM(IBGH))
                IH1(I) = IBGH
                IH2(I) = IBGH + 1
                IBGH = IBGH + 2
            END DO

            KBGH = KKMG(ILEV + 1, 1)             ! at higher level
            DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
                FKDW(K) = (ZMPM(K) - ZMPM(KBGH)) / (ZMPM(KBGH + 1) - ZMPM(KBGH))
                KH1(K) = KBGH
                KH2(K) = KBGH + 1
                KBGH = KBGH + 2
            END DO

        END DO

!====== calculate prolongation coefficients
        DO ILEV = 1, NLEV

            IBGL = IIMG(ILEV - 1, 1)             ! at lower level
            DO I = IIMG(ILEV, 1), IIMG(ILEV, 2), 2
                FIUP(I) = (XMPM(I) - XMPM(IBGL - 1)) / (XMPM(IBGL) - XMPM(IBGL - 1))
                FIUP(I + 1) = (XMPM(I + 1) - XMPM(IBGL)) / (XMPM(IBGL + 1) - XMPM(IBGL))
                IL1(I) = IBGL - 1
                IL2(I) = IBGL
                IL1(I + 1) = IBGL
                IL2(I + 1) = IBGL + 1
                IBGL = IBGL + 1
            END DO
            I = IIMG(ILEV, 1)
            FIUP(I) = 1.0_8                      ! Neumann B.C.
            I = IIMG(ILEV, 2)
            FIUP(I) = 0.0_8                      ! Neumann B.C.

            KBGL = KKMG(ILEV - 1, 1)             ! at lower level
            DO K = KKMG(ILEV, 1), KKMG(ILEV, 2), 2
                FKUP(K) = (ZMPM(K) - ZMPM(KBGL - 1)) / (ZMPM(KBGL) - ZMPM(KBGL - 1))
                FKUP(K + 1) = (ZMPM(K + 1) - ZMPM(KBGL)) / (ZMPM(KBGL + 1) - ZMPM(KBGL))
                KL1(K) = KBGL - 1
                KL2(K) = KBGL
                KL1(K + 1) = KBGL
                KL2(K + 1) = KBGL + 1
                KBGL = KBGL + 1
            END DO
            K = KKMG(ILEV, 1)
            FKUP(K) = 1.0_8                      ! Neumann B.C.
            K = KKMG(ILEV, 2)
            FKUP(K) = 0.0_8                      ! Neumann B.C.

        END DO

!       INTRODUCE DEFAULT BOUNDARY INDEX MAPPING
        DO ILEV = NLEV, 0, -1
            KBG = KKMG(ILEV, 1)
            KEND = KKMG(ILEV, 2)
            DO K = KBG, KEND
                KPM(K) = K + 1
                KMM(K) = K - 1
            END DO
            
            IBG = IIMG(ILEV, 1)
            IEND = IIMG(ILEV, 2)
            DO I = IBG, IEND
                IPM(I) = I + 1
                IMM(I) = I - 1
            END DO
        END DO

!       FOR X-PERIODICITY
        IF (XPRDIC == 1) THEN
            DO ILEV = NLEV, 0, -1
                IBG = IIMG(ILEV, 1)
                IEND = IIMG(ILEV, 2)
                IPM(IEND) = IBG
                IMM(IBG) = IEND
            END DO

            DO ILEV = NLEV - 1, 0, -1
                IBG = IIMG(ILEV, 1)
                IEND = IIMG(ILEV, 2)
                ISP = 2**(NLEV - ILEV)
                VDX_IBG = XMPM(IBG) - X(1) + X(N1) - XMPM(IEND)
                VDX_IEND = XMPM(IBG) - X(1) + X(N1) - XMPM(IEND)
                SDX_IBG = X(1 + ISP) - X(1)
                SDX_IEND = X(N1) - X(N1 - ISP)
                AW(IBG) = 1.0_8 / (VDX_IBG * SDX_IBG)
                AE(IEND) = 1.0_8 / (VDX_IEND * SDX_IEND)
            END DO

            DO ILEV = 1, NLEV
                IBG = IIMG(ILEV, 1)
                IEND = IIMG(ILEV, 2)
                IBGL = IIMG(ILEV - 1, 1)
                IENDL = IIMG(ILEV - 1, 2)
                DX = XMPM(IBGL) - X(1) + X(N1) - XMPM(IENDL)
                DX1 = XMPM(IEND) - XMPM(IENDL)
                DX2 = X(N1) - XMPM(IENDL) + XMPM(IBG) - X(1)
                FIUP(IBG) = DX2 / DX
                FIUP(IEND) = DX1 / DX
                IL1(IBG) = IENDL
                IL2(IBG) = IBGL
                IL1(IEND) = IENDL
                IL2(IEND) = IBGL
            END DO
        END IF

!       FOR Z-PERIODICITY
        IF (ZPRDIC == 1) THEN
            DO ILEV = NLEV, 0, -1
                KBG = KKMG(ILEV, 1)
                KEND = KKMG(ILEV, 2)
                KPM(KEND) = KBG
                KMM(KBG) = KEND
            END DO

            DO ILEV = NLEV - 1, 0, -1
                KBG = KKMG(ILEV, 1)
                KEND = KKMG(ILEV, 2)
                KSP = 2**(NLEV - ILEV)
                VDZ_KBG = ZMPM(KBG) - Z(1) + Z(N3) - ZMPM(KEND)
                VDZ_KEND = ZMPM(KBG) - Z(1) + Z(N3) - ZMPM(KEND)
                SDZ_KBG = Z(1 + KSP) - Z(1)
                SDZ_KEND = Z(N3) - Z(N3 - KSP)
                AB(KBG) = 1.0_8 / (VDZ_KBG * SDZ_KBG)
                AF(KEND) = 1.0_8 / (VDZ_KEND * SDZ_KEND)
            END DO

            DO ILEV = 1, NLEV
                KBG = KKMG(ILEV, 1)
                KEND = KKMG(ILEV, 2)
                KBGL = KKMG(ILEV - 1, 1)
                KENDL = KKMG(ILEV - 1, 2)
                DZ = ZMPM(KBGL) - Z(1) + Z(N3) - ZMPM(KENDL)
                DZ1 = ZMPM(KEND) - ZMPM(KENDL)
                DZ2 = Z(N3) - ZMPM(KENDL) + ZMPM(KBG) - Z(1)
                FKUP(KBG) = DZ2 / DZ
                FKUP(KEND) = DZ1 / DZ
                KL1(KBG) = KENDL
                KL2(KBG) = KBGL
                KL1(KEND) = KENDL
                KL2(KEND) = KBGL
            END DO
        END IF

        CALL TRIDCOEF

      RETURN
      END SUBROUTINE MGCOEF


!======================================================================
        SUBROUTINE TRIDCOEF
!======================================================================
!       compute tridiagonal matrix coefficients
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8) :: I, J, K
       REAL(8)    :: ACC
       INTEGER(8) :: ILEV

        DO ILEV = 0, NLEV
            DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
                DO I = IIMG(ILEV, 1), IIMG(ILEV, 2)
                    ACC = -1.0_8 * (AW(I) + AE(I) + AS(1) + AN(1) + AB(K) + AF(K))
                    BET(I, 1, K) = 1.0_8 / ACC
                END DO
            END DO

            DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
                DO J = 2, N2M
                    DO I = IIMG(ILEV, 1), IIMG(ILEV, 2)
                        ACC = -1.0_8 * (AW(I) + AE(I) + AS(J) + AN(J) + AB(K) + AF(K))
                        GAM(I, J, K) = AN(J - 1) * BET(I, J - 1, K)
                        BET(I, J, K) = 1.0_8 / (ACC - AS(J) * GAM(I, J, K))
                    END DO
                END DO
            END DO
        END DO

        RETURN
        END SUBROUTINE TRIDCOEF