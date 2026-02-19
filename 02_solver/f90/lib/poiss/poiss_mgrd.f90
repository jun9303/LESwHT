!=======================================================================
!
!     MAIN SOLVER OF POISSON EQUATION
!
!     AX=b
!     A: coefficient of discretized poisson equation
!     X: PHI: Pseudo pressure, Output of subroutine POISSON
!     b: DIVGSUM: OUTPUT of subroutine DIVGS.
!
!     Main algorithm: Preconditioned BICGSTAB (Biconjugate gradient stabilized method)
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
!
!=======================================================================
!     MAIN SOLVER OF POISSON EQUATION (3D Multigrid BiCGSTAB)
!=======================================================================
      SUBROUTINE POISINIT
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8) :: I, J, K, ILEV

       CALL MGRD_ALLO

       CALL COEFMG3

      OPEN(77, FILE='../output/ftr/mgrdresiduemax.dat')

      RETURN
      END SUBROUTINE POISINIT

!======================================================================
      SUBROUTINE POISSON(PHI, DIVGSUM)
!======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8) :: I, J, K, I_CG, IMAX, JMAX, KMAX
       REAL(8)    :: RHO_CG, RHO_OL, ALP_CG, W_CG, BET_CG, RV_TMP, ACC
       REAL(8)    :: PHIREF, RESM_CG, BTTT, UPPP
       
       REAL(8), DIMENSION(0:N1, 0:N2, 0:N3) :: PHI, DIVGSUM
       REAL(8), DIMENSION(0:N1, 0:N2, 0:N3) :: P_CG, S_CG, Y_CG, Z_CG, V_CG
       REAL(8), DIMENSION(N1, N2, N3)       :: T_CG, RES_CG

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
                                     -AW(I) * PHI(I - 1, J, K) - AE(I) * PHI(I + 1, J, K) &
                                     -AS(J) * PHI(I, J - 1, K) - AN(J) * PHI(I, J + 1, K) &
                                     -AB(K) * PHI(I, J, KMV(K)) - AF(K) * PHI(I, J, KPV(K)))
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
                       V_CG(I, J, K) = AW(I) * Y_CG(I - 1, J, K) + AE(I) * Y_CG(I + 1, J, K) &
                                     + AS(J) * Y_CG(I, J - 1, K) + AN(J) * Y_CG(I, J + 1, K) &
                                     + AB(K) * Y_CG(I, J, KMV(K)) + AF(K) * Y_CG(I, J, KPV(K)) &
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
                       T_CG(I, J, K) = AW(I) * Z_CG(I - 1, J, K) + AE(I) * Z_CG(I + 1, J, K) &
                                     + AS(J) * Z_CG(I, J - 1, K) + AN(J) * Z_CG(I, J + 1, K) &
                                     + AB(K) * Z_CG(I, J, KMV(K)) + AF(K) * Z_CG(I, J, KPV(K)) &
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
 102       FORMAT(I4, ES15.7, 3I4)

           IF (RESM_CG < TEST1) THEN
               WRITE(*, 300) I_CG, RESM_CG * DTCONST
               EXIT
           END IF
       END DO

       IF (RESM_CG >= TEST1) THEN
           PRINT*, '=== MUTLGRID : NOT CONVERGED ==='
           WRITE(*, 300) I_CG, RESM_CG * DTCONST
       END IF
 300   FORMAT('ICYC=', I10, '  RESMAX=', ES25.12)

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
 103   FORMAT(2F13.5, I5, ES15.7)

      RETURN
      END SUBROUTINE POISSON

!=======================================================================
       SUBROUTINE MG3D(PC, RHS, TEST, IOLDV)
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8) :: I, J, K, II, ILEV
       REAL(8)    :: TEST, SUMRES
       INTEGER(8) :: IOLDV
       REAL(8)    :: PC(0:N1, 0:N2, 0:N3), RHS(0:N1, 0:N2, 0:N3)
       REAL(8)    :: RESD(0:N1MD, N2M, N3MD), GGII(0:N1MD, 0:N2, N3MD)

       II = 0
       RESD = 0.0_8
       GGII = 0.0_8

       CALL TOPLEVEL3(0_8, PC, RHS, SUMRES, IOLDV, TEST, RESD, GGII)
       IF(SUMRES < TEST) GOTO 2000

       DO II = 1, 10
        DO ILEV = NLEV - 1, 1, -1
         CALL RELAX3(ILEV, 0.0_8, 1_8, RESD, GGII)
         CALL GODOWN3(ILEV, RESD, GGII)
        END DO

        CALL RELAX3(0_8, 0.0_8, NBLI, RESD, GGII)

        DO ILEV = 0, NLEV - 2
         CALL GOUP3(ILEV, RESD, GGII)
         CALL RELAX3(ILEV + 1, 1.0_8, 1_8, RESD, GGII)
        END DO

        CALL TOPLEVEL3(1_8, PC, RHS, SUMRES, 1_8, TEST, RESD, GGII)
        IF(SUMRES < TEST) GOTO 2000
       END DO

 2000  CONTINUE

      RETURN
      END SUBROUTINE MG3D

!=======================================================================
      SUBROUTINE TOPLEVEL3(ID, PC, RHS, SUMRES, IOLDV, TEST, RESD, GGII)
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
        INTEGER(8) :: I, J, K, IC, KC
        INTEGER(8) :: ID, IOLDV
        REAL(8)    :: SUMRES, TEST, TT, ACC
        REAL(8)    :: PC(0:N1, 0:N2, 0:N3), RHS(0:N1, 0:N2, 0:N3)
        REAL(8)    :: RESD(0:N1MD, N2M, N3MD), GGII(0:N1MD, 0:N2, N3MD)

       IF(ID == 1) THEN
        KC = 0
        DO K = KKMG(NLEV - 1, 1), KKMG(NLEV - 1, 2)
        KC = KC + 2
        IC = 0
        DO I = IIMG(NLEV - 1, 1), IIMG(NLEV - 1, 2) - 1
        IC = IC + 2
        DO J = 1, N2M
         PC(IC, J, KC) = PC(IC, J, KC) + &
             CORZ1(K) * (COI1(I) * GGII(I, J, K) + COI2(I) * GGII(I + 1, J, K)) + &
             CORZ2(K) * (COI1(I) * GGII(I, J, K + 1) + COI2(I) * GGII(I + 1, J, K + 1))
        END DO
        END DO

        I = IIMG(NLEV - 1, 2)
        IC = IC + 2
        DO J = 1, N2M
         PC(IC, J, KC) = PC(IC, J, KC) + CORZ1(K) * COI1(I) * GGII(I, J, K) + CORZ2(K) * COI1(I) * GGII(I, J, K + 1)
        END DO
        END DO
       END IF

      DO K = 1, N3M
      DO J = 1, N2M
      DO I = 1, N1M
       IF (IOLDV == 1) THEN
           ACC = AW(I) + AE(I) + AS(J) + AN(J) + AB(K) + AF(K)
           GGII(I, J, K) = (RHS(I, J, K) + AW(I) * PC(I - 1, J, K) + AE(I) * PC(I + 1, J, K) &
                                         + AS(J) * PC(I, J - 1, K) + AN(J) * PC(I, J + 1, K) &
                                         + AB(K) * PC(I, J, KMV(K)) + AF(K) * PC(I, J, KPV(K))) / ACC
       ELSE
           ACC = AW(I) + AE(I) + AS(J) + AN(J) + AB(K) + AF(K)
           GGII(I, J, K) = RHS(I, J, K) / ACC
       END IF
      END DO
      END DO
      END DO

      SUMRES = 0.0_8

      DO K = 1, N3M
      DO J = 1, N2M
      DO I = 1, N1M
       ACC = AW(I) + AE(I) + AS(J) + AN(J) + AB(K) + AF(K)
       RESD(I, J, K) = RHS(I, J, K) - AW(I) * GGII(I - 1, J, K) - AE(I) * GGII(I + 1, J, K) &
                                    - AS(J) * GGII(I, J - 1, K) - AN(J) * GGII(I, J + 1, K) &
                                    - AB(K) * GGII(I, J, KMV(K)) - AF(K) * GGII(I, J, KPV(K)) &
                                    - ACC * GGII(I, J, K)
       SUMRES = MAX(SUMRES, ABS(RESD(I, J, K)))
       PC(I, J, K) = GGII(I, J, K)
      END DO
      END DO
      END DO

      IF(SUMRES < TEST) RETURN
      IF(ID == 2) RETURN

      KC = -1
      DO K = KKMG(NLEV - 1, 1), KKMG(NLEV - 1, 2)
      KC = KC + 2
      IC = -1
      DO I = IIMG(NLEV - 1, 1), IIMG(NLEV - 1, 2)
      IC = IC + 2
      DO J = 1, N2M
      RESD(I, J, K) = CORZ1(K) * (RESD(IC, J, KC) * COR1(I) + RESD(IC + 1, J, KC) * COR2(I)) + &
                      CORZ2(K) * (RESD(IC, J, KC + 1) * COR1(I) + RESD(IC + 1, J, KC + 1) * COR2(I))
      END DO
      END DO
      END DO

      RETURN
      END SUBROUTINE TOPLEVEL3

!=======================================================================
      SUBROUTINE GOUP3(ILEV, RESD, GGII)
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
        INTEGER(8) :: I, J, K, ILEV, IBGH, KBGH
        REAL(8)    :: RESD(0:N1MD, N2M, N3MD), GGII(0:N1MD, 0:N2, N3MD)

      KBGH = KKMG(ILEV + 1, 1) - 1
      DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
      KBGH = KBGH + 2
      IBGH = IIMG(ILEV + 1, 1) - 1
      DO I = IIMG(ILEV, 1), IIMG(ILEV, 2) - 1
      IBGH = IBGH + 2
      DO J = 1, N2M
       GGII(IBGH, J, KBGH) = GGII(IBGH, J, KBGH) + &
           COIZ1(K) * (COI1(I) * GGII(I, J, K) + COI2(I) * GGII(I + 1, J, K)) + &
           COIZ2(K) * (COI1(I) * GGII(I, J, K + 1) + COI2(I) * GGII(I + 1, J, K + 1))
      END DO
      END DO

      I = IIMG(ILEV, 2)
      IBGH = IBGH + 2
      DO J = 1, N2M
       GGII(IBGH, J, KBGH) = GGII(IBGH, J, KBGH) + COIZ1(K) * COI1(I) * GGII(I, J, K) + COIZ2(K) * COI1(I) * GGII(I, J, K + 1)
      END DO
      END DO

      RETURN
      END SUBROUTINE GOUP3

!=======================================================================
      SUBROUTINE GODOWN3(ILEV, RESD, GGII)
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
        INTEGER(8) :: I, J, K, ILEV, IBG, KBG
        REAL(8)    :: ACC
        REAL(8)    :: RESD(0:N1MD, N2M, N3MD), GGII(0:N1MD, 0:N2, N3MD)

      DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
      DO J = 1, N2M
      DO I = IIMG(ILEV, 1), IIMG(ILEV, 2)
      ACC = AW(I) + AE(I) + AS(J) + AN(J) + AB(K) + AF(K)
      RESD(I, J, K) = RESD(I, J, K) - AW(I) * GGII(I - 1, J, K) - AE(I) * GGII(I + 1, J, K) &
                                    - AS(J) * GGII(I, J - 1, K) - AN(J) * GGII(I, J + 1, K) &
                                    - AB(K) * GGII(I, J, K - 1) - AF(K) * GGII(I, J, K + 1) &
                                    - ACC * GGII(I, J, K)
      END DO
      END DO
      END DO

      KBG = KKMG(ILEV, 1) - 2
      DO K = KKMG(ILEV - 1, 1), KKMG(ILEV - 1, 2)
      KBG = KBG + 2
      IBG = IIMG(ILEV, 1) - 2
      DO I = IIMG(ILEV - 1, 1), IIMG(ILEV - 1, 2)
      IBG = IBG + 2
      DO J = 1, N2M
      RESD(I, J, K) = CORZ1(K) * (RESD(IBG, J, KBG) * COR1(I) + RESD(IBG + 1, J, KBG) * COR2(I)) + &
                      CORZ2(K) * (RESD(IBG, J, KBG + 1) * COR1(I) + RESD(IBG + 1, J, KBG + 1) * COR2(I))
      END DO
      END DO
      END DO

      RETURN
      END SUBROUTINE GODOWN3

!=======================================================================
      SUBROUTINE RELAX3(ILEV, OLDV, IITER, RESD, GGII)
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
        INTEGER(8) :: I, J, K, ILEV, IITER, II
        REAL(8)    :: RESD(0:N1MD, N2M, N3MD), GGII(0:N1MD, 0:N2, N3MD)
        REAL(8)    :: OLDV, ACC

      DO II = 1, IITER
       DO K = KKMG(ILEV, 1), KKMG(ILEV, 2)
       DO J = 1, N2M
       DO I = IIMG(ILEV, 1), IIMG(ILEV, 2)
       ACC = AW(I) + AE(I) + AS(J) + AN(J) + AB(K) + AF(K)
       GGII(I, J, K) = (RESD(I, J, K) + OLDV * (AW(I) * GGII(I - 1, J, K) + AE(I) * GGII(I + 1, J, K) &
                                              + AS(J) * GGII(I, J - 1, K) + AN(J) * GGII(I, J + 1, K) &
                                              + AB(K) * GGII(I, J, K - 1) + AF(K) * GGII(I, J, K + 1))) / ACC
       END DO
       END DO
       END DO
      END DO

      RETURN
      END SUBROUTINE RELAX3

!=======================================================================
      SUBROUTINE COEFMG3
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8) :: I, J, K, ILEV
       INTEGER(8) :: MINROW, IBG, IEND, ISP, IC, IBGH
       INTEGER(8) :: MINROWZ, KBG, KEND, KSP, KC, KBGH
       REAL(8)    :: IWZ(N1MD), IEZ(N1MD), KBZ(N3MD), KFZ(N3MD)
       REAL(8)    :: XMPM(0:N1MD), ZMPM(0:N3MD)

       MINROW = N1M / (2**NLEV)
       XMPM = 0.0_8

       LEVHALF = NINT(NLEV / 2.0_8)
       IIMG(NLEV, 1) = 1
       IIMG(NLEV, 2) = N1M
       IIMG(NLEV, 3) = N1M

       DO I = 1, N1M
        XMPM(I) = XMP(I)
        IWZ(I) = 1.0_8 - 1.0_8 / FLOAT(I)
        IEZ(I) = 1.0_8 - FLOAT(I) / FLOAT(N1M)
       END DO

      DO ILEV = NLEV - 1, 0, -1
       IIMG(ILEV, 1) = IIMG(ILEV + 1, 1) + IIMG(ILEV + 1, 3)
       IIMG(ILEV, 3) = MINROW * (2**ILEV)
       IIMG(ILEV, 2) = IIMG(ILEV, 1) + IIMG(ILEV, 3) - 1

       IBG  = IIMG(ILEV, 1)
       IEND = IIMG(ILEV, 2)
       ISP  = 2**(NLEV - ILEV)

       IC = 0
       DO I = IBG, IEND
        IC = IC + 1
        XMPM(I) = 0.5_8 * (X(IC * ISP + 1) + X((IC - 1) * ISP + 1))
        IWZ(I)  = 1.0_8 - FLOAT(IBG) / FLOAT(I)
        IEZ(I)  = 1.0_8 - FLOAT(I) / FLOAT(IEND)
       END DO

       IC = 0
       DO I = IBG, IEND
        IC = IC + 1
        AW(I) = IWZ(I) / ((XMPM(I) - XMPM(I - 1)) * (X(IC * ISP + 1) - X((IC - 1) * ISP + 1)))
        AE(I) = IEZ(I) / ((XMPM(I + 1) - XMPM(I)) * (X(IC * ISP + 1) - X((IC - 1) * ISP + 1)))
       END DO
      END DO

       DO ILEV = NLEV, 1, -1
        IBGH = IIMG(ILEV, 1)
        DO I = IIMG(ILEV - 1, 1), IIMG(ILEV - 1, 2)
        COR1(I) = (XMPM(IBGH + 1) - XMPM(I)) / (XMPM(IBGH + 1) - XMPM(IBGH))
        COR2(I) = 1.0_8 - COR1(I)
        IBGH = IBGH + 2
        END DO
       END DO

       DO ILEV = 0, NLEV - 1
        IBGH = IIMG(ILEV + 1, 1) + 1
        DO I = IIMG(ILEV, 1), IIMG(ILEV, 2) - 1
        COI1(I) = (XMPM(I + 1) - XMPM(IBGH)) / (XMPM(I + 1) - XMPM(I))
        COI2(I) = 1.0_8 - COI1(I)
        IBGH = IBGH + 2
        END DO
        I = IIMG(ILEV, 2)
        COI1(I) = 1.0_8
       END DO

       MINROWZ = N3M / (2**NLEV)
       ZMPM = 0.0_8

       KKMG(NLEV, 1) = 1
       KKMG(NLEV, 2) = N3M
       KKMG(NLEV, 3) = N3M

       DO K = 1, N3M
        ZMPM(K) = ZMP(K)
        KBZ(K)  = 1.0_8 - 1.0_8 / FLOAT(K)
        KFZ(K)  = 1.0_8 - FLOAT(K) / FLOAT(N3M)
       END DO

      DO ILEV = NLEV - 1, 0, -1
       KKMG(ILEV, 1) = KKMG(ILEV + 1, 1) + KKMG(ILEV + 1, 3)
       KKMG(ILEV, 3) = MINROWZ * (2**ILEV)
       KKMG(ILEV, 2) = KKMG(ILEV, 1) + KKMG(ILEV, 3) - 1

       KBG  = KKMG(ILEV, 1)
       KEND = KKMG(ILEV, 2)
       KSP  = 2**(NLEV - ILEV)

       KC = 0
       DO K = KBG, KEND
        KC = KC + 1
        ZMPM(K) = 0.5_8 * (Z(KC * KSP + 1) + Z((KC - 1) * KSP + 1))
        KBZ(K)  = 1.0_8 - FLOAT(KBG) / FLOAT(K)
        KFZ(K)  = 1.0_8 - FLOAT(K) / FLOAT(KEND)
       END DO

       KC = 0
       DO K = KBG, KEND
        KC = KC + 1
        AB(K) = KBZ(K) / ((ZMPM(K) - ZMPM(K - 1)) * (Z(KC * KSP + 1) - Z((KC - 1) * KSP + 1)))
        AF(K) = KFZ(K) / ((ZMPM(K + 1) - ZMPM(K)) * (Z(KC * KSP + 1) - Z((KC - 1) * KSP + 1)))
       END DO
      END DO

       DO ILEV = NLEV, 1, -1
        KBGH = KKMG(ILEV, 1)
        DO K = KKMG(ILEV - 1, 1), KKMG(ILEV - 1, 2)
        CORZ1(K) = (ZMPM(KBGH + 1) - ZMPM(K)) / (ZMPM(KBGH + 1) - ZMPM(KBGH))
        CORZ2(K) = 1.0_8 - CORZ1(K)
        KBGH = KBGH + 2
        END DO
       END DO

       DO ILEV = 0, NLEV - 1
        KBGH = KKMG(ILEV + 1, 1) + 1
        DO K = KKMG(ILEV, 1), KKMG(ILEV, 2) - 1
        COIZ1(K) = (ZMPM(K + 1) - ZMPM(KBGH)) / (ZMPM(K + 1) - ZMPM(K))
        COIZ2(K) = 1.0_8 - COIZ1(K)
        KBGH = KBGH + 2
        END DO
        K = KKMG(ILEV, 2)
        COIZ1(K) = 1.0_8
       END DO

      RETURN
      END SUBROUTINE COEFMG3