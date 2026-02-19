!=======================================================================
!
!     MAIN SOLVER OF POISSON EQUATION
!
!     AX=b
!     A: coefficient of discretized poisson equation
!     X: PHI: Pseudo pressure, Output of subroutine POISSON
!     b: DIVGSUM: OUTPUT of subroutine DIVGS.
!
!     x & z direction: Fourier transform
!     y-direction: TDMA
!
!     AK3,AK1: matrix coefficient (modified wavenumber)
!     N3MH,N1MH: The number of wavenumber index
!
!     Apr. 2010, J. Lee
!     Jun. 2017, J. Park
!     Feb. 2026, S. Lee (Several fixes)
!
!=======================================================================
      SUBROUTINE POISSON(PHI, DIVGSUM)
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE

       INTEGER(8)  :: I, J, K, JJP
       REAL(8)     :: PHI(0:N1,0:N2,0:N3), DIVGSUM(0:N1,0:N2,0:N3)
       REAL(8)     :: CN1, CN3
       REAL(8)     :: AJCREF, AJMREF, AJPREF
       COMPLEX(8)  :: CCAP(N3MH,N1,N2)
       COMPLEX(8)  :: CRHSREF, PHREF

       REAL(8), DIMENSION(:,:),     ALLOCATABLE :: AJC, AJM, AJP
       COMPLEX(8), DIMENSION(:,:),  ALLOCATABLE :: ZZZ, CRHS
       COMPLEX(8), DIMENSION(:),    ALLOCATABLE :: ZZZZ, ZZZZ_B, XXXX, XXXX_B

      CN1 = 1.0_8 / FLOAT(N1M)
      CN3 = 1.0_8 / FLOAT(N3M)

! --- FORWARD FOURIER TRANSFORM
!!$OMP PARALLEL private(ZZZ, ZZZZ, XXXX, XXXX_B, ZZZZ_B)
      ALLOCATE(ZZZ(N3M, N1M))
      ALLOCATE(ZZZZ(N3M), ZZZZ_B(N3M*2))
      ALLOCATE(XXXX(N1M), XXXX_B(N1M*2))
      CALL ZFFT1D(XXXX, N1M, 0, XXXX_B)
      CALL ZFFT1D(ZZZZ, N3M, 0, ZZZZ_B)
!!$OMP DO
      DO J = 1, N2M
          DO K = 1, N3M
              DO I = 1, N1M
                  ZZZ(K, I) = DIVGSUM(I, J, K)
              END DO
          END DO

          DO I = 1, N1M
              CALL ZFFT1D(ZZZ(1, I), N3M, -1, ZZZZ_B)
          END DO

          DO K = 1, N3MH
              DO I = 1, N1M
                  XXXX(I) = ZZZ(K, I)
              END DO
              CALL ZFFT1D(XXXX, N1M, -1, XXXX_B)
              DO I = 1, N1M
                  CCAP(K, I, J) = XXXX(I) 
              END DO
          END DO
      END DO
!!$OMP END DO
      DEALLOCATE(ZZZ, ZZZZ, ZZZZ_B, XXXX, XXXX_B)
!!$OMP END PARALLEL

! --- SOLVE TDMA MATRIX
!!$OMP PARALLEL private(AJM, AJP, AJC, CRHS, CRHSREF, AJCREF, AJMREF, AJPREF, PHREF)
      ALLOCATE(CRHS(N2, N1))
      ALLOCATE(AJM(N2, N1), AJP(N2, N1), AJC(N2, N1))
!!$OMP DO
      DO K = 1, N3MH
          DO I = 1, N1M
              DO J = 1, N2M
                  JJP = JPV(J)
                  AJM(J, I) = F2FYI(J) * C2CYI(J) * (1.0_8 - FIXJL(J))
                  AJP(J, I) = F2FYI(J) * C2CYI(JJP) * (1.0_8 - FIXJU(J))
                  AJC(J, I) = -((AJM(J, I) + AJP(J, I) + AK1(I) + AK3(K)) * &
                                (1.0_8 - FIXJL(J)) * (1.0_8 - FIXJU(J)) + &
                                (F2FYI(1) * C2CYI(2) + AK1(I) + AK3(K)) * FIXJL(J) + &
                                (F2FYI(N2M) * C2CYI(N2M) + AK1(I) + AK3(K)) * FIXJU(J))
                  CRHS(J, I) = CCAP(K, I, J)
              END DO
          END DO

          IF (K == 1) THEN
              CRHSREF = CRHS(N2M, 1)
              AJCREF = AJC(N2M, 1)
              AJMREF = AJM(N2M, 1)
              AJPREF = AJP(N2M, 1)
              CRHS(N2M, 1) = 0.0_8
              AJC(N2M, 1)  = 1.0_8
              AJM(N2M, 1)  = 0.0_8
              AJP(N2M, 1)  = 0.0_8
          END IF

          CALL CTRDIAG(AJM, AJC, AJP, CRHS, 1_8, N2M, CRHS, N1M)

          IF (K == 1) THEN
              PHREF = (-AJMREF * CRHS(N2M-1, 1) + CRHSREF) / AJCREF
              DO J = 1, N2M
                  CRHS(J, 1) = CRHS(J, 1) - PHREF
              END DO
              CRHS(N2M, 1) = 0.0_8
          END IF

          DO I = 1, N1M
              DO J = 1, N2M
                  CCAP(K, I, J) = CRHS(J, I)
              END DO
          END DO
      END DO
!!$OMP END DO
      DEALLOCATE(CRHS, AJM, AJP, AJC)
!!$OMP END PARALLEL

! --- INVERSE FOURIER TRANSFORM
!!$OMP PARALLEL private(ZZZ, ZZZZ, XXXX, XXXX_B, ZZZZ_B)
      ALLOCATE(ZZZ(N3M, N1M))
      ALLOCATE(ZZZZ(N3M), ZZZZ_B(N3M*2))
      ALLOCATE(XXXX(N1M), XXXX_B(N1M*2))
      CALL ZFFT1D(XXXX, N1M, 0, XXXX_B)
      CALL ZFFT1D(ZZZZ, N3M, 0, ZZZZ_B)
!!$OMP DO
      DO J = 1, N2M
          DO K = 1, N3MH
              DO I = 1, N1M
                  XXXX(I) = CCAP(K, I, J)
              END DO
              CALL ZFFT1D(XXXX, N1M, 1, XXXX_B)
              DO I = 1, N1M
                  ZZZ(K, I) = XXXX(I)
              END DO
          END DO

          DO I = 1, N1M
              DO K = N3MH + 1, N3M
                  ZZZ(K, I) = CONJG(ZZZ(N3M + 2 - K, I))
              END DO
          END DO

          DO I = 1, N1M
              CALL ZFFT1D(ZZZ(1, I), N3M, 1, ZZZZ_B)
              DO K = 1, N3M
                  PHI(I, J, K) = REAL(ZZZ(K, I), 8)
              END DO
          END DO
      END DO
!!$OMP END DO
      DEALLOCATE(ZZZ, ZZZZ, ZZZZ_B, XXXX, XXXX_B)
!!$OMP END PARALLEL

      RETURN
      END SUBROUTINE POISSON

!=======================================================================
      SUBROUTINE POISINIT
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8)  :: I, J, K, KK
       REAL(8)     :: PI
       REAL(8)     :: SDZIS, SDXIS

       CALL FTFT_ALLO

! --- DEFINE MODIFIED WAVENUMBERS
       PI = 2.0_8 * ASIN(1.0_8)

      DO K = 1, N3MH
          AI3(K) = FLOAT(K - 1) * 2.0_8 * PI
      END DO
      AI3(1) = 0.0_8
      
      DO I = 1, N1M
          AI1(I) = FLOAT(I - 1) * 2.0_8 * PI
      END DO
      AI1(1) = 0.0_8

      SDZIS = F2FZI(1)
      SDXIS = F2FXI(1)

      DO KK = 1, N3MH
          AK3(KK) = 2.0_8 * (1.0_8 - COS(AI3(KK) / FLOAT(N3M))) * SDZIS * SDZIS
      END DO
      
      DO KK = 1, N1MH
          AK1(KK) = 2.0_8 * (1.0_8 - COS(AI1(KK) / FLOAT(N1M))) * SDXIS * SDXIS
      END DO
      
      DO KK = N1M, N1MH + 1, -1
          AK1(KK) = AK1(N1M + 2 - KK)
      END DO

      RETURN
      END SUBROUTINE POISINIT

!=======================================================================
      SUBROUTINE CTRDIAG(A, B, C, R, NI, NF, UU, MF)
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER(8)  :: I, J, K
       INTEGER(8)  :: NI, NF, MF
       REAL(8)     :: A(N2, N1), B(N2, N1), C(N2, N1)
       COMPLEX(8)  :: R(N2, N1), UU(N2, N1)

       ALLOCATE(GAM(N2, N1, 1))
       ALLOCATE(BET(N1, 1, 1))

      DO I = 1, MF
          BET(I, 1, 1) = 1.0_8 / B(NI, I)
          UU(NI, I) = R(NI, I) * BET(I, 1, 1)
      END DO

      DO I = 1, MF
          DO J = NI + 1, NF
              GAM(J, I, 1) = C(J - 1, I) * BET(I, 1, 1)
              BET(I, 1, 1) = 1.0_8 / (B(J, I) - A(J, I) * GAM(J, I, 1))
              UU(J, I) = (R(J, I) - A(J, I) * UU(J - 1, I)) * BET(I, 1, 1)
          END DO
      END DO

      DO I = 1, MF
          DO J = NF - 1, NI, -1
              UU(J, I) = UU(J, I) - GAM(J + 1, I, 1) * UU(J + 1, I)
          END DO
      END DO

      DEALLOCATE(GAM, BET)

      RETURN
      END SUBROUTINE CTRDIAG