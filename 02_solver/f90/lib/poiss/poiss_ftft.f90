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
!
!=======================================================================
      SUBROUTINE POISSON(PHI,DIVGSUM)
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE

       INTEGER*8  :: I,J,K,JJP
       REAL*8      :: PHI(0:N1,0:N2,0:N3),DIVGSUM(0:N1,0:N2,0:N3)
       REAL*8      :: CN1,CN3
       REAL*8      :: AJCREF,AJMREF,AJPREF
       COMPLEX*16  :: CCAP(N3MH,N1,N2)
       COMPLEX*16  :: CRHSREF,PHREF

       REAL*8, DIMENSION (:,:),     ALLOCATABLE :: AJC,AJM,AJP
       COMPLEX*16, DIMENSION (:,:), ALLOCATABLE :: ZZZ,CRHS
       COMPLEX*16, DIMENSION (:),   ALLOCATABLE :: ZZZZ,ZZZZ_B,XXXX,XXXX_B

      CN1 = 1./FLOAT(N1M)
      CN3 = 1./FLOAT(N3M)

!     FORWARD FOURIER TRANSFORM
!!$OMP PARALLEL  &
!!$OMP private(ZZZ,ZZZZ,XXXX,XXXX_B,ZZZZ_B)
      allocate(ZZZ(N3M,N1M))
      allocate(ZZZZ(N3M),ZZZZ_B(N3M*2))
      allocate(XXXX(N1M),XXXX_B(N1M*2))
      CALL ZFFT1D(XXXX,N1M, 0,XXXX_B)
      CALL ZFFT1D(ZZZZ,N3M, 0,ZZZZ_B)
!!$OMP DO
      DO J=1,N2M

      DO K=1,N3M
      DO I=1,N1M
      ZZZ(K,I)=DIVGSUM(I,J,K)
      ENDDO
      ENDDO

      DO I=1,N1M
      CALL ZFFT1D(ZZZ(1,I),N3M,-1,ZZZZ_B)
      ENDDO

      DO K=1,N3MH
      DO I=1,N1M
      XXXX(I)=ZZZ(K,I)
      ENDDO
      CALL ZFFT1D(XXXX,N1M,-1,XXXX_B)
      DO I=1,N1M
      CCAP(K,I,J)=XXXX(I)!*CN1*CN3
      ENDDO
      ENDDO

      ENDDO
!!$OMP END DO
      deallocate(ZZZ,ZZZZ,ZZZZ_B,XXXX,XXXX_B)
!!$OMP END PARALLEL

!!!!!!!!!     SOLVE TDMA MATRIX
!!$OMP PARALLEL  &
!!$OMP private(AJM,AJP,AJC,CRHS)  &
!!$OMP private(CRHSREF,AJCREF,AJMREF,AJPREF,PHREF)
      allocate(CRHS(N2,N1))
      allocate(AJM(N2,N1),AJP(N2,N1),AJC(N2,N1))
!!$OMP DO
      DO K=1,N3MH
      DO I=1,N1M
      DO J=1,N2M
      JJP=JPV(J)
      AJM(J,I)=F2FYI(J)*C2CYI(J)*(1.-FIXJL(J))
      AJP(J,I)=F2FYI(J)*C2CYI(JJP)*(1.-FIXJU(J))
      AJC(J,I)=-((AJM(J,I)+AJP(J,I)+AK1(I)+AK3(K))             &
                *(1.-FIXJL(J))*(1.-FIXJU(J))                   &
                +(F2FYI(1)*C2CYI(2)+AK1(I)+AK3(K))*FIXJL(J)      &
                +(F2FYI(N2M)*C2CYI(N2M)+AK1(I)+AK3(K))*FIXJU(J))
      CRHS(J,I)=CCAP(K,I,J)
      ENDDO
      ENDDO
! Wavenumber index = 1: mean phi.
! in the poisson equation, only difference of phi is important.
! therefore, we should choose the reference phi.
! we choose referece phi as the phi at the upper wall.
      IF (K.EQ.1) THEN
      CRHSREF=CRHS(N2M,1)
      AJCREF=AJC(N2M,1)
      AJMREF=AJM(N2M,1)
      AJPREF=AJP(N2M,1)
      CRHS(N2M,1)= 0.
      AJC(N2M,1) = 1.
      AJM(N2M,1) = 0.
      AJP(N2M,1) = 0.
      ENDIF

      CALL CTRDIAG(AJM,AJC,AJP,CRHS,1,N2M,CRHS,N1M)

      IF (K.EQ.1) THEN
      PHREF=(-AJMREF*CRHS(N2M-1,1)+CRHSREF)/AJCREF
      DO J=1,N2M
      CRHS(J,1)=CRHS(J,1)-PHREF
      ENDDO
      CRHS(N2M,1)=0.
      ENDIF

      DO I=1,N1M
      DO J=1,N2M
      CCAP(K,I,J)=CRHS(J,I)
      ENDDO
      ENDDO
      ENDDO
!!$OMP END DO
      deallocate(CRHS)
      deallocate(AJM,AJP,AJC)
!!$OMP END PARALLEL


!     INVERSE FOURIER TRANSFORM
!!$OMP PARALLEL  &
!!$OMP private(ZZZ,ZZZZ,XXXX,XXXX_B,ZZZZ_B)
      ALLOCATE(ZZZ(N3M,N1M))
      ALLOCATE(ZZZZ(N3M),ZZZZ_B(N3M*2))
      ALLOCATE(XXXX(N1M),XXXX_B(N1M*2))
      CALL ZFFT1D(XXXX,N1M, 0,XXXX_B)
      CALL ZFFT1D(ZZZZ,N3M, 0,ZZZZ_B)
!!$OMP DO
      DO J=1,N2M

      DO K=1,N3MH
      DO I=1,N1M
       XXXX(I)=CCAP(K,I,J)
      ENDDO
      CALL ZFFT1D(XXXX,N1M,1,XXXX_B)
      DO I=1,N1M
       ZZZ(K,I)=XXXX(I)
      ENDDO
      ENDDO
      DO I=1,N1M
      DO K=N3MH+1,N3M
       ZZZ(K,I)=CONJG(ZZZ(N3M+2-K,I))
      ENDDO
      ENDDO

      DO I=1,N1M
       CALL ZFFT1D(ZZZ(1,I),N3M,1,ZZZZ_B)
      DO K=1,N3M
       PHI(I,J,K) = REAL(ZZZ(K,I))
      ENDDO
      ENDDO

      ENDDO
!!$OMP END DO
      DEALLOCATE(ZZZ,ZZZZ,ZZZZ_B,XXXX,XXXX_B)
!!$OMP END PARALLEL

      RETURN
      END
!=======================================================================
      SUBROUTINE POISINIT
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER*8  :: I,J,K,KK
       REAL*8  :: PI
       REAL*8  :: SDZIS,SDXIS

       CALL FTFT_ALLO

!     DEFINE MODIFIED WAVENUMBERS
       PI = 2.*ASIN(1.)

      DO 16 K=1,N3MH
   16 AI3(K)= (K-1)*2.*PI
      AI3(1)= 0.
      DO 17 I=1,N1M
   17 AI1(I)= (I-1)*2.*PI
      AI1(1)= 0.

      SDZIS=F2FZI(1)
      SDXIS=F2FXI(1)

      DO 2 KK=1,N3MH
    2 AK3(KK)=2.*(1.-COS(AI3(KK)/N3M))*SDZIS*SDZIS
      DO 3 KK=1,N1MH
    3 AK1(KK)=2.*(1.-COS(AI1(KK)/N1M))*SDXIS*SDXIS
      DO 4 KK=N1M,N1MH+1,-1
    4 AK1(KK)=AK1(N1M+2-KK)

      RETURN
      END

!=======================================================================
      SUBROUTINE CTRDIAG(A,B,C,R,NI,NF,UU,MF)
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER*8     :: I,J,K
       INTEGER*8     :: NI,NF,MF
       REAL*8        :: A(N2,N1),B(N2,N1),C(N2,N1)
       COMPLEX*16    :: R(N2,N1),UU(N2,N1)

       ALLOCATE(GAM(N2,N1,1))
       ALLOCATE(BET(N1,1,1))

      DO 10 I=1,MF
       BET(I,1,1)= 1./B(NI,I)
       UU(NI,I)=R(NI,I)*BET(I,1,1)
   10 CONTINUE
      DO 21 I=1,MF
      DO 11 J=NI+1,NF
       GAM(J,I,1)=C(J-1,I)*BET(I,1,1)
       BET(I,1,1)=1./(B(J,I)-A(J,I)*GAM(J,I,1))
       UU(J,I)=(R(J,I)-A(J,I)*UU(J-1,I))*BET(I,1,1)
   11 CONTINUE
   21 CONTINUE
      DO 22 I=1,MF
      DO 12 J=NF-1,NI,-1
       UU(J,I)=UU(J,I)-GAM(J+1,I,1)*UU(J+1,I)
   12 CONTINUE
   22 CONTINUE

      DEALLOCATE(GAM,BET)

      RETURN
      END
