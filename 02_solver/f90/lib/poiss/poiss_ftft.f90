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

!
!     FFTE: A FAST FOURIER TRANSFORM PACKAGE
!
!     (C) COPYRIGHT SOFTWARE, 2000-2004, 2008-2011, ALL RIGHTS RESERVED
!                BY
!         DAISUKE TAKAHASHI
!         FACULTY OF ENGINEERING, INFORMATION AND SYSTEMS
!         UNIVERSITY OF TSUKUBA
!         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
!         E-MAIL: daisuke@cs.tsukuba.ac.jp
!
!
!     1-D COMPLEX FFT ROUTINE
!
!     FORTRAN77 SOURCE PROGRAM
!
!     CALL ZFFT1D(A,N,IOPT,B)
!
!     A(N) IS COMPLEX INPUT/OUTPUT VECTOR (COMPLEX*16)
!     B(N*2) IS WORK/COEFFICIENT VECTOR (COMPLEX*16)
!     N IS THE LENGTH OF THE TRANSFORMS (INTEGER*4)
!       -----------------------------------
!         N = (2**IP) * (3**IQ) * (5**IR)
!       -----------------------------------
!     IOPT = 0 FOR INITIALIZING THE COEFFICIENTS (INTEGER*4)
!          = -1 FOR FORWARD TRANSFORM
!          = +1 FOR INVERSE TRANSFORM
!
!     WRITTEN BY DAISUKE TAKAHASHI
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE ZFFT1D(A,N,IOPT,B)

      IMPLICIT REAL*8 (A-H,O-Z)
!     FFTE: A FAST FOURIER TRANSFORM PACKAGE
!
!     (C) COPYRIGHT SOFTWARE, 2000-2004, 2008-2014, ALL RIGHTS RESERVED
!                BY
!         DAISUKE TAKAHASHI
!         FACULTY OF ENGINEERING, INFORMATION AND SYSTEMS
!         UNIVERSITY OF TSUKUBA
!         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
!         E-MAIL: daisuke@cs.tsukuba.ac.jp
!
!     HEADER FILE FOR PARAMETERS
!
!     FORTRAN77 SOURCE PROGRAM
!
!     WRITTEN BY DAISUKE TAKAHASHI
!
! The maximum supported number of processors is 65536.
      PARAMETER (MAXNPU=65536)
! The maximum supported 2-D transform length is 65536.
      PARAMETER (NDA2=65536)
! The maximum supported 3-D transform length is 4096.
      PARAMETER (NDA3=4096)
! The parameter NBLK is a blocking parameter.
      PARAMETER (NBLK=16)
! The parameter NB is a blocking parameter for NVIDIA GPUs.
      parameter (NB=128)
! The parameter NP is a padding parameter to avoid cache conflicts in
! the FFT routines.
      PARAMETER (NP=8)
! Size of L2 cache
      PARAMETER (L2SIZE=2097152)
      COMPLEX*16 A(*),B(*)
      COMPLEX*16 C((NDA2+NP)*NBLK),D(NDA2)
      COMPLEX*16 WX(NDA2),WY(NDA2)
      DIMENSION IP(3),LNX(3),LNY(3)
      SAVE WX,WY
!
      CALL FACTOR(N,IP)
!
      IF (IOPT .EQ. 1) THEN
!$OMP PARALLEL DO
!DIR$ VECTOR ALIGNED
        DO 10 I=1,N
          A(I)=DCONJG(A(I))
   10   CONTINUE
      END IF
!
      IF (N .LE. (L2SIZE/16)/3) THEN
        IF (IOPT .EQ. 0) THEN
          CALL SETTBL(B(N+1),N)
          RETURN
        END IF
!
        CALL FFT235(A,B,B(N+1),N,IP)
      ELSE
        CALL GETNXNY(N,NX,NY)
        CALL FACTOR(NX,LNX)
        CALL FACTOR(NY,LNY)
!
        IF (IOPT .EQ. 0) THEN
          CALL SETTBL(WX,NX)
          CALL SETTBL(WY,NY)
          CALL SETTBL2(B(N+1),NX,NY)
          RETURN
        END IF
!
!$OMP PARALLEL PRIVATE(C,D)
        CALL ZFFT1D0(A,A,B,C,C,D,WX,WY,B(N+1),NX,NY,LNX,LNY)
!$OMP END PARALLEL
      END IF
!
      IF (IOPT .EQ. 1) THEN
        DN=1.0D0/DBLE(N)
!$OMP PARALLEL DO
!DIR$ VECTOR ALIGNED
        DO 20 I=1,N
          A(I)=DCONJG(A(I))*DN
   20   CONTINUE
      END IF
      RETURN
      END
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE ZFFT1D0(A,AYX,B,CX,CY,D,WX,WY,W,NX,NY,LNX,LNY)
      IMPLICIT REAL*8 (A-H,O-Z)

!     FFTE: A FAST FOURIER TRANSFORM PACKAGE
!
!     (C) COPYRIGHT SOFTWARE, 2000-2004, 2008-2014, ALL RIGHTS RESERVED
!                BY
!         DAISUKE TAKAHASHI
!         FACULTY OF ENGINEERING, INFORMATION AND SYSTEMS
!         UNIVERSITY OF TSUKUBA
!         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
!         E-MAIL: daisuke@cs.tsukuba.ac.jp
!
!
!     HEADER FILE FOR PARAMETERS
!
!     FORTRAN77 SOURCE PROGRAM
!
!     WRITTEN BY DAISUKE TAKAHASHI
!
! The maximum supported number of processors is 65536.
      PARAMETER (MAXNPU=65536)
! The maximum supported 2-D transform length is 65536.
      PARAMETER (NDA2=65536)
! The maximum supported 3-D transform length is 4096.
      PARAMETER (NDA3=4096)
! The parameter NBLK is a blocking parameter.
      PARAMETER (NBLK=16)
! The parameter NB is a blocking parameter for NVIDIA GPUs.
      parameter (NB=128)
! The parameter NP is a padding parameter to avoid cache conflicts in
! the FFT routines.
      PARAMETER (NP=8)
! Size of L2 cache
      PARAMETER (L2SIZE=2097152)


      COMPLEX*16 A(NX,*),AYX(NY,*),B(NX,*)
      COMPLEX*16 CX(NX+NP,*),CY(NY+NP,*),D(*)
      COMPLEX*16 WX(*),WY(*),W(NX,*)
      DIMENSION LNX(*),LNY(*)
!
!$OMP DO
      DO 70 II=1,NX,NBLK
        DO 30 JJ=1,NY,NBLK
          DO 20 I=II,MIN0(II+NBLK-1,NX)
!DIR$ VECTOR ALIGNED
            DO 10 J=JJ,MIN0(JJ+NBLK-1,NY)
              CY(J,I-II+1)=A(I,J)
   10       CONTINUE
   20     CONTINUE
   30   CONTINUE
        DO 40 I=II,MIN0(II+NBLK-1,NX)
          CALL FFT235(CY(1,I-II+1),D,WY,NY,LNY)
   40   CONTINUE
        DO 60 J=1,NY
!DIR$ VECTOR ALIGNED
          DO 50 I=II,MIN0(II+NBLK-1,NX)
            B(I,J)=CY(J,I-II+1)
   50     CONTINUE
   60   CONTINUE
   70 CONTINUE
!$OMP DO
      DO 120 JJ=1,NY,NBLK
        DO 90 J=JJ,MIN0(JJ+NBLK-1,NY)
!DIR$ VECTOR ALIGNED
          DO 80 I=1,NX
            CX(I,J-JJ+1)=B(I,J)*W(I,J)
   80     CONTINUE
          CALL FFT235(CX(1,J-JJ+1),D,WX,NX,LNX)
   90   CONTINUE
        DO 110 I=1,NX
!DIR$ VECTOR ALIGNED
          DO 100 J=JJ,MIN0(JJ+NBLK-1,NY)
            AYX(J,I)=CX(I,J-JJ+1)
  100     CONTINUE
  110   CONTINUE
  120 CONTINUE
      RETURN
      END
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!
!     FFTE: A FAST FOURIER TRANSFORM PACKAGE
!
!     (C) COPYRIGHT SOFTWARE, 2000-2004, 2008-2014, ALL RIGHTS RESERVED
!                BY
!         DAISUKE TAKAHASHI
!         FACULTY OF ENGINEERING, INFORMATION AND SYSTEMS
!         UNIVERSITY OF TSUKUBA
!         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
!         E-MAIL: daisuke@cs.tsukuba.ac.jp
!
!
!     FACTORIZATION ROUTINE
!
!     FORTRAN77 SOURCE PROGRAM
!
!     WRITTEN BY DAISUKE TAKAHASHI
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE FACTOR(N,IP)
      DIMENSION IP(*)
!
      IP(1)=0
      IP(2)=0
      IP(3)=0
      N2=N
      IF (MOD(N,2) .NE. 0 .AND. MOD(N,3) .NE. 0 .AND. MOD(N,5) .NE. 0) RETURN
   10 IF (N2 .LE. 1) RETURN
      IF (MOD(N2,2) .EQ. 0) THEN
        IP(1)=IP(1)+1
        N2=N2/2
        GO TO 10
      ELSE IF (MOD(N2,3) .EQ. 0) THEN
        IP(2)=IP(2)+1
        N2=N2/3
        GO TO 10
      ELSE IF (MOD(N2,5) .EQ. 0) THEN
        IP(3)=IP(3)+1
        N2=N2/5
        GO TO 10
      END IF
      RETURN
      END
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE GETNXNY(N,NX,NY)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IP(3),LNX(3),LNY(3)
!
      ISQRTN=IDINT(DSQRT(DBLE(N)))
      CALL FACTOR(N,IP)
      DO 10 I=1,3
        LNX(I)=0
   10 CONTINUE
      IRES=ISQRTN
      DO 40 K=0,(IP(3)+1)/2
        DO 30 J=0,(IP(2)+1)/2
          DO 20 I=0,(IP(1)+1)/2
            NX=(2**I)*(3**J)*(5**K)
            IF (NX .LE. ISQRTN) THEN
              IRES2=ISQRTN-NX
              IF (IRES2 .LT. IRES) THEN
                LNX(1)=I
                LNX(2)=J
                LNX(3)=K
                IRES=IRES2
              END IF
            END IF
   20     CONTINUE
   30   CONTINUE
   40 CONTINUE
      DO 50 I=1,3
        LNY(I)=IP(I)-LNX(I)
   50 CONTINUE
      NX=(2**LNX(1))*(3**LNX(2))*(5**LNX(3))
      NY=(2**LNY(1))*(3**LNY(2))*(5**LNY(3))
      RETURN
      END
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!     FFTE: A FAST FOURIER TRANSFORM PACKAGE
!
!     (C) COPYRIGHT SOFTWARE, 2000-2004, 2008-2014, ALL RIGHTS RESERVED
!                BY
!         DAISUKE TAKAHASHI
!         FACULTY OF ENGINEERING, INFORMATION AND SYSTEMS
!         UNIVERSITY OF TSUKUBA
!         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
!         E-MAIL: daisuke@cs.tsukuba.ac.jp
!
!
!     RADIX-2, 3, 4, 5 AND 8 FFT ROUTINE
!
!     FORTRAN77 SOURCE PROGRAM
!
!     WRITTEN BY DAISUKE TAKAHASHI
!
      SUBROUTINE FFT235(A,B,W,N,IP)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*),W(*)
      DIMENSION IP(*)
!
      IF (IP(1) .NE. 1) THEN
        KP4=2-MOD(IP(1)+2,3)
        KP8=(IP(1)-KP4)/3
      ELSE
        KP4=0
        KP8=0
      END IF
!
      KEY=1
      J=1
      L=N
      M=1
      DO 10 K=1,KP8
        L=L/8
        IF (L .GE. 2) THEN
          IF (KEY .GE. 0) THEN
            CALL FFT8(A,B,W(J),M,L)
          ELSE
            CALL FFT8(B,A,W(J),M,L)
          END IF
          KEY=-KEY
        ELSE
          IF (KEY .GE. 0) THEN
            CALL FFT8(A,A,W(J),M,L)
          ELSE
            CALL FFT8(B,A,W(J),M,L)
          END IF
        END IF
        M=M*8
        J=J+L*7
   10 CONTINUE
      DO 20 K=1,IP(3)
        L=L/5
        IF (L .GE. 2) THEN
          IF (KEY .GE. 0) THEN
            CALL FFT5(A,B,W(J),M,L)
          ELSE
            CALL FFT5(B,A,W(J),M,L)
          END IF
          KEY=-KEY
        ELSE
          IF (KEY .GE. 0) THEN
            CALL FFT5(A,A,W(J),M,L)
          ELSE
            CALL FFT5(B,A,W(J),M,L)
          END IF
        END IF
        M=M*5
        J=J+L*4
   20 CONTINUE
      DO 30 K=1,KP4
        L=L/4
        IF (L .GE. 2) THEN
          IF (KEY .GE. 0) THEN
            CALL FFT4(A,B,W(J),M,L)
          ELSE
            CALL FFT4(B,A,W(J),M,L)
          END IF
          KEY=-KEY
        ELSE
          IF (KEY .GE. 0) THEN
            CALL FFT4(A,A,W(J),M,L)
          ELSE
            CALL FFT4(B,A,W(J),M,L)
          END IF
        END IF
        M=M*4
        J=J+L*3
   30 CONTINUE
      DO 40 K=1,IP(2)
        L=L/3
        IF (L .GE. 2) THEN
          IF (KEY .GE. 0) THEN
            CALL FFT3(A,B,W(J),M,L)
          ELSE
            CALL FFT3(B,A,W(J),M,L)
          END IF
          KEY=-KEY
        ELSE
          IF (KEY .GE. 0) THEN
            CALL FFT3(A,A,W(J),M,L)
          ELSE
            CALL FFT3(B,A,W(J),M,L)
          END IF
        END IF
        M=M*3
        J=J+L*2
   40 CONTINUE
      IF (IP(1) .EQ. 1) THEN
        IF (KEY .GE. 0) THEN
          CALL FFT2(A,A,M)
        ELSE
          CALL FFT2(B,A,M)
        END IF
      END IF
      RETURN
      END
      SUBROUTINE FFT3(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*),W(*)
!
      IF (M .EQ. 1) THEN
        CALL FFT3A(A,B,W,L)
      ELSE
        CALL FFT3B(A,B,W,M,L)
      END IF
      RETURN
      END
      SUBROUTINE FFT4(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*),W(*)
!
      IF (M .EQ. 1) THEN
        CALL FFT4A(A,B,W,L)
      ELSE
        CALL FFT4B(A,B,W,M,L)
      END IF
      RETURN
      END
      SUBROUTINE FFT5(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*),W(*)
!
      IF (M .EQ. 1) THEN
        CALL FFT5A(A,B,W,L)
      ELSE
        CALL FFT5B(A,B,W,M,L)
      END IF
      RETURN
      END
      SUBROUTINE FFT8(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*),W(*)
!
      IF (M .EQ. 1) THEN
        CALL FFT8A(A,B,W,L)
      ELSE
        CALL FFT8B(A,B,W,M,L)
      END IF
      RETURN
      END
      SUBROUTINE SETTBL(W,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 W(*)
      DIMENSION IP(3)
!
      CALL FACTOR(N,IP)
!
      IF (IP(1) .NE. 1) THEN
        KP4=2-MOD(IP(1)+2,3)
        KP8=(IP(1)-KP4)/3
      ELSE
        KP4=0
        KP8=0
      END IF
!
      J=1
      L=N
      DO 10 K=1,KP8
        L=L/8
        CALL SETTBL0(W(J),8,L)
        J=J+L*7
   10 CONTINUE
      DO 20 K=1,IP(3)
        L=L/5
        CALL SETTBL0(W(J),5,L)
        J=J+L*4
   20 CONTINUE
      DO 30 K=1,KP4
        L=L/4
        CALL SETTBL0(W(J),4,L)
        J=J+L*3
   30 CONTINUE
      DO 40 K=1,IP(2)
        L=L/3
        CALL SETTBL0(W(J),3,L)
        J=J+L*2
   40 CONTINUE
      RETURN
      END
      SUBROUTINE SETTBL0(W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 W(M-1,*)
!
      PI2=8.0D0*DATAN(1.0D0)
      PX=-PI2/(DBLE(M)*DBLE(L))
      DO 20 J=1,L
!DIR$ VECTOR ALIGNED
        DO 10 I=1,M-1
          TEMP=PX*DBLE(I)*DBLE(J-1)
          W(I,J)=DCMPLX(DCOS(TEMP),DSIN(TEMP))
   10   CONTINUE
   20 CONTINUE
      RETURN
      END
      SUBROUTINE SETTBL2(W,NX,NY)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 W(NX,*)
!
      PI2=8.0D0*DATAN(1.0D0)
      PX=-PI2/(DBLE(NX)*DBLE(NY))
!$OMP PARALLEL DO PRIVATE(TEMP)
      DO 20 J=1,NY
!DIR$ VECTOR ALIGNED
        DO 10 I=1,NX
          TEMP=PX*DBLE(I-1)*DBLE(J-1)
          W(I,J)=DCMPLX(DCOS(TEMP),DSIN(TEMP))
   10   CONTINUE
   20 CONTINUE
      RETURN
      END


!
!     FFTE: A FAST FOURIER TRANSFORM PACKAGE
!
!     (C) COPYRIGHT SOFTWARE, 2000-2004, 2008-2011, ALL RIGHTS RESERVED
!                BY
!         DAISUKE TAKAHASHI
!         FACULTY OF ENGINEERING, INFORMATION AND SYSTEMS
!         UNIVERSITY OF TSUKUBA
!         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
!         E-MAIL: daisuke@cs.tsukuba.ac.jp
!
!
!     RADIX-2, 3, 4, 5 AND 8 FFT KERNEL ROUTINE
!
!     FORTRAN77 SOURCE PROGRAM
!
!     WRITTEN BY DAISUKE TAKAHASHI
!
      SUBROUTINE FFT2(A,B,M)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(M,*),B(M,*)
      COMPLEX*16 C0,C1
!
!DIR$ VECTOR ALIGNED
      DO 10 I=1,M
        C0=A(I,1)
        C1=A(I,2)
        B(I,1)=C0+C1
        B(I,2)=C0-C1
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FFT3A(A,B,W,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(L,*),B(3,*),W(2,*)
      COMPLEX*16 C0,C1,C2,D0,D1,D2,W1,W2
      DATA C31/0.86602540378443865D0/C32/0.5D0/
!
!DIR$ VECTOR ALIGNED
      DO 10 J=1,L
        W1=W(1,J)
        W2=W(2,J)
        C0=A(J,1)
        C1=A(J,2)
        C2=A(J,3)
        D0=C1+C2
        D1=C0-C32*D0
        D2=(0.0D0,-1.0D0)*C31*(C1-C2)
        B(1,J)=C0+D0
        B(2,J)=W1*(D1+D2)
        B(3,J)=W2*(D1-D2)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FFT3B(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(M,L,*),B(M,3,*),W(2,*)
      COMPLEX*16 C0,C1,C2,D0,D1,D2,W1,W2
      DATA C31/0.86602540378443865D0/C32/0.5D0/
!
!DIR$ VECTOR ALIGNED
      DO 10 I=1,M
        C0=A(I,1,1)
        C1=A(I,1,2)
        C2=A(I,1,3)
        D0=C1+C2
        D1=C0-C32*D0
        D2=(0.0D0,-1.0D0)*C31*(C1-C2)
        B(I,1,1)=C0+D0
        B(I,2,1)=D1+D2
        B(I,3,1)=D1-D2
   10 CONTINUE
      DO 30 J=2,L
        W1=W(1,J)
        W2=W(2,J)
!DIR$ VECTOR ALIGNED
        DO 20 I=1,M
          C0=A(I,J,1)
          C1=A(I,J,2)
          C2=A(I,J,3)
          D0=C1+C2
          D1=C0-C32*D0
          D2=(0.0D0,-1.0D0)*C31*(C1-C2)
          B(I,1,J)=C0+D0
          B(I,2,J)=W1*(D1+D2)
          B(I,3,J)=W2*(D1-D2)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
      SUBROUTINE FFT4A(A,B,W,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(L,*),B(4,*),W(3,*)
      COMPLEX*16 C0,C1,C2,C3,D0,D1,D2,D3,W1,W2,W3
!
!DIR$ VECTOR ALIGNED
      DO 10 J=1,L
        W1=W(1,J)
        W2=W(2,J)
        W3=W(3,J)
        C0=A(J,1)
        C1=A(J,2)
        C2=A(J,3)
        C3=A(J,4)
        D0=C0+C2
        D1=C0-C2
        D2=C1+C3
        D3=(0.0D0,-1.0D0)*(C1-C3)
        B(1,J)=D0+D2
        B(2,J)=W1*(D1+D3)
        B(3,J)=W2*(D0-D2)
        B(4,J)=W3*(D1-D3)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FFT4B(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(M,L,*),B(M,4,*),W(3,*)
      COMPLEX*16 C0,C1,C2,C3,D0,D1,D2,D3,W1,W2,W3
!
!DIR$ VECTOR ALIGNED
      DO 10 I=1,M
        C0=A(I,1,1)
        C1=A(I,1,2)
        C2=A(I,1,3)
        C3=A(I,1,4)
        D0=C0+C2
        D1=C0-C2
        D2=C1+C3
        D3=(0.0D0,-1.0D0)*(C1-C3)
        B(I,1,1)=D0+D2
        B(I,2,1)=D1+D3
        B(I,3,1)=D0-D2
        B(I,4,1)=D1-D3
   10 CONTINUE
      DO 30 J=2,L
        W1=W(1,J)
        W2=W(2,J)
        W3=W(3,J)
!DIR$ VECTOR ALIGNED
        DO 20 I=1,M
          C0=A(I,J,1)
          C1=A(I,J,2)
          C2=A(I,J,3)
          C3=A(I,J,4)
          D0=C0+C2
          D1=C0-C2
          D2=C1+C3
          D3=(0.0D0,-1.0D0)*(C1-C3)
          B(I,1,J)=D0+D2
          B(I,2,J)=W1*(D1+D3)
          B(I,3,J)=W2*(D0-D2)
          B(I,4,J)=W3*(D1-D3)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
      SUBROUTINE FFT5A(A,B,W,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(L,*),B(5,*),W(4,*)
      COMPLEX*16 C0,C1,C2,C3,C4,D0,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10
      COMPLEX*16 W1,W2,W3,W4
      DATA C51/0.95105651629515357D0/C52/0.61803398874989485D0/   &
           C53/0.55901699437494742D0/C54/0.25D0/
!
!DIR$ VECTOR ALIGNED
      DO 10 J=1,L
        W1=W(1,J)
        W2=W(2,J)
        W3=W(3,J)
        W4=W(4,J)
        C0=A(J,1)
        C1=A(J,2)
        C2=A(J,3)
        C3=A(J,4)
        C4=A(J,5)
        D0=C1+C4
        D1=C2+C3
        D2=C51*(C1-C4)
        D3=C51*(C2-C3)
        D4=D0+D1
        D5=C53*(D0-D1)
        D6=C0-C54*D4
        D7=D6+D5
        D8=D6-D5
        D9=(0.0D0,-1.0D0)*(D2+C52*D3)
        D10=(0.0D0,-1.0D0)*(C52*D2-D3)
        B(1,J)=C0+D4
        B(2,J)=W1*(D7+D9)
        B(3,J)=W2*(D8+D10)
        B(4,J)=W3*(D8-D10)
        B(5,J)=W4*(D7-D9)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FFT5B(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(M,L,*),B(M,5,*),W(4,*)
      COMPLEX*16 C0,C1,C2,C3,C4,D0,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10
      COMPLEX*16 W1,W2,W3,W4
      DATA C51/0.95105651629515357D0/C52/0.61803398874989485D0/    &
           C53/0.55901699437494742D0/C54/0.25D0/
!
!DIR$ VECTOR ALIGNED
      DO 10 I=1,M
        C0=A(I,1,1)
        C1=A(I,1,2)
        C2=A(I,1,3)
        C3=A(I,1,4)
        C4=A(I,1,5)
        D0=C1+C4
        D1=C2+C3
        D2=C51*(C1-C4)
        D3=C51*(C2-C3)
        D4=D0+D1
        D5=C53*(D0-D1)
        D6=C0-C54*D4
        D7=D6+D5
        D8=D6-D5
        D9=(0.0D0,-1.0D0)*(D2+C52*D3)
        D10=(0.0D0,-1.0D0)*(C52*D2-D3)
        B(I,1,1)=C0+D4
        B(I,2,1)=D7+D9
        B(I,3,1)=D8+D10
        B(I,4,1)=D8-D10
        B(I,5,1)=D7-D9
   10 CONTINUE
      DO 30 J=2,L
        W1=W(1,J)
        W2=W(2,J)
        W3=W(3,J)
        W4=W(4,J)
!DIR$ VECTOR ALIGNED
        DO 20 I=1,M
          C0=A(I,J,1)
          C1=A(I,J,2)
          C2=A(I,J,3)
          C3=A(I,J,4)
          C4=A(I,J,5)
          D0=C1+C4
          D1=C2+C3
          D2=C51*(C1-C4)
          D3=C51*(C2-C3)
          D4=D0+D1
          D5=C53*(D0-D1)
          D6=C0-C54*D4
          D7=D6+D5
          D8=D6-D5
          D9=(0.0D0,-1.0D0)*(D2+C52*D3)
          D10=(0.0D0,-1.0D0)*(C52*D2-D3)
          B(I,1,J)=C0+D4
          B(I,2,J)=W1*(D7+D9)
          B(I,3,J)=W2*(D8+D10)
          B(I,4,J)=W3*(D8-D10)
          B(I,5,J)=W4*(D7-D9)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
      SUBROUTINE FFT8A(A,B,W,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(L,*),B(8,*),W(7,*)
      COMPLEX*16 C0,C1,C2,C3,C4,C5,C6,C7,D0,D1,D2,D3,D4,D5,D6,D7
      COMPLEX*16 E0,E1,E2,E3,E4,E5,E6,E7,E8,E9,W1,W2,W3,W4,W5,W6,W7
      DATA C81/0.70710678118654752D0/
!
!DIR$ VECTOR ALIGNED
      DO 10 J=1,L
        W1=W(1,J)
        W2=W(2,J)
        W3=W(3,J)
        W4=W(4,J)
        W5=W(5,J)
        W6=W(6,J)
        W7=W(7,J)
        C0=A(J,1)
        C1=A(J,2)
        C2=A(J,3)
        C3=A(J,4)
        C4=A(J,5)
        C5=A(J,6)
        C6=A(J,7)
        C7=A(J,8)
        D0=C0+C4
        D1=C0-C4
        D2=C2+C6
        D3=(0.0D0,-1.0D0)*(C2-C6)
        D4=C1+C5
        D5=C1-C5
        D6=C3+C7
        D7=C3-C7
        E0=D0+D2
        E1=D0-D2
        E2=D4+D6
        E3=(0.0D0,-1.0D0)*(D4-D6)
        E4=C81*(D5-D7)
        E5=(0.0D0,-1.0D0)*C81*(D5+D7)
        E6=D1+E4
        E7=D1-E4
        E8=D3+E5
        E9=D3-E5
        B(1,J)=E0+E2
        B(2,J)=W1*(E6+E8)
        B(3,J)=W2*(E1+E3)
        B(4,J)=W3*(E7-E9)
        B(5,J)=W4*(E0-E2)
        B(6,J)=W5*(E7+E9)
        B(7,J)=W6*(E1-E3)
        B(8,J)=W7*(E6-E8)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FFT8B(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(M,L,*),B(M,8,*),W(7,*)
      COMPLEX*16 C0,C1,C2,C3,C4,C5,C6,C7,D0,D1,D2,D3,D4,D5,D6,D7
      COMPLEX*16 E0,E1,E2,E3,E4,E5,E6,E7,E8,E9,W1,W2,W3,W4,W5,W6,W7
      DATA C81/0.70710678118654752D0/
!
!DIR$ VECTOR ALIGNED
      DO 10 I=1,M
        C0=A(I,1,1)
        C1=A(I,1,2)
        C2=A(I,1,3)
        C3=A(I,1,4)
        C4=A(I,1,5)
        C5=A(I,1,6)
        C6=A(I,1,7)
        C7=A(I,1,8)
        D0=C0+C4
        D1=C0-C4
        D2=C2+C6
        D3=(0.0D0,-1.0D0)*(C2-C6)
        D4=C1+C5
        D5=C1-C5
        D6=C3+C7
        D7=C3-C7
        E0=D0+D2
        E1=D0-D2
        E2=D4+D6
        E3=(0.0D0,-1.0D0)*(D4-D6)
        E4=C81*(D5-D7)
        E5=(0.0D0,-1.0D0)*C81*(D5+D7)
        E6=D1+E4
        E7=D1-E4
        E8=D3+E5
        E9=D3-E5
        B(I,1,1)=E0+E2
        B(I,2,1)=E6+E8
        B(I,3,1)=E1+E3
        B(I,4,1)=E7-E9
        B(I,5,1)=E0-E2
        B(I,6,1)=E7+E9
        B(I,7,1)=E1-E3
        B(I,8,1)=E6-E8
   10 CONTINUE
      DO 30 J=2,L
        W1=W(1,J)
        W2=W(2,J)
        W3=W(3,J)
        W4=W(4,J)
        W5=W(5,J)
        W6=W(6,J)
        W7=W(7,J)
!DIR$ VECTOR ALIGNED
        DO 20 I=1,M
          C0=A(I,J,1)
          C1=A(I,J,2)
          C2=A(I,J,3)
          C3=A(I,J,4)
          C4=A(I,J,5)
          C5=A(I,J,6)
          C6=A(I,J,7)
          C7=A(I,J,8)
          D0=C0+C4
          D1=C0-C4
          D2=C2+C6
          D3=(0.0D0,-1.0D0)*(C2-C6)
          D4=C1+C5
          D5=C1-C5
          D6=C3+C7
          D7=C3-C7
          E0=D0+D2
          E1=D0-D2
          E2=D4+D6
          E3=(0.0D0,-1.0D0)*(D4-D6)
          E4=C81*(D5-D7)
          E5=(0.0D0,-1.0D0)*C81*(D5+D7)
          E6=D1+E4
          E7=D1-E4
          E8=D3+E5
          E9=D3-E5
          B(I,1,J)=E0+E2
          B(I,2,J)=W1*(E6+E8)
          B(I,3,J)=W2*(E1+E3)
          B(I,4,J)=W3*(E7-E9)
          B(I,5,J)=W4*(E0-E2)
          B(I,6,J)=W5*(E7+E9)
          B(I,7,J)=W6*(E1-E3)
          B(I,8,J)=W7*(E6-E8)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END

!
!     FFTE: A FAST FOURIER TRANSFORM PACKAGE
!
!     (C) COPYRIGHT SOFTWARE, 2000-2004, 2008-2014, ALL RIGHTS RESERVED
!                BY
!         DAISUKE TAKAHASHI
!         FACULTY OF ENGINEERING, INFORMATION AND SYSTEMS
!         UNIVERSITY OF TSUKUBA
!         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
!         E-MAIL: daisuke@cs.tsukuba.ac.jp
!
!
!     RADIX-2, 3, 4, 5 AND 8 MULTIPLE FFT ROUTINE
!
!     FORTRAN77 SOURCE PROGRAM
!
!     WRITTEN BY DAISUKE TAKAHASHI
!
      SUBROUTINE MFFT235A(A,B,W,NS,N,IP)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*),W(*)
      DIMENSION IP(*)
!
      IF (IP(1) .NE. 1) THEN
        KP4=2-MOD(IP(1)+2,3)
        KP8=(IP(1)-KP4)/3
      ELSE
        KP4=0
        KP8=0
      END IF
!
      KEY=1
      J=1
      L=N
      M=1
      DO 10 K=1,KP8
        L=L/8
        IF (L .GE. 2) THEN
          IF (KEY .GE. 0) THEN
            CALL FFT8B(A,B,W(J),NS*M,L)
          ELSE
            CALL FFT8B(B,A,W(J),NS*M,L)
          END IF
          KEY=-KEY
        ELSE
          IF (KEY .GE. 0) THEN
            CALL FFT8B(A,A,W(J),NS*M,L)
          ELSE
            CALL FFT8B(B,A,W(J),NS*M,L)
          END IF
        END IF
        M=M*8
        J=J+L*7
   10 CONTINUE
      DO 20 K=1,IP(3)
        L=L/5
        IF (L .GE. 2) THEN
          IF (KEY .GE. 0) THEN
            CALL FFT5B(A,B,W(J),NS*M,L)
          ELSE
            CALL FFT5B(B,A,W(J),NS*M,L)
          END IF
          KEY=-KEY
        ELSE
          IF (KEY .GE. 0) THEN
            CALL FFT5B(A,A,W(J),NS*M,L)
          ELSE
            CALL FFT5B(B,A,W(J),NS*M,L)
          END IF
        END IF
        M=M*5
        J=J+L*4
   20 CONTINUE
      DO 30 K=1,KP4
        L=L/4
        IF (L .GE. 2) THEN
          IF (KEY .GE. 0) THEN
            CALL FFT4B(A,B,W(J),NS*M,L)
          ELSE
            CALL FFT4B(B,A,W(J),NS*M,L)
          END IF
          KEY=-KEY
        ELSE
          IF (KEY .GE. 0) THEN
            CALL FFT4B(A,A,W(J),NS*M,L)
          ELSE
            CALL FFT4B(B,A,W(J),NS*M,L)
          END IF
        END IF
        M=M*4
        J=J+L*3
   30 CONTINUE
      DO 40 K=1,IP(2)
        L=L/3
        IF (L .GE. 2) THEN
          IF (KEY .GE. 0) THEN
            CALL FFT3B(A,B,W(J),NS*M,L)
          ELSE
            CALL FFT3B(B,A,W(J),NS*M,L)
          END IF
          KEY=-KEY
        ELSE
          IF (KEY .GE. 0) THEN
            CALL FFT3B(A,A,W(J),NS*M,L)
          ELSE
            CALL FFT3B(B,A,W(J),NS*M,L)
          END IF
        END IF
        M=M*3
        J=J+L*2
   40 CONTINUE
      IF (IP(1) .EQ. 1) THEN
        IF (KEY .GE. 0) THEN
          CALL FFT2(A,A,NS*M)
        ELSE
          CALL FFT2(B,A,NS*M)
        END IF
      END IF
      RETURN
      END
      SUBROUTINE MFFT235B(A,B,W,NS,N,IP)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*),W(*)
      DIMENSION IP(*)
!
      IF (IP(1) .NE. 1) THEN
        KP4=2-MOD(IP(1)+2,3)
        KP8=(IP(1)-KP4)/3
      ELSE
        KP4=0
        KP8=0
      END IF
!
      KEY=1
      J=1
      L=N
      M=1
      DO 10 K=1,KP8
        L=L/8
        IF (L .GE. 2) THEN
          IF (KEY .GE. 0) THEN
            CALL FFT8(A,B,W(J),NS*M,L)
          ELSE
            CALL FFT8(B,A,W(J),NS*M,L)
          END IF
          KEY=-KEY
        ELSE
          IF (KEY .GE. 0) THEN
            CALL FFT8(A,B,W(J),NS*M,L)
          ELSE
            CALL FFT8(B,B,W(J),NS*M,L)
          END IF
        END IF
        M=M*8
        J=J+L*7
   10 CONTINUE
      DO 20 K=1,IP(3)
        L=L/5
        IF (L .GE. 2) THEN
          IF (KEY .GE. 0) THEN
            CALL FFT5(A,B,W(J),NS*M,L)
          ELSE
            CALL FFT5(B,A,W(J),NS*M,L)
          END IF
          KEY=-KEY
        ELSE
          IF (KEY .GE. 0) THEN
            CALL FFT5(A,B,W(J),NS*M,L)
          ELSE
            CALL FFT5(B,B,W(J),NS*M,L)
          END IF
        END IF
        M=M*5
        J=J+L*4
   20 CONTINUE
      DO 30 K=1,KP4
        L=L/4
        IF (L .GE. 2) THEN
          IF (KEY .GE. 0) THEN
            CALL FFT4(A,B,W(J),NS*M,L)
          ELSE
            CALL FFT4(B,A,W(J),NS*M,L)
          END IF
          KEY=-KEY
        ELSE
          IF (KEY .GE. 0) THEN
            CALL FFT4(A,B,W(J),NS*M,L)
          ELSE
            CALL FFT4(B,B,W(J),NS*M,L)
          END IF
        END IF
        M=M*4
        J=J+L*3
   30 CONTINUE
      DO 40 K=1,IP(2)
        L=L/3
        IF (L .GE. 2) THEN
          IF (KEY .GE. 0) THEN
            CALL FFT3(A,B,W(J),NS*M,L)
          ELSE
            CALL FFT3(B,A,W(J),NS*M,L)
          END IF
          KEY=-KEY
        ELSE
          IF (KEY .GE. 0) THEN
            CALL FFT3(A,B,W(J),NS*M,L)
          ELSE
            CALL FFT3(B,B,W(J),NS*M,L)
          END IF
        END IF
        M=M*3
        J=J+L*2
   40 CONTINUE
      IF (IP(1) .EQ. 1) THEN
        IF (KEY .GE. 0) THEN
          CALL FFT2(A,B,NS*M)
        ELSE
          CALL FFT2(B,B,NS*M)
        END IF
      END IF
      RETURN
      END
      SUBROUTINE ZTRANS(A,B,NX,NY)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*)
      DIMENSION LNX(3),LNY(3)
!
      CALL FACTOR(NX,LNX)
      CALL FACTOR(NY,LNY)
!
      IF (NX .EQ. 1 .OR. NY .EQ. 1) THEN
!DIR$ VECTOR ALIGNED
        DO 10 I=1,NX*NY
          B(I)=A(I)
   10   CONTINUE
        RETURN
      END IF
!
      IF (LNX(1)+LNY(1) .LE. 1) THEN
        CALL ZTRANSA(A,B,NX,NY)
      ELSE
        CALL ZTRANSB(A,B,NX,NY)
      END IF
      RETURN
      END
      SUBROUTINE ZTRANSA(A,B,NX,NY)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(NX,*),B(NY,*)
!
      DO 20 I=1,NX
!DIR$ VECTOR ALIGNED
        DO 10 J=1,NY
          B(J,I)=A(I,J)
   10   CONTINUE
   20 CONTINUE
      RETURN
      END
      SUBROUTINE ZTRANSB(A,B,NX,NY)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(NX,*),B(NY,*)
!
      IF (NY .GE. NX) THEN
        DO 20 I=0,NX-1
!DIR$ VECTOR ALIGNED
          DO 10 J=1,NX-I
            B(J,I+J)=A(I+J,J)
   10     CONTINUE
   20   CONTINUE
        DO 40 I=1,NY-NX
!DIR$ VECTOR ALIGNED
          DO 30 J=1,NX
            B(I+J,J)=A(J,I+J)
   30     CONTINUE
   40   CONTINUE
        DO 60 I=NY-NX+1,NY-1
!DIR$ VECTOR ALIGNED
          DO 50 J=1,NY-I
            B(I+J,J)=A(J,I+J)
   50     CONTINUE
   60   CONTINUE
      ELSE
        DO 80 I=0,NY-1
!DIR$ VECTOR ALIGNED
          DO 70 J=1,NY-I
            B(I+J,J)=A(J,I+J)
   70     CONTINUE
   80   CONTINUE
        DO 100 I=1,NX-NY
!DIR$ VECTOR ALIGNED
          DO 90 J=1,NY
            B(J,I+J)=A(I+J,J)
   90     CONTINUE
  100   CONTINUE
        DO 120 I=NX-NY+1,NX-1
!DIR$ VECTOR ALIGNED
          DO 110 J=1,NX-I
            B(J,I+J)=A(I+J,J)
  110     CONTINUE
  120   CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE ZTRANS2(A,B,NX,NY,NZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(NX,NY,*),B(NY,NX,*)
!
      DO 10 K=1,NZ
        CALL ZTRANS(A(1,1,K),B(1,1,K),NX,NY)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE MZTRANS(A,B,NS,NX,NY)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(NS,NX,*),B(NS,NY,*)
!
      IF (NS .EQ. 1) THEN
        CALL ZTRANS(A(1,1,1),B(1,1,1),NX,NY)
      ELSE
        DO 30 I=1,NX
          DO 20 J=1,NY
!DIR$ VECTOR ALIGNED
            DO 10 K=1,NS
              B(K,J,I)=A(K,I,J)
   10       CONTINUE
   20     CONTINUE
   30   CONTINUE
      END IF
      RETURN
      END
