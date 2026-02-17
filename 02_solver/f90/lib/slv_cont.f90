!=======================================================================
      SUBROUTINE DIVGS(DIVGSUM)
!=======================================================================
!
!     RHS OF POISSON EQUATION
!
!     Options
!       MASSON = 0, IMOVINGON = 0: w/o IBM body
!       MASSON = 1, IMOVINGON = 0: w/  stationary body
!       MASSON = 1, IMOVINGON = 1: w/  moving body
!
!     Reference for IBM
!     Kim, Kim & Choi, 2001, An immersed-boundary finite-volume method
!       for simulations of flow in complex geometries, J. Comp. Phys.,
!       171, 132-150.
!
!     Variables
!       DIVGSUM: divergence of velocity
!       QMASS: mass source/sink for IBM body
!       UBD,VBD,WBD: translational velocity of body (IMOVINGON=1)
!                    These velocities are defined in lica_movingIBM.f90.
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,QMASS,IFC,JFC,KFC
      IMPLICIT NONE
      INTEGER*8    :: I,J,K
      INTEGER*8    :: II,JJ,KK,IMM,JMM,KMM,N
      REAL*8       :: DIVG1,DIVG2,DIVG3
      REAL*8       :: DIVGSUM(0:N1,0:N2,0:N3)
      REAL*8       :: USTMP,VSTMP,WSTMP
      REAL*8       :: UBD,VBD,WBD

!############################################################### W/O IBM
!$OMP PARALLEL DO  &
!$OMP private(DIVG1,DIVG2,DIVG3)
      DO 73 K=1,N3M
      DO 73 J=1,N2M
      DO 73 I=1,N1M
        DIVG1=(U(I+1,J,K)-U(I,J,K))*F2FXI(I)
        DIVG2=(V(I,J+1,K)-V(I,J,K))*F2FYI(J)
        DIVG3=(W(I,J,K+1)-W(I,J,K))*F2FZI(K)
        DIVGSUM(I,J,K)=(DIVG1+DIVG2+DIVG3)
   73 CONTINUE

!############################################### W/  STATIONALY IBM body
      IF (MASSON.EQ.1 .AND. IMOVINGON.EQ.0) THEN

      QMASS = 0.

!$OMP PARALLEL DO  &
!$OMP private(II,JJ,KK)
      DO N=1,NBODY(1)
        II=IFC(N,1)
        JJ=JFC(N,1)
        KK=KFC(N,1)
        QMASS(II,JJ,KK) = QMASS(II,JJ,KK)-U(II,JJ,KK)*F2FXI(II)
        DIVGSUM(II,JJ,KK) = DIVGSUM(II,JJ,KK)+U(II,JJ,KK)*F2FXI(II)
      ENDDO
!$OMP PARALLEL DO  &
!$OMP private(II,JJ,KK,IMM)
      DO N=1,NBODY(1)
        II=IFC(N,1)
        JJ=JFC(N,1)
        KK=KFC(N,1)
        IMM=IMV(II)
        QMASS(IMM,JJ,KK)=QMASS(IMM,JJ,KK)+U(II,JJ,KK)*F2FXI(IMM)
        DIVGSUM(IMM,JJ,KK)=DIVGSUM(IMM,JJ,KK)-U(II,JJ,KK)*F2FXI(IMM)
      ENDDO

!$OMP PARALLEL DO  &
!$OMP private(II,JJ,KK)
      DO N=1,NBODY(2)
        II=IFC(N,2)
        JJ=JFC(N,2)
        KK=KFC(N,2)
        QMASS(II,JJ,KK)=QMASS(II,JJ,KK)-V(II,JJ,KK)*F2FYI(JJ)
        DIVGSUM(II,JJ,KK)=DIVGSUM(II,JJ,KK)+V(II,JJ,KK)*F2FYI(JJ)
      ENDDO
!$OMP PARALLEL DO  &
!$OMP private(II,JJ,KK,JMM)
      DO N=1,NBODY(2)
        II=IFC(N,2)
        JJ=JFC(N,2)
        KK=KFC(N,2)
        JMM=JMV(JJ)
        QMASS(II,JMM,KK)=QMASS(II,JMM,KK)+V(II,JJ,KK)*F2FYI(JMM)
        DIVGSUM(II,JMM,KK)=DIVGSUM(II,JMM,KK)-V(II,JJ,KK)*F2FYI(JMM)
      ENDDO

!$OMP PARALLEL DO  &
!$OMP private(II,JJ,KK)
      DO N=1,NBODY(3)
        II=IFC(N,3)
        JJ=JFC(N,3)
        KK=KFC(N,3)
        QMASS(II,JJ,KK)=QMASS(II,JJ,KK)-W(II,JJ,KK)*F2FZI(KK)
        DIVGSUM(II,JJ,KK)=DIVGSUM(II,JJ,KK)+W(II,JJ,KK)*F2FZI(KK)
      ENDDO
!$OMP PARALLEL DO  &
!$OMP private(II,JJ,KK,KMM)
      DO N=1,NBODY(3)
        II=IFC(N,3)
        JJ=JFC(N,3)
        KK=KFC(N,3)
        KMM=KMV(KK)
        QMASS(II,JJ,KMM)=QMASS(II,JJ,KMM)+W(II,JJ,KK)*F2FZI(KMM)
        DIVGSUM(II,JJ,KMM)=DIVGSUM(II,JJ,KMM)-W(II,JJ,KK)*F2FZI(KMM)
      ENDDO

! !################################################### W/  MOVING IBM body
!       ELSEIF (MASSON.EQ.1 .AND. IMOVINGON.EQ.1) THEN

!       QMASS = 0.

! !$OMP PARALLEL DO  &
! !$OMP private(II,JJ,KK,USTMP)
!       DO N=1,NBODY(1)
!         II=IFC(N,1)
!         JJ=JFC(N,1)
!         KK=KFC(N,1)

!         USTMP=U(II,JJ,KK)-UBD(X(II),YMP(JJ),ZMP(KK))

!         QMASS (II,JJ,KK) = QMASS  (II,JJ,KK)-USTMP*F2FXI(II)
!         DIVGSUM(II,JJ,KK)= DIVGSUM(II,JJ,KK)+USTMP*F2FXI(II)
!       ENDDO
! !$OMP PARALLEL DO  &
! !$OMP private(II,JJ,KK,IMM,USTMP)
!       DO N=1,NBODY(1)
!         II=IFC(N,1)
!         JJ=JFC(N,1)
!         KK=KFC(N,1)
!         IMM=IMV(II)

!         USTMP=U(II,JJ,KK)-UBD(X(II),YMP(JJ),ZMP(KK))

!         QMASS  (IMM,JJ,KK)=QMASS  (IMM,JJ,KK)+USTMP*F2FXI(IMM)
!         DIVGSUM(IMM,JJ,KK)=DIVGSUM(IMM,JJ,KK)-USTMP*F2FXI(IMM)
!       ENDDO

! !$OMP PARALLEL DO  &
! !$OMP private(II,JJ,KK,VSTMP)
!       DO N=1,NBODY(2)
!         II=IFC(N,2)
!         JJ=JFC(N,2)
!         KK=KFC(N,2)

!         VSTMP=V(II,JJ,KK)-VBD(XMP(II),Y(JJ),ZMP(KK))

!         QMASS  (II,JJ,KK)=QMASS  (II,JJ,KK)-VSTMP*F2FYI(JJ)
!         DIVGSUM(II,JJ,KK)=DIVGSUM(II,JJ,KK)+VSTMP*F2FYI(JJ)
!       ENDDO
! !$OMP PARALLEL DO  &
! !$OMP private(II,JJ,KK,JMM,VSTMP)
!       DO N=1,NBODY(2)
!         II=IFC(N,2)
!         JJ=JFC(N,2)
!         KK=KFC(N,2)
!         JMM=JMV(JJ)

!         VSTMP=V(II,JJ,KK)-VBD(XMP(II),Y(JJ),ZMP(KK))

!         QMASS  (II,JMM,KK)=QMASS  (II,JMM,KK)+VSTMP*F2FYI(JMM)
!         DIVGSUM(II,JMM,KK)=DIVGSUM(II,JMM,KK)-VSTMP*F2FYI(JMM)
!       ENDDO

! !$OMP PARALLEL DO  &
! !$OMP private(II,JJ,KK,WSTMP)
!       DO N=1,NBODY(3)
!         II=IFC(N,3)
!         JJ=JFC(N,3)
!         KK=KFC(N,3)

!         WSTMP=W(II,JJ,KK)-WBD(XMP(II),YMP(JJ),Z(KK))

!         QMASS  (II,JJ,KK)=QMASS  (II,JJ,KK)-WSTMP*F2FZI(KK)
!         DIVGSUM(II,JJ,KK)=DIVGSUM(II,JJ,KK)+WSTMP*F2FZI(KK)
!       ENDDO
! !$OMP PARALLEL DO  &
! !$OMP private(II,JJ,KK,KMM,WSTMP)
!       DO N=1,NBODY(3)
!         II=IFC(N,3)
!         JJ=JFC(N,3)
!         KK=KFC(N,3)
!         KMM=KMV(KK)

!         WSTMP=W(II,JJ,KK)-WBD(XMP(II),YMP(JJ),Z(KK))

!         QMASS  (II,JJ,KMM)=QMASS  (II,JJ,KMM)+WSTMP*F2FZI(KMM)
!         DIVGSUM(II,JJ,KMM)=DIVGSUM(II,JJ,KMM)-WSTMP*F2FZI(KMM)
!       ENDDO

      ENDIF

!$OMP PARALLEL DO
      DO 80 K=1,N3M
      DO 80 J=1,N2M
      DO 80 I=1,N1M
      DIVGSUM(I,J,K)=DIVGSUM(I,J,K)*DTCONSTI
   80 CONTINUE

      RETURN
      END
!=======================================================================