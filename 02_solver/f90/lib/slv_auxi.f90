!=======================================================================
       SUBROUTINE STEPINIT
!=======================================================================
      USE MOD_COMMON
      USE MOD_FLOWARRAY
      IMPLICIT NONE
      REAL*8      ::  DTCFL
      INTEGER*8   ::  ICFL,JCFL,KCFL

      NTIME = NTIME + 1
      M = NTIME
      FCVAVG = 0.
      IF (ICH .EQ. 1) PMIAVG = 0.

      CALL CFL(CFLMAX,ICFL,JCFL,KCFL)         ! CALCULATE CFL NUMBER
      
      IF (IDTOPT.NE.0 .AND. CFLMAX.NE.0.) THEN
         DTCFL = DMIN1(DT*CFLFAC/CFLMAX, DT*(0.80+0.20*CFLFAC/CFLMAX))
         IF (IDTOPT.EQ.1) DT = DTCFL
      ENDIF
      
      IF (CFLMAX.GT.(CFLFAC*1.1)) THEN
         PRINT*,' '
         WRITE(*,310) NTIME,CFLMAX,TIME
      ELSE
         PRINT*,' '
         WRITE(*,320) NTIME,TIME,DT
      ENDIF
 310  FORMAT(I15,'   CFL NUMBER EXCEED GIVEN CFL LIMIT :',ES18.5,' AT ',F12.5)
 320  FORMAT('--------------------------',I6,'  TIME=',F10.5,'  DT=',F12.8)

      RETURN
      END SUBROUTINE STEPINIT
!=======================================================================
!=======================================================================
       SUBROUTINE CFL(CFLM,ICFL,JCFL,KCFL)
!=======================================================================
!
!     Calculate the maximum CFL number of flow field
!     (http://en.wikipedia.org/wiki/Courant-Friedrichs-Lewy_condition)
!
!     CFL#=U_i*DT/DX_i
!     CFLI=(Uc/DX+Vc/DY+Wc/DZ)*DT  (c for cell center)
!     CFLMPT: index where the CFL# is maximum
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W
      IMPLICIT NONE
      REAL*8       :: CFLM
      REAL*8       :: CFLI(0:3)
      INTEGER*8    :: ICFL,JCFL,KCFL
      INTEGER*8    :: I,J,K

      CFLM = 0.

!$OMP PARALLEL DO private(CFLI) reduction(MAX:CFLM)
      DO K=1,N3M
        DO J=1,N2M
          DO I=1,N1M
            CFLI(1)=ABS(U(I,J,K)+U(I+1,J,K))*0.5*F2FXI(I)
            CFLI(2)=ABS(V(I,J,K)+V(I,J+1,K))*0.5*F2FYI(J)
            CFLI(3)=ABS(W(I,J,K)+W(I,J,K+1))*0.5*F2FZI(K)
            CFLI(0)=(CFLI(1)+CFLI(2)+CFLI(3))*DT
            IF (CFLI(0) .GE. CFLM) THEN
              CFLM=CFLI(0)
              ICFL=I
              JCFL=J
              KCFL=K
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CFLMAX = CFLM
      
      IF(NTIME .NE. 1) THEN
        WRITE(*,149) CFLMAX,XMP(ICFL),YMP(JCFL),ZMP(KCFL),ICFL,JCFL,KCFL
 149    FORMAT('CFLMAX =  ',F10.7,' @ ',3F10.4,' , ',3I5)
      ENDIF

      RETURN
      END SUBROUTINE CFL
!=======================================================================
!=======================================================================
       SUBROUTINE SUBSTEPINIT(SUBSTEP)
!=======================================================================
      USE MOD_COMMON
      USE MOD_FLOWARRAY
      IMPLICIT NONE
      INTEGER*8    :: SUBSTEP        ! RK3 SUBSTEP = 1, 2, 3

      ALPHA    = 0.5*(GAMMA(SUBSTEP)+RO(SUBSTEP))
      DTCONST  = DT *(GAMMA(SUBSTEP)+RO(SUBSTEP))
      DTCONSTI = 1./DTCONST
      TEST1    = RESID1*DTCONSTI              ! POISS. CONVG. CRITERION
      ACOEF    = ALPHA*DT/RE
      ACOEFI   = 1./ACOEF
      PMIAVG   = PMIAVG + 2.*ALPHA*PMI(0)
      SUBDT    = SUBDT + DTCONST               ! SUBTIME FOR RK3 METHOD
      MSUB     = SUBSTEP

      RETURN
      END SUBROUTINE SUBSTEPINIT
!=======================================================================
!=======================================================================
       SUBROUTINE QVOLCALC(QQ)
!=======================================================================
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U
      IMPLICIT NONE
      REAL*8       :: QQ,FUNCBODY
      INTEGER*8    :: I,J,K

      QQ = 0.

!$OMP PARALLEL DO reduction(+:QQ)
      DO K=1,N3M
        DO J=1,N2M
          DO I=1,N1M
            IF (FUNCBODY(X(I),YMP(J),ZMP(K)).GE.1.E-10) THEN 
              QQ = QQ + U(I,J,K)*C2CX(I)*F2FY(J)*F2FZ(K)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE QVOLCALC
!=======================================================================
!=======================================================================
       SUBROUTINE TVOLCALC(QQ)
!=======================================================================
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : T,CSTAR
      IMPLICIT NONE
      REAL*8       :: QQ,FUNCBODY
      INTEGER*8    :: I,J,K

      QQ = 0.

!$OMP PARALLEL DO reduction(+:QQ)
      DO K=1,N3M
        DO J=1,N2M
          DO I=1,N1M
            QQ = QQ + T(I,J,K) / CSTAR(I,J,K) * F2FX(I)*F2FY(J)*F2FZ(K)
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE TVOLCALC
!=======================================================================
!=======================================================================
       SUBROUTINE QVOLCORR
!=======================================================================
      USE MOD_COMMON
      USE MOD_FLOWARRAY
      IMPLICIT NONE
      REAL*8                   :: FUNCBODY
      REAL*8                   :: FLOWVOL
      INTEGER*8                :: I,J,K

      FLOWVOL = 0.

!$OMP PARALLEL DO reduction(+:FLOWVOL)
      DO K=1,N3M
        DO J=1,N2M
          DO I=1,N1M
            IF (FUNCBODY(X(I),YMP(J),ZMP(K)).GE.1.E-10) THEN 
              FLOWVOL = FLOWVOL + C2CX(I)*F2FY(J)*F2FZ(K)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      PHCAP = (QVOL(2) - QVOL(1)) / FLOWVOL

      RETURN
      END SUBROUTINE QVOLCORR
!=======================================================================
!=======================================================================
       SUBROUTINE TVOLCORR
!=======================================================================
      USE MOD_COMMON
      USE MOD_FLOWARRAY
      IMPLICIT NONE
      REAL*8                   :: CVOL
      INTEGER*8                :: I,J,K

      CVOL = 0.

!$OMP PARALLEL DO reduction(+:CVOL)
      DO K=1,N3M
        DO J=1,N2M
          DO I=1,N1M
            CVOL = CVOL + 1. / CSTAR(I,J,K) * F2FX(I)*F2FY(J)*F2FZ(K)
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      THCAP = (TVOL(2) - TVOL(1)) / CVOL 

      RETURN
      END SUBROUTINE TVOLCORR
!=======================================================================
!=======================================================================
      SUBROUTINE UCALC(PHI)
!=======================================================================
!
!     Calculate velocity (u_i) from u_i hat
!     u_i hat is derived from pseudo-pressure, phi.
!
!     Definition of pseudo-pressure, phi, is as follows:
!     ({u_i}^k - u_i hat)/(2.*alpha_k*dt) = - d({phi}^k)/d(x_i)
!
!     From the definition of phi, following equation is derived:
!     {u_i}^k = u_i hat - 2.*a_k*DT*d({phi}^k)/d(x_i)
!
!-----------------------------------------------------------------------
       USE MOD_COMMON
       USE MOD_FLOWARRAY, ONLY : U,V,W,T
       IMPLICIT NONE
       INTEGER*8     :: I,J,K
       REAL*8        :: IDUM,FUNCBODY
       INTEGER*8     :: IM,KM
       REAL*8        :: PHI(0:N1,0:N2,0:N3),FLXCR

!$OMP PARALLEL DO &
!$OMP private(IM)
      DO 21 K=1,N3M
      DO 21 J=1,N2M
      DO 21 I=I_BGPX,N1M                  ! I=1,N1 => BOUNDARY
         IM=IMV(I)
         U(I,J,K)=U(I,J,K)                                   &
                 -DTCONST*(PHI(I,J,K)-PHI(IM,J,K))*C2CXI(I)  
         IF(ICH.EQ.1) THEN
           IF (FUNCBODY(X(I),YMP(J),ZMP(K)).GE.1.E-10) THEN 
              U(I,J,K)=U(I,J,K)-PHCAP
           ENDIF
         ENDIF
   21 CONTINUE

!$OMP PARALLEL DO
      DO 31 K=1,N3M
      DO 31 J=2,N2M                       ! J=1,N2 => BOUNDARY
      DO 31 I=1,N1M
         V(I,J,K)=V(I,J,K)-DTCONST*(PHI(I,J,K)-PHI(I,J-1,K))*C2CYI(J)
   31 CONTINUE

!$OMP PARALLEL DO &
!$OMP private(KM)
      DO 41 K=K_BGPZ,N3M                  ! K=0,N3 => BOUNDARY
         KM=KMV(K)
      DO 41 J=1,N2M
      DO 41 I=1,N1M
         W(I,J,K)=W(I,J,K)-DTCONST*(PHI(I,J,K)-PHI(I,J,KM))*C2CZI(K)
   41 CONTINUE

      IDUM = IHIST
      FLXCR = 0.

! !$OMP PARALLEL DO
!       DO K=1,N3M
!       DO J=1,N2M
!       CALL RANDOM_NUMBER(IDUM)
!       U(1,J,K)=UBULK_I + UBULK_I*EPS_PTR*((IDUM * 2. ) - 1.)/50.
!       FLXCR = FLXCR + (U(1,J,K)-UBULK_I)
!       ENDDO
!       ENDDO

! !$OMP PARALLEL DO
!       DO K=1,N3M
!       DO J=1,N2M
!       U(1,J,K) = U(1,J,K) - FLXCR / (N2M*N3M)
!       ENDDO
!       ENDDO

!       FLXCR = 0.

! !$OMP PARALLEL DO
!       DO K=1,N3M
!       DO J=2,N2M
!       CALL RANDOM_NUMBER(IDUM)
!       V(1,J,K)=0. + EPS_PTR*((IDUM * 2. ) - 1.)/50.
!       FLXCR = FLXCR + (V(1,J,K)-0.)
!       ENDDO
!       ENDDO

! !$OMP PARALLEL DO
!       DO K=1,N3M
!       DO J=2,N2M
!       V(1,J,K) = V(1,J,K) - FLXCR / ((N2M-1)*N3M)
!       ENDDO
!       ENDDO

!     FLXCR = 0.

! !$OMP PARALLEL DO
!       DO K=K_BGPZ,N3M
!       DO J=1,N2M
!       CALL RANDOM_NUMBER(IDUM)
!       W(1,J,K)=0. + EPS_PTR*(IDUM*.3557 - .3557/2.)
!       FLXCR = FLXCR + (W(1,J,K)-0.)
!       ENDDO
!       ENDDO

! !$OMP PARALLEL DO
!       DO K=K_BGPZ,N3M
!       DO J=1,N2M
!       W(1,J,K) = W(1,J,K) - FLXCR / (N2M * (N3M-K_BGPZ))
!       ENDDO
!       ENDDO


!!!!!!     NEUMANN & DIRICHLET BOUNDARY CONDITION
      IF (JUT.EQ.1) THEN
!$OMP PARALLEL DO
         DO K=0,N3
         DO I=1,N1
            U(I,N2,K)=U(I,N2M,K)
         ENDDO
         ENDDO
      ENDIF
      IF (JWT.EQ.1) THEN
!$OMP PARALLEL DO
         DO K=1,N3
         DO I=0,N1
            W(I,N2,K)=W(I,N2M,K)
         ENDDO
         ENDDO
      ENDIF
      IF (JUB.EQ.1) THEN
!$OMP PARALLEL DO
         DO K=0,N3
         DO I=1,N1
            U(I,0,K)=U(I,1,K)
         ENDDO
         ENDDO
      ENDIF
      IF (JWB.EQ.1) THEN
!$OMP PARALLEL DO
         DO K=1,N3
         DO I=0,N1
            W(I,0,K)=W(I,1,K)
         ENDDO
         ENDDO
      ENDIF
      IF (KUT.EQ.1) THEN
!$OMP PARALLEL DO
         DO J=0,N2
         DO I=1,N1
            U(I,J,N3)=U(I,J,N3M)
         ENDDO
         ENDDO
      ENDIF
      IF (KVT.EQ.1) THEN
!$OMP PARALLEL DO
         DO J=1,N2
         DO I=0,N1
            V(I,J,N3)=V(I,J,N3M)
         ENDDO
         ENDDO
      ENDIF
      IF (KUB.EQ.1) THEN
!$OMP PARALLEL DO
         DO J=0,N2
         DO I=1,N1
            U(I,J,0)=U(I,J,1)
         ENDDO
         ENDDO
      ENDIF
      IF (KVB.EQ.1) THEN
!$OMP PARALLEL DO
         DO J=1,N2
         DO I=0,N1
            V(I,J,0)=V(I,J,1)
         ENDDO
         ENDDO
      ENDIF

!     Z PERIODICITY
      IF (ZPRDIC .EQ. 1) THEN
!$OMP PARALLEL DO
      DO J=0,N2
      DO I=1,N1
         U(I,J,0) =U(I,J,N3M)
         U(I,J,N3)=U(I,J,1)
      ENDDO
      ENDDO

!$OMP PARALLEL DO
      DO J=1,N2
      DO I=0,N1
         V(I,J,0) =V(I,J,N3M)
         V(I,J,N3)=V(I,J,1)
      ENDDO
      ENDDO

!$OMP PARALLEL DO
      DO J=0,N2
      DO I=0,N1
         W(I,J,0) =W(I,J,N3M)
         W(I,J,N3)=W(I,J,1)
      ENDDO
      ENDDO
      ENDIF

!     X PERIODICITY
      IF (XPRDIC .EQ. 1) THEN
!$OMP PARALLEL DO
      DO K=0,N3
      DO J=0,N2
         U(0 ,J,K)=U(N1M,J,K)
         U(N1,J,K)=U(1  ,J,K)
      ENDDO
      ENDDO

!$OMP PARALLEL DO
      DO K=0,N3
      DO J=1,N2
         V(0 ,J,K)=V(N1M,J,K)
         V(N1,J,K)=V(1  ,J,K)
      ENDDO
      ENDDO

!$OMP PARALLEL DO
      DO K=1,N3
      DO J=0,N2
         W(0 ,J,K)=W(N1M,J,K)
         W(N1,J,K)=W(1  ,J,K)
      ENDDO
      ENDDO
      ENDIF

      RETURN
      END

!=======================================================================
      SUBROUTINE PCALC(PHI,DIVGSUM)
!=======================================================================
!
!     Calculate pressure from DIVGSUM & PHI
!
!     ({u_i}^k-u_i hat)/(2.*a_k*DT)
!                             = - [ d(p^k)/d(x_i) - d(p^(k-1))/d(x_i) ]
!                               + 0.5*[L({u_i}^k) - L(u_i hat)]
!                             = - d({phi}^k)/d(x_i)
!
!     By the definition of phi, final equation is derived as follows:
!       p^k = p^(k-1) + {phi}^k - (a_k*DT/Re)*(d2({phi^k})/d(x_j)d(x_j))
!
!-----------------------------------------------------------------------
       USE MOD_COMMON
       USE MOD_FLOWARRAY, ONLY : U,V,W,P
       IMPLICIT NONE
       INTEGER*8     :: I,J,K
       REAL*8        :: PHIREF
       REAL*8        :: PHI(0:N1,0:N2,0:N3),DIVGSUM(0:N1,0:N2,0:N3)
       REAL*8        :: FUNCBODY


!$OMP PARALLEL DO
      DO 92 K=1,N3M
      DO 92 J=1,N2M
      DO 92 I=1,N1M
      P(I,J,K)=P(I,J,K)+PHI(I,J,K)-0.5*DTCONST*DIVGSUM(I,J,K)*1./RE
   92 CONTINUE


!     SET THE AVERAGE PHI AT THE UPPER WALL TO BE ZERO.
      PHIREF = 0.
!$OMP PARALLEL DO &
!$OMP reduction(+:PHIREF)
      DO 80 K=1,N3M
      DO 80 J=1,N2M
      PHIREF=PHI(N1M,J,K)*F2FY(J)*F2FZ(K)+PHIREF
   80 CONTINUE
      PHIREF=PHIREF/(YL*ZL)

!$OMP PARALLEL DO
      DO 93 K=1,N3M
      DO 93 J=1,N2M
      DO 93 I=1,N1M
       P(I,J,K)=P(I,J,K)-PHIREF
   93 CONTINUE

!     Z PERIODICITY
      IF (ZPRDIC .EQ. 1) THEN
!$OMP PARALLEL DO
      DO J=1,N2M
      DO I=1,N1M
         P(I,J,0) =P(I,J,N3M)
         P(I,J,N3)=P(I,J,1)
      ENDDO
      ENDDO
      ENDIF

!     X PERIODICITY
      IF (XPRDIC .EQ. 1) THEN
!$OMP PARALLEL DO
      DO K=1,N3M
      DO J=1,N2M
         P(0 ,J,K)=P(N1M,J,K)
         P(N1,J,K)=P(1  ,J,K)
      ENDDO
      ENDDO
      ENDIF

      RETURN
      END

!=======================================================================
      SUBROUTINE TCALC
!=======================================================================
       USE MOD_COMMON
       USE MOD_FLOWARRAY, ONLY : T,CSTAR
       IMPLICIT NONE
       INTEGER*8     :: I,J,K
       REAL*8        :: FUNCBODY,HFLUX,TTEMP

       IF (ICH .EQ. 1) THEN
!$OMP PARALLEL DO
         DO 31 K=1,N3M
         DO 31 J=1,N2M                       
         DO 31 I=1,N1M
            T(I,J,K)=T(I,J,K)-THCAP
   31    CONTINUE
       ENDIF

!        HFLUX = 1.

!        TTEMP = T(1,0,1)

! !$OMP PARALLEL DO
!          DO K=0,N3
!          DO I=0,N1
!            T(I,0,K)=T(I,1,K)+HFLUX*C2CY(1)
!            T(I,N2,K)=T(I,N2M,K)+HFLUX*C2CY(N2)
!          ENDDO
!          ENDDO

!          T = T - (T(1,0,1) - TTEMP)

      RETURN
      END

!=======================================================================
      SUBROUTINE LAGFORCE
!=======================================================================
!
!     Material derivate of velocity, one of major contribution of the
!       force on a body
!
!     Reference
!       Lee et al., 2011, Sources of spurious force oscillations
!       from an immersed boundary method for moving-body problems,
!       J. Comp. Phys., 230, 2677-2695.
!
!     Variables
!       DUDTA,DVDTA,DWDTA: Volume integration of material derivative of
!                          velocity over a body
!       DUDTR: Time derivative of velocity
!       RK3XO,RK3XOO: Convective terms at each RK3 sub-step
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY :  RK3XO,RK3YO,RK3ZO,DUDTR  &
                                ,RK3XOO,RK3YOO,RK3ZOO,IFC,JFC,KFC
      IMPLICIT NONE
      INTEGER*8     :: I,J,K,N,L,II,JJ,KK

      IF(MSUB .EQ. 1) THEN
        DUDTA = 0.
        DVDTA = 0.
        DWDTA = 0.
      ENDIF

!$OMP PARALLEL DO private(II,JJ,KK)
      DO L=1,3
      DO N=1,NBODY(L)
        II=IFC(N,L)
        JJ=JFC(N,L)
        KK=KFC(N,L)
        !Intermediate information
        IF( L .EQ. 1) THEN
          DUDTA=DUDTA+C2CX(II)*F2FY(JJ)*F2FZ(KK)*(DUDTR(N,L)-              &
                (GAMMA(MSUB)*RK3XO(II,JJ,KK)+RO(MSUB)*RK3XOO(II,JJ,KK)))
        ELSEIF( L .EQ. 2) THEN
          DVDTA=DVDTA+F2FX(II)*C2CY(JJ)*F2FZ(KK)*(DUDTR(N,L)-              &
                (GAMMA(MSUB)*RK3YO(II,JJ,KK)+RO(MSUB)*RK3YOO(II,JJ,KK)))
        ELSEIF( L .EQ. 3) THEN
          DWDTA=DWDTA+F2FX(II)*F2FY(JJ)*C2CZ(KK)*(DUDTR(N,L)-              &
                (GAMMA(MSUB)*RK3ZO(II,JJ,KK)+RO(MSUB)*RK3ZOO(II,JJ,KK)))
        ENDIF

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE LAGFORCE
!=======================================================================