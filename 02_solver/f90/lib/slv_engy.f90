!=======================================================================
        SUBROUTINE RHSNLHS_T
!=======================================================================
       USE MOD_COMMON
       IMPLICIT NONE

       IF (ILES .EQ. 1) CALL RHSSGS_T    ! SubGrid-Scale computations (LES case)
       CALL RHS_T
       IF ((IBMON .NE. 0) .AND. ICONJG .EQ. 0) CALL RHS_IBM_T  ! IB method computation
       ! IF ((IBMON .NE. 0) .AND. ICONJG .EQ. 1) CALL RHS_IBM_CONJG_T  ! CONJUGATE H.TRANS

       CALL CONVBC_T                     
       IF (ICH .NE. 1) CALL RHSINCORPBC_T

       CALL LHS_T

       ! WHEN DO BLOWING/SUCTION RETRV_UVW SHOULD BE TURNED ON
       CALL RETRV_T
       ! CALL PRDIC_ADJ_T

       RETURN
       END
!=======================================================================
!=======================================================================
       SUBROUTINE RHS_T
!=======================================================================
!
!     Computing intermediate velocity, u hat, step in delta form
!
!     variables in common:
!           x, y, z         : coordinate direction
!           u, v, w         : velocity for x, y, z direction
!           n, s, e, w, c, f: + & - for each x, y, z direction
!           AN, AL, TAL     : Non-linear, linear, turbulent (SGS)
!
!     RHS1(x,y,z,L):
!           RHS term consists of Non-linear/Linear/SGS terms
!           Components for IB method will be added in RHS_IBM subroutine
!
!-----------------------------------------------------------------------
       USE MOD_COMMON
       USE MOD_FLOWARRAY, ONLY : U,V,W,P,T,ALSGS,ALSGS1,RHS1   &
                                ,RK3TO,RK3TOO,NWALL_DVM,CSTAR,KSTAR
       IMPLICIT NONE
       INTEGER*8   :: I,J,K
       INTEGER*8   :: IPLUS,IMINUS,JPLUS,JMINUS,KPLUS,KMINUS
       REAL*8      :: OMEGA

!------------ variables for T (temperature)
       REAL*8      :: TE,TW,TN,TS,TC,TF
       REAL*8      :: ANT1,ANT2,ANT3,RK3T
       REAL*8      :: ALT1,ALT2,ALT3,ALT4,ALT5,ALT6,ALTX,ALTY,ALTZ,ALT
       REAL*8      :: TALTX,TALTY,TALTZ


!-----RHS1 calculation for T -----------------
!$OMP PARALLEL DO  &
!$OMP private(TE,TW,TN,TS,TC,TF)   &
!$OMP private(ANT1,ANT2,ANT3,RK3T) &
!$OMP private(ALT1,ALT2,ALT3,ALT4,ALT5,ALT6,ALTX,ALTY,ALTZ,ALT) &
!$OMP private(TALTX,TALTY,TALTZ)
       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M

       KPLUS=KPV(K)
       KMINUS=KMV(K)
       JPLUS=JPV(J)
       JMINUS=JMV(J)
       IPLUS=IPV(I)
       IMINUS=IMV(I)

       TE=0.5*(F2FX(I)*T(IPLUS,J,K)+F2FX(IPLUS)*T(I,J,K))*C2CXI(IPLUS) &
          *(1.-FIXIU(I))+T(IPLUS,J,K)*FIXIU(I)
       TW=0.5*(F2FX(IMINUS)*T(I,J,K)+F2FX(I)*T(IMINUS,J,K))*C2CXI(I)   &
          *(1.-FIXIL(I))+T(IMINUS,J,K)*FIXIL(I)

       TN=0.5*(F2FY(J)*T(I,JPLUS,K)+F2FY(JPLUS)*T(I,J,K))*C2CYI(JPLUS) &
          *(1.-FIXJU(J))+T(I,JPLUS,K)*FIXJU(J)
       TS=0.5*(F2FY(JMINUS)*T(I,J,K)+F2FY(J)*T(I,JMINUS,K))*C2CYI(J)   &
          *(1.-FIXJL(J))+T(I,JMINUS,K)*FIXJL(J)

       TC=0.5*(F2FZ(K)*T(I,J,KPLUS)+F2FZ(KPLUS)*T(I,J,K))*C2CZI(KPLUS) &
          *(1.-FIXKU(K))+T(I,J,KPLUS)*FIXKU(K)
       TF=0.5*(F2FZ(KMINUS)*T(I,J,K)+F2FZ(K)*T(I,J,KMINUS))*C2CZI(K)   &
          *(1.-FIXKL(K))+T(I,J,KMINUS)*FIXKL(K)

       ANT1=(TE*U(IPLUS,J,K)-TW*U(I,J,K))*F2FXI(I)
       ANT2=(TN*V(I,JPLUS,K)-TS*V(I,J,K))*F2FYI(J)
       ANT3=(TC*W(I,J,KPLUS)-TF*W(I,J,K))*F2FZI(K)

       IF ((ICONJG .EQ. 1) .AND. (NWALL_DVM(I,J,K) .EQ. 0)) THEN
         OMEGA = 0.
       ELSE
         OMEGA = 1.
       ENDIF

       RK3T=-OMEGA*(ANT1+ANT2+ANT3)              ! Non-linear term at k-substep

       ALT1=(T(IPLUS,J,K)-T(I,J,K))*C2CXI(IPLUS)
       ALT2=(T(I,J,K)-T(IMINUS,J,K))*C2CXI(I)
       ALT3=(T(I,JPLUS,K)-T(I,J,K))*C2CYI(JPLUS)
       ALT4=(T(I,J,K)-T(I,JMINUS,K))*C2CYI(J)  
       ALT5=(T(I,J,KPLUS)-T(I,J,K))*C2CZI(KPLUS)
       ALT6=(T(I,J,K)-T(I,J,KMINUS))*C2CZI(K)   

       ALTX=(ALT1-ALT2)*F2FXI(I)
       ALTY=(ALT3-ALT4)*F2FYI(J)
       ALTZ=(ALT5-ALT6)*F2FZI(K)

       IF (ICONJG .EQ. 1) THEN
         ALTX = CSTAR(I,J,K)*(KSTAR(I,J,K,1)*ALT1-KSTAR(I,J,K,2)*ALT2)*F2FXI(I)
         ALTY = CSTAR(I,J,K)*(KSTAR(I,J,K,3)*ALT3-KSTAR(I,J,K,4)*ALT4)*F2FYI(J)
         ALTZ = CSTAR(I,J,K)*(KSTAR(I,J,K,5)*ALT5-KSTAR(I,J,K,6)*ALT6)*F2FZI(K)
       ENDIF

       ALT=1./(RE*PR)*(ALTX+ALTY+ALTZ)           ! Linear terms at k-substep
       
!-----LES
      IF (ILES.EQ.1) THEN
       TALTX=F2FXI(I)*                                           &
             (ALSGS1(IPLUS,J,K,1)*C2CXI(IPLUS)*(T(IPLUS,J,K)-T(I,J,K)) &
             -ALSGS1(I,J,K,1)*C2CXI(I)*(T(I,J,K)-T(I-1,J,K))) 
       TALTY=F2FYI(J)*                                           &
             (ALSGS1(I,JPLUS,K,2)*C2CYI(JPLUS)*(T(I,JPLUS,K)-T(I,J,K)) &
             -ALSGS1(I,J,K,2)*C2CYI(J)*(T(I,J,K)-T(I,J-1,K)))
       TALTZ=F2FZI(K)*                                           &
             (ALSGS1(I,J,KPLUS,3)*C2CZI(KPLUS)*(T(I,J,KPLUS)-T(I,J,K)) &
             -ALSGS1(I,J,K,3)*C2CZI(K)*(T(I,J,K)-T(I,J,K-1)))
       ALT=ALT+FLOAT(ILES)*(TALTX+TALTY+TALTZ)
       RK3T=RK3T+FLOAT(ILES)*RHS1(I,J,K,4)
      ENDIF
!-----LES

       RHS1(I,J,K,4)=DT                                    &
               *( GAMMA(MSUB)*RK3T+RO(MSUB)*RK3TO(I,J,K)   &
               +2.*ALPHA*ALT )
       RK3TOO(I,J,K)=RK3TO(I,J,K)
       RK3TO(I,J,K)=RK3T

       ENDDO
       ENDDO
       ENDDO

       RETURN
       END
!=======================================================================
!=======================================================================
      SUBROUTINE RHS_IBM_T
!=======================================================================
!
!     Calculate momentum forcing
!
!     Option
!       IMOVINGON = 0, Stationary body => UBODY,VBODY,WBODY for translational vel.
!       IMOVINGON = 1, Moving body     => UBD,VBD,WBD in lica_cylinder.f90
!
!     Variables
!       UTARG,VTARG,WTARG: Target velocity to satisfy no-slip b.c.
!       FCV   : Momentum forcing from the target velocity
!       FCVAVG: Averaging forcing values    to calculate the force on a body
!       DUDTR : Time derivative of velocity to calculate the force on a body
!               Ref. Lee et al., 2011, Sources of spurious force
!               oscillations from an immersed boundary method for moving
!               -body problems, J. Comp. Phys., 230, 2677-2695.
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : RHS1,U,V,W,T,IFC,JFC,KFC,INTPINDX,GEOMFAC,&
                                FCV,FCVAVG,DUDTR
      IMPLICIT NONE
      INTEGER*8 :: I,J,K,N,II,JJ,KK,IP,JP,KP,IPP,JPP,KPP
      REAL*8    :: TBODY,DTI
      REAL*8    :: TTARG
      REAL*8    :: RHSTMP(N1M,N2M,N3M)

!---- Compute target velocities & forcing values at forcing points
      TBODY = 0.
      DTI   = 1./DT

!$OMP PARALLEL DO
      DO K= 1,N3M
      DO J= 1,N2M
      DO I= 1,N1M
        RHSTMP(I,J,K)=RHS1(I,J,K,4)
      ENDDO
      ENDDO
      ENDDO


!$OMP PARALLEL DO &
!$OMP private(II,JJ,KK,IP,JP,KP,IPP,JPP,KPP,TTARG,TBODY)
      DO N=1,NINTP(4)
         II =IFC(N,4)
         JJ =JFC(N,4)
         KK =KFC(N,4)
         IP =IFC(N,4)+INTPINDX(N,4,1)
         JP =JFC(N,4)+INTPINDX(N,4,2)
         KP =KFC(N,4)+INTPINDX(N,4,3)
         IPP=IFC(N,4)+INTPINDX(N,4,1)*2
         JPP=JFC(N,4)+INTPINDX(N,4,2)*2
         KPP=KFC(N,4)+INTPINDX(N,4,3)*2

         IF ((IPP.GE.N1).OR.(IPP.LE.0)) CALL REINDEX_I(IP,IPP) ! AT SLV_MMTM LIB
         IF ((KPP.GE.N3).OR.(KPP.LE.0)) CALL REINDEX_K(KP,KPP) ! AT SLV_MMTM LIB
         TTARG=GEOMFAC(N,4,0,0,0)* TBODY                                &
              +GEOMFAC(N,4,0,0,1)*(T(II ,JJ ,KP )+RHSTMP(II ,JJ ,KP ))  &
              +GEOMFAC(N,4,0,0,2)*(T(II ,JJ ,KPP)+RHSTMP(II ,JJ ,KPP))  &
              +GEOMFAC(N,4,0,1,0)*(T(II ,JP ,KK )+RHSTMP(II ,JP ,KK ))  &
              +GEOMFAC(N,4,0,1,1)*(T(II ,JP ,KP )+RHSTMP(II ,JP ,KP ))  &
              +GEOMFAC(N,4,0,1,2)*(T(II ,JP ,KPP)+RHSTMP(II ,JP ,KPP))  &
              +GEOMFAC(N,4,0,2,0)*(T(II ,JPP,KK )+RHSTMP(II ,JPP,KK ))  &
              +GEOMFAC(N,4,0,2,1)*(T(II ,JPP,KP )+RHSTMP(II ,JPP,KP ))  &
              +GEOMFAC(N,4,0,2,2)*(T(II ,JPP,KPP)+RHSTMP(II ,JPP,KPP))  &
              +GEOMFAC(N,4,1,0,0)*(T(IP ,JJ ,KK )+RHSTMP(IP ,JJ ,KK ))  &
              +GEOMFAC(N,4,1,0,1)*(T(IP ,JJ ,KP )+RHSTMP(IP ,JJ ,KP ))  &
              +GEOMFAC(N,4,1,0,2)*(T(IP ,JJ ,KPP)+RHSTMP(IP ,JJ ,KPP))  &
              +GEOMFAC(N,4,1,1,0)*(T(IP ,JP ,KK )+RHSTMP(IP ,JP ,KK ))  &
              +GEOMFAC(N,4,1,1,1)*(T(IP ,JP ,KP )+RHSTMP(IP ,JP ,KP ))  &
              +GEOMFAC(N,4,1,1,2)*(T(IP ,JP ,KPP)+RHSTMP(IP ,JP ,KPP))  &
              +GEOMFAC(N,4,1,2,0)*(T(IP ,JPP,KK )+RHSTMP(IP ,JPP,KK ))  &
              +GEOMFAC(N,4,1,2,1)*(T(IP ,JPP,KP )+RHSTMP(IP ,JPP,KP ))  &
              +GEOMFAC(N,4,1,2,2)*(T(IP ,JPP,KPP)+RHSTMP(IP ,JPP,KPP))  &
              +GEOMFAC(N,4,2,0,0)*(T(IPP,JJ ,KK )+RHSTMP(IPP,JJ ,KK ))  &
              +GEOMFAC(N,4,2,0,1)*(T(IPP,JJ ,KP )+RHSTMP(IPP,JJ ,KP ))  &
              +GEOMFAC(N,4,2,0,2)*(T(IPP,JJ ,KPP)+RHSTMP(IPP,JJ ,KPP))  &
              +GEOMFAC(N,4,2,1,0)*(T(IPP,JP ,KK )+RHSTMP(IPP,JP ,KK ))  &
              +GEOMFAC(N,4,2,1,1)*(T(IPP,JP ,KP )+RHSTMP(IPP,JP ,KP ))  &
              +GEOMFAC(N,4,2,1,2)*(T(IPP,JP ,KPP)+RHSTMP(IPP,JP ,KPP))  &
              +GEOMFAC(N,4,2,2,0)*(T(IPP,JPP,KK )+RHSTMP(IPP,JPP,KK ))  &
              +GEOMFAC(N,4,2,2,1)*(T(IPP,JPP,KP )+RHSTMP(IPP,JPP,KP ))  &
              +GEOMFAC(N,4,2,2,2)*(T(IPP,JPP,KPP)+RHSTMP(IPP,JPP,KPP))
         FCV(N,4)=(TTARG-T(II,JJ,KK)-RHSTMP(II,JJ,KK))*DTI
         RHS1(II,JJ,KK,4)=TTARG-T(II,JJ,KK)
      ENDDO

!$OMP PARALLEL DO &
!$OMP private(II,JJ,KK,TTARG,TBODY)
      DO N=NINTP(4)+1,NBODY(4)
         II=IFC(N,4)
         JJ=JFC(N,4)
         KK=KFC(N,4)
         TTARG=TBODY
         FCV(N,4)=(TTARG-T(II,JJ,KK)-RHS1(II,JJ,KK,4))*DTI
         RHS1(II,JJ,KK,4)=TTARG-T(II,JJ,KK)
      ENDDO

!-----compute average forcing values during RK3 steps
!$OMP PARALLEL DO
      DO N=1,NBODY(4)
         FCVAVG(N,4)=FCVAVG(N,4)+FCV(N,4)
      ENDDO

      RETURN
      END
!=======================================================================
!=======================================================================
      SUBROUTINE RHS_IBM_CONJG_T
!=======================================================================
!
!     Calculate momentum forcing
!
!     Option
!       IMOVINGON = 0, Stationary body => UBODY,VBODY,WBODY for translational vel.
!       IMOVINGON = 1, Moving body     => UBD,VBD,WBD in lica_cylinder.f90
!
!     Variables
!       UTARG,VTARG,WTARG: Target velocity to satisfy no-slip b.c.
!       FCV   : Momentum forcing from the target velocity
!       FCVAVG: Averaging forcing values    to calculate the force on a body
!       DUDTR : Time derivative of velocity to calculate the force on a body
!               Ref. Lee et al., 2011, Sources of spurious force
!               oscillations from an immersed boundary method for moving
!               -body problems, J. Comp. Phys., 230, 2677-2695.
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : RHS1,U,V,W,T,IFC,JFC,KFC,INTPINDX,GEOMFAC,&
                                FCV,FCVAVG,DUDTR,NWALL_DVM
      IMPLICIT NONE
      INTEGER*8 :: I,J,K,N,KPLUS,KMINUS,JPLUS,JMINUS,IPLUS,IMINUS
      REAL*8    :: TP,TM,HEATS
      REAL*8    :: RHSTMP(N1M,N2M,N3M)

!---- Compute target velocities & forcing values at forcing points

!$OMP PARALLEL DO &
!$OMP private(TP,TM,HEATS,KPLUS,KMINUS,JPLUS,JMINUS,IPLUS,IMINUS)
      DO K= 1,N3M 
      DO J= 1,N2M
      DO I= 1,N1M
        IF (NWALL_DVM(I,J,K).EQ.0) THEN

          HEATS = 0.

          KPLUS=KPV(K)
          KMINUS=KMV(K)
          JPLUS=JPV(J)
          JMINUS=JMV(J)
          IPLUS=IPV(I)
          IMINUS=IMV(I)

          TP=0.5*(F2FX(I)*T(IPLUS,J,K)+F2FX(IPLUS)*T(I,J,K))*C2CXI(IPLUS) &
             *(1.-FIXIU(I))+T(IPLUS,J,K)*FIXIU(I)
          TM=0.5*(F2FX(IMINUS)*T(I,J,K)+F2FX(I)*T(IMINUS,J,K))*C2CXI(I)   &
             *(1.-FIXIL(I))+T(IMINUS,J,K)*FIXIL(I)
          HEATS = HEATS + (TP*U(IPLUS,J,K) - TM*U(I,J,K)) * F2FXI(I)

          TP=0.5*(F2FY(J)*T(I,JPLUS,K)+F2FY(JPLUS)*T(I,J,K))*C2CYI(JPLUS) &
             *(1.-FIXJU(J))+T(I,JPLUS,K)*FIXJU(J)
          TM=0.5*(F2FY(JMINUS)*T(I,J,K)+F2FY(J)*T(I,JMINUS,K))*C2CYI(J)   &
             *(1.-FIXJL(J))+T(I,JMINUS,K)*FIXJL(J)
          HEATS = HEATS + (TP*V(I,JPLUS,K) - TM*V(I,J,K)) * F2FYI(J)

          TP=0.5*(F2FZ(K)*T(I,J,KPLUS)+F2FZ(KPLUS)*T(I,J,K))*C2CZI(KPLUS) &
             *(1.-FIXKU(K))+T(I,J,KPLUS)*FIXKU(K)
          TM=0.5*(F2FZ(KMINUS)*T(I,J,K)+F2FZ(K)*T(I,J,KMINUS))*C2CZI(K)   &
             *(1.-FIXKL(K))+T(I,J,KMINUS)*FIXKL(K)
          HEATS = HEATS + (TP*W(I,J,KPLUS) - TM*W(I,J,K)) * F2FZI(K)

          RHS1(I,J,K,4) = RHS1(I,J,K,4) + 2.*ALPHA*DT*HEATS
        ENDIF
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END
!=======================================================================
!=======================================================================
        SUBROUTINE CONVBC_T
!=======================================================================
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,T
      IMPLICIT NONE
      INTEGER*8     :: I,J,K
      REAL*8      :: QIN,QOUT,QRATIO
      REAL*8      :: UBAR,UCOEF,VWCOEF

      IF (ICH .NE. 1) THEN

!-----CALCULATE INFLUX
      QIN=0.
!$OMP PARALLEL DO &
!$OMP REDUCTION(+:QIN)
      DO K=1,N3M
        DO J=1,N2M
          QIN=U(1,J,K)*F2FY(J)*F2FZ(K)+QIN
        ENDDO
      ENDDO
!$OMP PARALLEL DO &
!$OMP REDUCTION(+:QIN)
      DO K=1,N3M
        DO I=1,N1M
          QIN=(V(I,1,K)-V(I,N2,K))*F2FX(I)*F2FZ(K)+QIN
        ENDDO
      ENDDO
!$OMP PARALLEL DO &
!$OMP REDUCTION(+:QIN)
      DO J=1,N2M
        DO I=1,N1M
         QIN=(W(I,J,1)-W(I,J,N3))*F2FX(I)*F2FY(J)+QIN
        ENDDO
      ENDDO

!-----CALCULATE CONVECTIVE VELOCITY Uc
      UBAR = 0.
!$OMP PARALLEL DO &
!$OMP REDUCTION(+:UBAR)
      DO K=1,N3M
        DO J=1,N2M
          UBAR=U(N1,J,K)*F2FY(J)*F2FZ(K)+UBAR
        ENDDO
      ENDDO
      UBAR  = UBAR/(YL*ZL)
      UCOEF = UBAR*DTCONST*F2FXI(N1M)
      VWCOEF= UBAR*DTCONST*C2CXI(N1)

      QOUT = 0.
!$OMP PARALLEL DO &
!$OMP REDUCTION(+:QOUT)
      DO K=1,N3M
        DO J=1,N2M
          UOUT(J,K)=U(N1,J,K)-UCOEF*(U(N1,J,K)-U(N1M,J,K))
          TOUT(J,K)=T(N1,J,K)-VWCOEF*(T(N1,J,K)-T(N1M,J,K))
          QOUT=UOUT(J,K)*F2FY(J)*F2FZ(K)+QOUT
        ENDDO
      ENDDO
      QRATIO=QIN/QOUT

!-----ADJUST BOUNDARY VELOCITY TO SATISFY GLOBAL MASS CONSERVATION
!$OMP PARALLEL DO
      DO K=1,N3M
        DO J=1,N2M
          TOUT(J,K)=TOUT(J,K)*QRATIO
          DTOUT(J,K)=TOUT(J,K)-T(N1,J,K)
        ENDDO
      ENDDO
   
      IF (YPRDIC .EQ. 1) THEN
!$OMP PARALLEL DO
        DO K=0,N3
          DTOUT(N2,K)=DTOUT(1,  K)
          DTOUT(0, K)=DTOUT(N2M,K)
        ENDDO
      ENDIF
     
      IF (ZPRDIC .EQ. 1) THEN
!$OMP PARALLEL DO
        DO J=0,N2
          DTOUT(J,N3)=DTOUT(J,  1)
          DTOUT(J, 0)=DTOUT(J,N3M)
        ENDDO
      ENDIF

      ELSE

      TOUT = 0.

      ENDIF

      RETURN
      END SUBROUTINE CONVBC_T
!=======================================================================
!=======================================================================
       SUBROUTINE RHSINCORPBC_T
!=======================================================================
       USE MOD_COMMON
       USE MOD_FLOWARRAY, ONLY : U,V,W,T,RHS1,ALSGS,CSTAR,KSTAR
       IMPLICIT NONE
       INTEGER*8     :: I,J,K
       REAL*8        :: CRE,CRE2,CS,KS

       IF (ILES .EQ. 1) THEN
         CRE=RE*PR
         CRE2=2.*RE*PR
       ELSE
         CRE=0.
         CRE2=0.
       ENDIF

!$OMP PARALLEL DO
       DO K=1,N3M
       DO J=1,N2M
       IF (ICONJG .EQ. 1) THEN
        CS = CSTAR(N1M,J,K)
        KS = (KSTAR(N1M,J,K,1)+KSTAR(N1M,J,K,2)+KSTAR(N1M,J,K,3)        &
             +KSTAR(N1M,J,K,4)+KSTAR(N1M,J,K,5)+KSTAR(N1M,J,K,6))/6.
       ELSE
        CS = 1.
        KS = 1.
       ENDIF

       RHS1(N1M,J,K,4)=RHS1(N1M,J,K,4)                                  &
                      -ACOEF*CS*KS                                      &
                      *CIU(N1M)*(1.+CRE2*ALSGS(N1M,J,K))*DTOUT(J,K)
       ENDDO
       ENDDO

       RETURN
       END SUBROUTINE RHSINCORPBC_T
!=======================================================================
!=======================================================================
      SUBROUTINE LHS_T
!=======================================================================
!
!     Calculate intermediate velocity, u_i hat through TDMA
!           In this routine, compute streamwise velocity (u hat)
!
!-----------------------------------------------------------------------
       USE MOD_COMMON
       USE MOD_FLOWARRAY, ONLY : U,V,W,T,RHS1,ALSGS,ALSGS1,CSTAR,KSTAR
       IMPLICIT NONE
       INTEGER*8     :: I,J,K
       REAL*8        :: CRE,CRE2,CS,KS1,KS2
       REAL*8, DIMENSION (:,:), ALLOCATABLE :: AI,BI,CI,GI
       REAL*8, DIMENSION (:,:), ALLOCATABLE :: AJ,BJ,CJ,GJ
       REAL*8, DIMENSION (:,:), ALLOCATABLE :: AK,BK,CK,GK

       CRE=FLOAT(ILES)*RE*PR
       CRE2=2.*FLOAT(ILES)*RE*PR

!=====ADI STARTS

      IF (N3M.EQ.1) GOTO 100
!-----Z-DIRECTION
!$OMP PARALLEL &
!$OMP private(AK,CK,BK,GK,CS,KS1,KS2)
      allocate(AK(N1,N3),BK(N1,N3),CK(N1,N3),GK(N1,N3))
!$OMP DO
      DO 131 J=1,N2M
      DO 141 K=1,N3M
      DO 141 I=1,N1M
       IF (ICONJG .EQ. 1) THEN
        CS = CSTAR(I,J,K)
        KS1 = KSTAR(I,J,K,6)
        KS2 = KSTAR(I,J,K,5)
       ELSE
        CS = 1.
        KS1 = 1.
        KS2 = 1.
       ENDIF
       AK(I,K)=AKUV(K)*(CS*KS1+CRE*ALSGS1(I,J,K,2)) 
       CK(I,K)=CKUV(K)*(CS*KS2+CRE*ALSGS1(I,J,K+1,2))
       BK(I,K)=ACOEFI*PR-AK(I,K)-CK(I,K)
       GK(I,K)=ACOEFI*PR*RHS1(I,J,K,4)
  141 CONTINUE

      IF (ZPRDIC .EQ. 0) THEN
         CALL TRDIAG3(AK,BK,CK,GK,GK,1,N3M,1,N1M)
      ELSE IF (ZPRDIC .EQ. 1) THEN
         CALL TRDIAG3P(AK,BK,CK,GK,1,N3M,1,N1M)  !z periodicity
      ENDIF

      DO 151 K=1,N3M
      DO 151 I=1,N1M
      RHS1(I,J,K,4)=GK(I,K)
  151 CONTINUE
  131 CONTINUE
!$OMP END DO
      deallocate(AK,BK,CK,GK)
!$OMP END PARALLEL

  100 CONTINUE

!$OMP PARALLEL  &
!$OMP private(AJ,CJ,BJ,GJ)  &
!$OMP private(AI,CI,BI,GI,CS,KS1,KS2)
      allocate(AI(N2,N1),BI(N2,N1),CI(N2,N1),GI(N2,N1))
      allocate(AJ(N1,N2),BJ(N1,N2),CJ(N1,N2),GJ(N1,N2))
!$OMP DO
      DO 40 K=1,N3M

!-----Y-DIRECTION
      DO 91 J=1,N2M
      DO 91 I=1,N1M
       IF (ICONJG .EQ. 1) THEN
        CS = CSTAR(I,J,K)
        KS1 = KSTAR(I,J,K,4)
        KS2 = KSTAR(I,J,K,3)
       ELSE
        CS = 1.
        KS1 = 1.
        KS2 = 1.
       ENDIF
      AJ(I,J)=AJUW(J)*(CS*KS1+CRE*ALSGS1(I,J,K,3))  
      CJ(I,J)=CJUW(J)*(CS*KS2+CRE*ALSGS1(I,J+1,K,3))
      BJ(I,J)=ACOEFI*PR-AJ(I,J)-CJ(I,J)
      GJ(I,J)=ACOEFI*PR*RHS1(I,J,K,4)
   91 CONTINUE

      CALL TRDIAG2(AJ,BJ,CJ,GJ,GJ,1,N2M,1,N1M)

!-----X-DIRECTION
      DO 51 I=1,N1M
      DO 51 J=1,N2M
       IF (ICONJG .EQ. 1) THEN
        CS = CSTAR(I,J,K)
        KS1 = KSTAR(I,J,K,2)
        KS2 = KSTAR(I,J,K,1)
       ELSE
        CS = 1.
        KS1 = 1.
        KS2 = 1.
       ENDIF
      AI(J,I)=AIVW(I)*(CS*KS1+CRE*ALSGS1(I,J,K,1))
      CI(J,I)=CIVW(I)*(CS*KS2+CRE*ALSGS1(I+1,J,K,1))
      BI(J,I)=ACOEFI*PR-AI(J,I)-CI(J,I)
      GI(J,I)=ACOEFI*PR*GJ(I,J)
   51 CONTINUE

      IF (XPRDIC .EQ. 0) THEN
         CALL TRDIAG1(AI,BI,CI,GI,GI,1,N1M,1,N2M)
      ELSE IF (XPRDIC .EQ. 1) THEN
         CALL TRDIAG1P(AI,BI,CI,GI,1,N1M,1,N2M)
      ENDIF

      DO 61 J=1,N2M
      DO 61 I=1,N1M
        T(I,J,K)=GI(J,I)+T(I,J,K)
   61 CONTINUE

   40 CONTINUE
!$OMP END DO
      deallocate(AI,BI,CI,GI)
      deallocate(AJ,BJ,CJ,GJ)
!$OMP END PARALLEL

      RETURN
      END
!=======================================================================
!=======================================================================
      SUBROUTINE RETRV_T
!=======================================================================
!
!     Adjust boundary velocity to satisfy global mass conservation
!           [retrieve uvw (at the exit boundary)]
!     See 'CONVBC' subroutine in lica_[bodyname, e.g. cylinder].f90 file
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,T
      IMPLICIT NONE
      INTEGER*8    :: I,J,K

!$OMP PARALLEL DO
      DO K=1,N3M
      DO J=1,N2M
      T(N1,J,K)=TOUT(J,K)
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE RETRV_T
!=======================================================================
!=======================================================================
      SUBROUTINE PRDIC_ADJ_T
!=======================================================================
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : T
      IMPLICIT NONE
      INTEGER*8    :: I,J,K

!     X PERIODICITY
      IF (XPRDIC .EQ. 1) THEN
!$OMP PARALLEL DO
      DO K=0,N3
      DO J=0,N2
         T(0 ,J,K)=T(N1M,J,K)
         T(N1,J,K)=T(1  ,J,K)
      ENDDO
      ENDDO
      ENDIF

!     Y PERIODICITY
      IF (YPRDIC .EQ. 1) THEN
!$OMP PARALLEL DO
      DO K=1,N3
      DO I=0,N1
         T(I ,0,K)=T(I,N2M,K)
         T(I,N2,K)=T(I,  1,K)
      ENDDO
      ENDDO
      ENDIF

!     Z PERIODICITY
      IF (ZPRDIC .EQ. 1) THEN
!$OMP PARALLEL DO
      DO J=0,N2
      DO I=1,N1
         T(I,J,0) =T(I,J,N3M)
         T(I,J,N3)=T(I,J,1)
      ENDDO
      ENDDO
      ENDIF

      RETURN
      END SUBROUTINE PRDIC_ADJ_T
!=======================================================================