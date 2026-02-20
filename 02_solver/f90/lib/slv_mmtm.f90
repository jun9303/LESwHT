!=======================================================================
       SUBROUTINE MEANPG
!=======================================================================
!
!     The mean pressure gradient (dP/dx) required to
!     keep the mass flow rate constant is determined
!     by integrating the wall-shear stresses at the channel walls.
!
!     When there is an IBM body inside the channel, 
!     the IBM forcing has to be included.
!     -> FCV(N,1) is considered (forcing in the streamwise direction)
!
!     PMI(0): overall mean pressure gradient
!     PMI(1): p.grad component at the y-bottom wall
!     PMI(2): p.grad component at the y-top    wall
!     PMI(3): p.grad component at the z-bottom wall
!     PMI(4): p.grad component at the z-top    wall
!
!-----------------------------------------------------------------------
!$    use omp_lib
      USE MOD_COMMON
      USE MOD_FLOWARRAY
      IMPLICIT NONE
      INTEGER*8         :: I,J,K,N
      REAL*8            :: TMP, VOLUME

      PMI = 0.

      IF (BC_YBTM .EQ. 0) THEN
      TMP = 0.
!$OMP PARALLEL DO private(TMP) reduction(+:PMI)
      DO K=1,N3M
        DO I=1,N1M
          TMP = (U(I,1,K) - 0.)*C2CYI(1)*C2CX(I)*F2FZ(K)
          PMI(1) = PMI(1) + TMP/RE
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ENDIF

      IF (BC_YTOP .EQ. 0) THEN
      TMP = 0.
!$OMP PARALLEL DO private(TMP) reduction(+:PMI)
      DO K=1,N3M
        DO I=1,N1M
          TMP = (U(I,N2M,K) - 0.)*C2CYI(N2)*C2CX(I)*F2FZ(K)
          PMI(2) = PMI(2) + TMP/RE
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ENDIF

      IF (BC_ZBTM .EQ. 0) THEN
      TMP = 0.
!$OMP PARALLEL DO private(TMP) reduction(+:PMI)
      DO J=1,N2M
        DO I=1,N1M
          TMP = (U(I,J,1) - 0.)*C2CZI(1)*C2CX(I)*F2FY(J)
          PMI(3) = PMI(3) + TMP/RE
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ENDIF

      IF (BC_ZTOP .EQ. 0) THEN
      TMP = 0.
!$OMP PARALLEL DO private(TMP) reduction(+:PMI)
      DO J=1,N2M
        DO I=1,N1M
          TMP = (U(I,J,N3M) - 0.)*C2CZI(N3)*C2CX(I)*F2FZ(J)
          PMI(4) = PMI(4) + TMP/RE
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ENDIF

      PMI(0) = PMI(1) + PMI(2) + PMI(3) + PMI(4)

      IF ((BC_YBTM .NE. 0) .AND. (BC_YTOP .NE. 0) .AND. &
          (BC_ZBTM .NE. 0) .AND. (BC_ZTOP .NE. 0)) THEN
        WRITE(*,*) ' TO USE ICH CONDITION, AT LEAST ONE OF Y,Z-WALLS'
        WRITE(*,*) ' SHOULD BE WALL SO THAT WALL-SHEAR STRESS EXISTS.'
        STOP
      ELSE
        TMP = 0.
!$OMP PARALLEL DO reduction(+:TMP,VOLUME)
        DO N=1,NBODY(1)
          TMP = TMP + FCV(N,1)*C2CX(IFC(N,1))*F2FY(JFC(N,1))&
                              *F2FZ(KFC(N,1))
        ENDDO
!$OMP END PARALLEL DO
        PMI(0) = PMI(0) + TMP
        VOLUME = XL*YL*ZL

        DO I = 0,4
          PMI(I) = PMI(I) / VOLUME
          ! Remove the fractional RK3 scaling so PMI evaluates to the true full-step force
          IF (ALPHA .GT. 0.0D0) THEN
             PMI(I) = PMI(I) / (2.0D0 * ALPHA)
          ENDIF
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE MEANPG
!=======================================================================
!=======================================================================
       SUBROUTINE LHSINIT
!=======================================================================
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8 :: IC,JC,KC
      INTEGER*8 :: IM,IP,JM,JP,KM,KP

      IF (XPRDIC .EQ. 1) THEN
        I_BGPX=1
      ELSE
        I_BGPX=2
      ENDIF

      IF (YPRDIC .EQ. 1) THEN
        J_BGPY=1
      ELSE 
        J_BGPY=2
      ENDIF

      IF (ZPRDIC .EQ. 1) THEN
        K_BGPZ=1
      ELSE 
        K_BGPZ=2
      ENDIF

!$OMP PARALLEL DO
      DO IC=I_BGPX,N1M
      IM=IMV(IC)
      AIU(IC)=-C2CXI(IC)*F2FXI(IM)
      CIU(IC)=-C2CXI(IC)*F2FXI(IC)
      ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
      DO IC=1,N1M
      IP=IPV(IC)
      AIVW(IC)=-C2CXI(IC)*F2FXI(IC)
      CIVW(IC)=-C2CXI(IP)*F2FXI(IC)
      ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
      DO JC=J_BGPY,N2M
      JM=JMV(JC)
      AJV(JC)=-C2CYI(JC)*F2FYI(JM)
      CJV(JC)=-C2CYI(JC)*F2FYI(JC)
      ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
      DO JC=1,N2M
      JP=JPV(JC)
      AJUW(JC)=-C2CYI(JC)*F2FYI(JC)
      CJUW(JC)=-C2CYI(JP)*F2FYI(JC)
      ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
      DO KC=K_BGPZ,N3M
      KM=KMV(KC)
      AKW(KC)=-C2CZI(KC)*F2FZI(KM)
      CKW(KC)=-C2CZI(KC)*F2FZI(KC)
      ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
      DO KC=1,N3M
      KP=KPV(KC)
      AKUV(KC)=-C2CZI(KC)*F2FZI(KC)
      CKUV(KC)=-C2CZI(KP)*F2FZI(KC)
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE LHSINIT
!=======================================================================
!=======================================================================
        SUBROUTINE RHSNLHS
!=======================================================================
!
!     Solving N-S equation to get intermediate velocity, u_i hat
!     Intermediate velocity, u_i hat: {u_i}^(k-1) -> u_i hat -> {u_i}^k
!           is based on the pressure at timestep k (p^k).
!           Indeed, u_i hat is not satisfied with continuity equation.
!           For details about u_i hat, see 'Fractional step method'.
!
!     RHS sub-grid scale part (LES) -> RHS -> RHS IB(Immersed boundary)
!     -> Computing convective boundary condition (exit condition)
!     -> Other boundary conditions -> LHS (Get u_i hat)
!
!     Schemes used in this code
!       Time : 3rd order Runge-Kutta (Convection/Non-linear term) +
!              2nd order Crank-Nicolson (Diffusion/Linear term)
!       Space: 2nd order Central-difference
!
!     Non-linear term: N(u_i)=  d(u_i*u_j)/d(x_j)
!     Linear term    : L(u_i)=( d(d(u_i))/(d(x_j)d(x_j)) )/Re
!
!
!     See the following papers for details:
!
!     Papers for Fractional step method (computation procedure):
!      Kim, J., Moin, P. & Moser, R. 1987 Turbulence statistics in
!        fully developed channel flow at low Reynolds number.
!        J. Fluid Mech. 177, 133.
!      Choi, H., Moin, P. & Kim, J. 1994 Active turbulence control for
!        drag reduction in wall-bounded flows. J. Fluid Mech. 262,
!        75-110.
!
!     Paper for Large eddy simulation (Dynamic global model):
!      Lee, J., Choi, H. & Park, N. 2010 Dynamic global model for
!        large eddy simulation of transient flow. Phys. Fluids 22,
!        075106.
!
!     Paper for Immersed Boundary(IB) method:
!      Kim, J., Kim, D. & Choi, H. 2001 An Immersed-Boundary Finite
!        Volume Method for Simulations of Flow in Complex Geometries.
!        J. Comput. Phys. 171 132-150.
!
!-----------------------------------------------------------------------
       USE MOD_COMMON
       IMPLICIT NONE

       IF (ILES .EQ. 1) CALL RHSSGS    ! SubGrid-Scale computations (LES case)
       CALL RHS
       IF (IBMON .NE. 0) CALL RHS_IBM    ! IB method computation
       CALL CONVBC                     ! defined in lica_[bodyname, e.g. cylinder].f90 file
       IF (ICH .NE. 1) CALL RHSINCORPBC

       CALL LHSU
       CALL LHSV
       IF (N3M.NE.1) CALL LHSW

       ! WHEN DO BLOWING/SUCTION RETRV_UVW SHOULD BE TURNED ON
       CALL RETRV_UVW
       CALL PRDIC_ADJ_UVW(0)

       RETURN
       END
!=======================================================================
       SUBROUTINE RHS
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
       USE MOD_FLOWARRAY, ONLY : U,V,W,P,T,NUSGS,NUSGS1,RHS1   &
                                ,RK3XO,RK3YO,RK3ZO,RK3XOO,RK3YOO,RK3ZOO
       IMPLICIT NONE
       INTEGER*8   :: I,J,K

!------------ variables for U (x, streamwise) velocity
       REAL*8      :: UN,US,UC,UF,VNX,VSX,WCX,WFX
       REAL*8      :: ANXE,ANXW,ANXX,ANX1,ANX2,RK3X
       REAL*8      :: ALX1,ALX2,ALX3,ALX4,ALX5,ALX6,ALXX,ALXY,ALXZ,ALX
       REAL*8      :: TALXX,TALXY,TALXZ
!------------ variables for V (y, traseverse) velocity
       REAL*8      :: VE,VW,VC,VF,UEY,UWY,WCY,WFY
       REAL*8      :: ANYN,ANYS,ANYY,ANY1,ANY2,RK3Y
       REAL*8      :: ALY1,ALY2,ALY3,ALY4,ALY5,ALY6,ALYX,ALYY,ALYZ,ALY
       REAL*8      :: TALYX,TALYY,TALYZ
!------------ variables for W (z, spanwise) velocity
       REAL*8      :: WE,WW,WN,WS,UEZ,UWZ,VNZ,VSZ
       REAL*8      :: ANZC,ANZF,ANZZ,ANZ1,ANZ2,RK3Z
       REAL*8      :: ALZ1,ALZ2,ALZ3,ALZ4,ALZ5,ALZ6,ALZX,ALZY,ALZZ,ALZ
       REAL*8      :: TALZX,TALZY,TALZZ

!-----RHS1 calculation for u momentum -----------------
!$OMP PARALLEL DO  &
!$OMP private(UN,US,UC,UF,VNX,VSX,WCX,WFX)   &
!$OMP private(ANXE,ANXW,ANXX,ANX1,ANX2,RK3X) &
!$OMP private(ALX1,ALX2,ALX3,ALX4,ALX5,ALX6) &
!$OMP private(ALXX,ALXY,ALXZ,ALX,TALXX,TALXY,TALXZ)
       DO K=1,N3M
       DO J=1,N2M
       DO I=I_BGPX,N1M

       UN=0.5*(F2FY(J+1)*U(I,J,K)+F2FY(J)*U(I,J+1,K))*C2CYI(J+1)    &
           *(1.-FIXJU(J))+U(I,N2,K)*FIXJU(J)
       US=0.5*(F2FY(J)*U(I,J-1,K)+F2FY(J-1)*U(I,J,K))*C2CYI(J)      &
           *(1.-FIXJL(J))+U(I,0,K)*FIXJL(J)
       UC=0.5*(F2FZ(K+1)*U(I,J,K)+F2FZ(K)*U(I,J,K+1))*C2CZI(K+1)    &
           *(1.-FIXKU(K))+U(I,J,N3)*FIXKU(K)
       UF=0.5*(F2FZ(K)*U(I,J,K-1)+F2FZ(K-1)*U(I,J,K))*C2CZI(K)      &
           *(1.-FIXKL(K))+U(I,J,0)*FIXKL(K)
       VNX=0.5*(F2FX(I-1)*V(I,J+1,K)+F2FX(I)*V(I-1,J+1,K))*C2CXI(I)
       VSX=0.5*(F2FX(I-1)*V(I,J,K)+F2FX(I)*V(I-1,J,K))*C2CXI(I)
       WCX=0.5*(F2FX(I-1)*W(I,J,K+1)+F2FX(I)*W(I-1,J,K+1))*C2CXI(I)
       WFX=0.5*(F2FX(I-1)*W(I,J,K)+F2FX(I)*W(I-1,J,K))*C2CXI(I)

       ANXE=0.5*(U(I+1,J,K)+U(I,J,K))
       ANXW=0.5*(U(I,J,K)+U(I-1,J,K))
       ANXX=(ANXE**2-ANXW**2)*C2CXI(I)
       ANX1=(UN*VNX-US*VSX)*F2FYI(J)
       ANX2=(UC*WCX-UF*WFX)*F2FZI(K)
       RK3X=-ANXX-ANX1-ANX2               ! Non-linear term at k-substep

       ALX1=(U(I+1,J,K)-U(I,J,K))*F2FXI(I)
       ALX2=(U(I,J,K)-U(I-1,J,K))*F2FXI(I-1)
       ALX3=(U(I,J+1,K)-U(I,J,K))*C2CYI(J+1)*(1.-FIXJU(J)*FLOAT(JUT))
       ALX4=(U(I,J,K)-U(I,J-1,K))*C2CYI(J)  *(1.-FIXJL(J)*FLOAT(JUB))
       ALX5=(U(I,J,K+1)-U(I,J,K))*C2CZI(K+1)*(1.-FIXKU(K)*FLOAT(KUT))
       ALX6=(U(I,J,K)-U(I,J,K-1))*C2CZI(K)  *(1.-FIXKL(K)*FLOAT(KUB))
       ALXX=(ALX1-ALX2)*C2CXI(I)
       ALXY=(ALX3-ALX4)*F2FYI(J)
       ALXZ=(ALX5-ALX6)*F2FZI(K)
       ALX=1./RE*(ALXX+ALXY+ALXZ)           ! Linear terms at k-substep

!-----LES
      IF (ILES.EQ.1) THEN
       TALXX=(NUSGS(I,J,K)*ALX1-NUSGS(I-1,J,K)*ALX2)*C2CXI(I)
       TALXY=(NUSGS1(I,J+1,K,3)*ALX3-NUSGS1(I,J,K,3)*ALX4)*F2FYI(J)
       TALXZ=(NUSGS1(I,J,K+1,2)*ALX5-NUSGS1(I,J,K,2)*ALX6)*F2FZI(K)
       ALX=ALX+FLOAT(ILES)*(2.*TALXX+TALXY+TALXZ)
       RK3X=RK3X+FLOAT(ILES)*RHS1(I,J,K,1)  ! RHS1 in here: RHSSGS(lica_sgs.f90)
      ENDIF
!-----LES

!-----HEAT TRANSFER
      IF ((IHTRANS.EQ.1).AND.(GRDIR.EQ.1)) THEN
        RK3X = RK3X + GR/RE**2.*T(I,J,K)
      ENDIF
!-----HEAT TRANSFER

       RHS1(I,J,K,1)=DT                                    &
               *( GAMMA(MSUB)*RK3X+RO(MSUB)*RK3XO(I,J,K)   &
               +2.*ALPHA*ALX                               &
               -2.*ALPHA*(P(I,J,K)-P(I-1,J,K))*C2CXI(I)    &
               -2.*ALPHA*PMI(0)*FLOAT(ICH) )
       RK3XOO(I,J,K)=RK3XO(I,J,K)
       RK3XO(I,J,K)=RK3X

       ENDDO
       ENDDO
       ENDDO

!-----RHS1 calculation for v momentum -----------------
!$OMP PARALLEL DO  &
!$OMP private(VE,VW,VC,VF,UEY,UWY,WCY,WFY)   &
!$OMP private(ANYN,ANYS,ANYY,ANY1,ANY2,RK3Y) &
!$OMP private(ALY1,ALY2,ALY3,ALY4,ALY5,ALY6) &
!$OMP private(ALYX,ALYY,ALYZ,ALY,TALYX,TALYY,TALYZ)
       DO K=1,N3M
       DO J=J_BGPY,N2M
       DO I=1,N1M

       VE=(0.5*(F2FX(I+1)*V(I,J,K)+F2FX(I)*V(I+1,J,K))*C2CXI(I+1))   &
          *(1.-FIXIU(I))+V(N1,J,K)*FIXIU(I)
       VW=(0.5*(F2FX(I)*V(I-1,J,K)+F2FX(I-1)*V(I,J,K))*C2CXI(I))     &
          *(1.-FIXIL(I))+V(0,J,K)*FIXIL(I)
       VC=(0.5*(F2FZ(K+1)*V(I,J,K)+F2FZ(K)*V(I,J,K+1))*C2CZI(K+1))   &
          *(1.-FIXKU(K))+V(I,J,N3)*FIXKU(K)
       VF=(0.5*(F2FZ(K)*V(I,J,K-1)+F2FZ(K-1)*V(I,J,K))*C2CZI(K))     &
          *(1.-FIXKL(K))+V(I,J,0)*FIXKL(K)
       UEY=0.5*(F2FY(J-1)*U(I+1,J,K)+F2FY(J)*U(I+1,J-1,K))*C2CYI(J)
       UWY=0.5*(F2FY(J-1)*U(I,J,K)+F2FY(J)*U(I,J-1,K))*C2CYI(J)
       WCY=0.5*(F2FY(J-1)*W(I,J,K+1)+F2FY(J)*W(I,J-1,K+1))*C2CYI(J)
       WFY=0.5*(F2FY(J-1)*W(I,J,K)+F2FY(J)*W(I,J-1,K))*C2CYI(J)

       ANYN=0.5*(V(I,J+1,K)+V(I,J,K))
       ANYS=0.5*(V(I,J,K)+V(I,J-1,K))
       ANYY=C2CYI(J)*(ANYN**2-ANYS**2)
       ANY1=(UEY*VE-UWY*VW)*F2FXI(I)
       ANY2=(VC*WCY-VF*WFY)*F2FZI(K)
       RK3Y=-ANYY-ANY1-ANY2                ! Ny term

       ALY1=(V(I+1,J,K)-V(I,J,K))*C2CXI(I+1)
       ALY2=(V(I,J,K)-V(I-1,J,K))*C2CXI(I)
       ALY3=(V(I,J+1,K)-V(I,J,K))*F2FYI(J)
       ALY4=(V(I,J,K)-V(I,J-1,K))*F2FYI(J-1)
       ALY5=(V(I,J,K+1)-V(I,J,K))*C2CZI(K+1)*(1.-FIXKU(K)*FLOAT(KVT))
       ALY6=(V(I,J,K)-V(I,J,K-1))*C2CZI(K) *(1.-FIXKL(K)*FLOAT(KVB))
       ALYX=(ALY1-ALY2)*F2FXI(I)
       ALYY=(ALY3-ALY4)*C2CYI(J)
       ALYZ=(ALY5-ALY6)*F2FZI(K)
       ALY=1./RE*(ALYX+ALYY+ALYZ)

!-----LES
      IF (ILES.EQ.1) THEN
       TALYX=(NUSGS1(I+1,J,K,3)*ALY1-NUSGS1(I,J,K,3)*ALY2)*F2FXI(I)
       TALYY=(NUSGS(I,J,K)*ALY3-NUSGS(I,J-1,K)*ALY4)*C2CYI(J)
       TALYZ=(NUSGS1(I,J,K+1,1)*ALY5-NUSGS1(I,J,K,1)*ALY6)*F2FZI(K)
       ALY=ALY+FLOAT(ILES)*(TALYX+2.*TALYY+TALYZ)
       RK3Y=RK3Y+FLOAT(ILES)*RHS1(I,J,K,2)
      ENDIF
!-----LES

!-----HEAT TRANSFER
      IF ((IHTRANS.EQ.1).AND.(GRDIR.EQ.2)) THEN
        RK3Y = RK3Y + GR/RE**2.*T(I,J,K)
      ENDIF
!-----HEAT TRANSFER

       RHS1(I,J,K,2)=DT                                    &
              *( GAMMA(MSUB)*RK3Y+RO(MSUB)*RK3YO(I,J,K)    &
                +2.*ALPHA*ALY                              &
                -2.*ALPHA*(P(I,J,K)-P(I,J-1,K))*C2CYI(J))
       RK3YOO(I,J,K)=RK3YO(I,J,K)  ! added term for Fy
       RK3YO(I,J,K)=RK3Y

       ENDDO
       ENDDO
       ENDDO

       IF (N3M.EQ.1) GOTO 100
!-----RHS1 calculation for w momentum -----------------
!$OMP PARALLEL DO &
!$OMP private(WE,WW,WN,WS,UEZ,UWZ,VNZ,VSZ)    &
!$OMP private(ANZC,ANZF,ANZZ,ANZ1,ANZ2,RK3Z)  &
!$OMP private(ALZ1,ALZ2,ALZ3,ALZ4,ALZ5,ALZ6)  &
!$OMP private(ALZX,ALZY,ALZZ,ALZ,TALZX,TALZY,TALZZ)
       DO K=K_BGPZ,N3M
       DO J=1,N2M
       DO I=1,N1M

       WE=0.5*(F2FX(I+1)*W(I,J,K)+F2FX(I)*W(I+1,J,K))*C2CXI(I+1)    &
           *(1.-FIXIU(I))+FIXIU(I)*W(N1,J,K)
       WW=0.5*(F2FX(I)*W(I-1,J,K)+F2FX(I-1)*W(I,J,K))*C2CXI(I)      &
           *(1.-FIXIL(I))+FIXIL(I)*W(0,J,K)
       WN=(0.5*(F2FY(J+1)*W(I,J,K)+F2FY(J)*W(I,J+1,K))*C2CYI(J+1))  &
           *(1.-FIXJU(J))+FIXJU(J)*W(I,N2,K)
       WS=(0.5*(F2FY(J)*W(I,J-1,K)+F2FY(J-1)*W(I,J,K))*C2CYI(J))    &
           *(1.-FIXJL(J))+FIXJL(J)*W(I,0,K)
       UEZ=0.5*(F2FZ(K-1)*U(I+1,J,K)+F2FZ(K)*U(I+1,J,K-1))*C2CZI(K)
       UWZ=0.5*(F2FZ(K-1)*U(I,J,K)+F2FZ(K)*U(I,J,K-1))*C2CZI(K)
       VNZ=0.5*(F2FZ(K-1)*V(I,J+1,K)+F2FZ(K)*V(I,J+1,K-1))*C2CZI(K)
       VSZ=0.5*(F2FZ(K-1)*V(I,J,K)+F2FZ(K)*V(I,J,K-1))*C2CZI(K)

       ANZC=0.5*(W(I,J,K+1)+W(I,J,K))
       ANZF=0.5*(W(I,J,K)+W(I,J,K-1))
       ANZZ=(ANZC**2-ANZF**2)*C2CZI(K)
       ANZ1=(UEZ*WE-UWZ*WW)*F2FXI(I)
       ANZ2=(VNZ*WN-VSZ*WS)*F2FYI(J)
       RK3Z=-ANZZ-ANZ1-ANZ2                        ! Nz term

       ALZ1=(W(I+1,J,K)-W(I,J,K))*C2CXI(I+1)
       ALZ2=(W(I,J,K)-W(I-1,J,K))*C2CXI(I)
       ALZ3=(W(I,J+1,K)-W(I,J,K))*C2CYI(J+1)*(1.-FIXJU(J)*FLOAT(JWT))
       ALZ4=(W(I,J,K)-W(I,J-1,K))*C2CYI(J)  *(1.-FIXJL(J)*FLOAT(JWB))
       ALZ5=(W(I,J,K+1)-W(I,J,K))*F2FZI(K)
       ALZ6=(W(I,J,K)-W(I,J,K-1))*F2FZI(K-1)
       ALZX=(ALZ1-ALZ2)*F2FXI(I)
       ALZY=(ALZ3-ALZ4)*F2FYI(J)
       ALZZ=(ALZ5-ALZ6)*C2CZI(K)
       ALZ=1./RE*(ALZX+ALZY+ALZZ)

!-----LES
      IF (ILES.EQ.1) THEN
       TALZX=(NUSGS1(I+1,J,K,2)*ALZ1-NUSGS1(I,J,K,2)*ALZ2)*F2FXI(I)
       TALZY=(NUSGS1(I,J+1,K,1)*ALZ3-NUSGS1(I,J,K,1)*ALZ4)*F2FYI(J)
       TALZZ=(NUSGS(I,J,K)*ALZ5-NUSGS(I,J,K-1)*ALZ6)*C2CZI(K)
       ALZ=ALZ+FLOAT(ILES)*(TALZX+TALZY+2.*TALZZ)
       RK3Z=RK3Z+FLOAT(ILES)*RHS1(I,J,K,3)
      ENDIF
!-----LES

!-----HEAT TRANSFER
      IF ((IHTRANS.EQ.1).AND.(GRDIR.EQ.3)) THEN
        RK3Z = RK3Z + GR/RE**2.*T(I,J,K)
      ENDIF
!-----HEAT TRANSFER

       RHS1(I,J,K,3)=DT                                   &
              *( GAMMA(MSUB)*RK3Z+RO(MSUB)*RK3ZO(I,J,K)   &
                +2.*ALPHA*ALZ                             &
                -2.*ALPHA*(P(I,J,K)-P(I,J,K-1))*C2CZI(K))
       RK3ZOO(I,J,K)=RK3ZO(I,J,K)      ! added term for Fz
       RK3ZO(I,J,K)=RK3Z

       ENDDO
       ENDDO
       ENDDO

  100  RETURN
       END
!=======================================================================
      SUBROUTINE RHS_IBM
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
      USE MOD_FLOWARRAY, ONLY : RHS1,U,V,W,IFC,JFC,KFC,INTPINDX,GEOMFAC,&
                                FCV,FCVAVG,DUDTR
      IMPLICIT NONE
      INTEGER*8 :: I,J,K,N,L,II,JJ,KK,IP,JP,KP,IPP,JPP,KPP
      REAL*8    :: UBODY,VBODY,WBODY,DTI
      REAL*8    :: UTARG,VTARG,WTARG
      REAL*8    :: RHSTMP(N1M,N2M,N3M,3)
      REAL*8    :: UBD,VBD,WBD

!---- Compute target velocities & forcing values at forcing points
      UBODY = 0.
      VBODY = 0.
      WBODY = 0.
      DTI   = 1./DT

!$OMP PARALLEL DO
      DO K= 1,N3M
      DO J= 1,N2M
      DO I= 1,N1M
      DO L= 1,3
        RHSTMP(I,J,K,L)=RHS1(I,J,K,L)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

!$OMP PARALLEL DO &
!$OMP private(II,JJ,KK,IP,JP,KP,IPP,JPP,KPP,UTARG,UBODY)
      DO N=1,NINTP(1)
         II =IFC(N,1)
         JJ =JFC(N,1)
         KK =KFC(N,1)
         IP =IFC(N,1)+INTPINDX(N,1,1)
         JP =JFC(N,1)+INTPINDX(N,1,2)
         KP =KFC(N,1)+INTPINDX(N,1,3)
         IPP=IFC(N,1)+INTPINDX(N,1,1)*2
         JPP=JFC(N,1)+INTPINDX(N,1,2)*2
         KPP=KFC(N,1)+INTPINDX(N,1,3)*2

         IF ((IPP.GE.N1).OR.(IPP.LE.0)) CALL REINDEX_I(IP,IPP)
         IF ((KPP.GE.N3).OR.(KPP.LE.0)) CALL REINDEX_K(KP,KPP)
         ! IF (IMOVINGON.EQ.1) UBODY=UBD(X(II),YMP(JJ),ZMP(KK))

         UTARG=GEOMFAC(N,1,0,0,0)* UBODY                                  &
              +GEOMFAC(N,1,0,0,1)*(U(II ,JJ ,KP )+RHSTMP(II ,JJ ,KP ,1))  &
              +GEOMFAC(N,1,0,0,2)*(U(II ,JJ ,KPP)+RHSTMP(II ,JJ ,KPP,1))  &
              +GEOMFAC(N,1,0,1,0)*(U(II ,JP ,KK )+RHSTMP(II ,JP ,KK ,1))  &
              +GEOMFAC(N,1,0,1,1)*(U(II ,JP ,KP )+RHSTMP(II ,JP ,KP ,1))  &
              +GEOMFAC(N,1,0,1,2)*(U(II ,JP ,KPP)+RHSTMP(II ,JP ,KPP,1))  &
              +GEOMFAC(N,1,0,2,0)*(U(II ,JPP,KK )+RHSTMP(II ,JPP,KK ,1))  &
              +GEOMFAC(N,1,0,2,1)*(U(II ,JPP,KP )+RHSTMP(II ,JPP,KP ,1))  &
              +GEOMFAC(N,1,0,2,2)*(U(II ,JPP,KPP)+RHSTMP(II ,JPP,KPP,1))  &
              +GEOMFAC(N,1,1,0,0)*(U(IP ,JJ ,KK )+RHSTMP(IP ,JJ ,KK ,1))  &
              +GEOMFAC(N,1,1,0,1)*(U(IP ,JJ ,KP )+RHSTMP(IP ,JJ ,KP ,1))  &
              +GEOMFAC(N,1,1,0,2)*(U(IP ,JJ ,KPP)+RHSTMP(IP ,JJ ,KPP,1))  &
              +GEOMFAC(N,1,1,1,0)*(U(IP ,JP ,KK )+RHSTMP(IP ,JP ,KK ,1))  &
              +GEOMFAC(N,1,1,1,1)*(U(IP ,JP ,KP )+RHSTMP(IP ,JP ,KP ,1))  &
              +GEOMFAC(N,1,1,1,2)*(U(IP ,JP ,KPP)+RHSTMP(IP ,JP ,KPP,1))  &
              +GEOMFAC(N,1,1,2,0)*(U(IP ,JPP,KK )+RHSTMP(IP ,JPP,KK ,1))  &
              +GEOMFAC(N,1,1,2,1)*(U(IP ,JPP,KP )+RHSTMP(IP ,JPP,KP ,1))  &
              +GEOMFAC(N,1,1,2,2)*(U(IP ,JPP,KPP)+RHSTMP(IP ,JPP,KPP,1))  &
              +GEOMFAC(N,1,2,0,0)*(U(IPP,JJ ,KK )+RHSTMP(IPP,JJ ,KK ,1))  &
              +GEOMFAC(N,1,2,0,1)*(U(IPP,JJ ,KP )+RHSTMP(IPP,JJ ,KP ,1))  &
              +GEOMFAC(N,1,2,0,2)*(U(IPP,JJ ,KPP)+RHSTMP(IPP,JJ ,KPP,1))  &
              +GEOMFAC(N,1,2,1,0)*(U(IPP,JP ,KK )+RHSTMP(IPP,JP ,KK ,1))  &
              +GEOMFAC(N,1,2,1,1)*(U(IPP,JP ,KP )+RHSTMP(IPP,JP ,KP ,1))  &
              +GEOMFAC(N,1,2,1,2)*(U(IPP,JP ,KPP)+RHSTMP(IPP,JP ,KPP,1))  &
              +GEOMFAC(N,1,2,2,0)*(U(IPP,JPP,KK )+RHSTMP(IPP,JPP,KK ,1))  &
              +GEOMFAC(N,1,2,2,1)*(U(IPP,JPP,KP )+RHSTMP(IPP,JPP,KP ,1))  &
              +GEOMFAC(N,1,2,2,2)*(U(IPP,JPP,KPP)+RHSTMP(IPP,JPP,KPP,1))
         FCV(N,1)=(UTARG-U(II,JJ,KK)-RHSTMP(II,JJ,KK,1))*DTI
         DUDTR(N,1)=(UTARG-U(II,JJ,KK))*DTI
         RHS1(II,JJ,KK,1)=UTARG-U(II,JJ,KK)
      ENDDO

!$OMP PARALLEL DO &
!$OMP private(II,JJ,KK,IP,JP,KP,IPP,JPP,KPP,VTARG,VBODY)
      DO N=1,NINTP(2)
         II =IFC(N,2)
         JJ =JFC(N,2)
         KK =KFC(N,2)
         IP =IFC(N,2)+INTPINDX(N,2,1)
         JP =JFC(N,2)+INTPINDX(N,2,2)
         KP =KFC(N,2)+INTPINDX(N,2,3)
         IPP=IFC(N,2)+INTPINDX(N,2,1)*2
         JPP=JFC(N,2)+INTPINDX(N,2,2)*2
         KPP=KFC(N,2)+INTPINDX(N,2,3)*2

         IF ((IPP.GE.N1).OR.(IPP.LE.0)) CALL REINDEX_I(IP,IPP)
         IF ((KPP.GE.N3).OR.(KPP.LE.0)) CALL REINDEX_K(KP,KPP)
         ! IF (IMOVINGON.EQ.1) VBODY=VBD(XMP(II),Y(JJ),ZMP(KK))

         VTARG=GEOMFAC(N,2,0,0,0)* VBODY                                  &
              +GEOMFAC(N,2,0,0,1)*(V(II ,JJ ,KP )+RHSTMP(II ,JJ ,KP ,2))  &
              +GEOMFAC(N,2,0,0,2)*(V(II ,JJ ,KPP)+RHSTMP(II ,JJ ,KPP,2))  &
              +GEOMFAC(N,2,0,1,0)*(V(II ,JP ,KK )+RHSTMP(II ,JP ,KK ,2))  &
              +GEOMFAC(N,2,0,1,1)*(V(II ,JP ,KP )+RHSTMP(II ,JP ,KP ,2))  &
              +GEOMFAC(N,2,0,1,2)*(V(II ,JP ,KPP)+RHSTMP(II ,JP ,KPP,2))  &
              +GEOMFAC(N,2,0,2,0)*(V(II ,JPP,KK )+RHSTMP(II ,JPP,KK ,2))  &
              +GEOMFAC(N,2,0,2,1)*(V(II ,JPP,KP )+RHSTMP(II ,JPP,KP ,2))  &
              +GEOMFAC(N,2,0,2,2)*(V(II ,JPP,KPP)+RHSTMP(II ,JPP,KPP,2))  &
              +GEOMFAC(N,2,1,0,0)*(V(IP ,JJ ,KK )+RHSTMP(IP ,JJ ,KK ,2))  &
              +GEOMFAC(N,2,1,0,1)*(V(IP ,JJ ,KP )+RHSTMP(IP ,JJ ,KP ,2))  &
              +GEOMFAC(N,2,1,0,2)*(V(IP ,JJ ,KPP)+RHSTMP(IP ,JJ ,KPP,2))  &
              +GEOMFAC(N,2,1,1,0)*(V(IP ,JP ,KK )+RHSTMP(IP ,JP ,KK ,2))  &
              +GEOMFAC(N,2,1,1,1)*(V(IP ,JP ,KP )+RHSTMP(IP ,JP ,KP ,2))  &
              +GEOMFAC(N,2,1,1,2)*(V(IP ,JP ,KPP)+RHSTMP(IP ,JP ,KPP,2))  &
              +GEOMFAC(N,2,1,2,0)*(V(IP ,JPP,KK )+RHSTMP(IP ,JPP,KK ,2))  &
              +GEOMFAC(N,2,1,2,1)*(V(IP ,JPP,KP )+RHSTMP(IP ,JPP,KP ,2))  &
              +GEOMFAC(N,2,1,2,2)*(V(IP ,JPP,KPP)+RHSTMP(IP ,JPP,KPP,2))  &
              +GEOMFAC(N,2,2,0,0)*(V(IPP,JJ ,KK )+RHSTMP(IPP,JJ ,KK ,2))  &
              +GEOMFAC(N,2,2,0,1)*(V(IPP,JJ ,KP )+RHSTMP(IPP,JJ ,KP ,2))  &
              +GEOMFAC(N,2,2,0,2)*(V(IPP,JJ ,KPP)+RHSTMP(IPP,JJ ,KPP,2))  &
              +GEOMFAC(N,2,2,1,0)*(V(IPP,JP ,KK )+RHSTMP(IPP,JP ,KK ,2))  &
              +GEOMFAC(N,2,2,1,1)*(V(IPP,JP ,KP )+RHSTMP(IPP,JP ,KP ,2))  &
              +GEOMFAC(N,2,2,1,2)*(V(IPP,JP ,KPP)+RHSTMP(IPP,JP ,KPP,2))  &
              +GEOMFAC(N,2,2,2,0)*(V(IPP,JPP,KK )+RHSTMP(IPP,JPP,KK ,2))  &
              +GEOMFAC(N,2,2,2,1)*(V(IPP,JPP,KP )+RHSTMP(IPP,JPP,KP ,2))  &
              +GEOMFAC(N,2,2,2,2)*(V(IPP,JPP,KPP)+RHSTMP(IPP,JPP,KPP,2))
         FCV(N,2)=(VTARG-V(II,JJ,KK)-RHSTMP(II,JJ,KK,2))*DTI
         DUDTR(N,2)=(VTARG-V(II,JJ,KK))*DTI
         RHS1(II,JJ,KK,2)=VTARG-V(II,JJ,KK)
      ENDDO

!$OMP PARALLEL DO &
!$OMP private(II,JJ,KK,IP,JP,KP,IPP,JPP,KPP,WTARG,WBODY)
      DO N=1,NINTP(3)
         II =IFC(N,3)
         JJ =JFC(N,3)
         KK =KFC(N,3)
         IP =IFC(N,3)+INTPINDX(N,3,1)
         JP =JFC(N,3)+INTPINDX(N,3,2)
         KP =KFC(N,3)+INTPINDX(N,3,3)
         IPP=IFC(N,3)+INTPINDX(N,3,1)*2
         JPP=JFC(N,3)+INTPINDX(N,3,2)*2
         KPP=KFC(N,3)+INTPINDX(N,3,3)*2

         IF ((IPP.GE.N1).OR.(IPP.LE.0)) CALL REINDEX_I(IP,IPP)
         IF ((KPP.GE.N3).OR.(KPP.LE.0)) CALL REINDEX_K(KP,KPP)
         ! IF (IMOVINGON.EQ.1) WBODY=WBD(XMP(II),YMP(JJ),Z(KK))

         WTARG=GEOMFAC(N,3,0,0,0)* WBODY                                  &
              +GEOMFAC(N,3,0,0,1)*(W(II ,JJ ,KP )+RHSTMP(II ,JJ ,KP ,3))  &
              +GEOMFAC(N,3,0,0,2)*(W(II ,JJ ,KPP)+RHSTMP(II ,JJ ,KPP,3))  &
              +GEOMFAC(N,3,0,1,0)*(W(II ,JP ,KK )+RHSTMP(II ,JP ,KK ,3))  &
              +GEOMFAC(N,3,0,1,1)*(W(II ,JP ,KP )+RHSTMP(II ,JP ,KP ,3))  &
              +GEOMFAC(N,3,0,1,2)*(W(II ,JP ,KPP)+RHSTMP(II ,JP ,KPP,3))  &
              +GEOMFAC(N,3,0,2,0)*(W(II ,JPP,KK )+RHSTMP(II ,JPP,KK ,3))  &
              +GEOMFAC(N,3,0,2,1)*(W(II ,JPP,KP )+RHSTMP(II ,JPP,KP ,3))  &
              +GEOMFAC(N,3,0,2,2)*(W(II ,JPP,KPP)+RHSTMP(II ,JPP,KPP,3))  &
              +GEOMFAC(N,3,1,0,0)*(W(IP ,JJ ,KK )+RHSTMP(IP ,JJ ,KK ,3))  &
              +GEOMFAC(N,3,1,0,1)*(W(IP ,JJ ,KP )+RHSTMP(IP ,JJ ,KP ,3))  &
              +GEOMFAC(N,3,1,0,2)*(W(IP ,JJ ,KPP)+RHSTMP(IP ,JJ ,KPP,3))  &
              +GEOMFAC(N,3,1,1,0)*(W(IP ,JP ,KK )+RHSTMP(IP ,JP ,KK ,3))  &
              +GEOMFAC(N,3,1,1,1)*(W(IP ,JP ,KP )+RHSTMP(IP ,JP ,KP ,3))  &
              +GEOMFAC(N,3,1,1,2)*(W(IP ,JP ,KPP)+RHSTMP(IP ,JP ,KPP,3))  &
              +GEOMFAC(N,3,1,2,0)*(W(IP ,JPP,KK )+RHSTMP(IP ,JPP,KK ,3))  &
              +GEOMFAC(N,3,1,2,1)*(W(IP ,JPP,KP )+RHSTMP(IP ,JPP,KP ,3))  &
              +GEOMFAC(N,3,1,2,2)*(W(IP ,JPP,KPP)+RHSTMP(IP ,JPP,KPP,3))  &
              +GEOMFAC(N,3,2,0,0)*(W(IPP,JJ ,KK )+RHSTMP(IPP,JJ ,KK ,3))  &
              +GEOMFAC(N,3,2,0,1)*(W(IPP,JJ ,KP )+RHSTMP(IPP,JJ ,KP ,3))  &
              +GEOMFAC(N,3,2,0,2)*(W(IPP,JJ ,KPP)+RHSTMP(IPP,JJ ,KPP,3))  &
              +GEOMFAC(N,3,2,1,0)*(W(IPP,JP ,KK )+RHSTMP(IPP,JP ,KK ,3))  &
              +GEOMFAC(N,3,2,1,1)*(W(IPP,JP ,KP )+RHSTMP(IPP,JP ,KP ,3))  &
              +GEOMFAC(N,3,2,1,2)*(W(IPP,JP ,KPP)+RHSTMP(IPP,JP ,KPP,3))  &
              +GEOMFAC(N,3,2,2,0)*(W(IPP,JPP,KK )+RHSTMP(IPP,JPP,KK ,3))  &
              +GEOMFAC(N,3,2,2,1)*(W(IPP,JPP,KP )+RHSTMP(IPP,JPP,KP ,3))  &
              +GEOMFAC(N,3,2,2,2)*(W(IPP,JPP,KPP)+RHSTMP(IPP,JPP,KPP,3))
         FCV(N,3)=(WTARG-W(II,JJ,KK)-RHSTMP(II,JJ,KK,3))*DTI
         DUDTR(N,3)=(WTARG-W(II,JJ,KK))*DTI
         RHS1(II,JJ,KK,3)=WTARG-W(II,JJ,KK)
      ENDDO

!*************************       CAUTION       *************************
!     FOLLOWING 3 DO-LOOPS REQUIRES THE OPTION OF
!     '-Wf "-pvctl vwork=stack"'
!     WHEN THE CODE IS COINTPINDXLED ON NEC MACHINE.
!     FOR DETAILED EXPLANATION OF THE REASON,
!     YOU CAN CONSULT "NEC PORTING GUIDE"
!     PROVIDED BY NEC SUPERCOMPUTIONG CENTER.
!***********************************************************************

!$OMP PARALLEL DO &
!$OMP private(II,JJ,KK,UTARG,UBODY)
      DO N=NINTP(1)+1,NBODY(1)
         II=IFC(N,1)
         JJ=JFC(N,1)
         KK=KFC(N,1)
         ! IF (IMOVINGON.EQ.1) UBODY=UBD(X(II),YMP(JJ),ZMP(KK))
         UTARG=UBODY
         FCV(N,1)=(UTARG-U(II,JJ,KK)-RHS1(II,JJ,KK,1))*DTI
         DUDTR(N,1)=(UTARG-U(II,JJ,KK))*DTI
         RHS1(II,JJ,KK,1)=UTARG-U(II,JJ,KK)
      ENDDO

!$OMP PARALLEL DO &
!$OMP private(II,JJ,KK,VTARG,VBODY)
      DO N=NINTP(2)+1,NBODY(2)
         II=IFC(N,2)
         JJ=JFC(N,2)
         KK=KFC(N,2)
         ! IF (IMOVINGON.EQ.1) VBODY=VBD(XMP(II),Y(JJ),ZMP(KK))
         VTARG=VBODY
         FCV(N,2)=(VTARG-V(II,JJ,KK)-RHS1(II,JJ,KK,2))*DTI
         DUDTR(N,2)=(VTARG-V(II,JJ,KK))*DTI
         RHS1(II,JJ,KK,2)=VTARG-V(II,JJ,KK)
      ENDDO

!$OMP PARALLEL DO &
!$OMP private(II,JJ,KK,WTARG,WBODY)
      DO N=NINTP(3)+1,NBODY(3)
         II=IFC(N,3)
         JJ=JFC(N,3)
         KK=KFC(N,3)
         ! IF (IMOVINGON.EQ.1) WBODY=WBD(XMP(II),YMP(JJ),Z(KK))
         WTARG=WBODY
         FCV(N,3)=(WTARG-W(II,JJ,KK)-RHS1(II,JJ,KK,3))*DTI
         DUDTR(N,3)=(WTARG-W(II,JJ,KK))*DTI
         RHS1(II,JJ,KK,3)=WTARG-W(II,JJ,KK)
      ENDDO

!-----compute average forcing values during RK3 steps
!$OMP PARALLEL DO
      DO L=1,3
      DO N=1,NBODY(L)
         FCVAVG(N,L)=FCVAVG(N,L)+FCV(N,L)
      ENDDO
      ENDDO

      RETURN
      END
!=======================================================================
      SUBROUTINE REINDEX_I(IP,IPP)
!=======================================================================
!
!     Assign appropriate indices to the variables used in IBM interpolation
!     in order to satisfy the periodic BC in x-direction.
!      IP =I +1
!      IPP=IP+1
!     *Generally, (IPP.GE.N1)~ situations occur in case of periodic BC
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8     :: IP,IPP

      IF (IPP.EQ.N1) THEN
       IPP = 1
      ELSEIF (IPP.EQ.N1+1) THEN
       IP  = 1
       IPP = 2
      ELSEIF (IPP.EQ.0) THEN
       IPP = N1M
      ELSEIF (IPP.EQ.-1) THEN
       IP  = N1M
       IPP = N1M-1
      ENDIF

      RETURN
      END SUBROUTINE REINDEX_I
!=======================================================================
!=======================================================================
      SUBROUTINE REINDEX_K(KP,KPP)
!=======================================================================
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8     :: KP,KPP

      IF (KPP.EQ.N3) THEN
       KPP = 1
      ELSEIF (KPP.EQ.N3+1) THEN
       KP  = 1
       KPP = 2
      ELSEIF (KPP.EQ.0) THEN
       KPP = N3M
      ELSEIF (KPP.EQ.-1) THEN
       KP  = N3M
       KPP = N3M-1
      ENDIF

      RETURN
      END SUBROUTINE REINDEX_K
!=======================================================================
!=======================================================================
        SUBROUTINE CONVBC
!=======================================================================
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W
      IMPLICIT NONE
      INTEGER*8   :: I,J,K
      REAL*8      :: QIN,QOUT,QRATIO
      REAL*8      :: UBAR,UCOEF,VWCOEF
      REAL*8      :: FUNCBODY

      IF (ICH .NE. 1) THEN

!-----CALCULATE INFLUX
      QIN=0.
!$OMP PARALLEL DO REDUCTION(+:QIN)
      DO K=1,N3M
        DO J=1,N2M
          QIN=U(1,J,K)*F2FY(J)*F2FZ(K)+QIN
        ENDDO
      ENDDO
!$OMP PARALLEL DO REDUCTION(+:QIN)
      DO K=1,N3M
        DO I=1,N1M
          QIN=(V(I,1,K)-V(I,N2,K))*F2FX(I)*F2FZ(K)+QIN
        ENDDO
      ENDDO
!$OMP PARALLEL DO REDUCTION(+:QIN)
      DO J=1,N2M
        DO I=1,N1M
         QIN=(W(I,J,1)-W(I,J,N3))*F2FX(I)*F2FY(J)+QIN
        ENDDO
      ENDDO

!-----CALCULATE CONVECTIVE VELOCITY Uc
      UBAR = 0.
!$OMP PARALLEL DO REDUCTION(+:UBAR)
      DO K=1,N3M
        DO J=1,N2M
          UBAR=U(N1,J,K)*F2FY(J)*F2FZ(K)+UBAR
        ENDDO
      ENDDO
      UBAR  = UBAR/(YL*ZL)
      UCOEF = UBAR*DTCONST*F2FXI(N1M)
      VWCOEF= UBAR*DTCONST*C2CXI(N1)

      QOUT = 0.
!$OMP PARALLEL DO REDUCTION(+:QOUT)
      DO K=1,N3M
        DO J=1,N2M
          IF (FUNCBODY(XMP(N1M),YMP(J),ZMP(K),TIME) .LE. 1.E-10 .AND. IMOVINGON .EQ. 0) THEN
             ! Solid Phase: Zero-gradient (preserves solid velocity = 0)
             UOUT(J,K) = U(N1M,J,K)
             VOUT(J,K) = V(N1M,J,K)
             WOUT(J,K) = W(N1M,J,K)
          ELSE
             ! Fluid Phase: Convective advection
             UOUT(J,K)=U(N1,J,K)-UCOEF*(U(N1,J,K)-U(N1M,J,K))
             VOUT(J,K)=V(N1,J,K)-VWCOEF*(V(N1,J,K)-V(N1M,J,K))
             WOUT(J,K)=W(N1,J,K)-VWCOEF*(W(N1,J,K)-W(N1M,J,K))
          ENDIF
          
          QOUT=UOUT(J,K)*F2FY(J)*F2FZ(K)+QOUT
        ENDDO
      ENDDO
      QRATIO=QIN/QOUT

!-----ADJUST BOUNDARY VELOCITY TO SATISFY GLOBAL MASS CONSERVATION
!$OMP PARALLEL DO
      DO K=1,N3M
        DO J=1,N2M
          ! Only scale the fluid phase by the mass flux error
          IF (FUNCBODY(XMP(N1M),YMP(J),ZMP(K),TIME) .GT. 1.E-10) THEN
             UOUT(J,K)=UOUT(J,K)*QRATIO
             VOUT(J,K)=VOUT(J,K)*QRATIO
             WOUT(J,K)=WOUT(J,K)*QRATIO
          ENDIF
          DUOUT(J,K)=UOUT(J,K)-U(N1,J,K)
          DVOUT(J,K)=VOUT(J,K)-V(N1,J,K)
          DWOUT(J,K)=WOUT(J,K)-W(N1,J,K)
        ENDDO
      ENDDO
   
      IF (YPRDIC .EQ. 1) THEN
!$OMP PARALLEL DO
        DO K=0,N3
          DUOUT(N2,K)=DUOUT(1,  K)
          DUOUT(0, K)=DUOUT(N2M,K)
          DVOUT(N2,K)=DVOUT(1,  K)
          DVOUT(0, K)=DVOUT(N2M,K)
          DWOUT(N2,K)=DWOUT(1,  K)
          DWOUT(0, K)=DWOUT(N2M,K)
        ENDDO
      ENDIF
     
      IF (ZPRDIC .EQ. 1) THEN
!$OMP PARALLEL DO
        DO J=0,N2
          DUOUT(J,N3)=DUOUT(J,  1)
          DUOUT(J, 0)=DUOUT(J,N3M)
          DVOUT(J,N3)=DVOUT(J,  1)
          DVOUT(J, 0)=DVOUT(J,N3M)
          DWOUT(J,N3)=DWOUT(J,  1)
          DWOUT(J, 0)=DWOUT(J,N3M)
        ENDDO
      ENDIF

      ELSE

      UOUT = 0.
      VOUT = 0.
      WOUT = 0.

      ENDIF

      RETURN
      END SUBROUTINE CONVBC
!=======================================================================
!=======================================================================
       SUBROUTINE RHSINCORPBC
!=======================================================================
!
!     Applying exit boundary condition & top boundary condition to RHS
!           [RHS incorporate B.C.]
!
!     ACOEF: main solver in this file
!     CIU, CIVW, CJV: subroutine LHSINIT in this file
!     DUOUT, DVOUT, DVTOP, DWOUT: subroutine CONVBC (lica_cylinder.f90)
!
!-----------------------------------------------------------------------
       USE MOD_COMMON
       USE MOD_FLOWARRAY, ONLY : U,V,W,RHS1,NUSGS,NUSGS1
       IMPLICIT NONE
       INTEGER*8     :: I,J,K
       REAL*8        :: CRE,CRE2

       IF (ILES .EQ. 1) THEN
         CRE=RE
         CRE2=2.*RE
       ELSE
         CRE=0.
         CRE2=0.
       ENDIF

!$OMP PARALLEL DO
       DO K=1,N3M
       DO J=1,N2M
       RHS1(N1M,J,K,1)=RHS1(N1M,J,K,1)                                  &
                      -ACOEF*CIU(N1M)*(1.+CRE2*NUSGS(N1M,J,K))*DUOUT(J,K)
       ENDDO
       ENDDO

!$OMP PARALLEL DO
       DO K=1,N3M
       DO J=2,N2M
       RHS1(N1M,J,K,2)=RHS1(N1M,J,K,2)                                  &
                      -ACOEF*CIVW(N1M)*(1.+CRE*NUSGS1(N1,J,K,3))*DVOUT(J,K)
       ENDDO
       ENDDO

!$OMP PARALLEL DO
       DO K=K_BGPZ,N3M
       DO J=1,N2M
       RHS1(N1M,J,K,3)=RHS1(N1M,J,K,3)                                  &
                      -ACOEF*CIVW(N1M)*(1.+CRE*NUSGS1(N1,J,K,2))*DWOUT(J,K)
       ENDDO
       ENDDO

       RETURN
       END SUBROUTINE RHSINCORPBC
!=======================================================================
!=======================================================================
      SUBROUTINE LHSU
!=======================================================================
!
!     Calculate intermediate velocity, u_i hat through TDMA
!           In this routine, compute streamwise velocity (u hat)
!
!-----------------------------------------------------------------------
       USE MOD_COMMON
       USE MOD_FLOWARRAY, ONLY : U,V,W,RHS1,NUSGS,NUSGS1
       IMPLICIT NONE
       INTEGER*8     :: I,J,K
       REAL*8        :: CRE,CRE2
       REAL*8, DIMENSION (:,:), ALLOCATABLE :: AI,BI,CI,GI
       REAL*8, DIMENSION (:,:), ALLOCATABLE :: AJ,BJ,CJ,GJ
       REAL*8, DIMENSION (:,:), ALLOCATABLE :: AK,BK,CK,GK

       CRE=FLOAT(ILES)*RE
       CRE2=2.*FLOAT(ILES)*RE

!=====ADI STARTS

      IF (N3M.EQ.1) GOTO 100
!-----Z-DIRECTION
!$OMP PARALLEL &
!$OMP private(AK,CK,BK,GK)
      allocate(AK(N1,N3),BK(N1,N3),CK(N1,N3),GK(N1,N3))
!$OMP DO
      DO 131 J=1,N2M
      DO 141 K=1,N3M
      DO 141 I=I_BGPX,N1M
       AK(I,K)=AKUV(K)*(1.+CRE*NUSGS1(I,J,K,2)) *(1.-FIXKL(K)*FLOAT(KUB))
       CK(I,K)=CKUV(K)*(1.+CRE*NUSGS1(I,J,K+1,2))*(1.-FIXKU(K)*FLOAT(KUT))
       BK(I,K)=ACOEFI-AK(I,K)-CK(I,K)
       GK(I,K)=ACOEFI*RHS1(I,J,K,1)
  141 CONTINUE

      IF (ZPRDIC .EQ. 0) THEN
         CALL TRDIAG3(AK,BK,CK,GK,GK,1,N3M,I_BGPX,N1M)
      ELSE IF (ZPRDIC .EQ. 1) THEN
         CALL TRDIAG3P(AK,BK,CK,GK,1,N3M,I_BGPX,N1M)  !z periodicity
      ENDIF

      DO 151 K=1,N3M
      DO 151 I=I_BGPX,N1M
      RHS1(I,J,K,1)=GK(I,K)
  151 CONTINUE
  131 CONTINUE
!$OMP END DO
      deallocate(AK,BK,CK,GK)
!$OMP END PARALLEL

  100 CONTINUE

!$OMP PARALLEL  &
!$OMP private(AJ,CJ,BJ,GJ)  &
!$OMP private(AI,CI,BI,GI)
      allocate(AI(N2,N1),BI(N2,N1),CI(N2,N1),GI(N2,N1))
      allocate(AJ(N1,N2),BJ(N1,N2),CJ(N1,N2),GJ(N1,N2))
!$OMP DO
      DO 40 K=1,N3M

!-----Y-DIRECTION
      DO 91 J=1,N2M
      DO 91 I=I_BGPX,N1M
      AJ(I,J)=AJUW(J)*(1.+CRE*NUSGS1(I,J,K,3))  *(1.-FIXJL(J)*FLOAT(JUB))
      CJ(I,J)=CJUW(J)*(1.+CRE*NUSGS1(I,J+1,K,3))*(1.-FIXJU(J)*FLOAT(JUT))
      BJ(I,J)=ACOEFI-AJ(I,J)-CJ(I,J)
      GJ(I,J)=ACOEFI*RHS1(I,J,K,1)
   91 CONTINUE

      CALL TRDIAG2(AJ,BJ,CJ,GJ,GJ,1,N2M,I_BGPX,N1M)

!-----X-DIRECTION
      DO 51 I=I_BGPX,N1M
      DO 51 J=1,N2M
      AI(J,I)=AIU(I)*(1.+CRE2*NUSGS(I-1,J,K))
      CI(J,I)=CIU(I)*(1.+CRE2*NUSGS(I,J,K))
      BI(J,I)=ACOEFI-AI(J,I)-CI(J,I)
      GI(J,I)=ACOEFI*GJ(I,J)
   51 CONTINUE

      IF (XPRDIC .EQ. 0) THEN
         CALL TRDIAG1(AI,BI,CI,GI,GI,I_BGPX,N1M,1,N2M)
      ELSE IF (XPRDIC .EQ. 1) THEN
         CALL TRDIAG1P(AI,BI,CI,GI,I_BGPX,N1M,1,N2M)
      ENDIF

      DO 61 J=1,N2M
      DO 61 I=I_BGPX,N1M
        U(I,J,K)=GI(J,I)+U(I,J,K)
   61 CONTINUE

   40 CONTINUE
!$OMP END DO
      deallocate(AI,BI,CI,GI)
      deallocate(AJ,BJ,CJ,GJ)
!$OMP END PARALLEL

      RETURN
      END

!=======================================================================
      SUBROUTINE LHSV
!=======================================================================
!
!     Calculate intermediate velocity, u_i hat through TDMA
!           In this routine, compute transverse velocity (v hat)
!
!-----------------------------------------------------------------------
       USE MOD_COMMON
       USE MOD_FLOWARRAY, ONLY : U,V,W,RHS1,NUSGS,NUSGS1
       IMPLICIT NONE
       INTEGER*8     :: I,J,K
       REAL*8        :: CRE,CRE2
       REAL*8, DIMENSION (:,:), ALLOCATABLE :: AI,BI,CI,GI
       REAL*8, DIMENSION (:,:), ALLOCATABLE :: AJ,BJ,CJ,GJ
       REAL*8, DIMENSION (:,:), ALLOCATABLE :: AK,BK,CK,GK

       CRE = FLOAT(ILES)*RE
       CRE2= 2.*FLOAT(ILES)*RE

!=====ADI STARTS

      IF (N3M.EQ.1) GOTO 100
!-----Z-DIRECTION
!$OMP PARALLEL &
!$OMP private(AK,CK,BK,GK)
      allocate(AK(N1,N3),BK(N1,N3),CK(N1,N3),GK(N1,N3))
!$OMP DO
      DO 131 J=2,N2M

      DO 141 K=1,N3M
      DO 141 I=1,N1M
      AK(I,K)=AKUV(K)*(1.+CRE*NUSGS1(I,J,K,1)) *(1.-FIXKL(K)*FLOAT(KVB))
      CK(I,K)=CKUV(K)*(1.+CRE*NUSGS1(I,J,K+1,1))*(1.-FIXKU(K)*FLOAT(KVT))
      BK(I,K)=ACOEFI-AK(I,K)-CK(I,K)
      GK(I,K)=ACOEFI*RHS1(I,J,K,2)
  141 CONTINUE

      IF (ZPRDIC .EQ. 0) THEN
         CALL TRDIAG3(AK,BK,CK,GK,GK,1,N3M,1,N1M)
      ELSE IF (ZPRDIC .EQ. 1) THEN
         CALL TRDIAG3P(AK,BK,CK,GK,1,N3M,1,N1M)
      ENDIF

      DO 151 K=1,N3M
      DO 151 I=1,N1M
      RHS1(I,J,K,2)=GK(I,K)
  151 CONTINUE
  131 CONTINUE
!$OMP END DO
      deallocate(AK,BK,CK,GK)
!$OMP END PARALLEL

  100 CONTINUE

!$OMP PARALLEL &
!$OMP private(AJ,CJ,BJ,GJ) &
!$OMP private(AI,CI,BI,GI)
      allocate(AI(N2,N1),BI(N2,N1),CI(N2,N1),GI(N2,N1))
      allocate(AJ(N1,N2),BJ(N1,N2),CJ(N1,N2),GJ(N1,N2))
!$OMP DO
      DO 40 K=1,N3M

!-----Y-DIRECTION
      DO 91 J=2,N2M
      DO 91 I=1,N1M
      AJ(I,J)=AJV(J)*(1.+CRE2*NUSGS(I,J-1,K))
      CJ(I,J)=CJV(J)*(1.+CRE2*NUSGS(I,J,K))
      BJ(I,J)=ACOEFI-AJ(I,J)-CJ(I,J)
      GJ(I,J)=ACOEFI*RHS1(I,J,K,2)
   91 CONTINUE

      CALL TRDIAG2(AJ,BJ,CJ,GJ,GJ,2,N2M,1,N1M)

!-----X-DIRECTION
      DO 51 I=1,N1M
      DO 51 J=2,N2M
      AI(J,I)=AIVW(I)*(1.+CRE*NUSGS1(I,J,K,3))
      CI(J,I)=CIVW(I)*(1.+CRE*NUSGS1(I+1,J,K,3))
      BI(J,I)=ACOEFI-AI(J,I)-CI(J,I)
      GI(J,I)=ACOEFI*GJ(I,J)
   51 CONTINUE

      IF (XPRDIC .EQ. 0) THEN
         CALL TRDIAG1(AI,BI,CI,GI,GI,1,N1M,2,N2M)
      ELSE IF (XPRDIC .EQ. 1) THEN
         CALL TRDIAG1P(AI,BI,CI,GI,1,N1M,2,N2M)
      ENDIF

      DO 61 J=2,N2M
      DO 61 I=1,N1M
      V(I,J,K)=GI(J,I)+V(I,J,K)
   61 CONTINUE

   40 CONTINUE
!$OMP END DO
      deallocate(AI,BI,CI,GI)
      deallocate(AJ,BJ,CJ,GJ)
!$OMP END PARALLEL

      RETURN
      END

!=======================================================================
      SUBROUTINE LHSW
!=======================================================================
!
!     Calculate intermediate velocity, u_i hat through TDMA
!           In this routine, compute spanwise velocity (w hat)
!
!-----------------------------------------------------------------------
       USE MOD_COMMON
       USE MOD_FLOWARRAY, ONLY : U,V,W,RHS1,NUSGS,NUSGS1
       IMPLICIT NONE
       INTEGER*8     :: I,J,K
       REAL*8        :: CRE,CRE2
       REAL*8, DIMENSION (:,:), ALLOCATABLE :: AI,BI,CI,GI
       REAL*8, DIMENSION (:,:), ALLOCATABLE :: AJ,BJ,CJ,GJ
       REAL*8, DIMENSION (:,:), ALLOCATABLE :: AK,BK,CK,GK

       CRE = FLOAT(ILES)*RE
       CRE2= 2.*FLOAT(ILES)*RE

!=====ADI STARTS

!-----Z-DIRECTION
!$OMP PARALLEL &
!$OMP private(AK,CK,BK,GK)
      allocate(AK(N1,N3),BK(N1,N3),CK(N1,N3),GK(N1,N3))
!$OMP DO
      DO 131 J=1,N2M

      DO 141 K=K_BGPZ,N3M
      DO 141 I=1,N1M
      AK(I,K)=AKW(K)*(1.+CRE2*NUSGS(I,J,K-1))
      CK(I,K)=CKW(K)*(1.+CRE2*NUSGS(I,J,K))
      BK(I,K)=ACOEFI-AK(I,K)-CK(I,K)
      GK(I,K)=ACOEFI*RHS1(I,J,K,3)
  141 CONTINUE

      IF (ZPRDIC .EQ. 0) THEN
         CALL TRDIAG3(AK,BK,CK,GK,GK,2,N3M,1,N1M)
      ELSE IF (ZPRDIC .EQ. 1) THEN
         CALL TRDIAG3P(AK,BK,CK,GK,1,N3M,1,N1M)
      ENDIF

      DO 151 K=K_BGPZ,N3M
      DO 151 I=1,N1M
      RHS1(I,J,K,3)=GK(I,K)
  151 CONTINUE

  131 CONTINUE
!$OMP END DO
      deallocate(AK,BK,CK,GK)
!$OMP END PARALLEL


!$OMP PARALLEL  &
!$OMP private(AJ,CJ,BJ,GJ)  &
!$OMP private(AI,CI,BI,GI)
      allocate(AI(N2,N1),BI(N2,N1),CI(N2,N1),GI(N2,N1))
      allocate(AJ(N1,N2),BJ(N1,N2),CJ(N1,N2),GJ(N1,N2))
!$OMP DO
      DO 40 K=K_BGPZ,N3M

!-----Y-DIRECTION
      DO 91 J=1,N2M
      DO 91 I=1,N1M
      AJ(I,J)=AJUW(J)*(1.+CRE*NUSGS1(I,J,K,1))  *(1.-FIXJL(J)*FLOAT(JWB))
      CJ(I,J)=CJUW(J)*(1.+CRE*NUSGS1(I,J+1,K,1))*(1.-FIXJU(J)*FLOAT(JWT))
      BJ(I,J)=ACOEFI-AJ(I,J)-CJ(I,J)
      GJ(I,J)=ACOEFI*RHS1(I,J,K,3)
   91 CONTINUE

      CALL TRDIAG2(AJ,BJ,CJ,GJ,GJ,1,N2M,1,N1M)

!-----X-DIRECTION
      DO 51 I=1,N1M
      DO 51 J=1,N2M
      AI(J,I)=AIVW(I)*(1.+CRE*NUSGS1(I,J,K,2))
      CI(J,I)=CIVW(I)*(1.+CRE*NUSGS1(I+1,J,K,2))
      BI(J,I)=ACOEFI-AI(J,I)-CI(J,I)
      GI(J,I)=ACOEFI*GJ(I,J)
   51 CONTINUE

      IF (XPRDIC .EQ. 0) THEN
         CALL TRDIAG1(AI,BI,CI,GI,GI,1,N1M,1,N2M)
      ELSE IF (XPRDIC .EQ. 1) THEN
         CALL TRDIAG1P(AI,BI,CI,GI,1,N1M,1,N2M)
      ENDIF

      DO 61 J=1,N2M
      DO 61 I=1,N1M
      W(I,J,K)=GI(J,I)+W(I,J,K)
   61 CONTINUE

   40 CONTINUE
!$OMP END DO
      deallocate(AI,BI,CI,GI)
      deallocate(AJ,BJ,CJ,GJ)
!$OMP END PARALLEL

      RETURN
      END

!=======================================================================
      SUBROUTINE TRDIAG1(A,B,C,R,UU,L1,L2,LL1,LL2)
!=======================================================================
!
!     SOLVE THE TRIDIAGONAL MATRIX (X-1 TRIDIAGONAL) WITH
!     DIFFERENT COEFFICIENTS in LHSU,LHSV,LHSW subroutines
!
!     Variables
!       A(),B(),C()  : Coefficients of tridiagonal matrix
!       R()          : RHS
!       UU()         : Solution
!       L1,L2  : 2nd index indicating x1-direction (TDMA)
!       LL2,LL2: 1st index indicating x2-direction
!
!       At each i,
!       |B_{L1}   C_{L1}   0        0               | = |R_{L1  }|
!       |A_{L1+1} B_{L1+1} C_{L1+1} 0               | = |R_{L1+1}|
!       |0        A_{L1+2} B_{L1+2} C_{L1+2}        | = |R_{L1+2}|
!       |0                                          | = |  :     |
!       |0                                          | = |  :     |
!       |0                             A_{L2} B_{L2}| = |R_{L2}  |
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8 :: I,J,L1,L2,LL1,LL2
      REAL*8    :: GAM(N2,N1),A(N2,N1),B(N2,N1),C(N2,N1),R(N2,N1),UU(N2,N1),BET(N2)

      DO 10 I=LL1,LL2
      BET(I)  = 1./B(I,L1)
      UU(I,L1)= R(I,L1)*BET(I)
   10 CONTINUE

      DO 20 J=L1+1,L2
      DO 20 I=LL1,LL2
      GAM(I,J)= C(I,J-1)*BET(I)
      BET(I)  = 1./(B(I,J)-A(I,J)*GAM(I,J))
      UU(I,J) = (R(I,J)-A(I,J)*UU(I,J-1))*BET(I)
   20 CONTINUE

      DO 30 J=L2-1,L1,-1
      DO 30 I=LL1,LL2
      UU(I,J) = UU(I,J)-GAM(I,J+1)*UU(I,J+1)
   30 CONTINUE

      RETURN
      END

!=======================================================================
      SUBROUTINE TRDIAG2(A,B,C,R,UU,L1,L2,LL1,LL2)
!=======================================================================
!
!     SOLVE THE TRIDIAGONAL MATRIX (X-2 TRIDIAGONAL) WITH
!     DIFFERENT COEFFICIENTS in LHSU,LHSV,LHSW subroutines
!
!     Variables
!       A(),B(),C()  : Coefficients of tridiagonal matrix
!       R()          : RHS
!       UU()         : Solution
!       L1,L2  : 2nd index indicating x2-direction (TDMA)
!       LL2,LL2: 1st index indicating x1-direction
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8 :: I,J
      INTEGER*8 :: L1,L2,LL1,LL2
      REAL*8    :: GAM(N1,N2),A(N1,N2),B(N1,N2),C(N1,N2),R(N1,N2),UU(N1,N2),BET(N1)

      DO 10 I=LL1,LL2
      BET(I)=1./B(I,L1)
      UU(I,L1)=R(I,L1)*BET(I)
   10 CONTINUE

      DO 20 J=L1+1,L2
      DO 20 I=LL1,LL2
      GAM(I,J)=C(I,J-1)*BET(I)
      BET(I)=1./(B(I,J)-A(I,J)*GAM(I,J))
      UU(I,J)=(R(I,J)-A(I,J)*UU(I,J-1))*BET(I)
   20 CONTINUE

      DO 30 J=L2-1,L1,-1
      DO 30 I=LL1,LL2
      UU(I,J)=UU(I,J)-GAM(I,J+1)*UU(I,J+1)
   30 CONTINUE

      RETURN
      END

!=======================================================================
      SUBROUTINE TRDIAG3(A,B,C,R,UU,L1,L2,LL1,LL2)
!=======================================================================
!
!     SOLVE THE TRIDIAGONAL MATRIX (X-3 TRIDIAGONAL) WITH
!     DIFFERENT COEFFICIENTS in LHSU,LHSV,LHSW subroutines
!
!     Variables
!       A(),B(),C()  : Coefficients of tridiagonal matrix
!       R()          : RHS
!       UU()         : Solution
!       L1,L2  : 2nd index indicating x3-direction (TDMA)
!       LL2,LL2: 1st index indicating x1-direction
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8 :: I,J
      INTEGER*8 :: L1,L2,LL1,LL2
      REAL*8    :: GAM(N1,N3),A(N1,N3),B(N1,N3),C(N1,N3),R(N1,N3),UU(N1,N3),BET(N1)

      DO 10 I=LL1,LL2
      BET(I)=1./B(I,L1)
      UU(I,L1)=R(I,L1)*BET(I)
   10 CONTINUE

      DO 20 J=L1+1,L2
      DO 20 I=LL1,LL2
      GAM(I,J)=C(I,J-1)*BET(I)
      BET(I)=1./(B(I,J)-A(I,J)*GAM(I,J))
      UU(I,J)=(R(I,J)-A(I,J)*UU(I,J-1))*BET(I)
   20 CONTINUE

      DO 30 J=L2-1,L1,-1
      DO 30 I=LL1,LL2
      UU(I,J)=UU(I,J)-GAM(I,J+1)*UU(I,J+1)
   30 CONTINUE

      RETURN
      END

!=======================================================================
      SUBROUTINE TRDIAG1P(A,B,C,F,J1,J2,L1,L2)
!=======================================================================
!
!     If XPRDIC = 1 (X periodicity),
!     INVERT A PERIODIC MATRIX (X-1 TRIDIAGONAL) WITH
!     DIFFERENT COEFFICIENTS in LHSU,LHSV,LHSW subroutines
!
!     Variables
!       A(),B(),C()  : Coefficients of tridiagonal matrix
!       R()          : RHS
!       UU()         : Solution
!       L1,L2  : 2nd index indicating x1-direction (periodic TDMA)
!       LL2,LL2: 1st index indicating x2-direction
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8 :: I,J,K,J1,J2,L1,L2,JA,JJ
      REAL*8    :: A(N2,N1),B(N2,N1),C(N2,N1),F(N2,N1)
      REAL*8    :: Q(N2,N1),S(N2,N1),QE(N2,N1),FN(N2),PN(N2)
      REAL*8    :: BINV

      JA=J1+1
      JJ=J1+J2
      DO 10 K=L1,L2
      BINV= 1./B(K,J1)
      Q(K,J1)=-C(K,J1)*BINV
      S(K,J1)=-A(K,J1)*BINV
      FN(K)=F(K,J2)
      F(K,J1)=F(K,J1)*BINV
   10 CONTINUE

!     FORWARD ELIMINATION SWEEP
      DO 20 J=JA,J2
      DO 20 K=L1,L2
      PN(K)=1./(B(K,J)+A(K,J)*Q(K,J-1))
      Q(K,J)=-C(K,J)*PN(K)
      S(K,J)=-A(K,J)*S(K,J-1)*PN(K)
      F(K,J)=(F(K,J)-A(K,J)*F(K,J-1))*PN(K)
   20 CONTINUE

!     BACKWARD PASS
      DO 30 K=L1,L2
      S(K,J2) = 1.
      QE(K,J2)= 0.
   30 CONTINUE
      DO 40 I=JA,J2
      J=JJ-I
      DO 40 K=L1,L2
      S(K,J)=S(K,J)+Q(K,J)*S(K,J+1)
      QE(K,J)=F(K,J)+Q(K,J)*QE(K,J+1)
   40 CONTINUE
      DO 50 K=L1,L2
      F(K,J2)=(FN(K)-C(K,J2)*QE(K,J1)-A(K,J2)*QE(K,J2-1)) &
             /(C(K,J2)*S(K,J1)+A(K,J2)*S(K,J2-1)+B(K,J2))
   50 CONTINUE

!     BACKWARD ELIMINATION PASS
      DO 60 I=JA,J2
      J=JJ-I
      DO 60 K=L1,L2
      F(K,J)=F(K,J2)*S(K,J)+QE(K,J)
   60 CONTINUE

      RETURN
      END

!=======================================================================
      SUBROUTINE TRDIAG3P(A,B,C,F,J1,J2,L1,L2)
!=======================================================================
!
!     If ZPRDIC = 1 (Z periodicity),
!     INVERT A PERIODIC MATRIX (X-3 TRIDIAGONAL) WITH
!     DIFFERENT COEFFICIENTS in LHSU,LHSV,LHSW subroutines
!
!     Variables
!       A(),B(),C()  : Coefficients of tridiagonal matrix
!       R()          : RHS
!       UU()         : Solution
!       L1,L2  : 2nd index indicating x3-direction (periodic TDMA)
!       LL2,LL2: 1st index indicating x2-direction
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8 :: I,J,K
      INTEGER*8 :: J1,J2,L1,L2,JA,JJ
      REAL*8    :: A(N1,N3),B(N1,N3),C(N1,N3),F(N1,N3)
      REAL*8    :: Q(N1,N3),S(N1,N3),QE(N1,N3),FN(N1),PN(N1)
      REAL*8    :: BINV

      JA=J1+1
      JJ=J1+J2
      DO 10 K=L1,L2
      BINV=1./B(K,J1)
      Q(K,J1)=-C(K,J1)*BINV
      S(K,J1)=-A(K,J1)*BINV
      FN(K)=F(K,J2)
      F(K,J1)=F(K,J1)*BINV
   10 CONTINUE

!     FORWARD ELIMINATION SWEEP
      DO 20 J=JA,J2
      DO 20 K=L1,L2
      PN(K)=1./(B(K,J)+A(K,J)*Q(K,J-1))
      Q(K,J)=-C(K,J)*PN(K)
      S(K,J)=-A(K,J)*S(K,J-1)*PN(K)
      F(K,J)=(F(K,J)-A(K,J)*F(K,J-1))*PN(K)
   20 CONTINUE

!     BACKWARD PASS
      DO 30 K=L1,L2
      S(K,J2)=1.
      QE(K,J2)=0.
   30 CONTINUE
      DO 40 I=JA,J2
      J=JJ-I
      DO 40 K=L1,L2
      S(K,J)=S(K,J)+Q(K,J)*S(K,J+1)
      QE(K,J)=F(K,J)+Q(K,J)*QE(K,J+1)
   40 CONTINUE
      DO 50 K=L1,L2
      F(K,J2)=(FN(K)-C(K,J2)*QE(K,J1)-A(K,J2)*QE(K,J2-1))  &
             /(C(K,J2)*S(K,J1)+A(K,J2)*S(K,J2-1)+B(K,J2))
   50 CONTINUE

!     BACKWARD ELIMINATION PASS
      DO 60 I=JA,J2
      J=JJ-I
      DO 60 K=L1,L2
      F(K,J)=F(K,J2)*S(K,J)+QE(K,J)
   60 CONTINUE

      RETURN
      END

!=======================================================================
      SUBROUTINE PTDIAG1A(A,B,C,D,E,F,UU,I1,IN,J1,JN)
!=======================================================================
!
!     SOLVE THE PENTADIAGONAL MATRIX (X-1 PENTADIAGONAL)
!
!     Variables
!       A(),B(),C(),D(),E(),F(): Coefficients of tridiagonal matrix
!       UU()   : RHS and solution
!       I1,IN  : 1st index indicating x2-direction
!       J1,JN  : 2nd index indicating x1-direction
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8 :: I,J,K
      INTEGER*8 :: I1,IN,J1,JN,J2,J3,JNM,JM,JMM
      REAL*8    ::  A(N2,N1),B(N2,N1),C(N2,N1),D(N2,N1),E(N2,N1)
      REAL*8    ::  O(N2,N1),Q(N2,N1),R(N2,N1),F(N2,N1),UU(N2,N1)
      REAL*8    ::  PDENO2(N2),PDENOJ(N2)

      J2=J1+1
      J3=J2+1
      JNM=JN-1
      DO 1 I=I1,IN
      O(I,J1)=-B(I,J1)/A(I,J1)
      Q(I,J1)=-C(I,J1)/A(I,J1)
      R(I,J1)=F(I,J1)/A(I,J1)
    1 CONTINUE

      DO 2 I=I1,IN
      PDENO2(I)=1./(A(I,J2)+D(I,J2)*O(I,J1))
      O(I,J2)=-(B(I,J2)+D(I,J2)*Q(I,J1))*PDENO2(I)
      Q(I,J2)=-C(I,J2)*PDENO2(I)
      R(I,J2)=(F(I,J2)-D(I,J2)*R(I,J1))*PDENO2(I)
    2 CONTINUE

      DO 10 J=J3,JN
      JM=J-1
      JMM=JM-1
      DO 10 I=I1,IN
      PDENOJ(I)=1./(A(I,J)+E(I,J)*Q(I,JMM)                            &
                  +(D(I,J)+E(I,J)*O(I,JMM))*O(I,JM))
      O(I,J)=-(B(I,J)+(D(I,J)+E(I,J)*O(I,JMM))*Q(I,JM))*PDENOJ(I)
      Q(I,J)=-C(I,J)*PDENOJ(I)
      R(I,J)=(F(I,J)-E(I,J)*R(I,JMM)-(D(I,J)+E(I,J)*O(I,JMM))*R(I,JM)) &
            *PDENOJ(I)
   10 CONTINUE

      DO 11 I=I1,IN
      UU(I,JN)=R(I,JN)
   11 CONTINUE

      DO 12 I=I1,IN
      UU(I,JNM)=O(I,JNM)*UU(I,JN)+R(I,JNM)
   12 CONTINUE

      DO 20 J=JNM-1,J1,-1
      DO 20 I=I1,IN
      UU(I,J)=O(I,J)*UU(I,J+1)+Q(I,J)*UU(I,J+2)+R(I,J)
   20 CONTINUE

      RETURN
      END
!=======================================================================
      SUBROUTINE RETRV_UVW
!=======================================================================
!
!     Adjust boundary velocity to satisfy global mass conservation
!           [retrieve uvw (at the exit boundary)]
!     See 'CONVBC' subroutine in lica_[bodyname, e.g. cylinder].f90 file
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W
      IMPLICIT NONE
      INTEGER*8    :: I,J,K

!$OMP PARALLEL DO
      DO K=1,N3M
      DO J=1,N2M
      U(N1,J,K)=UOUT(J,K)
      ENDDO
      ENDDO
!$OMP PARALLEL DO
      DO K=1,N3M
      DO J=2,N2M
      V(N1,J,K)=VOUT(J,K)
      ENDDO
      ENDDO
!$OMP PARALLEL DO
      DO K=K_BGPZ,N3M
      DO J=1,N2M
      W(N1,J,K)=WOUT(J,K)
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE RETRV_UVW
!=======================================================================
!=======================================================================
      SUBROUTINE PRDIC_ADJ_UVW(ISP)
!=======================================================================
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,P
      IMPLICIT NONE
      INTEGER*8    :: ISP
      INTEGER*8    :: I,J,K

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
      IF (ISP .NE. 0) THEN
!$OMP PARALLEL DO
        DO K=1,N3M
        DO J=1,N2M
           P(0 ,J,K)=P(N1M,J,K)
           P(N1,J,K)=P(1  ,J,K)
        ENDDO
        ENDDO
      ENDIF
      ENDIF

!     Y PERIODICITY
      IF (YPRDIC .EQ. 1) THEN
!$OMP PARALLEL DO
      DO K=1,N3
      DO I=0,N1
         U(I ,0,K)=U(I,N2M,K)
         U(I,N2,K)=U(I,  1,K)
      ENDDO
      ENDDO
!$OMP PARALLEL DO
      DO K=0,N3
      DO I=0,N1
         V(I ,0,K)=V(I,N2M,K)
         V(I,N2,K)=V(I,  1,K)
      ENDDO
      ENDDO
!$OMP PARALLEL DO
      DO K=0,N3
      DO I=1,N1
         W(I ,0,K)=W(I,N2M,K)
         W(I,N2,K)=W(I,  1,K)
      ENDDO
      ENDDO
      IF (ISP .NE. 0) THEN
!$OMP PARALLEL DO
        DO K=1,N3M
        DO I=1,N1M
           P(I ,0,K)=P(I,N2M,K)
           P(I,N2,K)=P(I,  1,K)
        ENDDO
        ENDDO
      ENDIF
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
         V(I,J,0 ) =V(I,J,N3M)
         V(I,J,N3) =V(I,J,1  )
      ENDDO
      ENDDO
!$OMP PARALLEL DO
      DO J=0,N2
      DO I=0,N1
         W(I,J,0) =W(I,J,N3M)
         W(I,J,N3)=W(I,J,1)
      ENDDO
      ENDDO
      IF (ISP .NE. 0) THEN
!$OMP PARALLEL DO
        DO J=1,N2M
        DO I=1,N1M
           P(I,J,0 )=P(I,J,N3M)
           P(I,J,N3)=P(I,J,1  )
        ENDDO
        ENDDO
      ENDIF
      ENDIF

      RETURN
      END SUBROUTINE PRDIC_ADJ_UVW
!=======================================================================