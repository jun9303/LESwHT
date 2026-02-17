!=======================================================================
!
!                     POISSON EQUATION SOLVER
!
!     THE POISSON EQUATION IS SOLVED USING A COMBINATION OF
!     A FOURIER TRANSFORM METHOD AND A MULTI-GRID METHOD
!
!     A SET OF TWO DIMENSIONAL HELMHOLTZ EQUATIONS
!     RESULTED FROM FOURIER TRANSFORM    IN Z DIRECTION
!     ARE SOLVED BY MULTI-GRID ITERATION IN X DIRECTION
!     INCORPORATED WITH TDMA INVERSION   IN Y DIRECTION
!
!     ALTHOUGH A UNIFORM GRID IN Z DIRECTION IS REQUIRED
!     TO UTILIZE FOURIER TRANSFORM, THIS SOLVER SHOWS ABOUT
!     30% FASTER CONVERGENCE SPEED THAN TWO DIMENSIONAL (X,Z)
!     MULTI-GRID SOLVER
!
!     IN ADDITION, GSOR METHOD INSTEAD OF MULTI-GRID METHOD
!     CAN BE SELECTED
!
!     TURBULENCE AND FLOW CONTROL LAB.
!     SEOUL NATIONAL UNIVERSITY
!     AUGUST 27, 2007
!     JUNGIL LEE
!     Ph. D. STUDENT
!
!     June 2017, J. Park: allocatable variables and f90
!
!=======================================================================
      SUBROUTINE POISINIT
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER*8   :: I,J,K,ILEV

       CALL Z_FT_ALLO

       ALLOCATE (AC (N1MD,N2M,N3MH))
       ALLOCATE (GAM(N1MD,N2M,N3MH))
       ALLOCATE (BET(N1MD,N2M,N3MH))
 
       CALL PMAT
!------MULTIGRID METHOD
       CALL COEFMG

       DO ILEV=0,NLEV
        DO K=1,N3MH
         DO I=IIMG(ILEV,1),IIMG(ILEV,2)
          BET(I,1,K)=1./AC(I,1,K)
         ENDDO
         DO J=2,N2M
         DO I=IIMG(ILEV,1),IIMG(ILEV,2)
          GAM(I,J,K)= AN(J-1)*BET(I,J-1,K)
          BET(I,J,K)= 1./(AC(I,J,K)-AS(J)*GAM(I,J,K))
         ENDDO
         ENDDO
        ENDDO
       ENDDO

      OPEN(77,FILE='../output/ftr/mgftresiduemax.dat')

      RETURN
      END
!=======================================================================
      SUBROUTINE PMAT
!=======================================================================
! --- CONSTRUCT MATRIX FOR POISSON EQ. ---------------------------------
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER*8   :: I,J,K,IP,JP

!-----FOR TOP LEVEL
      DO 10 I=1,N1M
       IP=IPV(I)
       AW(I)=(1.-FIXIL(I))*C2CXI(I)*F2FXI(I)
       AE(I)=(1.-FIXIU(I))*C2CXI(IP)*F2FXI(I)
   10 CONTINUE

      DO 20 J=1,N2M
       JP=JPV(J)
       AS(J)=(1.-FIXJL(J))*C2CYI(J)*F2FYI(J)
       AN(J)=(1.-FIXJU(J))*C2CYI(JP)*F2FYI(J)
   20 CONTINUE

      DO 30 J=1,N2M
      DO 30 I=1,N1M
       AC(I,J,1)=-1.*(AW(I)+AE(I)+AS(J)+AN(J))
   30 CONTINUE

      N3MH=N3M/2+1
      IF(N3M .GT. 1) CALL MWAVENUMBER     ! INIT. MODIFIED WAVE #.

      DO 40 K=2,N3MH
      DO 40 J=1,N2M
      DO 40 I=1,N1M
      AC(I,J,K)=AC(I,J,1)-AK3(K)
   40 CONTINUE

      RETURN
      END
!=======================================================================
      SUBROUTINE MWAVENUMBER
!=======================================================================
! --- MODIFIED WAVE NUMBER DEFINITION
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER*8   :: K
       REAL*8      :: PI

      PI = ACOS(-1.)

      DO K=1,N3MH
       AK3(K)= 2.*(1.-COS(2.*FLOAT(K-1)*PI/FLOAT(N3M)))*F2FZI(1)*F2FZI(1)
      ENDDO

      RETURN
      END
!=======================================================================
      SUBROUTINE POISSON(PHI,DIVGSUM)
!=======================================================================
! MAIN SOLVER OF POISSON EQUATION
! AX=b
! A: coefficient of discretized poisson equation
! X: PHI: Pseudo pressure, Output of subroutine POISSON
! b: DIVGSUM: OUTPUT of subroutine DIVGS.
! x-direction: Gauss-Seidal + Multigrid
! z-direction: Fourier transform
! y-direction: TDMA
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER*8   :: I,J,K,KKK
       REAL*8      :: PHI(0:N1,0:N2,0:N3),DIVGSUM(0:N1,0:N2,0:N3)
       COMPLEX*16  :: CCAP(N1,N2,N3MH)
       COMPLEX*16, dimension (:,:), allocatable :: XXX,CP,TMP
       COMPLEX*16, dimension (:),   allocatable :: XXXX,XXXX_B

       REAL*8      :: TEST

      IF(N3M.EQ.1) THEN
       allocate (CP(0:N1,0:N2),TMP(N1,N2))
       CP= 0.
!$OMP PARALLEL DO
       DO J=1,N2M
       DO I=1,N1M
        TMP(I,J)=DIVGSUM(I,J,1)
       ENDDO
       ENDDO

       CALL MG2D(CP,TMP,1,TEST1,0.)

!$OMP PARALLEL DO
       DO J=1,N2M
       DO I=1,N1M
        PHI(I,J,1)=CP(I,J)
       ENDDO
       ENDDO
       deallocate(CP,TMP)
       GOTO 400
      ENDIF


! --- DO THE FORWARD FFT
!$OMP PARALLEL &
!$OMP private(XXX,XXXX,XXXX_B)
      allocate(XXX(N3M,N1M))
      allocate(XXXX(N3M),XXXX_B(N3M*2))
      CALL ZFFT1D(XXXX,N3M, 0,XXXX_B)
!$OMP DO
      DO 100 J=1,N2M
        DO K=1,N3M
        DO I=1,N1M
         XXX(K,I)=DIVGSUM(I,J,K)
        ENDDO
        ENDDO
        
        DO I=1,N1M
         CALL ZFFT1D(XXX(1,I),N3M,-1,XXXX_B)
        END DO
        
        DO K=1,N3MH
        DO I=1,N1M
         CCAP(I,J,K)=XXX(K,I)
        ENDDO
        ENDDO
  100 CONTINUE
!$OMP END DO
      deallocate(XXX,XXXX,XXXX_B)
!$OMP END PARALLEL


! --- SOLVE A SET OF POISSON EQS.
       TEST=TEST1/FLOAT(N3MH)*FLOAT(N3M)*0.8   ! convergence criteria
!$OMP PARALLEL &
!$OMP private(CP)
      allocate(CP(0:N1,0:N2))
!$OMP DO
      DO 200 KKK=1,N3MH
        CP= 0.
        IF (KKK.LE.IMGSOR) THEN    
         CALL MG2D(CP,CCAP(1,1,KKK),KKK,TEST,0.)
        ELSE
         CALL GSOR2D(CP,CCAP(1,1,KKK),KKK,TEST,0.)
        ENDIF
! Multigrid acceleration is efficient and effective at low wavenumber.
! For the high wavenumber, GS+SOR can be used. 
! The number of operation of (GS+SOR) is smaller than that of (GS+MG). 
        
        DO J=1,N2M
        DO I=1,N1M
          CCAP(I,J,KKK)=CP(I,J)
        ENDDO
        ENDDO
  200 CONTINUE
!$OMP END DO
      deallocate(CP)
!$OMP END PARALLEL


! --- DO THE INVERSE FFT
!$OMP PARALLEL &
!$OMP private(XXX,XXXX,XXXX_B)
      allocate(XXX(N3M,N1M))
      allocate(XXXX(N3M),XXXX_B(N3M*2))
      CALL ZFFT1D(XXXX,N3M, 0,XXXX_B)
!$OMP DO
      DO 300 J=1,N2M
        DO K=1,N3MH
        DO I=1,N1M
          XXX(K,I)= CCAP(I,J,K)
        ENDDO
        ENDDO

        DO K=N3MH+1,N3M
        DO I=1,N1M
         XXX(K,I) = CONJG(CCAP(I,J,N3M+2-K))
        ENDDO
        ENDDO
        
        DO I=1,N1M
          CALL ZFFT1D(XXX(1,I),N3M, 1,XXXX_B)
        ENDDO
        
        DO K=1,N3M
        DO I=1,N1M
          PHI(I,J,K)= REAL(XXX(K,I))
        ENDDO
        ENDDO
  300 CONTINUE
!$OMP END DO
      deallocate(XXX,XXXX,XXXX_B)
!$OMP END PARALLEL

  400 CONTINUE

      RETURN
      END

!=======================================================================
       SUBROUTINE MG2D(PC,RHS,KV,TEST,OLDV)
!=======================================================================
!     MULTIGRID ENVIRONMENT VARIABLES
!     PC  : SOLUTION OF THE MG2D.
!     RHS : RHS OF THE POISSON EQUATION FOR EACH WAVENUMBER
!     KV  : wavenumber index
!     TEST    : CONDITION FOR CONVERGENCE
!     OLDV    : 0; INITALIZE TO ZERO, 1; USE PREVIOUS SOLUTION
!     NLEV    : TOTAL CELS = (MINROW1+MINROW2)*(2**NLEV)
!     LEVHALF : HALF OF NLEV
!     IWC     : FLAG FOR USING W-CYCLE ( 1 IF YOU WANT TO USE W-CYCLE)
!     NBLI    : NUMBER OF BASE LEVEL ITERATION
!     MGITR   : NUMBER OF MAXIMUM ITERATIONS
!     IIMG(ILEV,1)  : START POINT OF INDEX AT ILEV LEVEL
!     IIMG(ILEV,2)  : END POINT OF INDEX AT ILEV LEVEL
!     IIMG(ILEV,3)  : NUMBER OF ROW AT ILEV LEVEL

!     Subroutines
!     1. TOPLEVEL: iteration at the highest level (original grid system)
!     2. RELAX   : smoothing + obtating residue ate each level
!     3. GODOWN  : restriction
!     4. GOUP    : elongation
!-----------------------------------------------------------------------
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER*8   :: I,J,KV,KK,ILEV
       REAL*8      :: TEST,SUMRES,OLDV
       COMPLEX*16  :: PC(0:N1,0:N2),RHS(N1,N2)
       COMPLEX*16  :: RESD(N1MD,N2M),GGII(0:N1MD,0:N2)

       KK = 0
       RESD = 0.
       GGII = 0.

       CALL TOPLEVEL(0,PC,RHS,SUMRES,OLDV,TEST,KV,RESD,GGII)
       IF(SUMRES .LT. TEST) GOTO 205


       DO 10 KK= 1,MGITR            ! main iteration
         DO ILEV=NLEV-1,1,-1 
           CALL RELAX(ILEV,0.,1,KV,RESD,GGII) 
           CALL GODOWN(ILEV,KV,RESD,GGII)
         ENDDO

         CALL RELAX(0,0.,NBLI,KV,RESD,GGII)

         DO ILEV= 0,NLEV-2
          CALL GOUP(ILEV,RESD,GGII)
          CALL RELAX(ILEV+1,1.,1,KV,RESD,GGII)
         ENDDO

         CALL TOPLEVEL(1,PC,RHS,SUMRES,1.,TEST,KV,RESD,GGII)

         IF(SUMRES .LT. TEST) GOTO 205

   10   CONTINUE
        PRINT *,'ITERATION LIMIT EXCEEDED.'
  205   CONTINUE

        IF ((KV.LE.3)) THEN
        WRITE(77,999) KV,KK,SUMRES*DTCONST!*float(n3mh)
        ENDIF

  999   FORMAT('KV=',I4,3X,'KK=',I4,3X,'RM=',ES24.16)

      RETURN
      END
!=======================================================================
      SUBROUTINE TOPLEVEL(ID,PC,RHS,SUMRES,OLDV,TEST,KV,RESD,GGII)    ! Zebra Version
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
        INTEGER*8   :: I,J,IC
        INTEGER*8   :: ID,KV
        REAL*8      :: SUMRES,TEST,OLDV
        COMPLEX*16  :: TT
        COMPLEX*16  :: PC(0:N1,0:N2),RHS(N1,N2)
        COMPLEX*16  :: RESD(N1MD,N2M),GGII(0:N1MD,0:N2)

      IF(ID .NE. 1) GOTO 25

!     INTERPOLATE & ADD
      IC= 0
      DO I=IIMG(NLEV-1,1),IIMG(NLEV-1,2)-1
      IC=IC+2
      DO J=1,N2M
      PC(IC,J)=PC(IC,J)+COI1(I)*GGII(I,J)+COI2(I)*GGII(I+1,J)
      ENDDO
      ENDDO

      I=IIMG(NLEV-1,2)
      IC=IC+2
      DO J=1,N2M
      PC(IC,J)=PC(IC,J)+COI1(I)*GGII(I,J)      ! use one point COI(~)=1.
      ENDDO

  25  CONTINUE

!  RELAX
      DO J=1,N2M
      DO I=1,N1M,2
       GGII(I,J)=RHS(I,J)-OLDV*(AW(I)*PC(I-1,J)+AE(I)*PC(I+1,J))
      ENDDO
      ENDDO

      CALL TRDIAG1M(GGII,GGII,1,N1M,1,N2M,KV)

      IF(KV .EQ. 1) THEN
        TT= GGII(1,N2M-1)
      ELSE
        TT= 0.
      END IF

      DO J=1,N2M
      DO I=2,N1M,2
        GGII(I-1,J)=GGII(I-1,J)-TT
        GGII(I,J)=RHS(I,J)-AW(I)*GGII(I-1,J)-AE(I)*(GGII(I+1,J)-TT)
      ENDDO
      ENDDO

      CALL TRDIAG1M(GGII,GGII,2,N1M,1,N2M,KV)

!  CALCULATE RESIDUAL
      SUMRES= 0.

      DO J=1,N2M
      DO I=1,N1M
      RESD(I,J)=RHS(I,J)-AW(I)*GGII(I-1,J)-AE(I)*GGII(I+1,J)  &
                   -AS(J)*GGII(I,J-1)-AN(J)*GGII(I,J+1)       &
                   -AC(I,J,KV)*GGII(I,J)
      SUMRES= AMAX1( SUMRES,ABS(RESD(I,J)) )
      PC(I,J)=GGII(I,J)
      ENDDO
      ENDDO

      IF(SUMRES .LT. TEST) GOTO 99
      IF(ID .EQ. 2) GOTO 99

!  RESTRICT
      IC=-1

      DO I=IIMG(NLEV-1,1),IIMG(NLEV-1,2)
      IC=IC+2
      DO J=1,N2M
      RESD(I,J)=RESD(IC,J)*COR1(I)+RESD(IC+1,J)*COR2(I)
      ENDDO
      ENDDO

  99  RETURN
      END
!=======================================================================
      SUBROUTINE TRDIAG1M(RR,UU,L1,L2,LL1,LL2,KV)
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
        INTEGER*8   :: I,J,L1,L2,LL1,LL2,KV
        COMPLEX*16  :: RR(0:N1MD,0:N2),UU(0:N1MD,0:N2)

      DO I=L1,L2,2
        UU(I,LL1)=RR(I,LL1)*BET(I,1,KV)
      ENDDO

      DO J=LL1+1,LL2
      DO I=L1,L2,2
        UU(I,J)=(RR(I,J)-AS(J)*UU(I,J-1))*BET(I,J,KV)
      ENDDO
      ENDDO
      DO J=LL2-1,LL1,-1
      DO I=L1,L2,2
        UU(I,J)=UU(I,J)-GAM(I,J+1,KV)*UU(I,J+1)
      ENDDO
      ENDDO

      RETURN
      END

!=======================================================================
      SUBROUTINE GOUP(ILEV,RESD,GGII)! INTERPOLATE RESIDUAL & ADD IT TO HIGH LEVEL
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
        INTEGER*8   :: I,J,ILEV,IBGH
        COMPLEX*16  :: RESD(N1MD,N2M),GGII(0:N1MD,0:N2)

      IBGH=IIMG(ILEV+1,1)-1

      DO I=IIMG(ILEV,1),IIMG(ILEV,2)-1
      IBGH=IBGH+2
      DO J=1,N2M
      GGII(IBGH,J)=GGII(IBGH,J)+COI1(I)*GGII(I,J)+COI2(I)*GGII(I+1,J)
      ENDDO
      ENDDO

      I=IIMG(ILEV,2)
      IBGH=IBGH+2
      DO J=1,N2M
      GGII(IBGH,J)=GGII(IBGH,J)+COI1(I)*GGII(I,J)           ! use one point
      ENDDO


      RETURN
      END

!=======================================================================
      SUBROUTINE GODOWN(ILEV,KV,RESD,GGII)        ! COMPUTE RESIDUAL & RESTRICT IT
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
        INTEGER*8   :: I,J,ILEV,KV,IBG
        COMPLEX*16  :: RESD(N1MD,N2M),GGII(0:N1MD,0:N2)

      DO J=1,N2M
      DO I=IIMG(ILEV,1),IIMG(ILEV,2)
      RESD(I,J)=RESD(I,J)-AW(I)*GGII(I-1,J)-AE(I)*GGII(I+1,J)   &
               -AS(J)*GGII(I,J-1)-AN(J)*GGII(I,J+1)-AC(I,J,KV)*GGII(I,J)
      ENDDO
      ENDDO

      IBG=IIMG(ILEV,1)-2

      DO I=IIMG(ILEV-1,1),IIMG(ILEV-1,2)
      IBG=IBG+2
      DO J=1,N2M
      RESD(I,J)=RESD(IBG,J)*COR1(I)+RESD(IBG+1,J)*COR2(I)
      ENDDO
      ENDDO

      RETURN
      END

!=======================================================================
      SUBROUTINE RELAX(ILEV,OLDV,IITER,KV,RESD,GGII)   ! Zebra Version
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
        INTEGER*8   :: I,J,ILEV,IITER,KV,II
        COMPLEX*16  :: RESD(N1MD,N2M),GGII(0:N1MD,0:N2)
        REAL*8      :: OLDV

      DO J=1,N2M
      DO I=IIMG(ILEV,1),IIMG(ILEV,2),2
      GGII(I,J)=RESD(I,J)-OLDV*(AW(I)*GGII(I-1,J)+AE(I)*GGII(I+1,J))
      ENDDO
      ENDDO

      CALL TRDIAG1M(GGII,GGII,IIMG(ILEV,1),IIMG(ILEV,2),1,N2M,KV)

      DO J=1,N2M
      DO I=IIMG(ILEV,1)+1,IIMG(ILEV,2),2
      GGII(I,J)=RESD(I,J)-AW(I)*GGII(I-1,J)-AE(I)*GGII(I+1,J)
      ENDDO
      ENDDO

      CALL TRDIAG1M(GGII,GGII,IIMG(ILEV,1)+1,IIMG(ILEV,2),1,N2M,KV)

      DO 50 II=1,IITER-1

      DO J=1,N2M
      DO I=IIMG(ILEV,1),IIMG(ILEV,2),2
      GGII(I,J)=RESD(I,J)-(AW(I)*GGII(I-1,J)+AE(I)*GGII(I+1,J))
      ENDDO
      ENDDO

      CALL TRDIAG1M(GGII,GGII,IIMG(ILEV,1),IIMG(ILEV,2),1,N2M,KV)

      DO J=1,N2M
      DO I=IIMG(ILEV,1)+1,IIMG(ILEV,2),2
      GGII(I,J)=RESD(I,J)-AW(I)*GGII(I-1,J)-AE(I)*GGII(I+1,J)
      ENDDO
      ENDDO

      CALL TRDIAG1M(GGII,GGII,IIMG(ILEV,1)+1,IIMG(ILEV,2),1,N2M,KV)

  50  CONTINUE
      RETURN
      END
!=======================================================================
      SUBROUTINE GSOR2D(U,RHS,KV,TEST,OLDV)    ! 1 EQ. TYPE
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER*8    :: I,J,KV
       REAL*8       :: TEST,OLDV
       COMPLEX*16   :: U(0:N1,0:N2),RHS(N1,N2)
       COMPLEX*16   :: GGII(0:N1MD,0:N2)
       INTEGER*8    :: KK
       REAL*8       :: WW,WW2,ERRMAX


      WW = WWSOR
      WW2= 1.-WW
      KK = 0

!  HALF RELAX
!  ----------

      DO J=1,N2M
      DO I=1,N1M,2
       GGII(I,J)=RHS(I,J)-OLDV*(AW(I)*U(I-1,J)+AE(I)*U(I+1,J))
      ENDDO
      ENDDO

      CALL TRDIAG1M(GGII,GGII,1,N1M,1,N2M,KV)

      DO J=1,N2M
      DO I=1,N1M,2
       U(I,J)=WW*GGII(I,J)+OLDV*WW2*U(I,J)
      ENDDO
      ENDDO

!  ANOTHER HALF
!  ------------


      DO J=1,N2M
      DO I=2,N1M,2
       GGII(I,J)=RHS(I,J)-(AW(I)*U(I-1,J)+AE(I)*U(I+1,J))
      ENDDO
      ENDDO

      CALL TRDIAG1M(GGII,GGII,2,N1M,1,N2M,KV)

      DO J=1,N2M
      DO I=2,N1M,2
       U(I,J)=WW*GGII(I,J)+OLDV*WW2*U(I,J)
      ENDDO
      ENDDO

      CALL RESID3(U,RHS,KV,ERRMAX)
      IF(ERRMAX .LT. TEST) GOTO 88

!  MAIN ITERATION
!  ==============
      DO 100 KK=1,MGITR

!  HALF RELAX
!  ----------

      DO J=1,N2M
      DO I=1,N1M,2
      GGII(I,J)=RHS(I,J)-(AW(I)*U(I-1,J)+AE(I)*U(I+1,J))
      ENDDO
      ENDDO

      CALL TRDIAG1M(GGII,GGII,1,N1M,1,N2M,KV)

      DO J=1,N2M
      DO I=1,N1M,2
       U(I,J)=WW*GGII(I,J)+WW2*U(I,J)
      ENDDO
      ENDDO

!  ANOTHER HALF
!  ------------

      DO J=1,N2M
      DO I=2,N1M,2
       GGII(I,J)=RHS(I,J)-(AW(I)*U(I-1,J)+AE(I)*U(I+1,J))
      ENDDO
      ENDDO

      CALL TRDIAG1M(GGII,GGII,2,N1M,1,N2M,KV)

      DO J=1,N2M
      DO I=2,N1M,2
       U(I,J)=WW*GGII(I,J)+WW2*U(I,J)
      ENDDO
      ENDDO

      CALL RESID3(U,RHS,KV,ERRMAX)
      WRITE(77,999) KV,kk,ERRMAX*dtconst*float(n3mh)
      IF(ERRMAX .LT. TEST) GOTO 88

100   CONTINUE

      PRINT *,'ITERATION LIMIT EXCEEDED.'

88    CONTINUE

999   FORMAT('KV=',I3,3X,'KK=',I3,3X,'RG=',ES24.16)

99    RETURN
      END

!=======================================================================
      SUBROUTINE RESID3(U,RHS,KV,ERRMAX)
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER*8    :: I,J,KV
       REAL*8       :: ERRMAX
       COMPLEX*16   :: U(0:N1,0:N2),RHS(N1,N2)
       COMPLEX*16   :: ERR

       ERRMAX= 0.

      DO J=1,N2M
      DO I=1,N1M
       ERR=RHS(I,J)-AW(I)*U(I-1,J)-AE(I)*U(I+1,J)   &
                   -AS(J)*U(I,J-1)-AN(J)*U(I,J+1)   &
                   -AC(I,J,KV)*U(I,J)
       ERRMAX= AMAX1( ERRMAX,ABS(ERR) )
      ENDDO
      ENDDO

      RETURN
      END
!=======================================================================
      SUBROUTINE COEFMG
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER*8   :: I,J,K,ILEV
       INTEGER*8   :: MINROW,IBG,IEND,ISP,IC,IBGH
       REAL*8      :: IWZ(N1MD),IEZ(N1MD)
       REAL*8      :: XMPM(0:N1MD)


       MINROW=N1M/(2**NLEV)
       XMPM = 0.

       LEVHALF= NINT(NLEV/2.)
       IIMG(NLEV,1)=1
       IIMG(NLEV,2)=N1M
       IIMG(NLEV,3)=N1M

       DO 50 I=1,N1M
        XMPM(I)=XMP(I)
        IWZ(I)=1-1/I                  ! 0 only if I=1
        IEZ(I)=1-I/N1M                ! 0 only if I=N1M
  50   CONTINUE



!  COMPUTE FOR LOWER LEVELS
      DO 100 ILEV=NLEV-1,0,-1

       IIMG(ILEV,1)=IIMG(ILEV+1,1)+IIMG(ILEV+1,3)
       IIMG(ILEV,3)=MINROW*(2**ILEV)
       IIMG(ILEV,2)=IIMG(ILEV,1)+IIMG(ILEV,3)-1

       IBG=IIMG(ILEV,1)
       IEND=IIMG(ILEV,2)

       ISP=2**(NLEV-ILEV)           ! width of one cell at low level

       IC= 0
       DO 110 I=IBG,IEND
        IC=IC+1
        XMPM(I)=0.5*(X(IC*ISP+1)+X((IC-1)*ISP+1))
        IWZ(I)=1-IBG/I                ! 0 only if I=IBG
        IEZ(I)=1-I/IEND               ! 0 only if I=IEND
 110   CONTINUE

       IC= 0
       DO 120 I=IBG,IEND
        IC=IC+1
        AW(I)=IWZ(I)/((XMPM(I)-XMPM(I-1))*(X(IC*ISP+1)-X((IC-1)*ISP+1)))
        AE(I)=IEZ(I)/((XMPM(I+1)-XMPM(I))*(X(IC*ISP+1)-X((IC-1)*ISP+1)))
 120   CONTINUE

       DO 130 J=1,N2M
       DO 130 I=IIMG(ILEV,1),IIMG(ILEV,2)
        AC(I,J,1)=-1.*(AW(I)+AE(I)+AS(J)+AN(J))
 130   CONTINUE

       DO 140 K=2,N3MH
       DO 140 J=1,N2M
       DO 140 I=IIMG(ILEV,1),IIMG(ILEV,2)
        AC(I,J,K)=AC(I,J,1)-AK3(K)
 140   CONTINUE

 100  CONTINUE

!  CALCULATE RESTRICTION COEFFS
       DO 145 ILEV=NLEV,1,-1
        IBGH=IIMG(ILEV,1)
        DO 145 I=IIMG(ILEV-1,1),IIMG(ILEV-1,2)
        COR1(I)=(XMPM(IBGH+1)-XMPM(I))/(XMPM(IBGH+1)-XMPM(IBGH))
        COR2(I)=1.-COR1(I)
        IBGH=IBGH+2
 145   CONTINUE

!  CALCULATE INTERPOLATION COEFFS
       DO 150 ILEV=0,NLEV-1
        IBGH=IIMG(ILEV+1,1)+1
        DO 160 I=IIMG(ILEV,1),IIMG(ILEV,2)-1
        COI1(I)=(XMPM(I+1)-XMPM(IBGH))/(XMPM(I+1)-XMPM(I))  ! * lower value
        COI2(I)=1.-COI1(I)
        IBGH=IBGH+2
 160   CONTINUE
       I=IIMG(ILEV,2)
       COI1(I)= 1.                 ! use only one lower point at upper wall

 150  CONTINUE

!  INITIALIZE WORKING ARRAY
!  ------------------------
!      GI=0.

       OPEN(20,FILE='../output/ftr/poiscoef.out')
       WRITE(20,20)0,0.,0.,0.,XMPM(I),0.,0.
       DO I=1,N1M
       WRITE(20,20)I,IWZ(I),IEZ(I),XMP(I),XMPM(I),AW(I),AE(I)
       ENDDO
       DO I=N1M+1,N1MD
       WRITE(20,20)I,IWZ(I),IEZ(I),0.,XMPM(I),AW(I),AE(I)
       ENDDO
       CLOSE(20)
   20 FORMAT(I10,6ES13.5)

      DO ILEV=NLEV,0,-1
      WRITE(*,*) 'IIMG(1',ILEV,')=',IIMG(ILEV,1)
      WRITE(*,*) 'IIMG(2',ILEV,')=',IIMG(ILEV,2)
      WRITE(*,*) 'IIMG(3',ILEV,')=',IIMG(ILEV,3)
      ENDDO

      RETURN
      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     FFTE: A FAST FOURIER TRANSFORM PACKAGE
!
!     (C) COPYRIGHT SOFTWARE, 2000-2004, ALL RIGHTS RESERVED
!                BY
!         DAISUKE TAKAHASHI
!         GRADUATE SCHOOL OF SYSTEMS AND INFORMATION ENGINEERING
!         UNIVERSITY OF TSUKUBA
!         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
!         E-MAIL: daisuke@cs.tsukuba.ac.jp
!
!
!     1-D COMPLEX FFT ROUTINE (FOR VECTOR MACHINES)
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
!
      SUBROUTINE ZFFT1D(A,N,IOPT,B)
      IMPLICIT REAL*8 (A-H,O-Z)

! The maximum supported 2-D transform length is 65536.
      PARAMETER (NDA2=65536)
! The maximum supported 3-D transform length is 4096.
      PARAMETER (NDA3=4096)
      PARAMETER (NDA4=256)
! The parameter NBLK is a blocking parameter.
      PARAMETER (NBLK=16)
!      PARAMETER (NBLK=8)  (for PentiumIII and Athlon)
!      PARAMETER (NBLK=16) (for Pentium4, Athlon XP, Opteron, Itanium
!                           and ItaniuN2)
! The parameter NP is a padding parameter to avoid cache conflicts in
! the FFT routines.
      PARAMETER (NP=8)
!      PARAMETER (NP=2) (for PentiumIII)
!      PARAMETER (NP=4) (for Athlon, Athlon XP, Opteron and Itanium)
!      PARAMETER (NP=8) (for Pentium4 and ItaniuN2)
! Size of L2 cache
      PARAMETER (L2SIZE=1048576)

      COMPLEX*16 A(*),B(*)
      COMPLEX*16 W1(NDA2/2+NP),W2(NDA2/2+NP)
      DIMENSION IP(3),IP1(3),IP2(3)
      SAVE W1,W2
!
      CALL FACTOR(N,IP)
!
      IF (IOPT .EQ. 1) THEN
        DO 10 I=1,N
          A(I)=CONJG(A(I))
   10   CONTINUE
      END IF
!
      DO 20 I=1,3
        IP1(I)=(IP(I)+1)/2
        IP2(I)=IP(I)-IP1(I)
   20 CONTINUE
      N1=(2**IP1(1))*(3**IP1(2))*(5**IP1(3))
      N2=(2**IP2(1))*(3**IP2(2))*(5**IP2(3))
!
      IF (IOPT .EQ. 0) THEN
        CALL SETTBL(W1,N1)
        CALL SETTBL(W2,N2)
        CALL SETTBL2(B(N+1),N2,N1)
        RETURN
      END IF
!
      CALL MFFT235A(A,B,W2,N1,N2,IP2)
      CALL ZTRANSMUL(A,B,B(N+1),N1,N2)
      CALL MFFT235B(B,A,W1,N2,N1,IP1)
!
      IF (IOPT .EQ. 1) THEN
        DN=1.0D0/DBLE(N)
        DO 30 I=1,N
          A(I)=CONJG(A(I))*DN
   30   CONTINUE
      END IF
      RETURN
      END

!
!     FFTE: A FAST FOURIER TRANSFORM PACKAGE
!
!     (C) COPYRIGHT SOFTWARE, 2000-2004, ALL RIGHTS RESERVED
!                BY
!         DAISUKE TAKAHASHI
!         GRADUATE SCHOOL OF SYSTEMS AND INFORMATION ENGINEERING
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
        J=J+L
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
        J=J+L
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
        J=J+L
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
        J=J+L
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
        J=J+L
   10 CONTINUE
      DO 20 K=1,IP(3)
        L=L/5
        CALL SETTBL0(W(J),5,L)
        J=J+L
   20 CONTINUE
      DO 30 K=1,KP4
        L=L/4
        CALL SETTBL0(W(J),4,L)
        J=J+L
   30 CONTINUE
      DO 40 K=1,IP(2)
        L=L/3
        CALL SETTBL0(W(J),3,L)
        J=J+L
   40 CONTINUE
      RETURN
      END
      SUBROUTINE SETTBL0(W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION W(2,*)
!
      PI2=8.0D0*DATAN(1.0D0)
      PX=-PI2/(DBLE(M)*DBLE(L))
!NODIR$ VECTOR ALIGNED
      DO 10 I=1,L
        W(1,I)=DCOS(PX*DBLE(I-1))
        W(2,I)=DSIN(PX*DBLE(I-1))
   10 CONTINUE
      RETURN
      END
      SUBROUTINE SETTBL2(W,N1,N2)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION W(2,N1,*)
!
      PI2=8.0D0*DATAN(1.0D0)
      PX=-PI2/(DBLE(N1)*DBLE(N2))
      DO 20 K=1,N2
!NODIR$ VECTOR ALIGNED
        DO 10 J=1,N1
          W(1,J,K)=DCOS(PX*DBLE(J-1)*DBLE(K-1))
          W(2,J,K)=DSIN(PX*DBLE(J-1)*DBLE(K-1))
   10   CONTINUE
   20 CONTINUE
      RETURN
      END
      SUBROUTINE FACTOR(N,IP)
      DIMENSION IP(*)
!
      IP(1)=0
      IP(2)=0
      IP(3)=0
      N2=N
      IF (MOD(N,2) .NE. 0 .AND. MOD(N,3) .NE. 0 .AND.   &
          MOD(N,5) .NE. 0) RETURN
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
      SUBROUTINE FACTOR8(N,IP)
      DIMENSION IP(*)
      INTEGER*8 N,N2
!
      IP(1)=0
      IP(2)=0
      IP(3)=0
      N2=N
      IF (MOD(N,2) .NE. 0 .AND. MOD(N,3) .NE. 0 .AND.  &
          MOD(N,5) .NE. 0) RETURN
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

!
!     FFTE: A FAST FOURIER TRANSFORM PACKAGE
!
!     (C) COPYRIGHT SOFTWARE, 2000-2004, ALL RIGHTS RESERVED
!                BY
!         DAISUKE TAKAHASHI
!         GRADUATE SCHOOL OF SYSTEMS AND INFORMATION ENGINEERING
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
        J=J+L
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
        J=J+L
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
        J=J+L
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
        J=J+L
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
        J=J+L
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
        J=J+L
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
        J=J+L
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
        J=J+L
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
      SUBROUTINE ZTRANS(A,B,N1,N2)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*)
      DIMENSION IP1(3),IP2(3)
!
      CALL FACTOR(N1,IP1)
      CALL FACTOR(N2,IP2)
!
      IF (N1 .EQ. 1 .OR. N2 .EQ. 1) THEN
        DO 10 I=1,N1*N2
          B(I)=A(I)
   10   CONTINUE
        RETURN
      END IF
!
      IF (IP1(1)+IP2(1) .LE. 1) THEN
        CALL ZTRANSA(A,B,N1,N2)
      ELSE
        CALL ZTRANSB(A,B,N1,N2)
      END IF
      RETURN
      END
      SUBROUTINE ZTRANSA(A,B,N1,N2)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(N1,*),B(N2,*)
!
      DO 20 I=1,N1
        DO 10 J=1,N2
          B(J,I)=A(I,J)
   10   CONTINUE
   20 CONTINUE
      RETURN
      END
      SUBROUTINE ZTRANSB(A,B,N1,N2)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(N1,*),B(N2,*)
!
      IF (N2 .GE. N1) THEN
        DO 20 I=0,N1-1
          DO 10 J=1,N1-I
            B(J,I+J)=A(I+J,J)
   10     CONTINUE
   20   CONTINUE
        DO 40 I=1,N2-N1
          DO 30 J=1,N1
            B(I+J,J)=A(J,I+J)
   30     CONTINUE
   40   CONTINUE
        DO 60 I=N2-N1+1,N2-1
          DO 50 J=1,N2-I
            B(I+J,J)=A(J,I+J)
   50     CONTINUE
   60   CONTINUE
      ELSE
        DO 80 I=0,N2-1
          DO 70 J=1,N2-I
            B(I+J,J)=A(J,I+J)
   70     CONTINUE
   80   CONTINUE
        DO 100 I=1,N1-N2
          DO 90 J=1,N2
            B(J,I+J)=A(I+J,J)
   90     CONTINUE
  100   CONTINUE
        DO 120 I=N1-N2+1,N1-1
          DO 110 J=1,N1-I
            B(J,I+J)=A(I+J,J)
  110     CONTINUE
  120   CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE ZTRANSMUL(A,B,W,N1,N2)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*),W(*)
      DIMENSION IP1(3),IP2(3)
!
      CALL FACTOR(N1,IP1)
      CALL FACTOR(N2,IP2)
!
      IF (N1 .EQ. 1 .OR. N2 .EQ. 1) THEN
        DO 10 I=1,N1*N2
          B(I)=A(I)*W(I)
   10   CONTINUE
        RETURN
      END IF
!
      IF (IP1(1)+IP2(1) .LE. 1) THEN
        CALL ZTRANSMULA(A,B,W,N1,N2)
      ELSE
        CALL ZTRANSMULB(A,B,W,N1,N2)
      END IF
      RETURN
      END
      SUBROUTINE ZTRANSMULA(A,B,W,N1,N2)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(N1,*),B(N2,*),W(N2,*)
!
      DO 20 I=1,N1
        DO 10 J=1,N2
          B(J,I)=A(I,J)*W(J,I)
   10   CONTINUE
   20 CONTINUE
      RETURN
      END
      SUBROUTINE ZTRANSMULB(A,B,W,N1,N2)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(N1,*),B(N2,*),W(N2,*)
!
      IF (N2 .GE. N1) THEN
        DO 20 I=0,N1-1
          DO 10 J=1,N1-I
            B(J,I+J)=A(I+J,J)*W(J,I+J)
   10     CONTINUE
   20   CONTINUE
        DO 40 I=1,N2-N1
          DO 30 J=1,N1
            B(I+J,J)=A(J,I+J)*W(I+J,J)
   30     CONTINUE
   40   CONTINUE
        DO 60 I=N2-N1+1,N2-1
          DO 50 J=1,N2-I
            B(I+J,J)=A(J,I+J)*W(I+J,J)
   50     CONTINUE
   60   CONTINUE
      ELSE
        DO 80 I=0,N2-1
          DO 70 J=1,N2-I
            B(I+J,J)=A(J,I+J)*W(I+J,J)
   70     CONTINUE
   80   CONTINUE
        DO 100 I=1,N1-N2
          DO 90 J=1,N2
            B(J,I+J)=A(I+J,J)*W(J,I+J)
   90     CONTINUE
  100   CONTINUE
        DO 120 I=N1-N2+1,N1-1
          DO 110 J=1,N1-I
            B(J,I+J)=A(I+J,J)*W(J,I+J)
  110     CONTINUE
  120   CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE MZTRANSA(A,B,NS,NY,NZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(NS,NY,*),B(NS,NZ,*)
!
      IF (NS .EQ. 1) THEN
        CALL ZTRANS(A(1,1,1),B(1,1,1),NY,NZ)
      ELSE
        DO 30 J=1,NY
          DO 20 K=1,NZ
            DO 10 I=1,NS
              B(I,K,J)=A(I,J,K)
   10       CONTINUE
   20     CONTINUE
   30   CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE MZTRANSB(A,B,NX,NY,NS)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(NX,NY,*),B(NY,NX,*)
!
      DO 10 I=1,NS
        CALL ZTRANS(A(1,1,I),B(1,1,I),NX,NY)
   10 CONTINUE
      RETURN
      END

!
!     FFTE: A FAST FOURIER TRANSFORM PACKAGE
!
!     (C) COPYRIGHT SOFTWARE, 2000-2004, ALL RIGHTS RESERVED
!                BY
!         DAISUKE TAKAHASHI
!         GRADUATE SCHOOL OF SYSTEMS AND INFORMATION ENGINEERING
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
      DIMENSION A(2,M,*),B(2,M,*)
!
!NODIR$ VECTOR ALIGNED
      DO 10 I=1,M
        X0=A(1,I,1)
        Y0=A(2,I,1)
        X1=A(1,I,2)
        Y1=A(2,I,2)
        B(1,I,1)=X0+X1
        B(2,I,1)=Y0+Y1
        B(1,I,2)=X0-X1
        B(2,I,2)=Y0-Y1
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FFT3A(A,B,W,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,L,*),B(2,3,*),W(2,*)
      DATA C31/0.86602540378443865D0/C32/0.5D0/
!
!NODIR$ VECTOR ALIGNED
      DO 10 J=1,L
        WR1=W(1,J)
        WI1=W(2,J)
        WR2=WR1*WR1-WI1*WI1
        WI2=WR1*WI1+WR1*WI1
        X0=A(1,J,2)+A(1,J,3)
        Y0=A(2,J,2)+A(2,J,3)
        X1=A(1,J,1)-C32*X0
        Y1=A(2,J,1)-C32*Y0
        X2=C31*(A(2,J,2)-A(2,J,3))
        Y2=C31*(A(1,J,3)-A(1,J,2))
        B(1,1,J)=A(1,J,1)+X0
        B(2,1,J)=A(2,J,1)+Y0
        B(1,2,J)=WR1*(X1+X2)-WI1*(Y1+Y2)
        B(2,2,J)=WR1*(Y1+Y2)+WI1*(X1+X2)
        B(1,3,J)=WR2*(X1-X2)-WI2*(Y1-Y2)
        B(2,3,J)=WR2*(Y1-Y2)+WI2*(X1-X2)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FFT3B(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,M,L,*),B(2,M,3,*),W(2,*)
      DATA C31/0.86602540378443865D0/C32/0.5D0/
!
!NODIR$ VECTOR ALIGNED
      DO 10 I=1,M
        X0=A(1,I,1,2)+A(1,I,1,3)
        Y0=A(2,I,1,2)+A(2,I,1,3)
        X1=A(1,I,1,1)-C32*X0
        Y1=A(2,I,1,1)-C32*Y0
        X2=C31*(A(2,I,1,2)-A(2,I,1,3))
        Y2=C31*(A(1,I,1,3)-A(1,I,1,2))
        B(1,I,1,1)=A(1,I,1,1)+X0
        B(2,I,1,1)=A(2,I,1,1)+Y0
        B(1,I,2,1)=X1+X2
        B(2,I,2,1)=Y1+Y2
        B(1,I,3,1)=X1-X2
        B(2,I,3,1)=Y1-Y2
   10 CONTINUE
      DO 30 J=2,L
        WR1=W(1,J)
        WI1=W(2,J)
        WR2=WR1*WR1-WI1*WI1
        WI2=WR1*WI1+WR1*WI1
!NODIR$ VECTOR ALIGNED
        DO 20 I=1,M
          X0=A(1,I,J,2)+A(1,I,J,3)
          Y0=A(2,I,J,2)+A(2,I,J,3)
          X1=A(1,I,J,1)-C32*X0
          Y1=A(2,I,J,1)-C32*Y0
          X2=C31*(A(2,I,J,2)-A(2,I,J,3))
          Y2=C31*(A(1,I,J,3)-A(1,I,J,2))
          B(1,I,1,J)=A(1,I,J,1)+X0
          B(2,I,1,J)=A(2,I,J,1)+Y0
          B(1,I,2,J)=WR1*(X1+X2)-WI1*(Y1+Y2)
          B(2,I,2,J)=WR1*(Y1+Y2)+WI1*(X1+X2)
          B(1,I,3,J)=WR2*(X1-X2)-WI2*(Y1-Y2)
          B(2,I,3,J)=WR2*(Y1-Y2)+WI2*(X1-X2)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
      SUBROUTINE FFT4A(A,B,W,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,L,*),B(2,4,*),W(2,*)
!
!NODIR$ VECTOR ALIGNED
      DO 10 J=1,L
        WR1=W(1,J)
        WI1=W(2,J)
        WR2=WR1*WR1-WI1*WI1
        WI2=WR1*WI1+WR1*WI1
        WR3=WR1*WR2-WI1*WI2
        WI3=WR1*WI2+WI1*WR2
        X0=A(1,J,1)+A(1,J,3)
        Y0=A(2,J,1)+A(2,J,3)
        X1=A(1,J,1)-A(1,J,3)
        Y1=A(2,J,1)-A(2,J,3)
        X2=A(1,J,2)+A(1,J,4)
        Y2=A(2,J,2)+A(2,J,4)
        X3=A(2,J,2)-A(2,J,4)
        Y3=A(1,J,4)-A(1,J,2)
        B(1,1,J)=X0+X2
        B(2,1,J)=Y0+Y2
        B(1,3,J)=WR2*(X0-X2)-WI2*(Y0-Y2)
        B(2,3,J)=WR2*(Y0-Y2)+WI2*(X0-X2)
        B(1,2,J)=WR1*(X1+X3)-WI1*(Y1+Y3)
        B(2,2,J)=WR1*(Y1+Y3)+WI1*(X1+X3)
        B(1,4,J)=WR3*(X1-X3)-WI3*(Y1-Y3)
        B(2,4,J)=WR3*(Y1-Y3)+WI3*(X1-X3)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FFT4B(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,M,L,*),B(2,M,4,*),W(2,*)
!
!NODIR$ VECTOR ALIGNED
      DO 10 I=1,M
        X0=A(1,I,1,1)+A(1,I,1,3)
        Y0=A(2,I,1,1)+A(2,I,1,3)
        X1=A(1,I,1,1)-A(1,I,1,3)
        Y1=A(2,I,1,1)-A(2,I,1,3)
        X2=A(1,I,1,2)+A(1,I,1,4)
        Y2=A(2,I,1,2)+A(2,I,1,4)
        X3=A(2,I,1,2)-A(2,I,1,4)
        Y3=A(1,I,1,4)-A(1,I,1,2)
        B(1,I,1,1)=X0+X2
        B(2,I,1,1)=Y0+Y2
        B(1,I,3,1)=X0-X2
        B(2,I,3,1)=Y0-Y2
        B(1,I,2,1)=X1+X3
        B(2,I,2,1)=Y1+Y3
        B(1,I,4,1)=X1-X3
        B(2,I,4,1)=Y1-Y3
   10 CONTINUE
      DO 30 J=2,L
        WR1=W(1,J)
        WI1=W(2,J)
        WR2=WR1*WR1-WI1*WI1
        WI2=WR1*WI1+WR1*WI1
        WR3=WR1*WR2-WI1*WI2
        WI3=WR1*WI2+WI1*WR2
!NODIR$ VECTOR ALIGNED
        DO 20 I=1,M
          X0=A(1,I,J,1)+A(1,I,J,3)
          Y0=A(2,I,J,1)+A(2,I,J,3)
          X1=A(1,I,J,1)-A(1,I,J,3)
          Y1=A(2,I,J,1)-A(2,I,J,3)
          X2=A(1,I,J,2)+A(1,I,J,4)
          Y2=A(2,I,J,2)+A(2,I,J,4)
          X3=A(2,I,J,2)-A(2,I,J,4)
          Y3=A(1,I,J,4)-A(1,I,J,2)
          B(1,I,1,J)=X0+X2
          B(2,I,1,J)=Y0+Y2
          B(1,I,3,J)=WR2*(X0-X2)-WI2*(Y0-Y2)
          B(2,I,3,J)=WR2*(Y0-Y2)+WI2*(X0-X2)
          B(1,I,2,J)=WR1*(X1+X3)-WI1*(Y1+Y3)
          B(2,I,2,J)=WR1*(Y1+Y3)+WI1*(X1+X3)
          B(1,I,4,J)=WR3*(X1-X3)-WI3*(Y1-Y3)
          B(2,I,4,J)=WR3*(Y1-Y3)+WI3*(X1-X3)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
      SUBROUTINE FFT5A(A,B,W,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,L,*),B(2,5,*),W(2,*)
      DATA C51/0.95105651629515357D0/C52/0.61803398874989485D0/   &
           C53/0.55901699437494742D0/C54/0.25D0/
!
!NODIR$ VECTOR ALIGNED
      DO 10 J=1,L
        WR1=W(1,J)
        WI1=W(2,J)
        WR2=WR1*WR1-WI1*WI1
        WI2=WR1*WI1+WR1*WI1
        WR3=WR1*WR2-WI1*WI2
        WI3=WR1*WI2+WI1*WR2
        WR4=WR2*WR2-WI2*WI2
        WI4=WR2*WI2+WR2*WI2
        X0=A(1,J,2)+A(1,J,5)
        Y0=A(2,J,2)+A(2,J,5)
        X1=A(1,J,3)+A(1,J,4)
        Y1=A(2,J,3)+A(2,J,4)
        X2=C51*(A(1,J,2)-A(1,J,5))
        Y2=C51*(A(2,J,2)-A(2,J,5))
        X3=C51*(A(1,J,3)-A(1,J,4))
        Y3=C51*(A(2,J,3)-A(2,J,4))
        X4=X0+X1
        Y4=Y0+Y1
        X5=C53*(X0-X1)
        Y5=C53*(Y0-Y1)
        X6=A(1,J,1)-C54*X4
        Y6=A(2,J,1)-C54*Y4
        X7=X6+X5
        Y7=Y6+Y5
        X8=X6-X5
        Y8=Y6-Y5
        X9=Y2+C52*Y3
        Y9=-X2-C52*X3
        X10=C52*Y2-Y3
        Y10=X3-C52*X2
        B(1,1,J)=A(1,J,1)+X4
        B(2,1,J)=A(2,J,1)+Y4
        B(1,2,J)=WR1*(X7+X9)-WI1*(Y7+Y9)
        B(2,2,J)=WR1*(Y7+Y9)+WI1*(X7+X9)
        B(1,3,J)=WR2*(X8+X10)-WI2*(Y8+Y10)
        B(2,3,J)=WR2*(Y8+Y10)+WI2*(X8+X10)
        B(1,4,J)=WR3*(X8-X10)-WI3*(Y8-Y10)
        B(2,4,J)=WR3*(Y8-Y10)+WI3*(X8-X10)
        B(1,5,J)=WR4*(X7-X9)-WI4*(Y7-Y9)
        B(2,5,J)=WR4*(Y7-Y9)+WI4*(X7-X9)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FFT5B(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,M,L,*),B(2,M,5,*),W(2,*)
      DATA C51/0.95105651629515357D0/C52/0.61803398874989485D0/   &
           C53/0.55901699437494742D0/C54/0.25D0/
!
!NODIR$ VECTOR ALIGNED
      DO 10 I=1,M
        X0=A(1,I,1,2)+A(1,I,1,5)
        Y0=A(2,I,1,2)+A(2,I,1,5)
        X1=A(1,I,1,3)+A(1,I,1,4)
        Y1=A(2,I,1,3)+A(2,I,1,4)
        X2=C51*(A(1,I,1,2)-A(1,I,1,5))
        Y2=C51*(A(2,I,1,2)-A(2,I,1,5))
        X3=C51*(A(1,I,1,3)-A(1,I,1,4))
        Y3=C51*(A(2,I,1,3)-A(2,I,1,4))
        X4=X0+X1
        Y4=Y0+Y1
        X5=C53*(X0-X1)
        Y5=C53*(Y0-Y1)
        X6=A(1,I,1,1)-C54*X4
        Y6=A(2,I,1,1)-C54*Y4
        X7=X6+X5
        Y7=Y6+Y5
        X8=X6-X5
        Y8=Y6-Y5
        X9=Y2+C52*Y3
        Y9=-X2-C52*X3
        X10=C52*Y2-Y3
        Y10=X3-C52*X2
        B(1,I,1,1)=A(1,I,1,1)+X4
        B(2,I,1,1)=A(2,I,1,1)+Y4
        B(1,I,2,1)=X7+X9
        B(2,I,2,1)=Y7+Y9
        B(1,I,3,1)=X8+X10
        B(2,I,3,1)=Y8+Y10
        B(1,I,4,1)=X8-X10
        B(2,I,4,1)=Y8-Y10
        B(1,I,5,1)=X7-X9
        B(2,I,5,1)=Y7-Y9
   10 CONTINUE
      DO 30 J=2,L
        WR1=W(1,J)
        WI1=W(2,J)
        WR2=WR1*WR1-WI1*WI1
        WI2=WR1*WI1+WR1*WI1
        WR3=WR1*WR2-WI1*WI2
        WI3=WR1*WI2+WI1*WR2
        WR4=WR2*WR2-WI2*WI2
        WI4=WR2*WI2+WR2*WI2
!NODIR$ VECTOR ALIGNED
        DO 20 I=1,M
          X0=A(1,I,J,2)+A(1,I,J,5)
          Y0=A(2,I,J,2)+A(2,I,J,5)
          X1=A(1,I,J,3)+A(1,I,J,4)
          Y1=A(2,I,J,3)+A(2,I,J,4)
          X2=C51*(A(1,I,J,2)-A(1,I,J,5))
          Y2=C51*(A(2,I,J,2)-A(2,I,J,5))
          X3=C51*(A(1,I,J,3)-A(1,I,J,4))
          Y3=C51*(A(2,I,J,3)-A(2,I,J,4))
          X4=X0+X1
          Y4=Y0+Y1
          X5=C53*(X0-X1)
          Y5=C53*(Y0-Y1)
          X6=A(1,I,J,1)-C54*X4
          Y6=A(2,I,J,1)-C54*Y4
          X7=X6+X5
          Y7=Y6+Y5
          X8=X6-X5
          Y8=Y6-Y5
          X9=Y2+C52*Y3
          Y9=-X2-C52*X3
          X10=C52*Y2-Y3
          Y10=X3-C52*X2
          B(1,I,1,J)=A(1,I,J,1)+X4
          B(2,I,1,J)=A(2,I,J,1)+Y4
          B(1,I,2,J)=WR1*(X7+X9)-WI1*(Y7+Y9)
          B(2,I,2,J)=WR1*(Y7+Y9)+WI1*(X7+X9)
          B(1,I,3,J)=WR2*(X8+X10)-WI2*(Y8+Y10)
          B(2,I,3,J)=WR2*(Y8+Y10)+WI2*(X8+X10)
          B(1,I,4,J)=WR3*(X8-X10)-WI3*(Y8-Y10)
          B(2,I,4,J)=WR3*(Y8-Y10)+WI3*(X8-X10)
          B(1,I,5,J)=WR4*(X7-X9)-WI4*(Y7-Y9)
          B(2,I,5,J)=WR4*(Y7-Y9)+WI4*(X7-X9)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
      SUBROUTINE FFT8A(A,B,W,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,L,*),B(2,8,*),W(2,*)
      DATA C81/0.70710678118654752D0/
!
!NODIR$ VECTOR ALIGNED
      DO 10 J=1,L
        WR1=W(1,J)
        WI1=W(2,J)
        WR2=WR1*WR1-WI1*WI1
        WI2=WR1*WI1+WR1*WI1
        WR3=WR1*WR2-WI1*WI2
        WI3=WR1*WI2+WI1*WR2
        WR4=WR2*WR2-WI2*WI2
        WI4=WR2*WI2+WR2*WI2
        WR5=WR2*WR3-WI2*WI3
        WI5=WR2*WI3+WI2*WR3
        WR6=WR3*WR3-WI3*WI3
        WI6=WR3*WI3+WR3*WI3
        WR7=WR3*WR4-WI3*WI4
        WI7=WR3*WI4+WI3*WR4
        X0=A(1,J,1)+A(1,J,5)
        Y0=A(2,J,1)+A(2,J,5)
        X1=A(1,J,1)-A(1,J,5)
        Y1=A(2,J,1)-A(2,J,5)
        X2=A(1,J,3)+A(1,J,7)
        Y2=A(2,J,3)+A(2,J,7)
        X3=A(2,J,3)-A(2,J,7)
        Y3=A(1,J,7)-A(1,J,3)
        U0=X0+X2
        V0=Y0+Y2
        U1=X0-X2
        V1=Y0-Y2
        X4=A(1,J,2)+A(1,J,6)
        Y4=A(2,J,2)+A(2,J,6)
        X5=A(1,J,2)-A(1,J,6)
        Y5=A(2,J,2)-A(2,J,6)
        X6=A(1,J,4)+A(1,J,8)
        Y6=A(2,J,4)+A(2,J,8)
        X7=A(1,J,4)-A(1,J,8)
        Y7=A(2,J,4)-A(2,J,8)
        U2=X4+X6
        V2=Y4+Y6
        U3=Y4-Y6
        V3=X6-X4
        B(1,1,J)=U0+U2
        B(2,1,J)=V0+V2
        B(1,5,J)=WR4*(U0-U2)-WI4*(V0-V2)
        B(2,5,J)=WR4*(V0-V2)+WI4*(U0-U2)
        B(1,3,J)=WR2*(U1+U3)-WI2*(V1+V3)
        B(2,3,J)=WR2*(V1+V3)+WI2*(U1+U3)
        B(1,7,J)=WR6*(U1-U3)-WI6*(V1-V3)
        B(2,7,J)=WR6*(V1-V3)+WI6*(U1-U3)
        U0=X1+C81*(X5-X7)
        V0=Y1+C81*(Y5-Y7)
        U1=X1-C81*(X5-X7)
        V1=Y1-C81*(Y5-Y7)
        U2=X3+C81*(Y5+Y7)
        V2=Y3-C81*(X5+X7)
        U3=X3-C81*(Y5+Y7)
        V3=Y3+C81*(X5+X7)
        B(1,2,J)=WR1*(U0+U2)-WI1*(V0+V2)
        B(2,2,J)=WR1*(V0+V2)+WI1*(U0+U2)
        B(1,6,J)=WR5*(U1+U3)-WI5*(V1+V3)
        B(2,6,J)=WR5*(V1+V3)+WI5*(U1+U3)
        B(1,4,J)=WR3*(U1-U3)-WI3*(V1-V3)
        B(2,4,J)=WR3*(V1-V3)+WI3*(U1-U3)
        B(1,8,J)=WR7*(U0-U2)-WI7*(V0-V2)
        B(2,8,J)=WR7*(V0-V2)+WI7*(U0-U2)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FFT8B(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,M,L,*),B(2,M,8,*),W(2,*)
      DATA C81/0.70710678118654752D0/
!
!NODIR$ VECTOR ALIGNED
      DO 10 I=1,M
        X0=A(1,I,1,1)+A(1,I,1,5)
        Y0=A(2,I,1,1)+A(2,I,1,5)
        X1=A(1,I,1,1)-A(1,I,1,5)
        Y1=A(2,I,1,1)-A(2,I,1,5)
        X2=A(1,I,1,3)+A(1,I,1,7)
        Y2=A(2,I,1,3)+A(2,I,1,7)
        X3=A(2,I,1,3)-A(2,I,1,7)
        Y3=A(1,I,1,7)-A(1,I,1,3)
        U0=X0+X2
        V0=Y0+Y2
        U1=X0-X2
        V1=Y0-Y2
        X4=A(1,I,1,2)+A(1,I,1,6)
        Y4=A(2,I,1,2)+A(2,I,1,6)
        X5=A(1,I,1,2)-A(1,I,1,6)
        Y5=A(2,I,1,2)-A(2,I,1,6)
        X6=A(1,I,1,4)+A(1,I,1,8)
        Y6=A(2,I,1,4)+A(2,I,1,8)
        X7=A(1,I,1,4)-A(1,I,1,8)
        Y7=A(2,I,1,4)-A(2,I,1,8)
        U2=X4+X6
        V2=Y4+Y6
        U3=Y4-Y6
        V3=X6-X4
        B(1,I,1,1)=U0+U2
        B(2,I,1,1)=V0+V2
        B(1,I,5,1)=U0-U2
        B(2,I,5,1)=V0-V2
        B(1,I,3,1)=U1+U3
        B(2,I,3,1)=V1+V3
        B(1,I,7,1)=U1-U3
        B(2,I,7,1)=V1-V3
        U0=X1+C81*(X5-X7)
        V0=Y1+C81*(Y5-Y7)
        U1=X1-C81*(X5-X7)
        V1=Y1-C81*(Y5-Y7)
        U2=X3+C81*(Y5+Y7)
        V2=Y3-C81*(X5+X7)
        U3=X3-C81*(Y5+Y7)
        V3=Y3+C81*(X5+X7)
        B(1,I,2,1)=U0+U2
        B(2,I,2,1)=V0+V2
        B(1,I,6,1)=U1+U3
        B(2,I,6,1)=V1+V3
        B(1,I,4,1)=U1-U3
        B(2,I,4,1)=V1-V3
        B(1,I,8,1)=U0-U2
        B(2,I,8,1)=V0-V2
   10 CONTINUE
      DO 30 J=2,L
        WR1=W(1,J)
        WI1=W(2,J)
        WR2=WR1*WR1-WI1*WI1
        WI2=WR1*WI1+WR1*WI1
        WR3=WR1*WR2-WI1*WI2
        WI3=WR1*WI2+WI1*WR2
        WR4=WR2*WR2-WI2*WI2
        WI4=WR2*WI2+WR2*WI2
        WR5=WR2*WR3-WI2*WI3
        WI5=WR2*WI3+WI2*WR3
        WR6=WR3*WR3-WI3*WI3
        WI6=WR3*WI3+WR3*WI3
        WR7=WR3*WR4-WI3*WI4
        WI7=WR3*WI4+WI3*WR4
!NODIR$ VECTOR ALIGNED
        DO 20 I=1,M
          X0=A(1,I,J,1)+A(1,I,J,5)
          Y0=A(2,I,J,1)+A(2,I,J,5)
          X1=A(1,I,J,1)-A(1,I,J,5)
          Y1=A(2,I,J,1)-A(2,I,J,5)
          X2=A(1,I,J,3)+A(1,I,J,7)
          Y2=A(2,I,J,3)+A(2,I,J,7)
          X3=A(2,I,J,3)-A(2,I,J,7)
          Y3=A(1,I,J,7)-A(1,I,J,3)
          U0=X0+X2
          V0=Y0+Y2
          U1=X0-X2
          V1=Y0-Y2
          X4=A(1,I,J,2)+A(1,I,J,6)
          Y4=A(2,I,J,2)+A(2,I,J,6)
          X5=A(1,I,J,2)-A(1,I,J,6)
          Y5=A(2,I,J,2)-A(2,I,J,6)
          X6=A(1,I,J,4)+A(1,I,J,8)
          Y6=A(2,I,J,4)+A(2,I,J,8)
          X7=A(1,I,J,4)-A(1,I,J,8)
          Y7=A(2,I,J,4)-A(2,I,J,8)
          U2=X4+X6
          V2=Y4+Y6
          U3=Y4-Y6
          V3=X6-X4
          B(1,I,1,J)=U0+U2
          B(2,I,1,J)=V0+V2
          B(1,I,5,J)=WR4*(U0-U2)-WI4*(V0-V2)
          B(2,I,5,J)=WR4*(V0-V2)+WI4*(U0-U2)
          B(1,I,3,J)=WR2*(U1+U3)-WI2*(V1+V3)
          B(2,I,3,J)=WR2*(V1+V3)+WI2*(U1+U3)
          B(1,I,7,J)=WR6*(U1-U3)-WI6*(V1-V3)
          B(2,I,7,J)=WR6*(V1-V3)+WI6*(U1-U3)
          U0=X1+C81*(X5-X7)
          V0=Y1+C81*(Y5-Y7)
          U1=X1-C81*(X5-X7)
          V1=Y1-C81*(Y5-Y7)
          U2=X3+C81*(Y5+Y7)
          V2=Y3-C81*(X5+X7)
          U3=X3-C81*(Y5+Y7)
          V3=Y3+C81*(X5+X7)
          B(1,I,2,J)=WR1*(U0+U2)-WI1*(V0+V2)
          B(2,I,2,J)=WR1*(V0+V2)+WI1*(U0+U2)
          B(1,I,6,J)=WR5*(U1+U3)-WI5*(V1+V3)
          B(2,I,6,J)=WR5*(V1+V3)+WI5*(U1+U3)
          B(1,I,4,J)=WR3*(U1-U3)-WI3*(V1-V3)
          B(2,I,4,J)=WR3*(V1-V3)+WI3*(U1-U3)
          B(1,I,8,J)=WR7*(U0-U2)-WI7*(V0-V2)
          B(2,I,8,J)=WR7*(V0-V2)+WI7*(U0-U2)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
