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
