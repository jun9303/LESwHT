!=======================================================================
!
!     MAIN SOLVER OF POISSON EQUATION
!
!     AX=b
!     A: coefficient of discretized poisson equation
!     X: PHI: Pseudo pressure, Output of subroutine POISSON
!     b: DIVGSUM: OUTPUT of subroutine DIVGS.
!
!     x-direction: Fourier transform
!     y-direction: TDMA
!     z-direction: MULTI-GRID iteration/GSOR method
!
!     Jun. 2017, J. Park   
!
!=======================================================================
      SUBROUTINE POISINIT
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER*8   :: I,J,K,ILEV

       CALL X_FT_ALLO

       ALLOCATE (AC (N3MD,N2M,N1MH))
       ALLOCATE (GAM(N3MD,N2M,N1MH))
       ALLOCATE (BET(N3MD,N2M,N1MH))
 
       CALL PMAT
!------MULTIGRID METHOD
       CALL COEFMG

       DO ILEV= 0,NLEV
        DO I= 1,N1MH
         DO K=KKMG(ILEV,1),KKMG(ILEV,2)
          BET(K,1,I)= 1./AC(K,1,I)
         ENDDO
         DO J= 2,N2M
         DO K= KKMG(ILEV,1),KKMG(ILEV,2)
          GAM(K,J,I)= AN(J-1)*BET(K,J-1,I)
          BET(K,J,I)= 1./(AC(K,J,I)-AS(J)*GAM(K,J,I))
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
       INTEGER*8   :: I,J,K,JP,KP

!-----FOR TOP LEVEL
      DO K= 1,N3M
        KP=KPV(K)
        AB(K)=(1.-FIXKL(K))*C2CZI(K )*F2FZI(K)
        AF(K)=(1.-FIXKU(K))*C2CZI(KP)*F2FZI(K)
      ENDDO

      DO J= 1,N2M
        JP=JPV(J)
        AS(J)=(1.-FIXJL(J))*C2CYI(J )*F2FYI(J)
        AN(J)=(1.-FIXJU(J))*C2CYI(JP)*F2FYI(J)
      ENDDO

      DO J= 1,N2M
      DO K= 1,N3M
      AC(K,J,1)= -1.*(AB(K)+AF(K)+AS(J)+AN(J))
      ENDDO
      ENDDO

      N1MH=N1M/2+1

      CALL MWAVENUMBER     ! INIT. MODIFIED WAVE #.

      DO 40 I= 2,N1MH
      DO 40 J= 1,N2M
      DO 40 K= 1,N3M
       AC(K,J,I)=AC(K,J,1)-AI3(I)
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
       INTEGER*8   :: I
       REAL*8      :: PI

       PI = ACOS(-1.)

      DO I= 1,N1MH
       AI3(I)= 2.*(1.-COS(2.*PI*FLOAT(I-1)/FLOAT(N1M)))*F2FXI(1)*F2FXI(1)
      ENDDO

      RETURN
      END

!=======================================================================
      SUBROUTINE POISSON(PHI,DIVGSUM)
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER*8   :: I,J,K,III
       REAL*8      :: PHI(0:N1,0:N2,0:N3),DIVGSUM(0:N1,0:N2,0:N3)
       COMPLEX*16  :: CCAP(N3,N2,N1MH)
       COMPLEX*16, dimension (:,:), allocatable :: XXX,CP
       COMPLEX*16, dimension (:),   allocatable :: XXXX,XXXX_B

       REAL*8      :: TEST,PHIREF

! --- DO THE FORWARD FFT
!$OMP PARALLEL &
!$OMP private(XXX,XXXX,XXXX_B)
      allocate(XXX(N1M,N3M))
      allocate(XXXX(N1M),XXXX_B(N1M*2))
      CALL ZFFT1D(XXXX,N1M, 0,XXXX_B)
!$OMP DO
      DO 100 J=1,N2M

      DO I=1,N1M
      DO K=1,N3M
       XXX(I,K)=DIVGSUM(I,J,K)
      ENDDO
      ENDDO


      DO K=1,N3M
       CALL ZFFT1D(XXX(1,K),N1M,-1,XXXX_B)
      END DO

      DO I=1,N1MH
      DO K=1,N3M
       CCAP(K,J,I)=XXX(I,K)
      ENDDO
      ENDDO
  100 CONTINUE
!$OMP END DO
      deallocate(XXX,XXXX,XXXX_B)
!$OMP END PARALLEL

! --- SOLVE A SET OF POISSON EQS.
       TEST=TEST1/FLOAT(N1MH)*FLOAT(N1M)*0.9

!$OMP PARALLEL &
!$OMP private(CP)
      allocate(CP(0:N3,0:N2))
!$OMP DO
      DO 200 III= 1,N1MH
        CP   = 0.

        IF (III.LE.IMGSOR) THEN
         CALL MG2D(CP,CCAP(1,1,III),III,TEST,0.)
        ELSE
         CALL GSOR2D(CP,CCAP(1,1,III),III,TEST,0.)
        ENDIF
        
        DO J=1,N2M
        DO K=1,N3M
        CCAP(K,J,III)=CP(K,J)
        ENDDO
        ENDDO
  200 CONTINUE
!$OMP END DO
      deallocate(CP)
!$OMP END PARALLEL
 

! --- DO THE INVERSE FFT
!$OMP PARALLEL &
!$OMP private(XXX,XXXX,XXXX_B)
      allocate(XXX(N1M,N3M))
      allocate(XXXX(N1M),XXXX_B(N1M*2))
      CALL ZFFT1D(XXXX,N1M, 0,XXXX_B)
!$OMP DO
      DO 300 J=1,N2M
        DO I=1,N1MH
        DO K=1,N3M
        XXX(I,K)= CCAP(K,J,I)
        ENDDO
        ENDDO
        
        DO I=N1MH+1,N1M
        DO K=1,N3M  
         XXX(I,K)= CONJG(CCAP(K,J,N1M+2-I))
        ENDDO
        ENDDO
        
        DO K=1,N3M
        CALL ZFFT1D(XXX(1,K),N1M, 1,XXXX_B)
        ENDDO

        DO I=1,N1M        
        DO K=1,N3M
        PHI(I,J,K)= REAL(XXX(I,K))
        ENDDO
        ENDDO
  300 CONTINUE
!$OMP END DO
      deallocate(XXX,XXXX,XXXX_B)
!$OMP END PARALLEL

  400 CONTINUE


      IF(ICH.EQ.1) THEN
!     SET THE AVERAGE PHI AT THE UPPER WALL TO BE ZERO.
      PHIREF = 0.
!$OMP PARALLEL DO &
!$OMP reduction(+:PHIREF)
      DO 80 I=1,N1M
      PHIREF= PHIREF+ PHI(I,N2M,1)*F2FX(I)
   80 CONTINUE
      PHIREF= PHIREF/XL

!$OMP PARALLEL DO
      DO 93 K=1,N3M
      DO 93 J=1,N2M
      DO 93 I=1,N1M
      PHI(I,J,K)=PHI(I,J,K)-PHIREF
   93 CONTINUE
      ENDIF

      RETURN
      END

!=======================================================================
       SUBROUTINE MG2D(PC,RHS,IV,TEST,OLDV)
!=======================================================================
!     IV      : WAVE NUMBER INDEX
!     MULTIGRID ENVIRONMENT VARIABLES
!     TEST    : CONDITION FOR CONVERGENCE
!     NLEV    : TOTAL CELS = (MINROW1+MINROW2)*(2**NLEV)
!     LEVHALF : HALF OF NLEV
!     IWC     : FLAG FOR USING W-CYCLE ( 1 IF YOU WANT TO USE W-CYCLE)
!     NBLI    : NUMBER OF BASE LEVEL ITERATION
!     MGITR   : NUMBER OF MAXIMUM ITERATIONS
!     IIMG(ILEV,1)  : START POINT OF INDEX AT ILEV LEVEL
!     IIMG(ILEV,2)  : END POINT OF INDEX AT ILEV LEVEL
!     IIMG(ILEV,3)  : NUMBER OF ROW AT ILEV LEVEL
!-----------------------------------------------------------------------
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER*8   :: K,J,IV,II,ILEV
       REAL*8      :: TEST,SUMRES,OLDV
       COMPLEX*16  :: PC(0:N3,0:N2),RHS(N3,N2)
       COMPLEX*16  :: RESD(N3MD,N2M),GGII(0:N3MD,0:N2)

       II= 0
       RESD = 0.
       GGII = 0.

       CALL TOPLEVEL(0,PC,RHS,SUMRES,OLDV,TEST,IV,RESD,GGII)
       IF(SUMRES .LT. TEST) GOTO 2000

       DO 1000 II= 1,MGITR            ! main iteration
        DO ILEV=NLEV-1,1,-1
         CALL RELAX(ILEV,0.,1,IV,RESD,GGII)
         CALL GODOWN(ILEV,IV,RESD,GGII)
        ENDDO

        CALL RELAX(0,0.,NBLI,IV,RESD,GGII)

        DO ILEV=0,NLEV-2
         CALL GOUP(ILEV,RESD,GGII)
         CALL RELAX(ILEV+1,1.,1,IV,RESD,GGII)
        ENDDO

        CALL TOPLEVEL(1,PC,RHS,SUMRES,1.,TEST,IV,RESD,GGII)

        IF(SUMRES .LT. TEST) GOTO 2000
 1000  CONTINUE
       WRITE(*,*) 'ITERATION LIMIT EXCEEDED.'

 2000  CONTINUE

        IF (IV .LE. 2) WRITE(77,201) ' MG, TIME = ',TIME,IV,II,SUMRES*DTCONST*FLOAT(N1MH)/FLOAT(N1M)

  201  FORMAT(A13,F13.5,' IV=',I5,' II=',I5,' RG=',ES16.8)

      RETURN
      END
!=======================================================================
      SUBROUTINE TOPLEVEL(ID,PC,RHS,SUMRES,OLDV,TEST,IV,RESD,GGII)    ! Zebra Version
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
        INTEGER*8   :: K,J,KC
        INTEGER*8   :: ID,IV
        REAL*8      :: SUMRES,TEST,OLDV
        COMPLEX*16  :: TT
        COMPLEX*16  :: PC(0:N3,0:N2),RHS(N3,N2)
        COMPLEX*16  :: RESD(N3MD,N2M),GGII(0:N3MD,0:N2)

      IF (ID .NE. 1) GOTO 25
!     INTERPOLATE & ADD
      KC= 0
      DO K=KKMG(NLEV-1,1),KKMG(NLEV-1,2)
        KC=KC+2
        DO J=1,N2M
          PC(KC,J)=PC(KC,J)+COI1(K)*GGII(K,J)+COI2(K)*GGII(KPM(K),J)
        ENDDO
      ENDDO
25    CONTINUE

!  RELAX
      DO J=1,N2M
      DO K=1,N3M,2
       GGII(K,J) = RHS(K,J)-OLDV*(AB(K)*PC(KMM(K),J)+AF(K)*PC(KPM(K),J))
      ENDDO
      ENDDO

      CALL TRDIAG1M(GGII,GGII,1,N3M,1,N2M,IV)

      IF(IV .EQ. 1) THEN
        TT= GGII(1,N2M)
        DO J=1,N2M
        DO K=1,N3M,2
          GGII(K,J)=GGII(K,J)-TT
        ENDDO
        ENDDO
      END IF

      DO J=1,N2M
      DO K=2,N3M,2
        GGII(K,J)=RHS(K,J)-(AB(K)*GGII(KMM(K),J)+AF(K)*(GGII(KPM(K),J)))
      ENDDO
      ENDDO

      CALL TRDIAG1M(GGII,GGII,2,N3M,1,N2M,IV)

!  CALCULATE RESIDUAL
      SUMRES= 0.

      DO 72 J=1,N2M
      DO 72 K=1,N3M
      RESD(K,J)=RHS(K,J)-AB(K)*GGII(KMM(K),J)-AF(K)*GGII(KPM(K),J)  &
                   -AS(J)*GGII(K,J-1)-AN(J)*GGII(K,J+1)       &
                   -AC(K,J,IV)*GGII(K,J)
      SUMRES= AMAX1( SUMRES,ABS(RESD(K,J)) )
      PC(K,J)=GGII(K,J)
72    CONTINUE

      IF(SUMRES .LT. TEST) GOTO 99
      IF(ID .EQ. 2) GOTO 99

!  RESTRICT
      KC=-1
      DO 101 K=KKMG(NLEV-1,1),KKMG(NLEV-1,2)
      KC=KC+2
      DO 101 J=1,N2M
      RESD(K,J)=RESD(KC,J)*COR1(K)+RESD(KC+1,J)*COR2(K)
101   CONTINUE

99    RETURN
      END
!=======================================================================
      SUBROUTINE TRDIAG1M(RR,UU,L1,L2,LL1,LL2,IV)
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
        INTEGER*8   :: K,J,L1,L2,LL1,LL2,IV
        COMPLEX*16  :: RR(0:N3MD,0:N2),UU(0:N3MD,0:N2)

      DO 10 K=L1,L2,2
        UU(K,LL1)=RR(K,LL1)*BET(K,1,IV)
10    CONTINUE

      DO 20 J=LL1+1,LL2
      DO 20 K=L1,L2,2
        UU(K,J)=(RR(K,J)-AS(J)*UU(K,J-1))*BET(K,J,IV)
20    CONTINUE
      DO 30 J=LL2-1,LL1,-1
      DO 30 K=L1,L2,2
        UU(K,J)=UU(K,J)-GAM(K,J+1,IV)*UU(K,J+1)
30    CONTINUE

      RETURN
      END

!=======================================================================
      SUBROUTINE GOUP(ILEV,RESD,GGII)! INTERPOLATE RESIDUAL & ADD IT TO HIGH LEVEL
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
        INTEGER*8   :: K,J,ILEV,KBGH
        COMPLEX*16  :: RESD(N3MD,N2M),GGII(0:N3MD,0:N2)

      KBGH=KKMG(ILEV+1,1)-1

      DO 21 K=KKMG(ILEV,1),KKMG(ILEV,2)
      KBGH=KBGH+2
      DO 21 J=1,N2M
      GGII(KBGH,J)=GGII(KBGH,J)+COI1(K)*GGII(K,J)+COI2(K)*GGII(KPM(K),J)
21    CONTINUE

      RETURN
      END

!=======================================================================
      SUBROUTINE GODOWN(ILEV,IV,RESD,GGII)        ! COMPUTE RESIDUAL & RESTRICT IT
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
        INTEGER*8   :: K,J,ILEV,IV,KBG
        COMPLEX*16  :: RESD(N3MD,N2M),GGII(0:N3MD,0:N2)

      DO 10 J=1,N2M
      DO 10 K=KKMG(ILEV,1),KKMG(ILEV,2)
      RESD(K,J)=RESD(K,J)-AB(K)*GGII(KMM(K),J)-AF(K)*GGII(KPM(K),J)   &
               -AS(J)*GGII(K,J-1)-AN(J)*GGII(K,J+1)-AC(K,J,IV)*GGII(K,J)
10    CONTINUE

      KBG=KKMG(ILEV,1)-2

      DO 21 K=KKMG(ILEV-1,1),KKMG(ILEV-1,2)
      KBG=KBG+2
      DO 21 J=1,N2M
      RESD(K,J)=RESD(KBG,J)*COR1(K)+RESD(KBG+1,J)*COR2(K)
21    CONTINUE

      RETURN
      END

!=======================================================================
      SUBROUTINE RELAX(ILEV,OLDV,IITER,IV,RESD,GGII)   ! Zebra Version
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
        INTEGER*8   :: K,J,ILEV,IITER,IV,KK
        COMPLEX*16  :: RESD(N3MD,N2M),GGII(0:N3MD,0:N2)
        REAL*8      :: OLDV

      DO 12 J=1,N2M
      DO 12 K=KKMG(ILEV,1),KKMG(ILEV,2),2
      GGII(K,J)=RESD(K,J)-OLDV*(AB(K)*GGII(KMM(K),J)+AF(K)*GGII(KPM(K),J))
12    CONTINUE

      CALL TRDIAG1M(GGII,GGII,KKMG(ILEV,1),KKMG(ILEV,2),1,N2M,IV)

      DO 32 J=1,N2M
      DO 32 K=KKMG(ILEV,1)+1,KKMG(ILEV,2),2
      GGII(K,J)=RESD(K,J)-AB(K)*GGII(KMM(K),J)-AF(K)*GGII(KPM(K),J)
32    CONTINUE

      CALL TRDIAG1M(GGII,GGII,KKMG(ILEV,1)+1,KKMG(ILEV,2),1,N2M,IV)

      DO 50 KK=1,IITER-1

      DO 102 J=1,N2M
      DO 102 K=KKMG(ILEV,1),KKMG(ILEV,2),2
      GGII(K,J)=RESD(K,J)-(AB(K)*GGII(KMM(K),J)+AF(K)*GGII(KPM(K),J))
102    CONTINUE

      CALL TRDIAG1M(GGII,GGII,KKMG(ILEV,1),KKMG(ILEV,2),1,N2M,IV)

      DO 302 J=1,N2M
      DO 302 K=KKMG(ILEV,1)+1,KKMG(ILEV,2),2
      GGII(K,J)=RESD(K,J)-AB(K)*GGII(KMM(K),J)-AF(K)*GGII(KPM(K),J)
302    CONTINUE

      CALL TRDIAG1M(GGII,GGII,KKMG(ILEV,1)+1,KKMG(ILEV,2),1,N2M,IV)

50    CONTINUE
      RETURN
      END
!=======================================================================
      SUBROUTINE GSOR2D(U,RHS,IV,TEST,OLDV)    ! 1 EQ. TYPE
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER*8    :: K,J,IV
       REAL*8       :: TEST,OLDV
       COMPLEX*16   :: U(0:N3,0:N2),RHS(N3,N2)
       COMPLEX*16   :: GGII(0:N3MD,0:N2)
       COMPLEX*16   :: TT
       INTEGER*8    :: II
       REAL*8       :: WW,WW2,ERRMAX

      GGII = 0.
      WW   = WWSOR
      WW2  = 1.-WW
      II   = 0

!  HALF RELAX
!  ----------
      DO J=1,N2M
      DO K=1,N3M,2
       GGII(K,J)=RHS(K,J)-OLDV*(AB(K)*U(KMM(K),J)+AF(K)*U(KPM(K),J))
      ENDDO
      ENDDO

      CALL TRDIAG1M(GGII,GGII,1,N3M,1,N2M,IV)

      DO J=1,N2M
      DO K=1,N3M,2
       U(K,J)=WW*GGII(K,J)+OLDV*WW2*U(K,J)
      ENDDO
      ENDDO

!  ANOTHER HALF
!  ------------
      IF(IV .EQ. 1) THEN
        TT=U(1,N2M)
        DO J=1,N2M
        DO K=1,N3M,2
          U(K,J)=U(K,J)-TT
        ENDDO
        ENDDO
      END IF

      DO J=1,N2M
      DO K=2,N3M,2
       GGII(K,J)=RHS(K,J)-(AB(K)*U(KMM(K),J)+AF(K)*U(KPM(K),J))
      ENDDO
      ENDDO

      CALL TRDIAG1M(GGII,GGII,2,N3M,1,N2M,IV)

      DO J=1,N2M
      DO K=2,N3M,2
       U(K,J)=WW*GGII(K,J)+OLDV*WW2*U(K,J)
      ENDDO
      ENDDO

      CALL RESID3(U,RHS,IV,ERRMAX)
      IF(ERRMAX .LT. TEST) GOTO 1000

!  MAIN ITERATION
!  ==============
      DO 100 II= 1,MGITR

!  HALF RELAX
!  ----------
      DO J=1,N2M
      DO K=1,N3M,2
       GGII(K,J)=RHS(K,J)-(AB(K)*U(KMM(K),J)+AF(K)*U(KPM(K),J))
      ENDDO
      ENDDO

      CALL TRDIAG1M(GGII,GGII,1,N3M,1,N2M,IV)

      DO J=1,N2M
      DO K=1,N3M,2
       U(K,J)=WW*GGII(K,J)+WW2*U(K,J)
      ENDDO
      ENDDO

!  ANOTHER HALF
!  ------------

       IF(IV .EQ. 1) THEN
          TT= U(1,N2M)
          DO J=1,N2M
          DO K=1,N3M,2
            U(K,J)=U(K,J)-TT
          ENDDO
          ENDDO
       END IF
 
      DO J=1,N2M
      DO K=2,N3M,2
       GGII(K,J)=RHS(K,J)-(AB(K)*U(KMM(K),J)+AF(K)*U(KPM(K),J))
      ENDDO
      ENDDO

      CALL TRDIAG1M(GGII,GGII,2,N3M,1,N2M,IV)

      DO J=1,N2M
      DO K=2,N3M,2
       U(K,J)=WW*GGII(K,J)+WW2*U(K,J)
      ENDDO
      ENDDO

      CALL RESID3(U,RHS,IV,ERRMAX)
      IF(ERRMAX .LT. TEST) GOTO 1000

  100 CONTINUE
      PRINT *,'ITERATION LIMIT EXCEEDED.'
      WRITE(77,201) 'SOR',IV,II,ERRMAX*DTCONST*FLOAT(N1MH)/FLOAT(N1M)
 1000 CONTINUE

201   FORMAT(A5,'  IV=',I5,'  II=',I5,'  RG=',ES23.15)

      RETURN
      END

!=======================================================================
      SUBROUTINE RESID3(U,RHS,IV,ERRMAX)
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER*8    :: K,J,IV
       REAL*8       :: ERRMAX
       COMPLEX*16   :: U(0:N3,0:N2),RHS(N3,N2)
       COMPLEX*16   :: ERR

       ERRMAX= 0.

      DO 72 J=1,N2M
      DO 72 K=1,N3M
       ERR=RHS(K,J)-AB(K)*U(KMM(K),J)-AF(K)*U(KPM(K),J)   &
                   -AS(J)*U(K,J-1)-AN(J)*U(K,J+1)   &
                   -AC(K,J,IV)*U(K,J)
       ERRMAX= AMAX1( ERRMAX,ABS(ERR) )
72    CONTINUE

      RETURN
      END
!=======================================================================
      SUBROUTINE COEFMG
!=======================================================================
       USE MOD_COMMON
       USE MOD_POISS
       IMPLICIT NONE
       INTEGER*8   :: I,J,K,ILEV
       INTEGER*8   :: MINROW,KBG,KEND,KSP,KC,KBGH
       REAL*8      :: KBZ(N3MD),KFZ(N3MD)
       REAL*8      :: ZMPM(0:N3MD)
       REAL*8      :: VDZ_KBG,VDZ_KEND,SDZ_KBG,SDZ_KEND


       MINROW=N3M/(2**NLEV)
       ZMPM = 0.

       LEVHALF= NINT(NLEV/2.)
       KKMG(NLEV,1)=1    ! START INDEX
       KKMG(NLEV,2)=N3M  ! END INDEXT
       KKMG(NLEV,3)=N3M  ! THE NUMBER OF POINTS AT THE ILEV

       DO K= 1,N3M
        ZMPM(K)= ZMP(K)
        KBZ(K)= 1-1/K                  ! 0 only if K=1
        KFZ(K)= 1-K/N3M                ! 0 only if K=N3M
       ENDDO


!  COMPUTE FOR LOWER LEVELS
      DO 100 ILEV=NLEV-1,0,-1

       KKMG(ILEV,1)=KKMG(ILEV+1,1)+KKMG(ILEV+1,3)
       KKMG(ILEV,3)=MINROW*(2**ILEV)
       KKMG(ILEV,2)=KKMG(ILEV,1)+KKMG(ILEV,3)-1

       KBG = KKMG(ILEV,1)
       KEND= KKMG(ILEV,2)

       KSP= 2**(NLEV-ILEV)           ! width of one cell at low level

       KC= 0
       DO K=KBG,KEND
        KC=KC+1
        ZMPM(K)= 0.5*(Z(KC*KSP+1)+Z((KC-1)*KSP+1))
        KBZ(K)= 1-KBG/K                  ! 0 onlyif K=KBG
        KFZ(K)= 1-K/KEND                ! 0 only if K=KEND
       ENDDO

       KC= 0
       DO K=KBG,KEND
        KC=KC+1
        AB(K)=KBZ(K)/((ZMPM(K)-ZMPM(K-1))*(Z(KC*KSP+1)-Z((KC-1)*KSP+1)))
        AF(K)=KFZ(K)/((ZMPM(K+1)-ZMPM(K))*(Z(KC*KSP+1)-Z((KC-1)*KSP+1)))
       ENDDO

       DO J= 1,N2M
       DO K= KKMG(ILEV,1),KKMG(ILEV,2)
        AC(K,J,1)=-1.*(AB(K)+AF(K)+AS(J)+AN(J))
       ENDDO
       ENDDO

       DO I= 2,N1MH
       DO J= 1,N2M
       DO K= KKMG(ILEV,1),KKMG(ILEV,2)
        AC(K,J,I)=AC(K,J,1)-AI3(I)
       ENDDO
       ENDDO
       ENDDO

 100  CONTINUE

!  CALCULATE RESTRICTION COEFFS
       DO 145 ILEV=NLEV,1,-1
        KBGH=KKMG(ILEV,1)
        DO 145 K=KKMG(ILEV-1,1),KKMG(ILEV-1,2)
        COR1(K)= (ZMPM(KBGH+1)-ZMPM(K))/(ZMPM(KBGH+1)-ZMPM(KBGH))
        COR2(K)= 1.-COR1(K)
        KBGH=KBGH+2
 145   CONTINUE

!  CALCULATE INTERPOLATION COEFFS
       DO 150 ILEV= 0,NLEV-1
        KBGH=KKMG(ILEV+1,1)+1
        DO 160 K=KKMG(ILEV,1),KKMG(ILEV,2)-1
        COI1(K)= (ZMPM(K+1)-ZMPM(KBGH))/(ZMPM(K+1)-ZMPM(K))  ! * lower value
        COI2(K)= 1.-COI1(K)
        KBGH=KBGH+2
 160   CONTINUE
       K=KKMG(ILEV,2)
       COI1(K)= 1.                 ! use only one lower point at upper wall
 150  CONTINUE


!===== FOR THE Z PERIODICIRY
!       INTRODUCE KPM & KMM
       IF (ZPRDIC.EQ.1) THEN
        DO ILEV=NLEV,0,-1
         KBG=KKMG(ILEV,1)
         KEND=KKMG(ILEV,2)
         DO K=KBG,KEND
         KPM(K)=K+1
         KMM(K)=K-1
         ENDDO
         KPM(KEND)=KBG
         KMM(KBG)=KEND
        ENDDO

        DO ILEV=NLEV-1,0,-1
         KBG=KKMG(ILEV,1)
         KEND=KKMG(ILEV,2)
         KSP= 2**(NLEV-ILEV)
         VDZ_KBG = ZMPM(KBG)-Z(1)+Z(N3)-ZMPM(KEND)
         VDZ_KEND= ZMPM(KBG)-Z(1)+Z(N3)-ZMPM(KEND)
         SDZ_KBG = Z(1+KSP)-Z(1)
         SDZ_KEND= Z(N3)-Z(N3-KSP)
         AB(KBG) = 1./(VDZ_KBG*SDZ_KBG)
         AF(KEND)= 1./(VDZ_KEND*SDZ_KEND)
        ENDDO

        DO ILEV=NLEV-1,0,-1
          DO J= 1,N2M
          DO K= KKMG(ILEV,1),KKMG(ILEV,2)
           AC(K,J,1)=-1.*(AB(K)+AF(K)+AS(J)+AN(J))
          ENDDO
          ENDDO
   
          DO I= 2,N1MH
          DO J= 1,N2M
          DO K= KKMG(ILEV,1),KKMG(ILEV,2)
            AC(K,J,I)=AC(K,J,1)-AI3(I)
          ENDDO
          ENDDO
          ENDDO
  
        ENDDO

!  CALCULATE INTERPOLATION COEFFS
       DO ILEV= 0,NLEV-1
        KBG=KKMG(ILEV,1)
        KEND=KKMG(ILEV,2)
        VDZ_KEND= ZMPM(KBG)-Z(1)+Z(N3)-ZMPM(KEND)
        COI2(KEND)= (ZMPM(KKMG(ILEV+1,2))-ZMPM(KEND))/VDZ_KEND
        COI1(KEND)= 1.-COI2(KEND)
       ENDDO

       ENDIF


        DO ILEV=NLEV,0,-1
        WRITE(*,*) 'IIMG(1',ILEV,')=',KKMG(ILEV,1)
        WRITE(*,*) 'IIMG(2',ILEV,')=',KKMG(ILEV,2)
        WRITE(*,*) 'IIMG(3',ILEV,')=',KKMG(ILEV,3)
        ENDDO


      RETURN
      END
