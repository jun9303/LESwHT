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

      OPEN(77,FILE='mgftresiduemax.dat')

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
