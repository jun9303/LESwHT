!=======================================================================
      SUBROUTINE PRINT_REAL_TIME
!=======================================================================
      IMPLICIT NONE
      CHARACTER*8 ::  DATE
      CHARACTER*10::  NOW
      CHARACTER*5 ::  ZONE
      INTEGER*8   ::  VALS(8)

      CALL DATE_AND_TIME(DATE,NOW,ZONE,VALS)

      WRITE(*,101) VALS(1),VALS(2),VALS(3),VALS(5),VALS(6),VALS(7)
      WRITE(*,*) ''
 101  FORMAT(' @ 'I0.4,'-',I0.2,'-',I0.2,' ',I0.2,':',I0.2,':',I0.2)

      RETURN
      END SUBROUTINE PRINT_REAL_TIME
!=======================================================================
      SUBROUTINE READSETTINGS
!=======================================================================
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8    :: N
      CHARACTER*10 :: DUMMY

      OPEN(10,FILE='settings.in')
      READ(10,*) DUMMY
      READ(10,*) RE,PR,GR,GRDIR,T_INF
      READ(10,*) DUMMY
      READ(10,*) IRESET,IREAD,IAVG,IPZERO,EPS_PTR,UBULK_I
      READ(10,*) DUMMY
      READ(10,*) NTST,NPRINT,NPRIAVG,NPIN
      READ(10,*) DUMMY
      READ(10,*) IDTOPT,DT_SIZE,CFLFAC
      READ(10,*) DUMMY
      READ(10,*) RESID1,NLEV,NBLI,IOLDV,MGITR,IMGSOR,WWSOR
      READ(10,*) DUMMY
      READ(10,*) ILES,INSMDL,ITEMDL,IDVMON,CSGSTS,CSGSHF,FILTER
      READ(10,*) DUMMY
      READ(10,*) IBMON,MASSON,IMOVINGON,IHTRANS
      READ(10,*) DUMMY
      READ(10,*) GRIDFILE
      READ(10,*) PREV_FLD
      READ(10,*) DUMMY
      READ(10,*) NTRACE,NTR
      IF (NTRACE .GT. 0) THEN
        DO N = 1, NTRACE
          READ(10,*) TRPTS(N,1),TRPTS(N,2),TRPTS(N,3)
        ENDDO
      ENDIF
      CLOSE(10)

      WRITE(*,*) '========= SETTINGS ========='
      WRITE(*,*) ''
      WRITE(*,101) RE,PR,GR
      WRITE(*,102) IRESET,IREAD,IAVG,IPZERO,EPS_PTR,UBULK_I
      WRITE(*,103) NTST,NPRINT,NPRIAVG,NPIN
      WRITE(*,104) IDTOPT,DT,CFLFAC
      WRITE(*,105) RESID1,NLEV,NBLI,IOLDV,MGITR,IMGSOR,WWSOR
      WRITE(*,106) ILES,INSMDL,ITEMDL,IDVMON,CSGSTS,CSGSHF,FILTER
      WRITE(*,107) IBMON,MASSON,IMOVINGON,IHTRANS
      WRITE(*,108) GRIDFILE
      WRITE(*,109) PREV_FLD
      WRITE(*,*) ''

      IF (NTRACE .GT. 0) THEN
        WRITE(*,*) '========= TRACE POSITION ========='
        WRITE(*,*) ''
        DO  N=1,NTRACE
          WRITE(*,110) N, TRPTS(N,1), TRPTS(N,2), TRPTS(N,3)
        ENDDO
      ENDIF

 101  FORMAT('  RE=',ES11.3,'  PR=',ES11.3,'  GR=',ES11.3)
 102  FORMAT('  IRESET=',I5,'  IREAD=',I5,'  IAVG=',I5,'  IPZERO=',I5,'  EPS_PTR=',F7.3,'  UBULK_I=',F7.3)
 103  FORMAT('  NTST=',I10,'  NPRINT=',I8,'  NPRIAVG=',I8,'  NPIN=',I5)
 104  FORMAT('  IDTOPT=',I5,'  DT=',ES13.5,'  CFLFAC=',F11.3)
 105  FORMAT('  RESID=',ES12.4,'  NLEV=',I5,'  NBLI=',I5,'  IOLDV=',I5,'  MGITR=',I5,'  IMGSOR=',I5,'  WWSOR=',F7.3)
 106  FORMAT('  ILES=',I2,'  INSMDL=',I2,'  ITEMDL=',I2,'  IDVMON=',I2,'  CSGSTS='F12.4,'  CSGSHF=',F12.4,'  IFILTER=',I5)
 107  FORMAT('  IBMON=',I2,'  MASSON=',I2,'  IMOVINGON=',I2,'  IHTRANS=',I2)
 108  FORMAT('  GRIDFILE=',A25)
 109  FORMAT('  PREV_FLD=',A25)
 110  FORMAT(I5,3I6)

      RETURN
      END SUBROUTINE READSETTINGS
!=======================================================================
!=======================================================================
      SUBROUTINE READBCS
!=======================================================================
      USE MOD_COMMON
      IMPLICIT NONE
      CHARACTER*10 :: DUMMY

      OPEN(10,FILE='boundary.in')
      READ(10,*) DUMMY
      READ(10,*) BC_YBTM,BC_YTOP
      READ(10,*) DUMMY
      READ(10,*) BC_ZBTM,BC_ZTOP
      READ(10,*) DUMMY
      READ(10,*) ICH,ICONJG
      CLOSE(10)

      OPEN(11,FILE='../output/ibmpre/ibmpre_prdic.bin')
      READ(11,*) XPRDIC,YPRDIC,ZPRDIC,IINTP
      CLOSE(11)

      IF (XPRDIC .EQ. 1) THEN
        BC_XBTM = 4 ! 4:   PERIODIC B.C.
        BC_XTOP = 4 
      ELSE
        BC_XBTM = 2 ! 2:   VELOCITY INLET B.C.
        BC_XTOP = 3 ! 3:   CONVECTIVE OUTLET B.C.
      ENDIF

      IF (YPRDIC .EQ. 1) THEN
        BC_YBTM = 4
        BC_YTOP = 4
      ENDIF

      IF (ZPRDIC .EQ. 1) THEN
        BC_ZBTM = 4
        BC_ZTOP = 4
      ENDIF

      WRITE(*,*) '========= BOUNDARY CONDITIONS ========='
      WRITE(*,*) ''
      WRITE(*,*) '(0: WALL; 1: FARF; 2: VELIN; 3: CONVOUT; 4: PRDIC)'
      WRITE(*,101) BC_XBTM,BC_XTOP
      WRITE(*,102) BC_YBTM,BC_YTOP
      WRITE(*,103) BC_ZBTM,BC_ZTOP
      WRITE(*,*) ''
      WRITE(*,104) ICH, ICONJG, IINTP
      WRITE(*,*) ''

 101  FORMAT('  BC_XBTM=',I2,'  BC_XTOP=',I2)
 102  FORMAT('  BC_YBTM=',I2,'  BC_YTOP=',I2)
 103  FORMAT('  BC_ZBTM=',I2,'  BC_ZTOP=',I2)
 104  FORMAT('  ICH=',I2,'  ICONJG=',I2,'  IINTP=',I2)

      RETURN
      END SUBROUTINE READBCS
!=======================================================================
!=======================================================================
      SUBROUTINE READGEOM
!=======================================================================
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8    :: I,J,K

      OPEN(10,FILE=GRIDFILE)
      READ(10,*) N1,N2,N3
      READ(10,*) N1M,N2M,N3M
      READ(10,*) XL,YL,ZL

      CALL ALLO(N1,N2,N3)

      READ(10,*) (X(I),I=1,N1)
      READ(10,*) (Y(J),J=1,N2)
      READ(10,*) (Z(K),K=1,N3)
      READ(10,*) (IPV(I),I=1,N1M)
      READ(10,*) (JPV(J),J=1,N2M)
      READ(10,*) (KPV(K),K=1,N3M)
      READ(10,*) (IMV(I),I=1,N1M)
      READ(10,*) (JMV(J),J=1,N2M)
      READ(10,*) (KMV(K),K=1,N3M)
      READ(10,*) (FIXIL(I),I=1,N1M)
      READ(10,*) (FIXJL(J),J=1,N2M)
      READ(10,*) (FIXKL(K),K=1,N3M)
      READ(10,*) (FIXIU(I),I=1,N1M)
      READ(10,*) (FIXJU(J),J=1,N2M)
      READ(10,*) (FIXKU(K),K=1,N3M)
      READ(10,*) (C2CX(I),I=0,N1)
      READ(10,*) (C2CY(J),J=0,N2)
      READ(10,*) (C2CZ(K),K=0,N3)
      READ(10,*) (C2CXI(I),I=0,N1)
      READ(10,*) (C2CYI(J),J=0,N2)
      READ(10,*) (C2CZI(K),K=0,N3)
      READ(10,*) (F2FX(I),I=0,N1)
      READ(10,*) (F2FY(J),J=0,N2)
      READ(10,*) (F2FZ(K),K=0,N3)
      READ(10,*) (F2FXI(I),I=0,N1)
      READ(10,*) (F2FYI(J),J=0,N2)
      READ(10,*) (F2FZI(K),K=0,N3)
      READ(10,*) (XMP(I),I=0,N1)
      READ(10,*) (YMP(J),J=0,N2)
      READ(10,*) (ZMP(K),K=0,N3)

      CLOSE(10)

      RETURN
      END SUBROUTINE READGEOM
!=======================================================================
!=======================================================================
      SUBROUTINE ALLOINIT
!=======================================================================
      USE MOD_COMMON
      USE MOD_FLOWARRAY
      IMPLICIT NONE
      INTEGER*8    :: N

      OPEN(11,FILE='../output/ibmpre/ibmpre_fcpts.bin')
      READ(11,*) NINTP(1), NINTP(2), NINTP(3)
      READ(11,*) NBODY(1), NBODY(2), NBODY(3)
      DO N = 1,3
        NBODY(N) = NBODY(N) + NINTP(N)
      ENDDO
      CLOSE(11)

      OPEN(11,FILE='../output/ibmpre/ibmpre_nutzero.bin')
      READ(11,*) NZERO
      CLOSE(11)

      IF (IHTRANS .EQ. 1) THEN
        OPEN(11,FILE='../output/ibmpre/ibmpre_fcpts_t.bin')
        READ(11,*) NINTP(4)
        READ(11,*) NBODY(4)
        NBODY(4) = NBODY(4) + NINTP(4)
        CLOSE(11)
      ENDIF

      RETURN
      END SUBROUTINE ALLOINIT
!=======================================================================
!=======================================================================
      SUBROUTINE ALLO_ARRAY
!=======================================================================
      USE MOD_COMMON
      USE MOD_FLOWARRAY
      IMPLICIT NONE

      IF (IHTRANS .EQ. 0) THEN
        CALL BASIC_ALLO()
      ELSE
        CALL THERMAL_ALLO()
      ENDIF

      IF (ILES .EQ. 1) THEN
        CALL LES_ALLO()
        IF (IHTRANS .EQ. 1) CALL LES_THERMAL_ALLO()
      ENDIF

      IF (ICONJG .EQ. 1) THEN
        CALL CONJG_ALLO()
      ENDIF

      IF (IAVG .EQ. 1) CALL AVG_ALLO()

      RETURN
      END SUBROUTINE ALLO_ARRAY
!=======================================================================
!=======================================================================
      SUBROUTINE IBMPREREAD
!=======================================================================
      USE MOD_COMMON
      USE MOD_FLOWARRAY
      IMPLICIT NONE
      INTEGER*8    :: N,L,DUMMY
      INTEGER*8    :: I,J,K,ITMP

      OPEN(11,FILE='../output/ibmpre/ibmpre_fcpts.bin')
      READ(11,*) DUMMY ! NINTP, ALREADY READ IN ALLOINIT SUBROUTINE
      READ(11,*) DUMMY ! NBODY, ALREADY READ IN ALLOINIT SUBROUTINE
      DO N = 1,NBODY(1)
        READ(11,*) IFC(N,1), JFC(N,1), KFC(N,1)
      ENDDO
      DO N = 1,NBODY(2)
        READ(11,*) IFC(N,2), JFC(N,2), KFC(N,2)
      ENDDO
      DO N = 1,NBODY(3)
        READ(11,*) IFC(N,3), JFC(N,3), KFC(N,3)
      ENDDO
      DO N = 1,NINTP(1)
        READ(11,*) INTPINDX(N,1,1), INTPINDX(N,1,2), INTPINDX(N,1,3)
      ENDDO
      DO N = 1,NINTP(2)
        READ(11,*) INTPINDX(N,2,1), INTPINDX(N,2,2), INTPINDX(N,2,3)
      ENDDO
      DO N = 1,NINTP(3)
        READ(11,*) INTPINDX(N,3,1), INTPINDX(N,3,2), INTPINDX(N,3,3)
      ENDDO
      DO N = 1,NINTP(1)
        READ(11,*) (((GEOMFAC(N,1,I,J,K),K=0,2),J=0,2),I=0,2)
      ENDDO
      DO N = 1,NINTP(2)
        READ(11,*) (((GEOMFAC(N,2,I,J,K),K=0,2),J=0,2),I=0,2)
      ENDDO
      DO N = 1,NINTP(3)
        READ(11,*) (((GEOMFAC(N,3,I,J,K),K=0,2),J=0,2),I=0,2)
      ENDDO
      CLOSE(11)

      IF (IHTRANS .EQ. 1) THEN
        OPEN(11,FILE='../output/ibmpre/ibmpre_fcpts_t.bin')
        READ(11,*) DUMMY ! NINTP, ALREADY READ IN ALLOINIT SUBROUTINE
        READ(11,*) DUMMY ! NBODY, ALREADY READ IN ALLOINIT SUBROUTINE
        DO N = 1,NBODY(4)
          READ(11,*) IFC(N,4), JFC(N,4), KFC(N,4)
        ENDDO
        DO N = 1,NINTP(4)
          READ(11,*) INTPINDX(N,4,1), INTPINDX(N,4,2), INTPINDX(N,4,3)
        ENDDO
        DO N = 1,NINTP(4)
          READ(11,*) (((GEOMFAC(N,4,I,J,K),K=0,2),J=0,2),I=0,2)
        ENDDO
      ENDIF

      WRITE(*,*) '========= LOADING IBM-PREPROCESSING DATA ========='
      WRITE(*,*) ''
      WRITE(*,35) NINTP(1),NINTP(2),NINTP(3)
      WRITE(*,36) NBODY(1),NBODY(2),NBODY(3)
      IF (IHTRANS .EQ. 1) WRITE(*,37) NINTP(4),NBODY(4)
      WRITE(*,*) ''
 35   FORMAT('  NINTP_U :',I12,'  NINTP_V :',I12,'  NINTP_W :',I12)
 36   FORMAT('  NBODY_U :',I12,'  NBODY_V :',I12,'  NBODY_W :',I12)
 37   FORMAT('  NINTP_T :',I12,'  NBODY_T :',I12)

      RETURN
      END SUBROUTINE IBMPREREAD
!=======================================================================
!=======================================================================
      SUBROUTINE NUTZEROREAD
!=======================================================================
      USE MOD_COMMON
      USE MOD_FLOWARRAY
      IMPLICIT NONE
      INTEGER*8     :: I,J,K
      INTEGER*8 :: N

      OPEN(14,FILE='../output/ibmpre/ibmpre_nutzero.bin')
      READ(14,*) NZERO
      READ(14,*) (INZ(N),N=1,NZERO)
      READ(14,*) (JNZ(N),N=1,NZERO)
      READ(14,*) (KNZ(N),N=1,NZERO)
      CLOSE(14)

      OPEN(15,FILE='../output/ibmpre/ibmpre_wallfdvm.bin')
      READ(15,*) ((( NWALL_DVM(I,J,K) ,I=1,N1M),J=1,N2M),K=1,N3M)
      CLOSE(15)

      WRITE(*,*) '========= LOADING LES-PREPROCESSING DATA ========='
      WRITE(*,*) ''
      WRITE(*,30) NZERO
      WRITE(*,*) ''
 30   FORMAT('  # OF NUTZERO PTS = ',I10)

      RETURN
      END SUBROUTINE NUTZEROREAD
!=======================================================================
!=======================================================================
      SUBROUTINE CONJGREAD
!=======================================================================
      USE MOD_COMMON
      USE MOD_FLOWARRAY
      IMPLICIT NONE
      INTEGER*8     :: I,J,K
      INTEGER*8     :: L

      OPEN(15,FILE='../output/ibmpre/ibmpre_conjg.bin')
      READ(15,*) ((( CSTAR(I,J,K) ,I=1,N1M),J=1,N2M),K=1,N3M)
      READ(15,*) (((( KSTAR(I,J,K,L) ,I=1,N1M),J=1,N2M),K=1,N3M),L=1,6)
      CLOSE(15)

      WRITE(*,*) '========= LOADING CONJUGATE H.TRANS DATA ========='
      WRITE(*,*) ''

      RETURN
      END SUBROUTINE CONJGREAD
!=======================================================================
!=======================================================================
       SUBROUTINE MAKEFLD
!=======================================================================
!$    use omp_lib
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,P,T,QMASS
      IMPLICIT NONE
      INTEGER*8     :: I,J,K
      REAL*8        :: FUNCBODY
      REAL*8        :: FLOWAREA(N1), FL, FL_S, ADJ

      IHIST = 0
      TIME = 0.

      U = 0.
      V = 0.
      W = 0.
      P = 0.

      FLOWAREA = 0.

      WRITE(*,*) '========= MAKING INITIAL FIELD ========='
      WRITE(*,*) ''

!$OMP PARALLEL DO
      DO I = 1,N1
        DO J = 0,N2
          DO K = 0,N3
             U(I,J,K) = UBULK_I
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      IF (BC_YBTM .EQ. 0) THEN
        DO I = 1,N1
          DO K = 0,N3
            U(I,0,K) = 0.
            U(I,1,K) = 0.
          ENDDO
        ENDDO
      ENDIF

      IF (BC_YTOP .EQ. 0) THEN
        DO I = 1,N1
          DO K = 0,N3
            U(I,N2M,K) = 0.
            U(I,N2,K) = 0.
          ENDDO
        ENDDO
      ENDIF

      IF (BC_ZBTM .EQ. 0) THEN
        DO I = 1,N1
          DO J = 0,N2
            U(I,J,0) = 0.
            U(I,J,1) = 0.
          ENDDO
        ENDDO
      ENDIF

      IF (BC_ZTOP .EQ. 0) THEN
        DO I = 1,N1
          DO J = 0,N2
            U(I,J,N3M) = 0.
            U(I,J,N3) = 0.
          ENDDO
        ENDDO
      ENDIF   

      IF (IHTRANS .EQ. 1) THEN
        IF (T_INF .EQ. 0) THEN
          T = -1.
        ELSE
          T = 1.
        ENDIF

!$OMP PARALLEL DO
        DO I = 0,N1
          DO J = 1,N2M
            DO K = 0,N3
              IF(FUNCBODY(XMP(I),YMP(J),ZMP(K),TIME).GE.1.E-10) THEN
                T(I,J,K) = 0.
              ENDIF
            ENDDO
          ENDDO
        ENDDO
!$OMP END PARALLEL DO

      ENDIF

      RETURN
      END SUBROUTINE MAKEFLD
!=======================================================================
!=======================================================================
       SUBROUTINE PREFLD
!=======================================================================
!$    use omp_lib
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,T,P
      IMPLICIT NONE
      INTEGER*8     :: I,J,K
      INTEGER*8     :: NN1,NN2,NN3
      REAL*8        :: RRE,PPR,GGR
      
      WRITE(*,*) '========= LOADING PREVIOUS FIELD ========='
      WRITE(*,*) ''

      OPEN(12,FILE=PREV_FLD,FORM='UNFORMATTED')
      READ(12) NN1, NN2, NN3, RRE, PPR, GGR

      IF ((NN1 .NE. N1) .OR. (NN2 .NE. N2) .OR. (NN3 .NE. N3)) THEN
        WRITE(*,*) ' --- PREVIOUS FIELD DOES NOT MATCH WITH THE GRID SYSTEM.'
        CLOSE(12)
        STOP

      ELSE

        READ(12) IHIST, M, TIME, DT
        READ(12) XPRDIC, YPRDIC, ZPRDIC
        READ(12) BC_XBTM, BC_XTOP, BC_YBTM, BC_YTOP, BC_ZBTM, BC_ZTOP
        READ(12) ICH, ICONJG
        READ(12) ((( U(I,J,K) ,I=1,N1),J=0,N2),K=0,N3)
        READ(12) ((( V(I,J,K) ,I=0,N1),J=1,N2),K=0,N3)
        READ(12) ((( W(I,J,K) ,I=0,N1),J=0,N2),K=1,N3)
        READ(12) ((( P(I,J,K) ,I=1,N1M),J=1,N2M),K=1,N3M)
        IF (IHTRANS .EQ. 1) THEN
          IF (T_INF .EQ. 0) THEN
            T = -1.
          ELSE
            T = 1.
          ENDIF
          READ(12) ((( T(I,J,K), I=1,N1M),J=1,N2M),K=1,N3M)
        ENDIF
        CLOSE(12)

        IF (IPZERO .EQ. 1) P = 0.

        WRITE(*,*) ''
        WRITE(*,*) ' ----------- PREVIOUS FIELD INFORMATION -----------'
        WRITE(*,100) PREV_FLD
        WRITE(*,101) RRE, PPR, GGR
        WRITE(*,102) NN1, NN2, NN3
        WRITE(*,103) IHIST, M, TIME, DT
        WRITE(*,*) ''
 100    FORMAT(' PREVIOUS FIELD LOCATION : ', A25)
 101    FORMAT(' RE = ',ES12.3,' PR = ',ES12.3,' GR = ',ES12.3)
 102    FORMAT(' N1 = ',I12,' N2 = ',I12,' N3 = ',I12)
 103    FORMAT(' IHIST = ',I9,' M = ',I9,' TIME = ',F10.5,' DT = ',F12.8)

        IF (XPRDIC .EQ. 1) THEN
!$OMP PARALLEL DO
          DO J=0,N2
          DO K=0,N3
             U(0 ,J,K)=U(N1M,J,K)
             U(N1,J,K)=U(1  ,J,K)
          ENDDO
          ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
          DO J=1,N2
          DO K=0,N3
             V(0 ,J,K)=V(N1M,J,K)
             V(N1,J,K)=V(1  ,J,K)
          ENDDO
          ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
          DO J=0,N2
          DO K=1,N3
             W(0 ,J,K)=W(N1M,J,K)
             W(N1,J,K)=W(1  ,J,K)
          ENDDO
          ENDDO
!$OMP END PARALLEL DO

          IF (IHTRANS .EQ. 1) THEN
!$OMP PARALLEL DO
           DO J=0,N2
            DO K=0,N3
               T(0 ,J,K)=T(N1M,J,K)
               T(N1,J,K)=T(1  ,J,K)
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
          ENDIF
        ENDIF            

        IF (YPRDIC .EQ. 1) THEN
!$OMP PARALLEL DO
          DO K=0,N3
          DO I=1,N1
             U(I,0,K) =U(I,N2M,K)
             U(I,N2,K)=U(I,1,K)
          ENDDO
          ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
          DO K=0,N3
          DO I=0,N1 
             V(I,0,K) =V(I,N2M,K)
             V(I,N2,K)=V(I,1,K)
          ENDDO
          ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
          DO K=1,N3
          DO I=0,N1
             W(I,0,K) =W(I,N2M,K)
             W(I,N2,K)=W(I,1,K)
          ENDDO
          ENDDO
!$OMP END PARALLEL DO

          IF (IHTRANS .EQ. 1) THEN
!$OMP PARALLEL DO
           DO K=0,N3
            DO I=0,N1
               T(I,0 ,K)=T(I,N2M,K)
               T(I,N2,K)=T(I,1  ,K)
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
          ENDIF

        ENDIF

        IF (ZPRDIC .EQ. 1) THEN
!$OMP PARALLEL DO
          DO I=1,N1
          DO J=0,N2
             U(I,J,0) =U(I,J,N3M)
             U(I,J,N3)=U(I,J,1)
          ENDDO
          ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
          DO I=0,N1
          DO J=1,N2
             V(I,J,0) =V(I,J,N3M)
             V(I,J,N3)=V(I,J,1)
          ENDDO
          ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
          DO I=0,N1
          DO J=0,N2
             W(I,J,0) =W(I,J,N3M)
             W(I,J,N3)=W(I,J,1)
          ENDDO
          ENDDO
!$OMP END PARALLEL DO

          IF (IHTRANS .EQ. 1) THEN
!$OMP PARALLEL DO
           DO I=0,N1
            DO J=0,N2
               T(I,J,0)=W(I,J,N3M)
               T(I,J,N3)=T(I,J,1)
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
          ENDIF

        ENDIF

      ENDIF

      RETURN
      END SUBROUTINE PREFLD
!=======================================================================
!=======================================================================
       SUBROUTINE ADDPERTURB
!=======================================================================
!$    use omp_lib
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,P,QMASS
      IMPLICIT NONE
      REAL*8        :: PTB_ARRAY(0:N1,0:N2,0:N3)
      INTEGER*8     :: I,J,K
      REAL*8        :: FUNCBODY
      REAL*8        :: FLOWAREA, PERTB_RATE, ADJ, FLOWRATE

      CALL RANDOM_NUMBER(PTB_ARRAY)
      PTB_ARRAY = (PTB_ARRAY - 1.) * 2. * EPS_PTR
      ! NOW PTB_ARRAY ~ UNIFORMDIST.(-EPS_PTR, +EPS_PTR)

      DO I = 2,N1M
        FLOWAREA = 0
        PERTB_RATE = 0.
!$OMP PARALLEL DO reduction(+:FLOWAREA, PERTB_RATE)
        DO J = 1,N2M
          DO K = 1,N3M
              IF (FUNCBODY(X(I),YMP(J),ZMP(K),TIME).GT.1.E-10) THEN
                U(I,J,K) = U(I,J,K) + PTB_ARRAY(I,J,K)
                FLOWAREA = FLOWAREA + F2FY(J) * F2FZ(K)
                PERTB_RATE = PERTB_RATE + PTB_ARRAY(I,J,K) * F2FY(J) * F2FZ(K)
              ENDIF
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
          ADJ = PERTB_RATE / FLOWAREA
          FLOWRATE = 0
!$OMP PARALLEL DO reduction(+:FLOWRATE)
        DO J = 1,N2M
          DO K = 1,N3M
              IF (FUNCBODY(X(I),YMP(J),ZMP(K),TIME).GT.1.E-10) THEN
                U(I,J,K) = U(I,J,K) - ADJ
                FLOWRATE = FLOWRATE + U(I,J,K) * F2FY(J) * F2FZ(K)
              ENDIF
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ! WRITE(*,*)FLOWRATE
      ENDDO

      CALL RANDOM_NUMBER(PTB_ARRAY)
      PTB_ARRAY = (PTB_ARRAY - 1.) * 2. * EPS_PTR
      ! NOW PTB_ARRAY ~ UNIFORMDIST.(-EPS_PTR, +EPS_PTR)

      DO J = 2,N2M
        FLOWAREA = 0
        PERTB_RATE = 0.
!$OMP PARALLEL DO reduction(+:FLOWAREA, PERTB_RATE)
        DO K = 1,N3M
          DO I = 1,N1M
              IF (FUNCBODY(XMP(I),Y(J),ZMP(K),TIME).GT.1.E-10) THEN
                V(I,J,K) = V(I,J,K) + PTB_ARRAY(I,J,K)
                FLOWAREA = FLOWAREA + F2FZ(K) * F2FX(I)
                PERTB_RATE = PERTB_RATE + PTB_ARRAY(I,J,K) * F2FZ(K) * F2FX(I)
              ENDIF
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
        ADJ = PERTB_RATE / FLOWAREA
        FLOWRATE = 0
!$OMP PARALLEL DO reduction(+:FLOWRATE)
        DO K = 1,N3M
          DO I = 1,N1M
              IF (FUNCBODY(XMP(I),Y(J),ZMP(K),TIME).GT.1.E-10) THEN
                V(I,J,K) = V(I,J,K) - ADJ
                FLOWRATE = FLOWRATE + V(I,J,K) * F2FZ(K) * F2FX(I)
              ENDIF
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ! WRITE(*,*)FLOWRATE
      ENDDO

      CALL RANDOM_NUMBER(PTB_ARRAY)
      PTB_ARRAY = (PTB_ARRAY - 1.) * 2. * EPS_PTR
      ! NOW PTB_ARRAY ~ UNIFORMDIST.(-EPS_PTR, +EPS_PTR)

      DO K = 3,N3M
        FLOWAREA = 0
        PERTB_RATE = 0.
!$OMP PARALLEL DO reduction(+:FLOWAREA, PERTB_RATE)
        DO I = 1,N1M
          DO J = 1,N2M
              IF (FUNCBODY(XMP(I),YMP(J),Z(K),TIME).GT.1.E-10) THEN
                W(I,J,K) = W(I,J,K) + PTB_ARRAY(I,J,K)
                FLOWAREA = FLOWAREA + F2FX(I) * F2FY(J)
                PERTB_RATE = PERTB_RATE + PTB_ARRAY(I,J,K) * F2FX(I) * F2FY(J)
              ENDIF
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
        ADJ = PERTB_RATE / FLOWAREA
        FLOWRATE = 0
!$OMP PARALLEL DO reduction(+:FLOWRATE)
        DO I = 1,N1M
          DO J = 1,N2M
              IF (FUNCBODY(XMP(I),YMP(J),Z(K),TIME).GT.1.E-10) THEN
                W(I,J,K) = W(I,J,K) - ADJ
                FLOWRATE = FLOWRATE + W(I,J,K) * F2FX(I) * F2FY(J)
              ENDIF
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ! WRITE(*,*)FLOWRATE
      ENDDO

      WRITE(*,*) '--- ADDING ISOTROPIC PERTURBATION FROM THE UNIFORM DIST.'
      WRITE(*,201) EPS_PTR*100.
      WRITE(*,*) ''
 201  FORMAT(' --- MAX. SIZE OF PERTURBATION = ', F5.1, '% OF U_BULK')

      RETURN
      END SUBROUTINE ADDPERTURB
!=======================================================================
!=======================================================================
      SUBROUTINE DTTIMEINIT
!=======================================================================
      USE MOD_COMMON
      IMPLICIT NONE

      IF (IDTOPT .EQ. 0) DT = DT_SIZE ! FOR CONSTANT USAGE PURPOSE
      IF ((IDTOPT .NE.0) .AND. (IREAD .EQ. 0)) DT = DT_SIZE 
                                      ! FOR INITAL USAGE PURPOSE
      IF (IRESET .EQ. 1) THEN
        IHIST = 0
        TIME  = 0.
      ENDIF
      NTIME = 0

      IF (IAVG .EQ. 1) THEN
        TIMEINIT = TIME
        IHISTINIT = IHIST
      ENDIF

      RETURN
      END SUBROUTINE DTTIMEINIT
!=======================================================================
!=======================================================================
       SUBROUTINE RK3COEFINIT
!=======================================================================
      USE MOD_COMMON
      IMPLICIT NONE

      GAMMA(1) = 8./15.
      GAMMA(2) = 5./12.
      GAMMA(3) = 3./4.
      RO(1) = 0.
      RO(2) = -17./60.
      RO(3) = -5./12.

      RETURN
      END SUBROUTINE RK3COEFINIT
!=======================================================================
!=======================================================================
      SUBROUTINE CONVERGENCE_CHECK
!=======================================================================
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,QMASS
      IMPLICIT NONE
      INTEGER*8     :: I,J,K
      REAL*8 :: DVG11,DVG12,DVG13,DVG1

      DVMAX = 0.

!$OMP PARALLEL DO                        &
!$OMP private(DVG11,DVG12,DVG13,DVG1)    &
!$OMP reduction(MAX:DVMAX)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
         DVG11=(U(I+1,J,K)-U(I,J,K))*F2FXI(I)
         DVG12=(V(I,J+1,K)-V(I,J,K))*F2FYI(J)
         DVG13=(W(I,J,K+1)-W(I,J,K))*F2FZI(K)
         IF (MASSON.EQ.1) THEN
          DVG1=DVG11+DVG12+DVG13-QMASS(I,J,K)
         ELSE
          DVG1=DVG11+DVG12+DVG13
         ENDIF
         DVMAX=AMAX1(ABS(DVG1),DVMAX)
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CONVERGENCE_CHECK
!=======================================================================
!=======================================================================
      SUBROUTINE MASSCHECK
!=======================================================================
!
!     Calculate maximum mass source/sink in IBM
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : QMASS
      IMPLICIT NONE
      INTEGER*8     :: I,J,K

      QMMAX = 0.

!$OMP PARALLEL DO &
!$OMP reduction(MAX:QMMAX)
      DO 100 K=1,N3M
      DO 100 J=1,N2M
      DO 100 I=1,N1M
      QMMAX = AMAX1(QMMAX,ABS(QMASS(I,J,K)))
 100  CONTINUE

      RETURN
      END SUBROUTINE MASSCHECK
!=======================================================================
      SUBROUTINE DEALLO_ARRAY
!=======================================================================
      USE MOD_COMMON
      USE MOD_FLOWARRAY
      IMPLICIT NONE

      IF (IHTRANS .EQ. 0) THEN
        CALL BASIC_DEALLO()
      ELSE
        CALL THERMAL_DEALLO()
      ENDIF

      IF (ILES .EQ. 1) THEN
        CALL LES_DEALLO()
        IF (IHTRANS .EQ. 1) CALL LES_THERMAL_DEALLO()
      ENDIF

      IF (IAVG .EQ. 1) CALL AVG_DEALLO()

      RETURN
      END SUBROUTINE DEALLO_ARRAY
!=======================================================================
