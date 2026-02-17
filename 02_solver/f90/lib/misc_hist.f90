!=======================================================================
       SUBROUTINE FTRFILES(IO)
!=======================================================================
!$    use omp_lib
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8    :: IO

      IF (IO.NE.0) THEN
        IF (IREAD.NE.1) THEN
          OPEN(2000,FILE='../output/ftr/fhist.dat')
          OPEN(2001,FILE='../output/ftr/fcdcl.dat')
          OPEN(2002,FILE='../output/ftr/fnusselt.dat')          
          IF (ILES.EQ.1) THEN
            OPEN(2003,FILE='../output/ftr/fles.dat')
            IF(IHTRANS.EQ.1) THEN
              OPEN(2004,FILE='../output/ftr/flest.dat')
            ENDIF
          ENDIF
          IF (NTRACE.NE.0) THEN
            OPEN(2005,FILE='../output/ftr/futrace.dat')
            OPEN(2006,FILE='../output/ftr/fvtrace.dat')
            OPEN(2007,FILE='../output/ftr/fwtrace.dat')
            OPEN(2008,FILE='../output/ftr/fptrace.dat')
            IF (IHTRANS.EQ.1) THEN
              OPEN(2009,FILE='../output/ftr/fttrace.dat')
            ENDIF
          ENDIF
          IF (ICH.EQ.1) THEN
            OPEN(2010,FILE='../output/ftr/fcmfr.dat')
          ENDIF
          OPEN(2011,FILE='../output/ftr/ftime.dat')
        ELSE
          OPEN(2000,FILE='../output/ftr/fhist.dat',&
            POSITION='APPEND')
          OPEN(2001,FILE='../output/ftr/fcdcl.dat',&
            POSITION='APPEND')
          OPEN(2002,FILE='../output/ftr/fnusselt.dat',&
            POSITION='APPEND')          
          IF (ILES.EQ.1) THEN
            OPEN(2003,FILE='../output/ftr/fles.dat',&
                POSITION='APPEND')
            IF (IHTRANS.EQ.1) THEN
              OPEN(2004,FILE='../output/ftr/flest.dat',&
                POSITION='APPEND')
            ENDIF
          ENDIF
          IF (NTRACE.NE.0) THEN
            OPEN(2005,FILE='../output/ftr/futrace.dat',&
                POSITION='APPEND')
            OPEN(2006,FILE='../output/ftr/fvtrace.dat',&
                POSITION='APPEND')
            OPEN(2007,FILE='../output/ftr/fwtrace.dat',&
                POSITION='APPEND')
            OPEN(2008,FILE='../output/ftr/fptrace.dat',&
                POSITION='APPEND')
            IF (IHTRANS.EQ.1) THEN
              OPEN(2009,FILE='../output/ftr/fttrace.dat',&
                   POSITION='APPEND')
            ENDIF
          ENDIF
          IF (ICH.EQ.1) THEN
            OPEN(2010,FILE='../output/ftr/fcmfr.dat',&
                POSITION='APPEND')
          ENDIF
          OPEN(2011,FILE='../output/ftr/ftime.dat',&
            POSITION='APPEND')
        ENDIF
      ELSE
        CLOSE(2000)
        CLOSE(2001)
        CLOSE(2002)
        IF (ILES.EQ.1) THEN
          CLOSE(2003)
          IF (IHTRANS.EQ.1) THEN
            CLOSE(2004)
          ENDIF
        ENDIF
        IF (NTRACE.NE.0) THEN
          CLOSE(2005)
          CLOSE(2006)
          CLOSE(2007)
          CLOSE(2008)
          IF (IHTRANS.EQ.1) THEN
            CLOSE(2009)
          ENDIF
        ENDIF
        IF (ICH.EQ.1) THEN
          CLOSE(2010)
        ENDIF
        CLOSE(2011)

      ENDIF

      RETURN
      END SUBROUTINE FTRFILES
!=======================================================================
!=======================================================================
      SUBROUTINE WRITEHISTORY
!=======================================================================
!
!     Trace the maximum CFL#, divergenceU, and QMASS(IBM).
!     The tracing interval is set in lica.in (NPIN)
!
!     *In case of channel flow,
!     - the flow rate and the mean pressure gradient is traced
!     - QFLUX should be constant, while PMI, PMIL, PMIU fluctuate
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY
      IMPLICIT NONE
      INTEGER*8     :: I,J,K

      QFLUX = 0.

!$OMP PARALLEL DO &
!$OMP reduction(+:QFLUX)
      DO K=1,N3M
      DO J=1,N2M
      QFLUX = QFLUX+U(1,J,K)*F2FY(J)*F2FZ(K)
      ENDDO
      ENDDO

      WRITE(2000,130) TIME,DT,CFLMAX,DVMAX,QMMAX
 130  FORMAT(F13.5,4ES15.7)

      IF (ICH.EQ.1) THEN
      WRITE(2010,140) TIME,QFLUX,-PMI
 140  FORMAT(F13.5,6ES20.12)
      ENDIF

      RETURN
      END SUBROUTINE WRITEHISTORY
!=======================================================================
!=======================================================================
      SUBROUTINE TRACER
!=======================================================================
!
!     This subroutine is called when MOD(M,NTPRINT)=0.
!     Velocity components and the pressure at the cell center of the
!     tracing cells are saved.
!
!     NTPRINT : determining the time-tracing interval
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,P
      IMPLICIT NONE
      INTEGER*8 :: N

      WRITE(2005,101)TIME,(0.5*(U(TRPTS(N,1),TRPTS(N,2),TRPTS(N,3))              &
                            +U(TRPTS(N,1)+1,TRPTS(N,2),TRPTS(N,3))),N=1,NTRACE)
      WRITE(2006,101)TIME,(0.5*(V(TRPTS(N,1),TRPTS(N,2),TRPTS(N,3))              &
                            +V(TRPTS(N,1),TRPTS(N,2)+1,TRPTS(N,3))),N=1,NTRACE)
      WRITE(2007,101)TIME,(0.5*(W(TRPTS(N,1),TRPTS(N,2),TRPTS(N,3))              &
                            +W(TRPTS(N,1),TRPTS(N,2),TRPTS(N,3)+1)),N=1,NTRACE)
      WRITE(2008,101)TIME,(P(TRPTS(N,1),TRPTS(N,2),TRPTS(N,3)),N=1,NTRACE)
  101 FORMAT(F15.7,10000ES15.7)

      RETURN
      END SUBROUTINE TRACER
!=======================================================================
!=======================================================================
      SUBROUTINE WRITEFIELD
!=======================================================================
!
!     Make an output file of instantaneous field, when MOD(NTIME,NPRINT)=0,
!     NPRINT : instantaneous field file printing interval
!     tfn1='fld': prefix for instantaneous flow field
!     tname     : instantaneous field file name ex) fld006100
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,T,P
      IMPLICIT NONE
      INTEGER*8     :: I,J,K
      INTEGER*8     :: IDUM,IDG1,IDG2,IDG3,IDG4,IDG5,IDG6
      REAL*8        :: DUM
      CHARACTER*25  :: TNAME
      CHARACTER*19  :: TFN1

      IDUM = 0
      DUM  = 0.

      TFN1='../output/field/fld'
      IDG1=IHIST/100000
      IDG2=(IHIST-IDG1*100000)/10000
      IDG3=(IHIST-IDG1*100000-IDG2*10000)/1000
      IDG4=(IHIST-IDG1*100000-IDG2*10000-IDG3*1000)/100
      IDG5=(IHIST-IDG1*100000-IDG2*10000-IDG3*1000-IDG4*100)/10
      IDG6=IHIST-IDG1*100000-IDG2*10000-IDG3*1000-IDG4*100-IDG5*10
      TNAME=TFN1//CHAR(IDG1+48)//CHAR(IDG2+48)//                       &
           CHAR(IDG3+48)//CHAR(IDG4+48)//CHAR(IDG5+48)//CHAR(IDG6+48)

      OPEN(NV,FILE=TNAME,FORM='UNFORMATTED')
      WRITE(NV) N1,N2,N3,RE,PR,GR
      WRITE(NV) IHIST,M,TIME,DT
      WRITE(NV) XPRDIC, YPRDIC, ZPRDIC
      WRITE(NV) BC_XBTM, BC_XTOP, BC_YBTM, BC_YTOP, BC_ZBTM, BC_ZTOP
      WRITE(NV) ICH, ICONJG
      WRITE(NV) ((( U(I,J,K) ,I=1,N1),J=0,N2),K=0,N3)
      WRITE(NV) ((( V(I,J,K) ,I=0,N1),J=1,N2),K=0,N3)
      WRITE(NV) ((( W(I,J,K) ,I=0,N1),J=0,N2),K=1,N3)
      WRITE(NV) ((( P(I,J,K) ,I=1,N1M),J=1,N2M),K=1,N3M)
      IF (IHTRANS .EQ. 1) WRITE(NV) ((( T(I,J,K), I=1,N1M),J=1,N2M),K=1,N3M)
      CLOSE(NV)

      NV=NV+1

      RETURN
      END SUBROUTINE WRITEFIELD
!=======================================================================
!=======================================================================
      SUBROUTINE WRITEFTRTIME
!=======================================================================
!
!     Real time measurements for individual processes
!     Printed at every computational step
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      REAL*8     ::   FTRTIME1,FTRTIME2,FTRTIME3,FTRTIME4

      FTRTIME1=(SGSTIME_E(3)-SGSTIME_B(3)+SGSTIME_E(2)-SGSTIME_B(2)+SGSTIME_E(1)-SGSTIME_B(1))
      FTRTIME2=(RHSNLHSTIME_E(3)-RHSNLHSTIME_B(3)+RHSNLHSTIME_E(2)-RHSNLHSTIME_B(2)+RHSNLHSTIME_E(1)-RHSNLHSTIME_B(1))
      FTRTIME3=(POISSTIME_E(3)-POISSTIME_B(3)+POISSTIME_E(2)-POISSTIME_B(2)+POISSTIME_E(1)-POISSTIME_B(1))
      FTRTIME4=TIME_END-TIME_BEGIN

      WRITE(*,201) FTRTIME1
      WRITE(*,202) FTRTIME2
      WRITE(*,203) FTRTIME3
      WRITE(*,204) FTRTIME4

 201  FORMAT('TIME FOR SGS    : ',F12.3,' SECONDS')
 202  FORMAT('TIME FOR RHSnLHS: ',F12.3,' SECONDS')
 203  FORMAT('TIME FOR POISSON: ',F12.3,' SECONDS')
 204  FORMAT('TIME OF OPERTION: ',F12.3,' SECONDS')
      WRITE(2011,206)TIME,FTRTIME1,FTRTIME2,FTRTIME3,FTRTIME4
 206  FORMAT(F13.5,5F12.4)

      RETURN
      END
!=======================================================================
!=======================================================================
      SUBROUTINE FIELD_AVG
!=======================================================================
!     SUBROUTINE FOR FIELD AVERAGING
!     AVERAGED VARIABLES ARE DEFINED AT CELL CENTER
!     VARIABLES ARE LINEARLY INTERPOLATED
!     FROM STAGGERED GRID STRUCTURE TO DETERMINE CELL CENTER VALUES
!
!     No spatial averaging.
!
!     Print an averaged flow-field file, when MOD(NTIME,NPRIAVG)=0,
!     satisfying the averaged flow-field file printing interval (NPRIAVG set in lica.in).
!     tfn1='fav': prefix for average flow field
!     tname     : average field file name  ex) fav100000-110000
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,P,T,UAVG,VAVG,WAVG,UIUJAVG,PAVG,&
                                P2AVG,TAVG,T2AVG,VORAVG,VOR2AVG,SSAVG
      IMPLICIT NONE
      INTEGER*8     :: I,J,K,L
      REAL*8        :: UCC,VCC,WCC
      REAL*8        :: VG11,VG12,VG13,VG21,VG22,VG23,VG31,VG32,VG33
      REAL*8        :: UP,UM,VP,VM,WP,WM,SR(6),SRSR
      REAL*8        :: WX,WY,WZ
      REAL*8        :: TIMEEND
      INTEGER*8     :: IHISTEND
      INTEGER*8     :: idg1,idg2,idg3,idg4,idg5,idg6
      CHARACTER*16  :: tname
      CHARACTER*3   :: tfn1
      CHARACTER*6   :: tfn2,tfn3
      CHARACTER*1   :: tfnh

      tfnh='-'

!$OMP PARALLEL DO &
!$OMP private(UCC,VCC,WCC)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
       UCC=0.5*(U(I,J,K)+U(I+1,J,K))
       VCC=0.5*(V(I,J,K)+V(I,J+1,K))
       WCC=0.5*(W(I,J,K)+W(I,J,K+1))
         UAVG(I,J,K)     =UAVG(I,J,K)     +UCC         *DT
         VAVG(I,J,K)     =VAVG(I,J,K)     +VCC         *DT
         WAVG(I,J,K)     =WAVG(I,J,K)     +WCC         *DT
         UIUJAVG(I,J,K,1)=UIUJAVG(I,J,K,1)+UCC**2.     *DT
         UIUJAVG(I,J,K,2)=UIUJAVG(I,J,K,2)+UCC*VCC     *DT
         UIUJAVG(I,J,K,3)=UIUJAVG(I,J,K,3)+UCC*WCC     *DT
         UIUJAVG(I,J,K,4)=UIUJAVG(I,J,K,4)+VCC**2.     *DT
         UIUJAVG(I,J,K,5)=UIUJAVG(I,J,K,5)+VCC*WCC     *DT
         UIUJAVG(I,J,K,6)=UIUJAVG(I,J,K,6)+WCC**2.     *DT
         IF(IHTRANS.EQ.1) THEN
         TAVG(I,J,K)     =TAVG(I,J,K)     +T(I,J,K)    *DT
         T2AVG(I,J,K)    =T2AVG(I,J,K)    +T(I,J,K)**2.*DT
         ENDIF
         PAVG(I,J,K)     =PAVG(I,J,K)     +P(I,J,K)    *DT
         P2AVG(I,J,K)    =P2AVG(I,J,K)    +P(I,J,K)**2.*DT
      ENDDO
      ENDDO
      ENDDO

!$OMP PARALLEL DO &
!$OMP private(VG11,VG12,VG13,VG21,VG22,VG23,VG31,VG32,VG33) &
!$OMP private(UP,UM,VP,VM,WP,WM,SR,SRSR) &
!$OMP private(WX,WY,WZ)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M

       VG11=F2FXI(I)*(U(I+1,J,K)-U(I,J,K))
       UP=C2CYI(J+1)*0.25                                      &
          *(F2FY(J+1)*(U(I,J,K)+U(I+1,J,K))                    &
           +F2FY(J)*(U(I,J+1,K)+U(I+1,J+1,K)))*(1.-FIXJU(J))   &
          +0.5*(U(I,N2,K)+U(I+1,N2,K))*FIXJU(J)
       UM=C2CYI(J)*0.25                                        &
          *(F2FY(J)*(U(I,J-1,K)+U(I+1,J-1,K))                  &
           +F2FY(J-1)*(U(I,J,K)+U(I+1,J,K)))*(1.-FIXJL(J))     &
          +0.5*(U(I,0,K)+U(I+1,0,K))*FIXJL(J)
       VG12=F2FYI(J)*(UP-UM)

       UP=C2CZI(K+1)*0.25                                      &
          *(F2FZ(K+1)*(U(I,J,K)+U(I+1,J,K))                    &
           +F2FZ(K)*(U(I,J,K+1)+U(I+1,J,K+1)))*(1.-FIXKU(K))   &
          +0.5*(U(I,J,N3)+U(I+1,J,N3))*FIXKU(K)
       UM=C2CZI(K)*0.25                                        &
          *(F2FZ(K)*(U(I,J,K-1)+U(I+1,J,K-1))                  &
           +F2FZ(K-1)*(U(I,J,K)+U(I+1,J,K)))*(1.-FIXKL(K))     &
          +0.5*(U(I,J,0)+U(I+1,J,0))*FIXKL(K)
       VG13=F2FZI(K)*(UP-UM)

       VP=C2CXI(I+1)*0.25                                      &
          *(F2FX(I+1)*(V(I,J,K)+V(I,J+1,K))                    &
           +F2FX(I)*(V(I+1,J,K)+V(I+1,J+1,K)))*(1.-FIXIU(I))   &
          +0.5*(V(N1,J,K)+V(N1,J+1,K))*FIXIU(I)
       VM=C2CXI(I)*0.25                                        &
          *(F2FX(I)*(V(I-1,J,K)+V(I-1,J+1,K))                  &
           +F2FX(I-1)*(V(I,J,K)+V(I,J+1,K)))*(1.-FIXIL(I))     &
          +0.5*(V(0,J,K)+V(0,J+1,K))*FIXIL(I)
       VG21=F2FXI(I)*(VP-VM)

       VG22=F2FYI(J)*(V(I,J+1,K)-V(I,J,K))

       VP=C2CZI(K+1)*0.25                                      &
          *(F2FZ(K+1)*(V(I,J,K)+V(I,J+1,K))                    &
           +F2FZ(K)*(V(I,J,K+1)+V(I,J+1,K+1)))*(1.-FIXKU(K))   &
          +0.5*(V(I,J,N3)+V(I,J+1,N3))*FIXKU(K)
       VM=C2CZI(K)*0.25                                        &
          *(F2FZ(K)*(V(I,J,K-1)+V(I,J+1,K-1))                  &
           +F2FZ(K-1)*(V(I,J,K)+V(I,J+1,K)))*(1.-FIXKL(K))     &
          +0.5*(V(I,J,0)+V(I,J+1,0))*FIXKL(K)
       VG23=F2FZI(K)*(VP-VM)

       WP=C2CXI(I+1)*0.25                                      &
          *(F2FX(I+1)*(W(I,J,K)+W(I,J,K+1))                    &
           +F2FX(I)*(W(I+1,J,K)+W(I+1,J,K+1)))*(1.-FIXIU(I))   &
          +0.5*(W(N1,J,K)+W(N1,J,K+1))*FIXIU(I)
       WM=C2CXI(I)*0.25                                        &
          *(F2FX(I)*(W(I-1,J,K)+W(I-1,J,K+1))                  &
           +F2FX(I-1)*(W(I,J,K)+W(I,J,K+1)))*(1.-FIXIL(I))     &
          +0.5*(W(0,J,K)+W(0,J,K+1))*FIXIL(I)
       VG31=F2FXI(I)*(WP-WM)

       WP=C2CYI(J+1)*0.25                                      &
          *(F2FY(J+1)*(W(I,J,K)+W(I,J,K+1))                    &
           +F2FY(J)*(W(I,J+1,K)+W(I,J+1,K+1)))*(1.-FIXJU(J))   &
          +0.5*(W(I,N2,K)+W(I,N2,K+1))*FIXJU(J)
       WM=C2CYI(J)*0.25                                        &
          *(F2FY(J)*(W(I,J-1,K)+W(I,J-1,K+1))                  &
           +F2FY(J-1)*(W(I,J,K)+W(I,J,K+1)))*(1.-FIXJL(J))     &
          +0.5*(W(I,0,K)+W(I,0,K+1))*FIXJL(J)
       VG32=F2FYI(J)*(WP-WM)

       VG33=F2FZI(K)*(W(I,J,K+1)-W(I,J,K))

       SR(1)= VG11
       SR(2)= 0.5*(VG12+VG21)
       SR(3)= 0.5*(VG13+VG31)
       SR(4)= VG22
       SR(5)= 0.5*(VG23+VG32)
       SR(6)= VG33

       WX=VG32-VG23
       WY=VG13-VG31
       WZ=VG21-VG12

       VORAVG(I,J,K,1)=VORAVG(I,J,K,1)+WX*DT
       VORAVG(I,J,K,2)=VORAVG(I,J,K,2)+WY*DT
       VORAVG(I,J,K,3)=VORAVG(I,J,K,3)+WZ*DT

       VOR2AVG(I,J,K,1)=VOR2AVG(I,J,K,1)+WX*WX*DT
       VOR2AVG(I,J,K,2)=VOR2AVG(I,J,K,2)+WX*WY*DT
       VOR2AVG(I,J,K,3)=VOR2AVG(I,J,K,3)+WX*WZ*DT
       VOR2AVG(I,J,K,4)=VOR2AVG(I,J,K,4)+WY*WY*DT
       VOR2AVG(I,J,K,5)=VOR2AVG(I,J,K,5)+WY*WZ*DT
       VOR2AVG(I,J,K,6)=VOR2AVG(I,J,K,6)+WZ*WZ*DT

       SRSR= 2.*SR(2)**2.+SR(1)**2.+2.*SR(3)**2.+SR(4)**2.+2.*SR(5)**2.+SR(6)**2.
       SSAVG(I,J,K)= SSAVG(I,J,K)+SRSR*DT

      ENDDO
      ENDDO
      ENDDO

      IF (MOD(NTIME,NPRIAVG).EQ.0) THEN

       TIMEEND =TIME
       IHISTEND=IHIST

       tfn1='fav'
       idg1=ihistinit/100000
       idg2=(ihistinit-idg1*100000)/10000
       idg3=(ihistinit-idg1*100000-idg2*10000)/1000
       idg4=(ihistinit-idg1*100000-idg2*10000-idg3*1000)/100
       idg5=(ihistinit-idg1*100000-idg2*10000-idg3*1000-idg4*100)/10
       idg6=ihistinit-idg1*100000-idg2*10000-idg3*1000-idg4*100-idg5*10
       tfn2=char(idg1+48)//char(idg2+48)//char(idg3+48)//char(idg4+48)//char(idg5+48)//char(idg6+48)
       idg1=ihist/100000
       idg2=(ihist-idg1*100000)/10000
       idg3=(ihist-idg1*100000-idg2*10000)/1000
       idg4=(ihist-idg1*100000-idg2*10000-idg3*1000)/100
       idg5=(ihist-idg1*100000-idg2*10000-idg3*1000-idg4*100)/10
       idg6=ihist-idg1*100000-idg2*10000-idg3*1000-idg4*100-idg5*10
       tfn3=char(idg1+48)//char(idg2+48)//char(idg3+48)//char(idg4+48)//char(idg5+48)//char(idg6+48)
       tname=tfn1//tfn2//tfnh//tfn3

       OPEN(NAV,FILE='../output/field_avg/'//tname,FORM='UNFORMATTED')
       WRITE(NAV)N1M,N2M,N3M,RE
       WRITE(NAV)TIMEINIT,TIMEEND,IHISTINIT,IHISTEND
       WRITE(NAV)((( UAVG(I,J,K)     ,I=1,N1M),J=1,N2M),K=1,N3M)
       WRITE(NAV)((( VAVG(I,J,K)     ,I=1,N1M),J=1,N2M),K=1,N3M)
       WRITE(NAV)((( WAVG(I,J,K)     ,I=1,N1M),J=1,N2M),K=1,N3M)
       WRITE(NAV)((((UIUJAVG(I,J,K,L),I=1,N1M),J=1,N2M),K=1,N3M),L=1,6)
       WRITE(NAV)((( PAVG(I,J,K)     ,I=1,N1M),J=1,N2M),K=1,N3M)
       WRITE(NAV)((( P2AVG(I,J,K)    ,I=1,N1M),J=1,N2M),K=1,N3M)
       WRITE(NAV)((((VORAVG(I,J,K,L) ,I=1,N1M),J=1,N2M),K=1,N3M),L=1,3)
       WRITE(NAV)((((VOR2AVG(I,J,K,L)   ,I=1,N1M),J=1,N2M),K=1,N3M),L=1,6)
       WRITE(NAV)((( SSAVG(I,J,K)    ,I=1,N1M),J=1,N2M),K=1,N3M)
       IF (IHTRANS .EQ. 1) THEN
         WRITE(NAV)((( TAVG(I,J,K)     ,I=1,N1M),J=1,N2M),K=1,N3M)
         WRITE(NAV)((( T2AVG(I,J,K)    ,I=1,N1M),J=1,N2M),K=1,N3M)
       ENDIF
       CLOSE(NAV)

       NAV      = NAV+1
       UAVG     = 0.
       VAVG     = 0.
       WAVG     = 0.
       UIUJAVG  = 0.
       PAVG     = 0.
       P2AVG    = 0.
       IF (IHTRANS.EQ.1) THEN
       TAVG     = 0.
       T2AVG    = 0.
       ENDIF
       VORAVG   = 0.
       VOR2AVG     = 0.
       SSAVG    = 0.
       TIMEINIT = TIME
       IHISTINIT= IHIST

      ENDIF

      RETURN
      END
!=======================================================================