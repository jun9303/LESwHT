!==================== MODULE FOR IBM PREPROCESSING======================
      MODULE MOD_IBMPRE
!---- WRITE THE INPUTS
!-----GEOMETRY
       INTEGER*8, PARAMETER :: M1 = 97,M2 = 577,M3 = 641
       INTEGER*8, PARAMETER :: M1M=M1-1,M2M=M2-1,M3M=M3-1
       INTEGER*8 :: N1,N1M,N2,N2M,N3,N3M
       REAL*8    :: XL,YL,ZL,X(0:M1),Y(0:M2),Z(0:M3),XMP(0:M1),YMP(0:M2),ZMP(0:M3)
       INTEGER*8 :: IPV(M1M),JPV(M2M),KPV(M3M),IMV(M1M),JMV(M2M),KMV(M3M)
       REAL*8    :: FIXIL(M1M),FIXIU(M1M),FIXJL(M2M),FIXJU(M2M),FIXKL(M3M),FIXKU(M3M)
       REAL*8    :: SDX(0:M1),SDY(0:M2),SDZ(0:M3),VDX(M1),VDY(M2),VDZ(M3),SSDX(0:M1),SSDY(0:M2),SSDZ(0:M3),VVDX(M1),VVDY(M2),VVDZ(M3)

!------INPUT FILES OPTION!
!     : Refer to the routine READINPUTFILE and post_inst.in file for detailed explanation

       INTEGER*8    :: IBMON,ITOT,IUVWP,IWXYZ,IALL,ISKIP,JSKIP,KSKIP
       INTEGER*8    :: IPZ,IPX,NIP,NJP,NKP,NFLD,IND_FILM
       INTEGER*8    :: ISTART,IEND,JSTART,JEND,KSTART,KEND
       REAL*8       :: XIP(200),YJP(200),ZKP(200)
       CHARACTER*20 :: FLDNAME(100)
       CHARACTER*9  :: ftailijk
       CHARACTER*5  :: tfn1,tfn2
       CHARACTER*9  :: tfn3
       CHARACTER*10 :: tname1,tname2
       CHARACTER*14 :: tname3

!------VARIABLES
       INTEGER*8   :: IP(200),JP(200),KP(200)
       REAL*8      :: UAVG(M1M,M2M,M3M),VAVG(M1M,M2M,M3M),WAVG(M1M,M2M,M3M),UIUJAVG(M1M,M2M,M3M,6)
       REAL*8      :: PAVG(M1M,M2M,M3M),P2AVG(M1M,M2M,M3M)
       REAL*8      :: VORAVG(M1M,M2M,M3M,3),WRMS(M1M,M2M,M3M,6),SSAVG(M1M,M2M,M3M)

!------TIME
       REAL*8      :: TIMEINIT,TIMEEND
       REAL*8      :: IHISTINIT,IHISTEND

      END MODULE
!=======================================================================
!========================== MAIN PROCEDURE =============================
!=======================================================================
      PROGRAM POST_2D
      USE MOD_IBMPRE
      IMPLICIT NONE
      INTEGER*8    :: N,L
      CHARACTER*20 :: DUMMY
      INTEGER*8    :: IPOINT,JPOINT,KPOINT

      IPZ = 1
      IPX = 1

       CALL READINPUTFILE                             ! Read 'post_inst.in' file
       CALL GEOM('../output/grid/grid.dat')           ! Section for grid geometry

       CALL FINDIPJPKP                                ! Find indices of planes on which 2D fields will be printed

       DO L= 1,NFLD
          WRITE(*,*) ' '
          WRITE(*,108)
          WRITE(*,*)'WORKING ON ',FLDNAME(L)
          CALL PREFLD(FLDNAME(L))                     ! Read an input field file
          CALL DATAINIT
          CALL MAKEFHEAD                              ! Make head of an output field file

!======== Work on a plane and write a 2D output field
          DO N= 1,NIP                                 ! Work on a x-plane and write a 2D output field
           WRITE(*,*)'   WORKING ON X= ',XMP(IP(N))
           IPOINT=IP(N)
           CALL MAKEFTAIL(IPOINT,1)                   ! Make tail of an output field file
           CALL OUTPUTI(IPOINT)                       ! Write an output field file (center velocity, pressure, vorticity, vorticity magnitude, and lambda2)
          ENDDO

          DO N= 1,NJP
           WRITE(*,*)'   WORKING ON Y= ',YMP(JP(N))
           JPOINT=JP(N)
           CALL MAKEFTAIL(JPOINT,2)
           CALL OUTPUTJ(JPOINT)
          ENDDO

          DO N= 1,NKP
           WRITE(*,*)'   WORKING ON Z= ',ZMP(KP(N))
           KPOINT=KP(N)
           CALL MAKEFTAIL(KPOINT,3)
           CALL OUTPUTK(KPOINT)
          ENDDO
!=========

!======== Work on the whole 3D field
          IF (ITOT.EQ.1) THEN
           WRITE(*,*)'   WORKING ON 3D FIELD OUTPUT '
           CALL OUTPUT_3D
          ENDIF

          WRITE(*,108)
          IND_FILM=IND_FILM+1
      ENDDO
  108 FORMAT('====================================================================')

      STOP
      END

!=======================================================================
      SUBROUTINE MAKEFTAIL(IDEX,LL)
!=======================================================================
!
!     Write a tail of output file name
!     ex) ftailijk = _i048.dat if a working plane is a x-plane with
!         I = 48 i.e. if writing a 2d output field file on I = 48
!
!=======================================================================
      USE MOD_IBMPRE
      IMPLICIT NONE
      INTEGER*8    :: IDEX,LL
      INTEGER*8    :: I1,I2,I3
      CHARACTER*1  :: NN1,NN2,NN3
      CHARACTER*2  :: IJK

        I3=IDEX/100
        I2=IDEX/10-I3*10
        I1=IDEX-I3*100-I2*10
        NN3= CHAR(I3+48)
        NN2= CHAR(I2+48)
        NN1= CHAR(I1+48)

        IF (LL.EQ.1) THEN
          IJK='_i'
        ELSEIF (LL.EQ.2) THEN
          IJK='_j'
        ELSEIF (LL.EQ.3) THEN
          IJK='_k'
        ENDIF

        ftailijk=IJK//NN3(1:1)//NN2(1:1)//NN1(1:1)//'.dat'

      RETURN
      END

!=======================================================================
      SUBROUTINE MAKEFHEAD
!=======================================================================
!
!     Write a head of output file name
!     ex) If IND_FILM is 5 (in the post_inst.in), head names for the first
!         input fields are
!            2dfm_00005     (: 2d output for a certain plane)
!            3dfm_00005     (: 3d output)
!            2dfm_uni_00005 (: 2d output on a uniform grid)
!
!=======================================================================
      USE MOD_IBMPRE
      IMPLICIT NONE
      INTEGER*8    :: idg1,idg2,idg3,idg4,idg5

      idg1=IND_FILM/10000
      idg2=(IND_FILM-idg1*10000)/1000
      idg3=(IND_FILM-idg1*10000-idg2*1000)/100
      idg4=(IND_FILM-idg1*10000-idg2*1000-idg3*100)/10
      idg5=(IND_FILM-idg1*10000-idg2*1000-idg3*100-idg4*10)

      tfn1='2dav_'
      tfn2='3dav_'

      tname1=tfn1//char(idg1+48)//char(idg2+48)//char(idg3+48)//char(idg4+48)//char(idg5+48)
      tname2=tfn2//char(idg1+48)//char(idg2+48)//char(idg3+48)//char(idg4+48)//char(idg5+48)


      RETURN
      END

!=======================================================================
        SUBROUTINE READINPUTFILE
!=======================================================================
!
!     Read 'post_inst.in' file
!
!     Explanation for the 'post_inst.in': (Open the 'post_inst.in' and
!     read it together)
!========================     IBM OPTION      ==========================
! IBMON (IMMERSED BOUNDARY METHOD ON/OFF) : 0= OFF, 1= ON
!========================   2D POST OPTION    ==========================
! IUNIGRID (UNIFROM GRID DATA) : 0= OFF, 1= ON
! XYZ_position : 2D CUT FIELDS
!   1. The number of positions (0= OFF)
!   2. Position of cell where one want to print a 2D field
!    ex)
!   X_position
!   2
!   0.1
!   0.2
!     -> 2D flow fields at (the nearest points from) x = 0.1, 0.2
!=====================   3D POST OPTION    =============================
! ITOT  (3D fields)                  : 0=OFF, 1=ON
! ILD2  (Lambda2 criterion fields)   : 0=OFF, 1=ON
! IUVWP (Velocity & pressure fields) : 0=OFF, 1=ON
! IWXYZ (Vorticity fields)           : 0=OFF, 1=ON
! IALL  (Velocity, pressure, vorticity, lambda2 criterion) : 0=OFF, 1=ON
!
! ISKIP,JSKIP,KSKIP     : Grid interval in X,Y and Z directions
!  ex)
! ISKIP   JSKIP   KSKIP
! 1         3       1
!  -> Print all data in X and Z directions, and skip 2 data for every
!     3 cells in Y direction
!
! ISTART,JSTART,KSTART  : Initial grid indices in X,Y and Z directions
! IEND,JEND,KEND        : Final grid indices in X,Y and Z directions
! ANIFLD (NFLD)         : Number of fields
! INDFILM               : Initial number of field
!  ex)
! ANIFLD  IND_FILM
! 3       5
! fld010000
! fld020000
! fld022000
!  -> you load 3(ANIFLD) files 'fld010000','fld020000',and 'fld022000'
!  -> after post-processing,
!  file number of 'fld010000' becomes '2dfm_00005 ... .dat' ... etc
!  file number of 'fld020000' becomes '2dfm_00006 ... .dat' ... etc
!  file number of 'fld022000' becomes '2dfm_00007 ... .dat' ... etc
!=======================================================================
       USE MOD_IBMPRE
       IMPLICIT NONE
      INTEGER*8    :: N,L
      CHARACTER*20 :: DUMMY

      OPEN(2,FILE='post_avg.in')
       READ(2,*) DUMMY
       READ(2,*) IBMON

       READ(2,*) DUMMY
       READ(2,*) DUMMY
       READ(2,*) NIP
       DO N=1,NIP
       READ(2,*) XIP(N)
       ENDDO
       READ(2,*) DUMMY
       READ(2,*) NJP
       DO N=1,NJP
       READ(2,*) YJP(N)
       ENDDO
       READ(2,*) DUMMY
       READ(2,*) NKP
       DO N=1,NKP
       READ(2,*) ZKP(N)
       ENDDO
       READ(2,*) DUMMY
       READ(2,*) DUMMY
       READ(2,*) ITOT
       READ(2,*) DUMMY
       READ(2,*) IUVWP,IWXYZ,IALL
       READ(2,*) DUMMY
       READ(2,*) ISKIP,JSKIP,KSKIP
       READ(2,*) DUMMY,DUMMY
       READ(2,*) ISTART,IEND
       READ(2,*) DUMMY,DUMMY
       READ(2,*) JSTART,JEND
       READ(2,*) DUMMY,DUMMY
       READ(2,*) KSTART,KEND
       READ(2,*) DUMMY,DUMMY
       READ(2,*) NFLD,IND_FILM
       DO L=1,NFLD
       READ(2,*) FLDNAME(L)
       ENDDO
      CLOSE(2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      WRITE(*,101) NIP
      WRITE(*,102) NJP
      WRITE(*,103) NKP
      DO N=1,NIP
      WRITE(*,104) N,XIP(N)
      ENDDO
      DO N=1,NJP
      WRITE(*,105) N,YJP(N)
      ENDDO
      DO N=1,NKP
      WRITE(*,106) N,ZKP(N)
      ENDDO
      WRITE(*,107) ITOT
      WRITE(*,108) IUVWP,IWXYZ,IALL
      WRITE(*,109) ISKIP,JSKIP,KSKIP
      WRITE(*,110) ISTART,IEND
      WRITE(*,111) JSTART,JEND
      WRITE(*,112) KSTART,KEND
      WRITE(*,113) NFLD,IND_FILM

  101 FORMAT('# OF XPOSITION  = ',I5)
  102 FORMAT('# OF YPOSITION  = ',I5)
  103 FORMAT('# OF ZPOSITION  = ',I5)
  104 FORMAT('X_POSITION ',I3,' : ',F13.5)
  105 FORMAT('Y_POSITION ',I3,' : ',F13.5)
  106 FORMAT('Z_POSITION ',I3,' : ',F13.5)
  107 FORMAT('ITOT   =',I5)
  108 FORMAT('ILD2   =',I5,'  IUVWP =',I5,'  IWXYZ =',I5,'  IWXYZ =',I5)
  109 FORMAT('ISKIP  =',I5,'  JSKIP =',I5,'  KSKIP =',I5)
  110 FORMAT('ISTART =',I5,'  IEND =',I5)
  111 FORMAT('JSTART =',I5,'  JEND =',I5)
  112 FORMAT('KSTART =',I5,'  KEND =',I5)
  113 FORMAT('NFLD   =',I5,'  IND_FILM =',I5)
  115 FORMAT('IUNI_GRID =',I5)

        RETURN
        END

!=======================================================================
      SUBROUTINE DATAINIT
!=======================================================================
!
!     Calaulate cell center velocities, vorticities and lambda 2 values
!
!=======================================================================
       USE MOD_IBMPRE
       IMPLICIT NONE
       INTEGER*8  :: I,J,K,L
       REAL*8     :: FUNCBODY

      ! Let U,V,W = 0 if a given point is inside the body
        DO K= 1,N3M
        DO J= 1,N2M
        DO I= 1,N1M
          IF (FUNCBODY(X(I),YMP(J),ZMP(K),0) .LE. 1.E-10) THEN
           UAVG(I,J,K)= 0.
          ENDIF
          IF (FUNCBODY(XMP(I),Y(J),ZMP(K),0) .LE. 1.E-10) THEN
           VAVG(I,J,K)= 0.
          ENDIF
          IF (FUNCBODY(XMP(I),YMP(J),Z(K),0) .LE. 1.E-10) THEN
           WAVG(I,J,K)= 0.
          ENDIF
        ENDDO
        ENDDO
        ENDDO

      IF (IBMON.EQ.1) THEN
       DO K= 1,N3M
       DO J= 1,N2M
       DO I= 1,N1M

       ! Modify flow variables at (XMP(I),YMP(J),ZMP(K)) if (XMP(I),YMP(J),ZMP(K)) is inside the body
       IF (FUNCBODY(XMP(I),YMP(J),ZMP(K),0) .LE. 1.E-10) THEN
         UAVG(I,J,K)      =  0.
         VAVG(I,J,K)      =  0.
         WAVG(I,J,K)      =  0.
         DO L= 1, 6
          IF (L .LE. 3) VORAVG(I,J,K,L) = 0.
          UIUJAVG(I,J,K,L) = 0.
          WRMS(I,J,K,L)    = 0.
         ENDDO
         PAVG(I,J,K)      =  1000.
         P2AVG(I,J,K)     =  1000.
         SSAVG(I,J,K)     =  0.
       ENDIF

       ENDDO
       ENDDO
       ENDDO
      ENDIF

      RETURN
      END
!=======================================================================
      SUBROUTINE OUTPUTI(IPOINT)
!=======================================================================
!
!     Write a 2D flow field on a x-plane (I = IPOINT)
!
!=======================================================================
      USE MOD_IBMPRE
      IMPLICIT NONE
      INTEGER*8    :: IPOINT
      INTEGER*8    :: I,J,K,L
      INTEGER*8    :: II,JJ,KK
      REAL*8       :: FUNCBODY,URMS,VRMS,WWRMS,VORMAG


      II=IPOINT

      OPEN(102,FILE='../output/favpost/'//tname1//ftailijk)
      WRITE(102,135)'variables= "z","y","u","v","w"'//&
                    ',"urms","vrms","wrms","vormag","p","uxuy","uxuz"'
      WRITE(102,*)'ZONE T="ZONE1" , I=',N3M,', J=',N2M,', F=POINT'
      DO J=1,N2M
      DO K=1,N3M
        URMS = (UIUJAVG(II,J,K,1) - UAVG(II,J,K)**2)**0.5
        VRMS = (UIUJAVG(II,J,K,4) - VAVG(II,J,K)**2)**0.5
        WWRMS = (UIUJAVG(II,J,K,6) - WAVG(II,J,K)**2)**0.5
        VORMAG = (VORAVG(II,J,K,1)**2+VORAVG(II,J,K,2)**2+VORAVG(II,J,K,3)**2)**0.5
       WRITE(102,101) ZMP(K),YMP(J),UAVG(II,J,K),VAVG(II,J,K),WAVG(II,J,K)   &
                     ,URMS,VRMS,WWRMS,VORMAG,PAVG(II,J,K),UIUJAVG(II,J,K,2),UIUJAVG(II,J,K,3)
      ENDDO
      ENDDO
      CLOSE(102)
  101 format(2F14.6,21ES13.5)
  135 format(a200)

      RETURN
      END
!=======================================================================
      SUBROUTINE OUTPUTJ(JPOINT)
!=======================================================================
!
!     Write a 2D flow field on a y-plane (J = JPOINT)
!
!=======================================================================
      USE MOD_IBMPRE
      IMPLICIT NONE
      INTEGER*8    :: JPOINT
      INTEGER*8    :: I,J,K,L,M
      INTEGER*8    :: II,JJ,KK
      REAL*8       :: FUNCBODY,URMS,VRMS,WWRMS,VORMAG

      JJ=JPOINT

      OPEN(102,FILE='../output/favpost/'//tname1//ftailijk)
      WRITE(102,135)'variables= "x","z","u","v","w"'//&
                    ',"urms","vrms","wrms","vormag","p"'
      WRITE(102,*)'ZONE T="ZONE1", I=',N1M*5,', J=',N3M,', F=POINT'
      DO K=1,N3M
      DO M=1,5
      DO I=1,N1M
        URMS = (UIUJAVG(I,JJ,K,1) - UAVG(I,JJ,K)**2)**0.5
        VRMS = (UIUJAVG(I,JJ,K,4) - VAVG(I,JJ,K)**2)**0.5
        WWRMS = (UIUJAVG(I,JJ,K,6) - WAVG(I,JJ,K)**2)**0.5
        VORMAG = (VORAVG(I,JJ,K,1)**2+VORAVG(I,JJ,K,2)**2+VORAVG(I,JJ,K,3)**2)**0.5
       WRITE(102,101) ((M-3)*(XMP(N1M)-XMP(1))+XMP(I))*0.015702,ZMP(K)*0.015702,UAVG(I,JJ,K)*0.8367,VAVG(I,JJ,K)*0.8367,WAVG(I,JJ,K)*0.8367   &
                     ,URMS*0.8367,VRMS*0.8367,WWRMS*0.8367,VORMAG*0.8367/0.015702,PAVG(I,JJ,K)
      ENDDO
      ENDDO
      ENDDO
      CLOSE(102)
  101 format(2F14.6,21ES13.5)
  135 format(a200)

      RETURN
      END

!=======================================================================
      SUBROUTINE OUTPUTK(KPOINT)
!=======================================================================
!
!     Write a 2D flow field on a z-plane (K = KPOINT)
!
!=======================================================================
      USE MOD_IBMPRE
      IMPLICIT NONE
      INTEGER*8    :: KPOINT
      INTEGER*8    :: I,J,K,L
      INTEGER*8    :: II,JJ,KK
      REAL*8       :: FUNCBODY,URMS,VRMS,WWRMS,VORMAG

      KK=KPOINT

      OPEN(102,FILE='../output/favpost/'//tname1//ftailijk)
      WRITE(102,135)'variables= "x","y","u","v","w"'//&
                    ',"urms","vrms","wrms","vormag","p"'
      WRITE(102,*)'ZONE T="ZONE1", I=',N1M,', J=',N2M,', F=POINT'
      DO J=1,N2M
      DO I=1,N1M
        URMS = (UIUJAVG(I,J,KK,1) - UAVG(I,J,KK)**2)**0.5
        VRMS = (UIUJAVG(I,J,KK,4) - VAVG(I,J,KK)**2)**0.5
        WWRMS = (UIUJAVG(I,J,KK,6) - WAVG(I,J,KK)**2)**0.5
        VORMAG = (VORAVG(I,J,KK,1)**2+VORAVG(I,J,KK,2)**2+VORAVG(I,J,KK,3)**2)**0.5
       WRITE(102,101) XMP(I)*0.015702,YMP(J)*0.015702,UAVG(I,J,KK)*0.8367,VAVG(I,J,KK)*0.8367,WAVG(I,J,KK)*0.8367   &
                     ,URMS*0.8367,VRMS*0.8367,WWRMS*0.8367,VORMAG*0.8367/0.015702,PAVG(I,J,KK)
      ENDDO
      ENDDO
      CLOSE(102)
  101 format(2F14.6,21ES13.5)
  135 format(a200)

      RETURN
      END

!=======================================================================
      SUBROUTINE OUTPUT_3D
!=======================================================================
!
!     Write a 3D output field file
!
!=======================================================================
      USE MOD_IBMPRE
      IMPLICIT NONE
      INTEGER*8    :: I,J,K,L,M
      INTEGER*8    :: INUM,JNUM,KNUM
      REAL*8       :: FUNCBODY,URMS,VRMS,WWRMS,VORMAG

       INUM = 0
       JNUM = 0
       KNUM = 0
       DO I=ISTART,IEND,ISKIP
        INUM=INUM+1
       ENDDO
       DO J=JSTART,JEND,JSKIP
        JNUM=JNUM+1
       ENDDO
       DO K=KSTART,KEND,KSKIP
         KNUM=KNUM+1
       ENDDO


!====== (When IUVWP = 1) Write UVWP values only
      IF (IUVWP.EQ.1) THEN
      OPEN(101,FILE='../output/favpost/'//tname2//'_uvwp.dat')
      WRITE(101,110)'variables= "x","y","z","u","v","w"'     //&
                    ',"urms","vrms","wrms","vormag","p"'
      WRITE(101,*)'ZONE T="ZONE1" , I=',INUM,' J=',JNUM,' K=',KNUM,' F=POINT'
      DO K=KSTART,KEND,KSKIP
      DO J=JSTART,JEND,JSKIP
      DO I=ISTART,IEND,ISKIP
        URMS = (UIUJAVG(I,J,K,1) - UAVG(I,J,K)**2)**0.5
        VRMS = (UIUJAVG(I,J,K,4) - VAVG(I,J,K)**2)**0.5
        WWRMS = (UIUJAVG(I,J,K,6) - WAVG(I,J,K)**2)**0.5
        VORMAG = (VORAVG(I,J,K,1)**2+VORAVG(I,J,K,2)**2+VORAVG(I,J,K,3)**2)**0.5
      WRITE(101,100) XMP(I)*0.015702,YMP(J)*0.015702,ZMP(K)*0.015702,UAVG(I,J,K)*0.8367  &
                    ,VAVG(I,J,K)*0.8367,WAVG(I,J,K)*0.8367,URMS*0.8367,VRMS*0.8367,WWRMS*0.8367,VORMAG*0.8367/0.015702,PAVG(I,J,K)
      ENDDO
      ENDDO
      ENDDO
      CLOSE(101)
  100 format(3F14.6,8ES13.5)
      ENDIF
!======

!====== (When IWXYZ = 1) Write WXYZ values only
      IF (IWXYZ.EQ.1) THEN
      OPEN(102,FILE='../output/favpost/'//tname2//'_wxyz.dat')
      WRITE(102,110)'variables= "x","y","z","wx","wy","wz","wxwx","wxwy"'//&
                    ',"wxwz","wywy","wywz","wzwz","ss"'
      WRITE(102,*)'ZONE T="ZONE1" , I=',INUM,' J=',JNUM,' K=',KNUM,' F=POINT'
      DO K=KSTART,KEND,KSKIP
      DO J=JSTART,JEND,JSKIP
      DO I=ISTART,IEND,ISKIP
      WRITE(102,103) XMP(I),YMP(J),ZMP(K),(VORAVG(I,J,K,L),L=1,3) &
                    ,(WRMS(I,J,K,L),L=1,6),SSAVG(I,J,K)
      ENDDO
      ENDDO
      ENDDO
      CLOSE(102)
  103 format(3F8.3,10ES13.5)
      ENDIF
!======

!====== (When IALL = 1) Write all
      IF (IALL.EQ.1) THEN
      OPEN(103,FILE='../output/favpost/'//tname2//'_all.dat')
      WRITE(103,110)'variables= "x","y","z",u","v","w","p","p2","wx","wy"'//&
                    ',"wz","uxux","uxuy","uxuz","uyuy","uyuz","uzuz"'  //&
                    ',"wxwx","wxwy","wxwz","wywy","wywz","wzwz","ss"'
      WRITE(103,*)'ZONE T="ZONE1" , I=',INUM,' J=',JNUM,' K=',KNUM,' F=POINT'
      DO K=KSTART,KEND,KSKIP
      DO J=JSTART,JEND,JSKIP
      DO I=ISTART,IEND,ISKIP
      WRITE(103,105) XMP(I),YMP(J),ZMP(K)                           &
              ,UAVG(I,J,K),VAVG(I,J,K),WAVG(I,J,K)                  &
              ,PAVG(I,J,K),P2AVG(I,J,K),(VORAVG(I,J,K,L),L=1,3)     &
              ,(UIUJAVG(I,J,K,L),L=1,6),(WRMS(I,J,K,L),L=1,6)       &
              ,SSAVG(I,J,K)
      ENDDO
      ENDDO
      ENDDO
      CLOSE(103)
  105 format(3F8.3,21ES13.5)
      ENDIF
!======

  110 FORMAT(A200)

      RETURN
      END

!=======================================================================
      SUBROUTINE PREFLD(fileprevel)
!=======================================================================
!
!     Read an input field file
!
!=======================================================================
       USE MOD_IBMPRE
       IMPLICIT NONE
       CHARACTER*20 :: fileprevel
       INTEGER*8    :: I,J,K,L,N
       INTEGER*8    :: NN1,NN2,NN3
       REAL*8       :: RRE,TIMETOTAL
       REAL*8       :: UAVGT(M1M,M2M,M3M),VAVGT(M1M,M2M,M3M),WAVGT(M1M,M2M,M3M),UIUJAVGT(M1M,M2M,M3M,6)
       REAL*8       :: PAVGT(M1M,M2M,M3M),P2AVGT(M1M,M2M,M3M)
       REAL*8       :: VORAVGT(M1M,M2M,M3M,3),WRMST(M1M,M2M,M3M,6),SSAVGT(M1M,M2M,M3M)

       UAVG= 0.
       VAVG= 0.
       WAVG= 0.
       PAVG= 0.
       UIUJAVG= 0.
       P2AVG= 0.
       VORAVG= 0.
       WRMS= 0.
       SSAVG= 0.

       TIMETOTAL = 0.

       UAVGT= 0.
       VAVGT= 0.
       WAVGT= 0.
       PAVGT= 0.
       UIUJAVGT= 0.
       P2AVGT= 0.
       VORAVGT= 0.
       WRMST= 0.
       SSAVGT= 0.

      DO N = 0,4
      IF (N.EQ.0) THEN
        OPEN(12,FILE='../output/fav/fav000000-005000',FORM='UNFORMATTED')
        WRITE(*,*) '11111111111'
      ELSEIF (N.EQ.1) THEN
        OPEN(12,FILE='../output/fav/fav005000-010000',FORM='UNFORMATTED')
        WRITE(*,*) '22222222222'
      ELSEIF (N.EQ.2) THEN
        OPEN(12,FILE='../output/fav/fav010000-015000',FORM='UNFORMATTED')
        WRITE(*,*) '33333333333'
      ELSEIF (N.EQ.3) THEN
        OPEN(12,FILE='../output/fav/fav015000-020000',FORM='UNFORMATTED')
        WRITE(*,*) '44444444444'
      ELSEIF (N.EQ.4) THEN
        OPEN(12,FILE='../output/fav/fav020000-025000',FORM='UNFORMATTED')
        WRITE(*,*) '55555555555'

      ENDIF
!     dum for future use
      READ(12)NN1,NN2,NN3,RRE
      READ(12)TIMEINIT,TIMEEND,IHISTINIT,IHISTEND
      READ(12)((( UAVGT(I,J,K)     ,I=1,N1M),J=1,N2M),K=1,N3M)
      READ(12)((( VAVGT(I,J,K)     ,I=1,N1M),J=1,N2M),K=1,N3M)
      READ(12)((( WAVGT(I,J,K)     ,I=1,N1M),J=1,N2M),K=1,N3M)
      READ(12)((((UIUJAVGT(I,J,K,L),I=1,N1M),J=1,N2M),K=1,N3M),L=1,6)
      READ(12)((( PAVGT(I,J,K)     ,I=1,N1M),J=1,N2M),K=1,N3M)
      READ(12)((( P2AVGT(I,J,K)    ,I=1,N1M),J=1,N2M),K=1,N3M)
      READ(12)((((VORAVGT(I,J,K,L) ,I=1,N1M),J=1,N2M),K=1,N3M),L=1,3)
      READ(12)((((WRMST(I,J,K,L)   ,I=1,N1M),J=1,N2M),K=1,N3M),L=1,6)
      READ(12)((( SSAVGT(I,J,K)    ,I=1,N1M),J=1,N2M),K=1,N3M)

      UAVG = UAVG + UAVGT 
      VAVG = VAVG + VAVGT  
      WAVG = WAVG + WAVGT 
      PAVG = PAVG + PAVGT 
      UIUJAVG = UIUJAVG + UIUJAVGT 
      P2AVG = P2AVG + P2AVGT 
      VORAVG = VORAVG + VORAVGT 
      WRMS = WRMS + WRMST 
      SSAVG = SSAVG +SSAVGT 

      TIMETOTAL = TIMETOTAL + (TIMEEND-TIMEINIT)
      ENDDO

      UAVG = UAVG / TIMETOTAL
      VAVG = VAVG / TIMETOTAL
      WAVG = WAVG / TIMETOTAL
      PAVG = PAVG / TIMETOTAL
      UIUJAVG = UIUJAVG / TIMETOTAL
      P2AVG = P2AVG / TIMETOTAL
      VORAVG = VORAVG / TIMETOTAL
      WRMS = WRMS / TIMETOTAL
      SSAVG = SSAVG / TIMETOTAL

      WRITE(*,100)
      WRITE(*,101)
      WRITE(*,102) fileprevel
      WRITE(*,103) RRE
      WRITE(*,104) NN1,NN2,NN3
      WRITE(*,105) TIMEINIT,TIMEEND,IHISTINIT,IHISTEND,TIMETOTAL
      WRITE(*,*)TIMETOTAL
      WRITE(*,*)UAVG(1,1,1)

  100 FORMAT('----------- INITIAL FIELD INFORMATION -----------')
  101 FORMAT('INITIAL FIELD      : READING DONE')
  102 FORMAT('INITIAL FIELD NAME : ',A30)
  103 FORMAT('RE=',ES12.3)
  104 FORMAT('N1=',I10,'  N2=',I10,'  N3=',I10)
  105 FORMAT('TI=',F12.5,'  TF=',F12.5,'  HISTI= ',I10,'  HISTF=',I10,'  TOTAL=',F12.5)

      RETURN
      END

!=======================================================================
!     Section for Geometry
!=======================================================================
      SUBROUTINE GEOM(gridfile)
!=======================================================================
!
!     Read 'grid.dat' and define various variables related to grid
!
!     GEOM
!       K_BGPZ,I_BGPX: intial index in x-,z-directions
!     READGRID
!       N1,N2,N3   : number of grid in x-,y-,z-directions
!       N1M,N2M,N3M: values obtained by subtracting 1 from N1,N2,N3
!       XL,YL,ZL   : domain size
!       X,Y,Z      : cell face positions in each direction
!     INDICES
!       IPV,JPV,KPV: I,J,K => X,Y,Z and P => +1, IPV(I)=I+1
!       IMV,JMV,KMV:                    M => -1, IMV(I)=I-1
!     INDXFIX
!       FIXIL,FIXJL,FIXKL: lower boundary indicators, which are only 1
!                          at the lower (I,J,K=1) boundary
!       FIXIU,FIXJU,FIXKU: upper boundary indicators, which are only 1
!                          at the upper (I,J,K=Max) boundary
!     MESHES
!       SDX,SDY,SDZ: distance between cell face,   SDX(I)=X(I+1)-X(I)
!       VDX,VDY,VDZ: distance between cell center, VDX(I)=(SDX(I)+SDX(I-1))/2
!       SSDX,SSDY,SSDZ: reciprocal of SDX,SDY,SDZ for spatial discretization
!       VVDX,VVDY,VVDZ: reciprocal of VDX,VDY,VDZ for spatial discretization
!     PHYSPOS
!       XMP,YMP,ZMP: cell center positions,        XMP(I)=X(I)+SDX(I)/2
!
!-----------------------------------------------------------------------
      USE MOD_IBMPRE
      IMPLICIT NONE
      INTEGER*8    :: I,K
      CHARACTER*30 :: gridfile

       OPEN(11,FILE=gridfile)
       CALL READGRID
       CLOSE(11)

       CALL INDICES
       CALL INDXFIX
       CALL MESHES
       CALL PHYSPOS

       IF (IPZ .EQ. 1) THEN
         FIXKL(1)  = 0.
         FIXKU(N3M)= 0.
         SDZ(0)    = SDZ(N3M)
         SDZ(N3)   = SDZ(1)
         SSDZ(0)   = SSDZ(N3M)
         SSDZ(N3)  = SSDZ(1)
         VDZ(1)    = 0.5*(SDZ(1)+SDZ(N3M))
         VDZ(N3)   = 0.5*(SDZ(1)+SDZ(N3M))
         VVDZ(1)   = 1./VDZ(1)
         VVDZ(N3)  = 1./VDZ(N3)
         KMV(1)    = N3M
         KPV(N3M)  = 1
       ENDIF

       IF (IPX .EQ. 1) THEN
         FIXIL(1)  = 0.
         FIXIU(N1M)= 0.
         SDX(0)    = SDX(N1M)
         SDX(N1)   = SDX(1)
         SSDX(0)   = SSDX(N1M)
         SSDX(N1)  = SSDX(1)
         VDX(1)    = 0.5*(SDX(1)+SDX(N1M))
         VDX(N1)   = 0.5*(SDX(1)+SDX(N1M))
         VVDX(1)   = 1./VDX(1)
         VVDX(N1)  = 1./VDX(N1)
         IMV(1)    = N1M
         IPV(N1M)  = 1
       ENDIF

      RETURN
      END
!=======================================================================
      SUBROUTINE READGRID
!=======================================================================
!
!     Read grid information from 'grid.dat'
!
!=======================================================================
      USE MOD_IBMPRE
      IMPLICIT NONE
      INTEGER*8    :: I,J,K

      READ (11,*)  N1,N2,N3
      WRITE(*,101) N1,N2,N3
      READ (11,*)  XL,YL,ZL
      WRITE(*,102) XL,YL,ZL
 101  FORMAT('NX=',I12,'  NY=',I12,'  NZ=',I12)
 102  FORMAT('XL=',F12.4,'  YL=',F12.4,'  ZL=',F12.4)

      N1M=N1-1
      N2M=N2-1
      N3M=N3-1

      IF ((N1.GT.M1) .OR. (N2.GT.M2) .OR. (N3.GT.M3)) THEN
         PRINT*, 'ARRAY SIZE CAN NOT HANDLE THIS GRID.'
         STOP
      END IF

      READ(11,*) (X(I),I=1,N1)
      READ(11,*) (Y(J),J=1,N2)
      READ(11,*) (Z(K),K=1,N3)

      RETURN
      END
!=======================================================================
      SUBROUTINE INDICES
!=======================================================================
!
!     Calculate adjacent grid's indices
!
!=======================================================================
      USE MOD_IBMPRE
      IMPLICIT NONE
      INTEGER*8    :: I,J,K

!-----STREAMWISE DIRECTION
      DO 10 I=1,N1M
      IPV(I)=I+1
      IMV(I)=I-1
   10 CONTINUE

!-----NORMAL DIRECTION
      DO 20 J=1,N2M
      JPV(J)=J+1
      JMV(J)=J-1
   20 CONTINUE

!-----SPANWISE DIRECTION
      DO 30 K=1,N3M
      KPV(K)=K+1
      KMV(K)=K-1
   30 CONTINUE

      RETURN
      END

!=======================================================================
      SUBROUTINE INDXFIX
!=======================================================================
!
!     Variables FIXxx(ex FIXIU) are 0 if a point is an interior point,
!     and 1 if a point is a boundary point (These variables are used
!     in the routine VORNLAMBDA2)
!
!=======================================================================
      USE MOD_IBMPRE
      IMPLICIT NONE
      INTEGER*8    :: I,J,K

      DO 10 I=1,N1M
      FIXIL(I) = 0.
      FIXIU(I) = 0.
   10 CONTINUE
      FIXIL(1)  = 1.
      FIXIU(N1M)= 1.

      DO 20 J=1,N2M
      FIXJL(J) = 0.
      FIXJU(J) = 0.
   20 CONTINUE
      FIXJL(1)  = 1.
      FIXJU(N2M)= 1.

      DO 30 K=1,N3M
      FIXKL(K) = 0.
      FIXKU(K) = 0.
   30 CONTINUE
      FIXKL(1)  = 1.
      FIXKU(N3M)= 1.

      RETURN
      END

!=======================================================================
      SUBROUTINE MESHES
!=======================================================================
!
!     Calculate side length (SDX) and cell center length (VDX) of each cell
!
!=======================================================================
      USE MOD_IBMPRE
      IMPLICIT NONE
      INTEGER*8    :: I,J,K

      DO 10 I=1,N1M
   10 SDX(I) = X(I+1)-X(I)
      DO 20 I=2,N1M
   20 VDX(I) = 0.5*(SDX(I)+SDX(I-1))
      VDX(1) = 0.5*SDX(1)             ! LEFT BOUNDED COND.
      VDX(N1)= 0.5*SDX(N1M)          ! RIGHT BOUNDED COND.

      DO 30 I=1,N1M
   30 SSDX(I)= 1./SDX(I)
      DO 40 I=1,N1
   40 VVDX(I)= 1./VDX(I)

      DO 15 J=1,N2M
   15 SDY(J) = Y(J+1)-Y(J)
      DO 25 J=2,N2M
   25 VDY(J) = 0.5*(SDY(J)+SDY(J-1))
      VDY(1) = 0.5*SDY(1)             ! LOWER WALL BOUNDED COND.
      VDY(N2)= 0.5*SDY(N2M)          ! UPPER WALL BOUNDED COND.

      DO 35 J=1,N2M
   35 SSDY(J)= 1./SDY(J)
      DO 45 J=1,N2
   45 VVDY(J)= 1./VDY(J)


      DO 17 K=1,N3M
   17 SDZ(K)=Z(K+1)-Z(K)
      DO 27 K=2,N3M
   27 VDZ(K) = 0.5*(SDZ(K)+SDZ(K-1))
      VDZ(1) = 0.5*SDZ(1)             ! LOWER WALL BOUNDED COND.
      VDZ(N3)= 0.5*SDZ(N3M)          ! UPPER WALL BOUNDED COND.

      DO 37 K=1,N3M
   37 SSDZ(K)=1./SDZ(K)
      DO 47 K=1,N3
   47 VVDZ(K)=1./VDZ(K)

!-----set by zero
      SDX(0) = 0.
      SDX(M1)= 0.
      SDY(0) = 0.
      SDY(M2)= 0.
      SDZ(0) = 0.
      SDZ(M3)= 0.

      RETURN
      END

!=======================================================================
      SUBROUTINE PHYSPOS
!=======================================================================
!
!     Calculate midpoint coordinates of each cell
!
!=======================================================================
      USE MOD_IBMPRE
      IMPLICIT NONE
      INTEGER*8    :: I,J,K

      DO I=1,N1M
        XMP(I) = X(I)+0.5*SDX(I)
      ENDDO
      XMP(N1)= X(N1)
      XMP(0) = X(1)

      DO J=1,N2M
        YMP(J) = Y(J)+0.5*SDY(J)
      ENDDO
      YMP(N2)= Y(N2)
      YMP(0) = Y(1)

      DO K=1,N3M
        ZMP(K) = Z(K)+0.5*SDZ(K)
      ENDDO
      ZMP(N3)= Z(N3)
      ZMP(0) = Z(1)

      RETURN
      END

!=======================================================================
        SUBROUTINE FINDIPJPKP
!=======================================================================
!
!     Find indices which correspond to input points (X_position, Y_position,
!     Z_position of post_inst.in file)
!
!=======================================================================
      USE MOD_IBMPRE
      IMPLICIT NONE
      INTEGER*8    :: I,J,K


      DO I=1,NIP
        DO J=N1M,1,-1
           IF (XIP(I) .GE. XMP(J)) THEN
            IP(I)=J
            GOTO 111
           ENDIF
        ENDDO
 111  CONTINUE
      ENDDO

      DO I=1,NJP
        DO J=N2M,1,-1
           IF (YJP(I) .GE. YMP(J)) THEN
            JP(I)=J
            GOTO 112
           ENDIF
        ENDDO
 112  CONTINUE
      ENDDO

      DO I=1,NKP
        DO J=N3M,1,-1
           IF (ZKP(I) .GE. ZMP(J)) THEN
            KP(I)=J
            GOTO 113
           ENDIF
        ENDDO
 113  CONTINUE
      ENDDO

        RETURN
        END
!=======================================================================
!     End of the Section for Geometry
!=======================================================================

!=======================================================================
      SUBROUTINE SORTER(N,RA)
!=======================================================================
!
!     Sort the input array RA with size N
!
!=======================================================================
      REAL RA(3),RRA
      INTEGER I,IR,J,L

      IF(N.LE.1) GOTO 30
      L=N/2+1
      IR=N
   10 CONTINUE
      IF(L .GT. 1) THEN
            L=L-1
            RRA=RA(L)
         ELSE
            RRA=RA(IR)
            RA(IR)=RA(1)
            IR=IR-1
            IF(IR.EQ.1)THEN
                 RA(1)=RRA
                 RETURN
            END IF
         END IF
      I=L
      J=L+L
   20 IF(J.LE.IR) THEN
         IF(J.LT.IR) THEN
              IF(RA(J).GT.RA(J+1))J=J+1
         END IF
         IF(RRA.GT.RA(J)) THEN
              RA(I)=RA(J)
              I=J
              J=J+J
         ELSE
              J=IR+1
         END IF
      GOTO 20
      END IF
      RA(I)=RRA

      GOTO 10
   30 RETURN
      END

!=======================================================================
      SUBROUTINE jacobi(a,n,np,d,v,nrot)
!=======================================================================
!
!     Compute eigenvalues of a matrix (from Numerical recipe in FORTRAN)
!
!=======================================================================
      INTEGER n,np,nrot,NMAX
      REAL a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=500)
      INTEGER i,ip,iq,j
      REAL c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.
11      continue
        v(ip,ip)=1.
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))

14        continue
15      continue
        if(sm.eq.0.)return
        if(i.lt.4)then
          tresh=0.2*sm/n**2
        else
          tresh=0.
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq)))) then
              a(ip,iq)=0.
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h

              else
                theta=0.5*h/a(ip,iq)
                t=1./(abs(theta)+sqrt(1.+theta**2))
                if(theta.lt.0.)t=-t
              endif
              c=1./sqrt(1+t**2)
              s=t*c
              tau=s/(1.+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)

                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)

                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.
23      continue
24    continue
      pause 'too many iterations in jacobi'
      return
      END

!=======================================================================
        FUNCTION FUNCBODY(XX,YY,ZZ,IS)
!=======================================================================
!
!     FUNCBODY < 0 IF (XX,YY,ZZ) is inside a moving cylinder
!     FUNCBODY > 0 IF (XX,YY,ZZ) is outside a moving cylinder
!
!=======================================================================
      USE MOD_IBMPRE, ONLY : N3M, ZMP, N2M, YMP
      IMPLICIT NONE
      REAL*8     :: XX,YY,ZZ,XTEMP,YTEMP,ZTEMP
      REAL*8     :: FUNCBODY
      INTEGER*8  :: IS,I,J,K
      INTEGER*8  :: NZP
      REAL*8     :: RIBS,RIBH
      REAL*8     :: LOWDEL
      REAL*8     :: ZBACUP,YBACUP

      RIBS = 0.11
      RIBH = 0.096

        FUNCBODY = +  YY -  6.849234/1.085481 * ZZ
        IF (ZZ.GE.ZMP(N3M)) FUNCBODY = -1.
        IF (YY.GE.YMP(N2M)) FUNCBODY = -1.

        ! FUNCBODY = 1.
        ! IF((ABS(YY).LE.0.58   ) .AND. (ABS(ZZ).LE.3.46   )) FUNCBODY = -1.
        ! IF((ABS(YY).LE.0.5426 ) .AND. (ABS(ZZ).LE.3.42375) .AND. (YY.GE.-1.0852/6.8475*ZZ)) FUNCBODY = 1.

  !       DO J = 1,7
  !       XTEMP = XX - .82972*(J-4.)
  !       YTEMP = YY + 1.0852 / 2.
  !       ZTEMP = ZZ - 1.38182

  !       CALL ROTATE_YAXIS(XTEMP,YTEMP,ZTEMP,+45.)

  !       IF ((ABS(XTEMP).LE.0.082972/2.).AND.((YTEMP.GE.0.).AND.(YTEMP.LE.0.082972))&
  !           .AND.(ABS(ZTEMP).LE. 3.87339/2.)) THEN
  !            FUNCBODY = -1.
  !       ENDIF

  !       ENDDO

		! DO J = 1,7
  !       XTEMP = XX + 0.82795/2. - .829415*(J-4.)
  !       YTEMP = YY
  !       ZTEMP = ZZ

  !       CALL ROTATE_XAXIS(XTEMP,YTEMP,ZTEMP,-9.00541)

  !       ZTEMP = ZTEMP + 1.38132

  !       CALL ROTATE_YAXIS(XTEMP,YTEMP,ZTEMP,-45.)

  !       IF ((ABS(XTEMP).LE.0.082972/2.).AND.((YTEMP.LE.0.05).AND.(YTEMP.GE.-0.082972))&
  !           .AND.(ABS(ZTEMP).LE. 3.87339/2.)) THEN
  !            FUNCBODY = -1.
  !       ENDIF

  !       ENDDO   

      RETURN
      END

!ADDED BY SANGJOON LEE
! **********************************************************************
       SUBROUTINE ROTATE_ZAXIS(XR,YR,ZR,THETA)
! **********************************************************************
!      ROTATE X AND Y COORDINATE TO COUNTERCLOCKWISE DIRECTION 

      IMPLICIT NONE
      REAL*8       :: XR,YR,ZR,THETA
      REAL*8       :: XTEMPO,YTEMPO

      XTEMPO=XR
      YTEMPO=YR

      XR=COS(THETA*3.14159265/180.)*XTEMPO-SIN(THETA*3.14159265/180.)&
         *YTEMPO
      YR=SIN(THETA*3.14159265/180.)*XTEMPO+COS(THETA*3.14159265/180.)&
         *YTEMPO
      ZR=ZR

      RETURN
      END

! **********************************************************************
       SUBROUTINE ROTATE_YAXIS(XR,YR,ZR,THETA)
! **********************************************************************
!      ROTATE X AND Z COORDINATE TO COUNTERCLOCKWISE DIRECTION 

      IMPLICIT NONE
      REAL*8       :: XR,YR,ZR,THETA
      REAL*8       :: XTEMPO,ZTEMPO

      XTEMPO=XR
      ZTEMPO=ZR

      XR=COS(THETA*3.14159265/180.)*XTEMPO-SIN(THETA*3.14159265/180.)&
         *ZTEMPO
      ZR=-SIN(THETA*3.14159265/180.)*XTEMPO-COS(THETA*3.14159265/180.)&
         *ZTEMPO
      YR=YR

      RETURN
      END

! **********************************************************************
       SUBROUTINE ROTATE_XAXIS(XR,YR,ZR,THETA)
! **********************************************************************
!      ROTATE Y AND Z COORDINATE TO COUNTERCLOCKWISE DIRECTION 

      IMPLICIT NONE
      REAL*8       :: XR,YR,ZR,THETA
      REAL*8       :: YTEMPO,ZTEMPO

      YTEMPO=YR
      ZTEMPO=ZR

      ZR=-COS(THETA*3.14159265/180.)*ZTEMPO+SIN(THETA*3.14159265/180.)&
         *YTEMPO
      YR=SIN(THETA*3.14159265/180.)*ZTEMPO+COS(THETA*3.14159265/180.)&
         *YTEMPO
      XR=XR

      RETURN
      END
