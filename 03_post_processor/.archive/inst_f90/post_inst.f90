!==========================================
MODULE MOD_POSTINST
!---- WRITE THE INPUTS
!-----GEOMETRY
  INTEGER*8, PARAMETER :: M1 = 641, M2 = 641, M3 = 641
  INTEGER*8, PARAMETER :: M1M = M1 - 1, M2M = M2 - 1, M3M = M3 - 1
  INTEGER*8 :: N1, N1M, N2, N2M, N3, N3M
  REAL*8    :: XL, YL, ZL, X(0:M1), Y(0:M2), Z(0:M3), XMP(0:M1), YMP(0:M2), ZMP(0:M3), RE, PR, GR
  INTEGER*8 :: IPV(M1M), JPV(M2M), KPV(M3M), IMV(M1M), JMV(M2M), KMV(M3M)
  REAL*8    :: FIXIL(M1M), FIXIU(M1M), FIXJL(M2M), FIXJU(M2M), FIXKL(M3M), FIXKU(M3M)
  REAL*8    :: SDX(0:M1),SDY(0:M2),SDZ(0:M3),VDX(M1),VDY(M2),VDZ(M3),SSDX(0:M1),SSDY(0:M2),SSDZ(0:M3),VVDX(M1),VVDY(M2),VVDZ(M3)
!------INPUT FILES OPTION!
      INTEGER*8    :: IBMON,IHTRANS,IOUTFMT,IUNIGRID,XPRDIC,YPRDIC,ZPRDIC,BC_XBTM, BC_XTOP, BC_YBTM, BC_YTOP, BC_ZBTM, BC_ZTOP,ITOT,ILD2,IUVWP,IWXYZ,IALL,ISKIP,JSKIP,KSKIP
  INTEGER*8    :: NIP, NJP, NKP, NFLD, IND_FILM, ICH, ITEMP
  INTEGER*8    :: ISTART, IEND, JSTART, JEND, KSTART, KEND
  REAL*8       :: XIP(20), YJP(20), ZKP(20)
  CHARACTER*100:: FLDNAME(100)
  CHARACTER*9  :: ftailijk
  CHARACTER*5  :: tfn1, tfn2
  CHARACTER*9  :: tfn3
  CHARACTER*10 :: tname1, tname2
  CHARACTER*14 :: tname3

!------VARIABLES
  INTEGER*8   :: IP(20), JP(20), KP(20), IUNI(M1M), JUNI(M2M), KUNI(M3M), INUMUNI, JNUMUNI, KNUMUNI
  REAL*8      :: XUNI(M1), YUNI(M2), ZUNI(M3), FACUNIX(M1M), FACUNIY(M2M), FACUNIZ(M3M)
  REAL*8      :: U(0:M1, 0:M2, 0:M3), V(0:M1, 0:M2, 0:M3), W(0:M1, 0:M2, 0:M3), P(0:M1, 0:M2, 0:M3), T(0:M1, 0:M2, 0:M3)
  REAL*8      :: UC(M1M, M2M, M3M), VC(M1M, M2M, M3M), WC(M1M, M2M, M3M)
  REAL*8      :: VLAMBDA2(M1M, M2M, M3M), VORX(0:M1, 0:M2, 0:M3), VORY(0:M1, 0:M2, 0:M3), VORZ(0:M1, 0:M2, 0:M3)

END MODULE
!=======================================================================
!========================== MAIN PROCEDURE =============================
!=======================================================================
PROGRAM POST_2D
  USE MOD_POSTINST
  IMPLICIT NONE
  INTEGER*8    :: N, L
  CHARACTER*20 :: DUMMY
  INTEGER*8    :: IPOINT, JPOINT, KPOINT

  XPRDIC = 0
  ZPRDIC = 1

  CALL READINPUTFILE
  IF (NFLD .GE. 1) THEN
    CALL READ_PERIODIC_FLAGS(FLDNAME(1))
  END IF
  CALL GEOM('../output/grid/grid.dat')
  CALL ADJUST_POST_BOUNDS

  CALL FINDIPJPKP
  IF (IUNIGRID .EQ. 1) CALL GRIDUNI

  DO L = 1, NFLD
    WRITE (*, *) ' '
    WRITE (*, 108)
    WRITE (*, *) 'WORKING ON ', FLDNAME(L)
    CALL PREFLD(FLDNAME(L))
    CALL VORNLAMBDA2
    CALL MAKEFHEAD

    DO N = 1, NIP
      WRITE (*, *) '   WORKING ON X= ', XMP(IP(N))
      IPOINT = IP(N)
      CALL MAKEFTAIL(IPOINT, 1)
      CALL OUTPUTI(IPOINT)
    END DO

    DO N = 1, NJP
      WRITE (*, *) '   WORKING ON Y= ', YMP(JP(N))
      JPOINT = JP(N)
      CALL MAKEFTAIL(JPOINT, 2)
      CALL OUTPUTJ(JPOINT)
    END DO

    DO N = 1, NKP
      WRITE (*, *) '   WORKING ON Z= ', ZMP(KP(N))
      KPOINT = KP(N)
      CALL MAKEFTAIL(KPOINT, 3)
      CALL OUTPUTK(KPOINT)
    END DO

    IF (ITOT .EQ. 1) THEN
      WRITE (*, *) '   WORKING ON 3D FIELD OUTPUT '
      CALL OUTPUT_3D
    END IF

    WRITE (*, 108)
    IND_FILM = IND_FILM + 1
  END DO
108 FORMAT('====================================================================')

  STOP
END

!=======================================================================
SUBROUTINE READ_PERIODIC_FLAGS(fileprevel)
  !=======================================================================
  USE MOD_POSTINST
  IMPLICIT NONE
  CHARACTER*20 :: fileprevel
  INTEGER*8    :: NN1, NN2, NN3, IHIST, M
  REAL*8       :: RRE, TPR, TGR, TIME, DT

  OPEN (12, FILE='../output/field/'//fileprevel, FORM='UNFORMATTED')
  READ (12) NN1, NN2, NN3, RRE, TPR, TGR
  READ (12) IHIST, M, TIME, DT
  READ (12) XPRDIC, YPRDIC, ZPRDIC
  CLOSE (12)

  RETURN
END

!=======================================================================
SUBROUTINE MAKEFTAIL(IDEX, LL)
!=======================================================================
  USE MOD_POSTINST
  IMPLICIT NONE
  INTEGER*8    :: IDEX, LL
  INTEGER*8    :: I1, I2, I3
  CHARACTER*1  :: NN1, NN2, NN3
  CHARACTER*2  :: IJK

  I3 = IDEX/100
  I2 = IDEX/10 - I3*10
  I1 = IDEX - I3*100 - I2*10
  NN3 = CHAR(I3 + 48)
  NN2 = CHAR(I2 + 48)
  NN1 = CHAR(I1 + 48)

  IF (LL .EQ. 1) THEN
    IJK = '_i'
  ELSEIF (LL .EQ. 2) THEN
    IJK = '_j'
  ELSEIF (LL .EQ. 3) THEN
    IJK = '_k'
  END IF

  ftailijk = IJK//NN3(1:1)//NN2(1:1)//NN1(1:1)//'.vtk'

  RETURN
END

!=======================================================================
SUBROUTINE MAKEFHEAD
!=======================================================================
  USE MOD_POSTINST
  IMPLICIT NONE
  INTEGER*8    :: idg1, idg2, idg3, idg4, idg5

  idg1 = IND_FILM/10000
  idg2 = (IND_FILM - idg1*10000)/1000
  idg3 = (IND_FILM - idg1*10000 - idg2*1000)/100
  idg4 = (IND_FILM - idg1*10000 - idg2*1000 - idg3*100)/10
  idg5 = (IND_FILM - idg1*10000 - idg2*1000 - idg3*100 - idg4*10)

  tfn1 = '2dfm_'
  tfn2 = '3dfm_'
  tfn3 = '2dfm_uni_'

  tname1 = tfn1//char(idg1 + 48)//char(idg2 + 48)//char(idg3 + 48)//char(idg4 + 48)//char(idg5 + 48)
  tname2 = tfn2//char(idg1 + 48)//char(idg2 + 48)//char(idg3 + 48)//char(idg4 + 48)//char(idg5 + 48)
  tname3 = tfn3//char(idg1 + 48)//char(idg2 + 48)//char(idg3 + 48)//char(idg4 + 48)//char(idg5 + 48)

  RETURN
END

!=======================================================================
SUBROUTINE READINPUTFILE
!=======================================================================
  USE MOD_POSTINST
  IMPLICIT NONE
  INTEGER*8    :: N, L
  CHARACTER*20 :: DUMMY

  OPEN (2, FILE='post_inst.in')
  READ (2, *) DUMMY
  READ (2, *) IBMON, IHTRANS
  READ (2, *) DUMMY
  READ (2, *) IOUTFMT

  READ (2, *) DUMMY
  READ (2, *) DUMMY
  READ (2, *) IUNIGRID
  READ (2, *) DUMMY
  READ (2, *) NIP
  DO N = 1, NIP
    READ (2, *) XIP(N)
  END DO
  READ (2, *) DUMMY
  READ (2, *) NJP
  DO N = 1, NJP
    READ (2, *) YJP(N)
  END DO
  READ (2, *) DUMMY
  READ (2, *) NKP
  DO N = 1, NKP
    READ (2, *) ZKP(N)
  END DO
  READ (2, *) DUMMY
  READ (2, *) DUMMY
  READ (2, *) ITOT
  READ (2, *) DUMMY
  READ (2, *) ILD2, IUVWP, IWXYZ, IALL
  READ (2, *) DUMMY
  READ (2, *) ISKIP, JSKIP, KSKIP
  READ (2, *) DUMMY, DUMMY
  READ (2, *) ISTART, IEND
  READ (2, *) DUMMY, DUMMY
  READ (2, *) JSTART, JEND
  READ (2, *) DUMMY, DUMMY
  READ (2, *) KSTART, KEND
  READ (2, *) DUMMY, DUMMY
  READ (2, *) NFLD, IND_FILM
  DO L = 1, NFLD
    READ (2, *) FLDNAME(L)
  END DO
  CLOSE (2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  WRITE (*, 101) NIP
  WRITE (*, 102) NJP
  WRITE (*, 103) NKP
  DO N = 1, NIP
    WRITE (*, 104) N, XIP(N)
  END DO
  DO N = 1, NJP
    WRITE (*, 105) N, YJP(N)
  END DO
  DO N = 1, NKP
    WRITE (*, 106) N, ZKP(N)
  END DO
  WRITE (*, 107) ITOT
  WRITE (*, 108) ILD2, IUVWP, IWXYZ, IALL
  WRITE (*, 109) ISKIP, JSKIP, KSKIP
  WRITE (*, 110) ISTART, IEND
  WRITE (*, 111) JSTART, JEND
  WRITE (*, 112) KSTART, KEND
  WRITE (*, 113) NFLD, IND_FILM
  WRITE (*, 115) IUNIGRID
  WRITE (*, 116) IOUTFMT

101 FORMAT('# OF XPOSITION  = ', I5)
102 FORMAT('# OF YPOSITION  = ', I5)
103 FORMAT('# OF ZPOSITION  = ', I5)
104 FORMAT('X_POSITION ', I3, ' : ', F13.5)
105 FORMAT('Y_POSITION ', I3, ' : ', F13.5)
106 FORMAT('Z_POSITION ', I3, ' : ', F13.5)
107 FORMAT('ITOT   =', I5)
108 FORMAT('ILD2   =', I5, '  IUVWP =', I5, '  IWXYZ =', I5, '  IALL  =', I5)
109 FORMAT('ISKIP  =', I5, '  JSKIP =', I5, '  KSKIP =', I5)
110 FORMAT('ISTART =', I5, '  IEND =', I5)
111 FORMAT('JSTART =', I5, '  JEND =', I5)
112 FORMAT('KSTART =', I5, '  KEND =', I5)
113 FORMAT('NFLD   =', I5, '  IND_FILM =', I5)
115 FORMAT('IUNI_GRID =', I5)
116 FORMAT('IOUTFMT =', I5, '  (0:VTK, 1:TEC, 2:BOTH)')

  RETURN
END

!=======================================================================
SUBROUTINE ADJUST_POST_BOUNDS
!=======================================================================
  USE MOD_POSTINST
  IMPLICIT NONE
  INTEGER*8 :: TMP

  IF (ISKIP .LE. 0) ISKIP = 1
  IF (JSKIP .LE. 0) JSKIP = 1
  IF (KSKIP .LE. 0) KSKIP = 1

  IF (ISTART .LE. 0) ISTART = 1
  IF (JSTART .LE. 0) JSTART = 1
  IF (KSTART .LE. 0) KSTART = 1

  IF (IEND .LE. 0) IEND = N1M
  IF (JEND .LE. 0) JEND = N2M
  IF (KEND .LE. 0) KEND = N3M

  ISTART = MAX(1_8, MIN(ISTART, N1M))
  IEND = MAX(1_8, MIN(IEND, N1M))
  JSTART = MAX(1_8, MIN(JSTART, N2M))
  JEND = MAX(1_8, MIN(JEND, N2M))
  KSTART = MAX(1_8, MIN(KSTART, N3M))
  KEND = MAX(1_8, MIN(KEND, N3M))

  IF (ISTART .GT. IEND) THEN
    TMP = ISTART
    ISTART = IEND
    IEND = TMP
  END IF
  IF (JSTART .GT. JEND) THEN
    TMP = JSTART
    JSTART = JEND
    JEND = TMP
  END IF
  IF (KSTART .GT. KEND) THEN
    TMP = KSTART
    KSTART = KEND
    KEND = TMP
  END IF

  WRITE (*, '(A,I5,A,I5)') 'AUTO ISTART =', ISTART, '  IEND =', IEND
  WRITE (*, '(A,I5,A,I5)') 'AUTO JSTART =', JSTART, '  JEND =', JEND
  WRITE (*, '(A,I5,A,I5)') 'AUTO KSTART =', KSTART, '  KEND =', KEND

  RETURN
END

!=======================================================================
SUBROUTINE VORNLAMBDA2
!=======================================================================
  USE MOD_POSTINST
  IMPLICIT NONE
  INTEGER*8  :: I, J, K, II, JJ
  REAL*8     :: FUNCBODY
  REAL*8     :: VG11, VG12, VG13, VG21, VG22, VG23, VG31, VG32, VG33
  REAL*8     :: UP, UM, VP, VM, WP, WM, WX, WY, WZ
  REAL*8     :: SS(3, 3), OM(3, 3), TMAT(M1M, M2M, 3, 3), DD(3), EVEC(3, 3), TMP(3, 3)
  INTEGER*8  :: ITER

  DO K = 1, N3M
    DO J = 1, N2M
      DO I = 1, N1M
        IF (FUNCBODY(X(I), YMP(J), ZMP(K), 0) .LE. 1.E-10) THEN
          U(I, J, K) = 0.
        END IF
        IF (FUNCBODY(XMP(I), Y(J), ZMP(K), 0) .LE. 1.E-10) THEN
          V(I, J, K) = 0.
        END IF
        IF (FUNCBODY(XMP(I), YMP(J), Z(K), 0) .LE. 1.E-10) THEN
          W(I, J, K) = 0.
        END IF
      END DO
    END DO
  END DO

  DO K = 1, N3M
    DO J = 1, N2M
      DO I = 1, N1M
        UC(I, J, K) = 0.5*(U(I, J, K) + U(I + 1, J, K))
        VC(I, J, K) = 0.5*(V(I, J, K) + V(I, J + 1, K))
        WC(I, J, K) = 0.5*(W(I, J, K) + W(I, J, K + 1))
      END DO
    END DO
  END DO

  DO K = 1, N3M
    DO J = 1, N2M
      DO I = 1, N1M

        VG11 = SSDX(I)*(U(I + 1, J, K) - U(I, J, K))

        UP = VVDY(J + 1)*0.25 &
             *(SDY(J + 1)*(U(I, J, K) + U(I + 1, J, K)) &
               + SDY(J)*(U(I, J + 1, K) + U(I + 1, J + 1, K)))*(1.-FIXJU(J)) &
             + 0.5*(U(I, N2, K) + U(I + 1, N2, K))*FIXJU(J)
        UM = VVDY(J)*0.25 &
             *(SDY(J)*(U(I, J - 1, K) + U(I + 1, J - 1, K)) &
               + SDY(J - 1)*(U(I, J, K) + U(I + 1, J, K)))*(1.-FIXJL(J)) &
             + 0.5*(U(I, 0, K) + U(I + 1, 0, K))*FIXJL(J)
        VG12 = SSDY(J)*(UP - UM)

        UP = VVDZ(K + 1)*0.25 &
             *(SDZ(K + 1)*(U(I, J, K) + U(I + 1, J, K)) &
               + SDZ(K)*(U(I, J, K + 1) + U(I + 1, J, K + 1)))*(1.-FIXKU(K)) &
             + 0.5*(U(I, J, N3) + U(I + 1, J, N3))*FIXKU(K)
        UM = VVDZ(K)*0.25 &
             *(SDZ(K)*(U(I, J, K - 1) + U(I + 1, J, K - 1)) &
               + SDZ(K - 1)*(U(I, J, K) + U(I + 1, J, K)))*(1.-FIXKL(K)) &
             + 0.5*(U(I, J, 0) + U(I + 1, J, 0))*FIXKL(K)
        VG13 = SSDZ(K)*(UP - UM)

        VP = VVDX(I + 1)*0.25 &
             *(SDX(I + 1)*(V(I, J, K) + V(I, J + 1, K)) &
               + SDX(I)*(V(I + 1, J, K) + V(I + 1, J + 1, K)))*(1.-FIXIU(I)) &
             + 0.5*(V(N1, J, K) + V(N1, J + 1, K))*FIXIU(I)
        VM = VVDX(I)*0.25 &
             *(SDX(I)*(V(I - 1, J, K) + V(I - 1, J + 1, K)) &
               + SDX(I - 1)*(V(I, J, K) + V(I, J + 1, K)))*(1.-FIXIL(I)) &
             + 0.5*(V(0, J, K) + V(0, J + 1, K))*FIXIL(I)
        VG21 = SSDX(I)*(VP - VM)

        VG22 = SSDY(J)*(V(I, J + 1, K) - V(I, J, K))

        VP = VVDZ(K + 1)*0.25 &
             *(SDZ(K + 1)*(V(I, J, K) + V(I, J + 1, K)) &
               + SDZ(K)*(V(I, J, K + 1) + V(I, J + 1, K + 1)))*(1.-FIXKU(K)) &
             + 0.5*(V(I, J, N3) + V(I, J + 1, N3))*FIXKU(K)
        VM = VVDZ(K)*0.25 &
             *(SDZ(K)*(V(I, J, K - 1) + V(I, J + 1, K - 1)) &
               + SDZ(K - 1)*(V(I, J, K) + V(I, J + 1, K)))*(1.-FIXKL(K)) &
             + 0.5*(V(I, J, 0) + V(I, J + 1, 0))*FIXKL(K)
        VG23 = SSDZ(K)*(VP - VM)

        WP = VVDX(I + 1)*0.25 &
             *(SDX(I + 1)*(W(I, J, K) + W(I, J, K + 1)) &
               + SDX(I)*(W(I + 1, J, K) + W(I + 1, J, K + 1)))*(1.-FIXIU(I)) &
             + 0.5*(W(N1, J, K) + W(N1, J, K + 1))*FIXIU(I)
        WM = VVDX(I)*0.25 &
             *(SDX(I)*(W(I - 1, J, K) + W(I - 1, J, K + 1)) &
               + SDX(I - 1)*(W(I, J, K) + W(I, J, K + 1)))*(1.-FIXIL(I)) &
             + 0.5*(W(0, J, K) + W(0, J, K + 1))*FIXIL(I)
        VG31 = SSDX(I)*(WP - WM)

        WP = VVDY(J + 1)*0.25 &
             *(SDY(J + 1)*(W(I, J, K) + W(I, J, K + 1)) &
               + SDY(J)*(W(I, J + 1, K) + W(I, J + 1, K + 1)))*(1.-FIXJU(J)) &
             + 0.5*(W(I, N2, K) + W(I, N2, K + 1))*FIXJU(J)
        WM = VVDY(J)*0.25 &
             *(SDY(J)*(W(I, J - 1, K) + W(I, J - 1, K + 1)) &
               + SDY(J - 1)*(W(I, J, K) + W(I, J, K + 1)))*(1.-FIXJL(J)) &
             + 0.5*(W(I, 0, K) + W(I, 0, K + 1))*FIXJL(J)

        VG32 = SSDY(J)*(WP - WM)

        VG33 = SSDZ(K)*(W(I, J, K + 1) - W(I, J, K))

        WX = VG32 - VG23
        WY = VG13 - VG31
        WZ = VG21 - VG12
        VORX(I, J, K) = WX
        VORY(I, J, K) = WY
        VORZ(I, J, K) = WZ

        SS(1, 1) = VG11
        SS(1, 2) = 0.5*(VG12 + VG21)
        SS(1, 3) = 0.5*(VG13 + VG31)
        SS(2, 1) = 0.5*(VG21 + VG12)
        SS(2, 2) = VG22
        SS(2, 3) = 0.5*(VG23 + VG32)
        SS(3, 1) = 0.5*(VG31 + VG13)
        SS(3, 2) = 0.5*(VG32 + VG23)
        SS(3, 3) = VG33
        OM(1, 1) = 0.
        OM(1, 2) = 0.5*(VG12 - VG21)
        OM(1, 3) = 0.5*(VG13 - VG31)
        OM(2, 1) = 0.5*(VG21 - VG12)
        OM(2, 2) = 0.
        OM(2, 3) = 0.5*(VG23 - VG32)
        OM(3, 1) = 0.5*(VG31 - VG13)
        OM(3, 2) = 0.5*(VG32 - VG23)
        OM(3, 3) = 0.

        TMAT(I, J, 1, 1) = SS(1, 1)*SS(1, 1) + SS(1, 2)*SS(2, 1) + SS(1, 3)*SS(3, 1) &
                           + OM(1, 1)*OM(1, 1) + OM(1, 2)*OM(2, 1) + OM(1, 3)*OM(3, 1)
        TMAT(I, J, 1, 2) = SS(1, 1)*SS(1, 2) + SS(1, 2)*SS(2, 2) + SS(1, 3)*SS(3, 2) &
                           + OM(1, 1)*OM(1, 2) + OM(1, 2)*OM(2, 2) + OM(1, 3)*OM(3, 2)
        TMAT(I, J, 1, 3) = SS(1, 1)*SS(1, 3) + SS(1, 2)*SS(2, 3) + SS(1, 3)*SS(3, 3) &
                           + OM(1, 1)*OM(1, 3) + OM(1, 2)*OM(2, 3) + OM(1, 3)*OM(3, 3)
        TMAT(I, J, 2, 2) = SS(2, 1)*SS(1, 2) + SS(2, 2)*SS(2, 2) + SS(2, 3)*SS(3, 2) &
                           + OM(2, 1)*OM(1, 2) + OM(2, 2)*OM(2, 2) + OM(2, 3)*OM(3, 2)
        TMAT(I, J, 2, 3) = SS(2, 1)*SS(1, 3) + SS(2, 2)*SS(2, 3) + SS(2, 3)*SS(3, 3) &
                           + OM(2, 1)*OM(1, 3) + OM(2, 2)*OM(2, 3) + OM(2, 3)*OM(3, 3)
        TMAT(I, J, 3, 3) = SS(3, 1)*SS(1, 3) + SS(3, 2)*SS(2, 3) + SS(3, 3)*SS(3, 3) &
                           + OM(3, 1)*OM(1, 3) + OM(3, 2)*OM(2, 3) + OM(3, 3)*OM(3, 3)
        TMAT(I, J, 2, 1) = TMAT(I, J, 1, 2)
        TMAT(I, J, 3, 1) = TMAT(I, J, 1, 3)
        TMAT(I, J, 3, 2) = TMAT(I, J, 2, 3)
      END DO
    END DO

    DO J = 1, N2M
      DO I = 1, N1M
        DO JJ = 1, 3
          DO II = 1, 3
            TMP(II, JJ) = TMAT(I, J, II, JJ)
          END DO
        END DO
        CALL JACOBI(TMP, 3, 3, DD, EVEC, ITER)
        CALL SORTER(3, DD)
        VLAMBDA2(I, J, K) = 0.5*(DD(2) - ABS(DD(2)))
      END DO
    END DO

  END DO

  IF (IBMON .EQ. 1) THEN
    DO K = 1, N3M
      DO J = 1, N2M
        DO I = 1, N1M

          IF (FUNCBODY(XMP(I), YMP(J), ZMP(K), 0) .LE. 1.E-10) THEN
            UC(I, J, K) = 0.
            VC(I, J, K) = 0.
            WC(I, J, K) = 0.
            VORX(I, J, K) = 0.
            VORY(I, J, K) = 0.
            VORZ(I, J, K) = 0.
            P(I, J, K) = 1000.
            VLAMBDA2(I, J, K) = 10.
          END IF

        END DO
      END DO
    END DO
  END IF

  RETURN
END
!=======================================================================
SUBROUTINE OUTPUTI(IPOINT)
!=======================================================================
  USE MOD_POSTINST
  IMPLICIT NONE
  INTEGER*8    :: IPOINT
  INTEGER*8    :: I, J, K
  INTEGER*8    :: II, JJ, KK
  REAL*8       :: VXX(M2M, M3M), WXX(M2M, M3M), VORMAG

  II = IPOINT

  IF (IOUTFMT .EQ. 1) THEN
    CALL OUTPUTI_TEC(IPOINT)
    RETURN
  END IF

  OPEN (102, FILE='../output/post_inst/'//tname1//ftailijk)
  WRITE (102, '(A)') '# vtk DataFile Version 3.0'
  WRITE (102, '(A)') 'LESwHT instantaneous x-plane'
  WRITE (102, '(A)') 'ASCII'
  WRITE (102, '(A)') 'DATASET STRUCTURED_GRID'
  WRITE (102, *) 'DIMENSIONS ', N3M, N2M, 1
  WRITE (102, *) 'POINTS ', N3M*N2M, ' double'
  DO J = 1, N2M
    DO K = 1, N3M
      WRITE (102, *) XMP(II), YMP(J), ZMP(K)
    END DO
  END DO
  WRITE (102, *) 'POINT_DATA ', N3M*N2M
  WRITE (102, '(A)') 'VECTORS velocity double'
  DO J = 1, N2M
    DO K = 1, N3M
      WRITE (102, *) UC(II, J, K), VC(II, J, K), WC(II, J, K)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS pressure double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO K = 1, N3M
      WRITE (102, *) P(II, J, K)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS vort_x double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO K = 1, N3M
      WRITE (102, *) VORX(II, J, K)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS vort_y double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO K = 1, N3M
      WRITE (102, *) VORY(II, J, K)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS vort_z double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO K = 1, N3M
      WRITE (102, *) VORZ(II, J, K)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS vort_mag double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO K = 1, N3M
      VORMAG = SQRT(VORX(II, J, K)**2 + VORY(II, J, K)**2 + VORZ(II, J, K)**2)
      WRITE (102, *) VORMAG
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS lambda2 double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO K = 1, N3M
      WRITE (102, *) VLAMBDA2(II, J, K)
    END DO
  END DO
  IF (IHTRANS .EQ. 1) THEN
    WRITE (102, '(A)') 'SCALARS temperature double 1'
    WRITE (102, '(A)') 'LOOKUP_TABLE default'
    DO J = 1, N2M
      DO K = 1, N3M
        WRITE (102, *) T(II, J, K)
      END DO
    END DO
  END IF
  CLOSE (102)

  IF (IUNIGRID .EQ. 1) THEN

    DO K = 1, N3M
      DO J = 1, N2M
        JJ = JUNI(J)
        KK = KUNI(K)

        VXX(J, K) = &
          FACUNIY(J)*FACUNIZ(K)*V(II, JJ + 1, KK + 1) &
          + FACUNIY(J)*(1.-FACUNIZ(K))*V(II, JJ + 1, KK) &
          + (1.-FACUNIY(J))*FACUNIZ(K)*V(II, JJ, KK + 1) &
          + (1.-FACUNIY(J))*(1.-FACUNIZ(K))*V(II, JJ, KK)

        WXX(J, K) = &
          FACUNIY(J)*FACUNIZ(K)*W(II, JJ + 1, KK + 1) &
          + FACUNIY(J)*(1.-FACUNIZ(K))*W(II, JJ + 1, KK) &
          + (1.-FACUNIY(J))*FACUNIZ(K)*W(II, JJ, KK + 1) &
          + (1.-FACUNIY(J))*(1.-FACUNIZ(K))*W(II, JJ, KK)
      END DO
    END DO

    OPEN (104, FILE='../output/post_inst/'//tname3//ftailijk)
    WRITE (104, '(A)') '# vtk DataFile Version 3.0'
    WRITE (104, '(A)') 'LESwHT instantaneous x-plane uniform'
    WRITE (104, '(A)') 'ASCII'
    WRITE (104, '(A)') 'DATASET STRUCTURED_GRID'
    WRITE (104, *) 'DIMENSIONS ', KNUMUNI, JNUMUNI, 1
    WRITE (104, *) 'POINTS ', KNUMUNI*JNUMUNI, ' double'
    DO J = 1, N2M, 2
      DO K = 1, N3M, 2
        WRITE (104, *) XMP(II), YUNI(J), ZUNI(K)
      END DO
    END DO
    WRITE (104, *) 'POINT_DATA ', KNUMUNI*JNUMUNI
    WRITE (104, '(A)') 'SCALARS v double 1'
    WRITE (104, '(A)') 'LOOKUP_TABLE default'
    DO J = 1, N2M, 2
      DO K = 1, N3M, 2
        WRITE (104, *) VXX(J, K)
      END DO
    END DO
    WRITE (104, '(A)') 'SCALARS w double 1'
    WRITE (104, '(A)') 'LOOKUP_TABLE default'
    DO J = 1, N2M, 2
      DO K = 1, N3M, 2
        WRITE (104, *) WXX(J, K)
      END DO
    END DO
    CLOSE (104)

  END IF

  IF (IOUTFMT .EQ. 2) CALL OUTPUTI_TEC(IPOINT)

  RETURN
END

!=======================================================================
SUBROUTINE OUTPUTJ(JPOINT)
  !=======================================================================
  USE MOD_POSTINST
  IMPLICIT NONE
  INTEGER*8    :: JPOINT
  INTEGER*8    :: I, J, K
  INTEGER*8    :: II, JJ, KK, IIP, KKP
  REAL*8       :: UXX(M1M, M3M), WXX(M1M, M3M), VORMAG

  JJ = JPOINT

  IF (IOUTFMT .EQ. 1) THEN
    CALL OUTPUTJ_TEC(JPOINT)
    RETURN
  END IF

  OPEN (102, FILE='../output/post_inst/'//tname1//ftailijk)
  WRITE (102, '(A)') '# vtk DataFile Version 3.0'
  WRITE (102, '(A)') 'LESwHT instantaneous y-plane'
  WRITE (102, '(A)') 'ASCII'
  WRITE (102, '(A)') 'DATASET STRUCTURED_GRID'
  WRITE (102, *) 'DIMENSIONS ', N1M, N3M, 1
  WRITE (102, *) 'POINTS ', N1M*N3M, ' double'
  DO K = 1, N3M
    DO I = 1, N1M
      WRITE (102, *) XMP(I), YMP(JJ), ZMP(K)
    END DO
  END DO
  WRITE (102, *) 'POINT_DATA ', N1M*N3M
  WRITE (102, '(A)') 'VECTORS velocity double'
  DO K = 1, N3M
    DO I = 1, N1M
      WRITE (102, *) UC(I, JJ, K), VC(I, JJ, K), WC(I, JJ, K)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS pressure double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO K = 1, N3M
    DO I = 1, N1M
      WRITE (102, *) P(I, JJ, K)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS vort_x double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO K = 1, N3M
    DO I = 1, N1M
      WRITE (102, *) VORX(I, JJ, K)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS vort_y double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO K = 1, N3M
    DO I = 1, N1M
      WRITE (102, *) VORY(I, JJ, K)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS vort_z double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO K = 1, N3M
    DO I = 1, N1M
      WRITE (102, *) VORZ(I, JJ, K)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS vort_mag double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO K = 1, N3M
    DO I = 1, N1M
      VORMAG = SQRT(VORX(I, JJ, K)**2 + VORY(I, JJ, K)**2 + VORZ(I, JJ, K)**2)
      WRITE (102, *) VORMAG
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS lambda2 double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO K = 1, N3M
    DO I = 1, N1M
      WRITE (102, *) VLAMBDA2(I, JJ, K)
    END DO
  END DO
  IF (IHTRANS .EQ. 1) THEN
    WRITE (102, '(A)') 'SCALARS temperature double 1'
    WRITE (102, '(A)') 'LOOKUP_TABLE default'
    DO K = 1, N3M
      DO I = 1, N1M
        WRITE (102, *) T(I, JJ, K)
      END DO
    END DO
  END IF
  CLOSE (102)

  IF (IUNIGRID .EQ. 1) THEN
    DO K = 1, N3M
      DO I = 1, N1M

        II = IUNI(I)
        KK = KUNI(K)
        IIP = II + 1
        KKP = KK + 1
        UXX(I, K) = &
          FACUNIX(I)*FACUNIZ(K)*U(IIP, JJ, KKP) &
          + FACUNIX(I)*(1.-FACUNIZ(K))*U(IIP, JJ, KK) &
          + (1.-FACUNIX(I))*FACUNIZ(K)*U(II, JJ, KKP) &
          + (1.-FACUNIX(I))*(1.-FACUNIZ(K))*U(II, JJ, KK)

        WXX(I, K) = &
          FACUNIX(I)*FACUNIZ(K)*W(IIP, JJ, KKP) &
          + FACUNIX(I)*(1.-FACUNIZ(K))*W(IIP, JJ, KK) &
          + (1.-FACUNIX(I))*FACUNIZ(K)*W(II, JJ, KKP) &
          + (1.-FACUNIX(I))*(1.-FACUNIZ(K))*W(II, JJ, KK)
      END DO
    END DO

    OPEN (104, FILE='../output/post_inst/'//tname3//ftailijk)
    WRITE (104, '(A)') '# vtk DataFile Version 3.0'
    WRITE (104, '(A)') 'LESwHT instantaneous y-plane uniform'
    WRITE (104, '(A)') 'ASCII'
    WRITE (104, '(A)') 'DATASET STRUCTURED_GRID'
    WRITE (104, *) 'DIMENSIONS ', INUMUNI, KNUMUNI, 1
    WRITE (104, *) 'POINTS ', INUMUNI*KNUMUNI, ' double'
    DO K = 1, N3M, 2
      DO I = 1, N1M, 2
        WRITE (104, *) XUNI(I), YMP(JJ), ZUNI(K)
      END DO
    END DO
    WRITE (104, *) 'POINT_DATA ', INUMUNI*KNUMUNI
    WRITE (104, '(A)') 'SCALARS u double 1'
    WRITE (104, '(A)') 'LOOKUP_TABLE default'
    DO K = 1, N3M, 2
      DO I = 1, N1M, 2
        WRITE (104, *) UXX(I, K)
      END DO
    END DO
    WRITE (104, '(A)') 'SCALARS w double 1'
    WRITE (104, '(A)') 'LOOKUP_TABLE default'
    DO K = 1, N3M, 2
      DO I = 1, N1M, 2
        WRITE (104, *) WXX(I, K)
      END DO
    END DO
    CLOSE (104)
  END IF

  IF (IOUTFMT .EQ. 2) CALL OUTPUTJ_TEC(JPOINT)

  RETURN
END

!=======================================================================
SUBROUTINE OUTPUTK(KPOINT)
  !=======================================================================
  USE MOD_POSTINST
  IMPLICIT NONE
  INTEGER*8    :: KPOINT
  INTEGER*8    :: I, J, K
  INTEGER*8    :: II, JJ, KK, IIP, JJP
  REAL*8       :: UXX(M1M, M2M), VXX(M1M, M2M), VORMAG

  KK = KPOINT

  IF (IOUTFMT .EQ. 1) THEN
    CALL OUTPUTK_TEC(KPOINT)
    RETURN
  END IF

  OPEN (102, FILE='../output/post_inst/'//tname1//ftailijk)
  WRITE (102, '(A)') '# vtk DataFile Version 3.0'
  WRITE (102, '(A)') 'LESwHT instantaneous z-plane'
  WRITE (102, '(A)') 'ASCII'
  WRITE (102, '(A)') 'DATASET STRUCTURED_GRID'
  WRITE (102, *) 'DIMENSIONS ', N1M, N2M, 1
  WRITE (102, *) 'POINTS ', N1M*N2M, ' double'
  DO J = 1, N2M
    DO I = 1, N1M
      WRITE (102, *) XMP(I), YMP(J), ZMP(KK)
    END DO
  END DO
  WRITE (102, *) 'POINT_DATA ', N1M*N2M
  WRITE (102, '(A)') 'VECTORS velocity double'
  DO J = 1, N2M
    DO I = 1, N1M
      WRITE (102, *) UC(I, J, KK), VC(I, J, KK), WC(I, J, KK)
    END DO
  END DO
  IF (IHTRANS .EQ. 1) THEN
    WRITE (102, '(A)') 'SCALARS temperature double 1'
    WRITE (102, '(A)') 'LOOKUP_TABLE default'
    DO J = 1, N2M
      DO I = 1, N1M
        WRITE (102, *) T(I, J, KK)
      END DO
    END DO
  END IF
  WRITE (102, '(A)') 'SCALARS pressure double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO I = 1, N1M
      WRITE (102, *) P(I, J, KK)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS vort_x double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO I = 1, N1M
      WRITE (102, *) VORX(I, J, KK)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS vort_y double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO I = 1, N1M
      WRITE (102, *) VORY(I, J, KK)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS vort_z double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO I = 1, N1M
      WRITE (102, *) VORZ(I, J, KK)
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS vort_mag double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO I = 1, N1M
      VORMAG = SQRT(VORX(I, J, KK)**2 + VORY(I, J, KK)**2 + VORZ(I, J, KK)**2)
      WRITE (102, *) VORMAG
    END DO
  END DO
  WRITE (102, '(A)') 'SCALARS lambda2 double 1'
  WRITE (102, '(A)') 'LOOKUP_TABLE default'
  DO J = 1, N2M
    DO I = 1, N1M
      WRITE (102, *) VLAMBDA2(I, J, KK)
    END DO
  END DO
  CLOSE (102)

  IF (IUNIGRID .EQ. 1) THEN

    DO J = 1, N2M
      DO I = 1, N1M

        II = IUNI(I)
        JJ = JUNI(J)
        IIP = II + 1
        JJP = JJ + 1

        UXX(I, J) = &
          FACUNIX(I)*FACUNIY(J)*U(IIP, JJP, KK) &
          + FACUNIX(I)*(1.-FACUNIY(J))*U(IIP, JJ, KK) &
          + (1.-FACUNIX(I))*FACUNIY(J)*U(II, JJP, KK) &
          + (1.-FACUNIX(I))*(1.-FACUNIY(J))*U(II, JJ, KK)

        VXX(I, J) = &
          FACUNIX(I)*FACUNIY(J)*V(IIP, JJP, KK) &
          + FACUNIX(I)*(1.-FACUNIY(J))*V(IIP, JJ, KK) &
          + (1.-FACUNIX(I))*FACUNIY(J)*V(II, JJP, KK) &
          + (1.-FACUNIX(I))*(1.-FACUNIY(J))*V(II, JJ, KK)
      END DO
    END DO

    OPEN (104, FILE='../output/post_inst/'//tname3//ftailijk)
    WRITE (104, '(A)') '# vtk DataFile Version 3.0'
    WRITE (104, '(A)') 'LESwHT instantaneous z-plane uniform'
    WRITE (104, '(A)') 'ASCII'
    WRITE (104, '(A)') 'DATASET STRUCTURED_GRID'
    WRITE (104, *) 'DIMENSIONS ', INUMUNI, JNUMUNI, 1
    WRITE (104, *) 'POINTS ', INUMUNI*JNUMUNI, ' double'
    DO J = 1, N2M, 2
      DO I = 1, N1M, 2
        WRITE (104, *) XUNI(I), YUNI(J), ZMP(KK)
      END DO
    END DO
    WRITE (104, *) 'POINT_DATA ', INUMUNI*JNUMUNI
    WRITE (104, '(A)') 'SCALARS u double 1'
    WRITE (104, '(A)') 'LOOKUP_TABLE default'
    DO J = 1, N2M, 2
      DO I = 1, N1M, 2
        WRITE (104, *) UXX(I, J)
      END DO
    END DO
    WRITE (104, '(A)') 'SCALARS v double 1'
    WRITE (104, '(A)') 'LOOKUP_TABLE default'
    DO J = 1, N2M, 2
      DO I = 1, N1M, 2
        WRITE (104, *) VXX(I, J)
      END DO
    END DO
    CLOSE (104)
  END IF

  IF (IOUTFMT .EQ. 2) CALL OUTPUTK_TEC(KPOINT)

  RETURN
END

!=======================================================================
SUBROUTINE OUTPUT_3D
  !=======================================================================
  USE MOD_POSTINST
  IMPLICIT NONE
  INTEGER*8    :: I, J, K
  INTEGER*8    :: INUM, JNUM, KNUM
  REAL*8       :: VORMAG

  INUM = 0
  JNUM = 0
  KNUM = 0
  DO I = ISTART, IEND, ISKIP
    INUM = INUM + 1
  END DO
  DO J = JSTART, JEND, JSKIP
    JNUM = JNUM + 1
  END DO
  DO K = KSTART, KEND, KSKIP
    KNUM = KNUM + 1
  END DO

  IF (IOUTFMT .EQ. 1) THEN
    CALL OUTPUT_3D_TEC(INUM, JNUM, KNUM)
    RETURN
  END IF

  IF (IOUTFMT .EQ. 2) CALL OUTPUT_3D_TEC(INUM, JNUM, KNUM)

  IF (ILD2 .EQ. 1) THEN
    OPEN (101, FILE='../output/post_inst/'//tname2//'_ld2.vtk')
    WRITE (101, '(A)') '# vtk DataFile Version 3.0'
    WRITE (101, '(A)') 'LESwHT instantaneous 3D lambda2'
    WRITE (101, '(A)') 'ASCII'
    WRITE (101, '(A)') 'DATASET STRUCTURED_GRID'
    WRITE (101, *) 'DIMENSIONS ', INUM, JNUM, KNUM
    WRITE (101, *) 'POINTS ', INUM*JNUM*KNUM, ' double'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (101, *) XMP(I), YMP(J), ZMP(K)
        END DO
      END DO
    END DO
    WRITE (101, *) 'POINT_DATA ', INUM*JNUM*KNUM
    WRITE (101, '(A)') 'SCALARS lambda2 double 1'
    WRITE (101, '(A)') 'LOOKUP_TABLE default'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (101, *) VLAMBDA2(I, J, K)
        END DO
      END DO
    END DO
    CLOSE (101)
  END IF

  IF (IUVWP .EQ. 1) THEN
    OPEN (102, FILE='../output/post_inst/'//tname2//'_uvwp.vtk')
    WRITE (102, '(A)') '# vtk DataFile Version 3.0'
    WRITE (102, '(A)') 'LESwHT instantaneous 3D uvwp'
    WRITE (102, '(A)') 'ASCII'
    WRITE (102, '(A)') 'DATASET STRUCTURED_GRID'
    WRITE (102, *) 'DIMENSIONS ', INUM, JNUM, KNUM
    WRITE (102, *) 'POINTS ', INUM*JNUM*KNUM, ' double'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (102, *) XMP(I), YMP(J), ZMP(K)
        END DO
      END DO
    END DO
    WRITE (102, *) 'POINT_DATA ', INUM*JNUM*KNUM
    WRITE (102, '(A)') 'VECTORS velocity double'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (102, *) UC(I, J, K), VC(I, J, K), WC(I, J, K)
        END DO
      END DO
    END DO
    WRITE (102, '(A)') 'SCALARS pressure double 1'
    WRITE (102, '(A)') 'LOOKUP_TABLE default'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (102, *) P(I, J, K)
        END DO
      END DO
    END DO
    IF (IHTRANS .EQ. 1) THEN
      WRITE (102, '(A)') 'SCALARS temperature double 1'
      WRITE (102, '(A)') 'LOOKUP_TABLE default'
      DO K = KSTART, KEND, KSKIP
        DO J = JSTART, JEND, JSKIP
          DO I = ISTART, IEND, ISKIP
            WRITE (102, *) T(I, J, K)
          END DO
        END DO
      END DO
    END IF
    CLOSE (102)
  END IF

  IF (IWXYZ .EQ. 1) THEN
    OPEN (103, FILE='../output/post_inst/'//tname2//'_vorxyz.vtk')
    WRITE (103, '(A)') '# vtk DataFile Version 3.0'
    WRITE (103, '(A)') 'LESwHT instantaneous 3D vorticity'
    WRITE (103, '(A)') 'ASCII'
    WRITE (103, '(A)') 'DATASET STRUCTURED_GRID'
    WRITE (103, *) 'DIMENSIONS ', INUM, JNUM, KNUM
    WRITE (103, *) 'POINTS ', INUM*JNUM*KNUM, ' double'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (103, *) XMP(I), YMP(J), ZMP(K)
        END DO
      END DO
    END DO
    WRITE (103, *) 'POINT_DATA ', INUM*JNUM*KNUM
    WRITE (103, '(A)') 'VECTORS vorticity double'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (103, *) VORX(I, J, K), VORY(I, J, K), VORZ(I, J, K)
        END DO
      END DO
    END DO
    WRITE (103, '(A)') 'SCALARS vort_mag double 1'
    WRITE (103, '(A)') 'LOOKUP_TABLE default'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          VORMAG = SQRT(VORX(I, J, K)**2 + VORY(I, J, K)**2 + VORZ(I, J, K)**2)
          WRITE (103, *) VORMAG
        END DO
      END DO
    END DO
    CLOSE (103)
  END IF

  IF (IALL .EQ. 1) THEN
    OPEN (103, FILE='../output/post_inst/'//tname2//'_all.vtk')
    WRITE (103, '(A)') '# vtk DataFile Version 3.0'
    WRITE (103, '(A)') 'LESwHT instantaneous 3D all'
    WRITE (103, '(A)') 'ASCII'
    WRITE (103, '(A)') 'DATASET STRUCTURED_GRID'
    WRITE (103, *) 'DIMENSIONS ', INUM, JNUM, KNUM
    WRITE (103, *) 'POINTS ', INUM*JNUM*KNUM, ' double'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (103, *) XMP(I), YMP(J), ZMP(K)
        END DO
      END DO
    END DO
    WRITE (103, *) 'POINT_DATA ', INUM*JNUM*KNUM
    WRITE (103, '(A)') 'VECTORS velocity double'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (103, *) UC(I, J, K), VC(I, J, K), WC(I, J, K)
        END DO
      END DO
    END DO
    WRITE (103, '(A)') 'SCALARS pressure double 1'
    WRITE (103, '(A)') 'LOOKUP_TABLE default'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (103, *) P(I, J, K)
        END DO
      END DO
    END DO
    WRITE (103, '(A)') 'VECTORS vorticity double'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (103, *) VORX(I, J, K), VORY(I, J, K), VORZ(I, J, K)
        END DO
      END DO
    END DO
    WRITE (103, '(A)') 'SCALARS vort_mag double 1'
    WRITE (103, '(A)') 'LOOKUP_TABLE default'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          VORMAG = SQRT(VORX(I, J, K)**2 + VORY(I, J, K)**2 + VORZ(I, J, K)**2)
          WRITE (103, *) VORMAG
        END DO
      END DO
    END DO
    WRITE (103, '(A)') 'SCALARS lambda2 double 1'
    WRITE (103, '(A)') 'LOOKUP_TABLE default'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (103, *) VLAMBDA2(I, J, K)
        END DO
      END DO
    END DO
    IF (IHTRANS .EQ. 1) THEN
      WRITE (103, '(A)') 'SCALARS temperature double 1'
      WRITE (103, '(A)') 'LOOKUP_TABLE default'
      DO K = KSTART, KEND, KSKIP
        DO J = JSTART, JEND, JSKIP
          DO I = ISTART, IEND, ISKIP
            WRITE (103, *) T(I, J, K)
          END DO
        END DO
      END DO
    END IF

    CLOSE (103)
  END IF

  RETURN
END

!=======================================================================
SUBROUTINE OUTPUTI_TEC(IPOINT)
!=======================================================================
  USE MOD_POSTINST
  IMPLICIT NONE
  INTEGER*8 :: IPOINT, I, J, K
  REAL*8 :: VORMAG
  CHARACTER*9 :: FTAILTEC

  FTAILTEC = FTAILIJK(1:5)//'.tec'
  OPEN (202, FILE='../output/post_inst/'//TNAME1//FTAILTEC)
  IF (IHTRANS .EQ. 1) THEN
    WRITE (202, '(A)') 'VARIABLES = "X","Y","Z","U","V","W","P","VORX","VORY","VORZ","VORMAG","LAMBDA2","T"'
  ELSE
    WRITE (202, '(A)') 'VARIABLES = "X","Y","Z","U","V","W","P","VORX","VORY","VORZ","VORMAG","LAMBDA2"'
  END IF
  WRITE (202, '(A,I0,A,I0,A)') 'ZONE I=', N3M, ', J=', N2M, ', F=POINT'
  DO J = 1, N2M
    DO K = 1, N3M
      VORMAG = SQRT(VORX(IPOINT, J, K)**2 + VORY(IPOINT, J, K)**2 + VORZ(IPOINT, J, K)**2)
      IF (IHTRANS .EQ. 1) THEN
        WRITE (202, *) XMP(IPOINT), YMP(J), ZMP(K), UC(IPOINT, J, K), VC(IPOINT, J, K), WC(IPOINT, J, K), P(IPOINT, J, K), &
          VORX(IPOINT, J, K), VORY(IPOINT, J, K), VORZ(IPOINT, J, K), VORMAG, VLAMBDA2(IPOINT, J, K), T(IPOINT, J, K)
      ELSE
        WRITE (202, *) XMP(IPOINT), YMP(J), ZMP(K), UC(IPOINT, J, K), VC(IPOINT, J, K), WC(IPOINT, J, K), P(IPOINT, J, K), &
          VORX(IPOINT, J, K), VORY(IPOINT, J, K), VORZ(IPOINT, J, K), VORMAG, VLAMBDA2(IPOINT, J, K)
      END IF
    END DO
  END DO
  CLOSE (202)
  RETURN
END

!=======================================================================
SUBROUTINE OUTPUTJ_TEC(JPOINT)
!=======================================================================
  USE MOD_POSTINST
  IMPLICIT NONE
  INTEGER*8 :: JPOINT, I, K
  REAL*8 :: VORMAG
  CHARACTER*9 :: FTAILTEC

  FTAILTEC = FTAILIJK(1:5)//'.tec'
  OPEN (202, FILE='../output/post_inst/'//TNAME1//FTAILTEC)
  IF (IHTRANS .EQ. 1) THEN
    WRITE (202, '(A)') 'VARIABLES = "X","Y","Z","U","V","W","P","VORX","VORY","VORZ","VORMAG","LAMBDA2","T"'
  ELSE
    WRITE (202, '(A)') 'VARIABLES = "X","Y","Z","U","V","W","P","VORX","VORY","VORZ","VORMAG","LAMBDA2"'
  END IF
  WRITE (202, '(A,I0,A,I0,A)') 'ZONE I=', N1M, ', J=', N3M, ', F=POINT'
  DO K = 1, N3M
    DO I = 1, N1M
      VORMAG = SQRT(VORX(I, JPOINT, K)**2 + VORY(I, JPOINT, K)**2 + VORZ(I, JPOINT, K)**2)
      IF (IHTRANS .EQ. 1) THEN
        WRITE (202, *) XMP(I), YMP(JPOINT), ZMP(K), UC(I, JPOINT, K), VC(I, JPOINT, K), WC(I, JPOINT, K), P(I, JPOINT, K), &
          VORX(I, JPOINT, K), VORY(I, JPOINT, K), VORZ(I, JPOINT, K), VORMAG, VLAMBDA2(I, JPOINT, K), T(I, JPOINT, K)
      ELSE
        WRITE (202, *) XMP(I), YMP(JPOINT), ZMP(K), UC(I, JPOINT, K), VC(I, JPOINT, K), WC(I, JPOINT, K), P(I, JPOINT, K), &
          VORX(I, JPOINT, K), VORY(I, JPOINT, K), VORZ(I, JPOINT, K), VORMAG, VLAMBDA2(I, JPOINT, K)
      END IF
    END DO
  END DO
  CLOSE (202)
  RETURN
END

!=======================================================================
SUBROUTINE OUTPUTK_TEC(KPOINT)
!=======================================================================
  USE MOD_POSTINST
  IMPLICIT NONE
  INTEGER*8 :: KPOINT, I, J
  REAL*8 :: VORMAG
  CHARACTER*9 :: FTAILTEC

  FTAILTEC = FTAILIJK(1:5)//'.tec'
  OPEN (202, FILE='../output/post_inst/'//TNAME1//FTAILTEC)
  IF (IHTRANS .EQ. 1) THEN
    WRITE (202, '(A)') 'VARIABLES = "X","Y","Z","U","V","W","P","VORX","VORY","VORZ","VORMAG","LAMBDA2","T"'
  ELSE
    WRITE (202, '(A)') 'VARIABLES = "X","Y","Z","U","V","W","P","VORX","VORY","VORZ","VORMAG","LAMBDA2"'
  END IF
  WRITE (202, '(A,I0,A,I0,A)') 'ZONE I=', N1M, ', J=', N2M, ', F=POINT'
  DO J = 1, N2M
    DO I = 1, N1M
      VORMAG = SQRT(VORX(I, J, KPOINT)**2 + VORY(I, J, KPOINT)**2 + VORZ(I, J, KPOINT)**2)
      IF (IHTRANS .EQ. 1) THEN
        WRITE (202, *) XMP(I), YMP(J), ZMP(KPOINT), UC(I, J, KPOINT), VC(I, J, KPOINT), WC(I, J, KPOINT), P(I, J, KPOINT), &
          VORX(I, J, KPOINT), VORY(I, J, KPOINT), VORZ(I, J, KPOINT), VORMAG, VLAMBDA2(I, J, KPOINT), T(I, J, KPOINT)
      ELSE
        WRITE (202, *) XMP(I), YMP(J), ZMP(KPOINT), UC(I, J, KPOINT), VC(I, J, KPOINT), WC(I, J, KPOINT), P(I, J, KPOINT), &
          VORX(I, J, KPOINT), VORY(I, J, KPOINT), VORZ(I, J, KPOINT), VORMAG, VLAMBDA2(I, J, KPOINT)
      END IF
    END DO
  END DO
  CLOSE (202)
  RETURN
END

!=======================================================================
SUBROUTINE OUTPUT_3D_TEC(INUM, JNUM, KNUM)
!=======================================================================
  USE MOD_POSTINST
  IMPLICIT NONE
  INTEGER*8 :: INUM, JNUM, KNUM, I, J, K
  REAL*8 :: VORMAG

  IF (ILD2 .EQ. 1) THEN
    OPEN (210, FILE='../output/post_inst/'//TNAME2//'_ld2.tec')
    WRITE (210, '(A)') 'VARIABLES = "X","Y","Z","LAMBDA2"'
    WRITE (210, '(A,I0,A,I0,A,I0,A)') 'ZONE I=', INUM, ', J=', JNUM, ', K=', KNUM, ', F=POINT'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          WRITE (210, *) XMP(I), YMP(J), ZMP(K), VLAMBDA2(I, J, K)
        END DO
      END DO
    END DO
    CLOSE (210)
  END IF

  IF (IUVWP .EQ. 1) THEN
    OPEN (211, FILE='../output/post_inst/'//TNAME2//'_uvwp.tec')
    IF (IHTRANS .EQ. 1) THEN
      WRITE (211, '(A)') 'VARIABLES = "X","Y","Z","U","V","W","P","T"'
    ELSE
      WRITE (211, '(A)') 'VARIABLES = "X","Y","Z","U","V","W","P"'
    END IF
    WRITE (211, '(A,I0,A,I0,A,I0,A)') 'ZONE I=', INUM, ', J=', JNUM, ', K=', KNUM, ', F=POINT'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          IF (IHTRANS .EQ. 1) THEN
            WRITE (211, *) XMP(I), YMP(J), ZMP(K), UC(I, J, K), VC(I, J, K), WC(I, J, K), P(I, J, K), T(I, J, K)
          ELSE
            WRITE (211, *) XMP(I), YMP(J), ZMP(K), UC(I, J, K), VC(I, J, K), WC(I, J, K), P(I, J, K)
          END IF
        END DO
      END DO
    END DO
    CLOSE (211)
  END IF

  IF (IWXYZ .EQ. 1) THEN
    OPEN (212, FILE='../output/post_inst/'//TNAME2//'_vorxyz.tec')
    WRITE (212, '(A)') 'VARIABLES = "X","Y","Z","VORX","VORY","VORZ","VORMAG"'
    WRITE (212, '(A,I0,A,I0,A,I0,A)') 'ZONE I=', INUM, ', J=', JNUM, ', K=', KNUM, ', F=POINT'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          VORMAG = SQRT(VORX(I, J, K)**2 + VORY(I, J, K)**2 + VORZ(I, J, K)**2)
          WRITE (212, *) XMP(I), YMP(J), ZMP(K), VORX(I, J, K), VORY(I, J, K), VORZ(I, J, K), VORMAG
        END DO
      END DO
    END DO
    CLOSE (212)
  END IF

  IF (IALL .EQ. 1) THEN
    OPEN (213, FILE='../output/post_inst/'//TNAME2//'_all.tec')
    IF (IHTRANS .EQ. 1) THEN
      WRITE (213, '(A)') 'VARIABLES = "X","Y","Z","U","V","W","P","VORX","VORY","VORZ","VORMAG","LAMBDA2","T"'
    ELSE
      WRITE (213, '(A)') 'VARIABLES = "X","Y","Z","U","V","W","P","VORX","VORY","VORZ","VORMAG","LAMBDA2"'
    END IF
    WRITE (213, '(A,I0,A,I0,A,I0,A)') 'ZONE I=', INUM, ', J=', JNUM, ', K=', KNUM, ', F=POINT'
    DO K = KSTART, KEND, KSKIP
      DO J = JSTART, JEND, JSKIP
        DO I = ISTART, IEND, ISKIP
          VORMAG = SQRT(VORX(I, J, K)**2 + VORY(I, J, K)**2 + VORZ(I, J, K)**2)
          IF (IHTRANS .EQ. 1) THEN
            WRITE (213, *) XMP(I), YMP(J), ZMP(K), UC(I, J, K), VC(I, J, K), WC(I, J, K), P(I, J, K), &
              VORX(I, J, K), VORY(I, J, K), VORZ(I, J, K), VORMAG, VLAMBDA2(I, J, K), T(I, J, K)
          ELSE
            WRITE (213, *) XMP(I), YMP(J), ZMP(K), UC(I, J, K), VC(I, J, K), WC(I, J, K), P(I, J, K), &
              VORX(I, J, K), VORY(I, J, K), VORZ(I, J, K), VORMAG, VLAMBDA2(I, J, K)
          END IF
        END DO
      END DO
    END DO
    CLOSE (213)
  END IF

  RETURN
END

!=======================================================================
SUBROUTINE PREFLD(fileprevel)
  !=======================================================================
  USE MOD_POSTINST
  IMPLICIT NONE
  CHARACTER*100 :: fileprevel
  INTEGER*8    :: I, J, K
  INTEGER*8    :: NN1, NN2, NN3, IHIST, M, IPSS, IXPRDIC, IZPRDIC, IDUM
  REAL*8       :: RRE, TIME, DT, PRA, DUM, GRA
  INTEGER*8    :: IOS
  CHARACTER*256:: FPATH

  U = 0.
  V = 0.
  W = 0.
  P = 0.
  T = 0.

  FPATH = '../output/field/'//TRIM(fileprevel)
  OPEN (12, FILE=TRIM(FPATH), FORM='UNFORMATTED', STATUS='OLD', ACTION='READ', IOSTAT=IOS)
  IF (IOS .NE. 0) THEN
    WRITE (*, *) 'ERROR: CANNOT OPEN INSTANTANEOUS FIELD FILE: ', TRIM(FPATH)
    WRITE (*, *) 'IOSTAT=', IOS
    STOP
  END IF

  !     dum for future use
  READ (12, IOSTAT=IOS) NN1, NN2, NN3, RRE, PRA, GRA
  IF (IOS .NE. 0) GOTO 901
  READ (12, IOSTAT=IOS) IHIST, M, TIME, DT
  IF (IOS .NE. 0) GOTO 902
  READ (12, IOSTAT=IOS) XPRDIC, YPRDIC, ZPRDIC
  IF (IOS .NE. 0) GOTO 903
  READ (12, IOSTAT=IOS) BC_XBTM, BC_XTOP, BC_YBTM, BC_YTOP, BC_ZBTM, BC_ZTOP
  IF (IOS .NE. 0) GOTO 904
  READ (12, IOSTAT=IOS) ICH, ITEMP
  IF (IOS .NE. 0) GOTO 905

  IF ((NN1 .NE. N1) .OR. (NN2 .NE. N2) .OR. (NN3 .NE. N3)) THEN
    WRITE (*, *) 'ERROR: FIELD/GRID SIZE MISMATCH.'
    WRITE (*, *) 'GRID N1,N2,N3 = ', N1, N2, N3
    WRITE (*, *) 'FIELD N1,N2,N3 = ', NN1, NN2, NN3
    WRITE (*, *) 'FILE = ', TRIM(FPATH)
    STOP
  END IF

  RE = RRE
  PR = PRA
  GR = GRA

  READ (12, IOSTAT=IOS) (((U(I, J, K), I=1, N1), J=0, N2), K=0, N3)
  IF (IOS .NE. 0) GOTO 906
  READ (12, IOSTAT=IOS) (((V(I, J, K), I=0, N1), J=1, N2), K=0, N3)
  IF (IOS .NE. 0) GOTO 907
  READ (12, IOSTAT=IOS) (((W(I, J, K), I=0, N1), J=0, N2), K=1, N3)
  IF (IOS .NE. 0) GOTO 908
  READ (12, IOSTAT=IOS) (((P(I, J, K), I=1, N1M), J=1, N2M), K=1, N3M)
  IF (IOS .NE. 0) GOTO 909
  IF (IHTRANS .EQ. 1) THEN
    READ (12, IOSTAT=IOS) (((T(I, J, K), I=1, N1M), J=1, N2M), K=1, N3M)
    IF (IOS .NE. 0) GOTO 910
  END IF

  CLOSE (12)

  WRITE (*, 100)
  WRITE (*, 101)
  WRITE (*, 102) fileprevel
  WRITE (*, 103) RE
  WRITE (*, 104) NN1, NN2, NN3
  WRITE (*, 105) IHIST, M, TIME, DT

100 FORMAT('----------- INITIAL FIELD INFORMATION -----------')
101 FORMAT('INITIAL FIELD      : READING DONE')
102 FORMAT('INITIAL FIELD NAME : ', A30)
103 FORMAT('RE=', ES12.3)
104 FORMAT('N1=', I10, '  N2=', I10, '  N3=', I10)
105 FORMAT('IHIST=', I8, '  M=', I10, '  TIME= ', F12.5, '  DT=', F12.5)

  !     Z PERIODICITY
  IF (ZPRDIC .EQ. 1) THEN
    DO J = 0, N2
      DO I = 1, N1
        U(I, J, 0) = U(I, J, N3M)
        U(I, J, N3) = U(I, J, 1)
      END DO
    END DO

    DO J = 1, N2
      DO I = 0, N1
        V(I, J, 0) = V(I, J, N3M)
        V(I, J, N3) = V(I, J, 1)
      END DO
    END DO

    DO J = 0, N2
      DO I = 0, N1
        W(I, J, 0) = W(I, J, N3M)
        W(I, J, N3) = W(I, J, 1)
      END DO
    END DO
  END IF

  !     X PERIODICITY
  IF (XPRDIC .EQ. 1) THEN
    DO K = 0, N3
      DO J = 0, N2
        U(0, J, K) = U(N1M, J, K)
        U(N1, J, K) = U(1, J, K)
      END DO
    END DO

    DO K = 0, N3
      DO J = 1, N2
        V(0, J, K) = V(N1M, J, K)
        V(N1, J, K) = V(1, J, K)
      END DO
    END DO

    DO K = 1, N3
      DO J = 0, N2
        W(0, J, K) = W(N1M, J, K)
        W(N1, J, K) = W(1, J, K)
      END DO
    END DO
  END IF

  !     Z PERIODICITY
  IF (ZPRDIC .EQ. 1) THEN
    DO J = 1, N2M
      DO I = 1, N1M
        P(I, J, 0) = P(I, J, N3M)
        P(I, J, N3) = P(I, J, 1)
        IF (IHTRANS .EQ. 1) THEN
          T(I, J, 0) = T(I, J, N3M)
          T(I, J, N3) = T(I, J, 1)
        END IF
      END DO
    END DO
  END IF

  !     X PERIODICITY
  IF (XPRDIC .EQ. 1) THEN
    DO K = 1, N3M
      DO J = 1, N2M
        P(0, J, K) = P(N1M, J, K)
        P(N1, J, K) = P(1, J, K)
        IF (IHTRANS .EQ. 1) THEN
          T(0, J, K) = T(N1M, J, K)
          T(N1, J, K) = T(1, J, K)
        END IF
      END DO
    END DO
  END IF

  RETURN
901 WRITE (*, *) 'ERROR READING HEADER-1 FROM ', TRIM(FPATH), ' IOSTAT=', IOS
  STOP
902 WRITE (*, *) 'ERROR READING HEADER-2 FROM ', TRIM(FPATH), ' IOSTAT=', IOS
  STOP
903 WRITE (*, *) 'ERROR READING PERIODIC FLAGS FROM ', TRIM(FPATH), ' IOSTAT=', IOS
  STOP
904 WRITE (*, *) 'ERROR READING BC FLAGS FROM ', TRIM(FPATH), ' IOSTAT=', IOS
  STOP
905 WRITE (*, *) 'ERROR READING ICH/ITEMP FROM ', TRIM(FPATH), ' IOSTAT=', IOS
  STOP
906 WRITE (*, *) 'ERROR READING U FIELD FROM ', TRIM(FPATH), ' IOSTAT=', IOS
  STOP
907 WRITE (*, *) 'ERROR READING V FIELD FROM ', TRIM(FPATH), ' IOSTAT=', IOS
  STOP
908 WRITE (*, *) 'ERROR READING W FIELD FROM ', TRIM(FPATH), ' IOSTAT=', IOS
  STOP
909 WRITE (*, *) 'ERROR READING P FIELD FROM ', TRIM(FPATH), ' IOSTAT=', IOS
  STOP
910 WRITE (*, *) 'ERROR READING T FIELD FROM ', TRIM(FPATH), ' IOSTAT=', IOS
  STOP
END

!******************************************************************
!       SECTION FOR GEOMETRY
!      - READ GRID FILE & COMPUTE VARIABLES ABOUT GEOMETRY
!      - GEOMETRY VARIABLES ARE LISTED IN 'GEOM.H'
!      - DATA READ FROM FILE ARE BLOCK 'DIM','COORD','SCALES',
!        'GEOMINPUT'
!      - COMPUTED VARIABLES ARE BLOCK 'GEOMETRY','INDEX','FIX',
!        'POSITION','VAR','VVAR'
!
!                               Revised by DONGJOO KIM
!                               Turbulence & Flow Control Lab.
!                               School of Mech. & Aerospace Eng.
!                               Seoul National Univ.
!                               9/27/1999
!******************************************************************
SUBROUTINE GEOM(gridfile)
  USE MOD_POSTINST
  IMPLICIT NONE
  INTEGER*8    :: I, K
  CHARACTER*100:: gridfile

  OPEN (11, FILE=gridfile)
  CALL READGRID
  CLOSE (11)

  CALL INDICES
  CALL INDXFIX
  CALL MESHES
  CALL PHYSPOS

  IF (ZPRDIC .EQ. 1) THEN
    FIXKL(1) = 0.
    FIXKU(N3M) = 0.
    SDZ(0) = SDZ(N3M)
    SDZ(N3) = SDZ(1)
    SSDZ(0) = SSDZ(N3M)
    SSDZ(N3) = SSDZ(1)
    VDZ(1) = 0.5*(SDZ(1) + SDZ(N3M))
    VDZ(N3) = 0.5*(SDZ(1) + SDZ(N3M))
    VVDZ(1) = 1./VDZ(1)
    VVDZ(N3) = 1./VDZ(N3)
    KMV(1) = N3M
    KPV(N3M) = 1
  END IF

  IF (XPRDIC .EQ. 1) THEN
    FIXIL(1) = 0.
    FIXIU(N1M) = 0.
    SDX(0) = SDX(N1M)
    SDX(N1) = SDX(1)
    SSDX(0) = SSDX(N1M)
    SSDX(N1) = SSDX(1)
    VDX(1) = 0.5*(SDX(1) + SDX(N1M))
    VDX(N1) = 0.5*(SDX(1) + SDX(N1M))
    VVDX(1) = 1./VDX(1)
    VVDX(N1) = 1./VDX(N1)
    IMV(1) = N1M
    IPV(N1M) = 1
  END IF

  RETURN
END
!=======================================================================
SUBROUTINE READGRID
!=======================================================================
  USE MOD_POSTINST
  IMPLICIT NONE
  INTEGER*8    :: I, J, K

  READ (11, *) N1, N2, N3
  WRITE (*, 101) N1, N2, N3
  READ (11, *) XL, YL, ZL
  WRITE (*, 102) XL, YL, ZL
101 FORMAT('NX=', I12, '  NY=', I12, '  NZ=', I12)
102 FORMAT('XL=', F12.4, '  YL=', F12.4, '  ZL=', F12.4)

  N1M = N1 - 1
  N2M = N2 - 1
  N3M = N3 - 1

  IF ((N1 .GT. M1) .OR. (N2 .GT. M2) .OR. (N3 .GT. M3)) THEN
    PRINT *, 'ARRAY SIZE CAN NOT HANDLE THIS GRID.'
    STOP
  END IF

  READ (11, *) (X(I), I=1, N1)
  READ (11, *) (Y(J), J=1, N2)
  READ (11, *) (Z(K), K=1, N3)

  RETURN
END
!=======================================================================
SUBROUTINE INDICES
!=======================================================================
  USE MOD_POSTINST
  IMPLICIT NONE
  INTEGER*8    :: I, J, K

!-----STREAMWISE DIRECTION
  DO I = 1, N1M
    IPV(I) = I + 1
    IMV(I) = I - 1
  END DO

!-----NORMAL DIRECTION
  DO J = 1, N2M
    JPV(J) = J + 1
    JMV(J) = J - 1
  END DO

!-----SPANWISE DIRECTION
  DO K = 1, N3M
    KPV(K) = K + 1
    KMV(K) = K - 1
  END DO

  RETURN
END

!=======================================================================
SUBROUTINE INDXFIX
!=======================================================================
  USE MOD_POSTINST
  IMPLICIT NONE
  INTEGER*8    :: I, J, K

  DO I = 1, N1M
    FIXIL(I) = 0.
    FIXIU(I) = 0.
  END DO
  FIXIL(1) = 1.
  FIXIU(N1M) = 1.

  DO J = 1, N2M
    FIXJL(J) = 0.
    FIXJU(J) = 0.
  END DO
  FIXJL(1) = 1.
  FIXJU(N2M) = 1.

  DO K = 1, N3M
    FIXKL(K) = 0.
    FIXKU(K) = 0.
  END DO
  FIXKL(1) = 1.
  FIXKU(N3M) = 1.

  RETURN
END
!=======================================================================
SUBROUTINE MESHES
!=======================================================================
  USE MOD_POSTINST
  IMPLICIT NONE
  INTEGER*8    :: I, J, K

  DO I = 1, N1M
    SDX(I) = X(I + 1) - X(I)
  END DO
  DO I = 2, N1M
    VDX(I) = 0.5*(SDX(I) + SDX(I - 1))
  END DO
  VDX(1) = 0.5*SDX(1)             ! LEFT BOUNDED COND.
  VDX(N1) = 0.5*SDX(N1M)          ! RIGHT BOUNDED COND.

  DO I = 1, N1M
    SSDX(I) = 1./SDX(I)
  END DO
  DO I = 1, N1
    VVDX(I) = 1./VDX(I)
  END DO

  DO J = 1, N2M
    SDY(J) = Y(J + 1) - Y(J)
  END DO
  DO J = 2, N2M
    VDY(J) = 0.5*(SDY(J) + SDY(J - 1))
  END DO
  VDY(1) = 0.5*SDY(1)             ! LOWER WALL BOUNDED COND.
  VDY(N2) = 0.5*SDY(N2M)          ! UPPER WALL BOUNDED COND.

  DO J = 1, N2M
    SSDY(J) = 1./SDY(J)
  END DO
  DO J = 1, N2
    VVDY(J) = 1./VDY(J)
  END DO

  DO K = 1, N3M
    SDZ(K) = Z(K + 1) - Z(K)
  END DO
  DO K = 2, N3M
    VDZ(K) = 0.5*(SDZ(K) + SDZ(K - 1))
  END DO
  VDZ(1) = 0.5*SDZ(1)             ! LOWER WALL BOUNDED COND.
  VDZ(N3) = 0.5*SDZ(N3M)          ! UPPER WALL BOUNDED COND.

  DO K = 1, N3M
    SSDZ(K) = 1./SDZ(K)
  END DO
  DO K = 1, N3
    VVDZ(K) = 1./VDZ(K)
  END DO

!-----set by zero
  SDX(0) = 0.
  SDX(M1) = 0.
  SDY(0) = 0.
  SDY(M2) = 0.
  SDZ(0) = 0.
  SDZ(M3) = 0.

  RETURN
END
!=======================================================================
SUBROUTINE PHYSPOS
!=======================================================================
  USE MOD_POSTINST
  IMPLICIT NONE
  INTEGER*8    :: I, J, K

  DO I = 1, N1M
    XMP(I) = X(I) + 0.5*SDX(I)
  END DO
  XMP(N1) = X(N1)
  XMP(0) = X(1)

  DO J = 1, N2M
    YMP(J) = Y(J) + 0.5*SDY(J)
  END DO
  YMP(N2) = Y(N2)
  YMP(0) = Y(1)

  DO K = 1, N3M
    ZMP(K) = Z(K) + 0.5*SDZ(K)
  END DO
  ZMP(N3) = Z(N3)
  ZMP(0) = Z(1)

  RETURN
END
!=======================================================================
SUBROUTINE GRIDUNI
!=======================================================================
  USE MOD_POSTINST
  IMPLICIT NONE
  INTEGER*8    :: I, J, K
  REAL*8       :: XINIT, XEND, YINIT, YEND, ZINIT, ZEND, UNIXL, UNIYL, UNIZL
  REAL*8       :: DXUNI, DYUNI, DZUNI

  INUMUNI = 0
  DO I = 1, N1M, 2
    INUMUNI = INUMUNI + 1
  END DO
  JNUMUNI = 0
  DO J = 1, N2M, 2
    JNUMUNI = JNUMUNI + 1
  END DO
  KNUMUNI = 0
  DO K = 1, N3M, 2
    KNUMUNI = KNUMUNI + 1
  END DO

  XINIT = X(1)
  XEND = X(N1)
  YINIT = Y(1)
  YEND = Y(N2)
  ZINIT = Z(1)
  ZEND = Z(N3)
  UNIXL = XEND - XINIT
  UNIYL = YEND - YINIT
  UNIZL = ZEND - ZINIT
  DXUNI = UNIXL/N1M
  DYUNI = UNIYL/N2M
  DZUNI = UNIZL/N3M

  WRITE (*, 11) DXUNI
  WRITE (*, 12) DYUNI
  WRITE (*, 13) DZUNI
11 FORMAT('DXUNI = ', F8.5)
12 FORMAT('DYUNI = ', F8.5)
13 FORMAT('DZUNI = ', F8.5)

  XUNI(1) = XINIT + 0.5*DXUNI
  YUNI(1) = YINIT + 0.5*DYUNI
  ZUNI(1) = ZINIT + 0.5*DZUNI

  DO I = 2, N1
    XUNI(I) = XUNI(I - 1) + DXUNI
  END DO

  DO J = 2, N2
    YUNI(J) = YUNI(J - 1) + DYUNI
  END DO

  DO K = 2, N3
    ZUNI(K) = ZUNI(K - 1) + DZUNI
  END DO

  DO I = 1, N1M
    DO J = N1, 1, -1
      IF (XUNI(I) .GE. X(J)) THEN
        FACUNIX(I) = (XUNI(I) - X(J))/(X(J + 1) - X(J))  ! FOR U(J+1)
        IUNI(I) = J
        EXIT
      END IF
    END DO
  END DO

  DO I = 1, N2M
    DO J = N2, 1, -1
      IF (YUNI(I) .GE. Y(J)) THEN
        FACUNIY(I) = (YUNI(I) - Y(J))/(Y(J + 1) - Y(J))  ! FOR U(J+1)
        JUNI(I) = J
        EXIT
      END IF
    END DO
  END DO

  DO I = 1, N3M
    DO J = N3, 1, -1
      IF (ZUNI(I) .GE. Z(J)) THEN
        FACUNIZ(I) = (ZUNI(I) - Z(J))/(Z(J + 1) - Z(J))  ! FOR U(J+1)
        KUNI(I) = J
        EXIT
      END IF
    END DO
  END DO

  RETURN
END

!=======================================================================
SUBROUTINE FINDIPJPKP
!=======================================================================
  USE MOD_POSTINST
  IMPLICIT NONE
  INTEGER*8    :: I, J, K

  DO I = 1, NIP
    DO J = N1M, 1, -1
      IF (XIP(I) .GE. XMP(J)) THEN
        IP(I) = J
        GOTO 111
      END IF
    END DO
111 CONTINUE
  END DO

  DO I = 1, NJP
    DO J = N2M, 1, -1
      IF (YJP(I) .GE. YMP(J)) THEN
        JP(I) = J
        GOTO 112
      END IF
    END DO
112 CONTINUE
  END DO

  DO I = 1, NKP
    DO J = N3M, 1, -1
      IF (ZKP(I) .GE. ZMP(J)) THEN
        KP(I) = J
        GOTO 113
      END IF
    END DO
113 CONTINUE
  END DO

  RETURN
END

!=======================================================================
SUBROUTINE SORTER(N, RA)
!=======================================================================
  REAL RA(3), RRA
  INTEGER I, IR, J, L

  IF (N .LE. 1) GOTO 30
  L = N/2 + 1
  IR = N
10 CONTINUE
  IF (L .GT. 1) THEN
    L = L - 1
    RRA = RA(L)
  ELSE
    RRA = RA(IR)
    RA(IR) = RA(1)
    IR = IR - 1
    IF (IR .EQ. 1) THEN
      RA(1) = RRA
      RETURN
    END IF
  END IF
  I = L
  J = L + L
20 IF (J .LE. IR) THEN
    IF (J .LT. IR) THEN
      IF (RA(J) .GT. RA(J + 1)) J = J + 1
    END IF
    IF (RRA .GT. RA(J)) THEN
      RA(I) = RA(J)
      I = J
      J = J + J
    ELSE
      J = IR + 1
    END IF
    GOTO 20
  END IF
  RA(I) = RRA

  GOTO 10
30 RETURN
END
!=======================================================================
!--------------------------------------------------------------
!  ROUTINE TO COMPUTE EIGENVALUES OF MATRIX
!    - ROUTINE FROM NUMERICAL RECIPE IN FORTRAN
!--------------------------------------------------------------
SUBROUTINE jacobi(a, n, np, d, v, nrot)
!=======================================================================
  INTEGER n, np, nrot, NMAX
  REAL a(np, np), d(np), v(np, np)
  PARAMETER(NMAX=500)
  INTEGER i, ip, iq, j
  REAL c, g, h, s, sm, t, tau, theta, tresh, b(NMAX), z(NMAX)
  do ip = 1, n
    do iq = 1, n
      v(ip, iq) = 0.
    end do
    v(ip, ip) = 1.
  end do
  do ip = 1, n
    b(ip) = a(ip, ip)
    d(ip) = b(ip)
    z(ip) = 0.
  end do
  nrot = 0
  do i = 1, 50
    sm = 0.
    do ip = 1, n - 1
      do iq = ip + 1, n
        sm = sm + abs(a(ip, iq))
      end do
    end do
    if (sm .eq. 0.) return
    if (i .lt. 4) then
      tresh = 0.2*sm/n**2
    else
      tresh = 0.
    end if
    do ip = 1, n - 1
      do iq = ip + 1, n
        g = 100.*abs(a(ip, iq))
        if ((i .gt. 4) .and. (abs(d(ip)) + g .eq. abs(d(ip))) .and. (abs(d(iq)) + g .eq. abs(d(iq)))) then
          a(ip, iq) = 0.
        else if (abs(a(ip, iq)) .gt. tresh) then
          h = d(iq) - d(ip)
          if (abs(h) + g .eq. abs(h)) then
            t = a(ip, iq)/h

          else
            theta = 0.5*h/a(ip, iq)
            t = 1./(abs(theta) + sqrt(1.+theta**2))
            if (theta .lt. 0.) t = -t
          end if
          c = 1./sqrt(1 + t**2)
          s = t*c
          tau = s/(1.+c)
          h = t*a(ip, iq)
          z(ip) = z(ip) - h
          z(iq) = z(iq) + h
          d(ip) = d(ip) - h
          d(iq) = d(iq) + h
          a(ip, iq) = 0.
          do j = 1, ip - 1
            g = a(j, ip)
            h = a(j, iq)

            a(j, ip) = g - s*(h + g*tau)
            a(j, iq) = h + s*(g - h*tau)
          end do
          do j = ip + 1, iq - 1
            g = a(ip, j)
            h = a(j, iq)
            a(ip, j) = g - s*(h + g*tau)
            a(j, iq) = h + s*(g - h*tau)
          end do
          do j = iq + 1, n
            g = a(ip, j)
            h = a(iq, j)
            a(ip, j) = g - s*(h + g*tau)
            a(iq, j) = h + s*(g - h*tau)
          end do
          do j = 1, n
            g = v(j, ip)

            h = v(j, iq)
            v(j, ip) = g - s*(h + g*tau)
            v(j, iq) = h + s*(g - h*tau)
          end do
          nrot = nrot + 1
        end if
      end do
    end do
    do ip = 1, n
      b(ip) = b(ip) + z(ip)
      d(ip) = b(ip)
      z(ip) = 0.
    end do
  end do
  pause 'too many iterations in jacobi'
  return
END

