!=======================================================================
!
!     Copyright(c) 2017 All rights reserved by
!     Haecheon Choi
!     Department of Mechanical & Aerospace Engineering
!     Seoul National University
!     E-mail: choi@snu.ac.kr
!     URL: http://tfc.snu.ac.kr
!
!     Noncommercial use with copyrighted mark
!
!     Q&A: licacodeforum@gmail.com
!
!=======================================================================
!
!     2nd copyright(c) 2018 All right reserved by Sangjoon Lee
!     2018.02.28. Modified for F2PY(fortran-to-python) usage
!                 Sangjoon Lee from Energy and Environmentl Flow Lab.
!                 Department of Mechanical & Aerospace Engineering
!                 Seoul National University
!     Only for in-lab use. Do not distribute for commercial use.
!
!=======================================================================
      SUBROUTINE FIND_INOUT(NX,NY,NZ,XCOORD,YCOORD,ZCOORD,NBODY,INOUT,T)
!=======================================================================
!
!     Variables
!       INOUT():
!         0 => IN THE BODY,  if FUNCBODY <= 1.E-10
!         1 => OUT OF THE BODY, if FUNCBODY >  1.E-10
!         1.E-10 instead of 0. to prevent wrong surface recognizaion
!         caused by round-off error
!       FUNCBODY(): Function or point data of body
!                   It is defined in funcbody.f90.
!
!-----------------------------------------------------------------------
!$    use omp_lib
      IMPLICIT NONE
      INTEGER*8, intent(in)  :: NX,NY,NZ
      REAL*8   , intent(in)  :: XCOORD(0:NX),YCOORD(0:NY),ZCOORD(0:NZ),T
      INTEGER*8, intent(out) :: INOUT(0:NX,0:NY,0:NZ)
      INTEGER*8, intent(out) :: NBODY

      INTEGER*8    :: I,J,K
      REAL*8       :: FUNCBODY

      NBODY = 0
      INOUT = 1

!$OMP PARALLEL DO reduction(+:NBODY)
      DO 20 K = 0,NZ
      DO 20 J = 0,NY
      DO 20 I = 0,NX
      IF (FUNCBODY(XCOORD(I),YCOORD(J),ZCOORD(K),T) .LE. 1.E-10) THEN
         NBODY=NBODY+1
         INOUT(I,J,K) = 0
      ENDIF
 20   CONTINUE
!$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE FIND_INOUT

!=======================================================================
      SUBROUTINE FINDBDY_INTP(DIR, NX, NY, NZ, NBODY, INOUT, UFIX_X,   &
                              UFIX_Y, UFIX_Z, LFIX_X, LFIX_Y, LFIX_Z,  &
                              NINTP, NINNER, FCP, INTPTYPE, INTPINDX)
!=======================================================================
!
!     Variables
!     DIR : selection for velocity components (1:u, 2:v, 3:w)
!     UFIX_i : Upper fixed pts indicator in i-dir. (1: end & non-period.)
!     LFIX_i : Lower fixed pts indicator in i-dir. (1: end & non-period.)
!     INOUT <-- FIND_INOUT
!     NBODY <-- total number of pts where FUNCBODY <= 1E-10 (In the body)
!               NBODY >= NINTP + NINNER (Used for array allocation)
!
!     NINTP : number of forcing points with interpolation
!     NINNER : number of forcing points except NINTP
!     FCP[n, -] : (I, J, K) index of the n-th forcing point
!     INTPTYPE[n, 1] : interpolation type (from 1 to 3)
!     INTPINDX[n, -] : Direction where the n-th forcing point faces
!                      (-1 : opposite i-dir. / +1 : i-dir.)
!
!-----------------------------------------------------------------------
!$    use omp_lib
      IMPLICIT NONE
      INTEGER*8, intent(in)  :: DIR,NX,NY,NZ,NBODY ! DIR - 1:X, 2:Y, 3:Z
      INTEGER*8, intent(in)  :: INOUT(0:NX,0:NY,0:NZ)
      INTEGER*8, intent(in)  :: UFIX_X(NX-1),UFIX_Y(NY-1),UFIX_Z(NZ-1)
      INTEGER*8, intent(in)  :: LFIX_X(NX-1),LFIX_Y(NY-1),LFIX_Z(NZ-1)
      INTEGER*8, intent(out) :: NINTP, NINNER, FCP(NBODY,3)
      INTEGER*8, intent(out) :: INTPTYPE(NBODY,1), INTPINDX(NBODY,3)
      
      INTEGER*8    :: I,IP,IM,J,JP,JM,K,KP,KM
      INTEGER*8    :: ISTART,JSTART,KSTART
      INTEGER*8    :: INOUT_SUM, FCP_TEMP(NBODY,3)

      ISTART = 1
      JSTART = 1
      KSTART = 1
      IF (DIR .EQ. 1) THEN
         ISTART = 2
      ELSE IF (DIR .EQ. 2) THEN
         JSTART = 2
      ELSE IF (DIR .EQ. 3) THEN
         KSTART = 2
      ENDIF

      NINTP = 0
      NINNER = 0

      DO 30 K = KSTART, NZ-1
         KM = K-1
         KP = K+1
      DO 30 J = JSTART, NY-1
         JM = J-1
         JP = J+1
      DO 30 I = ISTART, NX-1
         IM = I-1
         IP = I+1

         INOUT_SUM = (ABS(INOUT(IP,J,K)-INOUT(IM,J,K)))*(1.-UFIX_X(I))*(1.-LFIX_X(I)) &
                    +(ABS(INOUT(I,JP,K)-INOUT(I,JM,K)))*(1.-UFIX_Y(J))*(1.-LFIX_Y(J)) &
                    +(ABS(INOUT(I,J,KP)-INOUT(I,J,KM)))*(1.-UFIX_Z(K))*(1.-LFIX_Z(K))

         IF ((INOUT(I,J,K) .EQ. 0) .AND. (INOUT_SUM .GT. 0)) THEN
            NINTP = NINTP + 1
            FCP(NINTP,1)=I
            FCP(NINTP,2)=J
            FCP(NINTP,3)=K
            INTPTYPE(NINTP,1)=INOUT_SUM
            INTPINDX(NINTP,1)=INOUT(IP,J,K)-INOUT(IM,J,K)
            INTPINDX(NINTP,2)=INOUT(I,JP,K)-INOUT(I,JM,K)
            INTPINDX(NINTP,3)=INOUT(I,J,KP)-INOUT(I,J,KM)
         ELSE IF ((INOUT(I,J,K) .EQ. 0) .AND. (INOUT_SUM .EQ. 0)) THEN
            NINNER = NINNER + 1
            FCP_TEMP(NINNER,1)=I
            FCP_TEMP(NINNER,2)=J
            FCP_TEMP(NINNER,3)=K
         ENDIF
 30   CONTINUE

!$OMP PARALLEL DO
      DO 40 I = 1, NINNER
         FCP(NINTP+I,1)=FCP_TEMP(I,1)
         FCP(NINTP+I,2)=FCP_TEMP(I,2)
         FCP(NINTP+I,3)=FCP_TEMP(I,3)
 40   CONTINUE
!$OMP END PARALLEL DO

      END SUBROUTINE FINDBDY_INTP

!=======================================================================
      SUBROUTINE GEOMFAC_PRESET(NI, I, IM, PRDIC, I_ADJ)
!=======================================================================
!
!     Variables
!     I : coordinates for i-gridlines
!     IM : center value between adjacent gridlines (.. X(j) | XM(j) | X(j+1) ..))
!     
!     I_ADJ : adjusted coordinates for calculating geometric factors
!             for interpolation
!             I_ADJ(_,1) for u, I_ADJ(_,2) for v, I_ADJ(_,3) for w.
!
!-----------------------------------------------------------------------
!$    use omp_lib
      IMPLICIT NONE
      INTEGER*8, intent(in)  :: NI
      REAL*8   , intent(in)  :: I(0:NI)
      REAL*8   , intent(in)  :: IM(0:NI)
      CHARACTER(len=3), intent(in) :: PRDIC
      REAL*8   , intent(out) :: I_ADJ(-1:NI+1,3)

      INTEGER*8    :: L

      I_ADJ = 0.

      DO 10 L = 1, NI
         I_ADJ(L,1) = I(L)
         I_ADJ(L,2) = IM(L)
         I_ADJ(L,3) = IM(L)
 10   CONTINUE
      I_ADJ(0,2) = IM(0)
      I_ADJ(0,3) = IM(0)

      IF (PRDIC .EQ. 'ON') THEN
         I_ADJ(0,1)=I(1)-(I(NI)-I(NI-1))
         I_ADJ(-1,1)=I(1)-(I(NI)-I(NI-2))
         I_ADJ(NI+1,1)=I(NI)+(I(2)-I(1))

         I_ADJ(0,2)=I(1)-0.5*(I(NI)-I(NI-1))
         I_ADJ(-1,2)=I(1)-(I(NI)-I(NI-1))-0.5*(I(NI-1)-I(NI-2))
         I_ADJ(NI,2)=I(NI)+0.5*(I(2)-I(1))
         I_ADJ(NI+1,2)=I(NI)+I(2)-I(1)+0.5*(I(3)-I(2))

         I_ADJ(0,3)=I_ADJ(0,2)
         I_ADJ(-1,3)=I_ADJ(-1,2)
         I_ADJ(NI,3)=I_ADJ(NI,2)
         I_ADJ(NI+1,3)=I_ADJ(NI+1,2)
      ENDIF

      END SUBROUTINE GEOMFAC_PRESET

!=======================================================================
      SUBROUTINE GEOMFAC_INTP(NX, NY, NZ, XPRE, YPRE, ZPRE, NINTP,&
                              NBODY, FCP, INTPINDX, GEOMFAC, T)
!=======================================================================
!
!     Variables
!     XYZPRE <-- GEOMFAC_PRESET
!     NINTP <-- total number of pts for interpolation
!       
!     FCP[n, -] : (I, J, K) index of the n-th forcing point 
!     INTPINDX[n, -] : Direction where the n-th forcing point faces
!                      (-1 : opposite i-dir. / +1 : i-dir.)
!     <-- FINDBDY_INTP
!
!-----------------------------------------------------------------------
!$    use omp_lib
      IMPLICIT NONE
      INTEGER*8, intent(in)  :: NX,NY,NZ,NINTP,NBODY
      REAL*8   , intent(in)  :: XPRE(-1:NX+1),YPRE(-1:NY+1),ZPRE(-1:NZ+1)
      INTEGER*8, intent(in)  :: FCP(NBODY,3),INTPINDX(NINTP,3)
      REAL*8   , intent(out) :: GEOMFAC(NINTP,3,3,3)

      INTEGER*8    :: L,M,N,I,J,K,INDIC
      REAL*8       :: X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,T
      REAL*8       :: XTEMP1,YTEMP1,ZTEMP1,XTEMP2,YTEMP2,ZTEMP2,FFS,FFE,FF1,FF2
      REAL*8       :: XX1,YY1,ZZ1,XX2,YY2,ZZ2,DDX,DDY,DDZ
      REAL*8       :: DX1,DX2,DX3,DY1,DY2,DY3,DZ1,DZ2,DZ3
      REAL*8       :: A0,B0,C0,A1,B1,C1
      REAL*8       :: FUNCBODY

!$OMP PARALLEL DO &
!$OMP private(L,M,N,I,J,K,INDIC,X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3) &
!$OMP private(XTEMP1,YTEMP1,ZTEMP1,XTEMP2,YTEMP2,ZTEMP2,FFS,FFE,FF1,FF2) &
!$OMP private(XX1,YY1,ZZ1,XX2,YY2,ZZ2,DDX,DDY,DDZ) &
!$OMP private(DX1,DX2,DX3,DY1,DY2,DY3,DZ1,DZ2,DZ3,A0,B0,C0,A1,B1,C1)
      DO 500 L = 1, NINTP
         X1 = XPRE(FCP(L,1))
         Y1 = YPRE(FCP(L,2))
         Z1 = ZPRE(FCP(L,3))

         X2 = XPRE(FCP(L,1) + INTPINDX(L,1))
         Y2 = YPRE(FCP(L,2) + INTPINDX(L,2))
         Z2 = ZPRE(FCP(L,3) + INTPINDX(L,3))

         INDIC = 0

         DO M = 1,3
            IF (M .EQ. 1) THEN
               XTEMP1 = X1
               YTEMP1 = Y1
               ZTEMP1 = Z1
               XTEMP2 = X2
               YTEMP2 = Y2
               ZTEMP2 = Z2
               FFS = FUNCBODY(X1,Y1,Z1,T)
               FFE = FUNCBODY(X2,Y2,Z2,T)
               IF (FFS*FFE .GT. 0.) THEN
                  DO I = 1,3
                     DO J = 1,3
                        DO K = 1,3
                           GEOMFAC(L,I,J,K) = 0.
                        ENDDO
                     ENDDO
                  ENDDO
                  GEOMFAC(L,1,1,1) = 1.
                  GOTO 45
               ENDIF
            ELSE
               XTEMP1 = XX1
               YTEMP1 = YY1
               ZTEMP1 = ZZ1
               XTEMP2 = XX2
               YTEMP2 = YY2
               ZTEMP2 = ZZ2
            ENDIF

            DO N = 0, 19
               DDX = XTEMP2 - XTEMP1
               DDY = YTEMP2 - YTEMP1
               DDZ = ZTEMP2 - ZTEMP1

               XX1 = XTEMP1+DDX*N/20.
               XX2 = XTEMP1+DDX*(N+1)/20.
               YY1 = YTEMP1+DDY*N/20.
               YY2 = YTEMP1+DDY*(N+1)/20.
               ZZ1 = ZTEMP1+DDZ*N/20.
               ZZ2 = ZTEMP1+DDZ*(N+1)/20.

               FF1 = FUNCBODY(XX1,YY1,ZZ1,T)
               FF2 = FUNCBODY(XX2,YY2,ZZ2,T)

               IF (FF1 .EQ. 0.) THEN
                  X0 = XX1
                  Y0 = YY1
                  Z0 = ZZ1
                  GOTO 33
               ELSEIF (FF2 .EQ. 0.) THEN
                  X0 = XX2
                  Y0 = YY2
                  Z0 = ZZ2
                  GOTO 33
               ELSEIF (FF1*FF2 .LT. 0.) THEN
                  X0 = .5*(XX1+XX2)
                  Y0 = .5*(YY1+YY2)
                  Z0 = .5*(ZZ1+ZZ2)
                  IF(M .EQ. 3) GOTO 33
                  GOTO 22
               ENDIF
            ENDDO

 22         CONTINUE
         ENDDO
 33      CONTINUE

         X3 = XPRE(FCP(L,1) + INTPINDX(L,1)*2)
         Y3 = YPRE(FCP(L,2) + INTPINDX(L,2)*2)
         Z3 = ZPRE(FCP(L,3) + INTPINDX(L,3)*2)

         DX1= ABS(X1-X0)
         DX2= ABS(X2-X0)
         DX3= ABS(X3-X0)
         DY1= ABS(Y1-Y0)
         DY2= ABS(Y2-Y0)
         DY3= ABS(Y3-Y0)
         DZ1= ABS(Z1-Z0)
         DZ2= ABS(Z2-Z0)
         DZ3= ABS(Z3-Z0)

         IF (INTPINDX(L,1) .EQ. 0) THEN
            A0=1.
            A1=1.
         ELSEIF (DX2 .GE. DX1) THEN
            A0=DX2/(DX1+DX2)
            A1=1.
         ELSEIF (DX2 .LT. DX1) THEN
            A0=0.5
            A1=(DX3-DX1)/(DX3-DX2)
         ENDIF
   
         IF (INTPINDX(L,2) .EQ. 0) THEN
            B0=1.
            B1=1.
         ELSEIF (DY2 .GE. DY1) THEN
            B0=DY2/(DY1+DY2)
            B1= 1.
         ELSEIF (DY2 .LT. DY1) THEN
            B0= 0.5
            B1=(DY3-DY1)/(DY3-DY2)
         ENDIF
   
         IF (INTPINDX(L,3) .EQ. 0) THEN
            C0=1.
            C1=1.
         ELSEIF (DZ2 .GE. DZ1) THEN
            C0=DZ2/(DZ1+DZ2)
            C1=1.
         ELSEIF (DZ2 .LT. DZ1) THEN
            C0=0.5
            C1=(DZ3-DZ1)/(DZ3-DZ2)
         ENDIF
         
         GEOMFAC(L,1,1,1)= 1./A0/B0/C0
         GEOMFAC(L,1,1,2)=-1./A0/B0/C0*(   A0)*(   B0)*(1.-C0)    &
                                    *                (   C1)
         GEOMFAC(L,1,1,3)=-1./A0/B0/C0*(   A0)*(   B0)*(1.-C0)    &
                                    *                (1.-C1)
         GEOMFAC(L,1,2,1)=-1./A0/B0/C0*(   A0)*(1.-B0)*(   C0)    &
                                    *        (   B1)
         GEOMFAC(L,1,2,2)=-1./A0/B0/C0*(   A0)*(1.-B0)*(1.-C0)    &
                                    *        (   B1)*(   C1)
         GEOMFAC(L,1,2,3)=-1./A0/B0/C0*(   A0)*(1.-B0)*(1.-C0)    &
                                    *        (   B1)*(1.-C1)
         GEOMFAC(L,1,3,1)=-1./A0/B0/C0*(   A0)*(1.-B0)*(   C0)    &
                                    *        (1.-B1)
         GEOMFAC(L,1,3,2)=-1./A0/B0/C0*(   A0)*(1.-B0)*(1.-C0)    &
                                    *        (1.-B1)*(   C1)
         GEOMFAC(L,1,3,3)=-1./A0/B0/C0*(   A0)*(1.-B0)*(1.-C0)    &
                                    *        (1.-B1)*(1.-C1)
         GEOMFAC(L,2,1,1)=-1./A0/B0/C0*(1.-A0)*(   B0)*(   C0)    &
                                    *(   A1)
         GEOMFAC(L,2,1,2)=-1./A0/B0/C0*(1.-A0)*(   B0)*(1.-C0)    &
                                    *(   A1)        *(   C1)
         GEOMFAC(L,2,1,3)=-1./A0/B0/C0*(1.-A0)*(   B0)*(1.-C0)    &
                                    *(   A1)        *(1.-C1)
         GEOMFAC(L,2,2,1)=-1./A0/B0/C0*(1.-A0)*(1.-B0)*(   C0)    &
                                    *(   A1)*(   B1)
         GEOMFAC(L,2,2,2)=-1./A0/B0/C0*(1.-A0)*(1.-B0)*(1.-C0)    &
                                    *(   A1)*(   B1)*(   C1)
         GEOMFAC(L,2,2,3)=-1./A0/B0/C0*(1.-A0)*(1.-B0)*(1.-C0)    &
                                    *(   A1)*(   B1)*(1.-C1)
         GEOMFAC(L,2,3,1)=-1./A0/B0/C0*(1.-A0)*(1.-B0)*(   C0)    &
                                    *(   A1)*(1.-B1)
         GEOMFAC(L,2,3,2)=-1./A0/B0/C0*(1.-A0)*(1.-B0)*(1.-C0)    &
                                    *(   A1)*(1.-B1)*(   C1)
         GEOMFAC(L,2,3,3)=-1./A0/B0/C0*(1.-A0)*(1.-B0)*(1.-C0)    &
                                    *(   A1)*(1.-B1)*(1.-C1)
         GEOMFAC(L,3,1,1)=-1./A0/B0/C0*(1.-A0)*(   B0)*(   C0)    &
                                    *(1.-A1)
         GEOMFAC(L,3,1,2)=-1./A0/B0/C0*(1.-A0)*(   B0)*(1.-C0)    &
                                    *(1.-A1)        *(   C1)
         GEOMFAC(L,3,1,3)=-1./A0/B0/C0*(1.-A0)*(   B0)*(1.-C0)    &
                                    *(1.-A1)        *(1.-C1)
         GEOMFAC(L,3,2,1)=-1./A0/B0/C0*(1.-A0)*(1.-B0)*(   C0)    &
                                    *(1.-A1)*(   B1)
         GEOMFAC(L,3,2,2)=-1./A0/B0/C0*(1.-A0)*(1.-B0)*(1.-C0)    &
                                    *(1.-A1)*(   B1)*(   C1)
         GEOMFAC(L,3,2,3)=-1./A0/B0/C0*(1.-A0)*(1.-B0)*(1.-C0)    &
                                    *(1.-A1)*(   B1)*(1.-C1)
         GEOMFAC(L,3,3,1)=-1./A0/B0/C0*(1.-A0)*(1.-B0)*(   C0)    &
                                    *(1.-A1)*(1.-B1)
         GEOMFAC(L,3,3,2)=-1./A0/B0/C0*(1.-A0)*(1.-B0)*(1.-C0)    &
                                    *(1.-A1)*(1.-B1)*(   C1)
         GEOMFAC(L,3,3,3)=-1./A0/B0/C0*(1.-A0)*(1.-B0)*(1.-C0)    &
                                    *(1.-A1)*(1.-B1)*(1.-C1)
 45   CONTINUE
 500  CONTINUE
!$OMP END PARALLEL DO

      END SUBROUTINE GEOMFAC_INTP

!=======================================================================
      SUBROUTINE FINDBDY_NOINTP(DIR, NX, NY, NZ, NBODY, INOUT, NINNER, &
                                FCP)
!=======================================================================
!
!     Variables
!     DIR : selection for velocity components (1:u, 2:v, 3:w)
!     INOUT <-- FIND_INOUT
!     NBODY <-- total number of pts where FUNCBODY <= 1E-10 (In the body)
!               NBODY >= NINNER (Used for array allocation)
!
!     NINNER : number of forcing points except NINTP
!     FCP[n, -] : (I, J, K) index of the n-th forcing point
!
!-----------------------------------------------------------------------
!$    use omp_lib
      IMPLICIT NONE
      INTEGER*8, intent(in)  :: DIR,NX,NY,NZ,NBODY ! DIR - 1:X, 2:Y, 3:Z
      INTEGER*8, intent(in)  :: INOUT(0:NX,0:NY,0:NZ)
      INTEGER*8, intent(out) :: NINNER, FCP(NBODY,3)

      INTEGER*8    :: I,J,K
      INTEGER*8    :: ISTART,JSTART,KSTART

      ISTART = 1
      JSTART = 1
      KSTART = 1
      IF (DIR .EQ. 1) THEN
         ISTART = 2
      ELSE IF (DIR .EQ. 2) THEN
         JSTART = 2
      ELSE IF (DIR .EQ. 3) THEN
         KSTART = 2
      ENDIF

      NINNER = 0

      DO 30 K = KSTART, NZ-1
      DO 30 J = JSTART, NY-1
      DO 30 I = ISTART, NX-1
         IF (INOUT(I,J,K) .EQ. 0) THEN
            NINNER = NINNER + 1
            FCP(NINNER,1)=I
            FCP(NINNER,2)=J
            FCP(NINNER,3)=K
         ENDIF
 30   CONTINUE

      END SUBROUTINE FINDBDY_NOINTP

!=======================================================================
      SUBROUTINE FIND_ZERO_NU_SGS(NX,NY,NZ,XM,YM,ZM,NZERO,INOUT,T)
!=======================================================================
!
!     Variables
!       INOUT():
!         0 => SGS EDDY VISCOSITY SHOULD BE 0
!         1 => SGS EDDY VISCOSITY IS NOT 0
!         1.E-10 instead of 0. to prevent wrong surface recognizaion
!         caused by round-off error
!       FUNCBODY(): Function or point data of body
!                   It is defined in funcbody.f90.
!
!-----------------------------------------------------------------------
!$    use omp_lib
      IMPLICIT NONE
      INTEGER*8, intent(in)  :: NX,NY,NZ
      REAL*8   , intent(in)  :: XM(0:NX),YM(0:NY),ZM(0:NZ)
      INTEGER*8, intent(out) :: INOUT(1:NX-1,1:NY-1,1:NZ-1)
      INTEGER*8, intent(out) :: NZERO

      INTEGER*8    :: I,J,K
      REAL*8       :: FUNCBODY, T

      NZERO = 0
      INOUT = 1

!$OMP PARALLEL DO reduction(+:NZERO)
      DO 20 K = 1,NZ-1
      DO 20 J = 1,NY-1
      DO 20 I = 1,NX-1
      IF ((FUNCBODY(XM(I  ),YM(J  ),ZM(K  ),T) .LE. 1.E-10) .OR. &
          (FUNCBODY(XM(I-1),YM(J  ),ZM(K  ),T) .LE. 1.E-10) .OR. &
          (FUNCBODY(XM(I+1),YM(J  ),ZM(K  ),T) .LE. 1.E-10) .OR. &
          (FUNCBODY(XM(I  ),YM(J-1),ZM(K  ),T) .LE. 1.E-10) .OR. &
          (FUNCBODY(XM(I  ),YM(J+1),ZM(K  ),T) .LE. 1.E-10) .OR. &
          (FUNCBODY(XM(I  ),YM(J  ),ZM(K-1),T) .LE. 1.E-10) .OR. &
          (FUNCBODY(XM(I  ),YM(J  ),ZM(K+1),T) .LE. 1.E-10)) THEN
         NZERO=NZERO+1
         INOUT(I,J,K) = 0
      ENDIF
 20   CONTINUE
!$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE FIND_ZERO_NU_SGS

!=======================================================================
      SUBROUTINE CONJG_INTP(NX,NY,NZ,CRATIO,KRATIO,XM,YM,ZM,X,Y,Z,     &
                            ISZERO,T,CSTAR,KSTAR)
!=======================================================================
!$    use omp_lib
      IMPLICIT NONE
      INTEGER*8, intent(in)  :: NX,NY,NZ
      REAL*8   , intent(in)  :: CRATIO,KRATIO
      REAL*8   , intent(in)  :: XM(0:NX),YM(0:NY),ZM(0:NZ)
      REAL*8   , intent(in)  :: X(0:NX),Y(0:NY),Z(0:NZ)
      INTEGER*8, intent(in)  :: ISZERO(1:NX-1,1:NY-1,1:NZ-1)
      REAL*8   , intent(out) :: CSTAR(1:NX-1,1:NY-1,1:NZ-1)
      REAL*8   , intent(out) :: KSTAR(1:NX-1,1:NY-1,1:NZ-1,6)

      INTEGER*8    :: I,J,K,L,SUBC
      REAL*8       :: FPTEMP
      REAL*8       :: AA
      REAL*8       :: FUNCBODY, FLUID_PORTION, K_SAM_RGN, K_INT_CHG, T

      SUBC = 1

!$OMP PARALLEL DO &
!$OMP private(FPTEMP,AA)
      DO 20 K = 1,NZ-1
      DO 20 J = 1,NY-1
      DO 20 I = 1,NX-1
         AA = FUNCBODY(XM(I  ),YM(J  ),ZM(K  ),T)   ! FROM SOLID REGION

         FPTEMP = FLUID_PORTION(X(I),X(I+1),Y(J),Y(J+1),Z(K),Z(K+1),T,SUBC)
         CSTAR(I,J,K) = (1.-FPTEMP)*CRATIO + FPTEMP*1.

         FPTEMP = FLUID_PORTION(XM(I),XM(I+1),Y(J),Y(J+1),Z(K),Z(K+1),T,SUBC)
         IF (AA*FUNCBODY(XM(I+1),YM(J  ),ZM(K  ),T) .GE. 0.) THEN
            KSTAR(I,J,K,1) = (1.-FPTEMP)*KRATIO + FPTEMP*1.
         ELSE
            KSTAR(I,J,K,1) = KRATIO / (KRATIO*FPTEMP + 1.*(1.-FPTEMP))
         ENDIF

         FPTEMP = FLUID_PORTION(XM(I-1),XM(I),Y(J),Y(J+1),Z(K),Z(K+1),T,SUBC)
         IF (AA*FUNCBODY(XM(I-1),YM(J  ),ZM(K  ),T) .GE. 0.) THEN
            KSTAR(I,J,K,2) = (1.-FPTEMP)*KRATIO + FPTEMP*1.
         ELSE
            KSTAR(I,J,K,2) = KRATIO / (KRATIO*FPTEMP + 1.*(1.-FPTEMP))
         ENDIF

         FPTEMP = FLUID_PORTION(X(I),X(I+1),YM(J),YM(J+1),Z(K),Z(K+1),T,SUBC)
         IF (AA*FUNCBODY(XM(I  ),YM(J+1),ZM(K  ),T) .GE. 0.) THEN
            KSTAR(I,J,K,3) = (1.-FPTEMP)*KRATIO + FPTEMP*1.
         ELSE
            KSTAR(I,J,K,3) = KRATIO / (KRATIO*FPTEMP + 1.*(1.-FPTEMP))
         ENDIF

         FPTEMP = FLUID_PORTION(X(I),X(I+1),YM(J-1),YM(J),Z(K),Z(K+1),T,SUBC)
         IF (AA*FUNCBODY(XM(I  ),YM(J-1),ZM(K  ),T) .GE. 0.) THEN
            KSTAR(I,J,K,4) = (1.-FPTEMP)*KRATIO + FPTEMP*1.
         ELSE
            KSTAR(I,J,K,4) = KRATIO / (KRATIO*FPTEMP + 1.*(1.-FPTEMP))
         ENDIF

         FPTEMP = FLUID_PORTION(X(I),X(I+1),Y(J),Y(J+1),ZM(K),ZM(K+1),T,SUBC)
         IF (AA*FUNCBODY(XM(I  ),YM(J  ),ZM(K+1),T) .GE. 0.) THEN
            KSTAR(I,J,K,5) = (1.-FPTEMP)*KRATIO + FPTEMP*1.
         ELSE
            KSTAR(I,J,K,5) = KRATIO / (KRATIO*FPTEMP + 1.*(1.-FPTEMP))
         ENDIF

         FPTEMP = FLUID_PORTION(X(I),X(I+1),Y(J),Y(J+1),ZM(K-1),ZM(K),T,SUBC)
         IF (AA*FUNCBODY(XM(I  ),YM(J  ),ZM(K-1),T) .GE. 0.) THEN
            KSTAR(I,J,K,6) = (1.-FPTEMP)*KRATIO + FPTEMP*1.
         ELSE
            KSTAR(I,J,K,6) = KRATIO / (KRATIO*FPTEMP + 1.*(1.-FPTEMP))
         ENDIF

 20   CONTINUE
!$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE CONJG_INTP
!=======================================================================
      FUNCTION FLUID_PORTION(X1,X2,Y1,Y2,Z1,Z2,T,DIV)
!=======================================================================
!$    use omp_lib
      IMPLICIT NONE
      REAL*8       :: X1,X2,Y1,Y2,Z1,Z2,T
      INTEGER*8    :: DIV
      REAL*8       :: FLUID_PORTION

      REAL*8       :: XTEMP,YTEMP,ZTEMP,FUNCBODY
      INTEGER*8    :: I,J,K,SOLIDCELL

      SOLIDCELL = 0

!$OMP PARALLEL DO reduction(+:SOLIDCELL) private(XTEMP,YTEMP,ZTEMP)
      DO I = 1,DIV
         DO J = 1,DIV
            DO K = 1,DIV
               XTEMP = X1 + (X2-X1)/(DIV*1.)*(2.*I-1.)/2.
               YTEMP = Y1 + (Y2-Y1)/(DIV*1.)*(2.*J-1.)/2.
               ZTEMP = Z1 + (Z2-Z1)/(DIV*1.)*(2.*K-1.)/2.
               IF (FUNCBODY(XTEMP,YTEMP,ZTEMP,T).LE.1.E-10) THEN
                 SOLIDCELL = SOLIDCELL + 1
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      FLUID_PORTION = 1. - SOLIDCELL/(DIV**3.)

      RETURN
      END FUNCTION FLUID_PORTION
