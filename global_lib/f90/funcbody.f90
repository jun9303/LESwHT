!=======================================================================
      FUNCTION FUNCBODY(X,Y,Z,T)
!=======================================================================
!
!     Function for an immersed body
!
!     Required condition to use secant method :
!        1. FUNCBODY(X,Y,Z)<0  inside the body
!           FUNCBODY(X,Y,Z)=0  at the body
!           FUNCBODY(X,Y,Z)>0  outside the body
!     Ex.
!        FUNCBODY=X**2+Y**2+Z**2-0.5**2       ! for sphere
!        FUNCBODY=X**2+Y**2-0.5**2            ! for infinite cylinder
!
!     Parameter T corresponds to non-dimensional time.
!     For static-body problems, T is set zero.
!-----------------------------------------------------------------------
      IMPLICIT NONE

      REAL*8     :: X,Y,Z
      REAL*8     :: T
      REAL*8     :: FUNCBODY

      IF (ABS(Y) .GE. 0.5) THEN
        FUNCBODY = -1.
      ELSE
        FUNCBODY = 1.
      ENDIF
 
!      FUNCBODY = X**2. + Y**2. - 0.5 **2.     
      RETURN
      END FUNCTION FUNCBODY

! COORD. ROTATING MODULES FOR ROTATIONERY BODY, MADE BY SANGJOON LEE
!=======================================================================
       SUBROUTINE ROTATE_ZAXIS(XR,YR,ZR,THETA)
!=======================================================================
!     ROTATE X AND Y COORDINATE TO COUNTERCLOCKWISE DIRECTION 
!     THETA IN DEGREE UNITS
      IMPLICIT NONE
      REAL*8   , intent(in)     :: THETA
      REAL*8   , intent(inout)  :: XR,YR,ZR
      REAL*8                    :: PI
      REAL*8                    :: XTEMPO,YTEMPO

      PI = ACOS(-1.)

      XTEMPO=XR
      YTEMPO=YR

      XR=COS(THETA*PI/180.)*XTEMPO-SIN(THETA*PI/180.)&
         *YTEMPO
      YR=SIN(THETA*PI/180.)*XTEMPO+COS(THETA*PI/180.)&
         *YTEMPO
      ZR=ZR

      RETURN
      END

!=======================================================================
       SUBROUTINE ROTATE_YAXIS(XR,YR,ZR,THETA)
!=======================================================================
!     ROTATE X AND Z COORDINATE TO COUNTERCLOCKWISE DIRECTION
!     THETA IN DEGREE UNITS
      IMPLICIT NONE
      REAL*8   , intent(in)     :: THETA
      REAL*8   , intent(inout)  :: XR,YR,ZR
      REAL*8                    :: PI
      REAL*8                    :: XTEMPO,ZTEMPO

      PI = ACOS(-1.)

      XTEMPO=XR
      ZTEMPO=ZR

      XR=COS(THETA*PI/180.)*XTEMPO-SIN(THETA*PI/180.)&
         *ZTEMPO
      ZR=-SIN(THETA*PI/180.)*XTEMPO-COS(THETA*PI/180.)&
         *ZTEMPO
      YR=YR

      RETURN
      END

!=======================================================================
       SUBROUTINE ROTATE_XAXIS(XR,YR,ZR,THETA)
!=======================================================================
!     ROTATE Y AND Z COORDINATE TO COUNTERCLOCKWISE DIRECTION
!     THETA IN DEGREE UNITS
      IMPLICIT NONE
      REAL*8   , intent(in)     :: THETA
      REAL*8   , intent(inout)  :: XR,YR,ZR
      REAL*8                    :: PI
      REAL*8                    :: YTEMPO,ZTEMPO

      PI = ACOS(-1.)
      
      YTEMPO=YR
      ZTEMPO=ZR

      ZR=-COS(THETA*PI/180.)*ZTEMPO+SIN(THETA*PI/180.)&
         *YTEMPO
      YR=SIN(THETA*PI/180.)*ZTEMPO+COS(THETA*PI/180.)&
         *YTEMPO
      XR=XR

      RETURN
      END
