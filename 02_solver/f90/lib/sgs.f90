!=======================================================================
      SUBROUTINE SGSFILTERINIT
!=======================================================================
!
!     Filter coefficients for the Simpson's rule
!     : For dynamic procedures to determine the model coefficient,
!       a test filter is applied based on the Simpson's rule.
!
!     CFX1, CFX2  : Filter coefficients in x-direction
!     CFY1, CFY2  : Filter coefficients in y-direction
!     CFZ1, CFZ2  : Filter coefficients in z-direction
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : CFX1, CFX2, CFY1, CFY2, CFZ1, CFZ2 
      IMPLICIT NONE
      INTEGER*8              :: I,J,K
      REAL*8                 :: H1,H2,H3,H4

      CFX1 = 0.
      CFX2 = 0.
      CFY1 = 0.
      CFY2 = 0.
      CFZ1 = 0.
      CFZ2 = 0.

      IF (XPRDIC .EQ. 1) THEN
      
        IF (N1M.EQ.1) GOTO 101

!$OMP PARALLEL DO private(H1,H2)
        DO I=1,N1M
          H1=C2CX(I)
          H2=C2CX(I+1)
          CFX1(I,-1)=(2.*H1-H2)/6./H1
          CFX1(I, 0)=((H1+H2)**2)/6./H1/H2
          CFX1(I, 1)=(2.*H2-H1)/6./H2
        ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO private(H1,H2,H3,H4)
        DO I=2,N1M-1
          H1=C2CX(I-1)
          H2=C2CX(I  )
          H3=C2CX(I+1)
          H4=C2CX(I+2)
          CFX2(I,-2)=(2.*H1-H2)/6./H1/2.
          CFX2(I,-1)=((H1+H2)**2)/6./H1/H2/2.
          CFX2(I, 0)=((2.*H2-H1)/6./H2+(2.*H3-H4)/6./H3)/2.
          CFX2(I, 1)=((H3+H4)**2)/6./H3/H4/2.
          CFX2(I, 2)=(2.*H4-H3)/6./H4/2.
        ENDDO
!$OMP END PARALLEL DO
        H1=C2CX(N1M)
        H2=C2CX(1)
        H3=C2CX(2)
        H4=C2CX(3)
        CFX2(1,-2)=(2.*H1-H2)/6./H1/2.
        CFX2(1,-1)=((H1+H2)**2)/6./H1/H2/2.
        CFX2(1, 0)=((2.*H2-H1)/6./H2+(2.*H3-H4)/6./H3)/2.
        CFX2(1, 1)=((H3+H4)**2)/6./H3/H4/2.
        CFX2(1, 2)=(2.*H4-H3)/6./H4/2.
        H1=C2CX(N1M-1)
        H2=C2CX(N1M)
        H3=C2CX(1)
        H4=C2CX(2)
        CFX2(N1M,-2)=(2.*H1-H2)/6./H1/2.
        CFX2(N1M,-1)=((H1+H2)**2)/6./H1/H2/2.
        CFX2(N1M, 0)=((2.*H2-H1)/6./H2+(2.*H3-H4)/6./H3)/2.
        CFX2(N1M, 1)=((H3+H4)**2)/6./H3/H4/2.
        CFX2(N1M, 2)=(2.*H4-H3)/6./H4/2.

      ELSE

        IF (N1M.EQ.1) GOTO 101

!$OMP PARALLEL DO private(H1,H2)
        DO I=1,N1M
          H1=C2CX(I)
          H2=C2CX(I+1)
          CFX1(I,-1)=(2.*H1-H2)/6./H1
          CFX1(I, 0)=((H1+H2)**2)/6./H1/H2
          CFX1(I, 1)=(2.*H2-H1)/6./H2
        ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO private(H1,H2,H3,H4)
        DO I=3,N1M-2
          H1=C2CX(I-1)
          H2=C2CX(I  )
          H3=C2CX(I+1)
          H4=C2CX(I+2)
          CFX2(I,-2)=(2.*H1-H2)/6./H1/2.
          CFX2(I,-1)=((H1+H2)**2)/6./H1/H2/2.
          CFX2(I, 0)=((2.*H2-H1)/6./H2+(2.*H3-H4)/6./H3)/2.
          CFX2(I, 1)=((H3+H4)**2)/6./H3/H4/2.
          CFX2(I, 2)=(2.*H4-H3)/6./H4/2.
        ENDDO
!$OMP END PARALLEL DO
        H1=C2CX(2)
        H2=C2CX(3)
        CFX2(2,-1)=(2.*H1-H2)/6./H1
        CFX2(2, 0)=((H1+H2)**2)/6./H1/H2
        CFX2(2, 1)=(2.*H2-H1)/6./H2
        H1=C2CX(N1M-1)
        H2=C2CX(N1M)
        CFX2(N1M-1,-1)=(2.*H1-H2)/6./H1
        CFX2(N1M-1, 0)=((H1+H2)**2)/6./H1/H2
        CFX2(N1M-1, 1)=(2.*H2-H1)/6./H2

      ENDIF

  101 CONTINUE


      IF (YPRDIC .EQ. 1) THEN
      
        IF (N2M.EQ.1) GOTO 102

!$OMP PARALLEL DO private(H1,H2)
        DO J=1,N2M
          H1=C2CY(J)
          H2=C2CY(J+1)
          CFY1(J,-1)=(2.*H1-H2)/6./H1
          CFY1(J, 0)=((H1+H2)**2)/6./H1/H2
          CFY1(J, 1)=(2.*H2-H1)/6./H2
        ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO private(H1,H2,H3,H4)
        DO J=2,N2M-1
          H1=C2CY(J-1)
          H2=C2CY(J  )
          H3=C2CY(J+1)
          H4=C2CY(J+2)
          CFY2(J,-2)=(2.*H1-H2)/6./H1/2.
          CFY2(J,-1)=((H1+H2)**2)/6./H1/H2/2.
          CFY2(J, 0)=((2.*H2-H1)/6./H2+(2.*H3-H4)/6./H3)/2.
          CFY2(J, 1)=((H3+H4)**2)/6./H3/H4/2.
          CFY2(J, 2)=(2.*H4-H3)/6./H4/2.
        ENDDO
!$OMP END PARALLEL DO
        H1=C2CY(N2M)
        H2=C2CY(1)
        H3=C2CY(2)
        H4=C2CY(3)
        CFY2(1,-2)=(2.*H1-H2)/6./H1/2.
        CFY2(1,-1)=((H1+H2)**2)/6./H1/H2/2.
        CFY2(1, 0)=((2.*H2-H1)/6./H2+(2.*H3-H4)/6./H3)/2.
        CFY2(1, 1)=((H3+H4)**2)/6./H3/H4/2.
        CFY2(1, 2)=(2.*H4-H3)/6./H4/2.
        H1=C2CY(N2M-1)
        H2=C2CY(N2M)
        H3=C2CY(1)
        H4=C2CY(2)
        CFY2(N2M,-2)=(2.*H1-H2)/6./H1/2.
        CFY2(N2M,-1)=((H1+H2)**2)/6./H1/H2/2.
        CFY2(N2M, 0)=((2.*H2-H1)/6./H2+(2.*H3-H4)/6./H3)/2.
        CFY2(N2M, 1)=((H3+H4)**2)/6./H3/H4/2.
        CFY2(N2M, 2)=(2.*H4-H3)/6./H4/2.

      ELSE

        IF (N2M.EQ.1) GOTO 102
!$OMP PARALLEL DO private(H1,H2)
        DO J=1,N2M
          H1=C2CY(J)
          H2=C2CY(J+1)
          CFY1(J,-1)=(2.*H1-H2)/6./H1
          CFY1(J, 0)=((H1+H2)**2)/6./H1/H2
          CFY1(J, 1)=(2.*H2-H1)/6./H2

        ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO private(H1,H2,H3,H4)
        DO J=3,N2M-2
          H1=C2CY(J-1)
          H2=C2CY(J  )
          H3=C2CY(J+1)
          H4=C2CY(J+2)
          CFY2(J,-2)=(2.*H1-H2)/6./H1/2.
          CFY2(J,-1)=((H1+H2)**2)/6./H1/H2/2.
          CFY2(J, 0)=((2.*H2-H1)/6./H2+(2.*H3-H4)/6./H3)/2.
          CFY2(J, 1)=((H3+H4)**2)/6./H3/H4/2.
          CFY2(J, 2)=(2.*H4-H3)/6./H4/2.
        ENDDO
!$OMP END PARALLEL DO
        H1=C2CY(2)
        H2=C2CY(3)
        CFY2(2,-1)=(2.*H1-H2)/6./H1
        CFY2(2, 0)=((H1+H2)**2)/6./H1/H2
        CFY2(2, 1)=(2.*H2-H1)/6./H2
        H1=C2CY(N2M-1)
        H2=C2CY(N2M)
        CFY2(N2M-1,-1)=(2.*H1-H2)/6./H1
        CFY2(N2M-1, 0)=((H1+H2)**2)/6./H1/H2
        CFY2(N2M-1, 1)=(2.*H2-H1)/6./H2
      ENDIF

  102 CONTINUE
  
      IF (ZPRDIC .EQ. 1) THEN

        IF (N3M.EQ.1) GOTO 103

!$OMP PARALLEL DO private(H1,H2)
        DO K=1,N3M
          H1=C2CZ(K)
          H2=C2CZ(K+1)
          CFZ1(K,-1)=(2.*H1-H2)/6./H1
          CFZ1(K, 0)=((H1+H2)**2)/6./H1/H2
          CFZ1(K, 1)=(2.*H2-H1)/6./H2
        ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO private(H1,H2,H3,H4)
        DO K=2,N3M-1
          H1=C2CZ(K-1)
          H2=C2CZ(K  )
          H3=C2CZ(K+1)
          H4=C2CZ(K+2)
          CFZ2(K,-2)=(2.*H1-H2)/6./H1/2.
          CFZ2(K,-1)=((H1+H2)**2)/6./H1/H2/2.
          CFZ2(K, 0)=((2.*H2-H1)/6./H2+(2.*H3-H4)/6./H3)/2.
          CFZ2(K, 1)=((H3+H4)**2)/6./H3/H4/2.
          CFZ2(K, 2)=(2.*H4-H3)/6./H4/2.
        ENDDO
!$OMP END PARALLEL DO
        H1=C2CZ(N3M)
        H2=C2CZ(1)
        H3=C2CZ(2)
        H4=C2CZ(3)
        CFZ2(1,-2)=(2.*H1-H2)/6./H1/2.
        CFZ2(1,-1)=((H1+H2)**2)/6./H1/H2/2.
        CFZ2(1, 0)=((2.*H2-H1)/6./H2+(2.*H3-H4)/6./H3)/2.
        CFZ2(1, 1)=((H3+H4)**2)/6./H3/H4/2.
        CFZ2(1, 2)=(2.*H4-H3)/6./H4/2.
        H1=C2CZ(N3M-1)
        H2=C2CZ(N3M)
        H3=C2CZ(1)
        H4=C2CZ(2)
        CFZ2(N3M,-2)=(2.*H1-H2)/6./H1/2.
        CFZ2(N3M,-1)=((H1+H2)**2)/6./H1/H2/2.
        CFZ2(N3M, 0)=((2.*H2-H1)/6./H2+(2.*H3-H4)/6./H3)/2.
        CFZ2(N3M, 1)=((H3+H4)**2)/6./H3/H4/2.
        CFZ2(N3M, 2)=(2.*H4-H3)/6./H4/2.

      ELSE

        IF (N3M.EQ.1) GOTO 103

!$OMP PARALLEL DO private(H1,H2)
        DO K=1,N3M
          H1=C2CZ(K)
          H2=C2CZ(K+1)
          CFZ1(K,-1)=(2.*H1-H2)/6./H1
          CFZ1(K, 0)=((H1+H2)**2)/6./H1/H2
          CFZ1(K, 1)=(2.*H2-H1)/6./H2
        ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO private(H1,H2,H3,H4)
        DO K=3,N3M-2
          H1=C2CZ(K-1)
          H2=C2CZ(K  )
          H3=C2CZ(K+1)
          H4=C2CZ(K+2)
          CFZ2(K,-2)=(2.*H1-H2)/6./H1/2.
          CFZ2(K,-1)=((H1+H2)**2)/6./H1/H2/2.
          CFZ2(K, 0)=((2.*H2-H1)/6./H2+(2.*H3-H4)/6./H3)/2.
          CFZ2(K, 1)=((H3+H4)**2)/6./H3/H4/2.
          CFZ2(K, 2)=(2.*H4-H3)/6./H4/2.
        ENDDO
!$OMP END PARALLEL DO
        H1=C2CZ(2)
        H2=C2CZ(3)
        CFZ2(2,-1)=(2.*H1-H2)/6./H1
        CFZ2(2, 0)=((H1+H2)**2)/6./H1/H2
        CFZ2(2, 1)=(2.*H2-H1)/6./H2
        H1=C2CZ(N3M-1)
        H2=C2CZ(N3M)
        CFZ2(N3M-1,-1)=(2.*H1-H2)/6./H1
        CFZ2(N3M-1, 0)=((H1+H2)**2)/6./H1/H2
        CFZ2(N3M-1, 1)=(2.*H2-H1)/6./H2

      ENDIF

  103 CONTINUE

      RETURN
      END SUBROUTINE SGSFILTERINIT
!=======================================================================
!=======================================================================
      SUBROUTINE SGSCALC
!=======================================================================
!
! CALCULATE SUBGRID-SCALE EDDY VISCOSITY, NUSGS
!
! INSMDL = 0 -> SMAGORINSKY MODEL (SMAGORINSKY,  1963. MON.WEATHER REV)
! ###### 2018.02.22. SMAGORINKSY MODEL IS NOT IMPLEMENTED YET
!
! INSMDL = 1 -> VREMAN MODEL      (VREMAN,       2004. PHYS.FLUIDS)
!
! IDVMON = OFF, CONSTANT COEFF. MODEL
! IDVMON = ON , DYNAMIC COEFF. PROCEDURE
! (DSLM MODEL: LILLY,       1991. PHYS.FLUIDS A)
! (DVMG MODEL: PARK ET AL., 2006. PHYS.FLUIDS)
!
! PLZ REFER TO "GERMANO IDENTITY" THEORY (GERMANO, 1992. J.FLUID MECH.)
!
!-----------------------------------------------------------------------
!$    use omp_lib
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : NUSGS
      IMPLICIT NONE
      IF (MSUB .EQ. 1) THEN
        IF (IREAD.EQ.0 .AND. NTIME.EQ.1) THEN
          NUSGS = 0D0
        ELSE
          IF (INSMDL .EQ. 0) THEN
            WRITE(*,*) 'SMAGORINSKY MODEL IS NOT IMPLEMENTED YET'
            WRITE(*,*) 'SGS.F90 LINE 317'            
            STOP
          ELSE IF (INSMDL .EQ. 1) THEN
            IF (IDVMON .EQ. 1) THEN
              CALL SGS_DVMG         ! DYNAMIC VREMAN MODEL W/ GLOBAL COEF
            ELSE
              CALL SGS_CVM          ! CONSTANT VREMAN MODEL   
            ENDIF
          ENDIF
        ENDIF

        CALL NUTZERO               ! SET NUSGS TO ZERO IN THE SOLID BODY
        CALL NUTINTERPOL           !              INTERPOLATION OF NUSGS
      ENDIF

      RETURN
      END SUBROUTINE SGSCALC
!=======================================================================
!=======================================================================
      SUBROUTINE SGS_CVM
!=======================================================================
!
!     NUSGS=CSGSTS*SQRT(BBB/AAA)
!     CSGSTS : VREMAN CONSTANT (= CSGSTS FROM MOD_COMMON)
!     BBB    : B_{11}*B_{22}-B_{12}*B_{12}+B_{11}*B_{33}
!             -B_{13}*B_{13}+B_{22}*B_{33}-B_{23}*B_{23}
!     AAA    : A_{ij}*A_{ij}
!     B_{ij} : SIGMA(DEL_{m}^2*A_{mi}*A_{mj}) where m=1, 2, 3
!     A_{ij} : DUj/DXi
!
!-----------------------------------------------------------------------  
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,NUSGS
      IMPLICIT NONE
      INTEGER*8 :: I,J,K,KPLUS,KMINUS,JPLUS,JMINUS,IPLUS,IMINUS
      REAL*8    :: VG11,VG12,VG13,VG21,VG22,VG23,VG31,VG32,VG33
      REAL*8    :: UP,UM,VP,VM,WP,WM
      REAL*8    :: A(9),B(6)
      REAL*8    :: AAA,BBB

!$OMP PARALLEL DO private(UP,UM,VP,VM,WP,WM)&
!$OMP private(A,B,AAA,BBB)
      DO K=1,N3M
        DO J=1,N2M
          DO I=1,N1M

            CALL VELGRAD(I,J,K,A,B)         ! B IS USED FOR DUMMY VARIABLE

            B(1)=F2FX(I)**2.*A(1)*A(1)  &
                +F2FY(J)**2.*A(4)*A(4)  &
                +F2FZ(K)**2.*A(7)*A(7)
            B(2)=F2FX(I)**2.*A(1)*A(2)  &
                +F2FY(J)**2.*A(4)*A(5)  &
                +F2FZ(K)**2.*A(7)*A(8)
            B(3)=F2FX(I)**2.*A(1)*A(3)  &
                +F2FY(J)**2.*A(4)*A(6)  &
                +F2FZ(K)**2.*A(7)*A(9)
            B(4)=F2FX(I)**2.*A(2)*A(2)  &
                +F2FY(J)**2.*A(5)*A(5)  &
                +F2FZ(K)**2.*A(8)*A(8)
            B(5)=F2FX(I)**2.*A(2)*A(3)  &
                +F2FY(J)**2.*A(5)*A(6)  &
                +F2FZ(K)**2.*A(8)*A(9)
            B(6)=F2FX(I)**2.*A(3)*A(3)  &
                +F2FY(J)**2.*A(6)*A(6)  &
                +F2FZ(K)**2.*A(9)*A(9)
      
            AAA=A(1)**2.+A(2)**2.+A(3)**2. &
               +A(4)**2.+A(5)**2.+A(6)**2. &
               +A(7)**2.+A(8)**2.+A(9)**2.
            BBB=B(1)*B(4)+B(4)*B(6)+B(6)*B(1)&
               -B(2)**2.-B(3)**2.-B(5)**2.
      
            BBB = DMAX1(BBB,1.0E-16)
   
            IF (AAA .EQ. 0.) THEN
              NUSGS(I,J,K) = 0D0   
            ELSE
              NUSGS(I,J,K)= CSGSTS*SQRT(BBB/AAA)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE SGS_CVM
!=======================================================================
!=======================================================================
      SUBROUTINE SGS_DVMG
!=======================================================================
!
!     NUSGS=CSGSTS*SQRT(BBB/AAA)
!     CSGSTS   : VREMAN CONSTANT, DETERMINED FROM LSM OF GERMANO IDENTITY
!     BBB    : B_{11}*B_{22}-B_{12}*B_{12}+B_{11}*B_{33}
!             -B_{13}*B_{13}+B_{22}*B_{33}-B_{23}*B_{23}
!     AAA    : A_{ij}*A_{ij}
!     B_{ij} : SIGMA(DEL_{m}^2*A_{mi}*A_{mj}) where m=1, 2, 3
!     A_{ij} : DUj/DXi
!
!----------------------------------------------------------------------- 
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,NUSGS,NWALL_DVM,AALP,LLIJ,MMIJ,UUI
      IMPLICIT NONE
      INTEGER*8     :: I,J,K
      REAL*8  :: LJMJ_V,MJMJ_V
      REAL*8  :: SDXF2,SDYF2,SDZF2,AAA,BBB,AMI,VOLUME
      REAL*8  :: B(6),STR(6),ALP(9)

      LJMJ_V = 0.
      MJMJ_V = 0.
      CSGSTS     = 0.

      ALLOCATE (AALP(N1M,N2M,N3M,9))
      ALLOCATE (LLIJ(N1M,N2M,N3M,6))
      ALLOCATE (MMIJ(N1M,N2M,N3M,6))
      ALLOCATE (UUI (N1M,N2M,N3M,3))

!$OMP PARALLEL DO&
!$OMP private(B,STR,ALP)&
!$OMP private(SDXF2,SDYF2,SDZF2,AAA,BBB)
      DO K=1,N3M
        DO J=1,N2M
          DO I=1,N1M
            CALL VELGRAD(I,J,K,ALP,STR)
            SDXF2=F2FX(I)**2.
            SDYF2=F2FY(J)**2.
            SDZF2=F2FZ(K)**2.
            B(1)=SDXF2*ALP(1)*ALP(1)  &
                +SDYF2*ALP(4)*ALP(4)  &
                +SDZF2*ALP(7)*ALP(7)
            B(2)=SDXF2*ALP(1)*ALP(2)  &
                +SDYF2*ALP(4)*ALP(5)  &
                +SDZF2*ALP(7)*ALP(8)
            B(3)=SDXF2*ALP(1)*ALP(3)  &
                +SDYF2*ALP(4)*ALP(6)  &
                +SDZF2*ALP(7)*ALP(9)
            B(4)=SDXF2*ALP(2)*ALP(2)  &
                +SDYF2*ALP(5)*ALP(5)  &
                +SDZF2*ALP(8)*ALP(8)
            B(5)=SDXF2*ALP(2)*ALP(3)  &
                +SDYF2*ALP(5)*ALP(6)  &
                +SDZF2*ALP(8)*ALP(9)
            B(6)=SDXF2*ALP(3)*ALP(3)  &
                +SDYF2*ALP(6)*ALP(6)  &
                +SDZF2*ALP(9)*ALP(9)
            AAA=ALP(1)**2.+ALP(2)**2.+ALP(3)**2. &
               +ALP(4)**2.+ALP(5)**2.+ALP(6)**2. &
               +ALP(7)**2.+ALP(8)**2.+ALP(9)**2.
            BBB=B(1)*B(4)+B(4)*B(6)+B(6)*B(1)  &
               -B(2)**2. -B(3)**2. -B(5)**2.
            BBB = DMAX1(BBB,1.0E-16)
            IF (AAA .EQ. 0.) THEN
              NUSGS(I,J,K) = 0D0   
            ELSE
              NUSGS(I,J,K)= SQRT(BBB/AAA)
            ENDIF
            MMIJ(I,J,K,1)=NUSGS(I,J,K)*STR(1)
            MMIJ(I,J,K,2)=NUSGS(I,J,K)*STR(2)
            MMIJ(I,J,K,3)=NUSGS(I,J,K)*STR(3)
            MMIJ(I,J,K,4)=NUSGS(I,J,K)*STR(4)
            MMIJ(I,J,K,5)=NUSGS(I,J,K)*STR(5)
            MMIJ(I,J,K,6)=NUSGS(I,J,K)*STR(6)
            AALP(I,J,K,1)=ALP(1)
            AALP(I,J,K,2)=ALP(2)
            AALP(I,J,K,3)=ALP(3)
            AALP(I,J,K,4)=ALP(4)
            AALP(I,J,K,5)=ALP(5)
            AALP(I,J,K,6)=ALP(6)
            AALP(I,J,K,7)=ALP(7)
            AALP(I,J,K,8)=ALP(8)
            AALP(I,J,K,9)=ALP(9)
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
      DO I = 1,6
        CALL TEST_FILTER(MMIJ(:,:,:,I),FILTER)
      ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
      DO J =1,9
        CALL TEST_FILTER(AALP(:,:,:,J),FILTER)
      ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO&
!$OMP private(B,STR,ALP)&
!$OMP private(SDXF2,SDYF2,SDZF2,AAA,BBB,AMI)
      DO K=1,N3M
        DO J=1,N2M
          DO I=1,N1M
            CALL VELGRAD(I,J,K,ALP,STR)
            SDXF2=(0.5*F2FX(I-1)*(1.-FIXIU(I))+F2FX(I)  &
                  +0.5*F2FX(I+1)*(1.-FIXIL(I)))**2.
            SDYF2=(0.5*F2FY(J-1)*(1.-FIXJU(J))+F2FY(J)  &
                  +0.5*F2FY(J+1)*(1.-FIXJL(J)))**2.
            SDZF2=(0.5*F2FZ(K-1)*(1.-FIXKU(K))+F2FZ(K)  &
                  +0.5*F2FZ(K+1)*(1.-FIXKL(K)))**2.
            B(1)=SDXF2*ALP(1)*ALP(1)&
                +SDYF2*ALP(4)*ALP(4)&
                +SDZF2*ALP(7)*ALP(7)
            B(2)=SDXF2*ALP(1)*ALP(2)&
                +SDYF2*ALP(4)*ALP(5)&
                +SDZF2*ALP(7)*ALP(8)
            B(3)=SDXF2*ALP(1)*ALP(3)&
                +SDYF2*ALP(4)*ALP(6)&
                +SDZF2*ALP(7)*ALP(9)
            B(4)=SDXF2*ALP(2)*ALP(2)&
                +SDYF2*ALP(5)*ALP(5)&
                +SDZF2*ALP(8)*ALP(8)
            B(5)=SDXF2*ALP(2)*ALP(3)&
                +SDYF2*ALP(5)*ALP(6)&
                +SDZF2*ALP(8)*ALP(9)
            B(6)=SDXF2*ALP(3)*ALP(3)&
                +SDYF2*ALP(6)*ALP(6)&
                +SDZF2*ALP(9)*ALP(9)
            AAA=ALP(1)**2.+ALP(2)**2.+ALP(3)**2. &
               +ALP(4)**2.+ALP(5)**2.+ALP(6)**2. &
               +ALP(7)**2.+ALP(8)**2.+ALP(9)**2.
            BBB=B(1)*B(4)+B(4)*B(6)+B(6)*B(1)  &
               -B(2)**2. -B(3)**2. -B(5)**2.
            BBB = DMAX1(BBB,1.0E-16)
            IF (AAA .EQ. 0.) THEN
              AMI = 0D0   
            ELSE
              AMI = SQRT(BBB/AAA)
            ENDIF
            MMIJ(I,J,K,1)=AMI*STR(1)-MMIJ(I,J,K,1)
            MMIJ(I,J,K,2)=AMI*STR(2)-MMIJ(I,J,K,2)
            MMIJ(I,J,K,3)=AMI*STR(3)-MMIJ(I,J,K,3)
            MMIJ(I,J,K,4)=AMI*STR(4)-MMIJ(I,J,K,4)
            MMIJ(I,J,K,5)=AMI*STR(5)-MMIJ(I,J,K,5)
            MMIJ(I,J,K,6)=AMI*STR(6)-MMIJ(I,J,K,6)
            UUI (I,J,K,1)=0.5*(U(I,J,K)+U(I+1,J,K))
            UUI (I,J,K,2)=0.5*(V(I,J,K)+V(I,J+1,K))
            UUI (I,J,K,3)=0.5*(W(I,J,K)+W(I,J,K+1))
            LLIJ(I,J,K,1)=UUI(I,J,K,1)*UUI(I,J,K,1)
            LLIJ(I,J,K,2)=UUI(I,J,K,1)*UUI(I,J,K,2)
            LLIJ(I,J,K,3)=UUI(I,J,K,1)*UUI(I,J,K,3)
            LLIJ(I,J,K,4)=UUI(I,J,K,2)*UUI(I,J,K,2)
            LLIJ(I,J,K,5)=UUI(I,J,K,2)*UUI(I,J,K,3)
            LLIJ(I,J,K,6)=UUI(I,J,K,3)*UUI(I,J,K,3)
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
      DO I = 1,3
        CALL TEST_FILTER(UUI(:,:,:,I),FILTER)
      ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
      DO J =1,6
        CALL TEST_FILTER(LLIJ(:,:,:,J),FILTER)
      ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO&
!$OMP reduction(+:LJMJ_V,MJMJ_V)
      DO K=1,N3M
        DO J=1,N2M
          DO I=1,N1M
            VOLUME = F2FX(I)*F2FY(J)*F2FZ(K)
            LLIJ(I,J,K,1)=LLIJ(I,J,K,1)-UUI(I,J,K,1)*UUI(I,J,K,1)
            LLIJ(I,J,K,2)=LLIJ(I,J,K,2)-UUI(I,J,K,1)*UUI(I,J,K,2)
            LLIJ(I,J,K,3)=LLIJ(I,J,K,3)-UUI(I,J,K,1)*UUI(I,J,K,3)
            LLIJ(I,J,K,4)=LLIJ(I,J,K,4)-UUI(I,J,K,2)*UUI(I,J,K,2)
            LLIJ(I,J,K,5)=LLIJ(I,J,K,5)-UUI(I,J,K,2)*UUI(I,J,K,3)
            LLIJ(I,J,K,6)=LLIJ(I,J,K,6)-UUI(I,J,K,3)*UUI(I,J,K,3)
            LJMJ_V=LJMJ_V                                &
                 +(2.*LLIJ(I,J,K,2)*MMIJ(I,J,K,2)+LLIJ(I,J,K,1)*MMIJ(I,J,K,1) &
                  +2.*LLIJ(I,J,K,3)*MMIJ(I,J,K,3)+LLIJ(I,J,K,4)*MMIJ(I,J,K,4) &
                  +2.*LLIJ(I,J,K,5)*MMIJ(I,J,K,5)+LLIJ(I,J,K,6)*MMIJ(I,J,K,6))&
                 *FLOAT(NWALL_DVM(I,J,K))*VOLUME
            MJMJ_V=MJMJ_V                                &
                   +(2.*MMIJ(I,J,K,2)**2.+MMIJ(I,J,K,1)**2.  &
                    +2.*MMIJ(I,J,K,3)**2.+MMIJ(I,J,K,4)**2.  &
                    +2.*MMIJ(I,J,K,5)**2.+MMIJ(I,J,K,6)**2.) &
                    *FLOAT(NWALL_DVM(I,J,K))*VOLUME
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      DEALLOCATE (AALP)
      DEALLOCATE (LLIJ)
      DEALLOCATE (MMIJ)
      DEALLOCATE (UUI)

      CSGSTS = DMAX1(0., -0.5*LJMJ_V/MJMJ_V)

!$OMP PARALLEL DO
      DO K=1,N3M
        DO J=1,N2M
          DO I=1,N1M
            NUSGS(I,J,K) = CSGSTS * NUSGS(I,J,K)
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE SGS_DVMG
!=======================================================================
!=======================================================================
      SUBROUTINE TEST_FILTER(A0,DIR)
!=======================================================================
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : CFX1,CFY1,CFZ1
      IMPLICIT NONE
      INTEGER*8    :: DIR
      REAL*8       :: A0(N1M,N2M,N3M)

      INTEGER*8   :: I,J,K
      INTEGER*8   :: XFILTER,YFILTER,ZFILTER
      REAL*8      :: TMP(N1M,N2M,N3M)

      IF (DIR .EQ. 1) THEN
        XFILTER = 1
        YFILTER = 0
        ZFILTER = 0
      ELSE IF (DIR .EQ. 2) THEN
        XFILTER = 0
        YFILTER = 0
        ZFILTER = 1
      ELSE IF (DIR .EQ. 3) THEN
        XFILTER = 1
        YFILTER = 0
        ZFILTER = 1
      ELSE
        WRITE(*,*) 'INVALID FILTER INPUT.'
        STOP
      ENDIF 
      
      TMP = A0

      IF (XFILTER .EQ. 1) THEN
        IF (XPRDIC.NE.1) THEN
          DO K=1,N3M
            DO J=1,N2M
              DO I=2,N1M-1
                TMP(I,J,K)=CFX1(I,-1)*TMP(I-1,J,K)  &
                          +CFX1(I, 0)*TMP(I  ,J,K)  &
                          +CFX1(I, 1)*TMP(I+1,J,K)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO J=1,N2M
            DO K=1,N3M
              DO I=2,N1M-1
                TMP(I,J,K)=CFX1(I,-1)*TMP(I-1,J,K)  &
                          +CFX1(I, 0)*TMP(I  ,J,K)  &
                          +CFX1(I, 1)*TMP(I+1,J,K)
              ENDDO
              TMP(1,J,K)=CFX1(I,-1)*TMP(N1M,J,K)  &
                        +CFX1(I, 0)*TMP(1  ,J,K)  &
                        +CFX1(I, 1)*TMP(2  ,J,K)
              TMP(N1M,J,K)=CFX1(I,-1)*TMP(N1M-1,J,K)  &
                          +CFX1(I, 0)*TMP(N1M  ,J,K)  &
                          +CFX1(I, 1)*TMP(1    ,J,K)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      IF (YFILTER .EQ. 1) THEN
        IF (YPRDIC.NE.1) THEN
          DO K=1,N3M
            DO J=2,N2M-1
              DO I=1,N1M
                TMP(I,J,K)=CFY1(J,-1)*TMP(I,J-1,K)  &
                          +CFY1(J, 0)*TMP(I,J  ,K)  &
                          +CFY1(J, 1)*TMP(I,J+1,K)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO K=1,N3M
            DO I=1,N1M
              DO J=2,N2M-1
                TMP(I,J,K)=CFY1(J,-1)*TMP(I,J-1,K)  &
                          +CFY1(J, 0)*TMP(I,J  ,K)  &
                          +CFY1(J, 1)*TMP(I,J+1,K)
              ENDDO
              TMP(I,1,K)=CFY1(J,-1)*TMP(I,N2M,K)  &
                        +CFY1(J, 0)*TMP(I,1  ,K)  &
                        +CFY1(J, 1)*TMP(I,2  ,K)
              TMP(I,N2M,K)=CFY1(J,-1)*TMP(I,N2M-1,K)  &
                          +CFY1(J, 0)*TMP(I,N2M  ,K)  &
                          +CFY1(J, 1)*TMP(I,1    ,K)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      IF (ZFILTER .EQ. 1) THEN
        IF (ZPRDIC.NE.1) THEN
          DO K=2,N3M-1
            DO J=1,N2M
              DO I=1,N1M
                TMP(I,J,K)=CFZ1(K,-1)*TMP(I,J,K-1)  &
                          +CFZ1(K, 0)*TMP(I,J,K  )  &
                          +CFZ1(K, 1)*TMP(I,J,K+1)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO I=1,N1M
            DO J=1,N2M
              DO K=2,N3M-1
                TMP(I,J,K)=CFZ1(K,-1)*TMP(I,J,K-1)  &
                          +CFZ1(K, 0)*TMP(I,J,K  )  &
                          +CFZ1(K, 1)*TMP(I,J,K+1)
              ENDDO
              TMP(I,J,1)=CFZ1(K,-1)*TMP(I,J,N3M)  &
                        +CFZ1(K, 0)*TMP(I,J,1)  &
                        +CFZ1(K, 1)*TMP(I,J,2)
              TMP(I,J,N3M)=CFZ1(K,-1)*TMP(I,J,N3M-1)  &
                          +CFZ1(K, 0)*TMP(I,J,N3M  )  &
                          +CFZ1(K, 1)*TMP(I,J,1    )
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      A0 = TMP

      RETURN
      END SUBROUTINE TEST_FILTER
!=======================================================================
!=======================================================================
      SUBROUTINE VELGRAD(I,J,K,A,SR)
!=======================================================================
!$    use omp_lib
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,NUSGS
      IMPLICIT NONE
      INTEGER*8    :: I,J,K
      REAL*8       :: A(9),SR(6)
      INTEGER*8    :: KPLUS,KMINUS,JPLUS,JMINUS,IPLUS,IMINUS
      REAL*8       :: VG11,VG12,VG13,VG21,VG22,VG23,VG31,VG32,VG33
      REAL*8       :: UP,UM,VP,VM,WP,WM

      KPLUS=KPV(K)
      KMINUS=KMV(K)
      JPLUS=JPV(J)
      JMINUS=JMV(J)
      IPLUS=IPV(I)
      IMINUS=IMV(I)

      VG11=F2FXI(I)*(U(IPLUS,J,K)-U(I,J,K))

      UP=C2CYI(JPLUS)*0.25                                       &
        *(F2FY(JPLUS)*(U(I,J,K)+U(IPLUS,J,K))                    &
         +F2FY(J)*(U(I,JPLUS,K)+U(IPLUS,JPLUS,K)))               &
        *(1.-FIXJU(J))+0.5*(U(I,N2,K)+U(IPLUS,N2,K))*FIXJU(J)
      UM=C2CYI(J)*0.25                                           &
        *(F2FY(J)*(U(I,JMINUS,K)+U(IPLUS,JMINUS,K))              &
          +F2FY(JMINUS)*(U(I,J,K)+U(IPLUS,J,K)))                 &
        *(1.-FIXJL(J))+0.5*(U(I,0,K)+U(IPLUS,0,K))*FIXJL(J)
      VG12=F2FYI(J)*(UP-UM)

      UP=C2CZI(KPLUS)*0.25                                       &
         *(F2FZ(KPLUS)*(U(I,J,K)+U(IPLUS,J,K))                   &
          +F2FZ(K)*(U(I,J,KPLUS)+U(IPLUS,J,KPLUS)))              &
         *(1.-FIXKU(K))+0.5*(U(I,J,N3)+U(IPLUS,J,N3))*FIXKU(K)
      UM=C2CZI(K)*0.25                                           &
         *(F2FZ(K)*(U(I,J,KMINUS)+U(IPLUS,J,KMINUS))             &
          +F2FZ(KMINUS)*(U(I,J,K)+U(IPLUS,J,K)))                 &
         *(1.-FIXKL(K))+0.5*(U(I,J,0)+U(IPLUS,J,0))*FIXKL(K)
      VG13=F2FZI(K)*(UP-UM)

      VP=C2CXI(IPLUS)*0.25                                       &
         *(F2FX(IPLUS)*(V(I,J,K)+V(I,JPLUS,K))                   &
          +F2FX(I)*(V(IPLUS,J,K)+V(IPLUS,JPLUS,K)))              &
         *(1.-FIXIU(I))+0.5*(V(N1,J,K)+V(N1,JPLUS,K))*FIXIU(I)
      VM=C2CXI(I)*0.25                                           &
         *(F2FX(I)*(V(IMINUS,J,K)+V(IMINUS,JPLUS,K))             &
          +F2FX(IMINUS)*(V(I,J,K)+V(I,JPLUS,K)))                 &
         *(1.-FIXIL(I))+0.5*(V(0,J,K)+V(0,JPLUS,K))*FIXIL(I)
      VG21=F2FXI(I)*(VP-VM)

      VG22=F2FYI(J)*(V(I,JPLUS,K)-V(I,J,K))

      VP=C2CZI(KPLUS)*0.25                                       &
         *(F2FZ(KPLUS)*(V(I,J,K)+V(I,JPLUS,K))                   &
          +F2FZ(K)*(V(I,J,KPLUS)+V(I,JPLUS,KPLUS)))              &
         *(1.-FIXKU(K))+0.5*(V(I,J,N3)+V(I,JPLUS,N3))*FIXKU(K)
      VM=C2CZI(K)*0.25                                           &
         *(F2FZ(K)*(V(I,J,KMINUS)+V(I,JPLUS,KMINUS))             &
          +F2FZ(KMINUS)*(V(I,J,K)+V(I,JPLUS,K)))                 &
         *(1.-FIXKL(K))+0.5*(V(I,J,0)+V(I,JPLUS,0))*FIXKL(K)
      VG23=F2FZI(K)*(VP-VM)
      
      WP=C2CXI(IPLUS)*0.25                                       &
         *(F2FX(IPLUS)*(W(I,J,K)+W(I,J,KPLUS))                   &
          +F2FX(I)*(W(IPLUS,J,K)+W(IPLUS,J,KPLUS)))              &
         *(1.-FIXIU(I))+0.5*(W(N1,J,K)+W(N1,J,KPLUS))*FIXIU(I)
      WM=C2CXI(I)*0.25                                           &
         *(F2FX(I)*(W(IMINUS,J,K)+W(IMINUS,J,KPLUS))             &
          +F2FX(IMINUS)*(W(I,J,K)+W(I,J,KPLUS)))                 &
         *(1.-FIXIL(I))+0.5*(W(0,J,K)+W(0,J,KPLUS))*FIXIL(I)
      VG31=F2FXI(I)*(WP-WM)
      
      WP=C2CYI(JPLUS)*0.25                                       &
         *(F2FY(JPLUS)*(W(I,J,K)+W(I,J,KPLUS))                   &
          +F2FY(J)*(W(I,JPLUS,K)+W(I,JPLUS,KPLUS)))              &
         *(1.-FIXJU(J))+0.5*(W(I,N2,K)+W(I,N2,KPLUS))*FIXJU(J)
      WM=C2CYI(J)*0.25                                           &
         *(F2FY(J)*(W(I,JMINUS,K)+W(I,JMINUS,KPLUS))             &
          +F2FY(JMINUS)*(W(I,J,K)+W(I,J,KPLUS)))                 &
          *(1.-FIXJL(J))+0.5*(W(I,0,K)+W(I,0,KPLUS))*FIXJL(J)
      VG32=F2FYI(J)*(WP-WM)
      
      VG33=F2FZI(K)*(W(I,J,KPLUS)-W(I,J,K))
      
      A(1)=VG11
      A(2)=VG21
      A(3)=VG31
      A(4)=VG12
      A(5)=VG22
      A(6)=VG32
      A(7)=VG13
      A(8)=VG23
      A(9)=VG33

      SR(1)= VG11
      SR(2)= 0.5*(VG12+VG21)
      SR(3)= 0.5*(VG13+VG31)
      SR(4)= VG22
      SR(5)= 0.5*(VG23+VG32)
      SR(6)= VG33

      RETURN
      END SUBROUTINE VELGRAD
!=======================================================================
!=======================================================================
      SUBROUTINE NUTZERO
!=======================================================================
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : NUSGS,INZ,JNZ,KNZ
      IMPLICIT NONE
      INTEGER*8  :: I,J,K
      INTEGER*8  :: IMAX,JMAX,KMAX,L
      REAL*8     :: NUSGSAVG, NUSGSMAX, AREA
      REAL*8     :: FUNCBODY

      IMAX = 0
      JMAX = 0
      KMAX = 0

!$OMP PARALLEL DO
        DO L=1,NZERO
          NUSGS(INZ(L),JNZ(L),KNZ(L)) = 0.
        ENDDO
!$OMP END PARALLEL DO

      NUSGSAVG = 0.
      AREA = 0.
      NUSGSMAX = 0.

!$OMP PARALLEL DO&
!$OMP reduction(+:NUSGSAVG, AREA)&
!$OMP reduction(MAX:NUSGSMAX)
      DO K=1,N3M
        DO J=1,N2M
          DO I=1,N1M
            IF (FUNCBODY(XMP(I),YMP(J),ZMP(K),TIME).GE.1.E-10) THEN
              NUSGSAVG = NUSGSAVG + NUSGS(I,J,K)*F2FX(I)*F2FY(J)*F2FZ(K)
              AREA = AREA + F2FX(I)*F2FY(J)*F2FZ(K)
              IF (NUSGS(I,J,K).GT.NUSGSMAX) THEN
                NUSGSMAX = NUSGS(I,J,K)
                IMAX = I
                JMAX = J
                KMAX = K
              ENDIF
            ELSE
              NUSGS(I,J,K) = 0.
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      NUSGSAVG = NUSGSAVG / AREA

      NUTAVG = NUSGSAVG * RE
      NUTMAX = NUSGSMAX * RE

      WRITE(*,99)  CSGSTS
      WRITE(*,100) NUTAVG
      WRITE(*,110) NUTMAX,XMP(IMAX),YMP(JMAX),ZMP(KMAX),IMAX,JMAX,KMAX
   99 FORMAT('CSGSTS =  ',ES10.3)
  100 FORMAT('NUTAVG =  ',ES10.3)
  110 FORMAT('NUTMAX =  ',ES10.3,' @ ',3F10.4,' , ',3I5)

      RETURN
      END SUBROUTINE NUTZERO
!=======================================================================
!=======================================================================
      SUBROUTINE NUTINTERPOL
!=======================================================================
!
!     INTERPOLATION OF NUT
!     MAKE FIELD OF NUSGS1(X,Y,Z)
!
!     AT THE NO SLIP WALL                    : NUT=0
!     VELOCITY INLET                         : D(NUT)/DX=0
!     AT THE CONVECTIVE WALL(OUTLET WALL)    : D(NUT)/DX=0
!     FAR FIELD WALL(SIED WALL)              : D(NUT)/DX=0
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : NUSGS,NUSGS1,INZ,JNZ,KNZ
      IMPLICIT NONE
      INTEGER*8  :: I,J,K,IM,JM,KM

      NUSGS1 = 0.

      IF (XPRDIC .EQ. 1) THEN
!$OMP PARALLEL DO
        DO K=0,N3
          DO J=0,N2
            NUSGS(0 ,J,K)=NUSGS(N1M,J,K)
            NUSGS(N1,J,K)=NUSGS(1  ,J,K)
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDIF

      IF (YPRDIC .EQ. 1) THEN
!$OMP PARALLEL DO
        DO K=0,N3
          DO I=0,N1
            NUSGS(I,0 ,K)=NUSGS(I,N2M,K)
            NUSGS(I,N2,K)=NUSGS(I,  1,K)
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDIF

      IF (ZPRDIC .EQ. 1) THEN
!$OMP PARALLEL DO
        DO J=0,N2
          DO I=0,N1
            NUSGS(I,J,0 )=NUSGS(I,J,N3M)
            NUSGS(I,J,N3)=NUSGS(I,J,1  )
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDIF

!$OMP PARALLEL DO
      DO K=K_BGPZ,N3M
        DO J=J_BGPY,N2M
          DO I=1,N1M
            NUSGS1(I,J,K,1)=C2CYI(J)*C2CZI(K)*0.25                          &
              *(F2FY(J-1)*(F2FZ(K-1)*NUSGS(I,J,K)+F2FZ(K)*NUSGS(I,J,K-1))   &
              +F2FY(J)*(F2FZ(K-1)*NUSGS(I,J-1,K)+F2FZ(K)*NUSGS(I,J-1,K-1)))
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
      DO K=K_BGPZ,N3M
        DO J=1,N2M
          DO I=I_BGPX,N1M
            NUSGS1(I,J,K,2)=C2CZI(K)*C2CXI(I)*0.25                          &
              *(F2FZ(K-1)*(F2FX(I-1)*NUSGS(I,J,K)+F2FX(I)*NUSGS(I-1,J,K))   &
              +F2FZ(K)*(F2FX(I-1)*NUSGS(I,J,K-1)+F2FX(I)*NUSGS(I-1,J,K-1)))
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
      DO K=1,N3M
        DO J=J_BGPY,N2M
          DO I=I_BGPX,N1M
            NUSGS1(I,J,K,3)=C2CXI(I)*C2CYI(J)*0.25                          &
              *(F2FX(I-1)*(F2FY(J-1)*NUSGS(I,J,K)+F2FY(J)*NUSGS(I,J-1,K))   &
              +F2FX(I)*(F2FY(J-1)*NUSGS(I-1,J,K)+F2FY(J)*NUSGS(I-1,J-1,K)))
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

!!$OMP PARALLEL DO private(JM,KM)
      DO J=J_BGPY,N2M
        DO K=K_BGPZ,N3M
          JM=JMV(J)
          KM=KMV(K)
          IF (XPRDIC .EQ. 1) THEN
            NUSGS1(N1,J,K,1) = NUSGS1(1,J,K,1)
            NUSGS1(N1,J,K,2) = NUSGS1(1,J,K,2)
            NUSGS1(N1,J,K,3) = NUSGS1(1,J,K,3)
          ELSE
            IF (BC_XBTM .EQ. 0) THEN
              NUSGS1(1,J,K,1) = 0.
              NUSGS1(1,J,K,2) = 0.
              NUSGS1(1,J,K,3) = 0.
            ELSE
              NUSGS1(1,J,K,2) = &
                C2CZI(K)*0.5*(F2FZ(K)*NUSGS(2,J,KM)+F2FZ(KM)*NUSGS(2,J,K))
              NUSGS1(1,J,K,3) = &
                C2CYI(J)*0.5*(F2FY(J)*NUSGS(2,JM,K)+F2FY(JM)*NUSGS(2,J,K))
            ENDIF
            IF (BC_XTOP .EQ. 0) THEN
              NUSGS1(N1,J,K,1) = 0.
              NUSGS1(N1,J,K,2) = 0.
              NUSGS1(N1,J,K,3) = 0.
            ELSE
              NUSGS1(N1,J,K,2) = &
                C2CZI(K)*0.5*(F2FZ(K)*NUSGS(N1M,J,KM)+F2FZ(KM)*NUSGS(N1M,J,K))
              NUSGS1(N1,J,K,3) = &
                C2CYI(J)*0.5*(F2FY(J)*NUSGS(N1M,JM,K)+F2FY(JM)*NUSGS(N1M,J,K))
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!!$OMP END PARALLEL DO

!!$OMP PARALLEL DO private(IM,KM)
      DO K=K_BGPZ,N3M
        DO I=I_BGPX,N1M
          IM=IMV(I)
          KM=KMV(K)
          IF (YPRDIC .EQ. 1) THEN
            NUSGS1(I,N2,K,1) = NUSGS1(I,1,K,1)
            NUSGS1(I,N2,K,2) = NUSGS1(I,1,K,2)
            NUSGS1(I,N2,K,3) = NUSGS1(I,1,K,3)
          ELSE
            IF (BC_YBTM .EQ. 0) THEN
              NUSGS1(I,1,K,1) = 0.
              NUSGS1(I,1,K,2) = 0.
              NUSGS1(I,1,K,3) = 0.
            ELSE
              NUSGS1(I,1,K,2) = &
                C2CZI(K)*0.5*(F2FZ(K)*NUSGS(I,2,KM)+F2FZ(KM)*NUSGS(I,2,K))
              NUSGS1(I,1,K,3) = &
                C2CXI(I)*0.5*(F2FX(I)*NUSGS(IM,2,K)+F2FX(IM)*NUSGS(I,2,K))
            ENDIF
            IF (BC_YTOP .EQ. 0) THEN
              NUSGS1(I,N2,K,1) = 0.
              NUSGS1(I,N2,K,2) = 0.
              NUSGS1(I,N2,K,3) = 0.
            ELSE
              NUSGS1(I,N2,K,2) = &
                C2CZI(K)*0.5*(F2FZ(K)*NUSGS(I,N2M,KM)+F2FZ(KM)*NUSGS(I,N2M,K))
              NUSGS1(I,N2,K,3) = &
                C2CXI(I)*0.5*(F2FX(I)*NUSGS(IM,N2M,K)+F2FX(IM)*NUSGS(I,N2M,K))
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!!$OMP END PARALLEL DO

!!$OMP PARALLEL DO private(IM,JM)
      DO I=I_BGPX,N1M
        DO J=J_BGPY,N2M
          IM=IMV(I)
          JM=JMV(J)
          IF (ZPRDIC .EQ. 1) THEN
            NUSGS1(I,J,N3,1) = NUSGS1(I,J,1,1)
            NUSGS1(I,J,N3,2) = NUSGS1(I,J,1,2)
            NUSGS1(I,J,N3,3) = NUSGS1(I,J,1,3)
          ELSE
            IF (BC_ZBTM .EQ. 0) THEN
              NUSGS1(I,J,1,1) = 0.
              NUSGS1(I,J,1,2) = 0.
              NUSGS1(I,J,1,3) = 0.
            ELSE
              NUSGS1(I,J,1,2) = &
                C2CYI(J)*0.5*(F2FY(J)*NUSGS(I,JM,2)+F2FY(JM)*NUSGS(I,J,2))
              NUSGS1(I,J,1,3) = &
                C2CXI(I)*0.5*(F2FX(I)*NUSGS(IM,J,2)+F2FX(IM)*NUSGS(I,J,2))
            ENDIF
            IF (BC_ZTOP .EQ. 0) THEN
              NUSGS1(I,J,N3,1) = 0.
              NUSGS1(I,J,N3,2) = 0.
              NUSGS1(I,J,N3,3) = 0.
            ELSE
              NUSGS1(I,J,N3,2) = &
                C2CYI(K)*0.5*(F2FY(J)*NUSGS(I,JM,N3M)+F2FY(JM)*NUSGS(I,J,N3M))
              NUSGS1(I,J,N3,3) = &
                C2CXI(I)*0.5*(F2FX(I)*NUSGS(IM,J,N3M)+F2FX(IM)*NUSGS(I,J,N3M))
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!!$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE NUTINTERPOL
!=======================================================================
!=======================================================================
      SUBROUTINE RHSSGS
!=======================================================================
!
!     D/DXj(NU_t*DUj/DXi)*(1-KroneckerDelta(ij)),
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,NUSGS1,RHS1
      IMPLICIT NONE
      INTEGER*8 :: I,J,K
      REAL*8    :: RHSU2,RHSU3,RHSV1,RHSV3,RHSW1,RHSW2

!$OMP PARALLEL DO&
!$OMP private(RHSU2,RHSU3)
      DO K=1,N3M
        DO J=1,N2M
          DO I=I_BGPX,N1M
            RHSU2=F2FYI(J)*C2CXI(I)                            &
                 *(NUSGS1(I,J+1,K,3)*(V(I,J+1,K)-V(I-1,J+1,K)) &
                  -NUSGS1(I,J,K,3)*(V(I,J,K)-V(I-1,J,K)))
            RHSU3=F2FZI(K)*C2CXI(I)                            &
                 *(NUSGS1(I,J,K+1,2)*(W(I,J,K+1)-W(I-1,J,K+1)) &
                  -NUSGS1(I,J,K,2)*(W(I,J,K)-W(I-1,J,K)))
            RHS1(I,J,K,1)=RHSU2+RHSU3
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO&
!$OMP private(RHSV1,RHSV3)
      DO K=1,N3M
        DO J=J_BGPY,N2M
          DO I=1,N1M
            RHSV1=F2FXI(I)*C2CYI(J)                              &
                  *(NUSGS1(I+1,J,K,3)*(U(I+1,J,K)-U(I+1,J-1,K))  &
                   -NUSGS1(I,J,K,3)*(U(I,J,K)-U(I,J-1,K)))
            RHSV3=F2FZI(K)*C2CYI(J)                              &
                  *(NUSGS1(I,J,K+1,1)*(W(I,J,K+1)-W(I,J-1,K+1))  &
                   -NUSGS1(I,J,K,1)*(W(I,J,K)-W(I,J-1,K)))
            RHS1(I,J,K,2)=RHSV1+RHSV3
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO&
!$OMP private(RHSW1,RHSW2)
      DO K=K_BGPZ,N3M
        DO J=1,N2M
          DO I=1,N1M
            RHSW1=F2FXI(I)*C2CZI(K)                              &
                   *(NUSGS1(I+1,J,K,2)*(U(I+1,J,K)-U(I+1,J,K-1)) &
                   -NUSGS1(I,J,K,2)*(U(I,J,K)-U(I,J,K-1)))
            RHSW2=F2FYI(J)*C2CZI(K)                              &
                   *(NUSGS1(I,J+1,K,1)*(V(I,J+1,K)-V(I,J+1,K-1)) &
                   -NUSGS1(I,J,K,1)*(V(I,J,K)-V(I,J,K-1)))
            RHS1(I,J,K,3)=RHSW1+RHSW2
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE RHSSGS
!=======================================================================
!=======================================================================
      SUBROUTINE SGSCALC_T
!=======================================================================
!
! CALCULATE SUBGRID-SCALE EDDY DIFFUSIVITY, ALSGS
!
! ITEMDL = 0 -> EDDY DIFFUSIVITY MODEL (MOIN ET AL.,  1991. PHYS.FLUIDS)
! 
! ITEMDL = 1 -> NONE
!
! IDVMON = OFF, CONSTANT COEFF. MODEL
! IDVMON = ON , DYNAMIC COEFF. PROCEDURE
! (DSLM MODEL: LILLY,       1991. PHYS.FLUIDS A)
! (DVMG MODEL: PARK ET AL., 2006. PHYS.FLUIDS)
!
! PLZ REFER TO "GERMANO IDENTITY" THEORY (GERMANO, 1992. J.FLUID MECH.)
!
!-----------------------------------------------------------------------
!$    use omp_lib
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : ALSGS
      IMPLICIT NONE
      IF (MSUB .EQ. 1) THEN
        
        IF (IREAD.EQ.0 .AND. NTIME.EQ.1) THEN
          ALSGS = 0D0
        ELSE
          IF (ITEMDL .EQ. 0) THEN
            IF (IDVMON .EQ. 1) THEN
             CALL SGS_EDM
            ELSE
             CALL SGS_CEDM
            ENDIF
          ELSE IF (ITEMDL .EQ. 1) THEN
            WRITE(*,*) 'OTHER MODEL IS NOT IMPLEMENTED IN THIS VERSION'
            WRITE(*,*) 'SGS.F90 LINE ***'            
            STOP
          ENDIF
        ENDIF

        CALL ALPZERO               ! SET ALSGS TO ZERO IN THE SOLID BODY
        CALL ALPINTERPOL           !              INTERPOLATION OF NUSGS
      ENDIF

      RETURN
      END SUBROUTINE SGSCALC_T
!=======================================================================
!=======================================================================
      SUBROUTINE SGS_CEDM
!=======================================================================
!
!     ALSGS=CSGSHF*SQRT(BBB/AAA)
!     CSGSHF : EDDY DIFFUSIVITY CONSTANT (= CSGSHF FROM MOD_COMMON)
!     DEL    : FILTER WIDTH CHARACTERIZED BY MEAN CELL WIDTH
!     SSS    : S_{ij}*S_{ij}
!     S_{ij} : STRAIN RATE (0.5 * [DUj/DXi + DUi/DXj]) 
!
!-----------------------------------------------------------------------  
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,T,ALSGS
      IMPLICIT NONE
      INTEGER*8 :: I,J,K,KPLUS,KMIALS,JPLUS,JMIALS,IPLUS,IMIALS
      REAL*8    :: VG11,VG12,VG13,VG21,VG22,VG23,VG31,VG32,VG33
      REAL*8    :: UP,UM,VP,VM,WP,WM
      REAL*8    :: D(9),S(6)
      REAL*8    :: DEL,SSS

!$OMP PARALLEL DO private(UP,UM,VP,VM,WP,WM)&
!$OMP private(D,S,DEL,SSS)
      DO K=1,N3M
        DO J=1,N2M
          DO I=1,N1M

            CALL VELGRAD(I,J,K,D,S)         

            D(1) = F2FX(I)
            D(2) = F2FY(J)
            D(3) = F2FZ(K)
      
            DEL=(D(1)*D(2)*D(3))**(1./3.)
            SSS=S(1)**2.+2*S(2)**2.      &
               +2*S(3)**2.+S(4)**2.      &
               +2*S(5)**2.+S(6)**2.
      
            SSS = DMAX1(SSS, 1.0E-16)
   
            ALSGS(I,J,K)= CSGSHF*DEL**2.*SQRT(2.*SSS)
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE SGS_CEDM
!=======================================================================
!=======================================================================
      SUBROUTINE SGS_EDM
!=======================================================================
!
!     ALSGS=CSGSHF*DEL**2.*SQRT(2.*SSS)
!           CSGSHF DETERMINED VIA GERMANO IDENTITY
!     CSGSHF : EDDY DIFFUSIVITY CONSTANT (= CSGSHF FROM MOD_COMMON)
!     DEL    : FILTER WIDTH CHARACTERIZED BY MEAN CELL WIDTH
!     SSS    : S_{ij}*S_{ij}
!     S_{ij} : STRAIN RATE (0.5 * [DUj/DXi + DUi/DXj]) 
!
!-----------------------------------------------------------------------  
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,T,ALSGS,NWALL_DVM,AALP,LLIJ,MMIJ,UUI
      IMPLICIT NONE
      INTEGER*8     :: I,J,K
      REAL*8  :: LJMJ_V,MJMJ_V
      REAL*8  :: DEL,SSS,AMI,VOLUME
      REAL*8  :: D(9),S(6),ALP(3)

      LJMJ_V = 0.
      MJMJ_V = 0.
      CSGSHF   = 0.

      ALLOCATE (AALP(N1M,N2M,N3M,6))
      ALLOCATE (LLIJ(N1M,N2M,N3M,3))
      ALLOCATE (MMIJ(N1M,N2M,N3M,3))
      ALLOCATE (UUI (N1M,N2M,N3M,3))

!$OMP PARALLEL DO&
!$OMP private(D,S,ALP)&
!$OMP private(DEL,SSS)
      DO K=1,N3M
        DO J=1,N2M
          DO I=1,N1M
            CALL VELGRAD(I,J,K,D,S) ! D IS DUMMY VARIABLE
            CALL TEMGRAD(I,J,K,ALP)
            
            D(1) = F2FX(I)
            D(2) = F2FY(J)
            D(3) = F2FZ(K)
      
            DEL=(D(1)*D(2)*D(3))**(1./3.)
            SSS=S(1)**2.+2*S(2)**2.      &
               +2*S(3)**2.+S(4)**2.      &
               +2*S(5)**2.+S(6)**2.
      
            SSS = DMAX1(SSS, 1.0E-16)
   
            ALSGS(I,J,K)= DEL**2.*SQRT(2.*SSS)
            
            MMIJ(I,J,K,1) = ALSGS(I,J,K)*ALP(1)
            MMIJ(I,J,K,2) = ALSGS(I,J,K)*ALP(2)
            MMIJ(I,J,K,3) = ALSGS(I,J,K)*ALP(3)
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
      DO I = 1,3
        CALL TEST_FILTER(MMIJ(:,:,:,I),FILTER)
      ENDDO
!$OMP END PARALLEL DO        

!$OMP PARALLEL DO&
!$OMP private(D,S,ALP)&
!$OMP private(DEL,SSS)
      DO K=1,N3M
        DO J=1,N2M
          DO I=1,N1M
            CALL VELGRAD(I,J,K,D,S) ! D IS DUMMY VARIABLE
            CALL TEMGRAD(I,J,K,ALP)
            
            D(1) = (0.5*F2FX(I-1)*(1.-FIXIU(I))+F2FX(I)  &
                  +0.5*F2FX(I+1)*(1.-FIXIL(I)))
            D(2) = (0.5*F2FY(J-1)*(1.-FIXJU(J))+F2FY(J)  &
                  +0.5*F2FY(J+1)*(1.-FIXJL(J)))
            D(3) = (0.5*F2FZ(K-1)*(1.-FIXKU(K))+F2FZ(K)  &
                  +0.5*F2FZ(K+1)*(1.-FIXKL(K)))
      
            DEL=(D(1)*D(2)*D(3))**(1./3.)
            SSS=S(1)**2.+2*S(2)**2.      &
               +2*S(3)**2.+S(4)**2.      &
               +2*S(5)**2.+S(6)**2.
      
            SSS = DMAX1(SSS,1.0E-16)
   
            AMI = DEL**2.*SQRT(2.* SSS)

            MMIJ(I,J,K,1) = AMI*ALP(1) - MMIJ(I,J,K,1)
            MMIJ(I,J,K,2) = AMI*ALP(2) - MMIJ(I,J,K,2)
            MMIJ(I,J,K,3) = AMI*ALP(3) - MMIJ(I,J,K,3)
            UUI (I,J,K,1)=0.5*(U(I,J,K)+U(I+1,J,K))
            UUI (I,J,K,2)=0.5*(V(I,J,K)+V(I,J+1,K))
            UUI (I,J,K,3)=0.5*(W(I,J,K)+W(I,J,K+1))
            AALP(I,J,K,1)=T(I,J,K)
            AALP(I,J,K,2)=T(I,J,K)
            AALP(I,J,K,3)=T(I,J,K)
            LLIJ(I,J,K,1)=UUI(I,J,K,1)*AALP(I,J,K,1)
            LLIJ(I,J,K,2)=UUI(I,J,K,2)*AALP(I,J,K,2)
            LLIJ(I,J,K,3)=UUI(I,J,K,3)*AALP(I,J,K,3)
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
      DO I = 1,3
        CALL TEST_FILTER(UUI(:,:,:,I),FILTER)
      ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
      DO I = 1,3
        CALL TEST_FILTER(AALP(:,:,:,I),FILTER)
      ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
      DO J =1,3
        CALL TEST_FILTER(LLIJ(:,:,:,J),FILTER)
      ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO&
!$OMP reduction(+:LJMJ_V,MJMJ_V)
      DO K=1,N3M
        DO J=1,N2M
          DO I=1,N1M
            VOLUME = F2FX(I)*F2FY(J)*F2FZ(K)
            LLIJ(I,J,K,1)=LLIJ(I,J,K,1)-UUI(I,J,K,1)*AALP(I,J,K,1)
            LLIJ(I,J,K,2)=LLIJ(I,J,K,2)-UUI(I,J,K,2)*AALP(I,J,K,2)
            LLIJ(I,J,K,3)=LLIJ(I,J,K,3)-UUI(I,J,K,3)*AALP(I,J,K,3)
            LJMJ_V=LJMJ_V                                &
                 +(LLIJ(I,J,K,1)*MMIJ(I,J,K,1) &
                  +LLIJ(I,J,K,2)*MMIJ(I,J,K,2) &
                  +LLIJ(I,J,K,3)*MMIJ(I,J,K,3))&
                 *FLOAT(NWALL_DVM(I,J,K))*VOLUME
            MJMJ_V=MJMJ_V                                &
                   +(MMIJ(I,J,K,1)*MMIJ(I,J,K,1) &
                    +MMIJ(I,J,K,2)*MMIJ(I,J,K,2) &
                    +MMIJ(I,J,K,3)*MMIJ(I,J,K,3))&
                   *FLOAT(NWALL_DVM(I,J,K))*VOLUME
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      DEALLOCATE (AALP)
      DEALLOCATE (LLIJ)
      DEALLOCATE (MMIJ)
      DEALLOCATE (UUI)

      CSGSHF = DMAX1(0., -LJMJ_V/MJMJ_V)

!$OMP PARALLEL DO
      DO K=1,N3M
        DO J=1,N2M
          DO I=1,N1M
            ALSGS(I,J,K) = CSGSHF * ALSGS(I,J,K)
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE SGS_EDM
!=======================================================================
      SUBROUTINE TEMGRAD(I,J,K,A)
!=======================================================================
!$    use omp_lib
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,T,ALSGS
      IMPLICIT NONE
      INTEGER*8    :: I,J,K
      REAL*8       :: A(3)
      INTEGER*8    :: KPLUS,KMINUS,JPLUS,JMINUS,IPLUS,IMINUS
      REAL*8       :: TG1,TG2,TG3
      REAL*8       :: TP,TM

      KPLUS=KPV(K)
      KMINUS=KMV(K)
      JPLUS=JPV(J)
      JMINUS=JMV(J)
      IPLUS=IPV(I)
      IMINUS=IMV(I)

      TP=0.5*(F2FX(I)*T(IPLUS,J,K)+F2FX(IPLUS)*T(I,J,K))*C2CXI(IPLUS) &
         *(1.-FIXIU(I))+T(IPLUS,J,K)*FIXIU(I)
      TM=0.5*(F2FX(IMINUS)*T(I,J,K)+F2FX(I)*T(IMINUS,J,K))*C2CXI(I)   &
         *(1.-FIXIL(I))+T(IMINUS,J,K)*FIXIL(I)

      TG1 = F2FXI(I)*(TP-TM)

      TP=0.5*(F2FY(J)*T(I,JPLUS,K)+F2FY(JPLUS)*T(I,J,K))*C2CYI(JPLUS) &
         *(1.-FIXJU(J))+T(I,JPLUS,K)*FIXJU(J)
      TM=0.5*(F2FY(JMINUS)*T(I,J,K)+F2FY(J)*T(I,JMINUS,K))*C2CYI(J)   &
         *(1.-FIXJL(J))+T(I,JMINUS,K)*FIXJL(J)

      TG2 = F2FYI(I)*(TP-TM)

      TP=0.5*(F2FZ(K)*T(I,J,KPLUS)+F2FZ(KPLUS)*T(I,J,K))*C2CZI(KPLUS) &
         *(1.-FIXKU(K))+T(I,J,KPLUS)*FIXKU(K)
      TM=0.5*(F2FZ(KMINUS)*T(I,J,K)+F2FZ(K)*T(I,J,KMINUS))*C2CZI(K)   &
         *(1.-FIXKL(K))+T(I,J,KMINUS)*FIXKL(K)

      TG3 = F2FZI(I)*(TP-TM)

      A(1)=TG1
      A(2)=TG2
      A(3)=TG3

      RETURN
      END SUBROUTINE TEMGRAD
!=======================================================================
!=======================================================================
      SUBROUTINE ALPZERO
!=======================================================================
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : ALSGS,INZ,JNZ,KNZ
      IMPLICIT NONE
      INTEGER*8  :: I,J,K
      INTEGER*8  :: IMAX,JMAX,KMAX,L
      REAL*8     :: ALSGSAVG, ALSGSMAX, AREA
      REAL*8     :: FUNCBODY

      IMAX = 0
      JMAX = 0
      KMAX = 0

!$OMP PARALLEL DO
        DO L=1,NZERO
          ALSGS(INZ(L),JNZ(L),KNZ(L)) = 0.
        ENDDO
!$OMP END PARALLEL DO

      ALSGSAVG = 0.
      AREA = 0.
      ALSGSMAX = 0.

!$OMP PARALLEL DO&
!$OMP reduction(+:ALSGSAVG, AREA)&
!$OMP reduction(MAX:ALSGSMAX)
      DO K=1,N3M
        DO J=1,N2M
          DO I=1,N1M
            IF (FUNCBODY(XMP(I),YMP(J),ZMP(K),TIME).GE.1.E-10) THEN
              ALSGSAVG = ALSGSAVG + ALSGS(I,J,K)*F2FX(I)*F2FY(J)*F2FZ(K)
              AREA = AREA + F2FX(I)*F2FY(J)*F2FZ(K)
              IF (ALSGS(I,J,K).GT.ALSGSMAX) THEN
                ALSGSMAX = ALSGS(I,J,K)
                IMAX = I
                JMAX = J
                KMAX = K
              ENDIF
            ELSE
              ALSGS(I,J,K) = 0.
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ALSGSAVG = ALSGSAVG / AREA

      ALPAVG = ALSGSAVG * RE * PR
      ALPMAX = ALSGSMAX * RE * PR

      WRITE(*,99)  CSGSHF
      WRITE(*,100) ALPAVG
      WRITE(*,110) ALPMAX,XMP(IMAX),YMP(JMAX),ZMP(KMAX),IMAX,JMAX,KMAX
   99 FORMAT('CSGSHF =  ',ES10.3)
  100 FORMAT('ALPAVG =  ',ES10.3)
  110 FORMAT('ALPMAX =  ',ES10.3,' @ ',3F10.4,' , ',3I5)

      RETURN
      END SUBROUTINE ALPZERO
!=======================================================================
!=======================================================================
      SUBROUTINE ALPINTERPOL
!=======================================================================
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : ALSGS,ALSGS1,INZ,JNZ,KNZ
      IMPLICIT NONE
      INTEGER*8  :: I,J,K,IM,JM,KM

      ALSGS1 = 0.

      IF (XPRDIC .EQ. 1) THEN
!$OMP PARALLEL DO
        DO K=0,N3
          DO J=0,N2
            ALSGS(0 ,J,K)=ALSGS(N1M,J,K)
            ALSGS(N1,J,K)=ALSGS(1  ,J,K)
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDIF

      IF (YPRDIC .EQ. 1) THEN
!$OMP PARALLEL DO
        DO K=0,N3
          DO I=0,N1
            ALSGS(I,0 ,K)=ALSGS(I,N2M,K)
            ALSGS(I,N2,K)=ALSGS(I,  1,K)
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDIF

      IF (ZPRDIC .EQ. 1) THEN
!$OMP PARALLEL DO
        DO J=0,N2
          DO I=0,N1
            ALSGS(I,J,0 )=ALSGS(I,J,N3M)
            ALSGS(I,J,N3)=ALSGS(I,J,1  )
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDIF

!$OMP PARALLEL DO
      DO K=1,N3M
        DO J=1,N2M
          DO I=1,N1M
            ALSGS1(I,J,K,1)=C2CYI(J)*C2CZI(K)*0.25                          &
              *(F2FY(J-1)*(F2FZ(K-1)*ALSGS(I,J,K)+F2FZ(K)*ALSGS(I,J,K-1))   &
              +F2FY(J)*(F2FZ(K-1)*ALSGS(I,J-1,K)+F2FZ(K)*ALSGS(I,J-1,K-1)))
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
      DO K=1,N3M
        DO J=1,N2M
          DO I=1,N1M
            ALSGS1(I,J,K,2)=C2CZI(K)*C2CXI(I)*0.25                          &
              *(F2FZ(K-1)*(F2FX(I-1)*ALSGS(I,J,K)+F2FX(I)*ALSGS(I-1,J,K))   &
              +F2FZ(K)*(F2FX(I-1)*ALSGS(I,J,K-1)+F2FX(I)*ALSGS(I-1,J,K-1)))
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
      DO K=1,N3M
        DO J=1,N2M
          DO I=1,N1M
            ALSGS1(I,J,K,3)=C2CXI(I)*C2CYI(J)*0.25                          &
              *(F2FX(I-1)*(F2FY(J-1)*ALSGS(I,J,K)+F2FY(J)*ALSGS(I,J-1,K))   &
              +F2FX(I)*(F2FY(J-1)*ALSGS(I-1,J,K)+F2FY(J)*ALSGS(I-1,J-1,K)))
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

!!$OMP PARALLEL DO private(JM,KM)
      DO J=1,N2M
        DO K=1,N3M
          JM=JMV(J)
          KM=KMV(K)
          IF (XPRDIC .EQ. 1) THEN
            ALSGS1(N1,J,K,1) = ALSGS1(1,J,K,1)
            ALSGS1(N1,J,K,2) = ALSGS1(1,J,K,2)
            ALSGS1(N1,J,K,3) = ALSGS1(1,J,K,3)
          ELSE
            IF (BC_XBTM .EQ. 0) THEN
              ALSGS1(1,J,K,1) = 0.
              ALSGS1(1,J,K,2) = 0.
              ALSGS1(1,J,K,3) = 0.
            ELSE
              ALSGS1(1,J,K,2) = &
                C2CZI(K)*0.5*(F2FZ(K)*ALSGS(2,J,KM)+F2FZ(KM)*ALSGS(2,J,K))
              ALSGS1(1,J,K,3) = &
                C2CYI(J)*0.5*(F2FY(J)*ALSGS(2,JM,K)+F2FY(JM)*ALSGS(2,J,K))
            ENDIF
            IF (BC_XTOP .EQ. 0) THEN
              ALSGS1(N1,J,K,1) = 0.
              ALSGS1(N1,J,K,2) = 0.
              ALSGS1(N1,J,K,3) = 0.
            ELSE
              ALSGS1(N1,J,K,2) = &
                C2CZI(K)*0.5*(F2FZ(K)*ALSGS(N1M,J,KM)+F2FZ(KM)*ALSGS(N1M,J,K))
              ALSGS1(N1,J,K,3) = &
                C2CYI(J)*0.5*(F2FY(J)*ALSGS(N1M,JM,K)+F2FY(JM)*ALSGS(N1M,J,K))
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!!$OMP END PARALLEL DO

!!$OMP PARALLEL DO private(IM,KM)
      DO K=1,N3M
        DO I=1,N1M
          IM=IMV(I)
          KM=KMV(K)
          IF (YPRDIC .EQ. 1) THEN
            ALSGS1(I,N2,K,1) = ALSGS1(I,1,K,1)
            ALSGS1(I,N2,K,2) = ALSGS1(I,1,K,2)
            ALSGS1(I,N2,K,3) = ALSGS1(I,1,K,3)
          ELSE
            IF (BC_YBTM .EQ. 0) THEN
              ALSGS1(I,1,K,1) = 0.
              ALSGS1(I,1,K,2) = 0.
              ALSGS1(I,1,K,3) = 0.
            ELSE
              ALSGS1(I,1,K,2) = &
                C2CZI(K)*0.5*(F2FZ(K)*ALSGS(I,2,KM)+F2FZ(KM)*ALSGS(I,2,K))
              ALSGS1(I,1,K,3) = &
                C2CXI(I)*0.5*(F2FX(I)*ALSGS(IM,2,K)+F2FX(IM)*ALSGS(I,2,K))
            ENDIF
            IF (BC_YTOP .EQ. 0) THEN
              ALSGS1(I,N2,K,1) = 0.
              ALSGS1(I,N2,K,2) = 0.
              ALSGS1(I,N2,K,3) = 0.
            ELSE
              ALSGS1(I,N2,K,2) = &
                C2CZI(K)*0.5*(F2FZ(K)*ALSGS(I,N2M,KM)+F2FZ(KM)*ALSGS(I,N2M,K))
              ALSGS1(I,N2,K,3) = &
                C2CXI(I)*0.5*(F2FX(I)*ALSGS(IM,N2M,K)+F2FX(IM)*ALSGS(I,N2M,K))
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!!$OMP END PARALLEL DO

!!$OMP PARALLEL DO private(IM,JM)
      DO I=1,N1M
        DO J=1,N2M
          IM=IMV(I)
          JM=JMV(J)
          IF (ZPRDIC .EQ. 1) THEN
            ALSGS1(I,J,N3,1) = ALSGS1(I,J,1,1)
            ALSGS1(I,J,N3,2) = ALSGS1(I,J,1,2)
            ALSGS1(I,J,N3,3) = ALSGS1(I,J,1,3)
          ELSE
            IF (BC_ZBTM .EQ. 0) THEN
              ALSGS1(I,J,1,1) = 0.
              ALSGS1(I,J,1,2) = 0.
              ALSGS1(I,J,1,3) = 0.
            ELSE
              ALSGS1(I,J,1,2) = &
                C2CYI(J)*0.5*(F2FY(J)*ALSGS(I,JM,2)+F2FY(JM)*ALSGS(I,J,2))
              ALSGS1(I,J,1,3) = &
                C2CXI(I)*0.5*(F2FX(I)*ALSGS(IM,J,2)+F2FX(IM)*ALSGS(I,J,2))
            ENDIF
            IF (BC_ZTOP .EQ. 0) THEN
              ALSGS1(I,J,N3,1) = 0.
              ALSGS1(I,J,N3,2) = 0.
              ALSGS1(I,J,N3,3) = 0.
            ELSE
              ALSGS1(I,J,N3,2) = &
                C2CYI(K)*0.5*(F2FY(J)*ALSGS(I,JM,N3M)+F2FY(JM)*ALSGS(I,J,N3M))
              ALSGS1(I,J,N3,3) = &
                C2CXI(I)*0.5*(F2FX(I)*ALSGS(IM,J,N3M)+F2FX(IM)*ALSGS(I,J,N3M))
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!!$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE ALPINTERPOL
!=======================================================================
!=======================================================================
      SUBROUTINE RHSSGS_T
!=======================================================================
!
!     Calculate non-linear SGS HF terms
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,T,ALSGS1,RHS1
      IMPLICIT NONE
      INTEGER*8 :: I,J,K
      REAL*8    :: RHST1,RHST2,RHST3

!$OMP PARALLEL DO&
!$OMP private(RHST1,RHST2,RHST3)
      DO K=1,N3M
        DO J=1,N2M
          DO I=1,N1M
            RHST1=0.
            RHST2=0.
            RHST3=0.
            RHS1(I,J,K,4)=RHST1+RHST2+RHST3
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE RHSSGS_T
!=======================================================================
