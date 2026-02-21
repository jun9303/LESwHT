!=======================================================================
      subroutine sgsfilterinit
!=======================================================================
!
!     FILTER COEFFICIENTS FOR THE SIMPSON'S RULE
!     : FOR DYNAMIC PROCEDURES TO DETERMINE THE MODEL COEFFICIENT,
!       A TEST FILTER IS APPLIED BASED ON THE SIMPSON'S RULE.
!
!     CFX1, CFX2  : FILTER COEFFICIENTS IN X-DIRECTION
!     CFY1, CFY2  : FILTER COEFFICIENTS IN Y-DIRECTION
!     CFZ1, CFZ2  : FILTER COEFFICIENTS IN Z-DIRECTION
!
!-----------------------------------------------------------------------
        use mod_common
        use mod_flowarray, only: cfx1, cfx2, cfy1, cfy2, cfz1, cfz2
        implicit none
        integer(8) :: i, j, k
        real(8) :: h1, h2, h3, h4

        cfx1 = 0.
        cfx2 = 0.
        cfy1 = 0.
        cfy2 = 0.
        cfz1 = 0.
        cfz2 = 0.

        if (xprdic .eq. 1) then

          if (n1m .eq. 1) goto 101

!$OMP PARALLEL DO PRIVATE(H1,H2)
          do i = 1, n1m
            h1 = c2cx(i)
            h2 = c2cx(i + 1)
            cfx1(i, -1) = (2.*h1 - h2) / 6./h1
            cfx1(i, 0) = ((h1 + h2)**2) / 6./h1 / h2
            cfx1(i, 1) = (2.*h2 - h1) / 6./h2
          end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(H1,H2,H3,H4)
          do i = 2, n1m - 1
            h1 = c2cx(i - 1)
            h2 = c2cx(i)
            h3 = c2cx(i + 1)
            h4 = c2cx(i + 2)
            cfx2(i, -2) = (2.*h1 - h2) / 6./h1 / 2.
            cfx2(i, -1) = ((h1 + h2)**2) / 6./h1 / h2 / 2.
            cfx2(i, 0) = ((2.*h2 - h1) / 6./h2 + (2.*h3 - h4) / 6./h3) / 2.
            cfx2(i, 1) = ((h3 + h4)**2) / 6./h3 / h4 / 2.
            cfx2(i, 2) = (2.*h4 - h3) / 6./h4 / 2.
          end do
!$OMP END PARALLEL DO
          h1 = c2cx(n1m)
          h2 = c2cx(1)
          h3 = c2cx(2)
          h4 = c2cx(3)
          cfx2(1, -2) = (2.*h1 - h2) / 6./h1 / 2.
          cfx2(1, -1) = ((h1 + h2)**2) / 6./h1 / h2 / 2.
          cfx2(1, 0) = ((2.*h2 - h1) / 6./h2 + (2.*h3 - h4) / 6./h3) / 2.
          cfx2(1, 1) = ((h3 + h4)**2) / 6./h3 / h4 / 2.
          cfx2(1, 2) = (2.*h4 - h3) / 6./h4 / 2.
          h1 = c2cx(n1m - 1)
          h2 = c2cx(n1m)
          h3 = c2cx(1)
          h4 = c2cx(2)
          cfx2(n1m, -2) = (2.*h1 - h2) / 6./h1 / 2.
          cfx2(n1m, -1) = ((h1 + h2)**2) / 6./h1 / h2 / 2.
          cfx2(n1m, 0) = ((2.*h2 - h1) / 6./h2 + (2.*h3 - h4) / 6./h3) / 2.
          cfx2(n1m, 1) = ((h3 + h4)**2) / 6./h3 / h4 / 2.
          cfx2(n1m, 2) = (2.*h4 - h3) / 6./h4 / 2.

        else

          if (n1m .eq. 1) goto 101

!$OMP PARALLEL DO PRIVATE(H1,H2)
          do i = 1, n1m
            h1 = c2cx(i)
            h2 = c2cx(i + 1)
            cfx1(i, -1) = (2.*h1 - h2) / 6./h1
            cfx1(i, 0) = ((h1 + h2)**2) / 6./h1 / h2
            cfx1(i, 1) = (2.*h2 - h1) / 6./h2
          end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(H1,H2,H3,H4)
          do i = 3, n1m - 2
            h1 = c2cx(i - 1)
            h2 = c2cx(i)
            h3 = c2cx(i + 1)
            h4 = c2cx(i + 2)
            cfx2(i, -2) = (2.*h1 - h2) / 6./h1 / 2.
            cfx2(i, -1) = ((h1 + h2)**2) / 6./h1 / h2 / 2.
            cfx2(i, 0) = ((2.*h2 - h1) / 6./h2 + (2.*h3 - h4) / 6./h3) / 2.
            cfx2(i, 1) = ((h3 + h4)**2) / 6./h3 / h4 / 2.
            cfx2(i, 2) = (2.*h4 - h3) / 6./h4 / 2.
          end do
!$OMP END PARALLEL DO
          h1 = c2cx(2)
          h2 = c2cx(3)
          cfx2(2, -1) = (2.*h1 - h2) / 6./h1
          cfx2(2, 0) = ((h1 + h2)**2) / 6./h1 / h2
          cfx2(2, 1) = (2.*h2 - h1) / 6./h2
          h1 = c2cx(n1m - 1)
          h2 = c2cx(n1m)
          cfx2(n1m - 1, -1) = (2.*h1 - h2) / 6./h1
          cfx2(n1m - 1, 0) = ((h1 + h2)**2) / 6./h1 / h2
          cfx2(n1m - 1, 1) = (2.*h2 - h1) / 6./h2

        end if

101     continue

        if (yprdic .eq. 1) then

          if (n2m .eq. 1) goto 102

!$OMP PARALLEL DO PRIVATE(H1,H2)
          do j = 1, n2m
            h1 = c2cy(j)
            h2 = c2cy(j + 1)
            cfy1(j, -1) = (2.*h1 - h2) / 6./h1
            cfy1(j, 0) = ((h1 + h2)**2) / 6./h1 / h2
            cfy1(j, 1) = (2.*h2 - h1) / 6./h2
          end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(H1,H2,H3,H4)
          do j = 2, n2m - 1
            h1 = c2cy(j - 1)
            h2 = c2cy(j)
            h3 = c2cy(j + 1)
            h4 = c2cy(j + 2)
            cfy2(j, -2) = (2.*h1 - h2) / 6./h1 / 2.
            cfy2(j, -1) = ((h1 + h2)**2) / 6./h1 / h2 / 2.
            cfy2(j, 0) = ((2.*h2 - h1) / 6./h2 + (2.*h3 - h4) / 6./h3) / 2.
            cfy2(j, 1) = ((h3 + h4)**2) / 6./h3 / h4 / 2.
            cfy2(j, 2) = (2.*h4 - h3) / 6./h4 / 2.
          end do
!$OMP END PARALLEL DO
          h1 = c2cy(n2m)
          h2 = c2cy(1)
          h3 = c2cy(2)
          h4 = c2cy(3)
          cfy2(1, -2) = (2.*h1 - h2) / 6./h1 / 2.
          cfy2(1, -1) = ((h1 + h2)**2) / 6./h1 / h2 / 2.
          cfy2(1, 0) = ((2.*h2 - h1) / 6./h2 + (2.*h3 - h4) / 6./h3) / 2.
          cfy2(1, 1) = ((h3 + h4)**2) / 6./h3 / h4 / 2.
          cfy2(1, 2) = (2.*h4 - h3) / 6./h4 / 2.
          h1 = c2cy(n2m - 1)
          h2 = c2cy(n2m)
          h3 = c2cy(1)
          h4 = c2cy(2)
          cfy2(n2m, -2) = (2.*h1 - h2) / 6./h1 / 2.
          cfy2(n2m, -1) = ((h1 + h2)**2) / 6./h1 / h2 / 2.
          cfy2(n2m, 0) = ((2.*h2 - h1) / 6./h2 + (2.*h3 - h4) / 6./h3) / 2.
          cfy2(n2m, 1) = ((h3 + h4)**2) / 6./h3 / h4 / 2.
          cfy2(n2m, 2) = (2.*h4 - h3) / 6./h4 / 2.

        else

          if (n2m .eq. 1) goto 102
!$OMP PARALLEL DO PRIVATE(H1,H2)
          do j = 1, n2m
            h1 = c2cy(j)
            h2 = c2cy(j + 1)
            cfy1(j, -1) = (2.*h1 - h2) / 6./h1
            cfy1(j, 0) = ((h1 + h2)**2) / 6./h1 / h2
            cfy1(j, 1) = (2.*h2 - h1) / 6./h2

          end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(H1,H2,H3,H4)
          do j = 3, n2m - 2
            h1 = c2cy(j - 1)
            h2 = c2cy(j)
            h3 = c2cy(j + 1)
            h4 = c2cy(j + 2)
            cfy2(j, -2) = (2.*h1 - h2) / 6./h1 / 2.
            cfy2(j, -1) = ((h1 + h2)**2) / 6./h1 / h2 / 2.
            cfy2(j, 0) = ((2.*h2 - h1) / 6./h2 + (2.*h3 - h4) / 6./h3) / 2.
            cfy2(j, 1) = ((h3 + h4)**2) / 6./h3 / h4 / 2.
            cfy2(j, 2) = (2.*h4 - h3) / 6./h4 / 2.
          end do
!$OMP END PARALLEL DO
          h1 = c2cy(2)
          h2 = c2cy(3)
          cfy2(2, -1) = (2.*h1 - h2) / 6./h1
          cfy2(2, 0) = ((h1 + h2)**2) / 6./h1 / h2
          cfy2(2, 1) = (2.*h2 - h1) / 6./h2
          h1 = c2cy(n2m - 1)
          h2 = c2cy(n2m)
          cfy2(n2m - 1, -1) = (2.*h1 - h2) / 6./h1
          cfy2(n2m - 1, 0) = ((h1 + h2)**2) / 6./h1 / h2
          cfy2(n2m - 1, 1) = (2.*h2 - h1) / 6./h2
        end if

102     continue

        if (zprdic .eq. 1) then

          if (n3m .eq. 1) goto 103

!$OMP PARALLEL DO PRIVATE(H1,H2)
          do k = 1, n3m
            h1 = c2cz(k)
            h2 = c2cz(k + 1)
            cfz1(k, -1) = (2.*h1 - h2) / 6./h1
            cfz1(k, 0) = ((h1 + h2)**2) / 6./h1 / h2
            cfz1(k, 1) = (2.*h2 - h1) / 6./h2
          end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(H1,H2,H3,H4)
          do k = 2, n3m - 1
            h1 = c2cz(k - 1)
            h2 = c2cz(k)
            h3 = c2cz(k + 1)
            h4 = c2cz(k + 2)
            cfz2(k, -2) = (2.*h1 - h2) / 6./h1 / 2.
            cfz2(k, -1) = ((h1 + h2)**2) / 6./h1 / h2 / 2.
            cfz2(k, 0) = ((2.*h2 - h1) / 6./h2 + (2.*h3 - h4) / 6./h3) / 2.
            cfz2(k, 1) = ((h3 + h4)**2) / 6./h3 / h4 / 2.
            cfz2(k, 2) = (2.*h4 - h3) / 6./h4 / 2.
          end do
!$OMP END PARALLEL DO
          h1 = c2cz(n3m)
          h2 = c2cz(1)
          h3 = c2cz(2)
          h4 = c2cz(3)
          cfz2(1, -2) = (2.*h1 - h2) / 6./h1 / 2.
          cfz2(1, -1) = ((h1 + h2)**2) / 6./h1 / h2 / 2.
          cfz2(1, 0) = ((2.*h2 - h1) / 6./h2 + (2.*h3 - h4) / 6./h3) / 2.
          cfz2(1, 1) = ((h3 + h4)**2) / 6./h3 / h4 / 2.
          cfz2(1, 2) = (2.*h4 - h3) / 6./h4 / 2.
          h1 = c2cz(n3m - 1)
          h2 = c2cz(n3m)
          h3 = c2cz(1)
          h4 = c2cz(2)
          cfz2(n3m, -2) = (2.*h1 - h2) / 6./h1 / 2.
          cfz2(n3m, -1) = ((h1 + h2)**2) / 6./h1 / h2 / 2.
          cfz2(n3m, 0) = ((2.*h2 - h1) / 6./h2 + (2.*h3 - h4) / 6./h3) / 2.
          cfz2(n3m, 1) = ((h3 + h4)**2) / 6./h3 / h4 / 2.
          cfz2(n3m, 2) = (2.*h4 - h3) / 6./h4 / 2.

        else

          if (n3m .eq. 1) goto 103

!$OMP PARALLEL DO PRIVATE(H1,H2)
          do k = 1, n3m
            h1 = c2cz(k)
            h2 = c2cz(k + 1)
            cfz1(k, -1) = (2.*h1 - h2) / 6./h1
            cfz1(k, 0) = ((h1 + h2)**2) / 6./h1 / h2
            cfz1(k, 1) = (2.*h2 - h1) / 6./h2
          end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(H1,H2,H3,H4)
          do k = 3, n3m - 2
            h1 = c2cz(k - 1)
            h2 = c2cz(k)
            h3 = c2cz(k + 1)
            h4 = c2cz(k + 2)
            cfz2(k, -2) = (2.*h1 - h2) / 6./h1 / 2.
            cfz2(k, -1) = ((h1 + h2)**2) / 6./h1 / h2 / 2.
            cfz2(k, 0) = ((2.*h2 - h1) / 6./h2 + (2.*h3 - h4) / 6./h3) / 2.
            cfz2(k, 1) = ((h3 + h4)**2) / 6./h3 / h4 / 2.
            cfz2(k, 2) = (2.*h4 - h3) / 6./h4 / 2.
          end do
!$OMP END PARALLEL DO
          h1 = c2cz(2)
          h2 = c2cz(3)
          cfz2(2, -1) = (2.*h1 - h2) / 6./h1
          cfz2(2, 0) = ((h1 + h2)**2) / 6./h1 / h2
          cfz2(2, 1) = (2.*h2 - h1) / 6./h2
          h1 = c2cz(n3m - 1)
          h2 = c2cz(n3m)
          cfz2(n3m - 1, -1) = (2.*h1 - h2) / 6./h1
          cfz2(n3m - 1, 0) = ((h1 + h2)**2) / 6./h1 / h2
          cfz2(n3m - 1, 1) = (2.*h2 - h1) / 6./h2

        end if

103     continue

        return
      end subroutine sgsfilterinit
!=======================================================================
!=======================================================================
      subroutine sgscalc
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
!$      USE OMP_LIB
        use mod_common
        use mod_flowarray, only: nusgs
        implicit none
        if (msub .eq. 1) then
          if (iread .eq. 0 .and. ntime .eq. 1) then
            nusgs = 0d0
          else
            if (insmdl .eq. 0) then
              write (*, *) 'SMAGORINSKY MODEL IS NOT IMPLEMENTED YET'
              write (*, *) 'SGS.F90 LINE 317'
              stop
            else if (insmdl .eq. 1) then
              if (idvmon .eq. 1) then
                call sgs_dvmg         ! DYNAMIC VREMAN MODEL W/ GLOBAL COEF
              else
                call sgs_cvm          ! CONSTANT VREMAN MODEL
              end if
            end if
          end if

          call nutzero               ! SET NUSGS TO ZERO IN THE SOLID BODY
          call nutinterpol           !              INTERPOLATION OF NUSGS
        end if

        return
      end subroutine sgscalc
!=======================================================================
!=======================================================================
      subroutine sgs_cvm
!=======================================================================
!
!     NUSGS=CSGSTS*SQRT(BBB/AAA)
!     CSGSTS : VREMAN CONSTANT (= CSGSTS FROM MOD_COMMON)
!     BBB    : B_{11}*B_{22}-B_{12}*B_{12}+B_{11}*B_{33}
!             -B_{13}*B_{13}+B_{22}*B_{33}-B_{23}*B_{23}
!     AAA    : A_{IJ}*A_{IJ}
!     B_{IJ} : SIGMA(DEL_{M}^2*A_{MI}*A_{MJ}) WHERE M=1, 2, 3
!     A_{IJ} : DUJ/DXI
!
!-----------------------------------------------------------------------
        use mod_common
        use mod_flowarray, only: u, v, w, nusgs
        implicit none
        integer(8) :: i, j, k, kplus, kminus, jplus, jminus, iplus, iminus
        real(8) :: vg11, vg12, vg13, vg21, vg22, vg23, vg31, vg32, vg33
        real(8) :: up, um, vp, vm, wp, wm
        real(8) :: a(9), b(6)
        real(8) :: aaa, bbb

!$OMP PARALLEL DO PRIVATE(UP,UM,VP,VM,WP,WM)&
!$OMP PRIVATE(A,B,AAA,BBB)
        do k = 1, n3m
          do j = 1, n2m
            do i = 1, n1m

              call velgrad(i, j, k, a, b)         ! B IS USED FOR DUMMY VARIABLE

              b(1) = f2fx(i)**2.*a(1) * a(1) &
                     + f2fy(j)**2.*a(4) * a(4) &
                     + f2fz(k)**2.*a(7) * a(7)
              b(2) = f2fx(i)**2.*a(1) * a(2) &
                     + f2fy(j)**2.*a(4) * a(5) &
                     + f2fz(k)**2.*a(7) * a(8)
              b(3) = f2fx(i)**2.*a(1) * a(3) &
                     + f2fy(j)**2.*a(4) * a(6) &
                     + f2fz(k)**2.*a(7) * a(9)
              b(4) = f2fx(i)**2.*a(2) * a(2) &
                     + f2fy(j)**2.*a(5) * a(5) &
                     + f2fz(k)**2.*a(8) * a(8)
              b(5) = f2fx(i)**2.*a(2) * a(3) &
                     + f2fy(j)**2.*a(5) * a(6) &
                     + f2fz(k)**2.*a(8) * a(9)
              b(6) = f2fx(i)**2.*a(3) * a(3) &
                     + f2fy(j)**2.*a(6) * a(6) &
                     + f2fz(k)**2.*a(9) * a(9)

              aaa = a(1)**2.+a(2)**2.+a(3)**2. &
                    +a(4)**2.+a(5)**2.+a(6)**2. &
                    +a(7)**2.+a(8)**2.+a(9)**2.
              bbb = b(1) * b(4) + b(4) * b(6) + b(6) * b(1) &
                    - b(2)**2.-b(3)**2.-b(5)**2.

              bbb = dmax1(bbb, 1.0e-16)

              if (aaa .eq. 0.) then
                nusgs(i, j, k) = 0d0
              else
                nusgs(i, j, k) = csgsts * sqrt(bbb / aaa)
              end if
            end do
          end do
        end do
!$OMP END PARALLEL DO

        return
      end subroutine sgs_cvm
!=======================================================================
!=======================================================================
      subroutine sgs_dvmg
!=======================================================================
!
!     NUSGS=CSGSTS*SQRT(BBB/AAA)
!     CSGSTS   : VREMAN CONSTANT, DETERMINED FROM LSM OF GERMANO IDENTITY
!     BBB    : B_{11}*B_{22}-B_{12}*B_{12}+B_{11}*B_{33}
!             -B_{13}*B_{13}+B_{22}*B_{33}-B_{23}*B_{23}
!     AAA    : A_{IJ}*A_{IJ}
!     B_{IJ} : SIGMA(DEL_{M}^2*A_{MI}*A_{MJ}) WHERE M=1, 2, 3
!     A_{IJ} : DUJ/DXI
!
!-----------------------------------------------------------------------
        use mod_common
        use mod_flowarray, only: u, v, w, nusgs, nwall_dvm, aalp, llij, mmij, uui
        implicit none
        integer(8) :: i, j, k
        real(8) :: ljmj_v, mjmj_v
        real(8) :: sdxf2, sdyf2, sdzf2, aaa, bbb, ami, volume
        real(8) :: b(6), str(6), alp(9)

        ljmj_v = 0.
        mjmj_v = 0.
        csgsts = 0.

        allocate (aalp(n1m, n2m, n3m, 9))
        allocate (llij(n1m, n2m, n3m, 6))
        allocate (mmij(n1m, n2m, n3m, 6))
        allocate (uui(n1m, n2m, n3m, 3))

!$OMP PARALLEL DO&
!$OMP PRIVATE(B,STR,ALP)&
!$OMP PRIVATE(SDXF2,SDYF2,SDZF2,AAA,BBB)
        do k = 1, n3m
          do j = 1, n2m
            do i = 1, n1m
              call velgrad(i, j, k, alp, str)
              sdxf2 = f2fx(i)**2.
              sdyf2 = f2fy(j)**2.
              sdzf2 = f2fz(k)**2.
              b(1) = sdxf2 * alp(1) * alp(1) &
                     + sdyf2 * alp(4) * alp(4) &
                     + sdzf2 * alp(7) * alp(7)
              b(2) = sdxf2 * alp(1) * alp(2) &
                     + sdyf2 * alp(4) * alp(5) &
                     + sdzf2 * alp(7) * alp(8)
              b(3) = sdxf2 * alp(1) * alp(3) &
                     + sdyf2 * alp(4) * alp(6) &
                     + sdzf2 * alp(7) * alp(9)
              b(4) = sdxf2 * alp(2) * alp(2) &
                     + sdyf2 * alp(5) * alp(5) &
                     + sdzf2 * alp(8) * alp(8)
              b(5) = sdxf2 * alp(2) * alp(3) &
                     + sdyf2 * alp(5) * alp(6) &
                     + sdzf2 * alp(8) * alp(9)
              b(6) = sdxf2 * alp(3) * alp(3) &
                     + sdyf2 * alp(6) * alp(6) &
                     + sdzf2 * alp(9) * alp(9)
              aaa = alp(1)**2.+alp(2)**2.+alp(3)**2. &
                    +alp(4)**2.+alp(5)**2.+alp(6)**2. &
                    +alp(7)**2.+alp(8)**2.+alp(9)**2.
              bbb = b(1) * b(4) + b(4) * b(6) + b(6) * b(1) &
                    - b(2)**2.-b(3)**2.-b(5)**2.
              bbb = dmax1(bbb, 1.0e-16)
              if (aaa .eq. 0.) then
                nusgs(i, j, k) = 0d0
              else
                nusgs(i, j, k) = sqrt(bbb / aaa)
              end if
              mmij(i, j, k, 1) = nusgs(i, j, k) * str(1)
              mmij(i, j, k, 2) = nusgs(i, j, k) * str(2)
              mmij(i, j, k, 3) = nusgs(i, j, k) * str(3)
              mmij(i, j, k, 4) = nusgs(i, j, k) * str(4)
              mmij(i, j, k, 5) = nusgs(i, j, k) * str(5)
              mmij(i, j, k, 6) = nusgs(i, j, k) * str(6)
              aalp(i, j, k, 1) = alp(1)
              aalp(i, j, k, 2) = alp(2)
              aalp(i, j, k, 3) = alp(3)
              aalp(i, j, k, 4) = alp(4)
              aalp(i, j, k, 5) = alp(5)
              aalp(i, j, k, 6) = alp(6)
              aalp(i, j, k, 7) = alp(7)
              aalp(i, j, k, 8) = alp(8)
              aalp(i, j, k, 9) = alp(9)
            end do
          end do
        end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
        do i = 1, 6
          call test_filter(mmij(:, :, :, i), filter)
        end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
        do j = 1, 9
          call test_filter(aalp(:, :, :, j), filter)
        end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO&
!$OMP PRIVATE(B,STR,ALP)&
!$OMP PRIVATE(SDXF2,SDYF2,SDZF2,AAA,BBB,AMI)
        do k = 1, n3m
          do j = 1, n2m
            do i = 1, n1m
              call velgrad(i, j, k, alp, str)
              sdxf2 = (0.5 * f2fx(i - 1) * (1.-fixiu(i)) + f2fx(i) &
                       + 0.5 * f2fx(i + 1) * (1.-fixil(i)))**2.
              sdyf2 = (0.5 * f2fy(j - 1) * (1.-fixju(j)) + f2fy(j) &
                       + 0.5 * f2fy(j + 1) * (1.-fixjl(j)))**2.
              sdzf2 = (0.5 * f2fz(k - 1) * (1.-fixku(k)) + f2fz(k) &
                       + 0.5 * f2fz(k + 1) * (1.-fixkl(k)))**2.
              b(1) = sdxf2 * alp(1) * alp(1) &
                     + sdyf2 * alp(4) * alp(4) &
                     + sdzf2 * alp(7) * alp(7)
              b(2) = sdxf2 * alp(1) * alp(2) &
                     + sdyf2 * alp(4) * alp(5) &
                     + sdzf2 * alp(7) * alp(8)
              b(3) = sdxf2 * alp(1) * alp(3) &
                     + sdyf2 * alp(4) * alp(6) &
                     + sdzf2 * alp(7) * alp(9)
              b(4) = sdxf2 * alp(2) * alp(2) &
                     + sdyf2 * alp(5) * alp(5) &
                     + sdzf2 * alp(8) * alp(8)
              b(5) = sdxf2 * alp(2) * alp(3) &
                     + sdyf2 * alp(5) * alp(6) &
                     + sdzf2 * alp(8) * alp(9)
              b(6) = sdxf2 * alp(3) * alp(3) &
                     + sdyf2 * alp(6) * alp(6) &
                     + sdzf2 * alp(9) * alp(9)
              aaa = alp(1)**2.+alp(2)**2.+alp(3)**2. &
                    +alp(4)**2.+alp(5)**2.+alp(6)**2. &
                    +alp(7)**2.+alp(8)**2.+alp(9)**2.
              bbb = b(1) * b(4) + b(4) * b(6) + b(6) * b(1) &
                    - b(2)**2.-b(3)**2.-b(5)**2.
              bbb = dmax1(bbb, 1.0e-16)
              if (aaa .eq. 0.) then
                ami = 0d0
              else
                ami = sqrt(bbb / aaa)
              end if
              mmij(i, j, k, 1) = ami * str(1) - mmij(i, j, k, 1)
              mmij(i, j, k, 2) = ami * str(2) - mmij(i, j, k, 2)
              mmij(i, j, k, 3) = ami * str(3) - mmij(i, j, k, 3)
              mmij(i, j, k, 4) = ami * str(4) - mmij(i, j, k, 4)
              mmij(i, j, k, 5) = ami * str(5) - mmij(i, j, k, 5)
              mmij(i, j, k, 6) = ami * str(6) - mmij(i, j, k, 6)
              uui(i, j, k, 1) = 0.5 * (u(i, j, k) + u(i + 1, j, k))
              uui(i, j, k, 2) = 0.5 * (v(i, j, k) + v(i, j + 1, k))
              uui(i, j, k, 3) = 0.5 * (w(i, j, k) + w(i, j, k + 1))
              llij(i, j, k, 1) = uui(i, j, k, 1) * uui(i, j, k, 1)
              llij(i, j, k, 2) = uui(i, j, k, 1) * uui(i, j, k, 2)
              llij(i, j, k, 3) = uui(i, j, k, 1) * uui(i, j, k, 3)
              llij(i, j, k, 4) = uui(i, j, k, 2) * uui(i, j, k, 2)
              llij(i, j, k, 5) = uui(i, j, k, 2) * uui(i, j, k, 3)
              llij(i, j, k, 6) = uui(i, j, k, 3) * uui(i, j, k, 3)
            end do
          end do
        end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
        do i = 1, 3
          call test_filter(uui(:, :, :, i), filter)
        end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
        do j = 1, 6
          call test_filter(llij(:, :, :, j), filter)
        end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO&
!$OMP REDUCTION(+:LJMJ_V,MJMJ_V)
        do k = 1, n3m
          do j = 1, n2m
            do i = 1, n1m
              volume = f2fx(i) * f2fy(j) * f2fz(k)
              llij(i, j, k, 1) = llij(i, j, k, 1) - uui(i, j, k, 1) * uui(i, j, k, 1)
              llij(i, j, k, 2) = llij(i, j, k, 2) - uui(i, j, k, 1) * uui(i, j, k, 2)
              llij(i, j, k, 3) = llij(i, j, k, 3) - uui(i, j, k, 1) * uui(i, j, k, 3)
              llij(i, j, k, 4) = llij(i, j, k, 4) - uui(i, j, k, 2) * uui(i, j, k, 2)
              llij(i, j, k, 5) = llij(i, j, k, 5) - uui(i, j, k, 2) * uui(i, j, k, 3)
              llij(i, j, k, 6) = llij(i, j, k, 6) - uui(i, j, k, 3) * uui(i, j, k, 3)
              ljmj_v = ljmj_v &
                       + (2.*llij(i, j, k, 2) * mmij(i, j, k, 2) + llij(i, j, k, 1) * mmij(i, j, k, 1) &
                          + 2.*llij(i, j, k, 3) * mmij(i, j, k, 3) + llij(i, j, k, 4) * mmij(i, j, k, 4) &
                          + 2.*llij(i, j, k, 5) * mmij(i, j, k, 5) + llij(i, j, k, 6) * mmij(i, j, k, 6)) &
                       * float(nwall_dvm(i, j, k)) * volume
              mjmj_v = mjmj_v &
                       + (2.*mmij(i, j, k, 2)**2.+mmij(i, j, k, 1)**2. &
                          +2.*mmij(i, j, k, 3)**2.+mmij(i, j, k, 4)**2. &
                          +2.*mmij(i, j, k, 5)**2.+mmij(i, j, k, 6)**2.) &
                       * float(nwall_dvm(i, j, k)) * volume
            end do
          end do
        end do
!$OMP END PARALLEL DO

        deallocate (aalp)
        deallocate (llij)
        deallocate (mmij)
        deallocate (uui)

        csgsts = dmax1(0., -0.5 * ljmj_v / mjmj_v)

!$OMP PARALLEL DO
        do k = 1, n3m
          do j = 1, n2m
            do i = 1, n1m
              nusgs(i, j, k) = csgsts * nusgs(i, j, k)
            end do
          end do
        end do
!$OMP END PARALLEL DO

        return
      end subroutine sgs_dvmg
!=======================================================================
!=======================================================================
      subroutine test_filter(a0, dir)
!=======================================================================
        use mod_common
        use mod_flowarray, only: cfx1, cfy1, cfz1
        implicit none
        integer(8) :: dir
        real(8) :: a0(n1m, n2m, n3m)

        integer(8) :: i, j, k
        integer(8) :: xfilter, yfilter, zfilter
        real(8) :: tmp(n1m, n2m, n3m)

        if (dir .eq. 1) then
          xfilter = 1
          yfilter = 0
          zfilter = 0
        else if (dir .eq. 2) then
          xfilter = 0
          yfilter = 0
          zfilter = 1
        else if (dir .eq. 3) then
          xfilter = 1
          yfilter = 0
          zfilter = 1
        else
          write (*, *) 'INVALID FILTER INPUT.'
          stop
        end if

        tmp = a0

        if (xfilter .eq. 1) then
          if (xprdic .ne. 1) then
            do k = 1, n3m
              do j = 1, n2m
                do i = 2, n1m - 1
                  tmp(i, j, k) = cfx1(i, -1) * tmp(i - 1, j, k) &
                                 + cfx1(i, 0) * tmp(i, j, k) &
                                 + cfx1(i, 1) * tmp(i + 1, j, k)
                end do
              end do
            end do
          else
            do j = 1, n2m
              do k = 1, n3m
                do i = 2, n1m - 1
                  tmp(i, j, k) = cfx1(i, -1) * tmp(i - 1, j, k) &
                                 + cfx1(i, 0) * tmp(i, j, k) &
                                 + cfx1(i, 1) * tmp(i + 1, j, k)
                end do
                tmp(1, j, k) = cfx1(i, -1) * tmp(n1m, j, k) &
                               + cfx1(i, 0) * tmp(1, j, k) &
                               + cfx1(i, 1) * tmp(2, j, k)
                tmp(n1m, j, k) = cfx1(i, -1) * tmp(n1m - 1, j, k) &
                                 + cfx1(i, 0) * tmp(n1m, j, k) &
                                 + cfx1(i, 1) * tmp(1, j, k)
              end do
            end do
          end if
        end if

        if (yfilter .eq. 1) then
          if (yprdic .ne. 1) then
            do k = 1, n3m
              do j = 2, n2m - 1
                do i = 1, n1m
                  tmp(i, j, k) = cfy1(j, -1) * tmp(i, j - 1, k) &
                                 + cfy1(j, 0) * tmp(i, j, k) &
                                 + cfy1(j, 1) * tmp(i, j + 1, k)
                end do
              end do
            end do
          else
            do k = 1, n3m
              do i = 1, n1m
                do j = 2, n2m - 1
                  tmp(i, j, k) = cfy1(j, -1) * tmp(i, j - 1, k) &
                                 + cfy1(j, 0) * tmp(i, j, k) &
                                 + cfy1(j, 1) * tmp(i, j + 1, k)
                end do
                tmp(i, 1, k) = cfy1(j, -1) * tmp(i, n2m, k) &
                               + cfy1(j, 0) * tmp(i, 1, k) &
                               + cfy1(j, 1) * tmp(i, 2, k)
                tmp(i, n2m, k) = cfy1(j, -1) * tmp(i, n2m - 1, k) &
                                 + cfy1(j, 0) * tmp(i, n2m, k) &
                                 + cfy1(j, 1) * tmp(i, 1, k)
              end do
            end do
          end if
        end if

        if (zfilter .eq. 1) then
          if (zprdic .ne. 1) then
            do k = 2, n3m - 1
              do j = 1, n2m
                do i = 1, n1m
                  tmp(i, j, k) = cfz1(k, -1) * tmp(i, j, k - 1) &
                                 + cfz1(k, 0) * tmp(i, j, k) &
                                 + cfz1(k, 1) * tmp(i, j, k + 1)
                end do
              end do
            end do
          else
            do i = 1, n1m
              do j = 1, n2m
                do k = 2, n3m - 1
                  tmp(i, j, k) = cfz1(k, -1) * tmp(i, j, k - 1) &
                                 + cfz1(k, 0) * tmp(i, j, k) &
                                 + cfz1(k, 1) * tmp(i, j, k + 1)
                end do
                tmp(i, j, 1) = cfz1(k, -1) * tmp(i, j, n3m) &
                               + cfz1(k, 0) * tmp(i, j, 1) &
                               + cfz1(k, 1) * tmp(i, j, 2)
                tmp(i, j, n3m) = cfz1(k, -1) * tmp(i, j, n3m - 1) &
                                 + cfz1(k, 0) * tmp(i, j, n3m) &
                                 + cfz1(k, 1) * tmp(i, j, 1)
              end do
            end do
          end if
        end if

        a0 = tmp

        return
      end subroutine test_filter
!=======================================================================
!=======================================================================
      subroutine velgrad(i, j, k, a, sr)
!=======================================================================
!$      USE OMP_LIB
        use mod_common
        use mod_flowarray, only: u, v, w, nusgs
        implicit none
        integer(8) :: i, j, k
        real(8) :: a(9), sr(6)
        integer(8) :: kplus, kminus, jplus, jminus, iplus, iminus
        real(8) :: vg11, vg12, vg13, vg21, vg22, vg23, vg31, vg32, vg33
        real(8) :: up, um, vp, vm, wp, wm

        kplus = kpv(k)
        kminus = kmv(k)
        jplus = jpv(j)
        jminus = jmv(j)
        iplus = ipv(i)
        iminus = imv(i)

        vg11 = f2fxi(i) * (u(iplus, j, k) - u(i, j, k))

        up = c2cyi(jplus) * 0.25 &
             * (f2fy(jplus) * (u(i, j, k) + u(iplus, j, k)) &
                + f2fy(j) * (u(i, jplus, k) + u(iplus, jplus, k))) &
             * (1.-fixju(j)) + 0.5 * (u(i, n2, k) + u(iplus, n2, k)) * fixju(j)
        um = c2cyi(j) * 0.25 &
             * (f2fy(j) * (u(i, jminus, k) + u(iplus, jminus, k)) &
                + f2fy(jminus) * (u(i, j, k) + u(iplus, j, k))) &
             * (1.-fixjl(j)) + 0.5 * (u(i, 0, k) + u(iplus, 0, k)) * fixjl(j)
        vg12 = f2fyi(j) * (up - um)

        up = c2czi(kplus) * 0.25 &
             * (f2fz(kplus) * (u(i, j, k) + u(iplus, j, k)) &
                + f2fz(k) * (u(i, j, kplus) + u(iplus, j, kplus))) &
             * (1.-fixku(k)) + 0.5 * (u(i, j, n3) + u(iplus, j, n3)) * fixku(k)
        um = c2czi(k) * 0.25 &
             * (f2fz(k) * (u(i, j, kminus) + u(iplus, j, kminus)) &
                + f2fz(kminus) * (u(i, j, k) + u(iplus, j, k))) &
             * (1.-fixkl(k)) + 0.5 * (u(i, j, 0) + u(iplus, j, 0)) * fixkl(k)
        vg13 = f2fzi(k) * (up - um)

        vp = c2cxi(iplus) * 0.25 &
             * (f2fx(iplus) * (v(i, j, k) + v(i, jplus, k)) &
                + f2fx(i) * (v(iplus, j, k) + v(iplus, jplus, k))) &
             * (1.-fixiu(i)) + 0.5 * (v(n1, j, k) + v(n1, jplus, k)) * fixiu(i)
        vm = c2cxi(i) * 0.25 &
             * (f2fx(i) * (v(iminus, j, k) + v(iminus, jplus, k)) &
                + f2fx(iminus) * (v(i, j, k) + v(i, jplus, k))) &
             * (1.-fixil(i)) + 0.5 * (v(0, j, k) + v(0, jplus, k)) * fixil(i)
        vg21 = f2fxi(i) * (vp - vm)

        vg22 = f2fyi(j) * (v(i, jplus, k) - v(i, j, k))

        vp = c2czi(kplus) * 0.25 &
             * (f2fz(kplus) * (v(i, j, k) + v(i, jplus, k)) &
                + f2fz(k) * (v(i, j, kplus) + v(i, jplus, kplus))) &
             * (1.-fixku(k)) + 0.5 * (v(i, j, n3) + v(i, jplus, n3)) * fixku(k)
        vm = c2czi(k) * 0.25 &
             * (f2fz(k) * (v(i, j, kminus) + v(i, jplus, kminus)) &
                + f2fz(kminus) * (v(i, j, k) + v(i, jplus, k))) &
             * (1.-fixkl(k)) + 0.5 * (v(i, j, 0) + v(i, jplus, 0)) * fixkl(k)
        vg23 = f2fzi(k) * (vp - vm)

        wp = c2cxi(iplus) * 0.25 &
             * (f2fx(iplus) * (w(i, j, k) + w(i, j, kplus)) &
                + f2fx(i) * (w(iplus, j, k) + w(iplus, j, kplus))) &
             * (1.-fixiu(i)) + 0.5 * (w(n1, j, k) + w(n1, j, kplus)) * fixiu(i)
        wm = c2cxi(i) * 0.25 &
             * (f2fx(i) * (w(iminus, j, k) + w(iminus, j, kplus)) &
                + f2fx(iminus) * (w(i, j, k) + w(i, j, kplus))) &
             * (1.-fixil(i)) + 0.5 * (w(0, j, k) + w(0, j, kplus)) * fixil(i)
        vg31 = f2fxi(i) * (wp - wm)

        wp = c2cyi(jplus) * 0.25 &
             * (f2fy(jplus) * (w(i, j, k) + w(i, j, kplus)) &
                + f2fy(j) * (w(i, jplus, k) + w(i, jplus, kplus))) &
             * (1.-fixju(j)) + 0.5 * (w(i, n2, k) + w(i, n2, kplus)) * fixju(j)
        wm = c2cyi(j) * 0.25 &
             * (f2fy(j) * (w(i, jminus, k) + w(i, jminus, kplus)) &
                + f2fy(jminus) * (w(i, j, k) + w(i, j, kplus))) &
             * (1.-fixjl(j)) + 0.5 * (w(i, 0, k) + w(i, 0, kplus)) * fixjl(j)
        vg32 = f2fyi(j) * (wp - wm)

        vg33 = f2fzi(k) * (w(i, j, kplus) - w(i, j, k))

        a(1) = vg11
        a(2) = vg21
        a(3) = vg31
        a(4) = vg12
        a(5) = vg22
        a(6) = vg32
        a(7) = vg13
        a(8) = vg23
        a(9) = vg33

        sr(1) = vg11
        sr(2) = 0.5 * (vg12 + vg21)
        sr(3) = 0.5 * (vg13 + vg31)
        sr(4) = vg22
        sr(5) = 0.5 * (vg23 + vg32)
        sr(6) = vg33

        return
      end subroutine velgrad
!=======================================================================
!=======================================================================
      subroutine nutzero
!=======================================================================
        use mod_common
        use mod_flowarray, only: nusgs, inz, jnz, knz
        implicit none
        integer(8) :: i, j, k
        integer(8) :: imax, jmax, kmax, l
        real(8) :: nusgsavg, nusgsmax, area
        real(8) :: funcbody

        imax = 0
        jmax = 0
        kmax = 0

!$OMP PARALLEL DO
        do l = 1, nzero
          nusgs(inz(l), jnz(l), knz(l)) = 0.
        end do
!$OMP END PARALLEL DO

        nusgsavg = 0.
        area = 0.
        nusgsmax = 0.

!$OMP PARALLEL DO&
!$OMP REDUCTION(+:NUSGSAVG, AREA)&
!$OMP REDUCTION(MAX:NUSGSMAX)
        do k = 1, n3m
          do j = 1, n2m
            do i = 1, n1m
              if (funcbody(xmp(i), ymp(j), zmp(k), time) .ge. 1.e-10) then
                nusgsavg = nusgsavg + nusgs(i, j, k) * f2fx(i) * f2fy(j) * f2fz(k)
                area = area + f2fx(i) * f2fy(j) * f2fz(k)
                if (nusgs(i, j, k) .gt. nusgsmax) then
                  nusgsmax = nusgs(i, j, k)
                  imax = i
                  jmax = j
                  kmax = k
                end if
              else
                nusgs(i, j, k) = 0.
              end if
            end do
          end do
        end do
!$OMP END PARALLEL DO

        nusgsavg = nusgsavg / area

        nutavg = nusgsavg * re
        nutmax = nusgsmax * re

        write (*, 99) csgsts
        write (*, 100) nutavg
        write (*, 110) nutmax, xmp(imax), ymp(jmax), zmp(kmax), imax, jmax, kmax
99      format('CSGSTS =  ', es10.3)
100     format('NUTAVG =  ', es10.3)
110     format('NUTMAX =  ', es10.3, ' @ ', 3f10.4, ' , ', 3i5)

        return
      end subroutine nutzero
!=======================================================================
!=======================================================================
      subroutine nutinterpol
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
        use mod_common
        use mod_flowarray, only: nusgs, nusgs1, inz, jnz, knz
        implicit none
        integer(8) :: i, j, k, im, jm, km

        nusgs1 = 0.

        if (xprdic .eq. 1) then
!$OMP PARALLEL DO
          do k = 0, n3
            do j = 0, n2
              nusgs(0, j, k) = nusgs(n1m, j, k)
              nusgs(n1, j, k) = nusgs(1, j, k)
            end do
          end do
!$OMP END PARALLEL DO
        end if

        if (yprdic .eq. 1) then
!$OMP PARALLEL DO
          do k = 0, n3
            do i = 0, n1
              nusgs(i, 0, k) = nusgs(i, n2m, k)
              nusgs(i, n2, k) = nusgs(i, 1, k)
            end do
          end do
!$OMP END PARALLEL DO
        end if

        if (zprdic .eq. 1) then
!$OMP PARALLEL DO
          do j = 0, n2
            do i = 0, n1
              nusgs(i, j, 0) = nusgs(i, j, n3m)
              nusgs(i, j, n3) = nusgs(i, j, 1)
            end do
          end do
!$OMP END PARALLEL DO
        end if

!$OMP PARALLEL DO
        do k = k_bgpz, n3m
          do j = j_bgpy, n2m
            do i = 1, n1m
              nusgs1(i, j, k, 1) = c2cyi(j) * c2czi(k) * 0.25 &
                                   * (f2fy(j - 1) * (f2fz(k - 1) * nusgs(i, j, k) + f2fz(k) * nusgs(i, j, k - 1)) &
                                      + f2fy(j) * (f2fz(k - 1) * nusgs(i, j - 1, k) + f2fz(k) * nusgs(i, j - 1, k - 1)))
            end do
          end do
        end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
        do k = k_bgpz, n3m
          do j = 1, n2m
            do i = i_bgpx, n1m
              nusgs1(i, j, k, 2) = c2czi(k) * c2cxi(i) * 0.25 &
                                   * (f2fz(k - 1) * (f2fx(i - 1) * nusgs(i, j, k) + f2fx(i) * nusgs(i - 1, j, k)) &
                                      + f2fz(k) * (f2fx(i - 1) * nusgs(i, j, k - 1) + f2fx(i) * nusgs(i - 1, j, k - 1)))
            end do
          end do
        end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
        do k = 1, n3m
          do j = j_bgpy, n2m
            do i = i_bgpx, n1m
              nusgs1(i, j, k, 3) = c2cxi(i) * c2cyi(j) * 0.25 &
                                   * (f2fx(i - 1) * (f2fy(j - 1) * nusgs(i, j, k) + f2fy(j) * nusgs(i, j - 1, k)) &
                                      + f2fx(i) * (f2fy(j - 1) * nusgs(i - 1, j, k) + f2fy(j) * nusgs(i - 1, j - 1, k)))
            end do
          end do
        end do
!$OMP END PARALLEL DO

!!$OMP PARALLEL DO PRIVATE(JM,KM)
        do j = j_bgpy, n2m
          do k = k_bgpz, n3m
            jm = jmv(j)
            km = kmv(k)
            if (xprdic .eq. 1) then
              nusgs1(n1, j, k, 1) = nusgs1(1, j, k, 1)
              nusgs1(n1, j, k, 2) = nusgs1(1, j, k, 2)
              nusgs1(n1, j, k, 3) = nusgs1(1, j, k, 3)
            else
              if (bc_xbtm .eq. 0) then
                nusgs1(1, j, k, 1) = 0.
                nusgs1(1, j, k, 2) = 0.
                nusgs1(1, j, k, 3) = 0.
              else
                nusgs1(1, j, k, 2) = &
                  c2czi(k) * 0.5 * (f2fz(k) * nusgs(2, j, km) + f2fz(km) * nusgs(2, j, k))
                nusgs1(1, j, k, 3) = &
                  c2cyi(j) * 0.5 * (f2fy(j) * nusgs(2, jm, k) + f2fy(jm) * nusgs(2, j, k))
              end if
              if (bc_xtop .eq. 0) then
                nusgs1(n1, j, k, 1) = 0.
                nusgs1(n1, j, k, 2) = 0.
                nusgs1(n1, j, k, 3) = 0.
              else
                nusgs1(n1, j, k, 2) = &
                  c2czi(k) * 0.5 * (f2fz(k) * nusgs(n1m, j, km) + f2fz(km) * nusgs(n1m, j, k))
                nusgs1(n1, j, k, 3) = &
                  c2cyi(j) * 0.5 * (f2fy(j) * nusgs(n1m, jm, k) + f2fy(jm) * nusgs(n1m, j, k))
              end if
            end if
          end do
        end do
!!$OMP END PARALLEL DO

!!$OMP PARALLEL DO PRIVATE(IM,KM)
        do k = k_bgpz, n3m
          do i = i_bgpx, n1m
            im = imv(i)
            km = kmv(k)
            if (yprdic .eq. 1) then
              nusgs1(i, n2, k, 1) = nusgs1(i, 1, k, 1)
              nusgs1(i, n2, k, 2) = nusgs1(i, 1, k, 2)
              nusgs1(i, n2, k, 3) = nusgs1(i, 1, k, 3)
            else
              if (bc_ybtm .eq. 0) then
                nusgs1(i, 1, k, 1) = 0.
                nusgs1(i, 1, k, 2) = 0.
                nusgs1(i, 1, k, 3) = 0.
              else
                nusgs1(i, 1, k, 2) = &
                  c2czi(k) * 0.5 * (f2fz(k) * nusgs(i, 2, km) + f2fz(km) * nusgs(i, 2, k))
                nusgs1(i, 1, k, 3) = &
                  c2cxi(i) * 0.5 * (f2fx(i) * nusgs(im, 2, k) + f2fx(im) * nusgs(i, 2, k))
              end if
              if (bc_ytop .eq. 0) then
                nusgs1(i, n2, k, 1) = 0.
                nusgs1(i, n2, k, 2) = 0.
                nusgs1(i, n2, k, 3) = 0.
              else
                nusgs1(i, n2, k, 2) = &
                  c2czi(k) * 0.5 * (f2fz(k) * nusgs(i, n2m, km) + f2fz(km) * nusgs(i, n2m, k))
                nusgs1(i, n2, k, 3) = &
                  c2cxi(i) * 0.5 * (f2fx(i) * nusgs(im, n2m, k) + f2fx(im) * nusgs(i, n2m, k))
              end if
            end if
          end do
        end do
!!$OMP END PARALLEL DO

!!$OMP PARALLEL DO PRIVATE(IM,JM)
        do i = i_bgpx, n1m
          do j = j_bgpy, n2m
            im = imv(i)
            jm = jmv(j)
            if (zprdic .eq. 1) then
              nusgs1(i, j, n3, 1) = nusgs1(i, j, 1, 1)
              nusgs1(i, j, n3, 2) = nusgs1(i, j, 1, 2)
              nusgs1(i, j, n3, 3) = nusgs1(i, j, 1, 3)
            else
              if (bc_zbtm .eq. 0) then
                nusgs1(i, j, 1, 1) = 0.
                nusgs1(i, j, 1, 2) = 0.
                nusgs1(i, j, 1, 3) = 0.
              else
                nusgs1(i, j, 1, 2) = &
                  c2cyi(j) * 0.5 * (f2fy(j) * nusgs(i, jm, 2) + f2fy(jm) * nusgs(i, j, 2))
                nusgs1(i, j, 1, 3) = &
                  c2cxi(i) * 0.5 * (f2fx(i) * nusgs(im, j, 2) + f2fx(im) * nusgs(i, j, 2))
              end if
              if (bc_ztop .eq. 0) then
                nusgs1(i, j, n3, 1) = 0.
                nusgs1(i, j, n3, 2) = 0.
                nusgs1(i, j, n3, 3) = 0.
              else
                nusgs1(i, j, n3, 2) = &
                  c2cyi(k) * 0.5 * (f2fy(j) * nusgs(i, jm, n3m) + f2fy(jm) * nusgs(i, j, n3m))
                nusgs1(i, j, n3, 3) = &
                  c2cxi(i) * 0.5 * (f2fx(i) * nusgs(im, j, n3m) + f2fx(im) * nusgs(i, j, n3m))
              end if
            end if
          end do
        end do
!!$OMP END PARALLEL DO

        return
      end subroutine nutinterpol
!=======================================================================
!=======================================================================
      subroutine rhssgs
!=======================================================================
!
!     D/DXJ(NU_T*DUJ/DXI)*(1-KRONECKERDELTA(IJ)),
!
!-----------------------------------------------------------------------
        use mod_common
        use mod_flowarray, only: u, v, w, nusgs1, rhs1
        implicit none
        integer(8) :: i, j, k
        real(8) :: rhsu2, rhsu3, rhsv1, rhsv3, rhsw1, rhsw2

!$OMP PARALLEL DO&
!$OMP PRIVATE(RHSU2,RHSU3)
        do k = 1, n3m
          do j = 1, n2m
            do i = i_bgpx, n1m
              rhsu2 = f2fyi(j) * c2cxi(i) &
                      * (nusgs1(i, j + 1, k, 3) * (v(i, j + 1, k) - v(i - 1, j + 1, k)) &
                         - nusgs1(i, j, k, 3) * (v(i, j, k) - v(i - 1, j, k)))
              rhsu3 = f2fzi(k) * c2cxi(i) &
                      * (nusgs1(i, j, k + 1, 2) * (w(i, j, k + 1) - w(i - 1, j, k + 1)) &
                         - nusgs1(i, j, k, 2) * (w(i, j, k) - w(i - 1, j, k)))
              rhs1(i, j, k, 1) = rhsu2 + rhsu3
            end do
          end do
        end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO&
!$OMP PRIVATE(RHSV1,RHSV3)
        do k = 1, n3m
          do j = j_bgpy, n2m
            do i = 1, n1m
              rhsv1 = f2fxi(i) * c2cyi(j) &
                      * (nusgs1(i + 1, j, k, 3) * (u(i + 1, j, k) - u(i + 1, j - 1, k)) &
                         - nusgs1(i, j, k, 3) * (u(i, j, k) - u(i, j - 1, k)))
              rhsv3 = f2fzi(k) * c2cyi(j) &
                      * (nusgs1(i, j, k + 1, 1) * (w(i, j, k + 1) - w(i, j - 1, k + 1)) &
                         - nusgs1(i, j, k, 1) * (w(i, j, k) - w(i, j - 1, k)))
              rhs1(i, j, k, 2) = rhsv1 + rhsv3
            end do
          end do
        end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO&
!$OMP PRIVATE(RHSW1,RHSW2)
        do k = k_bgpz, n3m
          do j = 1, n2m
            do i = 1, n1m
              rhsw1 = f2fxi(i) * c2czi(k) &
                      * (nusgs1(i + 1, j, k, 2) * (u(i + 1, j, k) - u(i + 1, j, k - 1)) &
                         - nusgs1(i, j, k, 2) * (u(i, j, k) - u(i, j, k - 1)))
              rhsw2 = f2fyi(j) * c2czi(k) &
                      * (nusgs1(i, j + 1, k, 1) * (v(i, j + 1, k) - v(i, j + 1, k - 1)) &
                         - nusgs1(i, j, k, 1) * (v(i, j, k) - v(i, j, k - 1)))
              rhs1(i, j, k, 3) = rhsw1 + rhsw2
            end do
          end do
        end do
!$OMP END PARALLEL DO

        return
      end subroutine rhssgs
!=======================================================================
!=======================================================================
      subroutine sgscalc_t
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
!$      USE OMP_LIB
        use mod_common
        use mod_flowarray, only: alsgs
        implicit none
        if (msub .eq. 1) then

          if (iread .eq. 0 .and. ntime .eq. 1) then
            alsgs = 0d0
          else
            if (itemdl .eq. 0) then
              if (idvmon .eq. 1) then
                call sgs_edm
              else
                call sgs_cedm
              end if
            else if (itemdl .eq. 1) then
              write (*, *) 'OTHER MODEL IS NOT IMPLEMENTED IN THIS VERSION'
              write (*, *) 'SGS.F90 LINE ***'
              stop
            end if
          end if

          call alpzero               ! SET ALSGS TO ZERO IN THE SOLID BODY
          call alpinterpol           !              INTERPOLATION OF NUSGS
        end if

        return
      end subroutine sgscalc_t
!=======================================================================
!=======================================================================
      subroutine sgs_cedm
!=======================================================================
!
!     ALSGS=CSGSHF*DEL**2.*SQRT(2.*SSS)
!     CSGSHF : EDDY DIFFUSIVITY CONSTANT (= CSGSHF FROM MOD_COMMON)
!
!     DEL    : FILTER WIDTH CHARACTERIZED BY MEAN CELL WIDTH
!     SSS    : S_{ij}*S_{ij}
!
!     S_{ij} : STRAIN RATE (0.5 * [DUj/DXi + DUi/DXj]) 
!
!-----------------------------------------------------------------------
      use mod_common
      use mod_flowarray, only : u,v,w,t,alsgs
      implicit none
      integer(8) :: i,j,k
      real(8)    :: d(9),s(6)
      real(8)    :: del,sss

!$OMP PARALLEL DO private(D,S,DEL,SSS)
      do k=1,n3m
        do j=1,n2m
          do i=1,n1m

            call velgrad(i,j,k,d,s)         

            d(1) = f2fx(i)
            d(2) = f2fy(j)
            d(3) = f2fz(k)
      
            del=(d(1)*d(2)*d(3))**(1./3.)
         
            sss=s(1)**2.+2.*s(2)**2.      &
               +2.*s(3)**2.+s(4)**2.      &
               +2.*s(5)**2.+s(6)**2.
            sss = dmax1(sss, 1.0e-16)
   
            alsgs(i,j,k)= csgshf*del**2.*sqrt(2.*sss)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

      return
      end subroutine sgs_cedm
!=======================================================================
!=======================================================================
      subroutine sgs_edm
!=======================================================================
!
!     ALSGS=CSGSHF*DEL**2.*SQRT(2.*SSS)
!           CSGSHF DETERMINED VIA GERMANO IDENTITY
!     CSGSHF : EDDY DIFFUSIVITY CONSTANT (= CSGSHF FROM MOD_COMMON)
!     DEL    : FILTER WIDTH CHARACTERIZED BY MEAN CELL WIDTH
!     SSS    : S_{IJ}*S_{IJ}
!     S_{IJ} : STRAIN RATE (0.5 * [DUJ/DXI + DUI/DXJ])
!
!-----------------------------------------------------------------------
        use mod_common
        use mod_flowarray, only: u, v, w, t, alsgs, nwall_dvm, aalp, llij, mmij, uui
        implicit none
        integer(8) :: i, j, k
        real(8) :: ljmj_v, mjmj_v
        real(8) :: del, sss, ami, volume
        real(8) :: d(9), s(6), alp(3)

        ljmj_v = 0.
        mjmj_v = 0.
        csgshf = 0.

        allocate (aalp(n1m, n2m, n3m, 6))
        allocate (llij(n1m, n2m, n3m, 3))
        allocate (mmij(n1m, n2m, n3m, 3))
        allocate (uui(n1m, n2m, n3m, 3))

!$OMP PARALLEL DO&
!$OMP PRIVATE(D,S,ALP)&
!$OMP PRIVATE(DEL,SSS)
        do k = 1, n3m
          do j = 1, n2m
            do i = 1, n1m
              call velgrad(i, j, k, d, s) ! D IS DUMMY VARIABLE
              call temgrad(i, j, k, alp)

              d(1) = f2fx(i)
              d(2) = f2fy(j)
              d(3) = f2fz(k)

              del = (d(1) * d(2) * d(3))**(1./3.)
              sss = s(1)**2.+2 * s(2)**2. &
                    +2 * s(3)**2.+s(4)**2. &
                    +2 * s(5)**2.+s(6)**2.

              sss = dmax1(sss, 1.0e-16)

              alsgs(i, j, k) = del**2.*sqrt(2.*sss)

              mmij(i, j, k, 1) = alsgs(i, j, k) * alp(1)
              mmij(i, j, k, 2) = alsgs(i, j, k) * alp(2)
              mmij(i, j, k, 3) = alsgs(i, j, k) * alp(3)
            end do
          end do
        end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
        do i = 1, 3
          call test_filter(mmij(:, :, :, i), filter)
        end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO&
!$OMP PRIVATE(D,S,ALP)&
!$OMP PRIVATE(DEL,SSS)
        do k = 1, n3m
          do j = 1, n2m
            do i = 1, n1m
              call velgrad(i, j, k, d, s) ! D IS DUMMY VARIABLE
              call temgrad(i, j, k, alp)

              d(1) = (0.5 * f2fx(i - 1) * (1.-fixiu(i)) + f2fx(i) &
                      + 0.5 * f2fx(i + 1) * (1.-fixil(i)))
              d(2) = (0.5 * f2fy(j - 1) * (1.-fixju(j)) + f2fy(j) &
                      + 0.5 * f2fy(j + 1) * (1.-fixjl(j)))
              d(3) = (0.5 * f2fz(k - 1) * (1.-fixku(k)) + f2fz(k) &
                      + 0.5 * f2fz(k + 1) * (1.-fixkl(k)))

              del = (d(1) * d(2) * d(3))**(1./3.)
              sss = s(1)**2.+2 * s(2)**2. &
                    +2 * s(3)**2.+s(4)**2. &
                    +2 * s(5)**2.+s(6)**2.

              sss = dmax1(sss, 1.0e-16)

              ami = del**2.*sqrt(2.*sss)

              mmij(i, j, k, 1) = ami * alp(1) - mmij(i, j, k, 1)
              mmij(i, j, k, 2) = ami * alp(2) - mmij(i, j, k, 2)
              mmij(i, j, k, 3) = ami * alp(3) - mmij(i, j, k, 3)
              uui(i, j, k, 1) = 0.5 * (u(i, j, k) + u(i + 1, j, k))
              uui(i, j, k, 2) = 0.5 * (v(i, j, k) + v(i, j + 1, k))
              uui(i, j, k, 3) = 0.5 * (w(i, j, k) + w(i, j, k + 1))
              aalp(i, j, k, 1) = t(i, j, k)
              aalp(i, j, k, 2) = t(i, j, k)
              aalp(i, j, k, 3) = t(i, j, k)
              llij(i, j, k, 1) = uui(i, j, k, 1) * aalp(i, j, k, 1)
              llij(i, j, k, 2) = uui(i, j, k, 2) * aalp(i, j, k, 2)
              llij(i, j, k, 3) = uui(i, j, k, 3) * aalp(i, j, k, 3)
            end do
          end do
        end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
        do i = 1, 3
          call test_filter(uui(:, :, :, i), filter)
        end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
        do i = 1, 3
          call test_filter(aalp(:, :, :, i), filter)
        end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
        do j = 1, 3
          call test_filter(llij(:, :, :, j), filter)
        end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO&
!$OMP REDUCTION(+:LJMJ_V,MJMJ_V)
        do k = 1, n3m
          do j = 1, n2m
            do i = 1, n1m
              volume = f2fx(i) * f2fy(j) * f2fz(k)
              llij(i, j, k, 1) = llij(i, j, k, 1) - uui(i, j, k, 1) * aalp(i, j, k, 1)
              llij(i, j, k, 2) = llij(i, j, k, 2) - uui(i, j, k, 2) * aalp(i, j, k, 2)
              llij(i, j, k, 3) = llij(i, j, k, 3) - uui(i, j, k, 3) * aalp(i, j, k, 3)
              ljmj_v = ljmj_v &
                       + (llij(i, j, k, 1) * mmij(i, j, k, 1) &
                          + llij(i, j, k, 2) * mmij(i, j, k, 2) &
                          + llij(i, j, k, 3) * mmij(i, j, k, 3)) &
                       * float(nwall_dvm(i, j, k)) * volume
              mjmj_v = mjmj_v &
                       + (mmij(i, j, k, 1) * mmij(i, j, k, 1) &
                          + mmij(i, j, k, 2) * mmij(i, j, k, 2) &
                          + mmij(i, j, k, 3) * mmij(i, j, k, 3)) &
                       * float(nwall_dvm(i, j, k)) * volume
            end do
          end do
        end do
!$OMP END PARALLEL DO

        deallocate (aalp)
        deallocate (llij)
        deallocate (mmij)
        deallocate (uui)

        csgshf = dmax1(0., -ljmj_v / mjmj_v)

!$OMP PARALLEL DO
        do k = 1, n3m
          do j = 1, n2m
            do i = 1, n1m
              alsgs(i, j, k) = csgshf * alsgs(i, j, k)
            end do
          end do
        end do
!$OMP END PARALLEL DO

        return
      end subroutine sgs_edm
!=======================================================================
      subroutine temgrad(i, j, k, a)
!=======================================================================
!$      USE OMP_LIB
        use mod_common
        use mod_flowarray, only: t
        implicit none
        integer(8) :: i, j, k
        real(8) :: a(3)
        integer(8) :: kplus, kminus, jplus, jminus, iplus, iminus
        real(8) :: tg1, tg2, tg3
        real(8) :: tp, tm

        kplus = kpv(k)
        kminus = kmv(k)
        jplus = jpv(j)
        jminus = jmv(j)
        iplus = ipv(i)
        iminus = imv(i)

        tp = 0.5 * (f2fx(i) * t(iplus, j, k) + f2fx(iplus) * t(i, j, k)) * c2cxi(iplus) &
             * (1.-fixiu(i)) + t(iplus, j, k) * fixiu(i)
        tm = 0.5 * (f2fx(iminus) * t(i, j, k) + f2fx(i) * t(iminus, j, k)) * c2cxi(i) &
             * (1.-fixil(i)) + t(iminus, j, k) * fixil(i)

        tg1 = f2fxi(i) * (tp - tm)

        tp = 0.5 * (f2fy(j) * t(i, jplus, k) + f2fy(jplus) * t(i, j, k)) * c2cyi(jplus) &
             * (1.-fixju(j)) + t(i, jplus, k) * fixju(j)
        tm = 0.5 * (f2fy(jminus) * t(i, j, k) + f2fy(j) * t(i, jminus, k)) * c2cyi(j) &
             * (1.-fixjl(j)) + t(i, jminus, k) * fixjl(j)

        tg2 = f2fyi(i) * (tp - tm)

        tp = 0.5 * (f2fz(k) * t(i, j, kplus) + f2fz(kplus) * t(i, j, k)) * c2czi(kplus) &
             * (1.-fixku(k)) + t(i, j, kplus) * fixku(k)
        tm = 0.5 * (f2fz(kminus) * t(i, j, k) + f2fz(k) * t(i, j, kminus)) * c2czi(k) &
             * (1.-fixkl(k)) + t(i, j, kminus) * fixkl(k)

        tg3 = f2fzi(i) * (tp - tm)

        a(1) = tg1
        a(2) = tg2
        a(3) = tg3

        return
      end subroutine temgrad
!=======================================================================
!=======================================================================
      subroutine alpzero
!=======================================================================
        use mod_common
        use mod_flowarray, only: alsgs, inz, jnz, knz
        implicit none
        integer(8) :: i, j, k
        integer(8) :: imax, jmax, kmax, l
        real(8) :: alsgsavg, alsgsmax, area
        real(8) :: funcbody

        imax = 0
        jmax = 0
        kmax = 0

!$OMP PARALLEL DO
        do l = 1, nzero
          alsgs(inz(l), jnz(l), knz(l)) = 0.
        end do
!$OMP END PARALLEL DO

        alsgsavg = 0.
        area = 0.
        alsgsmax = 0.

!$OMP PARALLEL DO&
!$OMP REDUCTION(+:ALSGSAVG, AREA)&
!$OMP REDUCTION(MAX:ALSGSMAX)
        do k = 1, n3m
          do j = 1, n2m
            do i = 1, n1m
              if (funcbody(xmp(i), ymp(j), zmp(k), time) .ge. 1.e-10) then
                alsgsavg = alsgsavg + alsgs(i, j, k) * f2fx(i) * f2fy(j) * f2fz(k)
                area = area + f2fx(i) * f2fy(j) * f2fz(k)
                if (alsgs(i, j, k) .gt. alsgsmax) then
                  alsgsmax = alsgs(i, j, k)
                  imax = i
                  jmax = j
                  kmax = k
                end if
              else
                alsgs(i, j, k) = 0.
              end if
            end do
          end do
        end do
!$OMP END PARALLEL DO

        alsgsavg = alsgsavg / area

        alpavg = alsgsavg * re * pr
        alpmax = alsgsmax * re * pr

        write (*, 99) csgshf
        write (*, 100) alpavg
        write (*, 110) alpmax, xmp(imax), ymp(jmax), zmp(kmax), imax, jmax, kmax
99      format('CSGSHF =  ', es10.3)
100     format('ALPAVG =  ', es10.3)
110     format('ALPMAX =  ', es10.3, ' @ ', 3f10.4, ' , ', 3i5)

        return
      end subroutine alpzero
!=======================================================================
!=======================================================================
      subroutine alpinterpol
!=======================================================================
        use mod_common
        use mod_flowarray, only: alsgs, alsgs1, inz, jnz, knz
        implicit none
        integer(8) :: i, j, k, im, jm, km

        alsgs1 = 0.

        if (xprdic .eq. 1) then
!$OMP PARALLEL DO
          do k = 0, n3
            do j = 0, n2
              alsgs(0, j, k) = alsgs(n1m, j, k)
              alsgs(n1, j, k) = alsgs(1, j, k)
            end do
          end do
!$OMP END PARALLEL DO
        end if

        if (yprdic .eq. 1) then
!$OMP PARALLEL DO
          do k = 0, n3
            do i = 0, n1
              alsgs(i, 0, k) = alsgs(i, n2m, k)
              alsgs(i, n2, k) = alsgs(i, 1, k)
            end do
          end do
!$OMP END PARALLEL DO
        end if

        if (zprdic .eq. 1) then
!$OMP PARALLEL DO
          do j = 0, n2
            do i = 0, n1
              alsgs(i, j, 0) = alsgs(i, j, n3m)
              alsgs(i, j, n3) = alsgs(i, j, 1)
            end do
          end do
!$OMP END PARALLEL DO
        end if

!$OMP PARALLEL DO
        do k = 1, n3m
          do j = 1, n2m
            do i = 1, n1m
              alsgs1(i, j, k, 1) = c2cyi(j) * c2czi(k) * 0.25 &
                                   * (f2fy(j - 1) * (f2fz(k - 1) * alsgs(i, j, k) + f2fz(k) * alsgs(i, j, k - 1)) &
                                      + f2fy(j) * (f2fz(k - 1) * alsgs(i, j - 1, k) + f2fz(k) * alsgs(i, j - 1, k - 1)))
            end do
          end do
        end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
        do k = 1, n3m
          do j = 1, n2m
            do i = 1, n1m
              alsgs1(i, j, k, 2) = c2czi(k) * c2cxi(i) * 0.25 &
                                   * (f2fz(k - 1) * (f2fx(i - 1) * alsgs(i, j, k) + f2fx(i) * alsgs(i - 1, j, k)) &
                                      + f2fz(k) * (f2fx(i - 1) * alsgs(i, j, k - 1) + f2fx(i) * alsgs(i - 1, j, k - 1)))
            end do
          end do
        end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
        do k = 1, n3m
          do j = 1, n2m
            do i = 1, n1m
              alsgs1(i, j, k, 3) = c2cxi(i) * c2cyi(j) * 0.25 &
                                   * (f2fx(i - 1) * (f2fy(j - 1) * alsgs(i, j, k) + f2fy(j) * alsgs(i, j - 1, k)) &
                                      + f2fx(i) * (f2fy(j - 1) * alsgs(i - 1, j, k) + f2fy(j) * alsgs(i - 1, j - 1, k)))
            end do
          end do
        end do
!$OMP END PARALLEL DO

!!$OMP PARALLEL DO PRIVATE(JM,KM)
        do j = 1, n2m
          do k = 1, n3m
            jm = jmv(j)
            km = kmv(k)
            if (xprdic .eq. 1) then
              alsgs1(n1, j, k, 1) = alsgs1(1, j, k, 1)
              alsgs1(n1, j, k, 2) = alsgs1(1, j, k, 2)
              alsgs1(n1, j, k, 3) = alsgs1(1, j, k, 3)
            else
              if (bc_xbtm .eq. 0) then
                alsgs1(1, j, k, 1) = 0.
                alsgs1(1, j, k, 2) = 0.
                alsgs1(1, j, k, 3) = 0.
              else
                alsgs1(1, j, k, 2) = &
                  c2czi(k) * 0.5 * (f2fz(k) * alsgs(2, j, km) + f2fz(km) * alsgs(2, j, k))
                alsgs1(1, j, k, 3) = &
                  c2cyi(j) * 0.5 * (f2fy(j) * alsgs(2, jm, k) + f2fy(jm) * alsgs(2, j, k))
              end if
              if (bc_xtop .eq. 0) then
                alsgs1(n1, j, k, 1) = 0.
                alsgs1(n1, j, k, 2) = 0.
                alsgs1(n1, j, k, 3) = 0.
              else
                alsgs1(n1, j, k, 2) = &
                  c2czi(k) * 0.5 * (f2fz(k) * alsgs(n1m, j, km) + f2fz(km) * alsgs(n1m, j, k))
                alsgs1(n1, j, k, 3) = &
                  c2cyi(j) * 0.5 * (f2fy(j) * alsgs(n1m, jm, k) + f2fy(jm) * alsgs(n1m, j, k))
              end if
            end if
          end do
        end do
!!$OMP END PARALLEL DO

!!$OMP PARALLEL DO PRIVATE(IM,KM)
        do k = 1, n3m
          do i = 1, n1m
            im = imv(i)
            km = kmv(k)
            if (yprdic .eq. 1) then
              alsgs1(i, n2, k, 1) = alsgs1(i, 1, k, 1)
              alsgs1(i, n2, k, 2) = alsgs1(i, 1, k, 2)
              alsgs1(i, n2, k, 3) = alsgs1(i, 1, k, 3)
            else
              if (bc_ybtm .eq. 0) then
                alsgs1(i, 1, k, 1) = 0.
                alsgs1(i, 1, k, 2) = 0.
                alsgs1(i, 1, k, 3) = 0.
              else
                alsgs1(i, 1, k, 2) = &
                  c2czi(k) * 0.5 * (f2fz(k) * alsgs(i, 2, km) + f2fz(km) * alsgs(i, 2, k))
                alsgs1(i, 1, k, 3) = &
                  c2cxi(i) * 0.5 * (f2fx(i) * alsgs(im, 2, k) + f2fx(im) * alsgs(i, 2, k))
              end if
              if (bc_ytop .eq. 0) then
                alsgs1(i, n2, k, 1) = 0.
                alsgs1(i, n2, k, 2) = 0.
                alsgs1(i, n2, k, 3) = 0.
              else
                alsgs1(i, n2, k, 2) = &
                  c2czi(k) * 0.5 * (f2fz(k) * alsgs(i, n2m, km) + f2fz(km) * alsgs(i, n2m, k))
                alsgs1(i, n2, k, 3) = &
                  c2cxi(i) * 0.5 * (f2fx(i) * alsgs(im, n2m, k) + f2fx(im) * alsgs(i, n2m, k))
              end if
            end if
          end do
        end do
!!$OMP END PARALLEL DO

!!$OMP PARALLEL DO PRIVATE(IM,JM)
        do i = 1, n1m
          do j = 1, n2m
            im = imv(i)
            jm = jmv(j)
            if (zprdic .eq. 1) then
              alsgs1(i, j, n3, 1) = alsgs1(i, j, 1, 1)
              alsgs1(i, j, n3, 2) = alsgs1(i, j, 1, 2)
              alsgs1(i, j, n3, 3) = alsgs1(i, j, 1, 3)
            else
              if (bc_zbtm .eq. 0) then
                alsgs1(i, j, 1, 1) = 0.
                alsgs1(i, j, 1, 2) = 0.
                alsgs1(i, j, 1, 3) = 0.
              else
                alsgs1(i, j, 1, 2) = &
                  c2cyi(j) * 0.5 * (f2fy(j) * alsgs(i, jm, 2) + f2fy(jm) * alsgs(i, j, 2))
                alsgs1(i, j, 1, 3) = &
                  c2cxi(i) * 0.5 * (f2fx(i) * alsgs(im, j, 2) + f2fx(im) * alsgs(i, j, 2))
              end if
              if (bc_ztop .eq. 0) then
                alsgs1(i, j, n3, 1) = 0.
                alsgs1(i, j, n3, 2) = 0.
                alsgs1(i, j, n3, 3) = 0.
              else
                alsgs1(i, j, n3, 2) = &
                  c2cyi(k) * 0.5 * (f2fy(j) * alsgs(i, jm, n3m) + f2fy(jm) * alsgs(i, j, n3m))
                alsgs1(i, j, n3, 3) = &
                  c2cxi(i) * 0.5 * (f2fx(i) * alsgs(im, j, n3m) + f2fx(im) * alsgs(i, j, n3m))
              end if
            end if
          end do
        end do
!!$OMP END PARALLEL DO

        return
      end subroutine alpinterpol
!=======================================================================
!=======================================================================
      subroutine rhssgs_t
!=======================================================================
!
!     CALCULATE NON-LINEAR SGS HF TERMS
!
!-----------------------------------------------------------------------
        use mod_common
        use mod_flowarray, only: u, v, w, t, alsgs1, rhs1
        implicit none
        integer(8) :: i, j, k
        real(8) :: rhst1, rhst2, rhst3

!$OMP PARALLEL DO&
!$OMP PRIVATE(RHST1,RHST2,RHST3)
        do k = 1, n3m
          do j = 1, n2m
            do i = 1, n1m
              rhst1 = 0.
              rhst2 = 0.
              rhst3 = 0.
              rhs1(i, j, k, 4) = rhst1 + rhst2 + rhst3
            end do
          end do
        end do
!$OMP END PARALLEL DO

        return
      end subroutine rhssgs_t
!=======================================================================
