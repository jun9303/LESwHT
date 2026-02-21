!=======================================================================
!
!     MAIN SOLVER OF POISSON EQUATION
!
!     AX=B
!     A: COEFFICIENT OF DISCRETIZED POISSON EQUATION
!     X: PHI: PSEUDO PRESSURE, OUTPUT OF SUBROUTINE POISSON
!     B: DIVGSUM: OUTPUT OF SUBROUTINE DIVGS.
!
!     X-DIRECTION: MULTI-GRID ITERATION/GSOR METHOD
!     Y-DIRECTION: TDMA
!     Z-DIRECTION: FOURIER TRANSFORM
!
!     JUN. 2017, J. PARK
!     FEB. 2026, S. LEE (SEVERAL BUG FIXES)
!
!=======================================================================
      subroutine poisinit
!=======================================================================
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, k, ilev

        call z_ft_allo

        allocate (ac(n1md, n2m, n3mh))
        allocate (gam(n1md, n2m, n3mh))
        allocate (bet(n1md, n2m, n3mh))

        call pmat
!------MULTIGRID METHOD
        call coefmg

        do ilev = 0, nlev
          do k = 1, n3mh
            do i = iimg(ilev, 1), iimg(ilev, 2)
              bet(i, 1, k) = 1.0_8 / ac(i, 1, k)
            end do
            do j = 2, n2m
              do i = iimg(ilev, 1), iimg(ilev, 2)
                gam(i, j, k) = an(j - 1) * bet(i, j - 1, k)
                bet(i, j, k) = 1.0_8 / (ac(i, j, k) - as(j) * gam(i, j, k))
              end do
            end do
          end do
        end do

        open (77, file='../output/ftr/mgftresiduemax.dat')

        return
      end subroutine poisinit

!=======================================================================
      subroutine pmat
!=======================================================================
! --- CONSTRUCT MATRIX FOR POISSON EQ. ---------------------------------
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, k, ip, jp

!-----FOR TOP LEVEL
        do i = 1, n1m
          ip = ipv(i)
          aw(i) = (1.0_8 - fixil(i)) * c2cxi(i) * f2fxi(i)
          ae(i) = (1.0_8 - fixiu(i)) * c2cxi(ip) * f2fxi(i)
        end do

        do j = 1, n2m
          jp = jpv(j)
          as(j) = (1.0_8 - fixjl(j)) * c2cyi(j) * f2fyi(j)
          an(j) = (1.0_8 - fixju(j)) * c2cyi(jp) * f2fyi(j)
        end do

        do j = 1, n2m
          do i = 1, n1m
            ac(i, j, 1) = -1.0_8 * (aw(i) + ae(i) + as(j) + an(j))
          end do
        end do

        n3mh = n3m / 2 + 1
        if (n3m > 1) call mwavenumber     ! INIT. MODIFIED WAVE #.

        do k = 2, n3mh
          do j = 1, n2m
            do i = 1, n1m
              ac(i, j, k) = ac(i, j, 1) - ak3(k)
            end do
          end do
        end do

        return
      end subroutine pmat

!=======================================================================
      subroutine mwavenumber
!=======================================================================
! --- MODIFIED WAVE NUMBER DEFINITION
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: k
        real(8) :: pi

        pi = acos(-1.0_8)

        do k = 1, n3mh
          ak3(k) = 2.0_8 * (1.0_8 - cos(2.0_8 * float(k - 1) * pi / float(n3m))) * f2fzi(1) * f2fzi(1)
        end do

        return
      end subroutine mwavenumber

!=======================================================================
      subroutine poisson(phi, divgsum)
!=======================================================================
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, k, kkk
        real(8) :: phi(0:n1, 0:n2, 0:n3), divgsum(0:n1, 0:n2, 0:n3)
        complex(8) :: ccap(n1, n2, n3mh)
        complex(8), dimension(:, :), allocatable :: xxx, cp, tmp
        complex(8), dimension(:), allocatable :: xxxx, xxxx_b

        real(8) :: test

        if (n3m == 1) then
          allocate (cp(0:n1, 0:n2), tmp(n1, n2))
          cp = 0.0_8
!$OMP PARALLEL DO
          do j = 1, n2m
            do i = 1, n1m
              tmp(i, j) = divgsum(i, j, 1)
            end do
          end do

          call mg2d(cp, tmp, 1_8, test1, 0.0_8)

!$OMP PARALLEL DO
          do j = 1, n2m
            do i = 1, n1m
              phi(i, j, 1) = real(cp(i, j), 8)
            end do
          end do
          deallocate (cp, tmp)
          return
        end if

! --- DO THE FORWARD FFT
!$OMP PARALLEL PRIVATE(XXX, XXXX, XXXX_B)
        allocate (xxx(n3m, n1m))
        allocate (xxxx(n3m), xxxx_b(n3m * 2))
        call zfft1d(xxxx, n3m, 0, xxxx_b)
!$OMP DO
        do j = 1, n2m
          do k = 1, n3m
            do i = 1, n1m
              xxx(k, i) = divgsum(i, j, k)
            end do
          end do

          do i = 1, n1m
            call zfft1d(xxx(1, i), n3m, -1, xxxx_b)
          end do

          do k = 1, n3mh
            do i = 1, n1m
              ccap(i, j, k) = xxx(k, i)
            end do
          end do
        end do
!$OMP END DO
        deallocate (xxx, xxxx, xxxx_b)
!$OMP END PARALLEL

! --- SOLVE A SET OF POISSON EQS.
        test = test1 / float(n3mh) * float(n3m) * 0.8_8   ! CONVERGENCE CRITERIA

!$OMP PARALLEL PRIVATE(CP)
        allocate (cp(0:n1, 0:n2))
!$OMP DO
        do kkk = 1, n3mh
          cp = 0.0_8
          if (kkk <= imgsor) then
            call mg2d(cp, ccap(1, 1, kkk), kkk, test, 0.0_8)
          else
            call gsor2d(cp, ccap(1, 1, kkk), kkk, test, 0.0_8)
          end if

          do j = 1, n2m
            do i = 1, n1m
              ccap(i, j, kkk) = cp(i, j)
            end do
          end do
        end do
!$OMP END DO
        deallocate (cp)
!$OMP END PARALLEL

! --- DO THE INVERSE FFT
!$OMP PARALLEL PRIVATE(XXX, XXXX, XXXX_B)
        allocate (xxx(n3m, n1m))
        allocate (xxxx(n3m), xxxx_b(n3m * 2))
        call zfft1d(xxxx, n3m, 0, xxxx_b)
!$OMP DO
        do j = 1, n2m
          do k = 1, n3mh
            do i = 1, n1m
              xxx(k, i) = ccap(i, j, k)
            end do
          end do

          do k = n3mh + 1, n3m
            do i = 1, n1m
              xxx(k, i) = conjg(ccap(i, j, n3m + 2 - k))
            end do
          end do

          do i = 1, n1m
            call zfft1d(xxx(1, i), n3m, 1, xxxx_b)
          end do

          do k = 1, n3m
            do i = 1, n1m
              phi(i, j, k) = real(xxx(k, i), 8)
            end do
          end do
        end do
!$OMP END DO
        deallocate (xxx, xxxx, xxxx_b)
!$OMP END PARALLEL

        return
      end subroutine poisson

!=======================================================================
      subroutine mg2d(pc, rhs, kv, test, oldv)
!=======================================================================
!     MULTIGRID ENVIRONMENT VARIABLES
!     PC  : SOLUTION OF THE MG2D.
!     RHS : RHS OF THE POISSON EQUATION FOR EACH WAVENUMBER
!     KV  : WAVENUMBER INDEX
!     TEST    : CONDITION FOR CONVERGENCE
!     OLDV    : 0; INITALIZE TO ZERO, 1; USE PREVIOUS SOLUTION
!     NLEV    : TOTAL CELS = (MINROW1+MINROW2)*(2**NLEV)
!     LEVHALF : HALF OF NLEV
!     IWC     : FLAG FOR USING W-CYCLE ( 1 IF YOU WANT TO USE W-CYCLE)
!     NBLI    : NUMBER OF BASE LEVEL ITERATION
!     MGITR   : NUMBER OF MAXIMUM ITERATIONS
!     IIMG(ILEV,1)  : START POINT OF INDEX AT ILEV LEVEL
!     IIMG(ILEV,2)  : END POINT OF INDEX AT ILEV LEVEL
!     IIMG(ILEV,3)  : NUMBER OF ROW AT ILEV LEVEL
!-----------------------------------------------------------------------
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, kv, kk, ilev
        real(8) :: test, sumres, oldv
        complex(8) :: pc(0:n1, 0:n2), rhs(n1, n2)
        complex(8) :: resd(n1md, n2m), ggii(0:n1md, 0:n2)

        kk = 0
        resd = 0.0_8
        ggii = 0.0_8

        call toplevel(0_8, pc, rhs, sumres, oldv, test, kv, resd, ggii)
        if (sumres < test) goto 205

        do kk = 1, mgitr            ! MAIN ITERATION
          do ilev = nlev - 1, 1, -1
            call relax(ilev, 0.0_8, 1_8, kv, resd, ggii)
            call godown(ilev, kv, resd, ggii)
          end do

          call relax(0_8, 0.0_8, nbli, kv, resd, ggii)

          do ilev = 0, nlev - 2
            call goup(ilev, resd, ggii)
            call relax(ilev + 1, 1.0_8, 1_8, kv, resd, ggii)
          end do

          call toplevel(1_8, pc, rhs, sumres, 1.0_8, test, kv, resd, ggii)

          if (sumres < test) goto 205
        end do

        print *, 'ITERATION LIMIT EXCEEDED.'

205     continue
        if (kv <= 3) then
          write (77, 999) kv, kk, sumres * dtconst
        end if

999     format('KV=', i4, 3x, 'KK=', i4, 3x, 'RM=', es24.16)

        return
      end subroutine mg2d

!=======================================================================
      subroutine toplevel(id, pc, rhs, sumres, oldv, test, kv, resd, ggii)    ! ZEBRA VERSION
!=======================================================================
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, ic
        integer(8) :: id, kv
        real(8) :: sumres, test, oldv
        complex(8) :: tt
        complex(8) :: pc(0:n1, 0:n2), rhs(n1, n2)
        complex(8) :: resd(n1md, n2m), ggii(0:n1md, 0:n2)

        if (id == 1) then
!         INTERPOLATE & ADD
          ic = 0
          do i = iimg(nlev - 1, 1), iimg(nlev - 1, 2) - 1
            ic = ic + 2
            do j = 1, n2m
              pc(ic, j) = pc(ic, j) + coi1(i) * ggii(i, j) + coi2(i) * ggii(ipm(i), j)
            end do
          end do

          i = iimg(nlev - 1, 2)
          ic = ic + 2
          do j = 1, n2m
            pc(ic, j) = pc(ic, j) + coi1(i) * ggii(i, j)      ! USE ONE POINT COI(~)=1.
          end do
        end if

!  RELAX
        do j = 1, n2m
          do i = 1, n1m, 2
            ggii(i, j) = rhs(i, j) - oldv * (aw(i) * pc(imm(i), j) + ae(i) * pc(ipm(i), j))
          end do
        end do

        call trdiag1m(ggii, ggii, 1_8, n1m, 1_8, n2m, kv)

        if (kv == 1) then
          tt = ggii(1, n2m - 1)
        else
          tt = 0.0_8
        end if

        do j = 1, n2m
          do i = 2, n1m, 2
            ggii(i - 1, j) = ggii(i - 1, j) - tt
            ggii(i, j) = rhs(i, j) - aw(i) * ggii(imm(i), j) - ae(i) * (ggii(ipm(i), j) - tt)
          end do
        end do

        call trdiag1m(ggii, ggii, 2_8, n1m, 1_8, n2m, kv)

!  CALCULATE RESIDUAL
        sumres = 0.0_8

        do j = 1, n2m
          do i = 1, n1m
            resd(i, j) = rhs(i, j) - aw(i) * ggii(imm(i), j) - ae(i) * ggii(ipm(i), j) &
                         - as(j) * ggii(i, j - 1) - an(j) * ggii(i, j + 1) &
                         - ac(i, j, kv) * ggii(i, j)
            sumres = max(sumres, abs(resd(i, j)))
            pc(i, j) = ggii(i, j)
          end do
        end do

        if (sumres < test .or. id == 2) return

!  RESTRICT
        ic = -1
        do i = iimg(nlev - 1, 1), iimg(nlev - 1, 2)
          ic = ic + 2
          do j = 1, n2m
            resd(i, j) = resd(ic, j) * cor1(i) + resd(ic + 1, j) * cor2(i)
          end do
        end do

        return
      end subroutine toplevel

!=======================================================================
      subroutine trdiag1m(rr, uu, l1, l2, ll1, ll2, kv)
!=======================================================================
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, l1, l2, ll1, ll2, kv
        complex(8) :: rr(0:n1md, 0:n2), uu(0:n1md, 0:n2)

        do i = l1, l2, 2
          uu(i, ll1) = rr(i, ll1) * bet(i, 1, kv)
        end do

        do j = ll1 + 1, ll2
          do i = l1, l2, 2
            uu(i, j) = (rr(i, j) - as(j) * uu(i, j - 1)) * bet(i, j, kv)
          end do
        end do

        do j = ll2 - 1, ll1, -1
          do i = l1, l2, 2
            uu(i, j) = uu(i, j) - gam(i, j + 1, kv) * uu(i, j + 1)
          end do
        end do

        return
      end subroutine trdiag1m

!=======================================================================
      subroutine goup(ilev, resd, ggii)! INTERPOLATE RESIDUAL & ADD IT TO HIGH LEVEL
!=======================================================================
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, ilev, ibgh
        complex(8) :: resd(n1md, n2m), ggii(0:n1md, 0:n2)

        ibgh = iimg(ilev + 1, 1) - 1

        do i = iimg(ilev, 1), iimg(ilev, 2) - 1
          ibgh = ibgh + 2
          do j = 1, n2m
            ggii(ibgh, j) = ggii(ibgh, j) + coi1(i) * ggii(i, j) + coi2(i) * ggii(ipm(i), j)
          end do
        end do

        i = iimg(ilev, 2)
        ibgh = ibgh + 2
        do j = 1, n2m
          ggii(ibgh, j) = ggii(ibgh, j) + coi1(i) * ggii(i, j)           ! USE ONE POINT
        end do

        return
      end subroutine goup

!=======================================================================
      subroutine godown(ilev, kv, resd, ggii)        ! COMPUTE RESIDUAL & RESTRICT IT
!=======================================================================
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, ilev, kv, ibg
        complex(8) :: resd(n1md, n2m), ggii(0:n1md, 0:n2)

        do j = 1, n2m
          do i = iimg(ilev, 1), iimg(ilev, 2)
            resd(i, j) = resd(i, j) - aw(i) * ggii(imm(i), j) - ae(i) * ggii(ipm(i), j) &
                         - as(j) * ggii(i, j - 1) - an(j) * ggii(i, j + 1) - ac(i, j, kv) * ggii(i, j)
          end do
        end do

        ibg = iimg(ilev, 1) - 2

        do i = iimg(ilev - 1, 1), iimg(ilev - 1, 2)
          ibg = ibg + 2
          do j = 1, n2m
            resd(i, j) = resd(ibg, j) * cor1(i) + resd(ibg + 1, j) * cor2(i)
          end do
        end do

        return
      end subroutine godown

!=======================================================================
      subroutine relax(ilev, oldv, iiter, kv, resd, ggii)   ! ZEBRA VERSION
!=======================================================================
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, ilev, iiter, kv, ii
        complex(8) :: resd(n1md, n2m), ggii(0:n1md, 0:n2)
        real(8) :: oldv

        do j = 1, n2m
          do i = iimg(ilev, 1), iimg(ilev, 2), 2
            ggii(i, j) = resd(i, j) - oldv * (aw(i) * ggii(imm(i), j) + ae(i) * ggii(ipm(i), j))
          end do
        end do

        call trdiag1m(ggii, ggii, iimg(ilev, 1), iimg(ilev, 2), 1_8, n2m, kv)

        do j = 1, n2m
          do i = iimg(ilev, 1) + 1, iimg(ilev, 2), 2
            ggii(i, j) = resd(i, j) - aw(i) * ggii(imm(i), j) - ae(i) * ggii(ipm(i), j)
          end do
        end do

        call trdiag1m(ggii, ggii, iimg(ilev, 1) + 1, iimg(ilev, 2), 1_8, n2m, kv)

        do ii = 1, iiter - 1

          do j = 1, n2m
            do i = iimg(ilev, 1), iimg(ilev, 2), 2
              ggii(i, j) = resd(i, j) - (aw(i) * ggii(imm(i), j) + ae(i) * ggii(ipm(i), j))
            end do
          end do

          call trdiag1m(ggii, ggii, iimg(ilev, 1), iimg(ilev, 2), 1_8, n2m, kv)

          do j = 1, n2m
            do i = iimg(ilev, 1) + 1, iimg(ilev, 2), 2
              ggii(i, j) = resd(i, j) - aw(i) * ggii(imm(i), j) - ae(i) * ggii(ipm(i), j)
            end do
          end do

          call trdiag1m(ggii, ggii, iimg(ilev, 1) + 1, iimg(ilev, 2), 1_8, n2m, kv)

        end do
        return
      end subroutine relax
!=======================================================================
      subroutine gsor2d(u, rhs, kv, test, oldv)    ! 1 EQ. TYPE
!=======================================================================
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, kv
        real(8) :: test, oldv
        complex(8) :: u(0:n1, 0:n2), rhs(n1, n2)
        complex(8) :: ggii(0:n1md, 0:n2)
        integer(8) :: kk
        real(8) :: ww, ww2, errmax

        ww = wwsor
        ww2 = 1.0_8 - ww
        kk = 0

!  HALF RELAX
!  ----------

        do j = 1, n2m
          do i = 1, n1m, 2
            ggii(i, j) = rhs(i, j) - oldv * (aw(i) * u(imm(i), j) + ae(i) * u(ipm(i), j))
          end do
        end do

        call trdiag1m(ggii, ggii, 1_8, n1m, 1_8, n2m, kv)

        do j = 1, n2m
          do i = 1, n1m, 2
            u(i, j) = ww * ggii(i, j) + oldv * ww2 * u(i, j)
          end do
        end do

!  ANOTHER HALF
!  ------------

        do j = 1, n2m
          do i = 2, n1m, 2
            ggii(i, j) = rhs(i, j) - (aw(i) * u(imm(i), j) + ae(i) * u(ipm(i), j))
          end do
        end do

        call trdiag1m(ggii, ggii, 2_8, n1m, 1_8, n2m, kv)

        do j = 1, n2m
          do i = 2, n1m, 2
            u(i, j) = ww * ggii(i, j) + oldv * ww2 * u(i, j)
          end do
        end do

        call resid3(u, rhs, kv, errmax)
        if (errmax < test) goto 88

!  MAIN ITERATION
!  ==============
        do kk = 1, mgitr

!  HALF RELAX
!  ----------

          do j = 1, n2m
            do i = 1, n1m, 2
              ggii(i, j) = rhs(i, j) - (aw(i) * u(imm(i), j) + ae(i) * u(ipm(i), j))
            end do
          end do

          call trdiag1m(ggii, ggii, 1_8, n1m, 1_8, n2m, kv)

          do j = 1, n2m
            do i = 1, n1m, 2
              u(i, j) = ww * ggii(i, j) + ww2 * u(i, j)
            end do
          end do

!  ANOTHER HALF
!  ------------

          do j = 1, n2m
            do i = 2, n1m, 2
              ggii(i, j) = rhs(i, j) - (aw(i) * u(imm(i), j) + ae(i) * u(ipm(i), j))
            end do
          end do

          call trdiag1m(ggii, ggii, 2_8, n1m, 1_8, n2m, kv)

          do j = 1, n2m
            do i = 2, n1m, 2
              u(i, j) = ww * ggii(i, j) + ww2 * u(i, j)
            end do
          end do

          call resid3(u, rhs, kv, errmax)
          write (77, 999) kv, kk, errmax * dtconst * float(n3mh)
          if (errmax < test) goto 88

        end do

        print *, 'ITERATION LIMIT EXCEEDED.'
88      continue

999     format('KV=', i3, 3x, 'KK=', i3, 3x, 'RG=', es24.16)

        return
      end subroutine gsor2d

!=======================================================================
      subroutine resid3(u, rhs, kv, errmax)
!=======================================================================
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, kv
        real(8) :: errmax
        complex(8) :: u(0:n1, 0:n2), rhs(n1, n2)
        complex(8) :: err

        errmax = 0.0_8

        do j = 1, n2m
          do i = 1, n1m
            err = rhs(i, j) - aw(i) * u(imm(i), j) - ae(i) * u(ipm(i), j) &
                  - as(j) * u(i, j - 1) - an(j) * u(i, j + 1) &
                  - ac(i, j, kv) * u(i, j)
            errmax = max(errmax, abs(err))
          end do
        end do

        return
      end subroutine resid3
!=======================================================================
      subroutine coefmg
!=======================================================================
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, k, ilev
        integer(8) :: minrow, ibg, iend, isp, ic, ibgh
        real(8) :: iwz(n1md), iez(n1md)
        real(8) :: xmpm(0:n1md)
        real(8) :: vdx_ibg, vdx_iend, sdx_ibg, sdx_iend

        minrow = n1m / (2**nlev)
        xmpm = 0.0_8

        levhalf = nint(nlev / 2.0_8)
        iimg(nlev, 1) = 1
        iimg(nlev, 2) = n1m
        iimg(nlev, 3) = n1m

        do i = 1, n1m
          xmpm(i) = xmp(i)
          iwz(i) = 1 - 1 / i                  ! 0 ONLY IF I=1
          iez(i) = 1 - i / n1m                ! 0 ONLY IF I=N1M
        end do

!  COMPUTE FOR LOWER LEVELS
        do ilev = nlev - 1, 0, -1

          iimg(ilev, 1) = iimg(ilev + 1, 1) + iimg(ilev + 1, 3)
          iimg(ilev, 3) = minrow * (2**ilev)
          iimg(ilev, 2) = iimg(ilev, 1) + iimg(ilev, 3) - 1

          ibg = iimg(ilev, 1)
          iend = iimg(ilev, 2)
          isp = 2**(nlev - ilev)           ! WIDTH OF ONE CELL AT LOW LEVEL

          ic = 0
          do i = ibg, iend
            ic = ic + 1
            xmpm(i) = 0.5_8 * (x(ic * isp + 1) + x((ic - 1) * isp + 1))
            iwz(i) = 1 - ibg / i                ! 0 ONLY IF I=IBG
            iez(i) = 1 - i / iend               ! 0 ONLY IF I=IEND
          end do

          ic = 0
          do i = ibg, iend
            ic = ic + 1
            aw(i) = iwz(i) / ((xmpm(i) - xmpm(i - 1)) * (x(ic * isp + 1) - x((ic - 1) * isp + 1)))
            ae(i) = iez(i) / ((xmpm(i + 1) - xmpm(i)) * (x(ic * isp + 1) - x((ic - 1) * isp + 1)))
          end do

          do j = 1, n2m
            do i = iimg(ilev, 1), iimg(ilev, 2)
              ac(i, j, 1) = -1.0_8 * (aw(i) + ae(i) + as(j) + an(j))
            end do
          end do

          do k = 2, n3mh
            do j = 1, n2m
              do i = iimg(ilev, 1), iimg(ilev, 2)
                ac(i, j, k) = ac(i, j, 1) - ak3(k)
              end do
            end do
          end do

        end do

!  CALCULATE RESTRICTION COEFFS
        do ilev = nlev, 1, -1
          ibgh = iimg(ilev, 1)
          do i = iimg(ilev - 1, 1), iimg(ilev - 1, 2)
            cor1(i) = (xmpm(ibgh + 1) - xmpm(i)) / (xmpm(ibgh + 1) - xmpm(ibgh))
            cor2(i) = 1.0_8 - cor1(i)
            ibgh = ibgh + 2
          end do
        end do

!  CALCULATE INTERPOLATION COEFFS
        do ilev = 0, nlev - 1
          ibgh = iimg(ilev + 1, 1) + 1
          do i = iimg(ilev, 1), iimg(ilev, 2) - 1
            coi1(i) = (xmpm(i + 1) - xmpm(ibgh)) / (xmpm(i + 1) - xmpm(i))  ! * LOWER VALUE
            coi2(i) = 1.0_8 - coi1(i)
            ibgh = ibgh + 2
          end do
          i = iimg(ilev, 2)
          coi1(i) = 1.0_8                 ! USE ONLY ONE LOWER POINT AT UPPER WALL
        end do

!===== FOR THE X PERIODICITY
!       INTRODUCE IPM & IMM
        ! ENSURE SAFE DEFAULT INITIALIZATION FOR NON-PERIODIC BOUNDARIES AS WELL
        do ilev = nlev, 0, -1
          ibg = iimg(ilev, 1)
          iend = iimg(ilev, 2)
          do i = ibg, iend
            ipm(i) = i + 1
            imm(i) = i - 1
          end do
        end do

        if (xprdic == 1) then
          do ilev = nlev, 0, -1
            ibg = iimg(ilev, 1)
            iend = iimg(ilev, 2)
            ipm(iend) = ibg
            imm(ibg) = iend
          end do

          do ilev = nlev - 1, 0, -1
            ibg = iimg(ilev, 1)
            iend = iimg(ilev, 2)
            isp = 2**(nlev - ilev)
            vdx_ibg = xmpm(ibg) - x(1) + x(n1) - xmpm(iend)
            vdx_iend = xmpm(ibg) - x(1) + x(n1) - xmpm(iend)
            sdx_ibg = x(1 + isp) - x(1)
            sdx_iend = x(n1) - x(n1 - isp)
            aw(ibg) = 1.0_8 / (vdx_ibg * sdx_ibg)
            ae(iend) = 1.0_8 / (vdx_iend * sdx_iend)
          end do

          do ilev = nlev - 1, 0, -1
            do j = 1, n2m
              do i = iimg(ilev, 1), iimg(ilev, 2)
                ac(i, j, 1) = -1.0_8 * (aw(i) + ae(i) + as(j) + an(j))
              end do
            end do

            do k = 2, n3mh
              do j = 1, n2m
                do i = iimg(ilev, 1), iimg(ilev, 2)
                  ac(i, j, k) = ac(i, j, 1) - ak3(k)
                end do
              end do
            end do
          end do

!  CALCULATE INTERPOLATION COEFFS
          do ilev = 0, nlev - 1
            ibg = iimg(ilev, 1)
            iend = iimg(ilev, 2)
            vdx_iend = xmpm(ibg) - x(1) + x(n1) - xmpm(iend)
            coi2(iend) = (xmpm(iimg(ilev + 1, 2)) - xmpm(iend)) / vdx_iend
            coi1(iend) = 1.0_8 - coi2(iend)
          end do

        end if

        do ilev = nlev, 0, -1
          write (*, *) 'IIMG(1', ilev, ')=', iimg(ilev, 1)
          write (*, *) 'IIMG(2', ilev, ')=', iimg(ilev, 2)
          write (*, *) 'IIMG(3', ilev, ')=', iimg(ilev, 3)
        end do

        return
      end subroutine coefmg
