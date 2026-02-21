!=======================================================================
!
!     MAIN SOLVER OF POISSON EQUATION
!
!     AX=B
!     A: COEFFICIENT OF DISCRETIZED POISSON EQUATION
!     X: PHI: PSEUDO PRESSURE, OUTPUT OF SUBROUTINE POISSON
!     B: DIVGSUM: OUTPUT OF SUBROUTINE DIVGS.
!
!     X-DIRECTION: FOURIER TRANSFORM
!     Y-DIRECTION: TDMA
!     Z-DIRECTION: MULTI-GRID ITERATION/GSOR METHOD
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

        call x_ft_allo

        allocate (ac(n3md, n2m, n1mh))
        allocate (gam(n3md, n2m, n1mh))
        allocate (bet(n3md, n2m, n1mh))

        call pmat
!------MULTIGRID METHOD
        call coefmg

        do ilev = 0, nlev
          do i = 1, n1mh
            do k = kkmg(ilev, 1), kkmg(ilev, 2)
              bet(k, 1, i) = 1.0_8 / ac(k, 1, i)
            end do
            do j = 2, n2m
              do k = kkmg(ilev, 1), kkmg(ilev, 2)
                gam(k, j, i) = an(j - 1) * bet(k, j - 1, i)
                bet(k, j, i) = 1.0_8 / (ac(k, j, i) - as(j) * gam(k, j, i))
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
        integer(8) :: i, j, k, jp, kp

!-----FOR TOP LEVEL
        do k = 1, n3m
          kp = kpv(k)
          ab(k) = (1.0_8 - fixkl(k)) * c2czi(k) * f2fzi(k)
          af(k) = (1.0_8 - fixku(k)) * c2czi(kp) * f2fzi(k)
        end do

        do j = 1, n2m
          jp = jpv(j)
          as(j) = (1.0_8 - fixjl(j)) * c2cyi(j) * f2fyi(j)
          an(j) = (1.0_8 - fixju(j)) * c2cyi(jp) * f2fyi(j)
        end do

        do j = 1, n2m
          do k = 1, n3m
            ac(k, j, 1) = -1.0_8 * (ab(k) + af(k) + as(j) + an(j))
          end do
        end do

        n1mh = n1m / 2 + 1

        call mwavenumber     ! INIT. MODIFIED WAVE #.

        do i = 2, n1mh
          do j = 1, n2m
            do k = 1, n3m
              ac(k, j, i) = ac(k, j, 1) - ai3(i)
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
        integer(8) :: i
        real(8) :: pi

        pi = acos(-1.0_8)

        do i = 1, n1mh
          ai3(i) = 2.0_8 * (1.0_8 - cos(2.0_8 * pi * float(i - 1) / float(n1m))) * f2fxi(1) * f2fxi(1)
        end do

        return
      end subroutine mwavenumber

!=======================================================================
      subroutine poisson(phi, divgsum)
!=======================================================================
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, k, iii
        real(8) :: phi(0:n1, 0:n2, 0:n3), divgsum(0:n1, 0:n2, 0:n3)
        complex(8) :: ccap(n3, n2, n1mh)
        complex(8), dimension(:, :), allocatable :: xxx, cp
        complex(8), dimension(:), allocatable :: xxxx, xxxx_b

        real(8) :: test, phiref

! --- DO THE FORWARD FFT
!$OMP PARALLEL PRIVATE(XXX, XXXX, XXXX_B)
        allocate (xxx(n1m, n3m))
        allocate (xxxx(n1m), xxxx_b(n1m * 2))
        call zfft1d(xxxx, n1m, 0, xxxx_b)
!$OMP DO
        do j = 1, n2m

          do i = 1, n1m
            do k = 1, n3m
              xxx(i, k) = divgsum(i, j, k)
            end do
          end do

          do k = 1, n3m
            call zfft1d(xxx(1, k), n1m, -1, xxxx_b)
          end do

          do i = 1, n1mh
            do k = 1, n3m
              ccap(k, j, i) = xxx(i, k)
            end do
          end do

        end do
!$OMP END DO
        deallocate (xxx, xxxx, xxxx_b)
!$OMP END PARALLEL

! --- SOLVE A SET OF POISSON EQS.
        test = test1 / float(n1mh) * float(n1m) * 0.9_8

!$OMP PARALLEL PRIVATE(CP)
        allocate (cp(0:n3, 0:n2))
!$OMP DO
        do iii = 1, n1mh
          cp = 0.0_8

          if (iii <= imgsor) then
            call mg2d(cp, ccap(1, 1, iii), iii, test, 0.0_8)
          else
            call gsor2d(cp, ccap(1, 1, iii), iii, test, 0.0_8)
          end if

          do j = 1, n2m
            do k = 1, n3m
              ccap(k, j, iii) = cp(k, j)
            end do
          end do
        end do
!$OMP END DO
        deallocate (cp)
!$OMP END PARALLEL

! --- DO THE INVERSE FFT
!$OMP PARALLEL PRIVATE(XXX, XXXX, XXXX_B)
        allocate (xxx(n1m, n3m))
        allocate (xxxx(n1m), xxxx_b(n1m * 2))
        call zfft1d(xxxx, n1m, 0, xxxx_b)
!$OMP DO
        do j = 1, n2m
          do i = 1, n1mh
            do k = 1, n3m
              xxx(i, k) = ccap(k, j, i)
            end do
          end do

          do i = n1mh + 1, n1m
            do k = 1, n3m
              xxx(i, k) = conjg(ccap(k, j, n1m + 2 - i))
            end do
          end do

          do k = 1, n3m
            call zfft1d(xxx(1, k), n1m, 1, xxxx_b)
          end do

          do i = 1, n1m
            do k = 1, n3m
              phi(i, j, k) = real(xxx(i, k), 8)
            end do
          end do
        end do
!$OMP END DO
        deallocate (xxx, xxxx, xxxx_b)
!$OMP END PARALLEL

        if (ich == 1) then
!     SET THE AVERAGE PHI AT THE UPPER WALL TO BE ZERO.
          phiref = 0.0_8
!$OMP PARALLEL DO REDUCTION(+:PHIREF)
          do i = 1, n1m
            phiref = phiref + phi(i, n2m, 1) * f2fx(i)
          end do
          phiref = phiref / xl

!$OMP PARALLEL DO
          do k = 1, n3m
            do j = 1, n2m
              do i = 1, n1m
                phi(i, j, k) = phi(i, j, k) - phiref
              end do
            end do
          end do
        end if

        return
      end subroutine poisson

!=======================================================================
      subroutine mg2d(pc, rhs, iv, test, oldv)
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
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: k, j, iv, ii, ilev
        real(8) :: test, sumres, oldv
        complex(8) :: pc(0:n3, 0:n2), rhs(n3, n2)
        complex(8) :: resd(n3md, n2m), ggii(0:n3md, 0:n2)

        ii = 0
        resd = 0.0_8
        ggii = 0.0_8

        call toplevel(0_8, pc, rhs, sumres, oldv, test, iv, resd, ggii)
        if (sumres < test) goto 2000

        do ii = 1, mgitr            ! MAIN ITERATION
          do ilev = nlev - 1, 1, -1
            call relax(ilev, 0.0_8, 1_8, iv, resd, ggii)
            call godown(ilev, iv, resd, ggii)
          end do

          call relax(0_8, 0.0_8, nbli, iv, resd, ggii)

          do ilev = 0, nlev - 2
            call goup(ilev, resd, ggii)
            call relax(ilev + 1, 1.0_8, 1_8, iv, resd, ggii)
          end do

          call toplevel(1_8, pc, rhs, sumres, 1.0_8, test, iv, resd, ggii)

          if (sumres < test) goto 2000
        end do

        write (*, *) 'ITERATION LIMIT EXCEEDED.'

2000    continue

        if (iv <= 2) write (77, 201) ' mg, time = ', time, iv, ii, sumres * dtconst * float(n1mh) / float(n1m)

201     format(a13, f13.5, ' IV=', i5, ' II=', i5, ' RG=', es16.8)

        return
      end subroutine mg2d
!=======================================================================
      subroutine toplevel(id, pc, rhs, sumres, oldv, test, iv, resd, ggii)    ! ZEBRA VERSION
!=======================================================================
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: k, j, kc
        integer(8) :: id, iv
        real(8) :: sumres, test, oldv
        complex(8) :: tt
        complex(8) :: pc(0:n3, 0:n2), rhs(n3, n2)
        complex(8) :: resd(n3md, n2m), ggii(0:n3md, 0:n2)

        if (id == 1) then
!         INTERPOLATE & ADD
          kc = 0
          do k = kkmg(nlev - 1, 1), kkmg(nlev - 1, 2)
            kc = kc + 2
            do j = 1, n2m
              pc(kc, j) = pc(kc, j) + coi1(k) * ggii(k, j) + coi2(k) * ggii(kpm(k), j)
            end do
          end do
        end if

!  RELAX
        do j = 1, n2m
          do k = 1, n3m, 2
            ggii(k, j) = rhs(k, j) - oldv * (ab(k) * pc(kmm(k), j) + af(k) * pc(kpm(k), j))
          end do
        end do

        call trdiag1m(ggii, ggii, 1_8, n3m, 1_8, n2m, iv)

        if (iv == 1) then
          tt = ggii(1, n2m)
          do j = 1, n2m
            do k = 1, n3m, 2
              ggii(k, j) = ggii(k, j) - tt
            end do
          end do
        end if

        do j = 1, n2m
          do k = 2, n3m, 2
            ggii(k, j) = rhs(k, j) - (ab(k) * ggii(kmm(k), j) + af(k) * (ggii(kpm(k), j)))
          end do
        end do

        call trdiag1m(ggii, ggii, 2_8, n3m, 1_8, n2m, iv)

!  CALCULATE RESIDUAL
        sumres = 0.0_8

        do j = 1, n2m
          do k = 1, n3m
            resd(k, j) = rhs(k, j) - ab(k) * ggii(kmm(k), j) - af(k) * ggii(kpm(k), j) &
                         - as(j) * ggii(k, j - 1) - an(j) * ggii(k, j + 1) &
                         - ac(k, j, iv) * ggii(k, j)
            sumres = max(sumres, abs(resd(k, j)))
            pc(k, j) = ggii(k, j)
          end do
        end do

        if (sumres < test) return
        if (id == 2) return

!  RESTRICT
        kc = -1
        do k = kkmg(nlev - 1, 1), kkmg(nlev - 1, 2)
          kc = kc + 2
          do j = 1, n2m
            resd(k, j) = resd(kc, j) * cor1(k) + resd(kc + 1, j) * cor2(k)
          end do
        end do

        return
      end subroutine toplevel
!=======================================================================
      subroutine trdiag1m(rr, uu, l1, l2, ll1, ll2, iv)
!=======================================================================
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: k, j, l1, l2, ll1, ll2, iv
        complex(8) :: rr(0:n3md, 0:n2), uu(0:n3md, 0:n2)

        do k = l1, l2, 2
          uu(k, ll1) = rr(k, ll1) * bet(k, 1, iv)
        end do

        do j = ll1 + 1, ll2
          do k = l1, l2, 2
            uu(k, j) = (rr(k, j) - as(j) * uu(k, j - 1)) * bet(k, j, iv)
          end do
        end do

        do j = ll2 - 1, ll1, -1
          do k = l1, l2, 2
            uu(k, j) = uu(k, j) - gam(k, j + 1, iv) * uu(k, j + 1)
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
        integer(8) :: k, j, ilev, kbgh
        complex(8) :: resd(n3md, n2m), ggii(0:n3md, 0:n2)

        kbgh = kkmg(ilev + 1, 1) - 1

        do k = kkmg(ilev, 1), kkmg(ilev, 2)
          kbgh = kbgh + 2
          do j = 1, n2m
            ggii(kbgh, j) = ggii(kbgh, j) + coi1(k) * ggii(k, j) + coi2(k) * ggii(kpm(k), j)
          end do
        end do

        return
      end subroutine goup

!=======================================================================
      subroutine godown(ilev, iv, resd, ggii)        ! COMPUTE RESIDUAL & RESTRICT IT
!=======================================================================
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: k, j, ilev, iv, kbg
        complex(8) :: resd(n3md, n2m), ggii(0:n3md, 0:n2)

        do j = 1, n2m
          do k = kkmg(ilev, 1), kkmg(ilev, 2)
            resd(k, j) = resd(k, j) - ab(k) * ggii(kmm(k), j) - af(k) * ggii(kpm(k), j) &
                         - as(j) * ggii(k, j - 1) - an(j) * ggii(k, j + 1) - ac(k, j, iv) * ggii(k, j)
          end do
        end do

        kbg = kkmg(ilev, 1) - 2

        do k = kkmg(ilev - 1, 1), kkmg(ilev - 1, 2)
          kbg = kbg + 2
          do j = 1, n2m
            resd(k, j) = resd(kbg, j) * cor1(k) + resd(kbg + 1, j) * cor2(k)
          end do
        end do

        return
      end subroutine godown

!=======================================================================
      subroutine relax(ilev, oldv, iiter, iv, resd, ggii)   ! ZEBRA VERSION
!=======================================================================
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: k, j, ilev, iiter, iv, kk
        complex(8) :: resd(n3md, n2m), ggii(0:n3md, 0:n2)
        real(8) :: oldv

        do j = 1, n2m
          do k = kkmg(ilev, 1), kkmg(ilev, 2), 2
            ggii(k, j) = resd(k, j) - oldv * (ab(k) * ggii(kmm(k), j) + af(k) * ggii(kpm(k), j))
          end do
        end do

        call trdiag1m(ggii, ggii, kkmg(ilev, 1), kkmg(ilev, 2), 1_8, n2m, iv)

        do j = 1, n2m
          do k = kkmg(ilev, 1) + 1, kkmg(ilev, 2), 2
            ggii(k, j) = resd(k, j) - ab(k) * ggii(kmm(k), j) - af(k) * ggii(kpm(k), j)
          end do
        end do

        call trdiag1m(ggii, ggii, kkmg(ilev, 1) + 1, kkmg(ilev, 2), 1_8, n2m, iv)

        do kk = 1, iiter - 1

          do j = 1, n2m
            do k = kkmg(ilev, 1), kkmg(ilev, 2), 2
              ggii(k, j) = resd(k, j) - (ab(k) * ggii(kmm(k), j) + af(k) * ggii(kpm(k), j))
            end do
          end do

          call trdiag1m(ggii, ggii, kkmg(ilev, 1), kkmg(ilev, 2), 1_8, n2m, iv)

          do j = 1, n2m
            do k = kkmg(ilev, 1) + 1, kkmg(ilev, 2), 2
              ggii(k, j) = resd(k, j) - ab(k) * ggii(kmm(k), j) - af(k) * ggii(kpm(k), j)
            end do
          end do

          call trdiag1m(ggii, ggii, kkmg(ilev, 1) + 1, kkmg(ilev, 2), 1_8, n2m, iv)

        end do
        return
      end subroutine relax
!=======================================================================
      subroutine gsor2d(u, rhs, iv, test, oldv)    ! 1 EQ. TYPE
!=======================================================================
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: k, j, iv
        real(8) :: test, oldv
        complex(8) :: u(0:n3, 0:n2), rhs(n3, n2)
        complex(8) :: ggii(0:n3md, 0:n2)
        complex(8) :: tt
        integer(8) :: ii
        real(8) :: ww, ww2, errmax

        ggii = 0.0_8
        ww = wwsor
        ww2 = 1.0_8 - ww
        ii = 0

!  HALF RELAX
!  ----------
        do j = 1, n2m
          do k = 1, n3m, 2
            ggii(k, j) = rhs(k, j) - oldv * (ab(k) * u(kmm(k), j) + af(k) * u(kpm(k), j))
          end do
        end do

        call trdiag1m(ggii, ggii, 1_8, n3m, 1_8, n2m, iv)

        do j = 1, n2m
          do k = 1, n3m, 2
            u(k, j) = ww * ggii(k, j) + oldv * ww2 * u(k, j)
          end do
        end do

!  ANOTHER HALF
!  ------------
        if (iv == 1) then
          tt = u(1, n2m)
          do j = 1, n2m
            do k = 1, n3m, 2
              u(k, j) = u(k, j) - tt
            end do
          end do
        end if

        do j = 1, n2m
          do k = 2, n3m, 2
            ggii(k, j) = rhs(k, j) - (ab(k) * u(kmm(k), j) + af(k) * u(kpm(k), j))
          end do
        end do

        call trdiag1m(ggii, ggii, 2_8, n3m, 1_8, n2m, iv)

        do j = 1, n2m
          do k = 2, n3m, 2
            u(k, j) = ww * ggii(k, j) + oldv * ww2 * u(k, j)
          end do
        end do

        call resid3(u, rhs, iv, errmax)
        if (errmax < test) goto 1000

!  MAIN ITERATION
!  ==============
        do ii = 1, mgitr

!  HALF RELAX
!  ----------
          do j = 1, n2m
            do k = 1, n3m, 2
              ggii(k, j) = rhs(k, j) - (ab(k) * u(kmm(k), j) + af(k) * u(kpm(k), j))
            end do
          end do

          call trdiag1m(ggii, ggii, 1_8, n3m, 1_8, n2m, iv)

          do j = 1, n2m
            do k = 1, n3m, 2
              u(k, j) = ww * ggii(k, j) + ww2 * u(k, j)
            end do
          end do

!  ANOTHER HALF
!  ------------

          if (iv == 1) then
            tt = u(1, n2m)
            do j = 1, n2m
              do k = 1, n3m, 2
                u(k, j) = u(k, j) - tt
              end do
            end do
          end if

          do j = 1, n2m
            do k = 2, n3m, 2
              ggii(k, j) = rhs(k, j) - (ab(k) * u(kmm(k), j) + af(k) * u(kpm(k), j))
            end do
          end do

          call trdiag1m(ggii, ggii, 2_8, n3m, 1_8, n2m, iv)

          do j = 1, n2m
            do k = 2, n3m, 2
              u(k, j) = ww * ggii(k, j) + ww2 * u(k, j)
            end do
          end do

          call resid3(u, rhs, iv, errmax)
          if (errmax < test) goto 1000

        end do

        print *, 'ITERATION LIMIT EXCEEDED.'
        write (77, 201) 'SOR', iv, ii, errmax * dtconst * float(n1mh) / float(n1m)
1000    continue

201     format(a5, '  IV=', i5, '  II=', i5, '  RG=', es23.15)

        return
      end subroutine gsor2d

!=======================================================================
      subroutine resid3(u, rhs, iv, errmax)
!=======================================================================
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: k, j, iv
        real(8) :: errmax
        complex(8) :: u(0:n3, 0:n2), rhs(n3, n2)
        complex(8) :: err

        errmax = 0.0_8

        do j = 1, n2m
          do k = 1, n3m
            err = rhs(k, j) - ab(k) * u(kmm(k), j) - af(k) * u(kpm(k), j) &
                  - as(j) * u(k, j - 1) - an(j) * u(k, j + 1) &
                  - ac(k, j, iv) * u(k, j)
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
        integer(8) :: minrow, kbg, kend, ksp, kc, kbgh
        real(8) :: kbz(n3md), kfz(n3md)
        real(8) :: zmpm(0:n3md)
        real(8) :: vdz_kbg, vdz_kend, sdz_kbg, sdz_kend

        minrow = n3m / (2**nlev)
        zmpm = 0.0_8

        levhalf = nint(nlev / 2.0_8)
        kkmg(nlev, 1) = 1    ! START INDEX
        kkmg(nlev, 2) = n3m  ! END INDEXT
        kkmg(nlev, 3) = n3m  ! THE NUMBER OF POINTS AT THE ILEV

        do k = 1, n3m
          zmpm(k) = zmp(k)
          kbz(k) = 1 - 1 / k                  ! 0 ONLY IF K=1
          kfz(k) = 1 - k / n3m                ! 0 ONLY IF K=N3M
        end do

!  COMPUTE FOR LOWER LEVELS
        do ilev = nlev - 1, 0, -1

          kkmg(ilev, 1) = kkmg(ilev + 1, 1) + kkmg(ilev + 1, 3)
          kkmg(ilev, 3) = minrow * (2**ilev)
          kkmg(ilev, 2) = kkmg(ilev, 1) + kkmg(ilev, 3) - 1

          kbg = kkmg(ilev, 1)
          kend = kkmg(ilev, 2)

          ksp = 2**(nlev - ilev)           ! WIDTH OF ONE CELL AT LOW LEVEL

          kc = 0
          do k = kbg, kend
            kc = kc + 1
            zmpm(k) = 0.5_8 * (z(kc * ksp + 1) + z((kc - 1) * ksp + 1))
            kbz(k) = 1 - kbg / k                  ! 0 ONLYIF K=KBG
            kfz(k) = 1 - k / kend                 ! 0 ONLY IF K=KEND
          end do

          kc = 0
          do k = kbg, kend
            kc = kc + 1
            ab(k) = kbz(k) / ((zmpm(k) - zmpm(k - 1)) * (z(kc * ksp + 1) - z((kc - 1) * ksp + 1)))
            af(k) = kfz(k) / ((zmpm(k + 1) - zmpm(k)) * (z(kc * ksp + 1) - z((kc - 1) * ksp + 1)))
          end do

          do j = 1, n2m
            do k = kkmg(ilev, 1), kkmg(ilev, 2)
              ac(k, j, 1) = -1.0_8 * (ab(k) + af(k) + as(j) + an(j))
            end do
          end do

          do i = 2, n1mh
            do j = 1, n2m
              do k = kkmg(ilev, 1), kkmg(ilev, 2)
                ac(k, j, i) = ac(k, j, 1) - ai3(i)
              end do
            end do
          end do

        end do

!  CALCULATE RESTRICTION COEFFS
        do ilev = nlev, 1, -1
          kbgh = kkmg(ilev, 1)
          do k = kkmg(ilev - 1, 1), kkmg(ilev - 1, 2)
            cor1(k) = (zmpm(kbgh + 1) - zmpm(k)) / (zmpm(kbgh + 1) - zmpm(kbgh))
            cor2(k) = 1.0_8 - cor1(k)
            kbgh = kbgh + 2
          end do
        end do

!  CALCULATE INTERPOLATION COEFFS
        do ilev = 0, nlev - 1
          kbgh = kkmg(ilev + 1, 1) + 1
          do k = kkmg(ilev, 1), kkmg(ilev, 2) - 1
            coi1(k) = (zmpm(k + 1) - zmpm(kbgh)) / (zmpm(k + 1) - zmpm(k))  ! * LOWER VALUE
            coi2(k) = 1.0_8 - coi1(k)
            kbgh = kbgh + 2
          end do
          k = kkmg(ilev, 2)
          coi1(k) = 1.0_8                 ! USE ONLY ONE LOWER POINT AT UPPER WALL
        end do

!===== FOR THE Z PERIODICIRY
!       INTRODUCE KPM & KMM
        if (zprdic == 1) then
          do ilev = nlev, 0, -1
            kbg = kkmg(ilev, 1)
            kend = kkmg(ilev, 2)
            do k = kbg, kend
              kpm(k) = k + 1
              kmm(k) = k - 1
            end do
            kpm(kend) = kbg
            kmm(kbg) = kend
          end do

          do ilev = nlev - 1, 0, -1
            kbg = kkmg(ilev, 1)
            kend = kkmg(ilev, 2)
            ksp = 2**(nlev - ilev)
            vdz_kbg = zmpm(kbg) - z(1) + z(n3) - zmpm(kend)
            vdz_kend = zmpm(kbg) - z(1) + z(n3) - zmpm(kend)
            sdz_kbg = z(1 + ksp) - z(1)
            sdz_kend = z(n3) - z(n3 - ksp)
            ab(kbg) = 1.0_8 / (vdz_kbg * sdz_kbg)
            af(kend) = 1.0_8 / (vdz_kend * sdz_kend)
          end do

          do ilev = nlev - 1, 0, -1
            do j = 1, n2m
              do k = kkmg(ilev, 1), kkmg(ilev, 2)
                ac(k, j, 1) = -1.0_8 * (ab(k) + af(k) + as(j) + an(j))
              end do
            end do

            do i = 2, n1mh
              do j = 1, n2m
                do k = kkmg(ilev, 1), kkmg(ilev, 2)
                  ac(k, j, i) = ac(k, j, 1) - ai3(i)
                end do
              end do
            end do
          end do

!  CALCULATE INTERPOLATION COEFFS
          do ilev = 0, nlev - 1
            kbg = kkmg(ilev, 1)
            kend = kkmg(ilev, 2)
            vdz_kend = zmpm(kbg) - z(1) + z(n3) - zmpm(kend)
            coi2(kend) = (zmpm(kkmg(ilev + 1, 2)) - zmpm(kend)) / vdz_kend
            coi1(kend) = 1.0_8 - coi2(kend)
          end do

        end if

        do ilev = nlev, 0, -1
          write (*, *) 'IIMG(1', ilev, ')=', kkmg(ilev, 1)
          write (*, *) 'IIMG(2', ilev, ')=', kkmg(ilev, 2)
          write (*, *) 'IIMG(3', ilev, ')=', kkmg(ilev, 3)
        end do

        return
      end subroutine coefmg
