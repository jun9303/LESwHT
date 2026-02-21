!=======================================================================
!
!     MAIN SOLVER OF POISSON EQUATION
!
!     AX=B
!     A: COEFFICIENT OF DISCRETIZED POISSON EQUATION
!     X: PHI: PSEUDO PRESSURE, OUTPUT OF SUBROUTINE POISSON
!     B: DIVGSUM: OUTPUT OF SUBROUTINE DIVGS.
!
!     X & Z DIRECTION: FOURIER TRANSFORM
!     Y-DIRECTION: TDMA
!
!     AK3,AK1: MATRIX COEFFICIENT (MODIFIED WAVENUMBER)
!     N3MH,N1MH: THE NUMBER OF WAVENUMBER INDEX
!
!     APR. 2010, J. LEE
!     JUN. 2017, J. PARK
!     FEB. 2026, S. LEE (SEVERAL FIXES)
!
!=======================================================================
      subroutine poisson(phi, divgsum)
!=======================================================================
        use mod_common
        use mod_poiss
        implicit none

        integer(8) :: i, j, k, jjp
        real(8) :: phi(0:n1, 0:n2, 0:n3), divgsum(0:n1, 0:n2, 0:n3)
        real(8) :: cn1, cn3
        real(8) :: ajcref, ajmref, ajpref
        complex(8) :: ccap(n3mh, n1, n2)
        complex(8) :: crhsref, phref

        real(8), dimension(:, :), allocatable :: ajc, ajm, ajp
        complex(8), dimension(:, :), allocatable :: zzz, crhs
        complex(8), dimension(:), allocatable :: zzzz, zzzz_b, xxxx, xxxx_b

        cn1 = 1.0_8 / float(n1m)
        cn3 = 1.0_8 / float(n3m)

! --- FORWARD FOURIER TRANSFORM
!$OMP PARALLEL PRIVATE(ZZZ, ZZZZ, XXXX, XXXX_B, ZZZZ_B)
        allocate (zzz(n3m, n1m))
        allocate (zzzz(n3m), zzzz_b(n3m * 2))
        allocate (xxxx(n1m), xxxx_b(n1m * 2))
        call zfft1d(xxxx, n1m, 0, xxxx_b)
        call zfft1d(zzzz, n3m, 0, zzzz_b)
!$OMP DO
        do j = 1, n2m
          do k = 1, n3m
            do i = 1, n1m
              zzz(k, i) = divgsum(i, j, k)
            end do
          end do

          do i = 1, n1m
            call zfft1d(zzz(1, i), n3m, -1, zzzz_b)
          end do

          do k = 1, n3mh
            do i = 1, n1m
              xxxx(i) = zzz(k, i)
            end do
            call zfft1d(xxxx, n1m, -1, xxxx_b)
            do i = 1, n1m
              ccap(k, i, j) = xxxx(i)
            end do
          end do
        end do
!$OMP END DO
        deallocate (zzz, zzzz, zzzz_b, xxxx, xxxx_b)
!$OMP END PARALLEL

! --- SOLVE TDMA MATRIX
!$OMP PARALLEL PRIVATE(AJM, AJP, AJC, CRHS, CRHSREF, AJCREF, AJMREF, AJPREF, PHREF)
        allocate (crhs(n2, n1))
        allocate (ajm(n2, n1), ajp(n2, n1), ajc(n2, n1))
!$OMP DO
        do k = 1, n3mh
          do i = 1, n1m
            do j = 1, n2m
              jjp = jpv(j)
              ajm(j, i) = f2fyi(j) * c2cyi(j) * (1.0_8 - fixjl(j))
              ajp(j, i) = f2fyi(j) * c2cyi(jjp) * (1.0_8 - fixju(j))
              ajc(j, i) = -((ajm(j, i) + ajp(j, i) + ak1(i) + ak3(k)) * &
                            (1.0_8 - fixjl(j)) * (1.0_8 - fixju(j)) + &
                            (f2fyi(1) * c2cyi(2) + ak1(i) + ak3(k)) * fixjl(j) + &
                            (f2fyi(n2m) * c2cyi(n2m) + ak1(i) + ak3(k)) * fixju(j))
              crhs(j, i) = ccap(k, i, j)
            end do
          end do

          if (k == 1) then
            crhsref = crhs(n2m, 1)
            ajcref = ajc(n2m, 1)
            ajmref = ajm(n2m, 1)
            ajpref = ajp(n2m, 1)
            crhs(n2m, 1) = 0.0_8
            ajc(n2m, 1) = 1.0_8
            ajm(n2m, 1) = 0.0_8
            ajp(n2m, 1) = 0.0_8
          end if

          call ctrdiag(ajm, ajc, ajp, crhs, 1_8, n2m, crhs, n1m)

          if (k == 1) then
            phref = (-ajmref * crhs(n2m - 1, 1) + crhsref) / ajcref
            do j = 1, n2m
              crhs(j, 1) = crhs(j, 1) - phref
            end do
            crhs(n2m, 1) = 0.0_8
          end if

          do i = 1, n1m
            do j = 1, n2m
              ccap(k, i, j) = crhs(j, i)
            end do
          end do
        end do
!$OMP END DO
        deallocate (crhs, ajm, ajp, ajc)
!$OMP END PARALLEL

! --- INVERSE FOURIER TRANSFORM
!$OMP PARALLEL PRIVATE(ZZZ, ZZZZ, XXXX, XXXX_B, ZZZZ_B)
        allocate (zzz(n3m, n1m))
        allocate (zzzz(n3m), zzzz_b(n3m * 2))
        allocate (xxxx(n1m), xxxx_b(n1m * 2))
        call zfft1d(xxxx, n1m, 0, xxxx_b)
        call zfft1d(zzzz, n3m, 0, zzzz_b)
!$OMP DO
        do j = 1, n2m
          do k = 1, n3mh
            do i = 1, n1m
              xxxx(i) = ccap(k, i, j)
            end do
            call zfft1d(xxxx, n1m, 1, xxxx_b)
            do i = 1, n1m
              zzz(k, i) = xxxx(i)
            end do
          end do

          do i = 1, n1m
            do k = n3mh + 1, n3m
              zzz(k, i) = conjg(zzz(n3m + 2 - k, i))
            end do
          end do

          do i = 1, n1m
            call zfft1d(zzz(1, i), n3m, 1, zzzz_b)
            do k = 1, n3m
              phi(i, j, k) = real(zzz(k, i), 8)
            end do
          end do
        end do
!$OMP END DO
        deallocate (zzz, zzzz, zzzz_b, xxxx, xxxx_b)
!$OMP END PARALLEL

        return
      end subroutine poisson

!=======================================================================
      subroutine poisinit
!=======================================================================
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, k, kk
        real(8) :: pi
        real(8) :: sdzis, sdxis

        call ftft_allo

! --- DEFINE MODIFIED WAVENUMBERS
        pi = 2.0_8 * asin(1.0_8)

        do k = 1, n3mh
          ai3(k) = float(k - 1) * 2.0_8 * pi
        end do
        ai3(1) = 0.0_8

        do i = 1, n1m
          ai1(i) = float(i - 1) * 2.0_8 * pi
        end do
        ai1(1) = 0.0_8

        sdzis = f2fzi(1)
        sdxis = f2fxi(1)

        do kk = 1, n3mh
          ak3(kk) = 2.0_8 * (1.0_8 - cos(ai3(kk) / float(n3m))) * sdzis * sdzis
        end do

        do kk = 1, n1mh
          ak1(kk) = 2.0_8 * (1.0_8 - cos(ai1(kk) / float(n1m))) * sdxis * sdxis
        end do

        do kk = n1m, n1mh + 1, -1
          ak1(kk) = ak1(n1m + 2 - kk)
        end do

        return
      end subroutine poisinit

!=======================================================================
      subroutine ctrdiag(a, b, c, r, ni, nf, uu, mf)
!=======================================================================
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j
        integer(8) :: ni, nf, mf
        real(8) :: a(n2, n1), b(n2, n1), c(n2, n1)
        complex(8) :: r(n2, n1), uu(n2, n1)
        real(8), dimension(:, :), allocatable :: gam_l
        real(8), dimension(:), allocatable :: bet_l

        allocate (gam_l(n2, n1))
        allocate (bet_l(n1))

        do i = 1, mf
          bet_l(i) = 1.0_8 / b(ni, i)
          uu(ni, i) = r(ni, i) * bet_l(i)
        end do

        do i = 1, mf
          do j = ni + 1, nf
            gam_l(j, i) = c(j - 1, i) * bet_l(i)
            bet_l(i) = 1.0_8 / (b(j, i) - a(j, i) * gam_l(j, i))
            uu(j, i) = (r(j, i) - a(j, i) * uu(j - 1, i)) * bet_l(i)
          end do
        end do

        do i = 1, mf
          do j = nf - 1, ni, -1
            uu(j, i) = uu(j, i) - gam_l(j + 1, i) * uu(j + 1, i)
          end do
        end do

        deallocate (gam_l, bet_l)

        return
      end subroutine ctrdiag
