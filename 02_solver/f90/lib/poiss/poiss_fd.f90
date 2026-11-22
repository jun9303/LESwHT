!=======================================================================
!
!     MAIN SOLVER OF POISSON EQUATION
!
!     AX=B
!     A: COEFFICIENT OF DISCRETIZED POISSON EQUATION
!     X: PHI: PSEUDO PRESSURE, OUTPUT OF SUBROUTINE POISSON
!     B: DIVGSUM: OUTPUT OF SUBROUTINE DIVGS.
!
!     DIRECT SOLVER (FAST DIAGONALIZATION), NO ITERATION
!     X-DIRECTION: EIGENVECTORS OF THE 1D OPERATOR (ANY SPACING; FOURIER IF PERIODIC)
!     Z-DIRECTION: EIGENVECTORS OF THE 1D OPERATOR (ANY SPACING; FOURIER IF PERIODIC)
!     Y-DIRECTION: TDMA, ONE FOR EVERY (X-MODE, Z-MODE) PAIR
!
!     THE 1D OPERATORS ARE THOSE OF POISS_Z_FT / POISS_X_FT (AW, AE / AB, AF / AS, AN).
!     MULTIPLIED BY THE CELL WIDTH THEY ARE SYMMETRIC, SO WITH SX = SQRT(F2FX) THE MATRIX
!     T = SX * L * SX^-1 IS SYMMETRIC AND L = SX^-1 * Q * LAMBDA * Q^T * SX (Q ORTHOGONAL).
!     A NON-PERIODIC DIRECTION IS TRIDIAGONAL (IMPLICIT QL); A PERIODIC ONE MUST BE UNIFORM
!     AND USES THE REAL FOURIER BASIS. EVERY LINE OF EVERY STEP IS INDEPENDENT.
!     THE MEAN MODE OF THE ALL-NEUMANN / PERIODIC PROBLEM IS SINGULAR: ITS LAST ROW IS
!     REPLACED BY PHI(N2M) = 0 (AS IN POISS_FTFT). THE CONSTANT IN PHI DOES NOT MATTER
!     (UCALC USES ITS GRADIENT, PCALC REMOVES ITS OUTLET MEAN).
!
!     NOV. 2026
!
!=======================================================================
      subroutine poisinit
!=======================================================================
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, jp

        if (yprdic .eq. 1 .or. n2m .lt. 2) then
          write (*, *) ' [poiss_fd] y must be non-periodic with at least two cells; use another POISS'
          stop
        end if

        call fd_allo

        do j = 1, n2m
          jp = jpv(j)
          as(j) = (1.0_8 - fixjl(j)) * c2cyi(j) * f2fyi(j)
          an(j) = (1.0_8 - fixju(j)) * c2cyi(jp) * f2fyi(j)
        end do

        call fd_modes(n1m, xprdic, f2fx, c2cxi, fixil, fixiu, ipv, qfx, lfx, kfx0)
        call fd_modes(n3m, zprdic, f2fz, c2czi, fixkl, fixku, kpv, qfz, lfz, kfz0)

        do i = 1, n1m
          sfx(i) = sqrt(f2fx(i))
        end do
        do i = 1, n3m
          sfz(i) = sqrt(f2fz(i))
        end do
        qfxt = transpose(qfx)
        qfzt = transpose(qfz)

        write (*, '(A,2I6,A,2ES12.4)') ' [poiss_fd] direct solver, x and z modes:', n1m, n3m, &
          '   largest |lambda| x, z:', maxval(abs(lfx)), maxval(abs(lfz))

        return
      end subroutine poisinit

!=======================================================================
      subroutine poisson(phi, divgsum)
!=======================================================================
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, k, kx, kz
        integer(8) :: adi_nblk, nbj, nbx, ib, j1, j2, k1, k2
        real(8) :: phi(0:n1, 0:n2, 0:n3), divgsum(0:n1, 0:n2, 0:n3)
        real(8), dimension(:), allocatable :: bt
        real(8), dimension(:, :), allocatable :: gm
        real(8) :: crref, phref

        nbj = adi_nblk(n3m, n2m)
        nbx = adi_nblk(n3m, n1m)

! --- SCALE AND TRANSFORM IN X: B(KX,J,K) = SUM_I QFX(I,KX) * SFX(I) * SFZ(K) * DIVGSUM(I,J,K)
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(J1, J2, I, J)
        do k = 1, n3m
          do ib = 1, nbj
            j1 = 1 + ((ib - 1) * n2m) / nbj
            j2 = (ib * n2m) / nbj
            do j = j1, j2
              do i = 1, n1m
                fda(i, j, k) = sfx(i) * sfz(k) * divgsum(i, j, k)
              end do
            end do
            fdb(:, j1:j2, k) = matmul(qfxt, fda(:, j1:j2, k))
          end do
        end do
!$OMP END PARALLEL DO

! --- TRANSFORM IN Z: C(KX,J,KZ) = SUM_K B(KX,J,K) * QFZ(K,KZ)
        if (n3m .gt. 1) then
!$OMP PARALLEL DO
          do j = 1, n2m
            fdc(:, j, :) = matmul(fdb(:, j, :), qfz)
          end do
!$OMP END PARALLEL DO
        else
!$OMP PARALLEL DO
          do j = 1, n2m
            fdc(:, j, 1) = fdb(:, j, 1) * qfz(1, 1)
          end do
!$OMP END PARALLEL DO
        end if

! --- SOLVE TDMA IN Y FOR EVERY (KX, KZ); LINES KX ARE VECTORIZED IN BLOCKS
!$OMP PARALLEL PRIVATE(BT, GM, K1, K2, KX, J, CRREF, PHREF)
        allocate (bt(n1m), gm(n1m, n2m))
!$OMP DO COLLAPSE(2)
        do kz = 1, n3m
          do ib = 1, nbx
            k1 = 1 + ((ib - 1) * n1m) / nbx
            k2 = (ib * n1m) / nbx
            if (kz .eq. kfz0 .and. kfx0 .ge. k1 .and. kfx0 .le. k2) crref = fdc(kfx0, n2m, kz)

            do kx = k1, k2
              bt(kx) = 1.0_8 / (-(as(1) + an(1)) + lfx(kx) + lfz(kz))
              fdc(kx, 1, kz) = fdc(kx, 1, kz) * bt(kx)
            end do

            do j = 2, n2m - 1
              do kx = k1, k2
                gm(kx, j) = an(j - 1) * bt(kx)
                bt(kx) = 1.0_8 / (-(as(j) + an(j)) + lfx(kx) + lfz(kz) - as(j) * gm(kx, j))
                fdc(kx, j, kz) = (fdc(kx, j, kz) - as(j) * fdc(kx, j - 1, kz)) * bt(kx)
              end do
            end do

            j = n2m
            do kx = k1, k2
              gm(kx, j) = an(j - 1) * bt(kx)
              if (kx .eq. kfx0 .and. kz .eq. kfz0) then
                fdc(kx, j, kz) = 0.0_8                    ! SINGULAR MODE: PHI(N2M) = 0
              else
                bt(kx) = 1.0_8 / (-(as(j) + an(j)) + lfx(kx) + lfz(kz) - as(j) * gm(kx, j))
                fdc(kx, j, kz) = (fdc(kx, j, kz) - as(j) * fdc(kx, j - 1, kz)) * bt(kx)
              end if
            end do

            do j = n2m - 1, 1, -1
              do kx = k1, k2
                fdc(kx, j, kz) = fdc(kx, j, kz) - gm(kx, j + 1) * fdc(kx, j + 1, kz)
              end do
            end do

!           SINGULAR MODE: SHIFT AS IN POISS_FTFT (THE LAST ROW OF THE ORIGINAL SYSTEM HOLDS)
            if (kz .eq. kfz0 .and. kfx0 .ge. k1 .and. kfx0 .le. k2) then
              phref = (-as(n2m) * fdc(kfx0, n2m - 1, kz) + crref) / (-(as(n2m) + an(n2m)))
              do j = 1, n2m
                fdc(kfx0, j, kz) = fdc(kfx0, j, kz) - phref
              end do
              fdc(kfx0, n2m, kz) = 0.0_8
            end if
          end do
        end do
!$OMP END DO
        deallocate (bt, gm)
!$OMP END PARALLEL

! --- INVERSE TRANSFORM IN Z
        if (n3m .gt. 1) then
!$OMP PARALLEL DO
          do j = 1, n2m
            fdb(:, j, :) = matmul(fdc(:, j, :), qfzt)
          end do
!$OMP END PARALLEL DO
        else
!$OMP PARALLEL DO
          do j = 1, n2m
            fdb(:, j, 1) = fdc(:, j, 1) * qfz(1, 1)
          end do
!$OMP END PARALLEL DO
        end if

! --- INVERSE TRANSFORM IN X AND UNSCALE
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(J1, J2, I, J)
        do k = 1, n3m
          do ib = 1, nbj
            j1 = 1 + ((ib - 1) * n2m) / nbj
            j2 = (ib * n2m) / nbj
            fda(:, j1:j2, k) = matmul(qfx, fdb(:, j1:j2, k))
            do j = j1, j2
              do i = 1, n1m
                phi(i, j, k) = fda(i, j, k) / (sfx(i) * sfz(k))
              end do
            end do
          end do
        end do
!$OMP END PARALLEL DO

        return
      end subroutine poisson

!=======================================================================
      subroutine fd_modes(n, iprd, f2f, c2ci, fixl, fixu, ipp, q, lam, k0)
!=======================================================================
!     EIGENVECTORS Q (COLUMNS) AND EIGENVALUES LAM OF ONE 1D OPERATOR
!       (L PHI)(I) = WEST(I) * (PHI(I-1) - PHI(I)) + EAST(I) * (PHI(I+1) - PHI(I))
!       WEST(I) = (1 - FIXL(I)) * C2CI(I) / F2F(I),  EAST(I) = (1 - FIXU(I)) * C2CI(IPP(I)) / F2F(I)
!     IN THE SYMMETRIC FORM T = S * L * S^-1, S = SQRT(F2F).  K0: INDEX OF THE ZERO MODE.
!-----------------------------------------------------------------------
        implicit none
        integer(8) :: n, iprd, k0
        real(8) :: f2f(0:*), c2ci(0:*), q(n, n), lam(n)
        integer(8) :: fixl(*), fixu(*), ipp(*)
        integer(8) :: i, m
        real(8) :: pi, h, west, east, dmin
        real(8), dimension(:), allocatable :: e

        pi = acos(-1.0_8)
        q = 0.0_8
        lam = 0.0_8

        if (n .eq. 1) then
          q(1, 1) = 1.0_8
          lam(1) = 0.0_8
          k0 = 1
          return
        end if

        if (iprd .eq. 1) then
!         PERIODIC: THE GRID MUST BE UNIFORM; REAL FOURIER BASIS
          h = f2f(1)
          do i = 2, n
            if (abs(f2f(i) - h) .gt. 1.0d-10 * h) then
              write (*, *) ' [poiss_fd] a periodic direction must have uniform spacing; use another POISS'
              stop
            end if
          end do
          do i = 1, n
            q(i, 1) = 1.0_8 / sqrt(dble(n))
          end do
          lam(1) = 0.0_8
          do m = 1, (n - 1) / 2
            do i = 1, n
              q(i, 2 * m) = sqrt(2.0_8 / dble(n)) * cos(2.0_8 * pi * dble(m) * dble(i - 1) / dble(n))
              q(i, 2 * m + 1) = sqrt(2.0_8 / dble(n)) * sin(2.0_8 * pi * dble(m) * dble(i - 1) / dble(n))
            end do
            lam(2 * m) = -4.0_8 / (h * h) * sin(pi * dble(m) / dble(n))**2
            lam(2 * m + 1) = lam(2 * m)
          end do
          if (mod(n, 2_8) .eq. 0) then
            do i = 1, n
              q(i, n) = dble(1 - 2 * mod(i - 1, 2_8)) / sqrt(dble(n))
            end do
            lam(n) = -4.0_8 / (h * h)
          end if
          k0 = 1
          return
        end if

!       NON-PERIODIC: SYMMETRIC TRIDIAGONAL T, IMPLICIT QL
        allocate (e(n))
        e = 0.0_8
        do i = 1, n
          west = (1.0_8 - fixl(i)) * c2ci(i) / f2f(i)
          east = (1.0_8 - fixu(i)) * c2ci(ipp(i)) / f2f(i)
          lam(i) = -(west + east)
          if (i .gt. 1) e(i) = (1.0_8 - fixl(i)) * c2ci(i) / sqrt(f2f(i) * f2f(i - 1))
          q(i, i) = 1.0_8
        end do
        call fd_tql(n, lam, e, q)
        deallocate (e)

!       THE NEUMANN OPERATOR HAS ONE ZERO EIGENVALUE: MAKE IT EXACT
        k0 = 1
        dmin = abs(lam(1))
        do i = 2, n
          if (abs(lam(i)) .lt. dmin) then
            dmin = abs(lam(i))
            k0 = i
          end if
        end do
        lam(k0) = 0.0_8

        return
      end subroutine fd_modes

!=======================================================================
      subroutine fd_tql(n, d, e, z)
!=======================================================================
!     IMPLICIT QL WITH SHIFTS FOR A SYMMETRIC TRIDIAGONAL MATRIX
!     D: DIAGONAL (OUT: EIGENVALUES), E(2:N): SUBDIAGONAL, Z: IDENTITY IN, EIGENVECTORS OUT
!-----------------------------------------------------------------------
        implicit none
        integer(8) :: n
        real(8) :: d(n), e(n), z(n, n)
        integer(8) :: i, k, l, m, iter
        real(8) :: b, c, dd, f, g, p, r, s
        logical :: under

        do i = 2, n
          e(i - 1) = e(i)
        end do
        e(n) = 0.0_8

        do l = 1, n
          iter = 0
          do
            do m = l, n - 1
              dd = abs(d(m)) + abs(d(m + 1))
              if (abs(e(m)) .le. epsilon(1.0_8) * dd) exit
            end do
            if (m .eq. l) exit
            iter = iter + 1
            if (iter .gt. 60) then
              write (*, *) ' [poiss_fd] QL did not converge'
              stop
            end if
            g = (d(l + 1) - d(l)) / (2.0_8 * e(l))
            r = hypot(g, 1.0_8)
            g = d(m) - d(l) + e(l) / (g + sign(r, g))
            s = 1.0_8
            c = 1.0_8
            p = 0.0_8
            under = .false.
            do i = m - 1, l, -1
              f = s * e(i)
              b = c * e(i)
              r = hypot(f, g)
              e(i + 1) = r
              if (r .eq. 0.0_8) then
                d(i + 1) = d(i + 1) - p
                e(m) = 0.0_8
                under = .true.
                exit
              end if
              s = f / r
              c = g / r
              g = d(i + 1) - p
              r = (d(i) - g) * s + 2.0_8 * c * b
              p = s * r
              d(i + 1) = g + p
              g = c * r - b
              do k = 1, n
                f = z(k, i + 1)
                z(k, i + 1) = s * z(k, i) + c * f
                z(k, i) = c * z(k, i) - s * f
              end do
            end do
            if (under) cycle
            d(l) = d(l) - p
            e(l) = g
            e(m) = 0.0_8
          end do
        end do

        return
      end subroutine fd_tql
