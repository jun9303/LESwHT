!=======================================================================
!
!     CODEBASE (LICA 2017 VERSION) BY
!     H. CHOI / DEPARTMENT OF MECHANICAL & AEROSPACE ENGINEERING
!     SEOUL NATIONAL UNIVERSITY
!
!=======================================================================
!
!     LESWHT (C) 2026 S. LEE (ORCID: 0000-0002-2063-6298)
!     2018.02.28. MODIFIED FOR F2PY (FORTRAN-TO-PYTHON) USAGE
!     2026.02.18. CODE MODERNIZATION
!
!=======================================================================
subroutine find_inout(nx, ny, nz, xcoord, ycoord, zcoord, nbody, inout, t)
!$ USE OMP_LIB
  implicit none

  integer(8), intent(in) :: nx, ny, nz
  real(8), intent(in) :: xcoord(0:nx), ycoord(0:ny), zcoord(0:nz), t
  integer(8), intent(out) :: inout(0:nx, 0:ny, 0:nz)
  integer(8), intent(out) :: nbody

  integer(8) :: i, j, k
  real(8) :: val
  ! EXTERNAL FUNCTION DECLARATION
  real(8), external :: funcbody

  nbody = 0
  inout = 1

  !$OMP PARALLEL DO REDUCTION(+:NBODY) PRIVATE(I, J, K, VAL)
  do k = 0, nz
    do j = 0, ny
      do i = 0, nx
        val = funcbody(xcoord(i), ycoord(j), zcoord(k), t)
        if (val .le. 1.0d-10) then
          nbody = nbody + 1
          inout(i, j, k) = 0
        end if
      end do
    end do
  end do
  !$OMP END PARALLEL DO

  return
end subroutine find_inout

!=======================================================================
subroutine findbdy_intp(dir, nx, ny, nz, nbody, inout, ufix_x, &
                        ufix_y, ufix_z, lfix_x, lfix_y, lfix_z, &
                        nintp, ninner, fcp, intptype, intpindx)
!$ USE OMP_LIB
  implicit none

  integer(8), intent(in) :: dir, nx, ny, nz, nbody
  integer(8), intent(in) :: inout(0:nx, 0:ny, 0:nz)
  integer(8), intent(in) :: ufix_x(nx - 1), ufix_y(ny - 1), ufix_z(nz - 1)
  integer(8), intent(in) :: lfix_x(nx - 1), lfix_y(ny - 1), lfix_z(nz - 1)

  integer(8), intent(out) :: nintp, ninner
  integer(8), intent(out) :: fcp(nbody, 3)
  integer(8), intent(out) :: intptype(nbody, 1), intpindx(nbody, 3)

  integer(8) :: i, j, k, ip, im, jp, jm, kp, km
  integer(8) :: istart, jstart, kstart
  integer(8) :: inout_sum

  ! TEMPORARY ARRAY FOR INNER POINTS TO AVOID RACE CONDITIONS OR COMPLEX SORTING
  integer(8), allocatable :: fcp_temp(:, :)

  allocate (fcp_temp(nbody, 3))

  istart = 1; jstart = 1; kstart = 1
  if (dir .eq. 1) istart = 2
  if (dir .eq. 2) jstart = 2
  if (dir .eq. 3) kstart = 2

  nintp = 0
  ninner = 0

  ! SERIAL LOOP TO CLASSIFY POINTS (COULD BE PARALLELIZED WITH OFFSET CALCULATION)
  do k = kstart, nz - 1
    km = k - 1; kp = k + 1
    do j = jstart, ny - 1
      jm = j - 1; jp = j + 1
      do i = istart, nx - 1
        im = i - 1; ip = i + 1

        inout_sum = (abs(inout(ip, j, k) - inout(im, j, k))) * (1 - ufix_x(i)) * (1 - lfix_x(i)) &
                    + (abs(inout(i, jp, k) - inout(i, jm, k))) * (1 - ufix_y(j)) * (1 - lfix_y(j)) &
                    + (abs(inout(i, j, kp) - inout(i, j, km))) * (1 - ufix_z(k)) * (1 - lfix_z(k))

        if ((inout(i, j, k) .eq. 0) .and. (inout_sum .gt. 0)) then
          nintp = nintp + 1
          fcp(nintp, 1) = i
          fcp(nintp, 2) = j
          fcp(nintp, 3) = k
          intptype(nintp, 1) = inout_sum
          intpindx(nintp, 1) = inout(ip, j, k) - inout(im, j, k)
          intpindx(nintp, 2) = inout(i, jp, k) - inout(i, jm, k)
          intpindx(nintp, 3) = inout(i, j, kp) - inout(i, j, km)
        else if ((inout(i, j, k) .eq. 0) .and. (inout_sum .eq. 0)) then
          ninner = ninner + 1
          fcp_temp(ninner, 1) = i
          fcp_temp(ninner, 2) = j
          fcp_temp(ninner, 3) = k
        end if
      end do
    end do
  end do

  ! COPY INNER POINTS TO THE END OF FCP
  !$OMP PARALLEL DO PRIVATE(I)
  do i = 1, ninner
    fcp(nintp + i, 1) = fcp_temp(i, 1)
    fcp(nintp + i, 2) = fcp_temp(i, 2)
    fcp(nintp + i, 3) = fcp_temp(i, 3)
  end do
  !$OMP END PARALLEL DO

  deallocate (fcp_temp)

end subroutine findbdy_intp

!=======================================================================
subroutine geomfac_preset(ni, icoord, im, prdic, i_adj)
!$ USE OMP_LIB
  implicit none

  integer(8), intent(in) :: ni
  real(8), intent(in) :: icoord(0:ni), im(0:ni)
  character(*), intent(in) :: prdic
  real(8), intent(out) :: i_adj(-1:ni + 1, 3)
  integer(8) :: l

  i_adj = 0.0d0

  do l = 0, ni
    i_adj(l, 1) = icoord(l)
    i_adj(l, 2) = im(l)
    i_adj(l, 3) = im(l)
  end do

  if (prdic .eq. 'on') then
    i_adj(-1, 1)   = icoord(0) - (icoord(ni) - icoord(ni-1))
    i_adj(ni+1, 1) = icoord(ni) + (icoord(1) - icoord(0))

    i_adj(0, 2)    = icoord(0) - 0.5d0*(icoord(ni) - icoord(ni-1))
    i_adj(-1, 2)   = i_adj(-1, 1) - 0.5d0*(icoord(ni-1) - icoord(ni-2))
        
    i_adj(ni+1, 2) = icoord(ni) + 0.5d0*(icoord(1) - icoord(0))

    i_adj(0, 3)    = i_adj(0, 2)
    i_adj(-1, 3)   = i_adj(-1, 2)
    i_adj(ni+1, 3) = i_adj(ni+1, 2)
  end if
end subroutine geomfac_preset

!=======================================================================
subroutine geomfac_intp(nx, ny, nz, xpre, ypre, zpre, nintp, &
                        nbody, fcp, intpindx, geomfac, t)
!$ USE OMP_LIB
  implicit none

  integer(8), intent(in) :: nx, ny, nz, nintp, nbody
  real(8), intent(in) :: xpre(-1:nx + 1), ypre(-1:ny + 1), zpre(-1:nz + 1)
  integer(8), intent(in) :: fcp(nbody, 3), intpindx(nintp, 3)
  real(8), intent(out) :: geomfac(nintp, 3, 3, 3)
  real(8), intent(in) :: t

  integer(8) :: l, m, n, i, j, k, indic
  real(8) :: x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3
  real(8) :: xtemp1, ytemp1, ztemp1, xtemp2, ytemp2, ztemp2
  real(8) :: ffs, ffe, ff1, ff2
  real(8) :: xx1, yy1, zz1, xx2, yy2, zz2
  real(8) :: ddx, ddy, ddz, dx1, dx2, dx3, dy1, dy2, dy3, dz1, dz2, dz3
  real(8) :: a0, b0, c0, a1, b1, c1

  real(8), external :: funcbody

  !$OMP PARALLEL DO PRIVATE(L, M, N, I, J, K, INDIC, X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3) &
  !$OMP PRIVATE(XTEMP1, YTEMP1, ZTEMP1, XTEMP2, YTEMP2, ZTEMP2, FFS, FFE, FF1, FF2) &
  !$OMP PRIVATE(XX1, YY1, ZZ1, XX2, YY2, ZZ2, DDX, DDY, DDZ) &
  !$OMP PRIVATE(DX1, DX2, DX3, DY1, DY2, DY3, DZ1, DZ2, DZ3, A0, B0, C0, A1, B1, C1)
  do l = 1, nintp
    x1 = xpre(fcp(l, 1))
    y1 = ypre(fcp(l, 2))
    z1 = zpre(fcp(l, 3))

    x2 = xpre(fcp(l, 1) + intpindx(l, 1))
    y2 = ypre(fcp(l, 2) + intpindx(l, 2))
    z2 = zpre(fcp(l, 3) + intpindx(l, 3))

    indic = 0
    x0 = x1; y0 = y1; z0 = z1 ! INITIALIZE TO AVOID WARNING

    do m = 1, 3
      if (m .eq. 1) then
        xtemp1 = x1; ytemp1 = y1; ztemp1 = z1
        xtemp2 = x2; ytemp2 = y2; ztemp2 = z2
        ffs = funcbody(x1, y1, z1, t)
        ffe = funcbody(x2, y2, z2, t)

        if (ffs * ffe .gt. 0.0d0) then
          ! SHOULD NOT HAPPEN IF SURFACE IS BETWEEN 1 AND 2, BUT FAIL-SAFE
          geomfac(l, :, :, :) = 0.0d0
          geomfac(l, 1, 1, 1) = 1.0d0
          goto 45
        end if
      else
        xtemp1 = xx1; ytemp1 = yy1; ztemp1 = zz1
        xtemp2 = xx2; ytemp2 = yy2; ztemp2 = zz2
      end if

      ! BISECT/LINEAR SEARCH FOR SURFACE
      do n = 0, 19
        ddx = xtemp2 - xtemp1
        ddy = ytemp2 - ytemp1
        ddz = ztemp2 - ztemp1

        xx1 = xtemp1 + ddx * dble(n) / 20.0d0
        xx2 = xtemp1 + ddx * dble(n + 1) / 20.0d0
        yy1 = ytemp1 + ddy * dble(n) / 20.0d0
        yy2 = ytemp1 + ddy * dble(n + 1) / 20.0d0
        zz1 = ztemp1 + ddz * dble(n) / 20.0d0
        zz2 = ztemp1 + ddz * dble(n + 1) / 20.0d0

        ff1 = funcbody(xx1, yy1, zz1, t)
        ff2 = funcbody(xx2, yy2, zz2, t)

        if (ff1 .eq. 0.0d0) then
          x0 = xx1; y0 = yy1; z0 = zz1
          goto 33
        else if (ff2 .eq. 0.0d0) then
          x0 = xx2; y0 = yy2; z0 = zz2
          goto 33
        else if (ff1 * ff2 .lt. 0.0d0) then
          x0 = 0.5d0 * (xx1 + xx2)
          y0 = 0.5d0 * (yy1 + yy2)
          z0 = 0.5d0 * (zz1 + zz2)
          if (m .eq. 3) goto 33
          goto 22
        end if
      end do
22    continue
    end do
33  continue

    x3 = xpre(fcp(l, 1) + intpindx(l, 1) * 2)
    y3 = ypre(fcp(l, 2) + intpindx(l, 2) * 2)
    z3 = zpre(fcp(l, 3) + intpindx(l, 3) * 2)

    dx1 = abs(x1 - x0); dx2 = abs(x2 - x0); dx3 = abs(x3 - x0)
    dy1 = abs(y1 - y0); dy2 = abs(y2 - y0); dy3 = abs(y3 - y0)
    dz1 = abs(z1 - z0); dz2 = abs(z2 - z0); dz3 = abs(z3 - z0)

    ! WEIGHT CALCULATION (UNCHANGED LOGIC)
    if (intpindx(l, 1) .eq. 0) then
      a0 = 1.0d0; a1 = 1.0d0
    else if (dx2 .ge. dx1) then
      a0 = dx2 / (dx1 + dx2)
      a1 = 1.0d0
    else
      a0 = 0.5d0
      a1 = (dx3 - dx1) / (dx3 - dx2)
    end if

    if (intpindx(l, 2) .eq. 0) then
      b0 = 1.0d0; b1 = 1.0d0
    else if (dy2 .ge. dy1) then
      b0 = dy2 / (dy1 + dy2)
      b1 = 1.0d0
    else
      b0 = 0.5d0
      b1 = (dy3 - dy1) / (dy3 - dy2)
    end if

    if (intpindx(l, 3) .eq. 0) then
      c0 = 1.0d0; c1 = 1.0d0
    else if (dz2 .ge. dz1) then
      c0 = dz2 / (dz1 + dz2)
      c1 = 1.0d0
    else
      c0 = 0.5d0
      c1 = (dz3 - dz1) / (dz3 - dz2)
    end if

    ! POPULATE GEOMFAC (DIRECT ASSIGNMENTS)
    geomfac(l, 1, 1, 1) = 1./(a0 * b0 * c0)
    geomfac(l, 1, 1, 2) = -1./(a0 * b0 * c0) * (a0) * (b0) * (1.-c0) * (c1)
    geomfac(l, 1, 1, 3) = -1./(a0 * b0 * c0) * (a0) * (b0) * (1.-c0) * (1.-c1)
    geomfac(l, 1, 2, 1) = -1./(a0 * b0 * c0) * (a0) * (1.-b0) * (c0) * (b1)
    geomfac(l, 1, 2, 2) = -1./(a0 * b0 * c0) * (a0) * (1.-b0) * (1.-c0) * (b1) * (c1)
    geomfac(l, 1, 2, 3) = -1./(a0 * b0 * c0) * (a0) * (1.-b0) * (1.-c0) * (b1) * (1.-c1)
    geomfac(l, 1, 3, 1) = -1./(a0 * b0 * c0) * (a0) * (1.-b0) * (c0) * (1.-b1)
    geomfac(l, 1, 3, 2) = -1./(a0 * b0 * c0) * (a0) * (1.-b0) * (1.-c0) * (1.-b1) * (c1)
    geomfac(l, 1, 3, 3) = -1./(a0 * b0 * c0) * (a0) * (1.-b0) * (1.-c0) * (1.-b1) * (1.-c1)

    geomfac(l, 2, 1, 1) = -1./(a0 * b0 * c0) * (1.-a0) * (b0) * (c0) * (a1)
    geomfac(l, 2, 1, 2) = -1./(a0 * b0 * c0) * (1.-a0) * (b0) * (1.-c0) * (a1) * (c1)
    geomfac(l, 2, 1, 3) = -1./(a0 * b0 * c0) * (1.-a0) * (b0) * (1.-c0) * (a1) * (1.-c1)
    geomfac(l, 2, 2, 1) = -1./(a0 * b0 * c0) * (1.-a0) * (1.-b0) * (c0) * (a1) * (b1)
    geomfac(l, 2, 2, 2) = -1./(a0 * b0 * c0) * (1.-a0) * (1.-b0) * (1.-c0) * (a1) * (b1) * (c1)
    geomfac(l, 2, 2, 3) = -1./(a0 * b0 * c0) * (1.-a0) * (1.-b0) * (1.-c0) * (a1) * (b1) * (1.-c1)
    geomfac(l, 2, 3, 1) = -1./(a0 * b0 * c0) * (1.-a0) * (1.-b0) * (c0) * (a1) * (1.-b1)
    geomfac(l, 2, 3, 2) = -1./(a0 * b0 * c0) * (1.-a0) * (1.-b0) * (1.-c0) * (a1) * (1.-b1) * (c1)
    geomfac(l, 2, 3, 3) = -1./(a0 * b0 * c0) * (1.-a0) * (1.-b0) * (1.-c0) * (a1) * (1.-b1) * (1.-c1)

    geomfac(l, 3, 1, 1) = -1./(a0 * b0 * c0) * (1.-a0) * (b0) * (c0) * (1.-a1)
    geomfac(l, 3, 1, 2) = -1./(a0 * b0 * c0) * (1.-a0) * (b0) * (1.-c0) * (1.-a1) * (c1)
    geomfac(l, 3, 1, 3) = -1./(a0 * b0 * c0) * (1.-a0) * (b0) * (1.-c0) * (1.-a1) * (1.-c1)
    geomfac(l, 3, 2, 1) = -1./(a0 * b0 * c0) * (1.-a0) * (1.-b0) * (c0) * (1.-a1) * (b1)
    geomfac(l, 3, 2, 2) = -1./(a0 * b0 * c0) * (1.-a0) * (1.-b0) * (1.-c0) * (1.-a1) * (b1) * (c1)
    geomfac(l, 3, 2, 3) = -1./(a0 * b0 * c0) * (1.-a0) * (1.-b0) * (1.-c0) * (1.-a1) * (b1) * (1.-c1)
    geomfac(l, 3, 3, 1) = -1./(a0 * b0 * c0) * (1.-a0) * (1.-b0) * (c0) * (1.-a1) * (1.-b1)
    geomfac(l, 3, 3, 2) = -1./(a0 * b0 * c0) * (1.-a0) * (1.-b0) * (1.-c0) * (1.-a1) * (1.-b1) * (c1)
    geomfac(l, 3, 3, 3) = -1./(a0 * b0 * c0) * (1.-a0) * (1.-b0) * (1.-c0) * (1.-a1) * (1.-b1) * (1.-c1)
45  continue
  end do
  !$OMP END PARALLEL DO

end subroutine geomfac_intp

!=======================================================================
subroutine findbdy_nointp(dir, nx, ny, nz, nbody, inout, ninner, fcp)
  implicit none
  integer(8), intent(in) :: dir, nx, ny, nz, nbody
  integer(8), intent(in) :: inout(0:nx, 0:ny, 0:nz)
  integer(8), intent(out) :: ninner
  integer(8), intent(out) :: fcp(nbody, 3)

  integer(8) :: i, j, k, istart, jstart, kstart

  istart = 1; jstart = 1; kstart = 1
  if (dir .eq. 1) istart = 2
  if (dir .eq. 2) jstart = 2
  if (dir .eq. 3) kstart = 2

  ninner = 0

  do k = kstart, nz - 1
    do j = jstart, ny - 1
      do i = istart, nx - 1
        if (inout(i, j, k) .eq. 0) then
          ninner = ninner + 1
          fcp(ninner, 1) = i
          fcp(ninner, 2) = j
          fcp(ninner, 3) = k
        end if
      end do
    end do
  end do

end subroutine findbdy_nointp

!=======================================================================
subroutine find_zero_nu_sgs(nx, ny, nz, xm, ym, zm, nzero, inout, t)
!$ USE OMP_LIB
  implicit none

  integer(8), intent(in) :: nx, ny, nz
  real(8), intent(in) :: xm(0:nx), ym(0:ny), zm(0:nz), t
  integer(8), intent(out) :: inout(1:nx - 1, 1:ny - 1, 1:nz - 1)
  integer(8), intent(out) :: nzero

  integer(8) :: i, j, k
  real(8), external :: funcbody

  nzero = 0
  inout = 1

  !$OMP PARALLEL DO REDUCTION(+:NZERO) PRIVATE(I, J, K)
  do k = 1, nz - 1
    do j = 1, ny - 1
      do i = 1, nx - 1
        ! CHECK 7-POINT STENCIL. IF ANY POINT IS INSIDE BODY, ZERO VISCOSITY.
        if ((funcbody(xm(i), ym(j), zm(k), t) .le. 1.0d-10) .or. &
            (funcbody(xm(i - 1), ym(j), zm(k), t) .le. 1.0d-10) .or. &
            (funcbody(xm(i + 1), ym(j), zm(k), t) .le. 1.0d-10) .or. &
            (funcbody(xm(i), ym(j - 1), zm(k), t) .le. 1.0d-10) .or. &
            (funcbody(xm(i), ym(j + 1), zm(k), t) .le. 1.0d-10) .or. &
            (funcbody(xm(i), ym(j), zm(k - 1), t) .le. 1.0d-10) .or. &
            (funcbody(xm(i), ym(j), zm(k + 1), t) .le. 1.0d-10)) then
          nzero = nzero + 1
          inout(i, j, k) = 0
        end if
      end do
    end do
  end do
  !$OMP END PARALLEL DO

end subroutine find_zero_nu_sgs

!=======================================================================
subroutine conjg_intp(nx, ny, nz, cratio, kratio, xm, ym, zm, x, y, z, &
                      iszero, t, cstar, kstar)
!$ USE OMP_LIB
  implicit none

  integer(8), intent(in) :: nx, ny, nz
  real(8), intent(in) :: cratio, kratio
  real(8), intent(in) :: xm(0:nx), ym(0:ny), zm(0:nz)
  real(8), intent(in) :: x(0:nx), y(0:ny), z(0:nz)
  integer(8), intent(in) :: iszero(1:nx - 1, 1:ny - 1, 1:nz - 1)
  real(8), intent(in) :: t

  real(8), intent(out) :: cstar(1:nx - 1, 1:ny - 1, 1:nz - 1)
  real(8), intent(out) :: kstar(1:nx - 1, 1:ny - 1, 1:nz - 1, 6)

  integer(8) :: i, j, k, subc
  real(8) :: fptemp, aa
  real(8), external :: funcbody, fluid_portion

  ! SUBDIVISION LEVEL FOR FLOOD FILL
  subc = 6  ! OPTIMIZED VALUE FROM LEE & HWANG (2019)

  !$OMP PARALLEL DO PRIVATE(I, J, K, FPTEMP, AA)
  do k = 1, nz - 1
    do j = 1, ny - 1
      do i = 1, nx - 1
        ! AA > 0 IF CENTER IS FLUID. AA <= 0 IF CENTER IS SOLID.
        aa = funcbody(xm(i), ym(j), zm(k), t)

        ! HEAT CAPACITY C*
        fptemp = fluid_portion(x(i), x(i + 1), y(j), y(j + 1), z(k), z(k + 1), t, subc)
        cstar(i, j, k) = (1.0d0 - fptemp) * cratio + fptemp * 1.0d0

        ! THERMAL CONDUCTIVITY K* (6 FACES)

        ! EAST FACE (I+1/2) - CENTER (XM) TO (XM+1)
        ! USING 'INTERIM' CELL LOGIC FROM PAPER
        fptemp = fluid_portion(xm(i), xm(i + 1), y(j), y(j + 1), z(k), z(k + 1), t, subc)
        if (aa * funcbody(xm(i + 1), ym(j), zm(k), t) .ge. 0.0d0) then
          ! SAME PHASE
          kstar(i, j, k, 1) = (1.0d0 - fptemp) * kratio + fptemp * 1.0d0
        else
          ! HARMONIC MEAN FOR INTERFACE
          kstar(i, j, k, 1) = kratio / (kratio * fptemp + 1.0d0 * (1.0d0 - fptemp))
        end if

        ! WEST FACE (I-1/2)
        fptemp = fluid_portion(xm(i - 1), xm(i), y(j), y(j + 1), z(k), z(k + 1), t, subc)
        if (aa * funcbody(xm(i - 1), ym(j), zm(k), t) .ge. 0.0d0) then
          kstar(i, j, k, 2) = (1.0d0 - fptemp) * kratio + fptemp * 1.0d0
        else
          kstar(i, j, k, 2) = kratio / (kratio * fptemp + 1.0d0 * (1.0d0 - fptemp))
        end if

        ! NORTH FACE (J+1/2)
        fptemp = fluid_portion(x(i), x(i + 1), ym(j), ym(j + 1), z(k), z(k + 1), t, subc)
        if (aa * funcbody(xm(i), ym(j + 1), zm(k), t) .ge. 0.0d0) then
          kstar(i, j, k, 3) = (1.0d0 - fptemp) * kratio + fptemp * 1.0d0
        else
          kstar(i, j, k, 3) = kratio / (kratio * fptemp + 1.0d0 * (1.0d0 - fptemp))
        end if

        ! SOUTH FACE (J-1/2)
        fptemp = fluid_portion(x(i), x(i + 1), ym(j - 1), ym(j), z(k), z(k + 1), t, subc)
        if (aa * funcbody(xm(i), ym(j - 1), zm(k), t) .ge. 0.0d0) then
          kstar(i, j, k, 4) = (1.0d0 - fptemp) * kratio + fptemp * 1.0d0
        else
          kstar(i, j, k, 4) = kratio / (kratio * fptemp + 1.0d0 * (1.0d0 - fptemp))
        end if

        ! TOP FACE (K+1/2)
        fptemp = fluid_portion(x(i), x(i + 1), y(j), y(j + 1), zm(k), zm(k + 1), t, subc)
        if (aa * funcbody(xm(i), ym(j), zm(k + 1), t) .ge. 0.0d0) then
          kstar(i, j, k, 5) = (1.0d0 - fptemp) * kratio + fptemp * 1.0d0
        else
          kstar(i, j, k, 5) = kratio / (kratio * fptemp + 1.0d0 * (1.0d0 - fptemp))
        end if

        ! BOTTOM FACE (K-1/2)
        fptemp = fluid_portion(x(i), x(i + 1), y(j), y(j + 1), zm(k - 1), zm(k), t, subc)
        if (aa * funcbody(xm(i), ym(j), zm(k - 1), t) .ge. 0.0d0) then
          kstar(i, j, k, 6) = (1.0d0 - fptemp) * kratio + fptemp * 1.0d0
        else
          kstar(i, j, k, 6) = kratio / (kratio * fptemp + 1.0d0 * (1.0d0 - fptemp))
        end if
      end do
    end do
  end do
  !$OMP END PARALLEL DO

end subroutine conjg_intp

!=======================================================================
function fluid_portion(x1, x2, y1, y2, z1, z2, t, div) result(val)
  implicit none
  real(8), intent(in) :: x1, x2, y1, y2, z1, z2, t
  integer(8), intent(in) :: div
  real(8) :: val

  ! LOCAL VARIABLES FOR FLOOD FILL
  real(8) :: dx, dy, dz, xc, yc, zc
  integer(8) :: i, j, k, ii, jj, kk
  integer(8) :: solid_count, start_phase

  ! AUTOMATIC ARRAYS (THREAD-SAFE ON STACK)
  ! PHASE: -1 (UNKNOWN), 0 (SOLID), 1 (FLUID)
  integer(8) :: phase(div, div, div)

  ! STACK FOR RECURSION (ITERATIVE IMPLEMENTATION)
  integer(8) :: stack(div * div * div, 3)
  integer(8) :: stack_ptr

  real(8), external :: funcbody

  dx = (x2 - x1) / dble(div)
  dy = (y2 - y1) / dble(div)
  dz = (z2 - z1) / dble(div)

  ! INITIALIZE
  phase = -1
  stack_ptr = 0
  solid_count = 0

  ! --- STEP 1: INITIATION (CHECK CORNER 1,1,1) ---
  xc = x1 + dx * 0.5d0
  yc = y1 + dy * 0.5d0
  zc = z1 + dz * 0.5d0

  if (funcbody(xc, yc, zc, t) .le. 1.0d-10) then
    start_phase = 0 ! SOLID
  else
    start_phase = 1 ! FLUID
  end if

  phase(1, 1, 1) = start_phase
  stack_ptr = stack_ptr + 1
  stack(stack_ptr, 1) = 1
  stack(stack_ptr, 2) = 1
  stack(stack_ptr, 3) = 1

  ! --- STEP 2: FLOOD FILL ---
  do while (stack_ptr .gt. 0)
    ! POP
    i = stack(stack_ptr, 1)
    j = stack(stack_ptr, 2)
    k = stack(stack_ptr, 3)
    stack_ptr = stack_ptr - 1

    ! CHECK 6 NEIGHBORS
    do kk = k - 1, k + 1
      do jj = j - 1, j + 1
        do ii = i - 1, i + 1
          ! CHECK MANHATTEN DISTANCE = 1 (NEIGHBORS ONLY, NO DIAGONALS, NO SELF)
          if (abs(ii - i) + abs(jj - j) + abs(kk - k) .ne. 1) cycle

          ! BOUNDARY CHECK
          if (ii .lt. 1 .or. ii .gt. div) cycle
          if (jj .lt. 1 .or. jj .gt. div) cycle
          if (kk .lt. 1 .or. kk .gt. div) cycle

          ! IF ALREADY VISITED/IDENTIFIED, SKIP
          if (phase(ii, jj, kk) .ne. -1) cycle

          ! IDENTIFY NEIGHBOR
          xc = x1 + (x2 - x1) / dble(div) * (dble(2 * ii - 1) / 2.0d0)
          yc = y1 + (y2 - y1) / dble(div) * (dble(2 * jj - 1) / 2.0d0)
          zc = z1 + (z2 - z1) / dble(div) * (dble(2 * kk - 1) / 2.0d0)

          if (funcbody(xc, yc, zc, t) .le. 1.0d-10) then
            phase(ii, jj, kk) = 0 ! SOLID
          else
            phase(ii, jj, kk) = 1 ! FLUID
          end if

          ! IF NEIGHBOR HAS SAME PHASE AS START, PUSH TO STACK (CONTINUE FLOOD)
          ! IF DIFFERENT, IT IS A BOUNDARY, SO WE STOP FLOODING THAT PATH.
          if (phase(ii, jj, kk) .eq. start_phase) then
            stack_ptr = stack_ptr + 1
            stack(stack_ptr, 1) = ii
            stack(stack_ptr, 2) = jj
            stack(stack_ptr, 3) = kk
          end if
        end do
      end do
    end do
  end do

  ! --- STEP 3: TERMINATION & COUNTING ---
  ! ANY SUB-CELL WITH PHASE == -1 IS INACCESSIBLE FROM (1,1,1).
  ! ASSUMING SINGLE INTERFACE, THESE MUST BE THE OPPOSITE PHASE OF START_PHASE.
  do k = 1, div
    do j = 1, div
      do i = 1, div
        if (phase(i, j, k) .eq. -1) then
          phase(i, j, k) = 1 - start_phase
        end if

        if (phase(i, j, k) .eq. 0) then
          solid_count = solid_count + 1
        end if
      end do
    end do
  end do

  val = 1.0d0 - dble(solid_count) / dble(div**3)

end function fluid_portion
