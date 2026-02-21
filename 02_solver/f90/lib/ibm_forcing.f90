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
subroutine findforcing()
  use mod_common
  use mod_flowarray
  implicit none

  ! LOCAL TEMPORARY VARIABLES
  integer(8) :: ninner_u, ninner_v, ninner_w, ninner_t
  integer(8), allocatable :: fcp_temp(:, :)

  ! ALLOCATE TEMPORARY ARRAY FOR FORCING POINTS (FCP)
  ! SIZE IT SAFELY TO TOTAL GRID POINTS TO HANDLE ANY BODY SIZE
  allocate (fcp_temp(n1 * n2 * n3, 3))

  !-------------------------------------------------------------------
  ! 1. U-VELOCITY COMPONENT (FACE X)
  !    GRID: X (FACE), YMP (CENTER), ZMP (CENTER)
  !-------------------------------------------------------------------
  call find_inout(n1, n2, n3, x, ymp, zmp, nbody(1), inout(:, :, :, 1), time)

  call findbdy_intp(1_8, n1, n2, n3, nbody(1), inout(:, :, :, 1), &
                    fixiu, fixju, fixku, fixil, fixjl, fixkl, &
                    nintp(1), ninner_u, fcp_temp, intptype(:, 1), intpindx(:, 1, :))

  ! COPY FCP DATA TO MODULE ARRAYS
  ifc(1:nbody(1), 1) = fcp_temp(1:nbody(1), 1)
  jfc(1:nbody(1), 1) = fcp_temp(1:nbody(1), 2)
  kfc(1:nbody(1), 1) = fcp_temp(1:nbody(1), 3)

  call geomfac_intp(n1, n2, n3, x, ymp, zmp, nintp(1), &
                    nbody(1), fcp_temp, intpindx(:, 1, :), geomfac(:, 1, :, :, :), time)

  !-------------------------------------------------------------------
  ! 2. V-VELOCITY COMPONENT (FACE Y)
  !    GRID: XMP (CENTER), Y (FACE), ZMP (CENTER)
  !-------------------------------------------------------------------
  call find_inout(n1, n2, n3, xmp, y, zmp, nbody(2), inout(:, :, :, 2), time)

  call findbdy_intp(2_8, n1, n2, n3, nbody(2), inout(:, :, :, 2), &
                    fixiu, fixju, fixku, fixil, fixjl, fixkl, &
                    nintp(2), ninner_v, fcp_temp, intptype(:, 2), intpindx(:, 2, :))

  ifc(1:nbody(2), 2) = fcp_temp(1:nbody(2), 1)
  jfc(1:nbody(2), 2) = fcp_temp(1:nbody(2), 2)
  kfc(1:nbody(2), 2) = fcp_temp(1:nbody(2), 3)

  call geomfac_intp(n1, n2, n3, xmp, y, zmp, nintp(2), &
                    nbody(2), fcp_temp, intpindx(:, 2, :), geomfac(:, 2, :, :, :), time)

  !-------------------------------------------------------------------
  ! 3. W-VELOCITY COMPONENT (FACE Z)
  !    GRID: XMP (CENTER), YMP (CENTER), Z (FACE)
  !-------------------------------------------------------------------
  call find_inout(n1, n2, n3, xmp, ymp, z, nbody(3), inout(:, :, :, 3), time)

  call findbdy_intp(3_8, n1, n2, n3, nbody(3), inout(:, :, :, 3), &
                    fixiu, fixju, fixku, fixil, fixjl, fixkl, &
                    nintp(3), ninner_w, fcp_temp, intptype(:, 3), intpindx(:, 3, :))

  ifc(1:nbody(3), 3) = fcp_temp(1:nbody(3), 1)
  jfc(1:nbody(3), 3) = fcp_temp(1:nbody(3), 2)
  kfc(1:nbody(3), 3) = fcp_temp(1:nbody(3), 3)

  call geomfac_intp(n1, n2, n3, xmp, ymp, z, nintp(3), &
                    nbody(3), fcp_temp, intpindx(:, 3, :), geomfac(:, 3, :, :, :), time)

  !-------------------------------------------------------------------
  ! 4. TEMPERATURE (SCALAR) - OPTIONAL
  !    GRID: XMP, YMP, ZMP (ALL CENTER)
  !-------------------------------------------------------------------
  if (ihtrans .eq. 1) then
    call find_inout(n1, n2, n3, xmp, ymp, zmp, nbody(4), inout(:, :, :, 4), time)

    call findbdy_intp(4_8, n1, n2, n3, nbody(4), inout(:, :, :, 4), &
                      fixiu, fixju, fixku, fixil, fixjl, fixkl, &
                      nintp(4), ninner_t, fcp_temp, intptype(:, 4), intpindx(:, 4, :))

    ifc(1:nbody(4), 4) = fcp_temp(1:nbody(4), 1)
    jfc(1:nbody(4), 4) = fcp_temp(1:nbody(4), 2)
    kfc(1:nbody(4), 4) = fcp_temp(1:nbody(4), 3)

    call geomfac_intp(n1, n2, n3, xmp, ymp, zmp, nintp(4), &
                      nbody(4), fcp_temp, intpindx(:, 4, :), geomfac(:, 4, :, :, :), time)
  end if

  deallocate (fcp_temp)

end subroutine findforcing
