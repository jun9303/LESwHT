!=======================================================================
!
!     CODEBASE (LICA 2017 VERSION) BY
!     H. CHOI / DEPARTMENT OF MECHANICAL & AEROSPACE ENGINEERING
!     SEOUL NATIONAL UNIVERSITY
!
!=======================================================================
!
!     LESWHT (C) 2026 S. LEE (ORCID: 0000-0002-2063-6298)
!     2018.02.28. UPDATED (DOI:10.1016/J.IJHEATMASSTRANSFER.2019.01.019)
!     2026.02.18. CODE MODERNIZATION
!
!=======================================================================
program solver
  use mod_common
  use mod_flowarray
  implicit none

  real(8) :: omp_get_wtime
  real(8) :: time_prev
  logical :: perturb_applied

  !===================================================================
  !     INITIALIZATION PHASE
  !===================================================================
  call print_real_time()                   ! AT MISC_INIT LIBRARY

  total_time_b = omp_get_wtime()           ! INTRINSIC SUBROUTINE
  call readsettings()                      ! AT MISC_INIT LIBRARY
  call readbcs()                           ! AT MISC_INIT LIBRARY
  call readgeom()                          ! AT MISC_INIT LIBRARY
  call alloinit()                          ! AT MISC_INIT LIBRARY
  call allo_array()                        ! AT MISC_INIT LIBRARY

  ! --- LOAD PRE-PROCESSED DATA ---
  if (ibmon .eq. 1) then
    call ibmpreread()                    ! AT MISC_INIT LIBRARY
    if (iles .eq. 1) then
      call nutzeroread()               ! AT MISC_INIT LIBRARY
    end if
    if (iconjg .eq. 1) call conjgread()  ! AT MISC_INIT LIBRARY
  end if

  if (iles .eq. 1) call sgsfilterinit()    ! AT SGS LIBRARY

  ! --- FIELD LOADING ---
  if (iread .eq. 1) then
    call prefld()                        ! AT MISC_INIT LIBRARY
  else
    call makefld()                       ! AT MISC_INIT LIBRARY
    call writefield()
  end if

  if (ich .eq. 1) call meanpg()             ! AT SLV_MMTM LIBRARY

  ! --- MISCELLANEOUS SETTINGS ---
  call dttimeinit()                        ! AT MISC_INIT LIBRARY
  call rk3coefinit()                       ! AT MISC_INIT LIBRARY
  call lhsinit()                           ! AT SLV_MMTM LIBRARY
  call poisinit()                          ! AT POISS LIBRARY
  call ftrfiles(1)

  cflmax = cflfac
  nv = 101
  nav = 3001
  perturb_applied = .false.
  cdavg_dur = 0.0d0
  cdavg_int = 0.0d0

  !===================================================================
  !     TIME-DEPENDENT CALCULATION (MAIN SOLVER)
  !===================================================================
  do while ((ntime .lt. nend) .and. (time .lt. tend))
    time_begin = omp_get_wtime()         ! INTRINSIC SUBROUTINE
    call stepinit()                      ! AT SLV_AUXI LIBRARY
    subdt = 0.0d0

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !     RK3 SUB-STEP LOOP
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    do msub = 1, 3
      call substepinit(msub)           ! AT SLV_AUXI LIBRARY

      if (ich .eq. 1) then
        call qvolcalc(qvol(1))       ! AT SLV_AUXI LIBRARY
      end if

      sgstime_b(msub) = omp_get_wtime()
      if (iles .eq. 1) then
        call sgscalc()               ! AT SGS LIBRARY
        if (ihtrans .eq. 1) call sgscalc_t() ! AT SGS LIBRARY
      end if
      sgstime_e(msub) = omp_get_wtime()

      if (imovingon .eq. 1) call findforcing() ! AT IBM_BODY (GLOBAL)

      rhsnlhstime_b(msub) = omp_get_wtime()
      call rhsnlhs()                   ! AT SLV_MMTM LIBRARY
      rhsnlhstime_e(msub) = omp_get_wtime()

      allocate (divgsum(0:n1, 0:n2, 0:n3))
      allocate (phi(0:n1, 0:n2, 0:n3))
      divgsum = 0.0d0
      phi = 0.0d0

      poisstime_b(msub) = omp_get_wtime()
      call divgs(divgsum)              ! AT SLV_CONT LIBRARY
      call poisson(phi, divgsum)       ! AT POISS LIBRARY
      poisstime_e(msub) = omp_get_wtime()

      if (ich .eq. 1) then
        call qvolcalc(qvol(2))       ! AT SLV_AUXI LIBRARY
        call qvolcorr()              ! AT SLV_AUXI LIBRARY
      end if

      call ucalc(phi)

      if (ich .eq. 1) then
        call qvolcalc(qvol(0))       ! AT SLV_AUXI LIBRARY
      end if

      call pcalc(phi, divgsum)

      if (ich .eq. 1) call meanpg()    ! AT SLV_AUXI LIBRARY
      if (ibmon .eq. 1) call lagforce()

      if (ihtrans .eq. 1) then
        if (ich .eq. 1) call tvolcalc(tvol(1)) ! AT SLV_AUXI LIBRARY
        call rhsnlhs_t()                       ! AT SLV_ENGY LIBRARY
        if (ich .eq. 1) call tvolcalc(tvol(2)) ! AT SLV_AUXI LIBRARY
        if (ich .eq. 1) call tvolcorr()        ! AT SLV_AUXI LIBRARY
        call tcalc()                           ! AT SLV_AUXI LIBRARY
        if (ich .eq. 1) call tvolcalc(tvol(0)) ! AT SLV_AUXI LIBRARY
      end if

      deallocate (divgsum, phi)
    end do
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    time_prev = time
    time = time + dt

    if ((iread .eq. 0) .and. (eps_ptr .ne. 0.0d0) .and. (.not. perturb_applied)) then
      if ((time_prev .lt. ptb_tst) .and. (time .ge. ptb_tst)) then
        call addperturb()
        call prdic_adj_uvw(0)
        perturb_applied = .true.
      end if
    end if

    call convergence_check()
    if ((ibmon .eq. 1) .and. (masson .eq. 1)) call masscheck()

    if (mod(ntime, npin) .eq. 0) then
      ihist = ihist + 1
      call writehistory()
    end if

    if ((ntrace .gt. 0) .and. (mod(m, ntr) .eq. 0)) call tracer()
    if (ibmon .eq. 1) then
      call draglift()
    end if
    if (ihtrans .eq. 1) call calc_boundary_heat_flux()

    if (mod(ntime, nprint) .eq. 0) call writefield()
    if (iavg .eq. 1) call field_avg()

    time_end = omp_get_wtime()
    call writeftrtime()

  end do

  !===================================================================
  !     FINISH AND CLEANUP
  !===================================================================
  call ftrfiles(0)

  if (mod(ntime, nprint) .ne. 0) call writefield()
  if (iavg .eq. 1) call field_avg_finalize()

  call cpu_time(total_time_e)

  call print_real_time()
  call deallo_array()

  stop
end program solver
