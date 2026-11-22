!=======================================================================
      subroutine print_real_time
!=======================================================================
        implicit none
        character*8 :: date
        character*10 :: now
        character*5 :: zone
        integer(8) :: vals(8)

        call date_and_time(date, now, zone, vals)

        write (*, 101) vals(1), vals(2), vals(3), vals(5), vals(6), vals(7)
        write (*, *) ''
101     format(' @ 'i0.4, '-', i0.2, '-', i0.2, ' ', i0.2, ':', i0.2, ':', i0.2)

        return
      end subroutine print_real_time
!=======================================================================
      subroutine readsettings
!=======================================================================
        use mod_common
        implicit none
        integer(8) :: n
        character*10 :: dummy
        character(len=256) :: sgsline
        integer :: ios

        open (10, file='settings.input')
        read (10, *) dummy
        read (10, *) re, pr, gr, grdir, tflu_i, tsol_i
        read (10, *) dummy
        read (10, *) ireset, iread, iavg, ipzero, eps_ptr, udrv_i
        read (10, *) dummy
        read (10, *) nend, nprint, npriavg, npin
        read (10, *) dummy
        read (10, *) tend, ptb_tst, avg_tst
        read (10, *) dummy
        read (10, *) idtopt, dt_size, cflfac
        read (10, *) dummy
        read (10, *) resid1, nlev, nbli, ioldv, mgitr, imgsor, wwsor
        read (10, *) dummy
        read (10, '(A)') sgsline
        ios = 0
        read (sgsline, *, iostat=ios) iles, insmdl, itemdl, idvmon, csgsts, csgshf, filter, idelta_sgs
        if (ios .ne. 0) then
          read (sgsline, *, iostat=ios) iles, insmdl, itemdl, idvmon, csgsts, csgshf, filter
          idelta_sgs = 1
        end if
        read (10, *) dummy
        read (10, *) ibmon, masson, imovingon, ihtrans, idecomp
!       THE MOVING-BODY PATH (FINDFORCING EACH STEP) IS NOT USABLE: IT PASSES ARRAYS WITH OTHER SHAPES
!       THAN THE IBM ROUTINES DECLARE (COORDINATES WITHOUT THE -1:N+1 RANGE, FCP/INTPINDX/GEOMFAC LEADING
!       DIMENSIONS) AND NO BODY-VELOCITY FUNCTIONS (UBD, VBD, WBD) ARE DEFINED. STOP HERE, NOT AT STEP 1.
        if (imovingon .eq. 1) then
          write (*, *) 'IMOVINGON = 1 IS NOT AVAILABLE IN THIS VERSION (SEE MISC_INIT.F90); SET IMOVINGON = 0.'
          stop
        end if
        read (10, *) dummy
        read (10, *) gridfile
        read (10, *) prev_fld
          call resolve_legacy_output_path(gridfile)
          call resolve_legacy_output_path(prev_fld)
        read (10, *) dummy
        read (10, *) ntrace, ntr
        if (ntrace .gt. 0) then
          do n = 1, ntrace
            read (10, *) trpts(n, 1), trpts(n, 2), trpts(n, 3)
          end do
        end if
        close (10)

        write (*, *) '========= SETTINGS ========='
        write (*, *) ''
        write (*, 101) re, pr, gr
        write (*, 112) tflu_i, tsol_i
        write (*, 102) ireset, iread, iavg, ipzero, eps_ptr, udrv_i
        write (*, 103) nend, nprint, npriavg, npin
        write (*, 111) tend, ptb_tst, avg_tst
        write (*, 104) idtopt, dt, cflfac
        write (*, 105) resid1, nlev, nbli, ioldv, mgitr, imgsor, wwsor
        write (*, 106) iles, insmdl, itemdl, idvmon, csgsts, csgshf, filter, idelta_sgs
        write (*, 107) ibmon, masson, imovingon, ihtrans, idecomp
        write (*, 108) gridfile
        write (*, 109) prev_fld
        write (*, *) ''

        if (ntrace .gt. 0) then
          write (*, *) '========= TRACE POSITION ========='
          write (*, *) ''
          do n = 1, ntrace
            write (*, 110) n, trpts(n, 1), trpts(n, 2), trpts(n, 3)
          end do
        end if

101     format('  RE=', es11.3, '  PR=', es11.3, '  GR=', es11.3)
112     format('  TFLU_I=', es13.5, '  TSOL_I=', es13.5)
102     format('  IRESET=', i5, '  IREAD=', i5, '  IAVG=', i5, '  IPZERO=', i5, '  EPS_PTR=', f7.3, '  UDRV_I=', f7.3)
103     format('  NEND=', i10, '  NPRINT=', i8, '  NPRIAVG=', i8, '  NPIN=', i5)
111     format('  TEND=', es13.5, '  PTB_TST=', es13.5, '  AVG_TST=', es13.5)
104     format('  IDTOPT=', i5, '  DT=', es13.5, '  CFLFAC=', f11.3)
105     format('  RESID=', es12.4, '  NLEV=', i5, '  NBLI=', i5, '  IOLDV=', i5, '  MGITR=', i5, '  IMGSOR=', i5, '  WWSOR=', f7.3)
106     format('  ILES=', i2, '  INSMDL=', i2, '  ITEMDL=', i2, '  IDVMON=', i2, '  CSGSTS='f12.4, '  CSGSHF=', f12.4, '  IFILTER=', i5, '  IDELTA_SGS=', i3)
107     format('  IBMON=', i2, '  MASSON=', i2, '  IMOVINGON=', i2, '  IHTRANS=', i2, '  IDECOMP=', i2)
  108     format('  GRIDFILE=', a)
  109     format('  PREV_FLD=', a)
110     format(i5, 3i6)

        return
      end subroutine readsettings
!=======================================================================
!=======================================================================
      subroutine readbcs
!=======================================================================
        use mod_common
        implicit none
        character*10 :: dummy

        open (10, file='boundary.input')
        read (10, *) dummy
        read (10, *) bc_ybtm, bc_ytop
        read (10, *) dummy
        read (10, *) bc_zbtm, bc_ztop
        read (10, *) dummy
        read (10, *) bc_t_ybtm, bc_t_ytop
        read (10, *) dummy
        read (10, *) val_t_ybtm, val_t_ytop
        read (10, *) dummy
        read (10, *) bc_t_zbtm, bc_t_ztop
        read (10, *) dummy
        read (10, *) val_t_zbtm, val_t_ztop
        read (10, *) dummy
        read (10, *) ich, iconjg
        close (10)

          open (11, file=trim(output_dir_ibmpre)//'/ibmpre_prdic.bin')
        read (11, *) xprdic, yprdic, zprdic, iintp
        close (11)

        if (xprdic .eq. 1) then
          bc_xbtm = 4 ! 4:   PERIODIC B.C.
          bc_xtop = 4
        else
          bc_xbtm = 2 ! 2:   VELOCITY INLET B.C.
          bc_xtop = 3 ! 3:   CONVECTIVE OUTLET B.C.
        end if

        if (yprdic .eq. 1) then
          bc_ybtm = 4
          bc_ytop = 4
        end if

        if (zprdic .eq. 1) then
          bc_zbtm = 4
          bc_ztop = 4
        end if

!       FAR FIELD (BC = 1): NO FLOW THROUGH THE BOUNDARY AND ZERO GRADIENT OF THE TANGENTIAL VELOCITY.
!       THESE FLAGS SWITCH THE TANGENTIAL GHOSTS (UCALC) AND THE BOUNDARY TERMS OF RHS AND LHSU/V/W
!       TO ZERO GRADIENT; WITHOUT THEM THE GHOSTS KEPT THEIR INITIAL VALUES.
        jub = 0
        jwb = 0
        jut = 0
        jwt = 0
        kub = 0
        kvb = 0
        kut = 0
        kvt = 0
        if (bc_ybtm .eq. 1) then
          jub = 1
          jwb = 1
        end if
        if (bc_ytop .eq. 1) then
          jut = 1
          jwt = 1
        end if
        if (bc_zbtm .eq. 1) then
          kub = 1
          kvb = 1
        end if
        if (bc_ztop .eq. 1) then
          kut = 1
          kvt = 1
        end if

!       CONSTANT FLOW RATE NEEDS A PERIODIC X: WITH AN INLET THE INLET VELOCITY SETS THE FLOW RATE, AND
!       CONVBC CLOSES THE OUTLET FOR ICH = 1 (NO SOLUTION OF THE PRESSURE EQUATION).
        if ((ich .eq. 1) .and. (xprdic .eq. 0)) then
          write (*, *) ' --- ICH=1 (CONSTANT FLOW RATE) REQUIRES X PERIODICITY. WITH AN INLET, USE ICH=0.'
          stop
        end if

        call check_idecomp_input

        write (*, *) '========= BOUNDARY CONDITIONS ========='
        write (*, *) ''
        write (*, *) '(0: WALL; 1: FARF; 2: VELIN; 3: CONVOUT; 4: PRDIC)'
        write (*, 101) bc_xbtm, bc_xtop
        write (*, 102) bc_ybtm, bc_ytop
        write (*, 103) bc_zbtm, bc_ztop
        write (*, *) ''
        write (*, 104) ich, iconjg, iintp
        write (*, *) ''

101     format('  BC_XBTM=', i2, '  BC_XTOP=', i2)
102     format('  BC_YBTM=', i2, '  BC_YTOP=', i2)
103     format('  BC_ZBTM=', i2, '  BC_ZTOP=', i2)
104     format('  ICH=', i2, '  ICONJG=', i2, '  IINTP=', i2)

        return
      end subroutine readbcs
!=======================================================================
!=======================================================================
      subroutine check_idecomp_input
!=======================================================================
        use mod_common
        implicit none
        real(8), parameter :: eps_bc = 1.0d-12

        if (idecomp .ne. 1) return

        if (xprdic .ne. 1) then
          write (*, *) ' --- IDECOMP=1 REQUIRES X PERIODICITY (XPRDIC=1).'
          stop
        end if

        if ((bc_ybtm .ne. 0) .or. (bc_ytop .ne. 0)) then
          write (*, *) ' --- IDECOMP=1 REQUIRES BOTH Y MOMENTUM BOUNDARIES TO BE WALLS.'
          stop
        end if

        if ((bc_t_ybtm .ne. 0) .or. (bc_t_ytop .ne. 0)) then
          write (*, *) ' --- IDECOMP=1 REQUIRES Y THERMAL BC TO BE DIRICHLET AT BOTH ENDS.'
          stop
        end if

        if ((dabs(val_t_ybtm) .gt. eps_bc) .or. (dabs(val_t_ytop) .gt. eps_bc)) then
          write (*, *) ' --- IDECOMP=1 REQUIRES HOMOGENEOUS Y THERMAL DIRICHLET: VAL_T_YBTM=VAL_T_YTOP=0.'
          stop
        end if

        if ((zprdic .eq. 0) .and. (bc_zbtm .eq. 0) .and. (bc_ztop .eq. 0)) then
          if ((bc_t_zbtm .ne. 0) .or. (bc_t_ztop .ne. 0)) then
            write (*, *) ' --- IDECOMP=1 WITH Z WALLS REQUIRES Z THERMAL BC TO BE DIRICHLET AT BOTH ENDS.'
            stop
          end if
          if ((dabs(val_t_zbtm) .gt. eps_bc) .or. (dabs(val_t_ztop) .gt. eps_bc)) then
            write (*, *) ' --- IDECOMP=1 WITH Z WALLS REQUIRES VAL_T_ZBTM=VAL_T_ZTOP=0.'
            stop
          end if
        end if

        return
      end subroutine check_idecomp_input
!=======================================================================
!=======================================================================
      subroutine readgeom
!=======================================================================
        use mod_common
        implicit none
        integer(8) :: i, j, k

        open (10, file=gridfile)
        read (10, *) n1, n2, n3
        read (10, *) n1m, n2m, n3m
        read (10, *) xl, yl, zl

        call allo(n1, n2, n3)

        read (10, *) (x(i), i=1, n1)
        read (10, *) (y(j), j=1, n2)
        read (10, *) (z(k), k=1, n3)
        read (10, *) (ipv(i), i=1, n1m)
        read (10, *) (jpv(j), j=1, n2m)
        read (10, *) (kpv(k), k=1, n3m)
        read (10, *) (imv(i), i=1, n1m)
        read (10, *) (jmv(j), j=1, n2m)
        read (10, *) (kmv(k), k=1, n3m)
        read (10, *) (fixil(i), i=1, n1m)
        read (10, *) (fixjl(j), j=1, n2m)
        read (10, *) (fixkl(k), k=1, n3m)
        read (10, *) (fixiu(i), i=1, n1m)
        read (10, *) (fixju(j), j=1, n2m)
        read (10, *) (fixku(k), k=1, n3m)
        read (10, *) (c2cx(i), i=0, n1)
        read (10, *) (c2cy(j), j=0, n2)
        read (10, *) (c2cz(k), k=0, n3)
        read (10, *) (c2cxi(i), i=0, n1)
        read (10, *) (c2cyi(j), j=0, n2)
        read (10, *) (c2czi(k), k=0, n3)
        read (10, *) (f2fx(i), i=0, n1)
        read (10, *) (f2fy(j), j=0, n2)
        read (10, *) (f2fz(k), k=0, n3)
        read (10, *) (f2fxi(i), i=0, n1)
        read (10, *) (f2fyi(j), j=0, n2)
        read (10, *) (f2fzi(k), k=0, n3)
        read (10, *) (xmp(i), i=0, n1)
        read (10, *) (ymp(j), j=0, n2)
        read (10, *) (zmp(k), k=0, n3)

        close (10)

        return
      end subroutine readgeom
!=======================================================================
!=======================================================================
      subroutine alloinit
!=======================================================================
        use mod_common
        use mod_flowarray
        implicit none
        integer(8) :: n

          open (11, file=trim(output_dir_ibmpre)//'/ibmpre_fcpts.bin')
        read (11, *) nintp(1), nintp(2), nintp(3)
        read (11, *) nbody(1), nbody(2), nbody(3)
        do n = 1, 3
          nbody(n) = nbody(n) + nintp(n)
        end do
        close (11)

        nzero = 0
        if (iles .eq. 1) then
          open (11, file=trim(output_dir_ibmpre)//'/ibmpre_nutzero.bin')
          read (11, *) nzero
          close (11)
        end if

        if (ihtrans .eq. 1) then
            open (11, file=trim(output_dir_ibmpre)//'/ibmpre_fcpts_t.bin')
          read (11, *) nintp(4)
          read (11, *) nbody(4)
          nbody(4) = nbody(4) + nintp(4)
          close (11)
        end if

        return
      end subroutine alloinit
!=======================================================================
!=======================================================================
      subroutine allo_array
!=======================================================================
        use mod_common
        use mod_flowarray
        implicit none

        if (ihtrans .eq. 0) then
          call basic_allo()
        else
          call thermal_allo()
        end if

        if (iles .eq. 1) then
          call les_allo()
          if (ihtrans .eq. 1) call les_thermal_allo()
        else
!         DNS: THE MOMENTUM AND ENERGY SOLVERS STILL READ THE EDDY VISCOSITY AND DIFFUSIVITY (MULTIPLIED
!         BY ZERO), SO THEY ARE ALLOCATED AND SET TO ZERO
          allocate (nusgs(0:n1, 0:n2, 0:n3), nusgs1(0:n1, 0:n2, 0:n3, 3))
          nusgs = 0.
          nusgs1 = 0.
          if (ihtrans .eq. 1) call les_thermal_allo()
        end if

!       CSTAR/KSTAR/OMEGA_DVM ARE USED BY THERMAL SOLVER PATHS EVEN WHEN ICONJG=0 (CSTAR = KSTAR = 1,
!       SET IN CONJG_ALLO, UNLESS ICONJG=1 LOADS IBMPRE_CONJG.BIN)
        if ((iconjg .eq. 1) .or. (ihtrans .eq. 1)) then
          call conjg_allo()
        end if

        if (iavg .eq. 1) call avg_allo()

        return
      end subroutine allo_array
!=======================================================================
!=======================================================================
      subroutine ibmpreread
!=======================================================================
        use mod_common
        use mod_flowarray
        implicit none
        integer(8) :: n, l, dummy
        integer(8) :: i, j, k, itmp

          open (11, file=trim(output_dir_ibmpre)//'/ibmpre_fcpts.bin')
        read (11, *) dummy ! NINTP, ALREADY READ IN ALLOINIT SUBROUTINE
        read (11, *) dummy ! NBODY, ALREADY READ IN ALLOINIT SUBROUTINE
        do n = 1, nbody(1)
          read (11, *) ifc(n, 1), jfc(n, 1), kfc(n, 1)
        end do
        do n = 1, nbody(2)
          read (11, *) ifc(n, 2), jfc(n, 2), kfc(n, 2)
        end do
        do n = 1, nbody(3)
          read (11, *) ifc(n, 3), jfc(n, 3), kfc(n, 3)
        end do
        do n = 1, nintp(1)
          read (11, *) intpindx(n, 1, 1), intpindx(n, 1, 2), intpindx(n, 1, 3)
        end do
        do n = 1, nintp(2)
          read (11, *) intpindx(n, 2, 1), intpindx(n, 2, 2), intpindx(n, 2, 3)
        end do
        do n = 1, nintp(3)
          read (11, *) intpindx(n, 3, 1), intpindx(n, 3, 2), intpindx(n, 3, 3)
        end do
        do n = 1, nintp(1)
          read (11, *) (((geomfac(n, 1, i, j, k), k=0, 2), j=0, 2), i=0, 2)
        end do
        do n = 1, nintp(2)
          read (11, *) (((geomfac(n, 2, i, j, k), k=0, 2), j=0, 2), i=0, 2)
        end do
        do n = 1, nintp(3)
          read (11, *) (((geomfac(n, 3, i, j, k), k=0, 2), j=0, 2), i=0, 2)
        end do
        close (11)

        if (ihtrans .eq. 1) then
            open (11, file=trim(output_dir_ibmpre)//'/ibmpre_fcpts_t.bin')
          read (11, *) dummy ! NINTP, ALREADY READ IN ALLOINIT SUBROUTINE
          read (11, *) dummy ! NBODY, ALREADY READ IN ALLOINIT SUBROUTINE
          do n = 1, nbody(4)
            read (11, *) ifc(n, 4), jfc(n, 4), kfc(n, 4)
          end do
          do n = 1, nintp(4)
            read (11, *) intpindx(n, 4, 1), intpindx(n, 4, 2), intpindx(n, 4, 3)
          end do
          do n = 1, nintp(4)
            read (11, *) (((geomfac(n, 4, i, j, k), k=0, 2), j=0, 2), i=0, 2)
          end do
        end if

        write (*, *) '========= LOADING IBM-PREPROCESSING DATA ========='
        write (*, *) ''
        write (*, 35) nintp(1), nintp(2), nintp(3)
        write (*, 36) nbody(1), nbody(2), nbody(3)
        if (ihtrans .eq. 1) write (*, 37) nintp(4), nbody(4)
        write (*, *) ''
35      format('  NINTP_U :', i12, '  NINTP_V :', i12, '  NINTP_W :', i12)
36      format('  NBODY_U :', i12, '  NBODY_V :', i12, '  NBODY_W :', i12)
37      format('  NINTP_T :', i12, '  NBODY_T :', i12)

        call check_stencils()

        return
      end subroutine ibmpreread
!=======================================================================
!=======================================================================
      subroutine check_stencils
!=======================================================================
!       AN INTERPOLATION STENCIL MUST STAY INSIDE A NON-PERIODIC DIRECTION (INDICES 1 .. N-1), OR THE
!       IBM FORCING READS OUTSIDE THE ARRAYS. FILES FROM AN OLDER PREPROCESSOR CAN HOLD SUCH STENCILS
!       WHEN A BODY IS TANGENT TO, OR WITHIN TWO CELLS OF, A DOMAIN END. STOP AND ASK FOR A RERUN.
        use mod_common
        use mod_flowarray, only: ifc, jfc, kfc, intpindx
        implicit none
        integer(8) :: n, l, lmax, nbad, ipp, jpp, kpp

        lmax = 3
        if (ihtrans .eq. 1) lmax = 4
        nbad = 0
        do l = 1, lmax
          do n = 1, nintp(l)
            ipp = ifc(n, l) + 2 * intpindx(n, l, 1)
            jpp = jfc(n, l) + 2 * intpindx(n, l, 2)
            kpp = kfc(n, l) + 2 * intpindx(n, l, 3)
            if ((xprdic .eq. 0) .and. ((ipp .lt. 1) .or. (ipp .gt. n1m))) nbad = nbad + 1
            if ((yprdic .eq. 0) .and. ((jpp .lt. 1) .or. (jpp .gt. n2m))) nbad = nbad + 1
            if ((zprdic .eq. 0) .and. ((kpp .lt. 1) .or. (kpp .gt. n3m))) nbad = nbad + 1
          end do
        end do

        if (nbad .gt. 0) then
          write (*, '(A,I10,A)') ' IBM DATA: ', nbad, ' INTERPOLATION STENCILS LEAVE A NON-PERIODIC DIRECTION.'
          write (*, *) ' RERUN 01_PRE_PROCESSOR (PREPROCESSING.PY) WITH THIS VERSION OF IBM_BODY.F90.'
          stop
        end if

        return
      end subroutine check_stencils
!=======================================================================
!=======================================================================
      subroutine nutzeroread
!=======================================================================
        use mod_common
        use mod_flowarray
        implicit none
        integer(8) :: i, j, k
        integer(8) :: n

          open (14, file=trim(output_dir_ibmpre)//'/ibmpre_nutzero.bin')
        read (14, *) nzero
        if (nzero .gt. 0) then
          read (14, *) (inz(n), n=1, nzero)
          read (14, *) (jnz(n), n=1, nzero)
          read (14, *) (knz(n), n=1, nzero)
        end if
        close (14)

          open (15, file=trim(output_dir_ibmpre)//'/ibmpre_wallfdvm.bin')
        read (15, *) (((nwall_dvm(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
        close (15)

        write (*, *) '========= LOADING LES-PREPROCESSING DATA ========='
        write (*, *) ''
        write (*, 30) nzero
        write (*, *) ''
30      format('  # OF NUTZERO PTS = ', i10)

        return
      end subroutine nutzeroread
!=======================================================================
!=======================================================================
      subroutine conjgread
!=======================================================================
        use mod_common
        use mod_flowarray
        implicit none
        integer(8) :: i, j, k
        integer(8) :: l
        logical :: has_omega_mask

          open (15, file=trim(output_dir_ibmpre)//'/ibmpre_conjg.bin')
        read (15, *) (((cstar(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
        read (15, *) ((((kstar(i, j, k, l), i=1, n1m), j=1, n2m), k=1, n3m), l=1, 6)
        close (15)

          inquire (file=trim(output_dir_ibmpre)//'/ibmpre_omega.bin', exist=has_omega_mask)
        if (has_omega_mask) then
            open (16, file=trim(output_dir_ibmpre)//'/ibmpre_omega.bin')
          read (16, *) (((omega_dvm(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
          close (16)
        else
          write (*, *) ' --- WARNING: ibmpre_omega.bin NOT FOUND. USING DEFAULT OMEGA=1.'
          omega_dvm = 1
        end if

        write (*, *) '========= LOADING CONJUGATE H.TRANS DATA ========='
        write (*, *) ''

        return
      end subroutine conjgread
!=======================================================================
!=======================================================================
      subroutine makefld
!=======================================================================
!$      USE OMP_LIB
        use mod_common
        use mod_flowarray, only: u, v, w, p, t, qmass
        implicit none
        integer(8) :: i, j, k
        real(8) :: funcbody
        real(8) :: flowarea(n1), fl, fl_s, adj
        real(8) :: xi, ylo, yhi, eta
        logical :: isfl(0:n2)
        real(8), parameter :: eps_pr = 1.0d-12

        ihist = 0
        time = 0.

        u = 0.
        v = 0.
        w = 0.
        p = 0.

        flowarea = 0.

        write (*, *) '========= MAKING INITIAL FIELD ========='
        write (*, *) ''

!       ICH = 0: UNIFORM UDRV_I IN THE FLUID. ICH = 1, 2: START FROM REST (SEE BELOW).
!$OMP PARALLEL DO PRIVATE(XI, J, K, ISFL, YLO, YHI, ETA)
        do i = 0, n1    ! START FROM 0 TO COVER THE INLET GHOST CELL
          ! X IS ALLOCATED 1:N1, SO THE GHOST FACE I = 0 IS PLACED ONE CELL BEFORE X(1)
          if (i .eq. 0) then
            xi = 2.0d0 * x(1) - x(2)
          else
            xi = x(i)
          end if
          do k = 0, n3
            do j = 0, n2
              isfl(j) = funcbody(xi, ymp(j), zmp(k), time) .ge. 1.e-10
            end do
            do j = 0, n2
              ! APPLY PROFILE ONLY IN THE FLUID DOMAIN
              if (isfl(j)) then
                if (ich .eq. 0) then
                  u(i, j, k) = udrv_i
                elseif ((ich .eq. 1) .or. (ich .eq. 2)) then
                  ! THIS BRANCH STARTS ICH = 1/2 RUNS FROM REST AND LETS THE DRIVING FORCE BUILD THE FLOW.
                  ! A PARABOLIC START WOULD BE 1.5 * UB * (1 - ETA**2) ACROSS THE FLUID GAP IN Y:
                  ! CALL YFLUIDGAP(ISFL, J, YLO, YHI); ETA = (YMP(J) - 0.5*(YLO + YHI)) / (0.5*(YHI - YLO))
                  u(i, j, k) = 1.d-12
                else
                  u(i, j, k) = udrv_i
                end if
              else
                u(i, j, k) = 0.0d0
              end if
            end do
          end do
        end do
!$OMP END PARALLEL DO

        if (bc_ybtm .eq. 0) then
          do i = 1, n1
            do k = 0, n3
              u(i, 0, k) = 0.
              u(i, 1, k) = 0.
            end do
          end do
        end if

        if (bc_ytop .eq. 0) then
          do i = 1, n1
            do k = 0, n3
              u(i, n2m, k) = 0.
              u(i, n2, k) = 0.
            end do
          end do
        end if

        if (bc_zbtm .eq. 0) then
          do i = 1, n1
            do j = 0, n2
              u(i, j, 0) = 0.
              u(i, j, 1) = 0.
            end do
          end do
        end if

        if (bc_ztop .eq. 0) then
          do i = 1, n1
            do j = 0, n2
              u(i, j, n3m) = 0.
              u(i, j, n3) = 0.
            end do
          end do
        end if

        if (ihtrans .eq. 1) then
          t = tsol_i

!$OMP PARALLEL DO
          do i = 0, n1
            do j = 1, n2m
              do k = 0, n3
                if (funcbody(xmp(i), ymp(j), zmp(k), time) .ge. 1.e-10) then
                  if (idecomp .eq. 1) then
                    if (dabs(pr) .gt. eps_pr) then
                      t(i, j, k) = u(i, j, k) / pr
                    else
                      t(i, j, k) = u(i, j, k)
                    end if
                  else
                    t(i, j, k) = tflu_i  ! FLUID DOMAIN
                  end if
                end if
              end do
            end do
          end do
!$OMP END PARALLEL DO

          ! --- APPLY INITIAL THERMAL BOUNDARY CONDITIONS ---
          if (xprdic .eq. 0) then
!$OMP PARALLEL DO
            do k = 0, n3
              do j = 0, n2
                if (funcbody(xmp(0), ymp(j), zmp(k), time) .ge. 1.e-10) then
                  t(0, j, k) = tflu_i
                else
                  t(0, j, k) = tsol_i
                end if
              end do
            end do
!$OMP END PARALLEL DO
          end if

          if (yprdic .eq. 0) then
!$OMP PARALLEL DO
            do k = 0, n3
              do i = 0, n1
                if (bc_t_ybtm .eq. 0) then
                  t(i, 0, k) = val_t_ybtm
                else if (bc_t_ybtm .eq. 1) then
                  t(i, 0, k) = t(i, 1, k) - val_t_ybtm / c2cyi(1)
                end if

                if (bc_t_ytop .eq. 0) then
                  t(i, n2, k) = val_t_ytop
                else if (bc_t_ytop .eq. 1) then
                  t(i, n2, k) = t(i, n2m, k) + val_t_ytop / c2cyi(n2)
                end if
              end do
            end do
!$OMP END PARALLEL DO
          end if

          if (zprdic .eq. 0) then
!$OMP PARALLEL DO
            do j = 0, n2
              do i = 0, n1
                if (bc_t_zbtm .eq. 0) then
                  t(i, j, 0) = val_t_zbtm
                else if (bc_t_zbtm .eq. 1) then
                  t(i, j, 0) = t(i, j, 1) - val_t_zbtm / c2czi(1)
                end if

                if (bc_t_ztop .eq. 0) then
                  t(i, j, n3) = val_t_ztop
                else if (bc_t_ztop .eq. 1) then
                  t(i, j, n3) = t(i, j, n3m) + val_t_ztop / c2czi(n3)
                end if
              end do
            end do
!$OMP END PARALLEL DO
          end if
          ! -------------------------------------------------
        end if

        return
      end subroutine makefld
!=======================================================================
!=======================================================================
      subroutine yfluidgap(isfl, j, ylo, yhi)
!=======================================================================
!       FLUID GAP IN Y THAT CONTAINS CELL J OF A COLUMN (ISFL = FLUID FLAGS OF ITS CELLS 0..N2):
!       YLO / YHI ARE THE LOWER FACE OF ITS FIRST FLUID CELL AND THE UPPER FACE OF ITS LAST ONE
!       (A BODY SURFACE ON THE GRID, OR A DOMAIN WALL). A GHOST ROW USES THE GAP OF THE ROW NEXT TO IT.
        use mod_common
        implicit none
        logical :: isfl(0:n2)
        integer(8) :: j, jlo, jhi
        real(8) :: ylo, yhi

        jlo = min(max(j, 1_8), n2m)
        jhi = jlo
        do while (jlo .gt. 1)
          if (.not. isfl(jlo - 1)) exit
          jlo = jlo - 1
        end do
        do while (jhi .lt. n2m)
          if (.not. isfl(jhi + 1)) exit
          jhi = jhi + 1
        end do
        ylo = y(jlo)
        yhi = y(jhi + 1)

        return
      end subroutine yfluidgap
!=======================================================================
!=======================================================================
      subroutine prefld
!=======================================================================
!$      USE OMP_LIB
        use mod_common
        use mod_flowarray, only: u, v, w, t, p
        implicit none
        integer(8) :: i, j, k
        integer(8) :: nn1, nn2, nn3
        real(8) :: rre, ppr, ggr

        write (*, *) '========= LOADING PREVIOUS FIELD ========='
        write (*, *) ''

        open (12, file=prev_fld, form='unformatted')
        read (12) nn1, nn2, nn3, rre, ppr, ggr

        if ((nn1 .ne. n1) .or. (nn2 .ne. n2) .or. (nn3 .ne. n3)) then
          write (*, *) ' --- PREVIOUS FIELD DOES NOT MATCH WITH THE GRID SYSTEM.'
          close (12)
          stop

        else

          read (12) ihist, m, time, dt
          read (12) xprdic, yprdic, zprdic
          read (12) bc_xbtm, bc_xtop, bc_ybtm, bc_ytop, bc_zbtm, bc_ztop
          read (12) ich, iconjg
          read (12) (((u(i, j, k), i=1, n1), j=0, n2), k=0, n3)
          read (12) (((v(i, j, k), i=0, n1), j=1, n2), k=0, n3)
          read (12) (((w(i, j, k), i=0, n1), j=0, n2), k=1, n3)
          read (12) (((p(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
          if (ihtrans .eq. 1) then
            t = tsol_i
            read (12) (((t(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
          end if
          close (12)

          if (ipzero .eq. 1) p = 0.

          write (*, *) ''
          write (*, *) ' ----------- PREVIOUS FIELD INFORMATION -----------'
          write (*, 100) prev_fld
          write (*, 101) rre, ppr, ggr
          write (*, 102) nn1, nn2, nn3
          write (*, 103) ihist, m, time, dt
          write (*, *) ''
  100       format(' PREVIOUS FIELD LOCATION : ', a)
101       format(' RE = ', es12.3, ' PR = ', es12.3, ' GR = ', es12.3)
102       format(' N1 = ', i12, ' N2 = ', i12, ' N3 = ', i12)
103       format(' IHIST = ', i9, ' M = ', i9, ' TIME = ', f10.5, ' DT = ', f12.8)

          if (xprdic .eq. 1) then
!$OMP PARALLEL DO
            do j = 0, n2
              do k = 0, n3
                u(0, j, k) = u(n1m, j, k)
                u(n1, j, k) = u(1, j, k)
              end do
            end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
            do j = 1, n2
              do k = 0, n3
                v(0, j, k) = v(n1m, j, k)
                v(n1, j, k) = v(1, j, k)
              end do
            end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
            do j = 0, n2
              do k = 1, n3
                w(0, j, k) = w(n1m, j, k)
                w(n1, j, k) = w(1, j, k)
              end do
            end do
!$OMP END PARALLEL DO

            if (ihtrans .eq. 1) then
!$OMP PARALLEL DO
              do j = 0, n2
                do k = 0, n3
                  t(0, j, k) = t(n1m, j, k)
                  t(n1, j, k) = t(1, j, k)
                end do
              end do
!$OMP END PARALLEL DO
            end if
          end if

          if (yprdic .eq. 1) then
!$OMP PARALLEL DO
            do k = 0, n3
              do i = 1, n1
                u(i, 0, k) = u(i, n2m, k)
                u(i, n2, k) = u(i, 1, k)
              end do
            end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
            do k = 0, n3
              do i = 0, n1
                v(i, 0, k) = v(i, n2m, k)
                v(i, n2, k) = v(i, 1, k)
              end do
            end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
            do k = 1, n3
              do i = 0, n1
                w(i, 0, k) = w(i, n2m, k)
                w(i, n2, k) = w(i, 1, k)
              end do
            end do
!$OMP END PARALLEL DO

            if (ihtrans .eq. 1) then
!$OMP PARALLEL DO
              do k = 0, n3
                do i = 0, n1
                  t(i, 0, k) = t(i, n2m, k)
                  t(i, n2, k) = t(i, 1, k)
                end do
              end do
!$OMP END PARALLEL DO
            end if

          end if

          if (zprdic .eq. 1) then
!$OMP PARALLEL DO
            do i = 1, n1
              do j = 0, n2
                u(i, j, 0) = u(i, j, n3m)
                u(i, j, n3) = u(i, j, 1)
              end do
            end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
            do i = 0, n1
              do j = 1, n2
                v(i, j, 0) = v(i, j, n3m)
                v(i, j, n3) = v(i, j, 1)
              end do
            end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
            do i = 0, n1
              do j = 0, n2
                w(i, j, 0) = w(i, j, n3m)
                w(i, j, n3) = w(i, j, 1)
              end do
            end do
!$OMP END PARALLEL DO

            if (ihtrans .eq. 1) then
!$OMP PARALLEL DO
              do i = 0, n1
                do j = 0, n2
                  t(i, j, 0) = t(i, j, n3m)
                  t(i, j, n3) = t(i, j, 1)
                end do
              end do
!$OMP END PARALLEL DO
            end if

          end if

        end if

        return
      end subroutine prefld
!=======================================================================
!=======================================================================
      subroutine addperturb
!=======================================================================
!$      USE OMP_LIB
        use mod_common
        use mod_flowarray, only: u, v, w
        implicit none
        real(8) :: ptb, flowarea, pertb_rate, adj
        real(8) :: flowarea_g, pertb_rate_g
        real(8) :: flowarea_g_k, pertb_rate_g_k
        real(8) :: flowarea_g_pln(max(n1m, n2m, n3m)), pertb_rate_g_pln(max(n1m, n2m, n3m))
        integer(8) :: i, j, k
        real(8) :: funcbody
        real(8), parameter :: ret = 180.0d0
        ! THE PERTURBATION SITS NEAR THE UPPER WALL OF THE FLUID GAP IN Y THAT CONTAINS THE CELL (SEE
        ! YFLUIDGAP): ETA IS THE GAP COORDINATE (-1 .. 1), WD THE WALL DISTANCE IN HALF-GAP UNITS. FOR A
        ! CHANNEL WITH WALLS AT Y = -1 AND 1, ETA = Y AND WD = 1 - Y.
        logical :: isfl(0:n2)
        real(8) :: ylo, yhi, eta, wd

        ! --- PERTURB U (GLOBAL MEAN CORRECTION) ---
        flowarea_g = 0.0d0
        pertb_rate_g = 0.0d0

!$OMP PARALLEL DO PRIVATE(PTB, FLOWAREA_G_K, PERTB_RATE_G_K, J, K, ISFL, YLO, YHI, ETA, WD)
        do i = 1, n1m
          flowarea_g_k = 0.0d0
          pertb_rate_g_k = 0.0d0
          do k = 1, n3m
            do j = 0, n2
              isfl(j) = funcbody(x(i), ymp(j), zmp(k), time) .gt. 1.e-10
            end do
            do j = 1, n2m
              if (isfl(j)) then
                flowarea_g_k = flowarea_g_k + f2fy(j) * f2fz(k)
                call yfluidgap(isfl, j, ylo, yhi)
                eta = (ymp(j) - 0.5d0 * (ylo + yhi)) / (0.5d0 * (yhi - ylo))
                wd = (yhi - ymp(j)) / (0.5d0 * (yhi - ylo))
                if ((eta .gt. .5d0) .and. (eta .lt. .99d0)) then
                  ptb = eps_ptr * udrv_i * cos(zmp(k)/zl*acos(-1.d0)) &
                        * wd * ret                                     &
                        * exp(-.01d0*(wd*ret)**2.d0 + .5d0)
                  u(i, j, k) = u(i, j, k) + ptb
                  pertb_rate_g_k = pertb_rate_g_k + ptb * f2fy(j) * f2fz(k)
                end if
              end if
            end do
          end do
          flowarea_g_pln(i) = flowarea_g_k
          pertb_rate_g_pln(i) = pertb_rate_g_k
        end do
!$OMP END PARALLEL DO
        do i = 1, n1m
          flowarea_g = flowarea_g + flowarea_g_pln(i)
          pertb_rate_g = pertb_rate_g + pertb_rate_g_pln(i)
        end do

        if (flowarea_g .gt. 0.0d0) then
          adj = pertb_rate_g / flowarea_g
!$OMP PARALLEL DO
          do i = 1, n1m
            do j = 1, n2m
              do k = 1, n3m
                if (funcbody(x(i), ymp(j), zmp(k), time) .gt. 1.e-10) then
                  u(i, j, k) = u(i, j, k) - adj
                end if
              end do
            end do
          end do
!$OMP END PARALLEL DO
        end if

        ! --- PERTURB W (GLOBAL MEAN CORRECTION) ---
        flowarea_g = 0.0d0
        pertb_rate_g = 0.0d0
!$OMP PARALLEL DO PRIVATE(PTB, FLOWAREA_G_K, PERTB_RATE_G_K, I, J, ISFL, YLO, YHI, ETA, WD)
        do k = 1, n3m
          flowarea_g_k = 0.0d0
          pertb_rate_g_k = 0.0d0
          do i = 1, n1m
            do j = 0, n2
              isfl(j) = funcbody(xmp(i), ymp(j), z(k), time) .gt. 1.e-10
            end do
            do j = 1, n2m
              if (isfl(j)) then
                flowarea_g_k = flowarea_g_k + f2fx(i) * f2fy(j)
                call yfluidgap(isfl, j, ylo, yhi)
                eta = (ymp(j) - 0.5d0 * (ylo + yhi)) / (0.5d0 * (yhi - ylo))
                wd = (yhi - ymp(j)) / (0.5d0 * (yhi - ylo))
                if ((eta .gt. .5d0) .and. (eta .lt. .99d0)) then
                  ptb = eps_ptr * sin(xmp(i)/xl*acos(-1.d0)) &
                        * wd * ret                           &
                        * exp(-.01d0*(wd*ret)**2.d0)
                  w(i, j, k) = w(i, j, k) + ptb
                  pertb_rate_g_k = pertb_rate_g_k + ptb * f2fx(i) * f2fy(j)
                end if
              end if
            end do
          end do
          flowarea_g_pln(k) = flowarea_g_k
          pertb_rate_g_pln(k) = pertb_rate_g_k
        end do
!$OMP END PARALLEL DO
        do k = 1, n3m
          flowarea_g = flowarea_g + flowarea_g_pln(k)
          pertb_rate_g = pertb_rate_g + pertb_rate_g_pln(k)
        end do

        if (flowarea_g .gt. 0.0d0) then
          adj = pertb_rate_g / flowarea_g
!$OMP PARALLEL DO
          do k = 1, n3m
            do i = 1, n1m
              do j = 1, n2m
                if (funcbody(xmp(i), ymp(j), z(k), time) .gt. 1.e-10) then
                  w(i, j, k) = w(i, j, k) - adj
                end if
              end do
            end do
          end do
!$OMP END PARALLEL DO
        end if

        write (*, *) '--- ADDING SINUSOIDAL PERTURBATION'
        write (*, 201) eps_ptr * 100.
        write (*, *) ''
201     format(' --- MAX. SIZE OF PERTURBATION = ', f5.1, '% OF U_BULK')

        return
      end subroutine addperturb
!=======================================================================
!=======================================================================
      subroutine addperturb_isotropic
!=======================================================================
!$      USE OMP_LIB
        use mod_common
        use mod_flowarray, only: u, v, w
        implicit none
        real(8) :: ptb_array(1:n1, 1:n2, 1:n3)
        real(8) :: flowarea_k, pertb_rate_k
        real(8) :: flowarea_pln(max(n1m, n2m, n3m)), pertb_rate_pln(max(n1m, n2m, n3m))
        integer(8) :: i, j, k
        real(8) :: funcbody
        real(8) :: flowarea, pertb_rate, adj

        ! --- PERTURB U ---
        call random_number(ptb_array)
        ptb_array = (ptb_array - .5d0) * 2.d0 * eps_ptr

        do i = 1, n1m
          flowarea = 0.0d0
          pertb_rate = 0.0d0
!$OMP PARALLEL DO PRIVATE(FLOWAREA_K, PERTB_RATE_K)
          do j = 1, n2m
            flowarea_k = 0.0d0
            pertb_rate_k = 0.0d0
            do k = 1, n3m
              if (funcbody(x(i), ymp(j), zmp(k), time) .gt. 1.e-10) then
                u(i, j, k) = u(i, j, k) + ptb_array(i, j, k)
                flowarea_k = flowarea_k + f2fy(j) * f2fz(k)
                pertb_rate_k = pertb_rate_k + ptb_array(i, j, k) * f2fy(j) * f2fz(k)
              end if
            end do
            flowarea_pln(j) = flowarea_k
            pertb_rate_pln(j) = pertb_rate_k
          end do
!$OMP END PARALLEL DO
          do j = 1, n2m
            flowarea = flowarea + flowarea_pln(j)
            pertb_rate = pertb_rate + pertb_rate_pln(j)
          end do
          if (flowarea .gt. 0.0d0) then
            adj = pertb_rate / flowarea
!$OMP PARALLEL DO
            do j = 1, n2m
              do k = 1, n3m
                if (funcbody(x(i), ymp(j), zmp(k), time) .gt. 1.e-10) then
                  u(i, j, k) = u(i, j, k) - adj
                end if
              end do
            end do
!$OMP END PARALLEL DO
          end if
        end do

        ! --- PERTURB V ---
        call random_number(ptb_array)
        ptb_array = (ptb_array - .5d0) * 2.d0 * eps_ptr

        do j = 1, n2m
          flowarea = 0.0d0
          pertb_rate = 0.0d0
!$OMP PARALLEL DO PRIVATE(FLOWAREA_K, PERTB_RATE_K)
          do k = 1, n3m
            flowarea_k = 0.0d0
            pertb_rate_k = 0.0d0
            do i = 1, n1m
               if (funcbody(xmp(i), y(j), zmp(k), time) .gt. 1.e-10) then
                v(i, j, k) = v(i, j, k) + ptb_array(i, j, k)
                flowarea_k = flowarea_k + f2fz(k) * f2fx(i)
                pertb_rate_k = pertb_rate_k + ptb_array(i, j, k) * f2fz(k) * f2fx(i)
              end if
            end do
            flowarea_pln(k) = flowarea_k
            pertb_rate_pln(k) = pertb_rate_k
          end do
!$OMP END PARALLEL DO
          do k = 1, n3m
            flowarea = flowarea + flowarea_pln(k)
            pertb_rate = pertb_rate + pertb_rate_pln(k)
          end do
          if (flowarea .gt. 0.0d0) then
            adj = pertb_rate / flowarea
!$OMP PARALLEL DO
            do k = 1, n3m
              do i = 1, n1m
                if (funcbody(xmp(i), y(j), zmp(k), time) .gt. 1.e-10) then
                  v(i, j, k) = v(i, j, k) - adj
                end if
              end do
            end do
!$OMP END PARALLEL DO
          end if
        end do

        ! --- PERTURB W ---
        call random_number(ptb_array)
        ptb_array = (ptb_array - .5d0) * 2.d0 * eps_ptr

        do k = 1, n3m
          flowarea = 0.0d0
          pertb_rate = 0.0d0
!$OMP PARALLEL DO PRIVATE(FLOWAREA_K, PERTB_RATE_K)
          do i = 1, n1m
            flowarea_k = 0.0d0
            pertb_rate_k = 0.0d0
            do j = 1, n2m
               if (funcbody(xmp(i), ymp(j), z(k), time) .gt. 1.e-10) then
                w(i, j, k) = w(i, j, k) + ptb_array(i, j, k)
                flowarea_k = flowarea_k + f2fx(i) * f2fy(j)
                pertb_rate_k = pertb_rate_k + ptb_array(i, j, k) * f2fx(i) * f2fy(j)
              end if
            end do
            flowarea_pln(i) = flowarea_k
            pertb_rate_pln(i) = pertb_rate_k
          end do
!$OMP END PARALLEL DO
          do i = 1, n1m
            flowarea = flowarea + flowarea_pln(i)
            pertb_rate = pertb_rate + pertb_rate_pln(i)
          end do
          if (flowarea .gt. 0.0d0) then
            adj = pertb_rate / flowarea
!$OMP PARALLEL DO
            do i = 1, n1m
              do j = 1, n2m
                if (funcbody(xmp(i), ymp(j), z(k), time) .gt. 1.e-10) then
                  w(i, j, k) = w(i, j, k) - adj
                end if
              end do
            end do
!$OMP END PARALLEL DO
          end if
        end do

        write (*, *) '--- ADDING ISOTROPIC PERTURBATION FROM THE UNIFORM DIST.'
        write (*, 201) eps_ptr * 100.
        write (*, *) ''
201     format(' --- MAX. SIZE OF PERTURBATION = ', f5.1, '% OF U_BULK')

        return
      end subroutine addperturb_isotropic
!=======================================================================
!=======================================================================
      subroutine dttimeinit
!=======================================================================
        use mod_common
        implicit none

        if (idtopt .eq. 0) dt = dt_size ! FOR CONSTANT USAGE PURPOSE
        if ((idtopt .ne. 0) .and. (iread .eq. 0)) dt = dt_size
        ! FOR INITAL USAGE PURPOSE
        if (ireset .eq. 1) then
          ihist = 0
          time = 0.
        end if
        ntime = 0

        if (iavg .eq. 1) then
          timeinit = avg_tst
          ihistinit = ihist
          npriavg_count = 0
          avg_started = .false.
          ihistavg_start = ihist
          if (iread .ne. 1) then
              open (2999, file=trim(output_dir_field_avg)//'/fav_manifest.dat', status='replace')
            close (2999)
          end if
        end if

        return
      end subroutine dttimeinit
!=======================================================================
!=======================================================================
      subroutine rk3coefinit
!=======================================================================
        use mod_common
        implicit none

        gamma(1) = 8./15.
        gamma(2) = 5./12.
        gamma(3) = 3./4.
        ro(1) = 0.
        ro(2) = -17./60.
        ro(3) = -5./12.

        return
      end subroutine rk3coefinit
!=======================================================================
!=======================================================================
      subroutine convergence_check
!=======================================================================
        use mod_common
        use mod_flowarray, only: u, v, w, qmass
        implicit none
        integer(8) :: i, j, k
        real(8) :: dvg11, dvg12, dvg13, dvg1

        dvmax = 0.

!$OMP PARALLEL DO                        &
!$OMP PRIVATE(DVG11,DVG12,DVG13,DVG1)    &
!$OMP REDUCTION(MAX:DVMAX)
        do k = 1, n3m
          do j = 1, n2m
            do i = 1, n1m
              dvg11 = (u(i + 1, j, k) - u(i, j, k)) * f2fxi(i)
              dvg12 = (v(i, j + 1, k) - v(i, j, k)) * f2fyi(j)
              dvg13 = (w(i, j, k + 1) - w(i, j, k)) * f2fzi(k)
              if (masson .eq. 1) then
                dvg1 = dvg11 + dvg12 + dvg13 - qmass(i, j, k)
              else
                dvg1 = dvg11 + dvg12 + dvg13
              end if
              dvmax = amax1(abs(dvg1), dvmax)
            end do
          end do
        end do

        return
      end subroutine convergence_check
!=======================================================================
!=======================================================================
      subroutine masscheck
!=======================================================================
!
!     CALCULATE MAXIMUM MASS SOURCE/SINK IN IBM
!
!-----------------------------------------------------------------------
        use mod_common
        use mod_flowarray, only: qmass
        implicit none
        integer(8) :: i, j, k

        qmmax = 0.

!$OMP PARALLEL DO &
!$OMP REDUCTION(MAX:QMMAX)
        do k = 1, n3m
          do j = 1, n2m
            do i = 1, n1m
              qmmax = amax1(qmmax, abs(qmass(i, j, k)))
            end do
          end do
        end do

        return
      end subroutine masscheck
!=======================================================================
      subroutine deallo_array
!=======================================================================
        use mod_common
        use mod_flowarray
        implicit none

        if (ihtrans .eq. 0) then
          call basic_deallo()
        else
          call thermal_deallo()
        end if

        if (iles .eq. 1) then
          call les_deallo()
          if (ihtrans .eq. 1) call les_thermal_deallo()
        else
          deallocate (nusgs, nusgs1)
          if (ihtrans .eq. 1) call les_thermal_deallo()
        end if

        if (allocated(cstar)) call conjg_deallo()

        if (iavg .eq. 1) call avg_deallo()

        return
      end subroutine deallo_array
!=======================================================================
