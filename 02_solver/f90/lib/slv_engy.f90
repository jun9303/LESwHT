!=======================================================================
        subroutine rhsnlhs_t
!=======================================================================
          use mod_common
          implicit none

          if (iles .eq. 1) call rhssgs_t
          call rhs_t
          if ((ibmon .ne. 0) .and. iconjg .eq. 0) call rhs_ibm_t

          ! BYPASS CONVECTIVE OUTLET FOR SCALING METHOD
          if (xprdic .eq. 0) then
            call convbc_t
            if (ich .ne. 1) call rhsincorpbc_t
          end if

          call lhs_t

          ! BYPASS CONVECTIVE RETRIEVAL FOR SCALING METHOD
          if (xprdic .eq. 0) call retrv_t

          call prdic_adj_t
          call wallbc_t

          return
        end
!=======================================================================
!=======================================================================
        subroutine rhs_t
!=======================================================================
!
!     COMPUTING INTERMEDIATE VELOCITY, U HAT, STEP IN DELTA FORM
!
!     VARIABLES IN COMMON:
!     X, Y, Z         : COORDINATE DIRECTION
!     U, V, W         : VELOCITY FOR X, Y, Z DIRECTION
!     N, S, E, W, C, F: + & - FOR EACH X, Y, Z DIRECTION
!     AN, AL, TAL     : NON-LINEAR, LINEAR, TURBULENT (SGS)
!
!     RHS1(X,Y,Z,L):
!           RHS TERM CONSISTS OF NON-LINEAR/LINEAR/SGS TERMS
!           COMPONENTS FOR IB METHOD WILL BE ADDED IN RHS_IBM SUBROUTINE
!
!-----------------------------------------------------------------------
          use mod_common
          use mod_flowarray, only: u, v, w, p, t, alsgs, alsgs1, rhs1 &
                                   , rk3to, rk3too, nwall_dvm, cstar, kstar
          implicit none
          integer(8) :: i, j, k
          integer(8) :: iplus, iminus, jplus, jminus, kplus, kminus
          real(8) :: omega

!------------ VARIABLES FOR T (TEMPERATURE)
          real(8) :: te, tw, tn, ts, tc, tf
          real(8) :: ant1, ant2, ant3, rk3t
          real(8) :: alt1, alt2, alt3, alt4, alt5, alt6, altx, alty, altz, alt
          real(8) :: taltx, talty, taltz

!-----RHS1 CALCULATION FOR T -----------------
!$OMP PARALLEL DO  &
!$OMP PRIVATE(TE,TW,TN,TS,TC,TF)   &
!$OMP PRIVATE(ANT1,ANT2,ANT3,RK3T) &
!$OMP PRIVATE(ALT1,ALT2,ALT3,ALT4,ALT5,ALT6,ALTX,ALTY,ALTZ,ALT) &
!$OMP PRIVATE(TALTX,TALTY,TALTZ,OMEGA)
          do k = 1, n3m
            do j = 1, n2m
              do i = 1, n1m

                kplus = kpv(k)
                kminus = kmv(k)
                jplus = jpv(j)
                jminus = jmv(j)
                iplus = ipv(i)
                iminus = imv(i)

                ! DECOUPLE PERIODIC WRAP-AROUND FOR TEMPERATURE SCALING
                if (xprdic .eq. 1) then
                  if (i .eq. n1m) iplus = n1
                  if (i .eq. 1) iminus = 0
                end if

                te = 0.5 * (f2fx(i) * t(iplus, j, k) + f2fx(iplus) * t(i, j, k)) * c2cxi(iplus) &
                     * (1.-fixiu(i)) + t(iplus, j, k) * fixiu(i)
                tw = 0.5 * (f2fx(iminus) * t(i, j, k) + f2fx(i) * t(iminus, j, k)) * c2cxi(i) &
                     * (1.-fixil(i)) + t(iminus, j, k) * fixil(i)

                tn = 0.5 * (f2fy(j) * t(i, jplus, k) + f2fy(jplus) * t(i, j, k)) * c2cyi(jplus) &
                     * (1.-fixju(j)) + t(i, jplus, k) * fixju(j)
                ts = 0.5 * (f2fy(jminus) * t(i, j, k) + f2fy(j) * t(i, jminus, k)) * c2cyi(j) &
                     * (1.-fixjl(j)) + t(i, jminus, k) * fixjl(j)

                tc = 0.5 * (f2fz(k) * t(i, j, kplus) + f2fz(kplus) * t(i, j, k)) * c2czi(kplus) &
                     * (1.-fixku(k)) + t(i, j, kplus) * fixku(k)
                tf = 0.5 * (f2fz(kminus) * t(i, j, k) + f2fz(k) * t(i, j, kminus)) * c2czi(k) &
                     * (1.-fixkl(k)) + t(i, j, kminus) * fixkl(k)

                ant1 = (te * u(iplus, j, k) - tw * u(i, j, k)) * f2fxi(i)
                ant2 = (tn * v(i, jplus, k) - ts * v(i, j, k)) * f2fyi(j)
                ant3 = (tc * w(i, j, kplus) - tf * w(i, j, k)) * f2fzi(k)

                ! OMEGA ACTS AS A SWITCH FOR THE CONVECTIVE TERM.
                ! FOR STATIONARY CONJUGATE HEAT TRANSFER (IMOVINGON = 0), VELOCITY INSIDE THE SOLID IS ZERO,
                ! SO WE FORCE OMEGA = 0. TO GUARANTEE PURE CONDUCTION AND SUPPRESS NUMERICAL NOISE.
                ! FOR A MOVING SOLID (IMOVINGON = 1), THE SOLID HAS A RIGID-BODY VELOCITY,
                ! AND THEREFORE ADVECTS ITS OWN TEMPERATURE. WE MUST LEAVE OMEGA = 1.
                if ((iconjg .eq. 1) .and. (nwall_dvm(i, j, k) .eq. 0) .and. (imovingon .eq. 0)) then
                  omega = 0.
                else
                  omega = 1.
                end if

                rk3t = -omega * (ant1 + ant2 + ant3)              ! NON-LINEAR TERM AT K-SUBSTEP

                alt1 = (t(iplus, j, k) - t(i, j, k)) * c2cxi(iplus)
                alt2 = (t(i, j, k) - t(iminus, j, k)) * c2cxi(i)
                alt3 = (t(i, jplus, k) - t(i, j, k)) * c2cyi(jplus)
                alt4 = (t(i, j, k) - t(i, jminus, k)) * c2cyi(j)
                alt5 = (t(i, j, kplus) - t(i, j, k)) * c2czi(kplus)
                alt6 = (t(i, j, k) - t(i, j, kminus)) * c2czi(k)

                altx = (alt1 - alt2) * f2fxi(i)
                alty = (alt3 - alt4) * f2fyi(j)
                altz = (alt5 - alt6) * f2fzi(k)

                if (iconjg .eq. 1) then
                  altx = cstar(i, j, k) * (kstar(i, j, k, 1) * alt1 - kstar(i, j, k, 2) * alt2) * f2fxi(i)
                  alty = cstar(i, j, k) * (kstar(i, j, k, 3) * alt3 - kstar(i, j, k, 4) * alt4) * f2fyi(j)
                  altz = cstar(i, j, k) * (kstar(i, j, k, 5) * alt5 - kstar(i, j, k, 6) * alt6) * f2fzi(k)
                end if

                alt = 1./(re * pr) * (altx + alty + altz)           ! LINEAR TERMS AT K-SUBSTEP

!-----LES
                if (iles .eq. 1) then
                  taltx = f2fxi(i) * &
                          (alsgs1(iplus, j, k, 1) * c2cxi(iplus) * (t(iplus, j, k) - t(i, j, k)) &
                           - alsgs1(i, j, k, 1) * c2cxi(i) * (t(i, j, k) - t(iminus, j, k)))
                  talty = f2fyi(j) * &
                          (alsgs1(i, jplus, k, 2) * c2cyi(jplus) * (t(i, jplus, k) - t(i, j, k)) &
                           - alsgs1(i, j, k, 2) * c2cyi(j) * (t(i, j, k) - t(i, jminus, k)))
                  taltz = f2fzi(k) * &
                          (alsgs1(i, j, kplus, 3) * c2czi(kplus) * (t(i, j, kplus) - t(i, j, k)) &
                           - alsgs1(i, j, k, 3) * c2czi(k) * (t(i, j, k) - t(i, j, kminus)))
                  alt = alt + float(iles) * (taltx + talty + taltz)
                  rk3t = rk3t + float(iles) * rhs1(i, j, k, 4)
                end if
!-----LES

                rhs1(i, j, k, 4) = dt &
                                   * (gamma(msub) * rk3t + ro(msub) * rk3to(i, j, k) &
                                      + 2.*alpha * alt)
                rk3too(i, j, k) = rk3to(i, j, k)
                rk3to(i, j, k) = rk3t

              end do
            end do
          end do

          return
        end
!=======================================================================
!=======================================================================
        subroutine rhs_ibm_t
!=======================================================================
!
!     CALCULATE MOMENTUM FORCING
!
!     OPTION
!     IMOVINGON = 0, STATIONARY BODY => UBODY,VBODY,WBODY FOR TRANSLATIONAL VEL.
!     IMOVINGON = 1, MOVING BODY     => UBD,VBD,WBD IN LICA_CYLINDER.F90
!
!     VARIABLES
!     UTARG,VTARG,WTARG: TARGET VELOCITY TO SATISFY NO-SLIP B.C.
!       FCV   : MOMENTUM FORCING FROM THE TARGET VELOCITY
!     FCVAVG: AVERAGING FORCING VALUES    TO CALCULATE THE FORCE ON A BODY
!     DUDTR : TIME DERIVATIVE OF VELOCITY TO CALCULATE THE FORCE ON A BODY
!               REF. LEE ET AL., 2011, SOURCES OF SPURIOUS FORCE
!               OSCILLATIONS FROM AN IMMERSED BOUNDARY METHOD FOR MOVING
!               -BODY PROBLEMS, J. COMP. PHYS., 230, 2677-2695.
!
!-----------------------------------------------------------------------
          use mod_common
          use mod_flowarray, only: rhs1, u, v, w, t, ifc, jfc, kfc, intpindx, geomfac, &
                                   fcv, fcvavg, dudtr
          implicit none
          integer(8) :: i, j, k, n, ii, jj, kk, ip, jp, kp, ipp, jpp, kpp
          real(8) :: tbody, dti
          real(8) :: ttarg
          real(8) :: rhstmp(n1m, n2m, n3m)

!---- COMPUTE TARGET VELOCITIES & FORCING VALUES AT FORCING POINTS
          tbody = 0.
          dti = 1./dt

!$OMP PARALLEL DO
          do k = 1, n3m
            do j = 1, n2m
              do i = 1, n1m
                rhstmp(i, j, k) = rhs1(i, j, k, 4)
              end do
            end do
          end do

!$OMP PARALLEL DO &
!$OMP PRIVATE(II,JJ,KK,IP,JP,KP,IPP,JPP,KPP,TTARG,TBODY)
          do n = 1, nintp(4)
            ii = ifc(n, 4)
            jj = jfc(n, 4)
            kk = kfc(n, 4)
            ip = ifc(n, 4) + intpindx(n, 4, 1)
            jp = jfc(n, 4) + intpindx(n, 4, 2)
            kp = kfc(n, 4) + intpindx(n, 4, 3)
            ipp = ifc(n, 4) + intpindx(n, 4, 1) * 2
            jpp = jfc(n, 4) + intpindx(n, 4, 2) * 2
            kpp = kfc(n, 4) + intpindx(n, 4, 3) * 2

            if ((ipp .ge. n1) .or. (ipp .le. 0)) call reindex_i(ip, ipp) ! AT SLV_MMTM LIB
            if ((kpp .ge. n3) .or. (kpp .le. 0)) call reindex_k(kp, kpp) ! AT SLV_MMTM LIB
            ttarg = geomfac(n, 4, 0, 0, 0) * tbody &
                    + geomfac(n, 4, 0, 0, 1) * (t(ii, jj, kp) + rhstmp(ii, jj, kp)) &
                    + geomfac(n, 4, 0, 0, 2) * (t(ii, jj, kpp) + rhstmp(ii, jj, kpp)) &
                    + geomfac(n, 4, 0, 1, 0) * (t(ii, jp, kk) + rhstmp(ii, jp, kk)) &
                    + geomfac(n, 4, 0, 1, 1) * (t(ii, jp, kp) + rhstmp(ii, jp, kp)) &
                    + geomfac(n, 4, 0, 1, 2) * (t(ii, jp, kpp) + rhstmp(ii, jp, kpp)) &
                    + geomfac(n, 4, 0, 2, 0) * (t(ii, jpp, kk) + rhstmp(ii, jpp, kk)) &
                    + geomfac(n, 4, 0, 2, 1) * (t(ii, jpp, kp) + rhstmp(ii, jpp, kp)) &
                    + geomfac(n, 4, 0, 2, 2) * (t(ii, jpp, kpp) + rhstmp(ii, jpp, kpp)) &
                    + geomfac(n, 4, 1, 0, 0) * (t(ip, jj, kk) + rhstmp(ip, jj, kk)) &
                    + geomfac(n, 4, 1, 0, 1) * (t(ip, jj, kp) + rhstmp(ip, jj, kp)) &
                    + geomfac(n, 4, 1, 0, 2) * (t(ip, jj, kpp) + rhstmp(ip, jj, kpp)) &
                    + geomfac(n, 4, 1, 1, 0) * (t(ip, jp, kk) + rhstmp(ip, jp, kk)) &
                    + geomfac(n, 4, 1, 1, 1) * (t(ip, jp, kp) + rhstmp(ip, jp, kp)) &
                    + geomfac(n, 4, 1, 1, 2) * (t(ip, jp, kpp) + rhstmp(ip, jp, kpp)) &
                    + geomfac(n, 4, 1, 2, 0) * (t(ip, jpp, kk) + rhstmp(ip, jpp, kk)) &
                    + geomfac(n, 4, 1, 2, 1) * (t(ip, jpp, kp) + rhstmp(ip, jpp, kp)) &
                    + geomfac(n, 4, 1, 2, 2) * (t(ip, jpp, kpp) + rhstmp(ip, jpp, kpp)) &
                    + geomfac(n, 4, 2, 0, 0) * (t(ipp, jj, kk) + rhstmp(ipp, jj, kk)) &
                    + geomfac(n, 4, 2, 0, 1) * (t(ipp, jj, kp) + rhstmp(ipp, jj, kp)) &
                    + geomfac(n, 4, 2, 0, 2) * (t(ipp, jj, kpp) + rhstmp(ipp, jj, kpp)) &
                    + geomfac(n, 4, 2, 1, 0) * (t(ipp, jp, kk) + rhstmp(ipp, jp, kk)) &
                    + geomfac(n, 4, 2, 1, 1) * (t(ipp, jp, kp) + rhstmp(ipp, jp, kp)) &
                    + geomfac(n, 4, 2, 1, 2) * (t(ipp, jp, kpp) + rhstmp(ipp, jp, kpp)) &
                    + geomfac(n, 4, 2, 2, 0) * (t(ipp, jpp, kk) + rhstmp(ipp, jpp, kk)) &
                    + geomfac(n, 4, 2, 2, 1) * (t(ipp, jpp, kp) + rhstmp(ipp, jpp, kp)) &
                    + geomfac(n, 4, 2, 2, 2) * (t(ipp, jpp, kpp) + rhstmp(ipp, jpp, kpp))
            fcv(n, 4) = (ttarg - t(ii, jj, kk) - rhstmp(ii, jj, kk)) * dti
            rhs1(ii, jj, kk, 4) = ttarg - t(ii, jj, kk)
          end do

!$OMP PARALLEL DO &
!$OMP PRIVATE(II,JJ,KK,TTARG,TBODY)
          do n = nintp(4) + 1, nbody(4)
            ii = ifc(n, 4)
            jj = jfc(n, 4)
            kk = kfc(n, 4)
            ttarg = tbody
            fcv(n, 4) = (ttarg - t(ii, jj, kk) - rhs1(ii, jj, kk, 4)) * dti
            rhs1(ii, jj, kk, 4) = ttarg - t(ii, jj, kk)
          end do

!-----COMPUTE AVERAGE FORCING VALUES DURING RK3 STEPS
!$OMP PARALLEL DO
          do n = 1, nbody(4)
            fcvavg(n, 4) = fcvavg(n, 4) + fcv(n, 4)
          end do

          return
        end
!=======================================================================
!=======================================================================
        subroutine convbc_t
!=======================================================================
          use mod_common
          use mod_flowarray, only: u, v, w, t
          implicit none
          integer(8) :: i, j, k
          real(8) :: qin, qout, qratio
          real(8) :: ubar, ucoef, vwcoef
          real(8) :: funcbody

          if (ich .ne. 1) then

!-----CALCULATE INFLUX
            qin = 0.
!$OMP PARALLEL DO REDUCTION(+:QIN)
            do k = 1, n3m
              do j = 1, n2m
                qin = u(1, j, k) * f2fy(j) * f2fz(k) + qin
              end do
            end do
!$OMP PARALLEL DO REDUCTION(+:QIN)
            do k = 1, n3m
              do i = 1, n1m
                qin = (v(i, 1, k) - v(i, n2, k)) * f2fx(i) * f2fz(k) + qin
              end do
            end do
!$OMP PARALLEL DO REDUCTION(+:QIN)
            do j = 1, n2m
              do i = 1, n1m
                qin = (w(i, j, 1) - w(i, j, n3)) * f2fx(i) * f2fy(j) + qin
              end do
            end do

!-----CALCULATE CONVECTIVE VELOCITY UC
            ubar = 0.
!$OMP PARALLEL DO REDUCTION(+:UBAR)
            do k = 1, n3m
              do j = 1, n2m
                ubar = u(n1, j, k) * f2fy(j) * f2fz(k) + ubar
              end do
            end do
            ubar = ubar / (yl * zl)
            ucoef = ubar * dtconst * f2fxi(n1m)
            vwcoef = ubar * dtconst * c2cxi(n1)

            qout = 0.
!$OMP PARALLEL DO REDUCTION(+:QOUT)
            do k = 1, n3m
              do j = 1, n2m
                uout(j, k) = u(n1, j, k) - ucoef * (u(n1, j, k) - u(n1m, j, k))

                if (funcbody(xmp(n1m), ymp(j), zmp(k), time) .le. 1.e-10 .and. imovingon .eq. 0) then
                  ! SOLID PHASE: ZERO-GRADIENT CONDUCTION BOUNDARY
                  tout(j, k) = t(n1m, j, k)
                else
                  ! FLUID PHASE: CONVECTIVE ADVECTION
                  tout(j, k) = t(n1, j, k) - vwcoef * (t(n1, j, k) - t(n1m, j, k))
                end if

                qout = uout(j, k) * f2fy(j) * f2fz(k) + qout
              end do
            end do
            qratio = qin / qout

!-----ADJUST BOUNDARY VELOCITY TO SATISFY GLOBAL MASS CONSERVATION
!$OMP PARALLEL DO
            do k = 1, n3m
              do j = 1, n2m
                ! ONLY SCALE THE FLUID PHASE BY THE MASS FLUX ERROR
                if (funcbody(xmp(n1m), ymp(j), zmp(k), time) .gt. 1.e-10) then
                  tout(j, k) = tout(j, k) * qratio
                end if
                dtout(j, k) = tout(j, k) - t(n1, j, k)
              end do
            end do

            if (yprdic .eq. 1) then
!$OMP PARALLEL DO
              do k = 0, n3
                dtout(n2, k) = dtout(1, k)
                dtout(0, k) = dtout(n2m, k)
              end do
            end if

            if (zprdic .eq. 1) then
!$OMP PARALLEL DO
              do j = 0, n2
                dtout(j, n3) = dtout(j, 1)
                dtout(j, 0) = dtout(j, n3m)
              end do
            end if

          else

            tout = 0.

          end if

          return
        end subroutine convbc_t
!=======================================================================
!=======================================================================
        subroutine rhsincorpbc_t
!=======================================================================
          use mod_common
          use mod_flowarray, only: u, v, w, t, rhs1, alsgs, cstar, kstar
          implicit none
          integer(8) :: i, j, k
          real(8) :: cre, cre2, cs, ks

          if (iles .eq. 1) then
            cre = re * pr
            cre2 = 2.*re * pr
          else
            cre = 0.
            cre2 = 0.
          end if

!$OMP PARALLEL DO
          do k = 1, n3m
            do j = 1, n2m
              if (iconjg .eq. 1) then
                cs = cstar(n1m, j, k)
                ks = (kstar(n1m, j, k, 1) + kstar(n1m, j, k, 2) + kstar(n1m, j, k, 3) &
                      + kstar(n1m, j, k, 4) + kstar(n1m, j, k, 5) + kstar(n1m, j, k, 6)) / 6.
              else
                cs = 1.
                ks = 1.
              end if

              rhs1(n1m, j, k, 4) = rhs1(n1m, j, k, 4) &
                                   - acoef * cs * ks &
                                   * ciu(n1m) * (1.+cre2 * alsgs(n1m, j, k)) * dtout(j, k)
            end do
          end do

          return
        end subroutine rhsincorpbc_t
!=======================================================================
!=======================================================================
        subroutine lhs_t
!=======================================================================
!
!     CALCULATE INTERMEDIATE VELOCITY, U_I HAT THROUGH TDMA
!           IN THIS ROUTINE, COMPUTE STREAMWISE VELOCITY (U HAT)
!
!-----------------------------------------------------------------------
          use mod_common
          use mod_flowarray, only: u, v, w, t, rhs1, alsgs, alsgs1, cstar, kstar
          implicit none
          integer(8) :: i, j, k
          real(8) :: cre, cre2, cs, ks1, ks2
          real(8), dimension(:, :), allocatable :: ai, bi, ci, gi
          real(8), dimension(:, :), allocatable :: aj, bj, cj, gj
          real(8), dimension(:, :), allocatable :: ak, bk, ck, gk

          cre = float(iles) * re * pr
          cre2 = 2.*float(iles) * re * pr

!=====ADI STARTS

          if (n3m .eq. 1) goto 100
!-----Z-DIRECTION
!$OMP PARALLEL &
!$OMP PRIVATE(AK,CK,BK,GK,CS,KS1,KS2)
          allocate (ak(n1, n3), bk(n1, n3), ck(n1, n3), gk(n1, n3))
!$OMP DO
          do j = 1, n2m
            do k = 1, n3m
              do i = 1, n1m
                if (iconjg .eq. 1) then
                  cs = cstar(i, j, k)
                  ks1 = kstar(i, j, k, 6)
                  ks2 = kstar(i, j, k, 5)
                else
                  cs = 1.
                  ks1 = 1.
                  ks2 = 1.
                end if
                ak(i, k) = akuv(k) * (cs * ks1 + cre * alsgs1(i, j, k, 2))
                ck(i, k) = ckuv(k) * (cs * ks2 + cre * alsgs1(i, j, k + 1, 2))
                bk(i, k) = acoefi * pr - ak(i, k) - ck(i, k)
                gk(i, k) = acoefi * pr * rhs1(i, j, k, 4)
              end do
            end do

            if (zprdic .eq. 0) then
              call trdiag3(ak, bk, ck, gk, gk, 1, n3m, 1, n1m)
            else if (zprdic .eq. 1) then
              call trdiag3p(ak, bk, ck, gk, 1, n3m, 1, n1m)  !Z PERIODICITY
            end if

            do k = 1, n3m
              do i = 1, n1m
                rhs1(i, j, k, 4) = gk(i, k)
              end do
            end do
          end do
!$OMP END DO
          deallocate (ak, bk, ck, gk)
!$OMP END PARALLEL

100       continue

!$OMP PARALLEL  &
!$OMP PRIVATE(AJ,CJ,BJ,GJ)  &
!$OMP PRIVATE(AI,CI,BI,GI,CS,KS1,KS2)
          allocate (ai(n2, n1), bi(n2, n1), ci(n2, n1), gi(n2, n1))
          allocate (aj(n1, n2), bj(n1, n2), cj(n1, n2), gj(n1, n2))
!$OMP DO
          do k = 1, n3m

!-----Y-DIRECTION
            do j = 1, n2m
              do i = 1, n1m
                if (iconjg .eq. 1) then
                  cs = cstar(i, j, k)
                  ks1 = kstar(i, j, k, 4)
                  ks2 = kstar(i, j, k, 3)
                else
                  cs = 1.
                  ks1 = 1.
                  ks2 = 1.
                end if
                aj(i, j) = ajuw(j) * (cs * ks1 + cre * alsgs1(i, j, k, 3))
                cj(i, j) = cjuw(j) * (cs * ks2 + cre * alsgs1(i, j + 1, k, 3))
                bj(i, j) = acoefi * pr - aj(i, j) - cj(i, j)
                gj(i, j) = acoefi * pr * rhs1(i, j, k, 4)
              end do
            end do

            call trdiag2(aj, bj, cj, gj, gj, 1, n2m, 1, n1m)

!-----X-DIRECTION
            do i = 1, n1m
              do j = 1, n2m
                if (iconjg .eq. 1) then
                  cs = cstar(i, j, k)
                  ks1 = kstar(i, j, k, 2)
                  ks2 = kstar(i, j, k, 1)
                else
                  cs = 1.
                  ks1 = 1.
                  ks2 = 1.
                end if
                ai(j, i) = aivw(i) * (cs * ks1 + cre * alsgs1(i, j, k, 1))
                ci(j, i) = civw(i) * (cs * ks2 + cre * alsgs1(i + 1, j, k, 1))
                bi(j, i) = acoefi * pr - ai(j, i) - ci(j, i)
                gi(j, i) = acoefi * pr * gj(i, j)
              end do
            end do

            call trdiag1(ai, bi, ci, gi, gi, 1, n1m, 1, n2m)

            do j = 1, n2m
              do i = 1, n1m
                t(i, j, k) = gi(j, i) + t(i, j, k)
              end do
            end do

          end do
!$OMP END DO
          deallocate (ai, bi, ci, gi)
          deallocate (aj, bj, cj, gj)
!$OMP END PARALLEL

          return
        end
!=======================================================================
!=======================================================================
        subroutine retrv_t
!=======================================================================
!
!     ADJUST BOUNDARY VELOCITY TO SATISFY GLOBAL MASS CONSERVATION
!           [RETRIEVE UVW (AT THE EXIT BOUNDARY)]
!     SEE 'CONVBC' SUBROUTINE IN LICA_[BODYNAME, E.G. CYLINDER].F90 FILE
!
!-----------------------------------------------------------------------
          use mod_common
          use mod_flowarray, only: u, v, w, t
          implicit none
          integer(8) :: i, j, k

!$OMP PARALLEL DO
          do k = 1, n3m
            do j = 1, n2m
              t(n1, j, k) = tout(j, k)
            end do
          end do

          return
        end subroutine retrv_t
!=======================================================================
!=======================================================================
        subroutine prdic_adj_t
!=======================================================================
          use mod_common
          use mod_flowarray, only: u, t
          implicit none
          integer(8) :: i, j, k
          real(8) :: t_b_out, m_dot
          real(8) :: t_wall, t_b_in, ratio
          real(8) :: funcbody

! X PERIODICITY WITH DIRICHLET SCALING (RECYCLING METHOD)
          if (xprdic .eq. 1) then

            if (t_inf .eq. 0) then
              t_wall = -1.0d0
            else
              t_wall = 1.0d0
            end if
            t_b_in = 0.0d0

            t_b_out = 0.0d0
            m_dot = 0.0d0

!$OMP PARALLEL DO REDUCTION(+:T_B_OUT, M_DOT)
            do k = 1, n3m
              do j = 1, n2m
                ! CHECK GEOMETRY AT OUTLET (N1M)
                if (funcbody(xmp(n1m), ymp(j), zmp(k), time) .ge. 1.e-10) then
                  m_dot = m_dot + u(n1m, j, k) * f2fy(j) * f2fz(k)
                  t_b_out = t_b_out + t(n1m, j, k) * u(n1m, j, k) * f2fy(j) * f2fz(k)
                end if
              end do
            end do
!$OMP END PARALLEL DO

            t_b_out = t_b_out / m_dot
            ratio = (t_wall - t_b_in) / (t_wall - t_b_out)

!$OMP PARALLEL DO
            do k = 0, n3
              do j = 0, n2
                ! 1. INLET GETS THE SCALED OUTLET PROFILE (SOURCE OF THE DEVELOPED SHAPE)
                t(0, j, k) = t_wall - (t_wall - t(n1m, j, k)) * ratio

                ! 2. OUTLET GHOST CELL USES INVERSELY SCALED INLET PROFILE TO PRESERVE DERIVATIVES
                t(n1, j, k) = t_wall - (t_wall - t(1, j, k)) / ratio
              end do
            end do
!$OMP END PARALLEL DO

          end if

! Y PERIODICITY
          if (yprdic .eq. 1) then
!$OMP PARALLEL DO
            do k = 1, n3
              do i = 0, n1
                t(i, 0, k) = t(i, n2m, k)
                t(i, n2, k) = t(i, 1, k)
              end do
            end do
          end if

! Z PERIODICITY
          if (zprdic .eq. 1) then
!$OMP PARALLEL DO
            do j = 0, n2
              do i = 1, n1
                t(i, j, 0) = t(i, j, n3m)
                t(i, j, n3) = t(i, j, 1)
              end do
            end do
          end if

          return
        end subroutine prdic_adj_t
!=======================================================================
!=======================================================================
        subroutine wallbc_t
!=======================================================================
          use mod_common
          use mod_flowarray, only: t
          implicit none
          integer(8) :: i, j, k
          real(8) :: t_solid
          real(8) :: funcbody

          if (t_inf .eq. 0) then
            t_solid = -1.0d0
          else
            t_solid = 1.0d0
          end if

          if (xprdic .eq. 0) then
!$OMP PARALLEL DO
            do k = 0, n3
              do j = 0, n2
                ! BOTTOM X (I=0): SEPARATE SOLID VS FLUID DIRICHLETS UNCONDITIONALLY
                if (funcbody(xmp(0), ymp(j), zmp(k), time) .le. 1.e-10) then
                  t(0, j, k) = t_solid   ! SOLID BODY
                else
                  t(0, j, k) = 0.0d0     ! FLUID
                end if

                ! TOP X (I=N1): CONVECTIVE OUTLET IS ACTIVELY MAPPED BY RETRV_T.
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

          return
        end subroutine wallbc_t
!=======================================================================
