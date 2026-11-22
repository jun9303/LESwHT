!=======================================================================
       subroutine stepinit
!=======================================================================
         use mod_common
         use mod_flowarray
         implicit none
         real(8) :: dtcfl, dt_evt, eps_t
         integer(8) :: icfl, jcfl, kcfl
         logical :: nonfin
         real(8), parameter :: cfl_blowup = 100.0d0

         ntime = ntime + 1
         m = ntime
         fcvavg = 0.
         if (ich .ne. 0) pmiavg = 0.

         if (idtopt .eq. 0) dt = dt_size

         call cfl(cflmax, icfl, jcfl, kcfl, nonfin) ! CALCULATE CFL NUMBER

!        STOP AFTER A BLOW-UP: A NON-FINITE VELOCITY, OR A CFL NUMBER ABOVE CFL_BLOWUP WITH THE LAST TIME STEP.
!        WITHOUT THIS THE CFL CONTROL SHRINKS DT TOWARD ZERO AND THE RUN GOES ON UNTIL NEND.
         if (nonfin .or. cflmax .gt. cfl_blowup) then
           write (*, *) ''
           if (nonfin) then
             write (*, '(A,I10,A,F12.5)') ' SOLUTION DIVERGED: NON-FINITE VELOCITY AFTER STEP', ntime - 1, &
               ', TIME', time
           else
             write (*, '(A,ES11.3,A,I10,A,F12.5,A,3I6)') ' SOLUTION DIVERGED: CFL =', cflmax, ' AFTER STEP', &
               ntime - 1, ', TIME', time, ', CELL', icfl, jcfl, kcfl
           end if
           write (*, *) 'CHECK DT_SIZE / CFLFAC, THE BOUNDARY CONDITIONS AND THE POISSON SOLVER OUTPUT ABOVE.'
           flush (6)
           error stop 'SOLUTION DIVERGED'
         end if

         if (idtopt .ne. 0 .and. cflmax .ne. 0.) then
           dtcfl = dmin1(dt * cflfac / cflmax, dt * (0.80 + 0.20 * cflfac / cflmax))
           if (idtopt .eq. 1) dt = dtcfl
         end if

         eps_t = 1.0d-12 * dmax1(1.0d0, dabs(time))
         dt_evt = dt

         if (ptb_tst .gt. (time + eps_t)) then
           dt_evt = dmin1(dt_evt, ptb_tst - time)
         end if
         if (avg_tst .gt. (time + eps_t)) then
           dt_evt = dmin1(dt_evt, avg_tst - time)
         end if
         if (tend .gt. (time + eps_t)) then
           dt_evt = dmin1(dt_evt, tend - time)
         end if

         if (dt_evt .lt. dt) dt = dt_evt

         if (cflmax .gt. (cflfac * 1.1)) then
           print *, ' '
           write (*, 310) ntime, cflmax, time
         else
           print *, ' '
           write (*, 320) ntime, time, dt
         end if
310      format(i15, '   CFL NUMBER EXCEED GIVEN CFL LIMIT :', es18.5, ' AT ', f12.5)
320      format('--------------------------', i6, '  TIME=', f10.5, '  DT=', f12.8)

         return
       end subroutine stepinit
!=======================================================================
!=======================================================================
       subroutine cfl(cflm, icfl, jcfl, kcfl, nonfin)
!=======================================================================
!
!     CALCULATE THE MAXIMUM CFL NUMBER OF FLOW FIELD
!     (HTTP://EN.WIKIPEDIA.ORG/WIKI/COURANT-FRIEDRICHS-LEWY_CONDITION)
!
!     CFL#=U_I*DT/DX_I
!     CFLI=(UC/DX+VC/DY+WC/DZ)*DT  (C FOR CELL CENTER)
!     CFLMPT: INDEX WHERE THE CFL# IS MAXIMUM
!     NONFIN: .TRUE. IF A CELL HAS A NON-FINITE (NAN OR INFINITE) CFL NUMBER, I.E. VELOCITY
!
!-----------------------------------------------------------------------
         use mod_common
         use mod_flowarray, only: u, v, w
         implicit none
         real(8) :: cflm
         real(8) :: cfli(0:3)
         integer(8) :: icfl, jcfl, kcfl
         logical :: nonfin
         integer(8) :: i, j, k
         real(8) :: cflt
         integer(8) :: it, jt, kt

         cflm = 0.
         icfl = 1
         jcfl = 1
         kcfl = 1
         nonfin = .false.

!     EACH THREAD FINDS THE MAXIMUM OF ITS CELLS; THE MERGE KEEPS THE SERIAL PICK
!     (LARGEST VALUE, AND THE LAST CELL IN K-J-I ORDER AMONG EQUAL VALUES).
!$OMP PARALLEL PRIVATE(CFLI, CFLT, IT, JT, KT)
         cflt = 0.
         it = 1
         jt = 1
         kt = 1
!$OMP DO COLLAPSE(2) REDUCTION(.OR.:NONFIN)
         do k = 1, n3m
           do j = 1, n2m
             do i = 1, n1m
               cfli(1) = abs(u(i, j, k) + u(i + 1, j, k)) * 0.5 * f2fxi(i)
               cfli(2) = abs(v(i, j, k) + v(i, j + 1, k)) * 0.5 * f2fyi(j)
               cfli(3) = abs(w(i, j, k) + w(i, j, k + 1)) * 0.5 * f2fzi(k)
               cfli(0) = (cfli(1) + cfli(2) + cfli(3)) * dt
               if (.not. (cfli(0) .le. huge(1.0d0))) nonfin = .true.
               if (cfli(0) .ge. cflt) then
                 cflt = cfli(0)
                 it = i
                 jt = j
                 kt = k
               end if
             end do
           end do
         end do
!$OMP END DO
!$OMP CRITICAL(CFL_MAX)
         if (cflt .gt. cflm .or. (cflt .eq. cflm .and. ((kt - kcfl) * n2m + (jt - jcfl)) * n1m + (it - icfl) .gt. 0)) then
           cflm = cflt
           icfl = it
           jcfl = jt
           kcfl = kt
         end if
!$OMP END CRITICAL(CFL_MAX)
!$OMP END PARALLEL

         cflmax = cflm

         if (ntime .ne. 1) then
           write (*, 149) cflmax, xmp(icfl), ymp(jcfl), zmp(kcfl), icfl, jcfl, kcfl
149        format('CFLMAX =  ', f10.7, ' @ ', 3f10.4, ' , ', 3i5)
         end if

         return
       end subroutine cfl
!=======================================================================
!=======================================================================
       subroutine substepinit(substep)
!=======================================================================
         use mod_common
         use mod_flowarray
         implicit none
         integer(8) :: substep        ! RK3 SUBSTEP = 1, 2, 3

         alpha = 0.5 * (gamma(substep) + ro(substep))
         dtconst = dt * (gamma(substep) + ro(substep))
         dtconsti = 1./dtconst
         test1 = resid1 * dtconsti               ! POISS. CONVG. CRITERION
         acoef = alpha * dt / re
         acoefi = 1./acoef
         pmiavg = pmiavg + pmi(0)*2.*ALPHA
         subdt = subdt + dtconst               ! SUBTIME FOR RK3 METHOD
         msub = substep

         return
       end subroutine substepinit
!=======================================================================
!=======================================================================
       subroutine qvolcalc(qq)
!=======================================================================
         use mod_common
         use mod_flowarray, only: u
         implicit none
         real(8) :: qq, funcbody
         real(8) :: qq_k, qq_pln(n3m)
         integer(8) :: i, j, k

         call body_values()
         qq = 0.

!$OMP PARALLEL DO PRIVATE(QQ_K)
         do k = 1, n3m
           qq_k = 0.0d0
           do j = 1, n2m
             do i = 1, n1m
               if (fbody_u(i, j, k) .ge. 1.e-10) then
                 qq_k = qq_k + u(i, j, k) * c2cx(i) * f2fy(j) * f2fz(k)
               end if
             end do
           end do
           qq_pln(k) = qq_k
         end do
!$OMP END PARALLEL DO
         do k = 1, n3m
           qq = qq + qq_pln(k)
         end do

         return
       end subroutine qvolcalc
!=======================================================================
!=======================================================================
       subroutine tvolcalc(qq)
!=======================================================================
         use mod_common
         use mod_flowarray, only: t, cstar
         implicit none
         real(8) :: qq, funcbody
         real(8) :: qq_k, qq_pln(n3m)
         integer(8) :: i, j, k

         qq = 0.

!$OMP PARALLEL DO PRIVATE(QQ_K)
         do k = 1, n3m
           qq_k = 0.0d0
           do j = 1, n2m
             do i = 1, n1m
               qq_k = qq_k + t(i, j, k) / cstar(i, j, k) * f2fx(i) * f2fy(j) * f2fz(k)
             end do
           end do
           qq_pln(k) = qq_k
         end do
!$OMP END PARALLEL DO
         do k = 1, n3m
           qq = qq + qq_pln(k)
         end do

         return
       end subroutine tvolcalc
!=======================================================================
!=======================================================================
       subroutine qvolcorr
!=======================================================================
         use mod_common
         use mod_flowarray
         implicit none
         real(8) :: funcbody
         real(8) :: flowvol
         integer(8) :: i, j, k

         flowvol = 0.

! !$OMP PARALLEL DO REDUCTION(+:FLOWVOL)
!          do k = 1, n3m
!            do j = 1, n2m
!              do i = 1, n1m
!                if (funcbody(x(i), ymp(j), zmp(k)) .ge. 1.e-10) then
!                  flowvol = flowvol + c2cx(i) * f2fy(j) * f2fz(k)
!                end if
!              end do
!            end do
!          end do
! !$OMP END PARALLEL DO
         ! phcap = (qvol(2) - qvol(1)) / flowvol

         phcap = 0.D0 ! DISABLE THE ARTIFICIAL FLOW RATE CORRECTION
                      ! INSTEAD, WE IMPLEMENT THE CONTROLLER IN MEANPG
         return
       end subroutine qvolcorr
!=======================================================================
!=======================================================================
       subroutine tvolcorr
!=======================================================================
         use mod_common
         use mod_flowarray
         implicit none
         real(8) :: cvol
         real(8) :: cvol_k, cvol_pln(n3m)
         integer(8) :: i, j, k

         thcap = 0.
         cvol = 0.

!$OMP PARALLEL DO PRIVATE(CVOL_K)
         do k = 1, n3m
           cvol_k = 0.0d0
           do j = 1, n2m
             do i = 1, n1m
               cvol_k = cvol_k + 1./cstar(i, j, k) * f2fx(i) * f2fy(j) * f2fz(k)
             end do
           end do
           cvol_pln(k) = cvol_k
         end do
!$OMP END PARALLEL DO
         do k = 1, n3m
           cvol = cvol + cvol_pln(k)
         end do

         thcap = (tvol(2) - tvol(1)) / cvol

         return
       end subroutine tvolcorr
!=======================================================================
!=======================================================================
       subroutine ucalc(phi)
!=======================================================================
!
!     CALCULATE VELOCITY (U_I) FROM U_I HAT
!     U_I HAT IS DERIVED FROM PSEUDO-PRESSURE, PHI.
!
!     DEFINITION OF PSEUDO-PRESSURE, PHI, IS AS FOLLOWS:
!     ({U_I}^K - U_I HAT)/(2.*ALPHA_K*DT) = - D({PHI}^K)/D(X_I)
!
!     FROM THE DEFINITION OF PHI, FOLLOWING EQUATION IS DERIVED:
!     {U_I}^K = U_I HAT - 2.*A_K*DT*D({PHI}^K)/D(X_I)
!
!-----------------------------------------------------------------------
         use mod_common
         use mod_flowarray, only: u, v, w, t
         implicit none
         integer(8) :: i, j, k
         real(8) :: idum, funcbody
         integer(8) :: im, km
         real(8) :: phi(0:n1, 0:n2, 0:n3), flxcr

!$OMP PARALLEL DO COLLAPSE(2) &
!$OMP PRIVATE(IM)
         do k = 1, n3m
           do j = 1, n2m
             do i = i_bgpx, n1m                  ! I=1,N1 => BOUNDARY
               im = imv(i)
               u(i, j, k) = u(i, j, k) &
                            - dtconst * (phi(i, j, k) - phi(im, j, k)) * c2cxi(i)
               ! DISABLED: PHCAP IS NOW 0.0 (SEE QCORR), AND WE USE THE ACTIVE CONTROLLER IN MEANPG  
               ! if (ich .eq. 1) then
               !   if (funcbody(x(i), ymp(j), zmp(k)) .ge. 1.e-10) then
               !     u(i, j, k) = u(i, j, k) - phcap
               !   end if
               ! end if
             end do
           end do
         end do

!$OMP PARALLEL DO COLLAPSE(2)
         do k = 1, n3m
           do j = 2, n2m                       ! J=1,N2 => BOUNDARY
             do i = 1, n1m
               v(i, j, k) = v(i, j, k) - dtconst * (phi(i, j, k) - phi(i, j - 1, k)) * c2cyi(j)
             end do
           end do
         end do

!$OMP PARALLEL DO COLLAPSE(2) &
!$OMP PRIVATE(KM)
         do k = k_bgpz, n3m                  ! K=0,N3 => BOUNDARY
           do j = 1, n2m
             km = kmv(k)
             do i = 1, n1m
               w(i, j, k) = w(i, j, k) - dtconst * (phi(i, j, k) - phi(i, j, km)) * c2czi(k)
             end do
           end do
         end do

         idum = ihist
         flxcr = 0.

!!!!!!     NEUMANN & DIRICHLET BOUNDARY CONDITION
         if (jut .eq. 1) then
!$OMP PARALLEL DO
           do k = 0, n3
             do i = 1, n1
               u(i, n2, k) = u(i, n2m, k)
             end do
           end do
         end if
         if (jwt .eq. 1) then
!$OMP PARALLEL DO
           do k = 1, n3
             do i = 0, n1
               w(i, n2, k) = w(i, n2m, k)
             end do
           end do
         end if
         if (jub .eq. 1) then
!$OMP PARALLEL DO
           do k = 0, n3
             do i = 1, n1
               u(i, 0, k) = u(i, 1, k)
             end do
           end do
         end if
         if (jwb .eq. 1) then
!$OMP PARALLEL DO
           do k = 1, n3
             do i = 0, n1
               w(i, 0, k) = w(i, 1, k)
             end do
           end do
         end if
         if (kut .eq. 1) then
!$OMP PARALLEL DO
           do j = 0, n2
             do i = 1, n1
               u(i, j, n3) = u(i, j, n3m)
             end do
           end do
         end if
         if (kvt .eq. 1) then
!$OMP PARALLEL DO
           do j = 1, n2
             do i = 0, n1
               v(i, j, n3) = v(i, j, n3m)
             end do
           end do
         end if
         if (kub .eq. 1) then
!$OMP PARALLEL DO
           do j = 0, n2
             do i = 1, n1
               u(i, j, 0) = u(i, j, 1)
             end do
           end do
         end if
         if (kvb .eq. 1) then
!$OMP PARALLEL DO
           do j = 1, n2
             do i = 0, n1
               v(i, j, 0) = v(i, j, 1)
             end do
           end do
         end if

!     Z PERIODICITY
         if (zprdic .eq. 1) then
!$OMP PARALLEL DO
           do j = 0, n2
             do i = 1, n1
               u(i, j, 0) = u(i, j, n3m)
               u(i, j, n3) = u(i, j, 1)
             end do
           end do

!$OMP PARALLEL DO
           do j = 1, n2
             do i = 0, n1
               v(i, j, 0) = v(i, j, n3m)
               v(i, j, n3) = v(i, j, 1)
             end do
           end do

!$OMP PARALLEL DO
           do j = 0, n2
             do i = 0, n1
               w(i, j, 0) = w(i, j, n3m)
               w(i, j, n3) = w(i, j, 1)
             end do
           end do
         end if

!     X PERIODICITY
         if (xprdic .eq. 1) then
!$OMP PARALLEL DO
           do k = 0, n3
             do j = 0, n2
               u(0, j, k) = u(n1m, j, k)
               u(n1, j, k) = u(1, j, k)
             end do
           end do

!$OMP PARALLEL DO
           do k = 0, n3
             do j = 1, n2
               v(0, j, k) = v(n1m, j, k)
               v(n1, j, k) = v(1, j, k)
             end do
           end do

!$OMP PARALLEL DO
           do k = 1, n3
             do j = 0, n2
               w(0, j, k) = w(n1m, j, k)
               w(n1, j, k) = w(1, j, k)
             end do
           end do
         end if

         return
       end

!=======================================================================
       subroutine pcalc(phi, divgsum)
!=======================================================================
!
!     CALCULATE PRESSURE FROM DIVGSUM & PHI
!
!     ({U_I}^K-U_I HAT)/(2.*A_K*DT)
!                             = - [ D(P^K)/D(X_I) - D(P^(K-1))/D(X_I) ]
!                               + 0.5*[L({U_I}^K) - L(U_I HAT)]
!                             = - D({PHI}^K)/D(X_I)
!
!     BY THE DEFINITION OF PHI, FINAL EQUATION IS DERIVED AS FOLLOWS:
!       P^K = P^(K-1) + {PHI}^K - (A_K*DT/RE)*(D2({PHI^K})/D(X_J)D(X_J))
!
!-----------------------------------------------------------------------
         use mod_common
         use mod_flowarray, only: u, v, w, p
         implicit none
         integer(8) :: i, j, k
         real(8) :: phiref
         real(8) :: phi(0:n1, 0:n2, 0:n3), divgsum(0:n1, 0:n2, 0:n3)
         real(8) :: funcbody

!$OMP PARALLEL DO COLLAPSE(2)
         do k = 1, n3m
           do j = 1, n2m
             do i = 1, n1m
               p(i, j, k) = p(i, j, k) + phi(i, j, k) - 0.5 * dtconst * divgsum(i, j, k) * 1./re
             end do
           end do
         end do

!     SET THE AVERAGE PHI AT THE UPPER WALL TO BE ZERO.
         phiref = 0.
         do k = 1, n3m
           do j = 1, n2m
             phiref = phi(n1m, j, k) * f2fy(j) * f2fz(k) + phiref
           end do
         end do
         phiref = phiref / (yl * zl)

!$OMP PARALLEL DO COLLAPSE(2)
         do k = 1, n3m
           do j = 1, n2m
             do i = 1, n1m
               p(i, j, k) = p(i, j, k) - phiref
             end do
           end do
         end do

!     Z PERIODICITY
         if (zprdic .eq. 1) then
!$OMP PARALLEL DO
           do j = 1, n2m
             do i = 1, n1m
               p(i, j, 0) = p(i, j, n3m)
               p(i, j, n3) = p(i, j, 1)
             end do
           end do
         end if

!     X PERIODICITY
         if (xprdic .eq. 1) then
!$OMP PARALLEL DO
           do k = 1, n3m
             do j = 1, n2m
               p(0, j, k) = p(n1m, j, k)
               p(n1, j, k) = p(1, j, k)
             end do
           end do
         end if

         return
       end

!=======================================================================
       subroutine tcalc
!=======================================================================
         use mod_common
         use mod_flowarray, only: t, cstar
         implicit none
         integer(8) :: i, j, k
         real(8) :: funcbody, hflux, ttemp

         if ((ich .eq. 1) .and. (xprdic .eq. 0)) then
!$OMP PARALLEL DO COLLAPSE(2)
           do k = 1, n3m
             do j = 1, n2m
               do i = 1, n1m
                 t(i, j, k) = t(i, j, k) - thcap
               end do
             end do
           end do
         end if

         return
       end

!=======================================================================
       subroutine lagforce
!=======================================================================
!
!     MATERIAL DERIVATE OF VELOCITY, ONE OF MAJOR CONTRIBUTION OF THE
!       FORCE ON A BODY
!
!     REFERENCE
!       LEE ET AL., 2011, SOURCES OF SPURIOUS FORCE OSCILLATIONS
!       FROM AN IMMERSED BOUNDARY METHOD FOR MOVING-BODY PROBLEMS,
!       J. COMP. PHYS., 230, 2677-2695.
!
!     VARIABLES
!       DUDTA,DVDTA,DWDTA: VOLUME INTEGRATION OF MATERIAL DERIVATIVE OF
!                          VELOCITY OVER A BODY
!       DUDTR: TIME DERIVATIVE OF VELOCITY
!       RK3XO,RK3XOO: CONVECTIVE TERMS AT EACH RK3 SUB-STEP
!
!-----------------------------------------------------------------------
         use mod_common
         use mod_flowarray, only: rk3xo, rk3yo, rk3zo, dudtr &
                                  , rk3xoo, rk3yoo, rk3zoo, ifc, jfc, kfc
         implicit none
         integer(8) :: i, j, k, n, l, ii, jj, kk, ib, nb
         integer(8), parameter :: lblk = 1024
         real(8) :: dsum, dsum_b(3, (max(nbody(1), nbody(2), nbody(3)) + lblk - 1) / lblk)

         if (msub .eq. 1) then
           dudta = 0.
           dvdta = 0.
           dwdta = 0.
         end if

!     FIXED BLOCKS OF LBLK POINTS, SO THE SUMS DO NOT DEPEND ON THE NUMBER OF THREADS
         nb = (max(nbody(1), nbody(2), nbody(3)) + lblk - 1) / lblk
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(II,JJ,KK,N,DSUM)
         do l = 1, 3
           do ib = 1, nb
             dsum = 0.
             do n = (ib - 1) * lblk + 1, min(ib * lblk, nbody(l))
               ii = ifc(n, l)
               jj = jfc(n, l)
               kk = kfc(n, l)
               if (l .eq. 1) then
                 dsum = dsum + c2cx(ii) * f2fy(jj) * f2fz(kk) * (dudtr(n, l) - &
                                                                 (gamma(msub) * rk3xo(ii, jj, kk) + ro(msub) * rk3xoo(ii, jj, kk)))
               elseif (l .eq. 2) then
                 dsum = dsum + f2fx(ii) * c2cy(jj) * f2fz(kk) * (dudtr(n, l) - &
                                                                 (gamma(msub) * rk3yo(ii, jj, kk) + ro(msub) * rk3yoo(ii, jj, kk)))
               elseif (l .eq. 3) then
                 dsum = dsum + f2fx(ii) * f2fy(jj) * c2cz(kk) * (dudtr(n, l) - &
                                                                 (gamma(msub) * rk3zo(ii, jj, kk) + ro(msub) * rk3zoo(ii, jj, kk)))
               end if
             end do
             dsum_b(l, ib) = dsum
           end do
         end do
!$OMP END PARALLEL DO
         do ib = 1, nb
           dudta = dudta + dsum_b(1, ib)
           dvdta = dvdta + dsum_b(2, ib)
           dwdta = dwdta + dsum_b(3, ib)
         end do

         return
       end subroutine lagforce
!=======================================================================
!=======================================================================
       subroutine heatflux
!=======================================================================
         use mod_common
         use mod_flowarray, only: t, u, v, w
         implicit none
         integer(8) :: i, j, k
         real(8) :: area, tface
         real(8) :: flux_xm, flux_xp, flux_ym, flux_yp, flux_zm, flux_zp
         real(8) :: heat_flux(3), heat_flux_avg(3)
         real(8) :: funcbody
         ! funcbody on the y and z boundary faces, stored once for a static geometry
         real(8), allocatable, save :: fby(:, :, :), fbz(:, :, :)

         call body_values()
         if (.not. allocated(fby) .or. imovingon .eq. 1) then
           if (.not. allocated(fby)) allocate (fby(1:n1m, 1:n3m, 2), fbz(1:n1m, 1:n2m, 2))
!$OMP PARALLEL DO COLLAPSE(2)
           do k = 1, n3m
             do i = 1, n1m
               fby(i, k, 1) = funcbody(xmp(i), y(1), zmp(k), time)
               fby(i, k, 2) = funcbody(xmp(i), y(n2), zmp(k), time)
             end do
           end do
!$OMP PARALLEL DO COLLAPSE(2)
           do j = 1, n2m
             do i = 1, n1m
               fbz(i, j, 1) = funcbody(xmp(i), ymp(j), z(1), time)
               fbz(i, j, 2) = funcbody(xmp(i), ymp(j), z(n3), time)
             end do
           end do
         end if

         flux_xm = 0.d0
         flux_xp = 0.d0
         flux_ym = 0.d0
         flux_yp = 0.d0
         flux_zm = 0.d0
         flux_zp = 0.d0

         ! Face sums stay serial, so the history is the same for any number of threads.
         do k = 1, n3m
           do j = 1, n2m
             area = f2fy(j) * f2fz(k)
             if (fbody_u(1, j, k) .ge. 1.d-10) then
               tface = 0.5d0 * (t(1, j, k) + t(0, j, k))
               flux_xm = flux_xm + u(1, j, k) * tface * area
             end if
             if (fbody_u(n1, j, k) .ge. 1.d-10) then
               tface = 0.5d0 * (t(n1m, j, k) + t(n1, j, k))
               flux_xp = flux_xp + u(n1, j, k) * tface * area
             end if
           end do
         end do

         do k = 1, n3m
           do i = 1, n1m
             area = f2fx(i) * f2fz(k)
             if (fby(i, k, 1) .ge. 1.d-10) then
               tface = 0.5d0 * (t(i, 1, k) + t(i, 0, k))
               flux_ym = flux_ym + v(i, 1, k) * tface * area
             end if
             if (fby(i, k, 2) .ge. 1.d-10) then
               tface = 0.5d0 * (t(i, n2m, k) + t(i, n2, k))
               flux_yp = flux_yp + v(i, n2, k) * tface * area
             end if
           end do
         end do

         do j = 1, n2m
           do i = 1, n1m
             area = f2fx(i) * f2fy(j)
             if (fbz(i, j, 1) .ge. 1.d-10) then
               tface = 0.5d0 * (t(i, j, 1) + t(i, j, 0))
               flux_zm = flux_zm + w(i, j, 1) * tface * area
             end if
             if (fbz(i, j, 2) .ge. 1.d-10) then
               tface = 0.5d0 * (t(i, j, n3m) + t(i, j, n3))
               flux_zp = flux_zp + w(i, j, n3) * tface * area
             end if
           end do
         end do

         heat_flux(1) = flux_xp - flux_xm
         heat_flux(2) = flux_yp - flux_ym
         heat_flux(3) = flux_zp - flux_zm

         heat_flux_avg = 0.0d0
         if (time .ge. avg_tst) then
           hfluxavg_int(1) = hfluxavg_int(1) + heat_flux(1) * dt
           hfluxavg_int(2) = hfluxavg_int(2) + heat_flux(2) * dt
           hfluxavg_int(3) = hfluxavg_int(3) + heat_flux(3) * dt
           hfluxavg_dur = hfluxavg_dur + dt
           if (hfluxavg_dur .gt. 0.0d0) then
             heat_flux_avg(1) = hfluxavg_int(1) / hfluxavg_dur
             heat_flux_avg(2) = hfluxavg_int(2) / hfluxavg_dur
             heat_flux_avg(3) = hfluxavg_int(3) / hfluxavg_dur
           end if
         end if

         if (mod(ntime, npin) .eq. 0) then
           write (2002, 110) time, heat_flux(1), heat_flux(2), heat_flux(3), &
             heat_flux_avg(1), heat_flux_avg(2), heat_flux_avg(3)
         end if
110      format(f13.5, 6es15.7)

         return
       end subroutine heatflux
!=======================================================================
!=======================================================================
       subroutine draglift
!=======================================================================
!     CALCULATES THE TOTAL X, Y, Z AERODYNAMIC FORCES ACTING ON THE
!     IMMERSED BOUNDARY (IB) BODIES.
!
!     REQUIRES DUDTA, DVDTA, DWDTA TO BE PRE-CALCULATED BY THE
!     LAGFORCE SUBROUTINE DURING THE RK3 SUBSTEPS.
!-----------------------------------------------------------------------
         use mod_common
         use mod_flowarray
         implicit none
         integer(8) :: i, j, k, n, l, ii, jj, kk
         real(8) :: cd(3), cdavg(3), vol_solid_geom, vol_cell
         real(8) :: dti, funcbody, tmp, vol_fluid
         real(8) :: dy1, dy2, dz1, dz2, grad_u, grad_v, grad_w
         ! fixed-order partial sums (per plane, per block of LBLK points), the same for any thread count
         integer(8) :: ib, nb
         integer(8), parameter :: lblk = 1024
         real(8) :: vs_k, vs_pln(n3m), csum, csum_b(3, (max(nbody(1), nbody(2), nbody(3)) + lblk - 1) / lblk)

         dti = 1.0d0 / dt
         call body_values()

         cd = 0.0d0
         vol_solid_geom = 0.0d0

         ! 1. ESTIMATE SOLID VOLUME FROM THE CURRENT GEOMETRY MASK.
!$OMP PARALLEL DO PRIVATE(VS_K)
         do k = 1, n3m
           vs_k = 0.0d0
           do j = 1, n2m
             do i = 1, n1m
               if (fbody_u(i, j, k) .lt. 1.d-10) then
                 vs_k = vs_k + c2cx(i) * f2fy(j) * f2fz(k)
               end if
             end do
           end do
           vs_pln(k) = vs_k
         end do
!$OMP END PARALLEL DO
         do k = 1, n3m
           vol_solid_geom = vol_solid_geom + vs_pln(k)
         end do

         vol_fluid = xl * yl * zl - vol_solid_geom
         if (vol_fluid .le. 0.0d0) vol_fluid = xl * yl * zl

         ! 2. INTEGRATE THE RAW IBM FORCING (FCVAVG) OVER THE ACTIVE IBM NODES
         nb = (max(nbody(1), nbody(2), nbody(3)) + lblk - 1) / lblk
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(N, II, JJ, KK, VOL_CELL, CSUM)
         do l = 1, 3
           do ib = 1, nb
             csum = 0.0d0
             do n = (ib - 1) * lblk + 1, min(ib * lblk, nbody(l))
               ii = ifc(n, l)
               jj = jfc(n, l)
               kk = kfc(n, l)

               vol_cell = c2cx(ii) * f2fy(jj) * f2fz(kk)

               ! RAW IBM FORCE (TERM 1 OF EQ. 11)
               csum = csum - fcvavg(n, l) * vol_cell
             end do
             csum_b(l, ib) = csum
           end do
         end do
!$OMP END PARALLEL DO
         do ib = 1, nb
           do l = 1, 3
             cd(l) = cd(l) + csum_b(l, ib)
           end do
         end do

         ! 3. ADD THE MATERIAL DERIVATIVE (TERM 2 OF EQ. 11)
         !    THESE ARE DIRECTLY SUPPLIED BY YOUR EXISTING LAGFORCE SUBROUTINE.
         !    FOR STATIONARY BODIES, DUDT ARE EFFECTIVELY ZERO, SO THIS MIGHT HAVE A NEGLIGIBLE CONTRIBUTION.
         cd(1) = cd(1) + dudta
         cd(2) = cd(2) + dvdta
         cd(3) = cd(3) + dwdta

         ! 4. ADD COMPUTATIONAL-WALL SHEAR CONTRIBUTIONS
         ! --- Y-BOTTOM WALL ---
         if (bc_ybtm .eq. 0) then
           tmp = 0.0d0
           dy1 = 1.0d0 / c2cyi(1)
           dy2 = dy1 + 1.0d0 / c2cyi(2)
           do k = 1, n3m
             do i = 1, n1m
               if (fbody_u(i, 1, k) .ge. 1.d-10) then
                 grad_u = (u(i, 1, k) * dy2**2 - u(i, 2, k) * dy1**2) / (dy1 * dy2 * (dy2 - dy1))
                 tmp = grad_u * c2cx(i) * f2fz(k)
                 cd(1) = cd(1) + tmp / re

                 grad_w = (w(i, 1, k) * dy2**2 - w(i, 2, k) * dy1**2) / (dy1 * dy2 * (dy2 - dy1))
                 tmp = grad_w * c2cx(i) * f2fz(k)
                 cd(3) = cd(3) + tmp / re
               end if
             end do
           end do
         end if

         ! --- Y-TOP WALL ---
         if (bc_ytop .eq. 0) then
           tmp = 0.0d0
           dy1 = 1.0d0 / c2cyi(n2)
           dy2 = dy1 + 1.0d0 / c2cyi(n2m)
           do k = 1, n3m
             do i = 1, n1m
               if (fbody_u(i, n2m, k) .ge. 1.d-10) then
                 grad_u = (u(i, n2m, k) * dy2**2 - u(i, n2m-1, k) * dy1**2) / (dy1 * dy2 * (dy2 - dy1))
                 tmp = grad_u * c2cx(i) * f2fz(k)
                 cd(1) = cd(1) + tmp / re

                 grad_w = (w(i, n2m, k) * dy2**2 - w(i, n2m-1, k) * dy1**2) / (dy1 * dy2 * (dy2 - dy1))
                 tmp = grad_w * c2cx(i) * f2fz(k)
                 cd(3) = cd(3) + tmp / re
               end if
             end do
           end do
         end if

         ! --- Z-BOTTOM WALL ---
         if (bc_zbtm .eq. 0) then
           tmp = 0.0d0
           dz1 = 1.0d0 / c2czi(1)
           dz2 = dz1 + 1.0d0 / c2czi(2)
           do j = 1, n2m
             do i = 1, n1m
               if (fbody_u(i, j, 1) .ge. 1.d-10) then
                 grad_u = (u(i, j, 1) * dz2**2 - u(i, j, 2) * dz1**2) / (dz1 * dz2 * (dz2 - dz1))
                 tmp = grad_u * c2cx(i) * f2fy(j)
                 cd(1) = cd(1) + tmp / re

                 grad_v = (v(i, j, 1) * dz2**2 - v(i, j, 2) * dz1**2) / (dz1 * dz2 * (dz2 - dz1))
                 tmp = grad_v * c2cx(i) * f2fy(j)
                 cd(2) = cd(2) + tmp / re
               end if
             end do
           end do
         end if

         ! --- Z-TOP WALL ---
         if (bc_ztop .eq. 0) then
           tmp = 0.0d0
           dz1 = 1.0d0 / c2czi(n3)
           dz2 = dz1 + 1.0d0 / c2czi(n3m)
           do j = 1, n2m
             do i = 1, n1m
               if (fbody_u(i, j, n3m) .ge. 1.d-10) then
                 grad_u = (u(i, j, n3m) * dz2**2 - u(i, j, n3m-1) * dz1**2) / (dz1 * dz2 * (dz2 - dz1))
                 tmp = grad_u * c2cx(i) * f2fy(j)
                 cd(1) = cd(1) + tmp / re

                 grad_v = (v(i, j, n3m) * dz2**2 - v(i, j, n3m-1) * dz1**2) / (dz1 * dz2 * (dz2 - dz1))
                 tmp = grad_v * c2cx(i) * f2fy(j)
                 cd(2) = cd(2) + tmp / re
               end if
             end do
           end do
         end if

         ! 5. NON-DIMENSIONALIZE TOTAL FORCE TO EVALUATE WALL SHEAR STRESS (TAU_W)
         cd(1) = cd(1) / (xl * yl * zl)
         cd(2) = cd(2) / (xl * yl * zl)
         cd(3) = cd(3) / (xl * yl * zl)

         cdavg = 0.0d0
         if (time .ge. avg_tst) then
           cdavg_int(1) = cdavg_int(1) + cd(1) * dt
           cdavg_int(2) = cdavg_int(2) + cd(2) * dt
           cdavg_int(3) = cdavg_int(3) + cd(3) * dt
           cdavg_dur = cdavg_dur + dt
           if (cdavg_dur .gt. 0.0d0) then
             cdavg(1) = cdavg_int(1) / cdavg_dur
             cdavg(2) = cdavg_int(2) / cdavg_dur
             cdavg(3) = cdavg_int(3) / cdavg_dur
           end if
         end if

         ! 6. OUTPUT TO HISTORY FILE
         if (mod(ntime, npin) .eq. 0) then
           write (2001, 110) time, cd(1), cd(2), cd(3), cdavg(1), cdavg(2), cdavg(3)
         end if
110      format(f13.5, 6es15.7)

         return
       end subroutine draglift
!=======================================================================
