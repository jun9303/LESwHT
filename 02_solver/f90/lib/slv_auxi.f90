!=======================================================================
       subroutine stepinit
!=======================================================================
         use mod_common
         use mod_flowarray
         implicit none
         real(8) :: dtcfl
         integer(8) :: icfl, jcfl, kcfl

         ntime = ntime + 1
         m = ntime
         fcvavg = 0.
         if (ich .eq. 1) pmiavg = 0.

         call cfl(cflmax, icfl, jcfl, kcfl)         ! CALCULATE CFL NUMBER

         if (idtopt .ne. 0 .and. cflmax .ne. 0.) then
           dtcfl = dmin1(dt * cflfac / cflmax, dt * (0.80 + 0.20 * cflfac / cflmax))
           if (idtopt .eq. 1) dt = dtcfl
         end if

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
       subroutine cfl(cflm, icfl, jcfl, kcfl)
!=======================================================================
!
!     CALCULATE THE MAXIMUM CFL NUMBER OF FLOW FIELD
!     (HTTP://EN.WIKIPEDIA.ORG/WIKI/COURANT-FRIEDRICHS-LEWY_CONDITION)
!
!     CFL#=U_I*DT/DX_I
!     CFLI=(UC/DX+VC/DY+WC/DZ)*DT  (C FOR CELL CENTER)
!     CFLMPT: INDEX WHERE THE CFL# IS MAXIMUM
!
!-----------------------------------------------------------------------
         use mod_common
         use mod_flowarray, only: u, v, w
         implicit none
         real(8) :: cflm
         real(8) :: cfli(0:3)
         integer(8) :: icfl, jcfl, kcfl
         integer(8) :: i, j, k

         cflm = 0.

!$OMP PARALLEL DO PRIVATE(CFLI) REDUCTION(MAX:CFLM)
         do k = 1, n3m
           do j = 1, n2m
             do i = 1, n1m
               cfli(1) = abs(u(i, j, k) + u(i + 1, j, k)) * 0.5 * f2fxi(i)
               cfli(2) = abs(v(i, j, k) + v(i, j + 1, k)) * 0.5 * f2fyi(j)
               cfli(3) = abs(w(i, j, k) + w(i, j, k + 1)) * 0.5 * f2fzi(k)
               cfli(0) = (cfli(1) + cfli(2) + cfli(3)) * dt
               if (cfli(0) .ge. cflm) then
                 cflm = cfli(0)
                 icfl = i
                 jcfl = j
                 kcfl = k
               end if
             end do
           end do
         end do
!$OMP END PARALLEL DO

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
         pmiavg = pmiavg + pmi(0) !*2.*ALPHA
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
         integer(8) :: i, j, k

         qq = 0.

!$OMP PARALLEL DO REDUCTION(+:QQ)
         do k = 1, n3m
           do j = 1, n2m
             do i = 1, n1m
               if (funcbody(x(i), ymp(j), zmp(k)) .ge. 1.e-10) then
                 qq = qq + u(i, j, k) * c2cx(i) * f2fy(j) * f2fz(k)
               end if
             end do
           end do
         end do
!$OMP END PARALLEL DO

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
         integer(8) :: i, j, k

         qq = 0.

!$OMP PARALLEL DO REDUCTION(+:QQ)
         do k = 1, n3m
           do j = 1, n2m
             do i = 1, n1m
               qq = qq + t(i, j, k) / cstar(i, j, k) * f2fx(i) * f2fy(j) * f2fz(k)
             end do
           end do
         end do
!$OMP END PARALLEL DO

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

!$OMP PARALLEL DO REDUCTION(+:FLOWVOL)
         do k = 1, n3m
           do j = 1, n2m
             do i = 1, n1m
               if (funcbody(x(i), ymp(j), zmp(k)) .ge. 1.e-10) then
                 flowvol = flowvol + c2cx(i) * f2fy(j) * f2fz(k)
               end if
             end do
           end do
         end do
!$OMP END PARALLEL DO

         ! PHCAP = (QVOL(2) - QVOL(1)) / FLOWVOL
         phcap = (qvol(2) - ubulk_i * flowvol) / flowvol

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
         integer(8) :: i, j, k

         thcap = 0.
         cvol = 0.

!$OMP PARALLEL DO REDUCTION(+:CVOL)
         do k = 1, n3m
           do j = 1, n2m
             do i = 1, n1m
               cvol = cvol + 1./cstar(i, j, k) * f2fx(i) * f2fy(j) * f2fz(k)
             end do
           end do
         end do
!$OMP END PARALLEL DO

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

!$OMP PARALLEL DO &
!$OMP PRIVATE(IM)
         do k = 1, n3m
           do j = 1, n2m
             do i = i_bgpx, n1m                  ! I=1,N1 => BOUNDARY
               im = imv(i)
               u(i, j, k) = u(i, j, k) &
                            - dtconst * (phi(i, j, k) - phi(im, j, k)) * c2cxi(i)
               if (ich .eq. 1) then
                 if (funcbody(x(i), ymp(j), zmp(k)) .ge. 1.e-10) then
                   u(i, j, k) = u(i, j, k) - phcap
                 end if
               end if
             end do
           end do
         end do

!$OMP PARALLEL DO
         do k = 1, n3m
           do j = 2, n2m                       ! J=1,N2 => BOUNDARY
             do i = 1, n1m
               v(i, j, k) = v(i, j, k) - dtconst * (phi(i, j, k) - phi(i, j - 1, k)) * c2cyi(j)
             end do
           end do
         end do

!$OMP PARALLEL DO &
!$OMP PRIVATE(KM)
         do k = k_bgpz, n3m                  ! K=0,N3 => BOUNDARY
           km = kmv(k)
           do j = 1, n2m
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

!$OMP PARALLEL DO
         do k = 1, n3m
           do j = 1, n2m
             do i = 1, n1m
               p(i, j, k) = p(i, j, k) + phi(i, j, k) - 0.5 * dtconst * divgsum(i, j, k) * 1./re
             end do
           end do
         end do

!     SET THE AVERAGE PHI AT THE UPPER WALL TO BE ZERO.
         phiref = 0.
!$OMP PARALLEL DO &
!$OMP REDUCTION(+:PHIREF)
         do k = 1, n3m
           do j = 1, n2m
             phiref = phi(n1m, j, k) * f2fy(j) * f2fz(k) + phiref
           end do
         end do
         phiref = phiref / (yl * zl)

!$OMP PARALLEL DO
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
!$OMP PARALLEL DO
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
         integer(8) :: i, j, k, n, l, ii, jj, kk

         if (msub .eq. 1) then
           dudta = 0.
           dvdta = 0.
           dwdta = 0.
         end if

!$OMP PARALLEL DO PRIVATE(II,JJ,KK)
         do l = 1, 3
           do n = 1, nbody(l)
             ii = ifc(n, l)
             jj = jfc(n, l)
             kk = kfc(n, l)
             !INTERMEDIATE INFORMATION
             if (l .eq. 1) then
               dudta = dudta + c2cx(ii) * f2fy(jj) * f2fz(kk) * (dudtr(n, l) - &
                                                                 (gamma(msub) * rk3xo(ii, jj, kk) + ro(msub) * rk3xoo(ii, jj, kk)))
             elseif (l .eq. 2) then
               dvdta = dvdta + f2fx(ii) * c2cy(jj) * f2fz(kk) * (dudtr(n, l) - &
                                                                 (gamma(msub) * rk3yo(ii, jj, kk) + ro(msub) * rk3yoo(ii, jj, kk)))
             elseif (l .eq. 3) then
               dwdta = dwdta + f2fx(ii) * f2fy(jj) * c2cz(kk) * (dudtr(n, l) - &
                                                                 (gamma(msub) * rk3zo(ii, jj, kk) + ro(msub) * rk3zoo(ii, jj, kk)))
             end if

           end do
         end do
!$OMP END PARALLEL DO

         return
       end subroutine lagforce
!=======================================================================
!=======================================================================
       subroutine calc_boundary_heat_flux
!=======================================================================
         use mod_common
         use mod_flowarray, only: t, kstar
         implicit none
         integer(8) :: i, k
         real(8) :: dtdy_btm, dtdy_top
         real(8) :: q_btm, q_top, q_total
         real(8) :: area, k_solid_btm, k_solid_top

         q_btm = 0.d0
         q_top = 0.d0

!$OMP PARALLEL DO PRIVATE(DTDY_BTM, DTDY_TOP, AREA, K_SOLID_BTM, K_SOLID_TOP) REDUCTION(+:Q_BTM, Q_TOP)
         do k = 1, n3m
           do i = 1, n1m
             area = f2fx(i) * f2fz(k)

             ! EXTRACT THE SOLID THERMAL CONDUCTIVITY AT THE BOUNDARIES
             k_solid_btm = kstar(i, 1, k, 4)   ! SOUTH FACE OF J=1
             k_solid_top = kstar(i, n2m, k, 3) ! NORTH FACE OF J=N2M

             ! --- BOTTOM WALL (J=0 TO J=1) ---
             ! GRADIENT: DT/DY AT THE BOTTOM WALL
             dtdy_btm = (t(i, 1, k) - t(i, 0, k)) * c2cyi(1)
             ! HEAT FLOWS UP (+Y) INTO THE DOMAIN. FOURIER'S LAW: Q = -K * DT/DY
             q_btm = q_btm - (k_solid_btm / (re * pr)) * dtdy_btm * area

             ! --- TOP WALL (J=N2M TO J=N2) ---
             ! GRADIENT: DT/DY AT THE TOP WALL
             dtdy_top = (t(i, n2, k) - t(i, n2m, k)) * c2cyi(n2)
             ! HEAT FLOWS DOWN (-Y) INTO THE DOMAIN, SO THE SIGN IS FLIPPED.
             q_top = q_top + (k_solid_top / (re * pr)) * dtdy_top * area

           end do
         end do
!$OMP END PARALLEL DO

         q_total = q_btm + q_top

         ! WRITE TO THE FNUSSELT.DAT TRACKER
         write (2002, 110) time, q_btm, q_top, q_total
110      format(f12.5, 3es15.6)

         return
       end subroutine calc_boundary_heat_flux
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
         real(8) :: cd(3), vol_solid_geom, vol_cell
         real(8) :: dti, funcbody

         dti = 1.0d0 / dt

         cd = 0.0d0
         vol_solid_geom = 0.0d0

         ! 1. CALCULATE THE EXACT GEOMETRIC SOLID VOLUME USING FUNCBODY.
         !    THIS PERFECTLY ISOLATES THE FLUID VOLUME FROM THE SOLID SLABS
         !    SO THE GLOBAL PMIAVG PENALTY CAN BE NEUTRALIZED.
!$OMP PARALLEL DO REDUCTION(+:VOL_SOLID_GEOM)
         do k = 1, n3m
           do j = 1, n2m
             do i = 1, n1m
               if (funcbody(x(i), ymp(j), zmp(k), time) .lt. 1.d-10) then
                 vol_solid_geom = vol_solid_geom + c2cx(i) * f2fy(j) * f2fz(k)
               end if
             end do
           end do
         end do
!$OMP END PARALLEL DO

         ! 2. INTEGRATE THE RAW IBM FORCING (FCVAVG) OVER THE ACTIVE IBM NODES
!$OMP PARALLEL DO PRIVATE(N, L, II, JJ, KK, VOL_CELL) REDUCTION(-:CD)
         do l = 1, 3
           do n = 1, nbody(l)
             ii = ifc(n, l)
             jj = jfc(n, l)
             kk = kfc(n, l)

             vol_cell = c2cx(ii) * f2fy(jj) * f2fz(kk)

             ! RAW IBM FORCE (TERM 1 OF EQ. 11)
             cd(l) = cd(l) - fcvavg(n, l) * vol_cell
           end do
         end do
!$OMP END PARALLEL DO

         ! 3. ADD THE MATERIAL DERIVATIVE (TERM 2 OF EQ. 11)
         !    THESE ARE DIRECTLY SUPPLIED BY YOUR EXISTING LAGFORCE SUBROUTINE.
         !    FOR STATIONARY BODIES, DUDT ARE EFFECTIVELY ZERO, SO THIS MIGHT HAVE A NEGLIGIBLE CONTRIBUTION.
         cd(1) = cd(1) + dudta
         cd(2) = cd(2) + dvdta
         cd(3) = cd(3) + dwdta

         ! 4. NON-DIMENSIONALIZE TOTAL FORCE TO EVALUATE WALL SHEAR STRESS (TAU_W)
         !    DIVIDE BY THE TOTAL WETTED SURFACE AREA: 2 WALLS * (L_X * L_Z)
         cd(1) = cd(1) / (xl * zl * 2.d0)
         cd(2) = cd(2) / (xl * zl * 2.d0)
         cd(3) = cd(3) / (xl * zl * 2.d0)

         ! 5. OUTPUT TO HISTORY FILE
         if (mod(ntime, npin) .eq. 0) then
           write (2001, 110) time, cd(1), cd(2), cd(3)
         end if
110      format(f13.5, 3es15.6)

         return
       end subroutine draglift
!=======================================================================
