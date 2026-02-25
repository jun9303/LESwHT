!=======================================================================
       subroutine meanpg
!=======================================================================
!
!     THE MEAN PRESSURE GRADIENT (DP/DX) REQUIRED TO
!     KEEP THE MASS FLOW RATE CONSTANT IS DETERMINED
!     BY INTEGRATING THE WALL-SHEAR STRESSES AT THE CHANNEL WALLS.
!
!     WHEN THERE IS AN IBM BODY INSIDE THE CHANNEL,
!     THE IBM FORCING HAS TO BE INCLUDED.
!     -> FCV(N,1) IS CONSIDERED (FORCING IN THE STREAMWISE DIRECTION)
!
!     PMI(0): OVERALL MEAN PRESSURE GRADIENT
!     PMI(1): P.GRAD COMPONENT AT THE Y-BOTTOM WALL
!     PMI(2): P.GRAD COMPONENT AT THE Y-TOP    WALL
!     PMI(3): P.GRAD COMPONENT AT THE Z-BOTTOM WALL
!     PMI(4): P.GRAD COMPONENT AT THE Z-TOP    WALL
!
!-----------------------------------------------------------------------
!$       USE OMP_LIB
         use mod_common
         use mod_flowarray
         implicit none
         integer(8) :: i, j, k, n
         real(8) :: tmp, volume

         ! --- Controller Variables ---
         real(8) :: flowvol, qvol_current, u_bulk_current, err_u
         real(8), save :: err_integral = 0.0d0
         real(8) :: Kp, Ki
         real(8) :: funcbody
         real(8) :: dy1, dy2, dz1, dz2, grad_u

         flowvol = 0.0d0
         qvol_current = 0.0d0
         u_bulk_current = 0.0d0

!$OMP PARALLEL DO REDUCTION(+:FLOWVOL, QVOL_CURRENT)
         do k = 1, n3m
           do j = 1, n2m
             do i = 1, n1m
               if (funcbody(x(i), ymp(j), zmp(k), time) .ge. 1.e-10) then
                 flowvol = flowvol + c2cx(i) * f2fy(j) * f2fz(k)
                 qvol_current = qvol_current + u(i, j, k) * c2cx(i) * f2fy(j) * f2fz(k)
               end if
             end do
           end do
         end do
!$OMP END PARALLEL DO

         if (flowvol .gt. 0.0d0) then
           u_bulk_current = qvol_current / flowvol
         else
           u_bulk_current = 0.0d0
         end if

         pmi = 0.

         if (bc_ybtm .eq. 0) then
           tmp = 0.
           dy1 = 1.0d0 / c2cyi(1)
           dy2 = dy1 + 1.0d0 / c2cyi(2)
!$OMP PARALLEL DO PRIVATE(GRAD_U) REDUCTION(+:TMP)
           do k = 1, n3m
             do i = 1, n1m
               grad_u = (u(i, 1, k) * dy2**2 - u(i, 2, k) * dy1**2) / (dy1 * dy2 * (dy2 - dy1))
               tmp = tmp + grad_u * c2cx(i) * f2fz(k) / re
             end do
           end do
!$OMP END PARALLEL DO
           pmi(1) = pmi(1) + tmp

         end if

         if (bc_ytop .eq. 0) then
           tmp = 0.
           dy1 = 1.0d0 / c2cyi(n2)
           dy2 = dy1 + 1.0d0 / c2cyi(n2m)
!$OMP PARALLEL DO PRIVATE(GRAD_U) REDUCTION(+:TMP)
           do k = 1, n3m
             do i = 1, n1m
               grad_u = (u(i, n2m, k) * dy2**2 - u(i, n2m-1, k) * dy1**2) / (dy1 * dy2 * (dy2 - dy1))
               tmp = tmp + grad_u * c2cx(i) * f2fz(k) / re
             end do
           end do
!$OMP END PARALLEL DO
           pmi(2) = pmi(2) + tmp

         end if

         if (bc_zbtm .eq. 0) then
           tmp = 0.
           dz1 = 1.0d0 / c2czi(1)
           dz2 = dz1 + 1.0d0 / c2czi(2)
!$OMP PARALLEL DO PRIVATE(GRAD_U) REDUCTION(+:TMP)
           do j = 1, n2m
             do i = 1, n1m
               grad_u = (u(i, j, 1) * dz2**2 - u(i, j, 2) * dz1**2) / (dz1 * dz2 * (dz2 - dz1))
               tmp = tmp + grad_u * c2cx(i) * f2fy(j) / re
             end do
           end do
!$OMP END PARALLEL DO
           pmi(3) = pmi(3) + tmp

         end if

         if (bc_ztop .eq. 0) then
           tmp = 0.
           dz1 = 1.0d0 / c2czi(n3)
           dz2 = dz1 + 1.0d0 / c2czi(n3m)
!$OMP PARALLEL DO PRIVATE(GRAD_U) REDUCTION(+:TMP)
           do j = 1, n2m
             do i = 1, n1m
               grad_u = (u(i, j, n3m) * dz2**2 - u(i, j, n3m-1) * dz1**2) / (dz1 * dz2 * (dz2 - dz1))
               tmp = tmp + grad_u * c2cx(i) * f2fy(j) / re
             end do
           end do
!$OMP END PARALLEL DO
           pmi(4) = pmi(4) + tmp

         end if

         pmi(0) = pmi(1) + pmi(2) + pmi(3) + pmi(4)

         if ((ich .ne. 0) .and. (bc_ybtm .ne. 0) .and. (bc_ytop .ne. 0) .and. &
             (bc_zbtm .ne. 0) .and. (bc_ztop .ne. 0)) then
           write (*, *) ' TO USE ICH=/=0, AT LEAST ONE OF Y,Z-WALLS MUST'
           write (*, *) ' BE WALL SO THAT WALL-SHEAR STRESS EXISTS.'
           stop
         else
           tmp = 0.
!$OMP PARALLEL DO REDUCTION(+:TMP)
           do n = 1, nbody(1)
             tmp = tmp + fcv(n, 1) * c2cx(ifc(n, 1)) * f2fy(jfc(n, 1)) &
                   * f2fz(kfc(n, 1))
           end do
!$OMP END PARALLEL DO
           pmi(0) = pmi(0) + tmp
           volume = xl * yl * zl

           do i = 0, 4
             pmi(i) = pmi(i) / volume
           end do
         end if

         if (ich .eq. 0) then
           pmi(0) = 0.0d0
         elseif (ich .eq. 1) then
           ! =========================================================
           ! ACTIVE PI CONTROLLER FOR CONSTANT FLOW RATE (CFR)
           ! =========================================================
           if (ntime .gt. 0 .and. dt .gt. 0.0d0) then
             err_u = udrv_i - u_bulk_current

             Kp = 1.0d0
             Ki = 0.1d0

             err_integral = err_integral + err_u * dtconst

             pmi(0) = pmi(0) - (Kp * err_u + Ki * err_integral)
           end if
         elseif (ich .eq. 2) then
           ! CONSTANT PRESSURE GRADIENT (CPG)
           ! RHS USES "-pmi(0)" IN X-MOMENTUM, SO NEGATIVE PMI DRIVES +X FLOW.
           pmi(0) = -udrv_i
         end if

         return
       end subroutine meanpg
!=======================================================================
!=======================================================================
       subroutine lhsinit
!=======================================================================
         use mod_common
         implicit none
         integer(8) :: ic, jc, kc
         integer(8) :: im, ip, jm, jp, km, kp

         if (xprdic .eq. 1) then
           i_bgpx = 1
         else
           i_bgpx = 2
         end if

         if (yprdic .eq. 1) then
           j_bgpy = 1
         else
           j_bgpy = 2
         end if

         if (zprdic .eq. 1) then
           k_bgpz = 1
         else
           k_bgpz = 2
         end if

!$OMP PARALLEL DO
         do ic = i_bgpx, n1m
           im = imv(ic)
           aiu(ic) = -c2cxi(ic) * f2fxi(im)
           ciu(ic) = -c2cxi(ic) * f2fxi(ic)
         end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
         do ic = 1, n1m
           ip = ipv(ic)
           aivw(ic) = -c2cxi(ic) * f2fxi(ic)
           civw(ic) = -c2cxi(ip) * f2fxi(ic)
         end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
         do jc = j_bgpy, n2m
           jm = jmv(jc)
           ajv(jc) = -c2cyi(jc) * f2fyi(jm)
           cjv(jc) = -c2cyi(jc) * f2fyi(jc)
         end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
         do jc = 1, n2m
           jp = jpv(jc)
           ajuw(jc) = -c2cyi(jc) * f2fyi(jc)
           cjuw(jc) = -c2cyi(jp) * f2fyi(jc)
         end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
         do kc = k_bgpz, n3m
           km = kmv(kc)
           akw(kc) = -c2czi(kc) * f2fzi(km)
           ckw(kc) = -c2czi(kc) * f2fzi(kc)
         end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
         do kc = 1, n3m
           kp = kpv(kc)
           akuv(kc) = -c2czi(kc) * f2fzi(kc)
           ckuv(kc) = -c2czi(kp) * f2fzi(kc)
         end do
!$OMP END PARALLEL DO

         return
       end subroutine lhsinit
!=======================================================================
!=======================================================================
       subroutine rhsnlhs
!=======================================================================
!
!     SOLVING N-S EQUATION TO GET INTERMEDIATE VELOCITY, U_I HAT
!     INTERMEDIATE VELOCITY, U_I HAT: {U_I}^(K-1) -> U_I HAT -> {U_I}^K
!           IS BASED ON THE PRESSURE AT TIMESTEP K (P^K).
!           INDEED, U_I HAT IS NOT SATISFIED WITH CONTINUITY EQUATION.
!           FOR DETAILS ABOUT U_I HAT, SEE 'FRACTIONAL STEP METHOD'.
!
!     RHS SUB-GRID SCALE PART (LES) -> RHS -> RHS IB(IMMERSED BOUNDARY)
!     -> COMPUTING CONVECTIVE BOUNDARY CONDITION (EXIT CONDITION)
!     -> OTHER BOUNDARY CONDITIONS -> LHS (GET U_I HAT)
!
!     SCHEMES USED IN THIS CODE
!       TIME : 3RD ORDER RUNGE-KUTTA (CONVECTION/NON-LINEAR TERM) +
!              2ND ORDER CRANK-NICOLSON (DIFFUSION/LINEAR TERM)
!       SPACE: 2ND ORDER CENTRAL-DIFFERENCE
!
!     NON-LINEAR TERM: N(U_I)=  D(U_I*U_J)/D(X_J)
!     LINEAR TERM    : L(U_I)=( D(D(U_I))/(D(X_J)D(X_J)) )/RE
!
!
!     SEE THE FOLLOWING PAPERS FOR DETAILS:
!
!     PAPERS FOR FRACTIONAL STEP METHOD (COMPUTATION PROCEDURE):
!      KIM, J., MOIN, P. & MOSER, R. 1987 TURBULENCE STATISTICS IN
!        FULLY DEVELOPED CHANNEL FLOW AT LOW REYNOLDS NUMBER.
!        J. FLUID MECH. 177, 133.
!      CHOI, H., MOIN, P. & KIM, J. 1994 ACTIVE TURBULENCE CONTROL FOR
!        DRAG REDUCTION IN WALL-BOUNDED FLOWS. J. FLUID MECH. 262,
!        75-110.
!
!     PAPER FOR LARGE EDDY SIMULATION (DYNAMIC GLOBAL MODEL):
!      LEE, J., CHOI, H. & PARK, N. 2010 DYNAMIC GLOBAL MODEL FOR
!        LARGE EDDY SIMULATION OF TRANSIENT FLOW. PHYS. FLUIDS 22,
!        075106.
!
!     PAPER FOR IMMERSED BOUNDARY(IB) METHOD:
!      KIM, J., KIM, D. & CHOI, H. 2001 AN IMMERSED-BOUNDARY FINITE
!        VOLUME METHOD FOR SIMULATIONS OF FLOW IN COMPLEX GEOMETRIES.
!        J. COMPUT. PHYS. 171 132-150.
!
!-----------------------------------------------------------------------
         use mod_common
         implicit none

         if (iles .eq. 1) call rhssgs    ! SUBGRID-SCALE COMPUTATIONS (LES CASE)
         call rhs
         if (ibmon .ne. 0) call rhs_ibm    ! IB METHOD COMPUTATION
         call convbc                     ! DEFINED IN LICA_[BODYNAME, E.G. CYLINDER].F90 FILE
         if (ich .ne. 1) call rhsincorpbc

         call lhsu
         call lhsv
         if (n3m .ne. 1) call lhsw

         ! WHEN DO BLOWING/SUCTION RETRV_UVW SHOULD BE TURNED ON
         call retrv_uvw
         call prdic_adj_uvw(0)

         return
       end
!=======================================================================
       subroutine rhs
!=======================================================================
!
!     COMPUTING INTERMEDIATE VELOCITY, U HAT, STEP IN DELTA FORM
!
!     VARIABLES IN COMMON:
!           X, Y, Z         : COORDINATE DIRECTION
!           U, V, W         : VELOCITY FOR X, Y, Z DIRECTION
!           N, S, E, W, C, F: + & - FOR EACH X, Y, Z DIRECTION
!           AN, AL, TAL     : NON-LINEAR, LINEAR, TURBULENT (SGS)
!
!     RHS1(X,Y,Z,L):
!           RHS TERM CONSISTS OF NON-LINEAR/LINEAR/SGS TERMS
!           COMPONENTS FOR IB METHOD WILL BE ADDED IN RHS_IBM SUBROUTINE
!
!-----------------------------------------------------------------------
         use mod_common
         use mod_flowarray, only: u, v, w, p, t, nusgs, nusgs1, rhs1 &
                                  , rk3xo, rk3yo, rk3zo, rk3xoo, rk3yoo, rk3zoo
         implicit none
         integer(8) :: i, j, k
         integer(8) :: ichfac

!------------ VARIABLES FOR U (X, STREAMWISE) VELOCITY
         real(8) :: un, us, uc, uf, vnx, vsx, wcx, wfx
         real(8) :: anxe, anxw, anxx, anx1, anx2, rk3x
         real(8) :: alx1, alx2, alx3, alx4, alx5, alx6, alxx, alxy, alxz, alx
         real(8) :: talxx, talxy, talxz
!------------ VARIABLES FOR V (Y, TRASEVERSE) VELOCITY
         real(8) :: ve, vw, vc, vf, uey, uwy, wcy, wfy
         real(8) :: anyn, anys, anyy, any1, any2, rk3y
         real(8) :: aly1, aly2, aly3, aly4, aly5, aly6, alyx, alyy, alyz, aly
         real(8) :: talyx, talyy, talyz
!------------ VARIABLES FOR W (Z, SPANWISE) VELOCITY
         real(8) :: we, ww, wn, ws, uez, uwz, vnz, vsz
         real(8) :: anzc, anzf, anzz, anz1, anz2, rk3z
         real(8) :: alz1, alz2, alz3, alz4, alz5, alz6, alzx, alzy, alzz, alz
         real(8) :: talzx, talzy, talzz

         if (ich .ne. 0) then
           ichfac = 1
         else
           ichfac = 0
         end if

!-----RHS1 CALCULATION FOR U MOMENTUM -----------------
!$OMP PARALLEL DO  &
!$OMP PRIVATE(UN,US,UC,UF,VNX,VSX,WCX,WFX)   &
!$OMP PRIVATE(ANXE,ANXW,ANXX,ANX1,ANX2,RK3X) &
!$OMP PRIVATE(ALX1,ALX2,ALX3,ALX4,ALX5,ALX6) &
!$OMP PRIVATE(ALXX,ALXY,ALXZ,ALX,TALXX,TALXY,TALXZ)
         do k = 1, n3m
           do j = 1, n2m
             do i = i_bgpx, n1m

               un = 0.5 * (f2fy(j + 1) * u(i, j, k) + f2fy(j) * u(i, j + 1, k)) * c2cyi(j + 1) &
                    * (1.-fixju(j)) + u(i, n2, k) * fixju(j)
               us = 0.5 * (f2fy(j) * u(i, j - 1, k) + f2fy(j - 1) * u(i, j, k)) * c2cyi(j) &
                    * (1.-fixjl(j)) + u(i, 0, k) * fixjl(j)
               uc = 0.5 * (f2fz(k + 1) * u(i, j, k) + f2fz(k) * u(i, j, k + 1)) * c2czi(k + 1) &
                    * (1.-fixku(k)) + u(i, j, n3) * fixku(k)
               uf = 0.5 * (f2fz(k) * u(i, j, k - 1) + f2fz(k - 1) * u(i, j, k)) * c2czi(k) &
                    * (1.-fixkl(k)) + u(i, j, 0) * fixkl(k)
               vnx = 0.5 * (f2fx(i - 1) * v(i, j + 1, k) + f2fx(i) * v(i - 1, j + 1, k)) * c2cxi(i)
               vsx = 0.5 * (f2fx(i - 1) * v(i, j, k) + f2fx(i) * v(i - 1, j, k)) * c2cxi(i)
               wcx = 0.5 * (f2fx(i - 1) * w(i, j, k + 1) + f2fx(i) * w(i - 1, j, k + 1)) * c2cxi(i)
               wfx = 0.5 * (f2fx(i - 1) * w(i, j, k) + f2fx(i) * w(i - 1, j, k)) * c2cxi(i)

               anxe = 0.5 * (u(i + 1, j, k) + u(i, j, k))
               anxw = 0.5 * (u(i, j, k) + u(i - 1, j, k))
               anxx = (anxe**2 - anxw**2) * c2cxi(i)
               anx1 = (un * vnx - us * vsx) * f2fyi(j)
               anx2 = (uc * wcx - uf * wfx) * f2fzi(k)
               rk3x = -anxx - anx1 - anx2               ! NON-LINEAR TERM AT K-SUBSTEP

               alx1 = (u(i + 1, j, k) - u(i, j, k)) * f2fxi(i)
               alx2 = (u(i, j, k) - u(i - 1, j, k)) * f2fxi(i - 1)
               alx3 = (u(i, j + 1, k) - u(i, j, k)) * c2cyi(j + 1) * (1.-fixju(j) * float(jut))
               alx4 = (u(i, j, k) - u(i, j - 1, k)) * c2cyi(j) * (1.-fixjl(j) * float(jub))
               alx5 = (u(i, j, k + 1) - u(i, j, k)) * c2czi(k + 1) * (1.-fixku(k) * float(kut))
               alx6 = (u(i, j, k) - u(i, j, k - 1)) * c2czi(k) * (1.-fixkl(k) * float(kub))
               alxx = (alx1 - alx2) * c2cxi(i)
               alxy = (alx3 - alx4) * f2fyi(j)
               alxz = (alx5 - alx6) * f2fzi(k)
               alx = 1./re * (alxx + alxy + alxz)           ! LINEAR TERMS AT K-SUBSTEP

!-----LES
               if (iles .eq. 1) then
                 talxx = (nusgs(i, j, k) * alx1 - nusgs(i - 1, j, k) * alx2) * c2cxi(i)
                 talxy = (nusgs1(i, j + 1, k, 3) * alx3 - nusgs1(i, j, k, 3) * alx4) * f2fyi(j)
                 talxz = (nusgs1(i, j, k + 1, 2) * alx5 - nusgs1(i, j, k, 2) * alx6) * f2fzi(k)
                 alx = alx + float(iles) * (2.*talxx + talxy + talxz)
                 rk3x = rk3x + float(iles) * rhs1(i, j, k, 1)  ! RHS1 IN HERE: RHSSGS(LICA_SGS.F90)
               end if
!-----LES

!-----HEAT TRANSFER
               if ((ihtrans .eq. 1) .and. (grdir .eq. 1)) then
                 rk3x = rk3x + gr / re**2.*t(i, j, k)
               end if
!-----HEAT TRANSFER

               rhs1(i, j, k, 1) = dt &
                                  * (gamma(msub) * rk3x + ro(msub) * rk3xo(i, j, k) &
                                     + 2.*alpha * alx &
                                     - 2.*alpha * (p(i, j, k) - p(i - 1, j, k)) * c2cxi(i) &
                             - 2.*alpha * pmi(0) * dble(ichfac))
               rk3xoo(i, j, k) = rk3xo(i, j, k)
               rk3xo(i, j, k) = rk3x

             end do
           end do
         end do

!-----RHS1 CALCULATION FOR V MOMENTUM -----------------
!$OMP PARALLEL DO  &
!$OMP PRIVATE(VE,VW,VC,VF,UEY,UWY,WCY,WFY)   &
!$OMP PRIVATE(ANYN,ANYS,ANYY,ANY1,ANY2,RK3Y) &
!$OMP PRIVATE(ALY1,ALY2,ALY3,ALY4,ALY5,ALY6) &
!$OMP PRIVATE(ALYX,ALYY,ALYZ,ALY,TALYX,TALYY,TALYZ)
         do k = 1, n3m
           do j = j_bgpy, n2m
             do i = 1, n1m

               ve = (0.5 * (f2fx(i + 1) * v(i, j, k) + f2fx(i) * v(i + 1, j, k)) * c2cxi(i + 1)) &
                    * (1.-fixiu(i)) + v(n1, j, k) * fixiu(i)
               vw = (0.5 * (f2fx(i) * v(i - 1, j, k) + f2fx(i - 1) * v(i, j, k)) * c2cxi(i)) &
                    * (1.-fixil(i)) + v(0, j, k) * fixil(i)
               vc = (0.5 * (f2fz(k + 1) * v(i, j, k) + f2fz(k) * v(i, j, k + 1)) * c2czi(k + 1)) &
                    * (1.-fixku(k)) + v(i, j, n3) * fixku(k)
               vf = (0.5 * (f2fz(k) * v(i, j, k - 1) + f2fz(k - 1) * v(i, j, k)) * c2czi(k)) &
                    * (1.-fixkl(k)) + v(i, j, 0) * fixkl(k)
               uey = 0.5 * (f2fy(j - 1) * u(i + 1, j, k) + f2fy(j) * u(i + 1, j - 1, k)) * c2cyi(j)
               uwy = 0.5 * (f2fy(j - 1) * u(i, j, k) + f2fy(j) * u(i, j - 1, k)) * c2cyi(j)
               wcy = 0.5 * (f2fy(j - 1) * w(i, j, k + 1) + f2fy(j) * w(i, j - 1, k + 1)) * c2cyi(j)
               wfy = 0.5 * (f2fy(j - 1) * w(i, j, k) + f2fy(j) * w(i, j - 1, k)) * c2cyi(j)

               anyn = 0.5 * (v(i, j + 1, k) + v(i, j, k))
               anys = 0.5 * (v(i, j, k) + v(i, j - 1, k))
               anyy = c2cyi(j) * (anyn**2 - anys**2)
               any1 = (uey * ve - uwy * vw) * f2fxi(i)
               any2 = (vc * wcy - vf * wfy) * f2fzi(k)
               rk3y = -anyy - any1 - any2                ! NY TERM

               aly1 = (v(i + 1, j, k) - v(i, j, k)) * c2cxi(i + 1)
               aly2 = (v(i, j, k) - v(i - 1, j, k)) * c2cxi(i)
               aly3 = (v(i, j + 1, k) - v(i, j, k)) * f2fyi(j)
               aly4 = (v(i, j, k) - v(i, j - 1, k)) * f2fyi(j - 1)
               aly5 = (v(i, j, k + 1) - v(i, j, k)) * c2czi(k + 1) * (1.-fixku(k) * float(kvt))
               aly6 = (v(i, j, k) - v(i, j, k - 1)) * c2czi(k) * (1.-fixkl(k) * float(kvb))
               alyx = (aly1 - aly2) * f2fxi(i)
               alyy = (aly3 - aly4) * c2cyi(j)
               alyz = (aly5 - aly6) * f2fzi(k)
               aly = 1./re * (alyx + alyy + alyz)

!-----LES
               if (iles .eq. 1) then
                 talyx = (nusgs1(i + 1, j, k, 3) * aly1 - nusgs1(i, j, k, 3) * aly2) * f2fxi(i)
                 talyy = (nusgs(i, j, k) * aly3 - nusgs(i, j - 1, k) * aly4) * c2cyi(j)
                 talyz = (nusgs1(i, j, k + 1, 1) * aly5 - nusgs1(i, j, k, 1) * aly6) * f2fzi(k)
                 aly = aly + float(iles) * (talyx + 2.*talyy + talyz)
                 rk3y = rk3y + float(iles) * rhs1(i, j, k, 2)
               end if
!-----LES

!-----HEAT TRANSFER
               if ((ihtrans .eq. 1) .and. (grdir .eq. 2)) then
                 rk3y = rk3y + gr / re**2.*t(i, j, k)
               end if
!-----HEAT TRANSFER

               rhs1(i, j, k, 2) = dt &
                                  * (gamma(msub) * rk3y + ro(msub) * rk3yo(i, j, k) &
                                     + 2.*alpha * aly &
                                     - 2.*alpha * (p(i, j, k) - p(i, j - 1, k)) * c2cyi(j))
               rk3yoo(i, j, k) = rk3yo(i, j, k)  ! ADDED TERM FOR FY
               rk3yo(i, j, k) = rk3y

             end do
           end do
         end do

         if (n3m .eq. 1) goto 100
!-----RHS1 CALCULATION FOR W MOMENTUM -----------------
!$OMP PARALLEL DO &
!$OMP PRIVATE(WE,WW,WN,WS,UEZ,UWZ,VNZ,VSZ)    &
!$OMP PRIVATE(ANZC,ANZF,ANZZ,ANZ1,ANZ2,RK3Z)  &
!$OMP PRIVATE(ALZ1,ALZ2,ALZ3,ALZ4,ALZ5,ALZ6)  &
!$OMP PRIVATE(ALZX,ALZY,ALZZ,ALZ,TALZX,TALZY,TALZZ)
         do k = k_bgpz, n3m
           do j = 1, n2m
             do i = 1, n1m

               we = 0.5 * (f2fx(i + 1) * w(i, j, k) + f2fx(i) * w(i + 1, j, k)) * c2cxi(i + 1) &
                    * (1.-fixiu(i)) + fixiu(i) * w(n1, j, k)
               ww = 0.5 * (f2fx(i) * w(i - 1, j, k) + f2fx(i - 1) * w(i, j, k)) * c2cxi(i) &
                    * (1.-fixil(i)) + fixil(i) * w(0, j, k)
               wn = (0.5 * (f2fy(j + 1) * w(i, j, k) + f2fy(j) * w(i, j + 1, k)) * c2cyi(j + 1)) &
                    * (1.-fixju(j)) + fixju(j) * w(i, n2, k)
               ws = (0.5 * (f2fy(j) * w(i, j - 1, k) + f2fy(j - 1) * w(i, j, k)) * c2cyi(j)) &
                    * (1.-fixjl(j)) + fixjl(j) * w(i, 0, k)
               uez = 0.5 * (f2fz(k - 1) * u(i + 1, j, k) + f2fz(k) * u(i + 1, j, k - 1)) * c2czi(k)
               uwz = 0.5 * (f2fz(k - 1) * u(i, j, k) + f2fz(k) * u(i, j, k - 1)) * c2czi(k)
               vnz = 0.5 * (f2fz(k - 1) * v(i, j + 1, k) + f2fz(k) * v(i, j + 1, k - 1)) * c2czi(k)
               vsz = 0.5 * (f2fz(k - 1) * v(i, j, k) + f2fz(k) * v(i, j, k - 1)) * c2czi(k)

               anzc = 0.5 * (w(i, j, k + 1) + w(i, j, k))
               anzf = 0.5 * (w(i, j, k) + w(i, j, k - 1))
               anzz = (anzc**2 - anzf**2) * c2czi(k)
               anz1 = (uez * we - uwz * ww) * f2fxi(i)
               anz2 = (vnz * wn - vsz * ws) * f2fyi(j)
               rk3z = -anzz - anz1 - anz2                        ! NZ TERM

               alz1 = (w(i + 1, j, k) - w(i, j, k)) * c2cxi(i + 1)
               alz2 = (w(i, j, k) - w(i - 1, j, k)) * c2cxi(i)
               alz3 = (w(i, j + 1, k) - w(i, j, k)) * c2cyi(j + 1) * (1.-fixju(j) * float(jwt))
               alz4 = (w(i, j, k) - w(i, j - 1, k)) * c2cyi(j) * (1.-fixjl(j) * float(jwb))
               alz5 = (w(i, j, k + 1) - w(i, j, k)) * f2fzi(k)
               alz6 = (w(i, j, k) - w(i, j, k - 1)) * f2fzi(k - 1)
               alzx = (alz1 - alz2) * f2fxi(i)
               alzy = (alz3 - alz4) * f2fyi(j)
               alzz = (alz5 - alz6) * c2czi(k)
               alz = 1./re * (alzx + alzy + alzz)

!-----LES
               if (iles .eq. 1) then
                 talzx = (nusgs1(i + 1, j, k, 2) * alz1 - nusgs1(i, j, k, 2) * alz2) * f2fxi(i)
                 talzy = (nusgs1(i, j + 1, k, 1) * alz3 - nusgs1(i, j, k, 1) * alz4) * f2fyi(j)
                 talzz = (nusgs(i, j, k) * alz5 - nusgs(i, j, k - 1) * alz6) * c2czi(k)
                 alz = alz + float(iles) * (talzx + talzy + 2.*talzz)
                 rk3z = rk3z + float(iles) * rhs1(i, j, k, 3)
               end if
!-----LES

!-----HEAT TRANSFER
               if ((ihtrans .eq. 1) .and. (grdir .eq. 3)) then
                 rk3z = rk3z + gr / re**2.*t(i, j, k)
               end if
!-----HEAT TRANSFER

               rhs1(i, j, k, 3) = dt &
                                  * (gamma(msub) * rk3z + ro(msub) * rk3zo(i, j, k) &
                                     + 2.*alpha * alz &
                                     - 2.*alpha * (p(i, j, k) - p(i, j, k - 1)) * c2czi(k))
               rk3zoo(i, j, k) = rk3zo(i, j, k)      ! ADDED TERM FOR FZ
               rk3zo(i, j, k) = rk3z

             end do
           end do
         end do

100      return
       end
!=======================================================================
       subroutine rhs_ibm
!=======================================================================
!
!     CALCULATE MOMENTUM FORCING
!
!     OPTION
!       IMOVINGON = 0, STATIONARY BODY => UBODY,VBODY,WBODY FOR TRANSLATIONAL VEL.
!       IMOVINGON = 1, MOVING BODY     => UBD,VBD,WBD IN LICA_CYLINDER.F90
!
!     VARIABLES
!       UTARG,VTARG,WTARG: TARGET VELOCITY TO SATISFY NO-SLIP B.C.
!       FCV   : MOMENTUM FORCING FROM THE TARGET VELOCITY
!       FCVAVG: AVERAGING FORCING VALUES    TO CALCULATE THE FORCE ON A BODY
!       DUDTR : TIME DERIVATIVE OF VELOCITY TO CALCULATE THE FORCE ON A BODY
!               REF. LEE ET AL., 2011, SOURCES OF SPURIOUS FORCE
!               OSCILLATIONS FROM AN IMMERSED BOUNDARY METHOD FOR MOVING
!               -BODY PROBLEMS, J. COMP. PHYS., 230, 2677-2695.
!
!-----------------------------------------------------------------------
         use mod_common
         use mod_flowarray, only: rhs1, u, v, w, ifc, jfc, kfc, intpindx, geomfac, &
                                  fcv, fcvavg, dudtr
         implicit none
         integer(8) :: i, j, k, n, l, ii, jj, kk, ip, jp, kp, ipp, jpp, kpp
         real(8) :: ubody, vbody, wbody, dti
         real(8) :: utarg, vtarg, wtarg
         real(8) :: rhstmp(n1m, n2m, n3m, 3)
         real(8) :: ubd, vbd, wbd

!---- COMPUTE TARGET VELOCITIES & FORCING VALUES AT FORCING POINTS
         ubody = 0.
         vbody = 0.
         wbody = 0.
         dti = 1./dt

!$OMP PARALLEL DO
         do k = 1, n3m
           do j = 1, n2m
             do i = 1, n1m
               do l = 1, 3
                 rhstmp(i, j, k, l) = rhs1(i, j, k, l)
               end do
             end do
           end do
         end do

!$OMP PARALLEL DO &
!$OMP PRIVATE(II,JJ,KK,IP,JP,KP,IPP,JPP,KPP,UTARG,UBODY)
         do n = 1, nintp(1)
           ii = ifc(n, 1)
           jj = jfc(n, 1)
           kk = kfc(n, 1)
           ip = ifc(n, 1) + intpindx(n, 1, 1)
           jp = jfc(n, 1) + intpindx(n, 1, 2)
           kp = kfc(n, 1) + intpindx(n, 1, 3)
           ipp = ifc(n, 1) + intpindx(n, 1, 1) * 2
           jpp = jfc(n, 1) + intpindx(n, 1, 2) * 2
           kpp = kfc(n, 1) + intpindx(n, 1, 3) * 2

           if ((ipp .ge. n1) .or. (ipp .le. 0)) call reindex_i(ip, ipp)
           if ((kpp .ge. n3) .or. (kpp .le. 0)) call reindex_k(kp, kpp)
           ! IF (IMOVINGON.EQ.1) UBODY=UBD(X(II),YMP(JJ),ZMP(KK))

           utarg = geomfac(n, 1, 0, 0, 0) * ubody &
                   + geomfac(n, 1, 0, 0, 1) * (u(ii, jj, kp) + rhstmp(ii, jj, kp, 1)) &
                   + geomfac(n, 1, 0, 0, 2) * (u(ii, jj, kpp) + rhstmp(ii, jj, kpp, 1)) &
                   + geomfac(n, 1, 0, 1, 0) * (u(ii, jp, kk) + rhstmp(ii, jp, kk, 1)) &
                   + geomfac(n, 1, 0, 1, 1) * (u(ii, jp, kp) + rhstmp(ii, jp, kp, 1)) &
                   + geomfac(n, 1, 0, 1, 2) * (u(ii, jp, kpp) + rhstmp(ii, jp, kpp, 1)) &
                   + geomfac(n, 1, 0, 2, 0) * (u(ii, jpp, kk) + rhstmp(ii, jpp, kk, 1)) &
                   + geomfac(n, 1, 0, 2, 1) * (u(ii, jpp, kp) + rhstmp(ii, jpp, kp, 1)) &
                   + geomfac(n, 1, 0, 2, 2) * (u(ii, jpp, kpp) + rhstmp(ii, jpp, kpp, 1)) &
                   + geomfac(n, 1, 1, 0, 0) * (u(ip, jj, kk) + rhstmp(ip, jj, kk, 1)) &
                   + geomfac(n, 1, 1, 0, 1) * (u(ip, jj, kp) + rhstmp(ip, jj, kp, 1)) &
                   + geomfac(n, 1, 1, 0, 2) * (u(ip, jj, kpp) + rhstmp(ip, jj, kpp, 1)) &
                   + geomfac(n, 1, 1, 1, 0) * (u(ip, jp, kk) + rhstmp(ip, jp, kk, 1)) &
                   + geomfac(n, 1, 1, 1, 1) * (u(ip, jp, kp) + rhstmp(ip, jp, kp, 1)) &
                   + geomfac(n, 1, 1, 1, 2) * (u(ip, jp, kpp) + rhstmp(ip, jp, kpp, 1)) &
                   + geomfac(n, 1, 1, 2, 0) * (u(ip, jpp, kk) + rhstmp(ip, jpp, kk, 1)) &
                   + geomfac(n, 1, 1, 2, 1) * (u(ip, jpp, kp) + rhstmp(ip, jpp, kp, 1)) &
                   + geomfac(n, 1, 1, 2, 2) * (u(ip, jpp, kpp) + rhstmp(ip, jpp, kpp, 1)) &
                   + geomfac(n, 1, 2, 0, 0) * (u(ipp, jj, kk) + rhstmp(ipp, jj, kk, 1)) &
                   + geomfac(n, 1, 2, 0, 1) * (u(ipp, jj, kp) + rhstmp(ipp, jj, kp, 1)) &
                   + geomfac(n, 1, 2, 0, 2) * (u(ipp, jj, kpp) + rhstmp(ipp, jj, kpp, 1)) &
                   + geomfac(n, 1, 2, 1, 0) * (u(ipp, jp, kk) + rhstmp(ipp, jp, kk, 1)) &
                   + geomfac(n, 1, 2, 1, 1) * (u(ipp, jp, kp) + rhstmp(ipp, jp, kp, 1)) &
                   + geomfac(n, 1, 2, 1, 2) * (u(ipp, jp, kpp) + rhstmp(ipp, jp, kpp, 1)) &
                   + geomfac(n, 1, 2, 2, 0) * (u(ipp, jpp, kk) + rhstmp(ipp, jpp, kk, 1)) &
                   + geomfac(n, 1, 2, 2, 1) * (u(ipp, jpp, kp) + rhstmp(ipp, jpp, kp, 1)) &
                   + geomfac(n, 1, 2, 2, 2) * (u(ipp, jpp, kpp) + rhstmp(ipp, jpp, kpp, 1))
           fcv(n, 1) = (utarg - u(ii, jj, kk) - rhstmp(ii, jj, kk, 1)) * dti
           dudtr(n, 1) = (utarg - u(ii, jj, kk)) * dti
           rhs1(ii, jj, kk, 1) = utarg - u(ii, jj, kk)
         end do

!$OMP PARALLEL DO &
!$OMP PRIVATE(II,JJ,KK,IP,JP,KP,IPP,JPP,KPP,VTARG,VBODY)
         do n = 1, nintp(2)
           ii = ifc(n, 2)
           jj = jfc(n, 2)
           kk = kfc(n, 2)
           ip = ifc(n, 2) + intpindx(n, 2, 1)
           jp = jfc(n, 2) + intpindx(n, 2, 2)
           kp = kfc(n, 2) + intpindx(n, 2, 3)
           ipp = ifc(n, 2) + intpindx(n, 2, 1) * 2
           jpp = jfc(n, 2) + intpindx(n, 2, 2) * 2
           kpp = kfc(n, 2) + intpindx(n, 2, 3) * 2

           if ((ipp .ge. n1) .or. (ipp .le. 0)) call reindex_i(ip, ipp)
           if ((kpp .ge. n3) .or. (kpp .le. 0)) call reindex_k(kp, kpp)
           ! IF (IMOVINGON.EQ.1) VBODY=VBD(XMP(II),Y(JJ),ZMP(KK))

           vtarg = geomfac(n, 2, 0, 0, 0) * vbody &
                   + geomfac(n, 2, 0, 0, 1) * (v(ii, jj, kp) + rhstmp(ii, jj, kp, 2)) &
                   + geomfac(n, 2, 0, 0, 2) * (v(ii, jj, kpp) + rhstmp(ii, jj, kpp, 2)) &
                   + geomfac(n, 2, 0, 1, 0) * (v(ii, jp, kk) + rhstmp(ii, jp, kk, 2)) &
                   + geomfac(n, 2, 0, 1, 1) * (v(ii, jp, kp) + rhstmp(ii, jp, kp, 2)) &
                   + geomfac(n, 2, 0, 1, 2) * (v(ii, jp, kpp) + rhstmp(ii, jp, kpp, 2)) &
                   + geomfac(n, 2, 0, 2, 0) * (v(ii, jpp, kk) + rhstmp(ii, jpp, kk, 2)) &
                   + geomfac(n, 2, 0, 2, 1) * (v(ii, jpp, kp) + rhstmp(ii, jpp, kp, 2)) &
                   + geomfac(n, 2, 0, 2, 2) * (v(ii, jpp, kpp) + rhstmp(ii, jpp, kpp, 2)) &
                   + geomfac(n, 2, 1, 0, 0) * (v(ip, jj, kk) + rhstmp(ip, jj, kk, 2)) &
                   + geomfac(n, 2, 1, 0, 1) * (v(ip, jj, kp) + rhstmp(ip, jj, kp, 2)) &
                   + geomfac(n, 2, 1, 0, 2) * (v(ip, jj, kpp) + rhstmp(ip, jj, kpp, 2)) &
                   + geomfac(n, 2, 1, 1, 0) * (v(ip, jp, kk) + rhstmp(ip, jp, kk, 2)) &
                   + geomfac(n, 2, 1, 1, 1) * (v(ip, jp, kp) + rhstmp(ip, jp, kp, 2)) &
                   + geomfac(n, 2, 1, 1, 2) * (v(ip, jp, kpp) + rhstmp(ip, jp, kpp, 2)) &
                   + geomfac(n, 2, 1, 2, 0) * (v(ip, jpp, kk) + rhstmp(ip, jpp, kk, 2)) &
                   + geomfac(n, 2, 1, 2, 1) * (v(ip, jpp, kp) + rhstmp(ip, jpp, kp, 2)) &
                   + geomfac(n, 2, 1, 2, 2) * (v(ip, jpp, kpp) + rhstmp(ip, jpp, kpp, 2)) &
                   + geomfac(n, 2, 2, 0, 0) * (v(ipp, jj, kk) + rhstmp(ipp, jj, kk, 2)) &
                   + geomfac(n, 2, 2, 0, 1) * (v(ipp, jj, kp) + rhstmp(ipp, jj, kp, 2)) &
                   + geomfac(n, 2, 2, 0, 2) * (v(ipp, jj, kpp) + rhstmp(ipp, jj, kpp, 2)) &
                   + geomfac(n, 2, 2, 1, 0) * (v(ipp, jp, kk) + rhstmp(ipp, jp, kk, 2)) &
                   + geomfac(n, 2, 2, 1, 1) * (v(ipp, jp, kp) + rhstmp(ipp, jp, kp, 2)) &
                   + geomfac(n, 2, 2, 1, 2) * (v(ipp, jp, kpp) + rhstmp(ipp, jp, kpp, 2)) &
                   + geomfac(n, 2, 2, 2, 0) * (v(ipp, jpp, kk) + rhstmp(ipp, jpp, kk, 2)) &
                   + geomfac(n, 2, 2, 2, 1) * (v(ipp, jpp, kp) + rhstmp(ipp, jpp, kp, 2)) &
                   + geomfac(n, 2, 2, 2, 2) * (v(ipp, jpp, kpp) + rhstmp(ipp, jpp, kpp, 2))
           fcv(n, 2) = (vtarg - v(ii, jj, kk) - rhstmp(ii, jj, kk, 2)) * dti
           dudtr(n, 2) = (vtarg - v(ii, jj, kk)) * dti
           rhs1(ii, jj, kk, 2) = vtarg - v(ii, jj, kk)
         end do

!$OMP PARALLEL DO &
!$OMP PRIVATE(II,JJ,KK,IP,JP,KP,IPP,JPP,KPP,WTARG,WBODY)
         do n = 1, nintp(3)
           ii = ifc(n, 3)
           jj = jfc(n, 3)
           kk = kfc(n, 3)
           ip = ifc(n, 3) + intpindx(n, 3, 1)
           jp = jfc(n, 3) + intpindx(n, 3, 2)
           kp = kfc(n, 3) + intpindx(n, 3, 3)
           ipp = ifc(n, 3) + intpindx(n, 3, 1) * 2
           jpp = jfc(n, 3) + intpindx(n, 3, 2) * 2
           kpp = kfc(n, 3) + intpindx(n, 3, 3) * 2

           if ((ipp .ge. n1) .or. (ipp .le. 0)) call reindex_i(ip, ipp)
           if ((kpp .ge. n3) .or. (kpp .le. 0)) call reindex_k(kp, kpp)
           ! IF (IMOVINGON.EQ.1) WBODY=WBD(XMP(II),YMP(JJ),Z(KK))

           wtarg = geomfac(n, 3, 0, 0, 0) * wbody &
                   + geomfac(n, 3, 0, 0, 1) * (w(ii, jj, kp) + rhstmp(ii, jj, kp, 3)) &
                   + geomfac(n, 3, 0, 0, 2) * (w(ii, jj, kpp) + rhstmp(ii, jj, kpp, 3)) &
                   + geomfac(n, 3, 0, 1, 0) * (w(ii, jp, kk) + rhstmp(ii, jp, kk, 3)) &
                   + geomfac(n, 3, 0, 1, 1) * (w(ii, jp, kp) + rhstmp(ii, jp, kp, 3)) &
                   + geomfac(n, 3, 0, 1, 2) * (w(ii, jp, kpp) + rhstmp(ii, jp, kpp, 3)) &
                   + geomfac(n, 3, 0, 2, 0) * (w(ii, jpp, kk) + rhstmp(ii, jpp, kk, 3)) &
                   + geomfac(n, 3, 0, 2, 1) * (w(ii, jpp, kp) + rhstmp(ii, jpp, kp, 3)) &
                   + geomfac(n, 3, 0, 2, 2) * (w(ii, jpp, kpp) + rhstmp(ii, jpp, kpp, 3)) &
                   + geomfac(n, 3, 1, 0, 0) * (w(ip, jj, kk) + rhstmp(ip, jj, kk, 3)) &
                   + geomfac(n, 3, 1, 0, 1) * (w(ip, jj, kp) + rhstmp(ip, jj, kp, 3)) &
                   + geomfac(n, 3, 1, 0, 2) * (w(ip, jj, kpp) + rhstmp(ip, jj, kpp, 3)) &
                   + geomfac(n, 3, 1, 1, 0) * (w(ip, jp, kk) + rhstmp(ip, jp, kk, 3)) &
                   + geomfac(n, 3, 1, 1, 1) * (w(ip, jp, kp) + rhstmp(ip, jp, kp, 3)) &
                   + geomfac(n, 3, 1, 1, 2) * (w(ip, jp, kpp) + rhstmp(ip, jp, kpp, 3)) &
                   + geomfac(n, 3, 1, 2, 0) * (w(ip, jpp, kk) + rhstmp(ip, jpp, kk, 3)) &
                   + geomfac(n, 3, 1, 2, 1) * (w(ip, jpp, kp) + rhstmp(ip, jpp, kp, 3)) &
                   + geomfac(n, 3, 1, 2, 2) * (w(ip, jpp, kpp) + rhstmp(ip, jpp, kpp, 3)) &
                   + geomfac(n, 3, 2, 0, 0) * (w(ipp, jj, kk) + rhstmp(ipp, jj, kk, 3)) &
                   + geomfac(n, 3, 2, 0, 1) * (w(ipp, jj, kp) + rhstmp(ipp, jj, kp, 3)) &
                   + geomfac(n, 3, 2, 0, 2) * (w(ipp, jj, kpp) + rhstmp(ipp, jj, kpp, 3)) &
                   + geomfac(n, 3, 2, 1, 0) * (w(ipp, jp, kk) + rhstmp(ipp, jp, kk, 3)) &
                   + geomfac(n, 3, 2, 1, 1) * (w(ipp, jp, kp) + rhstmp(ipp, jp, kp, 3)) &
                   + geomfac(n, 3, 2, 1, 2) * (w(ipp, jp, kpp) + rhstmp(ipp, jp, kpp, 3)) &
                   + geomfac(n, 3, 2, 2, 0) * (w(ipp, jpp, kk) + rhstmp(ipp, jpp, kk, 3)) &
                   + geomfac(n, 3, 2, 2, 1) * (w(ipp, jpp, kp) + rhstmp(ipp, jpp, kp, 3)) &
                   + geomfac(n, 3, 2, 2, 2) * (w(ipp, jpp, kpp) + rhstmp(ipp, jpp, kpp, 3))
           fcv(n, 3) = (wtarg - w(ii, jj, kk) - rhstmp(ii, jj, kk, 3)) * dti
           dudtr(n, 3) = (wtarg - w(ii, jj, kk)) * dti
           rhs1(ii, jj, kk, 3) = wtarg - w(ii, jj, kk)
         end do

!*************************       CAUTION       *************************
!     FOLLOWING 3 DO-LOOPS REQUIRES THE OPTION OF
!     '-WF "-PVCTL VWORK=STACK"'
!     WHEN THE CODE IS COINTPINDXLED ON NEC MACHINE.
!     FOR DETAILED EXPLANATION OF THE REASON,
!     YOU CAN CONSULT "NEC PORTING GUIDE"
!     PROVIDED BY NEC SUPERCOMPUTIONG CENTER.
!***********************************************************************

!$OMP PARALLEL DO &
!$OMP PRIVATE(II,JJ,KK,UTARG,UBODY)
         do n = nintp(1) + 1, nbody(1)
           ii = ifc(n, 1)
           jj = jfc(n, 1)
           kk = kfc(n, 1)
           ! IF (IMOVINGON.EQ.1) UBODY=UBD(X(II),YMP(JJ),ZMP(KK))
           utarg = ubody
           fcv(n, 1) = (utarg - u(ii, jj, kk) - rhs1(ii, jj, kk, 1)) * dti
           dudtr(n, 1) = (utarg - u(ii, jj, kk)) * dti
           rhs1(ii, jj, kk, 1) = utarg - u(ii, jj, kk)
         end do

!$OMP PARALLEL DO &
!$OMP PRIVATE(II,JJ,KK,VTARG,VBODY)
         do n = nintp(2) + 1, nbody(2)
           ii = ifc(n, 2)
           jj = jfc(n, 2)
           kk = kfc(n, 2)
           ! IF (IMOVINGON.EQ.1) VBODY=VBD(XMP(II),Y(JJ),ZMP(KK))
           vtarg = vbody
           fcv(n, 2) = (vtarg - v(ii, jj, kk) - rhs1(ii, jj, kk, 2)) * dti
           dudtr(n, 2) = (vtarg - v(ii, jj, kk)) * dti
           rhs1(ii, jj, kk, 2) = vtarg - v(ii, jj, kk)
         end do

!$OMP PARALLEL DO &
!$OMP PRIVATE(II,JJ,KK,WTARG,WBODY)
         do n = nintp(3) + 1, nbody(3)
           ii = ifc(n, 3)
           jj = jfc(n, 3)
           kk = kfc(n, 3)
           ! IF (IMOVINGON.EQ.1) WBODY=WBD(XMP(II),YMP(JJ),Z(KK))
           wtarg = wbody
           fcv(n, 3) = (wtarg - w(ii, jj, kk) - rhs1(ii, jj, kk, 3)) * dti
           dudtr(n, 3) = (wtarg - w(ii, jj, kk)) * dti
           rhs1(ii, jj, kk, 3) = wtarg - w(ii, jj, kk)
         end do

!-----COMPUTE AVERAGE FORCING VALUES DURING RK3 STEPS
!$OMP PARALLEL DO
         do l = 1, 3
           do n = 1, nbody(l)
             fcvavg(n, l) = fcvavg(n, l) + fcv(n, l)
           end do
         end do

         return
       end
!=======================================================================
       subroutine reindex_i(ip, ipp)
!=======================================================================
!
!     ASSIGN APPROPRIATE INDICES TO THE VARIABLES USED IN IBM INTERPOLATION
!     IN ORDER TO SATISFY THE PERIODIC BC IN X-DIRECTION.
!      IP =I +1
!      IPP=IP+1
!     *GENERALLY, (IPP.GE.N1)~ SITUATIONS OCCUR IN CASE OF PERIODIC BC
!
!-----------------------------------------------------------------------
         use mod_common
         implicit none
         integer(8) :: ip, ipp

         if (ipp .eq. n1) then
           ipp = 1
         elseif (ipp .eq. n1 + 1) then
           ip = 1
           ipp = 2
         elseif (ipp .eq. 0) then
           ipp = n1m
         elseif (ipp .eq. -1) then
           ip = n1m
           ipp = n1m - 1
         end if

         return
       end subroutine reindex_i
!=======================================================================
!=======================================================================
       subroutine reindex_k(kp, kpp)
!=======================================================================
         use mod_common
         implicit none
         integer(8) :: kp, kpp

         if (kpp .eq. n3) then
           kpp = 1
         elseif (kpp .eq. n3 + 1) then
           kp = 1
           kpp = 2
         elseif (kpp .eq. 0) then
           kpp = n3m
         elseif (kpp .eq. -1) then
           kp = n3m
           kpp = n3m - 1
         end if

         return
       end subroutine reindex_k
!=======================================================================
!=======================================================================
       subroutine convbc
!=======================================================================
         use mod_common
         use mod_flowarray, only: u, v, w
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
               if (funcbody(xmp(n1m), ymp(j), zmp(k), time) .le. 1.e-10 .and. imovingon .eq. 0) then
                 ! SOLID PHASE: ZERO-GRADIENT (PRESERVES SOLID VELOCITY = 0)
                 uout(j, k) = u(n1m, j, k)
                 vout(j, k) = v(n1m, j, k)
                 wout(j, k) = w(n1m, j, k)
               else
                 ! FLUID PHASE: CONVECTIVE ADVECTION
                 uout(j, k) = u(n1, j, k) - ucoef * (u(n1, j, k) - u(n1m, j, k))
                 vout(j, k) = v(n1, j, k) - vwcoef * (v(n1, j, k) - v(n1m, j, k))
                 wout(j, k) = w(n1, j, k) - vwcoef * (w(n1, j, k) - w(n1m, j, k))
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
                 uout(j, k) = uout(j, k) * qratio
                 vout(j, k) = vout(j, k) * qratio
                 wout(j, k) = wout(j, k) * qratio
               end if
               duout(j, k) = uout(j, k) - u(n1, j, k)
               dvout(j, k) = vout(j, k) - v(n1, j, k)
               dwout(j, k) = wout(j, k) - w(n1, j, k)
             end do
           end do

           if (yprdic .eq. 1) then
!$OMP PARALLEL DO
             do k = 0, n3
               duout(n2, k) = duout(1, k)
               duout(0, k) = duout(n2m, k)
               dvout(n2, k) = dvout(1, k)
               dvout(0, k) = dvout(n2m, k)
               dwout(n2, k) = dwout(1, k)
               dwout(0, k) = dwout(n2m, k)
             end do
           end if

           if (zprdic .eq. 1) then
!$OMP PARALLEL DO
             do j = 0, n2
               duout(j, n3) = duout(j, 1)
               duout(j, 0) = duout(j, n3m)
               dvout(j, n3) = dvout(j, 1)
               dvout(j, 0) = dvout(j, n3m)
               dwout(j, n3) = dwout(j, 1)
               dwout(j, 0) = dwout(j, n3m)
             end do
           end if

         else

           uout = 0.
           vout = 0.
           wout = 0.

         end if

         return
       end subroutine convbc
!=======================================================================
!=======================================================================
       subroutine rhsincorpbc
!=======================================================================
!
!     APPLYING EXIT BOUNDARY CONDITION & TOP BOUNDARY CONDITION TO RHS
!           [RHS INCORPORATE B.C.]
!
!     ACOEF: MAIN SOLVER IN THIS FILE
!     CIU, CIVW, CJV: SUBROUTINE LHSINIT IN THIS FILE
!     DUOUT, DVOUT, DVTOP, DWOUT: SUBROUTINE CONVBC (LICA_CYLINDER.F90)
!
!-----------------------------------------------------------------------
         use mod_common
         use mod_flowarray, only: u, v, w, rhs1, nusgs, nusgs1
         implicit none
         integer(8) :: i, j, k
         real(8) :: cre, cre2

         if (iles .eq. 1) then
           cre = re
           cre2 = 2.*re
         else
           cre = 0.
           cre2 = 0.
         end if

!$OMP PARALLEL DO
         do k = 1, n3m
           do j = 1, n2m
             rhs1(n1m, j, k, 1) = rhs1(n1m, j, k, 1) &
                                  - acoef * ciu(n1m) * (1.+cre2 * nusgs(n1m, j, k)) * duout(j, k)
           end do
         end do

!$OMP PARALLEL DO
         do k = 1, n3m
           do j = 2, n2m
             rhs1(n1m, j, k, 2) = rhs1(n1m, j, k, 2) &
                                  - acoef * civw(n1m) * (1.+cre * nusgs1(n1, j, k, 3)) * dvout(j, k)
           end do
         end do

!$OMP PARALLEL DO
         do k = k_bgpz, n3m
           do j = 1, n2m
             rhs1(n1m, j, k, 3) = rhs1(n1m, j, k, 3) &
                                  - acoef * civw(n1m) * (1.+cre * nusgs1(n1, j, k, 2)) * dwout(j, k)
           end do
         end do

         return
       end subroutine rhsincorpbc
!=======================================================================
!=======================================================================
       subroutine lhsu
!=======================================================================
!
!     CALCULATE INTERMEDIATE VELOCITY, U_I HAT THROUGH TDMA
!           IN THIS ROUTINE, COMPUTE STREAMWISE VELOCITY (U HAT)
!
!-----------------------------------------------------------------------
         use mod_common
         use mod_flowarray, only: u, v, w, rhs1, nusgs, nusgs1
         implicit none
         integer(8) :: i, j, k
         real(8) :: cre, cre2
         real(8), dimension(:, :), allocatable :: ai, bi, ci, gi
         real(8), dimension(:, :), allocatable :: aj, bj, cj, gj
         real(8), dimension(:, :), allocatable :: ak, bk, ck, gk

         cre = float(iles) * re
         cre2 = 2.*float(iles) * re

!=====ADI STARTS

         if (n3m .eq. 1) goto 100
!-----Z-DIRECTION
!$OMP PARALLEL &
!$OMP PRIVATE(AK,CK,BK,GK)
         allocate (ak(n1, n3), bk(n1, n3), ck(n1, n3), gk(n1, n3))
!$OMP DO
         do j = 1, n2m
           do k = 1, n3m
             do i = i_bgpx, n1m
               ak(i, k) = akuv(k) * (1.+cre * nusgs1(i, j, k, 2)) * (1.-fixkl(k) * float(kub))
               ck(i, k) = ckuv(k) * (1.+cre * nusgs1(i, j, k + 1, 2)) * (1.-fixku(k) * float(kut))
               bk(i, k) = acoefi - ak(i, k) - ck(i, k)
               gk(i, k) = acoefi * rhs1(i, j, k, 1)
             end do
           end do

           if (zprdic .eq. 0) then
             call trdiag3(ak, bk, ck, gk, gk, 1, n3m, i_bgpx, n1m)
           else if (zprdic .eq. 1) then
             call trdiag3p(ak, bk, ck, gk, 1, n3m, i_bgpx, n1m)  !Z PERIODICITY
           end if

           do k = 1, n3m
             do i = i_bgpx, n1m
               rhs1(i, j, k, 1) = gk(i, k)
             end do
           end do
         end do
!$OMP END DO
         deallocate (ak, bk, ck, gk)
!$OMP END PARALLEL

100      continue

!$OMP PARALLEL  &
!$OMP PRIVATE(AJ,CJ,BJ,GJ)  &
!$OMP PRIVATE(AI,CI,BI,GI)
         allocate (ai(n2, n1), bi(n2, n1), ci(n2, n1), gi(n2, n1))
         allocate (aj(n1, n2), bj(n1, n2), cj(n1, n2), gj(n1, n2))
!$OMP DO
         do k = 1, n3m

!-----Y-DIRECTION
           do j = 1, n2m
             do i = i_bgpx, n1m
               aj(i, j) = ajuw(j) * (1.+cre * nusgs1(i, j, k, 3)) * (1.-fixjl(j) * float(jub))
               cj(i, j) = cjuw(j) * (1.+cre * nusgs1(i, j + 1, k, 3)) * (1.-fixju(j) * float(jut))
               bj(i, j) = acoefi - aj(i, j) - cj(i, j)
               gj(i, j) = acoefi * rhs1(i, j, k, 1)
             end do
           end do

           call trdiag2(aj, bj, cj, gj, gj, 1, n2m, i_bgpx, n1m)

!-----X-DIRECTION
           do i = i_bgpx, n1m
             do j = 1, n2m
               ai(j, i) = aiu(i) * (1.+cre2 * nusgs(i - 1, j, k))
               ci(j, i) = ciu(i) * (1.+cre2 * nusgs(i, j, k))
               bi(j, i) = acoefi - ai(j, i) - ci(j, i)
               gi(j, i) = acoefi * gj(i, j)
             end do
           end do

           if (xprdic .eq. 0) then
             call trdiag1(ai, bi, ci, gi, gi, i_bgpx, n1m, 1, n2m)
           else if (xprdic .eq. 1) then
             call trdiag1p(ai, bi, ci, gi, i_bgpx, n1m, 1, n2m)
           end if

           do j = 1, n2m
             do i = i_bgpx, n1m
               u(i, j, k) = gi(j, i) + u(i, j, k)
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
       subroutine lhsv
!=======================================================================
!
!     CALCULATE INTERMEDIATE VELOCITY, U_I HAT THROUGH TDMA
!           IN THIS ROUTINE, COMPUTE TRANSVERSE VELOCITY (V HAT)
!
!-----------------------------------------------------------------------
         use mod_common
         use mod_flowarray, only: u, v, w, rhs1, nusgs, nusgs1
         implicit none
         integer(8) :: i, j, k
         real(8) :: cre, cre2
         real(8), dimension(:, :), allocatable :: ai, bi, ci, gi
         real(8), dimension(:, :), allocatable :: aj, bj, cj, gj
         real(8), dimension(:, :), allocatable :: ak, bk, ck, gk

         cre = float(iles) * re
         cre2 = 2.*float(iles) * re

!=====ADI STARTS

         if (n3m .eq. 1) goto 100
!-----Z-DIRECTION
!$OMP PARALLEL &
!$OMP PRIVATE(AK,CK,BK,GK)
         allocate (ak(n1, n3), bk(n1, n3), ck(n1, n3), gk(n1, n3))
!$OMP DO
         do j = 2, n2m

           do k = 1, n3m
             do i = 1, n1m
               ak(i, k) = akuv(k) * (1.+cre * nusgs1(i, j, k, 1)) * (1.-fixkl(k) * float(kvb))
               ck(i, k) = ckuv(k) * (1.+cre * nusgs1(i, j, k + 1, 1)) * (1.-fixku(k) * float(kvt))
               bk(i, k) = acoefi - ak(i, k) - ck(i, k)
               gk(i, k) = acoefi * rhs1(i, j, k, 2)
             end do
           end do

           if (zprdic .eq. 0) then
             call trdiag3(ak, bk, ck, gk, gk, 1, n3m, 1, n1m)
           else if (zprdic .eq. 1) then
             call trdiag3p(ak, bk, ck, gk, 1, n3m, 1, n1m)
           end if

           do k = 1, n3m
             do i = 1, n1m
               rhs1(i, j, k, 2) = gk(i, k)
             end do
           end do
         end do
!$OMP END DO
         deallocate (ak, bk, ck, gk)
!$OMP END PARALLEL

100      continue

!$OMP PARALLEL &
!$OMP PRIVATE(AJ,CJ,BJ,GJ) &
!$OMP PRIVATE(AI,CI,BI,GI)
         allocate (ai(n2, n1), bi(n2, n1), ci(n2, n1), gi(n2, n1))
         allocate (aj(n1, n2), bj(n1, n2), cj(n1, n2), gj(n1, n2))
!$OMP DO
         do k = 1, n3m

!-----Y-DIRECTION
           do j = 2, n2m
             do i = 1, n1m
               aj(i, j) = ajv(j) * (1.+cre2 * nusgs(i, j - 1, k))
               cj(i, j) = cjv(j) * (1.+cre2 * nusgs(i, j, k))
               bj(i, j) = acoefi - aj(i, j) - cj(i, j)
               gj(i, j) = acoefi * rhs1(i, j, k, 2)
             end do
           end do

           call trdiag2(aj, bj, cj, gj, gj, 2, n2m, 1, n1m)

!-----X-DIRECTION
           do i = 1, n1m
             do j = 2, n2m
               ai(j, i) = aivw(i) * (1.+cre * nusgs1(i, j, k, 3))
               ci(j, i) = civw(i) * (1.+cre * nusgs1(i + 1, j, k, 3))
               bi(j, i) = acoefi - ai(j, i) - ci(j, i)
               gi(j, i) = acoefi * gj(i, j)
             end do
           end do

           if (xprdic .eq. 0) then
             call trdiag1(ai, bi, ci, gi, gi, 1, n1m, 2, n2m)
           else if (xprdic .eq. 1) then
             call trdiag1p(ai, bi, ci, gi, 1, n1m, 2, n2m)
           end if

           do j = 2, n2m
             do i = 1, n1m
               v(i, j, k) = gi(j, i) + v(i, j, k)
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
       subroutine lhsw
!=======================================================================
!
!     CALCULATE INTERMEDIATE VELOCITY, U_I HAT THROUGH TDMA
!           IN THIS ROUTINE, COMPUTE SPANWISE VELOCITY (W HAT)
!
!-----------------------------------------------------------------------
         use mod_common
         use mod_flowarray, only: u, v, w, rhs1, nusgs, nusgs1
         implicit none
         integer(8) :: i, j, k
         real(8) :: cre, cre2
         real(8), dimension(:, :), allocatable :: ai, bi, ci, gi
         real(8), dimension(:, :), allocatable :: aj, bj, cj, gj
         real(8), dimension(:, :), allocatable :: ak, bk, ck, gk

         cre = float(iles) * re
         cre2 = 2.*float(iles) * re

!=====ADI STARTS

!-----Z-DIRECTION
!$OMP PARALLEL &
!$OMP PRIVATE(AK,CK,BK,GK)
         allocate (ak(n1, n3), bk(n1, n3), ck(n1, n3), gk(n1, n3))
!$OMP DO
         do j = 1, n2m

           do k = k_bgpz, n3m
             do i = 1, n1m
               ak(i, k) = akw(k) * (1.+cre2 * nusgs(i, j, k - 1))
               ck(i, k) = ckw(k) * (1.+cre2 * nusgs(i, j, k))
               bk(i, k) = acoefi - ak(i, k) - ck(i, k)
               gk(i, k) = acoefi * rhs1(i, j, k, 3)
             end do
           end do

           if (zprdic .eq. 0) then
             call trdiag3(ak, bk, ck, gk, gk, 2, n3m, 1, n1m)
           else if (zprdic .eq. 1) then
             call trdiag3p(ak, bk, ck, gk, 1, n3m, 1, n1m)
           end if

           do k = k_bgpz, n3m
             do i = 1, n1m
               rhs1(i, j, k, 3) = gk(i, k)
             end do
           end do

         end do
!$OMP END DO
         deallocate (ak, bk, ck, gk)
!$OMP END PARALLEL

!$OMP PARALLEL  &
!$OMP PRIVATE(AJ,CJ,BJ,GJ)  &
!$OMP PRIVATE(AI,CI,BI,GI)
         allocate (ai(n2, n1), bi(n2, n1), ci(n2, n1), gi(n2, n1))
         allocate (aj(n1, n2), bj(n1, n2), cj(n1, n2), gj(n1, n2))
!$OMP DO
         do k = k_bgpz, n3m

!-----Y-DIRECTION
           do j = 1, n2m
             do i = 1, n1m
               aj(i, j) = ajuw(j) * (1.+cre * nusgs1(i, j, k, 1)) * (1.-fixjl(j) * float(jwb))
               cj(i, j) = cjuw(j) * (1.+cre * nusgs1(i, j + 1, k, 1)) * (1.-fixju(j) * float(jwt))
               bj(i, j) = acoefi - aj(i, j) - cj(i, j)
               gj(i, j) = acoefi * rhs1(i, j, k, 3)
             end do
           end do

           call trdiag2(aj, bj, cj, gj, gj, 1, n2m, 1, n1m)

!-----X-DIRECTION
           do i = 1, n1m
             do j = 1, n2m
               ai(j, i) = aivw(i) * (1.+cre * nusgs1(i, j, k, 2))
               ci(j, i) = civw(i) * (1.+cre * nusgs1(i + 1, j, k, 2))
               bi(j, i) = acoefi - ai(j, i) - ci(j, i)
               gi(j, i) = acoefi * gj(i, j)
             end do
           end do

           if (xprdic .eq. 0) then
             call trdiag1(ai, bi, ci, gi, gi, 1, n1m, 1, n2m)
           else if (xprdic .eq. 1) then
             call trdiag1p(ai, bi, ci, gi, 1, n1m, 1, n2m)
           end if

           do j = 1, n2m
             do i = 1, n1m
               w(i, j, k) = gi(j, i) + w(i, j, k)
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
       subroutine trdiag1(a, b, c, r, uu, l1, l2, ll1, ll2)
!=======================================================================
!
!     SOLVE THE TRIDIAGONAL MATRIX (X-1 TRIDIAGONAL) WITH
!     DIFFERENT COEFFICIENTS IN LHSU,LHSV,LHSW SUBROUTINES
!
!     VARIABLES
!       A(),B(),C()  : COEFFICIENTS OF TRIDIAGONAL MATRIX
!       R()          : RHS
!       UU()         : SOLUTION
!       L1,L2  : 2ND INDEX INDICATING X1-DIRECTION (TDMA)
!       LL2,LL2: 1ST INDEX INDICATING X2-DIRECTION
!
!       AT EACH I,
!       |B_{L1}   C_{L1}   0        0               | = |R_{L1  }|
!       |A_{L1+1} B_{L1+1} C_{L1+1} 0               | = |R_{L1+1}|
!       |0        A_{L1+2} B_{L1+2} C_{L1+2}        | = |R_{L1+2}|
!       |0                                          | = |  :     |
!       |0                                          | = |  :     |
!       |0                             A_{L2} B_{L2}| = |R_{L2}  |
!
!-----------------------------------------------------------------------
         use mod_common
         implicit none
         integer(8) :: i, j, l1, l2, ll1, ll2
         real(8) :: gam(n2, n1), a(n2, n1), b(n2, n1), c(n2, n1), r(n2, n1), uu(n2, n1), bet(n2)

         do i = ll1, ll2
           bet(i) = 1./b(i, l1)
           uu(i, l1) = r(i, l1) * bet(i)
         end do

         do j = l1 + 1, l2
           do i = ll1, ll2
             gam(i, j) = c(i, j - 1) * bet(i)
             bet(i) = 1./(b(i, j) - a(i, j) * gam(i, j))
             uu(i, j) = (r(i, j) - a(i, j) * uu(i, j - 1)) * bet(i)
           end do
         end do

         do j = l2 - 1, l1, -1
           do i = ll1, ll2
             uu(i, j) = uu(i, j) - gam(i, j + 1) * uu(i, j + 1)
           end do
         end do

         return
       end

!=======================================================================
       subroutine trdiag2(a, b, c, r, uu, l1, l2, ll1, ll2)
!=======================================================================
!
!     SOLVE THE TRIDIAGONAL MATRIX (X-2 TRIDIAGONAL) WITH
!     DIFFERENT COEFFICIENTS IN LHSU,LHSV,LHSW SUBROUTINES
!
!     VARIABLES
!       A(),B(),C()  : COEFFICIENTS OF TRIDIAGONAL MATRIX
!       R()          : RHS
!       UU()         : SOLUTION
!       L1,L2  : 2ND INDEX INDICATING X2-DIRECTION (TDMA)
!       LL2,LL2: 1ST INDEX INDICATING X1-DIRECTION
!
!-----------------------------------------------------------------------
         use mod_common
         implicit none
         integer(8) :: i, j
         integer(8) :: l1, l2, ll1, ll2
         real(8) :: gam(n1, n2), a(n1, n2), b(n1, n2), c(n1, n2), r(n1, n2), uu(n1, n2), bet(n1)

         do i = ll1, ll2
           bet(i) = 1./b(i, l1)
           uu(i, l1) = r(i, l1) * bet(i)
         end do

         do j = l1 + 1, l2
           do i = ll1, ll2
             gam(i, j) = c(i, j - 1) * bet(i)
             bet(i) = 1./(b(i, j) - a(i, j) * gam(i, j))
             uu(i, j) = (r(i, j) - a(i, j) * uu(i, j - 1)) * bet(i)
           end do
         end do

         do j = l2 - 1, l1, -1
           do i = ll1, ll2
             uu(i, j) = uu(i, j) - gam(i, j + 1) * uu(i, j + 1)
           end do
         end do

         return
       end

!=======================================================================
       subroutine trdiag3(a, b, c, r, uu, l1, l2, ll1, ll2)
!=======================================================================
!
!     SOLVE THE TRIDIAGONAL MATRIX (X-3 TRIDIAGONAL) WITH
!     DIFFERENT COEFFICIENTS IN LHSU,LHSV,LHSW SUBROUTINES
!
!     VARIABLES
!       A(),B(),C()  : COEFFICIENTS OF TRIDIAGONAL MATRIX
!       R()          : RHS
!       UU()         : SOLUTION
!       L1,L2  : 2ND INDEX INDICATING X3-DIRECTION (TDMA)
!       LL2,LL2: 1ST INDEX INDICATING X1-DIRECTION
!
!-----------------------------------------------------------------------
         use mod_common
         implicit none
         integer(8) :: i, j
         integer(8) :: l1, l2, ll1, ll2
         real(8) :: gam(n1, n3), a(n1, n3), b(n1, n3), c(n1, n3), r(n1, n3), uu(n1, n3), bet(n1)

         do i = ll1, ll2
           bet(i) = 1./b(i, l1)
           uu(i, l1) = r(i, l1) * bet(i)
         end do

         do j = l1 + 1, l2
           do i = ll1, ll2
             gam(i, j) = c(i, j - 1) * bet(i)
             bet(i) = 1./(b(i, j) - a(i, j) * gam(i, j))
             uu(i, j) = (r(i, j) - a(i, j) * uu(i, j - 1)) * bet(i)
           end do
         end do

         do j = l2 - 1, l1, -1
           do i = ll1, ll2
             uu(i, j) = uu(i, j) - gam(i, j + 1) * uu(i, j + 1)
           end do
         end do

         return
       end

!=======================================================================
       subroutine trdiag1p(a, b, c, f, j1, j2, l1, l2)
!=======================================================================
!
!     IF XPRDIC = 1 (X PERIODICITY),
!     INVERT A PERIODIC MATRIX (X-1 TRIDIAGONAL) WITH
!     DIFFERENT COEFFICIENTS IN LHSU,LHSV,LHSW SUBROUTINES
!
!     VARIABLES
!       A(),B(),C()  : COEFFICIENTS OF TRIDIAGONAL MATRIX
!       R()          : RHS
!       UU()         : SOLUTION
!       L1,L2  : 2ND INDEX INDICATING X1-DIRECTION (PERIODIC TDMA)
!       LL2,LL2: 1ST INDEX INDICATING X2-DIRECTION
!
!-----------------------------------------------------------------------
         use mod_common
         implicit none
         integer(8) :: i, j, k, j1, j2, l1, l2, ja, jj
         real(8) :: a(n2, n1), b(n2, n1), c(n2, n1), f(n2, n1)
         real(8) :: q(n2, n1), s(n2, n1), qe(n2, n1), fn(n2), pn(n2)
         real(8) :: binv

         ja = j1 + 1
         jj = j1 + j2
         do k = l1, l2
           binv = 1./b(k, j1)
           q(k, j1) = -c(k, j1) * binv
           s(k, j1) = -a(k, j1) * binv
           fn(k) = f(k, j2)
           f(k, j1) = f(k, j1) * binv
         end do

!     FORWARD ELIMINATION SWEEP
         do j = ja, j2
           do k = l1, l2
             pn(k) = 1./(b(k, j) + a(k, j) * q(k, j - 1))
             q(k, j) = -c(k, j) * pn(k)
             s(k, j) = -a(k, j) * s(k, j - 1) * pn(k)
             f(k, j) = (f(k, j) - a(k, j) * f(k, j - 1)) * pn(k)
           end do
         end do

!     BACKWARD PASS
         do k = l1, l2
           s(k, j2) = 1.
           qe(k, j2) = 0.
         end do
         do i = ja, j2
           j = jj - i
           do k = l1, l2
             s(k, j) = s(k, j) + q(k, j) * s(k, j + 1)
             qe(k, j) = f(k, j) + q(k, j) * qe(k, j + 1)
           end do
         end do
         do k = l1, l2
           f(k, j2) = (fn(k) - c(k, j2) * qe(k, j1) - a(k, j2) * qe(k, j2 - 1)) &
                      / (c(k, j2) * s(k, j1) + a(k, j2) * s(k, j2 - 1) + b(k, j2))
         end do

!     BACKWARD ELIMINATION PASS
         do i = ja, j2
           j = jj - i
           do k = l1, l2
             f(k, j) = f(k, j2) * s(k, j) + qe(k, j)
           end do
         end do

         return
       end

!=======================================================================
       subroutine trdiag3p(a, b, c, f, j1, j2, l1, l2)
!=======================================================================
!
!     IF ZPRDIC = 1 (Z PERIODICITY),
!     INVERT A PERIODIC MATRIX (X-3 TRIDIAGONAL) WITH
!     DIFFERENT COEFFICIENTS IN LHSU,LHSV,LHSW SUBROUTINES
!
!     VARIABLES
!       A(),B(),C()  : COEFFICIENTS OF TRIDIAGONAL MATRIX
!       R()          : RHS
!       UU()         : SOLUTION
!       L1,L2  : 2ND INDEX INDICATING X3-DIRECTION (PERIODIC TDMA)
!       LL2,LL2: 1ST INDEX INDICATING X2-DIRECTION
!
!-----------------------------------------------------------------------
         use mod_common
         implicit none
         integer(8) :: i, j, k
         integer(8) :: j1, j2, l1, l2, ja, jj
         real(8) :: a(n1, n3), b(n1, n3), c(n1, n3), f(n1, n3)
         real(8) :: q(n1, n3), s(n1, n3), qe(n1, n3), fn(n1), pn(n1)
         real(8) :: binv

         ja = j1 + 1
         jj = j1 + j2
         do k = l1, l2
           binv = 1./b(k, j1)
           q(k, j1) = -c(k, j1) * binv
           s(k, j1) = -a(k, j1) * binv
           fn(k) = f(k, j2)
           f(k, j1) = f(k, j1) * binv
         end do

!     FORWARD ELIMINATION SWEEP
         do j = ja, j2
           do k = l1, l2
             pn(k) = 1./(b(k, j) + a(k, j) * q(k, j - 1))
             q(k, j) = -c(k, j) * pn(k)
             s(k, j) = -a(k, j) * s(k, j - 1) * pn(k)
             f(k, j) = (f(k, j) - a(k, j) * f(k, j - 1)) * pn(k)
           end do
         end do

!     BACKWARD PASS
         do k = l1, l2
           s(k, j2) = 1.
           qe(k, j2) = 0.
         end do
         do i = ja, j2
           j = jj - i
           do k = l1, l2
             s(k, j) = s(k, j) + q(k, j) * s(k, j + 1)
             qe(k, j) = f(k, j) + q(k, j) * qe(k, j + 1)
           end do
         end do
         do k = l1, l2
           f(k, j2) = (fn(k) - c(k, j2) * qe(k, j1) - a(k, j2) * qe(k, j2 - 1)) &
                      / (c(k, j2) * s(k, j1) + a(k, j2) * s(k, j2 - 1) + b(k, j2))
         end do

!     BACKWARD ELIMINATION PASS
         do i = ja, j2
           j = jj - i
           do k = l1, l2
             f(k, j) = f(k, j2) * s(k, j) + qe(k, j)
           end do
         end do

         return
       end

!=======================================================================
       subroutine ptdiag1a(a, b, c, d, e, f, uu, i1, in, j1, jn)
!=======================================================================
!
!     SOLVE THE PENTADIAGONAL MATRIX (X-1 PENTADIAGONAL)
!
!     VARIABLES
!       A(),B(),C(),D(),E(),F(): COEFFICIENTS OF TRIDIAGONAL MATRIX
!       UU()   : RHS AND SOLUTION
!       I1,IN  : 1ST INDEX INDICATING X2-DIRECTION
!       J1,JN  : 2ND INDEX INDICATING X1-DIRECTION
!
!-----------------------------------------------------------------------
         use mod_common
         implicit none
         integer(8) :: i, j, k
         integer(8) :: i1, in, j1, jn, j2, j3, jnm, jm, jmm
         real(8) :: a(n2, n1), b(n2, n1), c(n2, n1), d(n2, n1), e(n2, n1)
         real(8) :: o(n2, n1), q(n2, n1), r(n2, n1), f(n2, n1), uu(n2, n1)
         real(8) :: pdeno2(n2), pdenoj(n2)

         j2 = j1 + 1
         j3 = j2 + 1
         jnm = jn - 1
         do i = i1, in
           o(i, j1) = -b(i, j1) / a(i, j1)
           q(i, j1) = -c(i, j1) / a(i, j1)
           r(i, j1) = f(i, j1) / a(i, j1)
         end do

         do i = i1, in
           pdeno2(i) = 1./(a(i, j2) + d(i, j2) * o(i, j1))
           o(i, j2) = -(b(i, j2) + d(i, j2) * q(i, j1)) * pdeno2(i)
           q(i, j2) = -c(i, j2) * pdeno2(i)
           r(i, j2) = (f(i, j2) - d(i, j2) * r(i, j1)) * pdeno2(i)
         end do

         do j = j3, jn
           jm = j - 1
           jmm = jm - 1
           do i = i1, in
             pdenoj(i) = 1./(a(i, j) + e(i, j) * q(i, jmm) &
                             + (d(i, j) + e(i, j) * o(i, jmm)) * o(i, jm))
             o(i, j) = -(b(i, j) + (d(i, j) + e(i, j) * o(i, jmm)) * q(i, jm)) * pdenoj(i)
             q(i, j) = -c(i, j) * pdenoj(i)
             r(i, j) = (f(i, j) - e(i, j) * r(i, jmm) - (d(i, j) + e(i, j) * o(i, jmm)) * r(i, jm)) &
                       * pdenoj(i)
           end do
         end do

         do i = i1, in
           uu(i, jn) = r(i, jn)
         end do

         do i = i1, in
           uu(i, jnm) = o(i, jnm) * uu(i, jn) + r(i, jnm)
         end do

         do j = jnm - 1, j1, -1
           do i = i1, in
             uu(i, j) = o(i, j) * uu(i, j + 1) + q(i, j) * uu(i, j + 2) + r(i, j)
           end do
         end do

         return
       end
!=======================================================================
       subroutine retrv_uvw
!=======================================================================
!
!     ADJUST BOUNDARY VELOCITY TO SATISFY GLOBAL MASS CONSERVATION
!           [RETRIEVE UVW (AT THE EXIT BOUNDARY)]
!     SEE 'CONVBC' SUBROUTINE IN LICA_[BODYNAME, E.G. CYLINDER].F90 FILE
!
!-----------------------------------------------------------------------
         use mod_common
         use mod_flowarray, only: u, v, w
         implicit none
         integer(8) :: i, j, k

!$OMP PARALLEL DO
         do k = 1, n3m
           do j = 1, n2m
             u(n1, j, k) = uout(j, k)
           end do
         end do
!$OMP PARALLEL DO
         do k = 1, n3m
           do j = 2, n2m
             v(n1, j, k) = vout(j, k)
           end do
         end do
!$OMP PARALLEL DO
         do k = k_bgpz, n3m
           do j = 1, n2m
             w(n1, j, k) = wout(j, k)
           end do
         end do

         return
       end subroutine retrv_uvw
!=======================================================================
!=======================================================================
       subroutine prdic_adj_uvw(isp)
!=======================================================================
         use mod_common
         use mod_flowarray, only: u, v, w, p
         implicit none
         integer(8) :: isp
         integer(8) :: i, j, k

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
           if (isp .ne. 0) then
!$OMP PARALLEL DO
             do k = 1, n3m
               do j = 1, n2m
                 p(0, j, k) = p(n1m, j, k)
                 p(n1, j, k) = p(1, j, k)
               end do
             end do
           end if
         end if

!     Y PERIODICITY
         if (yprdic .eq. 1) then
!$OMP PARALLEL DO
           do k = 1, n3
             do i = 0, n1
               u(i, 0, k) = u(i, n2m, k)
               u(i, n2, k) = u(i, 1, k)
             end do
           end do
!$OMP PARALLEL DO
           do k = 0, n3
             do i = 0, n1
               v(i, 0, k) = v(i, n2m, k)
               v(i, n2, k) = v(i, 1, k)
             end do
           end do
!$OMP PARALLEL DO
           do k = 0, n3
             do i = 1, n1
               w(i, 0, k) = w(i, n2m, k)
               w(i, n2, k) = w(i, 1, k)
             end do
           end do
           if (isp .ne. 0) then
!$OMP PARALLEL DO
             do k = 1, n3m
               do i = 1, n1m
                 p(i, 0, k) = p(i, n2m, k)
                 p(i, n2, k) = p(i, 1, k)
               end do
             end do
           end if
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
           if (isp .ne. 0) then
!$OMP PARALLEL DO
             do j = 1, n2m
               do i = 1, n1m
                 p(i, j, 0) = p(i, j, n3m)
                 p(i, j, n3) = p(i, j, 1)
               end do
             end do
           end if
         end if

         return
       end subroutine prdic_adj_uvw
!=======================================================================
