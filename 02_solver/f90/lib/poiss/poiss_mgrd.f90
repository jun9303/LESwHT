!=======================================================================
!
!     MAIN SOLVER OF POISSON EQUATION
!
!     AX=B
!     A: COEFFICIENT OF DISCRETIZED POISSON EQUATION
!     X: PHI: PSEUDO PRESSURE, OUTPUT OF SUBROUTINE POISSON
!     B: DIVGSUM: OUTPUT OF SUBROUTINE DIVGS.
!
!     MAIN ALGORITHM: PRECONDITIONED BICGSTAB
!     [BICONJUGATE GRADIENT STABILIZED METHOD]
!     ALGORITHM: HTTPS://EN.WIKIPEDIA.ORG/WIKI/BICONJUGATE_GRADIENT_STABILIZED_METHOD
!
!     PRECONDITION: GAUSS-SEIDAL WITH MULTIGRID ACCELERATION
!                   X & Z DIRECTION: GAUSS-SEIDAL + MULTIGRID
!                   Y DIRECTION: TDMA
!                   ONLY ONE V-CYCLE FOR THE PRECONDITION
!
!
!     SEPT. 1998, S. KANG: MG3D
!     JUN. 2017,  J. PARK: BICG AND F90
!     FEB. 2026, S. LEE: PERIODIC DOMAIN SUPPORT AND SEVERAL BUG FIXES
!
!======================================================================
      subroutine poisson(phi, divgsum)
!======================================================================
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, k
        integer(8) :: i_cg, imax, jmax, kmax
        real(8) :: rho_cg, rho_ol, alp_cg, w_cg, bet_cg, rv_tmp, acc
        real(8) :: phiref, resm_cg, bttt, uppp

        real(8), dimension(0:n1, 0:n2, 0:n3) :: phi
        real(8), dimension(0:n1, 0:n2, 0:n3) :: divgsum
        real(8), dimension(0:n1, 0:n2, 0:n3) :: p_cg
        real(8), dimension(0:n1, 0:n2, 0:n3) :: s_cg
        real(8), dimension(0:n1, 0:n2, 0:n3) :: y_cg
        real(8), dimension(0:n1, 0:n2, 0:n3) :: z_cg
        real(8), dimension(0:n1, 0:n2, 0:n3) :: v_cg
        real(8), dimension(n1, n2, n3) :: t_cg
        real(8), dimension(n1, n2, n3) :: res_cg

        p_cg = 0.0_8
        s_cg = 0.0_8
        y_cg = 0.0_8
        z_cg = 0.0_8
        v_cg = 0.0_8
        t_cg = 0.0_8

        if (ioldv == 0) then
          phi = 0.0_8
!$OMP PARALLEL DO PRIVATE(I, J)
          do k = 1, n3m
            do j = 1, n2m
              do i = 1, n1m
                res_cg(i, j, k) = divgsum(i, j, k)
              end do
            end do
          end do
        else
!$OMP PARALLEL DO PRIVATE(I, J, ACC)
          do k = 1, n3m
            do j = 1, n2m
              do i = 1, n1m
                acc = -1.0_8 * (aw(i) + ae(i) + as(j) + an(j) + ab(k) + af(k))
                res_cg(i, j, k) = divgsum(i, j, k) + (-acc * phi(i, j, k) &
                                                      - aw(i) * phi(imm(i), j, k) - ae(i) * phi(ipm(i), j, k) &
                                                      - as(j) * phi(i, j - 1, k) - an(j) * phi(i, j + 1, k) &
                                                      - ab(k) * phi(i, j, kmm(k)) - af(k) * phi(i, j, kpm(k)))
              end do
            end do
          end do
        end if

        rho_ol = 1.0_8
        alp_cg = 1.0_8
        w_cg = 1.0_8

        do i_cg = 1, mgitr

          rho_cg = 0.0_8
!$OMP PARALLEL DO PRIVATE(I, J) REDUCTION(+:RHO_CG)
          do k = 1, n3m
            do j = 1, n2m
              do i = 1, n1m
                rho_cg = rho_cg + divgsum(i, j, k) * res_cg(i, j, k)
              end do
            end do
          end do

          bet_cg = rho_cg / rho_ol * alp_cg / w_cg

!$OMP PARALLEL DO PRIVATE(I, J)
          do k = 1, n3m
            do j = 1, n2m
              do i = 1, n1m
                p_cg(i, j, k) = res_cg(i, j, k) + bet_cg * (p_cg(i, j, k) - w_cg * v_cg(i, j, k))
              end do
            end do
          end do

          call mg3d(y_cg, p_cg, test1, ioldv)

!$OMP PARALLEL DO PRIVATE(I, J, ACC)
          do k = 1, n3m
            do j = 1, n2m
              do i = 1, n1m
                acc = -1.0_8 * (aw(i) + ae(i) + as(j) + an(j) + ab(k) + af(k))
                v_cg(i, j, k) = aw(i) * y_cg(imm(i), j, k) + ae(i) * y_cg(ipm(i), j, k) &
                                + as(j) * y_cg(i, j - 1, k) + an(j) * y_cg(i, j + 1, k) &
                                + ab(k) * y_cg(i, j, kmm(k)) + af(k) * y_cg(i, j, kpm(k)) &
                                + acc * y_cg(i, j, k)
              end do
            end do
          end do

          rv_tmp = 0.0_8
!$OMP PARALLEL DO PRIVATE(I, J) REDUCTION(+:RV_TMP)
          do k = 1, n3m
            do j = 1, n2m
              do i = 1, n1m
                rv_tmp = rv_tmp + divgsum(i, j, k) * v_cg(i, j, k)
              end do
            end do
          end do

          alp_cg = rho_cg / rv_tmp

!$OMP PARALLEL DO PRIVATE(I, J)
          do k = 1, n3m
            do j = 1, n2m
              do i = 1, n1m
                s_cg(i, j, k) = res_cg(i, j, k) - alp_cg * v_cg(i, j, k)
              end do
            end do
          end do

          call mg3d(z_cg, s_cg, test1, ioldv)

!$OMP PARALLEL DO PRIVATE(I, J, ACC)
          do k = 1, n3m
            do j = 1, n2m
              do i = 1, n1m
                acc = -1.0_8 * (aw(i) + ae(i) + as(j) + an(j) + ab(k) + af(k))
                t_cg(i, j, k) = aw(i) * z_cg(imm(i), j, k) + ae(i) * z_cg(ipm(i), j, k) &
                                + as(j) * z_cg(i, j - 1, k) + an(j) * z_cg(i, j + 1, k) &
                                + ab(k) * z_cg(i, j, kmm(k)) + af(k) * z_cg(i, j, kpm(k)) &
                                + acc * z_cg(i, j, k)
              end do
            end do
          end do

          bttt = 0.0_8
          uppp = 0.0_8
!$OMP PARALLEL DO PRIVATE(I, J) REDUCTION(+:UPPP) REDUCTION(+:BTTT)
          do k = 1, n3m
            do j = 1, n2m
              do i = 1, n1m
                uppp = uppp + t_cg(i, j, k) * s_cg(i, j, k)
                bttt = bttt + t_cg(i, j, k)**2.0_8
              end do
            end do
          end do
          w_cg = uppp / bttt

!$OMP PARALLEL DO PRIVATE(I, J)
          do k = 1, n3m
            do j = 1, n2m
              do i = 1, n1m
                phi(i, j, k) = phi(i, j, k) + alp_cg * y_cg(i, j, k) + w_cg * z_cg(i, j, k)
              end do
            end do
          end do

          resm_cg = 0.0_8
!$OMP PARALLEL DO PRIVATE(I, J) REDUCTION(MAX:RESM_CG)
          do k = 1, n3m
            do j = 1, n2m
              do i = 1, n1m
                res_cg(i, j, k) = s_cg(i, j, k) - w_cg * t_cg(i, j, k)
                resm_cg = max(resm_cg, abs(res_cg(i, j, k)))
                if (resm_cg == abs(res_cg(i, j, k))) then
                  imax = i
                  jmax = j
                  kmax = k
                end if
              end do
            end do
          end do
          rho_ol = rho_cg

          write (77, 102) i_cg, resm_cg * dtconst, imax, jmax, kmax
102       format(i4, es15.7, 3i4)

          if (resm_cg < test1) then
            write (*, 300) i_cg, resm_cg * dtconst
            goto 1000
          end if

        end do

        print *, '=== MUTLGRID : NOT CONVERGED ==='
        write (*, 300) i_cg, resm_cg * dtconst
300     format('ICYC=', i10, '  RESMAX=', es25.12)

1000    continue

        phiref = phi(1, n2m, 1)
!$OMP PARALLEL DO PRIVATE(I, J)
        do k = 1, n3m
          do j = 1, n2m
            do i = 1, n1m
              phi(i, j, k) = phi(i, j, k) - phiref
            end do
          end do
        end do

        write (78, 103) time, dt, i_cg, resm_cg * dtconst
103     format(2f13.5, i5, es15.7)

        return
      end subroutine poisson

!======================================================================
      subroutine mg3d(sol, rhs, epsil, iold)
!======================================================================
!       CONTAIN OVERALL MULTIGRID PROCEDURES
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, k
        real(8) :: sol(0:n1, 0:n2, 0:n3), rhs(0:n1, 0:n2, 0:n3), epsil
        real(8) :: resid(n1md, n2m, n3md), gi(0:n1md, 0:n2, 0:n3md)
        integer(8) :: iold, icyc, ilev
        real(8) :: resm

        call toplevel(gi, resid, rhs, resm, iold)

        if (nlev == 0) goto 777
        do ilev = nlev - 1, 1, -1
          call relax(gi, resid, ilev, 1_8, 0_8)
          call rescal(gi, resid, ilev)
          call godown(resid, ilev - 1)
        end do

        call relax(gi, resid, 0_8, nbli, 0_8)

        do ilev = 1, nlev - 1
          call goup(gi, ilev)
          call relax(gi, resid, ilev, 1_8, 1_8)
        end do

        call goup(gi, nlev)
777     call toplevel(gi, resid, rhs, resm, 1_8)

!$OMP PARALLEL DO PRIVATE(I, J)
        do k = 1, n3m
          do j = 1, n2m
            do i = 1, n1m
              sol(i, j, k) = gi(i, j, k)
            end do
          end do
        end do

        return
      end subroutine mg3d

!======================================================================
      subroutine toplevel(gi, resid, rhs, resm, iold)   ! ZEBRA VERSION
!======================================================================
!       SOLVE THE POISSON EQUATION WITH ZEBRA GS AT THE TOP LEVEL
!       WITH ONLY 1 ITERATION
!       IOLD = 0 (GODOWN), 1 (GOUP)
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, k
        integer(8) :: istart, iold
        real(8) :: rhs(0:n1, 0:n2, 0:n3)
        real(8) :: resid(n1md, n2m, n3md), gi(0:n1md, 0:n2, 0:n3md)
        real(8) :: acc, resm

!------ MAKE RHS OF ZEBRA GS (ODD LINE)
!$OMP PARALLEL DO PRIVATE(ISTART, I, J)
        do k = 1, n3m
          istart = 1 + mod(kpm(k), 2_8)
          do j = 1, n2m
            do i = istart, n1m, 2
              gi(i, j, k) = rhs(i, j, k) - iold * &
                            (aw(i) * gi(imm(i), j, k) + ae(i) * gi(ipm(i), j, k) &
                             + ab(k) * gi(i, j, kmm(k)) + af(k) * gi(i, j, kpm(k)))
            end do
          end do
        end do

!------ SOLVE TDMA
!$OMP PARALLEL DO PRIVATE(ISTART)
        do k = 1, n3m
          istart = 1 + mod(kpm(k), 2_8)
          do i = istart, n1m, 2
            gi(i, 1, k) = gi(i, 1, k) * bet(i, 1, k)
          end do
        end do

!$OMP PARALLEL DO PRIVATE(ISTART)
        do k = 1, n3m
          istart = 1 + mod(kpm(k), 2_8)
          do j = 2, n2m
            do i = istart, n1m, 2
              gi(i, j, k) = (gi(i, j, k) - as(j) * gi(i, j - 1, k)) * bet(i, j, k)
            end do
          end do
        end do

!$OMP PARALLEL DO PRIVATE(ISTART)
        do k = 1, n3m
          istart = 1 + mod(kpm(k), 2_8)
          do j = n2m - 1, 1, -1
            do i = istart, n1m, 2
              gi(i, j, k) = gi(i, j, k) - gam(i, j + 1, k) * gi(i, j + 1, k)
            end do
          end do
        end do

!------ MAKE RHS OF ZEBRA GS (EVEN LINE)
!$OMP PARALLEL DO PRIVATE(ISTART, I, J)
        do k = 1, n3m
          istart = 1 + mod(k, 2_8)
          do j = 1, n2m
            do i = istart, n1m, 2
              gi(i, j, k) = rhs(i, j, k) &
                            - (aw(i) * gi(imm(i), j, k) + ae(i) * gi(ipm(i), j, k) &
                               + ab(k) * gi(i, j, kmm(k)) + af(k) * gi(i, j, kpm(k)))
            end do
          end do
        end do

!------ SOLVE TDMA
!$OMP PARALLEL DO PRIVATE(ISTART)
        do k = 1, n3m
          istart = 1 + mod(k, 2_8)
          do i = istart, n1m, 2
            gi(i, 1, k) = gi(i, 1, k) * bet(i, 1, k)
          end do
        end do

!$OMP PARALLEL DO PRIVATE(ISTART)
        do k = 1, n3m
          istart = 1 + mod(k, 2_8)
          do j = 2, n2m
            do i = istart, n1m, 2
              gi(i, j, k) = (gi(i, j, k) - as(j) * gi(i, j - 1, k)) * bet(i, j, k)
            end do
          end do
        end do

!$OMP PARALLEL DO PRIVATE(ISTART)
        do k = 1, n3m
          istart = 1 + mod(k, 2_8)
          do j = n2m - 1, 1, -1
            do i = istart, n1m, 2
              gi(i, j, k) = gi(i, j, k) - gam(i, j + 1, k) * gi(i, j + 1, k)
            end do
          end do
        end do

!------ CALCULATE RESIDUAL

        resm = 0.0_8
!$OMP PARALLEL DO PRIVATE(I, J, ACC) REDUCTION(MAX:RESM)
        do k = 1, n3m
          do j = 1, n2m
            do i = 1, n1m
              acc = -1.0_8 * (aw(i) + ae(i) + as(j) + an(j) + ab(k) + af(k))
              resid(i, j, k) = rhs(i, j, k) - acc * gi(i, j, k) &
                               - aw(i) * gi(imm(i), j, k) - ae(i) * gi(ipm(i), j, k) &
                               - as(j) * gi(i, j - 1, k) - an(j) * gi(i, j + 1, k) &
                               - ab(k) * gi(i, j, kmm(k)) - af(k) * gi(i, j, kpm(k))
              resm = max(resm, abs(resid(i, j, k)))
            end do
          end do
        end do

        if (nlev == 0) goto 97

101     format(es15.8, 3i6)

!------ RESTRICT RESIDUAL TO A LOWER LEVEL
!$OMP PARALLEL DO PRIVATE(I, J)
        do k = kkmg(nlev - 1, 1), kkmg(nlev - 1, 2)
          do j = 1, n2m
            do i = iimg(nlev - 1, 1), iimg(nlev - 1, 2)
              resid(i, j, k) = (1.0_8 - fidw(i)) * (1.0_8 - fkdw(k)) * resid(ih1(i), j, kh1(k)) &
                               + fidw(i) * (1.0_8 - fkdw(k)) * resid(ih2(i), j, kh1(k)) &
                               + (1.0_8 - fidw(i)) * fkdw(k) * resid(ih1(i), j, kh2(k)) &
                               + fidw(i) * fkdw(k) * resid(ih2(i), j, kh2(k))
            end do
          end do
        end do

97      return
      end subroutine toplevel

!======================================================================
      subroutine rescal(gi, resid, ilev)
!======================================================================
!      CALCULATE RESIDUAL AT EACH LEVEL
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, k
        integer(8) :: ilev
        real(8) :: resid(n1md, n2m, n3md), gi(0:n1md, 0:n2, 0:n3md)
        real(8) :: acc

!$OMP PARALLEL DO PRIVATE(I, J, ACC)
        do k = kkmg(ilev, 1), kkmg(ilev, 2)
          do j = 1, n2m
            do i = iimg(ilev, 1), iimg(ilev, 2)
              acc = -1.0_8 * (aw(i) + ae(i) + as(j) + an(j) + ab(k) + af(k))
              resid(i, j, k) = resid(i, j, k) - acc * gi(i, j, k) &
                               - aw(i) * gi(imm(i), j, k) - ae(i) * gi(ipm(i), j, k) &
                               - as(j) * gi(i, j - 1, k) - an(j) * gi(i, j + 1, k) &
                               - ab(k) * gi(i, j, kmm(k)) - af(k) * gi(i, j, kpm(k))
            end do
          end do
        end do

        return
      end subroutine rescal

!======================================================================
      subroutine godown(resid, ilev)
!======================================================================
!      RESTRICT RESIDUAL TO A LOWER LEVEL
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, k
        integer(8) :: ilev
        real(8) :: resid(n1md, n2m, n3md)

!$OMP PARALLEL DO PRIVATE(I, J)
        do k = kkmg(ilev, 1), kkmg(ilev, 2)
          do j = 1, n2m
            do i = iimg(ilev, 1), iimg(ilev, 2)
              resid(i, j, k) = (1.0_8 - fidw(i)) * (1.0_8 - fkdw(k)) * resid(ih1(i), j, kh1(k)) &
                               + fidw(i) * (1.0_8 - fkdw(k)) * resid(ih2(i), j, kh1(k)) &
                               + (1.0_8 - fidw(i)) * fkdw(k) * resid(ih1(i), j, kh2(k)) &
                               + fidw(i) * fkdw(k) * resid(ih2(i), j, kh2(k))
            end do
          end do
        end do

        return
      end subroutine godown

!======================================================================
      subroutine goup(gi, ilev)
!======================================================================
!       INTERPOLATE RESIDUAL & ADD IT TO A HIGHER LEVEL
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, k
        integer(8) :: ilev, istart
        real(8) :: gi(0:n1md, 0:n2, 0:n3md)

!$OMP PARALLEL DO PRIVATE(ISTART, I, J)
        do k = kkmg(ilev, 1), kkmg(ilev, 2)
          istart = iimg(ilev, 1) + mod(k, 2_8)
          do j = 1, n2m
            do i = istart, iimg(ilev, 2), 2
              gi(i, j, k) = gi(i, j, k) &
                            + (1.0_8 - fiup(i)) * (1.0_8 - fkup(k)) * gi(il1(i), j, kl1(k)) &
                            + fiup(i) * (1.0_8 - fkup(k)) * gi(il2(i), j, kl1(k)) &
                            + (1.0_8 - fiup(i)) * fkup(k) * gi(il1(i), j, kl2(k)) &
                            + fiup(i) * fkup(k) * gi(il2(i), j, kl2(k))
            end do
          end do
        end do

        return
      end subroutine goup

!======================================================================
      subroutine relax(gi, resid, ilev, iter, iold)   ! ZEBRA VERSION
!======================================================================
!       SOLVE THE POISSON EQUATION WITH ZEBRA GS
!       IOLD = 0 (GODOWN), 1 (GOUP)
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, k
        integer(8) :: ilev, iter, iold, istart, ii
        real(8) :: resid(n1md, n2m, n3md), gi(0:n1md, 0:n2, 0:n3md)

!====== 1ST ITERATION
!------ MAKE RHS OF ZEBRA GS (ODD LINE)
!$OMP PARALLEL DO PRIVATE(ISTART, I, J)
        do k = kkmg(ilev, 1), kkmg(ilev, 2)
          istart = iimg(ilev, 1) + mod(kpm(k), 2_8)
          do j = 1, n2m
            do i = istart, iimg(ilev, 2), 2
              gi(i, j, k) = resid(i, j, k) - iold * &
                            (aw(i) * gi(imm(i), j, k) + ae(i) * gi(ipm(i), j, k) &
                             + ab(k) * gi(i, j, kmm(k)) + af(k) * gi(i, j, kpm(k)))
            end do
          end do
        end do

!------ SOLVE TDMA
!$OMP PARALLEL DO PRIVATE(ISTART)
        do k = kkmg(ilev, 1), kkmg(ilev, 2)
          istart = iimg(ilev, 1) + mod(kpm(k), 2_8)
          do i = istart, iimg(ilev, 2), 2
            gi(i, 1, k) = gi(i, 1, k) * bet(i, 1, k)
          end do
        end do

!$OMP PARALLEL DO PRIVATE(ISTART)
        do k = kkmg(ilev, 1), kkmg(ilev, 2)
          istart = iimg(ilev, 1) + mod(kpm(k), 2_8)
          do j = 2, n2m
            do i = istart, iimg(ilev, 2), 2
              gi(i, j, k) = (gi(i, j, k) - as(j) * gi(i, j - 1, k)) * bet(i, j, k)
            end do
          end do
        end do

!$OMP PARALLEL DO PRIVATE(ISTART)
        do k = kkmg(ilev, 1), kkmg(ilev, 2)
          istart = iimg(ilev, 1) + mod(kpm(k), 2_8)
          do j = n2m - 1, 1, -1
            do i = istart, iimg(ilev, 2), 2
              gi(i, j, k) = gi(i, j, k) - gam(i, j + 1, k) * gi(i, j + 1, k)
            end do
          end do
        end do

!------ MAKE RHS OF ZEBRA GS (EVEN LINE)
!$OMP PARALLEL DO PRIVATE(ISTART, I, J)
        do k = kkmg(ilev, 1), kkmg(ilev, 2)
          istart = iimg(ilev, 1) + mod(k, 2_8)
          do j = 1, n2m
            do i = istart, iimg(ilev, 2), 2
              gi(i, j, k) = resid(i, j, k) &
                            - (aw(i) * gi(imm(i), j, k) + ae(i) * gi(ipm(i), j, k) &
                               + ab(k) * gi(i, j, kmm(k)) + af(k) * gi(i, j, kpm(k)))
            end do
          end do
        end do

!------ SOLVE TDMA
!$OMP PARALLEL DO PRIVATE(ISTART)
        do k = kkmg(ilev, 1), kkmg(ilev, 2)
          istart = iimg(ilev, 1) + mod(k, 2_8)
          do i = istart, iimg(ilev, 2), 2
            gi(i, 1, k) = gi(i, 1, k) * bet(i, 1, k)
          end do
        end do

!$OMP PARALLEL DO PRIVATE(ISTART)
        do k = kkmg(ilev, 1), kkmg(ilev, 2)
          istart = iimg(ilev, 1) + mod(k, 2_8)
          do j = 2, n2m
            do i = istart, iimg(ilev, 2), 2
              gi(i, j, k) = (gi(i, j, k) - as(j) * gi(i, j - 1, k)) * bet(i, j, k)
            end do
          end do
        end do

!$OMP PARALLEL DO PRIVATE(ISTART)
        do k = kkmg(ilev, 1), kkmg(ilev, 2)
          istart = iimg(ilev, 1) + mod(k, 2_8)
          do j = n2m - 1, 1, -1
            do i = istart, iimg(ilev, 2), 2
              gi(i, j, k) = gi(i, j, k) - gam(i, j + 1, k) * gi(i, j + 1, k)
            end do
          end do
        end do

!====== REPEAT PREVIOUS PROCEDURES
        do ii = 1, iter - 1

!------ MAKE RHS OF ZEBRA GS (ODD LINE)
!$OMP PARALLEL DO PRIVATE(ISTART, I, J)
          do k = kkmg(ilev, 1), kkmg(ilev, 2)
            istart = iimg(ilev, 1) + mod(kpm(k), 2_8)
            do j = 1, n2m
              do i = istart, iimg(ilev, 2), 2
                gi(i, j, k) = resid(i, j, k) &
                              - (aw(i) * gi(imm(i), j, k) + ae(i) * gi(ipm(i), j, k) &
                                 + ab(k) * gi(i, j, kmm(k)) + af(k) * gi(i, j, kpm(k)))
              end do
            end do
          end do

!------ SOLVE TDMA
!$OMP PARALLEL DO PRIVATE(ISTART)
          do k = kkmg(ilev, 1), kkmg(ilev, 2)
            istart = iimg(ilev, 1) + mod(kpm(k), 2_8)
            do i = istart, iimg(ilev, 2), 2
              gi(i, 1, k) = gi(i, 1, k) * bet(i, 1, k)
            end do
          end do

!$OMP PARALLEL DO PRIVATE(ISTART)
          do k = kkmg(ilev, 1), kkmg(ilev, 2)
            istart = iimg(ilev, 1) + mod(kpm(k), 2_8)
            do j = 2, n2m
              do i = istart, iimg(ilev, 2), 2
                gi(i, j, k) = (gi(i, j, k) - as(j) * gi(i, j - 1, k)) * bet(i, j, k)
              end do
            end do
          end do

!$OMP PARALLEL DO PRIVATE(ISTART)
          do k = kkmg(ilev, 1), kkmg(ilev, 2)
            istart = iimg(ilev, 1) + mod(kpm(k), 2_8)
            do j = n2m - 1, 1, -1
              do i = istart, iimg(ilev, 2), 2
                gi(i, j, k) = gi(i, j, k) - gam(i, j + 1, k) * gi(i, j + 1, k)
              end do
            end do
          end do

!------ MAKE RHS OF ZEBRA GS (EVEN LINE)
!$OMP PARALLEL DO PRIVATE(ISTART, I, J)
          do k = kkmg(ilev, 1), kkmg(ilev, 2)
            istart = iimg(ilev, 1) + mod(k, 2_8)
            do j = 1, n2m
              do i = istart, iimg(ilev, 2), 2
                gi(i, j, k) = resid(i, j, k) &
                              - (aw(i) * gi(imm(i), j, k) + ae(i) * gi(ipm(i), j, k) &
                                 + ab(k) * gi(i, j, kmm(k)) + af(k) * gi(i, j, kpm(k)))
              end do
            end do
          end do

!------ SOLVE TDMA
!$OMP PARALLEL DO PRIVATE(ISTART)
          do k = kkmg(ilev, 1), kkmg(ilev, 2)
            istart = iimg(ilev, 1) + mod(k, 2_8)
            do i = istart, iimg(ilev, 2), 2
              gi(i, 1, k) = gi(i, 1, k) * bet(i, 1, k)
            end do
          end do

!$OMP PARALLEL DO PRIVATE(ISTART)
          do k = kkmg(ilev, 1), kkmg(ilev, 2)
            istart = iimg(ilev, 1) + mod(k, 2_8)
            do j = 2, n2m
              do i = istart, iimg(ilev, 2), 2
                gi(i, j, k) = (gi(i, j, k) - as(j) * gi(i, j - 1, k)) * bet(i, j, k)
              end do
            end do
          end do

!$OMP PARALLEL DO PRIVATE(ISTART)
          do k = kkmg(ilev, 1), kkmg(ilev, 2)
            istart = iimg(ilev, 1) + mod(k, 2_8)
            do j = n2m - 1, 1, -1
              do i = istart, iimg(ilev, 2), 2
                gi(i, j, k) = gi(i, j, k) - gam(i, j + 1, k) * gi(i, j + 1, k)
              end do
            end do
          end do

        end do

        return
      end subroutine relax

!======================================================================
      subroutine poisinit  ! CALCULATE COEFS FOR POISSON EQ.
!======================================================================
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, k
        integer(8) :: ip, jp, kp

        call mgrd_allo

        allocate (gam(n1md, n2m, n3md))
        allocate (bet(n1md, n2m, n3md))

!-----COMPUTE COEFS OF POISSON EQ.
        do i = 1, n1m
          ip = ipv(i)
          aw(i) = (1.0_8 - fixil(i)) * c2cxi(i) * f2fxi(i)
          ae(i) = (1.0_8 - fixiu(i)) * c2cxi(ip) * f2fxi(i)
        end do

        do j = 1, n2m
          jp = jpv(j)
          as(j) = (1.0_8 - fixjl(j)) * c2cyi(j) * f2fyi(j)
          an(j) = (1.0_8 - fixju(j)) * c2cyi(jp) * f2fyi(j)
        end do

        do k = 1, n3m
          kp = kpv(k)
          ab(k) = (1.0_8 - fixkl(k)) * c2czi(k) * f2fzi(k)
          af(k) = (1.0_8 - fixku(k)) * c2czi(kp) * f2fzi(k)
        end do

        open (77, file='../output/ftr/poiss_itr.dat')
        open (78, file='../output/ftr/ftrpoittr.dat')

        call mgcoef

        return
      end subroutine poisinit

!======================================================================
      subroutine mgcoef  ! CALCULATE COEFS FOR POISSON EQ.
!======================================================================
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, k
        real(8) :: xmpm(0:n1md), zmpm(0:n3md)
        integer(8) :: iwest(n1md), ieast(n1md), kback(n3md), kforw(n3md)
        integer(8) :: ilev
        integer(8) :: ic, ibg, isp, iend, kc, kbg, ksp, kend, ibgh, kbgh, ibgl, kbgl, kendl, iendl
        real(8) :: vdz_kbg, vdz_kend, sdz_kbg, sdz_kend, dz, dz1, dz2
        real(8) :: vdx_ibg, vdx_iend, sdx_ibg, sdx_iend, dx, dx1, dx2

!====== CALCULATE INDICES AT EACH LEVEL
!------ IIMG(),KKMG(),IKMG()
        iimg(nlev, 1) = 1
        iimg(nlev, 2) = n1m
        iimg(nlev, 3) = n1m
        kkmg(nlev, 1) = 1
        kkmg(nlev, 2) = n3m
        kkmg(nlev, 3) = n3m

        do ilev = nlev - 1, 0, -1
          iimg(ilev, 1) = iimg(ilev + 1, 1) + iimg(ilev + 1, 3)
          iimg(ilev, 3) = (n1m / (2**nlev)) * (2**ilev)
          iimg(ilev, 2) = iimg(ilev, 1) + iimg(ilev, 3) - 1

          kkmg(ilev, 1) = kkmg(ilev + 1, 1) + kkmg(ilev + 1, 3)
          kkmg(ilev, 3) = (n3m / (2**nlev)) * (2**ilev)
          kkmg(ilev, 2) = kkmg(ilev, 1) + kkmg(ilev, 3) - 1
        end do

!====== COMPUTE FOR THE FINEST GRID
        do i = 1, n1m
          xmpm(i) = xmp(i)
          iwest(i) = 1 - 1 / i          ! 1 FOR I > 1
          ieast(i) = 1 - i / n1m        ! 1 FOR I < N1M
        end do

        do k = 1, n3m
          zmpm(k) = zmp(k)
          kback(k) = 1 - 1 / k          ! 1 FOR K > 1
          kforw(k) = 1 - k / n3m        ! 1 FOR K < N3M
        end do

!====== COMPUTE FOR COARSE GRIDS
        do ilev = nlev - 1, 0, -1

          ibg = iimg(ilev, 1)
          iend = iimg(ilev, 2)
          isp = 2**(nlev - ilev)

          kbg = kkmg(ilev, 1)
          kend = kkmg(ilev, 2)
          ksp = 2**(nlev - ilev)

          ic = 0
          do i = ibg, iend
            ic = ic + 1
            xmpm(i) = 0.5_8 * (x((ic - 1) * isp + 1) + x(ic * isp + 1))
            iwest(i) = 1 - ibg / i        ! 1 FOR I > IBG
            ieast(i) = 1 - i / iend       ! 1 FOR I < IEND
          end do

          kc = 0
          do k = kbg, kend
            kc = kc + 1
            zmpm(k) = 0.5_8 * (z((kc - 1) * ksp + 1) + z(kc * ksp + 1))
            kback(k) = 1 - kbg / k        ! 1 FOR K > KBG
            kforw(k) = 1 - k / kend       ! 1 FOR K < KEND
          end do

!------ CALCULATE POISSON COEFFICIENTS FOR COARSE GRIDS
          ic = 0
          do i = ibg, iend
            ic = ic + 1
            aw(i) = iwest(i) * 1.0_8 / ((xmpm(i) - xmpm(i - 1)) * (x(ic * isp + 1) - x((ic - 1) * isp + 1)))
            ae(i) = ieast(i) * 1.0_8 / ((xmpm(i + 1) - xmpm(i)) * (x(ic * isp + 1) - x((ic - 1) * isp + 1)))
          end do

          kc = 0
          do k = kbg, kend
            kc = kc + 1
            ab(k) = kback(k) * 1.0_8 / ((zmpm(k) - zmpm(k - 1)) * (z(kc * ksp + 1) - z((kc - 1) * ksp + 1)))
            af(k) = kforw(k) * 1.0_8 / ((zmpm(k + 1) - zmpm(k)) * (z(kc * ksp + 1) - z((kc - 1) * ksp + 1)))
          end do

        end do

!====== CALCULATE RESTRICTION COEFFICIENTS
        do ilev = 0, nlev - 1

          ibgh = iimg(ilev + 1, 1)             ! AT HIGHER LEVEL
          do i = iimg(ilev, 1), iimg(ilev, 2)
            fidw(i) = (xmpm(i) - xmpm(ibgh)) / (xmpm(ibgh + 1) - xmpm(ibgh))
            ih1(i) = ibgh
            ih2(i) = ibgh + 1
            ibgh = ibgh + 2
          end do

          kbgh = kkmg(ilev + 1, 1)             ! AT HIGHER LEVEL
          do k = kkmg(ilev, 1), kkmg(ilev, 2)
            fkdw(k) = (zmpm(k) - zmpm(kbgh)) / (zmpm(kbgh + 1) - zmpm(kbgh))
            kh1(k) = kbgh
            kh2(k) = kbgh + 1
            kbgh = kbgh + 2
          end do

        end do

!====== CALCULATE PROLONGATION COEFFICIENTS
        do ilev = 1, nlev

          ibgl = iimg(ilev - 1, 1)             ! AT LOWER LEVEL
          do i = iimg(ilev, 1), iimg(ilev, 2), 2
            fiup(i) = (xmpm(i) - xmpm(ibgl - 1)) / (xmpm(ibgl) - xmpm(ibgl - 1))
            fiup(i + 1) = (xmpm(i + 1) - xmpm(ibgl)) / (xmpm(ibgl + 1) - xmpm(ibgl))
            il1(i) = ibgl - 1
            il2(i) = ibgl
            il1(i + 1) = ibgl
            il2(i + 1) = ibgl + 1
            ibgl = ibgl + 1
          end do
          i = iimg(ilev, 1)
          fiup(i) = 1.0_8                      ! NEUMANN B.C.
          i = iimg(ilev, 2)
          fiup(i) = 0.0_8                      ! NEUMANN B.C.

          kbgl = kkmg(ilev - 1, 1)             ! AT LOWER LEVEL
          do k = kkmg(ilev, 1), kkmg(ilev, 2), 2
            fkup(k) = (zmpm(k) - zmpm(kbgl - 1)) / (zmpm(kbgl) - zmpm(kbgl - 1))
            fkup(k + 1) = (zmpm(k + 1) - zmpm(kbgl)) / (zmpm(kbgl + 1) - zmpm(kbgl))
            kl1(k) = kbgl - 1
            kl2(k) = kbgl
            kl1(k + 1) = kbgl
            kl2(k + 1) = kbgl + 1
            kbgl = kbgl + 1
          end do
          k = kkmg(ilev, 1)
          fkup(k) = 1.0_8                      ! NEUMANN B.C.
          k = kkmg(ilev, 2)
          fkup(k) = 0.0_8                      ! NEUMANN B.C.

        end do

!       INTRODUCE DEFAULT BOUNDARY INDEX MAPPING
        do ilev = nlev, 0, -1
          kbg = kkmg(ilev, 1)
          kend = kkmg(ilev, 2)
          do k = kbg, kend
            kpm(k) = k + 1
            kmm(k) = k - 1
          end do

          ibg = iimg(ilev, 1)
          iend = iimg(ilev, 2)
          do i = ibg, iend
            ipm(i) = i + 1
            imm(i) = i - 1
          end do
        end do

!       FOR X-PERIODICITY
        if (xprdic == 1) then
          do ilev = nlev, 0, -1
            ibg = iimg(ilev, 1)
            iend = iimg(ilev, 2)
            ipm(iend) = ibg
            imm(ibg) = iend
          end do

          do ilev = nlev - 1, 0, -1
            ibg = iimg(ilev, 1)
            iend = iimg(ilev, 2)
            isp = 2**(nlev - ilev)
            vdx_ibg = xmpm(ibg) - x(1) + x(n1) - xmpm(iend)
            vdx_iend = xmpm(ibg) - x(1) + x(n1) - xmpm(iend)
            sdx_ibg = x(1 + isp) - x(1)
            sdx_iend = x(n1) - x(n1 - isp)
            aw(ibg) = 1.0_8 / (vdx_ibg * sdx_ibg)
            ae(iend) = 1.0_8 / (vdx_iend * sdx_iend)
          end do

          do ilev = 1, nlev
            ibg = iimg(ilev, 1)
            iend = iimg(ilev, 2)
            ibgl = iimg(ilev - 1, 1)
            iendl = iimg(ilev - 1, 2)
            dx = xmpm(ibgl) - x(1) + x(n1) - xmpm(iendl)
            dx1 = xmpm(iend) - xmpm(iendl)
            dx2 = x(n1) - xmpm(iendl) + xmpm(ibg) - x(1)
            fiup(ibg) = dx2 / dx
            fiup(iend) = dx1 / dx
            il1(ibg) = iendl
            il2(ibg) = ibgl
            il1(iend) = iendl
            il2(iend) = ibgl
          end do
        end if

!       FOR Z-PERIODICITY
        if (zprdic == 1) then
          do ilev = nlev, 0, -1
            kbg = kkmg(ilev, 1)
            kend = kkmg(ilev, 2)
            kpm(kend) = kbg
            kmm(kbg) = kend
          end do

          do ilev = nlev - 1, 0, -1
            kbg = kkmg(ilev, 1)
            kend = kkmg(ilev, 2)
            ksp = 2**(nlev - ilev)
            vdz_kbg = zmpm(kbg) - z(1) + z(n3) - zmpm(kend)
            vdz_kend = zmpm(kbg) - z(1) + z(n3) - zmpm(kend)
            sdz_kbg = z(1 + ksp) - z(1)
            sdz_kend = z(n3) - z(n3 - ksp)
            ab(kbg) = 1.0_8 / (vdz_kbg * sdz_kbg)
            af(kend) = 1.0_8 / (vdz_kend * sdz_kend)
          end do

          do ilev = 1, nlev
            kbg = kkmg(ilev, 1)
            kend = kkmg(ilev, 2)
            kbgl = kkmg(ilev - 1, 1)
            kendl = kkmg(ilev - 1, 2)
            dz = zmpm(kbgl) - z(1) + z(n3) - zmpm(kendl)
            dz1 = zmpm(kend) - zmpm(kendl)
            dz2 = z(n3) - zmpm(kendl) + zmpm(kbg) - z(1)
            fkup(kbg) = dz2 / dz
            fkup(kend) = dz1 / dz
            kl1(kbg) = kendl
            kl2(kbg) = kbgl
            kl1(kend) = kendl
            kl2(kend) = kbgl
          end do
        end if

        call tridcoef

        return
      end subroutine mgcoef

!======================================================================
      subroutine tridcoef
!======================================================================
!       COMPUTE TRIDIAGONAL MATRIX COEFFICIENTS
        use mod_common
        use mod_poiss
        implicit none
        integer(8) :: i, j, k
        real(8) :: acc
        integer(8) :: ilev

        do ilev = 0, nlev
          do k = kkmg(ilev, 1), kkmg(ilev, 2)
            do i = iimg(ilev, 1), iimg(ilev, 2)
              acc = -1.0_8 * (aw(i) + ae(i) + as(1) + an(1) + ab(k) + af(k))
              bet(i, 1, k) = 1.0_8 / acc
            end do
          end do

          do k = kkmg(ilev, 1), kkmg(ilev, 2)
            do j = 2, n2m
              do i = iimg(ilev, 1), iimg(ilev, 2)
                acc = -1.0_8 * (aw(i) + ae(i) + as(j) + an(j) + ab(k) + af(k))
                gam(i, j, k) = an(j - 1) * bet(i, j - 1, k)
                bet(i, j, k) = 1.0_8 / (acc - as(j) * gam(i, j, k))
              end do
            end do
          end do
        end do

        return
      end subroutine tridcoef
