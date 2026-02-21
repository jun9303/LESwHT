!=======================================================================
      subroutine divgs(divgsum)
!=======================================================================
!
!     RHS OF POISSON EQUATION
!
!     OPTIONS
!       MASSON = 0, IMOVINGON = 0: W/O IBM BODY
!       MASSON = 1, IMOVINGON = 0: W/  STATIONARY BODY
!       MASSON = 1, IMOVINGON = 1: W/  MOVING BODY
!
!     REFERENCE FOR IBM
!     KIM, KIM & CHOI, 2001, AN IMMERSED-BOUNDARY FINITE-VOLUME METHOD
!       FOR SIMULATIONS OF FLOW IN COMPLEX GEOMETRIES, J. COMP. PHYS.,
!       171, 132-150.
!
!     VARIABLES
!       DIVGSUM: DIVERGENCE OF VELOCITY
!       QMASS: MASS SOURCE/SINK FOR IBM BODY
!       UBD,VBD,WBD: TRANSLATIONAL VELOCITY OF BODY (IMOVINGON=1)
!                    THESE VELOCITIES ARE DEFINED IN LICA_MOVINGIBM.F90.
!
!-----------------------------------------------------------------------
        use mod_common
        use mod_flowarray, only: u, v, w, qmass, ifc, jfc, kfc
        implicit none
        integer(8) :: i, j, k
        integer(8) :: ii, jj, kk, imm, jmm, kmm, n
        real(8) :: divg1, divg2, divg3
        real(8) :: divgsum(0:n1, 0:n2, 0:n3)
        real(8) :: ustmp, vstmp, wstmp
        real(8) :: ubd, vbd, wbd

!############################################################### W/O IBM
!$OMP PARALLEL DO  &
!$OMP PRIVATE(DIVG1,DIVG2,DIVG3)
        do k = 1, n3m
          do j = 1, n2m
            do i = 1, n1m
              divg1 = (u(i + 1, j, k) - u(i, j, k)) * f2fxi(i)
              divg2 = (v(i, j + 1, k) - v(i, j, k)) * f2fyi(j)
              divg3 = (w(i, j, k + 1) - w(i, j, k)) * f2fzi(k)
              divgsum(i, j, k) = (divg1 + divg2 + divg3)
            end do
          end do
        end do

!############################################### W/  STATIONALY IBM BODY
        if (masson .eq. 1 .and. imovingon .eq. 0) then

          qmass = 0.

!$OMP PARALLEL DO  &
!$OMP PRIVATE(II,JJ,KK)
          do n = 1, nbody(1)
            ii = ifc(n, 1)
            jj = jfc(n, 1)
            kk = kfc(n, 1)
            qmass(ii, jj, kk) = qmass(ii, jj, kk) - u(ii, jj, kk) * f2fxi(ii)
            divgsum(ii, jj, kk) = divgsum(ii, jj, kk) + u(ii, jj, kk) * f2fxi(ii)
          end do
!$OMP PARALLEL DO  &
!$OMP PRIVATE(II,JJ,KK,IMM)
          do n = 1, nbody(1)
            ii = ifc(n, 1)
            jj = jfc(n, 1)
            kk = kfc(n, 1)
            imm = imv(ii)
            qmass(imm, jj, kk) = qmass(imm, jj, kk) + u(ii, jj, kk) * f2fxi(imm)
            divgsum(imm, jj, kk) = divgsum(imm, jj, kk) - u(ii, jj, kk) * f2fxi(imm)
          end do

!$OMP PARALLEL DO  &
!$OMP PRIVATE(II,JJ,KK)
          do n = 1, nbody(2)
            ii = ifc(n, 2)
            jj = jfc(n, 2)
            kk = kfc(n, 2)
            qmass(ii, jj, kk) = qmass(ii, jj, kk) - v(ii, jj, kk) * f2fyi(jj)
            divgsum(ii, jj, kk) = divgsum(ii, jj, kk) + v(ii, jj, kk) * f2fyi(jj)
          end do
!$OMP PARALLEL DO  &
!$OMP PRIVATE(II,JJ,KK,JMM)
          do n = 1, nbody(2)
            ii = ifc(n, 2)
            jj = jfc(n, 2)
            kk = kfc(n, 2)
            jmm = jmv(jj)
            qmass(ii, jmm, kk) = qmass(ii, jmm, kk) + v(ii, jj, kk) * f2fyi(jmm)
            divgsum(ii, jmm, kk) = divgsum(ii, jmm, kk) - v(ii, jj, kk) * f2fyi(jmm)
          end do

!$OMP PARALLEL DO  &
!$OMP PRIVATE(II,JJ,KK)
          do n = 1, nbody(3)
            ii = ifc(n, 3)
            jj = jfc(n, 3)
            kk = kfc(n, 3)
            qmass(ii, jj, kk) = qmass(ii, jj, kk) - w(ii, jj, kk) * f2fzi(kk)
            divgsum(ii, jj, kk) = divgsum(ii, jj, kk) + w(ii, jj, kk) * f2fzi(kk)
          end do
!$OMP PARALLEL DO  &
!$OMP PRIVATE(II,JJ,KK,KMM)
          do n = 1, nbody(3)
            ii = ifc(n, 3)
            jj = jfc(n, 3)
            kk = kfc(n, 3)
            kmm = kmv(kk)
            qmass(ii, jj, kmm) = qmass(ii, jj, kmm) + w(ii, jj, kk) * f2fzi(kmm)
            divgsum(ii, jj, kmm) = divgsum(ii, jj, kmm) - w(ii, jj, kk) * f2fzi(kmm)
          end do

! !################################################### W/  MOVING IBM BODY
!       ELSEIF (MASSON.EQ.1 .AND. IMOVINGON.EQ.1) THEN

!       QMASS = 0.

! !$OMP PARALLEL DO  &
! !$OMP PRIVATE(II,JJ,KK,USTMP)
!       DO N=1,NBODY(1)
!         II=IFC(N,1)
!         JJ=JFC(N,1)
!         KK=KFC(N,1)

!         USTMP=U(II,JJ,KK)-UBD(X(II),YMP(JJ),ZMP(KK))

!         QMASS (II,JJ,KK) = QMASS  (II,JJ,KK)-USTMP*F2FXI(II)
!         DIVGSUM(II,JJ,KK)= DIVGSUM(II,JJ,KK)+USTMP*F2FXI(II)
!       END DO
! !$OMP PARALLEL DO  &
! !$OMP PRIVATE(II,JJ,KK,IMM,USTMP)
!       DO N=1,NBODY(1)
!         II=IFC(N,1)
!         JJ=JFC(N,1)
!         KK=KFC(N,1)
!         IMM=IMV(II)

!         USTMP=U(II,JJ,KK)-UBD(X(II),YMP(JJ),ZMP(KK))

!         QMASS  (IMM,JJ,KK)=QMASS  (IMM,JJ,KK)+USTMP*F2FXI(IMM)
!         DIVGSUM(IMM,JJ,KK)=DIVGSUM(IMM,JJ,KK)-USTMP*F2FXI(IMM)
!       END DO

! !$OMP PARALLEL DO  &
! !$OMP PRIVATE(II,JJ,KK,VSTMP)
!       DO N=1,NBODY(2)
!         II=IFC(N,2)
!         JJ=JFC(N,2)
!         KK=KFC(N,2)

!         VSTMP=V(II,JJ,KK)-VBD(XMP(II),Y(JJ),ZMP(KK))

!         QMASS  (II,JJ,KK)=QMASS  (II,JJ,KK)-VSTMP*F2FYI(JJ)
!         DIVGSUM(II,JJ,KK)=DIVGSUM(II,JJ,KK)+VSTMP*F2FYI(JJ)
!       END DO
! !$OMP PARALLEL DO  &
! !$OMP PRIVATE(II,JJ,KK,JMM,VSTMP)
!       DO N=1,NBODY(2)
!         II=IFC(N,2)
!         JJ=JFC(N,2)
!         KK=KFC(N,2)
!         JMM=JMV(JJ)

!         VSTMP=V(II,JJ,KK)-VBD(XMP(II),Y(JJ),ZMP(KK))

!         QMASS  (II,JMM,KK)=QMASS  (II,JMM,KK)+VSTMP*F2FYI(JMM)
!         DIVGSUM(II,JMM,KK)=DIVGSUM(II,JMM,KK)-VSTMP*F2FYI(JMM)
!       END DO

! !$OMP PARALLEL DO  &
! !$OMP PRIVATE(II,JJ,KK,WSTMP)
!       DO N=1,NBODY(3)
!         II=IFC(N,3)
!         JJ=JFC(N,3)
!         KK=KFC(N,3)

!         WSTMP=W(II,JJ,KK)-WBD(XMP(II),YMP(JJ),Z(KK))

!         QMASS  (II,JJ,KK)=QMASS  (II,JJ,KK)-WSTMP*F2FZI(KK)
!         DIVGSUM(II,JJ,KK)=DIVGSUM(II,JJ,KK)+WSTMP*F2FZI(KK)
!       END DO
! !$OMP PARALLEL DO  &
! !$OMP PRIVATE(II,JJ,KK,KMM,WSTMP)
!       DO N=1,NBODY(3)
!         II=IFC(N,3)
!         JJ=JFC(N,3)
!         KK=KFC(N,3)
!         KMM=KMV(KK)

!         WSTMP=W(II,JJ,KK)-WBD(XMP(II),YMP(JJ),Z(KK))

!         QMASS  (II,JJ,KMM)=QMASS  (II,JJ,KMM)+WSTMP*F2FZI(KMM)
!         DIVGSUM(II,JJ,KMM)=DIVGSUM(II,JJ,KMM)-WSTMP*F2FZI(KMM)
!       END DO

        end if

!$OMP PARALLEL DO
        do k = 1, n3m
          do j = 1, n2m
            do i = 1, n1m
              divgsum(i, j, k) = divgsum(i, j, k) * dtconsti
            end do
          end do
        end do

        return
      end
!=======================================================================
