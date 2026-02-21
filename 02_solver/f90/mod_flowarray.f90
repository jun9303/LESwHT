module mod_flowarray
!!!!!!!!!!!!!!!!!!!!! BASIC VARIABLES
  real(8), dimension(:, :, :), allocatable :: u, v, w, p, t
  real(8), dimension(:, :, :, :), allocatable :: rhs1
  real(8), dimension(:, :, :), allocatable :: rk3xo, rk3yo, rk3zo, rk3to
!     U,V,W,P           : VELOCITY AND PRESSURE
!     RHS1              : EXPLICITLY TREATED TERMS IN NAVIER-STOKES EQUATION
!     RK3XO,RK3YO,RK3ZO : PREVIOUS STEP INFORMATION FOR RUNGE-KUTTA METHOD (RK3)

!!!!!!!!!!!!!!!!!!!!! IBM VARIABLES
  integer(8), dimension(:, :), allocatable :: ifc, jfc, kfc
  integer(8), dimension(:, :, :), allocatable :: intpindx
  integer(8), dimension(:, :, :, :), allocatable :: inout
  integer(8), dimension(:, :), allocatable :: intptype
  real(8), dimension(:, :), allocatable :: fcv, fcvavg
  real(8), dimension(:, :), allocatable :: dudtr
  real(8), dimension(:, :, :, :, :), allocatable :: geomfac
  real(8), dimension(:, :, :), allocatable :: qmass
  real(8), dimension(:, :, :), allocatable :: divgsum, phi
!     IFC,JFC,KFC : INDEX OF INTERPOLATION POINT FOR IBM
!     MPI         : DISTINGUISH THE SOLID-TO-FLUID DIRECTION
!     INOUT       : 1 => OUTSIDE POINT, 0 => INSIDE POINT
!     JTYPE       : THE TYPE OF INTERPOLATION (1: 1D, 2: 2D, 3: 3D)
!     FCV,FCVAVG  : INSTANTANEOUS AND AVERAGE FORCES OBTAINED FROM THE MOMENTUM FORCING IN IBM
!     DUDTR       : INERTIA CONTRIBUTION IN IBM FORCING
!     GFI         : GEOETRIC FACTOR FOR IBM INTERPOLATION
!     QMASS       : MASS SOURCE AND SINK

!!!!!!!!!!!!!!!!!!!!! LES VARIABLES
  integer(8), dimension(:), allocatable :: inz, jnz, knz, itz, jtz, ktz
  integer(8), dimension(:, :, :), allocatable :: nwall_dvm, nheat_dvm
  real(8), dimension(:, :, :), allocatable :: nusgs, alsgs
  real(8), dimension(:, :, :, :), allocatable :: nusgs1, alsgs1
  real(8), dimension(:, :), allocatable :: cfx1, cfx2, cfy1, cfy2, cfz1, cfz2
  real(8), dimension(:, :, :, :), allocatable :: aalp, llij, mmij, uui
!     INZ,JNZ,KNZ           : INDEX FOR ZERO TURBULENT VISCOSITY POINT
!     NWALL_DVM             : ZERO TURBULENT VISCOSITY IN THE IB BODY
!     NUSGS,NUSGS1          : EDDY VISCOSITY AT EACH VELOCITY CONTROL VOLUME
!     ALSGS,ALSGS1          : EDDY DIFFUSIVITY AT EACH VELOCITY CONTROL VOLUME
!     CFX1,CFX2,CFZP1,CFZP2 : COEFFICIENT FOR LES FILTERING
!     AALP,LLIJ,MMIJ,UUI    : DYNAMIC COEFFICIENT MODELING IN SGS

!!!!!!!!!!!!!!!!!!!!! CONJUGATE HEAT TRANSFER
  real(8), dimension(:, :, :), allocatable :: cstar
  real(8), dimension(:, :, :, :), allocatable :: kstar

!!!!!!!!!!!!!!!!!!!!! SEMI_IMPLICIT
  real(8), dimension(:, :, :), allocatable :: rk3xoo, rk3yoo, rk3zoo, rk3too
!     RK3XOO,RK3YOO,RK3ZOO : CONVECTION TERM FROM INERTIA CONTRIBUTION IN IBM FORCING

!!!!!!!!!!!!!!!!!!!!! FILED AVG. VARIABLES
  real(8), dimension(:, :, :), allocatable :: uavg, vavg, wavg, pavg, tavg
  real(8), dimension(:, :, :), allocatable :: p2avg, t2avg, ssavg
  real(8), dimension(:, :, :, :), allocatable :: uiujavg, voravg, vor2avg

!     UAVG,VAVG,WAVG,PAVG,TAVG : AVERAGED VELOCITY, PRESSURE & TEMPERATURE
!     P2AVG,T2AVG,SSAVG        : AVERAGED PRESSURE & TEMPERATURE FLUCTUATIONS, AVERAGED MAGNITUDE OF STRAIN
!     UIUJAVG,VORAVG,VOR2AVG   : AVERAGED REYNOLDS STRESS, AVERAGED ROOT-MEAN-SQUARE OF VORTICITY

contains

!=======================================================================
  subroutine basic_allo
!=======================================================================
    use mod_common, only: n1, n2, n3, n1m, n2m, n3m, nintp, nbody
    implicit none
    integer(8) :: mintp, mbody

    mintp = maxval(nintp)
    mbody = maxval(nbody)

    allocate (u(0:n1, 0:n2, 0:n3))
    allocate (v(0:n1, 0:n2, 0:n3))
    allocate (w(0:n1, 0:n2, 0:n3))
    allocate (p(0:n1, 0:n2, 0:n3))
    allocate (rhs1(0:n1, 0:n2, 0:n3, 3))
    allocate (rk3xo(n1m, n2m, n3m))
    allocate (rk3yo(n1m, n2m, n3m))
    allocate (rk3zo(n1m, n2m, n3m))
    allocate (rk3xoo(n1m, n2m, n3m))
    allocate (rk3yoo(n1m, n2m, n3m))
    allocate (rk3zoo(n1m, n2m, n3m))
    allocate (inout(0:n1, 0:n2, 0:n3, 3))
    allocate (ifc(mbody, 3), jfc(mbody, 3), kfc(mbody, 3))
    allocate (fcv(mbody, 3), fcvavg(mbody, 3))
    allocate (intptype(mintp, 3))
    allocate (intpindx(mintp, 3, 3))
    allocate (geomfac(mintp, 3, 0:2, 0:2, 0:2))
    allocate (qmass(n1m, n2m, n3m))
    allocate (dudtr(mbody, 3))

    u = 0.
    v = 0.
    w = 0.
    p = 0.
    rhs1 = 0.
    rk3xo = 0.
    rk3yo = 0.
    rk3zo = 0.
    rk3xoo = 0.
    rk3yoo = 0.
    rk3zoo = 0.
    inout = 0
    ifc = 0
    jfc = 0
    kfc = 0
    fcv = 0.
    fcvavg = 0.
    intptype = 0
    intpindx = 0
    geomfac = 0.
    qmass = 0.
    dudtr = 0.

    return
  end subroutine basic_allo
!=======================================================================
  subroutine basic_deallo
!=======================================================================
    implicit none

    deallocate (u, v, w, p)
    deallocate (rhs1, rk3xo, rk3yo, rk3zo)
    deallocate (rk3xoo, rk3yoo, rk3zoo)
    deallocate (inout, ifc, jfc, kfc)
    deallocate (fcv, fcvavg)
    deallocate (intptype, intpindx, geomfac)
    deallocate (qmass, dudtr)

    return
  end subroutine basic_deallo
!=======================================================================
  subroutine thermal_allo
!=======================================================================
    use mod_common, only: n1, n2, n3, n1m, n2m, n3m, nintp, nbody
    implicit none
    integer(8) :: mintp, mbody

    mintp = maxval(nintp)
    mbody = maxval(nbody)

    allocate (u(0:n1, 0:n2, 0:n3))
    allocate (v(0:n1, 0:n2, 0:n3))
    allocate (w(0:n1, 0:n2, 0:n3))
    allocate (p(0:n1, 0:n2, 0:n3))
    allocate (t(0:n1, 0:n2, 0:n3))
    allocate (rhs1(0:n1, 0:n2, 0:n3, 4))
    allocate (rk3xo(n1m, n2m, n3m))
    allocate (rk3yo(n1m, n2m, n3m))
    allocate (rk3zo(n1m, n2m, n3m))
    allocate (rk3to(n1m, n2m, n3m))
    allocate (rk3xoo(n1m, n2m, n3m))
    allocate (rk3yoo(n1m, n2m, n3m))
    allocate (rk3zoo(n1m, n2m, n3m))
    allocate (rk3too(n1m, n2m, n3m))
    allocate (inout(0:n1, 0:n2, 0:n3, 4))
    allocate (ifc(mbody, 4), jfc(mbody, 4), kfc(mbody, 4))
    allocate (fcv(mbody, 4), fcvavg(mbody, 4))
    allocate (intptype(mintp, 4))
    allocate (intpindx(mintp, 4, 3))
    allocate (geomfac(mintp, 4, 0:2, 0:2, 0:2))
    allocate (qmass(n1m, n2m, n3m))
    allocate (dudtr(mbody, 3))

    u = 0.
    v = 0.
    w = 0.
    p = 0.
    t = 0.
    rhs1 = 0.
    rk3xo = 0.
    rk3yo = 0.
    rk3zo = 0.
    rk3xoo = 0.
    rk3yoo = 0.
    rk3zoo = 0.
    inout = 0
    ifc = 0
    jfc = 0
    kfc = 0
    fcv = 0.
    fcvavg = 0.
    intptype = 0
    intpindx = 0
    geomfac = 0.
    qmass = 0.
    dudtr = 0.

    return
  end subroutine thermal_allo
!=======================================================================
  subroutine thermal_deallo
!=======================================================================
    implicit none

    deallocate (u, v, w, p, t)
    deallocate (rhs1, rk3xo, rk3yo, rk3zo)
    deallocate (rk3xoo, rk3yoo, rk3zoo)
    deallocate (inout, ifc, jfc, kfc)
    deallocate (fcv, fcvavg)
    deallocate (intptype, intpindx, geomfac)
    deallocate (qmass, dudtr)

    return
  end subroutine thermal_deallo
!=======================================================================
  subroutine les_allo
!=======================================================================
    use mod_common, only: n1, n2, n3, n1m, n2m, n3m, nzero
    implicit none

    allocate (inz(nzero))
    allocate (jnz(nzero))
    allocate (knz(nzero))
    allocate (nwall_dvm(n1m, n2m, n3m))
    allocate (nusgs(0:n1, 0:n2, 0:n3))
    allocate (nusgs1(0:n1, 0:n2, 0:n3, 3))
    allocate (cfx1(n1m, -1:1))
    allocate (cfx2(n1m, -2:2))
    allocate (cfy1(n2m, -1:1))
    allocate (cfy2(n2m, -2:2))
    allocate (cfz1(n3m, -1:1))
    allocate (cfz2(n3m, -2:2))

    inz = 0
    jnz = 0
    knz = 0
    nwall_dvm = 1
    nusgs = 0.
    nusgs1 = 0.
    cfx1 = 0.
    cfx2 = 0.
    cfy1 = 0.
    cfy2 = 0.
    cfz1 = 0.
    cfz2 = 0.

    return
  end subroutine les_allo
!=======================================================================
  subroutine les_deallo
!=======================================================================
    implicit none

    deallocate (inz, jnz, knz, nwall_dvm)
    deallocate (nusgs, nusgs1)
    deallocate (cfx1, cfx2, cfy1, cfy2, cfz1, cfz2)

    return
  end subroutine les_deallo
!=======================================================================
  subroutine les_thermal_allo
!=======================================================================
    use mod_common, only: n1, n2, n3, n1m, n2m, n3m
    implicit none

    allocate (alsgs(0:n1, 0:n2, 0:n3))
    allocate (alsgs1(0:n1, 0:n2, 0:n3, 3))

    alsgs = 0.
    alsgs1 = 0.

    return
  end subroutine les_thermal_allo
!=======================================================================
  subroutine les_thermal_deallo
!=======================================================================
    use mod_common, only: n1, n2, n3, n1m, n2m, n3m
    implicit none

    deallocate (alsgs, alsgs1)

    return
  end subroutine les_thermal_deallo
!=======================================================================
  subroutine avg_allo
!=======================================================================
    use mod_common, only: n1, n2, n3, n1m, n2m, n3m
    implicit none

    allocate (uavg(n1m, n2m, n3m))
    allocate (vavg(n1m, n2m, n3m))
    allocate (wavg(n1m, n2m, n3m))
    allocate (pavg(n1m, n2m, n3m))
    allocate (tavg(n1m, n2m, n3m))
    allocate (p2avg(n1m, n2m, n3m))
    allocate (t2avg(n1m, n2m, n3m))
    allocate (ssavg(n1m, n2m, n3m))
    allocate (uiujavg(n1m, n2m, n3m, 6))
    allocate (voravg(n1m, n2m, n3m, 3))
    allocate (vor2avg(n1m, n2m, n3m, 6))

    uavg = 0.
    vavg = 0.
    wavg = 0.
    pavg = 0.
    tavg = 0.
    p2avg = 0.
    t2avg = 0.
    ssavg = 0.
    uiujavg = 0.
    voravg = 0.
    vor2avg = 0.

    return
  end subroutine avg_allo
!=======================================================================
  subroutine avg_deallo
!=======================================================================
    implicit none

    deallocate (uavg, vavg, wavg, pavg, tavg)
    deallocate (p2avg, t2avg, ssavg)
    deallocate (uiujavg, voravg, vor2avg)

    return
  end subroutine avg_deallo
!=======================================================================
!=======================================================================
  subroutine conjg_allo
!=======================================================================
    use mod_common, only: n1, n2, n3, n1m, n2m, n3m
    implicit none

    allocate (cstar(n1m, n2m, n3m))
    allocate (kstar(n1m, n2m, n3m, 6))

    cstar = 1.
    kstar = 1.

    return
  end subroutine conjg_allo
!=======================================================================
  subroutine conjg_deallo
!=======================================================================
    implicit none

    deallocate (cstar, kstar)

    return
  end subroutine conjg_deallo
!=======================================================================
end module mod_flowarray
