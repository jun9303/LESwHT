module mod_poiss
!======================= PARAMETERS ==================================
  integer(8), parameter :: mlev = 20

  real(8), dimension(:), allocatable :: aw, ae, as, an, ab, af
  real(8), dimension(:), allocatable :: fidw, fkdw, fiup, fkup
  integer(8) :: iimg(0:mlev, 3), kkmg(0:mlev, 3)
  integer(8), dimension(:), allocatable :: ih1, ih2, kh1, kh2, il1, il2, kl1, kl2
  real(8), dimension(:), allocatable :: coi1, coi2, cor1, cor2
  real(8), dimension(:), allocatable :: ai3, ak3, ai1, ak1
  integer(8) :: levhalf
  integer(8), dimension(:), allocatable :: ipm, imm, kpm, kmm
  real(8), dimension(:, :, :), allocatable :: ac, gam, bet

! MLEV: MAXIMUM LEVEL OF MULTIGRID
! AW,AE,AS,AN,AB,AF: COEEFICIENT FOR THE ITERATIVE METHOD (AX=B)
! FIDW,FKDW, FIUP, FKUP: COEEFICIENT FOR THE MULTIGRID RESTRICTION AND PROLONGATION
! IIMG,KKMG: START AND END INDEX FOR EACH LEVEL OF THE GRID SYSTEM
! IH1,IH2,KH1,KH2: INDEX OF RESIDUE TO BE RESTRICTION
! IL1,IL2,KL1,KL2: INDEX OF RHS TO BE PROLONGATION
! COI1,COI2,COR1,COR2: COEEFICIENT FOR THE MULTIGRID RESTRICTION AND PROLONGATION
! AK3,AI3,AK1,AI1: MODIFIED WAVENUMBER
! IPM,IMM,KPM,KMM: PLUS-MINUS INDEX OF X & Z DIECTION FOR CONSIDERING PERIODICITY
! AC,GAM,BET: COEFICIENT FOR THE SOLVER OF TRIDIGONAL MATRIX

contains
!=======================================================================
  subroutine x_ft_allo
!=======================================================================
    use mod_common
    implicit none

    allocate (ab(n3md))
    allocate (af(n3md))
    allocate (as(n2m))
    allocate (an(n2m))
    allocate (coi1(n3md))
    allocate (coi2(n3md))
    allocate (cor1(n3md))
    allocate (cor2(n3md))
    allocate (ai3(n1mh))
    allocate (kpm(n3md), kmm(n3md))

    ab = 0.
    af = 0.
    as = 0.
    an = 0.
    coi1 = 0.
    coi2 = 0.
    cor1 = 0.
    cor2 = 0.
    ai3 = 0.

    return
  end subroutine x_ft_allo
!=======================================================================
!=======================================================================
  subroutine z_ft_allo
!=======================================================================
    use mod_common
    implicit none

    allocate (aw(n1md))
    allocate (ae(n1md))
    allocate (as(n2m))
    allocate (an(n2m))
    allocate (coi1(n1md))
    allocate (coi2(n1md))
    allocate (cor1(n1md))
    allocate (cor2(n1md))
    allocate (ak3(n3mh))
    allocate (ipm(n1md), imm(n1md))

    aw = 0.
    ae = 0.
    as = 0.
    an = 0.
    coi1 = 0.
    coi2 = 0.
    cor1 = 0.
    cor2 = 0.
    ak3 = 0.

    return
  end subroutine z_ft_allo
!=======================================================================
!=======================================================================
  subroutine mgrd_allo
!=======================================================================
    use mod_common
    implicit none

    allocate (aw(n1md))
    allocate (ae(n1md))
    allocate (as(n2m))
    allocate (an(n2m))
    allocate (ab(n3md))
    allocate (af(n3md))
    allocate (fidw(n1md))
    allocate (fkdw(n3md))
    allocate (fiup(n1md))
    allocate (fkup(n3md))
    allocate (ih1(n1md))
    allocate (ih2(n1md))
    allocate (kh1(n3md))
    allocate (kh2(n3md))
    allocate (il1(n1md))
    allocate (il2(n1md))
    allocate (kl1(n3md))
    allocate (kl2(n3md))
    allocate (ipm(n1md), imm(n1md))
    allocate (kpm(n3md), kmm(n3md))

    aw = 0.
    ae = 0.
    as = 0.
    an = 0.
    ab = 0.
    af = 0.
    fidw = 0.
    fkdw = 0.
    fiup = 0.
    fkup = 0.
    ih1 = 0.
    ih2 = 0.
    kh1 = 0.
    kh2 = 0.
    il1 = 0.
    il2 = 0.
    kl1 = 0.
    kl2 = 0.

    return
  end subroutine mgrd_allo
!=======================================================================
  subroutine ftft_allo
!=======================================================================
    use mod_common
    implicit none

    allocate (ai3(n3mh))
    allocate (ai1(n1))
    allocate (ak3(n3mh))
    allocate (ak1(n1))

    ai3 = 0.
    ai1 = 0.
    ak3 = 0.
    ak1 = 0.

    return
  end subroutine ftft_allo
!=======================================================================
end module mod_poiss
