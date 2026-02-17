MODULE MOD_POISS
!======================= PARAMETERS ==================================
      INTEGER*8, PARAMETER :: MLEV = 20

      REAL*8, DIMENSION(:), ALLOCATABLE :: AW,AE,AS,AN,AB,AF
      REAL*8, DIMENSION(:), ALLOCATABLE :: FIDW,FKDW,FIUP,FKUP
      INTEGER*8 :: IIMG(0:MLEV,3),KKMG(0:MLEV,3)
      INTEGER*8, DIMENSION(:), ALLOCATABLE :: IH1,IH2,KH1,KH2,IL1,IL2,KL1,KL2
      REAL*8, DIMENSION(:), ALLOCATABLE :: COI1,COI2,COR1,COR2
      REAL*8, DIMENSION(:), ALLOCATABLE :: AI3,AK3,AI1,AK1
      INTEGER*8 :: LEVHALF
      INTEGER*8, DIMENSION(:), ALLOCATABLE :: IPM,IMM,KPM,KMM
      REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: AC,GAM,BET

! MLEV: MAXIMUM LEVEL OF MULTIGRID 
! AW,AE,AS,AN,AB,AF: COEEFICIENT FOR THE ITERATIVE METHOD (AX=b)
! FIDW,FKDW, FIUP, FKUP: COEEFICIENT FOR THE MULTIGRID RESTRICTION AND PROLONGATION
! IIMG,KKMG: START AND END INDEX FOR EACH LEVEL OF THE GRID SYSTEM
! IH1,IH2,KH1,KH2: INDEX OF RESIDUE TO BE RESTRICTION
! IL1,IL2,KL1,KL2: INDEX OF RHS TO BE PROLONGATION
! COI1,COI2,COR1,COR2: COEEFICIENT FOR THE MULTIGRID RESTRICTION AND PROLONGATION
! AK3,AI3,AK1,AI1: MODIFIED WAVENUMBER
! IPM,IMM,KPM,KMM: PLUS-MINUS INDEX OF X & Z DIECTION FOR CONSIDERING PERIODICITY
! AC,GAM,BET: COEFICIENT FOR THE SOLVER OF TRIDIGONAL MATRIX

CONTAINS
!=======================================================================
      SUBROUTINE X_FT_ALLO
!=======================================================================
      USE MOD_COMMON
      IMPLICIT NONE

      ALLOCATE(AB(N3MD))
      ALLOCATE(AF(N3MD))
      ALLOCATE(AS(N2M))
      ALLOCATE(AN(N2M))
      ALLOCATE(COI1(N3MD))
      ALLOCATE(COI2(N3MD))
      ALLOCATE(COR1(N3MD))
      ALLOCATE(COR2(N3MD))
      ALLOCATE(AI3(N1MH))
      ALLOCATE(KPM(N3MD),KMM(N3MD))

      AB = 0.
      AF = 0.
      AS = 0.
      AN = 0.
      COI1 = 0.
      COI2 = 0.
      COR1 = 0.
      COR2 = 0.
      AI3 = 0.
      
      RETURN
      END SUBROUTINE X_FT_ALLO
!=======================================================================
!=======================================================================
      SUBROUTINE Z_FT_ALLO
!=======================================================================
      USE MOD_COMMON
      IMPLICIT NONE

      ALLOCATE(AW(N1MD))
      ALLOCATE(AE(N1MD))
      ALLOCATE(AS(N2M))
      ALLOCATE(AN(N2M))
      ALLOCATE(COI1(N1MD))
      ALLOCATE(COI2(N1MD))
      ALLOCATE(COR1(N1MD))
      ALLOCATE(COR2(N1MD))
      ALLOCATE(AK3(N3MH))
      ALLOCATE(IPM(N1MD),IMM(N1MD))

      AW = 0.
      AE = 0.
      AS = 0.
      AN = 0.
      COI1 = 0.
      COI2 = 0.
      COR1 = 0.
      COR2 = 0.
      AK3 = 0.

      RETURN
      END SUBROUTINE Z_FT_ALLO
!=======================================================================
!=======================================================================
      SUBROUTINE MGRD_ALLO
!=======================================================================
      USE MOD_COMMON
      IMPLICIT NONE

      ALLOCATE(AW(N1MD))
      ALLOCATE(AE(N1MD))
      ALLOCATE(AS(N2M))
      ALLOCATE(AN(N2M))
      ALLOCATE(AB(N3MD))
      ALLOCATE(AF(N3MD))
      ALLOCATE(FIDW(N1MD))
      ALLOCATE(FKDW(N3MD))
      ALLOCATE(FIUP(N1MD))
      ALLOCATE(FKUP(N3MD))
      ALLOCATE(IH1(N1MD))
      ALLOCATE(IH2(N1MD))
      ALLOCATE(KH1(N3MD))
      ALLOCATE(KH2(N3MD))
      ALLOCATE(IL1(N1MD))
      ALLOCATE(IL2(N1MD))
      ALLOCATE(KL1(N3MD))
      ALLOCATE(KL2(N3MD))
      ALLOCATE(IPM(N1MD),IMM(N1MD))
      ALLOCATE(KPM(N3MD),KMM(N3MD))

      AW = 0.
      AE = 0.
      AS = 0.
      AN = 0.
      AB = 0.
      AF = 0.
      FIDW = 0.
      FKDW = 0.
      FIUP = 0.
      FKUP = 0.
      IH1 = 0.
      IH2 = 0.
      KH1 = 0.
      KH2 = 0.
      IL1 = 0.
      IL2 = 0.
      KL1 = 0.
      KL2 = 0.

      RETURN
      END SUBROUTINE MGRD_ALLO
!=======================================================================
      SUBROUTINE FTFT_ALLO
!=======================================================================
      USE MOD_COMMON
      IMPLICIT NONE

      ALLOCATE(AI3(N3MH))
      ALLOCATE(AI1(N1))
      ALLOCATE(AK3(N3MH))
      ALLOCATE(AK1(N1))

      AI3 = 0.
      AI1 = 0.
      AK3 = 0.
      AK1 = 0.

      RETURN
      END SUBROUTINE FTFT_ALLO
!=======================================================================
END MODULE MOD_POISS
