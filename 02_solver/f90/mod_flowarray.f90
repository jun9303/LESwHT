MODULE MOD_FLOWARRAY
!!!!!!!!!!!!!!!!!!!!! BASIC VARIABLES 
      REAL*8, DIMENSION (:,:,:),     ALLOCATABLE :: U,V,W,P,T
      REAL*8, DIMENSION (:,:,:,:),   ALLOCATABLE :: RHS1
      REAL*8, DIMENSION (:,:,:),     ALLOCATABLE :: RK3XO,RK3YO,RK3ZO,RK3TO
!     U,V,W,P           : Velocity and Pressure
!     RHS1              : Explicitly treated terms in Navier-Stokes equation
!     RK3XO,RK3YO,RK3ZO : Previous step information for Runge-Kutta method (RK3)

!!!!!!!!!!!!!!!!!!!!! IBM VARIABLES
      INTEGER*8, DIMENSION (:,:),    ALLOCATABLE :: IFC,JFC,KFC
      INTEGER*8, DIMENSION (:,:,:),  ALLOCATABLE :: INTPINDX
      INTEGER*8, DIMENSION (:,:,:,:),ALLOCATABLE :: INOUT
      INTEGER*8, DIMENSION (:,:),    ALLOCATABLE :: INTPTYPE
      REAL*8, DIMENSION (:,:),       ALLOCATABLE :: FCV,FCVAVG
      REAL*8, DIMENSION (:,:),       ALLOCATABLE :: DUDTR
      REAL*8, DIMENSION (:,:,:,:,:), ALLOCATABLE :: GEOMFAC
      REAL*8, DIMENSION (:,:,:),     ALLOCATABLE :: QMASS
      REAL*8, DIMENSION (:,:,:),     ALLOCATABLE :: DIVGSUM, PHI
!     IFC,JFC,KFC : Index of interpolation point for IBM
!     MPI         : Distinguish the solid-to-fluid direction
!     INOUT       : 1 => Outside point, 0 => Inside point
!     JTYPE       : The type of interpolation (1: 1D, 2: 2D, 3: 3D)
!     FCV,FCVAVG  : Instantaneous and average forces obtained from the momentum forcing in IBM
!     DUDTR       : Inertia contribution in IBM forcing
!     GFI         : Geoetric factor for IBM interpolation
!     QMASS       : Mass source and sink

!!!!!!!!!!!!!!!!!!!!! LES VARIABLES
      INTEGER*8, DIMENSION (:),      ALLOCATABLE :: INZ,JNZ,KNZ,ITZ,JTZ,KTZ
      INTEGER*8, DIMENSION (:,:,:),  ALLOCATABLE :: NWALL_DVM,NHEAT_DVM
      REAL*8, DIMENSION (:,:,:),     ALLOCATABLE :: NUSGS,ALSGS
      REAL*8, DIMENSION (:,:,:,:),   ALLOCATABLE :: NUSGS1,ALSGS1
      REAL*8, DIMENSION (:,:),       ALLOCATABLE :: CFX1,CFX2,CFY1,CFY2,CFZ1,CFZ2
      REAL*8, DIMENSION (:,:,:,:),   ALLOCATABLE :: AALP,LLIJ,MMIJ,UUI
!     INZ,JNZ,KNZ           : Index for zero turbulent viscosity point
!     NWALL_DVM             : Zero turbulent viscosity in the IB body
!     NUSGS,NUSGS1          : Eddy viscosity at each velocity control volume
!     ALSGS,ALSGS1          : Eddy diffusivity at each velocity control volume
!     CFX1,CFX2,CFZP1,CFZP2 : Coefficient for LES filtering
!     AALP,LLIJ,MMIJ,UUI    : Dynamic coefficient modeling in SGS

!!!!!!!!!!!!!!!!!!!!! CONJUGATE HEAT TRANSFER
      REAL*8, DIMENSION (:,:,:),     ALLOCATABLE :: CSTAR
      REAL*8, DIMENSION (:,:,:,:),   ALLOCATABLE :: KSTAR

!!!!!!!!!!!!!!!!!!!!! SEMI_IMPLICIT
      REAL*8, DIMENSION (:,:,:),     ALLOCATABLE :: RK3XOO,RK3YOO,RK3ZOO,RK3TOO
!     RK3XOO,RK3YOO,RK3ZOO : Convection term from inertia contribution in IBM forcing

!!!!!!!!!!!!!!!!!!!!! FILED AVG. VARIABLES 
      REAL*8, DIMENSION (:,:,:),     ALLOCATABLE :: UAVG,VAVG,WAVG,PAVG,TAVG
      REAL*8, DIMENSION (:,:,:),     ALLOCATABLE :: P2AVG,T2AVG,SSAVG
      REAL*8, DIMENSION (:,:,:,:),   ALLOCATABLE :: UIUJAVG,VORAVG,VOR2AVG

!     UAVG,VAVG,WAVG,PAVG,TAVG : AVERAGED VELOCITY, PRESSURE & TEMPERATURE
!     P2AVG,T2AVG,SSAVG        : AVERAGED PRESSURE & TEMPERATURE FLUCTUATIONS, AVERAGED MAGNITUDE OF STRAIN
!     UIUJAVG,VORAVG,VOR2AVG   : AVERAGED REYNOLDS STRESS, AVERAGED ROOT-MEAN-SQUARE OF VORTICITY

CONTAINS

!=======================================================================
      SUBROUTINE BASIC_ALLO
!=======================================================================
      USE MOD_COMMON, ONLY : N1,N2,N3,N1M,N2M,N3M,NINTP,NBODY
      IMPLICIT NONE
      INTEGER*8         :: MINTP,MBODY

      MINTP = MAXVAL(NINTP)
      MBODY = MAXVAL(NBODY)

      ALLOCATE( U(0:N1,0:N2,0:N3))
      ALLOCATE( V(0:N1,0:N2,0:N3))
      ALLOCATE( W(0:N1,0:N2,0:N3))
      ALLOCATE( P(0:N1,0:N2,0:N3))
      ALLOCATE( RHS1(0:N1,0:N2,0:N3,3))
      ALLOCATE( RK3XO(N1M,N2M,N3M))
      ALLOCATE( RK3YO(N1M,N2M,N3M))
      ALLOCATE( RK3ZO(N1M,N2M,N3M))
      ALLOCATE(RK3XOO(N1M,N2M,N3M))
      ALLOCATE(RK3YOO(N1M,N2M,N3M))
      ALLOCATE(RK3ZOO(N1M,N2M,N3M))
      ALLOCATE(INOUT(0:N1,0:N2,0:N3,3))
      ALLOCATE(IFC(MBODY,3),JFC(MBODY,3),KFC(MBODY,3))
      ALLOCATE(FCV(MBODY,3),FCVAVG(MBODY,3))
      ALLOCATE(INTPTYPE(MINTP,3))
      ALLOCATE(INTPINDX(MINTP,3,3))
      ALLOCATE(GEOMFAC(MINTP,3,0:2,0:2,0:2))
      ALLOCATE(QMASS(N1M,N2M,N3M))
      ALLOCATE(DUDTR(MBODY,3))

      U = 0.
      V = 0.
      W = 0.
      P = 0.
      RHS1 = 0.
      RK3XO = 0.
      RK3YO = 0.
      RK3ZO = 0.
      RK3XOO = 0.
      RK3YOO = 0.
      RK3ZOO = 0.
      INOUT = 0
      IFC = 0
      JFC = 0
      KFC = 0
      FCV = 0.
      FCVAVG = 0.
      INTPTYPE = 0
      INTPINDX = 0
      GEOMFAC = 0.
      QMASS = 0.
      DUDTR = 0.

      RETURN
      END SUBROUTINE BASIC_ALLO
!=======================================================================
      SUBROUTINE BASIC_DEALLO
!=======================================================================
      IMPLICIT NONE

      DEALLOCATE(U,V,W,P)
      DEALLOCATE(RHS1,RK3XO,RK3YO,RK3ZO)
      DEALLOCATE(RK3XOO,RK3YOO,RK3ZOO)
      DEALLOCATE(INOUT,IFC,JFC,KFC)
      DEALLOCATE(FCV,FCVAVG)
      DEALLOCATE(INTPTYPE,INTPINDX,GEOMFAC)
      DEALLOCATE(QMASS,DUDTR)

      RETURN
      END SUBROUTINE BASIC_DEALLO
!=======================================================================
      SUBROUTINE THERMAL_ALLO
!=======================================================================
      USE MOD_COMMON, ONLY : N1,N2,N3,N1M,N2M,N3M,NINTP,NBODY
      IMPLICIT NONE
      INTEGER*8         :: MINTP,MBODY

      MINTP = MAXVAL(NINTP)
      MBODY = MAXVAL(NBODY)
      
      ALLOCATE( U(0:N1,0:N2,0:N3))
      ALLOCATE( V(0:N1,0:N2,0:N3))
      ALLOCATE( W(0:N1,0:N2,0:N3))
      ALLOCATE( P(0:N1,0:N2,0:N3))
      ALLOCATE( T(0:N1,0:N2,0:N3))
      ALLOCATE( RHS1(0:N1,0:N2,0:N3,4))
      ALLOCATE( RK3XO(N1M,N2M,N3M))
      ALLOCATE( RK3YO(N1M,N2M,N3M))
      ALLOCATE( RK3ZO(N1M,N2M,N3M))
      ALLOCATE( RK3TO(N1M,N2M,N3M))
      ALLOCATE(RK3XOO(N1M,N2M,N3M))
      ALLOCATE(RK3YOO(N1M,N2M,N3M))
      ALLOCATE(RK3ZOO(N1M,N2M,N3M))
      ALLOCATE(RK3TOO(N1M,N2M,N3M))
      ALLOCATE(INOUT(0:N1,0:N2,0:N3,4))
      ALLOCATE(IFC(MBODY,4),JFC(MBODY,4),KFC(MBODY,4))
      ALLOCATE(FCV(MBODY,4),FCVAVG(MBODY,4))
      ALLOCATE(INTPTYPE(MINTP,4))
      ALLOCATE(INTPINDX(MINTP,4,3))
      ALLOCATE(GEOMFAC(MINTP,4,0:2,0:2,0:2))
      ALLOCATE(QMASS(N1M,N2M,N3M))
      ALLOCATE(DUDTR(MBODY,3))

      U = 0.
      V = 0.
      W = 0.
      P = 0.
      T = 0.
      RHS1 = 0.
      RK3XO = 0.
      RK3YO = 0.
      RK3ZO = 0.
      RK3XOO = 0.
      RK3YOO = 0.
      RK3ZOO = 0.
      INOUT = 0
      IFC = 0
      JFC = 0
      KFC = 0
      FCV = 0.
      FCVAVG = 0.
      INTPTYPE = 0
      INTPINDX = 0
      GEOMFAC = 0.
      QMASS = 0.
      DUDTR = 0.

      RETURN
      END SUBROUTINE THERMAL_ALLO
!=======================================================================
      SUBROUTINE THERMAL_DEALLO
!=======================================================================
      IMPLICIT NONE
      
      DEALLOCATE(U,V,W,P,T)
      DEALLOCATE(RHS1,RK3XO,RK3YO,RK3ZO)
      DEALLOCATE(RK3XOO,RK3YOO,RK3ZOO)
      DEALLOCATE(INOUT,IFC,JFC,KFC)
      DEALLOCATE(FCV,FCVAVG)
      DEALLOCATE(INTPTYPE,INTPINDX,GEOMFAC)
      DEALLOCATE(QMASS,DUDTR)

      RETURN
      END SUBROUTINE THERMAL_DEALLO
!=======================================================================
      SUBROUTINE LES_ALLO
!=======================================================================
      USE MOD_COMMON, ONLY : N1,N2,N3,N1M,N2M,N3M,NZERO
      IMPLICIT NONE

      ALLOCATE(INZ(NZERO))
      ALLOCATE(JNZ(NZERO))
      ALLOCATE(KNZ(NZERO))
      ALLOCATE(NWALL_DVM(N1M,N2M,N3M))
      ALLOCATE(NUSGS (0:N1,0:N2,0:N3))
      ALLOCATE(NUSGS1(0:N1,0:N2,0:N3,3))
      ALLOCATE(CFX1(N1M,-1:1))
      ALLOCATE(CFX2(N1M,-2:2))
      ALLOCATE(CFY1(N2M,-1:1))
      ALLOCATE(CFY2(N2M,-2:2))
      ALLOCATE(CFZ1(N3M,-1:1))
      ALLOCATE(CFZ2(N3M,-2:2))

      INZ = 0
      JNZ = 0
      KNZ = 0
      NWALL_DVM = 1
      NUSGS = 0.
      NUSGS1 = 0.
      CFX1 = 0.
      CFX2 = 0.
      CFY1 = 0.
      CFY2 = 0.
      CFZ1 = 0.
      CFZ2 = 0.

      RETURN
      END SUBROUTINE LES_ALLO
!=======================================================================
      SUBROUTINE LES_DEALLO
!=======================================================================
      IMPLICIT NONE

      DEALLOCATE(INZ,JNZ,KNZ,NWALL_DVM)
      DEALLOCATE(NUSGS,NUSGS1)
      DEALLOCATE(CFX1,CFX2,CFY1,CFY2,CFZ1,CFZ2)

      RETURN
      END SUBROUTINE LES_DEALLO
!=======================================================================
      SUBROUTINE LES_THERMAL_ALLO
!=======================================================================
      USE MOD_COMMON, ONLY : N1,N2,N3,N1M,N2M,N3M
      IMPLICIT NONE

      ALLOCATE(ALSGS (0:N1,0:N2,0:N3))
      ALLOCATE(ALSGS1(0:N1,0:N2,0:N3,3))

      ALSGS = 0.
      ALSGS1 = 0.

      RETURN
      END SUBROUTINE LES_THERMAL_ALLO
!=======================================================================
      SUBROUTINE LES_THERMAL_DEALLO
!=======================================================================
      USE MOD_COMMON, ONLY : N1,N2,N3,N1M,N2M,N3M
      IMPLICIT NONE

      DEALLOCATE(ALSGS,ALSGS1)

      RETURN
      END SUBROUTINE LES_THERMAL_DEALLO
!=======================================================================
      SUBROUTINE AVG_ALLO
!=======================================================================
      USE MOD_COMMON, ONLY : N1,N2,N3,N1M,N2M,N3M
      IMPLICIT NONE

      ALLOCATE(UAVG(N1M,N2M,N3M))
      ALLOCATE(VAVG(N1M,N2M,N3M))
      ALLOCATE(WAVG(N1M,N2M,N3M))
      ALLOCATE(PAVG(N1M,N2M,N3M))
      ALLOCATE(TAVG(N1M,N2M,N3M))
      ALLOCATE(P2AVG(N1M,N2M,N3M))
      ALLOCATE(T2AVG(N1M,N2M,N3M))
      ALLOCATE(SSAVG(N1M,N2M,N3M))
      ALLOCATE(UIUJAVG(N1M,N2M,N3M,6))
      ALLOCATE(VORAVG(N1M,N2M,N3M,3))
      ALLOCATE(VOR2AVG(N1M,N2M,N3M,6))

      UAVG = 0.
      VAVG = 0.
      WAVG = 0.
      PAVG = 0.
      TAVG = 0.
      P2AVG = 0.
      T2AVG = 0.
      SSAVG = 0.
      UIUJAVG = 0.
      VORAVG = 0.
      VOR2AVG = 0.

      RETURN
      END SUBROUTINE AVG_ALLO
!=======================================================================
      SUBROUTINE AVG_DEALLO
!=======================================================================
      IMPLICIT NONE

      DEALLOCATE(UAVG,VAVG,WAVG,PAVG,TAVG)
      DEALLOCATE(P2AVG,T2AVG,SSAVG)
      DEALLOCATE(UIUJAVG,VORAVG,VOR2AVG)

      RETURN
      END SUBROUTINE AVG_DEALLO
!=======================================================================
!=======================================================================
      SUBROUTINE CONJG_ALLO
!=======================================================================
      USE MOD_COMMON, ONLY : N1,N2,N3,N1M,N2M,N3M
      IMPLICIT NONE

      ALLOCATE(CSTAR(N1M,N2M,N3M))
      ALLOCATE(KSTAR(N1M,N2M,N3M,6))

      CSTAR = 1.
      KSTAR = 1.

      RETURN
      END SUBROUTINE CONJG_ALLO
!=======================================================================
      SUBROUTINE CONJG_DEALLO
!=======================================================================
      IMPLICIT NONE

      DEALLOCATE(CSTAR,KSTAR)

      RETURN
      END SUBROUTINE CONJG_DEALLO
!=======================================================================
END MODULE MOD_FLOWARRAY