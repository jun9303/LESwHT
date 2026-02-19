MODULE MOD_COMMON
!!!!!!!!!!!!!!!!!!!!! BASIC INPUTS (SETTINGS.INPUT)
      REAL*8                 :: RE, PR, GR, EPS_PTR, UBULK_I, DT_SIZE,      &
                                CFLFAC, RESID1, WWSOR, CSGSTS, CSGSHF
      INTEGER*8              :: IRESET, IREAD, IAVG, IPZERO, NTST, NPRINT,  &
                                NPRIAVG, NPIN, IDTOPT, NLEV, MODE, NBLI,    &
                                IOLDV, MGITR, IMGSOR, ILES, INSMDL, ITEMDL, &
                                IDVMON, FILTER, IBMON, MASSON, GRDIR, T_INF,&
                                IMOVINGON, IHTRANS, NTRACE, NTR, TRPTS(64,3)
      CHARACTER*25           :: GRIDFILE
      CHARACTER*25           :: PREV_FLD

!!!!!!!!!!!!!!!!!!!!! BOUNDARY CONDITIONS (BOUNDARY.INPUT & IBMPRE_PRDIC.BIN)
      INTEGER*8              :: XPRDIC, YPRDIC, ZPRDIC, IINTP
      INTEGER*8              :: BC_XBTM, BC_XTOP, BC_YBTM, BC_YTOP,        &
                                BC_ZBTM, BC_ZTOP
      INTEGER*8              :: BC_T_YBTM, BC_T_YTOP, BC_T_ZBTM, BC_T_ZTOP
      REAL*8                 :: VAL_T_YBTM, VAL_T_YTOP, VAL_T_ZBTM, VAL_T_ZTOP
      INTEGER*8              :: I_BGPX, J_BGPY, K_BGPZ
      REAL*8, DIMENSION(:,:), ALLOCATABLE  :: DUOUT,DVOUT,DWOUT,DTOUT
      REAL*8, DIMENSION(:,:), ALLOCATABLE  :: UOUT,VOUT,WOUT,TOUT
      INTEGER*8              :: JUT,JUB,JWT,JWB,KUT,KUB,KVT,KVB

!!!!!!!!!!!!!!!!!!!!! GRID & MESH & DOMAIN (GRID.DAT)
      INTEGER*8              :: N1,N1M,N1MD,N1MH,N2,N2M,N2MD,N2MH,N3,N3M,N3MD,N3MH
      REAL*8                 :: XL,YL,ZL
      REAL*8, DIMENSION(:), ALLOCATABLE    :: X,Y,Z,XMP,YMP,ZMP
      INTEGER*8, DIMENSION(:), ALLOCATABLE :: IPV,JPV,KPV,IMV,JMV,KMV
      INTEGER*8, DIMENSION(:), ALLOCATABLE :: FIXIL,FIXIU,FIXJL,FIXJU,FIXKL,FIXKU
      REAL*8, DIMENSION(:), ALLOCATABLE    :: F2FX,F2FY,F2FZ,F2FXI,F2FYI,F2FZI
      REAL*8, DIMENSION(:), ALLOCATABLE    :: C2CX,C2CY,C2CZ,C2CXI,C2CYI,C2CZI
!     N_,N_M,N_MD,N_MH                    : THE NUMBER OF GRID POINTS FOR X,Y,Z DIRECTIONS
!     XL,YL,ZL                            : DOMAIN SIZES
!     X,Y,Z                               : POSITIONS OF GRID POINTS
!     IPV,JPV,KPV,IMV,JMV,KMV             : NEXT AND PREVIOUS GRID INDEX
!     FIXIL,FIXIU,FIXJL,FIXJU,FIXKL,FIXKU : TREAT DOMAIN BOUNDARY
!     F2FX,F2FY,F2FZ                      : GRID SIZES FROM CELL FACE TO FACE
!     F2FXI,F2FYI,F2FZI                   : INVERSE OF GRID SIZES (SDX,SDY,SDZ)
!     C2CX,C2CY,C2CZ                      : GRID SIZES FROM CELL CENTER TO CENTER
!     C2CXI,C2CYI,C2CZI                   : INVERSE OF GRID SIZES (VDX,VDY,VDZ)
!     XMP,YMP,ZMP                         : POSITION OF CELL CENTER POINTS

!!!!!!!!!!!!!!!!!!!!! IBM & LES VARIABLES (GRID.DAT)
      INTEGER*8              :: NINTP(4), NBODY(4)
      INTEGER*8              :: NZERO
      REAL*8                 :: NUTAVG, NUTMAX, ALPAVG, ALPMAX
!     NINTP, NBODY : THE NUMBER OF INTERPOLATION AND INNER POINTS FOR IBM
!     NZERO : NUMBER OF ZERO SGS EDDY VISCOSITY POINT
!     TZERO : NUMBER OF ZERO SGS EDDY DIFFUSIVITY POINT

!!!!!!!!!!!!!!!!!!!!! N-S EQUATION
      INTEGER*8              :: NTIME, IHIST, M, MSUB, NV, NAV
      REAL*8                 :: DT, TIME, SUBDT
      REAL*8                 :: DVMAX, CFLMAX, QMMAX
      REAL*8                 :: GAMMA(3)
      REAL*8                 :: RO(3)
      REAL*8                 :: ALPHA
      REAL*8                 :: DTCONST,DTCONSTI,TEST1,ACOEF,ACOEFI
      REAL*8, DIMENSION(:),  ALLOCATABLE    :: AIU,CIU,AIVW,CIVW, &
                                               AJV,CJV,AJUW,CJUW, &
                                               AKW,CKW,AKUV,CKUV
!     NTIME,IHIST  : TIME STEP INDEX FOR CURRENT SIMULATION, TOTAL TIME STEP INDEX
!     M,MSUB       : TIME STEP INDEX (SAME TO NTIME), RUNGE-KUTTA SUB-TIME STEP INDEX
!     NV,NAV       : FIELD ADRESS
!     DT,TIME      : TIME STEP SIZE, TIME
!     DVMAX,CFLMAX,QMMAX : MAXIMUM VELOCITY DIVERGENCE, MAXIMUM CFL NUMBER,
!                          MAXIMUM MASS SOURCE/SINK FOR IBM
!     ALPHA,GAMMA(3),RO(3) : RK3 COEFFICIENTS
!     DTCONST,DTCONSTI,TEST1,ACOEF,ACOEFI : DT FOR EACH RK3 STEP, INVERSE OF DTCONST
!                                         : POISSON CONVERGENCE CRITERION
!     AIU(:),CIU(:),AIVW(:),CIVW(:),  : COEFFICIENT FOR ADI METHOD TO SOLVE LHS
!     AJV(:),CJV(:),AJUW(:),CJUW(:),
!     AKW(:),CKW(:),AKUV(:),CKUV(:)

!!!!!!!!!!!!!!!!!!!!! CONSTANT FLOW RATE, WALL TEMPERATURE CONDITION
      INTEGER*8              :: ICH,ICONJG
      REAL*8                 :: PMI(0:4), PMIAVG, QVOL(0:2), TVOL(0:2), QFLUX, PHCAP, THCAP
!     CMFR : CONSTANT MEAN FLOW RATE OPTION
!     PMI,PMIAVG,QVOL,QFLUX,PHCAP    : DP/DX (MEAN PRESSURE GRADIENT),
!                                      AVERAGED DP/DX
!                                      VOLUME FLOW RATE(CURRENT, PREV, INTERMEDIATE)
!                                      VOLUME FLUX
!                                      CORRECTION TERM FOR VOLUME FLOW RATE CONSTANT SIMULATION

!!!!!!!!!!!!!!!!!!!!! CALCULATION-TIME MEASUREMENTS
      REAL*8                :: TOTAL_TIME_B,TOTAL_TIME_E,                          &
                               TIME_BEGIN,TIME_END,SGSTIME_B(3),SGSTIME_E(3),      &
                               RHSNLHSTIME_B(3),RHSNLHSTIME_E(3),                  &
                               POISSTIME_B(3),POISSTIME_E(3)

!!!!!!!!!!!!!!!!!!!!! NON-DIMENSIONAL QUANTITIES
      REAL*8                :: FORCESUM(3),FORCESUMA(3)
      REAL*8                :: DUDTA,DVDTA,DWDTA
!     FORCESUM(3),FORCESUMA(3) : FORCE OBTAINED FROM THE MOMENTUM FORCING IN IBM
!     DUDTA,DVDTA,DWDTA        : INERTIA CONTRIBUTION IN IBM FORCING

!!!!!!!!!!!!!!!!!!!!! FIELD AVERAGE
      REAL*8                :: TIMEINIT
      INTEGER*8             :: IHISTINIT
!     TIMEINIT,TIMEINITZ   : START AND END TIME FOR AVERAGE FIELD
!     IHISTINIT,IHISTINITZ : START AND END TIME STEP FOR AVERAGE FIELD

CONTAINS
!=======================================================================
SUBROUTINE ALLO(NN1,NN2,NN3)
!=======================================================================
      IMPLICIT NONE
      INTEGER*8, INTENT(IN) :: NN1,NN2,NN3

      N1MD = N1M*2
      N2MD = N2M*2
      N3MD = N3M*2
      N1MH = N1M/2+1
      N2MH = N2M/2+1
      N3MH = N3M/2+1

      ALLOCATE(X(N1))
      ALLOCATE(Y(N2))
      ALLOCATE(Z(N3))
      ALLOCATE(XMP(0:N1))
      ALLOCATE(YMP(0:N2))
      ALLOCATE(ZMP(0:N3))

      ALLOCATE(IPV(N1M))
      ALLOCATE(IMV(N1M))
      ALLOCATE(JPV(N2M))
      ALLOCATE(JMV(N2M))
      ALLOCATE(KPV(N3M))
      ALLOCATE(KMV(N3M))

      ALLOCATE(FIXIL(N1M))
      ALLOCATE(FIXIU(N1M))
      ALLOCATE(FIXJL(N2M))
      ALLOCATE(FIXJU(N2M))
      ALLOCATE(FIXKL(N3M))
      ALLOCATE(FIXKU(N3M))

      ALLOCATE(F2FX(0:N1))
      ALLOCATE(F2FY(0:N2))
      ALLOCATE(F2FZ(0:N3))
      ALLOCATE(F2FXI(0:N1))
      ALLOCATE(F2FYI(0:N2))
      ALLOCATE(F2FZI(0:N3))

      ALLOCATE(C2CX(0:N1))
      ALLOCATE(C2CY(0:N2))
      ALLOCATE(C2CZ(0:N3))
      ALLOCATE(C2CXI(0:N1))
      ALLOCATE(C2CYI(0:N2))
      ALLOCATE(C2CZI(0:N3))

      ALLOCATE(AIU(N1))
      ALLOCATE(CIU(N1))
      ALLOCATE(AIVW(N1))
      ALLOCATE(CIVW(N1))
      ALLOCATE(AJV(N2))
      ALLOCATE(CJV(N2))
      ALLOCATE(AJUW(N2))
      ALLOCATE(CJUW(N2))
      ALLOCATE(AKW(N3))
      ALLOCATE(CKW(N3))
      ALLOCATE(AKUV(N3))
      ALLOCATE(CKUV(N3))

      ALLOCATE(DUOUT(0:N2,0:N3))
      ALLOCATE(DVOUT(0:N2,0:N3))
      ALLOCATE(DWOUT(0:N2,0:N3))
      ALLOCATE(DTOUT(0:N2,0:N3))

      ALLOCATE(UOUT(N2M,N3M))
      ALLOCATE(VOUT(N2M,N3M))
      ALLOCATE(WOUT(N2M,N3M))
      ALLOCATE(TOUT(N2M,N3M))

END SUBROUTINE ALLO

END MODULE MOD_COMMON