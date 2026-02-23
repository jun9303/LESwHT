module mod_common
!!!!!!!!!!!!!!!!!!!!! BASIC INPUTS (SETTINGS.INPUT)
  real(8) :: re, pr, gr, eps_ptr, ubulk_i, dt_size, &
             cflfac, resid1, wwsor, csgsts, csgshf
  integer(8) :: ireset, iread, iavg, ipzero, nend, nprint, &
                npriavg, npin, idtopt, nlev, mode, nbli, &
                ioldv, mgitr, imgsor, iles, insmdl, itemdl, &
                idvmon, filter, ibmon, masson, grdir, t_inf, &
                imovingon, ihtrans, ntrace, ntr, trpts(64, 3)
  character*25 :: gridfile
  character*25 :: prev_fld

!!!!!!!!!!!!!!!!!!!!! BOUNDARY CONDITIONS (BOUNDARY.INPUT & IBMPRE_PRDIC.BIN)
  integer(8) :: xprdic, yprdic, zprdic, iintp
  integer(8) :: bc_xbtm, bc_xtop, bc_ybtm, bc_ytop, &
                bc_zbtm, bc_ztop
  integer(8) :: bc_t_ybtm, bc_t_ytop, bc_t_zbtm, bc_t_ztop
  real(8) :: val_t_ybtm, val_t_ytop, val_t_zbtm, val_t_ztop
  integer(8) :: i_bgpx, j_bgpy, k_bgpz
  real(8), dimension(:, :), allocatable :: duout, dvout, dwout, dtout
  real(8), dimension(:, :), allocatable :: uout, vout, wout, tout
  integer(8) :: jut, jub, jwt, jwb, kut, kub, kvt, kvb

!!!!!!!!!!!!!!!!!!!!! GRID & MESH & DOMAIN (GRID.DAT)
  integer(8) :: n1, n1m, n1md, n1mh, n2, n2m, n2md, n2mh, n3, n3m, n3md, n3mh
  real(8) :: xl, yl, zl
  real(8), dimension(:), allocatable :: x, y, z, xmp, ymp, zmp
  integer(8), dimension(:), allocatable :: ipv, jpv, kpv, imv, jmv, kmv
  integer(8), dimension(:), allocatable :: fixil, fixiu, fixjl, fixju, fixkl, fixku
  real(8), dimension(:), allocatable :: f2fx, f2fy, f2fz, f2fxi, f2fyi, f2fzi
  real(8), dimension(:), allocatable :: c2cx, c2cy, c2cz, c2cxi, c2cyi, c2czi
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
  integer(8) :: nintp(4), nbody(4)
  integer(8) :: nzero
  real(8) :: nutavg, nutmax, alpavg, alpmax
!     NINTP, NBODY : THE NUMBER OF INTERPOLATION AND INNER POINTS FOR IBM
!     NZERO : NUMBER OF ZERO SGS EDDY VISCOSITY POINT
!     TZERO : NUMBER OF ZERO SGS EDDY DIFFUSIVITY POINT

!!!!!!!!!!!!!!!!!!!!! N-S EQUATION
  integer(8) :: ntime, ihist, m, msub, nv, nav
  real(8) :: dt, time, subdt, tend
  real(8) :: dvmax, cflmax, qmmax
  real(8) :: gamma(3)
  real(8) :: ro(3)
  real(8) :: alpha
  real(8) :: dtconst, dtconsti, test1, acoef, acoefi
  real(8), dimension(:), allocatable :: aiu, ciu, aivw, civw, &
                                        ajv, cjv, ajuw, cjuw, &
                                        akw, ckw, akuv, ckuv
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
  integer(8) :: ich, iconjg
  real(8) :: pmi(0:4), pmiavg, qvol(0:2), tvol(0:2), qflux, phcap, thcap
!     CMFR : CONSTANT MEAN FLOW RATE OPTION
!     PMI,PMIAVG,QVOL,QFLUX,PHCAP    : DP/DX (MEAN PRESSURE GRADIENT),
!                                      AVERAGED DP/DX
!                                      VOLUME FLOW RATE(CURRENT, PREV, INTERMEDIATE)
!                                      VOLUME FLUX
!                                      CORRECTION TERM FOR VOLUME FLOW RATE CONSTANT SIMULATION

!!!!!!!!!!!!!!!!!!!!! CALCULATION-TIME MEASUREMENTS
  real(8) :: total_time_b, total_time_e, &
             time_begin, time_end, sgstime_b(3), sgstime_e(3), &
             rhsnlhstime_b(3), rhsnlhstime_e(3), &
             poisstime_b(3), poisstime_e(3)

!!!!!!!!!!!!!!!!!!!!! NON-DIMENSIONAL QUANTITIES
  real(8) :: forcesum(3), forcesuma(3)
  real(8) :: dudta, dvdta, dwdta
  real(8) :: ptb_tst
  real(8) :: avg_tst, cdavg_dur, cdavg_int(3)
  integer(8) :: npriavg_count, ihistavg_start
  logical :: avg_started
!     FORCESUM(3),FORCESUMA(3) : FORCE OBTAINED FROM THE MOMENTUM FORCING IN IBM
!     DUDTA,DVDTA,DWDTA        : INERTIA CONTRIBUTION IN IBM FORCING

!!!!!!!!!!!!!!!!!!!!! FIELD AVERAGE
  real(8) :: timeinit
  integer(8) :: ihistinit
!     TIMEINIT,TIMEINITZ   : START AND END TIME FOR AVERAGE FIELD
!     IHISTINIT,IHISTINITZ : START AND END TIME STEP FOR AVERAGE FIELD

contains
!=======================================================================
  subroutine allo(nn1, nn2, nn3)
!=======================================================================
    implicit none
    integer(8), intent(in) :: nn1, nn2, nn3

    n1md = n1m * 2
    n2md = n2m * 2
    n3md = n3m * 2
    n1mh = n1m / 2 + 1
    n2mh = n2m / 2 + 1
    n3mh = n3m / 2 + 1

    allocate (x(n1))
    allocate (y(n2))
    allocate (z(n3))
    allocate (xmp(0:n1))
    allocate (ymp(0:n2))
    allocate (zmp(0:n3))

    allocate (ipv(n1m))
    allocate (imv(n1m))
    allocate (jpv(n2m))
    allocate (jmv(n2m))
    allocate (kpv(n3m))
    allocate (kmv(n3m))

    allocate (fixil(n1m))
    allocate (fixiu(n1m))
    allocate (fixjl(n2m))
    allocate (fixju(n2m))
    allocate (fixkl(n3m))
    allocate (fixku(n3m))

    allocate (f2fx(0:n1))
    allocate (f2fy(0:n2))
    allocate (f2fz(0:n3))
    allocate (f2fxi(0:n1))
    allocate (f2fyi(0:n2))
    allocate (f2fzi(0:n3))

    allocate (c2cx(0:n1))
    allocate (c2cy(0:n2))
    allocate (c2cz(0:n3))
    allocate (c2cxi(0:n1))
    allocate (c2cyi(0:n2))
    allocate (c2czi(0:n3))

    allocate (aiu(n1))
    allocate (ciu(n1))
    allocate (aivw(n1))
    allocate (civw(n1))
    allocate (ajv(n2))
    allocate (cjv(n2))
    allocate (ajuw(n2))
    allocate (cjuw(n2))
    allocate (akw(n3))
    allocate (ckw(n3))
    allocate (akuv(n3))
    allocate (ckuv(n3))

    allocate (duout(0:n2, 0:n3))
    allocate (dvout(0:n2, 0:n3))
    allocate (dwout(0:n2, 0:n3))
    allocate (dtout(0:n2, 0:n3))

    allocate (uout(n2m, n3m))
    allocate (vout(n2m, n3m))
    allocate (wout(n2m, n3m))
    allocate (tout(n2m, n3m))

  end subroutine allo

end module mod_common
