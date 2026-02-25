"""
post_inst.py — Instantaneous-field postprocessor for LESwHT.

Python port of 03_post_inst_processor/post_inst.f90.

Usage (from the 03_post_processor/inst/ directory):
    python post_inst.py

Input:  post_inst.input   (control parameters)
    ../../output/grid/grid.dat
    ../../output/field/<fldname>   (Fortran unformatted binary)

Output: ../../output/post_inst/<tname>     (VTK and/or Tecplot ASCII)

Requires: numpy, scipy (for eigvalsh), and lib_post.py in the parent dir.
"""

from __future__ import annotations
import os
import sys
import struct
import numpy as np

# Allow importing lib_post from ../lib and f2py extension from ..
_THIS_DIR = os.path.dirname(__file__)
_POST_ROOT = os.path.abspath(os.path.join(_THIS_DIR, '..'))
_LIB_DIR = os.path.join(_POST_ROOT, 'lib')
sys.path.insert(0, _LIB_DIR)
sys.path.insert(0, _POST_ROOT)
from lib_post import (
    GridData,
    read_input_inst, read_grid, adjust_post_bounds, find_ip_jp_kp,
    compute_griduni, make_ftail, make_fhead,
    _vtk_header, write_vtk_scalar, write_vtk_vector,
)

try:
    from scipy.linalg import eigvalsh as _eigvalsh
    _HAVE_SCIPY = True
except ImportError:
    _HAVE_SCIPY = False


# ===========================================================================
# Fortran unformatted binary reader
# ===========================================================================

def _read_frec(fh) -> bytes:
    """Read one Fortran unformatted sequential record and return raw bytes."""
    raw = fh.read(4)
    if len(raw) < 4:
        raise EOFError("Unexpected end of file reading record marker")
    nbytes = struct.unpack('<i', raw)[0]
    data = fh.read(nbytes)
    if len(data) != nbytes:
        raise EOFError("Short record read")
    tail = fh.read(4)
    if len(tail) == 4 and struct.unpack('<i', tail)[0] != nbytes:
        raise ValueError("Trailing record marker mismatch")
    return data


def _rec_to_array(raw: bytes, dtype) -> np.ndarray:
    return np.frombuffer(raw, dtype=dtype)


# ===========================================================================
# Read periodic flags from the first field file header
# ===========================================================================

def read_periodic_flags(fldname: str, field_dir: str) -> tuple[int, int, int]:
    """
    Peek at the first field file to read XPRDIC, YPRDIC, ZPRDIC flags.
    Returns (XPRDIC, YPRDIC, ZPRDIC).
    """
    fpath = os.path.join(field_dir, fldname)
    with open(fpath, 'rb') as fh:
        _read_frec(fh)          # skip header-1 (NN1 NN2 NN3 RE PRA GRA)
        _read_frec(fh)          # skip header-2 (IHIST M TIME DT)
        rec = _read_frec(fh)    # XPRDIC YPRDIC ZPRDIC
        vals = _rec_to_array(rec, np.int64)
        return int(vals[0]), int(vals[1]), int(vals[2])


# ===========================================================================
# Read instantaneous field file
# ===========================================================================

def read_inst_field(
        fldname: str,
        g: GridData,
        ihtrans: int,
        xprdic: int,
        zprdic: int,
        field_dir: str = '../../output/field',
) -> tuple:
    """
    Read a Fortran unformatted instantaneous field file.

    Returns
    -------
    U, V, W : ndarray shape (N1+1, N2+1, N3+1)
    P, T    : ndarray shape (N1+1, N2+1, N3+1)
    Re, Pr, Gr : float scalars
    """
    N1, N2, N3 = g.N1, g.N2, g.N3
    N1M, N2M, N3M = g.N1M, g.N2M, g.N3M

    U = np.zeros((N1 + 1, N2 + 1, N3 + 1))
    V = np.zeros((N1 + 1, N2 + 1, N3 + 1))
    W = np.zeros((N1 + 1, N2 + 1, N3 + 1))
    P = np.zeros((N1 + 1, N2 + 1, N3 + 1))
    T = np.zeros((N1 + 1, N2 + 1, N3 + 1))

    fpath = os.path.join(field_dir, fldname)
    print(f"  Opening: {fpath}")

    with open(fpath, 'rb') as fh:
        # --- Record 1: NN1 NN2 NN3 RRE PRA GRA ---
        rec = _read_frec(fh)
        h1 = _rec_to_array(rec[:24], np.int64)   # 3 × int64
        h1f = _rec_to_array(rec[24:], np.float64) # 3 × float64
        NN1, NN2, NN3 = int(h1[0]), int(h1[1]), int(h1[2])
        Re, Pr, Gr = float(h1f[0]), float(h1f[1]), float(h1f[2])

        # --- Record 2: IHIST M TIME DT ---
        rec = _read_frec(fh)
        r2i = _rec_to_array(rec[:16], np.int64)
        r2f = _rec_to_array(rec[16:], np.float64)
        IHIST, M_step = int(r2i[0]), int(r2i[1])
        TIME, DT = float(r2f[0]), float(r2f[1])

        # --- Record 3: XPRDIC YPRDIC ZPRDIC ---
        rec = _read_frec(fh)
        prd = _rec_to_array(rec, np.int64)
        xprdic_f, yprdic_f, zprdic_f = int(prd[0]), int(prd[1]), int(prd[2])

        # --- Record 4: BC flags ---
        _read_frec(fh)

        # --- Record 5: ICH ITEMP ---
        _read_frec(fh)

        if NN1 != N1 or NN2 != N2 or NN3 != N3:
            raise ValueError(
                f"Grid/field size mismatch: grid ({N1},{N2},{N3}) "
                f"vs field ({NN1},{NN2},{NN3})"
            )

        # --- Record 6: U(I=1:N1, J=0:N2, K=0:N3) ---
        # I is the fastest-varying index (Fortran column-major)
        rec = _read_frec(fh)
        u_flat = _rec_to_array(rec, np.float64)
        # reshape: innermost=I (size N1), J=0..N2 (size N2+1), K=0..N3 (size N3+1)
        u_3d = u_flat.reshape(N3 + 1, N2 + 1, N1, order='C')
        # u_3d[k, j, i_0indexed] = U(i_0indexed+1, j, k)
        U[1:N1 + 1, 0:N2 + 1, 0:N3 + 1] = u_3d.transpose(2, 1, 0)

        # --- Record 7: V(I=0:N1, J=1:N2, K=0:N3) ---
        rec = _read_frec(fh)
        v_flat = _rec_to_array(rec, np.float64)
        v_3d = v_flat.reshape(N3 + 1, N2, N1 + 1, order='C')
        V[0:N1 + 1, 1:N2 + 1, 0:N3 + 1] = v_3d.transpose(2, 1, 0)

        # --- Record 8: W(I=0:N1, J=0:N2, K=1:N3) ---
        rec = _read_frec(fh)
        w_flat = _rec_to_array(rec, np.float64)
        w_3d = w_flat.reshape(N3, N2 + 1, N1 + 1, order='C')
        W[0:N1 + 1, 0:N2 + 1, 1:N3 + 1] = w_3d.transpose(2, 1, 0)

        # --- Record 9: P(I=1:N1M, J=1:N2M, K=1:N3M) ---
        rec = _read_frec(fh)
        p_flat = _rec_to_array(rec, np.float64)
        p_3d = p_flat.reshape(N3M, N2M, N1M, order='C')
        P[1:N1M + 1, 1:N2M + 1, 1:N3M + 1] = p_3d.transpose(2, 1, 0)

        # --- Record 10 (optional): T ---
        if ihtrans == 1:
            rec = _read_frec(fh)
            t_flat = _rec_to_array(rec, np.float64)
            t_3d = t_flat.reshape(N3M, N2M, N1M, order='C')
            T[1:N1M + 1, 1:N2M + 1, 1:N3M + 1] = t_3d.transpose(2, 1, 0)

    print("----------- INITIAL FIELD INFORMATION -----------")
    print("INITIAL FIELD      : READING DONE")
    print(f"INITIAL FIELD NAME : {fldname}")
    print(f"RE={Re:.6e}")
    print(f"N1={NN1:10d}  N2={NN2:10d}  N3={NN3:10d}")
    print(f"IHIST={IHIST:8d}  M={M_step:10d}  TIME= {TIME:12.5f}  DT={DT:12.5f}")

    # --- Apply periodicity ---
    if zprdic_f == 1 or zprdic == 1:
        # U(I, J, 0) = U(I, J, N3M); U(I, J, N3) = U(I, J, 1)
        U[1:N1 + 1, 0:N2 + 1, 0]  = U[1:N1 + 1, 0:N2 + 1, N3M]
        U[1:N1 + 1, 0:N2 + 1, N3] = U[1:N1 + 1, 0:N2 + 1, 1]
        V[0:N1 + 1, 1:N2 + 1, 0]  = V[0:N1 + 1, 1:N2 + 1, N3M]
        V[0:N1 + 1, 1:N2 + 1, N3] = V[0:N1 + 1, 1:N2 + 1, 1]
        W[0:N1 + 1, 0:N2 + 1, 0]  = W[0:N1 + 1, 0:N2 + 1, N3M]
        W[0:N1 + 1, 0:N2 + 1, N3] = W[0:N1 + 1, 0:N2 + 1, 1]
        P[1:N1M + 1, 1:N2M + 1, 0]  = P[1:N1M + 1, 1:N2M + 1, N3M]
        P[1:N1M + 1, 1:N2M + 1, N3] = P[1:N1M + 1, 1:N2M + 1, 1]
        if ihtrans == 1:
            T[1:N1M + 1, 1:N2M + 1, 0]  = T[1:N1M + 1, 1:N2M + 1, N3M]
            T[1:N1M + 1, 1:N2M + 1, N3] = T[1:N1M + 1, 1:N2M + 1, 1]

    if xprdic_f == 1 or xprdic == 1:
        U[0,  0:N2 + 1, 0:N3 + 1] = U[N1M, 0:N2 + 1, 0:N3 + 1]
        U[N1, 0:N2 + 1, 0:N3 + 1] = U[1,   0:N2 + 1, 0:N3 + 1]
        V[0,  1:N2 + 1, 0:N3 + 1] = V[N1M, 1:N2 + 1, 0:N3 + 1]
        V[N1, 1:N2 + 1, 0:N3 + 1] = V[1,   1:N2 + 1, 0:N3 + 1]
        W[0,  0:N2 + 1, 1:N3 + 1] = W[N1M, 0:N2 + 1, 1:N3 + 1]
        W[N1, 0:N2 + 1, 1:N3 + 1] = W[1,   0:N2 + 1, 1:N3 + 1]
        P[0,  1:N2M + 1, 1:N3M + 1] = P[N1M, 1:N2M + 1, 1:N3M + 1]
        P[N1, 1:N2M + 1, 1:N3M + 1] = P[1,   1:N2M + 1, 1:N3M + 1]
        if ihtrans == 1:
            T[0,  1:N2M + 1, 1:N3M + 1] = T[N1M, 1:N2M + 1, 1:N3M + 1]
            T[N1, 1:N2M + 1, 1:N3M + 1] = T[1,   1:N2M + 1, 1:N3M + 1]

    return U, V, W, P, T, Re, Pr, Gr


# ===========================================================================
# Vorticity and λ₂ computation
# ===========================================================================

def compute_vorn_lambda2(
        g: GridData,
        U: np.ndarray,
        V: np.ndarray,
        W: np.ndarray,
    P: np.ndarray,
        ibmon: int,
) -> tuple:
    """
    Compute:
      UC, VC, WC  – cell-centre velocities      (shape N1M×N2M×N3M, 1-indexed)
      VORX,Y,Z    – vorticity components         (same shape)
      VLAMBDA2    – λ₂ criterion                 (same shape)

    Returns 1-indexed arrays with unused index-0 set to zero.
    """
    N1, N2, N3 = g.N1, g.N2, g.N3
    N1M, N2M, N3M = g.N1M, g.N2M, g.N3M

    # ------------------------------------------------------------------ #
    # Cell-centre velocities                                              #
    #   UC(i,j,k) = 0.5*(U(i,j,k) + U(i+1,j,k))  etc.                  #
    # ------------------------------------------------------------------ #
    I = slice(1, N1M + 1); J = slice(1, N2M + 1); K = slice(1, N3M + 1)
    Ip = slice(2, N1 + 1); Jp = slice(2, N2 + 1); Kp = slice(2, N3 + 1)
    Im = slice(0, N1M);    Jm = slice(0, N2M);    Km = slice(0, N3M)

    UC = np.zeros((N1 + 1, N2 + 1, N3 + 1))
    VC = np.zeros((N1 + 1, N2 + 1, N3 + 1))
    WC = np.zeros((N1 + 1, N2 + 1, N3 + 1))
    UC[I, J, K] = 0.5 * (U[I, J, K] + U[Ip, J, K])
    VC[I, J, K] = 0.5 * (V[I, J, K] + V[I, Jp, K])
    WC[I, J, K] = 0.5 * (W[I, J, K] + W[I, J, Kp])

    # ------------------------------------------------------------------ #
    # Velocity gradient tensor at cell centres (vectorised)               #
    # ------------------------------------------------------------------ #
    SSDX_I  = g.SSDX[I][:, None, None]
    SSDY_J  = g.SSDY[J][None, :, None]
    SSDZ_K  = g.SSDZ[K][None, None, :]
    VVDX_I  = g.VVDX[I][:, None, None]
    VVDX_Ip = g.VVDX[2:N1 + 1][:, None, None]
    VVDY_J  = g.VVDY[J][None, :, None]
    VVDY_Jp = g.VVDY[2:N2 + 1][None, :, None]
    VVDZ_K  = g.VVDZ[K][None, None, :]
    VVDZ_Kp = g.VVDZ[2:N3 + 1][None, None, :]
    SDX_Im  = g.SDX[Im][:, None, None]
    SDX_I   = g.SDX[I][:, None, None]
    SDX_Ip  = g.SDX[2:N1 + 1][:, None, None]
    SDY_Jm  = g.SDY[Jm][None, :, None]
    SDY_J   = g.SDY[J][None, :, None]
    SDY_Jp  = g.SDY[2:N2 + 1][None, :, None]
    SDZ_Km  = g.SDZ[Km][None, None, :]
    SDZ_K   = g.SDZ[K][None, None, :]
    SDZ_Kp  = g.SDZ[2:N3 + 1][None, None, :]
    FIXIL_I = g.FIXIL[I][:, None, None]
    FIXIU_I = g.FIXIU[I][:, None, None]
    FIXJL_J = g.FIXJL[J][None, :, None]
    FIXJU_J = g.FIXJU[J][None, :, None]
    FIXKL_K = g.FIXKL[K][None, None, :]
    FIXKU_K = g.FIXKU[K][None, None, :]

    # VG11 = dU/dx
    VG11 = SSDX_I * (U[Ip, J, K] - U[I, J, K])

    # VG22 = dV/dy
    VG22 = SSDY_J * (V[I, Jp, K] - V[I, J, K])

    # VG33 = dW/dz
    VG33 = SSDZ_K * (W[I, J, Kp] - W[I, J, K])

    # VG12 = dU/dy  (interpolate U to cell-centre across j faces)
    UI   = U[I, J, K]  + U[Ip, J, K]          # (N1M, N2M, N3M)
    UIjp = U[I, Jp, K] + U[Ip, Jp, K]
    UIjm = U[I, Jm, K] + U[Ip, Jm, K]
    UI_N2 = (U[I, N2, K] + U[Ip, N2, K])[:, None, :]  # boundary
    UI_0y = (U[I, 0, K]  + U[Ip, 0, K])[:, None, :]
    UP12 = (VVDY_Jp * 0.25 * (SDY_Jp * UI + SDY_J * UIjp) * (1 - FIXJU_J)
            + 0.5 * UI_N2 * FIXJU_J)
    UM12 = (VVDY_J  * 0.25 * (SDY_J * UIjm + SDY_Jm * UI) * (1 - FIXJL_J)
            + 0.5 * UI_0y * FIXJL_J)
    VG12 = SSDY_J * (UP12 - UM12)

    # VG13 = dU/dz
    UIkp = U[I, J, Kp] + U[Ip, J, Kp]
    UIkm = U[I, J, Km] + U[Ip, J, Km]
    UI_N3 = (U[I, J, N3] + U[Ip, J, N3])[:, :, None]
    UI_0z = (U[I, J, 0]  + U[Ip, J, 0])[:, :, None]
    UP13 = (VVDZ_Kp * 0.25 * (SDZ_Kp * UI + SDZ_K * UIkp) * (1 - FIXKU_K)
            + 0.5 * UI_N3 * FIXKU_K)
    UM13 = (VVDZ_K  * 0.25 * (SDZ_K * UIkm + SDZ_Km * UI) * (1 - FIXKL_K)
            + 0.5 * UI_0z * FIXKL_K)
    VG13 = SSDZ_K * (UP13 - UM13)

    # VG21 = dV/dx
    VI   = V[I, J, K]  + V[I, Jp, K]
    VIip = V[Ip, J, K] + V[Ip, Jp, K]
    VIim = V[Im, J, K] + V[Im, Jp, K]
    VI_N1 = (V[N1, J, K] + V[N1, Jp, K])[None, :, :]
    VI_0x = (V[0,  J, K] + V[0,  Jp, K])[None, :, :]
    UP21 = (VVDX_Ip * 0.25 * (SDX_Ip * VI + SDX_I * VIip) * (1 - FIXIU_I)
            + 0.5 * VI_N1 * FIXIU_I)
    UM21 = (VVDX_I  * 0.25 * (SDX_I * VIim + SDX_Im * VI) * (1 - FIXIL_I)
            + 0.5 * VI_0x * FIXIL_I)
    VG21 = SSDX_I * (UP21 - UM21)

    # VG23 = dV/dz
    VIkp = V[I, J, Kp] + V[I, Jp, Kp]
    VIkm = V[I, J, Km] + V[I, Jp, Km]
    VI_N3 = (V[I, J, N3] + V[I, Jp, N3])[:, :, None]
    VI_0z = (V[I, J, 0]  + V[I, Jp, 0])[:, :, None]
    UP23 = (VVDZ_Kp * 0.25 * (SDZ_Kp * VI + SDZ_K * VIkp) * (1 - FIXKU_K)
            + 0.5 * VI_N3 * FIXKU_K)
    UM23 = (VVDZ_K  * 0.25 * (SDZ_K * VIkm + SDZ_Km * VI) * (1 - FIXKL_K)
            + 0.5 * VI_0z * FIXKL_K)
    VG23 = SSDZ_K * (UP23 - UM23)

    # VG31 = dW/dx
    WI   = W[I, J, K]  + W[I, J, Kp]
    WIip = W[Ip, J, K] + W[Ip, J, Kp]
    WIim = W[Im, J, K] + W[Im, J, Kp]
    WI_N1 = (W[N1, J, K] + W[N1, J, Kp])[None, :, :]
    WI_0x = (W[0,  J, K] + W[0,  J, Kp])[None, :, :]
    UP31 = (VVDX_Ip * 0.25 * (SDX_Ip * WI + SDX_I * WIip) * (1 - FIXIU_I)
            + 0.5 * WI_N1 * FIXIU_I)
    UM31 = (VVDX_I  * 0.25 * (SDX_I * WIim + SDX_Im * WI) * (1 - FIXIL_I)
            + 0.5 * WI_0x * FIXIL_I)
    VG31 = SSDX_I * (UP31 - UM31)

    # VG32 = dW/dy
    WIjp = W[I, Jp, K] + W[I, Jp, Kp]
    WIjm = W[I, Jm, K] + W[I, Jm, Kp]
    WI_N2 = (W[I, N2, K] + W[I, N2, Kp])[:, None, :]
    WI_0y = (W[I, 0, K]  + W[I, 0, Kp])[:, None, :]
    UP32 = (VVDY_Jp * 0.25 * (SDY_Jp * WI + SDY_J * WIjp) * (1 - FIXJU_J)
            + 0.5 * WI_N2 * FIXJU_J)
    UM32 = (VVDY_J  * 0.25 * (SDY_J * WIjm + SDY_Jm * WI) * (1 - FIXJL_J)
            + 0.5 * WI_0y * FIXJL_J)
    VG32 = SSDY_J * (UP32 - UM32)

    # ------------------------------------------------------------------ #
    # Vorticity: WX = VG32-VG23, WY = VG13-VG31, WZ = VG21-VG12          #
    # ------------------------------------------------------------------ #
    VORX = np.zeros((N1 + 1, N2 + 1, N3 + 1))
    VORY = np.zeros((N1 + 1, N2 + 1, N3 + 1))
    VORZ = np.zeros((N1 + 1, N2 + 1, N3 + 1))
    VORX[I, J, K] = VG32 - VG23
    VORY[I, J, K] = VG13 - VG31
    VORZ[I, J, K] = VG21 - VG12

    # ------------------------------------------------------------------ #
    # Lambda-2 criterion (vectorised eigenvalue computation)              #
    # ------------------------------------------------------------------ #
    VLAMBDA2 = np.zeros((N1 + 1, N2 + 1, N3 + 1))

    # Symmetric S and antisymmetric OM parts of velocity gradient
    # S[a,b] = 0.5*(VGab + VGba),  OM[a,b] = 0.5*(VGab - VGba)
    # M = S²+OM² — symmetric 3×3

    # Vectorised: pack into (N1M*N2M*N3M) × 3 × 3 array, compute eigvals
    nc = N1M * N2M * N3M
    # Flatten all cells
    def _f(a): return a.reshape(nc)
    vg = {(1,1): _f(VG11), (1,2): _f(VG12), (1,3): _f(VG13),
          (2,1): _f(VG21), (2,2): _f(VG22), (2,3): _f(VG23),
          (3,1): _f(VG31), (3,2): _f(VG32), (3,3): _f(VG33)}

    def S(a, b):  return 0.5 * (vg[a, b] + vg[b, a])
    def OM(a, b): return 0.5 * (vg[a, b] - vg[b, a])

    # Build M = S²+OM² (symmetric) element by element
    M = np.zeros((nc, 3, 3))
    for a in range(3):
        for b in range(3):
            aa, bb = a + 1, b + 1
            for c in range(3):
                cc = c + 1
                M[:, a, b] += S(aa, cc) * S(cc, bb) + OM(aa, cc) * OM(cc, bb)

    # Eigenvalues of symmetric M
    if _HAVE_SCIPY:
        from scipy.linalg import eigvalsh as _ev
        eigs = np.array([_ev(M[n]) for n in range(nc)])  # (nc, 3)
    else:
        eigs = np.linalg.eigvalsh(M)                      # (nc, 3)

    eigs.sort(axis=1)   # ascending
    # Lambda2 = 2nd eigenvalue (index 1) if negative, else 0
    lam2 = np.where(eigs[:, 1] < 0, eigs[:, 1], 0.0)
    VLAMBDA2[I, J, K] = lam2.reshape(N1M, N2M, N3M)

    return UC, VC, WC, VORX, VORY, VORZ, VLAMBDA2


# ===========================================================================
# VTK plane output helpers
# ===========================================================================

def _iter_jk(N2M, N3M):
    """Iterate (j, k) in the order Fortran writes: outer j, inner k."""
    for j in range(1, N2M + 1):
        for k in range(1, N3M + 1):
            yield j, k


def _iter_ik(N1M, N3M):
    for k in range(1, N3M + 1):
        for i in range(1, N1M + 1):
            yield i, k


def _iter_ij(N1M, N2M):
    for j in range(1, N2M + 1):
        for i in range(1, N1M + 1):
            yield i, j


def output_i_vtk(
        ipoint: int, g: GridData,
        UC, VC, WC, P, VORX, VORY, VORZ, VLAMBDA2,
        tname1: str, ftailijk: str,
        ihtrans: int, T,
        outdir: str = '../../output/post_inst',
) -> None:
    """Write 2D x-plane VTK file."""
    N2M, N3M = g.N2M, g.N3M
    ii = ipoint
    fpath = os.path.join(outdir, tname1 + ftailijk)

    with open(fpath, 'w') as fh:
        _vtk_header(fh, 'LESwHT instantaneous x-plane',
                    (N3M, N2M, 1), N3M * N2M)
        for j, k in _iter_jk(N2M, N3M):
            fh.write(f" {g.XMP[ii]:.16e} {g.YMP[j]:.16e} {g.ZMP[k]:.16e}\n")
        fh.write(f"POINT_DATA {N3M * N2M}\n")
        write_vtk_vector(fh, 'velocity',
                         [UC[ii, j, k] for j, k in _iter_jk(N2M, N3M)],
                         [VC[ii, j, k] for j, k in _iter_jk(N2M, N3M)],
                         [WC[ii, j, k] for j, k in _iter_jk(N2M, N3M)])
        write_vtk_scalar(fh, 'pressure',
                         [P[ii, j, k] for j, k in _iter_jk(N2M, N3M)])
        write_vtk_scalar(fh, 'vort_x',
                         [VORX[ii, j, k] for j, k in _iter_jk(N2M, N3M)])
        write_vtk_scalar(fh, 'vort_y',
                         [VORY[ii, j, k] for j, k in _iter_jk(N2M, N3M)])
        write_vtk_scalar(fh, 'vort_z',
                         [VORZ[ii, j, k] for j, k in _iter_jk(N2M, N3M)])
        write_vtk_scalar(fh, 'vort_mag',
                         [np.sqrt(VORX[ii,j,k]**2 + VORY[ii,j,k]**2 + VORZ[ii,j,k]**2)
                          for j, k in _iter_jk(N2M, N3M)])
        write_vtk_scalar(fh, 'lambda2',
                         [VLAMBDA2[ii, j, k] for j, k in _iter_jk(N2M, N3M)])
        if ihtrans == 1:
            write_vtk_scalar(fh, 'temperature',
                             [T[ii, j, k] for j, k in _iter_jk(N2M, N3M)])
    print(f"  Wrote {fpath}")


def output_j_vtk(
        jpoint: int, g: GridData,
        UC, VC, WC, P, VORX, VORY, VORZ, VLAMBDA2,
        tname1: str, ftailijk: str,
        ihtrans: int, T,
        outdir: str = '../../output/post_inst',
) -> None:
    """Write 2D y-plane VTK file."""
    N1M, N3M = g.N1M, g.N3M
    jj = jpoint
    fpath = os.path.join(outdir, tname1 + ftailijk)

    with open(fpath, 'w') as fh:
        _vtk_header(fh, 'LESwHT instantaneous y-plane',
                    (N1M, N3M, 1), N1M * N3M)
        for i, k in _iter_ik(N1M, N3M):
            fh.write(f" {g.XMP[i]:.16e} {g.YMP[jj]:.16e} {g.ZMP[k]:.16e}\n")
        fh.write(f"POINT_DATA {N1M * N3M}\n")
        write_vtk_vector(fh, 'velocity',
                         [UC[i, jj, k] for i, k in _iter_ik(N1M, N3M)],
                         [VC[i, jj, k] for i, k in _iter_ik(N1M, N3M)],
                         [WC[i, jj, k] for i, k in _iter_ik(N1M, N3M)])
        write_vtk_scalar(fh, 'pressure',
                         [P[i, jj, k] for i, k in _iter_ik(N1M, N3M)])
        write_vtk_scalar(fh, 'vort_x',
                         [VORX[i, jj, k] for i, k in _iter_ik(N1M, N3M)])
        write_vtk_scalar(fh, 'vort_y',
                         [VORY[i, jj, k] for i, k in _iter_ik(N1M, N3M)])
        write_vtk_scalar(fh, 'vort_z',
                         [VORZ[i, jj, k] for i, k in _iter_ik(N1M, N3M)])
        write_vtk_scalar(fh, 'vort_mag',
                         [np.sqrt(VORX[i,jj,k]**2 + VORY[i,jj,k]**2 + VORZ[i,jj,k]**2)
                          for i, k in _iter_ik(N1M, N3M)])
        write_vtk_scalar(fh, 'lambda2',
                         [VLAMBDA2[i, jj, k] for i, k in _iter_ik(N1M, N3M)])
        if ihtrans == 1:
            write_vtk_scalar(fh, 'temperature',
                             [T[i, jj, k] for i, k in _iter_ik(N1M, N3M)])
    print(f"  Wrote {fpath}")


def output_k_vtk(
        kpoint: int, g: GridData,
        UC, VC, WC, P, VORX, VORY, VORZ, VLAMBDA2,
        tname1: str, ftailijk: str,
        ihtrans: int, T,
        outdir: str = '../../output/post_inst',
) -> None:
    """Write 2D z-plane VTK file."""
    N1M, N2M = g.N1M, g.N2M
    kk = kpoint
    fpath = os.path.join(outdir, tname1 + ftailijk)

    with open(fpath, 'w') as fh:
        _vtk_header(fh, 'LESwHT instantaneous z-plane',
                    (N1M, N2M, 1), N1M * N2M)
        for i, j in _iter_ij(N1M, N2M):
            fh.write(f" {g.XMP[i]:.16e} {g.YMP[j]:.16e} {g.ZMP[kk]:.16e}\n")
        fh.write(f"POINT_DATA {N1M * N2M}\n")
        write_vtk_vector(fh, 'velocity',
                         [UC[i, j, kk] for i, j in _iter_ij(N1M, N2M)],
                         [VC[i, j, kk] for i, j in _iter_ij(N1M, N2M)],
                         [WC[i, j, kk] for i, j in _iter_ij(N1M, N2M)])
        if ihtrans == 1:
            write_vtk_scalar(fh, 'temperature',
                             [T[i, j, kk] for i, j in _iter_ij(N1M, N2M)])
        write_vtk_scalar(fh, 'pressure',
                         [P[i, j, kk] for i, j in _iter_ij(N1M, N2M)])
        write_vtk_scalar(fh, 'vort_x',
                         [VORX[i, j, kk] for i, j in _iter_ij(N1M, N2M)])
        write_vtk_scalar(fh, 'vort_y',
                         [VORY[i, j, kk] for i, j in _iter_ij(N1M, N2M)])
        write_vtk_scalar(fh, 'vort_z',
                         [VORZ[i, j, kk] for i, j in _iter_ij(N1M, N2M)])
        write_vtk_scalar(fh, 'vort_mag',
                         [np.sqrt(VORX[i,j,kk]**2 + VORY[i,j,kk]**2 + VORZ[i,j,kk]**2)
                          for i, j in _iter_ij(N1M, N2M)])
        write_vtk_scalar(fh, 'lambda2',
                         [VLAMBDA2[i, j, kk] for i, j in _iter_ij(N1M, N2M)])
    print(f"  Wrote {fpath}")


# ===========================================================================
# Uniform-grid 2D output helpers
# ===========================================================================

def _interp2_jk(g, arr2d_jk, II):
    """Bilinear interpolation from non-uniform arr(i=II, J, K) to uniform JK."""
    N2M, N3M = g.N2M, g.N3M
    out = np.zeros((N2M + 1, N3M + 1))
    for j in range(1, N2M + 1, 2):
        for k in range(1, N3M + 1, 2):
            jj = g.JUNI[j]; kk = g.KUNI[k]
            fy = g.FACUNIY[j]; fz = g.FACUNIZ[k]
            out[j, k] = (fy * fz * arr2d_jk[II, jj + 1, kk + 1]
                         + fy * (1 - fz) * arr2d_jk[II, jj + 1, kk]
                         + (1 - fy) * fz * arr2d_jk[II, jj,     kk + 1]
                         + (1 - fy) * (1 - fz) * arr2d_jk[II, jj, kk])
    return out


def _interp2_ik(g, arr3d, JJ, use_V_or_U='U'):
    N1M, N3M = g.N1M, g.N3M
    out = np.zeros((N1M + 1, N3M + 1))
    for i in range(1, N1M + 1, 2):
        for k in range(1, N3M + 1, 2):
            ii = g.IUNI[i]; kk = g.KUNI[k]
            fx = g.FACUNIX[i]; fz = g.FACUNIZ[k]
            out[i, k] = (fx * fz * arr3d[ii + 1, JJ, kk + 1]
                         + fx * (1 - fz) * arr3d[ii + 1, JJ, kk]
                         + (1 - fx) * fz * arr3d[ii, JJ, kk + 1]
                         + (1 - fx) * (1 - fz) * arr3d[ii, JJ, kk])
    return out


def _interp2_ij(g, arr3d, KK):
    N1M, N2M = g.N1M, g.N2M
    out = np.zeros((N1M + 1, N2M + 1))
    for i in range(1, N1M + 1, 2):
        for j in range(1, N2M + 1, 2):
            ii = g.IUNI[i]; jj = g.JUNI[j]
            fx = g.FACUNIX[i]; fy = g.FACUNIY[j]
            out[i, j] = (fx * fy * arr3d[ii + 1, jj + 1, KK]
                         + fx * (1 - fy) * arr3d[ii + 1, jj, KK]
                         + (1 - fx) * fy * arr3d[ii, jj + 1, KK]
                         + (1 - fx) * (1 - fy) * arr3d[ii, jj, KK])
    return out


def output_i_uni_vtk(
        ipoint: int, g: GridData,
        U, V, W,
        tname3: str, ftailijk: str,
        outdir: str = '../../output/post_inst',
) -> None:
    """Write uniform-grid x-plane VTK (V, W interpolated)."""
    N2M, N3M = g.N2M, g.N3M
    II = ipoint
    VXX = _interp2_jk(g, V, II)
    WXX = _interp2_jk(g, W, II)
    fpath = os.path.join(outdir, tname3 + ftailijk)

    with open(fpath, 'w') as fh:
        _vtk_header(fh, 'LESwHT instantaneous x-plane uniform',
                    (g.KNUMUNI, g.JNUMUNI, 1), g.KNUMUNI * g.JNUMUNI)
        for j in range(1, N2M + 1, 2):
            for k in range(1, N3M + 1, 2):
                fh.write(f" {g.XMP[II]:.16e} {g.YUNI[j]:.16e} {g.ZUNI[k]:.16e}\n")
        fh.write(f"POINT_DATA {g.KNUMUNI * g.JNUMUNI}\n")
        write_vtk_scalar(fh, 'v',
                         [VXX[j, k] for j in range(1, N2M + 1, 2)
                          for k in range(1, N3M + 1, 2)])
        write_vtk_scalar(fh, 'w',
                         [WXX[j, k] for j in range(1, N2M + 1, 2)
                          for k in range(1, N3M + 1, 2)])
    print(f"  Wrote {fpath}")


def output_j_uni_vtk(
        jpoint: int, g: GridData,
        U, W,
        tname3: str, ftailijk: str,
        outdir: str = '../../output/post_inst',
) -> None:
    """Write uniform-grid y-plane VTK (U, W interpolated)."""
    N1M, N3M = g.N1M, g.N3M
    JJ = jpoint
    UXX = _interp2_ik(g, U, JJ)
    WXX = _interp2_ik(g, W, JJ)
    fpath = os.path.join(outdir, tname3 + ftailijk)

    with open(fpath, 'w') as fh:
        _vtk_header(fh, 'LESwHT instantaneous y-plane uniform',
                    (g.INUMUNI, g.KNUMUNI, 1), g.INUMUNI * g.KNUMUNI)
        for k in range(1, N3M + 1, 2):
            for i in range(1, N1M + 1, 2):
                fh.write(f" {g.XUNI[i]:.16e} {g.YMP[JJ]:.16e} {g.ZUNI[k]:.16e}\n")
        fh.write(f"POINT_DATA {g.INUMUNI * g.KNUMUNI}\n")
        write_vtk_scalar(fh, 'u',
                         [UXX[i, k] for k in range(1, N3M + 1, 2)
                          for i in range(1, N1M + 1, 2)])
        write_vtk_scalar(fh, 'w',
                         [WXX[i, k] for k in range(1, N3M + 1, 2)
                          for i in range(1, N1M + 1, 2)])
    print(f"  Wrote {fpath}")


def output_k_uni_vtk(
        kpoint: int, g: GridData,
        U, V,
        tname3: str, ftailijk: str,
        outdir: str = '../../output/post_inst',
) -> None:
    """Write uniform-grid z-plane VTK (U, V interpolated)."""
    N1M, N2M = g.N1M, g.N2M
    KK = kpoint
    UXX = _interp2_ij(g, U, KK)
    VXX = _interp2_ij(g, V, KK)
    fpath = os.path.join(outdir, tname3 + ftailijk)

    with open(fpath, 'w') as fh:
        _vtk_header(fh, 'LESwHT instantaneous z-plane uniform',
                    (g.INUMUNI, g.JNUMUNI, 1), g.INUMUNI * g.JNUMUNI)
        for j in range(1, N2M + 1, 2):
            for i in range(1, N1M + 1, 2):
                fh.write(f" {g.XUNI[i]:.16e} {g.YUNI[j]:.16e} {g.ZMP[KK]:.16e}\n")
        fh.write(f"POINT_DATA {g.INUMUNI * g.JNUMUNI}\n")
        write_vtk_scalar(fh, 'u',
                         [UXX[i, j] for j in range(1, N2M + 1, 2)
                          for i in range(1, N1M + 1, 2)])
        write_vtk_scalar(fh, 'v',
                         [VXX[i, j] for j in range(1, N2M + 1, 2)
                          for i in range(1, N1M + 1, 2)])
    print(f"  Wrote {fpath}")


# ===========================================================================
# Tecplot ASCII output
# ===========================================================================

def output_i_tec(ipoint, g, UC, VC, WC, P, VORX, VORY, VORZ, VLAMBDA2,
                 tname1, ftailijk, ihtrans, T, outdir='../../output/post_inst'):
    ftail = ftailijk[:5] + '.tec'
    fpath = os.path.join(outdir, tname1 + ftail)
    N2M, N3M = g.N2M, g.N3M
    ii = ipoint
    with open(fpath, 'w') as fh:
        if ihtrans == 1:
            fh.write('VARIABLES = "X","Y","Z","U","V","W","P","VORX","VORY","VORZ","VORMAG","LAMBDA2","T"\n')
        else:
            fh.write('VARIABLES = "X","Y","Z","U","V","W","P","VORX","VORY","VORZ","VORMAG","LAMBDA2"\n')
        fh.write(f'ZONE I={N3M}, J={N2M}, F=POINT\n')
        for j in range(1, N2M + 1):
            for k in range(1, N3M + 1):
                vm = np.sqrt(VORX[ii,j,k]**2 + VORY[ii,j,k]**2 + VORZ[ii,j,k]**2)
                row = (f" {g.XMP[ii]:.16e} {g.YMP[j]:.16e} {g.ZMP[k]:.16e}"
                       f" {UC[ii,j,k]:.16e} {VC[ii,j,k]:.16e} {WC[ii,j,k]:.16e}"
                       f" {P[ii,j,k]:.16e}"
                       f" {VORX[ii,j,k]:.16e} {VORY[ii,j,k]:.16e} {VORZ[ii,j,k]:.16e}"
                       f" {vm:.16e} {VLAMBDA2[ii,j,k]:.16e}")
                if ihtrans == 1:
                    row += f" {T[ii,j,k]:.16e}"
                fh.write(row + '\n')
    print(f"  Wrote {fpath}")


def output_j_tec(jpoint, g, UC, VC, WC, P, VORX, VORY, VORZ, VLAMBDA2,
                 tname1, ftailijk, ihtrans, T, outdir='../../output/post_inst'):
    ftail = ftailijk[:5] + '.tec'
    fpath = os.path.join(outdir, tname1 + ftail)
    N1M, N3M = g.N1M, g.N3M
    jj = jpoint
    with open(fpath, 'w') as fh:
        if ihtrans == 1:
            fh.write('VARIABLES = "X","Y","Z","U","V","W","P","VORX","VORY","VORZ","VORMAG","LAMBDA2","T"\n')
        else:
            fh.write('VARIABLES = "X","Y","Z","U","V","W","P","VORX","VORY","VORZ","VORMAG","LAMBDA2"\n')
        fh.write(f'ZONE I={N1M}, J={N3M}, F=POINT\n')
        for k in range(1, N3M + 1):
            for i in range(1, N1M + 1):
                vm = np.sqrt(VORX[i,jj,k]**2 + VORY[i,jj,k]**2 + VORZ[i,jj,k]**2)
                row = (f" {g.XMP[i]:.16e} {g.YMP[jj]:.16e} {g.ZMP[k]:.16e}"
                       f" {UC[i,jj,k]:.16e} {VC[i,jj,k]:.16e} {WC[i,jj,k]:.16e}"
                       f" {P[i,jj,k]:.16e}"
                       f" {VORX[i,jj,k]:.16e} {VORY[i,jj,k]:.16e} {VORZ[i,jj,k]:.16e}"
                       f" {vm:.16e} {VLAMBDA2[i,jj,k]:.16e}")
                if ihtrans == 1:
                    row += f" {T[i,jj,k]:.16e}"
                fh.write(row + '\n')
    print(f"  Wrote {fpath}")


def output_k_tec(kpoint, g, UC, VC, WC, P, VORX, VORY, VORZ, VLAMBDA2,
                 tname1, ftailijk, ihtrans, T, outdir='../../output/post_inst'):
    ftail = ftailijk[:5] + '.tec'
    fpath = os.path.join(outdir, tname1 + ftail)
    N1M, N2M = g.N1M, g.N2M
    kk = kpoint
    with open(fpath, 'w') as fh:
        if ihtrans == 1:
            fh.write('VARIABLES = "X","Y","Z","U","V","W","P","VORX","VORY","VORZ","VORMAG","LAMBDA2","T"\n')
        else:
            fh.write('VARIABLES = "X","Y","Z","U","V","W","P","VORX","VORY","VORZ","VORMAG","LAMBDA2"\n')
        fh.write(f'ZONE I={N1M}, J={N2M}, F=POINT\n')
        for j in range(1, N2M + 1):
            for i in range(1, N1M + 1):
                vm = np.sqrt(VORX[i,j,kk]**2 + VORY[i,j,kk]**2 + VORZ[i,j,kk]**2)
                row = (f" {g.XMP[i]:.16e} {g.YMP[j]:.16e} {g.ZMP[kk]:.16e}"
                       f" {UC[i,j,kk]:.16e} {VC[i,j,kk]:.16e} {WC[i,j,kk]:.16e}"
                       f" {P[i,j,kk]:.16e}"
                       f" {VORX[i,j,kk]:.16e} {VORY[i,j,kk]:.16e} {VORZ[i,j,kk]:.16e}"
                       f" {vm:.16e} {VLAMBDA2[i,j,kk]:.16e}")
                if ihtrans == 1:
                    row += f" {T[i,j,kk]:.16e}"
                fh.write(row + '\n')
    print(f"  Wrote {fpath}")


# ===========================================================================
# 3D output
# ===========================================================================

def output_3d(
        g: GridData,
        UC, VC, WC, P, VORX, VORY, VORZ, VLAMBDA2,
        tname2: str,
        params: dict,
        ihtrans: int, T,
        ioutfmt: int,
        outdir: str = '../../output/post_inst',
) -> None:
    """Write 3D output files according to ILD2, IUVWP, IWXYZ, IALL flags."""
    is_, ie = params['ISTART'], params['IEND']
    js_, je = params['JSTART'], params['JEND']
    ks_, ke = params['KSTART'], params['KEND']
    ids_, jds_, kds_ = params['ISKIP'], params['JSKIP'], params['KSKIP']

    I_rng = range(is_, ie + 1, ids_)
    J_rng = range(js_, je + 1, jds_)
    K_rng = range(ks_, ke + 1, kds_)
    INUM = len(I_rng); JNUM = len(J_rng); KNUM = len(K_rng)
    npts = INUM * JNUM * KNUM

    if ioutfmt in (1, 2):
        output_3d_tec(g, UC, VC, WC, P, VORX, VORY, VORZ, VLAMBDA2,
                      tname2, params, ihtrans, T, I_rng, J_rng, K_rng,
                      INUM, JNUM, KNUM, outdir)

    if ioutfmt in (0, 2):
        # Shared point coordinate writer
        def _pts(fh):
            for k in K_rng:
                for j in J_rng:
                    for i in I_rng:
                        fh.write(f" {g.XMP[i]:.16e} {g.YMP[j]:.16e} {g.ZMP[k]:.16e}\n")

        if params['ILD2'] == 1:
            fpath = os.path.join(outdir, tname2 + '_ld2.vtk')
            with open(fpath, 'w') as fh:
                _vtk_header(fh, 'LESwHT instantaneous 3D lambda2',
                            (INUM, JNUM, KNUM), npts)
                _pts(fh)
                fh.write(f"POINT_DATA {npts}\n")
                write_vtk_scalar(fh, 'lambda2',
                                 [VLAMBDA2[i, j, k] for k in K_rng
                                  for j in J_rng for i in I_rng])
            print(f"  Wrote {fpath}")

        if params['IUVWP'] == 1:
            fpath = os.path.join(outdir, tname2 + '_uvwp.vtk')
            with open(fpath, 'w') as fh:
                _vtk_header(fh, 'LESwHT instantaneous 3D uvwp',
                            (INUM, JNUM, KNUM), npts)
                _pts(fh)
                fh.write(f"POINT_DATA {npts}\n")
                write_vtk_vector(fh, 'velocity',
                                 [UC[i,j,k] for k in K_rng for j in J_rng for i in I_rng],
                                 [VC[i,j,k] for k in K_rng for j in J_rng for i in I_rng],
                                 [WC[i,j,k] for k in K_rng for j in J_rng for i in I_rng])
                write_vtk_scalar(fh, 'pressure',
                                 [P[i,j,k] for k in K_rng for j in J_rng for i in I_rng])
                if ihtrans == 1:
                    write_vtk_scalar(fh, 'temperature',
                                     [T[i,j,k] for k in K_rng for j in J_rng for i in I_rng])
            print(f"  Wrote {fpath}")

        if params['IWXYZ'] == 1:
            fpath = os.path.join(outdir, tname2 + '_vorxyz.vtk')
            with open(fpath, 'w') as fh:
                _vtk_header(fh, 'LESwHT instantaneous 3D vorticity',
                            (INUM, JNUM, KNUM), npts)
                _pts(fh)
                fh.write(f"POINT_DATA {npts}\n")
                write_vtk_vector(fh, 'vorticity',
                                 [VORX[i,j,k] for k in K_rng for j in J_rng for i in I_rng],
                                 [VORY[i,j,k] for k in K_rng for j in J_rng for i in I_rng],
                                 [VORZ[i,j,k] for k in K_rng for j in J_rng for i in I_rng])
                write_vtk_scalar(fh, 'vort_mag',
                                 [np.sqrt(VORX[i,j,k]**2 + VORY[i,j,k]**2 + VORZ[i,j,k]**2)
                                  for k in K_rng for j in J_rng for i in I_rng])
            print(f"  Wrote {fpath}")

        if params['IALL'] == 1:
            fpath = os.path.join(outdir, tname2 + '_all.vtk')
            with open(fpath, 'w') as fh:
                _vtk_header(fh, 'LESwHT instantaneous 3D all',
                            (INUM, JNUM, KNUM), npts)
                _pts(fh)
                fh.write(f"POINT_DATA {npts}\n")
                write_vtk_vector(fh, 'velocity',
                                 [UC[i,j,k] for k in K_rng for j in J_rng for i in I_rng],
                                 [VC[i,j,k] for k in K_rng for j in J_rng for i in I_rng],
                                 [WC[i,j,k] for k in K_rng for j in J_rng for i in I_rng])
                write_vtk_scalar(fh, 'pressure',
                                 [P[i,j,k] for k in K_rng for j in J_rng for i in I_rng])
                write_vtk_vector(fh, 'vorticity',
                                 [VORX[i,j,k] for k in K_rng for j in J_rng for i in I_rng],
                                 [VORY[i,j,k] for k in K_rng for j in J_rng for i in I_rng],
                                 [VORZ[i,j,k] for k in K_rng for j in J_rng for i in I_rng])
                write_vtk_scalar(fh, 'vort_mag',
                                 [np.sqrt(VORX[i,j,k]**2 + VORY[i,j,k]**2 + VORZ[i,j,k]**2)
                                  for k in K_rng for j in J_rng for i in I_rng])
                write_vtk_scalar(fh, 'lambda2',
                                 [VLAMBDA2[i,j,k] for k in K_rng for j in J_rng for i in I_rng])
                if ihtrans == 1:
                    write_vtk_scalar(fh, 'temperature',
                                     [T[i,j,k] for k in K_rng for j in J_rng for i in I_rng])
            print(f"  Wrote {fpath}")


def output_3d_tec(
        g, UC, VC, WC, P, VORX, VORY, VORZ, VLAMBDA2,
        tname2, params, ihtrans, T,
        I_rng, J_rng, K_rng, INUM, JNUM, KNUM,
        outdir='../../output/post_inst',
) -> None:
    npts_tag = f'ZONE I={INUM}, J={JNUM}, K={KNUM}, F=POINT\n'

    if params['ILD2'] == 1:
        fpath = os.path.join(outdir, tname2 + '_ld2.tec')
        with open(fpath, 'w') as fh:
            fh.write('VARIABLES = "X","Y","Z","LAMBDA2"\n')
            fh.write(npts_tag)
            for k in K_rng:
                for j in J_rng:
                    for i in I_rng:
                        fh.write(f" {g.XMP[i]:.16e} {g.YMP[j]:.16e} {g.ZMP[k]:.16e}"
                                 f" {VLAMBDA2[i,j,k]:.16e}\n")
        print(f"  Wrote {fpath}")

    if params['IUVWP'] == 1:
        fpath = os.path.join(outdir, tname2 + '_uvwp.tec')
        with open(fpath, 'w') as fh:
            if ihtrans == 1:
                fh.write('VARIABLES = "X","Y","Z","U","V","W","P","T"\n')
            else:
                fh.write('VARIABLES = "X","Y","Z","U","V","W","P"\n')
            fh.write(npts_tag)
            for k in K_rng:
                for j in J_rng:
                    for i in I_rng:
                        row = (f" {g.XMP[i]:.16e} {g.YMP[j]:.16e} {g.ZMP[k]:.16e}"
                               f" {UC[i,j,k]:.16e} {VC[i,j,k]:.16e} {WC[i,j,k]:.16e}"
                               f" {P[i,j,k]:.16e}")
                        if ihtrans == 1:
                            row += f" {T[i,j,k]:.16e}"
                        fh.write(row + '\n')
        print(f"  Wrote {fpath}")

    if params['IWXYZ'] == 1:
        fpath = os.path.join(outdir, tname2 + '_vorxyz.tec')
        with open(fpath, 'w') as fh:
            fh.write('VARIABLES = "X","Y","Z","VORX","VORY","VORZ","VORMAG"\n')
            fh.write(npts_tag)
            for k in K_rng:
                for j in J_rng:
                    for i in I_rng:
                        vm = np.sqrt(VORX[i,j,k]**2 + VORY[i,j,k]**2 + VORZ[i,j,k]**2)
                        fh.write(f" {g.XMP[i]:.16e} {g.YMP[j]:.16e} {g.ZMP[k]:.16e}"
                                 f" {VORX[i,j,k]:.16e} {VORY[i,j,k]:.16e}"
                                 f" {VORZ[i,j,k]:.16e} {vm:.16e}\n")
        print(f"  Wrote {fpath}")

    if params['IALL'] == 1:
        fpath = os.path.join(outdir, tname2 + '_all.tec')
        with open(fpath, 'w') as fh:
            if ihtrans == 1:
                fh.write('VARIABLES = "X","Y","Z","U","V","W","P","VORX","VORY","VORZ","VORMAG","LAMBDA2","T"\n')
            else:
                fh.write('VARIABLES = "X","Y","Z","U","V","W","P","VORX","VORY","VORZ","VORMAG","LAMBDA2"\n')
            fh.write(npts_tag)
            for k in K_rng:
                for j in J_rng:
                    for i in I_rng:
                        vm = np.sqrt(VORX[i,j,k]**2 + VORY[i,j,k]**2 + VORZ[i,j,k]**2)
                        row = (f" {g.XMP[i]:.16e} {g.YMP[j]:.16e} {g.ZMP[k]:.16e}"
                               f" {UC[i,j,k]:.16e} {VC[i,j,k]:.16e} {WC[i,j,k]:.16e}"
                               f" {P[i,j,k]:.16e}"
                               f" {VORX[i,j,k]:.16e} {VORY[i,j,k]:.16e} {VORZ[i,j,k]:.16e}"
                               f" {vm:.16e} {VLAMBDA2[i,j,k]:.16e}")
                        if ihtrans == 1:
                            row += f" {T[i,j,k]:.16e}"
                        fh.write(row + '\n')
        print(f"  Wrote {fpath}")


# ===========================================================================
# Main
# ===========================================================================

def main():
    SEP = '=' * 68

    # --- Read input file ---
    params = read_input_inst('post_inst.input')
    IBMON    = params['IBMON']
    IHTRANS  = params['IHTRANS']
    IOUTFMT  = params['IOUTFMT']
    IUNIGRID = params['IUNIGRID']

    # --- Read periodic flags from first field file (if available) ---
    XPRDIC = 0; YPRDIC = 0; ZPRDIC = 1   # defaults
    if params['NFLD'] >= 1:
        try:
            XPRDIC, YPRDIC, ZPRDIC = read_periodic_flags(
                params['FLDNAME'][0], '../../output/field')
        except Exception as e:
            print(f"  Warning: could not read periodic flags: {e}")

    # --- Grid geometry ---
    g = read_grid('../../output/grid/grid.dat',
                  periodic_z=(ZPRDIC == 1),
                  periodic_x=(XPRDIC == 1))
    adjust_post_bounds(g, params)

    # --- Plane index lookup ---
    IP, JP, KP = find_ip_jp_kp(g, params)
    if IUNIGRID == 1:
        compute_griduni(g)

    OUTDIR = '../../output/post_inst'
    os.makedirs(OUTDIR, exist_ok=True)

    # --- Loop over field files ---
    for L, fldname in enumerate(params['FLDNAME']):
        print()
        print(SEP)
        print(f" WORKING ON  {fldname}")

        U, V, W, P, T, Re, Pr, Gr = read_inst_field(
            fldname, g, IHTRANS, XPRDIC, ZPRDIC)

        UC, VC, WC, VORX, VORY, VORZ, VLAMBDA2 = compute_vorn_lambda2(
            g, U, V, W, P, IBMON)

        tname1, tname2, tname3 = make_fhead(params['IND_FILM'], mode='inst')

        # --- X-plane cuts ---
        for n in range(params['NIP']):
            ipt = IP[n]
            print(f"   WORKING ON X= {g.XMP[ipt]:.6f}")
            ftail = make_ftail(ipt, 1)
            if IOUTFMT in (0, 2):
                output_i_vtk(ipt, g, UC, VC, WC, P, VORX, VORY, VORZ,
                             VLAMBDA2, tname1, ftail, IHTRANS, T, OUTDIR)
            if IOUTFMT in (1, 2):
                output_i_tec(ipt, g, UC, VC, WC, P, VORX, VORY, VORZ,
                             VLAMBDA2, tname1, ftail, IHTRANS, T, OUTDIR)
            if IUNIGRID == 1:
                output_i_uni_vtk(ipt, g, U, V, W, tname3, ftail, OUTDIR)

        # --- Y-plane cuts ---
        for n in range(params['NJP']):
            jpt = JP[n]
            print(f"   WORKING ON Y= {g.YMP[jpt]:.6f}")
            ftail = make_ftail(jpt, 2)
            if IOUTFMT in (0, 2):
                output_j_vtk(jpt, g, UC, VC, WC, P, VORX, VORY, VORZ,
                             VLAMBDA2, tname1, ftail, IHTRANS, T, OUTDIR)
            if IOUTFMT in (1, 2):
                output_j_tec(jpt, g, UC, VC, WC, P, VORX, VORY, VORZ,
                             VLAMBDA2, tname1, ftail, IHTRANS, T, OUTDIR)
            if IUNIGRID == 1:
                output_j_uni_vtk(jpt, g, U, W, tname3, ftail, OUTDIR)

        # --- Z-plane cuts ---
        for n in range(params['NKP']):
            kpt = KP[n]
            print(f"   WORKING ON Z= {g.ZMP[kpt]:.6f}")
            ftail = make_ftail(kpt, 3)
            if IOUTFMT in (0, 2):
                output_k_vtk(kpt, g, UC, VC, WC, P, VORX, VORY, VORZ,
                             VLAMBDA2, tname1, ftail, IHTRANS, T, OUTDIR)
            if IOUTFMT in (1, 2):
                output_k_tec(kpt, g, UC, VC, WC, P, VORX, VORY, VORZ,
                             VLAMBDA2, tname1, ftail, IHTRANS, T, OUTDIR)
            if IUNIGRID == 1:
                output_k_uni_vtk(kpt, g, U, V, tname3, ftail, OUTDIR)

        # --- Full 3D output ---
        if params['ITOT'] == 1:
            print("   WORKING ON 3D FIELD OUTPUT")
            output_3d(g, UC, VC, WC, P, VORX, VORY, VORZ, VLAMBDA2,
                      tname2, params, IHTRANS, T, IOUTFMT, OUTDIR)

        print(SEP)
        params['IND_FILM'] += 1


if __name__ == '__main__':
    main()
