"""
post_avg.py — Time-averaged field postprocessor for LESwHT.

Python port of 04_post_avg_processor/post_avg.f90.

Usage (from the 03_post_processor/avg/ directory):
    python post_avg.py

Input:  post_avg.input
    ../../output/grid/grid.dat
    ../../output/field_avg/<fldname>   (Fortran unformatted binary)

Output: ../../output/post_avg/<tname>     (VTK and/or Tecplot ASCII)

Requires: numpy and lib_post.py in the parent directory.
"""

from __future__ import annotations
import os
import sys
import struct
import numpy as np

_THIS_DIR = os.path.dirname(__file__)
_POST_ROOT = os.path.abspath(os.path.join(_THIS_DIR, '..'))
_LIB_DIR = os.path.join(_POST_ROOT, 'lib')
sys.path.insert(0, _LIB_DIR)
sys.path.insert(0, _POST_ROOT)
from lib_post import (
    GridData,
    read_input_avg, read_grid, adjust_post_bounds, find_ip_jp_kp,
    compute_griduni, make_ftail, make_fhead,
    _vtk_header, write_vtk_scalar, write_vtk_vector,
)


# ===========================================================================
# Fortran unformatted binary helpers (same helper as post_inst.py)
# ===========================================================================

def _read_frec(fh) -> bytes:
    raw = fh.read(4)
    if len(raw) < 4:
        raise EOFError("Unexpected end of file")
    nbytes = struct.unpack('<i', raw)[0]
    data   = fh.read(nbytes)
    fh.read(4)   # trailing marker (not validated for performance)
    return data


def _flt(raw: bytes) -> np.ndarray:
    return np.frombuffer(raw, dtype=np.float64)


def _decode_avg_record2(raw: bytes) -> tuple[float, float, float, float]:
    """
    Decode averaged-field record-2 metadata.

    Supports:
      - New encoding: [TIMEINIT(f8), TIMEEND(f8), IHISTINIT(f8), IHISTEND(f8)]
      - Legacy encoding: [TIMEINIT(f8), TIMEEND(f8), IHISTINIT(i8), IHISTEND(i8)]
    """
    if len(raw) != 32:
        raise ValueError(f"Unexpected record-2 length: {len(raw)} bytes")

    # First two values are always float64 times.
    tvals = np.frombuffer(raw[:16], dtype=np.float64)
    TIMEINIT = float(tvals[0])
    TIMEEND = float(tvals[1])

    # Prefer all-float decode when it looks valid.
    tail_as_f = np.frombuffer(raw[16:], dtype=np.float64)
    ih0_f, ih1_f = float(tail_as_f[0]), float(tail_as_f[1])
    if np.isfinite(ih0_f) and np.isfinite(ih1_f) and abs(ih0_f) > 1e-100 and abs(ih1_f) > 1e-100:
        return TIMEINIT, TIMEEND, ih0_f, ih1_f

    # Fallback: legacy mixed encoding where tail is int64.
    tail_as_i = np.frombuffer(raw[16:], dtype=np.int64)
    return TIMEINIT, TIMEEND, float(tail_as_i[0]), float(tail_as_i[1])


def _blank(g: GridData) -> np.ndarray:
    return np.zeros((g.N1 + 1, g.N2 + 1, g.N3 + 1))


def _blank_m(g: GridData) -> np.ndarray:
    """Cell-centre sized zero array (indices 1..NiM)."""
    return np.zeros((g.N1 + 1, g.N2 + 1, g.N3 + 1))


# ===========================================================================
# Read averaged field file
# ===========================================================================

def read_avg_field(
        fldname: str,
        g: GridData,
        ihtrans: int,
        field_dir: str = '../../output/field_avg',
) -> dict:
    """
    Read a Fortran unformatted averaged field file.

    Returns a dict with keys:
        UAVG, VAVG, WAVG           – mean velocities  (face-staggered)
        UIUJAVG                    – list of 6 Reynolds-stress components
                                     order: uu, uv, uw, vv, vw, ww
        PAVG, P2AVG                – mean pressure and pressure variance
        VORAVG                     – list of 3 mean vorticity components
        VOR2AVG                    – list of 6 vorticity-correlation components
        SSAVG                      – mean SGS viscosity (or strain-rate proxy)
        TAVG, T2AVG                – mean temperature / variance (ihtrans=1 only)
        Re                         – Reynolds number
        TIMEINIT, TIMEEND          – averaging window boundaries
        IHISTINIT, IHISTEND        – corresponding step indices
    """
    N1M, N2M, N3M = g.N1M, g.N2M, g.N3M

    UAVG  = _blank_m(g); VAVG  = _blank_m(g); WAVG  = _blank_m(g)
    UIUJAVG = [_blank_m(g) for _ in range(6)]
    PAVG  = _blank_m(g); P2AVG = _blank_m(g)
    VORAVG  = [_blank_m(g) for _ in range(3)]
    VOR2AVG = [_blank_m(g) for _ in range(6)]
    SSAVG   = _blank_m(g)
    TAVG  = _blank_m(g); T2AVG = _blank_m(g)

    fpath = os.path.join(field_dir, fldname)
    print(f"  Opening: {fpath}")

    with open(fpath, 'rb') as fh:
        # --- Record 1: NN1(i8) NN2(i8) NN3(i8) RRE(f8) ---
        rec = _read_frec(fh)
        h1i = np.frombuffer(rec[:24], dtype=np.int64)
        h1f = np.frombuffer(rec[24:32], dtype=np.float64)
        NN1, NN2, NN3 = int(h1i[0]), int(h1i[1]), int(h1i[2])
        Re = float(h1f[0])

        # --- Record 2: TIMEINIT(f8) TIMEEND(f8) IHISTINIT(f8) IHISTEND(f8) ---
        rec = _read_frec(fh)
        TIMEINIT, TIMEEND, IHISTINIT, IHISTEND = _decode_avg_record2(rec)

        DT = TIMEEND - TIMEINIT
        if abs(DT) < 1e-30:
            print("  WARNING: averaging window DT is near zero; setting DT=1")
            DT = 1.0

        if NN1 != N1M or NN2 != N2M or NN3 != N3M:
            raise ValueError(
                f"Grid/field size mismatch: grid ({N1M},{N2M},{N3M}) "
                f"vs field ({NN1},{NN2},{NN3})"
            )

        print("----------- AVERAGED FIELD INFORMATION -----------")
        print(f"RE={Re:.6e}  TIMEINIT={TIMEINIT:.5f}  TIMEEND={TIMEEND:.5f}")
        print(f"IHISTINIT={IHISTINIT:.0f}  IHISTEND={IHISTEND:.0f}")

        def _read_cell():
            """Read cell-centre field: I=1..N1M, J=1..N2M, K=1..N3M."""
            rec = _read_frec(fh)
            a = _flt(rec).reshape(N3M, N2M, N1M, order='C')
            out = _blank_m(g)
            out[1:N1M + 1, 1:N2M + 1, 1:N3M + 1] = a.transpose(2, 1, 0)
            return out / DT

        def _read_cell_packed(ncomp: int):
            """Read packed cell-centre components from one Fortran record.

            Legacy Fortran writes arrays as ((((A(i,j,k,l),i),j),k),l), i.e.
            one unformatted record containing ncomp blocks of size N1M*N2M*N3M.
            """
            rec = _read_frec(fh)
            a = _flt(rec).reshape(ncomp, N3M, N2M, N1M, order='C')
            outs = [_blank_m(g) for _ in range(ncomp)]
            for n in range(ncomp):
                outs[n][1:N1M + 1, 1:N2M + 1, 1:N3M + 1] = a[n].transpose(2, 1, 0)
                outs[n] /= DT
            return outs

        # --- Records 3–13 (and optionally 14–15 for temperature) ---
        UAVG  = _read_cell()   # Rec 3
        VAVG  = _read_cell()   # Rec 4
        WAVG  = _read_cell()   # Rec 5

        UIUJAVG = _read_cell_packed(6)  # Rec 6: uu,uv,uw,vv,vw,ww

        PAVG   = _read_cell()    # Rec 7
        P2AVG  = _read_cell()    # Rec 8

        VORAVG = _read_cell_packed(3)   # Rec 9: wx,wy,wz

        VOR2AVG = _read_cell_packed(6)  # Rec 10: wxwx,wxwy,wxwz,wywy,wywz,wzwz

        SSAVG  = _read_cell()    # Rec 11

        if ihtrans == 1:
            TAVG  = _read_cell() # Rec 12
            T2AVG = _read_cell() # Rec 13

    return dict(
        UAVG=UAVG, VAVG=VAVG, WAVG=WAVG,
        UIUJAVG=UIUJAVG,
        PAVG=PAVG, P2AVG=P2AVG,
        VORAVG=VORAVG, VOR2AVG=VOR2AVG,
        SSAVG=SSAVG,
        TAVG=TAVG, T2AVG=T2AVG,
        Re=Re,
        TIMEINIT=TIMEINIT, TIMEEND=TIMEEND,
        IHISTINIT=IHISTINIT, IHISTEND=IHISTEND,
    )


# ===========================================================================
# DATAINIT — zero IBM interior cells in averaged fields
# ===========================================================================

def datainit(g: GridData, avgs: dict, ibmon: int) -> None:
    """No-op: raw averaged fields are written without IBM/body masking."""
    return


# ===========================================================================
# Cell-centre velocity mapping for output
# ===========================================================================

def avg_cellcenter_velocities(g, avgs):
    """Averaged files store velocity on output cells; use values directly."""
    return avgs['UAVG'].copy(), avgs['VAVG'].copy(), avgs['WAVG'].copy()


# ===========================================================================
# Iteration helpers (same as post_inst)
# ===========================================================================

def _iter_jk(N2M, N3M):
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


# ===========================================================================
# VTK plane output helpers
# ===========================================================================

def _write_avg_scalars_vtk(
        fh, ii, jj, kk, avgs, g,
        ihtrans, mode='i',
        iter_pts=None,
) -> None:
    """Write all averaged scalar/vector fields for one VTK 2D plane."""
    UIUJAVG = avgs['UIUJAVG']
    VORAVG  = avgs['VORAVG']
    VOR2AVG = avgs['VOR2AVG']
    UC, VC, WC = avgs['_UC'], avgs['_VC'], avgs['_WC']

    def _vals(arr):
        if mode == 'i':
            return [arr[ii, j, k] for j, k in iter_pts]
        elif mode == 'j':
            return [arr[i, jj, k] for i, k in iter_pts]
        else:
            return [arr[i, j, kk] for i, j in iter_pts]

    write_vtk_vector(fh, 'mean_velocity', _vals(UC), _vals(VC), _vals(WC))
    write_vtk_scalar(fh, 'pavg',   _vals(avgs['PAVG']))
    write_vtk_scalar(fh, 'p2avg',  _vals(avgs['P2AVG']))
    write_vtk_scalar(fh, 'uxux',   _vals(UIUJAVG[0]))
    write_vtk_scalar(fh, 'uxuy',   _vals(UIUJAVG[1]))
    write_vtk_scalar(fh, 'uxuz',   _vals(UIUJAVG[2]))
    write_vtk_scalar(fh, 'uyuy',   _vals(UIUJAVG[3]))
    write_vtk_scalar(fh, 'uyuz',   _vals(UIUJAVG[4]))
    write_vtk_scalar(fh, 'uzuz',   _vals(UIUJAVG[5]))
    write_vtk_vector(fh, 'mean_vorticity',
                     _vals(VORAVG[0]), _vals(VORAVG[1]), _vals(VORAVG[2]))
    write_vtk_scalar(fh, 'wxwx',   _vals(VOR2AVG[0]))
    write_vtk_scalar(fh, 'wxwy',   _vals(VOR2AVG[1]))
    write_vtk_scalar(fh, 'wxwz',   _vals(VOR2AVG[2]))
    write_vtk_scalar(fh, 'wywy',   _vals(VOR2AVG[3]))
    write_vtk_scalar(fh, 'wywz',   _vals(VOR2AVG[4]))
    write_vtk_scalar(fh, 'wzwz',   _vals(VOR2AVG[5]))
    write_vtk_scalar(fh, 'ssavg',  _vals(avgs['SSAVG']))
    if ihtrans == 1:
        write_vtk_scalar(fh, 'tavg',  _vals(avgs['TAVG']))
        write_vtk_scalar(fh, 't2avg', _vals(avgs['T2AVG']))


def output_i_vtk(ipoint, g, avgs, tname1, ftailijk, ihtrans,
                 outdir='../../output/post_avg'):
    N2M, N3M = g.N2M, g.N3M
    ii = ipoint
    pts = list(_iter_jk(N2M, N3M))
    npts = N2M * N3M
    fpath = os.path.join(outdir, tname1 + ftailijk)
    with open(fpath, 'w') as fh:
        _vtk_header(fh, 'LESwHT averaged x-plane', (N3M, N2M, 1), npts)
        for j, k in pts:
            fh.write(f" {g.XMP[ii]:.16e} {g.YMP[j]:.16e} {g.ZMP[k]:.16e}\n")
        fh.write(f"POINT_DATA {npts}\n")
        _write_avg_scalars_vtk(fh, ii, None, None, avgs, g,
                               ihtrans, mode='i', iter_pts=pts)
    print(f"  Wrote {fpath}")


def output_j_vtk(jpoint, g, avgs, tname1, ftailijk, ihtrans,
                 outdir='../../output/post_avg'):
    N1M, N3M = g.N1M, g.N3M
    jj = jpoint
    pts = list(_iter_ik(N1M, N3M))
    npts = N1M * N3M
    fpath = os.path.join(outdir, tname1 + ftailijk)
    with open(fpath, 'w') as fh:
        _vtk_header(fh, 'LESwHT averaged y-plane', (N1M, N3M, 1), npts)
        for i, k in pts:
            fh.write(f" {g.XMP[i]:.16e} {g.YMP[jj]:.16e} {g.ZMP[k]:.16e}\n")
        fh.write(f"POINT_DATA {npts}\n")
        _write_avg_scalars_vtk(fh, None, jj, None, avgs, g,
                               ihtrans, mode='j', iter_pts=pts)
    print(f"  Wrote {fpath}")


def output_k_vtk(kpoint, g, avgs, tname1, ftailijk, ihtrans,
                 outdir='../../output/post_avg'):
    N1M, N2M = g.N1M, g.N2M
    kk = kpoint
    pts = list(_iter_ij(N1M, N2M))
    npts = N1M * N2M
    fpath = os.path.join(outdir, tname1 + ftailijk)
    with open(fpath, 'w') as fh:
        _vtk_header(fh, 'LESwHT averaged z-plane', (N1M, N2M, 1), npts)
        for i, j in pts:
            fh.write(f" {g.XMP[i]:.16e} {g.YMP[j]:.16e} {g.ZMP[kk]:.16e}\n")
        fh.write(f"POINT_DATA {npts}\n")
        _write_avg_scalars_vtk(fh, None, None, kk, avgs, g,
                               ihtrans, mode='k', iter_pts=pts)
    print(f"  Wrote {fpath}")


# ===========================================================================
# Tecplot ASCII output
# ===========================================================================

def _var_header_tec(ihtrans: int) -> str:
    base = ('"X","Y","Z",'
            '"U","V","W",'
            '"PAVG","P2AVG",'
            '"UXUX","UXUY","UXUZ","UYUY","UYUZ","UZUZ",'
            '"VORX","VORY","VORZ",'
            '"WXWX","WXWY","WXWZ","WYWY","WYWZ","WZWZ",'
            '"SSAVG"')
    if ihtrans == 1:
        base += ',"TAVG","T2AVG"'
    return f'VARIABLES = {base}\n'


def _row_avg_tec(xi, yi, zi, avgs, ii, jj, kk, mode, ihtrans):
    if mode == 'i':
        i_, j_, k_ = ii, jj, kk
    elif mode == 'j':
        i_, j_, k_ = ii, jj, kk
    else:
        i_, j_, k_ = ii, jj, kk
    uc = avgs['_UC'][i_, j_, k_]; vc = avgs['_VC'][i_, j_, k_]; wc = avgs['_WC'][i_, j_, k_]
    row = (f" {xi:.16e} {yi:.16e} {zi:.16e}"
           f" {uc:.16e} {vc:.16e} {wc:.16e}"
           f" {avgs['PAVG'][i_,j_,k_]:.16e} {avgs['P2AVG'][i_,j_,k_]:.16e}"
           f" {avgs['UIUJAVG'][0][i_,j_,k_]:.16e}"
           f" {avgs['UIUJAVG'][1][i_,j_,k_]:.16e}"
           f" {avgs['UIUJAVG'][2][i_,j_,k_]:.16e}"
           f" {avgs['UIUJAVG'][3][i_,j_,k_]:.16e}"
           f" {avgs['UIUJAVG'][4][i_,j_,k_]:.16e}"
           f" {avgs['UIUJAVG'][5][i_,j_,k_]:.16e}"
           f" {avgs['VORAVG'][0][i_,j_,k_]:.16e}"
           f" {avgs['VORAVG'][1][i_,j_,k_]:.16e}"
           f" {avgs['VORAVG'][2][i_,j_,k_]:.16e}"
           f" {avgs['VOR2AVG'][0][i_,j_,k_]:.16e}"
           f" {avgs['VOR2AVG'][1][i_,j_,k_]:.16e}"
           f" {avgs['VOR2AVG'][2][i_,j_,k_]:.16e}"
           f" {avgs['VOR2AVG'][3][i_,j_,k_]:.16e}"
           f" {avgs['VOR2AVG'][4][i_,j_,k_]:.16e}"
           f" {avgs['VOR2AVG'][5][i_,j_,k_]:.16e}"
           f" {avgs['SSAVG'][i_,j_,k_]:.16e}")
    if ihtrans == 1:
        row += f" {avgs['TAVG'][i_,j_,k_]:.16e} {avgs['T2AVG'][i_,j_,k_]:.16e}"
    return row


def output_i_tec(ipoint, g, avgs, tname1, ftailijk, ihtrans,
                 outdir='../../output/post_avg'):
    ftail = ftailijk[:5] + '.tec'
    fpath = os.path.join(outdir, tname1 + ftail)
    N2M, N3M = g.N2M, g.N3M
    ii = ipoint
    with open(fpath, 'w') as fh:
        fh.write(_var_header_tec(ihtrans))
        fh.write(f'ZONE I={N3M}, J={N2M}, F=POINT\n')
        for j in range(1, N2M + 1):
            for k in range(1, N3M + 1):
                fh.write(_row_avg_tec(g.XMP[ii], g.YMP[j], g.ZMP[k],
                                      avgs, ii, j, k, 'i', ihtrans) + '\n')
    print(f"  Wrote {fpath}")


def output_j_tec(jpoint, g, avgs, tname1, ftailijk, ihtrans,
                 outdir='../../output/post_avg'):
    ftail = ftailijk[:5] + '.tec'
    fpath = os.path.join(outdir, tname1 + ftail)
    N1M, N3M = g.N1M, g.N3M
    jj = jpoint
    with open(fpath, 'w') as fh:
        fh.write(_var_header_tec(ihtrans))
        fh.write(f'ZONE I={N1M}, J={N3M}, F=POINT\n')
        for k in range(1, N3M + 1):
            for i in range(1, N1M + 1):
                fh.write(_row_avg_tec(g.XMP[i], g.YMP[jj], g.ZMP[k],
                                      avgs, i, jj, k, 'j', ihtrans) + '\n')
    print(f"  Wrote {fpath}")


def output_k_tec(kpoint, g, avgs, tname1, ftailijk, ihtrans,
                 outdir='../../output/post_avg'):
    ftail = ftailijk[:5] + '.tec'
    fpath = os.path.join(outdir, tname1 + ftail)
    N1M, N2M = g.N1M, g.N2M
    kk = kpoint
    with open(fpath, 'w') as fh:
        fh.write(_var_header_tec(ihtrans))
        fh.write(f'ZONE I={N1M}, J={N2M}, F=POINT\n')
        for j in range(1, N2M + 1):
            for i in range(1, N1M + 1):
                fh.write(_row_avg_tec(g.XMP[i], g.YMP[j], g.ZMP[kk],
                                      avgs, i, j, kk, 'k', ihtrans) + '\n')
    print(f"  Wrote {fpath}")


# ===========================================================================
# 3D output
# ===========================================================================

def output_3d(
        g: GridData,
        avgs: dict,
        tname2: str,
        params: dict,
        ihtrans: int,
        ioutfmt: int,
        outdir: str = '../../output/post_avg',
) -> None:
    is_, ie = params['ISTART'], params['IEND']
    js_, je = params['JSTART'], params['JEND']
    ks_, ke = params['KSTART'], params['KEND']
    ids_, jds_, kds_ = params['ISKIP'], params['JSKIP'], params['KSKIP']

    I_rng = range(is_, ie + 1, ids_)
    J_rng = range(js_, je + 1, jds_)
    K_rng = range(ks_, ke + 1, kds_)
    INUM = len(I_rng); JNUM = len(J_rng); KNUM = len(K_rng)
    npts  = INUM * JNUM * KNUM

    UC = avgs['_UC']; VC = avgs['_VC']; WC = avgs['_WC']
    UIUJAVG = avgs['UIUJAVG']
    VORAVG  = avgs['VORAVG']
    VOR2AVG = avgs['VOR2AVG']

    if ioutfmt in (1, 2):
        output_3d_tec(g, avgs, tname2, I_rng, J_rng, K_rng,
                      INUM, JNUM, KNUM, ihtrans, params, outdir)

    if ioutfmt in (0, 2):
        def _pts(fh):
            for k in K_rng:
                for j in J_rng:
                    for i in I_rng:
                        fh.write(f" {g.XMP[i]:.16e} {g.YMP[j]:.16e} {g.ZMP[k]:.16e}\n")

        if params['ISS'] == 1:
            fpath = os.path.join(outdir, tname2 + '_ss.vtk')
            with open(fpath, 'w') as fh:
                _vtk_header(fh, 'LESwHT averaged 3D ISS (from ssavg)', (INUM, JNUM, KNUM), npts)
                _pts(fh)
                fh.write(f"POINT_DATA {npts}\n")
                write_vtk_scalar(fh, 'iss',
                                 [avgs['SSAVG'][i,j,k] for k in K_rng for j in J_rng for i in I_rng])
            print(f"  Wrote {fpath}")

        if params['IUVWP'] == 1:
            fpath = os.path.join(outdir, tname2 + '_uvwp.vtk')
            with open(fpath, 'w') as fh:
                _vtk_header(fh, 'LESwHT averaged 3D uvwp', (INUM, JNUM, KNUM), npts)
                _pts(fh)
                fh.write(f"POINT_DATA {npts}\n")
                write_vtk_vector(fh, 'mean_velocity',
                                 [UC[i,j,k] for k in K_rng for j in J_rng for i in I_rng],
                                 [VC[i,j,k] for k in K_rng for j in J_rng for i in I_rng],
                                 [WC[i,j,k] for k in K_rng for j in J_rng for i in I_rng])
                for lbl, idx in [('uxux', 0), ('uxuy', 1), ('uxuz', 2),
                                  ('uyuy', 3), ('uyuz', 4), ('uzuz', 5)]:
                    write_vtk_scalar(fh, lbl,
                                     [UIUJAVG[idx][i,j,k]
                                      for k in K_rng for j in J_rng for i in I_rng])
            print(f"  Wrote {fpath}")

        if params['IWXYZ'] == 1:
            fpath = os.path.join(outdir, tname2 + '_vorticity.vtk')
            with open(fpath, 'w') as fh:
                _vtk_header(fh, 'LESwHT averaged 3D vorticity', (INUM, JNUM, KNUM), npts)
                _pts(fh)
                fh.write(f"POINT_DATA {npts}\n")
                write_vtk_vector(fh, 'mean_vorticity',
                                 [VORAVG[0][i,j,k] for k in K_rng for j in J_rng for i in I_rng],
                                 [VORAVG[1][i,j,k] for k in K_rng for j in J_rng for i in I_rng],
                                 [VORAVG[2][i,j,k] for k in K_rng for j in J_rng for i in I_rng])
                for lbl, idx in [('wxwx', 0), ('wxwy', 1), ('wxwz', 2),
                                  ('wywy', 3), ('wywz', 4), ('wzwz', 5)]:
                    write_vtk_scalar(fh, lbl,
                                     [VOR2AVG[idx][i,j,k]
                                      for k in K_rng for j in J_rng for i in I_rng])
            print(f"  Wrote {fpath}")

        if params['IALL'] == 1:
            fpath = os.path.join(outdir, tname2 + '_all.vtk')
            with open(fpath, 'w') as fh:
                _vtk_header(fh, 'LESwHT averaged 3D all', (INUM, JNUM, KNUM), npts)
                _pts(fh)
                fh.write(f"POINT_DATA {npts}\n")
                write_vtk_vector(fh, 'mean_velocity',
                                 [UC[i,j,k] for k in K_rng for j in J_rng for i in I_rng],
                                 [VC[i,j,k] for k in K_rng for j in J_rng for i in I_rng],
                                 [WC[i,j,k] for k in K_rng for j in J_rng for i in I_rng])
                write_vtk_scalar(fh, 'pavg',
                                 [avgs['PAVG'][i,j,k] for k in K_rng for j in J_rng for i in I_rng])
                for lbl, idx in [('uxux', 0), ('uxuy', 1), ('uxuz', 2),
                                  ('uyuy', 3), ('uyuz', 4), ('uzuz', 5)]:
                    write_vtk_scalar(fh, lbl,
                                     [UIUJAVG[idx][i,j,k]
                                      for k in K_rng for j in J_rng for i in I_rng])
                write_vtk_vector(fh, 'mean_vorticity',
                                 [VORAVG[0][i,j,k] for k in K_rng for j in J_rng for i in I_rng],
                                 [VORAVG[1][i,j,k] for k in K_rng for j in J_rng for i in I_rng],
                                 [VORAVG[2][i,j,k] for k in K_rng for j in J_rng for i in I_rng])
                write_vtk_scalar(fh, 'ssavg',
                                 [avgs['SSAVG'][i,j,k] for k in K_rng for j in J_rng for i in I_rng])
                if ihtrans == 1:
                    write_vtk_scalar(fh, 'tavg',
                                     [avgs['TAVG'][i,j,k] for k in K_rng for j in J_rng for i in I_rng])
                    write_vtk_scalar(fh, 't2avg',
                                     [avgs['T2AVG'][i,j,k] for k in K_rng for j in J_rng for i in I_rng])
            print(f"  Wrote {fpath}")


def output_3d_tec(
        g, avgs, tname2, I_rng, J_rng, K_rng,
        INUM, JNUM, KNUM, ihtrans, params, outdir,
) -> None:
    npts_tag = f'ZONE I={INUM}, J={JNUM}, K={KNUM}, F=POINT\n'
    UC = avgs['_UC']; VC = avgs['_VC']; WC = avgs['_WC']
    VORAVG  = avgs['VORAVG']
    UIUJAVG = avgs['UIUJAVG']
    VOR2AVG = avgs['VOR2AVG']

    if params['ISS'] == 1:
        fpath = os.path.join(outdir, tname2 + '_ss.tec')
        with open(fpath, 'w') as fh:
            hdr = '"X","Y","Z","ISS"'
            fh.write(f'VARIABLES = {hdr}\n')
            fh.write(npts_tag)
            for k in K_rng:
                for j in J_rng:
                    for i in I_rng:
                        row = (f" {g.XMP[i]:.16e} {g.YMP[j]:.16e} {g.ZMP[k]:.16e}"
                               f" {avgs['SSAVG'][i,j,k]:.16e}")
                        fh.write(row + '\n')
        print(f"  Wrote {fpath}")

    if params['IUVWP'] == 1:
        fpath = os.path.join(outdir, tname2 + '_uvwp.tec')
        with open(fpath, 'w') as fh:
            fh.write('VARIABLES = "X","Y","Z","U","V","W","UXUX","UXUY","UXUZ","UYUY","UYUZ","UZUZ"\n')
            fh.write(npts_tag)
            for k in K_rng:
                for j in J_rng:
                    for i in I_rng:
                        row = (f" {g.XMP[i]:.16e} {g.YMP[j]:.16e} {g.ZMP[k]:.16e}"
                               f" {UC[i,j,k]:.16e} {VC[i,j,k]:.16e} {WC[i,j,k]:.16e}")
                        for n in range(6):
                            row += f" {UIUJAVG[n][i,j,k]:.16e}"
                        fh.write(row + '\n')
        print(f"  Wrote {fpath}")

    if params['IWXYZ'] == 1:
        fpath = os.path.join(outdir, tname2 + '_vorticity.tec')
        with open(fpath, 'w') as fh:
            fh.write('VARIABLES = "X","Y","Z","VORX","VORY","VORZ","WXWX","WXWY","WXWZ","WYWY","WYWZ","WZWZ"\n')
            fh.write(npts_tag)
            for k in K_rng:
                for j in J_rng:
                    for i in I_rng:
                        row = f" {g.XMP[i]:.16e} {g.YMP[j]:.16e} {g.ZMP[k]:.16e}"
                        for n in range(3):
                            row += f" {VORAVG[n][i,j,k]:.16e}"
                        for n in range(6):
                            row += f" {VOR2AVG[n][i,j,k]:.16e}"
                        fh.write(row + '\n')
        print(f"  Wrote {fpath}")

    if params['IALL'] == 1:
        fpath = os.path.join(outdir, tname2 + '_all.tec')
        with open(fpath, 'w') as fh:
            fh.write(_var_header_tec(ihtrans))
            fh.write(npts_tag)
            for k in K_rng:
                for j in J_rng:
                    for i in I_rng:
                        row = _row_avg_tec(g.XMP[i], g.YMP[j], g.ZMP[k],
                                           avgs, i, j, k, 'k', ihtrans)
                        fh.write(row + '\n')
        print(f"  Wrote {fpath}")


# ===========================================================================
# Main
# ===========================================================================

def main():
    SEP = '=' * 68

    params   = read_input_avg('post_avg.input')
    IBMON    = params['IBMON']
    IHTRANS  = params['IHTRANS']
    IOUTFMT  = params['IOUTFMT']
    IUNIGRID = params['IUNIGRID']

    g = read_grid('../../output/grid/grid.dat')
    adjust_post_bounds(g, params)
    IP, JP, KP = find_ip_jp_kp(g, params)
    if IUNIGRID == 1:
        compute_griduni(g)

    OUTDIR = '../../output/post_avg'
    os.makedirs(OUTDIR, exist_ok=True)

    for L, fldname in enumerate(params['FLDNAME']):
        print()
        print(SEP)
        print(f" WORKING ON  {fldname}")

        avgs = read_avg_field(fldname, g, IHTRANS)
        datainit(g, avgs, IBMON)

        # Direct mapping from averaged velocity fields to output cell centres
        UC, VC, WC = avg_cellcenter_velocities(g, avgs)
        avgs['_UC'] = UC; avgs['_VC'] = VC; avgs['_WC'] = WC

        tname1, tname2, tname3 = make_fhead(params['IND_FILM'], mode='avg')

        for n in range(params['NIP']):
            ipt = IP[n]
            print(f"   WORKING ON X= {g.XMP[ipt]:.6f}")
            ftail = make_ftail(ipt, 1)
            if IOUTFMT in (0, 2):
                output_i_vtk(ipt, g, avgs, tname1, ftail, IHTRANS, OUTDIR)
            if IOUTFMT in (1, 2):
                output_i_tec(ipt, g, avgs, tname1, ftail, IHTRANS, OUTDIR)

        for n in range(params['NJP']):
            jpt = JP[n]
            print(f"   WORKING ON Y= {g.YMP[jpt]:.6f}")
            ftail = make_ftail(jpt, 2)
            if IOUTFMT in (0, 2):
                output_j_vtk(jpt, g, avgs, tname1, ftail, IHTRANS, OUTDIR)
            if IOUTFMT in (1, 2):
                output_j_tec(jpt, g, avgs, tname1, ftail, IHTRANS, OUTDIR)

        for n in range(params['NKP']):
            kpt = KP[n]
            print(f"   WORKING ON Z= {g.ZMP[kpt]:.6f}")
            ftail = make_ftail(kpt, 3)
            if IOUTFMT in (0, 2):
                output_k_vtk(kpt, g, avgs, tname1, ftail, IHTRANS, OUTDIR)
            if IOUTFMT in (1, 2):
                output_k_tec(kpt, g, avgs, tname1, ftail, IHTRANS, OUTDIR)

        if params['ITOT'] == 1:
            print("   WORKING ON 3D AVERAGED FIELD OUTPUT")
            output_3d(g, avgs, tname2, params, IHTRANS, IOUTFMT, OUTDIR)

        print(SEP)
        params['IND_FILM'] += 1


if __name__ == '__main__':
    main()
