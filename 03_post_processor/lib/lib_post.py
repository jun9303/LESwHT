"""
lib_post.py – Shared postprocessing utilities for LESwHT.

All numpy arrays use **1-based indexing** (index 0 is reserved for ghost /
boundary values) to mirror the Fortran convention directly.  An array that
Fortran declares as `X(0:M1)` is represented here as a NumPy array of length
M1+1, and `X[i]` == Fortran `X(I)`.

Requires: numpy, scipy
"""

from __future__ import annotations
import os
import sys
import struct
import numpy as np
from dataclasses import dataclass, field
from typing import Optional

def funcbody(x: float, y: float, z: float, t: float = 0.0) -> float:
    """Pure-Python placeholder body function (always fluid region)."""
    return 1.0


# ---------------------------------------------------------------------------
# Input-file parser
# ---------------------------------------------------------------------------

class _LineReader:
    """Read a Fortran-style input file line by line, stripping comments."""

    def __init__(self, path: str):
        with open(path) as fh:
            raw = fh.readlines()
        self._lines: list[list[str]] = []
        for ln in raw:
            # Strip inline Fortran comment
            if '!' in ln:
                ln = ln[:ln.index('!')]
            ln = ln.strip()
            if ln:
                self._lines.append(ln.split())
        self._pos = 0

    def skip(self) -> None:
        """Consume one line (dummy header)."""
        self._pos += 1

    def tokens(self) -> list[str]:
        """Return tokens of the current line and advance."""
        if self._pos >= len(self._lines):
            return []
        t = self._lines[self._pos]
        self._pos += 1
        return t

    def ints(self) -> list[int]:
        return [int(v) for v in self.tokens()]

    def floats(self) -> list[float]:
        return [float(v) for v in self.tokens()]


def read_input_inst(path: str) -> dict:
    """
    Parse `post_inst.in` and return a parameter dictionary.
    Keys mirror the Fortran variable names (upper-case).
    """
    r = _LineReader(path)

    r.skip()                                    # IBMON  IHTRANS  (header)
    ibmon, ihtrans = r.ints()
    r.skip()                                    # OUTPUT_FORMAT   (header)
    ioutfmt = r.ints()[0]

    r.skip()                                    # 2D_POST_OPTION  (header)
    r.skip()                                    # IUNIGRID        (header)
    iunigrid = r.ints()[0]

    r.skip()                                    # X_position      (header)
    nip = r.ints()[0]
    xip = [r.floats()[0] for _ in range(nip)]

    r.skip()                                    # Y_position      (header)
    njp = r.ints()[0]
    yjp = [r.floats()[0] for _ in range(njp)]

    r.skip()                                    # Z_position      (header)
    nkp = r.ints()[0]
    zkp = [r.floats()[0] for _ in range(nkp)]

    r.skip()                                    # 3D_POST_OPTION  (header)
    r.skip()                                    # ITOT            (header)
    itot = r.ints()[0]
    r.skip()                                    # ILD2 IUVWP ...  (header)
    ild2, iuvwp, iwxyz, iall = r.ints()
    r.skip()                                    # ISKIP JSKIP KSKIP (header)
    iskip, jskip, kskip = r.ints()
    r.skip()                                    # ISTART IEND header
    istart, iend = r.ints()
    r.skip()                                    # JSTART JEND header
    jstart, jend = r.ints()
    r.skip()                                    # KSTART KEND header
    kstart, kend = r.ints()
    r.skip()                                    # ANIFLD IND_FILM header
    nfld, ind_film = r.ints()
    fldname = [r.tokens()[0] for _ in range(nfld)]

    params = dict(
        IBMON=ibmon, IHTRANS=ihtrans, IOUTFMT=ioutfmt,
        IUNIGRID=iunigrid,
        NIP=nip, XIP=xip, NJP=njp, YJP=yjp, NKP=nkp, ZKP=zkp,
        ITOT=itot, ILD2=ild2, IUVWP=iuvwp, IWXYZ=iwxyz, IALL=iall,
        ISKIP=iskip, JSKIP=jskip, KSKIP=kskip,
        ISTART=istart, IEND=iend,
        JSTART=jstart, JEND=jend,
        KSTART=kstart, KEND=kend,
        NFLD=nfld, IND_FILM=ind_film,
        FLDNAME=fldname,
    )

    _print_params(params)
    return params


def read_input_avg(path: str) -> dict:
    """
    Parse `post_avg.in` and return a parameter dictionary.
    """
    r = _LineReader(path)

    r.skip()
    ibmon, ihtrans = r.ints()
    r.skip()
    ioutfmt = r.ints()[0]

    r.skip()                                    # 2D_POST_OPTION
    r.skip()                                    # IUNIGRID header
    iunigrid = r.ints()[0]

    r.skip()                                    # X_position
    nip = r.ints()[0]
    xip = [r.floats()[0] for _ in range(nip)]

    r.skip()                                    # Y_position
    njp = r.ints()[0]
    yjp = [r.floats()[0] for _ in range(njp)]

    r.skip()                                    # Z_position
    nkp = r.ints()[0]
    zkp = [r.floats()[0] for _ in range(nkp)]

    r.skip()                                    # 3D_POST_OPTION
    r.skip()                                    # ITOT header
    itot = r.ints()[0]
    r.skip()                                    # ISS/ILD2 … header
    iss, iuvwp, iwxyz, iall = r.ints()
    r.skip()                                    # ISKIP … header
    iskip, jskip, kskip = r.ints()
    r.skip()                                    # ISTART IEND header
    istart, iend = r.ints()
    r.skip()                                    # JSTART JEND header
    jstart, jend = r.ints()
    r.skip()                                    # KSTART KEND header
    kstart, kend = r.ints()
    r.skip()                                    # ANIFLD IND_FILM header
    nfld, ind_film = r.ints()
    fldname = [r.tokens()[0] for _ in range(nfld)]

    params = dict(
        IBMON=ibmon, IHTRANS=ihtrans, IOUTFMT=ioutfmt,
        IUNIGRID=iunigrid,
        NIP=nip, XIP=xip, NJP=njp, YJP=yjp, NKP=nkp, ZKP=zkp,
        ITOT=itot, ISS=iss, ILD2=iss, IUVWP=iuvwp, IWXYZ=iwxyz, IALL=iall,
        ISKIP=iskip, JSKIP=jskip, KSKIP=kskip,
        ISTART=istart, IEND=iend,
        JSTART=jstart, JEND=jend,
        KSTART=kstart, KEND=kend,
        NFLD=nfld, IND_FILM=ind_film,
        FLDNAME=fldname,
    )

    _print_params(params)
    return params


def _print_params(p: dict) -> None:
    print(f"# OF XPOSITION  = {p['NIP']:5d}")
    print(f"# OF YPOSITION  = {p['NJP']:5d}")
    print(f"# OF ZPOSITION  = {p['NKP']:5d}")
    for n, x in enumerate(p['XIP'], 1):
        print(f"X_POSITION {n:3d} : {x:13.5f}")
    for n, y in enumerate(p['YJP'], 1):
        print(f"Y_POSITION {n:3d} : {y:13.5f}")
    for n, z in enumerate(p['ZKP'], 1):
        print(f"Z_POSITION {n:3d} : {z:13.5f}")
    print(f"ITOT   ={p['ITOT']:5d}")
    crit_key = 'ISS' if 'ISS' in p else 'ILD2'
    print(f"{crit_key:<6}={p[crit_key]:5d}  IUVWP ={p['IUVWP']:5d}  IWXYZ ={p['IWXYZ']:5d}  IALL  ={p['IALL']:5d}")
    print(f"ISKIP  ={p['ISKIP']:5d}  JSKIP ={p['JSKIP']:5d}  KSKIP ={p['KSKIP']:5d}")
    print(f"ISTART ={p['ISTART']:5d}  IEND ={p['IEND']:5d}")
    print(f"JSTART ={p['JSTART']:5d}  JEND ={p['JEND']:5d}")
    print(f"KSTART ={p['KSTART']:5d}  KEND ={p['KEND']:5d}")
    print(f"NFLD   ={p['NFLD']:5d}  IND_FILM ={p['IND_FILM']:5d}")
    print(f"IUNI_GRID ={p['IUNIGRID']:5d}")
    print(f"IOUTFMT ={p['IOUTFMT']:5d}  (0:VTK, 1:TEC, 2:BOTH)")


# ---------------------------------------------------------------------------
# Grid data container
# ---------------------------------------------------------------------------

@dataclass
class GridData:
    """Holds all grid arrays.  Index conventions follow Fortran (1-based)."""
    N1: int = 0; N2: int = 0; N3: int = 0
    N1M: int = 0; N2M: int = 0; N3M: int = 0
    XL: float = 0.0; YL: float = 0.0; ZL: float = 0.0

    # Face coordinates (1-indexed, size N+1; index 0 = ghost boundary)
    X: np.ndarray = field(default_factory=lambda: np.zeros(1))
    Y: np.ndarray = field(default_factory=lambda: np.zeros(1))
    Z: np.ndarray = field(default_factory=lambda: np.zeros(1))

    # Cell-centre coordinates
    XMP: np.ndarray = field(default_factory=lambda: np.zeros(1))
    YMP: np.ndarray = field(default_factory=lambda: np.zeros(1))
    ZMP: np.ndarray = field(default_factory=lambda: np.zeros(1))

    # Adjacent-index maps
    IPV: np.ndarray = field(default_factory=lambda: np.zeros(1, dtype=np.int64))
    IMV: np.ndarray = field(default_factory=lambda: np.zeros(1, dtype=np.int64))
    JPV: np.ndarray = field(default_factory=lambda: np.zeros(1, dtype=np.int64))
    JMV: np.ndarray = field(default_factory=lambda: np.zeros(1, dtype=np.int64))
    KPV: np.ndarray = field(default_factory=lambda: np.zeros(1, dtype=np.int64))
    KMV: np.ndarray = field(default_factory=lambda: np.zeros(1, dtype=np.int64))

    # Boundary indicators (1 at boundary cell, 0 elsewhere)
    FIXIL: np.ndarray = field(default_factory=lambda: np.zeros(1))
    FIXIU: np.ndarray = field(default_factory=lambda: np.zeros(1))
    FIXJL: np.ndarray = field(default_factory=lambda: np.zeros(1))
    FIXJU: np.ndarray = field(default_factory=lambda: np.zeros(1))
    FIXKL: np.ndarray = field(default_factory=lambda: np.zeros(1))
    FIXKU: np.ndarray = field(default_factory=lambda: np.zeros(1))

    # Cell-face / cell-centre spacing and their reciprocals
    SDX: np.ndarray = field(default_factory=lambda: np.zeros(1))
    SDY: np.ndarray = field(default_factory=lambda: np.zeros(1))
    SDZ: np.ndarray = field(default_factory=lambda: np.zeros(1))
    VDX: np.ndarray = field(default_factory=lambda: np.zeros(1))
    VDY: np.ndarray = field(default_factory=lambda: np.zeros(1))
    VDZ: np.ndarray = field(default_factory=lambda: np.zeros(1))
    SSDX: np.ndarray = field(default_factory=lambda: np.zeros(1))
    SSDY: np.ndarray = field(default_factory=lambda: np.zeros(1))
    SSDZ: np.ndarray = field(default_factory=lambda: np.zeros(1))
    VVDX: np.ndarray = field(default_factory=lambda: np.zeros(1))
    VVDY: np.ndarray = field(default_factory=lambda: np.zeros(1))
    VVDZ: np.ndarray = field(default_factory=lambda: np.zeros(1))

    # Uniform-grid interpolation helpers (set by compute_griduni)
    INUMUNI: int = 0; JNUMUNI: int = 0; KNUMUNI: int = 0
    XUNI: np.ndarray = field(default_factory=lambda: np.zeros(1))
    YUNI: np.ndarray = field(default_factory=lambda: np.zeros(1))
    ZUNI: np.ndarray = field(default_factory=lambda: np.zeros(1))
    IUNI: np.ndarray = field(default_factory=lambda: np.zeros(1, dtype=np.int64))
    JUNI: np.ndarray = field(default_factory=lambda: np.zeros(1, dtype=np.int64))
    KUNI: np.ndarray = field(default_factory=lambda: np.zeros(1, dtype=np.int64))
    FACUNIX: np.ndarray = field(default_factory=lambda: np.zeros(1))
    FACUNIY: np.ndarray = field(default_factory=lambda: np.zeros(1))
    FACUNIZ: np.ndarray = field(default_factory=lambda: np.zeros(1))


# ---------------------------------------------------------------------------
# Grid reading and metric computation
# ---------------------------------------------------------------------------

_OPEN_UNIT = None   # module-level file handle for READGRID (kept simple)


def read_grid(gridfile: str, periodic_z: bool = False, periodic_x: bool = False) -> GridData:
    """
    Read `grid.dat`, compute all grid metrics, and return a GridData object.

    periodic_z / periodic_x mirror the Fortran ZPRDIC / XPRDIC flags
    (equivalent to IPZ=1 / IPX=1  in post_avg, or set via field header in
    post_inst).
    """
    g = GridData()

    with open(gridfile) as fh:
        _readgrid(fh, g)

    _indices(g)
    _indxfix(g)
    _meshes(g)
    _physpos(g)

    if periodic_z:
        _apply_z_periodicity(g)
    if periodic_x:
        _apply_x_periodicity(g)

    return g


def _readgrid(fh, g: GridData) -> None:
    tokens = fh.readline().split()
    g.N1, g.N2, g.N3 = int(tokens[0]), int(tokens[1]), int(tokens[2])
    print(f"NX={g.N1:12d}  NY={g.N2:12d}  NZ={g.N3:12d}")
    tokens = fh.readline().split()
    g.XL, g.YL, g.ZL = float(tokens[0]), float(tokens[1]), float(tokens[2])
    print(f"XL={g.XL:12.4f}  YL={g.YL:12.4f}  ZL={g.ZL:12.4f}")

    g.N1M = g.N1 - 1
    g.N2M = g.N2 - 1
    g.N3M = g.N3 - 1

    # Allocate with extra index-0 slot (ghost boundary)
    g.X = np.zeros(g.N1 + 1)
    g.Y = np.zeros(g.N2 + 1)
    g.Z = np.zeros(g.N3 + 1)

    # Read face coordinates; Fortran reads X(I), I=1..N1 etc.
    x_raw = []
    while len(x_raw) < g.N1:
        x_raw.extend(fh.readline().split())
    g.X[1:g.N1 + 1] = [float(v) for v in x_raw[:g.N1]]

    y_raw = []
    while len(y_raw) < g.N2:
        y_raw.extend(fh.readline().split())
    g.Y[1:g.N2 + 1] = [float(v) for v in y_raw[:g.N2]]

    z_raw = []
    while len(z_raw) < g.N3:
        z_raw.extend(fh.readline().split())
    g.Z[1:g.N3 + 1] = [float(v) for v in z_raw[:g.N3]]


def _indices(g: GridData) -> None:
    N1M, N2M, N3M = g.N1M, g.N2M, g.N3M
    g.IPV = np.zeros(N1M + 1, dtype=np.int64)
    g.IMV = np.zeros(N1M + 1, dtype=np.int64)
    g.JPV = np.zeros(N2M + 1, dtype=np.int64)
    g.JMV = np.zeros(N2M + 1, dtype=np.int64)
    g.KPV = np.zeros(N3M + 1, dtype=np.int64)
    g.KMV = np.zeros(N3M + 1, dtype=np.int64)
    for i in range(1, N1M + 1):
        g.IPV[i] = i + 1
        g.IMV[i] = i - 1
    for j in range(1, N2M + 1):
        g.JPV[j] = j + 1
        g.JMV[j] = j - 1
    for k in range(1, N3M + 1):
        g.KPV[k] = k + 1
        g.KMV[k] = k - 1


def _indxfix(g: GridData) -> None:
    N1M, N2M, N3M = g.N1M, g.N2M, g.N3M
    g.FIXIL = np.zeros(N1M + 1); g.FIXIU = np.zeros(N1M + 1)
    g.FIXJL = np.zeros(N2M + 1); g.FIXJU = np.zeros(N2M + 1)
    g.FIXKL = np.zeros(N3M + 1); g.FIXKU = np.zeros(N3M + 1)
    g.FIXIL[1]   = 1.0;  g.FIXIU[N1M] = 1.0
    g.FIXJL[1]   = 1.0;  g.FIXJU[N2M] = 1.0
    g.FIXKL[1]   = 1.0;  g.FIXKU[N3M] = 1.0


def _meshes(g: GridData) -> None:
    N1, N2, N3 = g.N1, g.N2, g.N3
    N1M, N2M, N3M = g.N1M, g.N2M, g.N3M

    g.SDX = np.zeros(N1 + 1); g.SDY = np.zeros(N2 + 1); g.SDZ = np.zeros(N3 + 1)
    g.VDX = np.zeros(N1 + 1); g.VDY = np.zeros(N2 + 1); g.VDZ = np.zeros(N3 + 1)
    g.SSDX = np.zeros(N1 + 1); g.SSDY = np.zeros(N2 + 1); g.SSDZ = np.zeros(N3 + 1)
    g.VVDX = np.zeros(N1 + 1); g.VVDY = np.zeros(N2 + 1); g.VVDZ = np.zeros(N3 + 1)

    # SDX[i] = X(i+1) - X(i)  for i=1..N1M
    g.SDX[1:N1M + 1] = g.X[2:N1 + 1] - g.X[1:N1M + 1]
    # VDX[i] = 0.5*(SDX[i]+SDX[i-1]) for i=2..N1M
    g.VDX[2:N1M + 1] = 0.5 * (g.SDX[2:N1M + 1] + g.SDX[1:N1M])
    g.VDX[1]  = 0.5 * g.SDX[1]
    g.VDX[N1] = 0.5 * g.SDX[N1M]
    g.SSDX[1:N1M + 1] = 1.0 / g.SDX[1:N1M + 1]
    g.VVDX[1:N1 + 1]  = 1.0 / g.VDX[1:N1 + 1]

    g.SDY[1:N2M + 1] = g.Y[2:N2 + 1] - g.Y[1:N2M + 1]
    g.VDY[2:N2M + 1] = 0.5 * (g.SDY[2:N2M + 1] + g.SDY[1:N2M])
    g.VDY[1]  = 0.5 * g.SDY[1]
    g.VDY[N2] = 0.5 * g.SDY[N2M]
    g.SSDY[1:N2M + 1] = 1.0 / g.SDY[1:N2M + 1]
    g.VVDY[1:N2 + 1]  = 1.0 / g.VDY[1:N2 + 1]

    g.SDZ[1:N3M + 1] = g.Z[2:N3 + 1] - g.Z[1:N3M + 1]
    g.VDZ[2:N3M + 1] = 0.5 * (g.SDZ[2:N3M + 1] + g.SDZ[1:N3M])
    g.VDZ[1]  = 0.5 * g.SDZ[1]
    g.VDZ[N3] = 0.5 * g.SDZ[N3M]
    g.SSDZ[1:N3M + 1] = 1.0 / g.SDZ[1:N3M + 1]
    g.VVDZ[1:N3 + 1]  = 1.0 / g.VDZ[1:N3 + 1]

    # Boundaries remain zero (SDX[0]=0, SDX[N1]=0, etc.)


def _physpos(g: GridData) -> None:
    N1, N2, N3 = g.N1, g.N2, g.N3
    N1M, N2M, N3M = g.N1M, g.N2M, g.N3M

    g.XMP = np.zeros(N1 + 1); g.YMP = np.zeros(N2 + 1); g.ZMP = np.zeros(N3 + 1)

    g.XMP[1:N1M + 1] = g.X[1:N1M + 1] + 0.5 * g.SDX[1:N1M + 1]
    g.XMP[N1] = g.X[N1];  g.XMP[0] = g.X[1]

    g.YMP[1:N2M + 1] = g.Y[1:N2M + 1] + 0.5 * g.SDY[1:N2M + 1]
    g.YMP[N2] = g.Y[N2];  g.YMP[0] = g.Y[1]

    g.ZMP[1:N3M + 1] = g.Z[1:N3M + 1] + 0.5 * g.SDZ[1:N3M + 1]
    g.ZMP[N3] = g.Z[N3];  g.ZMP[0] = g.Z[1]


def _apply_z_periodicity(g: GridData) -> None:
    N3, N3M = g.N3, g.N3M
    g.FIXKL[1]   = 0.0
    g.FIXKU[N3M] = 0.0
    g.SDZ[0]  = g.SDZ[N3M];  g.SDZ[N3]  = g.SDZ[1]
    g.SSDZ[0] = g.SSDZ[N3M]; g.SSDZ[N3] = g.SSDZ[1]
    g.VDZ[1]  = 0.5 * (g.SDZ[1] + g.SDZ[N3M])
    g.VDZ[N3] = g.VDZ[1]
    g.VVDZ[1]  = 1.0 / g.VDZ[1]
    g.VVDZ[N3] = g.VVDZ[1]
    g.KMV[1]   = N3M
    g.KPV[N3M] = 1


def _apply_x_periodicity(g: GridData) -> None:
    N1, N1M = g.N1, g.N1M
    g.FIXIL[1]   = 0.0
    g.FIXIU[N1M] = 0.0
    g.SDX[0]  = g.SDX[N1M];  g.SDX[N1]  = g.SDX[1]
    g.SSDX[0] = g.SSDX[N1M]; g.SSDX[N1] = g.SSDX[1]
    g.VDX[1]  = 0.5 * (g.SDX[1] + g.SDX[N1M])
    g.VDX[N1] = g.VDX[1]
    g.VVDX[1]  = 1.0 / g.VDX[1]
    g.VVDX[N1] = g.VVDX[1]
    g.IMV[1]   = N1M
    g.IPV[N1M] = 1


# ---------------------------------------------------------------------------
# Adjust output bounds
# ---------------------------------------------------------------------------

def adjust_post_bounds(g: GridData, params: dict) -> None:
    """Clamp and validate ISTART/IEND etc., mirroring ADJUST_POST_BOUNDS."""
    N1M, N2M, N3M = g.N1M, g.N2M, g.N3M

    for key in ('ISKIP', 'JSKIP', 'KSKIP'):
        if params[key] <= 0:
            params[key] = 1

    for key in ('ISTART', 'JSTART', 'KSTART'):
        if params[key] <= 0:
            params[key] = 1

    if params['IEND'] <= 0: params['IEND'] = N1M
    if params['JEND'] <= 0: params['JEND'] = N2M
    if params['KEND'] <= 0: params['KEND'] = N3M

    params['ISTART'] = max(1, min(params['ISTART'], N1M))
    params['IEND']   = max(1, min(params['IEND'],   N1M))
    params['JSTART'] = max(1, min(params['JSTART'], N2M))
    params['JEND']   = max(1, min(params['JEND'],   N2M))
    params['KSTART'] = max(1, min(params['KSTART'], N3M))
    params['KEND']   = max(1, min(params['KEND'],   N3M))

    for s, e in (('ISTART', 'IEND'), ('JSTART', 'JEND'), ('KSTART', 'KEND')):
        if params[s] > params[e]:
            params[s], params[e] = params[e], params[s]

    print(f"AUTO ISTART ={params['ISTART']:5d}  IEND ={params['IEND']:5d}")
    print(f"AUTO JSTART ={params['JSTART']:5d}  JEND ={params['JEND']:5d}")
    print(f"AUTO KSTART ={params['KSTART']:5d}  KEND ={params['KEND']:5d}")


# ---------------------------------------------------------------------------
# Find I/J/K plane indices
# ---------------------------------------------------------------------------

def find_ip_jp_kp(g: GridData, params: dict) -> tuple[list, list, list]:
    """Return IP, JP, KP index lists (1-based) nearest to the requested positions."""
    IP = [1] * params['NIP']
    JP = [1] * params['NJP']
    KP = [1] * params['NKP']

    for n, xi in enumerate(params['XIP']):
        for j in range(g.N1M, 0, -1):
            if xi >= g.XMP[j]:
                IP[n] = j; break

    for n, yj in enumerate(params['YJP']):
        for j in range(g.N2M, 0, -1):
            if yj >= g.YMP[j]:
                JP[n] = j; break

    for n, zk in enumerate(params['ZKP']):
        for j in range(g.N3M, 0, -1):
            if zk >= g.ZMP[j]:
                KP[n] = j; break

    return IP, JP, KP


# ---------------------------------------------------------------------------
# Uniform-grid interpolation setup
# ---------------------------------------------------------------------------

def compute_griduni(g: GridData) -> None:
    """
    Build the uniform-grid over-sampling arrays (GRIDUNI subroutine).
    Adds results directly to *g*.
    """
    N1M, N2M, N3M = g.N1M, g.N2M, g.N3M
    N1, N2, N3    = g.N1, g.N2, g.N3

    g.INUMUNI = len(range(1, N1M + 1, 2))
    g.JNUMUNI = len(range(1, N2M + 1, 2))
    g.KNUMUNI = len(range(1, N3M + 1, 2))

    DXUNI = (g.X[N1] - g.X[1]) / N1M
    DYUNI = (g.Y[N2] - g.Y[1]) / N2M
    DZUNI = (g.Z[N3] - g.Z[1]) / N3M
    print(f"DXUNI = {DXUNI:8.5f}")
    print(f"DYUNI = {DYUNI:8.5f}")
    print(f"DZUNI = {DZUNI:8.5f}")

    g.XUNI = np.zeros(N1 + 1); g.YUNI = np.zeros(N2 + 1); g.ZUNI = np.zeros(N3 + 1)
    g.XUNI[1] = g.X[1] + 0.5 * DXUNI
    g.YUNI[1] = g.Y[1] + 0.5 * DYUNI
    g.ZUNI[1] = g.Z[1] + 0.5 * DZUNI
    for i in range(2, N1 + 1): g.XUNI[i] = g.XUNI[i - 1] + DXUNI
    for j in range(2, N2 + 1): g.YUNI[j] = g.YUNI[j - 1] + DYUNI
    for k in range(2, N3 + 1): g.ZUNI[k] = g.ZUNI[k - 1] + DZUNI

    g.IUNI = np.zeros(N1M + 1, dtype=np.int64)
    g.JUNI = np.zeros(N2M + 1, dtype=np.int64)
    g.KUNI = np.zeros(N3M + 1, dtype=np.int64)
    g.FACUNIX = np.zeros(N1M + 1)
    g.FACUNIY = np.zeros(N2M + 1)
    g.FACUNIZ = np.zeros(N3M + 1)

    # --- X interpolation indices (post_avg uses XMP; post_inst uses X node) ---
    # We use the cell-centre version (matching post_avg GRIDUNI):
    for i in range(1, N1M + 1):
        xu = g.XUNI[i]
        if xu <= g.XMP[1]:
            g.IUNI[i] = 1; g.FACUNIX[i] = 0.0
        elif xu >= g.XMP[N1M]:
            g.IUNI[i] = N1M - 1; g.FACUNIX[i] = 1.0
        else:
            for j in range(N1M - 1, 0, -1):
                if xu >= g.XMP[j]:
                    g.IUNI[i] = j
                    g.FACUNIX[i] = (xu - g.XMP[j]) / (g.XMP[j + 1] - g.XMP[j])
                    break

    for i in range(1, N2M + 1):
        yu = g.YUNI[i]
        if yu <= g.YMP[1]:
            g.JUNI[i] = 1; g.FACUNIY[i] = 0.0
        elif yu >= g.YMP[N2M]:
            g.JUNI[i] = N2M - 1; g.FACUNIY[i] = 1.0
        else:
            for j in range(N2M - 1, 0, -1):
                if yu >= g.YMP[j]:
                    g.JUNI[i] = j
                    g.FACUNIY[i] = (yu - g.YMP[j]) / (g.YMP[j + 1] - g.YMP[j])
                    break

    for i in range(1, N3M + 1):
        zu = g.ZUNI[i]
        if zu <= g.ZMP[1]:
            g.KUNI[i] = 1; g.FACUNIZ[i] = 0.0
        elif zu >= g.ZMP[N3M]:
            g.KUNI[i] = N3M - 1; g.FACUNIZ[i] = 1.0
        else:
            for j in range(N3M - 1, 0, -1):
                if zu >= g.ZMP[j]:
                    g.KUNI[i] = j
                    g.FACUNIZ[i] = (zu - g.ZMP[j]) / (g.ZMP[j + 1] - g.ZMP[j])
                    break


# ---------------------------------------------------------------------------
# Output filename helpers
# ---------------------------------------------------------------------------

def make_ftail(idex: int, ll: int) -> str:
    """
    Build the plane-file tail string (e.g. `_i048.vtk`).
    ll=1 → x-plane ('_i'), ll=2 → y-plane ('_j'), ll=3 → z-plane ('_k').
    """
    i3 = idex // 100
    i2 = (idex // 10) - i3 * 10
    i1 = idex - i3 * 100 - i2 * 10
    ijk = {1: '_i', 2: '_j', 3: '_k'}[ll]
    return f"{ijk}{i3}{i2}{i1}.vtk"


def make_fhead(ind_film: int, mode: str = 'inst') -> tuple[str, str, str]:
    """
    Return (tname1, tname2, tname3) for the given film index.
    mode='inst' → prefixes 2dfm_, 3dfm_, 2dfm_uni_
    mode='avg'  → prefixes 2dav_, 3dav_, 2dav_uni_
    """
    s = f"{ind_film:05d}"
    if mode == 'inst':
        return f"2dfm_{s}", f"3dfm_{s}", f"2dfm_uni_{s}"
    else:
        return f"2dav_{s}", f"3dav_{s}", f"2dav_uni_{s}"


# ---------------------------------------------------------------------------
# VTK header helpers
# ---------------------------------------------------------------------------

def _vtk_header(fh, title: str, dims: tuple, npts: int) -> None:
    fh.write("# vtk DataFile Version 3.0\n")
    fh.write(f"{title}\n")
    fh.write("ASCII\n")
    fh.write("DATASET STRUCTURED_GRID\n")
    fh.write(f"DIMENSIONS {dims[0]} {dims[1]} {dims[2]}\n")
    fh.write(f"POINTS {npts} double\n")


def write_vtk_scalar(fh, name: str, data_flat) -> None:
    fh.write(f"SCALARS {name} double 1\n")
    fh.write("LOOKUP_TABLE default\n")
    for v in data_flat:
        fh.write(f" {v:.16e}\n")


def write_vtk_vector(fh, name: str, u_flat, v_flat, w_flat) -> None:
    fh.write(f"VECTORS {name} double\n")
    for u, v, w in zip(u_flat, v_flat, w_flat):
        fh.write(f" {u:.16e} {v:.16e} {w:.16e}\n")
