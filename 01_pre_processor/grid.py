#!/bin/env python
#-*- coding:utf-8 -*-
#
# GRID GENERATION PROGRAM FOR IMMERSED BOUNDARY METHOD
# CODED BY SANGJOON LEE 2026
# IN-LAB USE ONLY. PLZ DO NOT DISTRIBUTE FOR COMMERCIAL USE.
#
import math
import os
import sys
from lib import lib_gridfunc, lib_gridvalid, lib_griddebug

OUTPUT_ROOT = os.environ.get('LESWHT_OUTPUT_ROOT', '../output')
GRID_DIR = os.path.join(OUTPUT_ROOT, 'grid')


def _streamwise_spacing(topo_idx):
    if topo_idx in [0, 1]:
        return 5.0 * (math.sqrt(3.0) / 2.0)
    if topo_idx in [2, 3, 4]:
        return 2.5 * (4.0 - math.sqrt(3.0))
    if topo_idx in [5, 6]:
        return 5.0
    return 5.0


def _is_235_smooth(n):
    if n < 1:
        return False
    for p in (2, 3, 5):
        while n % p == 0:
            n //= p
    return n == 1


def _nearest_235_smooth(target, min_n=2):
    n = max(min_n, int(round(target)))
    if _is_235_smooth(n):
        return n

    up = n
    while not _is_235_smooth(up):
        up += 1

    down = n
    while down >= min_n and not _is_235_smooth(down):
        down -= 1

    if down < min_n:
        return up

    if (up - target) < (target - down):
        return up
    return down


def _read_latest_case_params():
    case_name = os.environ.get('LESWHT_CASE_NAME', '').strip()
    if not case_name:
        return None

    parts = case_name.split('_')
    if len(parts) < 4:
        return None

    try:
        topo = int(parts[0])
        depth_wu = float(parts[1])
        scale = float(parts[2])
        stretch = float(parts[3])
    except ValueError:
        return None

    return {
        'topology': topo,
        'depth_wu': depth_wu,
        'scale': scale,
        'stretch': stretch,
    }


def _solve_hypertan_gamma(ncells, seg_length, first_delta, factor1=1.0):
    # Match lib_gridfunc.hypertan discrete expression for the first interval.
    if ncells <= 1:
        raise ValueError('ncells for hypertan must be greater than 1.')
    if first_delta <= 0.0 or seg_length <= 0.0:
        raise ValueError('invalid hypertan spacing request.')

    p1 = 1.0 / float(ncells)

    def delta1(gamma):
        if gamma < 1.0e-12:
            return seg_length / float(ncells)
        num = 1.0 - math.tanh(gamma * (factor1 - p1)) / math.tanh(gamma * factor1)
        return seg_length * factor1 * num

    low = 0.0
    high = 1.0

    d_low = delta1(low)
    if abs(d_low - first_delta) < 1.0e-14:
        return 0.0

    # For factor1=1.0, delta1 decreases with gamma. Need target <= delta1(0).
    if first_delta > d_low + 1.0e-12:
        raise ValueError(
            'Requested first delta is larger than uniform spacing for chosen hypertan cell count.'
        )

    d_high = delta1(high)
    while d_high > first_delta and high < 200.0:
        high *= 2.0
        d_high = delta1(high)

    if d_high > first_delta:
        raise ValueError('Failed to bracket hypertan gamma for first-delta target.')

    for _ in range(80):
        mid = 0.5 * (low + high)
        d_mid = delta1(mid)
        if d_mid > first_delta:
            low = mid
        else:
            high = mid

    return 0.5 * (low + high)


def _hypertan_deltas(ncells, seg_length, gamma, factor1=1.0):
    vals = []
    for i in range(ncells + 1):
        prop = i / float(ncells)
        y = seg_length * factor1 * (
            1.0 - math.tanh(gamma * (factor1 - prop)) / math.tanh(gamma * factor1)
        )
        vals.append(y)
    return [vals[i + 1] - vals[i] for i in range(ncells)]


def _build_case_gridlines(params):
    re_tau = 180.0
    dy_wu = 1.0 / re_tau
    dx_target = 10.0 / re_tau

    topo = params['topology']
    depth_wu = params['depth_wu']
    scale = params['scale']
    stretch = params['stretch']

    # Streamwise domain: half upstream pattern + full middle pattern + half downstream pattern.
    # This is equivalent to [-pattern_len_x, +pattern_len_x] with total length 2*pattern_len_x.
    pattern_len_x = _streamwise_spacing(topo) * stretch * scale
    x_min = -1.0 * pattern_len_x
    x_max = 1.0 * pattern_len_x

    z_min = -2.5
    z_max = 2.5

    cell_x = _nearest_235_smooth((x_max - x_min) / dx_target, min_n=10)
    cell_z = _nearest_235_smooth((z_max - z_min) / dx_target, min_n=10)
    N_x = cell_x + 1
    N_z = cell_z + 1

    # Y domain specification:
    # add 3 bumper cells below -depth and 3 bumper cells above y=2.0,
    # each with 0.5 wall-unit spacing. Physical top wall location remains y=2.0.
    depth_cells = max(0, int(round(depth_wu)))
    y_dimple = -depth_cells * dy_wu
    y_bot = y_dimple - 3.0 * dy_wu
    y_top = 2.0 + 3.0 * dy_wu

    gridlines = []

    gridlines.append({
        'dir': 'X', 'coord_i': x_min, 'coord_f': x_max,
        'index_i': 1, 'index_f': N_x,
        'grid_opt': 'U', 'factor1': 0.0, 'factor2': 0.0,
    })

    gridlines.append({
        'dir': 'Z', 'coord_i': z_min, 'coord_f': z_max,
        'index_i': 1, 'index_f': N_z,
        'grid_opt': 'U', 'factor1': 0.0, 'factor2': 0.0,
    })

    y_idx = 1

    # Segment A: 3 uniform cells below -depth.
    gridlines.append({
        'dir': 'Y', 'coord_i': y_bot, 'coord_f': y_dimple,
        'index_i': y_idx, 'index_f': y_idx + 3,
        'grid_opt': 'U', 'factor1': 0.0, 'factor2': 0.0,
    })
    y_idx += 3

    if depth_cells > 0:
        # Segment B: uniform 0.5 wall-unit spacing from -depth to 0.
        gridlines.append({
            'dir': 'Y', 'coord_i': y_dimple, 'coord_f': 0.0,
            'index_i': y_idx, 'index_f': y_idx + depth_cells,
            'grid_opt': 'U', 'factor1': 0.0, 'factor2': 0.0,
        })
        y_idx += depth_cells

    # Segment C: symmetric hypertangent in [0,2] (factor1=0.5).
    # Constrain first and last spacing to 0.5 wall unit and allow center spacing up to ~10 wall units.
    # This applies even for depth=0 so the wall-normal grid does not collapse to uniform.
    max_mid_wu = 12.0
    max_mid_delta = max_mid_wu / re_tau
    sym_factor = 0.5

    n_hyper = None
    gamma = None
    for cand in range(40, 360):
        g = _solve_hypertan_gamma(cand, 2.0, dy_wu, factor1=sym_factor)
        deltas = _hypertan_deltas(cand, 2.0, g, factor1=sym_factor)
        if max(deltas) <= max_mid_delta:
            n_hyper = cand
            gamma = g
            break

    if n_hyper is None:
        n_hyper = 359
        gamma = _solve_hypertan_gamma(n_hyper, 2.0, dy_wu, factor1=sym_factor)

    gridlines.append({
        'dir': 'Y', 'coord_i': 0.0, 'coord_f': 2.0,
        'index_i': y_idx, 'index_f': y_idx + n_hyper,
        'grid_opt': 'H', 'factor1': sym_factor, 'factor2': gamma,
    })
    y_idx += n_hyper

    # Segment D: 3 uniform bumper cells above +2.
    gridlines.append({
        'dir': 'Y', 'coord_i': 2.0, 'coord_f': y_top,
        'index_i': y_idx, 'index_f': y_idx + 3,
        'grid_opt': 'U', 'factor1': 0.0, 'factor2': 0.0,
    })
    y_idx += 3

    N_y = y_idx

    return {
        'N_x': N_x,
        'N_y': N_y,
        'N_z': N_z,
        'L_x': x_max - x_min,
        'L_y': y_top - y_bot,
        'L_z': z_max - z_min,
        'gridlines': gridlines,
        'debugopt': {
            'Midpoints': 'ON',
            'XY_plane_grid': 'ON',
            'YZ_plane_grid': 'ON',
            'ZX_plane_grid': 'ON',
            'dx_plot': 'ON',
            'dy_plot': 'ON',
            'dz_plot': 'ON',
        }
    }


def _read_legacy_grid_input():
    with open('grid.input', 'r', encoding='utf-8') as file:
        _ = file.readline()
        inputs = file.readline().split()
        N_x = int(inputs[0]); N_y = int(inputs[1]); N_z = int(inputs[2])
        L_x = float(inputs[3]); L_y = float(inputs[4]); L_z = float(inputs[5])
        _ = file.readline()

        gridlines = []
        while True:
            inputs = file.readline().split()
            if '-----' in inputs:
                break
            gridlines.append({
                'dir': str(inputs[0]),
                'coord_i': float(inputs[1]),
                'coord_f': float(inputs[2]),
                'index_i': int(inputs[3]),
                'index_f': int(inputs[4]),
                'grid_opt': str(inputs[5]),
                'factor1': float(inputs[6]),
                'factor2': float(inputs[7]),
            })

        _ = file.readline(); _ = file.readline()
        inputs = file.readline().split()
        debugopt = {
            'Midpoints': str(inputs[0]),
            'XY_plane_grid': str(inputs[1]),
            'YZ_plane_grid': str(inputs[2]),
            'ZX_plane_grid': str(inputs[3]),
            'dx_plot': str(inputs[4]),
            'dy_plot': str(inputs[5]),
            'dz_plot': str(inputs[6]),
        }

    return {
        'N_x': N_x, 'N_y': N_y, 'N_z': N_z,
        'L_x': L_x, 'L_y': L_y, 'L_z': L_z,
        'gridlines': gridlines, 'debugopt': debugopt,
    }


def _validate_gridlines(N_x, N_y, N_z, L_x, L_y, L_z, gridlines):
    gridlines_x = []
    gridlines_y = []
    gridlines_z = []

    i = 0
    for gridline in gridlines:
        i += 1
        if gridline['dir'] in ['X', 'x']:
            gridlines_x.append(gridline)
        elif gridline['dir'] in ['Y', 'y']:
            gridlines_y.append(gridline)
        elif gridline['dir'] in ['Z', 'z']:
            gridlines_z.append(gridline)
        else:
            print('[Error] Unavailable grid direction at %d. check the direction(X,Y,Z) again.' % (i))
            sys.exit(1)

    error = []
    lib_gridvalid.is_valid(gridlines_x, N_x, L_x, error)
    lib_gridvalid.is_valid(gridlines_y, N_y, L_y, error)
    lib_gridvalid.is_valid(gridlines_z, N_z, L_z, error)

    if error:
        for err in error:
            print(err)
        sys.exit(1)


def _emit_grid(N_x, N_y, N_z, L_x, L_y, L_z, gridlines, debugopt):
    x_coord = {}
    y_coord = {}
    z_coord = {}

    i = 0
    for gridline in gridlines:
        i += 1
        if gridline['grid_opt'] in ['U', 'u']:
            lib_gridfunc.uniform(x_coord, y_coord, z_coord, gridline)
        elif gridline['grid_opt'] in ['G', 'g']:
            lib_gridfunc.geometric(x_coord, y_coord, z_coord, gridline)
        elif gridline['grid_opt'] in ['H', 'h']:
            lib_gridfunc.hypertan(x_coord, y_coord, z_coord, gridline)
        else:
            print('[Error] Unavailable grid option at line %d. only U, G, H are recognized.' % (i))
            sys.exit(1)

    os.makedirs(GRID_DIR, exist_ok=True)
    with open(os.path.join(GRID_DIR, 'grid.dat'), 'w', encoding='utf-8') as file:
        file.write('%23d %22d %22d\n' % (N_x, N_y, N_z))
        file.write('%23.15f %22.15f %22.15f\n' % (L_x, L_y, L_z))
        for i in range(1, N_x + 1):
            file.write('%23.15f' % (x_coord[i]))
        file.write('\n')
        for i in range(1, N_y + 1):
            file.write('%23.15f' % (y_coord[i]))
        file.write('\n')
        for i in range(1, N_z + 1):
            file.write('%23.15f' % (z_coord[i]))

    if debugopt['Midpoints'] == 'ON':
        lib_griddebug.midpoints(x_coord, N_x, 'x', GRID_DIR + '/')
        lib_griddebug.midpoints(y_coord, N_y, 'y', GRID_DIR + '/')
        lib_griddebug.midpoints(z_coord, N_z, 'z', GRID_DIR + '/')
    if debugopt['XY_plane_grid'] == 'ON':
        lib_griddebug.plane_grid(x_coord, y_coord, N_x, N_y, 'xy', GRID_DIR + '/')
    if debugopt['YZ_plane_grid'] == 'ON':
        lib_griddebug.plane_grid(y_coord, z_coord, N_y, N_z, 'yz', GRID_DIR + '/')
    if debugopt['ZX_plane_grid'] == 'ON':
        lib_griddebug.plane_grid(z_coord, x_coord, N_z, N_x, 'zx', GRID_DIR + '/')
    if debugopt['dx_plot'] == 'ON':
        lib_griddebug.deltaplot(x_coord, N_x, 'x', GRID_DIR + '/')
    if debugopt['dy_plot'] == 'ON':
        lib_griddebug.deltaplot(y_coord, N_y, 'y', GRID_DIR + '/')
    if debugopt['dz_plot'] == 'ON':
        lib_griddebug.deltaplot(z_coord, N_z, 'z', GRID_DIR + '/')


def main():
    case_params = _read_latest_case_params()
    if case_params is not None:
        cfg = _build_case_gridlines(case_params)
        print('Grid mode: case-driven (from LESWHT_CASE_NAME).')
        print(
            'Case params: topo={topology}, depth={depth_wu}, scale={scale}, stretch={stretch}'.format(
                **case_params
            )
        )
    else:
        cfg = _read_legacy_grid_input()
        print('Grid mode: legacy grid.input.')

    _validate_gridlines(
        cfg['N_x'], cfg['N_y'], cfg['N_z'],
        cfg['L_x'], cfg['L_y'], cfg['L_z'],
        cfg['gridlines'],
    )

    _emit_grid(
        cfg['N_x'], cfg['N_y'], cfg['N_z'],
        cfg['L_x'], cfg['L_y'], cfg['L_z'],
        cfg['gridlines'], cfg['debugopt'],
    )


if __name__ == '__main__':
    main()
