import numpy as np
from global_lib.lib_funcbody import funcbody

def find_inout(nx, ny, nz, xcoord, ycoord, zcoord, t):
    """
    Determines if grid points are inside (0) or outside (1) the body.
    Vectorized implementation of FIND_INOUT.
    """
    inout = np.ones((nx + 1, ny + 1, nz + 1), dtype=int)
    nbody = 0

    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                val = funcbody(xcoord[i], ycoord[j], zcoord[k], t)
                if val <= 1.0e-10:
                    nbody += 1
                    inout[i, j, k] = 0
                    
    return inout, nbody

def findbdy_intp(direction, nx, ny, nz, inout, ufix_x, ufix_y, ufix_z, lfix_x, lfix_y, lfix_z):
    """
    Finds forcing points near the boundary for IBM interpolation.
    Correctly accepts separate arrays for X, Y, Z boundary fixations.
    """
    istart, jstart, kstart = 0, 0, 0
    if direction == 1: istart = 1
    if direction == 2: jstart = 1
    if direction == 3: kstart = 1
    
    nintp = 0
    ninner = 0
    
    fcp_intp_list = []
    type_list = []
    indx_list = []
    fcp_inner_list = []

    for k in range(kstart, nz):
        for j in range(jstart, ny):
            for i in range(istart, nx):
                
                im, ip = i - 1, i + 1
                jm, jp = j - 1, j + 1
                km, kp = k - 1, k + 1
                
                # Corrected logic: Use specific arrays for specific indices
                term_x = abs(inout[ip,j,k] - inout[im,j,k]) * (1 - ufix_x[i]) * (1 - lfix_x[i])
                term_y = abs(inout[i,jp,k] - inout[i,jm,k]) * (1 - ufix_y[j]) * (1 - lfix_y[j])
                term_z = abs(inout[i,j,kp] - inout[i,j,km]) * (1 - ufix_z[k]) * (1 - lfix_z[k])
                
                inout_sum = int(term_x + term_y + term_z)
                
                if inout[i, j, k] == 0 and inout_sum > 0:
                    nintp += 1
                    fcp_intp_list.append([i, j, k])
                    type_list.append(inout_sum)
                    
                    idx_x = inout[ip,j,k] - inout[im,j,k]
                    idx_y = inout[i,jp,k] - inout[i,jm,k]
                    idx_z = inout[i,j,kp] - inout[i,j,km]
                    indx_list.append([idx_x, idx_y, idx_z])
                    
                elif inout[i, j, k] == 0 and inout_sum == 0:
                    ninner += 1
                    fcp_inner_list.append([i, j, k])

    fcp = np.array(fcp_intp_list + fcp_inner_list, dtype=int)
    intptype = np.array(type_list, dtype=int)
    intpindx = np.array(indx_list, dtype=int)
    
    return nintp, ninner, fcp, intptype, intpindx

def geomfac_preset(ni, i_arr, im_arr, prdic):
    """Prepares coordinate arrays including ghost cells."""
    i_adj = np.zeros((ni + 3, 3))
    
    # Map input 0:NI to output 1:NI+1
    i_adj[2:ni+2, 0] = i_arr[1:ni+1]
    i_adj[2:ni+2, 1] = im_arr[1:ni+1]
    i_adj[2:ni+2, 2] = im_arr[1:ni+1]
    
    # Boundary conditions
    i_adj[1, 1] = im_arr[0]
    i_adj[1, 2] = im_arr[0]
    
    if prdic == 'ON':
        i_adj[1, 0] = i_arr[1] - (i_arr[ni] - i_arr[ni-1])
        i_adj[0, 0] = i_arr[1] - (i_arr[ni] - i_arr[ni-2])
        i_adj[ni+2, 0] = i_arr[ni] + (i_arr[2] - i_arr[1])

        i_adj[1, 1] = i_arr[1] - 0.5 * (i_arr[ni] - i_arr[ni-1])
        i_adj[0, 1] = i_arr[1] - (i_arr[ni] - i_arr[ni-1]) - 0.5 * (i_arr[ni-1] - i_arr[ni-2])
        i_adj[ni+1, 1] = i_arr[ni] + 0.5 * (i_arr[2] - i_arr[1])
        i_adj[ni+2, 1] = i_arr[ni] + i_arr[2] - i_arr[1] + 0.5 * (i_arr[3] - i_arr[2])

        i_adj[0:ni+3, 2] = i_adj[0:ni+3, 1]
        
    return i_adj

def geomfac_intp(nx, ny, nz, xpre, ypre, zpre, nintp, fcp, intpindx, t):
    """Calculates 27-point geometric interpolation factors."""
    geomfac = np.zeros((nintp, 3, 3, 3))

    for L in range(nintp):
        ix, iy, iz = fcp[L]
        x1, y1, z1 = xpre[ix+1], ypre[iy+1], zpre[iz+1]
        idx_x, idx_y, idx_z = intpindx[L]
        x2, y2, z2 = xpre[ix+idx_x+1], ypre[iy+idx_y+1], zpre[iz+idx_z+1]
        
        xx1, yy1, zz1 = 0., 0., 0.
        xx2, yy2, zz2 = 0., 0., 0.
        x0, y0, z0 = 0., 0., 0.
        skip_to_end = False
        
        for m in range(3):
            if m == 0:
                xt1, yt1, zt1 = x1, y1, z1
                xt2, yt2, zt2 = x2, y2, z2
                ffs = funcbody(x1, y1, z1, t)
                ffe = funcbody(x2, y2, z2, t)
                if ffs * ffe > 0:
                    geomfac[L, 0, 0, 0] = 1.0
                    skip_to_end = True
                    break
            else:
                xt1, yt1, zt1 = xx1, yy1, zz1
                xt2, yt2, zt2 = xx2, yy2, zz2

            found_root = False
            ddx, ddy, ddz = xt2 - xt1, yt2 - yt1, zt2 - zt1
            for n in range(20):
                nx1 = xt1 + ddx * n / 20.0
                nx2 = xt1 + ddx * (n + 1) / 20.0
                ny1 = yt1 + ddy * n / 20.0
                ny2 = yt1 + ddy * (n + 1) / 20.0
                nz1 = zt1 + ddz * n / 20.0
                nz2 = zt1 + ddz * (n + 1) / 20.0
                
                xx1, yy1, zz1 = nx1, ny1, nz1
                xx2, yy2, zz2 = nx2, ny2, nz2
                
                ff1 = funcbody(nx1, ny1, nz1, t)
                ff2 = funcbody(nx2, ny2, nz2, t)
                if ff1 == 0.: x0, y0, z0 = nx1, ny1, nz1; found_root = True; break
                elif ff2 == 0.: x0, y0, z0 = nx2, ny2, nz2; found_root = True; break
                elif ff1 * ff2 < 0.:
                    x0 = 0.5 * (nx1 + nx2); y0 = 0.5 * (ny1 + ny2); z0 = 0.5 * (nz1 + nz2)
                    if m == 2: found_root = True
                    break
            if found_root: break

        if skip_to_end: continue

        x3 = xpre[ix + idx_x * 2 + 1]
        y3 = ypre[iy + idx_y * 2 + 1]
        z3 = zpre[iz + idx_z * 2 + 1]
        dx1, dx2, dx3 = abs(x1-x0), abs(x2-x0), abs(x3-x0)
        dy1, dy2, dy3 = abs(y1-y0), abs(y2-y0), abs(y3-y0)
        dz1, dz2, dz3 = abs(z1-z0), abs(z2-z0), abs(z3-z0)

        if idx_x == 0: a0, a1 = 1., 1.
        elif dx2 >= dx1: a0, a1 = dx2/(dx1+dx2), 1.
        else: a0, a1 = 0.5, (dx3-dx1)/(dx3-dx2)
            
        if idx_y == 0: b0, b1 = 1., 1.
        elif dy2 >= dy1: b0, b1 = dy2/(dy1+dy2), 1.
        else: b0, b1 = 0.5, (dy3-dy1)/(dy3-dy2)

        if idx_z == 0: c0, c1 = 1., 1.
        elif dz2 >= dz1: c0, c1 = dz2/(dz1+dz2), 1.
        else: c0, c1 = 0.5, (dz3-dz1)/(dz3-dz2)
            
        inv = -1. / (a0 * b0 * c0)
        geomfac[L,0,0,0] = 1./(a0*b0*c0)
        geomfac[L,0,0,1] = inv * a0 * b0 * (1.-c0) * c1
        geomfac[L,0,0,2] = inv * a0 * b0 * (1.-c0) * (1.-c1)
        geomfac[L,0,1,0] = inv * a0 * (1.-b0) * c0 * b1
        geomfac[L,0,1,1] = inv * a0 * (1.-b0) * (1.-c0) * b1 * c1
        geomfac[L,0,1,2] = inv * a0 * (1.-b0) * (1.-c0) * b1 * (1.-c1)
        geomfac[L,0,2,0] = inv * a0 * (1.-b0) * c0 * (1.-b1)
        geomfac[L,0,2,1] = inv * a0 * (1.-b0) * (1.-c0) * (1.-b1) * c1
        geomfac[L,0,2,2] = inv * a0 * (1.-b0) * (1.-c0) * (1.-b1) * (1.-c1)

        geomfac[L,1,0,0] = inv * (1.-a0) * b0 * c0 * a1
        geomfac[L,1,0,1] = inv * (1.-a0) * b0 * (1.-c0) * a1 * c1
        geomfac[L,1,0,2] = inv * (1.-a0) * b0 * (1.-c0) * a1 * (1.-c1)
        geomfac[L,1,1,0] = inv * (1.-a0) * (1.-b0) * c0 * a1 * b1
        geomfac[L,1,1,1] = inv * (1.-a0) * (1.-b0) * (1.-c0) * a1 * b1 * c1
        geomfac[L,1,1,2] = inv * (1.-a0) * (1.-b0) * (1.-c0) * a1 * b1 * (1.-c1)
        geomfac[L,1,2,0] = inv * (1.-a0) * (1.-b0) * c0 * a1 * (1.-b1)
        geomfac[L,1,2,1] = inv * (1.-a0) * (1.-b0) * (1.-c0) * a1 * (1.-b1) * c1
        geomfac[L,1,2,2] = inv * (1.-a0) * (1.-b0) * (1.-c0) * a1 * (1.-b1) * (1.-c1)

        geomfac[L,2,0,0] = inv * (1.-a0) * b0 * c0 * (1.-a1)
        geomfac[L,2,0,1] = inv * (1.-a0) * b0 * (1.-c0) * (1.-a1) * c1
        geomfac[L,2,0,2] = inv * (1.-a0) * b0 * (1.-c0) * (1.-a1) * (1.-c1)
        geomfac[L,2,1,0] = inv * (1.-a0) * (1.-b0) * c0 * (1.-a1) * b1
        geomfac[L,2,1,1] = inv * (1.-a0) * (1.-b0) * (1.-c0) * (1.-a1) * b1 * c1
        geomfac[L,2,1,2] = inv * (1.-a0) * (1.-b0) * (1.-c0) * (1.-a1) * b1 * (1.-c1)
        geomfac[L,2,2,0] = inv * (1.-a0) * (1.-b0) * c0 * (1.-a1) * (1.-b1)
        geomfac[L,2,2,1] = inv * (1.-a0) * (1.-b0) * (1.-c0) * (1.-a1) * (1.-b1) * c1
        geomfac[L,2,2,2] = inv * (1.-a0) * (1.-b0) * (1.-c0) * (1.-a1) * (1.-b1) * (1.-c1)

    return geomfac

def findbdy_nointp(direction, nx, ny, nz, inout):
    """Finds boundary points without interpolation."""
    istart, jstart, kstart = 0, 0, 0
    if direction == 1: istart = 1
    if direction == 2: jstart = 1
    if direction == 3: kstart = 1
    
    fcp_list = []
    for k in range(kstart, nz):
        for j in range(jstart, ny):
            for i in range(istart, nx):
                if inout[i, j, k] == 0:
                    fcp_list.append([i, j, k])
    
    ninner = len(fcp_list)
    # Return 0x3 array even if empty to prevent shape errors
    if ninner == 0:
        return 0, np.zeros((0,3), dtype=int)
    return ninner, np.array(fcp_list, dtype=int)

def find_zero_nu_sgs(nx, ny, nz, xm, ym, zm, t):
    """Identifies cells near the wall where SGS viscosity should be zero."""
    inout = np.ones((nx, ny, nz), dtype=int)
    nzero = 0
    
    for k in range(1, nz):
        for j in range(1, ny):
            for i in range(1, nx):
                vals = [funcbody(xm[i], ym[j], zm[k], t),
                        funcbody(xm[i-1], ym[j], zm[k], t), funcbody(xm[i+1], ym[j], zm[k], t),
                        funcbody(xm[i], ym[j-1], zm[k], t), funcbody(xm[i], ym[j+1], zm[k], t),
                        funcbody(xm[i], ym[j], zm[k-1], t), funcbody(xm[i], ym[j], zm[k+1], t)]
                if any(v <= 1.0e-10 for v in vals):
                    nzero += 1
                    inout[i, j, k] = 0
    return inout, nzero

def fluid_portion(x1, x2, y1, y2, z1, z2, t, div):
    """Modified Flood-Fill Algorithm (Lee & Hwang, 2019)."""
    dx = (x2 - x1) / div
    dy = (y2 - y1) / div
    dz = (z2 - z1) / div
    sub_status = np.zeros((div, div, div), dtype=np.int8)
    
    start_x = x1 + dx * 0.5
    start_y = y1 + dy * 0.5
    start_z = z1 + dz * 0.5
    start_is_solid = (funcbody(start_x, start_y, start_z, t) <= 1.0e-10)
    
    stack = [(0, 0, 0)]
    sub_status[0, 0, 0] = 1 
    
    while stack:
        ci, cj, ck = stack.pop()
        neighbors = [(ci+1, cj, ck), (ci-1, cj, ck), (ci, cj+1, ck), (ci, cj-1, ck), (ci, cj, ck+1), (ci, cj, ck-1)]
        for ni, nj, nk in neighbors:
            if 0 <= ni < div and 0 <= nj < div and 0 <= nk < div:
                if sub_status[ni, nj, nk] == 0:
                    nx = x1 + dx * (ni + 0.5)
                    ny = y1 + dy * (nj + 0.5)
                    nz = z1 + dz * (nk + 0.5)
                    is_solid = (funcbody(nx, ny, nz, t) <= 1.0e-10)
                    if is_solid == start_is_solid:
                        sub_status[ni, nj, nk] = 1
                        stack.append((ni, nj, nk))
                    else:
                        sub_status[ni, nj, nk] = 2

    fluid_count = 0
    if not start_is_solid:
        fluid_count = np.sum(sub_status == 1)
    else:
        fluid_count = np.sum(sub_status == 2) + np.sum(sub_status == 0)

    return float(fluid_count) / (div**3)

def conjg_intp(nx, ny, nz, cratio, kratio, xm, ym, zm, x, y, z, iszero, t):
    """Calculates CSTAR and KSTAR using div=6 (Optimal)."""
    cstar = np.zeros((nx, ny, nz))
    kstar = np.zeros((nx, ny, nz, 6))
    subc = 6
    
    for k in range(1, nz):
        for j in range(1, ny):
            for i in range(1, nx):
                aa = funcbody(xm[i], ym[j], zm[k], t)
                
                fptemp = fluid_portion(x[i], x[i+1], y[j], y[j+1], z[k], z[k+1], t, subc)
                cstar[i,j,k] = (1. - fptemp) * cratio + fptemp * 1.0
                
                # Check 6 faces
                # 1. (+x)
                fptemp = fluid_portion(xm[i], xm[i+1], y[j], y[j+1], z[k], z[k+1], t, subc)
                val = funcbody(xm[i+1], ym[j], zm[k], t)
                kstar[i,j,k,0] = (1.-fptemp)*kratio + fptemp if aa*val >= 0 else kratio/(kratio*fptemp + (1.-fptemp))
                
                # 2. (-x)
                fptemp = fluid_portion(xm[i-1], xm[i], y[j], y[j+1], z[k], z[k+1], t, subc)
                val = funcbody(xm[i-1], ym[j], zm[k], t)
                kstar[i,j,k,1] = (1.-fptemp)*kratio + fptemp if aa*val >= 0 else kratio/(kratio*fptemp + (1.-fptemp))
                
                # 3. (+y)
                fptemp = fluid_portion(x[i], x[i+1], ym[j], ym[j+1], z[k], z[k+1], t, subc)
                val = funcbody(xm[i], ym[j+1], zm[k], t)
                kstar[i,j,k,2] = (1.-fptemp)*kratio + fptemp if aa*val >= 0 else kratio/(kratio*fptemp + (1.-fptemp))

                # 4. (-y)
                fptemp = fluid_portion(x[i], x[i+1], ym[j-1], ym[j], z[k], z[k+1], t, subc)
                val = funcbody(xm[i], ym[j-1], zm[k], t)
                kstar[i,j,k,3] = (1.-fptemp)*kratio + fptemp if aa*val >= 0 else kratio/(kratio*fptemp + (1.-fptemp))

                # 5. (+z)
                fptemp = fluid_portion(x[i], x[i+1], y[j], y[j+1], zm[k], zm[k+1], t, subc)
                val = funcbody(xm[i], ym[j], zm[k+1], t)
                kstar[i,j,k,4] = (1.-fptemp)*kratio + fptemp if aa*val >= 0 else kratio/(kratio*fptemp + (1.-fptemp))

                # 6. (-z)
                fptemp = fluid_portion(x[i], x[i+1], y[j], y[j+1], zm[k-1], zm[k], t, subc)
                val = funcbody(xm[i], ym[j], zm[k-1], t)
                kstar[i,j,k,5] = (1.-fptemp)*kratio + fptemp if aa*val >= 0 else kratio/(kratio*fptemp + (1.-fptemp))

    return cstar, kstar