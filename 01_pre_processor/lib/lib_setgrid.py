def readcoord(file):
    nxyz = file.readline().split()
    nx = int(nxyz[0]); ny = int(nxyz[1]); nz = int(nxyz[2])
    lxyz = file.readline().split()
    lx = float(lxyz[0]); ly = float(lxyz[1]); lz = float(lxyz[2])
    Xs = file.readline().split(); Xs.append(float(Xs[-1]))
    for i in range(nx+1):
        Xs[i] = float(Xs[i])
    Ys = file.readline().split(); Ys.append(float(Ys[-1]))
    for j in range(ny+1):
        Ys[j] = float(Ys[j])
    Zs = file.readline().split(); Zs.append(float(Zs[-1]))
    for k in range(nz+1):
        Zs[k] = float(Zs[k])

    for i in range(2,nx+1):
        Xs[-i] = Xs[-(i+1)]
    for j in range(2,ny+1):
        Ys[-j] = Ys[-(j+1)]
    for k in range(2,nz+1):
        Zs[-k] = Zs[-(k+1)]

    return [nx, ny, nz, lx, ly, lz, Xs, Ys, Zs]

def setindex(gridinfo, prdic):
    cellnum = len(gridinfo) - 2
    for i in range(1,cellnum+1):
        gridinfo[i].update({'idxplus':i+1,'idxminus':i-1})
        if i == 1:
            gridinfo[i].update({'lowerfix':1, 'upperfix':0})
        elif i == cellnum:
            gridinfo[i].update({'lowerfix':0, 'upperfix':1})
        else:
            gridinfo[i].update({'lowerfix':0, 'upperfix':0})
    if prdic == 'ON':
        gridinfo[1].update({'lowerfix':0,'idxminus':cellnum})
        gridinfo[cellnum].update({'upperfix':0,'idxplus':1})

def setmesh(gridinfo, coord, prdic):
    cellnum = len(gridinfo) - 2
    if cellnum != 1:
        for i in range(1,cellnum+1):
            gridinfo[i].update({'f2fd':coord[i+1] - coord[i]}) #distance between cell surface, f2fd(i) = xyz(i+1) - xyz(i)
            if i != 1:
                gridinfo[i].update({'c2cd':0.5*(gridinfo[i]['f2fd']+gridinfo[i-1]['f2fd'])})#distance between cell center, c2cd(i) = (xyz(i) - xyz(i-1))/2

        gridinfo[1].update({'c2cd':0.5*gridinfo[1]['f2fd']})
        gridinfo[cellnum+1].update({'c2cd':0.5*gridinfo[cellnum]['f2fd']})

        for i in range(1,cellnum+1):
            gridinfo[i].update({'1/f2fd':1/gridinfo[i]['f2fd']})
            gridinfo[i].update({'1/c2cd':1/gridinfo[i]['c2cd']})
        gridinfo[cellnum+1].update({'1/c2cd':1/gridinfo[cellnum+1]['c2cd']})
    else:
        coord[2] = coord[1] + 1.
        gridinfo[1].update({'f2fd':1.0,'c2cd':1.0,'1/f2fd':1.0,'1/c2cd':1.0})
        gridinfo[2].update({'f2fd':1.0,'c2cd':1.0,'1/f2fd':1.0,'1/c2cd':1.0})

    gridinfo[0].update({'f2fd':0, 'c2cd':0, '1/c2cd':0, '1/f2fd':0})
    gridinfo[cellnum+1].update({'f2fd':0, '1/f2fd':0})

    if prdic == 'ON':
        coord[0] = coord[1] - gridinfo[cellnum]['f2fd']
        gridinfo[0].update({'f2fd':gridinfo[cellnum]['f2fd'],'1/f2fd':gridinfo[cellnum]['1/f2fd']})
        gridinfo[1].update({'c2cd':0.5*(gridinfo[1]['f2fd']+gridinfo[cellnum]['f2fd'])})
        gridinfo[cellnum+1].update({'f2fd':gridinfo[1]['f2fd'],'1/f2fd':gridinfo[1]['1/f2fd']})
        gridinfo[cellnum+1].update({'c2cd':0.5*(gridinfo[1]['f2fd']+gridinfo[cellnum]['f2fd'])})
        for i in range(1,cellnum+2):
            gridinfo[i].update({'1/c2cd':1/gridinfo[i]['c2cd']})

    for i in range(0, cellnum+2):
        gridinfo[i].update({'coord':coord[i]})

def physcenter(gridinfo, coord, prdic):
    cellnum = len(gridinfo) - 2
    for i in range(1,cellnum+1):
        gridinfo[i].update({'center':coord[i]+0.5*gridinfo[i]['f2fd']})
    gridinfo[cellnum+1].update({'center':coord[cellnum+1]})
    gridinfo[0].update({'center':coord[0]})

def setgrid(griddir, xprdic, yprdic, zprdic):
    file = open(griddir, 'r')
    [N_x, N_y, N_z, L_x, L_y, L_z, X, Y, Z] = readcoord(file)
    file.close()

    result = dict()

    Cell_x = N_x - 1; result['N_x'] = N_x; result['Cell_x'] = Cell_x
    Cell_y = N_y - 1; result['N_y'] = N_y; result['Cell_y'] = Cell_y
    Cell_z = N_z - 1; result['N_z'] = N_z; result['Cell_z'] = Cell_z
    result['L_x'] = L_x; result['L_y'] = L_y; result['L_z'] = L_z

    grid_info_x = list()
    for i in range(N_x + 1):
        grid_info_x.append(dict())
    grid_info_y = list()
    for j in range(N_y + 1):
        grid_info_y.append(dict())
    grid_info_z = list()
    for k in range(N_z + 1):
        grid_info_z.append(dict())

    setindex(grid_info_x, xprdic); setmesh(grid_info_x, X, xprdic); physcenter(grid_info_x, X, xprdic)
    setindex(grid_info_y, yprdic); setmesh(grid_info_y, Y, yprdic); physcenter(grid_info_y, Y, yprdic)
    setindex(grid_info_z, zprdic); setmesh(grid_info_z, Z, zprdic); physcenter(grid_info_z, Z, zprdic)

    result['grid_info_x'] = grid_info_x
    result['grid_info_y'] = grid_info_y
    result['grid_info_z'] = grid_info_z

    return result