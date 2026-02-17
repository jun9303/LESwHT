# Grid variables setting module
from .lib_geometry import readcoord, setindex, setmesh, physcenter
import math

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