#!/usr/bin/env python3
#-*- coding:utf-8 -*-
#
# GRID GENERATION PROGRAM FOR IMMERSED BOUNDARY METHOD
# CODED BY SANGJOON LEE
# 
#

import sys
import os
from local_lib import lib_gridfunc, lib_gridvalid, lib_griddebug

# Ensure output directory exists
if not os.path.exists('../output/grid'):
    os.makedirs('../output/grid')

# Read the input file
with open('grid.input', 'r') as file:
    dummyline = file.readline()
    inputs = file.readline().split()
    N_x = int(inputs[0]); N_y = int(inputs[1]); N_z = int(inputs[2])
    L_x = float(inputs[3]); L_y = float(inputs[4]); L_z = float(inputs[5])
    dummyline = file.readline()

    gridlines = list()

    while True:
        inputs = file.readline().split()
        if not inputs or '-----' in inputs:
            break
        gridlines.append({'dir' : str(inputs[0]), 'coord_i' : float(inputs[1]), 'coord_f' : float(inputs[2]), \
                          'index_i' : int(inputs[3]), 'index_f' : int(inputs[4]), 'grid_opt' : str(inputs[5]), \
                          'factor1' : float(inputs[6]), 'factor2' : float(inputs[7])})

    dummyline = file.readline(); dummyline = file.readline()
    inputs = file.readline().split()
    debugopt = {'Midpoints' : str(inputs[0]), 'XY_plane_grid' : str(inputs[1]), 'YZ_plane_grid' : str(inputs[2]), \
                'ZX_plane_grid' : str(inputs[3]), 'dx_plot' : str(inputs[4]), 'dy_plot' : str(inputs[5]), 'dz_plot' : str(inputs[6])}

# Check the validity of the input
gridlines_x = [line for line in gridlines if line['dir'] in ['X','x']]
gridlines_y = [line for line in gridlines if line['dir'] in ['Y','y']]
gridlines_z = [line for line in gridlines if line['dir'] in ['Z','z']]

error = list()
lib_gridvalid.isValid(gridlines_x, N_x, L_x, error)
lib_gridvalid.isValid(gridlines_y, N_y, L_y, error)
lib_gridvalid.isValid(gridlines_z, N_z, L_z, error)

if error:
    for err in error:
        print(err)
    sys.exit(1)

# Generate Grid
x_coord = dict(); y_coord = dict(); z_coord = dict() 

for i, gridline in enumerate(gridlines, 1):
    if gridline['grid_opt'] in ['U','u']:
        lib_gridfunc.uniform(x_coord,y_coord,z_coord,gridline)
    elif gridline['grid_opt'] in ['G','g']:
        lib_gridfunc.geometric(x_coord,y_coord,z_coord,gridline)
    elif gridline['grid_opt'] in ['H','h']:
        lib_gridfunc.hypertan(x_coord,y_coord,z_coord,gridline)
    else:
        print(f'[Error] Unavailable grid option at line {i}. only U, G, H are recognized.')
        sys.exit(1)

# Generate the output 'grid.dat'
with open('../output/grid/grid.dat', 'w') as file:
    file.write('%23d %22d %22d\n' %(N_x, N_y, N_z))
    file.write('%23.15f %22.15f %22.15f\n' %(L_x, L_y, L_z))
    for i in range(1,N_x+1):
        file.write('%23.15f' %(x_coord[i]))
    file.write('\n')
    for i in range(1,N_y+1):
        file.write('%23.15f' %(y_coord[i]))
    file.write('\n')
    for i in range(1,N_z+1):
        file.write('%23.15f' %(z_coord[i]))

# Debugging
if debugopt['Midpoints'] == 'ON':
    lib_griddebug.midpoints(x_coord,N_x, 'x', '../output/grid/')
    lib_griddebug.midpoints(y_coord,N_y, 'y', '../output/grid/')
    lib_griddebug.midpoints(z_coord,N_z, 'z', '../output/grid/')
if debugopt['XY_plane_grid'] == 'ON':
    lib_griddebug.plane_grid(x_coord,y_coord,N_x,N_y, 'xy')
if debugopt['YZ_plane_grid'] == 'ON':
    lib_griddebug.plane_grid(y_coord,z_coord,N_y,N_z, 'yz')
if debugopt['ZX_plane_grid'] == 'ON':
    lib_griddebug.plane_grid(z_coord,x_coord,N_z,N_x, 'zx')
if debugopt['dx_plot'] == 'ON':
    lib_griddebug.deltaplot(x_coord,N_x, 'x')
if debugopt['dy_plot'] == 'ON':
    lib_griddebug.deltaplot(y_coord,N_y, 'y')
if debugopt['dz_plot'] == 'ON':
    lib_griddebug.deltaplot(z_coord,N_z, 'z')