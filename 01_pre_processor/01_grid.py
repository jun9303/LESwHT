#!/bin/env python
#-*- coding:utf-8 -*-
#
# GRID GENERATION PROGRAM FOR IMMERSED BOUNDARY METHOD
# CODED BY SANGJOON LEE, RESEARCHER
# FROM SNU ENERGY AND ENVIRONMENTAL FLOW LABORATORY
# REFERENCE. LICA SOURCE BY SNU-TFC
# IN-LAB USE ONLY. PLZ DO NOT DISTRIBUTE FOR COMMERCIAL USE.
#
import sys
from local_lib import lib_gridfunc, lib_gridvalid, lib_griddebug

# Read the input file
file = open('grid.input', 'r')

dummyline = file.readline()
inputs = file.readline().split()
N_x = int(inputs[0]); N_y = int(inputs[1]); N_z = int(inputs[2])
L_x = float(inputs[3]); L_y = float(inputs[4]); L_z = float(inputs[5])
dummyline = file.readline()

gridlines = list()

while True:
    inputs = file.readline().split()
    if '-----' in inputs:
        break
    gridlines.append({'dir' : str(inputs[0]), 'coord_i' : float(inputs[1]), 'coord_f' : float(inputs[2]), \
                      'index_i' : int(inputs[3]), 'index_f' : int(inputs[4]), 'grid_opt' : str(inputs[5]), \
                      'factor1' : float(inputs[6]), 'factor2' : float(inputs[7])})

setnumber = len(gridlines)

dummyline = file.readline(); dummyline = file.readline()
inputs = file.readline().split()
debugopt = {'Midpoints' : str(inputs[0]), 'XY_plane_grid' : str(inputs[1]), 'YZ_plane_grid' : str(inputs[2]), \
            'ZX_plane_grid' : str(inputs[3]), 'dx_plot' : str(inputs[4]), 'dy_plot' : str(inputs[5]), 'dz_plot' : str(inputs[6])}

del(dummyline, inputs)
file.close()

# Check the validity of the input
gridlines_x = list(); gridlines_y = list(); gridlines_z = list()
i = 0
for gridline in gridlines:
    i = i+1
    if gridline['dir'] in ['X','x']:
        gridlines_x.append(gridline)
    elif gridline['dir'] in ['Y','y']:
        gridlines_y.append(gridline)
    elif gridline['dir'] in ['Z','z']:
        gridlines_z.append(gridline)
    else:
        print('[Error] Unavailable grid direction at %d. check the direction(X,Y,Z) again.' %(i))
        sys.exit(1)

error = list()
lib_gridvalid.is_valid(gridlines_x, N_x, L_x, error)
lib_gridvalid.is_valid(gridlines_y, N_y, L_y, error)
lib_gridvalid.is_valid(gridlines_z, N_z, L_z, error)

if error:
    for err in error:
        print(err)
    sys.exit(1)

del(error, gridlines_x, gridlines_y, gridlines_z)

x_coord = dict(); y_coord = dict(); z_coord = dict() # grid intervals in x,y,z directions 
i = 0

# Call the grid functions from the 'gridfunc' library.
for gridline in gridlines:
    i = i+1
    if gridline['grid_opt'] in ['U','u']: # unifrom gridlines
        lib_gridfunc.uniform(x_coord,y_coord,z_coord,gridline)
    elif gridline['grid_opt'] in ['G','g']: # geometric progression gridlines
        lib_gridfunc.geometric(x_coord,y_coord,z_coord,gridline)
    elif gridline['grid_opt'] in ['H','h']: # hyperbolic tangent gridlines
        lib_gridfunc.hypertan(x_coord,y_coord,z_coord,gridline)
    else:
        print('[Error] Unavailable grid option at line %d. only U, G, H are recognized.' %(i))
        sys.exit(1)

# Generate the output 'grid.dat' at /output/grid/ folder
file = open('../output/grid/grid.dat', 'w')
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
file.close()

# Debugging
if debugopt['Midpoints'] == 'ON':
    lib_griddebug.midpoints(x_coord,N_x, 'x', '../output/grid/')
    lib_griddebug.midpoints(y_coord,N_y, 'y', '../output/grid/')
    lib_griddebug.midpoints(z_coord,N_z, 'z', '../output/grid/')
if debugopt['XY_plane_grid'] == 'ON':
    lib_griddebug.plane_grid(x_coord,y_coord,N_x,N_y, 'xy', '../output/grid/')
if debugopt['YZ_plane_grid'] == 'ON':
    lib_griddebug.plane_grid(y_coord,z_coord,N_y,N_z, 'yz', '../output/grid/')
if debugopt['ZX_plane_grid'] == 'ON':
    lib_griddebug.plane_grid(z_coord,x_coord,N_z,N_x, 'zx', '../output/grid/')
if debugopt['dx_plot'] == 'ON':
    lib_griddebug.deltaplot(x_coord,N_x, 'x', '../output/grid/')
if debugopt['dy_plot'] == 'ON':
    lib_griddebug.deltaplot(y_coord,N_y, 'y', '../output/grid/')
if debugopt['dz_plot'] == 'ON':
    lib_griddebug.deltaplot(z_coord,N_z, 'z', '../output/grid/')
