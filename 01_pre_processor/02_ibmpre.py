#!/usr/bin/env python3
#-*- coding:utf-8 -*-
#
# PRE-DETERMINATION OF IMMERSED BODIES 
# CODED BY SANGJOON LEE
# 
#

import os, sys
import numpy as np
from local_lib import lib_ibmpre

# Add global_lib to path
sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(__file__))))
from global_lib import lib_setgrid
from global_lib import lib_ibm_body # New Python Native Library

# Read the input file
with open('ibmpre.input', 'r') as file:
    dummyline = file.readline()
    inputs = file.readline().split()
    LESSGS = str(inputs[0])
    XPRDIC = str(inputs[1])
    YPRDIC = 'OFF'
    ZPRDIC = str(inputs[2])
    IBMINT = str(inputs[3])
    HTRNFR = str(inputs[4])
    dummyline = file.readline()
    inputs = file.readline().split()
    CONJGHTRANS = str(inputs[0])
    CRATIO = float(inputs[1])
    KRATIO = float(inputs[2])
    dummyline = file.readline()
    inputs = file.readline().split()
    BDMODE = int(inputs[0])
    dummyline = file.readline(); dummyline = file.readline()
    inputs = file.readline().split()
    debugopt = {'U_surf_3D' : str(inputs[0]), 'V_surf_3D' : str(inputs[1]), 'W_surf_3D' : str(inputs[2])}

# Create output folder if not exists
if not os.path.exists('../output/ibmpre'):
    os.makedirs('../output/ibmpre')

with open('../output/ibmpre/ibmpre_prdic.bin', 'w') as file:
    XPRDIC_ = 1 if XPRDIC == 'ON' else 0
    YPRDIC_ = 1 if YPRDIC == 'ON' else 0
    ZPRDIC_ = 1 if ZPRDIC == 'ON' else 0
    IBMINT_ = 1 if IBMINT == 'ON' else 0
    file.write('%d %d %d %d' %(XPRDIC_, YPRDIC_, ZPRDIC_, IBMINT_))

print('Checking input options ----------')
print('LESSGS = %s, X_PERIODIC = %s, Y_PERIODIC = %s, Z_PERIODIC = %s, IBM_INTERPOLATION = %s\n' \
      %(LESSGS, XPRDIC, YPRDIC, ZPRDIC, IBMINT))

print('Body recognition ---------')
if BDMODE != 0:
    print('Error: Only DIRECT CODING (0) is supported in this version.')
    sys.exit(1)
print('Solid bodies will be recognized from the lib_funcbody module.')

# Call the grid.dat file
grid = lib_setgrid.setgrid('../output/grid/grid.dat', XPRDIC, YPRDIC, ZPRDIC)
lib_ibmpre.grid_preprocessing_data(grid)

N_x = grid['N_x']; N_y = grid['N_y']; N_z = grid['N_z']
X = np.array([grid['grid_info_x'][i]['coord'] for i in range(N_x+1)])
Y = np.array([grid['grid_info_y'][i]['coord'] for i in range(N_y+1)])
Z = np.array([grid['grid_info_z'][i]['coord'] for i in range(N_z+1)])

XM = np.array([grid['grid_info_x'][i]['center'] for i in range(N_x+1)])
YM = np.array([grid['grid_info_y'][i]['center'] for i in range(N_y+1)])
ZM = np.array([grid['grid_info_z'][i]['center'] for i in range(N_z+1)])

NBODY = {'u':0, 'v':0, 'w':0}
INOUT = {}

# 1. FIND INOUT (UPDATED CALLS)
# Note: lib_ibm_body.find_inout returns (inout_array, nbody_count)
# Also note: Explicitly pass dimensions N_x, N_y, N_z
INOUT['u'], NBODY['u'] = lib_ibm_body.find_inout(N_x, N_y, N_z, X, YM, ZM, 0.0)
INOUT['v'], NBODY['v'] = lib_ibm_body.find_inout(N_x, N_y, N_z, XM, Y, ZM, 0.0)
INOUT['w'], NBODY['w'] = lib_ibm_body.find_inout(N_x, N_y, N_z, XM, YM, Z, 0.0)

print('\n# of body pts for U = %d' %(NBODY['u']))
print('# of body pts for V = %d' %(NBODY['v']))
print('# of body pts for W = %d' %(NBODY['w']))

if HTRNFR == 'ON':
    INOUT['t'], NBODY['t'] = lib_ibm_body.find_inout(N_x, N_y, N_z, XM, YM, ZM, 0.0)
    print('# of body pts for T = %d' %(NBODY['t']))

# Initialize Data Structures
NINNER = {'u':0, 'v':0, 'w':0}
FCP = {}
NINTP = {'u':0, 'v':0, 'w':0}
INTPTYPE = {}
INTPINDX = {}
GEOMFAC = {}

if IBMINT == 'ON':
    print('\n*** INTERPOLATION MODE ON ***')
    
    # Prepare Boundary Fix Arrays
    # Python list comprehension for speed and clarity
    Ufix_X = np.array([grid['grid_info_x'][i]['upperfix'] for i in range(N_x+1)])
    Lfix_X = np.array([grid['grid_info_x'][i]['lowerfix'] for i in range(N_x+1)])
    Ufix_Y = np.array([grid['grid_info_y'][i]['upperfix'] for i in range(N_y+1)])
    Lfix_Y = np.array([grid['grid_info_y'][i]['lowerfix'] for i in range(N_y+1)])
    Ufix_Z = np.array([grid['grid_info_z'][i]['upperfix'] for i in range(N_z+1)])
    Lfix_Z = np.array([grid['grid_info_z'][i]['lowerfix'] for i in range(N_z+1)])

    # 2. FIND BOUNDARY INTERPOLATION POINTS (UPDATED CALLS)
    # Pass all fix arrays explicitly
    NINTP['u'], NINNER['u'], FCP['u'], INTPTYPE['u'], INTPINDX['u'] = \
        lib_ibm_body.findbdy_intp(1, N_x, N_y, N_z, INOUT['u'], Ufix_X, Ufix_Y, Ufix_Z, Lfix_X, Lfix_Y, Lfix_Z)

    NINTP['v'], NINNER['v'], FCP['v'], INTPTYPE['v'], INTPINDX['v'] = \
        lib_ibm_body.findbdy_intp(2, N_x, N_y, N_z, INOUT['v'], Ufix_X, Ufix_Y, Ufix_Z, Lfix_X, Lfix_Y, Lfix_Z)

    NINTP['w'], NINNER['w'], FCP['w'], INTPTYPE['w'], INTPINDX['w'] = \
        lib_ibm_body.findbdy_intp(3, N_x, N_y, N_z, INOUT['w'], Ufix_X, Ufix_Y, Ufix_Z, Lfix_X, Lfix_Y, Lfix_Z)

    if HTRNFR == 'ON':
        NINTP['t'], NINNER['t'], FCP['t'], INTPTYPE['t'], INTPINDX['t'] = \
            lib_ibm_body.findbdy_intp(4, N_x, N_y, N_z, INOUT['t'], Ufix_X, Ufix_Y, Ufix_Z, Lfix_X, Lfix_Y, Lfix_Z)

    print('\n# of intp forcing pts for U = %d' %(NINTP['u']))
    print('# of intp forcing pts for V = %d' %(NINTP['v']))
    print('# of intp forcing pts for W = %d' %(NINTP['w']))
    print('\n# of inner forcing pts for U = %d' %(NINNER['u']))
    print('# of inner forcing pts for V = %d' %(NINNER['v']))
    print('# of inner forcing pts for W = %d' %(NINNER['w']))

    # 3. GEOMETRIC FACTORS
    XX = lib_ibm_body.geomfac_preset(N_x, X, XM, XPRDIC)
    YY = lib_ibm_body.geomfac_preset(N_y, Y, YM, YPRDIC)
    ZZ = lib_ibm_body.geomfac_preset(N_z, Z, ZM, ZPRDIC)

    # Note: Pass slices (column 0 for X/U-loc, etc)
    GEOMFAC['u'] = lib_ibm_body.geomfac_intp(N_x, N_y, N_z, XX[:,0], YY[:,1], ZZ[:,2], NINTP['u'], FCP['u'], INTPINDX['u'], 0.0)
    GEOMFAC['v'] = lib_ibm_body.geomfac_intp(N_x, N_y, N_z, XX[:,1], YY[:,2], ZZ[:,0], NINTP['v'], FCP['v'], INTPINDX['v'], 0.0)
    GEOMFAC['w'] = lib_ibm_body.geomfac_intp(N_x, N_y, N_z, XX[:,2], YY[:,0], ZZ[:,1], NINTP['w'], FCP['w'], INTPINDX['w'], 0.0)

    if HTRNFR == 'ON':
        GEOMFAC['t'] = lib_ibm_body.geomfac_intp(N_x, N_y, N_z, XX[:,2], YY[:,2], ZZ[:,2], NINTP['t'], FCP['t'], INTPINDX['t'], 0.0)

    print('\n*** GEOMETRIC FACTOR CALCULATION FINISHED ***')
    lib_ibmpre.GFIDebug(GEOMFAC, INTPTYPE, INTPINDX, FCP, X, Y, Z)

elif IBMINT == 'OFF':
    print('\n*** INTERPOLATION MODE OFF ***')
    
    # 2b. FIND BODY NO INTERPOLATION
    NINNER['u'], FCP['u'] = lib_ibm_body.findbdy_nointp(1, N_x, N_y, N_z, INOUT['u'])
    NINNER['v'], FCP['v'] = lib_ibm_body.findbdy_nointp(2, N_x, N_y, N_z, INOUT['v'])
    NINNER['w'], FCP['w'] = lib_ibm_body.findbdy_nointp(3, N_x, N_y, N_z, INOUT['w'])

    print('\n# of forcing pts for U = %d' %(NINNER['u']))
    print('# of forcing pts for V = %d' %(NINNER['v']))
    print('# of forcing pts for W = %d' %(NINNER['w']))

    if HTRNFR == 'ON':
        NINNER['t'], FCP['t'] = lib_ibm_body.findbdy_nointp(4, N_x, N_y, N_z, INOUT['t'])

# Write output files
lib_ibmpre.ibm_preprocessing_data(NINTP, NINNER, FCP, INTPINDX, GEOMFAC)
if HTRNFR == 'ON':
    lib_ibmpre.ibm_preprocessing_data_htransfer(NINTP, NINNER, FCP, INTPINDX, GEOMFAC)

if LESSGS == 'ON':
    print('\n*** SGS-LES MODE ON ***')
    
    # 4. FIND ZERO SGS
    ISZERO, NZERO = lib_ibm_body.find_zero_nu_sgs(N_x, N_y, N_z, XM, YM, ZM, 0.0)
    print('\n# of zero-eddy-viscosity pts = %d' %(NZERO))
    lib_ibmpre.les_preprocessing_data(NZERO, ISZERO)

    if CONJGHTRANS == 'ON':
        print('\n*** CONJG-HTRANS MODE ON ***')
        # 5. CONJUGATE HEAT TRANSFER
        CSTAR, KSTAR = lib_ibm_body.conjg_intp(N_x, N_y, N_z, CRATIO, KRATIO, XM, YM, ZM, X, Y, Z, ISZERO, 0.0)
        lib_ibmpre.conjg_preprocessing_data(CSTAR, KSTAR)

print('\n***** PRE-PROCESSING FINISHED. *****')

# Debugging
if debugopt['U_surf_3D'] == 'ON':
    lib_ibmpre.surf3D(X, YM, ZM, INOUT['u'], 'u')
if debugopt['V_surf_3D'] == 'ON':
    lib_ibmpre.surf3D(XM, Y, ZM, INOUT['v'], 'v')
if debugopt['W_surf_3D'] == 'ON':
    lib_ibmpre.surf3D(XM, YM, Z, INOUT['w'], 'w')