#!/bin/env python
#-*- coding:utf-8 -*-
#
# PRE-DETERMINATION OF IMMERSED BODIES 
# CODED BY SANGJOON LEE, RESEARCHER
# FROM SNU ENERGY AND ENVIRONMENTAL FLOW LABORATORY
# REFERENCE. LICA SOURCE BY SNU-TFC
# IN-LAB USE ONLY. PLZ DO NOT DISTRIBUTE FOR COMMERCIAL USE.
#
import os, sys, math
import numpy as np
from local_lib import lib_ibmpre

sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(__file__))))
from global_lib import lib_setgrid
from global_lib import lib_ibm_body

# Read the input file
file = open('ibmpre.input', 'r')

dummyline = file.readline()
inputs = file.readline().split()
LESSGS = str(inputs[0])
XPRDIC = str(inputs[1])
YPRDIC = 'OFF'             # YPRDIC NOT AVAILABLE IN THIS CODE
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

del(dummyline, inputs)
file.close()

file = open('../output/ibmpre/ibmpre_prdic.bin', 'w')

if XPRDIC == 'ON':
    XPRDIC_ = 1
else:
    XPRDIC_ = 0
if YPRDIC == 'ON':
    YPRDIC_ = 1
else:
    YPRDIC_ = 0
if ZPRDIC == 'ON':
    ZPRDIC_ = 1
else:
    ZPRDIC_ = 0
if IBMINT == 'ON':
    IBMINT_ = 1
else:
    IBMINT_ = 0

file.write('%d %d %d %d' %(XPRDIC_, YPRDIC_, ZPRDIC_, IBMINT_))

file.close()

# Check the input
print('Checking input options ----------')
print('LESSGS = %s, X_PERIODIC = %s, Y_PERIODICPRAN= %s, Z_PERIODIC = %s, IBM_INTERPOLATION = %s\n' \
      %(LESSGS, XPRDIC, YPRDIC, ZPRDIC, IBMINT))

print('Body recognization ---------')
if BDMODE == False:
    print('Solid bodies will be recognized from the funcbody module.')
else:
    # print('Solid bodies will be extracted from the STL file.')
    print('Point extraction from STL is under development.... Ver.180108')
    sys.exit(1)

# Call the grid.dat file and define geometric variables for calculation
grid = lib_setgrid.setgrid('../output/grid/grid.dat', XPRDIC, YPRDIC, ZPRDIC)
# grid['N_i'] : number of gridlines in i-direction
# grid['Cell_i'] : number of cells in i-direction = grid['N-i']-1
# grid['L_i'] : computational domain length in i-direction
# grid['grid_info_i'] : grid information on i-direction
#                    [j]['f2fd'] : face-2-face distance at j-th cell
#                    [j]['1/f2fd'] : inverse of face-2-face distance
#                    [j]['c2cd'] : center-2-center distance between (j-1)-th cell and j-th cell
#                    [j]['1/c2cd'] : inverse of center-2-center distance
#                    [j]['coord'] : lower face of j-th cell or coordinate of j-th gridline
#                    [j]['center'] : center point of j-th cell
#                    [j]['idxplus'] : (index+1) = j+1
#                    [j]['idxminus'] : (index-1) = j-1
#                    [j]['upperfix'] : 1 if upper end without periodic condition
#                    [j]['lowerfix'] : 1 if lower end without periodic condition

lib_ibmpre.grid_preprocessing_data(grid)

N_x = grid['N_x']; N_y = grid['N_y']; N_z = grid['N_z']
Cell_x = grid['Cell_x']; Cell_y = grid['Cell_y']; Cell_z = grid['Cell_z']
L_x = grid['L_x']; L_y = grid['L_y']; L_z = grid['L_z']

X = list(); Y = list(); Z = list()
XM = list(); YM = list(); ZM = list()

for i in range(0,N_x+1):
    X.append(grid['grid_info_x'][i]['coord'])
    XM.append(grid['grid_info_x'][i]['center'])
for j in range(0,N_y+1):
    Y.append(grid['grid_info_y'][j]['coord'])
    YM.append(grid['grid_info_y'][j]['center'])
for k in range(0,N_z+1):
    Z.append(grid['grid_info_z'][k]['coord'])
    ZM.append(grid['grid_info_z'][k]['center'])

# Find inner body pts
X = np.array(X, order='F'); Y = np.array(Y, order='F'); Z = np.array(Z, order='F');
XM = np.array(XM, order='F'); YM = np.array(YM, order='F'); ZM = np.array(ZM, order='F');

NBODY = {'u':0, 'v':0, 'w':0}                             # Number of pts in the body defined
INOUT = {'u':np.empty((N_x+1,N_y+1,N_z+1), dtype=int, order='F'), \
         'v':np.empty((N_x+1,N_y+1,N_z+1), dtype=int, order='F'), \
         'w':np.empty((N_x+1,N_y+1,N_z+1), dtype=int, order='F')}     # Whether pts is in the body or not

[NBODY['u'], INOUT['u']] = lib_ibm_body.find_inout(X,YM,ZM,.0)
[NBODY['v'], INOUT['v']] = lib_ibm_body.find_inout(XM,Y,ZM,.0)
[NBODY['w'], INOUT['w']] = lib_ibm_body.find_inout(XM,YM,Z,.0)

print('\n# of body pts for U = %d' %(NBODY['u']))
print('# of body pts for V = %d' %(NBODY['v']))
print('# of body pts for W = %d' %(NBODY['w']))

if HTRNFR == 'ON':
    NBODY.update({ 't':0 })
    INOUT.update({ 't':np.empty((N_x+1,N_y+1,N_z+1), dtype=int, order='F')})
    [NBODY['t'], INOUT['t']] = lib_ibm_body.find_inout(XM,YM,ZM,.0)
    print('# of body pts for T = %d' %(NBODY['t']))

NINNER = {'u':0, 'v':0, 'w':0}                                # Number of inner pts
FCP = {'u':np.empty([NBODY['u'],3], dtype=int, order='F'), \
       'v':np.empty([NBODY['v'],3], dtype=int, order='F'), \
       'w':np.empty([NBODY['w'],3], dtype=int, order='F')}    # Index of the forcing pts

NINTP = {'u':0, 'v':0, 'w':0}                                 # Number of bdy pts for interpolation
INTPTYPE = {'u':np.empty([NBODY['u'],1], dtype=int, order='F'), \
            'v':np.empty([NBODY['v'],1], dtype=int, order='F'), \
            'w':np.empty([NBODY['w'],1], dtype=int, order='F')}    # Type of the intp. forcing pts
            # (0 : inner, 1 : face, 2 : edge, 3 : single volume)
INTPINDX = {'u':np.empty([NBODY['u'],3], dtype=int, order='F'), \
            'v':np.empty([NBODY['v'],3], dtype=int, order='F'), \
            'w':np.empty([NBODY['w'],3], dtype=int, order='F')}    # Intp. forcing pts direction indicator

GEOMFAC = {'u':np.empty([NBODY['u'],3,3,3], order='F'), \
           'v':np.empty([NBODY['v'],3,3,3], order='F'), \
           'w':np.empty([NBODY['w'],3,3,3], order='F')}   # Geometric factor for interpolation

if HTRNFR == 'ON':
    FCP.update({ 't':np.empty([NBODY['t'],3], dtype=int, order='F') })
    NINTP.update({ 't':0 })
    INTPTYPE.update({ 't':np.empty([NBODY['t'],1], dtype=int, order='F') })
    INTPINDX.update({ 't':np.empty([NBODY['t'],3], dtype=int, order='F') })
    GEOMFAC.update({ 't':np.empty([NBODY['t'],3,3,3], order='F') })

if IBMINT == 'ON':
    # Find pts for interpolation (located at the surface of the body)
    print('\n*** INTERPOLATION MODE ON ***')
    Ufix_X = list(); Ufix_Y = list(); Ufix_Z = list()
    Lfix_X = list(); Lfix_Y = list(); Lfix_Z = list()
    for i in range(1,N_x):
        Ufix_X.append(grid['grid_info_x'][i]['upperfix'])
        Lfix_X.append(grid['grid_info_x'][i]['lowerfix'])
    for j in range(1,N_y):
        Ufix_Y.append(grid['grid_info_y'][j]['upperfix'])
        Lfix_Y.append(grid['grid_info_y'][j]['lowerfix'])
    for k in range(1,N_z):
        Ufix_Z.append(grid['grid_info_z'][k]['upperfix'])
        Lfix_Z.append(grid['grid_info_z'][k]['lowerfix'])

    [NINTP['u'], NINNER['u'], FCP['u'], INTPTYPE['u'], INTPINDX['u']] = \
         lib_ibm_body.findbdy_intp(1, NBODY['u'], INOUT['u'], Ufix_X, Ufix_Y, Ufix_Z, Lfix_X, Lfix_Y, Lfix_Z)
    [NINTP['v'], NINNER['v'], FCP['v'], INTPTYPE['v'], INTPINDX['v']] = \
         lib_ibm_body.findbdy_intp(2, NBODY['v'], INOUT['v'], Ufix_X, Ufix_Y, Ufix_Z, Lfix_X, Lfix_Y, Lfix_Z)
    [NINTP['w'], NINNER['w'], FCP['w'], INTPTYPE['w'], INTPINDX['w']] = \
         lib_ibm_body.findbdy_intp(3, NBODY['w'], INOUT['w'], Ufix_X, Ufix_Y, Ufix_Z, Lfix_X, Lfix_Y, Lfix_Z)

    INTPTYPE['u'] = INTPTYPE['u'][0:NINTP['u'],:].copy()    # Cut the unused array to NINTP
    INTPTYPE['v'] = INTPTYPE['v'][0:NINTP['v'],:].copy()    # (Unused) = (NBODY) - (NINTP)
    INTPTYPE['w'] = INTPTYPE['w'][0:NINTP['w'],:].copy()

    INTPINDX['u'] = INTPINDX['u'][0:NINTP['u'],:].copy()    # Cut the unused array to NINTP
    INTPINDX['v'] = INTPINDX['v'][0:NINTP['v'],:].copy()    # (Unused) = (NBODY) - (NINTP)
    INTPINDX['w'] = INTPINDX['w'][0:NINTP['w'],:].copy()

    GEOMFAC['u'] = GEOMFAC['u'][0:NINTP['u'],:,:,:].copy()  # Cut the unused array to NINTP
    GEOMFAC['v'] = GEOMFAC['v'][0:NINTP['v'],:,:,:].copy()  # (Unused) = (NBODY) - (NINTP)
    GEOMFAC['w'] = GEOMFAC['w'][0:NINTP['w'],:,:,:].copy()

    FCP['u'] = FCP['u'][0:NINTP['u']+NINNER['u'],:].copy()  # Cut the unused array to (NINTP+NINNER)
    FCP['v'] = FCP['v'][0:NINTP['v']+NINNER['v'],:].copy()
    FCP['w'] = FCP['w'][0:NINTP['w']+NINNER['w'],:].copy()

    if HTRNFR == 'ON':
        [NINTP['t'], NINNER['t'], FCP['t'], INTPTYPE['t'], INTPINDX['t']] = \
             lib_ibm_body.findbdy_intp(4, NBODY['t'], INOUT['t'], Ufix_X, Ufix_Y, Ufix_Z, Lfix_X, Lfix_Y, Lfix_Z)

        INTPTYPE['t'] = INTPTYPE['t'][0:NINTP['t'],:].copy()
        INTPINDX['t'] = INTPINDX['t'][0:NINTP['t'],:].copy()
        GEOMFAC['t'] = GEOMFAC['t'][0:NINTP['t'],:,:,:].copy()
        FCP['t'] = FCP['t'][0:NINTP['t']+NINNER['t'],:].copy()

    print('\n# of intp forcing pts for U = %d' %(NINTP['u']))
    print('# of intp forcing pts for V = %d' %(NINTP['v']))
    print('# of intp forcing pts for W = %d' %(NINTP['w']))
    if HTRNFR == 'ON':
        print('# of intp forcing pts for T = %d' %(NINTP['t']))
    print('\n# of inner forcing pts for U = %d' %(NINNER['u']))
    print('# of inner forcing pts for V = %d' %(NINNER['v']))
    print('# of inner forcing pts for W = %d' %(NINNER['w']))
    if HTRNFR == 'ON':
        print('# of inner forcing pts for T = %d' %(NINNER['t']))

    XX = lib_ibm_body.geomfac_preset(X, XM, XPRDIC)
    YY = lib_ibm_body.geomfac_preset(Y, YM, YPRDIC)
    ZZ = lib_ibm_body.geomfac_preset(Z, ZM, ZPRDIC)

    GEOMFAC['u'] = lib_ibm_body.geomfac_intp(XX[:,0], YY[:,1], ZZ[:,2], FCP['u'], INTPINDX['u'], .0)
    GEOMFAC['v'] = lib_ibm_body.geomfac_intp(XX[:,1], YY[:,2], ZZ[:,0], FCP['v'], INTPINDX['v'], .0)
    GEOMFAC['w'] = lib_ibm_body.geomfac_intp(XX[:,2], YY[:,0], ZZ[:,1], FCP['w'], INTPINDX['w'], .0)

    if HTRNFR == 'ON':
        GEOMFAC['t'] = lib_ibm_body.geomfac_intp(XX[:,2], YY[:,2], ZZ[:,2], FCP['t'], INTPINDX['t'], .0)

    print('\n*** GEOMETRIC FACTOR CALCULATION FINISHED ***')
    lib_ibmpre.GFIDebug(GEOMFAC, INTPTYPE, INTPINDX, FCP, X, Y, Z)

elif IBMINT == 'OFF':
    # interpolation scheme is off
    print('\n*** INTERPOLATION MODE OFF ***')

    [NINNER['u'], FCP['u']] = lib_ibm_body.findbdy_nointp(1, NBODY['u'], INOUT['u'])
    [NINNER['v'], FCP['v']] = lib_ibm_body.findbdy_nointp(2, NBODY['v'], INOUT['v'])
    [NINNER['w'], FCP['w']] = lib_ibm_body.findbdy_nointp(3, NBODY['w'], INOUT['w'])

    FCP['u'] = FCP['u'][0:NINNER['u'],:].copy()
    FCP['v'] = FCP['v'][0:NINNER['v'],:].copy()
    FCP['w'] = FCP['w'][0:NINNER['w'],:].copy()

    print('\n# of forcing pts for U = %d' %(NINNER['u']))
    print('# of forcing pts for V = %d' %(NINNER['v']))
    print('# of forcing pts for W = %d' %(NINNER['w']))

    if HTRNFR == 'ON':
        [NINNER['t'], FCP['t']] = lib_ibm_body.findbdy_nointp(4, NBODY['t'], INOUT['t'])
        FCP['t'] = FCP['t'][0:NINNER['t'],:].copy()    
        print('# of forcing pts for T = %d' %(NINNER['t']))

else:
    print('Wrong IBMINT input(ON/OFF only). Plz check again.')
    sys.exit(1)

lib_ibmpre.ibm_preprocessing_data(NINTP,NINNER,FCP,INTPINDX,GEOMFAC)
if HTRNFR == 'ON':
	lib_ibmpre.ibm_preprocessing_data_htransfer(NINTP,NINNER,FCP,INTPINDX,GEOMFAC)
np.set_printoptions(threshold=np.nan)

if LESSGS == 'ON':

    print('\n*** SGS-LES MODE ON ***')

    NZERO = 0        # Number of pts where eddy visosity should be zero
    ISZERO = np.empty((N_x-1,N_y-1,N_z-1), dtype=int, order='F')

    [NZERO, ISZERO] = lib_ibm_body.find_zero_nu_sgs(XM,YM,ZM,0)

    print('\n# of zero-eddy-visocisy pts = %d' %(NZERO))

    lib_ibmpre.les_preprocessing_data(NZERO,ISZERO)

    if CONJGHTRANS == 'ON':
    
        print('\n*** CONJG-HTRANS MODE ON ***')
    
        CSTAR = np.empty((N_x-1,N_y-1,N_z-1), dtype=float, order='F')
        KSTAR = np.empty((N_x-1,N_y-1,N_z-1,6), dtype=float, order='F')
    
        [CSTAR, KSTAR] = lib_ibm_body.conjg_intp(CRATIO,KRATIO,XM,YM,ZM,X,Y,Z,ISZERO,0)
        # np.set_printoptions(threshold=np.nan)
        # print(KSTAR[:,:,1,5])
        lib_ibmpre.conjg_preprocessing_data(CSTAR,KSTAR)

print('\n***** PRE-PROCESSING FINISHED. *****')

# Debugging options
if debugopt['U_surf_3D'] == 'ON':
    print('\nDebug -- Printing surface forcing pts for U')
    lib_ibmpre.surf3D(X,YM,ZM,INOUT['u'], 'u')
if debugopt['V_surf_3D'] == 'ON':
    print('Debug -- Printing surface forcing pts for V')
    lib_ibmpre.surf3D(XM,Y,ZM,INOUT['v'], 'v')
if debugopt['W_surf_3D'] == 'ON':
    print('Debug -- Printing surface forcing pts for W')
    lib_ibmpre.surf3D(XM,YM,Z,INOUT['w'], 'w')

