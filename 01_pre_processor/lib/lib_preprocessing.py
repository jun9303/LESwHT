import numpy as np

def GFIDebug(geomfac, intptype, intpindx, fcp, x, y, z):
    gfimax = [0, 0, 0]
    gfimax_indx = [0, 0, 0]
    gfi_ltone = 0

    max_type1 = 0; type1 = 0 
    max_type2 = 0; type2 = 0
    max_type3 = 0; type3 = 0

    temp = geomfac['u'].copy(); temp = np.absolute(temp); temp[:,0,0,0] = 0
    temp2 = temp.flatten()
    if temp2.shape[0] is not 0:
        gfimax[0] = max(np.max(temp2),gfimax[0]); del(temp2)
        gfimax_indx[0]=np.where(temp >= gfimax[0])[0]
        gfi_ltone = gfi_ltone + len(np.where(temp > 1.0)[0])
    
    temp = geomfac['v'].copy(); temp = np.absolute(temp); temp[:,0,0,0] = 0
    temp2 = temp.flatten()
    if temp2.shape[0] is not 0:
        gfimax[1] = max(np.max(temp2),gfimax[1]); del(temp2)
        gfimax_indx[1]=np.where(temp >= gfimax[1])[0]
        gfi_ltone = gfi_ltone + len(np.where(temp > 1.0)[0])
    
    temp = geomfac['w'].copy(); temp = np.absolute(temp); temp[:,0,0,0] = 0
    temp2 = temp.flatten()
    if temp2.shape[0] is not 0:
        gfimax[2] = max(np.max(temp2),gfimax[2]); del(temp2)
        gfimax_indx[2]=np.where(temp >= gfimax[2])[0]
        gfi_ltone = gfi_ltone + len(np.where(temp > 1.0)[0])

    for i in range(3):
        if gfimax[i] == max(gfimax):
            for indx in gfimax_indx[i]:
                if i == 0:
                    if intptype['u'][indx] == 1:
                        max_type1 = max_type1 + 1
                    elif intptype['u'][indx] == 2:
                        max_type2 = max_type2 + 1
                    elif intptype['u'][indx] == 3:
                        max_type3 = max_type3 + 1
                elif i == 1:
                    if intptype['v'][indx] == 1:
                        max_type1 = max_type1 + 1
                    elif intptype['v'][indx] == 2:
                        max_type2 = max_type2 + 1
                    elif intptype['v'][indx] == 3:
                        max_type3 = max_type3 + 1
                else:
                    if intptype['w'][indx] == 1:
                        max_type1 = max_type1 + 1
                    elif intptype['w'][indx] == 2:
                        max_type2 = max_type2 + 1
                    elif intptype['w'][indx] == 3:
                        max_type3 = max_type3 + 1

    for i in intptype['u'].flatten():
            if i == 1:
                type1 = type1 + 1
            elif i == 2:
                type2 = type2 + 1
            elif i == 3:
                type3 = type3 + 1

    for j in intptype['v'].flatten():
            if j == 1:
                type1 = type1 + 1
            elif j == 2:
                type2 = type2 + 1
            elif j == 3:
                type3 = type3 + 1

    for k in intptype['w'].flatten():
            if k == 1:
                type1 = type1 + 1
            elif k == 2:
                type2 = type2 + 1
            elif k == 3:
                type3 = type3 + 1

    print('\nGeometric factor_Max                = %10.3f' %(max(gfimax)))
    print('# of GeomFac_Max(intptype = 1)      = %10d' %(max_type1))
    print('# of GeomFac_Max(intptype = 2)      = %10d' %(max_type2))
    print('# of GeomFac_Max(intptype = 3)      = %10d' %(max_type3))
    print('# of GeomFac larger than 1          = %10d' %(gfi_ltone))
    print('# of intptype equals to 1           = %10d' %(type1))
    print('# of intptype equals to 2           = %10d' %(type2))
    print('# of intptype equals to 3           = %10d' %(type3))

def surf3D(x_coord, y_coord, z_coord, INOUT, filename):
    file = open('../output/ibmpre/surf_body_3D_%s.dat' %(filename), 'w')

    inout_ori = INOUT[1:INOUT.shape[0]-1, 1:INOUT.shape[1]-1, 1:INOUT.shape[2]-1] 
    inout_ip  = INOUT[0:INOUT.shape[0]-2, 1:INOUT.shape[1]-1, 1:INOUT.shape[2]-1]
    inout_im  = INOUT[2:INOUT.shape[0]  , 1:INOUT.shape[1]-1, 1:INOUT.shape[2]-1]
    inout_jp  = INOUT[1:INOUT.shape[0]-1, 0:INOUT.shape[1]-2, 1:INOUT.shape[2]-1]
    inout_jm  = INOUT[1:INOUT.shape[0]-1, 2:INOUT.shape[1]  , 1:INOUT.shape[2]-1]
    inout_kp  = INOUT[1:INOUT.shape[0]-1, 1:INOUT.shape[1]-1, 0:INOUT.shape[2]-2]
    inout_km  = INOUT[1:INOUT.shape[0]-1, 1:INOUT.shape[1]-1, 2:INOUT.shape[2]  ]

    surf = np.absolute(inout_ip - inout_im) + np.absolute(inout_jp - inout_jm) + np.absolute(inout_kp - inout_km)

    for i in range(surf.shape[0]):
        for j in range(surf.shape[1]):
            for k in range(surf.shape[2]):
                if (INOUT[i+1][j+1][k+1] == 0) and (surf[i][j][k]):
                    file.write('%23.15f%23.15f%23.15f\n' %(x_coord[i], y_coord[j], z_coord[k]))

    file.close()

def grid_preprocessing_data(grid):
    file = open('../output/grid/grid.bin', 'w')
    file.write('%d %d %d\n' %(grid['N_x'], grid['N_y'], grid['N_z']))
    file.write('%d %d %d\n' %(grid['Cell_x'], grid['Cell_y'], grid['Cell_z']))
    file.write('%f %f %f\n' %(grid['L_x'], grid['L_y'], grid['L_z']))

    for tmp in [grid['grid_info_x'][i]['coord'] for i in range(1,grid['N_x']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')
    for tmp in [grid['grid_info_y'][j]['coord'] for j in range(1,grid['N_y']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')
    for tmp in [grid['grid_info_z'][k]['coord'] for k in range(1,grid['N_z']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')

    for tmp in [grid['grid_info_x'][i]['idxplus'] for i in range(1,grid['Cell_x']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')
    for tmp in [grid['grid_info_y'][j]['idxplus'] for j in range(1,grid['Cell_y']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')
    for tmp in [grid['grid_info_z'][k]['idxplus'] for k in range(1,grid['Cell_z']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')

    for tmp in [grid['grid_info_x'][i]['idxminus'] for i in range(1,grid['Cell_x']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')
    for tmp in [grid['grid_info_y'][j]['idxminus'] for j in range(1,grid['Cell_y']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')
    for tmp in [grid['grid_info_z'][k]['idxminus'] for k in range(1,grid['Cell_z']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')

    for tmp in [grid['grid_info_x'][i]['lowerfix'] for i in range(1,grid['Cell_x']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')
    for tmp in [grid['grid_info_y'][j]['lowerfix'] for j in range(1,grid['Cell_y']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')
    for tmp in [grid['grid_info_z'][k]['lowerfix'] for k in range(1,grid['Cell_z']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')

    for tmp in [grid['grid_info_x'][i]['upperfix'] for i in range(1,grid['Cell_x']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')
    for tmp in [grid['grid_info_y'][j]['upperfix'] for j in range(1,grid['Cell_y']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')
    for tmp in [grid['grid_info_z'][k]['upperfix'] for k in range(1,grid['Cell_z']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')

    for tmp in [grid['grid_info_x'][i]['c2cd'] for i in range(grid['N_x']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')
    for tmp in [grid['grid_info_y'][j]['c2cd'] for j in range(grid['N_y']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')
    for tmp in [grid['grid_info_z'][k]['c2cd'] for k in range(grid['N_z']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')

    for tmp in [grid['grid_info_x'][i]['1/c2cd'] for i in range(grid['N_x']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')
    for tmp in [grid['grid_info_y'][j]['1/c2cd'] for j in range(grid['N_y']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')
    for tmp in [grid['grid_info_z'][k]['1/c2cd'] for k in range(grid['N_z']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')

    for tmp in [grid['grid_info_x'][i]['f2fd'] for i in range(grid['N_x']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')
    for tmp in [grid['grid_info_y'][j]['f2fd'] for j in range(grid['N_y']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')
    for tmp in [grid['grid_info_z'][k]['f2fd'] for k in range(grid['N_z']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')

    for tmp in [grid['grid_info_x'][i]['1/f2fd'] for i in range(grid['N_x']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')
    for tmp in [grid['grid_info_y'][j]['1/f2fd'] for j in range(grid['N_y']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')
    for tmp in [grid['grid_info_z'][k]['1/f2fd'] for k in range(grid['N_z']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')

    for tmp in [grid['grid_info_x'][i]['center'] for i in range(grid['N_x']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')
    for tmp in [grid['grid_info_y'][j]['center'] for j in range(grid['N_y']+1)]:
        file.write(' %f ' %(tmp))
    file.write('\n')
    for tmp in [grid['grid_info_z'][k]['center'] for k in range(grid['N_z']+1)]:
        file.write(' %f ' %(tmp))

    file.close()

def write_ibmpre_placeholders(grid, lessgs='OFF', htrnfr='OFF', conjg='OFF'):
    print('\n*** WRITING PLACEHOLDER IBM-PREPROCESSING DATA (SAFE DEFAULTS) ***')

    # 1) IBM forcing-point file (mandatory for solver alloinit/ibmpreread)
    with open('../output/ibmpre/ibmpre_fcpts.bin', 'w') as file:
        file.write('0 0 0\n')
        file.write('0 0 0\n')

    # 2) Thermal IBM forcing-point file (mandatory when ihtrans==1)
    if htrnfr == 'ON':
        with open('../output/ibmpre/ibmpre_fcpts_t.bin', 'w') as file:
            file.write('0\n')
            file.write('0\n')

    # 3) LES zero-SGS files (mandatory when iles==1)
    if lessgs == 'ON':
        n1m = grid['Cell_x']
        n2m = grid['Cell_y']
        n3m = grid['Cell_z']

        with open('../output/ibmpre/ibmpre_nutzero.bin', 'w') as file:
            file.write('0\n\n\n\n')

        with open('../output/ibmpre/ibmpre_wallfdvm.bin', 'w') as file:
            # Fluid everywhere: wall-distance mask defaults to 1
            for _k in range(n3m):
                for _j in range(n2m):
                    for _i in range(n1m):
                        file.write(' 1 ')

    # 4) Conjugate heat-transfer fields (mandatory when iconjg==1 and read)
    if conjg == 'ON':
        n1m = grid['Cell_x']
        n2m = grid['Cell_y']
        n3m = grid['Cell_z']

        with open('../output/ibmpre/ibmpre_conjg.bin', 'w') as file:
            # cstar
            for _k in range(n3m):
                for _j in range(n2m):
                    for _i in range(n1m):
                        file.write(' 1.0 ')
            file.write('\n')

            # kstar (6 faces)
            for _l in range(6):
                for _k in range(n3m):
                    for _j in range(n2m):
                        for _i in range(n1m):
                            file.write(' 1.0 ')

def ibm_preprocessing_data(nintp,ninner,fcp,intpindx,geomfac):
    print('\n*** WRITING IMMERSED-BODY DATA INTO ../output/ibmpre/ibmpre_fcpts.bin ***')
    file = open('../output/ibmpre/ibmpre_fcpts.bin', 'w')

    file.write('%d %d %d\n' %(nintp['u'], nintp['v'], nintp['w']))
    file.write('%d %d %d\n' %(ninner['u'], ninner['v'], ninner['w']))

    nu = nintp['u'] + ninner['u']
    nv = nintp['v'] + ninner['v']
    nw = nintp['w'] + ninner['w']

    if nu+nv+nw >= 5000000:
        print('Take a cup of coffee. Plz do not finish the process ...')
    for i in range(nu):
        file.write('%d %d %d\n' %(fcp['u'][i,0], fcp['u'][i,1], fcp['u'][i,2]))
    for j in range(nv):
        file.write('%d %d %d\n' %(fcp['v'][j,0], fcp['v'][j,1], fcp['v'][j,2]))
    for k in range(nw):
        file.write('%d %d %d\n' %(fcp['w'][k,0], fcp['w'][k,1], fcp['w'][k,2]))

    print('--- FORCING PTS (INTP + NON-INTP) DATA RECORD DONE')

    nu = nintp['u']
    nv = nintp['v']
    nw = nintp['w']

    for i in range(nu):
        file.write('%d %d %d\n' %(intpindx['u'][i,0], intpindx['u'][i,1], intpindx['u'][i,2]))
    for j in range(nv):
        file.write('%d %d %d\n' %(intpindx['v'][j,0], intpindx['v'][j,1], intpindx['v'][j,2]))
    for k in range(nw):
        file.write('%d %d %d\n' %(intpindx['w'][k,0], intpindx['w'][k,1], intpindx['w'][k,2]))
    print('--- FCP DIRECTION FOR INTERPOLATION DATA RECORD DONE')

    for i in range(nu):
        for l in range(3):
            for m in range(3):
                for n in range(3):
                    file.write(' %15.8f ' %(geomfac['u'][i,l,m,n]))
        file.write('\n')
    for j in range(nv):
        for l in range(3):
            for m in range(3):
                for n in range(3):
                    file.write(' %15.8f ' %(geomfac['v'][j,l,m,n]))
        file.write('\n')
    for k in range(nw):
        for l in range(3):
            for m in range(3):
                for n in range(3):
                    file.write(' %15.8f ' %(geomfac['w'][k,l,m,n]))
        file.write('\n')
    print('--- GEOMETRIC FACTOR DATA RECORD DONE')

    file.close()

def ibm_preprocessing_data_htransfer(nintp,ninner,fcp,intpindx,geomfac):
    file = open('../output/ibmpre/ibmpre_fcpts_t.bin', 'w')

    file.write('%d\n' %(nintp['t']))
    file.write('%d\n' %(ninner['t']))

    nt = nintp['t'] + ninner['t']

    for t in range(nt):
        file.write('%d %d %d\n' %(fcp['t'][t,0], fcp['t'][t,1], fcp['t'][t,2]))

    nt = nintp['t']

    for t in range(nt):
        file.write('%d %d %d\n' %(intpindx['t'][t,0], intpindx['t'][t,1], intpindx['t'][t,2]))

    for k in range(nt):
        for l in range(3):
            for m in range(3):
                for n in range(3):
                    file.write(' %15.8f ' %(geomfac['t'][k,l,m,n]))
        file.write('\n')
    print('--- HEAT TRANSFER DATA RECORD DONE')

    file.close()


def les_preprocessing_data(nzero,iszero):
    print('\n*** WRITING LES-SGS_ZERO DATA INTO ../output/ibmpre/ibmpre_nutzero.bin ***')
    file = open('../output/ibmpre/ibmpre_nutzero.bin', 'w')

    file.write('%d \n' %(nzero))

    indices = list()
    
    for i in range(iszero.shape[0]):
        for j in range(iszero.shape[1]):
            for k in range(iszero.shape[2]):
                if iszero[i][j][k] == 0:
                    indices.append([i,j,k])

    for index in indices:
        file.write(' %d ' %(index[0]))
    file.write('\n')
    for index in indices:
        file.write(' %d ' %(index[1]))
    file.write('\n')
    for index in indices:
        file.write(' %d ' %(index[2]))

    file.close()

    print('*** WRITING LES-SGS_ZERO DATA INTO ../output/ibmpre/ibmpre_wallfdvm.bin ***')
    file = open('../output/ibmpre/ibmpre_wallfdvm.bin', 'w')

    for k in range(iszero.shape[2]):
        for j in range(iszero.shape[1]):
            for i in range(iszero.shape[0]):
                file.write(' %d ' %(iszero[i][j][k]))

    file.close()

def conjg_preprocessing_data(cstar,kstar):
    print('\n*** WRITING CONJUGATE_HTRANS DATA INTO ../output/ibmpre/ibmpre_conjg.bin ***')
    file = open('../output/ibmpre/ibmpre_conjg.bin', 'w')

    for k in range(cstar.shape[2]):
        for j in range(cstar.shape[1]):
            for i in range(cstar.shape[0]):
                file.write(' %f ' %(cstar[i][j][k]))

    file.write('\n')

    for l in range(kstar.shape[3]):
        for k in range(kstar.shape[2]):
            for j in range(kstar.shape[1]):
                for i in range(kstar.shape[0]):
                    file.write(' %f ' %(kstar[i][j][k][l]))

    file.close()