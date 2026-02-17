def midpoints(coord, N, linename, output_path):
    """Write the midpoint information for given grid coordinates."""
    file = open(output_path + '/midpts_%s_plot.csv' %(linename), 'w')
    for i in range(1,N):
        file.write('%23.15f' %((coord[i]+coord[i+1])/2))
    file.close()

def deltaplot(coord, N, linename, output_path):
    """Write the spacing (delta) information for given grid coordinates."""
    file = open(output_path + '/delta_%s_plot.csv' %(linename), 'w')
    file.write('%-5s,%-13s,%-13s\n' %('i','x','dx'))
    for i in range(1,N):
        file.write('%5d,%13.5f,%13.5f' %(i, coord[i], coord[i+1]-coord[i]))
        if not i==N:
            file.write('\n')
    file.close()
    
def plane_grid(i_coord, j_coord, N_i, N_j, planename, output_path):
    """Generate the plane grid file -- compatible with tec360"""
    file = open(output_path + '/grid_plane_%s.tec' %(planename), 'w')
    file.write('VARIABLES="%s","%s"\n' %(planename[0], planename[1]))
    file.write('ZONE I=%d,J=%d,F=POINT\n' %(N_i, N_j))
    for j in range(1,N_j+1):
        for i in range(1,N_i+1):
            file.write('%13.5f %12.5f' %(i_coord[i], j_coord[j]))
            if not ((i==N_i+1) and (j==N_j+1)):
                file.write('\n')
    file.close()
