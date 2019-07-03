### Alec Petersen
#
# interpolation.py
#
# q_at_p -- function for interpolation a fluid quantity (velocity, vorticity, etc) at a sub-grid
# particle location
#
# grid_ptv -- function for binning particle velocities (using the mean) onto a regular-spaced grid
###
# imports
import numpy as np
from scipy import interpolate
import scipy.stats as stats

def q_at_p(points, q, dstep):
    """
    Inputs:

    points --- (n,2) array of particle positions in pixel units
    q      --- 2d array of fluid field values
    dstep  --- int, pixel spacing between piv vectors

    Returns:
    q_at_ps --- n,1 array of interpolated values at points
    """
    X,Y = np.meshgrid(dstep * np.arange(0,q.shape[0]) + dstep,dstep * np.arange(0,q.shape[1]) + dstep)
    xi=X.ravel()/dstep
    yi=Y.ravel()/dstep
    z = q.ravel()
    xi=xi[~np.isnan(z)]
    yi=yi[~np.isnan(z)]
    xi=xi.reshape(len(xi),1)
    yi=yi.reshape(len(yi),1)
    z = z[~np.isnan(z)]
    q_at_p = interpolate.griddata(np.concatenate((xi,yi),axis=1),z,points/dstep,method='linear')
    return q_at_p



def grid_ptv(vel_data, grid_dim,ix,iy):
    """
    inputs: vel_data --- numpy array of [x,y,velocity]
            grid_dim --- dimension of grid wanted when structuring the data
                         can be an array [row,col] or a scalar, which results in
                         a square grid of size grid_dim x grid_dim
    outputs: vel_grid: a 2d array of gridded horizontal velocities
    """

    if hasattr(grid_dim,"__len__"):
        vel_data[:,0] = np.floor(vel_data[:,0]/(ix/grid_dim[0]))
        vel_data[:,1] = np.floor(vel_data[:,1]/(iy/grid_dim[1]))
        vel_grid = stats.binned_statistic_2d(vel_data[:,1],vel_data[:,0],vel_data[:,2],
                            statistic='count',
                            bins=[np.arange(0,grid_dim[0]+1), np.arange(0,grid_dim[1]+1)])

    else:
        vel_data[:,0] = np.floor(vel_data[:,0]/(ix/grid_dim))
        vel_data[:,1] = np.floor(vel_data[:,1]/(iy/grid_dim))
        vel_grid = stats.binned_statistic_2d(vel_data[:,1],vel_data[:,0],vel_data[:,2],
                                           statistic='mean',bins=np.arange(0,grid_dim+1))

    return vel_grid.statistic
