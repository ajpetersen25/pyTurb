### Alec Petersen
#
# rdf.py -- calculating the radial distance function for inertial particles in turbulence
#
# rdf_pair_2d -- calculates mirrored RDF
#
# calc_rdf -- function for calling rdf_pair_2d and parallelizing frames on multiple cores for speed
#
###

# Imports
import numpy as np
import multiprocessing
import psutil
import itertools

def rdf_pair_2d(params):
    # takes the PTV (pos, vel) info for a single frame, along with the RDF bins, and 
    # vertical & horizontal spatial limits of the field-of-view and calculates the 2D RDF
    # also uses mirroring to extend the spatial extent of RDF
    
    frame,ptv,rvec,xlim,ylim = params
    pframe = ptv[frame]
    pos = pframe[:,1:3] #extract particle positions
    avg_pd = pos.shape[0]/(xlim*ylim) # average particle density
    m = len(rvec)
    r12_array = np.array([],dtype=np.int64).reshape(pos.shape[0],0) #initialize distance array
    rbin = np.zeros((m-1,2))
    rbin[:,0] = rvec[1:] - np.diff(rvec)/2
    
    # mirror particle fields in all directions
    LL = np.array([0,ylim]-pos - [0,ylim])
    L = np.array([-pos[:,0],pos[:,1]]).transpose()
    UL = np.array([L[:,0],ylim*2-L[:,1]]).transpose()
    U = np.array([pos[:,0],ylim*2-pos[:,1]]).transpose()
    D = np.array([pos[:,0],-pos[:,1]]).transpose()
    R = np.array([xlim*2-pos[:,0],pos[:,1]]).transpose()
    LR = np.array([R[:,0],-R[:,1]]).transpose()
    UR = np.array([R[:,0],ylim*2-R[:,1]]).transpose()
    # calculate distances between all particles
    for p in range(0,len(pos)):
        r12 = np.linalg.norm(pos[p,:]-pos,axis=1,keepdims=True)
        r12_LL = np.linalg.norm(pos[p,:]-LL, axis=1,keepdims=True)
        r12_UL = np.linalg.norm(pos[p,:]-UL, axis=1,keepdims=True)
        r12_L = np.linalg.norm(pos[p,:]-L, axis=1,keepdims=True)
        r12_U = np.linalg.norm(pos[p,:]-U, axis=1,keepdims=True)  
        r12_D = np.linalg.norm(pos[p,:]-D, axis=1,keepdims=True)
        r12_R = np.linalg.norm(pos[p,:]-R, axis=1,keepdims=True)
        r12_LR = np.linalg.norm(pos[p,:]-LR, axis=1,keepdims=True)  
        r12_UR = np.linalg.norm(pos[p,:]-UR, axis=1,keepdims=True)  
        r12_array = np.concatenate((r12_array,r12,r12_LL,r12_UL,r12_L,r12_U,
                               r12_D,r12_R,r12_LR,r12_UR),axis=1)
    # bin distances
    for j in range(0,m-1):
        rbin[j,1] = np.sum(np.bitwise_and(r12_array.flatten() > rvec[j],
                r12_array.flatten()<rvec[j+1]))/(np.pi*(rvec[j+1]**2-rvec[j]**2)*avg_pd*len(pos))

    return rbin
    
def calc_rdf(ptv,xlim,ylim,rvec,cores=1):
    """
    Inputs:
           ptv:        List of (N,4) numpy arrays        list of arrays of ptv data (xpos, ypos, xvel, yvel) 
           xlim:       float                             spatial extent of FOV in x-direction
           ylim:       float                             spatial extent of FOV in y-direction
           rvec:       numpy array                       array of rdf bins 
           cores:      int                               number of cores to use for parallel processing
    Outputs:
           results:    (N,2) numpy array                 first col is the rdf bins vector, second contains the rdf values
    
    """

    #if __name__ == 'calc_rdf':
    param1 = range(0,len(ptv)*2,2)
    f_tot =len(range(0,len(ptv.ptv)*2,2))
    param2 = ptv
    param3 = rvec
    param4 = xlim
    param5 = ylim
    objList = zip(param1,itertools.repeat(param2,times=f_tot),
                  itertools.repeat(param3,times=f_tot),
                  itertools.repeat(param4,times=f_tot),
                  itertools.repeat(param5,times=f_tot))
        
    avail_cores = cores
    pool = multiprocessing.Pool(processes=avail_cores)

    instantaneous_results = pool.map(rdf_pair_2d,objList)
    pool.close()
    pool.join()
    results = np.nanmean(instantaneous_results,axis=0)

    return results
