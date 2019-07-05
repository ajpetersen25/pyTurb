### Alec Petersen
# pstats.py
#
# Functions for calculating inertial particle velocity statistics from PTV data
#
###

import numpy as np
from pyTurb.flow import flowstats as fs

def mean_ptv(ptv):
  """
  A function to calculate the ensemble average of ptv velocities from a list
  of numpy arrays
  Inputs:
        ptv:      list of numpy arrays        Each numpy array is of size (N,4) where N is the
                                              number of particles in that PTV snapshot and each 
                                              column represents x position, y position, x velocity
                                              and y velocity respectively
  Outputs:
        avg_vels:  (2,) numpy array            first entry is the ensembled averaged x velocity over
                                              all particles while the second is that of the y velocity
                                              
  """
  xlist = []
  ylist = []
  for i in range(0,len(ptv)):
    xvels = ptv[:,2]
    yvels = ptv[:,3]
    xlist.extend(xvels)
    ylist.extend(yvels)
  avg_vels np.array([np.mean(xvels), np.mean(yvels)])
  return avg_vels

def rms_ptv(ptv):
    """
    A function to calculate the Root Mean Square of ptv velocities from a list
    of numpy arrays
    Inputs:
          ptv:      list of numpy arrays        Each numpy array is of size (N,4) where N is the
                                                number of particles in that PTV snapshot and each 
                                                column represents x position, y position, x velocity
                                                and y velocity respectively
    Outputs:
          rms_vels:  (2,) numpy array           first entry is the Root Mean Sqaure x velocity over
                                                all particles while the second is that of the y velocity

    """
    xlist = []
    ylist = []
    for i in range(0,len(ptv)):
      xvels = ptv[:,2]
      yvels = ptv[:,3]
      xlist.extend(xvels)
      ylist.extend(yvels)
    rms_vels np.array([fs.rms(np.asarray(xvels)), fs.rms(np.asarray(yvels))])
    return rms_vels
