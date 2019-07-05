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
          xrms:  float                          rms velocity in the x-direction
          yrms:  float                          rms velocity in the y-direction                                 

    """
    xlist = []
    ylist = []
    for i in range(0,len(ptv)):
      xvels = ptv[:,2]
      yvels = ptv[:,3]
      xlist.extend(xvels)
      ylist.extend(yvels)
    return fs.rms(np.asarray(xvels)), fs.rms(np.asarray(yvels))
  
  def cal_vol_frac(ptv,dp,FOV_vol):
    """
    Estimates the solid volume fraction of the multiphase flow by calculating the volume of
    particles in each ptv frame, comparing to the field-of-view of the experiment and averaging
    over all frames
    Inputs:
          ptv:      list of numpy arrays        Each numpy array is of size (N,4) where N is the
                                                number of particles in that PTV snapshot and each 
                                                column represents x position, y position, x velocity
                                                and y velocity respectively
          dp:       float                       expected particle diameter
          FOV_vol:  float                       volume of the experimental domain
    Outputs:
          ppf:      (N,) numpy array            particles-per-frame
          Phi_v:    float                       estimated solid-volume fraction
    """
                                            
    ppf = np.zeros(len(ptv))
    for i in range(0,len(ptv)):
      ppf[i] = ptv[i].shape[0]
    avg_ppf = np.mean(ppf)
    Vp_tot = ((dp**3)*3.14159265359/6)*avg_ppf
    Phi_v = Vp_tot/FOV_vol
                                            
                                            
    
