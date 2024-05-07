# -*- coding: utf-8 -*-


# =============================================================================
# Drive the Tubeplot .py file
# =============================================================================

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import scipy as sc
import math
import matplotlib
import copy as cp
import matplotlib.pyplot as plt
import first_der as fd
import second_der as sd
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LightSource
from mpl_toolkits.mplot3d import art3d
import tubeplot as TB

# =============================================================================
# Initialize
# =============================================================================
N = 60
r = np.zeros([N,3])
a = np.ones([N])*1
tTime = 50
Omega3 = np.ones([N])*.2
RadPoints = 10 
TubeColor = [0.1, 0.9, 0.9]
StripeColor = [0.9, 0.9, 0.1]
opacity = 0.9

for i in range(0,tTime):
    theta = np.linspace(0+(2*np.pi*i/tTime),4*np.pi+(2*np.pi*i/tTime),N)
    r[:,0] = 0.2*np.cos(theta)
    r[:,1] = theta
    r[:,2] = np.zeros([N])
    
    [ax, pc] = TB.tubeplot(r,Omega3,a,RadPoints,TubeColor,StripeColor,opacity)
    
    minX = np.min(r[:,0])-np.max(a)
    maxX = np.max(r[:,0])+np.max(a)
    minY = np.min(r[:,1])-np.max(a)
    maxY = np.max(r[:,1])+np.max(a)
    minZ = np.min(r[:,2])-np.max(a)*2
    maxZ = np.max(r[:,2])+np.max(a)*2
    
    diffX = maxX-minX
    diffY = maxY-minY
    diffZ = maxZ-minZ
    
    absDiff = np.max([abs(diffX), abs(diffY), abs(diffZ)])
    
    ax.set_xlim3d(minX, minX+absDiff)
    ax.set_ylim3d(minY, minY+absDiff)
    ax.set_zlim3d(minZ, minZ+absDiff)
   
    # Hide grid lines
    ax.grid(False)
    
    # Hide axes ticks
    
    # ax.set_xticks([])
    # ax.set_yticks([])
    # ax.set_zticks([])
    
    # Set view angle
    ax.view_init(-135, 35)
    
    # Set background color (R,G,B, transparency)
    # ax.xaxis.set_pane_color((.8,.8,.8, 0.0))
    # ax.yaxis.set_pane_color((.8,.8,.8, 0.0))
    # ax.zaxis.set_pane_color((.8,.8,.8, 0.0))

    plt.show()
    

