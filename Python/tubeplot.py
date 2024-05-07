# -*- coding: utf-8 -*-


# =============================================================================
# Transferring tubeplot code to Python
# =============================================================================
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import matplotlib.pyplot as plt
import first_der as fd
import second_der as sd
from mpl_toolkits.mplot3d import art3d

# =============================================================================
# Initial inputs
# =============================================================================

# J = 100
# RadPoints = 20

# minX = 0
# maxX = 9

# Omega3 = np.ones(J)*0.1
# a = np.ones(J)*1.0
# TubeColor = [0.1, 0.5, 0.5] # RGB
# StripeColor = [0.9, 0.1, 0.1] # RGB
# transparency = 0.8
# x2 = np.linspace(minX, maxX, num=J, endpoint=True, retstep=False, dtype=float, axis=0).reshape(J,1)
# y2 = np.linspace(minX, maxX, num=J, endpoint=True, retstep=False, dtype=float, axis=0).reshape(J,1)
# z2 = np.linspace(minX, maxX, num=J, endpoint=True, retstep=False, dtype=float, axis=0).reshape(J,1)

# r = np.concatenate((x2,y2,z2), axis=1, out=None, dtype=None, casting="same_kind")
# =============================================================================
# Comment out all variables above this
# =============================================================================


   # plots a 3D tube about the curve defined by
    # the vector r.  a is the user defined radius
    # of the tube and RadPoints is the number of
    # circumferential patches used to make the tube.

def tubeplot(r,Omega3,a,RadPoints,TubeColor,StripeColor,transparency):
    
    N = len(r)
    # =============================================================================
    # Initialize variables
    # =============================================================================
    e1 = np.zeros([N,3])
    e2 = np.zeros([N,3])
    e3 = np.zeros([N,3])
    De3Ds = np.zeros([N,3])
    DrDs = np.zeros([N,3])
    e32 = np.zeros([N-1,3])
    e1_half = np.zeros([N-1,3])
    e2_half = np.zeros([N-1,3])
    om1 = np.zeros([N-1])
    om2 = np.zeros([N-1])
    om1_half = np.zeros([N-1])
    om2_half = np.zeros([N-1])
    MT = np.zeros([RadPoints,4],dtype=int)
    M = np.zeros([RadPoints,4],dtype=int)
    m1 = np.zeros(N)
    m2 = np.zeros(N)
    curv_half = np.zeros([N-1,3])
    vertices = np.zeros([N*RadPoints,3])
    fsv = np.zeros([RadPoints*N,3])
    
    # =============================================================================
    # START OF CODE - find tangent vector
    # =============================================================================
    
    dr = np.sqrt(np.sum(np.square(r[0:-1,:]-r[1:,:]), axis=1,dtype=float)).reshape(N-1,1)
    
    L1 = fd.first_der(N,dr)       
    L2 = sd.second_der(N,dr)
    
    L1 = L1.toarray()
    L2 = L2.toarray()
    
    DrDs[:,0] = np.dot(L1,r[:,0])
    DrDs[:,1] = np.dot(L1,r[:,1])
    DrDs[:,2] = np.dot(L1,r[:,2])
    
    g = np.sqrt(np.sum(np.square(DrDs),axis=1))
    
    e3[:,0] = DrDs[:,0]/g
    e3[:,1] = DrDs[:,1]/g
    e3[:,2] = DrDs[:,2]/g
    
    # =============================================================================
    # Define e3 and Omega3 at half grids
    # =============================================================================
    rDiff = r[1:N,:]-r[0:N-1,:]
    
    e32[:,0] = rDiff[:,0]/dr[:,0]
    e32[:,1] = rDiff[:,1]/dr[:,0]
    e32[:,2] = rDiff[:,2]/dr[:,0]
    
    om32 = 0.5*(Omega3[0:N-1]+Omega3[1:N])
    g2 = np.sqrt(np.sum(np.square(e32),axis=1))
    
    e32[:,0] = e32[:,0]/g2
    e32[:,1] = e32[:,1]/g2
    e32[:,2] = e32[:,2]/g2
    
    # =============================================================================
    # Define e1(1) and e2(1)
    # =============================================================================
    if e3[0,2]>=1-np.finfo(float).eps:
        e1[0,0] = 1
        e1[0,1] = 0
        e1[0,2] = 0
    else:
        e1[0,0] = -e3[0,0]*e3[0,2]/(np.sqrt(1-np.square(e3[0,2])))
        e1[0,1] = -e3[0,1]*e3[0,2]/(np.sqrt(1-np.square(e3[0,2])))
        e1[0,2] = (np.sqrt(1-np.square(e3[0,2])))  
    
    e2[0,0] = e3[0,1]*e1[0,2] - e3[0,2]*e1[0,1]
    e2[0,1] = e3[0,2]*e1[0,0] - e3[0,0]*e1[0,2]
    e2[0,2] = e3[0,0]*e1[0,1] - e3[0,1]*e1[0,0]
    
    # =============================================================================
    # Find curvatures 
    # =============================================================================
    
    De3Ds[:,0] = np.dot(L2,r[:,0])
    De3Ds[:,1] = np.dot(L2,r[:,1])
    De3Ds[:,2] = np.dot(L2,r[:,2])
    
    # =============================================================================
    # Define half grid curvatures
    # =============================================================================
    
    curv_half[:,0] = 0.5*(De3Ds[0:N-1,0]+De3Ds[1:N,0])
    curv_half[:,1] = 0.5*(De3Ds[0:N-1,1]+De3Ds[1:N,1])
    curv_half[:,2] = 0.5*(De3Ds[0:N-1,2]+De3Ds[1:N,2])
    
    # =============================================================================
    # Integrate e1 and e2
    # =============================================================================
    
    for m in range(0, N-1):
        om1[m] = -(De3Ds[m,0]*e2[m,0] + De3Ds[m,1]*e2[m,1] + De3Ds[m,2]*e2[m,2])
        om2[m] =  (De3Ds[m,0]*e1[m,0] + De3Ds[m,1]*e1[m,1] + De3Ds[m,2]*e1[m,2])
    
    # =============================================================================
    # Integrate e1 and e2 half step
    # =============================================================================
        # e1_half[m,:] = e1[m,:] + np.dot(0.5*dr[m],Omega3[m]*e2[m,:]-om2[m]*e3[m,:])
    
        e1_half[m,:] = e1[m,:] + 0.5*dr[m]*(Omega3[m]*e2[m,:]-om2[m]*e3[m,:])
        e2_half[m,:] = e2[m,:] + 0.5*dr[m]*(om1[m]*e3[m,:]-Omega3[m]*e1[m,:])
    
        m1[m] = np.sqrt(np.sum(np.square(e1_half[m,:])))
        m2[m] = np.sqrt(np.sum(np.square(e2_half[m,:])))                
                    
        e1_half[m,0] = e1_half[m,0]/m1[m]
        e1_half[m,1] = e1_half[m,1]/m1[m]
        e1_half[m,2] = e1_half[m,2]/m1[m]
        
        e2_half[m,0] = e2_half[m,0]/m2[m]
        e2_half[m,1] = e2_half[m,1]/m2[m]
        e2_half[m,2] = e2_half[m,2]/m2[m]
        
    # =============================================================================
    # Define half grid omegas
    # =============================================================================
        om1_half[m] = -(curv_half[m,0]*e2_half[m,0]) - (curv_half[m,1]*e2_half[m,1]) - (curv_half[m,2]*e2_half[m,2])
        om2_half[m] = -(curv_half[m,0]*e1_half[m,0]) - (curv_half[m,1]*e1_half[m,1]) - (curv_half[m,2]*e1_half[m,2])
            
    # =============================================================================
    # Integrate full step
    # =============================================================================
        e1[m+1,:] = e1[m,:] + dr[m]*(om32[m]*e2_half[m,:] - om2_half[m]*e32[m,:])
        e2[m+1,:] = e2[m,:] + dr[m]*(om1_half[m]*e32[m,:] - om32[m]*e1_half[m,:])
    
        m1[m+1] = np.sqrt(np.sum(np.square(e1[m+1])))
        m2[m+1] = np.sqrt(np.sum(np.square(e2[m+1])))
    
        e1[m+1,0] = e1[m+1,0]/m1[m+1]
        e1[m+1,1] = e1[m+1,1]/m1[m+1]
        e1[m+1,2] = e1[m+1,2]/m1[m+1]
    
        e2[m+1,0] = e2[m+1,0]/m2[m+1]
        e2[m+1,1] = e2[m+1,1]/m2[m+1]
        e2[m+1,2] = e2[m+1,2]/m2[m+1]
        
        
    # =============================================================================
    # Define angle about the curve
    # =============================================================================
    theta = np.linspace(0,2*np.pi,RadPoints)
    
    # =============================================================================
    # Define patches
    # =============================================================================
    
    for m in range(0, N):
        for i in range(0,RadPoints):
            for j in range(0,3):
                vertices[(RadPoints*m)+i,j] = r[m,j] + a[m]*(np.cos(theta[i])*e1[m,j] + np.sin(theta[i])*e2[m,j])
    
    # =============================================================================
    # Define connectivity matrix
    # =============================================================================
    
    vec1 = np.linspace(0,RadPoints-1,RadPoints,dtype=int)
    vec2 = np.roll(vec1, -1, axis=0)
    
    for m in range(0,RadPoints):
        MT[m,0] = vec1[m]
        MT[m,1] = vec2[m]
        MT[m,3] = vec1[m] + RadPoints
        MT[m,2] = MT[m,1] + RadPoints
        
    for i in range(2,N):
        for m in range(0,RadPoints):
            M[m,0] = (i-1)*RadPoints + vec1[m];
            M[m,1] = (i-1)*RadPoints + vec2[m];
            M[m,3] = M[m,0] + RadPoints;
            M[m,2] = M[m,1] + RadPoints;
    
        MT = np.append(MT, M, axis = 0)
        
    for k in range(1,N+1):
        for l in range(0,RadPoints-2):
            fsv[RadPoints*(k-1)+l,:] = TubeColor
        fsv[(RadPoints)*k-1,:] = StripeColor
        fsv[(RadPoints)*k-2,:] = StripeColor
    
    
    # =============================================================================
    # Plot it
    # =============================================================================
       
    fig = plt.figure(figsize=(20.0,20.0))
    ax = fig.add_subplot(projection="3d")
    pc = art3d.Poly3DCollection(vertices[MT], facecolors=fsv, edgecolor="black", alpha=transparency)
    ax.add_collection(pc)

    return ax , pc

