# -*- coding: utf-8 -*-

import numpy as np
# import math
import scipy as sc

def second_der(LATT,dx):
    
    ## creates a matrix that calculates the second order
    ## second derivative with unequal spacings dx

    N = LATT

    ##  make the big matrix; make it sparse 
    
    ia=np.zeros([3*N],dtype=int) #i index assigned to zero
    ja=np.zeros([3*N],dtype=int)   #j index assigned to zero
    za=np.zeros([3*N])   #value index assigned to zero
    
    # =============================================================================
    # Diagonal piece
    # =============================================================================
    
    ia[0]= 0
    ja[0]= 0
    za[0]=(2.0/dx[0]/(dx[0]+dx[1]))

    
    ia[1:N-1]= np.linspace(1,N-2,num=N-2)
    ja[1:N-1]= np.linspace(1,N-2,num=N-2)
    za[1:N-1]= -2/dx[1:N-1,0]/dx[0:N-2,0]
    
    ia[N-1]= N-1
    ja[N-1]= N-1
    za[N-1]= (2.0/dx[N-2,0]/(dx[N-2,0]+dx[N-3,0]))
   
    # =============================================================================
    # 1st lower diagonal
    # =============================================================================
      
    ia[N:2*N-2]= np.linspace(1,N-2,num=N-2)	
    ja[N:2*N-2]= np.linspace(0,N-3,num=N-2)	
    za[N:2*N-2]= 2/dx[0:N-2,0]/(dx[0:N-2,0]+dx[1:N-1,0])
    
    ia[2*N-2]= N-1	
    ja[2*N-2]= N-2	
    za[2*N-2]= -(2/dx[N-2]/dx[N-3])
       
    # =============================================================================
    # 2nd lower diagonal
    # =============================================================================
    
    ia[2*N-1]= N-1	
    ja[2*N-1]= N-3	
    za[2*N-1]= 2/(dx[N-2]+dx[N-3])/dx[N-3]

    # =============================================================================
    # 1st upper diagonal
    # =============================================================================
    ia[2*N]= 0
    ja[2*N]= 1
    za[2*N]= -2/dx[0]/dx[1]  
    
    ia[2*N+1:3*N-1]= np.linspace(1, N-2, num = N-2)
    ja[2*N+1:3*N-1]= np.linspace(2, N-1, num = N-2)	
    za[2*N+1:3*N-1]= 2/dx[1:N-1,0]/(dx[0:N-2,0]+dx[1:N-1,0])
    
    # =============================================================================
    # 2nd upper diagonal
    # =============================================================================
    ia[3*N-1]= 0
    ja[3*N-1]= 2
    za[3*N-1]=  2/(dx[0]+dx[1])/dx[1] 

    # za = np.round(za, 12)
    # =============================================================================
    # Make it sparse
    # =============================================================================
    Ld = sc.sparse.csr_matrix((za,(ia,ja)), dtype=np.float)
    Ld.eliminate_zeros()
    return Ld