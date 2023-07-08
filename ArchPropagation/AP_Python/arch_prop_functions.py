# -*- coding: utf-8 -*-
"""

@author: colej
arch_prop_functions.py developed from previous MATLAB 

Functions called in other Archictecture Propagation Files

-:func:`connectivity`:  Used to generate a connectivity matrix for a 
                        network of brain components.
-:func:`modular`:       Generates a modular network of size n.
-:func:`HMPerturbMatrix`:Takes HM matrix and produces a matrix starting from
                        modular matrix w/ intra-module connection probability p 
-:func:`ComputeSine`:   Takes in A and E & produces a vector S which is the 
                        angles by which eigenvectors of A rotate under 
                        perturbation by E
-:func:`RotationPairs`: takes in an eigenvalue rotation vector A & returns 
                        index value corresponding to the last eigenvalue in the
                        community eigenspace    
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy import linalg as LA

def ComputeSine(A,E,k):
    """
    % This function takes in A and E, and produces a vector S which is the
    % angles by which eigenvectors of A rotate under perturbation by E
    % with k large eigenvectors to retain
    % Keep track of versios here: 
    % Date: Version 1: 1 Sept2015
    % Author: Somwrita Sarkar, translated to Python by Cole Jetton

    """
    #computer eigenvectors for A
    wA,vA = LA.eig(A)
    wAE,vAE = LA.eig(A+E)
    #note, eigs() in matlab allows to specify top x amount of eigenvectors you 
    #get based on magnitude, but LA.eig() already orders it by magntitude so
    #this works just as well
    
    #now compute sine between eigenvectors of A and eigenvectors of A+E
    S = np.ones(k) #initialize
    for i in range(0,k):
        c_i = np.dot(vA[:,i],vAE[:,i]/(LA.norm(vA[:,i])*LA.norm(vAE[:,i])))
        s_i = np.sqrt(1-c_i**2)
        S[i] = np.arcsin(s_i)*180/np.pi
        
    
    return S

def ComputeSineReverse(A,E,k):
    """
    % This function takes in A and E, and produces a vector S which is the
    % angles by which eigenvectors of A rotate under perturbation by E
    % with k smallest eigenvectors to retain
    % Keep track of versios here: 
    % Date: Version 1: 1 Sept2015
    % Author: Somwrita Sarkar, translated to Python by Cole Jetton

    """
    #computer eigenvectors for A
    wA,vA = LA.eig(A)
    wAE,vAE = LA.eig(A+E)
    #note, eigs() in matlab allows to specify top x amount of eigenvectors you 
    #get based on magnitude, but LA.eig() already orders it by magntitude so
    #this works just as well
    
    #to get the smallest, simply reverse the order
    vA=np.fliplr(vA); vAE=np.fliplr(vAE)
    
    #now compute sine between eigenvectors of A and eigenvectors of A+E
    S = np.ones(k) #initialize
    for i in range(0,k):
        c_i = np.dot(vA[:,i],vAE[:,i]/(LA.norm(vA[:,i])*LA.norm(vAE[:,i])))
        s_i = np.sqrt(1-c_i**2)
        S[i] = np.arcsin(s_i)*180/np.pi
        
    
    return S


def connectivity(n,f_type,params):
    #NOTE! may change since not all 
    """

 This function is used to generate a connectivity matrix
 for a network of brain components.

  Author: Richard Gray, School of Physics, The University of Sydney,
		    N.S.W 2006, Australia
  Using:  Matlab 7.0
  
  Adapted to python by Cole Jetton, Oregonstate University

Returns a connectivity matrix of size n.
The input type is the type of connectivity used.
param are parameters dependent on type

	type, params split based on if_else statement
	----

 	1a	Random
		param1 = Probability of connection

	1b  Random symmetric (Network is symmetric)
		param1 = Probability of connection

	1c	Random partially symmetric
		param1 = Probability of connection
		param2 = Probability of connection being symmetric

	2a	Small World (Newman/Watts/Monasson  no rewiring)
		param1 = degree of underlying regular network
		param2 = probability of shortcut

	2b	Small World (Watts/Strogatz  rewiring)
		param1 = degree of underlying regular network
		param2 = probability of rewiring

      2c  Small World (ours rewiring)
		param1 = degree of underlying regular network
		param2 = probability of rewiring

	3	Scale Free, future work

	4a	Sigmoid distribution
       param1 = minimum probability
       param2 = steepness of transition
       param3 = position where transition occurs

	5a  Randomly distributed with fixed number of connections
       No self conections. Used by Sporns et al
       param1= Number of connections

   5b  Same as 5a but self connections allowed.
       Used by modular.m
		param1= Number of connections

   5c  Randomly distributed with fixed number of connections but
       each node has the same indegree. No self conections. Used by Sporns et al
       param1= Number of connections

   6a  Modular network. Using modular2.
       param1 = number of modules
       param2 = probability of connection in module
       param3 = size of modules
       param4 = number of intermodule connections
       param5 = probability connection exponent
                pe is used to scale pc

   6b  Partially Symmetric Modular network. Using modular2symm.
       param1 = number of modules
       param2 = probability of connection in module
       param3 = size of modules
       param4 = number of intermodule connections
       param5 = probability connection exponent
                pe is used to scale pc
       param6 = probability of symmetric connection
    """
    
    #initialize

    
    if f_type == '1a':
        temp=np.random.rand(n,n)<params
        cij = pd.DataFrame(temp.astype(int))
    else:
        print('please type correct option or specify correct length')
    
    return cij

def HMPerturbMatrix(A,n,m,p,q,pw,pu):
    """ 
    % This function takes in HM matrix A and produces a matrix S starting from
    % a modular matrix with intra-module connection probability p and inter-module
    % connection probability q then rewired such that
    % S(A==0)=1 with probability pw (probability of wiring)
    % S(A==1)=-1 with probability pd (probability of unwiring)
    % A is assumed square and the parameters n,m,r should be correct for A
    
    % This is used to perturb a HM network within or outside a module.
    % To perturb intra-module set q = 0
    % To perturb inter-module set p = 0
    % 
    % Keep track of versions here: 
    % Date: Version 1: 24 November 2015
    """   
    
   

    S = modular(n,m,p,int(n/m),q,0.5)
    
    #where there is an edge in A pertubation may only be 0 or -1
    vu = np.where(A+S.to_numpy()==2)
    Dis = np.random.rand(1,len(vu[0])) #initialize pairs to connect
    for i in range(0,len(Dis)):

        if Dis[0,i]>=pu: #do a probability test to create perturbation or not
            S.iloc[vu[0][i],vu[1][i]]=0 #make no perturbation
        else:
            S.iloc[vu[0][i],vu[1][i]]=-1 #make perturbation
    
    #where there is no edge in A perturbation may only be 0 or 1        
    vw = np.where(A-S.to_numpy()==-1)
    Con = np.random.rand(1,len(vw[0]))
    for i in range(0,len(Dis)):
        if Con[0,i]>=pw:  #do a probability test to create perturbation or not
            S.iloc[vw[0][i],vw[1][i]]=0 #make no perturbation
        else:
            S.iloc[vw[0][i],vw[1][i]]=-1 #make perturbation
    
    #format based on symmetry and set diagonal to zero
    S = np.triu(S.to_numpy())
    S = S+S.T
    S[(np.arange(0,n),np.arange(0,n))]=0
    
    return  S


def idealModular(n,m):
    """
    % Function idealModular generates perfect modular networks with fully
    % connected equally sized modules. 
    % n: number of nodes in the network
    % m: number of modules
    % size of module: n/m needs to be an integer
    """ 
    
    
    return G


def modular(n,nm,pm,km,pc,pe):
    """
   This function generates a modular network of size n.

  Inputs: n size of network
          nm number of modules
          pm probability of connection in module
          km size of modules
          pc probability of intermodule connections
          pe probability connection exponent
          pe is used to scale pc
   Outputs: cij modular network

	Author: Richard Gray, School of Physics, The University of Sydney,
		    N.S.W 2006, Australia
   Using:  Matlab 7.0

   Modifications: (Somwrita Sarkar, 2011): Instead of pe being used to
   scale pc, pc is intermodule connection probablity, and p sets the rate
   of decay for pc. Further, instead of pc*pm being the intermodule, use
   pc directly - for cleaner analytics (values of eigenvalues). 

           5/03/2007
           
   Translation to python: Cole Jetton, Oregon State University, 2023         
   
   Easily generalized to case where inter module connections are created
   with probability that varies with level.
    """

    #memory initialization
    cij = pd.DataFrame(np.zeros((n,n)))
    levels = np.log2(nm); 
    
    #do some checks
    if n != nm*km:
        print('Parameters n,m, and km incompatible')
        return
    if (levels-round(levels)) != 0.0:
        print('Number of Modules not a power of 2')
        return
   
    """
    # First generate the randomly connected modules. This is done
    # by using connectivity to make a random network of the correct 
    # size and then moving each module to the correct position in cij.    
    """
    
    for a in range(1,nm+1):

        module =connectivity(km,'1a',pm) #create the module
        q = np.arange(km*(a-1),km*a) 
        cij.iloc[q,q] = module #put module in correct position

    """
    # Now add the connections between modules. Again use connectivity
    # to generate a random matrix of the correct size and probability of
    # connection. Then move to the correct position in cij. As loop over the 
    # level number the size of the intermodule matrix increases while
    # the number of the connections remains constant. 
    # One intermodule connection is added to ensure that the network
    # is strongly connected (if pc is too small there is a chance no
    # connections are added.)
    """
    impc = pc
    
    for a in range(1,int(levels)+1):
        #below the diagonals
        for b in range(1,int(nm/(2**a))+1):
            imsize = km*(2**(a-1))

            intermod=connectivity(imsize,'1a',impc)#create intermodule connection matrix
            #move into position 

            
            qr = np.arange(2**(a-1)*km + (b-1)*2**a*km, 2**a*km + (b-1)*2**a*km) #2^(a-1)*km+1+(b-1)*2^a*km:2^a*km+(b-1)*2^a*km
            qc = np.arange( (b-1)*2**a*km, 2**(a-1)*km + (b-1)*2**a*km)  #(b-1)*2^a*km:2^(a-1)*km   +(b-1)*2^a*km
            cij.iloc[qr,qc] = intermod


        #above the diagonal
        for b in range(1,int(nm/(2**a))+1):
            imsize = km*(2**(a-1))
            #create intermodule connection matrix
            intermod=connectivity(imsize,'1a',impc)
            #move into position
            qc = np.arange(2**(a-1)*km + (b-1)*2**a*km, 2**a*km + (b-1)*2**a*km) #2^(a-1)*km+1+(b-1)*2^a*km:2^a*km+(b-1)*2^a*km
            qr = np.arange((b-1)*2**a*km, 2**(a-1)*km + (b-1)*2**a*km) #(b-1)*2^a*km:2^(a-1)*km   +(b-1)*2^a*km
            cij.iloc[qr,qc] = intermod
                

        impc=impc*pe #change probability
    
    return cij


def ModularPerturbMatrix(A,m,p,q,pw,pu):
    """ 
    % This function takes in A and produces a matrix S starting from
% a modular matrix with intra-module connection probability p and inter-module
% connection probability q then rewired such that
    % S(A==0)=1 with probability pw (probability of wiring)
    % S(A==1)=-1 with probability pd (probability of unwiring)
    % A is assumed square and the parameters n,m,r should be correct for A
    
    % This is used to perturb a HM network within or outside a module.
    % To perturb intra-module set q = 0
    % To perturb inter-module set p = 0
    % 
    % Keep track of versions here: 
    % Date: Version 1: 24 November 2015
    % Original MATLAB Author: Andy Dong
    % Converted to Python: Cole Jetton, 2023
    """   
    
    n=A.shape[0] #get dimension of A for modular creation

    S = modular(n,m,p,int(n/m),q,1)
    
    #where there is an edge in A pertubation may only be 0 or -1
    vu = np.where(A+S.to_numpy()==2)
    Dis = np.random.rand(1,len(vu[0])) #initialize pairs to connect
    for i in range(0,len(Dis)):

        if Dis[0,i]>=pu: #do a probability test to create perturbation or not
            S.iloc[vu[0][i],vu[1][i]]=0 #make no perturbation
        else:
            S.iloc[vu[0][i],vu[1][i]]=-1 #make perturbation
    
    #where there is no edge in A perturbation may only be 0 or 1        
    vw = np.where(A-S.to_numpy()==-1)
    Con = np.random.rand(1,len(vw[0]))
    for i in range(0,len(Dis)):
        if Con[0,i]>=pw:  #do a probability test to create perturbation or not
            S.iloc[vw[0][i],vw[1][i]]=0 #make no perturbation
        else:
            S.iloc[vw[0][i],vw[1][i]]=-1 #make perturbation
    
    #format based on symmetry and set diagonal to zero
    S = np.triu(S.to_numpy())
    S = S+S.T
    S[(np.arange(0,n),np.arange(0,n))]=0
    
    return  S

def PerturbMatrix(A,p,pw,pu):
    """ 
    % This function takes in HM matrix A and produces a matrix S starting from
    % a modular matrix with intra-module connection probability p and inter-module
    % connection probability q then rewired such that
    % S(A==0)=1 with probability pw (probability of wiring)
    % S(A==1)=-1 with probability pd (probability of unwiring)
    % A is assumed square and the parameters n,m,r should be correct for A
    
    % This is used to perturb a HM network within or outside a module.
    % To perturb intra-module set q = 0
    % To perturb inter-module set p = 0
    % 
    % Keep track of versions here: 
    % Date: Version 1: 24 November 2015
    """   
    
   
    n= A.shape[0] #get shape of A for calculations
    
    S = RandNetwork(n,p)
    
    #where there is an edge in A pertubation may only be 0 or -1
    vu = np.where(A+S.to_numpy()==2)
    Dis = np.random.rand(1,len(vu[0])) #initialize pairs to connect
    for i in range(0,len(Dis)):

        if Dis[0,i]>=pu: #do a probability test to create perturbation or not
            S.iloc[vu[0][i],vu[1][i]]=0 #make no perturbation
        else:
            S.iloc[vu[0][i],vu[1][i]]=-1 #make perturbation
    
    #where there is no edge in A perturbation may only be 0 or 1        
    vw = np.where(A-S.to_numpy()==-1)
    Con = np.random.rand(1,len(vw[0]))
    for i in range(0,len(Dis)):
        if Con[0,i]>=pw:  #do a probability test to create perturbation or not
            S.iloc[vw[0][i],vw[1][i]]=0 #make no perturbation
        else:
            S.iloc[vw[0][i],vw[1][i]]=-1 #make perturbation
    
    #format based on symmetry and set diagonal to zero
    S = np.triu(S.to_numpy())
    S = S+S.T
    S[(np.arange(0,n),np.arange(0,n))]=0
    
    return  S

def RandNetwork(n,p):
    #makes a random network similar to connectivity but symmetric
    A = np.random.rand(n,n)<p
    A = np.triu(A)
    A = A+A.T
    A[(np.arange(0,n),np.arange(0,n))]=0
    return A

def RotationPairs(A):
    """
    This function takes in an eigenvalue rotation vector A
     and returns the index value corresponding to the last eigenvalue
     in the community eigenspace
     See 	arXiv:1510.07064 [physics.soc-ph]
     
    Keep track of versions here: 
    Date: Version 1: 24 November 2015
    Author: Andy Dong\
    Adapted to Python by Cole Jetton, July 2023

    """    
    S=np.array([]) #generate empty S
    peak = False #has peak been reached?
    trough = False #has a trough been reached?
    noise =180/256 #allow for some noise in eigenvalues
    
    if A[0] > 180/16: #unsual situation for ideal modular network
        peak = True
    
    for i in range(0,len(A)-1):
        if not(peak):
            if A[i+1]<A[i] and A[i+1]+noise < A[i]:
                peak = True
                
        else:
            if A[i+1]>A[i] and (A[i+1]-noise > A[i] and A[i+1] - noise > A[0]):
                trough = True


        if trough:
            S=i
            break                  

    if type(S)!=int: #if the eigenvalues increase monotonically into the bulk       
        diff = np.zeros((len(A)-1))
        for i in range(0,len(A)-1):
            diff[i] = A[i+1]-A[i]
        S = np.argmax(diff) + 1   #find the largest jump to the bulk eigenspace 
    
    
    return S