# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 11:58:42 2023

@author: colej
"""
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from arch_prop_functions import *
import matplotlib.pyplot as plt
import time
import csv
from numpy import linalg as LA
from infinite_regress_functions import *
import scipy


#constants
k = 100 #number of eigenvectors to retain
g = 30 #number of computation loops


# %%

G = nx.read_edgelist("cj610_edgelist.txt", nodetype=int)

CJ610_DSM = nx.to_numpy_array(G)


#DV = list(G.degree())

NodeDegree = [val for (node, val) in G.degree()]

MaxNodeIndex = np.argsort(NodeDegree)[-10:]

MinNodeIndex = np.argsort(NodeDegree)[0:5]

"""
Max:
# (307) Front Frame to Compressor casing nut 17: 1056
# (2005) Compressor rotor stage 7 blade 87: 303
# (1751) Compressor rotor stage 6 blade 87: 295
# (1498) Compressor rotor stage 5 blade 87: 294
# (996) Compressor rotor stage 3 blade 80: 271

Min: Note, there's a bunch that have a degree of 1, so these are arbitrary-ish
# (36) Front frame ring shim 14
# (68) Front frame ring oversize insert 4
# (2699) Combustion lining cowl to dome ID rivet 12
# (66) Front frame ring oversize insert 2
# (65) Front frame ring oversize insert 1    
"""


# %%

def PerturbComponent(A,c,pw,pu):
    """ 
    This function takes in A and c (component being changed), and produces a matrix S starting from
    a random matrix with probability p and then rewired such that
    S(A==0)=1 with probability pw (probability of wiring)
    S(A==1)=-1 with probability pd (probability of unwiring)
    A is assumed square
     
    Keep track of versions here: 
    Date: Version 1: 10 November 2015
    Author: Andy Dong
    Adapted to Python by Cole Jetton, 2023
    """   
    
   
    n= A.shape[0] #get shape of A for calculations
    
    S = pd.DataFrame(np.zeros((n,n)))
    Sc = S.to_numpy()
    #Where there is an edge in A with component c perturbation may only be 0 or -1
    vu = np.where(A[c-1,:]+Sc[c-1,:]==1)
    Dis = np.random.rand(1,len(vu[0])) #initialize pairs to connect

    for i in range(0,Dis.shape[1]):

        if Dis[0,i]>=pu: #do a probability test to create perturbation or not
            S.iloc[c-1,vu[0][i]]=-1 #make perturbation
        else:
            S.iloc[c-1,vu[0][i]]=0 #don't make perturbation
    
    #Where there is no edge in A with component c perturbation may only be 0 or 1
    Sc = S.to_numpy()      
    vw = np.where(A[c-1,:]-Sc[c-1,:]==0)

    Con = np.random.rand(1,len(vw[0]))
    for i in range(0,Con.shape[1]):

        if Con[0,i]>=pw:  #do a probability test to create perturbation or not
            S.iloc[c-1,vw[0][i]]=1 #make perturbation
        else:
            S.iloc[c-1,vw[0][i]]=0 #don't make perturbation
    
    #format based on symmetry and set diagonal to zero
    S = np.triu(S.to_numpy())
    S = S+S.T
    S[(np.arange(0,n),np.arange(0,n))]=0
    
    return  S

# %%
def CompareEigKendall(A, p_A, k):
    # Calculates the Eigenvector rotation value and k statistic for analysis    
    n = A.shape[0]
    e = np.ones((1,n))    
    w1, v1 = np.linalg.eig(A)
    w2, v2 = np.linalg.eig(A + p_A)

    #take first eigen for comparison 
    q1 = v1[:,0]
    q2 = v2[:,0]
    #print(q2)    
    S = np.ones(k) #initialize

    for i in range(0, k):
        # sine calculations
        c_i = np.dot(v1[:, i], v2[:,i]/(np.linalg.norm(v1[:, i])*np.linalg.norm(v2[:, i])))
        s_i = np.sqrt(1-c_i**2)
        S[i] = np.arcsin(s_i)*180/np.pi
                        
    y1 = np.real(np.divide(k*q1,np.sum(q1*e)))
    y2 = np.real(np.divide(k*q2,np.sum(q2*e)))

    #print('k check',y1, y2, S)
    resy = scipy.stats.kendalltau(y1,y2)
    resouty = resy.statistic

    
    return S, resouty

# %%

#MAIN KNOB: Component for analysis on this run
n = 65

a = pd.DataFrame(index=range(1),columns = range(g))
x = a.copy()
y = a.copy()
z = a.copy()
aK = a.copy()
xK = a.copy()
yK = a.copy()
zK = a.copy()

starttime = time.time()
for i in range(0, g):
    Er1 = PerturbComponent(CJ610_DSM ,n,.1,.1)
    Er2 = PerturbComponent(CJ610_DSM ,n,.5,.5)
    Er3 = PerturbComponent(CJ610_DSM ,n,.9,.1)
    Er4 = PerturbComponent(CJ610_DSM ,n,.1,.9)
    a.iloc[0,i], aK.iloc[0,i] =CompareEigKendall(CJ610_DSM ,Er1,k)
    x.iloc[0,i], xK.iloc[0,i] =CompareEigKendall(CJ610_DSM ,Er2,k)
    y.iloc[0,i], yK.iloc[0,i] =CompareEigKendall(CJ610_DSM ,Er3,k)
    z.iloc[0,i], zK.iloc[0,i]  =CompareEigKendall(CJ610_DSM ,Er4,k)
    
print(time.time()-starttime)

# %%

#find mean from tests
agrandmean = np.mean(a[0])
xgrandmean = np.mean(x[0])
ygrandmean = np.mean(y[0])
zgrandmean = np.mean(z[0])

#plotting results
figs, axs = plt.subplots(2,2)
axs[0,0].scatter(np.arange(0,k),agrandmean,marker='.')
axs[0,0].set_title('Perturbed pw=.1 pu=.1')

axs[0,1].scatter(np.arange(0,k),xgrandmean,marker='.')
axs[0,1].set_title('Perturbed pw=0.5, pu=0.5')

axs[1,0].scatter(np.arange(0,k),ygrandmean,marker='.')
axs[1,0].set_title('Perturbed pw=0.9, pu=0.1')

axs[1,1].scatter(np.arange(0,k),zgrandmean,marker='.')
axs[1,1].set_title('Perturbed pw=0.1, pu=0.9')


# %% eigenvalue bulk eigenspace calculations



# Number of eigenvalues in community eigenspace
nace = RotationPairs(agrandmean)
# Number of eigenvalues in bulk eigenspace
nabe = k - nace
# Mean angle community eigenspace
mace = np.mean(agrandmean[0:nace])
# Mean angle bulk eigenspace
mabe = np.mean(agrandmean[nace + 1 : k])

# Number of eigenvalues in community eigenspace
nxce = RotationPairs(xgrandmean)
# Number of eigenvalues in bulk eigenspace
nxbe = k - nxce
# Mean angle community eigenspace
mxce = np.mean(xgrandmean[0:nxce])
# Mean angle bulk eigenspace
mxbe = np.mean(xgrandmean[nxce + 1 : k])

# Number of eigenvalues in community eigenspace
nyce = RotationPairs(ygrandmean)
# Number of eigenvalues in bulk eigenspace
nybe = k - nyce
# Mean angle community eigenspace
myce = np.mean(ygrandmean[0:nyce])
# Mean angle bulk eigenspace
mybe = np.mean(ygrandmean[nyce + 1 : k])

# Number of eigenvalues in community eigenspace
nzce = RotationPairs(zgrandmean)
# Number of eigenvalues in bulk eigenspace
nzbe = k - nzce
# Mean angle community eigenspace
mzce = np.mean(zgrandmean[0:nzce])
# Mean angle bulk eigenspace
mzbe = np.mean(zgrandmean[nzce + 1 : k])


#'%u & %.1f & %.1f & %.1f & %u & %.1f & %.1f & %.1f & %u & %.1f & %.1f & %.1f & %.1f',nace, ...
# mace, mabe, mabe-mace, nxce, mxce, mxbe, mxbe-mxce, nyce, myce, mybe, mybe-myce,
# max([mabe-mace,mxbe-mxce,mybe-myce])-min([mabe-mace,mxbe-mxce,mybe-myce]));
print(
    "%d" % nace,
    "& %0.1f" % mace,
    "& %0.1f" % mabe,
    "& %0.1f" % (mabe - mace),
    "& %d" % nxce,
    "& %0.1f" % mxce,
    "& %0.1f" % mxbe,
    "& %0.1f" % (mxbe - mxce),
    "& %d" % nyce,
    "& %0.1f" % myce,
    "& %0.1f" % mybe,
    "& %0.1f" % (mybe - myce),
    "& %d" % nzce,
    "& %0.1f" % mzce,
    "& %0.1f" % mzbe,
    "& %0.1f" % (mzbe - mzce),
    "& %0.1f"
    % (
        max(mabe - mace, mxbe - mxce, mybe - myce, mzbe-mzce)
        - min(mabe - mace, mxbe - mxce, mybe - myce, mzbe-mzce)
    ),
)

print('Don\'t forever to copy in the kendall stuff!')