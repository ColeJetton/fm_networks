# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 19:32:58 2023

@author: colej

Pratt and Whitney DSM analysis
Imports PW engine DSM
Performs Infinite Regress Analysis
Performs arch propigation

Graphs results at this point, maybe we can save stuff
"""
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from arch_prop_functions import *
import matplotlib.pyplot as plt
import time

# %%  OLD (for reference) import and turn into dsm:
G = nx.read_edgelist("pw_edgelist.txt", nodetype=int)

mapping = {
    1: "FAN1",
    2: "FAN2",
    3: "FAN3",
    4: "FAN4",
    5: "FAN5",
    6: "FAN6",
    7: "FAN7",
    
    8: "LPC1",
    9: "LPC2",
    10: "LPC3",
    11: "LPC4",
    12: "LPC5",
    13: "LPC6",
    14: "LPC7",
    
    15: "HPC1",
    16: "HPC2",
    17: "HPC3",
    18: "HPC4",
    19: "HPC5",
    20: "HPC6",
    21: "HPC7",
    
    22: "CC1",
    23: "CC2",
    24: "CC3",
    25: "CC4",
    26: "CC5",
    
    27: "HPT1",
    28: "HPT2",
    29: "HPT3",
    30: "HPT4",
    31: "HPT5",
    
    32: "LPT1",
    33: "LPT2",
    34: "LPT3",
    35: "LPT4",
    36: "LPT5",
    37: "LPT6",
    
    38: "MC1",
    39: "MC2",
    40: "MC3",
    41: "MC4",
    42: "MC5",
    43: "MC6",
    44: "MC7",
    
    45: "EC1",
    46: "EC2",
    47: "EC3",
    48: "EC4",
    49: "EC5",
    50: "EC6",
    51: "EC7",
    52: "EC8",
    53: "EC9",
    54: "EC10",
}
G = nx.relabel_nodes(G, mapping)

PW_DSM = nx.to_numpy_array(G)  # node names are eliminated in this form


# %% 
from numpy import genfromtxt
from numpy import linalg as LA
PW_DSM = genfromtxt('PW_DSM.csv',delimiter=',')
PW_DSM[0][0]=0 #issue with first index as nan

# %% 
"""
def Infinite_Regress(A):
    n= A.shape[0]
    
    return
"""

n = PW_DSM.shape[0]
PW_DSM_ones = np.copy(PW_DSM); PW_DSM_ones[PW_DSM_ones== 2] =1
e1 = sum(PW_DSM_ones.T)
e2 = sum(PW_DSM_ones)
e = np.ones((1,54))
#_vals, e_vecs = LA.eig(CJ610_DSM.T)

#find eigenvector and take the first one (corresponds to largest eigenvalue)
v1, p1 = LA.eig(PW_DSM.T);p = p1[:,0]#p = np.real(p[:,0]) 

v2, q1 = LA.eig(PW_DSM);q= q1[:,0]#q = np.real(q[:,0].T)

#x = (n*p)/(p.T*e), y = (n*q)/(q.T*e)

#x = np.divide(n*p,np.multiply(p.T,e1));x=np.real(x)
#y = np.divide(n*q,np.multiply(q.T,e2));y=np.real(y)

x = np.divide(n*p,np.sum(p*e));x=np.real(x)
y = np.divide(n*q,np.sum(q*e));y=np.real(y)

# %%

names = list(mapping.values())

plt.scatter(x, y)
plt.yscale("log")
plt.xscale("log")

for i in range(0,54):
    plt.annotate(names[i],(x[i],y[i]))
