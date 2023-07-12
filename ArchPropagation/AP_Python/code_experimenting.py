# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 19:41:37 2023

@author: colej
"""
import numpy as np
import pandas as pd
n= 8
import matplotlib.pyplot as plt
from numpy import linalg as LA
from arch_prop_functions import *

p = np.array([0.05, 0.025, 0.01, 0.005, 0.025, 0.001])
# Number of eigenvectors to retain
k=100;
# Number of computation loops
g=100;
a = pd.DataFrame(index=range(len(p)),columns = range(g))

b = np.random.rand(10,10)
b = np.array([[1, 2, 3],[4, 5, 6],[7,8,9]])
w,v=LA.eig(b)


def HMNetwork(n,m,r,p,q):
    """
    This function takes produces a matrix A which is a hierarchically modular
    network of n nodes, m modules per level, and r levels where p is
    the edge probability within a module and q is the connectivity off the
    hierarchy
     
    Keep track of versions here: 
    Date: Version 1: 9 October 2015
    Author: Andy Dong
    Adapted to Python by Cole Jetton, 2023
    """
    # Initial lowest level hierarchy module
    # Number of nodes in lowest level hierarchy
    Lnodes = int(n/(2*2**(r-1)))
    A=modular(Lnodes,m,p,int(Lnodes/m),q,1)
    A=A.to_numpy()
    
    for i in range(0,r):
        P = RandNetwork(int(Lnodes*2**i),p*q**(i+1)) #altered from original based on python vs MATLAB indexing schemes
        #update A with a terrible way to do it
        A = np.vstack((np.hstack((A,P)),np.hstack((P,A))))
    

    A = np.triu(A)
    A = A+A.T
    n = A.shape[0]
    A[(np.arange(0,n),np.arange(0,n))]=0
    return A

test = HMNetwork(16,2,2,0.5,1)