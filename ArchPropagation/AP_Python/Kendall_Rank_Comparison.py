# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 18:17:39 2023

For this test, there are 4 

Needs:
    Modular
    PerturbMatrix
    ModularPerturbMatrix
    
@author: colej
"""

import numpy as np
import pandas as pd
from arch_prop_functions import modular, HMPerturbMatrix, PerturbMatrix, RotationPairs
import matplotlib.pyplot as plt
import time
import os
import scipy
import openpyxl

beginning = time.time()
# reset random number generator for reproducability (same as MATLAB)
np.random.seed(1337)
# p = np.array([0.05, 0.025, 0.01, 0.005, 0.025, 0.001])
#perturbation probability for test
p = np.array([0.05, 0.025, 0.01, 0.005, 0.0025, 0.001])
# Number of eigenvectors to retain
k = 100
# Number of computation loops
g = 100


# Number of components in the system
nodes =1024
modules = 8

#intra modular and intermodular
p1d = 0.5
p1s = 0.1
p2d = 0.5
p2s = 0.1


def TurnSymmetric(A, n):
    A = np.triu(A.to_numpy())
    A = A + A.T
    A[(np.arange(0, n), np.arange(0, n))] = 0
    return A
# %%  reference systems
#FD: Fully Dense
#FS: Fully Sparse
#HD: Hierarchical Dense
#HS: Hierarchical Sparse


#fully dense and fully sparse
FD = modular(nodes, 1, p1d, nodes, 1, 1)
FS = modular(nodes, 1, p1s, nodes, 1,1)
FD = TurnSymmetric(FD, nodes)
FS = TurnSymmetric(FS, nodes)




# heirachical dense and sparse, no cascading dependencies
HD = modular(nodes, 8,  p1d, int(nodes/modules), p2d, 0.5)
HS = modular(nodes, 8,  p1s, int(nodes/modules), p2s, 0.5)
HD = TurnSymmetric(HD, nodes)
HS = TurnSymmetric(HS, nodes)

# show

figs, axs = plt.subplots(2, 2)
axs[0, 0].imshow(FD)
axs[0, 1].imshow(FS)
axs[1, 0].imshow(HD)
axs[1, 1].imshow(HS)


# %% Shared Functions


def CalcAndStore(name, mat, F_Eig_df, Eig_df, IRk_df, k, g):
    # inputs data from each calculation and stores it in an excel file
    # Calculate community, bulk, and mean angle in both eigenspaces
    
    Eig_np = Eig_df.to_numpy()
    Eig_gm = np.mean(Eig_np)
    
    IRk_np = IRk_df.to_numpy()
    IRk_gm = np.mean(IRk_np, axis = 0)
    IRk_gs = np.std(IRk_np.astype(np.float64), axis = 0)
    

    mean_np = np.vstack([IRk_gm, IRk_gs, Eig_gm])

    mean_df = pd.DataFrame(mean_np, index = ['IRk Mean','IRk Std','Eig mean'], columns = range(k))
    
    nce = RotationPairs(Eig_gm)
    nbe = k - nce
    mce = np.mean(Eig_gm[0:nce])
    mbe = np.mean(Eig_gm[nce + 1 : k])
    
    
    
    Eig_Calc_Sheet = pd.DataFrame([nce, nbe, mce, mbe], index = \
                                  ['N Eig V in com eig space', \
                                   'N Eig V in bulk eig space', \
                                       'Mean Angle com eig space', \
                                           'Mean angle bulk eig space'], \
                                      columns = ['Data'])

    writer = pd.ExcelWriter(name+'.xlsx', engine = 'xlsxwriter')    
    Eig_df.to_excel(writer, sheet_name='Eig Data')
    F_Eig_df.to_excel(writer, sheet_name='First Eig')
    IRk_df.to_excel(writer, sheet_name='IRk Data')
    mean_df.to_excel(writer, sheet_name = 'Means')
    Eig_Calc_Sheet.to_excel(writer, sheet_name = 'Calcs')
    writer.close()    
    
    #plot and save
    figs, axs = plt.subplots(1, 2, figsize=(10, 5) )
    axs[0].imshow(mat)
    axs[0].set_xticks([])
    axs[0].set_yticks([])
    axs[0].set_title("Unperturbed System")
    axs[1].scatter(np.arange(0,k), Eig_gm)
    
    
    plt.savefig(name + '.png')
    plt.clf()
    return


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

    
    return S, S[0], resouty

# %% Function for Hierarchical and Random Perturbation Tests


def HMPerturbTests(name, mat, p, nodes, modules, k, g):
    """
    

    Parameters
    ----------
    name : File name prefix for image and excel
    mat : Original, unperturbed matrix
    p : Array of perturbation probabilities
    nodes : Number of nodes 
    modules : Number of modules
    k : Number of eigenvectors to retain
    g : Number of iterations for the sample size

    Returns
    -------
    None.

    """
    
    begin = time.time()
    
    EI_df = pd.DataFrame(index=range(len(p)), columns = range(g))
    F_EI_df = EI_df.copy()
    IR_df = EI_df.copy()
    

    for i in range(0, len(p)):
        print("working on perturb prob", p[i])
        for j in range(0,g):
            Er = HMPerturbMatrix(mat, nodes, modules, 0, p[i], 0.5, 0.5)
            EI_df.iloc[i,j], F_EI_df.iloc[i,j], IR_df.iloc[i,j] =  CompareEigKendall(mat, Er, k)
    

    CalcAndStore(name, mat, F_EI_df, EI_df, IR_df, k, g)
    
    print("Calc for ", name,'done in', time.time()-begin,'s')
    return

def RandPerturbTests(name, mat, p, nodes, modules, k, g):
    """
    

    Parameters
    ----------
    name : File name prefix for image and excel
    mat : Original, unperturbed matrix
    p : Array of perturbation probabilities
    nodes : Number of nodes 
    modules : Number of modules
    k : Number of eigenvectors to retain
    g : Number of iterations for the sample size

    Returns
    -------
    None.

    """
    
    begin = time.time()
    
    EI_df = pd.DataFrame(index=range(len(p)), columns = range(g))
    F_EI_df = EI_df.copy()
    IR_df = EI_df.copy()
    

    for i in range(0, len(p)):
        print("working on perturb prob", p[i])
        for j in range(0,g):
            Er = PerturbMatrix(mat, p[i], 0.5, 0.5)
            EI_df.iloc[i,j], F_EI_df.iloc[i,j], IR_df.iloc[i,j] =  CompareEigKendall(mat, Er, k)
    

    CalcAndStore(name, mat, F_EI_df, EI_df, IR_df, k, g)
    
    print("Calc for ", name,'done in', time.time()-begin,'s')
    return
# %% Calculations
#HMPerturbTests('FD_HM', FD, p, nodes, modules, k, g)
#RandPerturbTests('FD_R', FD, p, nodes, modules, k, g)

#HMPerturbTests('FS_HM', FS, p, nodes, modules, k, g)
#RandPerturbTests('FS_R', FS, p, nodes, modules, k, g)

#HMPerturbTests('HD_HM', HD, p, nodes, modules, k, g)
#RandPerturbTests('HD_R', HD, p, nodes, modules, k, g)

#HMPerturbTests('HS_HM', HS, p, nodes, modules, k, g)
#RandPerturbTests('HS_R', HS, p, nodes, modules, k, g)
