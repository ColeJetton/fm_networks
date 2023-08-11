# -*- coding: utf-8 -*-
"""
Created on Sun Aug  6 13:02:10 2023

@author: colej
"""
import numpy as np
from numpy import linalg as LA
from scipy.stats import wilcoxon,kendalltau

def ComputeSine_and_InfReg(A,E,k):
    """
    % This function takes in A and E, and produces a vector S which is the
    % angles by which eigenvectors of A rotate under perturbation by E
    % with k large eigenvectors to retain
    % it also produces IR_u and IR_p which is the percentage of
    values above 1 for the infinite regress marking "high "
    % Keep track of versios here: 
    % Date: Version 1: 1 Sept2015
    % Author: Somwrita Sarkar, translated to Python by Cole Jetton

    """
    n = A.shape[0]  # number of components
    e = np.ones((1, n))
    
    #computer eigenvectors for A
    wA,vA = LA.eig(A)
    wAE,vAE = LA.eig(A+E)
    #note, eigs() in matlab allows to specify top x amount of eigenvectors you 
    #get based on magnitude, but LA.eig() already orders it by magntitude so
    #this works just as well
    
    #now compute sine between eigenvectors of A and eigenvectors of A+E
    S = np.ones(k) #initialize
    IR_u = np.ones(k)
    IR_p = np.ones(k)
    for i in range(0, k):
        # sine calculations
        c_i = np.dot(vA[:, i], vAE[:,i]/(LA.norm(vA[:, i])*LA.norm(vAE[:, i])))
        s_i = np.sqrt(1-c_i**2)
        S[i] = np.arcsin(s_i)*180/np.pi
        # infinite regress percentage calculations
        x_u = np.real(np.divide(n*vA[:, i], np.sum(vA[:, i]*e)))
        x_p = np.real(np.divide(n*vAE[:, i], np.sum(vAE[:, i]*e)))        
        IR_u[i] = np.count_nonzero(x_u>1)/n
        IR_p[i] = np.count_nonzero(x_p>1)/n

    return S, IR_u, IR_p
#sus_score = np.count_nonzero(x<1)/PW_DSM_s.shape[1]

def AssymetricRandNetwork(n,p):
    
    
    A = np.random.rand(n,n)<p

    A[(np.arange(0,n),np.arange(0,n))]=0
    return A.astype(int)


def Compare(A,p_A,k):
    
    n = A.shape[0]
    e = np.ones((1,n))
    
    w1, v1 = np.linalg.eig(A)
    w2, v2 = np.linalg.eig(A.T)

    w3, v3 = np.linalg.eig(p_A)
    w4, v4 = np.linalg.eig(p_A.T)

    #take first eigen for comparison (note, may look at whole thing in a bit)
    p1 = v2[:,0]
    q1 = v1[:,0]

    p2 = v4[:,0]
    q2 = v3[:,0]

    
    S = np.ones(k) #initialize

    for i in range(0, k):
        # sine calculations
        c_i = np.dot(v1[:, i], v3[:,i]/(np.linalg.norm(v1[:, i])*np.linalg.norm(v3[:, i])))
        s_i = np.sqrt(1-c_i**2)
        S[i] = np.arcsin(s_i)*180/np.pi
            
            
    x1 = np.real(np.divide(k*p1,np.sum(p1*e)))
    y1 = np.real(np.divide(k*q1,np.sum(q1*e)))

    x2 = np.real(np.divide(k*p2,np.sum(p2*e)))
    y2 = np.real(np.divide(k*q2,np.sum(q2*e)))

    xabs = np.mean(np.abs(x1-x2))
    yabs = np.mean(np.abs(y1-y2))
    #        IR_u[i] = np.count_nonzero(x_u>1)/n    

    xpd = (np.count_nonzero(x2>1)/n) - (np.count_nonzero(x1>1)/n)
    ypd = (np.count_nonzero(y2>1)/n) - (np.count_nonzero(y1>1)/n)
    

    resx = kendalltau(x1,x2)
    resoutx = resx.statistic

    resy = kendalltau(y1,y2)
    resouty = resy.statistic

    
    return xabs, yabs, xpd, ypd, S[0], np.mean(S), resoutx, resouty


"""
o_x = x.argsort()
r_x = o_x.argsort()

o_y = y.argsort()
r_y = o_y.argsort()

plt.figure(1)
plt.scatter(r_x, r_y)
"""