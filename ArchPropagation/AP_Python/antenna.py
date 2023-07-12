# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 19:07:30 2023

@author: colej
"""

#load required packages
import numpy as np
import pandas as pd
from arch_prop_functions import *
import matplotlib.pyplot as plt
from numpy import linalg as LA
import time

#constants:
k = 16 #number of eigenvectors to retain
g = 100 #number of computation loops

#load connectivity and convert to component DSM
directlikelihood = pd.read_excel('antennacpm.xlsx', \
                                 sheet_name = 'directlikelihood', \
                                     header=None)
    
b = directlikelihood
b = np.triu(b.to_numpy())
b[np.where(b>0)]=1
b = b+b.T
nodes = b.shape[0]
b[(np.arange(0,nodes),np.arange(0,nodes))]=0

#initialize memory
a = pd.DataFrame(index=range(1),columns = range(g))
x = a.copy()
y = a.copy()
z= a.copy()


#create perturbationation matrix
#engines is component 3
#egine auxiliaries ic omponent 10
#auxiliary electrics is component 13
#cabling 
Er1=AntennaPerturbMatrix(b,3,.1,.1)
Er2=AntennaPerturbMatrix(b,3,.5,.5)
Er3=AntennaPerturbMatrix(b,3,.9,.1)
Er4=AntennaPerturbMatrix(b,3,.1,.9)

for i in range(0,g):
    Er1=AntennaPerturbMatrix(b,3,.1,.1)
    Er2=AntennaPerturbMatrix(b,3,.5,.5)
    Er3=AntennaPerturbMatrix(b,3,.9,.1)
    Er4=AntennaPerturbMatrix(b,3,.1,.9)
    a.iloc[0,i] =ComputeSine(b,Er1,k)
    x.iloc[0,i] =ComputeSine(b,Er2,k)
    y.iloc[0,i] =ComputeSine(b,Er3,k)
    z.iloc[0,i] =ComputeSine(b,Er4,k)

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

