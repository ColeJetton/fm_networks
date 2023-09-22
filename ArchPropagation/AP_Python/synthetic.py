# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 09:28:15 2023

@author: colej

synthetic.py: Python translation of synthetic.m written by ___

Note: this is computationally intensive due to translation from MATLAB,
the loops for large sample size, and my own limited skills as a programmer
"""
# import required packages and functions from support file
import numpy as np
import pandas as pd
from arch_prop_functions import *
import matplotlib.pyplot as plt
import time

beginning = time.time()
# reset random number generator for reproducability (same as MATLAB)
np.random.seed(1337)

# Number of components in the system
nodes = 256
# modular perturbation probability
# p = np.array([0.05, 0.025, 0.01, 0.005, 0.025, 0.001])
p = np.array([0.05, 0.025])#, 0.01, 0.005, 0.025, 0.001])
# Number of eigenvectors to retain
k = 100
# Number of computation loops
g = 100


# %%
# memory initialization
Pt = pd.DataFrame(index=range(1), columns=range(len(p)))
a = pd.DataFrame(index=range(len(p)), columns=range(g))
x = pd.DataFrame(index=range(len(p)), columns=range(g))
y = pd.DataFrame(index=range(len(p)), columns=range(g))
am = pd.DataFrame(np.zeros((len(p), k)))
xm = pd.DataFrame(np.zeros((len(p), k)))
ym = pd.DataFrame(np.zeros((len(p), k)))

# %%reference system
modules = 8

b = modular(nodes, modules, 0.9, int(nodes / modules), 0.5, 0.5)
b = np.triu(b.to_numpy())
b = (
    b + b.T
)  # +5*np.identity(nodes) #add values to diagonal to filter for removal

b[
    (np.arange(0, nodes), np.arange(0, nodes))
] = 0  # turn diagonal values into zeros (this should be simpler to do but idk)

# plt.imshow(b)
# plt.show()


#%% Doing the perturbations\

for i in range(0, len(p)):
    print("Working on Perturbation Probability", p[i])
    for j in range(0, g):
        Er1 = HMPerturbMatrix(b, nodes, modules, 0, p[i], 0.5, 0.5)
        Er2 = HMPerturbMatrix(b, nodes, modules, 0, p[i], 0.9, 0.1)
        Er3 = HMPerturbMatrix(b, nodes, modules, 0, p[i], 0.1, 0.9)
        a.iloc[i, j] = ComputeSine(b, Er1, k)
        x.iloc[i, j] = ComputeSine(b, Er2, k)
        y.iloc[i, j] = ComputeSine(b, Er3, k)

print("Calculations elapsed in ", time.time() - beginning, " seconds")

#%% Calculating mean perturbations per value of p
# a.iloc[0,:] = np.arange(0,100)
a = a.to_numpy()
x = x.to_numpy()
y = y.to_numpy()

# I don't think there's any point to this code because it's not used elsewhere but is in the original
for i in range(0, len(p)):
    am.iloc[i, :] = np.mean(a[i, :])
    xm.iloc[i, :] = np.mean(x[i, :])
    ym.iloc[i, :] = np.mean(y[i, :])


agrandmean = np.mean(a)
xgrandmean = np.mean(x)
ygrandmean = np.mean(y)

#%% Plotting the results
figs, axs = plt.subplots(2, 2)
axs[0, 0].imshow(b)
axs[0, 0].set_xticks([])
axs[0, 0].set_yticks([])
axs[0, 0].set_title("Unperturbed System")

axs[0, 1].scatter(np.arange(0, k), agrandmean, marker=".")
axs[0, 1].set_title("Perturbed pw=0.5, pu=0.5")

axs[1, 0].scatter(np.arange(0, k), xgrandmean, marker=".")
axs[1, 0].set_title("Perturbed pw=0.9, pu=0.1")

axs[1, 1].scatter(np.arange(0, k), ygrandmean, marker=".")
axs[1, 1].set_title("Perturbed pw=0.1, pu=0.9")

#%% Calculating and printing eigenvalues

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
    "& %0.1f"
    % (
        max(mabe - mace, mxbe - mxce, mybe - myce)
        - min(mabe - mace, mxbe - mxce, mybe - myce)
    ),
)
