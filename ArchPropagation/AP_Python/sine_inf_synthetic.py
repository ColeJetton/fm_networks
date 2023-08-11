# -*- coding: utf-8 -*-
"""
Experimentation code with comparing DSM rotation data to  infinte
regress calculations
@author: colej
"""

# import required packages and functions from support file
import numpy as np
import pandas as pd
from arch_prop_functions import modular, HMPerturbMatrix
from sine_inf_comparison_functions import ComputeSine_and_InfReg
import matplotlib.pyplot as plt
import time

beginning = time.time()
# reset random number generator for reproducability (same as MATLAB)
np.random.seed(1337)

# Number of components in the system
nodes = 256
# modular perturbation probability
p = np.array([0.05, 0.025, 0.01, 0.005, 0.0025, 0.001])
# Number of eigenvectors to retain
k = 100
# Number of computation loops
g = 50


# %%
# memory initialization
Pt = pd.DataFrame(index=range(1), columns=range(len(p)))

S1 = pd.DataFrame(index=range(len(p)), columns=range(g))
S2 = S1.copy()
S3 = S1.copy()

IRu1 = S1.copy()
IRu2 = S1.copy()
IRu3 = S1.copy()

IRp1 = S1.copy()
IRp2 = S1.copy()
IRp3 = S1.copy()

S1mean = pd.DataFrame(np.zeros((len(p), k)))
S2mean = S1mean.copy()
S3mean = S1mean.copy()

IRu1mean = S1mean.copy()
IRu2mean = S1mean.copy()
IRu2mean = S1mean.copy()

IRp1mean = S1mean.copy()
IRp2mean = S1mean.copy()
IRp3mean = S1mean.copy()

# %% synthetic reference system
modules = 8

b = modular(nodes, modules, 0.1, int(nodes / modules), 0.9, 0.1)
b = np.triu(b.to_numpy())
b = b + b.T  #add values to diagonal to filter for removal

b[(np.arange(0, nodes), np.arange(0, nodes))] = 0  
plt.imshow(b)

# %%


for i in range(0, len(p)):
    print("Working on Perturbation Probability", p[i])
    for j in range(0, g):
        Er1 = HMPerturbMatrix(b, nodes, modules, 0, p[i], 0.2, 0.8)
        Er2 = HMPerturbMatrix(b, nodes, modules, 0, p[i], 0.3, 0.7)
        Er3 = HMPerturbMatrix(b, nodes, modules, 0, p[i], 0.4, 0.6)
        S1.iloc[i, j], IRu1.iloc[i,j], IRp1.iloc[i,j] = ComputeSine_and_InfReg(b, Er1, k)
        S2.iloc[i, j], IRu2.iloc[i,j], IRp2.iloc[i,j] = ComputeSine_and_InfReg(b, Er2, k)
        S3.iloc[i, j], IRu3.iloc[i,j], IRp3.iloc[i,j] = ComputeSine_and_InfReg(b, Er3, k)

print("Calculations elapsed in ", time.time() - beginning, " seconds")


# %% calculating means

S1 = S1.to_numpy()
S2 = S2.to_numpy()
S3 = S3.to_numpy()

IRu1 = IRu1.to_numpy()
IRu2 = IRu2.to_numpy()
IRu3 = IRu3.to_numpy()

IRp1 = IRp1.to_numpy()
IRp2 = IRp2.to_numpy()
IRp3 = IRp3.to_numpy()
# %%
S1_mean = np.mean(S1); S1_gm = np.mean(S1_mean)
S2_mean = np.mean(S2); S2_gm = np.mean(S2_mean)
S3_mean = np.mean(S3); S3_gm = np.mean(S3_mean)

IRu1_mean = np.mean(IRu1); IRu1_gm = np.mean(IRu1_mean)
IRu2_mean = np.mean(IRu2); IRu2_gm = np.mean(IRu2_mean)
IRu3_mean = np.mean(IRu3); IRu3_gm = np.mean(IRu3_mean)

IRp1_mean = np.mean(IRp1); IRp1_gm = np.mean(IRp1_mean)
IRp2_mean = np.mean(IRp2); IRp2_gm = np.mean(IRp2_mean)
IRp3_mean = np.mean(IRp3); IRp3_gm = np.mean(IRp3_mean)

print('1', S1_gm, IRu1_gm, IRp1_gm, IRp1_gm-IRu1_gm)
print('2', S2_gm, IRu2_gm, IRp2_gm, IRp2_gm-IRu2_gm)
print('3', S3_gm, IRu3_gm, IRp3_gm, IRp3_gm-IRu3_gm)
# %%


"""
# %% plotting???

figs, axs = plt.subplots(3, 4)

axs[0,0].scatter(np.arange(0, k), S1_mean, marker=".")
axs[0,0].set_title("Eigenvector Rotation Angle Perturbed pw=0.5, pu=0.5")
axs[0,0].set_xticks([])

axs[0,1].scatter(np.arange(0, k), IRu1_mean, marker=".")
axs[0,1].set_title("Symmetric Infinite Regress Percentage above 1 (unperturbed matrix)")
axs[0,1].set_xticks([])

axs[0,2].scatter(np.arange(0, k), IRp1_mean, marker=".")
axs[0,2].set_title("Symmetric Infinite Regress Percentage above 1 (perturbed matrix)")
axs[0,2].set_xticks([])


axs[0,3].scatter(np.arange(0, k), IRp1_mean - IRu1_mean, marker=".")
axs[0,3].set_title("Symmetric Infinite Regress Percentage above 1 (perturbed matrix)")
axs[0,3].set_xticks([])


axs[1,0].scatter(np.arange(0, k), S2_mean, marker=".")
axs[1,0].set_title("Eigenvector Rotation Angle Perturbed pw=0.5, pu=0.5")
axs[1,0].set_xticks([])

axs[1,1].scatter(np.arange(0, k), IRu2_mean, marker=".")
axs[1,1].set_title("Symmetric Infinite Regress Percentage above 1 (unperturbed matrix)")
axs[1,1].set_xticks([])

axs[1,2].scatter(np.arange(0, k), IRp2_mean, marker=".")
axs[1,2].set_title("Symmetric Infinite Regress Percentage above 1 (perturbed matrix)")
axs[1,2].set_xticks([])

axs[1,3].scatter(np.arange(0, k), IRp2_mean - IRu3_mean, marker=".")
axs[1,3].set_title("Symmetric Infinite Regress Percentage above 1 (perturbed matrix)")
axs[1,3].set_xticks([])


axs[2,0].scatter(np.arange(0, k), S3_mean, marker=".")
axs[2,0].set_title("Eigenvector Rotation Angle Perturbed pw=0.5, pu=0.5")
axs[2,0].set_xticks([])

axs[2,1].scatter(np.arange(0, k), IRu3_mean, marker=".")
axs[2,1].set_title("Symmetric Infinite Regress Percentage above 1 (unperturbed matrix)")
axs[2,1].set_xticks([])

axs[2,2].scatter(np.arange(0, k), IRp3_mean, marker=".")
axs[2,2].set_title("Symmetric Infinite Regress Percentage above 1 (perturbed matrix)")
axs[2,2].set_xticks([])

axs[2,3].scatter(np.arange(0, k), IRp3_mean - IRu3_mean, marker=".")
axs[2,3].set_title("Symmetric Infinite Regress Percentage above 1 (perturbed matrix)")
axs[2,3].set_xticks([])


# %% calculating variance????

S1_grand_var = np.var(S1)
S2_mean = np.mean(S2)
S3_mean = np.mean(S3)

IRu1_grand_var = np.var(IRu1)
IRu2_mean = np.mean(IRu2)
IRu3_mean = np.mean(IRu3)

IRp1_grand_var = np.var(IRp1)
IRp2_mean = np.mean(IRp2)
IRp3_mean = np.mean(IRp3)

# %% plotting var
figs, axs = plt.subplots(3, 1)

axs[0].scatter(np.arange(0, k), S1_grand_var, marker=".")
axs[0].set_title("Eigenvector Rotation Angle Perturbed pw=0.5, pu=0.5")

axs[1].scatter(np.arange(0, k), IRu1_grand_var, marker=".")
axs[1].set_title("Symmetric Infinite Regress Percentage above 1 (unperturbed matrix)")

axs[2].scatter(np.arange(0, k), IRp1_grand_var, marker=".")
axs[2].set_title("Symmetric Infinite Regress Percentage above 1 (perturbed matrix)")

# %% Plotting change???

figs, axs = plt.subplots(2, 1)

axs[0].scatter(np.arange(0, k), S1_mean, marker=".")
axs[0].set_title("Eigenvector Rotation Angle Perturbed pw=0.5, pu=0.5")
axs[0].set_xticks([])

axs[1].scatter(np.arange(0, k), IRu1_mean - IRp1_mean, marker=".")
axs[1].set_title("Symmetric Infinite Regress Percentage above 1 (unperturbed matrix)")

"""
# %% mean differences????
