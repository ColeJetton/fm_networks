# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 13:48:17 2023

@author: colej
"""

import numpy as np

from arch_prop_functions import RandNetwork, modular, HMPerturbMatrix
from sine_inf_comparison_functions import AssymetricRandNetwork, Compare

import matplotlib.pyplot as plt

PW_DSM = np.genfromtxt('PW_DSM.csv',delimiter=',')
PW_DSM[0][0]=0 #issue with first index as nan

PW_DSM_s  = PW_DSM +PW_DSM.T

PW_DSM_s[np.where(PW_DSM_s>1)]=1

#p = [0, 0.025, 0.05, 0.00.1, 0.15, 0.2, 0.25]

p = np.arange(0,0.5,0.025 )

x_d = np.ones((len(p),len(p)))
y_d = np.copy(x_d)
xpd = np.copy(x_d)
ypd = np.copy(x_d)
S = np.copy(x_d)
Sm = np.copy(x_d)
Ktx = np.copy(x_d)
Kty = np.copy(x_d)

PW_DSM[np.where(PW_DSM > 1)] = 1


modules = 4
nodes = 256


b = modular(nodes, modules, 0.9, int(nodes / modules), 0.5, 0.5)
b = np.triu(b.to_numpy())
b = (
    b + b.T
)  # +5*np.identity(nodes) #add values to diagonal to filter for removal

b[
    (np.arange(0, nodes), np.arange(0, nodes))
] = 0  # turn diagonal values into zeros (this should be simpler to do but idk)


for i in range(0, len(p)):
    for j in range(0,len(p)):
        
         #add = AssymetricRandNetwork(PW_DSM.shape[0],p[i])#RandNetwork(25,0.1)
         #rem = AssymetricRandNetwork(PW_DSM.shape[0],p[j])

         #add = RandNetwork(PW_DSM.shape[0], p[i])
         #rem = RandNetwork(PW_DSM.shape[0], p[j])

         #p_PW_DSM = PW_DSM_s + add - rem
         #p_PW_DSM = PW_DSM + add - rem

         #p_PW_DSM[np.where(p_PW_DSM<0)] = 0 

         #p_PW_DSM[np.where(p_PW_DSM > 1)] = 1 
         add = AssymetricRandNetwork(b.shape[0],p[i])#RandNetwork(25,0.1)
         rem = AssymetricRandNetwork(b.shape[0],p[j])
         
         p_b = b + add - rem
         p_b[np.where(p_b<0)]=0
         p_b[np.where(p_b>1)]=1
         #x_d[i,j],y_d[i,j], S[i,j], Sm[i,j] = Compare(PW_DSM_s,p_PW_DSM, PW_DSM.shape[0])
         #x_d[i,j],y_d[i,j], xpd[i,j], ypd[i,j], S[i,j], Sm[i,j], Ktx[i,j],Kty[i,j] = Compare(PW_DSM,p_PW_DSM, PW_DSM.shape[0])
         x_d[i,j],y_d[i,j], xpd[i,j], ypd[i,j], S[i,j], Sm[i,j], Ktx[i,j],Kty[i,j] = Compare(b,p_b, 100)

"""
fig, ax = plt.subplots(2,2)

ax[0,0].imshow(x_d)
ax[0,0].set_title("mean diff in influence")
ax[0,0].set_xticks([])
ax[0,0].set_yticks([])
plt.colorbar(fig, ax = ax[0,0])

ax[0,1].imshow(y_d)
ax[0,1].set_title("mean diff in susceptability")
ax[0,1].set_xticks([])
ax[0,1].set_yticks([])

ax[1,0].imshow(S)
ax[1,0].set_title("angle largest eigen")
ax[1,0].set_xticks([])
ax[1,0].set_yticks([])


ax[1,1].imshow(Sm)
ax[1,1].set_title("mean angle eigen")
ax[1,1].set_xticks([])
ax[1,1].set_yticks([])
"""
# %%
q = 0
data = np.array([x_d,y_d,xpd,ypd,S,Sm,Ktx,Kty])
titles = ["mu abs diff in inf of each component", \
          "mu abs diff in sus of each component", \
          "mu change in above 1 inf", "mu change above 1 sus", \
          "angle largest eig", "mu angle eig", \
              "Krcc influence", "Krcc Suscept"]

fig, ax = plt.subplots(2,4)
for i in range(2):
    for j in range(4):
        im = ax[i,j].imshow(data[q].T, origin = 'lower')
        plt.colorbar(im, ax = ax[i,j])
        ax[i,j].set_title(titles[q])
        ax[i,j].set_xticks([])
        ax[i,j].set_yticks([])
        #ax[i,j].set_xlabel('addition percentage 0->0.25')
        #ax[i,j].set_ylabel('subtraction percentage 0->0.25')
        q = q+1
        
plt.show
#%%
"""
for i in range(2):
    for j in range(2):
        plt.colorbar(fig, ax = ax[i,j])
"""
#%%
#fig2, ax2 = plt.subplots(1,2)
#ax2[0].imshow(PW_DSM_s)
#ax2[1].imshow(p_PW_DSM)


"""
add = AssymetricRandNetwork(PW_DSM.shape[0],0.5)#RandNetwork(25,0.1)
rem = AssymetricRandNetwork(PW_DSM.shape[0],0.5)

p_PW_DSM = PW_DSM + add - rem

p_PW_DSM[np.where(p_PW_DSM<0)] = 0 

p_PW_DSM[np.where(p_PW_DSM > 1)] = 1 

PW_DSM[np.where(PW_DSM > 1)] = 1

e = np.ones((1,PW_DSM.shape[0]))

# %% keeping it here for manual checks
# back practice but I'm tired

w1, v1 = np.linalg.eig(PW_DSM)
w2, v2 = np.linalg.eig(PW_DSM.T)

w3, v3 = np.linalg.eig(p_PW_DSM)
w4, v4 = np.linalg.eig(p_PW_DSM.T)

#take first eigen for comparison (note, may look at whole thing in a bit)
p1 = v2[:,0]
q1 = v1[:,0]

p2 = v4[:,0]
q2 = v3[:,0]

#c = np.dot(q1,q2/(np.linalg.norm(q1)*np.linalg.norm(q2)))
#s = np.sqrt(1-c**2)
#S = np.arcsin(s)*180/np.pi

k = PW_DSM.shape[0]
S = np.ones(k) #initialize

for i in range(0, k):
    # sine calculations
    c_i = np.dot(v1[:, i], v3[:,i]/(np.linalg.norm(v1[:, i])*np.linalg.norm(v3[:, i])))
    s_i = np.sqrt(1-c_i**2)
    S[i] = np.arcsin(s_i)*180/np.pi\
        
        
x1 = np.real(np.divide(k*p1,np.sum(p1*e)))
y1 = np.real(np.divide(k*q1,np.sum(q1*e)))

x2 = np.real(np.divide(k*p2,np.sum(p2*e)))
y2 = np.real(np.divide(k*q2,np.sum(q2*e)))

xabs = np.mean(np.abs(x1-x2))
yabs = np.mean(np.abs(y1-y2))


"""
