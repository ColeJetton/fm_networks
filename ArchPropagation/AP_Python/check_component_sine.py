# -*- coding: utf-8 -*-
"""
one component at a time infinte regress to 
@author: colej
"""
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from arch_prop_functions import *
import matplotlib.pyplot as plt
import time
from numpy import genfromtxt
from numpy import linalg as LA

from infinite_regress_functions import *
PW_DSM = genfromtxt('PW_DSM.csv',delimiter=',')
PW_DSM[0][0]=0 #issue with first index as nan

PW_DSM_s  = PW_DSM +PW_DSM.T

PW_DSM_s[np.where(PW_DSM_s>1)]=1


x,y = Infinite_Regress(PW_DSM_s)

sine_data = np.ones((PW_DSM_s.shape[0]))

for i in range(0, PW_DSM_s.shape[1]):
    PW_DSM_cr = np.copy(PW_DSM_s) # copy
    PW_DSM_cr[:,i]=0; PW_DSM_cr[0,i]=0# remove a component to test sine change
    sine_data[i] = np.mean(ComputeSine(PW_DSM_cr,PW_DSM_s,PW_DSM_s.shape[0]))
    
x = x[~np.isnan(sine_data)]


sine_data = sine_data[~np.isnan(sine_data)]

plt.scatter(sine_data,x)

# %%
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
names = list(mapping.values())
for i in range(0,54):
    plt.annotate(names[i],(sine_data[i],x[i]))



# %%

IR_data = np.ones((PW_DSM_s.shape[0],PW_DSM_s.shape[0]))
sine_data = np.ones((PW_DSM_s.shape[0]))

for i in range(0, PW_DSM_s.shape[1]):
    PW_DSM_cr = np.copy(PW_DSM_s) # copy
    PW_DSM_cr[:,i]=0; PW_DSM_cr[0,i]=0# remove a component to test sine change
    sine_data[i] = np.mean(ComputeSine(PW_DSM_cr,PW_DSM_s,PW_DSM_s.shape[0]))
    x,y = Infinite_Regress(PW_DSM_cr)
    IR_data[i,:] = x

IR_data = np.mean(IR_data, axis=0)



IR_data = IR_data[~np.isnan(sine_data)]
sine_data = sine_data[~np.isnan(sine_data)]
plt.scatter(sine_data,IR_data)