# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 19:47:52 2023

@author: colej

CJ610 DSM analysis
Imports CJ610 engine edgelist
Turns into DSM
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
import csv
from numpy import linalg as LA



# %%

G = nx.read_edgelist("cj610_edgelist.txt", nodetype=int)

CJ610_DSM = nx.to_numpy_array(G)

#later, read in csv to 

# %%

n= CJ610_DSM.shape[0]

e_vals, e_vecs = LA.eig(CJ610_DSM.T)

# %%
e_vecs = np.real(e_vecs)