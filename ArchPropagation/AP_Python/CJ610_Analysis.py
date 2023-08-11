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

NOTE: It takes 34 seconds to do the infinite regress function, so about 
15-20 seconds per eigenvector calculation.

meaning..
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
from infinite_regress_functions import *


# %%

G = nx.read_edgelist("cj610_edgelist.txt", nodetype=int)

CJ610_DSM = nx.to_numpy_array(G)

# time testing for eigenvector
c_t = time.time()


x,y = Infinite_Regress(CJ610_DSM)

c_f = time.time()

print(c_f - c_t, "seconds")