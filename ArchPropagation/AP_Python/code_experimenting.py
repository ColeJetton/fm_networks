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
p = np.array([0.05, 0.025, 0.01, 0.005, 0.025, 0.001])
# Number of eigenvectors to retain
k=100;
# Number of computation loops
g=100;
a = pd.DataFrame(index=range(len(p)),columns = range(g))

b = np.random.rand(10,10)
b = np.array([[1, 2, 3],[4, 5, 6],[7,8,9]])
w,v=LA.eig(b)
