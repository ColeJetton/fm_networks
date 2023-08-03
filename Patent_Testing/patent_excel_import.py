# -*- coding: utf-8 -*-
"""
@author: colej

purpose of this is to experiment with the upload and processing of component
to function mapping.
"""

import pandas as pd
import numpy as np

# import csv and remove nan/blank spaces
test_csv = pd.read_csv(
    "Sample Functional Model Data - Component x Function Table.csv"
)

test_csv = test_csv.fillna(0)

# remove extra columns and rows
col_index = []
for col in test_csv.columns:
    if "Unnamed" in col:
        col_index.append(col)


test_csv = test_csv.drop(columns=col_index)

row_index = []
for i in test_csv.index:
    if test_csv.iloc[i][0] == 0:
        row_index.append(i)

test_csv = test_csv.drop(index=row_index)

# %% creating the network-of-networks representation
test_numpy = test_csv.to_numpy()

CF = test_numpy[2:,3:].astype(int)

C = np.matmul(CF,CF.T)
C[np.where(C>1)]=1

import matplotlib.pyplot as plt

plt.figure(1)
plt.imshow(C)


plt.figure(2)
plt.imshow(CF)


plt.figure(3)
plt.imshow(np.vstack((C,CF.T)))