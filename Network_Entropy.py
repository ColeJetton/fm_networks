# -*- coding: utf-8 -*-
"""
Shannon Entropy In Network Calculation
@author: colej
Adapted from: 
https://stackoverflow.com/questions/70858169/networkx-entropy-of-subgraphs-generated-from-detected-communities
"""

import networkx as nx
import numpy as np
import scipy


def degree_distribution(G):
    vk = dict(G.degree())
    vk = list(vk.values()) #get degree values
    maxk = np.max(vk)
    # mink = np.min(min)
    kvalues = np.arange(0, maxk+1)
    Pk = np.zeros(maxk+1) #P(k)
    
    for k in vk:
        Pk[k] = Pk[k]+1 
    Pk = Pk/sum(Pk)
    
    return kvalues, Pk


def shannon_entropy(G):
    k, Pk = degree_distribution(G)
    H = 0
    for p in Pk:
        if p > 0:
            H = H - p*np.log(p)
            
    return H

G = nx.karate_club_graph()

#G = nx.davis_southern_women_graph()

test = shannon_entropy(G)

t_k, t_Pk = degree_distribution(G)

# %% other way of doing it, just to make it easier

def network_entropy(G):
    vk = dict(G.degree())
    vk = list(vk.values()) #get degree values
    maxk = np.max(vk)
    Pk = np.zeros(maxk+1) #P(k)
    
    for k in vk:
        Pk[k] = Pk[k]+1 
    
    Pk = Pk/sum(Pk)
    H = scipy.stats.entropy(Pk)
    
    return H
    
test2 = network_entropy(G)