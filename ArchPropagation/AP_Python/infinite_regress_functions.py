# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 16:31:00 2023

@author: colej

Improvement/Generalization of Infinite Regress Analysis
"""
import numpy as np
from numpy import linalg as LA


def Infinite_Regress(A):
    """

    Parameters
    ----------
    A : Assymetric Design Structure Matrix

    Returns
    -------
    x: Influence Score of a  component
    y: Susceptability Score of a component

    """
    n = A.shape[0]  # number of components
    e = np.ones((1, n))

    # find eigenvectors of the DSM and return the largests for calculations
    v1, p = LA.eig(A.T)
    v2, q = LA.eig(A)
    p = p[:, 0]
    q = q[:, 0]

    # calculate influence and susceptability scores
    x = np.divide(n*p, np.sum(p*e))
    y = np.divide(n*q, np.sum(q*e))
    x = np.real(x)
    y = np.real(y)
    return x, y
