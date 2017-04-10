#!/usr/bin/env python
# coding: utf-8
# ---------------------------------#
"""
File : electrostatic_equilibrium.py
Author : Calandra J.
Description : 4ième TD algorithme numerique
"""
# ---------------------------------#

import numpy as np

# ---------------------------------#

def Grad_E(X):
    n, m = np.shape(X)
    assert(m == 1)
    res = np.zeros((n,1))
    for i in range(n):
        s = 0
        for j in range(n):
            if (i != j):
                s += 1./(X[i,0] - X[j,0])
        res[i, 0] = 1./(X[i,0]+1) + 1./(X[i,0]-1) + s
    return res
