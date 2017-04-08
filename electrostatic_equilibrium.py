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

def Jacobian_Delta_E(X):
    n, m = np.shape(X)
    assert(m == 1)
    res = np.mat(np.zeros((n, n)))
    for i in range(n):
        for j in range(n):
            if i != j:
                res[i, j] = -1.0/((abs(x[i, 0] - x[j, 0])**2)*np.log(10.0))
            else:
                res[i, j] -= 1.0/(abs(x[i, 0]+1)**2) + 1.0/(abs(x[i, 0]-1)**2)
                for k in range(0,n):
                    if k != i:
                        res[i, j] -= 1.0/(abs(x[i, 0] - x[k, 0])**2)
                        res[i, j] *= 1.0/np.log(10.0)
            return res
