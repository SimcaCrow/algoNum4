#!/usr/bin/env python
# coding: utf-8
# ---------------------------------#
"""
File : newton.py
Author : Guichard A.
Description : 4ième TD algorithme numerique
"""
# ---------------------------------#

import numpy as np

# ---------------------------------#
"""
Newton Raphson Implementation
"""
def newton_raphson(f, J, U0, N, epsilon):
    U = np.copy(U0)
    for i in range(N):
        fu = f(U)
        na = np.linalg.norm(fu)
        if (na < epsilon):
            return U
        ju = J(U)
        V = np.linalg.lstsq(ju, -fu)[0] 
        U = U + V
    print("ERROR : precision not reached")
    return U

# ---------------------------------#
"""
Newton Raphson Implementation with Backtracking
"""
def newton_raphson_back(f, J, U0, N, epsilon):
    U = np.copy(U0)
    for i in range(N):
        fu = f(U)
        na = np.linalg.norm(fu)
        if (na < epsilon):
            return U
        ju = J(U)
        V = np.linalg.lstsq(ju,-fu)[0]
        
        if(np.linalg.norm(f(U+V)) - na >= 0):
            V = (1.0/3.0) * V
        U = U + V
    print("ERROR : precision not reached")
    return U

# ---------------------------------#
