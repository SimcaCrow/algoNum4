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
Question 1
H(U)*V = -F(V)
"""
def Newton_Raphson(f, J, U0, N, epsilon):
    U = np.copy(U0)
    for i in range(N):
        fu = f(U)
        na = np.linalg.norm(fu)
        if (na < epsilon):
            return U
        ju = J(U)
        V = np.linalg.lstsq(ju, -fu)[0] 
        U = U + V
    print("ERREUR : precision non atteinte")
    return U

def Newton_Raphson_back(f, J, U0, N, epsilon):
    U = np.copy(U0)
    for i in range(N):
        fu = f(U)
        na = np.linalg.norm(fu)
        if (na < epsilon):
            return U
        ju = J(U)
        V = np.linalg.lstsq(ju,-fu)[0]
        
        st = 1.0
        while (np.linalg.norm(f(U+st*V)) >= na + 0.001):
            st *= 1.0/3.0
        U = U + st*V
    print("ERREUR : precision non atteinte")
    return U
