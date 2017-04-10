#!/usr/bin/env python
# coding: utf-8
# ---------------------------------#
"""
File : electrostatic_equilibrium.py
Author : Calandra J.
Description : 4ième TD algorithme numerique
"""
# ---------------------------------#

from newton import *
import numpy as np

# ---------------------------------#

# Calcul du gradient de l'energie
# X : un vecteur de positions
# renvoie un vecteur
def grad_E(X):
    n, m = np.shape(X)
    assert(m == 1)
    res = np.zeros((n,1))
    for i in range(n):
        s = 0
        for j in range(n):
            if (i != j):
                s += 1./(X[i,0] - X[j,0])
        res[i, 0] = 1./X[i,0]+1 + 1./X[i,0]-1 + s
    return res

# Calcul de la position d'équilibre pour
# n charges
# retourne le vecteur des positions
def pos_equilibre(n):
    U0 = np.zeros((n,0))
    N = 1000
    eps = 10**-5
    U = Newton_Raphson_back(grad_E, J, U0, N, eps)
    return U
