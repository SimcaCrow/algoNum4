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
import matplotlib.pyplot as plt
import numpy.polynomial.legendre as L

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
def pos_equilibre(X):
    N = 1000
    eps = 10**-5
    U = Newton_Raphson_back(grad_E, JACOBIEN, X, N, eps) #TODO : rajouter le jacobien
    return U

def echange(A,i,j):
    tmp = A[i]
    A[i] = A[j]
    A[j] = tmp
    
def miroir(A):
    n = A.size
    for i in range(n/2):
        echange(A,i,n-i-1)
    return A

# Plot le polynome de legendre avec la couleur et le label pour le rang n
# Ainsi que les points d'equilibre d'energie
def add_plot(X, lbl, clr):
    n = X.size
    R = pos_equilibre(X)
    z = np.zeros(n)
    plt.plot(R, z, type='o', color=clr)
    
    c = [0]*(n+2)
    c[n+1] = 1
    
    d = L.legder(c)
    P = L.leg2poly(d)
    
    P = miroir(P)
    Poly = np.poly1d(P)
    x = np.linspace(-1,1,100)
    y = Poly(x)
    plt.plot(x, y, label=lbl, color=clr)

"""
Notes pour le rapport :
Quand tu lances le test.py, il y a les tests de newton-Raphson puis
ceux de ce fichier. Ca te sort un graphique avec les points d'equilibre pour
1,2,3 et 4 charges et les courbes des polynomes de Legendre. (Normalement) les*
polynomes de Legendre s'annulent aux memes endroit que les positions 
d'equilibre des charges
"""
