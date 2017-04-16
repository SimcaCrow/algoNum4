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
"""
Calculation of the energy gradient
X : a position vector
Return a vector
"""
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

# ---------------------------------#
"""
Jacobian calculation of the energy gradient
Terms i,j are calculated
X : a position vector
"""

def jacobian_grad_E_ij(X, i, j):
    N = X.size
    if ( i == j ):
        jac__ij = 0
        jac__ij -= 1.0/((X[i,0] + 1)*(X[i,0] + 1))
        jac__ij -= 1.0/((X[i,0] - 1)*(X[i,0] - 1))
        for k in range (N):
            if ( k != i ):
                jac__ij += 1.0/((X[i,0] - X[k,0])*(X[i,0] - X[k,0]))
    else:
        jac__ij = -1.0/((X[i,0] - X[j,0])*(X[i,0] - X[j,0]))
    return jac__ij

def jacobian_grad_E(X):
    N = X.size
    J = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            J[i,j] = jacobian_grad_E_ij(X, i, j)
    return J

 

# ---------------------------------#
"""
Calculation of the equilibrium position for n charges
Return a position vector
"""
def pos_equilibrium(X):
    N = 1000
    eps = 10**-5
    U = newton_raphson(grad_E, jacobian_grad_E, X, N, eps) #TODO : add jacobian
    return U

def exchange(A,i,j):
    tmp = A[i]
    A[i] = A[j]
    A[j] = tmp
    
def mirror(A):
    n = A.size
    for i in range(n/2):
        exchange(A,i,n-i-1)
    return A

# ---------------------------------#
"""
Plot Legendre polynomials with color and label for rank n
And equilibrium points of energy
"""
def add_plot(X, lbl, clr, type='o'):
    Peq = pos_equilibrium(X)
    n= Peq.size
    for i in range(n):
        plt.plot(Peq[i,0],0,type, color=clr)
    
    c = [0]*(n+2)
    c[n+1] = 1
    
    d = L.legder(c)
    P = L.leg2poly(d)

    P = mirror(P)
    Poly = np.poly1d(P)
    x = np.linspace(-1,1,100)
    y = Poly(x)
    plt.plot(x, y, label=lbl, color=clr)

# ---------------------------------#

