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

def Newton_Raphson(f, J, U0, N, epsilon):
    U = np.copy(U0)
    for i in range(N):
        if np.all(abs(f(U)) < epsilon):
            return U
        fu = f(U)
        ju = J(U)
        tmp = np.linalg.lstsq(ju, -fu)
        V = tmp[0]
        U = U+V
    assert(False and "precision non atteinte")

def Newton_Raphson_back(f, J, x0, N, epsilon):
    x = x0
    for i in range(N):
        print(i)
        va = f(x)
        na = np.linalg.norm(va)
        if na < epsilon:
            return x
        dv = J(x)
        dx = np.linalg.lstsq(dv,-va)[0]
        st = 1.0
        while np.linalg.norm(f(x+st*dx)) >= na+0.001:
            st *= 2.0/3.0
        x += st*dx
    assert(False and "precision non atteinte")
