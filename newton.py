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
