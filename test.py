#!/usr/bin/env python
# coding: utf-8
# ---------------------------------#
"""
File : test.py
Author : Calandra J., Guichard A.
Description : 4ième TD algorithme numerique
"""
# ---------------------------------#

import numpy as np
import matplotlib.pyplot as plt
from newton import *

# ---------------------------------#
# Test
# ---------------------------------#
def Newton_Raphson_Test():
    f = lambda x:np.sin(x)
    J = lambda x:np.cos(x)
    f3 = lambda x:np.array([x[0]*x[0]-2,x[1]*x[1]-x[1]-1])
    J3 = lambda x:np.array([[2*x[0],0],[0,2*x[1]-1]])
    
    def f2(x):
        res = np.mat(np.zeros((2,1)))
        res[0,0]=x[0,0]**2+x[1,0]
        res[1,0]=x[0,0]-x[1,0]**2
        return res
    
    def J2(x):
        res = np.mat(np.zeros((2,2)))
        res[0,1]=1
        res[0,0]=2*x[0,0]
        res[1,0]=-2*x[1,0]
        res[1,1]=1
        return res
    
    x= np.mat(np.zeros((1,1)))
    x[0,0] = 1
    
    x2=np.mat(np.zeros((2,1)))
    x2[0,0]=9
    x2[1,0]=7
    
    x3=np.array([[2],[2]])
    """
    print(Newton_Raphson(f,J,x,1000,10**-10))
    
    print(Newton_Raphson(f2,J2,x2,1000,10**-10))

    print(Newton_Raphson(f3, J3,  x3, 1000, 10**-10))
    print(Newton_Raphson_back(f3, J3,  x3, 1000, 10**-10))
    """
# ---------------------------------#
if __name__ == '__main__':

    print("Test de Newton_Raphson...")
    Newton_Raphson_Test()
    print("Test OK")
# ---------------------------------#
