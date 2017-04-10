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
from electrostatic_equilibrium import *

# ---------------------------------#
# Test
# ---------------------------------#
def Newton_Raphson_Test():
    pres = 10**-10
    
    f = lambda x:np.sin(x)
    J = lambda x:np.cos(x)
   
    f2 = lambda x:np.array([x[0]*x[0]-2,x[1]*x[1]-x[1]-1])
    J2 = lambda x:np.array([[2*x[0],0],[0,2*x[1]-1]])
        
    x= np.mat([1.])
    print x
    x2=np.array([[2],[2]])

    print("\nTest de Newton_Raphson 1")
    U = Newton_Raphson(f, J, x, 1000, pres)
    assert(np.linalg.norm(f(U)) < pres)
    print("Test OK")
    
    print("\nTest de Newton_Raphson 2")
    U = Newton_Raphson(f2, J2, x2, 1000, pres)
    assert(np.linalg.norm(f2(U)) < pres)
    print("Test OK")

    print("\nTest de Newton_Raphson_back")
    print("avec les memes param que pour le test precedent")
    U = Newton_Raphson_back(f2, J2, x2, 1000, pres)
    assert(np.linalg.norm(f2(U)) < pres)
    print("Test OK")

def Elec_Equ_Test():
    A = np.matrix([[0.2]])
    
    B = np.matrix([[0.4],
                   [-0.4],
                   [0.5],
                   [0.6]])
    
    C = np.matrix([[0.5],
                   [0.6]])
    
    W = np.matrix([[0.4],
                   [-0.5],
                   [0.7]])
    
    U = np.matrix([[0.1],
                   [-0.5],
                   [0.3],
                   [-0.9],
                   [-0.8]])


    
    x = np.linspace(-1,1)
    plt.plot([-1,1], [0,0], 'k')
    plt.axis([-1,1,-15,15])
    plt.show()
    
# ---------------------------------#
if __name__ == '__main__':

    Newton_Raphson_Test()
    Elec_Equ_Test()
# ---------------------------------#
