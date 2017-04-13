#!/usr/bin/env python
# coding: utf-8
# ---------------------------------#
"""
File : test.py
Author : Calandra J., El-Habr C., Guichard A.
Description : 4ième TD algorithme numerique
"""
# ---------------------------------#

import numpy as np
import matplotlib.pyplot as plt
from newton import *
from lagrange import *
from electrostatic_equilibrium import *

# ---------------------------------#
# Tests
# ---------------------------------#
"""
Newton Raphson Test
"""
def newton_raphson_test():
    pres = 10**-10
    
    f = lambda x:np.sin(x)
    J = lambda x:np.cos(x)
    
    f2 = lambda x:np.array([x[0]*x[0]-2,x[1]*x[1]-x[1]-1])
    J2 = lambda x:np.array([[2*x[0],0],[0,2*x[1]-1]])
        
    x= np.mat([1.])
    print("x =", x)
    x2=np.array([[2],[2]])

    print("\nTest Newton_Raphson 1")
    U = newton_raphson(f, J, x, 1000, pres)
    assert(np.linalg.norm(f(U)) < pres)
    print("Test OK")
    
    print("\nTest Newton_Raphson 2")
    U = newton_raphson(f2, J2, x2, 1000, pres)
    assert(np.linalg.norm(f2(U)) < pres)
    print("Test OK")

    print("\nTest Newton_Raphson_back")
    print("With same parameters than previously")
    U = newton_raphson_back(f2, J2, x2, 1000, pres)
    assert(np.linalg.norm(f2(U)) < pres)
    print("Test OK\n")

# ---------------------------------#
"""
Electrostatic Equilibrium Test
"""
def elec_equ_test():
    A = np.matrix([[0.2]])

    B = np.matrix([[0.5],
                   [0.6]])
    
    C = np.matrix([[0.4],
                   [-0.5],
                   [0.7]])

    D = np.matrix([[0.4],
                   [-0.4],
                   [0.5],
                   [0.6]])

    add_plot(A, 'n=1', 'g')
    add_plot(B, 'n=2', 'b')
    add_plot(C, 'n=3', 'y')
    add_plot(D, 'n=4', 'r')

    plt.plot([-1,1], [0,0], 'k')
    plt.axis([-1,1,-15,15])
    plt.legend()
    plt.title("Legendre polynomials and equilibrium positions")
    plt.show()

# ---------------------------------#
"""
Display 3 forces for a given k, v, v0
"""
def test_force(k, v, v0):
    print("Elastic force for k = ", k, " : ", elastic_force(v, k))
    print("Centrifugal force for k = ", k, "and v0 = ", v0, " : ", centrifugal_force(v, k, v0))
    print("Gravitational force for k = ", k, "and v0 = ", v0, " : ", gravitational_force(v, k, v0), "\n")

# ---------------------------------#
"""
Display jacobians for a given k, v, v0
"""
def test_jac(k, v, v0):
    print("Elastic Jacobian for k = ", k, " : ", elastic_jac(k))
    print("Centrifugal Jacobian for k = ", k, "and v0 = ", v0, " : ", centrifugal_jac(k))
    print("Gravitational Jacobian for k = ", k, "and v0 = ", v0, " : ", gravitational_jac(v, k, v0), "\n")

# ---------------------------------#
"""
Test all calculation functions on forces
"""
def forces_test():
    k = 1
    u = np.array([1.5, 0])
    v0 = np.array([0,0])

    print("Display forces examples")
    test_force(k, u, v0)
    print("Display Jacobian examples")
    test_jac(k, u, v0)

    print("Test forces sum")
    force = np.array([1.00565457, 0])
    own_force = add_force(u)
    print("Expected f(U)=", force, "Our f(U) =", own_force)
    np.testing.assert_almost_equal(force, own_force,8) # Erreur d'arrondi machine au dela de 8 decimales
    print("Test OK\n")

# ---------------------------------#
"""
Test calculation functions on jacobian
"""
def jacobians_test():
    print("Test Jacobian sum")
    u = np.array([1.5, 0])
    jac = np.array([[1.75259259, 0], [0, 0.6237037]])
    own_jac = add_jac(u)
    print("Expected df(U)=", jac, "Our df(U) =", own_jac)
    np.testing.assert_almost_equal(jac, own_jac, 8)
    print("Test OK\n")

# ---------------------------------#
if __name__ == '__main__':

    newton_raphson_test()
    # elec_equ_test()
    forces_test()
    jacobians_test()
    lagrangian_points()

# ---------------------------------#
