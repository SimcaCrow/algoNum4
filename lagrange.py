#!/usr/bin/env python
# coding: utf-8
# ---------------------------------#
"""
File : lagrange.py
Author : El-Habr C.
Description : 4ième TD algorithme numerique
"""
# ---------------------------------#

import numpy as np
from math import *
from newton import  *
import matplotlib.pyplot as plt

# ---------------------------------#

# Three Forces
elastic_force = lambda v, k: [-k * v[0], v[1]]

centrifugal_force =  lambda v, k, v0 : k * (v - v0)

gravitational_force = lambda v, k, v0 : -k * (v - v0) / pow(np.linalg.norm(v - v0), 3)

# ---------------------------------#

# Jacobians
elastic_jac = lambda k : [[-k,0], [0, 1]]

centrifugal_jac = lambda k: [[k,0], [0, k]]

gravitational_jac = lambda v, k, v0 : [[-k * ((v[1] - v0[1])**2 - 2*(v[0]-v0[0])**2) / (((v[0]-v0[0])**2 + (v[1]-v0[1])**2)**(5/2)),
                                       (3*k * (v[0]-v0[0]) * (v[1]-v0[1])) / (((v[0]-v0[0])**2 + (v[1]-v0[1])**2)**(5/2))],
                                       [(3*k * (v[0]-v0[0]) * (v[1]-v0[1])) / (((v[0]-v0[0])**2 + (v[1]-v0[1])**2)**(5/2)),
                                       -k * ((v[0] - v0[0])**2 - 2*(v[1]-v0[1])**2) / (((v[0]-v0[0])**2 + (v[1]-v0[1])**2)**(5/2))]]

# ---------------------------------#
"""
Adding 3-force test
"""
def add_force(u):
    f1 = gravitational_force(u, 1, np.array([0, 0]))
    f2 = gravitational_force(u, 0.01, np.array([1, 0]))
    f3 = centrifugal_force(u, 1, np.array([0.0099009900990099, 0]))
    return f1 + f2 + f3

# ---------------------------------#
"""
Adding 3-Jacobian test
"""
def add_jac(u):
    j1 = gravitational_jac(u, 1, np.array([0, 0]))
    j2 = gravitational_jac(u, 0.01, np.array([1, 0]))
    j3 = centrifugal_jac(1)
    return np.array(j1) + np.array(j2) + np.array(j3)

# ---------------------------------#
"""
Compute each equilibrium point to obtain Lagrangian points
"""
def lagrangian_points():
    pres = 10 ** -10
    n = 1000

    L1 = np.array([0.75, 0])
    L2 = np.array([1.15, 0])
    L3 = np.array([-1, 0])
    L4 = np.array([0.5, 0.75])
    L5 = np.array([0.5, -0.75])

    tab_points = [L1, L2, L3, L4, L5]
    names = ["L1", "L2", "L3", "L4", "L5"]
    equilibrium_point = [0,0,0,0,0]

    for i in range(len(tab_points)):
        equilibrium_point[i] = newton_raphson_back(add_force, add_jac, tab_points[i], n, pres)

    lagrangian_points_plot(equilibrium_point, names)

# ---------------------------------#
"""
Plot Lagrangian points on a graph
"""
def lagrangian_points_plot(points, names):
    plt.plot(0,0, 'yo')
    plt.text(-0.05, 0.1, 'Sun')
    plt.plot(1,0, 'bo')
    plt.text(0.9, 0.1, 'Earth')

    for i in range(len(points)):
        plt.plot(points[i][0], points[i][1], 'ro')
        plt.text(points[i][0]-0.05, points[i][1]-0.15, names[i])

    ax = plt.gca()
    circle = plt.Circle((0, 0), 1, color='k', fill=False)
    ax.add_artist(circle)
    plt.axis([-1.5, 1.5, -1.5, 1.5])
    plt.title("Lagrangian Points")
    plt.show()

# ---------------------------------#
