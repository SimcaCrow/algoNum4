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

# ---------------------------------#

def elastic_force(k):
    return lambda v: [-k * v[0], v[1]]

# ---------------------------------#

def centrifugal_force(k, v0):
    return lambda v : [k * (v[0] - v0[0]), k * (v[1] - v0[1])]

# ---------------------------------#

def gravitational_force(k, v0):
    x = -k * ((v[0]-v0[0]) / pow(((v[0]-v0[0])**2 + (v[1]-v0[1])**2), 3/2))
    y = -k * ((v[1]-v0[1]) / pow(((v[0]-v0[0])**2 + (v[1]-v0[1])**2), 3/2))
    return lambda v : [x,y]

# ---------------------------------#

def elastic_jac(k):
    return lambda v : [[-k,0], [0, 1]]

# ---------------------------------#

def centrifugal_jac(k):
    return lambda v : [[k,0], [0, k]]

# ---------------------------------#

def gravitational_jac(k, v0):
    xx = -k * ((v[1] - v0[1])**2 - 2*(v[0]-v0[0])**2) / (((v[0]-v0[0])**2 + (v[1]-v0[1])**2)**(5/2))
    yx = (3*k * (v[0]-v0[0]) * (v[1]-v0[1])) / (((v[0]-v0[0])**2 + (v[1]-v0[1])**2)**(5/2))
    xy = (3*k * (v[0]-v0[0]) * (v[1]-v0[1])) / (((v[0]-v0[0])**2 + (v[1]-v0[1])**2)**(5/2))
    yy = -k * ((v[0] - v0[0])**2 - 2*(v[1]-v0[1])**2) / (((v[0]-v0[0])**2 + (v[1]-v0[1])**2)**(5/2))
    return lambda v: [[xx, yx],[xy, yy]]

# ---------------------------------#

def newton_on_force(k, v0):
    u = [1.5,0]
    f = [1.00565457,0]

# ---------------------------------#

def start_test(name):
    print("...", name, "\n#" + '-' * 25 + "#")

# ---------------------------------#

def end_test(name):
    print("#" + '-' * 25 + "#\n", name ,"... Ok\n")

# ---------------------------------#

def test_force(k, v, v0):
    start_test(test_force.__name__)
    print("Elastic force for k = ", k, " : ", elastic_force(k)(v))
    print("Centrifugal force for k = ", k, "and v0 = ", v0, " : ", centrifugal_force(k, v0)(v))
    print("Gravitational force for k = ", k, "and v0 = ", v0, " : ", gravitational_force(k, v0)(v))
    end_test(test_force.__name__)

# ---------------------------------#

def test_jac(k, v, v0):
    start_test(test_jac.__name__)
    print("Elastic Jacobian for k = ", k, " : ", elastic_jac(k)(v))
    print("Centrifugal Jacobian for k = ", k, "and v0 = ", v0, " : ", centrifugal_jac(k)(v))
    print("Gravitational Jacobian for k = ", k, "and v0 = ", v0, " : ", gravitational_jac(k, v0)(v))
    end_test(test_jac.__name__)

# ---------------------------------#
if __name__ == '__main__':

    k = 1
    v =  np.array([1, 2]) # FIX ME
    v0 = np.array([0,0])
    test_force(k, v, v0)
    test_jac(k, v, v0)

    k = 0.01
    v = np.array([1, 2])
    v0 = np.array([1, 0])
    test_force(k, v, v0)
    test_jac(k, v, v0)

# ---------------------------------#
