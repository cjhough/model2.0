#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 16:30:45 2017

@author: carly
"""
import numpy as np
import matplotlib.pyplot as plt

def gain_pl_sqrt(x):
    """gain function"""
    y = np.sqrt((x>0)*x)
    return y

U = np.random.uniform(0,1,(2,2))
A = np.random.uniform(0,1,(2,1))
Z = np.random.uniform(0,1,(2,1))

#
dt = 1
sigma = 2.3
alpha = 10
beta = 2.6
gamma = 1000
tau_u = 0.01
tau_a = 1000

S1 = 1
S2 = 1

u1 = []
u2 = []
a1 = []
a2 = []
z1 = []
z2 = []

Z[0] += -dt*Z[0] +np.sqrt(dt)*sigma*np.random.randn()
Z[1] += -dt*Z[1] +np.sqrt(dt)*sigma*np.random.randn()

U[0][1] = U[1] +dt*(gain_pl_sqrt(S1 - beta*A[1]*U[1][0] + alpha*U[0]))/tau_u
U[1][1] = U[1][0] +dt*(gain_pl_sqrt(S2 - beta*A[0]*U[0][0] + alpha*U[1]))/tau_u
U[0][0]= U[0][1]
U[1][0]= U[1][1]



A[0] += dt*(1 - A[0] - gamma*A[0]*U[0][1])/tau_a
A[1] += dt*(1 - A[1] - gamma*A[1]*U[1][1])/tau_a

u1.append(U[0][0])
u2.append(U[1][0])
a1.append(A[0])
a2.append(A[1])
z1.append(Z[0])
z2.append(Z[1])

