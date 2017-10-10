#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 16:30:45 2017

@author: carly
"""
import numpy as np

def gain_pl_sqrt(x):
    """gain function"""
    y = np.sqrt((x>0)*x)
    return y

U = np.random.uniform(0,1,(2,2))
A = np.random.uniform(0,1,(2,1))
Z = np.random.uniform(0,1,(2,1))

p = {}
p['beta']=1.42528
p['alpha']=1
p['gamma'] = 2.07703
p['tau_a']=1248.0
p['tau_u']=10
p['sigma']=0.0419821
p['mu']=0.9
p['ks']=0.1918590
p['kc']=0.00253604
p['kb'] = 0.109369


u1 = []
u2 = []
a1 = []
a2 = []
z1 = []
z2 = []

#full
Z[0,0] = Z[0,0] -dt*(Z[0,0] + np.sqrt(dt)*sigma*np.random.randn())
Z[1,0] = Z[1,0] -dt*(Z[1,0] + np.sqrt(dt)*sigma*np.random.randn())

U[0,1] = U[0,0]+dt*(-U[0,0]+gain_pl_sqrt(kb+S1 - beta*A[0,0]*U[1,0] +alpha*U[0,0] + Z[0,0]*(1+mu*U[0,0])))/tau_u
U[1,1] = U[1,0]+dt*(-U[1,0]+gain_pl_sqrt(kb+S2 - beta*A[1,0]*U[0,0] +alpha*U[1,0] + Z[1,0]*(1+mu*U[1,0])))/tau_u
U[0,0] = U[0,1]
U[1,0] = U[1,1]

A[0,0] = A[0,0] + dt*(1 - A[0,0] - gamma*A[0,0]*U[0,0])/tau_a
A[1,0] = A[1,0] + dt*(1 - A[1,0] - gamma*A[1,0]*U[1,0])/tau_a

u1.append(U[0,0])
u2.append(U[1,0])
a1.append(A[0,0])
a2.append(A[1,0])
z1.append(Z[0,0])
z2.append(Z[1,0])

#fast version
#Z[0,0] = Z[0,0] -dt*(Z[0,0] + np.sqrt(dt)*sigma*np.random.randn())
#Z[1,0] = Z[1,0] -dt*(Z[1,0] + np.sqrt(dt)*sigma*np.random.randn())

#U[0,1] = gain_pl_sqrt(kb+S1 - beta*A[0,0]*U[1,0] +alpha*U[0,0] + Z[0,0]*(1+mu*U[0,0]))/tau_u
#U[1,1] = gain_pl_sqrt(kb+S2 - beta*A[1,0]*U[0,0] +alpha*U[0,1] + Z[1,0]*(1+mu*U[1,0]))/tau_u
#U[0,0]= U[0,1]
#U[1,0]= U[1,1]

#A[0,0] = A[0,0] + dt*(1 - A[0,0] - gamma*A[0,0]*U[0,0])/tau_a
#A[1,0] = A[1,0] + dt*(1 - A[1,0] - gamma*A[1,0]*U[1,0])/tau_a

#u1.append(U[0,0])
#u2.append(U[1,0])
#a1.append(A[0,0])
#a2.append(A[1,0])
#z1.append(Z[0,0])
#z2.append(Z[1,0])

