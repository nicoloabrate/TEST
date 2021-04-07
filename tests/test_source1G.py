#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 13:00:55 2021

@author: nabrate
"""
import sys
sys.path.append('/opt/programs_nabrate/mycodes')
import numpy as np
from TEST.geometry import Slab
import TEST.models.NeutronTransportEquation as NTE
from TEST.models.EigenProblem import eigenproblem
from TEST.models.SourceProblem import sourceproblem
import matplotlib.pyplot as plt

nev = 1
M = 100
N = 1
G = 1
bc = 'Mark'
H = 10
xlayers = [0, H]
# define geometry and mesh
myslab = Slab(M, xlayers, ['tmp2'], [bc], G, N, 'FD')
myslabD = Slab(M, xlayers, ['tmp2'], ['zero'], G, 0, 'FD')
myslab2 = Slab(M, xlayers, ['tmp2'], [bc], G, N+1, 'FD')

myDiff = NTE.Diffusion(myslabD, steady=True, fmt='csc', BC=True)
myPN = NTE.PN(myslab, N, steady=True, fmt='csc', BC=True)
mySN = NTE.SN(myslab2, N+1, steady=True, fmt='csc', BC=True)

# transform SN operators into PN operators to check consistency
S0 = 100
mysrc = lambda x: S0 # *(x>=H/4 and x<=3*H/4)
sD = sourceproblem(myDiff, 'static', myslabD, mysrc)
sP = sourceproblem(myPN, 'static', myslab, mysrc)
sS = sourceproblem(mySN, 'static', myslab2, mysrc)
sD.solve()
sP.solve()
sS.solve()

L = myslabD.regions['tmp2'].getxs('DiffLength')
XSA = myslabD.regions['tmp2'].getxs('Abs')
a = H/4
b = H/2
A = S0/XSA/(1/np.tanh((b-a)/L)*np.cosh(a/L)+np.sinh(a/L))
C = -A*np.cosh(a/L)/np.sinh((b-a)/L)
diff_anly = lambda x: (x>=0 and x<=H/4)*A*np.sinh(x/L)+(x>=H/4 and x<=H/2)*(S0/XSA+C*np.cosh((b-x)/L))
tmp = []
for x in np.linspace(0, H/2, 1000):
    tmp.append(diff_anly(x)[0][0]*20)
# source problem
# Q = np.zeros(((N+1)*myslab2.nS,))
# Q[0:myslab2.nS] = 1
phiSN = sS.solution.flux
totphiSN = phiSN[0:myslab2.nS]+np.flipud(phiSN[myslab2.nS:2*myslab2.nS])
# plt.plot(myslabD.mesh, sD.solution.flux, label='Diff')
plt.plot(myslab.mesh, sP.solution.flux[0:myslab.nS], label='PN')
plt.plot(myslab2.mesh, totphiSN, label='SN')
# plt.plot(np.linspace(0, H/2, 1000), tmp, label='Diff. analytical')
plt.legend()