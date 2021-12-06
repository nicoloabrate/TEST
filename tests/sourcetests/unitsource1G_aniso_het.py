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

H = 5
b = 1
xlayers = [-H, H]
M = -30 # [-10, -5, -10] -b, b,
N = 1
G = 1
bc = 'Mark'
mat2 = 'FullAbsorption'
mat1 = 'FullAbsorption' # 'Absorber_0L'  #
# xlayers = [-H, H]
# define geometry and mesh   [mat2, mat1, mat2]
myslabP = Slab(M, xlayers, mat1, [bc], G, N, 'FD')
myslabS = Slab(M, xlayers, mat1, [bc], G, N+1, 'FD')

myPN = NTE.PN(myslabP, N, steady=True, fmt='csc', BC=True)
mySN = NTE.SN(myslabS, N+1, steady=True, fmt='csc', BC=True)

# transform SN operators into PN operators to check consistency
S0 = 1
mysrc = lambda x, mu: S0*mu>0 # *(x>=-b and x<=b)
sP = sourceproblem(myPN, 'static', myslabP, mysrc)
sS = sourceproblem(mySN, 'static', myslabS, mysrc)

sP.solve()
sS.solve()

phiP1 = sP.solution.get(group=1, moment=0, mode=0)
phiS2 = sS.solution.get(group=1, moment=0, mode=0)

plt.plot(myslabP.mesh, phiP1, label='P1')
plt.plot(myslabS.mesh, phiS2, label='S2', linestyle='--')
# plt.plot(myslabSFD.mesh, phiSN1_FD, label='SN-FD', linestyle=':')
plt.legend()
plt.show()

# check maxima
print('P1 {:3f}'.format(phiP1.max()))
print('S2 {:3f}'.format(phiS2.max()))


sS.solution.xplot(1, angle=0, marker='s', label='S mu1')
sS.solution.xplot(1, angle=1, marker='s', label='S -mu1')
sP.solution.xplot(1, angle=0, label='P mu1', ls='--', marker='o')
sP.solution.xplot(1, angle=1, label='P -mu1', ls='-.', marker='o')
plt.legend()
plt.show()
# get operators
# Lp1=myPN.L.todense()
# Sp1=myPN.S.todense()
# Fp1=myPN.F.todense()
# Rp1=myPN.R.todense()
# S0p1=myPN.S0.todense()
# F0p1=myPN.F0.todense()

Ls2=mySN.L.todense()
Ss2=mySN.S.todense()
Fs2=mySN.F.todense()
Rs2=mySN.R.todense()
S0s2=mySN.S0.todense()
F0s2=mySN.F0.todense()
