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
xlayers = [-H, -b, b, H]
M = 50 # [-10, -5, -10]
N = 1
G = 1
bc = 'Mark'
mat2 = 'FullAbsorption'
mat1 = 'Absorber_0L'  #'FullAbsorption' #
# xlayers = [-H, H]
# define geometry and mesh   [mat2, mat1, mat2]
myslabP = Slab(M, xlayers, [mat2, mat1, mat2], [bc], G, N, 'FD')
myslabD = Slab(M, xlayers, [mat2, mat1, mat2], ['zero'], G, 0, 'FD')
myslabS = Slab(M, xlayers, [mat2, mat1, mat2], [bc], G, N+1, 'FD')

myDiff = NTE.Diffusion(myslabD, steady=True, fmt='csc', BC=True)
myPN = NTE.PN(myslabP, N, steady=True, fmt='csc', BC=True)
mySN = NTE.SN(myslabS, N+1, steady=True, fmt='csc', BC=True)

# transform SN operators into PN operators to check consistency
S0 = 1
mysrc = lambda x: S0*(x>=-b and x<=b)
sD = sourceproblem(myDiff, 'static', myslabD, mysrc)
sP = sourceproblem(myPN, 'static', myslabP, mysrc)
sS = sourceproblem(mySN, 'static', myslabS, mysrc)
sD.solve()
sP.solve()
sS.solve()

phiD = sD.solution.get(group=1, moment=0, mode=0)
phiP1 = sP.solution.get(group=1, moment=0, mode=0)
phiS2 = sS.solution.get(group=1, moment=0, mode=0)
# phiSN1_FD = sS_FD.solution.get(group=1, moment=0, mode=0)
# phiSN2 = sS.solution.get(group=2, moment=0, mode=0)
# phiPN2 = sP.solution.get(group=2, moment=0, mode=0)

plt.plot(myslabD.mesh, phiD, label='Diffusion')
plt.plot(myslabP.mesh, phiP1, label='P1')
plt.plot(myslabS.mesh, phiS2, label='S2', linestyle='--')
# plt.plot(myslabSFD.mesh, phiSN1_FD, label='SN-FD', linestyle=':')
plt.legend()
plt.show()

# check maxima
print('D  {:3f}'.format(phiD.max()))
print('P1 {:3f}'.format(phiP1.max()))
print('S2 {:3f}'.format(phiS2.max()))

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