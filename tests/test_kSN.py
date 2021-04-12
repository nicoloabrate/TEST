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
import matplotlib.pyplot as plt

nev = 1
M = -50
N = 1
G = 1
bc = 'Mark'
H = 0.1
xlayers = [0, H]
# define geometry and mesh
myslabP = Slab(M, xlayers, ['tmp'], [bc], G, N, 'FD')
# myslabD = Slab(M, xlayers, ['tmp'], ['zero'], G, 0, 'FD')
myslabS = Slab(M, xlayers, ['tmp'], [bc], G, N+1, 'FV')

# myDiff = NTE.Diffusion(myslabD, steady=True, fmt='csc', BC=True)
myPN = NTE.PN(myslabP, N, steady=True, fmt='csc', BC=True)
mySN = NTE.SN(myslabS, N+1, steady=True, fmt='csc', BC=True)

# transform SN operators into PN operators to check consistency
# sD = eigenproblem(myDiff, 'kappa', myslabD)
sP = eigenproblem(myPN, 'kappa', myslabP)
sS = eigenproblem(mySN, 'kappa', myslabS)
# sD.solve(algo='eig')
sP.solve()
sS.solve()

print('Fundamental k-eigenvalue:')
# print('Diffusion: {}'.format(sD.solution.eigvals[0]))
print('P1: {}'.format(sP.solution.eigvals[0]))
print('S2: {}'.format(sS.solution.eigvals[0]))

phiSN = sS.solution.eigvect
totphiSN = phiSN[0:myslabS.nS-1]+np.flipud(phiSN[myslabS.nS-1:2*myslabS.nS-1])

plt.plot(myslabS.stag_mesh,  phiSN[0:myslabS.nS-1, 0], label='phi1 SN-FV')
plt.plot(myslabS.stag_mesh, np.flipud(phiSN[myslabS.nS-1:2*myslabS.nS-1, 0]), label='phi2 SN-FV')
plt.legend()
plt.show()

# plt.plot(myslabD.mesh, sD.solution.eigvect[0:myslabS.nS, 0], label='Diffusion')
plt.plot(myslabP.mesh, sP.solution.eigvect[0:myslabS.nS, 0], label='P1')
plt.plot(myslabS.stag_mesh, totphiSN[:, 0], label='SN-FV')
plt.legend()
plt.show()

