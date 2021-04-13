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
M = -20
N = 6
G = 2
bc = 'Mark'
mat = 'Pu239a'
H = 1
xlayers = [0, H]
# define geometry and mesh
myslabP = Slab(M, xlayers, [mat], [bc], G, N, 'FD')
myslabS = Slab(M, xlayers, [mat], [bc], G, N+1, 'FV')

# myDiff = NTE.Diffusion(myslabD, steady=True, fmt='csc', BC=True)
myPN = NTE.PN(myslabP, N, steady=True, fmt='csc', BC=True)
mySN = NTE.SN(myslabS, N+1, steady=True, fmt='csc', BC=True)

# transform SN operators into PN operators to check consistency
# sD = eigenproblem(myDiff, 'kappa', myslabD)
sP = eigenproblem(myPN, 'kappa', myslabP)
sS = eigenproblem(mySN, 'kappa', myslabS)
# sD.solve(algo='eig')
sP.solve(normalisation='peaktotalflux')
sS.solve(normalisation='peaktotalflux')

print('Fundamental k-eigenvalue:')
# print('Diffusion: {}'.format(sD.solution.eigvals[0]))
print('P{}: {}'.format(N, sP.solution.eigvals[0]))
print('S{}: {}'.format(N+1, sS.solution.eigvals[0]))
print('diff: {} [pcm]'.format(1E5*abs(sS.solution.eigvals[0]-sP.solution.eigvals[0]).real))

phiSN1 = sS.solution.get(group=1, moment=0, mode=0)
phiPN1 = sP.solution.get(group=1, moment=0, mode=0)
# phiSN2 = sS.solution.get(group=2, moment=0, mode=0)
# phiPN2 = sP.solution.get(group=2, moment=0, mode=0)

# plt.plot(myslabD.mesh, sD.solution.eigvect[0:myslabS.nS, 0], label='Diffusion')
plt.plot(myslabP.mesh, phiPN1, label='P1')
plt.plot(myslabS.mesh, phiSN1, label='SN-FV', linestyle='--')
plt.legend()
plt.show()

# plt.plot(myslabP.mesh, phiPN2, label='P1')
# plt.plot(myslabS.mesh, phiSN2, label='SN-FV', linestyle='--')
# plt.legend()
# plt.show()

