#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 13:00:55 2021

@author: nabrate
"""
import sys
sys.path.append('/opt/programs_nabrate/mycodes')
sys.path.append('C:\\Users\\39346\\Documenti\\mycodes\\TEST')
from TEST.geometry import Slab
import TEST.models.NeutronTransportEquation as NTE
from TEST.models.EigenProblem import eigenproblem
import matplotlib.pyplot as plt

algo = 'SLEPc'
M = -10
N = 1
G = 1
bc = 'Mark'
mat2 = 'Pu239_0L_noS'
mat1 = 'Absorber_0L'
H = 1.2
b = 0.9
xlayers = [-H, 0, b]
# define geometry and mesh
myslabP = Slab(M, xlayers, [mat1, mat2], [bc], G, N, 'FD')
myslabSFV = Slab(M, xlayers, [mat1, mat2], [bc], G, N+1, 'FV')  #[mat1, mat2, mat1]
myslabSFD = Slab(M, xlayers, [mat1, mat2], [bc], G, N+1, 'FD')

# myDiff = NTE.Diffusion(myslabD, steady=True, fmt='csc', BC=True)
myPN = NTE.PN(myslabP, N, steady=True, fmt='csc', BC=True)
mySNFV = NTE.SN(myslabSFV, N+1, steady=True, fmt='csc', BC=True)
mySNFD = NTE.SN(myslabSFD, N+1, steady=True, fmt='csc', BC=True)

# transform SN operators into PN operators to check consistency
sP = eigenproblem(myPN, 'kappa', myslabP)
sS_FV = eigenproblem(mySNFV, 'kappa', myslabSFV)
# sS_FD = eigenproblem(mySNFD, 'kappa', myslabSFD)

# sD.solve(algo='eig')
sP.solve(normalisation='peaktotalflux', algo=algo)
sS_FV.solve(normalisation='peaktotalflux', algo=algo)
# sS_FD.solve(normalisation='peaktotalflux', algo=algo)

print('Fundamental k-eigenvalue:')
# print('Diffusion: {}'.format(sD.solution.eigvals[0]))
print('P{}: {}'.format(N, sP.solution.eigvals[0]))
print('S{} FV: {}'.format(N+1, sS_FV.solution.eigvals[0]))
print('diff: {} [pcm]'.format(1E5*abs(sS_FV.solution.eigvals[0]-sP.solution.eigvals[0]).real))

# print('S{} FD: {}'.format(N+1, sS_FD.solution.eigvals[0]))
# print('diff: {} [pcm]'.format(1E5*abs(sS_FD.solution.eigvals[0]-sP.solution.eigvals[0]).real))

phiSN1_FV = sS_FV.solution.get(group=1, moment=0, mode=0)
# phiSN1_FD = sS_FD.solution.get(group=1, moment=0, mode=0)
phiPN1 = sP.solution.get(group=1, moment=0, mode=0)
# phiSN2 = sS.solution.get(group=2, moment=0, mode=0)
# phiPN2 = sP.solution.get(group=2, moment=0, mode=0)

# plt.plot(myslabD.mesh, sD.solution.eigvect[0:myslabS.nS, 0], label='Diffusion')
plt.plot(myslabP.mesh, phiPN1, label='P1')
plt.plot(myslabSFV.mesh, phiSN1_FV, label='SN-FV', linestyle='--')
# plt.plot(myslabSFD.mesh, phiSN1_FD, label='SN-FD', linestyle=':')
plt.legend()
plt.show()

# plt.plot(myslabP.mesh, phiPN2, label='P1')
# plt.plot(myslabS.mesh, phiSN2, label='SN-FV', linestyle='--')
# plt.legend()
# plt.show()
